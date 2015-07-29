#include "stdafx.h"
#include "RandomFunction.hpp"
#include "Coloring.h"
namespace po = boost::program_options;

void linspace(const float start, const float end, 
  const int count, float* result)
{
  float dx = (end - start) / (count - 1);
  for (size_t i = 0; i < count; ++i)
  {
    result[i] = start + dx*i;
  }
}

float poly(const float x, const float *const args)
{
  return x*(args[0] * x + args[1]) + args[2];
}

// returns true if normalization succeeded and false if it failed
bool normalizeImage(float *image, const size_t pixelCount)
{
  vector<float> r(pixelCount);
  vector<float> g(pixelCount);
  vector<float> b(pixelCount);

  accumulator_set<float, stats<tag::mean, tag::variance,
          tag::min, tag::max>> accR, accG, accB;

  for (size_t i = 0; i < pixelCount; i++)
  {
    accR(image[i * 3 + 0]);
    accG(image[i * 3 + 1]);
    accB(image[i * 3 + 2]);
  }
      
  float meanR = mean(accR), meanG = mean(accG), meanB = mean(accB);
  float varR = variance(accR), varG = variance(accG), varB = variance(accB);
  float minR = boost::accumulators::min(accR), minG = boost::accumulators::min(accG), minB = boost::accumulators::min(accB);
  float maxR = boost::accumulators::max(accR), maxG = boost::accumulators::max(accG), maxB = boost::accumulators::max(accB);

  float rangeR = maxR - minR;
  float rangeG = maxG - minG;
  float rangeB = maxB - minB;

  cout << "mean: " << meanR << " " << meanG << " " << meanB << endl;
  cout << "var : " << varR << " " << varG << " " << varB << endl;
  cout << "min : " << minR << " " << minG << " " << minB << endl;
  cout << "max : " << maxR << " " << maxG << " " << maxB << endl;

  float all[] = { meanR, meanG, meanB, varR, varG, varB };
  for (auto x : all)
    if (std::isnan(x) || std::isinf(x)) return false;

  // now transform the image
  for (size_t i = 0; i < pixelCount; i++)
  {
    image[i * 3 + 0] = (image[i * 3 + 0] - meanR) / (varR / 3);
    image[i * 3 + 1] = (image[i * 3 + 1] - meanG) / (varG / 3);
    image[i * 3 + 2] = (image[i * 3 + 2] - meanB) / (varB / 3);
  }

  // postprocessing, you really feel like it
  /*for (size_t i = 0; i < pixelCount*3; i++)
  {
    image[i] *= image[i];
  }*/

  return true;
}

inline uint8_t constrain(float x)
{
  return fmaxf(0.f, fminf(x, 255.f));
}

int equivalence_check()
{
  assert(unaryFunctions.size() == unaryVectorFunctions.size());
  assert(binaryFunctions.size() == binaryVectorFunctions.size());

  // ensure that scalar and vector random functions are equivalent
  constexpr int count = 101;
  float f[count];
  linspace(-M_PI, M_PI, count, f);
  float result[count*count], result2[count*count];

  RandomFunction rf(2, 8, false);
  rf.Eval(count, count, f, f, result, false);
  rf.Eval(count, count, f, f, result2, true);

  for (size_t i = 0; i < count; i++)
  {
    if (result[i] != result2[i])
    {
      cout << "failed at " << i << " difference = " << abs(result[i]-result2[i]) << endl;
    }
  }
  return 0;
}

static random_device rd{};

struct image_rendering_options
{
  float xmin = -M_PI;
  float xmax = M_PI;
  float ymin = -M_PI;
  float ymax = M_PI;
  const unsigned int dimensions = 2;
  unsigned int iterations = 2;
  unsigned int depth = 6;
  boost::optional<unsigned> seed;
  int w = 1024;
  int h = 1024;
  int threads = 1;
  unsigned int images_to_generate = 1; // per process

  friend class boost::serialization::access;

  template<class Ar>
  inline void serialize(Ar& ar, const unsigned int fv)
  {
    ar & xmin;
    ar & xmax;
    ar & ymin;
    ar & ymax;
    ar & dimensions;
    ar & iterations;
    ar & depth;
    ar & seed;
    ar & w;
    ar & h;
    ar & threads;
    ar & images_to_generate;
  }
};

image_rendering_options parse_options(int ac, char* av[])
{
  vector<double> renderTimes; // not here

  image_rendering_options iro;

  po::options_description desc("generator options");
  desc.add_options()
    ("help", "show some help")
    ("n", po::value<unsigned int>(&iro.images_to_generate), "number of images to generate")
    ("i", po::value<unsigned int>(&iro.iterations), "number of iterations")
    ("d", po::value<unsigned int>(&iro.depth), "function depth")
    ("z", po::value<int>(), "size of image to generate (w/h)")
    ("f", po::value<vector<float>>(), "frame (xmin xmax ymin ymax) in which to generate")
    ("s", po::value<string>(), "rng seed")
    ("t", po::value<int>(&iro.threads), "number of threads");

  po::variables_map vm;
  try
  {
    po::store(po::parse_command_line(ac, av, desc,
      po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
  }
  catch (const exception& e)
  {
    cout << "Cannot parse command-line arguments: " << e.what() << endl;
    throw;
  }
  po::notify(vm);

  if (vm.count("z")) iro.w = iro.h = vm["z"].as<int>();

  if (vm.count("t"))
  {
    //omp_set_num_threads(threads);
    // this is now specified later
  }

  if (vm.count("f"))
  {
    auto vals = vm["f"].as<vector<float>>();
    if (vals.size() == 4)
    {
      iro.xmin = vals[0];
      iro.xmax = vals[1];
      iro.ymin = vals[2];
      iro.ymax = vals[3];
    }
  }

  if (vm.count("help"))
  {
    cout << desc << endl;
    iro.images_to_generate = 0;
  }

  // the seed is used to initialize the RNG
  if (vm.count("s"))
  {
    iro.seed = (unsigned)boost::lexical_cast<unsigned>(vm["s"].as<string>());
  }

  return iro;
}

void render(image_rendering_options& iro) {
  vector<double> renderTimes;


    //#pragma omp parallel for
    for (int i = 0; i < iro.images_to_generate; ++i) {
    everything:
      auto seed = iro.seed.is_initialized() ? iro.seed.value() : rd();
      srand(seed);

    sanity:
      #ifdef USEMKL // does vector routine checks
      equivalence_check();
      //return;
      #endif

      cout << endl << "Render started with seed " << seed << endl;
      auto start_time = high_resolution_clock::now();
      RandomFunction rf(iro.dimensions, iro.depth, false);

      auto w = iro.w;
      auto h = iro.h;

      float *x = new float[w];
      float *y = new float[h];
      linspace(iro.xmin, iro.xmax, w, x);
      linspace(iro.ymin, iro.ymax, h, y);

      cout << rf << endl;
      float *f = new float[w * h];

      uint8_t *image = new uint8_t[w * h * 3];
      float *intermediate = new float[w * h * 3];
      memset(image, 0, w * h * 3 * sizeof(uint8_t));

      int it_count = iro.iterations;
      while (it_count-- > 0) {
        rf.Eval(w, h, x, y, f, false);

        PolynomialColoring ca(3);

        // prepare intermediate values
        for (size_t j = 0; j < h; ++j) {
          for (size_t i = 0; i < w; ++i) {
            int k = j * w + i;
            int k3 = 3 * k;
            ca.Value(f[k], intermediate[k3], intermediate[k3 + 1], intermediate[k3 + 2]);
          }
        }
      };

      cout << "Calculations done." << endl;

      // normalize intermediate values using chosen algos
      if (!normalizeImage(intermediate, w * h)) {
        cout << "Generation failed, re-doing this!" << endl;
        // if this fails, try again with same params
        goto everything;
      }

      for (size_t j = 0; j < h; ++j) {
        for (size_t i = 0; i < w; ++i) {
          int k = j * w + i;
          int k3 = 3 * k;
          image[k3 + 0] = constrain(128.f * (intermediate[k3 + 0] + 1.f));
          image[k3 + 1] = constrain(128.f * (intermediate[k3 + 1] + 1.f));
          image[k3 + 2] = constrain(128.f * (intermediate[k3 + 2] + 1.f));
        }
      }

      // send the image over to the recipient
      auto filename = string("./samples/") + boost::lexical_cast<string>(seed) + ".png";
      lodepng_encode24_file(filename.c_str(), image, w, h);

      delete[] image;
      delete[] intermediate;
      delete[] f;
      delete[] x;
      delete[] y;

      auto end_time = high_resolution_clock::now();
      auto duration = duration_cast<milliseconds>(end_time - start_time).count();
      renderTimes.push_back(duration);
      cout << "Render took " << duration << "msec." << endl;
    }

  // compute average render time
  double sum = 0.0;
  uint32_t renderCount = renderTimes.size();
  for (auto t : renderTimes) sum += t;
  cout << "Average rendering time is " << (sum / renderCount) << "msec." << endl;
}

int main(int ac, char* av[])
{
  auto options = parse_options(ac, av);
  options.images_to_generate = 50;

  render(options);

  return 0;
}