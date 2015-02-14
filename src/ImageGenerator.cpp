#include "stdafx.h"
#include "lodepng.h"
#include "functions.h"
#include "RandomFunction.hpp"

void linspace(float start, float end, int count, float* result)
{
  float dx = (end - start) / (count-1);
  for (int i = 0; i < count; ++i)
  {
    result[i] = start + dx*i;
  }
}

inline unsigned char getPixelColor(float value)
{
  float c = value * 255.0;
  return static_cast<unsigned char>(std::max(0.0f, std::min(c,255.0f)));
}

float poly(const float x, const float *const args)
{
  return x*(args[0]*x*x + args[1]*x + args[2]) + args[3];
}

boost::tuple<float,float,float,float> meanAndVariance(float* data, int count)
{
  accumulator_set<float, stats<tag::mean, tag::variance, tag::min, tag::max>> acc;
  for_each(data, data + count, [&](float d) { acc(d); });
  float avg = mean(acc);
  float var = variance(acc);
  float mn = boost::accumulators::min(acc);
  float mx = boost::accumulators::max(acc);
  return boost::tuple<float,float,float,float>(avg,var,mn,mx);
}

void normalizeImageData(float* image, int count)
{
  int width = count / 3;
  float* r = new float[width];
  float* g = new float[width];
  float* b = new float[width];

  for (int i = 0; i < (width); ++i)
  {
    r[i] = image[i * 3];
    g[i] = image[i * 3 + 1];
    b[i] = image[i * 3 + 2];
  }

  auto mvr = meanAndVariance(r, width);
  auto mvg = meanAndVariance(g, width);
  auto mvb = meanAndVariance(b, width);
  auto minr = get<2>(mvr);
  auto ming = get<2>(mvg);
  auto minb = get<2>(mvb);
  auto rr = get<3>(mvr) - minr;
  auto rg = get<3>(mvg) - ming;
  auto rb = get<3>(mvb) - minb;

  cout << "Min: " << minr << " " << ming << " " << minb << endl;
  cout << "Range: " << rr << " " << rg << " " << rb << endl;

  parallel_for_each(r, r + width, [=](float& d) {
    d = (d - minr) / rr;
  });

  parallel_for_each(g, g + width, [=](float& d) {
    d = (d - ming) / rg;
  });

  parallel_for_each(b, b+ width, [=](float& d) {
    d = (d - minb) / rb;
  });

  // here we introduce a bit of a saturation bias

  for (int i = 0; i < (width); ++i)
  {
    image[i * 3] = r[i];
    image[i * 3 + 1] = g[i];
    image[i * 3 + 2] = b[i];
  }

  delete[] r;
  delete[] g;
  delete[] b;
}

void generateData(float* data, int w, int h, float* x, float* y, int dimensions, 
                  RandomFunction& rf, const vector<float>& coeff, bool warp)
{
  #pragma omp parallel for
  for (int j = 0; j < h; ++j)
  {
    for (int i = 0; i < w; ++i)
    {
      vector<float> params;
      params.push_back(x[i] * coeff[0]);
      params.push_back(y[j] * coeff[1]);
      
      int n = j*w+i;
      data[n] = rf.Eval(params);
    }
  }
}

void generateDataVector(float* data, int w, int h, float* x, float* y, int dimensions,
  RandomFunction& rf, const vector<float>& coeff, bool warp)
{
  int pts = w*h;

  // we need a meshgrid of x-y values
  float *xx = new float[pts];
  float *yy = new float[pts];
  for (size_t j = 0; j < h; ++j)
  {
    for (size_t i = 0; i < w; ++i)
    {
      xx[j*w + i] = x[i];
      yy[j*w + i] = y[i];
    }
  }

  rf.Eval(&pts, xx, yy, data);

  delete[] xx;
  delete[] yy;
}

// module load intel/mpi
// mpicc -cxx=icpc
// cat /proc/meminfo
// cat /proc/cpuinfo
int w, h, depth, iterations, dimensions;

void generate();
void generate_vector();

int main(int argc, char* argv[])
{
  if (argc < 3)
  {
    w = 720;
    h = 720;
    depth = 8;
    iterations = 2;
    dimensions = 2;
  } 
  else
  {
    w = atoi(argv[1]);
    h = atoi(argv[2]);
    depth = atoi(argv[3]);
    iterations = atoi(argv[4]);
    dimensions = atoi(argv[5]);
    wcout << L"Requested image is " << w << "x" << h << endl;
  }

  generate();
}

void generate_vector()
{
  unsigned seed;

start:
  seed = static_cast<unsigned>(time(nullptr));
  //seed = 1413732061;
  srand(seed);
  float *data = new float[w*h];
  float *rgbData = new float[w*h * 3];
  std::fill_n(rgbData, w*h * 3, 0.f);
  unsigned char* imageData = new unsigned char[w*h * 3];
  std::fill_n(imageData, w*h * 3, 0);

  float *x = new float[w];
  float *y = new float[h];

  auto xRange = make_pair(-1, 1);
  auto yRange = make_pair(-1, 1);

  linspace(xRange.first, xRange.second, w, x);
  linspace(yRange.first, yRange.second, h, y);

  RandomFunction rf(dimensions, depth);

  vector<float> coeffs(dimensions);
  fill_n(coeffs.begin(), dimensions, 1.0f);

  float adjustmentFactor = 1.0;

  for (int iteration = 0; iteration < iterations; ++iteration)
  {
    cout << "Iteration " << iteration << endl;

    float rco[] = { randfloat2, randfloat2, randfloat2 };
    float gco[] = { randfloat2, randfloat2, randfloat2 };
    float bco[] = { randfloat2, randfloat2, randfloat2 };

    // grab some random data
    generateDataVector(data, w, h, x, y, dimensions, rf, coeffs, false);

    //normalizeImageData(data, w*h);

    // augment the rgb values
    #pragma omp parallel for
    for (int j = 0; j < h; ++j)
    {
      for (int i = 0; i < w; ++i)
      {
        int n = j*w + i;
        int n3 = 3 * n;
        float z = data[n];
        rgbData[n3 + 0] += adjustmentFactor * poly(z, rco);
        rgbData[n3 + 1] += adjustmentFactor * poly(z, gco);
        rgbData[n3 + 2] += adjustmentFactor * poly(z, bco);
      }
    }

    // change the data and try again
    /*for (int i = 0; i < dimensions; ++i)
    coeffs[i] = randfloat;*/

    // change the factor
    adjustmentFactor = randfloat;
  }
  normalizeImageData(rgbData, w*h * 3);
  auto mv = meanAndVariance(rgbData, w*h * 3);
  cout << "Postprocessed mean = " << mv.get<0>() << " var = " << mv.get<1>() << endl;

#pragma omp parallel for
  for (int j = 0; j < h; ++j)
  {
    for (int i = 0; i < w; ++i)
    {
      for (int k = 0; k < 3; ++k)
      {
        int n = 3 * (j*w + i) + k;
        imageData[n] = getPixelColor(rgbData[n]);
      }
    }
  }

  lodepng_encode24_file("a.png", imageData, w, h);
  wcout << L"Image rendering complete, seed=" << seed << endl;

  delete[] x;
  delete[] y;
  delete[] data;
  delete[] imageData;
  delete[] rgbData;

  getchar();
  goto start;
}


void generate()
{
  unsigned seed;

start:
  seed = static_cast<unsigned>(time(nullptr));
  //seed = 1413732061;
  srand(seed);
  float *data = new float[w*h];
  float *rgbData = new float[w*h*3];
  std::fill_n(rgbData,w*h*3,0.f);
  unsigned char* imageData = new unsigned char[w*h*3];
  std::fill_n(imageData,w*h*3,0);

  float *x = new float[w];
  float *y = new float[h];

  auto xRange = make_pair(-1, 1);
  auto yRange = make_pair(1, 1);

  linspace(xRange.first,xRange.second,w,x);
  linspace(yRange.first,yRange.second,h,y);

  RandomFunction rf(dimensions, depth);

  vector<float> coeffs(dimensions);
  fill_n(coeffs.begin(),dimensions,1.0f);

  float adjustmentFactor = 1.0;
  for (int iteration = 0; iteration < iterations; ++iteration)
  {
    cout << "Iteration " << iteration << endl;

    float rco[] = { randfloat2, randfloat2, randfloat2 };
    float gco[] = { randfloat2, randfloat2, randfloat2 };
    float bco[] = { randfloat2, randfloat2, randfloat2 };

    // grab some random data
    generateData(data,w,h,x,y,dimensions,rf,coeffs, false);

    //normalizeImageData(data, w*h);

    // augment the rgb values
    #pragma omp parallel for
    for (int j = 0; j < h; ++j)
    {
      for (int i = 0; i < w; ++i)
      {
        int n = j*w + i;
        int n3 = 3 * n;
        float z = data[n];
        rgbData[n3 + 0] += adjustmentFactor * poly(z, rco);
        rgbData[n3 + 1] += adjustmentFactor * poly(z, gco);
        rgbData[n3 + 2] += adjustmentFactor * poly(z, bco);
      }
    }

    // change the data and try again
    /*for (int i = 0; i < dimensions; ++i)
      coeffs[i] = randfloat;*/

    // change the factor
    adjustmentFactor = randfloat;
  }
  normalizeImageData(rgbData,w*h*3);
  auto mv = meanAndVariance(rgbData, w*h*3);
  cout << "Postprocessed mean = " << mv.get<0>() << " var = " << mv.get<1>() << endl;

  #pragma omp parallel for
  for (int j = 0; j < h; ++j)
  {
    for (int i = 0; i < w; ++i)
    {
      for (int k = 0; k < 3; ++k)
      {
        int n = 3*(j*w+i)+k;
        imageData[n] = getPixelColor(rgbData[n]);
      }
    }
  }

  lodepng_encode24_file("a.png", imageData, w, h);
  wcout << L"Image rendering complete, seed=" << seed << endl;

  delete[] x;
  delete[] y;
  delete[] data;
  delete[] imageData;
  delete[] rgbData;

  getchar();
  goto start;
}
