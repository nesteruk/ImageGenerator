#include "stdafx.h"
#include "../RandomFunction.hpp"
#include "Coloring.h"

void linspace(const float start, const float end, 
  const int count, float* result)
{
  float dx = (end - start) / (count - 1);
  for (size_t i = 0; i < count; ++i)
  {
    result[i] = start + dx*i;
  }
}

float poly(const float x, const float const* args)
{
  return x*(args[0] * x + args[1]) + args[2];
}

void normalizeImage(float *image, const size_t pixelCount)
{
  vector<float> r(pixelCount);
  vector<float> g(pixelCount);
  vector<float> b(pixelCount);

  accumulator_set<float, stats<tag::mean, tag::variance, tag::min, tag::max>> accR, accG, accB;

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

  // now transform the image
  for (size_t i = 0; i < pixelCount; i++)
  {
    image[i * 3 + 0] = (image[i * 3 + 0] - meanR) / (varR/3);
    image[i * 3 + 1] = (image[i * 3 + 1] - meanG) / (varG/3);
    image[i * 3 + 2] = (image[i * 3 + 2] - meanB) / (varB/3);
  }

  // postprocessing, you really feel like it
  /*for (size_t i = 0; i < pixelCount*3; i++)
  {
    image[i] *= image[i];
  }*/
}

inline uint8_t constrain(float x)
{
  return fmaxf(0.f, fminf(x, 255.f));
}

int main(int, char* argv[])
{
  float xmin = -M_PI, xmax = M_PI;
  float ymin = -M_PI, ymax = M_PI;
  const int w = 720;
  const int h = 720;
  const int count = 2; // x and y
  int imagesToGenerate = 100;
  random_device rd;
everything:
  //unsigned seed = rd();
  auto seed = (unsigned)time(nullptr);
  //auto seed = 3584265032;
  srand(seed);
  int depth = 8;
  int iterations = 1;
start:
  cout << "Render started with seed " << seed << endl;
  RandomFunction rf(count, depth);
regen:
  float *x = new float[w];
  float *y = new float[h];
  linspace(xmin, xmax, w, x);
  linspace(ymin, ymax, h, y);
  
  cout << rf << endl;
  float *f = new float[w * h];

  uint8_t *image = new uint8_t[w * h * 3];
  float *intermediate = new float[w * h * 3];
  memset(image, 0, w*h * 3 * sizeof(uint8_t));

  int it_count = iterations;
  while (it_count --> 0) 
  {
    #pragma omp parallel for
    for (size_t j = 0; j < h; ++j)
    {
      for (size_t i = 0; i < w; ++i)
      {
        f[j*w + i] = rf.Eval(x[i], y[j]);
      }
    }

    PolynomialColoring ca(3);
        
    // prepare intermediate values
    for (size_t j = 0; j < h; ++j)
    {
      for (size_t i = 0; i < w; ++i)
      {
        int k = j*w + i;
        int k3 = 3 * k;
        ca.Value(f[k], intermediate[k3], intermediate[k3 + 1], intermediate[k3 + 2]);
      }
    }
  };

  // normalize intermediate values using chosen algos
  normalizeImage(intermediate, w*h);

  for (size_t j = 0; j < h; ++j)
  {
    for (size_t i = 0; i < w; ++i)
    {
      int k = j*w + i;
      int k3 = 3 * k;
      image[k3 + 0] = constrain(128.f * (intermediate[k3 + 0] + 1.f));
      image[k3 + 1] = constrain(128.f * (intermediate[k3 + 1] + 1.f));
      image[k3 + 2] = constrain(128.f * (intermediate[k3 + 2] + 1.f));
    }
  }

  auto filename = string(".\\samples\\") + boost::lexical_cast<string>(seed) + ".png";
  lodepng_encode24_file(filename.c_str(), image, w, h);

  delete[] image;
  delete[] intermediate;
  delete[] f;
  delete[] x;
  delete[] y;
  
  if (--imagesToGenerate <= 0)
  {
    cout << "Render complete" << endl;

    cout << "0. Go again" << endl;
    cout << "1. Change xy scope" << endl;
    cout << "2. Change depth" << endl;
    cout << "3. Recolor" << endl;
    cout << "4. Change iterations" << endl;

    int input;
    cin >> input;

    if (input == 0)
    {
      goto everything;
    }
    else if (input == 1)
    {
      cout << "Enter xmin xmax ymin ymax" << endl;

      cin >> xmin >> xmax >> ymin >> ymax;
      goto regen;
    }
    else if (input == 2)
    {
      cout << "Enter new depth: ";
      cin >> depth;
      goto start;
    }
    else if (input == 3)
    {
      goto regen;
    }
    else if (input == 4)
    {
      cout << "Enter no. of iterations: ";
      cin >> iterations;
      goto regen;
    }
  }
  else goto everything;

	return 0;
}

