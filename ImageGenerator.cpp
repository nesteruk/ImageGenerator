#include "lodepng.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <string>
#include <cmath>
#include <ctime>
using namespace std;

#include <boost/math/special_functions.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/bind.hpp>
#include <boost/tuple/tuple.hpp>
using namespace boost::accumulators;
using namespace boost::math;

#include <ppl.h>
using namespace concurrency;

#define UNARY 0
#define BINARY 1
#define TERMINALVALUE 2
#define TERMINALINDEX 3

#define UNARYCOUNT 3
#define BINARYCOUNT 3

#define randfloat ((float)rand() / (float)RAND_MAX)
#define randfloat2 (randfloat * 2.f - 1.f)

#ifndef nullptr
#define nullptr NULL
#endif

float Unary(float x, int n)
{
  switch (n)
  {
  case 0:
    return sin(x);
  case 1:
    return cos(x);
  case 2:
    return tanh(x);
  }
}

float Binary(float x, float y, int n)
{
  switch (n)
  {
  case 0:
    return x+y;
  case 1:
    return x-y;
  case 2:
    return x*y;
  }
}

class RandomFunction
{
  int type;
  int index;
  float value;
  RandomFunction* lhs;
  RandomFunction* rhs;
public:
  RandomFunction(int count, int depth)
  {
    lhs = rhs = NULL;
    if (depth == 0)
    {
      if (rand() % 50 > 10)
      {
        index = rand() % count;
        type = TERMINALINDEX;
        wcout << L"x[" << index << L"]";
      } 
      else
      {
        value = randfloat;
        type = TERMINALVALUE;
        wcout << value << endl;
      }
    } 
    else 
    {
      int r = rand() % 100;
      if(r > 80)
      {
        lhs = new RandomFunction(count, depth-1);
        type = UNARY;
        index = rand() % UNARYCOUNT;
      }
      else
      {
        lhs = new RandomFunction(count, depth-1);
        rhs = new RandomFunction(count, depth-1);
        type = BINARY;
        index = rand() % BINARYCOUNT;
      }
    }
  }

  ~RandomFunction()
  {
    if (lhs != nullptr) delete lhs;
    if (rhs != nullptr) delete rhs;
  }

  float Eval(float input)
  {
    float d[] = { input };
    return Eval(d);
  }

  float Eval(float params[])
  {
    switch (type)
    {
    case UNARY:
      return Unary(lhs->Eval(params), index);
    case BINARY:
      float a,b;
      a = lhs->Eval(params);
      b = rhs->Eval(params);
      return Binary(a, b, index);
    case TERMINALINDEX:
      return params[index];
    default:
      return value;
    }
  }
};

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
  return (unsigned char) std::max(0.0f, std::min(c,255.0f));
}

float poly(float x, float* args)
{
  return args[0]*x*x*x + args[1]*x*x + args[2]*x;
}

boost::tuple<float,float> meanAndSD(float* data, int count)
{
  accumulator_set<float, stats<tag::mean, tag::variance>> acc;
  for_each(data, data+count, [&](float d) { acc(d); });
  float avg = mean(acc);
  float sd = sqrt(variance(acc));
  assert(!isnan(avg) && !isinf(avg));
  assert(!isnan(sd) && !isinf(sd));
  return boost::tuple<float,float>(avg,sd);
}

void normalizeImageData(float* image, int count)
{
  auto msd = meanAndSD(image,count);
  const float avg = msd.get<0>();
  const float sd = msd.get<1>();
  const float shift = randfloat2;

  parallel_for_each(image, image+count, [=](float& d) { 
    d = 0.5 + (d - avg) / (sd * (3.0+shift)); 
  });
}



// module load intel/mpi
// mpicc -cxx=icpc
// cat /proc/meminfo
// cat /proc/cpuinfo
int main(int argc, char* argv[])
{
  unsigned seed = (unsigned)time(NULL);
  srand(seed);

  int w, h, depth, iterations;
  if (argc < 3)
  {
    w = 1920;
    h = 1080;
    depth = 8;
    iterations = 4;
  } 
  else
  {
    w = atoi(argv[1]);
    h = atoi(argv[2]);
    depth = atoi(argv[3]);
    iterations = atoi(argv[4]);
    wcout << L"Requested image is " << w << "x" << h << endl;
  }

start:
  float *data = new float[w*h];
  float *rgbData = new float[w*h*3];
  std::fill_n(rgbData,w*h*3,0.f);
  unsigned char* imageData = new unsigned char[w*h*3];
  std::fill_n(imageData,w*h*3,0);

  float *x = new float[w];
  float *y = new float[h];
  linspace(-10.f,10.f,w,x);
  linspace(-10.f,10.f,h,y);

  RandomFunction rf(2,depth);

  #pragma omp parallel for
  for (int j = 0; j < h; ++j)
  {
    wcout << j << endl;
    for (int i = 0; i < w; ++i)
    {
      int n = j*w+i;
      float xy[] = { x[i], y[j] };
      data[n] = rf.Eval(xy);
      assert(!isnan(data[n]) && !isinf(data[n]));
    }
  }

  normalizeImageData(data,w*h);


  for (int iteration = 0; iteration < iterations; ++iteration)
  {
    wcout << L"Iteration " << iteration << endl;

    float rco[] = { randfloat2, randfloat2, randfloat2 };
    float gco[] = { randfloat2, randfloat2, randfloat2 };
    float bco[] = { randfloat2, randfloat2, randfloat2 };
    float weight = 1.0 / (float)(1 << iteration);

    #pragma omp parallel for
    for (int j = 0; j < h; ++j)
    {
      for (int i = 0; i < w; ++i)
      {
        int n = j*w+i;
        int n3 = 3*n;
        float x = data[n];
        assert(!isnan(x) && !isinf(x));
        rgbData[n3] += 2*weight * poly(x, rco);
        rgbData[n3+1] += weight * poly(x, gco);
        rgbData[n3+2] += weight * poly(x, bco);
      }
    }
  }
  
  //normalizeImageData(rgbData,w*h*3);

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