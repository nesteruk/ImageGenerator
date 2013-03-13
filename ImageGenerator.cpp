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
using namespace boost::accumulators;
using namespace boost::math;

#define UNARY 0
#define BINARY 1
#define TERMINALVALUE 2
#define TERMINALINDEX 3

#define UNARYCOUNT 8
#define BINARYCOUNT 5

#define randfloat ((float)rand() / (float)RAND_MAX)

float Unary(float x, int n)
{
  switch (n)
  {
  case 0:
    return sin(2.0f*M_PI*x);
  case 1:
    return cos(2.0f*M_PI*x);
  case 2:
    return (exp(x)-1)/(M_E-1);
  case 3:
    return sqrt(x);
  case 4:
    return x*x;
  case 5:
    return abs(x-0.5);
  case 6:
    return sqrt(x+1) / M_SQRT2;
  default:
    return x;
  }
}

float Binary(float x, float y, int n)
{
  switch (n)
  {
  case 0:
    return 0.5*(x+y);
  case 1:
    return 0.5*(x-y+1);
  case 2:
    return 0.5*(1+tanh(4*x-2));
  case 3:
    return log(x+1);
  default:
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
      if (rand() % 50 > 1)
      {
        index = rand() % count;
        type = TERMINALINDEX;
      } 
      else
      {
        value = randfloat;
        type = TERMINALVALUE;
      }
    } 
    else 
    {
      int r = rand() % 100;
      if(r > 50)
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
  return args[0]*x*x + args[1]*x + args[2];
}

std::tuple<float,float> meanAndSD(float* data, int count)
{
  accumulator_set<float, stats<tag::mean, tag::variance>> acc;
  for_each(data, data+count, [&](float d) { acc(d); });
  float avg = mean(acc);
  float sd = sqrt(variance(acc));
  return std::tuple<float,float>(avg,sd);
}

void normalizeImageData(float* image, int count)
{
  auto msd = meanAndSD(image,count);
  for_each(image, image+count, [=](float& d) { d = 0.5 + (d - std::get<0>(msd)) / (std::get<1>(msd) * 3.0); });
}



// module load intel/mpi
// mpicc -cxx=icpc
// cat /proc/meminfo
// cat /proc/cpuinfo
int main(int argc, char* argv[])
{
  srand((unsigned)time(NULL));

  int w, h, depth, iterations;
  if (argc < 3)
  {
    w = 1024;
    h = 768;
    depth = 10;
    iterations = 3;
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

  float rco[] = { randfloat, randfloat, randfloat };
  float gco[] = { randfloat, randfloat, randfloat };
  float bco[] = { randfloat, randfloat, randfloat };

  float *data = new float[w*h];
  float *rgbData = new float[w*h*3];
  std::fill_n(rgbData,w*h*3,0.0);
  unsigned char* imageData = new unsigned char[w*h*3];
  std::fill_n(imageData,w*h*3,0);

  float *x = new float[w];
  float *y = new float[h];
  linspace(0.f,1.f,w,x);
  linspace(0.f,1.f,h,y);

  RandomFunction rf(3,depth);
  RandomFunction rf2(3,depth);
  RandomFunction rf3(3,depth);

  #pragma omp parallel for
  for (int j = 0; j < h; ++j)
  {
    for (int i = 0; i < w; ++i)
    {
      int n = j*w+i;
      int n3 = 3*n;
      float xy[] = { x[i], y[j] };
      rgbData[n3] += rf.Eval(xy);
      rgbData[n3+1] += rf2.Eval(xy);
      rgbData[n3+2] += rf3.Eval(xy);
    }
  }

  normalizeImageData(rgbData,w*h*3);

  for (int i = 0; i < w*h*3; ++i)
    imageData[i] = getPixelColor(rgbData[i]);

  wcout << L"Calculation finished" << endl;
  lodepng_encode24_file("a.png", imageData, w, h);

  delete[] x;
  delete[] y;
  delete[] data;
  delete[] imageData;
  delete[] rgbData;
}