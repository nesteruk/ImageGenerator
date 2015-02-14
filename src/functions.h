#pragma once

typedef function<float(float)> unaryFunction;
typedef function<float(float,float)> binaryFunction;

inline vector<unaryFunction> initUnaryFunctions()
{
  vector<unaryFunction> result { sinf, cosf, expf, sinhf, coshf, tanhf };
  result.push_back([](float x) { return x*x; });
  result.push_back([](float x) { return x*x*x; });
  return result;
}

inline vector<function<float(float,float)>> initBinaryFunctions()
{
  vector<binaryFunction> result;
  result.push_back([](float a, float b) { return a + b; });
  result.push_back([](float a, float b) { return a - b; });
  result.push_back([](float a, float b) { return a * b; });

  /*
  result.push_back([](float a, float b) { return a / b; });
  result.push_back([](float a, float b) { return sin(a*b); });
  result.push_back([](float a, float b) { return cos(a*b); });
  result.push_back([](float a, float b) { return powf(a,b); });
  result.push_back([](float a, float b) { return std::max(a,b); });
  result.push_back([](float a, float b) { return std::min(a,b); });*/
  return result;
}

const vector<unaryFunction> unaryFunctions = initUnaryFunctions();
const vector<binaryFunction> binaryFunctions = initBinaryFunctions();

typedef function<void(int*, float*, float*)> unaryVectorFunction;
typedef function<void(int*, float*, float*, float*)> binaryVectorFunction;

inline vector<unaryVectorFunction> initUnaryVectorFunctions()
{
  vector<unaryVectorFunction> result;
  result.push_back(vssin);
  result.push_back(vscos);
  result.push_back(vssqrt);
  result.push_back(vstanh);
  return result;
}

inline vector<binaryVectorFunction> initBinaryVectorFunctions()
{
  vector<binaryVectorFunction> result;
  result.push_back([](int* n, float* a, float* b, float* c)
  {
    #pragma omp parallel for
    #pragma simd
    for (int i = 0; i < *n; ++i)
      c[i] = a[i] + b[i];
  });
  result.push_back([](int* n, float* a, float* b, float* c)
  {
    #pragma omp parallel for
    #pragma simd
    for (int i = 0; i < *n; ++i)
      c[i] = a[i] - b[i];
  });
  result.push_back([](int* n, float* a, float* b, float* c)
  {
    #pragma omp parallel for
    #pragma simd
    for (int i = 0; i < *n; ++i)
      c[i] = a[i] * b[i];
  });
  return result;
}

const vector<unaryVectorFunction> unaryVectorFunctions = initUnaryVectorFunctions();
const vector<binaryVectorFunction> binaryVectorFunctions = initBinaryVectorFunctions();