#pragma once

typedef function<float(const float)> unaryFunction;
typedef function<float(const float,const float)> binaryFunction;

inline vector<unaryFunction> initUnaryFunctions()
{
  vector<unaryFunction> result { 
    sinf, cosf, expf, 
    logf, 
    sinhf, coshf, tanhf, 
    fabsf 
  };
  result.push_back([](float x) { return -x; });
  result.push_back([](float x) { return x*x; });
  result.push_back([](float x) { return x*x*x; });
  return result;
}

inline vector<function<float(const float,const float)>> initBinaryFunctions()
{
  vector<binaryFunction> result;
  result.push_back([](const float a, const float b) { return a + b; });
  result.push_back([](const float a, const float b) { return a - b; });
  result.push_back([](const float a, const float b) { return a * b; });
    
  //result.push_back([](float a, float b) { return a / b; });
  result.push_back([](const float a, const float b) { return sin(a*b); });
  //result.push_back([](float a, float b) { return cos(a*b); });
  //result.push_back([](float a, float b) { return tanh(a*b); });
  //result.push_back([](float a, float b) { return powf(a,b); });
  //result.push_back([](float a, float b) { return std::max(a,b); });
  //result.push_back([](float a, float b) { return std::min(a,b); });
  return result;
}

const vector<unaryFunction> unaryFunctions = initUnaryFunctions();
const vector<binaryFunction> binaryFunctions = initBinaryFunctions();

typedef function<void(const int*, const float*, float*)> unaryVectorFunction;
typedef function<void(const int*, const float*, const float*, float*)> binaryVectorFunction;

#ifdef USEMKL

inline vector<unaryVectorFunction> initUnaryVectorFunctions()
{
  vector<unaryVectorFunction> result;
  result.push_back(vssin);
  result.push_back(vscos);
  result.push_back(vsexp);
  //result.push_back(vslog10);
  result.push_back(vssinh);
  result.push_back(vscosh);
  result.push_back(vstanh);
  result.push_back(vsabs);
  //result.push_back(vssqrt);
  //result.push_back(vstanh);
  return result;
}

inline vector<binaryVectorFunction> initBinaryVectorFunctions()
{
  vector<binaryVectorFunction> result;
  result.push_back(vsadd);
  result.push_back(vssub);
  result.push_back(vsmul);

  // this is what a compound function (in this case, sin(a*b) looks like)
  result.push_back([](const int* n, const float* a, const float* b, float* c)
  {
    float* temp = new float[*n];
    vsmul(n, a, b, temp);
    vssin(n, temp, c);
    delete[] temp;
  });

  return result;
}

#endif

const vector<unaryVectorFunction> unaryVectorFunctions
#ifdef USEMKL
	= initUnaryVectorFunctions()
#endif
;
const vector<binaryVectorFunction> binaryVectorFunctions
#ifdef USEMKL
	= initBinaryVectorFunctions()
#endif
;