#pragma once

#include "functions.h"

#define randfloat ((float)rand() / (float)RAND_MAX)
#define randfloat2 ((2*randfloat) - 1.f)

#define UNARY 0
#define BINARY 1
#define TERMINALINDEX 2
#define TERMINALVALUE 3

class RandomFunction
{
  int type;
  int index;
  float value;
  RandomFunction *lhs;
  RandomFunction *rhs;
  size_t unaryFunctionCount;
  size_t binaryFunctionCount;
  bool vectorize;
public:
  RandomFunction(int count, int depth, bool vectorize = true)
  {
    unaryFunctionCount = vectorize ? unaryVectorFunctions.size() : unaryFunctions.size();
    binaryFunctionCount = vectorize ? binaryVectorFunctions.size() : binaryFunctions.size();

    lhs = nullptr;
    rhs = nullptr;
    int r = rand() % 100;
    if (depth < 2)
    {
      if (r < 10)
      {
        value = randfloat2;
        type = TERMINALVALUE;
      }
      else 
      {
        index = rand() % count;
        type = TERMINALINDEX;
      }
    } 
    else 
    {
      if (r < 2)
      {
        index = rand() % count;
        type = TERMINALINDEX;
      }
      else if(r > 50 && depth < 4 )
      {
        do {
          delete lhs;
          lhs = new RandomFunction(count, depth - 1, vectorize);
          type = UNARY;
          index = rand() % unaryFunctionCount;
        } while (lhs->type == TERMINALVALUE);
      }
      else
      {
        do {
          delete lhs;
          delete rhs;
          lhs = new RandomFunction(count, depth - 1, vectorize);
          rhs = new RandomFunction(count, depth - 1, vectorize);
        } while ((*lhs == *rhs) || 
          (lhs->type == TERMINALVALUE && rhs->type == TERMINALVALUE));
        type = BINARY;
        index = rand() % binaryFunctionCount;
      }
    }
  }

  ~RandomFunction()
  {
    delete lhs;
    delete rhs;
  }

  bool operator==(const RandomFunction& other)
  {
    if (type != other.type) return false;

    switch (type)
    {
    case UNARY:
      return index == other.index && *lhs == *other.lhs;
    case BINARY:
      return index == other.index && *lhs == *other.lhs && *rhs == *other.rhs;
    case TERMINALVALUE:
      return value == other.value;
    case TERMINALINDEX:
      return index == other.index;
    }
    return false;
  }

  float Eval(float input)
  {
    vector<float> d;
    d.push_back(input);
    return Eval(d);
  }

  float Eval(float x, float y)
  {
    vector<float> d;
    d.push_back(x);
    d.push_back(y);
    return Eval(d);
  }

  friend ostream& operator<<(ostream& os, const RandomFunction& f)
  {
    switch (f.type)
    {
    case UNARY:
      os << "U" << f.index << "(" << *f.lhs << ")";
      break;
    case BINARY:
      os << "B" << f.index << "(" << *f.lhs << "," << *f.rhs << ")";
      break;
    case TERMINALVALUE:
      os << f.value;
      break;
    case TERMINALINDEX:
      switch (f.index)
      {
      case 0: os << "x"; break;
      case 1: os << "y"; break;
      }
    }
    return os;
  }

  void Eval(int w, int h, float* x, float* y, float* result, bool vectorize)
  {
    if (vectorize)
      EvalVector(w, h, x, y, result);
    else {
      #pragma omp parallel for
      for (size_t j = 0; j < h; ++j)
      {
        #pragma omp parallel for
        for (size_t i = 0; i < w; ++i)
        {
          result[j*w + i] = Eval(x[i], y[j]);
        }
      }
    }
  }

private:
  void EvalVector(int w, int h, float* x, float* y, float* result)
  {
    int count = w * h;
    
    // construct a mesh grid of x and y values
    float *xx = new float[count];
    float *yy = new float[count];
    for (int i = 0; i < h; ++i)
    {
      for (int j = 0; j < w; ++j)
      {
        int pos = w*i + j;
        xx[pos] = x[j];
        yy[pos] = y[i];
      }
    }

    EvalVector(&count, xx, yy, result);

    delete[] xx;
    delete[] yy;
  }

  // assume # of dimensions == 2
  void EvalVector(int* count, float* x, float* y, float* result)
  {
    float *xx = new float[*count];
    float *yy = new float[*count];
    switch (type)
    {
    case UNARY:
      lhs->EvalVector(count, x, y, xx);
      unaryVectorFunctions[index](count, xx, result);
      break;
    case BINARY:
      lhs->EvalVector(count, x, y, xx);
      rhs->EvalVector(count, x, y, yy);
      binaryVectorFunctions[index](count, xx, yy, result);
      break;
    case TERMINALVALUE:
      for (size_t i = 0; i < *count; ++i)
        result[i] = value;
      break;
    case TERMINALINDEX:
      for (size_t i = 0; i < *count; ++i)
        result[i] = index == 0 ? x[i] : y[i];
      break;
    }
    delete[] xx;
    delete[] yy;
  }

  float Eval(vector<float>& params) const
  {
    float result;
    switch (type)
    {
    case UNARY:
      result = unaryFunctions[index](lhs->Eval(params));
      break;
    case BINARY:
      result = binaryFunctions[index](lhs->Eval(params), rhs->Eval(params));
      break;
    case TERMINALVALUE:
      result = value;
      break;
    case TERMINALINDEX:
      result = params[index];
      break;
    }
    return result;
  }
};