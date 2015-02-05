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
public:
  RandomFunction(int count, int depth)
  {
    lhs = nullptr;
    rhs = nullptr;
    int r = rand() % 100;
    if (depth < 3)
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
      else if(r > 50 /*&& depth < 4*/)
      {
        do {
          delete lhs;
          lhs = new RandomFunction(count, depth - 1);
          type = UNARY;
          index = rand() % unaryFunctions.size();
        } while (lhs->type == TERMINALVALUE);
      }
      else
      {
        do {
          delete lhs;
          delete rhs;
          lhs = new RandomFunction(count, depth - 1);
          rhs = new RandomFunction(count, depth - 1);
        } while ((*lhs == *rhs) || 
          (lhs->type == TERMINALVALUE && rhs->type == TERMINALVALUE));
        type = BINARY;
        index = rand() % binaryFunctions.size();
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
      os << "U(" << *f.lhs << ")";
      break;
    case BINARY:
      os << "B(" << *f.lhs << "," << *f.rhs << ")";
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

  void Eval(int *count, float* x, float* y, float* result)
  {
    float *xx = new float[*count];
    float *yy = new float[*count];
    switch (type)
    {
    case UNARY:
      lhs->Eval(count, x, y, xx);
      unaryVectorFunctions[index](count, xx, result);
      break;
    case BINARY:
      lhs->Eval(count, x, y, xx);
      rhs->Eval(count, x, y, yy);
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