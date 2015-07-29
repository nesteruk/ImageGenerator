#pragma once

#include "ColoringAlgorithm.h"
#include <boost/math/tools/rational.hpp>

struct PolynomialColoring : ColoringAlgorithm
{
  uint8_t terms;
  float *r, *g, *b;

  explicit PolynomialColoring(const uint8_t terms) : terms(terms)
  {
    r = new float[terms];
    g = new float[terms];
    b = new float[terms];
    for (size_t i = 0; i < terms; ++i)
    {
      r[i] = randfloat;
      g[i] = randfloat;
      b[i] = randfloat;
    }
  }

  

  virtual void Value(const float input, float &r, float &g, float &b) const override
  {
    r += tools::evaluate_polynomial(this->r, input, terms);
    g += tools::evaluate_polynomial(this->g, input, terms);
    b += tools::evaluate_polynomial(this->b, input, terms);
  }

  virtual ~PolynomialColoring()
  {
    delete[] r;
    delete[] g;
    delete[] b;
  }
};