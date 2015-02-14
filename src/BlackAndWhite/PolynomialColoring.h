#pragma once

#include "ColoringAlgorithm.h"
#include <boost/math/tools/rational.hpp>

struct PolynomialColoring : public ColoringAlgorithm
{
  uint8_t terms;
  float *r, *g, *b;
  PolynomialColoring(const uint8_t terms) : terms(terms)
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

  virtual ~PolynomialColoring()
  {
    delete[] r;
    delete[] g;
    delete[] b;
  }

  virtual void Value(float input, float &r, float &g, float &b) override
  {
    r += tools::evaluate_polynomial(this->r, input, terms);
    g += tools::evaluate_polynomial(this->g, input, terms);
    b += tools::evaluate_polynomial(this->b, input, terms);
  }
};