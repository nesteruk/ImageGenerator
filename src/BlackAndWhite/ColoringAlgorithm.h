#pragma once

struct ColoringAlgorithm
{
  virtual void Value(float input, float &r, float &g, float &b) = 0;
};