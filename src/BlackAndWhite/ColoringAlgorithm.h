#pragma once

struct ColoringAlgorithm
{
  virtual void Value(float input, float &r, float &g, float &b) = 0;
};

struct NoColoring : ColoringAlgorithm
{
  NoColoring(int _) {}
  virtual void Value(float input, float &r, float &g, float &b) override 
  {
    r = g = b = input;
  }
};