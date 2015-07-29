#pragma once

struct ColoringAlgorithm
{
protected:
  ~ColoringAlgorithm()  {  }

public:
  virtual void Value(const float input, float &r, float &g, float &b) const = 0;
};

struct NoColoring : ColoringAlgorithm
{
protected:
  ~NoColoring()  {  }

public:
  explicit NoColoring(int _) {}
  virtual void Value(const float input, float &r, float &g, float &b) const override 
  {
    r = g = b = input;
  }
};