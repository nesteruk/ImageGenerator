#include <vector>
#include <cmath>
#include <random>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <chrono>
using namespace std;
using namespace chrono;

#include "functions.h"
#include "RandomFunction.hpp"

int main()
{
    RandomFunction rf{2, 8, false};

    float z = 0;

    const int w = 1024, h  = 1024;

    float* x = new float[w];
    float* y = new float[h];
    float* result = new float[w*h];

    auto now = high_resolution_clock::now();
    rf.Eval(w,h,x,y,result, false);
    auto end = high_resolution_clock::now();

    cout << result[rand() % (w*h)] << endl;
    cout << "time taken " << duration_cast<milliseconds>(end - now).count() << "msec" << endl;

    delete[]x;
    delete[]y;
    delete[]result;
}