#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <string>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstdlib>
using namespace std;

#include <boost/math/special_functions.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/bind.hpp>
#include <boost/tuple/tuple.hpp>
using namespace boost::accumulators;
using namespace boost::math;

#include <tbb/parallel_for_each.h>
using namespace tbb;

#include <ipp.h>
#include <mkl.h>