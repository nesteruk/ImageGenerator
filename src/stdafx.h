#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdint>
#include <ctime>
#include <functional>
#include <vector>
#include <random>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <chrono>
using namespace std;
using namespace chrono;

// M_PI... not on GCC, apparently
#ifndef M_PI
constexpr auto M_PI = 3.14159265358979323846;
#endif

#include "lodepng.h"

#include <boost/math/special_functions.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/bind.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
//#include <boost/mpi.hpp>
#include <boost/optional.hpp>
#include <boost/serialization/optional.hpp>
using namespace boost::accumulators;
using namespace boost::math;

#ifdef USEMKL
#include <mkl.h>
#endif
#include <omp.h>
//#include <mpi.h>