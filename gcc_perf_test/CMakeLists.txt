cmake_minimum_required(VERSION 3.2)
project(gcc_perf_test)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

set(SOURCE_FILES main.cpp)

add_executable(gcc_perf_test ${SOURCE_FILES})