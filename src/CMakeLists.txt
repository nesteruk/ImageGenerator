cmake_minimum_required(VERSION 3.2)
project(src)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")
set(Boost_INCLUDE_DIR /cygdrive/c/boost)
set(Boost_LIBRARY_DIR /cygdrive/c/boost/stage/lib)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -O3")

find_package(OpenMP)
find_package(Boost 1.56 COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

set(SOURCE_FILES
    stdafx.cpp
    imagegen.cpp
    lodepng.cpp
    RandomFunction.hpp
)

add_executable(ImageGenCygwin ${SOURCE_FILES})
target_link_libraries(ImageGenCygwin ${Boost_LIBRARIES})