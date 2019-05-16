cmake_minimum_required(VERSION 3.14)
project(SimplexSolver)

set(CMAKE_CXX_STANDARD 17)

add_executable(SimplexSolver main.cpp simplex_algorithm.cpp simplex_algorithm.h simplex_utils.cpp simplex_utils.h)