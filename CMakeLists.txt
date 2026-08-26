cmake_minimum_required(VERSION 3.10)
project(EULER_SOLVER_1D)

# must use C++17 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# multiple cpp source files
set(SOURCES ./src/EulerSolver.cpp ./src/main.cpp)
# now generate executable "euler_solver.exe"
add_executable(euler_solver ${SOURCES})