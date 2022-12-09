cmake_minimum_required(VERSION 3.16.3)

project(PIC_Project)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O1")

add_executable(test.exe main.cpp Wave_Equation/Wave_Equation.cpp)

