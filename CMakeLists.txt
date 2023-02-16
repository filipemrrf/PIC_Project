cmake_minimum_required(VERSION 3.16.3)

project(PIC_Project)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O1")

add_executable(main main.cpp Equations.cpp Core.cpp)

add_executable(test test_drive.cpp)