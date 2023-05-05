cmake_minimum_required(VERSION 3.16.3)

project(PIC_Project)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O1")

add_executable(main Main_Code/main.cpp Main_Code/Core.cpp Main_Code/Equations.cpp Main_Code/Output.cpp)