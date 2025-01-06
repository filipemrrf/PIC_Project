/**
 * @file main.cpp
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief Reads the arguments from the command line and solves the specified equations with given parameters and initial conditions
 * @version 4.0
 * @date 2024-11-14
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include <fstream>
#include <iostream>
#include <string>
#include <string.h>

#include "Core.h"
#include "Equations.h"
#include "Output.h"
#include "Read_Par.h"

int main(int argc, char** argv){
    // Declares the string to store the parameter file name
    std::string INPUT;

    // Checks if the command line arguments are valid
    if(argc == 1){
        std::cout << "ERROR: No parameter file given" << std::endl;
        exit(-1);
    }
    if(argc > 2){
        std::cout << "ERROR: Too many arguments given" << std::endl;
        exit(-1);
    }

    // Reads the name of the parameter file
    INPUT = argv[1];

    // Declares the variables that will be read from the parameter file
    rh_sideFunc* rh_side = nullptr; // Equation to be solved
    int Acc = 0; // Accuracy order
    BoundaryFunc* boundary = nullptr; // Boundary conditions
    double* params_rh_side = nullptr; // Parameters for the equation

    double tmax = 0.0; // Time until the system is to be solved

    double* x = nullptr; // Space array
    double* IC = nullptr; // Initial conditions for the different fields
    int N_IC = 0; // Size of the initial conditions array
    int* N_Ghosts = nullptr; // Number of ghost points
    double step_x = 0.0; // Space step

    OutputFunc** output = nullptr; // Array of output functions
    int N_output = 0; // Size of the output functions array
    double** params_output = nullptr; // Parameters for the output functions
    std::string* out_filenames = nullptr; // Array of filenames where the output will be saved
    int write_con = 1; // Writting conditions for the output

    double cfl = 0.25; // CFL coefficient
    double diss = 0.0; // Dissipation coefficient


    // Reads the parameter file
    read_par_file(INPUT, &rh_side, &Acc, &boundary, &params_rh_side, &tmax, &x, &IC, &N_IC, &N_Ghosts, &step_x, &output, &N_output, &params_output, &out_filenames, 
    &write_con, &cfl, &diss);


    // Solves the equation specifies with the given parameters
    Runge_Kutta_4(rh_side, x, IC, N_IC, N_Ghosts, step_x, Acc, boundary, params_rh_side, output, N_output, params_output, out_filenames, write_con, tmax, 
    step_x*cfl, diss);

    // Deletes the memory allocated for the initial conditions
    delete[] params_rh_side;
    delete[] x;
    delete[] IC;
    delete[] N_Ghosts;
    for (int i = 0; i < N_output; ++i)
        delete[] params_output[i];
    delete[] output;
    delete[] params_output;
    delete[] out_filenames;
    
    return 0;
}