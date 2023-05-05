/**
 * @file Output.h
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief Defines the output functions for the main code
 * @version 1.0
 * @date 2023-04-19
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#include "Output.h"

void Output_Solution(std::fstream* FILE, double* u, int N, int N_Ghosts, double time, double* params){
    // Writes the time stepthat is being saved
    *FILE << "\"Time = " << time << std::endl;

    // Saves the solution to disk (ignoring ghost points)
    for(int i = 0; i < (N/((int) params[0]) - 2*N_Ghosts); ++i)
        *FILE << i*params[2] << " " << u[i + ((int) params[1])*(N/((int) params[0])) + N_Ghosts] << std::endl;

    *FILE << std::endl;
}

void Debug(std::fstream* FILE, double* u, int N, int N_Ghosts, double time, double* params){
    // Writes the time step that is being saved
    *FILE << "\"Time = " << time << std::endl;

    double* debug = new double[N];
    Even_0_Boundary(u, N, 2, 4);
    //Second_Derivative_4th_Order(u, debug, N, params[2], 2);
    //KO_Dissipation_4th_Order(u, derivative, N, params[2], 2, 0.02);

    // Saves the solution to disk (ignoring ghost points)
    for(int i = 0; i < (N/((int) params[0]) /*- 2*N_Ghosts*/); ++i){
        *FILE << (i-N_Ghosts)*params[2] << " " << debug[i + N/2] << std::endl;
    }

    *FILE << std::endl;

    delete[] debug;
}