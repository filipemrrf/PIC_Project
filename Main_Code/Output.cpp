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

void Hamiltonian_Constraint(std::fstream* FILE, double* u, int N, int N_Ghosts, double time, double* params){
    // Writes the time step that is being saved
    *FILE << "\"Time = " << time << std::endl;

    // Populates the ghost points
    //Even_0_Boundary(u, N, 9, Acc);

    // Defines pointers to make the manipulation of the state vector easier
    double* A = u;
    double* DA = &(u[N/9]);
    double* KA = &(u[(2*N)/9]);

    double* B = &(u[(3*N)/9]);
    double* DB = &(u[(4*N)/9]);
    double* KB = &(u[(5*N)/9]);

    double* lambda =&(u[(6*N)/9]);

    double* alpha = &(u[(7*N)/9]);
    double* Dalpha = &(u[(8*N)/9]);

    // Allocates memory for the derivative of DB and for the constraint
    double* dr_DB = new double[N/9];
    double* Hamiltonean = new double[N/9];

    // Calculates the derivative of DB
    //if()
        //Second_Derivative_2nd_Order(DB, dr_DB, N/9, params[], 1);
    //else if()
        //Second_Derivative_2nd_Order(DB, dr_DB, N/9, params[], 1);


    // Calculates the Hamiltonean constraint value
    for(int i = 0; i < N/9; ++i)
        //Hamiltonean[i] = -dr_DB[i] - lambda[i]/
    // Saves the value of the Hamiltonian constraint to disk (ignoring ghost points)


    // Frees the memory allocated
    delete[] dr_DB;
    delete[] Hamiltonean;
}

void Debug(std::fstream* FILE, double* u, int N, int N_Ghosts, double time, double* params){
    // Writes the time step that is being saved
    *FILE << "\"Time = " << time << std::endl;

    double* debug = new double[N];
    Even_0_Boundary(u, N, 2, 2);
    Second_Derivative_2nd_Order(u, debug, N, params[2], 2);
    //KO_Dissipation_4th_Order(u, derivative, N, params[2], 2, 0.02);

    // Saves the solution to disk (ignoring ghost points)
    for(int i = 0; i < (N/((int) params[0])); ++i){
        *FILE << (i-N_Ghosts)*params[2] << " " << debug[i /*+ (8*N)/9*/] << std::endl;
    }

    *FILE << std::endl;

    delete[] debug;
}