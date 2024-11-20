/**
 * @file Output.h
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief Declares the output functions for the main code
 * @version 3.0
 * @date 2024-11-14
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef OUTPUT
#define OUTPUT

#include <fstream>

#include "Core.h"

/**
 * @brief Outputs the solution of the equation
 * 
 * @param FILE Pointer to the file to be written
 * @param x Array with the spatial points
 * @param u Array with the solution of the equation
 * @param N Size of the u array x
 * @param N_Ghosts Number of ghost points on each side of the array
 * @param t Current time
 * @param params Additional parameters for the output function
 */
void Output_Solution(std::fstream* FILE, double* x, double* u, int N, int* N_Ghosts, double t, double* params);

/**
 * @brief Outputs the right hand side of the equation
 * 
 * @param FILE Pointer to the file to be written
 * @param x Array with the spatial points
 * @param u Array with the solution of the equation
 * @param N Size of the u array x
 * @param N_Ghosts Number of ghost points on each side of the array
 * @param t Current time
 * @param params Additional parameters for the output function
 */
void Output_RHS(std::fstream* FILE, double* x, double* u, int N, int* N_Ghosts, double t, double* params);

/**
 * @brief Outputs the dissipation operator of the equation
 * 
 * @param FILE Pointer to the file to be written
 * @param x Array with the spatial points
 * @param u Array with the solution of the equation
 * @param N Size of the u array x
 * @param N_Ghosts Number of ghost points on each side of the array
 * @param t Current time
 * @param params Additional parameters for the output function
 */
void Output_Dissipation(std::fstream* FILE, double* x, double* u, int N, int* N_Ghosts, double t, double* params);

/**
 * @brief Outputs the first derivative of the solution of the equation
 * 
 * @param FILE Pointer to the file to be written
 * @param x Array with the spatial points
 * @param u Array with the solution of the equation
 * @param N Size of the u array x
 * @param N_Ghosts Number of ghost points on each side of the array
 * @param t Current time
 * @param params Additional parameters for the output function
 */
void Output_1st_Derivative(std::fstream* FILE, double* x, double* u, int N, int* N_Ghosts, double t, double* params);

/**
 * @brief Outputs the second derivative of the solution of the equation
 * 
 * @param FILE Pointer to the file to be written
 * @param x Array with the spatial points
 * @param u Array with the solution of the equation
 * @param N Size of the u array x
 * @param N_Ghosts Number of ghost points on each side of the array
 * @param t Current time
 * @param params Additional parameters for the output function
 */
void Output_2nd_Derivative(std::fstream* FILE, double* x, double* u, int N, int* N_Ghosts, double t, double* params);

/**
 * @brief Outputs the reduction constraint of the solution of the equation
 * 
 * @param FILE Pointer to the file to be written
 * @param x Array with the spatial points
 * @param u Array with the solution of the equation
 * @param N Size of the u array x
 * @param N_Ghosts Number of ghost points on each side of the array
 * @param t Current time
 * @param params Additional parameters for the output function
 */
void Output_Constraint(std::fstream* FILE, double* x, double* u, int N, int* N_Ghosts, double t, double* params);


//void Reduction_Constraints_ADM(std::fstream* FILE, double* x, double* u, int N, int* N_Ghosts, double t, double* params);
//void Hamiltonian_Constraint_ADM(std::fstream* FILE, double* x, double* u, int N, int* N_Ghosts, double t, double* params);
//void Momentum_Constraint_ADM(std::fstream* FILE, double* x, double* u, int N, int* N_Ghosts, double t, double* params);

#endif