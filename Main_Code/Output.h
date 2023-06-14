/**
 * @file Output.h
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief Declares the output functions for the main code
 * @version 2.0
 * @date 2023-06-13
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef OUTPUT
#define OUTPUT

#include <fstream>

#include "Core.h"

/**
 * @brief Writes the solution to the ODE system for a given time step
 * 
 * @param FILE (Pointer to the) Pointer to the file that will be written
 * @param u Array with the solution that will be written
 * @param N Size of the solution array
 * @param N_Ghosts Number of ghost points on each side of de array
 * @param t Time the solution corresponds to
 * @param params Array containing {N_Fields, Field_Select, step_x}, where N_N_Fields contains the number of fields to be solved, Field_Select the field that will 
 *                  be output and step_x is the space step
 */
void Output_Solution(std::fstream* FILE, double* u, int N, int N_ghosts, double time, double* params);

/**
 * @brief Writes the reduction constraint value at a given point in time
 * 
 * @param FILE (Pointer to the) Pointer to the file that will be written
 * @param u Array with the solution that will be written
 * @param N Size of the solution array
 * @param N_Ghosts Number of ghost points on each side of de array
 * @param t Time the solution corresponds to
 * @param params Array containing {step_x, Acc, Field_Select}, where step_x is the space step, Acc the accuracy the code is being run at, and Field_Select the field 
 *                  that will be output
 */
void Reduction_Constraints(std::fstream* FILE, double* u, int N, int N_Ghosts, double time, double* params);

/**
 * @brief Writes the hamiltonian constraint value at a given point in time
 * 
 * @param FILE (Pointer to the) Pointer to the file that will be written
 * @param u Array with the solution that will be written
 * @param N Size of the solution array
 * @param N_Ghosts Number of ghost points on each side of de array
 * @param t Time the solution corresponds to
 * @param params Array containing {step_x, Acc}, where step_x is the space step, and Acc the accuracy the code is being run at
 */
void Hamiltonian_Constraint(std::fstream* FILE, double* u, int N, int N_Ghosts, double time, double* params);

/**
 * @brief Writes the momentum constraint value at a given point in time
 * 
 * @param FILE (Pointer to the) Pointer to the file that will be written
 * @param u Array with the solution that will be written
 * @param N Size of the solution array
 * @param N_Ghosts Number of ghost points on each side of de array
 * @param t Time the solution corresponds to
 * @param params Array containing {step_x, Acc}, where step_x is the space step, and Acc the accuracy the code is being run at
 */
void Momentum_Constraint(std::fstream* FILE, double* u, int N, int N_Ghosts, double time, double* params);

/**
 * @brief Output function used to help debug the code
 * 
 * @param FILE (Pointer to the) Pointer to the file that will be written
 * @param u Array with the solution that will be written
 * @param N Size of the solution array
 * @param N_Ghosts Number of ghost points on each side of de array
 * @param t Time the solution corresponds to
 * @param params Array containing {}
 */
void Debug(std::fstream* FILE, double* u, int N, int N_Ghosts, double time, double* params);

#endif