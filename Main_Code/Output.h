/**
 * @file Output.h
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief Declares the output functions for the main code
 * @version 1.0
 * @date 2023-04-19
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
 * @param params Array containing {N_Eqs, Eq_Select, step_x}, where N_Eqs contains the number of equations in the system, Eq_Select the equation that will be output 
 *                  and step_x the space step
 */
void Output_Solution(std::fstream* FILE, double* u, int N, int N_ghosts, double t, double* params);

void Debug(std::fstream* FILE, double* u, int N, int N_Ghosts, double time, double* params);

#endif