/**
 * @file Core.h
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief Declares the core algorithms to solve differential equations
 * @version 1.2
 * @date 2023-03-10
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef CORE
#define CORE

#include <fstream>
#include <iomanip>
#include <limits>
#include <string>

/**
 * @brief Solves a system of ODE's using the Runge-Kutta 4 method
 * 
 * @param f System of equations to be solved
 * @param IC Intial conditions
 * @param N Number of elements of the initial conditions array
 * @param N_Eq Number of equations in the ODE system
 * @param step_x Space step
 * @param params Parameters for the system of equations
 * @param tmax Time until the the system is to be solved
 * @param step_t Time step
 * @param filename Name of the file which will store the data
 * @param w Multiples of the iterations that will be saved to disk
 */
void Runge_Kutta_4(void (*f)(double* u, int N, double step_x, double* params), double* IC, int N, int N_Eq, double step_x, double* params,
    double tmax, double step_t, std::string filename, int w);



#endif