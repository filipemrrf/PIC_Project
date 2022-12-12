/**
 * @file Wave_Equation.h
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief Wave_Equation solver and output saver declaration
 * @version 0.4
 * @date 2022-12-3
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef WAVE_EQUATION
#define WAVE_EQUATION

#include <fstream>
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
 */
void Runge_Kutta_4(void (*f)(double* u, int N, double step_x, void* params), double* IC, int N, int N_Eq, double step_x, void* params, double tmax, double step_t, std::string filename);

/**
 * @brief Transforms an array using the wave equation (separated into a 2 ODE system)
 * 
 * @param u Array to be transformed by the equation (in the form u = {Phi_0, Phi_1, ..., Phi_n, Pi_0, Pi_1, ..., Pi_n})
 * @param N Number of elements of array u
 * @param step_x Space step
 * @param params Parameters for the equation: params = {c}, where c is the speed of the wave
 */
void Wave_Equation(double* u, int N, double step_x, void* params);

#endif