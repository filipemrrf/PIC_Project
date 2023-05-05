/**
 * @file Equations.h
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief Declares the equations that are to be solved
 * @version 2.0
 * @date 2023-03-23
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef EQUATIONS
#define EQUATIONS

#include <cmath>

#include "Core.h"

/**
 * @brief Transforms an array using the wave equation (separated into a 2 ODE system)
 * 
 * @param u Array to be transformed by the equation (in the form u = {Phi_0, Phi_1, ..., Phi_n, Pi_0, Pi_1, ..., Pi_n})
 * @param N Number of elements of array u
 * @param Acc Accuracy order of the equation
 * @param boundary Boundary conditions to be imposed
 * @param params Parameters for the equation: params = {c, diss, step_x}, where c is the speed of the wave, diss is the dissipation coeficient and
 *                  step_X is the spatial step       
 */
void Wave_Equation(double* u, int N, int Acc, BoundaryFunc* boundary, double* params);

/**
 * @brief Transforms an array using the non linear wave equation (separated into a 2 ODE system) with a non linear coefficient n
 * 
 * @param u Array to be transformed by the equation (in the form u = {Phi_0, Phi_1, ..., Phi_n, Pi_0, Pi_1, ..., Pi_n})
 * @param N Number of elements of array u
 * @param step_x Space step
 * @param params Parameters for the equation: params = {c, n, diss, step_x}, where c is the speed of the wave, n is the power of the non linear coeficient of the wave,
 *                  diss is the dissipation coeddicient and step_x is the spatial step
 */
void Non_Linear_Wave_Equation(double* u, int N, int Acc, BoundaryFunc* boundary, double* params);

/**
 * @brief Transforms an array using the wave equation (separated into a 2 ODE system) in spherical symmetry
 * 
 * @param u Array to be transformed by the equation (in the form u = {Phi_0, Phi_1, ..., Phi_n, Pi_0, Pi_1, ..., Pi_n})
 * @param N Number of elements of array u
 * @param step_x Space step
 * @param params Parameters for the equation: params = {c}, where c is the speed of the wave
 */
void Spherical_Wave_Equation(double* u, int N, int Acc, BoundaryFunc* boundary, double* params);

#endif