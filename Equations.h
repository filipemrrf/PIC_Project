/**
 * @file Equations.h
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief Declares the equations that are to be solved
 * @version 1.2
 * @date 2022-12-28
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef EQUATIONS
#define EQUATIONS

#include <cmath>

/**
 * @brief Transforms an array using the wave equation (separated into a 2 ODE system)
 * 
 * @param u Array to be transformed by the equation (in the form u = {Phi_0, Phi_1, ..., Phi_n, Pi_0, Pi_1, ..., Pi_n})
 * @param N Number of elements of array u
 * @param step_x Space step
 * @param params Parameters for the equation: params = {c}, where c is the speed of the wave
 */
void Wave_Equation(double* u, int N, double step_x, double* params);

/**
 * @brief Transforms an array using the wave equation (separated into a 2 ODE system)
 * 
 * @param u Array to be transformed by the equation (in the form u = {Phi_0, Phi_1, ..., Phi_n, Pi_0, Pi_1, ..., Pi_n})
 * @param N Number of elements of array u
 * @param step_x Space step
 * @param params Parameters for the equation: params = {c, d}, where c is the speed of the wave and d is the dissipation coefficient
 */
void Wave_Equation_Dissipation(double* u, int N, double step_x, double* params);

/**
 * @brief Transforms an array using the 4th order accurate wave equation (separated into a 2 ODE system)
 * 
 * @param u Array to be transformed by the equation (in the form u = {Phi_0, Phi_1, ..., Phi_n, Pi_0, Pi_1, ..., Pi_n})
 * @param N Number of elements of array u
 * @param step_x Space step
 * @param params Parameters for the equation: params = {c}, where c is the speed of the wave
 */
void Wave_Equation_4th_Order(double* u, int N, double step_x, double* params);

/**
 * @brief Transforms an array using the non linear wave equation (separated into a 2 ODE system) with a non linear coefficient n
 * 
 * @param u Array to be transformed by the equation (in the form u = {Phi_0, Phi_1, ..., Phi_n, Pi_0, Pi_1, ..., Pi_n})
 * @param N Number of elements of array u
 * @param step_x Space step
 * @param params Parameters for the equation: params = {c, n}, where c is the speed of the wave and n is the power of the non linear coeficient of the wave
 */
void Non_Linear_Wave_Equation(double* u, int N, double step_x, double* params);

/**
 * @brief Transforms an array using the non linear wave equation (separated into a 2 ODE system) with a non linear coefficient n and artificial (Kreiss-Oliger) dissipation
 * 
 * @param u Array to be transformed by the equation (in the form u = {Phi_0, Phi_1, ..., Phi_n, Pi_0, Pi_1, ..., Pi_n})
 * @param N Number of elements of array u
 * @param step_x Space step
 * @param params Parameters for the equation: params = {c, n, d}, where c is the speed of the wave, n is the power of the non linear coeficient of the wave and d is the dissipation coefficient
 */
void Non_Linear_Wave_Equation_Dissipation(double* u, int N, double step_x, double* params);

/**
 * @brief Transforms an array using the wave equation (separated into a 2 ODE system) in spherical symmetry
 * 
 * @param u Array to be transformed by the equation (in the form u = {Phi_0, Phi_1, ..., Phi_n, Pi_0, Pi_1, ..., Pi_n})
 * @param N Number of elements of array u
 * @param step_x Space step
 * @param params Parameters for the equation: params = {c}, where c is the speed of the wave
 */
void Spherical_Wave_Equation(double* u, int N, double step_x, double* params);


void Spherical_Wave_Equation_Dissipation(double* u, int N, double step_x, double* params);

#endif