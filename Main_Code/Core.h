/**
 * @file Core.h
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief Declares the core algorithms to solve differential equations
 * @version 3.0
 * @date 2023-03-23
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef CORE
#define CORE

#include <fstream>
#include <iomanip>
#include <limits>
#include <string>

#include <iostream>

/* Defines the function types that will be referenced */
typedef void (BoundaryFunc)(double* u, int N, int N_Var, int Acc);
typedef void (OutputFunc)(std::fstream* FILE, double* u, int N, int N_ghosts, double t, double* params);
typedef void (rh_sideFunc)(double* u, int N, int Acc, BoundaryFunc* boundary, double* params);

/* Declares numerical derivatives */
/**
 * @brief Calculates the 1st numerical derivative of an array with 2nd order accuracy
 * 
 * @param u Array to calculate the derivative of
 * @param du Array to save the derivative to
 * @param N Number of points in the array
 * @param step_x Space step 
 * @param N_Vars Number of variables in the array
 */
void First_Derivative_2nd_Order(double* u, double* du, int N, double step_x, int N_Vars);

/**
 * @brief Calculates the 1st numerical derivative of an array with 4th order accuracy
 * 
 * @param u Array to calculate the derivative of
 * @param du Array to save the derivative to
 * @param N Number of points in the array
 * @param step_x Space step 
 * @param N_Vars Number of variables in the array
 */
void First_Derivative_4th_Order(double* u, double* du, int N, double step_x, int N_Vars);

/**
 * @brief Calculates the 2nd numerical derivative of an array with 2nd order accuracy
 * 
 * @param u Array to calculate the derivative of
 * @param du Array to save the derivative to
 * @param N Number of points in the array
 * @param step_x Space step 
 * @param N_Vars Number of variables in the array
 */
void Second_Derivative_2nd_Order(double* u, double* du, int N, double step_x, int N_Vars);

/**
 * @brief Calculates the 2nd numerical derivative of an array with 4th order accuracy
 * 
 * @param u Array to calculate the derivative of
 * @param du Array to save the derivative to
 * @param N Number of points in the array
 * @param step_x Space step 
 * @param N_Vars Number of variables in the array
 */
void Second_Derivative_4th_Order(double* u, double* du, int N, double step_x, int N_Vars);


/* Declares dissipation operators */
/**
 * @brief Calculates the Kreiss-Oliger dissipation coeficient of the array compatible with a 2nd order accurate code
 * 
 * @param u Array to calculate the KO dissipation of
 * @param du Array to save the KO dissipation to
 * @param N Size of the array
 * @param step_x Space step
 * @param N_Vars Number of variables in the array
 * @param dissipation Dissipation coefficient
 */
void KO_Dissipation_4th_Order(double* u, double* du, int N, double step_x, int N_Vars, double dissipation);

/**
 * @brief Calculates the Kreiss-Oliger dissipation coeficient of the array compatible with a 4th order accurate code
 * 
 * @param u Array to calculate the KO dissipation of
 * @param du Array to save the KO dissipation to
 * @param N Size of the array
 * @param step_x Space step
 * @param N_Vars Number of variables in the array
 * @param dissipation Dissipation coefficient
 */
void KO_Dissipation_6th_Order(double* u, double* du, int N, double step_x, int N_Vars, double dissipation);


/* Declares boundary conditions */

void Even_Constant_Boundary(double* u, int N, int N_Var, int Acc);

void Odd_Constant_Boundary(double* u, int N, int N_Var, int Acc);

/**
 * @brief Populates ghost points using parity of the IC in the left and with 0 on the right
 * 
 * @param u Array whose ghost points will be populated
 * @param N Size of the array
 * @param N_Var Number of variables in the array
 * @param Acc Accuracy order being used
 */
void Even_0_Boundary(double* u, int N, int N_Var, int Acc);

/**
 * @brief Populates ghost points applying periodicity to the IC
 * 
 * @param u Array whose ghost points will be populated
 * @param N Size of the array
 * @param N_Var Number of variables in the array
 * @param Acc Accuracy order being used
 */
void Periodic_Boundary(double* u, int N, int N_Var, int Acc);


/* Declares equation solver methods */
/**
 * @brief Solves a system of ODE's using the Runge-Kutta 4 method
 * 
 * @param rh_side System of equations to be solved
 * @param IC Intial conditions
 * @param N_IC Size of the initial conditions array
 * @param params_rh_side Parameters for the system of equations
 * @param output Array ou output functions that will be called
 * @param N_output Size of the output function array
 * @param params_output Parameters for the output functions
 * @param filename Array with the name of the output functions
 * @param write_con Array containing information on whether or not to output
 * @param tmax Time until the the system is to be solved
 * @param step_t Time step
 */
void Runge_Kutta_4(rh_sideFunc* rh_side, double* IC, int N_IC, int Acc, BoundaryFunc* boundary, double* params_rh_side, 
    OutputFunc** output, int N_output, double** params_output, std::string* out_filenames, int* write_con, double tmax, double step_t);

#endif