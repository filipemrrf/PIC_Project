/**
 * @file Core.h
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief Declares the core algorithms to solve differential equations
 * @version 4.0
 * @date 2024-11-17
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef CORE
#define CORE

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>

#include <iostream>

/* Defines the function types that will be referenced */
typedef void (BoundaryFunc)(double* u, int N, int N_Var, int* N_Ghosts);
typedef void (OutputFunc)(std::fstream* FILE, double* x, double* u, int N, int* N_Ghosts, double t, double* params);
typedef void (rh_sideFunc)(double* x, double* u, int N, int* N_Ghosts, double step_x, int Acc, BoundaryFunc* boundary, double* params, double dissipation);

/* Declares difference operators */
/**
 * @brief Calculates the first derivative of a function using a 2nd order central difference scheme
 * 
 * @param u Array with the function values to be differentiated
 * @param du Array where the derivative will be stored
 * @param N Number of points in the array
 * @param step_x Step size between points
 * @param N_Ghosts Number of ghost points at the boundaries
 * @param N_Vars Number of fields in the array
 */
void First_Derivative_2nd_Order(double* u, double* du, int N, double step_x, int* N_Ghosts, int N_Vars);

/**
 * @brief Calculates the first derivative of a function using a 4th order central difference scheme
 * 
 * @param u Array with the function values to be differentiated
 * @param du Array where the derivative will be stored
 * @param N Number of points in the array
 * @param step_x Step size between points
 * @param N_Ghosts Number of ghost points at the boundaries
 * @param N_Vars Number of fields in the array
 */
void First_Derivative_4th_Order(double* u, double* du, int N, double step_x, int* N_Ghosts, int N_Vars);

/**
 * @brief Calculates the second derivative of a function using a 2nd order central difference scheme
 * 
 * @param u Array with the function values to be differentiated
 * @param du Array where the derivative will be stored
 * @param N Number of points in the array
 * @param step_x Step size between points
 * @param N_Ghosts Number of ghost points at the boundaries
 * @param N_Vars Number of fields in the array
 */
void Second_Derivative_2nd_Order(double* u, double* du, int N, double step_x, int* N_Ghosts, int N_Vars);

/**
 * @brief Calculates the second derivative of a function using a 4th order central difference scheme
 * 
 * @param u Array with the function values to be differentiated
 * @param du Array where the derivative will be stored
 * @param N Number of points in the array
 * @param step_x Step size between points
 * @param N_Ghosts Number of ghost points at the boundaries
 * @param N_Vars Number of fields in the array
 */
void Second_Derivative_4th_Order(double* u, double* du, int N, double step_x, int* N_Ghosts, int N_Vars);

/**
 * @brief Applies the Evans method to calculate the derivative of a function
 * 
 * @param u Array with the function values to be differentiated
 * @param du Array where the derivative will be stored
 * @param N Number of points in the array
 * @param step_x Step size between points
 * @param N_Ghosts Number of ghost points at the boundaries
 * @param N_Vars Number of fields in the array
 * @param p Exponent of the Evans method
 */
void Evans_Method(double* x, double* u, double* du, int N, double step_x, int* N_Ghosts, int N_Vars, double p);


/* Declares dissipation operators */
/**
 * @brief Calculates the Kreiss-Oliger dissipation of a function using a 4th order scheme
 * 
 * @param u Array with the values of the fields
 * @param du Array where the dissipation will be stored
 * @param N Number of points in the array
 * @param step_x Step size between points
 * @param N_Ghosts Number of ghost points at the boundaries
 * @param N_Vars Number of fields in the array
 * @param dissipation Dissipation coefficient
 */
void KO_Dissipation_4th_Order(double* u, double* du, int N, double step_x, int* N_Ghosts, int N_Vars, double dissipation);

/**
 * @brief Calculates the Kreiss-Oliger dissipation of a function using a 6th order scheme
 * 
 * @param u Array with the values of the fields
 * @param du Array where the dissipation will be stored
 * @param N Number of points in the array
 * @param step_x Step size between points
 * @param N_Ghosts Number of ghost points at the boundaries
 * @param N_Vars Number of fields in the array
 * @param dissipation Dissipation coefficient
 */
void KO_Dissipation_6th_Order(double* u, double* du, int N, double step_x, int* N_Ghosts, int N_Vars, double dissipation);


/* Declares boundary conditions */
/**
 * @brief Applies the even boundary condition to the left  boundary and constant boundary condition to the right boundary
 * 
 * @param u Array with the values of the fields
 * @param N Size of the array
 * @param N_Var Number of fields in the array
 * @param N_Ghosts Number of ghost points at each boundary
 */
void Even_Constant_Boundary(double* u, int N, int N_Var, int* N_Ghosts);

/**
 * @brief Applies the odd boundary condition to the left  boundary and constant boundary condition to the right boundary
 * 
 * @param u Array with the values of the fields
 * @param N Size of the array
 * @param N_Var Number of fields in the array
 * @param N_Ghosts Number of ghost points at each boundary
 */
void Odd_Constant_Boundary(double* u, int N, int N_Var, int* N_Ghosts);

/**
 * @brief Applies the even boundary condition to the left boundary and sets the ghost points on the right boundary to 0
 * 
 * @param u Array with the values of the fields
 * @param N Size of the array
 * @param N_Var Number of fields in the array
 * @param N_Ghosts Number of ghost points at each boundary
 */
void Even_0_Boundary(double* u, int N, int N_Var, int* N_Ghosts);

/**
 * @brief Applies the even boundary condition to the left boundary and does extrapolation to the right boundary
 * 
 * @param u Array with the values of the fields
 * @param N Size of the array
 * @param N_Var Number of fields in the array
 * @param N_Ghosts Number of ghost points at each boundary
 */
void Even_Extrapolation_Boundary(double* u, int N, int N_Var, int* N_Ghosts);

/**
 * @brief Applies the odd boundary condition to the left boundary and does extrapollation for the right boundary
 * 
 * @param u Array with the values of the fields
 * @param N Size of the array
 * @param N_Var Number of fields in the array
 * @param N_Ghosts Number of ghost points at each boundary
 */
void Odd_Extrapolation_Boundary(double* u, int N, int N_Var, int* N_Ghosts);

/**
 * @brief Applies constant boundary conditions to both boundaries
 * 
 * @param u Array with the values of the fields
 * @param N Size of the array
 * @param N_Var Number of fields in the array
 * @param N_Ghosts Number of ghost points at each boundary
 */
void Constant_Constant_Boundary(double* u, int N, int N_Var, int* N_Ghosts);

/**
 * @brief Applies periodic boundary conditions to the boundaries
 * 
 * @param u Array with the values of the fields
 * @param N Size of the array
 * @param N_Var Number of fields in the array
 * @param N_Ghosts Number of ghost points at each boundary
 */
void Periodic_Boundary(double* u, int N, int N_Var, int* N_Ghosts);

/**
 * @brief Applies poison to the boundaries
 * 
 * @param u Array with the values of the fields
 * @param N Size of the array
 * @param N_Var Number of fields in the array
 * @param N_Ghosts Number of ghost points at each boundary
 */
void Poison(double* u, int N, int N_Var, int* N_Ghosts);


/* Declares integration methods */
/**
 * @brief 
 * 
 * @param rh_side Pointer to the function that calculates the right-hand side of the differential equation
 * @param x Spatial grid
 * @param IC Initial conditions array
 * @param N_IC Size of the initial conditions array
 * @param N_Ghosts Number of ghost points at each boundary
 * @param step_x Spatial step size
 * @param Acc Accuracy of the difference operators
 * @param boundary Pointer to the function that applies the boundary conditions
 * @param params_rh_side Array with the parameters of the right-hand side function
 * @param output Pointer to the array of functions that write the output
 * @param N_output Size of the output array
 * @param params_output Array with the parameters of the output functions
 * @param out_filenames Array with the filenames of the output files
 * @param write_con Write condition
 * @param tmax Total time of the simulation
 * @param step_t Time step size
 * @param dissipation Dissipation coefficient
 */
void Runge_Kutta_4(rh_sideFunc* rh_side, double* x, double* IC, int N_IC, int* N_Ghosts, double step_x, int Acc, BoundaryFunc* boundary, double* params_rh_side, 
    OutputFunc** output, int N_output, double** params_output, std::string* out_filenames, int write_con, double tmax, double step_t, double dissipation);

#endif