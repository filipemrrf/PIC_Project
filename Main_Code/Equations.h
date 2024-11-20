/**
 * @file Equations.h
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief Declares the equations that are to be solved
 * @version 4.0
 * @date 2024-11-14
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef EQUATIONS
#define EQUATIONS

#include <cmath>

#include "Core.h"

/**
 * @brief 2nd order in space and 1st order in time wave equation
 * 
 * @param x Array with the spatial points
 * @param u Array with the values of the fields
 * @param N Number of points in the fields array
 * @param N_Ghosts Number of ghost points in each boundary
 * @param step_x Spatial step
 * @param Acc Accuracy order of the derivatives
 * @param boundary Boundary conditions function
 * @param params Array with the parameters of the equation
 * @param dissipation Dissipation parameter
 */
void Wave_Equation(double* x, double* u, int N, int* N_Ghosts, double step_x, int Acc, BoundaryFunc* boundary, double* params, double dissipation);

/**
 * @brief 1st order in space and 1st order in time wave equation in compactified hyperboloidal coordinates
 * 
 * @param x Array with the spatial points
 * @param u Array with the values of the fields
 * @param N Number of points in the fields array
 * @param N_Ghosts Number of ghost points in each boundary
 * @param step_x Spatial step
 * @param Acc Accuracy order of the derivatives
 * @param boundary Boundary conditions function
 * @param params Array with the parameters of the equation
 * @param dissipation Dissipation parameter
 */
void Compact_Wave_Equation(double* x, double* u, int N, int* N_Ghosts, double step_x, int Acc, BoundaryFunc* boundary, double* params, double dissipation);

/**
 * @brief 2nd order in space and 1st order in time wave equation in spherical coordinates
 * 
 * @param x Array with the spatial points
 * @param u Array with the values of the fields
 * @param N Number of points in the fields array
 * @param N_Ghosts Number of ghost points in each boundary
 * @param step_x Spatial step
 * @param Acc Accuracy order of the derivatives
 * @param boundary Boundary conditions function
 * @param params Array with the parameters of the equation
 * @param dissipation Dissipation parameter
 */
void Spherical_Wave_Equation(double* x, double* u, int N, int* N_Ghosts, double step_x, int Acc, BoundaryFunc* boundary, double* params, double dissipation);

/**
 * @brief 1st order in space and 1st order in time wave equation in spherical coordinates
 * 
 * @param x Array with the spatial points
 * @param u Array with the values of the fields
 * @param N Number of points in the fields array
 * @param N_Ghosts Number of ghost points in each boundary
 * @param step_x Spatial step
 * @param Acc Accuracy order of the derivatives
 * @param boundary Boundary conditions function
 * @param params Array with the parameters of the equation
 * @param dissipation Dissipation parameter
 */
void Spherical_Reduced_Wave_Equation(double* x, double* u, int N, int* N_Ghosts, double step_x, int Acc, BoundaryFunc* boundary, double* params, double dissipation);

/**
 * @brief 1st order in space and 1st order in time wave equation in spherical compactified hyperboloidal coordinates
 * 
 * @param x Array with the spatial points
 * @param u Array with the values of the fields
 * @param N Number of points in the fields array
 * @param N_Ghosts Number of ghost points in each boundary
 * @param step_x Spatial step
 * @param Acc Accuracy order of the derivatives
 * @param boundary Boundary conditions function
 * @param params Array with the parameters of the equation
 * @param dissipation Dissipation parameter
 */
void Spherical_Compact_Wave_Equation(double* x, double* u, int N, int* N_Ghosts, double step_x, int Acc, BoundaryFunc* boundary, double* params, double dissipation);

/**
 * @brief 2nd order in space and 1st order in time non linear wave equation
 * 
 * @param x Array with the spatial points
 * @param u Array with the values of the fields
 * @param N Number of points in the fields array
 * @param N_Ghosts Number of ghost points in each boundary
 * @param step_x Spatial step
 * @param Acc Accuracy order of the derivatives
 * @param boundary Boundary conditions function
 * @param params Array with the parameters of the equation
 * @param dissipation Dissipation parameter
 */
void Non_Linear_Wave_Equation(double* x, double* u, int N, int* N_Ghosts, double step_x, int Acc, BoundaryFunc* boundary, double* params, double dissipation);

/**
 * @brief 2nd order in space and 1st order in time non linear wave equation in spherical coordinates
 * 
 * @param x Array with the spatial points
 * @param u Array with the values of the fields
 * @param N Number of points in the fields array
 * @param N_Ghosts Number of ghost points in each boundary
 * @param step_x Spatial step
 * @param Acc Accuracy order of the derivatives
 * @param boundary Boundary conditions function
 * @param params Array with the parameters of the equation
 * @param dissipation Dissipation parameter
 */
void Non_Linear_Spherical_Wave_Equation(double* x, double* u, int N, int* N_Ghosts, double step_x, int Acc, BoundaryFunc* boundary, double* params, double dissipation);

/**
 * @brief 2nd order in space and 1st order in time ADM evolution equations
 * 
 * @param x Array with the spatial points
 * @param u Array with the values of the fields
 * @param N Number of points in the fields array
 * @param N_Ghosts Number of ghost points in each boundary
 * @param step_x Spatial step
 * @param Acc Accuracy order of the derivatives
 * @param boundary Boundary conditions function
 * @param params Array with the parameters of the equation
 * @param dissipation Dissipation parameter
 */
void ADM_Evolution(double* x, double* u, int N, int* N_Ghosts, double step_x, int Acc, BoundaryFunc* boundary, double* params, double dissipation);

#endif