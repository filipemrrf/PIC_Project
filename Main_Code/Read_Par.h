/**
 * @file Read_Par.h
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief Delcares the functions that read the parameter file and store the values in the given pointers
 * @version 1.0
 * @date 2024-11-14
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef READ_PAR
#define READ_PAR

#include <sstream>
#include <string>
#include <string.h>

#include "Core.h"
#include "Equations.h"
#include "Output.h"

/**
 * @brief Reads the parameter file and stores the values in the given pointers
 * 
 * @param parfilename Name of the parameter file
 * @param rh_side Pointer to the equation to be solved
 * @param Acc Pointer to the accuracy order
 * @param boundary Pointer to the boundary conditions
 * @param params_rh_side Pointer to the parameters for the equation
 * @param tmax Pointer to the time until the system is to be solved
 * @param x Pointer to the space array
 * @param IC Pointer to the initial conditions for the different fields
 * @param N_IC Pointer to the size of the initial conditions array
 * @param N_Ghosts Pointer to thearray with the number of ghost points
 * @param step_x Pointer to the space step
 * @param output Pointer to the array of output functions
 * @param N_output Pointer to the size of the output functions array
 * @param params_output Pointer to the parameters for the output functions
 * @param out_filenames Pointer to the array of filenames where the output will be saved
 * @param write_con Pointer to the writting condition for the output
 * @param cfl Pointer to the CFL coefficient
 * @param diss Pointer to the dissipation coefficient
 */
void read_par_file(std::string parfilename, rh_sideFunc** rh_side, int* Acc, BoundaryFunc** boundary, double** params_rh_side, double* tmax, double ** x, double** IC,
    int* N_IC, int** N_Ghosts, double* step_x, OutputFunc*** output, int* N_output, double*** params_output, std::string** out_filenames, 
    int* write_con, double* cfl, double* diss);

/**
 * @brief Reads the equation to be solved and returns the pointer to the function
 * 
 * @param eq_name String with the name of the equation to be solved
 * @param N_Vars Number of fields of the equation
 * @return Pointer to the function that calculates the right-hand side of the equation
 */
rh_sideFunc* read_eq(std::string eq_name, int* N_Vars);

/**
 * @brief Reads the accuracy to be used and returns the integer value
 * 
 * @param acc String with the name of the accuracy to be used
 * @return Integer value of the accuracy to be used
 */
int read_acc(std::string acc);

/**
 * @brief Reads the boundary conditions to be used and returns the pointer to the function
 * 
 * @param bound String with the name of the boundary conditions to be used
 * @return Pointer to the function that calculates the boundary conditions
 */
BoundaryFunc* read_bound(std::string bound);

/**
 * @brief Reads the output functions to be used and returns the pointer to the function
 * 
 * @param data String with the name of the output function to be used
 * @param output Pointer to the array of output functions
 * @param params_output Pointer to the array of parameters for the output functions
 * @param out_filenames Pointer to the array of filenames where the output will be saved
 * @param N_Vars Number of fields of the equation
 * @param index Current index of the output array
 * @param step_x Space step
 * @param dissipation Dissipation coefficient
 */
void read_output(std::string data, OutputFunc*** output, double*** params_output, std::string** out_filenames, double N_Vars, int* index, double step_x, double dissipation);

#endif