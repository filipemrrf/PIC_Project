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
#include <vector>

#include <gsl/gsl_const_mksa.h>

/**
 * @brief Solves the wave equation with the given parameters
 * 
 * @param x0 Initial consitions
 * @param step_x Space step
 * @param tmax Time the equation is to be solved until
 * @param c Speed of the wave
 * @param step_t Time step
 * @return std::vector<std::vector<std::vector<double>>> Solution of the equation (for each time, returns a vector with the function value and another with its derivative)
 */
std::vector<std::vector<std::vector<double>>> Solve_Wave_Eq(std::vector<std::vector<double>> x0, double step_x, double tmax, double c = {GSL_CONST_MKSA_SPEED_OF_LIGHT}, double step_t = 1e-3);

/**
 * @brief Saves the solution to the equation in a way muninn can display
 * 
 * @param Sol Solution to the equation
 * @param step_x Space step
 * @param filename Name of the file that will be saved
 * @param step_t Time step
 */
void Save_Sol(std::vector<std::vector<std::vector<double>>> Sol, double step_x, std::string filename, double step_t = 1e-3);

#endif