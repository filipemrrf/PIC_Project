/**
 * @file ODE_System.cpp
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief ODE_System class definition
 * @version 1.0
 * @date 2022-11-30
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "ODE_System.h"

/**
 * @brief Construct a new ode system solver::ode system solver object
 * 
 * @param X0 initial conditions
 * @param step_X step for space
 * @param step_T step for time
 */
ODE_System::ODE_System(const std::vector<std::vector<double>> &X0, double step_X, double step_T){
    x0 = X0;
    step_x = step_X;
    step_t = step_T;
}

/**
 * @brief Set the initial conditions on a preexisting ODE_System_Solver object
 * 
 * @param X0 initial conditions
 * @param step_X step for the positions 
 */
void ODE_System::Set_IC(const std::vector<std::vector<double>> &X0, double step_X){
    x0 = X0;
    step_x = step_X;
}

/**
 * @brief Set the time step on a preexisting ODE_System_Solver object
 * 
 * @param step_T 
 */
void ODE_System::Set_Step_T(double step_T){
    step_t = step_T;
}

/**
 * @brief Returns the solution to the system of equations
 * 
 * @return std::vector<std::vector<double>> solutions of the ODE system
 */
std::vector<std::vector<std::vector<double>>> ODE_System::Get_Sol(){
    return Sol;
}

/**
 * @brief Destroy the ode system solver object
 * 
 */
ODE_System::~ODE_System(){
    x0.clear();
    Sol.clear();
}