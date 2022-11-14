/**
 * @file Wave_Equation.cpp
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief Wave_Equation class definition
 * @version 0.1
 * @date 2022-11-14
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "Wave_Equation.h"

/**
 * @brief Construct a new Wave_Equation::Wave_Equation object
 * 
 * @param X0 initial conditions for the function
 * @param XL0 initial coditions for the function's derivative
 * @param step_X step for the positions
 * @param step_T step for the time
 * @param C wave propagation speed
 */
Wave_Equation::Wave_Equation(std::vector<double> X0, std::vector<double> XL0, double step_X, double step_T = 1e-3, double C = 299792458){
    x0 = X0;
    xl0 = XL0;
    step_x =step_X;
    step_t = step_T;
    c = C;
}

/**
 * @brief Construct a new Wave_Equation::Wave_Equation object
 * 
 * @param X0 initial conditions for the function
 * @param XL0 initial coditions for the function's derivative
 * @param xmax length of the domain
 * @param step_X step for the positions
 * @param step_T step for the time
 * @param C wave propagation speed 
 */
Wave_Equation::Wave_Equation(std::function<double(double)> X0, std::function<double(double)> XL0, double xmax, double step_X = 1e-3, double step_T = 1e-3, double C = 299792458){
    for(int i = 0; (double)i*step_X < xmax; ++i){
        x0.push_back(X0((double)i*step_X));
        xl0.push_back(XL0((double)i*step_X));
    }

    step_x = step_X;
    step_t = step_T;
    c = C;
}

/**
 * @brief Set the initial conditions on a preexisting Wave_Equation::Wave_Equation object
 * 
 * @param X0 initial conditions for the function
 * @param XL0 initial coditions for the function's derivative
 * @param step_X step for the positions 
 */
void Wave_Equation::Set_IC(std::vector<double> X0, std::vector<double> XL0, double step_X){
    x0 = X0;
    xl0 = XL0;
    step_x = step_X;
}

/**
 * @brief Set the initial conditions on a preexisting Wave_Equation::Wave_Equation object
 * 
 * @param X0 initial conditions for the function
 * @param XL0 initial coditions for the function's derivative
 * @param xmax length of the domain
 * @param step_X step for the positions
 */
void Wave_Equation::Set_IC(std::function<double(double)> X0, std::function<double(double)> XL0, double xmax, double step_X = 1e-3){
    for(int i = 0; (double)i*step_X < xmax; ++i){
        x0.push_back(X0((double)i*step_X));
        xl0.push_back(XL0((double)i*step_X));
    }

    step_x = step_X;
}

/**
 * @brief Set the time step on a preexisting Wave_Equation::Wave_Equation object
 * 
 * @param step_T new time step
 */
void Wave_Equation::Set_Step(double step_T){
    step_t = step_T;
}

/**
 * @brief Set the wave propagation speed on a preexisting Wave_Equation::Wave_Equation object
 * 
 * @param C new wave propagation speed
 */
void Wave_Equation::Set_C(double C){
    c = C;
}

/**
 * @brief Solves the wave function with the parameters previously defined up to time t=tmax
 * 
 * @param tmax final time of the solution
 */
void Wave_Equation::Solve(double tmax){
    //initialization of a auxiliary variables
    std::vector<double> aux; //aux vector
    double K = c*step_t/step_x;
    K = K*K;

    //adds t=0 to the solution
    sol.push_back(x0);

    //adds t=delta t to the solution
    for(int i = 0; i < x0.size(); ++i)
        aux.push_back(x0[i]+xl0[i]*step_t);
    
    sol.push_back(aux);
    aux.clear();

    //loops through the time
    for(int i = 2; (double)i*step_t < tmax; ++i){
        //loops through the positions
        for(int j = 0; j < x0.size(); ++j)
            aux.push_back(2.0*sol[i-1][j]-sol[i-2][j]+K*(sol[i-1][j+1]-2.0*sol[i-1][j]+sol[i-1][j-1])); //Does not work! The borders are not accounted for

        sol.push_back(aux);
        aux.clear();
    }
}

/**
 * @brief Destroy the Wave_Equation::Wave_Equation object
 * 
 */
Wave_Equation::~Wave_Equation(){
    x0.clear();
    xl0.clear();
    sol.clear();
}