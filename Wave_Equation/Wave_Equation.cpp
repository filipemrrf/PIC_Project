/**
 * @file Wave_Equation.cpp
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief Wave_Equation class definition
 * @version 0.3
 * @date 2022-11-30
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "Wave_Equation.h"

/**
 * @brief Construct a new Wave_Equation::Wave_Equation object
 * 
 * @param X0 initial conditions
 * @param step_X step for the positions
 * @param step_T step for the time
 * @param C wave propagation speed
 */
Wave_Equation::Wave_Equation(const std::vector<std::vector<double>> &X0, double step_X, double step_T, double C){
    x0 = X0;

    step_x = step_X;
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
Wave_Equation::Wave_Equation(const std::function<double(double)> &X0, const std::function<double(double)> &XL0, double xmax, double step_X, double step_T, double C){
    for(int i = 0; (double)i*step_X < xmax; ++i){
        x0[0].push_back(X0((double)i*step_X));
        x0[1].push_back(XL0((double)i*step_X));
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
void Wave_Equation::Set_IC(const std::vector<std::vector<double>> &X0, double step_X){
    x0 = X0;

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
void Wave_Equation::Set_IC(std::function<double(double)> X0, std::function<double(double)> XL0, double xmax, double step_X){
    for(int i = 0; (double)i*step_X < xmax; ++i){
        x0[0].push_back(X0((double)i*step_X));
        x0[1].push_back(XL0((double)i*step_X));
    }

    step_x = step_X;
}

/**
 * @brief Sets the time-step on a preexisting Wave_Equation::Wave_Equation object
 * 
 * @param step_T 
 */
void Wave_Equation::Set_Step_t(double step_T){
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
 * @brief Solves the wave function with the parameters previously defined up to time t=tmax, using Runge-Kutta 4
 * 
 * @param tmax final time of the solution
 */
void Wave_Equation::Solve(double tmax){
    //Calculates a constant that will be used multiple times throughout the loops
    double k = c/step_x;
    k *= k*step_t;

    //Declares the variables to store the 4 slopes of the rk4 method
    double K1, K2, K3, K4;

    //Loops through the time
    for(int t = 0; (double)t*step_t < tmax; ++t){
        //Creates auxiliary vectors to store the solutions
        std::vector<std::vector<double>> aux;

        //Loops through the positions
        for(int x = 0; x < x0[0].size()-1; ++x){
                //1st ODE
                //Calculates the 4 slopes necessary for the rk4 method

                //Calculates the 1st slope
                K1 = step_t*x0[1][x];


                //Calculates the index of the point where the 2nd slope will be calculated
                int xK2 = x + K1/(2.0*step_x);

                while(xK2 > (x0[1].size() - 2))
                    xK2 -= (x0[1].size() + 2);

                //Calculates the 2nd slope
                K2 = step_t*(x0[1][xK2] + 0.5*step_t*x0[1][xK2]);


                //Calculates the index of the point where the 3rd slope will be calculated
                int xK3 = x + K2/(2.0*step_x);

                while(xK3 > (x0[1].size() - 2))
                    xK3 -= (x0[1].size() + 2);
                    
                //Calculates the 3rd slope
                K3 = step_t*(x0[1][xK3] + 0.5*step_t*x0[1][xK3]);


                //Calculates the index of the point where the 4th slope will be calculated
                int xK4 = x + K3/step_x;

                while(xK4 > (x0[1].size() - 2))
                    xK4 -= (x0[1].size() + 2);

                //Calculates the 4th slope
                K4 = k*(x0[1][xK4] + 0.5*step_t*x0[1][xK4]);


                //Saves the solution to the solution vector
                aux[0].push_back(x0[1][x] + (K1 + 2*K2 + 2* K3 + K4)/6.0);


                //2nd ODE
                //Calculates the 4 slopes necessary for the rk4 method

                //Checks if we are calculating for the first point and calculates the first slope accordingly
                if(x == 0)
                    K1 = k*(x0[0][1] - 2*x0[0][0] + x0[0][x0[0].size()-2]);
                else
                    K1 = k*(x0[0][x+1] - 2*x0[0][x] + x0[0][x-1]);


                //Calculates the index of the point where the 2nd slope will be calculated
                xK2 = x + K1/(2.0*step_x);

                while(xK2 > (x0[0].size() - 2))
                    xK2 -= (x0[0].size() + 2);

                //Calculates the value of the function in the point where the slope will be calculated
                std::vector<double> yK2;

                for(int i = -1; i <= 1; ++i){
                    if(xK2+i == 0)
                        yK2.push_back(x0[0][0] + 0.5*k*(x0[0][1] - 2*x0[0][0] + x0[0][x0[0].size()-2]));
                    else
                        yK2.push_back(x0[0][xK2+i] + 0.5*k*(x0[0][xK2+i+1] - 2*x0[0][xK2+i] + x0[0][xK2+i-1]));
                }

                //Calculates the 2nd slope
                K2 = k*(yK2[2] - 2*yK2[1] + yK2[0]);


                //Calculates the index of the point where the 3rd slope will be calculated
                xK3 = x + K2/(2.0*step_x);

                while(xK3 > (x0[0].size() - 2))
                    xK3 -= (x0[0].size() + 2);

                //Calculates the value of the function in the point where the slope will be calculated
                std::vector<double> yK3;

                for(int i = -1; i <= 1; ++i){
                    if(xK3+i == 0)
                        yK3.push_back(x0[0][0] + 0.5*k*(x0[0][1] - 2*x0[0][0] + x0[0][x0[0].size()-2]));
                    else
                        yK3.push_back(x0[0][xK3+i] + 0.5*k*(x0[0][xK3+i+1] - 2*x0[0][xK3+i] + x0[0][xK3+i-1]));
                }

                //Calculates the 3rd slope
                K3 = k*(yK3[2] - 2*yK3[1] + yK3[0]);


                //Calculates the index of the point where the 4th slope will be calculated
                xK4 = x + K3/step_x;

                while(xK4 > (x0[0].size() - 2))
                    xK4 -= (x0[0].size() + 2);

                //Calculates the value of the function in the point where the slope will be calculated
                std::vector<double> yK4;

                for(int i = -1; i <= 1; ++i){
                    if(xK4+i == 0)
                        yK4.push_back(x0[0][0] + k*(x0[0][1] - 2*x0[0][0] + x0[0][x0[0].size()-2]));
                    else
                        yK4.push_back(x0[0][xK4+i] + k*(x0[0][xK4+i+1] - 2*x0[0][xK4+i] + x0[0][xK4+i-1]));
                }

                //Calculates the 4th slope
                K4 = k*(yK4[2] - 2*yK4[1] + yK4[0]);


                //Saves the solution to the solution vector
                aux[1].push_back(x0[0][x] + (K1 + 2*K2 + 2* K3 + K4)/6.0);

        }

        //Mantains the symmetry by copying the first point
        aux[0].push_back(aux[0][0]);
        aux[1].push_back(aux[1][0]);

        //Saves the solution to the solution vector
        Sol.push_back(aux);

        //clears the aux vector (not needed but just to make sure nothing goes wrong)
        aux.clear();
    }
}


void Wave_Equation::Save_Sol(std::string filename){
    //Opens the file
    std::fstream F(filename);

    //Writtes the header of the file
    F << "\"Output of X axis, variable ana.CSI, rank 0, grid 0, type CU, orientation: CE" << std::endl;

    for(int i = 0; i < Sol.size(); ++i){
        //Writtes the time in the file
        F << "\" Time = " << step_t*i << std::endl;
        for(int j = 0; j < Sol[0][0].size(); ++j)
            //Writtes the position and the function value
            F << j*step_x << " " << Sol[i][0][j] << std::endl;

        F << std::endl;
    }

    //Closes the file
    F.close();
}