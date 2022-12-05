/**
 * @file Wave_Equation.h
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief Wave_Equation solver and output saver definition
 * @version 0.4
 * @date 2022-12-3
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "Wave_Equation.h"

std::vector<std::vector<std::vector<double>>> Solve_Wave_Eq(std::vector<std::vector<double>> x0, double step_x, double tmax, double c, double step_t){
    //Calculates a constant that will be used multiple times throughout the loops
    double k = c/step_x;
    k *= k*step_t;

    //Declares the variables to store the 4 slopes of the rk4 method
    double K1, K2, K3, K4;

    //Declares the solution vector and adds the initial conditions to the solution vector
    std::vector<std::vector<std::vector<double>>> Sol = {x0};

    //Loops through the time
    for(int t = 0; (double)t*step_t < tmax; ++t){
        //Creates auxiliary vectors to store the solutions
        std::vector<std::vector<double>> aux;
        std::vector<double> aux_0, aux_1;

        //Loops through the positions
        for(int x = 0; x < x0[0].size()-1; ++x){
                //1st ODE
                //Calculates the 4 slopes necessary for the rk4 method

                //Calculates the 1st slope
                K1 = step_t*x0[1][x];


                //Calculates the index of the point where the 2nd slope will be calculated
                int xK2 = x + K1/(2.0*step_x);

                if(xK2 < 0){
                    xK2 %= (x0[1].size() - 1);
                    xK2 += (x0[1].size() - 1);
                }

                if(xK2 > (x0[1].size() - 1))
                    xK2 %= (x0[1].size() - 1);


                /*if(xK2 == (x0[1].size() - 1))
                    xK2 = 0;*/

                //Calculates the 2nd slope
                K2 = step_t*(x0[1][xK2] + 0.5*step_t*x0[1][xK2]);


                //Calculates the index of the point where the 3rd slope will be calculated
                int xK3 = x + K2/(2.0*step_x);

                if(xK3 < 0){
                    xK3 %= (x0[1].size() - 1);
                    xK3 += (x0[1].size() - 1);
                }

                if(xK3 > (x0[1].size() - 1))
                    xK3 %= (x0[1].size() - 1);
                    
                //Calculates the 3rd slope
                K3 = step_t*(x0[1][xK3] + 0.5*step_t*x0[1][xK3]);


                //Calculates the index of the point where the 4th slope will be calculated
                int xK4 = x + K3/*/step_x*/;

                if(xK4 < 0){
                    xK4 %= (x0[1].size() - 1);
                    xK4 += (x0[1].size() - 1);
                }

                if(xK4 > (x0[1].size() - 1))
                    xK4 %= (x0[1].size() - 1);

                //Calculates the 4th slope
                K4 = k*(x0[1][xK4] + 0.5*step_t*x0[1][xK4]);


                //Saves the solution to the solution vector
                aux_0.push_back(x0[1][x] + (K1 + 2*K2 + 2* K3 + K4)/6.0);


                //2nd ODE
                //Calculates the 4 slopes necessary for the rk4 method

                //Checks if we are calculating for the first point and calculates the first slope accordingly
                if(x == 0)
                    K1 = k*(x0[0][1] - 2*x0[0][0] + x0[0][x0[0].size()-2]);
                else
                    K1 = k*(x0[0][x+1] - 2*x0[0][x] + x0[0][x-1]);


                //Calculates the index of the point where the 2nd slope will be calculated
                xK2 = x + K1/(2.0*step_x);

                if(xK2 < 0){
                    xK2 %= (x0[0].size() - 1);
                    xK2 += (x0[0].size() - 1);
                }

                if(xK2 > (x0[0].size() - 1))
                    xK2 %= (x0[0].size() - 1);

                //Calculates the value of the function in the point where the slope will be calculated
                std::vector<double> yK2;

                if(xK2 == 0){
                    yK2.push_back(x0[0][x0[0].size()-1] + 0.5*k*(x0[0][x0[0].size()-2] - 2*x0[0][0] + x0[0][x0[0].size()-3]));
                    yK2.push_back(x0[0][0] + 0.5*k*(x0[0][1] - 2*x0[0][0] + x0[0][x0[0].size()-2]));
                    yK2.push_back(x0[0][2] + 0.5*k*(x0[0][1] - 2*x0[0][1] + x0[0][0]));
                }
                else
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

                if(xK3 < 0){
                    xK3 %= (x0[0].size() - 1);
                    xK3 += (x0[0].size() - 1);
                }

                if(xK3 > (x0[0].size() - 1))
                    xK3 %= (x0[0].size() - 1);

                //Calculates the value of the function in the point where the slope will be calculated
                std::vector<double> yK3;

                if(xK3 == 0){
                    yK3.push_back(x0[0][x0[0].size()-1] + 0.5*k*(x0[0][x0[0].size()-2] - 2*x0[0][0] + x0[0][x0[0].size()-3]));
                    yK3.push_back(x0[0][0] + 0.5*k*(x0[0][1] - 2*x0[0][0] + x0[0][x0[0].size()-2]));
                    yK3.push_back(x0[0][2] + 0.5*k*(x0[0][1] - 2*x0[0][1] + x0[0][0]));
                }
                else
                    for(int i = -1; i <= 1; ++i){
                        if(xK3+i == 0)
                            yK3.push_back(x0[0][0] + 0.5*k*(x0[0][1] - 2*x0[0][0] + x0[0][x0[0].size()-2]));
                        else
                            yK3.push_back(x0[0][xK3+i] + 0.5*k*(x0[0][xK3+i+1] - 2*x0[0][xK3+i] + x0[0][xK3+i-1]));
                    }

                //Calculates the 3rd slope
                K3 = k*(yK3[2] - 2*yK3[1] + yK3[0]);


                //Calculates the index of the point where the 4th slope will be calculated
                xK4 = x + K3/*/step_x*/;

                if(xK4 < 0){
                    xK4 %= (x0[0].size() - 1);
                    xK4 += (x0[0].size() - 1);
                }

                if(xK4 > (x0[0].size() - 1))
                    xK4 %= (x0[0].size() - 1);

                //Calculates the value of the function in the point where the slope will be calculated
                std::vector<double> yK4;

                if(xK4 == 0){
                    yK4.push_back(x0[0][x0[0].size()-1] + 0.5*k*(x0[0][x0[0].size()-2] - 2*x0[0][0] + x0[0][x0[0].size()-3]));
                    yK4.push_back(x0[0][0] + 0.5*k*(x0[0][1] - 2*x0[0][0] + x0[0][x0[0].size()-2]));
                    yK4.push_back(x0[0][2] + 0.5*k*(x0[0][1] - 2*x0[0][1] + x0[0][0]));
                }
                else
                    for(int i = -1; i <= 1; ++i){
                        if(xK4+i == 0)
                            yK4.push_back(x0[0][0] + k*(x0[0][1] - 2*x0[0][0] + x0[0][x0[0].size()-2]));
                        else
                            yK4.push_back(x0[0][xK4+i] + k*(x0[0][xK4+i+1] - 2*x0[0][xK4+i] + x0[0][xK4+i-1]));
                    }

                //Calculates the 4th slope
                K4 = k*(yK4[2] - 2*yK4[1] + yK4[0]);


                //Saves the solution to the solution vector
                aux_1.push_back(x0[0][x] + (K1 + 2*K2 + 2* K3 + K4)/6.0);

        }

        //Mantains the symmetry by copying the first point
        aux_0.push_back(aux_0[0]);
        aux_1.push_back(aux_1[0]);

        //Saves the solution to the solution vector
        aux.push_back(aux_0);
        aux.push_back(aux_1);
        Sol.push_back(aux);

        //clears the aux vector (not needed but just to make sure nothing goes wrong)
        aux_0.clear();
        aux_1.clear();
        aux.clear();
    }
    
    return Sol;
}

void Save_Sol(std::vector<std::vector<std::vector<double>>> Sol, double step_x, std::string filename, double step_t){
    //Opens the file
    std::fstream F;
    F.open(filename, std::fstream::out);

    //Writtes the header of the file
    F << "\"Output of X axis, variable ana.CSI, rank 0, grid 0, type CU, orientation: CE" << std::endl;

    for(int i = 0; i < Sol.size(); ++i){
        //Writtes the time in the file
        F << "\"Time = " << step_t*i << std::endl;
        for(int j = 0; j < Sol[0][0].size(); ++j)
            //Writtes the position and the function value
            F << j*step_x << " " << Sol[i][0][j] << std::endl;

        F << std::endl;
    }

    F << "\"Continue from checkpoint" << std::endl;
    F << "\"Continue from checkpoint" << std::endl;

    F.close();
}