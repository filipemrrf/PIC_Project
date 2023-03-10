/**
 * @file Core.h
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief Defines the core algorithms to solve differential equations (declared in Core.h)
 * @version 1.2
 * @date 2023-03-10
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "Core.h"

void Runge_Kutta_4(void (*f)(double* u, int N, double step_x, double* params), double* IC, int N, int N_Eq, double step_x, double* params,
    double tmax, double step_t, std::string filename, int w){

    //Opens the file the solution is to be written to
    std::fstream FILE;
    FILE.open(filename, std::fstream::out);

    //Makes sure the values are printed to full precision
    FILE.precision(std::numeric_limits<double>::max_digits10 + 2);

    //Writes the initial conditions to the file
    FILE << "\"Time = 0.0" << std::endl;
    for(int i = 0; i < N/N_Eq; ++i){
        FILE << step_x*i << " ";

        for(int j = 0; j < N_Eq; ++j)
            FILE << IC[j*N/N_Eq + i] << " ";
        
        FILE << std::endl;
    }

    FILE << std::endl;

    //Allocates memory for the auxiliary arrays
    double* prev_state = new double[N];
    double* temp_state = new double[N];

    //Copies the IC to the previous state array
    for(int i = 0; i < N; ++i)
        prev_state[i] = IC[i];

    //Allocates memory for the 4 slopes in every point
    double* K1 = new double[N];
    double* K2 = new double[N];
    double* K3 = new double[N];
    double* K4 = new double[N];

    //Loops through the time
    for(int t = 1; ((double) t)*step_t <= tmax; ++t){
        //Copies the IC array to the temporary state array
        for(int i = 0; i < N; ++i)
            temp_state[i] = prev_state[i];

        //Calculates the slope K1 for every point
        f(temp_state, N, step_x, params);

        for(int i = 0; i < N; ++i)
            K1[i] = temp_state[i]*step_t;


        //Calculates where the slope K2 is to be calculated at for every point
        for(int i = 0; i < N; ++i)
            temp_state[i] = prev_state[i] + 0.5*K1[i];
        
        //Calculates the slope K2 for every point
        f(temp_state, N, step_x, params);

        for(int i = 0; i < N; ++i)
            K2[i] = temp_state[i]*step_t;

        //Calculates where the slope K3 is to be calculated at for every point
        for(int i = 0; i < N; ++i)
            temp_state[i] = prev_state[i] + 0.5*K2[i];
        
        //Calculates the slope K3 for every point
        f(temp_state, N, step_x, params);

        for(int i = 0; i < N; ++i)
            K3[i] = temp_state[i]*step_t;


        //Calculates where the slope K4 is to be calculated at for every point
        for(int i = 0; i < N; ++i)
            temp_state[i] = prev_state[i] + K3[i];
        
        //Calculates the slope K4 for every point
        f(temp_state, N, step_x, params);

        for(int i = 0; i < N; ++i)
            K4[i] = temp_state[i]*step_t;

        //Calculates the final value of the time evolution
        for(int i = 0; i < N; ++i)
                prev_state[i] = prev_state[i] + (K1[i] + 2.0*K2[i] + 2.0*K3[i] + K4[i])/6.0;

        //Decides if this timestep is to be saved to disk
        if((t%w) == 0){
            //Saves the timestep to disk
            FILE << "\"Time = " << (double)t*step_t << std::endl;

            for(int i = 0; i < N/N_Eq; ++i){
                FILE << step_x*i << " ";

                for(int j = 0; j < N_Eq; ++j)
                    FILE << prev_state[j*N/N_Eq + i] << " ";
                
                FILE << std::endl;
            }

            FILE << std::endl;
        }
    }

    //Closes the file
    FILE.close();

    //Frees the memory allocated for the 4 slopes in every point and the auxiliary vector
    delete[] K1;
    delete[] K2;
    delete[] K3;
    delete[] K4;
    delete[] prev_state;
    delete[] temp_state;
}