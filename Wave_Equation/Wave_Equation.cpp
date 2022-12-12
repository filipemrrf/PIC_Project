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

void Runge_Kutta_4(void (*f)(double* u, int N, double step_x, void* params), double* IC, int N, int N_Eq, double step_x, void* params, double tmax, double step_t, std::string filename){
    //Opens the file the solution is to be written to
    std::fstream FILE;
    FILE.open(filename, std::fstream::out);

    //Writes the initial conditions to the file
    FILE << "\"Time = 0.0" << std::endl;
    for(int i = 0; i < N/N_Eq; ++i){
        FILE << step_x*i << " ";

        for(int j = 0; j < N_Eq; ++j)
            FILE << IC[j*N/N_Eq + i] << " ";
        
        FILE << std::endl;
    }

    FILE << std::endl;

    //Copies the initial conditions to an auxiliary array
    double* aux = new double[N];

    for(int i = 0; i < N; ++i)
        aux[i] = IC[i];

    //Allocates memory for the 4 slopes in every point
    double* K1 = new double[N];
    double* K2 = new double[N];
    double* K3 = new double[N];
    double* K4 = new double[N];

    //Loops through the time
    for(int t = 0; ((double) t)*step_t <= tmax; ++t){
        //Copies the IC array to K1
        for(int i = 0; i < N; ++i)
            K1[i] = aux[i];

        //Calculates the slope K1 for every point
        f(K1, N, step_x, params);

        for(int i = 0; i < N; ++i)
            K1[i] *= step_t;


        //Calculates where the slope K2 is to be calculated at for every point
        for(int i = 0; i < N; ++i)
            K2[i] = aux[i] + 0.5*K1[i];
        
        //Calculates the slope K2 for every point
        f(K2, N, step_x, params);

        for(int i = 0; i < N; ++i)
            K2[i] *= step_t;

        //Calculates where the slope K3 is to be calculated at for every point
        for(int i = 0; i < N; ++i)
            K3[i] = aux[i] + 0.5*K2[i];
        
        //Calculates the slope K3 for every point
        f(K3, N, step_x, params);

        for(int i = 0; i < N; ++i)
            K3[i] *= step_t;


        //Calculates where the slope K4 is to be calculated at for every point
        for(int i = 0; i < N; ++i)
            K4[i] = aux[i] + K3[i];
        
        //Calculates the slope K4 for every point
        f(K4, N, step_x, params);

        for(int i = 0; i < N; ++i)
            K4[i] *= step_t;

        //Saves the solution to the ODE system
        FILE << "\"Time = " << (double)t*step_t << std::endl;

        for(int i = 0; i < N; ++i)
            aux[i] = aux[i] + (K1[i] + 2.0*K2[i] + 2.0*K3[i] + K4[i])/6.0;

        for(int i = 0; i < N/N_Eq; ++i){
            FILE << step_x*i << " ";

            for(int j = 0; j < N_Eq; ++j)
                FILE << aux[j*N/N_Eq + i] << " ";
            
            FILE << std::endl;
        }

        FILE << std::endl;
    }

    //Closes the file
    FILE.close();

    //Frees the memory allocated for the 4 slopes in every point and the auxiliary vector
    delete[] K1;
    delete[] K2;
    delete[] K3;
    delete[] K4;
    delete[] aux;
}

void Wave_Equation(double* u, int N, double step_x, void* params){
    //Allocates memory for the transformed array
    double* Du = new double[N];

    //Declarates auxiliary pointers for easier readability of the function
    double* Phi0 = u;
    double* Pi0 = &(u[N/2]);

    double* Phi1 = Du;
    double* Pi1 = &(Du[N/2]);

    //Calculates the auxiliary value k = (c/step_x)^2
    double k = *((double*) params)/step_x;
    k *= k;

    //Transforms the array
    for(int i = 0; i < N/2; ++i){
        Phi1[i] = Pi0[i];

        if((i == 0) || (i == (N/2 - 1)))
            Pi1[i] = k*(Phi0[1] - 2*Phi0[0] + Phi0[N/2-2]);
        else
            Pi1[i] = k*(Phi0[i+1] - 2*Phi0[i] + Phi0[i-1]);
    }

    //Copies the transformed array to the original one
    for(int i = 0; i < N; ++i)
        u[i] = Du[i];

    //Frees the memory allocated for the transformed array
    delete[] Du;
}