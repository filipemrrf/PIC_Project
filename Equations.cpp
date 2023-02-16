/**
 * @file Equations.h
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief Defines the equations that are to be solved (declared in Equations.h)
 * @version 1.2
 * @date 2022-12-28
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "Equations.h"

void Wave_Equation(double* u, int N, double step_x, double* params){
    //Allocates memory for the transformed array
    double* Du = new double[N];

    //Declarates auxiliary pointers for easier readability of the function
    double* Phi0 = u;
    double* Pi0 = &(u[N/2]);

    double* Phi1 = Du;
    double* Pi1 = &(Du[N/2]);

    //Calculates the auxiliary value k = (c/step_x)^2
    double k = *(params)/step_x;
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

void Wave_Equation_Dissipation(double* u, int N, double step_x, double* params){
    //Allocates memory for the transformed array
    double* Du = new double[N];

    //Declarates auxiliary pointers for easier readability of the function
    double* Phi0 = u;
    double* Pi0 = &(u[N/2]);

    double* Phi1 = Du;
    double* Pi1 = &(Du[N/2]);

    //Calculates the auxiliary value k = (c/step_x)^2
    double k = *(params)/step_x;
    k *= k;

    //Transforms the array
    for(int i = 0; i < N/2; ++i){
        Phi1[i] = Pi0[i];

        if((i == 0) || (i == (N/2-1))){
            Pi1[i] = k*(Phi0[1] - 2*Phi0[0] + Phi0[N/2-1]);
            Pi1[i] -= params[1]*(Phi0[0] - Phi0[1])*(Phi0[0] - Phi0[1])*(Phi0[0] - Phi0[N/2-1])*(Phi0[0] - Phi0[N/2-1])/(16.0*step_x);
        }
        else{
            Pi1[i] = k*(Phi0[i+1] - 2*Phi0[i] + Phi0[i-1]);
            Pi1[i] -= params[1]*(Phi0[i] - Phi0[i+1])*(Phi0[i] - Phi0[i+1])*(Phi0[i] - Phi0[i-1])*(Phi0[i] - Phi0[i-1])/(16.0*step_x);
        }
    }

    //Copies the transformed array to the original one
    for(int i = 0; i < N; ++i)
        u[i] = Du[i];

    //Frees the memory allocated for the transformed array
    delete[] Du;
}

void Wave_Equation_4th_Order(double* u, int N, double step_x, double* params){
    //Allocates memory for the transformed array
    double* Du = new double[N];

    //Declarates auxiliary pointers for easier readability of the function
    double* Phi0 = u;
    double* Pi0 = &(u[N/2]);

    double* Phi1 = Du;
    double* Pi1 = &(Du[N/2]);

    //Calculates the auxiliary value k = c^2/(12.0*(step_x)^2)
    double k = *(params)/step_x;
    k *= k/12.0;

    //Transforms the array
    for(int i = 0; i < N/2; ++i){
        Phi1[i] = Pi0[i];

        if((i == 0) || (i == (N/2 - 1)))
            Pi1[i] = k*(-Phi0[2] + 16.0*Phi0[1] - 30.0*Phi0[0] + 16.0*Phi0[N/2-2] - Phi0[N/2-3]);
        if(i == 1)
            Pi1[1] = k*(-Phi0[3] + 16.0*Phi0[2] - 30.0*Phi0[1] + 16.0*Phi0[0] - Phi0[N/2-2]);
        if(i == (N/2 - 2))
            Pi1[N/2-2] = k*(-Phi0[1] + 16.0*Phi0[0] - 30.0*Phi0[N/2-2] + 16.0*Phi0[N/2-3] - Phi0[N/2-4]);
        if((i != 0) && (i != 1) && (i != (N/2-2)) && (i != N/2-1))
            Pi1[i] = k*(-Phi0[i+2] + 16.0*Phi0[i+1] - 30.0*Phi0[i] + 16.0*Phi0[i-1] - Phi0[i-2]);
    }

    //Copies the transformed array to the original one
    for(int i = 0; i < N; ++i)
        u[i] = Du[i];

    //Frees the memory allocated for the transformed array
    delete[] Du;
}

void Non_Linear_Wave_Equation(double* u, int N, double step_x, double* params){
    //Allocates memory for the transformed array
    double* Du = new double[N];

    //Declarates auxiliary pointers for easier readability of the function
    double* Phi0 = u;
    double* Pi0 = &(u[N/2]);

    double* Phi1 = Du;
    double* Pi1 = &(Du[N/2]);

    //Calculates the auxiliary value k = (c/step_x)^2
    double k = (params[0])/step_x;
    k *= k;

    //Transforms the array
    for(int i = 0; i < N/2; ++i){
        Phi1[i] = Pi0[i];

        if((i == 0) || (i == (N/2 - 1)))
            Pi1[i] = k*(Phi0[1] - 2*Phi0[0] + Phi0[N/2-2]) + pow(Phi0[0], params[1]);
        else
            Pi1[i] = k*(Phi0[i+1] - 2*Phi0[i] + Phi0[i-1]) + pow(Phi0[i], params[1]);
    }

    //Copies the transformed array to the original one
    for(int i = 0; i < N; ++i)
        u[i] = Du[i];

    //Frees the memory allocated for the transformed array
    delete[] Du;
}

void Non_Linear_Wave_Equation_Dissipation(double* u, int N, double step_x, double* params){
    //Allocates memory for the transformed array
    double* Du = new double[N];

    //Declarates auxiliary pointers for easier readability of the function
    double* Phi0 = u;
    double* Pi0 = &(u[N/2]);

    double* Phi1 = Du;
    double* Pi1 = &(Du[N/2]);

    //Calculates the auxiliary value k = (c/step_x)^2
    double k = (params[0])/step_x;
    k *= k;

    //Transforms the array
    for(int i = 0; i < N/2; ++i){
        Phi1[i] = Pi0[i];

        if((i == 0) || (i == (N/2 - 1))){
            Pi1[i] = k*(Phi0[1] - 2*Phi0[0] + Phi0[N/2-2]) + pow(Phi0[0], params[1]);
            Pi1[i] -= (Phi0[0] - Phi0[1])*(Phi0[0] - Phi0[1])*(Phi0[0] - Phi0[N/2-2])*(Phi0[0] - Phi0[N/2-2])/(16.0*step_x);
        }
        else{
            Pi1[i] = k*(Phi0[i+1] - 2*Phi0[i] + Phi0[i-1]) + pow(Phi0[i], params[1]);
            Pi1[i] -= (Phi0[i] - Phi0[i+1])*(Phi0[i] - Phi0[i+1])*(Phi0[i] - Phi0[i-1])*(Phi0[i] - Phi0[i-1])/(16.0*step_x);
        }
    }

    //Copies the transformed array to the original one
    for(int i = 0; i < N; ++i)
        u[i] = Du[i];

    //Frees the memory allocated for the transformed array
    delete[] Du;
}