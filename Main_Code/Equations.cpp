/**
 * @file Equations.h
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief Defines the equations that are to be solved (declared in Equations.h)
 * @version 2.0
 * @date 2023-03-23
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#include "Equations.h"

void Wave_Equation(double* u, int N, int Acc, BoundaryFunc* boundary, double* params){
    // Populates the ghost points
    boundary(u, N, 2, Acc);

    // Allocates memory for the 2nd derivative and the KO dissipation
    double* d2u = new double[N];
    double* dissipation = new double[N];

    // Checks the accuracy to be used
    if(Acc == 2){
        // Calculates the derivative of u
        Second_Derivative_2nd_Order(u, d2u, N, params[2], 2);

        // Calculates the artificial dissipation of u
        KO_Dissipation_4th_Order(u, dissipation, N, params[2], 2, params[1]);
    }
    if(Acc == 4){
        // Calculates the derivative of u
        Second_Derivative_4th_Order(u, d2u, N, params[2], 2);

        // Calculates the artificial dissipation of u
        KO_Dissipation_6th_Order(u, dissipation, N, params[2], 2, params[1]);
    }

    // Allocates memory for the transformed array
    double* udot = new double[N];

    // Declarates auxiliary pointers for easier readability of the function
    double* Phi0 = u;
    double* d2Phi0 = d2u;
    double* dissPhi0 = dissipation;

    double* Pi0 = &(u[N/2]);
    double* d2Pi0 = &(d2u[N/2]);
    double* dissPi0 = &(dissipation[N/2]);

    double* Phi1 = udot;
    double* Pi1 = &(udot[N/2]);

    // Transforms the array
    for(int i = 0; i < N/2; ++i){
        Phi1[i] = Pi0[i] + dissPhi0[i];
        Pi1[i] = params[0]*params[0]*d2Phi0[i] + dissPi0[i];
    }

    // Copies the transformed array to the received one
    for(int i = 0; i < N; ++i)
        u[i] = udot[i];

    // Deletes the allocated memory
    delete[] udot;
    delete[] d2u;
    delete[] dissipation;
}

void Non_Linear_Wave_Equation(double* u, int N, int Acc, BoundaryFunc* boundary, double* params){
    // Populates the ghost points
    boundary(u, N, 2, Acc);

    // Allocates memory for the 2nd derivative and the KO dissipation
    double* d2u = new double[N];
    double* dissipation = new double[N];

    // Checks the accuracy to be used
    if(Acc == 2){
        // Calculates the derivative of u
        Second_Derivative_2nd_Order(u, d2u, N, params[3], 2);

        // Calculates the artificial dissipation of u
        KO_Dissipation_4th_Order(u, dissipation, N, params[3], 2, params[2]);
    }
    if(Acc == 4){
        // Calculates the derivative of u
        Second_Derivative_4th_Order(u, d2u, N, params[3], 2);

        // Calculates the artificial dissipation of u
        KO_Dissipation_6th_Order(u, dissipation, N, params[3], 2, params[2]);
    }

    // Allocates memory for the transformed array
    double* udot = new double[N];

    // Declarates auxiliary pointers for easier readability of the function
    double* Phi0 = u;
    double* d2Phi0 = d2u;
    double* dissPhi0 = dissipation;

    double* Pi0 = &(u[N/2]);
    double* d2Pi0 = &(d2u[N/2]);
    double* dissPi0 = &(dissipation[N/2]);

    double* Phi1 = udot;
    double* Pi1 = &(udot[N/2]);

    // Transforms the array
    for(int i = 0; i < N/2; ++i){
        Phi1[i] = Pi0[i] + dissPhi0[i];
        Pi1[i] = params[0]*params[0]*(d2Phi0[i] + pow(Phi0[i], params[1])) + dissPi0[i];
    }

    // Copies the transformed array to the received one
    for(int i = 0; i < N; ++i)
        u[i] = udot[i];

    // Deletes the allocated memory
    delete[] udot;
    delete[] d2u;
    delete[] dissipation;
}

void Spherical_Wave_Equation(double* u, int N, int Acc, BoundaryFunc* boundary, double* params){
    // Populates the ghost points
    boundary(u, N, 2, Acc);

    // Allocates memory for the 1st and 2nd derivatives and for the KO dissipation
    double* du = new double[N];
    double* d2u = new double[N];
    double* dissipation = new double[N];

    // Checks the accuracy to be used
    if(Acc == 2){
        // Calculates the derivatives of u
        First_Derivative_2nd_Order(u, du, N, params[2], 2);
        Second_Derivative_2nd_Order(u, d2u, N, params[2], 2);

        // Calculates the artificial dissipation of u
        KO_Dissipation_4th_Order(u, dissipation, N, params[2], 2, params[1]);
    }
    if(Acc == 4){
        // Calculates the derivatives of u
        First_Derivative_4th_Order(u, du, N, params[2], 2);
        Second_Derivative_4th_Order(u, d2u, N, params[2], 2);

        // Calculates the artificial dissipation of u
        KO_Dissipation_6th_Order(u, dissipation, N, params[2], 2, params[1]);
    }

    // Allocates memory for the transformed array
    double* udot = new double[N];

    // Declarates auxiliary pointers for easier readability of the function
    double* Phi0 = u;
    double* d2Phi0 = d2u;
    double* dissPhi0 = dissipation;

    double* Pi0 = &(u[N/2]);
    double* d2Pi0 = &(d2u[N/2]);
    double* dissPi0 = &(dissipation[N/2]);

    double* Phi1 = udot;
    double* Pi1 = &(udot[N/2]);

    // Saves the number of ghost points
    int N_Ghosts = Acc/2 + 1;

    // Transforms the origin
    Phi1[N_Ghosts] = Pi0[N_Ghosts] + dissPhi0[N_Ghosts];
    Pi1[N_Ghosts] = 3.0*params[0]*params[0]*d2u[N_Ghosts] + dissPi0[N_Ghosts];

    // Transforms the rest of the array
    for(int i = N_Ghosts + 1 ; i < N/2 - N_Ghosts; ++i){
        Phi1[i] = Pi0[i] + dissPhi0[i];
        Pi1[i] = params[0]*params[0]*((2.0/((i-N_Ghosts)*params[2]))*du[i] + d2u[i]) - dissPi0[i];
    }

    // Copies the transformed array to the received one
    for(int i = 0; i < N; ++i)
        u[i] = udot[i];

    // Deletes the allocated memory
    delete[] udot;
    delete[] d2u;
    delete[] dissipation;
}