/**
 * @file Equations.h
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief Defines the equations that are to be solved (declared in Equations.h)
 * @version 3.0
 * @date 2023-06-13
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
        Pi1[i] = params[0]*params[0]*((2.0/((i-N_Ghosts)*params[2]))*du[i] + d2u[i]) + dissPi0[i];
    }

    // Copies the transformed array to the received one
    for(int i = 0; i < N; ++i)
        u[i] = udot[i];

    // Deletes the allocated memory
    delete[] udot;
    delete[] d2u;
    delete[] dissipation;
}

void Non_Linear_Spherical_Wave_Equation(double* u, int N, int Acc, BoundaryFunc* boundary, double* params){
    // Populates the ghost points
    boundary(u, N, 2, Acc);

    // Allocates memory for the 1st and 2nd derivatives and for the KO dissipation
    double* du = new double[N];
    double* d2u = new double[N];
    double* dissipation = new double[N];

    // Checks the accuracy to be used
    if(Acc == 2){
        // Calculates the derivatives of u
        First_Derivative_2nd_Order(u, du, N, params[3], 2);
        Second_Derivative_2nd_Order(u, d2u, N, params[3], 2);

        // Calculates the artificial dissipation of u
        KO_Dissipation_4th_Order(u, dissipation, N, params[3], 2, params[2]);
    }
    if(Acc == 4){
        // Calculates the derivatives of u
        First_Derivative_4th_Order(u, du, N, params[3], 2);
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

    // Saves the number of ghost points
    int N_Ghosts = Acc/2 + 1;

    // Transforms the origin
    Phi1[N_Ghosts] = Pi0[N_Ghosts] + dissPhi0[N_Ghosts];
    Pi1[N_Ghosts] = 3.0*params[0]*params[0]*(d2Phi0[N_Ghosts] + pow(Phi0[N_Ghosts], params[1])) + dissPi0[N_Ghosts];

    // Transforms the rest of the array
    for(int i = N_Ghosts + 1 ; i < N/2 - N_Ghosts; ++i){
        Phi1[i] = Pi0[i] + dissPhi0[i];
        Pi1[i] = params[0]*params[0]*((2.0/((i-N_Ghosts)*params[2]))*du[i] + d2Phi0[i] + pow(Phi0[N_Ghosts], params[1])) + dissPi0[i];
    }

    // Copies the transformed array to the received one
    for(int i = 0; i < N; ++i)
        u[i] = udot[i];

    // Deletes the allocated memory
    delete[] udot;
    delete[] d2u;
    delete[] dissipation;
}

void ADM_Evolution(double* u, int N, int Acc, BoundaryFunc* boundary, double* params){
    // Calculates the number of ghost points
    int N_ghosts = Acc/2 + 1;

    // Defines pointers to make the manipulation of the state vector easier
    double* A = u;
    double* DA = &(u[N/9]);
    double* KA = &(u[(2*N)/9]);

    double* B = &(u[(3*N)/9]);
    double* DB = &(u[(4*N)/9]);
    double* KB = &(u[(5*N)/9]);

    double* lambda =&(u[(6*N)/9]);

    double* alpha = &(u[(7*N)/9]);
    double* Dalpha = &(u[(8*N)/9]);

    // Populates the ghost points
    Even_Constant_Boundary(A, N/9, 1, Acc);
    Odd_Constant_Boundary(DA, N/9, 1, Acc);
    Even_Constant_Boundary(KA, N/9, 1, Acc);

    Even_Constant_Boundary(B, N/9, 1, Acc);
    Odd_Constant_Boundary(DB, N/9, 1, Acc);
    Even_Constant_Boundary(KB, N/9, 1, Acc);

    Odd_Constant_Boundary(lambda, N/9, 1, Acc);

    Even_Constant_Boundary(alpha, N/9, 1, Acc);
    Odd_Constant_Boundary(Dalpha, N/9, 1, Acc);

    // Allocates memory for the derivatives needed and the KO dissipation
    double* dr_KA = new double[N/9];
    double* dr_KB = new double[N/9];
    double* dr_DA = new double[N/9];
    double* dr_DB = new double[N/9];
    double* dr_lambda = new double[N/9];
    double* dr_Dalpha = new double[N/9];
    double* dissipation = new double[N];

    // Defines pointers to make the manipulation of the dissipation array easier
    double* diss_A = dissipation;
    double* diss_DA = &(dissipation[N/9]);
    double* diss_KA = &(dissipation[(2*N)/9]);

    double* diss_B = &(dissipation[(3*N)/9]);
    double* diss_DB = &(dissipation[(4*N)/9]);
    double* diss_KB = &(dissipation[(5*N)/9]);

    double* diss_lambda = &(dissipation[(6*N)/9]);

    double* diss_alpha = &(dissipation[(7*N)/9]);
    double* diss_Dalpha = &(dissipation[(8*N)/9]);

    // Checks the accuracy to be used
    if(Acc == 2){
        // Calculates the required derivative of the fields
        First_Derivative_2nd_Order(KA, dr_KA, N/9, params[1], 1);
        First_Derivative_2nd_Order(KB, dr_KB, N/9, params[1], 1);
        First_Derivative_2nd_Order(DA, dr_DA, N/9, params[1], 1);
        First_Derivative_2nd_Order(DB, dr_DB, N/9, params[1], 1);
        First_Derivative_2nd_Order(lambda, dr_lambda, N/9, params[1], 1);
        First_Derivative_2nd_Order(Dalpha, dr_Dalpha, N/9, params[1], 1);

        // Calculates the artificial dissipation of u
        KO_Dissipation_4th_Order(u, dissipation, (7*N)/9, params[1], 7, params[0]);
    }
    if(Acc == 4){
        // Calculates the required derivative of the fields
        First_Derivative_4th_Order(KA, dr_KA, N/9, params[1], 1);
        First_Derivative_4th_Order(KB, dr_KB, N/9, params[1], 1);
        First_Derivative_4th_Order(DA, dr_DA, N/9, params[1], 1);
        First_Derivative_4th_Order(DB, dr_DB, N/9, params[1], 1);
        First_Derivative_4th_Order(lambda, dr_lambda, N/9, params[1], 1);
        First_Derivative_4th_Order(Dalpha, dr_Dalpha, N/9, params[1], 1);

        // Calculates the artificial dissipation of u
        KO_Dissipation_6th_Order(u, dissipation, (7*N)/9, params[1], 7, params[0]);
    }

    // Allocates memory for the transformed array
    double* udot = new double[N];

    // Defines pointers to make the manipulation of the transformed state vector easier
    double* A_dot = udot;
    double* DA_dot = &(udot[N/9]);
    double* KA_dot = &(udot[(2*N)/9]);

    double* B_dot = &(udot[(3*N)/9]);
    double* DB_dot = &(udot[(4*N)/9]);
    double* KB_dot = &(udot[(5*N)/9]);

    double* lambda_dot = &(udot[(6*N)/9]);

    double* alpha_dot = &(udot[(7*N)/9]);
    double* Dalpha_dot = &(udot[(8*N)/9]);


    // Transforms the origin
    // Calculates the value of the fields with the evolution equations
    A_dot[N_ghosts] = -2.0*alpha[N_ghosts]*A[N_ghosts]*KA[N_ghosts] + diss_A[N_ghosts];
    B_dot[N_ghosts] = -2.0*alpha[N_ghosts]*B[N_ghosts]*KB[N_ghosts] + diss_B[N_ghosts];

    // Calculates the value of the fields with the evolution equations
    DA_dot[N_ghosts] = -2.0*alpha[N_ghosts]*(KA[N_ghosts]*Dalpha[N_ghosts] + dr_KA[N_ghosts]) + diss_DA[N_ghosts];
    DB_dot[N_ghosts] = -2.0*alpha[N_ghosts]*(KB[N_ghosts]*Dalpha[N_ghosts] + dr_KB[N_ghosts]) + diss_DB[N_ghosts];


    // Calculates the value of the fields with the evolution equations
    KA_dot[N_ghosts] = -(alpha[N_ghosts]/A[N_ghosts])*(dr_Dalpha[N_ghosts] + dr_DB[N_ghosts] + Dalpha[N_ghosts]*Dalpha[N_ghosts] - 0.5*Dalpha[N_ghosts]*DA[N_ghosts] 
        + 0.5*DB[N_ghosts]*DB[N_ghosts] - 0.5*DA[N_ghosts]*DB[N_ghosts] - A[N_ghosts]*KA[N_ghosts]*(KA[N_ghosts] + 2.0*KB[N_ghosts]) - dr_DA[N_ghosts] 
        + 2.0*dr_DB[N_ghosts]) + diss_KA[N_ghosts];
    KB_dot[N_ghosts] = -(alpha[N_ghosts]/(2.0*A[N_ghosts]))*(dr_DB[N_ghosts] + Dalpha[N_ghosts]*DB[N_ghosts] + DB[N_ghosts]*DB[N_ghosts] - 0.5*DA[N_ghosts]*DB[N_ghosts] 
        - dr_DA[N_ghosts] + 2.0*dr_Dalpha[N_ghosts] + 4.0*dr_DB[N_ghosts] + 2.0*dr_lambda[N_ghosts]) + alpha[N_ghosts]*KB[N_ghosts]*(KA[N_ghosts] + 2.0*KB[N_ghosts]) 
        + diss_KB[N_ghosts];


    lambda_dot[N_ghosts] = 2.0*alpha[N_ghosts]*A[N_ghosts]/B[N_ghosts]*(dr_KB[N_ghosts] - 0.5*DB[N_ghosts]*(KA[N_ghosts]-KB[N_ghosts])) + diss_lambda[N_ghosts];

    alpha_dot[N_ghosts] = 0;
    Dalpha_dot[N_ghosts] = 0;

    // Transforms the array
    for(int i = N_ghosts+1; i < N/9-N_ghosts-1; ++i){
        // Calculates 1/r
        double inv_r = (1.0/((i-(Acc/2+1))*params[1]));

        A_dot[i] = -2.0*alpha[i]*A[i]*KA[i] + diss_A[i];
        B_dot[i] = -2.0*alpha[i]*B[i]*KB[i] + diss_B[i];

        DA_dot[i] = -2.0*alpha[i]*(KA[i]*Dalpha[i] + dr_KA[i]) + diss_DA[i];
        DB_dot[i] = -2.0*alpha[i]*(KB[i]*Dalpha[i] + dr_KB[i]) + diss_DB[i];

        KA_dot[i] = -(alpha[i]/A[i])*(dr_Dalpha[i] + dr_DB[i] + Dalpha[i]*Dalpha[i] - 0.5*Dalpha[i]*DA[i] + 0.5*DB[i]*DB[i] - 0.5*DA[i]*DB[i] -
            A[i]*KA[i]*(KA[i] + 2.0*KB[i]) - inv_r*(DA[i] - 2.0*DB[i])) + diss_KA[i];
        KB_dot[i] = -(alpha[i]/(2.0*A[i]))*(dr_DB[i] + Dalpha[i]*DB[i] + DB[i]*DB[i] - 0.5*DA[i]*DB[i] - inv_r*(DA[i] - 2.0*Dalpha[i] - 
            4.0*DB[i]) + 2.0*inv_r*lambda[i]) + alpha[i]*KB[i]*(KA[i] + 2.0*KB[i]) + diss_KB[i];
    
        lambda_dot[i] = 2.0*alpha[i]*A[i]/B[i]*(dr_KB[i] - 0.5*DB[i]*(KA[i]-KB[i])) + diss_lambda[i];

        alpha_dot[i] = 0;//-2.0*alpha[i]*(KA[i]+2.0*KB[i]);
        Dalpha_dot[i] = 0; //-2.0*(dr_KA[i] + 2.0*dr_KB[i]);
    }

    // Makes the equation not evolve "at infinity"
    for(int i = N/9 - N_ghosts - 1; i < N_ghosts + 1; ++i){
        A_dot[i] = 0.0;
        DA_dot[i] = 0.0;
        KA_dot[i] = 0.0;
        B_dot[i] = 0.0;
        DB_dot[i] = 0.0;
        KB_dot[i] = 0.0;
        lambda_dot[i] = 0.0;
        alpha_dot[i] = 0.0;
        Dalpha_dot[i] = 0.0;
    }
    // Copies the transformed array to the received one
    for(int i = 0; i < N; ++i)
        u[i] = udot[i];

    // Frees the allocated memory
    delete[] udot;
    delete[] dr_KA;
    delete[] dr_KB;
    delete[] dr_DA;
    delete[] dr_DB;
    delete[] dr_lambda;
    delete[] dr_Dalpha;
    delete[] dissipation;
}
