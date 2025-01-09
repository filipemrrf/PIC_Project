/**
 * @file Equations.h
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief Defines the equations that are to be solved (declared in Equations.h)
 * @version 4.0
 * @date 2024-11-18
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "Equations.h"

void Wave_Equation(double* x, double* u, int N, int* N_Ghosts, double step_x, int Acc, BoundaryFunc* boundary, double* params, double diss){
    // Populates the ghost points
    boundary(u, N, 2, N_Ghosts);

    // Allocates memory for the 2nd derivative and the KO dissipation
    double* d2u = new double[N];
    double* dissipation = new double[N];

    // Checks the accuracy to be used
    if(Acc == 2){
        // Calculates the derivative of u
        Second_Derivative_2nd_Order(u, d2u, N, step_x, N_Ghosts, 2);

        // Calculates the artificial dissipation of u
        KO_Dissipation_4th_Order(u, dissipation, N, step_x, N_Ghosts, 2, diss);
    }
    if(Acc == 4){
        // Calculates the derivative of u
        Second_Derivative_4th_Order(u, d2u, N, step_x, N_Ghosts, 2);

        // Calculates the artificial dissipation of u
        KO_Dissipation_6th_Order(u, dissipation, N, step_x, N_Ghosts, 2, diss);
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

void Compact_Wave_Equation(double* x, double* u, int N, int* N_Ghosts, double step_x, int Acc, BoundaryFunc* boundary, double* params, double diss){
    // Populates the ghost points
    boundary(u, N, 5, N_Ghosts);

    // Allocates memory for the 2nd derivative and the KO dissipation
    double* drho_u = new double[3*N/5];
    double* dissipation = new double[3*N/5];

    // Calculates the derivatives of u
    First_Derivative_2nd_Order(u, drho_u, 3*N/5, step_x, N_Ghosts, 3);

    // Calculates the artificial dissipation of u
    KO_Dissipation_4th_Order(u, dissipation, 3*N/5, step_x, N_Ghosts, 3, diss);


    // Allocates memory for the transformed array
    double* udot = new double[3*N/5]();

    // Declarates auxiliary pointers for easier readability of the function
    double* Psi0 = u;
    double* drho_Psi0 = drho_u;
    double* dissPsi0 = dissipation;

    double* Phi0 = &(u[N/5]);
    double* drho_Phi0 = &(drho_u[N/5]);
    double* dissPhi0 = &(dissipation[N/5]);

    double* Pi0 = &(u[(2*N)/5]);
    double* drho_Pi0 = &(drho_u[(2*N)/5]);
    double* dissPi0 = &(dissipation[(2*N)/5]);

    double* H = &(u[(3*N)/5]);
    double* A = &(u[(4*N)/5]);
    
    double* Psi1 = udot;
    double* Phi1 = &(udot[N/5]);
    double* Pi1 = &(udot[(2*N)/5]);


    // Transforms the array
    for(int i = N_Ghosts[0]; i < N/5 - N_Ghosts[1]; ++i){
        Psi1[i] = -Pi0[i] + dissPsi0[i];
        Phi1[i] = A[i]*(H[i]*drho_Phi0[i] + drho_Pi0[i]) + dissPhi0[i];
        Pi1[i] = A[i]*(H[i]*drho_Pi0[i] + drho_Phi0[i]) + dissPi0[i];
    }

    // Copies the transformed array to the received one
    for(int i = 0; i < 3*N/5; ++i)
        u[i] = udot[i];

    for(int i = 3*N/5; i < N; ++i)
        u[i] = 0;

    // Deletes the allocated memory
    delete[] udot;
    delete[] drho_u;
    delete[] dissipation;
}

void Spherical_Wave_Equation(double* x, double* u, int N, int* N_Ghosts, double step_x, int Acc, BoundaryFunc* boundary, double* params, double diss){
    // Populates the ghost points
    boundary(u, N, 2, N_Ghosts);

    // Allocates memory for the 1st and 2nd derivatives and for the KO dissipation
    double* du = new double[N];
    double* d2u = new double[N];
    double* dissipation = new double[N];

    // Checks the accuracy to be used
    if(Acc == 2){
        // Calculates the derivatives of u
        First_Derivative_2nd_Order(u, du, N, step_x, N_Ghosts, 2);
        Second_Derivative_2nd_Order(u, d2u, N, step_x, N_Ghosts, 2);

        // Calculates the artificial dissipation of u
        KO_Dissipation_4th_Order(u, dissipation, N, step_x, N_Ghosts, 2, diss);
    }
    if(Acc == 4){
        // Calculates the derivatives of u
        First_Derivative_4th_Order(u, du, N, step_x, N_Ghosts, 2);
        Second_Derivative_4th_Order(u, d2u, N, step_x, N_Ghosts, 2);

        // Calculates the artificial dissipation of u
        KO_Dissipation_6th_Order(u, dissipation, N, step_x, N_Ghosts, 2, diss);
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

    // Transforms the origin
    Phi1[N_Ghosts[0]] = Pi0[N_Ghosts[0]] + dissPhi0[N_Ghosts[0]];
    Pi1[N_Ghosts[0]] = 3.0*params[0]*params[0]*d2u[N_Ghosts[0]] + dissPi0[N_Ghosts[0]];

    // Transforms the rest of the array
    for(int i = N_Ghosts[0] + 1 ; i < N/2 - N_Ghosts[0]; ++i){
        Phi1[i] = Pi0[i] + dissPhi0[i];
        Pi1[i] = params[0]*params[0]*((2.0/((i-N_Ghosts[0])*params[2]))*du[i] + d2u[i]) + dissPi0[i];
    }

    // Copies the transformed array to the received one
    for(int i = 0; i < N; ++i)
        u[i] = udot[i];

    // Deletes the allocated memory
    delete[] udot;
    delete[] d2u;
    delete[] dissipation;
}

void Spherical_Reduced_Wave_Equation(double* x, double* u, int N, int* N_Ghosts, double step_x, int Acc, BoundaryFunc* boundary, double* params, double diss){
    // Populates the ghost points
    Even_Extrapolation_Boundary(u, N/3, 1, N_Ghosts);
    Odd_Extrapolation_Boundary(&(u[N/3]), N/3, 1, N_Ghosts);
    Even_Extrapolation_Boundary(&(u[(2*N)/3]), N/3, 1, N_Ghosts);

    // Allocates memory for the 2nd derivative and the KO dissipation
    double* drho_u = new double[N];
    double* dissipation = new double[N];
    double* evans = new double[N];

    // Calculates the derivatives of u
    First_Derivative_2nd_Order(u, drho_u, N, step_x, N_Ghosts, 3);

    // Calculates the artificial dissipation of u
    KO_Dissipation_4th_Order(u, dissipation, N, step_x, N_Ghosts, 3, diss);

    // Allocates memory for the transformed array
    double* udot = new double[N];

    // Declarates auxiliary pointers for easier readability of the function
    double* Psi0 = u;
    double* drho_Psi0 = drho_u;
    double* dissPsi0 = dissipation;

    double* Phi0 = &(u[N/3]);
    double* drho_Phi0 = &(drho_u[N/3]);
    double* dissPhi0 = &(dissipation[N/3]);

    double* Pi0 = &(u[(2*N)/3]);
    double* drho_Pi0 = &(drho_u[(2*N)/3]);
    double* dissPi0 = &(dissipation[(2*N)/3]);

    double* Psi1 = udot;
    double* Phi1 = &(udot[N/3]);
    double* Pi1 = &(udot[(2*N)/3]);

    // Calculates the Evans Method
    Evans_Method(x, Phi0, evans, N/3, step_x, N_Ghosts, 1, 2.0);

    // Transforms the array
    for(int i = N_Ghosts[0]; i < N/3 - N_Ghosts[1]; ++i){

        Psi1[i] = -Pi0[i] + dissPsi0[i];
        Phi1[i] = -drho_Pi0[i] + dissPhi0[i];
        Pi1[i] = -evans[i] + (2.0*x[i]*Phi0[i])/(1+x[i]*x[i]) + 3.0*Psi0[i]/pow(x[i]*x[i]+1, 2) + dissPi0[i];

    }

    // Copies the transformed array to the received one
    for(int i = 0; i < N; ++i)
        u[i] = udot[i];

    for(int i = N; i < N; ++i)
        u[i] = 0;

    // Deletes the allocated memory
    delete[] udot;
    delete[] drho_u;
    delete[] dissipation;
    delete[] evans;
}

void Spherical_Compact_Wave_Equation(double* x, double* u, int N, int* N_Ghosts, double step_x, int Acc, BoundaryFunc* boundary, double* params, double diss){
    // Populates the ghost points
    Even_Extrapolation_Boundary(u, N/7, 1, N_Ghosts);
    Odd_Extrapolation_Boundary(&(u[N/7]), N/7, 1, N_Ghosts);
    Even_Extrapolation_Boundary(&(u[(2*N)/7]), N/7, 1, N_Ghosts);

    // Allocates memory for the 2nd derivative and the KO dissipation
    double* drho_u = new double[3*N/7]();
    double* dissipation = new double[3*N/7]();
    double* evans = new double[3*N/7]();

    // Calculates the derivatives of u
    First_Derivative_2nd_Order(u, drho_u, 3*N/7, step_x, N_Ghosts, 3);

    // Calculates the artificial dissipation of u
    KO_Dissipation_4th_Order(u, dissipation, 3*N/7, step_x, N_Ghosts, 3, diss);

    // Allocates memory for the transformed array
    double* udot = new double[3*N/7]();

    // Declarates auxiliary pointers for easier readability of the function
    double* Psi0 = u;
    double* drho_Psi0 = drho_u;
    double* dissPsi0 = dissipation;

    double* Phi0 = &(u[N/7]);
    double* drho_Phi0 = &(drho_u[N/7]);
    double* dissPhi0 = &(dissipation[N/7]);

    double* Pi0 = &(u[(2*N)/7]);
    double* drho_Pi0 = &(drho_u[(2*N)/7]);
    double* dissPi0 = &(dissipation[(2*N)/7]);

    double* H = &(u[(3*N)/7]);
    double* Omega = &(u[(4*N)/7]);
    double* L = &(u[(5*N)/7]);
    double* B = &(u[(6*N)/7]);
    
    double* Psi1 = udot;
    double* Phi1 = &(udot[N/7]);
    double* Pi1 = &(udot[(2*N)/7]);

    // Calculates the Evans Method
    Evans_Method(x, Phi0, evans, N/7, step_x, N_Ghosts, 1, 2.0);

    // Transforms the array
    for(int i = N_Ghosts[0]; i < N/7 - N_Ghosts[1]; ++i){
        Psi1[i] = -Pi0[i] + dissPsi0[i];
        Phi1[i] = B[i]*(pow((pow(x[i], 2) + pow(Omega[i], 2)), 2)*(drho_Phi0[i]*H[i] + drho_Pi0[i]) + H[i]*L[i]*Omega[i]*(2.0*x[i]*Phi0[i] - 3.0*Omega[i]*Psi0[i] -
            pow(Omega[i], 2)*drho_Phi0[i]) + H[i]*L[i]*pow(Omega[i], 3)*evans[i]) + dissPhi0[i];
        Pi1[i] = B[i]*(pow((pow(x[i], 2) + pow(Omega[i], 2)), 2)*(drho_Pi0[i]*H[i] + drho_Phi0[i]) + L[i]*Omega[i]*(2.0*x[i]*Phi0[i] - 3.0*Omega[i]*Psi0[i] - 
            pow(Omega[i], 2)*drho_Phi0[i]) + L[i]*pow(Omega[i], 3)*evans[i]) + dissPi0[i];
    }

    // Copies the transformed array to the received one
    for(int i = 0; i < 3*N/7; ++i)
        u[i] = udot[i];

    for(int i = 3*N/7; i < N; ++i)
        u[i] = 0;

    // Deletes the allocated memory
    delete[] udot;
    delete[] drho_u;
    delete[] dissipation;
    delete[] evans;
}

void Power_Non_Linear_Spherical_Compact_Wave_Equation(double* x, double* u, int N, int* N_Ghosts, double step_x, int Acc, BoundaryFunc* boundary, double* params, double diss){
    // Sets the non-linear parameter
    double power = params[0];

    // Populates the ghost points
    Even_Extrapolation_Boundary(u, N/8, 1, N_Ghosts);
    Odd_Extrapolation_Boundary(&(u[N/8]), N/8, 1, N_Ghosts);
    Even_Extrapolation_Boundary(&(u[(2*N)/8]), N/8, 1, N_Ghosts);

    // Allocates memory for the 2nd derivative and the KO dissipation
    double* drho_u = new double[3*N/8]();
    double* dissipation = new double[3*N/8]();
    double* evans = new double[3*N/8]();

    // Calculates the derivatives of u
    First_Derivative_2nd_Order(u, drho_u, 3*N/8, step_x, N_Ghosts, 3);

    // Calculates the artificial dissipation of u
    KO_Dissipation_4th_Order(u, dissipation, 3*N/8, step_x, N_Ghosts, 3, diss);

    // Allocates memory for the transformed array
    double* udot = new double[3*N/8]();

    // Declarates auxiliary pointers for easier readability of the function
    double* Psi0 = u;
    double* drho_Psi0 = drho_u;
    double* dissPsi0 = dissipation;

    double* Phi0 = &(u[N/8]);
    double* drho_Phi0 = &(drho_u[N/8]);
    double* dissPhi0 = &(dissipation[N/8]);

    double* Pi0 = &(u[(2*N)/8]);
    double* drho_Pi0 = &(drho_u[(2*N)/8]);
    double* dissPi0 = &(dissipation[(2*N)/8]);

    double* H = &(u[(3*N)/8]);
    double* Omega = &(u[(4*N)/8]);
    double* L = &(u[(5*N)/8]);
    double* A = &(u[(6*N)/8]);
    double* B = &(u[(7*N)/8]);
    
    double* Psi1 = udot;
    double* Phi1 = &(udot[N/8]);
    double* Pi1 = &(udot[(2*N)/8]);

    // Calculates the Evans Method
    Evans_Method(x, Phi0, evans, N/8, step_x, N_Ghosts, 1, 2.0);

    // Transforms the array
    for(int i = N_Ghosts[0]; i < N/8 - N_Ghosts[1]; ++i){
        Psi1[i] = -Pi0[i] + dissPsi0[i];
        Phi1[i] = B[i]*(pow((pow(x[i], 2) + pow(Omega[i], 2)), 2)*(drho_Phi0[i]*H[i] + drho_Pi0[i]) + H[i]*L[i]*Omega[i]*(2.0*x[i]*Phi0[i] - 3.0*Omega[i]*Psi0[i] -
            pow(Omega[i], 2)*drho_Phi0[i]) + H[i]*L[i]*pow(Omega[i], 3)*evans[i]) + H[i]*L[i]*A[i]/(x[i]*x[i] + Omega[i]*Omega[i])*pow(Psi0[i], power) + dissPhi0[i];
        Pi1[i] = B[i]*(pow((pow(x[i], 2) + pow(Omega[i], 2)), 2)*(drho_Pi0[i]*H[i] + drho_Phi0[i]) + L[i]*Omega[i]*(2.0*x[i]*Phi0[i] - 3.0*Omega[i]*Psi0[i] - 
            pow(Omega[i], 2)*drho_Phi0[i]) + L[i]*pow(Omega[i], 3)*evans[i]) + L[i]*A[i]/(x[i]*x[i] + Omega[i]*Omega[i])*pow(Psi0[i], power) + dissPi0[i];
    }

    // Copies the transformed array to the received one
    for(int i = 0; i < 3*N/8; ++i)
        u[i] = udot[i];

    for(int i = 3*N/8; i < N; ++i)
        u[i] = 0;

    // Deletes the allocated memory
    delete[] udot;
    delete[] drho_u;
    delete[] dissipation;
    delete[] evans;
}

void Non_Linear_Wave_Equation(double* x, double* u, int N, int* N_Ghosts, double step_x, int Acc, BoundaryFunc* boundary, double* params, double diss){
    // Populates the ghost points
    boundary(u, N, 2, N_Ghosts);

    // Allocates memory for the 2nd derivative and the KO dissipation
    double* d2u = new double[N];
    double* dissipation = new double[N];

    // Checks the accuracy to be used
    if(Acc == 2){
        // Calculates the derivative of u
        Second_Derivative_2nd_Order(u, d2u, N, step_x, N_Ghosts, 2);

        // Calculates the artificial dissipation of u
        KO_Dissipation_4th_Order(u, dissipation, N, step_x, N_Ghosts, 2, diss);
    }
    if(Acc == 4){
        // Calculates the derivative of u
        Second_Derivative_4th_Order(u, d2u, N, step_x, N_Ghosts, 2);

        // Calculates the artificial dissipation of u
        KO_Dissipation_6th_Order(u, dissipation, N, step_x, N_Ghosts, 2, diss);
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

void Non_Linear_Spherical_Wave_Equation(double* x, double* u, int N, int* N_Ghosts, double step_x, int Acc, BoundaryFunc* boundary, double* params, double diss){
    // Populates the ghost points
    boundary(u, N, 2, N_Ghosts);

    // Allocates memory for the 1st and 2nd derivatives and for the KO dissipation
    double* du = new double[N];
    double* d2u = new double[N];
    double* dissipation = new double[N];

    // Checks the accuracy to be used
    if(Acc == 2){
        // Calculates the derivatives of u
        First_Derivative_2nd_Order(u, du, N, step_x, N_Ghosts, 2);
        Second_Derivative_2nd_Order(u, d2u, N, step_x, N_Ghosts, 2);

        // Calculates the artificial dissipation of u
        KO_Dissipation_4th_Order(u, dissipation, N, step_x, N_Ghosts, 2, diss);
    }
    if(Acc == 4){
        // Calculates the derivatives of u
        First_Derivative_4th_Order(u, du, N, step_x, N_Ghosts, 2);
        Second_Derivative_4th_Order(u, d2u, N, step_x, N_Ghosts, 2);

        // Calculates the artificial dissipation of u
        KO_Dissipation_6th_Order(u, dissipation, N, step_x, N_Ghosts, 2, diss);
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


    // Transforms the origin
    Phi1[N_Ghosts[0]] = Pi0[N_Ghosts[0]] + dissPhi0[N_Ghosts[0]];
    Pi1[N_Ghosts[0]] = 3.0*params[0]*params[0]*(d2Phi0[N_Ghosts[0]] + pow(Phi0[N_Ghosts[0]], params[1])) + dissPi0[N_Ghosts[0]];

    // Transforms the rest of the array
    for(int i = N_Ghosts[0] + 1 ; i < N/2 - N_Ghosts[0]; ++i){
        Phi1[i] = Pi0[i] + dissPhi0[i];
        Pi1[i] = params[0]*params[0]*((2.0/((i-N_Ghosts[0])*params[2]))*du[i] + d2Phi0[i] + pow(Phi0[N_Ghosts[0]], params[1])) + dissPi0[i];
    }

    // Copies the transformed array to the received one
    for(int i = 0; i < N; ++i)
        u[i] = udot[i];

    // Deletes the allocated memory
    delete[] udot;
    delete[] d2u;
    delete[] dissipation;
}

void ADM_Evolution(double* x, double* u, int N, int* N_Ghosts, double step_x, int Acc, BoundaryFunc* boundary, double* params, double diss){
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
    Even_Constant_Boundary(A, N/9, 1, N_Ghosts);
    Odd_Constant_Boundary(DA, N/9, 1, N_Ghosts);
    Even_Constant_Boundary(KA, N/9, 1, N_Ghosts);

    Even_Constant_Boundary(B, N/9, 1, N_Ghosts);
    Odd_Constant_Boundary(DB, N/9, 1, N_Ghosts);
    Even_Constant_Boundary(KB, N/9, 1, N_Ghosts);

    Odd_Constant_Boundary(lambda, N/9, 1, N_Ghosts);

    Even_Constant_Boundary(alpha, N/9, 1, N_Ghosts);
    Odd_Constant_Boundary(Dalpha, N/9, 1, N_Ghosts);

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
        First_Derivative_2nd_Order(KA, dr_KA, N/9, step_x, N_Ghosts, 1);
        First_Derivative_2nd_Order(KB, dr_KB, N/9, step_x, N_Ghosts, 1);
        First_Derivative_2nd_Order(DA, dr_DA, N/9, step_x, N_Ghosts, 1);
        First_Derivative_2nd_Order(DB, dr_DB, N/9, step_x, N_Ghosts, 1);
        First_Derivative_2nd_Order(lambda, dr_lambda, N/9, step_x, N_Ghosts, 1);
        First_Derivative_2nd_Order(Dalpha, dr_Dalpha, N/9, step_x, N_Ghosts, 1);

        // Calculates the artificial dissipation of u
        KO_Dissipation_4th_Order(u, dissipation, (7*N)/9, step_x, N_Ghosts, 7, diss);
    }
    if(Acc == 4){
        // Calculates the required derivative of the fields
        First_Derivative_4th_Order(KA, dr_KA, N/9, step_x, N_Ghosts, 1);
        First_Derivative_4th_Order(KB, dr_KB, N/9, step_x, N_Ghosts, 1);
        First_Derivative_4th_Order(DA, dr_DA, N/9, step_x, N_Ghosts, 1);
        First_Derivative_4th_Order(DB, dr_DB, N/9, step_x, N_Ghosts, 1);
        First_Derivative_4th_Order(lambda, dr_lambda, N/9, step_x, N_Ghosts, 1);
        First_Derivative_4th_Order(Dalpha, dr_Dalpha, N/9, step_x, N_Ghosts, 1);

        // Calculates the artificial dissipation of u
        KO_Dissipation_6th_Order(u, dissipation, (7*N)/9, step_x, N_Ghosts, 7, diss);
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
    A_dot[N_Ghosts[0]] = -2.0*alpha[N_Ghosts[0]]*A[N_Ghosts[0]]*KA[N_Ghosts[0]] + diss_A[N_Ghosts[0]];
    B_dot[N_Ghosts[0]] = -2.0*alpha[N_Ghosts[0]]*B[N_Ghosts[0]]*KB[N_Ghosts[0]] + diss_B[N_Ghosts[0]];

    // Calculates the value of the fields with the evolution equations
    DA_dot[N_Ghosts[0]] = -2.0*alpha[N_Ghosts[0]]*(KA[N_Ghosts[0]]*Dalpha[N_Ghosts[0]] + dr_KA[N_Ghosts[0]]) + diss_DA[N_Ghosts[0]];
    DB_dot[N_Ghosts[0]] = -2.0*alpha[N_Ghosts[0]]*(KB[N_Ghosts[0]]*Dalpha[N_Ghosts[0]] + dr_KB[N_Ghosts[0]]) + diss_DB[N_Ghosts[0]];


    // Calculates the value of the fields with the evolution equations
    KA_dot[N_Ghosts[0]] = -(alpha[N_Ghosts[0]]/A[N_Ghosts[0]])*(dr_Dalpha[N_Ghosts[0]] + dr_DB[N_Ghosts[0]] + Dalpha[N_Ghosts[0]]*Dalpha[N_Ghosts[0]] - 0.5*Dalpha[N_Ghosts[0]]*DA[N_Ghosts[0]] 
        + 0.5*DB[N_Ghosts[0]]*DB[N_Ghosts[0]] - 0.5*DA[N_Ghosts[0]]*DB[N_Ghosts[0]] - A[N_Ghosts[0]]*KA[N_Ghosts[0]]*(KA[N_Ghosts[0]] + 2.0*KB[N_Ghosts[0]]) - dr_DA[N_Ghosts[0]] 
        + 2.0*dr_DB[N_Ghosts[0]]) + diss_KA[N_Ghosts[0]];
    KB_dot[N_Ghosts[0]] = -(alpha[N_Ghosts[0]]/(2.0*A[N_Ghosts[0]]))*(dr_DB[N_Ghosts[0]] + Dalpha[N_Ghosts[0]]*DB[N_Ghosts[0]] + DB[N_Ghosts[0]]*DB[N_Ghosts[0]] - 0.5*DA[N_Ghosts[0]]*DB[N_Ghosts[0]] 
        - dr_DA[N_Ghosts[0]] + 2.0*dr_Dalpha[N_Ghosts[0]] + 4.0*dr_DB[N_Ghosts[0]] + 2.0*dr_lambda[N_Ghosts[0]]) + alpha[N_Ghosts[0]]*KB[N_Ghosts[0]]*(KA[N_Ghosts[0]] + 2.0*KB[N_Ghosts[0]]) 
        + diss_KB[N_Ghosts[0]];


    lambda_dot[N_Ghosts[0]] = 2.0*alpha[N_Ghosts[0]]*A[N_Ghosts[0]]/B[N_Ghosts[0]]*(dr_KB[N_Ghosts[0]] - 0.5*DB[N_Ghosts[0]]*(KA[N_Ghosts[0]]-KB[N_Ghosts[0]])) + diss_lambda[N_Ghosts[0]];

    alpha_dot[N_Ghosts[0]] = 0;
    Dalpha_dot[N_Ghosts[0]] = 0;

    // Transforms the array
    for(int i = N_Ghosts[0]+1; i < N/9-N_Ghosts[1]-1; ++i){
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
    for(int i = N/9 - N_Ghosts[1] - 1; i < N_Ghosts[1] + 1; ++i){
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