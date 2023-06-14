/**
 * @file Output.h
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief Defines the output functions for the main code
 * @version 1.0
 * @date 2023-04-19
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#include "Output.h"

void Output_Solution(std::fstream* FILE, double* u, int N, int N_Ghosts, double time, double* params){
    // Writes the time stepthat is being saved
    *FILE << "\"Time = " << time << std::endl;

    // Saves the solution to disk (ignoring ghost points)
    for(int i = 0; i < (N/((int) params[0]) - 2*N_Ghosts); ++i)
        *FILE << i*params[2] << " " << u[i + ((int) params[1])*(N/((int) params[0])) + N_Ghosts] << std::endl;

    *FILE << std::endl;
}

void Hamiltonian_Constraint(std::fstream* FILE, double* u, int N, int N_Ghosts, double time, double* params){
    // Writes the time step that is being saved
    *FILE << "\"Time = " << time << std::endl;

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
    Even_Constant_Boundary(A, N/9, 1, params[1]);
    Odd_Constant_Boundary(DA, N/9, 1, params[1]);
    Even_Constant_Boundary(KA, N/9, 1, params[1]);

    Even_Constant_Boundary(B, N/9, 1, params[1]);
    Odd_Constant_Boundary(DB, N/9, 1, params[1]);
    Even_Constant_Boundary(KB, N/9, 1, params[1]);

    Odd_Constant_Boundary(lambda, N/9, 1, params[1]);

    Even_Constant_Boundary(alpha, N/9, 1, params[1]);
    Odd_Constant_Boundary(Dalpha, N/9, 1, params[1]);

    // Allocates memory for the derivative of DB and for the constraint
    double* dr_DA = new double[N/9];
    double* dr_DB = new double[N/9];
    double* dr_lambda = new double[N/9];
    double* Hamiltonean = new double[N/9];

    // Calculates the derivatives necessary for the calculation
    Second_Derivative_2nd_Order(DA, dr_DB, N/9, params[0], 1);
    Second_Derivative_2nd_Order(DB, dr_DB, N/9, params[0], 1);
    Second_Derivative_2nd_Order(lambda, dr_lambda, N/9, params[0], 1);

    // Calculates the Hamiltonean constraint value
    Hamiltonean[N_Ghosts] = -dr_DB[N_Ghosts] - dr_lambda[N_Ghosts] + A[N_Ghosts]*KB[N_Ghosts]*(2.0*KA[N_Ghosts] + KB[N_Ghosts]) + 
        dr_DA[N_Ghosts] - 3.0*dr_DB[N_Ghosts] + 0.5*DA[N_Ghosts]*DB[N_Ghosts] - 0.75*DB[N_Ghosts]*DB[N_Ghosts];

    for(int i = N_Ghosts+1; i < N/9 - N_Ghosts; ++i){
        double inv_r = (1.0/((i-(params[1]/2+1))*params[0]));
        Hamiltonean[i] = -dr_DB[i] - lambda[i]*inv_r + A[i]*KB[i]*(2.0*KA[i] + KB[i]) + inv_r*(DA[i] - 3.0*DB[i]) + 0.5*DA[i]*DB[i] - 0.75*DB[i]*DB[i];
    }

    // Saves the value of the Hamiltonian constraint to disk (ignoring ghost points)
    for(int i = N_Ghosts; i < N/9 - N_Ghosts; ++i)
        *FILE << (i-N_Ghosts)*params[0] << " " << Hamiltonean[i] << std::endl;

    *FILE << std::endl;

    // Frees the memory allocated
    delete[] dr_DA;
    delete[] dr_DB;
    delete[] dr_lambda;
    delete[] Hamiltonean;
}

void Reduction_Constraints(std::fstream* FILE, double* u, int N, int N_Ghosts, double time, double* params){
    // Writes the time step that is being saved
    *FILE << "\"Time = " << time << std::endl;

    // Defines pointers to make the manipulation of the state vector easier
    double* var = &(u[((int) params[2]*3)*(N/9)]);
    double* Dvar = &(u[(1 + (int) params[2]*3)*(N/9)]);

    // Populates the ghost points
    Even_Constant_Boundary(var, N/9, 1, params[1]);
    Odd_Constant_Boundary(Dvar, N/9, 1, params[1]);

    // Allocates memory for the derivative of DB and for the constraint
    double* dr_var = new double[N/9];

    // Calculates the derivatives necessary for the calculation
    First_Derivative_2nd_Order(var, dr_var, N/9, params[0], 1);

    // Saves the value of the Hamiltonian constraint to disk (ignoring ghost points)
    for(int i = N_Ghosts; i < N/9 - N_Ghosts; ++i)
        *FILE << (i-N_Ghosts)*params[0] << " " << dr_var[i]/var[i] - Dvar[i] << std::endl;

    *FILE << std::endl;

    // Frees the memory allocated
    delete[] dr_var;
}

void Momentum_Constraint(std::fstream* FILE, double* u, int N, int N_Ghosts, double time, double* params){
    // Writes the time step that is being saved
    *FILE << "\"Time = " << time << std::endl;

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
    Even_Constant_Boundary(A, N/9, 1, params[1]);
    Odd_Constant_Boundary(DA, N/9, 1, params[1]);
    Even_Constant_Boundary(KA, N/9, 1, params[1]);

    Even_Constant_Boundary(B, N/9, 1, params[1]);
    Odd_Constant_Boundary(DB, N/9, 1, params[1]);
    Even_Constant_Boundary(KB, N/9, 1, params[1]);

    Odd_Constant_Boundary(lambda, N/9, 1, params[1]);

    Even_Constant_Boundary(alpha, N/9, 1, params[1]);
    Odd_Constant_Boundary(Dalpha, N/9, 1, params[1]);

    // Allocates memory for the derivative of DB and for the constraint
    double* dr_KA = new double[N/9];
    double* dr_KB = new double[N/9];
    double* Momentum = new double[N/9];

    // Calculates the derivative of DB
    Second_Derivative_2nd_Order(KA, dr_KA, N/9, params[0], 1);
    Second_Derivative_2nd_Order(KB, dr_KB, N/9, params[0], 1);

    // Calculates the Momentum constraint value
    Momentum[N_Ghosts] = -dr_KB[N_Ghosts] + 0.5*DB[N_Ghosts]*(KA[N_Ghosts] - KB[N_Ghosts]) + dr_KA[N_Ghosts] - dr_KB[N_Ghosts];

    for(int i = N_Ghosts+1; i < N/9 - N_Ghosts; ++i){
        double inv_r = (1.0/((i-(params[1]/2+1))*params[0]));
        Momentum[i] = -dr_KB[i] + (KA[i] - KB[i])*(inv_r + 0.5*DB[i]);
    }

    // Saves the value of the Hamiltonian constraint to disk (ignoring ghost points)
    for(int i = N_Ghosts; i < N/9 - N_Ghosts; ++i)
        *FILE << (i-N_Ghosts)*params[0] << " " << Momentum[i] << std::endl;

    *FILE << std::endl;

    // Frees the memory allocated
    delete[] dr_KA;
    delete[] dr_KB;
    delete[] Momentum;
}

void Debug(std::fstream* FILE, double* u, int N, int N_Ghosts, double time, double* params){
    // Writes the time step that is being saved
    *FILE << "\"Time = " << time << std::endl;

    // Defines pointers to make the manipulation of the state vector easier
    double* A = u;
    double* KA = &(u[N/4]);
    double* B = &(u[(2*N)/4]);
    double* KB = &(u[(3*N)/4]);

    // Allocates memory for the derivatives needed and the KO dissipation
    double* dr_A = new double[N/4];
    double* dr2_A = new double[N/4];
    double* dr_B = new double[N/4];
    double* dr2_B = new double[N/4];

    // Populates the ghost points
    Even_Constant_Boundary(u, N, 4, 2);

    // Calculates the required derivative of the fields
    First_Derivative_2nd_Order(A, dr_A, N/4, params[2], 1);
    Second_Derivative_2nd_Order(A, dr2_A, N/4, params[2], 1);
    First_Derivative_2nd_Order(B, dr_B, N/4, params[2], 1);
    Second_Derivative_2nd_Order(B, dr2_B, N/4, params[2], 1);

    // Saves the solution to disk (ignoring ghost points)
    for(int i = N_Ghosts; i < (N/params[0]) - N_Ghosts; ++i){
        double r = (i-N_Ghosts)*params[2];

        *FILE << r << " " << dr2_A[i] << std::endl;
    }

    *FILE << std::endl;

    delete[] dr_A;
    delete[] dr2_A;
    delete[] dr_B;
    delete[] dr2_B;
}