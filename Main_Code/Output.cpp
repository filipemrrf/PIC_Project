/**
 * @file Output.h
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief Defines the output functions for the main code
 * @version 3.0
 * @date 2024-11-21
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "Output.h"
#include "Equations.h"

void Output_Solution(std::fstream* FILE, double* x, double* u, int N, int* N_Ghosts, double t, double* params){
    // Defines the parameters for the output
    int N_Fields = (int) params[0];
    int Field_Select = (int) params[1];

    // Calculates the size of the fields
    int N_Fields_Size = N/N_Fields;

    // Writes the time step that is being saved
    *FILE << "\"Time = " << t << std::endl;

    // Saves the solution to disk (ignoring ghost points)
    for(int i = N_Ghosts[0]; i < (N_Fields_Size - N_Ghosts[1]); ++i)
        *FILE << x[i] << " " << u[i + Field_Select*N_Fields_Size] << std::endl;

    *FILE << std::endl;
}

void Output_RHS(std::fstream* FILE, double* x, double* u, int N, int* N_Ghosts, double t, double* params){
    // Defines the parameters for the output
    int N_Fields = (int) params[0];
    int Field_Select = (int) params[1];
    double step_x = params[2];
    double dissipation = params[3];

    // Calculates the size of the fields
    int N_Fields_Size = N/N_Fields;

    // Writes the time step that is being saved
    *FILE << "\"Time = " << t << std::endl;

    // Allocates memory for the RHS
    double* RHS = new double[N]();

    // Copies the state data to the RHS array
    for(int i = 0; i < N; ++i)
        RHS[i] = u[i];

    // Calculates the RHS
    Compact_Wave_Equation(x, RHS, N, N_Ghosts, step_x, 2, Poison, params, dissipation); /*Needs some fixing*/

    // Saves the solution to disk (ignoring ghost points)
    for(int i = N_Ghosts[0]; i < (N_Fields_Size - N_Ghosts[1]); ++i)
        *FILE << x[i] << " " << RHS[i + Field_Select*N_Fields_Size] << std::endl;

    *FILE << std::endl;

    // Frees the memory allocated
    delete[] RHS;
}

void Output_Dissipation(std::fstream* FILE, double* x, double* u, int N, int* N_Ghosts, double t, double* params){
    // Defines the parameters for the output
    int N_Fields = (int) params[0];
    int Field_Select = (int) params[1];
    double step_x = params[2];
    double dissipation = params[3];

    // Calculates the size of the fields
    int N_Fields_Size = N/N_Fields;

    // Writes the time step that is being saved
    *FILE << "\"Time = " << t << std::endl;

    // Allocates memory for the dissipation
    double* diss = new double[N]();

    // Calculates the dissipation
    KO_Dissipation_4th_Order(u, diss, N, step_x, N_Ghosts, N_Fields, dissipation);

    // Saves the solution to disk (ignoring ghost points)
    for(int i = N_Ghosts[0]; i < (N_Fields_Size - N_Ghosts[1]); ++i)
        *FILE << x[i] << " " << diss[i + Field_Select*N_Fields_Size] << std::endl;

    *FILE << std::endl;

    // Frees the memory allocated
    delete[] diss;
}

void Output_1st_Derivative(std::fstream* FILE, double* x, double* u, int N, int* N_Ghosts, double t, double* params){
    // Defines the parameters for the output
    int N_Fields = (int) params[0];
    int Field_Select = (int) params[1];
    double step_x = params[2];

    // Calculates the size of the fields
    int N_Fields_Size = N/N_Fields;

    // Writes the time step that is being saved
    *FILE << "\"Time = " << t << std::endl;

    // Allocates memory for the 1st derivative
    double* dr_u = new double[N]();

    // Calculates the 1st derivative
    First_Derivative_2nd_Order(u, dr_u, N, step_x, N_Ghosts, N_Fields);

    // Saves the solution to disk (ignoring ghost points)
    for(int i = N_Ghosts[0]; i < (N_Fields_Size - N_Ghosts[1]); ++i)
        *FILE << x[i] << " " << dr_u[i + Field_Select*N_Fields_Size] << std::endl;

    *FILE << std::endl;

    // Frees the memory allocated
    delete[] dr_u;
}

void Output_2nd_Derivative(std::fstream* FILE, double* x, double* u, int N, int* N_Ghosts, double t, double* params){
    // Defines the parameters for the output
    int N_Fields = (int) params[0];
    int Field_Select = (int) params[1];
    double step_x = params[2];

    // Calculates the size of the fields
    int N_Fields_Size = N/N_Fields;

    // Writes the time step that is being saved
    *FILE << "\"Time = " << t << std::endl;

    // Allocates memory for the 1st derivative
    double* dr2_u = new double[N]();

    // Calculates the 1st derivative
    Second_Derivative_2nd_Order(u, dr2_u, N, step_x, N_Ghosts, N_Fields);

    // Saves the solution to disk (ignoring ghost points)
    for(int i = N_Ghosts[0]; i < (N_Fields_Size - N_Ghosts[1]); ++i)
        *FILE << x[i] << " " << dr2_u[i + Field_Select*N_Fields_Size] << std::endl;

    *FILE << std::endl;

    // Frees the memory allocated
    delete[] dr2_u;
}

void Output_Constraint(std::fstream* FILE, double* x, double* u, int N, int* N_Ghosts, double t, double* params){
    // Defines the parameters for the output
    double step_x = params[0];
    int N_Fields = (int) params[1];
    
    // Populates the ghost points
    Even_Extrapolation_Boundary(u, N/N_Fields, 1, N_Ghosts);
    Odd_Extrapolation_Boundary(&(u[N/N_Fields]), N/N_Fields, 1, N_Ghosts);
    Even_Extrapolation_Boundary(&(u[(2*N)/N_Fields]), N/N_Fields, 1, N_Ghosts);   

    // Declarates auxiliary pointers for easier readability of the function
    double* Psi = u;
    double* Phi = &(u[N/N_Fields]);
    double* Pi = &(u[(2*N)/N_Fields]);
    
    // Calculates the 1st derivative of psi
    double* drho_Psi = new double[N]();

    First_Derivative_2nd_Order(Psi, drho_Psi, N/N_Fields, step_x, N_Ghosts, 1);

    // Writes the time step that is being saved
    *FILE << "\"Time = " << t << std::endl;

    if(N_Fields == 5){
        double* H = &(u[(3*N)/5]);
        double* A = &(u[(4*N)/5]);

        // Saves the solution to disk (ignoring ghost points)
        for(int i = N_Ghosts[0]; i < N/5 - N_Ghosts[1]; ++i)
            *FILE << x[i] << " " << H[i]*Pi[i] + A[i]*(H[i]*H[i]-1.0)*drho_Psi[i] - Phi[i] << std::endl;
    }

    if(N_Fields >= 7){
        double* H = &(u[(3*N)/N_Fields]);
        double* Omega = &(u[(4*N)/N_Fields]);
        double* L = &(u[(5*N)/N_Fields]);;

        // Saves the solution to disk (ignoring ghost points)
        for(int i = N_Ghosts[0]; i < N/N_Fields - N_Ghosts[1]; ++i)
            *FILE << x[i] << " " << H[i]*Pi[i] + Omega[i]*Omega[i]*drho_Psi[i]/L[i] - Phi[i] << std::endl;
    }

    *FILE << std::endl;

    // Frees the memory allocated
    delete[] drho_Psi;
}

/*void Hamiltonian_Constraint(std::fstream* FILE, double* u, int N, int N_Ghosts, double time, double* params){
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
}*/

/*void Reduction_Constraints(std::fstream* FILE, double* u, int N, int N_Ghosts, double time, double* params){
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
}*/

/*void Momentum_Constraint(std::fstream* FILE, double* u, int N, int N_Ghosts, double time, double* params){
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
}*/

/*void Debug(std::fstream* FILE, double* u, int N, int N_Ghosts, double time, double* params){
    // Writes the time step that is being saved
    *FILE << "\"Time = " << time << std::endl;

    // Populates the ghost points
    Poison(u, N, 5, 2);

    // Allocates memory for the 2nd derivative and the KO dissipation
    double* drho_u = new double[3*N/5];
    double* dissipation = new double[3*N/5];
    
    // Calculates the derivatives of u
    First_Derivative_2nd_Order(u, drho_u, 3*N/5, params[1], 3);

    // Calculates the artificial dissipation of u
    KO_Dissipation_4th_Order(u, dissipation, 3*N/5, params[1], 3, params[0]);

    // Corrects the dissipation
    for(int i = 0; i < 3; ++i){
        for(int j = 0; j < N_Ghosts+ 2; ++j){
            dissipation[i*(N/5) + j] = 0;
            dissipation[(i+1)*(N/5) - j - 1] = 0;
        }
    }


    // Corrects the derivative at the borders
    for(int i = 0; i < 3; ++i){
       drho_u[i*(N/5) + N_Ghosts] = (u[i*(N/5) + N_Ghosts+3] - 4.0*u[i*(N/5) + N_Ghosts+2] + 7.0*u[i*(N/5) + N_Ghosts+1] - 4.0*u[i*(N/5) + N_Ghosts])/(2.0*params[1]);
       drho_u[(i+1)*(N/5) - N_Ghosts-1] = (4.0*u[(i+1)*(N/5)- N_Ghosts - 1] - 7.0*u[(i+1)*(N/5) - N_Ghosts- 2] + 4.0*u[(i+1)*(N/5)- N_Ghosts - 3] - u[(i+1)*(N/5)- N_Ghosts - 4])/(2.0*params[1]);
    }

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
    for(int i = N_Ghosts; i < N/5 - N_Ghosts; ++i){
        Psi1[i] = -Pi0[i] + dissPsi0[i];
        Phi1[i] = A[i]*(H[i]*drho_Phi0[i] + drho_Pi0[i] - drho_Psi0[i]) + dissPhi0[i];
        Pi1[i] = A[i]*(H[i]*(drho_Pi0[i] - drho_Psi0[i]) + drho_Phi0[i]) + dissPi0[i];
    }

    // Saves the solution to disk (ignoring ghost points)
    for(int i = N_Ghosts; i < (N/5)-N_Ghosts; ++i){
        double r = (i-N_Ghosts)*params[2] -1;

        *FILE << r << " " << Pi1[i] << std::endl;
    }

    *FILE << std::endl;

    // Deletes the allocated memory
    delete[] udot;
    delete[] drho_u;
    delete[] dissipation;
}*/