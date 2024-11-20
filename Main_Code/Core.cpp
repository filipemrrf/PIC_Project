/**
 * @file Core.h
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief Defines the core algorithms to solve differential equations (declared in Core.h)
 * @version 4.0
 * @date 2024-11-18
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "Core.h"

void First_Derivative_2nd_Order(double* u, double* du, int N, double step_x, int* N_Ghosts, int N_Vars){
    // Calculates 1/(2 step_x) so it doesn't need to calculate it every iteration
    double aux = 1.0/(2.0*step_x);

    // Calculates the number of points each variable has, so it doesn't need to calculate it more than once
    int n_var = N/N_Vars;
    
    // Calculates the derivative for all variables
    for(int i = 0; i < N_Vars; ++i){
        // Calculates the derivative for all the central points
        for(int j = 1; j < n_var-1; ++j)
            du[i*n_var + j] = aux*(u[i*n_var + j+1] - u[i*n_var + j-1]);

        // Chooses how to fill the remaining points
        if(N_Ghosts[0] == 0)
            // Calculates the derivative for the first point
            du[i*n_var] = aux*(u[i*n_var + 3] - 4.0*u[i*n_var + 2] + 7.0*u[i*n_var + 1] - 4.0*u[i*n_var]);
        else
            // Sets the derivative at the ghost points to 0
            for(int j = 0; j < N_Ghosts[0]; ++j)
                du[i*n_var + j] = 0.0;

        if(N_Ghosts[1] == 0)
            // Calculates the derivative for the last point
            du[i*n_var + n_var-1] = aux*(4.0*u[i*n_var + n_var-1] - 7.0*u[i*n_var + n_var-2] + 4.0*u[i*n_var + n_var-3] - u[i*n_var + n_var-4]);
        else
            // Sets the derivative at the ghost points to 0
            for(int j = 0; j < N_Ghosts[1]; ++j)
                du[(i+1)*n_var - 1 - j] = 0.0;
    }
}

void First_Derivative_4th_Order(double* u, double* du, int N, double step_x, int* N_Ghosts, int N_Vars){
    // Calculates 1/(12 step_x) so it doesn't need to calculate it every iteration
    double aux = 1.0/(12.0*step_x);
    
    // Calculates the number of points each variable has, so it doesn't need to calculate it more than once
    int n_var = N/N_Vars;

    // Calculates the derivative for all variables
    for(int i = 0; i < N_Vars; ++i){
        // Calculates the derivative for all non ghost points
        for(int j = 2; j < n_var - 2; ++j)
            du[i*n_var + j] = aux*(-u[i*n_var + j+2] + 8.0*u[i*n_var + j+1] - 8.0*u[i*n_var + j-1] + u[i*n_var + j-2]);

        // Sets the derivative at the ghost points to 0
        for(int j = 0; j < std::max(N_Ghosts[0], 2); ++j)
            du[i*n_var + j] = 0.0;
        
        for(int j = 0; j < std::max(N_Ghosts[1], 2); ++j)
            du[(i+1)*n_var - 1 - j] = 0.0;
    }
}

void Second_Derivative_2nd_Order(double* u, double* du, int N, double step_x, int* N_Ghosts, int N_Vars){
    // Calculates 1/step_x^2 so it doesn't need to calculate it every iteration
    double aux = 1.0/(step_x*step_x);

    // Calculates the number of points each variable has, so it doesn't need to calculate it more than once
    int n_var = N/N_Vars;

    // Calculates the derivative for all variables
    for(int i = 0; i < N_Vars; ++i){
        // Calculates the derivative for all non ghost points
        for(int j = 1; j < n_var - 1; ++j)
            du[i*n_var + j] = aux*(u[i*n_var + j+1] - 2.0*u[i*n_var + j] + u[i*n_var + j-1]);

        // Sets the derivative at the ghost points to 0
        for(int j = 0; j < std::max(N_Ghosts[0], 1); ++j)
            du[i*n_var + j] = 0.0;
        
        for(int j = 0; j < std::max(N_Ghosts[1], 1); ++j)
            du[(i+1)*n_var - 1 - j] = 0.0;
    }
}

void Second_Derivative_4th_Order(double* u, double* du, int N, double step_x, int* N_Ghosts, int N_Vars){
    // Calculates 1/(12 step_x^2) so it doesn't need to calculate it every iteration
    double aux = 1.0/(12.0*step_x*step_x);

    // Calculates the number of points each variable has, so it doesn't need to calculate it more than once
    int n_var = N/N_Vars;

    // Calculates the derivative for all variables
    for(int i = 0; i < N_Vars; ++i){
        // Calculates the derivative for all non ghost points
        for(int j = 2; j < n_var - 2; ++j)
            du[i*n_var + j] = aux*(-u[i*n_var + j+2] + 16.0*u[i*n_var + j+1] - 30.0*u[i*n_var + j] + 16.0*u[i*n_var + j-1] - u[i*n_var + j-2]);

        // Sets the derivative at the ghost points to 0
        for(int j = 0; j < std::max(N_Ghosts[0], 2); ++j)
            du[i*n_var + j] = 0.0;

        for(int j = 0; j < std::max(N_Ghosts[1], 2); ++j)
            du[(i+1)*n_var - 1 - j] = 0.0;
    }
}

void Evans_Method(double* x, double* u, double* du, int N, double step_x, int* N_Ghosts, int N_Vars, double p){
    // Calculates the number of points each variable has
    int nvar = N/N_Vars;

    // Calculates the derivative for all variables
    for(int i = 0; i < N_Vars; ++i){
        for(int j = 1; j < nvar-1; ++j)
            du[i*nvar + j] = (p + 1.0)*(pow(x[i*nvar + j+1], p)*u[i*nvar + j+1] - 
            pow(x[i*nvar + j-1], p)*u[i*nvar + j-1])/(pow(x[i*nvar + j+1], p+1) - 
            pow(x[i*nvar + j-1], p+1));
        
        // Sets the derivative at the ghost points to poison
        for(int j = 0; j < std::max(N_Ghosts[0], 1); ++j)
            du[i*nvar + j] = 666.0;
        
        for(int j = 0; j < std::max(N_Ghosts[1], 1); ++j)
            du[(i+1)*nvar - 1 - j] = 666.0;

    }
}

void KO_Dissipation_4th_Order(double* u, double* du, int N, double step_x, int* N_Ghosts, int N_Vars, double dissipation){
    // Calculates 1/(16 step_x) so it doesn't need to calculate it every iteration
    double aux = -dissipation/(16.0*step_x);

    // Calculates the number of points each variable has, so it doesn't need to calculate it more than once
    int n_var = N/N_Vars;

    // Goes through each variable
    for(int i = 0; i < N_Vars; ++i){
        // Calculates the dissipation for all non ghost points
        for(int j = 2; j < n_var-2; ++j){
            int idx = i*n_var + j;
            du[idx] = aux*(u[idx + 2] - 4.0*(u[idx + 1]) + 6.0*(u[idx]) - 4.0*(u[idx - 1]) + u[idx - 2]);
        }
        
        // Sets the dissipation at the ghost points to 0
        for(int j = 0; j < std::max(N_Ghosts[0], 2); ++j)
            du[i*n_var + j] = 0.0;
        
        for(int j = 0; j < std::max(N_Ghosts[1], 2); ++j)
            du[(i+1)*n_var - 1 - j] = 0.0;
    }
}

void KO_Dissipation_6th_Order(double* u, double* du, int N, double step_x, int* N_Ghosts, int N_Vars, double dissipation){
    // Calculates 1/(64 step_x) so it doesn't need to calculate it every iteration
    double aux = dissipation/(64.0*step_x);

    // Calculates the number of points each variable has, so it doesn't need to calculate it more than once
    int n_var = N/N_Vars;

    // Goes through each variable
    for(int i = 0; i < N_Vars; ++i){
        // Calculates the dissipation for all non ghost points
        for(int j = 3; j < n_var-3; ++j)
            du[i*n_var + j] = aux*(u[i*n_var + j+3] - 6.0*u[i*n_var + j+2] + 15.0*u[i*n_var + j+1] - 20.0*u[i*n_var + j] + 15.0*u[i*n_var + j-1] 
                - 6.0*u[i*n_var + j-2] + u[i*n_var + j-3]);

        // Sets the dissipation at the ghost points to 0
        for(int j = 0; j < std::max(N_Ghosts[0], 2); ++j)
            du[i*n_var + j] = 0.0;;
        
        for(int j = 0; j < std::max(N_Ghosts[1], 2); ++j)
            du[(i+1)*n_var - 1 - j] = 0.0;
    }
}

void Periodic_Boundary(double* u, int N, int N_Var, int* N_Ghosts){
    // Calculates the number of points for each variable
    int n_var = N/N_Var;

    // Loops through the multiple variables
    for(int i = 0; i < N_Var; ++i){
        // Loops through the ghost points and populates them
        for(int j = 0; j < N_Ghosts[0]; ++j)
            u[i*n_var + j] = u[(i+1)*n_var - (N_Ghosts[0] + N_Ghosts[1]) + j];
        
        for(int j = 0; j < N_Ghosts[1]; ++j)
            u[(i+1)*n_var - N_Ghosts[1] + j] = u[i*n_var + N_Ghosts[0] + j];
    }
}

void Even_0_Boundary(double* u, int N, int N_Var, int* N_Ghosts){
    // Calculates the number of points for each variable
    int n_var = N/N_Var;

    // Loops through the multiple variables
    for(int i = 0; i < N_Var; ++i){
        // Loops through the ghost points and populates them
        for(int j = 1; j <= N_Ghosts[0]; ++j)
            u[i*n_var + N_Ghosts[0]-j] = u[i*n_var + N_Ghosts[0]+j];

        for(int j = 1; j <= N_Ghosts[1]; ++j)
            u[(i+1)*n_var - N_Ghosts[1] + j-1] = 0.0;
    }
}

void Even_Constant_Boundary(double* u, int N, int N_Var, int* N_Ghosts){
    // Calculates the number of points for each variable
    int n_var = N/N_Var;

    // Loops through the multiple variables
    for(int i = 0; i < N_Var; ++i){
        // Loops through the ghost points and populates them
        for(int j = 1; j <= N_Ghosts[0]; ++j)
            u[i*n_var + N_Ghosts[0]-j] = u[i*n_var + N_Ghosts[0]+j];

        for(int j = 1; j <= N_Ghosts[1]; ++j)
            u[(i+1)*n_var - N_Ghosts[1] + j-1] = u[(i+1)*n_var - N_Ghosts[1] - 2];
    }
}

void Constant_Constant_Boundary(double* u, int N, int N_Var, int* N_Ghosts){
    // Calculates the number of points for each variable
    int n_var = N/N_Var;

    // Loops through the multiple variables
    for(int i = 0; i < N_Var; ++i){
        // Loops through the ghost points and populates them
        for(int j = 1; j <= N_Ghosts[0]; ++j)
            u[i*n_var + N_Ghosts[0]-j] = u[i*n_var + N_Ghosts[0]];

        for(int j = 1; j <= N_Ghosts[1]; ++j)
            u[(i+1)*n_var - N_Ghosts[1] + j-1] = u[(i+1)*n_var - N_Ghosts[1] - 2];
    }
}

void Odd_Constant_Boundary(double* u, int N, int N_Var, int* N_Ghosts){
    // Calculates the number of points for each variable
    int n_var = N/N_Var;

    // Loops through the multiple variables
    for(int i = 0; i < N_Var; ++i){
        // Loops through the ghost points and populates them
        for(int j = 1; j <= N_Ghosts[0]; ++j)
            u[i*n_var + N_Ghosts[0]-j] = -u[i*n_var + N_Ghosts[1]+j];

        for(int j = 1; j <= N_Ghosts[1]; ++j)
            u[(i+1)*n_var - N_Ghosts[1] + j-1] = u[(i+1)*n_var - N_Ghosts[1] - 2];
    }
}

void Even_Extrapolation_Boundary(double* u, int N, int N_Var, int* N_Ghosts){
    // Calculates the number of points for each variable
    int n_var = N/N_Var;

    // Calculates the position of the first ghost point on the right boundary
    int ghost = n_var - N_Ghosts[1];

    // Loops through the multiple variables
    for(int i = 0; i < N_Var; ++i){
        // Loops through the ghost points and populates them
        for(int j = 1; j <= N_Ghosts[0]; ++j)
            u[i*n_var + N_Ghosts[0]-j] = u[i*n_var + N_Ghosts[0]+j];

        for(int j = 0; j < N_Ghosts[1]; ++j)
            u[i*n_var + ghost + j] = 4.0*u[i*n_var + ghost + j-1] - 6.0*u[i*n_var + ghost + j-2] + 4.0*u[i*n_var + ghost + j-3] - u[i*n_var + ghost + j-4];
    }
}

void Odd_Extrapolation_Boundary(double* u, int N, int N_Var, int* N_Ghosts){
    // Calculates the number of points for each variable
    int n_var = N/N_Var;

    // Calculates the position of the first ghost point on the right boundary
    int ghost = n_var - N_Ghosts[1];

    // Loops through the multiple variables
    for(int i = 0; i < N_Var; ++i){
        // Loops through the ghost points and populates them
        for(int j = 1; j <= N_Ghosts[0]; ++j)
            u[i*n_var + N_Ghosts[0]-j] = -u[i*n_var + N_Ghosts[1]+j];

        for(int j = 0; j < N_Ghosts[1]; ++j)
            u[i*n_var + ghost + j] = 4.0*u[i*n_var + ghost + j-1] - 6.0*u[i*n_var + ghost + j-2] + 4.0*u[i*n_var + ghost + j-3] - u[i*n_var + ghost + j-4];
    }
}

void Poison(double* u, int N, int N_Var, int* N_Ghosts){
    // Calculates the number of points for each variable
    int n_var = N/N_Var;

    // Loops through the multiple variables
    for(int i = 0; i < N_Var; ++i){
        // Loops through the ghost points and populates them
        for(int j = 1; j <= N_Ghosts[0]; ++j){
            u[i*n_var + N_Ghosts[0]-j] = 999;
        }
        for(int j = 1; j <= N_Ghosts[1]; ++j){
            u[(i+1)*n_var - N_Ghosts[1] + j-1] = 999;
        }
    }
}

void Runge_Kutta_4(rh_sideFunc* rh_side, double* x, double* IC, int N_IC, int* N_Ghosts, double step_x, int Acc, BoundaryFunc* boundary, double* params_rh_side, 
    OutputFunc** output, int N_output, double** params_output, std::string* out_filenames, int write_con, double tmax, double step_t, double dissipation){

    // Allocates memory for the array of pointers to the output files
    std::fstream* OUT = new std::fstream[N_output];

    // Opens the output files and writes the IC in them
    for(int i = 0; i < N_output; ++i){
        OUT[i].open(out_filenames[i], std::fstream::out);
        OUT[i].precision(std::numeric_limits<double>::max_digits10 + 2);

        output[i](&(OUT[i]), x, IC, N_IC, N_Ghosts, 0.0, params_output[i]);
    }

    // Allocates memory for the auxiliary arrays
    double* prev_state = new double[N_IC];
    double* temp_state = new double[N_IC];

    // Copies the IC to the previous state array
    for(int i = 0; i < N_IC; ++i)
        prev_state[i] = IC[i];

    // Allocates memory for the 4 slopes in every point
    double* K1 = new double[N_IC];
    double* K2 = new double[N_IC];
    double* K3 = new double[N_IC];
    double* K4 = new double[N_IC];

    // Loops through the time
    for(int t = 1; ((double) t)*step_t <= tmax; ++t){
        //Copies the IC array to the temporary state array
        for(int i = 0; i < N_IC; ++i)
            temp_state[i] = prev_state[i];

        // Calculates the slope K1 for every point
        rh_side(x, temp_state, N_IC, N_Ghosts, step_x, Acc, boundary, params_rh_side, dissipation);

        for(int i = 0; i < N_IC; ++i)
            K1[i] = temp_state[i]*step_t;


        // Calculates where the slope K2 is to be calculated at for every point
        for(int i = 0; i < N_IC; ++i)
            temp_state[i] = prev_state[i] + 0.5*K1[i];
        
        // Calculates the slope K2 for every point
        rh_side(x, temp_state, N_IC, N_Ghosts, step_x, Acc, boundary, params_rh_side, dissipation);

        for(int i = 0; i < N_IC; ++i)
            K2[i] = temp_state[i]*step_t;


        // Calculates where the slope K3 is to be calculated at for every point
        for(int i = 0; i < N_IC; ++i)
            temp_state[i] = prev_state[i] + 0.5*K2[i];
        
        // Calculates the slope K3 for every point
        rh_side(x, temp_state, N_IC, N_Ghosts, step_x, Acc, boundary, params_rh_side, dissipation);

        for(int i = 0; i < N_IC; ++i)
            K3[i] = temp_state[i]*step_t;


        // Calculates where the slope K4 is to be calculated at for every point
        for(int i = 0; i < N_IC; ++i)
            temp_state[i] = prev_state[i] + K3[i];
        
        // Calculates the slope K4 for every point
        rh_side(x, temp_state, N_IC, N_Ghosts, step_x, Acc, boundary, params_rh_side, dissipation);

        for(int i = 0; i < N_IC; ++i)
            K4[i] = temp_state[i]*step_t;


        // Calculates the final value of the time evolution
        for(int i = 0; i < N_IC; ++i)
                prev_state[i] = prev_state[i] + (K1[i] + 2.0*K2[i] + 2.0*K3[i] + K4[i])/6.0;

        // Decides if this timestep is to be saved to disk
        for(int i = 0; i < N_output; ++i)
            if((t%write_con) == 0)
                output[i](&(OUT[i]), x, prev_state, N_IC, N_Ghosts, t*step_t, params_output[i]);
    }


    // Closes the output files
    for(int i = 0; i < N_output; ++i)
        OUT[i].close();

    // Frees the memory allocated
    delete[] OUT;
    delete[] K1;
    delete[] K2;
    delete[] K3;
    delete[] K4;
    delete[] prev_state;
    delete[] temp_state;
}