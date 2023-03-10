/**
 * @file main.cpp
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief Reads the arguments from the command line and solves the specified equations with given parameters and initial conditions
 * @version 2.1
 * @date 2023-01-09
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "Core.h"
#include "Equations.h"

#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>
#include <string.h>

void read_par_file(std::string parfilename, void (**Eq)(double* u, int N, double step_x, double* params), int* N_Eqs, double** params, 
double** IC, double* step_x, int* N, std::string* outfilename, double* cfl, double* tmax, double* write){
    //Reads the parameters file
    //Opens the file
    std::fstream FILE;
    FILE.open(parfilename, std::fstream::in);

    //Checks if the parameter file is open
    if(!FILE.is_open()){
        std::cout << "ERROR: Parameters file could not be opened" << std::endl;
        exit(-1);
    }

    //Declares 2 buffers to read the file
    std::string buffer1, buffer2;

    //Declares a string to save the name of the IC file
    std::string IC_filename;

    //Reads the first 3 lines to get the number of equations of the ODE system, the space step, and the number of points
    while(getline(FILE, buffer1)){
        //Declares a stringstream to read and split the buffer
        std::stringstream stream(buffer1);

        getline(stream, buffer2, ' ');

        //Sets up the equation to be solved
        if(!buffer2.compare("#Eq:")){
            getline(stream, buffer2, ' ');

            if(!buffer2.compare("simple_wave")){
                *Eq = &Wave_Equation;
                *N_Eqs = 2;
            }

            if(!buffer2.compare("simple_wave_4th_order")){
                *Eq = &Wave_Equation_4th_Order;
                *N_Eqs = 2;
            }

            if(!buffer2.compare("non_linear_wave")){
                *Eq = &Non_Linear_Wave_Equation;
                *N_Eqs = 2;
            }

            if(!buffer2.compare("simple_wave_dissipation")){
                *Eq = &Wave_Equation_Dissipation;
                *N_Eqs = 2;
            }

            if(!buffer2.compare("non_linear_wave_dissipation")){
                *Eq = &Non_Linear_Wave_Equation_Dissipation;
                *N_Eqs = 2;
            }

            if(!buffer2.compare("spherical_wave")){
                *Eq = &Spherical_Wave_Equation;
                *N_Eqs = 2;
            }

            if(!buffer2.compare("spherical_wave_dissipation")){
                *Eq = &Spherical_Wave_Equation_Dissipation;
                *N_Eqs = 2;
            }
        }

        //Sets up the parameters needed for the equation from the file
        if(!buffer2.compare("#params:")){
            getline(stream, buffer2, ' ');

            (*params) = new double[stoi(buffer2)];

            for(int i = 0; getline(stream, buffer2, ' '); ++i)
                (*params)[i] = atof(buffer2.c_str());
        }
            
        //Sets up the name of the initial conditions file
        if(!buffer2.compare("#IC:")){
            getline(stream, buffer2, ' ');

            IC_filename = buffer2;
        }
        
        //Sets up the space step from the file
        if(!buffer2.compare("#step_x:")){
            getline(stream, buffer2, ' ');

            *step_x = atof(buffer2.c_str());
        }

        //Sets up the number of points from the file
        if(!buffer2.compare("#NPoints:")){
            getline(stream, buffer2, ' ');

            *N = stoi(buffer2);
        }

        //Sets up the name of the output file
        if(!buffer2.compare("#FN:")){
            getline(stream, buffer2, ' ');

            *outfilename = buffer2;
        }

        //Sets up the cfl constant
        if(!buffer2.compare("#CFL:")){
            getline(stream, buffer2, ' ');

            *cfl = atof(buffer2.c_str());
        }

        //Sets the time until the system is to be solved
        if(!buffer2.compare("#T:")){
            getline(stream, buffer2, ' ');

            *tmax = atof(buffer2.c_str());
        }

        //Sets the iterations in which the data will be written to disk
        if(!buffer2.compare("#W:")){
            getline(stream, buffer2, ' ');

            *write = stoi(buffer2);
        }
    }

    //Allocates memory for the initial conditions
    (*IC) = new double[(*N)*(*N_Eqs)];

    //Opens the IC file
    std::fstream IC_FILE;
    IC_FILE.open(IC_filename, std::fstream::in);

    //Checks if the file is open
    if(!IC_FILE.is_open()){
        std::cout << "ERROR: IC file could not be opened" << std::endl;
        exit(-1);
    }

    //Reads the initial conditions and saves the to the array
    for(int i = 0; getline(IC_FILE, buffer1); ++i){
        std::stringstream stream(buffer1);

        for(int j = 0; j < *N_Eqs; ++j){
            getline(stream, buffer2, ' ');
            (*IC)[i + j*(*N)] = atof(buffer2.c_str());
        }
    }

    //Closes both files
    FILE.close();
    IC_FILE.close();
}

int main(int argc, char** argv){
    //Declares the string to store the parameter file name
    std::string INPUT;

    //checks if the command line arguments are valid
    if(argc == 1){
        std::cout << "ERROR: No parameter file given" << std::endl;
        exit(-1);
    }
    if(argc > 2){
        std::cout << "ERROR: Too many arguments given" << std::endl;
        exit(-1);
    }

    //Reads command line arguments
    //if the command -h is used, display help
    if(!strcmp(argv[1], "-h")){
        std::cout << "this is helpful" << std::endl;
        exit(0);
    }

    //reads the name of the parameter file
    INPUT = argv[1];

    //Declares the variables that will be read from the parameter file
    void (*equation)(double* u, int N, double step_x, double* params);
    int N_Eqs = 1;
    double* params = NULL;
    double* u0 = NULL;
    double step_x = 0;
    int NPoints = 0;
    std::string outfilename = "Output.dat";
    double cfl = 0.5;
    double tmax = 1;
    double write = 1;

    //Reads the parameter file
    read_par_file(INPUT, &equation, &N_Eqs, &params, &u0, &step_x, &NPoints, &outfilename, &cfl, &tmax, &write);

    //Solves the equation specifies with the given parameters
    Runge_Kutta_4(equation, u0, N_Eqs*NPoints, N_Eqs, step_x, params, tmax, cfl*step_x, outfilename, write);

    //Deletes the memory allocated for the initial conditions
    delete[] u0;
    delete[] params;

    return 0;
}