/**
 * @file main.cpp
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief 
 * @version 0.1
 * @date 2022-12-12
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "Wave_Equation/Wave_Equation.h"

#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>
#include <string.h>

int main(int argc, char** argv){
    //Sets up the boolean variables to check if the program has enough information to run
    bool IC = false;
    bool EQ = false;

    //Declares the string to store the input file name
    std::string INPUT;

    //Declares the function pointer that will save the equation system to solve
    void (*equation)(double* u, int N, double step_x, double* params);
    int N_Eqs;

    //Sets up the string to store the name of the output file
    std::string filename = "Output.dat";

    //Reads command line arguments
    for(int i = 0; i < argc; ++i){
        if(argv[i][0] == '-'){
            //Using -IC INPUT_FILE.dat sets the initial conditions
            if((argv[i][1] == 'I') && (argv[i][2] == 'C') && (argv[i][3] == '\0')){
                INPUT = argv[i+1];
                IC = true;
            }

            //Using -EQ equation chooses the equation to be solved
            if((argv[i][1] == 'E') && (argv[i][2] == 'Q') && (argv[i][3] == '\0')){
                if(strcmp(argv[i+1], "simple_wave") == 0){
                    equation = &Wave_Equation;
                    N_Eqs = 2;
                }

                EQ = true;
            }

            //Using -FN filename.dat sets the name of the output file
            if((argv[i][1] == 'F') && (argv[i][2] == 'N') && (argv[i][3] == '\0'))
                filename = argv[i+1];

            //Using -h opens the help menu
            if((argv[i][1] == 'h') && (argv[i][2] == '\0')){
                
            }
        }
    }

    //Checks if the minimum required parameters were set
    if(!IC || !EQ){
        std::cout << "ERROR: The program must be called with the -IC and the -EQ arguments. Use the argument -h for help" << std::endl;
        exit(-1);
    }

    //Declares the variables that will be read from the file
    double step_x = 0;
    int NPoints = 0;
    double* params = NULL;

    //Reads the Initial conditions file
    //Opens the file
    std::fstream FILE;
    FILE.open(INPUT, std::fstream::in);

    //Checks if the file is open
    if(!FILE.is_open()){
        std::cout << "ERROR: Input file could not be opened" << std::endl;
        exit(-1);
    }

    //Declares 2 buffers to read the file
    std::string buffer1, buffer2;

    //Reads the first 3 lines to get the number of equations of the ODE system, the space step, and the number of points
    for(int i = 0; i < 4; ++i){
        getline(FILE, buffer1, '\n');

        //Declares a stringstream to read and split the buffer
        std::stringstream stream(buffer1);

        getline(stream, buffer2, ' ');

        //Checks if the initial conditions provided are compatible with the equation selected
        if(!buffer2.compare("#NEq:")){
            getline(stream, buffer2, ' ');

            if(N_Eqs != stoi(buffer2)){
                std::cout << "ERROR: Initial conditions provided are not compatible with the equation specified" << std::endl;
                exit(-1);
            }
        }
            
        //Sets up the space step from the file
        if(!buffer2.compare("#step_x:")){
            getline(stream, buffer2, ' ');

            step_x = stod(buffer2);
        }

        //Sets up the number of points from the file
        if(!buffer2.compare("#NPoints:")){
            getline(stream, buffer2, ' ');

            NPoints = stoi(buffer2);
        }

        //Sets up the parameters needed for the equation from the file
        if(!buffer2.compare("#pars:")){
            getline(stream, buffer2, ' ');

            params = new double[stoi(buffer2)];

            for(int i = 0; getline(stream, buffer2, ' '); ++i)
                params[i] = stod(buffer2);
        }
    }

    //Allocates memory for the initial conditions
    double* u0 = new double[(NPoints*N_Eqs)];

    //Saves the initial conditions to the array
    for(int i = 0; getline(FILE, buffer1); ++i){
        std::stringstream stream(buffer1);

        for(int j = 0; j < N_Eqs; ++j){
            getline(stream, buffer2, ' ');
            u0[i + j*NPoints] = stod(buffer2);
        }
    }

    //Closes the file
    FILE.close();

    //Solves the equation specifies with the given parameters
    Runge_Kutta_4(equation, u0, N_Eqs*NPoints, N_Eqs, step_x, params, 1, 0.5*step_x, filename);

    //Deletes the memory allocated for the initial conditions
    delete[] u0;
    delete[] params;

    return 0;
}