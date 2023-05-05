/**
 * @file main.cpp
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief Reads the arguments from the command line and solves the specified equations with given parameters and initial conditions
 * @version 3.0
 * @date 2023-04-19
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>
#include <string.h>

#include "Core.h"
#include "Equations.h"
#include "Output.h"

/**
 * @brief Reads the parameter file and sets up the necessary variables to run the program
 * 
 * @param parfilename Name of the parameter file
 * @param rh_side (Pointer to the) Equation to be solved
 * @param Acc (Pointer to the) Accuracy order being run
 * @param boundary (Pointer to the) Boundary condition to be used
 * @param params_rh_side (Pointer to the) Parameters of the equation to be solved
 * @param tmax (Pointer to the) Time until the the system is to be solved
 * @param IC (Pointer to the) Initial conditions array
 * @param N_IC (Pointer to the) Size of the initial conditions array
 * @param step_x (Pointer to the) Space step
 * @param output (Pointer to the) Array of output functions that will be executed
 * @param N_output (Pointer to the) Size of the array of output functions
 * @param params_output (Pointer to the) Array of arrays of parameters to the output functions
 * @param out_filenames (Pointer to the) Array of filenames where the output will be saved to
 * @param write_con (Pointer to the) Array of writting conditions for the output
 * @param cfl (Pointer to the) CFL coefficient
 */
void read_par_file(std::string parfilename, rh_sideFunc** rh_side, int* Acc, BoundaryFunc** boundary, double** params_rh_side, double* tmax, double** IC,
    int* N_IC, double* step_x, OutputFunc*** output, int* N_output, double*** params_output, std::string** out_filenames, int** write_con, double* cfl){
    // Reads the parameters file
    // Opens the file
    std::fstream FILE;
    FILE.open(parfilename, std::fstream::in);

    // Checks if the parameter file is open
    if(!FILE.is_open()){
        std::cout << "ERROR: Parameters file could not be opened" << std::endl;
        exit(-1);
    }

    // Declares 2 buffers to read the file
    std::string buffer1, buffer2;

    // Declares a stringstream to read and split the buffer
    std::stringstream stream;

    // Declares a string to save the name of the IC file
    std::string* IC_filename;

    // Declares auxilary variables
    int N_Vars = 0;
    int N_Eq_pars = 0;
    double dissipation = 0.0;

    // Reads one line of the file (ignoring comments and blank lines)
    do{
        getline(FILE, buffer1);
        stream << buffer1;  
        getline(stream, buffer2, ' ');
    }while((buffer2.empty()) || (buffer2[0] == '#'));

    // Sets up the ODE system to be solved
    if(!buffer2.compare("Eq:")){
        getline(stream, buffer2, ' ');

        if(!buffer2.compare("simple_wave")){
            *rh_side = &Wave_Equation;
            N_Vars = 2;
        }

        if(!buffer2.compare("non_linear_simple_wave")){
            *rh_side = &Non_Linear_Wave_Equation;
            N_Vars = 2;
        }

        if(!buffer2.compare("spherical_wave")){
            *rh_side = &Spherical_Wave_Equation;
            N_Vars = 2;
        }
    }
    else{
        std::cout << "ERROR: Equation information not given" << std::endl;
        exit(-1);
    }

    // Reads one line of the file (ignoring comments and blank lines)
    do{
        // Clears the stringstream
        stream.str("");
        stream.clear(); 

        getline(FILE, buffer1);
        stream << buffer1;  
        getline(stream, buffer2, ' ');
    }while((buffer2.empty()) || (buffer2[0] == '#'));

    // Sets up the accuracy of the ODE to be solved
    if(!buffer2.compare("Acc:")){
        getline(stream, buffer2, ' ');

        if(!buffer2.compare("2")){
            *Acc = 2;
        }
        else if(!buffer2.compare("4")){
            *Acc = 4;
        }
        else{
            std::cout << "ERROR: Accuracy provided not supported" << std::endl;
            exit(-1);
        }
    }
    else{
        std::cout << "ERROR: Accuracy information not given" << std::endl;
        exit(-1);
    }

    // Reads one line of the file (ignoring comments and blank lines)
    do{
        // Clears the stringstream
        stream.str("");
        stream.clear(); 

        getline(FILE, buffer1);
        stream << buffer1;  
        getline(stream, buffer2, ' ');
    }while((buffer2.empty()) || (buffer2[0] == '#'));

    // Sets up the boundary conditions of the ODE to be solved
    if(!buffer2.compare("Bound:")){
        getline(stream, buffer2, ' ');

        if(!buffer2.compare("periodic"))
            *boundary = &Periodic_Boundary;
        
        if(!buffer2.compare("even_0"))
            *boundary = &Even_0_Boundary;
    }
    else{
        std::cout << "ERROR: Accuracy information not given" << std::endl;
        exit(-1);
    }
    
    // Reads one line of the file (ignoring comments and blank lines)
    do{
        // Clears the stringstream
        stream.str("");
        stream.clear(); 

        getline(FILE, buffer1);
        stream << buffer1;  
        getline(stream, buffer2, ' ');
    }while((buffer2.empty()) || (buffer2[0] == '#'));

    // Sets up the parameters needed for the equation from the file
    if(!buffer2.compare("params:")){
        getline(stream, buffer2, ' ');

        N_Eq_pars = stoi(buffer2);

        (*params_rh_side) = new double[N_Eq_pars + 2];

        for(int i = 0; getline(stream, buffer2, ' '); ++i)
            (*params_rh_side)[i] = atof(buffer2.c_str());

        (*params_rh_side)[N_Eq_pars] = dissipation;
    }
    else{
        std::cout << "ERROR: Equation parameters not given" << std::endl;
        exit(-1);
    }

    // Reads one line of the file (ignoring comments and blank lines)
    do{
        // Clears the stringstream
        stream.str("");
        stream.clear(); 

        getline(FILE, buffer1);
        stream << buffer1;  
        getline(stream, buffer2, ' ');
    }while((buffer2.empty()) || (buffer2[0] == '#'));


    // Sets the time until the system is to be solved
    if(!buffer2.compare("tmax:")){
        getline(stream, buffer2, ' ');

        *tmax = atof(buffer2.c_str());
    }
    else{
        std::cout << "ERROR: Time until the system is to be solved not given" << std::endl;
        exit(-1);
    }

    // Reads one line of the file (ignoring comments and blank lines)
    do{
        // Clears the stringstream
        stream.str("");
        stream.clear(); 

        getline(FILE, buffer1);
        stream << buffer1;  
        getline(stream, buffer2, ' ');
    }while((buffer2.empty()) || (buffer2[0] == '#'));

    // Sets the name of the IC file
    if(!buffer2.compare("IC:")){
        // Allocates memory for the filenames arrays
        IC_filename = new std::string[N_Vars];

        for(int i = 0; getline(stream, buffer2, ' '); ++i)
            IC_filename[i] = buffer2;
    }
    else{
        std::cout << "ERROR: Initial conditions not given" << std::endl;
        exit(-1);
    }

    // Reads one line of the file (ignoring comments and blank lines)
    do{
        // Clears the stringstream
        stream.str("");
        stream.clear(); 

        getline(FILE, buffer1);
        stream << buffer1;  
        getline(stream, buffer2, ' ');
    }while((buffer2.empty()) || (buffer2[0] == '#'));

    // Sets the number of points of the IC
    if(!buffer2.compare("NPoints:")){
        getline(stream, buffer2, ' ');

        *N_IC = N_Vars*(stoi(buffer2)+2*(((*Acc)/2)+1));

        // Allocates memory for the initial conditions
        *IC = new double[*N_IC];
    }
    else{
        std::cout << "ERROR: Number of points of the initial conditions not given" << std::endl;
        exit(-1);
    }

    // Reads one line of the file (ignoring comments and blank lines)
    do{
        // Clears the stringstream
        stream.str("");
        stream.clear(); 

        getline(FILE, buffer1);
        stream << buffer1;  
        getline(stream, buffer2, ' ');
    }while((buffer2.empty()) || (buffer2[0] == '#'));

    // Sets the space step of the IC
    if(!buffer2.compare("step_x:")){
        getline(stream, buffer2, ' ');

        *step_x = atof(buffer2.c_str());
        (*params_rh_side)[N_Eq_pars + 1] = *step_x;
    }
    else{
        std::cout << "ERROR: Space step not given" << std::endl;
        exit(-1);
    }

    // Reads one line of the file (ignoring comments and blank lines)
    do{
        // Clears the stringstream
        stream.str("");
        stream.clear(); 

        getline(FILE, buffer1);
        stream << buffer1;  
        getline(stream, buffer2, ' ');
    }while((buffer2.empty()) || (buffer2[0] == '#'));

    // Sets the number of files that will be output
    if(!buffer2.compare("NOut:")){
        getline(stream, buffer2, ' ');

        *N_output = stoi(buffer2);

        *params_output = new double*[*N_output];
        *out_filenames = new std::string[*N_output];
        *write_con = new int[*N_output];
        *output = new OutputFunc*[*N_output];
    }
    else{
        std::cout << "ERROR: Number of output files not given" << std::endl;
        exit(-1);
    }

    // Reads one line of the file (ignoring comments and blank lines)
    do{
        // Clears the stringstream
        stream.str("");
        stream.clear(); 

        getline(FILE, buffer1);
        stream << buffer1;  
        getline(stream, buffer2, ' ');
    }while((buffer2.empty()) || (buffer2[0] == '#'));

    // Sets the ouput functions that will bw called
    if(!buffer2.compare("Out_Type:")){
        for(int i = 0; i < *N_output; ++i){
            getline(stream, buffer2, ' ');

            if(!buffer2.compare("solution")){
                for(int j = 0; j < N_Vars; ++j){            
                    (*params_output)[i+j] = new double[3];
                    (*params_output)[i+j][0] = (double) N_Vars;
                    (*params_output)[i+j][1] = (double) j;
                    (*params_output)[i+j][2] = (double) (*step_x);
                
                    (*output)[i+j] = &Output_Solution;
                }

                i += N_Vars;
            }

            if(!buffer2.compare("debug")){          
                    (*params_output)[i] = new double[3];
                    (*params_output)[i][0] = (double) N_Vars;
                    (*params_output)[i][1] = (double) 0;
                    (*params_output)[i][2] = (double) (*step_x);
                
                    (*output)[i] = &Debug;
            }
        }
    }
    else{
        std::cout << "ERROR: Output type not given" << std::endl;
        exit(-1);
    }

    // Reads one line of the file (ignoring comments and blank lines)
    do{
        // Clears the stringstream
        stream.str("");
        stream.clear(); 

        getline(FILE, buffer1);
        stream << buffer1;  
        getline(stream, buffer2, ' ');
    }while((buffer2.empty()) || (buffer2[0] == '#'));


    // Sets the filenames of the output files
    if(!buffer2.compare("Out_Filename:")){
        for(int i = 0; i < *N_output; ++i){
            getline(stream, buffer2, ' ');

            (*out_filenames)[i] = buffer2;
        }
    }
    else{
        std::cout << "ERROR: Output file names not given" << std::endl;
        exit(-1);
    }

    // Reads one line of the file (ignoring comments and blank lines)
    do{
        // Clears the stringstream
        stream.str("");
        stream.clear(); 

        getline(FILE, buffer1);
        stream << buffer1;  
        getline(stream, buffer2, ' ');
    }while((buffer2.empty()) || (buffer2[0] == '#'));


    // Sets the writing condition of the files
    if(!buffer2.compare("W:")){
        for(int i = 0; i < *N_output; ++i){
            for(int i = 0; getline(stream, buffer2, ' '); ++i)
                (*write_con)[i] = atof(buffer2.c_str());
        }
    }
    else{
        std::cout << "ERROR: Output file names not given" << std::endl;
        exit(-1);
    }


    // Reads the remaining (optional) arguments in any order
    while(getline(FILE, buffer1)){
        // Declares a stringstream to read and split the buffer
        std::stringstream stream(buffer1);

        getline(stream, buffer2, ' ');

        // Sets up the cfl coefficient
        if(!buffer2.compare("CFL:")){
            getline(stream, buffer2, ' ');

            *cfl = atof(buffer2.c_str());
        }

        // Sets up the dissipation
        if(!buffer2.compare("dissipation:")){
            getline(stream, buffer2, ' ');

            dissipation = atof(buffer2.c_str());
            (*params_rh_side)[N_Eq_pars] = dissipation;
        }
    }

    for(int i = 0; i < N_Vars; ++i){
        // Opens the IC file
        std::fstream IC_FILE;
        IC_FILE.open(IC_filename[i], std::fstream::in);

        // Checks if the file is open
        if(!IC_FILE.is_open()){
            std::cout << "ERROR: IC file could not be opened" << std::endl;
            exit(-1);
        }

        // Reads the initial conditions and saves the to the array
        for(int j = 0; getline(IC_FILE, buffer1); ++j){
            std::stringstream stream(buffer1);
            getline(stream, buffer2, ' ');
            getline(stream, buffer2, ' ');
            (*IC)[j + i*((*N_IC)/N_Vars) + (*Acc)/2 + 1] = atof(buffer2.c_str());
        }

        // Closes both files
        FILE.close();
        IC_FILE.close();
    }

    delete[] IC_filename;
}

int main(int argc, char** argv){
    // Declares the string to store the parameter file name
    std::string INPUT;

    // Checks if the command line arguments are valid
    if(argc == 1){
        std::cout << "ERROR: No parameter file given" << std::endl;
        exit(-1);
    }
    if(argc > 2){
        std::cout << "ERROR: Too many arguments given" << std::endl;
        exit(-1);
    }

    // If the command -h is used, display help
    if(!strcmp(argv[1], "-h")){
        std::cout << "this is helpful" << std::endl;
        exit(0);
    }

    // Reads the name of the parameter file
    INPUT = argv[1];

    // Declares the variables that will be read from the parameter file
    rh_sideFunc* rh_side = nullptr;
    int Acc = 0; 
    BoundaryFunc* boundary = nullptr;
    double* params_rh_side = nullptr;
    double tmax = 0.0;
    double* IC = nullptr;
    int N_IC = 0;
    double step_x = 0.0;
    OutputFunc** output = nullptr;
    int N_output = 0;
    double** params_output = nullptr;
    std::string* out_filenames = nullptr;
    int* write_con = nullptr;
    double cfl = 0.25;

    // Reads the parameter file
    read_par_file(INPUT, &rh_side, &Acc, &boundary, &params_rh_side, &tmax, &IC, &N_IC, &step_x, &output, &N_output, &params_output, &out_filenames, &write_con, &cfl);

    // Solves the equation specifies with the given parameters
    Runge_Kutta_4(rh_side, IC, N_IC, Acc, boundary, params_rh_side, output, N_output, params_output, out_filenames, write_con, tmax, step_x*cfl);

    // Deletes the memory allocated for the initial conditions
    delete[] params_rh_side;
    delete[] IC;
    delete[] params_output;
    delete[] out_filenames;
    delete[] write_con;

    return 0;
}