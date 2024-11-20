/**
 * @file Read_Par.cpp
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief Defines the functions that read the parameter file and store the values in the given pointers (declared in Read_Par.h)
 * @version 1.0
 * @date 2024-11-19
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "Read_Par.h"

void read_par_file(std::string parfilename, rh_sideFunc** rh_side, int* Acc, BoundaryFunc** boundary, double** params_rh_side, double* tmax, double ** x, double** IC,
    int* N_IC, int** N_Ghosts, double* step_x, OutputFunc*** output, int* N_output, double*** params_output, std::string** out_filenames, 
    int* write_con, double* cfl, double* diss){
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
    int N_IC_Files = 0;
    int N_Eq_pars = 0;
    std::string folder = "";


    // Reads the parameters file and sets up the necessary variables
    while(getline(FILE, buffer1)){
        // Declares a stringstream to read and split the buffer
        std::stringstream stream(buffer1);

        // Reads the first word of the line
        getline(stream, buffer2, ' ');

        // Ignores comments and blank lines
        if ((buffer2.empty()) || (buffer2[0] == '#'))
            continue;

        // Sets up the ODE system to be solved
        if(!buffer2.compare("Eq:")){
            getline(stream, buffer2, ' ');

            *rh_side = read_eq(buffer2, &N_Vars);
        }

        // Sets up the accuracy of the ODE to be solved
        if(!buffer2.compare("Acc:")){
            getline(stream, buffer2, ' ');

            *Acc = read_acc(buffer2);
        }

        if(!buffer2.compare("Bound:")){
            getline(stream, buffer2, ' ');

            *boundary = read_bound(buffer2);
        }

        // Sets up the parameters needed for the equation from the file
        if(!buffer2.compare("Params:")){
            getline(stream, buffer2, ' ');

            N_Eq_pars = stoi(buffer2);
            (*params_rh_side) = new double[N_Eq_pars];

            // Reads the parameters and saves them to memory
            for(int i = 0; getline(stream, buffer2, ' '); ++i)
                (*params_rh_side)[i] = atof(buffer2.c_str());
        }

        // Sets the time until the system is to be solved
        if(!buffer2.compare("TMax:")){
            getline(stream, buffer2, ' ');

            *tmax = atof(buffer2.c_str());
        }

        // Sets the name of the IC files
        if(!buffer2.compare("IC:")){
            //Sets the size of the IC array
            getline(stream, buffer2, ' ');

            N_IC_Files = stoi(buffer2);

            // Allocates memory for the filenames arrays
            IC_filename = new std::string[N_IC_Files];

            for(int i = 0; getline(stream, buffer2, ' '); ++i)
                IC_filename[i] = buffer2;
        }

        // Sets the number of points of the IC
        if(!buffer2.compare("NPoints:")){
            getline(stream, buffer2, ' ');

            *N_IC = stoi(buffer2);
        }

        // Sets the number of ghost points
        if(!buffer2.compare("NGhosts:")){
            getline(stream, buffer2, ' ');

            *N_Ghosts = new int[2];
            
            for(int i = 0; getline(stream, buffer2, ' ') || (i < 2); ++i)
                (*N_Ghosts)[i] = stoi(buffer2);
        }

        // Sets the space step of the IC
        if(!buffer2.compare("Step_x:")){
            getline(stream, buffer2, ' ');

            *step_x = atof(buffer2.c_str());
        }

        // Sets the output functions that will be called
        if(!buffer2.compare("Out:")){
            getline(stream, buffer2, ' ');

            *N_output = stoi(buffer2);

            *params_output = new double*[*N_output];
            *out_filenames = new std::string[*N_output];
            *output = new OutputFunc*[*N_output];

            for(int i = 0; getline(stream, buffer2, ' '); )
                read_output(buffer2, output, params_output, out_filenames, N_Vars, &i, *step_x, *diss);
        }

        // Sets the folder where the output will be saved
        if(!buffer2.compare("Folder:")){
            getline(stream, buffer2, ' ');

            folder = buffer2;
        }

        // Sets the writing condition of the files
        if(!buffer2.compare("W:")){
            getline(stream, buffer2, ' ');

            (*write_con) = atoi(buffer2.c_str());
        }

        // Sets up the cfl coefficient
        if(!buffer2.compare("CFL:")){
            getline(stream, buffer2, ' ');

            *cfl = atof(buffer2.c_str());
        }

        // Sets up the dissipation
        if(!buffer2.compare("Dissipation:")){
            getline(stream, buffer2, ' ');

            *diss = atof(buffer2.c_str());
        }
    };

    // Checks if the necessary inputs were given
    if ((*rh_side == nullptr) || (*Acc == 0) || (*boundary == nullptr) || (*tmax == 0.0) || (IC_filename == nullptr) || (*N_IC == 0) || (N_Ghosts == nullptr) || 
    (*step_x == 0.0) || (*output == nullptr) || (*N_output == 0)){
        std::cout << "ERROR: Necessary inputs not given" << std::endl;
        exit(-1);
    }

    //Checks if the number of IC files is the same as the number of variables
    if(N_IC_Files != N_Vars){
        std::cout << "ERROR: Number of IC files is different from the number of variables" << std::endl;
        exit(-1);
    }

    // Fixes the file names
    for(int i = 0; i < *N_output; ++i)
        (*out_filenames)[i] = folder + "/" + (*out_filenames)[i];

    // Corrects the number of points of the IC
    *N_IC = (*N_IC + (*N_Ghosts)[0] + (*N_Ghosts)[1])*N_Vars;

    // Allocates memory for the initial conditions
    *x = new double[(*N_IC)/N_Vars];
    *IC = new double[*N_IC];

    // Saves the IC to memory
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

            // Reads the positions array and saves it to memory
            if (i == 0)
                (*x)[j+(*N_Ghosts)[0]] = atof(buffer2.c_str());

            getline(stream, buffer2, ' ');
            (*IC)[j + i*((*N_IC)/N_Vars) + (*N_Ghosts)[0]] = atof(buffer2.c_str());
        }

        // Closes the IC file
        IC_FILE.close();
    }

    // Closes the parameter file
    FILE.close();

    // Sets the x for the ghost points
    for(int i = 1; i <= (*N_Ghosts)[0]; ++i)
        (*x)[(*N_Ghosts)[0]-i] = (*x)[(*N_Ghosts)[0]] - ((double) i)*(*step_x);
    for(int i = 1; i <= (*N_Ghosts)[1]; ++i)
        (*x)[(*N_IC)/N_Vars - (*N_Ghosts)[1] - 1 + i] = (*x)[(*N_IC)/N_Vars - (*N_Ghosts)[1] - 1] + ((double) i)*(*step_x);

    // Deletes the memory allocated for the IC filenames
    delete[] IC_filename;
}

rh_sideFunc* read_eq(std::string eq_name, int* N_Vars){
    if(!eq_name.compare("simple_wave")){
        *N_Vars = 2;
        return &Wave_Equation;
    }

    else if(!eq_name.compare("non_linear_simple_wave")){
        *N_Vars = 2;
        return &Non_Linear_Wave_Equation;
    }

    else if(!eq_name.compare("spherical_wave")){
        *N_Vars = 2;
        return &Spherical_Wave_Equation;
    }

    else if(!eq_name.compare("spherical_reduced_wave")){
        *N_Vars = 3;
        return &Spherical_Reduced_Wave_Equation;
    }
    
    else if(!eq_name.compare("non_linear_spherical_wave")){
        *N_Vars = 2;
        return &Non_Linear_Spherical_Wave_Equation;
    }

    else if(!eq_name.compare("adm_evolution")){
        *N_Vars = 9;
        return &ADM_Evolution;
    }

    else if(!eq_name.compare("compact_wave_equation")){
        *N_Vars = 5;
        return &Compact_Wave_Equation;
    }

    else if(!eq_name.compare("spherical_compact_wave_equation")){
        *N_Vars = 7;
        return &Spherical_Compact_Wave_Equation;
    }

    else{
        std::cout << "ERROR: Equation provided not supported" << std::endl;
        exit(-1);
    }
}

int read_acc(std::string acc){
    if(!acc.compare("2")){
        return 2;
    }
    else if(!acc.compare("4")){
        return 4;
    }
    else{
        std::cout << "ERROR: Accuracy provided not supported" << std::endl;
        exit(-1);
    }
}

BoundaryFunc* read_bound(std::string bound){
    if(!bound.compare("even_constant"))
        return &Even_Constant_Boundary;

    else if(!bound.compare("periodic"))
        return &Periodic_Boundary;
    
    else if(!bound.compare("even_0"))
        return &Even_0_Boundary;
    
    else if(!bound.compare("constant_constant"))
        return &Constant_Constant_Boundary;

    else if(!bound.compare("poison"))
        return &Poison;

    else{
        std::cout << "ERROR: Boundary condition provided not supported" << std::endl;
        exit(-1);
    }
}

void read_output(std::string data, OutputFunc*** output, double*** params_output, std::string** out_filenames, double N_Vars, int* index, double step_x, double dissipation){
    // Declares an auxilary variable
    int aux = *index;

    if(!data.compare("solution")){
        for(int j = 0; j < N_Vars; ++j){
            (*params_output)[aux+j] = new double[2];
            (*params_output)[aux+j][0] = (double) N_Vars;
            (*params_output)[aux+j][1] = (double) j;
        
            (*output)[aux+j] = &Output_Solution;

            (*out_filenames)[aux+j] = "Field_" + std::to_string(j+1) + ".dat";
        }

        *index += N_Vars;
    }

    if(!data.compare("rhs")){
        for(int j = 0; j < N_Vars; ++j){
            (*params_output)[aux+j] = new double[4];
            (*params_output)[aux+j][0] = (double) N_Vars;
            (*params_output)[aux+j][1] = (double) j;
            (*params_output)[aux+j][2] = (double) step_x;
            (*params_output)[aux+j][3] = (double) dissipation;
        
            (*output)[aux+j] = &Output_RHS;

            (*out_filenames)[aux+j] = "RHS_Field_" + std::to_string(j+1) + ".dat";
        }

        *index += N_Vars;
    }

    if(!data.compare("dissipation")){
        for(int j = 0; j < N_Vars; ++j){
            (*params_output)[aux+j] = new double[4];
            (*params_output)[aux+j][0] = (double) N_Vars;
            (*params_output)[aux+j][1] = (double) j;
            (*params_output)[aux+j][2] = (double) step_x;
            (*params_output)[aux+j][3] = (double) dissipation;
        
            (*output)[aux+j] = &Output_Dissipation;

            (*out_filenames)[aux+j] = "Dissipation_Field_" + std::to_string(j+1) + ".dat";
        }

        *index += N_Vars;
    }

    if(!data.compare("1st_derivative")){
        for(int j = 0; j < N_Vars; ++j){
            (*params_output)[aux+j] = new double[3];
            (*params_output)[aux+j][0] = (double) N_Vars;
            (*params_output)[aux+j][1] = (double) j;
            (*params_output)[aux+j][2] = (double) step_x;
        
            (*output)[aux+j] = &Output_1st_Derivative;

            (*out_filenames)[aux+j] = "1st_Derivative_Field_" + std::to_string(j+1) + ".dat";
        }

        *index += N_Vars;
    }

    if(!data.compare("2nd_derivative")){
        for(int j = 0; j < N_Vars; ++j){
            (*params_output)[aux+j] = new double[3];
            (*params_output)[aux+j][0] = (double) N_Vars;
            (*params_output)[aux+j][1] = (double) j;
            (*params_output)[aux+j][2] = (double) step_x;
        
            (*output)[aux+j] = &Output_2nd_Derivative;

            (*out_filenames)[aux+j] = "2nd_Derivative_Field_" + std::to_string(j+1) + ".dat";
        }

        *index += N_Vars;
    }    

    if(!data.compare("constraint")){
        (*params_output)[aux] = new double[2];
        (*params_output)[aux][0] = (double) step_x;
        (*params_output)[aux][1] = (double) N_Vars;
    
        (*output)[aux] = &Output_Constraint;

        (*out_filenames)[aux] = "Constraint.dat";

        *index += 1;
    }

    /*if(!buffer2.compare("hamiltonian")){ 
        (*params_output)[i] = new double[2];
        (*params_output)[i][0] = (double) (*step_x);
        (*params_output)[i][1] = (double) (*Acc);
    
        (*output)[i] = &Hamiltonian_Constraint;
        

        i++;
    }

    if(!buffer2.compare("momentum")){          
        (*params_output)[i] = new double[2];
        (*params_output)[i][0] = (double) (*step_x);
        (*params_output)[i][1] = (double) (*Acc);
    
        (*output)[i] = &Momentum_Constraint;

        i++;
    }

    if(!buffer2.compare("reduction")){         
        for(int j = 0; j < 2; ++j){
            (*params_output)[i+j] = new double[3];
            (*params_output)[i+j][0] = (double) (*step_x);
            (*params_output)[i+j][1] = (double) (*Acc);
            (*params_output)[i+j][2] = (double) j;
        
            (*output)[i+j] = &Reduction_Constraints;
        }

        i += 2;
    }

    if(!buffer2.compare("debug")){          
        (*params_output)[i] = new double[3];
        (*params_output)[i][0] = (double) N_Vars;
        (*params_output)[i][1] = (double) 0;
        (*params_output)[i][2] = (double) (*step_x);
    
        (*output)[i] = &Debug;

        i++;
    }*/
}