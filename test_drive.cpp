#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <math.h>
#include <sstream>
#include <string>
#include <string.h>

void Test_Equation(double* u, int N, double step_x, double* params){
    //Allocates memory for the transformed array
    double* Du = new double[N];

    //Declarates auxiliary pointers for easier readability of the function
    double* Phi0 = u;
    double* Pi0 = &(u[N/2]);

    double* Phi1 = Du;
    double* Pi1 = &(Du[N/2]);

    //Calculates the auxiliary value k = (c/step_x)^2
    double k = *(params)/step_x;
    k *= k;

    //Transforms the array
    for(int i = 0; i < N/2; ++i){
        Phi1[i] = Pi0[i];

        if(i == 0)
            Pi1[0] = 6.0*k*(Phi0[1]-Phi0[0]);
        else if(i == (N/2-1))
            Pi1[N/2-1] = 0.0;
        else
            Pi1[i] = k*((Phi0[i+1] - 2.0*Phi0[i] + Phi0[i-1]) + (1.0/((double)i))*(Phi0[i+1] - Phi0[i-1]));
    }

    //Copies the transformed array to the original one
    for(int i = 0; i < N; ++i)
        u[i] = Du[i];

    //Frees the memory allocated for the transformed array
    delete[] Du;
}

void read_par_file(std::string parfilename, void (**Eq)(double* u, int N, double step_x, double* params), int* N_Eqs, double** params, double** IC, double* step_x, int* NPoints, std::string* outfilename, double* cfl, double* tmax, double* write){
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

            if(!buffer2.compare("test")){
                *Eq = &Test_Equation;
                *N_Eqs = 2;
            }
            else{
                std::cout << "ERROR: Parameter file is not a test" << std::endl;
                exit(-1);
            }
        }

        //Sets up the parameters needed for the equation from the file
        if(!buffer2.compare("#params:")){
            getline(stream, buffer2, ' ');

            (*params) = new double[stoi(buffer2)];

            for(int i = 0; getline(stream, buffer2, ' '); ++i)
                (*params)[i] = stod(buffer2);
        }
            
        //Sets up the name of the initial conditions file
        if(!buffer2.compare("#IC:")){
            getline(stream, buffer2, ' ');

            IC_filename = buffer2;
        }
        
        //Sets up the space step from the file
        if(!buffer2.compare("#step_x:")){
            getline(stream, buffer2, ' ');

            *step_x = stod(buffer2);
        }

        //Sets up the number of points from the file
        if(!buffer2.compare("#NPoints:")){
            getline(stream, buffer2, ' ');

            *NPoints = stoi(buffer2);
        }

        //Sets up the name of the output file
        if(!buffer2.compare("#FN:")){
            getline(stream, buffer2, ' ');

            *outfilename = buffer2;
        }

        //Sets up the cfl constant
        if(!buffer2.compare("#CFL:")){
            getline(stream, buffer2, ' ');

            *cfl = stod(buffer2);
        }

        //Sets the time until the system is to be solved
        if(!buffer2.compare("#T:")){
            getline(stream, buffer2, ' ');

            *tmax = stod(buffer2);
        }

        //Sets the iterations in which the data will be written to disk
        if(!buffer2.compare("#W:")){
            getline(stream, buffer2, ' ');

            *write = stoi(buffer2);
        }
    }

    //Allocates memory for the initial conditions
    (*IC) = new double[(*NPoints)*(*N_Eqs)];

    //Opens the IC file
    std::fstream IC_FILE;
    IC_FILE.open(IC_filename, std::fstream::in);

    //Checks if the file is open
    if(!IC_FILE.is_open()){
        std::cout << "ERROR: IC file could not be opened" << std::endl;
        exit(-1);
    }

    //Reads the initial conditions and saves the to the array
    for(int i = 0; i < *NPoints; ++i){
        getline(IC_FILE, buffer1);
        std::stringstream stream(buffer1);

        for(int j = 0; j < *N_Eqs; ++j){
            getline(stream, buffer2, ' ');
            (*IC)[i + j*(*NPoints)] = stod(buffer2);
        }
    }

    //Closes both files
    FILE.close();
    IC_FILE.close();
}

void Runge_Kutta_4(void (*f)(double* u, int N, double step_x, double* params), double* IC, int N, int N_Eq, double step_x, double* params, double tmax, double step_t, std::string filename, int w){
    //Opens the file the solution is to be written to
    std::fstream FILE;
    FILE.open(filename, std::fstream::out);

    //Makes sure the values are printed to full precision
    FILE.precision(std::numeric_limits<double>::max_digits10 + 2);

    //Writes the initial conditions to the file
    FILE << "\"Time = 0.0" << std::endl;
    for(int i = 0; i < N/N_Eq; ++i){
        FILE << step_x*i << " ";

        for(int j = 0; j < N_Eq; ++j)
            FILE << IC[j*N/N_Eq + i] << " ";
        
        FILE << std::endl;
    }

    FILE << std::endl;

    //Copies the initial conditions to an auxiliary array
    double* aux = new double[N];

    for(int i = 0; i < N; ++i)
        aux[i] = IC[i];

    //Allocates memory for the 4 slopes in every point
    double* K1 = new double[N];
    double* K2 = new double[N];
    double* K3 = new double[N];
    double* K4 = new double[N];

    //Loops through the time
    for(int t = 1; ((double) t)*step_t <= tmax; ++t){
        //Copies the IC array to K1
        for(int i = 0; i < N; ++i)
            K1[i] = aux[i];

        //Calculates the slope K1 for every point
        f(K1, N, step_x, params);

        for(int i = 0; i < N; ++i)
            K1[i] *= step_t;


        //Calculates where the slope K2 is to be calculated at for every point
        for(int i = 0; i < N; ++i)
            K2[i] = aux[i] + 0.5*K1[i];
        
        //Calculates the slope K2 for every point
        f(K2, N, step_x, params);

        for(int i = 0; i < N; ++i)
            K2[i] *= step_t;

        //Calculates where the slope K3 is to be calculated at for every point
        for(int i = 0; i < N; ++i)
            K3[i] = aux[i] + 0.5*K2[i];
        
        //Calculates the slope K3 for every point
        f(K3, N, step_x, params);

        for(int i = 0; i < N; ++i)
            K3[i] *= step_t;


        //Calculates where the slope K4 is to be calculated at for every point
        for(int i = 0; i < N; ++i)
            K4[i] = aux[i] + K3[i];
        
        //Calculates the slope K4 for every point
        f(K4, N, step_x, params);

        for(int i = 0; i < N; ++i)
            K4[i] *= step_t;

        //Calculates the final value of the time evolution
        for(int i = 0; i < N; ++i)
                aux[i] = aux[i] + (K1[i] + 2.0*K2[i] + 2.0*K3[i] + K4[i])/6.0;

        //Decides if this timestep is to be saved to disk
        if((t%w) == 0){
            //Saves the timestep to disk
            FILE << "\"Time = " << (double)t*step_t << std::endl;

            for(int i = 0; i < N/N_Eq; ++i){
                FILE << step_x*i << " ";

                for(int j = 0; j < N_Eq; ++j)
                    FILE << aux[j*N/N_Eq + i] << " ";
                
                FILE << std::endl;
            }

            FILE << std::endl;
        }
    }

    //Closes the file
    FILE.close();

    //Frees the memory allocated for the 4 slopes in every point and the auxiliary vector
    delete[] K1;
    delete[] K2;
    delete[] K3;
    delete[] K4;
    delete[] aux;
}

void save_test(std::string filename, double* u, int NPoints, double step_x){
    std::fstream OUT;
    OUT.open(filename, std::fstream::out);

    OUT << "\"Time = 0.0\n";

    for(int i = 0; i < NPoints; ++i){
        OUT << step_x*i << " " << u[i] << std::endl;
    }

    OUT.close();
}

int main(){
    //Declares the variables that will be read from the file
    void (*equation)(double* u, int N, double step_x, double* params);
    int N_Eqs = 1;
    double* params = NULL;
    double* u0 = NULL;
    double step_x = 0;
    int NPoints = 0;
    std::string out_filename = "Output.dat";
    double cfl = 0.5;
    double tmax = 1;
    double write = 1;

    //Reads the parameter file
    read_par_file("Parameters.txt", &equation, &N_Eqs, &params, &u0, &step_x, &NPoints, &out_filename, &cfl, &tmax, &write);

    save_test(out_filename, u0, NPoints, step_x);

    Runge_Kutta_4(equation, u0, NPoints*N_Eqs, N_Eqs, step_x, params, tmax, cfl*step_x, out_filename, write);

    //Frees the memory of the allocated array
    delete[] u0;
    delete[] params;

    return 0;
}