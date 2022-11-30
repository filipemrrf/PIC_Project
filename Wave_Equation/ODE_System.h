/**
 * @file ODE_System.h
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief ODE_System class declaration - This class is the base for the classes solving the Wave Equation and the Advection Equation
 * @version 1.0
 * @date 2022-11-30
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef ODE_SYSTEM
#define ODE_SYSTEM

#include <string>
#include <vector>

class ODE_System{
    public:
		ODE_System() = default;
		ODE_System(const std::vector<std::vector<double>> &x0, double step_x, double step_t = 1e-3);

		void Set_IC(const std::vector<std::vector<double>> &x0, double step_x);

		void Set_Step_T(double step_t);

		std::vector<std::vector<std::vector<double>>> Get_Sol();

		//Solves the system of equations using Runge-Kutta 4 (to be defined by the parent class)
        virtual std::vector<std::vector<std::vector<double>>> rk4(double tmax);

		//Saves the solution to a file (to be defined by the parent class)
		virtual void Save_Sol(std::string filename);

		~ODE_System();

    protected:
		//Vector of initial conditions in the form: (x(t=0), x'(t=0), x''(t=0), ...)
		std::vector<std::vector<double>> x0; 
	
		//Step
		double step_t {1e-3}; //Time step
		double step_x; //Space step

		//Solution
		std::vector<std::vector<std::vector<double>>> Sol;
};

#endif