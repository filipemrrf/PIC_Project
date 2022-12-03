/**
 * @file Wave_Equation.h
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief Wave_Equation class declaration - This class solves the wave equation with the given parameters
 * @version 0.3
 * @date 2022-11-30
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef WAVE_EQUATION
#define WAVE_EQUATION

#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include <gsl/gsl_const_mksa.h>

#include "ODE_System.h"

class Wave_Equation : private ODE_System{
	public:
		Wave_Equation() = default;
		Wave_Equation(const std::vector<std::vector<double>> &x0, double step_x, double step_t = 1e-3, double c = GSL_CONST_MKSA_SPEED_OF_LIGHT);
		Wave_Equation(const std::function<double(double)> &x0, const std::function<double(double)> &xl0, double xmax, double step_x = 1e-3, double step_t = 1e-3, double c = GSL_CONST_MKSA_SPEED_OF_LIGHT);
		
		void Set_IC(const std::vector<std::vector<double>> &x0, double step_x);
		void Set_IC(std::function<double(double)> x0, std::function<double(double)> xl0, double xmax, double step_x = 1e-3);
		
		void Set_Step_t(double step_t);
		
		void Set_C(double c);
		
		void Solve(double tmax);

		void Save_Sol(std::string filename);

		~Wave_Equation(){};
		
	private:		
		//Wave propagation speed
		double c {GSL_CONST_MKSA_SPEED_OF_LIGHT};
};

#endif