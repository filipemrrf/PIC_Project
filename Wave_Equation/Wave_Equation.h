/**
 * @file Wave_Equation.h
 * @author Filipe Ficalho (filipe.ficalho@tecnico.ulisboa.pt)
 * @brief Wave_Equation class declaration - This class solves the wave equation with the given parameters
 * @version 0.1
 * @date 2022-11-14
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <functional>
#include <vector>

class Wave_Equation {
	public:
		//constructors
		Wave_Equation() = default;
		Wave_Equation(std::vector<double> x0, std::vector<double> xl0, double step_x, double step_t = 1e-3, double c = 299792458);
		Wave_Equation(std::function<double(double)> x0, std::function<double(double)> xl0, double xmax, double step_x = 1e-3, double step_t = 1e-3, double c = 299792458);
		
		//setters
		//initial conditions
		void Set_IC(std::vector<double> x0, std::vector<double> xl0, double step_x);
		void Set_IC(std::function<double(double)> x0, std::function<double(double)> xl0, double xmax, double step_x = 1e-3);
		
		//step
		void Set_Step(double step_t);
		
		//wave propagation speed
		void Set_C(double c);
		
		//solve wave equation
		void Solve(double tmax);

		//destructor
		~Wave_Equation();
		
	private:
		//Initial conditions
		std::vector<double> x0; //x(t=0);
		std::vector<double> xl0; //x'(t=0):
		
		//Step
		double step_t {1e-3};
		double step_x {1e-3};
		
		//Wave propagation speed
		double c {299792458};
		
		//Solution
		std::vector<std::vector<double>> sol;
};