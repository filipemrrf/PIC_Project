In order to [[Numericaly Solving the Wave Equation|solve the wave equation numerically]], we need to define the data-types that will be used to store the required variables. We will opt for:

- Initial conditions -> std::vector\<double\>
- step ($\Delta x$ and $\Delta t$) -> double
- Solution -> std::vector\<std::vector\>

Even though using arrays instead of std::vectors would speed up the code by a bit, they need more careful usage and probably the time saved in the runtime would be spent debugging the code. Because of this, vectors will be used (at least in a first stage).

We also need constructors and destructors for the class, as well as getters and setters for the parameters and results.

This way, the class should have this header:

>class Wave_Equation {
>	public:
>		//constructors
>		Wave_Equation();
>		Wave_Equation(std::vector\<double\> x0, std::vector\<double\> xl0);
>		Wave_Equation(std::function\<double(double)\> x0, std::function\<double(double)\> xl0);
>		Wave_Equation(std::vector\<double\> x0, std::vector\<double\> xl0, double step_t = 1e-3, double step_x = 1e-3, double c = 299 792 458);
>		Wave_Equation(std::function\<double(double)\> x0, std::function\<double(double)\> xl0, double step_t = 1e-3, double step_x = 1e-3, double c = 299 792 458);
>		
>		//setters
>		//initial conditions
>		Set_IC(std::vector\<double\> x0, std::vector\<double\> xl0);
>		Set_IC(std::function\<double(double)\> x0, std::function\<double(double)\> xl0);
>		
>		//step
>		Set_Step(double step_t, double step_x);
>		
>		//wave propagation speed
>		Set_C(double c);
>		
>		//destructor
>		~Wave_Equation();
>		
>	private:
>		//Initial conditions
>		std::vector\<double\> x0; //x(t=0);
>		std::vector\<double\> xl0; //x'(t=0):
>		
>		//Step
>		double step_t {1e-3};
>		double step_x {1e-3};
>		
>		//Wave propagation speed
>		double c {299 792 458};
>		
>		//Solution
>		std::vector\<std::vector\<double\>\> sol;
>}