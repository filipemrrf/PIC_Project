In order to [[Numerically Solving the Einstein Equations in 3+1 Dimensions with Spherical Symmetry|solve the Einstein equations numerically]], we need to start developing our code so it can solve hyperbolic PDE's. To accomplish that goal, we will start by solving the wave equation in 1+1 dimensions with Cauchy boundary conditions and periodic boundaries. That is, solving
$$\frac{\partial^2 u(t,x)}{\partial t^2} = c^2 \frac{\partial^2 u(t,x)}{\partial x^2}$$

while having the initial conditions
$$u(0,x)=f(x)$$
$$\frac{\partial u(0,x)}{\partial t}=g(x)$$
and imposing the boundary condition
$$u(t,x) = u(t,x+l),$$
where $l$ is the length of our system (which was set to 1).

To do that, we separate the 2nd order PDE into a system of 2 1st order ODE's in order to time and use the method of finite differences to get the [[Numerical Derivatives|numerical derivatives]] of the positions.

Doing the separation of the wave equation into a system of ODE's we get
$$\frac{\partial^2 u(t,x)}{\partial t^2} = c^2 \frac{\partial^2 u(t,x)}{\partial x^2} \Leftrightarrow$$
$$\Leftrightarrow \left\{ \begin{array}{@{}l@{}} \frac{\partial u(t,x)}{\partial t} = \Pi (t,x) \\ \frac{\partial \Pi(t,x)}{\partial t} = \Phi(t,x) = c^2 \frac{\partial^2 u (t,x)}{\partial x^2} \end{array} \right.\, $$


Now, we can use the [[Runge-Kutta 4]] method to solve the ODE system.


We will be solving the equation this way for initial conditions with 50, 100, 200, 400 and 800 points, such as
$$\Phi(0,x) = sin(2\pi x)$$
$$\Pi(0,x) = 2\pi \: cos(2\pi x)$$
We will also use a CFL coefficient of 0.25 as to not break the [[The CFL condition|CFL Condition]].

After solving the equation, we check the [[Convergence of a Numerical Method|convergence]] of the solutions, obtaining the following results for 2nd and 4th order finite differencing respectively:

![[simple_wave-2nd_order-norm_convergence.png]]
Norm convergence of the 2nd order accurate in space wave equation

![[simple_wave-4th_order-norm_convergence.png]]
Norm convergence of the 4th order accurate in space wave equation


As we can observe, both accuracies show a clean convergence, as it is almost a straight line and the norm convergence matches the accuracy order of the finite differencing.