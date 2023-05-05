The next step towards our [[Numerically Solving the Einstein Equations in 3+1 Dimensions with Spherical Symmetry||goal]] is to solve a Non Linear Wave Equation, since we need [[Kreiss-Oliger Dissipation||artificial dissipation]] for both this problem and the Einstein Equations, as they are non linear.

First, we will start by solving the equation without dissipation so we can add it later into a fully functional code, as to reduce the number of sources of problems in the code.

Similarly as before, we will be solving the non linear wave equation in 1+1 dimensions with Cauchy boundary conditions and periodic boundaries. This is, solving
$$\frac{\partial^2 u(t,x)}{\partial t^2} = c^2 \left( \frac{\partial^2 u(t,x)}{\partial x^2} + u^m \right)$$
while having the initial conditions
$$u(0,x)=f(x)$$
$$\frac{\partial u(0,x)}{\partial t}=g(x)$$
and imposing the boundary condition
$$u(t,x) = u(t,x+l),$$
where $l$ is the length of our system (which was set to 1).

Using the same method we used for the [[Numerically Solving the Wave Equation in 1+1 Dimensions with Periodic Boundary Conditions||wave equation]], we get the following system of ODE's
$$\left\{ \begin{array}{@{}l@{}} \frac{\partial u(t,x)}{\partial t} = \Pi (t,x) \\ \frac{\partial \Pi(t,x)}{\partial t} = \Phi(t,x) = c^2 \left( \frac{\partial^2 u (t,x)}{\partial x^2} + u^m \right) \end{array} \right.\, $$

[[Runge-Kutta 4||Solving the equation]] with non linear coefficient $m=2$ and initial conditions with 50, 100, 200, 400 and 800 points, such as
$$\Phi(0,x) = sin(2\pi x)$$
$$\Pi(0,x) = 2\pi \: cos(2\pi x)$$
We will also use a [[The CFL condition||CFL coefficient]] of 0.25.

Checking the [[Convergence of a Numerical Method|convergence]] of the solutions, we obtain the following results for 2nd and 4th order finite differencing respectively:

![[non_linear_simple_wave-2nd_order-without_dissipation-norm_convergence.png]]
Norm convergence of the 2nd order accurate in space non linear wave equation with non linear power $m = 2$ without artificial dissipation

![[non_linear_simple_wave-4th_order-without_dissipation-norm_convergence.png]]
Norm convergence of the 4th order accurate in space non linear wave equation with non linear power $m = 2$ without artificial dissipation

As we can see, we get a clean convergence in both plots up to a little after $t = 4$. This happens because the solutions to the non linear wave equation diverge, since the non linear term adds energy to the system faster than it can dissipate it by spreading the wave out. By our convergence plots, we can see that we get a clean convergence almost until it clearly diverges. We can conclude that our numerical solution is trustworthy until very near the analytical blow up of the solution.


Now that we have a working non linear wave equation, we will then add [[Kreiss-Oliger Dissipation||artificial KO dissipation]]. This way, our system of equations becomes
$$\left\{ \begin{array}{@{}l@{}} \frac{\partial u(t,x)}{\partial t} = \Pi (t,x) + Q u(t,x) \\ \frac{\partial \Pi(t,x)}{\partial t} = \Phi(t,x) + Q \Pi(t,x) = c^2 \left( \frac{\partial^2 u (t,x)}{\partial x^2} + u^m \right) + Q \Pi(t,x) \end{array} \right.\, ,$$
where $Q$ is the artificial dissipation operator.

Using the same parameters and the same initial conditions as before, we get the following results for 2nd and 4th order finite differencing respectively:

![[non_linear_simple_wave-2nd_order-with_dissipation-norm_convergence.png]]
Norm convergence of the 2nd order accurate in space non linear wave equation with non linear power $m = 2$ with artificial dissipation

![[non_linear_simple_wave-4th_order-with_dissipation-norm_convergence.png]]
Norm convergence of the 4th order accurate in space non linear wave equation with non linear power $m = 2$ with artificial dissipation

As we can see, the artificial dissipation (even though not strictly necessary for this equation) made the code converge in a slightly cleaner way near the singularity. However, it made it converge slightly worse in the first few iterations. However, that is not very relevant since it is still a good convergence and doesn't deviate too much from the desired norm convergence value.