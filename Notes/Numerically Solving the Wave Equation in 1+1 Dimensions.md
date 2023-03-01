Our objective is to numerically solve the wave equation with Cauchy boundary conditions. That is
$$\frac{\partial^2 u(t,x)}{\partial t^2} = c^2 \frac{\partial^2 u(t,x)}{\partial x^2}$$
$$u(0,x)=f(x)$$
$$\frac{\partial u(0,x)}{\partial t}=g(x)$$
To do that, we will separate the 2nd order PDE into a system of 2 1st order ODE's in order to time and use the method of finite differences for the positions. We will also be solving it in a circle, which means that in a data vector v of size n, v\[0\] = v\[n-1\].


Doing the separation of the wave equation into a system of ODE's we get
$$\frac{\partial^2 u(t,x)}{\partial t^2} = c^2 \frac{\partial^2 u(t,x)}{\partial x^2} \Leftrightarrow$$
$$\Leftrightarrow \left\{ \begin{array}{@{}l@{}} \frac{\partial u(t,x)}{\partial t} = \Pi (t,x) \\ \frac{\partial \Pi(t,x)}{\partial t} = c^2 \frac{\partial^2 u (t,x)}{\partial x^2} \end{array} \right.\, $$

Doing this, we get
$$\left\{ \begin{array}{@{}l@{}} \frac{\partial u(t,x)}{\partial t} = \Pi (t,x) \\ \frac{\partial \Pi(t,x)}{\partial t} = \frac{c^2}{ \Delta x^2} \left(u(t,x+\Delta x) - 2 u(t,x) + u(t,x-\Delta x)\right) \end{array} \right.\, $$


Now, we can use the [[Runge-Kutta 4]] to solve the ODE system.

Solving the equation this way for initial conditions with 50, 100, 200, 400 and 800 points, then using the data acquired, we check the [[Norm Convergence|norm convergence]] of the solutions, obtaining the following results:


This system is only second order accurate. However, we can make it 4th order accurate.


$$\left\{ \begin{array}{@{}l@{}} \frac{\partial u(t,x)}{\partial t} = \Pi (t,x) \\ \frac{\partial \Pi(t,x)}{\partial t} = \frac{c^2\left(-u(t,x+2\Delta x) + 16u(t,x+\Delta x) - 30u(t,x) + 16u(t,x-\Delta x) - u(t,x-\Delta x)\right)}{12 \Delta x^2} \end{array} \right.\, $$