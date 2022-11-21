Our objective is to numerically solve the wave equation with Cauchy boundary conditions. That is
$$\frac{\partial^2 u(t,x)}{\partial t^2} = c^2 \frac{\partial^2 u(t,x)}{\partial x^2}$$
$$u(0,x)=f(x)$$
$$\frac{\partial u(0,x)}{\partial t}=g(x)$$
To do that, we will separate the 2nd order PDE into a system of 2 1st order ODE's in order to time and use the finite difference technique for the position. We will also be solving it in a circle, which means that in a data vector v of size n, v\[0\] = v\[n-1\].


Doing the separation of the wave equation into a system of ODE's we get
$$\frac{\partial^2 u(t,x)}{\partial t^2} = c^2 \frac{\partial^2 u(t,x)}{\partial x^2} \Leftrightarrow$$
$$\Leftrightarrow \left\{ \begin{array}{@{}l@{}} \frac{\partial u(t,x)}{\partial t} = \Pi (t,x) \\ \frac{\partial \Pi (t,x)}{\partial t} = c^2\frac{\partial^2 u(t,x)}{\partial x^2}\end{array} \right.\, $$


Using the Tayor series to approximate small displacements in $x$ while keeping time contant, we get
$$u(t, x+\Delta x) = u(t,x) + \frac{\partial u(t,x)}{\partial x}\Delta x + \frac{1}{2}\frac{\partial^2u(t,x)}{\partial x^2} \Delta x^2+O(\Delta x^3)$$
$$u(t,x-\Delta x) = u(t,x) - \frac{\partial u(t,x)}{\partial x}\Delta x + \frac{1}{2}\frac{\partial^2u(t,x)}{\partial x^2} \Delta x^2-O(\Delta x^3)$$
By adding both expressions, we get
$$u(t,x+\Delta x)+u(t,x-\Delta x) = 2u(t,x)+\frac{\partial^2u(t,x)}{\partial x^2}\Delta x^2 +O(\Delta x^4)\Leftrightarrow$$
$$\Leftrightarrow \frac{\partial^2u(t,x)}{\partial x^2} \approx \frac{u(t,x+\Delta x)-2u(t,x)+u(t,x-\Delta x)}{\Delta x^2}$$


Substituting this into the ODE system, we get
$$\left\{ \begin{array}{@{}l@{}} \frac{d u(t,x)}{d t} = \Pi (t,x) \\ \frac{d \Pi (t,x)}{d t} = \left( \frac{c}{\Delta x} \right)^2 \left( u(t,x+\Delta x)-2u(t,x)+u(t,x-\Delta x) \right) \end{array} \right.\, $$


Now, we can use the [[Runge-Kutta method]] to solve the ODE system.