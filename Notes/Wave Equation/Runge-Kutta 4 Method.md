We are trying to solve an ODE using the Runge-Kutta 4 method. The equation we are trying to solve has the form:
$$\frac{\partial \Pi(t,x)}{\partial t} = f(t,x)$$

In the case of the [[Numericaly Solving the Wave Equation|Wave Equation]], we have
$$\left\{ \begin{array}{@{}l@{}} \frac{\partial u(t,x)}{\partial t} = \Pi (t,x) \\ \frac{\partial \Pi(t,x)}{\partial t} = \frac{c^2}{ \Delta x^2} \left(u(t,x+\Delta x) - 2 u(t,x) + u(t,x-\Delta x)\right) \end{array} \right.\, $$


This method uses the slopes at 4 different points to calculate the solution to our ODE.

First we calculate the slope at the point we want to evolve in time according to the equation. This slope is calculated as:
$$K_1 = \Delta t \cdot f(t,x)$$

In the wave equation, this is translated as
$$\left\{ \begin{array}{@{}l@{}} K_{1_1} = \Delta t \cdot \Pi (t,x) \\ K_{1_2} = \Delta t \frac{c^2}{ \Delta x^2} \left(u(t,x+\Delta x) - 2 u(t,x) + u(t,x-\Delta x)\right) \end{array} \right.\, $$


Now, using this first slope, we evolve the function in time for half a step to find $f(t+\Delta t, x)$