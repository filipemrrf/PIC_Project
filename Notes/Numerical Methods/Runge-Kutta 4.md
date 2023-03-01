The Runge-Kutta 4 method s used to solve ODE's in the form:
$$\frac{\partial \Phi(t,x)}{\partial t} = \Pi(t,x)$$
However, it may also be used to solve ODE systems in the form:
$$\left\{ \begin{array}{@{}l@{}} \frac{\partial \Phi(t,x)}{\partial t} = \Pi (t,x) \\ \frac{\partial \Pi(t,x)}{\partial t} = f(t,x) \end{array} \right.\, $$

The [[Numerically Solving the Wave Equation in 1+1 Dimensions|Wave Equation]] is an example of an ODE system in that form:
$$\left\{ \begin{array}{@{}l@{}} \frac{\partial \Phi(t,x)}{\partial t} = \Pi (t,x) \\ \frac{\partial \Pi(t,x)}{\partial t} = \frac{\partial^2 \Phi(t,x)}{\partial x^2} \end{array} \right.\, $$


This method uses the slopes at 4 different points to calculate the solution to our ODE.

First we calculate the slope at the point we want to evolve in time according to the equation. This slope is calculated as:
$$K_1(t,x) = \Delta t \cdot \Pi(t,x)$$

Now, using this first slope, we evolve the function in time for half a step and then we calculate the second slope at this point:
$$K_2(t,x) = \Delta t \cdot \Pi(t+\Delta t/2, y_t +K_1(t,x)/2)$$

Coming back to the first point, we evolve the function again but using the slope $K_2$. Here, we calculate the slope $K_3$ as:
$$K_3(t,x) = \Delta t \cdot \Pi(t+\Delta t/2, y_t +K_2(t,x)/2)$$

Then, we come back again to the first point and evolve the function a full time step using the slope $K_4$:
$$K_4(t,x) = \Delta t \cdot \Pi(t+\Delta t, y_t +K_3(t,x)$$

Finally, we use all four slopes ($K_1$, $K_2$, $K_3$ and $K_4$) to compute a more accurate estimate to the function evolution:
$$y_{t+\Delta t}(x) = y_t(x) + \frac{1}{6}\left(K_1(t,x) + 2 \cdot K_2(t,x) + 2 \cdot K_3(t,x) + K_4(t,x)\right)$$

![[Images/Runge_Kutta_Scheme.png]]

