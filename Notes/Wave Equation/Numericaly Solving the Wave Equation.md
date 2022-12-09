Our objective is to numerically solve the wave equation with Cauchy boundary conditions. That is
$$\frac{\partial^2 u(t,x)}{\partial t^2} = c^2 \frac{\partial^2 u(t,x)}{\partial x^2}$$
$$u(0,x)=f(x)$$
$$\frac{\partial u(0,x)}{\partial t}=g(x)$$
To do that, we will separate the 2nd order PDE into a system of 2 1st order ODE's in order to time and use the method of finite differences for the positions. We will also be solving it in a circle, which means that in a data vector v of size n, v\[0\] = v\[n-1\].


Doing the separation of the wave equation into a system of ODE's we get
$$\frac{\partial^2 u(t,x)}{\partial t^2} = c^2 \frac{\partial^2 u(t,x)}{\partial x^2} \Leftrightarrow$$
$$\Leftrightarrow \left\{ \begin{array}{@{}l@{}} \frac{\partial u(t,x)}{\partial t} = \Pi (t,x) \\ \frac{\partial \Pi(t,x)}{\partial t} = c^2 \frac{\partial^2 u (t,x)}{\partial x^2} \end{array} \right.\, $$


Now, we need to use finite differences in order to rewrite the space derivative.

First, we discretize the spatial part of our plane into a grid with N cells with side $\Delta x$, covering this axis of our plane. This discretization is edge centered.

Expanding $u(t,x)$ around $x$, its immediate neighbors ($x+\Delta x$ and $x - \Delta x$), and the ones next to them ($x + 2\Delta x$ and $x - 2 \Delta x$), we get 
$$ u(t, x+\Delta x) = u(x) + \Delta x\frac{\partial u(t, x)}{\partial x} + \frac{1}{2}\Delta x^2\frac{\partial^2 u(t, x)}{\partial x^2} + \mathcal{O}(\Delta x^3)$$
$$u(t,x) = u(t,x)$$
$$ u(t, x+\Delta x) = u(x) - \Delta x\frac{\partial u(t, x)}{\partial x} + \frac{1}{2}\Delta x^2\frac{\partial^2 u(t, x)}{\partial x^2} + \mathcal{O}(\Delta x^3)$$

We can now write these equations in the matrix equation
$$ \textbf{u} = \textbf{M} \cdot \textbf{d} + \textbf{e} $$

In this equation, $\textbf{u} = (u(t, x+\Delta x), u(t,x), u(t,x-\Delta x))$ is a vector containing the function values, $\textbf{d} = (u(t,x), \frac{\partial u(t,x)}{\partial x}, \frac{\partial^2 u(t,x)}{\partial x^2})$ is a vector containing the derivatives on the grid point x, $\textbf{e}$ is the error of our approximation and $\textbf{M}$ is
$$\textbf{M} = \left( \begin{matrix} 
1 & \Delta x & \Delta x^2 / 2 \\
1 & 0 & 0 \\
1 & -\Delta x & \Delta x^2/2
\end{matrix} \right)$$

Solving the matrix equation for $\textbf{d}$, we have
$$ \textbf{d} = \textbf{M}^{-1} \cdot \textbf{f} - \textbf{M}^{-1} \cdot \textbf{e} $$
where the matrix $\textbf{M}^{-1}$  is given by
$$\textbf{M}^{-1} = \left( \begin{matrix} 
0 & 1 & 0 \\
1/(2\Delta x) & 0 & -1/(2\Delta x) \\
1/\Delta x^2 & -2/\Delta x^2 & 1/\Delta x^2
\end{matrix} \right)$$

The 3rd line of this equation gives us an expression for the 2nd derivative that we can plug into our ODE system
$$\frac{\partial^2 u(t,x)}{\partial x^2} = \frac{u(t,x+\Delta x) - 2u(t,x) +u(t, x-\Delta x))}{\Delta x^2}+\mathcal{O}(\Delta x^2)$$

Doing this, we get
$$\left\{ \begin{array}{@{}l@{}} \frac{\partial u(t,x)}{\partial t} = \Pi (t,x) \\ \frac{\partial \Pi(t,x)}{\partial t} = \frac{c^2}{ \Delta x^2} \left(u(t,x+\Delta x) - 2 u(t,x) + u(t,x-\Delta x)\right) \end{array} \right.\, $$


Now, we can use the [[Runge-Kutta 4 Method]] to solve the ODE system.