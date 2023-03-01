Many ODE's require the calculation of the function's derivative, making a method to calculate numerical derivatives necessary. One example of such an equation is the [[Numerically Solving the Wave Equation in 1+1 Dimensions|Wave Equation]].

To calculate a numerical derivative, we need to use finite differences, discretizing our axis into a grid with N cells with side $\Delta x$, covering the whole axis. This discretization is edge centered.

Expanding $u(x)$ around $x$ and its immediate neighbors ($x+\Delta x$ and $x - \Delta x$), we get 
$$ u(x+\Delta x) = u(x) + \Delta x\frac{\partial u(x)}{\partial x} + \frac{1}{2}\Delta x^2\frac{\partial^2 u(x)}{\partial x^2} + \mathcal{O}(\Delta x^3)$$
$$u(x) = u(x)$$
$$ u(x-\Delta x) = u(x) - \Delta x\frac{\partial u(x)}{\partial x} + \frac{1}{2}\Delta x^2\frac{\partial^2 u(x)}{\partial x^2} + \mathcal{O}(\Delta x^3)$$

We can now write these equations in the matrix equation
$$ \textbf{u} = \textbf{M} \cdot \textbf{d} + \textbf{e} $$

In this equation, $\textbf{u} = (u(x+\Delta x), u(x), u(x-\Delta x))$ is a vector containing the function values, $\textbf{d} = (u(x), \frac{\partial u(x)}{\partial x}, \frac{\partial^2 u(x)}{\partial x^2})$ is a vector containing the derivatives on the grid point x, $\textbf{e}$ is the error of our approximation and $\textbf{M}$ is
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

The 2nd and 3rd line of this equation give us 2nd order accurate expressions for the 1st and 2nd derivatives
$$\frac{\partial u(x)}{\partial x} \approx \frac{u(x+\Delta x) - u(x-\Delta x))}{2\Delta x}$$
$$\frac{\partial^2 u(x)}{\partial x^2} \approx \frac{u(x+\Delta x) - 2u(x) +u(x-\Delta x))}{\Delta x^2}$$

We can make a 4th order accurate one by expanding $u(x)$ around $x$, its immediate neighbors ($x+\Delta x$ and $x - \Delta x$), and the ones next to them ($x + 2\Delta x$ and $x - 2 \Delta x$) up to the 4th order. By doing this, we get 
$$\begin{array}{@{}@{}} u(x+2\Delta x) = u(x) + 2\Delta x\frac{\partial u(x)}{\partial x} + 2\Delta x^2\frac{\partial^2 u(x)}{\partial x^2} + \\ + \frac{4}{3} \Delta x^3 \frac{\partial^3 u(x)}{\partial x^3} + \frac{2}{3} \Delta x^4 \frac{\partial^4 u(x)}{\partial x^4}  + \mathcal{O}(\Delta x^5) \end{array}$$
$$\begin{array}{@{}@{}} u(x+\Delta x) = u(x) + \Delta x\frac{\partial u(x)}{\partial x} + \frac{1}{2}\Delta x^2\frac{\partial^2 u(x)}{\partial x^2} + \\ + \frac{1}{6} \Delta x^3 \frac{\partial^3 u(x)}{\partial x^3} + \frac{1}{24} \Delta x^4 \frac{\partial^4 u(x)}{\partial x^4} +\mathcal{O}(\Delta x^5) \end{array}$$
$$u(x) = u(x)$$
$$\begin{array}{@{}@{}} u(x-\Delta x) = u(x) - \Delta x\frac{\partial u(x)}{\partial x} + \frac{1}{2}\Delta x^2\frac{\partial^2 u(x)}{\partial x^2} - \\ - \frac{1}{6} \Delta x^3 \frac{\partial^3 u(x)}{\partial x^3} + \frac{1}{24} \Delta x^4 \frac{\partial^4 u(x)}{\partial x^4} + \mathcal{O}(\Delta x^5) \end{array}$$
$$\begin{array}{@{}@{}} u(x-2\Delta x) = u(x) - 2\Delta x\frac{\partial u(x)}{\partial x} + 2\Delta x^2\frac{\partial^2 u(x)}{\partial x^2} - \\ - \frac{4}{3} \Delta x^3 \frac{\partial^3 u(x)}{\partial x^3} + \frac{2}{3} \Delta x^4 \frac{\partial^4 u(x)}{\partial x^4} + \mathcal{O}(\Delta x^5) \end{array}$$

We can write this equations in the same matrix equation as before, but in this case $\textbf{u} = (u(x+2\Delta x), u(x+\Delta x), u(x), u(x-\Delta x), u(x-2\Delta x))$, $\textbf{d} = (u(x), \frac{\partial u(x)}{\partial x}, \frac{\partial^2 u(x)}{\partial x^2}, \frac{\partial^3 u(x)}{\partial x^3}, \frac{\partial^4 u(x)}{\partial x^4})$, $\textbf{e}$ is the error of our approximation and $\textbf{M}$ is
$$\textbf{M} = \left( \begin{matrix} 
1 & 2\Delta x & 2\Delta x^2 & 4 \Delta x^3/3 & 2\Delta x^4/3\\
1 & \Delta x & \Delta x^2/2 & \Delta x^3/6 & \Delta x^4/24\\
1 & 0 & 0 & 0 & 0\\
1 & -\Delta x & \Delta x^2/2 & -\Delta x^3/6 & \Delta x^4/24\\
1 & -2\Delta x & 2\Delta x^2 & -4 \Delta x^3/3 & 2\Delta x^4/3
\end{matrix} \right)$$

Inverting the matrix, we get
$$M^{-1} = \left( \begin{matrix} 
0 & 0 & 1 & 0 & 0\\
-1/(12\Delta x) & 2/(3\Delta x) & 0 & -2/(3\Delta x) & 1/(12\Delta x)\\
-1/(12\Delta x^2) & 4/(3\Delta x^2) & -5/(2\Delta x^2) & 4/(3\Delta x^2) & -1/(12\Delta x^2)\\
1/(2\Delta x^3) & -1/\Delta x^3 & 0 & 1/\Delta x^3 & -1/(2\Delta x^3)\\
1/\Delta x^4 & -4/\Delta x^4 & 6/\Delta x^4 & -4/\Delta x^4 & 1/\Delta x^4
\end{matrix} \right),$$
which using the matrix equation gives us a fourth order accurate expression for the second spatial derivative:
$$\frac{\partial^2 u(x)}{\partial x^2} \approx \frac{-u(x+2\Delta x) + 16u(x+\Delta x) - 30u(x) + 16u(x-\Delta x) - u(x-2\Delta x)}{12 \Delta x^2}$$