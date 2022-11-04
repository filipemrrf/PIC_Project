Our objective is to numerically solve the wave equation with Cauchy boundary conditions. That is
$$\frac{\partial^2 u(t,x)}{\partial t^2} = c^2 \frac{\partial^2 u(t,x)}{\partial x^2}$$
$$u(0,x)=f(x)$$
$$\frac{\partial u(0,x)}{\partial t}=g(x)$$
To do that, we will use the finite difference technique.


Using the Tayor series to approximate small displacements in $x$ while keeping time contant, we get
$$u(t, x+\Delta x) = u(t,x) + \frac{\partial u(t,x)}{\partial x}\Delta x + \frac{1}{2}\frac{\partial^2u(t,x)}{\partial x^2} \Delta x^2+O(\Delta x^3)$$
$$u(t,x-\Delta x) = u(t,x) - \frac{\partial u(t,x)}{\partial x}\Delta x + \frac{1}{2}\frac{\partial^2u(t,x)}{\partial x^2} \Delta x^2-O(\Delta x^3)$$
By adding both expressions, we get
$$u(t,x+\Delta x)+u(t,x-\Delta x) = 2u(t,x)+\frac{\partial^2u(t,x)}{\partial x^2}\Delta x^2 \Leftrightarrow$$
$$\Leftrightarrow \frac{\partial^2u(t,x)}{\partial x^2} = \frac{u(t,x+\Delta x)-2u(t,x)+u(t,x-\Delta x)}{\Delta x^2}$$

Doing the same procedure for the time derivative (keeping the position constant), we get
$$\frac{\partial^2u(t,x)}{\partial t^2} = \frac{u(t+\Delta t,x)-2u(t,x)+u(t-\Delta t,x)}{\Delta t^2}$$

This way, we can writte the wave equation as
$$\frac{\partial^2 u}{\partial t^2} = c^2 \frac{\partial^2 u}{\partial x^2} \Leftrightarrow$$
$$\Leftrightarrow \frac{u(t+\Delta t,x)-2u(t,x)+u(t-\Delta t,x)}{\Delta t^2} = c^2 \frac{u(t,x+\Delta x)-2u(t,x)+u(t,x-\Delta x)}{\Delta x^2} \Leftrightarrow$$
$$\Leftrightarrow u(t+\Delta t,x) = 2u(t,x)-u(t-\Delta t,x)+c^2\frac{\Delta t^2}{\Delta x^2}(u(t, x+\Delta x)-2u(t,x)+u(t, x-\Delta x))$$


To solve the equation, we also need to discretize the boundary contition for the derivative. By expanding the first derivatives in a Taylor series
$$u(t+\Delta t,x) = u(t,x) + \frac{\partial u(t,x)}{\partial t}\Delta t+O(\Delta t^2) \approx u(t,x) + \frac{\partial u(t,x)}{\partial t}\Delta t \Leftrightarrow$$
$$\Leftrightarrow \frac{\partial u(t,x)}{\partial t} = \frac{u(t+\Delta t,x)-u(t,x)}{\Delta t}$$
By setting the boundary condition $\frac{\partial u(0,x)}{\partial t} = g(x)$, we get
$$g(x) = \frac{u(\Delta t,x)-u(0,x)}{\Delta t} \Leftrightarrow $$
$$\Leftrightarrow u(\Delta t,x) = u(0,x)+g(x)\Delta t$$