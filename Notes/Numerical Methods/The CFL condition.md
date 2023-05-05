The Courant–Friedrichs–Lewy (CFL) condition is a necessary condition for the stability of explicit finite difference schemes for the numerical solution of partial differential equations (PDEs). It is named after Richard Courant, Kurt Friedrichs, and Hans Lewy, who first derived it in the 1920s.

This condition relates the length of the time step to a function of the interval lengths of each spatial coordinate and of the maximum speed that information can travel in the physical space. If a point in our solution travels more than one space step in one time step, our method may become unstable. To solve this we need to either increase our space step our decrease our time step.

This condition can be expressed as
$$CFL = c \frac{\Delta t}{\Delta x} \leq 1$$

![[CFL Condition.png]]