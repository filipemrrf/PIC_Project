To do a [[Pointwise Convergence]] analysis of our solutions, we compare 2 different solutions at the same points in space and time, using an intermediate solution. This gives us a function of space and time that tells us how close our solutions are to each other at a given point $(t,x)$.

If we have three solutions $f_1 (t,x)$, $f_2 (t,x)$ and $f_3 (t,x)$ with increasing resolution, we can do a pointwise convergence by comparing the functions $g_{12} (t,x) = f_2 (t,x) - f_1 (t,x)$ and $g_{23} (t,x) = f_3 (t,x) - f_2 (t,x)$. If, for example, $g_{2,3}(t,x) \approx 2^n \: g_{12}(t,x)$, our method converges in $n$th order.

This analysis also makes it easy to identify mistakes in the implementation of our method (or in the method itself), as we can see if there is a particular region of our solutions that is diverging.