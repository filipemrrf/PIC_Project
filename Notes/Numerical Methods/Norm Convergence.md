To do a [[Norm Convergence]], we define the norm $L^2$ of the distance between our solutions $\Delta u(t,x)$, as
$$||\Delta u||(t) = \sqrt{\int_{x_i}^{x_f} (\Delta u(t,x))^2 \,dx}$$

However, as our solutions are discrete, we have to write our integral as a sum, getting
$$||\Delta u||(t) = \sqrt{\int_{x_i}^{x_f} (\Delta u(t,x))^2 \,dx} \approx \sqrt{\Delta x \sum_{x} (\Delta u(t,x))^2}$$

Even though this expression is true for cartesian coordinates, if we are working in spherical coordinates with spherical symmetry, we need to do a change in coordinates in our expression, getting
$$||\Delta u||(t) = \sqrt{4 \pi \int_{x_i}^{x_f} r^2(\Delta u(t,r))^2 \,dr} \approx \sqrt{4 \pi \Delta r \sum_{r} r^2(\Delta u(t,r))^2}$$

By comparing the norm of the distance between solutions, we can tell how quick our method converges.