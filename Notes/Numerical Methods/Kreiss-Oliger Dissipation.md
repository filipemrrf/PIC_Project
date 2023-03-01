When working with non-linear systems of equations like the equations of general relativity, many stable numerical methods may become slightly unstable because of effects of coefficients that depend on dynamical variables and because of lower order terms.

In order to solve these problems, we can add dissipative terms to the finite difference operators that act as "filters" that preferentially damp modes with wavelength similar to the grid spacing. These high frequency modes are already unresolved and are therefore severely affected by truncation errors, and they are also frequently the source of instabilities. Because of that, we are often better off dissipating them away as otherwise they might become unstable and ruin the simulation.

The standard way of adding artificial dissipation to a finite difference approximation is known as [[Kreiss-Oliger Dissipation]]. 

The general formula for this kind of dissipation for a scheme with accuracy $2r-2$ is:
$$Q = \sigma (-1)^{r-1} h^{2r-1}(D_+)^r \rho (D_-)^r /2^{2r},$$
where $D_+$ and $D_-$ are the forward and backward finite difference operators respectively,  $\sigma$ regulates the strength of the dissipation and $\rho$ is a weighting function typically set to $1$ in the interior but may go to $0$ at the boundary.


When our code is 2nd order accurate, the dissipation operator is written as:
$$Q = -\sigma h^3 (D_+)^2 \rho (D_-)^2/16$$

Applying this operator to a vector $u$ and using $\rho$ as a constant, we get
$$Q u = (-\sigma h^3 (D_+)^2 \rho (D_-)^2/16)u = $$
$$= - \frac{\sigma h^3 \rho}{16} (D_+)^2 (D_-)^2 u =$$
$$= - \frac{\sigma h^2 \rho}{16} (D_+)^2 D_-(u_{i} - u_{i-1}) = $$
$$= - \frac{\sigma h \rho}{16} (D_+)^2 (u_i - 2 u_{i-1} + u_{i-2}) =$$
$$= - \frac{\sigma \rho}{16} D_+ (u_{i+1}-3u_i+3u_{i-1}-u_{i-2}) =$$
$$=- \frac{\sigma \rho}{16h} (u_{i+2}-4u_{i+1}+6u_i-4u_{i-1}+u_{i-2})$$


For a 4th order accurate code, our operator takes the form:
$$Q = \sigma h^5(D_+)^3 \rho (D_-)^3 /64$$

Applying this operator to a vector $u$ like before, we get
$$Qu = (\sigma h^5(D_+)^3 \rho (D_-)^3 /64)u =$$
$$= - \frac{\sigma h^5 \rho}{64} (D_+)^3 (D_-)^3 u =$$
$$= - \frac{\sigma h^4 \rho}{64} (D_+)^3 (D_-)^2 (u_{i} - u_{i-1}) =$$
$$= - \frac{\sigma h^3 \rho}{64} (D_+)^3 D_- (u_i - 2 u_{i-1} + u_{i-2}) =$$
$$= - \frac{\sigma h^2 \rho}{64} (D_+)^3(u_i-3u_{i-1}+3u_{i-2}-u_{i-3}) =$$
$$= - \frac{\sigma h \rho}{64} (D_+)^2 (u_{i+1}-4u_{i}+6u_{i-1}-4u_{i-2}+u_{i-3}) =$$
$$= - \frac{\sigma \rho}{64} D_+ (u_{i+2} - 5 u_{i+1}+10u_i-10u_{i-1}+5u_{i-2} -u_{i-3}) =$$
$$= - \frac{\sigma \rho}{64h}(u_{i+3}-6u_{i+2}+15u_{i+1}-20u_i+15u_{i-1}-6u_{i-2}+u_{i-3})$$