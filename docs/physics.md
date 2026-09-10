# The 2-D Euler Equations

The model simulates the 2-D inviscid Euler equations for stratified fluid
dynamics, which are defined as follows:

<img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\frac{\partial}{\partial&space;t}\left[\begin{array}{c}&space;\rho\\&space;\rho&space;u\\&space;\rho&space;w\\&space;\rho\theta&space;\end{array}\right]&plus;\frac{\partial}{\partial&space;x}\left[\begin{array}{c}&space;\rho&space;u\\&space;\rho&space;u^{2}&plus;p\\&space;\rho&space;uw\\&space;\rho&space;u\theta&space;\end{array}\right]&plus;\frac{\partial}{\partial&space;z}\left[\begin{array}{c}&space;\rho&space;w\\&space;\rho&space;wu\\&space;\rho&space;w^{2}&plus;p\\&space;\rho&space;w\theta&space;\end{array}\right]=\left[\begin{array}{c}&space;0\\&space;0\\&space;-\rho&space;g\\&space;0&space;\end{array}\right]" title="\large \frac{\partial}{\partial t}\left[\begin{array}{c} \rho\\ \rho u\\ \rho w\\ \rho\theta \end{array}\right]+\frac{\partial}{\partial x}\left[\begin{array}{c} \rho u\\ \rho u^{2}+p\\ \rho uw\\ \rho u\theta \end{array}\right]+\frac{\partial}{\partial z}\left[\begin{array}{c} \rho w\\ \rho wu\\ \rho w^{2}+p\\ \rho w\theta \end{array}\right]=\left[\begin{array}{c} 0\\ 0\\ -\rho g\\ 0 \end{array}\right]" />

<img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\rho_{H}=-\frac{1}{g}\frac{\partial&space;p}{\partial&space;z}" title="\large \rho_{H}=-\frac{1}{g}\frac{\partial p}{\partial z}" />

where

 - <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\rho" title="\large \rho" /> is density,
 - u, and w are winds in the x-, and z-directions, respectively,
 - <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\theta" title="\large \theta" /> is potential temperature related to temperature, T, by
 - <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\theta=T\left(P_{0}/P\right)^{R_{d}/c_{p}}" title="\large \theta=T\left(P_{0}/P\right)^{R_{d}/c_{p}}" />,
 where <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;P_{0}=10^{5}\,\text{Pa}" title="\large P_{0}=10^{5}\,\text{Pa}" />,
 is the surface pressure,
 - <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;g=9.8\,\mbox{m}\,\mbox{s}^{-2}" title="\large g=9.8\,\mbox{m}\,\mbox{s}^{-2}" /> is acceleration due to gravity,
 - <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;p=C_{0}\left(\rho\theta\right)^{\gamma}" title="\large p=C_{0}\left(\rho\theta\right)^{\gamma}" /> is the pressure as determined by an alternative form of the ideal gas equation of state,
 - <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;C_{0}=R_{d}^{\gamma}p_{0}^{-R_{d}/c_{v}}" title="\large C_{0}=R_{d}^{\gamma}p_{0}^{-R_{d}/c_{v}}" />,
 - <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;R_{d}=287\,\mbox{J}\,\mbox{kg}^{-1}\,\mbox{K}^{-1}" title="\large R_{d}=287\,\mbox{J}\,\mbox{kg}^{-1}\,\mbox{K}^{-1}" /> is the dry gas constant
 - <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\gamma=c_{p}/c_{v}" title="\large \gamma=c_{p}/c_{v}" />,
 - <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;c_{p}=1004\,\mbox{J}\,\mbox{kg}^{-1}\,\mbox{K}^{-1}" title="\large c_{p}=1004\,\mbox{J}\,\mbox{kg}^{-1}\,\mbox{K}^{-1}" /> is specific heat at constant pressure, and
 - <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;c_{v}=717\,\mbox{J}\,\mbox{kg}^{-1}\,\mbox{K}^{-1}" title="\large c_{v}=717\,\mbox{J}\,\mbox{kg}^{-1}\,\mbox{K}^{-1}" /> is specific heat at constant volume.

This can be cast in a more convenient form as:

<img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\frac{\partial\mathbf{q}}{\partial&space;t}&plus;\frac{\partial\mathbf{f}}{\partial&space;x}&plus;\frac{\partial\mathbf{h}}{\partial&space;z}=\mathbf{s}" title="\large \frac{\partial\mathbf{q}}{\partial t}+\frac{\partial\mathbf{f}}{\partial x}+\frac{\partial\mathbf{h}}{\partial z}=\mathbf{s}" />

where a bold font represents a vector quantity.

## Maintaining Hydrostatic Balance

The flows this code simulates are relatively small perturbations off of a “hydrostatic” balance, which balances gravity with a difference in pressure:

<img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\frac{dp}{dz}=-\rho&space;g" title="\large \frac{dp}{dz}=-\rho g" />

Because small violations of this balance lead to significant noise in the vertical momentum, it's best not to try to directly reconstruct this balance but rather to only reconstruct the perturbations. Therefore, hydrostasis is subtracted from the equations to give:

<img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\frac{\partial}{\partial&space;t}\left[\begin{array}{c}&space;\rho^{\prime}\\&space;\rho&space;u\\&space;\rho&space;w\\&space;\left(\rho\theta\right)^{\prime}&space;\end{array}\right]&plus;\frac{\partial}{\partial&space;x}\left[\begin{array}{c}&space;\rho&space;u\\&space;\rho&space;u^{2}&plus;p\\&space;\rho&space;uw\\&space;\rho&space;u\theta&space;\end{array}\right]&plus;\frac{\partial}{\partial&space;z}\left[\begin{array}{c}&space;\rho&space;w\\&space;\rho&space;wu\\&space;\rho&space;w^{2}&plus;p^{\prime}\\&space;\rho&space;w\theta&space;\end{array}\right]=\left[\begin{array}{c}&space;0\\&space;0\\&space;-\rho^{\prime}g\\&space;0&space;\end{array}\right]" title="\large \frac{\partial}{\partial t}\left[\begin{array}{c} \rho^{\prime}\\ \rho u\\ \rho w\\ \left(\rho\theta\right)^{\prime} \end{array}\right]+\frac{\partial}{\partial x}\left[\begin{array}{c} \rho u\\ \rho u^{2}+p\\ \rho uw\\ \rho u\theta \end{array}\right]+\frac{\partial}{\partial z}\left[\begin{array}{c} \rho w\\ \rho wu\\ \rho w^{2}+p^{\prime}\\ \rho w\theta \end{array}\right]=\left[\begin{array}{c} 0\\ 0\\ -\rho^{\prime}g\\ 0 \end{array}\right]" />

where a “prime” quantity represents that variable with the hydrostatic background state subtracted off (not a spatial derivative).

## Dimensional Splitting

This equation is solved using dimensional splitting for simplicity and speed. The equations are split into x- and z-direction solves that are, respectively:

<img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;x:\,\,\,\,\,\,\,\,\,\,\frac{\partial\mathbf{q}}{\partial&space;t}&plus;\frac{\partial\mathbf{f}}{\partial&space;x}=\mathbf{0}" title="\large x:\,\,\,\,\,\,\,\,\,\,\frac{\partial\mathbf{q}}{\partial t}+\frac{\partial\mathbf{f}}{\partial x}=\mathbf{0}" />

<img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;z:\,\,\,\,\,\,\,\,\,\,\frac{\partial\mathbf{q}}{\partial&space;t}&plus;\frac{\partial\mathbf{h}}{\partial&space;x}=\mathbf{s}" title="\large z:\,\,\,\,\,\,\,\,\,\,\frac{\partial\mathbf{q}}{\partial t}+\frac{\partial\mathbf{h}}{\partial x}=\mathbf{s}" />

Each time step, the order in which the dimensions are solved is reversed, giving second-order accuracy overall.

## Finite-Volume Spatial Discretization

A Finite-Volume discretization is used in which the PDE in a given dimension is integrated over a cell domain, <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\Omega_{i}\in\left[x_{i-1/2},x_{i&plus;1/2}\right]" title="\large \Omega_{i}\in\left[x_{i-1/2},x_{i+1/2}\right]" />, where <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;x_{i\pm1/2}=x_{i}\pm\Delta&space;x" title="\large x_{i\pm1/2}=x_{i}\pm\Delta x" />, <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;x_{i}" title="\large x_{i}" /> is the cell center, and <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\Delta&space;x" title="\large \Delta x" /> is the width of the cell. The integration is the same in the z-direction. Using the Gauss divergence theorem, this turns the equation into (using the z-direction as an example):

<img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\frac{\partial\overline{\mathbf{q}}_{i,k}}{\partial&space;t}=-\frac{\mathbf{h}_{i,k&plus;1/2}-\mathbf{h}_{i,k-1/2}}{\Delta&space;z}&plus;\overline{\mathbf{s}}_{i,k}" title="\large \frac{\partial\overline{\mathbf{q}}_{i,k}}{\partial t}=-\frac{\mathbf{h}_{i,k+1/2}-\mathbf{h}_{i,k-1/2}}{\Delta z}+\overline{\mathbf{s}}_{i,k}" />

where <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\overline{\mathbf{q}}_{i,k}" title="\large \overline{\mathbf{q}}_{i,k}" /> and <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\overline{\mathbf{s}}_{i,k}" title="\large \overline{\mathbf{s}}_{i,k}" /> are the cell-average of the fluid state and source term over the cell of index `i,k`.

To compute the update one needs the flux vector at the cell interfaces and the cell-averaged source term. To compute the flux vector at interfaces, fourth-order-accurate polynomial interpolation is used using the four cell averages surrounding the cell interface in question.

## Runge-Kutta Time Integration

So far the PDEs have been discretized in space but are still continuous in time. To integrate in time, we use a simple three-stage, linearly third-order-accurate Runge-Kutta integrator. It is solved as follows:

<img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\mathbf{q}^{\star}=\mathbf{q}^{n}&plus;\frac{\Delta&space;t}{3}RHS\left(\mathbf{q}^{n}\right)" title="\large \mathbf{q}^{\star}=\mathbf{q}^{n}+\frac{\Delta t}{3}RHS\left(\mathbf{q}^{n}\right)" />

<img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\mathbf{q}^{\star\star}=\mathbf{q}^{n}&plus;\frac{\Delta&space;t}{2}RHS\left(\mathbf{q}^{\star}\right)" title="\large \mathbf{q}^{\star\star}=\mathbf{q}^{n}+\frac{\Delta t}{2}RHS\left(\mathbf{q}^{\star}\right)" />

<img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\mathbf{q}^{n&plus;1}=\mathbf{q}^{n}&plus;\Delta&space;tRHS\left(\mathbf{q}^{\star\star}\right)" title="\large \mathbf{q}^{n+1}=\mathbf{q}^{n}+\Delta tRHS\left(\mathbf{q}^{\star\star}\right)" />

When it comes to time step stability, it is assumed a maximum speed of propagation of <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;450\,\text{m}\,\text{s}^{-1}" />, which basically means that the maximum wind speed is assumed to be <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;100\,\text{m}\,\text{s}^{-1}" />, which is a safe assumption.
The CFL value is set to 1.5.

## Hyper-viscosity

The centered fourth-order discretization is unstable for non-linear equations and requires extra dissipation to damp out small-wavelength energy that would otherwise blow up the simulation. This damping is accomplished with a scale-selective fourth-order so-called “hyper”-viscosity that is defined as:

<img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\frac{\partial\mathbf{q}}{\partial&space;t}&plus;\frac{\partial}{\partial&space;x}\left(-\kappa\frac{\partial^{3}\mathbf{q}}{\partial&space;x^{3}}\right)=\mathbf{0}" title="\large \frac{\partial\mathbf{q}}{\partial t}+\frac{\partial}{\partial x}\left(-\kappa\frac{\partial^{3}\mathbf{q}}{\partial x^{3}}\right)=\mathbf{0}" />

and this is also solved with the Finite-Volume method just like above. The hyperviscosity constant is defined as:

<img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\kappa=-\beta\left(\Delta&space;x\right)^{4}2^{-4}\left(\Delta&space;t\right)^{-1}" title="\large \kappa=-\beta\left(\Delta x\right)^{4}2^{-4}\left(\Delta t\right)^{-1}" />

where <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\beta\in\left[0,1\right]" title="\large \beta\in\left[0,1\right]" /> is a user-defined parameter to control the strength of the diffusion, where a higher value gives more diffusion. The parameter <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\beta" title="\large \beta" /> is not sensitive to the grid spacing, and it seems that <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;\beta=0.25" title="\large \beta=0.25" /> is generally enough to get rid of <img src="https://latex.codecogs.com/svg.latex?\dpi{300}&space;\large&space;2\Delta&space;x" title="\large 2\Delta x" /> noise contamination.

