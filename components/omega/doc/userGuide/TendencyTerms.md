(omega-user-tend-terms)=

# Tendency Terms

Tendencies for updating state variables are computed using functors, which are
class objects that can be called like functions. These functors are used within
Kokkos parallel loops to compute the full tendency arrays. The following
tendency terms are currently implemented:
| Name | Description |
| ---------------- | ---------------- |
| PseudoThicknessFluxDivOnCell | divergence of pseudo-thickness flux, defined at cell center
| PotentialVortHAdvOnEdge | horizontal advection of potential vorticity, defined on edges
| KEGradOnEdge | gradient of kinetic energy, defined on edges
| SSHGradOnEdge | gradient of sea-surface height, multiplied by gravitational acceleration, defined on edges
| VelocityDiffusionOnEdge | Laplacian horizontal mixing, defined on edges
| VelocityHyperDiffOnEdge | biharmonic horizontal mixing, defined on edges
| TracerHorzAdvOnCell | horizontal advection of thickness-weighted tracers
| TracerHighOrderHorzAdvOnCell | second order horizontal advection of thickness-weighted tracers
| TracerDiffOnCell | horizontal diffusion of thickness-weighted tracers
| TracerHyperDiffOnCell | biharmonic horizontal mixing of thickness-weighted tracers
| WindForcingOnEdge | forcing by wind stress, defined on edges
| BottomDragOnEdge | bottom drag, defined on edges
| SurfaceTracerRestoringOnCell | surface tracer restoring, defined on cells

Among the internal data stored by each functor is a `bool` which can enable or
disable the contribution of that particular term to the tendency. These flags
need to be provided in the configuration file:
```yaml
Tendencies:
   ThicknessFluxTendencyEnable: true
```

Some functors have member variables that can be set by the user by the
configuration file. The following are all the user-configurable parameters for
the currently available tendency terms:
| Term | Parameter | Description
| ------------ | ------------ | ------------ |
| PseudoThicknessFluxDivOnCell | ThicknessFluxTendencyEnable | enable/disable term
| PotentialVortHAdvOnEdge | PVTendencyEnable | enable/disable term
| KEGradOnEdge | KETendencyEnable | enable/disable term
| SSHGradONEdge | SSHTendencyEnable | enable/disable term
| VelocityDiffusionOnEdge | VelDiffTendencyEnable | enable/disable term
| | ViscDel2 | horizontal viscosity
| VelocityHyperDiffOnEdge | VelHyperDiffTendencyEnable | enable/disable term
| | ViscDel4 | horizontal biharmonic mixing coefficient for normal velocity
| | DivFactor | scale factor for the divergence term
| TracerHorzAdvOnCell | TracerHorzAdvTendencyEnable | enable/disable term
| | HorzTracerFluxOrder | 1 for standard linear advection
| TracerHighOrderHorzAdvOnCell | TracerHorzAdvTendencyEnable | enable/disable term
| | HorzTracerFluxOrder | 2 for second order advection algorithm
| TracerDiffOnCell | TracerDiffTendencyEnable | enable/disable term
| | EddyDiff2 | horizontal diffusion coefficient
| TracerHyperDiffOnCell | TracerHyperDiffTendencyEnable | enable/disable term
| | EddyDiff4 | biharmonic horizontal mixing coeffienct for tracers
| WindForcingOnEdge | WindForcingTendencyEnable | enable/disable term
| BottomDragOnEdge | BottomDragTendencyEnable | enable/disable term
| | BottomDragCoeff | bottom drag coefficient
| SurfaceTracerRestoringOnCell | SurfaceTracerRestoringEnable | enable/disable term

## Second Order Horizontal Advection Algorithm

The horizontal advection is done independently within each ocean layer
so a two dimensional scheme is required to advect mixing ratios.
The second order horizontal advection scheme is described in the paper

William C. Skamarock and Almut Gassmann,
"Conservative Transport Schemes for Spherical Geodesic Grids: High-Order Flux Operators for ODE-Based Time Integration"
Monthly Weather Review, Vol 139: Issue 9, pp 2962–2975, 2011. doi.org/10.1175/MWR-D-10-05056.1

and only a summary of a few equations from that paper are reproduced here.

The conservative form of the scalar transport equation
within a fluid is given in equation (1) of Skamarock and Gassmann,
$$
 \frac{\partial(\rho\psi)}{\partial t} = -\nabla\cdot\mathbf{V}\rho\psi
$$
where $\rho$ is the fluid density, $\psi$ is the mixing ratio, and $\mathbf{V}$ is
the fluid velocity. This is simplified by assuming that $\rho$ is constant and cancels out.
The advection scheme for a finite volume method involves approximating this equation along an element
edge at a certain time step and multiplying by the time step length to get an approximate scalar transport
from element to element during the time step.
The left hand side of this equation is a time derivative and for Omega this is handled by standard Runge—Kutta methods and
not discussed here but can be found in the Skamarock and Gassmann paper.
The right hand side of this equation is the spatial derivative that needs to be evaluated by the Omega horizontal advection scheme.
Omega ocean uses unstructured Voronoi meshes and a
finite volume scheme with $\psi$ defined on cell centers. A first order approximation
to the spatial derivative of the left hand termacross a single edge of an element would be the difference of $\psi$ at
element centroids divided by the distance from cell centroid to cell centroid.  This is expressed in equation (5) of the paper.
For edge $i$ between elements $i+1/2$ and $i-1/2$,
$$
 \frac{\partial(u\psi_i)}{\partial x} = \frac{1}{\Delta x}[F_{i+1/2}(u\psi) - F_{i-1/2}(u\psi)] + O(\Delta x^2).
$$
Where $u$ is $\mathbf{V}\cdot\mathbf{n}$ where $\mathbf{n}$ is the unit normal along the edge being evaluated and $F$ is the evaluation
of $\psi$ for the elements on each side of the edge being evaluated.
This is used when the TracerHorzAdvOnCell user option is active and HorzTracerFluxOrder is 1.
To make the comparison to higher order methods explicit, this first order method is equivalent to defining $\psi$ as a linear function,
$\psi = c_0 + c_x x + c_y y$, between the two elements sharing the edge $i$ and taking the directed derivative along the edge.


For higher order flux calculations, higher order derivatives of $\psi$ are needed and for the user option of
TracerHorzAdvOnCell along with HorzTracerFluxOrder set to 2, a quadratic approximation of $\psi$ is used,
$$
\psi = c_0 + c_x x + c_y y + c_{xx} x^2 + c_{xy} xy + c_{yy} y^2,
$$
defined for each edge in the mesh. The six coefficients in this quadratic equation are calculated from
a least squares fit of $\psi$ in a neighborhood of elements around an edge. The
least squares fit is defined in equation (12) of the paper,
$$
 \mathbf{f = (P^T P)^{-1}P^T s = Bs}
$$
where the vector $\mathbf{f} = [c_0,c_x,c_y,c_{xx},c_{xy},c_{yy}]$ and $\mathbf{P}$ is an $m\times 6$ matrix with each row,
$(1,x_i,y_i,x_i^2,x_iy_i,y_i^2)$, based on the cell center coordinates where there is one row for each of the $m$ elements in the neighborhood.
This results in a $6\times m$ matrix $\mathbf{B}$. The values of $\mathbf{s}$ are the mixing ratios,
$\mathbf{s}=[\psi_0, \psi_1,...,\psi_m]$, of the elements in the neighborhood at the time step being evaluated.
Note that $\mathbf{B}$ depends on the mesh geometry only, so only needs to
be computed once for each edge during initialization. Since only the second order terms are used to determine the
second order derivatives needed for the higher order flux corrections to the lower order linear flux given above
and the derivatives are taken along a line connecting
cell centers, the amount of data that needs to be stored is minimized and $\mathbf{B}$ reduces to a vector.

The second order derivatives of $\psi$ are then combined as shown in equation (11) of Skamarock and Gassmann to get
a second order horizontal transport algorithm.

The computation of these second order fluxes is complicated by having to do these computations on a sphere.
Figure 2 from the Skamarock and Gassmann paper describs the slight modifications needed to compute the entries
of the $\mathbf{P}$ matrix for a spherical geometry where the distances need to be measured as arc lengths along
the surface.

In practice, the convergence of this higher order algorithm on a sphere is slightly below
second order. For the advection of a $\psi$ in the shape of a cosine bell
being advected around a sphere and compared with the analytic solution, the spatial convergence for grid refinement
is about order 1.7 as shown in  {numref}`tracer-higher-order-convergence`:

```{figure} images/higher_order_tracer_convergence_on_sphere.jpeg
:name:  tracer-higher-order-convergence
:align: center
:width: 600 px
Tracer higer order convergence example of a cosine bell advected on a sphere showing an order 1.71 convergence rate
```

## See Also

Additional information on forcing (currently wind forcing and surface tracer
restoring) is detailed in [](omega-user-forcing).
