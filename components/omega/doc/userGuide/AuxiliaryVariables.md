(omega-user-aux-vars)=

# Auxiliary Variables

Omega needs to compute a couple of auxiliary variables at every model time step
in order to advance the model state. Some auxiliary variables can be computed in
different ways, and the user can specify that in the input configuration file.
For example, the value of thickness used in advective flux can be centered or upwind.
```yaml
    Advection:
       FluxThicknessType: 'Center'
```
Auxiliary variables are also available for output.

The following auxiliary variables are currently available:
| Name | Description |
| ------------------- | ------- |
| KineticEnergyCell | kinetic energy of horizontal velocity on cells
| VelocityDivCell | divergence of horizontal velocity
| FluxPseudoThickEdge | pseudo-thickness used for fluxes through edges. May be centered, upwinded, or a combination of the two
| MeanPseudoThickEdge | pseudo-thickness averaged from cell center to edges
| RelVortVertex | curl of horizontal velocity, defined at vertices
| NormRelVortVertex | curl of horizontal velocity divided by pseudo-thickness
| NormPlanetVortVertex | earth's rotational rate (Coriolis parameter, f) divided by pseudo-thickness
| NormRelVortEdge | curl of horizontal velocity divided by pseudo-thickness, averaged from vertices to edges
| NormPlanetVortEdge | earth's rotational rate (Coriolis parameter, f) divided by pseudo-thickness, averaged from vertices to edges
| VelDel2Edge | laplacian of horizontal velocity on edges
| VelDel2DivCell | divergence of laplacian of horizontal velocity on cells
| VelDel2RelVortVertex | curl of laplacian of horizontal velocity on cells
| HTracersEdge | thickness-weighted tracers used for fluxes through edges. May be centered, upwinded or a combination of the two
| Del2TracersCell | laplacian of tracers on cells
| SurfTracerRestoringDiffsCell | surface tracer restoring differences on cells
| TracersMonthlySurfClimoCell | monthly climatology values to restore to for surface tracer on cells
| DenMixLayerDepth | density-threshold mixed-layer depth on cells


## Kinetic energy on cells

In [Ringler et al. (2010)](https://www.sciencedirect.com/science/article/pii/S0021999109006780), the cell-centered kinetic energy ($K_i$) is a geometry-weighted combination of the squared edge-normal velocities surrounding cell ($i$), constructed so that its discrete gradient enters the vector-invariant momentum equation as part of the Bernoulli-gradient term, ($-\nabla(K_i+\Phi_i)$). Although only one velocity component is stored at each C-grid edge, the differently oriented edges collectively represent the two-dimensional velocity magnitude; for uniform flow on an isotropic cell, the construction recovers ($K_i=\tfrac12|\mathbf{u}|^2$). Importantly, this definition is chosen for algebraic compatibility with the edge-based kinetic-energy norm and the discrete continuity equation, enabling the nondissipative momentum terms to conserve total energy to within time-discretization error rather than providing an arbitrary pointwise reconstruction of the velocity magnitude.

It should be noted that this is the irrotational kinetic energy and is not the total kinetic energy. In principle, this implies that that $K_i$ could be an underestimate of the total kinetic energy. Even for a regular hexagon in an irrotational constant flow, this formulation underestimates the total kinetic energy by 13-16\% depeneding on the flow orientation. However, it has been shown with MPAS-Ocean at standard resolution (Icos30) that $K_i$ *exceeds* $K = 0.5 (u_{cell}^2 + v_{cell}^2)$, likely due to velocity noise at the grid scale that is filtered in the course of reconstructing velocities at cell-centers. This offers a justification for employing $K_i$ in other terms of the momentum equation such as bottom drag.
