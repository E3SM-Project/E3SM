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
| FluxLayerThickEdge | layer thickness used for fluxes through edges. May be centered, upwinded, or a combination of the two
| MeanLayerThickEdge | layer thickness averaged from cell center to edges
| RelVortVertex | curl of horizontal velocity, defined at vertices
| NormRelVortVertex | curl of horizontal velocity divided by layer thickness
| NormPlanetVortVertex | earth's rotational rate (Coriolis parameter, f) divided by layer thickness
| NormRelVortEdge | curl of horizontal velocity divided by layer thickness, averaged from vertices to edges
| NormPlanetVortEdge | earth's rotational rate (Coriolis parameter, f) divided by layer thickness, averaged from vertices to edges
| VelDel2Edge | laplacian of horizontal velocity on edges
| VelDel2DivCell | divergence of laplacian of horizontal velocity on cells
| VelDel2RelVortVertex | curl of laplacian of horizontal velocity on cells
| HTracersEdge | thickness-weighted tracers used for fluxes through edges. May be centered, upwinded or a combination of the two
| Del2TracersCell | laplacian of tracers on cells
