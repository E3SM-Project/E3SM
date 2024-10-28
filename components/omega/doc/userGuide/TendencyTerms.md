(omega-user-tend-terms)=

# Tendency Terms

Tendencies for updating state variables are computed using functors, which are
class objects that can be called like functions. These functors are used within
Kokkos parallel loops to compute the full tendency arrays. The following
tendency terms are currently implemented:
| Name | Description |
| ---------------- | ---------------- |
| ThicknessFluxDivOnCell | divergence of layer thickness flux, defined at cell center
| PotentialVortHAdvOnEdge | horizontal advection of potential vorticity, defined on edges
| KEGradOnEdge | gradient of kinetic energy, defined on edges
| SSHGradOnEdge | gradient of sea-surface height, multiplied by gravitational acceleration, defined on edges
| VelocityDiffusionOnEdge | Laplacian horizontal mixing, defined on edges
| VelocityHyperDiffOnEdge | biharmonic horizontal mixing, defined on edges
| TracerHorzAdvOnCell | horizontal advection of thickness-weighted tracers
| TracerDiffOnCell | horizontal diffusion of thickness-weighted tracers
| TracerHyperDiffOnCell | biharmonic horizontal mixing of thickness-weighted tracers

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
| ThicknessFluxDivOnCell | ThicknessFluxTendencyEnable | enable/disable term
| PotentialVortHAdvOnEdge | PVTendencyEnable | enable/disable term
| KEGradOnEdge | KETendencyEnable | enable/disable term
| SSHGradONEdge | SSHTendencyEnable | enable/disable term
| VelocityDiffusionOnEdge | VelDiffTendencyEnable | enable/disable term
| | ViscDel2 | horizontal viscosity
| VelocityHyperDiffOnEdge | VelHyperDiffTendencyEnable | enable/disable term
| | ViscDel4 | horizontal biharmonic mixing coefficient for normal velocity
| TracerHorzAdvOnCell | TracerHorzAdvTendencyEnable | enable/disable term
| TracerDiffOnCell | TracerDiffTendencyEnable | enable/disable term
| | EddyDiff2 | horizontal diffusion coefficient
| TracerHyperDiffOnCell | TracerHyperDiffTendencyEnable | enable/disable term
| | EddyDiff4 | biharmonic horizontal mixing coeffienct for tracers
