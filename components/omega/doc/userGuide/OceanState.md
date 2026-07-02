(omega-user-state)=

## Ocean State

The `OceanState` class provides a container for the non-tracer prognostic variables in Omega:
`NormalVelocity` and `PseudoThickness`.
Upon creation of an `OceanState` instance, these variables are allocated and registered with the
IO infrastructure as part of the `State` and `Restart` field groups, so they are read from the
initial-condition file and written to restart files.
The class contains a method to update the time levels for the state variables between timesteps.
This involves a halo update, time level index update, and updating the `IOFields` data references.

The relative pressure at the top of the ocean column, `SurfacePressure`, is owned by the
[vertical coordinate](omega-user-vert-coord) rather than the `OceanState`, since it is the top
boundary condition of the pressure field. It is still registered in the `State` and `Restart`
field groups, so it continues to be written to restart files and read from the initial-condition
and restart files when present; if it is absent, it defaults to zero.
