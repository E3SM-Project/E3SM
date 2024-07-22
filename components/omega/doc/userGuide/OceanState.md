(omega-user-state)=

## Ocean State

The `OceanState` class provides a container for the non-tracer prognostic variables in OMEGA, namely `normalVelocity` and `layerThickness`.
Upon creation of a `OceanState` instance, these variables are allocated and registered with the IO infrastructure.
The class contains a method to update the time levels for the state variables between timesteps.
This involves a halo update, pointer swap, and updating the `IOFields` data references.
