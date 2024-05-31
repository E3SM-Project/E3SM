(omega-dev-ocean-state)=

## Ocean State

The OMEGA `OceanState` class
A state object is created by the `init` method, which assumes that `Decomp` and `HorzMesh` have
already been initialized.
```c++
Err = OMEGA::Decomp::init();
Err = OMEGA::HorzMesh::init();
Err = OMEGA::OceanState::init();
```
The constructor:
```c++
OceanState(const std::string &Name, ///< [in] Name for mesh
           HorzMesh *Mesh,          ///< [in] Horizontal mesh
           Decomp *MeshDecomp,      ///< [in] Decomp for Mesh
           Halo *MeshHalo_,         ///< [in] Halo for Mesh
           const int NVertLevels_,  ///< [in] Number of vertical levels
           const int NTimeLevels_   ///< [in] Number of time levels
);
```
allocates the `normalVelocity` and `layerThickness` arrays for a given number of time levels.
The current time level is then registered with the IO infrastructure.

After initialization, the default mesh object can be retrieved via:
```
OMEGA::OceanState *State = OMEGA::OceanState::getDefault();
```

The `OceanState` is meant to be a container that allows the non-tracer prognostic variables to be
passed to the PDE solver routines for tendency and auxiliary variable calculation. For example
```
void AuxiliaryState::compute(OceanState *State) const {
  ...
}
```

For now, member variables that are host arrays have variable names are appended with an
`H`.  Array variable names not ending in `H` are device arrays.  For a given time level,
host to device array is performed via:
```c++
State->copyToDevice(TimeLevel);
```
and a copy from device to host is performed by:
```c++
State->copyToHost(TimeLevel);
```
Eventually, the host arrays will be eliminated when `IO` and `Halo` are extended to
handle host <-> device transfers.   

A time level update to advance the solution at the end of a model timestep is done by:
```c++
State->updateTimeLevels();
```
This shifts the time levels to the previous index with the last index being shifted to the
highest index to be overwritten in the next timestep. A halo exchange is also performed on
these arrays and the IOFields data pointer is attached to the current time level.

The state arrays are deallocated by the `OceanState::clear()` method, which is
necessary before calling `Kokkos::finalize`.
