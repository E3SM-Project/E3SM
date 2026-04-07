(omega-dev-forcing)=

# Forcing

This page describes design and implementation details for forcing-related
pathways in Omega, currently this includes:

- Wind forcing
- Surface tracer restoring

## Wind forcing design

### Wind forcing data flow

1. External fields provide:
   - `WindStressZonal`
   - `WindStressMeridional`
2. Auxiliary-state compute builds `NormalStressEdge` from those fields.
3. Tendency term applies wind-stress forcing to edge-normal velocity tendency.

### Wind forcing key classes/components

- `WindForcingAuxVars`
  - Stores wind-stress cell fields and computed `NormalStressEdge`
  - Applies configured interpolation choice (`InterpType`)
- `AuxiliaryState::computeMomAux`
  - Calls `WindForcingAuxVars::computeVarsOnEdge`
- `WindForcingOnEdge` tendency term
  - Adds contribution proportional to normal stress and inverse layer
    thickness in the surface layer

### Wind forcing config coupling

- `Omega.WindStress.InterpType`
  - mapped to `InterpCellToEdgeOption`
- `Omega.Tendencies.WindForcingTendencyEnable`
  - gates execution of wind forcing tendency kernel

## Surface tracer restoring design

### Surface tracer restoring data flow

1. External fields provide target values: `TracersMonthlySurfClimoCell` (values and units should match the state variables)
2. Auxiliary-state compute forms restoring differences: `SurfTracerRestoringDiffsCell = target - tracer_surface`
3. Tendency term applies restoring only at surface layer and only for tracers selected from `SurfaceRestoring.TracersToRestore`.

### Surface tracer restoring key classes/components

- `SurfTracerRestAuxVars`
  - Inputs: `TracersMonthlySurfClimoCell`, tracer state array
  - Output: `SurfTracerRestoringDiffsCell`
  - Uses `MinLayerCell` to select surface layer index
- `SurfaceTracerRestoringOnCell` tendency term
  - Applies `PistonVelocity * SurfTracerRestoringDiffsCell` at surface
- `Tendencies`
  - Parses `SurfaceRestoring.TracersToRestore` and resolves tracer indices
  - Builds `TracerIdsToRestore` and `NTracersToRestore`
  - Applies tracer-selection logic at call site in
    `computeTracerTendenciesOnly`
  - Aborts if restoring is enabled but no tracer IDs are available

### Surface tracer restoring config coupling

- `Omega.SurfaceRestoring.PistonVelocity`
  - tendency scaling
- `Omega.SurfaceRestoring.TracersToRestore`
  - tracer-level enable list used to build `TracerIdsToRestore`
- `Omega.Tendencies.SurfaceTracerRestoringEnable`
  - gates restoring tendency execution

## Notes

- If a tracer is not listed in `TracersToRestore`, no restoring tendency is applied to that tracer.
- If restoring is enabled but no tracer IDs are available at tendency compute-time, Omega aborts with an error.
- It is assumed that the incoming `TracersMonthlySurfClimoCell` fields (values and units) match the Omega state variables (i.e. conservative temperature and absolute salinity for Teos-10). If not, a pre-processing conversion should be implemented.
- Surface tracer restoring is active everywhere if enabled. A flag to turn it off under sea ice will need to be added in later development if this feature is desired.
- Unlike MPAS-Ocean, a `MaxDiff` clamping is not applied here. This check should instead be implemented in Ocean Validate when that is available.
- A global scaling to ensure zero-sum has not been implemented for the surface tracer restoring, but should be added in later development.
- At this stage, there is no temporal interpolation applied to the restoring targets, the raw `TracersMonthlySurfClimoCell` snapshot is used.
