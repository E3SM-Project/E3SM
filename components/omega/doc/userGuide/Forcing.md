(omega-user-forcing)=

# Forcing

This page documents the user-facing configuration and behavior for current forcing in Omega:

- Wind forcing
- Surface tracer restoring

## Wind forcing

Wind forcing adds momentum tendency from surface wind stress.

### Wind forcing configuration

Wind forcing behavior is controlled by two configuration blocks:

```yaml
Omega:
  WindStress:
    InterpType: Isotropic

  Tendencies:
    WindForcingTendencyEnable: true
```

- `WindStress.InterpType`
  - `Isotropic`: isotropic cell-to-edge interpolation for wind stress
  - `Anisotropic`: anisotropic interpolation option
- `Tendencies.WindForcingTendencyEnable`: switch to enable wind forcing tendency

### Required input fields

Wind forcing uses auxiliary wind-stress fields:

- `WindStressZonal`
- `WindStressMeridional`

These are used to form edge-normal stress (`NormalStressEdge`) that enters
momentum tendencies.

## Surface tracer restoring

Surface tracer restoring applies a piston-velocity tendency, or damping, at the ocean
surface for selected tracers. This is implemented to mitigate drifts in chosen tracers
(most often salinity) by nudging the model's simulated tracer values towards observed climatological values.
This process prevents oceanic regimes from shifting away from reality due to errors in surface freshwater
forcing (in the case of salinity restoring). Currently, it is applied everywhere when enabled.

### Surface tracer restoring configuration

Surface tracer restoring is controlled by two configuration blocks:

```yaml
Omega:
  SurfaceRestoring:
    TracersToRestore: [Temperature, Salinity]
    PistonVelocity: 1.585e-5
    MaxDiff: 100.0

  Tendencies:
    SurfaceTracerRestoringEnable: true
```

- `TracersToRestore`: list of tracer names that restoring is applied to
- `PistonVelocity`: restoring rate coefficient
- `MaxDiff`: cap on restoring difference magnitude
- `SurfaceTracerRestoringEnable`: switch to enable surface tracer restoring

When restoring is enabled, Omega resolves `TracersToRestore` into an internal
list of tracer IDs and applies restoring only to tracers in that list.

### Restoring target fields

Surface restoring uses auxiliary fields:

- `TracersMonthlySurfClimoCell`: restoring target climatological values
- `SurfTracerRestoringDiffsCell`: computed target-minus-state differences

The restoring tendency is computed at the surface layer only and is limited by
`MaxDiff`.

## Notes

- If a tracer is not listed in `TracersToRestore`, no restoring tendency is
  applied to that tracer.
- If surface restoring is enabled but no tracer IDs are available at tendency
  compute-time, Omega aborts with an error.
- `MaxDiff` must be positive. A runtime check will error out if not.
