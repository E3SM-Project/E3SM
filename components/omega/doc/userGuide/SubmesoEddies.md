(omega-user-submeso-eddies)=

# Submesoscale Eddy Parameterization

Omega includes an optional submesoscale mixed layer instability (MLI)
parameterization through the `SubmesoEddies` class. When enabled, the model
computes an eddy-induced transport velocity on edges and adds it to the normal
transport velocity used by thickness and tracer advection to represent the transport due to the unresolved submesoscale eddies.

The current implementation follows the Fox-Kemper et al. (2011) (FK11) closure
described in the
{ref}`omega-design-submesoscale-eddies` design document.

## Configuration

Configure the parameterization in the `Submeso` section of the YAML input:

```yaml
Submeso:
  Enable: false
  Tau: 172800.0
  Ce: 0.08
  LfMin: 1.0e3
  DsMax: 100.0e3
```

- `Enable`: turns the parameterization on/off.
- `Tau`: MLI timescale parameter (s), used to prevent division by zero at the equator.
- `Ce`: nondimensional efficiency coefficient.
- `LfMin`: minimum frontal width limiter (m).
- `DsMax`: maximum grid-length limiter used in the closure (m).

## Diagnostics

When enabled, the following fields are available in the `Submeso` field group:

- `GradBuoyEdgeInterface` (s^-2): buoyancy gradient on edge interfaces.
- `EddyVelocity` (m/s): eddy-induced transport velocity.
