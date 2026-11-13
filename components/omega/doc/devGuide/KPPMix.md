(omega-dev-kppmix)=

# KPP Boundary Layer Mixing

This guide maps Omega KPP implementation details to runtime behavior and code
locations. It complements the design page by focusing on concrete APIs,
call flow, and developer test strategy.

## Implementation Overview

The Omega implementation of KPP mostly follows directly from the MPAS-Ocean implementation.
Notably, by default it does not match diffusivity and viscosity from interior mixing
(below the ocean surface boundary layer) sources at the base of the boundary
layer. Instead these separate sources are added directly to the KPP diagnosed
diffusivity and viscosity.  Matching can be enabled via the `MatchTechnique` parameter to
`MatchBoth` in the Omega yaml file.  Boundary layer depth is computed as the depth
where the bulk Richardson number exceeds a critical value.  It is then smoothed horizontally.
KPP diffusivity, viscosity, and non local fluxes are computed based on this boundary layer.  Unlike MPAS-Ocean, the boundary layer depth, vertical vicosity and vertical diffusivity are calculated at the beginning of a time step to ensure consistency of the nonlocal and local parts of the KPP scheme.

Omega KPP is implemented in `KPPMix` as a singleton with two major compute
phases:

1. OSBL depth diagnosis (`computeOSBLDepth`)
2. Coefficient/profile construction (`computeMixingCoefficients`)

Main class/API surface is in `src/ocn/KPPMix.h` and implementation is in
`src/ocn/KPPMix.cpp`.

### Depth convention

All KPP depths are measured downward from the free surface, not from the geoid.
`VertCoord::GeomZInterface` and `GeomZMid` are geometric heights relative to
`z = 0`, with `GeomZInterface(ICell, MinLayerCell)` equal to
`VertCoord::SshCell`. KPP therefore forms depths as
`SshCell(ICell) - GeomZ...(ICell, K)`. Layer thicknesses are differences of
geometric heights and are unaffected by the sea surface height. After the
Richardson-number search and any horizontal smoothing, the OSBL depth is clamped
to be no shallower than the bottom of the first active layer.

## KPP Notation

KPP splits the water column at the ocean surface boundary layer (OSBL) depth `h`. Inside the OSBL, diffusivity and viscosity are prescribed via a cubic shape function
scaled by `h` and a turbulent velocity scale; below it, only interior mixing
 applies (e.g., shear instability driven mixing).

| Symbol | Code name | Units | Meaning |
| --- | --- | --- | --- |
| `h` | `HOSBL`, `OSBLDepth` | m | OSBL depth below the free surface |
| `d` | `ZDepth`, `ZCenter` | m | depth below the free surface |
| `sigma` | `Sigma` | - | normalized depth in `[-1,0]`, Omega sign convention: 0 at the surface, -1 at the OBL base |
| `sigma_mu` | `SigmaMu` | - | `-Sigma`, in `[0,1]`; the CVMix/Large et al. convention |
| `u*` | `UStar` | m/s | surface friction velocity from wind stress |
| `B_0` | `BuoyFlux` | m^2/s^3 | surface buoyancy flux; **negative is destabilizing** (convection) |
| `L` | `LMoninObukhov` | m | Monin-Obukhov length, `u*^3 / (kappa B_0)` |
| `zeta` | `Zeta` | - | stability coordinate `d/L`; negative is unstable |
| `kappa` | `VonKar` | - | von Karman constant |
| `epsilon` | `SurfaceLayerExtent` | - | surface layer as a fraction of `h` |
| `w_m`, `w_s` | `WMTurb`, `WSTurb` | m/s | turbulent velocity scales for momentum and scalars |
| `Vt^2` | `Vt2`, `UnresolvedShear` | m^2/s^2 | unresolved turbulent shear, Large et al. (1994) Eq. 23 |
| `Ri_b` | `RiBulk`, `BulkRichardsonNumber` | - | bulk Richardson number, Eq. 21 |
| `G(sigma)` | `kppShape*` | - | non-dimensional profile shapes |
| `gamma` | `VertNonLocalFlux` | - | non-local (counter-gradient) tracer flux coefficient |

The OSBL depth is the shallowest `d` at which `Ri_b(d)` reaches the critical
value; the search is done per cell in `computeOSBLDepth`, and the crossing depth
is refined by a quadratic fit through the three nearest cell-center `Ri_b`
values.

Two properties of that search are worth knowing before editing it. The
surface-layer reference averages span the top `epsilon*d`, which grows
monotonically with the trial depth, so the running sums are carried across trial
depths in a two-pointer scan rather than rebuilt from the surface at each one;
resetting them inside the loop would make the search `O(K^2)`. The loop also
stops at the crossing, since only the crossing level and the two above it feed
the refinement. `DebugDiagnostics` suppresses that early exit so the `Ri_b`
diagnostic profiles are filled over the whole column; it has no effect on
`OSBLDepth` or on any mixing coefficient.

### Shape and stability functions

The non-dimensional functions live in `src/ocn/KPPConstants.h` in namespace
`OMEGA::KPP`. All shape functions take `Sigma` in `[-1,0]` and convert to
`SigmaMu` internally, so callers never flip signs:

- `kppShapeMomentum`, `kppShapeScalar` -- gradient shapes for viscosity and
  diffusivity (`SimpleShapes`)
- `kppShapeMatched` -- gradient shape that additionally reaches a prescribed
  value at the OBL base, used by `MatchBoth`
- the non-local flux has no shape function of its own; it reuses whichever
  scalar shape is in effect, scaled by `C_s` instead of `h*w_s`
- `kppPhiInvMomentum`, `kppPhiInvScalar` -- **inverse** Monin-Obukhov stability
  functions `phi^-1(zeta)`; they already return the reciprocal, so do not
  invert them again at the call site

### Constants and defaults

`src/ocn/KPPConstants.h` is the single authoritative source for KPP default
values. The runtime-configurable members of `KPPMix` (`SurfaceLayerExtent`, the
two ice-fraction thresholds and `MinimumOSBLUnderSeaIce`) are initialized from
those constants rather than from inline literals, so a default is changed in
exactly one place.

Per-thread edge scratch arrays in `computeOSBLDepth` are sized from
`HorzMesh::MaxEdgesBound`, the shared compile-time bound on edges per cell;
KPP does not define its own maximum.

## Runtime Call Flow

### Tendency coupling

KPP coupling into tendencies occurs through:

- `KPPMix::computeKPPFields(...)`

`computeKPPFields(...)` assembles required inputs:

- the current tracer and normal-velocity state views
- temperature and salinity tracer indices
- the current EOS and surface forcing instances
- the configured non-local tracer-flux policy

It then calls `KPPMix::update(...)`. KPPMix owns preparation of its derived
state: EOS fields, surface-referenced potential density, edge tangential
velocity, surface friction velocity, surface buoyancy flux, and temporary ice
fraction. Temperature and salinity are zero-copy subviews of the canonical
tracer array, while surface tracer fluxes remain owned by the forcing state.
`KPPMix::update(...)` validates required non-local tracer fluxes and performs
the OSBL-depth and mixing-coefficient calculation.

```c++
KPPInstance->update(...)
```

### Time stepper interaction

KPP is evaluated exactly once per time step by every Omega time stepper
(`src/timeStepping/RungeKutta4Stepper.cpp`,
`src/timeStepping/RungeKutta2Stepper.cpp`,
`src/timeStepping/ForwardBackwardStepper.cpp`) through two hooks in
`src/timeStepping/TimeStepper.cpp`:

1. **Start-of-step compute**: `TimeStepper::updateKPPFields(...)` fetches the
   current tracer array and calls `KPPMix::computeKPPFields(...)`. Each
   stepper calls it once at the top of `doStep()`, after `prescribeState` /
   `prescribeVelocity` and before the first tendency evaluation:

   - **RungeKutta4Stepper**: in the `Stage == 0` branch, before the first
     `computeAllTendencies(...)`.
   - **RungeKutta2Stepper**: after the initial `prescribeState(...)`.
   - **ForwardBackwardStepper**: after `prescribeVelocity(...)`.

2. **End-of-step application**:
   `TimeStepper::applyImplicitVerticalMixing(...)` is called by every
   stepper's `doStep()` immediately after `State->updateTimeLevels()`. It
   recomputes auxiliary state and calls `VertMix::VertMixImplicit(...)` if
   enabled. It does **not** recompute KPP; the fields from the start of the
   step are reused.

The KPP fields are consumed in two places within the step:

- `Tendencies::computeTracerTendenciesOnly(...)` adds the non-local tracer
  tendency from `KPPMix::VertNonLocalFlux` at every stage.
- `VertMix::computeVertMix(...)`, called from `VertMixImplicit(...)`, merges
  `KPPMix::VertDiff` / `VertVisc` into the final coefficients using
  `KPPMix::OSBLDepthIndex`.

Both therefore see the same boundary-layer depth and the same shape function,
so the non-local flux and the diffusivity it is paired with cannot become
inconsistent. KPP diagnostics written for a step describe the state at the
beginning of that step.

Note that `VertMix::VertMixImplicit(...)` recomputes Brunt-Vaisala frequency
itself before calling `computeVertMix(...)`, so shear and convective mixing
still use end-of-step stratification; only the KPP contribution is lagged.

## Configuration Mapping

KPP reads configuration from the `VertMix: KPP` subgroup during `KPPMix::init`.
Important keys and class members:

- `Enable` -> `Enabled`
- `CriticalBulkRichardsonNumber` -> `CriticalRichardson`
- `SurfaceLayerExtent` -> `SurfaceLayerExtent`
- `MatchTechnique` -> `MatchTechnique` (a `KPPMatchType` enum, not a string)
- `InterpType2` -> `InterpType2Str`
- `UseEnhancedDiffusion` -> `UseEnhancedDiffusion`
- `UseLangmuirTurbulence` -> `UseLangmuirTurbulence`
- `IceFractionThresholdForLangmuir` -> `IceFractionThresholdForLangmuir`
- `IceFractionThresholdForMinimumOSBL` -> `IceFractionThresholdForMinimumOSBL`
- `MinimumOSBLUnderSeaIce` -> `MinimumOSBLUnderSeaIce`
- `DebugDiagnostics` -> `DebugDiagnostics`

See [User KPP guide](../userGuide/KPPMix.md) for defaults and runnable examples.

## Output and Diagnostic Fields

`KPPMix::defineFields()` registers KPP outputs for I/O. Frequently used outputs:

- `OSBLDepth`
- `VertNonLocalFlux`
- `BulkRichardsonNumber`
- `BulkRichardsonShear`
- `UnresolvedShear`
- `BuoyancyJump`
- `TurbulentVelocityScale`
- `PotentialDensity`
- `SurfaceFrictionVelocity`
- `SurfaceBuoyancyFlux`

These can be enabled in output stream contents to diagnose OSBL and profile
behavior in experiments.

## Developer Notes

- `MatchTechnique` accepts only `SimpleShapes` and `MatchBoth`. Any other value
  aborts at init rather than falling back, so a typo cannot silently change the
  scheme being run.
- `MatchBoth` matches the interior coefficient *value* at the OBL base. The
  shape derivative there is zero, so the gradient is not yet matched despite the
  name.
- `MatchBoth` needs interior coefficients to be passed in; without them
  `ShapeAtBase` is zero and it degenerates exactly to `SimpleShapes`.
- `Interptype2` accepts `LMD94`, `Linear`, `Quadratic`, and `Cubic`, but `LMD94` is strongly recommended.  For `MatchBoth` other options can result in negative diffusivities and viscosities.
- When `DebugDiagnostics` is enabled in debug builds, targeted diagnostic
  logging is available; behavior is compile/build-mode aware.

## Testing Strategy

### Code-level checks

1. Verify KPP initialization with explicit and default YAML keys.
2. Verify the `computeKPPFields` path executes with KPP enabled and
   early-returns when disabled.
3. Verify per-stepper sequencing: exactly one KPP evaluation per step for all
   three steppers, occurring before the first tendency evaluation and reused
   by the end-of-step implicit vertical mixing. The `KPP:computeKPPFields`
   Pacer region can be used to confirm the call count.
4. Run the KPP single column tests, comparing to prior baselines, MPAS-Ocean,
   and analytic solutions where available

### Diagnostics-based checks

1. Output `OSBLDepth`, `BulkRichardsonNumber`, and `VertDiff`.
2. Confirm OSBL depth and coefficient evolution under changing forcing.
3. Validate optional outputs (`VertNonLocalFlux`) only when enabled.

### Regression checks

Run CTests and appropriate OMEGA regression workflows documented in:
[Testing Code](./Testing.md).
