(omega-kppmix)=

# KPP Boundary Layer Mixing

This page explains how to enable, configure, and use Omega K-Profile
Parameterization (KPP) boundary layer mixing in runs.  The implementation follows
directly from the MPAS-Ocean implementation.

Related pages:
- KPP design/theory: [Design KPP document](../design/KPPMix.md)
- KPP implementation details: [Developer KPP document](../devGuide/KPPMix.md)
- Broader vertical mixing options: [Vertical Mixing Coefficients](./VerticalMixingCoeff.md)

## What KPP Provides

KPP computes:

- Ocean surface boundary layer depth (`OSBLDepth`)
- Vertical viscosity (`VertVisc`)
- Vertical diffusivity (`VertDiff`)
- Non-local tracer flux profile (`VertNonLocalFlux`)

It uses a bulk Richardson depth search followed by profile-based coefficient
construction. The non-local flux is a standard part of the KPP formulation and
is included by default.

## How KPP Is Used in Time Stepping

KPP is connected to all three Omega time steppers: Forward-Backward, RungeKutta2, and RungeKutta4. For
whichever stepper is active, KPP is computed **once per time step**, at the
start of the step (in contrast to MPAS-Ocean), before any tendency is evaluated.
Boundary-layer depth, viscosity, diffusivity, and the non-local flux profile
are then held fixed for the rest of the step and are used by:

- the non-local tracer tendency at every internal stage of the stepper, and
- the implicit vertical mixing solve applied after the time levels are
  advanced.

Because both use the same KPP fields, the non-local flux and the diffusivity
it is paired with are always consistent with one another.

The KPP diagnostics written to output for a step therefore describe the ocean
state at the **beginning** of that step, not the updated state at the end of
it. This is the same convention used by the MPAS-Ocean split-explicit
stepper.

## Configuration

KPP settings are under `VertMix: KPP` in `omega.yml`.

### Example

```yaml
VertMix:
  KPP:
    Enable: true
    CriticalBulkRichardsonNumber: 0.25
    SurfaceLayerExtent: 0.1
    MatchTechnique: SimpleShapes
    InterpType2: LMD94
    UseEnhancedDiffusion: true
    UseOSBLSmoothing: true
    UseLangmuirTurbulence: true
    IceFractionThresholdForLangmuir: 0.05
    IceFractionThresholdForMinimumOSBL: 0.15
    MinimumOSBLUnderSeaIce: 5.0
    DebugDiagnostics: false
```

### Key Options

| Key | Meaning | Typical default |
|---|---|---|
| `Enable` | Enable KPP mixing | `true` |
| `CriticalBulkRichardsonNumber` | OBL depth criterion threshold | `0.25` |
| `SurfaceLayerExtent` | Surface layer thickness as a fraction of the OBL depth ($\epsilon$ in Large et al. 1994) | `0.1` |
| `MatchTechnique` | How the K profile meets interior mixing at the OBL base: `SimpleShapes` or `MatchBoth` | `SimpleShapes` |
| `InterpType2` | Interpolation type used near OBL matching/base logic | `LMD94` |
| `UseEnhancedDiffusion` | Enable enhanced diffusion treatment near OBL base | `true` |
| `UseOSBLSmoothing` | Apply horizontal smoothing to the OSBL depth | `true` |
| `UseLangmuirTurbulence` | Apply Langmuir enhancement to the turbulent velocity scale | `true` |
| `IceFractionThresholdForLangmuir` | Above this ice fraction, disable Langmuir enhancement | `0.05` |
| `IceFractionThresholdForMinimumOSBL` | Above this ice fraction, enforce minimum OSBL depth | `0.15` |
| `MinimumOSBLUnderSeaIce` | Minimum OSBL depth under sea ice (m) | `5.0` |
| `DebugDiagnostics` | Enable additional KPP diagnostics/logging in debug workflows, and extend the Ri diagnostic profiles below the boundary layer base | `false` |

Independent of those options, KPP keeps the diagnosed OSBL depth at least as deep
as the bottom of the first active layer.

## Output and Diagnostics

To diagnose KPP, include KPP fields in output stream contents. Common fields:

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

The boundary layer search stops at the first level where the bulk Richardson
number reaches its critical value, so `BulkRichardsonNumber`,
`BulkRichardsonShear`, `UnresolvedShear` and `BuoyancyJump` are zero below the
boundary layer base. Set `DebugDiagnostics: true` to compute and output the
full water column profile of these four fields; it does not change
`OSBLDepth` or any mixing coefficient.

## Typical Workflow

1. Enable KPP and set baseline options in `omega.yml`.
2. Run a short case.
3. Inspect `OSBLDepth` and coefficient fields.
4. If needed, tune `CriticalBulkRichardsonNumber`, `MatchTechnique`, and
   `InterpType2`.
5. Re-run and compare diagnostics.

## Practical Notes

- Non-local tracer transport is required for physically correct KPP boundary
  layer tracer fluxes. It is applied through the
  `Tendencies: KPPNonLocalTracerFluxTendencyEnable` flag, which is on by default
  and should be left enabled; see
  [Tendency Terms](./TendencyTerms.md).
- Use `DebugDiagnostics` sparingly for troubleshooting targeted cases.
- When studying sea-ice regions, review minimum-OBL and ice-threshold options.
