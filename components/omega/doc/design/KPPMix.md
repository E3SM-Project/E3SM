(omega-design-kppmix)=
# KPP Boundary Layer Mixing

**Table of Contents**
1. [Overview](#1-overview)
2. [Requirements](#2-requirements)
3. [Algorithmic Formulation](#3-algorithmic-formulation)
4. [Design](#4-design)
5. [Verification and Testing](#5-verification-and-testing)

## 1 Overview

This document describes the Omega implementation of the K Profile Parameterization
(KPP) ocean boundary layer mixing. KPP computes a boundary-layer depth, vertical
viscosity, vertical diffusivity, and a non-local tracer flux shape implemented outside
the implicit vertical mixing routine.  The implementation follows that in MPAS-Ocean and
uses direct ports of the functions defined in the [CVMix](https://github.com/CVMix/CVMix-src)
version of KPP.

The implementation is in `KPPMix` and is integrated with the Omega tendency and
RK2, RK4, and Forward-backward stepping routines. Relative to broad vertical mixing documentation, this
page focuses specifically on KPP theory, algorithmic choices, and testing.

## 2 Requirements

### 2.1 Requirement: Boundary-layer depth from bulk Richardson criterion

Following [Large et al (1994)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/94RG01872),
the OBL depth must be diagnosed from a bulk Richardson criterion so that
mixing depth responds to evolving stratification, shear, and surface forcing. It also
must include a unresolved turbulent shear contribution.

### 2.2 Requirement: Coefficients must be computable in parallel over columns

The KPP implementation must operate over many columns in parallel using Omega
array/kernels, rather than serial single-column calls.

### 2.3 Requirement: Compatible with additive vertical-mixing framework

KPP viscosity/diffusivity fields must be compatible with existing Omega vertical
mixing infrastructure so that other chosen vertical mixing sources can be merged
with KPP.

### 2.4 Desired: Non-local flux and profile matching controls

KPP will support a non-local tracer flux from LMD94 and include configurable
viscosity/diffusivity matching at the base of the boundary layer.

## 3 Algorithmic Formulation

The implementation follows a two-stage KPP structure.

### 3.1 Stage 1: OBL depth search

For each water column, OBL depth $h$ is diagnosed by searching downward until
bulk Richardson number reaches a critical value:

$$
Ri_b(z) = \frac{\Delta b(z)\, z}{|\Delta \mathbf{U}(z)|^2 + V_t^2(z)}
$$

with threshold

$$
Ri_b(h) = Ri_{crit}.
$$

Here, $\Delta b$ is buoyancy jump relative to a surface layer average,
$|\Delta \mathbf{U}|^2$ is shear contribution again computed relative to the
surface layer average, and $V_t^2$ is unresolved shear.

When the bulk Richardson number falls between model layers, quadratic interpolation
is utilized to find the depth.  In addition the boundary layer depth is constrained
to fall between a configurable minimum OBL under sea ice and a maximum set by the water
column depth.

The boundary layer depth search is a growing inner loop.  The outer loop iterates over all model layers.
The current model layer is set as a boundary layer depth candidate and $Ri_b$ is calculated for all
model layers shallower than the current depth.  If any layer in the inner loop has an $Ri_b*StopOBL$ that
exceeds $Ri_{crit}$ the loop terminates.

### 3.2 Stage 2: KPP coefficients and non-local flux

Given a diagnosed $h$, KPP computes viscosity and diffusivity coefficients at the top of every
Omega cell except the surface and the bottom using shape functions
in normalized depth $\sigma = -d/h$, where $d$ is the depth relative to the sea surface height, not the physical depth:

$$
K_m(\sigma) = h\, w_m(\sigma)\, M_1(\sigma),
$$

$$
K_s(\sigma) = h\, w_s(\sigma)\, S_1(\sigma),
$$

where $w_m$ and $w_s$ are turbulent velocity scales from Monin-Obukhov style
stability functions, see Appendix B of [Large et al, 1994](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/94RG01872). $M_1$ and $S_1$ are shape functions.
The generic form of the shape function is given by

$$
X(\sigma) = c_1 \sigma^3 + c_2 \sigma^2 + c_3 \sigma + c_4
$$

The coefficients are determined by various conditions, e.g., zero viscosity and diffusivity at the
surface, assuming a linear reduction of the turbulent flux with distance from the surface in the
surface layer.  As in MPAS-Ocean, we include two options to determine the final coefficients.  The original
version of KPP matches predicted viscosities and diffusivities to those predicted by other schemes
(e.g., shear instability driven mixing) and a second option where viscosities and diffusivities are
instead additive.  In the latter case, the shape function greatly simplifies to $X(\sigma) = \sigma(1-\sigma)^2.
When diffusivities and viscosities are matched at the OSBL base, only KPP is active in the OSBL.  When
diffusivities and viscosities are not matched, we assume that convective diffusivities are not applied
in the OSBL.  KPP and convective mixing target the same physical processes in the OSBL, so including
both is double counting.

For either shape function, enhanced diffusivity can be included near the boundary layer base.  This can smooth boundary layer deepening in time.

The non-local tracer flux uses the scalar shape function, scaled by the constant
$C_s$ from Eq. (20) of Large et al. (1994) rather than by $h\, w_s$:

$$
\gamma_s(\sigma) = C_s\, S_1(\sigma).
$$

Here $S_1$ is always the unmatched scalar shape $\sigma(1-\sigma)^2$, regardless
of the `MatchTechnique` setting. Matching is a property of the diffusivity
profile only: the matched shape is non-zero at $\sigma = -1$ by construction, so
reusing it for $\gamma_s$ would leave a finite non-local flux at the boundary
layer base that drops discontinuously to zero immediately below it. CVMix draws
the same distinction, exposing the non-local shape as a separate choice from the
matching option.



## 4 Design

### 4.1 Data types and parameters

#### 4.1.1 Parameters

KPP is configured from the `VertMix: KPP` YAML group. Key parameters include:

- `Enable`
- `CriticalBulkRichardsonNumber`
- `MatchTechnique` (`SimpleShapes` or `MatchBoth`)
- `InterpType2` (`LMD94`, `Linear`, `Quadratic`, `Cubic`)
- `UseEnhancedDiffusion`
- `IceFractionThresholdForLangmuir`
- `IceFractionThresholdForMinimumOSBL`
- `MinimumOSBLUnderSeaIce`
- `DebugDiagnostics`

Defaults and usage examples are documented in the user guide page:
[KPP in the User Guide](../userGuide/KPPMix.md).

#### 4.1.2 Class/data structure

`KPPMix` is a singleton that owns persistent output fields and KPP-only
workspaces, including:

- `OSBLDepth`, `OSBLDepthIndex`
- `VertDiff`, `VertVisc`
- `VertNonLocalFlux`
- diagnostics such as `BulkRichardsonNumber`, `BulkRichardsonShear`,
  `UnresolvedShear`, `BuoyancyJump`, and `TurbulentVelocityScale`
- surface-referenced pressure, reconstructed edge tangential velocity, and
   temporary ice fraction

### 4.2 Methods

The production interface is `update(...)`, which receives canonical tracer,
velocity, EOS, and forcing state, prepares KPP-derived state, validates
required non-local tracer fluxes, and computes the KPP coefficients. It uses
zero-copy temperature and salinity tracer subviews; the forcing state retains
ownership of surface tracer fluxes.

The lower-level numerical stage is:

```c++
void computeKPPMix(const Array2DReal &PotentialDensity,
                   const Array2DReal &NormalVelocity,
                   const Array2DReal &TangentialVelocity,
                   const Array1DReal &SurfaceFrictionVelocity,
                   const Array1DReal &SurfaceBuoyancyFlux,
                   const Array2DReal &BruntVaisalaFreqSq,
                   const Array1DReal &IceFraction,
                   const Array1DReal &WindSpeed10m = Array1DReal());
```

Internal stages:
- `computeOSBLDepth(...)`
- `computeMixingCoefficients(...)`, using surface forcing and optional
   interior mixing coefficients

### 4.3 Time stepper coupling behavior

KPP is coupled to all three Omega time steppers -- Forward-Backward,
RungeKutta2, and RungeKutta4. For every stepper, KPP is evaluated exactly
once per time step, at the start of the step on the state at time $n$, before
any tendency is evaluated. The resulting boundary-layer depth, viscosity,
diffusivity, and non-local flux profile are then held fixed for the remainder
of the step.

This design differs from MPAS-Ocean, where the boundary layer depth, diffusivity, and viscosity are computed at the end of the time step and the non local flux is applied on the following timestep.  The primary advantage of this new approach is:

1. **Consistency.** The non-local flux $\gamma_s$ applied in the tracer
   tendency at each stage and the diffusivity $K_s$ used by the end-of-step
   implicit vertical mixing solve are derived from the same OBL depth and the
   same shape function $S_1(\sigma)$. Recomputing KPP at each stage would
   pair a stage-dependent $\gamma_s$ with a different $K_s$, breaking the
   correspondence described in section 3.2.

Because KPP is evaluated before the step advances, the KPP diagnostics
written for a step describe the state at the beginning of that step. The
coefficients are lagged relative to the state during the implicit solve at the end of the
step. Full call-flow detail per stepper is described for developers and
users in:

- [Developer KPP workflow](../devGuide/KPPMix.md)
- [User runtime notes](../userGuide/KPPMix.md)

## 5 Verification and Testing

### 5.1 Unit-level checks

The KPP unit-test driver is [KPPMixTest.cpp](../../test/ocn/KPPMixTest.cpp).
It supports focused test groups as well as an `all` group that runs all 28
unit tests:

- `profiles` (10 tests) checks the stability-function (used for the turbulent velocity scales) branches and continuity,
   simple shape and matched-shape functions, Langmuir utilities, OSBL depth clamping
   and indexing, turbulent velocity scales, and non-local flux scaling.
- `osbl` (7 tests) checks bulk-Richardson threshold crossing and quadratic interpolation,
   nonuniform and horizontally varying layer pseudo thickness, SSH-offset invariance to confirm positive definite geometric depths,
   partial columns, and area-weighted OSBL
   smoothing.
- `vmix` (8 tests) checks wind-only and convection-only profiles, stable and
   zero-forcing behavior, unmatched non-local profiles, `MatchBoth` interior
   coefficient matching, enhanced diffusion at the OSBL base, valid vertical
   domain edges, and invalid wet bounds.
- `integration` (3 tests) checks configuration parameters, the enabled
   end-to-end KPP call and diagnostics, and preservation of existing output
   fields when KPP is disabled.

The driver also provides three separate expected-failure cases:
`config-gradient`, `config-unsupported`, and `config-parabolic`. These inject
unsupported `MatchTechnique` values and verify that `KPPMix::init()` rejects
them; they are not included in the `all` group. The test registration and
group dispatch are shown at the end of [KPPMixTest.cpp](../../test/ocn/KPPMixTest.cpp).

Together, these tests verify:

- OBL depth search threshold crossing, interpolation, diagnostics, and edge
   handling
- KPP profile, matching, enhanced-diffusion, and non-local-flux invariants
- Langmuir, sea-ice, smoothing, stable-forcing, and zero-forcing behavior
- KPP's separation from shared `VertMix` background mixing and its handling of
   supplied interior coefficients
- Correct enabled, disabled, and invalid-configuration behavior

Tests cover requirements: 2.1, 2.2, 2.3, 2.4.

### 5.2 Polaris testing

The single column test case in Polaris has been extended to include a new set of regimes
to exercise KPP across a set of stratifications and forcings.  These tests are primarily
drawn from [Wagner et al., 2025](https://doi.org/10.1029/2024MS004522) and
[Van Roekel et al., 2018](https://doi.org/10.1029/2018MS001336), hereafter W25 and VR18 respectively.

The new tasks include

| Test | Description |
| :-------: | :---------- |
| Wind-only, no Coriolis. | Tests mechanical shear-driven OSBL growth.<br>Can be compared to solution from [Kato and Phillips](https://doi.org/10.1017/S0022112069000784) |
| Free convection (cooling) | VR18 FC case including a comparison to the analytic solution in Eq. F11 |
| Mixed forcing convection | Wind, Evaporation, and cooling (all destabilizing) |
| Nonlocal supression | Stable heating with evaporation (VR18 CEW case).<br>Tests suppression of nonlocal forcing under stable buoyancy forcing |
| Langmuir | Uses a fixed Stokes drift to test langmuir turbulence enabled and disabled |
| Sea ice | Wind forcing with a sea-ice fraction.  Tests OSBL clamping under sea ice. |
| Convection with Evaporation | VR18 FCE case, cooling with evaporation.<br>Tests the salinity contribution to surface buoyancy forcing. |
| Cooling with ML | VR18 FCML case with surface cooling on an existing ML. |
| Strong convection and stratification | W25 case showing KPP creating heating in the presence of cooling<br>This is expected behavior |

These tests include conservation and physics checks for each regime.  These tests also include plots of critical quantities like the bulk Richardson number,
ocean surface boundary layer depth, vertical diffusivity and viscosity, and the nonlocal flux.

Global stand-alone Omega cases, forced by annual averaged ERA-5 net surface heat, freshwater, and
momentum fluxes provides a qualitative assessment of KPP behavior.
