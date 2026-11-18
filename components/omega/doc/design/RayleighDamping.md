(omega-design-rayleigh-damping)=
# Rayleigh Damping

**Table of Contents**
1. [Overview](#1-overview)
2. [Requirements](#2-requirements)
3. [Algorithmic Formulation](#3-algorithmic-formulation)
4. [Design](#4-design)
5. [Verification and Testing](#5-verification-and-testing)

## 1 Overview

An ocean initial condition assembled from observational climatologies is not in
dynamical balance with the model's discrete equations. When such a state is used
to start a simulation, the adjustment radiates large-amplitude barotropic and
inertia-gravity waves. These waves carry velocities far larger than the
climatological circulation, and they force the model to run with a very small
time step (or to fail outright on the CFL condition) until they have propagated
away or dissipated.

The standard remedy, used routinely to produce E3SM ocean initial conditions
with MPAS-Ocean, is a *dynamic adjustment* spin-up: the model is run for a short
period with a strong linear drag on momentum, the so-called Rayleigh damping
term, which removes the transient wave energy on a prescribed time scale. The
spin-up proceeds as a sequence of adjustment segments, each run as a separate
job. Across the sequence the damping coefficient is reduced and the time step is
increased, until the damping is small enough to be turned off and the time step
has reached its production value. At that point the state is dynamically
balanced and can be integrated forward.

The staging of the time step is as important as the staging of the coefficient.
Running the entire sequence at the small initial time step would take far longer
than necessary, and it would leave the stability of the production time step
untested until the adjustment was over. Each segment therefore serves double
duty: it damps transients, and it demonstrates that the state tolerates the
next, longer time step.

Omega currently has no such capability, which is a blocker for generating E3SM
initial conditions (see
[Omega issue #495](https://github.com/E3SM-Project/Omega/issues/495)). This
document describes the addition of Rayleigh damping to Omega's momentum
equation, following the formulation used by MPAS-Ocean but expressed in Omega's
non-Boussinesq, pseudo-thickness-weighted form and configured with Omega naming
conventions.

The scope here is deliberately narrow: a single, spatially uniform damping
coefficient applied to the normal velocity at all active layers, treated
implicitly in time. Two variants present in MPAS-Ocean are explicitly *not*
carried over:

- `config_Rayleigh_damping_depth_variable`, which divides the coefficient by the
  total column thickness. This is used by the subgrid wetting-and-drying work
  rather than by dynamic adjustment, and its coefficient has different units
  (m s$^{-1}$ rather than s$^{-1}$), which would make a single config option
  ambiguous.
- `config_Rayleigh_bottom_damping_coeff`, which applies the damping only in the
  bottom layer. That term is physically a linear bottom friction, not Rayleigh
  damping, and belongs with Omega's bottom drag rather than here. If it is
  needed later, it should be added as a new `Type` under the existing
  `BottomDragTendency` config group.

## 2 Requirements

### 2.1 Requirement: Linear damping of momentum with a configurable time scale

Omega shall provide a term that damps the normal velocity toward zero at a rate
set by a single, run-time configurable inverse time scale $c_R$ with units of
s$^{-1}$.

Dynamic adjustment requires control over how fast transient wave energy is
removed. A damping time scale of a few hours ($c_R \sim 10^{-4}$ s$^{-1}$) is
short enough to suppress barotropic and near-inertial transients within a few
days of simulation, while a subsequent segment with a weaker coefficient allows
the slower baroclinic adjustment to proceed without being over-damped.

### 2.2 Requirement: Damping is applied to all active layers of the momentum equation only

The damping shall act on the normal velocity in every active layer between the
top and bottom of each edge column, and shall not act on pseudo-thickness or on
tracers.

The transients being removed are momentum modes. Damping thickness or tracers
would relax the model toward a state that is not the intended initial condition
and would violate conservation of mass and tracer content. Applying the damping
over the whole column rather than only near the bottom removes barotropic and
baroclinic transients at comparable rates, which is the behavior needed for
adjustment.

### 2.3 Requirement: Unconditional stability for large damping coefficients

The discrete damping term shall be stable and monotone for arbitrarily large
values of $c_R \Delta t$, and shall not produce sign reversal or growth of the
velocity.

The whole point of the term is to allow the time steps to be as *large* as
possible during adjustment. An explicit treatment would impose
$c_R \Delta t < 1$ for monotone decay, which couples the very quantity being
tuned ($c_R$) to the quantity being maximized ($\Delta t$) and defeats the
purpose. An implicit treatment removes that coupling entirely.

### 2.4 Requirement: Disabled by default with no effect on existing results

Rayleigh damping shall be off by default, and when off, model results shall be
bit-for-bit identical to results without the feature.

Rayleigh damping is a spin-up tool, not part of the physics of a production
simulation. Leaving it enabled in a science run would artificially dissipate the
circulation. Bit-for-bit neutrality when disabled keeps existing test baselines
valid.

### 2.5 Requirement: Configuration errors are detected at initialization

If Rayleigh damping is enabled in a configuration where it cannot be applied,
Omega shall abort during initialization with an informative message rather than
silently ignoring the setting.

A silently dropped damping term would produce an adjustment run that appears to
succeed but yields an unbalanced initial condition, and the failure would only
be discovered much later in the E3SM initialization workflow.

### 2.6 Desired: Time-varying damping coefficient

It is desirable for the damping coefficient to decrease smoothly over the course
of an adjustment segment, for example linearly or exponentially between a
starting and an ending value.

The value of such a ramp is limited, and it is worth being precise about what it
would and would not buy. It would *not* replace the sequence of adjustment jobs:
as described in Section 1, each segment also raises the time step, and Omega has
no mechanism for changing the time step mid-run and no intention of adding one.
The sequence is therefore intrinsic to the workflow regardless of how $c_R$ is
specified. What a ramp would offer is the removal of the abrupt change in $c_R$
at each segment boundary, replacing a step with a continuous transition. That is
a nice-to-have rather than a motivating capability, which is why it is listed as
desired rather than required.

The design accommodates it at no cost: the coefficient is a plain scalar member
of the setup functor, so a ramp can later be implemented by updating that scalar
once per time step, with no change to the tridiagonal setup or to the config
group structure.

## 3 Algorithmic Formulation

### 3.1 Continuous form

Rayleigh damping adds a linear drag to the momentum equation. In terms of the
normal velocity at an edge, $u_e$, the added tendency is

$$
\left. \frac{\partial u_e}{\partial t} \right|_{\mathrm{Rayleigh}} = -c_R u_e ,
$$

where $c_R \ge 0$ has units of s$^{-1}$ and $1 / c_R$ is the e-folding time
scale of the damping. In isolation this gives pure exponential decay,
$u_e(t) = u_e(0) e^{-c_R t}$, at every wavenumber and every vertical mode. The
term is scale-selective in neither space nor time: it damps the mean circulation
along with the transients. That is acceptable for adjustment precisely because
the adjustment run is short and the transients are much larger in amplitude than
the circulation being preserved.

The term is not energy conserving by construction. It removes kinetic energy at
a rate

$$
\left. \frac{\partial}{\partial t} \frac{u_e^2}{2} \right|_{\mathrm{Rayleigh}}
= -c_R u_e^2 ,
$$

and that energy is not returned to any other reservoir in the model.

### 3.2 Layered, pseudo-thickness-weighted form

Omega's discrete momentum equation is non-Boussinesq and is advanced in
pseudo-thickness-weighted form, with $\left[\tilde{h}_k\right]_e$ the
pseudo-thickness of layer $k$ interpolated to edge $e$. Grouping the Rayleigh
term with the vertical momentum stress divergence, the portion of the momentum
equation that is treated implicitly is

$$
\left[\tilde{h}_k\right]_e \frac{\partial u_{e,k}}{\partial t}
= \left[ \tilde{\nu}_v \frac{\partial u}{\partial \tilde{z}} \right]_{e,k}
- \left[ \tilde{\nu}_v \frac{\partial u}{\partial \tilde{z}} \right]_{e,k+1}
- c_R \left[\tilde{h}_k\right]_e u_{e,k} .
$$

Because the damping rate is independent of thickness, the Rayleigh term picks up
the same $\left[\tilde{h}_k\right]_e$ factor as the time derivative. This is the
non-Boussinesq counterpart of the MPAS-Ocean form, in which the equivalent
system is divided through by the layer thickness and the term appears simply as
$c_R$ on the diagonal.

### 3.3 Discrete implicit form

Omega already solves the vertical momentum stress divergence implicitly with a
backward-Euler step, using the diffusion-form tridiagonal solver
`TriDiagDiffSolver` (see the [Tridiagonal Solver](TridiagonalSolver) design
document). The system solved at each edge for the updated normal velocity
$Y_k = u^{n+1}_{e,k}$ is

$$
-G_{k-1} Y_{k-1} + \left( G_{k-1} + G_k + H_k \right) Y_k - G_k Y_{k+1} = X_k ,
$$

where $G_k$ holds the off-diagonal viscous coupling between layers $k$ and
$k+1$, $H_k$ holds contributions that appear only on the diagonal, and $X_k$ is
the right-hand side. Without damping,

$$
H_k = \left[\tilde{h}_k\right]_e
+ \delta_{k,k_{\mathrm{bot}}} \, \Delta t \, C_D \, \frac{\rho}{\rho_0} \left| \mathbf{u}_e \right| ,
\qquad
X_k = \left[\tilde{h}_k\right]_e u^n_{e,k} ,
$$

with the second term in $H_k$ the optional implicit bottom drag, present only in
the bottom layer.

Backward-Euler treatment of the Rayleigh term contributes
$\Delta t \, c_R \left[\tilde{h}_k\right]_e u^{n+1}_{e,k}$ to the left-hand side,
so the *only* change to the system is an additional diagonal contribution in
every active layer:

$$
H_k = \left[\tilde{h}_k\right]_e \left( 1 + \Delta t \, c_R \right)
+ \delta_{k,k_{\mathrm{bot}}} \, \Delta t \, C_D \, \frac{\rho}{\rho_0} \left| \mathbf{u}_e \right| .
$$

Neither $G_k$ nor $X_k$ is modified, and the matrix keeps its symmetric,
diagonally dominant structure, so the existing solver applies unchanged.

In the limit of no vertical viscosity and no bottom drag, the system decouples
and reduces to

$$
u^{n+1}_{e,k} = \frac{u^n_{e,k}}{1 + \Delta t \, c_R} ,
$$

which is the analytic solution of the damping equation accurate to first order
in $\Delta t$, is strictly decreasing in magnitude, preserves sign, and is
bounded for any $\Delta t \, c_R \ge 0$. This satisfies the stability
requirement in Section 2.3, and is the relation the unit test in Section 5.1
checks. It is also identical to the single-active-layer branch of MPAS-Ocean's
`ocn_vel_vmix_tend_implicit_rayleigh`, after accounting for the fact that
MPAS-Ocean divides its system by layer thickness.

### 3.4 Typical coefficient values

For reference, the coefficients used in MPAS-Ocean dynamic adjustment sequences
are of order $10^{-4}$ s$^{-1}$ (a damping time scale of about 2.8 hours) for
the first, most aggressive segment, decreasing by roughly an order of magnitude
over several segments before damping is switched off. Omega's default value
follows this convention, but the feature is disabled by default so the default
value is inert.

## 4 Design

The design adds no new class and no new field. Omega already assembles the
implicit velocity solve in one place, the `VelVertMixSetupOnEdge` functor in
`VertMix.h`, and that functor already carries an optional diagonal contribution
for implicit bottom drag. Rayleigh damping is a second such contribution, so it
is added the same way: two new public members on the functor, a few lines in
`operator()`, config reading in `VertMix::init()`, and a consistency check in
`Tendencies`.

### 4.1 Data types and parameters

#### 4.1.1 Parameters

A new top-level `RayleighDamping` group is added to the Omega config, placed
after the `VertMix` group in `components/omega/configs/Default.yml`:

```yaml
  RayleighDamping:
    Enable: false
    DampingCoeff: 1.0e-4
```

- `Enable` (logical, default `false`). Turns Rayleigh damping on. When `false`,
  the damping term is not evaluated at all and results are bit-for-bit
  unchanged, satisfying requirement 2.4.
- `DampingCoeff` (real, default `1.0e-4`, units s$^{-1}$). The coefficient
  $c_R$. Any non-negative value is allowed; a negative value is a configuration
  error and aborts at initialization. The default corresponds to a damping time
  scale of about 2.8 hours, the typical value for the first segment of a dynamic
  adjustment, but has no effect unless `Enable` is `true`.

A top-level group is used rather than a subgroup of `Tendencies` because
Rayleigh damping is a self-contained spin-up capability rather than one of the
physical tendency terms, and because a top-level group is the natural place to
add the ramp parameters of requirement 2.6 later.

Note that `DampingCoeff` is read only at initialization. Changing it therefore
requires a new job segment. This is not a limitation in practice: the time step
is also fixed for the duration of a run, so an adjustment sequence must be
staged as separate jobs in any case.

Because the term is applied inside the implicit velocity vertical mixing solve,
enabling it requires the existing option

```yaml
  Tendencies:
    VelVertMixTendencyEnable: true
```

which is the default. This is the same constraint that already applies to
`BottomDragTendency` with `Mode: Implicit`.

#### 4.1.2 Class/structs/data types

No new class is introduced. Two public members are added to
`VelVertMixSetupOnEdge` in `components/omega/src/ocn/VertMix.h`, alongside the
existing bottom drag members:

```c++
class VelVertMixSetupOnEdge {
 public:
   bool Enabled;
   bool ImplicitBottomDragEnabled;
   Real BottomDragCoeff;
   bool RayleighDampingEnabled; ///< Enable Rayleigh damping flag
   Real RayleighDampingCoeff;   ///< Rayleigh damping coefficient (s^-1)
   Real LocRhoSw;
   ...
};
```

Both are initialized in the constructor initializer list in `VertMix.cpp` to
`false` and `0.0_Real` respectively, so that a `VertMix` instance created
without reading a config performs no damping.

### 4.2 Methods

No new public method is added. Four existing pieces of code change.

#### 4.2.1 `VelVertMixSetupOnEdge::operator()`

The functor already computes the diagonal contribution `H` for layer `K` at edge
`IEdge`. The damping contribution is added immediately after the
pseudo-thickness is assigned to `H` and before the bottom drag contribution:

```c++
      H = PseudoThickEdgeK;

      // Rayleigh damping: a linear drag on momentum applied in every active
      // layer. Treated implicitly, so it is unconditionally stable for any
      // DT * RayleighDampingCoeff.
      if (RayleighDampingEnabled) {
         H += DT * RayleighDampingCoeff * PseudoThickEdgeK;
      }
```

The functor is only invoked for layers between `KMin` and `KMax` at each edge;
layers outside that range are filled with the inactive-column values
(`G = 0`, `H = 1`, `X = 0`) by the caller in `VertMix::applyVelVertMixImplicit`.
The damping therefore automatically applies to exactly the active layers,
satisfying requirement 2.2, and needs no additional masking. The off-diagonal
term `G` and the right-hand side `X` are untouched.

Guarding the addition with `RayleighDampingEnabled` rather than relying on a
zero coefficient mirrors the existing bottom drag pattern, avoids the extra
work in the common case, and makes the disabled path trivially bit-for-bit.

#### 4.2.2 `VertMix::init()`

`VertMix::init()` already reads the `VertMix` config group and its `Background`,
`Convective`, and `Shear` subgroups. It gains a block that reads the new
top-level group, following the same error-handling and logging idiom:

```c++
   /// Get RayleighDamping group from Omega config
   Config RayleighConfig("RayleighDamping");
   Err += OmegaConfig->get(RayleighConfig);
   CHECK_ERROR_ABORT(
       Err, "VertMix::init: RayleighDamping group not found in Config");

   Err += RayleighConfig.get(
       "Enable", DefVertMix->VelVertMixSetup.RayleighDampingEnabled);
   CHECK_ERROR_ABORT(Err, "VertMix::init: Parameter RayleighDamping:Enable not "
                          "found in RayleighConfig");

   if (!DefVertMix->VelVertMixSetup.RayleighDampingEnabled) {
      LOG_INFO("VertMix::init: Rayleigh damping is disabled.");
   } else {
      Err += RayleighConfig.get(
          "DampingCoeff", DefVertMix->VelVertMixSetup.RayleighDampingCoeff);
      CHECK_ERROR_ABORT(Err, "VertMix::init: Parameter "
                             "RayleighDamping:DampingCoeff not found in "
                             "RayleighConfig");

      if (DefVertMix->VelVertMixSetup.RayleighDampingCoeff < 0.0_Real) {
         ABORT_ERROR("VertMix::init: RayleighDamping:DampingCoeff must be "
                     "non-negative but got {}",
                     DefVertMix->VelVertMixSetup.RayleighDampingCoeff);
      }

      LOG_INFO("VertMix::init: Rayleigh damping is enabled with coefficient "
               "{} s^-1.",
               DefVertMix->VelVertMixSetup.RayleighDampingCoeff);
   }
```

`VertMix::init()` is called from `OceanInit.cpp` before `Tendencies::init()`,
so the flag is set before the consistency check below runs.

#### 4.2.3 `Tendencies` configuration check

`Tendencies` reads `VelVertMixTendencyEnable` into
`VMix->VelVertMixSetup.Enabled` and already aborts if implicit bottom drag is
requested without it. An analogous check for Rayleigh damping is added
immediately after, satisfying requirement 2.5:

```c++
   if (this->VMix->VelVertMixSetup.RayleighDampingEnabled &&
       !this->VMix->VelVertMixSetup.Enabled) {
      ABORT_ERROR("Tendencies: RayleighDamping Enable requires "
                  "VelVertMixTendencyEnable to be true");
   }
```

#### 4.2.4 Documentation

The `\ConfigInput` block in the header comment of `VertMix.h` is extended with
the new group, and the user guide gains a short section describing the option
and its use in dynamic adjustment. A new entry `design/RayleighDamping` is added
to the design document list in `components/omega/doc/index.md`.

### 4.3 Items explicitly unaffected

- **Halo exchanges.** The damping is a pointwise operation on each edge column
  and introduces no new horizontal stencil, so no additional halo exchange is
  needed. The implicit solve is already performed over all edges including
  halos.
- **Restart and I/O.** No new prognostic or diagnostic field is created, so the
  restart contents and all I/O streams are unchanged. A restart written during
  an adjustment segment is read normally by a subsequent segment with a
  different coefficient.
- **Explicit tendency path.** `Tendencies::computeVelocityTendencies` and the
  functors in `TendencyTerms.h` are untouched.
- **Time steppers.** The implicit vertical mixing step is applied by the time
  steppers through `VertMix::VertMixImplicit`, so all time steppers that call it
  pick up the damping with no change.
- **Performance.** The added cost is one multiply-add per active edge layer
  inside a kernel that is already memory-bound, and only when the feature is
  enabled.

## 5 Verification and Testing

### 5.1 Test: single-layer analytic decay

In `components/omega/test/ocn/VertMixTest.cpp`, add a test that configures a
column with vertical viscosity and bottom drag set to zero, Rayleigh damping
enabled, and a uniform nonzero initial normal velocity. After a single call to
`VertMix::applyVelVertMixImplicit`, the normal velocity in every active layer at
every edge must equal $u^n / (1 + \Delta t \, c_R)$ to within round-off of the
working precision.

This is the exact analytic solution of the discrete system in the decoupled
limit, so the test verifies both the presence of the term and the correctness of
its pseudo-thickness weighting: an unweighted or doubly weighted term would give
a thickness-dependent answer and fail.

List which requirements it tests:
  - tests requirement 2.1
  - tests requirement 2.2

### 5.2 Test: disabled damping is bit-for-bit neutral

Run the existing `VertMixTest` cases with `Enable: false` and confirm the
velocity and tracer results are bit-for-bit identical to the values expected
before this feature was added. Separately, confirm that with `Enable: true` and
`DampingCoeff: 0.0` the results are also bit-for-bit identical, which verifies
that the added diagonal contribution vanishes exactly rather than approximately.

List which requirements it tests:
  - tests requirement 2.4

### 5.3 Test: stability at large coefficient

Repeat the configuration of Section 5.1 with $c_R \Delta t = 10^3$, far outside
the range an explicit scheme could tolerate. The resulting velocity must be
reduced in magnitude by the expected factor, must retain the sign of the initial
velocity, and must be finite. The test is repeated with both signs of the
initial velocity.

List which requirements it tests:
  - tests requirement 2.3

### 5.4 Test: combined with vertical viscosity and bottom drag

Configure a multi-layer column with nonuniform pseudo-thickness, nonzero
vertical viscosity, implicit bottom drag enabled, and Rayleigh damping enabled.
Assemble the same tridiagonal system independently on the host with a
straightforward Thomas-algorithm implementation in the test, and require the
solver result to match to round-off.

This verifies that the damping contribution is added to the diagonal only, that
it does not perturb the off-diagonal coupling or the right-hand side, and that
it composes correctly with the bottom drag term that shares the same diagonal.

List which requirements it tests:
  - tests requirement 2.1
  - tests requirement 2.2

### 5.5 Test: configuration error handling

Confirm that a configuration with `RayleighDamping: Enable: true` and
`Tendencies: VelVertMixTendencyEnable: false` aborts during initialization with
the expected message, and that a negative `DampingCoeff` aborts during
`VertMix::init()`.

List which requirements it tests:
  - tests requirement 2.5

### 5.6 Test: dynamic adjustment in a realistic configuration

As a system-level check, run a global Omega configuration started from a
climatological initial condition, first without damping and then with
`DampingCoeff: 1.0e-4`. The damped run must show the maximum normal velocity and
the globally averaged kinetic energy decaying over the first few simulated days,
and must remain stable at a time step at which the undamped run either fails or
requires a substantially smaller step. This is the end-to-end demonstration that
the feature serves its purpose, and is the natural place for a Polaris dynamic
adjustment task.

List which requirements it tests:
  - tests requirement 2.1
  - tests requirement 2.3
