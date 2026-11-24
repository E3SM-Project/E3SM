(split-time-stepping)=
# Split Time Stepping

<!--- use table of contents if desired for longer documents  -->
**Table of Contents**

1. [Overview](#1-overview)
2. [Requirements](#2-requirements)
3. [Algorithmic Formulation](#3-algorithmic-formulation)
4. [Design](#4-design)
5. [Verification and Testing](#5-verification-and-testing)


## 1. Overview

To enhance computational efficiency by allowing longer time steps, ocean models require split barotropic-baroclinic time-stepping methods. The implementation described here is based on the approach of [Higdon (2005)](https://www.sciencedirect.com/science/article/pii/S0021999104005236) and the MPAS-Ocean `split_explicit` scheme, with modifications for Omega's non-Boussinesq pseudo-height $\tilde{z}$ vertical coordinate. The pseudo-height coordinate and pseudo thickness $\tilde h$ are defined by Eqs. {eq}`pseudo-height` and {eq}`def-pseudo-thickness`, respectively, in the {ref}`Omega V1 governing equations <omega-design-governing-eqns-omega1>` design document.

Omega currently implements two closely related steppers: `SplitExplicitRK2`, which explicitly subcycles the barotropic mode, and `UnsplitRK2`, which uses the same outer predictor-corrector framework without a barotropic split. The only implemented barotropic algorithm is `Predictor-Corrector`; implicit barotropic stepping and an Adams--Bashforth split-explicit method are future extensions.

The split-explicit method involves the following sequence:

- Decompose the velocity into barotropic and baroclinic components.
- Advance the baroclinic velocity with the long baroclinic time step and form the depth-averaged forcing $\overline{G}$ Eq. {eq}`split-stage1-gbar` that drives the barotropic momentum equation, Eq. {eq}`split-barotropic-momentum`.
- Subcycle the barotropic velocity with short explicit barotropic time steps.
- Construct a corrected transport velocity, update pseudo thickness and tracers, and reconstruct the full velocity.

## 2. Requirements

### 2.1 Requirement: A split time-stepping method in the pseudo-height coordinate

The algorithm is based on Section 2.3 of [Higdon (2005)](https://www.sciencedirect.com/science/article/pii/S0021999104005236), with modifications to accommodate the $\tilde{z}$-coordinate variables in the non-Boussinesq framework. It accepts the outer time step, requested barotropic time step, barotropic algorithm, number of outer predictor-corrector iterations, number of baroclinic Coriolis iterations, and an option to recompute the velocity split. The `UnsplitRK2` variant advances the full velocity ${\bf u}_{e,k}$ instead of the baroclinic velocity ${\bf u}'_{e,k}$, so the barotropic subcycle is redundant and is bypassed. Its stable time step is then set by the surface-gravity-wave CFL condition rather than by the internal-wave CFL condition, so `UnsplitRK2` is intended for verification and for idealized configurations rather than for production simulations.

### 2.2 Requirement: Stable time integration for long-term high-resolution simulations

Stability constrains the maximum allowable timestep, which in turn affects computational cost. The implemented split-explicit time stepping methods must allow for reasonably long timesteps while preventing numerical instabilities, such as those arising from internal gravity waves and barotropic modes, which are particularly important in global-scale and high-resolution ocean modeling. At a minimum, the time-stepping approach used in Omega V2 should accommodate the same timestep sizes as MPAS-Ocean for both the baroclinic and barotropic subsystems since Omega V2 is non-Boussinesq but hydrostatic.

### 2.3 Requirement: Modularization of baroclinic and barotropic time-stepping methods

Modularity ensures ease of testing and future-proofing of the Omega V2 codebase. Implementing a modular design enables mix-and-match time-stepping schemes of the baroclinic and barotropic subsystems, straightforward integration of alternative time-stepping schemes, and easier maintenance by separating the baroclinic and barotropic time-stepping codes, thereby enhancing flexibility.

## 3. Algorithmic Formulation

### 3.1 Barotropic (external) and baroclinic (internal) mode splitting

The split-explicit method separates ocean velocity into depth-integrated barotropic and depth-dependent baroclinic components. This separation allows the computationally expensive slow baroclinic modes to be advanced with a long time step, while the computationally inexpensive fast barotropic mode, which carries the surface gravity waves that impose the most restrictive CFL limit, is subcycled with short time steps. The benefit of the split is computational efficiency rather than accuracy. Mode splitting introduces a splitting error, since the two subsystems are coupled only once per baroclinic step. That error is controlled by the outer predictor-corrector iteration, by averaging the barotropic transport over the subcycles, and by the transport-velocity correction of Section 3.2.4, which restores consistency between the time-averaged barotropic transport and the layered thickness flux.

The layered discrete governing equations for Omega V2 are described in the {ref}`Omega V1 governing equations <omega-design-governing-eqns-omega1>` design document. The mass, tracer, and velocity equations used as the starting point for the mode split are summarized below.

**Mass:**

$$
\frac{\partial \tilde{h}_{i,k}}{\partial t}
+ \nabla \cdot \left( [\tilde{h}_{k}]_e {\bf u}_{e,k} \right)
+ \left[ \tilde{W}_{tr} \right]^{\text{top}}_{k}
- \left[ \tilde{W}_{tr} \right]^{\text{top}}_{k+1}
= 0 .
$$ (split-discrete-mass)

**Tracer:**

$$
\frac{\partial \tilde{h}_{i,k}\varphi_{i,k}}{\partial t}
+ \nabla \cdot \left( [\tilde{h}_{i,k}\varphi_{i,k}]_e {\bf u}_{e,k} \right)
+ \left\{
\left[\varphi \tilde{W}_{tr}\right]^{\text{top}}_{k}
-
\left[\varphi \tilde{W}_{tr}\right]^{\text{top}}_{k+1}
\right\}
= [D_h^\varphi]_{i,k} - [D_v^\varphi]_{i,k},
$$ (split-discrete-tracer)

where

$$
[D_h^\varphi]_{i,k}
= \kappa_{2,e}\nabla^2\varphi_{i,k}
- \kappa_{4,e}\nabla^4\varphi_{i,k},
$$ (split-tracer-horizontal-diffusion)

and

$$
[D_v^\varphi]_{i,k}
= [\tilde{\kappa}_{v}]_{i,k}
\left[\frac{\partial h\varphi}{\partial \tilde{z}}\right]_{i,k}
-
[\tilde{\kappa}_{v}]_{i,k+1}
\left[\frac{\partial h\varphi}{\partial \tilde{z}}\right]_{i,k+1} .
$$ (split-tracer-vertical-diffusion)

**Velocity:**

$$
\frac{\partial {\bf u}_{e,k}}{\partial t}
+ \left[ {\bf k} \cdot \nabla \times {\bf u}_{e,k} + f_v \right]_e {\bf u}^{\perp}_{e,k}
+ [\nabla K]_e
+ \frac{1}{[\tilde{h}_{i,k}]_e}
\left[\tilde{W}_{tr}\frac{\partial U}{\partial \tilde{z}}\right]_{e,k}
= - [(\alpha \nabla p + \nabla \Phi)]_{e,k}
+ [D_h^{\bf u}]_{e,k}
- [D_v^{\bf u}]_{e,k},
$$ (split-discrete-velocity)

where

$$
[D_h^{\bf u}]_{e,k}
= \nu_{2,e}\nabla^2 {\bf u}_{e,k}
- \nu_{4,e}\nabla^4 {\bf u}_{e,k},
$$ (split-velocity-horizontal-diffusion)

and

$$
[D_v^{\bf u}]_{e,k}
= \frac{1}{[\tilde{h}_{i,k}]_e}
\left\{
[\tilde{\nu}_{v}]^{\text{top}}_{e,k}
\left[\frac{\partial {\bf u}}{\partial \tilde{z}}\right]_{e,k}
-
[\tilde{\nu}_{v}]^{\text{top}}_{e,k+1}
\left[\frac{\partial {\bf u}}{\partial \tilde{z}}\right]_{e,k+1}
\right\} .
$$ (split-velocity-vertical-diffusion)

Define the barotropic velocity, baroclinic velocity, barotropic pressure, and barotropic pressure anomaly as follows.

**Barotropic velocity:**

$$
\overline{{\bf u}}
\equiv
\frac{1}{\tilde{H}}
\sum_{k=0}^{K_{\max}} \tilde{h}_{k}{\bf u}_{k}.
$$ (split-barotropic-velocity)

**Baroclinic velocity:**

$$
{\bf u}'_k \equiv {\bf u}_k - \overline{{\bf u}}.
$$ (split-baroclinic-velocity)

**Barotropic pressure:**

$$
B \equiv p^{\text{floor}} - p^{\text{surf}}.
$$ (split-barotropic-pressure)

**Barotropic pressure anomaly:**

$$
B' \equiv B - \rho_0 g b,
$$ (split-barotropic-pressure-anomaly)

so that

$$
B = B' + \rho_0 g b,
$$ (split-barotropic-pressure-reconstruction)

and

$$
p^{\text{floor}}
= p^{\text{surf}} + B' + \rho_0 g b,
$$ (split-floor-pressure-reconstruction)

where $b \equiv -z^{\text{floor}}$ is the geometric bottom depth.

Define the column-integrated pseudo thickness and geometric thickness as follows.

**Column-integrated pseudo thickness**, or mass-equivalent column depth:

$$
\tilde{H}
\equiv
\sum_{k=0}^{K_{\max}} \tilde{h}_{k}
=
\frac{1}{\rho_0 g}\sum_{k=0}^{K_{\max}}
\left(p_{k+1}^{\text{top}} - p_k^{\text{top}}\right)
=
\frac{1}{\rho_0 g}\left(p^{\text{floor}}-p^{\text{surf}}\right)
=
\frac{B}{\rho_0 g}
=
\frac{1}{\rho_0 g}\left(B' + \rho_0 g b\right).
$$ (split-column-pseudo-thickness)

**Column-integrated geometric thickness**, or geometric column depth:

$$
H
=
\sum_{k=0}^{K_{\max}} h_k
=
\sum_{k=0}^{K_{\max}} \frac{\rho_0}{\rho_k}\tilde{h}_k
=
\rho_0 \sum_{k=0}^{K_{\max}} \alpha_k \tilde{h}_k
=
\rho_0 S,
$$ (split-column-geometric-thickness)

where

$$
S = \sum_{k=0}^{K_{\max}} \alpha_k \tilde{h}_k
$$ (split-column-specific-volume)

is the column-integrated specific volume.

Sea-surface height can be diagnosed by

$$
\eta = \rho_0 S - b.
$$ (split-diagnostic-ssh)

The sea-surface height is decomposed into two components:

$$
\eta = \eta^{\text{mass}} + \eta^{\text{steric}},
$$ (split-ssh-decomposition)

where

$$
\eta^{\text{mass}} = \frac{B'}{\rho_0 g},
\qquad
\eta^{\text{steric}} = \rho_0 S - \frac{B}{\rho_0 g}.
$$ (split-mass-and-steric-ssh)

This expression is equivalent to Eq. (10.60) in [Griffies (2012)](https://mom-ocean.github.io/assets/pdfs/MOM5_manual.pdf) and [Madec et al. (2015)](https://epic.awi.de/id/eprint/39698/1/NEMO_book_v6039.pdf), but formulated using $S$ and $B$. The term $\eta^{\text{steric}}$ represents the change in sea-surface height resulting from the non-Boussinesq steric effect, while $\eta^{\text{mass}}$ is the mass-related component of sea-surface height that is advanced in time in the barotropic system as the fast process.

Under the Boussinesq approximation, the steric height $\eta^{\text{steric}}$ reduces to zero:

$$
\eta^{\text{steric}} = H - \tilde{H};
\qquad
\rho \rightarrow \rho_0;
\qquad
H - H = 0.
$$ (split-boussinesq-steric-limit)

A common expression for steric height, following [Griffies (2012)](https://mom-ocean.github.io/assets/pdfs/MOM5_manual.pdf) and [Madec et al. (2015)](https://epic.awi.de/id/eprint/39698/1/NEMO_book_v6039.pdf), is

$$
\eta^{\text{steric}}
= -\int_{-H}^{\eta}
\left(\frac{\rho - \rho_0}{\rho_0}\right) dz .
$$ (split-steric-height-griffies-madec)

The derivation below demonstrates that this definition of $\eta^{\text{steric}}$ is identical to the definition employed in the split system of Omega:

$$
\begin{aligned}
\eta^{\text{steric}}
&= -\int_{-H}^{\eta}
\left(\frac{\rho - \rho_0}{\rho_0}\right) dz \\
&= \int_{-H}^{\eta}
\left(1 - \frac{\rho}{\rho_0}\right) dz \\
&= \int_{-H}^{\eta} dz
- \int_{-H}^{\eta}\frac{\rho}{\rho_0} dz \\
&= \int_{-H}^{\eta} dz
- \int_{-H}^{\eta} d\tilde{z} \\
&= H - \tilde{H} \\
&= \rho_0 S - \frac{B}{\rho_0 g} .
\end{aligned}
$$ (split-steric-height-derivation)

The fourth line uses the definition of the pseudo-height coordinate, $\rho\,dz=\rho_0\,d\tilde{z}$, given in the {ref}`Omega V1 governing equations <omega-design-governing-eqns-omega1>` design document, and the fifth line uses Eqs. {eq}`split-column-pseudo-thickness` and {eq}`split-column-geometric-thickness`.

The barotropic system is obtained from the layered equations of Section 3.1 by summing the mass equation over the column and by taking the pseudo-thickness-weighted vertical average of the momentum equation. Summing Eq. {eq}`split-discrete-mass` over $k$ cancels the vertical transport terms between adjacent interfaces, leaving only the surface value, which is the freshwater flux $Q$:

$$
\frac{\partial \tilde{H}_i}{\partial t}
+ \left[\nabla \cdot \sum_{k=0}^{K_{\max}}
  \left([\tilde{h}_{k}]_e {\bf u}_{e,k}\right)\right]_i
= -Q_i .
$$ (split-column-integrated-mass)

To advance the column mass in the non-Boussinesq system, the barotropic prognostic variable is taken to be the pressure anomaly $B'$ rather than $\tilde{H}$ itself. Replacing the column transport by $[\tilde{H}_i]_e\overline{{\bf u}}_e$ using Eq. {eq}`split-barotropic-velocity`, and $\tilde{H}$ by $(B'+\rho_0 g b)/(\rho_0 g)$ using Eq. {eq}`split-column-pseudo-thickness`, converts this into an equation for the barotropic pressure anomaly. The time-independent term $\rho_0 g b$ drops out of the time derivative, and multiplying through by $\rho_0 g$ gives the barotropic continuity equation below.

Averaging Eq. {eq}`split-discrete-velocity` over the column in the same pseudo-thickness-weighted sense separates the pressure-gradient term into a depth-mean part, $-[\overline{\alpha}_i]_e[\nabla B_i']_e$, which is the fast surface-pressure-gradient term retained explicitly in the barotropic system, and a remainder. The planetary Coriolis term is also retained explicitly, because rotation modifies the surface gravity waves on the barotropic time scale. All remaining depth-averaged terms are collected into a single forcing $\overline{G}_e$ that is diagnosed once per outer iteration and held fixed through the barotropic subcycle. Subtracting the same $\overline{G}_e$ from the layered momentum equation gives the baroclinic momentum equation, so that the sum of the two subsystems recovers Eq. {eq}`split-discrete-velocity`.

The barotropic continuity and momentum equations are written as follows.

**Barotropic continuity equation:**

$$
\frac{\partial B_i'}{\partial t}
+ \left[
\nabla \cdot \left( [(B_i' + \rho_0 g b_i)]_e \overline{{\bf u}}_e \right)
\right]_i
= -\rho_0 g Q_i .
$$ (split-barotropic-continuity)

**Barotropic momentum equation:**

$$
\frac{\partial \overline{{\bf u}}_e}{\partial t}
+ f_e \overline{{\bf u}}^{\perp}_e
= -[\overline{\alpha}_i]_e [\nabla B_i']_e
+ \overline{G}_e .
$$ (split-barotropic-momentum)

**Baroclinic momentum equation:**

$$
\frac{\partial {\bf u}'_{e,k}}{\partial t}
= -[f_v]_e {\bf u}_{e,k}^{\prime\perp}
+ \Gamma_{e,k}
+ [\overline{\alpha}_i]_e [\nabla B_i']_e
- \overline{G}_e,
$$ (split-baroclinic-momentum)

where

$$
\Gamma_{e,k}
\equiv
-[\nabla K]_e
-(\alpha \nabla p + \nabla \Phi)_{e,k}
-
[ {\bf k} \cdot \nabla \times {\bf u}_{e,k}]_e {\bf u}^{\perp}_{e,k}
-
\frac{1}{[\tilde{h}_{i,k}]_e}
\left[\tilde{W}_{tr}\frac{\partial U}{\partial \tilde{z}}\right]_{e,k}
+ [D_h^{\bf u}]_{e,k}
+ [F_s^{\bf u}]_{e,k}
+ [F_b^{\bf u}]_{e,k}.
$$ (split-gamma-definition)

Here, $\overline{G}_e$ includes all remaining terms in the barotropic equation.

Each term of $\Gamma$ is applied only when its own tendency switch is enabled,
and each is evaluated from the *full* normal velocity of the working state
rather than from the baroclinic velocity alone. $[F_s^{\bf u}]$ and
$[F_b^{\bf u}]$ are the surface stress forcing and the explicit bottom drag,
which the implementation includes in the baroclinic tendency alongside the
terms of Eq. {eq}`split-discrete-velocity`. A user-registered custom velocity
tendency, when present, is also added here. $Q$ is the surface freshwater flux.

The vertical mixing term $[D_v^{\bf u}]$ does not appear in $\Gamma$. Vertical
mixing of momentum and tracers is applied implicitly once per slow step, after
the outer iterations and the time-level rotation have completed; see
Section 3.2.6. Consequently $\Gamma$ contains only the terms integrated
explicitly by the outer scheme.

### 3.2 Split-explicit time stepping algorithm

The mode-splitting algorithm follows the MPAS-Ocean split-explicit framework. It advances baroclinic velocity over the slow time step, explicitly subcycles the barotropic system, and uses the time-averaged barotropic transport to correct the velocity used for pseudo-thickness and tracer transport. The implementation calls this option `SplitExplicitRK2`.

`SplitExplicitRK2` uses two outer iterations by default. The first evaluates momentum at the time-$n$ state and constructs a midpoint state. The second evaluates momentum at that midpoint and retains a full-step conservative pseudo-thickness and tracer update. Each outer iteration passes the full slow time step $\Delta t$ to all three stages; it does not use $\Delta t/\text{NTimeStepIteration}$. The midpoint is produced by the one-half coefficients in the baroclinic, pseudo-thickness, and tracer updates. Setting `NTimeStepIteration` to one is permitted but leaves only the predictor and reduces the outer method to forward Euler.

An `SE-AB2` method is not implemented. It can be added later using the same barotropic-stepper interface and additional previous-tendency storage.


#### 3.2.1 Initialization

The default outer-iteration count is

$$
\text{NTimeStepIteration} = 2.
$$ (split-ntimestepiter-rk2)

Compute `NBtrSubcycle` as

$$
\text{NBtrSubcycles}
=
\max\left(1,
\left\lceil
\frac{\text{TimeStep}}{\text{BtrTimeStep}}
\right\rceil\right).
$$ (split-nbtrsubcycle)

`BtrTimeStep` is therefore an upper bound requested by the user. The effective
barotropic step is

$$
\Delta t_{\mathrm{bt}}
= \frac{\Delta t}{\text{NBtrSubcycles}}.
$$ (split-effective-btr-timestep)

Compute the barotropic velocity:

$$
\overline{{\bf u}}_e
\equiv
\frac{1}{[\tilde{H}_i]_e}
\sum_{k=0}^{K_{\max}}
[\tilde{h}_{i,k}]_e {\bf u}_{e,k}.
$$ (split-initial-barotropic-velocity)

On startup from initial conditions, the barotropic velocity of
Eq. {eq}`split-initial-barotropic-velocity` is formed with the edge pseudo
thickness $[\tilde{h}_{i,k}]_e$ taken as the arithmetic mean of the two
neighboring cell values, and the result is masked with the edge mask of the top
active layer so that dry edges carry no barotropic velocity. For a restart,
`NormalBarotropicVelocity` and `BarotropicPressureAnomaly` are preserved, while
`NormalBaroclinicVelocity` is made consistent with the full and barotropic
velocities read from the restart file through
${\bf u}'={\bf u}-\overline{\bf u}$.

Compute the baroclinic velocity:

$$
{\bf u}'_{e,k} = {\bf u}_{e,k} - \overline{{\bf u}}_e.
$$ (split-initial-baroclinic-velocity)

Before initializing the pressure split, the equation of state and the momentum
vertical auxiliary fields are computed from the time-level-0 state and tracers
through `AuxiliaryState::computeMomVertAux`. That call computes, in order, the
interface and midpoint pressures together with the column-integrated pseudo
thickness $\tilde H$, the specific volume, the depth-integrated and depth-mean
specific volume, the column-integrated geometric thickness $H=\rho_0 S$, the
geometric height, and the target thickness.

Compute the barotropic pressure:

$$
B_i = p_{i,K_{\max}+1} - p_{i}^{\text {surf}}.
$$ (split-initial-barotropic-pressure)

Here $p_{i,0}$ is the surface pressure and $p_{i,K_{\max}+1}$ is the pressure at
the interface below the bottom-most active layer. $B$ is stored once as a
diagnostic in the split-explicit scratch space and is not used again during time
stepping; the prognostic quantity is the anomaly $B'$.

Compute the barotropic pressure anomaly:

$$
B_i' = B_i - \rho_0 g b_i.
$$ (split-initial-barotropic-pressure-anomaly)

Time level 1 is then initialized from time level 0 for normal velocity,
baroclinic velocity, barotropic velocity, pseudo thickness, and barotropic
pressure anomaly; tracers are copied at the start of `doStep`. If
`ReinitSplitVelocity` is true, the velocity split is recomputed from
`NormalVelocity` at the beginning of each slow step, including a new depth mean
for `NormalBarotropicVelocity`. When it is false, the existing barotropic
velocity is retained and only baroclinic velocity is reconstructed from the
other two velocity fields.

Throughout Section 3.2 an asterisk marks a *provisional*, or working, value: the
most recent estimate of a field produced by the outer predictor-corrector
iteration and stored at time level 1. It is distinct from the time-$n$ value at
time level 0, which is held fixed for the whole slow step, and from the final
$n+1$ value. All tendencies and auxiliary variables are evaluated from these
starred quantities. At the start of a slow step the working state is a copy of
the time-$n$ state,

$$
{\bf u}^{*}_{e,k} = {\bf u}^n_{e,k},
\qquad
\tilde{W}^{*}_{i,k} = \tilde{W}^{n}_{i,k},
\qquad
\tilde{h}^{*}_{i,k} = \tilde{h}^{n}_{i,k},
\qquad
\varphi^{*}_{i,k} = \varphi^{n}_{i,k},
\qquad
p^{*}_{i,k} = p^{n}_{i,k},
\qquad
\tilde{z}^{*}_{i,k} = \tilde{z}^{n}_{i,k},
\quad \text{etc.}
$$ (split-initial-provisional-variables)

#### 3.2.2 Stage 1: Baroclinic velocity advance with long time step

This stage advances the baroclinic velocity $u'$ with the long time step and computes the barotropic forcing term $\overline{G}$.

Compute the baroclinic forcing plus the barotropic pressure-gradient contribution:

$$
\Gamma^*_{e,k} + [\overline{\alpha}^{*}_{i}]_e [\nabla B_i^{\prime *}]_e .
$$ (split-stage1-baroclinic-forcing)

Both terms are evaluated from the working state: the pseudo thickness, full
normal velocity, barotropic velocity, and barotropic pressure anomaly are all
taken at the working time level, and the momentum auxiliary variables are
recomputed from that state first. The barotropic pressure-anomaly gradient term
is gated on `SplitFactor`, because it exists solely to cancel the barotropic part of the
full pressure gradient and must be present whenever the mode split is active.
The horizontal-advection term uses the relative vorticity alone, with the
planetary contribution removed, so that the Coriolis acceleration can be
iterated separately below.

Compute the column-integrated pseudo thickness:

$$
\tilde{H}^{*}_i
= \sum_{k=0}^{K} \tilde{h}^{*}_{i,k}.
$$ (split-stage1-column-pseudo-thickness)

The non-Coriolis tendency is saved in `BaseVelocityTend`. Compute the Coriolis
term using a fixed-point treatment. For
$j = 0, \ldots, \text{NBclCoriolisIteration}-1$, with default value 2, restore
the saved tendency after the first iteration, add the Coriolis acceleration
computed from the latest working baroclinic velocity, recompute $\overline G$,
and update the baroclinic velocity from the fixed time-$n$ state:

$$
{\bf u}_{e,k}^{\prime\perp *}
=
\sum_{e'\in ECP(e)} \tilde{E}_{e,e'} f_{e'} {\bf u}_{e',k}^{\prime *}.
$$ (split-stage1-baroclinic-coriolis)

The Coriolis kernel shares its enable switch with the potential-vorticity
horizontal-advection term: when that tendency is disabled, the Coriolis
acceleration is skipped as well, in both the baroclinic and barotropic systems.

Before each update, form the provisional full-step baroclinic velocity

$$
{\bf u}_{e,k}^{\prime\mathrm{prov}}
= {\bf u}_{e,k}^{\prime n}
+ \Delta t\,{\cal R}_{e,k}^{*},
$$ (split-stage1-provisional-baroclinic-velocity)

where ${\cal R}^{*}$ is the current velocity tendency including the Coriolis
term. The implemented barotropic forcing is the pseudo-thickness-weighted
column mean of this provisional velocity divided by $\Delta t$:

$$
\overline G_e^{*}
= \frac{1}{\Delta t}
  \frac{\displaystyle\sum_k [\tilde h_k^*]_e
        {\bf u}_{e,k}^{\prime\mathrm{prov}}}
       {\displaystyle\sum_k [\tilde h_k^*]_e}.
$$ (split-stage1-gbar)

Here $[\tilde h_k^*]_e$ is the arithmetic mean of the two neighboring cell
pseudo thicknesses at the working time level, and ${\bf u}^{\prime n}$ is the
baroclinic velocity at time level 0, held fixed through the Coriolis iteration.
The result is stored in `BarotropicForcing` and is the $\overline G$ used
unchanged by every barotropic subcycle step in Stage 2.

The forcing is subtracted from every active layer of ${\cal R}^{*}$, ensuring
that the tendency passed to the baroclinic update has its depth-mean component
removed. The working baroclinic velocity is then advanced to the midpoint:

$$
{\bf u}_{e,k}^{\prime n+1/2}
=
{\bf u}_{e,k}^{\prime n}
+ \frac{\Delta t}{2}
\left(
-[f_v]_e {\bf u}_{e,k}^{\prime\perp *}
+ \Gamma^*_{e,k}
+ [\overline{\alpha}^{*}_{i}]_e [\nabla B_i^{\prime *}]_e
- \overline G_e^{*}
\right).
$$ (split-stage1-baroclinic-advance)

For `UnsplitRK2`, set

$$
\overline{G}^{*}_e = 0
$$ (split-stage1-gbar-unsplit)

and use one Coriolis iteration. After Stage 1, the implementation exchanges
only the working baroclinic-velocity halo.

#### 3.2.3 Stage 2: Barotropic velocity advance, explicitly subcycled

This stage advances $B'$ and $\overline{\bf {u}}$ as a coupled system through $2M$
predictor-corrector subcycle steps, where $M=\text{NBtrSubcycles}$ and
$\Delta t_{\mathrm{bt}}=\Delta t/M$. As in the MPAS-Ocean algorithm, the
subcycle steps span an extended $2\Delta t$ barotropic averaging window. For
`UnsplitRK2`, this stage is skipped, $\overline u=0$, and $u=u'$.

The discrete barotropic continuity equation is

$$
B_i^{\prime n+1}
= B_i^{\prime n}
- \Delta t
\left[
\nabla \cdot \left( [(B_i^{\prime n}+\rho_0 g b_i)]_e \overline{{\bf u}}_e^n \right)
\right]_i
- \Delta t\,\rho_0 g Q_i^n.
$$ (split-stage2-discrete-continuity)

The discrete barotropic momentum update is

$$
\overline{{\bf u}}^{n+1}_e
= \overline{{\bf u}}^{n}_e
+ \Delta t
\left(
- f_e \overline{{\bf u}}^{\perp n}_e
- [\overline{\alpha}_i^{*}\nabla B_i^{\prime n}]_e
+ \overline{G}^{*}_e
\right).
$$ (split-stage2-discrete-momentum)

Throughout this section a hat marks a quantity that lives inside the barotropic
subcycle and is updated on every subcycle step, as distinct from the unhatted
outer-iteration quantities carried by Stages 1 and 3.

At the start of an outer iteration, initialize the subcycle variables and
accumulators:

$$
\hat{\overline{{\bf u}}}^{n}_e = \overline{{\bf u}}^{n}_e,
\qquad
\hat{B}^{\prime n}_i = B_i^{\prime n},
\qquad
\overline{{\bf u}}_e^{\mathrm{sum}}=\overline{{\bf u}}_e^n,
\qquad F_e^{\mathrm{sum}}=0.
$$ (split-stage2-initialization)

The barotropic mass flux is not formed directly from
$(B'+\rho_0gb){\overline {\bf{u}}}$. Instead, the implementation uses the most recent
baroclinic edge pseudo thickness and a pressure-anomaly correction relative to
the provisional outer-iteration state $B'^{*}$:

$$
{\cal P}_e({\cal B}')
= \rho_0g\sum_k[\tilde h_k^*]_e
  + [{\cal B}'-B'^{*}]_e.
$$ (split-stage2-effective-pressure)

${\cal P}_e$ is written as a function of ${\cal B}'$ because it is evaluated
more than once per subcycle step with different arguments. ${\cal B}'$ is a placeholder
for whichever barotropic pressure anomaly is passed in; it is not a stored field
of its own. The $B'$ predictor, the $B'$ corrector, and the transport
accumulation each form a barotropic mass flux, so each evaluates ${\cal P}_e$,
passing a different combination of the three subcycle work arrays, `Cur`, `Pre`,
and `Cor`, which hold the anomaly at the start of the step, the predictor
output, and the corrector output; the combinations are written out in the
equations below.

Only the bracketed difference in Eq. {eq}`split-stage2-effective-pressure`
depends on ${\cal B}'$. The column sum $\sum_k[\tilde h_k^*]_e$ is identical in
all three evaluations: it is the flux pseudo thickness on edges produced by the
pseudo-thickness auxiliary state, so it already carries the configured centered
or upwind edge reconstruction, summed over the column. It is diagnosed once from
the provisional state before the subcycle loop begins and is held fixed for all
$2M$ subcycles. The barotropic fluxes formed with ${\cal P}_e$ do not update the
pseudo thickness, so the column sum does not change as the subcycles proceed.

The reference $B'^{*}$ is the provisional barotropic pressure anomaly at the
working time level on entry to Stage 2, in the same sense the asterisk carries
elsewhere in Section 3.2: $B'^n$ on the first outer iteration and the
$\rho_0 g(\tilde H^{*}-b)$ of Eq. {eq}`split-reset-barotropic-pressure-anomaly`
left by the previous outer iteration afterward. It is likewise held fixed and is
only overwritten by the final corrected anomaly after the last pass, so
${\cal P}_e$ measures the departure of the subcycled anomaly from the baroclinic
column mass. This is what enforces the [Hallberg and Adcroft (2009)](https://adcroft.github.io/assets/pdf/hallberg_adcroft_OM_2009.pdf) consistency
described in Section 3.2.4. Note that $B'^{*}$ is unhatted and so is an
outer-iteration quantity, whereas the hatted $\hat B^{\prime *}$ appearing in the
predictor-corrector equations below is the subcycle predictor output.

The edge value of the anomaly correction uses the same centered or upwind choice
selected for pseudo-thickness fluxes. With the upwind choice the upwind
direction is undefined where the normal velocity is exactly zero, and the
implementation then takes the larger of the two neighboring cell values.

The depth-mean specific volume $[\overline{\alpha}_i]_e$ is
`Eos::DepthMeanSpecificVolume` averaged to the edge. It is computed from the
working state during Stage 1 and is not updated during subcycling.

`Predictor-Corrector` is the only barotropic time stepper currently implemented, so
the subcycle steps below are specific to it. Section 4.2.3 gives the interface an
alternative barotropic stepper would implement.

For each predictor-corrector subcycle step $m=0,\ldots,2M-1$, apply the updates
below. Both velocity updates are multiplied by the edge mask of the top active
layer, and edges with no active layers are set to zero.

**$\overline{\bf {u}}$ predictor:**

$$
[\hat{\overline{{\bf u}}}^{*}_e]^{n+(m+1)/M}
=
[\hat{\overline{{\bf u}}}_e]^{n+m/M}
+
\frac{\Delta t}{M}
\left(
-f_e[\hat{\overline{{\bf u}}}_e^{\perp}]^{n+m/M}
-[\overline{\alpha}_i]^{*}_e[\nabla \hat{B}_i']^{n+m/M}_e
+ \overline{G}_e
\right).
$$ (split-stage2-u-predictor)

**$B'$ predictor:**

$$
[F^{*}_e]^{m+1}
=
{\cal P}_e\left([\hat{B}_i']^{n+m/M}\right)
\left(
(1-\gamma_1)[\hat{\overline{{\bf u}}}_e]^{n+m/M}
+
\gamma_1[\hat{\overline{{\bf u}}}^{*}_e]^{n+(m+1)/M}
\right),
$$ (split-stage2-b-predictor-flux)

and

$$
[\hat{B}_i^{\prime *}]^{n+(m+1)/M}
=
[\hat{B}_i']^{n+m/M}
+
\frac{\Delta t}{M}
[-\nabla \cdot [F^{*}_e]^{m+1}]_i.
$$ (split-stage2-b-predictor)


**$\overline{\bf {u}}$ corrector:**

$$
[\hat{\overline{{\bf u}}}_e]^{n+(m+1)/M}
=
[\hat{\overline{{\bf u}}}_e]^{n+m/M}
+
\frac{\Delta t}{M}
\left(
-f_e[\hat{\overline{{\bf u}}}_e^{\perp *}]^{n+(m+1)/M}
-
[\overline{\alpha}_i]^{*}_e
\nabla
\left(
(1-\gamma_2)[\hat{B}_i']^{n+m/M}
+
\gamma_2[\hat{B}_i^{\prime *}]^{n+(m+1)/M}
\right)_e
+
\overline{G}_e
\right).
$$ (split-stage2-u-corrector)

**$B'$ corrector:**

$$
[F_e]^{m+1}
=
{\cal P}_e\left(
(1-\gamma_2)[\hat{B}_i']^{n+m/M}
+
\gamma_2[\hat{B}_i^{\prime *}]^{n+(m+1)/M}
\right)
\left(
(1-\gamma_3)[\hat{\overline{{\bf u}}}_e]^{n+m/M}
+
\gamma_3[\hat{\overline{{\bf u}}}^{*}_e]^{n+(m+1)/M}
\right),
$$ (split-stage2-b-corrector-flux)

and

$$
[\hat{B}_i^{\prime \mathrm{cor}}]^{n+(m+1)/M}
=
[\hat{B}_i']^{n+m/M}
+
\frac{\Delta t}{M}
[-\nabla \cdot [F_e]^{m+1}]_i.
$$ (split-stage2-b-corrector)

On *hatted* subcycle quantities the asterisk denotes a different work array
depending on where it appears, so the mapping to the code is worth stating
explicitly. Starred hatted pressures are always the predictor output `Pre`.
Starred hatted velocities are the predictor output `Pre` in the $\overline {\bf{u}}$
predictor, in the $B'$ predictor flux, and in the Coriolis term of the
$\overline {\bf {u}}$ corrector; they are the corrector output `Cor` in the $B'$
corrector flux and in the transport accumulation. On unhatted quantities such as
$B'^{*}$, $\tilde h^{*}$, and $\overline{\alpha}^{*}$, the asterisk keeps its
Section 3.2 meaning of a provisional outer-iteration value, which is frozen for
the whole of Stage 2. The coefficients are

$$
\gamma_1=0.5333,\qquad \gamma_2=0.5333,\qquad \gamma_3=1.
$$ (split-stage2-gamma-values)

Accumulate the barotropic velocity and flux during subcycling:

$$
\overline{{\bf u}}_e^{\mathrm{sum}}
=
\sum_{m=0}^{2M}
[\hat{\overline{{\bf u}}}_e]^{n+m/M},
$$ (split-stage2-barotropic-velocity-accumulate)

and

$$
F_e^{\mathrm{sum}}
=
\sum_{m=0}^{2M-1}
[F_e]^{(m+1)/M}.
$$ (split-stage2-barotropic-flux-accumulate)

Compute the time average after subcycling:

$$
\overline{ {\bar {\bf u}}}_e^{\text{bt}}
=
\frac{1}{2M+1}
\overline{{\bf u}}_e^{\mathrm{sum}}
$$ (split-stage2-barotropic-velocity-average)

and

$$
\overline{F}_e^{\text{bt}}
=
\frac{1}{2M}F_e^{\mathrm{sum}}
$$ (split-stage2-barotropic-flux-average)

The flux accumulated for transport uses the corrected velocity and evaluates
${\cal P}_e$ with $(1-\gamma_2)\hat B^{\prime n+m/M}
+\gamma_2\hat B^{\prime\mathrm{cor}}$, which differs from the pressure used in
the $B'$ corrector flux above, where the second term is the predictor anomaly.
Both accumulations run over owned edges only.

#### 3.2.4 Barotropic-baroclinic coupling and barotropic pressure consistency

For the mode-split consistency of the barotropic pressure anomaly between $B'$ from the barotropic mode and $\rho_0 g(\tilde{H}-b)$ from the baroclinic mode, Omega follows the scheme from [Hallberg and Adcroft (2009)](https://adcroft.github.io/assets/pdf/hallberg_adcroft_OM_2009.pdf), as implemented in MPAS-Ocean.

The barotropic update of $B'$ is given by

$$
\frac{B_i^{\prime n+1}-B_i^{\prime n}}{\Delta t}
+
\left[
\nabla \cdot
\overline{\left([(B_i' + \rho_0 g b_i)]_e \overline{{\bf u}}_e\right)}^{\text{bt}}
\right]_i
= -\rho_0 g Q_i .
$$ (split-btr-update-consistency)

Here, $\overline{\varphi}^{\text{bt}}$ denotes a time-averaged quantity from the barotropic subcycles, and $n$ indicates the baroclinic time step.

The velocity correction ${\bf u}^{\text{co}}$ is written as

$$
{\bf u}_{e,k}^{\text{co}}
=
\left\{
\frac{\overline{F}_e^{\text{bt}}}{\rho_0g}
-
\sum_{k=0}^{K}
[\tilde{h}_i^{*}]_{e,k}
\left(
\overline{{\bar {\bf u}}}_e^{\text{bt}}
+{\bf u}_{e,k}^{\prime n+0.5}
\right)
\right\}
\bigg/
[\tilde{H}_i^{*}]_e ,
$$ (split-velocity-correction)

where $\overline{F}_e^{\text{bt}}$ is the time-averaged pressure transport
from Stage 2. The same correction is added uniformly to every active layer on
an edge. Here $[\tilde h_i^*]_{e,k}$ is the arithmetic mean of the two
neighboring cell pseudo thicknesses at the working time level, not the flux
pseudo thickness used by ${\cal P}_e$ in Stage 2.

The same kernel that forms ${\bf u}^{\text{co}}$ also writes the uncorrected sum
$\overline{{\bar{\bf u}}}^{\text{bt}}_e + {\bf u}^{\prime n+0.5}_{e,k}$ into the
working `NormalVelocity`, which is what the Stage-3 auxiliary variables and
tendencies see. On the final outer iteration this value is replaced later by the
$n+1$ reconstruction of Eq. {eq}`split-final-full-velocity`.

The transport velocity $u^{\text{tr}}$ is defined as

$$
{\bf u}^{\text{tr}}_{e,k}
=
\overline{{\bar {\bf u}}}_e^{\text{bt}}
+{\bf u}_{e,k}^{\prime n+0.5}
+{\bf u}_{e,k}^{\text{co}}.
$$ (split-transport-velocity)

The transport velocity is used to compute vertical transport velocity and horizontal transport for both pseudo thickness and tracers.

For the unsplit algorithm, the above processes are skipped except that

$$
{\bf u}_{e,k}^{\text{tr}}
= {\bf u}_{e,k}^{\prime n+0.5},
$$ (split-unsplit-transport-velocity)

where

$$
{\bf u}_{e,k}^{\prime n+0.5} = {\bf u}_{e,k}^{n+0.5}.
$$ (split-unsplit-midpoint-relation)

#### 3.2.5 Stage 3: Update tracers and diagnostics

Compute the pseudo-thickness and tracer auxiliary variables from
${\bf u}_{e,k}^{\text{tr}}$ through
`AuxiliaryState::computePseudoThicknessTracerAux`, which also refreshes the
vertical momentum auxiliary variables, and therefore the pressure, specific
volume, column thicknesses, geometric height, and target thickness, from the
working pseudo thickness and tracers.

Compute $\tilde{W}_{i,k}^{*}$ using ${\bf u}_{e,k}^{\text{tr}}$.

Compute pseudo thickness tendencies using ${\bf u}_{e,k}^{\text{tr}}$:

$$
{\cal T}_{\tilde h}
=
-
\nabla \cdot
\left([\tilde{h}_{k}^{*}]_e {\bf u}_{e,k}^{\text{tr}}\right)
-
\left(
[\tilde{W}_{tr}]^{\text{top}}_{k}
-
[\tilde{W}_{tr}]^{\text{top}}_{k+1}
\right).
$$ (split-stage3-pseudo-thickness-update)

Compute tracer tendencies using ${\bf u}_{e,k}^{\text{tr}}$:

$$
{\cal T}_{\varphi}
=
-
\nabla \cdot
\left([
\tilde{h}_{i,k}^{*}]_e
[\varphi_{i,k}^{*}]_e
{\bf u}_{e,k}^{\text{tr}}
\right)
-
\left(
[\varphi^{*}\tilde{W}_{tr}]^{\text{top}}_{i,k}
-
[\varphi^{*}\tilde{W}_{tr}]^{\text{top}}_{i,k+1}
\right).
$$ (split-stage3-tracer-update)

On a non-final outer iteration, pseudo thickness is advanced to the midpoint
with $\Delta t/2$:

$$
\tilde h^{*}
=
\tilde h^n+\frac{\Delta t}{2}{\cal T}_{\tilde h}.
$$ (split-reset-provisional-thickness)

This is the same midpoint state as the MPAS-Ocean form
$\tilde h^{*}=\tfrac{1}{2}(\tilde h^{n}+\tilde h^{n+1})$, because
$\tilde h^{n+1}=\tilde h^{n}+\Delta t\,{\cal T}_{\tilde h}$. Writing it as a
half-step tendency update lets the same kernel that performs the full-step update
on the final iteration be reused on the intermediate iterations, and avoids
storing $\tilde h^{n+1}$ and averaging it in an extra loop.

For each tracer, the implementation first forms the conservative
provisional end concentration:

$$
\varphi^{\mathrm{end}}
=\frac{\tilde h^n\varphi^n+\Delta t\,{\cal T}_\varphi}
       {\tilde h^n+\Delta t\,{\cal T}_{\tilde h}},
$$ (split-stage3-provisional-tracer)

and stores $\varphi^*=(\varphi^n+\varphi^{\mathrm{end}})/2$. On the final
outer iteration, it retains the full-step conservative updates shown above.

#### 3.2.6 Reset variables

The reconstruction of the full velocity and the reset of $\tilde H$ and $B'$ are
performed by a single routine at the end of Stage 3 on every outer iteration.
Only the velocity reconstruction differs between an intermediate iteration and
the final one; the column pseudo thickness and the barotropic pressure anomaly
are recomputed the same way in both cases. For unsplit stepping the routine
reconstructs the velocity and returns without touching $B'$.

If iterating, reset the provisional variables as follows:

$$
{\bf u}^{\prime *} = {\bf u}^{\prime n+0.5} \quad \text{from Stage 1},
$$ (split-reset-baroclinic-velocity)

$$
\overline{{\bf u}}^{*} = \overline{{\bar {\bf u}}}^{\text{bt}} \quad \text{from Stage 2},
$$ (split-reset-barotropic-velocity)

$$
{\bf u}^{*} = \overline{{\bf u}}^{*} + u^{\prime *},
$$ (split-reset-full-velocity)

with pseudo thickness and tracer concentration reset by
Eqs. {eq}`split-reset-provisional-thickness` and
{eq}`split-stage3-provisional-tracer` of Section 3.2.5. The column totals then
follow:

$$
\tilde{H}^{*} = \sum_{k=0}^{K} \tilde{h}_k^{*},
$$ (split-reset-column-pseudo-thickness)

and

$$
B^{\prime *} = \rho_0 g\left(\tilde{H}^{*}-b\right).
$$ (split-reset-barotropic-pressure-anomaly)

The full working velocity is the sum of the working split velocities. The state
and tracer halos are exchanged before the next outer iteration.

After the final iteration, Stage 1 still stores the working baroclinic velocity
at the midpoint:

$$
{\bf u}^{\prime *}={\bf u}^{\prime n+1/2}.
$$ (split-final-baroclinic-velocity)

$$
\overline{{\bf u}}^{n+1} = \overline{{\bar {\bf u}}}^{\text{bt}} \quad \text{from Stage 2},
$$ (split-final-barotropic-velocity)

$$
{\bf u}^{n+1}
= \overline{{\bf u}}^{n+1}
  +2{\bf u}^{\prime n+1/2}-{\bf u}^{\prime n},
$$ (split-final-full-velocity)

$$
\tilde{h}^{n+1}_{i,k}
\quad\text{and}\quad
\varphi^{n+1}_{i,k}
\quad \text{are the full-step values retained from Stage 3},
$$ (split-final-thickness-tracer)

$$
\tilde{H}^{n+1} = \sum_{k=0}^{K} \tilde{h}_k^{n+1},
$$ (split-final-column-pseudo-thickness)

and

$$
B^{\prime n+1} = \rho_0 g\left(\tilde{H}^{n+1}-b\right).
$$ (split-final-barotropic-pressure-anomaly)

The state and tracer time levels are rotated once, making the completed
working state the new time level 0. Kinetic energy, velocity divergence, and the
other kinetic auxiliary variables are then recomputed from the rotated $n+1$
velocity so that history output and state validation see diagnostics consistent
with the reconstructed velocity. Implicit vertical mixing is applied next if
either velocity or tracer vertical mixing is enabled, followed by a state and
tracer halo exchange. Finally the state is validated and the clock and
persistent step count advance.

### 3.3 Unsplit time stepping algorithm

`UnsplitRK2` uses the same `SplitExplicitRK2Stepper` implementation with
`SplitFactor=0`. The stored baroclinic velocity is set equal to full velocity,
barotropic velocity is zero, barotropic forcing is zero, and Stage 2 is
skipped. Barotropic configuration keys are ignored.

#### 3.3.1 Initialization

Compute the pressure $p$ and the other vertical momentum auxiliary variables.

Set the velocity split trivially, at both time levels and for both a startup and
a restart:

$$
\overline{{\bf u}}_e = 0,
\qquad
{\bf u}'_{e,k} = {\bf u}_{e,k}.
$$ (unsplit-velocity-split)

Barotropic pressure and pressure anomaly are not initialized, since $B'$ is
never referenced when `SplitFactor` is zero.

Prepare variables before the first iteration:

$$
{\bf u}^{*}_{e,k} = {\bf u}^n_{e,k},
\qquad
\tilde{W}^{*}_{i,k} = \tilde{W}^{n}_{i,k},
\qquad
\tilde{h}^{*}_{i,k} = \tilde{h}^{n}_{i,k},
\qquad
\varphi^{*}_{i,k} = \varphi^{n}_{i,k},
\qquad
p^{*}_{i,k} = p^{n}_{i,k},
\qquad
\tilde{z}^{*}_{i,k} = \tilde{z}^{n}_{i,k},
\quad \text{etc.}
$$ (unsplit-initial-provisional-variables)

#### 3.3.2 Stage 1: Velocity advance

This stage advances the full velocity $u$.

Compute $\Gamma_{e,k}^{*}$:

$$
\Gamma_{e,k}^{*}
=
-[{\bf k}\cdot\nabla\times {\bf u}_{e,k}]_e {\bf u}_{e,k}^{\perp}
-[\nabla K]_e
-
\frac{1}{[\tilde{h}_{i,k}]_e}
\left[\tilde{W}_{tr}\frac{\partial U}{\partial \tilde{z}}\right]_{e,k}
-(\alpha\nabla p+\nabla\Phi)_{e,k}
+[D_h^{\bf u}]_{e,k}
+[F_s^{\bf u}]_{e,k}
+[F_b^{\bf u}]_{e,k}.
$$ (unsplit-gamma-definition)

As in Section 3.1, the vertical mixing term is absent because it is applied
implicitly after the outer iterations, and the barotropic pressure-anomaly
gradient term of Eq. {eq}`split-stage1-baroclinic-forcing` is absent because
`SplitFactor` is zero.

Compute the column-integrated pseudo thickness:

$$
\tilde{H}_i^{*}
= \sum_{k=0}^{K} \tilde{h}^{*}_{i,k}.
$$ (unsplit-column-pseudo-thickness)

The unsplit configuration forces `NBclCoriolisIteration=1`. Compute the
Coriolis term from the current working full velocity:

$$
{\bf u}_{e,k}^{\perp *}
=
\sum_{e'\in ECP(e)} \tilde{E}_{e,e'} f_{e'} {\bf u}_{e',k}^{*}.
$$ (unsplit-coriolis)

Set

$$
\overline{G}_e^{*} = 0,
$$ (unsplit-gbar)

so the barotropic forcing is neither computed nor subtracted from the tendency.
The shared update routine then advances the working velocity directly to the
midpoint with the half coefficient:

$$
{\bf u}_{e,k}^{n+0.5}
=
{\bf u}_{e,k}^{n}
+
\frac{\Delta t}{2}
\left(
-[f_v]_e {\bf u}_{e,k}^{\perp *}
+
\Gamma_{e,k}^{*}
\right).
$$ (unsplit-velocity-advance)

The full-step velocity ${\bf u}^{n+1}$ is not formed here; it is recovered from
the midpoint at the end of the final outer iteration by
Eq. {eq}`unsplit-final-velocity`.

#### 3.3.3 Stage 2: Barotropic velocity advance, explicitly subcycled

For the unsplit time stepper, $\overline{{\bf u}}=0$. This stage is skipped.

#### 3.3.4 Stage 3: Update tracers and diagnostics

Compute ${\bf u}_{e,k}^{\text{tr}}$:

$$
{\bf u}_{e,k}^{\text{tr}}
= {\bf u}_{e,k}^{n+0.5}.
$$ (unsplit-transport-velocity)

There is no barotropic mode to reconcile, so the velocity correction
${\bf u}^{\text{co}}$ of Eq. {eq}`split-velocity-correction` is identically
zero and the transport velocity is simply the midpoint velocity from Stage 1.

From this point Stage 3 is identical to the split case, because it depends on
the mode split only through ${\bf u}^{\text{tr}}$. The same routine computes
the pseudo-thickness and tracer auxiliary variables, diagnoses
$\tilde{W}_{tr}$, and advances pseudo thickness and tracers, applying the
half-step midpoint construction on a non-final outer iteration and retaining
the full-step conservative update on the final one. Section 3.2.5 gives the
equations.

#### 3.3.5 Reset variables

Only the velocity handling differs from the split case, since the mode split
enters Section 3.2.6 solely through the velocity. If iterating, the working full
velocity is simply the Stage 1 midpoint,

$$
{\bf u}^{*} = {\bf u}^{n+0.5} \quad \text{from Stage 1},
$$ (unsplit-reset-velocity)

and after the final iteration the physical full-step velocity is reconstructed
from it with no barotropic contribution:

$$
{\bf u}^{n+1}=2{\bf u}^{n+1/2}-{\bf u}^{n}.
$$ (unsplit-final-velocity)

Everything else follows Section 3.2.6 unchanged: the pseudo-thickness and tracer
resets, the column totals $\tilde H^{*}$ and $\tilde H^{n+1}$, the working-state
and tracer halo exchange between outer iterations, and the closing sequence of
time-level rotation, kinetic diagnostics, implicit vertical mixing, state
validation, and clock advancement. Because `SplitFactor` is zero, the routine
skips the $B'$ reset of
Eq. {eq}`split-reset-barotropic-pressure-anomaly`.

## 4. Design

The split-explicit time stepping implementation consists of a baroclinic stepper with a barotropic subcycle stage.  The same
stepper also supports an unsplit mode by setting the split forcing factor to zero and skipping the barotropic subcycling stage.

The implementation is spread across the following files.

| File | Contents |
| --- | --- |
| `timeStepping/SplitExplicitTypes.h` | `SplitExplicitConfig` and `SplitExplicitScratch` |
| `timeStepping/SplitExplicitInit.{h,cpp}` | Config reading, scratch allocation, velocity split and recombination, barotropic pressure initialization |
| `timeStepping/SplitExplicitRK2Stepper.{h,cpp}` | Outer iteration, Stage 1, Stage 3, transport velocity, time-level management |
| `timeStepping/SplitExplicitBarotropicPCStepper.{h,cpp}` | Stage 2 predictor-corrector barotropic subcycle |
| `timeStepping/TimeStepper.{h,cpp}` | New stepper types and the `initializeStateFromInput` hook |
| `ocn/OceanState.{h,cpp}` | The three new split prognostic fields |
| `ocn/Tendencies.{h,cpp}` | Baroclinic velocity tendency and standalone Coriolis acceleration |
| `ocn/TendencyTerms.h` | `CoriolisAccelerationOnEdge` and the split-specific operator overloads |
| `ocn/AuxiliaryState.{h,cpp}` | `computePseudoThicknessTracerAux` and transport-velocity overloads |
| `ocn/Eos.{h,cpp}` | Depth-integrated and depth-mean specific volume |
| `ocn/VertCoord.{h,cpp}` | Column-integrated pseudo and geometric thickness |
| `ocn/OceanInit.cpp` | Calls `initializeStateFromInput` after input is read |

### 4.1 Data types and parameters

The slow `TimeStepper` and `TimeStep` options are read from `TimeIntegration`.
Split-specific options are read from its optional `ModeSplitShare` subgroup.
Missing keys retain the defaults in `SplitExplicitConfig`.

#### 4.1.1 Parameters

- `TimeStepper`: `SplitExplicitRK2` selects mode-split stepping and
  `UnsplitRK2` selects its unsplit counterpart. These are the only recognized
  split-stepper names.
- `NTimeStepIteration`: Number of outer predictor-corrector iterations. The
  default is 2 and values must be positive. Every iteration uses the full
  `TimeStep`; one iteration reduces the outer scheme to forward Euler.
- `NBclCoriolisIteration`: Number of centered Coriolis iterations used in the
  baroclinic velocity update.  The default is 2 for split-explicit stepping and
  1 for unsplit stepping.
- `BtrTimeStep`: Requested maximum barotropic subcycle time step. It must be
  positive and is ignored by `UnsplitRK2`.
- `BtrTimeStepper`: Barotropic subcycle algorithm.  The initial implementation
  supports `Predictor-Corrector`.
- `SplitFactor`: Internal factor multiplying the barotropic pressure-anomaly
  contribution to the baroclinic velocity tendency.  It is 1 for
  `SplitExplicitRK2` and 0 for unsplit stepping.
- `ReinitSplitVelocity`: If true, recompute barotropic and baroclinic velocity
  from `NormalVelocity` at the beginning of every slow step. The default is
  false.
- `NBtrSubcycles`: Internal count computed as
  `max(1, ceil(TimeStep / BtrTimeStep))` when the stepper is constructed, and
  only for `SplitExplicitRK2`. `TimeStepper::changeTimeStep` does not currently
  recompute this count, so a stepper whose time step is changed after
  construction keeps its original subcycle count and shrinks
  $\Delta t_{\mathrm{bt}}$ proportionally.

If `ModeSplitShare` or `BtrTimeStep` is omitted from a user configuration, the
requested barotropic step defaults to the slow `TimeStep`, giving one nominal
subcycle and two predictor-corrector subcycle steps. `configs/Default.yml` supplies the
whole subgroup, so the fallback applies only when a run configuration replaces
rather than extends the defaults.

The `ModeSplitShare` settings shipped in `configs/Default.yml`, together with a
`SplitExplicitRK2` selection, are:

```yaml
TimeIntegration:
  TimeStepper: SplitExplicitRK2
  TimeStep: 0000_00:10:00
  ModeSplitShare:
    BtrTimeStepper: Predictor-Corrector
    BtrTimeStep: 0000_00:00:20
    NTimeStepIteration: 2
    NBclCoriolisIteration: 2
    ReinitSplitVelocity: false
```

#### 4.1.2 Class/structs/data types

The types below are grouped into those the split-explicit time stepper introduces
and those that already exist in Omega. For the latter, only the additions this
time stepper requires are listed; no existing member is removed or renamed, and the
other time steppers are unaffected.

**New types**

- `SplitExplicitConfig`: Holds split-explicit configuration, including
  barotropic stepper selection, subcycle count, time-step iteration count,
  Coriolis iteration count, and `SplitFactor`.
- `SplitExplicitScratch`: Holds temporary arrays for barotropic subcycling and
  baroclinic updates:
  - `NormalBarotropicVelocitySubcycle{Cur,Pre,Cor}` and
    `BarotropicPressureAnomalySubcycle{Cur,Pre,Cor}`, the three subcycle work arrays
    per field described in Section 3.2.3.
  - `BarotropicPressure`, the diagnostic $B$ formed at initialization.
  - `BarotropicForcing`, the $\overline G$ produced by Stage 1 and consumed by
    Stage 2.
  - `BarotropicFlux`, the time-averaged pressure transport
    $\overline F^{\mathrm{bt}}$ produced by Stage 2.
  - `BaroclinicPseudoThicknessEdge`, the column sum of flux pseudo thickness on
    edges used by ${\cal P}_e$.
  - `BaseVelocityTend`, the saved non-Coriolis baroclinic velocity tendency.
  - `NormalTransportVelocity`, the corrected transport velocity
    ${\bf u}^{\mathrm{tr}}$ used by Stage 3.

  All arrays are sized over the full halo extent and zeroed at allocation. The
  scratch struct is a `mutable` member of the stepper, since `doStep` is
  `const`.
- `SplitExplicitInit`: Provides initialization utilities for reading
  split-explicit options, allocating scratch arrays, splitting full normal
  velocity into barotropic and baroclinic parts, combining the split velocity,
  and initializing barotropic pressure and pressure anomaly.
- `SplitExplicitRK2Stepper`: Implements the RK2 time step, including the
  baroclinic velocity stage, optional barotropic subcycle stage, thickness and
  tracer stage, halo exchanges, and time-level rotation.
- `SplitExplicitBarotropicPCStepper`: Implements the predictor-corrector
  barotropic stage interface.  This class owns the details of the barotropic
  algorithm so that `SplitExplicitRK2Stepper` can call a generic stage-2
  function.

**Additions to existing classes**

- `TimeStepper`: Gains the `SplitExplicitRK2` and `UnsplitRK2` enumerators, the
  matching strings in `getTimeStepperFromStr`, the construction cases in
  `TimeStepper::create`, and a virtual `initializeStateFromInput` whose base
  implementation is an empty no-op so that the other steppers are unaffected.
- `OceanState`: Stores the additional split-explicit prognostic fields:
  `NormalBaroclinicVelocity` (edges, layered), `NormalBarotropicVelocity`
  (edges, single level), and `BarotropicPressureAnomaly` (cells, single level).
  All three are added to the `Restart` field group but not to the `State` group,
  so they are written to and read from restart files without appearing in default
  history output.
- `Tendencies`: Exposes `computeCoriolisAccelerationOnEdge` in layered and
  single-layer forms, so that the Coriolis term can be applied on
  its own inside the Stage-1 iteration and by the barotropic subcycle. Both
  forms return without acting when the potential-vorticity tendency is disabled.
- `TendencyTerms`: Adds a `CoriolisAccelerationOnEdge` operator with layered and
  single-layer call operators; an overload of `PotentialVortHAdvOnEdge` that
  advects with relative vorticity only, leaving the planetary part to the
  Coriolis operator.
- `AuxiliaryState`: Adds `computePseudoThicknessTracerAux`, which computes just
  the pseudo-thickness and tracer auxiliary variables needed by Stage 3, in
  time-level and explicit-velocity-array forms.
- `Eos`: Stores `DepthIntegSpecificVolume`, the column integral of specific
  volume weighted by pseudo thickness, and `DepthMeanSpecificVolume`, that
  integral divided by the column-integrated pseudo thickness. Both are computed
  by `computeDepthIntegratedSpecificVolume` and are registered as `Eos` fields.
  `DepthMeanSpecificVolume` is the $\overline{\alpha}$ used by both the
  barotropic momentum equation and the barotropic pressure-anomaly contribution
  to the baroclinic velocity tendency.
- `VertCoord`: Stores `TotalPseudoThickness`, the column sum $\tilde H$, and
  `TotalGeometricThickness`, the column geometric depth $H=\rho_0 S$ of
  Eq. {eq}`split-column-geometric-thickness`.
- `OceanInit`: `initStateForTimeStepper` calls
  `TimeStepper::initializeStateFromInput` after the initial-state or restart
  read and after `initUpdateHaloAndHostArrays`, then exchanges the state halo and
  copies to the host.

### 4.2 Methods

The implemented `doStep` control flow is shown below. Pacer timer
instrumentation, and the communicator used only by its timing barriers, are
omitted for readability.

```c++
void SplitExplicitRK2Stepper::doStep(OceanState *State,
                                     TimeInstant &SimTime) const {

   const int CurLevel  = 0;
   const int NextLevel = 1;

   Array3DReal CurTracerArray  = Tracers::getAll(CurLevel);
   Array3DReal NextTracerArray = Tracers::getAll(NextLevel);

   // Initialize NextLevel from CurLevel, optionally recomputing the split
   initializeNextState(State, CurLevel, NextLevel, SEConfig.SplitFactor,
                       SEConfig.ReinitSplitVelocity);
   deepCopy(NextTracerArray, CurTracerArray);

   const TimeInstant StageTime = SimTime;
   for (I4 TimeStepIteration = 0;
        TimeStepIteration < SEConfig.NTimeStepIteration; ++TimeStepIteration) {

      const bool FinalIteration =
          TimeStepIteration + 1 == SEConfig.NTimeStepIteration;

      // The first iteration evaluates the momentum right-hand side at the n
      // state copied into NextLevel; every later iteration sees the midpoint
      // state left by its predecessor, so time-dependent terms are sampled at
      // n+1/2 to keep the predictor-corrector second order.
      const TimeInstant VelStageTime =
          TimeStepIteration == 0 ? StageTime : StageTime + 0.5 * TimeStep;

      // Stage 1: Baroclinic velocity advance, with long time step
      doBaroclinicVelocityUpdate(State, NextTracerArray, CurLevel, NextLevel,
                                 VelStageTime, TimeStep);

      MeshHalo->exchangeFullArrayHalo(
          State->getNormalBaroclinicVelocity(NextLevel), OnEdge);

      if (SEConfig.SplitFactor != 0._Real) {
         // Stage 2: Barotropic velocity advance, explicitly subcycled
         doBarotropicVelocityUpdate(State, CurLevel, NextLevel,
                                    StageTime + 0.5 * TimeStep, TimeStep);
      }

      // Physical total velocity and the corrected transport velocity
      computeTransportVelocity(State, NextLevel);

      // Stage 3: Update thickness, tracers, other diagnostics
      doThicknessTracerUpdate(State, CurTracerArray, NextTracerArray, CurLevel,
                              NextLevel, StageTime, TimeStep, FinalIteration);

      if (!FinalIteration) {
         State->exchangeHalo(NextLevel);
         MeshHalo->exchangeFullArrayHalo(NextTracerArray, OnCell);
      }
   }

   State->updateTimeLevels();
   Tracers::updateTimeLevels();

   // Refresh the kinetic auxiliary variables from the completed n+1 velocity,
   // an inline parallelForOuter over KineticAux::computeVarsOnCell

   // Apply implicit vertical mixing, then re-exchange the affected halos
   VertMix *VMix = VertMix::getInstance();
   if (VMix->VelVertMixSetup.Enabled or VMix->TracerVertMixSetup.Enabled) {
      VMix->VertMixImplicit(State, AuxState, Tracers::getAll(CurLevel),
                            Tracers::getNumTracers(), CurLevel);
      State->exchangeHalo(CurLevel);
      Tracers::exchangeHalo(CurLevel);
   }

   validateOceanState(State, AuxState, VertCoord::getDefault(), CurLevel);

   StepClock->advance();
   SimTime = StepClock->getCurrentTime();
   ++StepCount;
}
```

The stepper also implements `finalizeInit`, which validates the tendency, mesh,
vertical coordinate, and halo pointers, allocates the scratch arrays, and binds
the barotropic stage callback for the configured `BtrTimeStepper`.

#### 4.2.1 Initialization

After the initial condition or restart state has been read, `OceanInit`
initializes the split-explicit state through
`TimeStepper::initializeStateFromInput`.  For `SplitExplicitRK2Stepper`, this
initialization performs the following operations:

1. Compute momentum vertical auxiliary variables needed for pressure.
2. For a non-restart `SplitExplicitRK2` run, compute barotropic pressure and
   pressure anomaly from pressure-interface, surface-pressure, and bottom-depth
   fields, then initialize all three barotropic subcycle pressure work arrays.
3. Initialize velocity split fields. For a non-restart run with
   `SplitExplicitRK2`, full normal velocity is split into barotropic and
   baroclinic components. For a restart run, barotropic velocity is preserved
   and baroclinic velocity is reconstructed as full minus barotropic velocity.
   For unsplit stepping, barotropic velocity is set to zero and baroclinic
   velocity is set equal to full normal velocity.
4. Initialize time level 1 from time level 0 so that the first time step starts
   from a consistent pair of time levels.

`OceanInit` then exchanges the time-level-0 state halo and copies it to the
host, so that the newly split velocity fields are consistent across ranks and on
the host before any output or time stepping.

At the beginning of every full time step, time level 1 is refreshed from time
level 0 for the state fields advanced by the split-explicit scheme, and tracer
time level 1 is copied from tracer time level 0.  The state at time level 0 is
kept fixed during the internal time-step iterations.  Time level 1 is updated
repeatedly and provides the most recent state for tendency computations.

This refresh is also where `ReinitSplitVelocity` acts. When it is true and the
mode split is active, the velocity split is recomputed from the time-level-0
`NormalVelocity`, including a new depth mean for `NormalBarotropicVelocity`.
When it is false, the stored barotropic velocity is kept and only the baroclinic
velocity is rebuilt as ${\bf u}-\overline{\bf u}$, which preserves the barotropic
state carried forward from the previous step.

#### 4.2.2 Stage 1: baroclinic velocity

Stage 1 computes the baroclinic velocity tendency at the current iteration
state and advances `NormalBaroclinicVelocity` by half of the full slow time
step. The first outer iteration samples at $t^n$ and every later iteration
samples time-dependent forcing at $t^{n+1/2}$. The non-Coriolis tendency is
saved in `SplitExplicitScratch` before the
centered Coriolis iteration begins.

The centered Coriolis treatment repeats `NBclCoriolisIteration` times.  Each
iteration restores the saved non-Coriolis tendency, adds the Coriolis
acceleration computed from the updated baroclinic velocity, and updates
`NormalBaroclinicVelocity` again from the fixed time-level-0 base state.  After
stage 1, only the updated baroclinic velocity halo is exchanged.

#### 4.2.3 Stage 2: barotropic subcycling

Stage 2 advances the barotropic velocity and barotropic pressure anomaly using
the configured barotropic stepper.  The RK2 stepper delegates this work through
a generic barotropic-stage callback so that additional barotropic algorithms can
be added without complicating the RK2 stage logic.

For `SplitExplicitRK2`, the barotropic stage is active and currently uses the
predictor-corrector barotropic stepper.  For unsplit stepping, `SplitFactor` is
zero and stage 2 is skipped. The predictor-corrector implementation owns three
work arrays (`Cur`, `Pre`, and `Cor`) for both barotropic velocity and pressure
anomaly. It exchanges the `Cur` work arrays once per subcycle step, evaluates the
remaining kernels over progressively smaller halo ranges so that no further halo
exchange is needed within a subcycle step, accumulates the corrected velocity and
the pressure transport, and exchanges the resulting time-averaged fields after
the last subcycle step.

The callback signature takes only the state, the two time levels, the stage
time, and the stage time step; the barotropic stepper reaches everything else it
needs, including the scratch arrays, the configuration, the mesh, the halo, the
vertical coordinate, and the equation of state, through the arguments the RK2
stepper binds into the lambda. Adding an implicit or Adams--Bashforth barotropic
algorithm therefore means adding an enumerator, a class implementing the same
`doBarotropicVelocityUpdate` entry point, and a branch in
`initBarotropicStepper`, with no change to the RK2 stage logic.

#### 4.2.4 Stage 3: thickness and tracers

Before Stage 3, `computeTransportVelocity` reconstructs the full normal
velocity from the split fields and adds a depth-uniform correction to a
separate `NormalTransportVelocity` scratch field. Stage 3 computes
pseudo-thickness and tracer auxiliary variables with this corrected transport,
diagnoses vertical pseudo velocity, and advances pseudo thickness and tracers.
These Stage-3 tendencies are sampled at $t^{n+1/2}$ on every outer iteration.

During intermediate time-step iterations, full normal velocity is reconstructed
as the sum of barotropic and baroclinic velocity at time level 1.  On the final
iteration, full normal velocity is reconstructed at time level $n+1$ using

$$
{\bf u}^{n+1} = {\overline {\bar {\bf u}}}^{\mathrm{bt}}
        + 2 {\bf u}'^{n+1/2}
        - {\bf u}'^{n}.
$$

This final reconstruction follows the MPAS-Ocean split-explicit update: the
baroclinic velocity stored at time level 1 is the midpoint value, but the full
normal velocity must be available at the final time level for output,
diagnostics, restart, and the next time step.

#### 4.2.5 Time-level management

The internal time-step iterations do not rotate Ocean state or tracer time
levels.  Instead, time level 0 remains the base state for RK2 updates and time
level 1 is updated in place.  Halo exchanges are performed on time level 1
between internal iterations.  After the final iteration, `updateTimeLevels` is
called once for the Ocean state and tracers, making the completed time-level-1
state the new time-level-0 state. `NormalBaroclinicVelocity`,
`NormalBarotropicVelocity`, and `BarotropicPressureAnomaly` participate in
`OceanState::exchangeHalo` and in restart I/O, but are registered only in the
`Restart` field group. The physical `NormalVelocity` and `PseudoThickness`
remain in the normal state/history group.

Because the three split fields are in the restart group, a restart carries the
barotropic velocity and pressure anomaly forward exactly, and Section 4.2.1
rebuilds only the baroclinic velocity from them. The alternative, recomputing
the split from `NormalVelocity` alone, would discard the barotropic pressure
anomaly and break restart reproducibility.

## 5. Verification and testing

### 5.1 Unit testing

`TimeStepperTest` exercises both new time-stepper types with a velocity-decay tendency. It initializes each stepper through
`initializeStateFromInput`, performs time-step refinement, and checks the
observed convergence rate. The current expected rates in the implementation
are first order for `SplitExplicitRK2` (tolerance 0.15) and second order for
`UnsplitRK2` (tolerance 0.1). This records the behavior of the current code and
does not claim second-order convergence for the complete split algorithm.

The same test creates `SplitExplicitRK2` and `UnsplitRK2` without a stop time,
verifies that no end alarm is present, calls `doStep` twice, and checks that the
persistent step count increments across calls. Test setup overrides
`BtrTimeStep` in the in-memory configuration so that the coarsest refinement
level runs ten barotropic subcycles per baroclinic step, exercising the subcycle
loop rather than the degenerate single-subcycle case. Because
`changeTimeStep` does not recompute `NBtrSubcycles`, the finer refinement level
runs the same ten subcycles at half the barotropic step.

The supporting operators added for the split scheme have their own unit tests:

- `TendencyTermsTest` adds `testCoriolisAccelerationOnEdge`, which checks both
  the layered and single-layer forms of `CoriolisAccelerationOnEdge` against an
  independently computed tangential reconstruction.
- `EosTest` adds `testDepthIntegratedSpecificVolume`, which checks the column
  integral against a constant specific volume and thickness.
- `VertCoordTest` checks `TotalPseudoThickness` after `computePressure` and
  `TotalGeometricThickness` after `computeTotalGeometricThickness` against
  analytic column values.

Additional focused unit tests should cover velocity split and recombination,
barotropic-pressure initialization, ceiling-based subcycle selection,
centered and upwind anomaly interpolation, the transport-velocity correction,
restart preservation of the three split prognostic fields, and rejection of halo
widths smaller than three.

### 5.2 Polaris tests

The following end-to-end verification cases are available in Polaris:

- Overflow
- Baroclinic channel
- Realistic global ocean
