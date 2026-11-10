(omega-design-pressure-grad-high-order)=
# Higher-Order Horizontal Pressure Gradient

**Table of Contents**
1. [Overview](#1-overview)
2. [Requirements](#2-requirements)
3. [Algorithmic Formulation](#3-algorithmic-formulation)
4. [Design](#4-design)
5. [Verification and Testing](#5-verification-and-testing)
6. [Supporting evidence](#6-supporting-evidence)

## 1 Overview

This document specifies a more accurate discretization of the horizontal pressure gradient force
(PGF) in Omega's momentum equation, delivered in two phases: a consistent second-order scheme, then
a fourth-order one. It is the "future design document" promised in
{ref}`Pressure Gradient <omega-design-pressure-grad>` (§2.3, §3.2), which commits Omega to a
high-order option "similar to
[Adcroft et al. 2008](https://doi.org/10.1016/j.ocemod.2008.02.001)" but defers the details. This
design fills in the second `PressureGradType` entry, which the `PressureGrad` class anticipated as a
placeholder. The existing centered, second-order Montgomery-potential scheme
(`PressureGradCentered`) is kept as a runtime option and as an independent cross-check; once Phase 1
had been measured in a realistic global configuration, the new scheme replaced it as the default
(§4.4).

That entry, and the functor behind it, are named **`FiniteVolume`** and
**`PressureGradFiniteVolume`** — not anything containing "high order" (§4.1.1). The choice is
deliberate. Phase 1 uses the same two-cell horizontal stencil as the centered scheme and is second
order, not higher order, so naming the class for an order it does not have would misdescribe what
ships first; and what Phase 1 does deliver is *consistency* (§1.1), which is a separate property from
order of accuracy. `FiniteVolume` refers to the layer-integrated control-volume form of §3.1, which
both phases share.

The PGF is the single most error-sensitive term in a layered ocean model. In the
non-Boussinesq, hydrostatic momentum equation
({ref}`layered momentum equation <omega-v1-momentum-eq>`) the horizontal acceleration
from pressure and gravity is

$$
-\left(\alpha \nabla p + \nabla \Phi \right)_{e,k},
$$

in which $\alpha\nabla p$ and $\nabla\Phi = g\nabla z + \nabla(\phi_{TP}+\phi_{SAL})$
are individually large and nearly cancel. The residual is the dynamically relevant
baroclinic signal. When layers are tilted in geometric height — as they always are in
a general ALE coordinate, and severely so for thin, steeply sloped layers near shelf
breaks and (future) ice-shelf cavities — truncation error in this cancellation appears
directly as spurious velocity (the classic "pressure gradient error"). Experience with
MPAS-Ocean and published results for other models show that an accurate PGF goes a long
way toward making thin, steeply sloped layers usable.

One property of Omega's vertical coordinate simplifies the *pressure* part of this problem,
and §3.1.1–§3.2 turn on it. Pseudo-height is defined as $\tilde z \equiv -p/(\rho_0 g)$ with
**no offset** ({ref}`omega-design-governing-eqns-omega1` §5.1, which considers adding one and
rejects it so that "$\tilde z$ varies identically to pressure"). A surface of constant
$\tilde z$ is therefore a surface of constant $p$ in every column and at every time, so
$\nabla_{\tilde z}\,p \equiv 0$: the pressure-gradient term vanishes when the gradient is taken
along a coordinate surface, and the entire horizontal PGF is the geopotential compared at fixed
pressure.

**This is a statement about coordinate surfaces, not about layer interfaces, and the distinction is
the whole subject of this document.** Layer interfaces sit at $\tilde z^{\text{top}}_{i,k}(x,y,t)$ —
wherever the ALE coordinate puts them. They are tilted, they move in time, and layer $k$ occupies a
*different pressure range* in every column. Nothing here assumes otherwise, and a scheme that did
would be worthless: it is exactly that difference in pressure range between neighbouring columns that
produces the pressure-gradient error this design exists to remove.

The two facts coexist because $\nabla_{\tilde z}$ differentiates along a coordinate surface, while the
layer interfaces enter the layer-integrated form of §3.1 as the *limits of integration*. Their
horizontal variation is what generates that form's two metric terms, through Leibniz' rule, and those
terms are individually large. What §3.1.1 shows is that four individually nonzero terms **sum** to
zero — not that any one of them vanishes, and not that the interfaces are level. Every quantity the
resulting scheme actually computes is a tilt term: if interfaces were isobaric, the cross-edge
interface pressure differences would vanish, the scheme would collapse to a plain difference of
`GeomZInterface`, and there would be nothing left to design.

The reduction accordingly holds for any ALE layering *of pseudo-height* — p-star, sigma-p,
z-tilde-star, or any hybrid — because all of them relocate interfaces within $\tilde z$ and none of
them redefines $\tilde z$. It would fail only if Omega replaced pseudo-height with a coordinate that
is not a function of pressure alone: geometric $z$, a true $z$-star, or a genuine isopycnal
coordinate. §3.8 states the scope precisely and records what would have to come back in that event.

This design is shaped by three competing principles:

1. **Accuracy.** The PGF must be markedly more accurate than the centered scheme,
   especially in the presence of sloping layers.
2. **Cost.** Omega uses TEOS-10 (a 75-term polynomial) for specific volume $\alpha$.
   Naively evaluating $\alpha$ at every quadrature point would dominate run time. The
   design must bound TEOS-10 evaluations.
3. **General ALE coordinate.** The scheme must remain accurate and stable for thin,
   steeply sloped layers. The PGF alone cannot guarantee this, but it is a necessary
   ingredient.

A guiding choice that follows from principles 2 and 3 is that the metric of success *for
simulations* is **low absolute PGF error at affordable (coarse) resolution**, not high
asymptotic order of convergence. Omega will be strongly resource-limited in how fine a
mesh it can run, so a scheme whose error is small at coarse resolution is preferred over
one with a steeper convergence slope but larger coarse-resolution error. Order of
accuracy is configurable and a fourth-order variant is targeted, but it is a means to the
end of coarse-resolution accuracy, not the objective itself.

The measured convergence rate plays a different but equally essential role: it is the
**verification** metric. A scheme that does not converge at its designed order is
mis-implemented, and no amount of favorable absolute error excuses that. The two metrics
therefore answer two different questions and both must be satisfied — *is the code correct?*
(convergence rate at the designed order) and *is it useful at the resolutions we can run?*
(absolute error).

### 1.1 Why this work is delivered in two phases

Realistic global configurations run with `PressureGradCentered` show spurious flow in the
bottom layer over sloping bathymetry, large enough to drive numerical instability. What is
missing there is not high-order accuracy — it is *consistency*. §3.9 shows that
`PressureGradCentered` is, exactly and identically, the **first-order** conversion of a
height difference taken at fixed layer index into one taken at fixed pressure. Its error is
therefore the truncation of that conversion, and it accumulates downward through the column
(§3.7.3). That error is present at second order and does not require a fourth-order scheme
to remove.

The work is therefore split so the consistency fix can be delivered and tested on its own:

- **Phase 1 — a consistent second-order scheme.** The exact fixed-pressure geopotential
  difference: the equation-of-state expansion about a state shared across each edge, a
  mean-preserving *linear* reconstruction of $\Theta$ and $S$ in pressure, and the remainder
  by which the exact fixed-pressure height difference exceeds the centered scheme's
  first-order estimate of it (§3.5). The horizontal operator is the same two-cell stencil the
  centered scheme already uses, so the scheme remains second order in the horizontal. What
  Phase 1 buys is the robustness property of Requirement 2.3: the pressure gradient is zero to
  machine precision for any resting ocean whose profile varies linearly with pressure, no
  matter how the layers tilt.
- **Phase 2 — fourth order.** A parabolic vertical reconstruction and a wider horizontal
  stencil, plus the option of a second-order equation-of-state expansion. This raises the
  order of accuracy; it does not change the robustness property, which Phase 1 establishes and
  Phase 2 inherits by construction rather than having to re-establish (§3.6.1).

The two phases share one implementation, one set of configuration options, and one test
suite; Phase 2 widens settings that Phase 1 puts in place. Section 3 marks which parts of
the formulation belong to which phase, and §4.5 records what each phase delivers and what it
depends on.

Phase 1 is implemented in Omega
([#493](https://github.com/E3SM-Project/Omega/pull/493), with the sea-floor anchor refined in
[#635](https://github.com/E3SM-Project/Omega/pull/635)) and tested in Polaris
([#675](https://github.com/E3SM-Project/polaris/pull/675),
[#681](https://github.com/E3SM-Project/polaris/pull/681) and
[#900](https://github.com/E3SM-Project/polaris/pull/900)). The implementation details in §4 describe
that code, and §4.5.2 says where Phase 2 attaches to it.

## 2 Requirements

### 2.1 Requirement: Higher accuracy than the centered scheme at affordable resolution

The new PGF must produce substantially lower absolute error than
`PressureGradCentered` at the coarse-to-moderate resolutions Omega can afford to run,
for representative stratified columns with horizontal gradients of temperature, salinity,
surface pressure, and coordinate slope. Absolute error at affordable resolution is what
determines whether the scheme is *useful*; it is complemented by the separate verification
requirement in §2.6, which determines whether it is *correct*.

**"Substantially" is a factor of 5–10, and that figure comes from measurement rather than from
taste.** On a realistic global resting state the centered scheme's bottom-layer truncation error is
$1.3$–$1.6\times$ the pressure gradient the model actually carries there, so a modest gain of order
$2\times$ would leave the error comparable to the signal. §5.3.2 gives the measurement and §6.6 the
per-variant ratios achieved.

### 2.2 Requirement: Bounded TEOS-10 cost

The number of TEOS-10 specific-volume evaluations performed by the PGF must be
independent of the reconstruction order and the horizontal stencil width, and bounded at
approximately one evaluation (plus its first derivatives) per cell per layer per time step —
comparable to the cost Omega already pays to compute the `Eos::SpecVol` field. The scheme
must not require evaluating the full TEOS-10 polynomial at sub-layer sample points, nor
integrating the TEOS-10 polynomial itself semi-analytically.

### 2.3 Requirement: Robustness for thin, steeply sloped layers via discrete hydrostatic consistency

The scheme must suppress spurious pressure-gradient accelerations when layers are thin and steeply
sloped relative to surfaces of constant pressure. Four things are required.

**2.3.1 No spurious PGF for water columns the scheme resolves exactly.** For any resting ocean
whose temperature and salinity profile the scheme's vertical reconstruction reproduces exactly, the
computed horizontal PGF must be **zero to machine precision** at every edge and layer — for any
tilt of the coordinate surfaces, any layer thickness, and any bathymetry. This is the core
robustness property: it says the scheme adds nothing of its own when there is nothing to add.

**2.3.2 The design must say which water columns those are.** "Water columns the scheme resolves
exactly" is not a fixed set; it depends on the reconstruction, and it widens as the reconstruction
order increases. Each phase of this design must state its own set explicitly (§3.7.3) so that the
guarantee in §2.3.1 can be tested rather than assumed.

**2.3.3 That set must be large enough to matter.** A scheme that is exact only for an ocean of
uniform temperature and salinity would satisfy §2.3.1 trivially and be of no practical use, since
no such ocean exists. The set must at minimum include profiles that vary **linearly with
pressure**, which captures both compressibility and smooth stratification to leading order.

**2.3.4 Everywhere else, substantially smaller error than the centered scheme.** For realistic
profiles, which no reconstruction reproduces exactly, some residual PGF is unavoidable and is
*correct* — the two neighboring columns genuinely do hold slightly different water (§3.7.1). The
requirement is that this residual be substantially smaller than `PressureGradCentered`'s at the
resolutions Omega can afford, and that it not be swamped by errors the scheme itself introduces
through an inconsistent geopotential or a truncated fixed-pressure comparison.

The comparison is made on **absolute error at the resolutions Omega can afford**, not on relative
convergence order. The two schemes converge at similar rates on smooth profiles, so an
order-of-convergence gate against the centered scheme would not express what this requirement is
for; §6.5 gives the measurements behind that choice and §5.1 the gate that implements it.

These properties must hold without reference to any background or reference profile, so that they
hold equally in the warm pool, at the poles, and over a seamount.

### 2.4 Requirement: Consistency with the rest of the model state

The reconstructions used by the PGF must be **mean-preserving** (reproduce each layer's
prognostic mean $\Theta_{i,k}$, $S_{i,k}$ exactly) and **anchored to `VertCoord`**
(use the model's diagnostic interface pressures and geometric heights). The PGF must
therefore see the same layer-mean state the tracer and thickness equations evolve, and the
same hydrostatic pressure the rest of the model uses. It may assume a smoother sub-layer
profile for integration, but it must not introduce a second, inconsistent representation
of $\Theta$, $S$, $p$, or $z$.

### 2.5 Requirement: Runtime-selectable, backward-compatible option

The new scheme must be selectable at runtime through the existing `PressureGrad`
configuration group and `PressureGradType` enum, and `Centered` must remain selectable through the
same key. Its sub-options (edge stencil width, vertical-reconstruction mode, quadrature) must be
configurable and optional, so that a configuration that omits them still runs. The two phases must
share one set of configuration keys, so that a configuration written for Phase 1 continues to work
unchanged when Phase 2 lands.

Which scheme is the default is not part of this requirement. `Centered` was the default when this
design was written; Phase 1 changed that once it had been measured (§4.4).

### 2.6 Requirement: Verified order of accuracy

The implemented scheme must converge, under refinement of a smooth manufactured or
quasi-analytic reference solution, at the order of accuracy it is configured for
(nominally fourth order for the targeted variant, second order when reduced to the
centered limit). This is a verification requirement: a measured slope that falls short of
the designed order indicates a defect in the implementation — a mis-weighted quadrature
point, an inconsistent reconstruction stencil, a dropped correction term — and must be
diagnosed rather than accepted, even when the absolute error at coarse resolution already
satisfies Requirement 2.1. The two requirements are independent and both are gating:
Requirement 2.1 establishes that the scheme is useful at the resolutions Omega can afford,
Requirement 2.6 establishes that it is the scheme the design describes.

### 2.7 Desired: Extensibility to sixth order and to tidal/geoid geopotential terms

The framework should accommodate a future higher-order variant (e.g. sixth order, added to
the `PressureGradType` enum when implemented) and the
tidal-potential and self-attraction-and-loading contributions to the geopotential
(supplied by `VertCoord`; see {ref}`Pressure Gradient <omega-design-pressure-grad>`
§2.4), without restructuring.

## 3 Algorithmic Formulation

### 3.1 Continuous target form

The scheme starts from the layer-integrated, finite-volume form of the pressure
and geopotential terms derived in
{ref}`the governing equations <omega-design-governing-eqns-omega1>`
(the {ref}`full volume-integral form <omega-v1-vh-momentum-reynolds2>`). Neglecting the
turbulent correlations and dropping the resolved-component notation, the layer tendency is
({ref}`Pressure Gradient <omega-design-pressure-grad>` §3.2):

$$
T^p &= - \int_A \int_{\tilde{z}_k^{\text{bot}}}^{\tilde{z}_k^{\text{top}}} \rho_0 \, \left( \nabla \Phi \right) \, d\tilde{z} \, dA \\
& - \int_{\partial A} \left( \int_{\tilde{z}_k^{\text{bot}}}^{\tilde{z}_k^{\text{top}}} \rho_0 \left(\alpha p \right) \, d\tilde{z} \right) dl \\
& + \int_A \rho_0 \left[ \alpha \, p \, \nabla \tilde{z}_k^{\text{top}} \right]_{\tilde{z} = \tilde{z}_k^{\text{top}}} \, dA \\
& - \int_A \rho_0 \left[ \alpha \, p \, \nabla \tilde{z}_k^{\text{bot}} \right]_{\tilde{z} = \tilde{z}_k^{\text{bot}}} \, dA \\
& + \int_A \int_{\tilde{z}_k^{\text{bot}}}^{\tilde{z}_k^{\text{top}}} \rho_0 \, p \, \nabla \alpha \, d\tilde{z} \, dA.
$$ (ho-target)

The five terms are: the geopotential (gravity) body force integrated over the layer volume;
the side-wall integral of $\alpha p$ (the pressure traction on the cell faces); two metric
terms accounting for the pressure traction on the sloping top and bottom layer interfaces;
and a $p\nabla\alpha$ term.

**Two of those terms differ from [](#ho-target) as inherited from
{ref}`omega-design-governing-eqns-omega1`, and the differences are corrections.** The fifth term is
**absent** there; omitting it leaves a scheme that does not converge at all. The signs on lines 3 and 4
are **opposite** to those printed there; with the inherited signs the cancellation of §3.1.1 doubles
rather than vanishing. Both follow from applying Leibniz' rule directly, immediately below, and neither
is a convention difference in how the traction normal or $\nabla\tilde z^{\text{top}}$ is defined.
`OmegaV1GoverningEqns.md` is not maintained against the code and should be corrected there rather than
overridden here; §6.2 states what is owed to it and §6.1 the supporting measurements.

The derivation is short enough to give in full, because everything in §3 rests on it.
Differentiating $\tilde z \equiv -p/(\rho_0 g)$ and applying the hydrostatic relation gives the mass
element $\rho\,dz = \rho_0\,d\tilde z$ ({ref}`omega-design-governing-eqns-omega1` §5). The layer
integral of the pointwise acceleration $-(\alpha\nabla p + \nabla\Phi)$ at fixed $\tilde z$ is
therefore
$-\int\rho_0(\alpha\nabla p + \nabla\Phi)\,d\tilde z$. Writing
$\alpha\nabla p = \nabla(\alpha p) - p\nabla\alpha$ and applying Leibniz' rule to the first
piece,

$$
\int_{\tilde z_k^{\text{bot}}}^{\tilde z_k^{\text{top}}} \nabla(\alpha p)\,d\tilde z
= \nabla \int_{\tilde z_k^{\text{bot}}}^{\tilde z_k^{\text{top}}} \alpha p \,d\tilde z
- (\alpha p)_{\text{top}}\,\nabla \tilde z_k^{\text{top}}
+ (\alpha p)_{\text{bot}}\,\nabla \tilde z_k^{\text{bot}},
$$ (leibniz)

gives exactly lines 2–5 of [](#ho-target), with the divergence theorem turning the first
term on the right into the side-wall line integral. This is the Adcroft et al. (2008)
finite-volume route — the net pressure force on each control volume obtained by integrating
in-situ pressure over its faces, rather than by forming the pointwise product
$\alpha\nabla p$ — and it is the form MOM6 uses successfully for the same problem.

#### 3.1.1 In pseudo-height the four pressure terms cancel identically

Lines 2–5 of [](#ho-target) are, by construction, nothing but
$-\int_A\int\rho_0\,\alpha\nabla_{\tilde z} p\,d\tilde z\,dA$ rewritten. In Omega's
coordinate that quantity is **identically zero**.

Pseudo-height is $\tilde z \equiv -p/(\rho_0 g)$ with no offset, per
{ref}`omega-design-governing-eqns-omega1` §5.1, which considers adding one and rejects it so
that "$\tilde z$ varies identically to pressure". The map $p \leftrightarrow \tilde z$ is
therefore *universal*: it is the same in every column, at every time, irrespective of surface
pressure. A surface of constant $\tilde z$ is a surface of constant $p$, so

$$
\nabla_{\tilde z}\, p \;\equiv\; 0,
$$ (grad-p-zero)

and lines 2–5 of [](#ho-target) sum to zero — not approximately, and not only for a resting
ocean, but as an algebraic identity for any state whatsoever. This is the standard result that
in an isobaric coordinate the horizontal pressure gradient force is $-\nabla_p\Phi$.

**What [](#grad-p-zero) does and does not say.** It says the coordinate *surfaces* are isobars.
It says nothing about the layer *interfaces*, which are at $\tilde z^{\text{top}}_{i,k}(x,y,t)$
and are tilted, moving, and at different pressures in different columns — as they must be in an
ALE coordinate. The horizontal variation of those limits of integration is precisely what
produces the two metric terms through [](#leibniz), and each of the four pressure terms is
individually large: measured on a two-column state at 50 m/km tilt, the side-wall term alone is
$0.55$ m s$^{-2}$. The cancellation is a property of their *sum*. Read the other way, it is a
statement that the four terms carry no information the geopotential term does not already carry,
so discretizing them can only add error — not a statement that the tilt they describe is absent.

Three things follow, and they set the shape of the rest of §3.

1. **The entire horizontal PGF is the geopotential compared at fixed pressure.** §3.2 states
   the resulting target. Comparing "at fixed pressure" is hard exactly because the interfaces
   are not at fixed pressure; the difficulty of the original problem is preserved, relocated
   into a single term.
2. **Adcroft's problem is not Omega's problem.** Adcroft et al. (2008) needed the face-integral
   form because their coordinate was *not* isobaric, so $\alpha\nabla_r p$ was a real term
   requiring careful treatment. Omega's coordinate is isobaric by definition, and remains so
   under ALE: p-star, sigma-p and any other layering Omega adopts move interfaces *within*
   $\tilde z$ without redefining $\tilde z$, so [](#grad-p-zero) is untouched. The remaining
   problem is not a force-balance problem at all — it is the problem of integrating the
   hydrostatic relation accurately enough that two columns' geopotentials can be differenced at
   a common pressure.
3. **Discretizing the four terms separately is worse than not discretizing them.** Because they
   must cancel, any discretization of them contributes only error — and the cancellation is
   delicate in a way that is not apparent from [](#ho-target). It requires the $p\nabla\alpha$
   term to be integrated over the mean of the *two columns' own* layer pressure ranges rather
   than over the edge-mean range, which is the choice an edge-centred implementation invites;
   mispaired, the residual is $O(\tilde h^0)$ and does not converge under vertical refinement at
   all, while remaining invisible to every gate in §5 — on the exact set the mispaired term is
   identically zero. The reduced target of §3.2 has no such requirement, because it has no terms
   to pair.

[](#grad-p-zero) has been confirmed to hold discretely as well as continuously: the four pressure
terms were measured to sum to $3\times10^{-15}$ of the size of the side-wall term alone, at 50 m/km
tilt, for every profile tried. §6.1 gives the numbers and the harness.

### 3.2 The exact target

Pressure is an exactly linear function of $\tilde z$,

$$
p(\tilde z) = p_{k}^{\text{top}} - \rho_0 g \left( \tilde z - \tilde z_k^{\text{top}} \right),
\qquad \tilde z_k^{\text{bot}} \le \tilde z \le \tilde z_k^{\text{top}},
$$ (p-linear)

anchored at the interface pressures $p_k^{\text{top}}$, $p_k^{\text{bot}}$ that `VertCoord`
already provides (`PressureInterface`). No reconstruction of $p$ within a layer is needed or
permitted — using [](#p-linear) is exactly the `VertCoord` anchoring required by
Requirement 2.4. Within a column, pressure is used from here on as the vertical integration
variable in place of $\tilde z$; the two are interchangeable by [](#p-linear), and pressure is
the variable in which two columns' profiles can be compared.

Combining [](#grad-p-zero) with [](#ho-target) and mass-weighting the layer average with
$\rho_0\,d\tilde z = -dp/g$, the exact edge-normal acceleration at edge $e$, layer $k$ is

$$
T^p_{e,k} = -\,\frac{g}{\Delta p_{e,k}} \int_{p_{e,k}^{\text{top}}}^{p_{e,k}^{\text{bot}}}
\left[\nabla_n z\right]_{p} \; dp
\;-\; \left[\nabla_n \left(\phi_{TP} + \phi_{SAL}\right)\right]_{e},
$$ (ho-exact)

in which $[\nabla_n z]_p$ is the edge-normal gradient of geometric height taken **at fixed
pressure**, $\Delta p_{e,k} = p^{\text{bot}}_{e,k} - p^{\text{top}}_{e,k}$ is the edge control
volume's pressure range, and the geometric height of column $i$ as a function of pressure is
the hydrostatic integral

$$
z_{i}(p) = z_{i}^{\text{anchor}} + \frac{1}{g}\int_{p}^{p_i^{\text{anchor}}}
\hat\alpha_{i}(p')\, dp',
$$ (z-of-p)

with $\hat\alpha_{i}(p)$ the reconstructed specific volume in column $i$ (§3.3–§3.4).
Equation [](#ho-exact) is the whole scheme. Everything that follows is either an ingredient of
$\hat\alpha$ (§3.3, §3.4), a statement about how [](#ho-exact) is evaluated without losing
precision (§3.5), or a statement about when it returns zero (§3.7).

Two properties of [](#ho-exact) are worth noticing immediately, because between them they do the work
that a four-term force balance would need several delicate arguments to secure.

- **The robustness property is structural, not a cancellation.** If the two columns' reconstructed
  profiles agree as functions of pressure, $\hat\alpha_L(p) \equiv \hat\alpha_R(p)$, then by
  [](#z-of-p) the two $z_i(p)$ differ by a constant fixed by the anchors, and $[\nabla_n z]_p$ is
  *pointwise* zero at every pressure. The integral in [](#ho-exact) is then zero regardless of
  the tilt, the layer thicknesses, the bathymetry, or the quadrature used to evaluate it. There
  are no terms to pair and no signs to get right. For a resting ocean the anchors agree as well,
  the sea surface tilting exactly so as to compensate a horizontal gradient in $p^{\text{surf}}$.
- **The vertical quadrature is not an accuracy knob.** $[\nabla_n z]_p$ varies across the layer
  only through $\frac{d}{dp}[\nabla_n z]_p = -\frac{1}{g}[\nabla_n\hat\alpha]_p$, the horizontal
  contrast in specific volume at fixed pressure — a quantity that is *small* and that vanishes
  identically on the exact set. The layer average in [](#ho-exact) is therefore a nearly-constant
  integrand, handled in closed form for a polynomial $\hat\alpha$. There is no quadrature setting
  (§4.1.1).

### 3.3 Reference-state expansion of the equation of state

This is the core device that satisfies Requirement 2.2. We **never integrate TEOS-10**.
Instead, for each cell $i$ and layer $k$ we expand $\alpha$ in a Taylor series about the
layer's reference state $(\Theta_{i,k}, S_{i,k}, p_{i,k}^{\text{mid}})$, where
$\Theta_{i,k}$ and $S_{i,k}$ are the prognostic layer means and $p_{i,k}^{\text{mid}}$ is
the mid-layer pressure from `VertCoord` (`PressureMid`):

$$
\alpha(\Theta, S, p) \approx \alpha_{0} + \alpha_{\Theta}\,(\Theta - \Theta_{i,k}) + \alpha_{S}\,(S - S_{i,k}) + \alpha_{p}\,(p - p_{i,k}^{\text{mid}}),
$$ (alpha-taylor)

with

$$
\alpha_0 = \alpha(\Theta_{i,k}, S_{i,k}, p_{i,k}^{\text{mid}}), \quad
\alpha_\Theta = \left.\frac{\partial \alpha}{\partial \Theta}\right|_0, \quad
\alpha_S = \left.\frac{\partial \alpha}{\partial S}\right|_0, \quad
\alpha_p = \left.\frac{\partial \alpha}{\partial p}\right|_0.
$$ (alpha-derivs)

The four coefficients in [](#alpha-derivs) are obtained from a **single** evaluation of the
TEOS-10 polynomial: the derivatives are analytic derivatives of that same polynomial, taken
at the same normalized state, so no second call to the equation of state occurs (they are the
same quantities computed by `gsw.specvol_first_derivatives`, which the Polaris reference
solution already uses). The arithmetic is not free — $\alpha_p$ reuses the pressure
coefficients assembled for $\alpha_0$, but $\alpha_\Theta$ and $\alpha_S$ need coefficient
sets of their own — though it shares the square root and the normalization, and it is the
*call count*, not the polynomial arithmetic, that Requirement 2.2 bounds. Indeed $\alpha_0$ is
exactly the `Eos::SpecVol` field Omega already computes, so only the three derivative fields
are additional work.

The first-order expansion [](#alpha-taylor) is what Phase 1 implements. A second-order expansion
(adding $\alpha_{\Theta\Theta}$, $\alpha_{SS}$, $\alpha_{\Theta S}$, $\alpha_{pp}$, …) is a Phase 2
option where stratification is strong; it adds derivative evaluations but still no in-integral
TEOS-10 calls. It has no configuration key yet.

#### 3.3.1 The expansion point must be shared across the edge

Expanding about the *cell-local* reference state $(\Theta_{i,k}, S_{i,k}, p_{i,k}^{\text{mid}})$
would quietly defeat the robustness property of §3.7. Two neighboring columns would then be using
two slightly different approximations to the *same* equation of state, taken at two different
points. Their specific-volume profiles $\hat\alpha_{L,k}(p)$ and $\hat\alpha_{R,k}(p)$ would then
disagree — by $O(\alpha_{pp}\,\Delta p^{\text{mid}}_e\, \tilde h)$ — even for a resting ocean the
reconstruction reproduces perfectly. By [](#z-of-p) the two columns' $z_i(p)$ would then differ by
more than a constant, $[\nabla_n z]_p$ would be nonzero, and the scheme would generate spurious
flow from nothing but its own EOS approximation. This choice remains load-bearing under the reduced
target of §3.2 — it is the *only* remaining way for the exactness of [](#ho-exact) to fail on a
resolved profile — and guard test (b) of §5.2 confirms it fires.

The design therefore **shares the expansion point across each edge**. The coefficients
[](#alpha-derivs) are computed once per cell per layer, exactly as described above (one TEOS-10
evaluation, Requirement 2.2), and are then averaged to the edge together with the reference state
itself:

$$
\bar\alpha_{0}^{e,k}, \; \bar\alpha_{\Theta}^{e,k}, \; \bar\alpha_{S}^{e,k}, \; \bar\alpha_{p}^{e,k},
\qquad
\left(\bar\Theta^{e,k}, \bar S^{e,k}, \bar p^{e,k}\right) = \tfrac12 \sum_{i \in CE(e)} \left(\Theta_{i,k}, S_{i,k}, p_{i,k}^{\text{mid}}\right),
$$ (edge-ref)

and **both** columns' contributions to edge $e$ are evaluated with this single set, so that the two
sides of the edge see one and the same equation of state.

This has a cost consequence worth recording. The column quantity built from $\hat\alpha$ — under §3.2
the hydrostatic integral [](#z-of-p) — is not a per-cell quantity that can be computed once and
differenced across each of the cell's edges. It depends on the edge through [](#edge-ref) and must be
evaluated per edge, twice. On a hexagonal TRiSK mesh this is roughly three times as many column evaluations as a
cell-based formulation. The TEOS-10 call count is unaffected, which is the cost Requirement 2.2
binds; the extra work is polynomial arithmetic on coefficients already in cache. §3.5 shows how to
arrange it so that only a *small* correction is accumulated per edge, which bounds both the cost
and the round-off.

Because $\Theta$, $S$ are reconstructed as low-order polynomials in pressure and the expansion
[](#alpha-taylor) is linear in $\Theta$, $S$ and $p$, $\hat\alpha$ is a **low-order polynomial in
pressure** and every integral of it required by §3.2 and §3.5 is available in closed form. No
TEOS-10 evaluation occurs inside any integral — which is what resolves the feasibility and cost
concerns with applying Adcroft-style analytic integration to TEOS-10.

### 3.4 Mean-preserving vertical reconstruction

Within layer $k$ of column $i$ we reconstruct the conservative temperature and absolute
salinity as deviations from the prognostic layer means,

$$
\Theta(p) = \Theta_{i,k} + \Theta'_{i,k}(p), \qquad
S(p) = S_{i,k} + S'_{i,k}(p),
$$ (vert-recon)

where the deviations are polynomials in pressure built from the neighboring layer means and
constrained to **integrate to zero over the layer**
($\int_{p_{i,k}^{\text{top}}}^{p_{i,k}^{\text{bot}}} \Theta'_{i,k}\, dp = 0$). Because pressure
and $\tilde z$ differ within a column only by a fixed linear map [](#p-linear), a polynomial in
one is a polynomial of the same degree in the other; pressure is used here because it is the
variable in which the profiles of two neighboring columns can be compared (§3.7).

The mean-preserving property (Requirement 2.4) guarantees the PGF uses the same layer-mean state
as the rest of the model; the reconstruction only supplies the sub-layer shape needed to integrate
$\alpha p$.

**Phase 1 uses linear deviations; Phase 2 uses parabolic (PPM-style) ones.** The degree sets the
scheme's exact set in §3.7.3 directly: linear deviations make the scheme exact for profiles that
vary linearly with pressure, parabolic ones for profiles that vary quadratically.

That correspondence only holds if the reconstruction actually reproduces those profiles, which
places two requirements on how the deviations are built:

- **The estimator must be exact on a non-uniform vertical grid.** Omega's layers are not of equal
  thickness, and a slope or curvature formula derived assuming equal thickness will not recover a
  linear profile exactly when they are not. The estimator must be built from the layer means and
  the actual interface pressures, and must return the exact slope (Phase 1) or slope and curvature
  (Phase 2) when the underlying profile is of that degree, for any distribution of layer
  thicknesses. This is testable directly and is part of the unit test in §5.2.
- **Any limiter must be inactive on such profiles.** Monotonicity limiters are not required by the
  PGF and are not applied by default here, since the reconstruction feeds an integral rather than
  an advective flux. If one is later added for robustness, it must leave smooth monotone data
  untouched, or it will break the cancellation of §3.7.2 precisely where the profile is well
  resolved.

**The Phase 1 estimator.** Layer $k$'s deviation is $\Theta'_{i,k}(p) = s_{i,k}\,(p -
p^{\text{mid}}_{i,k})$, with the slope $s_{i,k}$ the difference of the two neighbouring layer means
divided by the difference of their `PressureMid`. The difference is centred, over layers $k-1$ and
$k+1$, in the interior; it is one-sided in the shallowest and deepest valid layers, and the slope is
zero in a column with a single valid layer. Because the mean of a linear function over a layer is
its value at the layer's exact pressure midpoint, the layer means of a profile linear in pressure
lie exactly on that line as a function of `PressureMid`, so this recovers the slope exactly for
any distribution of layer thicknesses. The one-sided branches matter: the deepest valid layer is
the partial cell that carries the whole signal where bathymetry steps. No limiter is applied. The
slopes are formed once per cell and layer and reused across each cell's edges.

The isopycnal-model assumption of constant properties within a layer is the degenerate choice
$\Theta'_{i,k} = S'_{i,k} = 0$. This is retained as a **verification-only** configuration — it is
what makes the reduction to the centered scheme in §3.9 possible, and it isolates compressibility
on its own in §5.2 — but it is not a supported production setting in either phase, because it
gives up the exactness for linearly varying profiles that Requirement 2.3.3 asks for.

### 3.5 Evaluating the fixed-pressure difference

Equation [](#ho-exact) is not evaluated by building each column's $z_i(p)$ and subtracting. It is
evaluated by **differencing the integrand first and integrating the difference**, which is what
makes the exactness of §3.7 structural rather than the outcome of a cancellation.

Write $\Delta_e f \equiv \sum_{i\in CE(e)} -n_{e,i}\,f_i$, so that $[\nabla_n f]_e = \Delta_e f/d_e$,
and let $\Delta_e z(p)$ be the fixed-pressure height difference [](#ho-exact) integrates.
Differentiating [](#z-of-p) gives the relation the whole scheme rests on:

$$
\frac{d}{dp}\,\Delta_e z(p) \;=\; -\,\frac{1}{g}\,\Delta_e \hat\alpha(p),
\qquad
\Delta_e \hat\alpha(p) \equiv \hat\alpha^{(e)}_R(p) - \hat\alpha^{(e)}_L(p),
$$ (dz-dp)

in which both columns are evaluated at **the same pressure** $p$. Two rules fix what that means,
and between them they are the whole of the robustness property:

- **One shared coefficient set at each pressure.** For $p$ in edge layer $k$, both columns use the
  edge-shared coefficients and expansion state [](#edge-ref) of *that edge layer* — not of each
  column's own layer $k$, which under tilt spans a different pressure range (§3.7.2, condition 1).
- **Each column supplies its own $\Theta(p)$, $S(p)$**, from whichever of *its own* layers contains
  $p$, using that layer's mean-preserving reconstruction (§3.4).

With those, the $\bar\alpha_0$ and $\bar\alpha_p\,(p-\bar p^{e})$ terms of [](#alpha-taylor) are
identical in the two columns and cancel in the difference, leaving

$$
\Delta_e \hat\alpha(p) \;=\; \bar\alpha_{\Theta}^{e,k}\,\Delta_e\Theta(p)
\;+\; \bar\alpha_{S}^{e,k}\,\Delta_e S(p).
$$ (dalpha)

**This is the central property of the scheme.** [](#dalpha) is a product of a coefficient with the
horizontal contrast in reconstructed $\Theta$ and $S$ *at matched pressure*, and that contrast is
**identically zero, pointwise**, whenever the two columns' reconstructions describe the same water
(§3.7.2, condition 1). Four consequences follow:

1. **Exactness does not depend on the coefficients' values, but does depend on their being
   shared.** Whatever $\bar\alpha_\Theta$, $\bar\alpha_S$ are, they multiply zero — so *which*
   shared set is used is an accuracy question, and selecting it by edge layer rather than by some
   other rule cannot break the robustness property. What is load-bearing is that **one** set
   multiply both columns: [](#dalpha) collapses to a coefficient times a contrast only then. Give
   each column its own expansion point and the $\bar\alpha_0$ and $\bar\alpha_p$ terms no longer
   cancel, and exactness is lost — measured at $10^{-5}$, which is guard (b) of §5.2. The
   edge-shared expansion point of §3.3.1 is therefore load-bearing, and it is load-bearing for
   *sharing* rather than for where the point sits.
2. **Exactness does not depend on the quadrature.** The integrand is zero at every point, so any
   rule integrates it to zero. Quadrature order is therefore an ordinary accuracy knob (§4.1.1),
   not a correctness requirement.
3. **Exactness does not depend on the interfaces lining up.** Where the two columns' layers are
   offset — at 50 m/km and 64 m layers they are offset by nearly three layer thicknesses, and
   overlap not at all — [](#dalpha) is still evaluated at matched pressure and is still zero.
4. **Compressibility drops out of the horizontal gradient entirely.** $\bar\alpha_p$ does not appear
   in [](#dalpha). This is correct physics rather than an approximation: if $\Theta$ and $S$ are
   horizontally uniform then so is $\alpha(p)$, and a horizontally uniform compressibility exerts no
   horizontal pressure gradient. It is *only* because `PressureGradCentered` compares at fixed layer
   index that compressibility appears in its error at all (§3.9).

#### 3.5.1 The column scan

$\Delta_e z$ is accumulated along each edge's column from an anchor at one end. §3.7.4 settles which
end, and it is the **sea floor**: the recurrence below therefore runs upward, from $k = K+1$ to
$k = 1$. Let $\bar q_k$ be the edge-layer interface pressures and $D_k \equiv \Delta_e z(\bar q_k)$.
Then

$$
D_{k+1} \;=\; D_k \;-\; \frac{1}{g}\int_{\bar q_k}^{\bar q_{k+1}} \Delta_e\hat\alpha(p)\; dp,
$$ (d-recurrence)

and the layer mean [](#ho-exact) needs is obtained from $D_{k+1}$ and the first moment of the same
integrand about the layer's top interface, $\int_{\bar q_k}^{\bar q_{k+1}} (p - \bar q_k)\,
\Delta_e\hat\alpha\,dp$. Both integrals use the same quadrature points, and at each point each
column's $\Theta$, $S$ come from the layer of that column containing the point.

**The anchor.** $K$ is the edge's deepest layer, the deepest layer valid in both columns, so
$\bar q_{K+1}$ lies at the shallower column's floor and, where `maxLevelCell` differs, at an interior
interface of the deeper one. The two columns' interfaces $K+1$ are in general at different
pressures, and $D_{K+1} = \Delta_e z(\bar q_{K+1})$ is the height difference there corrected to a
common pressure,

$$
D_{K+1} \;=\; \Delta_e Z_{K+1} \;+\; \frac{1}{g}\sum_{i \in CE(e)} -n_{e,i} \int_{\bar q_{K+1}}^{q_{i,K+1}} \hat\alpha^{(e)}_i(p)\,dp ,
$$ (anchor)

with both short integrals evaluated exactly by the scan's Gauss rule, each spanning half the
cross-edge pressure difference at that interface, and both vanishing where the two columns'
interface pressures agree. Anchored at the surface instead, the same expression applies at $k = 1$
with the sea-surface height difference and the surface pressures.

**The shift to a common pressure must stay within one layer.** The short integrals use one
first-order expansion, shared by the two layers that meet at the anchor (§3.3.1), and it is accurate
across a pressure difference of about a layer. That covers a flat floor, and most of a sloping one
cut by partial bottom cells beneath level interior layers: there the two columns' interfaces $K+1$
are within a layer of each other, except where a thin bottom cell lies entirely above the deeper
column's layer $K$. Terrain-following layers over a steep slope put one floor
several layers deeper than the other, and [](#anchor) then stretches the expansion across the whole
offset. The resulting error grows with the cube of the offset and is the same at every depth.
Measured in Omega on a resting state linear in pressure with the floors 6.4 layers apart, it leaves
the scheme at 0.45 of `PressureGradCentered`'s error (§5.2).

The comparison is therefore made where the two columns are within a layer of each other. Let column
$s$ be the one whose interface $K+1$ is at the lower pressure and column $d$ the other, and let $j
\le K+1$ be the interface of column $d$ at or just below $q_{s,K+1}$, the lower interface of the
layer of $d$ that contains $q_{s,K+1}$. Each column contributes one interface, $(Z^\ast_i,
q^\ast_i) = (Z_{s,K+1}, q_{s,K+1})$ or $(Z_{d,j}, q_{d,j})$, and with $q^\ast = \tfrac12(q^\ast_s +
q^\ast_d)$

$$
D(q^\ast) \;=\; \Delta_e Z^\ast \;+\; \frac{1}{g}\sum_{i \in CE(e)} -n_{e,i} \int_{q^\ast}^{q^\ast_i} \hat\alpha^{(a)}_i(p)\,dp,
\qquad
D_{K+1} \;=\; D(q^\ast) \;-\; \frac{1}{g}\int_{q^\ast}^{\bar q_{K+1}} \Delta_e\hat\alpha(p)\,dp ,
$$ (anchor-one-layer)

in which $\hat\alpha^{(a)}$ uses one expansion averaged from the two layers that meet at the
comparison, layer $K$ of column $s$ and layer $j-1$ of column $d$, and the second integral is the
matched-pressure integrand [](#dalpha) with edge layer $K$'s shared expansion. That integral is an
increment of the same kind as every step of [](#d-recurrence): it vanishes on the exact set, so it
carries the anchor down to $\bar q_{K+1}$ without disturbing exactness. Where $j = K+1$,
$q^\ast = \bar q_{K+1}$, the second integral has zero width, and [](#anchor-one-layer) is
[](#anchor) bit for bit. The rule changes answers only where the floors are more than a layer
apart, and on the 6.4-layer state above it brings the scheme from 0.45 to 0.014 of
`PressureGradCentered`.

What it cannot remove is the difference between the two columns' `VertCoord` heights. $Z_{s,K+1}$
and $Z_{d,j}$ are each `VertCoord`'s midpoint-rule sum over its own column's layers (§3.7.4), and
two different layer partitions give sums that differ at $O(\tilde h^2)$. That residual converges
with layer thickness rather than growing with the slope. It is why Requirement 2.3.1 is met to
machine precision only where the two columns' floors are at the same depth and the same pressure.

$D_{K+1}$ is **computed, not assumed**: it is whatever the model's geometric heights and interface
pressures imply, evaluated at a common pressure. Like every other quantity in the scan it is a
fixed-pressure height difference, and it is zero exactly when the two columns' heights agree at a
common pressure there — which is what a state at rest means. A state that is only approximately at
rest carries a real gradient there, and the scheme reports it. The anchor is therefore not a place
where the scheme assumes anything about the state; it is the deepest instance of the same
comparison [](#d-recurrence) makes at every other interface.

**Which end of the column the anchor sits at is not a conditioning choice**, though it was described
as one here before it was measured; §3.7.4 gives the argument and the trade-off. What follows is the
conditioning half of it.
`VertCoord::computeGeomZHeight` sets $Z_{i,\,\text{bot}} = -\text{BottomGeomDepth}_i$ and accumulates
$\Delta z = \rho_0\,\alpha_{i,k}\,\tilde h_{i,k}$ **upward**: the bathymetry is prescribed and the
sea-surface height is derived, so

$$
\Delta_e Z_{\text{surf}} = -\Delta_e H + \frac{1}{g}\,\Delta_e\!\left[\sum_k \alpha_{i,k}\,\Delta p_{i,k}\right].
$$ (ssh-difference)

The bracket is not an error. It is how the model defines geometric height, and the PGF must see the
model's own $z$ — a scheme that returns zero for the state the model actually holds is what
Requirement 2.3.1 asks for, so using the derived sea-surface height is self-consistent and correct.
What [](#ssh-difference) does show is that $\Delta_e Z_{\text{surf}}$ is a difference of two
column-length accumulations yielding a small result, where $-\Delta_e H$ is exact input and vanishes
identically for a flat floor. The sea-floor anchor is therefore better conditioned in floating point
*where the floor is unstepped*, at the cost of a larger common-pressure excursion, a partial cell to
handle where `maxLevelCell` differs, and — where the floor does step — a cancellation of order the
step. §3.7.4 records the trade-off and what is and is not established about it.

**Every quantity in the scan is small.** $D_k$ is the fixed-pressure height difference — of order
$10^{-1}$ m for a realistic baroclinic column and zero for a resting one, against the $10^{2}$ m
height *differences at fixed layer index* that a scheme comparing at fixed index has to form and
cancel. The increments in [](#d-recurrence) are smaller still, being integrals of a horizontal
contrast.

**The anchor is the exception, and only at the sea floor.** Where the bathymetry steps,
$\Delta_e Z$ at the anchor is of order the part of the step that falls within the deepest shared
layer, which [](#anchor-one-layer) keeps below a layer thickness, and the two short integrals of
[](#anchor) span it and cancel it. Measured on a two-column state with a stepped floor and a
profile linear in pressure: at a floor gradient of 1 m/km the anchor forms and cancels $4.0$ m to
leave $1.9\times10^{-12}$ m; at 25 m/km, $100$ m to leave $2.6\times10^{-8}$; at 200 m/km,
$205.6$ m to leave $1.7\times10^{-7}$ — nine to twelve decimal digits consumed. Anchored at the
surface the same states form nothing at all: $\Delta_e Z$ there is zero and, both columns sharing a
surface pressure, the correction is a zero-width integral. §3.7.4 records what this does to the
choice of end.

**`VertCoord`'s geometric height is used once**, for $\Delta_e Z^\ast$ in [](#anchor-one-layer): one
value of `GeomZInterface` per column, at the interface that column contributes to the anchor. It is
not accumulated and not re-derived, so Requirement 2.4 is met by the PGF not constructing a
geometric height at all: it constructs a *difference*, which is not a height field.

### 3.6 The edge operator

The discrete tendency lives at edge $e$ and is the edge-normal projection of the PGF. Every
horizontal difference in §3.5 is the two-cell TRiSK operator

$$
\left[\nabla_n f\right]_{e} = \frac{1}{d_e}\sum_{i \in CE(e)} -n_{e,i}\, f_{i}^{(e)},
$$ (edge-grad)

applied to column quantities $f_i^{(e)}$ evaluated with the equation-of-state coefficients shared
across edge $e$ [](#edge-ref). **Phase 1 uses this operator as written** — the same two-cell stencil
`PressureGradCentered` uses, so the horizontal accuracy is unchanged from today. Everything Phase 1
gains is in the vertical and in the fixed-pressure comparison.

Nothing is assumed constant across the edge: the horizontal contrast in $\Theta$ and $S$ between the
two columns enters at full strength through the two columns' $\hat\alpha^{(e)}_i(p)$ and hence
through $z_i(p)$. What Phase 1 does *not* do is raise the order of the horizontal difference, which
is Phase 2's job.

#### 3.6.1 Raising the horizontal order (Phase 2)

**The robustness property of §3.7 places no constraint on the wide stencil.** This is worth stating
explicitly, because the opposite would be true of a four-term force-balance formulation, where a
quantity interpolated to the edge from four or more cells belongs to no control volume and so has no
reason to balance. That would force the stencil to be a weighted sum of two-column pair
contributions, and settling its form on an unstructured, variable-resolution mesh would be a
prerequisite for Phase 2. Under [](#ho-exact) it is not.

The reason is that the quantity being differenced is $z_i(p)$ — a *per-column function of pressure*,
not a nonlinear product formed at the edge. For a resting ocean in the exact set, [](#z-of-p) makes
$z_i(p)$ the same function of $p$ in every cell of the neighbourhood, up to the constant set by each
column's anchor, and for a resting ocean those constants agree as well. **Any** consistent
edge-normal gradient operator — any operator that annihilates a constant, which every consistent
gradient operator does — therefore returns zero. The machine-precision property is inherited by
construction, and there is no pair-decomposition requirement.

What remains for Phase 2 is an ordinary accuracy question rather than a structural one:

- The **stencil width and weights** are chosen for order of accuracy on Omega's unstructured,
  variable-resolution TRiSK mesh, drawing on the wider edge neighbourhoods used for high-order
  tracer reconstruction (cf. {ref}`omega-design-governing-eqns-omega1` §10 and White & Adcroft
  (2008)). Assumption A1 (§3.7.6) still applies: the achieved order must be measured, not assumed.
- The **expansion point** of [](#alpha-taylor) must be shared across whatever set of cells the
  stencil couples, for the same reason §3.3.1 gives at a two-cell edge. Sharing one point across the
  whole neighbourhood is the simplest option and is sufficient for exactness; building the stencil as
  a weighted sum of pairs, each with its own shared point, keeps the expansion error smaller and is
  therefore preferable on assumption-A2 grounds. Both work; the choice is now a trade-off to be
  measured rather than a constraint imposed by the robustness property.

Phase 1 does not depend on either answer.

### 3.7 Where the cancellation is exact, and where it is not

This section makes Requirement 2.3 concrete: it explains why the obvious version of the robustness
property cannot be delivered, states the three conditions the implementation must meet, lists the
water columns that cancel exactly under each phase, and records what the design is still assuming
and must therefore test. It does *not* rely on subtracting a background profile.

#### 3.7.1 Why "uniform $\Theta$ and $S$" is the wrong condition

The natural way to state the robustness property is to require the scheme to return zero whenever
the layer means $\Theta_{i,k}$, $S_{i,k}$ are the same in both cells of an edge, for any tilt and
any vertical profile. **That requirement cannot be met, and asking for it would send the
implementation chasing an impossible target.**

When coordinate surfaces tilt relative to surfaces of constant pressure, layer $k$ in one column and
layer $k$ in its neighbor cover *different pressure ranges*. Each column builds its own
specific-volume profile from its own layer means on its own set of interfaces. Unless the vertical
reconstruction happens to reproduce the true profile exactly, the two columns therefore describe
slightly *different water*. Real water columns that differ do exert a real pressure gradient on one
another. A scheme that returned zero there would be hiding an error, not avoiding one.

The same argument runs the other way, and this is the more useful half: take a perfectly uniform
ocean whose temperature and salinity depend only on pressure, and sample it onto tilted layers. The
resulting layer means are **different in the two columns**, because the layers average over
different pressure ranges. So "uniform layer means" is neither necessary nor sufficient for the
cancellation we want. The condition that matters is a property of the reconstructed profiles, not of
the layer means, and it is stated next.

#### 3.7.2 Three conditions the implementation must meet

For each edge and layer:

1. **Both columns describe the same water, as a function of pressure, over the whole column.** The
   two columns must evaluate $\hat\alpha$ as one and the same function of $p$ — not merely within a
   layer, but at **every pressure from the anchor down**, since [](#z-of-p) integrates from the
   surface and [](#d-recurrence) accumulates from it.

   **The per-layer version of this condition is not sufficient, and the difference is the whole
   reason for [](#dalpha).** Requiring the two columns to agree only "over the pressure range the
   layer spans" is satisfiable while the scheme still misses exactness by three orders of magnitude,
   because [](#z-of-p) integrates from the surface and [](#d-recurrence) accumulates from it: two
   column integrals taken over different pressure ranges, with different layer partitions, have
   nothing forcing them to agree. Differencing the integrand at matched pressure makes the condition
   hold pointwise and by construction, and it is verifiable pointwise — [](#dalpha) is zero at a
   quadrature point or it is not. §6.5 gives the measurements.

2. **The two columns are differenced at a common pressure, before integration.** This is
   [](#dz-dp)–[](#dalpha). Integrating each column separately and subtracting satisfies condition 1
   only to the order of whatever shift is used to reconcile the two, and reintroduces the
   large-number cancellation §3.7.5 is about.
3. **The anchor is shifted to a common pressure too.** [](#anchor-one-layer), not the raw
   $\Delta_e Z$ at the anchor. Under any coordinate tilt, and wherever the bathymetry slopes, the
   two columns' interfaces at the sea floor sit at different pressures, so the comparison condition
   2 imposes in the interior must be made at the bottom of the column as well. This is the deepest
   instance of condition 2 rather than a separate idea, but it is stated separately because it is a
   distinct place in the code and fails in a distinguishable way — a depth-independent offset in
   $D_k$ (§5.2, guard (d)). The shift must also stay within one layer (§3.5.1), or the anchor's
   first-order expansion is stretched across the floor offset and leaves a depth-independent offset
   of its own.

Conditions 2 and 3 are about the discretization alone. Condition 1 is about
whether the *reconstruction* reproduces the true profile: it holds exactly for profiles in the exact
set of §3.7.3, and fails by however much the reconstruction misses when it does not. Meeting all
three gives a PGF that is **zero to machine precision**, for any tilt and any layer thickness,
wherever the anchor's inputs are exact: the two columns' floors at the same depth and pressure.
Elsewhere the anchor carries the difference between the two columns' `VertCoord` heights (§3.5.1,
§3.7.4), which is small and converges with layer thickness but is not round-off.

What makes this tractable, and is worth carrying away, is that condition 1 is now checked *pointwise*
rather than assembled: [](#dalpha) is zero at each quadrature point or it is not. There is no
cancellation between terms to get right, and consequently no way to satisfy the resting-state gate
while being subtly wrong — which the previous formulation permitted (§6.3).

#### 3.7.3 Which water columns cancel exactly

Condition 1 holds whenever the true vertical profile is one the reconstruction reproduces exactly:
two mean-preserving reconstructions of such a profile, built on *different* sets of interfaces from
*different* layer means, recover the identical profile anyway. Together with the edge-shared
expansion point [](#edge-ref) this gives, for each phase:

| Vertical profile of $\Theta$, $S$ | `Centered` | Phase 1 (linear reconstruction) | Phase 2 (parabolic reconstruction) |
|---|---|---|---|
| Uniform ($\alpha$ constant) | exact | exact | exact |
| Uniform $\Theta$, $S$; compressibility only | $O(\tilde h)$ | **exact** | **exact** |
| Linear in pressure | $O(\tilde h)$ | **exact** | **exact** |
| Quadratic in pressure | $O(\tilde h)$ | $O(\tilde h^2)$ | **exact** |
| Realistic (general smooth) | $O(\tilde h)$ | $O(\tilde h^2)$ | $O(\tilde h^3)$ |
| Horizontal structure in $\Theta$, $S$ | $O(d_e^2)$ | $O(d_e^2)$ | $O(d_e^4)$ |

"Exact" means zero to machine precision, given exact anchor inputs: the two columns' floors at the
same depth and pressure (§3.7.2). Where they are not, every row also carries the anchor residual of
§3.5.1. Orders are in layer thickness $\tilde h$ at fixed tilt, and in cell spacing $d_e$ for the
last row. The $\tilde h$ entries carry a factor set by how strongly the profile curves and by how
far the interfaces are displaced between the two columns, so they shrink as either shrinks.

Two consequences deserve emphasis:

- Phase 1 is exact for profiles that are **linear in pressure**, not merely for vertically uniform
  ones. That distinction is the whole reason Phase 1 is worth doing on its own: a vertically uniform
  ocean does not exist, but linear-in-pressure captures compressibility exactly and smooth
  stratification to leading order, which is most of what a real column looks like over a single
  layer. It satisfies Requirement 2.3.3.
- `PressureGradCentered` carries a first-order-in-tilt error for *every* stratified profile,
  including a horizontally uniform one. Its order in layer thickness is profile-dependent (§6.5), so
  the advantage Phase 1 offers is best stated as an absolute error ratio at a given resolution rather
  than as a difference in convergence rate.

That last point is worth showing, since it identifies the error Phase 1 targets. Take specific
volume uniform in the horizontal within each layer, and let the coordinate merely redistribute
thickness between the two columns ($\sum_j \Delta \tilde h_j = 0$). The centered scheme's
edge-normal Montgomery-potential difference is then

$$
\Delta \bar M_k = g \rho_0 \sum_{j>k} \left( \alpha_j - \alpha_k \right) \Delta \tilde h_j,
$$ (centered-error)

which is nonzero whenever $\alpha$ varies with depth and interfaces are displaced, and which
**accumulates downward**: the deeper the layer, the more layers contribute to the sum. The vertical
contrast in *in-situ* $\alpha$ is dominated by compressibility (a few percent over the full ocean
depth, against a few tenths of a percent from $\Theta$ and $S$), and that is exactly the part
Phase 1 cancels.

`PressureGradCentered` is exactly first order in the **coordinate tilt**: measured fitted exponent
$1.0000$ at three vertical resolutions on the Polaris resting-state configurations, with the
`ztilde_gradient` variant converging at $\approx 1.1$ in horizontal resolution once the bottom layer
is included. Its order in **layer thickness** is profile-dependent and is not reliably first — it
reaches second order and beyond on smooth resting profiles (§6.5) — so the tilt exponent, not the
thickness order, is the robust statement. The corresponding absolute errors reach
$2\times10^{-5}\ \mathrm{m\,s^{-2}}$ at a coordinate tilt of 50 m/km with 256 m layers — the order of
the bottom-layer error seen in realistic global configurations. Whether this downward accumulation is in fact what produces the bottom-layer noise
seen in realistic global configurations is a plausible diagnosis, not an established one; it is
carried as A4 in §3.7.6 and tested in §5.3.

**The tilt exponent is a statement about `PressureGradCentered` and must not be used to compare the
two schemes.** Measured in Omega on the Polaris resting-state variants, `PressureGradFiniteVolume`
returns a fitted exponent of $1.000$ as well — the *same* value — while its absolute error is one to
three orders of magnitude smaller: $1.750\times10^{-10}$ against $2.602\times10^{-8}\ \mathrm{m\,
s^{-2}}$ at a tilt of 0.05 m/km, and $2.022\times10^{-7}$ against $2.069\times10^{-5}$ at 50 m/km,
both at 4 km spacing and 256 m layers.

That is not a failure of the scheme and it is not a coincidence. On a Polaris-initialized state the
finite-volume residual is entirely the anchor of [](#anchor) (§5.2), and the anchor is *linear in
the tilt* because the cross-edge difference it corrects to a common pressure is. So both schemes
give $q = 1$ for unrelated reasons: the centered scheme because $\alpha$ is effectively piecewise
constant in pressure within a layer, the finite-volume scheme because its one remaining error term
inherits the tilt's own linearity.

**Reading $q$ alone would therefore report "no better than centered" for a scheme that is up to
500 times more accurate.** For comparing the schemes the absolute level is the discriminating
measurement, and it is what §5.1's accuracy gate uses. The exponent retains its original purpose —
calibrating the centered baseline, and detecting a scheme that has become *insensitive* to tilt —
but it no longer separates the two.

#### 3.7.4 What `VertCoord` supplies, and what it does not

`VertCoord` is the source of `PressureInterface`, `PressureMid`, `PseudoThickness` and the two
values of `GeomZInterface` that [](#anchor-one-layer) differences at the sea floor. It is **not**
the source of the height the scheme
integrates: [](#d-recurrence) builds $\Delta_e z$ from the reconstruction alone, so the question of
whether `VertCoord`'s $z$ is built from the same $\hat\alpha$ the PGF uses does not arise.

Two questions that would otherwise arise therefore do not.

**The quadrature question.**
`VertCoord::computeGeomZHeight` accumulates $\Delta z = \rho_0\,\alpha_{i,k}\,\tilde h_{i,k}$ —
apparently a midpoint rule. For the **linearized** $\hat\alpha$ of §3.3 it is the exact layer
integral:

$$
\frac{1}{g}\int_{p^{\text{top}}_{i,k}}^{p^{\text{bot}}_{i,k}} \hat\alpha_{i,k}(p)\,dp
= \frac{\alpha_{0}\,\Delta p_{i,k}}{g} = \rho_0\,\alpha_{i,k}\,\tilde h_{i,k},
$$ (z-increment-exact)

because $\int \Theta'\,dp = \int S'\,dp = 0$ by the mean-preserving constraint and
$\int (p - p^{\text{mid}})\,dp = 0$ because `PressureMid` is the exact arithmetic midpoint of the two
interface pressures.

**That is a statement about $\hat\alpha$, not about $\alpha$, and the distinction matters.**
`VertCoord` does not accumulate $\hat\alpha$. It uses the exact TEOS-10 $\alpha$ evaluated at the
layer-mean state $(\Theta_{i,k}, S_{i,k}, p^{\text{mid}}_{i,k})$, and $\alpha$ is nonlinear in
pressure, so its layer average is not its value at the midpoint state. `VertCoord`'s increment is
therefore a second-order-accurate approximation of the true layer integral, not an exact one.
Measured, the gap accumulates to $1.1\times10^{-3}$ m over a column.

**No change to `VertCoord` and no answer-changing baseline step is required** — but not because
[](#z-increment-exact) makes its rule exact, which it does not. It is because the scheme does not
accumulate `VertCoord`'s $z$ at all. It uses it once, as the difference [](#anchor-one-layer) takes
at the sea floor, where over a flat floor with matching `maxLevelCell` that difference is
$-\Delta_e H$, exact prescribed input, and no accumulation enters it. Where the two columns' floors
differ the anchor does span a partial accumulation, over the deeper column's layers below the
interface it contributes, and the gap enters at the level of the difference between the two
columns' partial sums. §5.1's `bathymetry_step` and the sloping-floor test of §5.2 are where that is
measured.

**The sharing question.** Whether a cell-based $z$ can carry an edge-dependent $\hat\alpha$ would
matter to a formulation that differenced two separately accumulated column integrals. It does not
arise for one that differences the integrand: there is no per-column integral to reconcile, and the
per-layer mismatch between `VertCoord`'s increment and the edge-shared profile's is not part of the
scheme.

**Two exact simplifications** follow from taking the edge control volume as the average of the two
columns' interface pressures, $p^{\text{top}}_{e,k} = \bar q_k$ and
$p^{\text{bot}}_{e,k} = \bar q_{k+1}$, and are worth using in the implementation:

- $p^{\text{mid}}_{e,k}$ is **exactly** the edge average of the two columns' own layer midpoints;
- $\Delta p_{e,k}$ is **exactly** the edge average of the two columns' own layer thicknesses.

Both follow from `PressureMid` being the exact arithmetic midpoint, and both were confirmed
numerically.

**The scan anchors at the sea floor, and the two ends fail in complementary ways.** The choice is a
trade-off rather than a settled preference, and neither end is uniformly better:

| | sea-floor anchor | surface anchor |
|---|---|---|
| coordinate tilts, flat floor | exact; $\Delta_e Z = -\Delta_e H$ is exact input | fails; the two columns' derived sea-surface heights differ by $7.1\times10^{-8}$ m, matching the difference of their midpoint-rule truncations to all figures measured |
| bathymetry steps | forms and cancels quantities of order the step within the deepest shared layer; $1.7\times10^{-7}$ m residual at 200 m/km, nine digits consumed | forms nothing; $\Delta_e Z$ at the surface is zero and the correction is a zero-width integral |

**Requirement 2.3.1 asks for machine precision at any tilt, thickness *and* bathymetry, and neither
end delivers all three.** The sea floor is chosen because it is exact on the case the requirement was
written for — a resting ocean on an unstepped floor at arbitrary coordinate tilt — and because its
stepped-floor residual is comparable to the surface anchor's unconditional one rather than worse.
That is a judgement, and it should be revisited if stepped bathymetry turns out to be the case that
matters; §5.1's `bathymetry_step` variants are where it is measured.

The stepped-floor residual is **not** quadrature error: refining from 2 to 16 Gauss points does not
move it at all, which is conclusive rather than suggestive, since $\hat\alpha$ is linear in pressure
and every rule integrates it exactly. Two sources remain, and only one has been settled.

- **The expansion stretched across the floor offset.** Where the floors are more than a layer apart,
  this dominates. It was confirmed on terrain-following layers, where it grows with the cube of the
  offset, and it is what [](#anchor-one-layer) removes (§3.5.1).
- **`VertCoord`'s midpoint-rule sums over different layers.** On `bathymetry_step` the floors are
  within a layer of each other at the anchor, so [](#anchor-one-layer) leaves those answers bit for
  bit unchanged. There the anchor still integrates the expansion across up to a layer, some
  $2\times10^{6}$ Pa at 200 m/km, and also differences the deeper column's midpoint-rule sum over
  the layers below the anchor, which the shallower column does not have. Which of the two dominates
  there is not measured. Using the exact equation of state in the anchor's short integrals was
  tried in a one-edge model of terrain-following layers and is not a clean fix: at gentle slopes it
  is worse, because the midpoint-rule difference dominates there. The second-order expansion of
  §3.3 has not been tried in the anchor.

The argument for the sea floor on a tilted coordinate is unchanged and follows.

**Why the sea floor wins there.** `VertCoord` builds $z$ upward from the bathymetry while
pressure is built downward from the surface, so the two accumulate round-off from opposite ends, and
[](#ssh-difference) favours the sea floor on conditioning grounds: $-\Delta_e H$ is exact input
where $\Delta_e Z_{\text{surf}}$ is the small residual of two column-length sums.

**Conditioning is not the deciding argument, though — exactness is.** `VertCoord` accumulates
$\Delta z = \rho_0\,\alpha_{i,k}\,\tilde h_{i,k}$ over each column's *own* layers, and that
increment is not the exact layer integral of the true $\alpha$ (previous subsection). The error
comes from the equation of state's curvature in pressure rather than from anything about the
profile, so it is present **even on the exact set**, and two columns with different layer partitions
accumulate different amounts of it. Their derived sea-surface heights therefore differ even when the
two columns hold identical water. Anchored at the surface that difference enters $D_1$ directly and
the machine-precision property of §2.3.1 is lost; anchored at the sea floor over a flat floor
$\Delta_e Z$ is exact input and vanishes identically, and none of it enters — because the scheme
integrates its own reconstruction rather than accumulating `VertCoord`'s height. The two ends are
therefore *not* equivalent, and the earlier claim that they agree to round-off does not hold, not
even on the exact set: the two columns' midpoint-rule errors differ there too, because their layer
partitions do.

One consequence is worth recording for whoever builds the reference implementation. A test harness
that anchors its geometric column at a *prescribed* sea surface and derives the bathymetry is
constructing the state from the opposite end to `VertCoord`, which prescribes the bathymetry and
derives the sea surface. The two agree only if the derived quantity is what is written to the initial
condition and the same $\Delta z$ increment is used on both sides — a round trip that is easy to
break silently and that nothing in this design would catch. That the existing Omega-vs-Polaris
comparison holds to $10^{-16}$–$10^{-13}$ on the centered scheme, which reads `GeomZInterface`
directly, is the evidence that it currently holds.

#### 3.7.5 Round-off in the deep ocean

A scheme that differences two column integrals of order $10^{2}$ m to obtain a residual of order
$10^{-3}$ m consumes roughly five significant digits before the physics appears, and roughly ten by
the time the tendency is formed. That is `PressureGradCentered`'s arithmetic, and it is the reason
this section exists.

**[](#dz-dp) removes the exposure rather than managing it — in the scan.** Every quantity in the
scan is a difference *before* it is an integral: [](#dalpha) is a horizontal contrast,
[](#d-recurrence) accumulates integrals of that contrast, and $D_k$ is the answer itself. No large
quantity is formed there, so nothing large has to cancel, and the scan's conditioning is set by the
size of the baroclinic signal rather than by the size of the hydrostatic terms. Measured, $D_k$
stays flat to $1.4\times10^{-15}$ m at every floor gradient tried, including where one column's
reconstruction is evaluated three layers below its own floor.

**The anchor is not covered by that argument.** An earlier version of this section claimed no large
quantity is formed *anywhere*; that is true of [](#d-recurrence) and false of [](#anchor) at the sea
floor, where a stepped bathymetry makes $\Delta_e Z$ of order the step and the short integrals
cancel it — nine to twelve digits, per the measurements in §3.5.1. The exposure the scheme removes
is the column-length one; a step-sized cancellation at a single interface remains.

Three consequences.

- **The ten-digit estimate does not apply to this scheme.** It applies to `PressureGradCentered`,
  which still forms and cancels $10^{2}$ m quantities, and to any formulation that reconciles two
  separately accumulated column integrals. Measured, restructuring the accumulation this way
  recovered about three decimal digits of headroom on a tilted coordinate — the largest increment
  being $\sim 10^{-3}$ of the height difference it replaces (§6.3).
- **A single-precision build passes §5.2.** Measured in Omega: the exactness gate returns
  $3.0\times10^{-9}$ m s$^{-2}$ in single precision against $3.1\times10^{-18}$ in double, which is
  $0.09\,\epsilon$ of the cancelling terms in one build and $0.05\,\epsilon$ in the other — the
  same result, scaled by the precision. `PressureGradCentered` on the same state is
  $3.4\times10^{-5}$ m s$^{-2}$ in single and $3.5\times10^{-5}$ in double, essentially unmoved,
  because its error is truncation rather than round-off. The scheme forms no large quantities, and
  the reduced precision costs it nothing beyond $\epsilon$ itself.
- **The perturbation form is not needed and is not specified.** It would have subtracted a local
  reference profile from $\hat\alpha$ before integrating and added its contribution back
  analytically. It addressed a problem that differencing first removes algebraically. If a
  single-precision measurement ever contradicts the paragraph above, this is the direction to look,
  but it should not be built speculatively.

#### 3.7.6 Assumptions that still need testing

Everything above concerns properties the implementation can be built to have. The following are
assumptions this design is making that only testing can confirm:

- **A1 — Accuracy at cell centers does not imply accuracy at edges.** Reconstructing $\Theta$ and
  $S$ accurately at each cell center does not by itself make $z_i(p)$, and hence
  $[\nabla_n z]_p$, accurate at the edge between them, because $\hat\alpha(\Theta,S,p)$ is a
  nonlinear function of the reconstructed quantities. The horizontal order claimed in §3.6 is
  therefore a design target to be confirmed by measured convergence (Requirement 2.6), not something
  that follows from the cell-centered reconstruction order. Note that A1 is now purely an *accuracy*
  assumption: the reformulation of §3.2 removes its consequences for the machine-precision property,
  which is why §3.6.1's stencil constraint could be dropped.
- **A2 — The EOS expansion is good enough across the edge.** [](#alpha-taylor) is expanded about the
  shared edge state [](#edge-ref), and its error grows with the horizontal contrast in $\Theta$ and
  $S$ across the edge. Whether the first-order expansion suffices in frontal regions, or whether the
  second-order option of §3.3 is needed, is an open question — though a narrower one than it was, since
  [](#dalpha) makes the coefficients multiply a quantity that is zero on the exact set, so A2 bears on
  accuracy only and not on robustness. Measured, the second-order remainder relative to $\alpha$ is
  $8\times10^{-6}$ on the hardest configuration available (a 12 °C contrast across a 4 km edge) and
  $\le 2\times10^{-9}$ on the resting-state variants (§6.3). The first-order expansion is not being
  worked near its limit, and `temperature_gradient` at coarse resolution is the configuration to watch.

  **A2 must be measured on the remainder itself, not inferred from a sweep of the contrast.** The
  contrast cannot be varied independently of the cell spacing in a two-column configuration, and
  A2's error is in any case confounded with the horizontal stencil's, which scales the same way and
  dominates. §5.1 gives the measurements behind both statements.
- **A3 — The residual outside the exact set is small enough in practice.** The $O(\tilde h^2)$ and
  $O(\tilde h^3)$ entries in §3.7.3 describe how the error *scales*; how large it actually is at the
  vertical resolutions Omega can afford is unknown.
- **A4 — PGF error is what is driving the observed instability.** Identifying [](#centered-error)
  with bottom-layer noise in realistic global runs is a plausible diagnosis, not a demonstrated one.
  A configuration whose profile falls in Phase 1's exact set gives a direct test: if spurious
  bottom-layer flow persists there, the cause lies elsewhere — most likely in the layer-mean
  treatment in the tracer and remapping operators (§4.3) — and Phase 1 will not cure it. See §5.3.

  **Measured since Phase 1 landed, and supported over the open ocean** (§5.3.2). Swapping
  `Centered` for `FiniteVolume` in a 10-day realistic global run cuts the area-mean bottom-layer
  kinetic energy sixteenfold and removes the grid-scale noise across the open ocean. What remains is
  mostly on shelves and in marginal seas, where A4 is not established.
- **A5 — The treatment at the top and bottom of the column is adequate.** [](#ho-exact) averages over
  the edge control volume's pressure range, which in general is not either column's own, so near the
  sea floor — where the two columns may have different `maxLevelCell`, and where a partial cell carries
  the whole signal for the Polaris `bathymetry_step` configuration — a column's deepest reconstruction
  must be evaluated slightly below its own floor.

  **A5 is an accuracy assumption, not a robustness one.** [](#dalpha) is evaluated at matched pressure
  whether or not the evaluation point lies inside the layer the reconstruction was built for, and on
  the exact set an extrapolated reconstruction still reproduces the true profile, so $\Delta_e\Theta(p)$
  is still zero and exactness survives. What extrapolation costs is accuracy off the exact set. The
  derivation harness could not measure it, having equal layer counts in both columns; §5.1's
  `bathymetry_step` variants have since measured it as benign (§6.5, §6.6).

  **The rule is to clamp to the outermost valid layer of the column and extrapolate that layer's
  reconstruction.** A pressure above a column's shallowest valid layer uses that layer; a pressure
  below its deepest valid layer uses that one. This is stated here rather than left to fall out of
  the implementation, and it is the rule the per-column pressure lookup applies at both ends.

### 3.8 Scope: the reduction and the vertical coordinate

The reduction of §3.1.1 rests on one property of Omega — that $\tilde z$ is a function of pressure
alone — so it is worth recording precisely how much freedom that leaves and where the boundary is.

**It holds for any ALE layering of pseudo-height.** p-star, sigma-p, z-tilde-star and hybrid schemes
all relocate layer interfaces within $\tilde z$, in space and in time, and none of them redefines
$\tilde z$. Interfaces may be arbitrarily tilted, arbitrarily thin, and at wholly different pressures
in neighbouring columns; [](#grad-p-zero) is untouched, because it is a statement about coordinate
surfaces (§3.1.1). Every term the scheme computes is a tilt term, so there is no configuration in
which the reduction quietly turns the scheme into a no-op.

**It fails only if Omega abandons pseudo-height.** A coordinate $r$ that is not a function of pressure
alone — geometric $z$, a true $z$-star, or a genuine isopycnal coordinate — gives $\nabla_r p \neq 0$,
and then all four pressure terms of [](#ho-target) become real work rather than a discretization of
zero. Should such a change ever be contemplated, [](#ho-target) as corrected here and the derivation in
§3.1 are the starting point, and §3.1.1's third consequence is the pitfall to carry forward with them.

**The geopotential body force** — the first line of [](#ho-target) — is the whole target. It is
$g\,\nabla z + \nabla(\phi_{TP}+\phi_{SAL})$; $z$ enters as `VertCoord`'s `GeomZInterface`, through
[](#centered-shift) in the centered scheme and only at the anchor [](#anchor-one-layer) in this
one, and the tidal-potential and self-attraction-and-loading contributions are differenced with
[](#edge-grad) (Requirement 2.7). They are placeholder arrays held by `PressureGrad`, zero until a
tidal forcing module supplies them, and are unaffected by anything in §3.5.

### 3.9 Relationship to the centered scheme

`PressureGradCentered` is **exactly** the first-order conversion of a height difference taken at
fixed layer index into one taken at fixed pressure. This is an algebraic identity about the
implemented code rather than an approximate correspondence, and it has been confirmed numerically to
$0.5\,\epsilon$ of the hydrostatic scale over fifty states, including the bottom partial cell and
configurations where the two columns' `maxLevelCell` differ (§6.2).

Let $q_{i,k}$ = `PressureInterface`, $Z_{i,k}$ = `GeomZInterface`, $\alpha_{i,k}$ = `SpecVol`, and

$$
\mathcal{S}_{e,k} \;=\; \tfrac{1}{2}\left(\Delta_e Z_{k} + \Delta_e Z_{k+1}\right)
\;+\; \frac{\bar\alpha_{e,k}}{2g}\left(\Delta_e q_{k} + \Delta_e q_{k+1}\right),
\qquad
\bar\alpha_{e,k} = \tfrac{1}{2}\!\!\sum_{i \in CE(e)}\!\! \alpha_{i,k}.
$$ (centered-shift)

Then the centered tendency is exactly

$$
T^{p,\text{ctr}}_{e,k} = -\frac{g}{d_e}\,\mathcal{S}_{e,k} - \left[\nabla_n(\phi_{TP}+\phi_{SAL})\right]_e .
$$ (centered-identity)

The proof is two lines: with
$M = gZ + \alpha q$, the identity $\Delta_e(\alpha q) = \bar\alpha\,\Delta_e q + \bar q\,\Delta_e\alpha$
splits the Montgomery difference, and the $\bar q\,\Delta_e\alpha$ halves cancel the scheme's
$\bar p\,\nabla\alpha$ term **term for term**, because
$\tfrac12(\bar q_{k} + \bar q_{k+1}) = \bar p^{\text{mid}}_{e,k}$ exactly — `PressureMid` is the exact
arithmetic midpoint (`VertCoord.cpp:1067`).

Three things this buys, and one it does not.

- **It explains the error the centered scheme makes.** Taylor-shifting each column's height from its
  own interface pressure to a common one gives
  $z_i(p) = Z_{i,k} - \tfrac{\alpha}{g}(p - q_{i,k}) + O((p-q_{i,k})^2)$, so $\mathcal{S}$ is the
  first-order truncation of what [](#ho-exact) asks for. `PressureGradCentered`'s entire
  Montgomery-potential apparatus is that conversion and nothing else, and its error is the
  conversion's truncation — which is why §3.7.3's first-order behaviour is what it is, and why
  compressibility appears in that error even though it cannot appear in the true horizontal gradient
  (§3.5, consequence 4).
- **It gives a regression test with a known answer.** §5.5 asserts [](#centered-shift) against
  `PressureGradCentered` directly. Because the two functors read the mesh, `VertCoord` and EOS state
  through independently written code, their agreement tests the shared upstream state — edge masks,
  interface indexing, `VertCoord` conventions — and not just the PGF arithmetic (§4.4).
- **It bounds the difference between the two schemes.** They agree to
  $O\!\left((\Delta_e p)^2\right)$, so a `FiniteVolume` answer far from `Centered` at small tilt
  indicates a defect in one of them rather than a property of the physics.

**What it does not give is a configuration switch that reduces one scheme to the other.**
[](#dz-dp) is not built by adding a correction to [](#centered-shift); it differences the integrand
and never forms $\mathcal{S}$ at all. There is therefore no `FiniteVolume` setting that recovers
`PressureGradCentered` bit-for-bit or to round-off, and §5.5 tests the identity above rather than a
reduced configuration — which is the more useful test in any case, since it can be asserted without
running the new scheme at all.

### 3.10 Per-step algorithm summary

Steps 1–3 are per cell and layer; step 4 is a per-edge column scan; step 5 is per edge and layer.
Phase differences are marked.

1. From `VertCoord`: read `PressureInterface`, `PressureMid` and `GeomZInterface`, the last only
   at the anchor (already computed diagnostically each step). `VertCoord` itself needs no change
   and no baseline step (§3.7.4).
2. Obtain $\alpha_0$ (= the existing `Eos::SpecVol` field) and the derivatives
   $\alpha_\Theta$, $\alpha_S$, $\alpha_p$ from one TEOS-10 evaluation ([](#alpha-derivs)), in
   `AuxiliaryState` before the tendencies are formed. $\alpha_p$ is not required by [](#dalpha) and
   is needed only for the anchor [](#anchor-one-layer); Phase 2 optionally adds second derivatives
   (§3.3).
3. Build the mean-preserving deviations $\Theta'$, $S'$ ([](#vert-recon)) — **linear in Phase 1,
   parabolic in Phase 2** — using the actual non-uniform interface pressures (§3.4).
4. For each edge, in one column scan anchored at the sea floor: for each edge layer, in parallel,
   form the shared coefficients [](#edge-ref) and evaluate [](#dalpha) at the layer's quadrature
   points — each column's $\Theta$, $S$ taken from *its own* layer containing that pressure — to
   form the layer's increment of [](#d-recurrence) and its first moment; then form the anchor
   [](#anchor-one-layer) and accumulate the increments upward from it to obtain $D_k$.
5. For each edge and layer: form the layer mean of $\Delta_e z$ from $D_{k+1}$ and the first moment,
   assemble [](#ho-exact) with the tidal/SAL difference, and accumulate into `Tend` with
   `EdgeMask`. **Phase 1** uses the two-cell operator [](#edge-grad); **Phase 2** uses the wider
   stencil of §3.6.1.

## 4 Design

The `FiniteVolume` scheme extends the existing `PressureGrad` class in
`components/omega/src/ocn/PGrad.{h,cpp}`. The class manager, creation/retrieval/removal
methods, named-instance map, and dispatch logic are unchanged; the work is to flesh out the
`PressureGradFiniteVolume` functor (a no-op placeholder before Phase 1), to add the kernels that
feed it, and to widen the inputs the manager hands to it. This section describes the code as Phase 1
implements it.

### 4.1 Data types and parameters

#### 4.1.1 Parameters

Configuration lives under the existing `PressureGrad` YAML group. The `PressureGradType`
enum in `PGrad.h` gives the new scheme a descriptive name and drops the unimplemented second
placeholder, leaving a commented stub for a future variant. A second enum names the vertical
reconstruction, with a commented stub for Phase 2's:

```c++
enum class PressureGradType {
   Centered,    ///< existing 2nd-order Montgomery scheme
   FiniteVolume ///< layer-integrated finite-volume scheme
   // , <FutureVariant>  ///< e.g. a 6th-order option, added when implemented
};

enum class PressureGradVertRecon {
   Linear ///< linear deviations
   // , PPM  ///< parabolic (PPM-style) deviations, added when implemented
};
```

New sub-options for `FiniteVolume`. Both phases use the same keys; Phase 2 adds values
rather than keys, so no configuration written for Phase 1 needs to change when Phase 2 lands.
`Default.yml` carries:

```yaml
    PressureGrad:
       PressureGradType: FiniteVolume     # Centered | FiniteVolume
       HorzOrder: 2                       # 2 = two-cell stencil (Phase 1); 4 = wide stencil (Phase 2)
       VerticalReconstruction: linear     # 'linear' (Phase 1) | 'ppm' (Phase 2)
       QuadraturePoints: 2                # per-edge-layer points for the matched-pressure integrand
```

- **Phase 1 default and target:** `HorzOrder: 2`, `VerticalReconstruction: 'linear'`,
  `QuadraturePoints: 2`.
- **Phase 2 target:** `HorzOrder: 4`, `VerticalReconstruction: 'ppm'`.

The three sub-options are optional: where a key is absent the Phase 1 default is used, so a
configuration written before they existed still runs. They are read and validated in the
`PressureGrad` constructor and cached on the `PressureGradFiniteVolume` functor. **Values that are
not implemented are rejected with an error rather than falling back** — today any `HorzOrder` other
than 2, and `VerticalReconstruction: 'ppm'` — so that a configuration requesting Phase 2 cannot
quietly run as Phase 1. `QuadraturePoints` must lie between 1 and 4 (`MaxPGradQuadPoints`). An
unrecognized `PressureGradType` is likewise a fatal error rather than a fallback to `Centered`,
which would otherwise turn a typo into a run that looks like a passing centered run.

`QuadraturePoints` sets the number of points at which [](#dalpha) is evaluated within each edge
layer. **It is an accuracy setting only.** Because the integrand is zero *pointwise* on the exact set
(§3.5, consequence 2), no quadrature rule can break the robustness property, and the choice therefore
trades cost against accuracy off the exact set with nothing else at stake. Two points is exact for
the piecewise-linear Phase 1 integrand within a sub-interval and is the default; more is worth
considering only where the two columns' interfaces are strongly offset, since the integrand is
piecewise linear with breakpoints at the union of the two columns' interfaces and a fixed rule does
not resolve those breaks.

`HorzOrder` selects the width of the edge stencil used by [](#edge-grad). It is purely an accuracy
setting: it carries no constraint from the robustness property (§3.6.1).

There is no setting that reduces the scheme to `PressureGradCentered`: [](#dz-dp) differences the
integrand rather than correcting the centered form, so no term can be switched off to recover it
(§3.9). `FiniteVolume` names the layer-integrated
control-volume form of §3.1 rather than a face-by-face assembly of the pressure traction, which is
the narrower sense Adcroft et al. (2008) use.

The configuration that isolates compressibility on its own is `'linear'` applied to a vertically uniform
$\Theta$, $S$ profile — the reconstruction slopes are then zero and $\alpha_p$ is retained — and it
needs no setting of its own. Tests that mean to isolate compressibility, including §5.2, select
`'linear'`.

#### 4.1.2 New EOS support

`PressureGradFiniteVolume` needs $\alpha$ together with its first derivatives. Note that
[](#dalpha) uses only $\alpha_\Theta$ and $\alpha_S$ — $\alpha_p$ cancels in the matched-pressure
difference and is required only by the anchor [](#anchor-one-layer). The three derivative arrays are provided
together regardless, since `Eos` computes them from one evaluation and `BruntVaisalaFreqSq` already
consumes $\alpha_p$. The `Eos` class
(Eos.h) exposed `computeSpecVol`, `computeSpecVolDisp`, and
`computeBruntVaisalaFreqSq`, but no specific-volume derivatives. One method and three device
arrays were added ahead of Phase 1, in [Omega#485](https://github.com/E3SM-Project/Omega/pull/485):

```c++
// On the Eos class
Array2DReal SpecVolDCt; ///< d(alpha)/d(ConservTemp) at cell centers
Array2DReal SpecVolDSa; ///< d(alpha)/d(AbsSalinity) at cell centers
Array2DReal SpecVolDP;  ///< d(alpha)/d(Pressure)    at cell centers

// Pressure is relative pressure in Pa; the derivatives are returned per degC,
// per (g/kg), and per Pa respectively, into the members above and SpecVol.
void computeSpecVolAndDerivs(const Array2DReal &ConservTemp,
                             const Array2DReal &AbsSalinity,
                             const Array2DReal &Pressure);
```

`Eos` owns the derivative arrays exactly as it owns `SpecVol`, `SpecVolDisplaced` and
`BruntVaisalaFreqSq`: they are allocated in the constructor and registered as `Field`s in the
`Eos` group, so they can be written to a stream and any number of consumers can read them
without each allocating its own copy. Their valid range must permit negative values —
$\alpha_S < 0$ everywhere, and $\alpha_\Theta < 0$ in cold fresh water — unlike `SpecVol`.

Because `computeSpecVolAndDerivs` fills `SpecVol` as well, it replaces rather than accompanies
a call to `computeSpecVol`; the two are kept separate so that the derivative arithmetic is
paid only where it is needed. `AuxiliaryState::computeMomVertAux` makes that choice: it calls
`computeSpecVolAndDerivs` when the default `PressureGrad` instance is `FiniteVolume` and
`computeSpecVol` otherwise, querying the default instance because the scheme is a model-wide choice
and `AuxiliaryState` holds no `PressureGrad` of its own.

The pressure gradient receives the equation of state as a `const Eos *`, so it cannot call any of
these compute methods, and the per-edge functions read only the four cached fields. That interface,
not a run-time count, is what keeps Requirement 2.2 in force; §5.6 gives the reason.

The TEOS-10 derivatives are obtained by differentiating the 75-term polynomial analytically.
The pressure derivative reuses the coefficients `calcPCoeffs` already assembles for
$\alpha$ itself and so is free; the $\Theta$ and $S_A$ derivatives need their own coefficient
sets, built from the same $s$ and $t$ (no additional square root). The marginal cost over
`computeSpecVol` is therefore roughly two extra coefficient assemblies and no additional
TEOS-10 evaluation, which is what Requirement 2.2 bounds. The linear and constant EOS
options supply trivial analytic derivatives.

#### 4.1.3 `PressureGradFiniteVolume` functor

The functor mirrors `PressureGradCentered` (cached mesh/coordinate arrays, `Enabled` flag,
`chunkStart`/`chunkLength` vertical iteration). Almost all of the work of forming the fixed-pressure
comparison happens before it runs, in the column scan below, so what it takes is that scan's two
output arrays rather than the state the scan consumed:

```c++
KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                const Array2DReal &PressureInterface,
                                const Array2DReal &DeltaZFixedP,
                                const Array2DReal &DeltaZMoment,
                                const Array1DReal &TidalPotential,
                                const Array1DReal &SelfAttractionLoading)
                                const;
```

(The existing centered signature is unchanged.) Additional cached members hold, in Phase 2, the
wide-stencil cell lists and weights. The functor implements step 5 of §3.10 per edge and vertical
chunk, accumulating into `Tend` with `EdgeMask`, exactly as the centered functor does.

Three aspects of the loop structure differ from the obvious implementation:

- **Column quantities built from $\hat\alpha$ are computed per edge, not cached per cell.** Both
  adjacent columns' contributions use equation-of-state coefficients formed at the edge
  [](#edge-ref), so the same cell yields a different result at each of its edges. The per-cell
  quantities that *can* be cached are the four EOS coefficients and the reconstruction slopes, which
  is where the TEOS-10 cost lives; the per-edge work is polynomial arithmetic on those cached
  values. This is what keeps Requirement 2.2 satisfied despite roughly three times as many column
  evaluations on a hexagonal mesh.
- **The column scan of §3.10 step 4 cannot live in this functor.** [](#d-recurrence) is a prefix sum
  up each column with edge-dependent coefficients, so it is not expressible as an independent
  per-vertical-chunk operation. `PressureGrad::computeColumnScan` computes it in a separate kernel
  filling three per-edge arrays,

  ```c++
  Array2DReal DeltaZIncr;    ///< (NEdgesSize, NVertLayers) — the recurrence's increment per layer
  Array2DReal DeltaZMoment;  ///< (NEdgesSize, NVertLayers) — its first moment about the top interface
  Array2DReal DeltaZFixedP;  ///< (NEdgesSize, NVertLayersP1) — D_k at edge interfaces
  ```

  as a `parallelForOuter` over edges, each edge a team that makes two passes over the edge's active
  layers. In the first, a `parallelForInner` over the layers forms each layer's shared expansion
  and evaluates the matched-pressure integrand at the layer's quadrature points, accumulating the
  increment of [](#d-recurrence) and its first moment from the same points. After a `teamBarrier`,
  the second forms the anchor [](#anchor-one-layer) and turns the increments into $D_k$ with a
  `parallelScanInner` upward from it, in the same shape as `VertCoord::computeGeomZHeight`. The
  functor then reads `DeltaZFixedP` and `DeltaZMoment` chunk-wise like any other input.

  $D_k$ is held at interfaces, so its extent is `NVertLayersP1`. The moment
  $\int (p - \bar q_k)\,\Delta_e\hat\alpha\,dp$ is what the layer mean of [](#ho-exact) needs
  alongside $D_{k+1}$; forming it in the same pass over the same quadrature points is what §3.5.1
  requires, and it is why the functor reads neither $\Theta$, $S$ nor the equation-of-state
  derivatives. Evaluating the integrand a second time in the functor to recover the moment would
  double the cost of the most expensive part of the scheme.

  The layer mean is then
  $D_{k+1} + \frac{1}{\Delta p_{e,k}}\int (p - \bar q_k)\,\Delta_e\hat\alpha\,dp$, taken from
  the layer's *lower* interface because that is the end nearer the anchor, with $\Delta p_{e,k}$ the
  edge average of the two columns' own layer thicknesses (§3.7.4).

  The reconstruction slopes are formed before the scan by `PressureGrad::computeReconSlopes`, a
  `parallelForOuter` over cells with a `parallelForInner` over each column's valid layers, into two
  per-cell arrays `ReconSlopeCt` and `ReconSlopeSa` (`NCellsSize`, `NVertLayers`).

  The cost is three edge-sized 2-D arrays, two cell-sized ones, and one column scan per edge per
  step. This is the one structural addition Phase 1 makes beyond the per-edge, per-chunk pattern the
  centered scheme uses. The arrays are allocated only when `FiniteVolume` is selected, so a
  `Centered` run pays nothing for them. Both kernel launchers are public members, which `nvcc`
  requires of a function containing a device lambda, and `getDeltaZFixedP` exposes $D_k$ so that
  tests can assert on it directly.
- **Each column's $\Theta$, $S$ are looked up by pressure, not by layer index.** At each quadrature
  point of edge layer $k$, [](#dalpha) needs the reconstruction of whichever of *that column's* layers
  contains the point, which under tilt is generally not layer $k$ (§3.7.2, condition 1).
  `findLayerForPress` walks from a hint, the edge layer index, up or down until the column's
  interfaces bracket the pressure, and clamps to the column's outermost valid layer outside it
  (§3.7.6, A5). Under tilt the walk is a few steps, and its result does not depend on the hint. It
  is a real difference from the obvious implementation, and getting it wrong reintroduces exactly the
  defect condition 1 warns about.

The arithmetic lives in two headers of inline functions, templated on the array type so that the
host-side reference assembly in the unit tests (§5.2) calls the same code the device kernels do:

| Header | Contents |
| ------ | -------- |
| `PGradRecon.h` | the reconstruction stencil and slope (`linearReconStencil`, `linearReconSlope`), its evaluation at a pressure (`linearReconEval`), and the per-column pressure lookup (`findLayerForPress`) |
| `PGradFiniteVolume.h` | the shared expansion (`PGradEdgeEos`, `buildEdgeEos`), one column's expanded $\hat\alpha$ for the anchor (`edgeSpecVol`), the matched-pressure integrand [](#dalpha) (`matchedPressSpecVolDiff`), and the Gauss-Legendre rule (`gaussLegendreRule`) |

`buildEdgeEos` takes a layer index for each cell. Within an edge layer the two are the same; the
anchor pairs layer $K$ of the shallower column with layer $j-1$ of the deeper one
([](#anchor-one-layer)).

The Phase 1 and Phase 2 code paths differ only in the reconstruction degree (§3.4) and in the width
of the edge stencil (§3.6.1). There is one functor, not two.

### 4.2 Methods

`PressureGrad::computePressureGrad` keeps its role of selecting the configured option, but
its input list grows so the `FiniteVolume` branch can reach $\Theta$, $S$, and the EOS
derivatives. The signature becomes:

```c++
void computePressureGrad(Array2DReal &Tend, const Array2DReal &PressureMid,
                         const Array2DReal &PressureInterface,
                         const Array2DReal &SpecVol,
                         const Array2DReal &GeomZInterface,
                         const Array2DReal &PseudoThick,
                         const Array2DReal &ConservTemp,        // new
                         const Array2DReal &AbsSalinity,        // new
                         const Eos *EqState) const;             // new: derivative fields
```

The `Centered` branch ignores the new arguments and is byte-for-byte unchanged; the
`FiniteVolume` branch runs `computeReconSlopes` and `computeColumnScan`, then dispatches to
`PressureGradFiniteVolume` through the same `parallelForOuter`/`parallelForInner` team pattern the
centered branch uses (PGrad.cpp). The call site in `Tendencies.cpp` is updated to pass $\Theta$,
$S$ (subviews of the tracer array) and the `Eos` instance, which it already references for
`SpecVol`.

Creation, retrieval, and removal (`init`, `create`, `get`, `getDefault`, `clear`, `erase`)
are unchanged. The constructor additionally reads and validates the `FiniteVolume` options
(§4.1.1), caches the mesh and coordinate arrays the column scan needs (`CellsOnEdge`,
`MinLayerCell`, `MaxLayerCell`), and allocates the scheme's working arrays when it is selected.
The Gauss-Legendre nodes and weights are formed inside the kernel from `QuadraturePoints`.
Accessors `getType`, `getHorzOrder`, `getVertRecon` and `getQuadraturePoints` report the
configured choice, which is how `AuxiliaryState` decides whether to compute the derivatives
(§4.1.2).

### 4.3 Consistency and follow-up work

This PGF is deliberately more accurate in the vertical than the rest of Omega's
layer-mean operators. This is sanctioned by {ref}`omega-design-governing-eqns-omega1` §9,
which singles out the PGF as **the** exception to the piecewise-constant assumption ("we will
ignore most of these [within-layer deviation] terms in Omega. The exception is the pressure
gradient force … appropriate for the *simple* pressure gradient force targeted for *early*
versions of Omega. This assumption will be revisited at a later date."). The mean-preserving,
`VertCoord`-anchored constraints (Requirement 2.4) keep the PGF consistent with the model
state at the layer-mean level. For completeness, the following terms would need higher-order
vertical treatment to make the **whole** model consistent with this PGF; they are explicitly
out of scope here and flagged as follow-up:

- **Vertical remapping and advection** ({ref}`omega-design-vert-adv` and the ALE
  remap): a mean-preserving PPM remap consistent with §3.4.
- **EOS evaluation point**: buoyancy and $N^2$ ({ref}`omega-design-eos`) are currently
  evaluated at layer mid; consistency with the reconstructed profile would evaluate them on
  the same reconstruction.
- **Dropped within-layer products** $\overline{\delta\varphi\,\delta u}$ in the tracer and
  momentum equations ({ref}`omega-design-governing-eqns-omega1` §9), which become non-negligible
  only at much higher resolution.

### 4.4 Retention of the centered implementation

`PressureGradCentered` (PGrad.h:25) is deliberately **kept as a separate functor** rather
than reimplemented as the lowest-order configuration of `PressureGradFiniteVolume`, even though
§3.9 shows the latter reduces to it. The redundancy is small — the centered functor is a
header-only, ~40-line loop body with no supporting machinery of its own — and it buys two
things that a single implementation cannot provide:

1. **An independent cross-check.** The two functors read the mesh, `VertCoord`, and EOS state
   through separately written code. Their agreement to round-off (§5.5) therefore tests the
   shared upstream state — edge masks, interface indexing, `VertCoord` conventions — and not
   just the PGF arithmetic. Collapsing them would make that comparison self-referential: a
   defect upstream of the order switch would appear identically in both limits and cancel.
2. **A fallback.** `FiniteVolume` is not a reduced form of `Centered` (§3.9), so no setting of
   the new scheme recovers the old answers. Keeping the centered functor keeps those answers one
   configuration key away, for comparison against earlier runs and for any configuration where
   the new scheme misbehaves, such as the steep sigma layers of §5.3.1.

**`FiniteVolume` is now the default.** Phase 1 made it so in the same change that introduced it,
after a 10-day realistic global run showed it removing the grid-scale bottom-layer noise that
`Centered` produces across the open ocean (§5.3.2). That was answer-changing for every
configuration with a horizontal density gradient: tasks that switch the pressure gradient off are
bit-for-bit, and tasks with no horizontal density gradient differ only at round-off.

Removing `PressureGradCentered` becomes reasonable once `FiniteVolume` has served out a period as
the default. That is deliberately left as **follow-up work**, not part of this design.

### 4.5 What each phase delivers and depends on

#### 4.5.1 Phase 1 — a consistent second-order scheme

**Delivers.** The exact fixed-pressure geopotential difference [](#ho-exact), evaluated by
differencing the integrand at matched pressure and integrating the difference (§3.5); the
equation-of-state expansion about a state shared across each edge (§3.3.1); mean-preserving linear
reconstruction of $\Theta$ and $S$ in pressure (§3.4); the per-edge column scan [](#d-recurrence) and
its sea-floor anchor, kept within one layer ([](#anchor-one-layer)). The result is a pressure
gradient that is zero to machine precision for any resting ocean whose profile varies linearly with
pressure, at any tilt and thickness, wherever the two columns' floors share depth and pressure
(§3.7.2, §3.7.3). Over a sloping or stepped floor a residual from `VertCoord`'s column sums remains,
small and convergent with layer thickness (§3.5.1). `FiniteVolume` is Omega's default (§4.4).

**Does not deliver.** Fourth-order accuracy. The horizontal operator is the same two-cell stencil
the centered scheme uses (§3.6), so horizontal truncation error is unchanged from
`PressureGradCentered`. Requirements 2.1 and 2.6 are met at second order only.

**Depends on.** No change to `VertCoord` and no baseline step (§3.7.4); the scheme does not
accumulate `VertCoord`'s geometric height at all, using `GeomZInterface` once, at the anchor
([](#anchor-one-layer)). All three conditions of §3.7.2 are properties of the PGF's own
discretization; none is a dependency on another module. The `Eos` derivatives it reads were added
separately, ahead of it (§4.1.2).

**Settled at implementation time.** The scan anchors at the **sea floor** (§3.7.4): anchoring at the
surface admits into $D_1$ the disagreement between the two columns' derived sea-surface heights,
which `VertCoord`'s midpoint rule produces even on the exact set. The anchor's shift to a common
pressure **never spans more than one layer** (§3.5.1); stretched across a multi-layer floor offset,
as terrain-following layers over a steep slope produce, the shift's first-order expansion left the
scheme at 0.45 of `PressureGradCentered`'s error, against 0.014 with the rule. The treatment at the
top and bottom of the column, assumption A5 of §3.7.6, is to **clamp to the outermost valid layer
and extrapolate its reconstruction**; measured, the scan stays flat to $1.4\times10^{-15}$ m even
where a column's reconstruction is evaluated three layers below its own floor, so the rule is
confirmed on the geometry it was most doubted on.

**Reopened by measurement.** Two things. First, the anchor end is a trade-off, not a settled
preference: at the sea floor a stepped bathymetry makes the anchor form and cancel quantities of
order the step within the deepest shared layer, which the surface anchor does not (§3.7.4).
Neither end gives machine precision at *every* tilt, thickness and bathymetry, so Requirement
2.3.1 is met where the two columns' floors share depth and pressure rather than unconditionally.
Second, on the seamount's sigma coordinate the scheme is unstable where `Centered` is not
(§5.3.1). It is not caused by the anchor, and its cause is unknown.

**Code and cost.** The `PressureGradFiniteVolume` functor; two kernels, `computeReconSlopes` and
`computeColumnScan`, with three per-edge and two per-cell working arrays (§4.1.3); two headers of
shared inline functions; no new TEOS-10 evaluations per cell and layer (Requirement 2.2), with
roughly three times as many column evaluations on a hexagonal mesh as a cell-based formulation
would need (§4.1.3). Most of the arithmetic is in the column scan, and §5.6 gives the measured
cost.

#### 4.5.2 Phase 2 — fourth order

**Delivers.** Parabolic vertical reconstruction (§3.4), widening the exact set to profiles that vary
quadratically with pressure; a wide horizontal stencil (§3.6.1); optionally a second-order
equation-of-state expansion (§3.3). Requirements 2.1 and 2.6 at fourth order.

**Must preserve.** Everything Phase 1 establishes. The machine-precision property is inherited by any
consistent edge-normal stencil (§3.6.1) rather than having to be re-established, so what Phase 2 must
not break is narrower than it appears: it must keep the expansion point shared across whatever cells
the stencil couples, and it must re-examine [](#z-increment-exact) before adopting the second-order
equation-of-state expansion, which would otherwise reopen the §3.7.4 question. The §5.2 gate is rerun
unchanged.

**No blocking open question.** The form of the wide stencil on Omega's unstructured,
variable-resolution TRiSK mesh is an ordinary accuracy-and-cost choice, settled by measurement against
Requirement 2.6, because exactness does not constrain it (§3.6.1).

One implementation question is new since Phase 1, though it does not block. The anchor
[](#anchor-one-layer) is defined for a pair of columns. A stencil built as a weighted sum of pairs
can reuse it, and the whole column scan, pair by pair. One built from a single
neighbourhood-wide expansion would need an anchor across several columns, whose floors may lie
several layers apart, and none exists. That favours the pair-sum form on implementation grounds as
well as on the assumption-A2 grounds §3.6.1 gives.

**Where Phase 2 attaches to the Phase 1 code.**

- *Configuration.* `VerticalReconstruction: 'ppm'` and `HorzOrder: 4` are already parsed and
  rejected in the `PressureGrad` constructor (§4.1.1); Phase 2 accepts them, and
  `PressureGradVertRecon` gains its `PPM` value.
- *Reconstruction.* Parabolic counterparts of the stencil, coefficient and evaluation functions in
  `PGradRecon.h`, exact for a quadratic profile on a non-uniform grid (§3.4), with per-cell
  coefficient arrays formed by `computeReconSlopes` alongside or in place of `ReconSlopeCt` and
  `ReconSlopeSa`. The column scan evaluates the reconstruction in three places — the integrand of
  each edge layer, the anchor's two short integrals, and the carry-down of
  [](#anchor-one-layer) — and all three must switch together, or the anchor will describe different
  water from the scan.
- *Quadrature.* With parabolic deviations and the first-order expansion, the integrand is quadratic
  within a sub-interval and its first moment cubic, so the two-point default still integrates both
  exactly there. A second-order expansion adds products of deviations, a quartic integrand and a
  quintic moment, which need three points. Exactness on the exact set is unaffected either way
  (§3.5, consequence 2).
- *Horizontal stencil.* The cell lists and weights are cached on the functor (§4.1.3).
- *Tests.* The reconstruction unit test gains a quadratic profile on non-uniform layers; the
  quadratic-in-pressure profile of §5.2 moves from the convergence group into the exact group; and
  the Polaris kernel that §5.1 compares against has to gain the same reconstruction, written from
  this document rather than from the C++, so that the comparison stays independent.
- *The sigma instability* (§5.3.1). Phase 2 should be run on the sigma seamount. A 700-fold smaller
  starting error grew at the same rate there, so higher accuracy alone is not expected to cure it.

#### 4.5.3 Suggested order of work

1. **Write the discrete form out in full, in the Polaris two-column harness, before writing C++.**
   Doing it there makes it executable — exactness becomes a pytest that runs in seconds on a login
   node — and doubles as the independent reference implementation §5.1 needs, which it would not be
   if written by reading the finished C++. *Substantially done*: the derivation of §3.1, the
   corrected [](#ho-target), the cancellation [](#grad-p-zero) and the numerical zero for a profile
   linear in pressure at large tilt have all been established this way. What remains is to re-express
   the confirmed result in the reduced form of §3.5, and to settle four things this design leaves to
   it:

   - the closed-form quadrature of [](#dalpha) and of the first moment [](#ho-exact) needs, together
     with the per-column pressure lookup that feeds them;
   - the anchor [](#anchor), including confirmation that it is zero for a resting ocean at the
     initialized sea-surface height;
   - which end of the column [](#d-recurrence) accumulates from (§3.7.4);
   - the rule used where a column's deepest reconstruction is evaluated below its own floor, and how
     much accuracy it costs off the exact set (assumption A5, §3.7.6).

2. Calibrate the seamount configuration with the centered scheme, which can be done immediately and
   in parallel with step 1. This establishes that the case exercises the pressure gradient at all —
   the coordinate must genuinely tilt, which a clipped p-star would prevent — and sizes the
   spurious velocity to expect. It is **not** the A4 diagnostic: §5.3 explains why that one cannot
   be run with `Centered` and has to wait for Phase 1.
3. Implement and verify Phase 1 against §5.1, §5.2, §5.3, and §5.5.
4. Take up Phase 2. Re-examine [](#z-increment-exact) before adopting the second-order
   equation-of-state expansion; that, and not the stencil, is the item that could reopen a settled
   question.

Steps 1–3 are done: step 1 in the Polaris `horiz_press_grad` kernel (polaris#675, with the anchor
rule in polaris#900), step 2 in the Polaris seamount tasks (polaris#681), and step 3 in Omega#493
and #635. Step 4 is next.

## 5 Verification and Testing

Testing reuses and extends the Polaris `horiz_press_grad` task family
(`polaris/tasks/ocean/horiz_press_grad`), which is already Omega-only and built around a
quasi-analytic, layer-mean TEOS-10 reference solution (`reference.py`, surface-anchored,
4-point Gauss). That analytic reference remains valid as "truth" for `FiniteVolume`,
so the principal changes are scheme selection and revised pass criteria, plus two full-model
acceptance tests.

The testing plan applies the two independent gates set out in §1: an **absolute error
tolerance at a representative coarse resolution** (Requirement 2.1), which establishes that
the scheme helps at the resolutions Omega can afford, and a **measured order of convergence**
(Requirement 2.6), which verifies that the implementation is the scheme this design
describes. Neither substitutes for the other, and both must pass.

Which tests gate which phase:

| Test | Phase 1 | Phase 2 |
|---|---|---|
| §5.1 Two-column convergence | gates, at second order | gates, at fourth order |
| §5.2 Machine-precision cancellation | gates | rerun unchanged; must still pass |
| §5.3 Seamount resting state | gates | rerun; error should drop further |
| §5.4 Overflow | gates | rerun |
| §5.5 Reduction to centered | gates | rerun unchanged |
| §5.6 Cost check | recorded | recorded; stencil width must not add EOS evaluations |

The A4 diagnostic within §5.3 gates nothing — it informs whether the work was prioritized correctly
rather than whether it is correct. It requires `FiniteVolume` and so can only be run *after* Phase 1,
not before it as earlier versions of this document assumed; §5.3 gives the reason.

### 5.1 Test: Two-column HPGA convergence (extend existing)

Extend the four existing variants — `temperature_gradient`, `salinity_gradient`,
`surface_pressure_gradient`, `ztilde_gradient` — to run the centered scheme alongside the new one:

- **Scheme selection.** *(Implemented and run; this bullet records what was built.)* Polaris adds a
  `PressureGrad` block to `forward.yaml` and parametrizes each task over the **two** configurations
  Phase 1 offers:
  - `centered` — the legacy `PressureGradCentered` functor, unchanged;
  - `finite_volume` — `HorzOrder: 2`, `VerticalReconstruction: linear`.

  An earlier version of this bullet listed a third, `finite_volume_phase2`, and described the last
  two as "the *same* implementation at two orders". That framing is retired: the Phase 2 values
  exist as keys but are **rejected with an error** (§4.1.1), so there is no third configuration to
  run until Phase 2 lands, and a configuration written for Phase 2 cannot quietly execute as
  Phase 1.

  Both schemes run over **one shared initial condition** — the Polaris `init` step is
  scheme-independent and writes both schemes' reference HPGA from the same state — so the
  comparison between them is made at an identical state rather than at two states that ought to
  match. That comparison is the clearest single measure of what Phase 1 buys. Running only the
  legacy functor would measure the convergence of code this design does not change and leave the
  new code unverified.

  The forward step still runs a single time step with only `PressureGradTendencyEnable: true`,
  reading the PGF acceleration from `NormalVelocityTend`.
- **Reference.** `reference.py`/`analysis.py` compare `NormalVelocityTend` against the
  layer-mean analytic HPGA (unchanged). For `FiniteVolume` the layer-mean comparison
  remains the correct target, since the scheme is a finite-volume, layer-averaged
  discretization.
- **Accuracy gate (new, Requirements 2.1 and 2.3.4):** at a representative coarse resolution (e.g.
  the coarsest in `horiz_resolutions`), the absolute RMS HPGA error vs. the reference must be below
  a tolerance, **and** the new scheme's RMS error must be below the centered RMS error at that same
  resolution. This gate applies to Phase 1, where it is the primary measure of value.

  **This is also the gate for Requirement 2.3.4.** Set the ratio per variant, at the coarse end of
  the sweep. **Measured in Omega** across the four variants, centered RMS over finite-volume RMS at
  4 km: $1.00\times$ (`temperature_gradient`), $1.00\times$ (`salinity_gradient`), $1.80\times$
  (`ztilde_gradient`), $2.19\times$ (`surface_pressure_gradient`).

  **The gate's default must therefore assert "not worse" rather than "better".** On
  `temperature_gradient` and `salinity_gradient` the two schemes agree to five significant figures,
  and that is correct rather than disappointing: those variants hold the coordinate level, so
  comparing at matched pressure and at matched layer index coincide, and the remaining error is the
  two-column representation of the tracer gradient, which both schemes share. A gate demanding an
  improvement would fail on a configuration where there is nothing to improve. Only the variants
  that tilt something raise it.

  Do **not** add a gate requiring the new scheme's convergence slope to exceed the centered
  scheme's. The reason is *not* that the two always converge alike — on `ztilde_gradient` they
  genuinely differ, and that is the headline result of Phase 1 (below). It is that they converge
  alike on three variants out of four, so such a gate fails where nothing is wrong.
- **Verification gate (Requirement 2.6):** the measured slope of RMS error vs. resolution,
  `omega_vs_reference_convergence_rate_*`, must fall within a band around the configured order
  of accuracy — nominally ~2 for `finite_volume` as Phase 1 ships it, and ~4 when Phase 2 lands. This band
  is retuned from its present values rather than loosened; a slope outside it fails the test and is
  treated as an implementation defect to be diagnosed, not as a tolerance to be widened.

  **The bands are *not* uniform across the four variants, and must be set per variant *and per
  scheme*, from measurement.** Measured in Omega on chrysalis with the bottom layer included:

  | variant | `centered` | `finite_volume` |
  | --- | --- | --- |
  | `temperature_gradient` | 1.631 | 1.632 |
  | `salinity_gradient` | 1.801 | 1.802 |
  | `surface_pressure_gradient` | 2.000 | 2.003 |
  | `ztilde_gradient` | **1.068** | **1.971** |

  A single "~2 for centered" band would fail three variants out of four, which is why the bands are
  per variant. Only `ztilde_gradient` needs a band of its own *per scheme*.

  **`ztilde_gradient` is the clearest single demonstration that Phase 1 works.** It holds $\Theta$
  and $S$ horizontally uniform and tilts only the coordinate, so the first-order term of
  [](#centered-error) is exposed with nothing else mixed in: the centered scheme is first order
  there (1.068) and the finite-volume scheme **restores second** (1.971). This is the resting-state
  behaviour of §3.7.3 and its cure, visible in a convergence rate rather than inferred.

  **One consequence for how the advantage is described.** Because the orders genuinely differ on
  this variant, its advantage *widens* under refinement rather than staying a constant factor:
  $1.80\times$ at 4 km against $14.8\times$ at 0.5 km. Elsewhere it is the constant factor this
  section previously assumed throughout — $2.19\times$ at both ends of the sweep on
  `surface_pressure_gradient`. Both statements are needed; neither generalizes.
- **Asymptotic range (Phase 2 implementation-time task):** it is not yet established that the
  existing `horiz_resolutions` sweep spans a range where a fourth-order slope is cleanly
  measurable — the sweep may be too coarse to have entered the asymptotic regime at its fine end,
  or fine enough that the reference solution's own quadrature error and roundoff contaminate the
  slope. Determining the usable window, and extending or tightening the sweep (and, if needed, the
  order of the Gauss quadrature in `reference.py`) so the designed order can be resolved, is part
  of implementing Phase 2. Phase 1's second-order slope is measurable on the existing sweep.
- **Consistency check (retained):** `omega_vs_polaris_rms_threshold` (~1e-10 m/s²) — Omega's
  forward output must still match the Python-computed HPGA, confirming the implementation
  matches the intended discretization.

  **Measured for the new scheme: $4\times10^{-19}$ to $5\times10^{-16}$ m s$^{-2}$** on the four
  gradient variants, against $1\times10^{-16}$ to $2\times10^{-15}$ for `centered` on the same
  states — six orders inside the threshold, and *better* than the legacy scheme's. Across all eight
  variants, including the resting-state ones of §6.6, it stays at or below $1.6\times10^{-13}$ on
  chrysalis and $5.4\times10^{-15}$ on pm-gpu. This is the strongest single piece of evidence
  that the scheme is implemented as designed, because the Python counterpart was written from this
  document rather than from the C++, so the two are independent realizations of the same
  specification.

  **That independence has to be maintained deliberately, and it nearly was not.** Reconciling the
  two implementations before the comparison was run found two places where they had silently
  diverged — the Python integrated with a 4-point Gauss rule where Omega defaults to 2
  (`QuadraturePoints`, §4.1.1), and it anchored the column scan at the sea surface where Omega
  anchors at the sea floor (§3.7.4). Neither shows up on the exact set, where the integrand is zero
  at every quadrature point and any rule integrates it to zero; both bite off it, which is exactly
  where this check operates. Had they not been reconciled first, the disagreement would have
  appeared as a single scalar per resolution with no indication which of nine discretization
  choices was responsible.

  The check does catch a single rule that differs. When the anchor rule of [](#anchor-one-layer)
  changed in Omega, the unchanged Python reference failed the comparison on both
  `hydrostatic_consistency` variants, at the 50 m/km tilts on 128 m and 64 m layers where the floors
  are more than a layer apart. With the rule mirrored in the reference it passes again, and that
  reference agrees with a separate one-edge model of the anchor to $2\times10^{-17}$ m s$^{-2}$ on
  terrain-following columns.
- **Assumption A2 is measured directly, not by a contrast sweep.** An earlier version of this
  section called for running `temperature_gradient` and `salinity_gradient` at several amplitudes of
  the horizontal contrast and requiring the error to grow no faster than linearly. **That test
  cannot work, for three independent reasons**, all established by measurement:

  - **The contrast axis is not independent of the resolution axis.** The two-column initialization
    interpolates at $x = \pm d_e/2$, so the state depends on the contrast $c$ and the spacing $d_e$
    only through their product. Verified exactly: contrast 2 at $d_e = 2$ km reproduces contrast 1
    at $d_e = 4$ km to every printed digit. A contrast sweep is therefore a *coarsening* sweep in
    disguise, which is why the convergence slope collapses from $1.46$ to $-0.04$ across it.
  - **A2 and the horizontal stencil's truncation error are confounded by construction.** Both scale
    with the cross-edge contrast, a comparison against the quasi-analytic reference cannot separate
    them, and the stencil term dominates by enough that `Centered` and `FiniteVolume` are
    indistinguishable — identical RMS to four significant figures at every contrast.
  - **The initialization's shape-preserving limiter is nonlinear in the node values**, so the state
    itself is not smooth in the contrast: the second difference of layer-mean salinity is
    $5.75\times10^{-2}$ g kg$^{-1}$ on the shipped profile against $7.1\times10^{-15}$ on a
    two-node profile where the limiter has no freedom. The nonlinearity saturates with contrast
    rather than scaling with it.

  **What replaces it is better.** A2 is the claim that the *second-order remainder* of
  [](#alpha-taylor) across an edge is small, and that remainder can be computed straight from the
  state without any reference comparison, so it is neither confounded with the stencil nor sensitive
  to how the state was built. Measured this way it is $8.0\times10^{-6}$ relative on
  `temperature_gradient` at 4 km — the hardest configuration in the family — with a companion check
  confirming it vanishes when the two columns coincide. **Pass:** the remainder stays below a stated
  fraction of $\alpha$ on the hardest configuration available. Should it ever fail, the second-order
  expansion of §3.3 becomes the default rather than an option.
- **Cfg keys.** New keys mirror the existing ones (`horiz_press_grad.cfg`): a coarse-resolution
  absolute tolerance, a `finite_volume_vs_centered` ratio gate, and per-scheme expected-rate bands.

**Covers:** Requirements 2.1, 2.2 (the bounded-EOS path is exercised), 2.4, 2.5, 2.6; the last row
of the §3.7.3 table (horizontal structure in $\Theta$, $S$); assumptions A1, A2, A3.

### 5.2 Unit test: Discrete hydrostatic consistency (exact resting state)

This is the primary test of the steep-layer robustness property (Requirement 2.3, §3.7). It follows
the table in §3.7.3 row by row rather than being a single pass/fail.

**Run the machine-precision gate on a state whose anchor inputs vanish identically.** The gate is
meant to test the scan, not the initialization, and the anchor is the one place where the state can
inject a residual the scheme is powerless to remove. At a sea-floor anchor that means the two
columns must share both a **floor depth** and a **bottom pressure**, so that $\Delta_e Z$ is zero and
the two short integrals of [](#anchor) are zero-width. At a surface anchor it means a horizontally
uniform surface pressure and equal sea-surface heights.

**This is a condition on how the state was built, and it does not come for free.** A state that
prescribes each column's floor depth and solves for the bottom pressure that reproduces it carries
`VertCoord`'s $O(\tilde h^2)$ truncation difference (§3.7.4) into the anchor's input, and then no
scheme however exact can return zero. Each column converges to the bottom pressure that makes *its
own* midpoint sum reach the prescribed depth, and under tilt the two columns' layer partitions, and
so their truncation errors, differ. Measured on Polaris-initialized states, which are built this
way, the sea-floor anchor is $-7.1\times10^{-8}$ m, equal to the difference of the two columns'
midpoint-rule truncations to four significant figures, where the surface anchor is zero to
round-off. That residual is the anchor and nothing else — every increment
of [](#d-recurrence) is still zero and $D_k$ is still flat down the column — so it does not
contradict the pointwise cancellation of [](#dalpha); it says the gate must be run on a state that
does not carry it.

Omega's own gate satisfies the condition by construction: it prescribes the pressure grid, giving
both columns the same surface and bottom pressure, and prescribes an equal `BottomGeomDepth` from
which `VertCoord` derives $z$ upward. Both anchor inputs are then identically zero and the measured
$3.1\times10^{-18}$ m s$^{-2}$ is a property of the scan alone. A run of the same scheme on a
Polaris-initialized state would be expected to return the anchor residual above rather than machine
zero, and that is the correct answer for that state. This costs nothing: surface-pressure gradients are exercised separately by §5.1,
where the reference solution supplies the expected value.

Where surface pressure *does* vary horizontally, the state carries a real fixed-pressure height
difference at the surface and the scheme is right to report it — **expect it, and assert it against a
computed value rather than against zero.** Two things make zero the wrong expectation there, and the
second is not removable:

- an initialization that balances the surface against a *reference* density $\rho_0$ rather than the
  in-situ density leaves a residual of relative size $(\alpha - 1/\rho_0)/\alpha$. In a pressure-based
  coordinate the in-situ balance is available in closed form — $z^{\text{surf}}$ is
  $-\tfrac1g\int\alpha\,dp$ over the surface pressure, a quadrature of a known integrand rather than
  the fixed-point solve a height-based Boussinesq model would need — so this part *can* be removed if
  a configuration wants to;
- [](#anchor) integrates $\hat\alpha^{(e)}$, the reconstruction, not the true $\alpha$. Exact
  cancellation at the surface would require the sea-surface height to be consistent with the
  *scheme's own linearized profile*, which is not something an initializer can reasonably be asked
  for. A residual at the equation-of-state linearization level therefore survives any physically
  correct initialization.

The second point is derived here and has not been measured; the estimate is that it is some two
orders below the first. Either way the conclusion for testing is the same, and it is why the first
paragraph exists: **do not put the exactness gate on a configuration with a varying surface
pressure.**

Set up a two-column (or seamount) configuration with the coordinate interfaces deliberately tilted
between the columns — including steep slopes and thin layers — and initialize the layer means as the
**exact layer averages** of a prescribed continuous profile $\Theta(p)$, $S(p)$. Under tilt those
averages come out *different in the two columns*, and that is the point: a configuration built by
copying identical layer means into both columns would not exercise the property at all (§3.7.1).
Three groups of profiles are run:

- **Profiles the scheme resolves exactly:** $\Theta$, $S$ linear in pressure, including the constant
  case, which isolates compressibility on its own. Run these with
  `VerticalReconstruction: 'linear'` (§4.1.1). **Pass:** the PGF is zero at every edge and layer to
  **machine precision** (double-precision builds; the threshold tracks `Real`'s epsilon and the size
  of the hydrostatic terms, not a physical tolerance). Phase 1 and Phase 2 must both pass.

  Report the tendency alongside `PressureGradCentered`'s on the same state, which is nonzero and
  tilt-dependent (§3.9). A pass is then a cancellation of ten orders of magnitude against a known
  nonzero reference rather than a comparison with zero, which is what makes it hard to pass by
  accident.

  Also assert the property one level down, at the integrand: [](#dalpha) must be zero **at every
  quadrature point**, not merely in the integral. This is the sharpest available form of the test,
  it localizes a failure to the layer and column that caused it, and it is the check that a
  per-layer-index sharing rule fails (§3.7.2, condition 1).

- **Profiles it does not:** $\Theta$, $S$ quadratic in pressure, then a realistic profile.
  **Pass:** the residual shrinks like $\tilde h^2$ (Phase 1) and $\tilde h^3$ (Phase 2) as the
  vertical grid is refined at fixed tilt, matching §3.7.3. A residual that does not shrink at the
  tabulated rate means one of the three conditions in §3.7.2 has been broken somewhere in the
  implementation; it is a bug to find, not a tolerance to widen.
- **Guard tests.** Rerun an exactly resolved profile with each of the following, and require each to
  *fail* the machine-precision check. Without them, a passing result could as easily come from a
  symmetry of the test setup, or from the scheme having become insensitive to tilt, as from the scheme
  being right.

  | guard | breaks | expected |
  |---|---|---|
  | (a) **tilt sensitivity** — assert the tendency and `PressureGradCentered`'s differ, and that the latter grows with tilt | nothing; runs on the passing configuration | must **pass**. A bug that zeroed the tilt response would satisfy every other check here perfectly |
  | (b) **cell-local expansion point** in place of [](#edge-ref) | condition 1 | fires. Measured in Omega at $3.4\times10^{-7}$ m s$^{-2}$ against $3.1\times10^{-18}$ with the shared point. Two expansion points mean the $\bar\alpha_0$ and $\bar\alpha_p$ terms no longer cancel in [](#dalpha) (§3.5, consequence 1) |
  | (c) **each column's $\Theta$, $S$ taken from its own layer $k$** rather than from the layer containing the pressure | condition 2 | **cannot be made to fire** — confirmed in Omega, which returns $3.4\times10^{-18}$ against $3.1\times10^{-18}$, indistinguishable. See the warning below |
  | (d) **anchor taken as the raw $\Delta_e Z$ alone**, dropping the short integrals of [](#anchor) | condition 3 | fires, but only where the two columns' end pressures differ; flat with depth, which distinguishes it from (b). Measured in Omega at $4.8\times10^{-3}$ m s$^{-2}$, with the offset in $D_k$ constant to sixteen digits down the column |

  Note that *which* shared coefficient set is used — selected by edge layer, by layer index, or
  otherwise — is deliberately **not** a guard: §3.5 consequence 1 says exactness cannot depend on it,
  and measurement confirms no effect (§6.5). A test asserting otherwise would be asserting something
  false.

  Guard (a) is the one to write first, because it is the only one that can fire on a configuration
  where every other check passes.

  **A warning that must not be lost, because it is the largest untested risk in the scheme.** Guard
  (c) **cannot be made to fire on any configuration in the Polaris family.** The exact-set variant's
  profile is a single line in pressure over the whole column, so every layer's mean-preserving
  reconstruction is that same line and looking up the wrong layer costs nothing; on the curved
  variant the two rules differ by less than a factor of two, which is not a usable discriminator
  either. **An implementation that looks up a column's state by layer index will therefore pass
  every exactness and accuracy check specified in this section.** The lookup must be pinned by
  direct property tests instead — that it returns a layer other than $k$ under tilt, and one whose
  interfaces actually bracket the pressure — and those tests are not optional. Whether an
  answer-level guard becomes available once the layer *integral* is formed is an open question.

  Guards must be checked against a deliberately broken configuration, not only against a passing
  one. A guard that cannot fire is worse than no guard, because it looks like protection.

  Guards (b), (c) and (d) need the scheme's rules switched off one at a time, which the shipped
  kernel cannot do. They therefore run against an assembly built in the test out of the same helper
  functions, with each rule selectable. That assembly is only trustworthy once it reproduces the
  kernel, and the fidelity check that has content is on a **curved** profile where both sides are
  large: agreeing at $10^{-18}$ on the exact set would also be satisfied by an assembly that
  computed nothing. Note what this does and does not establish — the guards show the three rules are
  load-bearing; they do not independently verify the arithmetic, since they share the helper
  functions with the kernel. That verification comes from §5.5 and from the Omega-versus-Polaris
  comparison of §5.1.
- **The anchor over a sloping floor.** The gate above deliberately gives both columns the same
  floor, so it cannot see the within-one-layer rule of [](#anchor-one-layer). A separate test does.
  It uses a resting state linear in pressure on 32 layers at 4 km spacing, whose floors are 8 MPa
  (6.4 layers) apart in pressure, a floor sloping at 200 m/km. Machine precision is not the
  expectation there, since the anchor carries the difference between the two columns' `VertCoord`
  sums (§3.5.1), so the gate is a ratio to `PressureGradCentered` on the same state. Measured in
  Omega, the rule gives 0.014 of `Centered`, against a gate of 0.05 (0.15 in single precision, where
  floor heights near 4000 m carry round-off of order $10^{-4}$ m). Its guard turns the rule off in
  the host assembly, shifting both floors to their mean pressure as [](#anchor) alone would. That
  gives 0.45 of `Centered` and must fail the gate by a wide margin, which it does. The host
  assembly's fidelity to the kernel is checked on the same state, to round-off of the height
  difference the anchor forms and cancels.

A separate and much smaller unit test covers the reconstruction estimator on its own (§3.4): given
layer means sampled from a profile of the reconstruction's own degree on a **deliberately
non-uniform** set of layer thicknesses, the recovered slope (Phase 1) or slope and curvature
(Phase 2) must match the exact values to round-off. This is worth testing separately because it is
the most likely place for the machine-precision gate above to fail, and it localizes the failure
immediately.

The test is also run in a single-precision build. §3.7.5 predicts it should *pass*, because
[](#dz-dp) forms no large quantities; running `PressureGradCentered` alongside, which does, is what
makes the result interpretable.

It is implemented as a fast C++ unit test and also run as a configuration of the Polaris
two-column task. In Omega it is `PGradTest.cpp`, registered as `PGRAD_TEST` and, built in single
precision, `PGRAD_SINGLE_PRECISION_TEST`. Its parts are the reconstruction estimator
(`testReconstruction`), the property tests that pin the pressure lookup (`testPressureLookup`), the
pointwise integrand (`testMatchedPressIntegrand`), the gate, convergence and guards above
(`testExactnessAndGuards`), the sloping-floor anchor (`testSlopingFloorAnchor`), the centered
identity of §5.5 (`testCenteredIdentity`) and the option parsing of §4.1.1 (`testPGradConfig`).

**Covers:** Requirement 2.3 in full; rows 1–5 of the §3.7.3 table; the round-off floor of §3.7.5;
the anchor rule of §3.5.1.

### 5.3 Test: Seamount resting state (steep-layer robustness)

Use the Polaris seamount task with layers tilted over a seamount and a stratification that is a
function of pressure alone — note that this gives layer means that *differ* between neighboring
columns wherever the layers tilt (§3.7.1), which is the situation the test is meant to create.
Integrate for a fixed period. **Pass:** the maximum spurious velocity stays below a threshold and
is substantially smaller for the new scheme than for the centered scheme. This is the dynamical
counterpart of §5.2 — the same resting state, now run forward in the full solver — and the direct
test of Requirement 2.3.

**Testing assumption A4 (§3.7.6).** Run the seamount case twice: once with a stratification the
scheme resolves exactly ($\Theta$, $S$ linear in pressure) and once with a realistic profile.
Spurious velocity that survives the first run cannot be PGF truncation error, since the PGF is zero
there to machine precision, and must come from somewhere else in the model — most likely the
layer-mean treatment in the tracer and remapping operators (§4.3).

**The diagnostic must be run with `FiniteVolume`, not with `PressureGradCentered`.** Earlier
versions of this section, of §4.5.3 step 2 and of the gate table above described it as a
*pre*-Phase-1 check that could be done immediately with the existing scheme. That was an error, and
§3.7.3's own table is what contradicts it: a profile linear in pressure is in Phase 1's exact set
and **not** in the centered scheme's, where it remains $O(\tilde h)$ with a coefficient
proportional to the coordinate tilt. Run with `Centered`, the linear-in-pressure seamount is not a
null configuration but merely a slightly cleaner non-null one — and a sigma coordinate over a
seamount is the largest tilt in the suite, so the contamination is at its worst exactly there. The
diagnostic would measure the thing it was designed to exclude, and its result would support neither
"A4 holds" nor "A4 fails".

The consequence for ordering is the opposite of what §4.5.3 assumed: A4 is a **post**-Phase-1 test,
and it became runnable for the first time when `FiniteVolume` landed. What can be run beforehand is
a *calibration* of the centered scheme's response to tilt, which is useful but answers a different
question; §5.3.2 is the measurement that does bear on A4 without needing the seamount at all.

#### 5.3.1 Measured: the scheme passes on near-level layers and is unstable on sigma

Run on the Polaris seamount for six days, with all explicit mixing off so that the pressure
gradient is what remains, under a linear equation of state and TEOS-10, each with an exponential
stratification and one linear in pressure. Maximum spurious velocity at day 6, in m s$^{-1}$, with
the day a run first overturned in parentheses:

| equation of state, stratification, coordinate | `Centered` | `FiniteVolume` |
| --- | --- | --- |
| linear, exponential, z-star | $8.4\times10^{-3}$ | $5.0\times10^{-4}$ |
| linear, linear in pressure, z-star | $1.3\times10^{-2}$ | $1.8\times10^{-10}$ |
| TEOS-10, exponential, z-star | $0.22$ (day 0.5) | $3.2\times10^{-4}$ |
| TEOS-10, linear in pressure, z-star | $0.14$ (day 1.9) | $1.3\times10^{-6}$ |
| linear, exponential, sigma | $1.1\times10^{-2}$ | $0.54$ (day 1.6) |
| linear, linear in pressure, sigma | $3.4\times10^{-5}$ | $0.67$ (day 4.0) |
| TEOS-10, exponential, sigma | $1.4\times10^{-2}$ | $1.4$ (day 1.4) |
| TEOS-10, linear in pressure, sigma | $4.1\times10^{-3}$ | $0.94$ (day 2.3) |

**On z-star, `FiniteVolume` passes §5.3's criterion by a wide margin**, and under TEOS-10 it stays
stable where `Centered` overturns. Extended to 20 days on the TEOS-10 trees it stays clean, with
kinetic energy near $10^{-10}$ m$^2$ s$^{-2}$ on the exponential profile and $3\times10^{-14}$ on
the linear-in-pressure one, while `Centered` overturns within two days. This is Requirement 2.3
under dynamics.

**On sigma it does not.** `FiniteVolume` overturns on all four sigma trees, between days 1.4 and 4,
where `Centered` on the identical configuration stays stable.

It is **not an accuracy failure**. At one hour the new scheme is still 1.6 to 10 times *more*
accurate than the centered one, so the problem is growth over time. Three things point at a
grid-scale numerical instability rather than a resolved response: the disturbance sits at a fixed
location on the flank, the same under both equations of state, so it follows the geometry rather
than the profile; its grid-scale content rises before any inversion appears; and the temperature
anomaly reaches a third of the profile's whole range.

**What has been ruled out.**

- *The sea-floor anchor.* The within-one-layer rule of §3.5.1 starts the linear-in-pressure sigma
  run about 260 times lower in kinetic energy, and an experimental anchor about 700 times lower.
  In both, the same mode then grows at the same rate, about 2.6 per day, to the same state. Also
  turning off extrapolation below a column's floor (§3.7.6, A5) makes it worse.
- *Unlimited tracer advection.* The monotonic tracer limiter was already on in these runs, which
  refutes the working hypothesis recorded here earlier, that an unlimited horizontal tracer
  advection admits a mode the more accurate pressure gradient no longer damps.
- *Vertical resolution.* At full seamount height 32, 64 and 128 layers grow at the same rate, and
  the 64- and 128-layer runs go NaN once the flow saturates. At half height it is delayed, not
  cured: the TEOS-10 linear-in-pressure tree overturns at day 18 with 32 layers and day 7 with 64.
  The instability tracks the seamount's slope, on the steep flank, not the layer count.

Because a 700-fold smaller starting error grows at the same rate, this looks like the scheme's
response on steep sigma layers rather than its accuracy, so Phase 2's higher order is not expected
to cure it by accuracy alone (§4.5.2).

**What this does and does not mean.** Sigma layers follow the bathymetry, which is the steepest
tilt against isobars in the whole suite and is not what Omega runs — p-star, like z-star, keeps
layers near level. So this does not qualify the Phase 1 result on the configurations Omega uses. It
does mean the scheme has a limit on steeply tilted layers that the centered scheme does not share,
that the limit is not understood, and that anything relying on sigma-like layering, including the
coordinate under an ice shelf, is untested territory.

#### 5.3.2 What a global resting state already shows, and what it does not

A measurement bearing on A4 exists and is recorded here so the assumption is argued from numbers
rather than from plausibility. **On its own it does not settle A4**; the end of this subsection says
what it leaves open, and what has since been measured with Phase 1 in hand.

**What was measured.** The centered scheme's PGF was evaluated **offline** on a realistic global
initial condition — `EC30to60E2r2` with horizontally uniform $\Theta$ and $S$
(`ocean.EC30to60E2r2.200908.teos10.uniformTS.nc`, by `measure_hpg_error.py`, both in the Polaris run
directory `test_20260727/hpg_uniform_ts/`). Uniform $\Theta$, $S$ as functions of pressure make this
an exact resting state on the real mesh and bathymetry, so **the entire computed PGF is error** —
§5.2's premise at global scale. Split into depth-mean (barotropic) and deviation (baroclinic), in
m s⁻²:

| | median | p99 |
| --- | --- | --- |
| all active edge-layers, total | $9.30\times10^{-7}$ | $1.13\times10^{-5}$ |
| all active edge-layers, baroclinic | $9.45\times10^{-8}$ | $1.81\times10^{-6}$ |
| deepest common layer, total | $2.70\times10^{-6}$ | $1.62\times10^{-5}$ |
| deepest common layer, baroclinic | $1.81\times10^{-6}$ | $1.40\times10^{-5}$ |

**Three things this establishes.**

1. **The error is concentrated in the bottom layer and is baroclinic there.** The interior error is
   ~90% barotropic, which a free-surface adjustment removes within hours and which therefore cannot
   sustain anything; the bottom layer is 67% baroclinic. That matches the bottom-trapped, 97–98%
   baroclinic character of the KE response that motivated this work.
2. **It is large enough to matter.** In the bottom layer the truncation error is $1.28\times$
   (median) and $1.55\times$ (p99) the pressure gradient the model actually carries there — and the
   carried gradient itself contains the error, so the ratio to the true signal is worse. **This sets
   the bar for Requirement 2.1**: a useful improvement is a factor of 5–10 at the slopes the global
   mesh actually has, not the ~2× a modest gain would give.
3. **Interior coordinate tilt is not the mechanism.** The global IC's interior coordinate is
   level to 2–18 Pa cross-edge, $2$–$5\times10^{-5}$ of a layer, while the smallest tilt in the
   Polaris `hydrostatic_consistency` sweep is $1.74\times10^{3}$ Pa — 100 to 1000 times more.
   Scaling the
   measured $q = 1$ line (§3.7.3) down to the global tilt gives $\sim1.6\times10^{-10}$ m s⁻², about
   6000× below the global interior error. **Whatever drives the bottom-layer error, it is not the
   tilting of interior coordinate surfaces** — which is why §5.2's tilt scan calibrates the scheme
   without predicting a global improvement.

**And one thing it calibrates.** Comparing like with like — baroclinic bottom-layer error, with the
Polaris `bathymetry_step` value scaled by $(199/256)^2$ for layer thickness — the two-column
configuration reproduces the global figure to within 20–35% across the p25–p75 slope range where
most of the mesh sits ($1.20\times$ at 1.7 m/km, $1.34\times$ at 4.5 m/km). It over-predicts steeply
sloped edges by 3–10× ($2.9\times$ at 11 m/km, $9.9\times$ at 30 m/km), presumably because an
isolated two-column step is harsher than a real slope where neighbouring cells share the descent. So
`bathymetry_step` is a quantitative proxy for the global bottom-layer error over the range that
dominates the mesh, and its tail should be read with that in mind.

**What this does not settle, and it is the whole of A4.** Every number above is the *error*, not its
dynamical consequence. That the PGF error exceeds the gradient the model carries is necessary for
A4 but not sufficient: it does not show this error is what *produces* the observed bottom-layer KE
rather than merely coexisting with it. The one forward-run observation available before Phase 1 —
the uniform-T/S state gives bottom KE $1.19\times$ that of the fully forced run — is suggestive and
no more.

**Measured since, with Phase 1 in hand.** Two forward-run results bear on A4 directly.

- *The realistic global run itself.* Polaris `realistic_global/EC30to60E2r2` was run for 10 days
  on chrysalis, with the task's defaults otherwise, once with each scheme. `FiniteVolume` removes
  the grid-scale noise in bottom-layer kinetic energy that `Centered` shows across the open ocean.
  What remains is mostly on shelves and in marginal seas. Bottom-layer KE, in m² s⁻²:

  | scheme | area mean, day 1 | day 5 | day 10 | 99th percentile, day 10 |
  | --- | --- | --- | --- | --- |
  | `Centered` | $1.30\times10^{-2}$ | $8.01\times10^{-3}$ | $7.44\times10^{-3}$ | $9.10\times10^{-2}$ |
  | `FiniteVolume` | $4.60\times10^{-4}$ | $4.86\times10^{-4}$ | $4.52\times10^{-4}$ | $7.88\times10^{-3}$ |

- *The seamount diagnostic of §5.3*, run with `FiniteVolume` as it must be. On z-star the
  linear-in-pressure run stays at kinetic energy near $3\times10^{-14}$ m² s⁻² for 20 days (§5.3.1),
  so nothing else in the model generates spurious flow there.

Together these support A4 over the open ocean: the PGF error was driving the bottom-layer kinetic
energy there, and Phase 1 removes it. They do not establish A4 on shelves and in marginal seas,
where some bottom-layer kinetic energy remains under `FiniteVolume`. The sigma seamount does not
bear on A4 either way, since there the new scheme's own instability dominates (§5.3.1).

**Covers:** Requirement 2.3 under dynamics; assumptions A3, A4.

### 5.4 Test: Overflow (full non-Boussinesq dynamics)

Use the Polaris overflow task to exercise the PGF within the full non-Boussinesq equations with
strongly sloping layers and active dynamics. Unlike the seamount task, the state here is neither
resting nor close to a profile the scheme resolves exactly, so this is the test of how the scheme
behaves where none of the exact-cancellation results apply and the residual of §2.3.4 is all there
is. **Pass:** the solution remains stable and the down-slope evolution agrees with the reference
behavior; spurious mixing and velocity attributable to PGF error are reduced relative to the
centered scheme.

This is also the most demanding test of assumption A2: the overflow front carries a large horizontal
temperature contrast across edges with steeply sloping layers, which is exactly where the shared
edge expansion point [](#edge-ref) is worked hardest. If the direct measurement of A2 in §5.1
indicates the second-order equation-of-state expansion is needed, this task is where that judgment
is confirmed under dynamics.

**So far only a smoke test has been run.** The 12-minute p-star overflow smoke test runs with
`FiniteVolume` and differs from `Centered` by $0.055$ m s$^{-1}$, against a maximum $|u|$ of about
$1.4$ m s$^{-1}$, and the overflow tasks in Polaris' `omega_pr` suite run to completion with it. The
comparison this section specifies is still outstanding.

**Covers:** Requirements 2.1 and 2.3 under realistic, coupled conditions; the last two rows of the
§3.7.3 table; assumptions A2, A3.

### 5.5 Test: Reduction to the centered scheme (permanent regression)

Assert [](#centered-identity) against `PressureGradCentered` directly, on the two-column test, across
a sweep of tilts. This is a unit test of the *identity* of §3.9, not of a reduced configuration of the
new scheme: [](#dz-dp) differences the integrand and never forms [](#centered-shift), so there is no
setting that recovers the centered scheme (§3.9, final paragraph).

It guards Requirement 2.5 and protects the centered option during refactoring, and it is cheap —
[](#centered-shift) is a few lines and needs none of the new machinery. Because the two functors read
the mesh, `VertCoord` and EOS state through independently written code, their agreement tests the
shared upstream state — edge masks, interface indexing, `VertCoord` conventions — and not just the PGF
arithmetic (§4.4). Measured, it holds to $0.5\,\epsilon$ of the hydrostatic scale over fifty states
(§6.2).

Expect agreement to round-off rather than bit-for-bit. In addition to the operation-order argument of
§4.4, `PressureMid` may be formed as $\tilde z_{\text{bot}} + \tfrac12\tilde h$ rather than as the mean
of the two interface pressures; those are algebraically identical but not bit-for-bit.

A second, weaker check is worth running alongside: the `FiniteVolume` and `Centered` tendencies must
agree to $O((\Delta_e p)^2)$ as the tilt is reduced. That is the only sense in which the two schemes
now reduce to one another, and a violation indicates a defect in one of them.

**Covers:** Requirement 2.5; §3.9.

### 5.6 Cost check

Requirement 2.2 bounds the number of TEOS-10 evaluations, and nothing above tests it — the scheme
could satisfy every accuracy gate while quietly calling the equation of state inside the column loop.
Two things guard against that:

- **Interface.** The pressure gradient performs **no** equation-of-state evaluations: the one per
  cell per layer the requirement allows is paid once by `AuxiliaryState`, before the tendency is
  formed, and the scheme works from the specific volume and derivatives that call leaves behind.
  `computePressureGrad` receives `Eos` as a `const Eos *`, so the compiler rejects a call to any of
  its compute methods. A run-time counter in those methods would add nothing: it could only see the
  host-side calls the `const` already forbids, while an evaluation inside the quadrature loop would
  be a device-side call to the TEOS-10 polynomial that no such counter sees. That remaining case is
  left to review and to the wall-time check below.
- **Wall time.** Record PGF kernel time relative to `PressureGradCentered` on a representative
  configuration, as a performance regression guard. The expected cost is dominated by the per-edge
  column scan (§4.1.3), roughly three times as many column evaluations as a cell-based formulation on
  a hexagonal mesh; a result far above that suggests the per-cell coefficients are being recomputed
  per edge rather than cached.

  **Measured.** In the 10-day `EC30to60E2r2` runs of §5.3.2 on chrysalis (Intel, 6 nodes), the
  pressure-gradient tendency took 9.0 s with `FiniteVolume` against 1.0 s with `Centered`, out of
  under 3 minutes per run. That was before the column kernels became team loops over each edge's
  active layers, with the two integrals fused into the scan kernel (§4.1.3). That change cut the
  `FiniteVolume` time by about a third on chrysalis (seamount 12.8 s to 8.1 s), and by 9–30% on the
  overflow and `QU240` cases on pm-gpu, Frontier and Aurora, CPU and GPU, with bit-for-bit answers
  throughout.

**Covers:** Requirement 2.2.

### 5.7 Coverage summary

| Requirement / assumption | Verified by |
|---|---|
| 2.1 Accuracy at affordable resolution | §5.1 accuracy gate; §5.3; §5.4 |
| 2.2 Bounded TEOS-10 cost | §5.6 |
| 2.3 Robustness for thin, sloped layers | §5.2 (machine precision); §5.3 (under dynamics) |
| 2.4 Consistency with model state | §5.1 (layer-mean reference); §5.2 (reconstruction unit test) |
| 2.5 Runtime-selectable, backward compatible | §5.1 (two configurations); §5.2 option parsing; §5.5 |
| 2.6 Verified order of accuracy | §5.1 verification gate |
| A1 Edge accuracy ≠ cell accuracy | §5.1 verification gate (Phase 2) |
| A2 EOS expansion adequate across an edge | §5.1, by direct measurement of the second-order remainder across the edge; §5.4 under dynamics |
| A3 Residual small enough in practice | §5.1 accuracy gate; §5.3; §5.4 |
| A4 PGF error causes the instability | §5.3 diagnostic, run with `FiniteVolume`; the realistic global run of §5.3.2 |
| A5 Top- and bottom-of-column treatment | §5.1 `bathymetry_step` and `surface_pressure_gradient` |
| §3.7.2 condition 1, shared at matched pressure | §5.2 guard (b), and the per-quadrature-point assertion |
| §3.7.2 condition 2, differenced before integration | §5.2 guard (c), which cannot fire; the pressure lookup's property tests |
| §3.7.2 condition 3, the anchor | §5.2 guard (d) |
| §3.5.1 the anchor's shift within one layer | §5.2 sloping-floor anchor test and its guard |
| §3.5 the scheme responds to tilt at all | §5.2 guard (a); §5.5 tilt sweep |
| §3.9 the centered identity | §5.5, as a standalone unit test |
| §3.7.5 Round-off floor | §5.2, run in both precisions, alongside `Centered` |

Requirement 2.7 (extensibility) is not testable directly; it is addressed by the configuration
design of §4.1.1 and by the phase structure of §4.5, which is itself the evidence that the framework
extends without restructuring.

## 6 Supporting evidence

Two things §3 asserts are not self-evident from the algebra alone: that the cancellation of §3.1.1
survives discretization, and that the form of [](#ho-target) used here is the correct one. This
section supplies the evidence for both.

### 6.1 Numerical confirmation of the reduction and the exact set

The claims below were established in a purpose-built two-column Python harness before any C++ was
written, and they are the reason §3.2 takes [](#ho-exact) as its target rather than assembling the
four pressure terms. The harness builds its state on the Polaris `horiz_press_grad` initialization:
interface pressures prescribed directly, $\Theta$ and $S$ sampled as `Init._interpolate_t_s` does,
equation-of-state coefficients from one `gsw.specvol_first_derivatives` call per cell per layer, $z$
accumulated as `VertCoord` does, geometry from `hydrostatic_consistency.cfg`
($d_e = 4$ km, $\tilde z_{\text{bot}} = -4096$ m). Where a truth value is needed it is
$-g\,\partial z/\partial x|_p$ by direct segment-wise Gauss integration of the true TEOS-10
$\alpha(p)$, sharing no code with the scheme under test.

| claim | §ref | measurement |
|---|---|---|
| the four pressure terms of [](#ho-target) sum to zero | §3.1.1 | $1.4$–$5.0\times10^{-15}$ m s$^{-2}$ against a side-wall term of $0.55$; relative $3\times10^{-15}$, at 50 m/km tilt, for profiles linear in $p$, quadratic in $p$, and with a horizontal temperature contrast |
| the $p\nabla\alpha$ term of [](#ho-target) is required | §3.1 | with it, error $1.45\times10^{-9} \to 2.26\times10^{-11}$ at rate $2.00$ per refinement; without it, error pinned at $4.55\times10^{-6}$ at rate $0.00$ — the scheme returns $9.0\times10^{-7}$ against a true $5.42\times10^{-6}$ |
| exactness on the exact set, independent of tilt | §3.7.3 | residual at the round-off floor, $9\times10^{-16}$ to $9\times10^{-15}$, **not growing** across 0.5–200 m/km at 16 and 64 layers. `Centered` on the same states is exactly first order in tilt |
| $O(\tilde h^2)$ off the exact set | §3.7.3 | quadratic-in-$p$ profile at 50 m/km: rates 2.21, 2.14, 2.09, 2.05 under vertical refinement |
| the edge-shared expansion point is load-bearing | §3.3.1 | replacing it with a cell-local one costs $2.2\times10^{-9}$ to $3.0\times10^{-5}$ across the tilt range — §5.2 guard (b) |
| the fixed-pressure comparison is the whole scheme, not a correction | §3.5 | omitting it costs $2.8\times10^{-3}$ to $1.1$ m s$^{-2}$ — three to five orders larger than any other guard |

Three limits on what this harness established. The first two have since been closed by the
Polaris task runs of §6.5 and §6.6.

- **Differing `maxLevelCell` was untested.** The harness gave both columns the same layer count, so
  the `bathymetry_step` case — where a bottom partial cell carries the whole signal — was not covered.
  This is assumption A5 of §3.7.6, since measured as benign.
- **The Polaris task machinery was untested against the new scheme.** The `horiz_press_grad` variants
  now run both schemes, so the interaction with p-star initialization and partial cells, and the
  `omega_vs_polaris` comparison for `FiniteVolume`, are measured (§5.1, §6.6).
- **The rows above were measured on an assembly of all four pressure terms**, not on the form §3.5
  specifies. Given that those terms cancel to $3\times10^{-15}$ relative, what the exactness and
  $O(\tilde h^2)$ rows measured *was* the geopotential term alone, so they carry over — but by that
  argument rather than by direct measurement. §6.3 records what happened when it was measured
  directly.

#### 6.1.1 Confirmation in Omega

The implementation was measured independently of the harness above, on a two-column state whose layer
means are the exact layer averages of one prescribed continuous profile, at 4 km spacing with 32
layers and interfaces offset by up to three layer thicknesses.

| claim | §ref | measurement |
|---|---|---|
| exactness on the exact set | §2.3.1, §3.7.3 | $3.1\times10^{-18}$ m s$^{-2}$, against `PressureGradCentered`'s $3.5\times10^{-5}$ on the same state — a ratio of $8.8\times10^{-14}$, and $0.05\,\epsilon$ of the terms that had to cancel |
| $D_k$ zero at every interface, not only in the layer mean | §3.5.1 | $1.4\times10^{-15}$ m, $0.06\,\epsilon$; larger in the upper column than the lower, i.e. growing away from the sea-floor anchor, which is the signature of an exact anchor plus an accumulating recurrence |
| the integrand is zero *pointwise* | §3.5, §3.7.2 | $5.8\times10^{-21}$ m$^3$ kg$^{-1}$ at every quadrature point on a profile linear in pressure, against $5.3\times10^{-10}$ on a quadratic one — the Polaris harness measured $5.3\times10^{-21}$ on the same check |
| $O(\tilde h^2)$ off the exact set | §3.7.3 | rates above 2 on quadratic and cubic profiles over 15, 30 and 60 layers at fixed tilt. The rates *rise* under refinement rather than settling, so the sweep is not yet in the asymptotic regime and the gate is a lower bound |
| the centered identity | §3.9, §5.5 | $0.49\,\epsilon$ of the hydrostatic scale over eight tilts, reproducing the harness' $0.5\,\epsilon$ through independently written code. The dropped-pressure-term guard is exactly linear in tilt, to six significant figures — §3.7.3's first-order-in-tilt statement falling out of a check written for another purpose |
| bounded equation-of-state cost | §2.2, §5.6 | by construction: the pressure gradient holds a `const Eos *` and reads only the four cached fields |
| single precision | §3.7.5 | the exactness gate passes; see §3.7.5 for the numbers |
| the anchor's shift kept within one layer | §3.5.1, §5.2 | on a separate state with floors 6.4 layers apart, profile linear in pressure: 0.014 of `PressureGradCentered`, against 0.45 with both floors shifted to their mean pressure |

### 6.2 Corrections owed to `OmegaV1GoverningEqns.md`

[](#ho-target) is inherited from {ref}`omega-design-governing-eqns-omega1`, which is not maintained
against the code. Three corrections are owed to it, all following from [](#leibniz), and they should
be made there rather than left as local overrides here.

1. **The $+\int_A\int\rho_0\,p\,\nabla\alpha\,d\tilde z\,dA$ term is missing.** Row 2 of §6.1 is what
   it costs: a scheme built on the four-term form without it does not converge at all. The omission is
   invisible in every resting-state configuration, because $\partial_x\alpha|_{\tilde z} = 0$ when
   $\Theta$ and $S$ are horizontally uniform, which is why it survived review.
2. **The signs on the two metric terms are reversed.** [](#leibniz) gives
   $+(\alpha p)_{\text{top}}\nabla\tilde z^{\text{top}}$ and
   $-(\alpha p)_{\text{bot}}\nabla\tilde z^{\text{bot}}$. With the signs as printed there, the
   cancellation of §3.1.1 doubles rather than vanishing.
3. **The small-slope interface normals are not needed.** Its §6 reaches [](#ho-target) through
   $\mathbf n^{\text{top}} \approx (-\nabla\tilde z^{\text{top}}, 1)$, unnormalized. The derivation in
   §3.1 here uses Leibniz' rule directly, so neither the corrected [](#ho-target) nor [](#ho-exact)
   depends on that approximation.

### 6.3 Measurements from the Polaris reference implementation

A second round of measurement, on the Polaris two-column harness at
`polaris/tasks/ocean/horiz_press_grad` and asserted by `pytest`, confirmed two claims of §3 and
falsified a third. It is the reason §3.5 differences the integrand rather than correcting
[](#centered-shift).

**Confirmed — the centered identity (§3.9).** [](#centered-shift) reproduces the Polaris centered
HPGA to at most $0.5\,\epsilon$ of the hydrostatic scale, over fifty states: both resting-state tilt
sweeps in full, and the four gradient variants at their coarsest and finest resolution. It holds in
the bottom partial cell and at tilts where the two columns' `maxLevelCell` differ. The check is not
vacuous — plausible mis-derivations miss by $700\times$ the tolerance (cell-local instead of
edge-averaged $\alpha$), $1.2\times10^{6}$ (top interface instead of the trapezoid) and
$2\times10^{9}$ (dropping the pressure term).

**Confirmed — accumulating small quantities recovers precision (§3.7.5).** The largest increment of
a cancellation-free accumulation, against the fixed-index height difference it replaces:
$1.1\times10^{-3}$ at 256 m layers and 50 m/km, $2.7\times10^{-4}$ at 64 m, $2.0\times10^{-4}$ on a
200 m/km bathymetry step — about three decimal digits of headroom on a tilted coordinate. The
gradient variants gain little only because their coordinate is flat and there is nothing to save.

**Falsified — building each column's height integral separately and reconciling the two.** On a
configuration inside the exact set, that formulation missed exactness by
$\approx 2\times10^{-3}$ of the centered scheme's answer, **flat across three decades of tilt**,
where the robustness property asks for round-off. A brute-force evaluation of [](#ho-exact) by direct
quadrature — sharing no code with the closed forms — missed by the same amount, which is what
identified the target rather than the algebra as the problem, and prompted [](#dz-dp).

**The mechanism of that failure was never established.** Two candidate explanations were proposed
and both were later excluded by measurement: sharing the equation-of-state expansion at fixed layer
index rather than at fixed pressure (no effect at all, §6.5), and extrapolating a reconstruction
beyond its own layer (benign, §6.5). The formulation is abandoned, so the question is left open
rather than pursued. What matters for the design is the positive result in §6.5: differencing the
integrand is exact pointwise, which the failed formulation never was.

Worth noting for scale: at 64 m layers and 50 m/km the two columns' layer $k$ are offset by 2.8 layer
thicknesses and **do not overlap in pressure at all**. Any formulation that pairs by layer index is
differencing unrelated water there.

### 6.4 Corrections carried over from the reference implementation

Five items, all confirmed numerically on the Polaris side and folded into §3–§5 above.

1. **[](#centered-error) reconciles only in the redistribution limit.** Summation by parts on a
   cancellation-free accumulation gives [](#centered-error) plus a $k$-independent anchor constant
   **only when $\sum_j \Delta_e\tilde h_j = 0$** — the two columns having equal total pseudo-thickness.
   A `z_tilde_bot` tilt gives them different pseudo-bottom depths, so this does not hold in
   `hydrostatic_consistency`, the configuration the design cites. [](#centered-error) should be read
   as the redistribution-limit statement its own derivation assumes.
2. **Sources of a correction are not individually small.** An earlier formulation required each
   component of its remainder to be "a small quantity". Measured, the largest was exactly the size of
   the term it corrected — necessarily, since it is the term that cancels it. The requirement that
   matters is that no component be *formed by subtracting large quantities*, which is a statement
   about how each is evaluated rather than about its magnitude. §3.5.1 states it that way.
3. **Two exact simplifications** hold when the edge control volume is the average of the two columns'
   interface pressures: the edge layer's mid-pressure is exactly the edge average of the two columns'
   own layer midpoints, and its thickness exactly the edge average of theirs. Both follow from
   `PressureMid` being the exact arithmetic midpoint. §3.7.4 records them.
4. **`gsw` units.** `gsw.specvol_first_derivatives(SA, CT, p)` takes pressure in **dbar** and returns
   the pressure derivative of specific volume per **Pa** — the convention of neither argument. The
   Omega unit test of §4.1.2 should assert this against a centred finite difference in Pa, and
   explicitly against being wrong by the factor $10^{4}$.
5. **Initializing $\Theta$ and $S$ as layer means.** Omega's prognostic $\Theta$, $S$ are layer means
   and §3.4 reconstructs them as mean-preserving, so a test initial condition should supply layer
   means rather than midpoint samples. Changing this moves the Polaris resting-state baselines by
   3–39% (and the gradient variants by $\le 0.02\%$), so the recorded `Centered` baseline and the
   Omega-vs-Polaris comparison both move with it. This is a property of the test configuration, not
   of the scheme, but it invalidates any comparison against the earlier baseline.

### 6.5 Measurements of the scheme as specified

The scheme of §3.5 was implemented in the Polaris two-column harness and measured. These are the
numbers Omega should expect to reproduce, and the two places where they contradict this document.

**Exactness (Requirement 2.3.1) — met.** On a resting profile exactly linear in pressure, across the
full tilt sweep at 256/128/64 m layers: [](#dalpha) is zero at every quadrature point to
$5.3\times10^{-21}$ m³ kg⁻¹, about $10^{-17}$ of specific volume; $D_k$ is zero at every edge
interface to $1.4\times10^{-15}$ m; and the assembled tendency is $2\times10^{-18}$ m s⁻².
`PressureGradCentered` on the same states ranges from $6.8\times10^{-10}$ to $2.6\times10^{-5}$
m s⁻². On a profile *outside* the exact set the same pointwise quantity is $2.9\times10^{-8}$, eleven
orders larger, so the check is not returning zero unconditionally.

Report **absolute** tendencies. At the smallest tilts the residual sits at $2\times10^{-14}$ rather
than $2\times10^{-18}$ m s⁻², not because the scheme is worse but because the floor is round-off of
the $\sim\!7\times10^{3}$ m hydrostatic scale.

**Accuracy off the exact set — the gain is real but modest, and narrows.** On a smooth curved resting
profile the scheme converges at 1.55–2.02, the $O(\tilde h^2)$ §3.7.3 predicts — but
`PressureGradCentered` converges at 1.88–2.93 from a larger starting value, so the advantage falls
from 6.5× at 256 m to 2.4× at 64 m. On a stepped bathymetry it is ~400×. This is what prompted the
restatement of Requirement 2.3.4, and expectations for Omega should be set from these numbers rather
than from the exact-set ones.

**The anchor (§3.5.1).** Exactly zero where the two columns share a surface pressure, which is why
§5.2 puts the exactness gate there. On the one Polaris configuration where they do not, whose sea
surface is initialized from a reference-density balance, it is $6.1\times10^{-4}$ m — the residual of
two 3.57 m terms, of relative size $(\alpha - 1/\rho_0)/\alpha$. That is a gradient the state genuinely
carries, not a defect: measured against a quasi-analytic reference there, the scheme gives an RMS
error of $3.9\times10^{-10}$ against `PressureGradCentered`'s $8.1\times10^{-10}$.

Note for anyone tempted to remove that residual by initializing against the in-situ density instead.
It would shrink but not vanish, because [](#anchor) integrates $\hat\alpha^{(e)}$ rather than the true
$\alpha$, so exact cancellation would require the sea surface to be consistent with the scheme's own
linearized profile. Estimated at the equation-of-state linearization level, some two orders below the
$6.1\times10^{-4}$ m above — derived, not measured.

**Assumption A5 — benign, as §3.7.6 predicted.** Clamping to the outermost valid layer and
extrapolating its reconstruction leaves exactness unaffected, including where the two columns reach
different `maxLevelCell` and where the deepest edge layer extends 65% of its thickness below the
shallower floor.

**Two claims of this document that measurement contradicted.**

1. **Selecting the shared coefficient set by layer index rather than by pressure has no effect**:
   $5.10\times10^{-21}$ against $5.33\times10^{-21}$. §3.5 consequence 1 explains why — whatever the
   set is, it multiplies a contrast that is zero — and §3.7.2 has been corrected accordingly. What
   *is* load-bearing is that one set multiply **both** columns: giving each its own breaks exactness
   at $10^{-5}$.
2. **`PressureGradCentered` is not first order in layer thickness on these profiles.** On a resting
   profile linear in pressure its error falls 3.9× and then 9.8× under successive halvings, checked
   at a tilt where both columns have identical valid-layer counts so that masking is excluded. §3.7.3
   and Requirement 2.3.4 have been corrected.

**Regression status.** `HPGAFiniteVolume` is written alongside the existing `HPGA`, and `HPGA` is
bit-identical to the previous reference run, so adding the scheme moves nothing already measured.

### 6.6 Confirmation from the Omega run

The eight Polaris `horiz_press_grad` variants were run on chrysalis with both schemes over shared
initial conditions. **Every §6.5 prediction was reproduced**, several to four or five significant
figures — the offline Polaris kernel predicted Omega's answers rather than merely resembling them,
which is what licenses using it to size gates before a run in future.

**Resting-state advantage, by variant.** The true HPGA is identically zero in all four, so both
schemes' output is pure error and the ratio is a clean statement of what Phase 1 buys. Smallest
ratio anywhere in each sweep:

| variant | what it tilts | profile | centered / finite-volume |
| --- | --- | --- | --- |
| `hydrostatic_consistency` | the coordinate | curved | $2.58\times$ |
| `hydrostatic_consistency_linear` | the coordinate | in the exact set | $7.71\times$ |
| `bathymetry_step` | the sea floor | curved | $1170\times$ |
| `bathymetry_step_linear` | the sea floor | in the exact set | $58828\times$ |

These are measured with the anchor of [](#anchor) alone. The within-one-layer rule of
[](#anchor-one-layer) leaves 220 of the 224 forward steps bit for bit unchanged, including every
`bathymetry_step` one. The four that change are the 50 m/km tilts of both `hydrostatic_consistency`
variants on 128 m and 64 m layers, where the floors are more than a layer apart; on
`hydrostatic_consistency_linear` the finite-volume error there falls 1.8 times (128 m) and 6.1 times
(64 m).

The ordering is the result, not an artefact of how the sweeps were chosen: the advantage is
smallest where the profile is unresolved and only the coordinate tilts, and largest where the
profile is resolved and the sea floor steps. **A single shared gate would be vacuous at one end and
unreachable at the other**, which is why Polaris sets this per variant.

**Read against Requirement 2.1's bar of 5–10×, three of the four clear it and `hydrostatic_consistency`
does not, at $2.58\times$ — and that is not a shortfall.** §5.3.2 measures the interior coordinate of
a realistic global initial condition as 100–1000× less tilted than the smallest tilt in that sweep,
which puts interior coordinate tilt some 6000× below the observed global interior error. So
`hydrostatic_consistency` is calibrating the scheme on a mechanism that is not the one the bar was set
from, and its ratio should not be read as a prediction of global improvement. The variants that *are*
quantitative proxies for the global bottom-layer error are the stepped-floor ones, and they clear the
bar by two to four orders of magnitude. §5.3.2's calibration paragraph gives the range over which
`bathymetry_step` tracks the global figure, and where it over-predicts.

**`bathymetry_step_linear` is new**, added because no existing variant is both inside the exact set
and stepped: `hydrostatic_consistency_linear` is inside the exact set but tilts only the coordinate,
while `bathymetry_step` has the right geometry but a curved profile whose truncation cannot be told
apart from the geometry's. Being both, it isolates what a bathymetry step costs the scheme.

(Of the two, the stepped floor is the one that matters globally: §5.3.2 measures the interior
coordinate of a realistic global IC as 100–1000× less tilted than the smallest tilt in the Polaris
sweep, putting its contribution some 6000× below the observed global interior error. That rules out
one mechanism; it does not establish A4, which asks whether PGF error *drives* the bottom-layer
instability and which remains undemonstrated. Nothing in this section rests on either.)

**What it establishes, and it is the assumption that was most doubted.** A5 (§3.7.6) — evaluating a
column's deepest reconstruction below its own floor — does **not** break the cancellation. $D_k$
stays flat to $1.4\times10^{-15}$ m at every sea-floor gradient in the sweep, including at 100 and
200 m/km where the two columns reach different `maxLevelCell` and the deepest edge layers extend two
and three layers below one column's floor. A5 is confirmed as an accuracy question rather than a
robustness one, on the geometry where it was least obvious.

**What remains there is the anchor**, not the scan: §5.2's construction condition and §3.7.5's
qualification about the cancellation over a bathymetry step both apply, and the residual on that
variant is $4.3\times10^{-10}$ m s⁻² against `PressureGradCentered`'s $9.1\times10^{-5}$.

**One reading to avoid** is recorded in §3.7.3: the fitted tilt exponent is $1.000$ for *both*
schemes, so it cannot be used to compare them. The absolute level is the discriminating measurement.
