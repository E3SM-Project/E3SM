(omega-design-governing-eqns-omega1)=
# Omega V1: Governing Equations

<!--
Add later, if it seems necessary. There is a toc on the left bar.
**Table of Contents**
1. [Overview](#1-overview)
2. [Requirements](#2-requirements)
3. [Algorithmic Formulation](#3-algorithmic-formulation)
4. [Design](#4-design)
5. [Verification and Testing](#5-verification-and-testing)
-->


## 1. Overview

This design document describes the governing equations for Omega, the Ocean Model for E3SM Global Applications. Overall, Omega is an unstructured-mesh ocean model based on TRiSK numerical methods ([Thuburn et al. 2009](https://www.sciencedirect.com/science/article/pii/S0021999109004434)) that is specifically designed for modern exascale computing architectures. The algorithms in Omega will be mostly identical to those in MPAS-Ocean, but it will be written in c++ rather than Fortran in order to take advantage of the Kokkos performance portability library to run on GPUs ([Trott et al. 2022](https://ieeexplore.ieee.org/document/9485033)). Significant differences between MPAS-Ocean and Omega are:

1. Omega is non-Boussinesq. This means that the full 3D density is used everywhere, and results in a mass-conserving model. MPAS-Ocean and POP were Boussinesq, so that a reference density $\rho_0$ is used in the pressure gradient term, and were therefore volume-conserving models. In Omega the layered mass-conservation equation is in terms of pseudo-thickness ($\tilde{h}=-\Delta p / \rho_0 g$). In MPAS-Ocean the simple thickness ($h=\Delta z$) is the prognostic volume variable (normalized by horizontal cell area).
1. Omega will use the updated equation of state TEOS10, while MPAS-Ocean used the Jackett-McDougall equation of state.

The planned versions of Omega are:

- **Omega-0: Shallow water equations with identical vertical layers and inactive tracers.** In his first version, there is no vertical transport or advection. The tracer equation is horizontal advection-diffusion, but tracers do not feed back to dynamics. Pressure gradient is simply gradient of sea surface height. Capability is similar to [Ringler et al. 2010](https://www.sciencedirect.com/science/article/pii/S0021999109006780)
- **Omega-1.0: Layered ocean, idealized, no surface fluxes.** This adds active temperature, salinity, and density as a function of pressure in the vertical. Vertical advection and diffusion terms are added to the momentum and tracer equations. An equation of state and simple vertical mixing, such as constant coefficient, are needed. Capability and testing are similar to [Petersen et al. 2015](http://www.sciencedirect.com/science/article/pii/S1463500314001796). Tests include overflow, internal gravity wave, baroclinic channel, seamount, and vertical merry-go-round.
- **Omega-1.1: Layered ocean, idealized, with surface fluxes.** Addition of simple vertical mixing scheme such as Pacanowski & Philander; nonlinear equation of state (TEOS10); tracer surface restoring to a constant field; constant-in-time wind forcing; and flux-corrected transport for horizontal advection. Testing will be with the baroclinic gyre and single column tests of surface fluxes and vertical mixing.
- **Omega-2.0: Coupled within E3SM, ocean only.** Ability to run C cases (active ocean only) within E3SM. Requires addition of E3SM coupling infrastructure; simple analysis (time-averaged output of mean, min, max); split baroclinic-barotropic time; global bounds checking on state. Testing and analysis similar to [Ringler et al. 2013](https://www.sciencedirect.com/science/article/pii/S1463500313000760), except Omega uses a non-Boussinesq formulation.
- **Omega-2.1: E3SM fully coupled** Ability to run G cases (active ocean and sea ice) and B cases (all components active) within E3SM. This will include: a full vertical mixing scheme, such as KPP; frazil ice formation in the ocean; and a submosescale parameterization. Omega 2.1 will mostly be run at eddy-resolving resolutions, and will not include parameterizations for lower resolutions such as Gent-McWilliams and Redi mixing.  Simulations will be compared to previous E3SM simulations, including [Caldwell et al. 2019](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2019MS001870) and [Petersen et al. 2019](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2018MS001373).

This document describes the governing equations for the layered ocean model, which are applicable for Omega-1.0 onwards. Specific terms, such as the pressure gradient, vertical mixing, and parameterizations, are left in their general form here but are described in more detail in other design documents.

## 2. Requirements

The requirements in the [Omega-0 design document](OmegaV0ShallowWater) still apply. Additional design requirements for Omega-1 are:

### Omega will be a hydrostatic, non-Boussinesq ocean model.
Omega will adopt a non-Boussinesq formulation, meaning it retains the full, spatially and temporally varying fluid density $\rho({\bf x}, t)$ in all governing equations. This ensures exact mass conservation, as opposed to volume conservation in Boussinesq models that approximate density as a constant reference value $\rho_0$ outside the pressure gradient. The equations are derived in terms of layer-integrated mass per unit area (see [layered equations](#Layer Equations)), and no approximation is made that filters out compressibility or density variations.

The model also assumes hydrostatic balance in the vertical momentum equation, which is a standard and well-justified simplification for large-scale geophysical flows. In such regimes, vertical accelerations are typically small compared to the vertical pressure gradient and gravitational forces. This assumption simplifies the dynamics and removes most sound waves while retaining fidelity for the mesoscale to planetary-scale ocean circulation that Omega is designed to simulate.

### Omega will use TEOS10 for the equation of state.
See additional [EOS design document](EOS)

### Omega-1.0 will add new terms for the pressure gradient, vertical mixing, and vertical advection.
See forthcoming design documents on the pressure gradient, [vertical mixing](VerticalMixingCoeff), and vertical advection.

## 3. Conservation Equations

### Control Volume Formulation

We begin with the continuous, control-volume form of the conservation equations. Consider an arbitrary control volume $V(t)$ with a bounding control surface $\partial V(t)$. The fluid has density $\rho({\bf x},t)$ and velocity ${\bf v}({\bf x},t)$, where $\mathbf{x}$ is the three-dimensional position vector. The control surface is moving at a velocity ${\bf v}_r({\bf x},t)$.  The vector ${\bf n}$ denotes the outward-facing unit normal vector on the control surface $\partial V(t)$, and is used in all surface integrals to represent the direction of fluxes across the boundary. The conservation equations are

**Mass:**

$$
\frac{d}{dt} \int_{V(t)} \rho({\bf x},t) \, dV
+ \int_{\partial V(t)}\rho({\bf x},t)\left({\bf v}({\bf x},t) - {\bf v}_r({\bf x},t) \right) \cdot {\bf n} \, dA
= 0
$$ (continuous-mass)

**Tracers:**

$$
\frac{d}{dt} \int_{V(t)} \rho({\bf x},t) \, \varphi({\bf x},t) \, dV
+ \int_{\partial V(t)}\rho({\bf x},t)\, \varphi({\bf x},t) \left({\bf v}({\bf x},t) - {\bf v}_r \right) \cdot {\bf n} \, dA
= 0
$$ (continuous-tracer)

**Momentum:**

$$
&\frac{d}{dt} \int_{V(t)} \rho({\bf x},t)\,  {\bf v}({\bf x},t) \, dV
+ \int_{\partial V(t)}\rho({\bf x},t)\, {\bf v}({\bf x},t) \left({\bf v}({\bf x},t) - {\bf v}_r \right) \cdot {\bf n} \, dA
\\ & \; \; \; =
\int_{V(t)} \rho({\bf x},t) \, {\bf b}({\bf x},t)\, dV
+ \int_{\partial V(t)} {\bf \tau}({\bf n},{\bf x},t) \, dA
$$ (continuous-momentum)

The operator $\frac{d}{dt}$ used here denotes the rate of change within a moving control volume, sometimes referred to as a Reynolds transport derivative. It differs from the partial derivative $\frac{\partial}{\partial t}$, which represents the local rate of change at a fixed point in space (Eulerian frame), and from the material derivative $\frac{D}{Dt}$, which follows an individual fluid parcel (Lagrangian frame). The use of $\frac{d}{dt}$ allows for conservation laws to be expressed in a general framework that includes both stationary and moving control volumes, consistent with the Reynolds transport theorem.

These equations are taken from [Kundu et al. 2016](https://doi.org/10.1016/C2012-0-00611-4) (4.5) for mass conservation and (4.17) for momentum conservation.
All notation is identical to Kundu, except that we use ${\bf v}$ for the three-dimensional velocity (${\bf u}$ will be used below for horizontal velocities) and ${\bf v}_r$ and $\partial V$ match the notation in [Ringler et al. 2013](https://www.sciencedirect.com/science/article/pii/S1463500313000760) Appendix A.2.


The tracer equation is simply mass conservation, where the conserved quantity is the tracer mass in the control volume $\rho \varphi$, as $\varphi$ is the tracer concentration per unit mass.
In all three equations, the first term is the change of the quantity within the control volume; the second term is the flux through the moving boundary.
If the control surface moves with the fluid in a Lagrangian fashion, then ${\bf v}_r={\bf v}$ and the second term (flux through the boundary) vanishes---there is no net mass, momentum, or tracer transport across the moving surface.


The momentum equation is an expression of Newton's second law and includes two types of external forces.
The first is the body force, represented here as ${\bf b}({\bf x}, t)$, which encompasses any volumetric force acting throughout the fluid, such as gravitational acceleration or the Coriolis force. In some contexts, body forces may be expressible as the gradient of a potential, ${\bf b} = -\nabla_{3D} \Phi$, but this is not assumed in general.
The second is the surface force or traction vector ${\bf \tau} = \sigma \cdot {\bf n}$, which acts on the boundary of the control volume and includes pressure and viscous stresses. These forces appear as surface integrals over the boundary and drive momentum exchange between adjacent fluid parcels or between the fluid and its environment.
The derivation of the momentum equation may also be found in [Leishman 2025](https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/conservation-of-momentum-momentum-equation/#chapter-260-section-2), Chapter 21, equation 10.

## 4. Hydrostatic Approximation

In large-scale geophysical flows, vertical accelerations are typically much smaller than the vertical pressure gradient and gravitational forces. As a result, we adopt the **hydrostatic approximation**, which simplifies the vertical momentum equation by assuming a leading-order balance between pressure and gravity:

$$
\frac{\partial p}{\partial z} = -\rho g
$$ (hydrostatic)

This balance allows us to express the pressure at any depth $z$ in terms of an integral over the density field above:

$$
p(x,y,z,t) = p^{\text{surf}}(x,y,t) + \int_{z}^{z^{\text{surf}}} \rho(x,y,z',t) g \, dz'
$$ (pressure)

Here, $p^{\text{surf}}$ is the pressure at the free surface $z^{\text{surf}}(x, y, t)$, $\rho(z')$ is the local fluid density and $z'$ is a dummy variable of integration. This relation holds pointwise in space and time and provides a foundation for defining vertical coordinates based on pressure.

The hydrostatic approximation also simplifies vertical pressure forces in control-volume integrations. For example, in a finite-volume cell bounded above and below by sloping surfaces, the vertical pressure force reduces to:

$$
- p|_{z = z^{\text{top}}} + p|_{z = z^{\text{bot}}}
$$ (vert-pressure-balance)

because the area correction from the slope and the projection of the pressure gradient into the vertical direction cancel to leading order.  [](#vert-pressure-balance) allows us to create simpler relations between pressure and pseudo-height in the following section.

This hydrostatic balance underlies the definition of **pseudo-height**, introduced next, and ensures consistency between the prognostic vertical coordinate and the model’s vertical force balance.

## 5. Pseudo-Height

In our non-Boussinesq, hydrostatic ocean model, we adopt **pseudo-height** $\tilde{z}$ as the prognostic vertical coordinate. This choice is motivated directly by the hydrostatic approximation introduced in [](#hydrostatic).
Under this balance, we define the pseudo-height as a normalized pressure:

$$
\tilde{z}(x, y, z, t) = -\frac{1}{\rho_0 g} \, p(x, y, z, t)
$$ (pseudo-height)

Here, $\rho_0$ is a constant reference density and $g$ is gravitational acceleration. Although $\tilde{z}$ is essentially a pressure coordinate, the normalization ensures it has units of length and retains an intuitive geometric interpretation.

The inverse relation, giving $\tilde{z}$ in terms of geometric height, follows from the hydrostatic pressure integral:

$$
\tilde{z}(x, y, z, t) = -\frac{1}{\rho_0 g} \left( p(x, y, z, t)^{\text{surf}} + \int_{z}^{z^{\text{surf}}} \rho(z') g \, dz' \right)
$$ (pseudo-heigh-from-z)

This expression shows that pseudo-height varies in space and time like pressure. It also makes clear that $z$ is no longer a prognostic variable but a derived quantity, recoverable via hydrostatic inversion.

Taking a differential, we relate pseudo-height and geometric height:

$$
d\tilde{z} = -\frac{1}{\rho_0 g} \, dp = \frac{\rho}{\rho_0} \, dz
$$ (d-z-tilde)

This motivates the definition of **pseudo-velocity** in the vertical direction:

$$
\tilde{w} \equiv \frac{d\tilde{z}}{dt} = \frac{\rho}{\rho_0} \, w
$$ (w-tilde)

In this formulation:
- $\tilde{w}$ is the vertical **mass flux per unit reference density**, a key variable in a non-Boussinesq framework.
- Vertical transport of mass, tracers, and momentum will use $\tilde{w}$ rather than the vertical velocity in height coordinates, $w$.

[Griffies 2018](https://doi.org/10.2307/j.ctv301gzg) p. 37  argues for the use of pseudo-velocities, which he calls the density-weighted velocity, for non-Boussinesq models.
Griffies recommends a value of $\rho_0=1035$ kg/m$^3$, following p. 47 of [Gill (1982)](https://doi.org/10.1016/S0074-6142(08)60028-5), because ocean density varies less than 2% from that value.

The use of a constant $\rho_0$ in defining pseudo-height does not imply the Boussinesq approximation. In Boussinesq models, $\rho$ is set to $\rho_0$ everywhere except in the buoyancy term (i.e., the vertical pressure gradient or gravitational forcing). Here, by contrast, we retain the full $\rho$ in all terms, and use $\rho_0$ only as a normalization constant—for example, so that $d\tilde{z} \approx dz$ when $\rho \approx \rho_0$. This preserves full mass conservation while making vertical units more intuitive.

This approach ensures consistency with mass conservation and simplifies the derivation of the discrete equations.

> **Note:** This definition does not invoke the Boussinesq approximation. The full, time- and space-dependent density $\rho(x, y, z, t)$ is retained in all conservation laws. The use of $\rho_0$ is solely for normalization.

### 5.1. Justification for Pseudo-Height Definition

Alternative definitions of pseudo-height could add an arbitrary offset, e.g., setting $\tilde{z} = 0$ at mean sea level. However, the adopted form:

$$
\tilde{z} = -\frac{1}{\rho_0 g} \, p
$$

was chosen for its practical advantages:
- It ensures that $\tilde{z}$ varies identically to pressure.
- Units of $\tilde{z}$, $\tilde{h}$, and $\tilde{w}$ are intuitive.
- It aligns layer thickness with pseudo-thickness: $\tilde{h} = \Delta \tilde{z} \propto \Delta p$
- $\tilde{z}^{\text{floor}} \propto p^{\text{bot}}$, which aids barotropic pressure gradient calculation.

This makes $\tilde{z}$ a natural coordinate for a mass-conserving hydrostatic model.

## 6. Horizontal and Vertical Flux Separation in Pseudo-Height Coordinates

In geophysical flows, the vertical and horizontal directions are treated differently due to rotation and stratification, which lead to distinct characteristic scales of motion. To reflect this, we separate horizontal and vertical fluxes explicitly in the governing equations.

We partition the control surface of a finite-volume cell into the side walls $\partial V^{\text{side}}$ (which are fixed in space and time) and the upper and lower pseudo-height surfaces $\partial V^{\text{top}}(t)$ and $\partial V^{\text{bot}}(t)$, which evolve in time.

As an example, we consider the tracer equation and drop explicit notation for spatial and temporal dependence for clarity. We write the control-volume form as:

$$
\frac{d}{dt} \int_{V(t)} \rho \, \varphi \, dV
& + \int_{\partial V^{\text{side}}} \rho \varphi ({\bf v} - {\bf v}_r) \cdot {\bf n} \, dA \\
& + \int_{\partial V^{\text{top}}(t)} \rho \varphi ({\bf v} - {\bf v}_r) \cdot {\bf n}^{\text{top}} \, dA
+ \int_{\partial V^{\text{bot}}(t)} \rho \varphi ({\bf v} - {\bf v}_r) \cdot {\bf n}^{\text{bot}} \, dA
= 0
$$ (tr-vh-split-pseudo)

We now express the top and bottom surfaces using the pseudo-height variable $\tilde{z}(x,y,t)$ rather than geometric height $z$. The unit normals are therefore:

$$
{\bf n}^{\text{top}} \approx (-\nabla \tilde{z}^{\text{top}}, 1), \quad
{\bf n}^{\text{bot}} \approx (\nabla \tilde{z}^{\text{bot}}, -1)
$$ (top-bot-normal-pseudo)

Here we apply a small-slope approximation, assuming $|\nabla \tilde{z}| \ll 1$, and omit normalization factors for clarity. This approximation retains leading-order effects of slope while simplifying the surface geometry.

To facilitate integration over a fixed horizontal domain, we introduce the following notation:
- $A$ is the horizontal footprint (in the $x$–$y$ plane) of the control volume $V(t)$.
- $dA$ is the horizontal area element.
- $\partial A$ is the boundary of $A$, and $dl$ is the line element along this boundary.
- $\mathbf{n}_\perp$ is the outward-pointing unit normal vector in the horizontal plane, defined on $\partial A$. It lies in the $x$–$y$ plane and is orthogonal to $dl$.

We now project the tracer equation [](#tr-vh-split-pseudo) onto a horizontal domain $A$ (the footprint of the control volume), over which the top and bottom boundaries vary in height. The side walls remain fixed in time and space, but not vertical extent.


The pseudo-height surfaces $\tilde{z}^{\text{top}}(x,y,t)$ and $\tilde{z}^{\text{bot}}(x,y,t)$ are not fixed in time, so their motion must be accounted for when computing vertical fluxes. We define $\tilde{w}_r$ as the pseudo-height velocity of the interface:

$$
\tilde{w}_{r} = \left.\frac{d\tilde{z}}{dt}\right |_{\text interface}
$$ (interface-velocity)

The **net vertical transport** through a moving surface is then

$$
\tilde{w}_{tr} = \tilde{w} - \tilde{w}_r
$$ (w-tilde-tr)

This gives the relative pseudo-mass flux through the interface per unit reference density. Vertical flux terms are written using this difference to ensure conservation in the presence of interface motion.

In the flux integrals below, we express all surface-normal transport using $\tilde{w}_{tr}$, i.e.,

$$
\rho_0 \varphi \left[ \tilde{w}_{tr} - \tilde{ u} \right]
= \rho_0 \varphi \left[ \tilde{w} - \tilde{w}_r - \tilde{ u} \right]
$$ (vertical-flux)

so that the vertical flux terms correctly account for both the local motion of the fluid and the motion of the interface itself.

The quantity $\tilde{w}$ in equation [](#w-tilde) is the pseudo-velocity. The second term, $\tilde{ u} = {\bf u} \cdot \nabla \tilde{z}^{\text{top}}$, captures the component of horizontal velocity advecting material across a sloping pseudo-height interface. Together, the expression in brackets represents the **net vertical transport** through a surface of constant $\tilde{z}$. A schematic of how $\tilde{w}$, $\tilde{w}_{tr}$, and $\tilde{u}$ relate to each other is shown in {numref}`pseudo-vel-schematic`:

```{figure} images/tilde_w.jpeg
:name:  pseudo-vel-schematic
:align: center
:width: 300 px

Schematic of pseudo-velocities.
```
Because both terms are scaled by $1/\rho_0$ in [](#vertical-flux), the full flux is multiplied by $\rho_0$ to recover a physical mass flux. This ensures dimensional consistency and physical equivalence to the traditional form $\rho w - \rho {\bf u} \cdot \nabla z$ expressed in geometric-height coordinates.

$\tilde{W}_{tr}$ will be introduced and defined in [](#notational-simplifications). Using this projection, we obtain:

$$
\frac{d}{dt} \int_A \int_{\tilde{z}^{\text{bot}}}^{\tilde{z}^{\text{top}}} \, \varphi \, d\tilde{z} \, dA
& + \int_{\partial A} \left( \int_{\tilde{z}^{\text{bot}}}^{\tilde{z}^{\text{top}}} \varphi \, {\bf u} \, d\tilde{z} \right) \cdot \mathbf{n}_\perp \, dl \\
& + \int_A \rho_0 \left\{\varphi \left[ \tilde{w}_{tr} - \tilde{ u} \right]\right\}_{\tilde{z} = \tilde{z}^{\text{top}}} dA \\
& - \int_A \rho_0 \left\{\varphi \left[ \tilde{w}_{tr} - \tilde{ u} \right]\right\}_{\tilde{z} = \tilde{z}^{\text{bot}}} dA
= 0
$$ (tr-vh-separation-pseudo)

This equation is structurally identical to its $z$-based counterpart but with all references to geometric height replaced by pseudo-height quantities.

This formulation allows vertical mass and tracer transport to be computed directly in terms of prognostic variables, without reconstructing $\rho$ or $z$ explicitly.

Similar expressions to [](#tr-vh-separation-pseudo) hold for mass and momentum:

**Mass:**

$$
\frac{d}{dt} \int_A \int_{\tilde{z}^{\text{bot}}}^{\tilde{z}^{\text{top}}} \, d\tilde{z} \, dA
& + \int_{\partial A} \left( \int_{\tilde{z}^{\text{bot}}}^{\tilde{z}^{\text{top}}} \, {\bf u} \, d\tilde{z} \right) \cdot \mathbf{n}_\perp \, dl \\
& + \int_A \left\{\left[ \tilde{w}_{tr} - \tilde{ u} \right]\right\}_{\tilde{z} = \tilde{z}^{\text{top}}} dA \\
& - \int_A \left\{\left[ \tilde{w}_{tr} - \tilde{ u} \right]\right\}_{\tilde{z} = \tilde{z}^{\text{bot}}} dA
= 0
$$ (vh-mass-pseudo)

**Tracer:**

$$
\frac{d}{dt} \int_A \int_{\tilde{z}^{\text{bot}}}^{\tilde{z}^{\text{top}}}  \varphi \, d\tilde{z} \, dA
& + \int_{\partial A} \left( \int_{\tilde{z}^{\text{bot}}}^{\tilde{z}^{\text{top}}}  \varphi \, {\bf u} \, d\tilde{z} \right) \cdot \mathbf{n}_\perp \, dl & \\
& + \int_A \left\{\varphi \left[\tilde{w}_{tr} - \tilde{ u} \right] \right\}_{\tilde{z} = \tilde{z}^{\text{top}}} \, dA & \\
& - \int_A \left\{\varphi \left[\tilde{w}_{tr} - \tilde{ u} \right] \right\}_{\tilde{z} = \tilde{z}^{\text{bot}}}\, dA
& = 0
$$ (vh-tracer-pseudo)

**Momentum:**

$$
\frac{d}{dt} \int_{A} \int_{\tilde{z}^{\text{bot}}}^{\tilde{z}^{\text{top}}} \, {\bf u} \, d\tilde{z} \, dA
&+
\int_{\partial A} \left( \int_{\tilde{z}^{\text{bot}}}^{\tilde{z}^{\text{top}}} \, {\bf u} \otimes {\bf u} \, d\tilde{z} \right) \cdot \mathbf{n}_\perp \, dl \\
&+
\int_A \, \left\{{\bf u} \left[ \tilde{w}_{tr} - \tilde{ u} \right]\right\}_{\tilde{z} = \tilde{z}^{\text{top}}} \, dA
-
\int_A \, \left\{{\bf u} \left[ \tilde{w}_{tr} - \tilde{ u} \right]\right\}_{\tilde{z} = \tilde{z}^{\text{bot}}} \, dA \\
&=
\int_A \int_{\tilde{z}^{\text{bot}}}^{\tilde{z}^{\text{top}}} {\bf b}_{\perp} \, d\tilde{z} \, dA
+
\int_{\partial A} \left( \int_{\tilde{z}^{\text{bot}}}^{\tilde{z}^{\text{top}}} \frac{{\bf \tau}_\perp}{\rho} \, d\tilde{z} \right) \mathbf{n}_\perp dl \\
&\quad
+ \int_A \left[ \frac{{\bf \tau}_\perp}{\rho} \right]_{\tilde{z} = \tilde{z}^{\text{top}}} \, dA
- \int_A \left[ \frac{{\bf \tau}_\perp}{\rho} \right]_{\tilde{z} = \tilde{z}^{\text{bot}}} \, dA
$$ (vh-momentum-pseudo)

Here, we have divided by $\rho_0$, since it occurs in every term.

In [](#vh-momentum-pseudo), $\mathbf{b}_\perp$ is defined similarly to $\mathbf{n}_\perp$, it is the body force normal to $dl$, and ${\bf \tau}_\perp$ is the horizontal component of the traction vector. These equations express horizontal and vertical fluxes naturally in terms of pseudo-height, enabling fully consistent discretization in the vertical coordinate used for prognostic evolution.

## 7. Vertical Discretization for the Layered Equations

The previous equation set [](#vh-mass-pseudo) to [](#vh-momentum-pseudo) is for an arbitrary layer bounded by $\tilde{z}^{\text{top}}$ above and $\tilde{z}^{\text{bot}}$ below. We now provide the details of the vertical discretization. The ocean is divided vertically into $K_{max}$ layers, with $k=0$ at the top and increasing downwards (opposite from $z$). Layer $k$ is bounded between $\tilde{z}_k^{\text{top}}$ above and $\tilde{z}_{k+1}^{\text{top}}$ below (i.e. $\tilde{z}_k^{\text{bot}} = \tilde{z}_{k+1}^{\text{top}}$).

The layer thickness of layer $k$, used in MPAS-Ocean, is

$$
h_k = \int_{z_{k}^{\text{bot}}}^{z_k^{\text{top}}} dz.
$$ (def-thickness)

In Omega we will use the pseudo-thickness,

$$
{\tilde h}_k(x,y,t)
&= \int_{{\tilde z}_{k}^{\text{bot}}}^{{\tilde z}_k^{\text{top}}} d\tilde{z} \\
&= \frac{1}{\rho_0} \int_{z_{k}^{\text{bot}}}^{z_k^{\text{top}}} \rho \, dz \\
&= \frac{1}{\rho_0 g} \left( p_{k}^{\text{bot}} - p_k^{\text{top}} \right)
$$ (def-pseudo-thickness)

which is the mass per unit area in the layer, normalized by $\rho_0$. This pseudo-thickness and layer-averaging will be used to express conservation laws in a mass-weighted coordinate system.  Pseudo-thickness, rather than geometric thickness will be the prognostic variable in Omega.

The density-weighted average of any variable $\varphi({\bf x},t)$ in layer $k$ is

$$
{\overline \varphi}^{\tilde{z}}_k(x,y,t) =
\frac{\int_{z_{k}^{\text{bot}}}^{z_k^{\text{top}}} \rho \varphi dz}
     {\int_{z_{k}^{\text{bot}}}^{z_k^{\text{top}}} \rho dz}
     =
\frac{\frac{1}{\rho_0}\int_{z_{k}^{\text{bot}}}^{z_k^{\text{top}}} \rho \varphi dz}
     {\frac{1}{\rho_0}\int_{z_{k}^{\text{bot}}}^{z_k^{\text{top}}} \rho dz}
     =
\frac{\int_{{\tilde z}_{k}^{\text{bot}}}^{{\tilde z}_k^{\text{top}}} \varphi d{\tilde z}}
     {\int_{{\tilde z}_{k}^{\text{bot}}}^{{\tilde z}_k^{\text{top}}} d{\tilde z}}
=
\frac{1}{{\tilde h}_k}\int_{{\tilde z}_{k}^{\text{bot}}}^{{\tilde z}_k^{\text{top}}} \varphi d{\tilde z}.
$$ (def-layer-average)

Rearranging, it is useful to note that

$$
\int_{\tilde{z}_{k}^{\text{bot}}}^{\tilde{z}_k^{\text{top}}} \varphi d\tilde{z}
 = {\tilde h}_k {\overline \varphi}^{\tilde{z}}_k.
$$ (h-phi)

This relation is frequently used in discretized fluxes and conservation equations to replace integrals with layer-mean quantities.

## 8. Decomposition and Averaging

### Decomposition

In this document, two decompositions will be critical.  The first decomposition is

$$
\varphi(x,y,z,t) = \overline{\varphi}^{\tilde{z}}_k(x,y,t) + \delta \varphi(x,y,z,t)
$$ (vertical-decomposition)

which is a density weighted vertical integral based upon [](#def-layer-average) and the deviation from this value within the layer.

The second decomposition is

$$
\varphi = \left<\varphi\right> + \varphi^\prime
$$ (reynolds-definition)

which is the traditional Reynolds' average and deviation from this value.

The fundamental ocean model equations are most often Reynolds' averaged to derive the sub gridscale stresses. Each variable is decomposed into an average over a sufficiently large sample of turbulent processes, which allows for a meaningful fluctuating component [[](#reynolds-definition)]. Commonly, the Reynolds' average is thought of as a time average, but an ensemble average is equally valid.  In fact, the ensemble and time averaging can be thought of as equivalent, given that a sufficiently long time average will effectively average over variations in the flow that could be thought of as an ensemble.  Given the functions are continuous, the averaging operator can be passed through derivatives and integrals without corrections and an average of products with a single perturbation quantity is zero.

The Reynolds' average is most commonly denoted by an overbar. However, to disambiguate from the definition of the bar as the vertical density weighted average, the Reynolds' average herein is denoted by $< . >$.

The Reynolds' approach is an attractive approach for Boussinesq ocean models since the fundamental equations do not include products of spatially variable density and tracer, pressure, or momentum.  When the ocean model is non Boussinesq, products of spatially varying density and other fields (e.g., tracer) arise and create difficulties, producing a term like $\left<\rho^\prime {\mathbf u}^\prime\right>$ which is difficult to parameterize. However, given the definition of our chosen pseudo-height, the density terms are wrapped up in $\tilde{z}$ and once again a Reynolds' approach can be cleanly used.

### Averaging

The quantity being averaged in [](#def-layer-average) is arbitrary. For example, this equation can apply equally to a Reynolds' averaged quantity, e.g.,

$$
\overline{\left<\varphi\right>}^{\tilde{z}}_k(x,y,t) =
\frac{1}{\tilde{h}_k} \int_{\tilde{z}_{k}^{\text{bot}}}^{\tilde{z}_k^{\text{top}}} \left<\varphi\right> d\tilde{z}.
$$ (def-layer-average-reynolds)

Additionally, for the decomposition related to vertical averaging, we will encounter terms such as ${\overline \varphi}^{\tilde{z}}_k \delta\varphi$ that need to be vertical averaged.  Using [](#def-layer-average) we have

$$
\overline{{\overline \varphi}^{\tilde{z}}_k \delta\varphi}^{\tilde{z}}_k &= \frac{\int_{{\tilde z}_{k}^{\text{bot}}}^{{\tilde z}_k^{\text{top}}} {\overline \varphi}^{\tilde{z}}_k \delta \varphi d{\tilde z}}
     {\int_{{\tilde z}_{k}^{\text{bot}}}^{{\tilde z}_k^{\text{top}}} d{\tilde z}} \\
     &= {\overline \varphi}^{\tilde{z}}_k \frac{\int_{{\tilde z}_{k}^{\text{bot}}}^{{\tilde z}_k^{\text{top}}} \delta \varphi d{\tilde z}}
     {\int_{{\tilde z}_{k}^{\text{bot}}}^{{\tilde z}_k^{\text{top}}} d{\tilde z}} \\
     &= 0,
$$ (delta-vert-average)

where the last equality is true by definition.

## 9. Layer Equations

### Tracer & Mass

We first Reynolds' average [](#vh-tracer-pseudo) and given the definition of the operator we can move the averaging past the derivatives and integrals without correction terms to yield

$$
\frac{d}{dt} \int_A \int_{\tilde{z}_k^{\text{bot}}}^{\tilde{z}_k^{\text{top}}} \left<\varphi\right> \, d\tilde{z} \, dA
& + \int_{\partial A} \left( \int_{\tilde{z}_k^{\text{bot}}}^{\tilde{z}_k^{\text{top}}} \left<\varphi \, {\bf u} \right> \, d\tilde{z} \right) \cdot {\bf n}_\perp \, dl & \\
& + \int_A \left<\left\{ \varphi \left[\tilde{w}_{tr} - \tilde{ u} \right] \right\}\right>_{\tilde{z} = \tilde{z}_k^{\text{top}}} \, dA & \\
& - \int_A \left<\left\{ \varphi \left[\tilde{w}_{tr} - \tilde{ u} \right] \right\}\right>_{\tilde{z} = \tilde{z}_k^{\text{bot}}}\, dA
& = 0.
$$ (Aintegral-tracer-first-reynolds)

Next, we use a Reynolds' decomposition on any terms involving products inside a Reynolds' average,

$$
\frac{d}{dt} \int_A \int_{\tilde{z}_{k}^{\text{bot}}}^{\tilde{z}_k^{\text{top}}} \left<\varphi \right> \, d\tilde{z} \, dA
& +
\int_{\partial A} \left( \int_{\tilde{z}_{k}^{\text{bot}}}^{\tilde{z}_k^{\text{top}}} \left(\left<{\varphi}\right>\left<{\bf u}\right> + \left<\varphi^{\prime}{\bf u}^{\prime}\right> \right) d\tilde{z} \right) \, \cdot \mathbf{n}_\perp \, dl & \\
& + \int_A \left\{ \left[ \left<\varphi\right>\left<\tilde{w}_{tr}\right> + \left<\varphi^{\prime} \tilde{w}_{tr}^\prime \right> - \left<\varphi\right> \left<\tilde{ u}\right> -\left<\varphi^{\prime} \tilde{ u}^{\prime}\right> \right] \right\}_{\tilde{z} = \tilde{z}_k^{\text{top}}} \, dA & \\
& - \int_A \left\{ \left[ \left<\varphi\right>\left<\tilde{w}_{tr}\right> + \left<\varphi^{\prime} \tilde{w}_{tr}^\prime \right> - \left<\varphi\right> \left<\tilde{ u}\right> - \left<\varphi^{\prime} \tilde{ u}^{\prime}\right> \right] \right\}_{\tilde{z} = \tilde{z}_k^{\text{bot}}}\, dA
& = 0.
$$ (Aintegral-tracer-reynolds)

We can rewrite the first two terms in [](#Aintegral-tracer-reynolds) using [](#def-layer-average-reynolds)

$$
\frac{d}{dt} \int_A \tilde{h}_k \, \overline{\left<\varphi\right>}^{\tilde{z}}_k \, dA
& + \int_{\partial A} \left( \tilde{h}_k \, \left[\overline{\left<\varphi\right> \left<{\bf u}\right> + \left<\varphi^{\prime}{\bf u}^{\prime}\right>}^{\tilde{z}}_k \right] \right) \cdot \mathbf{n}_\perp \, dl & \\
& + \int_A \left\{ \left[ \left<\varphi\right>\left<\tilde{w}_{tr}\right> + \left<\varphi^{\prime} \tilde{w}_{tr}^\prime \right> - \left<\varphi\right> \left<\tilde{ u}\right> -\left<\varphi^{\prime} \tilde{ u}^{\prime}\right> \right] \right\}_{\tilde{z} = \tilde{z}_k^{\text{top}}} \, dA & \\
& - \int_A \left\{ \left[ \left<\varphi\right>\left<\tilde{w}_{tr}\right> + \left<\varphi^{\prime} \tilde{w}_{tr}^\prime \right> - \left<\varphi\right> \left<\tilde{ u}\right> - \left<\varphi^{\prime} \tilde{ u}^{\prime}\right> \right] \right\}_{\tilde{z} = \tilde{z}_k^{\text{bot}}}\, dA
& = 0.
$$ (Aintegral-tracer)

The second term is expanded utilizing [](#vertical-decomposition)

$$
\frac{d}{dt} \int_A \tilde{h}_k \, \overline{\left<\varphi\right>}^{\tilde{z}}_k \, dA
& + \int_{\partial A} \left( \tilde{h}_k \, \overline{\left(\overline{\left<\varphi\right>}^{\tilde{z}}_k + \delta \varphi\right)\left(\overline{\left<{\bf u}\right>}^{\tilde{z}}_k + \delta {\bf u}\right)}^{\tilde{z}}_k + \tilde{h}_k \overline{\left<\varphi^{\prime}{\bf u}^{\prime}\right>}^{\tilde{z}}_k \right) \cdot \mathbf{n}_\perp \, dl & \\
& + \int_A \left\{ \left[ \left<\varphi\right>\left<\tilde{w}_{tr}\right> + \left<\varphi^{\prime} \tilde{w}_{tr}^\prime \right> - \left<\varphi\right> \left<\tilde{ u}\right> -\left<\varphi^{\prime} \tilde{ u}^{\prime}\right> \right] \right\}_{\tilde{z} = \tilde{z}_k^{\text{top}}} \, dA & \\
& - \int_A \left\{ \left[ \left<\varphi\right>\left<\tilde{w}_{tr}\right> + \left<\varphi^{\prime} \tilde{w}_{tr}^\prime \right> - \left<\varphi\right> \left<\tilde{ u}\right> - \left<\varphi^{\prime} \tilde{ u}^{\prime}\right> \right] \right\}_{\tilde{z} = \tilde{z}_k^{\text{bot}}}\, dA
& = 0.
$$ (Aintegral-tracer2)

And then is simplified using [](#delta-vert-average)

$$
\frac{d}{dt} \int_A \tilde{h}_k \, \overline{\left<\varphi\right>}^{\tilde{z}}_k \, dA
& + \int_{\partial A} \left[ \tilde{h}_k \, \left(\overline{\left<\varphi\right>}^{\tilde{z}}_k \overline{\left<{\bf u}\right>}^{\tilde{z}}_k + \overline{\delta \varphi \delta {\bf u}}^{\tilde{z}}_k + \overline{\left<\varphi^{\prime}{\bf u}^{\prime}\right>}^{\tilde{z}}_k \right)\right] \cdot \mathbf{n}_\perp \, dl & \\
& + \int_A \left\{ \left[ \left<\varphi\right>\left<\tilde{w}_{tr}\right> + \left<\varphi^{\prime} \tilde{w}_{tr}^\prime \right> - \left<\varphi\right> \left<\tilde{ u}\right> -\left<\varphi^{\prime} \tilde{ u}^{\prime}\right> \right] \right\}_{\tilde{z} = \tilde{z}_k^{\text{top}}} \, dA & \\
& - \int_A \left\{ \left[ \left<\varphi\right>\left<\tilde{w}_{tr}\right> + \left<\varphi^{\prime} \tilde{w}_{tr}^\prime \right> - \left<\varphi\right> \left<\tilde{ u}\right> - \left<\varphi^{\prime} \tilde{ u}^{\prime}\right> \right] \right\}_{\tilde{z} = \tilde{z}_k^{\text{bot}}}\, dA
& = 0.
$$

Given that the horizontal area is not changing in time, we can move the time derivative through the first integral.  We also invoke Green's Theorem for the second integral to convert that from a surface integral to an area integral.  With these changes, the equation becomes

$$
\int_A  \bigl\{ \frac{\partial \tilde{h}_k \, \overline{\left<\varphi\right>}^{\tilde{z}}_k}{\partial t}
& + \nabla \cdot \left[ \tilde{h}_k \, \left(\overline{\left<\varphi\right>}^{\tilde{z}}_k \overline{\left<{\bf u}\right>}^{\tilde{z}}_k + \overline{\delta \varphi \delta {\bf u}}^{\tilde{z}}_k + \overline{\left<\varphi^{\prime}{\bf u}^{\prime}\right>}^{\tilde{z}}_k \right)\right] & \\
& + \left\{ \left[ \left<\varphi\right>\left<\tilde{w}_{tr}\right> + \left<\varphi^{\prime} \tilde{w}_{tr}^\prime \right> - \left<\varphi\right> \left<\tilde{ u}\right> - \left<\varphi^{\prime} \tilde{ u}^{\prime}\right> \right] \right\}_{\tilde{z} = \tilde{z}_k^{\text{top}}} & \\
& - \left\{ \left[ \left<\varphi\right>\left<\tilde{w}_{tr}\right> + \left<\varphi^{\prime} \tilde{w}_{tr}^\prime \right> - \left<\varphi\right> \left<\tilde{ u}\right> - \left<\varphi^{\prime} \tilde{ u}^{\prime}\right> \right] \right\}_{\tilde{z} = \tilde{z}_k^{\text{bot}}} \bigr\}  \, dA
& = 0.
$$ (reynolds-tracer-final)


A few comments on the last two lines of [](#reynolds-tracer-final).  If only the first two terms of the integrals are retained, this represents the vertical advective and turbulent fluxes for a finite volume framework.  The second two terms represent the projection of the horizontal advective and turbulent fluxes into the direction normal to the sloping surface.  When the coordinate surfaces are flat $\nabla \tilde{z}^{\text{top}} = \nabla \tilde{z}^{\text{bot}} = 0$ and the vertical fluxes are aligned with the normal to the coordinate surface reducing the equation to the expected z-coordinate, finite volume, formulation.

Finally, given that the area integral operates on the entire equation, it is valid for any area and we can take the limit as $A \rightarrow 0$ to arrive at the final tracer equation

**Tracer:**

$$
\frac{\partial {\tilde h}_k \overline{\left<\varphi\right>}^{\tilde{z}}_k  }{\partial t}  & + \nabla \cdot \left({\tilde h}_k \overline{\left<\varphi\right>}^{\tilde{z}}_k \overline{\left<{\bf u}\right>}^{\tilde{z}}_k\right) \\
& + \left\{ \left[ \left<\varphi\right> \left<{\tilde w}_{tr}\right> - \left<\varphi\right>\left<\tilde{ u}\right> \right]_{{\tilde z}={\tilde z}_k^{\text{top}}} - \left[ \left<\varphi\right> \left<{\tilde w}_{tr}\right> - \left<\varphi\right>\left<\tilde{ u}\right> \right]_{{\tilde z}={\tilde z}_k^{\text{bot}}} \right\} \\
& = - \nabla \cdot \left({\tilde h}_k \overline{\left<\varphi^{\prime} {\bf u}^{\prime}\right>}^{\tilde{z}}_k \, + {\tilde h}_k \overline{\delta \varphi \delta {\bf u}}^{\tilde{z}}_k \right) \\
& - \left\{ \left[\left< \varphi^\prime {\tilde w}_{tr}^{\prime} \right> - \left< \varphi^\prime \tilde{ u}^{\prime}\right> \right]_{{\tilde z}={\tilde z}_{k}^{\text{top}}} - \left[\left< \varphi^\prime {\tilde w}_{tr}^{\prime} \right> - \left< \varphi^\prime \tilde{ u}^{\prime}\right> \right]_{{\tilde z}={\tilde z}_{k}^{\text{bot}}}\right\}.
$$ (layer-tracer)

A few notes on the layer-averaged tracer equation. In this complete form, it includes three types of fluctuating quantities that must be dealt with: (1) the layer averaged, Reynolds' averaged horizontal turbulent flux $\left( \overline{\left<\varphi^{\prime} {\bf u}^{\prime}\right>}^{\tilde{z}}_k \right)$, (2) the Reynolds' average vertical turbulent flux $\left( \left< \varphi^\prime {\tilde w}_{tr}^{\prime} \right> - \left< \varphi^\prime \tilde{u}^{\prime}\right> \right)$, and (3) the layer-averaged product of deviations from the layer-integrated variables $\left(\overline{\delta \varphi \delta {\bf u}}^{\tilde{z}}_k \right)$. The details of the first two quantities will be discussed later in this document and follow-on design documents. The terms involving perturbations from the layer integrated quantity are necessary to extend beyond piecewise constant represenation of variables. In this equation, variables with no overline are the full field variable at the interfaces.

The mass equation is identical to the tracer equation with $\varphi=1$.

**Mass:**

$$
\frac{\partial {\tilde h}_k }{\partial t}
+ \nabla \cdot \left({\tilde h}_k \overline{\left<{\bf u}\right>}^{\tilde{z}}_k\right)
+ \left[ \left<{\tilde w}_{tr}\right> - \left<\tilde{ u}\right> \right]_{{\tilde z}={\tilde z}_k^{\text{top}}}
- \left[ \left<{\tilde w}_{tr}\right> - \left<\tilde{ u}\right> \right]_{{\tilde z}={\tilde z}_{k}^{\text{bot}}}
= 0.
$$ (layer-mass)

### Pressure Gradient Force

Before deriving the layered momentum equation, we collect the pressure and gravitational contributions into a single **pressure gradient force (PGF)** in our pseudo-height framework. This subsection establishes the continuous form of the horizontal PGF per unit mass and serves as the target for the subsequent finite-volume, layer-averaged derivation.

Recall from [](vh-momentum-pseudo) that the horizontal momentum equation is written in terms of

- a **body force per unit mass** $\mathbf{b}_\perp$, and
- a **traction** vector $\boldsymbol{\tau}$ that enters the weak form as $\boldsymbol{\tau}/\rho$,

so that all terms in the equation have the units of acceleration.

In the finite-volume formulation, pressure acts on the boundary of a control volume $\mathcal{R}$ through the traction
$\boldsymbol{\tau}_p = -p\,\mathbf{n}$, where $\mathbf{n}$ is the outward unit normal. The net pressure force on $\mathcal{R}$ is therefore

$$
\mathbf{F}_p(\mathcal{R})
=
\int_{\partial\mathcal{R}} \boldsymbol{\tau}_p\,dS
=
-\int_{\partial\mathcal{R}} p\,\mathbf{n}\,dS .
$$ (pgf-pressure-traction)

Using the divergence theorem, this surface integral may be written equivalently as a volume integral of a force per unit volume,

$$
-\int_{\partial\mathcal{R}} p\,\mathbf{n}\,dS
=
-\int_{\mathcal{R}} \nabla p \, dV .
$$ (pgf-pressure-volume)

Interpreting $-\nabla p$ as a force per unit volume, the corresponding pressure acceleration (body force per unit mass) is

$$
\mathbf{b}_{p}
=
-\frac{1}{\rho}\,\nabla p .
$$ (pgf-pressure-accel)

Note that one should **not** apply the divergence theorem directly to $\boldsymbol{\tau}_p/\rho$; the divergence theorem applies to the physical traction $\boldsymbol{\tau}_p$ (force per unit area). The division by $\rho$ to obtain an acceleration is performed after forming the pressure force.

Gravity acts as a body force per unit mass $-\nabla\Phi$, so the combined pressure–gravity acceleration is

$$
\mathbf{b}_{\text{PGF},\perp}
=
-\,\frac{1}{\rho}\,\nabla p \;-\; \nabla \Phi .
$$ (pgf-continuous)

In Omega, the operator $\nabla$ denotes the **horizontal gradient along layers**, following the convention adopted throughout this section. It should not be interpreted as a gradient taken at constant geometric height $z$ or constant pseudo-height $\tilde{z}$.

Under the hydrostatic approximation [](hydrostatic), the form of the horizontal pressure–gravity force in [](pgf-continuous) is valid for *any* choice of vertical coordinate, provided the horizontal gradient is taken consistently with that coordinate. In particular, although the horizontal gradient of $\Phi$ vanishes in $z$-level coordinates in the absence of tides, self-attraction, and loading, the combined expression in [](pgf-continuous) remains the correct horizontal pressure–gravity force when written using the appropriate along-layer gradient operator.

When inserted into the weak, finite-volume form of the momentum equation, $\mathbf{b}_{\text{PGF},\perp}$ represents the combined effect of pressure traction and gravitational body force. In the next subsection, we will layer-average [](pgf-continuous) over the pseudo-height thickness of each layer and combine it with the advective and Coriolis terms to obtain the layered horizontal momentum equation. The resulting horizontal pressure-gradient force in Omega must be a consistent finite-volume discretization of [](pgf-continuous), together with additional metric terms associated with sloping layer interfaces.

### Momentum

We now derive the horizontal momentum equation in our non-Boussinesq, hydrostatic framework, following the same finite-volume approach used for mass and tracer conservation.

As established in the preceding **Pressure Gradient Force** subsection, the combined horizontal pressure–gravity acceleration is
$-\,\rho^{-1}\nabla p - \nabla \Phi$ when written using the appropriate along-layer horizontal gradient operator $\nabla$; in the finite-volume derivation below we retain the same physics by treating gravity as a body force and pressure as a surface traction.

We specify the body forces, including the pressure term, as:

$$
\mathbf{b}_{\perp} = - \left({\bf f} \times \mathbf{u} + \frac{1}{\rho} \nabla p + \nabla \Phi \right).
$$ (body-forces)

The first term on the right hand side is the **Coriolis force**, where $ \mathbf{f} $ is the vector Coriolis vector (e.g., $ f \hat{\mathbf{z}} $ on the sphere).

The second and third terms are $\mathbf{b}_{\text{PGF},\perp}$ from the previous subsection.

Substituting (body-forces) and (surface-forces) into [](#vh-momentum-pseudo), we arrive at:

$$
\frac{d}{dt} \int_A \int_{\tilde{z}_k^{\text{bot}}}^{\tilde{z}_k^{\text{top}}} {\bf u} \, d\tilde{z} \, dA
& + \int_{\partial A} \left( \int_{\tilde{z}_k^{\text{bot}}}^{\tilde{z}_k^{\text{top}}} {\bf u} \otimes {\bf u} \, d\tilde{z} \right) \cdot \mathbf{n}_\perp \, dl \\
& + \int_A \, {\bf u} \left[ \tilde{w}_{tr} - \tilde{\bf u} \right]_{\tilde{z} = \tilde{z}_k^{\text{top}}} \, dA
- \int_A \, {\bf u} \left[ \tilde{w}_{tr} - \tilde{\bf u} \right]_{\tilde{z} = \tilde{z}_k^{\text{bot}}} \, dA \\
& = -\int_A \int_{\tilde{z}_k^{\text{bot}}}^{\tilde{z}_k^{\text{top}}} \left({\bf f} \times \mathbf{u} + \frac{1}{\rho} \nabla p + \nabla \Phi \right) \, d\tilde{z} \, dA.
$$ (vh-momentum-forces)

As with the tracer derivation, we next Reynolds' average [](#vh-momentum-forces),

$$
\frac{d}{dt} \int_A \int_{\tilde{z}_k^{\text{bot}}}^{\tilde{z}_k^{\text{top}}} \left< {\bf u} \right> \, d\tilde{z} \, dA
& + \int_{\partial A} \left( \int_{\tilde{z}_k^{\text{bot}}}^{\tilde{z}_k^{\text{top}}} \left< {\bf u} \otimes {\bf u} \right> \, d\tilde{z} \right) \cdot \mathbf{n}_\perp \, dl \\
& + \int_A \,\left< {\bf u} \left[ \tilde{w}_{tr} - \tilde{\bf u} \right]_{\tilde{z} = \tilde{z}_k^{\text{top}}} \right> \, dA
- \int_A \, \left< {\bf u} \left[ \tilde{w}_{tr} - \tilde{\bf u} \right]_{\tilde{z} = \tilde{z}_k^{\text{bot}}} \right> \, dA \\
& = -\int_A \int_{\tilde{z}_k^{\text{bot}}}^{\tilde{z}_k^{\text{top}}} \, \left< {\bf f} \times \mathbf{u} + \frac{1}{\rho} \nabla p + \nabla \Phi \right> \, d\tilde{z} \, dA.
$$ (vh-momentum-reynolds1)

Here we have also moved the Reynolds' average through the spatial integrals given the properties of the averaging.  Next we do a Reynolds' decomposition, this yields

$$
\frac{d}{dt} \int_A \int_{\tilde{z}_k^{\text{bot}}}^{\tilde{z}_k^{\text{top}}} \left< {\bf u} \right> \, d\tilde{z} \, dA
& + \int_{\partial A} \left( \int_{\tilde{z}_k^{\text{bot}}}^{\tilde{z}_k^{\text{top}}} \left( \left< {\bf u} \right> \otimes \left< {\bf u} \right> + \left< {\bf u}^\prime \otimes {\bf u}^\prime \right> \right) \, d\tilde{z} \right) \cdot \mathbf{n}_\perp dl \, dl \\
& + \int_A \, \left[ \left(\left< {\bf u} \right> \left< \tilde{w}_{tr} \right> + \left<{\bf u}^\prime \tilde{w}_{tr}^\prime \right> \right) - \left(\left<{\bf u}\right> \left<\tilde{\bf u}\right> + \left< {\bf u}^\prime \tilde{\bf u}^\prime \right> \right) \right]_{\tilde{z} = \tilde{z}_k^{\text{top}}} \, dA \\
& - \int_A \, \left[ \left(\left< {\bf u} \right> \left< \tilde{w}_{tr} \right> + \left<{\bf u}^\prime \tilde{w}_{tr}^\prime \right> \right) - \left(\left<{\bf u}\right> \left<\tilde{\bf u}\right> + \left< {\bf u}^\prime \tilde{\bf u}^\prime \right> \right) \right]_{\tilde{z} = \tilde{z}_k^{\text{bot}}} \, dA  \\
& = - \int_A \int_{\tilde{z}_k^{\text{bot}}}^{\tilde{z}_k^{\text{top}}} \, \left( {\bf f} \times \left<\mathbf{u}\right> + \left< \alpha \right> \nabla \left< p \right> + \left< \alpha' \nabla p' \right> + \nabla \left<\Phi\right> \right) \, d\tilde{z} \, dA.
$$ (vh-momentum-reynolds2)

In [](#vh-momentum-reynolds2), we have also used $\alpha = \frac{1}{\rho}$ for notation conciseness.  The definition of the layer average [](#def-layer-average-reynolds) is now utilized on terms with vertical integrals to yield

$$
\frac{d}{dt} \int_A \tilde{h}_k \overline{\left< {\bf u} \right>}^{\tilde{z}}_k  \, dA
& + \int_{\partial A} \tilde{h}_k \left( \overline{\left< {\bf u} \right> \otimes \left< {\bf u} \right>}^{\tilde{z}}_k + \overline{\left< {\bf u}^\prime \otimes {\bf u}^\prime \right>}^{\tilde{z}}_k  \right) \cdot {\bf n}_\perp \, dl \\
& + \int_A \, \left[ \left(\left< {\bf u} \right> \left< \tilde{w}_{tr} \right> + \left<{\bf u}^\prime \tilde{w}_{tr}^\prime \right> \right) - \left(\left<{\bf u}\right> \left<\tilde{\bf u}\right> + \left< {\bf u}^\prime \tilde{\bf u}^\prime \right> \right) \right]_{\tilde{z} = \tilde{z}_k^{\text{top}}} \, dA \\
& - \int_A \, \left[ \left(\left< {\bf u} \right> \left< \tilde{w}_{tr} \right> + \left<{\bf u}^\prime \tilde{w}_{tr}^\prime \right> \right) - \left(\left<{\bf u}\right> \left<\tilde{\bf u}\right> + \left< {\bf u}^\prime \tilde{\bf u}^\prime \right> \right) \right]_{\tilde{z} = \tilde{z}_k^{\text{bot}}} \, dA  \\
& = - \int_A \, \tilde{h}_k \overline{\left( {\bf f} \times \left<\mathbf{u}\right> + \left< \alpha \right> \nabla \left< p \right> + \left< \alpha' \nabla p' \right> + \nabla \left<\Phi\right> \right)}^{\tilde{z}}_k \, dA.
$$ (vh-momentum-reynolds-lay-avg)

The next step is to decompose vertical averages of products.  However we retain the PGF terms as a product because approximating them as products of vertical averages can lead to large inaccuracies. The result is

$$
\frac{d}{dt} \int_{A} \tilde{h}_k \overline{\left< {\bf u} \right>}^{\tilde{z}}_k  \, dA
& + \int_{\partial A} \tilde{h}_k \left( \overline{\left< {\bf u} \right>}^{\tilde{z}}_k \otimes \overline{\left< {\bf u} \right>}^{\tilde{z}}_k + \overline{\delta {\bf u} \otimes \delta {\bf u}}^{\tilde{z}}_k + \overline{\left< {\bf u}^\prime \otimes {\bf u}^\prime \right>}^{\tilde{z}}_k  \right) \cdot {\bf n}_\perp \, dl \\
& + \int_A \, \left[ \left(\left< {\bf u} \right> \left< \tilde{w}_{tr} \right> + \left<{\bf u}^\prime \tilde{w}_{tr}^\prime \right> \right) - \left(\left<{\bf u}\right> \left<\tilde{\bf u}\right> + \left< {\bf u}^\prime \tilde{\bf u}^\prime \right> \right) \right]_{\tilde{z} = \tilde{z}_k^{\text{top}}} \, dA \\
& - \int_A \, \left[ \left(\left< {\bf u} \right> \left< \tilde{w}_{tr} \right> + \left<{\bf u}^\prime \tilde{w}_{tr}^\prime \right> \right) - \left(\left<{\bf u}\right> \left<\tilde{\bf u}\right> + \left< {\bf u}^\prime \tilde{\bf u}^\prime \right> \right) \right]_{\tilde{z} = \tilde{z}_k^{\text{bot}}} \, dA  \\
& = - \int_A \, \tilde{h}_k \overline{\left( {\bf f} \times \left<\mathbf{u}\right> + \left< \alpha \right> \nabla \left< p \right> + \left< \alpha' \nabla p' \right> + \nabla \left<\Phi\right> \right)}^{\tilde{z}}_k \, dA.
$$ (vh-momentum-reynolds-lay-avg2)

We next use the Divergence theorem on the surface integrals and combine terms into fewer integrals on each side of the equation.

$$
\int_A \Bigl\{ \frac{\partial\tilde{h}_k \overline{\left< {\bf u} \right>}^{\tilde{z}}_k}{\partial t}
& + \nabla \cdot \left[ \tilde{h}_k \left(\overline{\left< {\bf u} \right>}^{\tilde{z}}_k \otimes \overline{\left< {\bf u} \right>}^{\tilde{z}}_k + \overline{\delta {\bf u} \otimes \delta {\bf u}}^{\tilde{z}}_k + \overline{\left< {\bf u}^\prime \otimes {\bf u}^\prime \right>}^{\tilde{z}}_k  \right) \right] \\
& + \left[ \left(\left< {\bf u} \right> \left< \tilde{w}_{tr} \right> + \left<{\bf u}^\prime \tilde{w}_{tr}^\prime \right> \right) - \left(\left<{\bf u}\right> \left<\tilde{\bf u}\right> + \left< {\bf u}^\prime \tilde{\bf u}^\prime \right> \right) \right]_{\tilde{z} = \tilde{z}_k^{\text{top}}} \\
& - \left[ \left(\left< {\bf u} \right> \left< \tilde{w}_{tr} \right> + \left<{\bf u}^\prime \tilde{w}_{tr}^\prime \right> \right) - \left(\left<{\bf u}\right> \left<\tilde{\bf u}\right> + \left< {\bf u}^\prime \tilde{\bf u}^\prime \right> \right) \right]_{\tilde{z} = \tilde{z}_k^{\text{bot}}} \Bigr\} \, dA \\
& = - \int_A \left\{ \, \tilde{h}_k \overline{\left( {\bf f} \times \left<\mathbf{u}\right> + \left< \alpha \right> \nabla \left< p \right> + \left< \alpha' \nabla p' \right> + \nabla \left<\Phi\right> \right)}^{\tilde{z}}_k \right\} \, dA.
$$ (vh-momentum-reynolds-lay-avg3)

Since the equation is fully inside the integral, the equation is true for any area and therefore we can write the layer averaged momentum equation as

$$
\frac{\partial\tilde{h}_k \overline{\left< {\bf u} \right>}^{\tilde{z}}_k}{\partial t}
& + \nabla \cdot \left[ \tilde{h}_k \left( \overline{\left< {\bf u} \right>}^{\tilde{z}}_k \otimes \overline{\left< {\bf u} \right>}^{\tilde{z}}_k + \overline{\delta {\bf u} \otimes \delta {\bf u}}^{\tilde{z}}_k + \overline{\left< {\bf u}^\prime \otimes {\bf u}^\prime \right>}^{\tilde{z}}_k  \right) \right] \\
& + \left[ \left(\left< {\bf u} \right> \left< \tilde{w}_{tr} \right> + \left<{\bf u}^\prime \tilde{w}_{tr}^\prime \right> \right) - \left(\left<{\bf u}\right> \left<\tilde{\bf u}\right> + \left< {\bf u}^\prime \tilde{\bf u}^\prime \right> \right) \right]_{\tilde{z} = \tilde{z}_k^{\text{top}}} \\
& - \left[ \left(\left< {\bf u} \right> \left< \tilde{w}_{tr} \right> + \left<{\bf u}^\prime \tilde{w}_{tr}^\prime \right> \right) - \left(\left<{\bf u}\right> \left<\tilde{\bf u}\right> + \left< {\bf u}^\prime \tilde{\bf u}^\prime \right> \right) \right]_{\tilde{z} = \tilde{z}_k^{\text{bot}}} \\
& = - \, \tilde{h}_k \overline{\left( {\bf f} \times \left<\mathbf{u}\right> + \left< \alpha \right> \nabla \left< p \right> + \left< \alpha' \nabla p' \right> + \nabla \left<\Phi\right> \right)}^{\tilde{z}}_k.
$$ (vh-momentum-v1)

The product rule is used on the first two terms of [](#vh-momentum-v1) and then we multiply [](#layer-mass) by $\overline{\mathbf{u}}^{\tilde{z}}_k$ and subtract it from [](#vh-momentum-v1).  This yields

$$
\tilde{h}_k \frac{\partial \overline{\left< {\bf u} \right>}^{\tilde{z}}_k}{\partial t}
& + \tilde{h}_k \overline{\left< {\bf u} \right>}^{\tilde{z}}_k \cdot \nabla \overline{\left< {\bf u} \right>}^{\tilde{z}}_k + \nabla \cdot \left[ \tilde{h}_k \left(\overline{\delta {\bf u} \otimes \delta {\bf u}}^{\tilde{z}}_k + \overline{\left< {\bf u}^\prime \otimes {\bf u}^\prime \right>}^{\tilde{z}}_k  \right) \right] \\
& + \left(\left<\mathbf{u}\right> - \overline{\left<\mathbf{u}\right>}^{\tilde{z}}_k\right) \left\{ \left[ \left<\tilde{w}_{tr}\right> - \left<\tilde{\bf u}\right> \right]_{\tilde{z} = \tilde{z}_k^{\text{top}}} - \left[ \left<\tilde{w}_{tr}\right> - \left<\tilde{\bf u}\right> \right]_{\tilde{z} = \tilde{z}_k^{\text{bot}}} \right\} \\
& + \left\{ \left[ \left<\mathbf{u}^\prime \tilde{w}_{tr}^\prime \right> - \left< \mathbf{u}^\prime \tilde{\bf u}^\prime \right> \right]_{\tilde{z} = \tilde{z}_k^{\text{top}}} - \left[ \left<\mathbf{u}^\prime \tilde{w}_{tr}^\prime \right> - \left< \mathbf{u}^\prime \tilde{\bf u}^\prime \right> \right]_{\tilde{z} = \tilde{z}_k^{\text{bot}}} \right\} \\
& = - \, \tilde{h}_k \overline{\left( {\bf f} \times \left<\mathbf{u}\right> + \left< \alpha \right> \nabla \left< p \right> + \left< \alpha' \nabla p' \right> + \nabla \left<\Phi\right> \right)}^{\tilde{z}}_k.
$$ (vh-momentum-v2)

The previous equation is divided by $\tilde{h}_k$,

$$
\frac{\partial \overline{\left< {\bf u} \right>}^{\tilde{z}}_k}{\partial t}
& + \overline{\left< {\bf u} \right>}^{\tilde{z}}_k \cdot \nabla \overline{\left< {\bf u} \right>}^{\tilde{z}}_k + \frac{1}{\tilde{h}_k} \nabla \cdot \left[ \tilde{h}_k \left(\overline{\delta {\bf u} \otimes \delta {\bf u}}^{\tilde{z}}_k + \overline{\left< {\bf u}^\prime \otimes {\bf u}^\prime \right>}^{\tilde{z}}_k  \right) \right] \\
& + \frac{1}{\tilde{h}_k} \left(\left<\mathbf{u}\right> - \overline{\left<\mathbf{u}\right>}^{\tilde{z}}_k\right) \left\{ \left[ \left<\tilde{w}_{tr}\right> - \left<\tilde{\bf u}\right> \right]_{\tilde{z} = \tilde{z}_k^{\text{top}}} - \left[ \left<\tilde{w}_{tr}\right> - \left<\tilde{\bf u}\right> \right]_{\tilde{z} = \tilde{z}_k^{\text{bot}}} \right\} \\
& + \frac{1}{\tilde{h}_k} \left\{ \left[ \left<\mathbf{u}^\prime \tilde{w}_{tr}^\prime \right> - \left< \mathbf{u}^\prime \tilde{\bf u}^\prime \right> \right]_{\tilde{z} = \tilde{z}_k^{\text{top}}} - \left[ \left<\mathbf{u}^\prime \tilde{w}_{tr}^\prime \right> - \left< \mathbf{u}^\prime \tilde{\bf u}^\prime \right> \right]_{\tilde{z} = \tilde{z}_k^{\text{bot}}} \right\} \\
& = - \, \overline{\left( {\bf f} \times \left<\mathbf{u}\right> + \left< \alpha \right> \nabla \left< p \right> + \left< \alpha' \nabla p' \right> + \nabla \left<\Phi\right> \right)}^{\tilde{z}}_k.
$$ (vh-momentum-v3)


The $\overline{\bf u}^{\tilde{z}}_k \cdot \nabla \overline{\bf u}^{\tilde{z}}_k$ may be replaced with the vector identity

$$
\begin{aligned}
\overline{\left<{\bf u}\right>}^{\tilde{z}}_k \cdot \nabla \overline{\left<{\bf u}\right>}^{\tilde{z}}_k
& = (\nabla \times \overline{\left<{\bf u}\right>}^{\tilde{z}}_k) \times \overline{\left<{\bf u}\right>}^{\tilde{z}}_k + \nabla \frac{\left|\overline{\left<{\bf u}\right>}^{\tilde{z}}_k\right|^2}{2} \\
& = \left( \hat{k} \cdot (\nabla_\perp\times \overline{\left<{\bf u}\right>}^{\tilde{z}}_k)\right)
\left( \hat{k} \times \overline{\left<{\bf u}\right>}^{\tilde{z}}_k \right) + \nabla_\perp\frac{\left|\overline{\left<{\bf u}\right>}^{\tilde{z}}_k\right|^2}{2} \\
& = \zeta {\overline{\left<{\bf u}\right>}^{\tilde{z}}_k}^{\perp} + \nabla K.
\end{aligned}
$$ (advection-identity)

where $\zeta$ is relative vorticity and $K$ is kinetic energy.  This step separates the horizontal advection into non-divergent and non-rotational components, which is useful in the final TRiSK formulation.

The Coriolis term in [](#vh-momentum-v3), when projected into a local coordinate system can be written as

$$
{\bf f} \times \left<\mathbf{u}\right> = f \left<\mathbf{u}\right>^\perp + 2 \Omega \tilde{w}_{tr} \cos \phi
$$ (coriolis-expansion)

We expect the second term above to be negligible for the ocean.  Ocean vertical velocities are largest ($\approx 10$ cm) in the high latitudes, where $\cos \phi$ is small.  Further, these large vertical velocities only arise in non-hydrostatic ocean models.  Therefore we expect this second term to always be small in Omega and will neglect it.

Defining the absolute vorticity $\zeta_a$ as $\zeta + f$, where $f$ is the Coriolis parameter, which results from [](#coriolis-expansion), [](#vh-momentum-v3) becomes,

$$
\frac{\partial \overline{\left< {\bf u} \right>}^{\tilde{z}}_k}{\partial t}
& + \zeta_a {\overline{\left<{\bf u}\right>}^{\tilde{z}}_k}^{\perp} + \nabla K + \frac{1}{\tilde{h}_k} \nabla \cdot \left[ \tilde{h}_k \left(\overline{\delta {\bf u} \otimes \delta {\bf u}}^{\tilde{z}}_k + \overline{\left< {\bf u}^\prime \otimes {\bf u}^\prime \right>}^{\tilde{z}}_k  \right) \right] \\
& + \frac{1}{\tilde{h}_k} \left\{ \left[\left(\left<\mathbf{u}\right> - \overline{\left<\mathbf{u}\right>}^{\tilde{z}}_k\right) \left\{\left<\tilde{w}_{tr}\right> - \left<\tilde{\bf u}\right> \right\} \right]_{\tilde{z} = \tilde{z}_k^{\text{top}}} - \left[  \left(\left<\mathbf{u}\right> - \overline{\left<\mathbf{u}\right>}^{\tilde{z}}_k\right) \left\{\left<\tilde{w}_{tr}\right> - \left<\tilde{\bf u}\right> \right\} \right]_{\tilde{z} = \tilde{z}_k^{\text{bot}}} \right\} \\
& + \frac{1}{\tilde{h}_k}  \left\{ \left[ \left<\mathbf{u}^\prime \tilde{w}_{tr}^\prime \right> - \left< \mathbf{u}^\prime \tilde{\bf u}^\prime \right> \right]_{\tilde{z} = \tilde{z}_k^{\text{top}}} - \left[ \left<\mathbf{u}^\prime \tilde{w}_{tr}^\prime \right> - \left< \mathbf{u}^\prime \tilde{\bf u}^\prime \right> \right]_{\tilde{z} = \tilde{z}_k^{\text{bot}}} \right\} \\
& = - \, \overline{\left< \alpha \right> \nabla \left< p \right> + \left< \alpha' \nabla p' \right> + \nabla \left<\Phi\right>}^{\tilde{z}}_k.
$$ (layer-momentum-final)

The $\delta X$ terms in the layered equations are vertical deviations from the vertical layer average of a given quantity $X$.  We will assume these deviations are small and products of these deviations are even smaller.  Thus we will ignore most of these terms in Omega.  The exception is the pressure gradient force.  At this time we will assume piecewise constant approximations of the vertically continuous system, which is appropriate for the simple pressure gradient force targeted for early versions of Omega.  This assumption will be revisited at a later date.  We also note that these terms could potentially serve as a bridge to multiscale fluxes, as resolution is increased, these $\delta X$ terms would get larger, but likely only for significantly higher resolution (e.g. 10s of meters).  However, these terms would have to be further analyzed and developed as these terms are only deviations from layer averages, not temporal averages as in the Reynolds' approach.

(notational-simplifications)=
### Notational simplifications

Throughout the rest of this document, we will

1. Drop the $< >$ notation around single variables and note that all variables are assumed to be Reynolds' averaged.  Turbulent correlations will retain the angle bracket notation.
2. Change the vertical density weighted average notation from $\overline{\varphi}^{\tilde{z}}_k$ to the more simple $\varphi_k$.  Thus, any subscript $k$ implies a layer average.
3. Define a total vertical velocity across the pseudo height surface as $\tilde{W}_{tr} \equiv \tilde{w}_{tr} - \tilde{u}$.  As a reminder $\tilde{u}$ is the projection of the normal velocity onto the normal vector to the pseudo height surface, in many cases this can be a very small correction to $\tilde{w}_{tr}$.

With these three simplifications, the final continuous equations become

**Tracer:**

$$
\frac{\partial {\tilde h}_k \varphi_k}{\partial t}  & + \nabla \cdot \left({\tilde h}_k \varphi_k {\bf u}_k\right) + \left\{ \left[ \varphi \tilde{W}_{tr} \right]_{{\tilde z}={\tilde z}_k^{\text{top}}} - \left[ \varphi \tilde{W}_{tr} \right]_{{\tilde z}={\tilde z}_{k}^{\text{bot}}} \right\} \\
& = - \nabla \cdot \left({\tilde h}_k \left<\varphi^{\prime} {\bf u}^{\prime}\right>_k\right) \\
& - \left\{ \left[\left< \varphi^\prime {\tilde w}_{tr}^{\prime} \right> - \left< \varphi^\prime \tilde{ u}^{\prime}\right> \right]_{{\tilde z}={\tilde z}_{k}^{\text{top}}} - \left[\left< \varphi^\prime {\tilde w}_{tr}^{\prime} \right> - \left< \varphi^\prime \tilde{ u}^{\prime}\right> \right]_{{\tilde z}={\tilde z}_{k}^{\text{bot}}}\right\}.
$$ (layer-tracer-final-simple)

**Mass:**

$$
\frac{\partial {\tilde h}_k }{\partial t}
+ \nabla \cdot \left({\tilde h}_k {\bf u}_k\right)
+ \left[ \tilde{W}_{tr} \right]_{{\tilde z}={\tilde z}_{k}^{\text{top}}}
- \left[ \tilde{W}_{tr} \right]_{{\tilde z}={\tilde z}_{k}^{\text{bot}}}
= 0.
$$ (layer-mass-final-simple)

**Velocity:**

$$
\frac{\partial {\bf u}_k}{\partial t}
& + \zeta_a {{\bf u}_k}^{\perp} + \nabla K + \frac{1}{\tilde{h}_k} \nabla \cdot \left[ \tilde{h}_k \left< {\bf u}^\prime \otimes {\bf u}^\prime \right>_k \right] \\
& + \frac{1}{\tilde{h}_k} \left\{ \left[\left(\mathbf{u} - \mathbf{u}_k\right) \left\{\tilde{W}_{tr} \right\} \right]_{\tilde{z} = \tilde{z}_k^{\text{top}}} - \left[  \left(\mathbf{u} - \mathbf{u}_k\right) \left\{\tilde{W}_{tr} \right\} \right]_{\tilde{z} = \tilde{z}_k^{\text{bot}}} \right\} \\
& + \frac{1}{\tilde{h}_k}  \left\{ \left[ \left<\mathbf{u}^\prime \tilde{w}_{tr}^\prime \right> - \left< \mathbf{u}^\prime \tilde{ u}^\prime \right> \right]_{\tilde{z} = \tilde{z}_k^{\text{top}}} - \left[ \left<\mathbf{u}^\prime \tilde{w}_{tr}^\prime \right> - \left< \mathbf{u}^\prime \tilde{ u}^\prime \right> \right]_{\tilde{z} = \tilde{z}_k^{\text{bot}}} \right\} \\
& = - \left(\alpha \nabla p + \nabla \Phi \right)_k.
$$ (layer-momentum-final-simple)

**XSAD: why do we drop $\left< \alpha' \nabla p' \right>$??**

In [](#layer-tracer-final-simple) and [](#layer-momentum-final-simple) we have not combined turbulent flux terms (e.g., $\left<\varphi^\prime \tilde{W}_{tr}^\prime \right>$) as the two terms that comprise this term ($\left< \varphi^\prime {\tilde w}_{tr}^{\prime} \right>$ and $\left< \varphi^\prime \tilde{ u}^{\prime}\right>$) are fundamentally different processes that must be modeled separately.

## 10. Discrete Equations

In the following equations, the subscripts $i$, $e$, and $v$ indicate cell, edge, and vertex locations and subscript $k$ is the layer.  Square brackets $[\cdot]_e$ and $[\cdot]_v$ are quantities that are interpolated to edge and vertex locations. For vector quantities, $u_{e,k}$ denotes the normal component at the center of the edge (where the normal component is the local $\mathbf{n}_\perp$ described in above sections), while $u_{e,k}^\perp$ denotes the tangential component. We have switched from $\varphi_{i,k}^{bot}$ to the identical $\varphi_{i,k+1}^{top}$ for all variables in order for the notation to match the array names in the code. It is important to note that any term without a subscript $k$ are quantities evaluated at that location and are ***not*** individual layer averages, but will be reconstructed as a function of neighboring layer average values.  For example, Eq (44) of [White and Adcroft (2008)](https://www.sciencedirect.com/science/article/pii/S0021999108002593) shows a third order reconstruction.  Finally, we note that in the equations below, it is assumed that $k$ increases away from the surface.

**Mass:**

$$
\frac{\partial {\tilde h}_{i,k} }{\partial t}
+ \nabla \cdot \left(\left[{\tilde h}_k\right]_e {\bf u}_{e,k}\right)
+ \left[ \tilde{W}_{tr} \right]_k^\text{top}
- \left[ \tilde{W}_{tr} \right]_{k+1}^\text{top}
= 0.
$$ (discrete-mass)

In this equation, mass source term ($Q$), such as sources like river runoff, sea ice freshwater fluxes, precipitation, and evaporation are bound up in the surface value of $\tilde{W}_{tr}^{\text{top}}$.  The mass sources will be normalized by $\rho_0$ to achieve consistent units.

**Tracer:**

$$
\frac{\partial {\tilde h}_{i,k} \varphi_{i,k}}{\partial t}  & + \nabla \cdot \left(\left[{\tilde h}_{i,k} \varphi_{i,k}\right]_e {\bf u}_{e,k}\right) + \left\{ \left[ \varphi \tilde{W}_{tr} \right]_k^\text{top} - \left[ \varphi \tilde{W}_{tr} \right]_{k+1}^\text{top} \right\} \\
& = - \nabla \cdot \left({\tilde h}_k \left<\varphi^{\prime} {\bf u}^{\prime}\right>_k \right) - \left\{ \left[\left< \varphi^\prime {\tilde w}_{tr}^{\prime} \right> - \left< \varphi^\prime \tilde{ u}^{\prime}\right> \right]_k^\text{top} - \left[\left< \varphi^\prime {\tilde w}_{tr}^{\prime} \right> - \left< \varphi^\prime \tilde{ u}^{\prime}\right> \right]_{k+1}^\text{top} \right\}.
$$ (discrete-tracer)

In the tracer equation, we note that surface fluxes (e.g. latent heat fluxes) will appear in the surface value of the vertical turbulent tracer flux (e.g., $\left[\left<\varphi^\prime \tilde{w}_{tr}^\prime\right>\right]_{k=0}^\text{top}$).  The surface fluxes will also include effects from volumetric surface fluxes that carry a non zero tracer value. The terms on the second line represent the small scale turbulence and the form used in Omega to represent these terms will be discussed in the next section.

**Velocity:**

Omega will only predict the layer average normal velocity, so we drop the bold face on the $u$ terms except for the product of primes, which is specified in the next section.

$$
\frac{\partial u_{e,k}}{\partial t}
& + \left[ {\bf k} \cdot \nabla \times u_{e,k} +f_v\right]_e\left(u_{e,k}^{\perp}\right) + \left[\nabla K\right]_e  \\
& + \frac{1}{\left[\tilde{h}_{i,k}\right]_e} \left\{ \left[\left(u - u_k\right) \left\{\tilde{W}_{tr} \right\} \right]_{e,k}^\text{top} - \left[  \left(u - u_k\right) \left\{\tilde{W}_{tr} \right\} \right]_{e,k+1}^\text{top} \right\} \\
& = - \left(\alpha \nabla p + \nabla \Phi \right)_{e,k} - \frac{1}{\left[\tilde{h}_{i,k}\right]_e} \nabla \cdot \left( \tilde{h}_k \left< {\bf u}^\prime \otimes {\bf u}^\prime \right>_k \right) \\
& - \frac{1}{\left[\tilde{h}_{i,k}\right]_e^\text{top}}  \left\{ \left[ \left<\mathbf{u}^\prime \tilde{w}_{tr}^\prime \right> - \left< \mathbf{u}^\prime \tilde{ u}^\prime \right> \right]_{e,k}^\text{top} - \left[ \left<\mathbf{u}^\prime \tilde{w}_{tr}^\prime \right> - \left< \mathbf{u}^\prime \tilde{ u}^\prime \right> \right]_{e,k+1}^\text{top} \right\}.
$$ (discrete-momentum)

**Diagnostic Relations:**

 $$
 p_{i,k} = p_{i}^{surf} + g\rho_0 \sum_{k'=1}^{k-1} \tilde{h}_{i,k'} + \frac{1}{2} g\rho_0 \tilde{h}_{i,k}
 $$ (discrete-pressure)

$$
\alpha_{i,k} = f_{eos}(p_{i,k},\Theta_{i,k},S_{i,k})
$$ (discrete-eos)

$$
z_{i,k}^{top} = z_{i}^{floor} + \rho_0 \sum_{k'=k}^{K_{max}} \alpha_{i,k'} \tilde{h}_{i,k'}
$$ (discrete-z)

We refer to these as the discrete equations, but time derivatives remain continuous. The time discretization is described in the [time stepping design document](TimeStepping.md). The velocity, mass-thickness, and tracers are solved prognostically using [](discrete-momentum), [](discrete-mass), [](discrete-tracer). At the new time, these variables are used to compute pressure [](discrete-pressure), specific volume [](discrete-eos), and z-locations [](discrete-z). Additional variables are computed diagnostically at the new time: $\mathbf{u}^{\perp}$, $K$, $\zeta_a$, $z^{mid}$, $\Phi$, etc. The initial geopotential is simply $\Phi=gz$, but additional gravitational terms may be added later.

The horizontal operators $\nabla$, $\nabla\cdot$, and $\nabla \times$ are now in their discrete form. In the TRiSK design, gradients ($\nabla$) map cell centers to edges; divergence ($\nabla \cdot$) maps edge quantities to cells; and curl ($\nabla \times$) maps edges to vertices. The exact form of operators and interpolation stencils remain the same as those given in {ref}`Omega-0 operator formulation <33-operator-formulation>` The discrete version of terms common with Omega-0, such as advection, potential vorticity, and $\nabla K$, can be found in {ref}`Omega-0 Momentum Terms <34-momentum-terms>` and {ref}`Omega-0 Thickness and Tracer Terms <35-thickness-and-tracer-terms>`.


## 11. Sub gridscale parameterizations

### Horizontal Velocity Dissipation

There are two terms related to horizontal momentum dissipation in [](#discrete-momentum) that need to be parameterized, $\left<\mathbf{u}^\prime \tilde{u}^\prime \right>_k$ and $\frac{1}{\left[\tilde{h}_{i,k}\right]_e} \nabla \cdot \left( \tilde{h}_k \left< {\bf u}^\prime \otimes {\bf u}^\prime \right>_k \right)$.  The former only arises from the layer integration in pseudo-height coordinates, we interpret this term as the projection of the horizontal momentum dissipation that crosses $\tilde{z}$ interfaces.  Given this, we discuss the form of the horizontal dissipation parameterization first and return to the second term in a later subsection.

As in MPAS-Ocean, parameterization of the horizontal momentum dissipation is through Laplacian or Biharmonic dissipation,

$$
\frac{1}{\left[\tilde{h}_{i,k}\right]_e} \nabla \cdot \left( \tilde{h}_k \left< {\bf u}^\prime \otimes {\bf u}^\prime \right>_k \right) =  \nu_{2,e} \nabla^2 u_{e,k} - \nu_{4,e} \nabla^4 u_{e,k}.
$$ (discrete-mom-diss)

Again, the quantities in [](#discrete-mom-diss) are layer averaged.  The gradient of $\tilde{h}$ is assumed to be small relative to the stress tensor to allow the utilization of traditional parameterization of the dissipation. In [](#discrete-mom-diss), the viscosities are defined as

$$
\nu_{2,e} = \nu_2 \times \text{meshScaling}_{2,e}
$$

and

$$
\nu_{4,e} = \nu_4 \times \text{meshScaling}_{4,e}
$$

Where the $\nu_4$ and $\nu_2$ are the constant Biharmonic and Laplacian coefficients specified in the YAML configuration file.  The meshScaling variable is constructed to appropriately scale the viscosity with changes in mesh spacing (for details, see [Hecht et al 2008](https://www.researchgate.net/profile/Mark-Petersen-2/publication/255570421_Lateral_Mixing_in_the_Eddying_Regime_and_a_New_Broad-Ranging_Formulation/links/55886e1c08ae347f9bda9f04/Lateral-Mixing-in-the-Eddying-Regime-and-a-New-Broad-Ranging-Formulation.pdf)).  This formulation implies that for regionally refined configurations, we cannot strictly pull the $\nu_{2,e}$ and $\nu_{4,e}$ variables through the gradient operator.  To do this, we assume that the spatial gradient of $\nu$ is small locally and we can treat them as constant and pass them through the gradient operator.

#### Laplacian dissipation (del2)

$$
 \nu_{2,e} \nabla^2 u_{e,k} = \nu_{2,e} \left( \nabla D_{i,k} - \nabla^{\perp} \zeta_{v,k} \right)
$$ (discrete-mom-del2)

where $D$ is divergence and $\zeta$ is relative vorticity. See {ref}`Omega V0 Section 3.4.4 <344-del2-momentum-dissipation>` for further details.

#### Biharmonic dissipation (del4)
As in {ref}`Omega V0 Section 3.4.5 <345-del4-momentum-dissipation>`, biharmonic momentum dissipation is computed with two applications of the Del2 operator above.

$$
 - \nu_{4,e} \nabla^4 u_{e,k}
= - \nu_{4,e} \nabla^2 \left( \nabla^2 u_{e,k} \right)
$$ (discrete-mom-del4)

### Velocity dissipation across a sloping $\tilde{z}$ surface

We interpret $\left<\mathbf{u}^\prime \tilde{u}^\prime \right>$ as the dissipation of momentum across the sloping $\tilde{z}$ surface.

$$
\left<\mathbf{u}^\prime \tilde{u}^\prime \right> = \left\{\left[\nu_{2,e} \left( \nabla \tilde{D}_{i} - \nabla^{\perp} \tilde{\zeta}_{v} \right)\right]_k - \left[\nu_{4,e} \nabla^2 \left( \nabla^2 \tilde{u}_{e,k} \right)\right]_k\right\}
$$ (discrete-mom-flux-sloping)

While it looks very similar to [](#discrete-mom-del2) - [](#discrete-mom-del4), there are a few critical differences.  First, the normal velocities in the divergence and relative vorticity in [](#discrete-mom-flux-sloping) are the reconstructed velocity at the top of the cell along an edge, not the layer average.  Second, the velocities in the divergence and relative vorticity are only the projection across the interface (hence the tilde on $D$ and $\zeta$), computed in a discrete sense following

$$
\left[\tilde{u}_{e}\right]_k = u_e \nabla \tilde{z}_{e,k}
$$

in this relation, we have moved the subscript $k$ off the variable itself to prevent confusion with the layer average.  With this definition, [](#discrete-mom-flux-sloping) goes to zero for flat layer surfaces.

#### Vertical velocity dissipation
The vertical turbulent momentum stress is most commonly parameterized as a down-gradient process, i.e.,

$$
\left[ \left<\mathbf{u}^\prime \tilde{w}_{tr}^\prime \right> \right]_{e,k} = -\frac{\nu_v \rho}{\rho_0} \left[\frac{\partial u}{\partial \tilde{z}}\right]_{e,k}
$$

Plugging this relation into the last part of [](#discrete-momentum)

$$
\frac{1}{\left[\tilde{h}_{i,k}\right]_{e,k}}  \left\{ \left[ \nu_v \left[\frac{\partial u}{\partial \tilde{z}}\right]_{e,k} \right]_{e,k} - \left[ \nu_v \left[\frac{\partial u}{\partial \tilde{z}}\right]_{e,k} \right]_{e,k+1}  \right\}
$$ (discrete-mom-vert-diff)

##### Vertical derivatives in a finite volume framework

Since, Omega will predict layer average quantities (e.g., $u_k$), it's not immediately clear that a traditional discretization is appropriate as there is a discontinuity in the predicted variables at layer interfaces.  To circumvent this problem, we turn to weak derivatives instead of traditional pointwise forms.  A weak derivative of $\varphi$ is defined as

$$
\int_{\tilde{z}_{k+1/2}^{\text{top}}}^{\tilde{z}_{k-1/2}^{\text{top}}} \varphi^\prime \psi d\tilde{z} = \varphi \psi - \int_{\tilde{z}_{k+1/2}^{\text{top}}}^{\tilde{z}_{k-1/2}^{\text{top}}} \varphi \psi^\prime d\tilde{z}.
$$ (weak-derivative-defn)

Here $\psi$ is a smooth test function.  If we choose $\psi$ such that it has compact support on the interval $(\tilde{z}_{k+1/2}^{\text{top}},\tilde{z}_{k-1/2}^{\text{top}})$ then the first term on the right hand side of [](#weak-derivative-defn) is zero.  Using the fact that Omega predicts layer averages, the integral on the right hand side can be broken into

$$
-\int_{\tilde{z}_{k+1/2}^{\text{top}}}^{\tilde{z}_{k-1/2}^{\text{top}}} \varphi \psi^\prime d\tilde{z} = - \int_{\tilde{z}_{k+1/2}^{\text{top}}}^{\tilde{z}_{k}^{\text{top}}} \varphi_{k+1} \psi(z)^\prime d\tilde{z} - \int_{\tilde{z}_{k}^{\text{top}}}^{\tilde{z}_{k-1/2}^{\text{top}}} \varphi_k \psi(z)^\prime d\tilde{z}.
$$

Since $\varphi$ is constant over the integral domain, it can be pulled outside the integral and we use the fundamental theorem to yield

$$
-\varphi_{k+1} (\psi(k)-\psi(k+1/2)) - \varphi_k (\psi(k-1/2) - \psi(k)).
$$

Since we have assumed $\psi$ has compact support within the interval but not including the endpoints, $\psi(k-1/2) = \psi(k+1/2) = 0$ and we are left with

$$
\int_{\tilde{z}_{k+1/2}^{\text{top}}}^{\tilde{z}_{k-1/2}^{\text{top}}} \varphi^\prime \psi d\tilde{z} = \psi(k) \left(\varphi_k - \varphi_{k+1}\right).
$$

Since we are able to pick any test function that satisfies our assumptions, we choose the dirac delta function and rewrite the previous equation as

$$
\varphi^\prime(z) = \left(\varphi_k - \varphi_{k+1}\right) \delta \left(z-\tilde{z}_{k}^{\text{top}}\right).
$$ (weak-derivative)

If [](#weak-derivative) is averaged between $\tilde{z}_{k-1/2}^{\text{top}}$ and $\tilde{z}_{k+1/2}^{\text{top}}$, we arrive at the expected form of the gradient centered on layer interfaces

$$
\left[\frac{\partial \varphi}{\partial \tilde{z}}\right]_{avg} = \frac{\left(\varphi_k - \varphi_{k+1}\right)}{0.5 \left(\tilde{h}_k + \tilde{h}_{k+1}\right)}
$$ (weak-deriv-final)

where $0.5 \left(\tilde{h}_k + \tilde{h}_{k+1}\right)$ is the distance between $\tilde{z}_{k-1/2}^{\text{top}}$ and $\tilde{z}_{k+1/2}^{\text{top}}$.  A similar derivation can be followed to compute gradients across layer centers.  This will form discrete derivatives in Omega.

With this, we can now fully discretize [](#discrete-mom-vert-diff) as

$$
-\frac{1}{\left[\tilde{h}_k\right]_e} \left\{ \left[\left<u^\prime \tilde{w}_{tr}^\prime\right> \right]_{e,k} - \left[\left<u^\prime \tilde{w}_{tr}^\prime\right> \right]_{e,k+1} \right\} = -\frac{1}{\left[\tilde{h}_k\right]_e} \left\{ \frac{\left(u_{e,k-1} - u_{e,k}\right)}{0.5 \left(\tilde{h}_{k-1} + \tilde{h}_{k}\right)} - \frac{\left(u_{e,k} - u_{e,k+1}\right)}{0.5 \left(\tilde{h}_k + \tilde{h}_{k+1}\right)} \right\}.
$$ (final-vert-vel-dissipation)

This form can be interfaced with the Omega [tridiagonal solver](TridiagonalSolver.md) routine.

### Forcing at the top and bottom of the ocean

The discretized momentum and tracer forcing appear as the surface value of the vertical turbulent fluxes of tracer and momentum.  We note that surface forcings (e.g. river runoff) could also extend beyond the first layer.  This will be discussed further in a later design document.  Omega also includes a ocean floor vertical turbulent flux of momentum.

#### Wind Forcing

The wind forcing is applied as a top boundary condition during implicit vertical mixing as

$$
\frac{\tau_{e}}{\rho_0 [ \tilde{h}_{i,k}]_e}
$$

where $\tau$ is the wind stress in [Pa]. Since the pseudo-thickness $\tilde{h}$ is in [m], this results in the desired units of [m s$^{-2}$] for a momentum tendency term.

#### Bottom Drag

Bottom Drag is applied as a bottom boundary condition during implicit vertical mixing as

$$
- C_D \frac{u_{e,k}\left|u_{e,k}\right|}{\rho_0[\alpha_{i,k}\tilde{h}_{i,k}]_e}.
$$ (discrete-mom-bottom)

The units of the term in the denominator is length [m], so that the full term has units of [m s$^{-2}$].

#### Rayleigh Drag

Rayleigh drag is a simple drag applied to every cell.  It is used to ensure stability during spin-up.

$$
- Ra \, u_{e,k}
$$ (discrete-mom-Ra)

#### Temperature, salinity, and freshwater forcing

Direct forcing of temperature, e.g. from latent or sensible heat fluxes take a form similar to MPAS-Ocean

$$
\frac{LHF}{C_p \rho_1}
$$

where $\rho_1$ is the density in the top layer of Omega. This gives units of [m K s$^-1$].

#### Freshwater forcing

Since Omega is a non-Boussinesq ocean, surface sources of water will be mass fluxes instead of being converted into thickness fluxes.  Similar to MPAS-Ocean, Omega will include an ability to spread certain fluxes (e.g., river runoff) over a depth specified in the YAML configuration file.

### Horizontal Tracer Diffusion

As with momentum dissipation, the horizontal tracer diffusion arises from the $\left<\mathbf{u}_k^\prime \varphi_{k}^\prime \right>$ and $\left< \tilde{u}^\prime \varphi^\prime \right>$.  As in MPAS-Ocean, the former term can be parameterized either as Laplacian or Biharmonic diffusion,

$$
D^\varphi_{i,k} =  \kappa_{2,e} \nabla^2 \varphi_{i,k} - \kappa_{4,e} \nabla^4 \varphi_{i,k}.
$$ (discrete-tracer-diff)

As in [](#discrete-mom-diss), we have defined $\kappa_{2,e} = \kappa_2 \times \text{meshScaling}_{4,e}$ and $\kappa_{2,e} = \kappa_2 \times \text{meshScaling}_{4,e}$, where $\kappa_2$ and $\kappa_4$ are the constant values specified in the YAML configuration file.  We also assume that $\varphi_{i,k} \nabla^2 \kappa_{2,e}$ is small relative to $\kappa_{2,e} \nabla^2 \varphi_{i,k}$.

#### Laplacian diffusion (del2)
The Laplacian may be written as the divergence of the gradient,

$$
\nabla \cdot \left(\tilde{h}_k \left<\varphi^\prime u^\prime \right>_k \right) = \nabla \cdot \left( \tilde{h}_{i,k} \kappa_{2,e} \nabla \varphi_{i,k} \right).
$$ (discrete-tracer-del2)

See {ref}`Omega V0 Section 3.5.2 <352-del2-tracer-diffusion>` for details of this calculation.

#### Biharmonic diffusion (del4)
The biharmonic is a Laplacian operator applied twice,

$$
 -  \nabla \cdot \left( \kappa_{4,e} \nabla
\left[
\nabla \cdot \left( \tilde{h}_{i,k} \nabla \varphi_{i,k} \right)
\right]
 \right).
$$ (discrete-tracer-del4)

Each of these operators are written as horizontal stencils in the {ref}`Omega V0 Operator Formulation Section <33-operator-formulation>`.  Again we note that the variables in these equations are the layer average.

#### Horizontal tracer diffusion across a sloping surface
As with horizontal momentum dissipation, there is a turbulent flux of tracer across a sloping $\tilde{z}$ interface.  We interpret the $\left< \tilde{u}^\prime \varphi^\prime \right>$ as the projection of the horizontal turbulent flux across the sloping interface.  The form of the diffusion is similar, taking Laplacian diffusion as an example

$$
 \nabla \cdot \left( \tilde{h}_{i} \kappa_{2,e} \nabla \varphi_{i} \right)_k.
$$

While this is similar in form, this uses the reconstruction at the top of the layer and not the layer averages directly as in [](#discrete-tracer-del2).

### Vertical tracer diffusion
The vertical tracer diffusion arises from the $\rho_0\left(\left[\left<\varphi^\prime \tilde{w}_{tr}^\prime \right> \right]_k - \left[\left<\varphi^\prime \tilde{w}_{tr}^\prime \right> \right]_{k+1} \right)$ term. Discretization of this term differs from the vertical velocity diffusion. In the tracer equation, the quantity being updated is the mass-weighted tracer $\tilde{h}_k \varphi$, while in the momentum equation, the variable is simply the velocity $\bf u$, without multiplication by pseudo-thickness $\tilde{h}_k$. This difference affects how vertical mixing is handled. For the momentum equation, the vertical mixing term includes an explicit $1/\tilde{h}_k$ factor, which makes it easy to rearrange the equation into a standard tridiganal form that can be solved implicitly ([](#final-vert-vel-dissipation)). However, in the tracer equation [](#layer-tracer-final-simple), because $\tilde{h}_k$ is included in the prognostic variable $\tilde{h}_k \varphi_k$ and not divided through as in the velocity equation.  Thus isolating $\varphi_k^{n+1}$ becomes more difficult, where $n$ is an index of timestep. In particular, if $\tilde{h}_k^{n+1}$ were treated implicitly, the left-hand side of [](#layer-tracer-final-simple) would involve the product of two unknowns $\tilde{h}_k^{n+1}$ and $\varphi_k^{n+1}$, resulting in a non-linear system. This nonlinearity prevents the formation of the standard tridiagonal matrix for $\varphi_k^{n+1}$, unlike in the vertical momentum diffusion. To address this, Omega follows the approach used in the MPAS-Ocean: the tracer update is split into two steps. First, a provisional tracer field is computed with all processes except vertical mixing (e.g. advection, horizontal mixing). Then, vertical tracer diffusion is applied in a separate implicit step using this provisional field.

During the tracer update, Omega first computes a provisional tracer field that excludes vertical turbulent mixing to obtain $\varphi^{n+1}$. This temporary tracer update is given by

$$
\varphi^{*} = \frac{\tilde{h}_k^{n} \varphi^{n} + \Delta t  (\text{Adv} + \text{Diff}_{\text{H}})}{\tilde{h}_k^{n+1}}.
$$ (provision-tracer-update)

where the term "$\text{Adv}$" represents resolved tracer advection and includes both horizontal transport and resolved vertical transport by the projected vertical velocity $\tilde{W}_{tr}$, and the "$\text{Diff}_{\text{H}}$" term represents horizontal tracer diffusion including the harmonic and biharmonic diffusion.
Once the provisional tracer field ($\varphi^{*}$) is formed, vertical diffusion is applied using an implicit tridiagonal solver. The non mass weighted tracer can now be used to compute the diffusive tendency and compute the final tracer at the next timestep. Given diffusion is done on the non mass weighted tracer, the diffusive term is written as

$$
\frac{\rho_0}{\tilde{h}_k} \left\{\left[\left<\varphi^\prime \tilde{w}_{tr}^\prime \right>\right]_k - \left[\left<\varphi^\prime \tilde{w}_{tr}^\prime \right>\right]_{k+1} \right\}
$$

Taking a traditional down gradient approximation and using [](#weak-deriv-final) the vertical turbulent flux can be applied using a tridiagonal solver ($\mathcal{L}$) in the
[tridiagonal solver](TridiagonalSolver) in the implicit vertical mixing step:
$$
\varphi^{n+1}=\mathcal{L}_\varphi^{-1} (\varphi^*)
$$
Additional enhancements, such as non local tracer fluxes will be discussed in future design documents.


### MPAS-Ocean Equations of Motion

The MPAS-Ocean layered formulation are provided here for reference. MPAS-Ocean solves for momentum, thickness, and tracers at layer $k$. These are continuous in the horizontal and discrete in the vertical.

$$
\frac{\partial {\bf u}_k}{\partial t}
+ \frac{1}{2}\nabla \left| {\bf u}_k \right|^2
+ ( {\bf k} \cdot \nabla \times {\bf u}_k) {\bf u}^\perp_k
+ f{\bf u}^{\perp}_k
+ \frac{w_k^{bot}{\bf u}_k^{bot} - w_k^{top}{\bf u}_k^{top}}{h_k}
=
- \frac{1}{\rho_0}\nabla p_k
- \frac{\rho g}{\rho_0}\nabla z^{mid}_k
+ \nu_h\nabla^2{\bf u}_k
+ \frac{\partial }{\partial z} \left( \nu_v \frac{\partial {\bf u}_k}{\partial z} \right),
$$ (mpaso-continuous-momentum)

$$
\frac{\partial h_k}{\partial t} + \nabla \cdot \left( h_k^e {\bf u}_k \right) + w_k^{bot} - w_k^{top} = 0,
$$ (mpaso-continuous-thickness)

$$
\frac{\partial h_k\varphi_k}{\partial t} + \nabla \cdot \left( h_k^e\varphi_k^e {\bf u}_k \right)
+ \varphi_k^{bot} w_k^{bot} - \varphi_k^{top} w_k^{top}
= \nabla\cdot\left(h_k^e \kappa_h \nabla\varphi_k \right)
+ h_k \frac{\partial }{\partial z} \left( \kappa_v \frac{\partial \varphi_k}{\partial z} \right).
$$ (mpaso-continuous-tracer)

The layer thickness $h$, vertical velocity $w$, pressure $p$, and tracer $\varphi$, are cell-centered quantities, while the horizontal velocity ${\bf u}$ and $e$ superscript are variables interpolated to the cell edges.


## 12. Variable Definitions

Table 1. Definition of variables. Geometric variables may be found in the {ref}`Omega V0 design document, Table 1 <32-variable-definitions>`

| symbol  | name   | units    | location | name in code | notes  |
|---------------------|-----------------------------|----------|-|---------|-------------------------------------------------------|
|$D_{i,k}$   | divergence | s$^{-1}$      | cell | Divergence  |$D=\nabla\cdot\bf u$ |
|$f_v$       | Coriolis parameter| s$^{-1}$      | vertex   | FVertex  |  $f = 2\Omega sin(\phi)$, $\Omega$ rotation rate, $\phi$ latitude|
|$f_{eos}$ | equation of state | -  | any | function call | |
|$g$ | gravitational acceleration | m s$^{-2}$ | constant  | Gravity |
|$\tilde{h}_{i,k}$ | pseudo-thickness | m | cell | LayerThickness | $\tilde{h} = (\rho/\rho_0) h$ |
|$h_{i,k}$ | geometric layer thickness | m | cell | GeometricThickness | |
|$k$ | vertical index |  |
|${\bf k}$ | vertical unit vector |  |
|$K_{min}$ | shallowest active layer |  |
|$K_{max}$ | deepest active layer |  |
|$K_{i,k}$  | kinetic energy    | m$^2$ s$^{-2}$  | cell     | KineticEnergyCell  |$K = \left\| {\bf u} \right\|^2 / 2$ |
|$p_{i,k}$ | pressure | Pa | cell | Pressure | see [](discrete-pressure) |
|$p^{floor}_i$ | bottom pressure | Pa | cell | PFloor | pressure at ocean floor
|$p^{surf}_i$ | surface pressure | Pa | cell | PSurface | due to atm. pressure, sea ice, ice shelves
|$q_{v,k}$ | potential vorticity         | m$^{-1}$ s$^{-1}$    | vertex   | PotentialVorticity  |$q = \left(\zeta+f\right)/\tilde{h}$ |
|$Ra$      | Rayleigh drag coefficient   | s$^{-1}$      | constant |   |  |
|$S_{i,k}$ | salinity | PSU | cell | Salinity | a tracer $\varphi$  |
|$t$       | time    | s        | none     |   |  |
|${\bf u}_k$   | velocity, vector form       | m s$^{-1}$      | - |   |  |
|$u_{e,k}$   | velocity, normal to edge      | m s$^{-1}$      | edge     | NormalVelocity  | |
|$u^\perp_{e,k}$   | velocity, tangential to edge      | m s$^{-1}$      | edge     | TangentialVelocity  |${\bf u}^\perp = {\bf k} \times {\bf u}$|
|$\alpha_{i,k}$ | specific volume | m$^3$ kg$^{-1}$ | cell  | SpecificVolume | $v = 1/\rho$ |
|$\tilde{u}_{i,k}$ | projection of normal velocity across a pseudo height surface | m s$^{-1}$ | cell | | |
|$\tilde{w}_{i,k}$ | vertical velocity across a pseudo height surface | m s$^{-1}$ | cell  | VerticalVelocity | volume transport per m$^2$ |
|$\tilde{w}_{tr\ i,k}$ | net vertical transport through a moving surface | m s$^{-1}$ | cell  | NetVertTransportVelocity | volume transport per m$^2$ |
|$\tilde{W}_{tr\ i,k}$ | total vertical velocity across a pseudo height surface | m s$^{-1}$ | cell | TotalVertTransportVelocity | $\tilde{W}_{tr} \equiv \tilde{w}_{tr} - \tilde{u}$  |
|$\tilde{z}$ | vertical coordinate (pseudo-height) | m | - | | positive upward |
|$\tilde{z}^{top}_{i,k}$ | layer top $\tilde{z}$-location | m | cell | ZTop | see [](discrete-z) |
|$\tilde{z}^{mid}_{i,k}$ | layer mid-depth $\tilde{z}$-location | m | cell | ZMid |
|$\tilde{z}^{surf}_{i}$ | ocean surface $\tilde{z}$-location | m | cell | ZSurface |
|$\tilde{z}^{floor}_{i}$ | ocean floor $\tilde{z}$-location | m | cell | ZFloor |
|${z}^{floor}_{i}$ | ocean floor geometric z-location | m | cell | GeometricZFloor | -bottomDepth from MPAS-Ocean |
|$\zeta_{v,k}$   | relative vorticity| s$^{-1}$      | vertex   |  RelativeVorticity |$\zeta={\bf k} \cdot \left( \nabla \times {\bf u}\right)$ |
|$\zeta_a$ | absolute vorticity ($\zeta + f$) | s$^{-1}$ | vertex | |
|$\Theta_{i,k}$ | conservative temperature | C | cell  | Temperature  | a tracer $\varphi$ |
|$\kappa_2$| tracer diffusion  | m$^2$ s$^{-1}$    | cell     |   |  |
|$\kappa_4$| biharmonic tracer diffusion | m$^4$ s$^{-1}$    | cell     |   |  |
|$\kappa_v$| vertical tracer diffusion | m$^2$ s$^{-1}$    | cell     |   |  |
|meshScaling  | variable that holds the scaling factor for biharmonic and laplacian mixing        | unitless   | edge     |   | |
|$\nu_{2,e}$   | horizontal del2 viscosity scaled by mesh resolution        | m$^2$ s$^{-1}$    | edge     |   | |
|$\nu_{4,e}$   | horizontal biharmonic (del4) viscosity scaled by mesh resolution       | m$^4$ s$^{-1}$    | edge     |   |  |
|$\nu_v$| vertical momentum diffusion | m$^2$ s$^{-1}$    | edge       |   |  |
|$\varphi_{i,k}$ | tracer | kg m$^{-3}$ or similar | cell | | e.g. $\Theta$, $S$ |
|$\rho_{i,k}$ | density | kg m$^{-3}$ | cell  | Density |
|$\rho_0$ | Reference density | kg m$^{-3}$ | |  constant |
|$\tau_i$ | wind stress | Pa=N m$^{-2}$ | edge |  SurfaceStress |
|$\Phi_{i,k}$ | geopotential| m$^2$ s$^{-2}$  | cell | Geopotential |$\partial \Phi / \partial z = g$ for gravity |
|$\omega$   | mass transport | kg s$^{-1}$ m$^{-2}$      | cell | VerticalTransport |$\omega=\rho_0 w$|


## 13. Verification and Testing

Capability and testing are similar to [Petersen et al. 2015](http://www.sciencedirect.com/science/article/pii/S1463500314001796). The following tests are in idealized domains and do not require surface fluxes or surface restoring. For the following tests to show results comparable to those published with other models, the full dynamic sequence of density, pressure, momentum, and advection must work correctly. The successful completion of the following tests is a validation of the primitive equation functions in Omega 1.0. All of the following tests may exercise a linear equation of state or the nonlinear TEOS10. The first four tests quantify the anomalous mixing caused by the numerical schemes. The first five are on cartesian planes with regular hexagon meshes.

### Lock Exchange (Optional)
The Lock Exchange is the simplest possible test of a primitive equation model. There is an analytic formulation for the wave propagation speed. It is listed as optional because the Overflow tests the same dynamics.
See [Petersen et al. 2015](http://www.sciencedirect.com/science/article/pii/S1463500314001796) and the compass `lock_exchange` case.

### Overflow
The Overflow test case adds bathymetry to the Lock Exchange. It is a particularly effective test of vertical mass and tracer advection, and vertical mixing. It is useful to compare different vertical coordinates, like level (z- or p-level) versus terrain-following (sigma).
See [Petersen et al. 2015](http://www.sciencedirect.com/science/article/pii/S1463500314001796) and the compass `overflow` case.


### Internal Gravity Wave
The internal gravity wave tests horizontal and vertical advection.
See [Petersen et al. 2015](http://www.sciencedirect.com/science/article/pii/S1463500314001796) and the `internal_wave` case in both compass and polaris.

### Baroclinic Channel
This is the first test to add the Coriolis force and uses a three-dimensional domain. It is designed to result in an eddying simulation at sufficiently high resolution. This tests the combination of Coriolis and pressure gradient forces that produce geostrophic balance, as well as horizontal advection and dissipation for numerical stability.
See [Petersen et al. 2015](http://www.sciencedirect.com/science/article/pii/S1463500314001796) and the `baroclinic_channel` case in both compass and polaris.

### Seamount with zero velocity.
This is a 3D domain with a seamount in the center, where temperature and salinity are stratified in the vertical and constant in the horizontal. The test is simply that an initial velocity field of zero remains zero. For z-level layers the velocity trivially remains zero because the horizontal pressure gradient is zero. For tilted layers, this is a test of the pressure gradient error and the velocity is never exactly zero. This is a common test for sigma-coordinate models like ROMS because the bottom layers are extremely tilted along the seamount, but it is a good test for any model with tilted layers. Omega will use slightly tilted layers in p-star mode (pressure layers oscillating with SSH) and severely tilted layers below ice shelves, just like MPAS-Ocean. See [Ezer et al. 2002](https://www.sciencedirect.com/science/article/pii/S1463500302000033), [Haidvogel et al. 1993](https://journals.ametsoc.org/view/journals/phoc/23/11/1520-0485_1993_023_2373_nsofaa_2_0_co_2.xml), [Shchepetkin and McWilliams 2003](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2001JC001047), and previous MPAS-Ocean [confluence page](https://acme-climate.atlassian.net/wiki/spaces/OCNICE/blog/2015/11/19/40501447/MPAS-O+Sigma+coordinate+test+sea+mount).

### Cosine Bell on the Sphere
This test uses a fixed horizontal velocity field to test horizontal tracer advection. It is repeated from [Omega-0 design document](OmegaV0ShallowWater) and is important to conduct again as we convert Omega to a layered primitive-equation model. See `cosine_bell` case in both compass and polaris.

### Merry-Go-Round
This is an exact test for horizontal and vertical tracer advection. A fixed velocity field is provided, and a tracer distribution is advected around a vertical plane. See the `merry_go_round` test in compass, and the results on the [merry-go-round pull request](https://github.com/MPAS-Dev/compass/pull/108) and [compass port pull request](https://github.com/MPAS-Dev/compass/pull/452).


## References
This section is for references without webpage links. These are mostly textbooks.

- Cushman‐Roisin, B., & Beckers, J.M. (2011). Introduction to Geophysical Fluid Dynamics: Physical and Numerical Aspects. Academic Press.
- Gill, A. E. (2016). Atmosphere—Ocean dynamics. Elsevier.
- Kundu, P.K., Cohen, I.M., Dowling D.R. (2016) Fluid Mechanics 6th Edition, Academic Press.
- Pedlosky, J. (1987). Geophysical Fluid Dynamics (Vol. 710). Springer.
- Vallis, G. K. (2017). Atmospheric and oceanic fluid dynamics. Cambridge University Press.
