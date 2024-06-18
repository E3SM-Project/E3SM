# Time Stepping

<!--- use table of contents if desired for longer documents  -->
**Table of Contents**
1. [Overview](#1-overview)
2. [Requirements](#2-requirements)
3. [Algorithmic Formulation](#3-algorithmic-formulation)
4. [Design](#4-design)
5. [Verification and Testing](#5-verification-and-testing)

## 1 Overview

The solution methods for a partial differential equation come in two categories: the spatial discretizations and the temporal discretizations. Here we present the methods of solving the time derivative, $\partial / \partial t$, in the Omega-0 shallow water equations. The spatial discretizations are presented in the [Omega V0: Shallow Water](OmegaV0ShallowWater.md) design document.

## 2 Requirements

### 2.1 Requirement: Accurate discretization of time derivative

The time discretization scheme must be able to advance the shallow water equations with at least second-order accuracy for standard shallow water test cases.

### 2.2 Requirement: Stability

The scheme must be stable for reasonably long time steps compared to other schemes.

### 2.3 Requirement: Performance

The time stepping scheme must follow the scaling and performance criteria outlined in the [Omega V0: Shallow Water](OmegaV0ShallowWater.md) design document. It will be tested with metrics for:

1. Single CPU throughput
1. Parallel CPU scalability to high node counts
1. Single GPU throughput

The first two items will be compared against MPAS-Ocean with the same test case and resolution.

### 2.4 Requirement: Modularity

The time stepping code must be written in such a way that other time stepping schemes can be chosen, and the code is organized in a clear and readable fashion.

### 2.5 Requirement: Conservation

The total volume (area integrated layer thickness $h$) and total tracer (area integrated $\phi$) must be conserved. This means that the totals should remain constant to machine precision in time.


## 3 Algorithmic Formulation

### 3.1 Forward-Backward Scheme

The primary time stepping scheme will be forward-backward. This scheme is commonly used in shallow water models and the barotropic mode of layered models. It suits the design requirements because it is second order, is stable for reasonably long time steps compared to competing schemes, and may be implemented in a performant manner.

Given the shallow water equations introduced in the [Omega V0: Shallow Water](OmegaV0ShallowWater.md) design document, we may consolidate the right-hand-side tendency terms to express them as

$$
\frac{\partial \boldsymbol{u}}{\partial t} = \mathcal{RHS}_u \left( \boldsymbol{u},h \right),
$$

$$
\frac{\partial h}{\partial t} = \mathcal{RHS}_h \left( \boldsymbol{u},h \right),
$$

$$
\frac{\partial h \phi}{\partial t} = \mathcal{RHS}_\phi \left( \boldsymbol{u},h,\phi \right).
$$

The forward-backward scheme is simply a single forward time step, but always using the most recent information available:

$$
h^{n+1} = h^n + \mathcal{RHS}_h \left( \boldsymbol{u}^n,h^n \right) \Delta t
$$

$$
\phi^{n+1} = \frac{1}{h^{n+1}} \left( \phi^n h^n + \mathcal{RHS}_\phi \left( \boldsymbol{u}^n,h^n,\phi^n \right) \Delta t \right)
$$

$$
\boldsymbol{u}^{n+1} = \boldsymbol{u}^n + \mathcal{RHS}_u \left( \boldsymbol{u}^n,h^{n+1} \right) \Delta t
$$
Here the variables are discretized in time only. Spatial discretization can be shown with a subscript if desired, using $u^n_e, h^n_i, \phi^n_i$ etc.
The time domain is discretized into steps $t_0, t_1, ... t_n, t_{n+1}$ where $t_{n+1} = t_n+\Delta t$. The superscript indicates the time on all variables. For example, $h^n = h\left( t_n \right)$.


Note that the tracer equation is thickness-weighted, but we divide by the thickness in order to obtain the tracer variable $\phi$, which has units of concentration.
One might argue that the tracer equation could use $h^{n+1}$ in the right-hand-side calculation, because that quantity is available from the previous computation. However, that would violate the consistency of the solution method of the thickness and tracer equations, and cause the globally-integrated tracer to not be conserved. Like the continuous equations, when we set $\phi=1$ the numerical method for $\phi h$ must reduce to the numerical method for $h$.


### 3.2 Other Schemes

If time permits, other schemes will be added and tested. This includes the four-stage Runge-Kutta (RK4) for direct comparison to MPAS-Ocean. Another possibility is second-order Adams-Bashforth (AB2) which was implemented in MPAS-Ocean. Including multiple time-stepping schemes allows developers to test the modularity requirement above.

## 4 Design

Design specifics will be added at a later time.

## 5 Verification and Testing

The timestepping is an intrinsic part of the Omega-0 model, so the verification and testing are identical to those in the [Omega V0: Shallow Water](OmegaV0ShallowWater.md) design document.
