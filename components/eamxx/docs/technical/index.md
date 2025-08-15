# EAMxx Technical Guide

The goal of this document is to describe the specific equations,
parameterizations, and numerical methods used in the current version of EAMxx.
Because our master-branch implementation changes every time we make a new
commit, this documentation will also evolve continuously. As such,
documentation for master should always be considered to be preliminary and
under construction. If you want trustworthy documentation, pull it from an
official model release.

## Overview

Currently, EAMxx is only configured for km-scale convection-permitting runs. In
order to provide scientifically-credible simulations at lower resolution,
parameterizations for the following processes would be needed:

1. Deep Convection
2. Gravity-Wave Drag
3. Energy Fixer

The only configuration of EAMxx that is currently implemented is the
convection-permitting version, commonly known as the Simple Cloud-Resolving
E3SM Atmosphere Model (SCREAM).

Processes in EAMxx-SCREAM are:

1. a non-hydrostatic version of the **spectral-element dynamical core**
used by other E3SM Atmosphere Model versions[@Taylor_et20] with
semi-Lagrangian tracer advection as described by
Bradley et al. (2022).[@Bradley_et22]
2. **turbulent mountain stress** is crudely parameterized following
Fiedler and Panofsky (1972) [@Fiedler_Panofsky72] to reduce excessive winds
around topography.
3. the **Simple Higher-Order Closure (SHOC)** parameterization from
Bogenschutz and Krueger (2013) [@Bogenschutz_Krueger13],
which handles turbulent diffusion, condensation/evaporation, and liquid cloud fraction.
4. an **all-or-nothing ice cloud fraction** parameterization that sets ice
cloud fraction to 100% whenever cloud ice mass $q_i$ is less than a
user-specified threshold set by default to 1e-5 kg/kg.
This scheme also sets the total cloud fraction (used by microphysics)
to the maximum of the liquid and ice cloud fraction.
5. the effects of aerosol are prescribed via the
**Simple Prescribed Aerosol (SPA)** scheme, which is very similar to
MACv2-SP.[@Stevens_et17]
6. the **P3 microphysics** scheme from Morrison and Milbrandt (2015)
[@Morrison_Milbrandt15] modified as described by Caldwell et al. (2021)
[@Caldwell_et21] to assume instantaneous liquid saturation adjustment for
consistency with SHOC.
7. **RTE/RRTMGP radiation** from Pincus et al. (2019) [@Pincus_et19]
rewritten in C++ for consistency and performance.
8. the **CFMIP Observation Simulator Package (COSP)** is also integrated into
EAMxx, but currently only the ISCCP output is enabled.

By default processes are called in this order, but which processes to include
and in what order is modifiable at run time.
After all atmospheric processes are called, output is written.
Surface components are then called before the next atmosphere step starts.
These processes are described in more detail in
Caldwell et al. (2021) [@Caldwell_et21].
As in EAM, dynamics operates on a spectral element grid and all other
processes use a finite-volume grid that divides each spectral element into 4 quadrilaterals.
This physics grid is described in Hannah et al. (2021) [@Hannah_et21].
