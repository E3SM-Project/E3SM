# EAMxx Technical Guide

The goal of this document is to describe the specific equations, parameterizations, and numerical methods used in the current version of EAMxx. Because our master-branch implementation changes every time we make a new commit, this documentation will also evolve continuously. As such, documentation for master should always be considered to be preliminary and under construction. If you want trustworthy documentation, pull it from an official model release. 

## Overview

Currently, EAMxx is only configured for km-scale convection-permitting runs. In order to provide scientifically-credible simulations at lower resolution, parameterizations for the following processes would be needed:

1. deep convection
2. gravity-wave drag
3. energy fixer

The only configuration of EAMxx that is currently implemented is the convection-permitting version, commonly known as the Simple Cloud-Resolving E3SM Atmosphere Model (SCREAM). Parameterizations in EAMxx-SCREAM are:

1. a non-hydrostatic version of the **spectral-element dynamical core** used by other E3SM Atmosphere Model versions [@Taylor_et20] with semi-Lagrangian tracer advection as described by [@Bradley_et22]
2. the **Simple Higher-Order Closure (SHOC)** parameterization from [@Bogenschutz_Krueger13], which handles turbulent diffusion, condensation/evaporation, and liquid cloud fraction
3. an **all-or-nothing ice cloud fraction** parameterization that sets ice cloud fraction to 100% whenever cloud ice mass q<sub>i</sub> is less than a user-specified threshold set by default to 1e-5 kg/kg. This scheme also sets the total cloud fraction (used by microphysics) to the maximum of the liquid and ice cloud fraction.
4. the **P3 microphysics** scheme from [@Morrison_Milbrandt15] modified as described by [@Caldwell_et21] to assume instantaneous liquid saturation adjustment for consistency with SHOC
5. **RTE/RRTMGP radiation** from [@Pincus_et19] rewritten in C++ for consistency and performance

These parameterizations are described in more detail in Caldwell et al 2021. The CFMIP Observation Simulator Package (COSP) is also integrated into EAMxx, but currently only the ISCCP output is enabled. As in EAM, dynamics operates on a spectral element grid and all other processes use a finite-volume grid that divides each spectral element into 4 quadrilaterals. This physics grid is described in [@Hannah_et21].
