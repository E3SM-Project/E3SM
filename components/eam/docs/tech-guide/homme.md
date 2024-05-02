# High-Order Methods Modeling Environment

## Overview

EAM uses the a dynamical core (dycore) from the High Order Method Modeling Environment (HOMME).[@taylor_compatible_2010,@guba_spectral_2014,@taylor_energy_2020] The EAM dycore solves the atmospheric primitive equations governing the evolution of velocity, density, pressure and temperature, as well as the transport of water species and related hydrometers, aerosols and other atmospheric constituents.  The governing equations are written in a vertically lagrangian terrain following mass coordinate.  They are discretized with second order finite differences in the radial direction and spectral finite elements in the horizontal (surface of the sphere) directions, and advanced in time with a 3rd order accurate 5 stage Runge-Kutta method.  Dissipation is added through the use of monotoncity contraints on some advectiion terms, explicitly added hyperviscosity, and a Laplacian-based sponge layer in the first few layers at the model top.  The transported species  makes use of an efficient interpolatory semi-Lagrangian method.  EAMv3 uses 80 layers in the vertical.  The use of the spectral finite element method allows EAMv3 to run on fully unstructured grids, including the cubed-sphere grid ([SE Atmosphere Grid Overview (EAM & CAM)](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/34113147)) which provides quasi-uniform resolution over the globe, and regionally refined meshes (RRM) which enhance horizontal resolution in regions of interest ([Library of Regionally-Refined Model (RRM) Grids](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/3690397775)).

## Namelist parameters

Many dynamical core parameters can not be changed independently.  For example, increasing the hyperviscosity coefficient may require reducing the hyperviscosity timestep.   Dycore timesteps are tuned for each resolution and the defaults are close to the CFL stability limit. For complete details, as well as their interactions, see [EAM's HOMME dycore](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/1044644202/EAM+s+HOMME+dycore).

[HOMME Namelist Parameters](../user-guide/namelist_parameters.md#homme)
