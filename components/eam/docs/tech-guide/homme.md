# High-Order Methods Modeling Environment

## Overview

EAM using the a dynamical core (dycore)  from the High Order Method Modeling Environment (HOMME).  The EAM dycore solves the atmospheric primitive equations governing the evolution of velocity, density, pressure and temperature, as well as the transport of water species and related hydrometers, aerosols and other atmospheric constituents.  The governing equations are written in a vertically lagrangian terrain following mass coordinate.  They are discretized with second order finite differences in the radial direction and spectral finite elements in the horizontal (surface of the sphere) directions, and advanced in time with a 3rd order accurate 5 stage Runge-Kutta method.   Dissipation is added through the use of monotoncity contraints on some advectiion terms, explicitly added hyperviscosity, and a Laplacian-based sponge layer in the first few layers at the model top.  The transported species  makes use of an efficient interpolatory semi-Lagrangian method.   EAMv3 uses 80 layers in the vertical.   The use of the spectral finite element method allows EAMv3 to run on fully unstructured grids, including the cubed-sphere grid ([SE Atmosphere Grid Overview (EAM & CAM)](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/34113147)) which provides quasi-uniform resolution over the globe, and regionally refined meshes (RRM) which enhance horizontal resolution in regions of interest ([Library of Regionally-Refined Model (RRM) Grids](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/3690397775)).

### References

- Taylor et al. (2020) [@taylor_energy_2020]
- Bradley et al. (2022) [@bradley_islet_2022]
- Guba et al. (2014) [@guba_spectral_2014]
- Taylor and Fournier (2010) [@taylor_compatible_2010]

## Namelist parameters

Many dynamical core parameters can not be changed independently.  For example, increasing the hyperviscosity coefficient may require reducing the hyperviscosity timestep.   Dycore timesteps are tuned for each resolution and the defaults are close to the CFL stability limit. For complete details, as well as their interactions, see [EAM's HOMME dycore](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/1044644202/EAM+s+HOMME+dycore).

| Parameter        | Description                                                                                 | Default value  |
| ---------------- | ------------------------------------------------------------------------------------------- | -------------- |
| `se_tstep`       | Main dycore timestep. Additional parameters control the hyper viscsosity, trancer and vertical remap timesteps, which are derived from se_tstep. <!-- markdownlint-disable MD033 --><br> units = seconds | Scales linearly with horizontal resolution. <br> NE30 default: `300` |
| `nu`             | Tensor hyperviscosity coefficient, independent of spatial resolution. <br> units = 1/s      | `3.4e-8` |
| `nu_top`         | Scalar viscosity at model top. <br> units = m^2/s                                           | Horizontal resolution dependent <br> NE30 default: `2.5e5` |
| `transport_alg`  | Select between semi-lagrangian and Eulerian based transport schemes                         | `12` = semi-lagranian method with monotinicity and mass preservation |
| `statefreq`      | print a varieity of dycore metrics to the atm.log file every “statefreq” timesteps          | `480`          |
| `vert_remap_alg` | Algorithm used to remap the vertically lagrangian levels back to the reference levels       | `10` = strict monotonicity applied on top of a 2nd order accurate PPM method  |
| `se_ftype`       | Controls how physics tendencies are applied.  0=”dribbled” in during dynamics timesteps.  1=”hard adjustment” after each physics timestep.  2=hybrid approach: hard adjustment for tracers, dribbled for remaining tendencies | `2`          |
