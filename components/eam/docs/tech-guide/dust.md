# Dust aerosol

## Overview

Dust-related processes are represented in the E3SM atmosphere and land model components. In E3SMv3, dust deposition will also be coupled with the sea ice and ocean biogeochemistry in the ocean model. Total emission fluxes of dust particles are calculated at each model time step. A new dust emission scheme following Kok et al. (2014)[@kok_improved_2014] is implemented to E3SMv3, replacing the Zender scheme (Zender et al., 2003)[@zender_mineral_2003] in E3SMv1 and v2 as the default option. The new dust emission scheme includes a time-varying soil erodibility factor for dust mobilization, and includes dust sources in high latitudes (e.g., >60 degree N). A manuscript by Feng et al. is in prep to document the performance of the new emission scheme on dust life cycle and radiative effects in E3SMv3. Dust aerosol is represented in the accumulation and coarse aerosol modes of the MAM4 module following emission. Other dust properties such as optical properties and size distribution at emission are documented in Feng et al. (2022).[@feng_global_2022]

## Namelist parameters

[Dust Namelist Parameters](../user-guide/namelist_parameters.md#dust-aerosol)
