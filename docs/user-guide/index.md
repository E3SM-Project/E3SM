# User Guide

E3SM is not just one climate model but a modeling system that allows
many different configurations of atmosphere, ocean, land and other
components with both full model and data model options. Also, the configurations of model
components can run at different resolutions.  Some configurations
can run easily on a laptop.  Other's require the most powerful
supercomputers.

Using the model requires fist deciding what configuration to use and creating
a case with that configuration.
The configuration options are managed by the Case Control System (CCS)
within the
[Community Infrastructure for Modeling the Earth](https://esmci.github.io/cime/versions/master/html/what_cime/index.html) (CIME).

Before reading the rest of this guide, you should become familiar with
cases, compsets and grids by reading the
[Case Control System Basic Usage](https://esmci.github.io/cime/versions/master/html/users_guide/index.html#case-control-system-part-1-basic-usage)

A [step-by-step guide](https://docs.e3sm.org/running-e3sm-guide/) for running E3SM with a run script
is available.

## Supported Coupled Compsets

A *fully coupled* compset is one which has active components for at least the atmosphere, ocean, land surface, ocean and
sea-ice all interacting.  The fully coupled compsets supported in this version of E3SM are:

| Compset Alias  |   TIME     | ATM    | LND | Sea ice     | Ocean      | River |
| -------- |         ----     | ---    | --- | ---------   | ------     | ----- |
|WCYCL1850          | 1850SOI | EAM    | ELM | MPAS-Seaice | MPAS-Ocean | MOSART |
|WCYCL1850-1pcCO2   | 1850SOI | EAM    | ELM | MPAS-Seaice | MPAS-Ocean | MOSART |
|WCYCL1850-4xCO2    | 1850SOI | EAM    | ELM | MPAS-Seaice | MPAS-Ocean | MOSART |
|WCYCL20TR          | 20TRSOI | EAM    | ELM | MPAS-Seaice | MPAS-Ocean | MOSART |
|WCYCL1950          | 1950SOI | EAM    | ELM | MPAS-Seaice | MPAS-Ocean | MOSART |

Additional compsets are available for the component models. See the User Guides for the components for more information.

## Supported resolution

Only one grid set is supported for the above compsets.

`ne30pg2_r05_IcoswISC30E3r5` - For this grid set, the atmosphere is on the ne30pg2 cubed-sphere mesh with approximately 100km
resolution, the land and river
are on a 0.5deg x 0.5deg structured grid, and the ocean and sea ice are on a hexagonal mesh dervied from the dual of a 30km
resolution icosohedral mesh with ice shelf cavities (wISC) around Antarctica.

