# User Guide

E3SM is not just one earth system model but a modeling system that allows
many different configurations of atmosphere, ocean, land and other
components with both full model and data model options. Also, the configurations of model
components can run at different resolutions.  Some configurations
can run easily on a laptop.  Other's require the most powerful
supercomputers.

Using the model requires first deciding what configuration to use and creating
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
sea-ice all interacting.  Each compset is associated with a specific forcing condition.
Coupled compsets in E3SM are developed  for three science-driven simulation campaigns:  **water cycle change and impacts**, **human-earth system feedbacks**, and **polar processes, sea-level rise and coastal impacts**. The standard coupled configurations -- which consist of prognostic atmosphere, land, river, ocean and sea-ice components -- form the base physical coupled system and are mainly designed for `water cycle change and impacts` simulation campaign.
Below list the standard configuration compsets supported in the current version of E3SM:

|Compset alias | Description |
|:-----------  |:----------- |
|`WCYCL1850` | Standard configuration with pre-industrial climatological forcings |
|`WCYCL1850-4xCO2` | Same as `WCYCL1850` except with abrupt (then persistent) 4xCO2 forcing. |
|`WCYCL1850-1pctCO2` | Same as `WCYCL1850` except with 1 percent per year increase of CO2 concentration |
|`WCYCL1950` | Standard configuration with perpetual 1950 forcings |
|`WCYCL20TR` | Standard configuration with prescribed transient forcings over the historical period (1850-2014) |
|`WCYCLSSP245` | Standard configuration with prescribed SSP-245 forcings |
|`WCYCLSSP370` | Standard configuration with prescribed SSP-370 forcings |
|`WCYCLSSP585` | Standard configuration with prescribed SSP-585 forcings |
| | **Compsets for single-forcings simulations of the historical period (1850-2014)** |
|`WCYCL20TR-GHG` | Greenhouse gases only configuration (`GHGs`)|
|`WCYCL20TR-aer` | Anthropogenic aerosols and precursors only configuration (`aer`)|
|`WCYCL20TR-xGHG-xaer` | All forcings except GHGs and aer (`solar irradiance, stratospheric ozone and volcanic emissions, land use land cover`) |
|`WCYCL20TR-xaer` | All forcings except aer (`GHGs, solar irradiance, stratospheric ozone and volcanic emissions, land use land cover`) |
|`WCYCL20TR-nat` | Natural-only configuration (`solar irradiance, stratospheric volcanic emissions, land use land cover`) |
|`WCYCL20TR-ozone` | Stratospheric ozone only configuration |
|`WCYCL20TR-lulc` | Land use land cover only configuration |
|`WCYCL20TR-volc` | Stratospheric volcanic emissions only configuration |

The table below specifies which forcing category is fixed at 1850 conditions and which are allowed to vary over the historical period
for each of the historical ("20TR") compsets including the single-forcing compsets.

|Compset alias  |  GHGs    | Aerosol & precursors | Oxidants | Ozone (CI & Linoz) | Volcano | Solar   |  Land Use & ndep/popdensa |
|:------------  |:-----:   | :---:                | :---:    | :---:              | :---:   | :---:   | :---:   |
|`WCYCL20TR`    | varying       | varying         | varying  | varying            | varying | varying | varying |
|`WCYCL20TR-GHG`| varying       | 1850            | 1850     | 1850               | 1850    | 1850    | 1850    |
|`WCYCL20TR-aer`| 1850          | varying         | varying  | 1850               | 1850    | 1850    | 1850    |
|`WCYCL20TR-xGHG-xaer`| 1850    | 1850            | 1850     | varying            | varying | varying | varying |
|`WCYCL20TR-xaer`| varying      | 1850            | 1850     | varying            | varying | varying | varying |
|`WCYCL20TR-nat`| 1850          | 1850            | 1850     | 1850               | varying | varying | 1850    |
|`WCYCL20TR-ozone`| 1850        | 1850            | 1850     | varying            | 1850    | 1850    | 1850    |
|`WCYCL20TR-lulc`| 1850         | 1850            | 1850     | 1850               | 1850    | 1850    | varying |
|`WCYCL20TR-volc`| 1850         | 1850            | 1850     | 1850               | varying | 1850    | 1850    |

<!-- markdownlint-disable-next-line MD033 -->
- <sub> *Volcano* refers to stratospheric volcanic SO2 emissions; *1850* for *Volcano* refers to background (average) stratospheric volcanic emissions used in pre-industrial control experiments</sub>
<!-- markdownlint-disable-next-line MD033 -->
- <sub> *Oxidants* always follow *Aerosol & precursors* for fixed or varying.</sub>

The compsets for the other two science simulation campaigns are being finalized, with additional components and/or features.
The compset naming follows the same convention, e.g., `CRYO1850` and `CRYO1850-4xCO2` are with prognostic ice-shelf melt fluxes for the `polar processes` simulation campaign.

Compsets are also available for standalone component model configurations, See the User Guides for the components for more information.

## Supported resolution

Currently two grid sets are supported for the above compsets, including a nominal low-resosluton confiuguration and one regionally refined mesh. Additional regionally refined meshes and a high-resolution grid will become available in the near future.

| Grid Alias  |  Description  |
| ----------- |  ------------ |
|ne30pg2_r05_IcoswISC30E3r5 | For this grid set, the atmosphere is on the ne30pg2 cubed-sphere mesh with approximately 100km resolution, the land and river are on a 0.5deg x 0.5deg structured grid, and the ocean and sea ice are on a hexagonal mesh dervied from the dual of a 30km resolution icosahedral mesh with ice shelf cavities (wISC) around Antarctica.|
|northamericax4v1pg2_r025_IcoswISC30E3r5 | The atmosphere for this grid set uses North America regionally refined mesh from about 110 km down to 25 km over the refined region. The land and river are on 0.25deg x 0.25deg structured grid. The ocean and sea ice are on the same icosahedral mesh as for `ne30pg2_r05_IcoswISC30E3r5`.|

## Input data

Inputdata for coupled compsets at component model levels are the same as for the standalone component configurations
for a given forcing scenario (e.g., `1850` for the pre-industrial period,  `20TR` for the historical period, `2010`
for present-day condition, and `SSPs` for Shared Socioeconomic Pathways).
Between the coupled compsets, the differences are in the prescribed solar forcing, volcanic emissions,
atmospheric forcing data, and land use and land cover. The required inputdata for the pre-industrial and the historical periods
as well as the present-day condition are described in [the EAM User's Guide](https://e3sm-project.github.io/E3SM/EAM) and
[the ELM's User's Guide](https://e3sm-project.github.io/E3SM/ELM). Below are the prescribed forcing data for the SSP scenarios.

### [SSP245 forcing data](ssp245-forcings.md)

### [SSP370 forcing data](ssp370-forcings.md)

### [SSP585 forcing data](ssp585-forcings.md)
