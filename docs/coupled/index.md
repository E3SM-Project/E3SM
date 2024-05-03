# E3SM Coupled Model Users's Guide

This User's Guide describes how to set up and run coupled E3SM

## Steps to build and run coupled E3SM

The basic steps for creating a case, building and running the E3SM are described in [CIME](https://github.com/ESMCI/cime)'s
[Case Control System](https://esmci.github.io/cime/versions/master/html/users_guide/index.html#case-control-system-part-1-basic-usage).
Running coupled E3SM or standalone component configurations follow the same procedure.
A coupled E3SM case can be created using a supported coupled compset along with a grid that together meets the simulation design.
In practice, the CIME-based procedure are often consolidated into a single run script with more readable settings for managing the simulation and
for better provenance. More details about the run script template and its settings are described in
[E3SM step-by-step guide](https://docs.e3sm.org/running-e3sm-guide/).

## Scientifically supported coupled compsets and grids

### Compsets

Coupled compsets consist of interactive prognostic Earth system components under specific forcing conditions.
Coupled compsets in E3SM are developed  for three science-driven simulation campaigns:

1) water cycle change and impacts, 2) human-earth system feedbacks, and 3) polar processes, sea-level rise and coastal impacts.

The standard coupled configurations -- which consist of active atmosphere, land, ocean and sea-ice components -- form the base physical
coupled system and are mainly designed for water cycle change and impacts simulation campaign.
The compsets for such configurations are listed below.

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

The E3SMv3 compsets for the other two science simulation campaigns are being finalized,
with additional components and/or features. The compset naming follows the same convention, e.g., `CRYO1850` and `CRYO1850-4xCO2` are with prognostic ice-shelf melt fluxes for polar processes simulation campaign.

### Grids

|Compset alias | Description |
|:-----------  |:----------- |

|`ne30pg2_r05_IcoswISC30E3r5` | E3SMv3 standard low-resolution tri-grid, with ne30pg2 (110 km) atmosphere, 0.5deg x 0.5deg land and river, and Icosahedral 30 km ocean/seaice mesh with ice shelves cavities (wISC), E3SMv3 (E3) revision r5 |
|`northamericax4v1pg2_r025_IcoswISC30E3r5` |  North America regionally refined (110 to 25 km) atmosphere, 0.25deg x 0.25deg land and river, and the same ocean/seaice mesh for the low-resolution tri-grid. |
|`ne120pg2_r025_IcoswISC30E3r5` |  High-resolution tri-grid, with ne120pg2 (25 km) atmosphere, 0.25deg x 0.25deg land and river. The ocean and seaice mesh to be updated with a variable 18 to 6 km mesh. |

## Input data

Inputdata for coupled compsets at component model levels are the same as for the standalone compnent configurations
for a given forcing scenario (e.g., `1850` for the pre-industrial period,  `20TR` for the historical period, `2010`
for present-day condition, and `SSPs` for Shared Socioeconomiuc Pathways of climate change scenarios).
Between the coupled compsets, the differences are in the prescribed solar forcing, volcanic emissions,
atmospheric forcing data, and land use and land cover. The required inputdata for the pre-industrial and the historical periods
as well as the present-day condition are described in [the EAM User's Guide](../../EAM/docs/user-guide/index.md) and
[the ELM's User's Guide](../../ELM/docs/user-guide/index.md). Below are the prescribed forcing data for SSP scenarios.

### [SSP245 forcing data](ssp245-forcings.md)

### [SSP370 forcing data](ssp370-forcings.md)

### [SSP585 forcing data](ssp585-forcings.md)
