# Zhang and McFarlane deep convection scheme

## Overview

The ZM scheme (Zhang and McFarlane 1995) used in E3SMv3 is a bulk mass flux-type scheme; it has three components: a trigger for convection initiation, a cloud model including both updrafts and downdrafts, and a closure. The original CAPE-based trigger for convection was replaced by a trigger function based on dynamic CAPE generation by Xie et al. (2019) (see dCAPE-ULL description below for more details). The closure predicts cloud base mass flux using dilute CAPE (Neale et al. 2008). The updraft model is a bulk entraining plume model. Both updrafts and downdrafts are assumed saturated, with downdraft mass flux at the downdraft initiation level set proportional to the updraft cloud base mass flux. The microphysical processes inside the updrafts are represented by a convective microphysics scheme (see ZM convective microphysics description below). An additional adjustment is made to cloud base mass flux to incorporate the effect of large-scale circulation (see mass flux adjustment description below).

### dCAPE-ULL

A notable update related to clouds and precipitation in E3SMv2 is the use of a new convective trigger function described by Xie et al. (2019) in ZM to improve the simulation of precipitation and its diurnal cycle. The new convective trigger named as dCAPE-ULL uses the dynamic CAPE (dCAPE) trigger developed by Xie and Zhang (2000) with an unrestricted air parcel launch level (ULL) approach used by Wang et al. (2015). It was designed to address the unrealistically strong coupling of convection to the surface heating in ZM that often results in unrealistically too active model convection during the day in summer season over lands and improve the model capability to capture mid-level convection for nocturnal precipitation. 

### ZM convective microphysics

The convective microphysics scheme is based on the work of Song and Zhang (2011) to improve the representation of microphysical processes in convective clouds and their interaction with aerosol and stratiform clouds in GCMs. It explicitly treats the mass mixing ratio and number concentration of five hydrometeor species (cloud water, cloud ice, rain, snow, and graupel). The scheme is linked to stratiform cloud microphysics parameterization through convective detrainment of cloud liquid/ice water content and droplet/crystal number concentration, and to aerosols through cloud droplet activation and ice nucleation processes. Previous evaluations of the scheme showed that it improved the simulation of convective cloud properties and cloud hydrological cycle (Song et al. 2012, Storer et al. 2015). The assessment demonstrates that the convective microphysics scheme not only significantly improves the simulation of tropical variability across multiple scales but also enhances the simulation of climatological mean states.

### Mass flux adjustment

The convective mass flux adjustment (MAdj) is designed to represent the dynamical effects of large-scale vertical motion on convection. With MAdj, convection is enhanced (suppressed) when there is large-scale ascending (descending) motion at the planetary boundary layer top. The coupling of convection with the large-scale circulation significantly improves the simulation of climate variability across multiple scales from diurnal cycle, convectively coupled equatorial waves, to Madden-Julian oscillations (Song et al. 2023).

### MCSP

Due to inadequate model resolution, organized mesoscale convection cannot be resolved in general circulation models and thus needs to be parameterized. The Multiscale Coherent Structure Parameterization (MCSP) aims at representing the dynamical and physical effects of organized mesoscale convection.

MCSP applies a sinusoidal baroclinic profile in the temperature, moisture, and momentum fields to represent the impact. Moncrieff et al. (2017) and Chen et al. (2021) have found that by adding MCSP, the both the representation of large-scale precipitation systems and the modes of variability from Tropical waves are improved.

## Namelist parameters

| ZM Parameters             | Description                                                       | Default value          |
| ------------------------- | ----------------------------------------------------------------- | ---------------------- |
| `zmconv_ke`               | Tunable evaporation efficiency in ZM deep convection scheme       | `2.5E-6`               |
| `zmconv_tau`              | Relaxation time in ZM deep convection scheme                      | `3600`                 |
| `zmconv_dmpdz`            | Parcel fractional mass entrainment rate                           | `-0.7E-3`              |
| `zmconv_alfa`             | Initial downdraft mass flux fraction                              | `0.14D0`               |
| `zmconv_tiedke_add`       | Temperature perturbation of an air parcel                         | `0.8D0`                |
| `zmconv_cape_cin`         | Number of negative buoyancy regions that are allowed              | `1`                    |

| dCAPE-ULL Parameters      | Description                                                       | Default value          |
| ------------------------- | ----------------------------------------------------------------- | ---------------------- |
| `zmconv_trigdcape_ull`    | DCAPE trigger along with unrestricted launching level for ZM deep convection scheme        | `.true.`  |
| `zmconv_trig_dcape_only`  | DCAPE only trigger for ZM deep convection scheme                  | `.false.` <br> If true, zmconv_trigdcape_ull must be false to use the dcape only trigger. |
| `zmconv_trig_ull_only`    | Use unrestricted launching level (ULL) only trigger for ZM deep convection scheme          | `.false.` <br> If true, zmconv_trigdcape_ull must be false to use the ull only trigger. |

| Conv. micro. Parameters   | Description                                                       | Default value          |
| ------------------------- | ----------------------------------------------------------------- | ---------------------- |
| `zmconv_microp`           | Convective microphysics option in ZM convection scheme            | `true`                 |
| `zmconv_auto_fac`         | Cloud droplet-rain autoconversion enhancement factor in the convective microphysics scheme | `7.0`     |
| `zmconv_accr_fac`         | Cloud droplet-rain accretion enhancement factor in the convective microphysics scheme      | `1.5`     |
| `zmconv_micro_dcs`        | Autoconversion size threshold for cloud ice to snow (m)           | `150.E-6`                 |

| Mass flux adj. Parameters | Description                                                       | Default value          |
| ------------------------- | ----------------------------------------------------------------- | ---------------------- |
| `zmconv_clos_dyn_adj`     | Apply mass flux adjustment to ZM convection scheme               | `true`                 |

| MCSP Parameters              | Description                                                       | Default value          |
| ---------------------------- | ----------------------------------------------------------------- | ---------------------- |
| `zmconv_mcsp_heat_coeff`     | MCSP heating coefficient                                          | `0.3`                  |
| `zmconv_mcsp_moisture_coeff` | MCSP moisture coefficient                                         | `0.0`                  |
| `zmconv_mcsp_uwind_coeff`    | MCSP zonal wind coefficient                                       | `0.0`                  |
| `zmconv_mcsp_vwind_coeff`    | MCSP meridional wind coefficient                                  | `0.0`                  |