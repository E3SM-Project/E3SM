# ARM diagnostics

## Overview

The ARM data-oriented metrics and diagnostics package (ARM Diags) was developed to facilitate the use of ARM data in earth system model evaluation and model intercomparison (Zhang et al., 2020)[@zhang_arm_2020]. It includes ARM data sets, compiled from multiple ARM data products, and a Python-based analysis toolkit for computation ad visualization. It also includes simulation data from models participating in CMIP, which allows modeling groups to compare a new, candidate version of their model to existing CMIP models. The ARM Diags has been applied in several model evaluation studies to help address a range of issues in earth system models (Zheng et al., 2023;[@zheng_assessment_2023] Emmenegger et al., 2022;[@emmenegger_evaluating_2022] Zhang et al., 2018[@zhang_causes_2018]). The Majority of ARM Diags sets are ported into E3SM Diags (Zhang et al., 2022)[@zhang_e3sm_2022] for routine evaluation of the model.

## To enable the use of ARM Diags

To enable using ARM Diags for a simulation, often, a new tape that output at high-frequency over limited-area (nearest grid box to supported ARM site) needs to be included in the namelist file, an example as follows:

```fortran
fincl7 = 'PS','Q','T','Z3','CLOUD','CONCLD','CLDICE',
   'CLDLIQ','FREQR','REI','REL','PRECT','TMQ','PRECC',
   'TREFHT','QREFHT','OMEGA','CLDTOT','LHFLX','SHFLX',
   'FLDS','FSDS','FLNS','FSNS','FLNSC','FSDSC','FSNSC',
   'AODVIS','AODABS','LS_FLXPRC','LS_FLXSNW',
   'LS_REFFRAIN','ZMFLXPRC','ZMFLXSNW','CCN1','CCN2',
   'CCN3','CCN4','CCN5','num_a1','num_a2','num_a3',
   'num_a4','so4_a1','so4_a2','so4_a3','AREL','TGCLDLWP',
   'AQRAIN','ANRAIN','FREQR','PRECL','RELHUM'
fincl7lonlat='262.5e_36.6n','203.4e_71.3n','147.4e_2.0s',
   '166.9e_0.5s','130.9e_12.4s','331.97e_39.09n'
```

Note that in this example fincl7 should set to write output at hourly (`nhtfrq = -1`). And here additional variables are included for ARM simulator analysis. The ARM site information is shown below:

```fortran
    "sgpc1": ["97.5W 36.4N Oklahoma ARM"],

    "nsac1": ["156.6W 71.3N Barrow ARM"],

    "twpc1": ["147.4E 2.1S Manus ARM"],

    "twpc2": ["166.9E 0.5S Nauru ARM"],

    "twpc3": ["130.9E 12.4S Darwin ARM"],

    "enac1": ["28.0E 39.1N Graciosa Island ARM"], 
```

## Diagnostics and metrics currently implemented in the ARM Diags

| Statistical Metrics       | Variables                                                       |  Time sampling     |
| ------------------------- | --------------------------------------------------------------- | -----------------  |
| Line plots and Taylor diagrams for annual cycle variability of each variable | Precipitation, column water vapor, surface energy budget components, near-surface temperature and specific humidity, surface pressure, total cloud fraction, and aerosol optical depth. | Monthly mean       |
| Contour and vertical profiles of annual cycle and diurnal cycle of cloud fraction | Vertical profiles of cloud fraction | Hourly mean       |
| Line and harmonic dial plots of diurnal cycle of precipitation | Surface precipitation rate | Hourly mean       |
| Probability density function (PDF) plots of precipitation rate | Surface precipitation rate | Hourly mean       |
| CCN Annual Cycles  | CCN number concentrations at 0.1%, 0.2%, 0.5% and 1.0% supersaturation levels | Hourly mean       |
| Aerosol Annual Cycles | Total aerosol number concentration | Hourly mean       |
| Aerosol Chemical Annual Cycles | Organic, sulfate, nitrate, ammonium, chloride mass concentration | Hourly mean       |

| Process-oriented metrics  | Variables                                                     |  Time sampling     |
| ------------------------- | ------------------------------------------------------------- | -----------------  |
| Convection Onset | 1. Surface precipitation rate <!-- markdownlint-disable MD033 --><br>  2. Column Precipitable Water Vapor | Hourly mean       |
| Aerosol-CCN Activation | 1. Total aerosol number concentration <br>  2. CCN number concentrations at different supersaturation levels (0.1%, 0.2%, 0.5% and 1.0) | Hourly mean       |
