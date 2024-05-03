# Cloud Feedback Model Intercomparison Project (CFMIP) Observation Simulator Package

## Overview

The Cloud Feedback Model Intercomparison Project (CFMIP) Observation Simulator Package (COSP; Bodas-Salcedo et al., 2011;[@bodas-salcedo_cosp_2011] Swales et al., 2018;[@swales_cloud_2018] Zhang et al., 2019;[@zhang_evaluation_2019] Zhang et al., 2024[@zhang_understanding_2024]) was developed to improve the consistency between model clouds and satellite observations. COSP contains several independent satellite simulators for better comparing model clouds with satellite measurements collected by the International Satellite Cloud Climatology Project (ISCCP), the Moderate Resolution Imaging Spectroradiometer (MODIS), the Multi-angle Imaging SpectroRadiometer (MISR), Cloud-Aerosol Lidar and Infrared Pathfinder Satellite Observation (CALIPSO), and CloudSat. The use of satellite simulators will not only make a fairer comparison between model clouds and satellite data but also allow a more in-depth analysis of clouds. For example, clouds can be assessed in terms of their optical properties and vertical location, which dictate their radiative effects.

## Namelist parameters

[COSP Namelist Parameters](../user-guide/namelist_parameters.md#cloud-feedback-model-intercomparison-project-cfmip-observation-simulator-package)

## To turn on COSP outputs

Run (add to the run script) the following command before running `./case.setup`

`./xmlchange --id CAM_CONFIG_OPTS --append --val='-cosp'`

| Related Outputs           | Description                                                       |
| ------------------------- | ----------------------------------------------------------------- |
| `FISCCP1_COSP`            | Grid-box fraction covered by each ISCCP D level cloud type        |
| `CLMODIS`                 | MODIS Cloud Area Fraction                                         |
| `CLD_MISR`                | Cloud Fraction from MISR Simulator                                |
| `CLDTOT_CAL`              | Calipso Total Cloud Fraction                                      |
| `CLDHGH_CAL`              | Calipso High-level Cloud Fraction                                 |
| `CLDMED_CAL`              | Calipso Mid-level Cloud Fraction                                  |
| `CLDLOW_CAL`              | Calipso Low-level Cloud Fraction                                  |
| `CLD_CAL_TMPLIQ`          | Calipso Liquid Cloud Fraction as a function of temperature        |
| `CLD_CAL_TMPICE`          | Calipso Ice Cloud Fraction as a function of temperature           |
