# Cloud Layers Unified By Binormals

## Overview

The Cloud Layers Unified By Binormals (CLUBB) parameterization is a parameterization of subgrid-scale turbulence and clouds.  It prognoses turbulent fluxes of heat, moisture, and momentum, and it diagnoses the liquid cloud fraction and liquid water mixing ratio.  To do so, it prognoses higher-order turbulence moments and closes those prognostic equations by use of an assumed double-Gaussian shape of the subgrid probability density function.  CLUBB operates throughout the troposphere, but it contributes especially to the planetary boundary layer and low-cloud regimes, including stratocumulus and shallow cumulus regimes.  

### References

- Larson (2022) [@larson_clubb-silhs_2022]
- Bogenschutz et al. (2018) [@bogenschutz_path_2018]
- Larson and Golaz (2005) [@larson_using_2005]
- Golaz et al. (2002) [@golaz_pdf-based_2002]

## Namelist parameters

| Parameter      | Description                                                                                 | Default value  |
| -------------- | ------------------------------------------------------------------------------------------- | -------------- |
| `gamma_coef`   | Width of vertical velocity within a Gaussian PDF component at low skewness                  | `0.12`         |
| `gamma_coefb`  | Width of vertical velocity within a Gaussian PDF component at high skewness                 | `0.28`         |
| `C8`           | Coefficient of damping of third moment of vertical velocity, w’3                            | `5.2`          |
| `C1`           | Coefficient of damping of second vertical moment of vertical velocity, w’2, at low skewness | `2.4`          |
| `C14`          | Coefficient of damping of second horizontal moments of vertical velocity, u’2 and v’2       | `2.0`          |
| `c_k10`        | Ratio of diffusivity of momentum to heat                                                    | `0.35`         |
