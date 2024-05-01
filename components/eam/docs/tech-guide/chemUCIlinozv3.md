# chemUCI and Linoz

## Overview

Atmospheric interactive chemistry is handled by chemUCI (in the troposphere) and Linoz v3 (in the stratosphere). chemUCI consists of 28 advected tracers for the O3-CH4-HOx-NOx-NMVOCs chemistry. Compared to E3SMv2, the E3SMv3 linearized stratospheric chemistry scheme (Linoz v3) expends the interactive species to include O3, N2O, NOy, and CH4. The boundary between stratosphere and troposphere adopts the e90 tropopause algorithm.

## Namelist parameters

| Parameter                    | Description                                                              | Default value*         |
| ---------------------------- | ------------------------------------------------------------------------ | ---------------------- |
| `airpl_emis_file`            | Aviation emission                                                        |                        |
| `chlorine_loading_file`      | Chlorine loading                                                         |                        |
| `chlorine_loading_fixed_ymd` |                                                                          |                        |
| `chlorine_loading_type`      |                                                                          |                        |
| `ext_frc_specifier`          | 3-D emissions                                                            |                        |
| `ext_frc_cycle_yr`           |                                                                          |                        |
| `ext_frc_type`               |                                                                          |                        |
| `srf_emis_specifier`         | Surface emissions                                                        |                        |
| `srf_emis_cycle_yr`          |                                                                          |                        |
| `srf_emis_type`              | Upper bound of mean raindrop diameter                                    |                        |
| `linoz_data_file`            | Linoz data file                                                          |                        |
| `linoz_data_cycle_yr`        |                                                                          |                        |
| `linoz_data_path`            |                                                                          |                        |
| `linoz_data_type`            |                                                                          |                        |
| `lght_no_prd_factor`         | Lightning NOx emission factor                                            | `5.0`                  |
| `fstrat_efold_list`          | Tracer (from troposphere) list with e-folding decay in the stratosphere  |                        |

* Many of these namelist parameters specify input data files. Check the `atm_in` file for examples or refer to the [Users' Guide](../user-guide/index.md).
