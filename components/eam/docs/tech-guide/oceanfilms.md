# OCEANFILMS

## Overview

E3SM (v1-v3) uses the OCEANFILMS (Organic Compounds from Ecosystems to Aerosols: Natural Films and Interfaces via Langmuir Molecular Surfactants) parameterization to represent sea spray organic aerosol emissions.  OCEANFILMS is a physically based model that links sea spray chemistry with ocean biogeochemistry using a Langmuir partitioning approach.  The underlying physical assumptions and parameterization are described in Burrows et al. (2014); the implementation in E3SM and impact on clouds and climate are documented in Burrows et al. (2022).

## Namelist parameters

| Parameter                 | Description                                                       | Default value          |
| ------------------------- | ----------------------------------------------------------------- | ---------------------- |
| `mam_mom_cycle_yr`       |                                                                    | `1`                    |
| `mam_mom_datapath` | Full pathname of the directory that contains the files specified in mam_mom_filelist  | `'atm/cam/chem/trop_mam/marine_BGC/'`                 |
| `mam_mom_filename`     | Filename of file that contains a sequence of filenames for prescribed marine organic matter ocean concentrations.  The filenames in this file are relative to the directory specified by mam_mom_datapath.| `'monthly_macromolecules_0.1deg_bilinear_latlon_year01_merge_date.nc'` |
| `mam_mom_rmfile`   | Remove the file containing prescribed aerosol deposition fluxes from local disk when no longer needed. | `FALSE`                |
| `mam_mom_specifier`     | Names of variables containing aerosol data in the prescribed aerosol datasets. | `'chla:CHL1','mpoly:TRUEPOLYC','mprot:TRUEPROTC','mlip:TRUELIPC'`                 |
| `mam_mom_datatype`       | Type of time interpolation for data in mam_mom files. Can be set to `'CYCLICAL'`, `'SERIAL'`, `'INTERP_MISSING_MONTHS'`, or `'FIXED'`. | `'CYCLICAL'`                |
| `mam_mom_cycle_yr`         | The  cycle year of the prescribed aerosol flux data if mam_mom_type is `'CYCLICAL'`. Format: YYYY   | `1`               |
| `mam_mom_fixed_ymd`        | The date at which the prescribed aerosol flux data is fixed if mam_mom_type is `'FIXED'`. Format: YYYYMMDD | `0`                |
| `mam_mom_fixed_tod`  | The time of day (seconds) corresponding to mam_mom_fixed_ymd at which the prescribed aerosol flux data is fixed if mam_mom_type is 'FIXED'. | `0`           |
| `mam_mom_bubble_thickness`   | Bubble film thickness (in m) for marine organic aerosol emission mechanism.  The physically reasonable range is approximately (0.1 - 1) x 10^ -6. | `0.1e-6`            |
| `mam_mom_mixing_state`              | Switch to select mixing state assumption in marine organic aerosol code. Currently implemented options: 0 : total external mixture, add to mass; 1 : total external mixture, replace mass; 2 : total internal mixture, add to mass; 3 : total internal mixture, replace mass. | `0` [Note: set to 3 in the atm_in namelist]        |
| `mam_mom_parameterization`     | Selection of alternate parameterizations for marine organic matter emissions.  Set fmoa=1 for Burrows et al., 2014 parameterization; fmoa=2 for Gantt et al. (2011, ACP) parameterization; fmoa=3 for simple parameterization based on Quinn et al., 2014; fmoa=4 for
Rinaldi et al. (JGR, 2013).* | `1`                 |
*Note: non-default values have not been carefully tested and may not work as expected.