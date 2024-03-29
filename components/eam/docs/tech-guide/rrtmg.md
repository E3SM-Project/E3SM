# Rapid Radiative Transfer Model for GCMs

## Overview

The calculation of radiative energy flux through the atmosphere is done using the RRTMG radiation package (Iacono et al., 2008; Mlawer et al., 1997). The details are consistent with the implementation in CAM5 described in Neale et al. 2010. Radiative fluxes are broadly split into shortwave and longwave and computed by separate codes. The shortwave solver uses the 2-stream approximation, while the longwave is an absorption/emission code. Both shortwave and longwave use the correlated k-distribution method for integration of fluxes. Subgrid cloud overlap is accounted for using the Monte Carlo Independent Column Approximation (MCICA; Pincus and Morcrette, 2003), assuming the cloudy portions of the column are maximally overlapped in vertically contiguous layers and randomly overlapped when two layers are separated by a completely clear layer. Cloud optics are parameterized as described in Neale et al.(2010).

## Namelist parameters

| Parameter                 | Description                                                       | Default value          |
| ------------------------- | ----------------------------------------------------------------- | ---------------------- |
| `iradsw`                  | Frequency for updating shortwave fluxes and heating rate; iradsw > 0 interpreted as number of timesteps, iradsw < 0 interpreted as hours; iradsw = 0 disables shortwave radiation entirely | `-1`                 |
| `iradlw`                  | Frequency for updating longwave fluxes and heating rate; iradlw > 0 interpreted as number of timesteps, iradlw < 0 interpreted as hours; iradlw = 0 disables longwave radiation entirely   | `-1`                 |
| `irad_always`             | Length of time in timesteps (irad_always > 0) or in hours (irad_always < 0) SW/LW radiation will be run continuously from the start of an initial or restart run                           | `0`                  |
| `use_rad_dt_cosz`         | If true, use the radiation dt for all cosz calculations; calculates solar zenith angle averaged over a time step. In default model solar zenith angle is held constant over time           | `.true.` <br> (set by namelist_defaults_eam.xml for default physics)   |
| `spectralflux`            | Calculate fluxes (up and down) per band                                                                                                                                                    | `.false.`            |
| `liqcldoptics`            | Choice of cloud optical property parameterization for liquid clouds. Valid options are ‘slingo’ or ‘gammadist’                                                                             | `gammadist`          |
| `icecldoptics`            | Choice of cloud optical property parameterization for ice clouds. Valid options are ‘ebertcurry’ or ‘mitchell’                                                                             | `mitchell`           |