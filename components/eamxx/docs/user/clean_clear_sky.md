# Clean- and clean-clear-sky diagnostics

In order to decompose the aerosol effective radiative forcing, additional diagnostic radiation calls are needed.
These extra diagnostics are optionally added to the main radiation call. The extra diagnostics are:

- Clean-clear-sky fluxes: the fluxes that would be present if there were neither aerosols nor clouds, and are calculated by adding an additional radiation call at the very beginning of the logic before the optics class is endowed with aerosol and cloud properties.
- Clean-sky fluxes: the fluxes that would be present if there were no aerosols, and are calculated by adding an additional radiation call after substantiating an additional optics class, but not endowing it with aerosol properties.

It was necessary to add an additional optics class because the original optics class is endowed with aerosols before clouds (in order to calculate the clear-sky fluxes).

## Example setup (current as of April 2024)

The extra calls are controlled by runtime flags `extra_clnclrsky_diag` and `extra_clnsky_diag` (they take either `true` or `false` as their values).

```shell
    ./atmchange extra_clnclrsky_diag=true
    ./atmchange extra_clnsky_diag=true
```

An example output file with the additional radiation diagnostics requested.

```yaml
%YAML 1.1
---
filename_prefix: monthly.outputs
Averaging Type: Average
Max Snapshots Per File: 1
Fields:
  Physics PG2:
    Field Names:
    #2D vars
    - SW_flux_up_at_model_top
    - SW_flux_dn_at_model_top
    - LW_flux_up_at_model_top
    - SW_clnclrsky_flux_up_at_model_top
    - LW_clnclrsky_flux_up_at_model_top
    - SW_clrsky_flux_up_at_model_top
    - LW_clrsky_flux_up_at_model_top
    - SW_clnsky_flux_up_at_model_top
    - LW_clnsky_flux_up_at_model_top
    - SW_flux_dn_at_model_bot
    - SW_clnclrsky_flux_dn_at_model_bot
    - SW_clrsky_flux_dn_at_model_bot
    - SW_clnsky_flux_dn_at_model_bot
    - SW_flux_up_at_model_bot
    - SW_clnclrsky_flux_up_at_model_bot
    - SW_clrsky_flux_up_at_model_bot
    - SW_clnsky_flux_up_at_model_bot
    - LW_flux_dn_at_model_bot
    - LW_clnclrsky_flux_dn_at_model_bot
    - LW_clrsky_flux_dn_at_model_bot
    - LW_clnsky_flux_dn_at_model_bot
    - LW_flux_up_at_model_bot
    - LongwaveCloudForcing
    - ShortwaveCloudForcing
    - ps
    - SeaLevelPressure
    - T_2m
    - qv_2m
    - surf_radiative_T
    - VapWaterPath
    - IceWaterPath
    - LiqWaterPath
    - RainWaterPath
    - ZonalVapFlux
    - MeridionalVapFlux
    - surf_evap
    - surf_sens_flux
    - surface_upward_latent_heat_flux
    - precip_liq_surf_mass_flux
    - precip_ice_surf_mass_flux
    - landfrac
    - ocnfrac
    - PotentialTemperature_at_700hPa
    - PotentialTemperature_at_1000hPa
    - omega_at_500hPa
    - RelativeHumidity_at_700hPa
output_control:
  Frequency: 1
  frequency_units: nmonths
  MPI Ranks in Filename: false

```
