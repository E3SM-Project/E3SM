# Clean- and clean-clear-sky diagnostics

In order to decompose the aerosol effective radiative forcing, additional diagnostic radiation calls are needed.
These extra diagnostics are optionally added to the main radiation call. The extra diagnostics are:

- Clean-clear-sky fluxes: the fluxes that would be present if there were neither aerosols nor clouds, and are calculated by adding an additional radiation call at the very beginning of the logic before the optics class is endowed with aerosol and cloud properties.
- Clean-sky fluxes: the fluxes that would be present if there were no aerosols, and are calculated by adding an additional radiation call after substantiating an additional optics class, but not endowing it with aerosol properties.

## Example setup (current as of April 2024)

The extra calls are controlled by runtime flags `extra_clnclrsky_diag` and `extra_clnsky_diag` (they take either `true` or `false` as their values).

```shell
    ./atmchange extra_clnclrsky_diag=true
    ./atmchange extra_clnsky_diag=true
```

Below is an example output file to output the extra (clean and clean-clear-sky) radiation diagnostics atop the atmosphere.

```yaml
%YAML 1.1
---
filename_prefix: monthly.outputs
Averaging Type: Average
Max Snapshots Per File: 1
Fields:
  Physics PG2:
    Field Names:
    - SW_clnclrsky_flux_up_at_model_top
    - LW_clnclrsky_flux_up_at_model_top
    - SW_clnsky_flux_up_at_model_top
    - LW_clnsky_flux_up_at_model_top
output_control:
  Frequency: 1
  frequency_units: nmonths
  MPI Ranks in Filename: false

```
