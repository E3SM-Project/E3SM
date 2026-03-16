# CFMIP Observation Simulator Package (COSP) in EAMxx

[COSP](https://github.com/CFMIP/COSPv2.0)
is partially implemented and supported in EAMxx.
Currently, minimal outputs from the ISCCP, MODIS, and MISR simulators have been enabled.

## Running with COSP

COSP is implemented as a diagnostic in EAMxx. It is automatically triggered
when any of its output fields are requested in an output YAML file. There is
no need to add COSP to `atm_procs_list` — simply request COSP fields in your
output configuration.

COSP requires radiation fields (`dtau067`, `dtau105`, `cldfrac_rad`) and
microphysics fields (`eff_radius_qc`, `eff_radius_qi`) as inputs. These are
produced by the RRTMGP and P3 processes, so those must be in your
`atm_procs_list`.

COSP runs at the frequency specified by the output stream that requests its
fields (e.g., hourly, daily). There are no separate frequency controls.

COSP uses 10 subcolumns by default for internal subcolumn sampling using
`SCOPS`/`PREC_SCOPS`. The subcolumn count is not currently user-configurable
through the output YAML; it defaults to 10 for all grid resolutions.

A minimal output configuration example:

```yaml
averaging_type: average
fields:
  physics_pg2:
    field_names:
    - isccp_cldtot
    - isccp_ctptau
    - modis_ctptau
    - misr_cthtau
max_snapshots_per_file: 1
filename_prefix: eamxx
output_control:
  frequency: 1
  frequency_units: ndays
```

## Available output fields

The following output fields are available:

- `isccp_cldtot`
      - total cloud area from ISCCP simulator
- `isccp_ctptau`
      - ISCCP-simulated cloud top pressure/optical depth joint histogram
- `modis_ctptau`
      - MODIS-simulated cloud top pressure/optical depth joint histogram
- `misr_cthtau`
      - MISR-simulated cloud top height/optical depth joint histogram
