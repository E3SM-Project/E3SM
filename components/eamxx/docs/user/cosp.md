# CFMIP Observation Simulator Package (COSP) in EAMxx

[COSP](https://github.com/CFMIP/COSPv2.0)
is partially implemented and supported in EAMxx.
Currently, minimal outputs from the ISCCP, MODIS, and MISR simulators have been enabled.

## Running with COSP

COSP is a *diagnostic*, not an atmosphere process. There is nothing to add to
`atm_procs_list`: asking for any of the fields it computes in an output stream is
what runs it. Those fields are

- `isccp_cldtot`
- `isccp_ctptau`
- `modis_ctptau`
- `misr_cthtau`

COSP computes all four in a single pass, so asking for one costs the same as
asking for all four, and asking for them across several output streams still runs
COSP only once per timestep.

!!! note "Changed in this version"

    COSP used to be an atmosphere process, turned on with
    `./atmchange physics::atm_procs_list="...,cosp"` and configured with
    `physics::cosp::cosp_frequency` / `cosp_frequency_units` / `cosp_subcolumns`.
    Those settings no longer exist. Output yamls that simply list the COSP field
    names keep working unchanged; only the `atm_procs_list` entry has to go.

    COSP now runs whenever an output stream evaluates it, rather than on its own
    `cosp_frequency`. For `INSTANT` output that is the stream's output frequency;
    for averaged output it is every timestep. Subcolumn sampling is fixed at 10
    subcolumns (`SCOPS`/`PREC_SCOPS`), the former default.

Output streams need to be added manually.
A minimal example:

```shell
./atmchange output_yaml_files=eamxx_daily_output.yaml
cat << EOF > eamxx_cosp_daily_output.yaml
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
EOF
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
