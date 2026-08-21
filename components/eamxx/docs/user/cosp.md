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

## Controlling how often COSP runs

COSP is expensive, and an averaged output stream would otherwise evaluate it
every timestep. It therefore throttles itself, with the same two knobs it had as
a process, given in the `diag_params` section of the output yaml:

```yaml
diag_params:
  Cosp:
    cosp_frequency: 1
    cosp_frequency_units: hours   # 'hours' or 'steps'
    cosp_subcolumns: 10
```

Between runs the COSP fields keep the values of the last run, exactly as they
did when COSP was a process. The defaults are `cosp_frequency: 1` and
`cosp_frequency_units: steps`, i.e. every timestep, so setting them is worth it
for anything but a short run.

`diag_params` is keyed by the *diagnostic* name (`Cosp`), not by the name of any
field it computes, since one COSP pass produces all four. For the same reason,
one COSP diagnostic is shared by every output stream that asks for its fields:
if more than one stream configures it, they must all configure it identically,
and it is an error if they do not.

!!! note "Changed in this version"

    COSP used to be an atmosphere process, turned on with
    `./atmchange physics::atm_procs_list="...,cosp"` and configured with
    `physics::cosp::cosp_frequency` / `cosp_frequency_units` / `cosp_subcolumns`.
    Output yamls that simply list the COSP field names keep working unchanged,
    and the `atm_procs_list` entry has to go. The three settings move to the
    `diag_params: Cosp:` section of the output yaml, keeping their names,
    meanings, and defaults; `cosp_frequency: 0` ("never run COSP") is the one
    thing that no longer means anything, and is now an error.

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
