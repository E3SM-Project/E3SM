# CFMIP Observation Simulator Package (COSP) in EAMxx

[COSP](https://github.com/CFMIP/COSPv2.0)
is partially implemented and supported in EAMxx.
Currently, minimal outputs from the ISCCP, MODIS, and MISR simulators have been enabled.

## Running with COSP

Turning COSP on simply requires adding the `cosp` process to `atm_procs_list`
via `atmchange` in a case directory:

```shell
./atmchange physics::atm_procs_list="mac_aero_mic,rrtmgp,cosp"
```

Additionally, the frequency at which COSP is run can be configured via `atmchange`:

```shell
./atmchange physics::cosp::cosp_frequency_units="steps"
./atmchange physics::cosp::cosp_frequency=1
```

COSP can be run with or without subcolumn sampling.
This is configured by changing the `cosp_subcolumns` namelist variable via `atmchange`.
A value of 1 implies *no* subcolumn sampling, while values greater than 1
specify the number
of subcolumns to use for subcolumn sampling (assuming maximum-random overlap).
E.g.,

```shell
./atmchange physics::cosp:cosp_subcolumns=1
```

would disable subcolumn sampling, while

```shell
./atmchange physics::cosp::cosp_subcolumns=10
```

would use 10 subcolumns for the COSP internal subcolumn sampling using `SCOPS`/`PREC_SCOPS`.
The default for high resolution cases (e.g., `ne1024`) should be to *not* use
subcolumns, while lower resolutions (e.g., `ne30`) should enable subcolumn sampling.

Output streams need to be added manually.
A minimal example:

```shell
./atmchange output_yaml_files=eamxx_daily_output.yaml
cat << EOF > eamxx_cosp_daily_output.yaml
Averaging Type: Average
Fields:
  Physics PG2:
    Field Names:
    - isccp_cldtot
    - isccp_ctptau
    - modis_ctptau
    - misr_cthtau
    - cosp_sunlit
Max Snapshots Per File: 1
filename_prefix: eamxx
output_control:
  Frequency: 1
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
- `cosp_sunlit`
      - sunlit flag aggregated at COSP frequency for renormalizing daytime averages

ISCCP, MODIS, and MISR outputs are valid only for daytime/sunlit columns
(to be consistent with available satellite retrievals).
In order to aggregate only daytime columns in time averages, these outputs are
multiplied by the sunlit flag (0 or 1) at each COSP calculation time.
Time averages of these quantities are then aggregated, along with the COSP
sunlit flag each time COSP is called.
In order to back out the daytime-only time averages from the outputs,
one needs to divide the output fields by `cosp_sunlit`.
E.g.,

```shell
isccp_ctptau = mean(isccp_ctptau) / mean(cosp_sunlit)
```
