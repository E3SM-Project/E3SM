# Nudging in EAMxx

Nudging is supported in EAMxx.
Currently, it is possible to nudge EAMxx to the output from a different EAMxx run or to reanalysis. Nudging data can be on your model grid or an arbitrary coarser grid. Inline interpolating of finer-grid nudging data to a coarser model resolution isn't implemented yet but may be in the future.

## Data to nudge towards

The user is expected to prepapre (and then use `atmchange` to point to) nudging data files that are compliant with EAMxx specification.
In practice, this means that the data files must contain variable names known to EAMxx (only U, V, T_mid, and qv are supported now).
The files can be specified with an explicit list or via pattern matching.
The files must contain an arbitary global attribute `case_t0`, and it is recommended to be the same as the time dimension unit (the files must be time-dimensioned).
Finally, the dimension order must be such that `lev` is the last dimension, so most likely, the user must transpose the dimensions.

## Pressure in the nudging data

Pressure can be explicitly provided in the nudging data as time-varying `p_mid` corresponding to the option `TIME_DEPENDENT_3D_PROFILE` for `source_pressure_type`.
Alternatively, the data can contain a time-invariant pressure variable `p_lev` corresponding to the option `TIME_DEPENDENT_3D_PROFILE` for `source_pressure_type`.

## Weighted nudging for RRM applications

In regionally refined model applications, it is possible to use weighted nudging, for example, to avoid nudging the refined region.
To achieve that, the user can use `atmchange` to set `use_nudging_weights` (boolean) and provide `nudging_weights_file` that has the weight to apply for nudging (for example, zeros in the refined region).
Currently, weighted nudging is only supported if the user provides the nudging data at the target grid.

## Example setup (current as of April 2024)

To enable nudging as a process, one must declare it in the `atm_procs_list` runtime parameter.

```shell
./atmchange physics::atm_procs_list="mac_aero_mic,rrtmgp,cosp,nudging"
```

The following options are needed to specify the nudging.

```shell
./atmchange nudging::nudging_filenames_patterns="/pathto/nudging_files/*.nc" # can provide file name explicitly here instead (or multiple patterns)
./atmchange nudging::source_pressure_type=TIME_DEPENDENT_3D_PROFILE # see section on pressure
./atmchange nudging::nudging_fields=U,V # can include qv, T_mid as well
./atmchange nudging::nudging_timescale=21600 # in seconds
```

To gain a deeper understanding of these parameters and options, please refer to code implementation of the nudging process.
