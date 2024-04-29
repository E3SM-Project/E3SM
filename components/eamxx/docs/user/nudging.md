# Nudging in EAMxx

Nudging is supported in EAMxx.
Currently, it is possible to nudge EAMxx to the output from a different EAMxx run or to reanalysis. Nudging data can be on your model grid or an arbitrary coarser grid. Inline interpolating of finer-grid nudging data to a coarser model resolution isn't implemented yet but may be in the future. 

## Example setup (current as of April 2024)

To enable nudging as a process, one must declare it in the `atm_procs_list` runtime parameter.

```shell
    ./atmchange physics::atm_procs_list="mac_aero_mic,rrtmgp,cosp,nudging"
```

The following options are needed to specify the nudging.

```shell
./atmchange nudging::nudging_filenames_patterns="/pathto/nudging_files/*.nc"
./atmchange nudging::source_pressure_type=TIME_DEPENDENT_3D_PROFILE
./atmchange nudging::nudging_fields=U,V
./atmchange nudging::nudging_timescale=21600
```

To gain a deeper understanding of these parameters and options, please refer to code implementation of the nudging process.
In particular, depending on the specific needs, one may utilize weighted nudging (for regionally refined model runs) or nudging from coarse data via specifying a map in the `atmchange` calls above.
