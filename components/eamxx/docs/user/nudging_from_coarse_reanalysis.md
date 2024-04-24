# Nudging from coarse data

Because EAMxx is designed to support ultra-high resolutions (in fact, that was the initial reason for its inception), it is not feasible to produce nudging data at the same resolution.
Instead, in EAMxx, it is possible to nudge from coarse data.
This is done by remapping the coarse data provided by the user to the runtime physics grid of EAMxx.
In order to enable nudging from coarse data, the user must provide nudging data at the coarse resolution desired and an appropriate     ncremap-compatible mapping file.

## Example setup (current as of April 2024)

A user can produce coarse nudging data from running EAMxx or EAM at a ne30pg2 or any other applicable resolution.
Additionally, several users in the E3SM projects have produced nudging data at the ne30pg2 resolution from the MERRA2 and ERA5 datasets. Then, to enable nudging as a process, one must declare it in the `atm_procs_list` runtime parameter.

```shell
    ./atmchange physics::atm_procs_list="mac_aero_mic,rrtmgp,cosp,nudging"
```

The following options are needed to specify the nudging.

```shell
    ./atmchange nudging::nudging_filenames_patterns="${PSCRATCH}/mo_ne45pg2_nudg_trim_merra/*.nc"
    ./atmchange nudging::source_pressure_type=TIME_DEPENDENT_3D_PROFILE
    ./atmchange nudging::nudging_fields=U,V
    ./atmchange nudging::nudging_timescale=21600 #  6-hr
    ./atmchange nudging::nudging_refine_remap_mapfile="${PSCRATCH}/map_ne45pg2_to_ne256pg2.trbilin.20240410.nc"
```

To gain a deeper understanding of these parameters and options, please refer to code implementation of the nudging process.
