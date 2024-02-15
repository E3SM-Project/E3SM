# Nudging from coarse data

Because EAMxx is designed to support ultra-high resolutions (in fact, that was the initial reason for its inception), it is not feasible to produce nudging data at the same resolution.
Instead, in EAMxx, it is possible to nudge from coarse data.
This is done by remapping the coarse data provided by the user to the runtime physics grid of EAMxx.
In order to enable nudging from coarse data, the user must provide nudging data at the coarse resolution desired and an appropriate     ncremap-compatible mapping file.

## Example setup

A user can produce coarse nudging data from running EAMxx or EAM at a ne30pg2 or any other applicable resolution.
Additionally, several users in the E3SM projects have produced nudging data at the ne30pg2 resolution from the MERRA2 and ERA5 datasets.
A limitation for now is that the nudging data must be provided explicitly, either as one file or as a list of files.
This can be problematic for long list of files, but we are working on a solution to this problem.

Let's say that the nudging data is provided as one file in the following path: `/path/to/nudging_data_ne4pg2_L72.nc`.
Then, a mapping file is provided as `/another/path/to/mapping_file_ne4pg2_to_ne120pg2.nc`.
Then if the physics grid is ne120pg2, the user must enable the nudging process, specify the nudging files, and provide the specifies the nudging data and a remap file.
In other words, the following options are needed:

```shell
./atmchange atm_procs_list=(sc_import,nudging,homme,physics,sc_export)
./atmchange nudging_fields=U,V
./atmchange nudging_filename=/path/to/nudging_data_ne4pg2_L72.nc
./atmchange nudging_refine_remap_mapfile=/another/path/to/mapping_file_ne4pg2_to_ne120pg2.nc
```
