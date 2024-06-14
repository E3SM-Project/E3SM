# Generate mapping files

[Back to Adding Support for New Grids](../adding-grid-support-step-by-step-guide.md)

In order to pass data between different components at runtime, a set of mapping files between each component is generated offline.

See [Recommended Mapping Procedures for E3SM Atmosphere Grids](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/178848194/Recommended+Mapping+Procedures+for+E3SM+Atmosphere+Grids) for a discussion of different remap algorithms and when to use each.

TempestRemap and ESMF are the backends that generate the mapping weights, but this is all nicely encapsulated using ncremap. Tempest is the preferred method for creating mapping files. ncremap will call TempestRemap or ESMF depending on the algorithm argument and input file types. If exodus files are provided (i.e. `*.g`) then TempestRemap commands will be used. The ESMF tools are adequate for making atmosphere-only-type component sets for E3SM, but this tool is less conservative than TempestRemap. If you are making grids for a coupled run, then TempestRemap should be used wherever possible. Currently, TempestRemap has trouble with masked grids, such those that are needed for land data generation, so ESMF is still required for certain tasks.

!!! NOTE
    This page is still under construction

## Activate the E3SM Unified Env

Perlmutter (NERSC):

```shell
source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
```

Chrysalis (LCRC):

```shell
source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_chrysalis.sh
```

For unsupported machines you may need to build the unified environment:

```shell
conda install -c conda-forge -c e3sm e3sm-unified
```

## Specify the Input Data Path

Perlmutter (NERSC):

```shell
```

Chrysalis (LCRC):

```shell
```

## Create Mapping Files

ncremap provides a convient 

We can use ncremap to generate ALL the needed mapping files between two grids, in this example the ne4 atmosphere and the oQU240 ocean grid (for the moment, we will put the land on the same grid as the atmosphere):

```shell
atm_grid_file=ne4.g
ocn_grid_file=/global/cfs/cdirs/e3sm/inputdata/cpl/gridmaps/oQU240/ocean.QU.240km.scrip.151209.nc
cd ${output_root} && ncremap -P mwf -s $ocn_grid_file -g $atm_grid_file --nm_src=oQU240 --nm_dst=ne4np4 --dt_sng=20181114
```
