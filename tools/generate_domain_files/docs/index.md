# Generating Domain Files

Domain files are needed at runtime by the coupler, data models, and land model.
The land model uses the mask to determine where to run and the coupler use the
land fraction to merge fluxes from multiple surface types to the atmosphere
above them.

Domain files are created from a conservative, monotone mapping file from the 
ocean grid (where the mask is defined) to the atmosphere grid.

a mapping files created in the 
previous step, using a tool provided with CIME in 
${e3sm_root}/cime/tools/mapping/gen_domain_files. 


## Environement

The new domain generation tool requires a few special packages, 
such as xarray, numba, and itertools.
These are all included in the E3SM unified environment:
https://e3sm.org/resources/tools/other-tools/e3sm-unified-environment/

Alternatively, a simple conda environment can be created with the following command:
```
conda create --name example_env --channel conda-forge xarray numpy numba scikit-learn netcdf4
```

## Map File Generation

The map file used to generate the domain files can be created a few different ways.
For a typical E3SM configuration we recommend using a conservative, monotone map.
Here is an example command that can be used to generate one (as of NCO version 5.2.2)

```
ncremap -5 -a traave --src_grd=${OCN_GRID} --dst_grd=${ATM_GRID} --map_file=${MAP_FILE}
```

## Generating Domain Files

Below is a typical example of how to invoke the domain generation tool from the command line:

```
NE=30
MAP_FILE=${MAP_FILE_ROOT}/map_oEC60to30v3_to_ne${NE}pg2_traave.20240313.nc
python generate_domain_files_E3SM.py -m ${MAP_FILE} -o oEC60to30v3 -l ne${NE}pg2 --date-stamp=9999 --output-root=${OUTPUT_ROOT}
```
