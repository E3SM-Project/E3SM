# Generating Domain Files

Domain files are needed at runtime by the coupler, data models, and land model. The land model uses the mask to determine where to run and the coupler use the land fraction to merge fluxes from multiple surface types to the atmosphere above them.

Domain files are created from a conservative, monotone mapping file from the ocean grid (where the mask is defined) to the atmosphere grid.

## Environment

The new domain generation tool requires a few special packages, such as xarray, numba, and itertools. These are all included in the E3SM unified environment:
<https://e3sm.org/resources/tools/other-tools/e3sm-unified-environment/>

Alternatively, a simple conda environment can be created with the following command:

```shell
conda create --name example_env --channel conda-forge xarray numpy numba scikit-learn netcdf4
```

## Map File Generation

The map file used to generate the domain files can be created a few different ways. For a typical E3SM configuration we recommend using a conservative, monotone map. Here is an example command that can be used to generate one (as of NCO version 5.2.2)

```shell
ncremap -5 -a traave --src_grd=${OCN_GRID} --dst_grd=${ATM_GRID} --map_file=${MAP_FILE}
```

Note that existing ocean grid files can be found in the inputdata repository: `inputdata/ocn/mpas-o/<ocn_grid_name>/`

The atmosphere grid file should be on the "pg2" grid. These grid files are easily generated with three TempestRemap commands as follows:

```shell
NE=30
GenerateCSMesh --alt --res ${NE} --file ${GRID_FILE_PATH}/ne${NE}.g
GenerateVolumetricMesh --in ${GRID_FILE_PATH}/ne${NE}.g --out ${GRID_FILE_PATH}/ne${NE}pg2.g --np 2 --uniform
ConvertMeshToSCRIP --in ${GRID_FILE_PATH}/ne${NE}pg2.g --out ${GRID_FILE_PATH}/ne${NE}pg2_scrip.nc
```

For RRM grids the last two commands would be used on the exodus file produced by [SQuadGen](https://github.com/ClimateGlobalChange/squadgen) (See the [Adding Support for New Grids](https://docs.e3sm.org/user-guides/adding-grid-support-overview.md) tutorial for more information.).

## Running the Domain Generation Tool

Below is a typical example of how to invoke the domain generation tool from the command line:

```shell
NE=30
MAP_FILE=${MAP_FILE_ROOT}/map_oEC60to30v3_to_ne${NE}pg2_traave.20240313.nc
python generate_domain_files_E3SM.py -m ${MAP_FILE} -o oEC60to30v3 -l ne${NE}pg2 --date-stamp=9999 --output-root=${OUTPUT_ROOT}
```
