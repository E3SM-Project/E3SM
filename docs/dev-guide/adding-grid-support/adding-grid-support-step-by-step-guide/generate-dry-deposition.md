# Generate a Dry Deposition File

Atmospheric dry deposition of aerosols at the surface depends on certain surface properties, such as soil type. In some cases these calculations can be handled in the land model and passed to the atmosphere through the coupler. However, with modal areosols this method is not adequate and we must recalculate these fields in the atmosphere (see subroutine `interp_map` in `components/eam/src/chemistry/mozart/mo_drydep.F90`).

For unstructured grids it was determined to create this offline interpolation tool rather than generalize the subroutine interp_map.

Be sure to activate the E3SM unified environment when performing the steps below.

## Map File Generation

The destination atmosphere grid file should be on the "pg2" grid. These grid files are easily generated with three TempestRemap commands as follows:

```shell
NE=30
GenerateCSMesh --alt --res ${NE} --file ${GRID_ROOT}/ne${NE}.g
GenerateVolumetricMesh --in ${GRID_ROOT}/ne${NE}.g --out ${GRID_ROOT}/ne${NE}pg2.g --np 2 --uniform
ConvertMeshToSCRIP --in ${GRID_ROOT}/ne${NE}pg2.g --out ${GRID_ROOT}/ne${NE}pg2_scrip.nc
```

The map file used to generate atmsrf files can be created a few different ways. For a typical E3SM configuration we recommend using a conservative, monotone map. Here is an example command that can be used to generate one (as of NCO version 5.2.2)

```shell
SRC_GRID=${DIN_LOC_ROOT}/../mapping/grids/1x1d.nc
DST_GRID=${GRID_ROOT}/ne${NE}pg2_scrip.nc
MAP_FILE=${MAP_ROOT}/map_1x1_to_ne${NE}pg2_traave.nc
ncremap -5 -a traave --src_grd=${SRC_GRID} --dst_grd=${DST_GRID} --map_file=${MAP_FILE}
```

For RRM grids the last two commands would be used on the exodus file produced by [SQuadGen](https://github.com/ClimateGlobalChange/squadgen) (See the [Adding Support for New Grids](https://docs.e3sm.org/user-guides/adding-grid-support-overview.md) tutorial for more information.).

## Generating a New Dry Desposition File (atmsrf)

```shell
VEGE_FILE=${DIN_LOC_ROOT}/atm/cam/chem/trop_mozart/dvel/regrid_vegetation.nc
SOIL_FILE=${DIN_LOC_ROOT}/atm/cam/chem/trop_mozart/dvel/clim_soilw.nc

python ~/E3SM/mkatmsrffile.py --map_file=${MAP_FILE} --vegetation_file=${VEGE_FILE} --soil_water_file=${SOIL_FILE} --output_root=${atmsrf_root} --dst_grid=ne${NE}pg2 --date-stamp=20240613
```

Back to step-by-step guide for [Adding Support for New Grids](../adding-grid-support-step-by-step-guide.md)
