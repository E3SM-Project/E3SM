# Generate a Topography File

Topography needs to be interpolated from a high resolution USGS file, and then doctored up a bit to allow the model to run stably with the new topography. The tool chain used to compute topography is documented in the following paper:  

[P.H. Lauritzen, J.T. Bacmeister, P.F. Callaghan, M. Taylor,  NCAR_Topo (v1.0): NCAR global model topography generation software for unstructured grids, Geosci. Model Dev., 8, 3975-3986, 2015.](https://www.geosci-model-dev.net/8/3975/2015/)

## Input Topography Data

Traditionally the topography generation for E3SM start with **USGS-topo-cube3000.nc**, which is a high-resolution topography dataset on a 3km cubed sphere grid derived from 1 km resolution source data.  This file is located in the [CESM inputdata server here](https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/atm/cam/hrtopo/).

For target resolutions of 3 km or finer it is recommended to use **USGS-topo-cube12000.nc**, which was created by Jishi Zhang in 2024. This file has a resolution of 750m created from a  500m/250m USGS GMTED2010 source DEM dataset (see [here](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/4189520033/800m+cubed+topo+generation+from+GMTED2010+15s+DEM) for more information).

## homme_tool

homme_tool is included in the repository (`components/homme/test/tool`). This builds parts of the HOMME code directly, which is preferred to creating an offline tool. This tool allows us to process the topography data and address the following requirements:

- The dycore needs surface geopotential (`phi_s`) at GLL nodes (np4)
- The physics needs various things at FV cell centers (pg2)
    - surface geopotential (`phi_s`)
    - sub-grid variance of topography (`SGH` and `SGH30`)
    - Land area fraction (`LANDFRAC`)
- The np4 to pg2 map of `phi_s` data must be equal to the FV `phi_s` data
- HOMME's smoothing operator must be used on `phi_s` for stability

## Step-by-Step Topography Generation

1. Set some helpful environmental variables used in commands below
    ```shell
    e3sm_root=<path to E3SM source>
    grid_root=<path to write grid files>
    map_root=<path to write map files>
    topo_root=<path to write new topo files>

    NE=<ne value for target grid>
    ```
    ```shell
    e3sm_root=/lustre/orion/cli115/proj-shared/hannah6/e3sm_scratch/tmp_clone
    grid_root=/lustre/orion/cli115/proj-shared/hannah6/HICCUP/files_grid
    map_root=/lustre/orion/cli115/proj-shared/hannah6/HICCUP/files_map
    topo_root=/lustre/orion/cli115/proj-shared/hannah6/HICCUP/files_topo

    NE=30

    e3sm_root=/pscratch/sd/w/whannah/e3sm_scratch/tmp_clone
    grid_root=/pscratch/sd/w/whannah/e3sm_scratch/files_grid
    map_root=/pscratch/sd/w/whannah/e3sm_scratch/files_map
    topo_root=/pscratch/sd/w/whannah/e3sm_scratch/files_topo

    e3sm_root=/lcrc/group/e3sm/ac.whannah/scratch/chrys/tmp_clone
    grid_root=/pscratch/sd/w/whannah/e3sm_scratch/files_grid
    map_root=/pscratch/sd/w/whannah/e3sm_scratch/files_map
    topo_root=/pscratch/sd/w/whannah/e3sm_scratch/files_topo
    ```

1. Set machine specific environment (this sets `DIN_LOC_ROOT`)
    ```shell
    cd ${e3sm_root}/components/homme
    eval $(${e3sm_root}/cime/CIME/Tools/get_case_env)
    ${e3sm_root}/cime/CIME/scripts/configure --macros-format Makefile --mpilib mpi-serial
    source .env_mach_specific.sh
    # above approach is too slow - will this work?
    eval $(${e3sm_root}/cime/CIME/Tools/get_case_env)
    ${e3sm_root}/cime/CIME/scripts/configure && source .env_mach_specific.sh
    ```

1. Create grid files for the input high res topo
    ```shell
    GenerateCSMesh --alt --res 3000  --file ${grid_root}/exodus_ne3000.g
    ConvertExodusToSCRIP --in ${grid_root}/exodus_ne3000.g  --out ${grid_root}/scrip_ne3000pg1.nc
    ```

1. Create grid files target EAM grid
    ```shell
    GenerateCSMesh --alt --res ${NE} --file ${grid_root}/exodus_ne${NE}.g
    GenerateVolumetricMesh --in ${grid_root}/exodus_ne${NE}.g --out ${grid_root}/exodus_ne${NE}pg2.g --np 2 --uniform
    ConvertMeshToSCRIP --in ${grid_root}/exodus_ne${NE}pg2.g --out ${grid_root}/scrip_ne${NE}pg2.nc
    ```

1. Create map file
    ```shell
    ncremap -a fv2se_flx -5 --src_grd=${grid_root}/scrip_ne3000pg1.nc  --dst_grd=${grid_root}/exodus_ne${NE}.g --map_file=${map_root}/map_ne3000pg1_to_ne${NE}np4.nc
    ```

1. Map high-res topo to target np4 grid
    ```shell
    map_file=${map_root}/map_ne3000pg1_to_ne${NE}np4.nc
    ncremap -m ${map_file} -i ${DIN_LOC_ROOT}/atm/cam/hrtopo/USGS-topo-cube3000.nc -o ${topo_root}/USGS-topo_ne${NE}np4.nc
    ```

1. Build homme_tool
    ```shell
    eval $(${e3sm_root}/cime/CIME/Tools/get_case_env)
    ${e3sm_root}/cime/CIME/scripts/configure && source .env_mach_specific.sh

    cd ${e3sm_root}/components/homme

    E3SM_MACH=$(${e3sm_root}/cime/CIME/Tools/list_e3sm_tests cime_tiny | grep ERS | cut -f 4 -d . | cut -f 1 -d _ )
    E3SM_COMP=$(${e3sm_root}/cime/CIME/Tools/list_e3sm_tests cime_tiny | grep ERS | cut -f 4 -d . | cut -f 2 -d _ )

    # mach_file=${e3sm_root}/components/homme/cmake/machineFiles/perlmutter-nocuda-gnu.cmake
    # mach_file=${e3sm_root}/components/homme/cmake/machineFiles/perlmutter-gnu.cmake
    # mach_file=${e3sm_root}/components/homme/cmake/machineFiles/crusher-gpumpi.cmake
    # mach_file=${e3sm_root}/components/homme/cmake/machineFiles/chrysalis.cmake
    mach_file=${e3sm_root}/cime_config/machines/cmake_macros/gnu_chrysalis.cmake
    # E3SM_MACH=$(${e3sm_root}/cime/CIME/Tools/list_e3sm_tests cime_tiny | grep ERS | cut -f 4 -d . | cut -f 1 -d _ )
    # E3SM_COMP=$(${e3sm_root}/cime/CIME/Tools/list_e3sm_tests cime_tiny | grep ERS | cut -f 4 -d . | cut -f 2 -d _ )
    # mach_file=${e3sm_root}/cime_config/machines/cmake_macros/${E3SM_COMP}_${E3SM_MACH}.cmake
    # mach_file=${e3sm_root}/components/homme/cmake/machineFiles/${E3SM_MACH}.cmake
    cmake -C ${mach_file} -DBUILD_HOMME_WITHOUT_PIOLIBRARY=OFF -DPREQX_PLEV=26 ${e3sm_root}/components/homme/
    make -j4 homme_tool
    ```

1. run homme_tool
    ```shell
    cat <<EOF > input.nl
    &ctl_nl
    ne = ${NE}
    mesh_file = "${grid_root}/exodus_ne${NE}.g"
    smooth_phis_p2filt = 0
    smooth_phis_numcycle = 12  # v3 uses 6 for less smoothing
    smooth_phis_nudt = 4e-16
    hypervis_scaling = 2
    se_ftype = 2 ! actually output NPHYS; overloaded use of ftype
    /
    &vert_nl
    /
    &analysis_nl
    tool = 'topo_pgn_to_smoothed'
    infilenames = '${topo_root}/USGS-topo_${NE}np4.nc', '${topo_root}/USGS-topo_${NE}np4_smoothed'
    /
    EOF

    mpirun -np 8 ~/E3SM/E3SM_SRC4/components/homme/src/tool/homme_tool < input.nl
    ```

1. Compute phi_s on the np4 grid

1. Use homme_tool to apply dycore specific smoothing on GLL grid.

1. Compute SGH, SGH30, LANDFRAC, and LANDM_COSLAT on the pg2 grid, using the pg2 phi_s data.

1. Append the GLL phi_s data to the output of step 4.





