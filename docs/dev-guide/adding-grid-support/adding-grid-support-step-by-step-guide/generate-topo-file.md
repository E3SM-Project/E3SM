# Generate a Topography File

<!-- disable certain linter checks here for more readable nested markdown  -->
<!-- markdownlint-disable  MD007 --> <!-- ul-indent -->
<!-- markdownlint-disable  MD033 --> <!-- no-inline-html -->

Topography needs to be interpolated from a high resolution dataset, and then doctored a bit to allow the model to run stably with the new topography. More information can be found in the following paper:  

[P.H. Lauritzen, J.T. Bacmeister, P.F. Callaghan, M. Taylor,  NCAR_Topo (v1.0): NCAR global model topography generation software for unstructured grids, Geosci. Model Dev., 8, 3975-3986, 2015.](https://www.geosci-model-dev.net/8/3975/2015/)

Typically, input topography data generation for E3SM starts with a high resolution source dataset (`USGS-topo-cube3000.nc`). This is a high-resolution topography dataset on a 3km cubed sphere grid derived from 1 km resolution source data.  This file is located in the [CESM inputdata server here](https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/atm/cam/hrtopo/).

For target resolutions of 3 km or finer it is recommended to use an even higher resolution source dataset (`USGS-topo-cube12000.nc`), which was created by Jishi Zhang in 2024. This file has a resolution of 750m created from a  500m/250m USGS GMTED2010 source DEM dataset (see [here](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/4189520033/800m+cubed+topo+generation+from+GMTED2010+15s+DEM) for more information).

For testing the topography workflow the mapping between the ne3000 data is too burdensome, so a ne90pg1 (i.e. 1-degree) version of this data was created to allow efficient testing. This file can be found on various supported machines at `${DIN_LOC_ROOT}/atm/cam/hrtopo/USGS-topo-cube90.nc`.

## Data Processing Requirements

The workflow to generate topography data for the atmosphere model must address the following requirements:

- The dycore needs surface geopotential (`phi_s`) at GLL nodes (np4)
- The physics needs various things at FV cell centers (pg2)
    - surface geopotential (`phi_s`)
    - sub-grid variance of topography (`SGH` and `SGH30`)
    - Land area fraction (`LANDFRAC`)
- The np4 to pg2 map of `phi_s` data must be equal to the FV `phi_s` data
- HOMME's smoothing operator must be used on `phi_s` for stability

## Topography Smoothing

Smoothing of the input surface geopotential (`phi_s`) is an essential step to ensure numerical stability of the atmospheric dynamics, but the smoothing much be done in a way that is consistent with the internal Laplacian used by the HOMME dycor. To accomplish this we use `homme_tool`, which is a standalone build of the HOMME dycor.

### Building homme_tool

!!! NOTE
    homme_tool is not routinely tested on all supported machines, and the build can be broken without anyone noticing. If you encounter problems building homme_tool please reach out on the e3sm_help slack channel and include a detailed description of the error and the commands you used to produce the error.

Building `homme_tool` is a critical preliminary step to the topography generation workflow described below. The build process requires the user to select the appropriate cmake file that contains machine-specific settings (see `mach_file` below).

```shell
# Set the machine specific environment
cd ${e3sm_root}/components/homme

# load the appropriate machine environment
eval $(${e3sm_root}/cime/CIME/Tools/get_case_env)

# Specify machine configuration file
mach_file=${e3sm_root}/components/homme/cmake/machineFiles/pm-cpu.cmake
# mach_file=${e3sm_root}/components/homme/cmake/machineFiles/chrysalis.cmake

cmake -C ${mach_file} \
-DBUILD_HOMME_THETA_KOKKOS=FALSE \
-DBUILD_HOMME_PREQX_KOKKOS=FALSE \
-DHOMME_ENABLE_COMPOSE=FALSE \
-DHOMME_BUILD_EXECS=FALSE \
-DBUILD_HOMME_TOOL=TRUE \
-DBUILD_HOMME_WITHOUT_PIOLIBRARY=FALSE \
-DPREQX_PLEV=26 \
${e3sm_root}/components/homme

make -j4 homme_tool
```

## Sub-Grid Topography Variations

Certain physics calculations in the atmosphere require a characterization of the unresolved topography. This is provided in the input data as `SGH` and `SGH30`. These quantities need to be calculated by the `cube_to_target` tool that is included in the E3SM source  code. Similar to `homme_tool`, this tool can often be broken when the user goes to use it. Luckily, this tool is much simpler with fewer dependencies. Nevertheless, it is good to make sure the tools builds successfully before creating a new topography file.

### Building cube_to_target

The following commands were working on both Perlmutter/NERSC and Chrysalis/LCRC machines as of 2025.

```shell
cd ${e3sm_root}/components/eam/tools/topo_tool/cube_to_target

eval $(${e3sm_root}/cime/CIME/Tools/get_case_env)

make
```

!!! NOTE
    You can safely ignore compiler warnings that look like `Warning: Rank mismatch...`. These are a result of how we interface with the netcdf fortran routines, but fixing the warnings by switching from `#include <netcdf.inc>` to `use netcdf` can lead to other problems on certain machines.

## Step-by-Step Topography Generation

!!! NOTE
    Copying and pasting the relevant steps below is an easy way to step through the blocks of commands. Alternatively, a batch script is provided below that can be used to execute all steps at once in a batch job.

1. ### **Activate the E3SM Unified Env**

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

1. ### **Specify Source and Target Resolution**

    We will environement variables to specify the source and target grid resolutions based on the "ne" value of the cubed sphere grid ("ne" is the number of elements along a cube edge). A Typical use case will be mapping data from `ne3000pg1` to a chosen target resolution, in this case `ne30`:

    ```shell
    NE_SRC=3000
    NE_DST=30
    ```

    For testing use a special `ne90pg1` dataset:

    ```shell
    NE_SRC=90
    NE_DST=4
    ```

    !!! NOTE
        For grids with regional refinement (RRM) there is no corresponding "ne" value - so a slightly modified workflow is needed. This entails modifying the grid file generation step, and then modifying topo and map file paths to accomodate the new RRM grid name.

1. ### **Specify File Paths**

    First we need to set some enviroment variables that point to various "root" directories where we will be writing and/or reading files. It

    ```shell
    e3sm_root=?    # path to E3SM source
    grid_root=?    # path to write grid files
    map_root=?     # path to write map files
    topo_root=?    # path to write new topo files
    DIN_LOC_ROOT=? # path to E3SM inputdata
    ```

    Example path settings:

    - Perlmutter (NERSC):

        ```shell
        e3sm_root=${SCRATCH}/tmp_e3sm_src # make sure this contains an up-to-date clone of E3SM
        grid_root=${SCRATCH}/e3sm_scratch/files_grid
        map_root=${SCRATCH}/e3sm_scratch/files_map
        topo_root=${SCRATCH}/e3sm_scratch/files_topo
        DIN_LOC_ROOT=/global/cfs/cdirs/e3sm/inputdata
        ```

    - Chrysalis (ANL/LCRC):

        ```shell
        SCRATCH=/lcrc/group/e3sm/${USER}/scratch/chrys
        e3sm_root=${SCRATCH}/tmp_e3sm_src # make sure this contains an up-to-date clone of E3SM
        grid_root=${SCRATCH}/files_grid
        map_root=${SCRATCH}/files_map
        topo_root=${SCRATCH}/files_topo
        DIN_LOC_ROOT=/lcrc/group/e3sm/data/inputdata
        ```

    Make sure the root directories exist:

    ```shell
    mkdir -p ${grid_root} ${map_root} ${topo_root}
    ls -ld ${grid_root} ${map_root} ${topo_root} ${e3sm_root} ${DIN_LOC_ROOT}
    ```

    Now specify all the files that we will need for, including map files and temporary topo data.

    ```shell
    timestamp=$(date +%Y%m%d)
    topo_file_0=${DIN_LOC_ROOT}/atm/cam/hrtopo/USGS-topo-cube${NE_SRC}.nc
    topo_file_1=${topo_root}/tmp_USGS-topo_ne${NE_DST}np4.nc
    topo_file_2=${topo_root}/tmp_USGS-topo_ne${NE_DST}np4_smoothedx6t.nc
    topo_file_3=${topo_root}/USGS-topo_ne${NE_DST}np4_smoothedx6t_${timestamp}.nc

    map_file_src_to_np4=${map_root}/map_ne${NE_SRC}pg1_to_ne${NE_DST}np4_fv2se_flx.nc
    ```

1. ### **Create Grid Files**

    ```shell
    # Grid for source high res topo
    GenerateCSMesh --alt --res ${NE_SRC}  --file ${grid_root}/exodus_ne${NE_SRC}.g
    ConvertMeshToSCRIP --in ${grid_root}/exodus_ne${NE_SRC}.g  --out ${grid_root}/scrip_ne${NE_SRC}pg1.nc

    # Grid for target EAM grid
    GenerateCSMesh --alt --res ${NE_DST} --file ${grid_root}/exodus_ne${NE_DST}.g
    GenerateVolumetricMesh --in ${grid_root}/exodus_ne${NE_DST}.g --out ${grid_root}/exodus_ne${NE_DST}pg2.g --np 2 --uniform
    ConvertMeshToSCRIP --in ${grid_root}/exodus_ne${NE_DST}pg2.g --out ${grid_root}/scrip_ne${NE_DST}pg2.nc
    ```

1. ### **Create Map Files**

    !!!WARNING
        This step can potentially take a very long time - several hours for each map file in some cases!
        The use of `--mpi_nbr=32` will leverage parallelization via `mbtemptest` to dramatically reduce the time to generate mapping files, but this may not work on all machines. Feel free to experiment with a different number of tasks to optimize for the machine you are using. You can also drop this argument entirely if it is causing problems, but be sure to plan for hours of execution time.

    !!!NOTE
        If you see an error from these `ncremap` commands that looks like: `srun: error: Job request does not match any supported policy.` and you are on a login node then you will need to either launch an interactive compute node or use a batch job script.

    ```shell
    # Create map from source to target np4 
    ncremap ${MAP_ARGS} -a fv2se_flx \
    --src_grd=${grid_root}/scrip_ne${NE_SRC}pg1.nc  \
    --dst_grd=${grid_root}/exodus_ne${NE_DST}.g \
    --map_file=${map_file_src_to_np4} \
    --tmp_dir=${map_root}
    ```

1. ### **Remap Topograpy**

    The first command here will essentially create a copy of the source data, which is needed later on when calculating SGH.

    ```shell
    # Map high-res topo to target np4 grid
    ncremap -m ${map_file_src_to_np4} -i ${topo_file_0} -o ${topo_file_1}

    # Compute phi_s on the target np4 grid
    ncap2 -O -s 'PHIS=terr*9.80616' ${topo_file_1} ${topo_file_1}

    # rename the column dimension to be "ncol"
    ncrename -d grid_size,ncol ${topo_file_1}
    ```

1. ### **Apply Smoothing**

    Make sure `homme_tool` has been built following the [instructions above](#building-homme_tool).

    ```shell
    cd ${e3sm_root}/components/homme
    eval $(${e3sm_root}/cime/CIME/Tools/get_case_env)

    cat <<EOF > input.nl
    &ctl_nl
    mesh_file = "${grid_root}/exodus_ne${NE_DST}.g"
    smooth_phis_p2filt = 0
    smooth_phis_numcycle = 6 ! v2/v3 uses 12/6 for more/less smoothing
    smooth_phis_nudt = 4e-16
    hypervis_scaling = 2
    se_ftype = 2 ! actually output NPHYS; overloaded use of ftype
    /
    &vert_nl
    /
    &analysis_nl
    tool = 'topo_pgn_to_smoothed'
    infilenames = '${topo_file_1}', '${topo_file_2}'
    /
    EOF

    mpirun -np 8 ${e3sm_root}/components/homme/src/tool/homme_tool < input.nl

    # rename output file to remove "1.nc" suffix
    mv ${topo_file_2}1.nc ${topo_file_2}
    ```

1. ### **Compute SGH**

    ```shell
    ${e3sm_root}/components/eam/tools/topo_tool/cube_to_target/cube_to_target \
      --target-grid ${grid_root}/scrip_ne${NE_DST}pg2.nc \
      --input-topography ${topo_file_0} \
      --smoothed-topography ${topo_file_2} \
      --output-topography ${topo_file_3}

    # Append the GLL phi_s data to the output of step 4.
    ncks -A ${topo_file_2} ${topo_file_3}
    ```

1. ### **Clean up Temporary Files**

    ```shell
    rm ${topo_root}/tmp_USGS-topo_ne${NE_DST}*
    ```

## Batch script to streamline all steps

Running through all the steps above can be tedious and time-consuming. The batch scripts below include all these steps as well as batch system directives. The only step that is omitted is building homme_tool, since its better to do this manually in case problems arise.

Here is a check list of things to do before submitting this script:

- Build `homme_tool`
- Build `cube_to_target`
- Update batch allocation code (i.e. `--account`)
- Update batch job wallclock time (i.e. `--time`)
- Update paths at the top of the batch script
- Comment out any sections that were completed in advance (i.e. grid and map file creation)

To submit the slurm batch job use `sbatch <script>`

<details>
    <summary>batch_topo_slurm_lcrc.sh</summary>
    ```shell
    #!/bin/bash
    #SBATCH --account=e3sm
    #SBATCH --job-name=generate_topo
    #SBATCH --output=slurm-%x-%j.out
    #SBATCH --time=24:00:00
    #SBATCH --nodes=1
    #SBATCH --mail-type=END,FAIL
    #---------------------------------------------------------------------------
    # Make sure all these lines are correct for the machine
    source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_chrysalis.sh
    # Specify source and target resolutions
    NE_SRC=3000 ; NE_DST=30
    # NE_SRC=90 ; NE_DST=30 # low-res grid combination for testing
    # Specify time stamp for creation date
    timestamp=$(date +%Y%m%d)
    # Specify E3SM source code path - preferably a fresh clone
    SCRATCH=/lcrc/group/e3sm/${USER}/scratch/chrys
    e3sm_root=${SCRATCH}/tmp_e3sm_src
    # Specify root paths
    grid_root=${SCRATCH}/files_grid
    map_root=${SCRATCH}/files_map
    topo_root=${SCRATCH}/files_topo
    DIN_LOC_ROOT=/lcrc/group/e3sm/data/inputdata
    # argument for ncremap to select TempestRemap or mbtempest backend
    MAP_ARGS=
    # MAP_ARGS+="--mpi_nbr=32"
    #---------------------------------------------------------------------------
    # Stop script execution on error
    set -e
    # ANSI color codes for highlighting terminal output
    RED='\033[0;31m' ; GRN='\033[0;32m' CYN='\033[0;36m' ; NC='\033[0m'
    # start timer for entire script
    start=`date +%s`
    #---------------------------------------------------------------------------
    # Specify topo file names - including temporary files that will be deleted
    topo_file_0=${DIN_LOC_ROOT}/atm/cam/hrtopo/USGS-topo-cube${NE_SRC}.nc
    topo_file_1=${topo_root}/tmp_USGS-topo_ne${NE_DST}np4.nc
    topo_file_2=${topo_root}/tmp_USGS-topo_ne${NE_DST}np4_smoothedx6t.nc
    topo_file_3=${topo_root}/USGS-topo_ne${NE_DST}np4_smoothedx6t_${timestamp}.nc
    # Specify map file name
    map_file_src_to_np4=${map_root}/map_ne${NE_SRC}pg1_to_ne${NE_DST}np4_fv2se_flx.nc
    #---------------------------------------------------------------------------  
    # print some useful things
    echo --------------------------------------------------------------------------------
    echo --------------------------------------------------------------------------------
    echo
    echo   NE_SRC              = $NE_SRC
    echo   NE_DST              = $NE_DST
    echo
    echo   e3sm_root           = $e3sm_root
    echo   grid_root           = $grid_root
    echo   map_root            = $map_root
    echo   topo_root           = $topo_root
    echo   DIN_LOC_ROOT        = $DIN_LOC_ROOT
    echo
    echo   topo_file_0         = $topo_file_0
    echo   topo_file_1         = $topo_file_1
    echo   topo_file_2         = $topo_file_2
    echo   topo_file_3         = $topo_file_3
    echo
    echo   map_file_src_to_np4 = $map_file_src_to_np4
    echo
    echo --------------------------------------------------------------------------------
    echo --------------------------------------------------------------------------------
    #---------------------------------------------------------------------------
    # Make sure paths exist
    mkdir -p ${grid_root} ${map_root} ${topo_root}
    if [ ! -d ${DIN_LOC_ROOT} ]; then echo -e ${RED}ERROR directory does not exist:${NC} ${DIN_LOC_ROOT} ; fi
    if [ ! -d ${e3sm_root}    ]; then echo -e ${RED}ERROR directory does not exist:${NC} ${e3sm_root} ; fi
    if [ ! -d ${grid_root}    ]; then echo -e ${RED}ERROR directory does not exist:${NC} ${grid_root} ; fi
    if [ ! -d ${map_root}     ]; then echo -e ${RED}ERROR directory does not exist:${NC} ${map_root} ; fi
    if [ ! -d ${topo_root}    ]; then echo -e ${RED}ERROR directory does not exist:${NC} ${topo_root} ; fi
    #---------------------------------------------------------------------------
    # set to echo commands
    set -x
    #---------------------------------------------------------------------------
    # Create grid for source high res topo
    GenerateCSMesh --alt --res ${NE_SRC} --file ${grid_root}/exodus_ne${NE_SRC}.g
    ConvertMeshToSCRIP --in ${grid_root}/exodus_ne${NE_SRC}.g  --out ${grid_root}/scrip_ne${NE_SRC}pg1.nc
    # Create grid for target EAM grid
    GenerateCSMesh --alt --res ${NE_DST} --file ${grid_root}/exodus_ne${NE_DST}.g
    GenerateVolumetricMesh --in ${grid_root}/exodus_ne${NE_DST}.g --out ${grid_root}/exodus_ne${NE_DST}pg2.g --np 2 --uniform
    ConvertMeshToSCRIP --in ${grid_root}/exodus_ne${NE_DST}pg2.g --out ${grid_root}/scrip_ne${NE_DST}pg2.nc
    #---------------------------------------------------------------------------
    # Create map from source to target np4
    time ncremap ${MAP_ARGS} -a fv2se_flx \
      --src_grd=${grid_root}/scrip_ne${NE_SRC}pg1.nc \
      --dst_grd=${grid_root}/exodus_ne${NE_DST}.g \
      --map_file=${map_file_src_to_np4} \
      --tmp_dir=${map_root}
    #---------------------------------------------------------------------------
    # Remap high-res topo to target np4 grid
    ncremap -m ${map_file_src_to_np4} -i ${topo_file_0} -o ${topo_file_1}
    # Compute phi_s on the target np4 grid
    ncap2 -O -s 'PHIS=terr*9.80616' ${topo_file_1} ${topo_file_1}
    # rename the column dimension to be "ncol"
    ncrename -d grid_size,ncol ${topo_file_1}
    #---------------------------------------------------------------------------
    # Apply Smoothing
    cd ${e3sm_root}/components/homme
    eval $(${e3sm_root}/cime/CIME/Tools/get_case_env)
    # Create namelist file for HOMME
    cat <<EOF > input.nl
    &ctl_nl
    mesh_file = "${grid_root}/exodus_ne${NE_DST}.g"
    smooth_phis_p2filt = 0
    smooth_phis_numcycle = 6 ! v2/v3 uses 12/6 for more/less smoothing
    smooth_phis_nudt = 4e-16
    hypervis_scaling = 2
    se_ftype = 2 ! actually output NPHYS; overloaded use of ftype
    /
    &vert_nl
    /
    &analysis_nl
    tool = 'topo_pgn_to_smoothed'
    infilenames = '${topo_file_1}', '${topo_file_2}'
    /
    EOF
    # run homme_tool for topography smoothing
    srun -n 8 ${e3sm_root}/components/homme/src/tool/homme_tool < input.nl
    # rename output file to remove "1.nc" suffix
    mv ${topo_file_2}1.nc ${topo_file_2}
    #---------------------------------------------------------------------------
    # Compute SGH with cube_to_target
    ${e3sm_root}/components/eam/tools/topo_tool/cube_to_target/cube_to_target \
      --target-grid ${grid_root}/scrip_ne${NE_DST}pg2.nc \
      --input-topography ${topo_file_0} \
      --smoothed-topography ${topo_file_2} \
      --output-topography ${topo_file_3}
    # Append the GLL phi_s data to the output
    ncks -A ${topo_file_2} ${topo_file_3}
    #---------------------------------------------------------------------------
    # Clean up Temporary Files
    rm ${topo_root}/tmp_USGS-topo_ne${NE_DST}*
    #---------------------------------------------------------------------------
    # stop echoing commands
    set +x
    #---------------------------------------------------------------------------
    # Check that final topo output file was created
    if [ ! -f ${topo_file_3} ]; then
      echo
      echo -e ${RED} Failed to create topography file - Errors ocurred ${NC}
      echo
    else
      echo
      echo -e ${GRN} Sucessfully created new topography file  ${NC}
      echo $topo_file_3
      echo
    fi
    #---------------------------------------------------------------------------
    # Indicate overall run time for this script
    end=`date +%s`
    runtime_sc=$(( end - start ))
    runtime_mn=$(( runtime_sc/60 ))
    runtime_hr=$(( runtime_mn/60 ))
    echo -e    ${CYN} overall runtime: ${NC} ${runtime_sc} seconds / ${runtime_mn} minutes / ${runtime_hr} hours
    echo
    #---------------------------------------------------------------------------
    ```
</details>

<details>
    <summary>batch_topo_slurm_nersc.sh</summary>
    ```shell
    #!/bin/bash
    #SBATCH --account=e3sm
    #SBATCH --constraint=cpu
    #SBATCH --qos=regular
    #SBATCH --job-name=generate_topo
    #SBATCH --output=slurm-%x-%j.out
    #SBATCH --time=24:00:00
    #SBATCH --nodes=1
    #SBATCH --mail-type=END,FAIL
    #---------------------------------------------------------------------------
    # Make sure all these lines are correct for the machine
    source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
    # Specify source and target resolutions
    NE_SRC=3000 ; NE_DST=30
    # NE_SRC=90 ; NE_DST=30 # low-res grid combination for testing
    # Specify time stamp for creation date
    timestamp=$(date +%Y%m%d)
    # Specify E3SM source code path - preferably a fresh clone
    e3sm_root=${SCRATCH}/tmp_e3sm_src
    # Specify root paths
    grid_root=${SCRATCH}/files_grid
    map_root=${SCRATCH}/files_map
    topo_root=${SCRATCH}/files_topo
    DIN_LOC_ROOT=/global/cfs/cdirs/e3sm/inputdata
    # argument for ncremap to select TempestRemap or mbtempest backend
    MAP_ARGS=
    # MAP_ARGS+="--mpi_nbr=32"
    #---------------------------------------------------------------------------
    # Stop script execution on error
    set -e
    # ANSI color codes for highlighting terminal output
    RED='\033[0;31m' ; GRN='\033[0;32m' CYN='\033[0;36m' ; NC='\033[0m'
    # start timer for entire script
    start=`date +%s`
    #---------------------------------------------------------------------------
    # Specify topo file names - including temporary files that will be deleted
    topo_file_0=${DIN_LOC_ROOT}/atm/cam/hrtopo/USGS-topo-cube${NE_SRC}.nc
    topo_file_1=${topo_root}/tmp_USGS-topo_ne${NE_DST}np4.nc
    topo_file_2=${topo_root}/tmp_USGS-topo_ne${NE_DST}np4_smoothedx6t.nc
    topo_file_3=${topo_root}/USGS-topo_ne${NE_DST}np4_smoothedx6t_${timestamp}.nc
    # Specify map file name
    map_file_src_to_np4=${map_root}/map_ne${NE_SRC}pg1_to_ne${NE_DST}np4_fv2se_flx.nc
    #---------------------------------------------------------------------------  
    # print some useful things
    echo --------------------------------------------------------------------------------
    echo --------------------------------------------------------------------------------
    echo
    echo   NE_SRC              = $NE_SRC
    echo   NE_DST              = $NE_DST
    echo
    echo   e3sm_root           = $e3sm_root
    echo   grid_root           = $grid_root
    echo   map_root            = $map_root
    echo   topo_root           = $topo_root
    echo   DIN_LOC_ROOT        = $DIN_LOC_ROOT
    echo
    echo   topo_file_0         = $topo_file_0
    echo   topo_file_1         = $topo_file_1
    echo   topo_file_2         = $topo_file_2
    echo   topo_file_3         = $topo_file_3
    echo
    echo   map_file_src_to_np4 = $map_file_src_to_np4
    echo
    echo --------------------------------------------------------------------------------
    echo --------------------------------------------------------------------------------
    #---------------------------------------------------------------------------
    # Make sure paths exist
    mkdir -p ${grid_root} ${map_root} ${topo_root}
    if [ ! -d ${DIN_LOC_ROOT} ]; then echo -e ${RED}ERROR directory does not exist:${NC} ${DIN_LOC_ROOT} ; fi
    if [ ! -d ${e3sm_root}    ]; then echo -e ${RED}ERROR directory does not exist:${NC} ${e3sm_root} ; fi
    if [ ! -d ${grid_root}    ]; then echo -e ${RED}ERROR directory does not exist:${NC} ${grid_root} ; fi
    if [ ! -d ${map_root}     ]; then echo -e ${RED}ERROR directory does not exist:${NC} ${map_root} ; fi
    if [ ! -d ${topo_root}    ]; then echo -e ${RED}ERROR directory does not exist:${NC} ${topo_root} ; fi
    #---------------------------------------------------------------------------
    # set to echo commands
    set -x
    #---------------------------------------------------------------------------
    # Create grid for source high res topo
    GenerateCSMesh --alt --res ${NE_SRC} --file ${grid_root}/exodus_ne${NE_SRC}.g
    ConvertMeshToSCRIP --in ${grid_root}/exodus_ne${NE_SRC}.g  --out ${grid_root}/scrip_ne${NE_SRC}pg1.nc
    # Create grid for target EAM grid
    GenerateCSMesh --alt --res ${NE_DST} --file ${grid_root}/exodus_ne${NE_DST}.g
    GenerateVolumetricMesh --in ${grid_root}/exodus_ne${NE_DST}.g --out ${grid_root}/exodus_ne${NE_DST}pg2.g --np 2 --uniform
    ConvertMeshToSCRIP --in ${grid_root}/exodus_ne${NE_DST}pg2.g --out ${grid_root}/scrip_ne${NE_DST}pg2.nc
    #---------------------------------------------------------------------------
    # Create map from source to target np4
    time ncremap ${MAP_ARGS} -a fv2se_flx \
      --src_grd=${grid_root}/scrip_ne${NE_SRC}pg1.nc \
      --dst_grd=${grid_root}/exodus_ne${NE_DST}.g \
      --map_file=${map_file_src_to_np4} \
      --tmp_dir=${map_root}
    #---------------------------------------------------------------------------
    # Remap high-res topo to target np4 grid
    ncremap -m ${map_file_src_to_np4} -i ${topo_file_0} -o ${topo_file_1}
    # Compute phi_s on the target np4 grid
    ncap2 -O -s 'PHIS=terr*9.80616' ${topo_file_1} ${topo_file_1}
    # rename the column dimension to be "ncol"
    ncrename -d grid_size,ncol ${topo_file_1}
    #---------------------------------------------------------------------------
    # Apply Smoothing
    cd ${e3sm_root}/components/homme
    eval $(${e3sm_root}/cime/CIME/Tools/get_case_env)
    # Create namelist file for HOMME
    cat <<EOF > input.nl
    &ctl_nl
    mesh_file = "${grid_root}/exodus_ne${NE_DST}.g"
    smooth_phis_p2filt = 0
    smooth_phis_numcycle = 6 ! v2/v3 uses 12/6 for more/less smoothing
    smooth_phis_nudt = 4e-16
    hypervis_scaling = 2
    se_ftype = 2 ! actually output NPHYS; overloaded use of ftype
    /
    &vert_nl
    /
    &analysis_nl
    tool = 'topo_pgn_to_smoothed'
    infilenames = '${topo_file_1}', '${topo_file_2}'
    /
    EOF
    # run homme_tool for topography smoothing
    srun -n 8 ${e3sm_root}/components/homme/src/tool/homme_tool < input.nl
    # rename output file to remove "1.nc" suffix
    mv ${topo_file_2}1.nc ${topo_file_2}
    #---------------------------------------------------------------------------
    # Compute SGH with cube_to_target
    ${e3sm_root}/components/eam/tools/topo_tool/cube_to_target/cube_to_target \
      --target-grid ${grid_root}/scrip_ne${NE_DST}pg2.nc \
      --input-topography ${topo_file_0} \
      --smoothed-topography ${topo_file_2} \
      --output-topography ${topo_file_3}
    # Append the GLL phi_s data to the output
    ncks -A ${topo_file_2} ${topo_file_3}
    #---------------------------------------------------------------------------
    # Clean up Temporary Files
    rm ${topo_root}/tmp_USGS-topo_ne${NE_DST}*
    #---------------------------------------------------------------------------
    # stop echoing commands
    set +x
    #---------------------------------------------------------------------------
    # Check that final topo output file was created
    if [ ! -f ${topo_file_3} ]; then
      echo
      echo -e ${RED} Failed to create topography file - Errors ocurred ${NC}
      echo
    else
      echo
      echo -e ${GRN} Sucessfully created new topography file  ${NC}
      echo $topo_file_3
      echo
    fi
    #---------------------------------------------------------------------------
    # Indicate overall run time for this script
    end=`date +%s`
    runtime_sc=$(( end - start ))
    runtime_mn=$(( runtime_sc/60 ))
    runtime_hr=$(( runtime_mn/60 ))
    echo -e    ${CYN} overall runtime: ${NC} ${runtime_sc} seconds / ${runtime_mn} minutes / ${runtime_hr} hours
    echo
    #---------------------------------------------------------------------------
    ```
</details>

Back to step-by-step guide for [Adding Support for New Grids](../adding-grid-support-step-by-step-guide.md)
