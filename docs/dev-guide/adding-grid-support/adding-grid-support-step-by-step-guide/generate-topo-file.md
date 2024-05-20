# Generate a Topography File

Topography needs to be interpolated from a high resolution USGS file, and then doctored up a bit to allow the model to run stably with the new topography. The tool chain used to compute topography is documented in the following paper:  

[P.H. Lauritzen, J.T. Bacmeister, P.F. Callaghan, M. Taylor,  NCAR_Topo (v1.0): NCAR global model topography generation software for unstructured grids, Geosci. Model Dev., 8, 3975-3986, 2015.](https://www.geosci-model-dev.net/8/3975/2015/)

## Input Topography Data

Traditionally the topography generation for E3SM start with **USGS-topo-cube3000.nc**, which is included in the [E3SM inputdata repository]. This is a high-resolution topography dataset on a 3km cubed sphere grid derived from 1 km resolution source data.  This file is located in the [CESM inputdata server here](https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/atm/cam/hrtopo/).

For target resolutions of 3 km or finer it is recommended to use **USGS-topo-cube12000.nc**, which was created by Jishi Zhang in 2024. This file has a resolution of 750m created from a  500m/250m USGS GMTED2010 source DEM dataset (see [here](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/4189520033/800m+cubed+topo+generation+from+GMTED2010+15s+DEM) for more information).

For testing the topography workflow the mapping between the ne3000 data is much too burdensome, so a ne90pg1 (i.e. 1-degree) version of this data was created to allow efficient testing. This file can be found in the inputdata repository at `${DIN_LOC_ROOT}/atm/cam/hrtopo/USGS-topo-cube90.nc`.

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

!!! NOTE
    homme_tool is not routinely tested on all supported machines, and the build can be broken without anyone noticing. If you encounter problems building homme_tool please reach out on the e3sm_help slack channel and include a detailed description of the error and the commands you used to produce the error.

<!-- disable certain linter checks here for more readable nested markdown  -->
<!-- markdownlint-disable  MD007 --> <!-- ul-indent -->
<!-- markdownlint-disable  MD033 --> <!-- no-inline-html -->

## Step-by-Step Topography Generation

1. **Set some helpful environmental variables used in commands below**

    ```shell
    e3sm_root=<path to E3SM source>
    grid_root=<path to write grid files>
    map_root=<path to write map files>
    topo_root=<path to write new topo files>
    DIN_LOC_ROOT=<path to E3SM inputdata>
    NE_SRC=<ne value for source grid>
    NE_DST=<ne value for target grid>
    ```

    Example path settings for Perlmutter (NERSC):

    ```shell
    grid_root=${SCRATCH}/e3sm_scratch/files_grid
    map_root=${SCRATCH}/e3sm_scratch/files_map
    topo_root=${SCRATCH}/e3sm_scratch/files_topo
    DIN_LOC_ROOT=/global/cfs/cdirs/e3sm/inputdata
    ```

    Example path settings for Chrysalis (ANL/LCRC):

    ```shell
    grid_root=/lcrc/group/e3sm/${USER}/scratch/chrys/files_grid
    map_root=/lcrc/group/e3sm/${USER}/scratch/chrys/files_map
    topo_root=/lcrc/group/e3sm/${USER}/scratch/chrys/files_topo
    DIN_LOC_ROOT=/lcrc/group/e3sm/data/inputdata
    ```

    Make sure the directories exist:

    ```shell
    mkdir -p ${grid_root} ${map_root} ${topo_root}
    ```

1. **Specify all topo and map file paths for consistency**

    ```shell
    # topo_file_0=${DIN_LOC_ROOT}/atm/cam/hrtopo/USGS-topo-cube${NE_SRC}.nc
    # topo_file_1=${topo_root}/USGS-topo_tmp_ne${NE_DST}np4.nc
    # topo_file_2=${topo_root}/USGS-topo_tmp_ne${NE_DST}np4_phis.nc
    # topo_file_3=${topo_root}/USGS-topo_tmp_ne${NE_DST}np4_smoothedx6t.nc
    # topo_file_4=${topo_root}/USGS-topo_tmp_ne${NE_DST}np4_smoothedx6t_ne${NE_SRC}pg1.nc
    # topo_file_5=${topo_root}/USGS-topo_tmp_ne${NE_DST}np4_smoothedx6t_ne${NE_SRC}pg1_anomalies.nc
    # topo_file_6=${topo_root}/USGS-topo_tmp_ne${NE_DST}np4_smoothedx6t_anomalies.nc
    # topo_file_7=${topo_root}/USGS-topo_tmp_ne${NE_DST}np4_smoothedx6t_sgh.nc
    # topo_file_8=${topo_root}/USGS-topo_ne${NE_DST}np4_smoothedx6t.nc


    topo_file_0=${DIN_LOC_ROOT}/atm/cam/hrtopo/USGS-topo-cube${NE_SRC}.nc
    topo_file_1=${topo_root}/USGS-topo_tmp_ne${NE_SRC}pg1.nc
    topo_file_2=${topo_root}/USGS-topo_tmp_ne${NE_DST}np4.nc
    topo_file_3=${topo_root}/USGS-topo_tmp_ne${NE_DST}np4_smoothedx6t.nc
    topo_file_4=${topo_root}/USGS-topo_tmp_ne${NE_DST}np4_smoothedx6t_ne${NE_SRC}pg1.nc
    topo_file_5=${topo_root}/USGS-topo_tmp_ne${NE_DST}np4_smoothedx6t_ne${NE_SRC}pg1_anomalies.nc
    topo_file_6=${topo_root}/USGS-topo_tmp_ne${NE_DST}np4_smoothedx6t_anomalies.nc
    topo_file_7=${topo_root}/USGS-topo_tmp_ne${NE_DST}np4_smoothedx6t_sgh.nc
    topo_file_8=${topo_root}/USGS-topo_ne${NE_DST}np4_smoothedx6t.nc

    map_file_src_to_np4=${map_root}/map_ne${NE_SRC}pg1_to_ne${NE_DST}np4_fv2se_flx.nc
    map_file_src_to_pg2=${map_root}/map_ne${NE_SRC}pg1_to_ne${NE_DST}pg2_traave.nc
    map_file_pg2_to_src=${map_root}/map_ne${NE_DST}pg2_to_ne${NE_SRC}pg1_traave.nc
    map_file_np4_to_pg2=${map_root}/map_ne${NE_DST}np4_to_ne${NE_DST}pg2_se2fv_flx.nc
    ```

1. **Create grid and map files**

    1. source the unified env

        Perlmutter (NERSC):

        ```shell
        source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
        ```

        Chrysalis (ANL/LCRC):

        ```shell
        source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_chrysalis.sh
        ```

    1. Create grid files

        ```shell
        # Grid for source high res topo
        GenerateCSMesh --alt --res ${NE_SRC}  --file ${grid_root}/exodus_ne${NE_SRC}.g
        ConvertMeshToSCRIP --in ${grid_root}/exodus_ne${NE_SRC}.g  --out ${grid_root}/scrip_ne${NE_SRC}pg1.nc
        
        # Grid for target EAM grid
        GenerateCSMesh --alt --res ${NE_DST} --file ${grid_root}/exodus_ne${NE_DST}.g
        GenerateVolumetricMesh --in ${grid_root}/exodus_ne${NE_DST}.g --out ${grid_root}/exodus_ne${NE_DST}pg2.g --np 2 --uniform
        ConvertMeshToSCRIP --in ${grid_root}/exodus_ne${NE_DST}pg2.g --out ${grid_root}/scrip_ne${NE_DST}pg2.nc
        ```

    1. Create map files

        !!!WARNING
            This can take a long time - several hours for each map file in some cases

        ```shell
        # from source to target np4 
        ncremap -a fv2se_flx -5 --src_grd=${grid_root}/scrip_ne${NE_SRC}pg1.nc  --dst_grd=${grid_root}/exodus_ne${NE_DST}.g --map_file=${map_file_src_to_np4}
        
        # from source to target pg2
        ncremap -a traave -5 --src_grd=${grid_root}/scrip_ne${NE_SRC}pg1.nc  --dst_grd=${grid_root}/scrip_ne${NE_DST}pg2.nc --map_file=${map_file_src_to_pg2}
        
        # from target to source (needed for calculating sub-grid anomalies on target grid)
        ncremap -a traave -5 --src_grd=${grid_root}/scrip_ne${NE_DST}pg2.nc  --dst_grd=${grid_root}/scrip_ne${NE_SRC}pg1.nc --map_file=${
            map_file_pg2_to_src}

        # from target np4 to target pg2
        ncremap -a se2fv_flx -5 --src_grd=${grid_root}/exodus_ne${NE_DST}.g  --dst_grd=${grid_root}/scrip_ne${NE_DST}pg2.nc --map_file=${map_file_np4_to_pg2}
        ```

1. **Create new topograpy data on target grid**

    ```shell
    # # Map high-res topo to target np4 grid
    # ncremap -m ${map_file_src_to_np4} -i ${topo_file_0} -o ${topo_file_1}

    # # Compute phi_s on the np4 grid
    # ncap2 -s 'PHIS=terr*9.80616' ${topo_file_1} ${topo_file_2}

    # # rename the column dimension to be "ncol"
    # ncrename -d grid_size,ncol ${topo_file_2}

    # Compute phi_s on the source np4 grid
    ncap2 -s 'PHIS=terr*9.80616' ${topo_file_0} ${topo_file_1}

    # rename the column dimension to be "ncol"
    ncrename -d grid_size,ncol ${topo_file_1}

    # Map high-res topo to target np4 grid
    ncremap -m ${map_file_src_to_np4} -i ${topo_file_1} -o ${topo_file_2}
    ```

1. **Use homme_tool to smooth topography**

    1. Build homme_tool

        The build process requires the user to select the appropriate cmake file that contains machine-specific settings.

        ```shell
        # Set the machine specific environment
        cd ${e3sm_root}/components/homme
        eval $(${e3sm_root}/cime/CIME/Tools/get_case_env)

        # Specify machine configuration file
        mach_file=${e3sm_root}/components/homme/cmake/machineFiles/perlmutter-gnu.cmake
        # mach_file=${e3sm_root}/components/homme/cmake/machineFiles/chrysalis.cmake

        cmake -C ${mach_file} \
        -DBUILD_HOMME_THETA_KOKKOS=FALSE \
        -DBUILD_HOMME_PREQX_KOKKOS=FALSE \
        -DHOMME_ENABLE_COMPOSE=FALSE \
        -DHOMME_BUILD_EXECS=FALSE \
        -DBUILD_HOMME_TOOL=TRUE \
        -DPREQX_PLEV=26 \
        ${e3sm_root}/components/homme

        make -j4 homme_tool
        ```

    1. run homme_tool
        <!-- ne = ${NE_DST} -->
        <!-- mesh_file = "${grid_root}/exodus_ne${NE_DST}.g" -->
        ```shell
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
        infilenames = '${topo_file_2}', '${topo_file_3}'
        /
        EOF

        mpirun -np 8 ${e3sm_root}/components/homme/src/tool/homme_tool < input.nl

        # rename output file to remove "1.nc" suffix
        mv ${topo_file_3}1.nc ${topo_file_3}
        ```

1. **Compute SGH and SGH30 on the pg2 grid, using the pg2 phi_s data**

    ```shell
    # Remap smoothed data back to source grid
    ncremap -v PHIS -m ${map_file_pg2_to_src} -i ${topo_file_3} -o ${topo_file_4}
    
    # Calculate anomalies on source grid
    ncdiff -O ${topo_file_1} ${topo_file_4} ${topo_file_5}
    
    # Square the anomaly values
    ncap2 -O -s 'PHIS_ANOM_SQ=PHIS^2' ${topo_file_5} ${topo_file_5}
    
    # Remap squared anomalies back to target pg2 grid
    ncremap -v PHIS_ANOM_SQ -m ${map_file_src_to_pg2} -i ${topo_file_5} -o ${topo_file_6}
    
    # Take the square root of the remapped to get standard deviation (SGH)
    ncap2 -O -s 'SGH=sqrt(PHIS_ANOM_SQ)' ${topo_file_6} ${topo_file_7}
    ```

1. **Put all quantities into the final topo file**

    ```shell
    
    # Create final topo file starting with smoothed PHIS data
    cp ${topo_file_3} ${topo_file_8}

    # Append the SGH data
    ncks -A -v SGH ${topo_file_7} ${topo_file_8}

    # rename GLL coordinate to ncol_g
    ncrename -d ncol,ncol_d ${topo_file_2}

    # Map np4 LANDFRAC and SGH30 to the target pg2 grid
    ncremap -v SGH30,LANDFRAC -m ${map_file_np4_to_pg2} -i ${topo_file_2} -o ${????}

    # Append the pg2 LANDFRAC and SGH30 data to final file
    ncks -A -v SGH30,LANDFRAC --hdr_pad=100000 ${?????} ${topo_file_8}
    
    ```

1. **Clean up temporary files**

    ```shell
    rm ${topo_root}/USGS-topo_tmp_*
    ```

## Batch script to streamline all steps

Running through all the steps above can be tedious and time-consuming. The batch script below includes all these steps as well as example Slurm batch directives for running on Perlmutter CPU nodes(NERSC). The only step that is omitted is building homme_tool, since its better to do this manually in case problems arise. 

Here's a check list of things to do before submitting this script:

- Build `homme_tool`
- Update allocation code (i.e. `--account`
- Update batch job wallclock time
- Update paths at the top of the batch script
- Comment out any sections that were completed in advance (i.e. grid file creation)

To submit the slurm batch job use `sbatch batch_topo_slurm.sh`

<details>
    <summary>batch_topo_slurm.sh</summary>
    ```shell
    #!/bin/bash
    #SBATCH --account=### PUT ALLOCATION CODE HERE ###
    #SBATCH --constraint=cpu
    #SBATCH --qos=regular
    #SBATCH --job-name=generate_topo
    #SBATCH --output=slurm-%x-%j.out
    #SBATCH --time=24:00:00
    #SBATCH --nodes=1
    #---------------------------------------------------------------------------
    # Stop execution on any error
    set -e
    #---------------------------------------------------------------------------
    e3sm_root=<path to E3SM source>
    grid_root=<path to write grid files>
    map_root=<path to write map files>
    topo_root=<path to write new topo files>
    DIN_LOC_ROOT=<path to E3SM inputdata>
    NE_SRC=<ne value for source grid>
    NE_DST=<ne value for target grid>
    #---------------------------------------------------------------------------
    topo_file_0=${DIN_LOC_ROOT}/atm/cam/hrtopo/USGS-topo-cube${NE_SRC}.nc
    topo_file_1=${topo_root}/USGS-topo_ne${NE_DST}np4.nc
    topo_file_2=${topo_root}/USGS-topo_ne${NE_DST}np4_phis.nc
    topo_file_3=${topo_root}/USGS-topo_ne${NE_DST}np4_smoothed.nc
    topo_file_4=${topo_root}/USGS-topo_ne${NE_DST}np4_smoothed_ne${NE_SRC}pg1.nc
    topo_file_5=${topo_root}/USGS-topo_ne${NE_DST}np4_smoothed_ne${NE_SRC}pg1_anomalies.nc
    topo_file_6=${topo_root}/USGS-topo_ne${NE_DST}np4_smoothed_anomalies.nc
    #---------------------------------------------------------------------------
    # print some useful things to the log file
    echo e3sm_root    = $e3sm_root
    echo grid_root    = $grid_root
    echo map_root     = $map_root
    echo topo_root    = $topo_root
    echo DIN_LOC_ROOT = $DIN_LOC_ROOT
    echo topo_file_0  = $topo_file_0
    echo topo_file_1  = $topo_file_1
    echo topo_file_2  = $topo_file_2
    echo topo_file_3  = $topo_file_3
    echo topo_file_4  = $topo_file_4
    echo topo_file_5  = $topo_file_5
    echo topo_file_6  = $topo_file_6
    #---------------------------------------------------------------------------
    source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
    #---------------------------------------------------------------------------
    # Create grid files for the input high res topo
    GenerateCSMesh --alt --res ${NE_SRC}  --file ${grid_root}/exodus_ne${NE_SRC}.g
    ConvertMeshToSCRIP --in ${grid_root}/exodus_ne${NE_SRC}.g  --out ${grid_root}/scrip_ne${NE_SRC}pg1.nc
    #---------------------------------------------------------------------------
    # Create grid files target EAM grid
    GenerateCSMesh --alt --res ${NE_DST} --file ${grid_root}/exodus_ne${NE_DST}.g
    GenerateVolumetricMesh --in ${grid_root}/exodus_ne${NE_DST}.g --out ${grid_root}/exodus_ne${NE_DST}pg2.g --np 2 --uniform
    ConvertMeshToSCRIP --in ${grid_root}/exodus_ne${NE_DST}pg2.g --out ${grid_root}/scrip_ne${NE_DST}pg2.nc
    #---------------------------------------------------------------------------
    # Create map files
    # from source to target np4 
    ncremap -a fv2se_flx -5 --src_grd=${grid_root}/scrip_ne${NE_SRC}pg1.nc  --dst_grd=${grid_root}/exodus_ne${NE_DST}.g --map_file=${map_file_src_to_np4}
    # from source to target pg2
    ncremap -a traave -5 --src_grd=${grid_root}/scrip_ne${NE_SRC}pg1.nc  --dst_grd=${grid_root}/scrip_ne${NE_DST}pg2.nc --map_file=${map_file_src_to_pg2}
    # from target to source (needed for calculating sub-grid anomalies on target grid)
    ncremap -a traave -5 --src_grd=${grid_root}/scrip_ne${NE_DST}pg2.nc  --dst_grd=${grid_root}/scrip_ne${NE_SRC}pg1.nc --map_file=${map_file_pg2_to_src}
    #---------------------------------------------------------------------------
    # Map high-res topo to target np4 grid
    ncremap -m ${map_file_src_to_np4} -i ${topo_file_0} -o ${topo_file_1}
    #---------------------------------------------------------------------------
    # Compute phi_s on the np4 grid
    ncap2 -s 'PHIS=terr*9.80616' ${topo_file_1} ${topo_file_2}
    #---------------------------------------------------------------------------
    # Use homme_tool to smooth topography
    cat <<EOF > input.nl
    &ctl_nl
    ne = ${NE_DST}
    mesh_file = "${grid_root}/exodus_ne${NE_DST}.g"
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
    infilenames = '${topo_file_2}', '${topo_file_3::-3}'
    /
    EOF
    mpirun -np 8 ${e3sm_root}/components/homme/src/tool/homme_tool < input.nl
    #---------------------------------------------------------------------------
    # Compute SGH and SGH30 on the pg2 grid, using the pg2 phi_s data
    ncremap -v PHIS -m ${map_file_pg2_to_src} -i ${topo_file_3} -o ${topo_file_4}
    ncks -A ${topo_file_0} ${topo_file_4}
    ncdiff ${topo_file_0} ${topo_file_4} ${topo_file_5}
    ncremap -m ${map_file_src_to_pg2} -i ${topo_file_0} -o ${topo_file_6}
    #---------------------------------------------------------------------------
    #---------------------------------------------------------------------------
    #---------------------------------------------------------------------------
    ```
</details>
