# Generate a Topography File

Topography needs to be interpolated from a high resolution USGS file, and then doctored up a bit to allow the model to run stably with the new topography. The tool chain used to compute topography is documented in the following paper:  

[P.H. Lauritzen, J.T. Bacmeister, P.F. Callaghan, M. Taylor,  NCAR_Topo (v1.0): NCAR global model topography generation software for unstructured grids, Geosci. Model Dev., 8, 3975-3986, 2015.](https://www.geosci-model-dev.net/8/3975/2015/)

## Input Topography Data

Traditionally the topography generation for E3SM start with **USGS-topo-cube3000.nc**, which is included in the [E3SM inputdata repository]. This is a high-resolution topography dataset on a 3km cubed sphere grid derived from 1 km resolution source data.  This file is located in the [CESM inputdata server here](https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/atm/cam/hrtopo/).

For target resolutions of 3 km or finer it is recommended to use **USGS-topo-cube12000.nc**, which was created by Jishi Zhang in 2024. This file has a resolution of 750m created from a  500m/250m USGS GMTED2010 source DEM dataset (see [here](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/4189520033/800m+cubed+topo+generation+from+GMTED2010+15s+DEM) for more information).

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

## Step-by-Step Topography Generation

1. **Set some helpful environmental variables used in commands below**
  ```shell
  e3sm_root=<path to E3SM source>
  grid_root=<path to write grid files>
  map_root=<path to write map files>
  topo_root=<path to write new topo files>
  DIN_LOC_ROOT=<path to E3SM inputdata>
  NE=<ne value for target grid>
  ```
  Example settings for Perlmutter (NERSC):
  ```shell
  grid_root=${SCRATCH}/e3sm_scratch/files_grid
  map_root=${SCRATCH}/e3sm_scratch/files_map
  topo_root=${SCRATCH}/e3sm_scratch/files_topo
  DIN_LOC_ROOT=/global/cfs/cdirs/e3sm/inputdata
  NE=30
  ```
  Example settings for Chrysalis (ANL/LCRC):
  ```shell
  grid_root=/lcrc/group/e3sm/${USER}/scratch/chrys/files_grid
  map_root=/lcrc/group/e3sm/${USER}/scratch/chrys/files_map
  topo_root=/lcrc/group/e3sm/${USER}/scratch/chrys/files_topo
  DIN_LOC_ROOT=/lcrc/group/e3sm/data/inputdata
  NE=30
  ```
  Make sure the directories exist:
  ```
  mkdir -p ${grid_root} ${map_root} ${topo_root}
  ```

1. **Specify all topo and map file paths for consistency**
  ```shell
  topo_file_0=${DIN_LOC_ROOT}/atm/cam/hrtopo/USGS-topo-cube3000.nc
  topo_file_1=${topo_root}/USGS-topo_ne${NE}np4.nc
  topo_file_2=${topo_root}/USGS-topo_ne${NE}np4_phis.nc
  topo_file_3=${topo_root}/USGS-topo_ne${NE}np4_smoothed.nc
  topo_file_4=${topo_root}/USGS-topo_ne${NE}np4_smoothed_ne3000pg1.nc
  topo_file_5=${topo_root}/USGS-topo_ne${NE}np4_smoothed_ne3000pg1_anomalies.nc
  topo_file_6=${topo_root}/USGS-topo_ne${NE}np4_smoothed_anomalies.nc
  map_file_src_to_np4=${map_root}/map_ne3000pg1_to_ne${NE}np4_fv2se_flx.nc
  map_file_src_to_pg2=${map_root}/map_ne3000pg1_to_ne${NE}pg2_traave.nc
  map_file_pg2_to_src=${map_root}/map_ne${NE}pg2_to_ne3000pg1_traave.nc
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

  1. Create grid files for the input high res topo
    ```shell
    GenerateCSMesh --alt --res 3000  --file ${grid_root}/exodus_ne3000.g
    ConvertMeshToSCRIP --in ${grid_root}/exodus_ne3000.g  --out ${grid_root}/scrip_ne3000pg1.nc
    ```

  1. Create grid files target EAM grid
    ```shell
    GenerateCSMesh --alt --res ${NE} --file ${grid_root}/exodus_ne${NE}.g
    GenerateVolumetricMesh --in ${grid_root}/exodus_ne${NE}.g --out ${grid_root}/exodus_ne${NE}pg2.g --np 2 --uniform
    ConvertMeshToSCRIP --in ${grid_root}/exodus_ne${NE}pg2.g --out ${grid_root}/scrip_ne${NE}pg2.nc
    ```

  1. Create map files
    
    !!!WARNING
      this can take a long time
    1. from source to target np4
      ```shell
      ncremap -a fv2se_flx -5 --src_grd=${grid_root}/scrip_ne3000pg1.nc  --dst_grd=${grid_root}/exodus_ne${NE}.g --map_file=${map_file_src_to_np4}
      ```
    1. from source to target pg2
      ```shell
      ncremap -a traave -5 --src_grd=${grid_root}/scrip_ne3000pg1.nc  --dst_grd=${grid_root}/exodus_ne${NE}.g --map_file=${map_file_src_to_pg2}
      ```
    1. from target to source (needed for calculating sub-grid anomalies on target grid)
      ```shell
      ncremap -a traave -5 --src_grd=${grid_root}/scrip_ne${NE}pg2.nc  --dst_grd=${grid_root}/scrip_ne3000pg1.nc --map_file=${map_file_pg2_to_src}
      ```

1. **Create new topograpy data on target grid**

  1. Map high-res topo to target np4 grid
    ```shell
    ncremap -m ${map_file_src_to_np4} -i ${topo_file_0} -o ${topo_file_1}
    ```

  1. Compute phi_s on the np4 grid
    ```shell
    ncap2 -s 'PHIS=terr*9.80616' ${topo_file_1} ${topo_file_2}
    ```

1. **Use homme_tool to smooth topography**

  1. Set the machine specific environment
    ```shell
    cd ${e3sm_root}/components/homme
    eval $(${e3sm_root}/cime/CIME/Tools/get_case_env)
    ```

  1. Build homme_tool
    
    The build process requires the user to select the appropriate cmake file that contains machine-specific settings.
    ```shell
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
    infilenames = '${topo_file_2}', '${topo_file_3::-3}'
    /
    EOF

    mpirun -np 8 ${e3sm_root}/components/homme/src/tool/homme_tool < input.nl
    ```

1. **Compute SGH and SGH30 on the pg2 grid, using the pg2 phi_s data**
  
  1. Remap smoothed topo to ne3000 grid
    ```shell
    ncremap -v PHIS -m ${map_file_pg2_to_src} -i ${topo_file_3} -o ${topo_file_4}
    ```

  1. Append unsmoothed ne3000 data
  ```shell
  ncks -A ${topo_file_0} ${topo_file_4}
  ```

  1. Calculate anomalies on ne3000 grid

    (Note that for ncdiff the operation is `file_3 = file_1 - file_2`)
    ```shell
    ncdiff ${topo_file_0} ${topo_file_4} ${topo_file_5}
    ```

  1. Remap anomalies back to target pg2 grid
    ```shell
    ncremap -m ${map_file_src_to_pg2} -i ${topo_file_0} -o ${topo_file_6}
    ```


1. **???? Append LANDFRAC and LANDM_COSLAT ????**

1. **Append the GLL phi_s data to the output of step 4**
  ```shell
  TOPO_FILE_2=${DATA_FILE_ROOT}/USGS-topo_ne${NE}np4_phis_x6t_tmp
  TOPO_FILE_3=${DATA_FILE_ROOT}/USGS-topo_ne${NE}np4_smoothed_x6tensor.nc
  ncks -A ${TOPO_FILE_2}1.nc ${TOPO_FILE_3}
  ```

## Batch scripts to streamline all steps

Running through all the steps above can be tedious and time-consuming. The batch script below includes all these steps as well as example Slurm batch directives for running on Perlmutter CPU nodes(NERSC). The only step that is omitted is building homme_tool, since its better to do this manually in case problems arise. Also, be sure to comment out any parts that have already been completed in advance (like creating the grid files).

To submit the slurm batch job use `sbatch batch_topo_slurm.sh`

<details open>
  <summary>batch_topo_slurm.sh</summary>
  ```shell
  #!/bin/bash
  #SBATCH --constraint=cpu
  #SBATCH --account=m3312
  #SBATCH -q regular
  #SBATCH --job-name=generate_map
  #SBATCH --output=~/E3SM/logs_slurm/slurm-%x-%j.out
  #SBATCH --time=24:00:00
  #SBATCH --nodes=1
  #SBATCH --mail-user=hannah6@llnl.gov
  #SBATCH --mail-type=END,FAIL
  
  e3sm_root=/pscratch/sd/w/whannah/e3sm_scratch/tmp_clone
  grid_root=/lcrc/group/e3sm/${USER}/scratch/chrys/files_grid
  map_root=/lcrc/group/e3sm/${USER}/scratch/chrys/files_map
  topo_root=/lcrc/group/e3sm/${USER}/scratch/chrys/files_topo
  DIN_LOC_ROOT=/lcrc/group/e3sm/data/inputdata

  NE=30

  topo_file_0=${DIN_LOC_ROOT}/atm/cam/hrtopo/USGS-topo-cube3000.nc
  topo_file_1=${topo_root}/USGS-topo_ne${NE}np4.nc
  topo_file_2=${topo_root}/USGS-topo_ne${NE}np4_phis.nc
  topo_file_3=${topo_root}/USGS-topo_ne${NE}np4_smoothed.nc
  topo_file_4=${topo_root}/USGS-topo_ne${NE}np4_smoothed_ne3000pg1.nc
  topo_file_5=${topo_root}/USGS-topo_ne${NE}np4_smoothed_ne3000pg1_anomalies.nc

  source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh

  # Create grid files for the input high res topo
  GenerateCSMesh --alt --res 3000  --file ${grid_root}/exodus_ne3000.g
  ConvertMeshToSCRIP --in ${grid_root}/exodus_ne3000.g  --out ${grid_root}/scrip_ne3000pg1.nc

  # Create grid files target EAM grid
  GenerateCSMesh --alt --res ${NE} --file ${grid_root}/exodus_ne${NE}.g
  GenerateVolumetricMesh --in ${grid_root}/exodus_ne${NE}.g --out ${grid_root}/exodus_ne${NE}pg2.g --np 2 --uniform
  ConvertMeshToSCRIP --in ${grid_root}/exodus_ne${NE}pg2.g --out ${grid_root}/scrip_ne${NE}pg2.nc

  # Create map file - source to target
  ncremap -a fv2se_flx -5 --src_grd=${grid_root}/scrip_ne3000pg1.nc  --dst_grd=${grid_root}/exodus_ne${NE}.g --map_file=${map_root}/map_ne3000pg1_to_ne${NE}np4.nc

  # Create map file - target to source
  ncremap -a traave -5 --src_grd=${grid_root}/scrip_ne${NE}pg2.nc  --dst_grd=${grid_root}/scrip_ne3000pg1.nc --map_file=${map_root}/map_ne${NE}pg2_to_ne3000pg1.nc

  # Map high-res topo to target np4 grid
  ncremap -m ${map_file_src_to_np4} -i ${topo_file_0} -o ${topo_file_1}

  # Compute phi_s on the np4 grid
  ncap2 -s 'PHIS=terr*9.80616' ${topo_file_1} ${topo_file_2}

  # Use homme_tool to smooth topography
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
  infilenames = '${topo_file_2}', '${topo_file_3::-3}'
  /
  EOF

  mpirun -np 8 ${e3sm_root}/components/homme/src/tool/homme_tool < input.nl

  # Compute SGH and SGH30 on the pg2 grid, using the pg2 phi_s data
  ncremap -v PHIS -m ${map_file_pg2_to_src} -i ${topo_file_3} -o ${topo_file_4}
  ncks -A ${topo_file_0} ${topo_file_4}
  ncdiff ${topo_file_0} ${topo_file_4} ${topo_file_5}
  ncremap -m ${map_file_src_to_pg2} -i ${topo_file_0} -o ${topo_file_6}
  

  ```
</details>
