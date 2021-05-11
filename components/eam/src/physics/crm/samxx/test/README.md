# Testing workflow

```bash
################################################################
## Only do this first time
################################################################
cd E3SM/components/cam/src/physics/crm/samxx/test/build
./download_data.sh

################################################################
# setup cub library if not using E3SM/externals/cub
################################################################
# cd ~
# git clone https://github.com/NVlabs/cub.git
# If you put cub somewhere else, please change summit_*.sh to point to that location

################################################################
## Do this every time:
################################################################
# switch to build directory
cd E3SM/components/cam/src/physics/crm/samxx/test/build

# source the environemtn variables
source summit_gpu.sh  # or any of summit_*.sh

# activate a conda environment that includes netcdf4 and numpy
# use this command to create one:
# conda create --name crm_test_env --channel conda-forge netcdf4 numpy
source activate crm_test_env

# clean an old build
./cmakeclean.sh

# configure
./cmakescript.sh crmdata_nx32_ny1_nz28_nxrad2_nyrad1.nc crmdata_nx8_ny8_nz28_nxrad2_nyrad2.nc

# build the executables (this is also done at the top of ./runtest.sh)
make -j

# start an interactive job before running:
# (also reactivate the conda env after the interactive job starts)
# slurm:  salloc -N 1  -t 02:00:00 --account=cli115
# LSF:    bsub -Is -P cli115 -W 2:00 -nnodes 1 -q batch -J crm_standalone /bin/bash


# run the test
./runtest.sh

# if not using an interactive job submit the batch script
bsub run_standalone_batch.sh

################################################################
################################################################

# to just rerun the data comparison use a command like this
printf "\n2D data comparison:\n" ; python nccmp.py fortran2d/fortran_output_000001.nc cpp2d/cpp_output_000001.nc 
printf "\n3D data comparison:\n" ; python nccmp.py fortran3d/fortran_output_000001.nc cpp3d/cpp_output_000001.nc

```


cd E3SM/components/cam/src/physics/crm/samxx/test/build
source summit_gpu.sh  # or any of summit_*.sh
./cmakescript crmdata_nx32_ny1_nz28_nxrad2_nyrad1.nc crmdata_nx8_ny8_nz28_nxrad2_nyrad2.nc
make -j
./runtest.sh
```


