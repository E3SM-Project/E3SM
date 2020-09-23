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
cd E3SM/components/cam/src/physics/crm/samxx/test/build
source summit_gpu.sh  # or any of summit_*.sh
# activate an environment that includes netcdf4
# use this command to create one:
# conda create --name samxx_env --channel conda-forge netcdf4
source activate samxx_env
./cmakescript.sh crmdata_nx32_ny1_nz28_nxrad2_nyrad1.nc crmdata_nx8_ny8_nz28_nxrad2_nyrad2.nc
make -j
# if not in interactive job, start one before running:
# salloc -N 1  -t 02:00:00 --account=<account>
./runtest.sh
```

