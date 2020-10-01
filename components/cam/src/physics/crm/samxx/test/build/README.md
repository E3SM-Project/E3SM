# Testing workflow

```bash
################################
## Only do this first time
################################
cd E3SM/components/cam/src/physics/crm/samxx/test/build
./download_data.sh
cd ~
git clone https://github.com/NVlabs/cub.git
# If you put cub somewhere else, please change summit_*.sh to point to that location

################################
## Do this every time:
################################
cd E3SM/components/cam/src/physics/crm/samxx/test/build
source summit_gpu.sh  # or any of summit_*.sh
./cmakescript crmdata_nx32_ny1_nz28_nxrad2_nyrad1.nc crmdata_nx8_ny8_nz28_nxrad2_nyrad2.nc
make -j
./runtest.sh
```

