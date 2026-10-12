# Creating an ELM surface dataset

The notes describe the steps in creating an ELM surface dataset at 0.5x0.5 resolution for 1950 on Perlmutter.

## 1. Load the appropriate modules

```bash
cd <e3ms-dir>
eval $(./cime/CIME/Tools/get_case_env)
```

## 2. Compile `mksurfdata_map`

```bash
cd components/elm/tools/mksurfdata_map/src/

make clean
export USER_LDFLAGS="-L$NETCDF_DIR/lib -lnetcdf -lnetcdff -lnetcdf_intel"
export USER_LDFLAGS=$USER_LDFLAGS" -L$HDF5_DIR/lib -lhdf5 -lhdf5_fortran -lhdf5_hl_intel -lhdf5hl_fortran_intel"

USER_FC=ifort LIB_NETCDF="`nc-config --flibs`" INC_NETCDF="`nf-config --includedir`" make VERBOSE=1
```

## Build the namelist

This step assumes that the resolution for which the new surface dataset is being created is a supported resolution.
If the surface dataset is being created for an unsupported resolution, 16 mapping files will have to be created to map the raw datasets
onto this unsupported resolution. The `namelist` file with default number of glaciers (equal to zero) can be generated as:

```bash
cd ../

RES=0.5x0.5
YR=1950
DIN_LOC_ROOT=/global/cfs/cdirs/e3sm/inputdata

./mksurfdata.pl -res $RES -years $YR -d -dinlc $DIN_LOC_ROOT
```

An example of generating the namelist for 0.25 deg (`r025`) resolution for 1980 with 10 glacier layers is as follows:

```bash
RES=r025
YR=1980
DIN_LOC_ROOT=/global/cfs/cdirs/e3sm/inputdata

./mksurfdata.pl -res $RES -years $YR -d -dinlc $DIN_LOC_ROOT -glc_nec 10
```

## Run `mksurfdata_map` via an interactive job

```bash
salloc --nodes 1 --qos interactive --time 01:00:00 --constraint cpu --account e3sm

srun -n 1 ./mksurfdata_map < namelist
```

## Modifications for use on Chrysalis

The general instructions above also apply to use on Chrysalis. To build the necessary executable on Chrysalis, 
use the following instructions.

```bash
# Chrysalis modules (intel-classic + openmpi, the E3SM default stack)
source /gpfs/fs1/soft/chrysalis/spack/opt/spack/linux-centos8-x86_64/gcc-9.3.0/lmod-8.3-5be73rg/lmod/lmod/init/sh
module load perl/5.32.0-bsnc6lt intel/20.0.4-kodw73g intel-mkl/2020.4.304-g2qaxzf \
            openmpi/4.1.6-2mm63n2 hdf5/1.10.7-4cghwvq netcdf-c/4.7.4-4qjdadt \
            netcdf-fortran/4.5.3-qozrykr parallel-netcdf/1.11.0-icrpxty

# assuming you are starting from the root level of your E3SM source code
cd components/elm/tools/mksurfdata_map/src
make clean

# NetCDF / HDF5 paths derived from the loaded modules
export INC_NETCDF=$(nf-config --includedir)
export MOD_NETCDF=$INC_NETCDF
export LIB_NETCDF=$(nf-config --prefix)/lib
export USER_LDFLAGS="-L$(nc-config --prefix)/lib -lnetcdf \
  -L$(dirname $(dirname $(which h5dump)))/lib -lhdf5_hl -lhdf5"

USER_FC=ifort make VERBOSE=1

# update data input directory path for LCRC  
DIN_LOC_ROOT=/lcrc/group/e3sm/data/inputdata
```


