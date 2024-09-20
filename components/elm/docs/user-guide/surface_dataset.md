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
