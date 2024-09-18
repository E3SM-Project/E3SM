# Creating an ELM initial condition file

An ELM initial condition (IC) file can be created by remapping an existing IC file from
one resolution to another using the `interpinic`, located at
`components/elm/tools/interpinic`. An ELM IC file is in the same format as an ELM restart file.
The composet of the remapped IC file will be the same as that of the input IC file.
So, for a new ELM SP-mode IC file, use an ELM input file corresponding to the SP-mode.

The steps involved in creating a new IC files are as follows:

1. Identifying an input ELM IC or restart file that will be remapped.
2. Obtaining an ELM restart file at the new resolution.
3. Compiling `interpinic` on the machine of interest.
4. Running `interpinic` to perform the interpolation.

The notes below provide an example of creating 1850 ELM IC file for the NARRM grid using E3SM v3 LR piControl from year = 0101. These notes are provided for Chrysalis.

## 1. Identification of the input ELM IC file

The identified input land condition file for this case is the following:

```bash
/lcrc/group/e3sm2/ac.golaz/E3SMv3/v3.LR.piControl/archive/rest/0101-01-01-00000/v3.LR.piControl.elm.r.0101-01-01-00000.nc
```

## 2. Obtaining an ELM restart file

Using an existing NARRM land IC and making a copy of it

```bash
cd components/elm/tools/interpinic

cp /lcrc/group/e3sm/data/inputdata/lnd/clm2/initdata_map/elmi.v3-NARRM.northamericax4v1pg2_r025_IcoswISC30E3r5.1870-01-01-00000.c20240704.nc \
elmi.v3-NARRM.northamericax4v1pg2_r025_IcoswISC30E3r5.1850-01-01-00000.c`date "+%Y%m%d"`.nc
```

## 3. Compiling `interpinic`

```bash
# Load relevant modules
cd <e3sm-dir>
eval $(./cime/CIME/Tools/get_case_env)

# change directory
cd components/elm/tools/interpinic/src

export USER_LDFLAGS="-L$NETCDF_C_DIR/lib -lnetcdf -L$NETCDF_F_DIR/lib -lnetcdff  -L$HDF5_DIR/lib -lhdf5"

USER_FC=ifort LIB_NETCDF="`nc-config --flibs`" INC_NETCDF="`nf-config --includedir`" make VERBOSE=1
```

## 4. Run `interpinic`

The `interpinic` can then be run via the following batch job (e.g., `remap.r025_RRSwISC6to18E3r4.1850.batch`) to generate the initial condition.

```bash
>cat remap.r025_RRSwISC6to18E3r4.1850.batch

#!/bin/sh
#SBATCH  --job-name=remap
#SBATCH  --nodes=1
#SBATCH  --exclusive
#SBATCH  --time 24:00:00
#SBATCH  -p slurm
#SBATCH  --account esmd

# Load relevant modules.
cd <e3sm-dir>
eval $(./cime/CIME/Tools/get_case_env/get_case_env)

# Change dir to `interpinic`
cd components/elm/tools/interpinic/src

srun -n 1 ./interpinic \
-i /lcrc/group/e3sm2/ac.golaz/E3SMv3/v3.LR.piControl/archive/rest/0101-01-01-00000/v3.LR.piControl.elm.r.0101-01-01-00000.nc \
-o elmi.v3-NARRM.northamericax4v1pg2_r025_IcoswISC30E3r5.1850-01-01-00000.c20240903.nc
```
