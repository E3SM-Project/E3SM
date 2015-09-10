#!/bin/bash -f

set -e

#CHANGE PATHS TO YOUR LOCALS
# Oksana:
HOMME=$HOME/homme7
MACH=$HOMME/cmake/machineFiles/climate.cmake  # CMAKE
wdir=$HOME/runhomme/sweqx
NCPU=46

# Mark:
HOMME=$HOME/codes/homme
MACH=$HOMME/cmake/machineFiles/redsky.cmake  # CMAKE
wdir=$HOME/scratch1/sweqx
NCPU=0


input=$HOMME/test/sw_conservative



test_case=swtc5
#nu=3.5e-8 #3.5e-8 for hv_scaling=3.2, 1.2e-6 for hv_scaling=4.0
nu=40e-8 # maches nu=1e15 for ne30
hvscaling=3.2 

# arm x8 grid:
#   48e-8 blows up.   theory:  unstable at 61e-8
#   40e-8 stable
#   8e-8 matches ne30

#tstep=50
#meshfile="$HOMME/test/mesh_refine/grids/grid_10_x8_iter10halo2.g"
#meshfile="'$HOMME/test/mesh_refine/grids/mountain_10_x8.g'"

#tstep=90
#meshfile="'none'  ne=30"

tstep=10
meshfile="'$HOMME/utils/CUBIT_scripts/exodus/arm_30_x8_lowconn.g'"



name=${test_case}

# old configure build system:
#cd $HOMME/build/sweqx
#./configure PLEV=1 NP=4 --with-netcdf=$NETCDF_PATH  --with-pnetcdf=$PNETCDF_PATH --enable-blas --enable-lapack
#	make  -j 4 depends
#	make clean
#rm -f sweqx
#make -j 4
#exe=$HOMME/build/sweqx/sweqx

# CMAKE 

cd $wdir
#rm -rf CMakeFiles CMakeCache.txt
#cmake -C $MACH -DSWEQX_PLEV=1  -DSWEQX_NP=4 $HOMME
make -j4 sweqx
exe=$wdir/src/sweqx/sweqx


cd $wdir
let sfreq=24*3600
sfreq=`echo "$sfreq / $tstep" | bc`

sed s/tstep.\*/"tstep = $tstep"/  $input/swtc5-tensor-hv.nl |\
    sed s/hypervis_scaling.\*/"hypervis_scaling = $hvscaling"/  |\
    sed s/nu=.\*/"nu= $nu"/  |\
    sed s/nu_s=.\*/"nu_s= $nu"/  |\
    sed s:mesh_file.\*:"mesh_file=$meshfile ":  |\
    sed s/statefreq.\*/"statefreq=$sfreq"/  \
    > input.nl

pwd
set echo
time mpirun -np $NCPU  $exe < input.nl | tee  sweq.out

mv -f sweq.mass $name.mass
mv -f sweq.out $name.out
mv -f movies/swtc51.nc movies/$name.nc


exit
