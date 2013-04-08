#!/bin/tcsh -f 


set wdir = ~/scratch1/sweqx
set src = ~/codes/homme/build/sweqx
set input = ~/codes/homme/test/sw_conservative


set NCPU = 1
if ( ${?PBS_NODEFILE} ) then
   set NCPU = `wc $PBS_NODEFILE | awk '{print $1}' - `
endif
echo using NCPU = $NCPU


#configure the model
cd $src
set configure = 0
if ( $configure ) then
  cd $src
  ./configure --with-netcdf=$NETCDF_PATH --with-pnetcdf=$PNETCDF_PATH NP=8 PLEV=1
  if ($status ) exit

  make clean
  make -j4 depends
  if ($status ) exit
endif

make -j2 sweqx
if ($status ) exit



cd $src
make -j3 sweqx
mkdir $wdir
cd $wdir
mkdir movies


echo running $input/swtc5-ne5-si.nl
mpirun -np $NCPU $src/sweqx < $input/swtc5-ne5-si.nl  
#mpirun -np $NCPU $src/sweqx < ~/codes/homme/test/swtc5/semi_imp.nl

