#!/bin/tcsh -f 


set wdir = ~/scratch1/sweqx
set src = ~/codes/homme/build.Linux
set input = ~/codes/homme/test/sw_conservative


set NCPU = 1
if ( ${?PBS_NODEFILE} ) then
   set NCPU = `wc $PBS_NODEFILE | awk '{print $1}' - `
endif
echo using NCPU = $NCPU

set params = $input/Params.inc
diff  $params  $src/../Params.inc
if ($status != 0) then
   echo "replacing Params.inc"
   cp $params $src/../Params.inc
endif


cd $src
make -j3 sweqx
mkdir $wdir
cd $wdir
mkdir movies


echo running $input/swtc5-ne5-si.nl
mpirun -np $NCPU $src/sweqx < $input/swtc5-ne5-si.nl  
#mpirun -np $NCPU $src/sweqx < ~/codes/homme/test/swtc5/semi_imp.nl

