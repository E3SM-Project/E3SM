#!/bin/tcsh -f 
#
#  Test case to check for conservation of energy and mass
#  also does a diff on min/max u,v,p, to detect roundoff level changes 
#
if ( Linux == `uname` ) then 
    set wdir = ~/scratch1/sweqx
    set src = ~/codes/homme/build.Linux
    set input = ~/codes/homme/test/sw_conservative
    #set params = $input/Params.inc.nonetcdf
    set params = $input/Params.inc
    set NCPU = 4
endif


diff  $params  $src/../Params.inc
if ($status != 0) then
   echo "replacing Params.inc"
   cp $params $src/../Params.inc
endif



cd $src
make sweqx
mkdir $wdir
cd $wdir
mkdir movies



echo $1
if ( $1 == 'makeref' ) then
   echo making new reference output files
    mpirun -np $NCPU $src/sweqx < $input/swtc5-ne2-dt6.nl  > $input/swtc5-ne2-dt6.out
    mpirun -np $NCPU $src/sweqx < $input/swtc5-ne2.nl   > $input/swtc5-ne2.out
    mpirun -np $NCPU $src/sweqx < $input/swtc5-ne5.nl   > $input/swtc5-ne5.out
    mpirun -np $NCPU $src/sweqx < $input/swtc6.nl       > $input/swtc6.out
   exit
endif


echo "****CASE5 NE=2 dt=120 ****" running $input/swtc5-ne2.nl
mpirun -np $NCPU $src/sweqx < $input/swtc5-ne2.nl >  swtc5-ne2.new.out
#mpirun -np $NCPU $src/sweqx < $input/swtc5-ne2.nl ; exit
#$src/sweqx < $input/swtc5-ne2.nl 

tail -19 swtc5-ne2.new.out | head -18     > /tmp/run.out
tail -19 $input/swtc5-ne2.out | head -18   > /tmp/ref.out
echo "=============running diff on:============" 
cat /tmp/run.out
echo "========================================="
diff /tmp/run.out /tmp/ref.out
if ( $status == 0 ) echo "no differences" 




echo "****CASE5 NE=2 dt=6 ****" running $input/swtc5-ne2-dt6.nl
#mpirun -np $NCPU $src/sweqx < $input/swtc5-ne2-dt6.nl ; exit
mpirun -np $NCPU $src/sweqx < $input/swtc5-ne2-dt6.nl  > swtc5-ne2-dt6.new.out

tail -19 swtc5-ne2-dt6.new.out | head -18   > /tmp/run.out
tail -19 $input/swtc5-ne2-dt6.out | head -18   > /tmp/ref.out
echo "=============running diff on:============" 
cat /tmp/run.out
echo "========================================="
diff /tmp/run.out /tmp/ref.out
if ( $status == 0 ) echo "no differences" 




echo "****CASE5 NE=5 dt=120 ****" running $input/swtc5-ne5.nl
mpirun -np $NCPU $src/sweqx < $input/swtc5-ne5.nl   > swtc5-ne5.new.out

tail -19 swtc5-ne5.new.out | head -18   > /tmp/run.out
tail -19 $input/swtc5-ne5.out | head -18   > /tmp/ref.out
echo "=============running diff on:============" 
cat /tmp/run.out
echo "========================================="
diff /tmp/run.out /tmp/ref.out
if ( $status == 0 ) echo "no differences" 



echo "****CASE6 NE=5 dt=120 ****"  running $input/swtc6.nl
mpirun -np $NCPU $src/sweqx < $input/swtc6.nl   > swtc6.new.out

tail -19 swtc6.new.out | head -18   > /tmp/run.out
tail -19 $input/swtc6.out | head -18   > /tmp/ref.out
echo "=============running diff on:============" 
cat /tmp/run.out
echo "========================================="
diff /tmp/run.out /tmp/ref.out
if ( $status == 0 ) echo "no differences" 


