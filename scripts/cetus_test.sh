#!/bin/csh -f

#setenv RUNJOB_ENVS "GPFS_COLLMEMPROF=1 GPFSMPIO_NAGG_PSET=16 ROMIO_HINTS=/home/pkcoff/public/romio_hints GPFSMPIO_BALANCECONTIG=1 GPFSMPIO_AGGMETHOD=2 PAMID_TYPED_ONESIDED=1 PAMID_RMA_PENDING=1M GPFSMPIO_BRIDGERINGAGG=1"
#setenv MPIEXEC_PREFLAGS " --envs GPFSMPIO_NAGG_PSET=16 "

ctest --verbose 
grep 'tests passed' ./pio2build.out | mail -s'cetus pio2 tests' jedwards@ucar.edu
