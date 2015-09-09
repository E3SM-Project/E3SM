#! /bin/csh -f

if (-d ${RUNDIR}) then
  cp $DIN_LOC_ROOT/wav/ww3/restart.ww3.init.seed $RUNDIR/restart.ww3
  cp $DIN_LOC_ROOT/wav/ww3/core2_G4_wns_30min_20000601_to_05.nc $RUNDIR/wind.ww3
  cp $DIN_LOC_ROOT/wav/ww3/G4L1.mod_def.ww3.121031 $RUNDIR/mod_def.ww3
endif

