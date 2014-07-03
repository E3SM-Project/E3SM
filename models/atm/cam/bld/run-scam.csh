#!/bin/csh -f

# Script to make a community-available SCAM run library for ARM standard IOP cases

# Lin Su and John Truesdale (please contact Lin Su (linsu@ucar.edu) for SCAM diagnostic package)).

set VERSION=12072012
echo "***** Version - $VERSION *****"

#########################################################################
### Set vars needed for this script code dir, case, data dir, mods, wrkdir
#########################################################################

set CAM_ROOT  =  /fs/cgd/csm/collections/cam5_3_00

set CSMDATA   = /fs/cgd/csm/inputdata

set USR_SRC =  /somedirectory/with/some/mods # Test mods dir.

###set IOP_CASES_TO_RUN = 'arm95 arm97 gateIII mpace sparticus togaII twp06'
set IOP_CASES_TO_RUN = 'arm97'

###set CAM_TIMESTEPS_TO_RUN = '60 300 600 900 1200'
set CAM_TIMESTEPS_TO_RUN = '1200'

### set CAM_LEVELS  # options are 26,27,30,60,90,120,150,180,210,240
                    # you must have initial condition files for number of levels
set CAM_LEVELS_TO_RUN = 30

set CASE = cam5_3_00

set WRKDIR    = /scratch/$USER && if (! -e  $WRKDIR) mkdir -p $WRKDIR

#########################################################################
### Select compiler+libraries env vars and set paths depending on machine.
#########################################################################

set FC_DIR=/usr/local/pgi/linux86-64
set USER_FC=pgf95
set NCHOME=/usr/local/netcdf-pgi
set DBUG = "-debug"

#########################################################################
### Shouldn't have to modify below here
#########################################################################

setenv INC_NETCDF ${NCHOME}/include
setenv LIB_NETCDF ${NCHOME}/lib
setenv NCARG_ROOT /contrib/ncarg
setenv PATH ${NCHOME}/bin:${FC_DIR}/bin:${NCARG_ROOT}/bin:${PATH}
setenv LD_LIBRARY_PATH ${FC_DIR}/lib:${LIB_NETCDF}:${LD_LIBRARY_PATH}

alias MATH 'set \!:1 = `echo "\!:3-$" | bc -l`'  # do not modify this

set runtypes="test"

#########################################################################
# NOTE: Below, set iopname, levarr, tarray.  can be more than one values, if so will loop 
#########################################################################

foreach iopname ($IOP_CASES_TO_RUN) # change this, see above (depending on case you want to simulate)

foreach tarray ($CAM_TIMESTEPS_TO_RUN)    # change this, the host model timestep 

foreach levarr ($CAM_LEVELS_TO_RUN)      # change this, number of levels to run

set EXPNAME={$CASE}_{$iopname}_L{$levarr}_T{$tarray}

#########################################################################
### Set some case specific parameters here
### Here the boundary layer cases use prescribed aerosols while the deep convection
### and mixed phase cases use prognostic aerosols.  This is because the boundary layer
### cases are so short that the aerosols do not have time to spin up.

if ($iopname == 'arm95' ||$iopname == 'arm97' ||$iopname == 'mpace' ||$iopname == 'twp06' ||$iopname == 'sparticus' ||$iopname == 'togaII' ||$iopname == 'gateIII' ||$iopname == 'IOPCASE') then
  set aero_mode = 'trop_mam3'
else
  set aero_mode = 'none'
endif

set SCAM_MODS = $WRKDIR/$CASE/mods    && if (! -e  $SCAM_MODS) mkdir -p $SCAM_MODS 
rm -rf $SCAM_MODS/*
/bin/cp -f $USR_SRC/* $SCAM_MODS

set BLDDIR    = $WRKDIR/$CASE/{$CASE}_bld_L${levarr}_${aero_mode}  && if (! -e  $BLDDIR) mkdir -p $BLDDIR
cd $BLDDIR

set IOPDESC = `grep IOP\: $CAM_ROOT/models/atm/cam/bld/namelist_files/use_cases/scam_${iopname}.xml`

echo ""
echo "***** $IOPDESC *****"
echo ""

##------------------------------------------------
## Configure for building
##------------------------------------------------
   
$CAM_ROOT/models/atm/cam/bld/configure -s -chem $aero_mode -nlev $levarr -dyn eul -res 64x128 -nospmd -nosmp -cppdefs -DDISABLE_TIMERS -scam -usr_src $SCAM_MODS -fc $USER_FC $DBUG -ldflags "-llapack -lblas -Mnobounds" -cice_nx 1 -cice_ny 1 #-microphys mg1.5

##--------------------------
## compile
##--------------------------

echo ""
echo " -- Compile"
echo ""
gmake -j >&! MAKE.out || echo "ERROR: Compile failed for' bld_${levarr}_${aero_mode} - exiting run_scam" && exit 1

#--------------------------
## Build the namelist with extra fields needed for scam diagnostics
##--------------------------

cat <<EOF >! tmp_namelistfile
&camexp 
    history_budget       = .true.,
    dtime                = $tarray,
    ndens                = 1,
    fincl1               = 'CLDST','ZMDLF','ZMDT','ZMDQ',
                          'ICWMR','ICIMR','FREQL','FREQI','LANDFRAC','CDNUMC','FICE','WSUB','CCN3','ICLDIWP',
                          'CDNUMC', 'AQSNOW',  'WSUB', 'CCN3', 'FREQI', 'FREQL', 'FREQR', 'FREQS', 'CLDLIQ', 'CLDICE', 
                          'FSDS', 'FLDS','AREL','AREI','NSNOW','QSNOW','DSNOW','AWNC','AWNI',  
                          'FLNT','FLNTC','FSNT','FSNTC','FSNS','FSNSC','FLNT','FLNTC','QRS','QRSC','QRL','QRLC',
                          'LWCF','SWCF', 'NCAI', 'NCAL', 'NIHF','NIDEP','NIIMM','NIMEY', 'ICLDIWP','ICLDTWP', 'CONCLD', 
                          'QCSEVAP', 'QISEVAP', 'QVRES', 'CMELIQ', 'CMEIOUT', 'EVAPPREC', 'EVAPSNOW', 'TAQ', 
                          'ICLMRCU', 'ICIMRCU' ,'ICWMRSH' ,'ICWMRDP', 'ICLMRTOT' , 'ICIMRTOT' , 'SH_CLD' ,  'DP_CLD', 
                          'LIQCLDF','ICECLDF', 'ICWMRST', 'ICIMRST', 'EFFLIQ', 'EFFICE','ADRAIN','ADSNOW'
/
EOF


$CAM_ROOT/models/atm/cam/bld/build-namelist -s -runtype startup -infile tmp_namelistfile -use_case scam_${iopname} -csmdata $CSMDATA \
    || echo "build-namelist failed" && exit 1

set RUNDIR    = $WRKDIR/$CASE/$EXPNAME                  && if (! -e  $RUNDIR) mkdir -p $RUNDIR
cd $RUNDIR

 ### RUN

cp -f $BLDDIR/docn.stream.txt $RUNDIR
cp -f $BLDDIR/*_in            $RUNDIR
cp -f $BLDDIR/cam             $RUNDIR

echo ""
echo " -- Running SCAM in $RUNDIR"
echo ""
./cam >&! scam_output.txt

end           #foreach iopname
end           #foreach tarray
end           #foreach levarr

exit 0
