#! /bin/csh -f 

if !(-d $CASEBUILD/camconf) mkdir $CASEBUILD/camconf

cd $CASEBUILD/camconf || exit -1

if ($BUILD_COMPLETE == 'FALSE') then

  #--------------------------------------------------------------------
  # Invoke cam configure - output will go in CASEBUILD/camconf
  #--------------------------------------------------------------------

  # Some settings for single column mode.
  if ($PTS_MODE == 'TRUE') then
     set scm  = "-scam -nosmp"
  else
     set scm  = ""
  endif

  if ($MPILIB == 'mpi-serial') then
     set spmd = "-nospmd"
  else
     set spmd = "-spmd"
  endif

  # The ocean component setting is only used by CAM to do attribute matching for
  # setting default tuning parameter values.  In SOM mode we want to use the same
  # tunings as the fully coupled B compset, so set the ocean component to pop2 in
  # that case.
  set ocn = $COMP_OCN
  if ($?DOCN_MODE) then
    if ($DOCN_MODE == 'som') set ocn = pop2
  endif

  if ($COMP_INTERFACE == 'MCT' ) set comp = mct
  if ($COMP_INTERFACE == 'ESMF') set comp = esmf

  $CODEROOT/atm/cam/bld/configure -s  -ccsm_seq -ice none -ocn $ocn  -comp_intf $comp \
      $scm $spmd -dyn $CAM_DYCORE -res $ATM_GRID  $CAM_CONFIG_OPTS || exit -1

else

  # Verify that we have a config_cache file.
  if !(-e $CASEBUILD/camconf/config_cache.xml) then
    echo "cam.buildnml.csh: Build is complete but config_cache.xml is missing."
    echo "Cannot run build-namelist; try cleaning build and building again."
    exit -1
  endif

endif

#--------------------------------------------------------------------
# Invoke cam build-namelist - output will go in $CASEBUILD/camconf
#--------------------------------------------------------------------

if ($RUN_STARTDATE =~ *-01-01* || $RUN_STARTDATE =~ *-09-01*) then
    set ignore = "-ignore_ic_year"
else
    set ignore = "-ignore_ic_date"
endif
if ($CAM_NML_USE_CASE == UNSET) then
    set usecase = " "
else
    set usecase = "-use_case $CAM_NML_USE_CASE"
endif
if ($PTS_MODE == 'TRUE') then 
  # setting for single column mode.
  set scmb = "scmlon=$PTS_LON scmlat=$PTS_LAT"
else
  set scmb = ""
endif

set default_atm_in_filename = "atm_in"
set inst_counter = 1
while ($inst_counter <= $NINST_ATM)

if ($NINST_ATM > 1) then
   set inst_string = $inst_counter
   if ($inst_counter <= 999) set inst_string = "0$inst_string"
   if ($inst_counter <=  99) set inst_string = "0$inst_string"
   if ($inst_counter <=   9) set inst_string = "0$inst_string"
   set inst_string = "_${inst_string}"    
else
   set inst_string = ""       
endif
set atm_in_filename = ${default_atm_in_filename}${inst_string}

# create camconf/cesm_namelist
set ncdata = ""
if ($RUN_TYPE == 'hybrid') then
  if (-e ${RUN_REFCASE}.cam${inst_string}.i.${RUN_REFDATE}-${RUN_REFTOD}.nc) then
    set ncdata = "ncdata='${RUN_REFCASE}.cam${inst_string}.i.${RUN_REFDATE}-${RUN_REFTOD}.nc'"
  else 
    set ncdata = "ncdata='${RUN_REFCASE}.cam.i.${RUN_REFDATE}-${RUN_REFTOD}.nc'"
  endif
endif
set cam_branch_file = ""
if ($RUN_TYPE == 'branch') then
  if (-e ${RUNDIR}/${RUN_REFCASE}.cam${inst_string}.r.${RUN_REFDATE}-${RUN_REFTOD}.nc) then    
    set cam_branch_file = "cam_branch_file='${RUNDIR}/${RUN_REFCASE}.cam${inst_string}.r.${RUN_REFDATE}-${RUN_REFTOD}.nc'"
  else
    set cam_branch_file = "cam_branch_file='${RUNDIR}/${RUN_REFCASE}.cam.r.${RUN_REFDATE}-${RUN_REFTOD}.nc'"
  endif
endif
set co2vmr = "co2vmr=${CCSM_CO2_PPMV}e-6"
set yyyymmdd = `echo $RUN_STARTDATE | sed s/-//g `
set start_ymd = "start_ymd = $yyyymmdd"
@ dtime  = ( 3600 * 24 ) / $ATM_NCPL
@ ntasks = $NTASKS_ATM / $NINST_ATM
set scmb = 

cat >! $CASEBUILD/camconf/cesm_namelist << EOF2
&cam_inparm
 dtime = $dtime
 $ncdata
 $cam_branch_file
 $co2vmr
 $start_ymd 
 $scmb
EOF2
if (-e $CASEROOT/user_nl_cam${inst_string}) then
  $UTILROOT/Tools/user_nl_add -user_nl_file $CASEROOT/user_nl_cam${inst_string} >> $CASEBUILD/camconf/cesm_namelist 
endif
cat >> $CASEBUILD/camconf/cesm_namelist << EOF2
/
EOF2
     
if (-e $CASEBUILD/cam.input_data_list) rm $CASEBUILD/cam.input_data_list

$CODEROOT/atm/cam/bld/build-namelist -infile $CASEBUILD/camconf/cesm_namelist \
    -csmdata $DIN_LOC_ROOT $ignore $usecase -inputdata $CASEBUILD/cam.input_data_list \
    -ntasks $ntasks -namelist "&atmexp $CAM_NAMELIST_OPTS /" || exit -1

if (-d ${RUNDIR}) then
   cp $CASEBUILD/camconf/atm_in ${RUNDIR}/$atm_in_filename || exit -2
   cp $CASEBUILD/camconf/drv_flds_in ${RUNDIR}/drv_flds_in || exit -2
endif

@ inst_counter = $inst_counter + 1

end


