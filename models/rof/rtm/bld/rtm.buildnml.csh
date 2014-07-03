#! /bin/csh -f

if !(-d $CASEBUILD/rtmconf) mkdir -p $CASEBUILD/rtmconf

#------------------------------
# Verify rof grid is supported

set check_grid = "fail"
if (${ROF_GRID} == "null")set check_grid = "OK"
if (${ROF_GRID} == "r05") set check_grid = "OK"
if (${ROF_GRID} == "r01") set check_grid = "OK"

if (${check_grid} != "OK") then
  echo "ROF_GRID=${ROF_GRID} not supported in rtm"
  echo "  rtm support on null (for single point runs), r05 and r01 ROF_GRIDs only"
  exit -2
endif

#------------------------------

set default_rof_in_filename = "rof_in"

set inst_counter = 1
while ($inst_counter <= $NINST_ROF)

if ($NINST_ROF > 1) then
   set inst_string = $inst_counter
   if ($inst_counter <= 999) set inst_string = "0$inst_string"
   if ($inst_counter <=  99) set inst_string = "0$inst_string"
   if ($inst_counter <=   9) set inst_string = "0$inst_string"
   set inst_string = "_${inst_string}"    
else
   set inst_string = ""       
endif
set rof_in_filename = ${default_rof_in_filename}${inst_string}

setenv INST_STRING $inst_string

cd $CASEBUILD/rtmconf  

if (-e $CASEBUILD/rtm.input_data_list) rm $CASEBUILD/rtm.input_data_list

# The following is for backwards compatibility when runoff restart data was on clm restart files
set finidat_rtm = ""
if ($RUN_TYPE == 'hybrid') then
  set finidat_rtm = "finidat_rtm ='${RUN_REFCASE}.rtm${inst_string}.r.${RUN_REFDATE}-${RUN_REFTOD}.nc'" 
  if ($GET_REFCASE == 'TRUE') then
    set refdir = "$DIN_LOC_ROOT/ccsm4_init/$RUN_REFCASE/$RUN_REFDATE"
    ls $refdir/*rtm* >& /dev/null
    if ( $status != 0 ) then
      set finidat_rtm = "finidat_rtm ='${RUN_REFCASE}.clm2.r.${RUN_REFDATE}-${RUN_REFTOD}.nc'" 
    endif
  endif
endif

set nrevsn_rtm = ""
if ($RUN_TYPE == 'branch') then
  set nrevsn_rtm = "${RUN_REFCASE}.rtm${inst_string}.r.${RUN_REFDATE}-${RUN_REFTOD}.nc" 
    if( -e "$RUNDIR/$nrevsn_rtm") then
      set nrevsn_rtm = "nrevsn_rtm='$nrevsn_rtm'"
    else
      set nrevsn_rtm = "${RUN_REFCASE}.rtm.r.${RUN_REFDATE}-${RUN_REFTOD}.nc" 
    endif
  if ($GET_REFCASE == 'TRUE') then
     set refdir = "$DIN_LOC_ROOT/ccsm4_init/$RUN_REFCASE/$RUN_REFDATE"
     ls $refdir/*rtm* >& /dev/null
     if ( $status != 0 ) then
       set nrevsn_rtm = "nrevsn_rtm ='${RUN_REFCASE}.clm2.r.${RUN_REFDATE}-${RUN_REFTOD}.nc'" 
     endif
  endif
endif

cat >! $CASEBUILD/rtmconf/cesm_namelist << EOF2
&rtm_inparm
 $finidat_rtm
 $nrevsn_rtm
 $RTM_NAMELIST_OPTS
EOF2
if (-e $CASEROOT/user_nl_rtm${inst_string}) then
  $UTILROOT/Tools/user_nl_add -user_nl_file $CASEROOT/user_nl_rtm${inst_string} >> $CASEBUILD/rtmconf/cesm_namelist  || exit -2
endif
cat >> $CASEBUILD/rtmconf/cesm_namelist << EOF2
/
EOF2

cd $CASEBUILD/rtmconf  
$CODEROOT/rof/rtm/bld/build-namelist \
    -infile $CASEBUILD/rtmconf/cesm_namelist \
    -caseroot $CASEROOT \
    -scriptsroot $SCRIPTSROOT \
    -inst_string "$inst_string" $RTM_BLDNML_OPTS || exit -4

if (-d ${RUNDIR}) then
  cp $CASEBUILD/rtmconf/rof_in ${RUNDIR}/$rof_in_filename || exit -2
endif

@ inst_counter = $inst_counter + 1

end


