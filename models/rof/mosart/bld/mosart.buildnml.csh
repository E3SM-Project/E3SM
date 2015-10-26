#! /bin/csh -f

if !(-d $CASEBUILD/mosartconf) mkdir -p $CASEBUILD/mosartconf

#------------------------------
# Verify rof grid is supported

set check_grid = "fail"
if (${ROF_GRID} == "null")set check_grid = "OK"
if (${ROF_GRID} == "r05") set check_grid = "OK"
if (${ROF_GRID} == "r01") set check_grid = "OK"
if (${ROF_GRID} == "NLDAS") set check_grid = "OK"

if (${check_grid} != "OK") then
  echo "ROF_GRID=${ROF_GRID} not supported in mosart"
  echo "  mosart support on null (for single point runs), r05 and r01 ROF_GRIDs only"
  exit -2
endif

#------------------------------

set default_mosart_in_filename = "mosart_in"

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
set mosart_in_filename = ${default_mosart_in_filename}${inst_string}

setenv INST_STRING $inst_string

cd $CASEBUILD/mosartconf  

if (-e $CASEBUILD/mosart.input_data_list) rm $CASEBUILD/mosart.input_data_list

# The following is for backwards compatibility when runoff restart data was on clm restart files
set finidat_rtm = ""
set nrevsn_rtm = ""
if (${ROF_GRID} != "null") then
if ($RUN_TYPE == 'hybrid' || $RUN_TYPE == "branch" ) then

  # set search directory
  if ($GET_REFCASE == 'TRUE') then
    set refdir = "$DIN_LOC_ROOT/ccsm4_init/$RUN_REFCASE/$RUN_REFDATE"
  else
    set refdir = "$RUNDIR"
  endif

  # search for clm or mosart files with instance or not
  set fncheck = "${RUN_REFCASE}.mosart${inst_string}.r.${RUN_REFDATE}-${RUN_REFTOD}.nc"
  if !(-e "$refdir/$fncheck") then
    set fncheck = "${RUN_REFCASE}.mosart.r.${RUN_REFDATE}-${RUN_REFTOD}.nc"
    if !(-e "$refdir/$fncheck") then
      set fncheck = "${RUN_REFCASE}.clm2${inst_string}.r.${RUN_REFDATE}-${RUN_REFTOD}.nc"
      if !(-e "$refdir/$fncheck") then
        set fncheck = "${RUN_REFCASE}.clm2.r.${RUN_REFDATE}-${RUN_REFTOD}.nc"
        if !(-e "$refdir/$fncheck") then
          echo "mosart.buildnml.csh could not find restart file for branch or hybrid start"
          echo "refdir is" $refdir
          exit -8
        endif
      endif
    endif
  endif

  # set the namelist variable needed
  if ($RUN_TYPE == "hybrid") then
    set finidat_rtm = "finidat_rtm = '$fncheck'"
  endif
  if ($RUN_TYPE == "branch") then
    set nrevsn_rtm = "nrevsn_rtm = '$refdir/$fncheck'"
  endif

endif
endif

cat >! $CASEBUILD/mosartconf/cesm_namelist << EOF2
&mosart_inparm
 $finidat_rtm
 $nrevsn_rtm
 $MOSART_NAMELIST_OPTS
EOF2
if (-e $CASEROOT/user_nl_mosart${inst_string}) then
  $UTILROOT/Tools/user_nl_add -user_nl_file $CASEROOT/user_nl_mosart${inst_string} >> $CASEBUILD/mosartconf/cesm_namelist  || exit -2
endif
cat >> $CASEBUILD/mosartconf/cesm_namelist << EOF2
/
EOF2

cd $CASEBUILD/mosartconf  
$CODEROOT/rof/mosart/bld/build-namelist \
    -infile $CASEBUILD/mosartconf/cesm_namelist \
    -caseroot $CASEROOT \
    -scriptsroot $SCRIPTSROOT \
    -inst_string "$inst_string" $MOSART_BLDNML_OPTS || exit -4

if (-d ${RUNDIR}) then
  cp $CASEBUILD/mosartconf/mosart_in ${RUNDIR}/$mosart_in_filename || exit -2
endif

@ inst_counter = $inst_counter + 1

end


