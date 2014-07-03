#! /bin/csh -f 

#------------------------------
# Verify rof grid is supported

set check_grid = "fail"
if (${ROF_GRID} == "rx1" ) set check_grid = "OK"
if (${ROF_GRID} == "null") set check_grid = "OK"

if (${check_grid} != "OK") then
  echo "ROF_GRID=${ROF_GRID} not supported in drof"
  echo "  drof support on rx1 ROF_GRID only"
  exit -2
endif

#------------------------------

if !(-d $CASEBUILD/drofconf) mkdir $CASEBUILD/drofconf
rm $CASEBUILD/drofconf/* >& /dev/null
cd $CASEBUILD/drofconf || exit -1

set default_rof_in_filename0 = "drof_in"
set default_rof_in_filename  = "drof_rof_in"

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
set rof_in_filename1 = ${default_rof_in_filename0}${inst_string}
set rof_in_filename2 = ${default_rof_in_filename}${inst_string}

if (-e $CASEROOT/user_nl_drof${inst_string}) then
  $CASEROOT/Tools/user_nlcreate -user_nl_file $CASEROOT/user_nl_drof${inst_string} \
    -namelist_name drof_inparm >! $CASEBUILD/drofconf/cesm_namelist
endif
$CODEROOT/rof/drof/bld/build-namelist \
    -infile $CASEBUILD/drofconf/cesm_namelist \
    -caseroot $CASEROOT \
    -scriptsroot $SCRIPTSROOT \
    -inst_string "$inst_string" || exit -3

if (-d ${RUNDIR}) then
   cp $CASEBUILD/drofconf/drof_in     ${RUNDIR}/$rof_in_filename1 || exit -5
   cp $CASEBUILD/drofconf/drof_rof_in ${RUNDIR}/$rof_in_filename2 || exit -7
   cp $CASEBUILD/drofconf/*.txt*      ${RUNDIR}/. >& /dev/null
endif

@ inst_counter = $inst_counter + 1

end


