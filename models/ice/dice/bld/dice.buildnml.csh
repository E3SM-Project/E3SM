#! /bin/csh -f 

if !(-d $CASEBUILD/diceconf) mkdir $CASEBUILD/diceconf
rm $CASEBUILD/diceconf/* >& /dev/null
cd $CASEBUILD/diceconf || exit -1

set default_ice_in_filename1 = "dice_in"
set default_ice_in_filename2 = "dice_ice_in"
set inst_counter = 1
while ($inst_counter <= $NINST_ICE)

if ($NINST_ICE > 1) then
   set inst_string = $inst_counter
   if ($inst_counter <= 999) set inst_string = "0$inst_string"
   if ($inst_counter <=  99) set inst_string = "0$inst_string"
   if ($inst_counter <=   9) set inst_string = "0$inst_string"
   set inst_string = "_${inst_string}"    
else
   set inst_string = ""       
endif
set ice_in_filename1 = ${default_ice_in_filename1}${inst_string}
set ice_in_filename2 = ${default_ice_in_filename2}${inst_string}

if (-e $CASEROOT/user_nl_dice${inst_string}) then
  $CASEROOT/Tools/user_nlcreate -user_nl_file $CASEROOT/user_nl_dice${inst_string} -namelist_name dice_inparm >! $CASEBUILD/diceconf/cesm_namelist 
endif

$CODEROOT/ice/dice/bld/build-namelist \
    -infile $CASEBUILD/diceconf/cesm_namelist \
    -caseroot $CASEROOT \
    -scriptsroot $SCRIPTSROOT \
    -inst_string "$inst_string" || exit -2
    
if (-d ${RUNDIR}) then
   cp $CASEBUILD/diceconf/dice_in     ${RUNDIR}/$ice_in_filename1 || exit -3
   cp $CASEBUILD/diceconf/dice_ice_in ${RUNDIR}/$ice_in_filename2 || exit -4
   cp $CASEBUILD/diceconf/*txt*       ${RUNDIR}/                  >& /dev/null
endif

@ inst_counter = $inst_counter + 1

end



