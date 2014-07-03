#! /bin/csh -f 

if !(-d $CASEBUILD/datmconf) mkdir $CASEBUILD/datmconf
rm $CASEBUILD/datmconf/* >& /dev/null
cd $CASEBUILD/datmconf || exit -1

set default_atm_in_filename1 = "datm_in"
set default_atm_in_filename2 = "datm_atm_in"
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
set atm_in_filename1 = ${default_atm_in_filename1}${inst_string}
set atm_in_filename2 = ${default_atm_in_filename2}${inst_string}

if (-e $CASEROOT/user_nl_datm${inst_string}) then
  $CASEROOT/Tools/user_nlcreate -user_nl_file $CASEROOT/user_nl_datm${inst_string} -namelist_name datm_inparm >> $CASEBUILD/datmconf/cesm_namelist 
endif

$CODEROOT/atm/datm/bld/build-namelist \
    -infile $CASEBUILD/datmconf/cesm_namelist \
    -caseroot $CASEROOT \
    -scriptsroot $SCRIPTSROOT \
    -user_xml_dir $CASEROOT/SourceMods/src.datm \
    -inst_string "$inst_string" || exit -2

if (-d ${RUNDIR}) then
   cp $CASEBUILD/datmconf/datm_in     ${RUNDIR}/$atm_in_filename1 || exit -2
   cp $CASEBUILD/datmconf/datm_atm_in ${RUNDIR}/$atm_in_filename2 || exit -3
   cp $CASEBUILD/datmconf/*txt*       ${RUNDIR}/                  >& /dev/null
endif

@ inst_counter = $inst_counter + 1

end



