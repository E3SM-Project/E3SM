#! /bin/csh -f 

if !(-d $CASEBUILD/dlndconf) mkdir $CASEBUILD/dlndconf
rm $CASEBUILD/dlndconf/* >& /dev/null
cd $CASEBUILD/dlndconf || exit -1

set default_lnd_in_filename0 = "dlnd_in"
set default_lnd_in_filename  = "dlnd_lnd_in"
set inst_counter = 1
while ($inst_counter <= $NINST_LND)

if ($NINST_LND > 1) then
   set inst_string = $inst_counter
   if ($inst_counter <= 999) set inst_string = "0$inst_string"
   if ($inst_counter <=  99) set inst_string = "0$inst_string"
   if ($inst_counter <=   9) set inst_string = "0$inst_string"
   set inst_string = "_${inst_string}"    
else
   set inst_string = ""       
endif
set lnd_in_filename1 = ${default_lnd_in_filename0}${inst_string}
set lnd_in_filename2 = ${default_lnd_in_filename}${inst_string}

if (-e $CASEROOT/user_nl_dlnd${inst_string}) then
  $CASEROOT/Tools/user_nlcreate -user_nl_file $CASEROOT/user_nl_dlnd${inst_string} \
    -namelist_name dlnd_inparm >! $CASEBUILD/dlndconf/cesm_namelist
endif
$CODEROOT/lnd/dlnd/bld/build-namelist \
    -infile $CASEBUILD/dlndconf/cesm_namelist \
    -caseroot $CASEROOT \
    -scriptsroot $SCRIPTSROOT \
    -inst_string "$inst_string" || exit -2

if (-d ${RUNDIR}) then
   cp $CASEBUILD/dlndconf/dlnd_in     ${RUNDIR}/$lnd_in_filename1 || exit -5
   cp $CASEBUILD/dlndconf/dlnd_lnd_in ${RUNDIR}/$lnd_in_filename2 || exit -6
   cp $CASEBUILD/dlndconf/*txt*       ${RUNDIR}/. >& /dev/null
endif

@ inst_counter = $inst_counter + 1

end


