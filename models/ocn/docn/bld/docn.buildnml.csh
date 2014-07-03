#! /bin/csh -f 

if !(-d $CASEBUILD/docnconf) mkdir $CASEBUILD/docnconf
rm $CASEBUILD/docnconf/* >& /dev/null
cd $CASEBUILD/docnconf || exit -1

set default_ocn_in_filename1 = "docn_in"
set default_ocn_in_filename2 = "docn_ocn_in"
set inst_counter = 1

while ($inst_counter <= $NINST_OCN)

if ($NINST_OCN > 1) then
   set inst_string = $inst_counter
   if ($inst_counter <= 999) set inst_string = "0$inst_string"
   if ($inst_counter <=  99) set inst_string = "0$inst_string"
   if ($inst_counter <=   9) set inst_string = "0$inst_string"
   set inst_string = "_${inst_string}"    
else
   set inst_string = ""       
endif
set ocn_in_filename1 = ${default_ocn_in_filename1}${inst_string}
set ocn_in_filename2 = ${default_ocn_in_filename2}${inst_string}

if (-e $CASEROOT/user_nl_docn${inst_string}) then
  $CASEROOT/Tools/user_nlcreate -user_nl_file $CASEROOT/user_nl_docn${inst_string} -namelist_name docn_inparm >! $CASEBUILD/docnconf/cesm_namelist 
endif
$CODEROOT/ocn/docn/bld/build-namelist \
    -infile $CASEBUILD/docnconf/cesm_namelist \
    -caseroot $CASEROOT \
    -scriptsroot $SCRIPTSROOT \
    -inst_string "$inst_string" || exit -2
    
if (-d ${RUNDIR}) then
   cp $CASEBUILD/docnconf/docn_in     ${RUNDIR}/$ocn_in_filename1 || exit -2
   cp $CASEBUILD/docnconf/docn_ocn_in ${RUNDIR}/$ocn_in_filename2 || exit -2
   foreach file (*txt*)
      cp $CASEBUILD/docnconf/$file   ${RUNDIR}/$file || exit -2
   end
endif

@ inst_counter = $inst_counter + 1

end



