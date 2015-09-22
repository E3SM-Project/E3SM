#! /bin/csh -f 

if !(-d $CASEBUILD/cismconf) mkdir $CASEBUILD/cismconf

cd $CASEBUILD/cismconf 
set default_cism_in_filename = "cism_in"
set default_cism_config_filename = "cism.config"
set inst_counter = 1
while ($inst_counter <= $NINST_GLC)

if ($NINST_GLC > 1) then
   set inst_string = $inst_counter
   if ($inst_counter <= 999) set inst_string = "0$inst_string"
   if ($inst_counter <=  99) set inst_string = "0$inst_string"
   if ($inst_counter <=   9) set inst_string = "0$inst_string"
   set inst_string = "_${inst_string}"    
else
   set inst_string = ""       
endif
set cism_in_filename = ${default_cism_in_filename}${inst_string}
set cism_config_filename = ${default_cism_config_filename}${inst_string}

if (-e $CASEROOT/user_nl_cism${inst_string}) then
  $UTILROOT/Tools/user_nlcreate -user_nl_file $CASEROOT/user_nl_cism${inst_string} \
    -namelist_name cism_inparm >! $CASEBUILD/cismconf/cesm_namelist || exit -1
endif
     
if (-e $CASEBUILD/cism.input_data_list) rm $CASEBUILD/cism.input_data_list

$CODEROOT/glc/cism/bld/build-namelist \
    -infile $CASEBUILD/cismconf/cesm_namelist \
    -caseroot $CASEROOT \
    -scriptsroot $SCRIPTSROOT \
    -inst_string "$inst_string" \
    -paramfile "$cism_config_filename" \
    -lnd_grid $LND_GRID -glc_grid $CISM_GRID \
    -cism_phys $CISM_PHYS || exit -2

if (-d ${RUNDIR}) then
   cp $CASEBUILD/cismconf/cism_in     ${RUNDIR}/$cism_in_filename || exit -3
   cp $CASEBUILD/cismconf/cism.config ${RUNDIR}/$cism_config_filename || exit -4
endif

@ inst_counter = $inst_counter + 1

end

if ($CISM_USE_TRILINOS == "TRUE") then
   set sourcemod_dir=$CASEROOT/SourceMods/src.cism
   if (-e ${sourcemod_dir}/trilinosOptions.xml) then
      cp ${sourcemod_dir}/trilinosOptions.xml ${RUNDIR} || exit -5
   else
      set trilinos_options_dir=$CODEROOT/glc/cism/bld/trilinosOptions
      set trilinos_file=trilinosOptions_${CISM_GRID}.xml
      if (-e ${trilinos_options_dir}/$trilinos_file) then
         cp ${trilinos_options_dir}/$trilinos_file ${RUNDIR}/trilinosOptions.xml || exit -6
      else
         echo "ERROR: no trilinosOptions file found in $trilinos_options_dir for CISM_GRID=$CISM_GRID"
	 exit -7
      endif
   endif
endif


