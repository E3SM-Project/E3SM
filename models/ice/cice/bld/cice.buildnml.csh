#! /bin/csh -f 

if !(-d $CASEBUILD/ciceconf) mkdir -p $CASEBUILD/ciceconf

# Invoke cice configure - output will go in $CASEBUILD/ciceconf 
cd $CASEBUILD/ciceconf || exit -1
$CODEROOT/ice/cice/bld/configure -hgrid $ICE_GRID -nx $ICE_NX -ny $ICE_NY -comp_intf $COMP_INTERFACE \
    -cice_mode $CICE_MODE -nodecomp $CICE_CONFIG_OPTS || exit -1

if ($CICE_AUTO_DECOMP == 'true') then
   @ ntasks = $NTASKS_ICE / $NINST_ICE
   set hgrid = $ICE_GRID
   if ($ICE_GRID == ar9v2) set hgrid = 'ar9v1'
   if ($ICE_GRID == ar9v4) set hgrid = 'ar9v3'
   cd $CASEBUILD
   set config = `./generate_cice_decomp.pl -res $hgrid -nx $ICE_NX -ny $ICE_NY -nproc $ntasks -thrds $NTHRDS_ICE -output all`
   cd $CASEROOT 
   if ($config[1] >= 0) then
      ./xmlchange -file env_build.xml -id CICE_BLCKX      -val $config[3] || exit -1
      ./xmlchange -file env_build.xml -id CICE_BLCKY      -val $config[4] || exit -2
      ./xmlchange -file env_build.xml -id CICE_MXBLCKS    -val $config[5] || exit -3
      ./xmlchange -file env_build.xml -id CICE_DECOMPTYPE -val $config[6] || exit -4
      ./xmlchange -file env_build.xml -id CICE_DECOMPSETTING -val $config[7] || exit -5
      source $CASEROOT/Tools/ccsm_getenv # need to do this since env_build.xml just changed
   else
      echo "ERROR configure: cice decomp not set for $ICE_GRID on $ntasks x $NTHRDS_ICE procs"
      exit -1
   endif
endif

#--------------------------------------------------------------------
# Loop over ice instances
#--------------------------------------------------------------------

cd $CASEBUILD/ciceconf || exit -1 
set default_ice_in_filename = "ice_in"

set inst_counter = 1
while ($inst_counter <= $NINST_ICE)

if ($NINST_ICE > 1) then
   set inst_string = $inst_counter
   if ($inst_counter <= 999) set inst_string = 0$inst_string
   if ($inst_counter <=  99) set inst_string = 0$inst_string
   if ($inst_counter <=   9) set inst_string = 0$inst_string
   set inst_string = "_${inst_string}"    
else
   set inst_string = ""       
endif
set ice_in_filename = ${default_ice_in_filename}${inst_string}

if (-e $CASEROOT/user_nl_cice${inst_string}) then
  ${CASEROOT}/Tools/user_nlcreate -user_nl_file $CASEROOT/user_nl_cice${inst_string} \
	-namelist_name cice_inparm >! $CASEBUILD/ciceconf/cesm_namelist 
endif

# Invoke cice build-namelist - output will go in $CASEBUILD/ciceconf
setenv INST_STRING $inst_string
$CODEROOT/ice/cice/bld/build-namelist \
    -infile $CASEBUILD/ciceconf/cesm_namelist \
    -inputdata $CASEBUILD/cice.input_data_list \
    -caseroot $CASEROOT \
    -scriptsroot $SCRIPTSROOT  \
    -inst_string "$inst_string" \
    -namelist "&cice $CICE_NAMELIST_OPTS/" -config config_cache.xml || exit -1

# Copy resolved namelist to $CASEROOT/CaseDocs and $RUNDIR
if (-d ${RUNDIR}) then
  cp $CASEBUILD/ciceconf/ice_in ${RUNDIR}/$ice_in_filename || exit -2
endif

@ inst_counter = $inst_counter + 1

end



