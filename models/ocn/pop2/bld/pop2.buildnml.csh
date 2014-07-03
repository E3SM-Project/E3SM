#! /bin/csh -f

if !(-d $EXEROOT/ocn/obj   ) mkdir -p $EXEROOT/ocn/obj    || exit 2
if !(-d $EXEROOT/ocn/source) mkdir -p $EXEROOT/ocn/source || exit 3 

if ($POP_AUTO_DECOMP == 'true') then
   @ ntasks = $NTASKS_OCN / $NINST_OCN
   cd $CASEBUILD
   set config = `env UTILROOT=$UTILROOT ./generate_pop_decomp.pl -res $OCN_GRID -nproc $ntasks -thrds $NTHRDS_OCN -output all`
   cd $CASEROOT 
   if ($config[1] >= 0) then
      ./xmlchange -file env_build.xml -id POP_BLCKX      -val $config[3]
      ./xmlchange -file env_build.xml -id POP_BLCKY      -val $config[4]
      ./xmlchange -file env_build.xml -id POP_MXBLCKS    -val $config[5]
      ./xmlchange -file env_build.xml -id POP_DECOMPTYPE -val $config[6]
      ./xmlchange -file env_build.xml -id POP_NX_BLOCKS  -val $config[7]
      ./xmlchange -file env_build.xml -id POP_NY_BLOCKS  -val $config[8]
      source $CASEROOT/Tools/ccsm_getenv # need to do this since env_build.xml just changed
   else
      echo "ERROR configure: pop decomp not set for $OCN_GRID on $ntasks x $NTHRDS_OCN procs"
      exit -1
   endif
endif

if !(-d $CASEBUILD/pop2conf) mkdir $CASEBUILD/pop2conf || exit 1
cd $CASEBUILD/pop2conf || exit -1

set default_ocn_in_filename = "pop2_in"
set inst_counter = 1
while ($inst_counter <= $NINST_OCN)

if ($NINST_OCN > 1) then
    set inst_string = $inst_counter
    if ($inst_counter <= 999) set inst_string = "0$inst_string"
    if ($inst_counter <=  99) set inst_string = "0$inst_string"
    if ($inst_counter <=   9) set inst_string = "0$inst_string"
    set inst_string = _${inst_string}
else
    set inst_string = ""
endif

set ocn_in_filename = ${default_ocn_in_filename}${inst_string}

if ($NINST_OCN > 1) then
   # pop rpointer name for multi-instance case
   if (! -e $RUNDIR/rpointer.ocn${inst_string}.ovf && -e $RUNDIR/rpointer.ocn.ovf) then
      cp $RUNDIR/rpointer.ocn.ovf $RUNDIR/rpointer.ocn${inst_string}.ovf
   endif
   if (! -e $RUNDIR/rpointer.ocn${inst_string}.restart && -e $RUNDIR/rpointer.ocn.restart) then
      cp $RUNDIR/rpointer.ocn.restart $RUNDIR/rpointer.ocn${inst_string}.restart
   endif
   if (! -e $RUNDIR/rpointer.ocn${inst_string}.tavg && -e $RUNDIR/rpointer.ocn.tavg) then
      cp $RUNDIR/rpointer.ocn.tavg $RUNDIR/rpointer.ocn${inst_string}.tavg
   endif
   if (! -e $RUNDIR/rpointer.ocn${inst_string} && -e $RUNDIR/rpointer.ocn) then
      cp $RUNDIR/rpointer.ocn $RUNDIR/rpointer.ocn${inst_string}
   endif
endif

# env variable declaring type of restart file (nc vs bin) is not in any xml files
# but is needed by pop's build-namelist; it comes from rpointer.ocn.restart, which
# is in $RUNDIR for continued runs, but is in $refdir for hybrid / branch runs
# that are not continuations
setenv RESTART_INPUT_TS_FMT 'bin'
if ($RUN_TYPE == startup && $CONTINUE_RUN == 'FALSE') then
  set check_pointer_file = 'FALSE'
else
  set check_pointer_file = 'TRUE'
endif

if (($GET_REFCASE == 'TRUE') && ($RUN_TYPE != 'startup') && ($CONTINUE_RUN == 'FALSE')) then
  # During prestage step, rpointer files are copied from refdir
  set refdir = "ccsm4_init/$RUN_REFCASE/$RUN_REFDATE"
  # note: as of now, running multiple instances out of a single-instance reference
  #       case will cause problems if you try to repeat the run. A fix for this
  #       requires the prestage step to know that multi-instance is running and
  #       that reference directory only has single instance rpointer files, as well
  #       as an additional check for $DIN_LOC_ROOT/$refdir/rpointer.ocn.restart
  #       below.
  set pointer_file = "$DIN_LOC_ROOT/$refdir/rpointer.ocn${inst_string}.restart"
  if (! -e $pointer_file) then
    set pointer_file = "$RUNDIR/rpointer.ocn${inst_string}.restart"
  endif
else
  set pointer_file = "$RUNDIR/rpointer.ocn${inst_string}.restart"
endif

if ($check_pointer_file == 'TRUE') then
  grep 'RESTART_FMT=' $pointer_file >&! /dev/null
  if ($status == 0) then
    echo "Getting init_ts_file_fmt from $pointer_file"
    setenv RESTART_INPUT_TS_FMT `grep RESTART_FMT= $pointer_file | cut -c13-15`
  endif
endif

if (-e $CASEROOT/user_nl_pop2${inst_string}) then
  $UTILROOT/Tools/user_nlcreate -user_nl_file $CASEROOT/user_nl_pop2${inst_string} \
	-namelist_name pop2_inparm >! $CASEBUILD/pop2conf/cesm_namelist 
endif

# Check to see if "-preview" flag should be passed
if ( $?PREVIEW_NML ) then
  set PREVIEW_FLAG = "-preview"
else
  set PREVIEW_FLAG = ""
endif

# Check to see if build-namelist exists in SourceMods
if (-e $CASEROOT/SourceMods/src.pop2/build-namelist) then
  set BLD_NML_DIR = $CASEROOT/SourceMods/src.pop2
  set CFG_FLAG = "-cfg_dir $CODEROOT/ocn/pop2/bld"
else
  set BLD_NML_DIR = $CODEROOT/ocn/pop2/bld
  set CFG_FLAG = ""
endif

$BLD_NML_DIR/build-namelist $CFG_FLAG $PREVIEW_FLAG \
    -infile $CASEBUILD/pop2conf/cesm_namelist \
    -caseroot $CASEROOT \
    -casebuild $CASEBUILD \
    -scriptsroot $SCRIPTSROOT \
    -inst_string "$inst_string" \
    -ocn_grid "$OCN_GRID" || exit -1  

if (-d ${RUNDIR}) then
  cp $CASEBUILD/pop2conf/pop2_in ${RUNDIR}/$ocn_in_filename || exit -2
endif

# pop rpointer name for multi-instance case
foreach suffix ( ovf restart tavg ) 
   if (! -e $RUNDIR/rpointer.ocn${inst_string}.suffix && -e $RUNDIR/rpointer.ocn.suffix ) then
      cp $RUNDIR/rpointer.ocn.$suffix $RUNDIR/rpointer.ocn${inst_string}.$suffix
   endif
end 
if (! -e $RUNDIR/rpointer.ocn${inst_string} && -e $RUNDIR/rpointer.ocn) then
   cp $RUNDIR/rpointer.ocn $RUNDIR/rpointer.ocn${inst_string}
endif

if (-f $RUNDIR/pop2_in${inst_string}) rm $RUNDIR/pop2_in${inst_string}
cp -fp $CASEBUILD/pop2conf/pop2_in                    ${RUNDIR}/pop2_in${inst_string}
cp -fp $CASEBUILD/pop2conf/${OCN_GRID}_tavg_contents  ${RUNDIR}/${OCN_GRID}_tavg_contents


@ inst_counter = $inst_counter + 1

end




