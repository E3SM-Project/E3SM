#! /bin/csh -fe

### This script was created 2015-11-15 by Philip Cameron-Smith (pjc@llnl.gov) and Peter Caldwell
### and incorporates some features originally from Hui Wan, Kai Zhang, and Balwinder Singh.
### Significant improvements from Michael Deakin and Chris Golaz.
###
set first_argument = $1
if ( $first_argument != '' ) then
 echo 'This script does everything needed to configure/compile/run a case. As such, it'
 echo 'provides complete provenance for each run and makes sharing with newbies easy. Future'
 echo 'users should make sure that everything required for a run is in this script, the ACME'
 echo 'git repo, or the inputdata svn repo.'
 echo '** This script accepts no arguments. Please edit the script as needed and resubmit without arguments. **'
 exit 5
endif

# NOTE: CIME 5 and git commands are not cwd agnostic, so compute the absolute paths, then cd to the directories as needed
set this_script_name = `basename $0`
set relative_dir = `dirname $0`
set case_scripts_dir = `pwd -P`/$case_name
set this_script_dir = `cd $relative_dir ; pwd -P`
set this_script_path = $this_script_dir/$this_script_name

###===================================================================
### THINGS USERS USUALLY CHANGE (SEE END OF SECTION FOR GUIDANCE)
###===================================================================

### BASIC INFO ABOUT RUN
set run_name       = run_acme.template
set job_name       = A_WCYCL1850_template
set compset        = A_WCYCL1850
set resolution     = ne4np4_oQU240
set machine        = default
set walltime       = 00:10
setenv project       fy150001

### SOURCE CODE OPTIONS
set fetch_code     = false
set acme_tag       = master
set tag_name       = ACME

### CUSTOM CASE_NAME
set case_name = ${tag_name}.${job_name}.${resolution}

### BUILD OPTIONS
set debug_compile  = false

### AUTOMATIC DELETION OPTIONS
set seconds_before_delete_source_dir = -1
set seconds_before_delete_case_dir   = 10
set seconds_before_delete_run_dir    = -1

### SUBMIT OPTIONS
set submit_run       = true
set debug_queue      = false

### PROCESSOR CONFIGURATION
set processor_config = S

### STARTUP TYPE
set model_start_type = initial
set restart_files_dir = none

### DIRECTORIES
set code_root_dir    = `cd $this_script_dir/..; pwd -P`

### LENGTH OF SIMULATION, RESTARTS, AND ARCHIVING
set stop_units       = ndays
set stop_num         = 5
set restart_units    = $stop_units
set restart_num      = $stop_num
set num_submits      = 1
set do_short_term_archiving      = false
set do_long_term_archiving       = false

### SIMULATION OPTIONS
set atm_output_freq              = -24
set records_per_atm_output_file  = 40

### COUPLER HISTORY FILES
set do_cpl_hist    = true
set cpl_hist_units = ndays
set cpl_hist_num   = 1

#==============================
#EXPLANATION FOR OPTIONS ABOVE:
#==============================

### BASIC INFO ABOUT RUN (1)

#run_name: the run will be named: ${tag_name}.${compset}.${resolution}.${run_name}.  run_name is to explain the
#    purpose of the run (e.g. run_name=ParallelPhysDyn) or just to ensure the run name is unique (e.g. run_name=test1).
#job_name: This is only used to name the job in the batch system. The problem is that batch systems often only
#    provide the first few letters of the job name when reporting on jobs inthe queue, which may not be enough
#    to distinguish simultaneous jobs.
#compset: indicates which model components and forcings to use. List choices by typing `create_newcase -list compsets`.
#    An (outdated?) list of options is available at http://www.cesm.ucar.edu/models/cesm1.0/cesm/cesm_doc_1_0_4/a3170.html
#resolution: Model resolution to use. Type `create_newcase -list grids` for a list of options or see
#    http://www.cesm.ucar.edu/models/cesm1.0/cesm/cesm_doc_1_0_4/a3714.htm. Note that ACME always uses ne30 or ne120 in the atmosphere.
#machine: what machine you are going to run on. This should be 'default' for most machines, and only changed for machines with multiple queues, such as cori.
#    Valid machine names can also be listed via `create_newcase -list machines`
#walltime: How long to reserve the nodes for. The format is HH:MM(:SS); ie 00:10 -> 10 minutes.
#    Setting this to 'default' has the script determine a reasonable value for most runs.
#project: what bank to charge for your run time. May not be needed on some machines.
#    NOTE: project must be an *environment* variable on some systems.

### SOURCE CODE OPTIONS (2)

#fetch_code: If True, downloads code from github. If False, code is assumed to exist already.
#    NOTE: it is assumed that you have access to the ACME git repository.  To get access, see:
#    https://acme-climate.atlassian.net/wiki/display/Docs/Installing+the+ACME+Model
#acme_tag: ACME tagname in github. Can be 'master', a git hash, a tag, or a branch name. Only used if fetch_code=True.
#    NOTE: If acme_tag=master or master_detached, then this script will provide the latest master version, but detach from the head,
#          to minimize the risk of a user pushing something to master.
#tag_name: Short name for the ACME branch used. If fetch_code=True, this is a shorter replacement for acme_tag
#    (which could be a very long hash!). Otherwise, this is a short name for the branch used. You can
#    choose TAG_NAME to be whatever you want.

### BUILD OPTIONS (3)

#debug_compile: If TRUE, then compile with debug flags
#    Compiling in debug mode will stop the run at the actual location an error occurs, and provide more helpful output.
#    However, it runs about 10 times slower, and is not bit-for-bit the same because some optimizations make tiny change to the
#    numerics.

### AUTOMATIC DELETION OPTIONS (4)

#seconds_before_delete_source_dir : If seconds_before_delete_source_dir>=0 and fetch_code=true, this script automatically deletes
#    the old source code directory after waiting seconds_before_delete_source_dir seconds (to give you the opportunity to cancel
#    by pressing ctrl-C). To turn off the deletion (default behavior), set $num_seconds_before_deleting_source_dir to be negative.
#seconds_before_delete_case_dir : Similar to above, but remove the old case_scripts directory. Since create_newcase dies whenever
#    the case_scripts directory exists, it only makes sense to use $seconds_before_delete_case_dir<0 if you want to be extra careful and
#    only delete the case_scripts directory manually.
#seconds_before_delete_run_dir : As above, but the old run directory will be deleted.  This makes for a clean start.

### SUBMIT OPTIONS (5)

#submit_run:  If True, then submit the batch job to start the simulation.
#debug_queue: If True, then use the debug queue, otherwise use the queue specified in the section on QUEUE OPTIONS.

### PROCESSOR CONFIGURATION (6)

#processor_config: Indicates what processor configuration to use.
#    1=single processor, S=small, M=medium, L=large, X1=very large, X2=very very large, CUSTOM=defined below.
#    The actual configurations for S,M,L,X1,X2 are dependent on the machine.

### STARTUP TYPE (7)

#model_start_type:  Specify how this script should initialize the model:  initial, continue, branch.
#    These options are not necessarily related to the CESM run_type options.
#    'initial' means the intial files will be copied into the run directory,
#    and the ACME run_type can be 'initial', 'hybrid', or 'restart', as specified by this script below.
#    'continue' will do a standard restart, and assumes the restart files are already in the run directory.
#    'branch' is almost the same, but will set RUN_TYPE='branch', and other options as specified by this script below.
#    NOTE: To continue an existing simulation, it may be easier to edit env_run and [case].run manually in the
#    case_scripts directory.  The biggest difference is that doing it with this script
#    may delete the previous case_scripts directory, and will provide a way to pass a simulation to someone else.

### DIRECTORIES (8)

#code_root_dir: The directory that contains (or will contain) the source code and other code files. (formerly $CCSMROOT)
#     If fetch_code=false, this is the location where the code already resides.
#     If fetch_code=true, this is where to put the code.

### LENGTH OF SIMULATION, RESTARTS, AND ARCHIVING (9)

#stop_units: The units for the length of run, eg nhours, ndays, nmonths, nyears.
#stop_num: The simulation length for each batch submission, in units of $stop_units.
#restart_units: The units for how often restart files are written, eg nhours, ndays, nmonths, nyears.
#restart_num: How often restart files are written, in units of $restart_units.
#num_submits: After a batch job finishes successfully, a new batch job will automatically be submitted to
#    continue the simulation.  num_submits is the total number of submissions.
#    After the first submission, the CONTINUE_RUN flage in env_run.xml will be changed to TRUE.
#do_short_term_archiving: If TRUE, then move simulation output to the archive directory in your scratch directory.
#do_long_term_archiving : If TRUE, then move simulation output from the short_term_archive to the local mass storage system.

### SIMULATION OPTIONS (10)

#atm_output_freq (the namelist variable is nhtfrq) : The frequency with which the atmosphere writes its output.
#    0=monthly, +N=every N timesteps,  -N=every N hours
#    For more information:   http://www.cesm.ucar.edu/models/atm-cam/docs/usersguide/node45.html
#records_per_atm_output_file (the namelist variable is mfilt):  The number of time records in each netCDF output file
#    from the atmosphere model. If atm_output_freq=0, then there will only be one time record per file.
#NOTE: If there will be more than one 'history tape' defined in the atm namelist, then
#    atm_output_freq and records_per_atm_output_file should be a comma-separated list of numbers
#    (specifying the option for each history tape).

### GENERAL NOTES:

# 1. capitalization doesn't matter on most of the variables above because we lowercase variables before using them.
# 2. most of the code below does things you probably never want to change. However, users will often make settings
#    in the "USER_NL" and "RUN CONFIGURATION OPTIONS" sections below.

### PROGRAMMING GUIDELINES
#
# +) The exit error numbers are sequential through the code:
#        0-099 are before create_newcase
#      100-199 are between create_newcase and cesm_setup
#      200-299 are between cesm_setup and case_scripts.build
#      300-399 are between case_scripts.build and case_scripts.submit
#      400-499 are after case_scripts.submit
#    If this script dies, then print out the exit code.
#    (in csh: use 'echo $status' immediatedly after the script exits)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#  END OF COMMON OPTIONS - you may need to change things below here to access advanced
#  capabilities, but if you do you should know what you're doing.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#===========================================
# DOCUMENT WHICH VERSION OF THIS SCRIPT IS BEING USED:
#===========================================
set script_ver = 3.0.2

echo ''
echo 'run_acme: ++++++++ run_acme starting ('`date`'), version '$script_ver' ++++++++'
echo ''

#===========================================
# DEFINE THINGS NEEDED LATER:
#===========================================

alias lowercase "echo \!:1 | tr '[A-Z]' '[a-z]'"  #make function which lowercases any strings passed to it.
alias uppercase "echo \!:1 | tr '[a-z]' '[A-Z]'"  #make function which uppercases any strings passed to it.

#===========================================
# DOWNLOAD SOURCE CODE IF NEEDED:
#===========================================

### NOTE: you must be setup with access to the ACME repository before you can clone the repository. For access, see
###       https://acme-climate.atlassian.net/wiki/display/Docs/Installing+the+ACME+Model

if ( `lowercase $fetch_code` == true ) then
  echo 'run_acme: Downloading code from the ACME git repository.'
  if ( -d $code_root_dir/$tag_name ) then
    if ( $seconds_before_delete_source_dir >= 0 ) then
      set num_seconds_until_delete = $seconds_before_delete_source_dir
      echo 'run_acme: Removing old code directory '$code_root_dir/$tag_name' in '$num_seconds_until_delete' seconds.'
      echo 'run_acme: To abort, press ctrl-C'
      while ( ${num_seconds_until_delete} > 0 )
         echo 'run_acme:  '${num_seconds_until_delete}'  seconds until deletion.'
         sleep 1
         @ num_seconds_until_delete = ${num_seconds_until_delete} - 1
      end
      rm -fr $code_root_dir/$tag_name
      echo 'run_acme: Deleted '$code_root_dir/$tag_name
    else
      echo 'run_acme: ERROR: Your branch tag already exists, so dying instead of overwriting.'
      echo '          You likely want to either set fetch_code=false, change $tag_name, or'
      echo '          change seconds_before_delete_source_dir.'
      echo '          Note: $fetch_code = '$fetch_code
      echo '                $code_root_dir/$tag_name = '$code_root_dir/$tag_name
      echo '                $seconds_before_delete_source_dir = '$seconds_before_delete_source_dir
      exit 20
    endif #$seconds_before_delete_source_dir >=0
  endif #$code_root_dir exists

  echo 'run_acme: Cloning repository into $tag_name = '$tag_name'  under $code_root_dir = '$code_root_dir
  mkdir -p $code_root_dir
  git clone git@github.com:ACME-Climate/ACME.git $code_root_dir/$tag_name     # This will put repository, with all code, in directory $tag_name
  ## Setup git hooks
  rm -rf $code_root_dir/$tag_name/.git/hooks
  git clone git@github.com:ACME-Climate/ACME-Hooks.git $code_root_dir/$tag_name/.git/hooks         # checkout with write permission.
#  git clone git://github.com/ACME-Climate/ACME-Hooks.git .git/hooks      # checkout read-only.
  cd $code_root_dir/$tag_name
  git config commit.template $code_root_dir/$tag_name/.git/hooks/commit.template
  ## Bring in MPAS ocean/ice repo
  git submodule update --init

  if ( `lowercase $acme_tag` == master ) then
    echo ''
    ##echo 'run_acme: Detaching from the master branch to avoid accidental changes to master by user.'
    ##git checkout --detach
    echo 'KLUDGE: git version on anvil (1.7.1) is too old to be able to detach'
    echo 'edison uses git version 1.8.5.6 and it can git checkout --detach'
  else
    echo ''
    echo 'run_acme: Checking out branch ${acme_tag} = '${acme_tag}
    git checkout ${acme_tag}
  endif

endif

echo ''
echo 'run_acme: $case_name        = '$case_name



#===========================================
# DELETE PREVIOUS DIRECTORIES (IF REQUESTED)
#============================================

### Remove existing case_scripts directory (so it doesn't have to be done manually every time)
### Note: This script causes create_newcase to generate a temporary directory (part of a workaround to put the case_name into the script names)
###       If something goes wrong, this temporary directory is sometimes left behind, so we need to delete it too.
### Note: To turn off the deletion, set $num_seconds_until_delete to be negative.
###       To delete immediately, set $num_seconds_until_delete to be zero.

if ( -d $case_scripts_dir ) then
  if ( ${seconds_before_delete_case_dir} >= 0 ) then
    set num_seconds_until_delete = $seconds_before_delete_case_dir
    echo ''
    echo 'run_acme: Removing old $case_scripts_dir directory for '${case_name}' in '${num_seconds_until_delete}' seconds.'
    echo 'run_acme: To abort, press ctrl-C'
    while ( ${num_seconds_until_delete} > 0 )
      echo 'run_acme:  '${num_seconds_until_delete}'  seconds until deletion.'
      sleep 1
      @ num_seconds_until_delete = ${num_seconds_until_delete} - 1
    end
 #  ls -ld $case_scripts_dir     # For testing this script.
    rm -fr $case_scripts_dir
    echo 'run_acme:  Deleted $case_scripts_dir directory for : '${case_name}
  else
    echo 'run_acme: WARNING: $case_scripts_dir='$case_scripts_dir' exists '
    echo '          and is not being removed because seconds_before_delete_case_dir<0.'
    echo '          But create_newcase always fails when the case directory exists, so this script will now abort.'
    echo '          To fix this, either delete the case_scripts directory manually, or change seconds_before_delete_case_dir'
    exit 35
  endif
endif

#=============================================================
# HANDLE STANDARD PROCESSOR CONFIGURATIONS
#=============================================================
# NOTE: Some standard PE configurations are available (S,M,L,X1,X2).
#       If the requested configuration is 1 or CUSTOM, then set to M here, and handle later.

set lower_case = `lowercase $processor_config`
switch ( $lower_case )
  case 's':
    set std_proc_configuration = 'S'
    breaksw
  case 'm':
    set std_proc_configuration = 'M'
    breaksw
  case 'l':
    set std_proc_configuration = 'L'
    breaksw
  case 'x1':
    set std_proc_configuration = 'X1'
    breaksw
  case 'x2':
    set std_proc_configuration = 'X2'
    breaksw
  case '1':
    set std_proc_configuration = 'M'
    breaksw
  case 'custom*':
    # Note: this is just a placeholder so create_newcase will work.
    #       The actual configuration should be set under 'CUSTOMIZE PROCESSOR CONFIGURATION'
    set std_proc_configuration = 'M'
    breaksw
  default:
    echo 'run_acme ERROR: $processor_config='$processor_config' is not recognized'
    exit 40
    breaksw
endsw

set cime_dir = ${code_root_dir}/${tag_name}/cime
set create_newcase_exe = $cime_dir/scripts/create_newcase
if ( -f ${create_newcase_exe} ) then
  set acme_exe = acme.exe
  set case_setup_exe = $case_scripts_dir/case.setup
  set case_build_exe = $case_scripts_dir/case.build
  set case_run_exe = $case_scripts_dir/case.run
  set case_submit_exe = $case_scripts_dir/case.submit
  set xmlchange_exe = $case_scripts_dir/xmlchange
  set xmlquery_exe = $case_scripts_dir/xmlquery
  set shortterm_archive_script = $case_scripts_dir/case.st_archive
  set longterm_archive_script = $case_scripts_dir/case.lt_archive
else                                                                   # No version of create_newcase found
  echo 'run_acme ERROR: ${create_newcase_exe} not found'
  echo '                This is most likely because fetch_code should be true.'
  echo '                At the moment, $fetch_code = '$fetch_code
  exit 45
endif

#=============================================================
# CREATE CASE_SCRIPTS DIRECTORY AND POPULATE WITH NEEDED FILES
#=============================================================

set configure_options = "--case ${case_name} --compset ${compset} --res ${resolution} --project ${project} --pecount ${std_proc_configuration}"

if ( `lowercase $machine` != default ) then
  set configure_options = "$configure_options --mach ${machine}"
endif

echo ''
echo 'run_acme: -------- Starting create_newcase --------'
echo ''

echo $create_newcase_exe $configure_options
$create_newcase_exe $configure_options
cd ${case_name}

echo ''
echo 'run_acme: -------- Finished create_newcase --------'
echo ''

echo `pwd`

if ( `lowercase $machine` == default ) then
  set machine = `$xmlquery_exe MACH --value`
endif

set run_root_dir = `$xmlquery_exe CIME_OUTPUT_ROOT --value`

#================================
# SET WALLTIME FOR CREATE_NEWCASE
#================================

if ( `lowercase $walltime` == default ) then
  if ( `lowercase $debug_queue` == true ) then
    set walltime = '0:30:00'
  else
    if ( `lowercase $machine` == 'cab' || `lowercase $machine` == 'sierra' ) then
      set walltime = '1:29:00'
    else
      set walltime = '2:00:00'
    endif
  endif
endif

# Allow the user to specify how long the job taks
$xmlchange_exe JOB_WALLCLOCK_TIME=$walltime

#NOTE: Details of the configuration setup by create_newcase are in $case_scripts_dir/env_case.xml, which should NOT be edited.
#      It will be used by cesm_setup (formerly 'configure -case').
#NOTE: To get verbose output from create_newcase, add '-v' to the argument list.

#============================================================
#COPY STUFF TO SourceMods IF NEEDED:
#============================================================
#note: SourceMods is a horrible thing to do from a provenance perspective
#      it is much better to make changes to the code and to put those changes
#      into a git branch, which you then check out to do your run. Nonetheless,
#      it is useful for debugging so here's an example of how to use it...

#echo 'KLUDGE: Putting streams.ocean in SourceMods'
#cp /global/u1/p/petercal/junk/streams.ocean $case_scripts_dir/SourceMods/src.mpas-o/

#============================================================
#MAKE GROUP PERMISSIONS to $PROJECT FOR THIS DIRECTORY
#============================================================
#this stuff, combined with the umask command above, makes
#stuff in $run_root_dir readable by everyone in acme group.

if ( `lowercase $machine` == anvil ) then
  set run_root_dir = `$xmlquery_exe --value RUNDIR`
  chgrp climate $run_root_dir
  chmod g+s     $run_root_dir
  chgrp climate $case_scripts_dir
endif

#both run_root_dir and case_scripts_dir are created by create_newcase,
#so run_root_dir group isn't inherited for case_scripts_dir

#============================================================
# COPY THIS SCRIPT TO THE CASE DIRECTORY TO ENSURE PROVENANCE
#============================================================

set script_provenance_dir  = $case_scripts_dir/run_script_provenance
set script_provenance_name = $this_script_name.`date +%F_%T_%Z`
mkdir -p $script_provenance_dir
cp -f $this_script_path $script_provenance_dir/$script_provenance_name

#========================================================
# CREATE LOGICAL LINKS BETWEEN RUN_ROOT & THIS_SCRIPT_DIR
#========================================================

#NOTE: This is to make it easy for the user to cd to the case directory
#NOTE: Starting the suffix wit 'a' helps to keep this near the script in ls
#      (but in practice the behavior depends on the LC_COLLATE system variable).

echo 'run_acme: Creating logical links to make navigating easier.'
echo '          Note: Beware of using ".." with the links, since the behavior of shell commands can vary.'

# Link in this_script_dir case_dir
set run_dir_link = $this_script_dir/$this_script_name=a_run_link

echo ${run_dir_link}

if ( -l $run_dir_link ) then
  rm -f $run_dir_link
endif
echo "run_root ${run_root_dir}"
echo "run_dir ${run_dir_link}"

ln -s $run_root_dir $run_dir_link

#================================================
# COPY AUTO_CHAIN_RUNS SCRIPT TO CASE_SCRIPTS_DIR
#================================================

echo 'run_acme: copy auto_chain script, in case it is needed'

set auto_chain_run_file = ./auto_chain_runs.$machine
if ( -fx ${this_script_dir}/${auto_chain_run_file}  ) then
  echo 'run_acme: Copying '${auto_chain_run_file}' to '${case_scripts_dir}
  cp ${this_script_dir}/${auto_chain_run_file} ${case_scripts_dir}/${auto_chain_run_file}
endif

#=============================================
# CUSTOMIZE PROCESSOR CONFIGURATION
# ============================================
#NOTE: Changes to the processor configuration should be done by an expert.  \
#      Not all possible options will work.

if ( `lowercase $processor_config` == '1' ) then

  # NOTE: xmlchange won't work with shell variables for the id, so we have to write it out in full.
  set ntasks = 1
  set nthrds = 1
  set sequential_or_concurrent = 'sequential'
  foreach ntasks_name ( NTASKS_ATM  NTASKS_LND  NTASKS_ICE  NTASKS_OCN  NTASKS_CPL  NTASKS_GLC  NTASKS_ROF  NTASKS_WAV )
    $xmlchange_exe --id $ntasks_name --val $ntasks
  end

  foreach nthrds_name ( NTHRDS_ATM  NTHRDS_LND  NTHRDS_ICE  NTHRDS_OCN  NTHRDS_CPL  NTHRDS_GLC  NTHRDS_ROF  NTHRDS_WAV )
    $xmlchange_exe --id $nthrds_name --val $nthrds
  end

  foreach rootpe_name ( ROOTPE_ATM  ROOTPE_LND  ROOTPE_ICE  ROOTPE_OCN  ROOTPE_CPL  ROOTPE_GLC  ROOTPE_ROF  ROOTPE_WAV )
    $xmlchange_exe --id $rootpe_name --val 0
  end

  foreach layout_name ( NINST_ATM_LAYOUT NINST_LND_LAYOUT NINST_ICE_LAYOUT NINST_OCN_LAYOUT NINST_GLC_LAYOUT NINST_ROF_LAYOUT NINST_WAV_LAYOUT )
    $xmlchange_exe --id $layout_name --val $sequential_or_concurrent
  end

else if ( `lowercase $processor_config` == 'custom173' ) then

  echo 'run_acme: Setting custom processor configuration, because $processor_config = '$processor_config
###   This space is to allow a custom processor configuration to be defined.
###   If your layout will be useful to other people, then please get it added to the standard
###   configurations in the ACME repository.

###   NOTE: It is shorter and more robust to implement the PE configuration changes using xmlchange
###         rather than embedding the whole env_mach_pes.xml file

# These were chosen based on the smallest pe layouts of provided input files
  $xmlchange_exe --id NTASKS_ATM --val 2700
  $xmlchange_exe --id NTHRDS_ATM --val 1
  $xmlchange_exe --id ROOTPE_ATM --val 0

  $xmlchange_exe --id NTASKS_LND --val 312
  $xmlchange_exe --id NTHRDS_LND --val 1
  $xmlchange_exe --id ROOTPE_LND --val 2400

  $xmlchange_exe --id NTASKS_ICE --val 2400
  $xmlchange_exe --id NTHRDS_ICE --val 1
  $xmlchange_exe --id ROOTPE_ICE --val 0

  $xmlchange_exe --id NTASKS_OCN --val 1440
  $xmlchange_exe --id NTHRDS_OCN --val 1
  $xmlchange_exe --id ROOTPE_OCN --val 2712

  $xmlchange_exe --id NTASKS_CPL --val 2400
  $xmlchange_exe --id NTHRDS_CPL --val 1
  $xmlchange_exe --id ROOTPE_CPL --val 0

  $xmlchange_exe --id NTASKS_GLC --val 312
  $xmlchange_exe --id NTHRDS_GLC --val 1
  $xmlchange_exe --id ROOTPE_GLC --val 2400

  $xmlchange_exe --id NTASKS_ROF --val 312
  $xmlchange_exe --id NTHRDS_ROF --val 1
  $xmlchange_exe --id ROOTPE_ROF --val 2400

  $xmlchange_exe --id NTASKS_WAV --val 2400
  $xmlchange_exe --id NTHRDS_WAV --val 1
  $xmlchange_exe --id ROOTPE_WAV --val 0

  set sequential_or_concurrent = 'concurrent'
# CIME 5.1 has a typo in it's env_mach_pes.xml such that NINST_ICE_LAYOUT can be 'concurrent' or 'equential' (note the missing s)
  $xmlchange_exe --id NINST_ICE_LAYOUT --val concurrent
#  foreach layout_name ( NINST_ATM_LAYOUT NINST_LND_LAYOUT NINST_ICE_LAYOUT NINST_OCN_LAYOUT NINST_GLC_LAYOUT NINST_ROF_LAYOUT NINST_WAV_LAYOUT )
  foreach layout_name ( NINST_ATM_LAYOUT NINST_LND_LAYOUT NINST_OCN_LAYOUT NINST_GLC_LAYOUT NINST_ROF_LAYOUT NINST_WAV_LAYOUT )
    $xmlchange_exe --id $layout_name --val $sequential_or_concurrent
  end

else if ( `lowercase $processor_config` == 'custom375' ) then

  echo 'run_acme: Setting custom processor configuration, because $processor_config = '$processor_config
###   This space is to allow a custom processor configuration to be defined.
###   If your layout will be useful to other people, then please get it added to the standard
###   configurations in the ACME repository.

###   NOTE: It is shorter and more robust to implement the PE configuration changes using xmlchange
###         rather than embedding the whole env_mach_pes.xml file

# These were chosen based on the smallest pe layouts of provided input files
  $xmlchange_exe --id NTASKS_ATM --val 5400
  $xmlchange_exe --id NTHRDS_ATM --val 1
  $xmlchange_exe --id ROOTPE_ATM --val 0

  $xmlchange_exe --id NTASKS_LND --val 600
  $xmlchange_exe --id NTHRDS_LND --val 1
  $xmlchange_exe --id ROOTPE_LND --val 4800

  $xmlchange_exe --id NTASKS_ICE --val 4800
  $xmlchange_exe --id NTHRDS_ICE --val 1
  $xmlchange_exe --id ROOTPE_ICE --val 0

  $xmlchange_exe --id NTASKS_OCN --val 3600
  $xmlchange_exe --id NTHRDS_OCN --val 1
  $xmlchange_exe --id ROOTPE_OCN --val 5400

  $xmlchange_exe --id NTASKS_CPL --val 4800
  $xmlchange_exe --id NTHRDS_CPL --val 1
  $xmlchange_exe --id ROOTPE_CPL --val 0

  $xmlchange_exe --id NTASKS_GLC --val 600
  $xmlchange_exe --id NTHRDS_GLC --val 1
  $xmlchange_exe --id ROOTPE_GLC --val 4800

  $xmlchange_exe --id NTASKS_ROF --val 600
  $xmlchange_exe --id NTHRDS_ROF --val 1
  $xmlchange_exe --id ROOTPE_ROF --val 4800

  $xmlchange_exe --id NTASKS_WAV --val 4800
  $xmlchange_exe --id NTHRDS_WAV --val 1
  $xmlchange_exe --id ROOTPE_WAV --val 0

  set sequential_or_concurrent = 'concurrent'
# CIME 5.1 has a typo in it's env_mach_pes.xml such that NINST_ICE_LAYOUT can be 'concurrent' or 'equential' (note the missing s)
  $xmlchange_exe --id NINST_ICE_LAYOUT --val concurrent
#  foreach layout_name ( NINST_ATM_LAYOUT NINST_LND_LAYOUT NINST_ICE_LAYOUT NINST_OCN_LAYOUT NINST_GLC_LAYOUT NINST_ROF_LAYOUT NINST_WAV_LAYOUT )
  foreach layout_name ( NINST_ATM_LAYOUT NINST_LND_LAYOUT NINST_OCN_LAYOUT NINST_GLC_LAYOUT NINST_ROF_LAYOUT NINST_WAV_LAYOUT )
    $xmlchange_exe --id $layout_name --val $sequential_or_concurrent
  end


### The following couple of lines are for when no custom configuration is set (eg, in the archived version)
#    echo 'run_acme: Custom processor configuration not defined.  Please edit this script.'
#    exit 150

endif

#============================================
# SET PARALLEL I/O (PIO) SETTINGS
#============================================
#Having bad PIO_NUMTASKS and PIO_STRIDE values can wreck performance
#See https://acme-climate.atlassian.net/wiki/display/ATM/How+to+Create+a+New+PE+Layout

#$xmlchange_exe -file env_run.xml -id PIO_NUMTASKS -val 128

#============================================
# SET MODEL INPUT DATA DIRECTORY
#============================================
# NOTE: This section was moved from later in script, because sometimes it is needed for cesm_setup.

# The model input data directory should default to the managed location for your system.
# However, if this does not work properly, or if you want to use your own data, then
# that should be setup here (before case_scripts.build because it checks the necessary files exist)

# NOTE: This code handles the case when the default location is wrong.
#       If you want to use your own files then this code will need to be modified.

# NOTE: For information on the ACME input data repository, see:
#       https://acme-climate.atlassian.net/wiki/display/WORKFLOW/ACME+Input+Data+Repository

#set input_data_dir = 'input_data_dir_NOT_SET'
#if ( $machine == 'cori' || $machine == 'edison' ) then
#  set input_data_dir = '/project/projectdirs/m411/ACME_inputdata'    # PJC-NERSC
## set input_data_dir = '/project/projectdirs/ccsm1/inputdata'        # NERSC
#else if ( $machine == 'titan' || $machine == 'eos' ) then
#  set input_data_dir = '/lustre/atlas/proj-shared/cli112/pjcs/ACME_inputdata'    # PJC-OLCF
#endif
#if ( -d  $input_data_dir ) then
#  $xmlchange_exe --id DIN_LOC_ROOT --val $input_data_dir
#else
#  echo 'run_acme ERROR: User specified input data directory does NOT exist.'
#  echo '                $input_data_dir = '$input_data_dir
#  exit 270
#endif

### The following command extracts and stores the input_data_dir in case it is needed for user edits to the namelist later.
### NOTE: The following line may be necessary if the $input_data_dir is not set above, and hence defaults to the ACME default.
#NOTE: the following line can fail for old versions of the code (like v0.3) because
#"-value" is a new option in xmlquery. If that happens, comment out the next line and
#hardcode in the appropriate DIN_LOC_ROOT value for your machine.
set input_data_dir = `$xmlquery_exe DIN_LOC_ROOT --value`

#============================================
# COMPONENT CONFIGURATION OPTIONS
#============================================

#NOTE:  This section is for making specific component configuration selections.
#NOTE:  The input_data directory is best set in the section for it above.
#NOTE:  Setting CAM_CONFIG_OPTS will REPLACE anything set by the build system.
#       To add on instead, add '-append' to the xmlchange command.
#NOTE:  CAM_NAMELIST_OPTS should NOT be used.  Instead, use the user_nl section after case_scripts.build

#$xmlchange_exe --id CAM_CONFIG_OPTS --val "-phys cam5 -chem linoz_mam3"

## Chris Golaz: build with COSP
#./xmlchange -file env_build.xml -id CAM_CONFIG_OPTS -append -val "-cosp"

#===========================
# SET THE PARTITION OF NODES
#===========================

if ( `lowercase $debug_queue` == true ) then
  if ( $machine == cab || $machine == sierra ) then
    $xmlchange_exe --id JOB_QUEUE --val 'pdebug'
  else
    $xmlchange_exe --id JOB_QUEUE --val 'debug'
  endif
endif

#============================================
# CONFIGURE
#============================================

#note configure -case turned into cesm_setup in cam5.2

echo ''
echo 'run_acme: -------- Starting case.setup --------'
echo ''

echo $case_setup_exe

$case_setup_exe --reset

echo ''
echo 'run_acme: -------- Finished case.setup  --------'
echo ''

#============================================
# SET BUILD OPTIONS
#============================================

if ( `uppercase $debug_compile` != 'TRUE' && `uppercase $debug_compile` != 'FALSE' ) then
  echo 'run_acme ERROR: $debug_compile can be true or false but is instead '$debug_compile
  exit 220
endif

if ( `lowercase $machine` == 'edison' && `uppercase $debug_compile` == 'TRUE' ) then
  echo 'run_acme ERROR: Edison currently has a compiler bug and crashes when compiling in debug mode (Nov 2015)'
  exit 222
endif

$xmlchange_exe --id DEBUG --val `uppercase $debug_compile`

#Modify/uncomment the next line to change the number of processors used to compile.
#$xmlchange_exe --id GMAKE_J --val 4

#=============================================
# CREATE NAMELIST MODIFICATION FILES (USER_NL)
#=============================================

# Append desired changes to the default namelists generated by the build process.
#
# NOTE: It doesn't matter which namelist an option is in for any given component, the system will sort it out.
# NOTE: inputdata directory ($input_data_dir) is set above (before cesm_setup).
# NOTE: The user_nl files need to be set before the build, because case_scripts.build checks whether input files exist.
# NOTE: $atm_output_freq and $records_per_atm_output_file are so commonly used, that they are set in the options at the top of this script.

cat <<EOF >> user_nl_cam
 nhtfrq = $atm_output_freq
 mfilt  = $records_per_atm_output_file
EOF

### NOTES ON COMMON NAMELIST OPTIONS ###

### ATMOSPHERE NAMELIST ###

#NHTFRQ : The frequency with which the atmosphere writes its output.
#    0=monthly, +N=every N timesteps,  -N=every N hours
#    For more information:   http://www.cesm.ucar.edu/models/atm-cam/docs/usersguide/node45.html
#MFILT : The number of time records in each netCDF output file from the atmosphere model.
#    If mfilt is 0, then there will only be one time record per file.
#NOTE:  nhtfrq and mfilt can be a comma-separated list of numbers, corresponding to the 'history tapes' defined in the namelist.

#============================================
# BUILD CODE
#============================================

#NOTE: This will either build the code (if needed and $old_executable=false) or copy an existing executable.

echo ''
echo 'run_acme: -------- Starting Build --------'
echo ''

echo ${case_build_exe}
${case_build_exe}

echo ''
echo 'run_acme: -------- Finished Build --------'
echo ''


#============================================
# QUEUE OPTIONS
#============================================
# Edit the default queue and batch job lengths.

#HINT: To change queue after run submitted, the following works on most machines:
#      qalter -lwalltime=00:29:00 <run_descriptor>
#      qalter -W queue=debug <run_descriptor>

#NOTE: we are currently not modifying the archiving scripts to run in debug queue when $debug_queue=true
#      under the assumption that if you're debugging you shouldn't be archiving.

#NOTE: there was 1 space btwn MSUB or PBS and the commands before cime and there are 2 spaces
#      in post-cime versions. This is fixed by " \( \)*" in the lines below. The "*" here means
#      "match zero or more of the expression before". The expression before is a single whitespace.
#      The "\(" and "\)" bit indicate to sed that the whitespace in between is the expression we
#      care about. The space before "\(" makes sure there is at least one whitespace after #MSUB.
#      Taken all together, this stuff matches lines of the form #MSUB<one or more whitespaces>-<stuff>.

set machine = `lowercase $machine`

### Only specially authorized people can use the special_acme qos on Cori or Edison. Don't uncomment unless you're one.
### if ( `lowercase $debug_queue` == false && ( $machine == core || $machine == edison ) ) then
###   sed -i /"#SBATCH \( \)*--partition"/a"#SBATCH  --qos=special_acme"  ${case_run_exe}
### endif

#============================================
# BATCH JOB OPTIONS
#============================================

# Set options for batch scripts (see above for queue and batch time, which are handled separately)

# NOTE: This also modifies the short-term and long-term archiving scripts.
# NOTE: We want the batch job log to go into a sub-directory of case_scripts (to avoid it getting clogged up)

mkdir -p batch_output      ### Make directory that stdout and stderr will go into.

if ( $machine == hopper ) then
    sed -i /"#PBS \( \)*-N"/c"#PBS  -N ${job_name}"                                ${case_run_exe}
    sed -i /"#PBS \( \)*-j oe"/a'#PBS  -o batch_output/${PBS_JOBNAME}.o${PBS_JOBID}' ${case_run_exe}
    sed -i /"#PBS \( \)*-N"/c"#PBS  -N ST+${job_name}"                             $shortterm_archive_script
    sed -i /"#PBS \( \)*-j oe"/a'#PBS  -o batch_output/${PBS_JOBNAME}.o${PBS_JOBID}' $shortterm_archive_script
    sed -i /"#PBS \( \)*-N"/c"#PBS  -N LT+${job_name}"                             $longterm_archive_script
    sed -i /"#PBS \( \)*-j oe"/a'#PBS  -o batch_output/${PBS_JOBNAME}.o${PBS_JOBID}' $longterm_archive_script

else if ( $machine == cori || $machine == edison ) then
    sed -i /"#SBATCH \( \)*--job-name"/c"#SBATCH  --job-name=${job_name}"                 ${case_run_exe}
    sed -i /"#SBATCH \( \)*--job-name"/a"#SBATCH  --account=${project}"                   ${case_run_exe}
    sed -i /"#SBATCH \( \)*--output"/c"#SBATCH  --output=batch_output/"${case_name}'.o%j' ${case_run_exe}

    sed -i /"#SBATCH \( \)*--job-name"/c"#SBATCH  --job-name=ST+${job_name}"                  $shortterm_archive_script
    sed -i /"#SBATCH \( \)*--job-name"/a"#SBATCH  --account=${project}"                       $shortterm_archive_script
    sed -i /"#SBATCH \( \)*--output"/c'#SBATCH  --output=batch_output/ST+'${case_name}'.o%j'  $shortterm_archive_script
    sed -i /"#SBATCH \( \)*--job-name"/c"#SBATCH  --job-name=LT+${job_name}"                  $longterm_archive_script
    sed -i /"#SBATCH \( \)*--job-name"/a"#SBATCH  --account=${project}"                       $longterm_archive_script
    sed -i /"#SBATCH \( \)*--output"/c'#SBATCH  --output=batch_output/LT+'${case_name}'.o%j'  $longterm_archive_script

else if ( $machine == titan || $machine == eos ) then
    sed -i /"#PBS \( \)*-N"/c"#PBS  -N ${job_name}"                                ${case_run_exe}
    sed -i /"#PBS \( \)*-A"/c"#PBS  -A ${project}"                                 ${case_run_exe}
    sed -i /"#PBS \( \)*-j oe"/a'#PBS  -o batch_output/${PBS_JOBNAME}.o${PBS_JOBID}' ${case_run_exe}

    sed -i /"#PBS \( \)*-N"/c"#PBS  -N ST+${job_name}"                             $shortterm_archive_script
    sed -i /"#PBS \( \)*-j oe"/a'#PBS  -o batch_output/${PBS_JOBNAME}.o${PBS_JOBID}' $shortterm_archive_script
    sed -i /"#PBS \( \)*-N"/c"#PBS  -N LT+${job_name}"                             $longterm_archive_script
    sed -i /"#PBS \( \)*-j oe"/a'#PBS  -o batch_output/${PBS_JOBNAME}.o${PBS_JOBID}' $longterm_archive_script

else
    echo 'run_acme WARNING: This script does not have batch directives for $machine='$machine
    echo '                  Assuming default ACME values.'
endif

#============================================
# SETUP SHORT AND LONG TERM ARCHIVING
#============================================

$xmlchange_exe --id DOUT_S    --val `uppercase $do_short_term_archiving`

$xmlchange_exe --id DOUT_L_MS --val `uppercase $do_long_term_archiving`

# DOUT_L_MSROOT is the directory in your account on the local mass storage system (typically an HPSS tape system)
$xmlchange_exe --id DOUT_L_MSROOT --val "ACME_simulation_output/${case_name}"


#============================================
# COUPLER HISTORY OUTPUT
#============================================

#$xmlchange_exe --id HIST_OPTION --val ndays
#$xmlchange_exe --id HIST_N      --val 1


#=======================================================
# SETUP SIMULATION LENGTH AND FREQUENCY OF RESTART FILES
#=======================================================

#SIMULATION LENGTH
$xmlchange_exe --id STOP_OPTION --val `lowercase $stop_units`
$xmlchange_exe --id STOP_N      --val $stop_num

#RESTART FREQUENCY
$xmlchange_exe --id REST_OPTION --val `lowercase $restart_units`
$xmlchange_exe --id REST_N      --val $restart_num

#COUPLER BUDGETS
$xmlchange_exe --id BUDGETS     --val `uppercase $do_cpl_hist`

#COUPLER HISTORY FILES
$xmlchange_exe --id HIST_OPTION --val `lowercase $cpl_hist_units`
$xmlchange_exe --id HIST_N      --val $cpl_hist_num

#============================================
# SETUP SIMULATION INITIALIZATION
#============================================

echo ''
echo 'run_acme: $model_start_type = '${model_start_type}'  (This is NOT necessarily related to RUN_TYPE)'

set model_start_type = `lowercase $model_start_type`
#-----------------------------------------------------------------------------------------------
# start_type = initial means start a new run from default or user-specified initial conditions
#-----------------------------------------------------------------------------------------------
if ( $model_start_type == 'initial' ) then
  ### 'initial' run: cobble together files, with RUN_TYPE= 'startup' or 'hybrid'.
  $xmlchange_exe --id RUN_TYPE --val "startup"
  $xmlchange_exe --id CONTINUE_RUN --val "FALSE"

  # if you want to use your own initial conditions, uncomment and fix up the lines below:
#  set initial_files_dir = $PROJWORK/cli107/sulfur_DOE_restarts/2deg_1850_0011-01-01-00000
#  cp -fpu $initial_files_dir/* ${case_run_dir}

#-----------------------------------------------------------------------------------------------
# start_type = continue means you've already run long enough to produce restart files and want to
# continue the run
#-----------------------------------------------------------------------------------------------
else if ( $model_start_type == 'continue' ) then

  ### This is a standard restart.

  $xmlchange_exe --id CONTINUE_RUN --val "TRUE"

#-----------------------------------------------------------------------------------------------
# start_type = branch means you want to continue a run, but in a new run directory and/or with
# recompiled code
#-----------------------------------------------------------------------------------------------
else if ( $model_start_type == 'branch' ) then

  ### Branch runs are the same as restarts, except that the history output can be changed
  ### (eg to add new variables or change output frequency).

  ### Branch runs are often used when trying to handle a complicated situation.
  ### Hence, it is likely that the user will need to customize this section.

  ### the next lines attempt to automatically extract all needed info for setting up the branch run.
  set rpointer_filename = "${restart_files_dir}/rpointer.drv"
  if ( ! -f $rpointer_filename ) then
    echo 'run_acme ERROR: ${rpointer_filename} does not exist. It is needed to extract RUN_REFDATE.'
    echo "                This may be because you should set model_start_type to 'initial' or 'continue' rather than 'branch'."
    echo '                ${rpointer_filename} = '{rpointer_filename}
    exit 370
  endif
  set restart_coupler_filename = `cat $rpointer_filename`
  set restart_case_name = ${restart_coupler_filename:r:r:r:r}         # Extract out the case name for the restart files.
  set restart_filedate = ${restart_coupler_filename:r:e:s/-00000//}   # Extract out the date (yyyy-mm-dd).
  echo 'run_acme: $restart_case_name = '$restart_case_name
  echo 'run_acme: $restart_filedate  = '$restart_filedate

  ### the next line gets the YYYY-MM of the month before the restart time. Needed for staging history files.
  set restart_prevdate = `date -d "${restart_filedate} - 1 month" +%Y-%m`

  echo 'run_acme: $restart_prevdate  = '$restart_prevdate

  echo 'run_acme: Copying stuff for branch run'
  cp ${restart_files_dir}/${restart_case_name}.cam.r.${restart_filedate}-00000.nc $case_run_dir
  cp ${restart_files_dir}/${restart_case_name}.cam.rs.${restart_filedate}-00000.nc $case_run_dir
  cp ${restart_files_dir}/${restart_case_name}.clm2.r.${restart_filedate}-00000.nc $case_run_dir
  cp ${restart_files_dir}/${restart_case_name}.clm2.rh0.${restart_filedate}-00000.nc $case_run_dir
  cp ${restart_files_dir}/${restart_case_name}.cpl.r.${restart_filedate}-00000.nc $case_run_dir
  cp ${restart_files_dir}/${restart_case_name}.mosart.r.${restart_filedate}-00000.nc $case_run_dir
  cp ${restart_files_dir}/${restart_case_name}.mosart.rh0.${restart_filedate}-00000.nc $case_run_dir
  cp ${restart_files_dir}/rst.mpas-cice.${restart_filedate}_00.00.00.nc $case_run_dir
  cp ${restart_files_dir}/rst.mpas-o.${restart_filedate}_00.00.00.nc $case_run_dir
  cp ${restart_files_dir}/${restart_case_name}.cam.h0.${restart_prevdate}.nc $case_run_dir
  cp ${restart_files_dir}/${restart_case_name}.mosart.h0.${restart_prevdate}.nc $case_run_dir
  cp ${restart_files_dir}/${restart_case_name}.clm2.h0.${restart_prevdate}.nc $case_run_dir
  cp ${restart_files_dir}/rpointer* $case_run_dir

  $xmlchange_exe --id RUN_TYPE --val "branch"
  $xmlchange_exe --id RUN_REFCASE --val $restart_case_name
  $xmlchange_exe --id RUN_REFDATE --val $restart_filedate    # Model date of restart file
  $xmlchange_exe --id CONTINUE_RUN --val "FALSE"
  $xmlchange_exe --id BRNCH_RETAIN_CASENAME --val "FALSE"  ## Only TRUE if you want to continue the run with the same name (risky)!!

else

  echo 'run_acme ERROR: $model_start_type = '${model_start_type}' is unrecognized.   Exiting.'
  exit 380

endif

#============================================
# RUN CONFIGURATION OPTIONS
#============================================

#NOTE:  This section is for making specific changes to the run options (ie env_run.xml).

#if ( $machine == cori ) then      ### fix pnetcdf problem on Cori. (github #593)
#  $xmlchange_exe --id PIO_TYPENAME  --val "netcdf"
#endif

#=================================================
# SUBMIT THE SIMULATION TO THE RUN QUEUE
#=================================================
#note: to run the model in the totalview debugger,
# cd $case_run_dir
# totalview srun -a -n <number of procs> -p <name of debug queue> ../bld/$acme_exe
# where you may need to change srun to the appropriate submit command for your system, etc.


echo ''
echo 'run_acme: -------- Starting Submission to Run Queue --------'
echo ''

if ( `lowercase $submit_run` == 'true' ) then
  if ( $num_submits == 1 ) then
    echo 'run_acme:          SUBMITTING A SINGLE JOB.'
    ${case_submit_exe}
  else if ( $num_submits <= 0 ) then
    echo 'run_acme: $num_submits <= 0 so NOT submitting a job.'
    echo '          $num_submits = '$num_submits
  else if ( `lowercase $debug_queue` == 'true' && $num_submits > 1 ) then
    echo 'run_acme WARNING: $num_submits > 1  and  $debug_queue = "TRUE"'
    echo '                  Submitting chained jobs to the debug queue is usually forbidden'
    echo '                  $num_submits = '$num_submits
    echo '                  SUBMITTING JUST A SINGLE JOB.'
    ${case_submit_exe}
  else if ( ! -x ./auto_chain_runs.$machine && $num_submits > 1 ) then
    echo 'run_acme WARNING: $num_submits > 1  but auto_chain_runs.$machine excutable cannot be found.'
    echo '                  $num_submits = '$num_submits
    echo '                  $machine     = '$machine
    echo '                  SUBMITTING JUST A SINGLE JOB.'
    ${case_submit_exe}
  else
    echo 'run_acme: executing 'auto_chain_runs.$machine
    echo '          $num_submits = '$num_submits
    echo '          $do_short_term_archiving = '`uppercase $do_short_term_archiving`
    echo '          $do_long_term_archiving  = '`uppercase $do_long_term_archiving`
    # To avoid the error checking in the ACME scripts, it is necessary to tell ACME the archiving is FALSE, and then implement it manually.
    $xmlchange_exe --id DOUT_S    --val 'FALSE'
    $xmlchange_exe --id DOUT_L_MS --val 'FALSE'
    ./auto_chain_runs.$machine  $num_submits -1 `uppercase $do_short_term_archiving` `uppercase $do_long_term_archiving` ${case_run_exe}
  endif
else
    echo 'run_acme: Run NOT submitted because $submit_run = '$submit_run
endif

echo ''
echo 'run_acme: -------- Finished Submission to Run Queue --------'
echo ''

#=================================================
# DO POST-SUBMISSION THINGS (IF ANY)
#=================================================

# Actions after the run submission go here.

echo ''
echo 'run_acme: ++++++++ run_acme Completed ('`date`') ++++++++'
echo ''

#**********************************************************************************
### --- end of script - there are no commands beyond here, just useful comments ---
#**********************************************************************************

### -------- Version information --------
# 1.0.0    2015-11-19    Initial version.  Tested on Titan. (PJC)
# 1.0.1    2015-11-19    Fixed bugs and added features  for Hopper. (PJC)
# 1.0.2    2015-11-19    Modified to conform with ACME script standards. PJC)
# 1.0.3    2015-11-23    Modified to include Peter's ideas (PMC)
# 1.0.4    2015-11-23    Additional modification based on discusion with Peter and Chris Golaz. (PJC)
# 1.0.5    2015-11-23    Tweaks for Titan (PJC)
# 1.0.6    2015-11-23    Fixed some error messages, plus some other minor tweaks. (PJC)
# 1.0.7    2015-11-24    Fixed bug for setting batch options (CIME adds an extra space than before). (PJC)
#                        Also, removed GMAKE_J from option list (left it as a comment for users to find).
# 1.0.8    2015-11-24    Merged old_executable stuff and changed to loop over xmlchange statements for
#                        single proc run (PMC)
# 1.0.9    2015-11-25    Added support for using pre-cime code versions, fixed some bugs (PMC)
# 1.0.10   2015-11-25    Cosmetic changes to the edited batch script, and improved comments. (PJC)
# 1.0.11   2015-11-25    Fixed bug with naming of st_archive and lt_archive scripts.  Also cosmetic improvements (PJC).
# 1.0.12   2015-11-25    changed name of variable orig_dir/dir_of_this_script to "this_script_dir" and removed options
#                        for old_executable=true and seconds_before_delete_case_dir<0 because they were provenance-unsafe.(PMC)
# 1.0.13   2015-11-25    Merged changes from PMC with cosmetic changes from PJC.
#                        Also, reactivated old_executable=true, because the model recompiles many files unnecessarily.  (PJC)
# 1.0.14   2015-11-25    Added custom PE configuration so the ACME pre-alpha code will work on Titan for Chris Golaz.   (PJC)
#                        Fixed $cime_space bug introduced in 1.0.10 (PJC)
# 1.0.15   2015-11-25    Fixed bug with old_executable=true (PJC)
# 1.0.16   2015-11-30    Added $machine to the case_name (PJC)
# 1.0.17   2015-11-30    Added date to filename when archiving this script (so previous version doesn't get overwritten) (PJC)
# 1.0.18   2015-11-30    Will now automatically use 'git checkout --detach' so users cannot alter master by accident (PJC)
# 1.0.19   ??            Added an option to set the directory for short term archiving.  Also fixed some comments. (PMC)
# 1.0.20   2015-12-10    Improved comments, especially for 'old_executable' option. (PJC)
# 1.0.21   2015-12-10    Modified so that the script names contain "$case_name" rather than "case_scripts".
#                        Create_newcase doesn't have the flexibility to do what we need, and the rest of the CESM scripts
#                        are designed to stop us doing what we want, so we had to defeat those protections, but
#                        we do this in a safe way that reinstates the protections. (PJC)
# 1.0.22   2015-12-11    Creates logical links so it is easy to move between this this_script_dir and run_root_dir. (PJC)
# 1.0.23   2015-12-11    Changed references to build_and_run_script to just run_script, for consistency and brevity. (PJC)
# 1.0.24   2015-12-11    The temp_case_scripts_dir is now handled like case_scripts_dir for checking and deletion.  (PJC)
# 1.0.25   2015-12-11    Can have separate name for batch scheduler, to help distinguish runs. (PJC)
# 1.0.26   2015-12-16    Can now handle Cori (NERSC), plus improved error messages.  (PJC)
# 1.0.27   2015-12-16    Partial implementation for Eos (OLCF), plus cosmetic changes.  (PJC)
# 1.0.28   2015-12-17    Fixed Cori batch options.  Improved an error message.  (PJC)
# 1.0.29   2015-12-21    Added line to extract inputdata_dir from XML files, so it is available if needed in user_nl files. (PJC)
# 1.0.30   2015-12-23    Changed run.output dir to batch_output -- purpose is clearer, and better for filename completion. (PJC)
#                        Added option to set PE configuration to sequential or concurrent for option '1'. (PJC)
# 1.0.31   2016-01-07    Moved up the location where the input_data_dir is set, so it is availble to cesm_setup.    (PJC)
#                        Checks case_name is 79 characters, or less, which is a requirement of the ACME scripts.
#                        Improved options for SLURM machines.
#                        Added numbers for the ordering of options at top of file (in preperation for reordering).
#                        Added xxdiff calls to fix known bugs in master> (need to generalize for other people)
# 1.0.32   2016-01-07    Converted inputdata_dir to input_data_dir for consistency.      (PJC)
#                        Cosmetic improvements.
# 1.0.33   2016-01-08    Changed default tag to master_detached to improve clarity. (PJC)
#                        Now sets up ACME git hooks when fetch_code=true.
# 1.0.33p  2016-01-08    Changed compset from A_B1850CN to A_B1850 (pre-acme script only).  (PJC)
#                        Added finidat = '' to user_nl_clm, which allows A_B1850 to run.
# 1.0.34   2016-01-12    Commented out the input_data_dir user configuration, so it defaults to the ACME settings.   (PJC)
# 1.0.35   2016-01-13    Improved an error message.   (PJC)
# 1.0.36   2016-01-21    Reordered options to better match workflow. (PJC)
# 1.2.0    2016-01-21    Set options to settings for release. (PJC)
# 1.2.1    2016-01-21    Reordered and refined comments to match new ordering of options. (PJC)
# 1.2.2    2016-01-21    The batch submission problem on Cori has been repaired on master (#598),
#                        so I have undone the workaround in this script. (PJC)
# 1.2.3    2016-01-26    Commented out some of the workarounds for ACME bugs that are no longer needed.  (PJC)
# 1.4.0    2016-03-23    A number of modifications to handle changes in machines and ACME. [version archived to ACME] (PJC)
# 1.4.1    2016-03-23    Modified to defaults for Cori (NERSC). (PJC)
# 1.4.2    2016-08-05    Replaced cime_space with pattern matching, added num_depends functionality for daisychained
#                        jobs, added code for submitting to qos=acme_special on Edison, added cpl_hist options, and
#                        improved support for sierra and cab at LLNL.(PMC)
# 1.4.3    2016-08-10    Improved support for branch runs (PMC)
# 1.4.4    2016-08-11    Added umask command to make run directory world-readable by default.
# 2.0.0    2016-08-10    Added capability to a chain of submissions using the script auto_chain_runs.$machine (PJC)
# 2.0.1    2016-09-13    Fixed num_resubmits undefined error.
#                        Generalized setting of group permissions for other machines. (PJC)
# 2.0.2    2016-09-13    Turned off short- and long-term archiving so auto_chain_runs script can do it manually. (PJC)
# 2.0.3    2016-09-14    Removed 'git --set-upstream' command, because it does not work on tags.  (PJC)
# 2.0.4    2016-09-14    Long term archiving not working in ACME, so turn it off and warn user.  (PJC)
# 3.0.0    2016-12-15    Initial update for CIME5. Change script names, don't move the case directory
#                        as it's broken by the update, use xmlchange to set up the custom PE Layout.
#                        Remove support for CIME2 and pre-cime versions.  (MD)
# 3.0.1    2017-01-26    Setup to run A_WCYCL1850S simulation at ne30 resolution.  (CG) 
# 3.0.2    2017-02-13    Activated logical links by default, and tweaked the default settings. (PJC)

# NOTE:  PJC = Philip Cameron-Smith,  PMC = Peter Caldwell, CG = Chris Golaz, MD = Michael Deakin

### ---------- Desired features still to be implemented ------------
# +) fetch_code = update   (pull in latest updates to branch)    (PJC)
# +) A way to run the testsuite.? (PJC)
# +) Reorder options at top to match workflow. (PJC)
# +) make the handling of lowercase consistent.  $machine may need to be special. (PJC)
# +) generalize xxdiff commands (for fixing known bugs) to work for other people  (PJC)
# +) Add a 'default' option, for which REST_OPTION='$STOP_OPTION' and REST_N='$STOP_N'.
#    This is important if the user subsequently edits STOP_OPTION or STOP_N.      (PJC)
# +) Add defaults for Edison. (PJC)
# +) triggering on $acme_tag = master_detached doesn't make sense.  Fix logic. (PJC)

###Example sed commands
#============================
###To delete a line
#sed -i /'fstrat_list'/d $namelists_dir/cam.buildnml.csh

### To replace part of a line
#sed -i s/"#PBS -q regular"/"#PBS -q debug"/ ${case_run_exe}

### To replace a whole line based on a partial match
#sed -i /"#PBS -N"/c"#PBS -N ${run_job_name}" ${case_run_exe}

### To add a new line:
# sed -i /"PBS -j oe"/a"#PBS -o batch_output/${PBS_JOBNAME}.o${PBS_JOBID}" ${case_run_exe}

