#! /bin/csh -fe
### This script was created 2015-11-15 by Philip Cameron-Smith (pjc@llnl.gov) and Peter Caldwell
### and incorporates some features originally from Hui Wan, Kai Zhang, and Balwinder Singh.
### Significant improvements from Michael Deakin and Chris Golaz.
###

###===================================================================
### THINGS USERS USUALLY CHANGE (SEE END OF SECTION FOR GUIDANCE)
###===================================================================

### BASIC INFO ABOUT RUN
set job_name       = A_WCYCL1850_template
set compset        = A_WCYCL1850
set resolution     = ne4_oQU240
set machine        = default
set walltime       = default
setenv project       default

### SOURCE CODE OPTIONS
set fetch_code     = false
set e3sm_tag       = master
set tag_name       = default

### CUSTOM CASE_NAME
set case_name = ${machine}.${tag_name}.${job_name}.${resolution}

### BUILD OPTIONS
set debug_compile  = false
set old_executable = false

### AUTOMATIC DELETION OPTIONS
set seconds_before_delete_source_dir = -1
set seconds_before_delete_case_dir   = 10
set seconds_before_delete_bld_dir    = -1
set seconds_before_delete_run_dir    = -1

### SUBMIT OPTIONS
set submit_run       = true
set debug_queue      = true

### PROCESSOR CONFIGURATION
set processor_config = S

### STARTUP TYPE
set model_start_type = initial
set restart_files_dir = none

### DIRECTORIES
set code_root_dir               = default
set e3sm_simulations_dir        = default
set case_build_dir              = default
set case_run_dir                = default
set short_term_archive_root_dir = default

### LENGTH OF SIMULATION, RESTARTS, AND ARCHIVING
set stop_units                  = ndays
set stop_num                    = 1
set restart_units               = $stop_units
set restart_num                 = $stop_num
set num_resubmits               = 0
set do_short_term_archiving     = false

### SIMULATION OPTIONS
set atm_output_freq             = -24
set records_per_atm_output_file = 40
set start_date                  = default

### COUPLER HISTORY FILES
set do_cpl_hist    = true
set cpl_hist_units = ndays
set cpl_hist_num   = 1

#==============================
#EXPLANATION FOR OPTIONS ABOVE:
#==============================

### BASIC INFO ABOUT RUN (1)

#job_name: This is only used to name the job in the batch system. The problem is that batch systems often only
#    provide the first few letters of the job name when reporting on jobs in the queue, which may not be enough
#    to distinguish simultaneous jobs.
#compset: indicates which model components and forcings to use. List choices by typing `create_newcase -list compsets`.
#    An (outdated?) list of options is available at http://www.cesm.ucar.edu/models/cesm1.0/cesm/cesm_doc_1_0_4/a3170.html
#resolution: Model resolution to use. Type `create_newcase -list grids` for a list of options or see
#    http://www.cesm.ucar.edu/models/cesm1.0/cesm/cesm_doc_1_0_4/a3714.htm. Note that E3SM always uses ne30 or ne120 in the atmosphere.
#machine: what machine you are going to run on. This should be 'default' for most machines, and only changed for machines with multiple queues, such as cori.
#    Valid machine names can also be listed via `create_newcase -list machines`
#walltime: How long to reserve the nodes for. The format is HH:MM(:SS); ie 00:10 -> 10 minutes.
#    Setting this to 'default' has the script determine a reasonable value for most runs.
#project: what bank to charge for your run time. May not be needed on some machines.
#    Setting this to 'default' has CIME determine what project to use
#    NOTE: project must be an *environment* variable on some systems.

### SOURCE CODE OPTIONS (2)

#fetch_code: If True, downloads code from github. If False, code is assumed to exist already.
#    NOTE: it is assumed that you have access to the E3SM git repository.  To get access, see:
#    https://acme-climate.atlassian.net/wiki/display/Docs/Installing+the+ACME+Model
#e3sm_tag: E3SM tagname in github. Can be 'master', a git hash, a tag, or a branch name. Only used if fetch_code=True.
#    NOTE: If e3sm_tag=master or master_detached, then this script will provide the latest master version, but detach from the head,
#          to minimize the risk of a user pushing something to master.
#tag_name: Short name for the E3SM branch used. If fetch_code=True, this is a shorter replacement for e3sm_tag
#    (which could be a very long hash!). Otherwise, this is a short name for the branch used. You can
#    choose TAG_NAME to be whatever you want.

### BUILD OPTIONS (3)

#debug_compile: If TRUE, then compile with debug flags
#    Compiling in debug mode will stop the run at the actual location an error occurs, and provide more helpful output.
#    However, it runs about 10 times slower, and is not bit-for-bit the same because some optimizations make tiny change to the
#    numerics.
#old_executable: If this is a path to an executable, then it is used instead of recompiling (it is copied across).
#    If TRUE then skip the build step entirely.
#    If FALSE then build a new executable (using any already compiled files). If you want a clean build then
#    set seconds_before_delete_bld_dir>=0.
#    NOTE: The executable that will be copied should be the same as would be created by compiling (for provenance).
#    NOTE: The path should either be an absolute path, or a path relative to the case_scripts directory.
#    NOTE: old_executable=true is a risk to provenance, so this feature may be removed in the future.
#          However the build system currently rebuilds a few files every time which takes several minutes.
#          When this gets fixed the cost of deleting this feature will be minor.


### AUTOMATIC DELETION OPTIONS (4)

#seconds_before_delete_source_dir : If seconds_before_delete_source_dir>=0 and fetch_code=true, this script automatically deletes
#    the old source code directory after waiting seconds_before_delete_source_dir seconds (to give you the opportunity to cancel
#    by pressing ctrl-C). To turn off the deletion (default behavior), set $num_seconds_before_deleting_source_dir to be negative.
#seconds_before_delete_case_dir : Similar to above, but remove the old case_scripts directory. Since create_newcase dies whenever
#    the case_scripts directory exists, it only makes sense to use $seconds_before_delete_case_dir<0 if you want to be extra careful and
#    only delete the case_scripts directory manually.
#seconds_before_delete_bld_dir : As above, but the old bld directory will be deleted.  This makes for a clean start.
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
#    and the E3SM run_type can be 'initial', 'hybrid', or 'restart', as specified by this script below.
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
#num_resubmits: After a batch job finishes successfully, a new batch job will automatically be submitted to
#    continue the simulation.  num_resubmits is the number of times to submit after initial completion.
#    After the first submission, the CONTINUE_RUN flage in env_run.xml will be changed to TRUE.
#do_short_term_archiving: If TRUE, then move simulation output to the archive directory in your scratch directory.

### SIMULATION OPTIONS (10)

#atm_output_freq (the namelist variable is nhtfrq) : The frequency with which the atmosphere writes its output.
#    0=monthly, +N=every N timesteps,  -N=every N hours
#    For more information:   http://www.cesm.ucar.edu/models/atm-cam/docs/usersguide/node45.html
#records_per_atm_output_file (the namelist variable is mfilt):  The number of time records in each netCDF output file
#    from the atmosphere model. If atm_output_freq=0, then there will only be one time record per file.
#NOTE: If there will be more than one 'history tape' defined in the atm namelist, then
#    atm_output_freq and records_per_atm_output_file should be a comma-separated list of numbers
#    (specifying the option for each history tape).
#start_date: The day that the simulation starts

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
# VERSION OF THIS SCRIPT
#===========================================
set script_ver = 3.0.22

#===========================================
# DEFINE ALIASES
#===========================================

alias lowercase "echo \!:1 | tr '[A-Z]' '[a-z]'"  #make function which lowercases any strings passed to it.
alias uppercase "echo \!:1 | tr '[a-z]' '[A-Z]'"  #make function which uppercases any strings passed to it.

alias e3sm_print 'echo run_e3sm: \!*'
alias e3sm_newline "echo ''"

#===========================================
# ALERT THE USER IF THEY TRY TO PASS ARGUMENTS
#===========================================
set first_argument = $1
if ( $first_argument != '' ) then
 echo 'This script does everything needed to configure/compile/run a case. As such, it'
 echo 'provides complete provenance for each run and makes sharing simulation configurations easy.'
 echo 'Users should make sure that everything required for a run is in this script, the E3SM'
 echo 'git repo, and/or the inputdata svn repo.'
 echo '** This script accepts no arguments. Please edit the script as needed and resubmit without arguments. **'
 exit 5
endif

e3sm_newline
e3sm_print '++++++++ run_e3sm starting ('`date`'), version '$script_ver' ++++++++'
e3sm_newline

#===========================================
# DETERMINE THE LOCATION AND NAME OF THE SCRIPT
#===========================================

# NOTE: CIME 5 and git commands are not cwd agnostic, so compute the absolute paths, then cd to the directories as needed
set this_script_name = `basename $0`
set relative_dir = `dirname $0`
set this_script_dir = `cd $relative_dir ; pwd -P`
set this_script_path = $this_script_dir/$this_script_name

#===========================================
# SETUP DEFAULTS
#===========================================

if ( `lowercase $code_root_dir` == default ) then
  set code_root_dir = `cd $this_script_dir/..; pwd -P`
endif

if ( `lowercase $tag_name` == default ) then
  set pwd_temp       = `pwd`
  set tag_name       = ${pwd_temp:t}
  e3sm_print '$tag_name = '$tag_name
endif

#===========================================
# BASIC ERROR CHECKING
#===========================================

set seconds_after_warning = 10

if ( `lowercase $old_executable` == true ) then
  if ( $seconds_before_delete_source_dir >= 0 ) then
    e3sm_newline
    e3sm_print 'ERROR: It is unlikely that you want to delete the source code and then use the existing compiled executable.'
    e3sm_print '       Hence, this script will abort to avoid making a mistake.'
    e3sm_print '       $seconds_before_delete_source_dir = '$seconds_before_delete_source_dir'      $old_executable = '$old_executable
    exit 11
  endif

  if ( $seconds_before_delete_bld_dir >= 0 ) then
    e3sm_newline
    e3sm_print 'ERROR: It makes no sense to delete the source-compiled code and then use the existing compiled executable.'
    e3sm_print '       Hence, this script will abort to avoid making a mistake.'
    e3sm_print '       $seconds_before_delete_bld_dir = '$seconds_before_delete_bld_dir'      $old_executable = '$old_executable
    exit 12
  endif
endif

if ( `lowercase $case_build_dir` == default && $seconds_before_delete_bld_dir >= 0 ) then
  e3sm_print 'ERROR: run_e3sm cannot delete the build directory when CIME chooses it'
  e3sm_print '       To remedy this, either set $case_build_dir to the path of the executables or disable deleting the directory'
  exit 14
endif

if ( `lowercase $case_run_dir` == default && $seconds_before_delete_run_dir >= 0 ) then
  e3sm_print 'ERROR: run_e3sm cannot delete the run directory when CIME chooses it'
  e3sm_print '       To remedy this, either set $case_run_dir to the path of the executables or disable deleting the directory'
  exit 15
endif

if ( `lowercase $debug_queue` == true && ( $num_resubmits >= 1 || `lowercase $do_short_term_archiving` == true ) ) then
  e3sm_print 'ERROR: Supercomputer centers generally do not allow job chaining in debug queues'
  e3sm_print '       You should either use a different queue, or submit a single job without archiving.'
  e3sm_print '       $debug_queue             = '$debug_queue
  e3sm_print '       $num_resubmits           = '$num_resubmits
  e3sm_print '       $do_short_term_archiving = '$do_short_term_archiving
  exit 16
endif

if ( $restart_num != 0 ) then
  @ remaining_periods = $stop_num - ( $stop_num / $restart_num ) * $restart_num
  if ( $num_resubmits >= 1 && ( $stop_units != $restart_units || $remaining_periods != 0 ) ) then
    e3sm_print 'WARNING: run length is not divisible by the restart write frequency, or the units differ.'
    e3sm_print 'If restart write frequency doesnt evenly divide the run length, restarts will simulate the same time period multiple times.'
    e3sm_print '         $stop_units        = '$stop_units
    e3sm_print '         $stop_num          = '$stop_num
    e3sm_print '         $restart_units     = '$restart_units
    e3sm_print '         $restart_num       = '$restart_num
    e3sm_print '         $remaining_periods = '$remaining_periods
    e3sm_print '         $num_resubmits     = '$num_resubmits
    sleep $seconds_after_warning
  endif
endif

#===========================================
# DOWNLOAD SOURCE CODE IF NEEDED:
#===========================================

### NOTE: you must be setup with access to the E3SM repository before you can clone the repository. For access, see
###       https://acme-climate.atlassian.net/wiki/display/Docs/Installing+the+ACME+Model

if ( `lowercase $fetch_code` == true ) then
  e3sm_print 'Downloading code from the E3SM git repository.'
  if ( -d $code_root_dir/$tag_name ) then
    if ( $seconds_before_delete_source_dir >= 0 ) then
      set num_seconds_until_delete = $seconds_before_delete_source_dir
      e3sm_print 'Removing old code directory '$code_root_dir/$tag_name' in '$num_seconds_until_delete' seconds.'
      e3sm_print 'To abort, press ctrl-C'
      while ( ${num_seconds_until_delete} > 0 )
         e3sm_print ' '${num_seconds_until_delete}'  seconds until deletion.'
         sleep 1
         @ num_seconds_until_delete = ${num_seconds_until_delete} - 1
      end
      rm -fr $code_root_dir/$tag_name
      e3sm_print 'Deleted '$code_root_dir/$tag_name
    else
      e3sm_print 'ERROR: Your branch tag already exists, so dying instead of overwriting.'
      e3sm_print '          You likely want to either set fetch_code=false, change $tag_name, or'
      e3sm_print '          change seconds_before_delete_source_dir.'
      e3sm_print '          Note: $fetch_code = '$fetch_code
      e3sm_print '                $code_root_dir/$tag_name = '$code_root_dir/$tag_name
      e3sm_print '                $seconds_before_delete_source_dir = '$seconds_before_delete_source_dir
      exit 20
    endif #$seconds_before_delete_source_dir >=0
  endif #$code_root_dir exists

  e3sm_print 'Cloning repository into $tag_name = '$tag_name'  under $code_root_dir = '$code_root_dir
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

  if ( `lowercase $e3sm_tag` == master ) then
    e3sm_newline
    ##e3sm_print 'Detaching from the master branch to avoid accidental changes to master by user.'
    ##git checkout --detach
    echo 'KLUDGE: git version on anvil (1.7.1) is too old to be able to detach'
    echo 'edison uses git version 1.8.5.6 and it can git checkout --detach'
  else
    e3sm_newline
    e3sm_print 'Checking out branch ${e3sm_tag} = '${e3sm_tag}
    git checkout ${e3sm_tag}
  endif

endif

e3sm_newline
e3sm_print '$case_name        = '$case_name

#============================================
# DETERMINE THE SCRATCH DIRECTORY TO USE
#============================================

if ( $e3sm_simulations_dir == default ) then
  ### NOTE: csh doesn't short-circuit; so we can't check whether $SCRATCH exists or whether it's empty in the same condition
  if ( ! $?SCRATCH ) then
    e3sm_newline
    e3sm_print 'WARNING: Performing science runs while storing run output in your home directory is likely to exceed your quota'
    e3sm_print '         To avoid any issues, set $e3sm_simulations_dir to a scratch filesystem'
    set e3sm_simulations_dir = ${HOME}/E3SM_simulations
  else
    ### Verify that $SCRATCH is not an empty string
    if ( "${SCRATCH}" == "" ) then
      set e3sm_simulations_dir = ${HOME}/E3SM_simulations
      e3sm_newline
      e3sm_print 'WARNING: Performing science runs while storing run output in your home directory is likely to exceed your quota'
      e3sm_print '         To avoid any issues, set $e3sm_simulations_dir to a scratch filesystem'
    else
      set e3sm_simulations_dir = ${SCRATCH}/E3SM_simulations
    endif
  endif
endif

#============================================
# DELETE PREVIOUS DIRECTORIES (IF REQUESTED)
#============================================
### Determine the case_scripts directory
### Remove existing case_scripts directory (so it doesn't have to be done manually every time)
### Note: This script causes create_newcase to generate a temporary directory (part of a workaround to put the case_name into the script names)
###       If something goes wrong, this temporary directory is sometimes left behind, so we need to delete it too.
### Note: To turn off the deletion, set $num_seconds_until_delete to be negative.
###       To delete immediately, set $num_seconds_until_delete to be zero.

set case_scripts_dir = ${e3sm_simulations_dir}/${case_name}/case_scripts

if ( -d $case_scripts_dir ) then
  if ( ${seconds_before_delete_case_dir} >= 0 ) then
    set num_seconds_until_delete = $seconds_before_delete_case_dir
    e3sm_newline
    e3sm_print 'Removing old $case_scripts_dir directory for '${case_name}' in '${num_seconds_until_delete}' seconds.'
    e3sm_print 'To abort, press ctrl-C'
    while ( ${num_seconds_until_delete} > 0 )
      e3sm_print ' '${num_seconds_until_delete}'  seconds until deletion.'
      sleep 1
      @ num_seconds_until_delete = ${num_seconds_until_delete} - 1
    end
    rm -fr $case_scripts_dir
    e3sm_print ' Deleted $case_scripts_dir directory for : '${case_name}
  else
    e3sm_print 'WARNING: $case_scripts_dir='$case_scripts_dir' exists '
    e3sm_print '         and is not being removed because seconds_before_delete_case_dir<0.'
    e3sm_print '         But create_newcase always fails when the case directory exists, so this script will now abort.'
    e3sm_print '         To fix this, either delete the case_scripts directory manually, or change seconds_before_delete_case_dir'
    exit 35
  endif
endif

### Remove existing build directory (to force a clean compile).  This is good for a new run, but not usually necessary while developing.

if ( `lowercase $case_build_dir` != default && -d $case_build_dir ) then
  if ( ${seconds_before_delete_bld_dir} >= 0 ) then
    set num_seconds_until_delete = $seconds_before_delete_bld_dir
    e3sm_newline
    e3sm_print 'Removing old $case_build_dir directory for '${case_name}' in '${num_seconds_until_delete}' seconds.'
    e3sm_print 'To abort, press ctrl-C'
    while ( ${num_seconds_until_delete} > 0 )
      e3sm_print ' '${num_seconds_until_delete}'  seconds until deletion.'
      sleep 1
      @ num_seconds_until_delete = ${num_seconds_until_delete} - 1
    end
    rm -fr $case_build_dir
    e3sm_print ' Deleted $case_build_dir directory for '${case_name}
  else
    e3sm_print 'NOTE: $case_build_dir='$case_build_dir' exists '
    e3sm_print '      and is not being removed because seconds_before_delete_bld_dir<0.'
  endif
endif

### Remove existing run directory (for a clean start).  This is good for a new run, but often not usually necessary while developing.

if ( `lowercase $case_run_dir` != default &&  -d $case_run_dir ) then
  if ( ${seconds_before_delete_run_dir} >= 0 ) then
    set num_seconds_until_delete = $seconds_before_delete_run_dir
    e3sm_newline
    e3sm_print 'Removing old $case_run_dir directory for '${case_name}' in '${num_seconds_until_delete}' seconds.'
    e3sm_print 'To abort, press ctrl-C'
    while ( ${num_seconds_until_delete} > 0 )
     e3sm_print ' '${num_seconds_until_delete}'  seconds until deletion.'
     sleep 1
     @ num_seconds_until_delete = ${num_seconds_until_delete} - 1
    end
    rm -fr $case_run_dir
    e3sm_print ' Deleted $case_run_dir directory for '${case_name}
  else
    e3sm_print 'NOTE: $case_run_dir='$case_run_dir' exists '
    e3sm_print '      and is not being removed because seconds_before_delete_run_dir<0.'
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
    e3sm_print 'ERROR: $processor_config='$processor_config' is not recognized'
    exit 40
    breaksw
endsw

#================================================================================
# MAKE FILES AND DIRECTORIES CREATED BY THIS SCRIPT READABLE BY EVERYONE IN GROUP
#================================================================================
# Note:  we also want to change the group id for the files.
#        But this can only be done once the run_root_dir has been created,
#        so it is done later.
umask 022

set cime_dir = ${code_root_dir}/${tag_name}/cime
set create_newcase_exe = $cime_dir/scripts/create_newcase
if ( -f ${create_newcase_exe} ) then
  set e3sm_exe = e3sm.exe
  set case_setup_exe = $case_scripts_dir/case.setup
  set case_build_exe = $case_scripts_dir/case.build
  set case_run_exe = $case_scripts_dir/.case.run
  set case_submit_exe = $case_scripts_dir/case.submit
  set xmlchange_exe = $case_scripts_dir/xmlchange
  set xmlquery_exe = $case_scripts_dir/xmlquery
  set shortterm_archive_script = $case_scripts_dir/case.st_archive
  set preview_namelists_exe = $case_scripts_dir/preview_namelists
else                                                                   # No version of create_newcase found
  e3sm_print 'ERROR: ${create_newcase_exe} not found'
  e3sm_print '       This is most likely because fetch_code should be true.'
  e3sm_print '       At the moment, $fetch_code = '$fetch_code
  exit 45
endif

#=============================================================
# DETERMINE THE OPTIONS FOR CREATE_NEWCASE
#=============================================================

set configure_options = "--case ${case_name} --compset ${compset} --script-root ${case_scripts_dir} --res ${resolution} --pecount ${std_proc_configuration} --handle-preexisting-dirs u"

if ( `lowercase $machine` != default ) then
  set configure_options = "${configure_options} --mach ${machine}"
endif

if ( `lowercase $case_build_dir` == default ) then
  set case_build_dir = ${e3sm_simulations_dir}/${case_name}/build
endif

if ( `lowercase $case_run_dir` == default ) then
  set case_run_dir = ${e3sm_simulations_dir}/${case_name}/run
endif

# Default group and permissions on NERSC can be a pain for sharing data
# within the project. This should take care of it. Create top level 
# directory and set the default group to 'acme', permissions for 
# group read access for top level and all files underneath (Chris Golaz).
if ( $machine == 'cori*' || $machine == 'edison' ) then
  mkdir -p ${e3sm_simulations_dir}/${case_name}
  cd ${e3sm_simulations_dir}
  chgrp acme ${case_name}
  chmod 750 ${case_name}
  chmod g+s ${case_name}
  setfacl -d -m g::rx ${case_name}
endif

mkdir -p ${case_build_dir}
set build_root = `cd ${case_build_dir}/..; pwd -P`
mkdir -p ${case_run_dir}
set run_root = `cd ${case_run_dir}/..; pwd -P`

if ( ${build_root} == ${run_root} ) then
  set configure_options = "${configure_options} --output-root ${build_root}"
endif

if ( `lowercase $project` == default ) then
  unsetenv project
else
  set configure_options = "${configure_options} --project ${project}"
endif

#=============================================================
# CREATE CASE_SCRIPTS DIRECTORY AND POPULATE WITH NEEDED FILES
#=============================================================

e3sm_newline
e3sm_print '-------- Starting create_newcase --------'
e3sm_newline

e3sm_print ${create_newcase_exe} ${configure_options}
${create_newcase_exe} ${configure_options}

cd ${case_scripts_dir}

e3sm_newline
e3sm_print '-------- Finished create_newcase --------'
e3sm_newline

#================================================
# UPDATE VARIABLES WHICH REQUIRE A CASE TO BE SET
#================================================

if ( `lowercase $machine` == default ) then
  set machine = `$xmlquery_exe MACH --value`
endif
# machine is a commonly used variable; so make certain it's lowercase
set machine = `lowercase $machine`

if ( `lowercase $case_build_dir` == default ) then
  set case_build_dir = ${e3sm_simulations_dir}/${case_name}/bld
endif
${xmlchange_exe} EXEROOT=${case_build_dir}

if ( `lowercase $case_run_dir` == default ) then
  set case_run_dir = ${case_scripts_dir}/${case_name}/run
endif
${xmlchange_exe} RUNDIR=${case_run_dir}

if ( ! $?project ) then
  setenv project `$xmlquery_exe PROJECT --subgroup case.run --value`
else
  if ( $project == "" ) then
    setenv project `$xmlquery_exe PROJECT --subgroup case.run --value`
  endif
endif

e3sm_print "Project used for submission: ${project}"

#================================
# SET WALLTIME FOR CREATE_NEWCASE
#================================

if ( `lowercase $walltime` == default ) then
  if ( `lowercase $debug_queue` == true ) then
    set walltime = '0:30:00'
  else
    if ( $machine == 'cab' || $machine == 'syrah' ) then
      set walltime = '1:29:00'
    else
      set walltime = '2:00:00'
    endif
  endif
endif

# Allow the user to specify how long the job taks
$xmlchange_exe JOB_WALLCLOCK_TIME=$walltime

#================================
# SET THE STARTDATE FOR THE SIMULATION
#================================

if ( `lowercase $start_date` != 'default' ) then
  $xmlchange_exe RUN_STARTDATE=$start_date
endif

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
#cp /global/u1/p/petercal/junk/streams.ocean $case_scripts_dir/SourceMods/src.mpaso/

#============================================================
# COPY THIS SCRIPT TO THE CASE DIRECTORY TO ENSURE PROVENANCE
#============================================================

set script_provenance_dir  = $case_scripts_dir/run_script_provenance
set script_provenance_name = $this_script_name.`date +%F_%T_%Z`
mkdir -p $script_provenance_dir
cp -f $this_script_path $script_provenance_dir/$script_provenance_name

#=============================================
# CUSTOMIZE PROCESSOR CONFIGURATION
# ============================================
#NOTE: Changes to the processor configuration should be done by an expert.  \
#      Not all possible options will work.

if ( `lowercase $processor_config` == '1' ) then

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

else if ( `lowercase $processor_config` == 'customknl' ) then

  e3sm_print 'using custom layout for cori-knl because $processor_config = '$processor_config

  ${xmlchange_exe} MAX_TASKS_PER_NODE="64"
  ${xmlchange_exe} PES_PER_NODE="256"

  ${xmlchange_exe} NTASKS_ATM="5400"
  ${xmlchange_exe} ROOTPE_ATM="0"

  ${xmlchange_exe} NTASKS_LND="320"
  ${xmlchange_exe} ROOTPE_LND="5120"

  ${xmlchange_exe} NTASKS_ICE="5120"
  ${xmlchange_exe} ROOTPE_ICE="0"

  ${xmlchange_exe} NTASKS_OCN="3840"
  ${xmlchange_exe} ROOTPE_OCN="5440"

  ${xmlchange_exe} NTASKS_CPL="5120"
  ${xmlchange_exe} ROOTPE_CPL="0"

  ${xmlchange_exe} NTASKS_GLC="320"
  ${xmlchange_exe} ROOTPE_GLC="5120"

  ${xmlchange_exe} NTASKS_ROF="320"
  ${xmlchange_exe} ROOTPE_ROF="5120"

  ${xmlchange_exe} NTASKS_WAV="5120"
  ${xmlchange_exe} ROOTPE_WAV="0"

  ${xmlchange_exe} NTHRDS_ATM="1"
  ${xmlchange_exe} NTHRDS_LND="1"
  ${xmlchange_exe} NTHRDS_ICE="1"
  ${xmlchange_exe} NTHRDS_OCN="1"
  ${xmlchange_exe} NTHRDS_CPL="1"
  ${xmlchange_exe} NTHRDS_GLC="1"
  ${xmlchange_exe} NTHRDS_ROF="1"
  ${xmlchange_exe} NTHRDS_WAV="1"

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

# NOTE: For information on the E3SM input data repository, see:
#       https://acme-climate.atlassian.net/wiki/display/WORKFLOW/ACME+Input+Data+Repository

#set input_data_dir = 'input_data_dir_NOT_SET'
#if ( $machine == 'cori*' || $machine == 'edison' ) then
#  set input_data_dir = '/project/projectdirs/m411/ACME_inputdata'    # PJC-NERSC
## set input_data_dir = '/project/projectdirs/ccsm1/inputdata'        # NERSC
#else if ( $machine == 'titan' || $machine == 'eos' ) then
#  set input_data_dir = '/lustre/atlas/proj-shared/cli112/pjcs/ACME_inputdata'    # PJC-OLCF
#endif
#if ( -d  $input_data_dir ) then
#  $xmlchange_exe --id DIN_LOC_ROOT --val $input_data_dir
#else
#  echo 'run_e3sm ERROR: User specified input data directory does NOT exist.'
#  echo '                $input_data_dir = '$input_data_dir
#  exit 270
#endif

### The following command extracts and stores the input_data_dir in case it is needed for user edits to the namelist later.
### NOTE: The following line may be necessary if the $input_data_dir is not set above, and hence defaults to the E3SM default.
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
#NOTE: The xmlchange will fail if CAM is not active, so test whether a data atmosphere (datm) is used.

if ( `$xmlquery_exe --value COMP_ATM` == 'datm'  ) then 
  e3sm_newline
  e3sm_print 'The specified configuration uses a data atmosphere, so cannot activate COSP simulator.'
  e3sm_newline
else
  e3sm_newline
  e3sm_print 'Configuring E3SM to use the COSP simulator.'
  e3sm_newline
  $xmlchange_exe --id CAM_CONFIG_OPTS --append --val='-cosp'
endif

#===========================
# SET THE PARTITION OF NODES
#===========================

if ( `lowercase $debug_queue` == true ) then
  if ( $machine == cab || $machine == sierra ) then
    $xmlchange_exe --id JOB_QUEUE --val 'pdebug'
  else if ($machine != sandiatoss3 && $machine != bebop && $machine != blues) then
    $xmlchange_exe --id JOB_QUEUE --val 'debug'
  endif
endif

#============================================
# CONFIGURE
#============================================

#note configure -case turned into cesm_setup in cam5.2

e3sm_newline
e3sm_print '-------- Starting case.setup --------'
e3sm_newline

e3sm_print ${case_setup_exe}

${case_setup_exe} --reset

e3sm_newline
e3sm_print '-------- Finished case.setup  --------'
e3sm_newline

#============================================================
#MAKE GROUP PERMISSIONS to $PROJECT FOR THIS DIRECTORY
#============================================================
#this stuff, combined with the umask command above, makes
#stuff in $run_root_dir readable by everyone in e3sm group.

set run_root_dir = `cd $case_run_dir/..; pwd -P`

#both run_root_dir and case_scripts_dir are created by create_newcase,
#so run_root_dir group isn't inherited for case_scripts_dir

#========================================================
# CREATE LOGICAL LINKS BETWEEN RUN_ROOT & THIS_SCRIPT_DIR
#========================================================

#NOTE: This is to make it easy for the user to cd to the case directory
#NOTE: Starting the suffix wit 'a' helps to keep this near the script in ls
#      (but in practice the behavior depends on the LC_COLLATE system variable).

e3sm_print 'Creating logical links to make navigating easier.'
e3sm_print 'Note: Beware of using ".." with the links, since the behavior of shell commands can vary.'

# Customizations from Chris Golaz
# Link in this_script_dir case_dir
set run_dir_link = $this_script_dir/$this_script_name=a_run_link

e3sm_print ${run_dir_link}

if ( -l $run_dir_link ) then
  rm -f $run_dir_link
endif
e3sm_print "run_root ${run_root_dir}"
e3sm_print "run_dir ${run_dir_link}"

ln -s $run_root_dir $run_dir_link

#============================================
# SET BUILD OPTIONS
#============================================

if ( `uppercase $debug_compile` != 'TRUE' && `uppercase $debug_compile` != 'FALSE' ) then
  e3sm_print 'ERROR: $debug_compile can be true or false but is instead '$debug_compile
  exit 220
endif

if ( $machine == 'edison' && `uppercase $debug_compile` == 'TRUE' ) then
  e3sm_print 'ERROR: Edison currently has a compiler bug and crashes when compiling in debug mode (Nov 2015)'
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

cat <<EOF >> user_nl_clm
! finidat=''
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

if ( `lowercase $old_executable` == false ) then
  e3sm_newline
  e3sm_print '-------- Starting Build --------'
  e3sm_newline

  e3sm_print ${case_build_exe}
  ${case_build_exe}

  e3sm_newline
  e3sm_print '-------- Finished Build --------'
  e3sm_newline
else if ( `lowercase $old_executable` == true ) then
  if ( -x $case_build_dir/$e3sm_exe ) then       #use executable previously generated for this case_name.
    e3sm_print 'Skipping build because $old_executable='$old_executable
    e3sm_newline
    #create_newcase sets BUILD_COMPLETE to FALSE. By using an old executable you're certifying
    #that you're sure the old executable is consistent with the new case... so be sure you're right!
    #NOTE: This is a risk to provenance, so this feature may be removed in the future [PJC].
    #      However the build system currently rebuilds several files every time which takes many minutes.
    #      When this gets fixed the cost of deleting this feature will be minor.
    #      (Also see comments for user options at top of this file.)
    e3sm_print 'WARNING: Setting BUILD_COMPLETE = TRUE.  This is a little risky, but trusting the user.'
    $xmlchange_exe --id BUILD_COMPLETE --val TRUE
  else
    e3sm_print 'ERROR: $old_executable='$old_executable' but no executable exists for this case.'
    e3sm_print '                Expected to find executable = '$case_build_dir/$e3sm_exe
    exit 297
  endif
else
  if ( -x $old_executable ) then #if absolute pathname exists and is executable.
    #create_newcase sets BUILD_COMPLETE to FALSE. By copying in an old executable you're certifying
    #that you're sure the old executable is consistent with the new case... so be sure you're right!
    #NOTE: This is a risk to provenance, so this feature may be removed in the future [PJC].
    #      However the build system currently rebuilds several files every time which takes many minutes.
    #      When this gets fixed the cost of deleting this feature will be minor.
    #      (Also see comments for user options at top of this file.)
    #
    #NOTE: The alternative solution is to set EXEROOT in env_build.xml.
    #      That is cleaner and quicker, but it means that the executable is outside this directory,
    #      which weakens provenance if this directory is captured for provenance.
    e3sm_print 'WARNING: Setting BUILD_COMPLETE = TRUE.  This is a little risky, but trusting the user.'
    $xmlchange_exe --id BUILD_COMPLETE --val TRUE
    cp -fp $old_executable $case_build_dir/
  else
    e3sm_print 'ERROR: $old_executable='$old_executable' does not exist or is not an executable file.'
    exit 297
  endif
endif

# Some user_nl settings won't be updated to *_in files under the run directory
# until namelists are built again.
# Call preview_namelists to make sure *_in and user_nl files are consistent.
$preview_namelists_exe

#============================================
# BATCH JOB OPTIONS
#============================================

# Set options for batch scripts (see above for queue and batch time, which are handled separately)

# NOTE: This also modifies the short-term archiving script.
# NOTE: We want the batch job log to go into a sub-directory of case_scripts (to avoid it getting clogged up)

# NOTE: we are currently not modifying the archiving scripts to run in debug queue when $debug_queue=true
#       under the assumption that if you're debugging you shouldn't be archiving.

# NOTE: there was 1 space between MSUB or PBS and the commands before cime and there are 2 spaces
#       in post-cime versions. This is fixed by " \( \)*" in the lines below. The "*" here means
#       "match zero or more of the expression before". The expression before is a single whitespace.
#       The "\(" and "\)" bit indicate to sed that the whitespace in between is the expression we
#       care about. The space before "\(" makes sure there is at least one whitespace after #MSUB.
#       Taken all together, this stuff matches lines of the form #MSUB<one or more whitespaces>-<stuff>.

mkdir -p batch_output      ### Make directory that stdout and stderr will go into.

if ( $machine =~ 'cori*' || $machine == edison ) then
    $xmlchange_exe --subgroup case.run BATCH_COMMAND_FLAGS="--job-name=${job_name} --output=batch_output/${case_name}.o%j"
    $xmlchange_exe --subgroup case.st_archive BATCH_COMMAND_FLAGS="--job-name=ST+${job_name} --output=batch_output/ST+${case_name}.o%j --account=${project}"
else if ( $machine == titan || $machine == eos ) then
    $xmlchange_exe --subgroup case.run BATCH_COMMAND_FLAGS="-N ${job_name} -A ${project} -o batch_output/${case_name}.o%j"
    $xmlchange_exe --subgroup case.st_archive BATCH_COMMAND_FLAGS="-N ST+${job_name} -o batch_output/ST+${case_name}.o%j"

else if ( $machine == anvil ) then
# Priority for Anvil
# For more information, see
# https://acme-climate.atlassian.net/wiki/pages/viewpage.action?pageId=98992379#Anvil:ACME'sdedicatednodeshostedonBlues.-Settingthejobpriority
# If default; use the default priority of qsub.py. Otherwise, should be from 0-5, where 0 is the highest priority
# Note that only one user at a time is allowed to use the highest (0) priority
# env_batch.xml must be modified by hand, as it doesn't conform to the entry-id format
# ${xmlchange_exe} batch_submit="qsub.py "
    set anvil_priority   = default
    if ( `lowercase ${anvil_priority}` != default ) then
	sed -i 's:qsub:qsub.py:g' env_batch.xml
	$xmlchange_exe BATCH_COMMAND_FLAGS="-W x=QOS:pri${anvil_priority}"
    endif

else
    e3sm_print 'WARNING: This script does not have batch directives for $machine='$machine
    e3sm_print '         Assuming default E3SM values.'
endif

#============================================
# QUEUE OPTIONS
#============================================
# Edit the default queue and batch job lengths.

# HINT: To change queue after run submitted, the following works on most machines:
#       qalter -lwalltime=00:29:00 <run_descriptor>
#       qalter -W queue=debug <run_descriptor>

### Only specially authorized people can use the special_e3sm qos on Cori or Edison. Don't uncomment unless you're one.
#if ( `lowercase $debug_queue` == false && $machine == edison ) then
#  set batch_options = `$xmlquery_exe BATCH_COMMAND_FLAGS --value`
#  $xmlchange_exe BATCH_COMMAND_FLAGS="${batch_options} --qos=special_e3sm"
#endif

#============================================
# SETUP SHORT TERM ARCHIVING
#============================================

$xmlchange_exe --id DOUT_S --val `uppercase $do_short_term_archiving`
if ( `lowercase $short_term_archive_root_dir` != default ) then
  $xmlchange_exe --id DOUT_S_ROOT --val $short_term_archive_root_dir
endif

set short_term_archive_root_dir = `$xmlquery_exe DOUT_S_ROOT --value`

#==============================
# SETUP PERMISSIONS FOR SHARING
#==============================

set group_list = `groups`
if ( "$group_list" =~ "*${project}*" ) then
  # Determine what command to use to setup permissions
  where setfacl > /dev/null
  if ( $? == 0 ) then
    # setfacl exists, but may not work depending on kernel configuration
    # So, verify it works, and if not, just use chgrp
    set group_perms = "setfacl -Rdm g:${project}:r-x"
    e3sm_print ${group_perms} ${case_run_dir}
    ${group_perms} ${case_run_dir}
    if ( $? != 0 ) then
      set group_perms = "chgrp ${project}"
    endif
    ${group_perms} ${case_run_dir}
    # Ensure chgrp works also, in case there's something we didn't anticipate
    if ( $? != 0 ) then
      unset group_perms
      e3sm_print "Could not make results accessible to the group, results must be shared manually"
    endif
  endif

  if ( $?group_perms ) then
    # Make results which have been archived accessible to other project members
    if (! -d ${short_term_archive_root_dir} )then
      mkdir -p ${short_term_archive_root_dir}
    endif
    e3sm_print ${group_perms} ${short_term_archive_root_dir}
    ${group_perms} ${short_term_archive_root_dir}

    e3sm_print "All project members have been given access to the results of this simulation"
  endif
else
  e3sm_print "${project} not recognized as a group, results must be shared manually"
endif

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

e3sm_newline
e3sm_print '$model_start_type = '${model_start_type}'  (This is NOT necessarily related to RUN_TYPE)'

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
    e3sm_print 'ERROR: ${rpointer_filename} does not exist. It is needed to extract RUN_REFDATE.'
    e3sm_print "       This may be because you should set model_start_type to 'initial' or 'continue' rather than 'branch'."
    e3sm_print '       ${rpointer_filename} = '{rpointer_filename}
    exit 370
  endif
  set restart_coupler_filename = `cat $rpointer_filename`
  set restart_case_name = ${restart_coupler_filename:r:r:r:r}         # Extract out the case name for the restart files.
  set restart_filedate = ${restart_coupler_filename:r:e:s/-00000//}   # Extract out the date (yyyy-mm-dd).
  e3sm_print '$restart_case_name = '$restart_case_name
  e3sm_print '$restart_filedate  = '$restart_filedate

  ### the next line gets the YYYY-MM of the month before the restart time. Needed for staging history files.
  ### NOTE: This is broken for cases that have run for less than a month
  set restart_prevdate = `date -d "${restart_filedate} - 1 month" +%Y-%m`

  e3sm_print '$restart_prevdate  = '$restart_prevdate

  e3sm_print 'Copying stuff for branch run'
  cp -s ${restart_files_dir}/${restart_case_name}.cam.r.${restart_filedate}-00000.nc $case_run_dir
  cp -s ${restart_files_dir}/${restart_case_name}.cam.rs.${restart_filedate}-00000.nc $case_run_dir
  cp -s ${restart_files_dir}/${restart_case_name}.clm2.r.${restart_filedate}-00000.nc $case_run_dir
  cp -s ${restart_files_dir}/${restart_case_name}.clm2.rh0.${restart_filedate}-00000.nc $case_run_dir
  cp -s ${restart_files_dir}/${restart_case_name}.cpl.r.${restart_filedate}-00000.nc $case_run_dir
  cp -s ${restart_files_dir}/${restart_case_name}.mosart.r.${restart_filedate}-00000.nc $case_run_dir
  cp -s ${restart_files_dir}/${restart_case_name}.mosart.rh0.${restart_filedate}-00000.nc $case_run_dir
  cp -s ${restart_files_dir}/mpascice.rst.${restart_filedate}_00000.nc $case_run_dir
  cp -s ${restart_files_dir}/mpaso.rst.${restart_filedate}_00000.nc $case_run_dir
  cp -s ${restart_files_dir}/../../atm/hist/${restart_case_name}.cam.h0.${restart_prevdate}.nc $case_run_dir
  cp -s ${restart_files_dir}/../../rof/hist/${restart_case_name}.mosart.h0.${restart_prevdate}.nc $case_run_dir
  cp -s ${restart_files_dir}/../../lnd/hist/${restart_case_name}.clm2.h0.${restart_prevdate}.nc $case_run_dir
  cp ${restart_files_dir}/rpointer* $case_run_dir

  $xmlchange_exe --id RUN_TYPE --val "branch"
  $xmlchange_exe --id RUN_REFCASE --val $restart_case_name
  $xmlchange_exe --id RUN_REFDATE --val $restart_filedate    # Model date of restart file
  $xmlchange_exe --id CONTINUE_RUN --val "FALSE"
  # Currently broken in CIME
  # Only uncomment this if you want to continue the run with the same name (risky)!!
  # $xmlchange_exe --id BRNCH_RETAIN_CASENAME --val "TRUE"

else

  e3sm_print 'ERROR: $model_start_type = '${model_start_type}' is unrecognized.   Exiting.'
  exit 380

endif

#============================================
# RUN CONFIGURATION OPTIONS
#============================================

#NOTE:  This section is for making specific changes to the run options (ie env_run.xml).

#if ( $machine == 'cori*' ) then      ### fix pnetcdf problem on Cori. (github #593)
#  $xmlchange_exe --id PIO_TYPENAME  --val "netcdf"
#endif

#=================================================
# SUBMIT THE SIMULATION TO THE RUN QUEUE
#=================================================
#note: to run the model in the totalview debugger,
# cd $case_run_dir
# totalview srun -a -n <number of procs> -p <name of debug queue> ../bld/$e3sm_exe
# where you may need to change srun to the appropriate submit command for your system, etc.


e3sm_newline
e3sm_print '-------- Starting Submission to Run Queue --------'
e3sm_newline

if ( ${num_resubmits} > 0 ) then
  ${xmlchange_exe} --id RESUBMIT --val ${num_resubmits}
  e3sm_print 'Setting number of resubmits to be '${num_resubmits}
  @ total_submits = ${num_resubmits} + 1
  e3sm_print 'This job will submit '${total_submits}' times after completion'
endif

if ( `lowercase $submit_run` == 'true' ) then
  e3sm_print '         SUBMITTING JOB:'
  e3sm_print ${case_submit_exe}
  ${case_submit_exe}
else
    e3sm_print 'Run NOT submitted because $submit_run = '$submit_run
endif

e3sm_newline
e3sm_print '-------- Finished Submission to Run Queue --------'
e3sm_newline

#=================================================
# DO POST-SUBMISSION THINGS (IF ANY)
#=================================================

# Actions after the run submission go here.

e3sm_newline
e3sm_print '++++++++ run_e3sm Completed ('`date`') ++++++++'
e3sm_newline

#**********************************************************************************
### --- end of script - there are no commands beyond here, just useful comments ---
#**********************************************************************************

### -------- Version information --------
# 1.0.0    2015-11-19    Initial version.  Tested on Titan. (PJC)
# 1.0.1    2015-11-19    Fixed bugs and added features  for Hopper. (PJC)
# 1.0.2    2015-11-19    Modified to conform with E3SM script standards. PJC)
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
# 1.0.14   2015-11-25    Added custom PE configuration so the E3SM pre-alpha code will work on Titan for Chris Golaz.   (PJC)
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
#                        Checks case_name is 79 characters, or less, which is a requirement of the E3SM scripts.
#                        Improved options for SLURM machines.
#                        Added numbers for the ordering of options at top of file (in preperation for reordering).
#                        Added xxdiff calls to fix known bugs in master> (need to generalize for other people)
# 1.0.32   2016-01-07    Converted inputdata_dir to input_data_dir for consistency.      (PJC)
#                        Cosmetic improvements.
# 1.0.33   2016-01-08    Changed default tag to master_detached to improve clarity. (PJC)
#                        Now sets up E3SM git hooks when fetch_code=true.
# 1.0.33p  2016-01-08    Changed compset from A_B1850CN to A_B1850 (pre-e3sm script only).  (PJC)
#                        Added finidat = '' to user_nl_clm, which allows A_B1850 to run.
# 1.0.34   2016-01-12    Commented out the input_data_dir user configuration, so it defaults to the E3SM settings.   (PJC)
# 1.0.35   2016-01-13    Improved an error message.   (PJC)
# 1.0.36   2016-01-21    Reordered options to better match workflow. (PJC)
# 1.2.0    2016-01-21    Set options to settings for release. (PJC)
# 1.2.1    2016-01-21    Reordered and refined comments to match new ordering of options. (PJC)
# 1.2.2    2016-01-21    The batch submission problem on Cori has been repaired on master (#598),
#                        so I have undone the workaround in this script. (PJC)
# 1.2.3    2016-01-26    Commented out some of the workarounds for E3SM bugs that are no longer needed.  (PJC)
# 1.4.0    2016-03-23    A number of modifications to handle changes in machines and E3SM. [version archived to E3SM] (PJC)
# 1.4.1    2016-03-23    Modified to defaults for Cori (NERSC). (PJC)
# 1.4.2    2016-08-05    Replaced cime_space with pattern matching, added num_depends functionality for daisychained
#                        jobs, added code for submitting to qos=e3sm_special on Edison, added cpl_hist options, and
#                        improved support for sierra and cab at LLNL.(PMC)
# 1.4.3    2016-08-10    Improved support for branch runs (PMC)
# 1.4.4    2016-08-11    Added umask command to make run directory world-readable by default.
# 2.0.0    2016-08-10    Added capability to a chain of submissions using the script auto_chain_runs.$machine (PJC)
# 2.0.1    2016-09-13    Fixed num_resubmits undefined error.
#                        Generalized setting of group permissions for other machines. (PJC)
# 2.0.2    2016-09-13    Turned off short- and long-term archiving so auto_chain_runs script can do it manually. (PJC)
# 2.0.3    2016-09-14    Removed 'git --set-upstream' command, because it does not work on tags.  (PJC)
# 2.0.4    2016-09-14    Long term archiving not working in E3SM, so turn it off and warn user.  (PJC)
# 3.0.0    2016-12-15    Initial update for CIME5. Change script names, don't move the case directory
#                        as it's broken by the update, use xmlchange to set up the custom PE Layout.
#                        Remove support for CIME2 and pre-cime versions.  (MD)
# 3.0.1    2017-01-26    Setup to run A_WCYCL1850S simulation at ne30 resolution.  (CG)
# 3.0.2    2017-02-13    Activated logical links by default, and tweaked the default settings. (PJC)
# 3.0.3    2017-03-24    Added cori-knl support, made walltime an input variable, changed umask to 022, and
#                        deleted run_name since it wasn't being used any more.
# 3.0.4    2017-03-31    Added version to E3SM repository. Working on using more defaults from CIME.
#                        Use 'print' and 'newline' for standardized output (MD)
# 3.0.5    2017-04-07    Restored functionality to delete of run and build directories, and reuse other builds.
#                        Merged in PMC's changes from 3.0.4. Enabled using CIME defaults for more functionality
#                        Renamed 'print' and 'newline' to 'e3sm_print' and 'e3sm_newline'
#                        to disambiguate them from system commands (MD)
# 3.0.6    2017-04-27    Implemented PJC's "hack" in a machine independent way to
#                        restore the run e3sm groups preferred directory structure
#                        Add a warning if the default output directory is in the users home
#                        Give project a default value; if used, CIME will determine the batch account to use
#                        Remove the warning about not running in interactive mode;
#                        use the new CIME option --handle-preexisting-dirs to avoid this potential error
#                        Fix the usage of xmlchange for the customknl configuration
#                        Set walltime to default to get more time on Edison (MD)
# 3.0.7    2017-05-22    Fix for the new CIME 5.3; use the --script-root option instead of PJC's "hack"
#                        Note that this breaks compatibility with older versions of CIME
#                        Also add a fix to reenable using the special e3sm qos queue on Edison (MD)
# 3.0.8    2017-05-24    Fixed minor bug when $machine contained a capital letter. Bug was introduced recently. (PJC)
# 3.0.9    2017-06-19    Fixed branch runs. Also removed sed commands for case.run and use --batch-args in case.submit (MD)
# 3.0.10   2017-06-14    To allow data-atm compsets to work, I added a test for CAM_CONFIG_OPTS. (PJC)
# 3.0.11   2017-07-14    Replace auto-chaining code with E3SM's resubmit feature. Also fix Edison's qos setting (again...) (MD)
# 3.0.12   2017-07-24    Supports setting the queue priority for anvil. Also move making machine lowercase up to clean some things up (MD)
# 3.0.13   2017-08-07    Verify that the number of periods between a restart evenly divides the number until the stop with the same units.
#                        Update the machine check for cori to account for cori-knl (MD)
# 3.0.14   2017-09-11    Add checks for blues and bebop when trying to use the debug queue. Mostly by Andy Salinger with assist from (MD)
# 3.0.15   2017-09-18    Removes long term archiving settings, as they no longer exist in CIME (MD)
# 3.0.16   2017-10-17    Brings in CGs changes to make branch runs faster and easier. Also adds the machine name to case_name
# 3.0.17   2017-10-31    Trivial bug fix for setting cosp (MD)
# 3.0.18   2017-12-07    Update cime script names which have been hidden (MD)
# 3.0.19   2017-12-07    Remove all references to ACME except for online links. Also updates the case.st_archive name again (MD)
# 3.0.20   2018-04-03    Add in a setfacl command for the run and short term archiving directory (MD)
#
# NOTE:  PJC = Philip Cameron-Smith,  PMC = Peter Caldwell, CG = Chris Golaz, MD = Michael Deakin

### ---------- Desired features still to be implemented ------------
# +) fetch_code = update   (pull in latest updates to branch)    (PJC)
# +) A way to run the testsuite.? (PJC)
# +) make the handling of lowercase consistent.  $machine may need to be special. (PJC)
# +) generalize xxdiff commands (for fixing known bugs) to work for other people  (PJC)
# +) Add a 'default' option, for which REST_OPTION='$STOP_OPTION' and REST_N='$STOP_N'.
#    This is important if the user subsequently edits STOP_OPTION or STOP_N.      (PJC)
# +) triggering on $e3sm_tag = master_detached doesn't make sense.  Fix logic. (PJC)
# +) run_root and run_root_dir are duplicative.  Also, move logical link creation before case.setup (PJC)
# +) change comments referring to cesm_setup to case.setup (name has changed). (PJC)

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
