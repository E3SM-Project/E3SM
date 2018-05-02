#!/bin/bash

#==============================================================================
# $Id$
# $URL$
#
# Sets up a 12-month(default) or nine time step run
#
#==============================================================================

#==============================================================================
# Usage subroutine
#==============================================================================

usage() {

  if [ $ThisScript = "ensemble.sh" ]; then
    
      echo "USAGE: $ThisScript --case CASE ---mach MACH [--project PROJECT_NUM] [--walltime WALLTIME] [--compiler COMPILER] [--compset COMPSET] [--res RES] [--nb] [--ns] [--uf] [--ensemble ENS_SIZE]"
      echo ''
      echo 'Sets up CESM cases for either an ensemble of runs or a size 3 test set (default).'\
       'Use pyCECT utilities to create an ensemble summary file or to evaluate the small' \
       'test set of runs against the ensemble.'
 else
      echo "USAGE: $ThisScript --case CASE ---mach MACH [--project PROJECT_NUM] [--pertlim PERTLIM] [--walltime WALLTIME] [--compiler COMPILER] [--compset COMPSET] [--res RES] [--nb] [--ns] [--uf]"
      echo ''
      echo 'Sets up a single CESM case (typically called via ensemble.sh).'
 fi
    echo ''
    echo 'Required flags:'
    echo '  --case <name>    Case name passed on to create_newcase'
    echo '  --mach <name>    Machine name passed on to create_newcase'
   echo ''
    echo 'Optional flags (+ all "--" options to create_newcase): '
    echo '  --project <num>  Project number to charge in job scripts'
  if [ $ThisScript = "single_run.sh" ]; then
    echo '  --pertlim <num>     Run CAM with non-zero pertlim'
  fi
  echo '  --walltime <hr:mn> Amount of walltime requested (default = 4:30, or 0:10 with --uf enabled)'
  echo '  --compiler <name>  Compiler to use (default = same as Machine default) '
  echo '  --compset <name>   Compset to use (default = F2000)'
  echo '  --res <name>       Resolution to run (default = f19_f19)'
  echo '  --uf               Enable ninth time step runs (ultra-fast mode) - otherwise the default is 12-month runs'
  if [ $ThisScript = "ensemble.sh" ]; then
    echo '  --nb               Disables auto building the root case of the ensemble.'
    echo '  --ns               Disables auto submitting any members of the ensemble.'
    echo '  --ensemble <size>  Build the ensemble (instead of building 3 cases with random pertlim values for verification),'
    echo '                     and specify the number of ensemble members to generate (e.g.: 151 for annual averages or 350 for ultra-fast mode)'
  else
    echo '  --nb               Disables building (and submitting) the single case.'
    echo '  --ns               Disables submitting the single case.'
  fi
  echo '  -h, --help         Output this usage information'
}


#==============================================================================
# Main Program
#==============================================================================

ThisScript=$(basename $0)
ThisDir=$(cd `dirname $0`; pwd -P )

# Default Values
RES="f19_f19"
COMPSET="F2000"
WallTime="4:30"
PERTLIM="0"
nobuild="off"
nosubmit="off"
UF=0
WallTimeUser=0

# Process input arguments
while [ $# -gt 0 ]; do
  case $1 in
    --case )
      NewCaseFlags="$NewCaseFlags $1 $2"
      CASE=$2
      CASENAME=$(basename $CASE)
      shift
    ;;
    --mach ) 
      NewCaseFlags="$NewCaseFlags $1 $2"
      MACH=$2
      shift
    ;;
    --pertlim )
      PERTLIM="$2"
      shift
    ;;
     --walltime )
      WallTimeUser=$2
      shift
    ;;
    --project )
      PROJECT="$2"
      NewCaseFlags="$NewCaseFlags $1 $2"
      shift
    ;;
    --compiler )
      NewCaseFlags="$NewCaseFlags $1 $2"
      shift
    ;;
    --compset )
      COMPSET=$2
      shift
    ;;
    --res )
      RES=$2
      shift
    ;;
    --nb )
       nobuild="on"
       nosubmit="on"
    ;;
    --ns )
       nosubmit="on"
    ;;
     --uf )
      UF=1
    ;;
    -h|--help )
      usage
      exit 0
    ;;
    #other flags to pass through for create_newcase	  
    --verbose )
       NewCaseFlags="$NewCaseFlags $1"
    ;;
    --silent )
       NewCaseFlags="$NewCaseFlags $1"
    ;;
    --test )
       NewCaseFlags="$NewCaseFlags $1"
    ;;
    --multi-driver )
       NewCaseFlags="$NewCaseFlags $1"
    ;;
    --pecount )
      NewCaseFlags="$NewCaseFlags $1 $2"
      shift
    ;;
    --nist )
      NewCaseFlags="$NewCaseFlags $1 $2"
      shift
    ;;
    --mpilib )
      NewCaseFlags="$NewCaseFlags $1 $2"
      shift
    ;;
    --pesfile )
      NewCaseFlags="$NewCaseFlags $1 $2"
      shift
    ;;
   --gridfile )
       NewCaseFlags="$NewCaseFlags $1 $2"
       shift
    ;;
   --srcroot )
       NewCaseFlags="$NewCaseFlags $1 $2"
       shift
    ;;
  --output-root )
       NewCaseFlags="$NewCaseFlags $1 $2"
       shift
    ;;
  --script-root )
       NewCaseFlags="$NewCaseFlags $1 $2"
       shift
    ;;
  --queue )
       NewCaseFlags="$NewCaseFlags $1 $2"
       shift
    ;;
  --user-mods-dir )
       NewCaseFlags="$NewCaseFlags $1 $2"
       shift
    ;;
  --input-dir )
       NewCaseFlags="$NewCaseFlags $1 $2"
       shift
    ;;
   # Ignore these flags, they will be used by ensemble.sh
    --ensemble )
      shift
    ;;
    * )
      echo "ERROR: invalid argument $1"
      echo ''
      usage
      exit 1
    ;;
  esac
  shift
done

# SCRIPTS_ROOT is only used for generating ensemble!
SCRIPTS_ROOT=$(cd `dirname $0`; cd ../../scripts; pwd )
echo "STATUS: (Running from: $SCRIPTS_ROOT)"

#if user didn't specify, use default walltime
if [ "$WallTimeUser" == "0" ]; then
    echo "STATUS: user did not specify walltime"
        if [ $UF -eq 1 ]; then
	echo "STATUS: UF run walltime default..."
	WallTime="0:10"
    fi
else
    WallTime=$WallTimeUser 
fi

#machine and case are not optional
if [ -z "$CASE" ]; then
  echo "Must specify --case argument"
  echo "Invoke $ThisScript -h for usage"
  exit 1
fi

if [ -z "$MACH" ]; then
  echo "Must specify --mach argument"
  echo "Invoke $ThisScript -h for usage"
  exit 1
fi

#res and compset are required for create_newcase, so specify here + walltime + --run_unsupported
NewCaseFlags="$NewCaseFlags --res $RES --compset $COMPSET --walltime $WallTime --run-unsupported"


cd $SCRIPTS_ROOT
echo "STATUS: Currently in $SCRIPTS_ROOT"
echo "STATUS: Flags for create_newcase are $NewCaseFlags"
./create_newcase $NewCaseFlags

cd $CASE
chmod u+w *

# Change env_build to ensure bit-for-bit, disable archiving, disable restart
# files, and only run for 1 month
if [ ! -e env_run.xml.orig ]; then
  cp env_run.xml env_run.xml.orig
fi

echo "STATUS: Adjusting env_run.xml...."
./xmlchange --file env_run.xml --id BFBFLAG --val TRUE
./xmlchange --file env_run.xml --id DOUT_S --val FALSE
./xmlchange --file env_run.xml --id REST_OPTION --val never
# Set to time steps if uf selected
if [ $UF -eq 1 ]; then
  ./xmlchange --file env_run.xml --id STOP_OPTION --val nsteps
  ./xmlchange --file env_run.xml --id STOP_N --val 9
else
  ./xmlchange --file env_run.xml --id STOP_OPTION --val nmonths
  ./xmlchange --file env_run.xml --id STOP_N --val 12
fi

echo "STATUS: Running setup....."
./case.setup || { echo "Error running case_setup, probably a bad NP value"; exit 1; }
chmod u+w user_nl_*

echo "STATUS: Adjusting user_nl_* files...."
# Only edit user_nl_cam if the file exists (otherwise not using CAM!)
if [ -e user_nl_cam ]; then
  # Have CAM output everything (single precision)
  # But not initial data...
  # Check if uf runs
  if [ $UF -eq 1 ]; then
    echo "avgflag_pertape = 'I'" >> user_nl_cam
    echo "nhtfrq  = 9" >> user_nl_cam
  else
    echo "avgflag_pertape = 'A'" >> user_nl_cam
    echo "nhtfrq  = -8760" >> user_nl_cam
  fi
  echo "inithist = 'NONE'" >> user_nl_cam
  if [ ! "$PERTLIM" == "0" ]; then
    echo "pertlim = $PERTLIM" >> user_nl_cam
  fi
fi

# Only edit user_nl_clm if the file exists (otherwise not using CLM!)
if [ -e user_nl_clm ]; then
  # Have CLM output everything (single precision)
  if [ $UF -eq 1 ]; then
    echo "hist_avgflag_pertape = 'I'" >> user_nl_clm
    echo "hist_nhtfrq  = 9" >> user_nl_clm
  else
    echo "hist_avgflag_pertape = 'A'" >> user_nl_clm
    echo "hist_nhtfrq  = -8760" >> user_nl_clm
  fi
fi

# Disable output from all other components
if [ -e user_nl_cice ]; then
  echo "histfreq = 'x','x','x','x','x'" >> user_nl_cice
fi
if [ -e user_nl_pop2 ]; then
  echo "n_tavg_streams = 1" >> user_nl_pop2
  echo "ldiag_bsf = .false." >> user_nl_pop2
  echo "ldiag_global_tracer_budgets = .false." >> user_nl_pop2
  echo "ldiag_velocity = .false." >> user_nl_pop2
  echo "diag_gm_bolus = .false." >> user_nl_pop2
  echo "ltavg_nino_diags_requested = .false." >> user_nl_pop2
  echo "moc_requested = .false." >> user_nl_pop2
  echo "n_heat_trans_requested = .false." >> user_nl_pop2
  echo "n_salt_trans_requested = .false." >> user_nl_pop2
  echo "tavg_freq_opt = 'once', 'never', 'never'" >> user_nl_pop2
  echo "tavg_file_freq_opt = 'once', 'never', 'never'" >> user_nl_pop2
  echo "diag_cfl_freq_opt = 'never'" >> user_nl_pop2
  echo "diag_global_freq_opt = 'never'" >> user_nl_pop2
  echo "diag_transp_freq_opt = 'never'" >> user_nl_pop2
  rm -rf SourceMods/src.pop2/gx1v6_tavg_contents
  touch  SourceMods/src.pop2/gx1v6_tavg_contents
fi

echo "STATUS: Updating namelists...."
./preview_namelists || { echo "Error running preview_namelists"; exit 1; }

echo "nobuild: $nobuild, nosubmit: $nosubmit"
# Build executable
if [ $nobuild != 'on' ]; then
    echo "STATUS: Building exe..."
    ./case.build || { echo "Error building!"; exit 1; }
fi

# Submit job to queue
if [ $nosubmit != 'on' ]; then
    ./case.submit
fi

