#!/bin/bash

#==============================================================================
# $Id$
# $URL$
#
# Sets up a 12-month ne30_ne30 run
#
#==============================================================================

#==============================================================================
# Usage subroutine
#==============================================================================

usage() {
  echo "USAGE: $ThisScript -case CASE -mach MACH [-p PERTLIM] [-mach_pes ENV_MACH_PES_FILE] [-np NP] [-npocn NPOCN] [-npice NPICE] [-nt NTHRDS] [-w WALLTIME] [-account ACCOUNT_NUM] [-compiler COMPILER] [-compset COMPSET -r RES] [-nb] [-ns] [-usr_nl_cam \"VARIABLE1=VALUE1 VARIABLE2=VALUE2\"] [-cam_src_mod \"FILE1 FILE2 ... FILEN\"]"
  echo ''
  echo 'Sets up a 1-year run to compare to an ensemble generated on a trusted'\
       'machine. Need to automate generation of NCL script to run comparison'
  echo ''
  echo 'Required flags:'
  echo '  -case           Case name passed on to create_newcase'
  echo '  -mach           Machine name passed on to create_newcase'
  echo ''
  echo 'Optional flags:'
  if [ $ThisScript = "single_run.sh" ]; then
    echo '  -p,        --pertlim     Run CAM with non-zero pertlim'
  fi
  echo '  -mach_pes  --mach_pes    Use specified env_mach_pes file'
  echo '  -np,       --numproc     Number of processors requested (default = Machines default)'
  echo '  -npocn,    --npocn       Number of processors requested for ocean (default = Machines default)'
  echo '  -npice,    --npice       Number of processors requested for ice (default = Machines default)'
  echo '                           Note: npice and npocn will = np if np is specified and other two are not'
  echo '  -nt,       --nthrds      Number of threads for each component (default = Machines default)'
  echo '  -w,        --walltime    Amount of walltime requested (default = 4:30)'
  echo '  -account,  --account     Account number to use in job scripts (default = machine default)'
  echo '  -compiler, --compiler    Compiler to use (passed on to create_newcase)'
  echo '  -compset,  --compset     Compset to use (default = BC5)'
  echo '  -res,      --res         Resolution to run (default = ne30_g16)'
  if [ $ThisScript = "ensemble.sh" ]; then
    echo '  -nb,       --nobuild     Disables building the root case of the ensemble.'
    echo '  -ns,       --nosubmit    Disables submitting any members of the ensemble.'
    echo '  -ensemble                Instead of building 3 cases with random pertlim values, build the 101 member ensemble'
    echo ' -test_suite               Flag to indicate this is run from the CESM test script (makes assumption about case directory)'
  else
    echo '  -nb,       --nobuild     Disables building the single case.'
    echo '  -ns,       --nosubmit    Disables submitting the single case.'
  fi
  echo '  -usr_nl_cam              Include specified variables setting in user_nl_cam'
  echo '  -cam_src_mod             Include specified files in CAM SourceMod dir'
  echo '  -h,        --help        Output this usage information'
}

#==============================================================================
# Subroutine to fix run script
#==============================================================================

fix_run_script() {

  if [ -z "$1" ]; then
    #FileName=$CASENAME
    FileName='case'
  else
    FileName=$1
  fi
  echo "Checking to see if anything needs to be changed in $FileName.run..."

  # Adjust walltime
  if [ ! -z "$WalltimePrefix" ]; then
    cp $FileName.run $FileName.tmp
    echo "Changing wall time in $FileName.run to $WallTime"
    sed s/"$WalltimePrefix".\*/"${WalltimePrefix}$WallTime"/ $FileName.tmp > $FileName.run
    rm -f $FileName.tmp
  else
    echo "WARNING: Do not know how to change requested walltime on $MACH,"\
         "Just using default."
  fi

  # Adjust account (if requested)
  if [ ! -z "$ACCOUNT" ]; then
    if [ -z "$AccountPrefix" ]; then
      echo "ERROR: can not set account on $MACH because \$AccountPrefix is not set."
      exit 1
    fi
    echo "Setting account in $FileName.run to $ACCOUNT"
    cp $FileName.run $FileName.tmp
    sed s/"$AccountPrefix".\*/"${AccountPrefix}$ACCOUNT"/ $FileName.tmp > $FileName.run
    rm -f $FileName.tmp
  fi

  # For Yellowstone, set ptile=16 instead of 32
#  if [ "$MACH" == "yellowstone" ]; then
#    cp $FileName.run $FileName.tmp
#    echo "Changing ptile in $FileName.run to 16"
#    sed s/"ptile".\*/"ptile=16]\""/ $FileName.tmp > $FileName.run
#    rm -f $FileName.tmp
#  fi

}

#==============================================================================
# Main Program
#==============================================================================

ThisScript=$(basename $0)
ThisDir=$(cd `dirname $0`; pwd -P )

# Default Values
RES="ne30_g16"
COMPSET="BC5"
CHANGE_NP=0
CHANGE_NPICE=0
CHANGE_NPOCN=0
# MNL note: at this time, no flag to change ROOTPE all to 0
#           (no reason to change from default layout!)
CHANGE_ROOTPE=0
CHANGE_NTHRDS=0
WallTime="4:30"
PERTLIM="0"
nobuild="off"
nosubmit="off"
usr_nl_cam_val=()
cam_src_files=("")
POP_OUT_ON=0
CICE_OUT_ON=0
RTM_OUT_ON=0
# Large ensemble?
LENS=0

# SCRIPTS_ROOT is only used for generating ensemble!
SCRIPTS_ROOT=$(cd `dirname $0`; cd ../../scripts; pwd )
echo "CESM command: create_newcase $NewCaseFlags"
echo "(Running from: $SCRIPTS_ROOT)"
# Process input arguments
while [ $# -gt 0 ]; do
  case $1 in
    -scripts_root )
      SCRIPTS_ROOT=$2
      shift
    ;;
    -LENS )
      LENS=1
      POP_OUT_ON=1
      CICE_OUT_ON=1
      RTM_OUT_ON=1
    ;;
    -case )
      NewCaseFlags="$NewCaseFlags $1 $2"
      CASE=$2
      CASENAME=$(basename $CASE)
      shift
    ;;
    -mach )
      NewCaseFlags="$NewCaseFlags $1 $2"
      MACH=$2
      shift
    ;;
    -p|--pertlim )
      PERTLIM="$2"
      shift
    ;;
    -np|--numproc )
      CHANGE_NP=1
      NP=$2
      shift
    ;;
    -npocn|--npocn )
      CHANGE_NPOCN=1
      NPOCN=$2
      shift
    ;;
    -npice|--npice )
      CHANGE_NPICE=1
      NPICE=$2
      shift
    ;;
    -nt|--nthrds )
      CHANGE_NTHRDS=1
      NTHRDS=$2
      shift
    ;;
    -w|--walltime )
      WallTime=$2
      shift
    ;;
    -account|--account )
      ACCOUNT="$2"
      shift
    ;;
    -compiler|--compiler )
      NewCaseFlags="$NewCaseFlags $1 $2"
      shift
    ;;
    -compset|--compset )
      COMPSET=$2
      shift
    ;;
    -res|--res )
      RES=$2
      shift
    ;;
    -nb|--nobuild )
       nobuild="on"
       nosubmit="on"
    ;;
    -ns|--nosubmit )
       nosubmit="on"
    ;;
    -mach_pes|--mach_pes)
      arg_in=$2
      env_mach_pes_file=$(cd `dirname $arg_in`; pwd -P )/`basename $arg_in`
      shift
    ;;
    -usr_nl_cam )
      usr_nl_cam_val=("$2")
      shift
    ;;
    -cam_src_mod|--cam_src_mod )
      cam_src_files=()
      for arg_in in $2; do
        cam_src_files+=($(cd `dirname $arg_in`; pwd -P )/`basename $arg_in`)
      done
      shift
    ;;
    -h )
      usage
      exit 0
    ;;
    -test_suite|-ensemble )
      # Ignore these flags, they will be passed from ensemble.sh
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

if [ -z "$CASE" ]; then
  echo "Must specify -case argument"
  echo "Invoke $ThisScript -h for usage"
  exit 1
fi

if [ -z "$MACH" ]; then
  echo "Must specify -mach argument"
  echo "Invoke $ThisScript -h for usage"
  exit 1
fi

# SET NPOCN and $NPICE if NP is changed
if [ $CHANGE_NP -eq 1 ]; then
  if [ $CHANGE_NPOCN -eq 0 ]; then
    CHANGE_NPOCN=1
    NPOCN=$NP
  fi
  if [ $CHANGE_NPICE -eq 0 ]; then
    CHANGE_NPICE=1
    NPICE=$NP
  fi
fi

# Need to set strings used by sed to update run script (machine-dependent)
case $MACH in
  "bluefire" )
    WalltimePrefix="-W "
    AccountPrefix="-P "
  ;;
  "yellowstone" )
    WalltimePrefix="-W "
    AccountPrefix="-P "
  ;;
  "janus" )
    WalltimePrefix="-l walltime="
    WallTime="$WallTime:00"
  ;;
  "frankfurt" )
    WalltimePrefix="-l walltime="
    WallTime="$WallTime:00"
  ;;
  "hopper" )
    WalltimePrefix="-l walltime="
    WallTime="$WallTime:00"
  ;;
  "mira" )
    WalltimePrefix="set wt = "
  ;;
  * )
    echo "WARNING: Nothing machine-specific has been configured for $MACH"
  ;;
esac

NewCaseFlags="$NewCaseFlags -res $RES -compset $COMPSET"
cd $SCRIPTS_ROOT
echo "Currently in $SCRIPTS_ROOT"
echo "Flags for create_newcase are $NewCaseFlags"
./create_newcase $NewCaseFlags

cd $CASE
chmod u+w *

if [ ! -e env_mach_pes.xml.orig ]; then
  cp env_mach_pes.xml env_mach_pes.xml.orig
fi

if [ -e $env_mach_pes_file ] && [ ! -z $env_mach_pes_file ]; then
  rm -f env_mach_pes.xml
  cp -v $env_mach_pes_file env_mach_pes.xml
elif [ ! -z $env_mach_pes_file ]; then
  echo "WARNING: can not find $env_mach_pes_file to replace env_mach_pes.xml"
else
  if [ $CHANGE_NP -eq 1 ]; then
    echo "Changing to run ATM, LND, ROF, CPL, WAV, and GLC on $NP tasks"
    # Change env_mach_pes to run on NP processes
    ./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val $NP
    ./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val $NP
    ./xmlchange -file env_mach_pes.xml -id NTASKS_CPL -val $NP
    ./xmlchange -file env_mach_pes.xml -id NTASKS_GLC -val $NP
    ./xmlchange -file env_mach_pes.xml -id NTASKS_ROF -val $NP
    ./xmlchange -file env_mach_pes.xml -id NTASKS_WAV -val $NP
  fi
  if [ $CHANGE_NPOCN -eq 1 ]; then
    echo "Changing to run OCN on $NPOCN tasks"
    ./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val $NPOCN
  fi
  if [ $CHANGE_NPICE -eq 1 ]; then
    echo "Changing to run ICE on $NPICE tasks"
    ./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val $NPICE
  fi

  if [ $CHANGE_ROOTPE -eq 1 ]; then
    # Set rootpe = 0 for all components
    ./xmlchange -file env_mach_pes.xml -id ROOTPE_ATM -val 0
    ./xmlchange -file env_mach_pes.xml -id ROOTPE_LND -val 0
    ./xmlchange -file env_mach_pes.xml -id ROOTPE_ROF -val 0
    ./xmlchange -file env_mach_pes.xml -id ROOTPE_ICE -val 0
    ./xmlchange -file env_mach_pes.xml -id ROOTPE_OCN -val 0
    ./xmlchange -file env_mach_pes.xml -id ROOTPE_CPL -val 0
  fi

  # Set nthrds = $NTHRDS for all components
  if [ $CHANGE_NTHRDS -eq 1 ]; then
    ./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM -val $NTHRDS
    ./xmlchange -file env_mach_pes.xml -id NTHRDS_LND -val $NTHRDS
    ./xmlchange -file env_mach_pes.xml -id NTHRDS_ROF -val $NTHRDS
    ./xmlchange -file env_mach_pes.xml -id NTHRDS_ICE -val $NTHRDS
    ./xmlchange -file env_mach_pes.xml -id NTHRDS_OCN -val $NTHRDS
    ./xmlchange -file env_mach_pes.xml -id NTHRDS_CPL -val $NTHRDS
    ./xmlchange -file env_mach_pes.xml -id NTHRDS_GLC -val $NTHRDS
    ./xmlchange -file env_mach_pes.xml -id NTHRDS_WAV -val $NTHRDS
  fi
fi

# Change env_build to ensure bit-for-bit, disable archiving, disable restart
# files, and only run for 1 month
if [ ! -e env_run.xml.orig ]; then
  cp env_run.xml env_run.xml.orig
fi

./xmlchange -file env_run.xml -id BFBFLAG -val TRUE
./xmlchange -file env_run.xml -id DOUT_S -val FALSE
./xmlchange -file env_run.xml -id STOP_OPTION -val nmonths
./xmlchange -file env_run.xml -id STOP_N -val 12
./xmlchange -file env_run.xml -id REST_OPTION -val never

# Unset LS_COLORS before configure runs
LS_COLORS=
./case.setup || { echo "Error running case_setup, probably a bad NP value"; exit 1; }
chmod u+w user_nl_*
if [ $LENS -eq 1 ]; then
  cp -f $ThisDir/user_nl_cam_LENS ./user_nl_cam
fi

# Only edit user_nl_cam if the file exists (otherwise not using CAM!)
if [ -e user_nl_cam ]; then
  # Have CAM output everything (single precision)
  # But not initial data...
  echo "avgflag_pertape = 'A'" >> user_nl_cam
  echo "nhtfrq  = -8760" >> user_nl_cam
  echo "inithist = 'NONE'" >> user_nl_cam
  #echo "ndens  = 1" >> user_nl_cam # Double precision
  if [ ! "$PERTLIM" == "0" ]; then
    echo "pertlim = $PERTLIM" >> user_nl_cam
  fi

  # Add user-defined variables to user_nl_cam
  # (If you did not use the -usr_nl_cam option then $usr_nl_cam_val is empty)
  for new_variable in $usr_nl_cam_val; do
    echo "Adding additional variable to user_nl_cam!"
    echo "$new_variable" >> user_nl_cam
  done
fi

# Only edit user_nl_cam if the file exists (otherwise not using CAM!)
if [ -e user_nl_clm ]; then
  # Have CLM output everything (single precision)
  echo "hist_avgflag_pertape = 'A'" >> user_nl_clm
  echo "hist_nhtfrq  = -8760" >> user_nl_clm

  # Disable output?
  #echo "hist_empty_htapes = .true." >> user_nl_clm
fi

# Disable output from all other components
if [ -e user_nl_cice ] && [ $CICE_OUT_ON -eq 0 ]; then
  echo "histfreq = 'x','x','x','x','x'" >> user_nl_cice
fi
if [ -e user_nl_pop2 ] && [ $POP_OUT_ON -eq 0 ]; then
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
#echo "do_rtm = .false." >> user_nl_rtm

./preview_namelists || { echo "Error running preview_namelists"; exit 1; }

fix_run_script

# Copy any files over to SourceMods
if [ -d SourceMods/src.cam ]; then
  # (If you did not use the -cam_src_mod option then $cam_src_files is empty)
  for file in ${cam_src_files[@]}; do
    if [ -e $file ]; then
      cp -fv $file SourceMods/src.cam/
    else
      echo "WARNING: $file not found!"
      echo "(Probably a relative path issue, currently in `pwd -P`)"
    fi
  done
fi

# IGNORE THIS FOR NOW!
if [ 0 == 1 ]; then
  if [ -z "$EnsembleAvgsFile" ]; then
    echo "WARNING: not specifying location of netCDF file w/ avgs and standard"\
         "deviations from ensemble on known good machine. Script will skip checks"
  else
    # Add NCL script to run (if ensemble values are stored somewhere)
    # Check to see if NCL exists
    NCLCHECK=`command -v ncl || echo NotFound`
    if [ "$NCLCHECK" == "NotFound" ]; then
      echo "Can not find ncl in PATH"
      echo "This tool uses ncl to compare output to an ensemble"
    else
      # Note: need NCL 6.1.0 or later!
      NCL_VER=`ncl -V`
      NCL_FULL=(`echo $NCL_VER | tr '.' ' '`)
      NCL_MAIN=${NCL_FULL[0]}
      NCL_SUB=${NCL_FULL[1]}
  #    echo "ncl check_variance.ncl | tee check_variance.output" >> $CASENAME.run
    fi
  #  cat > check_variance.ncl << EOF
  #; Update this for annual averages!
  #EOF
  fi
fi

echo "nobuild: $nobuild, nosubmit: $nosubmit"
# Build executable
if [ $nobuild != 'on' ]; then
    if [ "$MACH" == "janus" ]; then
      ssh janus-compile1 "cd $PWD; ./$CASENAME.build" || { echo "Error building!"; exit 1; }
    else
      ./case.build || { echo "Error building!"; exit 1; }
    fi
fi

# Submit job to queue
if [ $nosubmit != 'on' ]; then
    if [ "$MACH" == "mira" ]; then
      ./case.run
    else
      ./case.submit
    fi
fi

