#!/bin/bash 

#==============================================================================
# $Id$
# $URL$
#
# Sets up and submits a 12-month  run; then calls create_clone 150 times, sets
# these up and submits them, for a total of 151 ensemble members
#
#==============================================================================

gen_random_numbers ()
{
  # For an even distribution, we don't want to allow numbers between max_rand
  # and 32767 (largest $RANDOM can be), otherwise distribution will be slightly
  # skewed towards smaller numbers: those smaller than 43 for first number,
  # smaller than 67 for second number, and smaller than 97 for final number.
  # Instead, we pick between 0 and floor(32768/N)*N-1. max_rand contains
  # floor(32768/N)*N-1 for N = 151, 150, and 149.
  max_rand=( 32723 32699 32669 )

  if [[ $1 -eq 1 ]]; then
     COUNT=349
  else
     COUNT=151
  fi
  for i in `seq 0 2`
  do
    # Want a random number between 0 and N-1, inclusive
    N=$(( $COUNT - i ))

    # Pick a random number
    tmp_rand=$((RANDOM))
    # if number is too big (see comment relating to max_rand), pick again
    while [ $tmp_rand -gt ${max_rand[$i]} ]; do
      tmp_rand=$((RANDOM))
    done

    # Store random number mod N (=> between 0 and N-1)
    rand_ints+=( $(( tmp_rand % N )) )
  done

  # We want 3 numbers in {0,150}, but no duplicates. So we pick one number in
  # {0,150}, one in {0,149}, and one in {0,148}. If the second number is smaller
  # than the first we keep it as is, otherwise we increment it by 1 (so instead
  # of picking R2 in {0,149}, we effectively are picking it in the set
  # {0,R1-1} U {R1+1,105}). Similarly, we want to pick R3 in the set
  # {0,R1-1} U (R1+1,R2-1} U {R2+1,150} (assuming R1 < R2, which is reasonable
  # because we know R1 != R2).

  # If first and second number are same, increment second number by 1
  if [ ${rand_ints[1]} -ge ${rand_ints[0]} ]; then
    rand_ints[1]=$((rand_ints[1]+1))
  fi

  # If the third number is larger than both the first and second, increment the
  # third number by 2
  # Otherwise, if the third number is only larger than one of the first two,
  # increment it by 1
  if [ ${rand_ints[2]} -ge ${rand_ints[0]} ] &&                               \
     [ ${rand_ints[2]} -ge ${rand_ints[1]} ]; then
    rand_ints[2]=$((rand_ints[2]+2))
  elif [ ${rand_ints[2]} -ge ${rand_ints[0]} ] ||                             \
       [ ${rand_ints[2]} -ge ${rand_ints[1]} ]; then
    rand_ints[2]=$((rand_ints[2]+1))
    # Check for cases where incrementing the third number by one makes it equal
    # to the number it was previously smaller than. In these cases, increment
    # by 1 again.
    if [ ${rand_ints[2]} -eq ${rand_ints[0]} ] ||                             \
       [ ${rand_ints[2]} -eq ${rand_ints[1]} ]; then
      rand_ints[2]=$((rand_ints[2]+1))
    fi
  fi
}


get_pertlim_uf ()
{
  i=$1
  if (( $i == 0 )); then
    ptlim=0
  else
    j=$(( 2*((i-1)/100) + 101 ))
    k=$(( (i-1) % 100 ))
    if [ $(( i % 2)) -ne 0 ]; then
      l=$(( j+(k/2)*18 ))
      ippt=$( printf "%0*d" 3 $l )
      ptlim="0.${ippt}d-13"
    else
      l=$(( j+((k-1)/2)*18 ))
      ippt=$( printf "%0*d" 3 $l )
      ptlim="-0.${ippt}d-13"
    fi
  fi
  echo $ptlim
}

create_cases ()
{
  # single_run.sh needs to be in same directory as ensemble.sh
  ThisDir=$(cd `dirname $0`; pwd -P )
  SingleRun="$ThisDir/single_run.sh"
  if [ ! -e $SingleRun ]; then
    echo "ERROR: can not find script to produce single run in $ThisDir"
    exit 1
  fi

  # If we're doing validation, we have to give the first run a random pertlim
  # and set the arguments appropriately, otherwise, no pertlim value needed for
  # the first case
  if [ $runtype = 'validation' ]; then
    firstpertlim=$( get_pertlim_uf ${rand_ints[0]} )
    . $SingleRun "$@" --pertlim $firstpertlim
    Status=$?
  else
    . $SingleRun "$@"
    Status=$?
  fi

  if [ "$Status" != "0" ]; then
    echo "Exit: $Status"
    exit $Status
  fi

  ######### END OF BUILDING ROOT CASE, NOW CLONING

  CASE_ROOT=$(dirname $CASE)
  CASE_PFX=$(basename $CASE .000) # It was checked above that the CASENAME suffix is "000"
  for i in `seq 1 $CLONECOUNT`; do
    iens=`/usr/bin/printf "%3.3d" $i`
    # If we're doing validations, set the pertlim based on
    # the random numbers previously generated, otherwise
    # just get a pertlim value in sequence.
    if [ $runtype = 'validation' ]; then
      PERTLIM=$(get_pertlim_uf ${rand_ints[$i]}  )
    else
      PERTLIM=$(get_pertlim_uf $i )
    fi

    CASE1_NAME=$CASE_PFX.$iens
    CASE1=$CASE_ROOT/$CASE1_NAME

    # Create clone
    cd $SCRIPTS_ROOT
    echo "=== SCRIPTS_ROOT ==="
    echo $SCRIPTS_ROOT
    #./create_clone --keepexe --case $CASE1 --clone $CASE # Copy $CASE to $CASE1
    $SCRIPTS_ROOT/create_clone --keepexe --case $CASE1 --clone $CASE # Copy $CASE to $CASE1

    # Get value for EXEROOT from $CASE
    # Note return string is "EXEROOT = $EXEROOT"
    if [ $test_suite = 'TRUE' ]; then
      cd $SCRIPTS_ROOT/$CASE
    else
      cd $CASE
    fi
    EXE=`./xmlquery EXEROOT -value`
    EXEROOT="$(echo -e "$EXE" | sed -e 's/[[:space:]]*$//')"

    # Edit env_build in cloned case
    if [ $test_suite = 'TRUE' ]; then
      cd $SCRIPTS_ROOT/$CASE1
    else
      cd $CASE1
    fi
	#echo "running ./xmlchange EXEROOT=\"$EXEROOT\""
    #./xmlchange EXEROOT="$EXEROOT"
    #echo "running ./xmlchange BUILD_COMPLETE=\"TRUE\""
    #./xmlchange BUILD_COMPLETE="TRUE"
    echo "running case.setup"
    ./case.setup

    # For validations, subsequent cloned cases will have the pertlim from the
    # parent case. We neet to remove the original case's pertlim before we set
    # the new value
    if [ $runtype = 'validation' ]; then
      sed  -i '/pertlim/d' user_nl_cam
    fi

    echo "mfilt = 2">>user_nl_cam
    # Change pertlim in clone
    echo "pertlim = $PERTLIM" >> user_nl_cam
    ./preview_namelists

    # Adjust walltime, account number and ptile in clone
    #fix_run_script $CASE1_NAME
    fix_run_script case

    # Only submit the cloned case if --nosubmit is off
    if [ $nosubmit != 'on' ]; then
      ./case.submit
    fi

  done
}

#==============================================================================
# Main Program
#==============================================================================

# Make sure SCRIPTS_ROOT is set
if [[ -z $SCRIPTS_ROOT ]]; then
  SCRIPTS_ROOT=`pwd`
fi

# The default runtype is 'validation', if the -ensemble option is set,
# we change this to 'ensemble'. We make two clones (to have three runs).
runtype="validation"
CLONECOUNT=2
UF=0

# Default assumes running from command line, if the -test_suite option is set
# then a couple of directories change
test_suite="FALSE"

# globals to disable building and submitting, this will be needed if and when
# this script gets used in the test suite.
nobuild="off"
nosubmit="off"

# Default ensemble count is 151, 0-150. This is changed to 2 when runtype is
# 'validation'

# Process input arguments
Args=("$@")
i=0
while [ $i -le ${#Args[@]} ]; do
  case ${Args[$i]} in
    --case )
      i=$((i+1))
      CASENAME=$(basename ${Args[$i]})
      if [ ! ${CASENAME##*.} == "000" ]; then
        echo "ERROR: When building an ensemble, \$CASE must end in 000!"
        exit 2
      fi
    ;;
    --test_suite )
      # Set test_suite to TRUE
      test_suite="TRUE"
    ;;
    --ensemble )
      i=$((i+1))
      # Set CLONECOUNT and runtype
      #CLONECOUNT=${Args[$i]}-1
      CLONECOUNT=$((${Args[$i]}-1))
      if [ $CLONECOUNT -gt 999 ]; then
        echo "ERROR: the number of ensemble member cannot be set to more than 999!"
        exit 2
      fi
      runtype="ensemble"
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
  esac
  i=$((i+1))
done

# If runtype is 'validation', print the three choices of pertlim to screen
if [ $runtype = 'validation' ]; then
    # Create empty array; this will be populated with 3 random integers
    # between 0 and 150, inclusive
    rand_ints=()

    gen_random_numbers $UF
fi

# Create all the cases
create_cases "$@"

if [ $runtype = 'validation' ]; then
  echo "---"
  echo "Set up three cases using the following pertlim values:"            \
     "$( get_pertlim_uf ${rand_ints[0]} ) $( get_pertlim_uf ${rand_ints[1]} )"  \
     "$( get_pertlim_uf ${rand_ints[2]} )"
fi

