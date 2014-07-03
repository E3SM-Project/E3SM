#!/bin/bash 

#==============================================================================
# $Id$
# $URL$
#
# Sets up and submits a 12-month  run; then calls create_clone 100 times, sets
# these up and submits them, for a total of 101 ensemble members
#
#==============================================================================

gen_random_numbers ()
{
  # For an even distribution, we don't want to allow numbers between max_rand
  # and 32767 (largest $RANDOM can be), otherwise distribution will be slightly
  # skewed towards smaller numbers: those smaller than 43 for first number,
  # smaller than 67 for second number, and smaller than 97 for final number.
  # Instead, we pick between 0 and floor(32768/N)*N-1. max_rand contains
  # floor(32768/N)*N-1 for N = 101, 100, and 99.
  max_rand=( 32723 32699 32669 )
    
  for i in `seq 0 2`
  do
    # Want a random number between 0 and N-1, inclusive
    N=$(( 101 - i ))
    
    # Pick a random number
    tmp_rand=$((RANDOM))
    # if number is too big (see comment relating to max_rand), pick again
    while [ $tmp_rand -gt ${max_rand[$i]} ]; do
      tmp_rand=$((RANDOM))
    done
    
    # Store random number mod N (=> between 0 and N-1)
    rand_ints+=( $(( tmp_rand % N )) )
  done
    
  # We want 3 numbers in {0,100}, but no duplicates. So we pick one number in
  # {0,100}, one in {0,99}, and one in {0,98}. If the second number is smaller
  # than the first we keep it as is, otherwise we increment it by 1 (so instead
  # of picking R2 in {0,99}, we effectively are picking it in the set
  # {0,R1-1} U {R1+1,100}). Similarly, we want to pick R3 in the set
  # {0,R1-1} U (R1+1,R2-1} U {R2+1,100} (assuming R1 < R2, which is reasonable
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

get_pertlim ()
{
  # Return a single pertlim value from an integer in the range of 0 - 100
  #  CASE   PERTLIM      CASE   PERTLIM
  #  ----   -------      ----   -------
  #   000   0.0
  #   001   1.0d-14       051  -1.0d-14
  #   002   1.1d-14       052  -1.1d-14
  #   003   1.2d-14       053  -1.2d-14
  #   004   1.3d-14       054  -1.3d-14
  #   005   1.4d-14       055  -1.4d-14
  #   006   1.5d-14       056  -1.5d-14
  #   007   1.6d-14       057  -1.6d-14
  #   008   1.7d-14       058  -1.7d-14
  #   009   1.8d-14       059  -1.8d-14
  #   010   1.9d-14       060  -1.9d-14
  #   011   2.0d-14       061  -2.0d-14
  #   012   2.1d-14       062  -2.1d-14
  #   013   2.2d-14       063  -2.2d-14
  #   014   2.3d-14       064  -2.3d-14
  #   015   2.4d-14       065  -2.4d-14
  #   016   2.5d-14       066  -2.5d-14
  #   017   2.6d-14       067  -2.6d-14
  #   018   2.7d-14       068  -2.7d-14
  #   019   2.8d-14       069  -2.8d-14
  #   020   2.9d-14       070  -2.9d-14
  #   021   3.0d-14       071  -3.0d-14
  #   022   3.1d-14       072  -3.1d-14
  #   023   3.2d-14       073  -3.2d-14
  #   024   3.3d-14       074  -3.3d-14
  #   025   3.4d-14       075  -3.4d-14
  #   026   3.5d-14       076  -3.5d-14
  #   027   3.6d-14       077  -3.6d-14
  #   028   3.7d-14       078  -3.7d-14
  #   029   3.8d-14       079  -3.8d-14
  #   030   3.9d-14       080  -3.9d-14
  #   031   4.0d-14       081  -4.0d-14
  #   032   4.1d-14       082  -4.1d-14
  #   033   4.2d-14       083  -4.2d-14
  #   034   4.3d-14       084  -4.3d-14
  #   035   4.4d-14       085  -4.4d-14
  #   036   4.5d-14       086  -4.5d-14
  #   037   4.6d-14       087  -4.6d-14
  #   038   4.7d-14       088  -4.7d-14
  #   039   4.8d-14       089  -4.8d-14
  #   040   4.9d-14       090  -4.9d-14
  #   041   5.0d-14       091  -5.0d-14
  #   042   5.1d-14       092  -5.1d-14
  #   043   5.2d-14       093  -5.2d-14
  #   044   5.3d-14       094  -5.3d-14
  #   045   5.4d-14       095  -5.4d-14
  #   046   5.5d-14       096  -5.5d-14
  #   047   5.6d-14       097  -5.6d-14
  #   048   5.7d-14       098  -5.7d-14
  #   049   5.8d-14       099  -5.8d-14
  #   050   5.9d-14       100  -5.9d-14

  i=$1
  if (( $i == 0 )); then 
    ptlim=0
  elif (( $i < 51 )); then
    let j=i+9
    ippt=`/usr/bin/printf "%2.2d" $j`
    ptlim="0.${ippt}d-13"
  else
    let j=i-41
    ippt=`/usr/bin/printf "%2.2d" $j`
    ptlim="-0.${ippt}d-13"
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
    firstpertlim=$( get_pertlim ${rand_ints[0]} )
    . $SingleRun "$@" -p $firstpertlim
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
      PERTLIM=$(get_pertlim ${rand_ints[$i]}  )
    else 
      PERTLIM=$(get_pertlim $i ) 
    fi

    CASE1_NAME=$CASE_PFX.$iens
    CASE1=$CASE_ROOT/$CASE1_NAME

    # Create clone
    cd $SCRIPTS_ROOT
    ./create_clone -case $CASE1 -clone $CASE # Copy $CASE to $CASE1

    # Get value for EXEROOT from $CASE
    # Note return string is "EXEROOT = $EXEROOT"
    if [ $test_suite = 'TRUE' ]; then
      cd $SCRIPTS_ROOT/$CASE
    else
      cd $CASE
    fi
    EXE=`./xmlquery EXEROOT -valonly`
    EXEROOT=${EXE#"EXEROOT = "}

    # Edit env_build in cloned case
    if [ $test_suite = 'TRUE' ]; then
      cd $SCRIPTS_ROOT/$CASE1
    else
      cd $CASE1
    fi
    ./xmlchange -file env_build.xml -id EXEROOT -val $EXEROOT
    ./xmlchange -file env_build.xml -id BUILD_COMPLETE -val TRUE
    ./cesm_setup

    # For validations, subsequent cloned cases will have the pertlim from the
    # parent case. We neet to remove the original case's pertlim before we set
    # the new value
    if [ $runtype = 'validation' ]; then
      sed  -i '/pertlim/d' user_nl_cam
    fi

    # Change pertlim in clone
    echo "pertlim = $PERTLIM" >> user_nl_cam
    ./preview_namelists

    # Adjust walltime, account number and ptile in clone
    fix_run_script $CASE1_NAME

    # Only submit the cloned case if --nosubmit is off
    if [ $nosubmit != 'on' ]; then
      ./$CASE1_NAME.submit
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

# Default assumes running from command line, if the -test_suite option is set
# then a couple of directories change
test_suite="FALSE"

# globals to disable building and submitting, this will be needed if and when 
# this script gets used in the test suite. 
nobuild="off"
nosubmit="off"

# Default ensemble count is 101, 0-100. This is changed to 2 when runtype is
# 'validation'

# Process input arguments
Args=("$@")
i=0
while [ $i -le ${#Args[@]} ]; do
  case ${Args[$i]} in
    -case )
      i=$((i+1))
      CASENAME=$(basename ${Args[$i]})
      if [ ! ${CASENAME##*.} == "000" ]; then
        echo "ERROR: When building an ensemble, \$CASE must end in 000!"
        exit 2
      fi
    ;;
    -test_suite )
      # Set test_suite to TRUE
      test_suite="TRUE"
    ;;
    -ensemble )
      # Set CLONECOUNT and runtype
      CLONECOUNT=100
      runtype="ensemble"
    ;;
    -nb|-nobuild )
      nobuild="on"         
      nosubmit="on"
    ;;
    -ns|-nosubmit )
      nosubmit="on"
    ;;
  esac
  i=$((i+1))
done

# If runtype is 'validation', print the three choices of pertlim to screen
if [ $runtype = 'validation' ]; then
    # Create empty array; this will be populated with 3 random integers
    # between 0 and 100, inclusive
    rand_ints=()

    gen_random_numbers
fi

# Create all the cases
create_cases "$@"

if [ $runtype = 'validation' ]; then
  echo "---"
  echo "Set up three cases using the following pertlim values:"            \
       "$( get_pertlim ${rand_ints[0]} ) $( get_pertlim ${rand_ints[1]} )"  \
       "$( get_pertlim ${rand_ints[2]} )"
fi

