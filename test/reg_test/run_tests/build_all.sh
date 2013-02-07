#!/bin/bash

# Read the tests names from a the test file
export TESTS
source ./test-list.in

#Grap the command line options
commandArgs="$@"

# common-setup.sh contains routines to build the tests
source ./common-setup.sh

# Set some variables
initBuild $commandArgs

# Here we loop through each tests using functions to peform each step
for thisTest in "${TESTS[@]}"
do

  resetVariables
  
  # File containing test specific data
  source $thisTest

  # Reconfigure/Rebuild or just build depending the previous state
  testConfigAndBuild

  # Copy test specific files to the build directory
  setupFileStructure

  # set up some MPI information
  setupMPI

  # on yellowstone (for now) create the lsf submission scripts
  if [[ `uname -n` =~ yslogin. ]]; then
    yellowstoneSetupLSF
  fi

done


