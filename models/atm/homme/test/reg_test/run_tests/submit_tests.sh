#!/bin/bash

# Get the build dir from cmake
HOMME_DIR=@Homme_Build_DIR@

# Source the following file to get the rest of the cmake variables
source ${HOMME_DIR}/tests/cmake_variables.sh

# The testing utilities
source ${HOMME_DIR}/tests/testing-utils.sh

cd $HOMME_TESTING_DIR

# Determine if we are creating the baseline
if [ "$1" == baseline ] ; then
  CREATE_BASELINE=true
else
  CREATE_BASELINE=false
fi

# Determine if we are submitting all of the jobs at once
if [ "$1" == all -o ${CREATE_BASELINE} == true ] ; then
  SUBMIT_ALL_AT_ONCE=true
else
  SUBMIT_ALL_AT_ONCE=false
fi

# Either read the list of tests or set the arguments
if [ "${SUBMIT_ALL_AT_ONCE}" == true ] ; then

  # If we are submitting all at once read the prepared list of tests
  if [ "${CREATE_BASELINE}" == true ] ; then
    source ${HOMME_DEFAULT_BASELINE_DIR}/submission-list.sh
  else
    source ${HOMME_DIR}/tests/submission-list.sh
  fi

else

  # Only running one test
  num_submissions=1
  subFile1=$1
  TEST_NAME=$2

  # To Do: make sure the above file exists

fi

if [ "${HOMME_QUEUING}" = TRUE -o "${HOMME_QUEUING}" = ON ]; then
  # Submit the tests to the queue
  submitTestsToQueue

  # Print a summary of the submissions
  printSubmissionSummary

  # Wait for the jobs to run through the queue
  queueWait
   
else
  runTestsStd
fi

if [ "${SUBMIT_ALL_AT_ONCE}" == true ] ; then

  # Do nothing for now

  # If baseline then move the netcdf output files to the baseline dir
  if [ ${CREATE_BASELINE} == true ] ; then
    echo "Creating baseline..."
    #moveBaseline
  fi

else

  ${HOMME_TESTING_DIR}/diff_output.sh ${TEST_NAME}

fi

