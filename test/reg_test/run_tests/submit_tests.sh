#!/bin/bash

# These variables are set by CMake
HOMME_DIR=@Homme_Build_DIR@
HOMME_TEST_RESULTS=@Homme_Results_DIR@

HOMME_NCRESULTS_DIR=@Homme_NCResults_DIR@

# The location of the tests directory
HOMME_TESTING_DIR=${HOMME_DIR}/tests
cd $HOMME_TESTING_DIR

# The "type" of submission (lsf, pbs, standard mpi etc.) for creating the executable scripts 
HOMME_Submission_Type=@Homme_Submission_Type@

# The cprnc Netcdf comparison tool
CPRNC_BINARY=@CPRNC_BINARY@

if [ "$1" == all ] ; then
  SUBMIT_ALL_AT_ONCE=true
else
  SUBMIT_ALL_AT_ONCE=false
fi

if [ "${SUBMIT_ALL_AT_ONCE}" == true ] ; then

  # The lists of tests to run
  source ${HOMME_DIR}/tests/submission-list.sh

else

  num_submissions=1
  subFile1=$1
  testName=$2

  # To Do: make sure the above file exists

fi

if [ "$HOMME_Submission_Type" = lsf -o "$HOMME_Submission_Type" = pbs ]; then
  HOMME_QUEUE_JOBS=true
else
  HOMME_QUEUE_JOBS=false
fi

# The testing utilities
source ${HOMME_DIR}/tests/testing-utils.sh

if [ "$HOMME_QUEUE_JOBS" = true ]; then
  # Submit the tests to the queue
  submitTestsToQueue

  # Print a summary of the submissions
  printSubmissionSummary

  # Wait for the jobs to run through the queue
  queueWait

  # Diff the output files with those saved in the repo
else
  runTestsStd
fi

# parse the stdout to grab only the relevant info
parseStdout

if [ "${SUBMIT_ALL_AT_ONCE}" == true ] ; then

  # Do nothing for now
  # Pass
  :

else

  ${HOMME_TESTING_DIR}/diff_output.sh ${testName}

fi

exit 0
