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

# The lists of tests to run
source ${HOMME_DIR}/tests/submission-list.sh

# The testing utilities
source ${HOMME_DIR}/tests/testing-utils.sh

if [ "$HOMME_Submission_Type" = lsf ]; then
  # Submit the tests to the queue
  submitTestsToLSF

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

