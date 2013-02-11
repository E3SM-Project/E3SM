#!/bin/bash

# These variables are set by CMake
HOMME_DIR=@Homme_Build_DIR@
HOMME_TEST_RESULTS=@Homme_Results_DIR@

# The location of the tests directory
HOMME_TESTING_DIR=${HOMME_DIR}/tests

# The lists of tests to run
source lsf-list.sh

# The testing utilities
source testing-utils.sh

# Submit the tests to the queue
submitTestsToLSF

# Print a summary of the submissions
printSubmissionSummary

# Wait for the jobs to run through the queue
queueWait

# Diff the output files with those saved in the repo
diffStdOut


