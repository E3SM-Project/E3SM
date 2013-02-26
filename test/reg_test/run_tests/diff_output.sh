#!/bin/bash

# These variables are set by CMake
HOMME_DIR=@Homme_Build_DIR@
HOMME_TEST_RESULTS=@Homme_Results_DIR@

HOMME_NC_RESULTS_DIR=/glade/scratch/jamroz/homme-results/yellowstone.intel

# The location of the tests directory
HOMME_TESTING_DIR=${HOMME_DIR}/tests

# The "type" of submission (lsf, pbs, standard mpi etc.) for creating the executable scripts 
HOMME_Submission_Type=@Homme_Submission_Type@

# The cprnc Netcdf comparison tool
CPRNC_BINARY=@CPRNC_BINARY@

# The lists of tests to run
source ${HOMME_DIR}/tests/submission-list.sh

# The testing utilities
source ${HOMME_DIR}/tests/testing-utils.sh

# Get the argument
TEST_NAME=$1
echo "TEST_NAME=${TEST_NAME}"

diffCprnc
