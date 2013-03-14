#!/bin/bash

# These variables are set by CMake
HOMME_DIR=@Homme_Build_DIR@
HOMME_TEST_RESULTS=@Homme_Results_DIR@

# The location of the Netcdf reference files (if they exist)
HOMME_NC_RESULTS_DIR=@NETCDF_RESULTS_DIR@

# The location of the tests directory
HOMME_TESTING_DIR=${HOMME_DIR}/tests
cd $HOMME_TESTING_DIR

# The "type" of submission (lsf, pbs, standard mpi etc.) for creating the executable scripts 
HOMME_Submission_Type=@Homme_Submission_Type@

# Whether to use cprnc to diff the Netcdf files
USE_CPRNC=@TEST_USING_CPRNC@

# The cprnc Netcdf comparison tool
CPRNC_BINARY=@CPRNC_BINARY@

# The cprnc Netcdf comparison tool
PYTHON_EXECUTABLE=@PYTHON_EXECUTABLE@

# The testing utilities
source ${HOMME_DIR}/tests/testing-utils.sh

# Get the argument
TEST_NAME=$1
echo "Test name = ${TEST_NAME}"

echo "Diffing the stdout of the run:"
diffStdout

echo "############################################################################"
echo "  The diff of the stdout has passed"
echo "############################################################################"

if [ "${USE_CPRNC}" == ON -o "${USE_CPRNC}" == TRUE ] ; then
  echo "Diffing the Netcdf output files"
  diffCprnc
  echo "############################################################################"
  echo "  The diff using CPRNC has passed"
  echo "############################################################################"

else
  echo "Not diffing the Netcdf output"
  echo "############################################################################"
fi

exit 0
