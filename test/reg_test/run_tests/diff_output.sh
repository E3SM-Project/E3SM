#!/bin/bash

# Get the build dir from cmake
HOMME_DIR=@Homme_Build_DIR@

# Source the following file to get the rest of the cmake variables
source ${HOMME_DIR}/tests/cmake_variables.sh

# Enter the testing dir
cd $HOMME_TESTING_DIR

# The "type" of submission (lsf, pbs, standard mpi etc.) for creating the executable scripts 
HOMME_Submission_Type=@HOMME_SUBMISSION_TYPE@

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

if [ "${USE_CPRNC}" == ON -o "${USE_CPRNC}" == TRUE ] ; then
  echo "Diffing the Netcdf output files"
  diffCprncOutput
  echo "############################################################################"
  echo "  The diff using CPRNC has passed"
  echo "############################################################################"

else
  echo "############################################################################"
  echo "Not diffing the Netcdf output"
  echo "############################################################################"
  echo "Diffing the stdout of the run:"
  echo "############################################################################"
  diffStdout
  echo "############################################################################"
  echo "  The diff of the stdout has passed"
  echo "############################################################################"



fi

exit 0
