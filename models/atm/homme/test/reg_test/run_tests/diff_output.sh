#!/bin/bash

# Get the build dir from cmake
HOMME_DIR=@Homme_Build_DIR@

# Source the following file to get the rest of the cmake variables
source ${HOMME_DIR}/tests/cmake_variables.sh

# Enter the testing dir
cd $HOMME_TESTING_DIR

# The testing utilities
source ${HOMME_DIR}/tests/testing-utils.sh

# Get the argument
TEST_NAME=$1
echo "Test name = ${TEST_NAME}"

echo "Diffing the Netcdf output files"
diffCprncOutput
echo "############################################################################"
echo "  The diff using CPRNC has passed"
echo "############################################################################"

exit 0
