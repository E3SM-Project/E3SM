#==============================================================================
#
#  This is the CTest script for generating test results for submission to the 
#  CTest Dashboard site: my.cdash.org.
#
#  Example originally stolen from:
#    http://www.vtk.org/Wiki/CTest:Using_CTEST_and_CDASH_without_CMAKE
#==============================================================================

#-------------------------------------------
#-- Get the common header information
#-------------------------------------------

list (APPEND CMAKE_MODULE_PATH ${CTEST_SCRIPT_DIRECTORY})
include (CTestScript-Header)

#-----------------------------------------------------------  
#-- Get build-specific information
#-----------------------------------------------------------  

## -- SRC Dir (where this script exists)
set (CTEST_SOURCE_DIRECTORY   "${CTEST_SCRIPT_DIRECTORY}/..")

# -----------------------------------------------------------  
# -- Run CTest- TESTING ONLY (Appended to existing TAG)
# -----------------------------------------------------------  

## -- Start
ctest_start("${CTEST_SCRIPT_ARG}" APPEND)

## -- TEST
ctest_test()
