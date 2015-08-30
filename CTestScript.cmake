#==============================================================================
#
#  This is the CTest script for PIO builds and submission to the CTest
#  Dashboard site: my.cdash.org.
#
#  Example originally stolen from:
#    http://www.vtk.org/Wiki/CTest:Using_CTEST_and_CDASH_without_CMAKE
#==============================================================================

#-------------------------------------------
#-- Get the common header information
#-------------------------------------------
set (CTEST_EXTRA_SCRIPT_PATH "${CTEST_SCRIPT_DIRECTORY}/ctest")
list (APPEND CMAKE_MODULE_PATH "${CTEST_EXTRA_SCRIPT_PATH}")
include (CTestScript-Header)

#-----------------------------------------------------------  
#-- Get build-specific information
#-----------------------------------------------------------  

## -- SRC Dir (where this script exists)
set (CTEST_SOURCE_DIRECTORY   "${CTEST_SCRIPT_DIRECTORY}")

## -- Empty the binary directory
ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})

## -- Add the CTest script directory to the module path
set (CTEST_ENVIRONMENT_SCRIPT "CTestEnvironment-${HOSTNAME_ID}")
set (CTEST_RUNCTEST_SCRIPT "${CTEST_EXTRA_SCRIPT_PATH}/runctest-${HOSTNAME_ID}.sh")

message ("CTEST_EXTRA_SCRIPT_PATH = ${CTEST_EXTRA_SCRIPT_PATH}")
message ("CTEST_SCRIPT_ARG = ${CTEST_SCRIPT_ARG}")

# -----------------------------------------------------------  
# -- Run CTest
# -----------------------------------------------------------  

## -- Start
message (" -- Start dashboard - ${CTEST_BUILD_NAME} --")
ctest_start("${CTEST_SCRIPT_ARG}")

## -- Update
message (" -- Update source - ${CTEST_BUILD_NAME} --")
set (CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")
ctest_update ()

## -- Configure 
message (" -- Configure build - ${CTEST_BUILD_NAME} --")
include (${CTEST_ENVIRONMENT_SCRIPT})
set (CTEST_CONFIGURE_COMMAND "${CMAKE_COMMAND} ${CTEST_CONFIGURE_OPTIONS} ${CTEST_SOURCE_DIRECTORY}")
ctest_configure ()

## -- BUILD
message (" -- Build - ${CTEST_BUILD_NAME} --")
set (CTEST_BUILD_COMMAND "${MAKE} tests")
ctest_build ()

## -- TEST
message (" -- Test - ${CTEST_BUILD_NAME} --")
execute_process (COMMAND "${CTEST_RUNCTEST_SCRIPT} ${CTEST_EXTRA_SCRIPT_PATH} ${CTEST_SCRIPT_ARG}"
                 WORKING_DIRECTORY ${CTEST_BINARY_DIRECTORY})

## -- SUBMIT
message (" -- Submit to dashboard - ${CTEST_BUILD_NAME} --")
ctest_submit ()

message (" -- Finished - ${CTEST_BUILD_NAME} --")
