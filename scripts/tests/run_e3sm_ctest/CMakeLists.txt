#
# Mini CMake project to test the run_e3sm script AND
# post the values on ACME-Climate CDash site.
#
# From any directory, do:
#   cmake <path_to_this_dir>
#   ctest -D Experimental
#
# Note: running ctest without -D will run test but not post to CDash site
#
# Questions to: Deakin (test script), Salinger (cmake)
#
cmake_minimum_required(VERSION 2.8)
set(BUILDNAME "run_e3sm_script_test" CACHE STRING "*Build Name* string on CDash")

# This script requires no language support, syntax has changed in CMake
if(${CMAKE_VERSION} VERSION_GREATER 3.3.0)
    set (LANG "LANGUAGES")
else()
    set (LANG "")
endif()

project(runE3smScriptTest ${LANG} NONE)

# Turn on testing and reporting to CDash -- see CTestConfig.cmake
ENABLE_TESTING()
INCLUDE(CTest)

# Get absolute path to dir where run_e3sm resides
set(RUN_E3SM_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../../.. CACHE PATH "Path to run_e3sm script")
# Create test script using configure_file to inject RUN_E3SM_DIR path
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/runE3smTest_template.sh
               ${CMAKE_CURRENT_BINARY_DIR}/runE3smTest.sh @ONLY)

add_test (run_e3sm_script_test runE3smTest.sh)
