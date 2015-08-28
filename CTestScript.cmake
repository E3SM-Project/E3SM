#==============================================================================
#
#  This is the CTest script for Nightly builds and submission to the CTest
#  Dashboard site: my.cdash.org.
#
#  Example originally stolen from:
#    http://www.vtk.org/Wiki/CTest:Using_CTEST_and_CDASH_without_CMAKE
#==============================================================================

#---------------------------------------
#-- Get the machine environment
#---------------------------------------

## -- Set hostname

find_program (HOSTNAME_CMD NAMES hostname)
execute_process (COMMAND ${HOSTNAME_CMD}
                 OUTPUT_VARIABLE HOSTNAME)
                 
set (CTEST_SITE "${HOSTNAME}")

## -- Set hostname ID (e.g., alcf, nwsc, nersc, ...)

# UCAR/NWSC Machines
if (CTEST_SITE MATCHES "^yslogin" OR
    CTEST_SITE MATCHES "^geyser" OR
    CTEST_SITE MATCHES "^caldera")
    set (CTEST_SITE_ID "nwsc")
# ALCF/Argonne Machines
elseif (CTEST_SITE MATCHES "^mira" OR
        CTEST_SITE MATCHES "^cetus" OR
        CTEST_SITE MATCHES "^vesta" OR
        CTEST_SITE MATCHES "^cooley")
    set (CTEST_SITE_ID "alcf")
# ALCF/Argonne Machines
elseif (CTEST_SITE MATCHES "^edison" OR
        CTEST_SITE MATCHES "^carver" OR
        CTEST_SITE MATCHES "^hopper")
    set (CTEST_SITE_ID "nersc")
else ()
    set (CTEST_SITE_ID "unknown")
endif ()

## -- Set site / build name

find_program (UNAME NAMES uname)
function (getuname name flag)
  execute_process (COMMAND ${UNAME} ${flag}
                   OUTPUT_VARIABLE res)
  string (STRIP "${res}" res)
  set (${name} ${res} PARENT_SCOPE)
endfunction ()

getuname (osname -s)
getuname (osrel -r)
getuname (cpu -m)

set(CTEST_BUILD_NAME "${osname}-${osrel}-${cpu}")

## -- Git command
find_program (CTEST_GIT_COMMAND NAMES git)

## -- make command
find_program (MAKE NAMES make)

#-----------------------------------------------------------  
#-- Get build-specific information
#-----------------------------------------------------------  

## -- CTest Dashboard Root Directory
if (DEFINED ENV{PIO_DASHBOARD_ROOT})
    set (CTEST_DASHBOARD_ROOT "$ENV{PIO_DASHBOARD_ROOT}")
else ()
    set (CTEST_DASHBOARD_ROOT "$ENV{HOME}/pio-dashboard")
endif ()

## -- SRC Dir (where this script exists)
set (CTEST_SOURCE_DIRECTORY   "${CTEST_SCRIPT_DIRECTORY}")

## -- BIN Dir (in-source build)  
set (CTEST_BINARY_DIRECTORY   "${CTEST_DASHBOARD_ROOT}/build-${CTEST_BUILD_NAME}")

## -- Empty the binary directory
ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})

## -- Add the CTest script directory to the module path
set (CTEST_EXTRA_SCRIPT_PATH "${CTEST_SOURCE_DIRECTORY}/ctest")
set (CTEST_ENVIRONMENT_SCRIPT "CTestEnvironment-${CTEST_SITE_ID}")
set (CTEST_RUNCTEST_SCRIPT "${CTEST_EXTRA_SCRIPT_PATH}/runctest-${CTEST_SITE_ID}.sh")
list (APPEND CMAKE_MODULE_PATH ${CTEST_EXTRA_SCRIPT_PATH})

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
execute_process (COMMAND ${CTEST_RUNCTEST_SCRIPT} ${CTEST_SCRIPT_ARG}
                 WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY})

## -- SUBMIT
message (" -- Submit to dashboard - ${CTEST_BUILD_NAME} --")
ctest_submit ()

message (" -- Finished - ${CTEST_BUILD_NAME} --")
