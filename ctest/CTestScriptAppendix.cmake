#==============================================================================
#
#  This is the CTest script for builds and submission to the CTest
#  Dashboard site: my.cdash.org.
#
#  This test is used to "append" to an existing tag (ctest run)
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
set (CTEST_SOURCE_DIRECTORY   "${CTEST_SCRIPT_DIRECTORY}/..")

## -- BIN Dir (in-source build)  
set (CTEST_BINARY_DIRECTORY   "${CTEST_DASHBOARD_ROOT}/build-${CTEST_BUILD_NAME}")

## -- Add the CTest script directory to the module path
set (CTEST_EXTRA_SCRIPT_PATH "${CTEST_SOURCE_DIRECTORY}/ctest")
set (CTEST_ENVIRONMENT_SCRIPT "CTestEnvironment-${CTEST_SITE_ID}")
set (CTEST_RUNCTEST_SCRIPT "${CTEST_EXTRA_SCRIPT_PATH}/runctest-${CTEST_SITE_ID}.sh")
list (APPEND CMAKE_MODULE_PATH ${CTEST_EXTRA_SCRIPT_PATH})

# -----------------------------------------------------------  
# -- Run CTest- TESTING ONLY
# -----------------------------------------------------------  

## -- Start
ctest_start("${CTEST_SCRIPT_ARG}" APPEND)

## -- TEST
execute_process (COMMAND ${CTEST_RUNCTEST_SCRIPT} ${CTEST_SCRIPT_ARG}
                 WORKING_DIRECTORY ${CTEST_BINARY_DIRECTORY})
