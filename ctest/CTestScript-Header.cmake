#==============================================================================
#
#  This is the common CTest script header for builds and submission to the 
#  CTest Dashboard site: my.cdash.org.
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
                 OUTPUT_VARIABLE HOSTNAME
                 OUTPUT_STRIP_TRAILING_WHITESPACE)

## -- Set hostname ID (e.g., alcf, nwsc, nersc, ...)

# UCAR/NWSC Machines
if (HOSTNAME MATCHES "^yslogin" OR
    HOSTNAME MATCHES "^geyser" OR
    HOSTNAME MATCHES "^caldera")
    set (HOSTNAME_ID "nwsc")
# ALCF/Argonne Machines
elseif (HOSTNAME MATCHES "^mira" OR
        HOSTNAME MATCHES "^cetus" OR
        HOSTNAME MATCHES "^vesta" OR
        HOSTNAME MATCHES "^cooley")
    set (HOSTNAME_ID "alcf")
# ALCF/Argonne Machines
elseif (HOSTNAME MATCHES "^edison" OR
        HOSTNAME MATCHES "^carver" OR
        HOSTNAME MATCHES "^hopper")
    set (HOSTNAME_ID "nersc")
else ()
    set (HOSTNAME_ID "unknown")
endif ()

## -- Set the site name

set (CTEST_SITE "${HOSTNAME_ID}")

## -- Set the build name

find_program (UNAME NAMES uname)
function (getuname name flag)
    execute_process (COMMAND ${UNAME} ${flag}
                     OUTPUT_VARIABLE res
                     OUTPUT_STRIP_TRAILING_WHITESPACE)
    set (${name} ${res} PARENT_SCOPE)
endfunction ()

getuname (osname -s)
getuname (osrel -r)
getuname (cpu -m)

set(CTEST_BUILD_NAME "${HOSTNAME}-${osname}-${osrel}-${cpu}")

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

## -- BIN Dir 
set (CTEST_BINARY_DIRECTORY   "${CTEST_DASHBOARD_ROOT}/build")

