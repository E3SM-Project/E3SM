# Example originally stolen from:
# http://www.vtk.org/Wiki/CTest:Using_CTEST_and_CDASH_without_CMAKE

# -----------------------------------------------------------  
# -- Get the machine environment
# -----------------------------------------------------------  

## -- Set hostname
## --------------------------

find_program (HOSTNAME_CMD NAMES hostname)
execute_process (COMMAND ${HOSTNAME_CMD}
                 OUTPUT_VARIABLE HOSTNAME)
                 
set (CTEST_SITE "${HOSTNAME}")

## -- Set hostname ID (e.g., alcf, nwsc, nersc, ...)
## --------------------------------------------------

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
## --------------------------

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
## ----------------
find_program (CTEST_GIT_COMMAND NAMES git)

## -- make command
## -----------------
find_program (MAKE NAMES make)

# -----------------------------------------------------------  
# -- build specific
# -----------------------------------------------------------  

## -- Get the CTest Dashboard Model
set (MODEL Experimental)
if (${CTEST_SCRIPT_ARG} MATCHES Nightly)
  set (MODEL Nightly)
endif ()

## -- Dashboard Root Dir
if (DEFINED ENV{PIO_DASHBOARD_ROOT})
    set (CTEST_DASHBOARD_ROOT "$ENV{PIO_DASHBOARD_ROOT}")
else ()
    set (CTEST_DASHBOARD_ROOT "$ENV{HOME}/pio-dashboard")
endif ()

## -- SRC Dir
set (CTEST_SOURCE_DIRECTORY   "${CTEST_DASHBOARD_ROOT}/src")

## -- BIN Dir                                            
set (CTEST_BINARY_DIRECTORY   "${CTEST_DASHBOARD_ROOT}/build-${CTEST_BUILD_NAME}")

## -- Empty the binary directory
ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})

# -----------------------------------------------------------  
# -- CTest Step Commands
# -----------------------------------------------------------  

## -- Checkout command

# -----------------------------------------------------------  
# -- CTest Step Commands
# -----------------------------------------------------------  

## -- Checkout command
if(NOT EXISTS ${CTEST_SOURCE_DIRECTORY})
    set (CTEST_GIT_URL "https://github.com/PARALLELIO/ParallelIO")
    set (CTEST_CHECKOUT_COMMAND "${CTEST_GIT_COMMAND} clone ${CTEST_GIT_URL} ${CTEST_SOURCE_DIRECTORY}")
endif()

## -- Update Command
set (CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")

## -- Configure Command
set (CTEST_CONFIGURE_COMMAND "${CMAKE_COMMAND} ${CTEST_BUILD_OPTIONS} ${CTEST_SOURCE_DIRECTORY}")

## -- Build Command
set (CTEST_BUILD_COMMAND "${MAKE} all tests")

# -----------------------------------------------------------  
# -- Run CTest
# -----------------------------------------------------------  

## -- Start
message (" -- Start dashboard ${MODEL} - ${CTEST_BUILD_NAME} --")
ctest_start("${MODEL}")

## -- Update
message (" -- Update ${MODEL} - ${CTEST_BUILD_NAME} --")
ctest_update ()

## -- Configure 
message (" -- Configure ${MODEL} - ${CTEST_BUILD_NAME} --")
ctest_configure ()

## -- BUILD
message (" -- Build ${MODEL} - ${CTEST_BUILD_NAME} --")
ctest_build ()

## -- TEST
message (" -- Test ${MODEL} - ${CTEST_BUILD_NAME} --")
#ctest_test ()
execute_process (COMMAND which ctest
                 OUTPUT_VARIABLE MY_CTEST_VAR)
message (" ***** CTest is ${MY_CTEST_VAR} *****")

## -- SUBMIT
#message (" -- Submit ${MODEL} - ${CTEST_BUILD_NAME} --")
#ctest_submit ()

message (" -- Finished ${MODEL}  - ${CTEST_BUILD_NAME} --")
