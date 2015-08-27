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

set(MODEL "Experimental")

## -- Dashboard Root Dir
if (DEFINED ENV{PIO_DASHBOARD_ROOT})
    set (CTEST_DASHBOARD_ROOT "$ENV{PIO_DASHBOARD_ROOT}")
else ()
    set (CTEST_DASHBOARD_ROOT "$ENV{HOME}/pio-dashboard")
endif ()

## -- SRC Dir
set (CTEST_SOURCE_DIRECTORY   "$ENV{HOME}/Development/Workspace/git/ParallelIO")

## -- BIN Dir                                            
set (CTEST_BINARY_DIRECTORY   "${CTEST_DASHBOARD_ROOT}/build-${CTEST_BUILD_NAME}")

## -- Build options
set (CTEST_BUILD_OPTIONS      "")

ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})

# -----------------------------------------------------------  
# -- commands
# -----------------------------------------------------------  

## -- Checkout command
if(NOT EXISTS ${CTEST_SOURCE_DIRECTORY})
    #set (CTEST_GIT_URL "https://github.com/PARALLELIO/ParallelIO")
    set (CTEST_GIT_URL "file://$ENV{HOME}/Development/Workspace/git/ParallelIO")
    set (CTEST_CHECKOUT_COMMAND "${CTEST_GIT_COMMAND} clone ${CTEST_GIT_URL} ${CTEST_SOURCE_DIRECTORY}")
endif()

## -- Update Command
set(CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")

## -- Configure Command
set(CTEST_CONFIGURE_COMMAND "${CMAKE_COMMAND} -DCMAKE_VERBOSE_MAKEFILE=1 ${CTEST_SOURCE_DIRECTORY}")

## -- Build Command
set(CTEST_BUILD_COMMAND "${MAKE} ${CTEST_BUILD_OPTIONS} all tests")

# -----------------------------------------------------------  
# -- Configure CTest
# -----------------------------------------------------------  

## -- CTest Config
#configure_file($ENV{HOME}/CTestConfiguration/CTestConfig.cmake  ${CTEST_SOURCE_DIRECTORY}/CTestConfig.cmake)

## -- CTest Custom
#configure_file($ENV{HOME}/CTestConfiguration/CTestCustom.cmake ${CTEST_BINARY_DIRECTORY}/CTestCustom.cmake)

## -- CTest Testfile
#configure_file($ENV{HOME}/CTestConfiguration/CTestTestfile.cmake ${CTEST_BINARY_DIRECTORY}/CTestTestfile.cmake)

## -- read CTestCustom.cmake file
#ctest_read_custom_files("${CTEST_BINARY_DIRECTORY}")

# -----------------------------------------------------------  
# -- Run CTest
# -----------------------------------------------------------  

## -- Start
message (" -- Start dashboard ${MODEL} - ${CTEST_BUILD_NAME} --")
ctest_start("${MODEL}")

## -- Update
#message (" -- Update ${MODEL} - ${CTEST_BUILD_NAME} --")
#ctest_update ()

## -- Configure 
message (" -- Configure ${MODEL} - ${CTEST_BUILD_NAME} --")
ctest_configure ()

## -- BUILD
message (" -- Build ${MODEL} - ${CTEST_BUILD_NAME} --")
ctest_build ()

## -- TEST
message (" -- Test ${MODEL} - ${CTEST_BUILD_NAME} --")
ctest_test ()

## -- SUBMIT
message (" -- Submit ${MODEL} - ${CTEST_BUILD_NAME} --")
ctest_submit ()

message (" -- Finished ${MODEL}  - ${CTEST_BUILD_NAME} --")
