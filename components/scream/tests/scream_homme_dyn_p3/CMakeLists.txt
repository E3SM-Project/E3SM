INCLUDE (ScreamUtils)

# Get or create the dynamics lib
#                 HOMME_TARGET   NP PLEV QSIZE_D USE_PIO USE_ENERGY 
CreateDynamicsLib("theta-l_kokkos"  4   72   35    FALSE   FALSE)

SET (NEED_LIBS p3 ${dynLibName} scream_control scream_share)

# Test atmosphere processes
CreateUnitTest(homme_p3 "homme_p3.cpp" "${NEED_LIBS}" EXE_ARGS "< namelist.nl" LABELS "dynamics;p3;driver")

# Copy yaml input file to run directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml COPYONLY)

# Set homme's test options, so that we can configure the namelist correctly
set (HOMME_TEST_NE 2)
set (HOMME_TEST_NDAYS 1)
set (HOMME_TEST_LIM 9)
set (HOMME_TEST_QSIZE 6)
set (HOMME_TEST_RSPLIT 3)
set (HOMME_TEST_MOISTURE dry)
set (HOMME_TEST_HVSCALING 0)
set (HOMME_TEST_NU 7e15)
set (HOMME_TEST_NUDIV 1e15)
set (HOMME_TEST_NUTOP 2.5e5)
set (HOMME_TEST_TIME_STEP 600)
set (HOMME_TEST_VCOORD_INT_FILE acme-72i.ascii)
set (HOMME_TEST_VCOORD_MID_FILE acme-72m.ascii)

# Configure the namelist into the test directory
CONFIGURE_FILE(${HOMME_SOURCE_DIR}/test/reg_test/namelists/theta.nl
               ${CMAKE_CURRENT_BINARY_DIR}/namelist.nl)

# Create movies and vcoord directory
FILE (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/movies)
FILE (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/vcoord)

# Copy vcoord files to build dir
CONFIGURE_FILE(${HOMME_SOURCE_DIR}/test/vcoord/acme-72i.ascii
               ${CMAKE_CURRENT_BINARY_DIR}/vcoord COPYONLY)
CONFIGURE_FILE(${HOMME_SOURCE_DIR}/test/vcoord/acme-72m.ascii
               ${CMAKE_CURRENT_BINARY_DIR}/vcoord COPYONLY)

# Copy p3 lookup tables to local data directory
FILE (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data)
CONFIGURE_FILE(${SCREAM_DATA_DIR}/p3_lookup_table_1.dat-v4
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
CONFIGURE_FILE(${SCREAM_DATA_DIR}/p3_lookup_table_2.dat-v4
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
CONFIGURE_FILE(${SCREAM_DATA_DIR}/p3_universal_constants.inp
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
