include (ScreamUtils)

# Get or create the dynamics lib
#                 HOMME_TARGET   NP PLEV QSIZE_D
CreateDynamicsLib("theta-l_kokkos"  4   72   35)

set (NEED_LIBS ${dynLibName} scream_control scream_share)

# Test Homme standalone
CreateUnitTest(homme_stand_alone "homme_stand_alone.cpp" "${NEED_LIBS}" EXE_ARGS "< namelist.nl" LABELS "dynamics;driver")

# Copy yaml input file to run directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml COPYONLY)

# Set homme's test options, so that we can configure the namelist correctly

# Discretization/algorithm settings
set (HOMME_TEST_NE 2)
set (HOMME_TEST_LIM 9)
set (HOMME_TEST_QSIZE 6)
set (HOMME_TEST_RSPLIT 3)
set (HOMME_TEST_NDAYS 1)
set (HOMME_TEST_TIME_STEP 600)
SET (HOMME_THETA_FORM 1)
SET (HOMME_TTYPE 5)

# Hyperviscosity settings
set (HOMME_TEST_HVSCALING 0)
SET (HOMME_TEST_HVS 1)
SET (HOMME_TEST_HVS_TOM 0)

set (HOMME_TEST_NU 7e15)
set (HOMME_TEST_NUDIV 1e15)
set (HOMME_TEST_NUTOP 2.5e5)

# Testcase settings
set (HOMME_TEST_MOISTURE dry)
SET (HOMME_THETA_HY_MODE true)  

# Vert coord settings
set (HOMME_TEST_VCOORD_INT_FILE acme-72i.ascii)
set (HOMME_TEST_VCOORD_MID_FILE acme-72m.ascii)

# Configure the namelist into the test directory
configure_file(${HOMME_SOURCE_DIR}/test/reg_test/namelists/theta.nl
               ${CMAKE_CURRENT_BINARY_DIR}/namelist.nl)

# Create vcoord directory
file (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/vcoord)

# Copy vcoord files to build dir
configure_file(${HOMME_SOURCE_DIR}/test/vcoord/acme-72i.ascii
               ${CMAKE_CURRENT_BINARY_DIR}/vcoord COPYONLY)
configure_file(${HOMME_SOURCE_DIR}/test/vcoord/acme-72m.ascii
               ${CMAKE_CURRENT_BINARY_DIR}/vcoord COPYONLY)
