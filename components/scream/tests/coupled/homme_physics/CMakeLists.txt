include (ScreamUtils)

# Needed for RRTMGP
find_library(NETCDF_C netcdf HINTS ${NetCDF_C_PATHS}/lib)
## Add CUDA flags for YAKL
if (CUDA_BUILD)
    set(ARCH "CUDA")
    set(YAKL_CXX_FLAGS "-D__USE_CUDA__ --expt-extended-lambda --expt-relaxed-constexpr ${CUDA_FLAGS} ${YAKL_CXX_FLAGS}")
endif()

# Get or create the dynamics lib
#                 HOMME_TARGET   NP PLEV QSIZE_D
CreateDynamicsLib("theta-l_kokkos"  4   72   35)

set (NEED_LIBS cld_fraction shoc p3 scream_rrtmgp rrtmgp ${NETCDF_C} ${dynLibName} scream_control scream_share physics_share yakl)

# Add source files
set (SRC shoc_cld_p3_rrtmgp.cpp ${SCREAM_BASE_DIR}/src/physics/rrtmgp/tests/rrtmgp_test_utils.cpp)

CreateUnitTest(homme_shoc_cld_p3_rrtmgp "homme_shoc_cld_p3_rrtmgp.cpp" "${NEED_LIBS}" EXE_ARGS "< namelist.nl" LABELS "dynamics;driver;shoc;cld;p3;rrtmgp;physics")
set_target_properties(homme_shoc_cld_p3_rrtmgp PROPERTIES COMPILE_FLAGS "${YAKL_CXX_FLAGS}")

# Copy yaml input file to run directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml COPYONLY)

# Set homme's test options, so that we can configure the namelist correctly

# Discretization/algorithm settings
set (HOMME_TEST_NE 2)
set (HOMME_TEST_LIM 9)
set (HOMME_TEST_QSIZE 1)
set (HOMME_TEST_RSPLIT 3)
set (HOMME_TEST_NMAX 24)
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
set (HOMME_TEST_MOISTURE notdry)
SET (HOMME_THETA_HY_MODE true)

# Vert coord settings
set (HOMME_TEST_VCOORD_INT_FILE acme-72i.ascii)
set (HOMME_TEST_VCOORD_MID_FILE acme-72m.ascii)

# Configure the namelist into the test directory
configure_file(${SCREAM_SOURCE_DIR}/src/dynamics/homme/tests/theta-nmax.nl
               ${CMAKE_CURRENT_BINARY_DIR}/namelist.nl)

# Copy initial conditions files
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/homme_shoc_cld_p3_rad_init_ne2np4.nc
               ${CMAKE_CURRENT_BINARY_DIR}/homme_shoc_cld_p3_rad_init_ne2np4.nc COPYONLY)

# Create vcoord directory
file (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/vcoord)

# Copy vcoord files to build dir
configure_file(${HOMME_SOURCE_DIR}/test/vcoord/acme-72i.ascii
               ${CMAKE_CURRENT_BINARY_DIR}/vcoord COPYONLY)
configure_file(${HOMME_SOURCE_DIR}/test/vcoord/acme-72m.ascii
               ${CMAKE_CURRENT_BINARY_DIR}/vcoord COPYONLY)

# P3 data files
CONFIGURE_FILE(${SCREAM_DATA_DIR}/p3_lookup_table_1.dat-v4.1.1
               ${CMAKE_CURRENT_BINARY_DIR}/data/p3_lookup_table_1.dat-v4.1.1 COPYONLY)
CONFIGURE_FILE(${SCREAM_DATA_DIR}/p3_lookup_table_2.dat-v4.1.1
               ${CMAKE_CURRENT_BINARY_DIR}/data/p3_lookup_table_2.dat-v4.1.1 COPYONLY)
# RRTMGP data files
CONFIGURE_FILE(${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
CONFIGURE_FILE(${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/rrtmgp/data/rrtmgp-data-lw-g256-2018-12-04.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
CONFIGURE_FILE(${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-sw.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
CONFIGURE_FILE(${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
