#include(ScreamUtils)
include(ScreamUtils)

# include(CsmShare)
# CreateCsmShareTarget()

set (ARRAY_SCORPIO_SRCS
  scorpio_tests.cpp
  scorpio_yaml_test.cpp
  scorpio_ad_test.cpp
  io.cpp
  restart.cpp
  ${SCREAM_BASE_DIR}/../cam/src/physics/cam/physics_utils.F90
)
set (SCREAM_LIBS
     scream_share
     scream_control
     scream_io
     p3
)

set(CSM_SHARE "${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/${NINST_VALUE}/csm_share")
set (SCREAM_INCLUDE
     ${CMAKE_BINARY_DIR}/cmake/atm
     ${INSTALL_SHAREDPATH}/include
     ${SCREAM_TPL_INCLUDE_DIRS}
     ${CSM_SHARE}
)

# Test atmosphere processes
CreateUnitTest(scorpio_interface "scorpio_tests.cpp" "${SCREAM_LIBS}"
  INCLUDE_DIRS ${SCREAM_INCLUDE}
  MPI_RANKS 1 4
)
CreateUnitTest(scorpio_yaml "scorpio_yaml_test.cpp" "${SCREAM_LIBS}"
  INCLUDE_DIRS ${SCREAM_INCLUDE}
  MPI_RANKS 1 4
)
CreateUnitTest(scorpio_ad "scorpio_ad_test.cpp" "${SCREAM_LIBS}"
  INCLUDE_DIRS ${SCREAM_INCLUDE}
  MPI_RANKS 1 4
)

CreateUnitTest(io_test "io.cpp" "${SCREAM_LIBS}"
  INCLUDE_DIRS ${SCREAM_INCLUDE}
)

CreateUnitTest(restart_test "restart.cpp" "${SCREAM_LIBS}"
  INCLUDE_DIRS ${SCREAM_INCLUDE}
)

configure_file(scorpio_output_baseline.nc scorpio_output_baseline.nc COPYONLY)
configure_file(scorpio_control.yaml scorpio_control.yaml COPYONLY)
configure_file(scorpio_output1.yaml scorpio_output1.yaml COPYONLY)
configure_file(scorpio_output2.yaml scorpio_output2.yaml COPYONLY)
configure_file(scorpio_ad_test.yaml scorpio_ad_test.yaml COPYONLY)
configure_file(scorpio_ad_out_ins.yaml scorpio_ad_out_ins.yaml COPYONLY)
configure_file(scorpio_ad_out_avg.yaml scorpio_ad_out_avg.yaml COPYONLY)

configure_file(io_test_instant.yaml io_test_instant.yaml COPYONLY)
configure_file(io_test_average.yaml io_test_average.yaml COPYONLY)
configure_file(io_test_max.yaml io_test_max.yaml COPYONLY)
configure_file(io_test_min.yaml io_test_min.yaml COPYONLY)
configure_file(io_test_restart.yaml io_test_restart.yaml COPYONLY)
# Copy p3 lookup tables to local data directory
FILE (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data)
CONFIGURE_FILE(${SCREAM_DATA_DIR}/p3_lookup_table_1.dat-v4
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
CONFIGURE_FILE(${SCREAM_DATA_DIR}/p3_lookup_table_2.dat-v4
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
CONFIGURE_FILE(${SCREAM_DATA_DIR}/p3_universal_constants.inp
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
