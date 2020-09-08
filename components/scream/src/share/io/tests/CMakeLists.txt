#include(ScreamUtils)
include(ScreamUtils)

include(CsmShare)
CreateCsmShareTarget()

set (ARRAY_SCORPIO_SRCS
  scorpio_tests.cpp
  scorpio_yaml_test.cpp
  ${SCREAM_BASE_DIR}/../cam/src/physics/cam/physics_utils.F90
)
set (SCREAM_LIBS
     scream_share
     scream_control
     io
     scream_piof
     scream_csm_share
     scream_pioc
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

configure_file(${SCREAM_DATA_DIR}/scorpio_output_baseline.nc scorpio_output_baseline.nc COPYONLY)
configure_file(${SCREAM_DATA_DIR}/scorpio_control.yaml scorpio_control.yaml COPYONLY)
configure_file(${SCREAM_DATA_DIR}/scorpio_output1.yaml scorpio_output1.yaml COPYONLY)
configure_file(${SCREAM_DATA_DIR}/scorpio_output2.yaml scorpio_output2.yaml COPYONLY)
