#include(ScreamUtils)
include(EkatCreateUnitTest)

set (ARRAY_SCORPIO_SRCS
  scorpio_tests.cpp
  ${SCREAM_BASE_DIR}/../cam/src/physics/cam/physics_utils.F90
)
set (SCREAM_LIBS
     scream_share
     scream_control
     atm
     ${SCREAM_TPL_LIBRARIES}
)

set(CSM_SHARE "${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/${NINST_VALUE}/csm_share")
set (SCREAM_INCLUDE
     ${CMAKE_BINARY_DIR}/cmake/atm
     ${INSTALL_SHAREDPATH}/include
     ${SCREAM_TPL_INCLUDE_DIRS}
     ${CSM_SHARE}
)
# Test atmosphere processes
EkatCreateUnitTest(scorpio_interface "${ARRAY_SCORPIO_SRCS}" 
  LIBS ${SCREAM_LIBS} csm_share
  LIBS_DIRS ${SCREAM_INCLUDE} ${CSM_SHARE}
  INCLUDE_DIRS ${SCREAM_INCLUDE}
  THREADS 1 4
)

configure_file(${SCREAM_DATA_DIR}/scorpio_output_baseline.nc scorpio_output_baseline.nc COPYONLY)

