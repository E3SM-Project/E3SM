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
message("ASD")
message("${ARRAY_SCORPIO_SRCS}")
message("${SCREAM_LIBS}")
message("${SCREAM_INCLUDE}")
# Test atmosphere processes
EkatCreateUnitTest(scorpio_interface "${ARRAY_SCORPIO_SRCS}" 
  LIBS ${SCREAM_LIBS} 
  LIBS_DIRS ${SCREAM_INCLUDE}
  INCLUDE_DIRS ${SCREAM_INCLUDE}
)
