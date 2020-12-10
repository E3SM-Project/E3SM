#include(ScreamUtils)
include(ScreamUtils)

set (ARRAY_SCORPIO_SRCS
  scorpio_tests.cpp
  ${SCREAM_BASE_DIR}/../eam/src/physics/cam/physics_utils.F90
)
set (SCREAM_LIBS
     scream_share
     scream_control
     atm
)

# Test atmosphere processes
CreateUnitTest(scorpio_interface "${ARRAY_SCORPIO_SRCS}" "${SCREAM_LIBS}"
  MPI_RANKS 1 4
)

configure_file(${SCREAM_DATA_DIR}/scorpio_output_baseline.nc scorpio_output_baseline.nc COPYONLY)
