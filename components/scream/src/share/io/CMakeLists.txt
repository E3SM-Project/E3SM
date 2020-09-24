set(SCREAM_SCORPIO_SRCS
  scream_scorpio_interface.F90
  scream_scorpio_interface.cpp
  scream_scorpio_interface_iso_c2f.F90
  ${SCREAM_BASE_DIR}/../cam/src/physics/cam/physics_utils.F90
)

set (SCREAM_LIBS
     scream_share
     scream_control
)

# Create or import scorpio targets
include (${SCREAM_SOURCE_DIR}/cmake/tpls/Scorpio.cmake)
CreateScorpioTarget(TRUE) # Note: the input is whether we want FLIB too. We want/need both, so look for both

set (SCREAM_CIME_LIBS
     pioc
     piof
)

# Create io lib
add_library(scream_io ${SCREAM_SCORPIO_SRCS})

# Set f90 compiled modules dir
set_target_properties(scream_io PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Link libraries
target_link_libraries(scream_io PRIVATE ${SCREAM_LIBS} ${SCREAM_CIME_LIBS})

if (NOT SCREAM_LIBS_ONLY)
  add_subdirectory(tests)
endif()
