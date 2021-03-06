set(SCREAM_SCORPIO_SRCS
  scorpio_input.cpp
  scream_scorpio_interface.F90
  scream_scorpio_interface.cpp
  scream_scorpio_interface_iso_c2f.F90
  scorpio_input.cpp
  scorpio_output.cpp
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
set_target_properties(scream_io PROPERTIES Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})
target_link_libraries(scream_io PUBLIC scream_share ${SCREAM_CIME_LIBS})
target_compile_options(scream_io PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${SCREAM_Fortran_FLAGS}>)

if (NOT SCREAM_LIBS_ONLY)
  add_subdirectory(tests)
endif()
