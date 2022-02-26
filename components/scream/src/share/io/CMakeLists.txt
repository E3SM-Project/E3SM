set(SCREAM_SCORPIO_SRCS
  scream_scorpio_interface.F90
  scream_scorpio_interface.cpp
  scream_scorpio_interface_iso_c2f.F90
  scream_output_manager.cpp
  scorpio_input.cpp
  scorpio_output.cpp
)

# Create or import scorpio targets
include (${SCREAM_BASE_DIR}/cmake/tpls/Scorpio.cmake)
CreateScorpioTargets()

set (SCREAM_CIME_LIBS
     pioc
     piof
)

# Create io lib
add_library(scream_io ${SCREAM_SCORPIO_SRCS})
set_target_properties(scream_io PROPERTIES Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})
target_link_libraries(scream_io PUBLIC scream_share ${SCREAM_CIME_LIBS})
target_compile_options(scream_io PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${SCREAM_Fortran_FLAGS}>)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
