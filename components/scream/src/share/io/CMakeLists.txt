set(SCREAM_SCORPIO_SRCS
  scream_scorpio_interface.F90
  scream_scorpio_interface.cpp
  scream_scorpio_interface_iso_c2f.F90
  scream_output_manager.cpp
  scorpio_input.cpp
  scorpio_output.cpp
  scream_io_utils.cpp
)

# Create io lib
add_library(scream_io ${SCREAM_SCORPIO_SRCS})
set_target_properties(scream_io PROPERTIES Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})
target_link_libraries(scream_io PUBLIC scream_share diagnostics piof)
target_compile_options(scream_io PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${SCREAM_Fortran_FLAGS}>)

if (SCREAM_CIME_BUILD)
  target_link_libraries(scream_io PUBLIC csm_share)
endif()

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
