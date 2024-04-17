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
set_target_properties(scream_io PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
)
target_include_directories(scream_io PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/modules)

if (DEFINED ENV{ADIOS2_ROOT})
  target_include_directories(scream_io PRIVATE $ENV{ADIOS2_ROOT}/include)
endif ()

if (SCORPIO_Fortran_INCLUDE_DIRS)
  target_include_directories(scream_io PUBLIC ${SCORPIO_Fortran_INCLUDE_DIRS})
endif ()

if (SCORPIO_C_INCLUDE_DIRS)
  target_include_directories(scream_io PUBLIC ${SCORPIO_C_INCLUDE_DIRS})
endif ()

target_link_libraries(scream_io PUBLIC scream_share piof pioc)

if (SCREAM_CIME_BUILD)
  target_link_libraries(scream_io PUBLIC csm_share)
endif()

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
