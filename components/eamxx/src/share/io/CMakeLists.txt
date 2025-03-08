#########################################
#       SCORPIO interface library       #
#########################################

# This small lib contains some interfaces to scorpio, which are formulated
# in terms of names rather than IDs. Internally, we store a list of
# small structs, holding handlers to PIO data and extra metadata, which
# allows to perform sanity checks and print more helpful messages

add_library(eamxx_scorpio_interface
  eamxx_scorpio_types.cpp
  eamxx_scorpio_interface.cpp
)
target_link_libraries(eamxx_scorpio_interface PUBLIC ekat)
target_link_libraries(eamxx_scorpio_interface PRIVATE pioc)
target_include_directories(eamxx_scorpio_interface PUBLIC
  ${SCREAM_BIN_DIR}/src   # For eamxx_config.h
)

if (DEFINED ENV{ADIOS2_ROOT})
  target_include_directories(eamxx_scorpio_interface PRIVATE $ENV{ADIOS2_ROOT}/include)
endif ()

if (SCORPIO_Fortran_INCLUDE_DIRS)
  target_include_directories(eamxx_scorpio_interface PUBLIC ${SCORPIO_Fortran_INCLUDE_DIRS})
endif ()

if (SCORPIO_C_INCLUDE_DIRS)
  target_include_directories(eamxx_scorpio_interface PUBLIC ${SCORPIO_C_INCLUDE_DIRS})
endif ()

if (SCREAM_CIME_BUILD)
  # Add interface to E3SM shr lib (to retrieve PIO subsystem info)
  target_sources (eamxx_scorpio_interface PRIVATE
    eamxx_shr_interface_c2f.F90
  )
  target_compile_definitions (eamxx_scorpio_interface PRIVATE SCREAM_CIME_BUILD)
  set_target_properties(eamxx_scorpio_interface PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
  )
  target_include_directories(eamxx_scorpio_interface PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/modules)
  target_link_libraries(eamxx_scorpio_interface PUBLIC csm_share)
endif()

##################################
#       EAMxx I/O library        #
##################################

# This library allows for a simple(r) interaction between EAMxx
# data structures (such as grids, fields, remappers,...) and
# the scorpio interface library

# Create io lib
add_library(scream_io
  eamxx_output_manager.cpp
  scorpio_input.cpp
  scorpio_scm_input.cpp
  scorpio_output.cpp
  eamxx_io_utils.cpp
)

target_link_libraries(scream_io PUBLIC scream_share eamxx_scorpio_interface diagnostics)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
