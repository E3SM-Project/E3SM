##################################
#       EAMxx I/O library        #
##################################

# This library allows for a simple(r) interaction between EAMxx
# data structures (such as grids, fields, remappers,...) and
# the scorpio interface library

add_library(eamxx_io
  eamxx_output_manager.cpp
  scorpio_output.cpp
  eamxx_io_utils.cpp
  eamxx_diag_spec.cpp
  eamxx_diag_tree.cpp
  eamxx_diag_bank.cpp
)

target_link_libraries(eamxx_io PUBLIC
  eamxx_data_managers
  eamxx_scorpio_interface
  eamxx_diagnostics)

# NOTE: PRIVATE on purpose. e3sm_diags_parser requires C++20, while E3SM builds
#       at C++17, so a PUBLIC link would push cxx_std_20 onto every consumer of
#       eamxx_io. eamxx_diag_spec.hpp does not include any edp header, so the
#       dependency can stay contained in this library.
target_link_libraries(eamxx_io PRIVATE e3sm_diags_parser)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
