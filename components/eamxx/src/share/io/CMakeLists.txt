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

# TODO: Until we are done with the refactor of eamxx libs, we still have some
# cyclic deps. Linking scream_io to eamxx_test_support seems to fix things for now
target_link_libraries(eamxx_test_support PUBLIC scream_io)
