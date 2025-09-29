add_library(eamxx_utils
  eamxx_bfbhash.cpp
  eamxx_time_stamp.cpp
  eamxx_timing.cpp
  eamxx_repro_sum_mod.F90
  eamxx_utils.cpp
)

set_target_properties(eamxx_utils PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
target_include_directories (eamxx_utils PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_link_libraries (eamxx_utils PUBLIC
  eamxx_core
  csm_share # for shr_reprosum_mod
  ekat::Logging
  ekat::Pack)

target_link_libraries (eamxx_utils PUBLIC gptl)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
