add_library(eamxx_atm_process
  atmosphere_process.cpp
  atmosphere_process_hash.cpp
  atmosphere_process_group.cpp
  atmosphere_process_dag.cpp
  atmosphere_diagnostic.cpp
)

target_link_libraries(eamxx_atm_process PUBLIC
  eamxx_core
  eamxx_utils
  eamxx_field
  eamxx_property_checks
  eamxx_data_managers)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
