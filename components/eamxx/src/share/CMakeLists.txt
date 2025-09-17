# Core types and configuration options
add_subdirectory(core)

# Small/generic utilities
add_subdirectory(util)

# Field structures
add_subdirectory(field)

set(SHARE_SRC
  atm_process/atmosphere_process.cpp
  atm_process/atmosphere_process_hash.cpp
  atm_process/atmosphere_process_group.cpp
  atm_process/atmosphere_process_dag.cpp
  atm_process/atmosphere_diagnostic.cpp
  atm_process/IOPDataManager.cpp
  grid/abstract_grid.cpp
  grid/grids_manager.cpp
  grid/grid_import_export.cpp
  grid/se_grid.cpp
  grid/point_grid.cpp
  grid/remap/abstract_remapper.cpp
  grid/remap/coarsening_remapper.cpp
  grid/remap/horiz_interp_remapper_base.cpp
  grid/remap/horiz_interp_remapper_data.cpp
  grid/remap/iop_remapper.cpp
  grid/remap/identity_remapper.cpp
  grid/remap/refining_remapper_p2p.cpp
  grid/remap/vertical_remapper.cpp
  manager/field_manager.cpp
  property_checks/property_check.cpp
  property_checks/field_nan_check.cpp
  property_checks/field_within_interval_check.cpp
  property_checks/mass_and_energy_conservation_check.cpp
  algorithm/eamxx_data_interpolation.cpp
  algorithm/eamxx_fv_phys_rrtmgp_active_gases_workaround.cpp
  algorithm/eamxx_time_interpolation.cpp
)

if (EAMXX_ENABLE_EXPERIMENTAL_CODE)
  list (APPEND
    SHARE_SRC
    grid/remap/refining_remapper_rma.cpp
  )
endif()

add_library(scream_share ${SHARE_SRC})

target_include_directories(scream_share PUBLIC
  ${SCREAM_SRC_DIR}
  ${SCREAM_BIN_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(scream_share PUBLIC
  eamxx_core
  eamxx_utils
  eamxx_field)

# Used in the data interpolation
target_link_libraries(scream_share PUBLIC stdc++fs)

target_link_libraries(scream_share PUBLIC ekat::AllLibs pioc)

# The "build_cf_dictionary" target downloads an XML file containing valid field
# names and aliases based on the CF conventions, and transforms it into a YAML
# file that we use to validate our own fields against it. This target is only
# available if CMake finds wget on your system.
find_program(WGET wget)
if (NOT WGET STREQUAL "WGET-NOTFOUND")
  set(CF_CONVENTIONS_XML_FILE cf-standard-name-table.xml)
  set(CF_CONVENTIONS_URL "http://cfconventions.org/Data/cf-standard-names/76/src/${CF_CONVENTIONS_XML_FILE}")
  add_custom_target(build_cf_dictionary
    ${WGET} ${CF_CONVENTIONS_URL}
    COMMAND ${PROJECT_SOURCE_DIR}/scripts/cf-xml-to-yaml ${CF_CONVENTIONS_XML_FILE}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/data
    COMMENT "Downloading CF conventions XML file and building YAML dictionary...")
endif()

# IO library
add_subdirectory(io)

if (NOT SCREAM_LIB_ONLY)
  # Create test_support lib
  add_library(scream_test_support
    grid/mesh_free_grids_manager.cpp
  )
  target_link_libraries(scream_test_support PUBLIC scream_share scream_io)
  target_include_directories (scream_test_support PUBLIC ${EKAT_SOURCE_DIR}/extern/Catch2/single_include)

  add_subdirectory(tests)
  add_subdirectory(io/tests)
endif()
