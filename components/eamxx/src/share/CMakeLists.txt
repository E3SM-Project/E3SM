include (EkatSetCompilerFlags)

set(SHARE_SRC
  scream_config.cpp
  scream_session.cpp
  atm_process/atmosphere_process.cpp
  atm_process/atmosphere_process_group.cpp
  atm_process/atmosphere_process_dag.cpp
  atm_process/atmosphere_diagnostic.cpp
  field/field_alloc_prop.cpp
  field/field_identifier.cpp
  field/field_header.cpp
  field/field_layout.cpp
  field/field_tracking.cpp
  field/field.cpp
  field/field_group.cpp
  field/field_manager.cpp
  grid/abstract_grid.cpp
  grid/grids_manager.cpp
  grid/se_grid.cpp
  grid/point_grid.cpp
  grid/remap/abstract_remapper.cpp
  grid/remap/coarsening_remapper.cpp
  grid/remap/vertical_remapper.cpp
  grid/remap/horizontal_remap_utility.cpp
  property_checks/property_check.cpp
  property_checks/field_nan_check.cpp
  property_checks/field_within_interval_check.cpp
  property_checks/mass_and_energy_column_conservation_check.cpp
  util/scream_time_stamp.cpp
  util/scream_timing.cpp
  util/scream_utils.cpp
)

add_library(scream_share ${SHARE_SRC})
set_target_properties(scream_share PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
)
target_include_directories(scream_share PUBLIC
  ${SCREAM_SRC_DIR}
  ${SCREAM_BIN_DIR}/src
  ${CMAKE_CURRENT_BINARY_DIR}/modules
)
target_link_libraries(scream_share PUBLIC ekat gptl)
target_compile_options(scream_share PUBLIC
  $<$<COMPILE_LANGUAGE:Fortran>:${SCREAM_Fortran_FLAGS}>
)

# This used to be in `scream_config.h/scream_config.f`, but we did accidentally remove
# the includes at least twice. To avoid this from happening again, we add this crucial
# macro as a compile definition of the scream_share, which is linked to all other
# scream libraries.
if (SCREAM_CIME_BUILD)
  target_compile_definitions (scream_share PUBLIC SCREAM_CIME_BUILD)
endif()
# We have some issues with RDC and repeated libraries in the link line.
# It's a known issue, and nvcc_wrapper has a flag for handling this.
if (Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE)
  if(CMAKE_VERSION VERSION_LESS "3.13.0")
    target_link_librarires(scream_share PUBLIC "--remove-duplicate-link-files")
  else()
    target_link_options(scream_share PUBLIC "--remove-duplicate-link-files")
  endif()
endif()
SetOmpFlags(scream_share Fortran)

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
    util/scream_test_session.cpp
  )
  target_link_libraries(scream_test_support PUBLIC scream_share scream_io)

  add_subdirectory(tests)
endif()
