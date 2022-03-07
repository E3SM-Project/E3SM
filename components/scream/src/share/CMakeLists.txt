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
  grid/mesh_free_grids_manager.cpp
  grid/se_grid.cpp
  grid/point_grid.cpp
  grid/user_provided_grids_manager.cpp
  grid/remap/abstract_remapper.cpp
  property_checks/property_check.cpp
  property_checks/field_nan_check.cpp
  property_checks/field_within_interval_check.cpp
  util/scream_test_session.cpp
  util/scream_time_stamp.cpp
)

add_library(scream_share ${SHARE_SRC})
target_include_directories(scream_share PUBLIC ${SCREAM_SRC_DIR} ${SCREAM_BIN_DIR}/src)
target_link_libraries(scream_share PUBLIC ekat)
set_target_properties(scream_share PROPERTIES Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})
target_compile_options(scream_share PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:${SCREAM_Fortran_FLAGS}>)
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

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()

add_subdirectory(io)
