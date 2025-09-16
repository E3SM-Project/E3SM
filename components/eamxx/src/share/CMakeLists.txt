# Core types and configuration options
add_subdirectory(core)

# Small/generic utilities
add_subdirectory(util)

set(SHARE_SRC
  atm_process/atmosphere_process.cpp
  atm_process/atmosphere_process_hash.cpp
  atm_process/atmosphere_process_group.cpp
  atm_process/atmosphere_process_dag.cpp
  atm_process/atmosphere_diagnostic.cpp
  atm_process/IOPDataManager.cpp
  field/field_alloc_prop.cpp
  field/field_identifier.cpp
  field/field_header.cpp
  field/field_layout.cpp
  field/field_tracking.cpp
  field/field.cpp
  field/field_group.cpp
  field/field_manager.cpp
  field/field_sync.cpp
  field/field_utils.cpp
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
  property_checks/property_check.cpp
  property_checks/field_nan_check.cpp
  property_checks/field_within_interval_check.cpp
  property_checks/mass_and_energy_conservation_check.cpp
  algorithm/eamxx_data_interpolation.cpp
  algorithm/eamxx_fv_phys_rrtmgp_active_gases_workaround.cpp
  algorithm/eamxx_time_interpolation.cpp
)

# Append ETI sources (I didn't do it above for clarity of reading)
list (APPEND SHARE_SRC
  field/eti/field_eti_deep_copy_double_device.cpp
  field/eti/field_eti_deep_copy_float_device.cpp
  field/eti/field_eti_deep_copy_int_device.cpp
  field/eti/field_eti_get_strided_view_double_device.cpp
  field/eti/field_eti_get_strided_view_float_device.cpp
  field/eti/field_eti_get_strided_view_int_device.cpp
  field/eti/field_eti_get_view_double_device.cpp
  field/eti/field_eti_get_view_float_device.cpp
  field/eti/field_eti_get_view_int_device.cpp
  field/eti/field_eti_update_double_device.cpp
  field/eti/field_eti_update_float_device.cpp
  field/eti/field_eti_update_impl_double_device.cpp
  field/eti/field_eti_update_impl_fill_double_device.cpp
  field/eti/field_eti_update_impl_fill_float_device.cpp
  field/eti/field_eti_update_impl_fill_int_device.cpp
  field/eti/field_eti_update_impl_float_device.cpp
  field/eti/field_eti_update_impl_int_device.cpp
  field/eti/field_eti_update_int_device.cpp
  field/eti/field_utils_eti_compute_mask_le_int.cpp
  field/eti/field_utils_eti_compute_mask_le_float.cpp
  field/eti/field_utils_eti_compute_mask_le_double.cpp
  field/eti/field_utils_eti_compute_mask_ne_int.cpp
  field/eti/field_utils_eti_compute_mask_ne_float.cpp
  field/eti/field_utils_eti_compute_mask_ne_double.cpp
  field/eti/field_eti_update_max_min_double_device.cpp
  field/eti/field_eti_update_max_min_float_device.cpp
  field/eti/field_eti_update_max_min_impl_double_device.cpp
  field/eti/field_eti_update_max_min_impl_fill_double_device.cpp
  field/eti/field_eti_update_max_min_impl_fill_float_device.cpp
  field/eti/field_eti_update_max_min_impl_fill_int_device.cpp
  field/eti/field_eti_update_max_min_impl_float_device.cpp
  field/eti/field_eti_update_max_min_impl_int_device.cpp
  field/eti/field_eti_update_max_min_int_device.cpp
)

if (EAMXX_ENABLE_GPU)
  # The HostOrDevice::Device is different from HostOrDevice::Host,
  # so we need to ETI also these fcns for HostOrDevice::Host
  # NOTE: we DON'T do ETI for the update/max/min utils, as they are
  #       only used on device (at least for now). OTOH, the methods
  #       to extract views are used in a few places
  list (APPEND SHARE_SRC
    field/eti/field_eti_get_strided_view_double_host.cpp
    field/eti/field_eti_get_strided_view_float_host.cpp
    field/eti/field_eti_get_strided_view_int_host.cpp
    field/eti/field_eti_get_view_double_host.cpp
    field/eti/field_eti_get_view_float_host.cpp
    field/eti/field_eti_get_view_int_host.cpp
  )
endif()

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
  eamxx_utils)

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
