include (EkatSetCompilerFlags)

set(SHARE_SRC
  eamxx_config.cpp
  eamxx_session.cpp
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
  property_checks/mass_and_energy_column_conservation_check.cpp
  util/eamxx_data_interpolation.cpp
  util/eamxx_fv_phys_rrtmgp_active_gases_workaround.cpp
  util/eamxx_time_interpolation.cpp
  util/eamxx_time_stamp.cpp
  util/eamxx_timing.cpp
  util/eamxx_utils.cpp
  util/eamxx_bfbhash.cpp
)

# Append ETI sources (I didn't do it above for clarity of reading)
list (APPEND SHARE_SRC
  field/eti/field_eti_deep_copy_double_device.cpp
  field/eti/field_eti_deep_copy_double_host.cpp
  field/eti/field_eti_deep_copy_float_device.cpp
  field/eti/field_eti_deep_copy_float_host.cpp
  field/eti/field_eti_deep_copy_int_device.cpp
  field/eti/field_eti_deep_copy_int_host.cpp
  field/eti/field_eti_get_strided_view_double_device.cpp
  field/eti/field_eti_get_strided_view_double_host.cpp
  field/eti/field_eti_get_strided_view_float_device.cpp
  field/eti/field_eti_get_strided_view_float_host.cpp
  field/eti/field_eti_get_strided_view_int_device.cpp
  field/eti/field_eti_get_strided_view_int_host.cpp
  field/eti/field_eti_get_view_double_device.cpp
  field/eti/field_eti_get_view_double_host.cpp
  field/eti/field_eti_get_view_float_device.cpp
  field/eti/field_eti_get_view_float_host.cpp
  field/eti/field_eti_get_view_int_device.cpp
  field/eti/field_eti_get_view_int_host.cpp
  field/eti/field_eti_update_double_device.cpp
  field/eti/field_eti_update_double_host.cpp
  field/eti/field_eti_update_float_device.cpp
  field/eti/field_eti_update_float_host.cpp
  field/eti/field_eti_update_impl_double_device.cpp
  field/eti/field_eti_update_impl_double_host.cpp
  field/eti/field_eti_update_impl_fill_double_device.cpp
  field/eti/field_eti_update_impl_fill_double_host.cpp
  field/eti/field_eti_update_impl_fill_float_device.cpp
  field/eti/field_eti_update_impl_fill_float_host.cpp
  field/eti/field_eti_update_impl_fill_int_device.cpp
  field/eti/field_eti_update_impl_fill_int_host.cpp
  field/eti/field_eti_update_impl_float_device.cpp
  field/eti/field_eti_update_impl_float_host.cpp
  field/eti/field_eti_update_impl_int_device.cpp
  field/eti/field_eti_update_impl_int_host.cpp
  field/eti/field_eti_update_int_device.cpp
  field/eti/field_eti_update_int_host.cpp
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
# Used in the data interpolation
target_link_libraries(scream_share PUBLIC stdc++fs)


if (GPTL_PATH)
  target_include_directories(scream_share PUBLIC ${GPTL_PATH})
else()
  target_include_directories(scream_share PUBLIC ${INSTALL_SHAREDPATH}/include)
endif()

target_link_libraries(scream_share PUBLIC ekat pioc)
target_compile_options(scream_share PUBLIC
  $<$<COMPILE_LANGUAGE:Fortran>:${SCREAM_Fortran_FLAGS}>
)

# This used to be in `eamxx_config.h/eamxx_config.f`, but we did accidentally remove
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
    util/eamxx_test_session.cpp
  )
  target_link_libraries(scream_test_support PUBLIC scream_share scream_io)

  add_subdirectory(tests)
endif()
