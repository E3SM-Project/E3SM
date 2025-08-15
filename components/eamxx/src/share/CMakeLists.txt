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
  field/field_sync.cpp
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
  util/eamxx_data_interpolation.cpp
  util/eamxx_fv_phys_rrtmgp_active_gases_workaround.cpp
  util/eamxx_time_interpolation.cpp
  util/eamxx_time_stamp.cpp
  util/eamxx_timing.cpp
  util/eamxx_utils.cpp
  util/eamxx_bfbhash.cpp
  util/eamxx_repro_sum_mod.F90
)

# Append ETI sources (I didn't do it above for clarity of reading)
list (APPEND SHARE_SRC
  field/eti/field_eti_deep_copy_double_host.cpp
  field/eti/field_eti_deep_copy_float_host.cpp
  field/eti/field_eti_deep_copy_int_host.cpp
  field/eti/field_eti_get_strided_view_double_host.cpp
  field/eti/field_eti_get_strided_view_float_host.cpp
  field/eti/field_eti_get_strided_view_int_host.cpp
  field/eti/field_eti_get_view_double_host.cpp
  field/eti/field_eti_get_view_float_host.cpp
  field/eti/field_eti_get_view_int_host.cpp
  field/eti/field_eti_update_double_host.cpp
  field/eti/field_eti_update_float_host.cpp
  field/eti/field_eti_update_impl_double_host.cpp
  field/eti/field_eti_update_impl_fill_double_host.cpp
  field/eti/field_eti_update_impl_fill_float_host.cpp
  field/eti/field_eti_update_impl_fill_int_host.cpp
  field/eti/field_eti_update_impl_float_host.cpp
  field/eti/field_eti_update_impl_int_host.cpp
  field/eti/field_eti_update_int_host.cpp
  field/eti/field_utils_eti_compute_mask_le_int.cpp
  field/eti/field_utils_eti_compute_mask_le_float.cpp
  field/eti/field_utils_eti_compute_mask_le_double.cpp
  field/eti/field_utils_eti_compute_mask_ne_int.cpp
  field/eti/field_utils_eti_compute_mask_ne_float.cpp
  field/eti/field_utils_eti_compute_mask_ne_double.cpp
  field/eti/field_eti_update_max_min_double_host.cpp
  field/eti/field_eti_update_max_min_float_host.cpp
  field/eti/field_eti_update_max_min_impl_double_host.cpp
  field/eti/field_eti_update_max_min_impl_fill_double_host.cpp
  field/eti/field_eti_update_max_min_impl_fill_float_host.cpp
  field/eti/field_eti_update_max_min_impl_fill_int_host.cpp
  field/eti/field_eti_update_max_min_impl_float_host.cpp
  field/eti/field_eti_update_max_min_impl_int_host.cpp
  field/eti/field_eti_update_max_min_int_host.cpp
)

if (EAMXX_ENABLE_GPU)
  # The HostOrDevice::Device is different from HostOrDevice::Host,
  # so we need to ETI also these fcns for HostOrDevice::Host
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
endif()

if (EAMXX_ENABLE_EXPERIMENTAL_CODE)
  list (APPEND
    SHARE_SRC
    grid/remap/refining_remapper_rma.cpp
  )
endif()

add_library(scream_share ${SHARE_SRC})
set_target_properties(scream_share PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
target_include_directories (scream_share PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(scream_share PUBLIC
  ${SCREAM_SRC_DIR}
  ${SCREAM_BIN_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# Used in the data interpolation
target_link_libraries(scream_share PUBLIC stdc++fs)

target_link_libraries(scream_share PUBLIC csm_share)

if (GPTL_PATH)
  target_include_directories(scream_share PUBLIC ${GPTL_PATH})
else()
  target_include_directories(scream_share PUBLIC ${INSTALL_SHAREDPATH}/include)
endif()

if (EAMXX_ENABLE_PYTHON)
  if (NOT Python_EXECUTABLE)
    message (FATAL_ERROR "You must set Python_EXECUTABLE to point to a valid python3 interpreter")
  endif()

  execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "import sys; assert sys.version_info >= (3,9), 'Python version must be >= 3.9'"
    RESULT_VARIABLE python_version_check_result
    ERROR_VARIABLE python_version_check_error
  )
  if (python_version_check_result)
      message(FATAL_ERROR "Python version check failed: ${python_version_check_error}")
  endif()

  # Need pybind11
  execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "import pybind11; print(pybind11.get_cmake_dir())"
    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE pybind11_ROOT)

  find_package(pybind11 REQUIRED HINTS ${pybind11_ROOT})

  # Add py sources and link pybind11
  target_sources (scream_share PUBLIC
    eamxx_pysession.cpp)

  target_link_libraries(scream_share PUBLIC pybind11::embed)
  target_compile_definitions (scream_share PUBLIC EAMXX_HAS_PYTHON)
endif()

target_link_libraries(scream_share PUBLIC ekat::AllLibs pioc)

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

# When dealing with subtle bugs in read/write of unit tests baselines,
# using ascii baselines can help inspect that the file indeed looks as expected
option (EAMXX_ASCII_BASELINES "Whether physics units tests should assume ascii (instead of binary) format for baselines" OFF)
if (EAMXX_ASCII_BASELINES)
  target_compile_definitions (scream_share PUBLIC EAMXX_ASCII_BASELINES)
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
  target_include_directories (scream_test_support PUBLIC ${EKAT_SOURCE_DIR}/extern/Catch2/single_include)

  add_subdirectory(tests)
endif()
