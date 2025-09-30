add_library(eamxx_field
  field_alloc_prop.cpp
  field_identifier.cpp
  field_header.cpp
  field_layout.cpp
  field_tracking.cpp
  field.cpp
  field_group.cpp
  field_sync.cpp
  field_utils.cpp
  eti/field_eti_deep_copy_double_device.cpp
  eti/field_eti_deep_copy_float_device.cpp
  eti/field_eti_deep_copy_int_device.cpp
  eti/field_eti_get_strided_view_double_device.cpp
  eti/field_eti_get_strided_view_float_device.cpp
  eti/field_eti_get_strided_view_int_device.cpp
  eti/field_eti_get_view_double_device.cpp
  eti/field_eti_get_view_float_device.cpp
  eti/field_eti_get_view_int_device.cpp
  eti/field_eti_update_double_device.cpp
  eti/field_eti_update_float_device.cpp
  eti/field_eti_update_impl_double_device.cpp
  eti/field_eti_update_impl_fill_double_device.cpp
  eti/field_eti_update_impl_fill_float_device.cpp
  eti/field_eti_update_impl_fill_int_device.cpp
  eti/field_eti_update_impl_float_device.cpp
  eti/field_eti_update_impl_int_device.cpp
  eti/field_eti_update_int_device.cpp
  eti/field_eti_update_max_min_double_device.cpp
  eti/field_eti_update_max_min_float_device.cpp
  eti/field_eti_update_max_min_impl_double_device.cpp
  eti/field_eti_update_max_min_impl_fill_double_device.cpp
  eti/field_eti_update_max_min_impl_fill_float_device.cpp
  eti/field_eti_update_max_min_impl_fill_int_device.cpp
  eti/field_eti_update_max_min_impl_float_device.cpp
  eti/field_eti_update_max_min_impl_int_device.cpp
  eti/field_eti_update_max_min_int_device.cpp
  eti/field_utils_eti_compute_mask_le_int.cpp
  eti/field_utils_eti_compute_mask_le_float.cpp
  eti/field_utils_eti_compute_mask_le_double.cpp
  eti/field_utils_eti_compute_mask_ne_int.cpp
  eti/field_utils_eti_compute_mask_ne_float.cpp
  eti/field_utils_eti_compute_mask_ne_double.cpp
)

if (EAMXX_ENABLE_GPU)
  # The HostOrDevice::Device is different from HostOrDevice::Host,
  # so we need to ETI also these fcns for HostOrDevice::Host
  # NOTE: we DON'T do ETI for the update/max/min utils, as they are
  #       only used on device (at least for now). OTOH, the methods
  #       to extract views are used in a few places
  target_sources (eamxx_field PRIVATE
    eti/field_eti_get_strided_view_double_host.cpp
    eti/field_eti_get_strided_view_float_host.cpp
    eti/field_eti_get_strided_view_int_host.cpp
    eti/field_eti_get_view_double_host.cpp
    eti/field_eti_get_view_float_host.cpp
    eti/field_eti_get_view_int_host.cpp
  )
endif()

target_link_libraries (eamxx_field PUBLIC
  eamxx_core
  eamxx_utils)

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
