set(SHARE_SRC
  atmosphere_process_group.cpp
  scream_assert.cpp
  scream_session.cpp
  field/field_alloc_prop.cpp
  field/field_identifier.cpp
  field/field_header.cpp
  field/field_layout.cpp
  field/field_tracking.cpp
  mpi/scream_comm.cpp
  scream_parameter_list.cpp
  scream_parse_yaml_file.cpp
  util/scream_utils.cpp
  util/scream_arch.cpp
  util/array_io.cpp
  util/array_io_mod.f90
  util/time_stamp.cpp
  util/scream_tridiag_bfb.cpp
  util/scream_tridiag_bfb_mod.F90
)

# So far this is pointless, since we are not installing the library
set(SHARE_HEADERS
  atmosphere_process.hpp
  atmosphere_process_group.hpp
  scream_assert.hpp
  field/field.hpp
  field/field_alloc_prop.hpp
  field/field_identifier.hpp
  field/field_header.hpp
  field/field_layout.hpp
  field/field_repository.hpp
  field/field_tag.hpp
  field/field_tracking.cpp
  field/field_utils.hpp
  grid/abstract_grid.hpp
  grid/default_grid.hpp
  grid/grid_utils.hpp
  grid/grids_manager.hpp
  grid/user_provided_grids_manager.hpp
  remap/abstract_remapper.hpp
  remap/identity_remapper.hpp
  remap/inverse_remapper.hpp
  scream_parameter_list.hpp
  scream_parse_yaml_file.hpp
  scream_session.hpp
  mpi/scream_comm.hpp
  util/scream_factory.hpp
  util/file_utils.hpp
  util/math_utils.hpp
  util/scream_utils.hpp
  util/scream_arch.hpp
  util/scream_kokkos.hpp
  util/scream_lin_interp.hpp
  util/scream_rational_constant.hpp
  util/scream_scaling_factor.hpp
  util/scream_units.hpp
  util/time_stamp.hpp
  scream_workspace.hpp
)

add_library(scream_share ${SHARE_SRC})
target_include_directories(scream_share PUBLIC ${SCREAM_INCLUDE_DIRS})
target_include_directories(scream_share SYSTEM PUBLIC ${SCREAM_TPL_INCLUDE_DIRS})
set_target_properties(scream_share PROPERTIES
  Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})
# link_directories(${SCREAM_TPL_LIBRARY_DIRS})
# target_link_libraries(scream_share PUBLIC ${SCREAM_TPL_LIBRARIES})

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
