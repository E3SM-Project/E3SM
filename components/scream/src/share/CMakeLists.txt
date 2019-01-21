set(SHARE_SRC
  scream_assert.cpp
  field/field_alloc_prop.cpp
  field/field_identifier.cpp
  field/field_header.cpp
  field/field_tracking.cpp
  scream_session.cpp
  mpi/scream_comm.cpp
  util/scream_utils.cpp
  util/array_io.cpp
  util/array_io_mod.f90
)

set(SHARE_HEADERS
  atmosphere_process.hpp
  scream_assert.hpp
  field/field.hpp
  field/field_alloc_prop.hpp
  field/field_identifier.hpp
  field/field_header.hpp
  field/field_repository.hpp
  field/field_tag.hpp
  field/field_tracking.cpp
  parameter_list.hpp
  scream_session.hpp
  mpi/scream_comm.hpp
  util/factory.hpp
  util/file_utils.hpp
  util/math_utils.hpp
  util/scream_utils.hpp
)

add_library(scream_share ${SHARE_SRC})
target_include_directories(scream_share PUBLIC ${SCREAM_INCLUDE_DIRS} ${SCREAM_TPL_INCLUDE_DIRS})
set_target_properties(scream_share PROPERTIES
  Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})
# link_directories(${SCREAM_TPL_LIBRARY_DIRS})
# target_link_libraries(scream_share PUBLIC ${SCREAM_TPL_LIBRARIES})

add_subdirectory(tests)
