set(EKAT_SRC
  mpi/scream_comm.cpp
  scream_assert.cpp
  scream_session.cpp
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

add_library(ekat ${EKAT_SRC})
target_include_directories(ekat PUBLIC ${SCREAM_INCLUDE_DIRS})
target_include_directories(ekat SYSTEM PUBLIC ${SCREAM_TPL_INCLUDE_DIRS})
target_link_libraries(ekat PUBLIC ${SCREAM_TPL_LIBRARIES})
set_target_properties(ekat PROPERTIES
  Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
