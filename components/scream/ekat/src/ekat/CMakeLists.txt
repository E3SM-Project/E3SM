set(EKAT_SOURCES
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
)

add_library(ekat ${EKAT_SOURCES})

target_include_directories(ekat PUBLIC ${EKAT_INCLUDE_DIRS})
target_include_directories(ekat SYSTEM PUBLIC ${EKAT_TPL_INCLUDE_DIRS})

target_compile_features(ekat PUBLIC cxx_std_11)

target_link_libraries(ekat PUBLIC ${EKAT_TPL_LIBRARIES})
set_target_properties(ekat PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ekat_f90_modules)
