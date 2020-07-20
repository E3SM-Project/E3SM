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

# TODO: The current version of kokkos we are using in e3sm has a nvcc_wrapper
#       which breaks the cmake toolchain features detection mechanisms.
#       If we update kokkos to, say, current develop (as of 06/2020),
#       the following line would work fine.
# target_compile_features(ekat PUBLIC cxx_std_11)

# Enable C Pre-Processor for Fortran.
# No reason to pollute C/C++ flags (CPP always runs), so add the flag only for Fortran files.
target_compile_options(ekat PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:-cpp>)

target_link_libraries(ekat PUBLIC ${EKAT_TPL_LIBRARIES})
set_target_properties(ekat PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ekat_f90_modules)
