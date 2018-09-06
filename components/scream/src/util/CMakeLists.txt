set(UTIL_SRC
  scream_util.cpp
  array_io.cpp
  array_io_mod.f90
)

set(UTIL_HEADERS
  scream_kokkos.hpp
  scream_util.hpp
)

add_library(scream_util ${UTIL_SRC})
target_include_directories(scream_util PUBLIC ${SCREAM_INCLUDE_DIRS})
set_target_properties(scream_util PROPERTIES
  Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})

add_subdirectory(tests)
