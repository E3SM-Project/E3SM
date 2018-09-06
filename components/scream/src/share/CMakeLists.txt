set(SHARE_SRC
  util/scream_util.cpp
  util/array_io.cpp
  util/array_io_mod.f90
)

set(SHARE_HEADERS
  util/scream_kokkos.hpp
  util/scream_util.hpp
)

add_library(scream_share ${SHARE_SRC})
target_include_directories(scream_share PUBLIC ${SCREAM_INCLUDE_DIRS})
set_target_properties(scream_share PROPERTIES
  Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})

add_subdirectory(tests)
