
set(SHOC_SRCS
  shoc.cpp
  shoc_f90.cpp
  shoc_functions_f90.cpp
  shoc_iso_c.f90
  shoc_iso_f.f90
  ${SCREAM_BASE_DIR}/../cam/src/physics/cam/shoc.F90
  atmosphere_macrophysics.cpp
  scream_shoc_interface.F90
)

set(SHOC_HEADERS
  shoc.hpp
  atmosphere_macrophysics.hpp
  scream_shoc_interface.hpp
)

# Add ETI source files if not on CUDA
if (NOT CUDA_BUILD)
  list(APPEND SHOC_SRCS
      )
endif()

add_library(shoc ${SHOC_SRCS})
target_include_directories(shoc PUBLIC ${SCREAM_INCLUDE_DIRS})
target_include_directories(shoc SYSTEM PUBLIC ${SCREAM_TPL_INCLUDE_DIRS})
set_target_properties(shoc PROPERTIES
  Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()

