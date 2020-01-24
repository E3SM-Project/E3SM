
set(SHOC_SRCS
  shoc.cpp
  atmosphere_macrophysics.cpp
  scream_shoc_interface.F90
  ${SCREAM_BASE_DIR}/../cam/src/physics/cam/shoc.F90
)

set(SHOC_HEADERS
  shoc.hpp
  atmosphere_macrophysics.hpp
  scream_shoc_interface.hpp
)

add_library(shoc ${SHOC_SRCS})
target_include_directories(shoc PUBLIC ${SCREAM_INCLUDE_DIRS} ${SCREAM_TPL_INCLUDE_DIRS})
set_target_properties(shoc PROPERTIES
  Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()

