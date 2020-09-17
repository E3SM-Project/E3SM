set(CIMEROOT ${SCREAM_BASE_DIR}/../../cime)
list(APPEND CMAKE_MODULE_PATH ${CIMEROOT}/src/CMake)

set(GENF90 ${CIMEROOT}/src/externals/genf90/genf90.pl)
set(ENABLE_GENF90 True)
include(genf90_utils)
include(Sourcelist_utils)

set(SHOC_SRCS
  shoc_f90.cpp
  shoc_functions_f90.cpp
  shoc_ic_cases.cpp
  shoc_iso_c.f90
  shoc_iso_f.f90
  ${SCREAM_BASE_DIR}/../cam/src/physics/cam/shoc.F90
  atmosphere_macrophysics.cpp
  scream_shoc_interface.F90
  shoc_diag_second_moments_srf.cpp
  shoc_diag_second_moments_ubycond.cpp
)

set(SHOC_HEADERS
  shoc.hpp
  atmosphere_macrophysics.hpp
  scream_shoc_interface.hpp
)

# Add ETI source files if not on CUDA
if (NOT CUDA_BUILD)
  list(APPEND SHOC_SRCS
    shoc_update_host_dse.cpp
    shoc_calc_shoc_varorcovar.cpp
    shoc_calc_shoc_vertflux.cpp
  )
endif()

add_library(shoc ${SHOC_SRCS})
target_include_directories(shoc PUBLIC ${SCREAM_INCLUDE_DIRS})
target_include_directories(shoc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../share)
target_include_directories(shoc SYSTEM PUBLIC ${SCREAM_TPL_INCLUDE_DIRS} ${CIMEROOT}/src/share/include)
set_target_properties(shoc PROPERTIES
  Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})
target_link_libraries(shoc physics_share scream_share ${SCREAM_TPL_LIBRARIES})

if (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()

