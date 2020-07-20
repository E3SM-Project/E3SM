set(CIMEROOT ${SCREAM_BASE_DIR}/../../cime)
list(APPEND CMAKE_MODULE_PATH ${CIMEROOT}/src/CMake)

set(GENF90 ${CIMEROOT}/src/externals/genf90/genf90.pl)
set(ENABLE_GENF90 True)
include(genf90_utils)
include(Sourcelist_utils)

set(PHYSICS_COMMON_SRCS
  physics_common_iso_f.F90
  physics_common.cpp
  physics_only_grids_manager.cpp
  ${SCREAM_BASE_DIR}/../cam/src/physics/cam/physics_utils.F90
  ${SCREAM_BASE_DIR}/../cam/src/physics/cam/scream_abortutils.F90
)

# Add ETI source files if not on CUDA
if (NOT CUDA_BUILD)
  list(APPEND P3_SRCS
  physics_saturation.cpp)
endif()

add_library(physics_common ${PHYSICS_COMMON_SRCS})
target_include_directories(physics_common PUBLIC ${SCREAM_INCLUDE_DIRS})
target_include_directories(physics_common SYSTEM PUBLIC ${SCREAM_TPL_INCLUDE_DIRS} ${CIMEROOT}/src/share/include)
set_target_properties(physics_common PROPERTIES
  Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})
target_link_libraries(physics_common scream_share ${SCREAM_TPL_LIBRARIES})

#if (NOT SCREAM_LIB_ONLY)
#  add_subdirectory(tests)
#endif()

