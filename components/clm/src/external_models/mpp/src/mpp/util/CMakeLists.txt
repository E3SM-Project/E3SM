set(MPP_UTIL_SOURCES
  ArrayDimThree.F90
  EOSWaterMod.F90
  MultiPhysicsProbConstants.F90
  mpp_abortutils.F90
  mpp_bounds.F90
  mpp_shr_log_mod.F90
  mpp_varcon.F90
  mpp_varctl.F90
  mpp_varpar.F90
  )

#include_directories(${CMAKE_BINARY_DIR}/src/stub_clm)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/dtypes)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/thermal)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/vsfm)

include(add_mpp_library)
add_mpp_library(mpp_util ${MPP_UTIL_SOURCES})
set_target_properties (mpp_util PROPERTIES COMPILE_DEFINITIONS "USE_PETSC_LIB")

set(MPP_LIBRARIES mpp_util;${MPP_LIBRARIES} PARENT_SCOPE)
set(MPP_LIBRARIES mpp_util;${MPP_LIBRARIES})


if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS mpp_util DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/mpp/util)
endif()

