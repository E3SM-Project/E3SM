set(MPP_TH_SOURCES
  MultiPhysicsProbTH.F90
  SystemOfEquationsTHType.F90
)

include_directories(${CMAKE_BINARY_DIR}/src/mpp/dtypes)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/util)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/vsfm)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/thermal-e)

include(add_mpp_library)
add_mpp_library(mpp_th ${MPP_TH_SOURCES})
set_target_properties (mpp_th PROPERTIES COMPILE_DEFINITIONS "USE_PETSC_LIB")

set(MPP_LIBRARIES mpp_th; ${MPP_LIBRARIES} PARENT_SCOPE)
set(MPP_LIBRARIES mpp_th; ${MPP_LIBRARIES})


if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS mpp_th DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/mpp/th)
endif()

