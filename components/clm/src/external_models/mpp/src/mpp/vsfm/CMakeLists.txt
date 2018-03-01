set(MPP_VSFM_SOURCES
  GoveqnRichardsODEPressureType.F90
  MultiPhysicsProbVSFM.F90
  PorosityFunctionMod.F90
  RichardsMod.F90
  RichardsODEPressureAuxType.F90
  RichardsODEPressureAuxMod.F90
  RichardsODEPressureConnAuxType.F90
  SaturationFunction.F90
  SystemOfEquationsVSFMAuxType.F90
  SystemOfEquationsVSFMType.F90  
  )

#include_directories(${CMAKE_BINARY_DIR}/src/stub_clm)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/dtypes)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/thermal)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/util)

include(add_mpp_library)
add_mpp_library(mpp_vsfm ${MPP_VSFM_SOURCES})
set_target_properties (mpp_vsfm PROPERTIES COMPILE_DEFINITIONS "USE_PETSC_LIB")

set(MPP_LIBRARIES mpp_vsfm;${MPP_LIBRARIES} PARENT_SCOPE)
set(MPP_LIBRARIES mpp_vsfm;${MPP_LIBRARIES})


if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS mpp_vsfm DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/mpp/vsfm)
endif()

