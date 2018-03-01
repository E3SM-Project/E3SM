set(MPP_THERMAL_E_SOURCES
  MultiPhysicsProbThermalEnthalpy.F90
  GoveqnThermalEnthalpySoilType.F90
  SystemOfEquationsThermalEnthalpyAuxType.F90
  SystemOfEquationsThermalEnthalpyAuxMod.F90
  SystemOfEquationsThermalEnthalpyType.F90
  ThermalEnthalpySoilAuxType.F90
  ThermalEnthalpySoilAuxMod.F90
  ThermalEnthalpyMod.F90
)

#include_directories(${CMAKE_BINARY_DIR}/src/stub_clm)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/dtypes)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/util)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/vsfm)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/thermal)

include(add_mpp_library)
add_mpp_library(mpp_thermal_e ${MPP_THERMAL_E_SOURCES})
set_target_properties (mpp_thermal_e PROPERTIES COMPILE_DEFINITIONS "USE_PETSC_LIB")

set(MPP_LIBRARIES mpp_thermal_e;${MPP_LIBRARIES} PARENT_SCOPE)
set(MPP_LIBRARIES mpp_thermal_e;${MPP_LIBRARIES})


if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS mpp_thermal_e DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/mpp/thermal)
endif()

