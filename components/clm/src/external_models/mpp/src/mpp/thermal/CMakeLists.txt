set(MPP_THERMAL_SOURCES
  GoveqnThermalKSPTemperatureSSWType.F90
  GoveqnThermalKSPTemperatureSnowType.F90
  GoveqnThermalKSPTemperatureSoilType.F90
  MultiPhysicsProbThermal.F90
  SystemOfEquationsThermalAuxMod.F90
  SystemOfEquationsThermalAuxType.F90
  SystemOfEquationsThermalType.F90
  ThermalKSPTemperatureBaseAuxType.F90
  ThermalKSPTemperatureSSWAuxMod.F90
  ThermalKSPTemperatureSSWAuxType.F90
  ThermalKSPTemperatureSnowAuxMod.F90
  ThermalKSPTemperatureSnowAuxType.F90
  ThermalKSPTemperatureSoilAuxMod.F90
  ThermalKSPTemperatureSoilAuxType.F90
  )

#include_directories(${CMAKE_BINARY_DIR}/src/stub_clm)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/dtypes)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/util)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/vsfm)

include(add_mpp_library)
add_mpp_library(mpp_thermal ${MPP_THERMAL_SOURCES})
set_target_properties (mpp_thermal PROPERTIES COMPILE_DEFINITIONS "USE_PETSC_LIB")

set(MPP_LIBRARIES mpp_thermal;${MPP_LIBRARIES} PARENT_SCOPE)
set(MPP_LIBRARIES mpp_thermal;${MPP_LIBRARIES})


if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS mpp_thermal DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/mpp/thermal)
endif()

