set(EMI_ELM_TYPES_STUB_SOURCES
  CanopyStateType.F90
  CNCarbonStateType.F90
  ChemStateType.F90
  EnergyFluxType.F90
  SoilHydrologyType.F90
  SoilStateType.F90
  TemperatureType.F90
  WaterStateType.F90
  WaterfluxType.F90
  atm2lndType.F90
  elm_instMod.F90
)

include_directories(${CMAKE_BINARY_DIR}/elm_stub/shr)
include_directories(${CMAKE_BINARY_DIR}/elm_stub/utils)
include_directories(${CMAKE_SOURCE_DIR}/elm_stub/utils)
include_directories(${CMAKE_BINARY_DIR}/elm_stub/elm_subgrid_types)

include(add_emi_library)
add_emi_library(emi_elm_types_stub ${EMI_ELM_TYPES_STUB_SOURCES})

set(EMI_LIBRARIES emi_elm_types_stub;${EMI_LIBRARIES} PARENT_SCOPE)
set(EMI_LIBRARIES emi_elm_types_stub;${EMI_LIBRARIES})


if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS emi_elm_types_stub DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()

