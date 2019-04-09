set(EMI_EMI_DATA_DEFINITION_SOURCES
  EMI_Atm2LndType_DataMod.F90
  EMI_CanopyStateType_DataMod.F90
  EMI_ChemStateType_DataMod.F90
  EMI_CNCarbonStateType_DataMod.F90
  EMI_EnergyFluxType_DataMod.F90
  EMI_SoilHydrologyType_DataMod.F90
  EMI_SoilStateType_DataMod.F90
  EMI_TemperatureType_DataMod.F90
  EMI_WaterFluxType_DataMod.F90
  EMI_WaterStateType_DataMod.F90
  EMI_DataMod.F90
)

include_directories(${CMAKE_BINARY_DIR}/elm_stub/shr)
include_directories(${CMAKE_BINARY_DIR}/elm_stub/utils)
include_directories(${CMAKE_BINARY_DIR}/constants)
include_directories(${CMAKE_BINARY_DIR}/emi_data_types)

include(add_emi_library)
add_emi_library(emi_emi_data_definition ${EMI_EMI_DATA_DEFINITION_SOURCES})

set(EMI_LIBRARIES emi_emi_data_definition;${EMI_LIBRARIES} PARENT_SCOPE)
set(EMI_LIBRARIES emi_emi_data_definition;${EMI_LIBRARIES})

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS emi_emi_data_definition DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()

