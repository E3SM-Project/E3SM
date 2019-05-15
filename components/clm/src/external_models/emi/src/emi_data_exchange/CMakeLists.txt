set(EMI_EMI_DATA_EXCHANGE_SOURCES
  EMI_Atm2LndType_ExchangeMod.F90
  EMI_CanopyStateType_ExchangeMod.F90
  EMI_ChemStateType_ExchangeMod.F90
  EMI_CNCarbonStateType_ExchangeMod.F90
  EMI_ColumnType_ExchangeMod.F90
  EMI_EnergyFluxType_ExchangeMod.F90
  EMI_Filter_ExchangeMod.F90
  EMI_Landunit_ExchangeMod.F90
  EMI_SoilHydrologyType_ExchangeMod.F90
  EMI_SoilStateType_ExchangeMod.F90
  EMI_TemperatureType_ExchangeMod.F90
  EMI_WaterFluxType_ExchangeMod.F90
  EMI_WaterStateType_ExchangeMod.F90
)

include_directories(${CMAKE_BINARY_DIR}/elm_stub/shr)
include_directories(${CMAKE_BINARY_DIR}/elm_stub/utils)
include_directories(${CMAKE_BINARY_DIR}/elm_stub/elm_types_stub)
include_directories(${CMAKE_BINARY_DIR}/elm_stub/elm_subgrid_types)
include_directories(${CMAKE_BINARY_DIR}/constants)
include_directories(${CMAKE_BINARY_DIR}/elm_subgrid_types)
include_directories(${CMAKE_BINARY_DIR}/emi_data_definition)
include_directories(${CMAKE_BINARY_DIR}/emi_data_types)

include(add_emi_library)
add_emi_library(emi_emi_data_exchange ${EMI_EMI_DATA_EXCHANGE_SOURCES})

set(EMI_LIBRARIES emi_emi_data_exchange;${EMI_LIBRARIES} PARENT_SCOPE)
set(EMI_LIBRARIES emi_emi_data_exchange;${EMI_LIBRARIES})

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS emi_emi_data_exchange DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()

