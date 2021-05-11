set(EMI_CONSTANTS_SOURCES
  EMI_Atm2LndType_Constants.F90
  EMI_CanopyStateType_Constants.F90
  EMI_ChemStateType_Constants.F90
  EMI_CNCarbonStateType_Constants.F90
  EMI_ColumnType_Constants.F90
  EMI_EnergyFluxType_Constants.F90
  EMI_Filter_Constants.F90
  EMI_Landunit_Constants.F90
  EMI_SoilHydrologyType_Constants.F90
  EMI_SoilStateType_Constants.F90
  EMI_TemperatureType_Constants.F90
  EMI_WaterFluxType_Constants.F90
  EMI_WaterStateType_Constants.F90
  ExternalModelConstants.F90 )

include(add_emi_library)
add_emi_library(emi_constants ${EMI_CONSTANTS_SOURCES})

set(EMI_LIBRARIES emi_constants;${EMI_LIBRARIES} PARENT_SCOPE)
set(EMI_LIBRARIES emi_constants;${EMI_LIBRARIES})

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS emi_constants DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()

