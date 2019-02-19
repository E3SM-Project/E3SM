set(EMI_ELM_SUBGRID_TYPES_SOURCES
  ColumnType.F90
  GridcellType.F90
  LandunitType.F90
  TopounitType.F90
  VegetationType.F90
)

include_directories(${CMAKE_BINARY_DIR}/elm_stub/shr)
include_directories(${CMAKE_BINARY_DIR}/elm_stub/utils)

include(add_emi_library)
add_emi_library(emi_elm_subgrid_types ${EMI_ELM_SUBGRID_TYPES_SOURCES})

set(EMI_LIBRARIES emi_elm_subgrid_types;${EMI_LIBRARIES} PARENT_SCOPE)
set(EMI_LIBRARIES emi_elm_subgrid_types;${EMI_LIBRARIES})


if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS emi_elm_subgrid_types DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()

