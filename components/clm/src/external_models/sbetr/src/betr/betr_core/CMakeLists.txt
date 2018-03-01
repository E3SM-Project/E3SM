set(BETR_CORE_SOURCES
  BeTRHistVarType.F90
  BeTRTracerType.F90
  TracerBaseType.F90
  TracerBoundaryCondType.F90
  TracerCoeffType.F90
  TracerFluxType.F90
  TracerStateType.F90
  TransportMod.F90
  BGCReactionsMod.F90
  PlantSoilBGCMod.F90
  )

include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_bgc)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_grid)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_dtype)

include(add_betr_library)
add_betr_library(betr_core ${BETR_CORE_SOURCES})

set(BETR_LIBRARIES betr_core;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES betr_core;${BETR_LIBRARIES})

add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS betr_core DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/betr/core)
endif()
