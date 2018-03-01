set(BETR_RXNS_SOURCES
  ReactionsFactory.F90
  MockBGCReactionsType.F90
  MockPlantSoilBGCType.F90
  DIOCPlantSoilBGCType.F90
  DIOCBGCReactionsType.F90
  H2OIsotopeBGCReactionsType.F90
  H2OIsotopePlantSoilBGCType.F90
  )

include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_bgc)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_dtype)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_core)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_para)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_grid)


include(add_betr_library)
add_betr_library(betr_rxns ${BETR_RXNS_SOURCES})


set(BETR_LIBRARIES betr_rxns;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES betr_rxns;${BETR_LIBRARIES})

link_libraries(${BETR_LIBRARIES};${BETR_TPLS})

#X#add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS betr_rxns DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/betr/rxns)
endif()
