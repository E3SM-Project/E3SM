set(BETR_GRID_SOURCES
  BeTR_ColumnType.F90
  BeTR_LandunitType.F90
  BeTR_PatchType.F90
  BeTR_landvarconType.F90
  BeTR_pftvarconType.F90
  )

include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)


include(add_betr_library)
add_betr_library(betr_grid ${BETR_GRID_SOURCES})

set(BETR_LIBRARIES betr_grid;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES betr_grid;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS betr_grid DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/betr/grid)
endif()
