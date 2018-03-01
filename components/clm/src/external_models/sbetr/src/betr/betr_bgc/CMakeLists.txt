set(BETR_BGC_SOURCES
  PlantNutKineticsMod.F90
  )

include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_grid)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_dtype)

include(add_betr_library)
add_betr_library(betr_bgc ${BETR_BGC_SOURCES})

set(BETR_LIBRARIES betr_bgc;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES betr_bgc;${BETR_LIBRARIES})


if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS betr_bgc DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/betr/bgc)
endif()
