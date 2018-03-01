set(BETR_PARA_SOURCES
  TracerParamSetMod.F90
  TracerParamSetWatIsoMod.F90
  TracerParamsMod.F90

  )

include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_bgc)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_grid)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_dtype)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_core)

include(add_betr_library)
add_betr_library(betr_para ${BETR_PARA_SOURCES})

set(BETR_LIBRARIES betr_para;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES betr_para;${BETR_LIBRARIES})

#add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS betr_para DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/betr/para)
endif()
