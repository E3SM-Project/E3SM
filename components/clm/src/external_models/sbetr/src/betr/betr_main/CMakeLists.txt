set(BETR_MAIN_SOURCES
  BetrBGCMod.F90
  TracerBalanceMod.F90
  )

include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_grid)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_dtype)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_core)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_para)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_rxns)


include(add_betr_library)
add_betr_library(betr_main ${BETR_MAIN_SOURCES})

set(BETR_LIBRARIES betr_main;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES betr_main;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS betr_main DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/betr/main)
endif()
