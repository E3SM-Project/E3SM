set(BGCFARM_UTIL_SOURCES
 BiogeoConType.F90
)

if (BETR_BGC EQUAL 1)
  include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
  include(add_betr_library)
  add_betr_library(bgcfarm_util ${BGCFARM_UTIL_SOURCES})

  set(BETR_LIBRARIES bgcfarm_util;${BETR_LIBRARIES} PARENT_SCOPE)
  set(BETR_LIBRARIES bgcfarm_util;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

  if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
    install(TARGETS bgcfarm_util DESTINATION lib)
    file(GLOB HEADERS *.h)
    install(FILES ${HEADERS} DESTINATION include/soil-farm/bgcfarm_util)
  endif()
else()
  include_directories(${CMAKE_BINARY_DIR}/src/betr_share/betr_util)

  include(add_bgcfarm_library)
  add_bgcfarm_library(bgcfarm_util ${BGCFARM_UTIL_SOURCES})

  set(BGCFARM_LIBRARIES bgcfarm_util;${BGCFARM_LIBRARIES} PARENT_SCOPE)
  set(BGCFARM_LIBRARIES bgcfarm_util;${BGCFARM_LIBRARIES})

#add_subdirectory(tests)

  if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
    install(TARGETS bgcfarm_util DESTINATION lib)
    file(GLOB HEADERS *.h)
    install(FILES ${HEADERS} DESTINATION include/bgcfarm/bgcfarm_util)
  endif()
endif()
