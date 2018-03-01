set(BETR_APP_UTIL_SOURCES
   ApplicationsFactory.F90
  )


include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_core)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/bgcfarm_util)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/CENT_ECACNP)
#include_directories(${CMAKE_BINARY_DIR}/src/Applications/soil-farm/CENT_ECADCNP)

include(add_betr_library)
add_betr_library(app_util ${BETR_APP_UTIL_SOURCES})


set(BETR_LIBRARIES app_util;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES app_util;${BETR_LIBRARIES})

link_libraries(${BETR_LIBRARIES};${BETR_TPLS})

#X#add_subdirectory(tests)


if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS app_util DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/betr/Applications/app_util)
endif()
