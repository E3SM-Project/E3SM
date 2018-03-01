set(STANDALONE_UTIL_SOURCES
  regression.F90
  )

include(add_mpp_library)
add_mpp_library(standalone_util ${STANDALONE_UTIL_SOURCES})
set_target_properties(standalone_util PROPERTIES COMPILE_DEFINITIONS "USE_PETSC_LIB")

set(MPP_LIBRARIES standalone_util; ${MPP_LIBRARIES} PARENT_SCOPE)
set(MPP_LIBRARIES standalone_util; ${MPP_LIBRARIES})

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS standalone_util DESTINATION lib)
  file(GLOB HEADERS *.mod)
  install(FILES ${HEADERS} DESTINATION include/driver/standalone/util)
endif()
