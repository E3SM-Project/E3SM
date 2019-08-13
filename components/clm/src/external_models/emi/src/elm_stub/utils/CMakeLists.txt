set(EMI_UTILS_SOURCES
  abortutils.F90
  clm_varcon.F90
  clm_varctl.F90
  clm_varpar.F90
  column_varcon.F90
  decompMod.F90
  landunit_varcon.F90
  spmdMod.F90
)

include_directories(${CMAKE_BINARY_DIR}/elm_stub/shr)

include(add_emi_library)
add_emi_library(emi_utils ${EMI_UTILS_SOURCES})

set(EMI_LIBRARIES emi_utils;${EMI_LIBRARIES} PARENT_SCOPE)
set(EMI_LIBRARIES emi_utils;${EMI_LIBRARIES})


if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS emi_utils DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include)
endif()

