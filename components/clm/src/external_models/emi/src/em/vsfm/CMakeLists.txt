set(EMI_EM_VSFM_SOURCES
  ExternalModelVSFMMod.F90
)

include_directories(${CMAKE_BINARY_DIR}/elm_stub/shr)
include_directories(${CMAKE_BINARY_DIR}/elm_stub/utils)
include_directories(${CMAKE_BINARY_DIR}/constants)
include_directories(${CMAKE_BINARY_DIR}/emi_data_types)
include_directories(${CMAKE_BINARY_DIR}/emi_data_definition)
include_directories(${CMAKE_BINARY_DIR}/em/base)

include(add_emi_library)
add_emi_library(emi_em_vsfm ${EMI_EM_VSFM_SOURCES})

set(EMI_LIBRARIES emi_em_vsfm;${EMI_LIBRARIES} PARENT_SCOPE)
set(EMI_LIBRARIES emi_em_vsfm;${EMI_LIBRARIES})


if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS emi_em_vsfm DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
