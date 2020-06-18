set(EMI_EMI_DATA_TYPES_SOURCES
  EMI_DataDimensionMod.F90
)

include_directories(${CMAKE_BINARY_DIR}/elm_stub/shr)
include_directories(${CMAKE_BINARY_DIR}/elm_stub/utils)
include_directories(${CMAKE_BINARY_DIR}/constants)

include(add_emi_library)
add_emi_library(emi_emi_data_types ${EMI_EMI_DATA_TYPES_SOURCES})

set(EMI_LIBRARIES emi_emi_data_types;${EMI_LIBRARIES} PARENT_SCOPE)
set(EMI_LIBRARIES emi_emi_data_types;${EMI_LIBRARIES})

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS emi_emi_data_types DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()

