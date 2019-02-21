set(EMI_EMI_SOURCES
  ExternalModelInterfaceMod.F90
)

include_directories(${CMAKE_BINARY_DIR}/elm_stub/shr)
include_directories(${CMAKE_BINARY_DIR}/elm_stub/utils)
include_directories(${CMAKE_BINARY_DIR}/elm_stub/elm_subgrid_types)
include_directories(${CMAKE_BINARY_DIR}/elm_stub/elm_types_stub)
include_directories(${CMAKE_BINARY_DIR}/constants)
include_directories(${CMAKE_BINARY_DIR}/emi_data_types)
include_directories(${CMAKE_BINARY_DIR}/emi_data_definition)
include_directories(${CMAKE_BINARY_DIR}/emi_data_exchange)
include_directories(${CMAKE_BINARY_DIR}/em/betr)
include_directories(${CMAKE_BINARY_DIR}/em/fates)
include_directories(${CMAKE_BINARY_DIR}/em/ptm)
include_directories(${CMAKE_BINARY_DIR}/em/stub)
include_directories(${CMAKE_BINARY_DIR}/em/vsfm)

include(add_emi_library)
add_emi_library(emi_emi ${EMI_EMI_SOURCES})

set(EMI_LIBRARIES emi_emi;${EMI_LIBRARIES} PARENT_SCOPE)
set(EMI_LIBRARIES emi_emi_data_types;${EMI_LIBRARIES})

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS emi_emi DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()
