set(EMI_EM_SOURCES
  ExternalModelBaseType.F90
)

include_directories(${CMAKE_BINARY_DIR}/elm_stub/shr)
include_directories(${CMAKE_BINARY_DIR}/elm_stub/utils)
include_directories(${CMAKE_BINARY_DIR}/constants)
include_directories(${CMAKE_BINARY_DIR}/emi_data_types)
include_directories(${CMAKE_BINARY_DIR}/emi_data_definition)

include(add_emi_library)
add_emi_library(emi_em ${EMI_EM_SOURCES})

set(EMI_LIBRARIES emi_em;${EMI_LIBRARIES} PARENT_SCOPE)
set(EMI_LIBRARIES emi_em;${EMI_LIBRARIES})

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS emi_em DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/)
endif()

add_subdirectory(betr)
add_subdirectory(fates)
add_subdirectory(ptm)
add_subdirectory(vsfm)
