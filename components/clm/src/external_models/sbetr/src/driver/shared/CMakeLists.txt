set(BETR_DRIVER_SHARED_SOURCES
  BeTRType.F90
  bncdio_pio.F90
  bhistFileMod.F90
  BeTRStatusSimType.F90
  BeTRSimulation.F90
  )

include_directories(${CMAKE_BINARY_DIR}/src/shr)
include_directories(${CMAKE_BINARY_DIR}/src/stub_clm)
include_directories(${CMAKE_BINARY_DIR}/src/Applications/app_util)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_util)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_math)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_dtype)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_bgc)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_core)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_para)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_main)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_grid)
include_directories(${CMAKE_BINARY_DIR}/src/betr/betr_rxns)

include(add_betr_library)
add_betr_library(betr_driver_shared ${BETR_DRIVER_SHARED_SOURCES})

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w")
if (SBETR)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DSBETR")
endif()
set(BETR_LIBRARIES betr_driver_shared;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES betr_driver_shared;${BETR_LIBRARIES})

#X#add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS betr_driver_shared DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/betr/driver/shared)
endif()
