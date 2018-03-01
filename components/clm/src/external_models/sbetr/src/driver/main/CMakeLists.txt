set(BETR_DRIVER_SOURCES
  BeTRSimulationFactory.F90
  sbetrDriverMod.F90
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
include_directories(${CMAKE_BINARY_DIR}/src/driver/shared)
include_directories(${CMAKE_BINARY_DIR}/src/driver/alm)
include_directories(${CMAKE_BINARY_DIR}/src/driver/clm)
include_directories(${CMAKE_BINARY_DIR}/src/driver/standalone)

include(add_betr_library)
add_betr_library(betr_driver_main ${BETR_DRIVER_SOURCES})

#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w")

set(BETR_LIBRARIES betr_driver_main;${BETR_LIBRARIES} PARENT_SCOPE)
set(BETR_LIBRARIES betr_driver_main;${BETR_LIBRARIES})

link_libraries(${BETR_LIBRARIES};${BETR_TPLS})
include(add_betr_executable)
add_betr_executable(sbetr sbetr.F90)

#X#add_subdirectory(tests)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS sbetr DESTINATION bin)
  file(GLOB HEADERS *.mod)
  install(FILES ${HEADERS} DESTINATION include/betr/driver/main)
endif()
