set(VSFM_DRIVER_SOURCES
  standalone_mpp.F90
  )

add_subdirectory(util)
add_subdirectory(vsfm)
add_subdirectory(thermal-e)

include_directories(${CMAKE_BINARY_DIR}/src/mpp/dtypes)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/th)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/thermal)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/thermal-e)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/vsfm)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/util)
include_directories(${CMAKE_BINARY_DIR}/src/driver/standalone/thermal-e)
include_directories(${CMAKE_BINARY_DIR}/src/driver/standalone/vsfm)

link_libraries(${MPP_LIBRARIES})
link_libraries(${PETSC_LIBRARIES})
link_libraries(${NETCDF_LIBRARY})
link_libraries(${NETCDF_F90_LIBRARIES})

include(add_mpp_executable)
add_mpp_executable(standalone_mpp standalone_mpp.F90)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS standalone_mpp DESTINATION bin)
  file(GLOB HEADERS *.mod)
  install(FILES ${HEADERS} DESTINATION include/mpp/driver)
endif()

