set(MPP_DRIVER_SOURCES
  mpp.F90
  )

include_directories(${CMAKE_BINARY_DIR}/src/mpp/dtypes)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/thermal)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/vsfm)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/util)

#include(add_mpp_library)
#add_mpp_library(betr_driver ${MPP_DRIVER_SOURCES})

#set(MPP_LIBRARIES mpp_driver;${MPP_LIBRARIES} PARENT_SCOPE)
#set(MPP_LIBRARIES mpp_driver;${MPP_LIBRARIES})

link_libraries(${MPP_LIBRARIES})
link_libraries(${PETSC_LIBRARIES})
link_libraries(${NETCDF_LIBRARY})
link_libraries(${NETCDF_F90_LIBRARIES})
include(add_mpp_executable)
add_mpp_executable(mpp mpp.F90)

message ("MPP_LIBRARIES = ${MPP_LIBRARIES} ")

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS mpp DESTINATION bin)
  file(GLOB HEADERS *.mod)
  install(FILES ${HEADERS} DESTINATION include/mpp/driver)
endif()

