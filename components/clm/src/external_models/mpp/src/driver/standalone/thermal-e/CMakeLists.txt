set(THERMALE_PROBLEM_SOURCES
  mass_and_heat_model_problem.F90
  heat_transport_1D_problem.F90
  th_manoli2014_problem.F90
  )

include_directories(${CMAKE_BINARY_DIR}/src/mpp/dtypes)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/th)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/thermal)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/thermal-e)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/vsfm)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/util)
include_directories(${CMAKE_BINARY_DIR}/src/driver/standalone/util)

link_libraries(${MPP_LIBRARIES})
link_libraries(${PETSC_LIBRARIES})
link_libraries(${NETCDF_LIBRARY})
link_libraries(${NETCDF_F90_LIBRARIES})

include(add_mpp_library)
add_mpp_library(thermale_problems ${THERMALE_PROBLEM_SOURCES})
set_target_properties (thermale_problems PROPERTIES COMPILE_DEFINITIONS "USE_PETSC_LIB")

set(MPP_LIBRARIES thermale_problems; ${MPP_LIBRARIES} PARENT_SCOPE)
set(MPP_LIBRARIES thermale_problems; ${MPP_LIBRARIES})

include(add_mpp_executable)
add_mpp_executable(heat_transport_1D heat_transport_1D.F90)
add_mpp_executable(heat_transport_1D_with_advection heat_transport_1D_with_advection.F90)
add_mpp_executable(mass_and_heat_model mass_and_heat_model.F90)
add_mpp_executable(th_manoli2014 th_manoli2014.F90)

if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS thermale_problems DESTINATION lib)

  install(TARGETS heat_transport_1D DESTINATION bin)
  install(TARGETS heat_transport_1D_with_advection DESTINATION bin)
  install(TARGETS mass_and_heat_model DESTINATION bin)
  install(TARGETS th_manoli2014 DESTINATION bin)

  file(GLOB HEADERS *.mod)
  install(FILES ${HEADERS} DESTINATION include/mpp/driver)
endif()

