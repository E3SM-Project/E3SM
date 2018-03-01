set(MPP_DTYPES_SOURCES
  ConditionType.F90
  ConnectionSetType.F90
  CouplingVariableType.F90
  GoverningEquationBaseType.F90
  MeshType.F90
  MultiPhysicsProbBaseType.F90
  SolverType.F90
  SystemOfEquationsBasePointerType.F90
  SystemOfEquationsBaseType.F90
  UnstructuredGridType.F90
 )

include_directories(${CMAKE_BINARY_DIR}/src/mpp/thermal)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/util)
include_directories(${CMAKE_BINARY_DIR}/src/mpp/vsfm)

include(add_mpp_library)
add_mpp_library(mpp_dtypes ${MPP_DTYPES_SOURCES})
set_target_properties (mpp_dtypes PROPERTIES COMPILE_DEFINITIONS "USE_PETSC_LIB")

set(MPP_LIBRARIES mpp_dtypes;${MPP_LIBRARIES} PARENT_SCOPE)
set(MPP_LIBRARIES mpp_dtypes;${MPP_LIBRARIES})


if (NOT CMAKE_INSTALL_PREFIX STREQUAL "INSTALL_DISABLED")
  install(TARGETS mpp_dtypes DESTINATION lib)
  file(GLOB HEADERS *.h)
  install(FILES ${HEADERS} DESTINATION include/mpp/dtypes)
endif()

