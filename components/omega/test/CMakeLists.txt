# Omega Unit Tests

##################
# Data type test
##################

set(_TestDataTypesName testDataTypes.exe)

add_executable(${_TestDataTypesName} base/DataTypesTest.cpp)

target_compile_options(
  ${_TestDataTypesName}
  PRIVATE
  "-I${OMEGA_SOURCE_DIR}/src/base"
  "-I${OMEGA_SOURCE_DIR}/src/infra"
  ${OMEGA_CXX_FLAGS}
)

target_link_options(
  ${_TestDataTypesName}
  PRIVATE
  ${OMEGA_LINK_OPTIONS}
)

target_link_libraries(${_TestDataTypesName} ${OMEGA_LIB_NAME} yakl metis)

add_test(
  NAME DATA_TYPES_TEST
  COMMAND ${MPI_EXEC} -n 1 -- ./${_TestDataTypesName}
)


##################
# Machine env test
##################

set(_TestMachEnvName testMachEnv.exe)

add_executable(${_TestMachEnvName} base/MachEnvTest.cpp)

target_compile_options(
  ${_TestMachEnvName}
  PRIVATE
  "-I${OMEGA_SOURCE_DIR}/src/base"
  "-I${OMEGA_SOURCE_DIR}/src/infra"
  ${OMEGA_CXX_FLAGS}
)

target_link_options(
  ${_TestMachEnvName}
  PRIVATE
  ${OMEGA_LINK_OPTIONS}
)

target_link_libraries(${_TestMachEnvName} ${OMEGA_LIB_NAME} yakl)

add_test(
  NAME MACHINE_ENV_TEST
  COMMAND ${MPI_EXEC} -n 8 -- ./${_TestMachEnvName}
)


##################
# Broadcast test
##################

set(_TestBroadcastName testBroadcast.exe)

# Add broadcast test
add_executable(${_TestBroadcastName} base/BroadcastTest.cpp)

target_compile_options(
  ${_TestBroadcastName}
  PRIVATE
  "-I${OMEGA_SOURCE_DIR}/src/base"
  "-I${OMEGA_SOURCE_DIR}/src/infra"
  ${OMEGA_CXX_FLAGS}
)

target_link_options(
  ${_TestBroadcastName}
  PRIVATE
  ${OMEGA_LINK_OPTIONS}
)

target_link_libraries(${_TestBroadcastName} ${OMEGA_LIB_NAME} yakl)

add_test(
  NAME BROADCAST_TEST
  COMMAND ${MPI_EXEC} -n 8 -- ./${_TestBroadcastName}
)


##################
# Logging test
##################

set(_TestLoggingName testLogging.exe)

add_executable(${_TestLoggingName} infra/LoggingTest.cpp)

target_compile_options(
  ${_TestLoggingName}
  PRIVATE
  "-I${OMEGA_SOURCE_DIR}/src/base"
  "-I${OMEGA_SOURCE_DIR}/src/infra"
  ${OMEGA_CXX_FLAGS}
)

target_link_options(
  ${_TestLoggingName}
  PRIVATE
  ${OMEGA_LINK_OPTIONS}
)

target_link_libraries(${_TestLoggingName} ${OMEGA_LIB_NAME} yakl)

add_test(NAME LOGGING_TEST COMMAND ./${_TestLoggingName})

#############
# Decomp test
#############

set(_TestDecompName testDecomp.exe)

# Add broadcast test
add_executable(${_TestDecompName} base/DecompTest.cpp)

target_include_directories(
  ${_TestDecompName}
  PRIVATE
  ${OMEGA_SOURCE_DIR}/src/base
  ${OMEGA_SOURCE_DIR}/src/infra
  ${PARMETIS_ROOT}/include
)

target_compile_options(
  ${_TestDecompName}
  PRIVATE
  ${OMEGA_CXX_FLAGS}
)

target_link_options(
  ${_TestDecompName}
  PRIVATE
  ${OMEGA_LINK_OPTIONS}
)

target_link_libraries(${_TestDecompName} ${OMEGA_LIB_NAME} yakl metis)

add_test(
  NAME DECOMP_TEST
  COMMAND ${MPI_EXEC} -n 8 -- ./${_TestDecompName}
)


##################
# YAKL test
##################

set(_TestYaklName testYakl.exe)

add_executable(
  ${_TestYaklName} ${E3SM_EXTERNALS_ROOT}/YAKL/unit/performance/performance.cpp
)

include(${E3SM_EXTERNALS_ROOT}/YAKL/yakl_utils.cmake)
yakl_process_target(${_TestYaklName})

# handles cuda cases
if (YAKL_ARCH STREQUAL "CUDA")
  set_target_properties(${_TestYaklName} PROPERTIES LINKER_LANGUAGE CXX)
  if (CMAKE_VERSION VERSION_GREATER "3.18.0")
    set_target_properties(${_TestYaklName} PROPERTIES CUDA_ARCHITECTURES OFF)
  endif()
endif()

add_test(NAME YAKL_TEST COMMAND ./${_TestYaklName})


##################
# test properties
##################

set_tests_properties(
  DATA_TYPES_TEST
  MACHINE_ENV_TEST
  BROADCAST_TEST
  LOGGING_TEST
  YAKL_TEST
  PROPERTIES FAIL_REGULAR_EXPRESSION "FAIL"
)
