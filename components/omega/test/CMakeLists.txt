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

target_link_libraries(${_TestDataTypesName} ${OMEGA_LIB_NAME} yakl parmetis metis)

if(GKlib_FOUND)
  target_link_libraries(${_TestDataTypesName} gklib)
endif()

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

target_link_libraries(${_TestBroadcastName} ${OMEGA_LIB_NAME} spdlog yakl)

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

target_link_libraries(${_TestLoggingName} ${OMEGA_LIB_NAME} spdlog yakl)

add_test(NAME LOGGING_TEST COMMAND ./${_TestLoggingName})

#############
# Decomp test
#############

set(_TestDecompName testDecomp.exe)

add_executable(${_TestDecompName} base/DecompTest.cpp)

target_include_directories(
  ${_TestDecompName}
  PRIVATE
  ${OMEGA_SOURCE_DIR}/src/base
  ${OMEGA_SOURCE_DIR}/src/infra
  ${Parmetis_INCLUDE_DIRS}
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

target_link_libraries(${_TestDecompName} ${OMEGA_LIB_NAME} spdlog yakl parmetis metis pioc)

if(GKlib_FOUND)
  target_link_libraries(${_TestDecompName} gklib)
endif()

add_test(
  NAME DECOMP_TEST
  COMMAND ${MPI_EXEC} -n 8 -- ./${_TestDecompName}
)

##################
# Halo test
##################

set(_TestHaloName testHalo.exe)

# Add halo test
add_executable(${_TestHaloName} base/HaloTest.cpp)

target_include_directories(
  ${_TestHaloName}
  PRIVATE
  ${OMEGA_SOURCE_DIR}/src/base
  ${OMEGA_SOURCE_DIR}/src/infra
  ${Parmetis_INCLUDE_DIRS}
)

target_compile_options(
  ${_TestHaloName}
  PRIVATE
  ${OMEGA_CXX_FLAGS}
)

target_link_libraries(${_TestHaloName} ${OMEGA_LIB_NAME} spdlog yakl parmetis metis pioc)

add_test(
  NAME HALO_TEST
  COMMAND ${MPI_EXEC} -n 8 -- ./${_TestHaloName}
)

################
# HorzMesh test
################

set(_TestHorzMeshName testHorzMesh.exe)

# Add broadcast test
add_executable(${_TestHorzMeshName} ocn/HorzMeshTest.cpp)

target_include_directories(
  ${_TestHorzMeshName}
  PRIVATE
  ${OMEGA_SOURCE_DIR}/src/base
  ${OMEGA_SOURCE_DIR}/src/infra
  ${OMEGA_SOURCE_DIR}/src/ocn
  ${Parmetis_INCLUDE_DIRS}
)

target_compile_options(
  ${_TestHorzMeshName}
  PRIVATE
  ${OMEGA_LINK_OPTIONS}
)

target_link_libraries(${_TestHorzMeshName} ${OMEGA_LIB_NAME} spdlog yakl parmetis metis pioc)

add_test(
  NAME HORZMESH_TEST
  COMMAND ${MPI_EXEC} -n 8 -- ./${_TestHorzMeshName}
)

#############
# IO test
#############

set(_TestIOName testIO.exe)

add_executable(${_TestIOName} base/IOTest.cpp)

target_include_directories(
   ${_TestIOName}
  PRIVATE
  ${OMEGA_SOURCE_DIR}/src/base
  ${OMEGA_SOURCE_DIR}/src/infra
  ${Parmetis_INCLUDE_DIRS}
)

target_compile_options(
   ${_TestIOName}
  PRIVATE
  ${OMEGA_CXX_FLAGS}
)

target_link_options(
   ${_TestIOName}
  PRIVATE
  ${OMEGA_LINK_OPTIONS}
)

target_link_libraries(${_TestIOName}     ${OMEGA_LIB_NAME} spdlog yakl parmetis metis pioc)

if(GKlib_FOUND)
   target_link_libraries(${_TestIOName} gklib)
endif()

add_test(
  NAME IO_TEST
  COMMAND ${MPI_EXEC} -n 8 -- ./${_TestIOName}
)

##################
# Config test
##################

set(_TestConfigName testConfig.exe)

# Add broadcast test
add_executable(${_TestConfigName} infra/ConfigTest.cpp)

target_include_directories(
  ${_TestConfigName}
  PRIVATE
  ${OMEGA_SOURCE_DIR}/src/base
  ${OMEGA_SOURCE_DIR}/src/infra
)

target_compile_options(
  ${_TestConfigName}
  PRIVATE
  ${OMEGA_CXX_FLAGS}
)

target_link_options(
  ${_TestConfigName}
  PRIVATE
  ${OMEGA_LINK_OPTIONS}
)

target_link_libraries(${_TestConfigName} ${OMEGA_LIB_NAME} yaml-cpp)

add_test(
  NAME CONFIG_TEST
  COMMAND ${MPI_EXEC} -n 8 -- ./${_TestConfigName}
)

##################
# Metadata test
##################

set(_TestMetadataName testMetadata.exe)

add_executable(${_TestMetadataName} infra/MetadataTest.cpp)

target_compile_options(
  ${_TestMetadataName}
  PRIVATE
  "-I${OMEGA_SOURCE_DIR}/src/base"
  "-I${OMEGA_SOURCE_DIR}/src/infra"
  ${OMEGA_CXX_FLAGS}
)

target_link_options(
  ${_TestMetadataName}
  PRIVATE
  ${OMEGA_LINK_OPTIONS}
)

target_link_libraries(${_TestMetadataName} ${OMEGA_LIB_NAME} yakl)

add_test(NAME METADATA_TEST COMMAND ./${_TestMetadataName})

##################
# IOField test
##################

set(_TestIOFieldName testIOField.exe)

add_executable(${_TestIOFieldName} infra/IOFieldTest.cpp)

target_compile_options(
   ${_TestIOFieldName}
  PRIVATE
  "-I${OMEGA_SOURCE_DIR}/src/base"
  "-I${OMEGA_SOURCE_DIR}/src/infra"
  ${OMEGA_CXX_FLAGS}
)

target_link_options(
   ${_TestIOFieldName}
  PRIVATE
  ${OMEGA_LINK_OPTIONS}
)

target_link_libraries(${_TestIOFieldName} ${OMEGA_LIB_NAME} yakl)

add_test(NAME IOFIELD_TEST
         COMMAND ${MPI_EXEC} -n 8 -- ./${_TestIOFieldName})

##################
# Time Manager test
##################

set(_TestTimeMgrName testTimeMgr.exe)

add_executable(${_TestTimeMgrName} infra/TimeMgrTest.cpp)

target_include_directories(
   ${_TestTimeMgrName}
  PRIVATE
  ${OMEGA_SOURCE_DIR}/src/base
  ${OMEGA_SOURCE_DIR}/src/infra
)

target_compile_options(
   ${_TestTimeMgrName}
  PRIVATE
  ${OMEGA_CXX_FLAGS}
)

target_link_options(
   ${_TestTimeMgrName}
  PRIVATE
  ${OMEGA_LINK_OPTIONS}
)

target_link_libraries(${_TestTimeMgrName}     ${OMEGA_LIB_NAME} spdlog yakl)

add_test(
  NAME TIMEMGR_TEST
  COMMAND ./${_TestTimeMgrName}
)

##################
# YAKL test
##################

set(_TestYaklName testYakl.exe)

add_executable(
  ${_TestYaklName} ${E3SM_EXTERNALS_ROOT}/YAKL/unit/ParForC/ParForC.cpp
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
  CONFIG_TEST
  DECOMP_TEST
  HALO_TEST
  HORZMESH_TEST
  IO_TEST
  METADATA_TEST
  IOFIELD_TEST
  PROPERTIES FAIL_REGULAR_EXPRESSION "FAIL" PASS_REGULAR_EXPRESSION "PASS"
)

set_tests_properties(
  YAKL_TEST
  PROPERTIES FAIL_REGULAR_EXPRESSION "ERROR"
)
