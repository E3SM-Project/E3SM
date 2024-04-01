# Omega Unit Tests

##################
# Data type test
##################

add_executable(testDataTypes.exe base/DataTypesTest.cpp)

target_link_libraries(testDataTypes.exe ${OMEGA_LIB_NAME})

add_test(NAME DATA_TYPES_TEST COMMAND ${MPI_EXEC} --cpu-bind=cores -n 1 -- ./testDataTypes.exe)

##################
# Machine env test
##################

add_executable(testMachEnv.exe base/MachEnvTest.cpp)

target_link_libraries(testMachEnv.exe ${OMEGA_LIB_NAME})

add_test(NAME MACHINE_ENV_TEST COMMAND ${MPI_EXEC} --cpu-bind=cores -n 8 -- ./testMachEnv.exe)


##################
# Broadcast test
##################

add_executable(testBroadcast.exe base/BroadcastTest.cpp)

target_link_libraries(testBroadcast.exe ${OMEGA_LIB_NAME})

add_test(NAME BROADCAST_TEST COMMAND ${MPI_EXEC} --cpu-bind=cores -n 8 -- ./testBroadcast.exe)


##################
# Logging test
##################

add_executable(testLogging.exe infra/LoggingTest.cpp)

target_link_libraries(testLogging.exe ${OMEGA_LIB_NAME})

add_test(NAME LOGGING_TEST COMMAND ./testLogging.exe)

#############
# Decomp test
#############

add_executable(testDecomp.exe base/DecompTest.cpp)

target_link_libraries(testDecomp.exe ${OMEGA_LIB_NAME})

add_test( NAME DECOMP_TEST COMMAND ${MPI_EXEC} -n 8 --cpu-bind=cores -- ./testDecomp.exe)

###################
## Halo test
###################
#
#add_executable(testHalo.exe base/HaloTest.cpp)
#
#target_link_libraries(testHalo.exe ${OMEGA_LIB_NAME})
#
#add_test(NAME HALO_TEST COMMAND ${MPI_EXEC} --cpu-bind=cores -n 8 -- ./testHalo.exe)

################
# HorzMesh test
################

add_executable(testHorzMesh.exe ocn/HorzMeshTest.cpp)

target_link_libraries(testHorzMesh.exe ${OMEGA_LIB_NAME})

add_test(NAME HORZMESH_TEST COMMAND ${MPI_EXEC} -n 8 --cpu-bind=cores -- ./testHorzMesh.exe)

#############
# IO test
#############

add_executable(testIO.exe base/IOTest.cpp)

target_link_libraries(testIO.exe ${OMEGA_LIB_NAME})

add_test(NAME IO_TEST COMMAND ${MPI_EXEC} -n 8 --cpu-bind=cores -- ./testIO.exe)

##################
# Config test
##################

add_executable(testConfig.exe infra/ConfigTest.cpp)

target_link_libraries(testConfig.exe ${OMEGA_LIB_NAME})

add_test(NAME CONFIG_TEST COMMAND ${MPI_EXEC} --cpu-bind=cores -n 8 -- ./testConfig.exe)

##################
# Metadata test
##################

add_executable(testMetadata.exe infra/MetadataTest.cpp)

target_link_libraries(testMetadata.exe ${OMEGA_LIB_NAME})

add_test(NAME METADATA_TEST COMMAND ./testMetadata.exe)

##################
# IOField test
##################

add_executable(testIOField.exe infra/IOFieldTest.cpp)

target_link_libraries(testIOField.exe ${OMEGA_LIB_NAME})

add_test(NAME IOFIELD_TEST COMMAND ${MPI_EXEC} --cpu-bind=cores -n 8 -- ./testIOField.exe)

##################
# Time Manager test
##################

add_executable(testTimeMgr.exe infra/TimeMgrTest.cpp)

target_link_libraries(testTimeMgr.exe ${OMEGA_LIB_NAME})

add_test(NAME TIMEMGR_TEST COMMAND ./testTimeMgr.exe)


##################
# Kokkos test
##################

add_executable(testKokkos.exe infra/OmegaKokkosTest.cpp)

target_link_libraries(testKokkos.exe ${OMEGA_LIB_NAME})

add_test(NAME KOKKOS_TEST COMMAND ./testKokkos.exe)

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
#  HALO_TEST
  HORZMESH_TEST
  IO_TEST
  METADATA_TEST
  TIMEMGR_TEST
  KOKKOS_TEST
  IOFIELD_TEST
  PROPERTIES FAIL_REGULAR_EXPRESSION "FAIL" PASS_REGULAR_EXPRESSION "PASS"
)
