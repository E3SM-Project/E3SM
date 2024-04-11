# Omega Unit Tests


#####################
# Omega Test Function
#####################

function(add_omega_test test_name exe_name source_files mpi_args)

  # Create the executable
  add_executable(${exe_name} ${source_files})

  # Link the library
  target_link_libraries(
    ${exe_name}
    PRIVATE
    ${OMEGA_LIB_NAME}
    OmegaLibFlags
  )

  # Add the test command
  add_test(
    NAME ${test_name} 
    COMMAND ${MPI_EXEC} ${mpi_args} -- ./${exe_name}
  )

endfunction()


##################
# Data type test
##################

add_omega_test(
    DATA_TYPES_TEST
    testDataTypes.exe
    base/DataTypesTest.cpp
    "-n 1;--cpu-bind=cores"
)

##################
# Machine env test
##################

add_omega_test(
    MACHINE_ENV_TEST
    testMachEnv.exe
    base/MachEnvTest.cpp
    "-n 8;--cpu-bind=cores"
)

##################
# Broadcast test
##################

add_omega_test(
    BROADCAST_TEST
    testBroadcast.exe
    base/BroadcastTest.cpp
    "-n 8;--cpu-bind=cores"
)

##################
# Logging test
##################

add_omega_test(
    LOGGING_TEST
    testLogging.exe
    infra/LoggingTest.cpp
    "-n 8;--cpu-bind=cores"
)

#############
# Decomp test
#############

add_omega_test(
    DECOMP_TEST
    testDecomp.exe
    base/DecompTest.cpp
    "-n 8;--cpu-bind=cores"
)

##################
# Halo test
##################

add_omega_test(
    HALO_TEST
    testHalo.exe
    base/HaloTest.cpp
    "-n 8;--cpu-bind=cores"
)

################
# HorzMesh test
################

add_omega_test(
    HORZMESH_TEST
    testHorzMesh.exe
    ocn/HorzMeshTest.cpp
    "-n 8;--cpu-bind=cores"
)

#############
# IO test
#############

add_omega_test(
    IO_TEST
    testIO.exe
    base/IOTest.cpp
    "-n 8;--cpu-bind=cores"
)

##################
# Config test
##################

add_omega_test(
    CONFIG_TEST
    testConfig.exe
    infra/ConfigTest.cpp
    "-n 8;--cpu-bind=cores"
)

##################
# Metadata test
##################

add_omega_test(
    METADATA_TEST
    testMetadata.exe
    infra/MetadataTest.cpp
    "-n 1;--cpu-bind=cores"
)

##################
# IOField test
##################

add_omega_test(
    IOFIELD_TEST
    testIOField.exe
    infra/IOFieldTest.cpp
    "-n 8;--cpu-bind=cores"
)

##################
# Time Manager test
##################

add_omega_test(
    TIMEMGR_TEST
    testTimeMgr.exe
    infra/TimeMgrTest.cpp
    "-n 1;--cpu-bind=cores"
)

##################
# Kokkos test
##################

add_omega_test(
    KOKKOS_TEST
    testKokkos.exe
    infra/OmegaKokkosTest.cpp
    "-n 1;--cpu-bind=cores"
)

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
  TIMEMGR_TEST
  KOKKOS_TEST
  IOFIELD_TEST
  PROPERTIES FAIL_REGULAR_EXPRESSION "FAIL" PASS_REGULAR_EXPRESSION "PASS"
)
