include (CheckPlatform)
include (LibCTestMPI)

check_platform ()

#==============================================================================
#  PREPARE FOR TESTING
#==============================================================================

file (COPY "./input.nl" 
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file (COPY "./not_netcdf.ieee" 
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

#==============================================================================
#  DEFINE THE TARGETS AND TESTS
#==============================================================================

set (SRCS basic_tests.F90
          driver.F90
          global_vars.F90
          ncdf_tests.F90
          nc_set_log_level2.c)

add_executable (pio_unit_test ${SRCS})
target_link_libraries (pio_unit_test piof pioc)

# Test Timeout (4 min = 240 sec)
set (DEFAULT_TEST_TIMEOUT 240)

add_mpi_test(pio_unit_test
    COMMANDS ${CMAKE_CURRENT_BINARY_DIR}/pio_unit_test
    NUMPROCS 4
    TIMEOUT ${DEFAULT_TEST_TIMEOUT})

if (PIO_HDF5_LOGGING) 
    target_compile_definitions (pio_unit_test 
        PRIVATE LOGGING)
endif ()
