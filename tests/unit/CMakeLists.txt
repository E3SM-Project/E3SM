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

add_mpi_test(piotest
    COMMANDS ${CMAKE_CURRENT_BINARY_DIR}/piotest
    NUMPROCS 4
    TIMEOUT ${DEFAULT_TEST_TIMEOUT})

if (PIO_HDF5_LOGGING) 
    target_compile_definitions (pio_unit_test 
        PRIVATE LOGGING)
endif ()
