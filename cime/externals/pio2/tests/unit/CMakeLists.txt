include (LibMPI)

# Compiler-specific compiler options
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c99")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -c99")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -c99")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c99")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c99")
endif()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")

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
	
add_executable (pio_unit_test EXCLUDE_FROM_ALL ${SRCS})
target_link_libraries (pio_unit_test piof)
if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    target_compile_options (pio_unit_test
        PRIVATE -ffree-line-length-none)
endif()

if (PIO_ENABLE_ASYNC AND NOT PIO_USE_MPI_SERIAL)
  add_executable (test_intercomm EXCLUDE_FROM_ALL test_intercomm.c)
  target_link_libraries (test_intercomm pioc)
  add_dependencies (tests test_intercomm)
  add_executable (test_darray_async EXCLUDE_FROM_ALL test_darray_async.c)
  target_link_libraries (test_darray_async pioc)
  add_dependencies (tests test_darray_async)
endif ()
add_executable (test_names EXCLUDE_FROM_ALL test_names.c)
target_link_libraries (test_names pioc)
add_executable (test_nc4 EXCLUDE_FROM_ALL test_nc4.c)
target_link_libraries (test_nc4 pioc)
add_executable (test_darray EXCLUDE_FROM_ALL test_darray.c)
target_link_libraries (test_darray pioc)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "NAG")
   set ( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mismatch_all" )
#    target_compile_options (gptl
#        PRIVATE -mismatch_all)
endif ()

add_dependencies (tests test_names)
add_dependencies (tests test_darray)
add_dependencies (tests test_nc4)
add_dependencies (tests pio_unit_test)

# Test Timeout (4 min = 240 sec)
set (DEFAULT_TEST_TIMEOUT 240)

if (PIO_USE_MPISERIAL)
    add_test(NAME test_names
        COMMAND test_names)
    add_test(NAME test_nc4
        COMMAND test_nc4)
    add_test(NAME pio_unit_test
        COMMAND pio_unit_test)
    set_tests_properties(pio_unit_test
        PROPERTIES TIMEOUT ${DEFAULT_TEST_TIMEOUT})
else ()
    add_mpi_test(test_darray
        EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/test_darray
        NUMPROCS 4
        TIMEOUT ${DEFAULT_TEST_TIMEOUT})
  if (PIO_ENABLE_ASYNC)
    add_mpi_test(test_intercomm
        EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/test_intercomm
        NUMPROCS 4
        TIMEOUT ${DEFAULT_TEST_TIMEOUT})
  endif ()
    add_mpi_test(test_names
        EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/test_names
        NUMPROCS 4
        TIMEOUT ${DEFAULT_TEST_TIMEOUT})
    add_mpi_test(test_nc4
        EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/test_nc4
        NUMPROCS 4
        TIMEOUT ${DEFAULT_TEST_TIMEOUT})
    add_mpi_test(pio_unit_test
        EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/pio_unit_test
        NUMPROCS 4
        TIMEOUT ${DEFAULT_TEST_TIMEOUT})
endif ()

if (PIO_HDF5_LOGGING) 
    target_compile_definitions (pio_unit_test 
        PUBLIC LOGGING)
endif ()
