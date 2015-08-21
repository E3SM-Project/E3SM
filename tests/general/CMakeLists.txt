include (LibMPI)

#==============================================================================
#  GENERATE TARGET SOURCES
#==============================================================================

SET(GENERATED_SRCS pio_file_simple_tests.F90
                   pio_init_finalize.F90
                   pio_fail.F90 
                   pio_file_fail.F90
                   ncdf_simple_tests.F90
                   ncdf_get_put.F90
                   ncdf_fail.F90 
                   pio_decomp_tests.F90
                   pio_decomp_fillval.F90)

foreach (SRC_FILE IN LISTS GENERATED_SRCS)
    add_custom_command (
        OUTPUT ${SRC_FILE}
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/util/pio_tf_f90gen.pl
                --annotate-source
                --out=${CMAKE_CURRENT_BINARY_DIR}/${SRC_FILE}
                ${CMAKE_CURRENT_SOURCE_DIR}/${SRC_FILE}.in
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${SRC_FILE}.in)
endforeach ()

#==============================================================================
#  DEFINE THE TARGETS AND TESTS
#==============================================================================

# Test Timeout (4 min = 240 sec)
set (DEFAULT_TEST_TIMEOUT 240)

#===== pio_init_finalize =====
add_executable (pio_init_finalize EXCLUDE_FROM_ALL 
    pio_init_finalize.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/util/pio_tutil.F90)
target_link_libraries (pio_init_finalize piof)
add_dependencies (tests pio_init_finalize)

add_mpi_test(init_finialize_1_proc
    EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/pio_init_finalize
    NUMPROCS 1
    TIMEOUT ${DEFAULT_TEST_TIMEOUT})
             
add_mpi_test(init_finialize_2_proc 
    EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/pio_init_finalize 
    NUMPROCS 2 
    TIMEOUT ${DEFAULT_TEST_TIMEOUT})
             
add_mpi_test(init_finalize_2_proc_with_args
    EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/pio_init_finalize
    ARGUMENTS --pio-tf-stride=2 --pio-tf-num-aggregators=2
    NUMPROCS 2 
    TIMEOUT ${DEFAULT_TEST_TIMEOUT})

#===== pio_file_simple_tests =====
add_executable (pio_file_simple_tests EXCLUDE_FROM_ALL
    pio_file_simple_tests.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/util/pio_tutil.F90)
target_link_libraries (pio_file_simple_tests piof)
add_dependencies (tests pio_file_simple_tests)

add_mpi_test(pio_file_simple_tests
    EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/pio_file_simple_tests
    NUMPROCS 2 
    TIMEOUT ${DEFAULT_TEST_TIMEOUT})
             
#===== pio_file_fail =====
add_executable (pio_file_fail EXCLUDE_FROM_ALL
    pio_file_fail.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/util/pio_tutil.F90)
target_link_libraries (pio_file_fail piof)
add_dependencies (tests pio_file_fail)

add_mpi_test(pio_file_fail
    EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/pio_file_fail
    NUMPROCS 2 
    TIMEOUT ${DEFAULT_TEST_TIMEOUT})

#===== ncdf_simple_tests =====
add_executable (ncdf_simple_tests EXCLUDE_FROM_ALL
    ncdf_simple_tests.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/util/pio_tutil.F90)
target_link_libraries (ncdf_simple_tests piof)
add_dependencies (tests ncdf_simple_tests)

add_mpi_test(ncdf_simple_tests
    EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ncdf_simple_tests
    NUMPROCS 2 
    TIMEOUT ${DEFAULT_TEST_TIMEOUT})

#===== ncdf_get_put =====
add_executable (ncdf_get_put EXCLUDE_FROM_ALL
    ncdf_get_put.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/util/pio_tutil.F90)
target_link_libraries (ncdf_get_put piof)
add_dependencies (tests ncdf_get_put)

add_mpi_test(ncdf_get_put_1proc
    EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ncdf_simple_tests
    NUMPROCS 1 
    TIMEOUT ${DEFAULT_TEST_TIMEOUT})
add_mpi_test(ncdf_get_put_2proc
    EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ncdf_simple_tests
    NUMPROCS 2 
    TIMEOUT ${DEFAULT_TEST_TIMEOUT})

#===== ncdf_fail =====
add_executable (ncdf_fail EXCLUDE_FROM_ALL
    ncdf_fail.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/util/pio_tutil.F90)
target_link_libraries (ncdf_fail piof)
add_dependencies (tests ncdf_fail)

add_mpi_test(ncdf_fail
    EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/ncdf_fail
    NUMPROCS 2 
    TIMEOUT ${DEFAULT_TEST_TIMEOUT})
    
#===== pio_decomp_tests =====
add_executable (pio_decomp_tests EXCLUDE_FROM_ALL
    pio_decomp_tests.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/util/pio_tutil.F90)
target_link_libraries (pio_decomp_tests piof)
add_dependencies (tests pio_decomp_tests)

add_mpi_test(pio_decomp_tests_1p
    EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/pio_decomp_tests
    NUMPROCS 1
    TIMEOUT ${DEFAULT_TEST_TIMEOUT})
add_mpi_test(pio_decomp_tests_2p
    EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/pio_decomp_tests
    NUMPROCS 2 
    TIMEOUT ${DEFAULT_TEST_TIMEOUT})
add_mpi_test(pio_decomp_tests_3p
    EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/pio_decomp_tests
    NUMPROCS 3 
    TIMEOUT ${DEFAULT_TEST_TIMEOUT})

add_mpi_test(pio_decomp_tests_4p_1agg
    EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/pio_decomp_tests 
    ARGUMENTS --pio-tf-num-aggregators=1
    NUMPROCS 4 
    TIMEOUT ${DEFAULT_TEST_TIMEOUT})
add_mpi_test(pio_decomp_tests_4p_2agg
    EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/pio_decomp_tests 
    ARGUMENTS --pio-tf-num-aggregators=2
    NUMPROCS 4 
    TIMEOUT ${DEFAULT_TEST_TIMEOUT})
add_mpi_test(pio_decomp_tests_4p_3agg
    EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/pio_decomp_tests
    ARGUMENTS --pio-tf-num-aggregators=3
    NUMPROCS 4 
    TIMEOUT ${DEFAULT_TEST_TIMEOUT})

add_mpi_test(pio_decomp_tests_4p_1iop
    EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/pio_decomp_tests 
    ARGUMENTS --pio-tf-num-io-tasks=1
    NUMPROCS 4 
    TIMEOUT ${DEFAULT_TEST_TIMEOUT})
add_mpi_test(pio_decomp_tests_4p_2iop
    EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/pio_decomp_tests 
    ARGUMENTS --pio-tf-num-io-tasks=2
    NUMPROCS 4 
    TIMEOUT ${DEFAULT_TEST_TIMEOUT})
add_mpi_test(pio_decomp_tests_4p_3iop
    EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/pio_decomp_tests 
    ARGUMENTS --pio-tf-num-io-tasks=3
    NUMPROCS 4 
    TIMEOUT ${DEFAULT_TEST_TIMEOUT})

add_mpi_test(pio_decomp_tests_4p_2iop_2str
    EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/pio_decomp_tests 
    ARGUMENTS --pio-tf-num-io-tasks=2 --pio-tf-stride=2
    NUMPROCS 4 
    TIMEOUT ${DEFAULT_TEST_TIMEOUT})
add_mpi_test(pio_decomp_tests_4p_2iop_1agg
    EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/pio_decomp_tests 
    ARGUMENTS --pio-tf-num-io-tasks=2 --pio-tf-num-aggregators=1
    NUMPROCS 4 
    TIMEOUT ${DEFAULT_TEST_TIMEOUT})

#===== pio_decomp_fillval =====
add_executable (pio_decomp_fillval EXCLUDE_FROM_ALL
    pio_decomp_fillval.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/util/pio_tutil.F90)
target_link_libraries (pio_decomp_fillval piof)
add_dependencies (tests pio_decomp_fillval)

add_mpi_test(pio_decomp_fillval
    EXECUTABLE ${CMAKE_CURRENT_BINARY_DIR}/pio_decomp_fillval
    NUMPROCS 4
    TIMEOUT ${DEFAULT_TEST_TIMEOUT})
