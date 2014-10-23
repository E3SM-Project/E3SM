###-------------------------------------------------------------------------###
### CMakeList.txt for tests of pio
###-------------------------------------------------------------------------###

# Adding PIO definitions - defined in the root directory
ADD_DEFINITIONS(${PIO_DEFINITIONS})

# Include PIO include and lib directories
INCLUDE_DIRECTORIES(${PIO_INCLUDE_DIRS})
LINK_DIRECTORIES(${PIO_LIB_DIR})

# Include test util
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/util)

# PIO Test util sources
SET(SRCS_PIO_TUTIL	${CMAKE_CURRENT_SOURCE_DIR}/util/pio_tutil.F90)

# PIO test source generator path
SET(PIO_TUTIL_GEN_PATH ${CMAKE_CURRENT_SOURCE_DIR}/util)

# Test source files that are generated (from template files)
SET(SRCS_GENERATED pio_file_simple_tests.F90 pio_init_finalize.F90 pio_fail.F90 
    pio_file_fail.F90 ncdf_simple_tests.F90 ncdf_get_put.F90 ncdf_fail.F90 
    pio_decomp_tests.F90)

# Generate test cases from template files
FOREACH(gen_src IN LISTS SRCS_GENERATED)
  ADD_CUSTOM_COMMAND(
    OUTPUT ${gen_src}
    COMMAND ${PIO_TUTIL_GEN_PATH}/pio_tf_f90gen.pl 
    --annotate-source --out=./${gen_src} ${CMAKE_CURRENT_SOURCE_DIR}/${gen_src}.in
    DEPENDS ${gen_src}.in
  )
ENDFOREACH()

if(${PIO_BUILD_TIMING} MATCHES "ON") 
  SET(TIMING_LINK_LIB timing)
endif()

###-------------------------------------------------------------------------###
# add_test usage:
#
# add_mpi_test("test name" "test exe" numPes)
#
###-------------------------------------------------------------------------###

###-------------------------------------------------------------------------###
#  arguments to binarys should be a CMake list not as a CMake string
#  deal with building the list and differences in platform mpiexe dispatch
#  in the function add_mpi_test (can be extended for other platforms
###-------------------------------------------------------------------------###

function( add_mpi_test _testName _testExe _testArgs _numProc )

    if (${PLATFORM} STREQUAL "yellowstone" )
        ###
        ### note: no space between -n and num_proc for mpirun.lsf on
        ### yellowstone
        ###
        set(MPIEXEC_NPF -n${_numProc})
        set(EXE_CMD ${EXECCA} ${MPIEXEC} ${_testExe} ${_testArgs} ${MPIEXEC_NPF})
    else()
        set(MPIEXEC_NPF ${MPIEXEC_NUMPROC_FLAG} ${_numProc})
        set(EXE_CMD ${MPIEXEC} ${MPIEXEC_NPF} ${_testExe} ${_testArgs} ${_testArgs}
    endif()
    add_test(NAME ${_testName} COMMAND ${EXE_CMD} )

endfunction(add_mpi_test)

###-------------------------------------------------------------------------###
# test pio_init_finalize
###-------------------------------------------------------------------------###

ADD_EXECUTABLE(pio_init_finalize pio_init_finalize.F90 ${SRCS_PIO_TUTIL})
TARGET_LINK_LIBRARIES(pio_init_finalize piof pioc ${TIMING_LINK_LIB} )

add_mpi_test("init_finialize_1_proc" "./pio_init_finalize" "" 1)
add_mpi_test("init_finialize_2_proc" "./pio_init_finalize" "" 2)
add_mpi_test("init_finalize_2_proc_with_args" "./pio_init_finalize" "--pio-tf-stride=2 --pio-tf-num-aggregators=2" 2)

###-------------------------------------------------------------------------###
# test pio_file_simple_tests
###-------------------------------------------------------------------------###

ADD_EXECUTABLE(pio_file_simple_tests pio_file_simple_tests.F90 ${SRCS_PIO_TUTIL})
TARGET_LINK_LIBRARIES(pio_file_simple_tests  piof pioc ${TIMING_LINK_LIB} )

add_mpi_test("pio_file_simple_tests" "./pio_file_simple_tests" "" 2)

###-------------------------------------------------------------------------###
# test pio_file_fail
###-------------------------------------------------------------------------###

ADD_EXECUTABLE(pio_file_fail pio_file_fail.F90 ${SRCS_PIO_TUTIL})
TARGET_LINK_LIBRARIES(pio_file_fail  piof pioc ${TIMING_LINK_LIB} )

add_mpi_test("pio_file_fail" "./pio_file_fail" "" 2)

###-------------------------------------------------------------------------###
# test ncdf_simple_tests
###-------------------------------------------------------------------------###

ADD_EXECUTABLE(ncdf_simple_tests ncdf_simple_tests.F90 ${SRCS_PIO_TUTIL})
TARGET_LINK_LIBRARIES(ncdf_simple_tests  piof pioc ${TIMING_LINK_LIB} )

add_mpi_test("ncdf_simple_tests" "./ncdf_simple_tests" "" 2)

###-------------------------------------------------------------------------###
# test ncdf_get_put
###-------------------------------------------------------------------------###

ADD_EXECUTABLE(ncdf_get_put ncdf_get_put.F90 ${SRCS_PIO_TUTIL})
TARGET_LINK_LIBRARIES(ncdf_get_put  piof pioc ${TIMING_LINK_LIB} )

add_mpi_test("ncdf_get_put_1proc" "./ncdf_simple_tests" "" 1)
add_mpi_test("ncdf_get_put_2proc" "./ncdf_simple_tests" "" 2)

###-------------------------------------------------------------------------###
# test ncdf_fail
###-------------------------------------------------------------------------###

ADD_EXECUTABLE(ncdf_fail ncdf_fail.F90 ${SRCS_PIO_TUTIL})
TARGET_LINK_LIBRARIES(ncdf_fail  piof pioc ${TIMING_LINK_LIB} )

add_mpi_test("ncdf_fail" "./ncdf_fail" "" 2)

###-------------------------------------------------------------------------###
# test pio_decomp_tests
###-------------------------------------------------------------------------###

ADD_EXECUTABLE(pio_decomp_tests pio_decomp_tests.F90 ${SRCS_PIO_TUTIL})
TARGET_LINK_LIBRARIES(pio_decomp_tests  piof pioc ${TIMING_LINK_LIB} )

add_mpi_test("pio_decomp_tests_1p" "./pio_decomp_tests" "" 1)
add_mpi_test("pio_decomp_tests_2p" "./pio_decomp_tests" "" 2)
add_mpi_test("pio_decomp_tests_3p" "./pio_decomp_tests" "" 3)

add_mpi_test("pio_decomp_tests_4p_1agg" "./pio_decomp_tests" "--pio-tf-num-aggregators=1" 4)
add_mpi_test("pio_decomp_tests_4p_2agg" "./pio_decomp_tests" "--pio-tf-num-aggregators=2" 4)
add_mpi_test("pio_decomp_tests_4p_3agg" "./pio_decomp_tests" "--pio-tf-num-aggregators=3" 4)

add_mpi_test("pio_decomp_tests_4p_1iop" "./pio_decomp_tests" "--pio-tf-num-io-tasks=1" 4)
add_mpi_test("pio_decomp_tests_4p_2iop" "./pio_decomp_tests" "--pio-tf-num-io-tasks=2" 4)
add_mpi_test("pio_decomp_tests_4p_3iop" "./pio_decomp_tests" "--pio-tf-num-io-tasks=3" 4)

add_mpi_test("pio_decomp_tests_4p_2iop_2str" "./pio_decomp_tests" "--pio-tf-num-io-tasks=2 --pio-tf-stride=2" 4)
add_mpi_test("pio_decomp_tests_4p_2iop_1agg" "./pio_decomp_tests" "--pio-tf-num-io-tasks=2 --pio-tf-num-aggregators=1" 4)
