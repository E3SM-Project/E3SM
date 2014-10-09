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
SET(SRCS_GENERATED pio_file_simple_tests.F90 pio_init_finalize.F90 pio_fail.F90 pio_file_fail.F90 ncdf_simple_tests.F90 ncdf_get_put.F90 ncdf_fail.F90 pio_decomp_tests.F90)

# Generate test cases from template files
FOREACH(gen_src IN LISTS SRCS_GENERATED)
  ADD_CUSTOM_COMMAND(
    OUTPUT ${gen_src}
    COMMAND ${PIO_TUTIL_GEN_PATH}/pio_tf_f90gen.pl --annotate-source --out=./${gen_src} ${CMAKE_CURRENT_SOURCE_DIR}/${gen_src}.in
    DEPENDS ${gen_src}.in
  )
ENDFOREACH()

if(${PIO_BUILD_TIMING} MATCHES "ON") 
  SET(TIMING_LINK_LIB timing)
endif()


ADD_EXECUTABLE(pio_init_finalize pio_init_finalize.F90 ${SRCS_PIO_TUTIL})
TARGET_LINK_LIBRARIES(pio_init_finalize piof pioc ${TIMING_LINK_LIB} )

ADD_TEST(init_finalize_1_proc ${MPIEXEC} ./pio_init_finalize ${MPIEXEC_NUMPROC_FLAG} 1)
ADD_TEST(init_finalize_2_proc ${MPIEXEC} ./pio_init_finalize ${MPIEXEC_NUMPROC_FLAG} 2)
ADD_TEST(init_finalize_2_proc_with_args ${MPIEXEC} ./pio_init_finalize --pio-tf-stride=2 --pio-tf-num-aggregators=2 ${MPIEXEC_NUMPROC_FLAG} 2)

ADD_EXECUTABLE(pio_file_simple_tests pio_file_simple_tests.F90 ${SRCS_PIO_TUTIL})
TARGET_LINK_LIBRARIES(pio_file_simple_tests  piof pioc ${TIMING_LINK_LIB} )

ADD_TEST(pio_file_simple_tests ${MPIEXEC} ./pio_file_simple_tests ${MPIEXEC_NUMPROC_FLAG} 2)

ADD_EXECUTABLE(pio_file_fail pio_file_fail.F90 ${SRCS_PIO_TUTIL})
TARGET_LINK_LIBRARIES(pio_file_fail  piof pioc ${TIMING_LINK_LIB} )

ADD_TEST(pio_file_fail ${MPIEXEC} ./pio_file_fail ${MPIEXEC_NUMPROC_FLAG} 2)

ADD_EXECUTABLE(ncdf_simple_tests ncdf_simple_tests.F90 ${SRCS_PIO_TUTIL})
TARGET_LINK_LIBRARIES(ncdf_simple_tests  piof pioc ${TIMING_LINK_LIB} )

ADD_TEST(ncdf_simple_tests ${MPIEXEC} ./ncdf_simple_tests ${MPIEXEC_NUMPROC_FLAG} 2)

ADD_EXECUTABLE(ncdf_get_put ncdf_get_put.F90 ${SRCS_PIO_TUTIL})
TARGET_LINK_LIBRARIES(ncdf_get_put  piof pioc ${TIMING_LINK_LIB} )

ADD_TEST(ncdf_get_put_1proc ${MPIEXEC} ./ncdf_get_put ${MPIEXEC_NUMPROC_FLAG} 1)
ADD_TEST(ncdf_get_put_2proc ${MPIEXEC} ./ncdf_get_put ${MPIEXEC_NUMPROC_FLAG} 2)

ADD_EXECUTABLE(ncdf_fail ncdf_fail.F90 ${SRCS_PIO_TUTIL})
TARGET_LINK_LIBRARIES(ncdf_fail  piof pioc ${TIMING_LINK_LIB} )

ADD_TEST(ncdf_fail ${MPIEXEC} ./ncdf_fail ${MPIEXEC_NUMPROC_FLAG} 2 )

ADD_EXECUTABLE(pio_decomp_tests pio_decomp_tests.F90 ${SRCS_PIO_TUTIL})
TARGET_LINK_LIBRARIES(pio_decomp_tests  piof pioc ${TIMING_LINK_LIB} )

ADD_TEST(pio_decomp_tests_1p ${MPIEXEC} ./pio_decomp_tests ${MPIEXEC_NUMPROC_FLAG} 1 )
ADD_TEST(pio_decomp_tests_2p ${MPIEXEC} ./pio_decomp_tests ${MPIEXEC_NUMPROC_FLAG} 2)
ADD_TEST(pio_decomp_tests_3p ${MPIEXEC} ./pio_decomp_tests ${MPIEXEC_NUMPROC_FLAG} 3)
ADD_TEST(pio_decomp_tests_4p_1agg ${MPIEXEC} ./pio_decomp_tests --pio-tf-num-aggregators=1 ${MPIEXEC_NUMPROC_FLAG} 4)
ADD_TEST(pio_decomp_tests_4p_2agg ${MPIEXEC} ./pio_decomp_tests --pio-tf-num-aggregators=2 ${MPIEXEC_NUMPROC_FLAG} 4)
ADD_TEST(pio_decomp_tests_4p_3agg ${MPIEXEC} ./pio_decomp_tests --pio-tf-num-aggregators=3 ${MPIEXEC_NUMPROC_FLAG} 4 )
ADD_TEST(pio_decomp_tests_4p_1iop ${MPIEXEC} ./pio_decomp_tests --pio-tf-num-io-tasks=1 {MPIEXEC_NUMPROC_FLAG} 4)
ADD_TEST(pio_decomp_tests_4p_2iop ${MPIEXEC} ./pio_decomp_tests --pio-tf-num-io-tasks=2 ${MPIEXEC_NUMPROC_FLAG} 4)
ADD_TEST(pio_decomp_tests_4p_3iop ${MPIEXEC} ./pio_decomp_tests --pio-tf-num-io-tasks=3 ${MPIEXEC_NUMPROC_FLAG} 4)
ADD_TEST(pio_decomp_tests_4p_2iop_2str ${MPIEXEC} ./pio_decomp_tests --pio-tf-num-io-tasks=2 --pio-tf-stride=2 ${MPIEXEC_NUMPROC_FLAG} 4)
ADD_TEST(pio_decomp_tests_4p_2iop_1agg ${MPIEXEC} ./pio_decomp_tests --pio-tf-num-io-tasks=2 --pio-tf-num-aggregators=1 ${MPIEXEC_NUMPROC_FLAG} 4 )
