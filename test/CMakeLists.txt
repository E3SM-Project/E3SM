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
# ADD_TEST(<test_name>
#          <queue, whichMpirun, numpes>
#          <exe+arguments>)
#
#
# ADD_TEST(pio_decomp_tests_4p_2iop_2str
#          ${MPIEXEC_ARG1} ${MPIEXEC_ARG2}
#          ./pio_decomp_tests --pio-tf-num-io-tasks=2 --pio-tf-stride=2)
#
###-------------------------------------------------------------------------###

###-------------------------------------------------------------------------###
# set -np <numProcs> for darwin and set -n <numProcs> for yellowstone.
# so this once for selected pe counts.  mpirun.lsf on yellowstone wants the 
# -n <proc number> at the end of the command and we provide mpirunLsfReorderArgs.py
# as a wrapper that does this reordering so we can keep this file relatively simple.
###-------------------------------------------------------------------------###

#if (${PLATFORM} STREQUAL "yellowstone" )
#    set(MPIEXEC_NPF1 -n1 )
#    set(MPIEXEC_NPF2 -n2 )
#    set(MPIEXEC_NPF3 -n3 )
#    set(MPIEXEC_NPF4 -n4 )
#    set (EXECCA "execca")
#else()
#    set(MPIEXEC_NPF1 ${MPIEXEC_NUMPROC_FLAG} 1)
#    set(MPIEXEC_NPF2 ${MPIEXEC_NUMPROC_FLAG} 2)
#    set(MPIEXEC_NPF3 ${MPIEXEC_NUMPROC_FLAG} 3)
#    set(MPIEXEC_NPF4 ${MPIEXEC_NUMPROC_FLAG} 4)
#endif()

###-------------------------------------------------------------------------###
#  arguments to binarys should be a CMake list not as a CMake string
###-------------------------------------------------------------------------###


function( add_mpi_test test_name test_exe num_proc )

    if (${PLATFORM} STREQUAL "yellowstone" )
        ###
        ### note: no space between -n and num_proc
        ###
        set(MPIEXEC_NPF -n${num_proc})
        set(EXE_CMD ${EXECCA} ${MPIEXEC} ${test_exe} ${MPIEXEC_NPF})
    else()
        set(MPIEXEC_NPF ${MPIEXEC_NUMPROC_FLAG} ${num_proc})
        set(EXE_CMD ${MPIEXEC} ${MPIEXEC_NPF} ${test_exe})
    endif()
    add_test(NAME ${test_name} COMMAND ${EXE_CMD} )

endfunction(add_mpi_test)

###-------------------------------------------------------------------------###
# test pio_init_finalize
###-------------------------------------------------------------------------###


ADD_EXECUTABLE(pio_init_finalize pio_init_finalize.F90 ${SRCS_PIO_TUTIL})
TARGET_LINK_LIBRARIES(pio_init_finalize piof pioc ${TIMING_LINK_LIB} )


#set (test_exe ./pio_init_finalize)
#ADD_TEST(init_finalize_1_proc
#         ${MPIEXEC} ${MPIEXEC_NPF1}
#         ${test_exe}  )


#add_test(NAME init_finalize_1_proc
    #COMMAND bsub -Is -q caldera -PP93300070 -W2:0 ${MPIEXEC_NPF1} ${test_exe} )




add_mpi_test("init_finialize_1_proc_func" "./pio_init_finalize" 1)

set (test_exe "./pio_init_finalize")
if (${PLATFORM} STREQUAL "yellowstone" )
    set(EXE_CMD ${EXECCA} ${MPIEXEC} ${test_exe} ${MPIEXEC_NPF1})
else()
    set(EXE_CMD ${MPIEXEC} ${MPIEXEC_NPF1} ${test_exe})
endif()
add_test(NAME init_finalize_1_proc COMMAND ${EXE_CMD} )






#add_test(NAME init_finalize_1_proc COMMAND ${EXECCA} ${MPIEXEC} ${test_exe} ${MPIEXEC_NPF1}  )

#ADD_TEST(init_finalize_2_proc
#         ${MPIEXEC} ${MPIEXEC_NPF2}
#         ${test_exe}  )
  
  

#add_test(NAME init_finalize_2_proc
#    COMMAND bsub -Is -q caldera -PP93300070 -W2:0 -n2 mpirun.lsf ./pio_init_finalize )

#if (${PLATFORM} STREQUAL "yellowstone" )
#   set (test_exe "./pio_init_finalize --pio-tf-stride=2 --pio-tf-num-aggregators=2")
#else()
#   set (test_exe ./pio_init_finalize --pio-tf-stride=2 --pio-tf-num-aggregators=2)
#endif()
#ADD_TEST(init_finalize_2_proc_with_args
#         ${MPIEXEC} ${MPIEXEC_NPF2}
#         ${test_exe}  )

###-------------------------------------------------------------------------###
# test pio_file_simple_tests
###-------------------------------------------------------------------------###

#ADD_EXECUTABLE(pio_file_simple_tests pio_file_simple_tests.F90 ${SRCS_PIO_TUTIL})
#TARGET_LINK_LIBRARIES(pio_file_simple_tests  piof pioc ${TIMING_LINK_LIB} )

#set (test_exe ./pio_file_simple_tests)
#ADD_TEST(pio_file_simple_tests
#         ${MPIEXEC} ${MPIEXEC_NPF2}
#         ${test_exe} )

###-------------------------------------------------------------------------###
# test pio_file_fai
###-------------------------------------------------------------------------###

#ADD_EXECUTABLE(pio_file_fail pio_file_fail.F90 ${SRCS_PIO_TUTIL})
#TARGET_LINK_LIBRARIES(pio_file_fail  piof pioc ${TIMING_LINK_LIB} )

#set (test_exe ./pio_file_fail)
#ADD_TEST(pio_file_fail
#         ${MPIEXEC} ${MPIEXEC_NPF2}
#         ${test_exe} )

###-------------------------------------------------------------------------###
# test ncdf_simple_tests
###-------------------------------------------------------------------------###

#ADD_EXECUTABLE(ncdf_simple_tests ncdf_simple_tests.F90 ${SRCS_PIO_TUTIL})
#TARGET_LINK_LIBRARIES(ncdf_simple_tests  piof pioc ${TIMING_LINK_LIB} )

#set (test_exe ./ncdf_simple_tests)
#ADD_TEST(ncdf_simple_tests
#        ${MPIEXEC} ${MPIEXEC_NPF2}
#        ${test_exe} )

###-------------------------------------------------------------------------###
# test ncdf_get_put
###-------------------------------------------------------------------------###

#ADD_EXECUTABLE(ncdf_get_put ncdf_get_put.F90 ${SRCS_PIO_TUTIL})
#TARGET_LINK_LIBRARIES(ncdf_get_put  piof pioc ${TIMING_LINK_LIB} )

#set (test_exe ./ncdf_get_put)
#ADD_TEST(ncdf_get_put_1proc
#         ${MPIEXEC} ${MPIEXEC_NPF1}
#         ${test_exe} )

#ADD_TEST(ncdf_get_put_2proc
#         ${MPIEXEC} ${MPIEXEC_NPF2}
#         ${test_exe} )

###-------------------------------------------------------------------------###
# test ncdf_fail
###-------------------------------------------------------------------------###

#ADD_EXECUTABLE(ncdf_fail ncdf_fail.F90 ${SRCS_PIO_TUTIL})
#TARGET_LINK_LIBRARIES(ncdf_fail  piof pioc ${TIMING_LINK_LIB} )

#set (test_exe ./ncdf_fail)
#ADD_TEST(ncdf_fail
#         ${MPIEXEC} ${MPIEXEC_NPF2}
#         ${test_exe} )

###-------------------------------------------------------------------------###
# test pio_decomp_tests
###-------------------------------------------------------------------------###

#ADD_EXECUTABLE(pio_decomp_tests pio_decomp_tests.F90 ${SRCS_PIO_TUTIL})
#TARGET_LINK_LIBRARIES(pio_decomp_tests  piof pioc ${TIMING_LINK_LIB} )

#set (test_exe ./pio_decomp_tests)
#ADD_TEST(pio_decomp_tests_1p
 #        ${MPIEXEC} ${MPIEXEC_NPF1}
 #        ${test_exe} )
#ADD_TEST(pio_decomp_tests_2p
#         ${MPIEXEC} ${MPIEXEC_NPF2}
#         ${test_exe} )
#ADD_TEST(pio_decomp_tests_3p
#         ${MPIEXEC} ${MPIEXEC_NPF3}
#         ${test_exe} )

#if (${PLATFORM} STREQUAL "yellowstone" )
#   set (test_exe "./pio_decomp_tests --pio-tf-num-aggregators=1")
#else()
#   set (test_exe ./pio_decomp_tests --pio-tf-num-aggregators=1)
#endif()
#ADD_TEST(pio_decomp_tests_4p_1agg
#         ${MPIEXEC} ${MPIEXEC_NPF4}
#         ${test_exe} )

#if (${PLATFORM} STREQUAL "yellowstone" )
#   set (test_exe "./pio_decomp_tests --pio-tf-num-aggregators=2")
#else()
#   set (test_exe ./pio_decomp_tests --pio-tf-num-aggregators=2)
#endif()
#ADD_TEST(pio_decomp_tests_4p_2agg
#         ${MPIEXEC} ${MPIEXEC_NPF4}
#         ${test_exe} )

#if (${PLATFORM} STREQUAL "yellowstone" )
#   set (test_exe "./pio_decomp_tests --pio-tf-num-aggregators=3")
#else()
#   set (test_exe ./pio_decomp_tests --pio-tf-num-aggregators=3)
#endif()
#ADD_TEST(pio_decomp_tests_4p_3agg
#         ${MPIEXEC} ${MPIEXEC_NPF4}
#         ${test_exe} )#

#if (${PLATFORM} STREQUAL "yellowstone" )
#   set (test_exe "./pio_decomp_tests --pio-tf-num-io-tasks=1")
#else()
#   set (test_exe ./pio_decomp_tests --pio-tf-num-io-tasks=1)
#endif()
#ADD_TEST(pio_decomp_tests_4p_1iop
#         ${MPIEXEC} ${MPIEXEC_NPF4}
#         ${test_exe} )

#if (${PLATFORM} STREQUAL "yellowstone" )
#   set (test_exe "./pio_decomp_tests --pio-tf-num-io-tasks=2")
#else()
#   set (test_exe ./pio_decomp_tests --pio-tf-num-io-tasks=2)
#endif()
#ADD_TEST(pio_decomp_tests_4p_2iop
#         ${MPIEXEC} ${MPIEXEC_NPF4}
#         ${test_exe})

#if (${PLATFORM} STREQUAL "yellowstone" )
#   set (test_exe "./pio_decomp_tests --pio-tf-num-io-tasks=3")
#else()
#   set (test_exe ./pio_decomp_tests --pio-tf-num-io-tasks=3)
#endif()
#ADD_TEST(pio_decomp_tests_4p_3iop
#         ${MPIEXEC} ${MPIEXEC_NPF4}
#         ${test_exe})

#if (${PLATFORM} STREQUAL "yellowstone" )
#   set (test_exe "./pio_decomp_tests --pio-tf-num-io-tasks=2 --pio-tf-stride=2")
#else()
#   set (test_exe ./pio_decomp_tests --pio-tf-num-io-tasks=2 --pio-tf-stride=2)
#endif()
#ADD_TEST(pio_decomp_tests_4p_2iop_2str
#         ${MPIEXEC} ${MPIEXEC_NPF4}
#         ${test_exe})

#if (${PLATFORM} STREQUAL "yellowstone" )
#   set (test_exe "./pio_decomp_tests --pio-tf-num-io-tasks=2 --pio-tf-num-aggregators=1")
#else()
#   set (test_exe ./pio_decomp_tests --pio-tf-num-io-tasks=2 --pio-tf-num-aggregators=1)
#endif()
#ADD_TEST(pio_decomp_tests_4p_2iop_1agg
#         ${MPIEXEC} ${MPIEXEC_NPF4}
#         ${test_exe})
