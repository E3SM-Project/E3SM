include (CMakeParseArguments)
include (CheckPlatform)
check_platform ()

#
# - Functions for parallel testing with CTest
#

#==============================================================================
# - Add a new parallel test
#
# Syntax:  add_mpi_test (<TESTNAME>
#                        COMMAND <command> <arg1> <arg2> ...
#                        NUMPROCS <num_procs>
#                        TIMEOUT <timeout>)
function (add_mpi_test TESTNAME)

    # Parse the input arguments
    set (options)
    set (oneValueArgs NUMPROCS TIMEOUT)
    set (multiValueArgs COMMAND)
    cmake_parse_arguments (${TESTNAME} "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    
    # Store parsed arguments for convenience
    set (exe_cmds ${${TESTNAME}_COMMAND})
    set (num_procs ${${TESTNAME}_NUMPROCS})
    set (timeout ${${TESTNAME}_TIMEOUT})
    
    # UCAR LSF execution
    if (PLATFORM STREQUAL "ucar" )
        ###
        ### note: no space between -n and num_procs for mpirun.lsf on
        ### yellowstone
        ###
        set (EXE_CMD execca mpirun.lsf ${exe_cmds} -n${num_procs})
        
    # Argonne COBALT execution
    elseif (PLATFORM STREQUAL "alcf" )
        ###
        ###
        #set(PIO_RUNJOB ${CMAKE_BINARY_DIR}/scripts/pio_runjob.sh)
        set (REQUIRED_OPTION --block \$ENV{COBALT_PARTNAME}) 
        set (RUNJOB_NPF --np ${num_procs})
        if (DEFINED ENV{BGQ_RUNJOB})
          set (RUNJOB $ENV{BGQ_RUNJOB})
        else()
          set (RUNJOB runjob)
        endif()
        set (EXE_CMD ${RUNJOB} ${RUNJOB_NPF} ${REQUIRED_OPTION}
                     ${MPIEXEC_PREFLAGS} : ${exe_cmds})
    
    # All others
    else()
        set(MPIEXEC_NPF ${MPIEXEC_NUMPROC_FLAG} ${num_procs})
        set(EXE_CMD ${MPIEXEC} ${MPIEXEC_NPF} ${MPIEXEC_PREFLAGS} ${exe_cmds})
    endif()
    
    # Add the test to CTest
    add_test(NAME ${TESTNAME} COMMAND ${EXE_CMD})
    
    # Adjust the test timeout
    set_tests_properties(${TESTNAME} PROPERTIES TIMEOUT ${timeout})

endfunction()
