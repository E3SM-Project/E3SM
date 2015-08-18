include (CMakeParseArguments)

#
# - Functions for parallel testing with CTest
#

#==============================================================================
# - Check the platform
#
# Defines the following variables:
#
#    PLATFORM - The name of the platform (alcf, ucar, nersc, etc)
#    PLATFORM_MPIEXEC - The MPIEXEC command name to use when launching tests
function (check_platform)

    if (NOT DEFINED PLATFORM)
        site_name (sitename)

        set (PLATFORM_HOSTNAME ${sitename}
             CACHE STRING "Platform host name")             

        # UCAR/NCAR Machines
        if (sitename MATCHES "^yslogin" OR
            sitename MATCHES "^geyser" OR
            sitename MATCHES "^caldera")

            set (PLATFORM "ucar"
                 CACHE STRING "Platform name")             
    
            set (PLATFORM_MPIEXEC "execca mpirun.lsf"
                 CACHE STRING "Platform MPI job launch command")

        # Argonne ALCF Machines
        elseif (sitename MATCHES "^mira" OR
                sitename MATCHES "^cetus" OR
                sitename MATCHES "^vesta" OR
                sitename MATCHES "^cooley")

            set (PLATFORM "alcf"
                 CACHE STRING "Platform name")             
    
            set (PLATFORM_MPIEXEC "qsub"
                 CACHE STRING "Platform MPI job launch command")

        # NERSC Machines
        elseif (sitename MATCHES "^edison" OR
                sitename MATCHES "^hopper" OR
                sitename MATCHES "^carver")

            set (PLATFORM "nersc"
                 CACHE STRING "Platform name")             
    
            set (PLATFORM_MPIEXEC ${MPIEXEC}
                 CACHE STRING "Platform MPI job launch command")

        # All other machines (depend upon FindMPI's MPIEXEC variable)
        else ()

            set (PLATFORM "unknown"
                 CACHE STRING "Platform name")             
    
            set (PLATFORM_MPIEXEC ${MPIEXEC}
                 CACHE STRING "Platform MPI job launch command")

        endif ()

    endif ()

endfunction ()

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
        set (EXE_CMD ${PLATFORM_MPIEXEC} ${exe_cmds} -n${num_procs})
        
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
        set(EXE_CMD ${PLATFORM_MPIEXEC} ${MPIEXEC_NPF} ${MPIEXEC_PREFLAGS} ${exe_cmds})
    endif()
    
    # Add the test to CTest
    add_test(NAME ${TESTNAME} COMMAND ${EXE_CMD})
    
    # Adjust the test timeout
    set_tests_properties(${TESTNAME} PROPERTIES TIMEOUT ${timeout})

endfunction()
