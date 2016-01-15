include (CMakeParseArguments)

#
# - Functions for parallel testing with CTest
#

#==============================================================================
# - Get the machine platform-specific 
#
# Syntax:  platform_name (RETURN_VARIABLE)
#
function (platform_name RETURN_VARIABLE)

    # Determine platform name from site name...
    site_name (SITENAME)

    # UCAR/NWSC Machines
    if (SITENAME MATCHES "^yslogin" OR
        SITENAME MATCHES "^geyser" OR
        SITENAME MATCHES "^caldera")
        
        set (${RETURN_VARIABLE} "nwsc" PARENT_SCOPE)
        
    # ALCF/Argonne Machines
    elseif (SITENAME MATCHES "^mira" OR
            SITENAME MATCHES "^cetus" OR
            SITENAME MATCHES "^vesta" OR
            SITENAME MATCHES "^cooley")
        
        set (${RETURN_VARIABLE} "alcf" PARENT_SCOPE)
        
    # ALCF/Argonne Machines
    elseif (SITENAME MATCHES "^edison" OR
        SITENAME MATCHES "^cori")
        
        set (${RETURN_VARIABLE} "nersc" PARENT_SCOPE)

    # NCSA Machine (Blue Waters)
    elseif (SITENAME MATCHES "^h2ologin")

        set (${RETURN_VARIABLE} "ncsa" PARENT_SCOPE)
        
    else ()

        set (${RETURN_VARIABLE} "unknown" PARENT_SCOPE)
    
    endif ()

endfunction ()

#==============================================================================
# - Add a new parallel test
#
# Syntax:  add_mpi_test (<TESTNAME>
#                        EXECUTABLE <command>
#                        ARGUMENTS <arg1> <arg2> ...
#                        NUMPROCS <num_procs>
#                        TIMEOUT <timeout>)
function (add_mpi_test TESTNAME)

    # Parse the input arguments
    set (options)
    set (oneValueArgs NUMPROCS TIMEOUT EXECUTABLE)
    set (multiValueArgs ARGUMENTS)
    cmake_parse_arguments (${TESTNAME} "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    
    # Store parsed arguments for convenience
    set (exec_file ${${TESTNAME}_EXECUTABLE})
    set (exec_args ${${TESTNAME}_ARGUMENTS})
    set (num_procs ${${TESTNAME}_NUMPROCS})
    set (timeout ${${TESTNAME}_TIMEOUT})
    
    # Get the platform name
    platform_name (PLATFORM)
    
    # Default ("unknown" platform) execution
    if (PLATFORM STREQUAL "unknown")

        # Run tests directly from the command line
        set(EXE_CMD ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${num_procs} 
                    ${MPIEXEC_PREFLAGS} ${exec_file} 
                    ${MPIEXEC_POSTFLAGS} ${exec_args})

    else ()
                        
        # Run tests from the platform-specific executable
        set (EXE_CMD ${CMAKE_SOURCE_DIR}/cmake/mpiexec.${PLATFORM} 
                     ${num_procs} ${exec_file} ${exec_args})
                     
    endif ()
    
    # Add the test to CTest
    add_test(NAME ${TESTNAME} COMMAND ${EXE_CMD})
    
    # Adjust the test timeout
    set_tests_properties(${TESTNAME} PROPERTIES TIMEOUT ${timeout})

endfunction()
