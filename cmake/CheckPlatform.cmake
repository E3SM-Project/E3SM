#
# - Functions to check if the system/machine is recognized platform
#   This should be called AFTER find_package(MPI) is called.

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

            set (PLATFORM_MPIEXEC_NPF "-n"
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

            set (PLATFORM_MPIEXEC_NPF "-n"
                 CACHE STRING "Platform MPI job launch command")

        # NERSC Machines
        elseif (sitename MATCHES "^edison" OR
                sitename MATCHES "^hopper" OR
                sitename MATCHES "^carver")

            set (PLATFORM "nersc"
                 CACHE STRING "Platform name")             
    
            set (PLATFORM_MPIEXEC ""
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