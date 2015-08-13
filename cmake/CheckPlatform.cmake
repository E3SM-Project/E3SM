#
# - Functions to check if the system/machine is recognized platform
#

function (check_platform)

    if (NOT DEFINED PLATFORM)
        site_name (sitename)

        # UCAR/NCAR Machines
        if (sitename MATCHES "^yslogin")
            set (PLATFORM "ucar" CACHE STRING "Yellowstone Platform")

        # Argonne ALCF Machines
        elseif (sitename MATCHES "^mira")
            set (PLATFORM "alcf" CACHE STRING "Mira Platform")
        elseif (sitename MATCHES "^cetus")
            set (PLATFORM "alcf" CACHE STRING "Cetus Platform")
        elseif (sitename MATCHES "^vesta")
            set (PLATFORM "alcf" CACHE STRING "Vesta Platform")
        elseif (sitename MATCHES "^cooley")
            set (PLATFORM "alcf" CACHE STRING "Cooley Platform")

        # NERSC Machines
        elseif (sitename MATCHES "^edison")
            set (PLATFORM "nersc" CACHE STRING "Edison Platform")
        elseif (sitename MATCHES "^hopper")
            set (PLATFORM "nersc" CACHE STRING "Hopper Platform")
        elseif (sitename MATCHES "^carver")
            set (PLATFORM "nersc" CACHE STRING "Carver Platform")

        endif ()

    endif ()

endfunction ()