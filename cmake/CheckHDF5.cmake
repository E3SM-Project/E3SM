#
# - Functions to check/validate the HDF5 package(s) that were found
#

#
# - Check HDF5_C
#
function (check_HDF5_C)

    if (NOT HDF5_C_CHECKED)

        set (HDF5_C_CHECKED TRUE CACHE BOOL "HDF5_C checked")
        
    endif ()
        
endfunction ()
