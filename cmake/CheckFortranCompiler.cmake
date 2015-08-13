#
# - Functions to check/validate the Fortran compiler capabilities
#

#
# - Check for CSizeOf capability
#
function (check_Fortran_csizeof)
    
    if (NOT DEFINED Fortran_CSIZEOF)
        
        message (STATUS "Checking whether Fortran compiler has c_sizeof")
        try_compile (COMPILE_RESULT
                     ${CMAKE_CURRENT_BINARY_DIR}/tryFortran_csizeof
                     SOURCES ${CMAKE_SOURCE_DIR}/cmake/TryFortran_csizeof.f90
                     OUTPUT_VARIABLE TryFortran_OUT)
        if (COMPILE_RESULT)
            message (STATUS "Checking whether Fortran compiler has c_sizeof - yes")
        else ()
            message (STATUS "Checking whether Fortran compiler has c_sizeof - no")
        endif ()
        
        set (Fortran_CSIZEOF ${COMPILE_RESULT}
             CACHE BOOL "Whether the Fortran compiler has c_sizeof")
     
     endif ()
            
endfunction ()
