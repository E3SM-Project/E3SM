#
# - Functions to check/validate the Fortran compiler capabilities
#

#
# - Check for CSizeOf capability
#
function (check_Fortran_csizeof SEARCH_DIRS)
    
    if (NOT DEFINED Fortran_CSIZEOF)
        
        message (STATUS "Checking whether Fortran compiler has c_sizeof")
        find_file (TryCSIZEOF_FILE
                   NAMES TryCSizeOf.f90
                   HINTS ${SEARCH_DIRS})
        if (TryCSIZEOF_FILE)
            try_compile (COMPILE_RESULT
                         ${CMAKE_CURRENT_BINARY_DIR}/tryFortran_csizeof
                         SOURCES ${TryCSIZEOF_FILE}
                         OUTPUT_VARIABLE TryFortran_OUT)
            if (COMPILE_RESULT)
                message (STATUS "Checking whether Fortran compiler has c_sizeof - yes")
            else ()
                message (STATUS "Checking whether Fortran compiler has c_sizeof - no")
            endif ()
            
            set (Fortran_CSIZEOF ${COMPILE_RESULT}
                 CACHE BOOL "Whether the Fortran compiler has c_sizeof")
                 
         else ()
            message (STATUS "Checking whether Fortran compiler has c_sizeof - failed")
         endif ()
     
     endif ()
            
endfunction ()
