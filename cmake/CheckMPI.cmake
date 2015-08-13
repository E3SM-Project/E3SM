#
# - Functions to check/validate the MPI package that was found
#

#
# - Check for MPI I/O capability
#
function (check_MPIIO SEARCH_DIRS)
    
    if (NOT DEFINED MPIIO_DETECTED)
        
        message (STATUS "Checking whether MPIIO is supported")
        find_file (TryMPIIO_FILE
                   NAMES TryMPIIO.f90
                   HINTS ${SEARCH_DIRS})
        if (TryMPIIO_FILE)
            try_compile (COMPILE_RESULT
                         ${CMAKE_CURRENT_BINARY_DIR}/tryMPIIO
                         SOURCES ${TryMPIIO_FILE}
                         OUTPUT_VARIABLE TryFortran_OUT)
            if (COMPILE_RESULT)
                message (STATUS "Checking whether MPIIO is supported - yes")
            else ()
                message (STATUS "Checking whether MPIIO is supported - no")
            endif ()
            
            set (MPIIO_DETECTED ${COMPILE_RESULT}
                 CACHE BOOL "Whether MPIIO is supported")
                 
        else ()
            message (STATUS "Checking whether MPIIO is supported - failed")
        endif ()
     
    endif ()
            
endfunction ()

#
# - Check for MPI Fortran Module
#
function (check_MPIMOD SEARCH_DIRS)
    
    if (NOT DEFINED MPIMOD_DETECTED)
        
        message (STATUS "Checking whether MPI Fortran module is supported")
        find_file (TryMPIMOD_FILE
                   NAMES TryMPIMOD.f90
                   HINTS ${SEARCH_DIRS})
        if (TryMPIMOD_FILE)
            try_compile (COMPILE_RESULT
                         ${CMAKE_CURRENT_BINARY_DIR}/tryMPIMOD
                         SOURCES ${TryMPIMOD_FILE}
                         OUTPUT_VARIABLE TryFortran_OUT)
            if (COMPILE_RESULT)
                message (STATUS "Checking whether MPI Fortran module is supported - yes")
            else ()
                message (STATUS "Checking whether MPI Fortran module is supported - no")
            endif ()
            
            set (MPIMOD_DETECTED ${COMPILE_RESULT}
                 CACHE BOOL "Whether MPI Fortran module is supported")
                 
        else ()
            message (STATUS "Checking whether MPI Fortran module is supported - failed")
        endif ()
     
    endif ()
            
endfunction ()