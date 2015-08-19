#==============================================================================
#  DEFINE THE TARGETS AND TESTS
#==============================================================================

set (SRC pioperformance.F90)

add_executable (pioperf ${SRC})

target_link_libraries (pioperf piof)

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    target_compile_options (pioperf
        PRIVATE -ffree-line-length-none)
endif()

if (PIO_HDF5_LOGGING) 
    target_compile_definitions (pioperf 
        PUBLIC LOGGING)
endif ()
