#==============================================================================
#  DEFINE THE TARGETS AND TESTS
#==============================================================================

set (SRC pioperformance.F90)

add_executable (pioperf ${SRC})

target_link_libraries (pioperf piof)

if (PIO_HDF5_LOGGING) 
    target_compile_definitions (pioperf 
        PUBLIC LOGGING)
endif ()
