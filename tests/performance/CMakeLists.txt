#==============================================================================
#  DEFINE THE TARGETS AND TESTS
#==============================================================================

set (SRC pioperformance.F90)

add_executable (pio_perf_test ${SRC})

target_link_libraries (pio_perf_test piof pioc gptl)

if (PIO_HDF5_LOGGING) 
    add_definitions (-DLOGGING)
endif ()

set_target_properties (pio_perf_test PROPERTIES
    INCLUDE_DIRECTORIES "${PIO_Fortran_INCLUDE_DIRS}"
    COMPILE_DEFINITIONS "${PIO_Fortran_DEFINITIONS}"
    COMPILE_FLAGS "${PIO_Fortran_COMPILE_FLAGS}")
