#==============================================================================
#  DEFINE THE TARGETS AND TESTS
#==============================================================================

set (SRC pioperformance.F90)

add_executable (pio_perf_test ${SRC})
target_link_libraries (pio_perf_test piof pioc gptl)

if (PIO_HDF5_LOGGING) 
    target_compile_definitions (pio_perf_test 
        PRIVATE LOGGING)
endif ()
