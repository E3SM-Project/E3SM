#==============================================================================
#  DEFINE THE TARGETS AND TESTS
#==============================================================================

set (SRC pioperformance.F90)

add_executable (pio_perf_test ${SRC})
include_directories(${PIO_Fortran_INCLUDE_DIRS})
target_link_libraries (pio_perf_test piof pioc gptl)

if (PIO_HDF5_LOGGING) 
    add_definitions (-DLOGGING)
endif ()
