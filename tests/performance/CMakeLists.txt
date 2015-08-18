#==============================================================================
#  DEFINE THE TARGETS AND TESTS
#==============================================================================

set (SRC pioperformance.F90)

add_executable (pio_perf_test ${SRC})

target_link_libraries (pio_perf_test piof pioc gptl)

include_directories(${PIO_Fortran_INCLUDE_DIRS})
add_definitions(${PIO_Fortran_DEFINITIONS})
list (APPEND CMAKE_Fortran_FLAGS ${PIO_Fortran_COMPILE_FLAGS})

if (PIO_HDF5_LOGGING) 
    add_definitions (-DLOGGING)
endif ()
