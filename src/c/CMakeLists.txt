cmake_minimum_required (VERSION 3.0.0)
project (PIOC C)

set (PIO_C_SRCS topology.c
               pio_nc.c
               pio_file.c
               pioc_support.c
               pio_lists.c
               pioc.c
               pioc_sc.c
               pio_spmd.c
               pio_rearrange.c
               pio_darray.c
               pio_put_nc.c 
               pio_get_nc.c
               bget.c)

add_library (pioc ${PIO_C_SRCS})

target_include_directories (pioc
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Compiler
target_compile_definitions (pioc
    PUBLIC CPRGNU
    PUBLIC LINUX)
target_compile_options (pioc
    PUBLIC -std=c99)
    
# GPTL Timing
if (ENABLE_TIMING)
    target_link_libraries (pioc gptl)
endif ()

# MPI
find_package (MPI)
target_include_directories (pioc 
    PUBLIC ${MPI_C_INCLUDE_PATH})
target_compile_definitions (pioc 
    PUBLIC USEMPIIO)
target_link_libraries (pioc
    ${MPI_C_LIBRARIES})

# NetCDF
target_include_directories (pioc 
    PUBLIC $ENV{HOME}/Software/installs/netcdf-c/include)
target_compile_definitions (pioc 
    PUBLIC _NETCDF
    PUBLIC _NETCDF4)
target_link_libraries (pioc
    $ENV{HOME}/Software/installs/netcdf-c/lib/libnetcdf.a)

# PnetCDF
target_include_directories (pioc 
    PUBLIC $ENV{HOME}/Software/installs/pnetcdf/include)
target_compile_definitions (pioc 
    PUBLIC _PNETCDF
    PUBLIC USE_PNETCDF_VARN
    PUBLIC USE_PNETCDF_VARN_ON_READ)
target_link_libraries (pioc
    $ENV{HOME}/Software/installs/pnetcdf/lib/libpnetcdf.a)
