add_library (pioc
    topology.c
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

target_include_directories (pioc
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_include_directories (pioc 
    PRIVATE $ENV{HOME}/Software/installs/netcdf-c/include)
target_compile_definitions (pioc 
    PRIVATE _NETCDF)
target_link_libraries (pioc
    $ENV{HOME}/Software/installs/netcdf-c/lib/libnetcdf.a)

target_include_directories (pioc 
    PRIVATE $ENV{HOME}/Software/installs/pnetcdf/include)
target_compile_definitions (pioc 
    PRIVATE _PNETCDF)
target_link_libraries (pioc
    $ENV{HOME}/Software/installs/pnetcdf/lib/libpnetcdf.a)
