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
target_link_libraries (pioc pnetcdf)
target_link_libraries (pioc netcdf)
