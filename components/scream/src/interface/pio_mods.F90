module pio_mods

    use piolib_mod, only : PIO_init, PIO_finalize, PIO_createfile, PIO_closefile, &
        PIO_initdecomp, PIO_freedecomp, PIO_syncfile, PIO_openfile, PIO_setframe
    use pio_types,  only : iosystem_desc_t, file_desc_t, &
        pio_noerr, PIO_iotype_netcdf, var_desc_t, io_desc_t, PIO_int, &
        pio_clobber, PIO_nowrite, PIO_unlimited, pio_global
    use pio_kinds,  only : PIO_OFFSET_KIND, i4 
    use pio_nf,     only : PIO_redef, PIO_def_dim, PIO_def_var, PIO_enddef
    use piodarray,  only : PIO_write_darray, PIO_read_darray 
    use pionfatt_mod, only : PIO_put_att   => put_att

end module pio_mods
