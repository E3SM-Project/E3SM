&piotest_nml
ltest_bin        = .true.
!ltest_bin_direct = .true. ! Ignored if PIO is not built with -DUSEMPIIO
ltest_netcdf     = .true. ! Ignored if PIO is not built with -D_NETCDF
ltest_netcdf4p   = .true. ! Ignored if PIO is not built with -D_NETCDF4
ltest_netcdf4c   = .true. ! Ignored if PIO is not built with -D_NETCDF4
ltest_pnetcdf    = .true. ! Ignored if PIO is not built with -D_PNETCDF
stride           = 3
/
