!
! NetCDF Fortran Test for version string
!

program nf90_version
    use netcdf
    implicit none
    print *, trim(nf90_inq_libvers())
end program nf90_version