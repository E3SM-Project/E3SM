!
! PnetCDF Fortran Test for version string
!

program pnf_version
    use pnetcdf
    implicit none
    print *, trim(nfmpi_inq_libvers())
end program pnf_version