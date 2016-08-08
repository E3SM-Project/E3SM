program freeform
  use pnetcdf
  integer ierr
  ierr = nfmpi_put_att(4, 1, 'fred', 7)
end program
