program TryNC4
  use netcdf
  integer :: ierr, fh, varid
  ierr = nf90_var_par_access(fh,varid,NF90_COLLECTIVE)
end program
