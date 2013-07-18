

program checkNCLink
  use netcdf
  implicit none

  integer, parameter :: i4 = selected_int_kind(6)
  integer(i4) :: fileHandle
  character(32) :: fileName
  integer :: nmode
  integer :: ierr

  fileName = "someDummyFile.nc"

  nmode = ior(nmode,NF90_NETCDF4)
  ierr = nf90_create(fileName,nmode,fileHandle)

end program
