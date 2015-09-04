subroutine handle_error (ret)
  implicit none

  integer ret

  include 'netcdf.inc'

  write(6,*) nf_strerror (ret)
  call abort
  stop 999
end subroutine handle_error
