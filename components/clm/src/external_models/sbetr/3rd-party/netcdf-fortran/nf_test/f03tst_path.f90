program ftst_path

! Tests new nf_inq_path function
! Mimics tests in C tst_files5.c code

  USE netcdf_f03

  implicit NONE

  character(len=*), parameter  :: FILE_NAME="f03tst_path.nc"

  integer                      :: path_len, ncid
  character(LEN=NF_MAX_NAME+1) :: path_in

  path_in   = REPEAT(" ", LEN(path_in))
  path_len  = 0

  print *,''
  print *,'*** Testing netcdf file functions.'
  print *,'*** Checking the new inq_path function'

! Test with classic mode nf90_create

  call check(nf_create(FILE_NAME, nf_classic_model, ncid))
  call check(nf_inq_path(ncid, path_len, path_in))

  if ((path_len /= LEN(FILE_NAME)) .OR. (FILE_NAME /= TRIM(path_in)))  &
    call check(-1)
  call check(nf_close(ncid))

  path_in=REPEAT(" ", LEN(path_in))
  path_len=0

! Test with classic mode nf_open

  call check(nf_open(FILE_NAME, nf_classic_model, ncid))
  call check(nf_inq_path(ncid, path_len, path_in))

  if ((path_len /= LEN(FILE_NAME)) .OR. (FILE_NAME /= TRIM(path_in)))  &
    call check(-1)
  call check(nf_close(ncid))

  path_in=REPEAT(" ", LEN(path_in))
  path_len=0


! Test with netcdf4 mode nf_create

  call check(nf_create(FILE_NAME, nf_netcdf4, ncid))
  call check(nf_inq_path(ncid, path_len, path_in))

  if ((path_len /= LEN(FILE_NAME)) .OR. (FILE_NAME /= TRIM(path_in)))  &
    call check(-1)
  call check(nf_close(ncid))

  path_in=REPEAT(" ", LEN(path_in))
  path_len=0

! Test with netcdf4 mode nf_open

  call check(nf_open(FILE_NAME, nf_netcdf4, ncid))
  call check(nf_inq_path(ncid, path_len, path_in))

  if ((path_len /= LEN(FILE_NAME)) .OR. (FILE_NAME /= TRIM(path_in)))  &
    call check(-1)
  call check(nf_close(ncid))

  path_in=REPEAT(" ", LEN(path_in))
  path_len=0

  Print *,'*** SUCCESS!'

contains
!     This subroutine handles errors by printing an error message and
!     exiting with a non-zero status.
  subroutine check(errcode)
    use netcdf
    implicit none
    integer, intent(in) :: errcode

    if(errcode /= nf_noerr) then
       print *, 'Error: ', trim(nf_strerror(errcode))
       stop 2
    endif
  end subroutine check

end program ftst_path
