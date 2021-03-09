subroutine addglobal (ncid, cmdline)
  use netcdf
  use error_messages, only : handle_ncerr  

  implicit none
  
!
! Input arguments
!
  integer ncid
  character*(*) cmdline
!
! Local workspace
!
  integer ret
  integer numchars
  integer values(8)
  integer hnum
  integer hlen

  character*8 date
  character*10 time
  character*5 zone
  character*18 datetime
  character*16 logname
  character*16 hostname
  character*1500 :: hist

  call date_and_time (date, time, zone, values)

  datetime(1:8) =        date(5:6) // '/' // date(7:8) // '/' // date(3:4)
  datetime(9:)  = ' ' // time(1:2) // ':' // time(3:4) // ':' // time(5:6) // ' '

  call getenv ('LOGNAME', logname)
  call getenv ('HOST', hostname)

  hlen = 0
  hist = ' '
  if (nf90_inquire_attribute (ncid, nf90_global, 'history', len=hlen, attnum=hnum) == nf90_noerr) then
    call handle_ncerr( nf90_get_att (ncid, nf90_global, 'history', hist),&
      'addglobal.F90:40')
  end if

  hist = trim (hist) // char(10) // datetime // trim (logname) // ':' // &
         trim (hostname) // ':' // trim (cmdline)
!
! Add 3 to account for 1st newline and colons between each of 2 trimmed strings
!
  hlen = hlen + len(datetime) + len_trim(logname) + len_trim(hostname) + &
                len_trim(cmdline) + 3

  if (hlen > len (hist)) then
    write(6,*)'Warning: history attribute too long: truncating'
    hlen = len (hist)
  end if

  numchars = len_trim (hist)
  ret = nf90_put_att (ncid, nf90_global, 'history',  hist)

  return
end subroutine addglobal
  
#ifdef UNICOSMP
subroutine getenv (name, val)
  implicit none

  character(len=*), intent(in) :: name
  character(len=*), intent(out) :: val

  integer :: lenname
  integer :: lenval
  integer :: ierr

  call pxfgetenv (name, lenname, val, lenval, ierr)
  if (ierr /= 0) then
    write(6,*)'getenv: ierr not 0:', ierr
    stop 999
  end if

  return
end subroutine getenv
#endif
