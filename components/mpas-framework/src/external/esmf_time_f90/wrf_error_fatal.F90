
subroutine wrf_error_fatal(msg)
! use shr_sys_mod, only: shr_sys_abort
  implicit none
  character(len=*), intent(in) :: msg
  write(6,*) 'wrf_error_fatal: ',trim(msg)
! call shr_sys_abort( msg )
end subroutine wrf_error_fatal

