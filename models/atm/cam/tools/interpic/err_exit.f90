subroutine err_exit (string)
  implicit none

  character*(*) string

  write(6,*) string
  stop 999
end subroutine err_exit
