
SUBROUTINE wrf_message( str )
  IMPLICIT NONE
  CHARACTER*(*) str
  write(6,*) 'wrf_message ',trim(str)
END SUBROUTINE wrf_message


SUBROUTINE wrf_error_fatal( str )
  IMPLICIT NONE
  CHARACTER*(*) str
  write(6,*) 'wrf_error_fatal ',trim(str)
  stop
END SUBROUTINE wrf_error_fatal



