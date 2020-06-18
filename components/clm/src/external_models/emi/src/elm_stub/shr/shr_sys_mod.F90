MODULE shr_sys_mod

   use shr_kind_mod  ! defines real & integer kinds
   use shr_log_mod, only: s_logunit => shr_log_Unit
   implicit none
! PUBLIC: Public interfaces

   private

   public :: shr_sys_abort   ! abort a program
   public :: shr_sys_flush   !flush
   public :: shr_sys_chdir   ! change current working dir
   public :: shr_sys_sleep   ! have program sleep for a while
   public :: shr_sys_system  ! make a system call
   contains

SUBROUTINE shr_sys_abort(string,rc)

   IMPLICIT none

   character(*)        ,optional :: string  ! error message string
   integer(SHR_KIND_IN),optional :: rc      ! error code

   !----- local -----
   integer(SHR_KIND_IN) :: ierr
   logical              :: flag

   !----- formats -----
   character(*),parameter :: subName =   '(shr_sys_abort) '
   character(*),parameter :: F00     = "('(shr_sys_abort) ',4a)"

!-------------------------------------------------------------------------------
! PURPOSE: consistent stopping mechanism
!-------------------------------------------------------------------------------

   if (present(string)) then
      if (len_trim(string) > 0) write(s_logunit,F00) 'ERROR: '//trim(string)
   end if
   write(s_logunit,F00) 'WARNING: calling shr_mpi_abort() and stopping'

   stop
   rc = 0

END SUBROUTINE shr_sys_abort

SUBROUTINE shr_sys_flush(iulog)
  implicit none
  integer, intent(in) :: iulog

  flush(iulog)

END SUBROUTINE shr_sys_flush



SUBROUTINE shr_sys_chdir(path, rcode)

   IMPLICIT none

   !----- arguments -----
   character(*)        ,intent(in)  :: path    ! chdir to this dir
   integer(SHR_KIND_IN),intent(out) :: rcode   ! return code

   !----- local -----
   integer(SHR_KIND_IN)             :: lenpath ! length of path
#if (defined AIX || (defined LINUX && !defined CPRGNU && !defined CPRNAG) || defined CPRINTEL)
   integer(SHR_KIND_IN),external    :: chdir   ! AIX system call
#endif

   !----- formats -----
   character(*),parameter :: subName =   '(shr_sys_chdir) '
   character(*),parameter :: F00     = "('(shr_sys_chdir) ',4a)"
   rcode = 0
   if (len(path) > 0) continue
end SUBROUTINE shr_sys_chdir



SUBROUTINE shr_sys_sleep(sec)

   IMPLICIT none

   !----- arguments -----
   real   (SHR_KIND_R8),intent(in) :: sec  ! number of seconds to sleep

   !----- local -----
   integer(SHR_KIND_IN) :: isec   ! integer number of seconds
   integer(SHR_KIND_IN) :: rcode  ! return code
   character(90)        :: str    ! system call string

   !----- formats -----
   character(*),parameter :: subName =   '(shr_sys_sleep) '
   character(*),parameter :: F00     = "('(shr_sys_sleep) ',4a)"
   character(*),parameter :: F10     = "('sleep ',i8 )"
   if (sec > 0.0) continue
end SUBROUTINE shr_sys_sleep


SUBROUTINE shr_sys_system(str,rcode)

   IMPLICIT none

   !----- arguments ---
   character(*)        ,intent(in)  :: str    ! system/shell command string
   integer(SHR_KIND_IN),intent(out) :: rcode  ! function return error code

   !----- functions -----
#if (defined LINUX && !defined CPRGNU)
   integer(SHR_KIND_IN),external    :: system ! function to envoke shell command
#endif

   !----- formats -----
   character(*),parameter :: subName =   '(shr_sys_system) '
   character(*),parameter :: F00     = "('(shr_sys_system) ',4a)"
   rcode = 0
   if (len(str) == 0) continue
end SUBROUTINE shr_sys_system   
end module shr_sys_mod
