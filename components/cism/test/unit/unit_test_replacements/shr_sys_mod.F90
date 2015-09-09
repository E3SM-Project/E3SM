! Trimmed-down version of shr_sys_mod including just what is needed for cism unit tests,
! in order to avoid dependencies

module shr_sys_mod

   use shr_kind_mod
   use shr_log_mod, only: s_loglev  => shr_log_Level
   use shr_log_mod, only: s_logunit => shr_log_Unit

   implicit none
   
   private

   public :: shr_sys_abort   ! abort a program
   public :: shr_sys_flush   ! flush an i/o buffer


contains
   
!===============================================================================
!===============================================================================

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

   call shr_sys_flush(s_logunit)
   if (len_trim(string) > 0) write(s_logunit,F00) 'ERROR: '//trim(string)
   write(s_logunit,F00) 'WARNING: calling shr_mpi_abort() and stopping'
   call shr_sys_flush(s_logunit)

! WJS (12-6-11): Removed some mpi-related stuff that is here in the real version of this
! subroutine

   call abort()
   stop

END SUBROUTINE shr_sys_abort

!===============================================================================
!===============================================================================

SUBROUTINE shr_sys_flush(unit)

   IMPLICIT none

   !----- arguments -----
   integer(SHR_KIND_IN) :: unit  ! flush output buffer for this unit

   !----- formats -----
   character(*),parameter :: subName =   '(shr_sys_flush) '
   character(*),parameter :: F00     = "('(shr_sys_flush) ',4a)"

!-------------------------------------------------------------------------------
! PURPOSE: an architecture independant system call
!-------------------------------------------------------------------------------

! WJS (12-6-11): I have reworked this from the real version, in order to allow
! reassonable behavior when the sysstem is not defined

#if (defined AIX)
   call flush_(unit)
#else
   call flush(unit)
#endif

END SUBROUTINE shr_sys_flush

!===============================================================================
!===============================================================================

end module shr_sys_mod
