module abortutils

!-----------------------------------------------------------------------
!BOP
! !MODULE: abortutils
!
! !DESCRIPTION:
! Abort the model for abnormal termination
!
! !USES:
  use clm_varctl, only : iulog
!
! !REVISION HISTORY:
! Author: CCM Core group
!
!EOP
!-----------------------------------------------------------------------

   private
   save

   public :: endrun

CONTAINS

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: endrun
!
! !INTERFACE:
   subroutine endrun(msg)
!
! !DESCRIPTION:
! Abort the model for abnormal termination
!
! Author: CCM Core group
!
!-----------------------------------------------------------------------
   use spmdMod,     only: mpicom
   use shr_sys_mod, only: shr_sys_flush
!
! !ARGUMENTS:
   implicit none
   character(len=*), intent(in), optional :: msg    ! string to be printed
!
! !REVISION HISTORY:
! Author: CCM Core group
!
!EOP
!-----------------------------------------------------------------------

   if (present (msg)) then
      write(iulog,*)'ENDRUN:', msg
   else
      write(iulog,*)'ENDRUN: called without a message string'
   end if

#if defined(AIX) && !defined(BGL) && !defined(BGP)
   close(5)    ! needed to prevent batch jobs from hanging in xl__trbk
   call xl__trbk()
#endif

   call shr_sys_flush(iulog)   ! Flush all output to standard output

   ! passing an argument of 1 to mpi_abort will lead to a STOPALL output
   ! error code of 257
   call mpi_abort (mpicom, 1)

end subroutine endrun

end module abortutils
