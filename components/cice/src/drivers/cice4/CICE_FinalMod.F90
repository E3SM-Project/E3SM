!=======================================================================
!
!BOP
!
! !MODULE: CICE_FinalMod - routines for final exit of CICE model
!
! !DESCRIPTION:
!
!  This module contains routines for the final exit of the CICE model,
!  including final output and clean exit from any message passing
!  environments and frameworks.
!
! !REVISION HISTORY:
!  SVN:$Id: CICE_FinalMod.F90 136 2008-06-27 17:07:19Z eclare $
!
!  authors: Philip W. Jones, LANL
!  2006: Converted to free source form (F90) by Elizabeth Hunke
!  2008: E. Hunke moved ESMF code to its own driver
!
! !INTERFACE:
!
      module CICE_FinalMod
!
! !USES:
!
      use ice_exit
      use ice_fileunits
      use ice_kinds_mod
      use ice_timers

      implicit none
      private
      save

! !PUBLIC MEMBER FUNCTIONS:

      public :: CICE_Finalize

!
!EOP
!
!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: CICE_Finalize - final exit of CICE model
!
! !DESCRIPTION:
!
!  This routine shuts down CICE by exiting all relevent environments.
!
! !REVISION HISTORY:
!
!  author same as module
!
! !INTERFACE:
!
      subroutine CICE_Finalize
!
!EOP
!BOC
!
   !-------------------------------------------------------------------
   ! stop timers and print timer info
   !-------------------------------------------------------------------

      call ice_timer_stop(timer_total)        ! stop timing entire run
      call ice_timer_print_all(stats=.false.) ! print timing information

!echmod      if (nu_diag /= 6) close (nu_diag) ! diagnostic output
      call release_all_fileunits

#ifndef coupled
      call end_run       ! quit MPI
#endif
!
!EOC
!
      end subroutine CICE_Finalize

!=======================================================================

      end module CICE_FinalMod

!=======================================================================
