!=======================================================================
!BOP
!
! !MODULE: CICE_ComponentMod - ESMF component module for CICE
!
! !DESCRIPTION:
!
!  This module contains the routines for making CICE a component,
!  particularly for the Earth System Modeling Framework (ESMF).
!  It primarily makes public a routine for registering the init,
!  finalize and run methods for CICE.
!
! !REVISION HISTORY:
!  SVN:$Id: CICE_ComponentMod.F90 56 2007-03-15 14:42:35Z dbailey $
!
! authors:  Philip W. Jones, LANL
!
! 2005: Introduced module for ESMF compliance
! 2006: Converted to free source form (F90) by Elizabeth Hunke
!
! !INTERFACE:
!
      module CICE_ComponentMod
!
! !USES:
!
      use ice_kinds_mod
#ifdef USE_ESMF
      use esmf
      use CICE_InitMod
      use CICE_RunMod
      use CICE_FinalMod
#endif

      implicit none
      private
      save

! !PUBLIC MEMBER FUNCTIONS:

#ifdef USE_ESMF
      public :: CICE_SetServices
#endif
!
!EOP
!

#ifdef USE_ESMF    
! entire module

!=======================================================================

      contains

!=======================================================================

!BOP
!
! !IROUTINE: CICE_SetServices - registers methods with ESMF
!
! !INTERFACE:
!
      subroutine CICE_SetServices(griddedComp, errorCode)

!
! !DESCRIPTION:
!
!     This routine registers the initialize, run and finalize methods
!     for the CICE model in the ESMF.
!
! !REVISION HISTORY:
!
!     authors:  same as module
!

!
! !INPUT/OUTPUT PARAMETERS:

      type (ESMF_GridComp) :: &
          griddedComp            ! ESMF gridded component to which
                                 !  services are assigned

      integer (int_kind) :: &
          errorCode              ! Returns an error code if any init fails

!
!EOP
!BOC
!

   !--------------------------------------------------------------------
   !  initialize return flag
   !--------------------------------------------------------------------

      errorCode = ESMF_Success

   !--------------------------------------------------------------------
   !  register init method
   !--------------------------------------------------------------------

      call ESMF_GridCompSetEntryPoint(griddedComp, ESMF_SETINIT, &
                                      CICE_Initialize, 0, errorCode)

      if (errorCode /= ESMF_Success) then
         ! add error message here or use ESMF error and logging facilities
         return
      endif

   !--------------------------------------------------------------------
   !  register run method
   !--------------------------------------------------------------------

      call ESMF_GridCompSetEntryPoint(griddedComp, ESMF_SETRUN, &
                                      CICE_Run, 0, errorCode)

      if (errorCode /= ESMF_Success) then
         ! add error message here or use ESMF error and logging facilities
         return
      endif

   !--------------------------------------------------------------------
   !  register finalize method
   !--------------------------------------------------------------------

      call ESMF_GridCompSetEntryPoint(griddedComp, ESMF_SETFINAL, &
                                      CICE_Finalize, 0, errorCode)

      if (errorCode /= ESMF_Success) then
         ! add error message here or use ESMF error and logging facilities
         return
      endif

!
!EOC
!

      end subroutine CICE_SetServices

#endif   
!ifdef USE_ESMF

!=======================================================================

      end module CICE_ComponentMod

!=======================================================================
