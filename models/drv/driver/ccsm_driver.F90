program ccsm_driver

!-------------------------------------------------------------------------------
!
! Purpose: Main program for NCAR CCSM4/cpl7. Can have different
!          land, sea-ice, and ocean models plugged in at compile-time.
!          These models can be either: stub, dead, data, or active
!          components or some combination of the above.
!
!               stub -------- Do nothing.
!               dead -------- Send analytic data back.
!               data -------- Send data back interpolated from input files.
!               active ------ Prognostically simulate the given component.
!
! Method: Call appropriate initialization, run (time-stepping), and 
!         finalization routines.
! 
!-------------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! share code & libs
   !----------------------------------------------------------------------------
   use shr_sys_mod,       only: shr_sys_abort
   use perf_mod
   use ESMF
   use ccsm_comp_mod

   implicit none

   !--------------------------------------------------------------------------
   ! Local Variables
   !--------------------------------------------------------------------------
   integer                    :: localrc
#ifdef USE_ESMF_LIB
   character(len=ESMF_MAXSTR) :: compName
   type(ESMF_CplComp)         :: drvcomp
#endif


   !--------------------------------------------------------------------------
   ! Setup and initialize the communications and logging.  
   !--------------------------------------------------------------------------
   call ccsm_pre_init1()

   !--------------------------------------------------------------------------
   ! Initialize ESMF.  This is done outside of the ESMF_INTERFACE ifdef
   ! because it is needed for the time manager, even if the ESMF_INTERFACE
   ! is not used.
   !--------------------------------------------------------------------------
   call ESMF_Initialize()

   !--------------------------------------------------------------------------
   ! Read in the configuration information and initialize the time manager.
   !--------------------------------------------------------------------------
   call ccsm_pre_init2()

#ifdef USE_ESMF_LIB

   !--------------------------------------------------------------------------
   ! Create the "Cap" component and set the services using the register
   ! routine.  Setting the services is where the initialize, run and 
   ! finalize routines for this component are set.
   !--------------------------------------------------------------------------
   compName = "CESM_Component"
   drvcomp = ESMF_CplCompCreate(name=compName, rc=localrc)
   if (localrc /= 0) call shr_sys_abort('failed to create CESM Component')

   call ESMF_CplCompSetServices(drvcomp, userRoutine=ccsm_comp_register, &
        rc=localrc)
   if (localrc /= 0) call shr_sys_abort('failed to set services for CESM Comp')

   !--------------------------------------------------------------------------
   ! Call the initialize, run and finalize routines registered with the
   ! cap component.
   !--------------------------------------------------------------------------

   call ESMF_CplCompInitialize(drvcomp, rc=localrc)
   if (localrc /= 0) call shr_sys_abort('failed to esmf initialize')

   call ESMF_CplCompRun(drvcomp, rc=localrc)
   if (localrc /= 0) call shr_sys_abort('failed to esmf run')

   call ESMF_CplCompFinalize(drvcomp, rc=localrc)
   if (localrc /= 0) call shr_sys_abort('failed to esmf finalize')

#else

   !--------------------------------------------------------------------------
   ! If ESMF is not defined, then just call the initialize, run and finalize
   ! routines directly.
   !--------------------------------------------------------------------------
   call ccsm_init()
   call ccsm_run()
   call ccsm_final()

#endif

   !--------------------------------------------------------------------------
   ! Clean-up
   !--------------------------------------------------------------------------
   call ESMF_Finalize( )


end program ccsm_driver
