program cesm_driver

!-------------------------------------------------------------------------------
!
! Purpose: Main program for CESM and ACME models.  Can have different
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
   use perf_mod
   use ESMF
   use cesm_comp_mod, only : cesm_pre_init1
   use cesm_comp_mod, only : cesm_pre_init2
   use cesm_comp_mod, only : cesm_init
   use cesm_comp_mod, only : cesm_run
   use cesm_comp_mod, only : cesm_final

   implicit none

   !--------------------------------------------------------------------------
   ! Setup and initialize the communications and logging.  
   !--------------------------------------------------------------------------
   call cesm_pre_init1()

   !--------------------------------------------------------------------------
   ! Initialize ESMF.  This is done outside of the ESMF_INTERFACE ifdef
   ! because it is needed for the time manager, even if the ESMF_INTERFACE
   ! is not used.
   !--------------------------------------------------------------------------
   call ESMF_Initialize()

   !--------------------------------------------------------------------------
   ! Read in the configuration information and initialize the time manager.
   !--------------------------------------------------------------------------
   call cesm_pre_init2()

   !--------------------------------------------------------------------------
   ! Call the initialize, run and finalize routines.
   !--------------------------------------------------------------------------
   call cesm_init()
   call cesm_run()
   call cesm_final()

   !--------------------------------------------------------------------------
   ! Clean-up
   !--------------------------------------------------------------------------
   call ESMF_Finalize( )

end program cesm_driver
