program cesm_driver

!-------------------------------------------------------------------------------
!
! Purpose: Main program for NCAR CESM  Can have different
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
   use cesm_comp_mod

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
   ! Timer initialization has to be after determination of the maximum number
   ! of threads used across all components, so called inside of 
   ! ccsm_pre_init2, as are t_startf and t_stopf for CPL:INIT and 
   ! cesm_pre_init2.
   !--------------------------------------------------------------------------
   call cesm_pre_init2()

   call t_startf('CPL:INIT')
   call t_adj_detailf(+1)
#ifdef USE_ESMF_LIB

   !--------------------------------------------------------------------------
   ! Create the "Cap" component and set the services using the register
   ! routine.  Setting the services is where the initialize, run and 
   ! finalize routines for this component are set.
   !--------------------------------------------------------------------------
   compName = "CESM_Component"
   drvcomp = ESMF_CplCompCreate(name=compName, rc=localrc)
   if (localrc /= 0) call shr_sys_abort('failed to create CESM Component')

   call ESMF_CplCompSetServices(drvcomp, userRoutine=cesm_comp_register, &
        rc=localrc)
   if (localrc /= 0) call shr_sys_abort('failed to set services for CESM Comp')

   !--------------------------------------------------------------------------
   ! Call the initialize, run and finalize routines registered with the
   ! cap component.
   !--------------------------------------------------------------------------

   call ESMF_CplCompInitialize(drvcomp, rc=localrc)
   if (localrc /= 0) call shr_sys_abort('failed to esmf initialize')

   call t_adj_detailf(-1)
   call t_stopf('CPL:INIT')

   call ESMF_CplCompRun(drvcomp, rc=localrc)
   if (localrc /= 0) call shr_sys_abort('failed to esmf run')

   call ESMF_CplCompFinalize(drvcomp, rc=localrc)
   if (localrc /= 0) call shr_sys_abort('failed to esmf finalize')

#else

   !--------------------------------------------------------------------------
   ! If ESMF is not defined, then just call the initialize, run and finalize
   ! routines directly.
   !--------------------------------------------------------------------------
   call cesm_init()

   call t_adj_detailf(-1)
   call t_stopf('CPL:INIT')

   call cesm_run()
   call cesm_final()

#endif

   !--------------------------------------------------------------------------
   ! Clean-up
   !--------------------------------------------------------------------------
   call ESMF_Finalize( )


end program cesm_driver
