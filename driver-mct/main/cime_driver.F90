program cime_driver

!-------------------------------------------------------------------------------
!
! Purpose: Main program for a CIME-driven model.  Can have different
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
   use perf_mod,      only : t_startf, t_adj_detailf, t_stopf
   use ESMF,          only : ESMF_Initialize, ESMF_Finalize
   use seq_comm_mct,  only : esmf_logfile_kind
   use cime_comp_mod, only : cime_pre_init1
   use cime_comp_mod, only : cime_pre_init2
   use cime_comp_mod, only : cime_init
   use cime_comp_mod, only : cime_run
   use cime_comp_mod, only : cime_final

   implicit none

   !--------------------------------------------------------------------------
   ! Setup and initialize the communications and logging.
   !--------------------------------------------------------------------------
   call cime_pre_init1()

   !--------------------------------------------------------------------------
   ! Initialize ESMF.  This is done outside of the ESMF_INTERFACE ifdef
   ! because it is needed for the time manager, even if the ESMF_INTERFACE
   ! is not used.
   !--------------------------------------------------------------------------
   call ESMF_Initialize(logkindflag=esmf_logfile_kind)

   !--------------------------------------------------------------------------
   ! Read in the configuration information and initialize the time manager.
   !--------------------------------------------------------------------------
   ! Timer initialization has to be after determination of the maximum number
   ! of threads used across all components, so called inside of
   ! cime_pre_init2, as are t_startf and t_stopf for CPL:INIT and
   ! cime_pre_init2.
   !--------------------------------------------------------------------------
   call cime_pre_init2()

   !--------------------------------------------------------------------------
   ! Call the initialize, run and finalize routines.
   !--------------------------------------------------------------------------

   call t_startf('CPL:INIT')
   call t_adj_detailf(+1)

   call cime_init()

   call t_adj_detailf(-1)
   call t_stopf('CPL:INIT')

   call cime_run()
   call cime_final()

   !--------------------------------------------------------------------------
   ! Clean-up
   !--------------------------------------------------------------------------
   call ESMF_Finalize( )

end program cime_driver
