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
  use shr_kind_mod,  only : r8 => SHR_KIND_R8
  use shr_kind_mod,  only : i8 => SHR_KIND_I8
  use shr_sys_mod,   only : shr_sys_irtc
  use perf_mod,      only : t_startf, t_adj_detailf, t_stopf, t_startstop_valsf
  use ESMF,          only : ESMF_Initialize, ESMF_Finalize
  use seq_comm_mct,  only : esmf_logfile_kind
  use cime_comp_mod, only : cime_pre_init1
  use cime_comp_mod, only : cime_pre_init2
  use cime_comp_mod, only : cime_init
  use cime_comp_mod, only : cime_run
  use cime_comp_mod, only : cime_final

  implicit none

  !--------------------------------------------------------------------------
  ! timing variables
  !--------------------------------------------------------------------------
  integer(i8) :: beg_count, end_count, irtc_rate
  real(r8)    :: cime_pre_init1_time, ESMF_Initialize_time, &
       cime_pre_init2_time, cime_init_time_adjustment

  !--------------------------------------------------------------------------
  ! Setup and initialize the communications and logging.
  !--------------------------------------------------------------------------
  beg_count = shr_sys_irtc(irtc_rate)

  call cime_pre_init1()

  end_count = shr_sys_irtc(irtc_rate)
  cime_pre_init1_time = real( (end_count-beg_count), r8)/real(irtc_rate, r8)

  !--------------------------------------------------------------------------
  ! Initialize ESMF.  This is done outside of the ESMF_INTERFACE ifdef
  ! because it is needed for the time manager, even if the ESMF_INTERFACE
  ! is not used.
  !--------------------------------------------------------------------------
  beg_count = shr_sys_irtc(irtc_rate)

  call ESMF_Initialize(logkindflag=esmf_logfile_kind)

  end_count = shr_sys_irtc(irtc_rate)
  ESMF_Initialize_time = real( (end_count-beg_count), r8)/real(irtc_rate, r8)

  !--------------------------------------------------------------------------
  ! Read in the configuration information and initialize the time manager.
  !--------------------------------------------------------------------------
  ! Timer initialization has to be after determination of the maximum number
  ! of threads used across all components, so called inside of
  ! cime_pre_init2, as are t_startf and t_stopf for CPL:INIT and
  ! cime_pre_init2.
  !--------------------------------------------------------------------------
  beg_count = shr_sys_irtc(irtc_rate)

  call cime_pre_init2()

  end_count = shr_sys_irtc(irtc_rate)
  cime_pre_init2_time = real( (end_count-beg_count), r8)/real(irtc_rate, r8)

  !--------------------------------------------------------------------------
  ! Call the initialize, run and finalize routines.
  !--------------------------------------------------------------------------

  call t_startf('CPL:INIT')
  call t_adj_detailf(+1)

  call t_startstop_valsf('CPL:cime_pre_init1',  walltime=cime_pre_init1_time)
  call t_startstop_valsf('CPL:ESMF_Initialize', walltime=ESMF_Initialize_time)
  call t_startstop_valsf('CPL:cime_pre_init2',  walltime=cime_pre_init2_time)

  call cime_init()

  call t_adj_detailf(-1)
  call t_stopf('CPL:INIT')

  cime_init_time_adjustment = cime_pre_init1_time  &
       + ESMF_Initialize_time &
       + cime_pre_init2_time
  call t_startstop_valsf('CPL:INIT',  walltime=cime_init_time_adjustment, &
       callcount=0)

  call cime_run()
  call cime_final()

  !--------------------------------------------------------------------------
  ! Clean-up
  !--------------------------------------------------------------------------
  call ESMF_Finalize( )

end program cime_driver
