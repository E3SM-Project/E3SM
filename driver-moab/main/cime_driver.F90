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
  use shr_kind_mod,  only : CS => SHR_KIND_CS
  use shr_sys_mod,   only : shr_sys_irtc, shr_sys_abort
  use perf_mod,      only : t_startf, t_adj_detailf, t_stopf, t_startstop_valsf
  use ESMF,          only : ESMF_Initialize, ESMF_Finalize
  use ESMF,          only : ESMF_LogKind_Flag, ESMF_LOGKIND_NONE
  use ESMF,          only : ESMF_LOGKIND_SINGLE, ESMF_LOGKIND_MULTI
#if (! defined(USE_ESMF_LIB) ) || (ESMF_VERSION_MAJOR > 7)
  use ESMF,          only : ESMF_LOGKIND_MULTI_ON_ERROR
#endif
  use cime_comp_mod, only : cime_pre_init1
  use cime_comp_mod, only : cime_pre_init2
  use cime_comp_mod, only : cime_init
  use cime_comp_mod, only : cime_run
  use cime_comp_mod, only : cime_final
  use cime_comp_mod, only : cime_pre_init2_lb
  use seq_comm_mct,  only : logunit

  implicit none

  !--------------------------------------------------------------------------
  ! timing variables
  !--------------------------------------------------------------------------
  integer(i8) :: beg_count, end_count, irtc_rate
  real(r8)    :: cime_pre_init1_time, ESMF_Initialize_time, &
                 cime_pre_init2_time, cime_pre_init2_time_adjustment, &
                 cime_init_time_adjustment

  !--------------------------------------------------------------------------
  ! For ESMF logging
  !--------------------------------------------------------------------------
  character(len=CS)       :: esmf_logfile_option
  type(ESMF_LogKind_Flag) :: esmf_logfile_kind

  !--------------------------------------------------------------------------
  ! Setup and initialize the communications and logging.
  !--------------------------------------------------------------------------
  beg_count = shr_sys_irtc(irtc_rate)

  call cime_pre_init1(esmf_logfile_option)

  end_count = shr_sys_irtc(irtc_rate)
  cime_pre_init1_time = real( (end_count-beg_count), r8)/real(irtc_rate, r8)

  !--------------------------------------------------------------------------
  ! Initialize ESMF.  This is done outside of the ESMF_INTERFACE ifdef
  ! because it is needed for the time manager, even if the ESMF_INTERFACE
  ! is not used.
  !--------------------------------------------------------------------------
  beg_count = shr_sys_irtc(irtc_rate)

  select case(esmf_logfile_option)
  case('ESMF_LOGKIND_SINGLE')
     esmf_logfile_kind = ESMF_LOGKIND_SINGLE
  case('ESMF_LOGKIND_MULTI')
     esmf_logfile_kind = ESMF_LOGKIND_MULTI
  case('ESMF_LOGKIND_MULTI_ON_ERROR')
#if (! defined(USE_ESMF_LIB) ) || (ESMF_VERSION_MAJOR > 7)
     esmf_logfile_kind = ESMF_LOGKIND_MULTI_ON_ERROR
#else
     write(logunit,*) 'ESMF library version being used: ', ESMF_VERSION_MAJOR
     call shr_sys_abort('CIME ERROR: invalid ESMF logfile kind for this ESMF library version: ' &
                        //trim(esmf_logfile_option) )
#endif
  case('ESMF_LOGKIND_NONE')
     esmf_logfile_kind = ESMF_LOGKIND_NONE
  case default
     call shr_sys_abort('CIME ERROR: invalid ESMF logfile kind '//trim(esmf_logfile_option))
  end select
  call ESMF_Initialize(logkindflag=esmf_logfile_kind)

  end_count = shr_sys_irtc(irtc_rate)
  ESMF_Initialize_time = real( (end_count-beg_count), r8)/real(irtc_rate, r8)

  !--------------------------------------------------------------------------
  ! Read in the configuration information and initialize the time manager.
  !--------------------------------------------------------------------------
  ! Timer initialization has to be after determination of the maximum number
  ! of threads used across all components, so called inside of
  ! cime_pre_init2, as are t_startf for CPL:INIT and t_startf and t_stopf for
  ! cime_pre_init2. t_startf/t_stopf are also called for cime_pre_init1 and
  ! cime_ESMF_initialize, to establish the correct parent/child
  ! relationships. t_stopf for CPL:INIT is called below.
  ! The wallclock times are later adjusted, as needed, for all of these
  ! timers, using the timings collected in cime_driver using shr_sys_irtc.
  !--------------------------------------------------------------------------
  beg_count = shr_sys_irtc(irtc_rate)

  call cime_pre_init2()

  end_count = shr_sys_irtc(irtc_rate)
  cime_pre_init2_time = real( (end_count-beg_count), r8)/real(irtc_rate, r8)

  !--------------------------------------------------------------------------
  ! Call the initialize routines
  !--------------------------------------------------------------------------
  call t_adj_detailf(+1)

  call cime_init()

  call t_adj_detailf(-1)

  !--------------------------------------------------------------------------
  ! Adjust initialization timers, as needed
  !--------------------------------------------------------------------------
  call t_startstop_valsf('CPL:cime_pre_init1',  walltime=cime_pre_init1_time, &
                         callcount = 0)
  call t_startstop_valsf('CPL:ESMF_Initialize', walltime=ESMF_Initialize_time, &
                         callcount = 0)

  ! Some of the time for cime_pre_init2 has already been recorded. Adjust for this.
  cime_pre_init2_time_adjustment = cime_pre_init2_time - cime_pre_init2_lb
  call t_startstop_valsf('CPL:cime_pre_init2',  walltime=cime_pre_init2_time_adjustment, &
                         callcount = 0)

  ! CPL:INIT timer started in cime_pre_init2
  call t_stopf('CPL:INIT')

  ! Add in the time that was not captured by the CPL:INIT timer (because it was
  ! started inside of cime_pre_init2 instead of before cime_pre_init1).
  cime_init_time_adjustment = cime_pre_init1_time  &
                            + ESMF_Initialize_time &
                            + cime_pre_init2_time_adjustment
  call t_startstop_valsf('CPL:INIT',  walltime=cime_init_time_adjustment, &
                         callcount=0)

  !--------------------------------------------------------------------------
  ! Call the run and finalize routines.
  !--------------------------------------------------------------------------
  call cime_run()
  call cime_final()

  !--------------------------------------------------------------------------
  ! Clean-up
  !--------------------------------------------------------------------------
  call ESMF_Finalize( )

end program cime_driver
