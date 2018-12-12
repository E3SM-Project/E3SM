module shr_nuopc_time_mod
  ! !USES:
  use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet, ESMF_GridCompSet
  use ESMF                  , only : ESMF_Clock, ESMF_ClockCreate, ESMF_ClockGet, ESMF_ClockSet
  use ESMF                  , only : ESMF_ClockAdvance
  use ESMF                  , only : ESMF_Alarm, ESMF_AlarmCreate, ESMF_AlarmGet
  use ESMF                  , only : ESMF_Calendar, ESMF_CalKind_Flag, ESMF_CalendarCreate
  use ESMF                  , only : ESMF_CALKIND_NOLEAP, ESMF_CALKIND_GREGORIAN
  use ESMF                  , only : ESMF_Time, ESMF_TimeGet, ESMF_TimeSet
  use ESMF                  , only : ESMF_TimeInterval, ESMF_TimeIntervalSet, ESMF_TimeIntervalGet
  use ESMF                  , only : ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_FAILURE
  use ESMF                  , only : ESMF_VM, ESMF_VMGet, ESMF_VMBroadcast
  use ESMF                  , only : operator(<), operator(/=), operator(+)
  use ESMF                  , only : operator(-), operator(*) , operator(>=)
  use ESMF                  , only : operator(<=), operator(>), operator(==)
  use NUOPC                 , only : NUOPC_CompAttributeGet
  use med_constants_mod     , only : dbug_flag => med_constants_dbug_flag
  use shr_nuopc_utils_mod   , only : shr_nuopc_abort
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr

  implicit none
  private    ! default private

  public  :: shr_nuopc_time_alarmInit  ! initialize an alarm
  public  :: shr_nuopc_time_clockInit  ! initialize driver clock
  public  :: shr_nuopc_time_set_component_stop_alarm

  private :: shr_nuopc_time_timeInit
  private :: shr_nuopc_time_date2ymd

  ! Clock and alarm options
  character(len=*), private, parameter :: &
       optNONE           = "none"      , &
       optNever          = "never"     , &
       optNSteps         = "nsteps"    , &
       optNStep          = "nstep"     , &
       optNSeconds       = "nseconds"  , &
       optNSecond        = "nsecond"   , &
       optNMinutes       = "nminutes"  , &
       optNMinute        = "nminute"   , &
       optNHours         = "nhours"    , &
       optNHour          = "nhour"     , &
       optNDays          = "ndays"     , &
       optNDay           = "nday"      , &
       optNMonths        = "nmonths"   , &
       optNMonth         = "nmonth"    , &
       optNYears         = "nyears"    , &
       optNYear          = "nyear"     , &
       optMonthly        = "monthly"   , &
       optYearly         = "yearly"    , &
       optDate           = "date"      , &
       optIfdays0        = "ifdays0"   , &
       optGLCCouplingPeriod = "glc_coupling_period"

  ! Module data
  integer, parameter          :: SecPerDay = 86400 ! Seconds per day
  character(len=*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine shr_nuopc_time_clockInit(ensemble_driver, esmdriver, logunit, rc)

    use med_constants_mod     , only : CL, CS
    use shr_file_mod          , only : shr_file_getUnit, shr_file_freeUnit
    use shr_cal_mod           , only : shr_cal_noleap, shr_cal_gregorian, shr_cal_calendarname

    ! input/output variables
    type(ESMF_GridComp)  :: ensemble_driver, esmdriver
    integer, intent(in)  :: logunit
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)        :: clock
    type(ESMF_VM)           :: vm
    type(ESMF_Time)         :: StartTime           ! Start time
    type(ESMF_Time)         :: RefTime             ! Reference time
    type(ESMF_Time)         :: CurrTime            ! Current time
    type(ESMF_Time)         :: StopTime            ! Stop time
    type(ESMF_Time)         :: StopTime1           ! Stop time
    type(ESMF_Time)         :: StopTime2           ! Stop time
    type(ESMF_Time)         :: Clocktime           ! Loop time
    type(ESMF_TimeInterval) :: TimeStep            ! Clock time-step
    type(ESMF_Calendar)     :: calendar            ! esmf calendar
    type(ESMF_CalKind_Flag) :: caltype             ! esmf calendar type
    type(ESMF_Alarm)        :: alarm_stop          ! alarm
    type(ESMF_Alarm)        :: alarm_datestop      ! alarm
    integer                 :: ref_ymd             ! Reference date (YYYYMMDD)
    integer                 :: ref_tod             ! Reference time of day (seconds)
    integer                 :: start_ymd           ! Start date (YYYYMMDD)
    integer                 :: start_tod           ! Start time of day (seconds)
    integer                 :: curr_ymd            ! Current ymd (YYYYMMDD)
    integer                 :: curr_tod            ! Current tod (seconds)
    integer                 :: stop_n              ! Number until stop
    integer                 :: stop_ymd            ! Stop date (YYYYMMDD)
    integer                 :: stop_tod            ! Stop time-of-day
    character(CS)           :: stop_option         ! Stop option units
    integer                 :: atm_cpl_dt          ! Atmosphere coupling interval
    integer                 :: lnd_cpl_dt          ! Land coupling interval
    integer                 :: ice_cpl_dt          ! Sea-Ice coupling interval
    integer                 :: ocn_cpl_dt          ! Ocean coupling interval
    integer                 :: glc_cpl_dt          ! Glc coupling interval
    integer                 :: rof_cpl_dt          ! Runoff coupling interval
    integer                 :: wav_cpl_dt          ! Wav coupling interval
    integer                 :: esp_cpl_dt          ! Esp coupling interval
    character(CS)           :: glc_avg_period      ! Glc avering coupling period
    logical                 :: read_restart
    character(len=CL)       :: restart_file
    character(len=CL)       :: restart_pfile
    character(len=CL)       :: cvalue
    integer                 :: dtime_drv           ! time-step to use
    integer                 :: yr, mon, day        ! Year, month, day as integers
    integer                 :: localPet            ! local pet in esm domain
    logical                 :: mastertask          ! true if mastertask in esm domain
    integer                 :: unitn               ! unit number
    integer                 :: ierr                ! Return code
    character(CL)           :: tmpstr              ! temporary
    character(CS)           :: calendar_name       ! Calendar name
    character(CS)           :: inst_suffix
    integer                 :: tmp(6)              ! Array for Broadcast
    integer                 :: dbrc
    logical                 :: isPresent
    character(len=*), parameter :: subname = '(shr_nuopc_time_clockInit): '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif

    call ESMF_GridCompGet(esmdriver, vm=vm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    ! We may want to get the ensemble_driver vm here instead so that
    ! files are read on global task 0 only instead of each esm member task 0
    call ESMF_VMGet(vm, localPet=localPet, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    mastertask = localPet == 0
    !---------------------------------------------------------------------------
    ! Create the driver calendar
    !---------------------------------------------------------------------------

    call NUOPC_CompAttributeGet(esmdriver, name="calendar", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    calendar_name = shr_cal_calendarName(cvalue)

    if ( trim(calendar_name) == trim(shr_cal_noleap)) then
       caltype = ESMF_CALKIND_NOLEAP
    else if ( trim(calendar_name) == trim(shr_cal_gregorian)) then
       caltype = ESMF_CALKIND_GREGORIAN
    else
       call ESMF_LogWrite(trim(subname)//': unrecognized ESMF calendar specified: '//&
            trim(calendar_name), ESMF_LOGMSG_INFO, rc=rc)
       rc = ESMF_FAILURE
       return
    end if

    call ESMF_LogWrite(trim(subname)//': driver calendar is : '// trim(calendar_name), &
         ESMF_LOGMSG_INFO, rc=rc)

    calendar = ESMF_CalendarCreate( name='CMEPS_'//trim(calendar_name), &
         calkindflag=caltype, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------------------------------------------
    ! Determine clock start time, reference time and current time
    !---------------------------------------------------------------------------

    curr_ymd = 0
    curr_tod = 0

    call NUOPC_CompAttributeGet(esmdriver, name="start_ymd", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) start_ymd
    call NUOPC_CompAttributeGet(esmdriver, name="start_tod", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) start_tod

    call NUOPC_CompAttributeGet(esmdriver, name="ref_ymd", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) ref_ymd
    call NUOPC_CompAttributeGet(esmdriver, name="ref_tod", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) ref_tod

    call NUOPC_CompAttributeGet(esmdriver, name='read_restart', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) read_restart

    if (read_restart) then

       call NUOPC_CompAttributeGet(esmdriver, name='restart_file', value=restart_file, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       !--- read rpointer if restart_file is set to str_undefined ---
       if (trim(restart_file) == 'str_undefined') then

          ! Error check on restart_pfile
          call NUOPC_CompAttributeGet(esmdriver, name="restart_pfile", value=restart_pfile, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          call NUOPC_CompAttributeGet(esmdriver, name="inst_suffix", isPresent=isPresent, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          if(isPresent) then
             call NUOPC_CompAttributeGet(esmdriver, name="inst_suffix", value=inst_suffix, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          else
             inst_suffix = ""
          endif
          if ( len_trim(restart_pfile) == 0 ) then
             rc = ESMF_FAILURE
             call ESMF_LogWrite(trim(subname)//' ERROR restart_pfile must be defined', &
                  ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=dbrc)
             return
          end if
          restart_pfile = trim(restart_pfile)//inst_suffix
          if (mastertask) then
             unitn = shr_file_getUnit()
             call ESMF_LogWrite(trim(subname)//" read rpointer file = "//trim(restart_pfile), &
                  ESMF_LOGMSG_INFO, rc=dbrc)
             open(unitn, file=restart_pfile, form='FORMATTED', status='old',iostat=ierr)
             if (ierr < 0) then
                rc = ESMF_FAILURE
                call ESMF_LogWrite(trim(subname)//' ERROR rpointer file open returns error', &
                     ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=dbrc)
                return
             end if
             read(unitn,'(a)', iostat=ierr) restart_file
             if (ierr < 0) then
                rc = ESMF_FAILURE
                call ESMF_LogWrite(trim(subname)//' ERROR rpointer file read returns error', &
                     ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=dbrc)
                return
             end if
             close(unitn)
             call shr_file_freeUnit( unitn )
             call ESMF_LogWrite(trim(subname)//" read driver restart from file = "//trim(restart_file), &
                  ESMF_LOGMSG_INFO, rc=dbrc)
          endif
       endif
       if (mastertask) then
          call shr_nuopc_time_read_restart_calendar_settings(restart_file, &
               start_ymd, start_tod, ref_ymd, ref_tod, curr_ymd, curr_tod)
       endif
       tmp(1) = start_ymd
       tmp(2) = start_tod
       tmp(3) = ref_ymd
       tmp(4) = ref_tod
       tmp(5) = curr_ymd
       tmp(6) = curr_tod
       call ESMF_VMBroadcast(vm, tmp, 6, 0, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       start_ymd = tmp(1)
       start_tod = tmp(2)
       ref_ymd = tmp(3)
       ref_tod = tmp(4)
       curr_ymd = tmp(5)
       curr_tod = tmp(6)
    end if

    if ( ref_ymd == 0 ) then
       ref_ymd = start_ymd
       ref_tod = start_tod
    endif
    if ( curr_ymd == 0 ) then
       curr_ymd = start_ymd
       curr_tod = start_tod
    endif

    ! Determine start time
    call shr_nuopc_time_date2ymd(start_ymd, yr, mon, day)
    call ESMF_TimeSet( StartTime, yy=yr, mm=mon, dd=day, s=start_tod, calendar=calendar, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if(mastertask .or. dbug_flag > 2) then
       write(tmpstr,'(i10)') start_ymd
       call ESMF_LogWrite(trim(subname)//': driver start_ymd: '// trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
       write(logunit,*)   trim(subname)//': driver start_ymd: '// trim(tmpstr)
       write(tmpstr,'(i10)') start_tod
       call ESMF_LogWrite(trim(subname)//': driver start_tod: '// trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
       write(logunit,*)   trim(subname)//': driver start_tod: '// trim(tmpstr)
    endif

    ! Determine reference time
    call shr_nuopc_time_date2ymd(ref_ymd, yr, mon, day)
    call ESMF_TimeSet( RefTime, yy=yr, mm=mon, dd=day, s=ref_tod, calendar=calendar, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if(mastertask .or. dbug_flag > 2) then
       write(tmpstr,'(i10)') ref_ymd
       call ESMF_LogWrite(trim(subname)//': driver ref_ymd: '// trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
       write(logunit,*)   trim(subname)//': driver ref_ymd: '// trim(tmpstr)
       write(tmpstr,'(i10)') ref_tod
       call ESMF_LogWrite(trim(subname)//': driver ref_tod: '// trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
       write(logunit,*)   trim(subname)//': driver ref_tod: '// trim(tmpstr)
    endif
    ! Determine current time
    call shr_nuopc_time_date2ymd(curr_ymd, yr, mon, day)
    call ESMF_TimeSet( CurrTime, yy=yr, mm=mon, dd=day, s=curr_tod, calendar=calendar, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if(mastertask .or. dbug_flag > 2) then
       write(tmpstr,'(i10)') curr_ymd
       call ESMF_LogWrite(trim(subname)//': driver curr_ymd: '// trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
       write(logunit,*)   trim(subname)//': driver curr_ymd: '// trim(tmpstr)
       write(tmpstr,'(i10)') curr_tod
       call ESMF_LogWrite(trim(subname)//': driver curr_tod: '// trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
       write(logunit,*)   trim(subname)//': driver curr_tod: '// trim(tmpstr)
    endif
    !---------------------------------------------------------------------------
    ! Determine driver clock timestep
    !---------------------------------------------------------------------------

    call NUOPC_CompAttributeGet(esmdriver, name="atm_cpl_dt", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) atm_cpl_dt

    call NUOPC_CompAttributeGet(esmdriver, name="lnd_cpl_dt", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) lnd_cpl_dt

    call NUOPC_CompAttributeGet(esmdriver, name="ice_cpl_dt", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) ice_cpl_dt

    call NUOPC_CompAttributeGet(esmdriver, name="ocn_cpl_dt", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) ocn_cpl_dt

    call NUOPC_CompAttributeGet(esmdriver, name="glc_cpl_dt", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) glc_cpl_dt

    call NUOPC_CompAttributeGet(esmdriver, name="rof_cpl_dt", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) rof_cpl_dt

    call NUOPC_CompAttributeGet(esmdriver, name="wav_cpl_dt", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) wav_cpl_dt

    call NUOPC_CompAttributeGet(esmdriver, name="glc_avg_period", value=glc_avg_period, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) glc_avg_period

    ! TODO: for now - this is not in the namelist_definition_drv.xml file
    ! call NUOPC_CompAttributeGet(esmdriver, name="esp_cpl_dt", value=cvalue, rc=rc)
    ! if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    ! read(cvalue,*) esp_cpl_dt
    esp_cpl_dt = 9999

    dtime_drv = 9999
    dtime_drv = min(dtime_drv, atm_cpl_dt)
    dtime_drv = min(dtime_drv, lnd_cpl_dt)
    dtime_drv = min(dtime_drv, ocn_cpl_dt)
    dtime_drv = min(dtime_drv, ice_cpl_dt)
    dtime_drv = min(dtime_drv, glc_cpl_dt)
    dtime_drv = min(dtime_drv, rof_cpl_dt)
    dtime_drv = min(dtime_drv, wav_cpl_dt)
    dtime_drv = min(dtime_drv, esp_cpl_dt)
    if(mastertask .or. dbug_flag > 2) then
       write(tmpstr,'(i10)') dtime_drv
       call ESMF_LogWrite(trim(subname)//': driver time interval is : '// trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
       write(logunit,*)   trim(subname)//': driver time interval is : '// trim(tmpstr)
    endif
    call ESMF_TimeIntervalSet( TimeStep, s=dtime_drv, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------------------------------------------
    ! Create the driver clock with an artificial stop time
    !---------------------------------------------------------------------------

    ! Create the clock
    clock = ESMF_ClockCreate(TimeStep, StartTime, refTime=RefTime, name='ESMF Driver Clock', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Advance the clock to the current time (in case of a restart)
    call ESMF_ClockGet(clock, currTime=clocktime, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    do while( clocktime < CurrTime)
       call ESMF_ClockAdvance( clock, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_ClockGet( clock, currTime=clocktime, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end do

    ! Set the driver gridded component clock to the created clock
    call ESMF_GridCompSet(esmdriver, clock=clock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !-------------------------------
    ! Set driver clock stop time
    !-------------------------------

    call NUOPC_CompAttributeGet(esmdriver, name="stop_option", value=stop_option, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompAttributeGet(esmdriver, name="stop_n", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) stop_n
    call NUOPC_CompAttributeGet(esmdriver, name="stop_ymd", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) stop_ymd
    call NUOPC_CompAttributeGet(esmdriver, name="stop_tod", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) stop_tod
    if ( stop_ymd < 0) then
       stop_ymd = 99990101
       stop_tod = 0
    endif
    if(mastertask .or. dbug_flag > 2) then
       write(tmpstr,'(i10)') stop_ymd
       call ESMF_LogWrite(trim(subname)//': driver stop_ymd: '// trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
       write(logunit,*)   trim(subname)//': driver stop_ymd: '// trim(tmpstr)
       write(tmpstr,'(i10)') stop_tod
       call ESMF_LogWrite(trim(subname)//': driver stop_tod: '// trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
       write(logunit,*)   trim(subname)//': driver stop_tod: '// trim(tmpstr)
    endif
    call shr_nuopc_time_alarmInit(clock, &
         alarm   = alarm_stop,           &
         option  = stop_option,          &
         opt_n   = stop_n,               &
         opt_ymd = stop_ymd,             &
         opt_tod = stop_tod,             &
         RefTime = CurrTime,             &
         alarmname = 'alarm_stop', rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_time_alarmInit(clock, &
         alarm     = alarm_datestop,     &
         option    = optDate,            &
         opt_ymd   = stop_ymd,           &
         opt_tod   = stop_tod,           &
         RefTime   = StartTime,          &
         alarmname = 'alarm_datestop', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_AlarmGet(alarm_stop, RingTime=StopTime1, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_AlarmGet(alarm_datestop, RingTime=StopTime2, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (StopTime2 < StopTime1) then
       StopTime = StopTime2
    else
       StopTime = StopTime1
    endif

    call ESMF_ClockSet(clock, StopTime=StopTime, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Create the ensemble driver clock
    TimeStep = StopTime-ClockTime
    clock = ESMF_ClockCreate(TimeStep, ClockTime, StopTime=StopTime, &
         refTime=RefTime, name='ESMF ensemble Driver Clock', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_GridCompSet(ensemble_driver, clock=clock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return



 end subroutine shr_nuopc_time_clockInit

 subroutine shr_nuopc_time_set_component_stop_alarm(gcomp, rc)
   use ESMF, only : ESMF_GridComp, ESMF_Alarm, ESMF_Clock, ESMF_ClockGet
   use ESMF, only : ESMF_AlarmSet
   use NUOPC, only : NUOPC_CompAttributeGet
   use NUOPC_Model, only : NUOPC_ModelGet
   type(ESMF_gridcomp) :: gcomp

   character(len=256)       :: stop_option    ! Stop option units
   integer                  :: stop_n         ! Number until stop interval
   integer                  :: stop_ymd       ! Stop date (YYYYMMDD)
   type(ESMF_ALARM)         :: stop_alarm
   character(len=256)       :: cvalue
   type(ESMF_Clock)         :: mclock
   type(ESMF_Time)          :: mCurrTime
   integer                  :: rc
   !----------------
   ! Stop alarm
   !----------------
   call NUOPC_CompAttributeGet(gcomp, name="stop_option", value=stop_option, rc=rc)
   if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

   call NUOPC_ModelGet(gcomp, modelClock=mclock, rc=rc)
   if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

   call ESMF_ClockGet(mclock, CurrTime=mCurrTime, rc=rc)
   if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

   call NUOPC_CompAttributeGet(gcomp, name="stop_n", value=cvalue, rc=rc)
   if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
   read(cvalue,*) stop_n

   call NUOPC_CompAttributeGet(gcomp, name="stop_ymd", value=cvalue, rc=rc)
   if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
   read(cvalue,*) stop_ymd
   call shr_nuopc_time_alarmInit(mclock, stop_alarm, stop_option, &
        opt_n   = stop_n,           &
        opt_ymd = stop_ymd,         &
        RefTime = mcurrTime,           &
        alarmname = 'alarm_stop', rc=rc)
   if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

   call ESMF_AlarmSet(stop_alarm, clock=mclock, rc=rc)
   if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
 end subroutine shr_nuopc_time_set_component_stop_alarm

!===============================================================================

  subroutine shr_nuopc_time_alarmInit( clock, alarm, option, &
       opt_n, opt_ymd, opt_tod, RefTime, alarmname, rc)

    ! !DESCRIPTION: Setup an alarm in a clock
    ! Notes: The ringtime sent to AlarmCreate MUST be the next alarm
    ! time.  If you send an arbitrary but proper ringtime from the
    ! past and the ring interval, the alarm will always go off on the
    ! next clock advance and this will cause serious problems.  Even
    ! if it makes sense to initialize an alarm with some reference
    ! time and the alarm interval, that reference time has to be
    ! advance forward to be >= the current time.  In the logic below
    ! we set an appropriate "NextAlarm" and then we make sure to
    ! advance it properly based on the ring interval.

    ! input/output variables
    type(ESMF_Clock)            , intent(inout) :: clock     ! clock
    type(ESMF_Alarm)            , intent(inout) :: alarm     ! alarm
    character(len=*)            , intent(in)    :: option    ! alarm option
    integer          , optional , intent(in)    :: opt_n     ! alarm freq
    integer          , optional , intent(in)    :: opt_ymd   ! alarm ymd
    integer          , optional , intent(in)    :: opt_tod   ! alarm tod (sec)
    type(ESMF_Time)  , optional , intent(in)    :: RefTime   ! ref time
    character(len=*) , optional , intent(in)    :: alarmname ! alarm name
    integer                     , intent(inout) :: rc        ! Return code

    ! local variables
    type(ESMF_Calendar)     :: cal                ! calendar
    integer                 :: lymd             ! local ymd
    integer                 :: ltod             ! local tod
    integer                 :: cyy,cmm,cdd,csec ! time info
    character(len=64)       :: lalarmname       ! local alarm name
    logical                 :: update_nextalarm ! update next alarm
    type(ESMF_Time)         :: CurrTime         ! Current Time
    type(ESMF_Time)         :: NextAlarm        ! Next restart alarm time
    type(ESMF_TimeInterval) :: AlarmInterval    ! Alarm interval
    integer                 :: sec
    character(len=*), parameter :: subname = '(shr_nuopc_time_alarmInit): '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    lalarmname = 'alarm_unknown'
    if (present(alarmname)) lalarmname = trim(alarmname)
    ltod = 0
    if (present(opt_tod)) ltod = opt_tod
    lymd = -1
    if (present(opt_ymd)) lymd = opt_ymd

    call ESMF_ClockGet(clock, CurrTime=CurrTime, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet(CurrTime, yy=cyy, mm=cmm, dd=cdd, s=csec, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! initial guess of next alarm, this will be updated below
    if (present(RefTime)) then
       NextAlarm = RefTime
    else
       NextAlarm = CurrTime
    endif

    ! Determine calendar
    call ESMF_ClockGet(clock, calendar=cal)

    ! Determine inputs for call to create alarm
    selectcase (trim(option))

    case (optNONE)
       call ESMF_TimeIntervalSet(AlarmInterval, yy=9999, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeSet( NextAlarm, yy=9999, mm=12, dd=1, s=0, calendar=cal, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       update_nextalarm  = .false.

    case (optNever)
       call ESMF_TimeIntervalSet(AlarmInterval, yy=9999, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeSet( NextAlarm, yy=9999, mm=12, dd=1, s=0, calendar=cal, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       update_nextalarm  = .false.

    case (optDate)
       if (.not. present(opt_ymd)) then
          call shr_nuopc_abort(subname//trim(option)//' requires opt_ymd')
       end if
       if (lymd < 0 .or. ltod < 0) then
          call shr_nuopc_abort(subname//trim(option)//'opt_ymd, opt_tod invalid')
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, yy=9999, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_time_timeInit(NextAlarm, lymd, cal, tod=ltod, desc="optDate")
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       update_nextalarm  = .false.

    case (optIfdays0)
       if (.not. present(opt_ymd)) then
          call shr_nuopc_abort(subname//trim(option)//' requires opt_ymd')
       end if
       if (.not.present(opt_n)) then
          call shr_nuopc_abort(subname//trim(option)//' requires opt_n')
       end if
       if (opt_n <= 0)  then
          call shr_nuopc_abort(subname//trim(option)//' invalid opt_n')
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeSet( NextAlarm, yy=cyy, mm=cmm, dd=opt_n, s=0, calendar=cal, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       update_nextalarm  = .true.

    case (optNSteps)
       if (.not.present(opt_n)) then
          call shr_nuopc_abort(subname//trim(option)//' requires opt_n')
       end if
       if (opt_n <= 0) then
          call shr_nuopc_abort(subname//trim(option)//' invalid opt_n')
       end if
       call ESMF_ClockGet(clock, TimeStep=AlarmInterval, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNStep)
       if (.not.present(opt_n)) call shr_nuopc_abort(subname//trim(option)//' requires opt_n')
       if (opt_n <= 0)  call shr_nuopc_abort(subname//trim(option)//' invalid opt_n')
       call ESMF_ClockGet(clock, TimeStep=AlarmInterval, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNSeconds)
       if (.not.present(opt_n)) then
          call shr_nuopc_abort(subname//trim(option)//' requires opt_n')
       end if
       if (opt_n <= 0) then
          call shr_nuopc_abort(subname//trim(option)//' invalid opt_n')
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, s=1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNSecond)
       if (.not.present(opt_n)) then
          call shr_nuopc_abort(subname//trim(option)//' requires opt_n')
       end if
       if (opt_n <= 0) then
          call shr_nuopc_abort(subname//trim(option)//' invalid opt_n')
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, s=1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNMinutes)
       call ESMF_TimeIntervalSet(AlarmInterval, s=60, rc=rc)
       if (.not.present(opt_n)) then
          call shr_nuopc_abort(subname//trim(option)//' requires opt_n')
       end if
       if (opt_n <= 0) then
          call shr_nuopc_abort(subname//trim(option)//' invalid opt_n')
       end if
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNMinute)
       if (.not.present(opt_n)) then
          call shr_nuopc_abort(subname//trim(option)//' requires opt_n')
       end if
       if (opt_n <= 0) then
          call shr_nuopc_abort(subname//trim(option)//' invalid opt_n')
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, s=60, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNHours)
       if (.not.present(opt_n)) then
          call shr_nuopc_abort(subname//trim(option)//' requires opt_n')
       end if
       if (opt_n <= 0) then
          call shr_nuopc_abort(subname//trim(option)//' invalid opt_n')
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, s=3600, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNHour)
       if (.not.present(opt_n)) then
          call shr_nuopc_abort(subname//trim(option)//' requires opt_n')
       end if
       if (opt_n <= 0) then
          call shr_nuopc_abort(subname//trim(option)//' invalid opt_n')
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, s=3600, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNDays)
       if (.not.present(opt_n)) then
          call shr_nuopc_abort(subname//trim(option)//' requires opt_n')
       end if
       if (opt_n <= 0) then
          call shr_nuopc_abort(subname//trim(option)//' invalid opt_n')
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, d=1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNDay)
       if (.not.present(opt_n)) then
          call shr_nuopc_abort(subname//trim(option)//' requires opt_n')
       end if
       if (opt_n <= 0) then
          call shr_nuopc_abort(subname//trim(option)//' invalid opt_n')
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, d=1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNMonths)
       if (.not.present(opt_n)) then
          call shr_nuopc_abort(subname//trim(option)//' requires opt_n')
       end if
       if (opt_n <= 0) then
          call shr_nuopc_abort(subname//trim(option)//' invalid opt_n')
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNMonth)
       if (.not.present(opt_n)) then
          call shr_nuopc_abort(subname//trim(option)//' requires opt_n')
       end if
       if (opt_n <= 0) then
          call shr_nuopc_abort(subname//trim(option)//' invalid opt_n')
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optMonthly)
       call ESMF_TimeIntervalSet(AlarmInterval, mm=1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeSet( NextAlarm, yy=cyy, mm=cmm, dd=1, s=0, calendar=cal, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       update_nextalarm  = .true.

    case (optNYears)
       if (.not.present(opt_n)) then
          call shr_nuopc_abort(subname//trim(option)//' requires opt_n')
       end if
       if (opt_n <= 0) then
          call shr_nuopc_abort(subname//trim(option)//' invalid opt_n')
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, yy=1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optNYear)
       if (.not.present(opt_n)) then
          call shr_nuopc_abort(subname//trim(option)//' requires opt_n')
       end if
       if (opt_n <= 0) then
          call shr_nuopc_abort(subname//trim(option)//' invalid opt_n')
       end if
       call ESMF_TimeIntervalSet(AlarmInterval, yy=1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       AlarmInterval = AlarmInterval * opt_n
       update_nextalarm  = .true.

    case (optYearly)
       call ESMF_TimeIntervalSet(AlarmInterval, yy=1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeSet( NextAlarm, yy=cyy, mm=1, dd=1, s=0, calendar=cal, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       update_nextalarm  = .true.

    case default
       call shr_nuopc_abort(subname//'unknown option '//trim(option))

    end select

    ! --------------------------------------------------------------------------------
    ! --- AlarmInterval and NextAlarm should be set ---
    ! --------------------------------------------------------------------------------

    ! --- advance Next Alarm so it won't ring on first timestep for
    ! --- most options above. go back one alarminterval just to be careful

    if (update_nextalarm) then
       NextAlarm = NextAlarm - AlarmInterval
       do while (NextAlarm <= CurrTime)
          NextAlarm = NextAlarm + AlarmInterval
       enddo
    endif

    alarm = ESMF_AlarmCreate( name=lalarmname, clock=clock, ringTime=NextAlarm, &
         ringInterval=AlarmInterval, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine shr_nuopc_time_alarmInit

  !===============================================================================

  subroutine shr_nuopc_time_timeInit( Time, ymd, cal, tod, desc, logunit )

    !  Create the ESMF_Time object corresponding to the given input time, given in
    !  YMD (Year Month Day) and TOD (Time-of-day) format.
    !  Set the time by an integer as YYYYMMDD and integer seconds in the day

    ! input/output parameters:
    type(ESMF_Time)     , intent(inout)        :: Time ! ESMF time
    integer             , intent(in)           :: ymd  ! year, month, day YYYYMMDD
    type(ESMF_Calendar) , intent(in)           :: cal  ! ESMF calendar
    integer             , intent(in), optional :: tod  ! time of day in seconds
    character(len=*)    , intent(in), optional :: desc ! description of time to set
    integer             , intent(in), optional :: logunit

    ! local variables
    integer                     :: yr, mon, day ! Year, month, day as integers
    integer                     :: ltod         ! local tod
    character(len=256)          :: ldesc        ! local desc
    integer                     :: rc           ! return code
    character(len=*), parameter :: subname = '(shr_nuopc_time_m_ETimeInit) '
    !-------------------------------------------------------------------------------

    ltod = 0
    if (present(tod)) ltod = tod
    ldesc = ''
    if (present(desc)) ldesc = desc

    if ( (ymd < 0) .or. (ltod < 0) .or. (ltod > SecPerDay) )then
       if (present(logunit)) then
          write(logunit,*) subname//': ERROR yymmdd is a negative number or '// &
               'time-of-day out of bounds', ymd, ltod
       end if
       call shr_nuopc_abort( subname//'ERROR: Bad input' )
    end if

    call shr_nuopc_time_date2ymd (ymd,yr,mon,day)

    call ESMF_TimeSet( Time, yy=yr, mm=mon, dd=day, s=ltod, calendar=cal, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine shr_nuopc_time_timeInit

  !===============================================================================

  subroutine shr_nuopc_time_date2ymd (date, year, month, day)

    ! input/output variables
    integer, intent(in)  :: date             ! coded-date (yyyymmdd)
    integer, intent(out) :: year,month,day   ! calendar year,month,day

    ! local variables
    integer :: tdate   ! temporary date
    character(*),parameter :: subName = "(shr_nuopc_time_date2ymd)"
    !-------------------------------------------------------------------------------

    tdate = abs(date)
    year = int(tdate/10000)
    if (date < 0) then
       year = -year
    end if
    month = int( mod(tdate,10000)/  100)
    day = mod(tdate,  100)

  end subroutine shr_nuopc_time_date2ymd

  subroutine shr_nuopc_time_read_restart_calendar_settings(restart_file, &
       start_ymd, start_tod, ref_ymd, ref_tod, curr_ymd, curr_tod)

    use netcdf                , only : nf90_open, nf90_nowrite, nf90_noerr
    use netcdf                , only : nf90_inq_varid, nf90_get_var, nf90_close
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO
    use med_constants_mod     , only : CL

    character(len=*), intent(in) :: restart_file
    integer, intent(out)         :: ref_ymd             ! Reference date (YYYYMMDD)
    integer, intent(out)         :: ref_tod             ! Reference time of day (seconds)
    integer, intent(out)         :: start_ymd           ! Start date (YYYYMMDD)
    integer, intent(out)         :: start_tod           ! Start time of day (seconds)
    integer, intent(out)         :: curr_ymd            ! Current ymd (YYYYMMDD)
    integer, intent(out)         :: curr_tod            ! Current tod (seconds)

    integer                 :: status, ncid, varid ! netcdf stuff
    integer                 :: dbrc                ! error codes
    character(CL)           :: tmpstr              ! temporary
    character(len=*), parameter :: subname = "(shr_nuopc_time_read_restart_calendar_settings)"

    ! use netcdf here since it's serial
    status = nf90_open(restart_file, NF90_NOWRITE, ncid)
    if (status /= nf90_NoErr) then
       print *,__FILE__,__LINE__,trim(restart_file)
       call shr_nuopc_abort(trim(subname)//' ERROR: nf90_open')
    endif
    status = nf90_inq_varid(ncid, 'start_ymd', varid)
    if (status /= nf90_NoErr) call shr_nuopc_abort(trim(subname)//' ERROR: nf90_inq_varid start_ymd')
    status = nf90_get_var(ncid, varid, start_ymd)
    if (status /= nf90_NoErr) call shr_nuopc_abort(trim(subname)//' ERROR: nf90_get_var start_ymd')
    status = nf90_inq_varid(ncid, 'start_tod', varid)
    if (status /= nf90_NoErr) call shr_nuopc_abort(trim(subname)//' ERROR: nf90_inq_varid start_tod')
    status = nf90_get_var(ncid, varid, start_tod)
    if (status /= nf90_NoErr) call shr_nuopc_abort(trim(subname)//' ERROR: nf90_get_var start_tod')
    status = nf90_inq_varid(ncid, 'ref_ymd', varid)
    if (status /= nf90_NoErr) call shr_nuopc_abort(trim(subname)//' ERROR: nf90_inq_varid ref_ymd')
    status = nf90_get_var(ncid, varid, ref_ymd)
    if (status /= nf90_NoErr) call shr_nuopc_abort(trim(subname)//' ERROR: nf90_get_var ref_ymd')
    status = nf90_inq_varid(ncid, 'ref_tod', varid)
    if (status /= nf90_NoErr) call shr_nuopc_abort(trim(subname)//' ERROR: nf90_inq_varid ref_tod')
    status = nf90_get_var(ncid, varid, ref_tod)
    if (status /= nf90_NoErr) call shr_nuopc_abort(trim(subname)//' ERROR: nf90_get_var ref_tod')
    status = nf90_inq_varid(ncid, 'curr_ymd', varid)
    if (status /= nf90_NoErr) call shr_nuopc_abort(trim(subname)//' ERROR: nf90_inq_varid curr_ymd')
    status = nf90_get_var(ncid, varid, curr_ymd)
    if (status /= nf90_NoErr) call shr_nuopc_abort(trim(subname)//' ERROR: nf90_get_var curr_ymd')
    status = nf90_inq_varid(ncid, 'curr_tod', varid)
    if (status /= nf90_NoErr) call shr_nuopc_abort(trim(subname)//' ERROR: nf90_inq_varid curr_tod')
    status = nf90_get_var(ncid, varid, curr_tod)
    if (status /= nf90_NoErr) call shr_nuopc_abort(trim(subname)//' ERROR: nf90_get_var curr_tod')
    status = nf90_close(ncid)
    if (status /= nf90_NoErr) call shr_nuopc_abort(trim(subname)//' ERROR: nf90_close')

    write(tmpstr,*) trim(subname)//" read start_ymd = ",start_ymd
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    write(tmpstr,*) trim(subname)//" read start_tod = ",start_tod
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    write(tmpstr,*) trim(subname)//" read ref_ymd   = ",ref_ymd
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    write(tmpstr,*) trim(subname)//" read ref_tod   = ",ref_tod
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    write(tmpstr,*) trim(subname)//" read curr_ymd  = ",curr_ymd
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
    write(tmpstr,*) trim(subname)//" read curr_tod  = ",curr_tod
    call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine shr_nuopc_time_read_restart_calendar_settings

end module shr_nuopc_time_mod
