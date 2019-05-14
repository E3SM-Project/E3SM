module med_phases_history_mod

  !-----------------------------------------------------------------------------
  ! Mediator Phases
  !-----------------------------------------------------------------------------

  use ESMF              , only : ESMF_Alarm
  use med_constants_mod , only : R8, CL, CS, I8

  implicit none
  private

  public :: med_phases_history_write

  type(ESMF_Alarm)        :: AlarmHist
  type(ESMF_Alarm)        :: AlarmHistAvg
  character(*), parameter :: u_FILE_u  = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine med_phases_history_write(gcomp, rc)

    ! Write mediator history file

    use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet
    use ESMF                  , only : ESMF_VM, ESMF_VMGet
    use ESMF                  , only : ESMF_Clock, ESMF_ClockGet, ESMF_ClockGetNextTime
    use ESMF                  , only : ESMF_Calendar
    use ESMF                  , only : ESMF_Time, ESMF_TimeGet
    use ESMF                  , only : ESMF_TimeInterval, ESMF_TimeIntervalGet
    use ESMF                  , only : ESMF_Alarm, ESMF_AlarmIsRinging, ESMF_AlarmRingerOff
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LOGMSG_ERROR
    use ESMF                  , only : ESMF_SUCCESS, ESMF_FAILURE
    use ESMF                  , only : operator(==), operator(-)
    use ESMF                  , only : ESMF_FieldBundleIsCreated, ESMF_MAXSTR, ESMF_ClockPrint, ESMF_AlarmIsCreated
    use NUOPC                 , only : NUOPC_CompAttributeGet
    use esmFlds               , only : compatm, complnd, compocn, compice, comprof, compglc, ncomps, compname
    use esmFlds               , only : fldListFr, fldListTo
    use shr_nuopc_utils_mod   , only : chkerr          => shr_nuopc_utils_ChkErr
    use shr_nuopc_methods_mod , only : FB_reset        => shr_nuopc_methods_FB_reset
    use shr_nuopc_methods_mod , only : FB_diagnose     => shr_nuopc_methods_FB_diagnose
    use shr_nuopc_methods_mod , only : FB_GetFldPtr    => shr_nuopc_methods_FB_GetFldPtr
    use shr_nuopc_methods_mod , only : FB_accum        => shr_nuopc_methods_FB_accum
    use shr_nuopc_methods_mod , only : State_GetScalar => shr_nuopc_methods_State_GetScalar
    use shr_nuopc_time_mod    , only : alarmInit       => shr_nuopc_time_alarmInit
    use med_constants_mod     , only : dbug_flag       => med_constants_dbug_flag
    use med_constants_mod     , only : SecPerDay       => med_constants_SecPerDay
    use med_map_mod           , only : med_map_FB_Regrid_Norm
    use med_internalstate_mod , only : InternalState, mastertask
    use med_io_mod            , only : med_io_write, med_io_wopen, med_io_enddef
    use med_io_mod            , only : med_io_close, med_io_date2yyyymmdd, med_io_sec2hms
    use med_io_mod            , only : med_io_ymd2date
    use perf_mod              , only : t_startf, t_stopf

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_VM)           :: vm
    type(ESMF_Clock)        :: clock
    type(ESMF_Time)         :: currtime
    type(ESMF_Time)         :: reftime
    type(ESMF_Time)         :: starttime
    type(ESMF_Time)         :: nexttime
    type(ESMF_TimeInterval) :: timediff       ! Used to calculate curr_time
    type(ESMF_Calendar)     :: calendar       ! calendar type
    character(len=64)       :: currtimestr
    character(len=64)       :: nexttimestr
    type(InternalState)     :: is_local
    character(CS)           :: histavg_option ! Histavg option units
    integer                 :: i,j,m,n,n1,ncnt
    integer                 :: start_ymd      ! Starting date YYYYMMDD
    integer                 :: start_tod      ! Starting time-of-day (s)
    integer                 :: nx,ny          ! global grid size
    integer                 :: yr,mon,day,sec ! time units
    real(r8)                :: rval           ! real tmp value
    real(r8)                :: dayssince      ! Time interval since reference time
    integer                 :: fk             ! index
    character(CL)           :: time_units     ! units of time variable
    character(CL)           :: case_name      ! case name
    character(CL)           :: hist_file      ! Local path to history filename
    character(CS)           :: cpl_inst_tag   ! instance tag
    character(CL)           :: cvalue         ! attribute string
    character(CL)           :: freq_option    ! freq_option setting (ndays, nsteps, etc)
    integer                 :: freq_n         ! freq_n setting relative to freq_option
    logical                 :: alarmIsOn      ! generic alarm flag
    real(r8)                :: tbnds(2)       ! CF1.0 time bounds
    logical                 :: whead,wdata    ! for writing restart/history cdf files
    integer                 :: dbrc
    integer                 :: iam
    logical                 :: isPresent
    logical,save            :: first_call = .true.
    character(len=*), parameter :: subname='(med_phases_history_write)'
    !---------------------------------------

    call t_startf('MED:'//subname)
    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! --- Get the communicator and localpet
    !---------------------------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localPet=iam, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! --- Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=case_name, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name='inst_suffix', isPresent=isPresent, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if(isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='inst_suffix', value=cpl_inst_tag, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       cpl_inst_tag = ""
    endif

    !---------------------------------------
    ! --- Get the clock info
    !---------------------------------------

    call ESMF_GridCompGet(gcomp, clock=clock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(clock, currtime=currtime, reftime=reftime, starttime=starttime, calendar=calendar, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGetNextTime(clock, nextTime=nexttime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet(currtime,yy=yr, mm=mon, dd=day, s=sec, rc=dbrc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    write(currtimestr,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',sec
    if (dbug_flag > 1) then
       call ESMF_LogWrite(trim(subname)//": currtime = "//trim(currtimestr), ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    call ESMF_TimeGet(nexttime,yy=yr, mm=mon, dd=day, s=sec, rc=dbrc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    write(nexttimestr,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',sec
    if (dbug_flag > 1) then
       call ESMF_LogWrite(trim(subname)//": nexttime = "//trim(nexttimestr), ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    timediff = nexttime - reftime
    call ESMF_TimeIntervalGet(timediff, d=day, s=sec, rc=rc)
    dayssince = day + sec/real(SecPerDay,R8)

    call ESMF_TimeGet(reftime, yy=yr, mm=mon, dd=day, s=sec, rc=dbrc)
    call med_io_ymd2date(yr,mon,day,start_ymd)
    start_tod = sec
    time_units = 'days since ' // trim(med_io_date2yyyymmdd(start_ymd)) // ' ' // med_io_sec2hms(start_tod, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! --- History Alarms
    !---------------------------------------

    if (.not. ESMF_AlarmIsCreated(AlarmHist, rc=rc)) then
       ! Set instantaneous history output alarm
       call NUOPC_CompAttributeGet(gcomp, name='history_option', value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       freq_option = cvalue

       call NUOPC_CompAttributeGet(gcomp, name='history_n', value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) freq_n

       call ESMF_LogWrite(trim(subname)//" init history alarm with option, n = "//&
            trim(freq_option)//","//trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)

       call alarmInit(clock, AlarmHist, option=freq_option, opt_n=freq_n, &
            RefTime=RefTime, alarmname='history', rc=rc)
    endif

    if (ESMF_AlarmIsRinging(AlarmHist, rc=rc)) then
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       alarmIsOn = .true.
       call ESMF_AlarmRingerOff( AlarmHist, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
#if DEBUG
       if (mastertask) then
          call ESMF_ClockPrint(clock, options="currTime", preString="-------->"//trim(subname)//&
               " history alarm for: ", rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
#endif
    else
       alarmisOn = .false.
    endif

    ! Set average history output alarm TODO: fix the following
    ! if (.not. ESMF_AlarmIsCreated(AlarmHistAvg, rc=rc)) then
    !    call NUOPC_CompAttributeGet(gcomp, name="histavg_option", value=histavg_option, rc=rc)
    !    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !    freq_option = cvalue
    !    call NUOPC_CompAttributeGet(gcomp, name="histavg_n", value=cvalue, rc=rc)
    !    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !    read(cvalue,*) freq_n
    !    call alarmInit(clock, AlarmHistAvg, option=freq_option, opt_n=freq_n, &
    !         RefTime=RefTime, alarmname='history_avg', rc=rc)
    ! end if
    ! if (ESMF_AlarmIsRinging(AlarmHistAvg, rc=rc)) then
    !    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    !    alarmIsOn = .true.
    !    call ESMF_AlarmRingerOff( AlarmHist, rc=rc )
    !    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! else
    !    alarmisOn = .false.
    ! endif

    !---------------------------------------
    ! --- History File
    ! Use nexttimestr rather than currtimestr here since that is the time at the end of
    ! the timestep and is preferred for history file names
    !---------------------------------------

    if (alarmIsOn) then
       write(hist_file,"(6a)") &
            trim(case_name), '.cpl',trim(cpl_inst_tag),'.hi.', trim(nexttimestr),'.nc'
       call ESMF_LogWrite(trim(subname)//": write "//trim(hist_file), ESMF_LOGMSG_INFO, rc=dbrc)
       call med_io_wopen(hist_file, vm, iam, clobber=.true.)

       do m = 1,2
          whead=.false.
          wdata=.false.
          if (m == 1) then
             whead=.true.
          elseif (m == 2) then
             wdata=.true.
             call med_io_enddef(hist_file)
          endif

          tbnds = dayssince

          call ESMF_LogWrite(trim(subname)//": time "//trim(time_units), ESMF_LOGMSG_INFO, rc=dbrc)
          if (tbnds(1) >= tbnds(2)) then
             call med_io_write(hist_file, iam, &
                  time_units=time_units, calendar=calendar, time_val=dayssince, &
                  whead=whead, wdata=wdata, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          else
             call med_io_write(hist_file, iam, &
                  time_units=time_units, calendar=calendar, time_val=dayssince, &
                  whead=whead, wdata=wdata, tbnds=tbnds, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          endif

          do n = 1,ncomps
             if (is_local%wrap%comp_present(n)) then
                if (ESMF_FieldBundleIsCreated(is_local%wrap%FBimp(n,n),rc=rc)) then
                   nx = is_local%wrap%nx(n)
                   ny = is_local%wrap%ny(n)
                   call med_io_write(hist_file, iam, is_local%wrap%FBimp(n,n), &
                       nx=nx, ny=ny, nt=1, whead=whead, wdata=wdata, pre=trim(compname(n))//'Imp', rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                endif
                if (ESMF_FieldBundleIsCreated(is_local%wrap%FBexp(n),rc=rc)) then
                   nx = is_local%wrap%nx(n)
                   ny = is_local%wrap%ny(n)
                   call med_io_write(hist_file, iam, is_local%wrap%FBexp(n), &
                       nx=nx, ny=ny, nt=1, whead=whead, wdata=wdata, pre=trim(compname(n))//'Exp', rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                endif
             endif
          enddo

       enddo

       call med_io_close(hist_file, iam, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    endif

    !---------------------------------------
    !--- clean up
    !---------------------------------------

    first_call = .false.

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    call t_stopf('MED:'//subname)

  end subroutine med_phases_history_write

  !===============================================================================

end module med_phases_history_mod
