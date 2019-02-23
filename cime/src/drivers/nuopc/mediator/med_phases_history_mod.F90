module med_phases_history_mod

  !-----------------------------------------------------------------------------
  ! Mediator Phases
  !-----------------------------------------------------------------------------

  use ESMF, only : ESMF_Alarm

  implicit none
  private

  character(*)      , parameter :: u_FILE_u  = __FILE__
  type(ESMF_Alarm)              :: AlarmHist
  type(ESMF_Alarm)              :: AlarmHistAvg

  public  :: med_phases_history_write

!===============================================================================
contains
!===============================================================================

  subroutine med_phases_history_write(gcomp, rc)

    ! Write mediator history file

    use ESMF                  , only : ESMF_GridComp, ESMF_VM, ESMF_Clock, ESMF_Time, ESMF_TimeInterval, ESMF_CalKind_Flag
    use ESMF                  , only : ESMF_GridCompGet, ESMF_ClockGet, ESMF_ClockGetNextTime, ESMF_VMGet, ESMF_TimeGet
    use ESMF                  , only : ESMF_TimeIntervalGet, ESMF_AlarmIsRinging, ESMF_AlarmRingerOff
    use ESMF                  , only : ESMF_CALKIND_GREGORIAN, ESMF_CALKIND_NOLEAP
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LOGMSG_ERROR, ESMF_SUCCESS, ESMF_FAILURE
    use ESMF                  , only : operator(==), operator(-)
    use ESMF                  , only : ESMF_FieldBundleIsCreated, ESMF_MAXSTR, ESMF_ClockPrint, ESMF_AlarmIsCreated
    use NUOPC                 , only : NUOPC_CompAttributeGet
    use shr_cal_mod           , only : shr_cal_ymd2date
    use esmFlds               , only : compatm, complnd, compocn, compice, comprof, compglc, ncomps, compname
    use esmFlds               , only : fldListFr, fldListTo
    use esmFlds               , only : fldListMed_aoflux_a, fldListMed_aoflux_o
    use shr_nuopc_scalars_mod , only : flds_scalar_index_nx, flds_scalar_index_ny
    use shr_nuopc_scalars_mod , only : flds_scalar_name, flds_scalar_num
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_reset
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_diagnose
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_GetFldPtr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_accum
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_GetScalar
    use shr_nuopc_time_mod    , only : shr_nuopc_time_alarmInit
    use med_constants_mod     , only : dbug_flag =>med_constants_dbug_flag
    use med_constants_mod     , only : SecPerDay =>med_constants_SecPerDay
    use med_constants_mod     , only : R8, CL, CS, IN
    use med_constants_mod     , only : med_constants_noleap, med_constants_gregorian
    use med_infodata_mod      , only : med_infodata, med_infodata_GetData
    use med_map_mod           , only : med_map_FB_Regrid_Norm
    use med_internalstate_mod , only : InternalState, mastertask
    use med_io_mod            , only : med_io_write, med_io_wopen, med_io_enddef
    use med_io_mod            , only : med_io_close, med_io_date2yyyymmdd
    use med_io_mod            , only : med_io_sec2hms
    use perf_mod              , only : t_startf, t_stopf
    ! Input/output variables
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
    type(ESMF_CalKind_Flag) :: calkindflag
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
    character(CL)           :: calendar       ! calendar type
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
    logical,save            :: first_call = .true.
    character(len=*), parameter :: subname='(med_phases_history_write)'
    logical :: isPresent

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
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localPet=iam, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! --- Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=case_name, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name='inst_suffix', isPresent=isPresent, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if(isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name='inst_suffix', value=cpl_inst_tag, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       cpl_inst_tag = ""
    endif
    !---------------------------------------
    ! --- Get the clock info
    !---------------------------------------

    call ESMF_GridCompGet(gcomp, clock=clock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(clock, currtime=currtime, reftime=reftime, starttime=starttime, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGetNextTime(clock, nextTime=nexttime, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(clock, calkindflag=calkindflag, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (calkindflag == ESMF_CALKIND_GREGORIAN) then
      calendar = med_constants_gregorian
    elseif (calkindflag == ESMF_CALKIND_NOLEAP) then
      calendar = med_constants_noleap
    else
      call ESMF_LogWrite(trim(subname)//' ERROR: calendar not supported', ESMF_LOGMSG_ERROR, rc=dbrc)
      rc=ESMF_Failure
      return
    endif

    call ESMF_TimeGet(currtime,yy=yr, mm=mon, dd=day, s=sec, rc=dbrc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    write(currtimestr,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',sec
    if (dbug_flag > 1) then
       call ESMF_LogWrite(trim(subname)//": currtime = "//trim(currtimestr), ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    call ESMF_TimeGet(nexttime,yy=yr, mm=mon, dd=day, s=sec, rc=dbrc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    write(nexttimestr,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',sec
    if (dbug_flag > 1) then
       call ESMF_LogWrite(trim(subname)//": nexttime = "//trim(nexttimestr), ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    timediff = nexttime - reftime
    call ESMF_TimeIntervalGet(timediff, d=day, s=sec, rc=rc)
    dayssince = day + sec/real(SecPerDay,R8)

    call ESMF_TimeGet(reftime, yy=yr, mm=mon, dd=day, s=sec, rc=dbrc)
    call shr_cal_ymd2date(yr,mon,day,start_ymd)
    start_tod = sec
    time_units = 'days since ' &
         // trim(med_io_date2yyyymmdd(start_ymd)) // ' ' // med_io_sec2hms(start_tod)

    !---------------------------------------
    ! --- History Alarms
    !---------------------------------------

    if (.not. ESMF_AlarmIsCreated(AlarmHist, rc=rc)) then
       ! Set instantaneous history output alarm
       call NUOPC_CompAttributeGet(gcomp, name='history_option', value=cvalue, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       freq_option = cvalue

       call NUOPC_CompAttributeGet(gcomp, name='history_n', value=cvalue, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) freq_n

       call ESMF_LogWrite(trim(subname)//" init history alarm with option, n = "//&
            trim(freq_option)//","//trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)

       call shr_nuopc_time_alarmInit(clock, AlarmHist, option=freq_option, opt_n=freq_n, &
            RefTime=RefTime, alarmname='history', rc=rc)
    endif

    if (ESMF_AlarmIsRinging(AlarmHist, rc=rc)) then
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       alarmIsOn = .true.
       call ESMF_AlarmRingerOff( AlarmHist, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
#if DEBUG
       if (mastertask) then
          call ESMF_ClockPrint(clock, options="currTime", preString="-------->"//trim(subname)//&
               " history alarm for: ", rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
#endif
    else
       alarmisOn = .false.
    endif

    ! Set average history output alarm TODO: fix the following
    ! if (.not. ESMF_AlarmIsCreated(AlarmHistAvg, rc=rc)) then
    !    call NUOPC_CompAttributeGet(gcomp, name="histavg_option", value=histavg_option, rc=rc)
    !    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    !    freq_option = cvalue
    !    call NUOPC_CompAttributeGet(gcomp, name="histavg_n", value=cvalue, rc=rc)
    !    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    !    read(cvalue,*) freq_n
    !    call shr_nuopc_time_alarmInit(clock, AlarmHistAvg, option=freq_option, opt_n=freq_n, &
    !         RefTime=RefTime, alarmname='history_avg', rc=rc)
    ! end if
    ! if (ESMF_AlarmIsRinging(AlarmHistAvg, rc=rc)) then
    !    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    !    alarmIsOn = .true.
    !    call ESMF_AlarmRingerOff( AlarmHist, rc=rc )
    !    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
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
                  time_units=time_units, time_cal=calendar, time_val=dayssince, &
                  whead=whead, wdata=wdata)
          else
             call med_io_write(hist_file, iam, &
                  time_units=time_units, time_cal=calendar, time_val=dayssince, &
                  whead=whead, wdata=wdata, tbnds=tbnds)
          endif

          do n = 1,ncomps
             if (is_local%wrap%comp_present(n)) then
                if (ESMF_FieldBundleIsCreated(is_local%wrap%FBimp(n,n),rc=rc)) then
                   call med_infodata_GetData(med_infodata, ncomp=n, nx=nx, ny=ny)
                   call med_io_write(hist_file, iam, is_local%wrap%FBimp(n,n), &
                       nx=nx, ny=ny, nt=1, whead=whead, wdata=wdata, pre=trim(compname(n))//'Imp', rc=rc)
                   if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                endif
                if (ESMF_FieldBundleIsCreated(is_local%wrap%FBexp(n),rc=rc)) then
                   call med_infodata_GetData(med_infodata, ncomp=n, nx=nx, ny=ny)
                   call med_io_write(hist_file, iam, is_local%wrap%FBexp(n), &
                       nx=nx, ny=ny, nt=1, whead=whead, wdata=wdata, pre=trim(compname(n))//'Exp', rc=rc)
                   if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                endif
             endif
          enddo

       enddo

       call med_io_close(hist_file, iam)

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

end module med_phases_history_mod
