module med_phases_history_mod

  !-----------------------------------------------------------------------------
  ! Mediator Phases
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use shr_kind_mod            , only : IN=>SHR_KIND_IN, R8=>SHR_KIND_R8
  use shr_kind_mod            , only : CL=>SHR_KIND_CL, CS=>SHR_KIND_CS
  use shr_cal_mod             , only : shr_cal_noleap, shr_cal_gregorian
  use shr_cal_mod             , only : shr_cal_ymd2date
  use seq_timemgr_mod         , only : seq_timemgr_AlarmInit, seq_timemgr_AlarmIsOn
  use seq_timemgr_mod         , only : seq_timemgr_AlarmSetOff
  use esmFlds                 , only : compatm, complnd, compocn, compice, comprof, compglc
  use esmFlds                 , only : ncomps, compname 
  use esmFlds                 , only : fldListFr, fldListTo
  use esmFlds                 , only : fldListMed_aoflux_a, fldListMed_aoflux_o
  use esmFlds                 , only : flds_scalar_name, flds_scalar_num
  use esmFlds                 , only : flds_scalar_index_nx, flds_scalar_index_ny
  use shr_nuopc_methods_mod   , only : shr_nuopc_methods_ChkErr
  use shr_nuopc_methods_mod   , only : shr_nuopc_methods_FB_reset
  use shr_nuopc_methods_mod   , only : shr_nuopc_methods_FB_diagnose
  use shr_nuopc_methods_mod   , only : shr_nuopc_methods_FB_GetFldPtr
  use shr_nuopc_methods_mod   , only : shr_nuopc_methods_FB_accum
  use shr_nuopc_methods_mod   , only : shr_nuopc_methods_State_GetScalar
  use med_constants_mod       , only : med_constants_dbug_flag
  use med_constants_mod       , only : med_constants_czero
  use med_infodata_mod        , only : med_infodata, med_infodata_GetData
  use med_merge_mod           , only : med_merge_auto
  use med_map_mod             , only : med_map_FB_Regrid_Norm 
  use med_internalstate_mod   , only : InternalState
  use med_io_mod              , only : med_io_write, med_io_wopen, med_io_enddef
  use med_io_mod              , only : med_io_close, med_io_date2yyyymmdd
  use med_io_mod              , only : med_io_sec2hms

  implicit none
  private

  integer           , parameter :: dbug_flag = med_constants_dbug_flag
  real(ESMF_KIND_R8), parameter :: czero     = med_constants_czero
  character(*)      , parameter :: u_FILE_u  = __FILE__
  character(len=ESMF_MAXSTR)    :: tmpstr
  integer                       :: dbrc
  logical                       :: mastertask
  integer, parameter            :: SecPerDay = 86400    ! Seconds per day
  type(ESMF_Alarm)              :: AlarmHist

  public  :: med_phases_history

!===============================================================================
contains
!===============================================================================

  subroutine med_phases_history(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Write mediator history file

    ! local variables
    type(InternalState)     :: is_local
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
    integer                 :: i,j,m,n,n1,ncnt
    logical,save            :: first_call = .true.
    integer(IN)             :: start_ymd      ! Starting date YYYYMMDD
    integer(IN)             :: start_tod      ! Starting time-of-day (s)
    integer(IN)             :: nx,ny          ! global grid size
    integer(IN)             :: yr,mon,day,sec ! time units
    real(r8)                :: dayssince      ! Time interval since reference time
    character(CL)           :: time_units     ! units of time variable
    character(CL)           :: calendar       ! calendar type
    character(CL)           :: case_name      ! case name
    character(CL)           :: hist_file      ! Local path to history filename
    character(CS)           :: cpl_inst_tag   ! instance tag
    character(CL)           :: cvalue         ! attribute string
    character(CL)           :: freq_option    ! freq_option setting (ndays, nsteps, etc)
    integer(IN)             :: freq_n         ! freq_n setting relative to freq_option
    logical                 :: alarmIsOn      ! generic alarm flag
    real(r8)                :: tbnds(2)       ! CF1.0 time bounds
    logical                 :: whead,wdata    ! for writing restart/history cdf files
    integer                 :: mpicom, iam 
    character(len=*), parameter :: subname='(med_phases_history)'
    !---------------------------------------

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! --- Get the communicator and localpet
    !---------------------------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, mpiCommunicator=mpicom, localPet=iam, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! --- Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=case_name, rc=rc)
    cpl_inst_tag = ''

    !---------------------------------------
    ! --- Get the clock info
    !---------------------------------------

    call ESMF_GridCompGet(gcomp, clock=clock)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(clock, currtime=currtime, reftime=reftime, starttime=starttime, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGetNextTime(clock, nextTime=nexttime, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(clock, calkindflag=calkindflag, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (calkindflag == ESMF_CALKIND_GREGORIAN) then
      calendar = shr_cal_gregorian
    elseif (calkindflag == ESMF_CALKIND_NOLEAP) then
      calendar = shr_cal_noleap
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

    call ESMF_ClockPrint(clock, options="currTime", preString="-------->"//trim(subname)//" mediating for: ", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    timediff = currtime - reftime
    call ESMF_TimeIntervalGet(timediff, d=day, s=sec, rc=rc)
    dayssince = day + sec/real(SecPerDay,R8)

    call ESMF_TimeGet(reftime, yy=yr, mm=mon, dd=day, s=sec, rc=dbrc)
    call shr_cal_ymd2date(yr,mon,day,start_ymd)
    start_tod = sec
    time_units = 'days since ' &
         // trim(med_io_date2yyyymmdd(start_ymd)) // ' ' // med_io_sec2hms(start_tod)

    !---------------------------------------
    ! --- History Alarm
    !---------------------------------------

    if (.not. ESMF_AlarmIsCreated(AlarmHist, rc=rc)) then
       call NUOPC_CompAttributeGet(gcomp, name='history_option', value=cvalue, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       freq_option = cvalue
       call NUOPC_CompAttributeGet(gcomp, name='history_n', value=cvalue, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) freq_n
       call ESMF_LogWrite(trim(subname)//" init history alarm with option, n = "//trim(freq_option)//","//trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)
       call seq_timemgr_alarmInit(clock, AlarmHist, option=freq_option, opt_n=freq_n, RefTime=RefTime, alarmname='history', rc=rc)
    endif

    alarmIsOn = seq_timemgr_alarmIsOn(clock, 'history', rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call seq_timemgr_AlarmSetOff(clock, 'history', rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! --- History File
    ! Use nexttimestr rather than currtimestr here since that is the time at the end of 
    ! the timestep and is preferred for history file names
    !---------------------------------------

    if (alarmIsOn) then

       write(hist_file,"(6a)") &
            !trim(case_name), '.cpl',trim(cpl_inst_tag),'.hi.', trim(currtimestr),'.nc'
            trim(case_name), '.cpl',trim(cpl_inst_tag),'.hi.', trim(nexttimestr),'.nc'

       call ESMF_LogWrite(trim(subname)//": write "//trim(hist_file), ESMF_LOGMSG_INFO, rc=dbrc)
       call med_io_wopen(hist_file, mpicom, iam, clobber=.true.)

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
          !------- tcx nov 2011 tbnds of same values causes problems in ferret                  
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
                   !write(tmpstr,*) subname,' nx,ny = ',trim(compname(n)),nx,ny
                   !call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

                   call med_io_write(hist_file, iam, is_local%wrap%FBimp(n,n), &
                       nx=nx, ny=ny, nt=1, whead=whead, wdata=wdata, pre=trim(compname(n))//'Imp', rc=rc)
                   if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                endif
                if (ESMF_FieldBundleIsCreated(is_local%wrap%FBexp(n),rc=rc)) then
                   call med_infodata_GetData(med_infodata, ncomp=n, nx=nx, ny=ny)
                   !write(tmpstr,*) subname,' nx,ny = ',trim(compname(n)),nx,ny
                   !call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)

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

  end subroutine med_phases_history

end module med_phases_history_mod
