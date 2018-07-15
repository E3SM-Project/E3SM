module med_phases_restart_mod

  !-----------------------------------------------------------------------------
  ! Write/Read mediator restart files
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use shr_file_mod            , only : shr_file_getUnit, shr_file_freeUnit
  use shr_mpi_mod             , only : shr_mpi_bcast
  use esmFlds                 , only : ncomps, compname, compocn 
  use shr_nuopc_methods_mod   , only : shr_nuopc_methods_ChkErr
  use shr_nuopc_methods_mod   , only : shr_nuopc_methods_FB_reset
  use shr_nuopc_methods_mod   , only : shr_nuopc_methods_FB_diagnose
  use shr_nuopc_methods_mod   , only : shr_nuopc_methods_FB_GetFldPtr
  use shr_nuopc_methods_mod   , only : shr_nuopc_methods_State_GetScalar
  use med_constants_mod       , only : med_constants_dbug_flag
  use med_constants_mod       , only : med_constants_czero
  use med_infodata_mod        , only : med_infodata, med_infodata_GetData
  use med_io_mod              , only : med_io_write, med_io_wopen, med_io_enddef
  use med_io_mod              , only : med_io_close, med_io_date2yyyymmdd
  use med_io_mod              , only : med_io_sec2hms, med_io_read
  use med_internalstate_mod   , only : InternalState
  !
  use seq_timemgr_mod         , only : seq_timemgr_AlarmIsOn
  use seq_timemgr_mod         , only : seq_timemgr_AlarmSetOff, seq_timemgr_alarm_restart

  implicit none
  private

  integer            , parameter :: dbug_flag = med_constants_dbug_flag
  real(ESMF_KIND_R8) , parameter :: czero     = med_constants_czero
  integer            , parameter :: SecPerDay = 86400    ! Seconds per day
  character(len=*)   , parameter :: sp_str = 'str_undefined'
  character(len=*)   , parameter :: shr_cal_noleap    = 'NO_LEAP'
  character(len=*)   , parameter :: shr_cal_gregorian = 'GREGORIAN'
  integer                        :: dbrc
  character(*)       , parameter :: u_FILE_u  = &
       __FILE__

  public  :: med_phases_restart_read
  public  :: med_phases_restart_write

!=================================================================================
contains
!=================================================================================

  subroutine med_phases_restart_write(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Carry out fast accumulation for the ocean

    ! local variables
    type(ESMF_VM)              :: vm
    type(ESMF_Clock)           :: clock
    type(ESMF_Time)            :: currtime, reftime, starttime, nexttime
    type(ESMF_TimeInterval)    :: timediff       ! Used to calculate curr_time
    type(ESMF_CalKind_Flag)    :: calkindflag
    character(len=64)          :: currtimestr, nexttimestr
    type(InternalState)        :: is_local
    integer                    :: i,j,m,n,n1,ncnt
    integer                    :: curr_ymd       ! Current date YYYYMMDD
    integer                    :: curr_tod       ! Current time-of-day (s)
    integer                    :: start_ymd      ! Starting date YYYYMMDD
    integer                    :: start_tod      ! Starting time-of-day (s)
    integer                    :: ref_ymd        ! Reference date YYYYMMDD
    integer                    :: ref_tod        ! Reference time-of-day (s)
    integer                    :: next_ymd       ! Starting date YYYYMMDD
    integer                    :: next_tod       ! Starting time-of-day (s)
    integer                    :: nx,ny          ! global grid size
    integer                    :: yr,mon,day,sec ! time units
    real(ESMF_KIND_R8)         :: dayssince      ! Time interval since reference time
    integer                    :: unitn          ! unit number
    character(ESMF_MAXSTR)     :: time_units     ! units of time variable
    character(ESMF_MAXSTR)     :: calendar       ! calendar type
    character(ESMF_MAXSTR)     :: case_name      ! case name
    character(ESMF_MAXSTR)     :: restart_file   ! Local path to restart filename
    character(ESMF_MAXSTR)     :: restart_pfile  ! Local path to restart pointer filename
    character(ESMF_MAXSTR)     :: cpl_inst_tag   ! instance tag
    character(ESMF_MAXSTR)     :: cvalue         ! attribute string
    character(ESMF_MAXSTR)     :: freq_option    ! freq_option setting (ndays, nsteps, etc)
    integer                    :: freq_n         ! freq_n setting relative to freq_option
    logical                    :: alarmIsOn      ! generic alarm flag
    real(ESMF_KIND_R8)         :: tbnds(2)       ! CF1.0 time bounds
    logical                    :: whead,wdata    ! for writing restart/restart cdf files
    integer                    :: iam,mpicom     ! vm stuff
    character(len=ESMF_MAXSTR) :: tmpstr
    character(len=*), parameter :: subname='(med_phases_restart_write)'
    !---------------------------------------

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! --- Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, mpiCommunicator=mpicom, localPet=iam, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=case_name, rc=rc)
    cpl_inst_tag = ''

    !---------------------------------------
    ! --- Get the clock info
    !---------------------------------------

    call ESMF_GridCompGet(gcomp, clock=clock)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! --- Restart Alarm
    !---------------------------------------

    alarmIsOn = seq_timemgr_alarmIsOn(clock, seq_timemgr_alarm_restart, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call seq_timemgr_AlarmSetOff(clock, seq_timemgr_alarm_restart, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (alarmIsOn) then

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

       timediff = nexttime - reftime
       call ESMF_TimeIntervalGet(timediff, d=day, s=sec, rc=rc)
       dayssince = day + sec/real(SecPerDay,ESMF_KIND_R8)

       call ESMF_TimeGet(reftime, yy=yr, mm=mon, dd=day, s=sec, rc=dbrc)
       call ymd2date(yr,mon,day,start_ymd)
       start_tod = sec
       time_units = 'days since ' &
          // trim(med_io_date2yyyymmdd(start_ymd)) // ' ' // med_io_sec2hms(start_tod)

       call ESMF_TimeGet(nexttime, yy=yr, mm=mon, dd=day, s=sec, rc=dbrc)
       call ymd2date(yr,mon,day,next_ymd)
       next_tod = sec

       call ESMF_TimeGet(reftime, yy=yr, mm=mon, dd=day, s=sec, rc=dbrc)
       call ymd2date(yr,mon,day,ref_ymd)
       ref_tod = sec

       call ESMF_TimeGet(currtime, yy=yr, mm=mon, dd=day, s=sec, rc=dbrc)
       call ymd2date(yr,mon,day,curr_ymd)
       curr_tod = sec

       !---------------------------------------
       ! --- Restart File
       ! Use nexttimestr rather than currtimestr here since that is the time at the end of 
       ! the timestep and is preferred for restart file names
       !---------------------------------------

       write(restart_file,"(6a)") &
           !trim(case_name), '.cpl',trim(cpl_inst_tag),'.r.', trim(currtimestr),'.nc'
            trim(case_name), '.cpl',trim(cpl_inst_tag),'.r.', trim(nexttimestr),'.nc'

       if (iam == 0) then
          call NUOPC_CompAttributeGet(gcomp, name='restart_pfile', value=restart_pfile, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          call ESMF_LogWrite(trim(subname)//" write rpointer file = "//trim(restart_pfile), ESMF_LOGMSG_INFO, rc=dbrc)
          unitn = shr_file_getUnit()
          open(unitn, file=restart_pfile, form='FORMATTED')
          write(unitn,'(a)') trim(restart_file)
          close(unitn)
          call shr_file_freeUnit( unitn )
       endif

       call ESMF_LogWrite(trim(subname)//": write "//trim(restart_file), ESMF_LOGMSG_INFO, rc=dbrc)
       call med_io_wopen(restart_file, mpicom, iam, clobber=.true.)

       do m = 1,2
          if (m == 1) then
             whead = .true.
             wdata = .false.
          else if (m == 2) then
             whead = .false.
             wdata = .true.
          endif
          if (wdata) then
             call med_io_enddef(restart_file)
          end if

          tbnds = dayssince
          call ESMF_LogWrite(trim(subname)//": time "//trim(time_units), ESMF_LOGMSG_INFO, rc=dbrc)
          if (tbnds(1) >= tbnds(2)) then
             call med_io_write(restart_file, iam=iam, &
                  time_units=time_units, time_cal=calendar, time_val=dayssince, &
                  whead=whead, wdata=wdata)
          else
             call med_io_write(restart_file, iam=iam, &
                  time_units=time_units, time_cal=calendar, time_val=dayssince, &
                  whead=whead, wdata=wdata, tbnds=tbnds)
          endif

          ! Write out next ymd/tod in place of curr ymd/tod because
          ! the currently the restart represents the time at end of
          ! the current timestep and that is where we want to start
          ! the next run.

          call med_io_write(restart_file, iam, start_ymd, 'seq_timemgr_start_ymd', whead=whead, wdata=wdata)
          call med_io_write(restart_file, iam, start_tod, 'seq_timemgr_start_tod', whead=whead, wdata=wdata)
          call med_io_write(restart_file, iam, ref_ymd  , 'seq_timemgr_ref_ymd'  , whead=whead, wdata=wdata)
          call med_io_write(restart_file, iam, ref_tod  , 'seq_timemgr_ref_tod'  , whead=whead, wdata=wdata)
          call med_io_write(restart_file, iam, next_ymd , 'seq_timemgr_curr_ymd' , whead=whead, wdata=wdata)
          call med_io_write(restart_file, iam, next_tod , 'seq_timemgr_curr_tod' , whead=whead, wdata=wdata)

          call med_io_write(restart_file, iam, is_local%wrap%FBExpAccumCnt, dname='ExpAccumCnt', &
               whead=whead, wdata=wdata)
          !if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
 
          do n = 1,ncomps
             if (is_local%wrap%comp_present(n)) then
                ! Write import field bundles
                if (ESMF_FieldBundleIsCreated(is_local%wrap%FBimp(n,n),rc=rc)) then
                   call med_infodata_GetData(med_infodata, ncomp=n, nx=nx, ny=ny)
                   !write(tmpstr,*) subname,' nx,ny = ',trim(compname(n)),nx,ny
                   !call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
                   call med_io_write(restart_file, iam, is_local%wrap%FBimp(n,n), &
                       nx=nx, ny=ny, nt=1, whead=whead, wdata=wdata, pre=trim(compname(n))//'Imp', rc=rc)
                   if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                endif

                ! Write fraction field bundles
                if (ESMF_FieldBundleIsCreated(is_local%wrap%FBfrac(n),rc=rc)) then
                   call med_infodata_GetData(med_infodata, ncomp=n, nx=nx, ny=ny)
                   !write(tmpstr,*) subname,' nx,ny = ',trim(compname(n)),nx,ny
                   !call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
                   call med_io_write(restart_file, iam, is_local%wrap%FBfrac(n), &
                       nx=nx, ny=ny, nt=1, whead=whead, wdata=wdata, pre=trim(compname(n))//'Frac', rc=rc)
                   if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                endif

                ! Write export accumulators 
                if (ESMF_FieldBundleIsCreated(is_local%wrap%FBExpAccum(n),rc=rc)) then
                   ! TODO: only write this out if actually have done accumulation
                   call med_infodata_GetData(med_infodata, ncomp=n, nx=nx, ny=ny)
                   !write(tmpstr,*) subname,' nx,ny = ',trim(compname(n)),nx,ny
                   !call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
                   call med_io_write(restart_file, iam,  is_local%wrap%FBExpAccum(n), &
                       nx=nx, ny=ny, nt=1, whead=whead, wdata=wdata, pre=trim(compname(n))//'ExpAccum', rc=rc)
                   if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                endif
             endif
          enddo

          ! Write ocn albedo field bundle (CESM only)
          ! if (ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_ocnalb_o,rc=rc)) then
          !    call med_infodata_GetData(med_infodata, ncomp=compocn, nx=nx, ny=ny)
          !    !write(tmpstr,*) subname,' nx,ny = ',trim(compname(n)),nx,ny
          !    !call ESMF_LogWrite(trim(tmpstr), ESMF_LOGMSG_INFO, rc=dbrc)
          !    call med_io_write(restart_file, iam, is_local%wrap%FBMed_ocnalb_o, &
          !         nx=nx, ny=ny, nt=1, whead=whead, wdata=wdata, pre='MedOcnAlb_o', rc=rc)
          !    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          ! end if

       enddo

       ! Close file
       call med_io_close(restart_file, iam)
    endif

    !---------------------------------------
    !--- clean up
    !---------------------------------------

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_phases_restart_write

  !===============================================================================
  subroutine med_phases_restart_read(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Carry out fast accumulation for the ocean

    ! local variables
    type(ESMF_VM)          :: vm
    type(ESMF_Clock)       :: clock
    type(ESMF_Time)        :: currtime
    character(len=64)      :: currtimestr
    type(InternalState)    :: is_local
    integer                :: i,j,m,n,n1,ncnt
    integer                :: ierr, unitn
    integer                :: yr,mon,day,sec ! time units
    integer                :: iam,mpicom     ! vm stuff
    character(ESMF_MAXSTR) :: case_name      ! case name
    character(ESMF_MAXSTR) :: restart_file   ! Local path to restart filename
    character(ESMF_MAXSTR) :: restart_pfile  ! Local path to restart pointer filename
    character(ESMF_MAXSTR) :: cpl_inst_tag   ! instance tag
    character(len=*), parameter :: subname='(med_phases_restart_read)'
    !---------------------------------------

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! --- Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, localPet=iam, mpiCommunicator=mpicom, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=case_name, rc=rc)
    cpl_inst_tag = ''

    !---------------------------------------
    ! --- Get the clock info
    !---------------------------------------

    call ESMF_GridCompGet(gcomp, clock=clock)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(clock, currtime=currtime, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet(currtime,yy=yr, mm=mon, dd=day, s=sec, rc=dbrc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    write(currtimestr,'(i4.4,a,i2.2,a,i2.2,a,i5.5)') yr,'-',mon,'-',day,'-',sec
    if (dbug_flag > 1) then
       call ESMF_LogWrite(trim(subname)//": currtime = "//trim(currtimestr), ESMF_LOGMSG_INFO, rc=dbrc)
    endif

    call ESMF_ClockPrint(clock, options="currTime", preString="-------->"//trim(subname)//" mediating for: ", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! --- Restart File
    !---------------------------------------

    ! Error check on restart_pfile                                                                                                      
    call NUOPC_CompAttributeGet(gcomp, name="restart_pfile", value=restart_pfile, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if ( len_trim(restart_pfile) == 0 ) then
       call ESMF_LogWrite(trim(subname)//' ERROR restart_pfile must be set', &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__, rc=rc)
       rc = ESMF_FAILURE
       return
    end if

    if (iam == 0) then
       call NUOPC_CompAttributeGet(gcomp, name='restart_file', value=restart_file, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       !--- read rpointer if restart_file is set to sp_str ---                                                                       

       if (trim(restart_file) == trim(sp_str)) then
          call NUOPC_CompAttributeGet(gcomp, name='restart_pfile', value=restart_file, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          unitn = shr_file_getUnit()
          call ESMF_LogWrite(trim(subname)//" read rpointer file = "//trim(restart_pfile), ESMF_LOGMSG_INFO, rc=dbrc)
          open(unitn, file=restart_pfile, form='FORMATTED', status='old',iostat=ierr)
          if (ierr < 0) then
             call ESMF_LogWrite(trim(subname)//' rpointer file open returns error', ESMF_LOGMSG_INFO, rc=dbrc)
             rc=ESMF_Failure
             return
          end if
          read(unitn,'(a)', iostat=ierr) restart_file
          if (ierr < 0) then
             call ESMF_LogWrite(trim(subname)//' rpointer file read returns error', ESMF_LOGMSG_INFO, rc=dbrc)
             rc=ESMF_Failure
             return
          end if
          close(unitn)
          call shr_file_freeUnit( unitn )
          call ESMF_LogWrite(trim(subname)//' restart file from rpointer = '//trim(restart_file), &
               ESMF_LOGMSG_INFO, rc=dbrc)
       endif
    endif
    call shr_mpi_bcast(restart_file, mpicom)

    call ESMF_LogWrite(trim(subname)//": read "//trim(restart_file), ESMF_LOGMSG_INFO, rc=dbrc)
 
    call med_io_read(restart_file, mpicom, iam, is_local%wrap%FBExpAccumCnt, dname='ExpAccumCnt')

    do n = 1,ncomps
       if (is_local%wrap%comp_present(n)) then
          ! Read import field bundle
          if (ESMF_FieldBundleIsCreated(is_local%wrap%FBimp(n,n),rc=rc)) then
             call med_io_read(restart_file, mpicom, iam, is_local%wrap%FBimp(n,n), &
                 pre=trim(compname(n))//'Imp', rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          endif

          ! Read import fractions
          if (ESMF_FieldBundleIsCreated(is_local%wrap%FBfrac(n),rc=rc)) then
             call med_io_read(restart_file, mpicom, iam, is_local%wrap%FBfrac(n), &
                 pre=trim(compname(n))//'Frac', rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          endif

          ! Read export field bundle accumulator
          if (ESMF_FieldBundleIsCreated(is_local%wrap%FBExpAccum(n),rc=rc)) then
             call med_io_read(restart_file, mpicom, iam, is_local%wrap%FBExpAccum(n), &
                 pre=trim(compname(n))//'ExpAccum', rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          endif
       endif
    enddo

    ! read ocn albedo field bundle (CESM only)
    ! if (ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_ocnalb_o,rc=rc)) then
    !    call med_io_read(restart_file, mpicom, iam, is_local%wrap%FBMed_ocnalb_o, &
    !         pre='MedOcnAlb_o', rc=rc) 
    !    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    ! end if

    !---------------------------------------
    !--- clean up
    !---------------------------------------

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_phases_restart_read


  !===============================================================================
  subroutine ymd2date(year,month,day,date)
    ! Converts  year, month, day to coded-date
    ! NOTE: this calendar has a year zero (but no day or month zero)

    integer,intent(in ) :: year,month,day  ! calendar year,month,day
    integer,intent(out) :: date            ! coded (yyyymmdd) calendar date

    ! local variables
    character(*),parameter :: subName = "(ymd2date)"
    !---------------------------------------

    date = abs(year)*10000 + month*100 + day  ! coded calendar date
    if (year < 0) date = -date
  end subroutine ymd2date

end module med_phases_restart_mod
