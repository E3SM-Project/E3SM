module ice_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for DICE
  !----------------------------------------------------------------------------

  use ESMF
  use NUOPC            , only : NUOPC_CompDerive, NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize
  use NUOPC            , only : NUOPC_CompAttributeGet, NUOPC_Advertise
  use NUOPC_Model      , only : model_routine_SS        => SetServices
  use NUOPC_Model      , only : model_label_Advance     => label_Advance
  use NUOPC_Model      , only : model_label_SetRunClock => label_SetRunClock
  use NUOPC_Model      , only : model_label_Finalize    => label_Finalize
  use NUOPC_Model      , only : NUOPC_ModelGet
  use shr_file_mod     , only : shr_file_getlogunit, shr_file_setlogunit
  use shr_kind_mod     , only : r8=>shr_kind_r8, i8=>shr_kind_i8, cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_const_mod    , only : shr_const_spval, shr_const_pi
  use shr_sys_mod      , only : shr_sys_abort
  use shr_cal_mod      , only : shr_cal_noleap, shr_cal_gregorian, shr_cal_ymd2date, shr_cal_ymd2julian
  use shr_mpi_mod      , only : shr_mpi_bcast
  use dshr_strdata_mod , only : shr_strdata_type, shr_strdata_readnml
  use dshr_methods_mod , only : chkerr, state_setscalar, state_diagnose, memcheck
  use dshr_methods_mod , only : set_component_logging, log_clock_advance
  use dshr_nuopc_mod   , only : dshr_advertise, dshr_model_initphase, dshr_set_runclock
  use dshr_nuopc_mod   , only : dshr_sdat_init
  use dshr_nuopc_mod   , only : dshr_restart_read, dshr_restart_write
  use dice_comp_mod    , only : dice_comp_advertise, dice_comp_realize, dice_comp_run
  use dice_comp_mod    , only : water  ! for restart
  use perf_mod         , only : t_startf, t_stopf, t_adj_detailf, t_barrierf

  implicit none
  private ! except

  public  :: SetServices

  private :: InitializeAdvertise
  private :: InitializeRealize
  private :: ModelAdvance
  private :: ModelFinalize

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  type(shr_strdata_type)   :: sdat
  character(len=CS)        :: flds_scalar_name = ''
  integer                  :: flds_scalar_num = 0
  integer                  :: flds_scalar_index_nx = 0
  integer                  :: flds_scalar_index_ny = 0
  type(ESMF_Mesh)          :: mesh                      ! model mesh
  integer                  :: compid                    ! mct comp id
  integer                  :: mpicom                    ! mpi communicator
  integer                  :: my_task                   ! my task in mpi communicator mpicom
  character(len=16)        :: inst_suffix = ""          ! char string associated with instance (ie. "_0001" or "")
  integer                  :: logunit                   ! logging unit number
  logical                  :: read_restart              ! start from restart
  character(*) , parameter :: nullstr = 'undefined'
  integer      , parameter :: master_task=0             ! task number of master task
  character(*) , parameter :: rpfile = 'rpointer.ice'
  character(*) , parameter :: modName =  "(ice_comp_nuopc)"

  ! constants
  real(R8), parameter :: pi  = shr_const_pi ! pi
  real(R8)            :: dt           ! real model timestep

  ! nuopc attributes
  logical       :: flds_i2o_per_cat   ! .true. if select per ice thickness
  character(CS) :: calendar           ! calendar name

  ! dice_in namelist input
  real(R8)      :: flux_swpf          ! short-wave penatration factor
  real(R8)      :: flux_Qmin          ! bound on melt rate
  logical       :: flux_Qacc          ! activates water accumulation/melt wrt Q
  real(R8)      :: flux_Qacc0         ! initial water accumulation value
  character(CL) :: restfilm = nullstr ! model restart file namelist
  character(CL) :: restfils = nullstr ! stream restart file namelist
  character(CL) :: rest_file          ! restart filename
  character(CL) :: rest_file_strm     ! restart filename for streams

  character(*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Local varaibles
    character(len=*),parameter  :: subname=trim(modName)//':(SetServices) '
    !--------------------------------

    rc = ESMF_SUCCESS
    ! the NUOPC gcomp component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=dshr_model_initphase, phase=0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p3"/), userRoutine=InitializeRealize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! attach specializing method(s)
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, specRoutine=ModelAdvance, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MethodRemove(gcomp, label=model_label_SetRunClock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_SetRunClock, specRoutine=dshr_set_runclock, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, specRoutine=ModelFinalize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine SetServices

  !===============================================================================

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    integer           :: inst_index ! number of current instance (ie. 1)
    character(len=CL) :: cvalue     ! temporary
    integer           :: shrlogunit ! original log unit
    character(len=CL) :: fileName   ! generic file name
    integer           :: nu         ! unit number
    integer           :: ierr       ! error code
    character(len=*),parameter  :: subname=trim(modName)//':(InitializeAdvertise) '
    !-------------------------------------------------------------------------------

    namelist / dice_nml / flux_swpf, flux_Qmin, flux_Qacc, flux_Qacc0, restfilm, restfils

    rc = ESMF_SUCCESS

    ! obtain flds_scalar values, mpi values and multi-instance values
    call dshr_advertise(gcomp, mpicom, my_task,  inst_index, inst_suffix, &
         flds_scalar_name, flds_scalar_num, flds_scalar_index_nx, flds_scalar_index_ny, rc)

    ! set logunit and set shr logging to my log file
    call set_component_logging(gcomp, my_task==master_task, logunit, shrlogunit, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! set input namelist filename
    filename = "dice_in"//trim(inst_suffix)

    ! Read dice_nml from filename
    if (my_task == master_task) then
       open (newunit=nu,file=trim(filename),status="old",action="read")
       read (nu,nml=dice_nml,iostat=ierr)
       close(nu)
       if (ierr > 0) then
          write(logunit,*) 'ERROR: reading input namelist, '//trim(filename)//' iostat=',ierr
          call shr_sys_abort(subName//': namelist read error '//trim(filename))
       end if
       write(logunit,*)' flux_swpf  = ',flux_swpf
       write(logunit,*)' flux_Qmin  = ',flux_Qmin
       write(logunit,*)' flux_Qacc  = ',flux_Qacc
       write(logunit,*)' flux_Qacc0 = ',flux_Qacc0
       write(logunit,*)' restfilm   = ',trim(restfilm)
       write(logunit,*)' restfils   = ',trim(restfils)
    endif
    call shr_mpi_bcast(flux_swpf ,mpicom,'flux_swpf')
    call shr_mpi_bcast(flux_Qmin ,mpicom,'flux_Qmin')
    call shr_mpi_bcast(flux_Qacc ,mpicom,'flux_Qacc')
    call shr_mpi_bcast(flux_Qacc0,mpicom,'flux_Qacc0')
    call shr_mpi_bcast(restfilm,mpicom,'restfilm')
    call shr_mpi_bcast(restfils,mpicom,'restfils')
    rest_file      = trim(restfilm)
    rest_file_strm = trim(restfils)

    ! Read shr_strdata_nml from filename
    ! Read sdat namelist (need to do this here in order to get the datamode value - which
    ! is needed or order to do the advertise phase
    call shr_strdata_readnml(sdat, trim(filename), mpicom=mpicom)

    ! Validate sdat%datamode
    if (trim(sdat%datamode) == 'NULL' .or. trim(sdat%datamode) == 'SSTDATA' .or. trim(sdat%datamode) == 'COPYALL') then
       if (my_task == master_task) write(logunit,*) ' dice datamode = ',trim(sdat%datamode)
    else
       call shr_sys_abort(' ERROR illegal dice datamode = '//trim(sdat%datamode))
    endif

    ! Advertise import and export fields
    if (sdat%datamode /= 'NULL') then
       call NUOPC_CompAttributeGet(gcomp, name='flds_i2o_per_cat', value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) flds_i2o_per_cat  ! module variable

       call dice_comp_advertise(importState, exportState, flds_scalar_name, flds_i2o_per_cat, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Reset shr logging to original values
    call shr_file_setLogUnit (shrlogunit)

  end subroutine InitializeAdvertise

  !===============================================================================

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_TimeInterval) :: TimeStep
    type(ESMF_Time)         :: currTime
    type(ESMF_Calendar)     :: esmf_calendar ! esmf calendar
    type(ESMF_CalKind_Flag) :: esmf_caltype  ! esmf calendar type
    integer                 :: current_ymd   ! model date
    integer                 :: current_year  ! model year
    integer                 :: current_mon   ! model month
    integer                 :: current_day   ! model day
    integer                 :: current_tod   ! model sec into model date
    real(R8)                :: cosarg        ! for setting ice temp pattern
    real(R8)                :: jday, jday0   ! elapsed day counters
    character(CL)           :: cvalue        ! temporary
    integer                 :: shrlogunit    ! original log unit
    integer                 :: n,k           ! generic counters
    logical                 :: scmMode       ! single column mode
    real(R8)                :: scmLat        ! single column lat
    real(R8)                :: scmLon        ! single column lon
    integer                 :: model_dt      ! integer model timestep
    integer, allocatable, target :: gindex(:)
    character(len=*), parameter :: F00   = "('ice_comp_nuopc: ')',8a)"
    character(len=*), parameter :: subname=trim(modName)//':(InitializeRealize) '
    !-------------------------------------------------------------------------------

    if (sdat%datamode == 'NULL') RETURN

    rc = ESMF_SUCCESS

    ! Reset shr logging to my log file
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_setLogUnit (logUnit)

    ! Create the data model mesh
    call NUOPC_CompAttributeGet(gcomp, name='mesh_ice', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (my_task == master_task) then
       write(logunit,*) ' obtaining mesh_ice from '//trim(cvalue)
    end if
    mesh = ESMF_MeshCreate(filename=trim(cvalue), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get compid (for mct)
    call NUOPC_CompAttributeGet(gcomp, name='MCTID', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) compid

    ! Set single column values
    call NUOPC_CompAttributeGet(gcomp, name='scmlon', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scmlon
    call NUOPC_CompAttributeGet(gcomp, name='scmlat', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scmlat
    call NUOPC_CompAttributeGet(gcomp, name='single_column', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scmMode

    ! Initialize sdat
    call t_startf('dice_strdata_init')
    call dshr_sdat_init(mpicom, compid, my_task, master_task, logunit, &
         scmmode, scmlon, scmlat, clock, mesh, 'dice', sdat, use_new=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (my_task == master_task) write(logunit,*) ' initialized SDAT'
    call t_stopf('dice_strdata_init')

    ! Realize the actively coupled fields, now that a mesh is established and
    ! initialize dfields data type (to map streams to export state fields)
    call dice_comp_realize(sdat, importState, exportState, flds_scalar_name, flds_scalar_num, mesh, &
         logunit, my_task==master_task, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Read restart if necessary
    call NUOPC_CompAttributeGet(gcomp, name='read_restart', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) read_restart
    if (read_restart) then
       call dshr_restart_read(rest_file, rest_file_strm, rpfile, inst_suffix, nullstr, &
            logunit, my_task, master_task, mpicom, sdat, fld=water, fldname='water')
    end if

    ! get the time to interpolate the stream data to
    call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(currTime, yy=current_year, mm=current_mon, dd=current_day, s=current_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(current_year, current_mon, current_day, current_ymd)

    ! get model timestep
    call ESMF_TimeIntervalGet( timeStep, s=model_dt, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    dt = model_dt * 1.0_r8

    ! get cosarg
    call ESMF_ClockGet( clock, calkindflag=esmf_caltype, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (esmf_caltype == ESMF_CALKIND_NOLEAP) then
       calendar = shr_cal_noleap
    else if (esmf_caltype == ESMF_CALKIND_GREGORIAN) then
       calendar = shr_cal_gregorian
    else
       call shr_sys_abort(" ERROR bad ESMF calendar name "//trim(calendar))
    end if
    call shr_cal_ymd2julian(0, current_mon, current_day, current_tod, jDay , calendar) ! julian day for model
    call shr_cal_ymd2julian(0, 9,           1,           0,           jDay0, calendar) ! julian day for Sept 1
    cosArg = 2.0_R8*pi*(jday - jday0)/365.0_R8

    ! run dice
    call dice_comp_run(mpicom, my_task, master_task, logunit, current_ymd, current_tod, sdat, &
         mesh, cosarg, flux_swpf, flux_Qmin, flux_Qacc, flux_Qacc0, dt, read_restart, rc)

    ! add scalars to export state
    call State_SetScalar(dble(sdat%nxg),flds_scalar_index_nx, exportState, flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_SetScalar(dble(sdat%nyg),flds_scalar_index_ny, exportState, flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! diagnostics
    call State_diagnose(exportState,subname//':ES',rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_file_setLogUnit (shrlogunit)

   end subroutine InitializeRealize

  !===============================================================================

  subroutine ModelAdvance(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_State)        :: importState, exportState
    type(ESMF_Clock)        :: clock
    type(ESMF_Alarm)        :: alarm
    type(ESMF_TimeInterval) :: timeStep
    type(ESMF_Time)         :: currTime, nextTime
    integer                 :: current_mon   ! model month
    integer                 :: current_day   ! model day
    integer                 :: current_tod   ! model sec into model date
    real(R8)                :: cosarg        ! for setting ice temp pattern
    real(R8)                :: jday, jday0   ! elapsed day counters
    integer                 :: shrlogunit    ! original log unit
    integer                 :: next_ymd      ! model date
    integer                 :: next_tod      ! model sec into model date
    integer                 :: yr            ! year
    integer                 :: mon           ! month
    integer                 :: day           ! day in month
    character(CL)           :: case_name     ! case name
    character(len=*),parameter :: subname=trim(modName)//':(ModelAdvance) '
    !-------------------------------------------------------------------------------

    if (sdat%datamode == 'NULL') RETURN

    rc = ESMF_SUCCESS

    call t_startf(subname)
    call memcheck(subname, 5, my_task == master_task)

    ! Reset shr logging to my log file
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_setLogUnit (logunit)

    ! Query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! For nuopc - the component clock is advanced at the end of the time interval
    ! For these to match for now - need to advance nuopc one timestep ahead for
    ! shr_strdata time interpolation
    call ESMF_ClockGet( clock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    nextTime = currTime + timeStep
    call ESMF_TimeGet( nextTime, yy=yr, mm=mon, dd=day, s=next_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yr, mon, day, next_ymd)

    ! Get cosarg
    call shr_cal_ymd2julian(0, mon, day, next_tod, jDay , calendar)    ! julian day for model
    call shr_cal_ymd2julian(0, 9,   1,   0,        jDay0, calendar)    ! julian day for Sept 1
    cosArg = 2.0_R8*pi*(jday - jday0)/365.0_R8

    ! Run dice
    call dice_comp_run(mpicom, my_task, master_task, logunit, next_ymd, next_tod, sdat, &
         mesh, cosarg, flux_swpf, flux_Qmin, flux_Qacc, flux_Qacc0, dt, read_restart, rc)

    ! Write_restart if alarm is ringing
    call ESMF_ClockGetAlarm(clock, alarmname='alarm_restart', alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call t_startf('dice_restart')
       call NUOPC_CompAttributeGet(gcomp, name='case_name', value=case_name, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call dshr_restart_write(rpfile, case_name, 'dice', inst_suffix, next_ymd, next_tod, &
            logunit, mpicom, my_task, master_task, sdat, fld=water, fldname='water')
       call t_stopf('dice_restart')
    endif

    ! Write diagnostics
    call State_diagnose(exportState,subname//':ES',rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (my_task == master_task) then
       call log_clock_advance(clock, 'dice', logunit, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    ! Reset shr logging to original values
    call shr_file_setLogUnit (shrlogunit)
    call t_stopf(subname)

  end subroutine ModelAdvance

  !===============================================================================

  subroutine ModelFinalize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (my_task == master_task) then
       write(logunit,*)
       write(logunit,*) 'dice : end of main integration loop'
       write(logunit,*)
    end if

  end subroutine ModelFinalize

end module ice_comp_nuopc
