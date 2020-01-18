module rof_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for DROF
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
  use shr_const_mod    , only : SHR_CONST_SPVAL
  use shr_sys_mod      , only : shr_sys_abort
  use shr_cal_mod      , only : shr_cal_ymd2date
  use shr_strdata_mod  , only : shr_strdata_type, shr_strdata_readnml
  use shr_mpi_mod      , only : shr_mpi_bcast
  use dshr_nuopc_mod   , only : fld_list_type, fldsMax, dshr_realize
  use dshr_nuopc_mod   , only : dshr_advertise, dshr_model_initphase, dshr_set_runclock
  use dshr_nuopc_mod   , only : dshr_sdat_init, dshr_check_mesh 
  use dshr_nuopc_mod   , only : dshr_restart_read, dshr_restart_write
  use dshr_methods_mod , only : chkerr, state_setscalar,  state_diagnose, alarmInit, memcheck
  use dshr_methods_mod , only : set_component_logging, log_clock_advance
  use drof_comp_mod    , only : drof_comp_advertise, drof_comp_dfields_init, drof_comp_run 
  use drof_comp_mod    , only : fldsfrRof, fldsFrRof_num
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

  type(shr_strdata_type) :: sdrof

  character(len=CS)      :: flds_scalar_name = ''
  integer                :: flds_scalar_num = 0
  integer                :: flds_scalar_index_nx = 0
  integer                :: flds_scalar_index_ny = 0

  integer                :: compid                    ! mct comp id
  integer                :: mpicom                    ! mpi communicator
  integer                :: my_task                   ! my task in mpi communicator mpicom
  integer                :: inst_index                ! number of current instance (ie. 1)
  character(len=16)      :: inst_suffix = ""          ! char string associated with instance (ie. "_0001" or "")
  integer                :: logunit                   ! logging unit number
  logical                :: rof_prognostic            ! flag
  character(CL)          :: case_name                 ! case name
  character(CL)          :: tmpstr                    ! tmp string
  character(len=80)      :: calendar                  ! calendar name

  character(*) , parameter :: nullstr = 'undefined'
  integer      , parameter :: master_task=0             ! task number of master task
  character(*) , parameter :: rpfile = 'rpointer.rof'
  character(*) , parameter :: modName =  "(rof_comp_nuopc)"
  character(*) , parameter :: u_FILE_u = &
       __FILE__

  logical       :: force_prognostic_true = .false. ! if true set prognostic true
  character(CL) :: restfilm = nullstr              ! model restart file namelist
  character(CL) :: restfils = nullstr              ! stream restart file namelist
  character(CL) :: rest_file                       ! restart filename
  character(CL) :: rest_file_strm                  ! restart filename for streams

  integer      , parameter :: debug_import = 0          ! if > 0 will diagnose import fields
  integer      , parameter :: debug_export = 0          ! if > 0 will diagnose export fields

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
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

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

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine SetServices

  !===============================================================================

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    character(len=CL)  :: cvalue
    integer            :: shrlogunit ! original log unit
    character(len=CL)  :: fileName   ! generic file name
    character(len=*),parameter  :: subname=trim(modName)//':(InitializeAdvertise) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! obtain flds_scalar values, mpi values and multi-instance values
    call dshr_advertise(gcomp, mpicom, my_task,  inst_index, inst_suffix, &
         flds_scalar_name, flds_scalar_num, flds_scalar_index_nx, flds_scalar_index_ny, rc)

    ! set logunit and set shr logging to my log file
    call set_component_logging(gcomp, my_task==master_task, logunit, shrlogunit, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Read input namelists and set present and prognostic flags
    filename = "drof_in"//trim(inst_suffix)
    call drof_read_namelists(filename)

    ! Read sdrof namelist (need to do this here in order to get the datamode value - which
    ! is needed or order to do the advertise phase
    call shr_strdata_readnml(sdrof, trim(filename), mpicom=mpicom)

    ! Validate sdrof datamode
    if (trim(sdrof%dataMode) == 'NULL' .or. trim(sdrof%dataMode) == 'COPYALL') then
       if (my_task == master_task) write(logunit,*) 'drof datamode = ',trim(sdrof%dataMode)
    else
       call shr_sys_abort(' ERROR illegal drof datamode = '//trim(sdrof%dataMode))
    end if

    if (trim(sdrof%datamode) /= 'NULL') then
       call drof_comp_advertise(importState, exportState, flds_scalar_name, force_prognostic_true, rc=rc)
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
     type(ESMF_Mesh) :: mesh
    type(ESMF_TIME) :: currTime
    integer         :: current_ymd  ! model date
    integer         :: current_year ! model year
    integer         :: current_mon  ! model month
    integer         :: current_day  ! model day
    integer         :: current_tod  ! model sec into model date
    character(CL)   :: cvalue       ! temporary
    integer         :: shrlogunit   ! original log unit
    integer         :: n,k          ! generic counters
    logical         :: scmMode      ! single column mode
    real(R8)        :: scmLat       ! single column lat
    real(R8)        :: scmLon       ! single column lon
    logical         :: read_restart
    character(CS)   :: model_name
    integer, allocatable, target :: gindex(:)
    character(len=*), parameter :: F00   = "('rof_comp_nuopc: ')',8a)"
    character(len=*), parameter :: subname=trim(modName)//':(InitializeRealize) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_setLogUnit (logUnit)

    !----------------------------------------------------------------------------
    ! Create the data model mesh
    !----------------------------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='mesh_rof', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (my_task == master_task) then
       write(logunit,*) ' obtaining mesh_rof from '//trim(cvalue)
    end if

    mesh = ESMF_MeshCreate(filename=trim(cvalue), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------------------------
    ! Initialize sdrof
    !----------------------------------------------------------------------------

    call t_startf('drof_strdata_init')
    call NUOPC_CompAttributeGet(gcomp, name='MCTID', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) compid

    scmmode = .false.
    scmlon = shr_const_spval 
    scmlat = shr_const_spval 
    model_name = 'drof'

    call dshr_sdat_init(mpicom, compid, my_task, master_task, logunit, &
       scmmode, scmlon, scmlat, clock, mesh, model_name, sdrof, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (my_task == master_task) write(logunit,*) ' initialized SDROF'
    call t_stopf('drof_strdata_init')

    !----------------------------------------------------------------------------
    ! Check that mesh lats and lons correspond to those on the input domain file
    !----------------------------------------------------------------------------

    call dshr_check_mesh(mesh, sdrof, 'drof', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! realize the actively coupled fields, now that a mesh is established
    !--------------------------------

    ! NUOPC_Realize "realizes" a previously advertised field in the importState and exportState
    ! by replacing the advertised fields with the newly created fields of the same name.
    call dshr_realize( state=ExportState, fldList=fldsFrRof, numflds=fldsFrRof_num, &
         flds_scalar_name=flds_scalar_name, flds_scalar_num=flds_scalar_num, &
         tag=subname//':drofExport', mesh=mesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Initialize dfields data type (to map streams to export state fields)
    !--------------------------------

    call drof_comp_dfields_init(sdrof, exportState, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Read restart if necessary
    !--------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='read_restart', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) read_restart
    if (read_restart) then
       call dshr_restart_read(rest_file, rest_file_strm, rpfile, inst_suffix, nullstr, &
            logunit, my_task, master_task, mpicom, sdrof)
    end if

    !--------------------------------
    ! Run drof to create export state
    !--------------------------------

    ! get the time to interpolate the stream data to
    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(currTime, yy=current_year, mm=current_mon, dd=current_day, s=current_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(current_year, current_mon, current_day, current_ymd)

    ! run drof
    call drof_comp_run(mpicom, my_task, master_task, logunit, current_ymd, current_tod, sdrof, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! add scalars to export state
    call State_SetScalar(dble(sdrof%nxg),flds_scalar_index_nx, exportState, flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_SetScalar(dble(sdrof%nyg),flds_scalar_index_ny, exportState, flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! diagnostics
    call State_diagnose(exportState,subname//':ES',rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_file_setLogUnit (shrlogunit)

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

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
    integer                 :: shrlogunit    ! original log unit
    logical                 :: write_restart ! write restart
    integer                 :: next_ymd      ! model date
    integer                 :: next_tod      ! model sec into model date
    integer                 :: yr            ! year
    integer                 :: mon           ! month
    integer                 :: day           ! day in month
    character(len=*),parameter :: subname=trim(modName)//':(ModelAdvance) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
    call memcheck(subname, 5, my_task == master_task)

    ! Reset shr logging to my log file
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_setLogUnit (logunit)

    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine if need to write restarts
    call ESMF_ClockGetAlarm(clock, alarmname='alarm_restart', alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       write_restart = .true.
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       write_restart = .false.
    endif

    ! For nuopc - the component clock is advanced at the end of the time interval
    ! For these to match for now - need to advance nuopc one timestep ahead for
    ! shr_strdata time interpolation
    call ESMF_ClockGet( clock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    nextTime = currTime + timeStep
    call ESMF_TimeGet( nextTime, yy=yr, mm=mon, dd=day, s=next_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yr, mon, day, next_ymd)

    ! run drof
    call drof_comp_run(mpicom, my_task, master_task, logunit, next_ymd, next_tod, sdrof, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! write_restart if alarm is ringing
    call ESMF_ClockGetAlarm(clock, alarmname='alarm_restart', alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call t_startf('drof_restart')
       call NUOPC_CompAttributeGet(gcomp, name='case_name', value=case_name, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call dshr_restart_write(rpfile, case_name, 'drof', inst_suffix, next_ymd, next_tod, &
            logunit, mpicom, my_task, master_task, sdrof)
       call t_stopf('drof_restart')
    endif

    ! write diagnostics
    call State_diagnose(exportState,subname//':ES',rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (my_task == master_task) then
       call log_clock_advance(clock, 'DROF', logunit, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif
    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

    call shr_file_setLogUnit (shrlogunit)

  end subroutine ModelAdvance

  !===============================================================================

  subroutine ModelFinalize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(*), parameter :: F00   = "('(drof_comp_final) ',8a)"
    character(*), parameter :: F91   = "('(drof_comp_final) ',73('-'))"
    character(len=*),parameter  :: subname=trim(modName)//':(ModelFinalize) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
    if (my_task == master_task) then
       write(logunit,F91)
       write(logunit,F00) 'drof : end of main integration loop'
       write(logunit,F91)
    end if
    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelFinalize

  !===============================================================================

  subroutine drof_read_namelists(filename)

    ! Read in non-sdat part of drof_in

    ! input/output parameters
    character(len=*), intent(in) :: filename ! input namelist filename

    ! local variables
    integer       :: nu     ! unit number
    integer       :: ierr   ! error code
    character(CL) :: decomp ! decomp strategy - not used for NUOPC - but still needed in namelist for now
    character(*), parameter :: subName = "(shr_drof_read_namelists) "
    !-------------------------------------------------------------------------------

    namelist / drof_nml / decomp, restfilm, restfils, force_prognostic_true

    restfilm   = trim(nullstr)
    restfils   = trim(nullstr)

    if (my_task == master_task) then
       open (newunit=nu,file=trim(filename),status="old",action="read")
       read (nu,nml=drof_nml,iostat=ierr)
       close(nu)
       if (ierr > 0) then
          write(logunit,*) 'ERROR: reading input namelist, '//trim(filename)//' iostat=',ierr
          call shr_sys_abort(subName//': namelist read error '//trim(filename))
       end if
       write(logunit,*)' restfilm   = ',trim(restfilm)
       write(logunit,*)' restfils   = ',trim(restfils)
       write(logunit,*)' force_prognostic_true = ',force_prognostic_true
    endif
    call shr_mpi_bcast(restfilm,mpicom,'restfilm')
    call shr_mpi_bcast(restfils,mpicom,'restfils')
    call shr_mpi_bcast(force_prognostic_true,mpicom,'force_prognostic_true')
    rest_file = trim(restfilm)
    rest_file_strm = trim(restfils)

  end subroutine drof_read_namelists

end module rof_comp_nuopc
