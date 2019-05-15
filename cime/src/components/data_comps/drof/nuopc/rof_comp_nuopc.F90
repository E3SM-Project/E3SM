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
  use shr_cal_mod      , only : shr_cal_noleap, shr_cal_gregorian, shr_cal_ymd2date
  use dshr_nuopc_mod   , only : fld_list_type, fldsMax, dshr_realize
  use dshr_nuopc_mod   , only : ModelInitPhase, ModelSetRunClock, ModelSetMetaData
  use dshr_methods_mod , only : chkerr, state_setscalar,  state_diagnose, alarmInit, memcheck
  use dshr_methods_mod , only : set_component_logging, get_component_instance, log_clock_advance
  use drof_shr_mod     , only : drof_shr_read_namelists
  use drof_comp_mod    , only : drof_comp_init, drof_comp_run, drof_comp_advertise
  use drof_comp_mod    , only : drof_comp_export

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

  character(len=CS)      :: flds_scalar_name = ''
  integer                :: flds_scalar_num = 0
  integer                :: flds_scalar_index_nx = 0
  integer                :: flds_scalar_index_ny = 0

  integer                :: fldsToRof_num = 0
  integer                :: fldsFrRof_num = 0
  type (fld_list_type)   :: fldsToRof(fldsMax)
  type (fld_list_type)   :: fldsFrRof(fldsMax)

  integer                :: compid                    ! mct comp id
  integer                :: mpicom                    ! mpi communicator
  integer                :: my_task                   ! my task in mpi communicator mpicom
  integer                :: inst_index                ! number of current instance (ie. 1)
  character(len=16)      :: inst_suffix = ""          ! char string associated with instance (ie. "_0001" or "")
  integer                :: logunit                   ! logging unit number
  integer    ,parameter  :: master_task=0             ! task number of master task
  logical                :: rof_prognostic            ! flag
  character(CL)          :: case_name                 ! case name
  character(CL)          :: tmpstr                    ! tmp string
  character(len=80)      :: calendar                  ! calendar name
  logical                :: use_esmf_metadata = .false.
  character(*),parameter :: modName =  "(rof_comp_nuopc)"
  integer, parameter     :: debug_import = 0          ! if > 0 will diagnose import fields
  integer, parameter     :: debug_export = 0          ! if > 0 will diagnose export fields
  character(*),parameter :: u_FILE_u = &
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
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! the NUOPC gcomp component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=ModelInitPhase, phase=0, rc=rc)
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
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_SetRunClock, specRoutine=ModelSetRunClock, rc=rc)
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
    logical            :: rof_present    ! flag
    type(ESMF_VM)      :: vm
    integer            :: lmpicom
    character(len=CL)  :: cvalue
    integer            :: n
    integer            :: ierr           ! error code
    integer            :: shrlogunit     ! original log unit
    character(len=CL)  :: diro
    character(len=CL)  :: logfile
    integer            :: localPet
    character(len=CL)  :: fileName    ! generic file name
    character(len=CL)  :: logmsg
    logical            :: isPresent, isSet
    character(len=*),parameter  :: subname=trim(modName)//':(InitializeAdvertise) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    !----------------------------------------------------------------------------
    ! generate local mpi comm
    !----------------------------------------------------------------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, mpiCommunicator=lmpicom, localPet=localPet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call mpi_comm_dup(lmpicom, mpicom, ierr)
    call mpi_comm_rank(mpicom, my_task, ierr)

    !----------------------------------------------------------------------------
    ! determine instance information
    !----------------------------------------------------------------------------

    call get_component_instance(gcomp, inst_suffix, inst_index, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------------------------
    ! set logunit and set shr logging to my log file
    !----------------------------------------------------------------------------

    call set_component_logging(gcomp, my_task==master_task, logunit, shrlogunit, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------------------------
    ! Read input namelists and set present and prognostic flags
    !----------------------------------------------------------------------------

    filename = "drof_in"//trim(inst_suffix)
    call drof_shr_read_namelists(filename, mpicom, my_task, master_task, &
         logunit, rof_present, rof_prognostic)

    !--------------------------------
    ! advertise export fields
    !--------------------------------

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       flds_scalar_name = trim(cvalue)
       call ESMF_LogWrite(trim(subname)//' flds_scalar_name = '//trim(flds_scalar_name), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldCount", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue, *) flds_scalar_num
       write(logmsg,*) flds_scalar_num
       call ESMF_LogWrite(trim(subname)//' flds_scalar_num = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNX", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_nx
       write(logmsg,*) flds_scalar_index_nx
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_nx = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxGridNY", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue,*) flds_scalar_index_ny
       write(logmsg,*) flds_scalar_index_ny
       call ESMF_LogWrite(trim(subname)//' : flds_scalar_index_ny = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    call drof_comp_advertise(importState, exportState, flds_scalar_name, &
         rof_present, rof_prognostic, &
         fldsFrRof_num, fldsFrRof, fldsToRof_num, fldsToRof, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    call shr_file_setLogUnit (shrlogunit)

  end subroutine InitializeAdvertise

  !===============================================================================

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Mesh)         :: Emesh
    type(ESMF_TIME)         :: currTime
    type(ESMF_TimeInterval) :: timeStep
    type(ESMF_Calendar)     :: esmf_calendar ! esmf calendar
    type(ESMF_CalKind_Flag) :: esmf_caltype  ! esmf calendar type
    character(CL)           :: cvalue
    integer                 :: shrlogunit    ! original log unit
    integer                 :: ierr          ! error code
    integer                 :: current_ymd   ! model date
    integer                 :: current_year  ! model year
    integer                 :: current_mon   ! model month
    integer                 :: current_day   ! model day
    integer                 :: current_tod   ! model sec into model date
    logical                 :: read_restart  ! start from restart
    integer                 :: nxg, nyg
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

    !--------------------------------
    ! Determine necessary config variables
    !--------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=case_name, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name='read_restart', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) read_restart

    call NUOPC_CompAttributeGet(gcomp, name='MCTID', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) compid

    !----------------------------------------------------------------------------
    ! Determine calendar info
    !----------------------------------------------------------------------------

    call ESMF_ClockGet( clock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet( currTime, yy=current_year, mm=current_mon, dd=current_day, s=current_tod, &
         calkindflag=esmf_caltype, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(current_year, current_mon, current_day, current_ymd)

    if (esmf_caltype == ESMF_CALKIND_NOLEAP) then
       calendar = shr_cal_noleap
    else if (esmf_caltype == ESMF_CALKIND_GREGORIAN) then
       calendar = shr_cal_gregorian
    else
       call ESMF_LogWrite(subname//" ERROR bad ESMF calendar name "//trim(calendar), ESMF_LOGMSG_ERROR)
       rc = ESMF_Failure
       return
    end if

    !--------------------------------
    ! Generate the mesh
    !--------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='mesh_rof', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (my_task == master_task) then
       write(logunit,*) " obtaining drof mesh from " // trim(cvalue)
    end if

    Emesh = ESMF_MeshCreate(filename=trim(cvalue), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------------------------
    ! Initialize model
    !----------------------------------------------------------------------------

    call drof_comp_init(mpicom, compid, my_task, master_task, &
         inst_suffix, logunit, read_restart, &
         current_ymd, current_tod, calendar, Emesh, nxg, nyg)

    !--------------------------------
    ! realize the actively coupled fields, now that a mesh is established
    ! NUOPC_Realize "realizes" a previously advertised field in the importState and exportState
    ! by replacing the advertised fields with the newly created fields of the same name.
    !--------------------------------

    call dshr_realize( &
         state=ExportState, &
         fldList=fldsFrRof, &
         numflds=fldsFrRof_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':drofExport',&
         mesh=Emesh, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! No import state for now

    !--------------------------------
    ! Pack export state
    ! Copy from r2x to exportState
    ! Set the coupling scalars
    !--------------------------------

    call drof_comp_export(exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call State_SetScalar(dble(nxg),flds_scalar_index_nx, exportState,  &
         flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call State_SetScalar(dble(nyg),flds_scalar_index_ny, exportState,  &
         flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (debug_export > 0) then
       call State_diagnose(exportState,subname//':ES',rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    call shr_file_setLogUnit (shrlogunit)

    if (use_esmf_metadata) then
       call ModelSetMetaData(gcomp, name='DROF', rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine InitializeRealize

  !===============================================================================

  subroutine ModelAdvance(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)        :: clock
    type(ESMF_Alarm)        :: alarm
    type(ESMF_Time)         :: currTime, nextTime
    type(ESMF_TimeInterval) :: timeStep
    type(ESMF_State)        :: importState, exportState
    integer                 :: shrlogunit    ! original log unit
    logical                 :: write_restart ! write restart
    logical                 :: read_restart  ! read restart
    integer                 :: yr            ! year
    integer                 :: mon           ! month
    integer                 :: day           ! day in month
    integer                 :: next_ymd      ! model date
    integer                 :: next_tod      ! model sec into model date
    character(len=*),parameter :: subname=trim(modName)//':(ModelAdvance) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
    call memcheck(subname, 5, my_task == master_task)

    !--------------------------------
    ! Reset shr logging to my log file
    !--------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_setLogUnit (logunit)

    !--------------------------------
    ! query the Component for its clock, importState and exportState
    !--------------------------------

    call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Unpack import state
    !--------------------------------

    if (rof_prognostic) then
       ! Do nothing for now
    end if

    !--------------------------------
    ! Run model
    !--------------------------------

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

    ! Run the model

    read_restart = .false.
    call drof_comp_run(mpicom, compid, my_task, master_task, &
         inst_suffix, logunit, read_restart, write_restart, &
         next_ymd, next_tod, case_name)

    !--------------------------------
    ! Pack export state
    !--------------------------------

    call drof_comp_export(exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (debug_export > 0) then
       call State_diagnose(exportState,subname//':ES',rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       if (my_task == master_task) then
          call log_clock_advance(clock, 'DROF', logunit, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    endif
    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

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

end module rof_comp_nuopc
