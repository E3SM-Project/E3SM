module ice_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for DICE
  !----------------------------------------------------------------------------

  use ESMF
  use NUOPC                 , only : NUOPC_CompDerive, NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize
  use NUOPC                 , only : NUOPC_CompAttributeGet, NUOPC_Advertise
  use NUOPC_Model           , only : model_routine_SS        => SetServices
  use NUOPC_Model           , only : model_label_Advance     => label_Advance
  use NUOPC_Model           , only : model_label_SetRunClock => label_SetRunClock
  use NUOPC_Model           , only : model_label_Finalize    => label_Finalize
  use NUOPC_Model           , only : NUOPC_ModelGet
  use med_constants_mod     , only : R8, CXX, CL, CS
  use med_constants_mod     , only : shr_log_Unit
  use med_constants_mod     , only : shr_file_getlogunit, shr_file_setlogunit
  use med_constants_mod     , only : shr_file_getloglevel, shr_file_setloglevel
  use med_constants_mod     , only : shr_file_setIO, shr_file_getUnit
  use med_constants_mod     , only : shr_cal_ymd2date, shr_cal_noleap, shr_cal_gregorian
  use shr_nuopc_scalars_mod , only : flds_scalar_name
  use shr_nuopc_scalars_mod , only : flds_scalar_num
  use shr_nuopc_scalars_mod , only : flds_scalar_index_nx
  use shr_nuopc_scalars_mod , only : flds_scalar_index_ny
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_Clock_TimePrint
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_SetScalar
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_Diagnose
  use shr_nuopc_grid_mod    , only : shr_nuopc_grid_ArrayToState
  use shr_nuopc_grid_mod    , only : shr_nuopc_grid_StateToArray
  use shr_const_mod         , only : SHR_CONST_SPVAL
  use shr_strdata_mod       , only : shr_strdata_type
  use shr_cal_mod           , only : shr_cal_ymd2julian
  use shr_const_mod         , only : shr_const_pi
  use dshr_nuopc_mod        , only : fld_list_type, fldsMax, fld_list_realize
  use dshr_nuopc_mod        , only : ModelInitPhase, ModelSetRunClock, ModelSetMetaData
  use dice_shr_mod          , only : dice_shr_read_namelists
  use dice_comp_mod         , only : dice_comp_init, dice_comp_run, dice_comp_advertise
  use mct_mod               , only : mct_Avect, mct_Avect_info


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

  integer                    :: fldsToIce_num = 0
  integer                    :: fldsFrIce_num = 0
  type (fld_list_type)       :: fldsToIce(fldsMax)
  type (fld_list_type)       :: fldsFrIce(fldsMax)
  type(shr_strdata_type)     :: SDICE
  type(mct_aVect)            :: x2i
  type(mct_aVect)            :: i2x
  integer                    :: compid                    ! mct comp id
  integer                    :: mpicom                    ! mpi communicator
  integer                    :: my_task                   ! my task in mpi communicator mpicom
  integer                    :: inst_index                ! number of current instance (ie. 1)
  character(len=16)          :: inst_name                 ! fullname of current instance (ie. "lnd_0001")
  character(len=16)          :: inst_suffix = ""          ! char string associated with instance (ie. "_0001" or "")
  integer                    :: logunit                   ! logging unit number
  integer, parameter         :: master_task=0             ! task number of master task
  logical                    :: read_restart              ! start from restart
  character(len=256)         :: case_name                 ! case name
  integer                    :: dbrc
  logical                    :: flds_i2o_per_cat          ! .true. if select per ice thickness
                                                          ! category fields are passed from ice to ocean
  character(len=80)          :: calendar                  ! calendar name
  integer                    :: modeldt                   ! integer timestep
  character(len=CXX)         :: flds_i2x = ''
  character(len=CXX)         :: flds_x2i = ''
  logical                    :: use_esmf_metadata = .false.
  real(R8)    ,parameter     :: pi  = shr_const_pi      ! pi
  character(*),parameter     :: modName =  "(ice_comp_nuopc)"
  integer, parameter         :: debug_import = 0          ! if > 0 will diagnose import fields
  integer, parameter         :: debug_export = 0          ! if > 0 will diagnose export fields
  character(*),parameter     :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter  :: subname=trim(modName)//':(SetServices) '

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    ! the NUOPC gcomp component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=ModelInitPhase, phase=0, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p3"/), userRoutine=InitializeRealize, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! attach specializing method(s)

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
         specRoutine=ModelAdvance, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MethodRemove(gcomp, label=model_label_SetRunClock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_SetRunClock, &
         specRoutine=ModelSetRunClock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
         specRoutine=ModelFinalize, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine SetServices

  !===============================================================================

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)
    use shr_nuopc_utils_mod, only : shr_nuopc_set_component_logging
    use shr_nuopc_utils_mod, only : shr_nuopc_get_component_instance

    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    logical            :: ice_present    ! flag
    logical            :: ice_prognostic ! flag
    type(ESMF_VM)      :: vm
    integer            :: lmpicom
    character(len=CL)  :: cvalue
    integer            :: n
    integer            :: ierr           ! error code
    integer            :: shrlogunit     ! original log unit
    integer            :: shrloglev      ! original log level
    logical            :: isPresent
    character(len=CL)  :: diro
    character(len=CL)  :: logfile
    integer            :: localPet
    character(len=CL)  :: fileName    ! generic file name
    character(len=*),parameter :: subname=trim(modName)//':(InitializeAdvertise) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    !----------------------------------------------------------------------------
    ! generate local mpi comm
    !----------------------------------------------------------------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, mpiCommunicator=lmpicom, localPet=localPet, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call mpi_comm_dup(lmpicom, mpicom, ierr)
    call mpi_comm_rank(mpicom, my_task, ierr)

    !----------------------------------------------------------------------------
    ! determine instance information
    !----------------------------------------------------------------------------

    call shr_nuopc_get_component_instance(gcomp, inst_suffix, inst_index)
    inst_name = "ICE"//trim(inst_suffix)

    !----------------------------------------------------------------------------
    ! set logunit and set shr logging to my log file
    !----------------------------------------------------------------------------

    call shr_nuopc_set_component_logging(gcomp, my_task==master_task, logunit, shrlogunit, shrloglev)

    !----------------------------------------------------------------------------
    ! Read input namelists and set present and prognostic flags
    !----------------------------------------------------------------------------

    filename = "dice_in"//trim(inst_suffix)
    call dice_shr_read_namelists(filename, mpicom, my_task, master_task, &
         logunit, SDICE, ice_present, ice_prognostic)

    !--------------------------------
    ! Advertise import and export fields
    !--------------------------------

    call dice_comp_advertise(importstate, exportState, &
       ice_present, ice_prognostic, &
       fldsFrIce_num, fldsFrIce, fldsToIce_num, fldsToIce, &
       flds_i2x, flds_x2i, rc)

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

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
    type(ESMF_Calendar)     :: esmf_calendar             ! esmf calendar
    type(ESMF_CalKind_Flag) :: esmf_caltype              ! esmf calendar type
    integer                 :: nx_global, ny_global      ! global sizes
    integer                 :: n                         ! index
    character(len=256)      :: cvalue                    ! tempoaray character string
    integer                 :: shrlogunit                ! original log unit
    integer                 :: shrloglev                 ! original log level
    logical                 :: scmMode = .false.         ! single column mode
    real(R8)                :: scmLat  = shr_const_SPVAL ! single column lat
    real(R8)                :: scmLon  = shr_const_SPVAL ! single column lon
    integer                 :: current_ymd               ! model date
    integer                 :: current_year              ! model year
    integer                 :: current_mon               ! model month
    integer                 :: current_day               ! model day
    integer                 :: current_tod               ! model sec into model date
    real(R8)                :: cosarg                    ! for setting ice temp pattern
    real(R8)                :: jday, jday0               ! elapsed day counters
    logical                 :: write_restart
    character(len=*),parameter :: subname=trim(modName)//':(InitializeRealize) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogLevel(max(shrloglev,1))
    call shr_file_setLogUnit (logUnit)

    !--------------------------------
    ! get config variables
    !--------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=case_name, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name='scmlon', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scmlon

    call NUOPC_CompAttributeGet(gcomp, name='scmlat', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scmlat

    call NUOPC_CompAttributeGet(gcomp, name='single_column', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scmMode

    call NUOPC_CompAttributeGet(gcomp, name='read_restart', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) read_restart

    call NUOPC_CompAttributeGet(gcomp, name='flds_i2o_per_cat', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_i2o_per_cat  ! module variable

    call NUOPC_CompAttributeGet(gcomp, name='MCTID', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) compid

    !----------------------------------------------------------------------------
    ! Determine calendar info
    !----------------------------------------------------------------------------

    call ESMF_ClockGet( clock, calkindflag=esmf_caltype, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (esmf_caltype == ESMF_CALKIND_NOLEAP) then
       calendar = shr_cal_noleap
    else if (esmf_caltype == ESMF_CALKIND_GREGORIAN) then
       calendar = shr_cal_gregorian
    else
       call ESMF_LogWrite(subname//" ERROR bad ESMF calendar name "//trim(calendar), &
            ESMF_LOGMSG_ERROR, rc=dbrc)
       rc = ESMF_Failure
       return
    end if

    !--------------------------------
    ! Generate the mesh
    !--------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='mesh_ice', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    Emesh = ESMF_MeshCreate(filename=trim(cvalue), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Initialize model
    !--------------------------------

    call dice_comp_init(x2i, i2x, &
         flds_x2i, flds_i2x, flds_i2o_per_cat, &
         SDICE, mpicom, compid, my_task, master_task, &
         inst_suffix, inst_name, logunit, read_restart, &
         scmMode, scmlat, scmlon, calendar, Emesh)

    !--------------------------------
    ! realize the actively coupled fields, now that a mesh is established
    ! NUOPC_Realize "realizes" a previously advertised field in the importState and exportState
    ! by replacing the advertised fields with the newly created fields of the same name.
    !--------------------------------

    call fld_list_realize( &
         state=ExportState, &
         fldList=fldsFrIce, &
         numflds=fldsFrIce_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':diceExport',&
         mesh=Emesh, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call fld_list_realize( &
         state=importState, &
         fldList=fldsToIce, &
         numflds=fldsToIce_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':diceImport',&
         mesh=Emesh, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------------------------
    ! Set initial ice state and pack export state
    !----------------------------------------------------------------------------

    call ESMF_ClockGet( clock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet( currTime, yy=current_year, mm=current_mon, dd=current_day, s=current_tod, &
         calkindflag=esmf_caltype, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(current_year, current_mon, current_day, current_ymd)

    call ESMF_TimeIntervalGet( timeStep, s=modeldt, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_cal_ymd2julian(0, current_mon, current_day, current_tod, jDay , calendar)    ! julian day for model
    call shr_cal_ymd2julian(0, 9,           1,           0,           jDay0, calendar)    ! julian day for Sept 1
    cosArg = 2.0_R8*pi*(jday - jday0)/365.0_R8

    write_restart = .false.
    call dice_comp_run(x2i, i2x, &
         flds_i2o_per_cat, SDICE, mpicom, my_task, master_task, &
         inst_suffix, logunit, read_restart, write_restart, &
         calendar, modeldt, current_ymd, current_tod, cosArg)

    ! Pack export state
    call shr_nuopc_grid_ArrayToState(i2x%rattr, flds_i2x, exportState, grid_option='mesh', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    nx_global = SDICE%nxg
    ny_global = SDICE%nyg
    call shr_nuopc_methods_State_SetScalar(dble(nx_global),flds_scalar_index_nx, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_State_SetScalar(dble(ny_global),flds_scalar_index_ny, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (debug_export > 0) then
       call shr_nuopc_methods_State_diagnose(exportState,subname//':ES',rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    call shr_file_setLogLevel(shrloglev)
    call shr_file_setLogUnit (shrlogunit)

    if (use_esmf_metadata) then
       call ModelSetMetaData(gcomp, name='DICE', rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine InitializeRealize

  !===============================================================================

  subroutine ModelAdvance(gcomp, rc)
    use shr_nuopc_utils_mod, only : shr_nuopc_memcheck, shr_nuopc_log_clock_advance
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)        :: clock
    type(ESMF_Alarm)        :: alarm
    type(ESMF_Time)         :: currTime, nextTime
    type(ESMF_TimeInterval) :: timeStep
    type(ESMF_State)        :: importState, exportState
    integer                 :: shrlogunit    ! original log unit
    integer                 :: shrloglev     ! original log level
    logical                 :: write_restart ! restart alarm is ringing
    logical                 :: read_restart  ! read restart flag
    integer                 :: next_ymd      ! model date
    integer                 :: next_tod      ! model sec into model date
    integer                 :: yr            ! year
    integer                 :: mon           ! month
    integer                 :: day           ! day in month
    real(R8)                :: cosarg        ! for setting ice temp pattern
    real(R8)                :: jday, jday0   ! elapsed day counters
    character(len=CL)       :: cvalue
    character(len=*),parameter  :: subname=trim(modName)//':(ModelAdvance) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)
    call shr_nuopc_memcheck(subname, 5, my_task==master_task)
    !--------------------------------
    ! Reset shr logging to my log file
    !--------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogLevel(max(shrloglev,1))
    call shr_file_setLogUnit (logunit)

    !--------------------------------
    ! query the Component for its clock, importState and exportState
    !--------------------------------

    call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (debug_import > 0 .and. my_task == master_task) then
       call shr_nuopc_methods_Clock_TimePrint(clock,subname//'clock',rc=rc)
    end if

    !--------------------------------
    ! Unpack import state
    !--------------------------------

    call shr_nuopc_grid_StateToArray(importState, x2i%rattr, flds_x2i, grid_option='mesh', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Run model
    !--------------------------------

    ! Determine if will write restart

    call ESMF_ClockGetAlarm(clock, alarmname='alarm_restart', alarm=alarm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       write_restart = .true.
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       write_restart = .false.
    endif

    ! For nuopc - the component clock is advanced at the end of the time interval
    ! For these to match for now - need to advance nuopc one timestep ahead for
    ! shr_strdata time interpolation

    call ESMF_ClockGet( clock, currTime=currTime, timeStep=timeStep, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    nextTime = currTime + timeStep
    call ESMF_TimeGet( nextTime, yy=yr, mm=mon, dd=day, s=next_tod, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yr, mon, day, next_ymd)

    call shr_cal_ymd2julian(0, mon, day, next_tod, jDay , calendar)    ! julian day for model
    call shr_cal_ymd2julian(0, 9,   1,   0,        jDay0, calendar)    ! julian day for Sept 1
    cosArg = 2.0_R8*pi*(jday - jday0)/365.0_R8

    ! Run dice

    read_restart = .false.
    call dice_comp_run(x2i, i2x, &
         flds_i2o_per_cat, SDICE, mpicom, my_task, master_task, &
         inst_suffix, logunit, read_restart, write_restart, &
         calendar, modeldt, next_ymd, next_tod, cosArg, case_name)

    !--------------------------------
    ! Pack export state
    !--------------------------------

    call shr_nuopc_grid_ArrayToState(i2x%rattr, flds_i2x, exportState, grid_option='mesh', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (debug_export > 0) then
       call shr_nuopc_methods_State_diagnose(exportState,subname//':ES',rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    if (my_task == master_task) then
       call shr_nuopc_log_clock_advance(clock, 'ICE', logunit)
    end if

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    call shr_file_setLogLevel(shrloglev)
    call shr_file_setLogUnit (shrlogunit)

  end subroutine ModelAdvance

  !===============================================================================

  subroutine ModelFinalize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(*), parameter :: F00   = "('(dice_comp_final) ',8a)"
    character(*), parameter :: F91   = "('(dice_comp_final) ',73('-'))"
    character(len=*),parameter  :: subname=trim(modName)//':(ModelFinalize) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (my_task == master_task) then
       write(logunit,F91)
       write(logunit,F00) ' dice: end of main integration loop'
       write(logunit,F91)
    end if
    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine ModelFinalize

end module ice_comp_nuopc
