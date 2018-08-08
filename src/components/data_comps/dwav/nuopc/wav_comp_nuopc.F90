module wav_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for DWAV
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
  use shr_nuopc_grid_mod    , only : shr_nuopc_grid_Meshinit
  use shr_nuopc_grid_mod    , only : shr_nuopc_grid_ArrayToState
  use shr_nuopc_grid_mod    , only : shr_nuopc_grid_StateToArray
  use shr_nuopc_time_mod    , only : shr_nuopc_time_alarmInit
  use shr_strdata_mod       , only : shr_strdata_type
  use dshr_nuopc_mod        , only : fld_list_type, fldsMax
  use dshr_nuopc_mod        , only : fld_list_add
  use dshr_nuopc_mod        , only : fld_list_realize
  use dshr_nuopc_mod        , only : ModelInitPhase
  use dshr_nuopc_mod        , only : ModelSetRunClock
  use dshr_nuopc_mod        , only : ModelSetMetaData

  use dwav_shr_mod , only: dwav_shr_read_namelists
  use dwav_comp_mod, only: dwav_comp_init, dwav_comp_run, dwav_comp_final
  use mct_mod

  implicit none

  public  :: SetServices

  private :: InitializeAdvertise
  private :: InitializeRealize
  private :: ModelAdvance
  private :: ModelFinalize

  private ! except

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  integer                    :: fldsToWav_num = 0
  integer                    :: fldsFrWav_num = 0
  type (fld_list_type)       :: fldsToWav(fldsMax)
  type (fld_list_type)       :: fldsFrWav(fldsMax)

  character(len=3)           :: myModelName = 'wav'       ! user defined model name
  type(shr_strdata_type)     :: SDWAV
  type(mct_gsMap), target    :: gsMap_target
  type(mct_gGrid), target    :: ggrid_target
  type(mct_gsMap), pointer   :: gsMap
  type(mct_gGrid), pointer   :: ggrid
  type(mct_aVect)            :: x2d
  type(mct_aVect)            :: d2x
  integer                    :: compid                    ! mct comp id
  integer                    :: mpicom                    ! mpi communicator
  integer                    :: my_task                   ! my task in mpi communicator mpicom
  integer                    :: inst_index                ! number of current instance (ie. 1)
  character(len=16)          :: inst_name                 ! fullname of current instance (ie. "wav_0001")
  character(len=16)          :: inst_suffix = ""          ! char string associated with instance (ie. "_0001" or "")
  integer                    :: logunit                   ! logging unit number
  integer    ,parameter      :: master_task=0             ! task number of master task
  integer                    :: localPet
  logical                    :: read_restart              ! start from restart
  logical                    :: wav_prognostic            ! flag
  logical                    :: unpack_import
  character(CL)              :: tmpstr                    ! tmp string
  integer                    :: dbrc
  integer, parameter         :: dbug = 10
  character(len=*),parameter :: grid_option = "mesh"      ! grid_de, grid_arb, grid_reg, mesh
  character(CXX)             :: flds_w2x = ''
  character(CXX)             :: flds_x2w = ''

  !----- formats -----
  character(*),parameter :: modName =  "(wav_comp_nuopc)"
  character(*),parameter :: u_FILE_u = __FILE__

  !===============================================================================
  contains
  !===============================================================================

  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter  :: subname=trim(modName)//':(SetServices) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

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

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, specRoutine=ModelAdvance, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_MethodRemove(gcomp, label=model_label_SetRunClock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_SetRunClock, specRoutine=ModelSetRunClock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, specRoutine=ModelFinalize, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine SetServices

  !===============================================================================

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    logical            :: wav_present ! flag
    type(ESMF_VM)      :: vm
    integer            :: lmpicom
    character(CL)      :: cvalue
    character(CS)      :: stdname, shortname
    logical            :: activefld
    integer            :: n,nflds
    integer            :: ierr        ! error code
    integer            :: shrlogunit  ! original log unit
    integer            :: shrloglev   ! original log level
    logical            :: isPresent
    character(len=512) :: diro
    character(len=512) :: logfile
    character(len=*),parameter  :: subname=trim(modName)//':(InitializeAdvertise) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    !----------------------------------------------------------------------------
    ! generate local mpi comm
    !----------------------------------------------------------------------------

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, mpiCommunicator=lmpicom, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call mpi_comm_dup(lmpicom, mpicom, ierr)
    call mpi_comm_rank(mpicom, my_task, ierr)

    !----------------------------------------------------------------------------
    ! determine instance information
    !----------------------------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name="inst_name", value=inst_name, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name="inst_index", value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) inst_index

    call ESMF_AttributeGet(gcomp, name="inst_suffix", isPresent=isPresent, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent) then
       call NUOPC_CompAttributeGet(gcomp, name="inst_suffix", value=inst_suffix, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       inst_suffix = ''
    end if

    !----------------------------------------------------------------------------
    ! set logunit and set shr logging to my log file
    !----------------------------------------------------------------------------

    if (my_task == master_task) then
       call NUOPC_CompAttributeGet(gcomp, name="diro", value=diro, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call NUOPC_CompAttributeGet(gcomp, name="logfile", value=logfile, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       logunit = shr_file_getUnit()
       open(logunit,file=trim(diro)//"/"//trim(logfile))
    else
       logUnit = 6
    endif

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogLevel(max(shrloglev,1))
    call shr_file_setLogUnit (logunit)

    !----------------------------------------------------------------------------
    ! Read input namelists and set present and prognostic flags
    !----------------------------------------------------------------------------

    call dwav_shr_read_namelists(mpicom, my_task, master_task, &
         inst_index, inst_suffix, inst_name,  &
         logunit, shrlogunit, SDWAV, wav_present, wav_prognostic)

    !--------------------------------
    ! advertise import and export fields
    !--------------------------------

    if (wav_present) then
       !-----------------
       ! export fields
       !-----------------

       call fld_list_add(fldsFrWav_num, fldsFrWav, trim(flds_scalar_name))
       call fld_list_add(fldsFrWav_num, fldsFrWav, 'Sw_lamult')
       call fld_list_add(fldsFrWav_num, fldsFrWav, 'Sw_ustokes')
       call fld_list_add(fldsFrWav_num, fldsFrWav, 'Sw_vstokes')

    end if

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    call shr_file_setLogLevel(shrloglev)
    call shr_file_setLogUnit (shrlogunit)

  end subroutine InitializeAdvertise

  !===============================================================================

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    character(ESMF_MAXSTR)   :: convCIM, purpComp
    type(ESMF_Grid)          :: Egrid
    type(ESMF_Mesh)          :: Emesh
    integer                  :: nx_global, ny_global
    type(ESMF_VM)            :: vm
    integer                  :: n
    character(CL)            :: cvalue
    integer                  :: shrlogunit                ! original log unit
    integer                  :: shrloglev                 ! original log level
    logical                  :: read_restart              ! start from restart
    integer                  :: ierr                      ! error code
    logical                  :: scmMode = .false.         ! single column mode
    real(R8)                 :: scmLat  = shr_const_SPVAL ! single column lat
    real(R8)                 :: scmLon  = shr_const_SPVAL ! single column lon
    logical                  :: connected                 ! is field connected?
    integer                  :: klon, klat
    integer                  :: lsize
    integer                  :: iam
    real(r8), pointer        :: lon(:),lat(:)
    integer , pointer        :: gindex(:)
    character(len=*),parameter :: subname=trim(modName)//':(InitializeRealize) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogLevel(max(shrloglev,1))
    call shr_file_setLogUnit (logUnit)

    !--------------------------------
    ! call dwav init routine
    !--------------------------------

    gsmap => gsmap_target
    ggrid => ggrid_target

    call NUOPC_CompAttributeGet(gcomp, name='read_restart', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) read_restart

    call NUOPC_CompAttributeGet(gcomp, name='MCTID', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) compid

    call dwav_comp_init(clock, x2d, d2x, &
         SDWAV, gsmap, ggrid, mpicom, compid, my_task, master_task, &
         inst_suffix, inst_name, logunit, read_restart, &
         flds_w2x, flds_x2w)

    !--------------------------------
    ! generate the grid or mesh from the gsmap and ggrid
    ! grid_option specifies grid or mesh
    !--------------------------------

    nx_global = SDWAV%nxg
    ny_global = SDWAV%nyg
    lsize = mct_gsMap_lsize(gsMap, mpicom)
    allocate(lon(lsize))
    allocate(lat(lsize))
    allocate(gindex(lsize))
    klat = mct_aVect_indexRA(ggrid%data, 'lat')
    klon = mct_aVect_indexRA(ggrid%data, 'lon')
    call mpi_comm_rank(mpicom, iam, ierr)
    call mct_gGrid_exportRattr(ggrid,'lon',lon,lsize)
    call mct_gGrid_exportRattr(ggrid,'lat',lat,lsize)
    call mct_gsMap_OrderedPoints(gsMap_target, iam, gindex)

    call shr_nuopc_grid_MeshInit(gcomp, nx_global, ny_global, mpicom, gindex, lon, lat, Emesh, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    deallocate(lon)
    deallocate(lat)
    deallocate(gindex)

    !--------------------------------
    ! realize the actively coupled fields, now that a grid or mesh is established
    ! Note: shr_nuopc_fldList_Realize does the following:
    ! 1) loops over all of the entries in fldsToWav and creates a field
    !    for each one via one of the following commands:
    !     field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R8, name=fldlist%shortname(n), rc=rc)
    !     field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=fldlist%shortname(n), meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    ! 2) realizes the field via the following call
    !     call NUOPC_Realize(state, field=field, rc=rc)
    !    where state is either importState or exportState
    !  NUOPC_Realize "realizes" a previously advertised field in the importState and exportState
    !  by replacing the advertised fields with the fields in fldsToWav of the same name.
    !--------------------------------

    call shr_nuopc_fldList_Realize(importState, fldListTo(compwav), flds_scalar_name, flds_scalar_num, &
         mesh=Emesh, tag=subname//':dwavImport', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_fldList_Realize(exportState, fldListFr(compwav), flds_scalar_name, flds_scalar_num, &
         mesh=Emesh, tag=subname//':dwavExport', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Pack export state
    ! Copy from d2x to exportState
    ! Set the coupling scalars
    !--------------------------------

    call shr_nuopc_grid_ArrayToState(d2x%rattr, flds_w2x, exportState, grid_option, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_State_SetScalar(dble(nx_global),flds_scalar_index_nx, exportState, mpicom, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_State_SetScalar(dble(ny_global),flds_scalar_index_ny, exportState, mpicom, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (dbug > 1) then
       if (my_task == master_task) then
          call mct_aVect_info(2, d2x, istr='initial diag'//':AV')
       end if
       call shr_nuopc_methods_State_diagnose(exportState,subname//':ES',rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

#ifdef USE_ESMF_METADATA
    convCIM  = "CIM"
    purpComp = "Model Component Simulation Description"
    call ESMF_AttributeAdd(comp, convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ShortName", "DWAV", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "LongName", "Climatological Wave Data Model", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "Description", &
         "The CIME data models perform the basic function of " // &
         "reading external data, modifying that data, and then " // &
         "sending it to the driver via coupling " // &
         "interfaces. The driver and other models have no " // &
         "fundamental knowledge of whether another component " // &
         "is fully active or just a data model.  In some cases, " // &
         "data models are prognostic and also receive and use " // &
         "some data sent by the driver to the data model.  But " // &
         "in most cases, the data models are not running " // &
         "prognostically and have no need to receive any data " // &
         "from the driver.", &
         convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ReleaseDate", "2010", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ModelType", "Wave", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "Name", "TBD", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "EmailAddress", "TBD", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ResponsiblePartyRole", "contact", convention=convCIM, purpose=purpComp, rc=rc)
#endif

    call shr_file_setLogLevel(shrloglev)
    call shr_file_setLogUnit (shrlogunit)

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine InitializeRealize

  !===============================================================================

  subroutine ModelAdvance(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)        :: clock
    type(ESMF_Alarm)        :: alarm
    type(ESMF_State)        :: importState, exportState
    integer                 :: shrlogunit    ! original log unit
    integer                 :: shrloglev     ! original log level
    character(CL)           :: case_name     ! case name
    logical                 :: write_restart ! write a restart
    integer                 :: currentYMD    ! model date
    integer                 :: currentTOD    ! model sec into model date
    integer                 :: nextYMD       ! model date
    integer                 :: nextTOD       ! model sec into model date
    type(ESMF_Time)         :: currTime, nextTime
    type(ESMF_TimeInterval) :: timeStep
    character(len=*),parameter :: subname=trim(modName)//':(ModelAdvance) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogLevel(max(shrloglev,1))
    call shr_file_setLogUnit (logunit)

    !--------------------------------
    ! query the Component for its clock, importState and exportState
    !--------------------------------

    call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug > 1) then
       if (my_task == master_task) then
          call shr_nuopc_methods_Clock_TimePrint(clock,subname//'clock',rc=rc)
       end if
    endif

    !--------------------------------
    ! Unpack import state
    !--------------------------------

    if (unpack_import) then
       call shr_nuopc_grid_StateToArray(importState, x2d%rattr, flds_x2w, grid_option, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !--------------------------------
    ! Run model
    !--------------------------------

    write_restart = seq_timemgr_alarmIsOn(clock, seq_timemgr_alarm_restart, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call seq_timemgr_AlarmSetOff(clock, seq_timemgr_alarm_restart, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! For nuopc - the component clock is advanced at the end of the time interval
    ! For these to match for now - need to advance nuopc one timestep ahead for
    ! shr_strdata time interpolation
    call ESMF_ClockGet( clock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    nextTime = currTime + timeStep
    call seq_timemgr_ETimeGet( currTime, ymd=CurrentYMD, tod=CurrentTOD )
    call seq_timemgr_ETimeGet( nextTime, ymd=NextYMD, tod=NextTOD )

    call dwav_comp_run(clock, x2d, d2x, &
       SDWAV, gsmap, ggrid, mpicom, compid, my_task, master_task, &
       inst_suffix, logunit, read_restart, write_restart, &
       nextYMD, nextTOD, case_name=case_name)

    !--------------------------------
    ! Pack export state
    !--------------------------------

    call shr_nuopc_grid_ArrayToState(d2x%rattr, flds_w2x, exportState, grid_option, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (dbug > 1) then
       if (my_task == master_task) then
          call mct_aVect_info(2, d2x, istr='initial diag'//':AV')
       end if
       call shr_nuopc_methods_State_diagnose(exportState,subname//':ES',rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    if (my_task == master_task) then
       call ESMF_ClockPrint(clock, options="currTime", preString="------>Advancing WAV from: ", rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_ClockPrint(clock, options="stopTime", preString="--------------------------------> to: ", rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

    call shr_file_setLogLevel(shrloglev)
    call shr_file_setLogUnit (shrlogunit)

  end subroutine ModelAdvance

  !===============================================================================

  subroutine ModelSetRunClock(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)         :: mclock, dclock
    type(ESMF_Time)          :: mcurrtime, dcurrtime
    type(ESMF_Time)          :: mstoptime
    type(ESMF_TimeInterval)  :: mtimestep, dtimestep
    character(len=128)       :: mtimestring, dtimestring
    type(ESMF_Alarm),pointer :: alarmList(:)
    type(ESMF_Alarm)         :: dalarm
    integer                  :: alarmcount, n
    character(len=*),parameter :: subname=trim(modName)//':(ModelSetRunClock) '

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(gcomp, driverClock=dclock, modelClock=mclock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(dclock, currTime=dcurrtime, timeStep=dtimestep, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGet(mclock, currTime=mcurrtime, timeStep=mtimestep, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! force model clock currtime and timestep to match driver and set stoptime
    !--------------------------------

    mstoptime = mcurrtime + dtimestep
    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! copy alarms from driver to model clock if model clock has no alarms (do this only once!)
    !--------------------------------

    call ESMF_ClockGetAlarmList(mclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmCount=alarmCount, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (alarmCount == 0) then
      call ESMF_ClockGetAlarmList(dclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmCount=alarmCount, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      allocate(alarmList(alarmCount))
      call ESMF_ClockGetAlarmList(dclock, alarmlistflag=ESMF_ALARMLIST_ALL, alarmList=alarmList, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

      do n = 1, alarmCount
         !call ESMF_AlarmPrint(alarmList(n), rc=rc)
         !if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
         dalarm = ESMF_AlarmCreate(alarmList(n), rc=rc)
         if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
         call ESMF_AlarmSet(dalarm, clock=mclock, rc=rc)
         if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      enddo

      deallocate(alarmList)
    endif

    !--------------------------------
    ! Advance model clock to trigger alarms then reset model clock back to currtime
    !--------------------------------

    call ESMF_ClockAdvance(mclock,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockSet(mclock, currTime=dcurrtime, timeStep=dtimestep, stopTime=mstoptime, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine ModelSetRunClock

  !===============================================================================

  subroutine ModelFinalize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(len=*),parameter  :: subname=trim(modName)//':(ModelFinalize) '

    !--------------------------------
    ! Finalize routine
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    call dwav_comp_final(my_task, master_task, logunit)

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine ModelFinalize

  !===============================================================================

end module wav_comp_nuopc
