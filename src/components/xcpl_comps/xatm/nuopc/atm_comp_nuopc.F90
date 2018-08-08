module atm_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for XATM
  !----------------------------------------------------------------------------

  use ESMF
  use NUOPC                 , only : NUOPC_CompDerive, NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize
  use NUOPC                 , only : NUOPC_CompAttributeGet, NUOPC_Advertise
  use NUOPC_Model           , only : model_routine_SS        => SetServices
  use NUOPC_Model           , only : model_label_Advance     => label_Advance
  use NUOPC_Model           , only : model_label_SetRunClock => label_SetRunClock
  use NUOPC_Model           , only : model_label_Finalize    => label_Finalize
  use NUOPC_Model           , only : NUOPC_ModelGet
  use med_constants_mod     , only : IN, R8, I8, CXX
  use med_constants_mod     , only : shr_log_Unit
  use med_constants_mod     , only : shr_file_getlogunit, shr_file_setlogunit
  use med_constants_mod     , only : shr_file_getloglevel, shr_file_setloglevel
  use med_constants_mod     , only : shr_file_setIO, shr_file_getUnit
  use shr_string_mod        , only : shr_string_listGetNum
  use shr_nuopc_scalars_mod , only : flds_scalar_name
  use shr_nuopc_scalars_mod , only : flds_scalar_num
  use shr_nuopc_scalars_mod , only : flds_scalar_index_nx
  use shr_nuopc_scalars_mod , only : flds_scalar_index_ny
  use shr_nuopc_scalars_mod , only : flds_scalar_index_nextsw_cday
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_Realize
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_Concat
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_Deactivate
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_Getnumflds
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_Getfldinfo
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_Clock_TimePrint
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_SetScalar
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_Diagnose
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_Print_FieldExchInfo
  use shr_nuopc_grid_mod    , only : shr_nuopc_grid_Meshinit
  use shr_nuopc_grid_mod    , only : shr_nuopc_grid_ArrayToState
  use shr_nuopc_grid_mod    , only : shr_nuopc_grid_StateToArray

  ! TODO: remove this
  use esmFlds               , only : fldListFr, fldListTo, compatm, compname

  use dead_data_mod , only : dead_grid_lat, dead_grid_lon, dead_grid_index
  use dead_nuopc_mod, only : dead_init_nuopc, dead_run_nuopc, dead_final_nuopc
  use dead_nuopc_mod, only : ModelSetRunClock

  implicit none
  private ! except

  public :: SetServices

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  real(r8), pointer          :: gbuf(:,:)            ! model info
  real(r8), pointer          :: lat(:)
  real(r8), pointer          :: lon(:)
  integer , allocatable      :: gindex(:)
  real(r8), allocatable      :: x2d(:,:)
  real(r8), allocatable      :: d2x(:,:)
  integer                    :: nflds_d2x
  integer                    :: nflds_x2d
  integer                    :: nxg                  ! global dim i-direction
  integer                    :: nyg                  ! global dim j-direction
  integer                    :: mpicom               ! mpi communicator
  integer                    :: my_task              ! my task in mpi communicator mpicom
  integer                    :: inst_index           ! number of current instance (ie. 1)
  character(len=16)          :: inst_name            ! fullname of current instance (ie. "lnd_0001")
  character(len=16)          :: inst_suffix = ""     ! char string associated with instance (ie. "_0001" or "")
  integer                    :: logunit              ! logging unit number
  integer    ,parameter      :: master_task=0        ! task number of master task
  character(len=*),parameter :: grid_option = "mesh" ! grid_de, grid_arb, grid_reg, mesh
  integer, parameter         :: dbug = 10
  integer                    :: dbrc
  logical                    :: atm_prognostic
  character(CXX)             :: flds_a2x = ''
  character(CXX)             :: flds_x2a = ''

  !----- formats -----
  character(*),parameter :: modName =  "(atm_comp_nuopc)"
  character(*),parameter :: u_FILE_u = __FILE__

  !===============================================================================
  contains
  !===============================================================================

  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter  :: subname=trim(modName)//':(SetServices) '

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    ! the NUOPC gcomp component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      userRoutine=InitializeP0, phase=0, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, phaseLabelList=(/"IPDv01p1"/), &
         userRoutine=InitializeAdvertise, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, phaseLabelList=(/"IPDv01p3"/), &
         userRoutine=InitializeRealize, rc=rc)
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

  !-----------------------------------------------------------------------------

  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries

    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, acceptStringList=(/"IPDv01p"/), rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine InitializeP0

  !===============================================================================

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_VM)      :: vm
    integer            :: lmpicom
    character(CL)      :: cvalue
    character(CS)      :: stdname, shortname
    logical            :: activefld
    integer            :: n,nflds
    integer            :: lsize       ! local array size
    integer            :: ierr        ! error code
    integer            :: shrlogunit  ! original log unit
    integer            :: shrloglev   ! original log level
    logical            :: isPresent
    character(len=512) :: diro
    character(len=512) :: logfile
    character(len=*),parameter :: subname=trim(modName)//':(InitializeAdvertise) '
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

    !--------------------------------
    ! create import and export field list
    !--------------------------------

    call shr_nuopc_fldList_Concat(fldListFr(compatm), fldListTo(compatm), flds_a2x, flds_x2a, flds_scalar_name)

    !----------------------------------------------------------------------------
    ! Initialize xatm
    !----------------------------------------------------------------------------

    call dead_init_nuopc('atm', mpicom, my_task, master_task, &
         inst_index, inst_suffix, inst_name, logunit, lsize, gbuf, nxg, nyg)

    nflds_d2x = shr_string_listGetNum(flds_a2x)
    nflds_x2d = shr_string_listGetNum(flds_x2a)

    allocate(gindex(lsize))
    allocate(lon(lsize))
    allocate(lat(lsize))
    allocate(d2x(nflds_d2x,lsize))
    allocate(x2d(nflds_x2d,lsize))

    gindex(:) = gbuf(:,dead_grid_index)
    lat(:)    = gbuf(:,dead_grid_lat)
    lon(:)    = gbuf(:,dead_grid_lon)
    d2x(:,:)  = 0._r8
    x2d(:,:)  = 0._r8

    if (nxg == 0 .and. nyg == 0) then
       atm_prognostic=.false.
    else
       atm_prognostic=.true.
    end if

    !--------------------------------
    ! advertise import and export fields
    !--------------------------------

    ! First deactivate fldListTo(compatm) if atm_prognostic is .false.
    if (.not. atm_prognostic) then
       call shr_nuopc_fldList_Deactivate(fldListTo(compatm), flds_scalar_name)
    end if

    nflds = shr_nuopc_fldList_Getnumflds(fldListFr(compatm))
    do n = 1,nflds
       call shr_nuopc_fldList_Getfldinfo(fldListFr(compatm), n, activefld, stdname, shortname)
       if (activefld) then
          call NUOPC_Advertise(exportState, standardName=stdname, shortname=shortname, name=shortname, &
               TransferOfferGeomObject='will provide', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_LogWrite(subname//':Fr_'//trim(compname(compatm))//': '//trim(shortname), ESMF_LOGMSG_INFO)
       end if
    end do

    nflds = shr_nuopc_fldList_Getnumflds(fldListTo(compatm))
    do n = 1,nflds
       call shr_nuopc_fldList_Getfldinfo(fldListTo(compatm), n, activefld, stdname, shortname)
       if (activefld) then
          call NUOPC_Advertise(importState, standardName=stdname, shortname=shortname, name=shortname, &
               TransferOfferGeomObject='will provide', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_LogWrite(subname//':To_'//trim(compname(compatm))//': '//trim(shortname), ESMF_LOGMSG_INFO)
       end if
    end do

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
    character(ESMF_MAXSTR) :: convCIM, purpComp
    type(ESMF_Mesh)        :: Emesh
    type(ESMF_Time)        :: nextTime
    real(r8)               :: nextsw_cday
    integer                :: shrlogunit                ! original log unit
    integer                :: shrloglev                 ! original log level
    type(ESMF_VM)          :: vm
    integer                :: n
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
    ! generate the mesh
    !--------------------------------

    call shr_nuopc_grid_MeshInit(gcomp, nxg, nyg, mpicom, gindex, lon, lat, Emesh, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! realize the actively coupled fields
    !--------------------------------

    call shr_nuopc_fldList_Realize(importState, fldListTo(compatm), flds_scalar_name, flds_scalar_num, &
         mesh=Emesh, tag=subname//':xatmImport', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_fldList_Realize(exportState, fldListFr(compatm), flds_scalar_name, flds_scalar_num, &
         mesh=Emesh, tag=subname//':xatmExport', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Pack export state
    ! Copy from d2x to exportState
    ! Set the coupling scalars
    !--------------------------------

    call shr_nuopc_grid_ArrayToState(d2x, flds_a2x, exportState, grid_option, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_State_SetScalar(dble(nxg),flds_scalar_index_nx, exportState, mpicom, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_State_SetScalar(dble(nyg),flds_scalar_index_ny, exportState, mpicom, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set time of next radiation computation

    call ESMF_ClockGetNextTime(clock, nextTime)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(nextTime, dayOfYear_r8=nextsw_cday)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_State_SetScalar(nextsw_cday, flds_scalar_index_nextsw_cday, exportState, mpicom, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (dbug > 1) then
       if (my_task == master_task) then
          call shr_nuopc_methods_Print_FieldExchInfo(flag=2, values=d2x, logunit=logunit, &
               fldlist=flds_a2x, nflds=nflds_d2x, istr="InitializeRealize: atm->mediator")
       end if
       call shr_nuopc_methods_State_diagnose(exportState,subname//':ES',rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

#ifdef USE_ESMF_METADATA
    convCIM  = "CIM"
    purpComp = "Model Component Simulation Description"
    call ESMF_AttributeAdd(comp,  convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ShortName", "XATM", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "LongName", "Atmosphere Dead Model", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "Description", &
         "The dead models stand in as test model for active components." // &
         "Coupling data is artificially generated ", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ReleaseDate", "2017", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ModelType", "Sea Ice", convention=convCIM, purpose=purpComp, rc=rc)
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
    type(ESMF_Clock) :: clock
    type(ESMF_Time)  :: nexttime
    type(ESMF_State) :: importState, exportState
    integer          :: CurrentYMD, CurrentTOD
    real(r8)         :: nextsw_cday
    integer          :: shrlogunit ! original log unit
    integer          :: shrloglev  ! original log level
    integer(IN)      :: n          ! index
    integer(IN)      :: nf         ! fields loop index
    integer(IN)      :: ki         ! index
    integer(IN)      :: lsize      ! size of AttrVect
    real(R8)         :: lat        ! latitude
    real(R8)         :: lon        ! longitude
    integer          :: nflds_x2d
    integer          :: nflds_d2x
    integer          :: ncomp
    character(len=*),parameter  :: subname=trim(modName)//':(ModelAdvance) '
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
      call shr_nuopc_methods_Clock_TimePrint(clock,subname//'clock',rc=rc)
    endif

    !--------------------------------
    ! unpack import state
    !--------------------------------

    call shr_nuopc_grid_StateToArray(importState, x2d, flds_x2a, grid_option, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Run model
    !--------------------------------

    call dead_run_nuopc('atm', d2x, gbuf, flds_d2x)

    !--------------------------------
    ! Pack export state
    !--------------------------------

    call shr_nuopc_grid_ArrayToState(d2x, flds_a2x, exportState, grid_option, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockGetNextTime(clock, nextTime)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(nextTime, dayOfYear_r8=nextsw_cday)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_State_SetScalar(nextsw_cday, flds_scalar_index_nextsw_cday, exportState, mpicom, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (dbug > 1) then
       if (my_task == master_task) then
          call shr_nuopc_methods_Print_FieldExchInfo(flag=2, values=d2x, logunit=logunit, &
               fldlist=flds_a2x, nflds=nflds_d2x, istr="ModelAdvance: atm->mediator")
       end if
       call shr_nuopc_methods_State_diagnose(exportState,subname//':ES',rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    call ESMF_ClockPrint(clock, options="currTime", preString="------>Advancing ATM from: ", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockPrint(clock, options="stopTime", preString="--------------------------------> to: ", rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_file_setLogLevel(shrloglev)
    call shr_file_setLogUnit (shrlogunit)

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine ModelAdvance

  !===============================================================================

  subroutine ModelFinalize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(len=*),parameter  :: subname=trim(modName)//':(ModelFinalize) '
    !-------------------------------------------------------------------------------

    !--------------------------------
    ! Finalize routine
    !--------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

    call dead_final_nuopc('atm', my_task, master_task, logunit)

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine ModelFinalize

end module atm_comp_nuopc
