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
  use med_constants_mod     , only : IN, R8, I8, CXX, CL, CS
  use med_constants_mod     , only : shr_log_Unit
  use med_constants_mod     , only : shr_file_getlogunit, shr_file_setlogunit
  use med_constants_mod     , only : shr_file_getloglevel, shr_file_setloglevel
  use med_constants_mod     , only : shr_file_setIO, shr_file_getUnit
  use shr_nuopc_scalars_mod , only : flds_scalar_name
  use shr_nuopc_scalars_mod , only : flds_scalar_num
  use shr_nuopc_scalars_mod , only : flds_scalar_index_nx
  use shr_nuopc_scalars_mod , only : flds_scalar_index_ny
  use shr_nuopc_scalars_mod , only : flds_scalar_index_nextsw_cday
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_Clock_TimePrint
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_SetScalar
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_Diagnose
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
  use shr_nuopc_grid_mod    , only : shr_nuopc_grid_Meshinit
  use dead_nuopc_mod        , only : dead_grid_lat, dead_grid_lon, dead_grid_index
  use dead_nuopc_mod        , only : dead_init_nuopc, dead_run_nuopc, dead_final_nuopc
  use dead_nuopc_mod        , only : fld_list_add, fld_list_realize, fldsMax, fld_list_type
  use dead_nuopc_mod        , only : state_getimport, state_setexport
  use dead_nuopc_mod        , only : ModelInitPhase, ModelSetRunClock, Print_FieldExchInfo
  use med_constants_mod, only : dbug => med_constants_dbug_flag

  implicit none
  private ! except

  public :: SetServices

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  integer                    :: fldsToAtm_num = 0
  integer                    :: fldsFrAtm_num = 0
  type (fld_list_type)       :: fldsToAtm(fldsMax)
  type (fld_list_type)       :: fldsFrAtm(fldsMax)
  real(r8), pointer          :: gbuf(:,:)            ! model info
  real(r8), pointer          :: lat(:)
  real(r8), pointer          :: lon(:)
  integer , allocatable      :: gindex(:)
  real(r8), allocatable      :: x2d(:,:)
  real(r8), allocatable      :: d2x(:,:)
  character(CXX)             :: flds_a2x = ''
  character(CXX)             :: flds_x2a = ''
  integer                    :: nxg                  ! global dim i-direction
  integer                    :: nyg                  ! global dim j-direction
  integer                    :: my_task              ! my task in mpi communicator
  integer                    :: inst_index           ! number of current instance (ie. 1)
  character(len=12)          :: inst_name            ! fullname of current instance (ie. "lnd_0001")
  character(len=5)          :: inst_suffix      ! char string associated with instance (ie. "_0001" or "")
  integer                    :: logunit              ! logging unit number
  logical :: mastertask
  logical                    :: atm_prognostic

  !----- formats -----
  character(*),parameter :: modName =  "(xatm_comp_nuopc)"
  character(*),parameter :: u_FILE_u = __FILE__


  !===============================================================================
  contains
  !===============================================================================

  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter  :: subname=trim(modName)//':(SetServices) '

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! the NUOPC gcomp component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=ModelInitPhase, phase=0, rc=rc)
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

    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

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
    type(ESMF_VM)      :: vm
    character(CL)      :: cvalue
    character(CS)      :: stdname
    integer            :: n
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
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=rc)

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localpet=my_task, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    mastertask = my_task==0
    !----------------------------------------------------------------------------
    ! determine instance information
    !----------------------------------------------------------------------------

    call shr_nuopc_get_component_instance(gcomp, inst_suffix, inst_index)
    inst_name = "ATM"//trim(inst_suffix)

    !----------------------------------------------------------------------------
    ! set logunit and set shr logging to my log file
    !----------------------------------------------------------------------------

    call shr_nuopc_set_component_logging(gcomp, mastertask, logunit, shrlogunit, shrloglev)

    !----------------------------------------------------------------------------
    ! Initialize xatm
    !----------------------------------------------------------------------------

    call dead_init_nuopc('atm', inst_suffix, logunit, lsize, gbuf, nxg, nyg)

    allocate(gindex(lsize))
    allocate(lon(lsize))
    allocate(lat(lsize))

    gindex(:) = gbuf(:,dead_grid_index)
    lat(:)    = gbuf(:,dead_grid_lat)
    lon(:)    = gbuf(:,dead_grid_lon)

    !--------------------------------
    ! advertise import and export fields
    !--------------------------------

    if (nxg /= 0 .and. nyg /= 0) then

       call fld_list_add(fldsFrAtm_num, fldsFrAtm, trim(flds_scalar_name))
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_topo'       , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_z'          , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_u'          , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_v'          , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_tbot'       , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_ptem'       , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_shum'       , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_pbot'       , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_dens'       , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Sa_pslv'       , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_rainc'    , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_rainl'    , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_snowc'    , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_snowl'    , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_lwdn'     , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_swndr'    , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_swvdr'    , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_swndf'    , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_swvdf'    , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_swnet'    , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_bcphidry' , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_bcphodry' , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_bcphiwet' , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_ocphidry' , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_ocphodry' , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_ocphiwet' , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_dstwet1'  , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_dstwet2'  , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_dstwet3'  , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_dstwet4'  , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_dstdry1'  , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_dstdry2'  , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_dstdry3'  , flds_concat=flds_a2x)
       call fld_list_add(fldsFrAtm_num, fldsFrAtm, 'Faxa_dstdry4'  , flds_concat=flds_a2x)

       do n = 1,fldsFrAtm_num
          if(mastertask) write(logunit,*)'Advertising From Xatm ',trim(fldsFrAtm(n)%stdname)
          call NUOPC_Advertise(exportState, standardName=fldsFrAtm(n)%stdname, &
               TransferOfferGeomObject='will provide', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end do

       call fld_list_add(fldsToAtm_num, fldsToAtm, trim(flds_scalar_name))
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sx_anidr'  , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sx_avsdf'  , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sx_anidf'  , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sx_avsdr'  , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sl_lfrac'  , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Si_ifrac'  , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'So_ofrac'  , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sx_tref'   , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sx_qref'   , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sx_t'      , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'So_t'      , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sl_fv'     , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sl_ram1'   , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sl_snowh'  , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Si_snowh'  , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'So_ssq'    , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'So_re'     , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Sx_u10'    , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Faxx_taux' , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Faxx_tauy' , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Faxx_lat'  , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Faxx_sen'  , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Faxx_lwup' , flds_concat=flds_x2a)
       call fld_list_add(fldsToAtm_num, fldsToAtm, 'Faxx_evap' , flds_concat=flds_x2a)

       do n = 1,fldsToAtm_num
          if(mastertask) write(logunit,*)'Advertising To Xatm',trim(fldsToAtm(n)%stdname)
          call NUOPC_Advertise(importState, standardName=fldsToAtm(n)%stdname, &
               TransferOfferGeomObject='will provide', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       enddo

       allocate(d2x(FldsFrAtm_num,lsize)); d2x(:,:)  = 0._r8
       allocate(x2d(FldsToAtm_num,lsize)); x2d(:,:)  = 0._r8
    end if

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    call shr_file_setLogLevel(shrloglev)
    call shr_file_setLogUnit (shrlogunit)

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=rc)

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
    integer                :: n
    integer                :: shrlogunit                ! original log unit
    integer                :: shrloglev                 ! original log level
    character(len=*),parameter :: subname=trim(modName)//':(InitializeRealize: xatm) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=rc)

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

    call shr_nuopc_grid_MeshInit(gcomp, nxg, nyg, gindex, lon, lat, Emesh, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! realize the actively coupled fields, now that a mesh is established
    ! NUOPC_Realize "realizes" a previously advertised field in the importState and exportState
    ! by replacing the advertised fields with the newly created fields of the same name.
    !--------------------------------

    call fld_list_realize( &
         state=exportState, &
         fldList=fldsFrAtm, &
         numflds=fldsFrAtm_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':datmExport',&
         mesh=Emesh, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call fld_list_realize( &
         state=importState, &
         fldList=fldsToAtm, &
         numflds=fldsToAtm_num, &
         flds_scalar_name=flds_scalar_name, &
         flds_scalar_num=flds_scalar_num, &
         tag=subname//':datmImport',&
         mesh=Emesh, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Pack export state
    ! Copy from d2x to exportState
    ! Set the coupling scalars
    !--------------------------------

    do n = 1, FldsFrAtm_num
       if (fldsFrAtm(n)%stdname /= flds_scalar_name) then
          call state_setexport(exportState, trim(fldsFrAtm(n)%stdname), d2x(n,:), rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end do

    call shr_nuopc_methods_State_SetScalar(dble(nxg),flds_scalar_index_nx, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_State_SetScalar(dble(nyg),flds_scalar_index_ny, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set time of next radiation computation

    call ESMF_ClockGetNextTime(clock, nextTime)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(nextTime, dayOfYear_r8=nextsw_cday)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_State_SetScalar(nextsw_cday, flds_scalar_index_nextsw_cday, exportState, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (dbug > 1) then
       if (mastertask) then
          call Print_FieldExchInfo(values=d2x, logunit=logunit, &
               fldlist=fldsFrAtm, nflds=fldsFrAtm_num, istr="InitializeRealize: atm->mediator")
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

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=rc)

  end subroutine InitializeRealize

  !===============================================================================

  subroutine ModelAdvance(gcomp, rc)
    use shr_nuopc_utils_mod, only : shr_nuopc_memcheck, shr_nuopc_log_clock_advance
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock) :: clock
    type(ESMF_Time)  :: nexttime
    type(ESMF_State) :: importState, exportState
    real(r8)         :: nextsw_cday
    integer          :: n
    integer          :: shrlogunit ! original log unit
    integer          :: shrloglev  ! original log level
    character(len=CL)      :: clockstr
    character(len=*),parameter  :: subname=trim(modName)//':(ModelAdvance) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=rc)
    call shr_nuopc_memcheck(subname, 3, mastertask)
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
    ! Unpack import state
    !--------------------------------

    do n = 1, FldsFrAtm_num
       if (fldsFrAtm(n)%stdname /= flds_scalar_name) then
          call state_getimport(importState, trim(fldsToAtm(n)%stdname), x2d(n,:), rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end do

    !--------------------------------
    ! Run model
    !--------------------------------

    call dead_run_nuopc('atm', d2x, gbuf, flds_a2x)

    !--------------------------------
    ! Pack export state
    !--------------------------------

    do n = 1, FldsFrAtm_num
       if (fldsFrAtm(n)%stdname /= flds_scalar_name) then
          call state_setexport(exportState, trim(fldsFrAtm(n)%stdname), d2x(n,:), rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end do

    call ESMF_ClockGetNextTime(clock, nextTime)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(nextTime, dayOfYear_r8=nextsw_cday)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_State_SetScalar(nextsw_cday, flds_scalar_index_nextsw_cday, exportState, &
          flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (dbug > 1) then
       if (mastertask) then
          call Print_FieldExchInfo(values=d2x, logunit=logunit, &
               fldlist=fldsFrAtm, nflds=fldsFrAtm_num, istr="ModelAdvance: atm->mediator")
       end if
       call shr_nuopc_methods_State_diagnose(exportState,subname//':ES',rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif
    if(mastertask) then
       call shr_nuopc_log_clock_advance(clock, 'ATM', logunit)
    endif
    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=rc)

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
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=rc)

    call dead_final_nuopc('atm', logunit)

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=rc)

  end subroutine ModelFinalize

end module atm_comp_nuopc
