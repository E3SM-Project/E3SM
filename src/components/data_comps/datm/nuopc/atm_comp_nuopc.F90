module atm_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for DATM
  !----------------------------------------------------------------------------

  use shr_kind_mod          , only : CXX => shr_kind_CXX
  use med_constants_mod     , only : IN, R8, I8
  use shr_log_mod           , only : shr_log_Unit
  use shr_cal_mod           , only : shr_cal_ymd2date, shr_cal_noleap, shr_cal_gregorian
  use shr_file_mod          , only : shr_file_getlogunit, shr_file_setlogunit
  use shr_file_mod          , only : shr_file_getloglevel, shr_file_setloglevel
  use shr_file_mod          , only : shr_file_setIO, shr_file_getUnit
  use esmFlds               , only : fldListFr, fldListTo, compatm, compname
  use esmFlds               , only : flds_scalar_name
  use esmFlds               , only : flds_scalar_num
  use esmFlds               , only : flds_scalar_index_nx
  use esmFlds               , only : flds_scalar_index_ny
  use esmFlds               , only : flds_scalar_index_nextsw_cday
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_Realize
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_Concat
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_Deactivate
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_Getnumflds
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldList_Getfldinfo
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_Clock_TimePrint
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_SetScalar
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_Diagnose
  use shr_nuopc_grid_mod    , only : shr_nuopc_grid_Meshinit
  use shr_nuopc_grid_mod    , only : shr_nuopc_grid_ArrayToState
  use shr_nuopc_grid_mod    , only : shr_nuopc_grid_StateToArray
  use shr_nuopc_time_mod    , only : shr_nuopc_time_alarmInit
  use shr_strdata_mod       , only : shr_strdata_type

  use ESMF
  use NUOPC
  use NUOPC_Model, &
    model_routine_SS      => SetServices, &
    model_label_Advance   => label_Advance, &
    model_label_SetRunClock => label_SetRunClock, &
    model_label_Finalize  => label_Finalize

  use datm_shr_mod , only: datm_shr_read_namelists
  use datm_shr_mod , only: iradsw         ! namelist input
  use datm_shr_mod , only: presaero       ! namelist input
  use datm_shr_mod , only: datm_shr_getNextRadCDay
  use datm_comp_mod, only: datm_comp_init, datm_comp_run, datm_comp_final
  use mct_mod

  implicit none
  private ! except

  public  :: SetServices

  private :: InitializeP0
  private :: InitializeAdvertise
  private :: InitializeRealize
  private :: ModelAdvance
  private :: ModelSetRunClock
  private :: ModelFinalize

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  character(len=80)          :: myModelName = 'atm'       ! user defined model name
  type(shr_strdata_type)     :: SDATM
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
  character(len=16)          :: inst_name                 ! fullname of current instance (ie. "lnd_0001")
  character(len=16)          :: inst_suffix = ""          ! char string associated with instance (ie. "_0001" or "")
  integer                    :: logunit                   ! logging unit number
  integer    ,parameter      :: master_task=0             ! task number of master task
  integer                    :: localPet
  logical                    :: atm_present               ! flag
  logical                    :: atm_prognostic            ! flag
  logical                    :: unpack_import
  character(len=256)         :: case_name                 ! case name
  character(len=256)         :: tmpstr                    ! tmp string
  integer                    :: dbrc
  integer, parameter         :: dbug = 10
  character(len=*),parameter :: grid_option = "mesh"      ! grid_de, grid_arb, grid_reg, mesh
  character(len=80)          :: calendar                  ! calendar name
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

    ! local variables
    character(len=*),parameter  :: subname=trim(modName)//':(ATMSetServices) '
    !-------------------------------------------------------------------------------

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

    !----------------------------------------------------------------------------
    ! Switch to IPDv01 by filtering all other phaseMap entries
    !----------------------------------------------------------------------------

    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
         acceptStringList=(/"IPDv01p"/), rc=rc)
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
    character(len=256) :: cvalue
    character(len=80)  :: stdname, shortname
    logical            :: activefld
    integer            :: n,nflds
    integer            :: ierr       ! error code
    integer            :: shrlogunit ! original log unit
    integer            :: shrloglev  ! original log level
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

    call ESMF_VMGet(vm, mpiCommunicator=lmpicom, localPet=localPet, rc=rc)
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

    call datm_shr_read_namelists(mpicom, my_task, master_task, &
         inst_index, inst_suffix, inst_name, &
         logunit, shrlogunit, SDATM, atm_present, atm_prognostic)

    ! NOTE: atm_present flag is not needed - since the run sequence
    ! will have no call to this routine for the atm_present flag being
    ! set to false (i.e. null mode) - only the atm_prognostic flag is
    ! needed below

    if (atm_prognostic) then
       unpack_import = .true.
    else
       unpack_import = .false.
    end if

    !--------------------------------
    ! create import and export field list needed by data models
    !--------------------------------

    call shr_nuopc_fldList_Concat(fldListFr(compatm), fldListTo(compatm), flds_a2x, flds_x2a, flds_scalar_name)

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
    character(ESMF_MAXSTR)  :: convCIM, purpComp
    type(ESMF_Grid)         :: Egrid
    type(ESMF_TIME)         :: currTime
    type(ESMF_TimeInterval) :: timeStep
    type(ESMF_Mesh)         :: Emesh
    type(ESMF_Calendar)     :: esmf_calendar             ! esmf calendar
    type(ESMF_CalKind_Flag) :: esmf_caltype              ! esmf calendar type
    type(ESMF_VM)           :: vm
    integer                 :: nx_global, ny_global
    integer                 :: n
    character(len=256)      :: cvalue
    integer                 :: shrlogunit                ! original log unit
    integer                 :: shrloglev                 ! original log level
    logical                 :: read_restart              ! start from restart
    integer                 :: ierr                      ! error code
    logical                 :: scmMode = .false.         ! single column mode
    real(R8)                :: scmLat  = shr_const_SPVAL ! single column lat
    real(R8)                :: scmLon  = shr_const_SPVAL ! single column lon
    integer                 :: current_ymd               ! model date
    integer                 :: current_year              ! model year
    integer                 :: current_mon               ! model month
    integer                 :: current_day               ! model day
    integer                 :: current_tod               ! model sec into model date
    integer(I8)             :: stepno                    ! step number
    integer                 :: modeldt                   ! integer timestep
    real(R8)                :: nextsw_cday               ! calendar of next atm sw
    logical                 :: connected                 ! is field connected?
    integer                 :: klon, klat
    integer                 :: lsize
    integer                 :: iam
    real(r8), pointer       :: lon(:),lat(:)
    integer , pointer       :: gindex(:)
    real(R8)                :: orbEccen                  ! orb eccentricity (unit-less)
    real(R8)                :: orbMvelpp                 ! orb moving vernal eq (radians)
    real(R8)                :: orbLambm0                 ! orb mean long of perhelion (radians)
    real(R8)                :: orbObliqr                 ! orb obliquity (radians)
    character(len=*) , parameter :: subname=trim(modName)//':(InitializeRealize) '
    !-------------------------------------------------------------------------------

    ! TODO: read_restart, scmlat, scmlon, orbeccen, orbmvelpp, orblambm0, orbobliqr needs to be obtained
    ! from the config attributes of the gridded component

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
    ! call datm init routine
    !--------------------------------

    gsmap => gsmap_target
    ggrid => ggrid_target

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

    call NUOPC_CompAttributeGet(gcomp, name='MCTID', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) compid

    !----------------------------------------------------------------------------
    ! Determine orbital values (these might change dynamically)
    !----------------------------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='orb_eccen', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orbEccen
    call NUOPC_CompAttributeGet(gcomp, name='orb_obliqr', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orbObliqr
    call NUOPC_CompAttributeGet(gcomp, name='orb_lambm0', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orbLambm0
    call NUOPC_CompAttributeGet(gcomp, name='orb_mvelpp', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orbMvelpp

    !----------------------------------------------------------------------------
    ! Determine calendar info
    !----------------------------------------------------------------------------

    call ESMF_ClockGet( clock, currTime=currTime, timeStep=timeStep, advanceCount=stepno, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimeGet( currTime, yy=current_year, mm=current_mon, dd=current_day, s=current_tod, &
         calkindflag=esmf_caltype, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(current_year, current_mon, current_day, current_ymd)

    if (esmf_caltype == ESMF_CALKIND_NOLEAP) then
       calendar = shr_cal_noleap
    else if (esmf_caltype == ESMF_CALKIND_GREGORIAN) then
       calendar = shr_cal_gregorian
    else
       call ESMF_LogWrite(subname//" ERROR bad ESMF calendar name "//trim(calendar), ESMF_LOGMSG_ERROR, rc=dbrc)
       rc = ESMF_Failure
       return
    end if

    call ESMF_TimeIntervalGet( timeStep, s=modeldt, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call datm_comp_init(clock, x2d, d2x, &
         flds_x2a, flds_a2x, &
         SDATM, gsmap, ggrid, mpicom, compid, my_task, master_task, &
         inst_suffix, inst_name, logunit, read_restart, &
         scmMode, scmlat, scmlon, &
         orbEccen, orbMvelpp, orbLambm0, orbObliqr, &
         calendar, modeldt, current_ymd, current_tod, current_mon)

    !----------------------------------------------------------------------------
    ! Set nextsw_cday
    !----------------------------------------------------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='read_restart', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) read_restart

    if (read_restart) then
       nextsw_cday = datm_shr_getNextRadCDay( current_ymd, current_tod, stepno, modeldt, iradsw, calendar )
    else
       ! For a startup run the nextsw_cday is just the current calendar day
       call ESMF_TimeGet( currTime, dayofyear_r8=nextsw_cday, rc=rc )
    endif

    !--------------------------------
    ! Generate the mesh
    !--------------------------------

    nx_global = SDATM%nxg
    ny_global = SDATM%nyg
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
    ! realize the actively coupled fields, now that a mesh is established
    ! Note: shr_nuopc_fldList_Realize does the following:
    ! 1) loops over all of the entries in fldListTo and fldListFr creates a field
    !    for each one via one of the following commands:
    !     field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R8, name=shortname, rc=rc)
    !     field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=shortname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    ! 2) realizes the field via the following call
    !     call NUOPC_Realize(state, field=field, rc=rc)
    !    where state is either importState or exportState
    !  NUOPC_Realize "realizes" a previously advertised field in the importState and exportState
    !  by replacing the advertised fields with the newly created fields of the same name.
    !--------------------------------

    call shr_nuopc_fldList_Realize(importState, fldListTo(compatm), flds_scalar_name, flds_scalar_num, &
         mesh=Emesh, tag=subname//':datmImport', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_fldList_Realize(exportState, fldListFr(compatm), flds_scalar_name, flds_scalar_num, &
         mesh=Emesh, tag=subname//':datmExport', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Pack export state
    ! Copy from d2x to exportState
    ! Set the coupling scalars
    !--------------------------------

    call shr_nuopc_grid_ArrayToState(d2x%rattr, flds_a2x, exportState, grid_option, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !TODO: do we still need the nx and ny scalars?
    call shr_nuopc_methods_State_SetScalar(dble(nx_global),flds_scalar_index_nx, exportState, mpicom, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_State_SetScalar(dble(ny_global),flds_scalar_index_ny, exportState, mpicom, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_State_SetScalar(nextsw_cday, flds_scalar_index_nextsw_cday, exportState, mpicom, &
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
    call ESMF_AttributeSet(comp, "ShortName", "DATM", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "LongName", "Climatological Atmosphere Data Model", convention=convCIM, purpose=purpComp, rc=rc)
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
    call ESMF_AttributeSet(comp, "ModelType", "Atmosphere", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "Name", "TBD", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "EmailAddress", "TBD", convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ResponsiblePartyRole", "contact", convention=convCIM, purpose=purpComp, rc=rc)
#endif

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------

    call shr_file_setLogLevel(shrloglev)
    call shr_file_setLogUnit (shrlogunit)


  end subroutine InitializeRealize

  !===============================================================================

  subroutine ModelAdvance(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)        :: clock
    type(ESMF_Time)         :: time
    type(ESMF_State)        :: importState, exportState
    type(ESMF_Alarm)        :: alarm
    integer                 :: shrlogunit    ! original log unit
    integer                 :: shrloglev     ! original log level
    real(r8)                :: nextsw_cday
    logical                 :: write_restart ! restart alarm is ringing
    integer                 :: nextymd       ! model date
    integer                 :: nexttod       ! model sec into model date
    integer                 :: yr            ! year
    integer                 :: mon           ! month
    integer                 :: day           ! day in month
    type(ESMF_Time)         :: currTime
    type(ESMF_Time)         :: nextTime
    type(ESMF_TimeInterval) :: timeStep
    integer(I8)             :: stepno        ! step number
    integer                 :: modeldt       ! model timestep
    real(R8)                :: orbEccen      ! orb eccentricity (unit-less)
    real(R8)                :: orbMvelpp     ! orb moving vernal eq (radians)
    real(R8)                :: orbLambm0     ! orb mean long of perhelion (radians)
    real(R8)                :: orbObliqr     ! orb obliquity (radians)
    character(len=256)      :: cvalue
    character(len=*),parameter  :: subname=trim(modName)//':(ModelAdvance) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (dbug > 5) call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)

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

    if (dbug > 1) then
       if (my_task == master_task) then
          call shr_nuopc_methods_Clock_TimePrint(clock,subname//'clock',rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    endif

    !--------------------------------
    ! Unpack import state
    !--------------------------------

    if (unpack_import) then
       call shr_nuopc_grid_StateToArray(importState, x2d%rattr, flds_x2a, grid_option, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !--------------------------------
    ! Run model
    !--------------------------------

    ! Get orbital parameters (these can be changed by the mediator)
    ! TODO: need to put in capability for these to be modified for variable orbitals
    call NUOPC_CompAttributeGet(gcomp, name='orb_eccen', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orbEccen
    call NUOPC_CompAttributeGet(gcomp, name='orb_obliqr', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orbObliqr
    call NUOPC_CompAttributeGet(gcomp, name='orb_lambm0', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orbLambm0
    call NUOPC_CompAttributeGet(gcomp, name='orb_mvelpp', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orbMvelpp

    ! Determine if need to write restarts

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

    call ESMF_ClockGet( clock, currTime=currTime, timeStep=timeStep, advanceCount=stepno, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    nextTime = currTime + timeStep
    call ESMF_TimeGet( nextTime, yy=yr, mm=mon, dd=day, s=nexttod, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yr, mon, day, nextymd)

    call ESMF_TimeIntervalGet( timeStep, s=modeldt, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Advance the model

    call datm_comp_run(clock, &
         x2a=x2d, &
         a2x=d2x, &
         SDATM=SDATM, &
         gsmap=gsmap, &
         ggrid=ggrid, &
         mpicom=mpicom, &
         compid=compid, &
         my_task=my_task, &
         master_task=master_task, &
         inst_suffix=inst_suffix, &
         logunit=logunit, &
         orbEccen=orbEccen, &
         orbMvelpp=orbMvelpp, &
         orbLambm0=orbLambm0, &
         orbObliqr=orbObliqr, &
         write_restart=write_restart, &
         target_ymd=nextYMD, &
         target_tod=nextTOD, &
         target_mon=mon, &
         calendar=calendar, &
         modeldt=modeldt, &
         case_name=case_name)

    ! Use nextYMD and nextTOD here since since the component - clock is advance at the END of the time interval
    nextsw_cday = datm_shr_getNextRadCDay( nextYMD, nextTOD, stepno, modeldt, iradsw, calendar )

    !--------------------------------
    ! Pack export state
    !--------------------------------

    call shr_nuopc_grid_ArrayToState(d2x%rattr, flds_a2x, exportState, grid_option, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_State_SetScalar(nextsw_cday, flds_scalar_index_nextsw_cday, exportState, mpicom, &
         flds_scalar_name, flds_scalar_num, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! diagnostics
    !--------------------------------

    if (dbug > 1) then
       if (my_task == master_task) then
          call mct_aVect_info(2, d2x, istr='run diag'//':AV', pe=localPet)
       end if
       call shr_nuopc_methods_State_diagnose(exportState,subname//':ES',rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    if (my_task == master_task) then
       call ESMF_ClockPrint(clock, options="currTime", preString="------>Advancing ATM from: ", rc=rc)
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
    character(len=256)       :: cvalue
    character(len=256)       :: restart_option       ! Restart option units
    integer                  :: restart_n            ! Number until restart interval
    integer                  :: restart_ymd          ! Restart date (YYYYMMDD)
    type(ESMF_ALARM)         :: restart_alarm
    integer                  :: first_time = .true.
    character(len=*),parameter :: subname=trim(modName)//':(ModelSetRunClock) '
    !-------------------------------------------------------------------------------

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
    ! set restart alarm
    !--------------------------------                                                                                 

    if (first_time) then
       call NUOPC_CompAttributeGet(gcomp, name="restart_option", value=restart_option, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp,  name="restart_n", value=cvalue, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_n

       call NUOPC_CompAttributeGet(gcomp, name="restart_ymd", value=cvalue, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) restart_ymd

       call shr_nuopc_time_alarmInit(mclock, restart_alarm, restart_option, &
            opt_n   = restart_n,           &
            opt_ymd = restart_ymd,         &
            RefTime = mcurrTime,           &
            alarmname = 'alarm_restart', rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call ESMF_AlarmSet(restart_alarm, clock=mclock, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       first_time = .false.
    end if

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
    implicit none
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

    call datm_comp_final(my_task, master_task, logunit)

    if (dbug > 5) call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine ModelFinalize

  !===============================================================================

end module atm_comp_nuopc
