module ocn_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for DOCN
  !----------------------------------------------------------------------------

  use ESMF
  use NUOPC            , only : NUOPC_CompDerive, NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize
  use NUOPC            , only : NUOPC_CompAttributeGet, NUOPC_Advertise
  use NUOPC_Model      , only : model_routine_SS        => SetServices
  use NUOPC_Model      , only : model_label_Advance     => label_Advance
  use NUOPC_Model      , only : model_label_SetRunClock => label_SetRunClock
  use NUOPC_Model      , only : model_label_Finalize    => label_Finalize
  use NUOPC_Model      , only : NUOPC_ModelGet
  use shr_kind_mod     , only : r8=>shr_kind_r8, i8=>shr_kind_i8, cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_sys_mod      , only : shr_sys_abort
  use shr_cal_mod      , only : shr_cal_ymd2date
  use shr_mpi_mod      , only : shr_mpi_bcast
  use dshr_methods_mod , only : dshr_state_diagnose, chkerr, memcheck
  use dshr_strdata_mod , only : shr_strdata_type, shr_strdata_advance, shr_strdata_init_from_xml
  use dshr_mod         , only : dshr_model_initphase, dshr_init, dshr_mesh_init
  use dshr_mod         , only : dshr_state_setscalar, dshr_set_runclock
  use dshr_dfield_mod  , only : dfield_type, dshr_dfield_add, dshr_dfield_copy
  use dshr_fldlist_mod , only : fldlist_type, dshr_fldlist_realize
  use perf_mod         , only : t_startf, t_stopf, t_adj_detailf, t_barrierf
  use pio

  ! Datamode specialized modules
  use docn_set_ofrac_mod           , only : docn_set_ofrac
  use docn_datamode_copyall_mod    , only : docn_datamode_copyall_advertise
  use docn_datamode_copyall_mod    , only : docn_datamode_copyall_init_pointers
  use docn_datamode_copyall_mod    , only : docn_datamode_copyall_advance
  use docn_datamode_copyall_mod    , only : docn_datamode_copyall_restart_read
  use docn_datamode_copyall_mod    , only : docn_datamode_copyall_restart_write
  use docn_datamode_iaf_mod        , only : docn_datamode_iaf_advertise
  use docn_datamode_iaf_mod        , only : docn_datamode_iaf_init_pointers
  use docn_datamode_iaf_mod        , only : docn_datamode_iaf_advance
  use docn_datamode_iaf_mod        , only : docn_datamode_iaf_restart_read
  use docn_datamode_iaf_mod        , only : docn_datamode_iaf_restart_write
  use docn_datamode_som_mod        , only : docn_datamode_som_advertise
  use docn_datamode_som_mod        , only : docn_datamode_som_init_pointers
  use docn_datamode_som_mod        , only : docn_datamode_som_advance
  use docn_datamode_som_mod        , only : docn_datamode_som_restart_read
  use docn_datamode_som_mod        , only : docn_datamode_som_restart_write
  use docn_datamode_aquaplanet_mod , only : docn_datamode_aquaplanet_advertise
  use docn_datamode_aquaplanet_mod , only : docn_datamode_aquaplanet_init_pointers
  use docn_datamode_aquaplanet_mod , only : docn_datamode_aquaplanet_advance

  implicit none
  private ! except

  public  :: SetServices

  private :: InitializeAdvertise
  private :: InitializeRealize
  private :: ModelAdvance
  private :: docn_comp_run
  private :: ModelFinalize

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  ! module variables common to all data models
  type(shr_strdata_type)       :: sdat
  type(ESMF_Mesh)              :: model_mesh
  character(CS)                :: flds_scalar_name = ''
  integer                      :: flds_scalar_num = 0
  integer                      :: flds_scalar_index_nx = 0
  integer                      :: flds_scalar_index_ny = 0
  integer                      :: compid           ! mct comp id
  integer                      :: mpicom           ! mpi communicator
  integer                      :: my_task          ! my task in mpi communicator mpicom
  logical                      :: masterproc       ! true of my_task == master_task
  character(len=16)            :: inst_suffix = "" ! char string associated with instance (ie. "_0001" or "")
  integer                      :: logunit          ! logging unit number
  logical                      :: restart_read     ! start from restart
  character(CL)                :: case_name
  character(*) , parameter     :: nullstr = 'undefined'

  ! docn_in namelist input
  character(CL)                :: xmlfilename = nullstr               ! filename to obtain namelist info from
  character(CL)                :: nlfilename = nullstr                ! filename to obtain namelist info from
  character(CL)                :: datamode = nullstr                  ! flags physics options wrt input data
  character(CL)                :: model_meshfile = nullstr            ! full pathname to model meshfile
  character(CL)                :: model_maskfile = nullstr            ! full pathname to obtain mask from
  character(CL)                :: model_createmesh_fromfile = nullstr ! full pathname to obtain mask from
  real(R8)                     :: sst_constant_value                  ! sst constant value
  integer                      :: aquap_option                        ! if aqua-planet mode, option to use
  character(CL)                :: restfilm = nullstr                  ! model restart file namelist
  integer                      :: nx_global
  integer                      :: ny_global

  ! linked lists
  type(fldList_type) , pointer :: fldsImport => null()
  type(fldList_type) , pointer :: fldsExport => null()
  type(dfield_type)  , pointer :: dfields    => null()

  ! constants
  logical                      :: aquaplanet = .false.
  integer      , parameter     :: master_task = 0                 ! task number of master task
  character(*) , parameter     :: module_name = "(ocn_comp_nuopc)"
  character(*) , parameter     :: modelname = 'docn'

  ! ocean fraction
  real(r8), allocatable        :: ocn_fraction(:)

  character(*) , parameter     :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Local varaibles
    character(len=*),parameter  :: subname=trim(module_name)//':(SetServices) '
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
    integer           :: inst_index         ! number of current instance (ie. 1)
    integer           :: nu                 ! unit number
    integer           :: ierr               ! error code
    logical           :: exists             ! check for file existence
    integer           :: n
    type(fldlist_type), pointer :: fldList
    character(len=*),parameter :: subname=trim(module_name)//':(InitializeAdvertise) '
    character(*)    ,parameter :: F00 = "('(ocn_comp_nuopc) ',8a)"
    character(*)    ,parameter :: F01 = "('(ocn_comp_nuopc) ',a,2x,i8)"
    character(*)    ,parameter :: F02 = "('(ocn_comp_nuopc) ',a,l6)"
    !-------------------------------------------------------------------------------

    namelist / docn_nml / datamode, model_meshfile, model_maskfile, model_createmesh_fromfile, &
         restfilm,  sst_constant_value, nx_global, ny_global

    rc = ESMF_SUCCESS

    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=case_name, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Obtain flds_scalar values, mpi values, multi-instance values and
    ! set logunit and set shr logging to my log file
    call dshr_init(gcomp, mpicom, my_task, inst_index, inst_suffix, &
         flds_scalar_name, flds_scalar_num, flds_scalar_index_nx, flds_scalar_index_ny, &
         logunit, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine logical masterproc
    masterproc = (my_task == master_task)

    ! Read docn_nml from nlfilename
    if (my_task == master_task) then
       nlfilename = "docn_in"//trim(inst_suffix)
       open (newunit=nu,file=trim(nlfilename),status="old",action="read")
       read (nu,nml=docn_nml,iostat=ierr)
       close(nu)
       if (ierr > 0) then
          write(logunit,F00) 'ERROR: reading input namelist, '//trim(nlfilename)//' iostat=',ierr
          call shr_sys_abort(subName//': namelist read error '//trim(nlfilename))
       end if

       ! write namelist input to standard out
       write(logunit,F00)' datamode = ',trim(datamode)
       if (model_createmesh_fromfile /= nullstr) then
          write(logunit,F00)' model_create_meshfile_fromfile = ',trim(model_createmesh_fromfile)
       else
          write(logunit,F00)' model_meshfile = ',trim(model_meshfile)
          write(logunit,F00)' model_maskfile = ',trim(model_maskfile)
       end if
       write(logunit,F01)' nx_global = ',nx_global
       write(logunit,F01)' ny_global = ',ny_global
       write(logunit,F00)' restfilm = ',trim(restfilm)

       ! check that files exists
       if (model_createmesh_fromfile /= nullstr) then
          inquire(file=trim(model_createmesh_fromfile), exist=exists)
          if (.not.exists) then
             write(logunit, *)' ERROR: model_createmesh_fromfile '//&
                  trim(model_createmesh_fromfile)//' does not exist'
             call shr_sys_abort(trim(subname)//' ERROR: model_createmesh_fromfile '//&
                  trim(model_createmesh_fromfile)//' does not exist')
          end if
       else
          inquire(file=trim(model_meshfile), exist=exists)
          if (.not.exists) then
             write(logunit, *)' ERROR: model_meshfile '//trim(model_meshfile)//' does not exist'
             call shr_sys_abort(trim(subname)//' ERROR: model_meshfile '//trim(model_meshfile)//' does not exist')
          end if
          inquire(file=trim(model_maskfile), exist=exists)
          if (.not.exists) then
             write(logunit, *)' ERROR: model_maskfile '//trim(model_maskfile)//' does not exist'
             call shr_sys_abort(trim(subname)//' ERROR: model_maskfile '//trim(model_maskfile)//' does not exist')
          end if
       end if
    endif

    ! Broadcast namelist input
    call shr_mpi_bcast(datamode                  , mpicom, 'datamode')
    call shr_mpi_bcast(model_meshfile            , mpicom, 'model_meshfile')
    call shr_mpi_bcast(model_maskfile            , mpicom, 'model_maskfile')
    call shr_mpi_bcast(model_createmesh_fromfile , mpicom, 'model_createmesh_fromfile')
    call shr_mpi_bcast(nx_global                 , mpicom, 'nx_global')
    call shr_mpi_bcast(ny_global                 , mpicom, 'ny_global')
    call shr_mpi_bcast(sst_constant_value        , mpicom, 'sst_constant_value')

    ! Special logic for prescribed aquaplanet
    if (datamode(1:9) == 'sst_aquap') then
       ! First determine the prescribed aquaplanet option
       if (len_trim(datamode) == 10) then
          read(datamode(10:10),'(i1)') aquap_option
       else if (len_trim(datamode) == 11) then
          read(datamode(10:11),'(i2)') aquap_option
       end if
       ! Now remove the index from the datamode value, to have a generic setting for later use
       datamode = "sst_aquap_analytic"
    end if

    ! Validate datamode
    if ( trim(datamode) == 'null'               .or. & ! does nothing
         trim(datamode) == 'sstdata'            .or. & ! read stream, no import data
         trim(datamode) == 'iaf'                .or. & ! read stream, needs import data?
         trim(datamode) == 'sst_aquap_file'     .or. & ! read stream, no import data
         trim(datamode) == 'som'                .or. & ! read stream, needs import data
         trim(datamode) == 'som_aquap'          .or. & ! read stream, needs import data
         trim(datamode) == 'sst_aquap_analytic' .or. & ! analytic, no streams, import or export data
         trim(datamode) == 'sst_aquap_constant' ) then ! analytic, no streams, import or export data
       ! success do nothing
    else
       call shr_sys_abort(' ERROR illegal docn datamode = '//trim(datamode))
    endif

    ! Advertise docn fields
    if (trim(datamode) /= 'NULL') then
       if (trim(datamode)=='sst_aquap_analytic' .or. trim(datamode)=='sst_aquap_constant') then
          aquaplanet = .true.
          call docn_datamode_aquaplanet_advertise(exportState, fldsExport, flds_scalar_name, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else if (trim(datamode(1:3)) == 'som') then
          call docn_datamode_som_advertise(importState, exportState, fldsImport, fldsExport, flds_scalar_name, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else if (trim(datamode) == 'sstdata' .or. trim(datamode) == 'sst_aquap_file') then
          call docn_datamode_copyall_advertise(exportState, fldsExport, flds_scalar_name, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else if (trim(datamode) == 'iaf') then
          call docn_datamode_iaf_advertise(importState, exportState, fldsImport, fldsExport, flds_scalar_name, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end if

  end subroutine InitializeAdvertise

  !===============================================================================
  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_TimeInterval)         :: TimeStep
    type(ESMF_Time)                 :: currTime
    integer                         :: current_ymd  ! model date
    integer                         :: current_year ! model year
    integer                         :: current_mon  ! model month
    integer                         :: current_day  ! model day
    integer                         :: current_tod  ! model sec into model date
    integer                         :: fieldcount
    type(ESMF_Field)                :: lfield
    character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
    integer                         :: n
    character(len=*), parameter :: subname=trim(module_name)//':(InitializeRealize) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Initialize model mesh, restart flag, compid, and logunit
    call t_startf('docn_strdata_init')
    call dshr_mesh_init(gcomp, compid, logunit, trim(modelname), nx_global, ny_global, &
         model_meshfile, model_maskfile, model_createmesh_fromfile, model_mesh, restart_read, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Initialize stream data type if not aqua planet
    if (.not. aquaplanet) then
       xmlfilename = trim(modelname)//'.streams'//trim(inst_suffix)//'.xml'
       call shr_strdata_init_from_xml(sdat, xmlfilename, model_mesh, clock, compid, logunit, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    call t_stopf('docn_strdata_init')

    ! Realize the actively coupled fields, now that a mesh is established and
    ! NUOPC_Realize "realizes" a previously advertised field in the importState and exportState
    ! by replacing the advertised fields with the newly created fields of the same name.
    call dshr_fldlist_realize( exportState, fldsExport, flds_scalar_name, flds_scalar_num, model_mesh, &
         subname//trim(modelname)//':Export', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_fldlist_realize( importState, fldsImport, flds_scalar_name, flds_scalar_num, model_mesh, &
         subname//trim(modelname)//':Import', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the time to interpolate the stream data to
    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(currTime, yy=current_year, mm=current_mon, dd=current_day, s=current_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(current_year, current_mon, current_day, current_ymd)

    ! Run docn
    call docn_comp_run(importState, exportState, clock, current_ymd, current_tod, restart_write=.false., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Add scalars to export state
    call dshr_state_SetScalar(dble(nx_global),flds_scalar_index_nx, exportState, flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_SetScalar(dble(ny_global),flds_scalar_index_ny, exportState, flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Diagnostics
    call dshr_state_diagnose(exportState,subname//':ES',rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

   end subroutine InitializeRealize

  !===============================================================================
  subroutine ModelAdvance(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_State)        :: importState, exportState
    type(ESMF_Clock)        :: clock
    type(ESMF_TimeInterval) :: timeStep
    type(ESMF_Time)         :: currTime, nextTime
    type(ESMF_Alarm)        :: alarm
    integer                 :: next_ymd      ! model date
    integer                 :: next_tod      ! model sec into model date
    integer                 :: yr            ! year
    integer                 :: mon           ! month
    integer                 :: day           ! day in month
    logical                 :: restart_write
    character(len=*),parameter :: subname=trim(module_name)//':(ModelAdvance) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call memcheck(subname, 5, my_task == master_task)

    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! For nuopc - the component clock is advanced at the end of the time interval
    ! Need to advance nuopc one timestep ahead for shr_strdata time interpolation
    call ESMF_ClockGet( clock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    nextTime = currTime + timeStep
    call ESMF_TimeGet( nextTime, yy=yr, mm=mon, dd=day, s=next_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yr, mon, day, next_ymd)

    ! determine if will write restart
    call ESMF_ClockGetAlarm(clock, alarmname='alarm_restart', alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       restart_write = .true.
    else
       restart_write = .false.
    end if

    ! run docn
    call docn_comp_run(importState, exportState, clock, next_ymd, next_tod, restart_write, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine ModelAdvance

  !===============================================================================
  subroutine docn_comp_run(importState, exportState, clock, target_ymd, target_tod, restart_write, rc)

    ! --------------------------
    ! advance docn
    ! --------------------------

    ! input/output variables:
    type(ESMF_Clock) , intent(in)    :: clock
    type(ESMF_State) , intent(inout) :: importState
    type(ESMF_State) , intent(inout) :: exportState
    integer          , intent(in)    :: target_ymd       ! model date
    integer          , intent(in)    :: target_tod       ! model sec into model date
    logical          , intent(in)    :: restart_write
    integer          , intent(out)   :: rc

    ! local variables
    logical :: first_time = .true.
    integer :: numOwnedElements
    character(*), parameter :: subName = "(docn_comp_run) "
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('DOCN_RUN')

    !--------------------
    ! First time initialization
    !--------------------

    if (first_time) then

       ! Determine ocn fraction
       call ESMF_MeshGet(model_mesh, numOwnedElements=numOwnedElements,  rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       allocate(ocn_fraction(numOwnedElements))
       call docn_set_ofrac(exportState, model_mesh, model_meshfile, model_maskfile, &
            nx_global, ny_global, compid, ocn_fraction, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    
       ! Initialize dfields
       call docn_init_dfields(importState, exportState, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Initialize datamode module ponters
       select case (trim(datamode))
       case('sstdata', 'sst_aquap_file')
          call docn_datamode_copyall_init_pointers(exportState, ocn_fraction, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       case('iaf')
          call docn_datamode_iaf_init_pointers(importState, exportState, ocn_fraction, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       case('som', 'som_aquap')
          call docn_datamode_som_init_pointers(importState, exportState, sdat, ocn_fraction, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       case('sst_aquap_analytic', 'sst_aquap_constant')
          call  docn_datamode_aquaplanet_init_pointers(exportState, ocn_fraction, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end select

       ! Read restart if needed
       if (restart_read) then
          select case (trim(datamode))
          case('sstdata', 'sst_aquap_file')
             call docn_datamode_copyall_restart_read(restfilm, inst_suffix, logunit, my_task, mpicom, sdat)
          case('iaf')
             call docn_datamode_iaf_restart_read(restfilm, inst_suffix, logunit, my_task, mpicom, sdat)
          case('som', 'som_aquap')
             call docn_datamode_som_restart_read(restfilm, inst_suffix, logunit, my_task, mpicom, sdat)
          end select
       end if

       ! Reset first_time
       first_time = .false.
    end if

    !--------------------
    ! Update export (and possibly import data model states)
    !--------------------

    ! Advance data model streams - time and spatially interpolate to model time and grid
    call t_barrierf('docn_BARRIER',mpicom)
    call t_startf('docn_strdata_advance')
    call shr_strdata_advance(sdat, target_ymd, target_tod, logunit, 'docn', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('docn_strdata_advance')

    ! Copy all fields from streams to export state as default
    ! This automatically will update the fields in the export state
    call t_barrierf('docn_dfield_copy_BARRIER', mpicom)
    call t_startf('docn_dfield_copy')
    if(.not. aquaplanet) then
       call dshr_dfield_copy(dfields, sdat, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif
    call t_stopf('docn_dfield_copy')

    ! Perform data mode specific calculations
    select case (trim(datamode))
    case('sstdata','sst_aquap_file')
       call  docn_datamode_copyall_advance(rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    case('iaf')
       call  docn_datamode_iaf_advance(rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    case('som','som_aquap')
       call docn_datamode_som_advance(importState, exportState, clock, restart_read, datamode, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    case('sst_aquap_analytic')
       call  docn_datamode_aquaplanet_advance(exportstate, model_mesh, sst_option=aquap_option, rc=rc)
       if (chkerr(rc,__line__,u_file_u)) return
    case('sst_aquap_constant')
       call  docn_datamode_aquaplanet_advance(exportState, model_mesh, sst_constant_value=sst_constant_value, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end select

    ! Write restarts if needed (no restarts for aquaplanet analytic or aquaplanet input file)
    if (restart_write) then
       select case (trim(datamode))
       case('sstdata','sst_aquap_file')
          call docn_datamode_copyall_restart_write(case_name, inst_suffix, target_ymd, target_tod, &
               logunit, mpicom, my_task, sdat)
       case('iaf')
          call docn_datamode_iaf_restart_write(case_name, inst_suffix, target_ymd, target_tod, &
               logunit, mpicom, my_task, sdat)
       case('som','som_aquap')
          call docn_datamode_som_restart_write(case_name, inst_suffix, target_ymd, target_tod, &
               logunit, mpicom, my_task, sdat)
       end select
    end if

    call t_stopf('DOCN_RUN')

    ! write diagnostics
    call dshr_state_diagnose(exportState,subname//':ES',rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  contains

    subroutine docn_init_dfields(importState, exportState, rc)
      ! -----------------------------
      ! Initialize dfields arrays
      ! -----------------------------

      ! input/output variables
      type(ESMF_State)       , intent(inout) :: importState
      type(ESMF_State)       , intent(inout) :: exportState
      integer                , intent(out)   :: rc

      ! local variables
      integer                         :: n
      integer                         :: fieldcount
      type(ESMF_Field)                :: lfield
      character(ESMF_MAXSTR) ,pointer :: lfieldnamelist(:)
      character(*), parameter   :: subName = "(docn_init_dfields) "
      !-------------------------------------------------------------------------------

      rc = ESMF_SUCCESS

      ! Initialize dfields data type (to map streams to export state fields)
      ! Create dfields linked list - used for copying stream fields to export state fields
      call ESMF_StateGet(exportState, itemCount=fieldCount, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      allocate(lfieldnamelist(fieldCount))
      call ESMF_StateGet(exportState, itemNameList=lfieldnamelist, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      do n = 1, fieldCount
         call ESMF_StateGet(exportState, itemName=trim(lfieldNameList(n)), field=lfield, rc=rc)
         if (chkerr(rc,__LINE__,u_FILE_u)) return
         if (trim(lfieldnamelist(n)) /= flds_scalar_name) then
            call dshr_dfield_add( dfields, sdat, trim(lfieldnamelist(n)), trim(lfieldnamelist(n)), exportState, &
                 logunit, masterproc, rc)
            if (chkerr(rc,__LINE__,u_FILE_u)) return
         end if
      end do
    end subroutine docn_init_dfields

  end subroutine docn_comp_run

  !===============================================================================
  subroutine ModelFinalize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    !-------------------------------------------------------------------------------
    rc = ESMF_SUCCESS
    if (my_task == master_task) then
       write(logunit,*)
       write(logunit,*) 'docn : end of main integration loop'
       write(logunit,*)
    end if
  end subroutine ModelFinalize

end module ocn_comp_nuopc
