module atm_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for DATM
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
  use shr_const_mod    , only : shr_const_cday
  use shr_sys_mod      , only : shr_sys_abort
  use shr_cal_mod      , only : shr_cal_ymd2date, shr_cal_ymd2julian, shr_cal_date2julian
  use shr_mpi_mod      , only : shr_mpi_bcast
  use dshr_methods_mod , only : dshr_state_diagnose, chkerr, memcheck
  use dshr_strdata_mod , only : shr_strdata_type, shr_strdata_init_from_xml, shr_strdata_advance
  use dshr_strdata_mod , only : shr_strdata_get_stream_pointer, shr_strdata_setOrbs
  use dshr_mod         , only : dshr_model_initphase, dshr_init
  use dshr_mod         , only : dshr_state_setscalar, dshr_set_runclock, dshr_log_clock_advance
  use dshr_mod         , only : dshr_restart_read, dshr_restart_write, dshr_mesh_init
  use dshr_mod         , only : dshr_orbital_init, dshr_orbital_update
  use dshr_dfield_mod  , only : dfield_type, dshr_dfield_add, dshr_dfield_copy
  use dshr_fldlist_mod , only : fldlist_type, dshr_fldlist_add, dshr_fldlist_realize
  use perf_mod         , only : t_startf, t_stopf, t_barrierf

  use datm_datamode_core2_mod   , only : datm_datamode_core2_advertise
  use datm_datamode_core2_mod   , only : datm_datamode_core2_init_pointers
  use datm_datamode_core2_mod   , only : datm_datamode_core2_advance
  use datm_datamode_core2_mod   , only : datm_datamode_core2_restart_write
  use datm_datamode_core2_mod   , only : datm_datamode_core2_restart_read
  use datm_datamode_jra_mod     , only : datm_datamode_jra_advertise
  use datm_datamode_jra_mod     , only : datm_datamode_jra_init_pointers
  use datm_datamode_jra_mod     , only : datm_datamode_jra_advance
  use datm_datamode_jra_mod     , only : datm_datamode_jra_restart_write
  use datm_datamode_jra_mod     , only : datm_datamode_jra_restart_read
  use datm_datamode_clmncep_mod , only : datm_datamode_clmncep_advertise
  use datm_datamode_clmncep_mod , only : datm_datamode_clmncep_init_pointers
  use datm_datamode_clmncep_mod , only : datm_datamode_clmncep_advance
  use datm_datamode_clmncep_mod , only : datm_datamode_clmncep_restart_write
  use datm_datamode_clmncep_mod , only : datm_datamode_clmncep_restart_read

  implicit none
  private ! except

  public  :: SetServices

  private :: InitializeAdvertise
  private :: InitializeRealize
  private :: ModelAdvance
  private :: datm_comp_run
  private :: ModelFinalize

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  type(shr_strdata_type)       :: sdat
  type(ESMF_Mesh)              :: model_mesh                ! model mesh
  character(len=128)           :: flds_scalar_name = ''
  integer                      :: flds_scalar_num = 0
  integer                      :: flds_scalar_index_nx = 0
  integer                      :: flds_scalar_index_ny = 0
  integer                      :: flds_scalar_index_nextsw_cday = 0
  integer                      :: compid                    ! component id (needed by pio)
  integer                      :: mpicom                    ! mpi communicator
  integer                      :: my_task                   ! my task in mpi communicator mpicom
  logical                      :: masterproc                ! true of my_task == master_task
  integer                      :: inst_index                ! number of current instance (ie. 1)
  character(len=16)            :: inst_suffix = ""          ! char string associated with instance (ie. "_0001" or "")
  integer                      :: logunit                   ! logging unit number
  logical                      :: restart_read              ! start from restart
  character(CL)                :: case_name     ! case name
  character(len=*) , parameter :: nullstr = 'undefined'

  ! datm_in namelist input
  character(CL)                :: nlfilename = nullstr                ! filename to obtain namelist info from
  character(CL)                :: xmlfilename = nullstr                ! filename to obtain namelist info from
  character(CL)                :: dataMode = nullstr                  ! flags physics options wrt input data
  character(CL)                :: model_meshfile = nullstr            ! full pathname to model meshfile
  character(CL)                :: model_maskfile = nullstr            ! full pathname to obtain mask from
  character(CL)                :: model_createmesh_fromfile = nullstr ! full pathname to obtain mask from
  integer                      :: iradsw = 0                          ! radiation interval (input namelist)
  character(CL)                :: factorFn_mesh = 'null'              ! file containing correction factors mesh
  character(CL)                :: factorFn_data = 'null'              ! file containing correction factors data
  logical                      :: presaero = .false.                  ! true => send valid prescribe aero fields to coupler
  character(CL)                :: bias_correct = nullstr              ! send bias correction fields to coupler (not used here)
  character(CL)                :: anomaly_forcing(8) = nullstr        ! send anomaly forcing fields to coupler (not used here)
  logical                      :: force_prognostic_true = .false.     ! if true set prognostic true
  logical                      :: wiso_datm = .false.                 ! expect isotopic forcing from file?
  character(CL)                :: restfilm = nullstr                  ! model restart file namelist
  integer                      :: nx_global
  integer                      :: ny_global
                                                                      ! config attribute intput
  ! linked lists
  type(fldList_type) , pointer :: fldsImport => null()
  type(fldList_type) , pointer :: fldsExport => null()
  type(dfield_type)  , pointer :: dfields    => null()

  ! constants
  logical                      :: flds_co2
  logical                      :: flds_wiso
  integer                      :: idt                                 ! integer model timestep
  logical                      :: diagnose_data = .true.
  integer          , parameter :: master_task   = 0                   ! task number of master task
  character(len=*) , parameter :: rpfile        = 'rpointer.atm'
  character(*)     , parameter :: modName       = "(atm_comp_nuopc)"

  character(*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(len=*),parameter  :: subname=trim(modName)//':(SetServices) '
    !-------------------------------------------------------------------------------

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
    character(len=CL) :: cvalue     ! temporary
    integer           :: nu         ! unit number
    integer           :: ierr       ! error code
    logical           :: exists     ! check for file existence
    logical           :: flds_co2a, flds_co2b, flds_co2c 
    character(len=*),parameter :: subname='(atm_comp_nuopc):(InitializeAdvertise) '
    character(*)    ,parameter :: F00 = "('(atm_comp_nuopc) ',8a)"
    character(*)    ,parameter :: F01 = "('(atm_comp_nuopc) ',a,2x,i8)"
    character(*)    ,parameter :: F02 = "('(atm_comp_nuopc) ',a,l6)"
    !-------------------------------------------------------------------------------

    namelist / datm_nml / datamode, model_meshfile, model_maskfile, model_createmesh_fromfile, &
         nx_global, ny_global, restfilm, &
         iradsw, factorFn_data, factorFn_mesh, presaero, bias_correct, &
         anomaly_forcing, force_prognostic_true, wiso_datm

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

    ! Read atm_nml from nlfilename
    if (my_task == master_task) then
       nlfilename = "datm_in"//trim(inst_suffix)
       open (newunit=nu,file=trim(nlfilename),status="old",action="read")
       read (nu,nml=datm_nml,iostat=ierr)
       close(nu)
       if (ierr > 0) then
          write(logunit,*) 'ERROR: reading input namelist, '//trim(nlfilename)//' iostat=',ierr
          call shr_sys_abort(subName//': namelist read error '//trim(nlfilename))
       end if
    end if

    call shr_mpi_bcast(datamode                  , mpicom, 'datamode')
    call shr_mpi_bcast(model_meshfile            , mpicom, 'model_meshfile')
    call shr_mpi_bcast(model_maskfile            , mpicom, 'model_maskfile')
    call shr_mpi_bcast(model_createmesh_fromfile , mpicom, 'model_createmesh_fromfile')
    call shr_mpi_bcast(nx_global                 , mpicom, 'nx_global')
    call shr_mpi_bcast(ny_global                 , mpicom, 'ny_global')
    call shr_mpi_bcast(iradsw                    , mpicom, 'iradsw')
    call shr_mpi_bcast(factorFn_data             , mpicom, 'factorFn_data')
    call shr_mpi_bcast(factorFn_mesh             , mpicom, 'factorFn_mesh')
    call shr_mpi_bcast(restfilm                  , mpicom, 'restfilm')
    call shr_mpi_bcast(presaero                  , mpicom, 'presaero')
    call shr_mpi_bcast(wiso_datm                 , mpicom, 'wiso_datm')
    call shr_mpi_bcast(force_prognostic_true     , mpicom, 'force_prognostic_true')

    ! write namelist input to standard out
    if (my_task == master_task) then
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
       write(logunit,F02)' force_prognostic_true = ',force_prognostic_true
       write(logunit,F01)' iradsw = ',iradsw
       write(logunit,F00)' factorFn_data = ',trim(factorFn_data)
       write(logunit,F00)' factorFn_mesh = ',trim(factorFn_mesh)
       write(logunit,F00)' restfilm = ',trim(restfilm)
       write(logunit,F02)' presaero  = ',presaero
       write(logunit,F02)' wiso_datm = ',wiso_datm
    end if

    ! check that files exists
    if (my_task == master_task) then
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

    ! Validate sdat datamode
    if (masterproc) write(logunit,*) ' datm datamode = ',trim(datamode)
    if ( trim(datamode) == 'NULL'         .or. &
         trim(datamode) == 'CORE2_NYF'    .or. &
         trim(datamode) == 'CORE2_IAF'    .or. &
         trim(datamode) == 'CORE_IAF_JRA' .or. &
         trim(datamode) == 'CLMNCEP') then
    else
       call shr_sys_abort(' ERROR illegal datm datamode = '//trim(datamode))
    endif

    if (datamode /= 'NULL') then
       call NUOPC_CompAttributeGet(gcomp, name='flds_co2a', value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) flds_co2a
       call ESMF_LogWrite('flds_co2a = '// trim(cvalue), ESMF_LOGMSG_INFO)

       call NUOPC_CompAttributeGet(gcomp, name='flds_co2b', value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) flds_co2b
       call ESMF_LogWrite('flds_co2b = '// trim(cvalue), ESMF_LOGMSG_INFO)

       call NUOPC_CompAttributeGet(gcomp, name='flds_co2c', value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) flds_co2c
       call ESMF_LogWrite('flds_co2c = '// trim(cvalue), ESMF_LOGMSG_INFO)

       if (flds_co2a .or. flds_co2b .or. flds_co2c) then
          flds_co2 = .true.
       else
          flds_co2 = .false.
       end if
       
       call NUOPC_CompAttributeGet(gcomp, name='flds_wiso', value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) flds_wiso
       call ESMF_LogWrite('flds_wiso = '// trim(cvalue), ESMF_LOGMSG_INFO)

       ! check that flds_wiso matches wiso_datm
       ! TODO: should not need wiso_datm for nuopc data models
       if (wiso_datm /= flds_wiso) then
          call shr_sys_abort(subName//': datm namelist wiso_datm must match nuopc attribute flds_wiso')
       end if
    end if

    ! Advertise datm fields
    select case (trim(datamode))
    case ('CORE2_NYF', 'CORE2_IAF')
       call datm_datamode_core2_advertise(exportState, fldsExport, flds_scalar_name, &
            flds_co2, flds_wiso, presaero, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    case ('CORE2_IAF_JRA')
       call datm_datamode_jra_advertise(exportState, fldsExport, flds_scalar_name, rc)       
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    case ('CLMNCEP')
       call datm_datamode_clmncep_advertise(exportState, fldsExport, flds_scalar_name, &
            flds_co2, flds_wiso, presaero, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end select

  end subroutine InitializeAdvertise

  !===============================================================================

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_TimeInterval) :: timeStep
    type(ESMF_TIME)         :: currTime
    integer                 :: current_ymd   ! model date
    integer                 :: current_year  ! model year
    integer                 :: current_mon   ! model month
    integer                 :: current_day   ! model day
    integer                 :: current_tod   ! model sec into model date
    integer(i8)             :: stepno        ! step number
    real(r8)                :: nextsw_cday   ! calendar of next atm sw
    character(CL)           :: cvalue        ! character string for input config
    real(R8)                :: orbEccen      ! orb eccentricity (unit-less)
    real(R8)                :: orbMvelpp     ! orb moving vernal eq (radians)
    real(R8)                :: orbLambm0     ! orb mean long of perhelion (radians)
    real(R8)                :: orbObliqr     ! orb obliquity (radians)
    logical                 :: isPresent, isSet
    character(len=*), parameter :: subname=trim(modName)//':(InitializeRealize) '
    !-------------------------------------------------------------------------------

    if (datamode == 'NULL') RETURN

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! Initialize mesh, restart flag, compid, and logunit
    call t_startf('datm_strdata_init')

    call dshr_mesh_init(gcomp, compid, logunit, 'atm', nx_global, ny_global, &
         model_meshfile, model_maskfile, model_createmesh_fromfile, model_mesh, restart_read, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Initialize stream data type
    xmlfilename = 'datm.streams'//trim(inst_suffix)//'.xml'
    call shr_strdata_init_from_xml(sdat, xmlfilename, model_mesh, clock, compid, logunit, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('datm_strdata_init')

    ! NUOPC_Realize "realizes" a previously advertised field in the importState and exportState
    ! by replacing the advertised fields with the newly created fields of the same name.
    call dshr_fldlist_realize( exportState, fldsExport, flds_scalar_name, flds_scalar_num, model_mesh, &
         subname//':datmExport', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_fldlist_realize( importState, fldsImport, flds_scalar_name, flds_scalar_num, model_mesh, &
         subname//':datmImport', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the time to interpolate the stream data to
    call ESMF_ClockGet( clock, currTime=currTime, timeStep=timeStep, advanceCount=stepno, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(currTime, yy=current_year, mm=current_mon, dd=current_day, s=current_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(current_year, current_mon, current_day, current_ymd)

    ! Get model timestep (idt is module variable)
    call ESMF_TimeIntervalGet( timeStep, s=idt, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Initialize and update orbital values
    call dshr_orbital_init(gcomp, logunit, my_task == master_task, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_orbital_update(clock, logunit, my_task == master_task, &
         orbEccen, orbObliqr, orbLambm0, orbMvelpp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Run datm
    call datm_comp_run(importstate, exportstate, current_ymd, current_tod, current_mon, &
         orbEccen, orbMvelpp, orbLambm0, orbObliqr, restart_write=.false., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Add scalars to export state
    call dshr_state_SetScalar(dble(nx_global),flds_scalar_index_nx, exportState, flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_SetScalar(dble(ny_global),flds_scalar_index_ny, exportState, flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxNextSwCday", value=cvalue, &
         isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read (cvalue,*) flds_scalar_index_nextsw_cday
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldIdxNextSwCday')
    endif

    nextsw_cday = getNextRadCDay( current_ymd, current_tod, stepno, idt, iradsw, sdat%model_calendar )
    call dshr_state_SetScalar(nextsw_cday, flds_scalar_index_nextsw_cday, exportState, flds_scalar_name, flds_scalar_num, rc)
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
    type(ESMF_Time)         :: time
    type(ESMF_Alarm)        :: alarm
    type(ESMF_Time)         :: currTime
    type(ESMF_Time)         :: nextTime
    type(ESMF_TimeInterval) :: timeStep
    real(r8)                :: nextsw_cday
    logical                 :: restart_write ! restart alarm is ringing
    integer                 :: next_ymd      ! model date
    integer                 :: next_tod      ! model sec into model date
    integer                 :: yr, mon, day  ! year, month, day
    integer(i8)             :: stepno        ! step number
    real(R8)                :: orbEccen      ! orb eccentricity (unit-less)
    real(R8)                :: orbMvelpp     ! orb moving vernal eq (radians)
    real(R8)                :: orbLambm0     ! orb mean long of perhelion (radians)
    real(R8)                :: orbObliqr     ! orb obliquity (radians)
    character(len=CL)       :: cvalue        ! temporary
    character(len=*),parameter  :: subname=trim(modName)//':(ModelAdvance) '
    !-------------------------------------------------------------------------------

    if (datamode == 'NULL') RETURN

    rc = ESMF_SUCCESS

    call t_startf(subname)
    call memcheck(subname, 5, my_task==master_task)

    ! Query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! For nuopc - the component clock is advanced at the end of the time interval
    ! For these to match for now - need to advance nuopc one timestep ahead for
    ! shr_strdata time interpolation
    call ESMF_ClockGet( clock, currTime=currTime, timeStep=timeStep, advanceCount=stepno, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    nextTime = currTime + timeStep
    call ESMF_TimeGet( nextTime, yy=yr, mm=mon, dd=day, s=next_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yr, mon, day, next_ymd)

    ! Update the orbital values
    call dshr_orbital_update(clock, logunit, my_task == master_task, &
         orbEccen, orbObliqr, orbLambm0, orbMvelpp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

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
    endif

    ! Run datm
    call t_startf('datm_run')
    call datm_comp_run(importstate, exportstate, next_ymd, next_tod, mon, &
         orbEccen, orbMvelpp, orbLambm0, orbObliqr, restart_write,  rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('datm_run')

    ! Update nextsw_cday for scalar data
    ! Use nextYMD and nextTOD here since since the component - clock is advance at the END of the time interval
    nextsw_cday = getNextRadCDay( next_ymd, next_tod, stepno, idt, iradsw, sdat%model_calendar )
    call dshr_state_SetScalar(nextsw_cday, flds_scalar_index_nextsw_cday, exportState, flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call t_stopf(subname)

  end subroutine ModelAdvance

  !===============================================================================
  subroutine datm_comp_run(importState, exportState, target_ymd, target_tod, target_mon, &
       orbEccen, orbMvelpp, orbLambm0, orbObliqr, restart_write, rc)

    ! ----------------------------------
    ! run method for datm model
    ! ----------------------------------

    ! input/output variables
    type(ESMF_State)       , intent(inout) :: importState
    type(ESMF_State)       , intent(inout) :: exportState
    integer                , intent(in)    :: target_ymd       ! model date
    integer                , intent(in)    :: target_tod       ! model sec into model date
    integer                , intent(in)    :: target_mon       ! model month
    real(R8)               , intent(in)    :: orbEccen         ! orb eccentricity (unit-less)
    real(R8)               , intent(in)    :: orbMvelpp        ! orb moving vernal eq (radians)
    real(R8)               , intent(in)    :: orbLambm0        ! orb mean long of perhelion (radians)
    real(R8)               , intent(in)    :: orbObliqr        ! orb obliquity (radians)
    logical                , intent(in)    :: restart_write
    integer                , intent(out)   :: rc

    ! local variables
    logical :: first_time = .true.
    character(*), parameter :: subName = '(datm_comp_run) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('DATM_RUN')

    !--------------------
    ! First time initialization
    !--------------------

    if (first_time) then

       ! Initialize dfields
       call datm_init_dfields(rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Initialize datamode module ponters
       select case (trim(datamode))
       case('CORE2_NYF','CORE2_IAF')
          call datm_datamode_core2_init_pointers(exportState, sdat, datamode, factorfn_mesh, factorfn_data, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       case('CORE_IAF_JRA')
          call datm_datamode_jra_init_pointers(exportState, sdat, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       case('CLMNCEP')
          call datm_datamode_clmncep_init_pointers(importState, exportState, sdat, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end select

       ! Read restart if needed
       if (restart_read) then
          select case (trim(datamode))
          case('CORE2_NYF','CORE2_IAF')
             call datm_datamode_core2_restart_read(restfilm, inst_suffix, logunit, my_task, mpicom, sdat)
          case('CORE_IAF_JRA')
             call datm_datamode_jra_restart_read(restfilm, inst_suffix, logunit, my_task, mpicom, sdat)
          case('CLMNCEP')
             call datm_datamode_clmncep_restart_read(restfilm, inst_suffix, logunit, my_task, mpicom, sdat)
          end select
       end if

       ! reset first_time
       first_time = .false.
    end if

    !--------------------
    ! Advance datm streams
    !--------------------

    ! set data needed for cosz t-interp method
    call shr_strdata_setOrbs(sdat, orbEccen, orbMvelpp, orbLambm0, orbObliqr, idt)

    ! time and spatially interpolate to model time and grid
    call t_barrierf('datm_BARRIER',mpicom)
    call t_startf('datm_strdata_advance')
    call shr_strdata_advance(sdat, target_ymd, target_tod, logunit, 'datm', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('datm_strdata_advance')

    ! copy all fields from streams to export state as default
    ! This automatically will update the fields in the export state
    call t_barrierf('datm_comp_dfield_copy_BARRIER', mpicom)
    call t_startf('datm_dfield_copy')
    call dshr_dfield_copy(dfields, sdat, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('datm_dfield_copy')

    ! Determine data model behavior based on the mode
    call t_startf('datm_datamode')
    select case (trim(datamode))
    case('CORE2_NYF','CORE2_IAF')
       call datm_datamode_core2_advance(exportstate, datamode, target_ymd, target_tod, target_mon, &
            sdat%model_calendar, factorfn_mesh, factorfn_data, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    case('CORE_IAF_JRA')
       call datm_datamode_jra_advance(exportstate, target_ymd, target_tod, sdat%model_calendar, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    case('CLMNCEP')
       call datm_datamode_clmncep_advance(importstate, exportstate, masterproc, logunit, mpicom,  rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end select

    ! Write restarts if needed
    if (restart_write) then
       select case (trim(datamode))
       case('CORE2_NYF','CORE2_IAF')
          call datm_datamode_core2_restart_write(case_name, inst_suffix, target_ymd, target_tod, &
               logunit, mpicom, my_task, sdat)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       case('CORE_IAF_JRA')
          call datm_datamode_jra_restart_write(case_name, inst_suffix, target_ymd, target_tod, &
               logunit, mpicom, my_task, sdat)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       case('CLMNCEP')
          call datm_datamode_clmncep_restart_write(case_name, inst_suffix, target_ymd, target_tod, &
               logunit, mpicom, my_task, sdat)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end select
    end if

    ! Log output for model date
    if (masterproc) write(logunit,*) 'atm : model date ', target_ymd, target_tod

    ! Diagnostics
    if (diagnose_data) then
       call dshr_state_diagnose(exportState, subname//':ES',rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    call t_stopf('datm_datamode')
    call t_stopf('DATM_RUN')

  !--------
  contains
  !--------

    subroutine datm_init_dfields(rc)
      ! -----------------------------
      ! Initialize dfields arrays 
      ! (for export fields that have a corresponding stream field)
      ! -----------------------------
      
      ! input/output parameters
      integer, intent(out)   :: rc

      ! local variables
      integer                         :: n
      character(CS)                   :: strm_flds3(3)
      character(CS)                   :: strm_flds4(4)
      integer                         :: rank
      integer                         :: fieldcount
      type(ESMF_Field)                :: lfield
      character(ESMF_MAXSTR) ,pointer :: lfieldnames(:)
      character(*), parameter   :: subName = "(datm_init_dfields) "
      !-------------------------------------------------------------------------------

      rc = ESMF_SUCCESS

      call ESMF_StateGet(exportState, itemCount=fieldCount, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      allocate(lfieldnames(fieldCount))
      call ESMF_StateGet(exportState, itemNameList=lfieldnames, rc=rc)
      if (chkerr(rc,__LINE__,u_FILE_u)) return
      do n = 1, fieldCount
         call ESMF_StateGet(exportState, itemName=trim(lfieldnames(n)), field=lfield, rc=rc)
         if (chkerr(rc,__LINE__,u_FILE_u)) return
         call ESMF_FieldGet(lfield, rank=rank, rc=rc)
         if (chkerr(rc,__LINE__,u_FILE_u)) return
         if (rank == 1) then
            call dshr_dfield_add( dfields, sdat, trim(lfieldnames(n)), trim(lfieldnames(n)), &
                 exportState, logunit, masterproc, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
         else if (rank == 2) then
            select case (trim(lfieldnames(n)))
            case('Faxa_bcph')
               strm_flds3 = (/'Faxa_bcphidry', 'Faxa_bcphodry', 'Faxa_bcphiwet'/)
               call dshr_dfield_add(dfields, sdat, trim(lfieldnames(n)), strm_flds3, exportState, logunit, masterproc, rc)
               if (ChkErr(rc,__LINE__,u_FILE_u)) return
            case('Faxa_ocph')
               strm_flds3 = (/'Faxa_ocphidry', 'Faxa_ocphodry', 'Faxa_ocphiwet'/)
               call dshr_dfield_add(dfields, sdat, trim(lfieldnames(n)), strm_flds3, exportState, logunit, masterproc, rc)
               if (ChkErr(rc,__LINE__,u_FILE_u)) return
            case('Faxa_dstwet')
               strm_flds4 = (/'Faxa_dstwet1', 'Faxa_dstwet2', 'Faxa_dstwet3', 'Faxa_dstwet4'/)
               call dshr_dfield_add(dfields, sdat, trim(lfieldnames(n)), strm_flds4, exportState, logunit, masterproc, rc)
               if (ChkErr(rc,__LINE__,u_FILE_u)) return
            case('Faxa_dstdry')
               strm_flds4 = (/'Faxa_dstdry1', 'Faxa_dstdry2', 'Faxa_dstdry3', 'Faxa_dstdry4'/)
               call dshr_dfield_add(dfields, sdat, trim(lfieldnames(n)), strm_flds4, exportState, logunit, masterproc, rc)
               if (ChkErr(rc,__LINE__,u_FILE_u)) return
            case('Faxa_rainc_wiso')
               strm_flds3 = (/'Faxa_rainc_16O', 'Faxa_rainc_18O', 'Faxa_rainc_HDO'/)
               call dshr_dfield_add(dfields, sdat, trim(lfieldnames(n)), strm_flds3, exportState, logunit, masterproc, rc)
               if (ChkErr(rc,__LINE__,u_FILE_u)) return
            case('Faxa_rainl_wiso')
               strm_flds3 = (/'Faxa_rainl_16O', 'Faxa_rainl_18O', 'Faxa_rainl_HDO'/)
               call dshr_dfield_add(dfields, sdat, trim(lfieldnames(n)), strm_flds3, exportState, logunit, masterproc, rc)
               if (ChkErr(rc,__LINE__,u_FILE_u)) return
            case('Faxa_snowc_wiso')
               strm_flds3 = (/'Faxa_snowc_16O', 'Faxa_snowc_18O', 'Faxa_snowc_HDO'/)
               call dshr_dfield_add(dfields, sdat, trim(lfieldnames(n)), strm_flds3, exportState, logunit, masterproc, rc)
               if (ChkErr(rc,__LINE__,u_FILE_u)) return
            case('Faxa_snowl_wiso')
               strm_flds3 = (/'Faxa_snowl_16O', 'Faxa_snowl_18O', 'Faxa_snowl_HDO'/)
               call dshr_dfield_add(dfields, sdat, trim(lfieldnames(n)), strm_flds3, exportState, logunit, masterproc, rc)
               if (ChkErr(rc,__LINE__,u_FILE_u)) return
            end select
         end if
      end do
    end subroutine datm_init_dfields

  end subroutine datm_comp_run

  !===============================================================================
  real(R8) function getNextRadCDay( ymd, tod, stepno, dtime, iradsw, calendar )

    !  Return the calendar day of the next radiation time-step.
    !  General Usage: nextswday = getNextRadCDay(curr_date)

    ! input/output variables
    integer    , intent(in)    :: ymd
    integer    , intent(in)    :: tod
    integer(i8), intent(in)    :: stepno
    integer    , intent(in)    :: dtime
    integer    , intent(in)    :: iradsw
    character(*),intent(in)    :: calendar

    ! local variables
    real(R8) :: nextsw_cday
    real(R8) :: julday
    integer  :: liradsw
    integer  :: yy,mm,dd
    character(*),parameter :: subName =  '(getNextRadCDay) '
    !-------------------------------------------------------------------------------

    liradsw = iradsw
    if (liradsw < 0) liradsw  = nint((-liradsw *3600._r8)/dtime)
    call shr_cal_date2julian(ymd,tod,julday,calendar)
    if (liradsw > 1) then
       if (mod(stepno+1,liradsw) == 0 .and. stepno > 0) then
          nextsw_cday = julday + 2*dtime/shr_const_cday
       else
          nextsw_cday = -1._r8
       end if
    else
       nextsw_cday = julday + dtime/shr_const_cday
    end if
    getNextRadCDay = nextsw_cday

  end function getNextRadCDay

  !===============================================================================
  subroutine ModelFinalize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (my_task == master_task) then
       write(logunit,*)
       write(logunit,*) 'datm : end of main integration loop'
       write(logunit,*)
    end if
  end subroutine ModelFinalize

end module atm_comp_nuopc
