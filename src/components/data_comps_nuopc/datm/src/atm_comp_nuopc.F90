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
  use shr_file_mod     , only : shr_file_getlogunit, shr_file_setlogunit
  use shr_kind_mod     , only : r8=>shr_kind_r8, i8=>shr_kind_i8, cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_const_mod    , only : shr_const_spval, shr_const_pi
  use shr_sys_mod      , only : shr_sys_abort
  use shr_cal_mod      , only : shr_cal_noleap, shr_cal_gregorian, shr_cal_ymd2date, shr_cal_ymd2julian
  use shr_mpi_mod      , only : shr_mpi_bcast
  use shr_orb_mod      , only : shr_orb_params, SHR_ORB_UNDEF_INT, SHR_ORB_UNDEF_REAL
  use dshr_strdata_mod , only : shr_strdata_type, shr_strdata_readnml
  use dshr_methods_mod , only : chkerr, state_setscalar,  state_diagnose, memcheck
  use dshr_methods_mod , only : set_component_logging, log_clock_advance
  use dshr_nuopc_mod   , only : dshr_advertise, dshr_model_initphase, dshr_set_runclock
  use dshr_nuopc_mod   , only : dshr_sdat_init
  use dshr_nuopc_mod   , only : dshr_restart_read, dshr_restart_write
  use dshr_nuopc_mod   , only : dshr_create_mesh_from_grid
  use datm_comp_mod    , only : datm_comp_advertise, datm_comp_realize, datm_comp_run
  use perf_mod         , only : t_startf, t_stopf, t_barrierf

  implicit none
  private ! except

  public  :: SetServices

  private :: InitializeAdvertise
  private :: InitializeRealize
  private :: ModelAdvance
  private :: ModelFinalize
  private :: OrbitalInit
  private :: OrbitalUpdate

  interface getNextRadCday
     module procedure getNextRadCDay_i8
     module procedure getNextRadCDay_i4
  end interface getNextRadCday

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  type(shr_strdata_type)       :: sdat
  character(len=128)           :: flds_scalar_name = ''
  integer                      :: flds_scalar_num = 0
  integer                      :: flds_scalar_index_nx = 0
  integer                      :: flds_scalar_index_ny = 0
  integer                      :: flds_scalar_index_nextsw_cday = 0
  type(ESMF_Mesh)              :: mesh
  integer                      :: compid                    ! mct comp id
  integer                      :: mpicom                    ! mpi communicator
  integer                      :: my_task                   ! my task in mpi communicator mpicom
  character(len=16)            :: inst_suffix = ""          ! char string associated with instance (ie. "_0001" or "")
  integer                      :: logunit                   ! logging unit number
  character(len=*) , parameter :: nullstr = 'undefined'
  integer          , parameter :: master_task = 0           ! task number of master task
  character(len=*) , parameter :: rpfile = 'rpointer.atm'
  character(*)     , parameter :: modName = "(atm_comp_nuopc)"

  integer                      :: iradsw          ! radiation interval (input namelist)
  integer                      :: idt             ! integer model timestep
  character(len=CL)            :: orb_mode        ! attribute - orbital mode (nuopc attribute)
  integer                      :: orb_iyear       ! attribute - orbital year (nuopc attribute)
  integer                      :: orb_iyear_align ! attribute - associated with model year (nuopc attribute)
  real(R8)                     :: orb_obliq       ! attribute - obliquity in degrees (nuopc attribute)
  real(R8)                     :: orb_mvelp       ! attribute - moving vernal equinox longitude (nuopc attribute)
  real(R8)                     :: orb_eccen       ! attribute and update-  orbital eccentricity (nuopc attribute)
  character(CL)                :: rest_file       ! restart filename
  character(CL)                :: rest_file_strm  ! restart filename for streams

  character(len=*) , parameter :: orb_fixed_year       = 'fixed_year'
  character(len=*) , parameter :: orb_variable_year    = 'variable_year'
  character(len=*) , parameter :: orb_fixed_parameters = 'fixed_parameters'

  logical :: diagnose_data = .false.

  character(*)     , parameter :: u_FILE_u = &
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
    integer           :: inst_index            ! number of current instance (ie. 1)
    character(len=CL) :: cvalue                ! temporary
    integer           :: shrlogunit            ! original log unit
    integer           :: nu                    ! unit number
    integer           :: ierr                  ! error code
    character(CL)     :: decomp                ! decomp strategy - not used for NUOPC - but still needed in namelist for now
    logical           :: flds_co2a             ! use case
    logical           :: flds_co2b             ! use case
    logical           :: flds_co2c             ! use case
    logical           :: flds_wiso             ! use case
    character(len=CL) :: fileName              ! generic file name
    logical           :: force_prognostic_true ! if true set prognostic true
    logical           :: wiso_datm = .false.   ! expect isotopic forcing from file?
    logical           :: presaero              ! true => send valid prescribe aero fields to coupler
    character(CL)     :: factorFn              ! file containing correction factors
    character(CL)     :: anomaly_forcing(8)    ! true => send anomaly forcing fields to coupler (not used here)
    character(CL)     :: bias_correct          ! true => send bias correction fields to coupler (not used here)
    logical           :: get_importdata        ! data is sent back to datm
    character(CL)     :: restfilm = nullstr    ! model restart file namelist
    character(CL)     :: restfils = nullstr    ! stream restart file namelist
    character(len=*),parameter :: subname=trim(modName)//':(InitializeAdvertise) '
    !-------------------------------------------------------------------------------

    namelist / datm_nml / iradsw, factorFn, restfilm, restfils, &
         presaero, bias_correct, anomaly_forcing, force_prognostic_true, wiso_datm

    rc = ESMF_SUCCESS

    ! Obtain flds_scalar values, mpi values and multi-instance values
    call dshr_advertise(gcomp, mpicom, my_task,  inst_index, inst_suffix, &
         flds_scalar_name, flds_scalar_num, flds_scalar_index_nx, flds_scalar_index_ny, rc)

    ! Set logunit and set shr logging to my log file
    call set_component_logging(gcomp, my_task==master_task, logunit, shrlogunit, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set input namelist filename
    filename = "datm_in"//trim(inst_suffix)

    ! Read atm_nml from filename
    iradsw = 0
    factorFn = 'null'
    restfilm = trim(nullstr)
    restfils = trim(nullstr)
    presaero = .false.
    force_prognostic_true = .false.
    if (my_task == master_task) then
       open (newunit=nu,file=trim(filename),status="old",action="read")
       read (nu,nml=datm_nml,iostat=ierr)
       close(nu)
       if (ierr > 0) then
          write(logunit,*) 'ERROR: reading input namelist, '//trim(filename)//' iostat=',ierr
          call shr_sys_abort(subName//': namelist read error '//trim(filename))
       end if
       write(logunit,*)' iradsw   = ',iradsw
       write(logunit,*)' factorFn = ',trim(factorFn)
       write(logunit,*)' restfilm = ',trim(restfilm)
       write(logunit,*)' restfils = ',trim(restfils)
       write(logunit,*)' presaero = ',presaero
       write(logunit,*)' force_prognostic_true = ',force_prognostic_true
       write(logunit,*)' wiso_datm   = ',wiso_datm
    endif
    call shr_mpi_bcast(iradsw                ,mpicom, 'iradsw')
    call shr_mpi_bcast(factorFn              ,mpicom, 'factorFn')
    call shr_mpi_bcast(restfilm              ,mpicom, 'restfilm')
    call shr_mpi_bcast(restfils              ,mpicom, 'restfils')
    call shr_mpi_bcast(presaero              ,mpicom, 'presaero')
    call shr_mpi_bcast(wiso_datm             ,mpicom, 'wiso_datm')
    call shr_mpi_bcast(force_prognostic_true ,mpicom, 'force_prognostic_true')
    rest_file = trim(restfilm)
    rest_file_strm = trim(restfils)

    ! Read shr_strdata_nml from filename
    ! Read sdat namelist (need to do this here in order to get the datamode value - which
    ! is needed or order to do the advertise phase
    call shr_strdata_readnml(sdat, trim(filename), mpicom=mpicom)

    ! Validate sdat%datamode
    if (trim(sdat%datamode) == 'NULL'      .or. trim(sdat%datamode) == 'CORE2_NYF' .or. &
        trim(sdat%datamode) == 'CORE2_IAF' .or. trim(sdat%datamode) == 'CORE_IAF_JRA' .or. &
        trim(sdat%datamode) == 'CLMNCEP'   .or. trim(sdat%datamode) == 'COPYALL'   ) then
       if (my_task == master_task) then
          write(logunit,*) ' datm sdat%datamode = ',trim(sdat%datamode)
       end if
    else
       call shr_sys_abort(' ERROR illegal datm sdat%datamode = '//trim(sdat%datamode))
    endif

    if (sdat%datamode /= 'NULL') then
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

       call NUOPC_CompAttributeGet(gcomp, name='flds_wiso', value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) flds_wiso
       call ESMF_LogWrite('flds_wiso = '// trim(cvalue), ESMF_LOGMSG_INFO)

       if (force_prognostic_true) then
          get_importdata = .true.
       else
          get_importdata = .false.
       endif

       ! check that flds_wiso matches wiso_datm
       ! TODO: should not need wiso_datm for nuopc data models
       if (wiso_datm /= flds_wiso) then
          call shr_sys_abort(subName//': datm namelist wiso_datm must match nuopc attribute flds_wiso')
       end if

       call datm_comp_advertise(importState, exportState, flds_scalar_name, &
            get_importdata, presaero, flds_wiso, flds_co2a, flds_co2b, flds_co2c, factorfn, rc=rc)
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
    character(CL)           :: filename
    type(ESMF_TimeInterval) :: timeStep
    type(ESMF_TIME)         :: currTime
    type(ESMF_Calendar)     :: esmf_calendar ! esmf calendar
    type(ESMF_CalKind_Flag) :: esmf_caltype  ! esmf calendar type
    integer                 :: current_ymd   ! model date
    integer                 :: current_year  ! model year
    integer                 :: current_mon   ! model month
    integer                 :: current_day   ! model day
    integer                 :: current_tod   ! model sec into model date
    integer(i8)             :: stepno        ! step number
    real(r8)                :: nextsw_cday   ! calendar of next atm sw
    character(len=256)      :: cvalue        ! character string for input config
    integer                 :: shrlogunit    ! original log unit
    logical                 :: read_restart  ! start from restart
    logical                 :: scmMode       ! single column mode
    real(R8)                :: scmLat        ! single column lat
    real(R8)                :: scmLon        ! single column lon
    real(R8)                :: orbEccen      ! orb eccentricity (unit-less)
    real(R8)                :: orbMvelpp     ! orb moving vernal eq (radians)
    real(R8)                :: orbLambm0     ! orb mean long of perhelion (radians)
    real(R8)                :: orbObliqr     ! orb obliquity (radians)
    logical                 :: isPresent, isSet
    character(len=*), parameter :: subname=trim(modName)//':(InitializeRealize) '
    !-------------------------------------------------------------------------------

    if (sdat%datamode == 'NULL') RETURN

    ! TODO: read_restart, scmlat, scmlon, orbeccen, orbmvelpp, orblambm0, orbobliqr needs to be obtained
    ! from the config attributes of the gridded component

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! Reset shr logging to my log file
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_setLogUnit (logUnit)

    ! Create the data model mesh
    call NUOPC_CompAttributeGet(gcomp, name='mesh_atm', value=filename, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (trim(filename) == 'create_mesh') then
       ! get the datm grid from the domain file
       call NUOPC_CompAttributeGet(gcomp, name='domain_atm', value=filename, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call dshr_create_mesh_from_grid(trim(filename), mesh, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       mesh = ESMF_MeshCreate(trim(filename), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    if (my_task == master_task) then
       write(logunit,*) trim(subname)// " obtaining datm mesh from " // trim(filename)
    end if

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
    call t_startf('datm_strdata_init')
    call dshr_sdat_init(mpicom, compid, my_task, master_task, logunit, &
         scmmode, scmlon, scmlat, clock, mesh, 'datm', sdat, use_new=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('datm_strdata_init')

    ! Realize the actively coupled fields, now that a mesh is established and
    ! initialize dfields data type (to map streams to export state fields)
    call datm_comp_realize(sdat, importState, exportState, flds_scalar_name, flds_scalar_num, mesh, &
         logunit, my_task==master_task, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Read restart if necessary
    call NUOPC_CompAttributeGet(gcomp, name='read_restart', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) read_restart
    if (read_restart) then
       call dshr_restart_read(rest_file, rest_file_strm, rpfile, inst_suffix, nullstr, &
            logunit, my_task, master_task, mpicom, sdat)
    end if

    ! Get the time to interpolate the stream data to
    call ESMF_ClockGet( clock, currTime=currTime, timeStep=timeStep, advanceCount=stepno, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(currTime, yy=current_year, mm=current_mon, dd=current_day, s=current_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(current_year, current_mon, current_day, current_ymd)

    ! Get model timestep
    call ESMF_TimeIntervalGet( timeStep, s=idt, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Initialize and update orbital values
    call OrbitalInit(gcomp, logunit, my_task == master_task, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call OrbitalUpdate(clock, logunit, my_task == master_task, orbEccen, orbObliqr, orbLambm0, orbMvelpp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Run datm
    call datm_comp_run(mpicom, compid, my_task==master_task, logunit, current_ymd, current_tod, sdat, &
         mesh, current_mon, orbEccen, orbMvelpp, orbLambm0, orbObliqr, idt, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Add scalars to export state
    call State_SetScalar(dble(sdat%nxg), flds_scalar_index_nx, exportState, flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_SetScalar(dble(sdat%nyg),flds_scalar_index_ny, exportState,  flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldIdxNextSwCday", value=cvalue, &
         isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read (cvalue,*) flds_scalar_index_nextsw_cday
    else
       call shr_sys_abort(subname//'Need to set attribute ScalarFieldIdxNextSwCday')
    endif
    nextsw_cday = getNextRadCDay( current_ymd, current_tod, stepno, idt, iradsw, sdat%calendar )
    call State_SetScalar(nextsw_cday, flds_scalar_index_nextsw_cday, exportState, flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Diagnostics
    if (diagnose_data) then
       call State_diagnose(exportState, subname//':ES',rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Reset shr logging to original values
    call shr_file_setLogUnit (shrlogunit)

  end subroutine InitializeRealize

  !===============================================================================

  subroutine ModelAdvance(gcomp, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)        :: clock
    type(ESMF_State)        :: importState, exportState
    type(ESMF_Time)         :: time
    type(ESMF_Alarm)        :: alarm
    type(ESMF_Time)         :: currTime
    type(ESMF_Time)         :: nextTime
    type(ESMF_TimeInterval) :: timeStep
    integer                 :: shrlogunit    ! original log unit
    real(r8)                :: nextsw_cday
    logical                 :: write_restart ! restart alarm is ringing
    integer                 :: next_ymd      ! model date
    integer                 :: next_tod      ! model sec into model date
    integer                 :: yr, mon, day  ! year, month, day
    integer(I8)             :: stepno        ! step number
    real(R8)                :: orbEccen      ! orb eccentricity (unit-less)
    real(R8)                :: orbMvelpp     ! orb moving vernal eq (radians)
    real(R8)                :: orbLambm0     ! orb mean long of perhelion (radians)
    real(R8)                :: orbObliqr     ! orb obliquity (radians)
    character(len=CL)       :: case_name     ! case name
    character(len=CL)       :: cvalue        ! temporary
    character(len=*),parameter  :: subname=trim(modName)//':(ModelAdvance) '
    !-------------------------------------------------------------------------------

    if (sdat%datamode == 'NULL') RETURN

    rc = ESMF_SUCCESS

    call t_startf(subname)
    call memcheck(subname, 5, my_task==master_task)

    ! Reset shr logging to my log file
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_setLogUnit (logunit)

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
    call OrbitalUpdate(clock, logunit, my_task == master_task, orbEccen, orbObliqr, orbLambm0, orbMvelpp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Run datm
    call t_startf('datm_run')
    call datm_comp_run(mpicom, compid, my_task==master_task, logunit, next_ymd, next_tod, sdat, &
         mesh, mon, orbEccen, orbMvelpp, orbLambm0, orbObliqr, idt, rc)
    call t_stopf('datm_run')

    ! Update nextsw_cday for scalar data
    ! Use nextYMD and nextTOD here since since the component - clock is advance at the END of the time interval
    nextsw_cday = getNextRadCDay( next_ymd, next_tod, stepno, idt, iradsw, sdat%calendar )
    call State_SetScalar(nextsw_cday, flds_scalar_index_nextsw_cday, exportState, flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Write_restart if alarm is ringing
    call ESMF_ClockGetAlarm(clock, alarmname='alarm_restart', alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call t_startf('datm_restart')
       call NUOPC_CompAttributeGet(gcomp, name='case_name', value=case_name, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call dshr_restart_write(rpfile, case_name, 'datm', inst_suffix, next_ymd, next_tod, &
            logunit, mpicom, my_task, master_task, sdat)
       call t_stopf('datm_restart')
    endif

    ! Diagnostics
    if (diagnose_data) then
       call State_diagnose(exportState,subname//':ES',rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

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
       write(logunit,*) 'datm : end of main integration loop'
       write(logunit,*)
    end if

  end subroutine ModelFinalize

  !===============================================================================

  subroutine OrbitalInit(gcomp, logunit, mastertask, rc)

    !----------------------------------------------------------
    ! Initialize orbital related values
    !----------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp)                 :: gcomp
    integer             , intent(in)    :: logunit
    logical             , intent(in)    :: mastertask
    integer             , intent(out)   :: rc              ! output error

    ! local variables
    character(len=CL) :: msgstr          ! temporary
    character(len=CL) :: cvalue          ! temporary
    character(len=*) , parameter :: subname = "(OrbitalInit)"
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Determine orbital attributes from input
    call NUOPC_CompAttributeGet(gcomp, name="orb_mode", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_mode
    call NUOPC_CompAttributeGet(gcomp, name="orb_iyear", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_iyear
    call NUOPC_CompAttributeGet(gcomp, name="orb_iyear_align", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_iyear_align
    call NUOPC_CompAttributeGet(gcomp, name="orb_obliq", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_obliq
    call NUOPC_CompAttributeGet(gcomp, name="orb_eccen", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_eccen
    call NUOPC_CompAttributeGet(gcomp, name="orb_mvelp", value=cvalue, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) orb_mvelp

    ! Error checks
    if (trim(orb_mode) == trim(orb_fixed_year)) then
       if (orb_iyear == SHR_ORB_UNDEF_INT) then
          if (mastertask) then
             write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(orb_mode)
             write(logunit,*) trim(subname),' ERROR: fixed_year settings = ',orb_iyear
             write (msgstr, *) ' ERROR: invalid settings for orb_mode '//trim(orb_mode)
          end if
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          return  ! bail out
       else
          orb_obliq = SHR_ORB_UNDEF_REAL
          orb_eccen = SHR_ORB_UNDEF_REAL
          orb_mvelp = SHR_ORB_UNDEF_REAL
       endif
    elseif (trim(orb_mode) == trim(orb_variable_year)) then
       if (orb_iyear == SHR_ORB_UNDEF_INT .or. orb_iyear_align == SHR_ORB_UNDEF_INT) then
          if (mastertask) then
             write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(orb_mode)
             write(logunit,*) trim(subname),' ERROR: variable_year settings = ',orb_iyear, orb_iyear_align
             write (msgstr, *) subname//' ERROR: invalid settings for orb_mode '//trim(orb_mode)
          end if
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          return  ! bail out
       else
          orb_obliq = SHR_ORB_UNDEF_REAL
          orb_eccen = SHR_ORB_UNDEF_REAL
          orb_mvelp = SHR_ORB_UNDEF_REAL
       endif
    elseif (trim(orb_mode) == trim(orb_fixed_parameters)) then
       !-- force orb_iyear to undef to make sure shr_orb_params works properly
       if (orb_eccen == SHR_ORB_UNDEF_REAL .or. orb_obliq == SHR_ORB_UNDEF_REAL .or. orb_mvelp == SHR_ORB_UNDEF_REAL) then
          if (mastertask) then
             write(logunit,*) trim(subname),' ERROR: invalid settings orb_mode =',trim(orb_mode)
             write(logunit,*) trim(subname),' ERROR: orb_eccen = ',orb_eccen
             write(logunit,*) trim(subname),' ERROR: orb_obliq = ',orb_obliq
             write(logunit,*) trim(subname),' ERROR: orb_mvelp = ',orb_mvelp
             write (msgstr, *) subname//' ERROR: invalid settings for orb_mode '//trim(orb_mode)
          end if
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          return  ! bail out
       else
          orb_iyear       = SHR_ORB_UNDEF_INT
          orb_iyear_align = SHR_ORB_UNDEF_INT
       endif
    else
       write (msgstr, *) subname//' ERROR: invalid orb_mode '//trim(orb_mode)
       call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
       rc = ESMF_FAILURE
       return  ! bail out
    endif

  end subroutine OrbitalInit

  !===============================================================================

  subroutine OrbitalUpdate(clock, logunit,  mastertask, eccen, obliqr, lambm0, mvelpp, rc)

    !----------------------------------------------------------
    ! Update orbital settings
    !----------------------------------------------------------

    ! input/output variables
    type(ESMF_Clock) , intent(in)    :: clock
    integer          , intent(in)    :: logunit
    logical          , intent(in)    :: mastertask
    real(R8)         , intent(inout) :: eccen  ! orbital eccentricity
    real(R8)         , intent(inout) :: obliqr ! Earths obliquity in rad
    real(R8)         , intent(inout) :: lambm0 ! Mean long of perihelion at vernal equinox (radians)
    real(R8)         , intent(inout) :: mvelpp ! moving vernal equinox longitude of perihelion plus pi (radians)
    integer          , intent(out)   :: rc     ! output error

    ! local variables
    type(ESMF_Time)   :: CurrTime ! current time
    integer           :: year     ! model year at current time
    integer           :: orb_year ! orbital year for current orbital computation
    character(len=CL) :: msgstr   ! temporary
    logical           :: lprint
    logical           :: first_time = .true.
    character(len=*) , parameter :: subname = "(OrbitalUpdate)"
    !-------------------------------------------

    if (trim(orb_mode) == trim(orb_variable_year)) then
       call ESMF_ClockGet(clock, CurrTime=CurrTime, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeGet(CurrTime, yy=year, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       orb_year = orb_iyear + (year - orb_iyear_align)
       lprint = mastertask
    else
       orb_year = orb_iyear
       if (first_time) then
          lprint = mastertask
          first_time = .false.
       else
          lprint = .false.
       end if
    end if

    eccen = orb_eccen
    call shr_orb_params(orb_year, eccen, orb_obliq, orb_mvelp, obliqr, lambm0, mvelpp, lprint)

    if ( eccen  == SHR_ORB_UNDEF_REAL .or. obliqr == SHR_ORB_UNDEF_REAL .or. &
         mvelpp == SHR_ORB_UNDEF_REAL .or. lambm0 == SHR_ORB_UNDEF_REAL) then
       write (msgstr, *) subname//' ERROR: orb params incorrect'
       call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
       return  ! bail out
    endif

  end subroutine OrbitalUpdate

  !===============================================================================
  real(R8) function getNextRadCDay_i8( ymd, tod, stepno, dtime, iradsw, calendar )

    !  Return the calendar day of the next radiation time-step.
    !  General Usage: nextswday = getNextRadCDay(curr_date)

    use shr_kind_mod   , only : r8=>shr_kind_r8, i8=>shr_kind_i8, cs=>shr_kind_cs, cl=>shr_kind_cl
    use shr_cal_mod    , only : shr_cal_date2julian
    use shr_const_mod  , only : shr_const_cday

    ! input/output variables
    integer    , intent(in)    :: ymd
    integer    , intent(in)    :: tod
    integer(I8), intent(in)    :: stepno
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
    getNextRadCDay_i8 = nextsw_cday

  end function getNextRadCDay_i8

  !===============================================================================
  real(R8) function getNextRadCDay_i4( ymd, tod, stepno, dtime, iradsw, calendar )

    !  Return the calendar day of the next radiation time-step.
    !  General Usage: nextswday = getNextRadCDay(curr_date)

    use shr_kind_mod   , only : r8=>shr_kind_r8, i8=>shr_kind_i8, cs=>shr_kind_cs, cl=>shr_kind_cl
    use shr_cal_mod    , only : shr_cal_date2julian
    use shr_const_mod  , only : shr_const_cday

    ! input/output variables
    integer    , intent(in)    :: ymd
    integer    , intent(in)    :: tod
    integer    , intent(in)    :: stepno
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
    getNextRadCDay_i4 = nextsw_cday

  end function getNextRadCDay_i4

end module atm_comp_nuopc
