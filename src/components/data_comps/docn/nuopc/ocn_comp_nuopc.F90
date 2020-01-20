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
  use shr_file_mod     , only : shr_file_getlogunit, shr_file_setlogunit
  use shr_kind_mod     , only : r8=>shr_kind_r8, i8=>shr_kind_i8, cl=>shr_kind_cl, cs=>shr_kind_cs
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
  use docn_comp_mod    , only : docn_comp_advertise, docn_comp_init, docn_comp_run 
  use docn_comp_mod    , only : fldsExport, fldsExport_num, fldsImport, fldsImport_num
  use docn_comp_mod    , only : somtp  ! for restart
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

  type(shr_strdata_type)   :: sdat
  character(len=CS)        :: flds_scalar_name = ''
  integer                  :: flds_scalar_num = 0
  integer                  :: flds_scalar_index_nx = 0
  integer                  :: flds_scalar_index_ny = 0
  type(ESMF_Mesh)          :: mesh                            ! model mesh
  integer                  :: compid                          ! mct comp id
  integer                  :: mpicom                          ! mpi communicator
  integer                  :: my_task                         ! my task in mpi communicator mpicom
  character(len=16)        :: inst_suffix = ""                ! char string associated with instance (ie. "_0001" or "")
  integer                  :: logunit                         ! logging unit number
  logical                  :: read_restart                    ! start from restart
  character(*) , parameter :: nullstr = 'undefined'
  integer      , parameter :: master_task=0                   ! task number of master task
  character(*) , parameter :: rpfile = 'rpointer.ocn'
  character(*) , parameter :: modName = "(ocn_comp_nuopc)"

  ! docn_in namelist input
  real(R8)                 :: sst_constant_value
  integer                  :: aquap_option
  character(CL)            :: rest_file                       ! restart filename
  character(CL)            :: rest_file_strm                  ! restart filename for streams

  character(*) , parameter :: u_FILE_u = &
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
    character(len=CL) :: cvalue             ! temporary
    integer           :: shrlogunit         ! original log unit
    character(len=CL) :: fileName           ! generic file name
    integer           :: nu                 ! unit number
    integer           :: ierr               ! error code
    character(CL)     :: decomp             ! decomp strategy - not used for NUOPC - but still needed in namelist for now
    character(CL)     :: restfilm = nullstr ! model restart file namelist
    character(CL)     :: restfils = nullstr ! stream restart file namelist
    logical           :: ocn_prognostic     ! true => ocn expects import data
    logical           :: force_prognostic_true = .false. ! if true set prognostic true
    character(len=*),parameter  :: subname=trim(modName)//':(InitializeAdvertise) '
    !-------------------------------------------------------------------------------

    namelist / docn_nml / decomp, restfilm, restfils, force_prognostic_true, sst_constant_value

    rc = ESMF_SUCCESS

    ! obtain flds_scalar values, mpi values and multi-instance values
    call dshr_advertise(gcomp, mpicom, my_task,  inst_index, inst_suffix, &
         flds_scalar_name, flds_scalar_num, flds_scalar_index_nx, flds_scalar_index_ny, rc)

    ! set logunit and set shr logging to my log file
    call set_component_logging(gcomp, my_task==master_task, logunit, shrlogunit, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! set input namelist filename
    filename = "docn_in"//trim(inst_suffix)

    ! Read docn_nml from filename
    if (my_task == master_task) then
       open (newunit=nu,file=trim(filename),status="old",action="read")
       read (nu,nml=docn_nml,iostat=ierr)
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
    call shr_mpi_bcast(force_prognostic_true, mpicom, 'force_prognostic_true')
    call shr_mpi_bcast(sst_constant_value   , mpicom, 'sst_constant_value')
    rest_file      = trim(restfilm)
    rest_file_strm = trim(restfils)

    ! Read shr_strdata_nml from filename
    ! Read sdat namelist (need to do this here in order to get the datamode value - which
    ! is needed or order to do the advertise phase
    call shr_strdata_readnml(sdat, trim(filename), mpicom=mpicom)

    ! Validate sdat%datamode
    if ( trim(sdat%datamode) == 'NULL'          .or. trim(sdat%datamode) == 'COPYALL' .or.&
         trim(sdat%datamode) == 'SSTDATA'       .or. &
         trim(sdat%datamode) == 'SST_AQUAPANAL' .or. trim(sdat%datamode) == 'SST_AQUAPFILE' .or. &
         trim(sdat%datamode) == 'IAF'           .or. &
         trim(sdat%datamode) == 'SOM'           .or. trim(sdat%datamode) == 'SOM_AQUAP') then
       if (my_task == master_task) then
          write(logunit,*) ' docn sdat%datamode = ',trim(sdat%datamode)
       end if
    else
       call shr_sys_abort(' ERROR illegal docn datamode = '//trim(sdat%datamode))
    endif

    if (trim(sdat%datamode) /= 'NULL') then
       ! determine if ocn will receive import data
       if ( force_prognostic_true .or. trim(sdat%datamode) == 'IAF' .or. &
            trim(sdat%datamode) == 'SOM' .or. trim(sdat%datamode) == 'SOM_AQUAP') then
          ocn_prognostic = .true.
       else
          ocn_prognostic = .false.
       end if

       ! Special logic for prescribed aquaplanet
       if (sdat%datamode(1:9) == 'SST_AQUAP' .and. trim(sdat%datamode) /= 'SST_AQUAPFILE') then
          ! First determine the prescribed aquaplanet option
          if (len_trim(sdat%datamode) == 10) then
             read(sdat%datamode(10:10),'(i1)') aquap_option
          else if (len_trim(sdat%datamode) == 11) then
             read(sdat%datamode(10:11),'(i2)') aquap_option
          end if
          ! Now remove the index from the sdat%datamode value, to have a generic setting for later use
          sdat%datamode = "SST_AQUAPANAL"
       end if

       call docn_comp_advertise(importState, exportState, flds_scalar_name, ocn_prognostic, rc=rc)
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
    type(ESMF_TimeInterval) :: TimeStep
    type(ESMF_Time)         :: currTime
    integer                 :: current_ymd  ! model date
    integer                 :: current_year ! model year
    integer                 :: current_mon  ! model month
    integer                 :: current_day  ! model day
    integer                 :: current_tod  ! model sec into model date
    character(CL)           :: cvalue       ! temporary
    integer                 :: shrlogunit   ! original log unit
    integer                 :: n,k          ! generic counters
    logical                 :: scmMode      ! single column mode
    real(R8)                :: scmLat       ! single column lat
    real(R8)                :: scmLon       ! single column lon
    character(CS)           :: model_name
    integer                 :: model_dt     ! integer model timestep
    real(R8)                :: dt           ! real model timestep
    logical                 :: reset_domain_mask ! true => reset the domain mask for aquaplanet
    integer, allocatable, target :: gindex(:)
    character(len=*), parameter :: F00   = "('ocn_comp_nuopc: ')',8a)"
    character(len=*), parameter :: subname=trim(modName)//':(InitializeRealize) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !--------------------------------
    ! Reset shr logging to my log file
    !--------------------------------

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_setLogUnit (logUnit)

    !--------------------------------
    ! Create the data model mesh
    !--------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='mesh_ocn', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (my_task == master_task) then
       write(logunit,*) ' obtaining mesh_ocn from '//trim(cvalue)
    end if

    mesh = ESMF_MeshCreate(filename=trim(cvalue), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Initialize sdat
    !--------------------------------

    call t_startf('docn_strdata_init')
    call NUOPC_CompAttributeGet(gcomp, name='MCTID', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) compid

    ! set single column values
    call NUOPC_CompAttributeGet(gcomp, name='scmlon', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scmlon
    call NUOPC_CompAttributeGet(gcomp, name='scmlat', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scmlat
    call NUOPC_CompAttributeGet(gcomp, name='single_column', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) scmMode

    ! for aquaplanet mode - need to reset the domain
    if (sdat%datamode == 'SST_AQUAPANAL' .or. sdat%datamode == 'SST_AQUAPFILE' .or. sdat%datamode == 'SOM_AQUAP') then
       reset_domain_mask = .true.       
    else
       reset_domain_mask= .false.
    end if

    ! determine the model name
    model_name = 'docn'

    call dshr_sdat_init(mpicom, compid, my_task, master_task, logunit, &
         scmmode, scmlon, scmlat, clock, mesh, model_name, sdat, reset_domain_mask=reset_domain_mask, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (my_task == master_task) write(logunit,*) ' initialized SDAT'
    call t_stopf('docn_strdata_init')

    !--------------------------------
    ! Check that mesh lats and lons correspond to those on the input domain file
    !--------------------------------

    call dshr_check_mesh(mesh, sdat, 'docn', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! realize the actively coupled fields, now that a mesh is established
    !--------------------------------

    ! NUOPC_Realize "realizes" a previously advertised field in the importState and exportState
    ! by replacing the advertised fields with the newly created fields of the same name.

    call dshr_realize( exportState, fldsExport, fldsExport_num, &
         flds_scalar_name, flds_scalar_num, mesh, tag=subname//':docnExport', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call dshr_realize( importState, fldsImport, fldsImport_num, &
         flds_scalar_name, flds_scalar_num, mesh, tag=subname//':docnImport', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Initialize dfields data type (to map streams to export state fields)
    !--------------------------------

    call docn_comp_init(sdat, importState, exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Read restart if necessary
    !--------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='read_restart', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) read_restart
    if (read_restart) then
       call dshr_restart_read(rest_file, rest_file_strm, rpfile, inst_suffix, nullstr, &
            logunit, my_task, master_task, mpicom, sdat, fld=somtp, fldname='somtp')
    end if

    !--------------------------------
    ! Run docn to create export state
    !--------------------------------

    ! get the time to interpolate the stream data to
    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(currTime, yy=current_year, mm=current_mon, dd=current_day, s=current_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(current_year, current_mon, current_day, current_ymd)

    ! get model timestep
    call ESMF_ClockGet( clock, timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeIntervalGet( timeStep, s=model_dt, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    dt = model_dt * 1.0_r8

    ! run docn
    call docn_comp_run(mpicom, my_task, master_task, logunit, current_ymd, current_tod, &
         sdat, mesh, sst_constant_value, dt, aquap_option, read_restart, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! add scalars to export state
    call State_SetScalar(dble(sdat%nxg),flds_scalar_index_nx, exportState, flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call State_SetScalar(dble(sdat%nyg),flds_scalar_index_ny, exportState, flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! diagnostics
    call State_diagnose(exportState,subname//':ES',rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_file_setLogUnit (shrlogunit)

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
    integer                 :: next_ymd      ! model date
    integer                 :: next_tod      ! model sec into model date
    integer                 :: yr            ! year
    integer                 :: mon           ! month
    integer                 :: day           ! day in month
    character(CL)           :: case_name     ! case name
    integer                 :: model_dt      ! integer model timestep
    real(R8)                :: dt            ! real model timestep
    character(len=*),parameter :: subname=trim(modName)//':(ModelAdvance) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call memcheck(subname, 5, my_task == master_task)

    ! Reset shr logging to my log file
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_setLogUnit (logunit)

    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! For nuopc - the component clock is advanced at the end of the time interval
    ! For these to match for now - need to advance nuopc one timestep ahead for
    ! shr_strdata time interpolation
    call ESMF_ClockGet( clock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    nextTime = currTime + timeStep
    call ESMF_TimeGet( nextTime, yy=yr, mm=mon, dd=day, s=next_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(yr, mon, day, next_ymd)

    ! run docn
    call docn_comp_run(mpicom, my_task, master_task, logunit, next_ymd, next_tod, &
         sdat, mesh, sst_constant_value, dt, aquap_option, read_restart, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! write_restart if alarm is ringing
    call ESMF_ClockGetAlarm(clock, alarmname='alarm_restart', alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call t_startf('docn_restart')
       call NUOPC_CompAttributeGet(gcomp, name='case_name', value=case_name, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call dshr_restart_write(rpfile, case_name, 'docn', inst_suffix, next_ymd, next_tod, &
            logunit, mpicom, my_task, master_task, sdat, fld=somtp, fldname='somtp')
       call t_stopf('docn_restart')
    endif

    ! write diagnostics
    call State_diagnose(exportState,subname//':ES',rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (my_task == master_task) then
       call log_clock_advance(clock, 'docn', logunit, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    call shr_file_setLogUnit (shrlogunit)

  end subroutine ModelAdvance

  !===============================================================================

  subroutine ModelFinalize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(*), parameter :: F00   = "('(docn_comp_final) ',8a)"
    character(*), parameter :: F91   = "('(docn_comp_final) ',73('-'))"
    character(len=*),parameter  :: subname=trim(modName)//':(ModelFinalize) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
    if (my_task == master_task) then
       write(logunit,F91)
       write(logunit,F00) 'docn : end of main integration loop'
       write(logunit,F91)
    end if
    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelFinalize

end module ocn_comp_nuopc

