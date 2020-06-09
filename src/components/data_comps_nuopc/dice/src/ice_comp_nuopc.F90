module ice_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for DICE
  !----------------------------------------------------------------------------

  use ESMF
  use NUOPC                , only : NUOPC_CompDerive, NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize
  use NUOPC                , only : NUOPC_CompAttributeGet, NUOPC_Advertise
  use NUOPC_Model          , only : model_routine_SS        => SetServices
  use NUOPC_Model          , only : model_label_Advance     => label_Advance
  use NUOPC_Model          , only : model_label_SetRunClock => label_SetRunClock
  use NUOPC_Model          , only : model_label_Finalize    => label_Finalize
  use NUOPC_Model          , only : NUOPC_ModelGet
  use shr_kind_mod         , only : r8=>shr_kind_r8, cxx=>shr_kind_cxx, cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_const_mod        , only : shr_const_pi
  use shr_sys_mod          , only : shr_sys_abort
  use shr_cal_mod          , only : shr_cal_ymd2date, shr_cal_ymd2julian
  use shr_mpi_mod          , only : shr_mpi_bcast
  use dshr_mod             , only : dshr_model_initphase, dshr_init, dshr_mesh_init
  use dshr_mod             , only : dshr_state_setscalar, dshr_set_runclock, dshr_log_clock_advance
  use dshr_methods_mod     , only : dshr_state_diagnose, chkerr, memcheck
  use dshr_strdata_mod     , only : shr_strdata_type, shr_strdata_init_from_xml, shr_strdata_advance
  use dshr_dfield_mod      , only : dfield_type, dshr_dfield_add, dshr_dfield_copy
  use dshr_fldlist_mod     , only : fldlist_type, dshr_fldlist_add, dshr_fldlist_realize
  use perf_mod             , only : t_startf, t_stopf, t_adj_detailf, t_barrierf

  use dice_datamode_ssmi_mod , only : dice_datamode_ssmi_advertise
  use dice_datamode_ssmi_mod , only : dice_datamode_ssmi_init_pointers
  use dice_datamode_ssmi_mod , only : dice_datamode_ssmi_advance
  use dice_datamode_ssmi_mod , only : dice_datamode_ssmi_restart_read
  use dice_datamode_ssmi_mod , only : dice_datamode_ssmi_restart_write

  implicit none
  private ! except

  public  :: SetServices

  private :: InitializeAdvertise
  private :: InitializeRealize
  private :: ModelAdvance
  private :: dice_comp_run
  private :: ModelFinalize

  !--------------------------------------------------------------------------
  ! Module data
  !--------------------------------------------------------------------------

  type(shr_strdata_type)       :: sdat
  type(ESMF_Mesh)              :: model_mesh
  character(len=CS)            :: flds_scalar_name = ''
  integer                      :: flds_scalar_num = 0
  integer                      :: flds_scalar_index_nx = 0
  integer                      :: flds_scalar_index_ny = 0
  integer                      :: compid                              ! mct comp id
  integer                      :: mpicom                              ! mpi communicator
  integer                      :: my_task                             ! my task in mpi communicator mpicom
  logical                      :: masterproc                          ! true of my_task == master_task
  character(len=16)            :: inst_suffix = ""                    ! char string associated with instance (ie. "_0001" or "")
  integer                      :: logunit                             ! logging unit number
  logical                      :: restart_read                        ! start from restart
  character(CL)                :: case_name     ! case name
  character(*) , parameter     :: nullstr = 'undefined'

  ! dice_in namelist input
  character(CL)                :: xmlfilename = nullstr               ! filename to obtain namelist info from
  character(CL)                :: nlfilename = nullstr                ! filename to obtain namelist info from
  character(CL)                :: dataMode                            ! flags physics options wrt input data
  character(CL)                :: model_meshfile = nullstr            ! full pathname to model meshfile
  character(CL)                :: model_maskfile = nullstr            ! full pathname to obtain mask from
  character(CL)                :: model_createmesh_fromfile = nullstr ! full pathname to obtain mask from
  real(R8)                     :: flux_swpf                           ! short-wave penatration factor
  real(R8)                     :: flux_Qmin                           ! bound on melt rate
  logical                      :: flux_Qacc                           ! activates water accumulation/melt wrt Q
  real(R8)                     :: flux_Qacc0                          ! initial water accumulation value
  character(CL)                :: restfilm = nullstr                  ! model restart file namelist
  integer                      :: nx_global
  integer                      :: ny_global

  ! linked lists
  type(fldList_type) , pointer :: fldsImport => null()
  type(fldList_type) , pointer :: fldsExport => null()
  type(dfield_type)  , pointer :: dfields    => null()

  ! constants
  logical                      :: flds_i2o_per_cat                    ! .true. if select per ice thickness
  real(R8)                     :: dt                                  ! real model timestep

  integer      , parameter     :: master_task=0                       ! task number of master task
  character(*) , parameter     :: modName =  "(ice_comp_nuopc)"
  character(*) , parameter     :: u_FILE_u = &
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
    integer           :: nu                 ! unit number
    integer           :: ierr               ! error code
    logical           :: exists             ! check for file existence
    character(len=*),parameter  :: subname=trim(modName)//':(InitializeAdvertise) '
    character(*)    ,parameter :: F00 = "('(ice_comp_nuopc) ',8a)"
    character(*)    ,parameter :: F01 = "('(ice_comp_nuopc) ',a,2x,i8)"
    character(*)    ,parameter :: F02 = "('(ice_comp_nuopc) ',a,l6)"
    character(*)    ,parameter :: F03 = "('(ice_comp_nuopc) ',a,d13.5)"
    !-------------------------------------------------------------------------------

    namelist / dice_nml / datamode, model_meshfile, model_maskfile, model_createmesh_fromfile, &
         restfilm, nx_global, ny_global, flux_swpf, flux_Qmin, flux_Qacc, flux_Qacc0

    rc = ESMF_SUCCESS

    call NUOPC_CompAttributeGet(gcomp, name='case_name', value=case_name, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Obtain flds_scalar values, mpi values, multi-instance values and
    ! set logunit and set shr logging to my log file
    call dshr_init(gcomp, mpicom, my_task, inst_index, inst_suffix, &
         flds_scalar_name, flds_scalar_num, flds_scalar_index_nx, flds_scalar_index_ny, &
         logunit, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! determine logical masterproc
    masterproc = (my_task == master_task)

    ! Read dice_nml from nlfilename
    if (my_task == master_task) then
       nlfilename = "dice_in"//trim(inst_suffix)
       open (newunit=nu,file=trim(nlfilename),status="old",action="read")
       read (nu,nml=dice_nml,iostat=ierr)
       close(nu)
       if (ierr > 0) then
          write(logunit,*) 'ERROR: reading input namelist, '//trim(nlfilename)//' iostat=',ierr
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
       write(logunit,F01)' nx_global  = ',nx_global
       write(logunit,F01)' ny_global  = ',ny_global
       write(logunit,F03)' flux_swpf  = ',flux_swpf
       write(logunit,F03)' flux_Qmin  = ',flux_Qmin
       write(logunit,F03)' flux_Qacc  = ',flux_Qacc
       write(logunit,F03)' flux_Qacc0 = ',flux_Qacc0
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

    ! broadcast namelist input
    call shr_mpi_bcast(datamode                  , mpicom, 'datamode')
    call shr_mpi_bcast(model_maskfile            , mpicom, 'model_maskfile')
    call shr_mpi_bcast(model_meshfile            , mpicom, 'model_meshfile')
    call shr_mpi_bcast(model_maskfile            , mpicom, 'model_maskfile')
    call shr_mpi_bcast(model_createmesh_fromfile , mpicom, 'model_createmesh_fromfile')
    call shr_mpi_bcast(nx_global                 , mpicom, 'nx_global')
    call shr_mpi_bcast(ny_global                 , mpicom, 'ny_global')
    call shr_mpi_bcast(restfilm                  , mpicom, 'restfilm')
    call shr_mpi_bcast(flux_swpf                 , mpicom, 'flux_swpf')
    call shr_mpi_bcast(flux_Qmin                 , mpicom, 'flux_Qmin')
    call shr_mpi_bcast(flux_Qacc                 , mpicom, 'flux_Qacc')
    call shr_mpi_bcast(flux_Qacc0                , mpicom, 'flux_Qacc0')

    ! Validate datamode
    if ( trim(datamode) == 'ssmi' .or. trim(datamode) == 'ssmi_iaf') then
       if (my_task == master_task) write(logunit,*) ' dice datamode = ',trim(datamode)
    else
       call shr_sys_abort(' ERROR illegal dice datamode = '//trim(datamode))
    endif

    ! Advertise import and export fields
    call NUOPC_CompAttributeGet(gcomp, name='flds_i2o_per_cat', value=cvalue, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flds_i2o_per_cat  ! module variable

    select case (trim(datamode))
    case('ssmi', 'ssmi_iaf')
       call dice_datamode_ssmi_advertise(importState, exportState, fldsimport, fldsexport, &
            flds_scalar_name, flds_i2o_per_cat, rc)
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
    type(ESMF_TimeInterval)     :: TimeStep
    type(ESMF_Time)             :: currTime
    integer                     :: current_ymd   ! model date
    integer                     :: current_year  ! model year
    integer                     :: current_mon   ! model month
    integer                     :: current_day   ! model day
    integer                     :: current_tod   ! model sec into model date
    real(R8)                    :: cosarg        ! for setting ice temp pattern
    real(R8)                    :: jday, jday0   ! elapsed day counters
    character(CL)               :: cvalue        ! temporary
    integer                     :: n,k           ! generic counters
    integer                     :: model_dt      ! integer model timestep
    character(len=*), parameter :: F00   = "('ice_comp_nuopc: ')',8a)"
    character(len=*), parameter :: subname=trim(modName)//':(InitializeRealize) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Initialize mesh, restart flag, compid, and logunit
    call t_startf('dice_strdata_init')
    call dshr_mesh_init(gcomp, compid, logunit, 'ice', nx_global, ny_global, &
         model_meshfile, model_maskfile, model_createmesh_fromfile, model_mesh, restart_read, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Initialize stream data type
    xmlfilename = 'dice.streams'//trim(inst_suffix)//'.xml'
    call shr_strdata_init_from_xml(sdat, xmlfilename, model_mesh, clock, compid, logunit, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('dice_strdata_init')

    ! NUOPC_Realize "realizes" a previously advertised field in the importState and exportState
    ! by replacing the advertised fields with the newly created fields of the same name.
    call dshr_fldlist_realize( exportState, fldsExport, flds_scalar_name, flds_scalar_num, model_mesh, &
         subname//':diceExport', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_fldlist_realize( importState, fldsImport, flds_scalar_name, flds_scalar_num, model_mesh, &
         subname//':diceImport', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the time to interpolate the stream data to
    call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(currTime, yy=current_year, mm=current_mon, dd=current_day, s=current_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(current_year, current_mon, current_day, current_ymd)

    ! Get model timestep
    call ESMF_TimeIntervalGet( timeStep, s=model_dt, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    dt = model_dt * 1.0_r8

    ! Get cosarg
    call shr_cal_ymd2julian(0, current_mon, current_day, current_tod, jDay , sdat%model_calendar) ! julian day for model
    call shr_cal_ymd2julian(0, 9,           1,           0,           jDay0, sdat%model_calendar) ! julian day for Sept 1
    cosArg = 2.0_R8*shr_const_pi*(jday - jday0)/365.0_R8

    ! Run dice
    call dice_comp_run(importState, exportState, current_ymd, current_tod, cosarg, restart_write=.false., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Add scalars to export state
    call dshr_state_SetScalar(dble(nx_global), flds_scalar_index_nx, exportState, flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_SetScalar(dble(ny_global), flds_scalar_index_ny, exportState, flds_scalar_name, flds_scalar_num, rc)
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
    type(ESMF_Alarm)        :: alarm
    type(ESMF_TimeInterval) :: timeStep
    type(ESMF_Time)         :: currTime, nextTime
    integer                 :: current_mon   ! model month
    integer                 :: current_day   ! model day
    integer                 :: current_tod   ! model sec into model date
    real(R8)                :: cosarg        ! for setting ice temp pattern
    real(R8)                :: jday, jday0   ! elapsed day counters
    integer                 :: next_ymd      ! model date
    integer                 :: next_tod      ! model sec into model date
    integer                 :: yr            ! year
    integer                 :: mon           ! month
    integer                 :: day           ! day in month
    logical                 :: restart_write
    character(len=*),parameter :: subname=trim(modName)//':(ModelAdvance) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call t_startf(subname)
    call memcheck(subname, 5, my_task == master_task)

    ! Query the Component for its clock, importState and exportState
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

    ! Get cosarg
    call shr_cal_ymd2julian(0, mon, day, next_tod, jDay , sdat%model_calendar) ! julian day for model
    call shr_cal_ymd2julian(0, 9,   1,   0,        jDay0, sdat%model_calendar) ! julian day for Sept 1
    cosArg = 2.0_R8*shr_const_pi*(jday - jday0)/365.0_R8

    ! Determine if will write restarts
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

    ! Run dice
    call dice_comp_run(importState, exportState, next_ymd, next_tod, cosarg, restart_write, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call t_stopf(subname)

  end subroutine ModelAdvance

  !===============================================================================
  subroutine dice_comp_run(importstate, exportstate, target_ymd, target_tod, cosarg, restart_write, rc)

    ! --------------------------
    ! advance dice
    ! --------------------------

    ! input/output variables:
    type(ESMF_State) , intent(inout) :: exportState
    type(ESMF_State) , intent(inout) :: importState
    integer          , intent(in)    :: target_ymd ! model date
    integer          , intent(in)    :: target_tod ! model sec into model date
    real(R8)         , intent(in)    :: cosarg     ! for setting ice temp pattern
    logical          , intent(in)    :: restart_write
    integer          , intent(out)   :: rc

    ! local variables
    logical :: first_time = .true.
    character(*), parameter :: subName = "(dice_comp_run) "
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('DICE_RUN')

    !--------------------
    ! first time initialization
    !--------------------

    if (first_time) then

       ! Initialize dfields with export state data that has corresponding stream field
       call dshr_dfield_add(dfields, sdat, state_fld='Si_ifrac', strm_fld='Si_ifrac', &
            state=exportState, logunit=logunit, masterproc=masterproc, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Initialize datamode module ponters
       select case (trim(datamode))
       case('ssmi', 'ssmi_iaf')
          call dice_datamode_ssmi_init_pointers(importState, exportState, sdat, flds_i2o_per_cat, rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end select

       ! read restart if needed
       if (restart_read) then
          select case (trim(datamode))
          case('ssmi', 'ssmi_iaf')
             call dice_datamode_ssmi_restart_read(restfilm, inst_suffix, logunit, my_task, mpicom, sdat)
          end select
       end if

       ! reset first_time
       first_time = .false.
    end if

    !--------------------
    ! advance dice streams
    !--------------------

    ! time and spatially interpolate to model time and grid
    call t_barrierf('dice_BARRIER',mpicom)
    call t_startf('dice_strdata_advance')
    call shr_strdata_advance(sdat, target_ymd, target_tod, logunit, 'dice', rc=rc)
    call t_stopf('dice_strdata_advance')

    !--------------------
    ! copy all fields from streams to export state as default
    !--------------------

    ! This automatically will update the fields in the export state
    call t_barrierf('dice_comp_dfield_copy_BARRIER', mpicom)
    call t_startf('dice_dfield_copy')
    call dshr_dfield_copy(dfields, sdat, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('dice_dfield_copy')

    !-------------------------------------------------
    ! Determine data model behavior based on the mode
    !-------------------------------------------------

    call t_startf('dice_datamode')

    ! Perform data mode specific calculations
    select case (trim(datamode))
    case ('ssmi', 'ssmi_iaf')
       call dice_datamode_ssmi_advance(exportState, importState, cosarg, flds_i2o_per_cat, &
            flux_swpf, flux_Qmin, flux_Qacc, flux_Qacc0, dt, logunit, restart_read, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end select

    ! Write restarts if needed
    if (restart_write) then
       select case (trim(datamode))
       case('ssmi', 'ssmi_iaf')
          call dice_datamode_ssmi_restart_write(case_name, inst_suffix, target_ymd, target_tod, &
               logunit, mpicom, my_task, sdat)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end select
    end if

    ! Write diagnostics
    call dshr_state_diagnose(exportState,subname//':ES',rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call t_stopf('dice_datamode')
    call t_stopf('DICE_RUN')

  end subroutine dice_comp_run

  !===============================================================================
  subroutine ModelFinalize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (my_task == master_task) then
       write(logunit,*)
       write(logunit,*) 'dice : end of main integration loop'
       write(logunit,*)
    end if
  end subroutine ModelFinalize

end module ice_comp_nuopc
