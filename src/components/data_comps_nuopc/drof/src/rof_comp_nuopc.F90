module rof_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for DROF
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
  use shr_const_mod    , only : SHR_CONST_SPVAL
  use shr_sys_mod      , only : shr_sys_abort
  use shr_cal_mod      , only : shr_cal_ymd2date
  use shr_mpi_mod      , only : shr_mpi_bcast
  use dshr_methods_mod , only : dshr_state_getfldptr, dshr_state_diagnose, chkerr, memcheck
  use dshr_strdata_mod , only : shr_strdata_type, shr_strdata_advance, shr_strdata_get_stream_domain
  use dshr_strdata_mod , only : shr_strdata_init_from_xml
  use dshr_mod         , only : dshr_model_initphase, dshr_init
  use dshr_mod         , only : dshr_state_setscalar, dshr_set_runclock, dshr_log_clock_advance
  use dshr_mod         , only : dshr_restart_read, dshr_restart_write, dshr_mesh_init
  use dshr_dfield_mod  , only : dfield_type, dshr_dfield_add, dshr_dfield_copy
  use dshr_fldlist_mod , only : fldlist_type, dshr_fldlist_add, dshr_fldlist_realize
  use perf_mod         , only : t_startf, t_stopf, t_adj_detailf, t_barrierf

  implicit none
  private ! except

  public  :: SetServices

  private :: InitializeAdvertise
  private :: InitializeRealize
  private :: ModelAdvance
  private :: ModelFinalize
  private :: drof_comp_advertise
  private :: drof_comp_realize
  private :: drof_comp_run

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  type(shr_strdata_type)       :: sdat
  type(ESMF_Mesh)              :: model_mesh                          ! model mesh
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
  logical                      :: read_restart
  character(*) , parameter     :: nullstr = 'undefined'
                                                                      ! drof_in namelist input
  character(CL)                :: xmlfilename = nullstr               ! filename to obtain namelist info from
  character(CL)                :: nlfilename = nullstr                ! filename to obtain namelist info from
  character(CL)                :: dataMode = nullstr                  ! flags physics options wrt input data
  character(CL)                :: model_meshfile = nullstr            ! full pathname to model meshfile
  character(CL)                :: model_maskfile = nullstr            ! full pathname to obtain mask from
  character(CL)                :: model_createmesh_fromfile = nullstr ! full pathname to obtain mask from
  logical                      :: force_prognostic_true = .false.     ! if true set prognostic true
  character(CL)                :: restfilm = nullstr                  ! model restart file namelist
  character(CL)                :: restfils = nullstr                  ! stream restart file namelist
  integer                      :: nx_global
  integer                      :: ny_global

                                                                      ! constants
  integer      , parameter     :: master_task=0                       ! task number of master task
  character(*) , parameter     :: rpfile = 'rpointer.rof'
  character(*) , parameter     :: modName =  "(rof_comp_nuopc)"

  ! linked lists
  type(fldList_type) , pointer :: fldsImport => null()
  type(fldList_type) , pointer :: fldsExport => null()
  type(dfield_type)  , pointer :: dfields    => null()

  ! module pointer arrays
  real(r8), pointer            :: Forr_rofl(:) => null()
  real(r8), pointer            :: Forr_rofi(:) => null()

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
    integer           :: inst_index ! number of current instance (ie. 1)
    character(len=CL) :: cvalue     ! temporary
    integer           :: nu         ! unit number
    integer           :: ierr       ! error code
    logical           :: exists     ! check for file existence
    character(len=*),parameter :: subname=trim(modName)//':(InitializeAdvertise) '
    character(*)    ,parameter :: F00 = "('(rof_comp_nuopc) ',8a)"
    character(*)    ,parameter :: F01 = "('(rof_comp_nuopc) ',a,2x,i8)"
    character(*)    ,parameter :: F02 = "('(rof_comp_nuopc) ',a,l6)"
    !-------------------------------------------------------------------------------

    namelist / drof_nml / datamode, model_meshfile, model_maskfile, model_createmesh_fromfile, &
         restfilm, restfils, force_prognostic_true, nx_global, ny_global

    rc = ESMF_SUCCESS

    ! Obtain flds_scalar values, mpi values, multi-instance values and
    ! set logunit and set shr logging to my log file
    call dshr_init(gcomp, mpicom, my_task, inst_index, inst_suffix, &
         flds_scalar_name, flds_scalar_num, flds_scalar_index_nx, flds_scalar_index_ny, &
         logunit, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! determine logical masterproc
    masterproc = (my_task == master_task)

    ! Read drof_nml from nlfilename
    if (masterproc) then
       nlfilename = "drof_in"//trim(inst_suffix)
       open (newunit=nu,file=trim(nlfilename),status="old",action="read")
       read (nu,nml=drof_nml,iostat=ierr)
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
       write(logunit,F01)' nx_global = ',nx_global
       write(logunit,F01)' ny_global = ',ny_global
       write(logunit,F00)' restfilm = ',trim(restfilm)
       write(logunit,F00)' restfils = ',trim(restfils)
       write(logunit,F02)' force_prognostic_true = ',force_prognostic_true

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
    call shr_mpi_bcast(model_meshfile            , mpicom, 'model_meshfile')
    call shr_mpi_bcast(model_maskfile            , mpicom, 'model_maskfile')
    call shr_mpi_bcast(model_createmesh_fromfile , mpicom, 'model_createmesh_fromfile')
    call shr_mpi_bcast(nx_global                 , mpicom, 'nx_global')
    call shr_mpi_bcast(ny_global                 , mpicom, 'ny_global')
    call shr_mpi_bcast(restfilm                  , mpicom, 'restfilm')
    call shr_mpi_bcast(restfils                  , mpicom, 'restfils')
    call shr_mpi_bcast(force_prognostic_true     , mpicom, 'force_prognostic_true')

    ! Validate datamode
    if (trim(datamode) == 'NULL' .or. trim(datamode) == 'COPYALL') then
       if (masterproc) write(logunit,*) 'drof datamode = ',trim(datamode)
    else
       call shr_sys_abort(' ERROR illegal drof datamode = '//trim(datamode))
    end if
    if (trim(datamode) /= 'NULL') then
       call drof_comp_advertise(importState, exportState, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
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
    type(ESMF_TIME) :: currTime
    integer         :: current_ymd  ! model date
    integer         :: current_year ! model year
    integer         :: current_mon  ! model month
    integer         :: current_day  ! model day
    integer         :: current_tod  ! model sec into model date
    character(CL)   :: cvalue       ! temporary
    integer         :: n,k          ! generic counters
    character(len=*), parameter :: F00   = "('rof_comp_nuopc: ')',8a)"
    character(len=*), parameter :: subname=trim(modName)//':(InitializeRealize) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (datamode == 'NULL') RETURN

    ! Initialize mesh, restart flag, compid, and logunit
    call t_startf('drof_strdata_init')
    call dshr_mesh_init(gcomp, compid, logunit, 'rof', nx_global, ny_global, &
         model_meshfile, model_maskfile, model_createmesh_fromfile,  model_mesh, read_restart, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Initialize stream data type
    xmlfilename = 'drof.streams'//trim(inst_suffix)//'.xml'
    call shr_strdata_init_from_xml(sdat, xmlfilename, model_mesh, clock, compid, logunit, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('drof_strdata_init')

    ! Realize the actively coupled fields, now that a mesh is established and
    ! initialize dfields data type (to map streams to export state fields)
    call drof_comp_realize(importState, exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Read restart if necessary
    if (read_restart) then
       call dshr_restart_read(restfilm, restfils, rpfile, inst_suffix, nullstr, logunit, my_task, mpicom, sdat)
    end if

    ! Get the time to interpolate the stream data to
    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(currTime, yy=current_year, mm=current_mon, dd=current_day, s=current_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(current_year, current_mon, current_day, current_ymd)

    ! Run drof
    call drof_comp_run(current_ymd, current_tod, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Add scalars to export state
    call dshr_state_SetScalar(dble(nx_global),flds_scalar_index_nx, exportState, flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_SetScalar(dble(ny_global),flds_scalar_index_ny, exportState, flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Diagnostics
    call dshr_state_diagnose(exportState, subname//':ES',rc=rc)
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
    integer                 :: next_ymd      ! model date
    integer                 :: next_tod      ! model sec into model date
    integer                 :: yr            ! year
    integer                 :: mon           ! month
    integer                 :: day           ! day in month
    character(CL)           :: case_name     ! case name
    character(len=*),parameter :: subname=trim(modName)//':(ModelAdvance) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call memcheck(subname, 5, masterproc)

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

    ! run drof
    call drof_comp_run(next_ymd, next_tod, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! write_restart if alarm is ringing
    call ESMF_ClockGetAlarm(clock, alarmname='alarm_restart', alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call t_startf('drof_restart')
       call NUOPC_CompAttributeGet(gcomp, name='case_name', value=case_name, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call dshr_restart_write(rpfile, case_name, 'drof', inst_suffix, next_ymd, next_tod, &
            logunit, mpicom, my_task, sdat)
       call t_stopf('drof_restart')
    endif

    ! write diagnostics
    call dshr_state_diagnose(exportState,subname//':ES',rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (masterproc) then
       call dshr_log_clock_advance(clock, 'drof', logunit, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

  end subroutine ModelAdvance

  !===============================================================================

  subroutine ModelFinalize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    if (masterproc) then
       write(logunit,*)
       write(logunit,*) 'drof : end of main integration loop'
       write(logunit,*)
    end if

  end subroutine ModelFinalize

  !===============================================================================

  subroutine drof_comp_advertise(importState, exportState, rc)

    ! determine export and import fields to advertise to mediator

    ! input/output arguments
    type(ESMF_State) , intent(inout) :: importState
    type(ESMF_State) , intent(inout) :: exportState
    integer          , intent(out)   :: rc

    ! local variables
    integer :: n
    type(fldlist_type), pointer :: fldList
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !-------------------
    ! Advrtise export fields
    !-------------------

    call dshr_fldList_add(fldsExport, trim(flds_scalar_name))
    call dshr_fldlist_add(fldsExport, "Forr_rofl")
    call dshr_fldlist_add(fldsExport, "Forr_rofi")

    fldlist => fldsExport ! the head of the linked list
    do while (associated(fldlist))
       call NUOPC_Advertise(exportState, standardName=fldlist%stdname, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite('(drof_comp_advertise): Fr_rof '//trim(fldList%stdname), ESMF_LOGMSG_INFO)
       fldList => fldList%next
    enddo

    ! currently there is no import state to drof

  end subroutine drof_comp_advertise

  !===============================================================================

  subroutine drof_comp_realize(importState, exportState, rc)

    !input/output variables
    type(ESMF_State)       , intent(inout) :: importState
    type(ESMF_State)       , intent(inout) :: exportState
    integer                , intent(out)   :: rc

    ! local variables
    character(*), parameter :: subName = "(drof_comp_realize) "
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! -------------------------------------
    ! NUOPC_Realize "realizes" a previously advertised field in the importState and exportState
    ! by replacing the advertised fields with the newly created fields of the same name.
    ! -------------------------------------

    call dshr_fldlist_realize( exportState, fldsExport, flds_scalar_name, flds_scalar_num, model_mesh, &
         subname//':drofExport', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! -------------------------------------
    ! Set pointers to exportState fields
    ! -------------------------------------

    call dshr_dfield_add(dfields, sdat, state_fld='Forr_rofl', strm_fld='rofl', &
         state=exportState, state_ptr=Forr_rofl, logunit=logunit, masterproc=masterproc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add(dfields, sdat, state_fld='Forr_rofi', strm_fld='rofi', &
         state=exportState, state_ptr=Forr_rofi, logunit=logunit, masterproc=masterproc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine drof_comp_realize

  !===============================================================================

  subroutine drof_comp_run(target_ymd, target_tod, rc)

    ! --------------------------
    ! advance drof
    ! --------------------------

    ! input/output variables:
    integer , intent(in)    :: target_ymd       ! model date
    integer , intent(in)    :: target_tod       ! model sec into model date
    integer , intent(out)   :: rc

    ! local variables
    integer :: n
    character(*), parameter :: subName = "(drof_comp_run) "
    !-------------------------------------------------------------------------------

    call t_startf('DROF_RUN')

    !--------------------
    ! advance drof streams
    !--------------------

    ! time and spatially interpolate to model time and grid
    call t_barrierf('drof_BARRIER',mpicom)
    call t_startf('drof_strdata_advance')
    call shr_strdata_advance(sdat, target_ymd, target_tod, logunit, 'drof', rc=rc)
    call t_stopf('drof_strdata_advance')

    !--------------------
    ! copy all fields from streams to export state as default
    !--------------------

    ! This automatically will update the fields in the export state
    call t_barrierf('drof_comp_dfield_copy_BARRIER', mpicom)
    call t_startf('drof_dfield_copy')
    call dshr_dfield_copy(dfields,  sdat, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('drof_dfield_copy')

    !-------------------------------------------------
    ! determine data model behavior based on the mode
    !-------------------------------------------------

    call t_startf('drof_datamode')
    select case (trim(datamode))
    case('COPYALL')
       ! zero out "special values" of export fields
       do n = 1, size(Forr_rofl)
          if (abs(Forr_rofl(n)) > 1.0e28) Forr_rofl(n) = 0.0_r8
          if (abs(Forr_rofi(n)) > 1.0e28) Forr_rofi(n) = 0.0_r8
       enddo
    end select
    call t_stopf('drof_datamode')

    call t_stopf('DROF_RUN')

  end subroutine drof_comp_run

end module rof_comp_nuopc
