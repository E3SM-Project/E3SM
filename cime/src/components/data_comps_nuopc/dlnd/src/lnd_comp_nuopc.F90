module lnd_comp_nuopc

  !----------------------------------------------------------------------------
  ! This is the NUOPC cap for DLND
  !----------------------------------------------------------------------------

  use ESMF
  use NUOPC             , only : NUOPC_CompDerive, NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize
  use NUOPC             , only : NUOPC_CompAttributeGet, NUOPC_Advertise
  use NUOPC_Model       , only : model_routine_SS        => SetServices
  use NUOPC_Model       , only : model_label_Advance     => label_Advance
  use NUOPC_Model       , only : model_label_SetRunClock => label_SetRunClock
  use NUOPC_Model       , only : model_label_Finalize    => label_Finalize
  use NUOPC_Model       , only : NUOPC_ModelGet
  use shr_file_mod      , only : shr_file_getlogunit, shr_file_setlogunit
  use shr_kind_mod      , only : r8=>shr_kind_r8, i8=>shr_kind_i8, cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_sys_mod       , only : shr_sys_abort
  use shr_cal_mod       , only : shr_cal_ymd2date
  use shr_mpi_mod       , only : shr_mpi_bcast
  use dshr_methods_mod  , only : dshr_state_getfldptr, dshr_state_diagnose, chkerr, memcheck
  use dshr_strdata_mod  , only : shr_strdata_type, shr_strdata_advance, shr_strdata_get_stream_domain
  use dshr_mod          , only : dshr_model_initphase, dshr_init, dshr_sdat_init 
  use dshr_mod          , only : dshr_state_setscalar
  use dshr_mod          , only : dshr_set_runclock, dshr_log_clock_advance
  use dshr_mod          , only : dshr_restart_read, dshr_restart_write
  use dshr_mod          , only : dshr_create_mesh_from_grid
  use dshr_dfield_mod   , only : dfield_type, dshr_dfield_add, dshr_dfield_copy
  use dshr_fldlist_mod  , only : fldlist_type, dshr_fldlist_add, dshr_fldlist_realize
  use glc_elevclass_mod , only : glc_elevclass_as_string, glc_elevclass_init
  use perf_mod          , only : t_startf, t_stopf, t_adj_detailf, t_barrierf

  implicit none
  private ! except

  public  :: SetServices

  private :: InitializeAdvertise
  private :: InitializeRealize
  private :: ModelAdvance
  private :: ModelFinalize
  private :: dlnd_comp_advertise
  private :: dlnd_comp_realize
  private :: dlnd_comp_run

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  type(shr_strdata_type)       :: sdat
  type(ESMF_Mesh)              :: mesh
  character(len=CS)            :: flds_scalar_name = ''
  integer                      :: flds_scalar_num = 0
  integer                      :: flds_scalar_index_nx = 0
  integer                      :: flds_scalar_index_ny = 0
  integer                      :: compid                          ! mct comp id
  integer                      :: mpicom                          ! mpi communicator
  integer                      :: my_task                         ! my task in mpi communicator mpicom
  logical                      :: masterproc                      ! true of my_task == master_task
  integer                      :: inst_index                      ! number of current instance (ie. 1)
  character(len=16)            :: inst_name                       ! fullname of current instance (ie. "lnd_0001")
  character(len=16)            :: inst_suffix = ""                ! char string associated with instance (ie. "_0001" or "")
  integer                      :: logunit                         ! logging unit number
  logical                      :: read_restart                    ! start from restart
  character(*) , parameter     :: nullstr = 'undefined'

  ! dlnd_in namelist input
  character(CL)                :: nlfilename                      ! filename to obtain namelist info from
  character(CL)                :: dataMode                        ! flags physics options wrt input data
  character(CL)                :: domain_fracname = 'undefined'   ! name of fraction field on first stream file
  logical                      :: force_prognostic_true = .false. ! if true set prognostic true
  character(CL)                :: restfilm = nullstr              ! model restart file namelist
  character(CL)                :: restfils = nullstr              ! stream restart file namelist

  ! linked lists
  type(fldList_type) , pointer :: fldsImport => null()
  type(fldList_type) , pointer :: fldsExport => null()
  type(dfield_type)  , pointer :: dfields    => null()

  ! module pointer arrays
  real(r8), pointer            :: lfrac(:)

  ! module constants
  integer                      :: glc_nec
  integer      , parameter     :: master_task=0                   ! task number of master task
  character(*) , parameter     :: rpfile = 'rpointer.lnd'
  character(*) , parameter     :: modName =  "(lnd_comp_nuopc)"

  character(*) , parameter     :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

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
    character(CL) :: cvalue
    integer       :: shrlogunit                      ! original log unit
    integer       :: nu                              ! unit number
    integer       :: ierr                            ! error code
    character(len=*), parameter :: subname=trim(modName)//':(InitializeAdvertise) '
    !-------------------------------------------------------------------------------

    namelist / dlnd_nml / datamode, restfilm, restfils, force_prognostic_true, domain_fracname

    rc = ESMF_SUCCESS

    ! Obtain flds_scalar values, mpi values, multi-instance values and  
    ! set logunit and set shr logging to my log file
    call dshr_init(gcomp, mpicom, my_task, inst_index, inst_suffix, &
         flds_scalar_name, flds_scalar_num, flds_scalar_index_nx, flds_scalar_index_ny, &
         logunit, shrlogunit, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine namelist filename
    nlfilename = "dlnd_in"//trim(inst_suffix)

    ! Determine logical masterproc
    masterproc = (my_task == master_task)

    ! Read dlnd_nml from nlfilename
    if (my_task == master_task) then
       open (newunit=nu, file=trim(nlfilename), status="old", action="read")
       read (nu,nml=dlnd_nml,iostat=ierr)
       close(nu)
       if (ierr > 0) then
          write(logunit,*) 'ERROR: reading input namelist, '//trim(nlfilename)//' iostat=',ierr
          call shr_sys_abort(subName//': namelist read error '//trim(nlfilename))
       end if
       write(logunit,*)' datamode   = ',datamode
       write(logunit,*)' restfilm   = ',trim(restfilm)
       write(logunit,*)' restfils   = ',trim(restfils)
       write(logunit,*)' force_prognostic_true = ',force_prognostic_true
       write(logunit,*)' domain_fracname = ',trim(domain_fracname)
    endif
    call shr_mpi_bcast(datamode ,mpicom, 'datamode')
    call shr_mpi_bcast(domain_fracname,mpicom,'domain_fracname')
    call shr_mpi_bcast(force_prognostic_true,mpicom,'force_prognostic_true')
    call shr_mpi_bcast(restfilm ,mpicom, 'restfilm')
    call shr_mpi_bcast(restfils ,mpicom, 'restfils')

    ! Validate sdat datamode
    if (trim(datamode) == 'NULL' .or. trim(datamode) == 'COPYALL') then
       if (my_task == master_task) write(logunit,*) 'dlnd datamode = ',trim(datamode)
    else
       call shr_sys_abort(' ERROR illegal dlnd datamode = '//trim(datamode))
    end if

    ! Advertise the export fields
    if (trim(datamode) /= 'NULL') then
       call NUOPC_CompAttributeGet(gcomp, name='glc_nec', value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) glc_nec
       call ESMF_LogWrite('glc_nec = '// trim(cvalue), ESMF_LOGMSG_INFO)

       call dlnd_comp_advertise(importState, exportState, rc=rc)
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
    type(ESMF_TIME) :: currTime
    integer         :: current_ymd  ! model date
    integer         :: current_year ! model year
    integer         :: current_mon  ! model month
    integer         :: current_day  ! model day
    integer         :: current_tod  ! model sec into model date
    character(CL)   :: cvalue       ! temporary
    integer         :: shrlogunit   ! original log unit
    character(len=*),parameter :: subname=trim(modName)//':(InitializeRealize) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (datamode == 'NULL') RETURN

    ! Reset shr logging to my log file
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_setLogUnit (logUnit)

    ! Initialize sdat
    call t_startf('dlnd_strdata_init')
    call dshr_sdat_init(gcomp, clock, nlfilename, compid, logunit, 'lnd', mesh, read_restart, sdat, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('dlnd_strdata_init')

    ! Realize the actively coupled fields, now that a mesh is established and
    ! initialize dfields data type (to map streams to export state fields)
    call dlnd_comp_realize(importState, exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Read restart if necessary
    if (read_restart) then
       call dshr_restart_read(restfilm, restfils, rpfile, inst_suffix, nullstr, &
            logunit, my_task, mpicom, sdat)
    end if

    ! get the time to interpolate the stream data to
    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(currTime, yy=current_year, mm=current_mon, dd=current_day, s=current_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(current_year, current_mon, current_day, current_ymd)

    ! Run dlnd to create export state
    call dlnd_comp_run(current_ymd, current_tod, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! add scalars to export state
    call dshr_state_SetScalar(dble(sdat%nxg),flds_scalar_index_nx, exportState, flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_SetScalar(dble(sdat%nyg),flds_scalar_index_ny, exportState, flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! diagnostics
    call dshr_state_diagnose(exportState,subname//':ES',rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Reset shr logging to original values
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
    character(len=*),parameter :: subname=trim(modName)//':(ModelAdvance) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call memcheck(subname, 5, my_task==master_task)

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

    ! run dlnd
    call dlnd_comp_run(next_ymd, next_tod, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! write_restart if alarm is ringing
    call ESMF_ClockGetAlarm(clock, alarmname='alarm_restart', alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call t_startf('dlnd_restart')
       call NUOPC_CompAttributeGet(gcomp, name='case_name', value=case_name, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call dshr_restart_write(rpfile, case_name, 'dlnd', inst_suffix, next_ymd, next_tod, &
            logunit, mpicom, my_task, sdat)
       call t_stopf('dlnd_restart')
    endif

    ! write diagnostics
    call dshr_state_diagnose(exportState,subname//':ES',rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (my_task == master_task) then
       call dshr_log_clock_advance(clock, 'DLND', logunit, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif
    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

    call shr_file_setLogUnit (shrlogunit)

  end subroutine ModelAdvance

  !===============================================================================
  subroutine ModelFinalize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(*), parameter :: F00   = "('(dlnd_comp_final) ',8a)"
    character(*), parameter :: F91   = "('(dlnd_comp_final) ',73('-'))"
    character(len=*),parameter  :: subname=trim(modName)//':(ModelFinalize) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
    if (my_task == master_task) then
       write(logunit,F91)
       write(logunit,F00) ' dlnd : end of main integration loop'
       write(logunit,F91)
    end if
    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO)

  end subroutine ModelFinalize

  !===============================================================================
  subroutine dlnd_comp_advertise(importState, exportState, rc)

    ! determine export and import fields to advertise to mediator

    ! input/output arguments
    type(ESMF_State)               :: importState
    type(ESMF_State)               :: exportState
    integer          , intent(out) :: rc

    ! local variables
    type(fldlist_type), pointer :: fldList
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call glc_elevclass_init(glc_nec)

    !-------------------
    ! Advertise export fields
    !-------------------

    call dshr_fldList_add(fldsExport, trim(flds_scalar_name))
    call dshr_fldlist_add(fldsExport, "Sl_lfrin")

    ! The following puts all of the elevation class fields as an
    ! undidstributed dimension in the export state field - index1 is bare land - and the total number of
    ! elevation classes not equal to bare land go from index2 -> glc_nec+1
    if (glc_nec > 0) then
       call dshr_fldList_add(fldsExport, 'Sl_tsrf_elev'  , ungridded_lbound=1, ungridded_ubound=glc_nec+1)
       call dshr_fldList_add(fldsExport, 'Sl_topo_elev'  , ungridded_lbound=1, ungridded_ubound=glc_nec+1)
       call dshr_fldList_add(fldsExport, 'Flgl_qice_elev', ungridded_lbound=1, ungridded_ubound=glc_nec+1)
    end if

    fldlist => fldsExport ! the head of the linked list
    do while (associated(fldlist))
       call NUOPC_Advertise(exportState, standardName=fldlist%stdname, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite('(dlnd_comp_advertise): Fr_lnd '//trim(fldList%stdname), ESMF_LOGMSG_INFO)
       fldList => fldList%next
    enddo

    ! TODO: Non snow fields that nead to be added if dlnd is in cplhist mode
    ! "Sl_t        " "Sl_tref     " "Sl_qref     " "Sl_avsdr    "
    ! "Sl_anidr    " "Sl_avsdf    " "Sl_anidf    " "Sl_snowh    "
    ! "Fall_taux   " "Fall_tauy   " "Fall_lat    " "Fall_sen    "
    ! "Fall_lwup   " "Fall_evap   " "Fall_swnet  " "Sl_landfrac "
    ! "Sl_fv       " "Sl_ram1     "
    ! "Fall_flxdst1" "Fall_flxdst2" "Fall_flxdst3" "Fall_flxdst4"

  end subroutine dlnd_comp_advertise

  !===============================================================================
  subroutine dlnd_comp_realize(importState, exportState, rc)

    ! input/output variables
    type(ESMF_State) , intent(inout) :: importState
    type(ESMF_State) , intent(inout) :: exportState
    integer          , intent(out)   :: rc

    ! local variables
    integer                    :: n
    character(len=2)           :: nec_str
    character(CS), allocatable :: strm_flds(:)
    character(*), parameter    :: subName = "(dlnd_comp_realize) "
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    ! -------------------------------------
    ! NUOPC_Realize "realizes" a previously advertised field in the importState and exportState
    ! by replacing the advertised fields with the newly created fields of the same name.
    ! -------------------------------------

    call dshr_fldlist_realize( exportState, fldsExport, flds_scalar_name, flds_scalar_num,  mesh, &
         subname//':dlndExport', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Set fractional land pointer in export state
    call dshr_state_getfldptr(exportState, fldname='Sl_lfrin', fldptr1=lfrac, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Obtain fractional land from first stream
    call shr_strdata_get_stream_domain(sdat, 1, mpicom, my_task, domain_fracname, lfrac) 

    ! Create stream-> export state mapping
    allocate(strm_flds(0:glc_nec))
    do n = 0,glc_nec
       nec_str = glc_elevclass_as_string(n)
       strm_flds(n) = 'tsrf' // trim(nec_str)
    end do
    call dshr_dfield_add(dfields, sdat, state_fld='Sl_tsrf_elev', strm_flds=strm_flds, state=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do n = 0,glc_nec
       nec_str = glc_elevclass_as_string(n)
       strm_flds(n) = 'topo' // trim(nec_str)
    end do
    call dshr_dfield_add(dfields, sdat, state_fld='Sl_topo_elev', strm_flds=strm_flds, state=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do n = 0,glc_nec
       nec_str = glc_elevclass_as_string(n)
       strm_flds(n) = 'qice' // trim(nec_str)
    end do
    call dshr_dfield_add(dfields, sdat, state_fld='Flgl_qice_elev', strm_flds=strm_flds, state=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine dlnd_comp_realize

  !===============================================================================
  subroutine dlnd_comp_run(target_ymd, target_tod, rc)

    ! --------------------------
    ! advance dlnd
    ! --------------------------

    ! input/output variables:
    integer , intent(in)    :: target_ymd       ! model date
    integer , intent(in)    :: target_tod       ! model sec into model date
    integer , intent(out)   :: rc
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('DLND_RUN')

    !--------------------
    ! advance dlnd streams
    !--------------------

    ! time and spatially interpolate to model time and grid
    call t_barrierf('dlnd_BARRIER',mpicom)
    call t_startf('dlnd_strdata_advance')
    call shr_strdata_advance(sdat, target_ymd, target_tod, mpicom, 'dlnd')
    call t_stopf('dlnd_strdata_advance')

    !--------------------
    ! copy all fields from streams to export state as default
    !--------------------

    ! This automatically will update the fields in the export state
    call t_barrierf('dlnd_comp_strdata_copy_BARRIER', mpicom)
    call t_startf('dlnd_strdata_copy')
    call dshr_dfield_copy(dfields, sdat, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('dlnd_strdata_copy')

    !-------------------------------------------------
    ! determine data model behavior based on the mode
    !-------------------------------------------------

    call t_startf('dlnd_datamode')
    select case (trim(datamode))
    case('COPYALL')
       ! do nothing extra
    end select
    call t_stopf('dlnd_datamode')

    call t_stopf('DLND_RUN')

  end subroutine dlnd_comp_run

end module lnd_comp_nuopc
