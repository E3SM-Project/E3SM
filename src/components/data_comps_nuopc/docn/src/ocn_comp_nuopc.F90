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
  use shr_mpi_mod      , only : shr_mpi_bcast
  use shr_frz_mod      , only : shr_frz_freezetemp
  use shr_const_mod    , only : shr_const_cpsw, shr_const_rhosw, shr_const_TkFrz
  use shr_const_mod    , only : shr_const_TkFrzSw, shr_const_latice, shr_const_ocn_ref_sal
  use shr_const_mod    , only : shr_const_zsrflyr, shr_const_pi
  use dshr_methods_mod , only : dshr_state_getfldptr, dshr_state_diagnose, chkerr, memcheck
  use dshr_strdata_mod , only : shr_strdata_type, shr_strdata_advance, shr_strdata_get_stream_domain
  use dshr_mod         , only : dshr_model_initphase, dshr_init, dshr_sdat_init
  use dshr_mod         , only : dshr_state_setscalar, dshr_set_runclock, dshr_log_clock_advance
  use dshr_mod         , only : dshr_restart_read, dshr_restart_write
  use dshr_mod         , only : dshr_create_mesh_from_grid
  use dshr_strdata_mod , only : shr_strdata_type, shr_strdata_advance
  use dshr_dfield_mod  , only : dfield_type, dshr_dfield_add, dshr_dfield_copy
  use dshr_fldlist_mod , only : fldlist_type, dshr_fldlist_add, dshr_fldlist_realize
  use perf_mod         , only : t_startf, t_stopf, t_adj_detailf, t_barrierf
  use pio

  implicit none
  private ! except

  public  :: SetServices

  private :: InitializeAdvertise
  private :: InitializeRealize
  private :: ModelAdvance
  private :: ModelFinalize
  private :: docn_comp_advertise
  private :: docn_comp_realize
  private :: docn_comp_run
  private :: docn_prescribed_sst

  !--------------------------------------------------------------------------
  ! Private module data
  !--------------------------------------------------------------------------

  type(shr_strdata_type)       :: sdat
  type(ESMF_Mesh)              :: mesh                            ! model mesh
  character(len=CS)            :: flds_scalar_name = ''
  integer                      :: flds_scalar_num = 0
  integer                      :: flds_scalar_index_nx = 0
  integer                      :: flds_scalar_index_ny = 0
  integer                      :: compid                          ! mct comp id
  integer                      :: mpicom                          ! mpi communicator
  integer                      :: my_task                         ! my task in mpi communicator mpicom
  logical                      :: masterproc                      ! true of my_task == master_task
  character(len=16)            :: inst_suffix = ""                ! char string associated with instance (ie. "_0001" or "")
  integer                      :: logunit                         ! logging unit number
  logical                      :: read_restart                    ! start from restart
  character(*) , parameter     :: nullstr = 'undefined'

  ! docn_in namelist input
  character(CL)                :: nlfilename                      ! filename to obtain namelist info from
  character(CL)                :: dataMode                        ! flags physics options wrt input data
  character(CL)                :: model_maskfile = nullstr        ! full pathname to obtain mask from
  real(R8)                     :: sst_constant_value
  integer                      :: aquap_option
  character(CL)                :: restfilm = nullstr              ! model restart file namelist
  character(CL)                :: restfils = nullstr              ! stream restart file namelist
  logical                      :: force_prognostic_true = .false. ! if true set prognostic true
  logical                      :: ocn_prognostic

  ! linked lists
  type(fldList_type) , pointer :: fldsImport => null()
  type(fldList_type) , pointer :: fldsExport => null()
  type(dfield_type)  , pointer :: dfields    => null()

  ! constants
  real(R8)                     :: dt                              ! real model timestep
  real(r8)     , parameter     :: cpsw    = shr_const_cpsw        ! specific heat of sea h2o ~ j/kg/k
  real(r8)     , parameter     :: rhosw   = shr_const_rhosw       ! density of sea water ~ kg/m^3
  real(r8)     , parameter     :: tkfrz   = shr_const_tkfrz       ! freezing point, fresh water (kelvin)
  real(r8)     , parameter     :: tkfrzsw = shr_const_tkfrzsw     ! freezing point, sea   water (kelvin)
  real(r8)     , parameter     :: latice  = shr_const_latice      ! latent heat of fusion
  real(r8)     , parameter     :: ocnsalt = shr_const_ocn_ref_sal ! ocean reference salinity
  integer      , parameter     :: master_task = 0                 ! task number of master task
  character(*) , parameter     :: rpfile = 'rpointer.ocn'
  character(*) , parameter     :: modName = "(ocn_comp_nuopc)"

  ! internal fields
  real(r8),         pointer :: xc(:), yc(:) ! mesh lats and lons - needed for aquaplanet analytical
  real(R8), public, pointer :: somtp(:)     ! SOM ocean temperature needed for restart

  ! export fields
  real(r8), pointer :: So_omask(:)  => null()    ! real ocean fraction sent to mediator
  real(r8), pointer :: So_t(:)      => null()
  real(r8), pointer :: So_s(:)      => null()
  real(r8), pointer :: So_u(:)      => null()
  real(r8), pointer :: So_v(:)      => null()
  real(r8), pointer :: So_dhdx(:)   => null()
  real(r8), pointer :: So_dhdy(:)   => null()
  real(r8), pointer :: So_fswpen(:) => null()
  real(r8), pointer :: Fioo_q(:)    => null()

  ! import  fields
  real(r8), pointer :: Foxx_swnet(:) => null()
  real(r8), pointer :: Foxx_lwup(:)  => null()
  real(r8), pointer :: Foxx_sen(:)   => null()
  real(r8), pointer :: Foxx_lat(:)   => null()
  real(r8), pointer :: Faxa_lwdn(:)  => null()
  real(r8), pointer :: Faxa_snow(:)  => null()
  real(r8), pointer :: Fioi_melth(:) => null()
  real(r8), pointer :: Foxx_rofi(:)  => null()

  ! internal stream type
  real(r8), pointer :: strm_h(:)    => null()
  real(r8), pointer :: strm_qbot(:) => null()

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
    integer           :: nu                 ! unit number
    integer           :: ierr               ! error code
    logical           :: exists             ! check for file existence  
    character(len=*),parameter  :: subname=trim(modName)//':(InitializeAdvertise) '
    !-------------------------------------------------------------------------------

    namelist / docn_nml / datamode, model_maskfile, &
         restfilm, restfils, force_prognostic_true, sst_constant_value

    rc = ESMF_SUCCESS

    ! Obtain flds_scalar values, mpi values, multi-instance values and
    ! set logunit and set shr logging to my log file
    call dshr_init(gcomp, mpicom, my_task, inst_index, inst_suffix, &
         flds_scalar_name, flds_scalar_num, flds_scalar_index_nx, flds_scalar_index_ny, &
         logunit, shrlogunit, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine namelist filename
    nlfilename = "docn_in"//trim(inst_suffix)

    ! Determine logical masterproc
    masterproc = (my_task == master_task)

    ! Read docn_nml from nlfilename
    if (my_task == master_task) then
       open (newunit=nu,file=trim(nlfilename),status="old",action="read")
       read (nu,nml=docn_nml,iostat=ierr)
       close(nu)
       if (ierr > 0) then
          write(logunit,*) 'ERROR: reading input namelist, '//trim(nlfilename)//' iostat=',ierr
          call shr_sys_abort(subName//': namelist read error '//trim(nlfilename))
       end if
       write(logunit,*)' restfilm   = ',trim(restfilm)
       write(logunit,*)' restfils   = ',trim(restfils)
       write(logunit,*)' force_prognostic_true = ',force_prognostic_true
       if (trim(model_maskfile) == nullstr) then
          write(logunit,*)' obtaining model mask from model mesh'
       else
          ! obtain model mask from model_maskfile
          inquire(file=trim(model_maskfile), exist=exists)
          if (.not.exists) then
             write(logunit, *)' ERROR: model_maskfile '//trim(model_maskfile)//' does not exist'
             call shr_sys_abort(trim(subname)//' ERROR: model_maskfile '//trim(model_maskfile)//' does not exist')
          else
             write(logunit,*)' obtaining model mask from ',trim(model_maskfile)
          end if
       end if
    endif
    call shr_mpi_bcast(datamode              , mpicom, 'datamode')
    call shr_mpi_bcast(model_maskfile        , mpicom, 'model_maskfile')
    call shr_mpi_bcast(restfilm              , mpicom, 'restfilm')
    call shr_mpi_bcast(restfils              , mpicom, 'restfils')
    call shr_mpi_bcast(force_prognostic_true , mpicom, 'force_prognostic_true')
    call shr_mpi_bcast(sst_constant_value    , mpicom, 'sst_constant_value')

    if (trim(datamode) /= 'NULL') then
       ! determine if ocn will receive import data
       if ( force_prognostic_true .or. trim(datamode) == 'IAF' .or. &
            trim(datamode) == 'SOM' .or. trim(datamode) == 'SOM_AQUAP') then
          ocn_prognostic = .true.
       else
          ocn_prognostic = .false.
       end if

       ! Special logic for prescribed aquaplanet
       if (datamode(1:9) == 'SST_AQUAP' .and. trim(datamode) /= 'SST_AQUAPFILE') then
          ! First determine the prescribed aquaplanet option
          if (len_trim(datamode) == 10) then
             read(datamode(10:10),'(i1)') aquap_option
          else if (len_trim(datamode) == 11) then
             read(datamode(10:11),'(i2)') aquap_option
          end if
          ! Now remove the index from the datamode value, to have a generic setting for later use
          datamode = "SST_AQUAPANAL"
       end if

       ! Validate datamode
       if ( trim(datamode) == 'NULL'          .or. trim(datamode) == 'COPYALL' .or.&
            trim(datamode) == 'SSTDATA'       .or. &
            trim(datamode) == 'SST_AQUAPANAL' .or. trim(datamode) == 'SST_AQUAPFILE' .or. &
            trim(datamode) == 'IAF'           .or. &
            trim(datamode) == 'SOM'           .or. trim(datamode) == 'SOM_AQUAP') then
          if (my_task == master_task) then
             write(logunit,*) ' docn datamode = ',trim(datamode)
          end if
       else
          call shr_sys_abort(' ERROR illegal docn datamode = '//trim(datamode))
       endif

       call docn_comp_advertise(importState, exportState, rc=rc)
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
    real(r8), allocatable   :: rmask(:)
    type(file_desc_t)       :: pioid
    type(var_desc_t)        :: varid
    type(io_desc_t)         :: pio_iodesc
    integer                 :: rcode
    logical                 :: reset_mask
    character(len=*), parameter :: subname=trim(modName)//':(InitializeRealize) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Reset shr logging to my log file
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_setLogUnit (logUnit)

    ! Initialize sdat and set the model domain mask in sdat if appropriate
    ! TODO: need a check that the mask file has the same grid as the model mesh
    call t_startf('docn_strdata_init')
    reset_mask = .false.
    if (datamode == 'SST_AQUAPANAL' .or. datamode == 'SST_AQUAPFILE' .or. datamode == 'SOM_AQUAP') then
       reset_mask = .true.
    end if
    if (trim(model_maskfile) /= nullstr) then 
       call dshr_sdat_init(gcomp, clock, nlfilename, compid, logunit, 'ocn', mesh, read_restart, sdat, &
            reset_mask=reset_mask, model_maskfile=model_maskfile, rc=rc)
    else
       call dshr_sdat_init(gcomp, clock, nlfilename, compid, logunit, 'ocn', mesh, read_restart, sdat, &
            reset_mask=reset_mask, rc=rc)
    end if
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('docn_strdata_init')

    ! Realize the actively coupled fields, now that a mesh is established and
    ! initialize dfields data type (to map streams to export state fields)
    call docn_comp_realize(importState, exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Read restart if necessary
    if (read_restart) then
       call dshr_restart_read(restfilm, restfils, rpfile, inst_suffix, nullstr, &
            logunit, my_task, mpicom, sdat, fld=somtp, fldname='somtp')
    end if

    ! Get the time to interpolate the stream data to
    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet(currTime, yy=current_year, mm=current_mon, dd=current_day, s=current_tod, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_cal_ymd2date(current_year, current_mon, current_day, current_ymd)

    ! Run docn
    call docn_comp_run(current_ymd, current_tod, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Add scalars to export state
    call dshr_state_SetScalar(dble(sdat%nxg),flds_scalar_index_nx, exportState, flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_SetScalar(dble(sdat%nyg),flds_scalar_index_ny, exportState, flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Diagnostics
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
    integer                 :: idt          ! integer model timestep
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

    ! Get model timestep
    call ESMF_ClockGet( clock, timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeIntervalGet( timeStep, s=idt, rc=rc )
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    dt = idt * 1.0_r8

    ! run docn
    call docn_comp_run(next_ymd, next_tod, rc=rc)
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
            logunit, mpicom, my_task, sdat, fld=somtp, fldname='somtp')
       call t_stopf('docn_restart')
    endif

    ! write diagnostics
    call dshr_state_diagnose(exportState,subname//':ES',rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (my_task == master_task) then
       call dshr_log_clock_advance(clock, 'docn', logunit, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    call shr_file_setLogUnit (shrlogunit)

  end subroutine ModelAdvance

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

  !===============================================================================

  subroutine docn_comp_advertise(importState, exportState, rc)

    ! --------------------------------------------------------------
    ! determine export and import fields to advertise to mediator
    ! --------------------------------------------------------------

    ! input/output arguments
    type(ESMF_State)     , intent(inout) :: importState
    type(ESMF_State)     , intent(inout) :: exportState
    integer              , intent(out)   :: rc

    ! local variables
    integer           :: n
    type(fldlist_type), pointer :: fldList
    !-------------------------------------------------------------------------------

    ! Advertise export fields

    call dshr_fldList_add(fldsExport, trim(flds_scalar_name))
    call dshr_fldList_add(fldsExport, 'So_omask'            )
    call dshr_fldList_add(fldsExport, 'So_t'                )
    call dshr_fldList_add(fldsExport, 'So_s'                )
    call dshr_fldList_add(fldsExport, 'So_u'                )
    call dshr_fldList_add(fldsExport, 'So_v'                )
    call dshr_fldList_add(fldsExport, 'So_dhdx'             )
    call dshr_fldList_add(fldsExport, 'So_dhdy'             )
    call dshr_fldList_add(fldsExport, 'Fioo_q'              )
    call dshr_fldList_add(fldsExport, 'So_fswpen'           )

    fldlist => fldsExport ! the head of the linked list
    do while (associated(fldlist))
       call NUOPC_Advertise(exportState, standardName=fldlist%stdname, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite('(docn_comp_advertise): Fr_ocn'//trim(fldList%stdname), ESMF_LOGMSG_INFO)
       fldList => fldList%next
    enddo

    ! Advertise import fields

    if (ocn_prognostic) then
       call dshr_fldList_add(fldsImport, trim(flds_scalar_name))
       call dshr_fldList_add(fldsImport, 'Foxx_swnet'          )
       call dshr_fldList_add(fldsImport, 'Foxx_lwup'           )
       call dshr_fldList_add(fldsImport, 'Foxx_sen'            )
       call dshr_fldList_add(fldsImport, 'Foxx_lat'            )
       call dshr_fldList_add(fldsImport, 'Faxa_lwdn'           )
       call dshr_fldList_add(fldsImport, 'Faxa_snow'           )
       call dshr_fldList_add(fldsImport, 'Fioi_melth'          )
       call dshr_fldList_add(fldsImport, 'Foxx_rofi'           )

       fldlist => fldsImport ! the head of the linked list
       do while (associated(fldlist))
          call NUOPC_Advertise(importState, standardName=fldlist%stdname, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_LogWrite('(docn_comp_advertise): Fr_ocn'//trim(fldList%stdname), ESMF_LOGMSG_INFO)
          fldList => fldList%next
       enddo
    end if

  end subroutine docn_comp_advertise

!===============================================================================

  subroutine docn_comp_realize(importState, exportState, rc)

    ! input/output parameters
    type(ESMF_State)       , intent(inout) :: importState
    type(ESMF_State)       , intent(inout) :: exportState
    integer                , intent(out)   :: rc

    ! local variables
    type(ESMF_StateItem_Flag) :: itemFlag
    integer                   :: n
    integer                   :: numOwnedElements   ! number of elements owned by this PET
    type(ESMF_DistGrid)       :: distGrid           ! mesh distGrid
    type(ESMF_Array)          :: elemMaskArray
    integer, pointer          :: imask(:)
    type(file_desc_t)         :: pioid
    type(var_desc_t)          :: varid
    type(io_desc_t)           :: pio_iodesc
    integer                   :: rcode
    character(*), parameter   :: subName = "(docn_comp_realize) "
    real(R8)    , parameter   :: &
         swp = 0.67_R8*(exp((-1._R8*shr_const_zsrflyr) /1.0_R8)) + 0.33_R8*exp((-1._R8*shr_const_zsrflyr)/17.0_R8)
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! -------------------------------------
    ! NUOPC_Realize "realizes" a previously advertised field in the importState and exportState
    ! by replacing the advertised fields with the newly created fields of the same name.
    ! -------------------------------------

    call dshr_fldlist_realize( exportState, fldsExport, flds_scalar_name, flds_scalar_num, mesh, &
         subname//':docnExport', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_fldlist_realize( importState, fldsImport, flds_scalar_name, flds_scalar_num, mesh, &
         subname//':docnImport', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! -------------------------------------
    ! Determine ocean fraction 
    ! -------------------------------------

    ! Set pointers to exportState fields that have no corresponding stream field
    call dshr_state_getfldptr(exportState, fldname='So_omask', fldptr1=So_omask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Obtain So_omask (the ocean fraction)
    if (trim(model_maskfile) /= nullstr) then 
       ! Read in the ocean fraction from the input namelist ocean mask file and assume 'frac' name on domain file
       rcode = pio_openfile(sdat%pio_subsystem, pioid, sdat%io_type, trim(model_maskfile), pio_nowrite)
       call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)
       rcode = pio_inq_varid(pioid, 'frac', varid) 
       call pio_initdecomp(sdat%pio_subsystem, pio_double, (/sdat%nxg, sdat%nyg/), sdat%gindex, pio_iodesc)
       call pio_read_darray(pioid, varid, pio_iodesc, So_omask, rcode)
       call pio_closefile(pioid)
       call pio_freedecomp(sdat%pio_subsystem, pio_iodesc)
    else
       ! Obtain the ocean fraction from the mask values in the ocean mesh file
       call ESMF_MeshGet(mesh, numOwnedElements=numOwnedElements, elementdistGrid=distGrid, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       allocate(imask(numOwnedElements))
       elemMaskArray = ESMF_ArrayCreate(distGrid, imask, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! the following call sets the varues of imask
       call ESMF_MeshGet(mesh, elemMaskArray=elemMaskArray, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! now set the fraction as just the real mask
       So_omask(:) = real(imask(:), kind=r8)
       deallocate(imask)
       call ESMF_ArrayDestroy(elemMaskArray, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! -------------------------------------
    ! Set pointers to exportState fields
    ! -------------------------------------

    call dshr_state_getfldptr(exportState, fldname='Fioo_q', fldptr1=Fioo_q, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    Fioo_q(:) = 0._r8

    ! Initialize export state data that has corresponding stream field
    call dshr_dfield_add(dfields, sdat, state_fld='So_t', strm_fld='t', &
         state=exportState, state_ptr=So_t, logunit=logunit, masterproc=masterproc, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add(dfields, sdat, state_fld='So_s', strm_fld='s', &
         state=exportState, state_ptr=So_s, logunit=logunit, masterproc=masterproc, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add(dfields, sdat,  state_fld='So_u', strm_fld='u', &
         state=exportState, state_ptr=So_u, logunit=logunit, masterproc=masterproc, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add(dfields, sdat,  state_fld='So_v', strm_fld='v', &
         state=exportState, state_ptr=So_v, logunit=logunit, masterproc=masterproc, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add(dfields, sdat, state_fld='So_dhdx', strm_fld='dhdx', &
         state=exportState, state_ptr=So_dhdx, logunit=logunit, masterproc=masterproc, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add(dfields, sdat, state_fld='So_dhdy', strm_fld='dhdy', &
         state=exportState, state_ptr=So_dhdy, logunit=logunit, masterproc=masterproc, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Initialize dfields stream fields that have no corresponding export fields
    call dshr_dfield_add(dfields, sdat,  strm_fld='qbot', strm_ptr=strm_qbot)
    call dshr_dfield_add(dfields, sdat,  strm_fld='h'   , strm_ptr=strm_h)

    ! For So_fswpen is only needed for diurnal cycle calculation of atm/ocn fluxes - and
    ! currently this is not implemented in cmeps
    call ESMF_StateGet(exportState, 'So_fswpen', itemFlag, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (itemFlag /= ESMF_STATEITEM_NOTFOUND) then
       call dshr_state_getfldptr(exportState, 'So_fswpen', fldptr1=So_fswpen, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       So_fswpen(:) = swp
    end if

    ! Initialize export state pointers to non-zero
    So_t(:) = TkFrz
    So_s(:) = ocnsalt

    ! Allocate memory for somtp
    allocate(somtp(sdat%lsize))

    ! -------------------------------------
    ! Set pointers to importState fields
    ! -------------------------------------

    if (ocn_prognostic) then
       call dshr_state_getfldptr(importState, 'Foxx_swnet' , fldptr1=Foxx_swnet , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call dshr_state_getfldptr(importState, 'Foxx_lwup'  , fldptr1=Foxx_lwup  , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call dshr_state_getfldptr(importState, 'Foxx_lwup'  , fldptr1=Foxx_lwup  , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call dshr_state_getfldptr(importState, 'Foxx_sen'   , fldptr1=Foxx_sen   , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call dshr_state_getfldptr(importState, 'Foxx_lat'   , fldptr1=Foxx_lat   , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call dshr_state_getfldptr(importState, 'Faxa_lwdn'  , fldptr1=Faxa_lwdn  , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call dshr_state_getfldptr(importState, 'Faxa_snow'  , fldptr1=Faxa_snow  , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call dshr_state_getfldptr(importState, 'Fioi_melth' , fldptr1=Fioi_melth , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call dshr_state_getfldptr(importState, 'Foxx_rofi'  , fldptr1=Foxx_rofi  , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

  end subroutine docn_comp_realize

!===============================================================================

  subroutine docn_comp_run(target_ymd, target_tod, rc)

    ! --------------------------
    ! advance docn
    ! --------------------------

    ! input/output variables:
    integer , intent(in)    :: target_ymd       ! model date
    integer , intent(in)    :: target_tod       ! model sec into model date
    integer , intent(out)   :: rc

    ! local variables
    integer               :: n, lsize
    logical               :: first_time = .true.
    integer               :: spatialDim         ! number of dimension in mesh
    integer               :: numOwnedElements   ! size of mesh
    real(r8), pointer     :: ownedElemCoords(:) ! mesh lat and lons
    real(r8), allocatable :: tfreeze(:)         ! SOM ocean freezing temperature
    character(*), parameter :: subName = "(docn_comp_run) "
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('DOCN_RUN')

    !--------------------
    ! advance docn streams
    !--------------------

    ! time and spatially interpolate to model time and grid
    call t_barrierf('docn_BARRIER',mpicom)
    call t_startf('docn_strdata_advance')
    call shr_strdata_advance(sdat, target_ymd, target_tod, mpicom, 'docn')
    call t_stopf('docn_strdata_advance')

    !--------------------
    ! copy all fields from streams to export state as default
    !--------------------

    ! This automatically will update the fields in the export state
    call t_barrierf('docn_dfield_copy_BARRIER', mpicom)
    call t_startf('docn_dfield_copy')
    call dshr_dfield_copy(dfields, sdat, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('docn_dfield_copy')

    !-------------------------------------------------
    ! Determine additional data model behavior based on the mode
    !-------------------------------------------------

    lsize = sdat%lsize

    call t_startf('docn_datamode')
    select case (trim(datamode))

    case('COPYALL')
       ! do nothing extra

    case('SSTDATA')
       So_t(:) = So_t(:) + TkFrz

    case('SST_AQUAPANAL')
       So_s(:)      = 0.0_r8
       if (associated(So_fswpen)) then
          So_fswpen(:) = 0.0_r8
       end if
       if (first_time) then
          call ESMF_MeshGet(mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          allocate(ownedElemCoords(spatialDim*numOwnedElements))
          allocate(xc(numOwnedElements), yc(numOwnedElements))
          call ESMF_MeshGet(mesh, ownedElemCoords=ownedElemCoords)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          do n = 1,numOwnedElements
             xc(n) = ownedElemCoords(2*n-1)
             yc(n) = ownedElemCoords(2*n)
          end do
       end if
       call docn_prescribed_sst(xc, yc, lsize, aquap_option, So_t)
       So_t(:) = So_t(:) + TkFrz

    case('SST_AQUAPFILE')
       So_s(:)      = 0.0_r8
       if (associated(So_fswpen)) then
          So_fswpen(:) = 0.0_r8
       end if
       So_t(:) = So_t(:) + TkFrz

    case('SST_AQUAP_CONSTANT')
       So_s(:)      = 0.0_r8
       if (associated(So_fswpen)) then
          So_fswpen(:) = 0.0_r8
       end if
       So_t(:) = sst_constant_value

    case('IAF')
       So_t(:) = So_t(:) + TkFrz

    case('SOM')
       if (first_time) then
          do n = 1,lsize
             if (.not. read_restart) then
                somtp(n) = So_t(n) + TkFrz
             endif
             So_t(n) = somtp(n)
             Fioo_q(n) = 0.0_R8
          enddo
       else
          allocate(tfreeze(lsize))
          tfreeze(:) = shr_frz_freezetemp(So_s(:)) + TkFrz
          do n = 1,lsize
             if (So_omask(n) /= 0._r8) then
                ! compute new temp (last term is latent by prec and roff)
                So_t(n) = somtp(n) +  &
                     ( Foxx_swnet(n) + Foxx_lwup(n) + Faxa_lwdn(n) + Foxx_sen(n) + Foxx_lat(n) + &
                       Fioi_melth(n) - strm_qbot(n) - (Faxa_snow(n)+Foxx_rofi(n))*latice) * dt/(cpsw*rhosw* strm_h(n))

                ! compute ice formed or melt potential
                Fioo_q(n) = (tfreeze(n) - So_t(n))*(cpsw*rhosw*strm_h(n))/dt ! ice formed q>0

                ! reset temp
                So_t(n)  = max(tfreeze(n),So_t(n))

                ! save somtp to restart file
                somtp(n) = So_t(n)
             endif
          end do
          deallocate(tfreeze)
       endif   ! first_time

    case('SOM_AQUAP')
       if (first_time) then
          do n = 1,lsize
             if (.not. read_restart) then
                somtp(n) = So_t(n) + TkFrz
             endif
             So_t(n) = somtp(n)
             Fioo_q(n) = 0.0_R8
          enddo
       else
          allocate(tfreeze(lsize))
          tfreeze(:) = shr_frz_freezetemp(So_s(:)) + TkFrz
          do n = 1,lsize
             ! compute new temp (last term is latent by prec and roff)
             So_t(n) = somtp(n) + &
                  ( Foxx_swnet(n) + Foxx_lwup(n) + Faxa_lwdn(n) + Foxx_sen(n) + Foxx_lat(n) + &
                    Fioi_melth(n) - strm_qbot(n) - (Faxa_snow(n)+Foxx_rofi(n))*latice) * dt/(cpsw*rhosw*strm_h(n))

             ! compute ice formed or melt potential
             Fioo_q(n) = (tfreeze(n) - So_t(n))*(cpsw*rhosw*strm_h(n))/dt  ! ice formed q>0

             ! save somtp on restart file
             somtp(n) = So_t(n)
          enddo
          deallocate(tfreeze)
       endif   ! first_time

    end select

    first_time= .false.

    call t_stopf('docn_datamode')
    call t_stopf('DOCN_RUN')

  end subroutine docn_comp_run

!===============================================================================

  subroutine docn_prescribed_sst(xc, yc, lsize, sst_option, sst)

    ! input/output variables
    real(R8)     , intent(in)    :: xc(:)  !degrees
    real(R8)     , intent(in)    :: yc(:)  !degrees
    integer      , intent(in)    :: lsize
    integer      , intent(in)    :: sst_option
    real(R8)     , intent(inout) :: sst(:)

    ! local variables
    integer  :: i
    real(r8) :: tmp, tmp1, pi
    real(r8) :: rlon(lsize), rlat(lsize)

    real(r8), parameter :: pio180 = SHR_CONST_PI/180._r8

    ! Parameters for zonally symmetric experiments
    real(r8), parameter ::   t0_max     = 27._r8
    real(r8), parameter ::   t0_min     = 0._r8
    real(r8), parameter ::   maxlat     = 60._r8*pio180
    real(r8), parameter ::   shift      = 5._r8*pio180
    real(r8), parameter ::   shift9     = 10._r8*pio180
    real(r8), parameter ::   shift10    = 15._r8*pio180

    ! Parameters for zonally asymmetric experiments
    real(r8), parameter ::   t0_max6    = 1._r8
    real(r8), parameter ::   t0_max7    = 3._r8
    real(r8), parameter ::   latcen     = 0._r8*pio180
    real(r8), parameter ::   loncen     = 0._r8*pio180
    real(r8), parameter ::   latrad6    = 15._r8*pio180
    real(r8), parameter ::   latrad8    = 30._r8*pio180
    real(r8), parameter ::   lonrad     = 30._r8*pio180
    !-------------------------------------------------------------------------------

    pi = SHR_CONST_PI

    ! convert xc and yc from degrees to radians

    rlon(:) = xc(:) * pio180
    rlat(:) = yc(:) * pio180

    ! Control
    if (sst_option < 1 .or. sst_option > 10) then
       call shr_sys_abort ('docn_prescribed_sst: ERROR: sst_option must be between 1 and 10')
    end if

    if (sst_option == 1 .or. sst_option == 6 .or. sst_option == 7 .or. sst_option == 8) then
       do i = 1,lsize
          if (abs(rlat(i)) > maxlat) then
             sst(i) = t0_min
          else
             tmp = sin(rlat(i)*pi*0.5_r8/maxlat)
             tmp = 1._r8 - tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          end if
       end do
    end if
    if (sst_option == 2) then ! Flat
       do i = 1,lsize
          if (abs(rlat(i)) > maxlat) then
             sst(i) = t0_min
          else
             tmp = sin(rlat(i)*pi*0.5_r8/maxlat)
             tmp = 1._r8 - tmp*tmp*tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          end if
       end do
    end if
    if (sst_option == 3) then ! Qobs
       do i = 1,lsize
          if (abs(rlat(i)) > maxlat) then
             sst(i) = t0_min
          else
             tmp = sin(rlat(i)*pi*0.5_r8/maxlat)
             tmp = (2._r8 - tmp*tmp*tmp*tmp - tmp*tmp)*0.5_r8
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          end if
       end do
    end if
    if (sst_option == 4) then ! Peaked
       do i = 1,lsize
          if (abs(rlat(i)) > maxlat) then
             sst(i) = t0_min
          else
             tmp = (maxlat - abs(rlat(i)))/maxlat
             tmp1 = 1._r8 - tmp
             sst(i) = t0_max*tmp + t0_min*tmp1
          end if
       end do
    end if
    if (sst_option == 5) then ! Control-5N
       do i = 1,lsize
          if (abs(rlat(i)) > maxlat) then
             sst(i) = t0_min
          else if (rlat(i) > shift) then
             tmp = sin((rlat(i)-shift)*pi*0.5_r8/(maxlat-shift))
             tmp = 1._r8 - tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          else
             tmp = sin((rlat(i)-shift)*pi*0.5_r8/(maxlat+shift))
             tmp = 1._r8 - tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          end if
       end do
    end if
    if (sst_option == 6) then ! 1KEQ
       do i = 1,lsize
          if (abs(rlat(i)-latcen) <= latrad6) then
             tmp1 = cos((rlat(i)-latcen)*pi*0.5_r8/latrad6)
             tmp1 = tmp1*tmp1
             tmp = abs(rlon(i)-loncen)
             tmp = min(tmp , 2._r8*pi-tmp)
             if(tmp <= lonrad) then
                tmp = cos(tmp*pi*0.5_r8/lonrad)
                tmp = tmp*tmp
                sst(i) = sst(i) + t0_max6*tmp*tmp1
             end if
          end if
       end do
    end if
    if (sst_option == 7) then ! 3KEQ
       do i = 1, lsize
          if (abs(rlat(i)-latcen) <= latrad6) then
             tmp1 = cos((rlat(i)-latcen)*pi*0.5_r8/latrad6)
             tmp1 = tmp1*tmp1
             tmp = abs(rlon(i)-loncen)
             tmp = min(tmp , 2._r8*pi-tmp)
             if (tmp <= lonrad) then
                tmp = cos(tmp*pi*0.5_r8/lonrad)
                tmp = tmp*tmp
                sst(i) = sst(i) + t0_max7*tmp*tmp1
             end if
          end if
       end do
    end if
    if (sst_option == 8) then ! 3KW1
       do i = 1, lsize
          if (abs(rlat(i)-latcen) <= latrad8) then
             tmp1 = cos((rlat(i)-latcen)*pi*0.5_r8/latrad8)
             tmp1 = tmp1*tmp1
             tmp = cos(rlon(i)-loncen)
             sst(i) = sst(i) + t0_max7*tmp*tmp1
          end if
       end do
    end if
    if (sst_option == 9) then ! Control-10N
       do i = 1, lsize
          if (abs(rlat(i)) > maxlat) then
             sst(i) = t0_min
          else if (rlat(i) > shift9) then
             tmp = sin((rlat(i)-shift9)*pi*0.5_r8/(maxlat-shift9))
             tmp = 1._r8 - tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          else
             tmp = sin((rlat(i)-shift9)*pi*0.5_r8/(maxlat+shift9))
             tmp = 1._r8 - tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          end if
       end do
    end if
    if (sst_option == 10) then ! Control-15N
       do i = 1, lsize
          if (abs(rlat(i)) > maxlat) then
             sst(i) = t0_min
          else if(rlat(i) > shift10) then
             tmp = sin((rlat(i)-shift10)*pi*0.5_r8/(maxlat-shift10))
             tmp = 1._r8 - tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          else
             tmp = sin((rlat(i)-shift10)*pi*0.5_r8/(maxlat+shift10))
             tmp = 1._r8 - tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          end if
       end do
    end if

  end subroutine docn_prescribed_sst

end module ocn_comp_nuopc
