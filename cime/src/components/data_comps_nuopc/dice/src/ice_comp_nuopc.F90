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
  use shr_file_mod         , only : shr_file_getlogunit, shr_file_setlogunit
  use shr_kind_mod         , only : r8=>shr_kind_r8, cxx=>shr_kind_cxx, cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_const_mod        , only : shr_const_pi, shr_const_spval, shr_const_tkfrz, shr_const_latice
  use shr_sys_mod          , only : shr_sys_abort
  use shr_cal_mod          , only : shr_cal_noleap, shr_cal_gregorian, shr_cal_ymd2date, shr_cal_ymd2julian
  use shr_mpi_mod          , only : shr_mpi_bcast
  use shr_frz_mod          , only : shr_frz_freezetemp
  use dshr_methods_mod     , only : dshr_state_getfldptr, dshr_state_diagnose, chkerr, memcheck
  use dshr_mod             , only : dshr_model_initphase, dshr_init, dshr_sdat_init 
  use dshr_mod             , only : dshr_state_setscalar, dshr_set_runclock, dshr_log_clock_advance
  use dshr_mod             , only : dshr_restart_read, dshr_restart_write
  use dshr_mod             , only : dshr_create_mesh_from_grid
  use dshr_strdata_mod     , only : shr_strdata_type, shr_strdata_advance
  use dshr_dfield_mod      , only : dfield_type, dshr_dfield_add, dshr_dfield_copy
  use dshr_fldlist_mod     , only : fldlist_type, dshr_fldlist_add, dshr_fldlist_realize
  use dice_flux_atmice_mod , only : dice_flux_atmice
  use perf_mod             , only : t_startf, t_stopf, t_adj_detailf, t_barrierf
  use pio

  implicit none
  private ! except

  public  :: SetServices

  private :: InitializeAdvertise
  private :: InitializeRealize
  private :: ModelAdvance
  private :: ModelFinalize
  private :: dice_comp_advertise
  private :: dice_comp_realize
  private :: dice_comp_run

  !--------------------------------------------------------------------------
  ! Module data
  !--------------------------------------------------------------------------

  type(shr_strdata_type)       :: sdat
  type(ESMF_Mesh)              :: mesh                        ! model mesh
  character(len=CS)            :: flds_scalar_name = ''
  integer                      :: flds_scalar_num = 0
  integer                      :: flds_scalar_index_nx = 0
  integer                      :: flds_scalar_index_ny = 0
  integer                      :: compid                      ! mct comp id
  integer                      :: mpicom                      ! mpi communicator
  integer                      :: my_task                     ! my task in mpi communicator mpicom
  logical                      :: masterproc                  ! true of my_task == master_task
  character(len=16)            :: inst_suffix = ""            ! char string associated with instance (ie. "_0001" or "")
  integer                      :: logunit                     ! logging unit number
  logical                      :: read_restart                ! start from restart
  character(*) , parameter     :: nullstr = 'undefined'

  ! dice_in namelist input
  character(CL)                :: nlfilename                  ! filename to obtain namelist info from
  character(CL)                :: dataMode                    ! flags physics options wrt input data
  character(CL)                :: model_maskfile = nullstr    ! full pathname to obtain mask from
  real(R8)                     :: flux_swpf                   ! short-wave penatration factor
  real(R8)                     :: flux_Qmin                   ! bound on melt rate
  logical                      :: flux_Qacc                   ! activates water accumulation/melt wrt Q
  real(R8)                     :: flux_Qacc0                  ! initial water accumulation value
  character(CL)                :: restfilm = nullstr          ! model restart file namelist
  character(CL)                :: restfils = nullstr          ! stream restart file namelist

  ! nuopc attributes
  logical                      :: flds_i2o_per_cat            ! .true. if select per ice thickness
  character(CS)                :: calendar                    ! calendar name

  ! linked lists
  type(fldList_type) , pointer :: fldsImport => null()
  type(fldList_type) , pointer :: fldsExport => null()
  type(dfield_type)  , pointer :: dfields    => null()

  ! constants
  real(R8)                     :: dt                          ! real model timestep
  real(R8)     , parameter     :: pi       = shr_const_pi     ! pi
  real(r8)     , parameter     :: spval    = shr_const_spval  ! flags invalid data
  real(r8)     , parameter     :: tFrz     = shr_const_tkfrz  ! temp of freezing
  real(r8)     , parameter     :: latice   = shr_const_latice ! latent heat of fusion
  real(r8)     , parameter     :: waterMax = 1000.0_r8        ! wrt iFrac comp & frazil ice (kg/m^2)
  integer      , parameter     :: master_task=0               ! task number of master task
  character(*) , parameter     :: rpfile = 'rpointer.ice'
  character(*) , parameter     :: modName =  "(ice_comp_nuopc)"

  ! surface albedo constants
  real(r8),parameter :: snwfrac = 0.286_r8 ! snow cover fraction ~ [0,1]
  real(r8),parameter :: as_nidf = 0.950_r8 ! albedo: snow,near-infr,diffuse
  real(r8),parameter :: as_vsdf = 0.700_r8 ! albedo: snow,visible  ,diffuse
  real(r8),parameter :: as_nidr = 0.960_r8 ! albedo: snow,near-infr,direct
  real(r8),parameter :: as_vsdr = 0.800_r8 ! albedo: snow,visible  ,direct
  real(r8),parameter :: ai_nidf = 0.700_r8 ! albedo: ice, near-infr,diffuse
  real(r8),parameter :: ai_vsdf = 0.500_r8 ! albedo: ice, visible  ,diffuse
  real(r8),parameter :: ai_nidr = 0.700_r8 ! albedo: ice, near-infr,direct
  real(r8),parameter :: ai_vsdr = 0.500_r8 ! albedo: ice, visible  ,direct
  real(r8),parameter :: ax_nidf = ai_nidf*(1.0_r8-snwfrac) + as_nidf*snwfrac
  real(r8),parameter :: ax_vsdf = ai_vsdf*(1.0_r8-snwfrac) + as_vsdf*snwfrac
  real(r8),parameter :: ax_nidr = ai_nidr*(1.0_r8-snwfrac) + as_nidr*snwfrac
  real(r8),parameter :: ax_vsdr = ai_vsdr*(1.0_r8-snwfrac) + as_vsdr*snwfrac

  ! restart fields
  real(r8), pointer, public :: water(:) => null()

  ! internal fields
  real(r8), pointer :: yc(:)      => null() ! mesh lats (degrees)
  integer , pointer :: imask(:)   => null()
  real(r8), pointer :: tfreeze(:) => null()
  !real(r8), pointer:: ifrac0(:)  => null()

  ! export fields
  real(r8), pointer ::  Si_imask(:)      => null()
  real(r8), pointer ::  Si_ifrac(:)      => null()
  real(r8), pointer ::  Si_t(:)          => null()
  real(r8), pointer ::  Si_tref(:)       => null()
  real(r8), pointer ::  Si_qref(:)       => null()
  real(r8), pointer ::  Si_avsdr(:)      => null()
  real(r8), pointer ::  Si_anidr(:)      => null()
  real(r8), pointer ::  Si_avsdf(:)      => null()
  real(r8), pointer ::  Si_anidf(:)      => null()
  real(r8), pointer ::  Faii_swnet(:)    => null()
  real(r8), pointer ::  Faii_sen(:)      => null()
  real(r8), pointer ::  Faii_lat(:)      => null()
  real(r8), pointer ::  Faii_lwup(:)     => null()
  real(r8), pointer ::  Faii_evap(:)     => null()
  real(r8), pointer ::  Faii_taux(:)     => null()
  real(r8), pointer ::  Faii_tauy(:)     => null()
  real(r8), pointer ::  Fioi_melth(:)    => null()
  real(r8), pointer ::  Fioi_meltw(:)    => null()
  real(r8), pointer ::  Fioi_swpen(:)    => null()
  real(r8), pointer ::  Fioi_taux(:)     => null()
  real(r8), pointer ::  Fioi_tauy(:)     => null()
  real(r8), pointer ::  Fioi_salt(:)     => null()
  real(r8), pointer ::  Fioi_bcpho(:)    => null()
  real(r8), pointer ::  Fioi_bcphi(:)    => null()
  real(r8), pointer ::  Fioi_flxdst(:)   => null()
  real(r8), pointer ::  Si_ifrac_n(:,:)  => null()
  real(r8), pointer ::  Fioi_swpen_ifrac_n(:,:) => null()

  ! import fields
  real(r8), pointer :: Faxa_swvdr(:)    => null()
  real(r8), pointer :: Faxa_swvdf(:)    => null()
  real(r8), pointer :: Faxa_swndr(:)    => null()
  real(r8), pointer :: Faxa_swndf(:)    => null()
  real(r8), pointer :: Fioo_q(:)        => null()
  real(r8), pointer :: Sa_z(:)          => null()
  real(r8), pointer :: Sa_u(:)          => null()
  real(r8), pointer :: Sa_v(:)          => null()
  real(r8), pointer :: Sa_ptem(:)       => null()
  real(r8), pointer :: Sa_shum(:)       => null()
  real(r8), pointer :: Sa_dens(:)       => null()
  real(r8), pointer :: Sa_tbot(:)       => null()
  real(r8), pointer :: So_s(:)          => null()
  real(r8), pointer :: Faxa_bcph(:,:)   => null()
  real(r8), pointer :: Faxa_ocph(:,:)   => null()
  real(r8), pointer :: Faxa_dstdry(:,:) => null()
  real(r8), pointer :: Faxa_dstwet(:,:) => null()

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
    integer           :: shrlogunit         ! original log unit
    integer           :: nu                 ! unit number
    integer           :: ierr               ! error code
    logical           :: exists             ! check for file existence  
    character(len=*),parameter  :: subname=trim(modName)//':(InitializeAdvertise) '
    !-------------------------------------------------------------------------------

    namelist / dice_nml / datamode, model_maskfile, restfilm, restfils, &
         flux_swpf, flux_Qmin, flux_Qacc, flux_Qacc0 

    rc = ESMF_SUCCESS

    ! Obtain flds_scalar values, mpi values, multi-instance values and  
    ! set logunit and set shr logging to my log file
    call dshr_init(gcomp, mpicom, my_task, inst_index, inst_suffix, &
         flds_scalar_name, flds_scalar_num, flds_scalar_index_nx, flds_scalar_index_ny, &
         logunit, shrlogunit, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine namelist filename
    nlfilename = "dice_in"//trim(inst_suffix)

    ! determine logical masterproc
    masterproc = (my_task == master_task)

    ! Read dice_nml from nlfilename
    if (my_task == master_task) then
       open (newunit=nu,file=trim(nlfilename),status="old",action="read")
       read (nu,nml=dice_nml,iostat=ierr)
       close(nu)
       if (ierr > 0) then
          write(logunit,*) 'ERROR: reading input namelist, '//trim(nlfilename)//' iostat=',ierr
          call shr_sys_abort(subName//': namelist read error '//trim(nlfilename))
       end if
       write(logunit,*)' flux_swpf  = ',flux_swpf
       write(logunit,*)' flux_Qmin  = ',flux_Qmin
       write(logunit,*)' flux_Qacc  = ',flux_Qacc
       write(logunit,*)' flux_Qacc0 = ',flux_Qacc0
       write(logunit,*)' restfilm   = ',trim(restfilm)
       write(logunit,*)' restfils   = ',trim(restfils)
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
    call shr_mpi_bcast(datamode       , mpicom, 'datamode')
    call shr_mpi_bcast(model_maskfile , mpicom, 'model_maskfile')
    call shr_mpi_bcast(restfilm       , mpicom, 'restfilm')
    call shr_mpi_bcast(restfils       , mpicom, 'restfils')
    call shr_mpi_bcast(flux_swpf      , mpicom, 'flux_swpf')
    call shr_mpi_bcast(flux_Qmin      , mpicom, 'flux_Qmin')
    call shr_mpi_bcast(flux_Qacc      , mpicom, 'flux_Qacc')
    call shr_mpi_bcast(flux_Qacc0     , mpicom, 'flux_Qacc0')

    ! Validate datamode
    if (trim(datamode) == 'NULL' .or. trim(datamode) == 'SSTDATA' .or. trim(datamode) == 'COPYALL') then
       if (my_task == master_task) write(logunit,*) ' dice datamode = ',trim(datamode)
    else
       call shr_sys_abort(' ERROR illegal dice datamode = '//trim(datamode))
    endif

    ! Advertise import and export fields
    if (datamode /= 'NULL') then
       call NUOPC_CompAttributeGet(gcomp, name='flds_i2o_per_cat', value=cvalue, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) flds_i2o_per_cat  ! module variable

       call dice_comp_advertise(importState, exportState, rc=rc)
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
    type(ESMF_TimeInterval)     :: TimeStep
    type(ESMF_Time)             :: currTime
    type(ESMF_Calendar)         :: esmf_calendar ! esmf calendar
    type(ESMF_CalKind_Flag)     :: esmf_caltype  ! esmf calendar type
    integer                     :: current_ymd   ! model date
    integer                     :: current_year  ! model year
    integer                     :: current_mon   ! model month
    integer                     :: current_day   ! model day
    integer                     :: current_tod   ! model sec into model date
    real(R8)                    :: cosarg        ! for setting ice temp pattern
    real(R8)                    :: jday, jday0   ! elapsed day counters
    character(CL)               :: cvalue        ! temporary
    integer                     :: shrlogunit    ! original log unit
    integer                     :: n,k           ! generic counters
    integer                     :: model_dt      ! integer model timestep
    character(len=*), parameter :: F00   = "('ice_comp_nuopc: ')',8a)"
    character(len=*), parameter :: subname=trim(modName)//':(InitializeRealize) '
    !-------------------------------------------------------------------------------

    if (datamode == 'NULL') RETURN

    rc = ESMF_SUCCESS

    ! Reset shr logging to my log file
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_setLogUnit (logUnit)

    ! Initialize sdat
    call t_startf('dice_strdata_init')
    if (trim(model_maskfile) /= nullstr) then 
       call dshr_sdat_init(gcomp, clock, nlfilename, compid, logunit, 'ice', mesh, read_restart, sdat, &
            model_maskfile=model_maskfile, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call dshr_sdat_init(gcomp, clock, nlfilename, compid, logunit, 'ice', mesh, read_restart, sdat, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    call t_stopf('dice_strdata_init')

    ! Realize the actively coupled fields, now that a mesh is established and
    ! initialize dfields data type (to map streams to export state fields)
    call dice_comp_realize(importState, exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Read restart if necessary
    if (read_restart) then
       call dshr_restart_read(restfilm, restfils, rpfile, inst_suffix, nullstr, &
            logunit, my_task, mpicom, sdat, fld=water, fldname='water')
    end if

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
    call shr_cal_ymd2julian(0, current_mon, current_day, current_tod, jDay , sdat%calendar) ! julian day for model
    call shr_cal_ymd2julian(0, 9,           1,           0,           jDay0, sdat%calendar) ! julian day for Sept 1
    cosArg = 2.0_R8*pi*(jday - jday0)/365.0_R8

    ! Run dice
    call dice_comp_run(current_ymd, current_tod, cosarg, rc)

    ! Add scalars to export state
    call dshr_state_SetScalar(dble(sdat%nxg),flds_scalar_index_nx, exportState, flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_SetScalar(dble(sdat%nyg),flds_scalar_index_ny, exportState, flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Diagnostics
    call dshr_state_diagnose(exportState,subname//':ES',rc=rc)
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
    integer                 :: current_mon   ! model month
    integer                 :: current_day   ! model day
    integer                 :: current_tod   ! model sec into model date
    real(R8)                :: cosarg        ! for setting ice temp pattern
    real(R8)                :: jday, jday0   ! elapsed day counters
    integer                 :: shrlogunit    ! original log unit
    integer                 :: next_ymd      ! model date
    integer                 :: next_tod      ! model sec into model date
    integer                 :: yr            ! year
    integer                 :: mon           ! month
    integer                 :: day           ! day in month
    character(CL)           :: case_name     ! case name
    character(len=*),parameter :: subname=trim(modName)//':(ModelAdvance) '
    !-------------------------------------------------------------------------------

    if (datamode == 'NULL') RETURN

    rc = ESMF_SUCCESS

    call t_startf(subname)
    call memcheck(subname, 5, my_task == master_task)

    ! Reset shr logging to my log file
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_setLogUnit (logunit)

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
    call shr_cal_ymd2julian(0, mon, day, next_tod, jDay , sdat%calendar)    ! julian day for model
    call shr_cal_ymd2julian(0, 9,   1,   0,        jDay0, sdat%calendar)    ! julian day for Sept 1
    cosArg = 2.0_R8*pi*(jday - jday0)/365.0_R8

    ! Run dice
    call dice_comp_run(next_ymd, next_tod, cosarg, rc)

    ! Write_restart if alarm is ringing
    call ESMF_ClockGetAlarm(clock, alarmname='alarm_restart', alarm=alarm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (ESMF_AlarmIsRinging(alarm, rc=rc)) then
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_AlarmRingerOff( alarm, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call t_startf('dice_restart')
       call NUOPC_CompAttributeGet(gcomp, name='case_name', value=case_name, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call dshr_restart_write(rpfile, case_name, 'dice', inst_suffix, next_ymd, next_tod, &
            logunit, mpicom, my_task, sdat, fld=water, fldname='water')
       call t_stopf('dice_restart')
    endif

    ! Write diagnostics
    call dshr_state_diagnose(exportState,subname//':ES',rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (my_task == master_task) then
       call dshr_log_clock_advance(clock, 'dice', logunit, rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

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
       write(logunit,*) 'dice : end of main integration loop'
       write(logunit,*)
    end if

  end subroutine ModelFinalize

  !===============================================================================
  subroutine dice_comp_advertise(importState, exportState, rc)

    ! --------------------------------------------------------------
    ! determine export and import fields to advertise to mediator
    ! --------------------------------------------------------------

    ! input/output arguments
    type(ESMF_State)     , intent(inout) :: importState
    type(ESMF_State)     , intent(inout) :: exportState
    integer              , intent(out)   :: rc

    ! local variables
    integer         :: n
    type(fldlist_type), pointer :: fldList
    !-------------------------------------------------------------------------------

    ! Advertise export fields

    call dshr_fldList_add(fldsExport , trim(flds_scalar_name))
    call dshr_fldList_add(fldsExport ,'Si_ifrac'    )
    call dshr_fldList_add(fldsExport ,'Si_imask'    )
    call dshr_fldList_add(fldsExport ,'Si_t'        )
    call dshr_fldList_add(fldsExport ,'Si_tref'     )
    call dshr_fldList_add(fldsExport ,'Si_qref'     )
    call dshr_fldList_add(fldsExport ,'Si_avsdr'    )
    call dshr_fldList_add(fldsExport ,'Si_anidr'    )
    call dshr_fldList_add(fldsExport ,'Si_avsdf'    )
    call dshr_fldList_add(fldsExport ,'Si_anidf'    )
    call dshr_fldList_add(fldsExport ,'Faii_swnet'  )
    call dshr_fldList_add(fldsExport ,'Faii_sen'    )
    call dshr_fldList_add(fldsExport ,'Faii_lat'    )
    call dshr_fldList_add(fldsExport ,'Faii_lwup'   )
    call dshr_fldList_add(fldsExport ,'Faii_evap'   )
    call dshr_fldList_add(fldsExport ,'Faii_taux'   )
    call dshr_fldList_add(fldsExport ,'Faii_tauy'   )
    call dshr_fldList_add(fldsExport ,'Fioi_melth'  )
    call dshr_fldList_add(fldsExport ,'Fioi_meltw'  )
    call dshr_fldList_add(fldsExport ,'Fioi_swpen'  )
    call dshr_fldList_add(fldsExport ,'Fioi_taux'   )
    call dshr_fldList_add(fldsExport ,'Fioi_tauy'   )
    call dshr_fldList_add(fldsExport ,'Fioi_salt'   )
    call dshr_fldList_add(fldsExport ,'Fioi_bcpho'  )
    call dshr_fldList_add(fldsExport ,'Fioi_bcphi'  )
    call dshr_fldList_add(fldsExport ,'Fioi_flxdst' )
    if (flds_i2o_per_cat) then
       call dshr_fldList_add(fldsExport, 'Si_ifrac_n'        , ungridded_lbound=1, ungridded_ubound=1)
       call dshr_fldList_add(fldsExport, 'Fioi_swpen_ifrac_n', ungridded_lbound=1, ungridded_ubound=1)
    end if

    fldlist => fldsExport ! the head of the linked list
    do while (associated(fldlist))
       call NUOPC_Advertise(exportState, standardName=fldlist%stdname, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite('(dice_comp_advertise): Fr_ice'//trim(fldList%stdname), ESMF_LOGMSG_INFO)
       fldList => fldList%next
    enddo

    ! Advertise import fields

    call dshr_fldList_add(fldsImport , trim(flds_scalar_name))
    call dshr_fldList_add(fldsImport, 'Faxa_swvdr' )
    call dshr_fldList_add(fldsImport, 'Faxa_swvdf' )
    call dshr_fldList_add(fldsImport, 'Faxa_swndr' )
    call dshr_fldList_add(fldsImport, 'Faxa_swndf' )
    call dshr_fldList_add(fldsImport, 'Fioo_q'     )
    call dshr_fldList_add(fldsImport, 'Sa_z'       )
    call dshr_fldList_add(fldsImport, 'Sa_u'       )
    call dshr_fldList_add(fldsImport, 'Sa_v'       )
    call dshr_fldList_add(fldsImport, 'Sa_ptem'    )
    call dshr_fldList_add(fldsImport, 'Sa_shum'    )
    call dshr_fldList_add(fldsImport, 'Sa_dens'    )
    call dshr_fldList_add(fldsImport, 'Sa_tbot'    )
    call dshr_fldList_add(fldsImport, 'So_s'       )
    call dshr_fldList_add(fldsImport, 'Faxa_bcph'  , ungridded_lbound=1, ungridded_ubound=3)
    call dshr_fldList_add(fldsImport, 'Faxa_ocph'  , ungridded_lbound=1, ungridded_ubound=3)
    call dshr_fldList_add(fldsImport, 'Faxa_dstdry', ungridded_lbound=1, ungridded_ubound=4)
    call dshr_fldList_add(fldsImport, 'Faxa_dstwet', ungridded_lbound=1, ungridded_ubound=4)

    fldlist => fldsImport ! the head of the linked list
    do while (associated(fldlist))
       call NUOPC_Advertise(importState, standardName=fldlist%stdname, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite('(dice_comp_advertise): Fr_ice'//trim(fldList%stdname), ESMF_LOGMSG_INFO)
       fldList => fldList%next
    enddo

  end subroutine dice_comp_advertise

!===============================================================================

  subroutine dice_comp_realize(importState, exportState, rc)

    ! input/output parameters
    type(ESMF_State)       , intent(inout) :: importState
    type(ESMF_State)       , intent(inout) :: exportState
    integer                , intent(out)   :: rc

    ! local variables
    integer                 :: n, lsize, kf
    type(file_desc_t)       :: pioid
    type(var_desc_t)        :: varid
    type(io_desc_t)         :: pio_iodesc
    integer                 :: numOwnedElements   ! number of elements owned by this PET
    type(ESMF_DistGrid)     :: distGrid           ! mesh distGrid
    type(ESMF_Array)        :: elemMaskArray
    integer                 :: rcode
    character(*), parameter :: subName = "(dice_comp_realize) "
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! -------------------------------------
    ! NUOPC_Realize "realizes" a previously advertised field in the importState and exportState
    ! by replacing the advertised fields with the newly created fields of the same name.
    ! -------------------------------------

    call dshr_fldlist_realize( exportState, fldsExport, flds_scalar_name, flds_scalar_num, mesh, &
         subname//':diceExport', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_fldlist_realize( importState, fldsImport, flds_scalar_name, flds_scalar_num, mesh, &
         subname//':diceImport', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! -------------------------------------
    ! Set pointers to exportState fields
    ! -------------------------------------

    ! Initialize dfields with export state data that has corresponding stream field
    call dshr_dfield_add(dfields, sdat, state_fld='Si_ifrac', strm_fld='ifrac', &
         state=exportState, state_ptr=Si_ifrac, logunit=logunit, masterproc=masterproc, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (flds_i2o_per_cat) then
       call dshr_state_getfldptr(exportState, fldname='Si_ifrac_n'        , fldptr2=Si_ifrac_n        , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call dshr_state_getfldptr(exportState, fldname='Fioi_swpen_ifrac_n', fldptr2=Fioi_swpen_ifrac_n, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Set pointers to exportState fields that have no corresponding stream field
    call dshr_state_getfldptr(exportState, fldname='Si_imask'    , fldptr1=Si_imask    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Si_t'        , fldptr1=Si_t        , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Si_tref'     , fldptr1=Si_tref     , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Si_qref'     , fldptr1=Si_qref     , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Si_avsdr'    , fldptr1=Si_avsdr    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Si_anidr'    , fldptr1=Si_anidr    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Si_avsdf'    , fldptr1=Si_avsdf    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Si_anidf'    , fldptr1=Si_anidf    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Faii_swnet'  , fldptr1=Faii_swnet  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Faii_sen'    , fldptr1=Faii_sen    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Faii_lat'    , fldptr1=Faii_lat    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Faii_lwup'   , fldptr1=Faii_lwup   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Faii_evap'   , fldptr1=Faii_evap   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Faii_taux'   , fldptr1=Faii_taux   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Faii_tauy'   , fldptr1=Faii_tauy   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Fioi_melth'  , fldptr1=Fioi_melth  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Fioi_meltw'  , fldptr1=Fioi_meltw  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Fioi_swpen'  , fldptr1=Fioi_swpen  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Fioi_taux'   , fldptr1=Fioi_taux   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Fioi_tauy'   , fldptr1=Fioi_tauy   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Fioi_salt'   , fldptr1=Fioi_salt   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Fioi_bcpho'  , fldptr1=Fioi_bcpho  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Fioi_bcphi'  , fldptr1=Fioi_bcphi  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, fldname='Fioi_flxdst' , fldptr1=Fioi_flxdst , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Set Si_imask (this corresponds to the ocean mask)
    allocate(imask(sdat%lsize))
    if (trim(model_maskfile) /= nullstr) then 
       ! Read in the ocean fraction from the input namelist ocean mask file and assume 'mask' name on domain file
       rcode = pio_openfile(sdat%pio_subsystem, pioid, sdat%io_type, trim(model_maskfile), pio_nowrite)
       call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)
       rcode = pio_inq_varid(pioid, 'mask', varid) 
       call pio_initdecomp(sdat%pio_subsystem, pio_int, (/sdat%nxg, sdat%nyg/), sdat%gindex, pio_iodesc)
       call pio_read_darray(pioid, varid, pio_iodesc, imask, rcode)
       call pio_closefile(pioid)
       call pio_freedecomp(sdat%pio_subsystem, pio_iodesc)
       ! now set the mask as just the real mask
       Si_imask(:) = real(imask(:), kind=r8)
    else
       ! Obtain the ice mask in the ice mesh file
       call ESMF_MeshGet(mesh, numOwnedElements=numOwnedElements, elementdistGrid=distGrid, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       elemMaskArray = ESMF_ArrayCreate(distGrid, imask, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! the following call sets the varues of imask
       call ESMF_MeshGet(mesh, elemMaskArray=elemMaskArray, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! now set the mask as just the real mask
       Si_imask(:) = real(imask(:), kind=r8)
       call ESMF_ArrayDestroy(elemMaskArray, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! -------------------------------------
    ! Set pointers to importState fields
    ! -------------------------------------

    call dshr_state_getfldptr(importState, fldname='Faxa_swvdr'  , fldptr1=Faxa_swvdr  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='Faxa_swvdf'  , fldptr1=Faxa_swvdf  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='Faxa_swndr'  , fldptr1=Faxa_swndr  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='Faxa_swndf'  , fldptr1=Faxa_swndf  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='Faxa_bcph'   , fldptr2=Faxa_bcph   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='Faxa_ocph'   , fldptr2=Faxa_ocph   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='Faxa_dstdry' , fldptr2=Faxa_dstdry , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='Faxa_dstwet' , fldptr2=Faxa_dstwet , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='Fioo_q'      , fldptr1=Fioo_q      , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='Sa_z'        , fldptr1=Sa_z        , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='Sa_u'        , fldptr1=Sa_u        , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='Sa_v'        , fldptr1=Sa_v        , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='Sa_ptem'     , fldptr1=Sa_ptem     , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='Sa_tbot'     , fldptr1=Sa_tbot     , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='Sa_shum'     , fldptr1=Sa_shum     , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='Sa_dens'     , fldptr1=Sa_dens     , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, fldname='So_s'        , fldptr1=So_s        , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! initialize import arrays
    ! On initial call, these are used for the first use in generating the export state
    ! and should have no impact on the solution!!

    Faxa_swvdr(:)    = 0._r8
    Faxa_swvdf(:)    = 0._r8
    Faxa_swndr(:)    = 0._r8
    Faxa_swndf(:)    = 0._r8
    Faxa_bcph(:,:)   = 0._r8
    Faxa_ocph(:,:)   = 0._r8
    Faxa_dstdry(:,:) = 0._r8
    Faxa_dstwet(:,:) = 0._r8
    Fioo_q(:)        = 0._r8
    Sa_z(:)          = 10.0_r8
    Sa_u(:)          = 5.0_r8
    Sa_v(:)          = 5.0_r8
    Sa_ptem(:)       = 260.0_r8
    Sa_tbot(:)       = 260.0_r8
    Sa_shum(:)       = 0.0014_r8
    Sa_dens(:)       = 1.3_r8
    So_s(:)          = 0._r8

    ! -------------------------------------
    ! allocate module arrays that are not part of import and export state
    ! -------------------------------------

    lsize = size(Si_ifrac)
    allocate(yc(lsize))
    allocate(water(lsize))
    allocate(tfreeze(lsize))

  end subroutine dice_comp_realize

!===============================================================================

  subroutine dice_comp_run(target_ymd, target_tod, cosarg, rc)

    ! --------------------------
    ! advance dice
    ! --------------------------

    ! input/output variables:
    integer  , intent(in)    :: target_ymd ! model date
    integer  , intent(in)    :: target_tod ! model sec into model date
    real(R8) , intent(in)    :: cosarg     ! for setting ice temp pattern
    integer  , intent(out)   :: rc

    ! local variables
    integer           :: n, lsize
    real(r8)          :: jday, jday0        ! elapsed day counters
    integer           :: spatialDim         ! number of dimension in mesh
    integer           :: numOwnedElements   ! size of mesh
    real(r8), pointer :: ownedElemCoords(:) ! mesh lat and lons
    real(r8)          :: qmeltall          ! q that would melt all accumulated water
    logical           :: first_time = .true.
    character(*), parameter :: subName = "(dice_comp_run) "
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('DICE_RUN')

    !--------------------
    ! advance dice streams
    !--------------------

    ! time and spatially interpolate to model time and grid
    call t_barrierf('dice_BARRIER',mpicom)
    call t_startf('dice_strdata_advance')
    call shr_strdata_advance(sdat, target_ymd, target_tod, mpicom, 'dice')
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

    lsize = size(Si_ifrac)

    call t_startf('dice_datamode')
    select case (trim(datamode))

    case('COPYALL')
       ! do nothing extra

    case('SSTDATA')
       if (first_time) then
          if (.not. read_restart) then
             do n = 1,lsize
                if (Si_ifrac(n) > 0.0_r8) then
                   water(n) = flux_Qacc0
                else
                   water(n) = 0.0_r8
                end if
             end do
             ! iFrac0 = iFrac  ! previous step's ice fraction
          endif

          call ESMF_MeshGet(mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          allocate(ownedElemCoords(spatialDim*numOwnedElements))
          call ESMF_MeshGet(mesh, ownedElemCoords=ownedElemCoords)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          do n = 1,numOwnedElements
             yc(n) = ownedElemCoords(2*n)
          end do
          deallocate(ownedElemCoords)
       end if

       tfreeze(:) = shr_frz_freezetemp(So_s(:)) + tFrz ! convert to Kelvin

       do n = 1,lsize

          !--- fix erroneous iFrac ---
          Si_ifrac(n) = min(1.0_r8,max(0.0_r8,Si_ifrac(n)))

          !--- fabricate ice surface T, fix erroneous iFrac ---
          if ( yc(n) > 0.0_r8) then
             Si_t(n) = 260.0_r8 + 10.0_r8*cos(cosArg)
          else
             Si_t(n) = 260.0_r8 - 10.0_r8*cos(cosArg)
          end if

          !--- set albedos (constant) ---
          Si_avsdr(n) = ax_vsdr
          Si_anidr(n) = ax_nidr
          Si_avsdf(n) = ax_vsdf
          Si_anidf(n) = ax_nidf

          !--- swnet is sent to cpl as a diagnostic quantity only ---
          !--- newly recv'd swdn goes with previously sent albedo ---
          !--- but albedos are (currently) time invariant         ---
          Faii_swnet(n)   = (1.0_r8 - Si_avsdr(n))*Faxa_swvdr(n) &
                          + (1.0_r8 - Si_anidr(n))*Faxa_swndr(n) &
                          + (1.0_r8 - Si_avsdf(n))*Faxa_swvdf(n) &
                          + (1.0_r8 - Si_anidf(n))*Faxa_swndf(n)

          !--- compute melt/freeze water balance, adjust iFrac  -------------
          if ( .not. flux_Qacc ) then ! Q accumulation option is OFF
             Fioi_melth(n) = min(Fioo_q(n),0.0_r8 )          ! q<0 => melt potential
             Fioi_melth(n) = max(Fioi_melth(n),Flux_Qmin   ) ! limit the melt rate
             Fioi_meltw(n) =    -Fioi_melth(n)/latice        ! corresponding water flux

          else                                 ! Q accumulation option is ON
             !--------------------------------------------------------------
             ! 1a) Q<0 & iFrac > 0  =>  infinite supply of water to melt
             ! 1b) Q<0 & iFrac = 0  =>  melt accumulated water only
             ! 2a) Q>0 & iFrac > 0  =>  zero-out accumulated water
             ! 2b) Q>0 & iFrac = 0  =>  accumulated water
             !--------------------------------------------------------------

             if ( Fioo_q(n) <  0.0_r8 ) then ! Q<0 => melt
                if (Si_ifrac(n) > 0.0_r8 ) then
                   Fioi_melth(n) = Si_ifrac(n)*max(Fioo_q(n),Flux_Qmin)
                   Fioi_meltw(n) =    -Fioi_melth(n)/latice
                   !  water(n) = < don't change this value >
                else
                   Qmeltall   = -water(n)*latice/dt
                   Fioi_melth(n) = max(Fioo_q(n), Qmeltall, Flux_Qmin )
                   Fioi_meltw(n) = -Fioi_melth(n)/latice
                   water(n) =  water(n) - Fioi_meltw(n)*dt
                end if
             else                       ! Q>0 => freeze
                if (Si_ifrac(n) > 0.0_r8 ) then
                   Fioi_melth(n) = 0.0_r8
                   Fioi_meltw(n) = 0.0_r8
                   water(n) = 0.0_r8
                else
                   Fioi_melth(n) = 0.0_r8
                   Fioi_meltw(n) = 0.0_r8
                   water(n) = water(n) + dt*Fioo_q(n)/latice
                end if
             end if

             if (water(n) < 1.0e-16_r8 ) water(n) = 0.0_r8

             !--- non-zero water => non-zero iFrac ---
             if (Si_ifrac(n) <= 0.0_r8  .and.  water(n) > 0.0_r8) then
                Si_ifrac(n) = min(1.0_r8,water(n)/waterMax)
                ! Si_t(n) = tfreeze(n)     ! T can be above freezing?!?
             end if

             !--- cpl multiplies Fioi_melth & Fioi_meltw by iFrac ---
             !--- divide by iFrac here => fixed quantity flux (not per area) ---
             if (Si_ifrac(n) > 0.0_r8) then
                Si_ifrac(n) = max( 0.01_r8, Si_ifrac(n)) ! min iFrac
                Fioi_melth(n) = Fioi_melth(n)/Si_ifrac(n)
                Fioi_meltw(n) = Fioi_meltw(n)/Si_ifrac(n)
             else
                Fioi_melth(n) = 0.0_r8
                Fioi_meltw(n) = 0.0_r8
             end if
          end if

          !--- modify T wrt iFrac: (iFrac -> 0) => (T -> tfreeze) ---
          Si_t(n) = tfreeze(n) + Si_ifrac(n)*(Si_t(n)-tfreeze(n))
       end do

       ! compute ice/ice surface fluxes
       call dice_flux_atmice( &
            imask     ,Sa_z      ,Sa_u      ,Sa_v     , &
            Sa_ptem   ,Sa_shum   ,Sa_dens   ,Sa_tbot  , &
            Si_t      ,Faii_sen  ,Faii_lat  ,Faii_lwup, &
            Faii_evap ,Faii_taux ,Faii_tauy ,Si_tref  , &
            Si_qref   ,logunit )

       ! compute ice/oce surface fluxes (except melth & meltw, see above)
       do n=1,lsize
          if (imask(n) == 0) then
             Fioi_swpen(n) = spval
             Fioi_melth(n) = spval
             Fioi_meltw(n) = spval
             Fioi_salt (n) = spval
             Fioi_taux(n)  = spval
             Fioi_tauy(n)  = spval
             Si_ifrac(n)   = 0.0_r8
          else
             !--- penetrating short wave ---
             Fioi_swpen(n) = max(0.0_r8, flux_swpf*Faii_swnet(n) ) ! must be non-negative

             !--- i/o surface stress ( = atm/ice stress) ---
             Fioi_taux(n) = Faii_taux(n)
             Fioi_tauy(n) = Faii_tauy(n)

             !--- salt flux ---
             Fioi_salt(n) = 0.0_r8
          end if
          ! !--- save ifrac for next timestep
          ! iFrac0(n) = Si_ifrac(n)
       end do

       ! Compute outgoing aerosol fluxes
       do n = 1,lsize
          Fioi_bcpho(n) = Faxa_bcph(2,n)
          Fioi_bcphi(n) = Faxa_bcph(1,n) + Faxa_bcph(3,n)
          Fioi_flxdst(n) =  Faxa_dstdry(1,n) + Faxa_dstwet(1,n) + Faxa_dstdry(2,n) + Faxa_dstwet(2,n) &
                          + Faxa_dstdry(3,n) + Faxa_dstwet(3,n) + Faxa_dstdry(4,n) + Faxa_dstwet(4,n)
       end do

    end select

    !-------------------------------------------------
    ! optional per thickness category fields
    !-------------------------------------------------

    if (flds_i2o_per_cat) then
       do n = 1,lsize
          Si_iFrac_n(1,n)         = Si_ifrac(n)
          Fioi_swpen_iFrac_n(1,n) = Fioi_swpen(n) * Si_ifrac(n)
       end do
    end if

    first_time = .false.

    call t_stopf('dice_datamode')
    call t_stopf('DICE_RUN')

  end subroutine dice_comp_run

end module ice_comp_nuopc
