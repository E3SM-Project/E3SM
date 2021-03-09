module lnd_comp_esmf

#ifdef ESMF_INTERFACE  
  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  !  Interface of the active land model component of E3SM the ELM (E3SM Land Model)
  !  with the main E3SM driver. This is a thin interface taking E3SM driver information
  !  in MCT (Model Coupling Toolkit) format and converting it to use by CLM and outputing
  !  if in ESMF (Earth System Modelling Framework) format.
  !
  ! !USES:
  use esmf
  use esmfshr_util_mod
  use shr_kind_mod      , only : r8 => shr_kind_r8, SHR_KIND_CL
  use shr_string_mod    , only : shr_string_listGetNum
  use abortutils        , only : endrun
  use domainMod         , only : ldomain
  use decompMod         , only : ldecomp, bounds_type, get_proc_bounds
  use elm_varctl        , only : iulog
  use elm_instMod       , only : lnd2atm_vars, atm2lnd_vars, lnd2glc_vars, glc2lnd_vars
  use elm_cpl_indices
  use lnd_import_export
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  SAVE
  private                              ! By default make data private
  !
  public :: lnd_register_esmf          ! register elm initial, run, final methods
  public :: lnd_init_esmf              ! elm initialization
  public :: lnd_run_esmf               ! elm run phase
  public :: lnd_final_esmf             ! elm finalization/cleanup
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: lnd_distgrid_esmf        ! Distribute elm grid
  private :: lnd_domain_esmf          ! Set the land model domain information
  !---------------------------------------------------------------------------

contains

  !---------------------------------------------------------------------------
  subroutine lnd_register_esmf(comp, rc)
    !
    ! !DESCRIPTION:
    ! Register the elm initial, run, and final phase methods with ESMF.
    !
    ! !ARGUMENTS:
    type(ESMF_GridComp)  :: comp  ! CLM grid component
    integer, intent(out) :: rc    ! return status
    !-----------------------------------------------------------------------
    rc = ESMF_SUCCESS
    ! Register the callback routines.

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
         lnd_init_esmf, phase=1, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, &
         lnd_run_esmf, phase=1, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, &
         lnd_final_esmf, phase=1, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

  end subroutine lnd_register_esmf

  !---------------------------------------------------------------------------
  subroutine lnd_init_esmf(comp, import_state, export_state, EClock, rc)
    !
    ! !DESCRIPTION:
    ! Initialize land surface model and obtain relevant atmospheric model arrays
    ! back from (i.e. albedos, surface temperature and snow cover over land).
    !
    ! !USES:
    use shr_file_mod     , only : shr_file_setLogUnit, shr_file_setLogLevel
    use shr_file_mod     , only : shr_file_getLogUnit, shr_file_getLogLevel
    use shr_file_mod     , only : shr_file_getUnit, shr_file_setIO
    use clm_time_manager , only : get_nstep, get_step_size, set_timemgr_init, set_nextsw_cday
    use elm_initializeMod, only : initialize1, initialize2, initialize3
    use elm_instMod      , only : lnd2atm_vars, lnd2glc_vars
    use elm_varctl       , only : finidat,single_column, elm_varctl_set, noland
    use elm_varctl       , only : inst_index, inst_suffix, inst_name
    use elm_varctl       , only : nsrStartup, nsrContinue, nsrBranch
    use elm_varorb       , only : eccen, obliqr, lambm0, mvelpp
    use controlMod       , only : control_setNL
    use spmdMod          , only : masterproc, spmd_init
    use seq_timemgr_mod  , only : seq_timemgr_EClockGetData
    use seq_infodata_mod , only : seq_infodata_start_type_cont
    use seq_infodata_mod , only : seq_infodata_start_type_brnch
    use seq_infodata_mod , only : seq_infodata_start_type_start
    use seq_comm_mct     , only : seq_comm_suffix, seq_comm_inst, seq_comm_name
    use seq_flds_mod
    !
    ! !ARGUMENTS:
    type(ESMF_GridComp)          :: comp            ! CLM gridded component
    type(ESMF_State)             :: import_state    ! CLM import state
    type(ESMF_State)             :: export_state    ! CLM export state
    type(ESMF_Clock)             :: EClock          ! ESMF synchronization clock
    integer, intent(out)         :: rc              ! Return code
    !
    ! !LOCAL VARIABLES:
    integer                      :: mpicom_lnd, mpicom_vm, gsize
    type(ESMF_ArraySpec)         :: arrayspec
    type(ESMF_DistGrid)          :: distgrid
    type(ESMF_Array)             :: dom, l2x, x2l
    type(ESMF_VM)                :: vm
    integer                      :: lsize                 ! size of attribute vector
    integer                      :: g,i,j                 ! indices
    integer                      :: dtime_sync            ! coupling time-step from the input synchronization clock
    integer                      :: dtime_elm             ! elm time-step
    logical                      :: exists                ! true if file exists
    real(r8)                     :: scmlat                ! single-column latitude
    real(r8)                     :: scmlon                ! single-column longitude
    real(r8)                     :: nextsw_cday           ! calday from clock of next radiation computation
    character(len=SHR_KIND_CL)   :: caseid                ! case identifier name
    character(len=SHR_KIND_CL)   :: ctitle                ! case description title
    character(len=SHR_KIND_CL)   :: starttype             ! start-type (startup, continue, branch, hybrid)
    character(len=SHR_KIND_CL)   :: calendar              ! calendar type name
    character(len=SHR_KIND_CL)   :: hostname              ! hostname of machine running on
    character(len=SHR_KIND_CL)   :: version               ! Model version
    character(len=SHR_KIND_CL)   :: username              ! user running the model
    integer                      :: nsrest                ! elm restart type
    integer                      :: ref_ymd               ! reference date (YYYYMMDD)
    integer                      :: ref_tod               ! reference time of day (sec)
    integer                      :: start_ymd             ! start date (YYYYMMDD)
    integer                      :: start_tod             ! start time of day (sec)
    integer                      :: stop_ymd              ! stop date (YYYYMMDD)
    integer                      :: stop_tod              ! stop time of day (sec)
    logical                      :: brnch_retain_casename ! flag if should retain the case name on a branch start type
    logical                      :: atm_aero              ! Flag if aerosol data sent from atm model
    integer                      :: lbnum                 ! input to memory diagnostic
    integer                      :: shrlogunit,shrloglev  ! old values for log unit and log level
    type(bounds_type)            :: bounds                ! bounds
    integer                      :: LNDID                 ! cesm ID value
    integer                      :: nfields
    real(R8), pointer            :: fptr(:, :)
    character(len=32), parameter :: sub = 'lnd_init_esmf'
    character(len=*),  parameter :: format = "('("//trim(sub)//") :',A)"
    character(ESMF_MAXSTR)       :: convCIM, purpComp
    !-----------------------------------------------------------------------

    ! Determine indices

    call elm_cpl_indices_set()

    rc = ESMF_SUCCESS

    ! duplicate the mpi communicator from the current VM 

    call ESMF_VMGetCurrent(vm, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_VMGet(vm, mpiCommunicator=mpicom_vm, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call MPI_Comm_dup(mpicom_vm, mpicom_lnd, rc)
    if(rc /= 0) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="ID", value=LNDID, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call spmd_init( mpicom_lnd, LNDID )

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_init_esmf:start::',lbnum)
    endif
#endif                      

    inst_name   = seq_comm_name(LNDID)
    inst_index  = seq_comm_inst(LNDID)
    inst_suffix = seq_comm_suffix(LNDID)

    ! Initialize io log unit

    call shr_file_getLogUnit (shrlogunit)
    if (masterproc) then
       inquire(file='lnd_modelio.nml'//trim(inst_suffix),exist=exists)
       if (exists) then
          iulog = shr_file_getUnit()
          call shr_file_setIO('lnd_modelio.nml'//trim(inst_suffix),iulog)
       end if
       write(iulog,format) "CLM land model initialization"
    else
       iulog = shrlogunit
    end if

    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)

    ! Use infodata to set orbital values

    call ESMF_AttributeGet(export_state, name="orb_eccen", value=eccen, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="orb_mvelpp", value=mvelpp, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="orb_lambm0", value=lambm0, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="orb_obliqr", value=obliqr, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Consistency check on namelist filename	

    call control_setNL("lnd_in"//trim(inst_suffix))

    ! Initialize elm
    ! initialize1 reads namelist, grid and surface data
    ! initialize2 performs rest of initialization    

    call seq_timemgr_EClockGetData(EClock,                          &
         start_ymd=start_ymd, start_tod=start_tod, ref_ymd=ref_ymd, &
         ref_tod=ref_tod, stop_ymd=stop_ymd, stop_tod=stop_tod, calendar=calendar )

    call ESMF_AttributeGet(export_state, name="case_name", value=caseid, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="case_desc", value=ctitle, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="single_column", value=single_column, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="scmlat", value=scmlat, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="scmlon", value=scmlon, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="brnch_retain_casename", value=brnch_retain_casename, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="start_type", value=starttype, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="model_version", value=version, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="hostname", value=hostname, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="username", value=username, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call set_timemgr_init( &
         calendar_in=calendar, start_ymd_in=start_ymd, start_tod_in=start_tod, &
         ref_ymd_in=ref_ymd, ref_tod_in=ref_tod, stop_ymd_in=stop_ymd,         &
         stop_tod_in=stop_tod)

    if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
       nsrest = nsrStartup
    else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
       nsrest = nsrContinue
    else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
       nsrest = nsrBranch
    else
       call endrun( sub//' ERROR: unknown starttype' )
    end if

    call elm_varctl_set(                                         &
         caseid_in=caseid, ctitle_in=ctitle,                     &
         brnch_retain_casename_in=brnch_retain_casename,         &
         single_column_in=single_column, scmlat_in=scmlat,       &
         scmlon_in=scmlon, nsrest_in=nsrest, version_in=version, &
         hostname_in=hostname, username_in=username )

    call initialize1( )

    ! If no land then exit out of initialization

    if ( noland) then
       call ESMF_AttributeSet(export_state, name="lnd_present", value=.false., rc=rc)
       if( rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       call ESMF_AttributeSet(export_state, name="lnd_prognostic", value=.false., rc=rc)
       if( rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    end if

    ! Determine if aerosol and dust deposition come from atmosphere component

    rc = ESMF_SUCCESS

    call ESMF_AttributeGet(export_state, name="atm_aero", value=atm_aero, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    if ( .not. atm_aero )then
       call endrun( sub//' ERROR: atmosphere model MUST send aerosols to CLM' )
    end if

    !-----------------------------------------
    ! Initialize distgrid
    !-----------------------------------------

    call get_proc_bounds(bounds)

    distgrid = lnd_distgrid_esmf(bounds, gsize)

    call ESMF_AttributeSet(export_state, name="gsize", value=gsize, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !-----------------------------------------
    !  Set arrayspec for dom, l2x and x2l
    !-----------------------------------------
    
    call ESMF_ArraySpecSet(arrayspec, rank=2, typekind=ESMF_TYPEKIND_R8, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !-----------------------------------------
    ! Create dom 
    !-----------------------------------------

    nfields = shr_string_listGetNum(trim(seq_flds_dom_fields))

    dom = ESMF_ArrayCreate(distgrid=distgrid, arrayspec=arrayspec, distgridToArrayMap=(/2/), &
         undistLBound=(/1/), undistUBound=(/nfields/), name="domain", rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(dom, name="mct_names", value=trim(seq_flds_dom_fields), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Set values of dom
    call lnd_domain_esmf(bounds, dom)

    !----------------------------------------- 
    !  Create l2x 
    !-----------------------------------------

    ! 1d undistributed index of fields, 2d is packed data

    nfields = shr_string_listGetNum(trim(seq_flds_l2x_fields))

    l2x = ESMF_ArrayCreate(distgrid=distgrid, arrayspec=arrayspec, distgridToArrayMap=(/2/), &
         undistLBound=(/1/), undistUBound=(/nfields/), name="d2x", rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(l2x, name="mct_names", value=trim(seq_flds_l2x_fields), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !----------------------------------------- 
    !  Create x2l 
    !-----------------------------------------

    nfields = shr_string_listGetNum(trim(seq_flds_x2l_fields))

    x2l = ESMF_ArrayCreate(distgrid=distgrid, arrayspec=arrayspec, distgridToArrayMap=(/2/), &
         undistLBound=(/1/), undistUBound=(/nfields/), name="x2d", rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(x2l, name="mct_names", value=trim(seq_flds_x2l_fields), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !----------------------------------------- 
    ! Add esmf arrays to import and export state 
    !-----------------------------------------

    call ESMF_StateAdd(export_state, (/dom/), rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_StateAdd(export_state, (/l2x/), rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_StateAdd(import_state, (/x2l/), rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Finish initializing elm

    call initialize2()
    call initialize3()

    ! Check that elm internal dtime aligns with elm coupling interval

    call seq_timemgr_EClockGetData(EClock, dtime=dtime_sync )
    dtime_elm = get_step_size()
    if(masterproc) write(iulog,*)'dtime_sync= ',dtime_sync,&
         ' dtime_elm= ',dtime_elm,' mod = ',mod(dtime_sync,dtime_elm)
    if (mod(dtime_sync,dtime_elm) /= 0) then
       write(iulog,*)'elm dtime ',dtime_elm,' and Eclock dtime ',dtime_sync,' never align'
       call endrun( sub//' ERROR: time out of sync' )
    end if

    ! Create land export state 

    call ESMF_ArrayGet(l2x, localDe=0, farrayPtr=fptr, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call lnd_export(bounds, lnd2atm_vars, lnd2glc_vars, fptr)

    ! Set land modes

    call ESMF_AttributeSet(export_state, name="lnd_prognostic", value=.true., rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(export_state, name="lnd_nx", value=ldomain%ni, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(export_state, name="lnd_ny", value=ldomain%nj, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Get nextsw_cday

    call ESMF_AttributeGet(export_state, name="nextsw_cday", value=nextsw_cday, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call set_nextsw_cday( nextsw_cday )

    ! Reset shr logging to original values

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

#if (defined _MEMTRACE)
    if(masterproc) then
       write(iulog,*) TRIM(Sub) // ':end::'
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_int_esmf:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

#ifdef USE_ESMF_METADATA
    convCIM  = "CIM"
    purpComp = "Model Component Simulation Description"

    call ESMF_AttributeAdd(comp,  &
         convention=convCIM, purpose=purpComp, rc=rc)

    call ESMF_AttributeSet(comp, "ShortName", "CLM", &
         convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "LongName", &
         "Community Land Model", &
         convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "Description", &
         "The Community Land Model version 4.0 is " // &
         "the land model used in the CESM1.0.  " // &
         "More information on the CLM project " // &
         "and access to previous CLM model versions and " // &
         "documentation can be found via the CLM Web Page.", &
         convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ReleaseDate", "2010", &
         convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ModelType", "Land", &
         convention=convCIM, purpose=purpComp, rc=rc)
#endif

  end subroutine lnd_init_esmf

  !---------------------------------------------------------------------------
  subroutine lnd_run_esmf(comp, import_state, export_state, EClock, rc)
    !
    ! !DESCRIPTION:
    ! Run elm model
    !
    ! !USES:
    use shr_file_mod      , only : shr_file_setLogUnit, shr_file_setLogLevel
    use shr_file_mod      , only : shr_file_getLogUnit, shr_file_getLogLevel
    use shr_orb_mod       , only : shr_orb_decl
    use elm_instMod       , only : lnd2atm_vars, atm2lnd_vars, lnd2glc_vars, glc2lnd_vars
    use elm_driver        , only : elm_drv
    use elm_varorb        , only : eccen, obliqr, lambm0, mvelpp
    use clm_time_manager  , only : get_curr_date, get_nstep, get_curr_calday, get_step_size
    use clm_time_manager  , only : advance_timestep, set_nextsw_cday,update_rad_dtime
    use seq_timemgr_mod   , only : seq_timemgr_EClockGetData, seq_timemgr_StopAlarmIsOn
    use seq_timemgr_mod   , only : seq_timemgr_RestartAlarmIsOn, seq_timemgr_EClockDateInSync
    use spmdMod           , only : masterproc, mpicom
    use perf_mod          , only : t_startf, t_stopf, t_barrierf
    !
    ! !ARGUMENTS:
    type(ESMF_GridComp)  :: comp            ! CLM gridded component
    type(ESMF_State)     :: import_state    ! CLM import state
    type(ESMF_State)     :: export_state    ! CLM export state
    type(ESMF_Clock)     :: EClock          ! ESMF synchronization clock
    integer, intent(out) :: rc              ! Return code
    !
    ! !LOCAL VARIABLES:
    type(ESMF_Array)  :: l2x, x2l, dom
    real(R8), pointer :: fptr(:, :)
    integer :: ymd_sync                   ! Sync date (YYYYMMDD)
    integer :: yr_sync                    ! Sync current year
    integer :: mon_sync                   ! Sync current month
    integer :: day_sync                   ! Sync current day
    integer :: tod_sync                   ! Sync current time of day (sec)
    integer :: ymd                        ! CLM current date (YYYYMMDD)
    integer :: yr                         ! CLM current year
    integer :: mon                        ! CLM current month
    integer :: day                        ! CLM current day
    integer :: tod                        ! CLM current time of day (sec)
    integer :: dtime                      ! time step increment (sec)
    integer :: nstep                      ! time step index
    logical :: rstwr_sync                 ! .true. ==> write restart file before returning
    logical :: rstwr                      ! .true. ==> write restart file before returning
    logical :: nlend_sync                 ! Flag signaling last time-step
    logical :: nlend                      ! .true. ==> last time-step
    logical :: dosend                     ! true => send data back to driver
    logical :: doalb                      ! .true. ==> do albedo calculation on this time step
    real(r8):: nextsw_cday                ! calday from clock of next radiation computation
    real(r8):: caldayp1                   ! elm calday plus dtime offset
    real(r8):: calday                     ! calendar day for nstep
    real(r8):: declin                     ! solar declination angle in radians for nstep
    real(r8):: declinp1                   ! solar declination angle in radians for nstep+1
    real(r8):: eccf                       ! earth orbit eccentricity factor
    integer :: shrlogunit,shrloglev       ! old values
    integer :: lbnum                      ! input to memory diagnostic
    integer :: g,i,ka                     ! counters
    real(r8):: recip                      ! recip
    logical :: glcrun_alarm               ! if true, sno data is averaged and sent to glc this step
    type(bounds_type) :: bounds           ! bounds
    logical,save :: first_call = .true.   ! first call work
    character(len=32)            :: rdate ! date char string for restart file names
    character(len=32), parameter :: sub = "lnd_run_esmf"
    !---------------------------------------------------------------------------

    call get_proc_bounds(bounds)

    call t_startf ('lc_lnd_run1')
    rc = ESMF_SUCCESS

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_run_esmf:start::',lbnum)
    endif
#endif

    ! Reset shr logging to my log file

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)

    ! Determine time of next atmospheric shortwave calculation

    call seq_timemgr_EClockGetData(EClock, &
         curr_ymd=ymd, curr_tod=tod_sync,  &
         curr_yr=yr_sync, curr_mon=mon_sync, curr_day=day_sync)
    call ESMF_AttributeGet(export_state, name="nextsw_cday", value=nextsw_cday, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call set_nextsw_cday( nextsw_cday )
    dtime = get_step_size()

    write(rdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr_sync,mon_sync,day_sync,tod_sync
    nlend_sync = seq_timemgr_StopAlarmIsOn( EClock )
    rstwr_sync = seq_timemgr_RestartAlarmIsOn( EClock )

    call t_stopf ('lc_lnd_run1')

    ! Map ESMF to CLM data type

    call t_startf ('lc_lnd_import')

    call ESMF_StateGet(import_state, itemName="x2d", array=x2l, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_ArrayGet(x2l, localDe=0, farrayPtr=fptr, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call lnd_import( bounds, fptr, atm2lnd_vars, glc2lnd_vars )

    call t_stopf ('lc_lnd_import')

    ! Use infodata to set orbital values if it was updated at run time

    call ESMF_AttributeGet(export_state, name="orb_eccen", value=eccen, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="orb_mvelpp", value=mvelpp, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="orb_lambm0", value=lambm0, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="orb_obliqr", value=obliqr, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Loop over time steps in coupling interval

    call t_startf ('lc_lnd_run2')

    dosend = .false.
    do while(.not. dosend)

       ! Determine if dosend
       ! When time is not updated at the beginning of the loop - then return only if
       ! are in sync with clock before time is updated

       call get_curr_date( yr, mon, day, tod )
       ymd = yr*10000 + mon*100 + day
       tod = tod
       dosend = (seq_timemgr_EClockDateInSync( EClock, ymd, tod))

       ! Determine doalb based on nextsw_cday sent from atm model

       nstep = get_nstep()
       caldayp1 = get_curr_calday(offset=dtime)
       doalb = abs(nextsw_cday- caldayp1) < 1.e-10_r8
       if (nstep == 0) then
          doalb = .false. 	
       else if (nstep == 1) then 
          doalb = (abs(nextsw_cday- caldayp1) < 1.e-10_r8) 
       else
          doalb = (nextsw_cday >= -0.5_r8) 
       end if
       call update_rad_dtime(doalb)

       ! Determine if time to write cam restart and stop

       rstwr = .false.
       if (rstwr_sync .and. dosend) rstwr = .true.
       nlend = .false.
       if (nlend_sync .and. dosend) nlend = .true.

       ! Run elm 

       call t_barrierf('sync_elm_run', mpicom)
       call t_startf ('elm_run')
       calday = get_curr_calday()
       call shr_orb_decl( calday     , eccen, mvelpp, lambm0, obliqr, declin  , eccf )
       call shr_orb_decl( nextsw_cday, eccen, mvelpp, lambm0, obliqr, declinp1, eccf )
       call elm_drv(doalb, nextsw_cday, declinp1, declin, rstwr, nlend, rdate)
       call t_stopf ('elm_run')

       ! Map CLM data type to MCT
       ! Reset landfrac on atmosphere grid to have the right domain

       call t_startf ('lc_lnd_export')
       call ESMF_StateGet(export_state, itemName="d2x", array=l2x, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       call ESMF_ArrayGet(l2x, localDe=0, farrayPtr=fptr, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       call lnd_export(bounds, lnd2atm_vars, lnd2glc_vars, fptr)
       call t_stopf ('lc_lnd_export')

       ! Advance elm time step

       call t_startf ('lc_elm2_adv_timestep')
       call advance_timestep()
       call t_stopf ('lc_elm2_adv_timestep')

    end do

    call t_stopf ('lc_lnd_run2')
    call t_startf('lc_lnd_run3')

    ! Check that internal clock is in sync with master clock

    call get_curr_date( yr, mon, day, tod, offset=-dtime )
    ymd = yr*10000 + mon*100 + day
    tod = tod
    if ( .not. seq_timemgr_EClockDateInSync( EClock, ymd, tod ) )then
       call seq_timemgr_EclockGetData( EClock, curr_ymd=ymd_sync, curr_tod=tod_sync )
       write(iulog,*)' elm ymd=',ymd     ,'  elm tod= ',tod
       write(iulog,*)'sync ymd=',ymd_sync,' sync tod= ',tod_sync
       call endrun( sub//":: CLM clock not in sync with Master Sync clock" )
    end if

    ! Reset shr logging to my original values

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','lnd_run_esmf:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

    first_call = .false.
    call t_stopf ('lc_lnd_run3')

  end subroutine lnd_run_esmf

  !---------------------------------------------------------------------------

  subroutine lnd_final_esmf(comp, import_state, export_state, EClock, rc)
    !
    ! !DESCRIPTION:
    ! Finalize land surface model
    !
    ! !USES:
    use elm_finalizeMod, only : final
    !
    ! !ARGUMENTS:
    type(ESMF_GridComp)  :: comp            ! CLM gridded component
    type(ESMF_State)     :: import_state    ! CLM import state
    type(ESMF_State)     :: export_state    ! CLM export state
    type(ESMF_Clock)     :: EClock          ! ESMF synchronization clock
    integer, intent(out) :: rc              ! Return code
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call final()

    ! Destroy ESMF objects
    call esmfshr_util_StateArrayDestroy(export_state,'domain',rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call esmfshr_util_StateArrayDestroy(export_state,'d2x',rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call esmfshr_util_StateArrayDestroy(import_state,'x2d',rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

  end subroutine lnd_final_esmf

  !---------------------------------------------------------------------------
  function lnd_DistGrid_esmf(bounds, gsize)
    !
    ! !DESCRIPTION:
    ! Setup distributed grid for CLM
    !
    ! !USES:
    use shr_kind_mod , only : r8 => shr_kind_r8
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds ! bounds
    integer, intent(out)          :: gsize  ! grid size
    !
    ! RETURN:
    type(ESMF_DistGrid) :: lnd_DistGrid_esmf  ! Resulting distributed grid
    !
    ! !LOCAL VARIABLES:
    integer,allocatable :: gindex(:)  ! grid indices
    integer :: n                      ! indices
    integer :: rc                     ! error code
    !---------------------------------------------------------------------------

    ! number the local grid

    allocate(gindex(bounds%begg:bounds%endg))
    do n = bounds%begg, bounds%endg
       gindex(n) = ldecomp%gdc2glo(n)
    end do
    gsize = ldomain%ni * ldomain%nj

    lnd_distgrid_esmf = ESMF_DistGridCreate(arbSeqIndexList=gindex, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    deallocate(gindex)

  end function lnd_DistGrid_esmf

  !---------------------------------------------------------------------------
  subroutine lnd_domain_esmf( bounds, dom)
    !
    ! !DESCRIPTION:
    ! Send the land model domain information to the coupler
    !
    ! !USES:
    use elm_varcon  , only : re
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)    :: bounds ! bounds
    type(ESMF_Array)  , intent(inout) :: dom    ! CLM domain data
    !
    ! !LOCAL VARIABLES:
    integer           :: g,i                         ! index
    integer           :: klon,klat,karea,kmask,kfrac ! domain fields
    real(R8), pointer :: fptr (:,:)
    integer           :: rc     ! return code
    !---------------------------------------------------------------------------

    ! Initialize domain type
    ! lat/lon in degrees,  area in radians^2, mask is 1 (land), 0 (non-land)
    ! Note that in addition land carries around landfrac for the purposes of domain checking

    call ESMF_ArrayGet(dom, localDe=0, farrayPtr=fptr, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Fill in correct values for domain components
    klon  = esmfshr_util_ArrayGetIndex(dom,'lon ',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    klat  = esmfshr_util_ArrayGetIndex(dom,'lat ',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    karea = esmfshr_util_ArrayGetIndex(dom,'area',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    kmask = esmfshr_util_ArrayGetIndex(dom,'mask',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    kfrac = esmfshr_util_ArrayGetIndex(dom,'frac',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Fill in correct values for domain components
    ! Note aream will be filled in in the atm-lnd mapper

    fptr(:,:) = -9999.0_R8
    fptr(kmask,:) = -0.0_R8
    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       fptr(klon, i)  = ldomain%lonc(g)
       fptr(klat, i)  = ldomain%latc(g)
       fptr(karea, i) = ldomain%area(g)/(re*re)
       fptr(kmask, i) = real(ldomain%mask(g), r8)
       fptr(kfrac, i) = real(ldomain%frac(g), r8)
    end do

  end subroutine lnd_domain_esmf
  !---------------------------------------------------------------------------

#endif

end module lnd_comp_esmf
