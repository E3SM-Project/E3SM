module rof_comp_esmf
  
#ifdef ESMF_INTERFACE
!===============================================================================
! !DESCRIPTION:
! Interface of the active river runoff componetn (RTM) model component of CESM 
! with the main CESM driver. This is a thin interface taking CESM driver information
! in MCT (Model Coupling Toolkit) format and converting it to use by RTM and outputing
! if in ESMF (Earth System Modelling Framework) format.
!
! !USES:
  use esmf
  use esmfshr_util_mod
  use shr_kind_mod     , only : r8 => shr_kind_r8, SHR_KIND_CL
  use shr_file_mod     , only : shr_file_setLogUnit, shr_file_setLogLevel
  use shr_file_mod     , only : shr_file_getUnit, shr_file_setIO
  use shr_const_mod    , only : SHR_CONST_REARTH
  use shr_string_mod   , only : shr_string_listGetNum
  use shr_sys_mod
  use seq_timemgr_mod  , only : seq_timemgr_EClockGetData, seq_timemgr_StopAlarmIsOn
  use seq_timemgr_mod  , only : seq_timemgr_RestartAlarmIsOn, seq_timemgr_EClockDateInSync
  use seq_infodata_mod , only : seq_infodata_start_type_cont
  use seq_infodata_mod , only : seq_infodata_start_type_brnch
  use seq_infodata_mod , only : seq_infodata_start_type_start
  use seq_comm_mct     , only : seq_comm_suffix, seq_comm_inst, seq_comm_name
  use seq_flds_mod
  use RunoffMod        , only : runoff
  use RtmVar           , only : rtmlon, rtmlat, ice_runoff, iulog
  use RtmVar           , only : nsrStartup, nsrContinue, nsrBranch
  use RtmVar           , only : inst_index, inst_suffix, inst_name, RtmVarSet
  use RtmVar           , only : rtm_active, flood_active
  use RtmSpmd          , only : masterproc, iam, RtmSpmdInit
  use RtmMod           , only : Rtmini, Rtmrun
  use RtmTimeManager   , only : timemgr_setup, get_curr_date, get_step_size, advance_timestep 
  use rtm_cpl_indices  , only : rtm_cpl_indices_set, nt_rtm, rtm_tracers
  use rtm_cpl_indices  , only : index_r2x_Forr_rofl, index_r2x_Forr_rofi, index_r2x_Flrr_flood
  use rtm_cpl_indices  , only : index_x2r_Flrl_rofl, index_x2r_Flrl_rofi, index_r2x_Flrr_volr
  use perf_mod         , only : t_startf, t_stopf, t_barrierf
  use rof_import_export
!
  implicit none
  SAVE
  private                              ! By default make data private
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: rof_register_esmf          ! register rtm initial, run, final methods
  public :: rof_init_esmf              ! rtm initialization
  public :: rof_run_esmf               ! rtm run phase
  public :: rof_final_esmf             ! rtm finalization/cleanup
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: rof_distgrid_esmf        ! Distribute river runoff grid
  private :: rof_domain_esmf          ! Set the river runoff model domain information
!
! ! PRIVATE DATA MEMBERS:
  real(r8), pointer :: totrunin(:,:)   ! cell tracer lnd forcing on rtm grid (mm/s)

!===============================================================================
contains
!===============================================================================

  subroutine rof_register_esmf(comp, rc)

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Register the rof initial, run, and final phase methods with ESMF.
    !
    ! !ARGUMENTS:
    type(ESMF_GridComp)  :: comp  ! CLM grid component
    integer, intent(out) :: rc    ! return status
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS
    ! Register the callback routines.

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
         rof_init_esmf, phase=1, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, &
         rof_run_esmf, phase=1, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, &
         rof_final_esmf, phase=1, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

  end subroutine rof_register_esmf

!===============================================================================

  subroutine rof_init_esmf(comp, import_state, export_state, EClock, rc)

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Initialize runoff model (rtm)
    !
    ! !ARGUMENTS:
    type(ESMF_GridComp)  :: comp            ! RTM gridded component
    type(ESMF_State)     :: import_state    ! RTM import state
    type(ESMF_State)     :: export_state    ! RTM export state
    type(ESMF_Clock)     :: EClock          ! ESMF synchronization clock
    integer, intent(out) :: rc              ! Return code
    !
    ! !LOCAL VARIABLES:
    integer                      :: mpicom_rof, mpicom_vm, gsize
    type(ESMF_ArraySpec)         :: arrayspec
    type(ESMF_DistGrid)          :: distgrid
    type(ESMF_Array)             :: dom, x2r, r2x
    type(ESMF_VM)                :: vm
    integer                      :: ROFID	          ! rof identifyer
    integer                      :: lsize                 ! size of attribute vector
    integer                      :: g,i,j,n               ! indices
    logical                      :: exists                ! true if file exists
    integer                      :: nsrest                ! restart type
    integer                      :: ref_ymd               ! reference date (YYYYMMDD)
    integer                      :: ref_tod               ! reference time of day (sec)
    integer                      :: start_ymd             ! start date (YYYYMMDD)
    integer                      :: start_tod             ! start time of day (sec)
    integer                      :: stop_ymd              ! stop date (YYYYMMDD)
    integer                      :: stop_tod              ! stop time of day (sec)
    logical                      :: brnch_retain_casename ! flag if should retain the case name on a branch start type
    integer                      :: lbnum                 ! input to memory diagnostic
    integer                      :: shrlogunit,shrloglev  ! old values for log unit and log level
    integer                      :: begr, endr            ! beg and end per-proc runoff indices
    character(len=SHR_KIND_CL)   :: caseid                ! case identifier name
    character(len=SHR_KIND_CL)   :: ctitle                ! case description title
    character(len=SHR_KIND_CL)   :: starttype             ! start-type (startup, continue, branch, hybrid)
    character(len=SHR_KIND_CL)   :: calendar              ! calendar type name
    character(len=SHR_KIND_CL)   :: hostname              ! hostname of machine running on
    character(len=SHR_KIND_CL)   :: version               ! Model version
    character(len=SHR_KIND_CL)   :: username              ! user running the model
    character(len=32), parameter :: sub = 'rof_init_esmf'
    character(len=*),  parameter :: format = "('("//trim(sub)//") :',A)"
    integer                      :: nfields
    real(R8), pointer            :: fptr(:,:)
    character(ESMF_MAXSTR)       :: convCIM, purpComp
    !-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Determine indices

    call rtm_cpl_indices_set()
 
    ! Duplicate the mpi communicator from the current VM 

    call ESMF_VMGetCurrent(vm, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_VMGet(vm, mpiCommunicator=mpicom_vm, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call MPI_Comm_dup(mpicom_vm, mpicom_rof, rc)
    if(rc /= 0) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Initialize rtm MPI communicator

    call ESMF_AttributeGet(export_state, name="ID", value=ROFID, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call RtmSpmdInit(mpicom_rof)

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','rof_init_esmf:start::',lbnum)
    endif
#endif                      

    ! Initialize io log unit

    inst_name   = seq_comm_name(ROFID)
    inst_index  = seq_comm_inst(ROFID)
    inst_suffix = seq_comm_suffix(ROFID)

    call shr_file_getLogUnit (shrlogunit)
    if (masterproc) then
       inquire(file='rof_modelio.nml'//trim(inst_suffix),exist=exists)
       if (exists) then
          iulog = shr_file_getUnit()
          call shr_file_setIO('rof_modelio.nml'//trim(inst_suffix),iulog)
       end if
       write(iulog,format) "RTM model initialization"
    else
       iulog = shrlogunit
    end if

    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)
    
    ! Initialize RTM

    call seq_timemgr_EClockGetData(EClock,                               &
                                   start_ymd=start_ymd,                  &
                                   start_tod=start_tod, ref_ymd=ref_ymd, &
                                   ref_tod=ref_tod, stop_ymd=stop_ymd,   &
                                   stop_tod=stop_tod,                    &
                                   calendar=calendar )

    call ESMF_AttributeGet(export_state, name="case_name", value=caseid, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="case_desc", value=ctitle, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="brnch_retain_casename", &
         value=brnch_retain_casename, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="start_type", value=starttype, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="model_version", value=version, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="hostname", value=hostname, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeGet(export_state, name="username", value=username, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    
    call timemgr_setup( calendar_in=calendar, &
                        start_ymd_in=start_ymd, start_tod_in=start_tod, &
                        ref_ymd_in=ref_ymd, ref_tod_in=ref_tod, &
                        stop_ymd_in=stop_ymd, stop_tod_in=stop_tod)

    if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
       nsrest = nsrStartup
    else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
       nsrest = nsrContinue
    else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
       nsrest = nsrBranch
    else
       call shr_sys_abort( sub//' ERROR: unknown starttype' )
    end if

    call RtmVarSet(caseid_in=caseid, ctitle_in=ctitle,             &
                   brnch_retain_casename_in=brnch_retain_casename, &
                   nsrest_in=nsrest, version_in=version,           &
                   hostname_in=hostname, username_in=username)

    ! Initialize rtm and determine if rtm will be active

    call Rtmini()

    if ( rtm_active ) then
       begr = runoff%begr
       endr = runoff%endr
       allocate(totrunin(begr:endr,nt_rtm))

       !-----------------------------------------
       ! Initialize distgrid
       !-----------------------------------------
       
       distgrid = rof_distgrid_esmf( gsize)

       call ESMF_AttributeSet(export_state, name="gsize", value=gsize, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       !-----------------------------------------
       !  Set arrayspec for dom, r2x and x2r
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
       call rof_domain_esmf( dom)

       !----------------------------------------- 
       !  Create r2x 
       !-----------------------------------------
       
       ! 1d undistributed index of fields, 2d is packed data

       nfields = shr_string_listGetNum(trim(seq_flds_r2x_fields))
       
       r2x = ESMF_ArrayCreate(distgrid=distgrid, arrayspec=arrayspec, distgridToArrayMap=(/2/), &
            undistLBound=(/1/), undistUBound=(/nfields/), name="d2x", rc=rc)
       if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       
       call ESMF_AttributeSet(r2x, name="mct_names", value=trim(seq_flds_r2x_fields), rc=rc)
       if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       !----------------------------------------- 
       !  Create x2r 
       !-----------------------------------------
       
       nfields = shr_string_listGetNum(trim(seq_flds_x2r_fields))

       x2r = ESMF_ArrayCreate(distgrid=distgrid, arrayspec=arrayspec, distgridToArrayMap=(/2/), &
            undistLBound=(/1/), undistUBound=(/nfields/), name="x2d", rc=rc)
       if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       call ESMF_AttributeSet(x2r, name="mct_names", value=trim(seq_flds_x2r_fields), rc=rc)
       if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       !----------------------------------------- 
       ! Add esmf arrays to import and export state 
       !-----------------------------------------

       call ESMF_StateAdd(export_state, (/dom/), rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       call ESMF_StateAdd(import_state, (/x2r/), rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       call ESMF_StateAdd(export_state, (/r2x/), rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       ! Create river runoff export state 

       call ESMF_ArrayGet(r2x, localDe=0, farrayPtr=fptr, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       call rof_export(fptr)

    endif  ! rtm_active

    call ESMF_AttributeSet(export_state, name="rof_present", value=rtm_active, rc=rc)
    if( rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(export_state, name="rofice_present", value=.false., rc=rc)
    if( rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(export_state, name="rof_prognostic", value=rtm_active, rc=rc)
    if( rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(export_state, name="rof_nx", value=rtmlon, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(export_state, name="rof_ny", value=rtmlat, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(export_state, name="flood_present", value=flood_active, rc=rc)
    if( rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Reset shr logging to original values

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

#if (defined _MEMTRACE)
    if(masterproc) then
       write(iulog,*) TRIM(Sub) // ':end::'
       lbnum=1
       call memmon_dump_fort('memmon.out','rof_int_esmf:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

#ifdef USE_ESMF_METADATA
    convCIM  = "CIM"
    purpComp = "Model Component Simulation Description"

    call ESMF_AttributeAdd(comp,  &
                           convention=convCIM, purpose=purpComp, rc=rc)

    call ESMF_AttributeSet(comp, "ShortName", "RTM", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "LongName", &
                           "River Runoff from RTM", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "Description", &
                  "The RTM component is " // &
                  "the river runoff model used in the CESM1.1.  " // &
                  "More information on the RTM project " // &
                  "and access to previous RTM model versions and " // &
                  "documentation can be found via the RTM Web Page.", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ReleaseDate", "2012", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ModelType", "Runoff", &
                           convention=convCIM, purpose=purpComp, rc=rc)

    !    call ESMF_AttributeSet(comp, "Name", "Sam Levis", &
    !                           convention=convCIM, purpose=purpComp, rc=rc)
    !    call ESMF_AttributeSet(comp, "EmailAddress", &
    !                           "slevis@ucar.edu", &
    !                           convention=convCIM, purpose=purpComp, rc=rc)
    !    call ESMF_AttributeSet(comp, "ResponsiblePartyRole", "contact", &
    !                           convention=convCIM, purpose=purpComp, rc=rc)
#endif

  end subroutine rof_init_esmf

!===============================================================================

  subroutine rof_run_esmf(comp, import_state, export_state, EClock, rc)

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Run rtm model
    !
    ! !ARGUMENTS:
    type(ESMF_GridComp)  :: comp            ! RTM gridded component
    type(ESMF_State)     :: import_state    ! RTM import state
    type(ESMF_State)     :: export_state    ! RTM export state
    type(ESMF_Clock)     :: EClock          ! ESMF synchronization clock
    integer, intent(out) :: rc              ! Return code
    !
    ! !LOCAL VARIABLES:
    type(ESMF_Array)  :: x2r, r2x, dom
    integer :: ymd_sync                   ! Sync date (YYYYMMDD)
    integer :: yr_sync                    ! Sync current year
    integer :: mon_sync                   ! Sync current month
    integer :: day_sync                   ! Sync current day
    integer :: tod_sync                   ! Sync current time of day (sec)
    integer :: ymd                        ! RTM current date (YYYYMMDD)
    integer :: yr                         ! RTM current year
    integer :: mon                        ! RTM current month
    integer :: day                        ! RTM current day
    integer :: tod                        ! RTM current time of day (sec)
    logical :: rstwr_sync                 ! .true. ==> write restart file before returning
    logical :: rstwr                      ! .true. ==> write restart file before returning
    logical :: nlend_sync                 ! Flag signaling last time-step
    logical :: nlend                      ! .true. ==> last time-step
    integer :: shrlogunit,shrloglev       ! old values
    integer :: lbnum                      ! input to memory diagnostic
    integer :: g,i,kf                     ! counters
    real(R8), pointer :: fptr(:,:)
    character(len=32)            :: rdate ! date char string for restart file names
    character(len=32), parameter :: sub = "rof_run_esmf"
    !---------------------------------------------------------------------------

    call t_startf ('lc_rof_run1')

    rc = ESMF_SUCCESS

#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','rof_run_esmf:start::',lbnum)
    endif
#endif

    ! Reset shr logging to my log file

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)

    call seq_timemgr_EClockGetData(EClock, &
         curr_ymd=ymd, curr_tod=tod_sync,  &
         curr_yr=yr_sync, curr_mon=mon_sync, curr_day=day_sync)

    ! Map ESMF to rtm input (rof) data type

    call t_startf ('lc_rof_import')

    call ESMF_StateGet(import_state, itemName="x2d", array=x2r, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_ArrayGet(x2r, localDe=0, farrayPtr=fptr, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call rof_import(fptr, totrunin=totrunin)

    call t_stopf ('lc_rof_import')

    ! Run rtm (input is totrunin, output is runoff%runoff)
    ! First advance rtm time step

    write(rdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr_sync,mon_sync,day_sync,tod_sync
    nlend = seq_timemgr_StopAlarmIsOn( EClock )
    rstwr = seq_timemgr_RestartAlarmIsOn( EClock )
    call advance_timestep()
    call Rtmrun(totrunin, rstwr, nlend, rdate)

    ! Map roff data to MCT datatype (input is runoff%runoff, output is r2x_r)
      
    call t_startf ('lc_rof_export')

    call ESMF_StateGet(export_state, itemName="d2x", array=r2x, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_ArrayGet(r2x, localDe=0, farrayPtr=fptr, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call rof_export(fptr)

    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call t_stopf ('lc_rof_export')
       
    ! Check that internal clock is in sync with master clock
    call get_curr_date( yr, mon, day, tod)
    ymd = yr*10000 + mon*100 + day
    tod = tod
    if ( .not. seq_timemgr_EClockDateInSync( EClock, ymd, tod ) )then
       call seq_timemgr_EclockGetData( EClock, curr_ymd=ymd_sync, curr_tod=tod_sync )
       write(iulog,*)' rtm ymd=',ymd     ,'  rtm tod= ',tod
       write(iulog,*)'sync ymd=',ymd_sync,' sync tod= ',tod_sync
       call shr_sys_abort( sub//":: RTM clock not in sync with Master Sync clock" )
    end if
    
    ! Reset shr logging to my original values

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
  
#if (defined _MEMTRACE)
    if(masterproc) then
       lbnum=1
       call memmon_dump_fort('memmon.out','rof_run_esmf:end::',lbnum)
       call memmon_reset_addr()
    endif
#endif

  end subroutine rof_run_esmf

!===============================================================================

  subroutine rof_final_esmf(comp, import_state, export_state, EClock, rc)

    !------------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Finalize rtm
    !
    ! !ARGUMENTS:
    type(ESMF_GridComp)  :: comp            ! RTM gridded component
    type(ESMF_State)     :: import_state    ! RTM import state
    type(ESMF_State)     :: export_state    ! RTM export state
    type(ESMF_Clock)     :: EClock          ! ESMF synchronization clock
    integer, intent(out) :: rc              ! Return code
    !---------------------------------------------------------------------------
    
    rc = ESMF_SUCCESS
    
    ! Destroy ESMF objects
    call esmfshr_util_StateArrayDestroy(export_state,'domain',rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    
    call esmfshr_util_StateArrayDestroy(export_state,'d2x',rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call esmfshr_util_StateArrayDestroy(import_state,'x2d',rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

  end subroutine rof_final_esmf

!===============================================================================

  function rof_DistGrid_esmf(gsize)

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Set the ESMF distributed grid for the runoff model
    !
    ! RETURN
    implicit none
    type(ESMF_DistGrid)  :: rof_DistGrid_esmf  
    !
    ! ARGUMENTS:
    integer, intent(out) :: gsize            ! Grid size
    !
    ! Local Variables
    integer,allocatable :: gindex(:)         ! indexing for runoff grid cells
    integer :: n, ni                         ! indices
    integer :: lsize                         ! size of runoff data and number of grid cells
    integer :: begr, endr                    ! beg, end runoff indices
    integer :: rc                            ! return code
    character(len=32), parameter :: sub = 'rof_DistGrid_esmf'
    !-------------------------------------------------------------------

    ! Build the rof grid numbering
    ! NOTE:  Numbering scheme is: West to East and South to North
    ! starting at south pole.  Should be the same as what's used in SCRIP
    
    begr  = runoff%begr
    endr  = runoff%endr
    lsize = runoff%lnumr
    gsize = rtmlon*rtmlat

    ! Check 
    ni = 0
    do n = begr,endr
       ni = ni + 1
       if (ni > lsize) then
          write(iulog,*) sub, ' : ERROR runoff count',n,ni,runoff%lnumr
          call shr_sys_abort( sub//' ERROR: runoff > expected' )
       endif
    end do
    if (ni /= lsize) then
       write(iulog,*) sub, ' : ERROR runoff total count',ni,runoff%lnumr
       call shr_sys_abort( sub//' ERROR: runoff not equal to expected' )
    endif

    allocate(gindex(lsize))
    ni = 0
    do n = begr,endr
       ni = ni + 1
       gindex(ni) = runoff%gindex(n)
    end do

    rof_distgrid_esmf = ESMF_DistGridCreate(arbSeqIndexList=gindex, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    deallocate(gindex)

  end function rof_DistGrid_esmf

!===============================================================================

  subroutine rof_domain_esmf( dom )

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    !
    ! Send the runoff model domain information to the coupler
    !
    ! !ARGUMENTS:
    type(ESMF_Array), intent(inout)     :: dom   ! Domain information
    !
    ! Local Variables
    integer                      :: n, ni                          ! index
    integer                      :: klon,klat,karea,kmask,kfrac    ! domain fields
    real(r8)                     :: re = SHR_CONST_REARTH*0.001_r8 ! radius of earth (km)
    real(R8), pointer            :: fptr(:,:)
    integer                      :: rc                             ! Return code
    character(len=32), parameter :: sub = 'rof_domain_esmf'
    !-------------------------------------------------------------------
    !
    ! Initialize domain type
    ! lat/lon in degrees,  area in radians^2
    ! 
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
    ni = 0
    do n = runoff%begr,runoff%endr
       ni = ni + 1
       fptr(klon , ni) = runoff%lonc(n)
       fptr(klat , ni) = runoff%latc(n)
       fptr(karea, ni) = runoff%area(n)*1.0e-6_r8/(re*re)
       fptr(kmask, ni) = 1.0_r8
       fptr(kfrac, ni) = 1.0_r8
    end do

  end subroutine rof_domain_esmf
 
#endif
end module rof_comp_esmf
