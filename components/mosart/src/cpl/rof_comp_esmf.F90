module rof_comp_esmf
  
#ifdef ESMF_INTERFACE
!===============================================================================
! Interface of the active river runoff componetn (mosart) model component of CESM 
! with the main CESM driver. This is a thin interface taking CESM driver information
! in MCT (Model Coupling Toolkit) format and converting it to use by RTM and outputing
! if in ESMF (Earth System Modelling Framework) format.
!
! !DESCRIPTION:
!
! !USES:
  use shr_kind_mod     , only : r8 => shr_kind_r8, SHR_KIND_CL
  use shr_file_mod     , only : shr_file_setLogUnit, shr_file_setLogLevel, &
                                shr_file_getLogUnit, shr_file_getLogLevel, &
                                shr_file_getUnit, shr_file_setIO
  use shr_const_mod    , only : SHR_CONST_REARTH
  use shr_sys_mod
  use seq_timemgr_mod  , only : seq_timemgr_EClockGetData, seq_timemgr_StopAlarmIsOn, &
                                seq_timemgr_RestartAlarmIsOn, seq_timemgr_EClockDateInSync
  use seq_infodata_mod , only : seq_infodata_start_type_cont, &
                                seq_infodata_start_type_brnch, &
                                seq_infodata_start_type_start
  use seq_comm_mct     , only : seq_comm_suffix, seq_comm_inst, seq_comm_name
  use seq_flds_mod
  use esmf
  use esmfshr_mod
  use RunoffMod        , only : rtmCTL, TRunoff
  use RtmVar           , only : rtmlon, rtmlat, ice_runoff, iulog, &
                                nsrStartup, nsrContinue, nsrBranch, & 
                                inst_index, inst_suffix, inst_name, RtmVarSet
  use RtmSpmd          , only : masterproc, iam, npes, RtmSpmdInit, ROFID
  use RtmMod           , only : Rtmini, Rtmrun
  use RtmTimeManager   , only : timemgr_setup, get_curr_date, get_step_size, advance_timestep 
  use rof_cpl_indices  , only : rof_cpl_indices_set, nt_rtm, rtm_tracers, &
                                index_r2x_Forr_rofl, index_r2x_Forr_rofi, &
                                index_x2r_Flrl_rofi, index_x2r_Flrl_rofsur, &
                                index_x2r_Flrl_rofgwl, index_x2r_Flrl_rofsub, &
                                index_x2r_Flrl_rofdto, &
                                index_r2x_Flrr_flood, &
                                index_r2x_Flrr_volr, index_r2x_Flrr_volrmch
  use perf_mod         , only : t_startf, t_stopf, t_barrierf
!
! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  SAVE
  private                              ! By default make data private
!
! !PUBLIC MEMBER FUNCTIONS:
!
  public :: rof_register_esmf          ! register rtm initial, run, final methods
  public :: rof_init_esmf              ! rtm initialization
  public :: rof_run_esmf               ! rtm run phase
  public :: rof_final_esmf             ! rtm finalization/cleanup
!
! ! PUBLIC DATA MEMBERS: None
!
! !REVISION HISTORY:
! Author:  Mariana Vertenstein
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: rof_DistGrid_esmf        ! Distribute river runoff grid
  private :: rof_domain_esmf          ! Set the river runoff model domain information
  private :: rof_import_esmf          ! Import data from the CESM coupler to the river runoff model
  private :: rof_export_esmf          ! Export the river runoff model data to the CESM coupler
!
! ! PRIVATE DATA MEMBERS:
!

!===============================================================================
contains
!===============================================================================

  subroutine rof_register_esmf(comp, rc)

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Register the rof initial, run, and final phase methods with ESMF.
    !
    ! !ARGUMENTS:
    implicit none
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
    ! Initialize runoff model (mosart)
    !
    ! !ARGUMENTS:
    implicit none
    type(ESMF_GridComp)  :: comp            ! MOSART gridded component
    type(ESMF_State)     :: import_state    ! MOSART import state
    type(ESMF_State)     :: export_state    ! MOSART export state
    type(ESMF_Clock)     :: EClock          ! ESMF synchronization clock
    integer, intent(out) :: rc              ! Return code
    !
    ! !LOCAL VARIABLES:
    logical              :: rof_prognostic           ! flag
    logical              :: flood_present            ! flag
    integer              :: mpicom_rof, mpicom_vm, gsize
    type(ESMF_DistGrid)  :: distgrid
    type(ESMF_Array)     :: dom, x2r, r2x
    type(ESMF_VM)        :: vm
    integer :: lsize                                 ! size of attribute vector
    integer :: g,i,j,n                               ! indices
    logical :: exists                                ! true if file exists
    integer :: nsrest                                ! restart type
    integer :: ref_ymd                               ! reference date (YYYYMMDD)
    integer :: ref_tod                               ! reference time of day (sec)
    integer :: start_ymd                             ! start date (YYYYMMDD)
    integer :: start_tod                             ! start time of day (sec)
    integer :: stop_ymd                              ! stop date (YYYYMMDD)
    integer :: stop_tod                              ! stop time of day (sec)
    logical :: brnch_retain_casename                 ! flag if should retain the case name on a branch start type
    integer :: lbnum                                 ! input to memory diagnostic
    integer :: shrlogunit,shrloglev                  ! old values for log unit and log level
    integer :: begr, endr                            ! beg and end per-proc runoff indices
    character(len=SHR_KIND_CL) :: caseid             ! case identifier name
    character(len=SHR_KIND_CL) :: ctitle             ! case description title
    character(len=SHR_KIND_CL) :: starttype          ! start-type (startup, continue, branch, hybrid)
    character(len=SHR_KIND_CL) :: calendar           ! calendar type name
    character(len=SHR_KIND_CL) :: hostname           ! hostname of machine running on
    character(len=SHR_KIND_CL) :: version            ! Model version
    character(len=SHR_KIND_CL) :: username           ! user running the model
    character(len=32), parameter :: sub = 'rof_init_esmf'
    character(len=*),  parameter :: format = "('("//trim(sub)//") :',A)"
    character(ESMF_MAXSTR) :: convCIM, purpComp
    !-----------------------------------------------------------------------

    ! Determine indices

    call rof_cpl_indices_set()

    rc = ESMF_SUCCESS
 
    ! Duplicate the mpi communicator from the current VM 

    call ESMF_VMGetCurrent(vm, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_VMGet(vm, mpiCommunicator=mpicom_vm, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call MPI_Comm_dup(mpicom_vm, mpicom_rof, rc)
    if(rc /= 0) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Initialize mosart MPI communicator

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
       write(iulog,format) "MOSART model initialization"
    else
       iulog = shrlogunit
    end if

    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)

    if (masterproc) then
       write(iulog,*) ' mosart npes    = ',npes
       write(iulog,*) ' mosart iam     = ',iam
       write(iulog,*) ' mosart ROFID   = ',ROFID
       write(iulog,*) ' mosart name    = ',trim(inst_name)
       write(iulog,*) ' mosart inst    = ',inst_index
       write(iulog,*) ' mosart suffix  = ',trim(inst_suffix)
    endif

    ! Initialize MOSART

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

    ! Initialize mosart and determine if mosart will be active

    call Rtmini(rtm_active=rof_prognostic,flood_active=flood_present)

    if ( rof_prognostic ) then
       begr = rtmCTL%begr
       endr = rtmCTL%endr

       ! Initialize rof distgrid and domain

       distgrid = rof_DistGrid_esmf(gsize, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       call ESMF_AttributeSet(export_state, name="gsize", value=gsize, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    
       dom = mct2esmf_init(distgrid, attname=seq_flds_dom_fields, name="domain", rc=rc)
       if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    
       call rof_domain_esmf( dom, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       ! Initialize rof import and export states

       r2x = mct2esmf_init(distgrid, attname=seq_flds_r2x_fields, name="d2x", rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       x2r = mct2esmf_init(distgrid, attname=seq_flds_x2r_fields, name="x2d", rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       call ESMF_StateAdd(export_state, (/dom/), rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       call ESMF_StateAdd(import_state, (/x2r/), rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       call ESMF_StateAdd(export_state, (/r2x/), rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       ! Create river runoff export state 

       call rof_export_esmf(r2x, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    endif  ! rof_prognostic

    call ESMF_AttributeSet(export_state, name="rof_present", value=rof_prognostic, rc=rc)
    if( rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(export_state, name="rof_prognostic", value=rof_prognostic, rc=rc)
    if( rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(export_state, name="rof_nx", value=rtmlon, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(export_state, name="rof_ny", value=rtmlat, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(export_state, name="flood_present", value=flood_present, rc=rc)
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

    call ESMF_AttributeSet(comp, "ShortName", "MOSART", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "LongName", &
                           "River Runoff from MOSART", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "Description", &
                  "The MOSART component is " // &
                  "the river runoff model used in the CESM1.1.  " // &
                  "More information on the MOSART project " // &
                  "and access to previous MOSART model versions and " // &
                  "documentation can be found via the MOSART Web Page.", &
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
    ! Run mosart model
    !
    ! !ARGUMENTS:
    implicit none
    type(ESMF_GridComp)  :: comp            ! MOSART gridded component
    type(ESMF_State)     :: import_state    ! MOSART import state
    type(ESMF_State)     :: export_state    ! MOSART export state
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
    integer :: ymd                        ! MOSART current date (YYYYMMDD)
    integer :: yr                         ! MOSART current year
    integer :: mon                        ! MOSART current month
    integer :: day                        ! MOSART current day
    integer :: tod                        ! MOSART current time of day (sec)
    logical :: rstwr_sync                 ! .true. ==> write restart file before returning
    logical :: rstwr                      ! .true. ==> write restart file before returning
    logical :: nlend_sync                 ! Flag signaling last time-step
    logical :: nlend                      ! .true. ==> last time-step
    integer :: shrlogunit,shrloglev       ! old values
    integer :: begg, endg                 ! Beginning and ending gridcell index numbers
    integer :: lbnum                      ! input to memory diagnostic
    integer :: g,i,kf                     ! counters
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

    ! Map ESMF to mosart input (rof) data type

    call t_startf ('lc_rof_import')
    call ESMF_StateGet(import_state, itemName="x2d", array=x2r, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call rof_import_esmf( x2r, rc=rc )
    call t_stopf ('lc_rof_import')

    ! Run mosart (input is *runin, output is rtmCTL%runoff)
    ! First advance mosart time step

    write(rdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr_sync,mon_sync,day_sync,tod_sync
    nlend = seq_timemgr_StopAlarmIsOn( EClock )
    rstwr = seq_timemgr_RestartAlarmIsOn( EClock )
    call advance_timestep()
    call Rtmrun(rstwr,nlend,rdate)

    ! Map roff data to MCT datatype (input is rtmCTL%runoff, output is r2x_r)
      
    call t_startf ('lc_rof_export')
    call ESMF_StateGet(export_state, itemName="d2x", array=r2x, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call rof_export_esmf( r2x, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call t_stopf ('lc_rof_export')
       
    call get_curr_date( yr, mon, day, tod)
    ymd = yr*10000 + mon*100 + day
    tod = tod
    if ( .not. seq_timemgr_EClockDateInSync( EClock, ymd, tod ) )then
       call seq_timemgr_EclockGetData( EClock, curr_ymd=ymd_sync, curr_tod=tod_sync )
       write(iulog,*)' mosart ymd=',ymd     ,'  mosart tod= ',tod
       write(iulog,*)'sync ymd=',ymd_sync,' sync tod= ',tod_sync
       call shr_sys_abort( sub//":: MOSART clock not in sync with Master Sync clock" )
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
    ! Finalize mosart
    !
    ! !ARGUMENTS:
    implicit none
    type(ESMF_GridComp)  :: comp            ! MOSART gridded component
    type(ESMF_State)     :: import_state    ! MOSART import state
    type(ESMF_State)     :: export_state    ! MOSART export state
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

  function rof_DistGrid_esmf(gsize, rc)

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Set the ESMF distributed grid for the runoff model
    !
    ! RETURN
    implicit none
    type(ESMF_DistGrid)  :: rof_DistGrid_esmf  
    !
    ! ARGUMENTS:
    integer, intent(out) :: gsize              ! Grid size
    integer, intent(out) :: rc                 ! Return code
    !
    ! Local Variables
    integer,allocatable :: gindex(:)         ! indexing for runoff grid cells
    integer :: n, ni                         ! indices
    integer :: lsize                         ! size of runoff data and number of grid cells
    integer :: begr, endr                    ! beg, end runoff indices
    integer :: ier                           ! error code
    character(len=32), parameter :: sub = 'rof_DistGrid_esmf'
    !-------------------------------------------------------------------

    ! Build the rof grid numbering
    ! NOTE:  Numbering scheme is: West to East and South to North
    ! starting at south pole.  Should be the same as what's used in SCRIP
    
    rc = ESMF_SUCCESS
    
    begr  = rtmCTL%begr
    endr  = rtmCTL%endr
    lsize = rtmCTL%lnumr
    gsize = rtmlon*rtmlat

    ! Check 
    ni = 0
    do n = begr,endr
       ni = ni + 1
       if (ni > lsize) then
          write(iulog,*) sub, ' : ERROR runoff count',n,ni,rtmCTL%lnumr
          call shr_sys_abort( sub//' ERROR: runoff > expected' )
       endif
    end do
    if (ni /= lsize) then
       write(iulog,*) sub, ' : ERROR runoff total count',ni,rtmCTL%lnumr
       call shr_sys_abort( sub//' ERROR: runoff not equal to expected' )
    endif

    ! Determine rof  distgrid
    allocate(gindex(lsize),stat=ier)
    ni = 0
    do n = begr,endr
       ni = ni + 1
       gindex(ni) = rtmCTL%gindex(n)
    end do
    rof_DistGrid_esmf = mct2esmf_init(gindex, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    deallocate(gindex)

  end function rof_DistGrid_esmf

!===============================================================================

  subroutine rof_domain_esmf( dom, rc )

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    !
    ! Send the runoff model domain information to the coupler
    !
    ! !ARGUMENTS:
    implicit none
    type(ESMF_Array), intent(inout)     :: dom   ! Domain information
    integer, intent(out)                :: rc    ! Return code
    !
    ! Local Variables
    integer  :: n, ni                          ! index
    integer  :: klon,klat,karea,kmask,kfrac    ! domain fields
    real(r8) :: re = SHR_CONST_REARTH*0.001_r8 ! radius of earth (km)
    real(R8), pointer :: fptr(:,:)
    character(len=32), parameter :: sub = 'rof_domain_esmf'
    !-------------------------------------------------------------------
    !
    ! Initialize domain type
    ! lat/lon in degrees,  area in radians^2
    ! 
    rc = ESMF_SUCCESS

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
    do n = rtmCTL%begr,rtmCTL%endr
       ni = ni + 1
       fptr(klon , ni) = rtmCTL%lonc(n)
       fptr(klat , ni) = rtmCTL%latc(n)
       fptr(karea, ni) = rtmCTL%area(n)*1.0e-6_r8/(re*re)
       fptr(kmask, ni) = 1.0_r8
       fptr(kfrac, ni) = 1.0_r8
    end do

  end subroutine rof_domain_esmf

!====================================================================================
 
  subroutine rof_import_esmf( x2r_array, rc)

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    ! Obtain the runoff input from the coupler
    ! convert from kg/m2s to m3/s
    !
    ! ARGUMENTS:
    implicit none
    type(ESMF_Array), intent(inout) :: x2r_array
    integer         , intent(out)   :: rc
    !
    ! LOCAL VARIABLES
    real(R8), pointer :: fptr(:, :)
    integer :: n2, n, nt, begr, endr, nliq, nfrz
    character(len=32), parameter :: sub = 'rof_import_mct'
    !---------------------------------------------------------------------------
    
    rc = ESMF_SUCCESS

    nliq = 0
    nfrz = 0
    do nt = 1,nt_rtm
       if (trim(rtm_tracers(nt)) == 'LIQ') then
          nliq = nt
       endif
       if (trim(rtm_tracers(nt)) == 'ICE') then
          nfrz = nt
       endif
    enddo
    if (nliq == 0 .or. nfrz == 0) then
       write(iulog,*) sub,': ERROR in rtm_tracers LIQ ICE ',nliq,nfrz,rtm_tracers
       call shr_sys_abort()
    endif

    call ESMF_ArrayGet(x2r_array, localDe=0, farrayPtr=fptr, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    
    begr = rtmCTL%begr
    endr = rtmCTL%endr
    do n = begr,endr
       n2 = n - begr + 1

       rtmCTL%qsur(n,nliq) = fptr(index_x2r_Flrl_rofsur,n2) * (rtmCTL%area(n)*0.001_r8)
       rtmCTL%qsub(n,nliq) = fptr(index_x2r_Flrl_rofsub,n2) * (rtmCTL%area(n)*0.001_r8)
       rtmCTL%qgwl(n,nliq) = fptr(index_x2r_Flrl_rofgwl,n2) * (rtmCTL%area(n)*0.001_r8)
       if (index_x2r_Flrl_rofdto > 0) then
          rtmCTL%qdto(n,nliq) = fptr(index_x2r_Flrl_rofdto,n2) * (rtmCTL%area(n)*0.001_r8)
       else
          rtmCTL%qdto(n,nliq) = 0.0_r8
       endif

       rtmCTL%qsur(n,nfrz) = fptr(index_x2r_Flrl_rofi,n2) * (rtmCTL%area(n)*0.001_r8)
       rtmCTL%qsub(n,nfrz) = 0.0_r8
       rtmCTL%qgwl(n,nfrz) = 0.0_r8
       rtmCTL%qdto(n,nfrz) = 0.0_r8

    enddo

  end subroutine rof_import_esmf

!====================================================================================

  subroutine rof_export_esmf( r2x_array, rc)

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Send the mosart export state to the CESM coupler
    ! convert from kg/m2s to m3/s
    !
    ! !ARGUMENTS:
    implicit none
    type(ESMF_Array), intent(inout)            :: r2x_array
    integer, intent(out)                       :: rc
    !
    ! Local variables
    integer :: ni, n, nt, nliq, nfrz
    real(R8), pointer :: fptr(:, :)
    logical,save :: first_time = .true.
    character(len=*), parameter :: sub = 'rof_export_esmf'
    !---------------------------------------------------------------------------
    
    rc = ESMF_SUCCESS

    nliq = 0
    nfrz = 0
    do nt = 1,nt_rtm
       if (trim(rtm_tracers(nt)) == 'LIQ') then
          nliq = nt
       endif
       if (trim(rtm_tracers(nt)) == 'ICE') then
          nfrz = nt
       endif
    enddo
    if (nliq == 0 .or. nfrz == 0) then
       write(iulog,*)'RtmUpdateInput: ERROR in rtm_tracers LIQ ICE ',nliq,nfrz,rtm_tracers
       call shr_sys_abort()
    endif

    call ESMF_ArrayGet(r2x_array, localDe=0, farrayPtr=fptr, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    fptr(:,:) = 0._r8

    if (first_time) then
       if (masterproc) then
       if ( ice_runoff )then
          write(iulog,*)'Snow capping will flow out in frozen river runoff'
       else
          write(iulog,*)'Snow capping will flow out in liquid river runoff'
       endif
       endif
       first_time = .false.
    end if

    ni = 0
    if ( ice_runoff )then
       ! separate liquid and ice runoff
       do n = rtmCTL%begr,rtmCTL%endr
          ni = ni + 1
          fptr(index_r2x_Forr_rofl,ni) =  rtmCTL%direct(n,nliq) / (rtmCTL%area(n)*0.001_r8)
          fptr(index_r2x_Forr_rofi,ni) =  rtmCTL%direct(n,nfrz) / (rtmCTL%area(n)*0.001_r8)
          if (rtmCTL%mask(n) >= 2) then
             fptr(index_r2x_Forr_rofl,ni) = fptr(index_r2x_Forr_rofl,ni) + &
                 rtmCTL%runoff(n,nliq) / (rtmCTL%area(n)*0.001_r8)
             fptr(index_r2x_Forr_rofi,ni) = fptr(index_r2x_Forr_rofi,ni) + &
                 rtmCTL%runoff(n,nfrz) / (rtmCTL%area(n)*0.001_r8)
             if (ni > rtmCTL%lnumr) then
                write(iulog,*) sub, ' : ERROR runoff count',n,ni
                call shr_sys_abort( sub//' : ERROR runoff > expected' )
             endif
          endif
       end do
    else
       ! liquid and ice runoff added to liquid runoff, ice runoff is zero
       do n = rtmCTL%begr,rtmCTL%endr
          ni = ni + 1
          fptr(index_r2x_Forr_rofl,ni) =  &
             (rtmCTL%direct(n,nfrz)+rtmCTL%direct(n,nliq)) / (rtmCTL%area(n)*0.001_r8)
          if (rtmCTL%mask(n) >= 2) then
             fptr(index_r2x_Forr_rofl,ni) = fptr(index_r2x_Forr_rofl,ni) +  &
               (rtmCTL%runoff(n,nfrz)+rtmCTL%runoff(n,nliq))/(rtmCTL%area(n)*0.001_r8)
             if (ni > rtmCTL%lnumr) then
                write(iulog,*) sub, ' : ERROR runoff count',n,ni
                call shr_sys_abort( sub//' : ERROR runoff > expected' )
             endif
          endif
       end do
    end if

    ! Flooding back to land, sign convention is positive in land->rof direction
    ! so if water is sent from rof to land, the flux must be negative.
    ni = 0
    do n = rtmCTL%begr, rtmCTL%endr
       ni = ni + 1
       fptr(index_r2x_Flrr_flood,ni) = -rtmCTL%flood(n)/(rtmCTL%area(n)*0.001_r8)
       fptr(index_r2x_Flrr_volr,ni)    = (Trunoff%wr(n,nliq) + Trunoff%wt(n,nliq)) / rtmCTL%area(n)
       fptr(index_r2x_Flrr_volrmch,ni) = Trunoff%wr(n,nliq) / rtmCTL%area(n)
    end do

  end subroutine rof_export_esmf

#endif
end module rof_comp_esmf
