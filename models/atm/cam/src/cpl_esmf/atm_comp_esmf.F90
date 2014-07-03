module atm_comp_esmf

  use pio              , only: file_desc_t, io_desc_t, var_desc_t, pio_double, pio_def_dim, &
                               pio_put_att, pio_enddef, pio_initdecomp, pio_read_darray, pio_freedecomp, &
                               pio_closefile, pio_write_darray, pio_def_var, pio_inq_varid, &
                               pio_noerr, pio_bcast_error, pio_internal_error, pio_seterrorhandling
  use esmf
  use esmfshr_mod
  use seq_flds_mod
  use seq_infodata_mod
  use seq_timemgr_mod

  use shr_kind_mod     , only: r8 => shr_kind_r8, cl=>shr_kind_cl
  use shr_file_mod     , only: shr_file_getunit, shr_file_freeunit, &
                               shr_file_setLogUnit, shr_file_setLogLevel, &
                               shr_file_getLogUnit, shr_file_getLogLevel, &
		               shr_file_setIO
  use shr_sys_mod      , only: shr_sys_flush, shr_sys_abort

  use cam_cpl_indices
  use cam_comp
  use cam_instance     , only: cam_instance_init, inst_suffix
  use cam_control_mod  , only: nsrest, adiabatic, ideal_phys, aqua_planet, eccen, obliqr, lambm0, mvelpp
  use radiation        , only: radiation_get, radiation_do, radiation_nextsw_cday
  use phys_grid        , only: get_ncols_p, get_gcol_all_p, & 
                               ngcols, get_gcol_p, get_rlat_all_p, &
	                       get_rlon_all_p, get_area_all_p
  use ppgrid           , only: pcols, begchunk, endchunk       
  use dyn_grid         , only: get_horiz_grid_dim_d
  use camsrfexch       , only: cam_out_t, cam_in_t     
  use cam_restart      , only: get_restcase, get_restartdir
  use cam_history      , only: outfld, ctitle
  use abortutils       , only: endrun
  use filenames        , only: interpret_filename_spec, caseid, brnch_retain_casename
#ifdef SPMD
  use spmd_utils       , only: spmdinit, masterproc, iam
  use mpishorthand     , only: mpicom
#else
  use spmd_utils       , only: spmdinit, masterproc, mpicom, iam
#endif
  use time_manager     , only: get_curr_calday, advance_timestep, get_curr_date, get_nstep, &
                               is_first_step, get_step_size, timemgr_init, timemgr_check_restart
  use ioFileMod             
  use perf_mod
  use cam_logfile      , only: iulog
  use co2_cycle        , only: c_i, co2_readFlux_ocn, co2_readFlux_fuel, co2_transport, &
                               co2_time_interp_ocn, co2_time_interp_fuel, data_flux_ocn, data_flux_fuel
  use physconst       ,  only: mwco2
  use runtime_opts     , only: read_namelist
  use phys_control     , only: cam_chempkg_is

!
! !PUBLIC TYPES:
  implicit none
  save
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: atm_register_esmf
  public :: atm_init_esmf
  public :: atm_run_esmf
  public :: atm_final_esmf

!--------------------------------------------------------------------------
! Private interfaces
!--------------------------------------------------------------------------

  private :: atm_DistGrid_esmf
  private :: atm_import_esmf
  private :: atm_export_esmf
  private :: atm_domain_esmf
  private :: atm_read_srfrest_esmf
  private :: atm_write_srfrest_esmf

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

  type(cam_in_t) , pointer :: cam_in(:)
  type(cam_out_t), pointer :: cam_out(:)

  type(ESMF_Array)   :: a2x_a_SNAP
  type(ESMF_Array)   :: a2x_a_SUM

  integer, parameter  :: nlen = 256     ! Length of character strings
  character(len=nlen) :: fname_srf_cam  ! surface restart filename
  character(len=nlen) :: pname_srf_cam  ! surface restart full pathname

  ! Filename specifier for restart surface file
  character(len=cl) :: rsfilename_spec_cam

!
! Time averaged counter for flux fields
!
  integer :: avg_count
!
! Time averaged flux fields
!  
  character(*), parameter :: a2x_avg_flds = "Faxa_rainc:Faxa_rainl:Faxa_snowc:Faxa_snowl"  
!
! Are all surface types present   
!
  logical :: lnd_present ! if true => land is present
  logical :: ocn_present ! if true => ocean is present

!
!================================================================================
CONTAINS
!================================================================================

  subroutine atm_register_esmf(comp, rc)
    implicit none
    type(ESMF_GridComp)  :: comp
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS
    ! Register the callback routines.

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
      atm_init_esmf, phase=1, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, &
      atm_run_esmf, phase=1, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, &
      atm_final_esmf, phase=1, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

  end subroutine

  subroutine atm_init_esmf(comp, import_state, export_state, EClock, rc)

    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    implicit none
    type(ESMF_GridComp)          :: comp
    type(ESMF_State)             :: import_state
    type(ESMF_State)             :: export_state
    type(ESMF_Clock)             :: EClock
    integer, intent(out)         :: rc
    !
    ! Locals
    !
    integer :: ATMID            ! cesm ID value
    integer :: mpicom_atm, mpicom_vm
    integer :: iradsw
    logical :: exists           ! true if file exists
    real(r8):: nextsw_cday      ! calendar of next atm shortwave
    integer :: stepno           ! time step			 
    integer :: dtime_sync       ! integer timestep size
    integer :: currentymd       ! current year-month-day
    integer :: dtime            ! time step increment (sec)
    integer :: atm_cpl_dt       ! driver atm coupling time step 
    integer :: nstep            ! CAM nstep
    real(r8):: caldayp1         ! CAM calendar day for for next cam time step
    integer :: dtime_cam        ! Time-step increment (sec)
    integer :: ymd              ! CAM current date (YYYYMMDD)
    integer :: yr               ! CAM current year
    integer :: mon              ! CAM current month
    integer :: day              ! CAM current day
    integer :: tod              ! CAM current time of day (sec)
    integer :: start_ymd        ! Start date (YYYYMMDD)
    integer :: start_tod        ! Start time of day (sec)
    integer :: ref_ymd          ! Reference date (YYYYMMDD)
    integer :: ref_tod          ! Reference time of day (sec)
    integer :: stop_ymd         ! Stop date (YYYYMMDD)
    integer :: stop_tod         ! Stop time of day (sec)
    logical :: perpetual_run    ! If in perpetual mode or not
    integer :: perpetual_ymd    ! Perpetual date (YYYYMMDD)
    logical :: single_column
    real(r8):: scmlat,scmlon
    integer :: shrlogunit,shrloglev ! old values
    logical :: first_time = .true.
    character(len=SHR_KIND_CS) :: calendar  ! Calendar type
    character(len=SHR_KIND_CS) :: starttype ! infodata start type
    integer :: lbnum
    type(ESMF_DistGrid)                   :: distgrid
    type(ESMF_Array)                      :: d2x, x2d, dom
    type(ESMF_VM)                         :: vm
    integer                               :: nflds, gsize
    integer :: hdim1_d, hdim2_d ! dimensions of rectangular horizontal grid
                                ! data structure, If 1D data structure, then
                                ! hdim2_d == 1.
    character(ESMF_MAXSTR) :: convCIM, purpComp
    character(len=64) :: filein ! Input namelist filename
    !-----------------------------------------------------------------------
    !
    ! Determine cdata points
    !
    rc = ESMF_SUCCESS
#if (defined _MEMTRACE)
    if(masterproc) then
      lbnum=1
      call memmon_dump_fort('memmon.out','atm_init_esmf:start::',lbnum)
    endif                      
#endif                         

    ! duplicate the mpi communicator from the current VM 
    call ESMF_VMGetCurrent(vm, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_VMGet(vm, mpiCommunicator=mpicom_vm, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call MPI_Comm_dup(mpicom_vm, mpicom_atm, rc)
    if(rc /= 0) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    if (first_time) then
       
       ! Initialize cam id

       call ESMF_AttributeGet(export_state, name="ID", value=ATMID, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       call cam_instance_init(ATMID)

       ! Set filename specifier for restart surface file
       ! (%c=caseid, $y=year, $m=month, $d=day, $s=seconds in day)
       rsfilename_spec_cam = '%c.cam' // trim(inst_suffix) // '.rs.%y-%m-%d-%s.nc' 

       ! Determine attribute vector indices

       call cam_cpl_indices_set()

       ! Redirect share output to cam log
       
       call spmdinit(mpicom_atm)
       
       if (masterproc) then
          inquire(file='atm_modelio.nml'//trim(inst_suffix), exist=exists)
          if (exists) then
             iulog = shr_file_getUnit()
             call shr_file_setIO('atm_modelio.nml'//trim(inst_suffix), iulog)
          endif
          write(iulog,*) "CAM atmosphere model initialization"
       endif
       
       call shr_file_getLogUnit (shrlogunit)
       call shr_file_getLogLevel(shrloglev)
       call shr_file_setLogUnit (iulog)
       ! 
       ! Consistency check                              
       !
       if (co2_readFlux_ocn .and. index_x2a_Faoo_fco2_ocn /= 0) then
          write(iulog,*)'error co2_readFlux_ocn and index_x2a_Faoo_fco2_ocn cannot both be active'
          call shr_sys_abort()
       end if
       ! 
       ! Get data from infodata object
       !
       call ESMF_AttributeGet(export_state, name="case_name", value=caseid, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       call ESMF_AttributeGet(export_state, name="case_desc", value=ctitle, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       call ESMF_AttributeGet(export_state, name="start_type", value=starttype, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       call ESMF_AttributeGet(export_state, name="atm_adiabatic", value=adiabatic, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       call ESMF_AttributeGet(export_state, name="atm_ideal_phys", value=ideal_phys, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       call ESMF_AttributeGet(export_state, name="aqua_planet", value=aqua_planet, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       call ESMF_AttributeGet(export_state, name="brnch_retain_casename", value=brnch_retain_casename, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       call ESMF_AttributeGet(export_state, name="single_column", value=single_column, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       call ESMF_AttributeGet(export_state, name="scmlat", value=scmlat, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       call ESMF_AttributeGet(export_state, name="scmlon", value=scmlon, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       call ESMF_AttributeGet(export_state, name="orb_eccen", value=eccen, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       call ESMF_AttributeGet(export_state, name="orb_mvelpp", value=mvelpp, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       call ESMF_AttributeGet(export_state, name="orb_lambm0", value=lambm0, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       call ESMF_AttributeGet(export_state, name="orb_obliqr", value=obliqr, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       call ESMF_AttributeGet(export_state, name="lnd_present", value=lnd_present, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       call ESMF_AttributeGet(export_state, name="ocn_present", value=ocn_present, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       call ESMF_AttributeGet(export_state, name="perpetual", value=perpetual_run, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       call ESMF_AttributeGet(export_state, name="perpetual_ymd", value=perpetual_ymd, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       !
       ! Get nsrest from startup type methods
       !
       if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
          nsrest = 0
       else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
          nsrest = 1
       else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
          nsrest = 3
       else
          write(iulog,*) 'atm_comp_esmf: ERROR: unknown starttype'
          call shr_sys_abort()
       end if
       !
       ! Initialize time manager.
       !
       call seq_timemgr_EClockGetData(EClock, &
                                      start_ymd=start_ymd, start_tod=start_tod, &
                                      ref_ymd=ref_ymd, ref_tod=ref_tod,         &
                                      stop_ymd=stop_ymd, stop_tod=stop_tod,     &
                                      calendar=calendar )
       !
       ! Read namelist
       !
       filein = "atm_in" // trim(inst_suffix)
       call read_namelist(single_column_in=single_column, scmlat_in=scmlat, &
            scmlon_in=scmlon, nlfilename_in=filein)
       !
       ! Initialize cam time manager
       !
       if ( nsrest == 0 )then
          call timemgr_init( calendar_in=calendar, start_ymd=start_ymd, &
                             start_tod=start_tod, ref_ymd=ref_ymd,      &
                             ref_tod=ref_tod, stop_ymd=stop_ymd,        &
                             stop_tod=stop_tod,                         &
                             perpetual_run=perpetual_run,               &
                             perpetual_ymd=perpetual_ymd )
       end if
       !
       ! First phase of cam initialization 
       ! Initialize mpicom_atm, allocate cam_in and cam_out and determine 
       ! atm decomposition (needed to initialize gsmap) 
       ! for an initial run, cam_in and cam_out are allocated in cam_initial
       ! for a restart/branch run, cam_in and cam_out are allocated in restart 
       ! Set defaults then override with user-specified input and initialize time manager
       ! Note that the following arguments are needed to cam_init for timemgr_restart only
       !
       call cam_init( cam_out, cam_in, mpicom_atm, &
                      start_ymd, start_tod, ref_ymd, ref_tod, stop_ymd, stop_tod, &
                      perpetual_run, perpetual_ymd, calendar)
       !
       ! Check consistency of restart time information with input clock
       !
       if (nsrest /= 0) then
          dtime_cam = get_step_size()
          call timemgr_check_restart( calendar, start_ymd, start_tod, ref_ymd, &
                                      ref_tod, dtime_cam, perpetual_run, perpetual_ymd)
       end if
       !
       ! Initialize Distgrid
       !
       distgrid = atm_DistGrid_esmf(gsize, rc) 
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       call ESMF_AttributeSet(export_state, name="gsize", value=gsize, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       ! Initialize Arrays

       dom = mct2esmf_init(distgrid, attname=seq_flds_dom_fields, name="domain", rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       call atm_domain_esmf(dom, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       d2x = mct2esmf_init(distgrid, attname=seq_flds_a2x_fields, name="d2x", rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       
       x2d = mct2esmf_init(distgrid, attname=seq_flds_x2a_fields, name="x2d", rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       call ESMF_StateAdd(export_state, (/dom/), rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       call ESMF_StateAdd(export_state, (/d2x/), rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    
       call ESMF_StateAdd(import_state, (/x2d/), rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       
       a2x_a_SNAP = mct2esmf_init(distgrid, attname=a2x_avg_flds, name="a2x_a_SNAP", rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       
       a2x_a_SUM = mct2esmf_init(distgrid, attname=a2x_avg_flds, name="a2x_a_SUM", rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       !
       ! Initialize averaging counter
       !
       avg_count = 0
       !
       ! Create initial atm export state
       !
       call atm_export_esmf( cam_out, d2x, rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       !
       ! Set flag to specify that an extra albedo calculation is to be done (i.e. specify active)
       !
       call ESMF_AttributeSet(export_state, name="atm_prognostic", value=.true., rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       call get_horiz_grid_dim_d(hdim1_d, hdim2_d)
       call ESMF_AttributeSet(export_state, name="atm_nx", value=hdim1_d, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       call ESMF_AttributeSet(export_state, name="atm_ny", value=hdim2_d, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       ! Set flag to indicate that CAM will provide carbon and dust deposition fluxes.
       ! This is now hardcoded to .true. since the ability of CICE to read these
       ! fluxes from a file has been removed.
       call ESMF_AttributeSet(export_state, name="atm_aero", value=.true., rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       !
       ! Set time step of radiation computation as the current calday
       ! This will only be used on the first timestep of an initial run
       !
       if (nsrest == 0) then
          nextsw_cday = get_curr_calday()
          call ESMF_AttributeSet(export_state, name="nextsw_cday", value=nextsw_cday, rc=rc)
          if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       end if
       
       ! End redirection of share output to cam log
       
       call shr_file_setLogUnit (shrlogunit)
       call shr_file_setLogLevel(shrloglev)

       first_time = .false.

    else
       
       ! For initial run, run cam radiation/clouds and return
       ! For restart run, read restart x2a_a
       ! Note - a2x_a is computed upon the completion of the previous run - cam_run1 is called
       ! only for the purposes of finishing the flux averaged calculation to compute a2x_a
       ! Note - cam_run1 is called on restart only to have cam internal state consistent with the 
       ! a2x_a state sent to the coupler

       ! Redirect share output to cam log

       call shr_file_getLogUnit (shrlogunit)
       call shr_file_getLogLevel(shrloglev)
       call shr_file_setLogUnit (iulog)

       call ESMF_StateGet(export_state, itemName="domain", array=dom, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       call ESMF_StateGet(export_state, itemName="d2x", array=d2x, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       call ESMF_StateGet(import_state, itemName="x2d", array=x2d, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       call seq_timemgr_EClockGetData(EClock,curr_ymd=CurrentYMD, StepNo=StepNo, dtime=DTime_Sync )
       if (StepNo == 0) then
          call atm_import_esmf( x2d, cam_in, rc )
          if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
          call cam_run1 ( cam_in, cam_out ) 
          call atm_export_esmf( cam_out, d2x, rc )
          if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       else
          call atm_read_srfrest_esmf( EClock, dom, x2d, d2x )
          call atm_import_esmf( x2d, cam_in, rc)
          if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
          call cam_run1 ( cam_in, cam_out ) 
       end if

       ! Compute time of next radiation computation, like in run method for exact restart

! tcx was
!       nextsw_cday = radiation_nextsw_cday() 

       call seq_timemgr_EClockGetData(Eclock,dtime=atm_cpl_dt)
       dtime = get_step_size()          
       nstep = get_nstep()
       if (nstep < 1 .or. dtime < atm_cpl_dt) then
          nextsw_cday = radiation_nextsw_cday() 
       else if (dtime == atm_cpl_dt) then
          caldayp1 = get_curr_calday(offset=int(dtime))
          nextsw_cday = radiation_nextsw_cday() 
          if (caldayp1 /= nextsw_cday) nextsw_cday = -1._r8
       else
          call shr_sys_abort('dtime must be less than or equal to atm_cpl_dt')
       end if
       call ESMF_AttributeSet(export_state, name="nextsw_cday", value=nextsw_cday, rc=rc)
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       ! End redirection of share output to cam log
       
       call shr_file_setLogUnit (shrlogunit)
       call shr_file_setLogLevel(shrloglev)
       
    end if

#if (defined _MEMTRACE )
    if(masterproc) then
      lbnum=1
      call memmon_dump_fort('memmon.out','atm_init_esmf:end::',lbnum)
      call memmon_reset_addr()
    endif
#endif

!    convCIM  = "CIM 1.0"
    convCIM  = "CIM"
    purpComp = "Model Component Simulation Description"

    call ESMF_AttributeAdd(comp,  &
                           convention=convCIM, purpose=purpComp, rc=rc)

    call ESMF_AttributeSet(comp, "ShortName", "CAM", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "LongName", &
                           "Community Atmosphere Model", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "Description", &
                  "Version 5.0 of the Community Atmosphere Model (CAM) " // &
                  "is the latest in a series of global atmosphere models " // &
                  "developed primarily at the National Center for " // &
                  "Atmospheric Research (NCAR).  CAM 5.0 includes " // &
                  "significant enhancements to the representation of " // &
                  "atmospheric processes resulting in a number of notable " // &
                  "improvements.  CAM 4.0 is also available in the CESM " // &
                  "1.0 release.  Development of the model was led by the " // &
                  "Atmosphere Model Working Group (AMWG).", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ReleaseDate", "2010", &
                           convention=convCIM, purpose=purpComp, rc=rc)
    call ESMF_AttributeSet(comp, "ModelType", "Atmosphere", &
                           convention=convCIM, purpose=purpComp, rc=rc)

!    call ESMF_AttributeSet(comp, "Name", "Cecile Hannay", &
!                           convention=convCIM, purpose=purpComp, rc=rc)
!    call ESMF_AttributeSet(comp, "EmailAddress", &
!                           "hannay@ucar.edu", &
!                           convention=convCIM, purpose=purpComp, rc=rc)
!    call ESMF_AttributeSet(comp, "ResponsiblePartyRole", "contact", &
!                           convention=convCIM, purpose=purpComp, rc=rc)

    call shr_sys_flush(iulog)

 end subroutine atm_init_esmf

!================================================================================

subroutine atm_run_esmf(comp, import_state, export_state, EClock, rc)

    !-----------------------------------------------------------------------
    !
    ! Uses
    !
    use time_manager,    only: advance_timestep, get_curr_date, get_curr_calday, &
	                       get_nstep, get_step_size
    use scamMod,         only: single_column
!   use iop,             only: scam_use_iop_srf
    use pmgrid,          only: plev, plevp
    use constituents,    only: pcnst
    use shr_sys_mod,     only: shr_sys_flush
    use chemistry,       only: chem_reset_fluxes

    implicit none
    ! 
    ! Arguments
    !
    type(ESMF_GridComp)          :: comp
    type(ESMF_State)             :: import_state
    type(ESMF_State)             :: export_state
    type(ESMF_Clock)             :: EClock
    integer, intent(out)         :: rc
    !
    ! Local variables
    !
    integer :: StepNo          ! time step			 
    integer :: DTime_Sync      ! integer timestep size
    integer :: CurrentYMD      ! current year-month-day
    integer :: iradsw          ! shortwave radation frequency (time steps) 
    logical :: dosend          ! true => send data back to driver
    integer :: dtime           ! time step increment (sec)
    integer :: atm_cpl_dt      ! driver atm coupling time step 
    integer :: ymd_sync        ! Sync date (YYYYMMDD)
    integer :: yr_sync         ! Sync current year
    integer :: mon_sync        ! Sync current month
    integer :: day_sync        ! Sync current day
    integer :: tod_sync        ! Sync current time of day (sec)
    integer :: ymd             ! CAM current date (YYYYMMDD)
    integer :: yr              ! CAM current year
    integer :: mon             ! CAM current month
    integer :: day             ! CAM current day
    integer :: tod             ! CAM current time of day (sec)
    integer :: nstep           ! CAM nstep
    integer :: shrlogunit,shrloglev ! old values
    real(r8):: caldayp1        ! CAM calendar day for for next cam time step
    real(r8):: nextsw_cday     ! calendar of next atm shortwave
    logical :: rstwr           ! .true. ==> write restart file before returning
    logical :: nlend           ! Flag signaling last time-step
    logical :: rstwr_sync      ! .true. ==> write restart file before returning
    logical :: nlend_sync      ! Flag signaling last time-step
    logical :: first_time = .true.
    character(len=*), parameter :: subname="atm_run_esmf"
    !-----------------------------------------------------------------------
    integer :: lbnum
    type(ESMF_Array) :: d2x, x2d, dom

    real(R8), pointer :: fptr(:,:)        ! pointer into    array data

    rc = ESMF_SUCCESS

#if (defined _MEMTRACE)
    if(masterproc) then
      lbnum=1
      call memmon_dump_fort('memmon.out',SubName //':start::',lbnum)
    endif
#endif

    ! Redirect share output to cam log
    
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iulog)
    
    ! Note that sync clock time should match cam time at end of time step/loop not beginning
    
    call seq_timemgr_EClockGetData(EClock,curr_ymd=ymd_sync,curr_tod=tod_sync, &
       curr_yr=yr_sync,curr_mon=mon_sync,curr_day=day_sync)

    nlend_sync = seq_timemgr_StopAlarmIsOn(EClock)
    rstwr_sync = seq_timemgr_RestartAlarmIsOn(EClock)
    
    ! Map input from Array to cam data structure
    call ESMF_StateGet(export_state, itemName="domain", array=dom, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_StateGet(export_state, itemName="d2x", array=d2x, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_StateGet(import_state, itemName="x2d", array=x2d, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !load orbital parameters
    call ESMF_AttributeGet(export_state, name="orb_eccen", value=eccen, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_AttributeGet(export_state, name="orb_mvelpp", value=mvelpp, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_AttributeGet(export_state, name="orb_lambm0", value=lambm0, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_AttributeGet(export_state, name="orb_obliqr", value=obliqr, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    
    call t_startf ('CAM_import')
    call atm_import_esmf(x2d, cam_in, rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call t_stopf  ('CAM_import')
    
    ! Cycle over all time steps in the atm coupling interval
    
    dosend = .false.
    do while (.not. dosend)
       
       ! (re)set surface fluxes of chem tracers here to MEGAN fluxes (from CLM)
       ! or to zero so that fluxes read from file can be added to MEGAN
       call ESMF_ArrayGet(x2d, localDe=0, farrayPtr=fptr, rc=rc)
       call chem_reset_fluxes( fptr, cam_in )

       ! Determine if dosend
       ! When time is not updated at the beginning of the loop - then return only if
       ! are in sync with clock before time is updated
       
       call get_curr_date( yr, mon, day, tod )
       ymd = yr*10000 + mon*100 + day
       tod = tod
       dosend = (seq_timemgr_EClockDateInSync( EClock, ymd, tod))
       
       ! Determine if time to write cam restart and stop
       
       rstwr = .false.
       if (rstwr_sync .and. dosend) rstwr = .true.
       nlend = .false.
       if (nlend_sync .and. dosend) nlend = .true.
       
       ! Single column specific input 
       
       if (single_column) then
          call scam_use_iop_srf( cam_in )
       endif

       ! Run CAM (run2, run3, run4)
       
       call t_startf ('CAM_run2')
       call cam_run2( cam_out, cam_in )
       call t_stopf  ('CAM_run2')

       call t_startf ('CAM_run3')
       call cam_run3( cam_out )
       call t_stopf  ('CAM_run3')
       
       call t_startf ('CAM_run4')
       call cam_run4( cam_out, cam_in, rstwr, nlend, &
            yr_spec=yr_sync, mon_spec=mon_sync, day_spec=day_sync, sec_spec=tod_sync)
       call t_stopf  ('CAM_run4')
       
       ! Advance cam time step 
       
       call t_startf ('CAM_adv_timestep')
       call advance_timestep()
       call t_stopf  ('CAM_adv_timestep')
       
       ! Run cam radiation/clouds (run1)
          
       call t_startf ('CAM_run1')
       call cam_run1 ( cam_in, cam_out ) 
       call t_stopf  ('CAM_run1')
       
       ! Map output from cam to Array data structures
       
       call t_startf ('CAM_export')
       call atm_export_esmf( cam_out, d2x, rc )
       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       call t_stopf ('CAM_export')
       
       ! Compute snapshot attribute vector for accumulation
       
       ! don't accumulate on first coupling freq ts1 and ts2
       ! for consistency with ccsm3 when flxave is off
       nstep = get_nstep()
       if (nstep <= 2) then
          call esmfshr_util_ArrayCopy( d2x, a2x_a_SUM ,rc)
          if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
          avg_count = 1
       else
          call esmfshr_util_ArrayCopy( d2x, a2x_a_SNAP ,rc)
          if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
          call esmfshr_util_ArraySum( a2x_a_SNAP, a2x_a_SUM ,rc)
          if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
          avg_count = avg_count + 1
       endif
       
    end do

    ! Finish accumulation of attribute vector and average and copy accumulation 
    ! field into output attribute vector
    
    call esmfshr_util_ArrayAvg(a2x_a_SUM, avg_count, rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call esmfshr_util_ArrayCopy(a2x_a_SUM, d2x, rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call esmfshr_util_ArrayZero( a2x_a_SUM, rc) 
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    avg_count = 0                   
    
    ! Get time of next radiation calculation - albedos will need to be 
    ! calculated by each surface model at this time
    
    call seq_timemgr_EClockGetData(Eclock,dtime=atm_cpl_dt)
    dtime = get_step_size()          
    if (dtime < atm_cpl_dt) then
       nextsw_cday = radiation_nextsw_cday() 
    else if (dtime == atm_cpl_dt) then
       caldayp1 = get_curr_calday(offset=int(dtime))
       nextsw_cday = radiation_nextsw_cday() 
       if (caldayp1 /= nextsw_cday) nextsw_cday = -1._r8
    else
       call shr_sys_abort('dtime must be less than or equal to atm_cpl_dt')
    end if
    call ESMF_AttributeSet(export_state, name="nextsw_cday", value=nextsw_cday, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    
    ! Write merged surface data restart file if appropriate
    
    if (rstwr_sync) then
       call atm_write_srfrest_esmf( dom, x2d, d2x, &
            yr_spec=yr_sync, mon_spec=mon_sync, day_spec=day_sync, sec_spec=tod_sync)
    end if
    
    ! Check for consistency of internal cam clock with master sync clock 
    
    dtime = get_step_size()
    call get_curr_date( yr, mon, day, tod, offset=-dtime )
    ymd = yr*10000 + mon*100 + day
    tod = tod
    if ( .not. seq_timemgr_EClockDateInSync( EClock, ymd, tod ) )then
       call seq_timemgr_EClockGetData(EClock, curr_ymd=ymd_sync, curr_tod=tod_sync )
       write(iulog,*)' cam ymd=',ymd     ,'  cam tod= ',tod
       write(iulog,*)'sync ymd=',ymd_sync,' sync tod= ',tod_sync
       call shr_sys_abort( subname//': CAM clock is not in sync with master Sync Clock' )
    end if
    
    ! End redirection of share output to cam log

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

#if (defined _MEMTRACE)
    if(masterproc) then
      lbnum=1
      call memmon_dump_fort('memmon.out',SubName //':end::',lbnum)
      call memmon_reset_addr()
    endif
#endif

  end subroutine atm_run_esmf

!================================================================================

  subroutine atm_final_esmf(comp, import_state, export_state, EClock, rc)

    !----- arguments -----
    implicit none
    type(ESMF_GridComp)          :: comp
    type(ESMF_State)             :: import_state
    type(ESMF_State)             :: export_state
    type(ESMF_Clock)             :: EClock
    integer, intent(out)         :: rc

    ! local
    type(ESMF_Array)                 :: a2x_a, x2a_a
    type(ESMF_DistGrid)              :: distgrid_ref

    !----------------------------------------------------------------------------
    ! Finalize routine 
    !----------------------------------------------------------------------------
    rc = ESMF_SUCCESS

    call cam_final( cam_out, cam_in )

    ! Destroy ESMF objects
    call esmfshr_util_StateArrayDestroy(export_state,"d2x",rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call esmfshr_util_StateArrayDestroy(export_state,"domain",rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call esmfshr_util_StateArrayDestroy(import_state,"x2d",rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

  end subroutine atm_final_esmf

!================================================================================
  type(ESMF_DistGrid) function atm_DistGrid_esmf(gsize, rc)
    use phys_grid, only : get_nlcols_p
    implicit none
    !-------------------------------------------------------------------
    !
    ! Arguments
    !
    integer, intent(out)            :: gsize
    integer, intent(out)            :: rc
    !
    ! Local variables
    !
    integer, allocatable :: gindex(:)
    integer :: i, n, c, ncols, sizebuf, nlcols
    !-------------------------------------------------------------------

    ! Build the atmosphere grid numbering for Array
    ! NOTE:  Numbering scheme is: West to East and South to North
    ! starting at south pole.  Should be the same as what's used in SCRIP
    
    ! Determine global seg map

    rc = ESMF_SUCCESS

    sizebuf=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       do i = 1,ncols
          sizebuf = sizebuf+1
       end do
    end do

    allocate(gindex(sizebuf))

    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       do i = 1,ncols
          n=n+1
          gindex(n) = get_gcol_p(c,i)
       end do
    end do

    gsize = ngcols
    atm_DistGrid_esmf = mct2esmf_init(gindex, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    deallocate(gindex)

  end function atm_DistGrid_esmf

!===============================================================================

  subroutine atm_import_esmf( x2a_a, cam_in ,rc)

    !-----------------------------------------------------------------------
    !
    ! Uses	
    !
    use dust_intr,     only: dust_idx1
#if (defined MODAL_AERO)
    use mo_chem_utls,  only: get_spc_ndx
#endif
    use shr_const_mod, only: shr_const_stebol
    use seq_drydep_mod,only: n_drydep

    implicit none
    !
    ! Arguments
    !
    type(ESMF_Array),   intent(inout) :: x2a_a
    type(cam_in_t),     intent(inout) :: cam_in(begchunk:endchunk)
    integer, intent(out)              :: rc
    !
    ! Local variables
    !		
    integer  :: i,lat,n,c,ig  ! indices
    integer  :: ncols         ! number of columns
    integer  :: dust_ndx
    logical, save :: first_time = .true.
#if (defined MODAL_AERO)
    integer, parameter:: ndst =2
    integer, target   :: spc_ndx(ndst)
#if (defined MODAL_AERO_7MODE)
    integer, pointer  :: dst_a5_ndx, dst_a7_ndx
#elif (defined MODAL_AERO_3MODE)
    integer, pointer  :: dst_a1_ndx, dst_a3_ndx
#endif
#endif
    real(R8), pointer :: fptr(:,:)        ! pointer into    array data

    !-----------------------------------------------------------------------
    !
    
    rc = ESMF_SUCCESS
    call ESMF_ArrayGet(x2a_a, localDe=0, farrayPtr=fptr, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

#if (defined MODAL_AERO)
#if (defined MODAL_AERO_7MODE)

    dst_a5_ndx => spc_ndx(1)
    dst_a7_ndx => spc_ndx(2)
    dst_a5_ndx = get_spc_ndx( 'dst_a5' )
    dst_a7_ndx = get_spc_ndx( 'dst_a7' )
#elif (defined MODAL_AERO_3MODE)
    dst_a1_ndx => spc_ndx(1)
    dst_a3_ndx => spc_ndx(2)
    dst_a1_ndx = get_spc_ndx( 'dst_a1' )
    dst_a3_ndx = get_spc_ndx( 'dst_a3' )
#endif
#endif

    ! ccsm sign convention is that fluxes are positive downward

    ig=1
    do c=begchunk,endchunk
       ncols = get_ncols_p(c)                                                 
       do i =1,ncols                                                               
          cam_in(c)%wsx(i)       = -fptr(index_x2a_Faxx_taux,ig)     
          cam_in(c)%wsy(i)       = -fptr(index_x2a_Faxx_tauy,ig)     
          cam_in(c)%lhf(i)       = -fptr(index_x2a_Faxx_lat, ig)     
          cam_in(c)%shf(i)       = -fptr(index_x2a_Faxx_sen, ig)     
          cam_in(c)%lwup(i)      = -fptr(index_x2a_Faxx_lwup,ig)    
          cam_in(c)%cflx(i,1)    = -fptr(index_x2a_Faxx_evap,ig)                
          cam_in(c)%asdir(i)     =  fptr(index_x2a_Sx_avsdr, ig)  
          cam_in(c)%aldir(i)     =  fptr(index_x2a_Sx_anidr, ig)  
          cam_in(c)%asdif(i)     =  fptr(index_x2a_Sx_avsdf, ig)  
          cam_in(c)%aldif(i)     =  fptr(index_x2a_Sx_anidf, ig)
          cam_in(c)%ts(i)        =  fptr(index_x2a_Sx_t,     ig)  
          cam_in(c)%sst(i)       =  fptr(index_x2a_So_t,     ig)             
          cam_in(c)%snowhland(i) =  fptr(index_x2a_Sl_snowh, ig)  
          cam_in(c)%snowhice(i)  =  fptr(index_x2a_Si_snowh, ig)  
          cam_in(c)%tref(i)      =  fptr(index_x2a_Sx_tref,  ig)  
          cam_in(c)%qref(i)      =  fptr(index_x2a_Sx_qref,  ig)
          cam_in(c)%u10(i)       =  fptr(index_x2a_Sx_u10,   ig)
          cam_in(c)%icefrac(i)   =  fptr(index_x2a_Sf_ifrac, ig)  
          cam_in(c)%ocnfrac(i)   =  fptr(index_x2a_Sf_ofrac, ig)
	  cam_in(c)%landfrac(i)  =  fptr(index_x2a_Sf_lfrac, ig)
          if ( associated(cam_in(c)%ram1) ) &
               cam_in(c)%ram1(i) =  fptr(index_x2a_Sl_ram1 , ig)
          if ( associated(cam_in(c)%fv) ) &
               cam_in(c)%fv(i)   =  fptr(index_x2a_Sl_fv   , ig)
          if ( associated(cam_in(c)%soilw) ) &
               cam_in(c)%soilw(i)   =  fptr(index_x2a_Sl_soilw, ig)
          dust_ndx = dust_idx1()
          ! check that dust constituents are actually in the simulation
          if (dust_ndx>0) then
#if (defined MODAL_AERO)
#if (defined MODAL_AERO_7MODE)
            cam_in(c)%cflx(i,dust_ndx   )  = 0.13_r8  &  ! 1st mode, based on Zender et al (2003) Table 1
#elif (defined MODAL_AERO_3MODE)
            cam_in(c)%cflx(i,dust_ndx   )  = 0.032_r8  &  ! 1st mode, based on Zender et al (2003) Table 1
#endif
                                           * (-fptr(index_x2a_Fall_flxdst1, ig) &
                                              -fptr(index_x2a_Fall_flxdst2, ig) &
                                              -fptr(index_x2a_Fall_flxdst3, ig) &
                                              -fptr(index_x2a_Fall_flxdst4, ig))
#if (defined MODAL_AERO_7MODE)
            cam_in(c)%cflx(i,dust_ndx-spc_ndx(1)+spc_ndx(2))  = 0.87_r8 &  ! 2nd mode
#elif (defined MODAL_AERO_3MODE)
            cam_in(c)%cflx(i,dust_ndx-spc_ndx(1)+spc_ndx(2))  = 0.968_r8 &  ! 2nd mode
#endif
                                           * (-fptr(index_x2a_Fall_flxdst1, ig) &
                                              -fptr(index_x2a_Fall_flxdst2, ig) &
                                              -fptr(index_x2a_Fall_flxdst3, ig) &
                                              -fptr(index_x2a_Fall_flxdst4, ig))
#else
	    cam_in(c)%cflx(i,dust_ndx   )  = -fptr(index_x2a_Fall_flxdst1, ig)
	    cam_in(c)%cflx(i,dust_ndx +1)  = -fptr(index_x2a_Fall_flxdst2, ig)
	    cam_in(c)%cflx(i,dust_ndx +2)  = -fptr(index_x2a_Fall_flxdst3, ig)
	    cam_in(c)%cflx(i,dust_ndx +3)  = -fptr(index_x2a_Fall_flxdst4, ig)
#endif
          endif

          ! dry dep velocities
          if ( index_x2a_Sl_ddvel/=0 .and. n_drydep>0 ) then
             cam_in(c)%depvel(i,:n_drydep) = &
                  fptr(index_x2a_Sl_ddvel:index_x2a_Sl_ddvel+n_drydep-1, ig)
          endif
          !
          ! fields needed to calculate water isotopes to ocean evaporation processes
          !
          cam_in(c)%ustar(i) = fptr(index_x2a_So_ustar,ig)
          cam_in(c)%re(i)    = fptr(index_x2a_So_re   ,ig)
          cam_in(c)%ssq(i)   = fptr(index_x2a_So_ssq  ,ig)
          !
          ! bgc scenarios
          !
          if (index_x2a_Fall_fco2_lnd /= 0) then
             cam_in(c)%fco2_lnd(i) = -fptr(index_x2a_Fall_fco2_lnd,ig)
          end if
          if (index_x2a_Faoo_fco2_ocn /= 0) then
             cam_in(c)%fco2_ocn(i) = -fptr(index_x2a_Faoo_fco2_ocn,ig)
          end if
          if (index_x2a_Faoo_fdms_ocn /= 0) then
             cam_in(c)%fdms(i)     = -fptr(index_x2a_Faoo_fdms_ocn,ig)
          end if

          ig=ig+1

       end do
    end do

    ! Get total co2 flux from components,
    ! Note - co2_transport determines if cam_in(c)%cflx(i,c_i(1:4)) is allocated

    if (co2_transport()) then

       ! Interpolate in time for flux data read in
       if (co2_readFlux_ocn) then
          call co2_time_interp_ocn
       end if
       if (co2_readFlux_fuel) then
          call co2_time_interp_fuel
       end if
       
       ! from ocn : data read in or from coupler or zero
       ! from fuel: data read in or zero
       ! from lnd : through coupler or zero
       do c=begchunk,endchunk
          ncols = get_ncols_p(c)                                                 
          do i=1,ncols                                                               
             
             ! all co2 fluxes in unit kgCO2/m2/s ! co2 flux from ocn 
             if (index_x2a_Faoo_fco2_ocn /= 0) then
                cam_in(c)%cflx(i,c_i(1)) = cam_in(c)%fco2_ocn(i)
             else if (co2_readFlux_ocn) then 
                ! convert from molesCO2/m2/s to kgCO2/m2/s
                cam_in(c)%cflx(i,c_i(1)) = &
                     -data_flux_ocn%co2flx(i,c)*(1._r8- cam_in(c)%landfrac(i)) &
                     *mwco2*1.0e-3_r8
             else
                cam_in(c)%cflx(i,c_i(1)) = 0._r8
             end if
             
             ! co2 flux from fossil fuel
             if (co2_readFlux_fuel) then
                cam_in(c)%cflx(i,c_i(2)) = data_flux_fuel%co2flx(i,c)
             else
                cam_in(c)%cflx(i,c_i(2)) = 0._r8
             end if
             
             ! co2 flux from land (cpl already multiplies flux by land fraction)
             if (index_x2a_Fall_fco2_lnd /= 0) then
                cam_in(c)%cflx(i,c_i(3)) = cam_in(c)%fco2_lnd(i)
             else
                cam_in(c)%cflx(i,c_i(3)) = 0._r8
             end if
             
             ! merged co2 flux
             cam_in(c)%cflx(i,c_i(4)) = cam_in(c)%cflx(i,c_i(1)) + &
                                        cam_in(c)%cflx(i,c_i(2)) + &
                                        cam_in(c)%cflx(i,c_i(3))
          end do
       end do
    end if
    !
    ! if first step, determine longwave up flux from the surface temperature 
    !
    if (first_time) then
       if (is_first_step()) then
          do c=begchunk, endchunk
             ncols = get_ncols_p(c)
             do i=1,ncols
                cam_in(c)%lwup(i) = shr_const_stebol*(cam_in(c)%ts(i)**4)
             end do
          end do
       end if
       first_time = .false.
    end if

  end subroutine atm_import_esmf

!===============================================================================

  subroutine atm_export_esmf( cam_out, a2x_a, rc )

    implicit none
    !-------------------------------------------------------------------
    !
    ! Arguments
    !
    type(cam_out_t),   intent(in)    :: cam_out(begchunk:endchunk) 
    type(ESMF_Array) , intent(inout) :: a2x_a
    integer, intent(out)             :: rc
    !
    ! Local variables
    !
    integer :: avsize, avnat
    integer :: i,m,c,n,ig       ! indices
    integer :: ncols            ! Number of columns
    real(R8), pointer :: fptr(:,:)        ! pointer into    array data
    !-----------------------------------------------------------------------

    ! Copy from component Arrays into chunk array data structure
    ! Rearrange data from chunk structure into lat-lon buffer and subsequently
    ! create attribute vector

    rc = ESMF_SUCCESS
    call ESMF_ArrayGet(a2x_a, localDe=0, farrayPtr=fptr, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ig=1
    do c=begchunk, endchunk
       ncols = get_ncols_p(c)
       do i=1,ncols
          fptr(index_a2x_Sa_pslv   ,ig) = cam_out(c)%psl(i)
          fptr(index_a2x_Sa_z      ,ig) = cam_out(c)%zbot(i)   
          fptr(index_a2x_Sa_u      ,ig) = cam_out(c)%ubot(i)   
          fptr(index_a2x_Sa_v      ,ig) = cam_out(c)%vbot(i)   
          fptr(index_a2x_Sa_tbot   ,ig) = cam_out(c)%tbot(i)   
          fptr(index_a2x_Sa_ptem   ,ig) = cam_out(c)%thbot(i)  
          fptr(index_a2x_Sa_pbot   ,ig) = cam_out(c)%pbot(i)   
          fptr(index_a2x_Sa_shum   ,ig) = cam_out(c)%qbot(i,1) 
	      fptr(index_a2x_Sa_dens   ,ig) = cam_out(c)%rho(i)
          fptr(index_a2x_Faxa_swnet,ig) = cam_out(c)%netsw(i)      
          fptr(index_a2x_Faxa_lwdn ,ig) = cam_out(c)%flwds(i)  
          fptr(index_a2x_Faxa_rainc,ig) = (cam_out(c)%precc(i)-cam_out(c)%precsc(i))*1000._r8
          fptr(index_a2x_Faxa_rainl,ig) = (cam_out(c)%precl(i)-cam_out(c)%precsl(i))*1000._r8
          fptr(index_a2x_Faxa_snowc,ig) = cam_out(c)%precsc(i)*1000._r8
          fptr(index_a2x_Faxa_snowl,ig) = cam_out(c)%precsl(i)*1000._r8
          fptr(index_a2x_Faxa_swndr,ig) = cam_out(c)%soll(i)   
          fptr(index_a2x_Faxa_swvdr,ig) = cam_out(c)%sols(i)   
          fptr(index_a2x_Faxa_swndf,ig) = cam_out(c)%solld(i)  
          fptr(index_a2x_Faxa_swvdf,ig) = cam_out(c)%solsd(i)  

          ! aerosol deposition fluxes
          fptr(index_a2x_Faxa_bcphidry,ig) = cam_out(c)%bcphidry(i)
          fptr(index_a2x_Faxa_bcphodry,ig) = cam_out(c)%bcphodry(i)
          fptr(index_a2x_Faxa_bcphiwet,ig) = cam_out(c)%bcphiwet(i)
          fptr(index_a2x_Faxa_ocphidry,ig) = cam_out(c)%ocphidry(i)
          fptr(index_a2x_Faxa_ocphodry,ig) = cam_out(c)%ocphodry(i)
          fptr(index_a2x_Faxa_ocphiwet,ig) = cam_out(c)%ocphiwet(i)
          fptr(index_a2x_Faxa_dstwet1,ig)  = cam_out(c)%dstwet1(i)
          fptr(index_a2x_Faxa_dstdry1,ig)  = cam_out(c)%dstdry1(i)
          fptr(index_a2x_Faxa_dstwet2,ig)  = cam_out(c)%dstwet2(i)
          fptr(index_a2x_Faxa_dstdry2,ig)  = cam_out(c)%dstdry2(i)
          fptr(index_a2x_Faxa_dstwet3,ig)  = cam_out(c)%dstwet3(i)
          fptr(index_a2x_Faxa_dstdry3,ig)  = cam_out(c)%dstdry3(i)
          fptr(index_a2x_Faxa_dstwet4,ig)  = cam_out(c)%dstwet4(i)
          fptr(index_a2x_Faxa_dstdry4,ig)  = cam_out(c)%dstdry4(i)

          if (index_a2x_Sa_co2prog /= 0) then
             fptr(index_a2x_Sa_co2prog,ig) = cam_out(c)%co2prog(i) ! atm prognostic co2
          end if
          if (index_a2x_Sa_co2diag /= 0) then
             fptr(index_a2x_Sa_co2diag,ig) = cam_out(c)%co2diag(i) ! atm diagnostic co2
          end if

          ig=ig+1
       end do
    end do
 
  end subroutine atm_export_esmf

!===============================================================================
  subroutine atm_domain_esmf( dom, rc)

    implicit none
    !-------------------------------------------------------------------
    !
    ! Arguments
    !
    type(ESMF_Array), intent(inout)     :: dom
    integer, intent(out)                :: rc
    !
    ! Local Variables
    !
    integer  :: n,i,c,ncols           ! indices	
    real(r8) :: lats(pcols)           ! array of chunk latitudes
    real(r8) :: lons(pcols)           ! array of chunk longitude
    real(r8) :: area(pcols)           ! area in radians squared for each grid point
    real(r8), parameter:: radtodeg = 180.0_r8/SHR_CONST_PI

    real(R8),    pointer    :: fptr (:,:)
    integer :: klon,klat,karea,kmask,kfrac ! domain fields

    !-------------------------------------------------------------------
    rc = ESMF_SUCCESS
    !-------------------------------------------------------------------

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

    fptr(:,:) = -9999.0_R8
    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, ncols, lats)
       call get_rlon_all_p(c, ncols, lons)
       call get_area_all_p(c, ncols, area)
       do i=1,ncols
          n = n+1
          fptr(klat,n) = lats(i)*radtodeg
          fptr(klon,n) = lons(i)*radtodeg
          fptr(karea,n) = area(i) 
          fptr(kmask,n) = 1.0_r8
          fptr(kfrac,n) = 1.0_r8
       end do
    end do

  end subroutine atm_domain_esmf

!===============================================================================

!
!===========================================================================================
!
  subroutine atm_read_srfrest_esmf( EClock, dom, x2a, a2x)
    use cam_pio_utils
    implicit none
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(ESMF_Clock),intent(in)    :: EClock
    type(ESMF_Array),intent(inout) :: dom
    type(ESMF_Array),intent(inout) :: x2a
    type(ESMF_Array),intent(inout) :: a2x
    ! 
    ! Local variables
    !
    integer         :: rcode,rc        ! return error code
    integer         :: yr_spec      ! Current year
    integer         :: mon_spec     ! Current month
    integer         :: day_spec     ! Current day
    integer         :: sec_spec     ! Current time of day (sec)
    real(R8),    pointer    :: fptr (:,:)
    integer, pointer :: dof(:)
    integer :: lnx, nf_x2a, nf_a2x, k
    real(r8), allocatable :: tmp(:)
    type(file_desc_t) :: file
    type(io_desc_t) :: iodesc
    type(var_desc_t) :: varid
    character(CL)    :: itemc       ! string converted to char

    !-----------------------------------------------------------------------

    ! Determine and open surface restart dataset

    
    call seq_timemgr_EClockGetData( EClock, curr_yr=yr_spec,curr_mon=mon_spec, &
         curr_day=day_spec, curr_tod=sec_spec ) 
    fname_srf_cam = interpret_filename_spec( rsfilename_spec_cam, case=get_restcase(), &
         yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )
    pname_srf_cam = trim(get_restartdir() )//fname_srf_cam
    call getfil(pname_srf_cam, fname_srf_cam)
    
    call cam_pio_openfile(File, fname_srf_cam, 0)

    call esmf2mct_init(dom,dof,rc)
    lnx = ngcols
    call pio_initdecomp(pio_subsystem, pio_double, (/lnx/), dof, iodesc)
    allocate(tmp(size(dof)))
    deallocate(dof)
    
    call esmfshr_util_ArrayGetSize(x2a,lsize1=nf_x2a)
    call esmfshr_util_ArrayGetSize(a2x,lsize1=nf_a2x)

    call ESMF_ArrayGet(x2a, localDe=0, farrayPtr=fptr, rc=rc)
    do k=1,nf_x2a
       call esmfshr_util_ArrayGetName(x2a,k,itemc)

       call pio_seterrorhandling(File, pio_bcast_error)
       rcode = pio_inq_varid(File,'x2a_'//trim(itemc) ,varid)
       if (rcode == pio_noerr) then
          call pio_read_darray(File, varid, iodesc, tmp, rcode)
          fptr(k,:) = tmp(:)
       else
         if (masterproc) then
             write(iulog,*)'srfrest warning: field ',trim(itemc),' is not on restart file'
             write(iulog,*)'for backwards compatibility will set it to 0'
          end if
          fptr(k,:) = 0._r8
       end if
       call pio_seterrorhandling(File, pio_internal_error)
    end do

    call ESMF_ArrayGet(a2x, localDe=0, farrayPtr=fptr, rc=rc)
    do k=1,nf_a2x
       call esmfshr_util_ArrayGetName(a2x,k,itemc)

       rcode = pio_inq_varid(File,'a2x_'//trim(itemc) ,varid)
       call pio_read_darray(File, varid, iodesc, tmp, rcode)
       fptr(k,:) = tmp(:)
    end do

    call pio_freedecomp(File,iodesc)
    call pio_closefile(File)
    deallocate(tmp)

  end subroutine atm_read_srfrest_esmf

!===========================================================================================
!
  subroutine atm_write_srfrest_esmf( dom, x2a, a2x, &
       yr_spec, mon_spec, day_spec, sec_spec)

    use cam_pio_utils
    implicit none

    !-----------------------------------------------------------------------

    ! Arguments

    type(ESMF_Array),intent(inout) :: dom
    type(ESMF_Array),intent(inout) :: x2a
    type(ESMF_Array),intent(inout) :: a2x
    integer        , intent(in) :: yr_spec         ! Simulation year
    integer        , intent(in) :: mon_spec        ! Simulation month
    integer        , intent(in) :: day_spec        ! Simulation day
    integer        , intent(in) :: sec_spec        ! Seconds into current simulation day

    ! Local variables

    integer         :: rcode,rc        ! return error code
    integer, pointer :: dof(:)
    integer :: nf_x2a, nf_a2x, lnx, dimid(1), k
    real(R8),    pointer    :: fptr (:,:)
    type(file_desc_t) :: file
    type(var_desc_t), pointer :: varid_x2a(:), varid_a2x(:)
    type(io_desc_t)  :: iodesc
    character(CL)    :: itemc       ! string converted to char

    !-----------------------------------------------------------------------

    ! Determine and open surface restart dataset

    fname_srf_cam = interpret_filename_spec( rsfilename_spec_cam, &
         yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )
    call cam_pio_createfile(File, fname_srf_cam, 0)

    call esmf2mct_init(dom,dof,rc)
    lnx = ngcols
    call pio_initdecomp(pio_subsystem, pio_double, (/lnx/), dof, iodesc)
    deallocate(dof)
    
    call esmfshr_util_ArrayGetSize(x2a,lsize1=nf_x2a)
    call esmfshr_util_ArrayGetSize(a2x,lsize1=nf_a2x)
    allocate(varid_x2a(nf_x2a))
    allocate(varid_a2x(nf_a2x))

    rcode = pio_def_dim(File,'x2a_nx',lnx,dimid(1))
    do k = 1,nf_x2a
       call esmfshr_util_ArrayGetName(x2a,k,itemc)
       rcode = pio_def_var(File,'x2a_'//trim(itemc),PIO_DOUBLE,dimid,varid_x2a(k))
       rcode = pio_put_att(File,varid_x2a(k),"_fillvalue",fillvalue)
    enddo

    rcode = pio_def_dim(File,'a2x_nx',lnx,dimid(1))
    do k = 1,nf_a2x
       call esmfshr_util_ArrayGetName(a2x,k,itemc)
       rcode = PIO_def_var(File,'a2x_'//trim(itemc),PIO_DOUBLE,dimid,varid_a2x(k))
       rcode = PIO_put_att(File,varid_a2x(k),"_fillvalue",fillvalue)
    enddo

    rcode = pio_enddef(File)  ! don't check return code, might be enddef already

    call ESMF_ArrayGet(x2a, localDe=0, farrayPtr=fptr, rc=rc)
    do k=1,nf_x2a
       call pio_write_darray(File, varid_x2a(k), iodesc, fptr(k,:), rcode)
    end do

    call ESMF_ArrayGet(a2x, localDe=0, farrayPtr=fptr, rc=rc)
    do k=1,nf_a2x
       call pio_write_darray(File, varid_a2x(k), iodesc, fptr(k,:), rcode)       
    end do

    deallocate(varid_x2a, varid_a2x)

    call pio_freedecomp(File,iodesc)
    call pio_closefile(file)

  end subroutine atm_write_srfrest_esmf

!================================================================================

end module atm_comp_esmf
