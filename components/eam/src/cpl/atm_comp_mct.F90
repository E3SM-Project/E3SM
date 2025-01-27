module atm_comp_mct

  use pio              , only: file_desc_t, io_desc_t, var_desc_t, pio_double, pio_def_dim, &
                               pio_put_att, pio_enddef, pio_initdecomp, pio_read_darray, pio_freedecomp, &
                               pio_closefile, pio_write_darray, pio_def_var, pio_inq_varid, &
	                       pio_noerr, pio_bcast_error, pio_internal_error, pio_seterrorhandling 
  use mct_mod
  use seq_comm_mct     , only: info_taskmap_comp
  use seq_cdata_mod
  use esmf

  use seq_flds_mod ! for seq_flds_x2a_fields, seq_flds_dom_fields, etc
  use seq_infodata_mod
  use seq_timemgr_mod

  use shr_kind_mod     , only: r8 => shr_kind_r8, cl=>shr_kind_cl, cxx=>shr_kind_cxx
  use shr_kind_mod     , only: cs => shr_kind_cs
  use shr_file_mod     , only: shr_file_getunit, shr_file_freeunit, &
                               shr_file_setLogUnit, shr_file_setLogLevel, &
                               shr_file_getLogUnit, shr_file_getLogLevel, &
		               shr_file_setIO
  use shr_sys_mod      , only: shr_sys_flush, shr_sys_abort
  use shr_taskmap_mod  , only: shr_taskmap_write

  use cam_cpl_indices
  ! it has atm_import, atm_export 
  use atm_import_export
  !  atm_export_moab is private here, atm_import_moab too

  use cam_comp
  use cam_instance     , only: cam_instance_init, inst_index, inst_suffix
  use cam_control_mod  , only: nsrest, aqua_planet, eccen, obliqr, lambm0, mvelpp
  use radiation        , only: radiation_do, radiation_nextsw_cday
  use phys_grid        , only: get_ncols_p, ngcols, get_gcol_p, get_rlat_all_p, &
	                       get_rlon_all_p, get_area_all_p
  use ppgrid           , only: pcols, begchunk, endchunk       
  use dyn_grid         , only: get_horiz_grid_dim_d
  use camsrfexch       , only: cam_out_t, cam_in_t     
  use cam_restart      , only: get_restcase, get_restartdir
  use cam_history      , only: outfld, ctitle
  use cam_abortutils       , only: endrun
  use filenames        , only: interpret_filename_spec, caseid, brnch_retain_casename, &
                               hostname, username, version
#ifdef SPMD
  use spmd_utils       , only: spmdinit, masterproc, iam, npes, nsmps, &
                               proc_smp_map
  use mpishorthand     , only: mpicom
#else
  use spmd_utils       , only: spmdinit, masterproc, mpicom, iam, npes, nsmps, &
                               proc_smp_map
#endif
  use time_manager     , only: get_curr_calday, advance_timestep, get_curr_date, get_nstep, &
                               get_step_size, timemgr_init, timemgr_check_restart
  use ioFileMod
  use perf_mod
  use cam_logfile      , only: iulog
  use co2_cycle        , only: co2_readFlux_ocn, co2_readFlux_fuel
  use runtime_opts     , only: read_namelist

  use iop_data_mod     , only: single_column,scmlat,scmlon,scm_multcols
  use lnd_infodata     , only: precip_downscaling_method !Precipitation downscaling method used in the land model

#ifdef HAVE_MOAB
  use seq_comm_mct     , only: mphaid ! atm physics grid id in MOAB, on atm pes
  use iso_c_binding 
  use seq_comm_mct,     only : num_moab_exports
#ifdef MOABCOMP
  use seq_comm_mct, only:  seq_comm_compare_mb_mct
#endif
#endif


!
! !PUBLIC TYPES:
  implicit none
  save
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: atm_init_mct
  public :: atm_run_mct
  public :: atm_final_mct

!--------------------------------------------------------------------------
! Private interfaces
!--------------------------------------------------------------------------

  private :: atm_setgsmap_mct
  private :: atm_domain_mct
  private :: atm_read_srfrest_mct
  private :: atm_write_srfrest_mct
#ifdef HAVE_MOAB
  private :: atm_read_srfrest_moab
  private :: atm_write_srfrest_moab
#endif

!--------------------------------------------------------------------------
! Private data
!--------------------------------------------------------------------------

  type(cam_in_t) , pointer :: cam_in(:)
  type(cam_out_t), pointer :: cam_out(:)

  integer, parameter  :: nlen = 256     ! Length of character strings
  character(len=nlen) :: fname_srf_cam  ! surface restart filename
  character(len=nlen) :: pname_srf_cam  ! surface restart full pathname
#ifdef HAVE_MOAB
  character(len=nlen) :: moab_fname_srf_cam  ! surface restart filename
  character(len=nlen) :: moab_pname_srf_cam  ! surface restart full pathname
#endif
  ! Filename specifier for restart surface file
  character(len=cl) :: rsfilename_spec_cam

  ! Are all surface types present   
  logical :: lnd_present ! if true => land is present
  logical :: ocn_present ! if true => ocean is present

  integer,                 pointer :: dof(:) ! needed for pio_init decomp for restarts
  type(seq_infodata_type), pointer :: infodata

#ifdef HAVE_MOAB
  ! to store all fields to be set in moab 
  integer , private :: mblsize, totalmbls, nsend, totalmbls_r, nrecv
  real(r8) , allocatable, private :: a2x_am(:,:) ! atm to coupler, on atm mesh, on atm component pes
  real(r8) , allocatable, private :: x2a_am(:,:) ! coupler to atm, on atm mesh, on atm component pes
  integer,                pointer :: global_ids(:) ! they could be dof(), but better maintain our own list

#ifdef MOABCOMP
  integer  :: mpicom_atm_moab ! used just for mpi-reducing the difference between moab tags and mct avs
  integer :: rank2
#endif

#endif 

!================================================================================
CONTAINS
!================================================================================

  subroutine atm_init_mct( EClock, cdata_a, x2a_a, a2x_a, NLFilename )

    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    type(ESMF_Clock),intent(inout)              :: EClock
    type(seq_cdata), intent(inout)              :: cdata_a
    type(mct_aVect), intent(inout)              :: x2a_a
    type(mct_aVect), intent(inout)              :: a2x_a   
    character(len=*), optional,   intent(IN)    :: NLFilename ! Namelist filename
    !
    ! Locals
    !
    type(mct_gsMap), pointer   :: gsMap_atm
    type(mct_gGrid), pointer   :: dom_a
    integer :: ATMID
    integer :: mpicom_atm
    logical :: no_taskmap_output ! if true, do not write out task-to-node mapping
    logical :: verbose_taskmap_output ! if true, use verbose task-to-node mapping format
    integer :: lsize
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
    integer :: shrlogunit,shrloglev ! old values
    logical :: first_time = .true.
    character(len=cs) :: calendar      ! Calendar type
    character(len=cs) :: starttype     ! infodata start type
    character(len=8)           :: c_inst_index  ! instance number
    character(len=8)           :: c_npes        ! number of pes
    integer :: lbnum
    integer :: hdim1_d, hdim2_d ! dimensions of rectangular horizontal grid
                                ! data structure, If 1D data structure, then
                                ! hdim2_d == 1.
    character(len=64) :: filein ! Input namelist filename

#ifdef MOABCOMP
    real(r8)                 :: difference
    type(mct_list) :: temp_list
    integer :: size_list, index_list, ent_type
    type(mct_string)    :: mctOStr  !
    character(CXX) ::tagname, mct_field, modelStr
#endif 

    !-----------------------------------------------------------------------
    !
    ! Determine cdata points
    !
#if (defined _MEMTRACE)
    if(masterproc) then
      lbnum=1
      call memmon_dump_fort('memmon.out','atm_init_mct:start::',lbnum)
    endif                      
#endif                         

    call seq_cdata_setptrs(cdata_a, ID=ATMID, mpicom=mpicom_atm, &
         gsMap=gsMap_atm, dom=dom_a, infodata=infodata)

#ifdef MOABCOMP
    mpicom_atm_moab = mpicom_atm ! just store it now, for later use
    call shr_mpi_commrank( mpicom_atm_moab, rank2 )
#endif 

    if (first_time) then
       
       call cam_instance_init(ATMID)

       ! Set filename specifier for restart surface file
       ! (%c=caseid, $y=year, $m=month, $d=day, $s=seconds in day)
       rsfilename_spec_cam = '%c.eam' // trim(inst_suffix) // '.rs.%y-%m-%d-%s.nc'

       ! Determine attribute vector indices

       call cam_cpl_indices_set()

       ! Initialize MPI for CAM

       call spmdinit(mpicom_atm, calc_proc_smp_map=.false.)
       
       ! Redirect share output to cam log
       
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

       ! Identify SMP nodes and process/SMP mapping for this instance
       ! (Assume that processor names are SMP node names on SMP clusters.)
       write(c_inst_index,'(i8)') inst_index

       if (info_taskmap_comp > 0) then

          no_taskmap_output = .false.

          if (info_taskmap_comp == 1) then
             verbose_taskmap_output = .false.
          else
             verbose_taskmap_output = .true.
          endif

          write(c_npes,'(i8)') npes

          if (masterproc) then
             write(iulog,'(/,3A)') &
                trim(adjustl(c_npes)), &
                ' pes participating in computation of CAM instance #', &
                trim(adjustl(c_inst_index))
             call shr_sys_flush(iulog)
          endif

       else

          no_taskmap_output = .true.
          verbose_taskmap_output = .false.

       endif

       call t_startf('shr_taskmap_write')
       call shr_taskmap_write(iulog, mpicom_atm,                    &
                              'ATM #'//trim(adjustl(c_inst_index)), &
                              verbose=verbose_taskmap_output,       &
                              no_output=no_taskmap_output,          &
                              save_nnodes=nsmps,                    &
                              save_task_node_map=proc_smp_map       )
       call shr_sys_flush(iulog)
       call t_stopf('shr_taskmap_write')

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
       call seq_infodata_GetData( infodata,                                           &
            case_name=caseid, case_desc=ctitle,                                       &
            start_type=starttype,                                                     &
            aqua_planet=aqua_planet,                                                  &
            brnch_retain_casename=brnch_retain_casename,                              &
            hostname=hostname, username=username, model_version=version,              &
            single_column=single_column, scmlat=scmlat, scmlon=scmlon,                &
            scm_multcols=scm_multcols,                                                &
            orb_eccen=eccen, orb_mvelpp=mvelpp, orb_lambm0=lambm0, orb_obliqr=obliqr, &
            lnd_present=lnd_present, ocn_present=ocn_present,                         &
            perpetual=perpetual_run, perpetual_ymd=perpetual_ymd)
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
          write(iulog,*) 'atm_comp_mct: ERROR: unknown starttype'
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
       call t_startf('read_namelist')
       filein = "atm_in" // trim(inst_suffix)
       call read_namelist(single_column_in=single_column, scmlat_in=scmlat, &
            scmlon_in=scmlon, scm_multcols_in=scm_multcols, nlfilename_in=filein)
       call t_stopf('read_namelist')
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
       call t_startf('cam_init')
       call cam_init( cam_out, cam_in, mpicom_atm, &
            start_ymd, start_tod, ref_ymd, ref_tod, stop_ymd, stop_tod, &
            perpetual_run, perpetual_ymd, calendar)
       call t_stopf('cam_init')
       !
       ! Check consistency of restart time information with input clock
       !
       if (nsrest /= 0) then
          dtime_cam = get_step_size()
          call timemgr_check_restart( calendar, start_ymd, start_tod, ref_ymd, &
               ref_tod, dtime_cam, perpetual_run, perpetual_ymd)
       end if
       !
       ! Initialize MCT gsMap, domain and attribute vectors (and dof)
       !
       call atm_SetgsMap_mct( mpicom_atm, ATMID, gsMap_atm )
       lsize = mct_gsMap_lsize(gsMap_atm, mpicom_atm)

       ! Set dof (module variable, needed for pio for restarts)
       call mct_gsmap_orderedpoints(gsmap_atm, iam, dof)
       !
       ! Initialize MCT domain and add data
       !
       call atm_domain_mct( lsize, gsMap_atm, dom_a )

       ! Initialize MCT attribute vectors
       !
       call mct_aVect_init(a2x_a, rList=seq_flds_a2x_fields, lsize=lsize)
       call mct_aVect_zero(a2x_a)
       
       call mct_aVect_init(x2a_a, rList=seq_flds_x2a_fields, lsize=lsize) 
       call mct_aVect_zero(x2a_a)
       !
       ! Create initial atm export state
       !
       call atm_export( cam_out, a2x_a%rattr )

#ifdef HAVE_MOAB
       ! Initialize MOAB physgrid mesh and add coordinate,mask data
       ! NOTE:  dynamics mesh is initialized in homme source
       ! as part of cam_init above
       call init_moab_atm_phys( cdata_a )
       ! initialize arrays to hold data for MOAB
       mblsize = lsize
       nsend = mct_avect_nRattr(a2x_a)
       totalmbls = mblsize * nsend   ! size of the double array
       allocate (a2x_am(mblsize, nsend) )

       nrecv = mct_avect_nRattr(x2a_a)
       totalmbls_r = mblsize * nrecv   ! size of the double array used to receive
       allocate (x2a_am(mblsize, nrecv) ) ! these will be received by moab tags, then used to set cam in surf data
       !
       ! Create initial atm export state inside moab
       !
       call atm_export_moab(Eclock, cam_out )
#endif
       !
       ! Set flag to specify that an extra albedo calculation is to be done (i.e. specify active)
       !
       call seq_infodata_PutData(infodata, atm_prognostic=.true.)
       call get_horiz_grid_dim_d(hdim1_d, hdim2_d)
       call seq_infodata_PutData(infodata, atm_nx=hdim1_d, atm_ny=hdim2_d)

       ! Set flag to indicate that CAM will provide carbon and dust deposition fluxes.
       ! This is now hardcoded to .true. since the ability of CICE to read these
       ! fluxes from a file has been removed.
       call seq_infodata_PutData(infodata, atm_aero=.true.)

       !
       ! Set time step of radiation computation as the current calday
       ! This will only be used on the first timestep of an initial run
       !
       if (nsrest == 0) then
          nextsw_cday = get_curr_calday()
          call seq_infodata_PutData( infodata, nextsw_cday=nextsw_cday )
       end if
       
       ! End redirection of share output to cam log
       
       call shr_file_setLogUnit (shrlogunit)
       call shr_file_setLogLevel(shrloglev)

       first_time = .false.

    else ! so here first_time == .false.
       
       ! For initial run, run cam radiation/clouds and return
       ! For restart run, read restart x2a_a
       ! Note - a2x_a is computed upon the completion of the previous run - cam_run1 is called
       ! only for the purposes of finishing the flux averaged calculation to compute a2x_a
       ! Note - cam_run1 is called on restart only to have cam internal state consistent with the 
       ! a2x_a state sent to the coupler

       !Obtain the precipitation downscaling method from the land model
       call seq_infodata_GetData( infodata,                                           &
            precip_downscaling_method=precip_downscaling_method )


       ! Redirect share output to cam log

       call shr_file_getLogUnit (shrlogunit)
       call shr_file_getLogLevel(shrloglev)
       call shr_file_setLogUnit (iulog)

       call seq_timemgr_EClockGetData(EClock,curr_ymd=CurrentYMD, StepNo=StepNo, dtime=DTime_Sync )
       if (StepNo == 0) then
#ifdef MOABCOMP
         ! loop over all fields in seq_flds_x2a_fields
          call mct_list_init(temp_list ,seq_flds_x2a_fields)
          size_list=mct_list_nitem (temp_list)
          ent_type = 0 ! entity type is vertex for phys atm
          if (rank2 .eq. 0) print *, num_moab_exports, trim(seq_flds_x2a_fields), ' atm import check'
          modelStr='atm init2'
          do index_list = 1, size_list
            call mct_list_get(mctOStr,index_list,temp_list)
            mct_field = mct_string_toChar(mctOStr)
            tagname= trim(mct_field)//C_NULL_CHAR
            call seq_comm_compare_mb_mct(modelStr, mpicom_atm_moab, x2a_a, mct_field,  mphaid, tagname, ent_type, difference)
          enddo
          call mct_list_clean(temp_list)

#endif
       ! so the cam import is before moab    
          call atm_import( x2a_a%rattr, cam_in )
#ifdef HAVE_MOAB
       ! move moab import after cam import, so moab takes precedence
          call atm_import_moab(Eclock, cam_in)
#endif    



          call t_startf('CAM_run1')
          call cam_run1 ( cam_in, cam_out ) 
          call t_stopf('CAM_run1')
        
          call atm_export( cam_out, a2x_a%rattr )
#ifdef HAVE_MOAB
          call atm_export_moab(Eclock, cam_out)
#endif   
       else ! if (StepNo != 0) then

          call t_startf('atm_read_srfrest_mct')
          call atm_read_srfrest_mct( EClock, x2a_a, a2x_a )
          call t_stopf('atm_read_srfrest_mct')
#ifdef HAVE_MOAB
          call atm_read_srfrest_moab ( EClock )
#endif

          ! Sent .true. as an optional argument so that restart_init is set to .true.  in atm_import
	      ! This will ensure BFB restarts whenever qneg4 updates fluxes on the restart time step
          call atm_import( x2a_a%rattr, cam_in, .true. )
#ifdef HAVE_MOAB
          call atm_import_moab(Eclock, cam_in, .true. )
#endif   

          call t_startf('cam_run1')
          call cam_run1 ( cam_in, cam_out ) 
          call t_stopf('cam_run1')
       end if

       ! Compute time of next radiation computation, like in run method for exact restart

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
       call seq_infodata_PutData( infodata, nextsw_cday=nextsw_cday ) 

       ! End redirection of share output to cam log
       
       call shr_file_setLogUnit (shrlogunit)
       call shr_file_setLogLevel(shrloglev)
       
    end if

#if (defined _MEMTRACE )
    if(masterproc) then
      lbnum=1
      call memmon_dump_fort('memmon.out','atm_init_mct:end::',lbnum)
      call memmon_reset_addr()
    endif
#endif
    
    call shr_sys_flush(iulog)

 end subroutine atm_init_mct

 !================================================================================

 subroutine atm_run_mct( EClock, cdata_a, x2a_a, a2x_a)

    !-----------------------------------------------------------------------
    use time_manager,    only: advance_timestep, get_curr_date, get_curr_calday, &
	                       get_step_size
   !use iop,             only: scam_use_iop_srf
    use pmgrid,          only: plev, plevp
    use constituents,    only: pcnst
    use shr_sys_mod,     only: shr_sys_flush

    ! 
    ! Arguments
    !
    type(ESMF_Clock)            ,intent(inout) :: EClock
    type(seq_cdata)             ,intent(inout) :: cdata_a
    type(mct_aVect)             ,intent(inout) :: x2a_a
    type(mct_aVect)             ,intent(inout) :: a2x_a
    !
    ! Local variables
    !
    integer :: lsize           ! size of attribute vector
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
    integer :: lbnum
    character(len=*), parameter :: subname="atm_run_mct"
    !-----------------------------------------------------------------------
#ifdef MOABCOMP
    real(r8)                 :: difference
    type(mct_list) :: temp_list
    integer :: size_list, index_list, ent_type
    type(mct_string)    :: mctOStr  !
    character(CXX) ::tagname, mct_field, modelStr
#endif 

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

    !load orbital parameters 

    call seq_infodata_GetData( infodata,                                           &
       orb_eccen=eccen, orb_mvelpp=mvelpp, orb_lambm0=lambm0, orb_obliqr=obliqr)

    nlend_sync = seq_timemgr_StopAlarmIsOn(EClock)
    rstwr_sync = seq_timemgr_RestartAlarmIsOn(EClock)

    ! Map input from mct to cam data structure

    call t_startf ('CAM_import')
! move moab import after regular atm import, so it would be in charge
    call atm_import( x2a_a%rattr, cam_in )
#ifdef HAVE_MOAB

#ifdef MOABCOMP
    ! loop over all fields in seq_flds_x2a_fields
    call mct_list_init(temp_list ,seq_flds_x2a_fields)
    size_list=mct_list_nitem (temp_list)
    ent_type = 0 ! entity type is vertex for phys atm
    if (rank2 .eq. 0) print *, num_moab_exports, trim(seq_flds_x2a_fields)
    modelStr ='atm run'
    do index_list = 1, size_list
      call mct_list_get(mctOStr,index_list,temp_list)
      mct_field = mct_string_toChar(mctOStr)
      tagname= trim(mct_field)//C_NULL_CHAR
      call seq_comm_compare_mb_mct(modelStr, mpicom_atm_moab, x2a_a, mct_field,  mphaid, tagname, ent_type, difference)
    enddo
    call mct_list_clean(temp_list)
#endif

     call atm_import_moab(Eclock, cam_in)
#endif
   
    call t_stopf  ('CAM_import')
    
    ! Cycle over all time steps in the atm coupling interval
    
    dosend = .false.
    do while (.not. dosend)

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
       
       ! Map output from cam to mct data structures
       
       call t_startf ('CAM_export')
       call atm_export( cam_out, a2x_a%rattr )
#ifdef HAVE_MOAB
    ! call method to set all seq_flds_a2x_fields  on phys grid point cloud;
    ! it will be moved then to Atm Spectral mesh on coupler ; just to show how to move it to atm spectral
    ! on coupler
       call atm_export_moab(Eclock, cam_out)

#endif
       call t_stopf ('CAM_export')
       
    end do



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

    call seq_infodata_PutData( infodata, nextsw_cday=nextsw_cday ) 
    
    ! Write merged surface data restart file if appropriate
    
    if (rstwr_sync) then
       call t_startf('atm_write_srfrest_mct')
       call atm_write_srfrest_mct( x2a_a, a2x_a, &
            yr_spec=yr_sync, mon_spec=mon_sync, day_spec=day_sync, sec_spec=tod_sync)
       call t_stopf('atm_write_srfrest_mct')
#ifdef HAVE_MOAB
       call atm_write_srfrest_moab(yr_spec=yr_sync, &
            mon_spec=mon_sync, day_spec=day_sync, sec_spec=tod_sync)
#endif
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

  end subroutine atm_run_mct

  !================================================================================

  subroutine atm_final_mct( EClock, cdata_a, x2a_a, a2x_a)

    type(ESMF_Clock)            ,intent(inout) :: EClock
    type(seq_cdata)             ,intent(inout) :: cdata_a
    type(mct_aVect)             ,intent(inout) :: x2a_a
    type(mct_aVect)             ,intent(inout) :: a2x_a

    call t_startf('cam_final')
    call cam_final( cam_out, cam_in )
    call t_stopf('cam_final')

  end subroutine atm_final_mct

  !================================================================================

  subroutine atm_SetgsMap_mct( mpicom_atm, ATMID, GSMap_atm )

    !-------------------------------------------------------------------
    use phys_grid, only : get_nlcols_p
    !
    ! Arguments
    !
    integer        , intent(in)  :: mpicom_atm
    integer        , intent(in)  :: ATMID
    type(mct_gsMap), intent(out) :: GSMap_atm
    !
    ! Local variables
    !
    integer, allocatable :: gindex(:)
    integer :: i, n, c, ncols, sizebuf, nlcols
    integer :: ier            ! error status
    !-------------------------------------------------------------------

    ! Build the atmosphere grid numbering for MCT
    ! NOTE:  Numbering scheme is: West to East and South to North
    ! starting at south pole.  Should be the same as what's used in SCRIP
    
    ! Determine global seg map

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

    nlcols = get_nlcols_p()
    call mct_gsMap_init( gsMap_atm, gindex, mpicom_atm, ATMID, nlcols, ngcols)

    deallocate(gindex)

  end subroutine atm_SetgsMap_mct

  !===============================================================================

  subroutine atm_domain_mct( lsize, gsMap_a, dom_a )

    !-------------------------------------------------------------------
    !
    ! Arguments
    !
    integer        , intent(in)   :: lsize
    type(mct_gsMap), intent(in)   :: gsMap_a
    type(mct_ggrid), intent(inout):: dom_a  
    !
    ! Local Variables
    !
    integer  :: n,i,c,ncols           ! indices	
    real(r8) :: lats(pcols)           ! array of chunk latitudes
    real(r8) :: lons(pcols)           ! array of chunk longitude
    real(r8) :: area(pcols)           ! area in radians squared for each grid point
    real(r8), pointer  :: data(:)     ! temporary
    integer , pointer  :: idata(:)    ! temporary
    real(r8), parameter:: radtodeg = 180.0_r8/SHR_CONST_PI
    !-------------------------------------------------------------------
    !
    ! Initialize mct atm domain
    !
    call mct_gGrid_init( GGrid=dom_a, CoordChars=trim(seq_flds_dom_coord), OtherChars=trim(seq_flds_dom_other), lsize=lsize )
    !
    ! Allocate memory
    !
    allocate(data(lsize))
    !
    ! Initialize attribute vector with special value,
    ! then deallocate storage pointed to by idata
    !
    call mct_gsMap_orderedPoints(gsMap_a, iam, idata)
    call mct_gGrid_importIAttr(dom_a,'GlobGridNum',idata,lsize)
    if (associated(idata)) deallocate(idata)
    !
    ! Determine domain (numbering scheme is: West to East and South to North to South pole)
    ! Initialize attribute vector with special value
    !
    data(:) = -9999.0_R8 
    call mct_gGrid_importRAttr(dom_a,"lat"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_a,"lon"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_a,"area" ,data,lsize) 
    call mct_gGrid_importRAttr(dom_a,"aream",data,lsize) 
    data(:) = 0.0_R8     
    call mct_gGrid_importRAttr(dom_a,"mask" ,data,lsize) 
    data(:) = 1.0_R8
    call mct_gGrid_importRAttr(dom_a,"frac" ,data,lsize)
    !
    ! Fill in correct values for domain components
    !
    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, ncols, lats)
       do i=1,ncols
          n = n+1
          data(n) = lats(i)*radtodeg
       end do
    end do
    call mct_gGrid_importRAttr(dom_a,"lat",data,lsize) 

    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       call get_rlon_all_p(c, ncols, lons)
       do i=1,ncols
          n = n+1
          data(n) = lons(i)*radtodeg
       end do
    end do
    call mct_gGrid_importRAttr(dom_a,"lon",data,lsize) 

    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       call get_area_all_p(c, ncols, area)
       do i=1,ncols
          n = n+1
          data(n) = area(i) 
       end do
    end do
    call mct_gGrid_importRAttr(dom_a,"area",data,lsize) 

    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       do i=1,ncols
          n = n+1
          data(n) = 1._r8 ! mask
       end do
    end do
    call mct_gGrid_importRAttr(dom_a,"mask"   ,data,lsize) 
    deallocate(data)

  end subroutine atm_domain_mct

#ifdef HAVE_MOAB
  !===========================================================================================

  subroutine atm_read_srfrest_moab( EClock )

   !-----------------------------------------------------------------------
   use cam_pio_utils, only: cam_pio_openfile, cam_pio_closefile, pio_subsystem
   use iMOAB, only:    iMOAB_SetDoubleTagStorage, iMOAB_WriteMesh
   !
   ! Arguments
   !
   type(ESMF_Clock),intent(inout) :: EClock
   ! 
   ! Local variables
   !
   integer               :: rcode        ! return error code
   integer               :: yr_spec      ! Current year
   integer               :: mon_spec     ! Current month
   integer               :: day_spec     ! Current day
   integer               :: sec_spec     ! Current time of day (sec)
   integer               :: k
   real(r8), allocatable :: tmp(:)
   type(file_desc_t)     :: file
   type(io_desc_t)       :: iodesc
   type(var_desc_t)      :: varid
   character(CL)         :: itemc       ! string converted to char
   type(mct_string)      :: mstring     ! mct char type
   character(CXX)        :: tagname

   type(mct_list) :: temp_list
   integer :: size_list, index_list, ent_type, ierr

   !-----------------------------------------------------------------------

   ! Determine and open surface restart dataset

   call seq_timemgr_EClockGetData( EClock, curr_yr=yr_spec,curr_mon=mon_spec, &
        curr_day=day_spec, curr_tod=sec_spec ) 
   fname_srf_cam = interpret_filename_spec( rsfilename_spec_cam, case=get_restcase(), &
        yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )
   moab_fname_srf_cam = 'moab_'//trim(fname_srf_cam)
   moab_pname_srf_cam = trim(get_restartdir() )//trim(moab_fname_srf_cam)
   call getfil(moab_pname_srf_cam, moab_fname_srf_cam)
   
   call cam_pio_openfile(File, moab_fname_srf_cam, 0)
   
   call pio_initdecomp(pio_subsystem, pio_double, (/ngcols/), global_ids, iodesc)
   
   allocate(tmp(size(global_ids)))
   
   call mct_list_init(temp_list, seq_flds_x2a_fields)
   size_list=mct_list_nitem (temp_list) ! it should be the same as nrecv

   do k=1,nrecv
      call mct_list_get(mstring, k, temp_list)
      itemc = mct_string_toChar(mstring)
      call pio_seterrorhandling(File, pio_bcast_error)
      rcode = pio_inq_varid(File,'x2a_'//trim(itemc) ,varid)
      call mct_string_clean(mstring)
      if (rcode == pio_noerr) then
         call pio_read_darray(File, varid, iodesc, tmp, rcode)
         x2a_am(:,k) = tmp(:)
      else
         if (masterproc) then
            write(iulog,*)'srfrest warning: field ',trim(itemc),' is not on restart file'
            write(iulog,*)'for backwards compatibility will set it to 0'
         end if
         x2a_am(:,k) = 0._r8
      end if
      call pio_seterrorhandling(File, pio_internal_error)
   end do

   tagname=trim(seq_flds_x2a_fields)//C_NULL_CHAR
    ent_type = 0 ! vertices, point cloud
    ierr = iMOAB_SetDoubleTagStorage ( mphaid, tagname, totalmbls_r , ent_type, x2a_am )
    if ( ierr > 0) then
      call endrun('Error: fail to set  seq_flds_a2x_fields for atm physgrid moab mesh in restart')
    endif

   
   call mct_list_clean(temp_list)

   call mct_list_init(temp_list, seq_flds_a2x_fields)

   do k=1,nsend
      call mct_list_get(mstring, k, temp_list)
      itemc = mct_string_toChar(mstring)
      rcode = pio_inq_varid(File,'a2x_'//trim(itemc) ,varid)
      call mct_string_clean(mstring)

      if (rcode == pio_noerr) then
         call pio_read_darray(File, varid, iodesc, tmp, rcode)
         a2x_am(:,k) = tmp(:)
      else
         if (masterproc) then
            write(iulog,*)'srfrest warning: field ',trim(itemc),' is not on restart file'
            write(iulog,*)'for backwards compatibility will set it to 0'
         end if
         a2x_am(:,k) = 0._r8
      endif
   end do
   call mct_list_clean(temp_list)
   tagname=trim(seq_flds_a2x_fields)//C_NULL_CHAR
    ent_type = 0 ! vertices, point cloud
    ierr = iMOAB_SetDoubleTagStorage ( mphaid, tagname, totalmbls , ent_type, a2x_am )
    if ( ierr > 0) then
      call endrun('Error: fail to set  seq_flds_a2x_fields for atm physgrid moab mesh in restart')
    endif
   
   ! write moab phys atm after reading restart surface file

   call pio_freedecomp(File,iodesc)
   call cam_pio_closefile(File)
   deallocate(tmp)

 end subroutine atm_read_srfrest_moab

 !===========================================================================================

 subroutine atm_write_srfrest_moab( yr_spec, mon_spec, day_spec, sec_spec )

   !-----------------------------------------------------------------------
   use cam_pio_utils, only: cam_pio_createfile, cam_pio_closefile, pio_subsystem
   use cam_pio_utils, only: cam_pio_openfile
   use cam_history_support, only: fillvalue
   use iMOAB, only:    iMOAB_GetDoubleTagStorage, iMOAB_WriteMesh
   !
   ! Arguments
   !
   integer        , intent(in) :: yr_spec         ! Simulation year
   integer        , intent(in) :: mon_spec        ! Simulation month
   integer        , intent(in) :: day_spec        ! Simulation day
   integer        , intent(in) :: sec_spec        ! Seconds into current simulation day
   !
   ! Local variables
   !
   integer                   :: rcode        ! return error code
   integer                   :: dimid(1), k
   type(file_desc_t)         :: file
   real(r8), allocatable     :: tmp(:)
   type(var_desc_t), pointer :: varid_x2a(:), varid_a2x(:)
   type(io_desc_t)           :: iodesc
   character(CL)             :: itemc       ! string converted to char

   type(mct_string)          :: mstring     ! mct char type
   character(CXX)            :: tagname

   type(mct_list) :: temp_list
   integer :: size_list, index_list, ent_type, ierr

   !-----------------------------------------------------------------------

   ! Determine and open surface restart dataset

      ! Determine and open surface restart dataset

   moab_fname_srf_cam = 'moab_'//trim(fname_srf_cam)
   
   call cam_pio_createfile(File, trim(moab_fname_srf_cam))
   if (masterproc) then
      write(iulog,*)'create file :', trim(moab_fname_srf_cam) 
   end if

   call pio_initdecomp(pio_subsystem, pio_double, (/ngcols/), global_ids, iodesc)
   allocate(tmp(size(global_ids)))
   
   rcode = pio_def_dim(File,'x2a_nx',ngcols,dimid(1))
   call mct_list_init(temp_list ,seq_flds_x2a_fields)
   size_list=mct_list_nitem (temp_list) ! it should be the same as nrecv
   allocate(varid_x2a(size_list))
   if (masterproc) then
      write(iulog,*)'size list:', size_list, seq_flds_x2a_fields 
   end if
    
   do k = 1, size_list
      call mct_list_get(mstring, k, temp_list)
      itemc = mct_string_toChar(mstring)
      rcode = pio_def_var(File,'x2a_'//trim(itemc),PIO_DOUBLE,dimid,varid_x2a(k))
      call mct_string_clean(mstring)
      rcode = pio_put_att(File,varid_x2a(k),"_fillvalue",fillvalue)
   enddo
   call mct_list_clean(temp_list)

   call mct_list_init(temp_list ,seq_flds_a2x_fields)
   size_list=mct_list_nitem (temp_list) ! it should be the same as nsend
   allocate(varid_a2x(size_list))
   
   rcode = pio_def_dim(File,'a2x_nx',ngcols,dimid(1))
   do k = 1,size_list
      call mct_list_get(mstring,k,temp_list)
      itemc = mct_string_toChar(mstring)
      rcode = PIO_def_var(File,'a2x_'//trim(itemc),PIO_DOUBLE,dimid,varid_a2x(k))
      call mct_string_clean(mstring)
      rcode = PIO_put_att(File,varid_a2x(k),"_fillvalue",fillvalue)
   enddo

   rcode = pio_enddef(File)  ! don't check return code, might be enddef already

! do we need to fill it up with values? 
   ! ccsm sign convention is that fluxes are positive downward
   tagname=trim(seq_flds_x2a_fields)//C_NULL_CHAR
   ent_type = 0 ! vertices, point cloud
   ierr = iMOAB_GetDoubleTagStorage ( mphaid, tagname, totalmbls_r , ent_type, x2a_am )
   if ( ierr > 0) then
     call endrun('Error: fail to get  seq_flds_x2a_fields for atm physgrid moab mesh for writing restart surface')
   endif

   tagname=trim(seq_flds_a2x_fields)//C_NULL_CHAR
   ent_type = 0 ! vertices, point cloud
   ierr = iMOAB_GetDoubleTagStorage ( mphaid, tagname, totalmbls , ent_type, a2x_am )
   if ( ierr > 0) then
     call endrun('Error: fail to get  seq_flds_a2x_fields for atm physgrid moab mesh for writing restart surface')
   endif

   do k=1,nrecv
      call pio_write_darray(File, varid_x2a(k), iodesc, x2a_am(:,k), rcode)
   end do

   do k=1,nsend
      call pio_write_darray(File, varid_a2x(k), iodesc, a2x_am(:,k), rcode)       
   end do

   deallocate(varid_x2a, varid_a2x)

   call pio_freedecomp(File,iodesc)
   call cam_pio_closefile(file)


 end subroutine atm_write_srfrest_moab
#endif

  !===========================================================================================

  subroutine atm_read_srfrest_mct( EClock, x2a_a, a2x_a)

    !-----------------------------------------------------------------------
    use cam_pio_utils, only: cam_pio_openfile, cam_pio_closefile, pio_subsystem
    !
    ! Arguments
    !
    type(ESMF_Clock),intent(inout) :: EClock
    type(mct_aVect), intent(inout) :: x2a_a
    type(mct_aVect), intent(inout) :: a2x_a
    ! 
    ! Local variables
    !
    integer               :: rcode        ! return error code
    integer               :: yr_spec      ! Current year
    integer               :: mon_spec     ! Current month
    integer               :: day_spec     ! Current day
    integer               :: sec_spec     ! Current time of day (sec)
    integer               :: nf_x2a, nf_a2x, k
    real(r8), allocatable :: tmp(:)
    type(file_desc_t)     :: file
    type(io_desc_t)       :: iodesc
    type(var_desc_t)      :: varid
    character(CL)         :: itemc       ! string converted to char
    type(mct_string)      :: mstring     ! mct char type
    !-----------------------------------------------------------------------

    ! Determine and open surface restart dataset

    call seq_timemgr_EClockGetData( EClock, curr_yr=yr_spec,curr_mon=mon_spec, &
         curr_day=day_spec, curr_tod=sec_spec ) 
    fname_srf_cam = interpret_filename_spec( rsfilename_spec_cam, case=get_restcase(), &
         yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )
    pname_srf_cam = trim(get_restartdir() )//fname_srf_cam
    call getfil(pname_srf_cam, fname_srf_cam)
    
    call cam_pio_openfile(File, fname_srf_cam, 0)
    call pio_initdecomp(pio_subsystem, pio_double, (/ngcols/), dof, iodesc)
    allocate(tmp(size(dof)))
    
    nf_x2a = mct_aVect_nRattr(x2a_a)
    do k=1,nf_x2a
       call mct_aVect_getRList(mstring,k,x2a_a)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)

       call pio_seterrorhandling(File, pio_bcast_error)
       rcode = pio_inq_varid(File,'x2a_'//trim(itemc) ,varid)
       if (rcode == pio_noerr) then
          call pio_read_darray(File, varid, iodesc, tmp, rcode)
          x2a_a%rattr(k,:) = tmp(:)
       else
       if (masterproc) then
            write(iulog,*)'srfrest warning: field ',trim(itemc),' is not on restart file'
            write(iulog,*)'for backwards compatibility will set it to 0'
          end if
          x2a_a%rattr(k,:) = 0._r8
       end if
       call pio_seterrorhandling(File, pio_internal_error)
    end do

    nf_a2x = mct_aVect_nRattr(a2x_a)
    do k=1,nf_a2x
       call mct_aVect_getRList(mstring,k,a2x_a)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)

       rcode = pio_inq_varid(File,'a2x_'//trim(itemc) ,varid)
       call pio_read_darray(File, varid, iodesc, tmp, rcode)
       a2x_a%rattr(k,:) = tmp(:)
    end do

    call pio_freedecomp(File,iodesc)
    call cam_pio_closefile(File)
    deallocate(tmp)

  end subroutine atm_read_srfrest_mct

  !===========================================================================================

  subroutine atm_write_srfrest_mct( x2a_a, a2x_a, & 
       yr_spec, mon_spec, day_spec, sec_spec)

    !-----------------------------------------------------------------------
    use cam_pio_utils, only: cam_pio_createfile, cam_pio_closefile, pio_subsystem
    use cam_history_support, only: fillvalue
    !
    ! Arguments
    !
    type(mct_aVect), intent(in) :: x2a_a
    type(mct_aVect), intent(in) :: a2x_a
    integer        , intent(in) :: yr_spec         ! Simulation year
    integer        , intent(in) :: mon_spec        ! Simulation month
    integer        , intent(in) :: day_spec        ! Simulation day
    integer        , intent(in) :: sec_spec        ! Seconds into current simulation day
    !
    ! Local variables
    !
    integer                   :: rcode        ! return error code
    integer                   :: nf_x2a, nf_a2x, dimid(1), k
    type(file_desc_t)         :: file
    type(var_desc_t), pointer :: varid_x2a(:), varid_a2x(:)
    type(io_desc_t)           :: iodesc
    character(CL)             :: itemc       ! string converted to char
    type(mct_string)          :: mstring     ! mct char type
    !-----------------------------------------------------------------------

    ! Determine and open surface restart dataset

    fname_srf_cam = interpret_filename_spec( rsfilename_spec_cam, &
         yr_spec=yr_spec, mon_spec=mon_spec, day_spec=day_spec, sec_spec= sec_spec )

    call cam_pio_createfile(File, fname_srf_cam)
    call pio_initdecomp(pio_subsystem, pio_double, (/ngcols/), dof, iodesc)

    nf_x2a = mct_aVect_nRattr(x2a_a)
    allocate(varid_x2a(nf_x2a))
    
    rcode = pio_def_dim(File,'x2a_nx',ngcols,dimid(1))
    do k = 1,nf_x2a
       call mct_aVect_getRList(mstring,k,x2a_a)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)
       rcode = pio_def_var(File,'x2a_'//trim(itemc),PIO_DOUBLE,dimid,varid_x2a(k))
       rcode = pio_put_att(File,varid_x2a(k),"_fillvalue",fillvalue)
    enddo

    nf_a2x = mct_aVect_nRattr(a2x_a)
    allocate(varid_a2x(nf_a2x))
    
    rcode = pio_def_dim(File,'a2x_nx',ngcols,dimid(1))
    do k = 1,nf_a2x
       call mct_aVect_getRList(mstring,k,a2x_a)
       itemc = mct_string_toChar(mstring)
       call mct_string_clean(mstring)
       rcode = PIO_def_var(File,'a2x_'//trim(itemc),PIO_DOUBLE,dimid,varid_a2x(k))
       rcode = PIO_put_att(File,varid_a2x(k),"_fillvalue",fillvalue)
    enddo

    rcode = pio_enddef(File)  ! don't check return code, might be enddef already


    do k=1,nf_x2a
       call pio_write_darray(File, varid_x2a(k), iodesc, x2a_a%rattr(k,:), rcode)
    end do

    do k=1,nf_a2x
       call pio_write_darray(File, varid_a2x(k), iodesc, a2x_a%rattr(k,:), rcode)       
    end do

    deallocate(varid_x2a, varid_a2x)

    call pio_freedecomp(File,iodesc)
    call cam_pio_closefile(file)


  end subroutine atm_write_srfrest_mct

#ifdef HAVE_MOAB
  subroutine init_moab_atm_phys( cdata_a )

    use shr_mpi_mod,       only: shr_mpi_commrank, shr_mpi_commsize
    use shr_const_mod, only: SHR_CONST_PI
!-------------------------------------------------------------------
    use phys_grid, only : get_nlcols_p ! used to det local size ?
    use iMOAB, only : iMOAB_RegisterApplication, iMOAB_CreateVertices, iMOAB_WriteMesh, &
      iMOAB_DefineTagStorage, iMOAB_SetIntTagStorage, iMOAB_SetDoubleTagStorage, &
      iMOAB_ResolveSharedEntities, iMOAB_UpdateMeshInfo

    type(seq_cdata), intent(in)              :: cdata_a


    integer :: ATMID
    integer :: mpicom_atm

    integer :: ATM_PHYS ! our numbering

    ! local variables to fill in data
    integer, dimension(:), allocatable :: vgids
    !  retrieve everything we need from mct
    ! number of vertices is the size of mct grid
    real(r8), dimension(:), allocatable :: moab_vert_coords  ! temporary
    real(r8), dimension(:), allocatable :: areavals
    ! r
    real(r8)   :: latv, lonv
    integer   dims, i, ilat, ilon, igdx, ierr, tagindex
    integer tagtype, numco, ent_type

    integer ::  n, c, ncols, nlcols

    real(r8) :: lats(pcols)           ! array of chunk latitudes pcol is defined in ppgrid
    real(r8) :: lons(pcols)           ! array of chunk longitude
    real(r8) :: area(pcols)           ! area in radians squared for each grid point
    integer , dimension(:), allocatable :: chunk_index(:)    ! temporary
    real(r8), parameter:: radtodeg = 180.0_r8/SHR_CONST_PI

    character*100 outfile, wopts
    character(CXX) :: tagname ! will store all seq_flds_a2x_fields 
    character*32 appname


    call seq_cdata_setptrs(cdata_a, ID=ATMID, mpicom=mpicom_atm, &
         infodata=infodata)

    appname="ATM_PHYS"//C_NULL_CHAR
    ATM_PHYS = 200 + ATMID !
    ierr = iMOAB_RegisterApplication(appname, mpicom_atm, ATM_PHYS, mphaid)
    if (ierr > 0 )  &
       call endrun('Error: cannot register moab app for atm physics')
    if(masterproc) then
       write(iulog,*) " "
       write(iulog,*) "register MOAB app:", trim(appname), "  mphaid=", mphaid
       write(iulog,*) " "
    endif

    ! first, determine the size of local vertices

    nlcols = get_nlcols_p()
    dims = 3 !
    allocate(vgids(nlcols))
    allocate(global_ids(nlcols))
    allocate(moab_vert_coords(nlcols*dims))
    allocate(areavals(nlcols))
    allocate(chunk_index(nlcols))
    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, ncols, lats) !
       call get_rlon_all_p(c, ncols, lons)
       call get_area_all_p(c, ncols, area)
       do i = 1,ncols
          n=n+1
          vgids(n) = get_gcol_p(c,i)
          global_ids(n) = vgids(n)
          latv = lats(i) ! these are in rads ?
          lonv = lons(i)
          moab_vert_coords(3*n-2)=COS(latv)*COS(lonv)
          moab_vert_coords(3*n-1)=COS(latv)*SIN(lonv)
          moab_vert_coords(3*n  )=SIN(latv)
          areavals(n) = area(i)
          chunk_index(n) = c ! this is just for us, to see the chunk
       end do
    end do


    ierr = iMOAB_CreateVertices(mphaid, nlcols*3, dims, moab_vert_coords)
    if (ierr > 0 )  &
      call endrun('Error: fail to create MOAB vertices in phys atm model')

    tagtype = 0  ! dense, integer
    numco = 1
    tagname='GLOBAL_ID'//C_NULL_CHAR
    ierr = iMOAB_DefineTagStorage(mphaid, tagname, tagtype, numco,  tagindex )
    if (ierr > 0 )  &
      call endrun('Error: fail to retrieve GLOBAL_ID tag ')

    ent_type = 0 ! vertex type
    ierr = iMOAB_SetIntTagStorage ( mphaid, tagname, nlcols , ent_type, vgids)
    if (ierr > 0 )  &
      call endrun('Error: fail to set GLOBAL_ID tag ')

    ierr = iMOAB_ResolveSharedEntities( mphaid, nlcols, vgids );
    if (ierr > 0 )  &
      call endrun('Error: fail to resolve shared entities')

    ierr = iMOAB_UpdateMeshInfo(mphaid)
    if (ierr > 0 )  &
      call endrun('Error: fail to update mesh info')

    !there are no shared entities, but we will set a special partition tag, in order to see the
    ! partitions ; it will be visible with a Pseudocolor plot in VisIt
    tagname='partition'//C_NULL_CHAR
    ierr = iMOAB_DefineTagStorage(mphaid, tagname, tagtype, numco,  tagindex )
    if (ierr > 0 )  &
      call endrun('Error: fail to create new partition tag ')

    vgids = iam
    ierr = iMOAB_SetIntTagStorage ( mphaid, tagname, nlcols , ent_type, vgids)
    if (ierr > 0 )  &
      call endrun('Error: fail to set partition tag ')

    ! chunk_index ; it will be visible with a Pseudocolor plot in VisIt
    tagname='chunk_id'//C_NULL_CHAR
    ierr = iMOAB_DefineTagStorage(mphaid, tagname, tagtype, numco,  tagindex )
    if (ierr > 0 )  &
      call endrun('Error: fail to create new chunk index tag ')

    ierr = iMOAB_SetIntTagStorage ( mphaid, tagname, nlcols , ent_type, chunk_index)
    if (ierr > 0 )  &
      call endrun('Error: fail to set chunk id tag tag ')

    ! use areavals for area, aream; define also dom fields 
    ! define all  on moab 

    tagname=trim(seq_flds_dom_fields)//C_NULL_CHAR !  mask is double too lat:lon:hgt:area:aream:mask:frac
    tagtype = 1 ! dense, double
    ierr = iMOAB_DefineTagStorage(mphaid, tagname, tagtype, numco,  tagindex )
    if (ierr > 0 )  &
      call endrun('Error: fail to create  tags from seq_flds_dom_fields ')

    tagname='area'//C_NULL_CHAR
    ierr = iMOAB_SetDoubleTagStorage ( mphaid, tagname, nlcols , ent_type, areavals)
    if (ierr > 0 )  &
      call endrun('Error: fail to set area tag ')

    areavals = 1._r8 ! double
    tagname='mask'//C_NULL_CHAR
    ierr = iMOAB_SetDoubleTagStorage ( mphaid, tagname, nlcols , ent_type, areavals)
    if (ierr > 0 )  &
      call endrun('Error: fail to set mask tag ')

    areavals = 1._r8 ! double
    
    ! set lat, lon, and frac tags at the same time, reusing the moab_vert_coords array already allocated
    ! we set all at the same time, so use the 1-d array carefully, with values interlaced/or not 
    n=0
    do c = begchunk, endchunk
       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, ncols, lats) !
       call get_rlon_all_p(c, ncols, lons)
       do i = 1,ncols
          n=n+1
          vgids(n) = get_gcol_p(c,i)
          latv = lats(i) ! these are in rads ?
          lonv = lons(i)
          moab_vert_coords(            n ) = lats(i) * radtodeg
          moab_vert_coords(  nlcols +  n ) = lons(i) * radtodeg
          moab_vert_coords( 2*nlcols + n )= 1._r8 ! this for fractions
       end do
    end do
    tagname = 'lat:lon:frac'//C_NULL_CHAR
    ierr = iMOAB_SetDoubleTagStorage ( mphaid, tagname, nlcols*3 , ent_type, moab_vert_coords)
    if (ierr > 0 )  &
      call endrun('Error: fail to set lat lon frac tag ')


#ifdef MOABDEBUG
    outfile = 'AtmPhys.h5m'//C_NULL_CHAR
    wopts   = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
    ierr = iMOAB_WriteMesh(mphaid, outfile, wopts)
    if (ierr > 0 )  &
      call endrun('Error: fail to write the atm phys mesh file')
#endif
   ! define fields seq_flds_a2x_fields 
    tagtype = 1  ! dense, double
    numco = 1 !  one value per vertex / entity
    tagname = trim(seq_flds_a2x_fields)//C_NULL_CHAR
    ierr = iMOAB_DefineTagStorage(mphaid, tagname, tagtype, numco,  tagindex )
    if ( ierr > 0) then
       call endrun('Error: fail to define seq_flds_a2x_fields for atm physgrid moab mesh')
    endif
    ! make sure this is defined too; it could have the same fields, but in different order, or really different 
    ! fields; need to make sure we have them
    tagname = trim(seq_flds_x2a_fields)//C_NULL_CHAR
    ierr = iMOAB_DefineTagStorage(mphaid, tagname, tagtype, numco,  tagindex )
    if ( ierr > 0) then
       call endrun('Error: fail to define seq_flds_x2a_fields for atm physgrid moab mesh')
    endif

    deallocate(moab_vert_coords)
    deallocate(vgids)
    deallocate(areavals)
    deallocate(chunk_index)

  end subroutine init_moab_atm_phys

  subroutine atm_export_moab(Eclock, cam_out)
  !-------------------------------------------------------------------
    use camsrfexch, only: cam_out_t
    use phys_grid , only: get_ncols_p, get_nlcols_p
    use ppgrid    , only: begchunk, endchunk
    use seq_comm_mct, only: mphaid ! imoab pid for atm physics
    use cam_abortutils       , only: endrun
    use iMOAB, only:  iMOAB_WriteMesh,  iMOAB_SetDoubleTagStorage
    use iso_c_binding 
    !
    ! Arguments
    !
    type(ESMF_Clock),intent(inout) :: EClock
    type(cam_out_t), intent(in)    :: cam_out(begchunk:endchunk)

    integer tagtype, numco, ent_type
    character(CXX) ::  tagname ! 

    integer ierr, c, nlcols, ig, i, ncols
    integer  :: cur_atm_stepno

#ifdef MOABDEBUG
    character*100 outfile, wopts, lnum
    integer, save :: local_count = 0
    character*100 lnum2
#endif 
    ! Copy from component arrays into chunk array data structure
    ! Rearrange data from chunk structure into lat-lon buffer and subsequently
    ! create double array for moab tags

    ig=1
    do c=begchunk, endchunk
       ncols = get_ncols_p(c)
       do i=1,ncols
          a2x_am(ig, index_a2x_Sa_pslv   ) = cam_out(c)%psl(i)
          a2x_am(ig, index_a2x_Sa_z      ) = cam_out(c)%zbot(i)   
          a2x_am(ig, index_a2x_Sa_u      ) = cam_out(c)%ubot(i)   
          a2x_am(ig, index_a2x_Sa_v      ) = cam_out(c)%vbot(i)   
          a2x_am(ig, index_a2x_Sa_tbot   ) = cam_out(c)%tbot(i)   
          a2x_am(ig, index_a2x_Sa_ptem   ) = cam_out(c)%thbot(i)  
          a2x_am(ig, index_a2x_Sa_pbot   ) = cam_out(c)%pbot(i)   
          a2x_am(ig, index_a2x_Sa_shum   ) = cam_out(c)%qbot(i,1) 
          a2x_am(ig, index_a2x_Sa_dens   ) = cam_out(c)%rho(i)
          if (trim(adjustl(precip_downscaling_method)) == "FNM") then
             !if the land model's precip downscaling method is FNM, export uovern to the coupler
             a2x_am(ig, index_a2x_Sa_uovern) = cam_out(c)%uovern(i)
          else
             a2x_am(ig, index_a2x_Sa_uovern) = 0._R8
          endif

          a2x_am(ig, index_a2x_Faxa_swnet) = cam_out(c)%netsw(i)      
          a2x_am(ig, index_a2x_Faxa_lwdn ) = cam_out(c)%flwds(i)  
          a2x_am(ig, index_a2x_Faxa_rainc) = (cam_out(c)%precc(i)-cam_out(c)%precsc(i))*1000._r8
          a2x_am(ig, index_a2x_Faxa_rainl) = (cam_out(c)%precl(i)-cam_out(c)%precsl(i))*1000._r8
          a2x_am(ig, index_a2x_Faxa_snowc) = cam_out(c)%precsc(i)*1000._r8
          a2x_am(ig, index_a2x_Faxa_snowl) = cam_out(c)%precsl(i)*1000._r8
          a2x_am(ig, index_a2x_Faxa_swndr) = cam_out(c)%soll(i)   
          a2x_am(ig, index_a2x_Faxa_swvdr) = cam_out(c)%sols(i)   
          a2x_am(ig, index_a2x_Faxa_swndf) = cam_out(c)%solld(i)  
          a2x_am(ig, index_a2x_Faxa_swvdf) = cam_out(c)%solsd(i)  

          ! aerosol deposition fluxes
          a2x_am(ig, index_a2x_Faxa_bcphidry) = cam_out(c)%bcphidry(i)
          a2x_am(ig, index_a2x_Faxa_bcphodry) = cam_out(c)%bcphodry(i)
          a2x_am(ig, index_a2x_Faxa_bcphiwet) = cam_out(c)%bcphiwet(i)
          a2x_am(ig, index_a2x_Faxa_ocphidry) = cam_out(c)%ocphidry(i)
          a2x_am(ig, index_a2x_Faxa_ocphodry) = cam_out(c)%ocphodry(i)
          a2x_am(ig, index_a2x_Faxa_ocphiwet) = cam_out(c)%ocphiwet(i)
          a2x_am(ig, index_a2x_Faxa_dstwet1)  = cam_out(c)%dstwet1(i)
          a2x_am(ig, index_a2x_Faxa_dstdry1)  = cam_out(c)%dstdry1(i)
          a2x_am(ig, index_a2x_Faxa_dstwet2)  = cam_out(c)%dstwet2(i)
          a2x_am(ig, index_a2x_Faxa_dstdry2)  = cam_out(c)%dstdry2(i)
          a2x_am(ig, index_a2x_Faxa_dstwet3)  = cam_out(c)%dstwet3(i)
          a2x_am(ig, index_a2x_Faxa_dstdry3)  = cam_out(c)%dstdry3(i)
          a2x_am(ig, index_a2x_Faxa_dstwet4)  = cam_out(c)%dstwet4(i)
          a2x_am(ig, index_a2x_Faxa_dstdry4)  = cam_out(c)%dstdry4(i)

          if (index_a2x_Sa_co2prog /= 0) then
             a2x_am(ig, index_a2x_Sa_co2prog) = cam_out(c)%co2prog(i) ! atm prognostic co2
          end if
          if (index_a2x_Sa_co2diag /= 0) then
             a2x_am(ig, index_a2x_Sa_co2diag) = cam_out(c)%co2diag(i) ! atm diagnostic co2
          end if

          ig=ig+1
       end do
    end do
    tagname=trim(seq_flds_a2x_fields)//C_NULL_CHAR
    ent_type = 0 ! vertices, point cloud
    ierr = iMOAB_SetDoubleTagStorage ( mphaid, tagname, totalmbls , ent_type, a2x_am )
    if ( ierr > 0) then
      call endrun('Error: fail to set  seq_flds_a2x_fields for atm physgrid moab mesh')
    endif
    call seq_timemgr_EClockGetData( EClock, stepno=cur_atm_stepno )
#ifdef MOABDEBUG
    write(lnum,"(I0.2)")cur_atm_stepno
    local_count = local_count + 1
    write(lnum2,"(I0.2)")local_count
    outfile = 'atm_export_'//trim(lnum)//'_'//trim(lnum2)//'.h5m'//C_NULL_CHAR
    wopts   = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
    ierr = iMOAB_WriteMesh(mphaid, outfile, wopts)
    if (ierr > 0 )  &
      call endrun('Error: fail to write the atm phys mesh file with data')
#endif
    

  end subroutine atm_export_moab

subroutine atm_import_moab(Eclock, cam_in, restart_init )

    !-----------------------------------------------------------------------
    use cam_cpl_indices
    use camsrfexch,     only: cam_in_t
    use phys_grid ,     only: get_ncols_p
    use ppgrid    ,     only: begchunk, endchunk       
    use shr_const_mod,  only: shr_const_stebol
    use seq_drydep_mod, only: n_drydep
    use co2_cycle     , only: c_i, co2_readFlux_ocn, co2_readFlux_fuel
    use co2_cycle     , only: co2_transport, co2_time_interp_ocn, co2_time_interp_fuel
    use co2_cycle     , only: data_flux_ocn, data_flux_fuel
    use physconst     , only: mwco2
    use time_manager  , only: is_first_step
    use iMOAB, only:  iMOAB_WriteMesh,  iMOAB_GetDoubleTagStorage
    use iso_c_binding
    !
    ! Arguments
    !
    ! real(r8)      , intent(in)    :: x2a_am(:,:) will be retrieved from moab tags, and used to set cam_in 
    type(ESMF_Clock),intent(inout) :: EClock
    type(cam_in_t), intent(inout) :: cam_in(begchunk:endchunk)
    logical, optional, intent(in) :: restart_init
    !
    ! Local variables
    !		
    integer            :: i,lat,n,c,ig  ! indices
    integer            :: ncols         ! number of columns
    logical, save      :: first_time = .true.
    integer, parameter :: ndst = 2
    integer, target    :: spc_ndx(ndst)
    integer, pointer   :: dst_a5_ndx, dst_a7_ndx
    integer, pointer   :: dst_a1_ndx, dst_a3_ndx
    logical :: overwrite_flds

    character(CXX) ::  tagname ! 
    integer  :: ent_type, ierr
    integer  :: cur_atm_stepno

    call seq_timemgr_EClockGetData( EClock, stepno=cur_atm_stepno )
    !-----------------------------------------------------------------------
    overwrite_flds = .true.
    ! don't overwrite fields if invoked during the initialization phase 
    ! of a 'continue' or 'branch' run type with data from .rs file
    if (present(restart_init)) overwrite_flds = .not. restart_init
    

    ! ccsm sign convention is that fluxes are positive downward
    tagname=trim(seq_flds_x2a_fields)//C_NULL_CHAR
    ent_type = 0 ! vertices, point cloud
    ierr = iMOAB_GetDoubleTagStorage ( mphaid, tagname, totalmbls_r , ent_type, x2a_am )
    if ( ierr > 0) then
      call endrun('Error: fail to get  seq_flds_a2x_fields for atm physgrid moab mesh')
    endif

    ig=1
    do c=begchunk,endchunk
       ncols = get_ncols_p(c) 

       ! initialize constituent surface fluxes to zero
       ! NOTE:overwrite_flds is .FALSE. for the first restart
       ! time step making cflx(:,1)=0.0 for the first restart time step.
       ! cflx(:,1) should not be zeroed out, start the second index of cflx from 2.
       ! cam_in(c)%cflx(:,2:) = 0._r8 
                                               
       do i =1,ncols                                                               
          if (overwrite_flds) then
             ! Prior to this change, "overwrite_flds" was always .true. therefore wsx and wsy were always updated.
             ! Now, overwrite_flds is .false. for the first time step of the restart run. Move wsx and wsy out of 
             ! this if-condition so that they are still updated everytime irrespective of the value of overwrite_flds.

             ! Move lhf to this if-block so that it is not overwritten to ensure BFB restarts when qneg4 correction 
             ! occurs at the restart time step
             ! Modified by Wuyin Lin
             cam_in(c)%shf(i)    = -x2a_am(ig,index_x2a_Faxx_sen)     
             cam_in(c)%cflx(i,1) = -x2a_am(ig,index_x2a_Faxx_evap)                
             cam_in(c)%lhf(i)    = -x2a_am(ig,index_x2a_Faxx_lat)     
          endif

          if (index_x2a_Faoo_h2otemp /= 0) then
             cam_in(c)%h2otemp(i) = -x2a_am(ig,index_x2a_Faoo_h2otemp)
          end if
           
          cam_in(c)%wsx(i)    = -x2a_am(ig,index_x2a_Faxx_taux)     
          cam_in(c)%wsy(i)    = -x2a_am(ig,index_x2a_Faxx_tauy)     
          cam_in(c)%lwup(i)      = -x2a_am(ig,index_x2a_Faxx_lwup)    
          cam_in(c)%asdir(i)     =  x2a_am(ig,index_x2a_Sx_avsdr)  
          cam_in(c)%aldir(i)     =  x2a_am(ig,index_x2a_Sx_anidr)  
          cam_in(c)%asdif(i)     =  x2a_am(ig,index_x2a_Sx_avsdf)  
          cam_in(c)%aldif(i)     =  x2a_am(ig,index_x2a_Sx_anidf)
          cam_in(c)%ts(i)        =  x2a_am(ig,index_x2a_Sx_t)  
          cam_in(c)%sst(i)       =  x2a_am(ig,index_x2a_So_t)             
          cam_in(c)%snowhland(i) =  x2a_am(ig,index_x2a_Sl_snowh)  
          cam_in(c)%snowhice(i)  =  x2a_am(ig,index_x2a_Si_snowh)  
          cam_in(c)%tref(i)      =  x2a_am(ig,index_x2a_Sx_tref)  
          cam_in(c)%qref(i)      =  x2a_am(ig,index_x2a_Sx_qref)
          cam_in(c)%u10(i)       =  x2a_am(ig,index_x2a_Sx_u10)
          cam_in(c)%icefrac(i)   =  x2a_am(ig,index_x2a_Sf_ifrac)  
          cam_in(c)%ocnfrac(i)   =  x2a_am(ig,index_x2a_Sf_ofrac)
          cam_in(c)%landfrac(i)  =  x2a_am(ig,index_x2a_Sf_lfrac)
          if ( associated(cam_in(c)%ram1) ) &
               cam_in(c)%ram1(i) =  x2a_am(ig,index_x2a_Sl_ram1 )
          if ( associated(cam_in(c)%fv) ) &
               cam_in(c)%fv(i)   =  x2a_am(ig,index_x2a_Sl_fv   )
          if ( associated(cam_in(c)%soilw) ) &
               cam_in(c)%soilw(i) =  x2a_am(ig,index_x2a_Sl_soilw)
          if ( associated(cam_in(c)%dstflx) ) then
             cam_in(c)%dstflx(i,1) = x2a_am(ig,index_x2a_Fall_flxdst1)
             cam_in(c)%dstflx(i,2) = x2a_am(ig,index_x2a_Fall_flxdst2)
             cam_in(c)%dstflx(i,3) = x2a_am(ig,index_x2a_Fall_flxdst3)
             cam_in(c)%dstflx(i,4) = x2a_am(ig,index_x2a_Fall_flxdst4)
          endif
          if ( associated(cam_in(c)%meganflx) ) then
             cam_in(c)%meganflx(i,1:shr_megan_mechcomps_n) = &
                  x2a_am(ig,index_x2a_Fall_flxvoc:index_x2a_Fall_flxvoc+shr_megan_mechcomps_n-1)
          endif

          ! dry dep velocities
          if ( index_x2a_Sl_ddvel/=0 .and. n_drydep>0 ) then
             cam_in(c)%depvel(i,:n_drydep) = &
                  x2a_am(ig,index_x2a_Sl_ddvel:index_x2a_Sl_ddvel+n_drydep-1)
          endif
          !
          ! fields needed to calculate water isotopes to ocean evaporation processes
          !
          cam_in(c)%ustar(i) = x2a_am(ig,index_x2a_So_ustar)
          cam_in(c)%re(i)    = x2a_am(ig,index_x2a_So_re   )
          cam_in(c)%ssq(i)   = x2a_am(ig,index_x2a_So_ssq  )
          !
          ! bgc scenarios
          !
          if (index_x2a_Fall_fco2_lnd /= 0) then
             cam_in(c)%fco2_lnd(i) = -x2a_am(ig,index_x2a_Fall_fco2_lnd)
          end if
          if (index_x2a_Faoo_fco2_ocn /= 0) then
             cam_in(c)%fco2_ocn(i) = -x2a_am(ig,index_x2a_Faoo_fco2_ocn)
          end if
          if (index_x2a_Faoo_fdms_ocn /= 0) then
             cam_in(c)%fdms(i)     = -x2a_am(ig,index_x2a_Faoo_fdms_ocn)
          end if

          ig=ig+1

       end do
    end do

    ! Get total co2 flux from components,
    ! Note - co2_transport determines if cam_in(c)%cflx(i,c_i(1:4)) is allocated

    if (co2_transport().and.overwrite_flds) then

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
! The below section involves a temporary workaround for fluxes from data (read in from a file)
! There is an issue with infld that does not allow time-varying 2D files to be read correctly.
! The work around involves adding a singleton 3rd dimension offline and reading the files as 
! 3D fields.  Once this issue is corrected, the old implementation can be reinstated.
! This is the case for both data_flux_ocn and data_flux_fuel
!++BEH  vvv old implementation vvv
!                cam_in(c)%cflx(i,c_i(1)) = &
!                     -data_flux_ocn%co2flx(i,c)*(1._r8- cam_in(c)%landfrac(i)) &
!                     *mwco2*1.0e-3_r8
!       ^^^ old implementation ^^^   ///    vvv new implementation vvv
                cam_in(c)%cflx(i,c_i(1)) = &
                     -data_flux_ocn%co2flx(i,1,c)*(1._r8- cam_in(c)%landfrac(i)) &
                     *mwco2*1.0e-3_r8
!--BEH  ^^^ new implementation ^^^
             else
                cam_in(c)%cflx(i,c_i(1)) = 0._r8
             end if
             
             ! co2 flux from fossil fuel
             if (co2_readFlux_fuel) then
!++BEH  vvv old implementation vvv
!                cam_in(c)%cflx(i,c_i(2)) = data_flux_fuel%co2flx(i,c)
!       ^^^ old implementation ^^^   ///    vvv new implementation vvv
                cam_in(c)%cflx(i,c_i(2)) = data_flux_fuel%co2flx(i,1,c)
!--BEH  ^^^ new implementation ^^^
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


  end subroutine atm_import_moab

! endif for HAVE_MOAB
#endif



end module atm_comp_mct
