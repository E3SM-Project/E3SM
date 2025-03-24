module ocn_comp_mct

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
! !MODULE: ocn_comp_mct
! !INTERFACE:

! !DESCRIPTION:
!  This is the main driver for the Model for Predication Across Scales Ocean Model (MPASO).
!
! !USES:
   use mct_mod
   use esmf
   use seq_flds_mod
   use seq_cdata_mod
   use seq_infodata_mod
   use seq_timemgr_mod
   use seq_comm_mct,      only : seq_comm_suffix, seq_comm_inst, seq_comm_name, info_taskmap_comp
   use shr_file_mod
   use shr_cal_mod,       only : shr_cal_date2ymd
   use shr_sys_mod
   use shr_taskmap_mod,   only : shr_taskmap_write
   use shr_pio_mod
   use perf_mod

   use mpaso_cpl_indices
   use mpaso_mct_vars

   use mpas_framework
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_stream_manager
   use mpas_kind_types
   use mpas_io_units
   use mpas_timekeeping
   use mpas_bootstrapping
   use mpas_dmpar
   use mpas_constants
   use mpas_log
#ifdef HAVE_MOAB
   use shr_kind_mod    , only: cxx => SHR_KIND_CXX
   use mpas_moabmesh
   use seq_comm_mct, only: MPOID
   use seq_comm_mct,      only: num_moab_exports
   use iMOAB, only: iMOAB_DefineTagStorage, iMOAB_SetDoubleTagStorage, iMOAB_GetMeshInfo
   use seq_comm_mct, only:  seq_comm_compare_mb_mct
#endif
   use iso_c_binding, only : c_char, c_loc, c_ptr, c_int
   use mpas_c_interfacing, only : mpas_f_to_c_string, mpas_c_to_f_string

   use ocn_core
   use ocn_core_interface
   use ocn_forcing
   use ocn_constants
   use ocn_analysis_driver
   use ocn_time_integration
   use ocn_time_average_coupled
   use ocn_frazil_forcing
   use ocn_equation_of_state
   use ocn_surface_land_ice_fluxes
   use ocn_diagnostics
   use ocn_diagnostics_variables
   use ocn_tracer_short_wave_absorption
   use ocn_tracer_short_wave_absorption_variable
   use ocn_tracer_ecosys
   use ocn_tracer_CFC
   use ocn_tracer_surface_restoring
   use ocn_gm
   use ocn_config
   use ocn_submesoscale_eddies
   use ocn_eddy_parameterization_helpers
   use ocn_scaled_dismf
!
! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  public :: ocn_init_mct
  public :: ocn_run_mct
  public :: ocn_final_mct
  SAVE
  private                              ! By default make data private

!
! ! PUBLIC DATA:
!
! !REVISION HISTORY:
! Author: Doug Jacobsen
!
!EOP
! !PRIVATE MODULE FUNCTIONS:
  private :: ocn_export_mct
  private :: ocn_import_mct
  private :: ocn_SetGSMap_mct
  private :: ocn_domain_mct
!
! !PRIVATE MODULE VARIABLES
#ifdef HAVE_MOAB
  private :: ocn_export_moab
  private :: ocn_import_moab
  integer , private :: mblsize, totalmbls, totalmbls_r
  real (kind=RKIND) , allocatable, private :: o2x_om(:,:)

  real (kind=RKIND) , allocatable, private :: x2o_om(:,:)
  integer :: mpicom_moab
#endif


  integer, private ::   &
      my_task

   integer  :: nsend, nrecv

   logical  :: ocn_c2_glctf           ! .true.  => ocn to glc thermal forcing coupling on

   character(len=StrKIND) :: runtype, coupleTimeStamp

   type(seq_infodata_type), pointer :: infodata
   type (iosystem_desc_t), pointer :: io_system
   integer :: ocnLogUnit ! unit number for ocn log

   !! MPASO Datatypes
   type (core_type), pointer :: corelist => null()
   type (domain_type), pointer:: domain
   integer :: itimestep, &  ! time step number for MPAS
              ocn_cpl_dt    ! length of coupling interval in seconds - set by coupler/ESMF

!=======================================================================

contains

!***********************************************************************
!BOP
!
! !IROUTINE: ocn_init_mct
!
! !INTERFACE:
  subroutine ocn_init_mct( EClock, cdata_o, x2o_o, o2x_o, NLFilename )!{{{
!
! !DESCRIPTION:
! Initialize MPASO
!
! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock), intent(inout) :: EClock
    type(seq_cdata), intent(inout) :: cdata_o
    type(mct_aVect), intent(inout) :: x2o_o, o2x_o
    character(len=*), optional, intent(in) :: NLFilename ! Namelist filename

! !REVISION HISTORY:
! Author: Doug Jacobsen
!EOP
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    type (domain_type), pointer :: domain_ptr

    integer :: OCNID, mpicom_o, lsize, start_ymd, start_tod, start_year, start_day,   &
       start_month, start_hour, iyear, mpas_o_cpl_dt, pio_iotype
    integer :: shrloglev, shrlogunit !< shr log level and log unit
    integer :: ocnPossibleErrUnit !< unit number to reserve for if a err log needs to be opened

    type(mct_gsMap), pointer :: gsMap_o

    type(mct_gGrid), pointer :: dom_o

    integer :: errorCode  ! error code
    integer :: ioerr
    character(len=StrKIND) :: ioerrstr, scratchFName

    character(len=StrKIND) :: cpl_seq_option

    type (MPAS_Time_Type) :: currTime
    integer :: iam, npes, ierr, ierr_local
    integer :: iyear0, imonth0
    character(len=StrKIND)  :: starttype          ! infodata start type
    character(len=StrKIND)  :: timeStamp
    character(len=StrKIND)  :: nml_filename
    character(len=16) :: inst_suffix
    integer :: inst_index
    integer :: lbnum

    type (MPAS_Time_Type) :: alarmStartTime
    type (MPAS_TimeInterval_Type) :: alarmTimeStep
    type (block_type), pointer :: block_ptr

    type (mpas_pool_type), pointer :: meshPool, statePool, &
                                      forcingPool, &
                                      averagePool, scratchPool

    logical :: exists
    logical :: verbose_taskmap_output                ! true then use verbose task-to-node mapping format

    character(kind=c_char), dimension(StrKIND+1) :: c_filename       ! StrKIND+1 for C null-termination character
    integer(kind=c_int) :: c_comm
    integer(kind=c_int) :: c_ierr
    type (c_ptr) :: mgr_p
    character(len=StrKIND) :: mesh_stream
    character(len=StrKIND) :: mesh_filename
    character(len=StrKIND) :: mesh_filename_temp
    character(len=StrKIND) :: ref_time_temp
    character(len=StrKIND) :: filename_interval_temp
    character(kind=c_char), dimension(StrKIND+1) :: c_mesh_stream
    character(kind=c_char), dimension(StrKIND+1) :: c_mesh_filename_temp
    character(kind=c_char), dimension(StrKIND+1) :: c_ref_time_temp
    character(kind=c_char), dimension(StrKIND+1) :: c_filename_interval_temp
    character(kind=c_char), dimension(StrKIND+1) :: c_iotype
    character(len=8) :: c_inst_index ! instance number
    character(len=8) :: c_npes       ! number of pes
    type (MPAS_Time_type) :: start_time
    type (MPAS_Time_type) :: ref_time
    type (MPAS_TimeInterval_type) :: filename_interval
    type (MPAS_TimeInterval_type) :: denInterval, remInterval, zeroInterval
    integer (kind=I8KIND) :: numDivs
    character(len=StrKIND) :: start_timestamp
    character(len=StrKIND) :: iotype
    logical :: streamsExists
    integer :: mesh_iotype
    logical :: ocn_c2_glcshelf

    logical, pointer :: tempLogicalConfig
    character(len=StrKIND), pointer :: tempCharConfig

    logical, pointer :: config_use_ecosysTracers
    logical, pointer :: config_use_CFCTracers
    logical, pointer :: config_use_activeTracers_surface_restoring
    logical, pointer :: config_use_surface_salinity_monthly_restoring
    logical, pointer :: config_scale_dismf_by_removed_ice_runoff
    character (len=StrKIND), pointer :: config_land_ice_flux_mode
    character (len=StrKIND), pointer :: config_glc_thermal_forcing_coupling_mode

    ! ssh coupling interval initialization
    integer, pointer :: index_avgZonalSSHGradient, index_avgMeridionalSSHGradient
    real (kind=RKIND), dimension(:,:), pointer :: avgSSHGradient
    real (kind=RKIND), pointer :: &
       runningMeanRemovedIceRunoff ! the area integrated, running mean of removed ice runoff from the ocean

#ifdef HAVE_MOAB
    character*100 outfile, wopts
    integer :: ierrmb, numco, tagtype, tagindex, ent_type
    character(CXX) :: tagname
#ifdef MOABCOMP
    real(r8)                 :: difference
    type(mct_list) :: temp_list
    integer :: size_list, index_list
    type(mct_string)    :: mctOStr  !
    character(CXX) :: mct_field, modelStr
#endif

#endif
    interface
       subroutine xml_stream_parser(xmlname, mgr_p, comm, ierr) bind(c)
          use iso_c_binding, only : c_char, c_ptr, c_int
          character(kind=c_char), dimension(*), intent(in) :: xmlname
          type (c_ptr), intent(inout) :: mgr_p
          integer(kind=c_int), intent(inout) :: comm
          integer(kind=c_int), intent(out) :: ierr
       end subroutine xml_stream_parser

       subroutine xml_stream_get_attributes(xmlname, streamname, comm, filename, ref_time, filename_interval, io_type, ierr) bind(c)
          use iso_c_binding, only : c_char, c_int
          character(kind=c_char), dimension(*), intent(in) :: xmlname
          character(kind=c_char), dimension(*), intent(in) :: streamname
          integer(kind=c_int), intent(inout) :: comm
          character(kind=c_char), dimension(*), intent(out) :: filename
          character(kind=c_char), dimension(*), intent(out) :: ref_time
          character(kind=c_char), dimension(*), intent(out) :: filename_interval
          character(kind=c_char), dimension(*), intent(out) :: io_type
          integer(kind=c_int), intent(out) :: ierr
       end subroutine xml_stream_get_attributes
    end interface

!-----------------------------------------------------------------------
!
!   setup mpaso data structures
!
!-----------------------------------------------------------------------
    allocate(corelist)
    nullify(corelist % next)

    allocate(corelist % domainlist)
    nullify(corelist % domainlist % next)

    domain => corelist % domainlist
    domain % core => corelist

    call mpas_allocate_domain(domain)

!-----------------------------------------------------------------------
!
!   set cdata pointers
!
!-----------------------------------------------------------------------
    errorCode = 0

    call seq_cdata_setptrs(cdata_o, ID=OCNID, mpicom=mpicom_o, &
         gsMap=gsMap_o, dom=dom_o, infodata=infodata)

    MPASO_MCT_OCNID =  OCNID
    MPASO_MCT_gsMap_o => gsMap_o
    MPASO_MCT_dom_o   => dom_o

    call MPI_comm_rank(mpicom_o,iam,ierr)
#if (defined _MEMTRACE)
    if(iam == 0) then
        lbnum=1
        call memmon_dump_fort('memmon.out','ocn_init_mct:start::',lbnum)
    endif
#endif

    ! Determine coupling type
    call seq_infodata_GetData(infodata, cpl_seq_option=cpl_seq_option)

    ! Determine if ocn to glc thermal forcing coupling is on
    call seq_infodata_GetData(infodata, ocn_c2_glctf=ocn_c2_glctf)

!-----------------------------------------------------------------------
!
!   initialize the model run
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!   first initializaiton phase of mpaso
!   initialize mpaso because grid information is needed for
!   creation of GSMap_ocn.
!   call mpaso initialization routines
!
!-----------------------------------------------------------------------

    call t_startf('mpaso_init1')

    io_system => shr_pio_getiosys(ocnid)

    pio_iotype = shr_pio_getiotype(ocnid)
!    call MPAS_io_set_iotype(domain % iocontext, pio_iotype)
    call MPAS_io_unset_iotype(domain % iocontext)

    ! ----------------
    ! Set up log file information
    ! ----------------
    inst_suffix = seq_comm_suffix(ocnID) ! a suffix to append to log file name
    ocnLogUnit = shr_file_getUnit() ! reserve unit number for log unit
    ocnPossibleErrUnit = shr_file_getUnit() ! reserve unit number for possible error log file

    ! Note: In following code, a file ocn_modelio.nml is queried to determine the log file name to use.
    ! This file is generated when a case is created.  The default file name in there is currently ocn.log.DATE-TIME.
    if (iam==0) then
       inquire(file='ocn_modelio.nml'//trim(inst_suffix),exist=exists)
       call shr_file_setio('ocn_modelio.nml'//trim(inst_suffix), ocnLogUnit)
    endif

    ! Store shr log unit and level so we can reassign them later when OCN control is complete
    call shr_file_getLogUnit(shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    ! Set the shr log unit to be the ocnLogUnit
    call shr_file_setLogUnit(ocnLogUnit)

    domain_ptr => domain

    ! Write to ocnLogUnit here, because the log module is not open yet.
    if (iam==0) write(ocnLogUnit,'(a,i6)') '=== Beginning ocn_init_mct: rank=',iam

    call mpas_framework_init_phase1(domain % dminfo, mpicom_o)

    ! Identify SMP nodes and process/SMP mapping for this instance
    ! (Assume that processor names are SMP node names on SMP clusters.)
    inst_index = seq_comm_inst(ocnID)
    write(c_inst_index,'(i8)') inst_index

    if (info_taskmap_comp > 0) then

       if (info_taskmap_comp == 1) then
          verbose_taskmap_output = .false.
       else
          verbose_taskmap_output = .true.
       endif

       call MPI_comm_size(mpicom_o,npes,ierr)
       write(c_npes,'(i8)') npes

       if (iam==0) then
          write(ocnLogUnit,'(/,3A)') &
             trim(adjustl(c_npes)), &
             ' pes participating in computation of MPAS-OCEAN instance #', &
             trim(adjustl(c_inst_index))
          call shr_sys_flush(ocnLogUnit)
       endif

       call t_startf("shr_taskmap_write")
       call shr_taskmap_write(ocnLogUnit, mpicom_o,                 &
                              'OCN #'//trim(adjustl(c_inst_index)), &
                              verbose=verbose_taskmap_output        )
       call t_stopf("shr_taskmap_write")

    endif

    ! Setup function pointers for MPASO core and domain types
    call ocn_setup_core(corelist)
    call ocn_setup_domain(domain_ptr)

    ! ===========
    ! Initialize log manager
    call mpas_log_init(domain % logInfo, domain, unitNumbers=(/ocnLogUnit, ocnPossibleErrUnit/), err=ierr)
    if ( ierr /= 0 ) then
       write(ocnLogUnit,*) 'ERROR: log init failed for core ' // trim(domain % core % coreName)
       call mpas_dmpar_abort(domain % dminfo)
    end if

    ! Set core specific options here
    ! Disable output from all but the master task for E3SM!
    ! (This overrides the default set by mpas_log_init based on MPAS_DEBUG setting.)
    if (iam /= 0) then
       domain % logInfo % outputLog % isActive = .false.
    endif

    ! After core has had a chance to modify log defaults, open the output log
    call mpas_log_open(err=ierr)
    if ( ierr /= 0 ) then
       write(ocnLogUnit,*) 'ERROR: log open failed for core ' // trim(domain % core % coreName)
       call mpas_dmpar_abort(domain % dminfo)
    end if
    ! ===========


    ! ----------
    ! Process namelist and streams files
    ! ----------
    ! Override the names of the stream and namelist files
    domain_ptr % namelist_filename = 'mpaso_in' // trim(inst_suffix)
    domain_ptr % streams_filename = 'streams.ocean' // trim(inst_suffix)

    ! Setup namelist variables, and read the namelist
    ierr = domain_ptr % core % setup_namelist(domain_ptr % configs, domain_ptr % namelist_filename, domain_ptr % dminfo)
    if ( ierr /= 0 ) then
       call mpas_log_write('Namelist setup failed for core ' // trim(domain % core % coreName), MPAS_LOG_CRIT)
    end if

    call mpas_framework_init_phase2(domain_ptr, io_system)

    ! Define package variables
    ierr = domain_ptr % core % define_packages(domain_ptr % packages)
    if ( ierr /= 0 ) then
       call mpas_log_write('Package definition failed for core ' // trim(domain % core % coreName), MPAS_LOG_CRIT)
    end if

    ! Setup packages (i.e. determine if they should be on or off)
    ierr = domain_ptr % core % setup_packages(domain_ptr % configs, domain_ptr % packages, domain_ptr % ioContext)
    if ( ierr /= 0 ) then
       call mpas_log_write('Package setup failed for core ' // trim(domain % core % coreName), MPAS_LOG_CRIT)
    end if

    ! Setup decompositions available for dimensions
    ierr = domain_ptr % core % setup_decompositions(domain_ptr % decompositions)
    if ( ierr /= 0 ) then
       call mpas_log_write('Decomposition setup failed for core ' // trim(domain % core % coreName), MPAS_LOG_CRIT)
    end if


    ! ----------
    ! Override namelist options based on start type
    ! ----------
    call mpaso_cpl_indices_set()
    call seq_infodata_GetData( infodata, start_type=starttype)
    if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
       runtype = "initial"
    else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
       runtype = "continue"
    else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
       runtype = "branch"
    else
       call mpas_log_write('ocn_comp_mct ERROR: infodata provided unknown starttype: ' // trim(starttype), MPAS_LOG_CRIT)
    end if
    if (runtype == "initial") then ! Start up run
        ! Turn off restart
        call mpas_pool_get_config(domain_ptr % configs, "config_do_restart", tempLogicalConfig)
        tempLogicalConfig = .false.

        ! Setup start time. Will be over written later when clocks are synchronized
        call mpas_pool_get_config(domain_ptr % configs, "config_start_time", tempCharConfig)
        tempCharConfig = trim(tempCharConfig) // "0:00:00"

        ! Setup run duration. Will be ignored in coupled run, since coupler defines how long the run is.
        call mpas_pool_get_config(domain_ptr % configs, "config_run_duration", tempCharConfig)
        tempCharConfig = "0001-00-00_00:00:00"
    else if (runtype == "continue" .or. runtype == "branch") then ! Restart run or branch run
        ! Turn on restart
        call mpas_pool_get_config(domain_ptr % configs, "config_do_restart", tempLogicalConfig)
        tempLogicalConfig = .true.

        ! Set start time to be read from file
        call mpas_pool_get_config(domain_ptr % configs, "config_start_time", tempCharConfig)
        tempCharConfig = "file"

        ! Setup run duration. Will be ignored in coupled run, since coupler defines how long the run is.
        call mpas_pool_get_config(domain_ptr % configs, "config_run_duration", tempCharConfig)
        tempCharConfig = "0001-00-00_00:00:00"
    end if

    ! Setup MPASO simulation clock
    ierr = domain_ptr % core % setup_clock(domain_ptr % clock, domain_ptr % configs)
    if ( ierr /= 0 ) then
       call mpas_log_write('Clock setup failed for core ' // trim(domain % core % coreName), MPAS_LOG_CRIT)
    end if

    call mpas_log_write('Reading streams configuration from file '//trim(domain % streams_filename))
    inquire(file=trim(domain % streams_filename), exist=streamsExists)

    if ( .not. streamsExists ) then
       call mpas_log_write('Streams file '//trim(domain % streams_filename)//' does not exist.', MPAS_LOG_CRIT)
    end if

    !
    ! Using information from the namelist, a graph.info file, and a file containing
    !    mesh fields, build halos and allocate blocks in the domain
    !
    ierr = domain_ptr % core % get_mesh_stream(domain_ptr % configs, mesh_stream)
    if ( ierr /= 0 ) then
       call mpas_log_write('Failed to find mesh stream for core ' // trim(domain % core % coreName), MPAS_LOG_CRIT)
    end if


    call mpas_f_to_c_string(domain % streams_filename, c_filename)
    call mpas_f_to_c_string(mesh_stream, c_mesh_stream)
    c_comm = domain_ptr % dminfo % comm
    call xml_stream_get_attributes(c_filename, c_mesh_stream, c_comm, &
                                   c_mesh_filename_temp, c_ref_time_temp, &
                                   c_filename_interval_temp, c_iotype, c_ierr)
    if (c_ierr /= 0) then
       call mpas_log_write('xml_stream_get_attributes failed.', MPAS_LOG_CRIT)
    end if
    call mpas_c_to_f_string(c_mesh_filename_temp, mesh_filename_temp)
    call mpas_c_to_f_string(c_ref_time_temp, ref_time_temp)
    call mpas_c_to_f_string(c_filename_interval_temp, filename_interval_temp)
    call mpas_c_to_f_string(c_iotype, iotype)

    if (trim(iotype) == 'pnetcdf') then
       mesh_iotype = MPAS_IO_PNETCDF
    else if (trim(iotype) == 'pnetcdf,cdf5') then
       mesh_iotype = MPAS_IO_PNETCDF5
    else if (trim(iotype) == 'netcdf') then
       mesh_iotype = MPAS_IO_NETCDF
    else if (trim(iotype) == 'netcdf4') then
       mesh_iotype = MPAS_IO_NETCDF4
    else
       mesh_iotype = MPAS_IO_PNETCDF
    end if

    start_time = mpas_get_clock_time(domain_ptr % clock, MPAS_START_TIME, ierr)
    if ( trim(ref_time_temp) == 'initial_time' ) then
        call mpas_get_time(start_time, dateTimeString=ref_time_temp, ierr=ierr)
    end if

    if ( trim(filename_interval_temp) == 'none' ) then
        call mpas_expand_string(ref_time_temp, -1, mesh_filename_temp, mesh_filename)
    else
        call mpas_set_time(ref_time, dateTimeString=ref_time_temp, ierr=ierr)
        call mpas_set_timeInterval(filename_interval, timeString=filename_interval_temp, ierr=ierr)
        call mpas_build_stream_filename(ref_time, start_time, filename_interval, mesh_filename_temp, -1, mesh_filename, ierr)
    end if

    ! Bootstrap framework (1). Here data structures are setup, but dimensions and arrays are not finalized.
    call mpas_log_write(' ** Attempting to bootstrap MPAS framework using stream: ' // trim(mesh_stream))
    call mpas_bootstrap_framework_phase1(domain_ptr, mesh_filename, mesh_iotype)

    !
    ! Set up run-time streams
    !
    call MPAS_stream_mgr_init(domain_ptr % streamManager, domain_ptr % ioContext, domain_ptr % clock, &
                              domain_ptr % blocklist % allFields, domain_ptr % packages, domain_ptr % blocklist % allStructs)

    call add_stream_attributes(domain_ptr)

    ! Setup all immutable streams for the core
    ierr = domain_ptr % core % setup_immutable_streams(domain_ptr % streamManager)
    if ( ierr /= 0 ) then
       call mpas_log_write('Immutable streams setup failed for core ' // trim(domain % core % coreName), MPAS_LOG_CRIT)
    end if

    ! Parse / read all streams configuration
    mgr_p = c_loc(domain_ptr % streamManager)
    call xml_stream_parser(c_filename, mgr_p, c_comm, c_ierr)
    if (c_ierr /= 0) then
       call mpas_log_write('xml_stream_parser failed.', MPAS_LOG_CRIT)
    end if

    my_task = domain_ptr % dminfo % my_proc_id

    !
    ! Finalize the setup of blocks and fields
    !
    call mpas_bootstrap_framework_phase2(domain_ptr)
    call t_stopf('mpaso_init1')

    call t_startf('mpaso_init2')
    ! Initialize the MPASO core
    ierr = domain_ptr % core % core_init(domain_ptr, timeStamp)
    if ( ierr /= 0 ) then
       call mpas_log_write('Core init failed for core ' // trim(domain % core % coreName), MPAS_LOG_CRIT)
    end if
    call t_stopf('mpaso_init2')


!-----------------------------------------------------------------------
!
!   initialize time-stamp information
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!   check for consistency of mpaso and sync clock initial time
!
!-----------------------------------------------------------------------

    if (runtype == 'initial') then
       call seq_timemgr_EClockGetData(EClock, ECurrTime=currTime % t)
       call mpas_set_clock_time(domain_ptr % clock, currTime, MPAS_START_TIME, ierr)
       call mpas_set_clock_time(domain_ptr % clock, currTime, MPAS_NOW, ierr)
    else if (runtype == 'continue' .or. runtype == 'branch') then
       call seq_timemgr_EClockGetData(EClock, ECurrTime=currTime % t)
       call mpas_set_clock_time(domain_ptr % clock, currTime, MPAS_START_TIME, ierr)
       call mpas_set_clock_time(domain_ptr % clock, currTime, MPAS_NOW, ierr)
    end if
    ! Doublecheck that clocks are synced here.
    ! (Note that if this is an initial run, a section below will advance the ocean clock
    ! by one coupling interval, at which point we expect the clocks to be OUT of sync.)
    call check_clocks_sync(domain % clock, Eclock, ierr)

!-----------------------------------------------------------------------
!
!   initialize MCT attribute vectors and indices
!
!-----------------------------------------------------------------------

    call t_startf ('mpaso_mct_init')

    call ocn_SetGSMap_mct( mpicom_o, OCNID, GSMap_o )
    lsize = mct_gsMap_lsize(gsMap_o, mpicom_o)

    ! Initialize mct ocn domain (needs ocn initialization info)
    call ocn_domain_mct( lsize, gsMap_o, dom_o )

    ! Inialize mct attribute vectors

    call mct_aVect_init(x2o_o, rList=seq_flds_x2o_fields, lsize=lsize)
    call mct_aVect_zero(x2o_o)

    call mct_aVect_init(o2x_o, rList=seq_flds_o2x_fields, lsize=lsize)
    call mct_aVect_zero(o2x_o)

    nsend = mct_avect_nRattr(o2x_o)
    nrecv = mct_avect_nRattr(x2o_o)
#ifdef HAVE_MOAB
    call init_moab_mpas(domain_ptr, OCNID, MPOID) ! should return MPOID ..
    call mpas_log_write('initialized MOAB MPAS ocean instance... ')

    ! initialize moab tag fields array
   mblsize = lsize
   totalmbls = mblsize * nsend ! size of the double array for exporting to coupler
   allocate (o2x_om(lsize, nsend) )
   o2x_om = 0._r8

  ! define tags according to the seq_flds_o2x_fields
   tagtype = 1  ! dense, double
   numco = 1 !  one value per cell / entity
   tagname = trim(seq_flds_o2x_fields)//C_NULL_CHAR
   ierrmb = iMOAB_DefineTagStorage(MPOID, tagname, tagtype, numco,  tagindex )
   if ( ierrmb == 1 ) then
      call mpas_log_write('cannot define tags for MOAB o2x fields ' // trim(seq_flds_o2x_fields), MPAS_LOG_ERR)
   endif
   !  set all tags to 0 in one method
   ent_type = 1! cells
   ierr = iMOAB_SetDoubleTagStorage ( MPOID, tagname, totalmbls , ent_type, o2x_om(1, 1) )
   if ( ierr /= 0 ) then
      write(ocnLogUnit,*) 'Fail to set MOAB fields '
   endif

  ! define tags according to the seq_flds_x2o_fields
   totalmbls_r = mblsize * nrecv ! size of the double array for importing
   allocate (x2o_om(lsize, nrecv) )
   x2o_om = 0._r8

   tagname = trim(seq_flds_x2o_fields)//C_NULL_CHAR
   ierrmb = iMOAB_DefineTagStorage(MPOID, tagname, tagtype, numco,  tagindex )
   if ( ierrmb == 1 ) then
      call mpas_log_write('cannot define tags for MOAB x2o fields ' // trim(seq_flds_x2o_fields), MPAS_LOG_ERR)
   endif
   !  set all tags to 0 in one method
   ierr = iMOAB_SetDoubleTagStorage ( MPOID, tagname, totalmbls_r , ent_type, x2o_om(1, 1) )
   if ( ierr /= 0 ) then
      write(ocnLogUnit,*) 'Fail to set MOAB fields '
   endif

   ! add domain tags
   tagname=trim(seq_flds_dom_fields)//C_NULL_CHAR
   tagtype = 1
   ierr = iMOAB_DefineTagStorage(MPOID, tagname, tagtype, numco,  tagindex )
   if ( ierr /= 0 ) then
      write(ocnLogUnit,*) 'Fail to set define dom fields '
   endif
   ent_type = 1 ! cells

   call ocn_domain_moab(MPOID)
#endif
!-----------------------------------------------------------------------
!
!   initialize necessary  coupling info
!
!-----------------------------------------------------------------------

    call seq_timemgr_EClockGetData(EClock, dtime=ocn_cpl_dt)
    call mpas_set_timeInterval(alarmTimeStep, S=ocn_cpl_dt, ierr=ierr)
    call mpas_get_timeInterval(alarmTimeStep, timeString=coupleTimeStamp, ierr=ierr)

    ! Verify the mpas time step fits into a coupling interval
    call mpas_pool_get_config(domain % configs, 'config_dt', tempCharConfig)
    call mpas_set_timeInterval(denInterval, timeString=tempCharConfig, ierr=ierr)
    call mpas_set_timeInterval(zeroInterval, S=0, ierr=ierr)
    call mpas_interval_division(start_time, alarmTimeStep, denInterval, numDivs, remInterval)

    ierr = 0

    if ( alarmTimeStep < denInterval ) then
       ierr = 1
    end if
    ierr_local = ierr
    call mpas_dmpar_max_int(domain % dminfo, ierr_local, ierr)

    if ( ierr == 1 ) then
       call mpas_log_write('Coupling interval is: ' // trim(coupleTimeStamp), MPAS_LOG_ERR)
       call mpas_log_write('        OCN Model time step is: ' // trim(tempCharConfig), MPAS_LOG_ERR)
       call mpas_log_write('        The model time step cannot be longer then the coupling interval', MPAS_LOG_ERR)
       call mpas_log_write('Model is not properly configured for coupling interval.', MPAS_LOG_CRIT)
    end if

    if ( remInterval > zeroInterval ) then
       ierr = 1
    end if

    ierr_local = ierr
    call mpas_dmpar_max_int(domain % dminfo, ierr_local, ierr)

    if ( ierr == 1 ) then
       call mpas_log_write('Coupling interval is: ' // trim(coupleTimeStamp), MPAS_LOG_ERR)
       call mpas_log_write('        OCN Model time step is: ' // trim(tempCharConfig), MPAS_LOG_ERR)
       call mpas_log_write('        These are not synchronized, so time steps will not match to coupling interval boundaries.', MPAS_LOG_ERR)
       call mpas_log_write('        Please reconfigure either the coupling interval or the time step.', MPAS_LOG_ERR)
       call mpas_log_write('Model is not properly configured for coupling interval.', MPAS_LOG_CRIT)
    end if

    ! set coupling alarm
    alarmStartTime = currTime
    call mpas_add_clock_alarm(domain % clock, coupleAlarmID, alarmStartTime, alarmTimeStep, ierr=ierr)
    call mpas_reset_clock_alarm(domain_ptr % clock, coupleAlarmID, ierr=ierr)

    ! Build forcing arrays.
    block_ptr => domain % blocklist
    do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
        call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)

        call ocn_forcing_build_fraction_absorbed_array(meshPool, statePool, forcingPool, ierr, 1)
        call ocn_frazil_forcing_build_arrays(domain, meshPool, forcingPool, statePool, ierr)
        call mpas_timer_start("land_ice_build_arrays", .false.)
        call ocn_surface_land_ice_fluxes_build_arrays(meshPool, forcingPool, scratchPool, &
                                                      statePool, err=ierr)
        call mpas_timer_stop("land_ice_build_arrays")

        block_ptr => block_ptr % next
    end do

    currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
    call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
    call mpas_log_write( 'Initial time '//trim(timeStamp))

    ! read initial data required for variable shortwave
    call mpas_timer_start('io_shortwave',.false.)
    call ocn_get_shortWaveData(domain_ptr % streamManager, domain, domain % clock, .true.)
    call mpas_timer_stop('io_shortwave')

    ! read initial data required for ecosys forcing
    call mpas_pool_get_config(domain % configs, 'config_use_ecosysTracers', config_use_ecosysTracers)
    if (config_use_ecosysTracers) then
        call mpas_timer_start('io_ecosys',.false.)
        call ocn_get_ecosysData(domain_ptr % streamManager, domain, domain % clock, .true.)
        call mpas_timer_stop('io_ecosys')
    endif

    ! read initial data required for CFC forcing
    call mpas_pool_get_config(domain % configs, 'config_use_CFCTracers', config_use_CFCTracers)
    if (config_use_CFCTracers) then
        call mpas_timer_start('io_CFC',.false.)
        call ocn_get_CFCData(domain_ptr % streamManager, domain, domain % clock, .true.)
        call mpas_timer_stop('io_CFC')
    endif

    ! read initial data required for monthly surface salinity restoring
    call mpas_pool_get_config(domain % configs, 'config_use_activeTracers_surface_restoring',  &
                                                 config_use_activeTracers_surface_restoring)
    call mpas_pool_get_config(domain % configs, 'config_use_surface_salinity_monthly_restoring',  &
                                                 config_use_surface_salinity_monthly_restoring)
    if (config_use_activeTracers_surface_restoring .and.  config_use_surface_salinity_monthly_restoring) then
        call mpas_timer_start('io_monthly_surface_salinity',.false.)
        call ocn_get_surfaceSalinityData(domain_ptr % streamManager, domain, domain % clock, .true.)
        call mpas_timer_stop('io_monthly_surface_salinity')
    endif

    ! Setup clock for initial runs
    if ( runtype == 'initial' ) then
       ! Advance clock one coupling interval to be in sync with the coupler clock.
       do while (.not. mpas_is_alarm_ringing(domain_ptr % clock, coupleAlarmID, ierr=ierr))
          itimestep = itimestep + 1
          call mpas_advance_clock(domain_ptr % clock)
       end do

       block_ptr => domain % blocklist
       do while(associated(block_ptr))
          call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
          call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)

          call ocn_time_average_coupled_init(forcingPool)
          call ocn_time_average_coupled_accumulate(statePool, forcingPool, 1)
          block_ptr => block_ptr % next
       end do

       ! initialize scaled data ice-shelf melt fluxes based on remove ice runoff
       call ocn_init_scaled_dismf(domain)
    end if

!-----------------------------------------------------------------------
!
!   initialize coupling variables
!   NOTE: could be moved to subroutine
!
!-----------------------------------------------------------------------

    if ( runtype == "initial") then
       block_ptr => domain % blocklist
       do while(associated(block_ptr))
          call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
          call mpas_pool_get_dimension(forcingPool, 'index_avgSSHGradientZonal', index_avgZonalSSHGradient)
          call mpas_pool_get_dimension(forcingPool, 'index_avgSSHGradientMeridional', index_avgMeridionalSSHGradient)
          call mpas_pool_get_array(forcingPool, 'avgSSHGradient', avgSSHGradient)
          block_ptr => block_ptr % next
       end do
    endif

!-----------------------------------------------------------------------
!
!   send intial state to driver
!
!-----------------------------------------------------------------------

    call ocn_export_mct(o2x_o, errorCode)
    if (errorCode /= 0) then
       call mpas_log_write('ERROR in ocn_export_mct', MPAS_LOG_CRIT)
    endif
#ifdef HAVE_MOAB
    ! initial state has to be exported to moab too
    call ocn_export_moab(EClock)
#endif
    ! Setup clock for initial runs
    if (runtype == "continue" .or. runtype == "branch" ) then
       block_ptr => domain % blocklist
       do while(associated(block_ptr))
          call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
          call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)

          call ocn_time_average_coupled_init(forcingPool)
          call ocn_time_average_coupled_accumulate(statePool, forcingPool, 1)
          block_ptr => block_ptr % next
       end do

      ! initialize scaled data ice-shelf melt fluxes based on remove ice runoff
       call ocn_init_scaled_dismf(domain)
    end if

    call t_stopf ('mpaso_mct_init')

    call mpas_pool_get_config(domain % configs, 'config_land_ice_flux_mode', config_land_ice_flux_mode)
    if ( trim(config_land_ice_flux_mode) == 'off' .or. &
         trim(config_land_ice_flux_mode) == 'pressure_only' .or. &
         trim(config_land_ice_flux_mode) == 'data' .or. &
         trim(config_land_ice_flux_mode) == 'standalone' ) then
       ocn_c2_glcshelf = .false.
    else if ( trim(config_land_ice_flux_mode) .eq. 'coupled' ) then
       ocn_c2_glcshelf = .true.
    else
       call mpas_log_write('ERROR: unknown land_ice_flux_mode: ' // trim(config_land_ice_flux_mode), MPAS_LOG_CRIT)
    end if
    call seq_infodata_PutData(infodata, ocn_prognostic=.true., ocnrof_prognostic=.true., &
                              ocn_c2_glcshelf=ocn_c2_glcshelf)

    call mpas_pool_get_config(domain % configs, 'config_scale_dismf_by_removed_ice_runoff', &
                              config_scale_dismf_by_removed_ice_runoff)
    if (config_scale_dismf_by_removed_ice_runoff) then
       ! independent of space so should be no need to loop over blocks
       block_ptr => domain % blocklist
       call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
       call MPAS_pool_get_array(forcingPool, "runningMeanRemovedIceRunoff", &
                                runningMeanRemovedIceRunoff)
       call seq_infodata_PutData(infodata, rmean_rmv_ice_runoff=runningMeanRemovedIceRunoff)
    end if

    call mpas_pool_get_config(domain % configs, 'config_glc_thermal_forcing_coupling_mode', config_glc_thermal_forcing_coupling_mode)
    if ( trim(config_glc_thermal_forcing_coupling_mode) == 'off' ) then
       call seq_infodata_PutData(infodata, ocn_c2_glctf=.false.)
    else if ( trim(config_glc_thermal_forcing_coupling_mode) == '2d' ) then
       call seq_infodata_PutData(infodata, ocn_c2_glctf=.true.)
    else
       call mpas_log_write('ERROR: unknown config_glc_thermal_forcing_coupling_mode: ' // &
               trim(config_glc_thermal_forcing_coupling_mode), MPAS_LOG_CRIT)
    end if

!-----------------------------------------------------------------------
!
!   get initial state from driver
!
!-----------------------------------------------------------------------

#ifdef HAVE_MOAB
#ifdef MOABCOMP
    ! loop over all fields in seq_flds_x2o_fields
    call mct_list_init(temp_list ,seq_flds_x2o_fields)
    size_list=mct_list_nitem (temp_list)
    ent_type = 1 ! entity type is cell for ocn
    modelStr='ocn init'
    mpicom_moab = mpicom_o ! save it for run method
    do index_list = 1, size_list
      call mct_list_get(mctOStr,index_list,temp_list)
      mct_field = mct_string_toChar(mctOStr)
      tagname= trim(mct_field)//C_NULL_CHAR
      call seq_comm_compare_mb_mct(modelStr, mpicom_moab, x2o_o, mct_field,  MPOID, tagname, ent_type, difference)
    enddo
    call mct_list_clean(temp_list)
#endif
#endif

    call ocn_import_mct(x2o_o, errorCode)
    if (errorCode /= 0) then
       call mpas_log_write('Error in ocn_import_mct', MPAS_LOG_CRIT)
    endif

#ifdef HAVE_MOAB

    call ocn_import_moab(Eclock, errorCode)
    if (errorCode /= 0) then
       call mpas_log_write('Error in ocn_import_moab', MPAS_LOG_CRIT)
    endif
#endif


    itimestep = 0

    ! Reset all output alarms, to prevent intial time step from writing any output, unless it's ringing.
    call mpas_stream_mgr_reset_alarms(domain_ptr % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)
    call mpas_stream_mgr_reset_alarms(domain_ptr % streamManager, direction=MPAS_STREAM_INPUT, ierr=ierr)

!----------------------------------------------------------------------------
!   Reset shr logging to original values
!----------------------------------------------------------------------------
    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

#if defined (_MEMTRACE)
    if(iam  == 0) then
        lbnum=1
        call memmon_dump_fort('memmon.out','ocn_init_mct:end::',lbnum)
        call memmon_reset_addr()
    endif
#endif

    call mpas_log_write('=== Completed ocn_init_mct ===', flushNow=.true.)

!-----------------------------------------------------------------------
!EOC

 end subroutine ocn_init_mct!}}}

!***********************************************************************
!BOP
!
! !IROUTINE: ocn_run_mct
!
! !INTERFACE:
  subroutine ocn_run_mct( EClock, cdata_o, x2o_o, o2x_o)!{{{

!
! !DESCRIPTION:
! Run MPASO for one coupling interval
!
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_Clock)            , intent(inout) :: EClock
      type(seq_cdata)             , intent(inout) :: cdata_o
      type(mct_aVect)             , intent(inout) :: x2o_o
      type(mct_aVect)             , intent(inout) :: o2x_o

!! !REVISION HISTORY:
!! Author: Doug Jacobsen
!!EOP
!!-----------------------------------------------------------------------
!!
!!  local variables
!!
!!-----------------------------------------------------------------------
      integer :: ymd, tod, ihour, iminute, isecond
      integer :: iyear, imonth, iday, curr_ymd, curr_tod
      integer :: shrloglev, shrlogunit
      real (kind=RKIND) :: dt, current_wallclock_time

      type (block_type), pointer :: block_ptr

      type (mpas_pool_type), pointer :: meshPool, statePool, forcingPool, averagePool, scratchPool

      type (MPAS_Time_Type) :: currTime
      type (domain_type), pointer :: domain_ptr
      character(len=StrKIND) :: timeStamp, streamName, WCstring
      type (MPAS_timeInterval_type) :: timeStep
      integer :: iam, ierr, streamDirection
      logical :: streamActive, debugOn
      logical, pointer :: config_write_output_on_startup
      logical, pointer :: config_use_ecosysTracers
      logical, pointer :: config_use_CFCTracers
      logical, pointer :: config_use_activeTracers_surface_restoring
      logical, pointer :: config_use_surface_salinity_monthly_restoring
      logical, pointer :: config_scale_dismf_by_removed_ice_runoff
      character (len=StrKIND), pointer :: config_restart_timestamp_name
      character (len=StrKIND), pointer :: config_sw_absorption_type

      ! Added for coupling interval initialization
      integer, pointer :: index_avgZonalSSHGradient, index_avgMeridionalSSHGradient
      real (kind=RKIND), dimension(:,:), pointer :: avgSSHGradient
      real (kind=RKIND), pointer :: &
         runningMeanRemovedIceRunoff ! the area integrated, running mean of removed ice runoff from the ocean

#ifdef HAVE_MOAB
#ifdef MOABCOMP
    real(r8)                 :: difference
    type(mct_list) :: temp_list
    integer :: size_list, index_list, ent_type
    type(mct_string)    :: mctOStr  !
    character(CXX) :: mct_field, modelStr, tagname
#endif
#endif

      iam = domain % dminfo % my_proc_id

      debugOn = .false.
#ifdef MPAS_DEBUG
      debugOn = .true.
#endif

      domain_ptr => domain

      ! Set MPAS Log module instance
      mpas_log_info => domain % logInfo
      if (debugOn) call mpas_log_write("=== Beginning ocn_run_mct ===")

      call mpas_pool_get_config(domain_ptr % configs, 'config_restart_timestamp_name', config_restart_timestamp_name)
      call mpas_pool_get_config(domain_ptr % configs, 'config_sw_absorption_type', config_sw_absorption_type)

      ! Setup log information.
      call shr_file_getLogUnit (shrlogunit)
      call shr_file_getLogLevel(shrloglev)
      call shr_file_setLogUnit (ocnLogUnit)

      timeStep = mpas_get_clock_timestep(domain_ptr % clock, ierr=ierr)
      call mpas_get_timeInterval(timeStep, dt=dt)
      call mpas_reset_clock_alarm(domain_ptr % clock, coupleAlarmID, ierr=ierr)


#ifdef HAVE_MOAB
#ifdef MOABCOMP
    ! loop over all fields in seq_flds_x2o_fields
    call mct_list_init(temp_list ,seq_flds_x2o_fields)
    size_list=mct_list_nitem (temp_list)
    ent_type = 1 ! entity type is cell for ocn
    modelStr='ocn run'
    do index_list = 1, size_list
      call mct_list_get(mctOStr,index_list,temp_list)
      mct_field = mct_string_toChar(mctOStr)
      tagname= trim(mct_field)//C_NULL_CHAR
      call seq_comm_compare_mb_mct(modelStr, mpicom_moab, x2o_o, mct_field,  MPOID, tagname, ent_type, difference)
    enddo
    call mct_list_clean(temp_list)
#endif
#endif

      ! Import state from coupler
      call ocn_import_mct(x2o_o, ierr)
      ! Import state from moab coupler
#ifdef HAVE_MOAB

      call ocn_import_moab(Eclock, ierr)
      if (ierr /= 0) then
         call mpas_log_write('Error in ocn_import_moab', MPAS_LOG_CRIT)
      endif
#endif


      ! Ensures MPAS AM write/compute startup steps are performed
      call ocn_analysis_compute_startup(domain_ptr, ierr)

      ! Handle writing initial state
      call mpas_pool_get_config(domain_ptr % configs, 'config_write_output_on_startup', config_write_output_on_startup)
      if (config_write_output_on_startup) then
         call mpas_stream_mgr_write(domain_ptr % streamManager, 'output', forceWriteNow=.true., ierr=ierr)

         ! Reset config to false, so we don't write the state every coupling interval.
         config_write_output_on_startup = .false.
      end if

      ! Initialize time average fields
      block_ptr => domain_ptr % blocklist
      do while(associated(block_ptr))
         call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
         call ocn_time_average_coupled_init(forcingPool)
         block_ptr => block_ptr % next
      end do

      ! During integration, time level 1 stores the model state at the beginning of the
      !   time step, and time level 2 stores the state advanced dt in time by timestep(...)
      ! This integration loop continues for a single coupling interval.
      do while (.not. mpas_is_alarm_ringing(domain_ptr % clock, coupleAlarmID, ierr=ierr))
         call mpas_stream_mgr_read(domain_ptr % streamManager, ierr=ierr)
         call mpas_stream_mgr_reset_alarms(domain_ptr % streamManager, direction=MPAS_STREAM_INPUT, ierr=ierr)

         ! Build forcing arrays.
         if (debugOn) call mpas_log_write('   Building forcing arrays')
         block_ptr => domain_ptr % blocklist
         do while(associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
            call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
            call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)

            call ocn_forcing_build_fraction_absorbed_array(meshPool, statePool, forcingPool, ierr, 1)
            call ocn_frazil_forcing_build_arrays(domain, meshPool, forcingPool, statePool, ierr)
            call mpas_timer_start("land_ice_build_arrays", .false.)
            call ocn_surface_land_ice_fluxes_build_arrays(meshPool, &
                                                          forcingPool, scratchPool, statePool, ierr)
            call mpas_timer_stop("land_ice_build_arrays")
            call ocn_surface_land_ice_fluxes_accumulate_fluxes(meshPool, forcingPool, &
                                                               statePool, dt, ierr)

            call ocn_eddy_compute_mixed_layer_depth(statePool, forcingPool)
            if (config_use_GM .or. config_submesoscale_enable) then
                call ocn_eddy_compute_buoyancy_gradient()
            end if

            if (config_submesoscale_enable) then
                call mpas_timer_start("submesoscale eddy velocity compute", .false.)
                call ocn_submesoscale_compute_velocity()
                call mpas_timer_stop("submesoscale eddy velocity compute")
            end if
            ! Compute normalGMBolusVelocity, relativeSlope and RediDiffVertCoef if respective flags are turned on
            if (config_use_Redi.or.config_use_GM) then
              call ocn_gm_compute_Bolus_velocity(statePool, meshPool, scratchPool, timeLevelIn=1)
            endif

            block_ptr => block_ptr % next
         end do
         if (debugOn) call mpas_log_write('   Finished building forcing arrays')

         ! Read forcing streams for next timestep
         if (iam==0.and.debugOn) then
            call mpas_log_write('   Reading forcing streams')
         endif

         ! read next time level data required for variable shortwave
         call mpas_timer_start('io_shortwave',.false.)
         call ocn_get_shortWaveData(domain % streamManager, domain, domain_ptr % clock, .false.)
         call mpas_timer_stop('io_shortwave')

         ! read next time level data required for ecosys forcing
         call mpas_pool_get_config(domain % configs, 'config_use_ecosysTracers', config_use_ecosysTracers)
         if (config_use_ecosysTracers) then
           call mpas_timer_start('io_ecosys',.false.)
           call ocn_get_ecosysData(domain % streamManager, domain, domain_ptr % clock, .false.)
           call mpas_timer_stop('io_ecosys')
         endif

         ! read next time level data required for CFC forcing
         call mpas_pool_get_config(domain % configs, 'config_use_CFCTracers', config_use_CFCTracers)
         if (config_use_CFCTracers) then
           call mpas_timer_start('io_CFC',.false.)
           call ocn_get_CFCData(domain % streamManager, domain, domain_ptr % clock, .false.)
           call mpas_timer_stop('io_CFC')
         endif

         ! read next time level data required for monthly surface salinity restoring
         call mpas_pool_get_config(domain % configs, 'config_use_activeTracers_surface_restoring',  &
                                                      config_use_activeTracers_surface_restoring)
         call mpas_pool_get_config(domain % configs, 'config_use_surface_salinity_monthly_restoring',  &
                                                      config_use_surface_salinity_monthly_restoring)
         if (config_use_activeTracers_surface_restoring .and.  config_use_surface_salinity_monthly_restoring) then
             if ( mpas_is_alarm_ringing(domain % clock, 'salinityDataReadAlarm', ierr=ierr) ) then
                 call mpas_reset_clock_alarm(domain % clock, 'salinityDataReadAlarm', ierr=ierr)
                 call mpas_timer_start('io_monthly_surface_salinity',.false.)
                 call ocn_get_surfaceSalinityData(domain % streamManager, domain, domain_ptr % clock, .false.)
                 call mpas_timer_stop('io_monthly_surface_salinity')
             endif
         endif
         if (iam==0.and.debugOn) then
            call mpas_log_write('   Finished reading forcing streams')
         endif

         itimestep = itimestep + 1
         call mpas_advance_clock(domain_ptr % clock)

         currTime = mpas_get_clock_time(domain_ptr % clock, MPAS_NOW, ierr)
         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
         ! write time stamp at every step
         call mpas_dmpar_get_time(current_wallclock_time)
         write (WCstring,'(F18.3)') current_wallclock_time
         call mpas_log_write(trim(timeStamp) // '  WC time:' // WCstring)

        ! pre-timestep analysis computation
        if (debugOn) call mpas_log_write('   Starting analysis precompute', masterOnly=.true.)
        call ocn_analysis_precompute(domain, ierr)
        if (debugOn) call mpas_log_write('   Finished analysis precompute', masterOnly=.true.)

         if (debugOn) call mpas_log_write('   Performing forward update')
         call mpas_timer_start("time integration", .false.)
         call ocn_timestep(domain, dt, timeStamp)
         call mpas_timer_stop("time integration")
         if (debugOn) call mpas_log_write('   Finished performing forward update')

         ! Move time level 2 fields back into time level 1 for next time step
         block_ptr => domain_ptr % blocklist
         do while(associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
            call mpas_pool_shift_time_levels(statePool)
            block_ptr => block_ptr % next
         end do

        ! update scaled data ice-shelf melt fluxes based on remove ice runoff
        call ocn_update_scaled_dismf(domain)

         if (debugOn) call mpas_log_write('   Computing analysis members')
         call ocn_analysis_compute(domain_ptr, ierr)
         if (iam==0.and.debugOn) then
            call mpas_log_write( '   Finished computing analysis members')
            call mpas_log_write('   Preparing analysis members for restart')
         endif
         call ocn_analysis_restart(domain_ptr, ierr)
         if (iam==0.and.debugOn) then
            call mpas_log_write('   Finished preparing analysis members for restart')
            call mpas_log_write('   Writing analysis member output')
         endif
         call ocn_analysis_write(domain_ptr, ierr)
         if (debugOn) call mpas_log_write('   Finished writing analysis member output')

         ! Reset any restart alarms to prevent restart files being written without the coupler requesting it.
         if (debugOn) call mpas_log_write('   Resetting restart alarms')
         call mpas_stream_mgr_begin_iteration(domain % streamManager)

         do while ( mpas_stream_mgr_get_next_stream( domain % streamManager, streamID=streamName, &
                    directionProperty=streamDirection, activeProperty=streamActive ) )

            if ( streamActive .and. streamDirection == MPAS_STREAM_INPUT_OUTPUT ) then
               call mpas_stream_mgr_reset_alarms(domain % streamManager, streamID=streamName, ierr=ierr)
            end if
         end do
         if (debugOn) call mpas_log_write('   Finished resetting restart alarms')

         if (debugOn) call mpas_log_write('   Writing output streams')
         call mpas_stream_mgr_write(domain % streamManager, streamID='output', ierr=ierr)
         call mpas_stream_mgr_reset_alarms(domain % streamManager, streamID='output', ierr=ierr)

         call mpas_stream_mgr_write(domain % streamManager, ierr=ierr)
         call mpas_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)
         if (iam==0.and.debugOn) then
            call mpas_log_write('   Finished writing output streams')
            call mpas_log_write('   Validating ocean state')
         endif

         call ocn_validate_state(domain, timeLevel=1)
         if (debugOn) call mpas_log_write('   Completed validating ocean state')

         if (debugOn) call mpas_log_write('Completed timestep '//trim(timeStamp))
      end do

      ! update coupled variables that are calculated on coupling intervals
      block_ptr => domain_ptr % blocklist
      do while(associated(block_ptr))
         call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
         call mpas_pool_get_dimension(forcingPool, 'index_avgSSHGradientZonal', index_avgZonalSSHGradient)
         call mpas_pool_get_dimension(forcingPool, 'index_avgSSHGradientMeridional', index_avgMeridionalSSHGradient)
         call mpas_pool_get_array(forcingPool, 'avgSSHGradient', avgSSHGradient)
         block_ptr => block_ptr % next
      end do

      ! Check if coupler wants us to write a restart file.
      ! We only write restart files at the end of a coupling interval
      if (seq_timemgr_RestartAlarmIsOn(EClock)) then
         if(trim(config_sw_absorption_type)=='ohlmann00') call ocn_shortwave_forcing_write_restart(domain_ptr)

         if (config_use_ecosysTracers) call ocn_ecosys_forcing_write_restart(domain_ptr)

         if (config_use_CFCTracers) call ocn_CFC_forcing_write_restart(domain_ptr)

         if (config_use_activeTracers_surface_restoring .and. config_use_surface_salinity_monthly_restoring)  &
             call ocn_salinity_restoring_forcing_write_restart(domain_ptr)
         ! Write a restart file, because the coupler asked for it.
         if (debugOn) call mpas_log_write('Writing restart streams')
         call mpas_stream_mgr_begin_iteration(domain % streamManager)

         do while ( mpas_stream_mgr_get_next_stream( domain % streamManager, streamID=streamName, &
                    directionProperty=streamDirection, activeProperty=streamActive ) )

            if ( streamActive .and. streamDirection == MPAS_STREAM_INPUT_OUTPUT ) then
               if (debugOn) call mpas_log_write('   Writing stream ' // trim(streamName))
               call mpas_stream_mgr_write(domain % streamManager, forceWriteNow=.true., streamID=streamName, ierr=ierr)
               if (debugOn) call mpas_log_write('   Finished writing stream ' // trim(streamName))
            end if
         end do

         if ( iam == 0 ) then
            open(22, file=config_restart_timestamp_name, form='formatted', status='replace')
            write(22, *) trim(timeStamp)
            close(22)
         end if

         if (debugOn) call mpas_log_write('Finished writing restart streams')
      end if

      ! Export state to coupler
      if (debugOn) call mpas_log_write('Exporting ocean state')
      call ocn_export_mct(o2x_o, ierr)
      if (debugOn) call mpas_log_write('Finished exporting ocean state')
#ifdef HAVE_MOAB
      call ocn_export_moab(EClock)
#endif
      call check_clocks_sync(domain % clock, Eclock, ierr)

      call mpas_pool_get_config(domain % configs, 'config_scale_dismf_by_removed_ice_runoff', &
                                config_scale_dismf_by_removed_ice_runoff)
      if (config_scale_dismf_by_removed_ice_runoff) then
         ! independent of space so should be no need to loop over blocks
         block_ptr => domain % blocklist
         call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
         call MPAS_pool_get_array(forcingPool, "runningMeanRemovedIceRunoff", &
                                  runningMeanRemovedIceRunoff)
         call seq_infodata_PutData(infodata, rmean_rmv_ice_runoff=runningMeanRemovedIceRunoff)
      end if

      ! Reset I/O logs
      call shr_file_setLogUnit (shrlogunit)
      call shr_file_setLogLevel(shrloglev)

      call mpas_log_write('=== Completed coupling interval in ocn_run_mct. ===', flushNow=.true.)


!-----------------------------------------------------------------------
!EOC

  end subroutine ocn_run_mct!}}}

!***********************************************************************
!BOP
!
! !IROUTINE: ocn_final_mct
!
! !INTERFACE:
  subroutine ocn_final_mct( EClock, cdata_o, x2o_o, o2x_o)!{{{

!
! !DESCRIPTION:
! Finalize MPASO
!
! !USES:
! !ARGUMENTS:
    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata_o
    type(mct_aVect)             , intent(inout) :: x2o_o
    type(mct_aVect)             , intent(inout) :: o2x_o
!
! !REVISION HISTORY:
! Author: Doug Jacobsen
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------


    integer :: shrloglev, shrlogunit
    integer :: ioerr
    character(len=StrKIND) :: ioerrstr

    type (domain_type), pointer :: domain_ptr

    integer :: ierr, iam

!-----------------------------------------------------------------------
    iam = domain % dminfo % my_proc_id

    domain_ptr => domain

    ! Set MPAS Log module instance
    mpas_log_info => domain % logInfo

    ! Setup I/O logs
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (ocnLogUnit)

    ! Finalize MPASO
    iErr = domain_ptr % core % core_finalize(domain_ptr)
    if ( iErr /= 0 ) then
       call mpas_log_write('Core finalize failed for core ' // trim(domain % core % coreName), MPAS_LOG_CRIT)
    end if

    call mpas_timer_write()

    call MPAS_stream_mgr_finalize(domain % streamManager)

    call mpas_log_finalize(iErr)
    if ( iErr /= 0 ) then
       write(ocnLogUnit,*) 'ERROR: log finalize failed for core ' // trim(domain % core % coreName)
       call mpas_dmpar_abort(domain % dminfo)
    end if

    call mpas_framework_finalize(domain % dminfo, domain, io_system)
#ifdef HAVE_MOAB
    ! deallocate moab fields array
     deallocate (o2x_om)
#endif
    ! Reset I/O logs
    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

  end subroutine ocn_final_mct!}}}

!***********************************************************************
!BOP
!IROUTINE: ocn_SetGSMap_mct
! !INTERFACE:

 subroutine ocn_SetGSMap_mct( mpicom_ocn, OCNID, gsMap_ocn )!{{{

! !DESCRIPTION:
!  This routine sets up the MPASO grid numbering for MCT
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

    implicit none
    integer        , intent(in)    :: mpicom_ocn
    integer        , intent(in)    :: OCNID
    type(mct_gsMap), intent(inout) :: gsMap_ocn

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    integer,allocatable :: gindex(:)

    integer :: i, n, lsize, gsize, ier

    type (block_type), pointer :: block_ptr
    type (mpas_pool_type), pointer :: meshPool

    integer, dimension(:), pointer :: indexToCellID

    integer, pointer :: nCellsSolve

    ! Loop over all cells in all blocks to determine total number.
    n = 0
    block_ptr => domain % blocklist
    do while(associated(block_ptr))
      call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

      n = n + nCellsSolve
      block_ptr => block_ptr % next
    end do

    ! Determine total number of cells across all processors
    lsize = n
    call mpas_dmpar_sum_int(domain % dminfo, lsize, gsize)
    allocate(gindex(lsize),stat=ier)

    ! Setup the mapping (gindex)
    n = 0
    block_ptr => domain % blocklist
    do while(associated(block_ptr))
      call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

      call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)

      do i = 1, nCellsSolve
        n = n + 1
        gindex(n) = indexToCellID(i)
      end do
      block_ptr => block_ptr % next
    end do

    ! Init the gsMap with gindex
    call mct_gsMap_init( gsMap_ocn, gindex, mpicom_ocn, OCNID, lsize, gsize )

    deallocate(gindex)

!-----------------------------------------------------------------------
!EOC

  end subroutine ocn_SetGSMap_mct!}}}

!***********************************************************************
!BOP
! !IROUTINE: ocn_domain_mct
! !INTERFACE:

 subroutine ocn_domain_mct( lsize, gsMap_o, dom_o )!{{{

! !DESCRIPTION:
!  This routine sets up the MCT domain for MPASO
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT/OUTPUT PARAMETERS:

    integer        , intent(in)    :: lsize
    type(mct_gsMap), intent(in)    :: gsMap_o
    type(mct_ggrid), intent(inout) :: dom_o

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
    integer, pointer :: idata(:)
    real(kind=RKIND), pointer :: data(:)
    real(kind=RKIND) :: r2d

    integer :: i, n, ier

    type (block_type), pointer :: block_ptr

    type (mpas_pool_type), pointer :: meshPool, forcingPool

    integer, pointer :: nCellsSolve

    real (kind=RKIND), dimension(:), pointer :: lonCell, latCell, areaCell

    integer, dimension(:), pointer :: landIceMask

    real (kind=RKIND), pointer :: sphere_radius

    r2d = 180.0/pii

!-------------------------------------------------------------------
!
!  initialize mct domain type
!
!-------------------------------------------------------------------

    call mct_gGrid_init( GGrid=dom_o, CoordChars=trim(seq_flds_dom_coord), &
       OtherChars=trim(seq_flds_dom_other), lsize=lsize )
    call mct_aVect_zero(dom_o%data)
    allocate(data(lsize))

!-------------------------------------------------------------------
!
! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
!
!-------------------------------------------------------------------

    call mct_gsMap_orderedPoints(gsMap_o, my_task, idata)
    call mct_gGrid_importIAttr(dom_o,'GlobGridNum',idata,lsize)

!-------------------------------------------------------------------
!
! Determine domain
! Initialize attribute vector with special value
!
!-------------------------------------------------------------------

    data(:) = -9999.0_R8
    call mct_gGrid_importRAttr(dom_o,"lat"  ,data,lsize)
    call mct_gGrid_importRAttr(dom_o,"lon"  ,data,lsize)
    call mct_gGrid_importRAttr(dom_o,"area" ,data,lsize)
    call mct_gGrid_importRAttr(dom_o,"aream",data,lsize)
    data(:) = 1.0_R8
    call mct_gGrid_importRAttr(dom_o,"mask",data,lsize)
    call mct_gGrid_importRAttr(dom_o,"frac",data,lsize)

!-------------------------------------------------------------------
!
! Fill in correct values for domain components
!
!-------------------------------------------------------------------

    n = 0
    block_ptr => domain % blocklist
    do while(associated(block_ptr))
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

       call mpas_pool_get_array(meshPool, 'lonCell', lonCell)

       do i = 1, nCellsSolve
          n = n + 1
          data(n) = lonCell(i) * r2d
       end do

       block_ptr => block_ptr % next
    end do
    call mct_gGrid_importRattr(dom_o,"lon",data,lsize)

    n = 0
    block_ptr => domain % blocklist
    do while(associated(block_ptr))
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

       call mpas_pool_get_array(meshPool, 'latCell', latCell)

       do i = 1, nCellsSolve
          n = n + 1
          data(n) = latCell(i) * r2d
       end do
       block_ptr => block_ptr % next
    end do
    call mct_gGrid_importRattr(dom_o,"lat",data,lsize)

    n = 0
    block_ptr => domain % blocklist
    do while(associated(block_ptr))
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

       call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

       call mpas_pool_get_config(meshPool, 'sphere_radius', sphere_radius)
       do i = 1, nCellsSolve
          n = n + 1
          data(n) = areaCell(i) / (sphere_radius * sphere_radius)
       end do
       block_ptr => block_ptr % next
    end do
    call mct_gGrid_importRattr(dom_o,"area",data,lsize)

    ! Build mask and frac based on landIceMask
    block_ptr => domain % blocklist
    do while (associated(block_ptr))
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)

       call mpas_pool_get_array(forcingPool, 'landIceMask', landIceMask)

       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

       if ( associated(landIceMask) ) then
          do i = 1, nCellsSolve
             if ( landIceMask(i) == 1 ) then
                 data(i) = 0.0_RKIND
             else
                 data(i) = 1.0_RKIND
             end if
          end do
       else
          do i = 1, nCellsSolve
             data(i) = 1.0_RKIND
          end do
       end if

       block_ptr => block_ptr % next
    end do

    call mct_gGrid_importRattr(dom_o,"mask",data,lsize)
    call mct_gGrid_importRattr(dom_o,"frac",data,lsize)

    deallocate(data)
    deallocate(idata)

!-----------------------------------------------------------------------
!EOC

  end subroutine ocn_domain_mct!}}}

#ifdef HAVE_MOAB
!***********************************************************************
!BOP
! !IROUTINE: ocn_domain_moab
! !INTERFACE:

 subroutine ocn_domain_moab( mbid )!{{{

  use iso_c_binding, only : C_NULL_CHAR
! !DESCRIPTION:
!  This routine sets up the MOAB domain-related info for MPASO
!
! !REVISION HISTORY:
!  2025-01-21: initial version is copy from ocn_domain_mct
!
! !INPUT/OUTPUT PARAMETERS:

    integer        , intent(in)    :: mbid

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
    real(kind=RKIND), pointer :: data(:)
    real(kind=RKIND) :: r2d

    integer :: i, n, ier, ierr
    character(CXX)           :: tagname

    type (block_type), pointer :: block_ptr

    type (mpas_pool_type), pointer :: meshPool, forcingPool

    integer, pointer :: nCellsSolve

    real (kind=RKIND), dimension(:), pointer :: lonCell, latCell, areaCell

    integer, dimension(:), pointer :: landIceMask

    real (kind=RKIND), pointer :: sphere_radius

    integer nvert(3), nvise(3), nbl(3), nsurf(3), nvisBC(3)
    integer arrsize, ent_type

    r2d = 180.0/pii

!-------------------------------------------------------------------
!
! Determine local grid size
!
!-------------------------------------------------------------------
    ierr  = iMOAB_GetMeshInfo ( mbid, nvert, nvise, nbl, nsurf, nvisBC )

    arrsize = nvise(1)

    allocate(data(arrsize))


!-------------------------------------------------------------------
!
! Fill in correct values for domain tags
!
!-------------------------------------------------------------------

    n = 0
    block_ptr => domain % blocklist
    do while(associated(block_ptr))
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

       call mpas_pool_get_array(meshPool, 'lonCell', lonCell)

       do i = 1, nCellsSolve
          n = n + 1
          data(n) = lonCell(i) * r2d
       end do

       block_ptr => block_ptr % next
    end do
    tagname = "lon"//C_NULL_CHAR
    ent_type=1
    ierr = iMOAB_SetDoubleTagStorage ( mbid, tagname, arrSize , ent_type, data)

    n = 0
    block_ptr => domain % blocklist
    do while(associated(block_ptr))
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

       call mpas_pool_get_array(meshPool, 'latCell', latCell)

       do i = 1, nCellsSolve
          n = n + 1
          data(n) = latCell(i) * r2d
       end do
       block_ptr => block_ptr % next
    end do
    tagname = "lat"//C_NULL_CHAR
    ent_type=1
    ierr = iMOAB_SetDoubleTagStorage ( mbid, tagname, arrSize , ent_type, data)

    n = 0
    block_ptr => domain % blocklist
    do while(associated(block_ptr))
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

       call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

       call mpas_pool_get_config(meshPool, 'sphere_radius', sphere_radius)
       do i = 1, nCellsSolve
          n = n + 1
          data(n) = areaCell(i) / (sphere_radius * sphere_radius)
       end do
       block_ptr => block_ptr % next
    end do
    tagname = "area"//C_NULL_CHAR
    ent_type=1
    ierr = iMOAB_SetDoubleTagStorage ( mbid, tagname, arrSize , ent_type, data)

    ! Build mask and frac based on landIceMask
    block_ptr => domain % blocklist
    do while (associated(block_ptr))
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)

       call mpas_pool_get_array(forcingPool, 'landIceMask', landIceMask)

       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

       if ( associated(landIceMask) ) then
          do i = 1, nCellsSolve
             if ( landIceMask(i) == 1 ) then
                 data(i) = 0.0_RKIND
             else
                 data(i) = 1.0_RKIND
             end if
          end do
       else
          do i = 1, nCellsSolve
             data(i) = 1.0_RKIND
          end do
       end if

       block_ptr => block_ptr % next
    end do
    tagname = "mask"//C_NULL_CHAR
    ent_type=1
    ierr = iMOAB_SetDoubleTagStorage ( mbid, tagname, arrSize , ent_type, data)
    tagname = "frac"//C_NULL_CHAR
    ent_type=1
    ierr = iMOAB_SetDoubleTagStorage ( mbid, tagname, arrSize , ent_type, data)

    deallocate(data)

!-----------------------------------------------------------------------
!EOC
  end subroutine ocn_domain_moab !}}}
#endif


!***********************************************************************
!BOP
! !IROUTINE: ocn_import_mct
! !INTERFACE:

 subroutine ocn_import_mct(x2o_o, errorCode)!{{{

! !DESCRIPTION:
!-----------------------------------------------------------------------
!  This routine receives message from cpl7 driver
!
!    The following fields are always received from the coupler:
!
!    o  taux   -- zonal wind stress (taux)                 (W/m2   )
!    o  tauy   -- meridonal wind stress (tauy)             (W/m2   )
!    o  snow   -- water flux due to snow                   (kg/m2/s)
!    o  rain   -- water flux due to rain                   (kg/m2/s)
!    o  evap   -- evaporation flux                         (kg/m2/s)
!    o  meltw  -- snow melt flux                           (kg/m2/s)
!    o  salt   -- salt                                     (kg(salt)/m2/s)
!    o  swnet  -- net short-wave heat flux                 (W/m2   )
!    o  sen    -- sensible heat flux                       (W/m2   )
!    o  lwup   -- longwave radiation (up)                  (W/m2   )
!    o  lwdn   -- longwave radiation (down)                (W/m2   )
!    o  melth  -- heat flux from snow&ice melt             (W/m2   )
!    o  ifrac  -- ice fraction                             (%)
!    o  rofl   -- river runoff flux                        (kg/m2/s)
!    o  rofi   -- ice runoff flux                          (kg/m2/s)
!    o  rofDIN -- DIN runoff flux                          (kg/m2/s)
!    o  rofDIP -- DIP runoff flux                          (kg/m2/s)
!    o  rofDON -- DON runoff flux                          (kg/m2/s)
!    o  rofDOP -- DOP runoff flux                          (kg/m2/s)
!    o  rofDOC -- DOC runoff flux                          (kg/m2/s)
!    o  rofPP  -- PP  runoff flux                          (kg/m2/s)
!    o  rofDSi -- DSi runoff flux                          (kg/m2/s)
!    o  rofPOC -- POC runoff flux                          (kg/m2/s)
!    o  rofPN  -- PN  runoff flux                          (kg/m2/s)
!    o  rofDIC -- DIC runoff flux                          (kg/m2/s)
!    o  rofFe  -- Fe  runoff flux                          (kg/m2/s)
!
!    The following fields are sometimes received from the coupler,
!      depending on model options:
!
!    o  pslv   -- atmospheric pressure at sea level        (Pa)
!    o  duu10n -- 10m wind speed squared                   (m^2/s^2)
!    o  co2prog-- bottom atm level prognostic co2
!    o  co2diag-- bottom atm level diagnostic co2
!
!-----------------------------------------------------------------------
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

    type(mct_aVect)   , intent(inout) :: x2o_o

! !OUTPUT PARAMETERS:

   integer, intent(out) :: &
      errorCode              ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   character (len=StrKIND) ::   &
      label,                 &
      message

   integer ::  &
      i,n

   integer, pointer :: indexIT, indexIS

   logical, pointer :: config_use_ecosysTracers, &
                       config_use_ecosysTracers_sea_ice_coupling, &
                       config_use_DMSTracers,    &
                       config_use_DMSTracers_sea_ice_coupling,  &
                       config_use_MacroMoleculesTracers,  &
                       config_use_MacroMoleculesTracers_sea_ice_coupling, &
                       config_use_CFCTracers,  &
                       config_remove_ais_river_runoff, &
                       config_remove_ais_ice_runoff, &
                       config_cvmix_kpp_use_theory_wave

   character(len=StrKIND), pointer :: config_ecosys_atm_co2_option, &
                                      config_ecosys_atm_alt_co2_option

   real (kind=RKIND), pointer :: config_ecosys_atm_co2_constant_value

   real (kind=RKIND), pointer :: config_density0

   type (block_type), pointer :: block_ptr

   type (mpas_pool_type), pointer :: meshPool,             &
                                     forcingPool,          &
                                     ecosysAuxiliary,      &
                                     ecosysSeaIceCoupling, &
                                     DMSSeaIceCoupling,    &
                                     MacroMoleculesSeaIceCoupling,      &
                                     CFCAuxiliary

   integer, pointer :: nCellsSolve

   type (field1DReal), pointer :: windStressZonalField, windStressMeridionalField, &
                                  latentHeatFluxField, sensibleHeatFluxField, &
                                  longWaveHeatFluxUpField, &
                                  longWaveHeatFluxDownField, &
                                  evaporationFluxField, seaIceHeatFluxField, icebergHeatFluxField, &
                                  snowFluxField, seaIceFreshWaterFluxField, icebergFreshWaterFluxField, &
                                  seaIceSalinityFluxField, &
                                  riverRunoffFluxField, iceRunoffFluxField, &
                                  removedRiverRunoffFluxField, removedIceRunoffFluxField, &
                                  shortWaveHeatFluxField, rainFluxField, &
                                  atmosphericPressureField, iceFractionField, &
                                  seaIcePressureField, windSpeedSquared10mField, &
                                  atmosphericCO2Field, atmosphericCO2_ALT_CO2Field,  &
                                  iceFluxDICField, &
                                  iceFluxDONField, &
                                  iceFluxNO3Field, &
                                  iceFluxSiO3Field, &
                                  iceFluxNH4Field, &
                                  iceFluxDMSField, &
                                  iceFluxDMSPField, &
                                  iceFluxDOCrField, &
                                  iceFluxFeParticulateField, &
                                  iceFluxFeDissolvedField, &
                                  iceFluxDustField, &
                                  riverFluxNO3Field,  &
                                  riverFluxPO4Field,  &
                                  riverFluxSiO3Field, &
                                  riverFluxDOCField,  &
                                  riverFluxDONField,  &
                                  riverFluxDOPField,  &
                                  riverFluxDICField,  &
                                  riverFluxALKField,  &
                                  riverFluxFeField,   &
                                  landIceFreshwaterFluxField, &
                                  landIceHeatFluxField, &
                                  landIceFractionField, &
                                  windSpeed10mField, &
                                 significantWaveHeightField, &
                                  peakWaveFrequencyField, &
                                  peakWaveDirectionField
                                  !landIcePressureField

   type (field2DReal), pointer :: iceFluxPhytoCField, &
                                  iceFluxDOCField

   type (field2DReal), pointer :: landIceInterfaceTracersField

   type (field2DReal), pointer :: stokesDriftZonalWavenumberField, &
                                  stokesDriftMeridionalWavenumberField

   real (kind=RKIND), dimension(:), pointer :: windStressZonal, windStressMeridional, &
                                               latentHeatFlux, sensibleHeatFlux, &
                                               longWaveHeatFluxUp, &
                                               longWaveHeatFluxDown, &
                                               evaporationFlux, seaIceHeatFlux, icebergHeatFlux, &
                                               snowFlux, seaIceFreshWaterFlux, icebergFreshWaterFlux, &
                                               seaIceSalinityFlux, &
                                               riverRunoffFlux, iceRunoffFlux, &
                                               removedRiverRunoffFlux, removedIceRunoffFlux, &
                                               shortWaveHeatFlux, rainFlux, &
                                               atmosphericPressure, iceFraction, &
                                               seaIcePressure, windSpeedSquared10m, &
                                               atmosphericCO2, atmosphericCO2_ALT_CO2,  &
                                               windSpeedSquared10mCFC, &
                                               iceFluxDIC,       &
                                               iceFluxDON, &
                                               iceFluxNO3, &
                                               iceFluxSiO3, &
                                               iceFluxNH4, &
                                               iceFluxDMS, &
                                               iceFluxDMSP, &
                                               iceFluxDOCr, &
                                               iceFluxFeParticulate, &
                                               iceFluxFeDissolved, &
                                               iceFluxDust, &
                                               riverFluxNO3,  &
                                               riverFluxPO4,  &
                                               riverFluxSiO3, &
                                               riverFluxDOC,  &
                                               riverFluxDON,  &
                                               riverFluxDOP,  &
                                               riverFluxDIC,  &
                                               riverFluxALK,  &
                                               riverFluxFe,   &
                                               landIceFreshwaterFlux, &
                                               landIceHeatFlux, &
                                               landIceFraction, &
                                               areaCell, &
                                               windSpeed10m, &
                                               significantWaveHeight, &
                                               peakWaveFrequency, &
                                               peakWaveDirection
                                               !landIcePressure

   real (kind=RKIND), dimension(:), pointer :: latCell

   real (kind=RKIND), dimension(:,:), pointer :: iceFluxPhytoC, &
                                                 iceFluxDOC, &
                                                 stokesDriftZonalWavenumber, &
                                                 stokesDriftMeridionalWavenumber

   real (kind=RKIND) :: removedRiverRunoffFluxThisProc, removedIceRunoffFluxThisProc
   real (kind=RKIND) :: removedRiverRunoffFluxReduced, removedIceRunoffFluxReduced

   real (kind=RKIND), dimension(:,:), pointer :: landIceInterfaceTracers

   real (kind=RKIND) :: riverFactor

!-----------------------------------------------------------------------
!
!  zero out padded cells
!
!-----------------------------------------------------------------------

   errorCode = 0

!-----------------------------------------------------------------------
!
!  unpack and distribute wind stress, then convert to correct units
!  and rotate components to local coordinates
!
!-----------------------------------------------------------------------

   ! get configure options
   call mpas_pool_get_config(domain % configs, 'config_use_ecosysTracers', config_use_ecosysTracers)
   call mpas_pool_get_config(domain % configs, 'config_use_DMSTracers', config_use_DMSTracers)
   call mpas_pool_get_config(domain % configs, 'config_use_MacroMoleculesTracers', config_use_MacroMoleculesTracers)
   call mpas_pool_get_config(domain % configs, 'config_density0', config_density0)
   call mpas_pool_get_config(domain % configs, 'config_use_ecosysTracers_sea_ice_coupling',  &
                                                config_use_ecosysTracers_sea_ice_coupling)
   call mpas_pool_get_config(domain % configs, 'config_use_DMSTracers_sea_ice_coupling',  &
                                                config_use_DMSTracers_sea_ice_coupling)
   call mpas_pool_get_config(domain % configs, 'config_use_MacroMoleculesTracers_sea_ice_coupling',  &
                                                config_use_MacroMoleculesTracers_sea_ice_coupling)
   call mpas_pool_get_config(domain % configs, 'config_use_CFCTracers', config_use_CFCTracers)
   call mpas_pool_get_config(domain % configs, 'config_remove_ais_river_runoff', config_remove_ais_river_runoff)
   call mpas_pool_get_config(domain % configs, 'config_remove_ais_ice_runoff', config_remove_ais_ice_runoff)
   call mpas_pool_get_config(domain % configs, 'config_cvmix_kpp_use_theory_wave', config_cvmix_kpp_use_theory_wave)

   n = 0
   removedRiverRunoffFluxThisProc = 0.0_RKIND
   removedIceRunoffFluxThisProc = 0.0_RKIND

   block_ptr => domain % blocklist
   do while(associated(block_ptr))
      call mpas_pool_get_subpool(block_ptr % structs, 'mesh',    meshPool)
      call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)

      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

      call mpas_pool_get_field(forcingPool, 'windStressZonal', windStressZonalField)
      call mpas_pool_get_field(forcingPool, 'windStressMeridional', windStressMeridionalField)
      call mpas_pool_get_field(forcingPool, 'latentHeatFlux', latentHeatFluxField)
      call mpas_pool_get_field(forcingPool, 'sensibleHeatFlux', sensibleHeatFluxField)
      call mpas_pool_get_field(forcingPool, 'longWaveHeatFluxUp', longWaveHeatFluxUpField)
      call mpas_pool_get_field(forcingPool, 'longWaveHeatFluxDown', longWaveHeatFluxDownField)
      call mpas_pool_get_field(forcingPool, 'evaporationFlux', evaporationFluxField)
      call mpas_pool_get_field(forcingPool, 'seaIceHeatFlux', seaIceHeatFluxField)
      call mpas_pool_get_field(forcingPool, 'icebergHeatFlux', icebergHeatFluxField)
      call mpas_pool_get_field(forcingPool, 'snowFlux', snowFluxField)
      call mpas_pool_get_field(forcingPool, 'seaIceFreshWaterFlux', seaIceFreshWaterFluxField)
      call mpas_pool_get_field(forcingPool, 'icebergFreshWaterFlux', icebergFreshWaterFluxField)
      call mpas_pool_get_field(forcingPool, 'seaIceSalinityFlux', seaIceSalinityFluxField)
      call mpas_pool_get_field(forcingPool, 'riverRunoffFlux', riverRunoffFluxField)
      call mpas_pool_get_field(forcingPool, 'iceRunoffFlux', iceRunoffFluxField)
      call mpas_pool_get_field(forcingPool, 'shortWaveHeatFlux', shortWaveHeatFluxField)
      call mpas_pool_get_field(forcingPool, 'rainFlux', rainFluxField)
      call mpas_pool_get_field(forcingPool, 'atmosphericPressure', atmosphericPressureField)
      call mpas_pool_get_field(forcingPool, 'seaIcePressure', seaIcePressureField)
      call mpas_pool_get_field(forcingPool, 'iceFraction', iceFractionField)
      call mpas_pool_get_field(forcingPool, 'iceRunoffFlux', iceRunoffFluxField)
      call mpas_pool_get_field(forcingPool, 'removedRiverRunoffFlux', removedRiverRunoffFluxField)
      call mpas_pool_get_field(forcingPool, 'removedIceRunoffFlux', removedIceRunoffFluxField)
      call mpas_pool_get_field(forcingPool, 'stokesDriftZonalWavenumber', stokesDriftZonalWavenumberField)
      call mpas_pool_get_field(forcingPool, 'stokesDriftMeridionalWavenumber', stokesDriftMeridionalWavenumberField)
      call mpas_pool_get_field(forcingPool, 'significantWaveHeight', significantWaveHeightField)
      call mpas_pool_get_field(forcingPool, 'peakWaveFrequency', peakWaveFrequencyField)
      call mpas_pool_get_field(forcingPool, 'peakWaveDirection', peakWaveDirectionField)

      call mpas_pool_get_field(forcingPool, 'landIceFreshwaterFlux', landIceFreshwaterFluxField)
      call mpas_pool_get_field(forcingPool, 'landIceHeatFlux', landIceHeatFluxField)
      call mpas_pool_get_field(forcingPool, 'landIceFraction', landIceFractionField)
      call mpas_pool_get_field(forcingPool, 'landIceInterfaceTracers', landIceInterfaceTracersField)

      call mpas_pool_get_field(forcingPool, 'windSpeed10m', windSpeed10mField)

      call mpas_pool_get_dimension(forcingPool, 'index_landIceInterfaceTemperature', indexIT)
      call mpas_pool_get_dimension(forcingPool, 'index_landIceInterfaceSalinity', indexIS)

   !call mpas_pool_get_field(forcingPool, 'landIcePressure', landIcePressureField)

      windStressZonal => windStressZonalField % array
      windStressMeridional => windStressMeridionalField % array
      latentHeatFlux => latentHeatFluxField % array
      sensibleHeatFlux => sensibleHeatFluxField % array
      longWaveHeatFluxUp => longWaveHeatFluxUpField % array
      longWaveHeatFluxDown => longWaveHeatFluxDownField % array
      evaporationFlux => evaporationFluxField % array
      seaIceHeatFlux => seaIceHeatFluxField % array
      icebergHeatFlux => icebergHeatFluxField % array
      snowFlux => snowFluxField % array
      seaIceFreshWaterFlux => seaIceFreshWaterFluxField % array
      icebergFreshWaterFlux => icebergFreshWaterFluxField % array
      seaIceSalinityFlux => seaIceSalinityFluxField % array
      riverRunoffFlux => riverRunoffFluxField % array
      iceRunoffFlux => iceRunoffFluxField % array
      shortWaveHeatFlux => shortWaveHeatFluxField % array
      rainFlux => rainFluxField % array
      atmosphericPressure => atmosphericPressureField % array
      seaIcePressure => seaIcePressureField % array
      iceFraction => iceFractionField % array
      iceRunoffFlux => iceRunoffFluxField % array
      removedRiverRunoffFlux => removedRiverRunoffFluxField % array
      removedIceRunoffFlux => removedIceRunoffFluxField % array
      landIceFreshwaterFlux => landIceFreshwaterFluxField % array
      landIceHeatFlux => landIceHeatFluxField % array
      landIceInterfaceTracers => landIceInterfaceTracersField % array
      landIceFraction => landIceFractionField % array
      windSpeed10m => windSpeed10mField % array
      stokesDriftZonalWavenumber => stokesDriftZonalWavenumberField % array
      stokesDriftMeridionalWavenumber => stokesDriftMeridionalWavenumberField % array
      significantWaveHeight => significantWaveHeightField % array
      peakWaveFrequency => peakWaveFrequencyField % array
      peakWaveDirection => peakWaveDirectionField % array
      !landIcePressure => landIcePressureField % array

      call mpas_pool_get_array(meshPool, 'latCell', latCell)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

      ! BGC fields
      if (config_use_ecosysTracers) then
         call mpas_pool_get_subpool(forcingPool, 'ecosysAuxiliary', ecosysAuxiliary)

         call mpas_pool_get_field(ecosysAuxiliary, 'windSpeedSquared10m', windSpeedSquared10mField)
         windSpeedSquared10m => windSpeedSquared10mField % array
         call mpas_pool_get_field(ecosysAuxiliary, 'atmosphericCO2', atmosphericCO2Field)
         atmosphericCO2 => atmosphericCO2Field % array
         call mpas_pool_get_field(ecosysAuxiliary, 'atmosphericCO2_ALT_CO2', atmosphericCO2_ALT_CO2Field)
         atmosphericCO2_ALT_CO2 => atmosphericCO2_ALT_CO2Field % array

         if (config_use_ecosysTracers_river_inputs_from_coupler) then
           call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxNO3' , riverFluxNO3Field)
           call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxPO4' , riverFluxPO4Field)
           call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxDON' , riverFluxDONField)
           call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxDOP' , riverFluxDOPField)
           call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxSiO3', riverFluxSiO3Field)
           call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxDOC' , riverFluxDOCField)
           call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxDIC' , riverFluxDICField)
           call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxALK' , riverFluxALKField)
           call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxFe'  , riverFluxFeField)
           riverFluxNO3  => riverFluxNO3Field  % array
           riverFluxPO4  => riverFluxPO4Field  % array
           riverFluxDON  => riverFluxDONField  % array
           riverFluxDOP  => riverFluxDOPField  % array
           riverFluxSiO3 => riverFluxSiO3Field % array
           riverFluxDOC  => riverFluxDOCField  % array
           riverFluxDIC  => riverFluxDICField  % array
           riverFluxALK  => riverFluxALKField  % array
           riverFluxFe   => riverFluxFeField   % array
         endif

         call mpas_pool_get_config(domain % configs, 'config_ecosys_atm_co2_option',   &
            config_ecosys_atm_co2_option)
         call mpas_pool_get_config(domain % configs, 'config_ecosys_atm_alt_co2_option',   &
            config_ecosys_atm_alt_co2_option)
         call mpas_pool_get_config(domain % configs, 'config_ecosys_atm_co2_constant_value',   &
            config_ecosys_atm_co2_constant_value)

         if (config_use_ecosysTracers_sea_ice_coupling) then
            call mpas_pool_get_subpool(forcingPool, 'ecosysSeaIceCoupling', ecosysSeaIceCoupling)

            call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxPhytoC', iceFluxPhytoCField)
            call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxDIC', iceFluxDICField)
            call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxNO3', iceFluxNO3Field)
            call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxSiO3', iceFluxSiO3Field)
            call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxNH4', iceFluxNH4Field)
            call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxDOCr', iceFluxDOCrField)
            call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxFeParticulate', iceFluxFeParticulateField)
            call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxFeDissolved', iceFluxFeDissolvedField)
            call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxDust', iceFluxDustField)
            call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxDOC', iceFluxDOCField)
            call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxDON', iceFluxDONField)

            iceFluxPhytoC => iceFluxPhytoCField % array
            iceFluxDIC => iceFluxDICField % array
            iceFluxNO3 => iceFluxNO3Field % array
            iceFluxSiO3 => iceFluxSiO3Field % array
            iceFluxNH4 => iceFluxNH4Field % array
            iceFluxDOCr => iceFluxDOCrField % array
            iceFluxFeParticulate => iceFluxFeParticulateField % array
            iceFluxFeDissolved => iceFluxFeDissolvedField % array
            iceFluxDust => iceFluxDustField % array
            iceFluxDOC => iceFluxDOCField % array
            iceFluxDON => iceFluxDONField % array
         endif
      endif
      if (config_use_DMSTracers .and. config_use_DMSTracers_sea_ice_coupling) then
         call mpas_pool_get_subpool(forcingPool, 'DMSSeaIceCoupling', DMSSeaIceCoupling)

         call mpas_pool_get_field(DMSSeaIceCoupling, 'iceFluxDMS', iceFluxDMSField)
         call mpas_pool_get_field(DMSSeaIceCoupling, 'iceFluxDMSP', iceFluxDMSPField)

         iceFluxDMS => iceFluxDMSField % array
         iceFluxDMSP => iceFluxDMSPField % array
      endif

      ! CFC fields
      if (config_use_CFCTracers) then
         call mpas_pool_get_subpool(forcingPool, 'CFCAuxiliary', CFCAuxiliary)
         call mpas_pool_get_field(CFCAuxiliary, 'windSpeedSquared10mCFC', windSpeedSquared10mField)
         windSpeedSquared10mCFC => windSpeedSquared10mField % array
      endif

      if (config_remove_ais_river_runoff) then
         ! Initialize this field
         removedRiverRunoffFlux(:) = 0.0_RKIND
      endif

      if (config_remove_ais_ice_runoff) then
         ! Initialize this field
         removedIceRunoffFlux(:) = 0.0_RKIND
      endif

      if (config_cvmix_kpp_use_theory_wave) then
         ! Initialize this field
         windSpeed10m(:) = 0.0_RKIND
      endif

      do i = 1, nCellsSolve
        n = n + 1
        if ( windStressZonalField % isActive ) then
           windStressZonal(i) = x2o_o % rAttr(index_x2o_Foxx_taux, n)
        end if
        if ( windStressMeridionalField % isActive ) then
           windStressMeridional(i) = x2o_o % rAttr(index_x2o_Foxx_tauy, n)
        end if

        if ( latentHeatFluxField % isActive ) then
           latentHeatFlux(i) = x2o_o % rAttr(index_x2o_Foxx_lat, n)
        end if
        if ( sensibleHeatFluxField % isActive ) then
           sensibleHeatFlux(i) = x2o_o % rAttr(index_x2o_Foxx_sen, n)
        end if
        if ( longWaveHeatFluxUpField % isActive ) then
           longWaveHeatFluxUp(i) = x2o_o % rAttr(index_x2o_Foxx_lwup, n)
        end if
        if ( longWaveHeatFluxDownField % isActive ) then
           longWaveHeatFluxDown(i) = x2o_o % rAttr(index_x2o_Faxa_lwdn, n)
        end if
        if ( evaporationFluxField % isActive ) then
           evaporationFlux(i) = x2o_o % rAttr(index_x2o_Foxx_evap, n)
        end if
        if ( seaIceHeatFluxField % isActive ) then
           seaIceHeatFlux(i) = x2o_o % rAttr(index_x2o_Fioi_melth, n)
        end if
        if ( icebergHeatFluxField % isActive ) then
           icebergHeatFlux(i) = x2o_o % rAttr(index_x2o_Fioi_bergh, n)
        end if
        if ( snowFluxField % isActive ) then
           snowFlux(i) = x2o_o % rAttr(index_x2o_Faxa_snow, n)
        end if
        if ( seaIceFreshWaterFluxField % isActive ) then
           seaIceFreshWaterFlux(i) = x2o_o % rAttr(index_x2o_Fioi_meltw, n)
        end if
        if ( icebergFreshWaterFluxField % isActive ) then
           icebergFreshWaterFlux(i) = x2o_o % rAttr(index_x2o_Fioi_bergw, n)
        end if
        if ( seaIceSalinityFluxField % isActive ) then
           seaIceSalinityFlux(i) = x2o_o % rAttr(index_x2o_Fioi_salt, n)
        end if
        if ( riverRunoffFluxField % isActive ) then
           riverRunoffFlux(i) = x2o_o % rAttr(index_x2o_Foxx_rofl, n)
           if (config_remove_ais_river_runoff) then
              if (latCell(i) < -1.04719666667_RKIND) then ! 60S in radians
                 removedRiverRunoffFlux(i) = riverRunoffFlux(i)
                 riverRunoffFlux(i) = 0.0_RKIND
                 removedRiverRunoffFluxThisProc = removedRiverRunoffFluxThisProc + removedRiverRunoffFlux(i)
               endif
           endif
        end if
        if ( iceRunoffFluxField % isActive ) then
           iceRunoffFlux(i) = x2o_o % rAttr(index_x2o_Foxx_rofi, n)
           if(iceRunoffFlux(n) < 0.0_RKIND) then
               call shr_sys_abort ('Error: incoming rofi_F is negative')
           end if
           if (config_remove_ais_ice_runoff) then
              if (latCell(i) < -0.99483767345_RKIND) then ! 57S in radians
                 removedIceRunoffFlux(i) = iceRunoffFlux(i)
                 iceRunoffFlux(i) = 0.0_RKIND
                 removedIceRunoffFluxThisProc = removedIceRunoffFluxThisProc + removedIceRunoffFlux(i)
              endif
           endif
        end if
        if ( shortWaveHeatFluxField % isActive ) then
           shortWaveHeatFlux(i) = max(x2o_o % rAttr(index_x2o_Foxx_swnet, n), 0.0_RKIND)
        end if

        if ( rainFluxField % isActive ) then
           rainFlux(i) = x2o_o % rAttr(index_x2o_Faxa_rain, n)
        end if
        if ( atmosphericPressureField % isActive ) then
           atmosphericPressure(i) = x2o_o % rAttr(index_x2o_Sa_pslv,   n)
        end if
        if ( seaIcePressureField % isActive ) then
           ! Set seaIcePressure to be limited to 5m of pressure
           seaIcePressure(i) =  min( x2o_o % rAttr(index_x2o_Si_bpress, n), config_density0 * gravity * 5.0_RKIND )
        end if
        if ( iceFractionField % isActive ) then
           iceFraction(i) = x2o_o % rAttr(index_x2o_Si_ifrac, n)
        end if

        if ( stokesDriftZonalWavenumberField % isActive ) then
           stokesDriftZonalWavenumber(1,i) = x2o_o % rAttr(index_x2o_Sw_ustokes_wavenumber_1, n)
           stokesDriftZonalWavenumber(2,i) = x2o_o % rAttr(index_x2o_Sw_ustokes_wavenumber_2, n)
           stokesDriftZonalWavenumber(3,i) = x2o_o % rAttr(index_x2o_Sw_ustokes_wavenumber_3, n)
           stokesDriftZonalWavenumber(4,i) = x2o_o % rAttr(index_x2o_Sw_ustokes_wavenumber_4, n)
           stokesDriftZonalWavenumber(5,i) = x2o_o % rAttr(index_x2o_Sw_ustokes_wavenumber_5, n)
           stokesDriftZonalWavenumber(6,i) = x2o_o % rAttr(index_x2o_Sw_ustokes_wavenumber_6, n)
        end if
        if ( stokesDriftMeridionalWavenumberField % isActive ) then
           stokesDriftMeridionalWavenumber(1,i) = x2o_o % rAttr(index_x2o_Sw_vstokes_wavenumber_1, n)
           stokesDriftMeridionalWavenumber(2,i) = x2o_o % rAttr(index_x2o_Sw_vstokes_wavenumber_2, n)
           stokesDriftMeridionalWavenumber(3,i) = x2o_o % rAttr(index_x2o_Sw_vstokes_wavenumber_3, n)
           stokesDriftMeridionalWavenumber(4,i) = x2o_o % rAttr(index_x2o_Sw_vstokes_wavenumber_4, n)
           stokesDriftMeridionalWavenumber(5,i) = x2o_o % rAttr(index_x2o_Sw_vstokes_wavenumber_5, n)
           stokesDriftMeridionalWavenumber(6,i) = x2o_o % rAttr(index_x2o_Sw_vstokes_wavenumber_6, n)
        end if
        if ( significantWaveHeightField % isActive ) then
           significantWaveHeight(i) = x2o_o % rAttr(index_x2o_Sw_Hs, n)
        end if
        if ( peakWaveFrequencyField % isActive ) then
           peakWaveFrequency(i) = x2o_o % rAttr(index_x2o_Sw_Fp, n)
        end if
        if ( peakWaveDirectionField % isActive ) then
           peakWaveDirection(i) = x2o_o % rAttr(index_x2o_Sw_Dp, n)
        end if

        if (config_cvmix_kpp_use_theory_wave) then
           if ( windSpeed10mField% isActive ) then
              windSpeed10m(i) = sqrt( x2o_o % rAttr(index_x2o_So_duu10n, n))
           end if
        endif

        if ( landIceFreshwaterFluxField % isActive ) then
           !landIceFreshwaterFlux(i) = x2o_o % rAttr(index_x2o_Fogx_qicelo, n)
        end if
        if ( landIceHeatFluxField % isActive ) then
           !landIceHeatFlux(i) = x2o_o % rAttr(index_x2o_Fogx_qiceho, n)
        end if
        if ( landIceInterfaceTracersField % isActive ) then
           !landIceInterfaceTracers(indexIT, i) = x2o_o % rAttr(index_x2o_Sg_blit, n)
           !landIceInterfaceTracers(indexIS, i) = x2o_o % rAttr(index_x2o_Sg_blis, n)
        end if
        if ( landIceFractionField % isActive ) then
           !landIceFraction(i) = x2o_o % rAttr(index_x2o_Sg_icemask, n)
        end if
        !if ( landIcePressureField % isActive ) then
           !landIcePressure(i) = x2o_o % rAttr(index_x2o_Sg_lithop, n)
        !end if

        ! BGC fields
        if (config_use_ecosysTracers) then
           if ( windSpeedSquared10mField % isActive ) then
              windSpeedSquared10m(i) = x2o_o % rAttr(index_x2o_So_duu10n, n)
           end if
           if ( atmosphericCO2Field % isActive ) then
              if ( config_ecosys_atm_co2_option == 'prognostic') then
                 atmosphericCO2(i) = x2o_o % rAttr(index_x2o_Sa_co2prog, n)
              else if ( config_ecosys_atm_co2_option == 'diagnostic') then
                 atmosphericCO2(i) = x2o_o % rAttr(index_x2o_Sa_co2diag, n)
              else if ( config_ecosys_atm_co2_option == 'bcrc') then
                 atmosphericCO2(i) = config_ecosys_atm_co2_constant_value
              else if ( config_ecosys_atm_co2_option == 'bcrd') then
                 atmosphericCO2(i) = config_ecosys_atm_co2_constant_value
              else if ( config_ecosys_atm_co2_option == 'bdrc') then
                 atmosphericCO2(i) = x2o_o % rAttr(index_x2o_Sa_co2diag, n)
              else if ( config_ecosys_atm_co2_option == 'bdrd') then
                 atmosphericCO2(i) = x2o_o % rAttr(index_x2o_Sa_co2diag, n)
              else
                 atmosphericCO2(i) = config_ecosys_atm_co2_constant_value
              end if
           end if
           if ( atmosphericCO2_ALT_CO2Field % isActive ) then
              if ( config_ecosys_atm_alt_co2_option == 'prognostic') then
                 atmosphericCO2_ALT_CO2(i) = x2o_o % rAttr(index_x2o_Sa_co2prog, n)
              else if ( config_ecosys_atm_alt_co2_option == 'diagnostic') then
                 atmosphericCO2_ALT_CO2(i) = x2o_o % rAttr(index_x2o_Sa_co2diag, n)
              else if ( config_ecosys_atm_alt_co2_option == 'bcrc') then
                 atmosphericCO2_ALT_CO2(i) = config_ecosys_atm_co2_constant_value
              else if ( config_ecosys_atm_alt_co2_option == 'bcrd') then
                 atmosphericCO2_ALT_CO2(i) = x2o_o % rAttr(index_x2o_Sa_co2diag, n)
              else if ( config_ecosys_atm_alt_co2_option == 'bdrc') then
                 atmosphericCO2_ALT_CO2(i) = config_ecosys_atm_co2_constant_value
              else if ( config_ecosys_atm_alt_co2_option == 'bdrd') then
                 atmosphericCO2_ALT_CO2(i) = x2o_o % rAttr(index_x2o_Sa_co2diag, n)
              else
                 atmosphericCO2_ALT_CO2(i) = config_ecosys_atm_co2_constant_value
              end if
           end if

           if (config_use_ecosysTracers_river_inputs_from_coupler) then
             riverFluxNO3(i)  = x2o_o % rAttr(index_x2o_Foxx_rofDIN, n)
             riverFluxPO4(i)  = x2o_o % rAttr(index_x2o_Foxx_rofDIP, n)
             riverFluxDON(i)  = x2o_o % rAttr(index_x2o_Foxx_rofDON, n)
             riverFluxDOP(i)  = x2o_o % rAttr(index_x2o_Foxx_rofDOP, n)
             riverFluxSiO3(i) = x2o_o % rAttr(index_x2o_Foxx_rofDSi, n)
             riverFluxDOC(i)  = x2o_o % rAttr(index_x2o_Foxx_rofDOC, n)
             riverFluxDIC(i)  = x2o_o % rAttr(index_x2o_Foxx_rofDIC, n)
             riverFluxFe(i)   = x2o_o % rAttr(index_x2o_Foxx_rofFe , n)

! convert from kgNutrient/(m2-s) to mmol/m3 m/s
             riverFactor = 1.e6_RKIND
             riverFluxNO3(i)  = riverFluxNO3(i)*riverFactor/14.007_RKIND
             riverFluxPO4(i)  = riverFluxPO4(i)*riverFactor/30.974_RKIND
             riverFluxDON(i)  = riverFluxDON(i)*riverFactor/14.007_RKIND
             riverFluxDOP(i)  = riverFluxDOP(i)*riverFactor/30.974_RKIND
             riverFluxSiO3(i) = riverFluxSiO3(i)*riverFactor/28.085_RKIND
             riverFluxDOC(i)  = riverFluxDOC(i)*riverFactor/12.001_RKIND
             riverFluxDIC(i)  = riverFluxDIC(i)*riverFactor/12.001_RKIND
             riverFluxFe(i)   = riverFluxFe(i)*riverFactor/55.845_RKIND

             riverFluxALK(i)  = riverFluxDIC(i)
           endif

           if (config_use_ecosysTracers_sea_ice_coupling) then
              if ( iceFluxPhytoCField % isActive ) then
                 iceFluxPhytoC(1,i) = x2o_o % rAttr(index_x2o_Fioi_algae1, n)
                 iceFluxPhytoC(2,i) = x2o_o % rAttr(index_x2o_Fioi_algae2, n)
                 iceFluxPhytoC(3,i) = x2o_o % rAttr(index_x2o_Fioi_algae3, n)
              endif
              if ( iceFluxDICField % isActive ) then
                 iceFluxDIC(i) = x2o_o % rAttr(index_x2o_Fioi_dic1, n)
              endif
              if ( iceFluxNO3Field % isActive ) then
                 iceFluxNO3(i) = x2o_o % rAttr(index_x2o_Fioi_no3, n)
              endif
              if ( iceFluxSiO3Field % isActive ) then
                 iceFluxSiO3(i) = x2o_o % rAttr(index_x2o_Fioi_sio3, n)
              endif
              if ( iceFluxNH4Field % isActive ) then
                 iceFluxNH4(i) = x2o_o % rAttr(index_x2o_Fioi_nh4, n)
              endif
              if ( iceFluxDOCrField % isActive ) then
                 iceFluxDOCr(i) = x2o_o % rAttr(index_x2o_Fioi_docr, n)
              endif
              if ( iceFluxFeParticulateField % isActive ) then
                 iceFluxFeParticulate(i) = x2o_o % rAttr(index_x2o_Fioi_fep1, n)
              endif
              if ( iceFluxFeDissolvedField % isActive ) then
                 iceFluxFeDissolved(i) = x2o_o % rAttr(index_x2o_Fioi_fed1, n)
              endif
              if ( iceFluxDustField % isActive ) then
                 iceFluxDust(i) = x2o_o % rAttr(index_x2o_Fioi_dust1, n)
              endif
              if ( iceFluxDOCField % isActive ) then
                 iceFluxDOC(1,i) = x2o_o % rAttr(index_x2o_Fioi_doc1, n)
                 iceFluxDOC(2,i) = x2o_o % rAttr(index_x2o_Fioi_doc2, n)
                 iceFluxDOC(3,i) = x2o_o % rAttr(index_x2o_Fioi_doc3, n)
              endif
              if ( iceFluxDONField % isActive ) then
                 iceFluxDON(i) = x2o_o % rAttr(index_x2o_Fioi_don1, n)
              endif
            endif
         endif
         if (config_use_DMSTracers .and. config_use_DMSTracers_sea_ice_coupling) then
           if ( iceFluxDMSField % isActive ) then
              iceFluxDMS(i) = x2o_o % rAttr(index_x2o_Fioi_dms, n)
           endif
           if ( iceFluxDMSPField % isActive ) then
              !JW TODO: dmspp? dmspd? the sum?
              iceFluxDMSP(i) = x2o_o % rAttr(index_x2o_Fioi_dmspp, n)
           endif
         endif

        ! CFC fields
        if (config_use_CFCTracers) then
           if ( windSpeedSquared10mField % isActive ) then
              windSpeedSquared10mCFC(i) = x2o_o % rAttr(index_x2o_So_duu10n, n)
           end if
        end if

      end do

      block_ptr => block_ptr % next
   end do

   call mpas_pool_get_subpool(domain % blocklist % structs, 'forcing', forcingPool)

   call mpas_pool_get_field(forcingPool, 'windStressZonal', windStressZonalField)
   call mpas_pool_get_field(forcingPool, 'windStressMeridional', windStressMeridionalField)
   call mpas_pool_get_field(forcingPool, 'latentHeatFlux', latentHeatFluxField)
   call mpas_pool_get_field(forcingPool, 'sensibleHeatFlux', sensibleHeatFluxField)
   call mpas_pool_get_field(forcingPool, 'longWaveHeatFluxUp', longWaveHeatFluxUpField)
   call mpas_pool_get_field(forcingPool, 'longWaveHeatFluxDown', longWaveHeatFluxDownField)
   call mpas_pool_get_field(forcingPool, 'evaporationFlux', evaporationFluxField)
   call mpas_pool_get_field(forcingPool, 'seaIceHeatFlux', seaIceHeatFluxField)
   call mpas_pool_get_field(forcingPool, 'icebergHeatFlux', icebergHeatFluxField)
   call mpas_pool_get_field(forcingPool, 'snowFlux', snowFluxField)
   call mpas_pool_get_field(forcingPool, 'seaIceFreshWaterFlux', seaIceFreshWaterFluxField)
   call mpas_pool_get_field(forcingPool, 'icebergFreshWaterFlux', icebergFreshWaterFluxField)
   call mpas_pool_get_field(forcingPool, 'seaIceSalinityFlux', seaIceSalinityFluxField)
   call mpas_pool_get_field(forcingPool, 'riverRunoffFlux', riverRunoffFluxField)
   call mpas_pool_get_field(forcingPool, 'iceRunoffFlux', iceRunoffFluxField)
   call mpas_pool_get_field(forcingPool, 'shortWaveHeatFlux', shortWaveHeatFluxField)
   call mpas_pool_get_field(forcingPool, 'rainFlux', rainFluxField)
   call mpas_pool_get_field(forcingPool, 'atmosphericPressure', atmosphericPressureField)
   call mpas_pool_get_field(forcingPool, 'seaIcePressure', seaIcePressureField)
   call mpas_pool_get_field(forcingPool, 'iceFraction', iceFractionField)
   call mpas_pool_get_field(forcingPool, 'stokesDriftZonalWavenumber', stokesDriftZonalWavenumberField)
   call mpas_pool_get_field(forcingPool, 'stokesDriftMeridionalWavenumber', stokesDriftMeridionalWavenumberField)
   call mpas_pool_get_field(forcingPool, 'significantWaveHeight', significantWaveHeightField)
   call mpas_pool_get_field(forcingPool, 'peakWaveFrequency', peakWaveFrequencyField)
   call mpas_pool_get_field(forcingPool, 'peakWaveDirection', peakWaveDirectionField)

   call mpas_pool_get_field(forcingPool, 'landIceFreshwaterFlux', landIceFreshwaterFluxField)
   call mpas_pool_get_field(forcingPool, 'landIceHeatFlux', landIceHeatFluxField)
   call mpas_pool_get_field(forcingPool, 'landIceFraction', landIceFractionField)
   call mpas_pool_get_field(forcingPool, 'landIceInterfaceTracers', landIceInterfaceTracersField)

   call mpas_pool_get_field(forcingPool, 'windSpeed10m', windSpeed10mField)

   call mpas_pool_get_dimension(forcingPool, 'index_landIceInterfaceTemperature', indexIT)
   call mpas_pool_get_dimension(forcingPool, 'index_landIceInterfaceSalinity', indexIS)
   !call mpas_pool_get_field(forcingPool, 'landIcePressure', landIcePressureField)

   ! BGC fields
   if (config_use_ecosysTracers) then
      call mpas_pool_get_subpool(forcingPool, 'ecosysAuxiliary', ecosysAuxiliary)

      call mpas_pool_get_field(ecosysAuxiliary, 'windSpeedSquared10m', windSpeedSquared10mField)
      call mpas_pool_get_field(ecosysAuxiliary, 'atmosphericCO2', atmosphericCO2Field)
      call mpas_pool_get_field(ecosysAuxiliary, 'atmosphericCO2_ALT_CO2', atmosphericCO2_ALT_CO2Field)

      if (config_use_ecosysTracers_river_inputs_from_coupler) then
        call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxNO3' , riverFluxNO3Field)
        call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxPO4' , riverFluxPO4Field)
        call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxDON' , riverFluxDONField)
        call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxDOP' , riverFluxDOPField)
        call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxSiO3', riverFluxSiO3Field)
        call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxDOC' , riverFluxDOCField)
        call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxDIC' , riverFluxDICField)
        call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxALK' , riverFluxALKField)
        call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxFe'  , riverFluxFeField)
      endif

      if (config_use_ecosysTracers_sea_ice_coupling) then
         call mpas_pool_get_subpool(forcingPool, 'ecosysSeaIceCoupling', ecosysSeaIceCoupling)

         call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxPhytoC', iceFluxPhytoCField)
         call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxDIC', iceFluxDICField)
         call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxNO3', iceFluxNO3Field)
         call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxSiO3', iceFluxSiO3Field)
         call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxNH4', iceFluxNH4Field)
         call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxDOCr', iceFluxDOCrField)
         call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxFeParticulate', iceFluxFeParticulateField)
         call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxFeDissolved', iceFluxFeDissolvedField)
         call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxDust', iceFluxDustField)
         call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxDOC', iceFluxDOCField)
         call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxDON', iceFluxDONField)
      endif
   endif
   if (config_use_DMSTracers .and. config_use_DMSTracers_sea_ice_coupling) then
      call mpas_pool_get_subpool(forcingPool, 'DMSSeaIceCoupling', DMSSeaIceCoupling)

      call mpas_pool_get_field(DMSSeaIceCoupling, 'iceFluxDMS', iceFluxDMSField)
      call mpas_pool_get_field(DMSSeaIceCoupling, 'iceFluxDMSP', iceFluxDMSPField)
   endif

   ! CFC fields
   if (config_use_CFCTracers) then
      call mpas_pool_get_subpool(forcingPool, 'CFCAuxiliary', CFCAuxiliary)
      call mpas_pool_get_field(CFCAuxiliary, 'windSpeedSquared10mCFC', windSpeedSquared10mField)
   endif

   if ( windStressMeridionalField % isActive ) then
      call mpas_dmpar_exch_halo_field(windStressMeridionalField)
   end if
   if ( windStressZonalField % isActive ) then
      call mpas_dmpar_exch_halo_field(windStressZonalField)
   end if
   if ( latentHeatFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(latentHeatFluxField)
   end if
   if ( sensibleHeatFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(sensibleHeatFluxField)
   end if
   if ( longWaveHeatFluxUpField % isActive ) then
      call mpas_dmpar_exch_halo_field(longWaveHeatFluxUpField)
   end if
   if ( longWaveHeatFluxDownField % isActive ) then
      call mpas_dmpar_exch_halo_field(longWaveHeatFluxDownField)
   end if
   if ( evaporationFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(evaporationFluxField)
   end if
   if ( seaIceHeatFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(seaIceHeatFluxField)
   end if
   if ( icebergHeatFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(icebergHeatFluxField)
   end if
   if ( snowFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(snowFluxField)
   end if
   if ( seaIceFreshWaterFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(seaIceFreshWaterFluxField)
   end if
   if ( icebergFreshWaterFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(icebergFreshWaterFluxField)
   end if
   if ( seaIceSalinityFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(seaIceSalinityFluxField)
   end if
   if ( riverRunoffFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(riverRunoffFluxField)
   end if
   if ( iceRunoffFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(iceRunoffFluxField)
   end if
   if ( shortWaveHeatFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(shortWaveHeatFluxField)
   end if
   if ( rainFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(rainFluxField)
   end if
   if ( atmosphericPressureField % isActive ) then
      call mpas_dmpar_exch_halo_field(atmosphericPressureField)
   end if
   if ( seaIcePressureField % isActive ) then
      call mpas_dmpar_exch_halo_field(seaIcePressureField)
   end if
   if ( iceFractionField % isActive ) then
      call mpas_dmpar_exch_halo_field(iceFractionField)
   end if
   if ( stokesDriftZonalWavenumberField % isActive ) then
      call mpas_dmpar_exch_halo_field(stokesDriftZonalWavenumberField)
   end if
   if ( stokesDriftMeridionalWavenumberField % isActive ) then
      call mpas_dmpar_exch_halo_field(stokesDriftMeridionalWavenumberField)
   end if
   if ( significantWaveHeightField % isActive ) then
      call mpas_dmpar_exch_halo_field(significantWaveHeightField)
   end if
   if ( peakWaveFrequencyField % isActive ) then
      call mpas_dmpar_exch_halo_field(peakWaveFrequencyField)
   end if
   if ( peakWaveDirectionField % isActive ) then
      call mpas_dmpar_exch_halo_field(peakWaveDirectionField)
   end if

   if ( landIceFreshwaterFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(landIceFreshwaterFluxField)
   end if
   if ( landIceHeatFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(landIceHeatFluxField)
   end if
   if ( landIceInterfaceTracersField % isActive ) then
      call mpas_dmpar_exch_halo_field(landIceInterfaceTracersField)
   end if
   if ( landIceFractionField % isActive ) then
      call mpas_dmpar_exch_halo_field(landIceFractionField)
   end if
!   if ( landIcePressureField % isActive ) then
!      call mpas_dmpar_exch_halo_field(landIcePressureField)
!   end if

   if ( windSpeed10mField % isActive ) then
      call mpas_dmpar_exch_halo_field(windSpeed10mField)
   end if

   ! BGC fields
   if (config_use_ecosysTracers) then

      if ( windSpeedSquared10mField % isActive ) then
         call mpas_dmpar_exch_halo_field(windSpeedSquared10mField)
      end if
      if ( atmosphericCO2Field % isActive ) then
         call mpas_dmpar_exch_halo_field(atmosphericCO2Field)
      end if
      if ( atmosphericCO2_ALT_CO2Field % isActive ) then
         call mpas_dmpar_exch_halo_field(atmosphericCO2_ALT_CO2Field)
      end if

      if (config_use_ecosysTracers_river_inputs_from_coupler) then
        if ( riverFluxNO3Field % isActive ) then
           call mpas_dmpar_exch_halo_field(riverFluxNO3Field)
        end if
        if ( riverFluxPO4Field % isActive ) then
           call mpas_dmpar_exch_halo_field(riverFluxPO4Field)
        end if
        if ( riverFluxDONField % isActive ) then
           call mpas_dmpar_exch_halo_field(riverFluxDONField)
        end if
        if ( riverFluxDOPField % isActive ) then
           call mpas_dmpar_exch_halo_field(riverFluxDOPField)
        end if
        if ( riverFluxSiO3Field % isActive ) then
           call mpas_dmpar_exch_halo_field(riverFluxSiO3Field)
        end if
        if ( riverFluxDOCField % isActive ) then
           call mpas_dmpar_exch_halo_field(riverFluxDOCField)
        end if
        if ( riverFluxDICField % isActive ) then
           call mpas_dmpar_exch_halo_field(riverFluxDICField)
        end if
        if ( riverFluxALKField % isActive ) then
           call mpas_dmpar_exch_halo_field(riverFluxALKField)
        end if
        if ( riverFluxFeField % isActive ) then
           call mpas_dmpar_exch_halo_field(riverFluxFeField)
        end if
      endif

      if (config_use_ecosysTracers_sea_ice_coupling) then
         if ( iceFluxPhytoCField % isActive ) then
            call mpas_dmpar_exch_halo_field(iceFluxPhytoCField)
         endif
         if ( iceFluxDICField % isActive ) then
            call mpas_dmpar_exch_halo_field(iceFluxDICField)
         endif
         if ( iceFluxNO3Field % isActive ) then
            call mpas_dmpar_exch_halo_field(iceFluxNO3Field)
         endif
         if ( iceFluxSiO3Field % isActive ) then
            call mpas_dmpar_exch_halo_field(iceFluxSiO3Field)
         endif
         if ( iceFluxNH4Field % isActive ) then
            call mpas_dmpar_exch_halo_field(iceFluxNH4Field)
         endif
         if ( iceFluxDOCrField % isActive ) then
            call mpas_dmpar_exch_halo_field(iceFluxDOCrField)
         endif
         if ( iceFluxFeParticulateField % isActive ) then
            call mpas_dmpar_exch_halo_field(iceFluxFeParticulateField)
         endif
         if ( iceFluxFeDissolvedField % isActive ) then
            call mpas_dmpar_exch_halo_field(iceFluxFeDissolvedField)
         endif
         if ( iceFluxDustField % isActive ) then
            call mpas_dmpar_exch_halo_field(iceFluxDustField)
         endif
         if ( iceFluxDOCField % isActive ) then
            call mpas_dmpar_exch_halo_field(iceFluxDOCField)
         endif
         if ( iceFluxDONField % isActive ) then
            call mpas_dmpar_exch_halo_field(iceFluxDONField)
         endif
      endif
   endif
   if (config_use_DMSTracers .and. config_use_DMSTracers_sea_ice_coupling) then
      if ( iceFluxDMSField % isActive ) then
         call mpas_dmpar_exch_halo_field(iceFluxDMSField)
      endif
      if ( iceFluxDMSPField % isActive ) then
         call mpas_dmpar_exch_halo_field(iceFluxDMSPField)
      endif
   endif

   ! CFC fields
   if (config_use_CFCTracers .and. .not. config_use_ecosysTracers) then
      if ( windSpeedSquared10mField % isActive ) then
         call mpas_dmpar_exch_halo_field(windSpeedSquared10mField)
      end if
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine ocn_import_mct!}}}
!***********************************************************************
!BOP
! !IROUTINE: ocn_export_mct
! !INTERFACE:

 subroutine ocn_export_mct(o2x_o, errorCode)   !{{{

! !DESCRIPTION:
!  This routine calls the routines necessary to send MPASO fields to
!  the CCSM cpl7 driver
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT/OUTPUT PARAMETERS:

   type(mct_aVect)   , intent(inout) :: o2x_o

! !OUTPUT PARAMETERS:

   integer, intent(out) :: &
      errorCode              ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer :: i, n
   integer, pointer :: nCellsSolve, index_temperatureSurfaceValue, index_salinitySurfaceValue, &
                       index_avgZonalSurfaceVelocity, index_avgMeridionalSurfaceVelocity, &
                       index_avgZonalSSHGradient, index_avgMeridionalSSHGradient

   type (block_type), pointer :: block_ptr

   type (mpas_pool_type), pointer :: meshPool,             &
                                     forcingPool,          &
                                     statePool,            &
                                     tracersPool,          &
                                     ecosysAuxiliary,      &
                                     ecosysSeaIceCoupling, &
                                     DMSSeaIceCoupling,    &
                                     MacroMoleculesSeaIceCoupling

   integer, dimension(:), pointer :: landIceMask

   real (kind=RKIND), dimension(:), pointer :: seaIceEnergy, accumulatedFrazilIceMass, frazilSurfacePressure, &
                                               avgTotalFreshWaterTemperatureFlux, &
                                               avgCO2_gas_flux, DMSFlux, surfaceUpwardCO2Flux, &
                                               avgOceanSurfaceDIC, &
                                               avgOceanSurfaceDON, &
                                               avgOceanSurfaceNO3, &
                                               avgOceanSurfaceSiO3, &
                                               avgOceanSurfaceNH4, &
                                               avgOceanSurfaceDMS, &
                                               avgOceanSurfaceDMSP, &
                                               avgOceanSurfaceDOCr, &
                                               avgOceanSurfaceDOCSemiLabile, &
                                               avgOceanSurfaceFeParticulate, &
                                               avgOceanSurfaceFeDissolved, &
                                               ssh, &
                                               avgLandIceFreshwaterFlux, &
                                               avgRemovedRiverRunoffFlux, &
                                               avgRemovedIceRunoffFlux, &
                                               avgLandIceHeatFlux, &
                                               avgRemovedIceRunoffHeatFlux, &
                                               avgThermalForcingAtCritDepth

   real (kind=RKIND), dimension(:,:), pointer :: avgTracersSurfaceValue, avgSurfaceVelocity, &
                                                 avgSSHGradient, avgOceanSurfacePhytoC, &
                                                 avgOceanSurfaceDOC, layerThickness

   real (kind=RKIND) :: surfaceFreezingTemp

   logical, pointer :: frazilIceActive,          &
                       config_remove_ais_river_runoff, &
                       config_remove_ais_ice_runoff, &
                       config_use_ecosysTracers, &
                       config_use_DMSTracers,    &
                       config_use_MacroMoleculesTracers,  &
                       config_use_ecosysTracers_sea_ice_coupling, &
                       config_use_DMSTracers_sea_ice_coupling,  &
                       config_use_MacroMoleculesTracers_sea_ice_coupling

   character (len=StrKIND), pointer :: config_land_ice_flux_mode
   character (len=StrKIND), pointer :: config_glc_thermal_forcing_coupling_mode

   logical :: keepFrazil


   ! get configure options
   call mpas_pool_get_package(domain % packages, 'frazilIceActive', frazilIceActive)
   call mpas_pool_get_config(domain % configs, 'config_use_ecosysTracers', config_use_ecosysTracers)
   call mpas_pool_get_config(domain % configs, 'config_land_ice_flux_mode', config_land_ice_flux_mode)
   call mpas_pool_get_config(domain % configs, 'config_remove_ais_river_runoff', config_remove_ais_river_runoff)
   call mpas_pool_get_config(domain % configs, 'config_remove_ais_ice_runoff', config_remove_ais_ice_runoff)
   call mpas_pool_get_config(domain % configs, 'config_glc_thermal_forcing_coupling_mode', &
                                                config_glc_thermal_forcing_coupling_mode)
   call mpas_pool_get_config(domain % configs, 'config_use_DMSTracers', config_use_DMSTracers)
   call mpas_pool_get_config(domain % configs, 'config_use_MacroMoleculesTracers', config_use_MacroMoleculesTracers)
   call mpas_pool_get_config(domain % configs, 'config_use_ecosysTracers_sea_ice_coupling',  &
                                                config_use_ecosysTracers_sea_ice_coupling)
   call mpas_pool_get_config(domain % configs, 'config_use_DMSTracers_sea_ice_coupling',  &
                                                config_use_DMSTracers_sea_ice_coupling)
   call mpas_pool_get_config(domain % configs, 'config_use_MacroMoleculesTracers_sea_ice_coupling',  &
                                                config_use_MacroMoleculesTracers_sea_ice_coupling)

   n = 0
   block_ptr => domain % blocklist
   do while(associated(block_ptr))
     call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
     call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
     call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)

     call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

     call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

     call mpas_pool_get_dimension(forcingPool, 'index_avgTemperatureSurfaceValue', index_temperatureSurfaceValue)
     call mpas_pool_get_dimension(forcingPool, 'index_avgSalinitySurfaceValue', index_salinitySurfaceValue)
     call mpas_pool_get_dimension(forcingPool, 'index_avgSurfaceVelocityZonal', index_avgZonalSurfaceVelocity)
     call mpas_pool_get_dimension(forcingPool, 'index_avgSurfaceVelocityMeridional', index_avgMeridionalSurfaceVelocity)
     call mpas_pool_get_dimension(forcingPool, 'index_avgSSHGradientZonal', index_avgZonalSSHGradient)
     call mpas_pool_get_dimension(forcingPool, 'index_avgSSHGradientMeridional', index_avgMeridionalSSHGradient)

     call mpas_pool_get_array(statePool, 'ssh', ssh, 1)
     call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)

     call mpas_pool_get_array(forcingPool, 'landIceMask', landIceMask)
     call mpas_pool_get_array(forcingPool, 'avgTracersSurfaceValue', avgTracersSurfaceValue)
     call mpas_pool_get_array(forcingPool, 'avgSurfaceVelocity', avgSurfaceVelocity)
     call mpas_pool_get_array(forcingPool, 'avgSSHGradient', avgSSHGradient)
     call mpas_pool_get_array(forcingPool, 'avgTotalFreshWaterTemperatureFlux', avgTotalFreshWaterTemperatureFlux)

     if ( frazilIceActive ) then
        call mpas_pool_get_array(forcingPool, 'seaIceEnergy', seaIceEnergy)
        call mpas_pool_get_array(forcingPool, 'frazilSurfacePressure', frazilSurfacePressure)
        call mpas_pool_get_array(statePool, 'accumulatedFrazilIceMass', accumulatedFrazilIceMass, 1)
     end if

     ! Cryo fields
     if (trim(config_land_ice_flux_mode) == 'standalone' .or. trim(config_land_ice_flux_mode) == 'data') then
        call mpas_pool_get_array(forcingPool, 'avgLandIceFreshwaterFlux', avgLandIceFreshwaterFlux)
        call mpas_pool_get_array(forcingPool, 'avgLandIceHeatFlux', avgLandIceHeatFlux)
     endif
     if (config_remove_ais_river_runoff) then
        call mpas_pool_get_array(forcingPool, 'avgRemovedRiverRunoffFlux', avgRemovedRiverRunoffFlux)
     endif
     if (config_remove_ais_ice_runoff) then
        call mpas_pool_get_array(forcingPool, 'avgRemovedIceRunoffFlux', avgRemovedIceRunoffFlux)
        call mpas_pool_get_array(forcingPool, 'avgRemovedIceRunoffHeatFlux', avgRemovedIceRunoffHeatFlux)
     endif
     if (trim(config_glc_thermal_forcing_coupling_mode) == '2d') then
        call mpas_pool_get_array(forcingPool, 'avgThermalForcingAtCritDepth', avgThermalForcingAtCritDepth)
     endif

     ! BGC fields
     if (config_use_ecosysTracers) then

        call mpas_pool_get_subpool(forcingPool, 'ecosysAuxiliary', ecosysAuxiliary)
        call mpas_pool_get_array(ecosysAuxiliary, 'avgCO2_gas_flux', avgCO2_gas_flux)

     end if

     if (config_use_ecosysTracers .and. config_use_ecosysTracers_sea_ice_coupling) then
        call mpas_pool_get_subpool(forcingPool, 'ecosysSeaIceCoupling', ecosysSeaIceCoupling)

        call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfacePhytoC', avgOceanSurfacePhytoC)
        call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceDIC', avgOceanSurfaceDIC)
        call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceNO3', avgOceanSurfaceNO3)
        call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceSiO3', avgOceanSurfaceSiO3)
        call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceNH4', avgOceanSurfaceNH4)
        call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceDOCr', avgOceanSurfaceDOCr)
        call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceDOCSemiLabile', avgOceanSurfaceDOCSemiLabile)
        call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceFeParticulate', avgOceanSurfaceFeParticulate)
        call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceFeDissolved', avgOceanSurfaceFeDissolved)
     endif
     if (config_use_DMSTracers .and. config_use_DMSTracers_sea_ice_coupling) then
        call mpas_pool_get_subpool(forcingPool, 'DMSSeaIceCoupling', DMSSeaIceCoupling)

        call mpas_pool_get_array(DMSSeaIceCoupling, 'avgOceanSurfaceDMS', avgOceanSurfaceDMS)
        call mpas_pool_get_array(DMSSeaIceCoupling, 'avgOceanSurfaceDMSP', avgOceanSurfaceDMSP)
     endif
     if (config_use_MacroMoleculesTracers .and. config_use_MacroMoleculesTracers_sea_ice_coupling) then
        call mpas_pool_get_subpool(forcingPool, 'MacroMoleculesSeaIceCoupling', MacroMoleculesSeaIceCoupling)

        call mpas_pool_get_array(MacroMoleculesSeaIceCoupling, 'avgOceanSurfaceDOC', avgOceanSurfaceDOC)
        call mpas_pool_get_array(MacroMoleculesSeaIceCoupling, 'avgOceanSurfaceDON', avgOceanSurfaceDON)
     endif
!    call mpas_pool_get_array(forcingPool, 'CO2Flux', CO2Flux)
!    call mpas_pool_get_array(forcingPool, 'DMSFlux', DMSFlux)
!    call mpas_pool_get_array(forcingPool, 'surfaceUpwardCO2Flux', surfaceUpwardCO2Flux)

     do i = 1, nCellsSolve
       n = n + 1

       o2x_o % rAttr(index_o2x_So_t, n) = avgTracersSurfaceValue(index_temperatureSurfaceValue, i)
       o2x_o % rAttr(index_o2x_So_s, n) = avgTracersSurfaceValue(index_salinitySurfaceValue, i)
       o2x_o % rAttr(index_o2x_So_u, n) = avgSurfaceVelocity(index_avgZonalSurfaceVelocity, i)
       o2x_o % rAttr(index_o2x_So_v, n) = avgSurfaceVelocity(index_avgMeridionalSurfaceVelocity, i)

       o2x_o % rAttr(index_o2x_So_ssh, n)  = ssh(i)
       o2x_o % rAttr(index_o2x_So_dhdx, n) = avgSSHGradient(index_avgZonalSSHGradient, i)
       o2x_o % rAttr(index_o2x_So_dhdy, n) = avgSSHGradient(index_avgMeridionalSSHGradient, i)

       o2x_o % rAttr(index_o2x_Faoo_h2otemp, n) = avgTotalFreshWaterTemperatureFlux(i) * rho_sw * cp_sw

       ! Cryo fields
       if (trim(config_land_ice_flux_mode) == 'standalone' .or. trim(config_land_ice_flux_mode) == 'data') then
          o2x_o % rAttr(index_o2x_Foxo_ismw, n)  = avgLandIceFreshwaterFlux(i)
          o2x_o % rAttr(index_o2x_Foxo_ismh, n)  = avgLandIceHeatFlux(i)
       endif
       if (config_remove_ais_river_runoff) then
          o2x_o % rAttr(index_o2x_Foxo_rrofl, n) = avgRemovedRiverRunoffFlux(i)
       endif
       if (config_remove_ais_ice_runoff) then
          o2x_o % rAttr(index_o2x_Foxo_rrofi, n) = avgRemovedIceRunoffFlux(i)
          o2x_o % rAttr(index_o2x_Foxo_rrofih, n) = avgRemovedIceRunoffHeatFlux(i)
       endif

       if ( frazilIceActive ) then
          ! negative when frazil ice can be melted
          keepFrazil = .true.
          if ( associated(landIceMask) ) then
             if ( landIceMask(i) == 1 ) then
                keepFrazil = .false.
             end if
          end if

          if ( keepFrazil ) then

             ! Calculate energy associated with frazil mass transfer to sea ice if frazil has accumulated
             if ( accumulatedFrazilIceMass(i) > 0.0_RKIND ) then

              seaIceEnergy(i) = accumulatedFrazilIceMass(i) * config_frazil_heat_of_fusion

             ! Otherwise calculate the melt potential where avgTracersSurfaceValue represents only the
             ! top layer of the ocean
             else

              surfaceFreezingTemp = ocn_freezing_temperature(salinity=avgTracersSurfaceValue(index_salinitySurfaceValue, i), &
                 pressure=0.0_RKIND,  inLandIceCavity=.false.)

              seaIceEnergy(i) = min(rho_sw*cp_sw*layerThickness(1, i)*( surfaceFreezingTemp + T0_Kelvin &
                              - avgTracersSurfaceValue(index_temperatureSurfaceValue, i) ), 0.0_RKIND )

             end if

             o2x_o % rAttr(index_o2x_Fioo_q, n)  = seaIceEnergy(i) / ocn_cpl_dt
             o2x_o % rAttr(index_o2x_Fioo_frazil, n) = accumulatedFrazilIceMass(i) / ocn_cpl_dt

          else

             o2x_o % rAttr(index_o2x_Fioo_q, n)  = 0.0_RKIND
             o2x_o % rAttr(index_o2x_Fioo_frazil, n) = 0.0_RKIND
             if (trim(config_land_ice_flux_mode) == 'standalone' .or. trim(config_land_ice_flux_mode) == 'data') then
                o2x_o % rAttr(index_o2x_Foxo_q_li, n) = accumulatedFrazilIceMass(i) * config_frazil_heat_of_fusion / ocn_cpl_dt
                o2x_o % rAttr(index_o2x_Foxo_frazil_li, n) = accumulatedFrazilIceMass(i) / ocn_cpl_dt
             endif

          end if

          ! Reset SeaIce Energy and Accumulated Frazil Ice
          seaIceEnergy(i) = 0.0_RKIND
          accumulatedFrazilIceMass(i) = 0.0_RKIND
          frazilSurfacePressure(i) = 0.0_RKIND
       end if

       ! BGC fields
       if (config_use_ecosysTracers .and. index_o2x_Faoo_fco2_ocn /= 0) then
          ! convert from mmolC/m2/s to kg CO2/m2/s
          o2x_o % rAttr(index_o2x_Faoo_fco2_ocn, n) = avgCO2_gas_flux(i)*44.e-6_RKIND
       endif
       if (config_use_ecosysTracers .and. config_use_ecosysTracers_sea_ice_coupling) then
          o2x_o % rAttr(index_o2x_So_algae1, n) = max(0.0_RKIND,avgOceanSurfacePhytoC(1,i))
          o2x_o % rAttr(index_o2x_So_algae2, n) = max(0.0_RKIND,avgOceanSurfacePhytoC(2,i))
          o2x_o % rAttr(index_o2x_So_algae3, n) = max(0.0_RKIND,avgOceanSurfacePhytoC(3,i))
          o2x_o % rAttr(index_o2x_So_dic1,   n) = max(0.0_RKIND,avgOceanSurfaceDIC(i))
          o2x_o % rAttr(index_o2x_So_doc1,   n) = max(0.0_RKIND,avgOceanSurfaceDOCSemiLabile(i))
          o2x_o % rAttr(index_o2x_So_doc2,   n) = max(0.0_RKIND,avgOceanSurfaceDOCSemiLabile(i))
          o2x_o % rAttr(index_o2x_So_doc3,   n) = max(0.0_RKIND,avgOceanSurfaceDOCSemiLabile(i))
          o2x_o % rAttr(index_o2x_So_don1,   n) = 0.0_RKIND
          o2x_o % rAttr(index_o2x_So_no3,    n) = max(0.0_RKIND,avgOceanSurfaceNO3(i))
          o2x_o % rAttr(index_o2x_So_sio3,   n) = max(0.0_RKIND,avgOceanSurfaceSiO3(i))
          o2x_o % rAttr(index_o2x_So_nh4,    n) = max(0.0_RKIND,avgOceanSurfaceNH4(i))
          o2x_o % rAttr(index_o2x_So_docr,   n) = max(0.0_RKIND,avgOceanSurfaceDOCr(i))
          o2x_o % rAttr(index_o2x_So_fep1,   n) = max(0.0_RKIND,avgOceanSurfaceFeParticulate(i))
          o2x_o % rAttr(index_o2x_So_fed1,   n) = max(0.0_RKIND,avgOceanSurfaceFeDissolved(i))
       endif
       if (config_use_DMSTracers .and. config_use_DMSTracers_sea_ice_coupling) then
          o2x_o % rAttr(index_o2x_So_dms,  n) = max(0.0_RKIND,avgOceanSurfaceDMS(i))
          o2x_o % rAttr(index_o2x_So_dmsp, n) = max(0.0_RKIND,avgOceanSurfaceDMSP(i))
       endif
       if (config_use_MacroMoleculesTracers .and. config_use_MacroMoleculesTracers_sea_ice_coupling) then
          o2x_o % rAttr(index_o2x_So_doc1, n) = max(0.0_RKIND,avgOceanSurfaceDOC(1,i))
          o2x_o % rAttr(index_o2x_So_doc2, n) = max(0.0_RKIND,avgOceanSurfaceDOC(2,i))
          o2x_o % rAttr(index_o2x_So_doc3, n) = max(0.0_RKIND,avgOceanSurfaceDOC(3,i))
          o2x_o % rAttr(index_o2x_So_don1, n) = 0.0_RKIND
       endif
!      o2x_o % rAttr(index_o2x_Faoo_fco2_ocn, n) = CO2Flux(i)
!      o2x_o % rAttr(index_o2x_Faoo_fdms_ocn, n) = DMSFlux(i)
!      o2x_o % rAttr(index_o2x_Faoo_fco2_ocn, n) = surfaceUpwardCO2Flux(i)

!JW       o2x_o % rAttr(index_o2x_So_blt, n) = landIceBoundaryLayerTemperature(i)
!JW       o2x_o % rAttr(index_o2x_So_bls, n) = landIceBoundaryLayerSalinity(i)
!JW       o2x_o % rAttr(index_o2x_So_htv, n) = landIceHeatTransferVelocity(i)
!JW       o2x_o % rAttr(index_o2x_So_stv, n) = landIceSaltTransferVelocity(i)

       if ( trim(config_land_ice_flux_mode) .eq. 'standalone' .or. &
            trim(config_land_ice_flux_mode) .eq. 'coupled'  ) then
          o2x_o % rAttr(index_o2x_So_blt, n) = landIceBoundaryLayerTracers(indexBLT,i)
          o2x_o % rAttr(index_o2x_So_bls, n) = landIceBoundaryLayerTracers(indexBLS,i)
          o2x_o % rAttr(index_o2x_So_htv, n) = landIceTracerTransferVelocities(indexHeatTrans,i)
          o2x_o % rAttr(index_o2x_So_stv, n) = landIceTracerTransferVelocities(indexSaltTrans,i)
          o2x_o % rAttr(index_o2x_So_rhoeff, n) = 0.0_RKIND
       endif
       if (trim(config_glc_thermal_forcing_coupling_mode) == '2d' .and. ocn_c2_glctf) then
          o2x_o % rAttr(index_o2x_So_tf2d, n) = avgThermalForcingAtCritDepth(i)
       endif


       !Fyke: test
       !write(stderrUnit,*) 'n=',n
       !write(stderrUnit,*) 'o2x_o % rAttr(index_o2x_So_blt, n)=',o2x_o % rAttr(index_o2x_So_blt, n)
       !write(stderrUnit,*) 'o2x_o % rAttr(index_o2x_So_bls, n)=',o2x_o % rAttr(index_o2x_So_bls, n)
       !write(stderrUnit,*) 'o2x_o % rAttr(index_o2x_So_htv, n)=',o2x_o % rAttr(index_o2x_So_htv, n)
       !write(stderrUnit,*) 'o2x_o % rAttr(index_o2x_So_stv, n)=',o2x_o % rAttr(index_o2x_So_stv, n)
       !write(stderrUnit,*) 'o2x_o % rAttr(index_o2x_So_rhoeff, n)=',o2x_o % rAttr(index_o2x_So_rhoeff, n)
       !o2x_o % rAttr(index_o2x_So_blt, n) = 0._r8
       !o2x_o % rAttr(index_o2x_So_bls, n) = 34.5_r8
       !o2x_o % rAttr(index_o2x_So_htv, n) = 1.e-4_r8
       !o2x_o % rAttr(index_o2x_So_stv, n) = 3.e-6_r8
       !o2x_o % rAttr(index_o2x_So_rhoeff, n) = 1000._r8*9.81_r8*918._r8 !lithostatic pressure of 1km of ice

     end do

     block_ptr => block_ptr % next
   end do

   errorCode = 0

!-----------------------------------------------------------------------
!EOC

 end subroutine ocn_export_mct!}}}

 subroutine check_clocks_sync(MPASclock, Eclock, err)
    type (MPAS_Clock_type), intent(in) :: MPASclock
    type(ESMF_Clock), intent(in)       :: EClock
    integer, intent(out)               :: err
    ! local variables
    type (MPAS_Time_Type) :: currTime
    integer :: iyear, imonth, iday, ihour, iminute, isecond
    integer :: ymd, tod
    integer :: curr_ymd, curr_tod
    integer :: err_tmp

    err = 0

    ! Check if clocks are in sync
    currTime = mpas_get_clock_time(MPASclock, MPAS_NOW, err_tmp)
    err = ior(err,err_tmp)
    call mpas_get_time(curr_time=currTime, YYYY=iyear, MM=imonth, DD=iday, H=ihour, M=iminute, S=isecond, ierr=err_tmp)
    err = ior(err,err_tmp)
    ymd = iyear * 10000 + imonth * 100 + iday
    tod = ihour * 3600 + iminute * 60 + isecond

    call seq_timemgr_EClockGetData(EClock, curr_ymd=curr_ymd, curr_tod=curr_tod)

    if (.not. seq_timemgr_EClockDateInSync( EClock, ymd, tod)) then
       call mpas_log_write('MPAS ymd=$i  MPAS tod=$i', MPAS_LOG_ERR, intArgs=(/ymd, tod/))
       call mpas_log_write('sync ymd=$i  sync tod=$i', MPAS_LOG_ERR, intArgs=(/curr_ymd, curr_tod/))
       call mpas_log_write('Internal mpas clock not in sync with sync clock', MPAS_LOG_ERR)
       call mpas_log_write("Internal MPAS clock not in sync with E3SM sync clock. Aborting...", MPAS_LOG_CRIT)
    end if

 end subroutine check_clocks_sync

 subroutine convert_seconds_to_timestamp(seconds, timeStamp)!{{{
   integer, intent(in) :: seconds
   character (len=StrKIND), intent(out) :: timeStamp
   real (kind=RKIND) :: secondsPerHour, secondsPerMinute, remaining
   integer :: minutes, hours, secondsLeft

   secondsPerHour = 3600
   secondsPerMinute = 60

   if(seconds < 0 .or. seconds > 86400) then
     secondsLeft = 00
     minutes = 00
     hours = 00
   else
     hours = int(seconds/secondsPerHour)
     remaining = seconds - real(hours) * secondsPerHour

     minutes = int(remaining/secondsPerMinute)
     remaining = remaining - real(minutes) * secondsPerMinute

     secondsLeft = int(remaining)
   end if

   write(timeStamp,"(a,i2.2,a,i2.2,a,i2.2)") "00_",hours,":",minutes,":",secondsLeft
   timeStamp = trim(timeStamp)

 end subroutine convert_seconds_to_timestamp!}}}

   subroutine add_stream_attributes(domain)!{{{

      type (domain_type), intent(inout) :: domain

      type (MPAS_Pool_iterator_type) :: itr
      integer, pointer :: intAtt
      logical, pointer :: logAtt
      character (len=StrKIND), pointer :: charAtt
      real (kind=RKIND), pointer :: realAtt
      character (len=StrKIND) :: histAtt
      character (len=StrKIND) :: caseid, model_version, hostname, username, &
                                 timeString
      character (len=8) :: curdate, curtime

      integer :: local_ierr

      if (domain % dminfo % nProcs < 10) then
          write(histAtt, '(A,I1,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
      else if (domain % dminfo % nProcs < 100) then
          write(histAtt, '(A,I2,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
      else if (domain % dminfo % nProcs < 1000) then
          write(histAtt, '(A,I3,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
      else if (domain % dminfo % nProcs < 10000) then
          write(histAtt, '(A,I4,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
      else if (domain % dminfo % nProcs < 100000) then
          write(histAtt, '(A,I5,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
      else
          write(histAtt, '(A,I6,A,A,A)') 'mpirun -n ', domain % dminfo % nProcs, ' ./', trim(domain % core % coreName), '_model'
      end if

      call seq_infodata_GetData(infodata,case_name=caseid, &
                                model_version=model_version, &
                                username=username, &
                                hostname=hostname)

      !Add E3SM run specific variables to MPAS output files
      call MPAS_stream_mgr_add_att(domain % streamManager, 'title', 'MPAS-Ocean output file information')
      call MPAS_stream_mgr_add_att(domain % streamManager, 'source', 'MPAS Ocean')
      call MPAS_stream_mgr_add_att(domain % streamManager, 'source_id', trim(model_version))
      call MPAS_stream_mgr_add_att(domain % streamManager, 'product', 'model-output')
      call MPAS_stream_mgr_add_att(domain % streamManager, 'realm', 'ocean')
      call MPAS_stream_mgr_add_att(domain % streamManager, 'case', trim(caseid))
      call MPAS_stream_mgr_add_att(domain % streamManager, 'username', trim(username))
      call MPAS_stream_mgr_add_att(domain % streamManager, 'hostname', trim(hostname))
      call MPAS_stream_mgr_add_att(domain % streamManager, 'git_version', trim(model_version))
      call datetime(curdate, curtime)
      timeString = 'created on ' // curdate // ' ' // curtime
      call MPAS_stream_mgr_add_att(domain % streamManager, 'history', trim(timeString))
      call MPAS_stream_mgr_add_att(domain % streamManager, 'Conventions', 'CF-1.7')
      call MPAS_stream_mgr_add_att(domain % streamManager, 'institution_id', 'E3SM-Project')
      call MPAS_stream_mgr_add_att(domain % streamManager, 'institution', &
          'LLNL (Lawrence Livermore National Laboratory, Livermore, CA 94550, USA); &
          &ANL (Argonne National Laboratory, Argonne, IL 60439, USA); &
          &BNL (Brookhaven National Laboratory, Upton, NY 11973, USA); &
          &LANL (Los Alamos National Laboratory, Los Alamos, NM 87545, USA); &
          &LBNL (Lawrence Berkeley National Laboratory, Berkeley, CA 94720, USA); &
          &ORNL (Oak Ridge National Laboratory, Oak Ridge, TN 37831, USA); &
          &PNNL (Pacific Northwest National Laboratory, Richland, WA 99352, USA); &
          &SNL (Sandia National Laboratories, Albuquerque, NM 87185, USA). &
          &Mailing address: LLNL Climate Program, c/o David C. Bader, &
          &Principal Investigator, L-103, 7000 East Avenue, Livermore, CA 94550, USA')
      call MPAS_stream_mgr_add_att(domain % streamManager, 'contact', &
          'e3sm-data-support@llnl.gov')

      call MPAS_stream_mgr_add_att(domain % streamManager, 'on_a_sphere', domain % on_a_sphere)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'sphere_radius', domain % sphere_radius)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'model_name', domain % core % modelName)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'core_name', domain % core % coreName)
      ! DWJ 10/01/2014: Eventually add the real history attribute, for now (due to length restrictions)
      ! add a shortened version.
      call MPAS_stream_mgr_add_att(domain % streamManager, 'parent_id', domain % parent_id)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'mesh_spec', domain % mesh_spec)

      call mpas_pool_begin_iteration(domain % configs)

      do while (mpas_pool_get_next_member(domain % configs, itr))

         if ( itr % memberType == MPAS_POOL_CONFIG) then

            if ( itr % dataType == MPAS_POOL_REAL ) then
               call mpas_pool_get_config(domain % configs, itr % memberName, realAtt)
               call MPAS_stream_mgr_add_att(domain % streamManager, itr % memberName, realAtt, ierr=local_ierr)
            else if ( itr % dataType == MPAS_POOL_INTEGER ) then
               call mpas_pool_get_config(domain % configs, itr % memberName, intAtt)
               call MPAS_stream_mgr_add_att(domain % streamManager, itr % memberName, intAtt, ierr=local_ierr)
            else if ( itr % dataType == MPAS_POOL_CHARACTER ) then
               call mpas_pool_get_config(domain % configs, itr % memberName, charAtt)
               call MPAS_stream_mgr_add_att(domain % streamManager, itr % memberName, charAtt, ierr=local_ierr)
            else if ( itr % dataType == MPAS_POOL_LOGICAL ) then
               call mpas_pool_get_config(domain % configs, itr % memberName, logAtt)
               if (logAtt) then
                  call MPAS_stream_mgr_add_att(domain % streamManager, itr % memberName, 'YES', ierr=local_ierr)
               else
                  call MPAS_stream_mgr_add_att(domain % streamManager, itr % memberName, 'NO', ierr=local_ierr)
               end if
            end if

          end if
      end do

   end subroutine add_stream_attributes!}}}

   subroutine datetime(cdate, ctime)!{{{

      implicit none

      character, intent(out) :: cdate*8, ctime*8

      integer, dimension(8) :: values
      character :: date*8, time*10, zone*5

      call date_and_time (date, time, zone, values)
      cdate(1:2) = date(5:6)
      cdate(3:3) = '/'
      cdate(4:5) = date(7:8)
      cdate(6:6) = '/'
      cdate(7:8) = date(3:4)
      ctime(1:2) = time(1:2)
      ctime(3:3) = ':'
      ctime(4:5) = time(3:4)
      ctime(6:6) = ':'
      ctime(7:8) = time(5:6)

   end subroutine datetime!}}}

#ifdef HAVE_MOAB


!***********************************************************************
!BOP
! !IROUTINE: ocn_import_moab
! !INTERFACE:

 subroutine ocn_import_moab(Eclock, errorCode)!{{{

! !DESCRIPTION:
!-----------------------------------------------------------------------
!  This routine receives message from moab driver
!
!    The following fields are always received from the driver:
!
!    o  taux   -- zonal wind stress (taux)                 (W/m2   )
!    o  tauy   -- meridonal wind stress (tauy)             (W/m2   )
!    o  snow   -- water flux due to snow                   (kg/m2/s)
!    o  rain   -- water flux due to rain                   (kg/m2/s)
!    o  evap   -- evaporation flux                         (kg/m2/s)
!    o  meltw  -- snow melt flux                           (kg/m2/s)
!    o  salt   -- salt                                     (kg(salt)/m2/s)
!    o  swnet  -- net short-wave heat flux                 (W/m2   )
!    o  sen    -- sensible heat flux                       (W/m2   )
!    o  lwup   -- longwave radiation (up)                  (W/m2   )
!    o  lwdn   -- longwave radiation (down)                (W/m2   )
!    o  melth  -- heat flux from snow&ice melt             (W/m2   )
!    o  ifrac  -- ice fraction                             (%)
!    o  rofl   -- river runoff flux                        (kg/m2/s)
!    o  rofi   -- ice runoff flux                          (kg/m2/s)
!    o  rofDIN -- DIN runoff flux                          (kg/m2/s)
!    o  rofDIP -- DIP runoff flux                          (kg/m2/s)
!    o  rofDON -- DON runoff flux                          (kg/m2/s)
!    o  rofDOP -- DOP runoff flux                          (kg/m2/s)
!    o  rofDOC -- DOC runoff flux                          (kg/m2/s)
!    o  rofPP  -- PP  runoff flux                          (kg/m2/s)
!    o  rofDSi -- DSi runoff flux                          (kg/m2/s)
!    o  rofPOC -- POC runoff flux                          (kg/m2/s)
!    o  rofPN  -- PN  runoff flux                          (kg/m2/s)
!    o  rofDIC -- DIC runoff flux                          (kg/m2/s)
!    o  rofFe  -- Fe  runoff flux                          (kg/m2/s)
!
!    The following fields are sometimes received from the coupler,
!      depending on model options:
!
!    o  pslv   -- atmospheric pressure at sea level        (Pa)
!    o  duu10n -- 10m wind speed squared                   (m^2/s^2)
!    o  co2prog-- bottom atm level prognostic co2
!    o  co2diag-- bottom atm level diagnostic co2
!
!-----------------------------------------------------------------------
!
! !REVISION HISTORY:
!  same as module
   use iMOAB,  only       : iMOAB_GetDoubleTagStorage, iMOAB_WriteMesh
! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock), intent(inout) :: EClock
   ! type(mct_aVect)   , intent(inout) :: x2o_o
   ! instead, we will get x2o_om from MPOID

! !OUTPUT PARAMETERS:

   integer, intent(out) :: &
      errorCode              ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   character (len=StrKIND) ::   &
      label,                 &
      message

   integer :: ent_type, ierr
   character(CXX) :: tagname

   integer ::  &
      i,n

   integer, pointer :: indexIT, indexIS

   logical, pointer :: config_use_ecosysTracers, &
                       config_use_ecosysTracers_sea_ice_coupling, &
                       config_use_DMSTracers,    &
                       config_use_DMSTracers_sea_ice_coupling,  &
                       config_use_MacroMoleculesTracers,  &
                       config_use_MacroMoleculesTracers_sea_ice_coupling, &
                       config_use_CFCTracers,  &
                       config_remove_ais_river_runoff, &
                       config_remove_ais_ice_runoff, &
                       config_cvmix_kpp_use_theory_wave

   character(len=StrKIND), pointer :: config_ecosys_atm_co2_option, &
                                      config_ecosys_atm_alt_co2_option

   real (kind=RKIND), pointer :: config_ecosys_atm_co2_constant_value

   real (kind=RKIND), pointer :: config_density0

   type (block_type), pointer :: block_ptr

   type (mpas_pool_type), pointer :: meshPool,             &
                                     forcingPool,          &
                                     ecosysAuxiliary,      &
                                     ecosysSeaIceCoupling, &
                                     DMSSeaIceCoupling,    &
                                     MacroMoleculesSeaIceCoupling,      &
                                     CFCAuxiliary

   integer, pointer :: nCellsSolve

   type (field1DReal), pointer :: windStressZonalField, windStressMeridionalField, &
                                  latentHeatFluxField, sensibleHeatFluxField, &
                                  longWaveHeatFluxUpField, &
                                  longWaveHeatFluxDownField, &
                                  evaporationFluxField, seaIceHeatFluxField, icebergHeatFluxField, &
                                  snowFluxField, seaIceFreshWaterFluxField, icebergFreshWaterFluxField, &
                                  seaIceSalinityFluxField, &
                                  riverRunoffFluxField, iceRunoffFluxField, &
                                  removedRiverRunoffFluxField, removedIceRunoffFluxField, &
                                  shortWaveHeatFluxField, rainFluxField, &
                                  atmosphericPressureField, iceFractionField, &
                                  seaIcePressureField, windSpeedSquared10mField, &
                                  atmosphericCO2Field, atmosphericCO2_ALT_CO2Field,  &
                                  iceFluxDICField, &
                                  iceFluxDONField, &
                                  iceFluxNO3Field, &
                                  iceFluxSiO3Field, &
                                  iceFluxNH4Field, &
                                  iceFluxDMSField, &
                                  iceFluxDMSPField, &
                                  iceFluxDOCrField, &
                                  iceFluxFeParticulateField, &
                                  iceFluxFeDissolvedField, &
                                  iceFluxDustField, &
                                  riverFluxNO3Field,  &
                                  riverFluxPO4Field,  &
                                  riverFluxSiO3Field, &
                                  riverFluxDOCField,  &
                                  riverFluxDONField,  &
                                  riverFluxDOPField,  &
                                  riverFluxDICField,  &
                                  riverFluxALKField,  &
                                  riverFluxFeField,   &
                                  landIceFreshwaterFluxField, &
                                  landIceHeatFluxField, &
                                  landIceFractionField, &
                                  windSpeed10mField, &
                                 significantWaveHeightField, &
                                  peakWaveFrequencyField, &
                                  peakWaveDirectionField
                                  !landIcePressureField

   type (field2DReal), pointer :: iceFluxPhytoCField, &
                                  iceFluxDOCField

   type (field2DReal), pointer :: landIceInterfaceTracersField

   type (field2DReal), pointer :: stokesDriftZonalWavenumberField, &
                                  stokesDriftMeridionalWavenumberField

   real (kind=RKIND), dimension(:), pointer :: windStressZonal, windStressMeridional, &
                                               latentHeatFlux, sensibleHeatFlux, &
                                               longWaveHeatFluxUp, &
                                               longWaveHeatFluxDown, &
                                               evaporationFlux, seaIceHeatFlux, icebergHeatFlux, &
                                               snowFlux, seaIceFreshWaterFlux, icebergFreshWaterFlux, &
                                               seaIceSalinityFlux, &
                                               riverRunoffFlux, iceRunoffFlux, &
                                               removedRiverRunoffFlux, removedIceRunoffFlux, &
                                               shortWaveHeatFlux, rainFlux, &
                                               atmosphericPressure, iceFraction, &
                                               seaIcePressure, windSpeedSquared10m, &
                                               atmosphericCO2, atmosphericCO2_ALT_CO2,  &
                                               windSpeedSquared10mCFC, &
                                               iceFluxDIC,       &
                                               iceFluxDON, &
                                               iceFluxNO3, &
                                               iceFluxSiO3, &
                                               iceFluxNH4, &
                                               iceFluxDMS, &
                                               iceFluxDMSP, &
                                               iceFluxDOCr, &
                                               iceFluxFeParticulate, &
                                               iceFluxFeDissolved, &
                                               iceFluxDust, &
                                               riverFluxNO3,  &
                                               riverFluxPO4,  &
                                               riverFluxSiO3, &
                                               riverFluxDOC,  &
                                               riverFluxDON,  &
                                               riverFluxDOP,  &
                                               riverFluxDIC,  &
                                               riverFluxALK,  &
                                               riverFluxFe,   &
                                               landIceFreshwaterFlux, &
                                               landIceHeatFlux, &
                                               landIceFraction, &
                                               areaCell, &
                                               windSpeed10m, &
                                               significantWaveHeight, &
                                               peakWaveFrequency, &
                                               peakWaveDirection
                                               !landIcePressure

   real (kind=RKIND), dimension(:), pointer :: latCell

   real (kind=RKIND), dimension(:,:), pointer :: iceFluxPhytoC, &
                                                 iceFluxDOC, &
                                                 stokesDriftZonalWavenumber, &
                                                 stokesDriftMeridionalWavenumber

   real (kind=RKIND) :: removedRiverRunoffFluxThisProc, removedIceRunoffFluxThisProc
   real (kind=RKIND) :: removedRiverRunoffFluxReduced, removedIceRunoffFluxReduced

   real (kind=RKIND), dimension(:,:), pointer :: landIceInterfaceTracers

   real (kind=RKIND) :: riverFactor

!-----------------------------------------------------------------------
!
!  zero out padded cells
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer :: cur_ocn_stepno
#ifdef MOABDEBUG
    character*100 outfile, wopts, lnum
#endif

    call seq_timemgr_EClockGetData( EClock, stepno=cur_ocn_stepno )
#ifdef MOABDEBUG
    write(lnum,"(I0.2)")cur_ocn_stepno
    outfile = 'ocn_import_'//trim(lnum)//'.h5m'//C_NULL_CHAR
    wopts   = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
    ierr = iMOAB_WriteMesh(MPOID, outfile, wopts)
    if (ierr > 0 )  then
       write(ocnLogUnit,*) 'Fail to write ocean state '
    endif
#endif

   errorCode = 0

!  get moab tags from MPOID
   ent_type = 1 ! cells
   !  get all tags in one method
   tagname = trim(seq_flds_x2o_fields)//C_NULL_CHAR
   ierr = iMOAB_GetDoubleTagStorage ( MPOID, tagname, totalmbls_r , ent_type, x2o_om(1, 1) )
   if ( ierr /= 0 ) then
      write(ocnLogUnit,*) 'Fail to get MOAB fields '
   endif
!-----------------------------------------------------------------------
!
!  unpack and distribute wind stress, then convert to correct units
!  and rotate components to local coordinates
!
!-----------------------------------------------------------------------

   ! get configure options
   call mpas_pool_get_config(domain % configs, 'config_use_ecosysTracers', config_use_ecosysTracers)
   call mpas_pool_get_config(domain % configs, 'config_use_DMSTracers', config_use_DMSTracers)
   call mpas_pool_get_config(domain % configs, 'config_use_MacroMoleculesTracers', config_use_MacroMoleculesTracers)
   call mpas_pool_get_config(domain % configs, 'config_density0', config_density0)
   call mpas_pool_get_config(domain % configs, 'config_use_ecosysTracers_sea_ice_coupling',  &
                                                config_use_ecosysTracers_sea_ice_coupling)
   call mpas_pool_get_config(domain % configs, 'config_use_DMSTracers_sea_ice_coupling',  &
                                                config_use_DMSTracers_sea_ice_coupling)
   call mpas_pool_get_config(domain % configs, 'config_use_MacroMoleculesTracers_sea_ice_coupling',  &
                                                config_use_MacroMoleculesTracers_sea_ice_coupling)
   call mpas_pool_get_config(domain % configs, 'config_use_CFCTracers', config_use_CFCTracers)
   call mpas_pool_get_config(domain % configs, 'config_remove_ais_river_runoff', config_remove_ais_river_runoff)
   call mpas_pool_get_config(domain % configs, 'config_remove_ais_ice_runoff', config_remove_ais_ice_runoff)
   call mpas_pool_get_config(domain % configs, 'config_cvmix_kpp_use_theory_wave', config_cvmix_kpp_use_theory_wave)

   n = 0
   removedRiverRunoffFluxThisProc = 0.0_RKIND
   removedIceRunoffFluxThisProc = 0.0_RKIND

   block_ptr => domain % blocklist
   do while(associated(block_ptr))
      call mpas_pool_get_subpool(block_ptr % structs, 'mesh',    meshPool)
      call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)

      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

      call mpas_pool_get_field(forcingPool, 'windStressZonal', windStressZonalField)
      call mpas_pool_get_field(forcingPool, 'windStressMeridional', windStressMeridionalField)
      call mpas_pool_get_field(forcingPool, 'latentHeatFlux', latentHeatFluxField)
      call mpas_pool_get_field(forcingPool, 'sensibleHeatFlux', sensibleHeatFluxField)
      call mpas_pool_get_field(forcingPool, 'longWaveHeatFluxUp', longWaveHeatFluxUpField)
      call mpas_pool_get_field(forcingPool, 'longWaveHeatFluxDown', longWaveHeatFluxDownField)
      call mpas_pool_get_field(forcingPool, 'evaporationFlux', evaporationFluxField)
      call mpas_pool_get_field(forcingPool, 'seaIceHeatFlux', seaIceHeatFluxField)
      call mpas_pool_get_field(forcingPool, 'icebergHeatFlux', icebergHeatFluxField)
      call mpas_pool_get_field(forcingPool, 'snowFlux', snowFluxField)
      call mpas_pool_get_field(forcingPool, 'seaIceFreshWaterFlux', seaIceFreshWaterFluxField)
      call mpas_pool_get_field(forcingPool, 'icebergFreshWaterFlux', icebergFreshWaterFluxField)
      call mpas_pool_get_field(forcingPool, 'seaIceSalinityFlux', seaIceSalinityFluxField)
      call mpas_pool_get_field(forcingPool, 'riverRunoffFlux', riverRunoffFluxField)
      call mpas_pool_get_field(forcingPool, 'iceRunoffFlux', iceRunoffFluxField)
      call mpas_pool_get_field(forcingPool, 'shortWaveHeatFlux', shortWaveHeatFluxField)
      call mpas_pool_get_field(forcingPool, 'rainFlux', rainFluxField)
      call mpas_pool_get_field(forcingPool, 'atmosphericPressure', atmosphericPressureField)
      call mpas_pool_get_field(forcingPool, 'seaIcePressure', seaIcePressureField)
      call mpas_pool_get_field(forcingPool, 'iceFraction', iceFractionField)
      call mpas_pool_get_field(forcingPool, 'iceRunoffFlux', iceRunoffFluxField)
      call mpas_pool_get_field(forcingPool, 'removedRiverRunoffFlux', removedRiverRunoffFluxField)
      call mpas_pool_get_field(forcingPool, 'removedIceRunoffFlux', removedIceRunoffFluxField)
      call mpas_pool_get_field(forcingPool, 'stokesDriftZonalWavenumber', stokesDriftZonalWavenumberField)
      call mpas_pool_get_field(forcingPool, 'stokesDriftMeridionalWavenumber', stokesDriftMeridionalWavenumberField)
      call mpas_pool_get_field(forcingPool, 'significantWaveHeight', significantWaveHeightField)
      call mpas_pool_get_field(forcingPool, 'peakWaveFrequency', peakWaveFrequencyField)
      call mpas_pool_get_field(forcingPool, 'peakWaveDirection', peakWaveDirectionField)

      call mpas_pool_get_field(forcingPool, 'landIceFreshwaterFlux', landIceFreshwaterFluxField)
      call mpas_pool_get_field(forcingPool, 'landIceHeatFlux', landIceHeatFluxField)
      call mpas_pool_get_field(forcingPool, 'landIceFraction', landIceFractionField)
      call mpas_pool_get_field(forcingPool, 'landIceInterfaceTracers', landIceInterfaceTracersField)

      call mpas_pool_get_field(forcingPool, 'windSpeed10m', windSpeed10mField)

      call mpas_pool_get_dimension(forcingPool, 'index_landIceInterfaceTemperature', indexIT)
      call mpas_pool_get_dimension(forcingPool, 'index_landIceInterfaceSalinity', indexIS)

   !call mpas_pool_get_field(forcingPool, 'landIcePressure', landIcePressureField)

      windStressZonal => windStressZonalField % array
      windStressMeridional => windStressMeridionalField % array
      latentHeatFlux => latentHeatFluxField % array
      sensibleHeatFlux => sensibleHeatFluxField % array
      longWaveHeatFluxUp => longWaveHeatFluxUpField % array
      longWaveHeatFluxDown => longWaveHeatFluxDownField % array
      evaporationFlux => evaporationFluxField % array
      seaIceHeatFlux => seaIceHeatFluxField % array
      icebergHeatFlux => icebergHeatFluxField % array
      snowFlux => snowFluxField % array
      seaIceFreshWaterFlux => seaIceFreshWaterFluxField % array
      icebergFreshWaterFlux => icebergFreshWaterFluxField % array
      seaIceSalinityFlux => seaIceSalinityFluxField % array
      riverRunoffFlux => riverRunoffFluxField % array
      iceRunoffFlux => iceRunoffFluxField % array
      shortWaveHeatFlux => shortWaveHeatFluxField % array
      rainFlux => rainFluxField % array
      atmosphericPressure => atmosphericPressureField % array
      seaIcePressure => seaIcePressureField % array
      iceFraction => iceFractionField % array
      iceRunoffFlux => iceRunoffFluxField % array
      removedRiverRunoffFlux => removedRiverRunoffFluxField % array
      removedIceRunoffFlux => removedIceRunoffFluxField % array
      landIceFreshwaterFlux => landIceFreshwaterFluxField % array
      landIceHeatFlux => landIceHeatFluxField % array
      landIceInterfaceTracers => landIceInterfaceTracersField % array
      landIceFraction => landIceFractionField % array
      windSpeed10m => windSpeed10mField % array
      stokesDriftZonalWavenumber => stokesDriftZonalWavenumberField % array
      stokesDriftMeridionalWavenumber => stokesDriftMeridionalWavenumberField % array
      significantWaveHeight => significantWaveHeightField % array
      peakWaveFrequency => peakWaveFrequencyField % array
      peakWaveDirection => peakWaveDirectionField % array
      !landIcePressure => landIcePressureField % array

      call mpas_pool_get_array(meshPool, 'latCell', latCell)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

      ! BGC fields
      if (config_use_ecosysTracers) then
         call mpas_pool_get_subpool(forcingPool, 'ecosysAuxiliary', ecosysAuxiliary)

         call mpas_pool_get_field(ecosysAuxiliary, 'windSpeedSquared10m', windSpeedSquared10mField)
         windSpeedSquared10m => windSpeedSquared10mField % array
         call mpas_pool_get_field(ecosysAuxiliary, 'atmosphericCO2', atmosphericCO2Field)
         atmosphericCO2 => atmosphericCO2Field % array
         call mpas_pool_get_field(ecosysAuxiliary, 'atmosphericCO2_ALT_CO2', atmosphericCO2_ALT_CO2Field)
         atmosphericCO2_ALT_CO2 => atmosphericCO2_ALT_CO2Field % array

         if (config_use_ecosysTracers_river_inputs_from_coupler) then
           call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxNO3' , riverFluxNO3Field)
           call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxPO4' , riverFluxPO4Field)
           call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxDON' , riverFluxDONField)
           call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxDOP' , riverFluxDOPField)
           call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxSiO3', riverFluxSiO3Field)
           call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxDOC' , riverFluxDOCField)
           call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxDIC' , riverFluxDICField)
           call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxALK' , riverFluxALKField)
           call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxFe'  , riverFluxFeField)
           riverFluxNO3  => riverFluxNO3Field  % array
           riverFluxPO4  => riverFluxPO4Field  % array
           riverFluxDON  => riverFluxDONField  % array
           riverFluxDOP  => riverFluxDOPField  % array
           riverFluxSiO3 => riverFluxSiO3Field % array
           riverFluxDOC  => riverFluxDOCField  % array
           riverFluxDIC  => riverFluxDICField  % array
           riverFluxALK  => riverFluxALKField  % array
           riverFluxFe   => riverFluxFeField   % array
         endif

         call mpas_pool_get_config(domain % configs, 'config_ecosys_atm_co2_option',   &
            config_ecosys_atm_co2_option)
         call mpas_pool_get_config(domain % configs, 'config_ecosys_atm_alt_co2_option',   &
            config_ecosys_atm_alt_co2_option)
         call mpas_pool_get_config(domain % configs, 'config_ecosys_atm_co2_constant_value',   &
            config_ecosys_atm_co2_constant_value)

         if (config_use_ecosysTracers_sea_ice_coupling) then
            call mpas_pool_get_subpool(forcingPool, 'ecosysSeaIceCoupling', ecosysSeaIceCoupling)

            call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxPhytoC', iceFluxPhytoCField)
            call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxDIC', iceFluxDICField)
            call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxNO3', iceFluxNO3Field)
            call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxSiO3', iceFluxSiO3Field)
            call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxNH4', iceFluxNH4Field)
            call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxDOCr', iceFluxDOCrField)
            call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxFeParticulate', iceFluxFeParticulateField)
            call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxFeDissolved', iceFluxFeDissolvedField)
            call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxDust', iceFluxDustField)
            call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxDOC', iceFluxDOCField)
            call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxDON', iceFluxDONField)

            iceFluxPhytoC => iceFluxPhytoCField % array
            iceFluxDIC => iceFluxDICField % array
            iceFluxNO3 => iceFluxNO3Field % array
            iceFluxSiO3 => iceFluxSiO3Field % array
            iceFluxNH4 => iceFluxNH4Field % array
            iceFluxDOCr => iceFluxDOCrField % array
            iceFluxFeParticulate => iceFluxFeParticulateField % array
            iceFluxFeDissolved => iceFluxFeDissolvedField % array
            iceFluxDust => iceFluxDustField % array
            iceFluxDOC => iceFluxDOCField % array
            iceFluxDON => iceFluxDONField % array
         endif
      endif
      if (config_use_DMSTracers .and. config_use_DMSTracers_sea_ice_coupling) then
         call mpas_pool_get_subpool(forcingPool, 'DMSSeaIceCoupling', DMSSeaIceCoupling)

         call mpas_pool_get_field(DMSSeaIceCoupling, 'iceFluxDMS', iceFluxDMSField)
         call mpas_pool_get_field(DMSSeaIceCoupling, 'iceFluxDMSP', iceFluxDMSPField)

         iceFluxDMS => iceFluxDMSField % array
         iceFluxDMSP => iceFluxDMSPField % array
      endif

      ! CFC fields
      if (config_use_CFCTracers) then
         call mpas_pool_get_subpool(forcingPool, 'CFCAuxiliary', CFCAuxiliary)
         call mpas_pool_get_field(CFCAuxiliary, 'windSpeedSquared10mCFC', windSpeedSquared10mField)
         windSpeedSquared10mCFC => windSpeedSquared10mField % array
      endif

      if (config_remove_ais_river_runoff) then
         ! Initialize this field
         removedRiverRunoffFlux(:) = 0.0_RKIND
      endif

      if (config_remove_ais_ice_runoff) then
         ! Initialize this field
         removedIceRunoffFlux(:) = 0.0_RKIND
      endif

      if (config_cvmix_kpp_use_theory_wave) then
         ! Initialize this field
         windSpeed10m(:) = 0.0_RKIND
      endif
! !  replace 'x2o_o % rAttr(' to 'x2o_om(n, '  and ', n)' with ')'
      do i = 1, nCellsSolve
        n = n + 1
        if ( windStressZonalField % isActive ) then
           windStressZonal(i) = x2o_om(n, index_x2o_Foxx_taux)
        end if
        if ( windStressMeridionalField % isActive ) then
           windStressMeridional(i) = x2o_om(n, index_x2o_Foxx_tauy)
        end if

        if ( latentHeatFluxField % isActive ) then
           latentHeatFlux(i) = x2o_om(n, index_x2o_Foxx_lat)
        end if
        if ( sensibleHeatFluxField % isActive ) then
           sensibleHeatFlux(i) = x2o_om(n, index_x2o_Foxx_sen)
        end if
        if ( longWaveHeatFluxUpField % isActive ) then
           longWaveHeatFluxUp(i) = x2o_om(n, index_x2o_Foxx_lwup)
        end if
        if ( longWaveHeatFluxDownField % isActive ) then
           longWaveHeatFluxDown(i) = x2o_om(n, index_x2o_Faxa_lwdn)
        end if
        if ( evaporationFluxField % isActive ) then
           evaporationFlux(i) = x2o_om(n, index_x2o_Foxx_evap)
        end if
        if ( seaIceHeatFluxField % isActive ) then
           seaIceHeatFlux(i) = x2o_om(n, index_x2o_Fioi_melth)
        end if
        if ( icebergHeatFluxField % isActive ) then
           icebergHeatFlux(i) = x2o_om(n, index_x2o_Fioi_bergh)
        end if
        if ( snowFluxField % isActive ) then
           snowFlux(i) = x2o_om(n, index_x2o_Faxa_snow)
        end if
        if ( seaIceFreshWaterFluxField % isActive ) then
           seaIceFreshWaterFlux(i) = x2o_om(n, index_x2o_Fioi_meltw)
        end if
        if ( icebergFreshWaterFluxField % isActive ) then
           icebergFreshWaterFlux(i) = x2o_om(n, index_x2o_Fioi_bergw)
        end if
        if ( seaIceSalinityFluxField % isActive ) then
           seaIceSalinityFlux(i) = x2o_om(n, index_x2o_Fioi_salt)
        end if
        if ( riverRunoffFluxField % isActive ) then
           riverRunoffFlux(i) = x2o_om(n, index_x2o_Foxx_rofl)
           if (config_remove_ais_river_runoff) then
              if (latCell(i) < -1.04719666667_RKIND) then ! 60S in radians
                 removedRiverRunoffFlux(i) = riverRunoffFlux(i)
                 riverRunoffFlux(i) = 0.0_RKIND
                 removedRiverRunoffFluxThisProc = removedRiverRunoffFluxThisProc + removedRiverRunoffFlux(i)
               endif
           endif
        end if
        if ( iceRunoffFluxField % isActive ) then
           iceRunoffFlux(i) = x2o_om(n, index_x2o_Foxx_rofi)
           if(iceRunoffFlux(n) < 0.0_RKIND) then
               call shr_sys_abort ('Error: incoming rofi_F is negative')
           end if
           if (config_remove_ais_ice_runoff) then
              if (latCell(i) < -0.99483767345_RKIND) then ! 57S in radians
                 removedIceRunoffFlux(i) = iceRunoffFlux(i)
                 iceRunoffFlux(i) = 0.0_RKIND
                 removedIceRunoffFluxThisProc = removedIceRunoffFluxThisProc + removedIceRunoffFlux(i)
              endif
           endif
        end if
        if ( shortWaveHeatFluxField % isActive ) then
           shortWaveHeatFlux(i) = max(x2o_om(n, index_x2o_Foxx_swnet), 0.0_RKIND)
        end if

        if ( rainFluxField % isActive ) then
           rainFlux(i) = x2o_om(n, index_x2o_Faxa_rain)
        end if
        if ( atmosphericPressureField % isActive ) then
           atmosphericPressure(i) = x2o_om(n, index_x2o_Sa_pslv)
        end if
        if ( seaIcePressureField % isActive ) then
           ! Set seaIcePressure to be limited to 5m of pressure
           seaIcePressure(i) =  min( x2o_om(n, index_x2o_Si_bpress), config_density0 * gravity * 5.0_RKIND )
        end if
        if ( iceFractionField % isActive ) then
           iceFraction(i) = x2o_om(n, index_x2o_Si_ifrac)
        end if

        if ( stokesDriftZonalWavenumberField % isActive ) then
           stokesDriftZonalWavenumber(1,i) = x2o_om(n, index_x2o_Sw_ustokes_wavenumber_1)
           stokesDriftZonalWavenumber(2,i) = x2o_om(n, index_x2o_Sw_ustokes_wavenumber_2)
           stokesDriftZonalWavenumber(3,i) = x2o_om(n, index_x2o_Sw_ustokes_wavenumber_3)
           stokesDriftZonalWavenumber(4,i) = x2o_om(n, index_x2o_Sw_ustokes_wavenumber_4)
           stokesDriftZonalWavenumber(5,i) = x2o_om(n, index_x2o_Sw_ustokes_wavenumber_5)
           stokesDriftZonalWavenumber(6,i) = x2o_om(n, index_x2o_Sw_ustokes_wavenumber_6)
        end if
        if ( stokesDriftMeridionalWavenumberField % isActive ) then
           stokesDriftMeridionalWavenumber(1,i) = x2o_om(n, index_x2o_Sw_vstokes_wavenumber_1)
           stokesDriftMeridionalWavenumber(2,i) = x2o_om(n, index_x2o_Sw_vstokes_wavenumber_2)
           stokesDriftMeridionalWavenumber(3,i) = x2o_om(n, index_x2o_Sw_vstokes_wavenumber_3)
           stokesDriftMeridionalWavenumber(4,i) = x2o_om(n, index_x2o_Sw_vstokes_wavenumber_4)
           stokesDriftMeridionalWavenumber(5,i) = x2o_om(n, index_x2o_Sw_vstokes_wavenumber_5)
           stokesDriftMeridionalWavenumber(6,i) = x2o_om(n, index_x2o_Sw_vstokes_wavenumber_6)
        end if
        if ( significantWaveHeightField % isActive ) then
           significantWaveHeight(i) = x2o_om(n, index_x2o_Sw_Hs)
        end if
        if ( peakWaveFrequencyField % isActive ) then
           peakWaveFrequency(i) = x2o_om(n, index_x2o_Sw_Fp)
        end if
        if ( peakWaveDirectionField % isActive ) then
           peakWaveDirection(i) = x2o_om(n, index_x2o_Sw_Dp)
        end if

        if (config_cvmix_kpp_use_theory_wave) then
           if ( windSpeed10mField% isActive ) then
              windSpeed10m(i) = sqrt( x2o_om(n, index_x2o_So_duu10n))
           end if
        endif

        if ( landIceFreshwaterFluxField % isActive ) then
           !landIceFreshwaterFlux(i) = x2o_om(n, index_x2o_Fogx_qicelo)
        end if
        if ( landIceHeatFluxField % isActive ) then
           !landIceHeatFlux(i) = x2o_om(n, index_x2o_Fogx_qiceho)
        end if
        if ( landIceInterfaceTracersField % isActive ) then
           !landIceInterfaceTracers(indexIT, i) = x2o_om(n, index_x2o_Sg_blit)
           !landIceInterfaceTracers(indexIS, i) = x2o_om(n, index_x2o_Sg_blis)
        end if
        if ( landIceFractionField % isActive ) then
           !landIceFraction(i) = x2o_om(n, index_x2o_Sg_icemask)
        end if
        !if ( landIcePressureField % isActive ) then
           !landIcePressure(i) = x2o_om(n, index_x2o_Sg_lithop)
        !end if

        ! BGC fields
        if (config_use_ecosysTracers) then
           if ( windSpeedSquared10mField % isActive ) then
              windSpeedSquared10m(i) = x2o_om(n, index_x2o_So_duu10n)
           end if
           if ( atmosphericCO2Field % isActive ) then
              if ( config_ecosys_atm_co2_option == 'prognostic') then
                 atmosphericCO2(i) = x2o_om(n, index_x2o_Sa_co2prog)
              else if ( config_ecosys_atm_co2_option == 'diagnostic') then
                 atmosphericCO2(i) = x2o_om(n, index_x2o_Sa_co2diag)
              else if ( config_ecosys_atm_co2_option == 'bcrc') then
                 atmosphericCO2(i) = config_ecosys_atm_co2_constant_value
              else if ( config_ecosys_atm_co2_option == 'bcrd') then
                 atmosphericCO2(i) = config_ecosys_atm_co2_constant_value
              else if ( config_ecosys_atm_co2_option == 'bdrc') then
                 atmosphericCO2(i) = x2o_om(n, index_x2o_Sa_co2diag)
              else if ( config_ecosys_atm_co2_option == 'bdrd') then
                 atmosphericCO2(i) = x2o_om(n, index_x2o_Sa_co2diag)
              else
                 atmosphericCO2(i) = config_ecosys_atm_co2_constant_value
              end if
           end if
           if ( atmosphericCO2_ALT_CO2Field % isActive ) then
              if ( config_ecosys_atm_alt_co2_option == 'prognostic') then
                 atmosphericCO2_ALT_CO2(i) = x2o_om(n, index_x2o_Sa_co2prog)
              else if ( config_ecosys_atm_alt_co2_option == 'diagnostic') then
                 atmosphericCO2_ALT_CO2(i) = x2o_om(n, index_x2o_Sa_co2diag)
              else if ( config_ecosys_atm_alt_co2_option == 'bcrc') then
                 atmosphericCO2_ALT_CO2(i) = config_ecosys_atm_co2_constant_value
              else if ( config_ecosys_atm_alt_co2_option == 'bcrd') then
                 atmosphericCO2_ALT_CO2(i) = x2o_om(n, index_x2o_Sa_co2diag)
              else if ( config_ecosys_atm_alt_co2_option == 'bdrc') then
                 atmosphericCO2_ALT_CO2(i) = config_ecosys_atm_co2_constant_value
              else if ( config_ecosys_atm_alt_co2_option == 'bdrd') then
                 atmosphericCO2_ALT_CO2(i) = x2o_om(n, index_x2o_Sa_co2diag)
              else
                 atmosphericCO2_ALT_CO2(i) = config_ecosys_atm_co2_constant_value
              end if
           end if

           if (config_use_ecosysTracers_river_inputs_from_coupler) then
             riverFluxNO3(i)  = x2o_om(n, index_x2o_Foxx_rofDIN)
             riverFluxPO4(i)  = x2o_om(n, index_x2o_Foxx_rofDIP)
             riverFluxDON(i)  = x2o_om(n, index_x2o_Foxx_rofDON)
             riverFluxDOP(i)  = x2o_om(n, index_x2o_Foxx_rofDOP)
             riverFluxSiO3(i) = x2o_om(n, index_x2o_Foxx_rofDSi)
             riverFluxDOC(i)  = x2o_om(n, index_x2o_Foxx_rofDOC)
             riverFluxDIC(i)  = x2o_om(n, index_x2o_Foxx_rofDIC)
             riverFluxFe(i)   = x2o_om(n, index_x2o_Foxx_rofFe )

! convert from kgNutrient/(m2-s) to mmol/m3 m/s
             riverFactor = 1.e6_RKIND
             riverFluxNO3(i)  = riverFluxNO3(i)*riverFactor/14.007_RKIND
             riverFluxPO4(i)  = riverFluxPO4(i)*riverFactor/30.974_RKIND
             riverFluxDON(i)  = riverFluxDON(i)*riverFactor/14.007_RKIND
             riverFluxDOP(i)  = riverFluxDOP(i)*riverFactor/30.974_RKIND
             riverFluxSiO3(i) = riverFluxSiO3(i)*riverFactor/28.085_RKIND
             riverFluxDOC(i)  = riverFluxDOC(i)*riverFactor/12.001_RKIND
             riverFluxDIC(i)  = riverFluxDIC(i)*riverFactor/12.001_RKIND
             riverFluxFe(i)   = riverFluxFe(i)*riverFactor/55.845_RKIND

             riverFluxALK(i)  = riverFluxDIC(i)
           endif

           if (config_use_ecosysTracers_sea_ice_coupling) then
              if ( iceFluxPhytoCField % isActive ) then
                 iceFluxPhytoC(1,i) = x2o_om(n, index_x2o_Fioi_algae1)
                 iceFluxPhytoC(2,i) = x2o_om(n, index_x2o_Fioi_algae2)
                 iceFluxPhytoC(3,i) = x2o_om(n, index_x2o_Fioi_algae3)
              endif
              if ( iceFluxDICField % isActive ) then
                 iceFluxDIC(i) = x2o_om(n, index_x2o_Fioi_dic1)
              endif
              if ( iceFluxNO3Field % isActive ) then
                 iceFluxNO3(i) = x2o_om(n, index_x2o_Fioi_no3)
              endif
              if ( iceFluxSiO3Field % isActive ) then
                 iceFluxSiO3(i) = x2o_om(n, index_x2o_Fioi_sio3)
              endif
              if ( iceFluxNH4Field % isActive ) then
                 iceFluxNH4(i) = x2o_om(n, index_x2o_Fioi_nh4)
              endif
              if ( iceFluxDOCrField % isActive ) then
                 iceFluxDOCr(i) = x2o_om(n, index_x2o_Fioi_docr)
              endif
              if ( iceFluxFeParticulateField % isActive ) then
                 iceFluxFeParticulate(i) = x2o_om(n, index_x2o_Fioi_fep1)
              endif
              if ( iceFluxFeDissolvedField % isActive ) then
                 iceFluxFeDissolved(i) = x2o_om(n, index_x2o_Fioi_fed1)
              endif
              if ( iceFluxDustField % isActive ) then
                 iceFluxDust(i) = x2o_om(n, index_x2o_Fioi_dust1)
              endif
              if ( iceFluxDOCField % isActive ) then
                 iceFluxDOC(1,i) = x2o_om(n, index_x2o_Fioi_doc1)
                 iceFluxDOC(2,i) = x2o_om(n, index_x2o_Fioi_doc2)
                 iceFluxDOC(3,i) = x2o_om(n, index_x2o_Fioi_doc3)
              endif
              if ( iceFluxDONField % isActive ) then
                 iceFluxDON(i) = x2o_om(n, index_x2o_Fioi_don1)
              endif
            endif
         endif
         if (config_use_DMSTracers .and. config_use_DMSTracers_sea_ice_coupling) then
           if ( iceFluxDMSField % isActive ) then
              iceFluxDMS(i) = x2o_om(n, index_x2o_Fioi_dms)
           endif
           if ( iceFluxDMSPField % isActive ) then
              !JW TODO: dmspp? dmspd? the sum?
              iceFluxDMSP(i) = x2o_om(n, index_x2o_Fioi_dmspp)
           endif
         endif

        ! CFC fields
        if (config_use_CFCTracers) then
           if ( windSpeedSquared10mField % isActive ) then
              windSpeedSquared10mCFC(i) = x2o_om(n, index_x2o_So_duu10n)
           end if
        end if

      end do

      block_ptr => block_ptr % next
   end do

   call mpas_pool_get_subpool(domain % blocklist % structs, 'forcing', forcingPool)

   call mpas_pool_get_field(forcingPool, 'windStressZonal', windStressZonalField)
   call mpas_pool_get_field(forcingPool, 'windStressMeridional', windStressMeridionalField)
   call mpas_pool_get_field(forcingPool, 'latentHeatFlux', latentHeatFluxField)
   call mpas_pool_get_field(forcingPool, 'sensibleHeatFlux', sensibleHeatFluxField)
   call mpas_pool_get_field(forcingPool, 'longWaveHeatFluxUp', longWaveHeatFluxUpField)
   call mpas_pool_get_field(forcingPool, 'longWaveHeatFluxDown', longWaveHeatFluxDownField)
   call mpas_pool_get_field(forcingPool, 'evaporationFlux', evaporationFluxField)
   call mpas_pool_get_field(forcingPool, 'seaIceHeatFlux', seaIceHeatFluxField)
   call mpas_pool_get_field(forcingPool, 'icebergHeatFlux', icebergHeatFluxField)
   call mpas_pool_get_field(forcingPool, 'snowFlux', snowFluxField)
   call mpas_pool_get_field(forcingPool, 'seaIceFreshWaterFlux', seaIceFreshWaterFluxField)
   call mpas_pool_get_field(forcingPool, 'icebergFreshWaterFlux', icebergFreshWaterFluxField)
   call mpas_pool_get_field(forcingPool, 'seaIceSalinityFlux', seaIceSalinityFluxField)
   call mpas_pool_get_field(forcingPool, 'riverRunoffFlux', riverRunoffFluxField)
   call mpas_pool_get_field(forcingPool, 'iceRunoffFlux', iceRunoffFluxField)
   call mpas_pool_get_field(forcingPool, 'shortWaveHeatFlux', shortWaveHeatFluxField)
   call mpas_pool_get_field(forcingPool, 'rainFlux', rainFluxField)
   call mpas_pool_get_field(forcingPool, 'atmosphericPressure', atmosphericPressureField)
   call mpas_pool_get_field(forcingPool, 'seaIcePressure', seaIcePressureField)
   call mpas_pool_get_field(forcingPool, 'iceFraction', iceFractionField)
   call mpas_pool_get_field(forcingPool, 'stokesDriftZonalWavenumber', stokesDriftZonalWavenumberField)
   call mpas_pool_get_field(forcingPool, 'stokesDriftMeridionalWavenumber', stokesDriftMeridionalWavenumberField)
   call mpas_pool_get_field(forcingPool, 'significantWaveHeight', significantWaveHeightField)
   call mpas_pool_get_field(forcingPool, 'peakWaveFrequency', peakWaveFrequencyField)
   call mpas_pool_get_field(forcingPool, 'peakWaveDirection', peakWaveDirectionField)

   call mpas_pool_get_field(forcingPool, 'landIceFreshwaterFlux', landIceFreshwaterFluxField)
   call mpas_pool_get_field(forcingPool, 'landIceHeatFlux', landIceHeatFluxField)
   call mpas_pool_get_field(forcingPool, 'landIceFraction', landIceFractionField)
   call mpas_pool_get_field(forcingPool, 'landIceInterfaceTracers', landIceInterfaceTracersField)

   call mpas_pool_get_field(forcingPool, 'windSpeed10m', windSpeed10mField)

   call mpas_pool_get_dimension(forcingPool, 'index_landIceInterfaceTemperature', indexIT)
   call mpas_pool_get_dimension(forcingPool, 'index_landIceInterfaceSalinity', indexIS)
   !call mpas_pool_get_field(forcingPool, 'landIcePressure', landIcePressureField)

   ! BGC fields
   if (config_use_ecosysTracers) then
      call mpas_pool_get_subpool(forcingPool, 'ecosysAuxiliary', ecosysAuxiliary)

      call mpas_pool_get_field(ecosysAuxiliary, 'windSpeedSquared10m', windSpeedSquared10mField)
      call mpas_pool_get_field(ecosysAuxiliary, 'atmosphericCO2', atmosphericCO2Field)
      call mpas_pool_get_field(ecosysAuxiliary, 'atmosphericCO2_ALT_CO2', atmosphericCO2_ALT_CO2Field)

      if (config_use_ecosysTracers_river_inputs_from_coupler) then
        call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxNO3' , riverFluxNO3Field)
        call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxPO4' , riverFluxPO4Field)
        call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxDON' , riverFluxDONField)
        call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxDOP' , riverFluxDOPField)
        call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxSiO3', riverFluxSiO3Field)
        call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxDOC' , riverFluxDOCField)
        call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxDIC' , riverFluxDICField)
        call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxALK' , riverFluxALKField)
        call mpas_pool_get_field(ecosysAuxiliary, 'riverFluxFe'  , riverFluxFeField)
      endif

      if (config_use_ecosysTracers_sea_ice_coupling) then
         call mpas_pool_get_subpool(forcingPool, 'ecosysSeaIceCoupling', ecosysSeaIceCoupling)

         call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxPhytoC', iceFluxPhytoCField)
         call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxDIC', iceFluxDICField)
         call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxNO3', iceFluxNO3Field)
         call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxSiO3', iceFluxSiO3Field)
         call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxNH4', iceFluxNH4Field)
         call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxDOCr', iceFluxDOCrField)
         call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxFeParticulate', iceFluxFeParticulateField)
         call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxFeDissolved', iceFluxFeDissolvedField)
         call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxDust', iceFluxDustField)
         call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxDOC', iceFluxDOCField)
         call mpas_pool_get_field(ecosysSeaIceCoupling, 'iceFluxDON', iceFluxDONField)
      endif
   endif
   if (config_use_DMSTracers .and. config_use_DMSTracers_sea_ice_coupling) then
      call mpas_pool_get_subpool(forcingPool, 'DMSSeaIceCoupling', DMSSeaIceCoupling)

      call mpas_pool_get_field(DMSSeaIceCoupling, 'iceFluxDMS', iceFluxDMSField)
      call mpas_pool_get_field(DMSSeaIceCoupling, 'iceFluxDMSP', iceFluxDMSPField)
   endif

   ! CFC fields
   if (config_use_CFCTracers) then
      call mpas_pool_get_subpool(forcingPool, 'CFCAuxiliary', CFCAuxiliary)
      call mpas_pool_get_field(CFCAuxiliary, 'windSpeedSquared10mCFC', windSpeedSquared10mField)
   endif

   if ( windStressMeridionalField % isActive ) then
      call mpas_dmpar_exch_halo_field(windStressMeridionalField)
   end if
   if ( windStressZonalField % isActive ) then
      call mpas_dmpar_exch_halo_field(windStressZonalField)
   end if
   if ( latentHeatFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(latentHeatFluxField)
   end if
   if ( sensibleHeatFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(sensibleHeatFluxField)
   end if
   if ( longWaveHeatFluxUpField % isActive ) then
      call mpas_dmpar_exch_halo_field(longWaveHeatFluxUpField)
   end if
   if ( longWaveHeatFluxDownField % isActive ) then
      call mpas_dmpar_exch_halo_field(longWaveHeatFluxDownField)
   end if
   if ( evaporationFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(evaporationFluxField)
   end if
   if ( seaIceHeatFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(seaIceHeatFluxField)
   end if
   if ( icebergHeatFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(icebergHeatFluxField)
   end if
   if ( snowFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(snowFluxField)
   end if
   if ( seaIceFreshWaterFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(seaIceFreshWaterFluxField)
   end if
   if ( icebergFreshWaterFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(icebergFreshWaterFluxField)
   end if
   if ( seaIceSalinityFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(seaIceSalinityFluxField)
   end if
   if ( riverRunoffFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(riverRunoffFluxField)
   end if
   if ( iceRunoffFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(iceRunoffFluxField)
   end if
   if ( shortWaveHeatFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(shortWaveHeatFluxField)
   end if
   if ( rainFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(rainFluxField)
   end if
   if ( atmosphericPressureField % isActive ) then
      call mpas_dmpar_exch_halo_field(atmosphericPressureField)
   end if
   if ( seaIcePressureField % isActive ) then
      call mpas_dmpar_exch_halo_field(seaIcePressureField)
   end if
   if ( iceFractionField % isActive ) then
      call mpas_dmpar_exch_halo_field(iceFractionField)
   end if
   if ( stokesDriftZonalWavenumberField % isActive ) then
      call mpas_dmpar_exch_halo_field(stokesDriftZonalWavenumberField)
   end if
   if ( stokesDriftMeridionalWavenumberField % isActive ) then
      call mpas_dmpar_exch_halo_field(stokesDriftMeridionalWavenumberField)
   end if
   if ( significantWaveHeightField % isActive ) then
      call mpas_dmpar_exch_halo_field(significantWaveHeightField)
   end if
   if ( peakWaveFrequencyField % isActive ) then
      call mpas_dmpar_exch_halo_field(peakWaveFrequencyField)
   end if
   if ( peakWaveDirectionField % isActive ) then
      call mpas_dmpar_exch_halo_field(peakWaveDirectionField)
   end if

   if ( landIceFreshwaterFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(landIceFreshwaterFluxField)
   end if
   if ( landIceHeatFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(landIceHeatFluxField)
   end if
   if ( landIceInterfaceTracersField % isActive ) then
      call mpas_dmpar_exch_halo_field(landIceInterfaceTracersField)
   end if
   if ( landIceFractionField % isActive ) then
      call mpas_dmpar_exch_halo_field(landIceFractionField)
   end if
!   if ( landIcePressureField % isActive ) then
!      call mpas_dmpar_exch_halo_field(landIcePressureField)
!   end if

   if ( windSpeed10mField % isActive ) then
      call mpas_dmpar_exch_halo_field(windSpeed10mField)
   end if

   ! BGC fields
   if (config_use_ecosysTracers) then

      if ( windSpeedSquared10mField % isActive ) then
         call mpas_dmpar_exch_halo_field(windSpeedSquared10mField)
      end if
      if ( atmosphericCO2Field % isActive ) then
         call mpas_dmpar_exch_halo_field(atmosphericCO2Field)
      end if
      if ( atmosphericCO2_ALT_CO2Field % isActive ) then
         call mpas_dmpar_exch_halo_field(atmosphericCO2_ALT_CO2Field)
      end if

      if (config_use_ecosysTracers_river_inputs_from_coupler) then
        if ( riverFluxNO3Field % isActive ) then
           call mpas_dmpar_exch_halo_field(riverFluxNO3Field)
        end if
        if ( riverFluxPO4Field % isActive ) then
           call mpas_dmpar_exch_halo_field(riverFluxPO4Field)
        end if
        if ( riverFluxDONField % isActive ) then
           call mpas_dmpar_exch_halo_field(riverFluxDONField)
        end if
        if ( riverFluxDOPField % isActive ) then
           call mpas_dmpar_exch_halo_field(riverFluxDOPField)
        end if
        if ( riverFluxSiO3Field % isActive ) then
           call mpas_dmpar_exch_halo_field(riverFluxSiO3Field)
        end if
        if ( riverFluxDOCField % isActive ) then
           call mpas_dmpar_exch_halo_field(riverFluxDOCField)
        end if
        if ( riverFluxDICField % isActive ) then
           call mpas_dmpar_exch_halo_field(riverFluxDICField)
        end if
        if ( riverFluxALKField % isActive ) then
           call mpas_dmpar_exch_halo_field(riverFluxALKField)
        end if
        if ( riverFluxFeField % isActive ) then
           call mpas_dmpar_exch_halo_field(riverFluxFeField)
        end if
      endif

      if (config_use_ecosysTracers_sea_ice_coupling) then
         if ( iceFluxPhytoCField % isActive ) then
            call mpas_dmpar_exch_halo_field(iceFluxPhytoCField)
         endif
         if ( iceFluxDICField % isActive ) then
            call mpas_dmpar_exch_halo_field(iceFluxDICField)
         endif
         if ( iceFluxNO3Field % isActive ) then
            call mpas_dmpar_exch_halo_field(iceFluxNO3Field)
         endif
         if ( iceFluxSiO3Field % isActive ) then
            call mpas_dmpar_exch_halo_field(iceFluxSiO3Field)
         endif
         if ( iceFluxNH4Field % isActive ) then
            call mpas_dmpar_exch_halo_field(iceFluxNH4Field)
         endif
         if ( iceFluxDOCrField % isActive ) then
            call mpas_dmpar_exch_halo_field(iceFluxDOCrField)
         endif
         if ( iceFluxFeParticulateField % isActive ) then
            call mpas_dmpar_exch_halo_field(iceFluxFeParticulateField)
         endif
         if ( iceFluxFeDissolvedField % isActive ) then
            call mpas_dmpar_exch_halo_field(iceFluxFeDissolvedField)
         endif
         if ( iceFluxDustField % isActive ) then
            call mpas_dmpar_exch_halo_field(iceFluxDustField)
         endif
         if ( iceFluxDOCField % isActive ) then
            call mpas_dmpar_exch_halo_field(iceFluxDOCField)
         endif
         if ( iceFluxDONField % isActive ) then
            call mpas_dmpar_exch_halo_field(iceFluxDONField)
         endif
      endif
   endif
   if (config_use_DMSTracers .and. config_use_DMSTracers_sea_ice_coupling) then
      if ( iceFluxDMSField % isActive ) then
         call mpas_dmpar_exch_halo_field(iceFluxDMSField)
      endif
      if ( iceFluxDMSPField % isActive ) then
         call mpas_dmpar_exch_halo_field(iceFluxDMSPField)
      endif
   endif

   ! CFC fields
   if (config_use_CFCTracers .and. .not. config_use_ecosysTracers) then
      if ( windSpeedSquared10mField % isActive ) then
         call mpas_dmpar_exch_halo_field(windSpeedSquared10mField)
      end if
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine ocn_import_moab!}}}
!***********************************************************************
!BOP
! !IROUTINE: ocn_export_moab
! !INTERFACE:

 subroutine ocn_export_moab(EClock)   !{{{

! !DESCRIPTION:
!  This routine calls the routines necessary to send MPASO fields to
!  the MOAB driver
!
! !REVISION HISTORY:
!  same as module
   use iMOAB,  only       : iMOAB_SetDoubleTagStorage, iMOAB_WriteMesh
! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock) , intent(inout) :: EClock    ! Input synchronization clock from driver

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   integer :: ent_type, ierr, cur_ocn_stepno
   character(len=100) :: outfile, wopts, localmeshfile, lnum
   character(CXX) :: tagname

   integer :: i, n
   integer, pointer :: nCellsSolve, index_temperatureSurfaceValue, index_salinitySurfaceValue, &
                       index_avgZonalSurfaceVelocity, index_avgMeridionalSurfaceVelocity, &
                       index_avgZonalSSHGradient, index_avgMeridionalSSHGradient

   type (block_type), pointer :: block_ptr

   type (mpas_pool_type), pointer :: meshPool,             &
                                     forcingPool,          &
                                     statePool,            &
                                     tracersPool,          &
                                     ecosysAuxiliary,      &
                                     ecosysSeaIceCoupling, &
                                     DMSSeaIceCoupling,    &
                                     MacroMoleculesSeaIceCoupling

   integer, dimension(:), pointer :: landIceMask

   real (kind=RKIND), dimension(:), pointer :: seaIceEnergy, accumulatedFrazilIceMass, frazilSurfacePressure, &
                                               avgTotalFreshWaterTemperatureFlux, &
                                               avgCO2_gas_flux, DMSFlux, surfaceUpwardCO2Flux, &
                                               avgOceanSurfaceDIC, &
                                               avgOceanSurfaceDON, &
                                               avgOceanSurfaceNO3, &
                                               avgOceanSurfaceSiO3, &
                                               avgOceanSurfaceNH4, &
                                               avgOceanSurfaceDMS, &
                                               avgOceanSurfaceDMSP, &
                                               avgOceanSurfaceDOCr, &
                                               avgOceanSurfaceDOCSemiLabile, &
                                               avgOceanSurfaceFeParticulate, &
                                               avgOceanSurfaceFeDissolved, &
                                               ssh, &
                                               avgLandIceFreshwaterFlux, &
                                               avgRemovedRiverRunoffFlux, &
                                               avgRemovedIceRunoffFlux, &
                                               avgLandIceHeatFlux, &
                                               avgRemovedIceRunoffHeatFlux

   real (kind=RKIND), dimension(:,:), pointer :: avgTracersSurfaceValue, avgSurfaceVelocity, &
                                                 avgSSHGradient, avgOceanSurfacePhytoC, &
                                                 avgOceanSurfaceDOC, layerThickness

   real (kind=RKIND) :: surfaceFreezingTemp

   logical, pointer :: frazilIceActive,          &
                       config_remove_ais_river_runoff, &
                       config_remove_ais_ice_runoff, &
                       config_use_ecosysTracers, &
                       config_use_DMSTracers,    &
                       config_use_MacroMoleculesTracers,  &
                       config_use_ecosysTracers_sea_ice_coupling, &
                       config_use_DMSTracers_sea_ice_coupling,  &
                       config_use_MacroMoleculesTracers_sea_ice_coupling

   character (len=StrKIND), pointer :: config_land_ice_flux_mode

   logical :: keepFrazil


   ! get configure options
   call mpas_pool_get_package(domain % packages, 'frazilIceActive', frazilIceActive)
   call mpas_pool_get_config(domain % configs, 'config_use_ecosysTracers', config_use_ecosysTracers)
   call mpas_pool_get_config(domain % configs, 'config_land_ice_flux_mode', config_land_ice_flux_mode)
   call mpas_pool_get_config(domain % configs, 'config_remove_ais_river_runoff', config_remove_ais_river_runoff)
   call mpas_pool_get_config(domain % configs, 'config_remove_ais_ice_runoff', config_remove_ais_ice_runoff)
   call mpas_pool_get_config(domain % configs, 'config_use_DMSTracers', config_use_DMSTracers)
   call mpas_pool_get_config(domain % configs, 'config_use_MacroMoleculesTracers', config_use_MacroMoleculesTracers)
   call mpas_pool_get_config(domain % configs, 'config_use_ecosysTracers_sea_ice_coupling',  &
                                                config_use_ecosysTracers_sea_ice_coupling)
   call mpas_pool_get_config(domain % configs, 'config_use_DMSTracers_sea_ice_coupling',  &
                                                config_use_DMSTracers_sea_ice_coupling)
   call mpas_pool_get_config(domain % configs, 'config_use_MacroMoleculesTracers_sea_ice_coupling',  &
                                                config_use_MacroMoleculesTracers_sea_ice_coupling)

   n = 0
   block_ptr => domain % blocklist
   do while(associated(block_ptr))
     call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
     call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
     call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)

     call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

     call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

     call mpas_pool_get_dimension(forcingPool, 'index_avgTemperatureSurfaceValue', index_temperatureSurfaceValue)
     call mpas_pool_get_dimension(forcingPool, 'index_avgSalinitySurfaceValue', index_salinitySurfaceValue)
     call mpas_pool_get_dimension(forcingPool, 'index_avgSurfaceVelocityZonal', index_avgZonalSurfaceVelocity)
     call mpas_pool_get_dimension(forcingPool, 'index_avgSurfaceVelocityMeridional', index_avgMeridionalSurfaceVelocity)
     call mpas_pool_get_dimension(forcingPool, 'index_avgSSHGradientZonal', index_avgZonalSSHGradient)
     call mpas_pool_get_dimension(forcingPool, 'index_avgSSHGradientMeridional', index_avgMeridionalSSHGradient)

     call mpas_pool_get_array(statePool, 'ssh', ssh, 1)
     call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)

     call mpas_pool_get_array(forcingPool, 'landIceMask', landIceMask)
     call mpas_pool_get_array(forcingPool, 'avgTracersSurfaceValue', avgTracersSurfaceValue)
     call mpas_pool_get_array(forcingPool, 'avgSurfaceVelocity', avgSurfaceVelocity)
     call mpas_pool_get_array(forcingPool, 'avgSSHGradient', avgSSHGradient)
     call mpas_pool_get_array(forcingPool, 'avgTotalFreshWaterTemperatureFlux', avgTotalFreshWaterTemperatureFlux)

     if ( frazilIceActive ) then
        call mpas_pool_get_array(forcingPool, 'seaIceEnergy', seaIceEnergy)
        call mpas_pool_get_array(forcingPool, 'frazilSurfacePressure', frazilSurfacePressure)
        call mpas_pool_get_array(statePool, 'accumulatedFrazilIceMass', accumulatedFrazilIceMass, 1)
     end if

     ! Cryo fields
     if (trim(config_land_ice_flux_mode) == 'standalone' .or. trim(config_land_ice_flux_mode) == 'data') then
        call mpas_pool_get_array(forcingPool, 'avgLandIceFreshwaterFlux', avgLandIceFreshwaterFlux)
        call mpas_pool_get_array(forcingPool, 'avgLandIceHeatFlux', avgLandIceHeatFlux)
     endif
     if (config_remove_ais_river_runoff) then
        call mpas_pool_get_array(forcingPool, 'avgRemovedRiverRunoffFlux', avgRemovedRiverRunoffFlux)
     endif
     if (config_remove_ais_ice_runoff) then
        call mpas_pool_get_array(forcingPool, 'avgRemovedIceRunoffFlux', avgRemovedIceRunoffFlux)
        call mpas_pool_get_array(forcingPool, 'avgRemovedIceRunoffHeatFlux', avgRemovedIceRunoffHeatFlux)
     endif

     ! BGC fields
     if (config_use_ecosysTracers) then

        call mpas_pool_get_subpool(forcingPool, 'ecosysAuxiliary', ecosysAuxiliary)
        call mpas_pool_get_array(ecosysAuxiliary, 'avgCO2_gas_flux', avgCO2_gas_flux)

     end if

     if (config_use_ecosysTracers .and. config_use_ecosysTracers_sea_ice_coupling) then
        call mpas_pool_get_subpool(forcingPool, 'ecosysSeaIceCoupling', ecosysSeaIceCoupling)

        call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfacePhytoC', avgOceanSurfacePhytoC)
        call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceDIC', avgOceanSurfaceDIC)
        call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceNO3', avgOceanSurfaceNO3)
        call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceSiO3', avgOceanSurfaceSiO3)
        call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceNH4', avgOceanSurfaceNH4)
        call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceDOCr', avgOceanSurfaceDOCr)
        call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceDOCSemiLabile', avgOceanSurfaceDOCSemiLabile)
        call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceFeParticulate', avgOceanSurfaceFeParticulate)
        call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceFeDissolved', avgOceanSurfaceFeDissolved)
     endif
     if (config_use_DMSTracers .and. config_use_DMSTracers_sea_ice_coupling) then
        call mpas_pool_get_subpool(forcingPool, 'DMSSeaIceCoupling', DMSSeaIceCoupling)

        call mpas_pool_get_array(DMSSeaIceCoupling, 'avgOceanSurfaceDMS', avgOceanSurfaceDMS)
        call mpas_pool_get_array(DMSSeaIceCoupling, 'avgOceanSurfaceDMSP', avgOceanSurfaceDMSP)
     endif
     if (config_use_MacroMoleculesTracers .and. config_use_MacroMoleculesTracers_sea_ice_coupling) then
        call mpas_pool_get_subpool(forcingPool, 'MacroMoleculesSeaIceCoupling', MacroMoleculesSeaIceCoupling)

        call mpas_pool_get_array(MacroMoleculesSeaIceCoupling, 'avgOceanSurfaceDOC', avgOceanSurfaceDOC)
        call mpas_pool_get_array(MacroMoleculesSeaIceCoupling, 'avgOceanSurfaceDON', avgOceanSurfaceDON)
     endif
!    call mpas_pool_get_array(forcingPool, 'CO2Flux', CO2Flux)
!    call mpas_pool_get_array(forcingPool, 'DMSFlux', DMSFlux)
!    call mpas_pool_get_array(forcingPool, 'surfaceUpwardCO2Flux', surfaceUpwardCO2Flux)

! replace 'o2x_o % rAttr(' with 'o2x_om(n, ' and ', n)' with ')' 
     do i = 1, nCellsSolve
       n = n + 1

       o2x_om(n, index_o2x_So_t) = avgTracersSurfaceValue(index_temperatureSurfaceValue, i)
       o2x_om(n, index_o2x_So_s) = avgTracersSurfaceValue(index_salinitySurfaceValue, i)
       o2x_om(n, index_o2x_So_u) = avgSurfaceVelocity(index_avgZonalSurfaceVelocity, i)
       o2x_om(n, index_o2x_So_v) = avgSurfaceVelocity(index_avgMeridionalSurfaceVelocity, i)

       o2x_om(n, index_o2x_So_ssh)  = ssh(i)
       o2x_om(n, index_o2x_So_dhdx) = avgSSHGradient(index_avgZonalSSHGradient, i)
       o2x_om(n, index_o2x_So_dhdy) = avgSSHGradient(index_avgMeridionalSSHGradient, i)

       o2x_om(n, index_o2x_Faoo_h2otemp) = avgTotalFreshWaterTemperatureFlux(i) * rho_sw * cp_sw

       ! Cryo fields
       if (trim(config_land_ice_flux_mode) == 'standalone' .or. trim(config_land_ice_flux_mode) == 'data') then
          o2x_om(n, index_o2x_Foxo_ismw)  = avgLandIceFreshwaterFlux(i)
          o2x_om(n, index_o2x_Foxo_ismh)  = avgLandIceHeatFlux(i)
       endif
       if (config_remove_ais_river_runoff) then
          o2x_om(n, index_o2x_Foxo_rrofl) = avgRemovedRiverRunoffFlux(i)
       endif
       if (config_remove_ais_ice_runoff) then
          o2x_om(n, index_o2x_Foxo_rrofi) = avgRemovedIceRunoffFlux(i)
          o2x_om(n, index_o2x_Foxo_rrofih) = avgRemovedIceRunoffHeatFlux(i)
       endif

       if ( frazilIceActive ) then
          ! negative when frazil ice can be melted
          keepFrazil = .true.
          if ( associated(landIceMask) ) then
             if ( landIceMask(i) == 1 ) then
                keepFrazil = .false.
             end if
          end if

          if ( keepFrazil ) then

             ! Calculate energy associated with frazil mass transfer to sea ice if frazil has accumulated
             if ( accumulatedFrazilIceMass(i) > 0.0_RKIND ) then

              seaIceEnergy(i) = accumulatedFrazilIceMass(i) * config_frazil_heat_of_fusion

             ! Otherwise calculate the melt potential where avgTracersSurfaceValue represents only the
             ! top layer of the ocean
             else

              surfaceFreezingTemp = ocn_freezing_temperature(salinity=avgTracersSurfaceValue(index_salinitySurfaceValue, i), &
                 pressure=0.0_RKIND,  inLandIceCavity=.false.)

              seaIceEnergy(i) = min(rho_sw*cp_sw*layerThickness(1, i)*( surfaceFreezingTemp + T0_Kelvin &
                              - avgTracersSurfaceValue(index_temperatureSurfaceValue, i) ), 0.0_RKIND )

             end if

             o2x_om(n, index_o2x_Fioo_q)  = seaIceEnergy(i) / ocn_cpl_dt
             o2x_om(n, index_o2x_Fioo_frazil) = accumulatedFrazilIceMass(i) / ocn_cpl_dt

          else

             o2x_om(n, index_o2x_Fioo_q)  = 0.0_RKIND
             o2x_om(n, index_o2x_Fioo_frazil) = 0.0_RKIND
             if (trim(config_land_ice_flux_mode) == 'standalone' .or. trim(config_land_ice_flux_mode) == 'data') then
                o2x_om(n, index_o2x_Foxo_q_li) = accumulatedFrazilIceMass(i) * config_frazil_heat_of_fusion / ocn_cpl_dt
                o2x_om(n, index_o2x_Foxo_frazil_li) = accumulatedFrazilIceMass(i) / ocn_cpl_dt
             endif

          end if

          ! Reset SeaIce Energy and Accumulated Frazil Ice
          seaIceEnergy(i) = 0.0_RKIND
          accumulatedFrazilIceMass(i) = 0.0_RKIND
          frazilSurfacePressure(i) = 0.0_RKIND
       end if

       ! BGC fields
       if (config_use_ecosysTracers .and. index_o2x_Faoo_fco2_ocn /= 0) then
          ! convert from mmolC/m2/s to kg CO2/m2/s
          o2x_om(n, index_o2x_Faoo_fco2_ocn) = avgCO2_gas_flux(i)*44.e-6_RKIND
       endif
       if (config_use_ecosysTracers .and. config_use_ecosysTracers_sea_ice_coupling) then
          o2x_om(n, index_o2x_So_algae1) = max(0.0_RKIND,avgOceanSurfacePhytoC(1,i))
          o2x_om(n, index_o2x_So_algae2) = max(0.0_RKIND,avgOceanSurfacePhytoC(2,i))
          o2x_om(n, index_o2x_So_algae3) = max(0.0_RKIND,avgOceanSurfacePhytoC(3,i))
          o2x_om(n, index_o2x_So_dic1) = max(0.0_RKIND,avgOceanSurfaceDIC(i))
          o2x_om(n, index_o2x_So_doc1) = max(0.0_RKIND,avgOceanSurfaceDOCSemiLabile(i))
          o2x_om(n, index_o2x_So_doc2) = max(0.0_RKIND,avgOceanSurfaceDOCSemiLabile(i))
          o2x_om(n, index_o2x_So_doc3) = max(0.0_RKIND,avgOceanSurfaceDOCSemiLabile(i))
          o2x_om(n, index_o2x_So_don1) = 0.0_RKIND
          o2x_om(n, index_o2x_So_no3) = max(0.0_RKIND,avgOceanSurfaceNO3(i))
          o2x_om(n, index_o2x_So_sio3) = max(0.0_RKIND,avgOceanSurfaceSiO3(i))
          o2x_om(n, index_o2x_So_nh4) = max(0.0_RKIND,avgOceanSurfaceNH4(i))
          o2x_om(n, index_o2x_So_docr) = max(0.0_RKIND,avgOceanSurfaceDOCr(i))
          o2x_om(n, index_o2x_So_fep1) = max(0.0_RKIND,avgOceanSurfaceFeParticulate(i))
          o2x_om(n, index_o2x_So_fed1) = max(0.0_RKIND,avgOceanSurfaceFeDissolved(i))
       endif
       if (config_use_DMSTracers .and. config_use_DMSTracers_sea_ice_coupling) then
          o2x_om(n, index_o2x_So_dms) = max(0.0_RKIND,avgOceanSurfaceDMS(i))
          o2x_om(n, index_o2x_So_dmsp) = max(0.0_RKIND,avgOceanSurfaceDMSP(i))
       endif
       if (config_use_MacroMoleculesTracers .and. config_use_MacroMoleculesTracers_sea_ice_coupling) then
          o2x_om(n, index_o2x_So_doc1) = max(0.0_RKIND,avgOceanSurfaceDOC(1,i))
          o2x_om(n, index_o2x_So_doc2) = max(0.0_RKIND,avgOceanSurfaceDOC(2,i))
          o2x_om(n, index_o2x_So_doc3) = max(0.0_RKIND,avgOceanSurfaceDOC(3,i))
          o2x_om(n, index_o2x_So_don1) = 0.0_RKIND
       endif
!      o2x_om(n, index_o2x_Faoo_fco2_ocn) = CO2Flux(i)
!      o2x_om(n, index_o2x_Faoo_fdms_ocn) = DMSFlux(i)
!      o2x_om(n, index_o2x_Faoo_fco2_ocn) = surfaceUpwardCO2Flux(i)

!JW       o2x_om(n, index_o2x_So_blt) = landIceBoundaryLayerTemperature(i)
!JW       o2x_om(n, index_o2x_So_bls) = landIceBoundaryLayerSalinity(i)
!JW       o2x_om(n, index_o2x_So_htv) = landIceHeatTransferVelocity(i)
!JW       o2x_om(n, index_o2x_So_stv) = landIceSaltTransferVelocity(i)

       if ( trim(config_land_ice_flux_mode) .eq. 'standalone' .or. &
            trim(config_land_ice_flux_mode) .eq. 'coupled'  ) then
          o2x_om(n, index_o2x_So_blt) = landIceBoundaryLayerTracers(indexBLT,i)
          o2x_om(n, index_o2x_So_bls) = landIceBoundaryLayerTracers(indexBLS,i)
          o2x_om(n, index_o2x_So_htv) = landIceTracerTransferVelocities(indexHeatTrans,i)
          o2x_om(n, index_o2x_So_stv) = landIceTracerTransferVelocities(indexSaltTrans,i)
          o2x_om(n, index_o2x_So_rhoeff) = 0.0_RKIND
       endif

       !Fyke: test
       !write(stderrUnit,*) 'n=',n
       !write(stderrUnit,*) 'o2x_om(n, index_o2x_So_blt)=',o2x_om(n, index_o2x_So_blt)
       !write(stderrUnit,*) 'o2x_om(n, index_o2x_So_bls)=',o2x_om(n, index_o2x_So_bls)
       !write(stderrUnit,*) 'o2x_om(n, index_o2x_So_htv)=',o2x_om(n, index_o2x_So_htv)
       !write(stderrUnit,*) 'o2x_om(n, index_o2x_So_stv)=',o2x_om(n, index_o2x_So_stv)
       !write(stderrUnit,*) 'o2x_om(n, index_o2x_So_rhoeff)=',o2x_om(n, index_o2x_So_rhoeff)
       !o2x_om(n, index_o2x_So_blt) = 0._r8
       !o2x_om(n, index_o2x_So_bls) = 34.5_r8
       !o2x_om(n, index_o2x_So_htv) = 1.e-4_r8
       !o2x_om(n, index_o2x_So_stv) = 3.e-6_r8
       !o2x_om(n, index_o2x_So_rhoeff) = 1000._r8*9.81_r8*918._r8 !lithostatic pressure of 1km of ice

     end do

     block_ptr => block_ptr % next
   end do
      ent_type = 1 ! cells
      !  set all tags in one method
      tagname = trim(seq_flds_o2x_fields)//C_NULL_CHAR
      ierr = iMOAB_SetDoubleTagStorage ( MPOID, tagname, totalmbls , ent_type, o2x_om(1, 1) )
      if ( ierr /= 0 ) then
         write(ocnLogUnit,*) 'Fail to set o2x MOAB fields '
      endif
   !-----------------------------------------------------------------------
   !EOC
      call seq_timemgr_EClockGetData( EClock, stepno=cur_ocn_stepno )
#ifdef MOABDEBUG
      write(lnum,"(I0.2)")cur_ocn_stepno
      outfile = 'ocn_export_'//trim(lnum)//'.h5m'//C_NULL_CHAR
      wopts   = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
      ierr = iMOAB_WriteMesh(MPOID, outfile, wopts)
#endif

!-----------------------------------------------------------------------
!EOC

 end subroutine ocn_export_moab!}}}
#endif
end module ocn_comp_mct

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
