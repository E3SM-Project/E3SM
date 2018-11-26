module ice_comp_mct

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
! !MODULE: ice_comp_mct
! !INTERFACE:

! !DESCRIPTION:
!  This is the main driver for the Model for Predication Across Scales SEAICE Model (MPASSI).
!
! !USES:
   ! CPL modules
   use mct_mod
   use esmf
   use seq_flds_mod
   use seq_cdata_mod
   use seq_infodata_mod
   use seq_timemgr_mod
   use seq_comm_mct,      only : seq_comm_suffix, seq_comm_inst, seq_comm_name
   use shr_file_mod 
   use shr_cal_mod,       only : shr_cal_date2ymd
   use shr_sys_mod
   use shr_taskmap_mod,   only : shr_taskmap_write
   use shr_pio_mod
   use perf_mod

   ! MPASSI driver modules
   use mpassi_cpl_indices
   use mpassi_mct_vars

   ! MPAS framework modules
   use mpas_framework
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_stream_manager
   use mpas_kind_types
   use mpas_io_units
   use mpas_timekeeping
   use mpas_bootstrapping
   use mpas_dmpar
   use mpas_log

   use iso_c_binding, only : c_char, c_loc, c_ptr, c_int
   use mpas_c_interfacing, only : mpas_f_to_c_string, mpas_c_to_f_string

   ! MPASSI modules
   use seaice_analysis_driver
   use seaice_column, only: seaice_column_reinitialize_fluxes, &
                            seaice_column_coupling_prep
   use seaice_constants
   use seaice_core
   use seaice_core_interface
   use seaice_forcing
   use seaice_initialize, only: seaice_init_post_clock_advance
   use seaice_mesh, only: seaice_latlon_vector_rotation_backward
   use seaice_time_integration
   use seaice_error, only: seaice_check_critical_error

   use ice_constants_colpkg, only: &
       Tffresh

!
! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  public :: ice_init_mct
  public :: ice_run_mct
  public :: ice_final_mct
  SAVE
  private                              ! By default make data private

!
! ! PUBLIC DATA:
!
! !REVISION HISTORY:
! Author: Jon Wolfe
!
!EOP
! !PRIVATE MODULE FUNCTIONS:
  private :: ice_export_mct
  private :: ice_import_mct
  private :: ice_SetGSMap_mct
  private :: ice_domain_mct
!
! !PRIVATE MODULE VARIABLES

   integer, private :: my_task

   integer  :: nsend, nrecv

   character(len=StrKIND) :: runtype, coupleTimeStamp

   type(seq_infodata_type), pointer :: infodata   
   type (iosystem_desc_t), pointer :: io_system 

   integer :: iceLogUnit ! unit number for ice log

   !! MPASSI Datatypes
   type (core_type), pointer :: corelist => null()
   type (domain_type), pointer :: domain
   integer :: itimestep, ice_cpl_dt

!=======================================================================

contains

!***********************************************************************
!BOP
!
! !IROUTINE: ice_init_mct
!
! !INTERFACE:
  subroutine ice_init_mct( EClock, cdata_i, x2i_i, i2x_i, NLFilename )!{{{
!
! !DESCRIPTION:
! Initialize MPASSI
!
! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock), intent(inout) :: EClock
    type(seq_cdata), intent(inout) :: cdata_i
    type(mct_aVect), intent(inout) :: x2i_i, i2x_i
    character(len=*), optional, intent(in) :: NLFilename ! Namelist filename
!
! !REVISION HISTORY:
! Author: Jon Wolfe
!EOP
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    integer :: ICEID !< task number
    integer :: mpicom_i, lsize, start_ymd, start_tod, start_year, start_day,   &
       start_month, start_hour, iyear, mpas_i_cpl_dt, pio_iotype
    integer :: shrloglev, shrlogunit
    integer :: icePossibleErrUnit !< unit number to reserve for if a err log needs to be opened

    type(mct_gsMap), pointer :: gsMap_i

    type(mct_gGrid), pointer :: dom_i

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
    type (block_type), pointer :: block

    type (MPAS_Pool_Type), pointer :: shortwave

    logical :: exists

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

    logical, pointer :: tempLogicalConfig
    character(len=StrKIND), pointer :: tempCharConfig
    real (kind=RKIND), pointer :: tempRealConfig

    real(kind=RKIND), pointer :: &
         dayOfNextShortwaveCalculation ! needed for CESM like coupled simulations

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
!   setup mpassi data structures
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

    call seq_cdata_setptrs(cdata_i, ID=ICEID, mpicom=mpicom_i, &
         gsMap=gsMap_i, dom=dom_i, infodata=infodata)

    MPASSI_MCT_ICEID =  ICEID  ! task number
    MPASSI_MCT_gsMap_i => gsMap_i
    MPASSI_MCT_dom_i   => dom_i

    call MPI_comm_rank(mpicom_i,iam,ierr)
#if (defined _MEMTRACE)
    if(iam == 0) then
        lbnum=1
        call memmon_dump_fort('memmon.out','ice_init_mct:start::',lbnum) 
    endif
#endif

!-----------------------------------------------------------------------
!
!   first initializaiton phase of mpassi
!   initialize mpassi because grid information is needed for
!   creation of GSMap_ice.
!   call mpassi initialization routines
!
!-----------------------------------------------------------------------

    call t_startf('mpassi_init')

    io_system => shr_pio_getiosys(iceid)

    pio_iotype = shr_pio_getiotype(iceid)
!    call MPAS_io_set_iotype(domain % iocontext, pio_iotype)
    call MPAS_io_unset_iotype(domain % iocontext)

    ! ----------------
    ! Set up log file information
    ! ----------------
    inst_suffix = seq_comm_suffix(iceID) ! a suffix to append to log file name
    iceLogUnit = shr_file_getUnit() ! reserve unit number for log unit
    icePossibleErrUnit = shr_file_getUnit() ! reserve unit number for possible error log file

    ! Note: In following code, a file ice_modelio.nml is queried to determine the log file name to use.
    ! This file is generated when a case is created.  The default file name in there is currently ice.log.DATE-TIME.
    if (iam==0) then
       inquire(file='ice_modelio.nml'//trim(inst_suffix),exist=exists)
       call shr_file_setio('ice_modelio.nml'//trim(inst_suffix), iceLogUnit)
    endif

    ! Store shr log unit and level so we can reassign them later when ICE control is complete
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    ! Set the shr log unit to be the iceLogUnit
    !  (I believe this is so that shr routines called from the ICE driver write any messages to the ICE log file.)
    call shr_file_setLogUnit(iceLogUnit)

    ! Write to iceLogUnit here, because the log module is not open yet.
    if (iam==0) write(iceLogUnit,'(a,i6)') '=== Beginning ice_init_mct: rank=',iam

    call mpas_framework_init_phase1(domain % dminfo, mpicom_i)

    ! Identify SMP nodes and process/SMP mapping for this instance
    ! (Assume that processor names are SMP node names on SMP clusters.)
    call MPI_comm_size(mpicom_i,npes,ierr)
    write(c_npes,'(i8)') npes

    inst_index = seq_comm_inst(iceID)
    write(c_inst_index,'(i8)') inst_index

    if (iam==0) then
       write(iceLogUnit,100) trim(adjustl(c_npes)), trim(adjustl(c_inst_index))
100    format(/,a,' pes participating in computation of MPAS-SEAICE instance #',a)
       call flush(iceLogUnit)
    endif

    call t_startf("shr_taskmap_write")
    call shr_taskmap_write(iceLogUnit, mpicom_i, 'ICE #'//trim(adjustl(c_inst_index)))
    call t_stopf("shr_taskmap_write")

    ! Setup function pointers for MPASSI core and domain types
    call seaice_setup_core(corelist)
    call seaice_setup_domain(domain)

    ! ===========
    ! Initialize log manager
    call mpas_log_init(domain % logInfo, domain, unitNumbers=(/iceLogUnit, icePossibleErrUnit/), err=ierr)
    if ( ierr /= 0 ) then
       write(iceLogUnit,*) 'ERROR: log init failed for core ' // trim(domain % core % coreName)
       call mpas_dmpar_abort(domain % dminfo)
    end if

    ! Set core specific options here
    ! Disable output from all but the master task for ACME!
    ! (This overrides the default set by mpas_log_init based on MPAS_DEBUG setting.)
    if (iam /= 0) then
       domain % logInfo % outputLog % isActive = .false.
    endif

    ! After core has had a chance to modify log defaults, open the output log
    call mpas_log_open(err=ierr)
    if ( ierr /= 0 ) then
       write(iceLogUnit,*) 'ERROR: log open failed for core ' // trim(domain % core % coreName)
       call mpas_dmpar_abort(domain % dminfo)
    end if
    ! ===========


    ! ----------
    ! Process namelist and streams files
    ! ----------
    ! Override the names of the stream and namelist files
    domain % namelist_filename = 'mpassi_in'
    domain % streams_filename = 'streams.seaice'

    ! Setup namelist variables, and read the namelist
    ierr = domain % core % setup_namelist(domain % configs, domain % namelist_filename, domain % dminfo)
    if ( ierr /= 0 ) then
       call mpas_log_write('Namelist setup failed for core ' // trim(domain % core % coreName), MPAS_LOG_CRIT)
    end if

    call mpas_framework_init_phase2(domain, io_system)

    ! Define package variables
    ierr = domain % core % define_packages(domain % packages)
    if ( ierr /= 0 ) then
       call mpas_log_write('Package definition failed for core ' // trim(domain % core % coreName), MPAS_LOG_CRIT)
    end if

    ! Setup packages (i.e. determine if they should be on or off)
    ierr = domain % core % setup_packages(domain % configs, domain % packages, domain % iocontext)
    if ( ierr /= 0 ) then
       call mpas_log_write('Package setup failed for core ' // trim(domain % core % coreName), MPAS_LOG_CRIT)
    end if

    ! Setup decompositions available for dimensions
    ierr = domain % core % setup_decompositions(domain % decompositions)
    if ( ierr /= 0 ) then
       call mpas_log_write('Decomposition setup failed for core ' // trim(domain % core % coreName), MPAS_LOG_CRIT)
    end if


    ! ----------
    ! Override namelist options based on start type
    ! ----------
    call seq_infodata_GetData( infodata, start_type=starttype)
    if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
       runtype = "initial"
    else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
       runtype = "continue"
    else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
       runtype = "branch"
    else
       call mpas_log_write('ice_comp_mct ERROR: infodata provided unknown starttype: ' // trim(starttype), MPAS_LOG_CRIT)
    end if
    if (runtype == "initial") then ! Start up run

        ! Turn off restart
        call mpas_pool_get_config(domain % configs, "config_do_restart", tempLogicalConfig)
        tempLogicalConfig = .false.

        ! Setup start time. Will be over written later when clocks are synchronized
        call mpas_pool_get_config(domain % configs, "config_start_time", tempCharConfig)
        tempCharConfig = trim(tempCharConfig) // "0:00:00"

        ! Setup run duration. Will be ignored in coupled run, since coupler defines how long the run is.
        call mpas_pool_get_config(domain % configs, "config_run_duration", tempCharConfig)
        tempCharConfig = "0001-00-00_00:00:00"

    else if (runtype == "continue" .or. runtype == "branch") then ! Restart run or branch run

        ! Turn on restart
        call mpas_pool_get_config(domain % configs, "config_do_restart", tempLogicalConfig)
        tempLogicalConfig = .true.

        ! Set start time to be read from file
        call mpas_pool_get_config(domain % configs, "config_start_time", tempCharConfig)
        tempCharConfig = "file"

        ! Setup run duration. Will be ignored in coupled run, since coupler defines how long the run is.
        call mpas_pool_get_config(domain % configs, "config_run_duration", tempCharConfig)
        tempCharConfig = "0001-00-00_00:00:00"

    end if

    ! Setup MPASSI simulation clock
    ierr = domain % core % setup_clock(domain % clock, domain % configs)
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
    ierr = domain % core % get_mesh_stream(domain % configs, mesh_stream)
    if ( ierr /= 0 ) then
       call mpas_log_write('Failed to find mesh stream for core ' // trim(domain % core % coreName), MPAS_LOG_CRIT)
    end if


    call mpas_f_to_c_string(domain % streams_filename, c_filename)
    call mpas_f_to_c_string(mesh_stream, c_mesh_stream)
    c_comm = domain % dminfo % comm
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

    start_time = mpas_get_clock_time(domain % clock, MPAS_START_TIME, ierr)
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
    call mpas_bootstrap_framework_phase1(domain, mesh_filename, mesh_iotype)

    !
    ! Set up run-time streams
    !
    call MPAS_stream_mgr_init(domain % streamManager, domain % ioContext, domain % clock, domain % blocklist % allFields, &
                              domain % packages, domain % blocklist % allStructs)

    call add_stream_attributes(domain)

    ! Setup all immutable streams for the core
    ierr = domain % core % setup_immutable_streams(domain % streamManager)
    if ( ierr /= 0 ) then
       call mpas_log_write('Immutable streams setup failed for core ' // trim(domain % core % coreName), MPAS_LOG_CRIT)
    end if

    ! Parse / read all streams configuration
    mgr_p = c_loc(domain % streamManager)
    call xml_stream_parser(c_filename, mgr_p, c_comm, c_ierr)
    if (c_ierr /= 0) then
       call mpas_log_write('xml_stream_parser failed.', MPAS_LOG_CRIT)
    end if

    my_task = domain % dminfo % my_proc_id

    !
    ! Finalize the setup of blocks and fields
    !
    call mpas_bootstrap_framework_phase2(domain)

    ! Determine coupling type
    call seq_infodata_GetData(infodata, cpl_seq_option=cpl_seq_option)

    ! Determine time of next atmospheric shortwave calculation
    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "shortwave", shortwave)
       call MPAS_pool_get_array(shortwave, "dayOfNextShortwaveCalculation", dayOfNextShortwaveCalculation)
       call seq_infodata_GetData(infodata, nextsw_cday=dayOfNextShortwaveCalculation )

       ! Set dayOfNextShortwaveCalculation to -1 for continue and branch runs.
       if (trim(runtype) /= 'initial') dayOfNextShortwaveCalculation = -1
       call mpas_log_write('dayOfNextShortwaveCalculation = $r', realArgs=(/dayOfNextShortwaveCalculation/))

       block => block % next
    end do

    ! Initialize the MPASSI core
    ierr = domain % core % core_init(domain, timeStamp)
    if ( ierr /= 0 ) then
       call mpas_log_write('Core init failed for core ' // trim(domain % core % coreName), MPAS_LOG_CRIT)
    end if

!-----------------------------------------------------------------------
!
!   initialize time-stamp information
!
!-----------------------------------------------------------------------
    call t_stopf ('mpassi_init')

!-----------------------------------------------------------------------
!
!   check for consistency of mpassi and sync clock initial time
!
!-----------------------------------------------------------------------

    if (runtype == 'initial') then
       call seq_timemgr_EClockGetData(EClock, ECurrTime=currTime % t)
       call mpas_set_clock_time(domain % clock, currTime, MPAS_START_TIME, ierr)
       call mpas_set_clock_time(domain % clock, currTime, MPAS_NOW, ierr)

       currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
       call mpas_get_time(curr_time=currTime, YYYY=iyear0, MM=imonth0, ierr=ierr)
       call seq_timemgr_EClockGetData(EClock, start_ymd=start_ymd, start_tod=start_tod)
       call shr_cal_date2ymd(start_ymd,start_year,start_month,start_day)

       if (iyear0 /= start_year) then
          call mpas_log_write('iyear0=$i does not match start_year=$i', MPAS_LOG_CRIT, intArgs=(/iyear0, start_year/))
       end if
       if (imonth0 /= start_month) then
          call mpas_log_write('imonth0=$i does not match start_month=$i', MPAS_LOG_CRIT, intArgs=(/imonth0, start_month/))
       end if
    else if (runtype == 'continue' .or. runtype == 'branch') then
       call seq_timemgr_EClockGetData(EClock, ECurrTime=currTime % t)
       call mpas_set_clock_time(domain % clock, currTime, MPAS_START_TIME, ierr)
       call mpas_set_clock_time(domain % clock, currTime, MPAS_NOW, ierr)
    end if

!-----------------------------------------------------------------------
!
!   initialize MCT attribute vectors and indices
!
!-----------------------------------------------------------------------

    call t_startf ('mpassi_mct_init')

    call mpassi_cpl_indices_set()

    call ice_SetGSMap_mct( mpicom_i, ICEID, GSMap_i )
    lsize = mct_gsMap_lsize(gsMap_i, mpicom_i)

    ! Initialize mct ice domain (needs ice initialization info)
    call ice_domain_mct( lsize, gsMap_i, dom_i )
    
    ! Inialize mct attribute vectors
    
    call mct_aVect_init(x2i_i, rList=seq_flds_x2i_fields, lsize=lsize)
    call mct_aVect_zero(x2i_i)
    
    call mct_aVect_init(i2x_i, rList=seq_flds_i2x_fields, lsize=lsize) 
    call mct_aVect_zero(i2x_i)
    
    nsend = mct_avect_nRattr(i2x_i)
    nrecv = mct_avect_nRattr(x2i_i)

!-----------------------------------------------------------------------
!
!   initialize necessary  coupling info
!
!-----------------------------------------------------------------------

    ! set ice_cpl_dt (coupling interval in integer seconds) from the EClock
    call seq_timemgr_EClockGetData(EClock, dtime=ice_cpl_dt)
    ! convert this dt to a time interval
    call mpas_set_timeInterval(alarmTimeStep, S=ice_cpl_dt, ierr=ierr)
    call mpas_get_timeInterval(alarmTimeStep, timeString=coupleTimeStamp, ierr=ierr)
    call mpas_log_write( 'Applying ICE coupling dt (s) of: $r', intArgs=(/ice_cpl_dt/))

    ! Verify the mpas time step fits into a coupling interval
    call mpas_pool_get_config(domain % configs, 'config_dt', tempRealConfig)
    call mpas_set_timeInterval(denInterval, dt=tempRealConfig, ierr=ierr)
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
       call mpas_log_write('ICE Model time step is: $r', MPAS_LOG_ERR, realArgs=(/tempRealConfig/))

       call mpas_log_write('The model time step cannot be longer then the coupling interval', MPAS_LOG_ERR)
       call mpas_log_write('Model is not properly configured for coupling interval.', MPAS_LOG_CRIT)
    end if

    if ( remInterval > zeroInterval ) then
       ierr = 1
    end if

    ierr_local = ierr
    call mpas_dmpar_max_int(domain % dminfo, ierr_local, ierr)

    if ( ierr == 1 ) then
       call mpas_log_write('Coupling interval is: ' // trim(coupleTimeStamp), MPAS_LOG_ERR)
       call mpas_log_write('        ICE Model time step is: $r', MPAS_LOG_ERR, realArgs=(/tempRealConfig/))
       call mpas_log_write('        These are not synchronized, so time steps will not match to coupling interval boundaries.', MPAS_LOG_ERR)
       call mpas_log_write('        Please reconfigure either the coupling interval or the time step.', MPAS_LOG_ERR)
       call mpas_log_write('Model is not properly configured for coupling interval.', MPAS_LOG_CRIT)
    end if

    ! set coupling alarm
    alarmStartTime = currTime
    call mpas_add_clock_alarm(domain % clock, coupleAlarmID, alarmStartTime, alarmTimeStep, ierr=ierr)
    call mpas_reset_clock_alarm(domain % clock, coupleAlarmID, ierr=ierr)

    ! Coupling prep
    call seaice_column_coupling_prep(domain)


!-----------------------------------------------------------------------
!
!   send intial state to driver
!
!-----------------------------------------------------------------------

    call ice_export_mct(i2x_i, errorCode)  
    if (errorCode /= 0) then
       call mpas_log_write('Error in ice_export_mct', MPAS_LOG_CRIT)
    endif

    call t_stopf ('mpassi_mct_init')

    call seq_infodata_PutData( infodata, ice_prognostic=.true.)

!-----------------------------------------------------------------------
!
!   get intial state from driver
!
!-----------------------------------------------------------------------

    call ice_import_mct(x2i_i, errorCode)  
    if (errorCode /= 0) then
       call mpas_log_write('Error in ice_import_mct', MPAS_LOG_CRIT)
    endif

    currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
    call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
    call mpas_log_write('Initial time ' // trim(timeStamp))

    itimestep = 0

    call seaice_analysis_compute_startup(domain, ierr) 

    ! Reset all output alarms, to prevent intial time step from writing any output, unless it's ringing.
    call mpas_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)
    call mpas_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_INPUT, ierr=ierr)

!----------------------------------------------------------------------------
!
!   Reset shr logging to original values
!
!----------------------------------------------------------------------------

    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

#if defined (_MEMTRACE)
    if(iam  == 0) then
        lbnum=1
        call memmon_dump_fort('memmon.out','ice_init_mct:end::',lbnum) 
        call memmon_reset_addr()
    endif
#endif
    call mpas_log_write('=== Completed ice_init_mct ===', flushNow=.true.)

!-----------------------------------------------------------------------
!EOC

 end subroutine ice_init_mct!}}}

!***********************************************************************
!BOP
!
! !IROUTINE: ice_run_mct
!
! !INTERFACE:
  subroutine ice_run_mct( EClock, cdata_i, x2i_i, i2x_i)!{{{

!
! !DESCRIPTION:
! Run MPASSI for one coupling interval
!
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_Clock)            , intent(inout) :: EClock
      type(seq_cdata)             , intent(inout) :: cdata_i
      type(mct_aVect)             , intent(inout) :: x2i_i
      type(mct_aVect)             , intent(inout) :: i2x_i

!! !REVISION HISTORY:
!! Author: Jon Wolfe
!!EOP
!!-----------------------------------------------------------------------
!!
!!  local variables
!!
!!-----------------------------------------------------------------------
! Variables related to ACME coupler
      integer :: ymd, tod, ihour, iminute, isecond
      integer :: iyear, imonth, iday, curr_ymd, curr_tod
      integer :: shrloglev, shrlogunit

! Variable related to MPASSI
      type (block_type), pointer :: block
      type (MPAS_Pool_type), pointer :: &
            shortwave

      real (kind=RKIND) :: current_wallclock_time
      type (MPAS_Time_Type) :: currTime
      character(len=StrKIND) :: timeStamp, streamName, WCstring
      type (MPAS_timeInterval_type) :: timeStep
      integer :: ierr, streamDirection, iam
      logical :: streamActive, debugOn
      logical, pointer :: config_write_output_on_startup
      logical, save :: first=.true.
      character (len=StrKIND), pointer :: config_restart_timestamp_name
      real(kind=RKIND), pointer :: &
            dayOfNextShortwaveCalculation ! needed for CESM like coupled simulations

      iam = domain % dminfo % my_proc_id

      debugOn = .false.
#ifdef MPAS_DEBUG
      debugOn = .true.
#endif
 
      ! Set MPAS Log module instance
      mpas_log_info => domain % logInfo
      if (debugOn) call mpas_log_write("=== Beginning ice_run_mct ===")

      call mpas_pool_get_config(domain % configs, 'config_restart_timestamp_name', config_restart_timestamp_name)

      ! Setup log information.
      call shr_file_getLogUnit (shrlogunit)
      call shr_file_getLogLevel(shrloglev)
      call shr_file_setLogUnit (iceLogUnit)

      ! reinitialize fluxes
      call seaice_column_reinitialize_fluxes(domain)

      ! Import state from coupler
      call ice_import_mct(x2i_i, ierr)

      ! Post coupling calls
      block => domain % blocklist
      do while (associated(block))

         ! Determine time of next atmospheric shortwave calculation
         call MPAS_pool_get_subpool(block % structs, "shortwave", shortwave)
         call MPAS_pool_get_array(shortwave, "dayOfNextShortwaveCalculation", dayOfNextShortwaveCalculation)
         call seq_infodata_GetData(infodata, nextsw_cday=dayOfNextShortwaveCalculation )

         ! perform post coupling operations
         call post_atmospheric_coupling(block)
         call post_oceanic_coupling(block)

         block => block % next
      end do

      ! reset coupler alarm before we start
      call mpas_reset_clock_alarm(domain % clock, coupleAlarmID, ierr=ierr)

     ! Get current time
     currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)

     timeStep = mpas_get_clock_timestep(domain % clock, ierr=ierr)
     call mpas_reset_clock_alarm(domain % clock, coupleAlarmID, ierr=ierr)

     itimestep = 0  ! We may want to initialize this in init and make it a module, save variable

     ! During integration, time level 1 stores the model state at the beginning of the
     !   time step, and time level 2 stores the state advanced dt in time by timestep(...)
     do while (.not. mpas_is_alarm_ringing(domain % clock, coupleAlarmID, ierr=ierr))
        call mpas_stream_mgr_read(domain % streamManager, ierr=ierr)
        call mpas_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_INPUT, ierr=ierr)

        itimestep = itimestep + 1
        call mpas_advance_clock(domain % clock)

        ! final initialization after clock advance 
        if (first) then
           call seaice_forcing_get(domain % streamManager, domain, domain % clock, .true.)

           call seaice_init_post_clock_advance(domain, domain % clock)
           first = .false.
        else
           ! forcing
           call seaice_forcing_get(domain % streamManager, domain, domain % clock, .false.)
        endif

        currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
        call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
        ! write time stamp on every step
        call mpas_dmpar_get_time(current_wallclock_time)
        write (WCstring,'(F18.3)') current_wallclock_time
        call mpas_log_write(trim(timeStamp) // '  WC time:' // WCstring)

        ! pre-timestep analysis computation
        if (debugOn) call mpas_log_write('   Starting analysis precompute', masterOnly=.true.)
        call seaice_analysis_precompute(domain, ierr)
        if (debugOn) call mpas_log_write('   Finished analysis precompute', masterOnly=.true.)

        if (debugOn) call mpas_log_write('   Starting forward update', masterOnly=.true.)
        call mpas_timer_start("time integration", .false.)
        ierr = 0
        call seaice_timestep(domain, domain % clock, itimestep)
        call mpas_timer_stop("time integration")
        if (debugOn) call mpas_log_write('   Finished forward update', masterOnly=.true.)

        ! update analysis members
        if (debugOn) call mpas_log_write('   Starting AM compute', masterOnly=.true.)
        call seaice_analysis_compute(domain, ierr)
        if (debugOn) call mpas_log_write('   Finished AM compute', masterOnly=.true.)
        if (debugOn) call mpas_log_write('   Starting AM restart', masterOnly=.true.)
        call seaice_analysis_restart(domain, ierr)
        if (debugOn) call mpas_log_write('   Finished AM restart', masterOnly=.true.)
        if (debugOn) call mpas_log_write('   Starting AM write', masterOnly=.true.)
        call seaice_analysis_write(domain, ierr)
        if (debugOn) call mpas_log_write('   Finished AM write', masterOnly=.true.)

        ! Reset the restart alarm to prevent restart files being written without the coupler requesting it.
        if (debugOn) call mpas_log_write('   Resetting restart stream alarms', masterOnly=.true.)
        call mpas_stream_mgr_begin_iteration(domain % streamManager)
        do while ( mpas_stream_mgr_get_next_stream(domain % streamManager, streamID=streamName, &
                   directionProperty=streamDirection, activeProperty=streamActive) )
           if ( streamActive .and. streamDirection == MPAS_STREAM_INPUT_OUTPUT ) then
              call mpas_stream_mgr_reset_alarms(domain % streamManager, streamID=streamName, ierr=ierr)
           end if
        end do
        if (debugOn) call mpas_log_write('   Resetting restart stream alarms complete', masterOnly=.true.)
        if (debugOn) call mpas_log_write('   Writing output streams', masterOnly=.true.)
        call mpas_stream_mgr_write(domain % streamManager, streamID='output', ierr=ierr)
        call mpas_stream_mgr_reset_alarms(domain % streamManager, streamID='output', ierr=ierr)

        call mpas_stream_mgr_write(domain % streamManager, ierr=ierr)
        call mpas_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)
        if (debugOn) call mpas_log_write('   Finished writing output streams', masterOnly=.true.)
        if (debugOn) call mpas_log_write('Completed timestep ' // trim(timeStamp))
     end do

      ! Check if coupler wants us to write a restart file.
      ! We only write restart files at the end of a coupling interval
      if (seq_timemgr_RestartAlarmIsOn(EClock)) then
         ! Write a restart file, because the coupler asked for it.
         if (debugOn) call mpas_log_write('Writing restart streams', masterOnly=.true.)
         call mpas_stream_mgr_begin_iteration(domain % streamManager)
         do while ( mpas_stream_mgr_get_next_stream(domain % streamManager, streamID=streamName, &
                    directionProperty=streamDirection, activeProperty=streamActive) )
            if ( streamActive .and. streamDirection == MPAS_STREAM_INPUT_OUTPUT ) then
               if (debugOn) call mpas_log_write('    Writing stream: ' // trim(streamName), masterOnly=.true.)
               call mpas_stream_mgr_write(domain % streamManager, forceWriteNow=.true., streamID=streamName, iErr=iErr)
               if (debugOn) call mpas_log_write('    Finished Writing stream: ' // trim(streamName), masterOnly=.true.)
            end if
         end do

         if ( iam == 0 ) then
            open(22, file=config_restart_timestamp_name, form='formatted', status='replace')
            write(22, *) trim(timeStamp)
            close(22)
         end if
         if (debugOn) call mpas_log_write('Finished writing restart streams', masterOnly=.true.)
      end if

      ! Export state to coupler
      if (debugOn) call mpas_log_write('Exporting state', masterOnly=.true.)
      call ice_export_mct(i2x_i, ierr)
      if (debugOn) call mpas_log_write('Finished exporting state', masterOnly=.true.)

      ! Check if clocks are in sync
!TODO?      call check_clocks_sync(domain % clock, Eclock, ierr)
      ! Check if clocks are in sync
      currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
      call mpas_get_time(curr_time=currTime, YYYY=iyear, MM=imonth, DD=iday, H=ihour, M=iminute, S=isecond, ierr=ierr)
      call seq_timemgr_EClockGetData(EClock, curr_ymd=curr_ymd, curr_tod=curr_tod)

      ymd = iyear * 10000 + imonth * 100 + iday
      tod = ihour * 3600 + iminute * 60 + isecond
      if (.not. seq_timemgr_EClockDateInSync( EClock, ymd, tod)) then
         call mpas_log_write('MPAS ymd=$i MPAS tod=$i', MPAS_LOG_ERR, intArgs=(/ymd,tod/))
         call mpas_log_write('sync ymd=$i sync tod=$i', MPAS_LOG_ERR, intArgs=(/curr_ymd, curr_tod/))
         call mpas_log_write('Internal mpas clock not in sync with sync clock', MPAS_LOG_CRIT)
      end if

      call mpas_log_write('=== Completed coupling interval in ice_run_mct. ===', flushNow=.true.)

      ! Reset I/O logs
      call shr_file_setLogUnit (shrlogunit)
      call shr_file_setLogLevel(shrloglev)

!-----------------------------------------------------------------------
!EOC

  end subroutine ice_run_mct!}}}

!***********************************************************************
!BOP
!
! !IROUTINE: ice_final_mct
!
! !INTERFACE:
  subroutine ice_final_mct( EClock, cdata_i, x2i_i, i2x_i)!{{{

!
! !DESCRIPTION:
! Finalize MPASSI
!
! !USES:
! !ARGUMENTS:
    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata_i
    type(mct_aVect)             , intent(inout) :: x2i_i
    type(mct_aVect)             , intent(inout) :: i2x_i
!
! !REVISION HISTORY:
! Author: Jon Wolfe
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------


    integer :: shrloglev, shrlogunit, iam
    integer :: ioerr
    character(len=StrKIND) :: ioerrstr

    integer :: ierr

!-----------------------------------------------------------------------
    iam = domain % dminfo % my_proc_id

    ! Set MPAS Log module instance
    mpas_log_info => domain % logInfo

    ! Setup I/O logs
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (iceLogUnit)

    ! Finalize MPASSI
    iErr = domain % core % core_finalize(domain)
    if ( iErr /= 0 ) then
       call mpas_log_write('Core finalize failed for core ' // trim(domain % core % coreName), MPAS_LOG_CRIT)
    end if

    call mpas_timer_write()
    call MPAS_stream_mgr_finalize(domain % streamManager)

    call mpas_log_finalize(iErr)
    if ( iErr /= 0 ) then
       write(iceLogUnit,*) 'ERROR: log finalize failed for core ' // trim(domain % core % coreName)
       call mpas_dmpar_abort(domain % dminfo)
    end if

    call mpas_framework_finalize(domain % dminfo, domain, io_system)

    ! Reset I/O logs
    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

  end subroutine ice_final_mct!}}}

!***********************************************************************
!BOP
!IROUTINE: ice_SetGSMap_mct
! !INTERFACE:

 subroutine ice_SetGSMap_mct( mpicom_ice, ICEID, gsMap_ice )!{{{

! !DESCRIPTION:
!  This routine sets up the MPASSI grid numbering for MCT
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

    implicit none
    integer        , intent(in)    :: mpicom_ice
    integer        , intent(in)    :: ICEID
    type(mct_gsMap), intent(inout) :: gsMap_ice

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
    call mct_gsMap_init( gsMap_ice, gindex, mpicom_ice, ICEID, lsize, gsize )

    deallocate(gindex)

!-----------------------------------------------------------------------
!EOC

  end subroutine ice_SetGSMap_mct!}}}

!***********************************************************************
!BOP
! !IROUTINE: ice_domain_mct
! !INTERFACE:

 subroutine ice_domain_mct( lsize, gsMap_i, dom_i )!{{{

! !DESCRIPTION:
!  This routine sets up the MCT domain for MPASSI
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT/OUTPUT PARAMETERS:

    integer        , intent(in)    :: lsize
    type(mct_gsMap), intent(in)    :: gsMap_i
    type(mct_ggrid), intent(inout) :: dom_i     

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

    type (mpas_pool_type), pointer :: meshPool, oceanCouplingPool

    integer, pointer :: nCellsSolve

    real (kind=RKIND), dimension(:), pointer :: lonCell, latCell, areaCell

    integer, dimension(:), pointer :: landIceMask

    real (kind=RKIND), pointer :: sphere_radius

    r2d = 180.0_RKIND/pii

!-------------------------------------------------------------------
!
!  initialize mct domain type
!
!-------------------------------------------------------------------

    call mct_gGrid_init( GGrid=dom_i, CoordChars=trim(seq_flds_dom_coord), &
       OtherChars=trim(seq_flds_dom_other), lsize=lsize )
    call mct_aVect_zero(dom_i%data)
    allocate(data(lsize))

!-------------------------------------------------------------------
!
! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
!
!-------------------------------------------------------------------

    call mct_gsMap_orderedPoints(gsMap_i, my_task, idata)
    call mct_gGrid_importIAttr(dom_i,'GlobGridNum',idata,lsize)

!-------------------------------------------------------------------
!
! Determine domain 
! Initialize attribute vector with special value
!
!-------------------------------------------------------------------

    data(:) = -9999.0_R8 
    call mct_gGrid_importRAttr(dom_i,"lat"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_i,"lon"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_i,"area" ,data,lsize) 
    call mct_gGrid_importRAttr(dom_i,"aream",data,lsize) 
    data(:) = 1.0_R8     
    call mct_gGrid_importRAttr(dom_i,"mask",data,lsize) 
    call mct_gGrid_importRAttr(dom_i,"frac",data,lsize) 

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
    call mct_gGrid_importRattr(dom_i,"lon",data,lsize) 

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
    call mct_gGrid_importRattr(dom_i,"lat",data,lsize) 

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
    call mct_gGrid_importRattr(dom_i,"area",data,lsize) 

    ! Build mask and frac based on landIceMask
    block_ptr => domain % blocklist
    do while ( associated(block_ptr) )
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'ocean_coupling',oceanCouplingPool)

       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

       call mpas_pool_get_array(oceanCouplingPool, 'landIceMask',landIceMask)

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

    call mct_gGrid_importRattr(dom_i,"mask",data,lsize) 
    call mct_gGrid_importRattr(dom_i,"frac",data,lsize) 

    deallocate(data)
    deallocate(idata)

!-----------------------------------------------------------------------
!EOC

  end subroutine ice_domain_mct!}}}


!***********************************************************************
!BOP
! !IROUTINE: ice_import_mct
! !INTERFACE:

 subroutine ice_import_mct(x2i_i, errorCode)!{{{

! !DESCRIPTION:
!-----------------------------------------------------------------------
!  This routine receives message from cpl7 driver
!
!    The following fields are always received from the coupler:
! 
!    o  t        -- ocn layer temperature
!    o  s        -- ocn salinity
!    o  u        -- ocn u velocity
!    o  v        -- ocn v velocity
!    o  z        -- bottom atm level height
!    o  u        -- bottom atm level zon wind
!    o  v        -- bottom atm level mer wind
!    o  tbot     -- bottom atm level temp
!    o  pbot     -- bottom atm level pressure
!    o  ptem     -- bottom atm level pot temp
!    o  shum     -- bottom atm level spec hum
!    o  dens     -- bottom atm level air den
!    o  dhdx     -- ocn surface slope, zonal
!    o  dhdy     -- ocn surface slope, meridional
!    o  lwdn     -- downward lw heat flux
!    o  rain     -- prec: liquid 
!    o  snow     -- prec: frozen 
!    o  swndr    -- sw: nir direct  downward
!    o  swvdr    -- sw: vis direct  downward
!    o  swndf    -- sw: nir diffuse downward
!    o  swvdf    -- sw: vis diffuse downward
!    o  swnet    -- sw: net
!    o  meltp    -- ocn melt heat  
!    o  frazil   -- ocn frazil production
!    o  bcphidry -- Black Carbon hydrophilic dry deposition flux
!    o  bcphodry -- Black Carbon hydrophobic dry deposition flux
!    o  bcphiwet -- Black Carbon hydrophilic wet deposition flux
!    o  ocphidry -- Organic Carbon hydrophilic dry deposition flux
!    o  ocphodry -- Organic Carbon hydrophobic dry deposition flux
!    o  ocphiwet -- Organic Carbon hydrophilic dry deposition flux
!    o  dstwet1  -- Size 1 dust -- wet deposition flux
!    o  dstwet2  -- Size 2 dust -- wet deposition flux
!    o  dstwet3  -- Size 3 dust -- wet deposition flux
!    o  dstwet4  -- Size 4 dust -- wet deposition flux
!    o  dstdry1  -- Size 1 dust -- dry deposition flux
!    o  dstdry2  -- Size 2 dust -- dry deposition flux
!    o  dstdry3  -- Size 3 dust -- dry deposition flux
!    o  dstdry4  -- Size 4 dust -- dry deposition flux
! 
!    The following fields are sometimes received from the coupler,
!      depending on model options:
! 
!    o  algae1   --
!    o  algae2   --
!    o  algae3   --
!    o  doc1     --
!    o  doc2     --
!    o  doc3     --
!    o  dic1     --
!    o  don1     --
!    o  no3      --
!    o  sio3     --
!    o  nh4      --
!    o  dms      --
!    o  dmsp     --
!    o  docr     --
!    o  fep1     --
!    o  fep2     --
!    o  fed1     --
!    o  fed2     --
!    o  zaer1    --
!    o  zaer2    --
!    o  zaer3    --
!    o  zaer4    --
!    o  zaer5    --
!    o  zaer6    --
! 
!-----------------------------------------------------------------------
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

    type(mct_aVect), intent(inout) :: x2i_i

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

   character (len=StrKIND) :: &
      label,                  &
      message
 
   integer ::  &
      i,n

   real (kind=RKIND) :: &
      frazilProduction, &
      meltingPotential, &
      freezingPotential

   type (block_type), pointer :: block_ptr

   type (mpas_pool_type), pointer :: &
      configs,        &
      meshPool,       &
      aerosols,       &
      atmosCoupling,  &
      oceanCoupling,  &
      biogeochemistry

   integer, pointer :: nCellsSolve

   logical, pointer :: &
      config_use_aerosols,              &
      config_use_modal_aerosols,        &
      config_use_zaerosols,             &
      config_use_column_biogeochemistry

   character(len=strKIND), pointer ::   &
      config_thermodynamics_type,       &
      config_ocean_surface_type

   type (field1DReal), pointer ::       &
      seaSurfaceTemperatureField,       &
      seaSurfaceSalinityField,          &
      seaFreezingTemperatureField,      &
      freezingMeltingPotentialField,    &
      uOceanVelocityField,              &
      vOceanVelocityField,              &
      seaSurfaceTiltUField,             &
      seaSurfaceTiltVField,             &
      airLevelHeightField,              &
      airPotentialTemperatureField,     &
      airTemperatureField,              &
      airSpecificHumidityField,         &
      airDensityField,                  &
      shortwaveVisibleDirectDownField,  &
      shortwaveVisibleDiffuseDownField, &
      shortwaveIRDirectDownField,       &
      shortwaveIRDiffuseDownField,      &
      longwaveDownField,                &
      rainfallRateField,                &
      snowfallRateField,                &
      uAirVelocityField,                &
      vAirVelocityField,                &
      oceanNitrateConcField,            &
      oceanSilicateConcField,           &
      oceanAmmoniumConcField,           &
      oceanDMSConcField,                &
      oceanDMSPConcField,               &
      oceanHumicsConcField

   type (field2DReal), pointer ::       &
      oceanAlgaeConcField,              &
      oceanDOCConcField,                &
      oceanDICConcField,                &
      oceanDONConcField,                &
      oceanParticulateIronConcField,    &
      oceanDissolvedIronConcField,      &
      oceanZAerosolConcField,           &
      atmosAerosolFluxField,            &
      atmosBlackCarbonFluxField,        &
      atmosDustFluxField

   real (kind=RKIND), dimension(:), pointer :: &
      seaSurfaceTemperature,       &
      seaSurfaceSalinity,          &
      seaFreezingTemperature,      &
      freezingMeltingPotential,    &
      uOceanVelocity,              &
      vOceanVelocity,              &
      seaSurfaceTiltU,             &
      seaSurfaceTiltV,             &
      airLevelHeight,              &
      airPotentialTemperature,     &
      airTemperature,              &
      airSpecificHumidity,         &
      airDensity,                  &
      shortwaveVisibleDirectDown,  &
      shortwaveVisibleDiffuseDown, &
      shortwaveIRDirectDown,       &
      shortwaveIRDiffuseDown,      &
      longwaveDown,                &
      rainfallRate,                &
      snowfallRate,                &
      uAirVelocity,                &
      vAirVelocity,                &
      oceanNitrateConc,            &
      oceanSilicateConc,           &
      oceanAmmoniumConc,           &
      oceanDMSConc,                &
      oceanDMSPConc,               &
      oceanHumicsConc,             &
      carbonToNitrogenRatioAlgae,  &
      carbonToNitrogenRatioDON

   real (kind=RKIND), dimension(:,:), pointer :: &
      oceanAlgaeConc,              &
      oceanDOCConc,                &
      oceanDICConc,                &
      oceanDONConc,                &
      oceanParticulateIronConc,    &
      oceanDissolvedIronConc,      &
      oceanZAerosolConc,           &
      atmosAerosolFlux,            &
      atmosBlackCarbonFlux,        &
      atmosDustFlux

!----------------------------------------------------------------------- 
!
!  zero out padded cells 
!
!-----------------------------------------------------------------------

   errorCode = 0

!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------

   n = 0
   block_ptr => domain % blocklist
   do while(associated(block_ptr))

      configs => block_ptr % configs
      call mpas_pool_get_config(configs, "config_thermodynamics_type", config_thermodynamics_type)
      call mpas_pool_get_config(configs, "config_ocean_surface_type", config_ocean_surface_type)
      call mpas_pool_get_config(configs, "config_use_aerosols", config_use_aerosols)
      call mpas_pool_get_config(configs, "config_use_modal_aerosols", config_use_modal_aerosols)
      call mpas_pool_get_config(configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)

      call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block_ptr % structs, 'ocean_coupling', oceanCoupling)
      call mpas_pool_get_subpool(block_ptr % structs, 'atmos_coupling', atmosCoupling)

      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

      call mpas_pool_get_array(oceanCoupling, 'seaSurfaceTemperature', seaSurfaceTemperature)
      call mpas_pool_get_array(oceanCoupling, 'seaSurfaceSalinity', seaSurfaceSalinity)
      call mpas_pool_get_array(oceanCoupling, 'seaFreezingTemperature', seaFreezingTemperature)
      call mpas_pool_get_array(oceanCoupling, 'freezingMeltingPotential', freezingMeltingPotential)
      call mpas_pool_get_array(oceanCoupling, 'uOceanVelocity', uOceanVelocity)
      call mpas_pool_get_array(oceanCoupling, 'vOceanVelocity', vOceanVelocity)
      call mpas_pool_get_array(oceanCoupling, 'seaSurfaceTiltU', seaSurfaceTiltU)
      call mpas_pool_get_array(oceanCoupling, 'seaSurfaceTiltV', seaSurfaceTiltV)

      call mpas_pool_get_array(atmosCoupling, 'airLevelHeight', airLevelHeight)
      call mpas_pool_get_array(atmosCoupling, 'airPotentialTemperature', airPotentialTemperature)
      call mpas_pool_get_array(atmosCoupling, 'airTemperature', airTemperature)
      call mpas_pool_get_array(atmosCoupling, 'airSpecificHumidity', airSpecificHumidity)
      call mpas_pool_get_array(atmosCoupling, 'airDensity', airDensity)
      call mpas_pool_get_array(atmosCoupling, 'shortwaveVisibleDirectDown', shortwaveVisibleDirectDown)
      call mpas_pool_get_array(atmosCoupling, 'shortwaveVisibleDiffuseDown', shortwaveVisibleDiffuseDown)
      call mpas_pool_get_array(atmosCoupling, 'shortwaveIRDirectDown', shortwaveIRDirectDown)
      call mpas_pool_get_array(atmosCoupling, 'shortwaveIRDiffuseDown', shortwaveIRDiffuseDown)
      call mpas_pool_get_array(atmosCoupling, 'longwaveDown', longwaveDown)
      call mpas_pool_get_array(atmosCoupling, 'rainfallRate', rainfallRate)
      call mpas_pool_get_array(atmosCoupling, 'snowfallRate', snowfallRate)
      call mpas_pool_get_array(atmosCoupling, 'uAirVelocity', uAirVelocity)
      call mpas_pool_get_array(atmosCoupling, 'vAirVelocity', vAirVelocity)

      if (config_use_aerosols) then
         call mpas_pool_get_subpool(block_ptr % structs, 'aerosols', aerosols)

         call mpas_pool_get_array(aerosols, "atmosAerosolFlux", atmosAerosolFlux)
      endif

      if (config_use_column_biogeochemistry) then
         call mpas_pool_get_config(configs, "config_use_zaerosols", config_use_zaerosols)
         call mpas_pool_get_subpool(block_ptr % structs, 'biogeochemistry', biogeochemistry)

         call mpas_pool_get_array(biogeochemistry, 'oceanAlgaeConc', oceanAlgaeConc)
         call mpas_pool_get_array(biogeochemistry, 'oceanDOCConc', oceanDOCConc)
         call mpas_pool_get_array(biogeochemistry, 'oceanDICConc', oceanDICConc)
         call mpas_pool_get_array(biogeochemistry, 'oceanDONConc', oceanDONConc)
         call mpas_pool_get_array(biogeochemistry, 'oceanNitrateConc', oceanNitrateConc)
         call mpas_pool_get_array(biogeochemistry, 'oceanSilicateConc', oceanSilicateConc)
         call mpas_pool_get_array(biogeochemistry, 'oceanAmmoniumConc', oceanAmmoniumConc)
         call mpas_pool_get_array(biogeochemistry, 'oceanDMSConc', oceanDMSConc)
         call mpas_pool_get_array(biogeochemistry, 'oceanDMSPConc', oceanDMSPConc)
         call mpas_pool_get_array(biogeochemistry, 'oceanHumicsConc', oceanHumicsConc)
         call mpas_pool_get_array(biogeochemistry, 'oceanParticulateIronConc', oceanParticulateIronConc)
         call mpas_pool_get_array(biogeochemistry, 'oceanDissolvedIronConc', oceanDissolvedIronConc)
         call mpas_pool_get_array(biogeochemistry, 'oceanZAerosolConc', oceanZAerosolConc)
         call mpas_pool_get_array(biogeochemistry, 'carbonToNitrogenRatioAlgae', carbonToNitrogenRatioAlgae)
         call mpas_pool_get_array(biogeochemistry, 'carbonToNitrogenRatioDON', carbonToNitrogenRatioDON)
         if (config_use_zaerosols) then
            call mpas_pool_get_array(biogeochemistry, "atmosBlackCarbonFlux", atmosBlackCarbonFlux)
            call mpas_pool_get_array(biogeochemistry, "atmosDustFlux", atmosDustFlux)
         endif
      endif

      do i = 1, nCellsSolve
        n = n + 1
        seaSurfaceTemperature(i)       = x2i_i % rAttr(index_x2i_So_t, n)
        seaSurfaceSalinity(i)          = x2i_i % rAttr(index_x2i_So_s, n)
        uOceanVelocity(i)              = x2i_i % rAttr(index_x2i_So_u, n)
        vOceanVelocity(i)              = x2i_i % rAttr(index_x2i_So_v, n)
        seaSurfaceTiltU(i)             = x2i_i % rAttr(index_x2i_So_dhdx, n)
        seaSurfaceTiltV(i)             = x2i_i % rAttr(index_x2i_So_dhdy, n)

        if (trim(config_ocean_surface_type) == "free") then ! free surface (MPAS-O)
        
           meltingPotential            = x2i_i % rAttr(index_x2i_Fioo_meltp, n)

           ! frazil calculation
           frazilProduction            = x2i_i % rAttr(index_x2i_Fioo_frazil, n)
           call freezing_potential(freezingPotential, frazilProduction, seaSurfaceSalinity(i), &
                                   config_thermodynamics_type)
           freezingMeltingPotential(i) = freezingPotential + meltingPotential

        else ! non-free surface (SOM)

           freezingMeltingPotential(i) = x2i_i % rAttr(index_x2i_Fioo_q, n)

        endif

        airLevelHeight(i)              = x2i_i % rAttr(index_x2i_Sa_z, n)
        airPotentialTemperature(i)     = x2i_i % rAttr(index_x2i_Sa_ptem, n)
        airTemperature(i)              = x2i_i % rAttr(index_x2i_Sa_tbot, n)
        airSpecificHumidity(i)         = x2i_i % rAttr(index_x2i_Sa_shum, n)
        airDensity(i)                  = x2i_i % rAttr(index_x2i_Sa_dens, n)
        shortwaveVisibleDirectDown(i)  = x2i_i % rAttr(index_x2i_Faxa_swvdr, n)
        shortwaveVisibleDiffuseDown(i) = x2i_i % rAttr(index_x2i_Faxa_swvdf, n)
        shortwaveIRDirectDown(i)       = x2i_i % rAttr(index_x2i_Faxa_swndr, n)
        shortwaveIRDiffuseDown(i)      = x2i_i % rAttr(index_x2i_Faxa_swndf, n)
        longwaveDown(i)                = x2i_i % rAttr(index_x2i_Faxa_lwdn, n)
        rainfallRate(i)                = x2i_i % rAttr(index_x2i_Faxa_rain, n)
        snowfallRate(i)                = x2i_i % rAttr(index_x2i_Faxa_snow, n)
        uAirVelocity(i)                = x2i_i % rAttr(index_x2i_Sa_u, n)
        vAirVelocity(i)                = x2i_i % rAttr(index_x2i_Sa_v, n)

        ! set aerosols, if configured
        if (config_use_aerosols) then
           if (config_use_modal_aerosols) then
              atmosAerosolFlux(1,i) = x2i_i % rAttr(index_x2i_Faxa_bcphodry, n) &
                                    + x2i_i % rAttr(index_x2i_Faxa_bcphidry, n)
              atmosAerosolFlux(2,i) = x2i_i % rAttr(index_x2i_Faxa_bcphiwet, n)
              ! combine all the dust into one category
              atmosAerosolFlux(3,i) = x2i_i % rAttr(index_x2i_Faxa_dstwet1, n) &
                                    + x2i_i % rAttr(index_x2i_Faxa_dstwet2, n) &
                                    + x2i_i % rAttr(index_x2i_Faxa_dstwet3, n) &
                                    + x2i_i % rAttr(index_x2i_Faxa_dstwet4, n) &
                                    + x2i_i % rAttr(index_x2i_Faxa_dstdry1, n) &
                                    + x2i_i % rAttr(index_x2i_Faxa_dstdry2, n) &
                                    + x2i_i % rAttr(index_x2i_Faxa_dstdry3, n) &
                                    + x2i_i % rAttr(index_x2i_Faxa_dstdry4, n)
           else
              atmosAerosolFlux(1,i) = x2i_i % rAttr(index_x2i_Faxa_bcphodry, n)
              atmosAerosolFlux(2,i) = x2i_i % rAttr(index_x2i_Faxa_bcphidry, n) &
                                    + x2i_i % rAttr(index_x2i_Faxa_bcphiwet, n)
              ! combine all the dust into one category
              atmosAerosolFlux(3,i) = x2i_i % rAttr(index_x2i_Faxa_dstwet1, n) &
                                    + x2i_i % rAttr(index_x2i_Faxa_dstwet2, n) &
                                    + x2i_i % rAttr(index_x2i_Faxa_dstwet3, n) &
                                    + x2i_i % rAttr(index_x2i_Faxa_dstwet4, n) &
                                    + x2i_i % rAttr(index_x2i_Faxa_dstdry1, n) &
                                    + x2i_i % rAttr(index_x2i_Faxa_dstdry2, n) &
                                    + x2i_i % rAttr(index_x2i_Faxa_dstdry3, n) &
                                    + x2i_i % rAttr(index_x2i_Faxa_dstdry4, n)
           endif
        endif

        ! import biogeochemistry fields, if configured
        if (config_use_column_biogeochemistry) then
           oceanAlgaeConc(1,i)           = x2i_i % rAttr(index_x2i_So_algae1, n)
           oceanAlgaeConc(2,i)           = x2i_i % rAttr(index_x2i_So_algae2, n)
           oceanAlgaeConc(3,i)           = x2i_i % rAttr(index_x2i_So_algae3, n)
           oceanDOCConc(1,i)             = x2i_i % rAttr(index_x2i_So_doc1, n)
           oceanDOCConc(2,i)             = x2i_i % rAttr(index_x2i_So_doc2, n)
           oceanDOCConc(3,i)             = 0.0_RKIND
           oceanDICConc(1,i)             = x2i_i % rAttr(index_x2i_So_dic1, n) !JW not used, set to 0?
           oceanDONConc(1,i)             = x2i_i % rAttr(index_x2i_So_don1, n)
           oceanNitrateConc(i)           = x2i_i % rAttr(index_x2i_So_no3, n)
           oceanSilicateConc(i)          = x2i_i % rAttr(index_x2i_So_sio3, n)
           oceanAmmoniumConc(i)          = x2i_i % rAttr(index_x2i_So_nh4, n)
           oceanDMSConc(i)               = x2i_i % rAttr(index_x2i_So_dms, n)
           oceanDMSPConc(i)              = x2i_i % rAttr(index_x2i_So_dmsp, n)
           oceanHumicsConc(i)            = x2i_i % rAttr(index_x2i_So_docr, n)
           oceanParticulateIronConc(1,i) = x2i_i % rAttr(index_x2i_So_fep1, n)
           oceanParticulateIronConc(2,i) = x2i_i % rAttr(index_x2i_So_fep2, n)
           oceanDissolvedIronConc(1,i)   = x2i_i % rAttr(index_x2i_So_fed1, n)
           oceanDissolvedIronConc(2,i)   = x2i_i % rAttr(index_x2i_So_fed2, n)
           oceanZAerosolConc(1,i)        = x2i_i % rAttr(index_x2i_So_zaer1, n) !JW not used, set to 0?
           oceanZAerosolConc(2,i)        = x2i_i % rAttr(index_x2i_So_zaer2, n) !JW not used, set to 0?
           oceanZAerosolConc(3,i)        = x2i_i % rAttr(index_x2i_So_zaer3, n) !JW not used, set to 0?
           oceanZAerosolConc(4,i)        = x2i_i % rAttr(index_x2i_So_zaer4, n) !JW not used, set to 0?
           oceanZAerosolConc(5,i)        = x2i_i % rAttr(index_x2i_So_zaer5, n) !JW not used, set to 0?
           oceanZAerosolConc(6,i)        = x2i_i % rAttr(index_x2i_So_zaer6, n) !JW not used, set to 0?
           ! set aerosols, if configured
           if (config_use_zaerosols) then
              if (config_use_modal_aerosols) then
                 atmosBlackCarbonFlux(1,i) = x2i_i % rAttr(index_x2i_Faxa_bcphodry, n) &
                                           + x2i_i % rAttr(index_x2i_Faxa_bcphidry, n)
                 atmosBlackCarbonFlux(2,i) = x2i_i % rAttr(index_x2i_Faxa_bcphiwet, n)
                 ! combine wet and dry dust
                 atmosDustFlux(1,i) = x2i_i % rAttr(index_x2i_Faxa_dstwet1, n) &
                                    + x2i_i % rAttr(index_x2i_Faxa_dstdry1, n)
                 atmosDustFlux(2,i) = x2i_i % rAttr(index_x2i_Faxa_dstwet2, n) &
                                    + x2i_i % rAttr(index_x2i_Faxa_dstdry2, n)
                 atmosDustFlux(3,i) = x2i_i % rAttr(index_x2i_Faxa_dstwet3, n) &
                                    + x2i_i % rAttr(index_x2i_Faxa_dstdry3, n)
                 atmosDustFlux(4,i) = x2i_i % rAttr(index_x2i_Faxa_dstwet4, n) &
                                    + x2i_i % rAttr(index_x2i_Faxa_dstdry4, n)
              else
                 atmosBlackCarbonFlux(1,i) = x2i_i % rAttr(index_x2i_Faxa_bcphodry, n)
                 atmosBlackCarbonFlux(2,i) = x2i_i % rAttr(index_x2i_Faxa_bcphidry, n) &
                                           + x2i_i % rAttr(index_x2i_Faxa_bcphiwet, n)
                 ! combine wet and dry dust
                 atmosDustFlux(1,i) = x2i_i % rAttr(index_x2i_Faxa_dstwet1, n) &
                                    + x2i_i % rAttr(index_x2i_Faxa_dstdry1, n)
                 atmosDustFlux(2,i) = x2i_i % rAttr(index_x2i_Faxa_dstwet2, n) &
                                    + x2i_i % rAttr(index_x2i_Faxa_dstdry2, n)
                 atmosDustFlux(3,i) = x2i_i % rAttr(index_x2i_Faxa_dstwet3, n) &
                                    + x2i_i % rAttr(index_x2i_Faxa_dstdry3, n)
                 atmosDustFlux(4,i) = x2i_i % rAttr(index_x2i_Faxa_dstwet4, n) &
                                    + x2i_i % rAttr(index_x2i_Faxa_dstdry4, n)
              endif
           endif
        endif
      end do

!----------------------------------------------------------------------- 
!
!  unit conversions and any manipulation of coupled fields
!
!-----------------------------------------------------------------------
      do i = 1, nCellsSolve
        seaSurfaceTemperature(i)  = seaSurfaceTemperature(i) - seaiceFreshWaterFreezingPoint
        seaFreezingTemperature(i) = -1.8_RKIND                 ! hardwired for NCOM
             !         Tf (i,j,iblk) = -depressT*sss(i,j,iblk)       ! freezing temp (C)
             !         Tf (i,j,iblk) = -depressT*max(sss(i,j,iblk),ice_ref_salinity)

        if (config_use_column_biogeochemistry) then
           ! convert from mmol C/m^3 to mmol N/m^3
           oceanAlgaeConc(1,i)           = oceanAlgaeConc(1,i) / carbonToNitrogenRatioAlgae(1)
           oceanAlgaeConc(2,i)           = oceanAlgaeConc(2,i) / carbonToNitrogenRatioAlgae(2)
           oceanAlgaeConc(3,i)           = oceanAlgaeConc(3,i) / carbonToNitrogenRatioAlgae(3)
           ! convert from mmol Fe/m^3 to umol Fe/m^3
           oceanParticulateIronConc(1,i) = oceanParticulateIronConc(1,i) * 1000._RKIND
           oceanParticulateIronConc(2,i) = oceanParticulateIronConc(2,i) * 1000._RKIND
           oceanDissolvedIronConc(1,i)   = oceanDissolvedIronConc(1,i)   * 1000._RKIND
           oceanDissolvedIronConc(2,i)   = oceanDissolvedIronConc(2,i)   * 1000._RKIND
        endif
      end do

      block_ptr => block_ptr % next
   end do

   call mpas_pool_get_subpool(domain % blocklist % structs, 'ocean_coupling', oceanCoupling)
   call mpas_pool_get_subpool(domain % blocklist % structs, 'atmos_coupling', atmosCoupling)

   call mpas_pool_get_field(oceanCoupling, 'seaSurfaceTemperature', seaSurfaceTemperatureField)
   call mpas_pool_get_field(oceanCoupling, 'seaSurfaceSalinity', seaSurfaceSalinityField)
   call mpas_pool_get_field(oceanCoupling, 'seaFreezingTemperature', seaFreezingTemperatureField)
   call mpas_pool_get_field(oceanCoupling, 'freezingMeltingPotential', freezingMeltingPotentialField)
   call mpas_pool_get_field(oceanCoupling, 'uOceanVelocity', uOceanVelocityField)
   call mpas_pool_get_field(oceanCoupling, 'vOceanVelocity', vOceanVelocityField)
   call mpas_pool_get_field(oceanCoupling, 'seaSurfaceTiltU', seaSurfaceTiltUField)
   call mpas_pool_get_field(oceanCoupling, 'seaSurfaceTiltV', seaSurfaceTiltVField)

   call mpas_pool_get_field(atmosCoupling, 'airLevelHeight', airLevelHeightField)
   call mpas_pool_get_field(atmosCoupling, 'airPotentialTemperature', airPotentialTemperatureField)
   call mpas_pool_get_field(atmosCoupling, 'airTemperature', airTemperatureField)
   call mpas_pool_get_field(atmosCoupling, 'airSpecificHumidity', airSpecificHumidityField)
   call mpas_pool_get_field(atmosCoupling, 'airDensity', airDensityField)
   call mpas_pool_get_field(atmosCoupling, 'shortwaveVisibleDirectDown', shortwaveVisibleDirectDownField)
   call mpas_pool_get_field(atmosCoupling, 'shortwaveVisibleDiffuseDown', shortwaveVisibleDiffuseDownField)
   call mpas_pool_get_field(atmosCoupling, 'shortwaveIRDirectDown', shortwaveIRDirectDownField)
   call mpas_pool_get_field(atmosCoupling, 'shortwaveIRDiffuseDown', shortwaveIRDiffuseDownField)
   call mpas_pool_get_field(atmosCoupling, 'longwaveDown', longwaveDownField)
   call mpas_pool_get_field(atmosCoupling, 'rainfallRate', rainfallRateField)
   call mpas_pool_get_field(atmosCoupling, 'snowfallRate', snowfallRateField)
   call mpas_pool_get_field(atmosCoupling, 'uAirVelocity', uAirVelocityField)
   call mpas_pool_get_field(atmosCoupling, 'vAirVelocity', vAirVelocityField)

   if (config_use_aerosols) then
      call mpas_pool_get_subpool(domain % blocklist % structs, 'aerosols', aerosols)

      call mpas_pool_get_field(aerosols, "atmosAerosolFlux", atmosAerosolFluxField)
   endif

   if (config_use_column_biogeochemistry) then
      call mpas_pool_get_subpool(domain % blocklist % structs, 'biogeochemistry', biogeochemistry)

      call mpas_pool_get_field(biogeochemistry, 'oceanAlgaeConc', oceanAlgaeConcField)
      call mpas_pool_get_field(biogeochemistry, 'oceanDOCConc', oceanDOCConcField)
      call mpas_pool_get_field(biogeochemistry, 'oceanDICConc', oceanDICConcField)
      call mpas_pool_get_field(biogeochemistry, 'oceanDONConc', oceanDONConcField)
      call mpas_pool_get_field(biogeochemistry, 'oceanNitrateConc', oceanNitrateConcField)
      call mpas_pool_get_field(biogeochemistry, 'oceanSilicateConc', oceanSilicateConcField)
      call mpas_pool_get_field(biogeochemistry, 'oceanAmmoniumConc', oceanAmmoniumConcField)
      call mpas_pool_get_field(biogeochemistry, 'oceanDMSConc', oceanDMSConcField)
      call mpas_pool_get_field(biogeochemistry, 'oceanDMSPConc', oceanDMSPConcField)
      call mpas_pool_get_field(biogeochemistry, 'oceanHumicsConc', oceanHumicsConcField)
      call mpas_pool_get_field(biogeochemistry, 'oceanParticulateIronConc', oceanParticulateIronConcField)
      call mpas_pool_get_field(biogeochemistry, 'oceanDissolvedIronConc', oceanDissolvedIronConcField)
      call mpas_pool_get_field(biogeochemistry, 'oceanZAerosolConc', oceanZAerosolConcField)
      if (config_use_zaerosols) then
         call mpas_pool_get_field(biogeochemistry, "atmosBlackCarbonFlux", atmosBlackCarbonFluxField)
         call mpas_pool_get_field(biogeochemistry, "atmosDustFlux", atmosDustFluxField)
      endif
   endif

   call mpas_dmpar_exch_halo_field(seaSurfaceTemperatureField)
   call mpas_dmpar_exch_halo_field(seaSurfaceSalinityField)
   call mpas_dmpar_exch_halo_field(seaFreezingTemperatureField)
   call mpas_dmpar_exch_halo_field(freezingMeltingPotentialField)
   call mpas_dmpar_exch_halo_field(uOceanVelocityField)
   call mpas_dmpar_exch_halo_field(vOceanVelocityField)
   call mpas_dmpar_exch_halo_field(seaSurfaceTiltUField)
   call mpas_dmpar_exch_halo_field(seaSurfaceTiltVField)

   call mpas_dmpar_exch_halo_field(airLevelHeightField)
   call mpas_dmpar_exch_halo_field(airPotentialTemperatureField)
   call mpas_dmpar_exch_halo_field(airTemperatureField)
   call mpas_dmpar_exch_halo_field(airSpecificHumidityField)
   call mpas_dmpar_exch_halo_field(airDensityField)
   call mpas_dmpar_exch_halo_field(shortwaveVisibleDirectDownField)
   call mpas_dmpar_exch_halo_field(shortwaveVisibleDiffuseDownField)
   call mpas_dmpar_exch_halo_field(shortwaveIRDirectDownField)
   call mpas_dmpar_exch_halo_field(shortwaveIRDiffuseDownField)
   call mpas_dmpar_exch_halo_field(longwaveDownField)
   call mpas_dmpar_exch_halo_field(rainfallRateField)
   call mpas_dmpar_exch_halo_field(snowfallRateField)
   call mpas_dmpar_exch_halo_field(uAirVelocityField)
   call mpas_dmpar_exch_halo_field(vAirVelocityField)

   if (config_use_aerosols) then
      call mpas_dmpar_exch_halo_field(atmosAerosolFluxField)
   endif

   if (config_use_column_biogeochemistry) then
      call mpas_dmpar_exch_halo_field(oceanAlgaeConcField)
      call mpas_dmpar_exch_halo_field(oceanDOCConcField)
      call mpas_dmpar_exch_halo_field(oceanDICConcField)
      call mpas_dmpar_exch_halo_field(oceanDONConcField)
      call mpas_dmpar_exch_halo_field(oceanNitrateConcField)
      call mpas_dmpar_exch_halo_field(oceanSilicateConcField)
      call mpas_dmpar_exch_halo_field(oceanAmmoniumConcField)
      call mpas_dmpar_exch_halo_field(oceanDMSConcField)
      call mpas_dmpar_exch_halo_field(oceanDMSPConcField)
      call mpas_dmpar_exch_halo_field(oceanHumicsConcField)
      call mpas_dmpar_exch_halo_field(oceanParticulateIronConcField)
      call mpas_dmpar_exch_halo_field(oceanDissolvedIronConcField)
      call mpas_dmpar_exch_halo_field(oceanZAerosolConcField)
      if (config_use_zaerosols) then
         call mpas_dmpar_exch_halo_field(atmosBlackCarbonFluxField)
         call mpas_dmpar_exch_halo_field(atmosDustFluxField)
      endif
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine ice_import_mct!}}}
!***********************************************************************
!BOP
! !IROUTINE: ice_export_mct
! !INTERFACE:

 subroutine ice_export_mct(i2x_i, errorCode)   !{{{

! !DESCRIPTION:
!  This routine calls the routines necessary to send MPASSI fields to
!  the CCSM cpl7 driver
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT/OUTPUT PARAMETERS:

   type(mct_aVect)   , intent(inout) :: i2x_i

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

   integer ::  &
      i, n

   real(kind=RKIND) :: &
      ailohi,          &
      Tsrf,            &
      tauxa,           &
      tauya,           &
      tauxo,           &
      tauyo,           &
      basalPressure

   type (block_type), pointer :: &
      block_ptr

   type (mpas_pool_type), pointer :: &
      configs,          &
      meshPool,         &
      tracersAggregate, &
      velocitySolver,   &
      shortwave,        &
      atmosCoupling,    &
      atmosFluxes,      &
      oceanFluxes,      &
      icebergFluxes,    &
      biogeochemistry

   integer, pointer :: &
      nCellsSolve

   logical, pointer :: &
      config_rotate_cartesian_grid,      &
      config_use_topo_meltponds,         &
      config_use_column_biogeochemistry, &
      config_use_column_shortwave,       &
      config_use_data_icebergs

   real(kind=RKIND), pointer :: &
      sphere_radius

   real (kind=RKIND), dimension(:), pointer :: &
      latCell,                     &
      lonCell,                     &
      xCell,                       &
      yCell,                       &
      zCell,                       &
      iceAreaCell,                 &
      iceVolumeCell,               &
      snowVolumeCell,              &
      pondDepthCell,               &
      pondLidThicknessCell,        &
      pondAreaCell,                &
      surfaceTemperatureCell,      &
      airStressCellU,              &
      airStressCellV,              &
      oceanStressCellU,            &
      oceanStressCellV,            &
      albedoVisibleDirectCell,     &
      albedoIRDirectCell,          &
      albedoVisibleDiffuseCell,    &
      albedoIRDiffuseCell,         &
      atmosReferenceSpeed10m,      &
      atmosReferenceTemperature2m, &
      atmosReferenceHumidity2m,    &
      latentHeatFlux,              &
      sensibleHeatFlux,            &
      longwaveUp,                  &
      evaporativeWaterFlux,        &
      absorbedShortwaveFlux,       &
      oceanHeatFlux,               &
      oceanShortwaveFlux,          &
      oceanFreshWaterFlux,         &
      oceanSaltFlux,               &
      bergFreshwaterFlux,          &
      bergLatentHeatFlux,          &
      oceanNitrateFlux,            &
      oceanSilicateFlux,           &
      oceanAmmoniumFlux,           &
      oceanDMSFlux,                &
      oceanDMSPpFlux,              &
      oceanDMSPdFlux,              &
      oceanHumicsFlux,             &
      carbonToNitrogenRatioAlgae,  &
      carbonToNitrogenRatioDON

   real (kind=RKIND), dimension(:,:), pointer :: &
      oceanAlgaeFlux,              &
      oceanDOCFlux,                &
      oceanDICFlux,                &
      oceanDONFlux,                &
      oceanParticulateIronFlux,    &
      oceanDissolvedIronFlux

!----------------------------------------------------------------------- 

   errorCode = 0
   n = 0
   i2x_i % rAttr(: ,:) = 0.0_RKIND
   block_ptr => domain % blocklist
   do while(associated(block_ptr))

      configs => block_ptr % configs
      call MPAS_pool_get_config(configs, "config_rotate_cartesian_grid", config_rotate_cartesian_grid)
      call MPAS_pool_get_config(configs, "config_use_topo_meltponds", config_use_topo_meltponds)
      call MPAS_pool_get_config(configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)
      call MPAS_pool_get_config(configs, "config_use_column_shortwave", config_use_column_shortwave)
      call MPAS_pool_get_config(configs, "config_use_data_icebergs", config_use_data_icebergs)

      call MPAS_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
      call MPAS_pool_get_subpool(block_ptr % structs, "tracers_aggregate", tracersAggregate)
      call MPAS_pool_get_subpool(block_ptr % structs, "velocity_solver", velocitySolver)
      call MPAS_pool_get_subpool(block_ptr % structs, "shortwave", shortwave)
      call MPAS_pool_get_subpool(block_ptr % structs, 'atmos_coupling', atmosCoupling)
      call MPAS_pool_get_subpool(block_ptr % structs, "atmos_fluxes", atmosFluxes)
      call MPAS_pool_get_subpool(block_ptr % structs, "ocean_fluxes", oceanFluxes)

      call MPAS_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call MPAS_pool_get_config(meshPool, "sphere_radius", sphere_radius)
      call MPAS_pool_get_array(meshPool, "latCell", latCell)
      call MPAS_pool_get_array(meshPool, "lonCell", lonCell)
      call MPAS_pool_get_array(meshPool, "xCell", xCell)
      call MPAS_pool_get_array(meshPool, "yCell", yCell)
      call MPAS_pool_get_array(meshPool, "zCell", zCell)

      call MPAS_pool_get_array(tracersAggregate, 'iceAreaCell', iceAreaCell)
      call MPAS_pool_get_array(tracersAggregate, 'iceVolumeCell', iceVolumeCell)
      call MPAS_pool_get_array(tracersAggregate, 'snowVolumeCell', snowVolumeCell)
      call MPAS_pool_get_array(tracersAggregate, 'pondDepthCell', pondDepthCell)
      call MPAS_pool_get_array(tracersAggregate, 'pondLidThicknessCell', pondLidThicknessCell)
      call MPAS_pool_get_array(tracersAggregate, 'pondAreaCell', pondAreaCell)
      call MPAS_pool_get_array(tracersAggregate, 'surfaceTemperatureCell', surfaceTemperatureCell)

      call MPAS_pool_get_array(velocitySolver, 'airStressCellU', airStressCellU)
      call MPAS_pool_get_array(velocitySolver, 'airStressCellV', airStressCellV)
      call MPAS_pool_get_array(velocitySolver, 'oceanStressCellU', oceanStressCellU)
      call MPAS_pool_get_array(velocitySolver, 'oceanStressCellV', oceanStressCellV)

      call MPAS_pool_get_array(shortwave, 'albedoVisibleDirectCell', albedoVisibleDirectCell)
      call MPAS_pool_get_array(shortwave, 'albedoIRDirectCell', albedoIRDirectCell)
      call MPAS_pool_get_array(shortwave, 'albedoVisibleDiffuseCell', albedoVisibleDiffuseCell)
      call MPAS_pool_get_array(shortwave, 'albedoIRDiffuseCell', albedoIRDiffuseCell)
      call MPAS_pool_get_array(shortwave, 'absorbedShortwaveFlux', absorbedShortwaveFlux)

      call MPAS_pool_get_array(atmosCoupling, 'atmosReferenceSpeed10m', atmosReferenceSpeed10m)
      call MPAS_pool_get_array(atmosCoupling, 'atmosReferenceTemperature2m', atmosReferenceTemperature2m)
      call MPAS_pool_get_array(atmosCoupling, 'atmosReferenceHumidity2m', atmosReferenceHumidity2m)

      call MPAS_pool_get_array(atmosFluxes, 'latentHeatFlux', latentHeatFlux)
      call MPAS_pool_get_array(atmosFluxes, 'sensibleHeatFlux', sensibleHeatFlux)
      call MPAS_pool_get_array(atmosFluxes, 'longwaveUp', longwaveUp)
      call MPAS_pool_get_array(atmosFluxes, 'evaporativeWaterFlux', evaporativeWaterFlux)

      call MPAS_pool_get_array(oceanFluxes, 'oceanHeatFlux', oceanHeatFlux)
      call MPAS_pool_get_array(oceanFluxes, 'oceanShortwaveFlux', oceanShortwaveFlux)
      call MPAS_pool_get_array(oceanFluxes, 'oceanFreshWaterFlux', oceanFreshWaterFlux)
      call MPAS_pool_get_array(oceanFluxes, 'oceanSaltFlux', oceanSaltFlux)

      if (config_use_data_icebergs) then
         call MPAS_pool_get_subpool(block_ptr % structs, "berg_fluxes", icebergFluxes)

         call MPAS_pool_get_array(icebergFluxes, "bergFreshwaterFlux", bergFreshwaterFlux)
         call MPAS_pool_get_array(icebergFluxes, "bergLatentHeatFlux", bergLatentHeatFlux)
      endif

      if (config_use_column_biogeochemistry) then
         call mpas_pool_get_subpool(block_ptr % structs, 'biogeochemistry', biogeochemistry)

         call mpas_pool_get_array(biogeochemistry, 'oceanAlgaeFlux', oceanAlgaeFlux)
         call mpas_pool_get_array(biogeochemistry, 'oceanDOCFlux', oceanDOCFlux)
         call mpas_pool_get_array(biogeochemistry, 'oceanDICFlux', oceanDICFlux)
         call mpas_pool_get_array(biogeochemistry, 'oceanDONFlux', oceanDONFlux)
         call mpas_pool_get_array(biogeochemistry, 'oceanNitrateFlux', oceanNitrateFlux)
         call mpas_pool_get_array(biogeochemistry, 'oceanSilicateFlux', oceanSilicateFlux)
         call mpas_pool_get_array(biogeochemistry, 'oceanAmmoniumFlux', oceanAmmoniumFlux)
         call mpas_pool_get_array(biogeochemistry, 'oceanDMSFlux', oceanDMSFlux)
         call mpas_pool_get_array(biogeochemistry, 'oceanDMSPpFlux', oceanDMSPpFlux)
         call mpas_pool_get_array(biogeochemistry, 'oceanDMSPdFlux', oceanDMSPdFlux)
         call mpas_pool_get_array(biogeochemistry, 'oceanHumicsFlux', oceanHumicsFlux)
         call mpas_pool_get_array(biogeochemistry, 'oceanParticulateIronFlux', oceanParticulateIronFlux)
         call mpas_pool_get_array(biogeochemistry, 'oceanDissolvedIronFlux', oceanDissolvedIronFlux)
         call mpas_pool_get_array(biogeochemistry, 'carbonToNitrogenRatioAlgae', carbonToNitrogenRatioAlgae)
         call mpas_pool_get_array(biogeochemistry, 'carbonToNitrogenRatioDON', carbonToNitrogenRatioDON)
      endif

      do i = 1, nCellsSolve
         n = n + 1

         ! ice fraction
         ailohi = min(iceAreaCell(i), 1.0_RKIND)

         !TODO: CICE has a check for ailohi < 0

         ! surface temperature
         Tsrf = seaiceFreshWaterFreezingPoint + surfaceTemperatureCell(i)

         ! basal pressure
         if ( ailohi > 0.0_RKIND ) then
            call basal_pressure(&
                 basalPressure,           &
                 iceVolumeCell(i),        &
                 snowVolumeCell(i),       &
                 pondDepthCell(i),        &
                 pondLidThicknessCell(i), &
                 pondAreaCell(i),         &
                 config_use_topo_meltponds)
         endif

         ! wind stress  (on T-grid:  convert to lat-lon)
         call seaice_latlon_vector_rotation_backward(&
              tauxa,             &
              tauya,             &
              airStressCellU(i), &
              airStressCellV(i), &
              latCell(i),        &
              lonCell(i),        &
              xCell(i),          &
              yCell(i),          &
              zCell(i),          &
              sphere_radius,     &
              config_rotate_cartesian_grid)

         ! ice/ocean stress (on POP T-grid:  convert to lat-lon)
         call seaice_latlon_vector_rotation_backward(&
              tauxo,                &
              tauyo,                &
              -oceanStressCellU(i), &
              -oceanStressCellV(i), &
              latCell(i),           &
              lonCell(i),           &
              xCell(i),             &
              yCell(i),             &
              zCell(i),             &
              sphere_radius,        &
              config_rotate_cartesian_grid)

         !-------states-------------------- 
         i2x_i % rAttr(index_i2x_Si_ifrac ,n)    = ailohi

         if (config_use_data_icebergs) then
            i2x_i % rAttr(index_i2x_Fioi_bergw,n) = bergFreshwaterFlux(i)
            i2x_i % rAttr(index_i2x_Fioi_bergh,n) = bergLatentHeatFlux(i)
         endif

         if ( ailohi > 0.0_RKIND ) then

            !-------states-------------------- 
            i2x_i % rAttr(index_i2x_Si_t     ,n)  = Tsrf
            i2x_i % rAttr(index_i2x_Si_bpress,n)  = basalPressure
            i2x_i % rAttr(index_i2x_Si_u10   ,n)  = atmosReferenceSpeed10m(i)
            i2x_i % rAttr(index_i2x_Si_tref  ,n)  = atmosReferenceTemperature2m(i)
            i2x_i % rAttr(index_i2x_Si_qref  ,n)  = atmosReferenceHumidity2m(i)
            i2x_i % rAttr(index_i2x_Si_snowh ,n)  = snowVolumeCell(i) / ailohi

            !--- a/i fluxes computed by ice
            i2x_i % rAttr(index_i2x_Faii_taux ,n) = tauxa
            i2x_i % rAttr(index_i2x_Faii_tauy ,n) = tauya
            i2x_i % rAttr(index_i2x_Faii_lat  ,n) = latentHeatFlux(i)
            i2x_i % rAttr(index_i2x_Faii_sen  ,n) = sensibleHeatFlux(i)
            i2x_i % rAttr(index_i2x_Faii_lwup ,n) = longwaveUp(i)
            i2x_i % rAttr(index_i2x_Faii_evap ,n) = evaporativeWaterFlux(i)
            i2x_i % rAttr(index_i2x_Faii_swnet,n) = absorbedShortwaveFlux(i)
            i2x_i % rAttr(index_i2x_Faii_evap ,n) = evaporativeWaterFlux(i)

            if (config_use_column_shortwave) then
               i2x_i % rAttr(index_i2x_Si_avsdr ,n)  = albedoVisibleDirectCell(i)
               i2x_i % rAttr(index_i2x_Si_anidr ,n)  = albedoIRDirectCell(i)
               i2x_i % rAttr(index_i2x_Si_avsdf ,n)  = albedoVisibleDiffuseCell(i)
               i2x_i % rAttr(index_i2x_Si_anidf ,n)  = albedoIRDiffuseCell(i)

               i2x_i % rAttr(index_i2x_Faii_swnet,n) = absorbedShortwaveFlux(i)
            endif

            !--- i/o fluxes computed by ice
            i2x_i % rAttr(index_i2x_Fioi_melth,n) = oceanHeatFlux(i)
            i2x_i % rAttr(index_i2x_Fioi_swpen,n) = oceanShortwaveFlux(i)
            i2x_i % rAttr(index_i2x_Fioi_meltw,n) = oceanFreshWaterFlux(i)
            i2x_i % rAttr(index_i2x_Fioi_salt ,n) = oceanSaltFlux(i)
            i2x_i % rAttr(index_i2x_Fioi_taux ,n) = tauxo
            i2x_i % rAttr(index_i2x_Fioi_tauy ,n) = tauyo

            ! export biogeochemistry fields, if configured
            if (config_use_column_biogeochemistry) then
               ! convert from mmol N/m^3 to mmol C/m^3
               i2x_i % rAttr(index_i2x_Fioi_algae1,n) = oceanAlgaeFlux(1,i) * carbonToNitrogenRatioAlgae(1)
               i2x_i % rAttr(index_i2x_Fioi_algae2,n) = oceanAlgaeFlux(2,i) * carbonToNitrogenRatioAlgae(2)
               i2x_i % rAttr(index_i2x_Fioi_algae3,n) = oceanAlgaeFlux(3,i) * carbonToNitrogenRatioAlgae(3)
               i2x_i % rAttr(index_i2x_Fioi_doc1  ,n) = oceanDOCFlux(1,i)
               i2x_i % rAttr(index_i2x_Fioi_doc2  ,n) = oceanDOCFlux(2,i)
               i2x_i % rAttr(index_i2x_Fioi_doc3  ,n) = oceanDOCFlux(3,i) !JW set to 0?
               i2x_i % rAttr(index_i2x_Fioi_dic1  ,n) = oceanDICFlux(1,i)
               i2x_i % rAttr(index_i2x_Fioi_don1  ,n) = oceanDONFlux(1,i)
               i2x_i % rAttr(index_i2x_Fioi_no3   ,n) = oceanNitrateFlux(i)
               i2x_i % rAttr(index_i2x_Fioi_sio3  ,n) = oceanSilicateFlux(i)
               i2x_i % rAttr(index_i2x_Fioi_nh4   ,n) = oceanAmmoniumFlux(i)
               i2x_i % rAttr(index_i2x_Fioi_dms   ,n) = oceanDMSFlux(i)
               i2x_i % rAttr(index_i2x_Fioi_dmspp ,n) = oceanDMSPpFlux(i)
               i2x_i % rAttr(index_i2x_Fioi_dmspd ,n) = oceanDMSPdFlux(i)
               i2x_i % rAttr(index_i2x_Fioi_docr  ,n) = oceanHumicsFlux(i)
               ! convert from umol Fe/m^3 to mmol Fe/m^3
               i2x_i % rAttr(index_i2x_Fioi_fep1  ,n) = oceanParticulateIronFlux(1,i) / 1000._RKIND
               i2x_i % rAttr(index_i2x_Fioi_fep2  ,n) = oceanParticulateIronFlux(2,i) / 1000._RKIND
               i2x_i % rAttr(index_i2x_Fioi_fed1  ,n) = oceanDissolvedIronFlux(1,i) / 1000._RKIND
               i2x_i % rAttr(index_i2x_Fioi_fed2  ,n) = oceanDissolvedIronFlux(2,i) / 1000._RKIND
            endif
         endif
      enddo

      block_ptr => block_ptr % next
   enddo

!-----------------------------------------------------------------------
!EOC

 end subroutine ice_export_mct!}}}

   subroutine add_stream_attributes(domain)!{{{

      type (domain_type), intent(inout) :: domain

      type (MPAS_Pool_iterator_type) :: itr
      integer, pointer :: intAtt
      logical, pointer :: logAtt
      character (len=StrKIND), pointer :: charAtt
      real (kind=RKIND), pointer :: realAtt
      character (len=StrKIND) :: histAtt

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
     

      call MPAS_stream_mgr_add_att(domain % streamManager, 'on_a_sphere', domain % on_a_sphere)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'sphere_radius', domain % sphere_radius)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'model_name', domain % core % modelName)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'core_name', domain % core % coreName)
      ! DWJ 10/01/2014: Eventually add the real history attribute, for now (due to length restrictions)
      ! add a shortened version.
!     call MPAS_stream_mgr_add_att(domain % streamManager, 'history', domain % history)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'history', histAtt)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'source', domain % core % source)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'Conventions', domain % core % Conventions)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'parent_id', domain % parent_id)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'mesh_spec', domain % mesh_spec)
      call MPAS_stream_mgr_add_att(domain % streamManager, 'git_version', domain % core % git_version)

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

!***********************************************************************
!BOP
!
! !IROUTINE: basal_pressure
!
! !INTERFACE:
   subroutine basal_pressure(basalPressure,  iceVolumeCell,  snowVolumeCell,  pondDepthCell, &
                             pondLidThicknessCell,  pondAreaCell, config_use_topo_meltponds)!{{{
!
! !DESCRIPTION:
! Calculate basal pressure for a cell
!

! !INPUT PARAMETERS:
      real (kind=RKIND), intent(in)  :: iceVolumeCell
      real (kind=RKIND), intent(in)  :: snowVolumeCell
      real (kind=RKIND), intent(in)  :: pondDepthCell
      real (kind=RKIND), intent(in)  :: pondLidThicknessCell
      real (kind=RKIND), intent(in)  :: pondAreaCell
      logical,           intent(in)  :: config_use_topo_meltponds

! !OUTPUT PARAMETERS:
      real (kind=RKIND), intent(out) :: basalPressure

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   real(kind=RKIND) :: &
      seaIceSpecificMass

   ! sea ice and snow specific mass
   seaIceSpecificMass = &
      iceVolumeCell  * seaiceDensityIce + &
      snowVolumeCell * seaiceDensitySnow

   ! only topo ponds have real pond volume
   if (config_use_topo_meltponds) then

      ! add pond specific weight
      seaIceSpecificMass = seaIceSpecificMass + &
         pondDepthCell        * pondAreaCell * seaiceDensitySeaWater + &
         pondLidThicknessCell * pondAreaCell * seaiceDensityIce

   endif ! config_use_topo_meltponds

   ! convert specific mass to pressure at sea ice base
   basalPressure = seaIceSpecificMass * seaiceGravity

!-----------------------------------------------------------------------
!EOC

   end subroutine basal_pressure!}}}

!***********************************************************************
!BOP
!
! !IROUTINE: freezing_potential
!
! !INTERFACE:
   subroutine freezing_potential(freezingPotential, frazilProduction, seaSurfaceSalinity, &
                                 config_thermodynamics_type)
!
! !DESCRIPTION:
! Calculate freezing potential for a cell
!
! !USES:
      use ice_mushy_physics, only:  &
         liquidus_temperature_mush, &
         enthalpy_mush

      use ice_colpkg_shared, only: &
         dSin0_frazil,             &
         phi_init

! !INPUT PARAMETERS:
      real (kind=RKIND),      intent(in)  :: frazilProduction
      real (kind=RKIND),      intent(in)  :: seaSurfaceSalinity
      character(len=strKIND), intent(in)  :: config_thermodynamics_type

! !OUTPUT PARAMETERS:
      real (kind=RKIND),      intent(out) :: freezingPotential

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
   real(kind=RKIND) :: &
      Si0new,          &
      Ti,              &
      qi0new,          &
      vi0new


   if (frazilProduction > 0.0_RKIND) then

      if (trim(config_thermodynamics_type) == "mushy") then  ! mushy
         if (seaSurfaceSalinity > 2.0_RKIND * dSin0_frazil) then
             Si0new = seaSurfaceSalinity - dSin0_frazil
         else
             Si0new = seaSurfaceSalinity**2 / (4.0_RKIND*dSin0_frazil)
         endif
         Ti = min(liquidus_temperature_mush(Si0new/phi_init), -0.1_RKIND)
         qi0new = enthalpy_mush(Ti, Si0new)
      else
         qi0new = -seaiceDensityIce*seaiceLatentHeatMelting
      endif    ! ktherm

      vi0new = max (frazilProduction, 0.0_RKIND) / seaiceDensityIce
      freezingPotential = -vi0new*qi0new ! note sign convention, qi < 0

   else

      freezingPotential = 0.0_RKIND

   endif

!-----------------------------------------------------------------------
!EOC

   end subroutine freezing_potential!}}}

end module ice_comp_mct

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
