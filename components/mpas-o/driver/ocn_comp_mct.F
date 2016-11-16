module ocn_comp_mct

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
! !MODULE: ocn_comp_mct
! !INTERFACE:

! !DESCRIPTION:
!  This is the main driver for the Model for Predication Across Scales Ocean Model (MPAS-O).
!
! !USES:
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
   use ocn_surface_land_ice_fluxes
   use ocn_diagnostics

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

  logical :: ldiag_cpl = .false.

  integer, private ::   &
      cpl_write_restart,   &! flag id for write restart
      cpl_write_history,   &! flag id for write history
      cpl_write_tavg,      &! flag id for write tavg      
      cpl_diag_global,     &! flag id for computing diagnostics
      cpl_diag_transp,     &! flag id for computing diagnostics
      my_task

   real (kind=RKIND) :: tlast_coupled

   integer  :: nsend, nrecv

   character(len=StrKIND) :: runtype, coupleTimeStamp

   type(seq_infodata_type), pointer :: infodata   
   type (iosystem_desc_t), pointer :: io_system 

   !! MPAS-O Datatypes
   type (core_type), save, pointer :: corelist => null()
   type (domain_type), save, pointer:: domain
   integer :: itimestep, ocn_cpl_dt
   integer :: ocnStdErrUnit, ocnStdOutUnit

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
! Initialize MPAS-O 
!
! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock), intent(inout) :: EClock
    type(seq_cdata), intent(inout) :: cdata_o
    type(mct_aVect), intent(inout) :: x2o_o, o2x_o
    character(len=*), optional, intent(in) :: NLFilename ! Namelist filename
!
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
       start_month, start_hour, iyear, mpas_o_cpl_dt, shrloglev, shrlogunit, stdout_shr, pio_iotype

    type(mct_gsMap), pointer :: gsMap_o

    type(mct_gGrid), pointer :: dom_o

    integer :: errorCode  ! error code

    character(len=StrKIND) :: cpl_seq_option

    type (MPAS_Time_Type) :: currTime
    integer :: iam, ierr, ierr_local
    integer :: iyear0, imonth0
    character(len=StrKIND)  :: starttype          ! infodata start type
    character(len=StrKIND)  :: timeStamp
    character(len=StrKIND)  :: nml_filename
    character(len=16) :: inst_suffix
    integer :: lbnum

    type (MPAS_Time_Type) :: alarmStartTime
    type (MPAS_TimeInterval_Type) :: alarmTimeStep
    type (block_type), pointer :: block_ptr

    type (mpas_pool_type), pointer :: meshPool, statePool, &
                                      forcingPool, diagnosticsPool, &
                                      averagePool, scratchPool

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
    type (MPAS_Time_type) :: start_time
    type (MPAS_Time_type) :: ref_time
    type (MPAS_TimeInterval_type) :: filename_interval
    type (MPAS_TimeInterval_type) :: denInterval, remInterval, zeroInterval
    integer :: numDivs
    character(len=StrKIND) :: start_timestamp
    character(len=StrKIND) :: iotype
    logical :: streamsExists
    integer :: mesh_iotype

    logical, pointer :: tempLogicalConfig
    character(len=StrKIND), pointer :: tempCharConfig

    real (kind=RKIND) :: dt
    type (MPAS_timeInterval_type) :: timeStep

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
!   setup mpas-o data structures
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

!-----------------------------------------------------------------------
!
!   initialize the model run 
!
!-----------------------------------------------------------------------

    call mpaso_cpl_indices_set()
    call seq_infodata_GetData( infodata, start_type=starttype)

    if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
       runtype = "initial"
    else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
       runtype = "continue"
    else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
       runtype = "branch"
    else
       write(stderrUnit, *) 'Start type is:', trim(starttype)
       call mpas_dmpar_global_abort(' ocn_comp_mct ERROR: unknown starttype')
    end if
   
!-----------------------------------------------------------------------
!
!   first initializaiton phase of mpas-o
!   initialize mpas-o because grid information is needed for
!   creation of GSMap_ocn.
!   call mpas-o initialization routines
!
!-----------------------------------------------------------------------
    inst_suffix = seq_comm_suffix(OCNID)

    call t_startf('mpas-o_init1')

    io_system => shr_pio_getiosys(ocnid)

    pio_iotype = shr_pio_getiotype(ocnid)
    call MPAS_io_set_iotype(domain % iocontext, pio_iotype)

    stdout_shr = shr_file_getUnit()
#ifndef MPAS_DEBUG
    if ( iam /= 0 ) then
       open(stdout_shr, status='SCRATCH')
    else
       inquire(file='ocn_modelio.nml'//trim(inst_suffix),exist=exists)
       if(exists) then
          call shr_file_setio('ocn_modelio.nml'//trim(inst_suffix),stdout_shr)
       endif
    end if
#else
    inquire(file='ocn_modelio.nml'//trim(inst_suffix),exist=exists)
    if(exists) then
       call shr_file_setio('ocn_modelio.nml'//trim(inst_suffix),stdout_shr)
    endif
#endif

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (stdout_shr)

    domain_ptr => domain
    ocnStdErrUnit = stdout_shr
    ocnStdOutUnit = stdout_shr
    call mpas_framework_init_phase1(domain_ptr % dminfo, mpicom_o, stdoutUnit_in=stdout_shr, stderrUnit_in=stdout_shr)

    ! Setup function pointers for MPAS-O core and domain types
    call ocn_setup_core(corelist)
    call ocn_setup_domain(domain_ptr)

    ! Override the names of the stream and namelist files
    domain_ptr % namelist_filename = 'mpas-o_in'
    domain_ptr % streams_filename = 'streams.ocean'

    ! Setup namelist variables, and read the namelist
    ierr = domain_ptr % core % setup_namelist(domain_ptr % configs, domain_ptr % namelist_filename, domain_ptr % dminfo)
    if ( ierr /= 0 ) then
       call mpas_dmpar_global_abort('ERROR: Namelist setup failed for core ' // trim(domain_ptr % core % coreName))
    end if

    call mpas_framework_init_phase2(domain_ptr, io_system)

    ! Define package variables
    ierr = domain_ptr % core % define_packages(domain_ptr % packages)
    if ( ierr /= 0 ) then
       call mpas_dmpar_global_abort('ERROR: Package definition failed for core ' // trim(domain_ptr % core % coreName))
    end if

    ! Setup packages (i.e. determine if they should be on or off)
    ierr = domain_ptr % core % setup_packages(domain_ptr % configs, domain_ptr % packages, domain_ptr % ioContext)
    if ( ierr /= 0 ) then
       call mpas_dmpar_global_abort('ERROR: Package setup failed for core ' // trim(domain_ptr % core % coreName))
    end if

    ! Setup decompositions available for dimensions
    ierr = domain_ptr % core % setup_decompositions(domain_ptr % decompositions)
    if ( ierr /= 0 ) then
       call mpas_dmpar_global_abort('ERROR: Decomposition setup failed for core ' // trim(domain_ptr % core % coreName))
    end if

    ! Override namelist options based on start type
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

    ! Setup MPAS-O simulation clock
    ierr = domain_ptr % core % setup_clock(domain_ptr % clock, domain_ptr % configs)
    if ( ierr /= 0 ) then
       call mpas_dmpar_global_abort('ERROR: Clock setup failed for core ' // trim(domain_ptr % core % coreName))
    end if

    write(stderrUnit,*) 'Reading streams configuration from file '//trim(domain_ptr % streams_filename)
    inquire(file=trim(domain_ptr % streams_filename), exist=streamsExists)

    if ( .not. streamsExists ) then
       call mpas_dmpar_global_abort('ERROR: Streams file '//trim(domain_ptr % streams_filename)//' does not exist.')
    end if

    !
    ! Using information from the namelist, a graph.info file, and a file containing
    !    mesh fields, build halos and allocate blocks in the domain
    !
    ierr = domain_ptr % core % get_mesh_stream(domain_ptr % configs, mesh_stream)
    if ( ierr /= 0 ) then
       call mpas_dmpar_global_abort('ERROR: Failed to find mesh stream for core ' // trim(domain_ptr % core % coreName))
    end if


    call mpas_f_to_c_string(domain % streams_filename, c_filename)
    call mpas_f_to_c_string(mesh_stream, c_mesh_stream)
    c_comm = domain_ptr % dminfo % comm
    call xml_stream_get_attributes(c_filename, c_mesh_stream, c_comm, &
                                   c_mesh_filename_temp, c_ref_time_temp, &
                                   c_filename_interval_temp, c_iotype, c_ierr)
    if (c_ierr /= 0) then
       call mpas_dmpar_abort(domain_ptr % dminfo)
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
    write(stderrUnit, *) ' ** Attempting to bootstrap MPAS framework using stream: ', trim(mesh_stream)
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
       call mpas_dmpar_global_abort('ERROR: Immutable streams setup failed for core ' // trim(domain_ptr % core % coreName))
    end if

    ! Parse / read all streams configuration
    mgr_p = c_loc(domain_ptr % streamManager)
    call xml_stream_parser(c_filename, mgr_p, c_comm, c_ierr)
    if (c_ierr /= 0) then
       call mpas_dmpar_abort(domain_ptr % dminfo)
    end if

    my_task = domain_ptr % dminfo % my_proc_id

    !
    ! Finalize the setup of blocks and fields
    !
    call mpas_bootstrap_framework_phase2(domain_ptr)
    call t_stopf('mpas-o_init1')

    call t_startf('mpas-o_init2')
    ! Initialize the MPAS-O core
    ierr = domain_ptr % core % core_init(domain_ptr, timeStamp)
    if ( ierr /= 0 ) then
       call mpas_dmpar_global_abort('ERROR: Core init failed for core ' // trim(domain_ptr % core % coreName))
    end if
    call t_stopf('mpas-o_init2')

!-----------------------------------------------------------------------
!
!   initialize time-stamp information
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!   check for consistency of mpas-o and sync clock initial time
!
!-----------------------------------------------------------------------

    if (runtype == 'initial') then
       call seq_timemgr_EClockGetData(EClock, ECurrTime=currTime % t)
       call mpas_set_clock_time(domain_ptr % clock, currTime, MPAS_START_TIME, ierr)
       call mpas_set_clock_time(domain_ptr % clock, currTime, MPAS_NOW, ierr)

       currTime = mpas_get_clock_time(domain_ptr % clock, MPAS_NOW, ierr)
       call mpas_get_time(curr_time=currTime, YYYY=iyear0, MM=imonth0, ierr=ierr)
       call seq_timemgr_EClockGetData(EClock, start_ymd=start_ymd, start_tod=start_tod)
       call shr_cal_date2ymd(start_ymd,start_year,start_month,start_day)

       if (iyear0 /= start_year) then
          write(shrlogunit, *) ' Iyear0 is: ', iyear0
          write(shrlogunit,*) ' start_year is: ', start_year
          call mpas_dmpar_global_abort(' iyear0 does not match start_year')
       end if
       if (imonth0 /= start_month) then
          write(shrlogunit, *) ' Imonth0 is: ', imonth0
          write(shrlogunit, *) ' start_month is: ', start_month
          call mpas_dmpar_global_abort(' imonth0 does not match start_month')
       end if
    else if (runtype == 'continue' .or. runtype == 'branch') then
       call seq_timemgr_EClockGetData(EClock, ECurrTime=currTime % t)
       call mpas_set_clock_time(domain_ptr % clock, currTime, MPAS_START_TIME, ierr)
       call mpas_set_clock_time(domain_ptr % clock, currTime, MPAS_NOW, ierr)
    end if

!-----------------------------------------------------------------------
!
!   initialize MCT attribute vectors and indices
!
!-----------------------------------------------------------------------

    call t_startf ('mpas-o_mct_init')

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
       write(stderrUnit, *) ' ERROR: Coupling interval is: ', trim(coupleTimeStamp)
       write(stderrUnit, *) '        OCN Model time step is: ', trim(tempCharConfig)

       write(stderrUnit, *) '        The model time step cannot be longer then the coupling interval'
       call mpas_dmpar_global_abort('ERROR: Model is not properly configured for coupling interval.')
    end if

    if ( remInterval > zeroInterval ) then
       ierr = 1
    end if

    ierr_local = ierr
    call mpas_dmpar_max_int(domain % dminfo, ierr_local, ierr)

    if ( ierr == 1 ) then
       write(stderrUnit, *) ' ERROR: Coupling interval is: ', trim(coupleTimeStamp)
       write(stderrUnit, *) '        OCN Model time step is: ', trim(tempCharConfig)

       write(stderrUnit, *) '        These are not synchronized, so time steps'
       write(stderrUnit, *) '        will not match to coupling interval boundaries.'
       write(stderrUnit, *) '        Please reconfigure either the coupling interval, '
       write(stderrUnit, *) '        or the time step.'
       call mpas_dmpar_global_abort('ERROR: Model is not properly configured for coupling interval.')
    end if

    ! set coupling alarm
    alarmStartTime = currTime
    call mpas_add_clock_alarm(domain % clock, coupleAlarmID, alarmStartTime, alarmTimeStep, ierr=ierr)
    call mpas_reset_clock_alarm(domain_ptr % clock, coupleAlarmID, ierr=ierr)

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
          call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)

          call ocn_time_average_coupled_init(forcingPool)
          call ocn_time_average_coupled_accumulate(diagnosticsPool, statePool, forcingPool, 1)
          block_ptr => block_ptr % next
       end do
    end if

!-----------------------------------------------------------------------
!
!   send intial state to driver
!
!-----------------------------------------------------------------------

    call ocn_export_mct(o2x_o, errorCode)  
    if (errorCode /= 0) then
       call mpas_dmpar_global_abort('ERROR in ocn_export_mct')
    endif

    ! Setup clock for initial runs
    if ( runtype == 'continue' ) then
       block_ptr => domain % blocklist
       do while(associated(block_ptr))
          call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
          call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)

          call ocn_time_average_coupled_init(forcingPool)
          call ocn_time_average_coupled_accumulate(diagnosticsPool, statePool, forcingPool, 1)
          block_ptr => block_ptr % next
       end do
    end if

    call t_stopf ('mpas-o_mct_init')

    call seq_infodata_PutData( infodata, ocn_prognostic=.true., ocnrof_prognostic=.true.)

!-----------------------------------------------------------------------
!
!   get intial state from driver
!
!-----------------------------------------------------------------------

    timeStep = mpas_get_clock_timestep(domain_ptr % clock, ierr=ierr)
    call mpas_get_timeInterval(timeStep, dt=dt)

    call ocn_import_mct(x2o_o, errorCode)  
    if (errorCode /= 0) then
       call mpas_dmpar_global_abort('ERROR in ocn_import_mct')
    endif

    ! Build forcing arrays.
    block_ptr => domain % blocklist
    do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
        call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)

        call ocn_forcing_build_fraction_absorbed_array(meshPool, statePool, diagnosticsPool, forcingPool, ierr, 1)
        call mpas_timer_start("land_ice_build_arrays", .false.)
        call ocn_surface_land_ice_fluxes_build_arrays(meshPool, diagnosticsPool, &
                                                      forcingPool, scratchPool, statePool, dt, ierr)
        call mpas_timer_stop("land_ice_build_arrays")

        call ocn_frazil_forcing_build_arrays(domain, meshPool, forcingPool, diagnosticsPool, statePool, ierr)

        block_ptr => block_ptr % next
    end do

    currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
    call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
    write(stderrUnit,*) 'Initial time ', trim(timeStamp)

    itimestep = 0

    ! Reset all output alarms, to prevent intial time step from writing any output, unless it's ringing.
    call mpas_stream_mgr_reset_alarms(domain_ptr % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)
    call mpas_stream_mgr_reset_alarms(domain_ptr % streamManager, direction=MPAS_STREAM_INPUT, ierr=ierr)

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
        call memmon_dump_fort('memmon.out','ocn_init_mct:end::',lbnum) 
        call memmon_reset_addr()
    endif
#endif
    call shr_sys_flush(stderrUnit)

#ifndef MPAS_DEBUG
    if ( iam /= 0 ) then
       close(stdout_shr)
    end if
#endif


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
! Run MPAS-O for one coupling interval
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
      real (kind=RKIND) :: dt
      type (block_type), pointer :: block_ptr

      type (mpas_pool_type), pointer :: meshPool, statePool, diagnosticsPool, forcingPool, averagePool, scratchPool

      type (MPAS_Time_Type) :: currTime
      type (domain_type), pointer :: domain_ptr
      character(len=StrKIND) :: timeStamp, streamName
      type (MPAS_timeInterval_type) :: timeStep
      integer :: ierr, streamDirection
      logical :: streamActive
      logical, pointer :: config_write_output_on_startup
      character (len=StrKIND), pointer :: config_restart_timestamp_name

#ifndef MPAS_DEBUG
    if ( domain % dminfo % my_proc_id /= 0 ) then
       open(ocnStdErrUnit, status='SCRATCH')
    end if
#endif


      call mpas_io_units_set_stdout(ocnStdOutUnit)
      call mpas_io_units_set_stderr(ocnStdErrUnit)

      domain_ptr => domain

      call mpas_pool_get_config(domain_ptr % configs, 'config_restart_timestamp_name', config_restart_timestamp_name)

      ! Setup log information.
      call shr_file_getLogUnit (shrlogunit)
      call shr_file_getLogLevel(shrloglev)
      call shr_file_setLogUnit (stdoutUnit)

      timeStep = mpas_get_clock_timestep(domain_ptr % clock, ierr=ierr)
      call mpas_get_timeInterval(timeStep, dt=dt)
      call mpas_reset_clock_alarm(domain_ptr % clock, coupleAlarmID, ierr=ierr)

      ! Import state from coupler
      call ocn_import_mct(x2o_o, ierr)

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

         itimestep = itimestep + 1
         call mpas_advance_clock(domain_ptr % clock)

         currTime = mpas_get_clock_time(domain_ptr % clock, MPAS_NOW, ierr)
         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
         write(stderrUnit,*) 'Doing timestep ', trim(timeStamp)
         call shr_sys_flush(stderrUnit)

         ! Build forcing arrays.
         write(stderrUnit,*) '   Building forcing arrays'
         call shr_sys_flush(stderrUnit)
         block_ptr => domain_ptr % blocklist
         do while(associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
            call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
            call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
            call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
            call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)

            call ocn_forcing_build_fraction_absorbed_array(meshPool, statePool, diagnosticsPool, forcingPool, ierr, 1)

            call mpas_timer_start("land_ice_build_arrays", .false.)
            call ocn_surface_land_ice_fluxes_build_arrays(meshPool, diagnosticsPool, &
                                                          forcingPool, scratchPool, statePool, dt, ierr)
            call mpas_timer_stop("land_ice_build_arrays")

            call ocn_frazil_forcing_build_arrays(domain, meshPool, forcingPool, diagnosticsPool, statePool, ierr)
            block_ptr => block_ptr % next
         end do
         write(stderrUnit,*) '   Finished building forcing arrays'
         call shr_sys_flush(stderrUnit)

         write(stderrUnit,*) '   Performing forward update'
         call shr_sys_flush(stderrUnit)
         call mpas_timer_start("time integration", .false.)
         !$omp parallel default(firstprivate) shared(domain, dt, timeStamp)
         call ocn_timestep(domain, dt, timeStamp)
         !$omp end parallel
         call mpas_timer_stop("time integration")
         write(stderrUnit,*) '   Finished performing forward update'
         call shr_sys_flush(stderrUnit)

         ! Move time level 2 fields back into time level 1 for next time step
         block_ptr => domain_ptr % blocklist
         do while(associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
            call mpas_pool_shift_time_levels(statePool)
            block_ptr => block_ptr % next
         end do

         write(stderrUnit,*) '   Computing analysis members'
         call shr_sys_flush(stderrUnit)
         call ocn_analysis_compute(domain_ptr, ierr) 
         write(stderrUnit,*) '   Finished computing analysis members'
         write(stderrUnit,*) '   Preparing analysis members for restart'
         call shr_sys_flush(stderrUnit)
         call ocn_analysis_restart(domain_ptr, ierr) 
         write(stderrUnit,*) '   Finished preparing analysis members for restart'
         write(stderrUnit,*) '   Writing analysis member output'
         call shr_sys_flush(stderrUnit)
         call ocn_analysis_write(domain_ptr, ierr)
         write(stderrUnit,*) '   Finished writing analysis member output'
         call shr_sys_flush(stderrUnit)

         ! Reset any restart alarms to prevent restart files being written without the coupler requesting it.
         write(stderrUnit,*) '   Resetting restart alarms'
         call shr_sys_flush(stderrUnit)
         call mpas_stream_mgr_begin_iteration(domain % streamManager)

         do while ( mpas_stream_mgr_get_next_stream( domain % streamManager, streamID=streamName, &
                    directionProperty=streamDirection, activeProperty=streamActive ) )

            if ( streamActive .and. streamDirection == MPAS_STREAM_INPUT_OUTPUT ) then
               call mpas_stream_mgr_reset_alarms(domain % streamManager, streamID=streamName, ierr=ierr)
            end if
         end do
         write(stderrUnit,*) '   Finished resetting restart alarms'
         call shr_sys_flush(stderrUnit)

         write(stderrUnit,*) '   Writing output streams'
         call shr_sys_flush(stderrUnit)
         call mpas_stream_mgr_write(domain % streamManager, streamID='output', ierr=ierr)
         call mpas_stream_mgr_reset_alarms(domain % streamManager, streamID='output', ierr=ierr)

         call mpas_stream_mgr_write(domain % streamManager, ierr=ierr)
         call mpas_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)
         write(stderrUnit,*) '   Finished writing output streams'
         call shr_sys_flush(stderrUnit)
         call ocn_validate_state(domain, timeLevel=1)
         write(stderrUnit,*) '   Validating ocean state'
         write(stderrUnit,*) '   Completed validating ocean state'
         write(stderrUnit,*) 'Completed timestep ', trim(timeStamp)

         call shr_sys_flush(stderrUnit)
      end do

      ! Check if coupler wants us to write a restart file.
      ! We only write restart files at the end of a coupling interval
      if (seq_timemgr_RestartAlarmIsOn(EClock)) then
         ! Write a restart file, because the coupler asked for it.
         write(stderrUnit,*) 'Writing restart streams'
         call shr_sys_flush(stderrUnit)
         call mpas_stream_mgr_begin_iteration(domain % streamManager)

         do while ( mpas_stream_mgr_get_next_stream( domain % streamManager, streamID=streamName, &
                    directionProperty=streamDirection, activeProperty=streamActive ) )

            if ( streamActive .and. streamDirection == MPAS_STREAM_INPUT_OUTPUT ) then
               write(stderrUnit,*) '   Writing stream ' // trim(streamName)
               call shr_sys_flush(stderrUnit)
               call mpas_stream_mgr_write(domain % streamManager, forceWriteNow=.true., streamID=streamName, ierr=ierr)
               write(stderrUnit,*) '   Finished writing stream ' // trim(streamName)
               call shr_sys_flush(stderrUnit)
            end if
         end do

         if ( domain_ptr % dminfo % my_proc_id == 0 ) then
            open(22, file=config_restart_timestamp_name, form='formatted', status='replace')
            write(22, *) trim(timeStamp)
            close(22)
         end if
         write(stderrUnit,*) 'Finished writing restart streams'
         call shr_sys_flush(stderrUnit)
      end if

      ! Export state to coupler
      write(stderrUnit,*) 'Exporting ocean state'
      call shr_sys_flush(stderrUnit)
      call ocn_export_mct(o2x_o, ierr)
      write(stderrUnit,*) 'Finished exporting ocean state'
      call shr_sys_flush(stderrUnit)

      ! Check if clocks are in sync
      currTime = mpas_get_clock_time(domain_ptr % clock, MPAS_NOW, ierr)
      call mpas_get_time(curr_time=currTime, YYYY=iyear, MM=imonth, DD=iday, H=ihour, M=iminute, S=isecond, ierr=ierr)
      call seq_timemgr_EClockGetData(EClock, curr_ymd=curr_ymd, curr_tod=curr_tod)
 
      ymd = iyear * 10000 + imonth * 100 + iday
      tod = ihour * 3600 + iminute * 60 + isecond
      if (.not. seq_timemgr_EClockDateInSync( EClock, ymd, tod)) then
         write(stderrUnit,*) 'MPAS ymd=',ymd,' MPAS tod=', tod
         write(stderrUnit,*) 'sync ymd=',curr_ymd,' sync tod=', curr_tod
         write(stderrUnit,*) 'Internal mpas clock not in sync with sync clock'
      end if

      ! Reset I/O logs
      call shr_file_setLogUnit (shrlogunit)
      call shr_file_setLogLevel(shrloglev)

#ifndef MPAS_DEBUG
    if ( domain % dminfo % my_proc_id /= 0 ) then
       close(ocnStdErrUnit)
    end if
#endif


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
! Finalize MPAS-O
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
    integer :: errorCode

    type (domain_type), pointer :: domain_ptr

    integer :: ierr, my_proc_id

!-----------------------------------------------------------------------
    my_proc_id = domain % dminfo % my_proc_id
#ifndef MPAS_DEBUG
    if ( my_proc_id /= 0 ) then
       open(ocnStdErrUnit, status='SCRATCH')
    end if
#endif

    call mpas_io_units_set_stdout(ocnStdOutUnit)
    call mpas_io_units_set_stderr(ocnStdErrUnit)

    domain_ptr => domain

    ! Setup I/O logs
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (stdoutUnit)

    ! Finalize MPAS-O
    iErr = domain_ptr % core % core_finalize(domain_ptr)
    if ( iErr /= 0 ) then
       call mpas_dmpar_global_abort('ERROR: Core finalize failed for core ' // trim(domain_ptr % core % coreName))
    end if

    call mpas_timer_write()
    call MPAS_stream_mgr_finalize(domain_ptr % streamManager)
    call mpas_framework_finalize(domain_ptr % dminfo, domain_ptr, io_system)

    ! Reset I/O logs
    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

#ifndef MPAS_DEBUG
    if ( my_proc_id /= 0 ) then
       close(ocnStdErrUnit)
    end if
#endif

  end subroutine ocn_final_mct!}}}

!***********************************************************************
!BOP
!IROUTINE: ocn_SetGSMap_mct
! !INTERFACE:

 subroutine ocn_SetGSMap_mct( mpicom_ocn, OCNID, gsMap_ocn )!{{{

! !DESCRIPTION:
!  This routine sets up the MPAS-O grid numbering for MCT
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
!  This routine sets up the MCT domain for MPAS-O
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
! 
!    The following fields are sometimes received from the coupler,
!      depending on model options:
! 
!    o  pbot   -- bottom atm pressure                      (Pa)
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

   logical, pointer :: config_use_ecosysTracers, &
                       config_use_ecosysTracers_sea_ice_coupling, &
                       config_use_DMSTracers,    &
                       config_use_DMSTracers_sea_ice_coupling,  &
                       config_use_MacroMoleculesTracers,  &
                       config_use_MacroMoleculesTracers_sea_ice_coupling

   real (kind=RKIND), pointer :: config_density0

   type (block_type), pointer :: block_ptr

   type (mpas_pool_type), pointer :: meshPool,             &
                                     forcingPool,          &
                                     ecosysAuxiliary,      &
                                     ecosysSeaIceCoupling, &
                                     DMSSeaIceCoupling,    &
                                     MacroMoleculesSeaIceCoupling

   integer, pointer :: nCellsSolve

   type (field1DReal), pointer :: windStressZonalField, windStressMeridionalField, &
                                  latentHeatFluxField, sensibleHeatFluxField, &
                                  longWaveHeatFluxUpField, &
                                  longWaveHeatFluxDownField, &
                                  evaporationFluxField, seaIceHeatFluxField, &
                                  snowFluxField, seaIceFreshWaterFluxField, &
                                  seaIceSalinityFluxField, &
                                  riverRunoffFluxField, iceRunoffFluxField, &
                                  shortWaveHeatFluxField, rainFluxField, &
                                  seaSurfacePressureField, iceFractionField, &
                                  windSpeedSquared10mField, &
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
                                  iceFluxDustField

   type (field2DReal), pointer :: iceFluxPhytoCField, &
                                  iceFluxDOCField

   real (kind=RKIND), dimension(:), pointer :: windStressZonal, windStressMeridional, &
                                               latentHeatFlux, sensibleHeatFlux, &
                                               longWaveHeatFluxUp, &
                                               longWaveHeatFluxDown, &
                                               evaporationFlux, seaIceHeatFlux, &
                                               snowFlux, seaIceFreshWaterFlux, &
                                               seaIceSalinityFlux, &
                                               riverRunoffFlux, iceRunoffFlux, &
                                               shortWaveHeatFlux, rainFlux, &
                                               seaSurfacePressure, iceFraction, &
                                               windSpeedSquared10m, &
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
                                               iceFluxDust
           

   real (kind=RKIND), dimension(:,:), pointer :: iceFluxPhytoC, &
                                                 iceFluxDOC
   
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

   n = 0
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
      call mpas_pool_get_field(forcingPool, 'snowFlux', snowFluxField)
      call mpas_pool_get_field(forcingPool, 'seaIceFreshWaterFlux', seaIceFreshWaterFluxField)
      call mpas_pool_get_field(forcingPool, 'seaIceSalinityFlux', seaIceSalinityFluxField)
      call mpas_pool_get_field(forcingPool, 'riverRunoffFlux', riverRunoffFluxField)
      call mpas_pool_get_field(forcingPool, 'iceRunoffFlux', iceRunoffFluxField)
      call mpas_pool_get_field(forcingPool, 'shortWaveHeatFlux', shortWaveHeatFluxField)
      call mpas_pool_get_field(forcingPool, 'rainFlux', rainFluxField)
      call mpas_pool_get_field(forcingPool, 'seaSurfacePressure', seaSurfacePressureField)
      call mpas_pool_get_field(forcingPool, 'iceFraction', iceFractionField)
      call mpas_pool_get_field(forcingPool, 'iceRunoffFlux', iceRunoffFluxField)

      windStressZonal => windStressZonalField % array
      windStressMeridional => windStressMeridionalField % array
      latentHeatFlux => latentHeatFluxField % array
      sensibleHeatFlux => sensibleHeatFluxField % array
      longWaveHeatFluxUp => longWaveHeatFluxUpField % array
      longWaveHeatFluxDown => longWaveHeatFluxDownField % array
      evaporationFlux => evaporationFluxField % array
      seaIceHeatFlux => seaIceHeatFluxField % array
      snowFlux => snowFluxField % array
      seaIceFreshWaterFlux => seaIceFreshWaterFluxField % array
      seaIceSalinityFlux => seaIceSalinityFluxField % array
      riverRunoffFlux => riverRunoffFluxField % array
      iceRunoffFlux => iceRunoffFluxField % array
      shortWaveHeatFlux => shortWaveHeatFluxField % array
      rainFlux => rainFluxField % array
      seaSurfacePressure => seaSurfacePressureField % array
      iceFraction => iceFractionField % array
      iceRunoffFlux => iceRunoffFluxField % array

      ! BGC fields
      if (config_use_ecosysTracers) then
         call mpas_pool_get_subpool(forcingPool, 'ecosysAuxiliary', ecosysAuxiliary)

         call mpas_pool_get_field(ecosysAuxiliary, 'windSpeedSquared10m', windSpeedSquared10mField)
         windSpeedSquared10m => windSpeedSquared10mField % array

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
        if ( snowFluxField % isActive ) then
           snowFlux(i) = x2o_o % rAttr(index_x2o_Faxa_snow, n)
        end if

        if ( seaIceFreshWaterFluxField % isActive ) then
           seaIceFreshWaterFlux(i) = x2o_o % rAttr(index_x2o_Fioi_meltw, n)
        end if
        if ( seaIceSalinityFluxField % isActive ) then
           seaIceSalinityFlux(i) = x2o_o % rAttr(index_x2o_Fioi_salt, n)
        end if
        if ( riverRunoffFluxField % isActive ) then
           riverRunoffFlux(i) = x2o_o % rAttr(index_x2o_Foxx_rofl, n)
        end if
        if ( iceRunoffFluxField % isActive ) then
           iceRunoffFlux(i) = x2o_o % rAttr(index_x2o_Foxx_rofi, n)
           if(iceRunoffFlux(n) < 0.0_RKIND) then
               call shr_sys_abort ('Error: incoming rofi_F is negative')
           end if
        end if

        if ( shortWaveHeatFluxField % isActive ) then
           shortWaveHeatFlux(i) = max(x2o_o % rAttr(index_x2o_Foxx_swnet, n), 0.0_RKIND)
        end if

        if ( rainFluxField % isActive ) then
           rainFlux(i) = x2o_o % rAttr(index_x2o_Faxa_rain, n)
        end if
        if ( seaSurfacePressureField % isActive ) then
           ! Set seaSurfacePressure to be the sum of atmospheric bottom pressure 
           ! and a limited version of sea ice pressure where it limits it to 5m of pressure
           seaSurfacePressure(i) = x2o_o % rAttr(index_x2o_Sa_pbot,   n) &
                                 + min( x2o_o % rAttr(index_x2o_Si_bpress, n), config_density0 * gravity * 5.0_RKIND )
        end if
        if ( iceFractionField % isActive ) then
           iceFraction(i) = x2o_o % rAttr(index_x2o_Si_ifrac, n)
        end if

        ! BGC fields
        if (config_use_ecosysTracers) then
           if ( windSpeedSquared10mField % isActive ) then
              windSpeedSquared10m(i) = x2o_o % rAttr(index_x2o_So_duu10n, n)
           end if

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
   call mpas_pool_get_field(forcingPool, 'snowFlux', snowFluxField)
   call mpas_pool_get_field(forcingPool, 'seaIceFreshWaterFlux', seaIceFreshWaterFluxField)
   call mpas_pool_get_field(forcingPool, 'seaIceSalinityFlux', seaIceSalinityFluxField)
   call mpas_pool_get_field(forcingPool, 'riverRunoffFlux', riverRunoffFluxField)
   call mpas_pool_get_field(forcingPool, 'iceRunoffFlux', iceRunoffFluxField)
   call mpas_pool_get_field(forcingPool, 'shortWaveHeatFlux', shortWaveHeatFluxField)
   call mpas_pool_get_field(forcingPool, 'rainFlux', rainFluxField)
   call mpas_pool_get_field(forcingPool, 'seaSurfacePressure', seaSurfacePressureField)
   call mpas_pool_get_field(forcingPool, 'iceFraction', iceFractionField)
   call mpas_pool_get_field(forcingPool, 'iceRunoffFlux', iceRunoffFluxField)

   ! BGC fields
   if (config_use_ecosysTracers) then
      call mpas_pool_get_subpool(forcingPool, 'ecosysAuxiliary', ecosysAuxiliary)

      call mpas_pool_get_field(ecosysAuxiliary, 'windSpeedSquared10m', windSpeedSquared10mField)

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
   if ( snowFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(snowFluxField)
   end if
   if ( seaIceFreshWaterFluxField % isActive ) then
      call mpas_dmpar_exch_halo_field(seaIceFreshWaterFluxField)
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
   if ( seaSurfacePressureField % isActive ) then
      call mpas_dmpar_exch_halo_field(seaSurfacePressureField)
   end if
   if ( iceFractionField % isActive ) then
      call mpas_dmpar_exch_halo_field(iceFractionField)
   end if

   ! BGC fields
   if (config_use_ecosysTracers) then

      if ( windSpeedSquared10mField % isActive ) then
         call mpas_dmpar_exch_halo_field(windSpeedSquared10mField)
      end if

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

!-----------------------------------------------------------------------
!EOC

 end subroutine ocn_import_mct!}}}
!***********************************************************************
!BOP
! !IROUTINE: ocn_export_mct
! !INTERFACE:

 subroutine ocn_export_mct(o2x_o, errorCode)   !{{{

! !DESCRIPTION:
!  This routine calls the routines necessary to send MPAS-O fields to
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
                       index_avgZonalSurfaceVelocity, index_avgMeridionalSurfaceVelocity,      &
                       index_avgZonalSSHGradient, index_avgMeridionalSSHGradient

   type (block_type), pointer :: block_ptr

   type (mpas_pool_type), pointer :: meshPool,             &
                                     forcingPool,          &
                                     statePool,            &
                                     tracersPool,          &
                                     ecosysSeaIceCoupling, &
                                     DMSSeaIceCoupling,    &
                                     MacroMoleculesSeaIceCoupling

   integer, dimension(:), pointer :: landIceMask

   real (kind=RKIND), dimension(:), pointer :: seaIceEnergy, accumulatedFrazilIceMass, frazilSurfacePressure, &
                                               CO2Flux, DMSFlux, surfaceUpwardCO2Flux, &
                                               avgOceanSurfaceDIC, &
                                               avgOceanSurfaceDON, &
                                               avgOceanSurfaceNO3, &
                                               avgOceanSurfaceSiO3, &
                                               avgOceanSurfaceNH4, &
                                               avgOceanSurfaceDMS, &
                                               avgOceanSurfaceDMSP, &
                                               avgOceanSurfaceDOCr, &
                                               avgOceanSurfaceFeParticulate, &
                                               avgOceanSurfaceFeDissolved
 
   real (kind=RKIND), dimension(:,:), pointer :: avgTracersSurfaceValue, avgSurfaceVelocity, avgSSHGradient, &
                                                 avgOceanSurfacePhytoC, &
                                                 avgOceanSurfaceDOC, layerThickness

   logical, pointer :: frazilIceActive,          &
                       config_use_ecosysTracers, &
                       config_use_DMSTracers,    &
                       config_use_MacroMoleculesTracers,  &
                       config_use_ecosysTracers_sea_ice_coupling, &
                       config_use_DMSTracers_sea_ice_coupling,  &
                       config_use_MacroMoleculesTracers_sea_ice_coupling

   logical :: keepFrazil

   errorcode = 0

   ! get configure options
   call mpas_pool_get_package(domain % packages, 'frazilIceActive', frazilIceActive)
   call mpas_pool_get_config(domain % configs, 'config_use_ecosysTracers', config_use_ecosysTracers)
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

     call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)

     call mpas_pool_get_array(forcingPool, 'landIceMask', landIceMask)
     call mpas_pool_get_array(forcingPool, 'avgTracersSurfaceValue', avgTracersSurfaceValue)
     call mpas_pool_get_array(forcingPool, 'avgSurfaceVelocity', avgSurfaceVelocity)
     call mpas_pool_get_array(forcingPool, 'avgSSHGradient', avgSSHGradient)
     if ( frazilIceActive ) then
        call mpas_pool_get_array(forcingPool, 'seaIceEnergy', seaIceEnergy)
        call mpas_pool_get_array(forcingPool, 'frazilSurfacePressure', frazilSurfacePressure)
        call mpas_pool_get_array(statePool, 'accumulatedFrazilIceMass', accumulatedFrazilIceMass, 1)
     end if

     ! BGC fields
     if (config_use_ecosysTracers .and. config_use_ecosysTracers_sea_ice_coupling) then
        call mpas_pool_get_subpool(forcingPool, 'ecosysSeaIceCoupling', ecosysSeaIceCoupling)

        call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfacePhytoC', avgOceanSurfacePhytoC)
        call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceDIC', avgOceanSurfaceDIC)
        call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceNO3', avgOceanSurfaceNO3)
        call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceSiO3', avgOceanSurfaceSiO3)
        call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceNH4', avgOceanSurfaceNH4)
        call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceDOCr', avgOceanSurfaceDOCr)
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

       o2x_o % rAttr(index_o2x_So_dhdx, n) = avgSSHGradient(index_avgZonalSSHGradient, i)
       o2x_o % rAttr(index_o2x_So_dhdy, n) = avgSSHGradient(index_avgMeridionalSSHGradient, i)
       if ( frazilIceActive ) then
          ! negative when frazil ice can be melted
          keepFrazil = .true.
          if ( associated(landIceMask) ) then
             if ( landIceMask(i) == 1 ) then
                keepFrazil = .false.
             end if
          end if

          if ( keepFrazil ) then
             seaIceEnergy(i) = min( rho_sw * cp_sw * layerThickness(1, i) * ( -1.8_RKIND + T0_Kelvin &
                             - avgTracersSurfaceValue(index_temperatureSurfaceValue, i) ), 0.0_RKIND )
             if ( accumulatedFrazilIceMass(i) > 0.0_RKIND ) seaIceEnergy(i) = 0

             o2x_o % rAttr(index_o2x_Fioo_meltp, n)  = seaIceEnergy(i) / ocn_cpl_dt
             o2x_o % rAttr(index_o2x_Fioo_frazil, n) = accumulatedFrazilIceMass(i) / ocn_cpl_dt
          else
             o2x_o % rAttr(index_o2x_Fioo_meltp, n)  = 0.0_RKIND
             o2x_o % rAttr(index_o2x_Fioo_frazil, n) = 0.0_RKIND
          end if

          ! Reset SeaIce Energy and Accumulated Frazil Ice
          seaIceEnergy(i) = 0.0_RKIND
          accumulatedFrazilIceMass(i) = 0.0_RKIND
          frazilSurfacePressure(i) = 0.0_RKIND
       end if

       ! BGC fields
       if (config_use_ecosysTracers .and. config_use_ecosysTracers_sea_ice_coupling) then
          o2x_o % rAttr(index_o2x_So_algae1, n) = max(0.0_RKIND,avgOceanSurfacePhytoC(1,i))
          o2x_o % rAttr(index_o2x_So_algae2, n) = max(0.0_RKIND,avgOceanSurfacePhytoC(2,i))
          o2x_o % rAttr(index_o2x_So_algae3, n) = max(0.0_RKIND,avgOceanSurfacePhytoC(3,i))
          o2x_o % rAttr(index_o2x_So_dic1,   n) = max(0.0_RKIND,avgOceanSurfaceDIC(i))
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
          o2x_o % rAttr(index_o2x_So_don1, n) = max(0.0_RKIND,avgOceanSurfaceDON(i))
       endif
!      o2x_o % rAttr(index_o2x_Faoo_fco2_ocn, n) = CO2Flux(i)
!      o2x_o % rAttr(index_o2x_Faoo_fdms_ocn, n) = DMSFlux(i)
!      o2x_o % rAttr(index_o2x_Faoo_fco2_ocn, n) = surfaceUpwardCO2Flux(i)

     end do

     block_ptr => block_ptr % next
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine ocn_export_mct!}}}

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


end module ocn_comp_mct

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
