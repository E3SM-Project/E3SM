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
   use ocn_time_average
   use ocn_time_average_coupled
   use ocn_global_diagnostics

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
    integer :: iam,ierr 
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
                                      averagePool

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
    character(len=StrKIND) :: start_timestamp
    character(len=StrKIND) :: iotype
    logical :: streamsExists
    integer :: mesh_iotype

    logical, pointer :: tempLogicalConfig
    character(len=StrKIND), pointer :: tempCharConfig

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

    call t_startf('mpas-o_init')

    io_system => shr_pio_getiosys(ocnid)

    pio_iotype = shr_pio_getiotype(ocnid)
    call MPAS_io_set_iotype(pio_iotype)

    stdout_shr = shr_file_getUnit()
#ifndef MPAS_DEBUG
    if ( iam /= 0 ) then
       open(stdout_shr, file='/dev/null', position='APPEND')
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
    call mpas_framework_init_phase1(domain_ptr % dminfo, mpicom_o, stdoutUnit_in=stdout_shr, stderrUnit_in=stdout_shr)

    ! Setup function pointers for MPAS-O core and domain types
    call ocn_setup_core(corelist)
    call ocn_setup_domain(domain_ptr)

    ! Override the names of the stream and namelist files
    domain_ptr % namelist_filename = 'mpaso_in'
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
    ierr = domain_ptr % core % setup_packages(domain_ptr % configs, domain_ptr % packages)
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

        ! Setup restingThickness
        call mpas_pool_get_config(domain_ptr % configs, "config_set_restingThickness_to_IC", tempLogicalConfig)
        tempLogicalConfig = .true.

        ! Setup PBCs
        call mpas_pool_get_config(domain_ptr % configs, "config_alter_ICs_for_pbcs", tempLogicalConfig)
        tempLogicalConfig = .true.

        ! Setup start time. Will be over written later when clocks are synchronized
        call mpas_pool_get_config(domain_ptr % configs, "config_start_time", tempCharConfig)
        tempCharConfig = "0001-01-01_00:00:00"

        ! Setup run duration. Will be ignored in coupled run, since coupler defines how long the run is.
        call mpas_pool_get_config(domain_ptr % configs, "config_run_duration", tempCharConfig)
        tempCharConfig = "0001-00-00_00:00:00"
    else if (runtype == "continue" .or. runtype == "branch") then ! Restart run or branch run
        ! Turn on restart
        call mpas_pool_get_config(domain_ptr % configs, "config_do_restart", tempLogicalConfig)
        tempLogicalConfig = .true.

        ! Turn off restingThickness setup
        call mpas_pool_get_config(domain_ptr % configs, "config_set_restingThickness_to_IC", tempLogicalConfig)
        tempLogicalConfig = .false.

        ! Turn off PBC setup
        call mpas_pool_get_config(domain_ptr % configs, "config_alter_ICs_for_pbcs", tempLogicalConfig)
        tempLogicalConfig = .false.

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
        call mpas_expand_string(ref_time_temp, mesh_filename_temp, mesh_filename)
    else
        call mpas_set_time(ref_time, dateTimeString=ref_time_temp, ierr=ierr)
        call mpas_set_timeInterval(filename_interval, timeString=filename_interval_temp, ierr=ierr)
        call mpas_build_stream_filename(ref_time, start_time, filename_interval, mesh_filename_temp, mesh_filename, ierr)
    end if

    ! Bootstrap framework (1). Here data structures are setup, but dimensions and arrays are not finalized.
    write(stderrUnit, *) ' ** Attempting to bootstrap MPAS framework using stream: ', trim(mesh_stream)
    call mpas_bootstrap_framework_phase1(domain_ptr, mesh_filename, mesh_iotype)

    !
    ! Set up run-time streams
    !
    call MPAS_stream_mgr_init(domain_ptr % streamManager, domain_ptr % clock, domain_ptr % blocklist % allFields, domain_ptr % packages, domain_ptr % blocklist % allStructs)

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

    call mpas_timer_start("total time")

    ! Initialize the MPAS-O core
    ierr = domain_ptr % core % core_init(domain_ptr, timeStamp)
    if ( ierr /= 0 ) then
       call mpas_dmpar_global_abort('ERROR: Core init failed for core ' // trim(domain_ptr % core % coreName))
    end if

!-----------------------------------------------------------------------
!
!   initialize time-stamp information
!
!-----------------------------------------------------------------------
    call t_stopf ('mpas-o_init')

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
    call convert_seconds_to_timestamp(ocn_cpl_dt, coupleTimeStamp)

    ! set coupling alarm
    call mpas_set_timeInterval(alarmTimeStep, timeString=coupleTimeStamp, ierr=ierr)
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
          call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)

          call ocn_time_average_coupled_init(forcingPool)
          call ocn_time_average_coupled_accumulate(diagnosticsPool, forcingPool)
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
          call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'average', averagePool)

          call ocn_time_average_init(averagePool)
          call ocn_time_average_coupled_init(forcingPool)
          call ocn_time_average_coupled_accumulate(diagnosticsPool, forcingPool)
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

    call ocn_import_mct(x2o_o, errorCode)  
    if (errorCode /= 0) then
       call mpas_dmpar_global_abort('ERROR in ocn_import_mct')
    endif

    ! Build forcing arrays.
    block_ptr => domain % blocklist
    do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
        call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)

        call ocn_forcing_build_arrays(meshPool, statePool, forcingPool, ierr, 1)
        call ocn_forcing_build_fraction_absorbed_array(meshPool, statePool, forcingPool, ierr, 1)
        block_ptr => block_ptr % next
    end do

    currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
    call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
    write(stderrUnit,*) 'Initial time ', trim(timeStamp)

    itimestep = 0

    block_ptr => domain % blocklist
    do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'average', averagePool)
     
        call ocn_time_average_init(averagePool)

        block_ptr => block_ptr % next
    end do

    call ocn_analysis_compute_startup(domain_ptr, ierr) 

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

      type (mpas_pool_type), pointer :: meshPool, statePool, forcingPool, averagePool

      type (MPAS_Time_Type) :: currTime
      type (domain_type), pointer :: domain_ptr
      character(len=StrKIND) :: timeStamp
      type (MPAS_timeInterval_type) :: timeStep
      integer :: ierr
      logical, pointer :: config_write_output_on_startup
      character (len=StrKIND), pointer :: config_restart_timestamp_name

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

         ! Build forcing arrays.
         block_ptr => domain_ptr % blocklist
         do while(associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
            call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)

            call ocn_forcing_build_arrays(meshPool, statePool, forcingPool, ierr, 1)
            call ocn_forcing_build_fraction_absorbed_array(meshPool, statePool, forcingPool, ierr, 1)
            block_ptr => block_ptr % next
         end do

         call mpas_timer_start("time integration", .false.)
         call ocn_timestep(domain, dt, timeStamp)
         if (mpas_is_alarm_ringing(domain % clock, statsAlarmID, ierr=ierr)) then
            call mpas_reset_clock_alarm(domain % clock, statsAlarmID, ierr=ierr)

            call mpas_timer_start("global diagnostics")
            call ocn_compute_global_diagnostics(domain, 2, itimestep, dt);
            call mpas_timer_stop("global diagnostics")
         end if

         call mpas_timer_stop("time integration")

         ! Move time level 2 fields back into time level 1 for next time step
         block_ptr => domain_ptr % blocklist
         do while(associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
            call mpas_pool_shift_time_levels(statePool)
            block_ptr => block_ptr % next
         end do

         call ocn_analysis_compute_w_alarms(domain_ptr, ierr) 
         call ocn_analysis_write(domain % streamManager, ierr)

         ! Reset the restart alarm to prevent restart files being written without the coupler requesting it.
         call mpas_stream_mgr_reset_alarms(domain % streamManager, streamID='restart', ierr=ierr)

         call mpas_stream_mgr_write(domain % streamManager, streamID='output', ierr=ierr)
         call mpas_stream_mgr_reset_alarms(domain % streamManager, streamID='output', ierr=ierr)

         call mpas_stream_mgr_write(domain % streamManager, ierr=ierr)
         call mpas_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)

         call shr_sys_flush(stderrUnit)
      end do

      ! Check if coupler wants us to write a restart file.
      ! We only write restart files at the end of a coupling interval
      if (seq_timemgr_RestartAlarmIsOn(EClock)) then
         open(22, file=config_restart_timestamp_name, form='formatted', status='replace')
         write(22, *) trim(timeStamp)
         close(22)

         ! Write a restart file, because the coupler asked for it.
         call mpas_stream_mgr_write(domain % streamManager, forceWriteNow=.true., streamID='restart', ierr=ierr)
      end if

      ! Normalize time averaged fields
      block_ptr => domain_ptr % blocklist
      do while(associated(block_ptr))
         call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
         call ocn_time_average_coupled_normalize(forcingPool)
         block_ptr => block_ptr % next
      end do

      ! Export state to coupler
      call ocn_export_mct(o2x_o, ierr)

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

    integer :: ierr

!-----------------------------------------------------------------------

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

    call mpas_timer_stop("total time")
    call mpas_timer_write()
    call MPAS_stream_mgr_finalize(domain_ptr % streamManager)
    call mpas_framework_finalize(domain_ptr % dminfo, domain_ptr, io_system)

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

    type (mpas_pool_type), pointer :: meshPool

    integer, pointer :: nCellsSolve

    real (kind=RKIND), dimension(:), pointer :: lonCell, latCell, areaCell

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

    data(:) = 1.0_RKIND ! No land cells in MPAS-O, only Ocean cells
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
!    o  pslv   -- sea-level pressure                       (Pa)
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

   type (block_type), pointer :: block_ptr

   type (mpas_pool_type), pointer :: meshPool, forcingPool

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
                                  squaredWindSpeed10MeterField

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
                                               squaredWindSpeed10Meter
   
   logical, pointer :: bulkForcingActive

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

   call mpas_pool_get_package(domain % packages, 'bulkForcingActive', bulkForcingActive)

   n = 0
   block_ptr => domain % blocklist
   do while(associated(block_ptr))
      call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)

      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

      if ( bulkForcingActive ) then
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
         call mpas_pool_get_field(forcingPool, 'squaredWindSpeed10Meter', squaredWindSpeed10MeterField)
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
         squaredWindSpeed10Meter => squaredWindSpeed10MeterField % array
         iceRunoffFlux => iceRunoffFluxField % array

         do i = 1, nCellsSolve
           n = n + 1
           windStressZonal(i) = x2o_o % rAttr(index_x2o_Foxx_taux, n)
           windStressMeridional(i) = x2o_o % rAttr(index_x2o_Foxx_tauy, n)

           latentHeatFlux(i) = x2o_o % rAttr(index_x2o_Foxx_lat, n)
           sensibleHeatFlux(i) = x2o_o % rAttr(index_x2o_Foxx_sen, n)
           longWaveHeatFluxUp(i) = x2o_o % rAttr(index_x2o_Foxx_lwup, n)
           longWaveHeatFluxDown(i) = x2o_o % rAttr(index_x2o_Faxa_lwdn, n)
           evaporationFlux(i) = x2o_o % rAttr(index_x2o_Foxx_evap, n)
           seaIceHeatFlux(i) = x2o_o % rAttr(index_x2o_Fioi_melth, n)
           snowFlux(i) = x2o_o % rAttr(index_x2o_Faxa_snow, n)

           seaIceFreshWaterFlux(i) = x2o_o % rAttr(index_x2o_Fioi_meltw, n)
           seaIceSalinityFlux(i) = x2o_o % rAttr(index_x2o_Fioi_salt, n)
           riverRunoffFlux(i) = x2o_o % rAttr(index_x2o_Foxx_rofl, n)
           iceRunoffFlux(i) = x2o_o % rAttr(index_x2o_Foxx_rofi, n)

           shortWaveHeatFlux(i) = max(x2o_o % rAttr(index_x2o_Foxx_swnet, n), 0.0_RKIND)

           rainFlux(i) = x2o_o % rAttr(index_x2o_Faxa_rain, n)
           seaSurfacePressure(i) = x2o_o % rAttr(index_x2o_Sa_pslv, n)
           iceFraction(i) = x2o_o % rAttr(index_x2o_Si_ifrac, n)

! mrp causes out-of-bounds index error in debug mode
!          prognosticCO2(i) = x2o_o % rAttr(index_x2o_Sa_co2prog, n)
!          diagnosticCO2(i) = x2o_o % rAttr(index_x2o_Sa_co2diag, n)
           squaredWindSpeed10Meter(i) = x2o_o % rAttr(index_x2o_So_duu10n, n)

           if(iceRunoffFlux(n) < 0.0_RKIND) then
               call shr_sys_abort ('Error: incoming rofi_F is negative')
           end if
         end do
      end if

      block_ptr => block_ptr % next
   end do

   call mpas_pool_get_subpool(domain % blocklist % structs, 'forcing', forcingPool)

   if ( bulkForcingActive ) then
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
      call mpas_pool_get_field(forcingPool, 'squaredWindSpeed10Meter', squaredWindSpeed10MeterField)
      call mpas_pool_get_field(forcingPool, 'iceRunoffFlux', iceRunoffFluxField)

      call mpas_dmpar_exch_halo_field(windStressMeridionalField)
      call mpas_dmpar_exch_halo_field(windStressZonalField)
      call mpas_dmpar_exch_halo_field(latentHeatFluxField)
      call mpas_dmpar_exch_halo_field(sensibleHeatFluxField)
      call mpas_dmpar_exch_halo_field(longWaveHeatFluxUpField)
      call mpas_dmpar_exch_halo_field(longWaveHeatFluxDownField)
      call mpas_dmpar_exch_halo_field(evaporationFluxField)
      call mpas_dmpar_exch_halo_field(seaIceHeatFluxField)
      call mpas_dmpar_exch_halo_field(snowFluxField)
      call mpas_dmpar_exch_halo_field(seaIceFreshWaterFluxField)
      call mpas_dmpar_exch_halo_field(seaIceSalinityFluxField)
      call mpas_dmpar_exch_halo_field(riverRunoffFluxField)
      call mpas_dmpar_exch_halo_field(iceRunoffFluxField)
      call mpas_dmpar_exch_halo_field(shortWaveHeatFluxField)
      call mpas_dmpar_exch_halo_field(rainFluxField)
      call mpas_dmpar_exch_halo_field(seaSurfacePressureField)
      call mpas_dmpar_exch_halo_field(iceFractionField)
!     call mpas_dmpar_exch_halo_field(prognosticCO2Field)
!     call mpas_dmpar_exch_halo_field(diagnosticCO2Field)
!     call mpas_dmpar_exch_halo_field(squaredWindSpeed10MeterField)
   end if

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
   integer, pointer :: nCellsSolve, index_temperatureSurfaceValue, index_salinitySurfaceValue, index_avgZonalSurfaceVelocity, index_avgMeridionalSurfaceVelocity, index_avgZonalSSHGradient, index_avgMeridionalSSHGradient

   type (block_type), pointer :: block_ptr

   type (mpas_pool_type), pointer :: meshPool, forcingPool

   real (kind=RKIND), dimension(:), pointer :: seaIceEnergy, CO2Flux, DMSFlux, surfaceUpwardCO2Flux
   real (kind=RKIND), dimension(:,:), pointer :: avgTracersSurfaceValue, avgSurfaceVelocity, avgSSHGradient

   logical, pointer :: frazilIceActive

   errorcode = 0

   call mpas_pool_get_package(domain % packages, 'frazilIceActive', frazilIceActive)

   n = 0
   block_ptr => domain % blocklist
   do while(associated(block_ptr))
     call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
     call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)

     call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

     call mpas_pool_get_dimension(forcingPool, 'index_avgTemperatureSurfaceValue', index_temperatureSurfaceValue)
     call mpas_pool_get_dimension(forcingPool, 'index_avgSalinitySurfaceValue', index_salinitySurfaceValue)
     call mpas_pool_get_dimension(forcingPool, 'index_avgSurfaceVelocityZonal', index_avgZonalSurfaceVelocity)
     call mpas_pool_get_dimension(forcingPool, 'index_avgSurfaceVelocityMeridional', index_avgMeridionalSurfaceVelocity)
     call mpas_pool_get_dimension(forcingPool, 'index_avgSSHGradientZonal', index_avgZonalSSHGradient)
     call mpas_pool_get_dimension(forcingPool, 'index_avgSSHGradientMeridional', index_avgMeridionalSSHGradient)

     call mpas_pool_get_array(forcingPool, 'avgTracersSurfaceValue', avgTracersSurfaceValue)
     call mpas_pool_get_array(forcingPool, 'avgSurfaceVelocity', avgSurfaceVelocity)
     call mpas_pool_get_array(forcingPool, 'avgSSHGradient', avgSSHGradient)
     if ( frazilIceActive ) then
        call mpas_pool_get_array(forcingPool, 'seaIceEnergy', seaIceEnergy)
     end if
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
          o2x_o % rAttr(index_o2x_Fioo_q, n) = seaIceEnergy(i) / ocn_cpl_dt

          ! Reset SeaIce Energy
          seaIceEnergy(i) = 0.0_RKIND
       end if

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
