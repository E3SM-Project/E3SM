module glc_comp_mct

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
! !MODULE: glc_comp_mct
! !INTERFACE:

! !DESCRIPTION:
!  This is the main driver for the Model for Predication Across Scales Land Ice Model (MPASLI).
!
! !REVISION HISTORY:
!
! !AUTHOR:  Matt Hoffman
!
! !USES:
  ! cpl modules
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

  ! glc driver modules
  use glc_cpl_indices
  use glc_mct_vars

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
  use mpas_constants

  use iso_c_binding, only : c_char, c_loc, c_ptr, c_int
  use mpas_c_interfacing, only : mpas_f_to_c_string, mpas_c_to_f_string

  ! MPASLI modules
  use li_core
  use li_core_interface


  implicit none
  save
  private   ! By default make data private

  ! Public interfaces
  public :: glc_init_mct
  public :: glc_run_mct
  public :: glc_final_mct

  ! Private data interfaces
  private :: glc_SetGSMap_mct
  private :: glc_domain_mct

! !PRIVATE MODULE VARIABLES

  logical :: ldiag_cpl = .false.

  integer, private :: my_task

  integer, private ::   &
      cpl_write_restart,   &! flag id for write restart
      cpl_write_history,   &! flag id for write history
      cpl_write_tavg,      &! flag id for write tavg
      cpl_diag_global,     &! flag id for computing diagnostics
      cpl_diag_transp       ! flag id for computing diagnostics

  integer  :: nsend, nrecv

  character(len=StrKIND) :: runtype  !< ACME run type: initial, continue, branch
  character(len=StrKIND) :: coupleTimeStamp

  type(seq_infodata_type), pointer :: infodata  ! coupler infodata
  type (iosystem_desc_t), pointer :: io_system  ! ?

  ! MPAS Datatypes
  !type (dm_info), pointer :: dminfo
  type (core_type), pointer :: corelist => null()
  type (domain_type), pointer:: domain
  integer :: itimestep, &  ! time step number for MPAS
             glc_cpl_dt    ! length of coupling interval in seconds - set by coupler/ESMF
  character (len=*), parameter :: coupleAlarmID = 'coupling'  ! the MPAS alarm for coupling intervals
  integer :: glcStdErrUnit, glcStdOutUnit

!=======================================================================

contains

!***********************************************************************
! ! ROUTINE: glc_init_mct
!
! ! INTERFACE:
  subroutine glc_init_mct( EClock, cdata_g, x2g_g, g2x_g, NLFilename )!{{{
!
! ! DESCRIPTION:
! Initialize MPASLI
!
! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock), intent(inout) :: EClock
    type(seq_cdata), intent(inout) :: cdata_g
    type(mct_aVect), intent(inout) :: x2g_g, g2x_g
    character(len=*), optional, intent(in) :: NLFilename ! Namelist filename
!
! !REVISION HISTORY:
! Author: Matt Hoffman

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    integer :: glcID, mpicom_g, lsize, start_ymd, start_tod, start_year, start_day,   &
       start_month, start_hour, iyear, mpas_li_cpl_dt, shrloglev, shrlogunit, stdout_shr, pio_iotype

    type(mct_gsMap), pointer :: gsMap_g

    type(mct_gGrid), pointer :: dom_g

    integer :: errorCode  ! error code
    integer :: err_tmp

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
    type (MPAS_TimeInterval_type) :: denInterval, remInterval, zeroInterval
    integer :: numDivs
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
!   setup mpasli data structures
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

    call seq_cdata_setptrs(cdata_g, ID=glcID, mpicom=mpicom_g, &
         gsMap=gsMap_g, dom=dom_g, infodata=infodata)

    MPASLI_MCT_glcID =  glcID
    MPASLI_MCT_gsMap_g => gsMap_g
    MPASLI_MCT_dom_g   => dom_g

    call MPI_comm_rank(mpicom_g,iam,ierr)
#if (defined _MEMTRACE)
     if(iam == 0) then
         lbnum=1
         call memmon_dump_fort('memmon.out','glc_init_mct:start::',lbnum)
     endif
#endif

     ! Determine coupling type
     call seq_infodata_GetData(infodata, cpl_seq_option=cpl_seq_option)

!-----------------------------------------------------------------------
!
!   initialize the model run
!
!-----------------------------------------------------------------------

    call glc_cpl_indices_set()
    call seq_infodata_GetData( infodata, start_type=starttype)

    if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
       runtype = "initial"
    else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
       runtype = "continue"
    else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
       runtype = "branch"
    else
       write(stderrUnit, *) 'Start type is:', trim(starttype)
       call mpas_dmpar_global_abort(' glc_comp_mct ERROR: unknown starttype')
    end if

!-----------------------------------------------------------------------
!
!   first initializaiton phase of mpasli
!   initialize mpasli because grid information is needed for
!   creation of GSMap_glc.
!   call mpasli initialization routines
!
!-----------------------------------------------------------------------
    inst_suffix = seq_comm_suffix(glcID)

    call t_startf('mpasli_init')

    io_system => shr_pio_getiosys(glcid)

    pio_iotype = shr_pio_getiotype(glcid)
    call MPAS_io_set_iotype(domain % iocontext, pio_iotype)

    stdout_shr = shr_file_getUnit()
#ifndef MPAS_DEBUG
    if ( iam /= 0 ) then
       open(stdout_shr, status='SCRATCH')
    else
       inquire(file='glc_modelio.nml'//trim(inst_suffix),exist=exists)
       if(exists) then
          call shr_file_setio('glc_modelio.nml'//trim(inst_suffix),stdout_shr)
       endif
    end if
#else
    inquire(file='glc_modelio.nml'//trim(inst_suffix),exist=exists)
    if(exists) then
       call shr_file_setio('glc_modelio.nml'//trim(inst_suffix),stdout_shr)
    endif
#endif

    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (stdout_shr)
    glcStdErrUnit = stdout_shr
    glcStdOutUnit = stdout_shr

    ! Write to glcStdErrUnit here, because stderrUnit is not correctly set until the mpas_framework_init_phase1 call
    write(glcStdErrUnit,*) '=== Beginning glc_init_mct.'

    call mpas_framework_init_phase1(domain % dminfo, mpicom_g, stdoutUnit_in=stdout_shr, stderrUnit_in=stdout_shr)

    ! Setup function pointers for MPASLI core and domain types
    call li_setup_core(corelist)
    call li_setup_domain(domain)

    ! Override the names of the stream and namelist files
    domain % namelist_filename = 'mpasli_in'
    domain % streams_filename = 'streams.landice'

    ! Setup namelist variables, and read the namelist
    ierr = domain % core % setup_namelist(domain % configs, domain % namelist_filename, domain % dminfo)
    if ( ierr /= 0 ) then
       call mpas_dmpar_global_abort('ERROR: Namelist setup failed for core ' // trim(domain % core % coreName))
    end if

    call mpas_framework_init_phase2(domain, io_system)

    ! Define package variables
    ierr = domain % core % define_packages(domain % packages)
    if ( ierr /= 0 ) then
       call mpas_dmpar_global_abort('ERROR: Package definition failed for core ' // trim(domain % core % coreName))
    end if

    ! Setup packages (i.e. determine if they should be on or off)
    ierr = domain % core % setup_packages(domain % configs, domain % packages, domain % iocontext)
    if ( ierr /= 0 ) then
       call mpas_dmpar_global_abort('ERROR: Package setup failed for core ' // trim(domain % core % coreName))
    end if

    ! Setup decompositions available for dimensions
    ierr = domain % core % setup_decompositions(domain % decompositions)
    if ( ierr /= 0 ) then
       call mpas_dmpar_global_abort('ERROR: Decomposition setup failed for core ' // trim(domain % core % coreName))
    end if

    ! Override namelist options based on start type
    if (runtype == "initial") then ! Start up run

        ! Turn off restart
        call mpas_pool_get_config(domain % configs, "config_do_restart", tempLogicalConfig)
        tempLogicalConfig = .false.

        ! Setup start time. Will be over written later when clocks are synchronized
        call mpas_pool_get_config(domain % configs, "config_start_time", tempCharConfig)
        tempCharConfig = "0001-01-01_00:00:00"

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

    ! Setup MPASLI simulation clock
    ierr = domain % core % setup_clock(domain % clock, domain % configs)
    if ( ierr /= 0 ) then
       call mpas_dmpar_global_abort('ERROR: Clock setup failed for core ' // trim(domain % core % coreName))
    end if

    write(stderrUnit,*) 'Reading streams configuration from file '//trim(domain % streams_filename)
    inquire(file=trim(domain % streams_filename), exist=streamsExists)

    if ( .not. streamsExists ) then
       call mpas_dmpar_global_abort('ERROR: Streams file '//trim(domain % streams_filename)//' does not exist.')
    end if

    !
    ! Using information from the namelist, a graph.info file, and a file containing
    !    mesh fields, build halos and allocate blocks in the domain
    !
    ierr = domain % core % get_mesh_stream(domain % configs, mesh_stream)
    if ( ierr /= 0 ) then
       call mpas_dmpar_global_abort('ERROR: Failed to find mesh stream for core ' // trim(domain % core % coreName))
    end if


    call mpas_f_to_c_string(domain % streams_filename, c_filename)
    call mpas_f_to_c_string(mesh_stream, c_mesh_stream)
    c_comm = domain % dminfo % comm
    call xml_stream_get_attributes(c_filename, c_mesh_stream, c_comm, &
                                   c_mesh_filename_temp, c_ref_time_temp, &
                                   c_filename_interval_temp, c_iotype, c_ierr)
    if (c_ierr /= 0) then
       call mpas_dmpar_abort(domain % dminfo)
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

    write(stderrUnit, *) ' ** Attempting to bootstrap MPAS framework using stream: ', trim(mesh_stream)

    ! Bootstrap framework (1). Here data structures are setup, but dimensions and arrays are not finalized.
    call mpas_bootstrap_framework_phase1(domain, mesh_filename, mesh_iotype)

    !
    ! Set up run-time streams
    !
    call MPAS_stream_mgr_init(domain % streamManager, domain % ioContext, domain % clock, domain % blocklist % allFields, domain % packages, domain % blocklist % allStructs)

    call add_stream_attributes(domain % streamManager, domain)

    ! Setup all immutable streams for the core
    ierr = domain % core % setup_immutable_streams(domain % streamManager)
    if ( ierr /= 0 ) then
       call mpas_dmpar_global_abort('ERROR: Immutable streams setup failed for core ' // trim(domain % core % coreName))
    end if

    ! Parse / read all streams configuration
    mgr_p = c_loc(domain % streamManager)
    call xml_stream_parser(c_filename, mgr_p, c_comm, c_ierr)
    if (c_ierr /= 0) then
       call mpas_dmpar_abort(domain % dminfo)
    end if

    my_task = domain % dminfo % my_proc_id

    !
    ! Finalize the setup of blocks and fields
    !
    call mpas_bootstrap_framework_phase2(domain)


    ! Initialize the MPASLI core
    ierr = domain % core % core_init(domain, timeStamp)
    if ( ierr /= 0 ) then
       call mpas_dmpar_global_abort('ERROR: Core init failed for core ' // trim(domain % core % coreName))
    end if

!-----------------------------------------------------------------------
!
!   initialize time-stamp information
!
!-----------------------------------------------------------------------
    call t_stopf ('mpasli_init')

!-----------------------------------------------------------------------
!
!   check for consistency of mpasli and sync clock initial time
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
       call mpas_set_clock_time(domain % clock, currTime, MPAS_START_TIME, ierr)
       call mpas_set_clock_time(domain % clock, currTime, MPAS_NOW, ierr)
    end if

!-----------------------------------------------------------------------
!
!   initialize MCT attribute vectors and indices
!
!-----------------------------------------------------------------------

    call t_startf ('mpasli_mct_init')

    call glc_SetGSMap_mct( mpicom_g, glcID, GSMap_g )
    lsize = mct_gsMap_lsize(gsMap_g, mpicom_g)
    write(shrlogunit, *) ' MPASLI gsMap lsize is: ', lsize

    ! Initialize mct glc domain (needs glc initialization info)
    call glc_domain_mct( lsize, gsMap_g, dom_g )

    ! Initialize mct attribute vectors

    call mct_aVect_init(x2g_g, rList=seq_flds_x2g_fields, lsize=lsize)
    call mct_aVect_zero(x2g_g)

    call mct_aVect_init(g2x_g, rList=seq_flds_g2x_fields, lsize=lsize)
    call mct_aVect_zero(g2x_g)

    nsend = mct_avect_nRattr(g2x_g)
    nrecv = mct_avect_nRattr(x2g_g)

!-----------------------------------------------------------------------
!
!   initialize necessary coupling info
!
!-----------------------------------------------------------------------

    ! set glc_cpl_dt (coupling interval in integer seconds) from the EClock
    call seq_timemgr_EClockGetData(EClock, dtime=glc_cpl_dt)
    ! convert this dt to a time interval
    call mpas_set_timeInterval(alarmTimeStep, S=glc_cpl_dt, ierr=err_tmp)
    call mpas_get_timeInterval(alarmTimeStep, timeString=coupleTimeStamp, ierr=ierr)
    ierr = ior(ierr,err_tmp)
    write (stderrUnit,*) 'Applying GLC coupling dt (s) of: ', glc_cpl_dt

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
       write(stderrUnit, *) '        GLC Model time step is: ', trim(tempCharConfig)

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
       write(stderrUnit, *) '        GLC Model time step is: ', trim(tempCharConfig)

       write(stderrUnit, *) '        These are not synchronized, so time steps'
       write(stderrUnit, *) '        will not match to coupling interval boundaries.'
       write(stderrUnit, *) '        Please reconfigure either the coupling interval, '
       write(stderrUnit, *) '        or the time step.'
       call mpas_dmpar_global_abort('ERROR: Model is not properly configured for coupling interval.')
    end if

    ! set coupling alarm
    alarmStartTime = currTime
    call mpas_add_clock_alarm(domain % clock, coupleAlarmID, alarmStartTime, alarmTimeStep, ierr=err_tmp)
    ierr = ior(ierr,err_tmp)
    call mpas_reset_clock_alarm(domain % clock, coupleAlarmID, ierr=err_tmp)
    ierr = ior(ierr,err_tmp)


    ! Setup clock for initial runs
    if ( runtype == 'initial' ) then
       block_ptr => domain % blocklist
       do while(associated(block_ptr))
!          call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
!          call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)

!          call ocn_time_average_coupled_init(forcingPool)
!          call ocn_time_average_coupled_accumulate(diagnosticsPool, forcingPool)
          block_ptr => block_ptr % next
       end do
    end if

    !call mpas_print_alarm(domain % clock, couplealarmID, ierr) ! DEBUG to see alarm details

!-----------------------------------------------------------------------
!
!   get initial state from driver
!
!-----------------------------------------------------------------------

! GLC CURRENTLY DOES NOT REQUIRE IMPORT OF THE INITIAL STATE
!    call glc_import_mct(x2g_g, errorCode)
!    if (errorCode /= 0) then
!       call mpas_dmpar_global_abort('ERROR in glc_import_mct')
!    endif

! MPAS-O builds forcing arrays, which we don't have
!    ! Build forcing arrays.
!    block_ptr => domain % blocklist
!    do while(associated(block_ptr))
!        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
!        call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
!        call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
!        call ocn_forcing_build_arrays(meshPool, statePool, forcingPool, ierr, 1)
!        call ocn_forcing_build_transmission_array(meshPool, statePool, forcingPool, ierr, 1)
!        block_ptr => block_ptr % next
!    end do

!-----------------------------------------------------------------------
!
!   Calculate initial state in MPASLI (velocity, upperSurface, masks, etc.)
!
!-----------------------------------------------------------------------
    err_tmp = li_core_initial_solve(domain)
    ierr = ior(ierr,err_tmp)  ! li_core_finalize would abort if there was an error, but being safe.


    itimestep = 0  ! Initialize this module, save counter variable

! MPAS-O handles an average pool here, which we don't have
!    block_ptr => domain % blocklist
!    do while(associated(block_ptr))
!      call mpas_pool_get_subpool(block_ptr % structs, 'average', averagePool)
!      call ocn_time_average_init(averagePool)
!      block_ptr => block_ptr % next
!    end do

!    call ocn_analysis_compute_startup(domain, stream_manager, ierr)

    ! Reset all output alarms, to prevent initial time step from writing any output, unless it's ringing.
    call mpas_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)
    call mpas_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_INPUT, ierr=ierr)

!-----------------------------------------------------------------------
!
!   send initial state to driver
!
!-----------------------------------------------------------------------

    call glc_export_mct(g2x_g, errorCode)
    if (errorCode /= 0) then
       call mpas_dmpar_global_abort('ERROR in glc_export_mct')
    endif

    ! Setup clock for initial runs
! Some MPAS-O ocean pools handled here, but we don't have them.
!    if ( runtype == 'continue' ) then
!      block_ptr => domain % blocklist
!      do while(associated(block_ptr))
!        call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
!        call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
!        call mpas_pool_get_subpool(block_ptr % structs, 'average', averagePool)
!        call ocn_time_average_init(averagePool)
!        call ocn_time_average_coupled_init(forcingPool)
!        call ocn_time_average_coupled_accumulate(diagnosticsPool, forcingPool)
!        block_ptr => block_ptr % next
!      end do
!    end if

    call t_stopf ('mpasli_mct_init')

    call seq_infodata_PutData( infodata, glc_prognostic=.true.)


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
        call memmon_dump_fort('memmon.out','glc_init_mct:end::',lbnum)
        call memmon_reset_addr()
    endif
#endif
    call shr_sys_flush(stderrUnit)

    call check_clocks_sync(domain % clock, Eclock, stderrUnit, err_tmp)
    ierr = ior(ierr,err_tmp)

    write(stderrUnit,*) '=== glc_init_mct complete.'

#ifndef MPAS_DEBUG
    if ( iam /= 0 ) then
       close(stdout_shr)
    end if
#endif

!-----------------------------------------------------------------------
!EOC

  end subroutine glc_init_mct!}}}



!***********************************************************************
!
! !IROUTINE: glc_run_mct
!
! !INTERFACE:
  subroutine glc_run_mct( EClock, cdata_g, x2g_g, g2x_g)!{{{

    ! MPASLI modules
    use li_time_integration
    use li_statistics

    implicit none
!
! !DESCRIPTION:
! Run MPASLI for one coupling interval
!
! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata_g
    type(mct_aVect)             , intent(inout) :: x2g_g
    type(mct_aVect)             , intent(inout) :: g2x_g

!! !REVISION HISTORY:
!! Author: Matt Hoffman

!!-----------------------------------------------------------------------
!!
!!  local variables
!!
!!-----------------------------------------------------------------------
! Variables related to ACME coupler
      integer :: shrloglev, shrlogunit

! Variable related to MPASLI
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: geometryPool
      integer, pointer :: config_stats_interval   !< interval (number oftimesteps) for writing stats
      logical, pointer :: config_do_restart, config_write_output_on_startup, config_write_stats_on_startup
      character(len=StrKIND), pointer :: config_restart_timestamp_name

      type (MPAS_Time_Type) :: currTime
      character(len=StrKIND) :: timeStamp, streamName
      integer :: err, err_tmp, globalErr, streamDirection, my_proc_id
      logical :: solveVelo, streamActive

      my_proc_id = domain % dminfo % my_proc_id

#ifndef MPAS_DEBUG
      if ( my_proc_id /= 0 ) then
         open(glcStdErrUnit, status='SCRATCH')
      end if
#endif

      call mpas_io_units_set_stdout(glcStdOutUnit)
      call mpas_io_units_set_stderr(glcStdErrUnit)

      err = 0
      err_tmp = 0
      globalErr = 0

      ! Setup log information.
      call shr_file_getLogUnit (shrlogunit)
      call shr_file_getLogLevel(shrloglev)
      call shr_file_setLogUnit (stdoutUnit)

      write(stdoutUnit,*) '=== Beginning glc_run_mct.'

      ! Get config options that will be needed
      call mpas_pool_get_config(domain % configs, 'config_restart_timestamp_name', config_restart_timestamp_name)
      call mpas_pool_get_config(domain % configs, 'config_stats_interval', config_stats_interval)

      ! reset coupler alarm before we start

      call mpas_reset_clock_alarm(domain % clock, coupleAlarmID, ierr=err_tmp)
      err = ior(err, err_tmp)

      ! Get current time
      currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, err_tmp)
      err = ior(err, err_tmp)

      ! Import state from coupler
      call glc_import_mct(x2g_g, err_tmp)
      err = ior(err,err_tmp)

      ! Initialize time average fields
!      block_ptr => domain % blocklist
!      do while(associated(block_ptr))
!        call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
!        call ocn_time_average_coupled_init(forcingPool)
!        block_ptr => block_ptr % next
!      end do


      ! During integration, time level 1 stores the model state to be solved.
      ! For a variables with a second time level, it is the previous value.
      do while (.not. mpas_is_alarm_ringing(domain % clock, coupleAlarmID, ierr=err_tmp))
         err = ior(err,err_tmp)

         itimestep = itimestep + 1

         write(stderrUnit,*) 'Starting timestep number ', iTimeStep
         call shr_sys_flush(stderrUnit)

         ! ===
         ! === Perform Timestep
         ! ===
         call mpas_timer_start("time integration")

         call li_timestep(domain, err_tmp)
         err = ior(err,err_tmp)

         currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, err_tmp)
         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=err_tmp)

         ! Write statistics at designated interval
         if (config_stats_interval > 0) then
            if (mod(itimestep, config_stats_interval) == 0) then
               call mpas_timer_start("compute_statistics")
               call li_compute_statistics(domain, itimestep)
               call mpas_timer_stop("compute_statistics")
            end if
         end if

         call mpas_timer_stop("time integration")

         ! ===
         ! === Read time-varying inputs, if present (i.e., forcing)
         ! ===
         ! This should happen at the end of the time step so that if we write out
         ! the forcing it is at the correct time level.
         ! For an explicit time-stepping method, we want the forcing to be at the
         ! *old* time when it is applied during time integration.  Reading it here
         ! will allow that.
         ! Finally, set whence to latest_before so we have piecewise-constant forcing.
         ! Could add, e.g., linear interpolation later.
! In coupled mode, we don't want to use forcing from a file - that should come from the coupler!
!         call mpas_stream_mgr_read(stream_manager, whence=MPAS_STREAM_LATEST_BEFORE, ierr=err_tmp)
!         err = ior(err, err_tmp)
!         call mpas_stream_mgr_reset_alarms(stream_manager, direction=MPAS_STREAM_INPUT, ierr=err_tmp)
!         err = ior(err, err_tmp)


         ! ===
         ! === Write Output and/or Restart, if needed
         ! ===
         call mpas_timer_start("write output")

         ! Reset any restart alarms to prevent restart files being written without the coupler requesting it.
         call mpas_stream_mgr_begin_iteration(domain % streamManager)

         do while ( mpas_stream_mgr_get_next_stream( domain % streamManager, streamID=streamName, &
                    directionProperty=streamDirection, activeProperty=streamActive ) )

            if ( streamActive .and. streamDirection == MPAS_STREAM_INPUT_OUTPUT ) then
               call mpas_stream_mgr_reset_alarms(domain % streamManager, streamID=streamName, ierr=err_tmp)
               err = ior(err,err_tmp)
            end if
         end do
         err = ior(err, err_tmp)

         ! These calls will handle ALL output streams that need writing.
         ! [Could add them individually, as the ocean does, if some other actions are needed when a
         ! specific alarm is ringing (e.g., global stats calculated only when output stream gets written)]
         call mpas_stream_mgr_write(domain % streamManager, ierr=err_tmp)
         err = ior(err, err_tmp)
         call mpas_stream_mgr_reset_alarms(domain % streamManager,   direction=MPAS_STREAM_OUTPUT, ierr=err_tmp)
         err = ior(err, err_tmp)
         call mpas_timer_stop("write output")

         ! Move time level 1 fields (current values) into time level 2 (old values) for next time step
         ! (for those fields with multiple time levels)
         block => domain % blocklist
         do while(associated(block))
            call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
            call mpas_pool_shift_time_levels(geometryPool)
            block => block % next
         end do

         ! Reset the alarm for checking for force setting of the adaptive timestep interval
         if (mpas_is_alarm_ringing(domain % clock, 'adaptiveTimestepForceInterval', ierr=err_tmp)) then
            err = ior(err, err_tmp)
            call mpas_reset_clock_alarm(domain % clock, 'adaptiveTimestepForceInterval', ierr=err_tmp)
            err = ior(err, err_tmp)
         endif
         err = ior(err, err_tmp)

         ! === error check and exit
         call mpas_dmpar_max_int(domain % dminfo, err, globalErr)  ! Find out if any blocks got an error
         if (globalErr > 0) then
             call mpas_dmpar_global_abort("An error has occurred in mpas_core_run. Aborting...")
         endif

         write(stderrUnit,*) 'Completed timestep number ', iTimeStep
         call shr_sys_flush(stderrUnit)
      end do
      err = ior(err, err_tmp)

      ! Move time levels back so current values are in time level 1 for any additional output/export
      block => domain % blocklist
      do while(associated(block))
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_shift_time_levels(geometryPool)
         block => block % next
      end do


      ! Check if coupler wants us to write a restart file.
      ! We only write restart files at the end of a coupling interval
      if (seq_timemgr_RestartAlarmIsOn(EClock)) then

         ! Write a restart file, because the coupler asked for it.
         call mpas_stream_mgr_begin_iteration(domain % streamManager)

         do while ( mpas_stream_mgr_get_next_stream( domain % streamManager, streamID=streamName, &
                    directionProperty=streamDirection, activeProperty=streamActive ) )

            if ( streamActive .and. streamDirection == MPAS_STREAM_INPUT_OUTPUT ) then
               call mpas_stream_mgr_write(domain % streamManager, forceWriteNow=.true., streamID=streamName, ierr=err_tmp)
               err = ior(err,err_tmp)
            end if
         end do

         if ( domain % dminfo % my_proc_id == 0 ) then
            open(22, file=config_restart_timestamp_name, form='formatted', status='replace')
            write(22, *) trim(timeStamp)
            close(22)
         end if
      end if

      ! Normalize time averaged fields
!      block_ptr => domain % blocklist
!      do while(associated(block_ptr))
!        call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
!        call ocn_time_average_coupled_normalize(forcingPool)
!        block_ptr => block_ptr % next
!      end do

      ! Export state to coupler
      call glc_export_mct(g2x_g, err_tmp)
      err = ior(err,err_tmp)

      ! Move time levels back the way MPAS will expect them on the next time step
      block => domain % blocklist
      do while(associated(block))
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_shift_time_levels(geometryPool)
         block => block % next
      end do

      call check_clocks_sync(domain % clock, Eclock, stderrUnit, err_tmp)
      err = ior(err,err_tmp)

      write(stdoutUnit,*) '=== Completed glc_run_mct.'

      ! Reset I/O logs
      call shr_file_setLogUnit (shrlogunit)
      call shr_file_setLogLevel(shrloglev)

#ifndef MPAS_DEBUG
      if ( my_proc_id /= 0 ) then
         close(glcStdErrUnit)
      end if
#endif

!-----------------------------------------------------------------------
  end subroutine glc_run_mct!}}}



!***********************************************************************
!
! !ROUTINE: glc_final_mct
!
! !INTERFACE:
  subroutine glc_final_mct( EClock, cdata_g, x2g_g, g2x_g)!{{{

    ! MPAS Framework modules
    use mpas_stream_manager, only : MPAS_stream_mgr_finalize
    ! MPAS Land Ice modules
    use li_velocity, only: li_velocity_finalize
!
! !DESCRIPTION:
! Finalize MPASLI
!
! !USES:
! !ARGUMENTS:
    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata_g
    type(mct_aVect)             , intent(inout) :: x2g_g
    type(mct_aVect)             , intent(inout) :: g2x_g
!
! !REVISION HISTORY:
! Author: Matt Hoffman

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
    integer :: shrloglev, shrlogunit, my_proc_id
    integer :: err, err_tmp, iErr

!-----------------------------------------------------------------------

    my_proc_id = domain % dminfo % my_proc_id

#ifndef MPAS_DEBUG
    if ( my_proc_id /= 0 ) then
       open(glcStdErrUnit, status='SCRATCH')
    end if
#endif

    call mpas_io_units_set_stdout(glcStdOutUnit)
    call mpas_io_units_set_stderr(glcStdErrUnit)

    ! Setup I/O logs
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (stdoutUnit)

    ! Finalize MPASLI
    iErr = domain % core % core_finalize(domain)
    if ( iErr /= 0 ) then
       call mpas_dmpar_global_abort('ERROR: Core finalize failed for core ' // trim(domain % core % coreName))
    end if

    call mpas_timer_write()
    call MPAS_stream_mgr_finalize(domain % streamManager)
    call mpas_framework_finalize(domain % dminfo, domain, io_system)

    ! Reset I/O logs
    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

#ifndef MPAS_DEBUG
    if ( my_proc_id /= 0 ) then
       close(glcStdErrUnit)
    end if
#endif

  end subroutine glc_final_mct!}}}



!***********************************************************************
! !ROUTINE: glc_SetGSMap_mct
! !INTERFACE:

 subroutine glc_SetGSMap_mct( mpicom_glc, GLCID, gsMap_glc )!{{{

   use mpas_dmpar

! !DESCRIPTION:
!  This routine sets up the MPASLI grid numbering for MCT
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

    implicit none
    integer        , intent(in)    :: mpicom_glc
    integer        , intent(in)    :: GLCID
    type(mct_gsMap), intent(inout) :: gsMap_glc

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    integer,allocatable :: &
      gindex(:)

    integer ::   &
      i,j, k, n, iblock, &
      lsize, gsize,   &
      ier

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
    call mct_gsMap_init( gsMap_glc, gindex, mpicom_glc, GLCID, lsize, gsize )

    deallocate(gindex)

!-----------------------------------------------------------------------
  end subroutine glc_SetGSMap_mct!}}}



!***********************************************************************
! !ROUTINE: glc_domain_mct
! !INTERFACE:

  subroutine glc_domain_mct( lsize, gsMap_g, dom_g )!{{{

! !DESCRIPTION:
!  This routine sets up the MCT domain for MPASLI
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT/OUTPUT PARAMETERS:

    use shr_const_mod, only: SHR_CONST_PI

    implicit none
    integer        , intent(in)    :: lsize
    type(mct_gsMap), intent(in)    :: gsMap_g
    type(mct_ggrid), intent(inout) :: dom_g

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    integer, pointer :: idata(:)

    real(kind=RKIND), pointer :: data(:)
    real(kind=RKIND) :: r2d

    integer :: i,j, k, n, ier

    type (block_type), pointer :: block_ptr

    type (mpas_pool_type), pointer :: meshPool

    integer, pointer :: nCellsSolve

    real (kind=RKIND), dimension(:), pointer :: lonCell, latCell, areaCell

    real (kind=RKIND), pointer :: sphere_radius

    r2d = 180.0_RKIND / SHR_CONST_PI

!-------------------------------------------------------------------
!
!  initialize mct domain type
!
!-------------------------------------------------------------------

    call mct_gGrid_init( GGrid=dom_g, CoordChars=trim(seq_flds_dom_coord), &
       OtherChars=trim(seq_flds_dom_other), lsize=lsize )
    call mct_aVect_zero(dom_g%data)
    allocate(data(lsize))

!-------------------------------------------------------------------
!
! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
!
!-------------------------------------------------------------------

    call mct_gsMap_orderedPoints(gsMap_g, my_task, idata)
    call mct_gGrid_importIAttr(dom_g,'GlobGridNum',idata,lsize)

!-------------------------------------------------------------------
!
! Determine domain
! Initialize attribute vector with special value
!
!-------------------------------------------------------------------

    data(:) = -9999.0_RKIND
    call mct_gGrid_importRAttr(dom_g,"lat"  ,data,lsize)
    call mct_gGrid_importRAttr(dom_g,"lon"  ,data,lsize)
    call mct_gGrid_importRAttr(dom_g,"area" ,data,lsize)
    call mct_gGrid_importRAttr(dom_g,"aream",data,lsize)
    data(:) = 1.0_RKIND
    call mct_gGrid_importRAttr(dom_g,"mask",data,lsize)
    call mct_gGrid_importRAttr(dom_g,"frac",data,lsize)

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
    !write(stderrUnit,*) 'lonrange',minval(data), maxval(data), n
    call mct_gGrid_importRattr(dom_g,"lon",data,lsize)

    n = 0
    block_ptr => domain % blocklist
    do while(associated(block_ptr))
      call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

      call mpas_pool_get_array(meshPool, 'latCell', latCell)

      !write(stderrUnit,*) 'lsize=',lsize,' nCellsSolve=', nCellsSolve

      do i = 1, nCellsSolve
        n = n + 1
        data(n) = latCell(i) * r2d
      end do
      block_ptr => block_ptr % next
    end do
    !write(stderrUnit,*) 'latrange',minval(data), maxval(data), n
    call mct_gGrid_importRattr(dom_g,"lat",data,lsize)

    n = 0
    block_ptr => domain % blocklist
    do while(associated(block_ptr))
      call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

      ! Note: for spherical meshes, we should get the sphere_radius!
      ! But currently only supporting planar meshes
      !call mpas_pool_get_config(meshPool, 'sphere_radius', sphere_radius)

      do i = 1, nCellsSolve
        n = n + 1
        data(n) = areaCell(i) * 0.24635127E-13_RKIND  ! This is the conversion factor to rad^2
      end do
      block_ptr => block_ptr % next
    end do
    !write(stderrUnit,*) 'arearange',minval(data), maxval(data), n
    call mct_gGrid_importRattr(dom_g,"area",data,lsize)

    ! For now, assume mask and frac are 1 everywhere.
    ! This may need to be changed in the future.
    data(:) = 1.0_RKIND
    call mct_gGrid_importRattr(dom_g,"mask",data,lsize)
    call mct_gGrid_importRattr(dom_g,"frac",data,lsize)

    deallocate(data)
    deallocate(idata)

!-----------------------------------------------------------------------
  end subroutine glc_domain_mct!}}}


!=================================================================================


   subroutine glc_import_mct(x2g_g, errorCode)

! !INPUT/OUTPUT PARAMETERS:

    type(mct_aVect)   , intent(inout) :: x2g_g

! !OUTPUT PARAMETERS:

   integer, intent(out) :: &
      errorCode              ! returned error code

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

   type (block_type), pointer :: block

   type (mpas_pool_type), pointer :: meshPool
   type (mpas_pool_type), pointer :: geometryPool
   type (mpas_pool_type), pointer :: thermalPool

   integer, pointer :: nCellsSolve

   real (kind=RKIND), dimension(:), pointer :: sfcMassBal,&
                                               basalOceanMassBal,&
                                               surfaceTemperature,&
                                               basalOceanHeatflx,&
                                               OceanDensity

    errorCode = 0

    n = 0
    block => domain % blocklist

    do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

         ! Get variables from pools
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         !call mpas_pool_get_subpool(block % structs, 'thermal', thermalPool)

         call mpas_pool_get_array(geometryPool, 'sfcMassBal', sfcMassBal)
         !call mpas_pool_get_array(geometryPool, 'basalOceanMassBal',basalOceanMassBal)
         !call mpas_pool_get_array(thermalPool, 'surfaceTemperature',surfaceTemperature)
         !call mpas_pool_get_array(thermalPool, 'basalOceanHeatflx',basalOceanHeatflx)
         !call mpas_pool_get_array(geometryPool, 'OceanDensity',OceanDensity)

         do i = 1, nCellsSolve
            n = n + 1
            sfcMassBal(i) = x2g_g % rAttr(index_x2g_Flgl_qice, n)
            !basalOceanMassBal(i) = x2g_g % rAttr(index_x2g_Fogo_qicel, n)
            !surfaceTemperature(i) = x2g_g % rAttr(index_x2g_Sl_tsrf, n)
            !basalOceanHeatflx(i) = x2g_g % rAttr(index_x2g_Fogo_qiceh, n)
            !OceanDensity (i) = x2g_g % rAttr(, n)
         end do

         block => block % next
    end do


!Need to loop over blocks
!Jer: Need to define local pointers to MPAS-LI arrays, get these arrays, and point to them with local pointers
!Jer: Need to assign local pointers to x2g arrays (as in CESM)
!Jer: then, need to re-get forcingPool fields and call maps_dmpar_exch_halo_field on them

  end subroutine glc_import_mct


!=================================================================================


  subroutine glc_export_mct(g2x_g, errorCode)

   use li_mask

   !-------------------------------------------------------------------

   type(mct_aVect) , intent(inout) :: g2x_g
   integer, intent(out) :: errorCode

   !-------------------------------------------------------------------
   integer ::  &
      i,n

   type (block_type), pointer :: block

   real (kind=RKIND), pointer :: config_sea_level

   type (mpas_pool_type), pointer :: meshPool
   type (mpas_pool_type), pointer :: geometryPool

   integer, pointer :: nCellsSolve

   real (kind=RKIND), dimension(:), pointer :: upperSurface
   integer, dimension(:), pointer :: cellMask
   !-------------------------------------------------------------------

    errorCode = 0

    n = 0
    block => domain % blocklist

    call mpas_pool_get_config(domain % configs, 'config_sea_level', config_sea_level)

    do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

         ! Get variables from pools
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)

         call mpas_pool_get_array(geometryPool, 'upperSurface', upperSurface)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)

         do i = 1, nCellsSolve
            n = n + 1
            g2x_g % rAttr(index_g2x_Sg_topo, n) = upperSurface(i)
            g2x_g % rAttr(index_g2x_Sg_ice_covered, n) = li_mask_is_ice_int(cellMask(i))
            ! Sg_icemask should be 1 in locations where GLC will use a SMB value
            ! This is cells with ice or bare land (no ocean cells)
            ! This can easily be checked by any cells where upperSurface is
            ! above sea level!
            if (upperSurface(i) > config_sea_level) then
               g2x_g % rAttr(index_g2x_Sg_icemask, n) = 1
            else
               g2x_g % rAttr(index_g2x_Sg_icemask, n) = 0
            endif

            ! icemask_coupled_fluxes for now should be the same as icemask
            ! If we add in the ability to prevent the ice sheet from evolving
            ! (i.e., to ignore these fluxes) then this should be set to 0 everywhere.
            g2x_g % rAttr(index_g2x_Sg_icemask_coupled_fluxes, n) = g2x_g % rAttr(index_g2x_Sg_icemask, n)

         end do

         block => block % next
    end do

  end subroutine glc_export_mct


!=================================================================================


  subroutine check_clocks_sync(MPASclock, Eclock, stderrUnit, err)
    type (MPAS_Clock_type), intent(in) :: MPASclock
    type(ESMF_Clock), intent(in)       :: EClock
    integer, intent(in)                :: stderrUnit
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
       write(stderrUnit,*) 'MPAS ymd=',ymd,' MPAS tod=', tod
       write(stderrUnit,*) 'sync ymd=',curr_ymd,' sync tod=', curr_tod
       write(stderrUnit,*) 'Internal mpas clock not in sync with sync clock'
       call mpas_dmpar_global_abort("Internal MPAS clock not in sync with ACME sync clock. Aborting...")
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

   subroutine add_stream_attributes(stream_manager, domain)!{{{

      use mpas_stream_manager, only : MPAS_stream_mgr_add_att

      implicit none

      type (MPAS_streamManager_type), intent(inout) :: stream_manager
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


      call MPAS_stream_mgr_add_att(stream_manager, 'on_a_sphere', domain % on_a_sphere)
      call MPAS_stream_mgr_add_att(stream_manager, 'sphere_radius', domain % sphere_radius)
      call MPAS_stream_mgr_add_att(stream_manager, 'model_name', domain % core % modelName)
      call MPAS_stream_mgr_add_att(stream_manager, 'core_name', domain % core % coreName)
      ! DWJ 10/01/2014: Eventually add the real history attribute, for now (due to length restrictions)
      ! add a shortened version.
!     call MPAS_stream_mgr_add_att(stream_manager, 'history', domain % history)
      call MPAS_stream_mgr_add_att(stream_manager, 'history', histAtt)
      call MPAS_stream_mgr_add_att(stream_manager, 'source', domain % core % source)
      call MPAS_stream_mgr_add_att(stream_manager, 'Conventions', domain % core % conventions)
      call MPAS_stream_mgr_add_att(stream_manager, 'parent_id', domain % parent_id)
      call MPAS_stream_mgr_add_att(stream_manager, 'mesh_spec', domain % mesh_spec)
      call MPAS_stream_mgr_add_att(stream_manager, 'git_version', domain % core % git_version)

      call mpas_pool_begin_iteration(domain % configs)

      do while (mpas_pool_get_next_member(domain % configs, itr))

         if ( itr % memberType == MPAS_POOL_CONFIG) then

            if ( itr % dataType == MPAS_POOL_REAL ) then
               call mpas_pool_get_config(domain % configs, itr % memberName, realAtt)
               call MPAS_stream_mgr_add_att(stream_manager, itr % memberName, realAtt, ierr=local_ierr)
            else if ( itr % dataType == MPAS_POOL_INTEGER ) then
               call mpas_pool_get_config(domain % configs, itr % memberName, intAtt)
               call MPAS_stream_mgr_add_att(stream_manager, itr % memberName, intAtt, ierr=local_ierr)
            else if ( itr % dataType == MPAS_POOL_CHARACTER ) then
               call mpas_pool_get_config(domain % configs, itr % memberName, charAtt)
               call MPAS_stream_mgr_add_att(stream_manager, itr % memberName, charAtt, ierr=local_ierr)
            else if ( itr % dataType == MPAS_POOL_LOGICAL ) then
               call mpas_pool_get_config(domain % configs, itr % memberName, logAtt)
               if (logAtt) then
                  call MPAS_stream_mgr_add_att(stream_manager, itr % memberName, 'YES', ierr=local_ierr)
               else
                  call MPAS_stream_mgr_add_att(stream_manager, itr % memberName, 'NO', ierr=local_ierr)
               end if
            end if

          end if
      end do

   end subroutine add_stream_attributes!}}}


end module glc_comp_mct

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
