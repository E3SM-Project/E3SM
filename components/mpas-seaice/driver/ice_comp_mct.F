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
   use seq_comm_mct,      only : seq_comm_suffix, seq_comm_inst, seq_comm_name, info_taskmap_comp
   use shr_file_mod
   use shr_cal_mod,       only : shr_cal_date2ymd
   use shr_sys_mod
   use shr_taskmap_mod,   only : shr_taskmap_write
   use shr_pio_mod
   use perf_mod
   use shr_strdata_mod,   only : shr_strdata_type

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
#ifdef HAVE_MOAB
   use mpas_moabmesh
   use seq_comm_mct, only: MPSIID
   use iMOAB, only: iMOAB_DefineTagStorage, iMOAB_SetDoubleTagStorage
   use shr_kind_mod    , only: cxx => SHR_KIND_CXX
   use seq_comm_mct, only:  seq_comm_compare_mb_mct
   use seq_comm_mct,      only: num_moab_exports
#endif

   use iso_c_binding, only : c_char, c_loc, c_ptr, c_int
   use mpas_c_interfacing, only : mpas_f_to_c_string, mpas_c_to_f_string

   ! MPASSI modules
   use seaice_analysis_driver
   use seaice_column, only: seaice_column_reinitialize_fluxes, & !colpkg
                            seaice_column_coupling_prep
   use seaice_icepack, only: seaice_icepack_reinitialize_fluxes, &
                             seaice_icepack_coupling_prep
   use seaice_constants, only: coupleAlarmID, &
                               seaiceFreshWaterFreezingPoint, &
                               seaiceDensityIce, &
                               seaiceDensitySnow, &
                               seaiceDensitySeaWater, &
                               seaiceGravity, &
                               seaiceLatentHeatMelting, &
                               seaiceReferenceSalinity, &
                               seaiceFrazilSalinityReduction, &
                               seaiceFrazilIcePorosity, &
                               pii

   use seaice_core
   use seaice_core_interface
   use seaice_forcing
   use seaice_initialize, only: seaice_init_post_clock_advance
   use seaice_mesh, only: seaice_latlon_vector_rotation_backward
   use seaice_time_integration
   use seaice_error, only: seaice_check_critical_error

   use ice_colpkg, only: colpkg_sea_freezing_temperature
   use icepack_intfc, only: icepack_sea_freezing_temperature

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
#ifdef HAVE_MOAB
  private :: ice_export_moab
  private :: ice_import_moab
  integer , private :: mblsize, totalmbls,totalmblr
  real (kind=RKIND) , allocatable, private :: i2x_im(:,:)
  real (kind=RKIND) , allocatable, private :: x2i_im(:,:)
  integer :: mpicom_moab  ! save it for runtime debug
#endif
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

   ! prescribed ice mode
   type(shr_strdata_type) :: &
        sdat ! prescribed data stream

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

    type (MPAS_Pool_Type), pointer :: &
         shortwave, &
         berg_forcing

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
    character(len=strKIND) :: caseid
    character(len=strKIND) :: model_version
    character(len=strKIND) :: username
    character(len=strKIND) :: hostname
    character(len=strKIND) :: institution
    character(len=strKIND) :: curdate
    character(len=strKIND) :: curtime
    character(len=strKIND) :: history
#ifdef HAVE_MOAB
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
    logical, pointer :: &
         tempLogicalConfig, &
         config_scale_dib_by_removed_ice_runoff
    character(len=StrKIND), pointer :: tempCharConfig
    real (kind=RKIND), pointer :: tempRealConfig

    real(kind=RKIND), pointer :: &
         dayOfNextShortwaveCalculation, & ! needed for CESM like coupled simulations
         runningMeanRemovedIceRunoff ! the area integrated, running mean of removed ice runoff from the ocean

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
    inst_index = seq_comm_inst(iceID)
    write(c_inst_index,'(i8)') inst_index

    if (info_taskmap_comp > 0) then

       if (info_taskmap_comp == 1) then
          verbose_taskmap_output = .false.
       else
          verbose_taskmap_output = .true.
       endif

       call MPI_comm_size(mpicom_i,npes,ierr)
       write(c_npes,'(i8)') npes

       if (iam==0) then
          write(iceLogUnit,'(/,3A)') &
             trim(adjustl(c_npes)), &
             ' pes participating in computation of MPAS-SEAICE instance #', &
             trim(adjustl(c_inst_index))
          call shr_sys_flush(iceLogUnit)
       endif

       call t_startf("shr_taskmap_write")
       call shr_taskmap_write(iceLogUnit, mpicom_i,                 &
                              'ICE #'//trim(adjustl(c_inst_index)), &
                              verbose=verbose_taskmap_output        )
       call t_stopf("shr_taskmap_write")

    endif

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
    ! Disable output from all but the master task for E3SM!
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
    domain % namelist_filename = 'mpassi_in' // trim(inst_suffix)
    domain % streams_filename = 'streams.seaice' // trim(inst_suffix)

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

        ! Turn on restarts
        call mpas_pool_get_config(domain % configs, "config_do_restart", tempLogicalConfig)
        tempLogicalConfig = .true.
        call mpas_pool_get_config(domain % configs, "config_do_restart_bgc", tempLogicalConfig)
        tempLogicalConfig = .true.
        call mpas_pool_get_config(domain % configs, "config_do_restart_hbrine", tempLogicalConfig)
        tempLogicalConfig = .true.
        call mpas_pool_get_config(domain % configs, "config_do_restart_snow_density", tempLogicalConfig)
        tempLogicalConfig = .true.
        call mpas_pool_get_config(domain % configs, "config_do_restart_snow_grain_radius", tempLogicalConfig)
        tempLogicalConfig = .true.
        call mpas_pool_get_config(domain % configs, "config_do_restart_zsalinity", tempLogicalConfig)
        tempLogicalConfig = .true.
        call mpas_pool_get_config(domain % configs, "config_do_restart_bgc", tempLogicalConfig)
        tempLogicalConfig = .true.


        ! Set start time to be read from file
        call mpas_pool_get_config(domain % configs, "config_start_time", tempCharConfig)
        tempCharConfig = "file"

        ! Setup run duration. Will be ignored in coupled run, since coupler defines how long the run is.
        call mpas_pool_get_config(domain % configs, "config_run_duration", tempCharConfig)
        tempCharConfig = "0001-00-00_00:00:00"

    end if

    call MPAS_pool_get_config(domain % configs, "config_scale_dib_by_removed_ice_runoff", &
                              config_scale_dib_by_removed_ice_runoff)

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

    ! add E3SM stream attributes to MPAS output files
    call seq_infodata_GetData(infodata, case_name=caseid, model_version=model_version, username=username, hostname=hostname)
    institution = &
         'LLNL (Lawrence Livermore National Laboratory, Livermore, CA 94550, USA); &
         &ANL (Argonne National Laboratory, Argonne, IL 60439, USA); &
         &BNL (Brookhaven National Laboratory, Upton, NY 11973, USA); &
         &LANL (Los Alamos National Laboratory, Los Alamos, NM 87545, USA); &
         &LBNL (Lawrence Berkeley National Laboratory, Berkeley, CA 94720, USA); &
         &ORNL (Oak Ridge National Laboratory, Oak Ridge, TN 37831, USA); &
         &PNNL (Pacific Northwest National Laboratory, Richland, WA 99352, USA); &
         &SNL (Sandia National Laboratories, Albuquerque, NM 87185, USA). &
         &Mailing address: LLNL Climate Program, c/o David C. Bader, &
         &Principal Investigator, L-103, 7000 East Avenue, Livermore, CA 94550, USA'
    call datetime(curdate, curtime)
    history = 'created on ' // trim(curdate) // ' ' // trim(curtime)

    call MPAS_stream_mgr_add_att(domain % streamManager, 'case', trim(caseid))
    call MPAS_stream_mgr_add_att(domain % streamManager, 'source_id', trim(model_version))
    call MPAS_stream_mgr_add_att(domain % streamManager, 'realm', 'seaIce')
    call MPAS_stream_mgr_add_att(domain % streamManager, 'product', 'model-output')
    call MPAS_stream_mgr_add_att(domain % streamManager, 'title', 'MPAS-Seaice History file information')
    call MPAS_stream_mgr_add_att(domain % streamManager, 'source', 'E3SM Sea Ice Model')
    call MPAS_stream_mgr_add_att(domain % streamManager, 'institution',trim(institution))
    call MPAS_stream_mgr_add_att(domain % streamManager, 'institution_id', 'E3SM-Project')
    call MPAS_stream_mgr_add_att(domain % streamManager, 'contact', 'e3sm-data-support@llnl.gov')
    call MPAS_stream_mgr_add_att(domain % streamManager, 'username', trim(username))
    call MPAS_stream_mgr_add_att(domain % streamManager, 'hostname', trim(hostname))

    ! Other attributes
    call add_stream_attributes(domain)

    ! overwrite MPAS attributes
    call MPAS_stream_mgr_add_att(domain % streamManager, 'history', trim(history))
    call MPAS_stream_mgr_add_att(domain % streamManager, 'Conventions', 'CF-1.7')
    call MPAS_stream_mgr_add_att(domain % streamManager, 'git_version', trim(model_version))

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

    if (config_scale_dib_by_removed_ice_runoff) then
       ! independent of space so should be no need to loop over blocks
       block => domain % blocklist
       call MPAS_pool_get_subpool(block % structs, "berg_forcing", berg_forcing)
       call MPAS_pool_get_array(berg_forcing, "runningMeanRemovedIceRunoff", &
                                runningMeanRemovedIceRunoff)
       call seq_infodata_GetData(infodata, rmean_rmv_ice_runoff=runningMeanRemovedIceRunoff )
    end if

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
#ifdef HAVE_MOAB
    call init_moab_mpas(domain, ICEID, MPSIID) ! should return MPSIID ..
    call mpas_log_write('initialized MOAB MPAS sea-ice instance... ')
#endif
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
#ifdef HAVE_MOAB
   ! initialize moab tag fields array
    mblsize = lsize
    totalmbls = mblsize * nsend ! size of the double array
    totalmblr = mblsize * nrecv ! size of the double array
    allocate (i2x_im(lsize, nsend) )
    allocate (x2i_im(lsize, nrecv) )
    i2x_im = 0._r8
    x2i_im = 0._r8
    ! define tags according to the seq_flds_i2x_fields
    ! also zero them out
    tagtype = 1  ! dense, double
    numco = 1 !  one value per cell / entity
    ent_type = 1 ! cells
    tagname = trim(seq_flds_i2x_fields)//C_NULL_CHAR
    ierrmb = iMOAB_DefineTagStorage(MPSIID, tagname, tagtype, numco,  tagindex )
    if ( ierrmb == 1 ) then
      call mpas_log_write('cannot define tags for MOAB i2x fields ' // trim(seq_flds_i2x_fields), MPAS_LOG_ERR)
    endif
    ierrmb = iMOAB_SetDoubleTagStorage ( MPSIID, tagname, totalmbls , ent_type, i2x_im(1, 1) )
    if ( ierrmb /= 0 ) then
       call mpas_log_write('cannot set tags for MOAB i2x fields to zero' // trim(seq_flds_i2x_fields), MPAS_LOG_ERR)
    endif
    tagname = trim(seq_flds_x2i_fields)//C_NULL_CHAR
    ierrmb = iMOAB_DefineTagStorage(MPSIID, tagname, tagtype, numco,  tagindex )
    if ( ierrmb == 1 ) then
      call mpas_log_write('cannot define tags for MOAB i2x fields ' // trim(seq_flds_x2i_fields), MPAS_LOG_ERR)
    endif
    ierrmb = iMOAB_SetDoubleTagStorage ( MPSIID, tagname, totalmblr , ent_type, x2i_im(1, 1) )
    if ( ierrmb /= 0 ) then
       call mpas_log_write('cannot set tags for MOAB x2i fields to zero' // trim(seq_flds_x2i_fields), MPAS_LOG_ERR)
    endif

    ! add domain tags
    tagname=trim(seq_flds_dom_fields)//C_NULL_CHAR
    tagtype = 1
    ierrmb = iMOAB_DefineTagStorage(MPSIID, tagname, tagtype, numco,  tagindex )
    if ( ierrmb /= 0 ) then
       call mpas_log_write('Fail to set define dom fields ',  MPAS_LOG_ERR)
    endif

    call ice_domain_moab(MPSIID)
#endif


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
    call MPAS_pool_get_config(domain % configs, "config_column_physics_type", tempCharConfig)
    if (trim(tempCharConfig) == "icepack") then
      call seaice_icepack_coupling_prep(domain)
    else if (trim(tempCharConfig) == "column_package") then
      call seaice_column_coupling_prep(domain)
    endif ! config_column_physics_type

    call MPAS_pool_get_config(domain % configs, "config_use_floe_size_distribution", tempLogicalConfig)
    if (tempLogicalConfig) then
      call mpas_log_write('FloeSizeDistribution coming online soon. Turn FSD off for now.', MPAS_LOG_CRIT)
    endif
!-----------------------------------------------------------------------
!
!   send intial state to driver
!
!-----------------------------------------------------------------------

    call ice_export_mct(i2x_i, errorCode)
    if (errorCode /= 0) then
       call mpas_log_write('Error in ice_export_mct', MPAS_LOG_CRIT)
    endif
#ifdef HAVE_MOAB
    call ice_export_moab(EClock)
#endif

    call t_stopf ('mpassi_mct_init')

    call seq_infodata_PutData( infodata, ice_prognostic=.true.)

!-----------------------------------------------------------------------
!
!   get intial state from driver
!
    call ice_import_mct(x2i_i, errorCode)
    if (errorCode /= 0) then
       call mpas_log_write('Error in ice_import_mct', MPAS_LOG_CRIT)
    endif
!-----------------------------------------------------------------------
#ifdef HAVE_MOAB

#ifdef MOABCOMP
    mpicom_moab = mpicom_i ! save it for run method
    ! loop over all fields in seq_flds_x2i_fields
    call mct_list_init(temp_list ,seq_flds_x2i_fields)
    size_list=mct_list_nitem (temp_list)
    ent_type = 1 ! entity type is cell for ice
    modelStr='ice init'

    do index_list = 1, size_list
      call mct_list_get(mctOStr,index_list,temp_list)
      mct_field = mct_string_toChar(mctOStr)
      tagname= trim(mct_field)//C_NULL_CHAR
      call seq_comm_compare_mb_mct(modelStr, mpicom_moab, x2i_i, mct_field,  MPSIID, tagname, ent_type, difference)
    enddo
    call mct_list_clean(temp_list)
#endif

    call ice_import_moab(Eclock)
#endif


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
!   Initialize ice prescribed mode
!
!----------------------------------------------------------------------------

    call ice_prescribed_init( domain, ICEID, gsmap_i, dom_i )

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
! !IROUTINE: ice_prescribed_init
!
! !INTERFACE:
  subroutine ice_prescribed_init( domain, compid, gsmap, dom )!{{{
!
! !DESCRIPTION:
! Initialize prescribed ice
!
! !USES:
    use shr_pio_mod, only: &
         shr_pio_getiotype, &
         shr_pio_getiosys
    use seq_comm_mct, only: &
         seq_comm_name
    use shr_cal_mod, only: &
         shr_cal_noleap, &
         shr_cal_gregorian
    use shr_strdata_mod, only: &
         shr_strdata_create
    use shr_kind_mod, only: &
         SHR_KIND_IN
    include 'mpif.h'
! !INPUT/OUTPUT PARAMETERS:
    type (domain_type), pointer :: &
         domain
    integer, intent(in) :: &
         compid
    type(mct_gsMap) :: &
         gsmap
    type(mct_gGrid) :: &
         dom
!
! !REVISION HISTORY:
! Author: Adrian K. Turner
!EOP
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    integer(SHR_KIND_IN), parameter :: nFilesMaximum = 400 ! max number of files

    integer, pointer :: stream_year_first   ! first year in stream to use
    integer, pointer :: stream_year_last    ! last year in stream to use
    integer, pointer :: model_year_align    ! align stream_year_first
                                   ! with this model year
    character(len=strKIND), pointer :: stream_fldVarName
    character(len=strKIND), pointer :: stream_fldFileNameIn
    character(len=strKIND) :: stream_fldFileName(nFilesMaximum)
    character(len=strKIND), pointer :: stream_domTvarName
    character(len=strKIND), pointer :: stream_domXvarName
    character(len=strKIND), pointer :: stream_domYvarName
    character(len=strKIND), pointer :: stream_domAreaName
    character(len=strKIND), pointer :: stream_domMaskName
    character(len=strKIND), pointer :: stream_domFileName
    character(len=strKIND), pointer :: stream_mapread
    logical, pointer :: stream_fill ! true if data fill required

    integer :: &
         ierr, &
         iFile, &
         nFiles

    character(len=8) :: &
         fillalgo

    character(len=16) :: &
         inst_name

    integer, pointer :: &
         nCellsSolve

    Integer :: &
         nCells

    character(len=strKIND), pointer :: &
         config_calendar_type

    character(len=strKIND) :: &
         calendar_type

    logical, pointer :: &
         config_use_prescribed_ice

    call MPAS_pool_get_config(domain % configs, "config_use_prescribed_ice", config_use_prescribed_ice)
    if (config_use_prescribed_ice) then

       call mpas_pool_get_config(domain % configs, "config_prescribed_ice_stream_year_first", stream_year_first)
       call mpas_pool_get_config(domain % configs, "config_prescribed_ice_stream_year_last", stream_year_last)
       call mpas_pool_get_config(domain % configs, "config_prescribed_ice_model_year_align", model_year_align)
       call mpas_pool_get_config(domain % configs, "config_prescribed_ice_stream_fldvarname", stream_fldVarName)
       call mpas_pool_get_config(domain % configs, "config_prescribed_ice_stream_fldfilename", stream_fldFileNameIn)
       call mpas_pool_get_config(domain % configs, "config_prescribed_ice_stream_domtvarname", stream_domTvarName)
       call mpas_pool_get_config(domain % configs, "config_prescribed_ice_stream_domxvarname", stream_domXvarName)
       call mpas_pool_get_config(domain % configs, "config_prescribed_ice_stream_domyvarname", stream_domYvarName)
       call mpas_pool_get_config(domain % configs, "config_prescribed_ice_stream_domareaname", stream_domAreaName)
       call mpas_pool_get_config(domain % configs, "config_prescribed_ice_stream_dommaskname", stream_domMaskName)
       call mpas_pool_get_config(domain % configs, "config_prescribed_ice_stream_domfilename", stream_domFileName)
       call mpas_pool_get_config(domain % configs, "config_prescribed_ice_stream_mapread", stream_mapread)
       call mpas_pool_get_config(domain % configs, "config_prescribed_ice_stream_fill", stream_fill)

       ! file array
       stream_fldFileName(1) = trim(stream_fldFileNameIn)
       do iFile = 2, nFilesMaximum
          stream_fldFileName(iFile) = ' '
       end do

       ! file number
       nFiles = 0
       do iFile = 1, nFilesMaximum
          if (stream_fldFileName(iFile) /= ' ') nFiles = nFiles + 1
       end do

       ! Read shr_strdata_nml namelist
       if (stream_fill) then
          fillalgo='nn'
       else
          fillalgo='none'
       endif

       call MPAS_log_write(' ')
       call MPAS_log_write('This is the prescribed ice coverage option.')
       call MPAS_log_write('  stream_year_first  = $i', intArgs=(/stream_year_first/))
       call MPAS_log_write('  stream_year_last   = $i', intArgs=(/stream_year_last/))
       call MPAS_log_write('  model_year_align   = $i', intArgs=(/model_year_align/))
       call MPAS_log_write('  stream_fldVarName  = '//trim(stream_fldVarName))
       do iFile = 1, nFiles
          call MPAS_log_write('  stream_fldFileName = $i: '//trim(stream_fldFileName(iFile)), intArgs=(/iFile/))
       end do
       call MPAS_log_write('  stream_domTvarName = '//trim(stream_domTvarName))
       call MPAS_log_write('  stream_domXvarName = '//trim(stream_domXvarName))
       call MPAS_log_write('  stream_domYvarName = '//trim(stream_domYvarName))
       call MPAS_log_write('  stream_domFileName = '//trim(stream_domFileName))
       call MPAS_log_write('  stream_domAreaName = '//trim(stream_domAreaName))
       call MPAS_log_write('  stream_domMaskName = '//trim(stream_domMaskName))
       call MPAS_log_write('  stream_mapread     = '//trim(stream_mapread))
       call MPAS_log_write('  stream_fillalgo    = '//trim(fillalgo))
       call MPAS_log_write(' ')

       inst_name = seq_comm_name(compid)

       call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nCellsSolve", nCellsSolve)
       call MPI_Allreduce(nCellsSolve, nCells, 1, MPI_Integer, MPI_SUM, domain % dminfo % comm, ierr)


       call MPAS_pool_get_config(domain % configs, "config_calendar_type", config_calendar_type)
       if (trim(config_calendar_type) == "gregorian") then
          calendar_type = shr_cal_gregorian
       else if (trim(config_calendar_type) == "noleap") then
          calendar_type = shr_cal_noleap
       else
          call MPAS_log_write("ice_prescribed_init: Unknown calendar type: "//trim(config_calendar_type))
       endif

       call shr_strdata_create( &
            sdat, &
            name = "prescribed_ice", &
            mpicom = domain % dminfo % comm, &
            compid = compid, &
            gsmap = gsmap, &
            ggrid = dom, &
            nxg = nCells, &
            nyg = 1, &
            yearFirst = stream_year_first, &
            yearLast = stream_year_last, &
            yearAlign = model_year_align, &
            offset = 0, &
            domFilePath = '', &
            domFileName = trim(stream_domFileName), &
            domTvarName = trim(stream_domTvarName), &
            domXvarName = trim(stream_domXvarName), &
            domYvarName = trim(stream_domYvarName), &
            domAreaName = trim(stream_domAreaName), &
            domMaskName = trim(stream_domMaskName), &
            filePath = '', &
            filename = stream_fldFileName(1:nFiles), &
            fldListFile = trim(stream_fldVarName), &
            fldListModel = trim(stream_fldVarName), &
            pio_subsystem = shr_pio_getiosys(inst_name), &
            pio_iotype = shr_pio_getiotype(inst_name), &
            fillalgo = trim(fillalgo), &
            calendar = trim(calendar_type), &
            mapread = trim(stream_mapread))

    endif ! config_use_prescribed_ice

  end subroutine ice_prescribed_init

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
! Variables related to E3SM coupler
      integer :: ymd, tod, ihour, iminute, isecond
      integer :: iyear, imonth, iday, curr_ymd, curr_tod
      integer :: shrloglev, shrlogunit

! Variable related to MPASSI
      type (block_type), pointer :: block
      type (MPAS_Pool_type), pointer :: &
            shortwave, &
            berg_forcing

      real (kind=RKIND) :: current_wallclock_time
      type (MPAS_Time_Type) :: currTime
      character(len=StrKIND) :: timeStamp, streamName, WCstring
      type (MPAS_timeInterval_type) :: timeStep
      integer :: ierr, streamDirection, iam
      logical :: streamActive, debugOn
      logical, pointer :: &
            config_write_output_on_startup, &
            config_scale_dib_by_removed_ice_runoff
      logical, save :: first=.true.
      character (len=StrKIND), pointer :: &
            config_restart_timestamp_name, &
            config_column_physics_type
      real(kind=RKIND), pointer :: &
            dayOfNextShortwaveCalculation, & ! needed for CESM like coupled simulations
            runningMeanRemovedIceRunoff ! the area integrated, running mean of removed ice runoff from the ocean

#ifdef MOABCOMP
      real(r8)                 :: difference
      type(mct_list) :: temp_list
      integer :: size_list, index_list, ent_type
      type(mct_string)    :: mctOStr  !
      character(CXX) :: mct_field, modelStr, tagname
#endif


      iam = domain % dminfo % my_proc_id

      debugOn = .false.
#ifdef MPAS_DEBUG
      debugOn = .true.
#endif

      ! Set MPAS Log module instance
      mpas_log_info => domain % logInfo
      if (debugOn) call mpas_log_write("=== Beginning ice_run_mct ===")

      call MPAS_pool_get_config(domain % configs, 'config_restart_timestamp_name', config_restart_timestamp_name)
      call MPAS_pool_get_config(domain % configs, "config_column_physics_type", config_column_physics_type)
      call MPAS_pool_get_config(domain % configs, "config_scale_dib_by_removed_ice_runoff", &
                                config_scale_dib_by_removed_ice_runoff)

      ! Setup log information.
      call shr_file_getLogUnit (shrlogunit)
      call shr_file_getLogLevel(shrloglev)
      call shr_file_setLogUnit (iceLogUnit)

      ! reinitialize fluxes
      if (trim(config_column_physics_type) == "icepack") then
         call seaice_icepack_reinitialize_fluxes(domain)
      else if (trim(config_column_physics_type) == "column_package") then
         call seaice_column_reinitialize_fluxes(domain)
      endif ! config_column_physics_type

      ! Import state from coupler
      call ice_import_mct(x2i_i, ierr)
#ifdef HAVE_MOAB
#ifdef MOABCOMP
    ! loop over all fields in seq_flds_x2i_fields
    call mct_list_init(temp_list ,seq_flds_x2i_fields)
    size_list=mct_list_nitem (temp_list)
    ent_type = 1 ! entity type is cell for ice
    modelStr='ice run'

    do index_list = 1, size_list
      call mct_list_get(mctOStr,index_list,temp_list)
      mct_field = mct_string_toChar(mctOStr)
      tagname= trim(mct_field)//C_NULL_CHAR
      call seq_comm_compare_mb_mct(modelStr, mpicom_moab, x2i_i, mct_field,  MPSIID, tagname, ent_type, difference)
    enddo
    call mct_list_clean(temp_list)
#endif
      call ice_import_moab(Eclock)
#endif

      ! Post coupling calls
      block => domain % blocklist

      if (config_scale_dib_by_removed_ice_runoff) then
         ! independent of space so should be no need to loop over blocks
         call MPAS_pool_get_subpool(block % structs, "berg_forcing", berg_forcing)
         call MPAS_pool_get_array(berg_forcing, "runningMeanRemovedIceRunoff", &
                                  runningMeanRemovedIceRunoff)
         call seq_infodata_GetData(infodata, rmean_rmv_ice_runoff=runningMeanRemovedIceRunoff )
      end if

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

        ! get prescribed ice coverage
        call ice_prescribed_run(domain, currTime)

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
         call seaice_forcing_write_restart_times(domain)
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
#ifdef HAVE_MOAB
      call ice_export_moab(EClock)
#endif
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
! !IROUTINE: ice_prescribed_run
!
! !INTERFACE:
  subroutine ice_prescribed_run( domain, currTime )!{{{
!
! !DESCRIPTION:
! Set ice coverage for prescribed ice mode
!
! !USES:
    use shr_strdata_mod, only: &
         shr_strdata_advance
! !INPUT/OUTPUT PARAMETERS:
    type (domain_type), pointer :: &
         domain
    type (MPAS_Time_Type), intent(in) :: &
         currTime
!
! !REVISION HISTORY:
! Author: Adrian K. Turner
!EOP
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

    logical, pointer :: &
         config_use_prescribed_ice

    integer :: &
         mDateIn, &
         secIn, &
         YYYY, MM, DD, DoY, H, M, S, &
         ierr

    integer :: &
         n, &
         iCell

    type (block_type), pointer :: &
         block_ptr

    type (mpas_pool_type), pointer :: &
         meshPool, &
         prescribedIcePool

    integer, pointer :: &
         nCellsSolve

    real(kind=RKIND), dimension(:), pointer :: &
         iceCoverage

    type (field1DReal), pointer :: &
         iceCoverageField

    call MPAS_pool_get_config(domain % configs, "config_use_prescribed_ice", config_use_prescribed_ice)
    if (config_use_prescribed_ice) then

       call MPAS_get_time(currTime, YYYY, MM, DD, DoY, H, M, S, ierr)
       mDateIn = YYYY * 10000 + MM * 100 + DD
       secIn = H * 3600 + M * 60 + S

       call shr_strdata_advance(&
            sdat, &
            mDateIn, &
            SecIn, &
            domain % dminfo % comm, &
            'mpassi_pice')

       n = 0
       block_ptr => domain % blocklist
       do while(associated(block_ptr))

          call MPAS_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
          call MPAS_pool_get_subpool(block_ptr % structs, 'prescribed_ice', prescribedIcePool)

          call MPAS_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

          call MPAS_pool_get_array(prescribedIcePool, 'iceCoverage', iceCoverage)

          iceCoverage(:) = 0.0_RKIND

          do iCell = 1, nCellsSolve
             n = n + 1
             iceCoverage(iCell) = sdat % avs(1) % rAttr(1,n)
          enddo ! iCell

          block_ptr => block_ptr % next
       end do

       call MPAS_pool_get_subpool(domain % blocklist % structs, 'prescribed_ice', prescribedIcePool)
       call MPAS_pool_get_field(prescribedIcePool, 'iceCoverage', iceCoverageField)
       call MPAS_dmpar_exch_halo_field(iceCoverageField)

    endif

  end subroutine ice_prescribed_run

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

#ifdef HAVE_MOAB
   ! deallocate moab fields array
    deallocate (i2x_im)
    deallocate (x2i_im)
#endif
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

#ifdef HAVE_MOAB
!***********************************************************************
!BOP
! !IROUTINE: ice_domain_moab
! !INTERFACE:

 subroutine ice_domain_moab( mbid )!{{{

  use iso_c_binding, only : C_NULL_CHAR
  use iMOAB, only: iMOAB_GetMeshInfo
! !DESCRIPTION:
!  This routine sets up the MOAB domain-related info for MPASSI
!
! !REVISION HISTORY:
!  2025-01-23: initial version is copy from ice_domain_mct
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

    type (mpas_pool_type), pointer :: meshPool, oceanCouplingPool

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
    tagname = "mask"//C_NULL_CHAR
    ent_type=1
    ierr = iMOAB_SetDoubleTagStorage ( mbid, tagname, arrSize , ent_type, data)
    tagname = "frac"//C_NULL_CHAR
    ent_type=1
    ierr = iMOAB_SetDoubleTagStorage ( mbid, tagname, arrSize , ent_type, data)

    deallocate(data)

!-----------------------------------------------------------------------
!EOC
  end subroutine ice_domain_moab !}}}
#endif

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
!    o  q        -- ocn frazil heat flux(+) / melt potential(-)
!    o  frazil   -- ocn frazil mass flux
!    o  bcphidry -- Black Carbon hydrophilic dry deposition flux
!    o  bcphodry -- Black Carbon hydrophobic dry deposition flux
!    o  bcphiwet -- Black Carbon hydrophilic wet deposition flux
!    o  ocphidry -- Organic Carbon hydrophilic dry deposition flux
!    o  ocphodry -- Organic Carbon hydrophobic dry deposition flux
!    o  ocphiwet -- Organic Carbon hydrophilic wet deposition flux
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
      frazilMassFlux, &
      frazilMassFluxRev

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
      config_use_column_biogeochemistry, &
      config_couple_biogeochemistry_fields

   character(len=strKIND), pointer ::   &
      config_column_physics_type,       &
      config_thermodynamics_type,       &
      config_ocean_surface_type

   type (field1DReal), pointer ::       &
      seaSurfaceTemperatureField,       &
      seaSurfaceSalinityField,          &
      seaFreezingTemperatureField,      &
      freezingMeltingPotentialField,    &
      frazilMassAdjustField,            &
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
      atmosDustFluxField,               &
      atmosWetDustFluxField,            &
      atmosDryDustFluxField

   real (kind=RKIND), dimension(:), pointer :: &
      seaSurfaceTemperature,       &
      seaSurfaceSalinity,          &
      seaFreezingTemperature,      &
      freezingMeltingPotential,    &
      frazilMassAdjust,            &
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
      carbonToNitrogenRatioDON,    &
      DOCPoolFractions

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
      atmosDustFlux,               &
      atmosWetDustFlux,            &
      atmosDryDustFlux

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
      call mpas_pool_get_config(configs, "config_column_physics_type", config_column_physics_type)
      call mpas_pool_get_config(configs, "config_thermodynamics_type", config_thermodynamics_type)
      call mpas_pool_get_config(configs, "config_ocean_surface_type", config_ocean_surface_type)
      call mpas_pool_get_config(configs, "config_use_aerosols", config_use_aerosols)
      call mpas_pool_get_config(configs, "config_use_modal_aerosols", config_use_modal_aerosols)
      call mpas_pool_get_config(configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)
      call mpas_pool_get_config(configs, "config_couple_biogeochemistry_fields", config_couple_biogeochemistry_fields)
      call mpas_pool_get_config(configs, "config_use_zaerosols", config_use_zaerosols)

      call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block_ptr % structs, 'ocean_coupling', oceanCoupling)
      call mpas_pool_get_subpool(block_ptr % structs, 'atmos_coupling', atmosCoupling)

      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

      call mpas_pool_get_array(oceanCoupling, 'seaSurfaceTemperature', seaSurfaceTemperature)
      call mpas_pool_get_array(oceanCoupling, 'seaSurfaceSalinity', seaSurfaceSalinity)
      call mpas_pool_get_array(oceanCoupling, 'seaFreezingTemperature', seaFreezingTemperature)
      call mpas_pool_get_array(oceanCoupling, 'freezingMeltingPotential', freezingMeltingPotential)
      call mpas_pool_get_array(oceanCoupling, 'frazilMassAdjust', frazilMassAdjust)
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

      call mpas_pool_get_subpool(block_ptr % structs, 'biogeochemistry', biogeochemistry)
      if (config_use_zaerosols) then
         call mpas_pool_get_array(biogeochemistry, 'atmosBlackCarbonFlux', atmosBlackCarbonFlux)
         call mpas_pool_get_array(biogeochemistry, 'atmosDustFlux', atmosDustFlux)
         call mpas_pool_get_array(biogeochemistry, "atmosWetDustFlux", atmosWetDustFlux)
         call mpas_pool_get_array(biogeochemistry, "atmosDryDustFlux", atmosDryDustFlux)
      endif
      if (config_use_column_biogeochemistry .and. config_couple_biogeochemistry_fields) then

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
         call mpas_pool_get_array(biogeochemistry, 'carbonToNitrogenRatioAlgae', carbonToNitrogenRatioAlgae)
         call mpas_pool_get_array(biogeochemistry, 'carbonToNitrogenRatioDON', carbonToNitrogenRatioDON)
         call mpas_pool_get_array(biogeochemistry, 'DOCPoolFractions', DOCPoolFractions)

      endif

      do i = 1, nCellsSolve
        n = n + 1
        seaSurfaceTemperature(i)       = x2i_i % rAttr(index_x2i_So_t, n)
        seaSurfaceSalinity(i)          = x2i_i % rAttr(index_x2i_So_s, n)

        if (trim(config_column_physics_type) == "icepack") then
           seaFreezingTemperature(i) = icepack_sea_freezing_temperature(seaSurfaceSalinity(i))
        else if (trim(config_column_physics_type) == "column_package") then
           seaFreezingTemperature(i) = colpkg_sea_freezing_temperature(seaSurfaceSalinity(i))
        endif ! config_column_physics_type

        uOceanVelocity(i)              = x2i_i % rAttr(index_x2i_So_u, n)
        vOceanVelocity(i)              = x2i_i % rAttr(index_x2i_So_v, n)
        seaSurfaceTiltU(i)             = x2i_i % rAttr(index_x2i_So_dhdx, n)
        seaSurfaceTiltV(i)             = x2i_i % rAttr(index_x2i_So_dhdy, n)

        if (trim(config_ocean_surface_type) == "free") then ! free surface (MPAS-O)

           ! freezingMeltingPotential(i) is the ocean energy associated with frazil formation
           ! when it is positive and frazilMassFlux is positive. Conversely, freezingMeltingPotential(i)
           ! is negative when there is the melting potential in which case frazilMassFlux is zero.

           freezingMeltingPotential(i) = x2i_i % rAttr(index_x2i_Fioo_q, n)

           frazilMassFlux              = x2i_i % rAttr(index_x2i_Fioo_frazil, n)

           ! Now determine the sea ice mass associated with the frazil heat flux given when
           ! freezingMeltingPotential(i) is positive. This produces a revised mass flux, given
           ! in frazilMassFluxRev for the given sea surface salinity. The resulting difference
           ! is assigned to frazilMassAdjust(i) which is exported to the ocean in the subsequent
           ! coupling step as a freshwater and salt flux. This step is required to balance mass
           ! and heat with the ocean.

           call frazil_mass(freezingMeltingPotential(i), frazilMassFluxRev, seaSurfaceSalinity(i))

           frazilMassAdjust(i) = frazilMassFlux-frazilMassFluxRev

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

        ! set aerosols, if configured
        if (config_use_zaerosols .or. config_use_column_biogeochemistry) then
           call mpas_pool_get_array(biogeochemistry, 'oceanZAerosolConc', oceanZAerosolConc)
           oceanZAerosolConc(1,i)        = 0.0_RKIND !x2i_i % rAttr(index_x2i_So_zaer1, n)
           oceanZAerosolConc(2,i)        = 0.0_RKIND !x2i_i % rAttr(index_x2i_So_zaer2, n)
           oceanZAerosolConc(3,i)        = 0.0_RKIND !x2i_i % rAttr(index_x2i_So_zaer3, n)
           oceanZAerosolConc(4,i)        = 0.0_RKIND !x2i_i % rAttr(index_x2i_So_zaer4, n)
           oceanZAerosolConc(5,i)        = 0.0_RKIND !x2i_i % rAttr(index_x2i_So_zaer5, n)
           oceanZAerosolConc(6,i)        = 0.0_RKIND !x2i_i % rAttr(index_x2i_So_zaer6, n)
        end if
        if (config_use_zaerosols) then
           if (config_use_modal_aerosols) then
              atmosBlackCarbonFlux(1,i) = x2i_i % rAttr(index_x2i_Faxa_bcphodry, n) &
                                        + x2i_i % rAttr(index_x2i_Faxa_bcphidry, n)
              atmosBlackCarbonFlux(2,i) = x2i_i % rAttr(index_x2i_Faxa_bcphiwet, n)
              !combine wet and dry dust
              atmosDustFlux(1,i) = x2i_i % rAttr(index_x2i_Faxa_dstwet1, n) &
                                 + x2i_i % rAttr(index_x2i_Faxa_dstdry1, n)
              atmosDustFlux(2,i) = x2i_i % rAttr(index_x2i_Faxa_dstwet2, n) &
                                 + x2i_i % rAttr(index_x2i_Faxa_dstdry2, n)
              atmosDustFlux(3,i) = x2i_i % rAttr(index_x2i_Faxa_dstwet3, n) &
                                 + x2i_i % rAttr(index_x2i_Faxa_dstdry3, n)
              atmosDustFlux(4,i) = x2i_i % rAttr(index_x2i_Faxa_dstwet4, n) &
                                 + x2i_i % rAttr(index_x2i_Faxa_dstdry4, n)

              atmosWetDustFlux(1,i) = x2i_i % rAttr(index_x2i_Faxa_dstwet1, n)
              atmosWetDustFlux(2,i) = x2i_i % rAttr(index_x2i_Faxa_dstwet2, n)
              atmosWetDustFlux(3,i) = x2i_i % rAttr(index_x2i_Faxa_dstwet3, n)
              atmosWetDustFlux(4,i) = x2i_i % rAttr(index_x2i_Faxa_dstwet4, n)

              atmosDryDustFlux(1,i) = x2i_i % rAttr(index_x2i_Faxa_dstdry1, n)
              atmosDryDustFlux(2,i) = x2i_i % rAttr(index_x2i_Faxa_dstdry2, n)
              atmosDryDustFlux(3,i) = x2i_i % rAttr(index_x2i_Faxa_dstdry3, n)
              atmosDryDustFlux(4,i) = x2i_i % rAttr(index_x2i_Faxa_dstdry4, n)
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

              atmosDryDustFlux(1,i) = x2i_i % rAttr(index_x2i_Faxa_dstdry1, n)
              atmosDryDustFlux(2,i) = x2i_i % rAttr(index_x2i_Faxa_dstdry2, n)
              atmosDryDustFlux(3,i) = x2i_i % rAttr(index_x2i_Faxa_dstdry3, n)
              atmosDryDustFlux(4,i) = x2i_i % rAttr(index_x2i_Faxa_dstdry4, n)

              atmosWetDustFlux(1,i) = x2i_i % rAttr(index_x2i_Faxa_dstwet1, n)
              atmosWetDustFlux(2,i) = x2i_i % rAttr(index_x2i_Faxa_dstwet2, n)
              atmosWetDustFlux(3,i) = x2i_i % rAttr(index_x2i_Faxa_dstwet3, n)
              atmosWetDustFlux(4,i) = x2i_i % rAttr(index_x2i_Faxa_dstwet4, n)
           endif
        endif
        ! import biogeochemistry fields, if configured
        if (config_use_column_biogeochemistry .and. config_couple_biogeochemistry_fields) then
           oceanAlgaeConc(1,i)           = x2i_i % rAttr(index_x2i_So_algae1, n)
           oceanAlgaeConc(2,i)           = x2i_i % rAttr(index_x2i_So_algae2, n)
           oceanAlgaeConc(3,i)           = 0.0_RKIND !x2i_i % rAttr(index_x2i_So_algae3, n)
           oceanDOCConc(1,i)             = x2i_i % rAttr(index_x2i_So_doc1, n) * DOCPoolFractions(1)
           oceanDOCConc(2,i)             = x2i_i % rAttr(index_x2i_So_doc2, n) * DOCPoolFractions(2)
           oceanDOCConc(3,i)             = 0.0_RKIND
           oceanDICConc(1,i)             = x2i_i % rAttr(index_x2i_So_dic1, n)
           oceanDONConc(1,i)             = x2i_i % rAttr(index_x2i_So_doc3, n) * DOCPoolFractions(3) &
              /carbonToNitrogenRatioDON(1)
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
        endif
      end do

!-----------------------------------------------------------------------
!
!  unit conversions and any manipulation of coupled fields
!
!-----------------------------------------------------------------------
      do i = 1, nCellsSolve

        seaSurfaceTemperature(i)  = seaSurfaceTemperature(i) - seaiceFreshWaterFreezingPoint

        if (config_use_column_biogeochemistry .and. config_couple_biogeochemistry_fields) then
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
   call mpas_pool_get_field(oceanCoupling, 'frazilMassAdjust', frazilMassAdjustField)
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

   if (config_use_zaerosols .or. config_use_column_biogeochemistry) &
      call mpas_pool_get_field(biogeochemistry, 'oceanZAerosolConc', oceanZAerosolConcField)

   if (config_use_zaerosols) then
      call mpas_pool_get_field(biogeochemistry, "atmosBlackCarbonFlux", atmosBlackCarbonFluxField)
      call mpas_pool_get_field(biogeochemistry, "atmosDustFlux", atmosDustFluxField)
      call mpas_pool_get_field(biogeochemistry, "atmosWetDustFlux", atmosWetDustFluxField)
      call mpas_pool_get_field(biogeochemistry, "atmosDryDustFlux", atmosDryDustFluxField)
   endif
   if (config_use_column_biogeochemistry .and. config_couple_biogeochemistry_fields) then
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
   endif

   call mpas_dmpar_exch_halo_field(seaSurfaceTemperatureField)
   call mpas_dmpar_exch_halo_field(seaSurfaceSalinityField)
   call mpas_dmpar_exch_halo_field(seaFreezingTemperatureField)
   call mpas_dmpar_exch_halo_field(freezingMeltingPotentialField)
   call mpas_dmpar_exch_halo_field(frazilMassAdjustField)
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
   if (config_use_zaerosols .or. config_use_column_biogeochemistry) &
      call mpas_dmpar_exch_halo_field(oceanZAerosolConcField)

   if (config_use_zaerosols) then
      call mpas_dmpar_exch_halo_field(atmosBlackCarbonFluxField)
      call mpas_dmpar_exch_halo_field(atmosDustFluxField)
      call mpas_dmpar_exch_halo_field(atmosDryDustFluxField)
      call mpas_dmpar_exch_halo_field(atmosWetDustFluxField)
   endif
   if (config_use_column_biogeochemistry .and. config_couple_biogeochemistry_fields) then
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
   endif

! REVISION HISTORY:
! Revised Andrew Roberts May 2021
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
      oceanCoupling,    &
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
      config_use_zaerosols,              &
      config_couple_biogeochemistry_fields, &
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
      frazilMassAdjust,            &
      bergFreshwaterFlux,          &
      bergLatentHeatFlux,          &
      oceanNitrateFlux,            &
      oceanSilicateFlux,           &
      oceanAmmoniumFlux,           &
      oceanDMSFlux,                &
      oceanDMSPpFlux,              &
      oceanDMSPdFlux,              &
      oceanHumicsFlux,             &
      oceanDustIronFlux,           &
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
      call MPAS_pool_get_config(configs, "config_use_zaerosols", config_use_zaerosols)
      call MPAS_pool_get_config(configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)
      call mpas_pool_get_config(configs, "config_couple_biogeochemistry_fields", config_couple_biogeochemistry_fields)
      call MPAS_pool_get_config(configs, "config_use_column_shortwave", config_use_column_shortwave)
      call MPAS_pool_get_config(configs, "config_use_data_icebergs", config_use_data_icebergs)

      call MPAS_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
      call MPAS_pool_get_subpool(block_ptr % structs, "tracers_aggregate", tracersAggregate)
      call MPAS_pool_get_subpool(block_ptr % structs, "velocity_solver", velocitySolver)
      call MPAS_pool_get_subpool(block_ptr % structs, "shortwave", shortwave)
      call MPAS_pool_get_subpool(block_ptr % structs, 'atmos_coupling', atmosCoupling)
      call MPAS_pool_get_subpool(block_ptr % structs, 'ocean_coupling', oceanCoupling)
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

      call MPAS_pool_get_array(oceanCoupling, 'frazilMassAdjust', frazilMassAdjust)

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

      if (config_use_column_biogeochemistry .and. config_couple_biogeochemistry_fields) then
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

      if (config_use_zaerosols .and. config_couple_biogeochemistry_fields) &
         call mpas_pool_get_array(biogeochemistry, 'oceanDustIronFlux', oceanDustIronFlux)

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
            i2x_i % rAttr(index_i2x_Si_u10withgusts,n) = atmosReferenceSpeed10m(i)
            i2x_i % rAttr(index_i2x_Si_tref  ,n)  = atmosReferenceTemperature2m(i)
            i2x_i % rAttr(index_i2x_Si_qref  ,n)  = atmosReferenceHumidity2m(i)
            i2x_i % rAttr(index_i2x_Si_snowh ,n)  = snowVolumeCell(i) / ailohi
            i2x_i % rAttr(index_i2x_Si_ithick,n)  = iceVolumeCell(i) / ailohi

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

            ! i/o fluxes computed by ice, as well as additional freshwater and salt calculated at the last
            ! coupling import and needed to grow sea ice from frazil passed from the ocean model in the
            ! field frazilMassAdjust.
            i2x_i % rAttr(index_i2x_Fioi_melth,n) = oceanHeatFlux(i)
            i2x_i % rAttr(index_i2x_Fioi_swpen,n) = oceanShortwaveFlux(i)
            i2x_i % rAttr(index_i2x_Fioi_meltw,n) = oceanFreshWaterFlux(i) + frazilMassAdjust(i)/ailohi
            i2x_i % rAttr(index_i2x_Fioi_salt ,n) = oceanSaltFlux(i) + seaiceReferenceSalinity*0.001_RKIND*frazilMassAdjust(i)/ailohi
            i2x_i % rAttr(index_i2x_Fioi_taux ,n) = tauxo
            i2x_i % rAttr(index_i2x_Fioi_tauy ,n) = tauyo

            ! export biogeochemistry fields, if configured
            if (config_use_column_biogeochemistry .and. config_couple_biogeochemistry_fields) then
               ! convert from mmol N/m^3 to mmol C/m^3
               i2x_i % rAttr(index_i2x_Fioi_algae1,n) = oceanAlgaeFlux(1,i) * carbonToNitrogenRatioAlgae(1)
               i2x_i % rAttr(index_i2x_Fioi_algae2,n) = oceanAlgaeFlux(2,i) * carbonToNitrogenRatioAlgae(2)
               i2x_i % rAttr(index_i2x_Fioi_algae3,n) = oceanAlgaeFlux(3,i) * carbonToNitrogenRatioAlgae(3)
               i2x_i % rAttr(index_i2x_Fioi_doc1  ,n) = oceanDOCFlux(1,i)
               i2x_i % rAttr(index_i2x_Fioi_doc2  ,n) = oceanDOCFlux(2,i)
               i2x_i % rAttr(index_i2x_Fioi_doc3  ,n) = oceanDONFlux(1,i) * carbonToNitrogenRatioDON(1)
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
            ! export dust, if configured kg/m2/s of dust
            if (config_use_zaerosols .and. config_couple_biogeochemistry_fields) &
               i2x_i % rAttr(index_i2x_Fioi_dust1 ,n) = oceanDustIronFlux(i)
         endif
      enddo

      block_ptr => block_ptr % next
   enddo

! REVISION HISTORY:
! Revised Andrew Roberts May 2021
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
! !IROUTINE: frazil_mass
!
! !INTERFACE
   subroutine frazil_mass(freezingPotential, frazilMassFlux, seaSurfaceSalinity)
!
! !DESCRIPTION:
! Calculate frazil mass based on on the sea surface salinity, and frazil heat flux
! from the ocean, otherwise referred to as to as the freeze-melt potential. When
! freezingPotential is positive, it gives the heat flux, according to the ocean model
! associated with the frazil mass passed from the ocean.  This function calculates
! the frazil mass based on the freezingPotential according to sea ice model thermodynamics,
!
! !USES:
      use ice_mushy_physics, only:  &
         liquidus_temperature_mush, &
         enthalpy_mush

      use icepack_intfc, only:  &
         icepack_liquidus_temperature, &
         icepack_enthalpy_mush

! !INPUT PARAMETERS:
      real (kind=RKIND),      intent(in)   :: freezingPotential
      real (kind=RKIND),      intent(in)   :: seaSurfaceSalinity

! !OUTPUT PARAMETERS:
      real (kind=RKIND),      intent(out)  :: frazilMassFlux

      character(len=strKIND), pointer   :: config_thermodynamics_type
      character(len=strKIND), pointer   :: config_column_physics_type

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

   call MPAS_pool_get_config(domain % configs, "config_thermodynamics_type", config_thermodynamics_type)
   call MPAS_pool_get_config(domain % configs, "config_column_physics_type", config_column_physics_type)

   if (freezingPotential > 0.0_RKIND) then

      if (trim(config_thermodynamics_type) == "mushy") then  ! mushy
         if (seaSurfaceSalinity > 2.0_RKIND * seaiceFrazilSalinityReduction) then
             Si0new = seaSurfaceSalinity - seaiceFrazilSalinityReduction
         else
             Si0new = seaSurfaceSalinity**2 / (4.0_RKIND*seaiceFrazilSalinityReduction)
         endif

         if (trim(config_column_physics_type) == "icepack") then
            Ti = icepack_liquidus_temperature(Si0new/seaiceFrazilIcePorosity)
            qi0new = icepack_enthalpy_mush(Ti, Si0new)
         else if (trim(config_column_physics_type) == "column_package") then
            Ti = liquidus_temperature_mush(Si0new/seaiceFrazilIcePorosity)
            qi0new = enthalpy_mush(Ti, Si0new)
         endif ! config_column_physics_type

      else
         qi0new = -seaiceDensityIce*seaiceLatentHeatMelting
      endif    ! ktherm

      frazilMassFlux = -freezingPotential*seaiceDensityIce/qi0new

   else

      frazilMassFlux = 0.0_RKIND

   endif

! REVISION HISTORY:
! Revised Andrew Roberts May 2021
!-----------------------------------------------------------------------
!EOC

   end subroutine frazil_mass!}}}

!***********************************************************************
!BOP
!
! !IROUTINE: datetime
!
! !INTERFACE:
   subroutine datetime(cdate, ctime)
!
! !DESCRIPTION:
! Calculate current date and time for metadata history
!
! !USES:

! !INPUT PARAMETERS:
      character(len=*), intent(out) :: cdate
      character(len=*), intent(out) :: ctime

! !OUTPUT PARAMETERS:
!EOP
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
      character(len=strKIND) :: date
      character(len=strKIND) :: time
      character(len=strKIND) :: zone
      integer, dimension(8) :: values

      call date_and_time (date, time, zone, values)

      cdate = ""
      cdate(1:2) = date(5:6)
      cdate(3:3) = '/'
      cdate(4:5) = date(7:8)
      cdate(6:6) = '/'
      cdate(7:8) = date(3:4)

      ctime = ""
      ctime(1:2) = time(1:2)
      ctime(3:3) = ':'
      ctime(4:5) = time(3:4)
      ctime(6:6) = ':'
      ctime(7:8) = time(5:6)

!-----------------------------------------------------------------------
!EOC

   end subroutine datetime!}}}

#ifdef HAVE_MOAB
   subroutine ice_export_moab(EClock)
!BOP

   ! !DESCRIPTION:
   !  This routine calls the routines necessary to send MPASSI fields to MOAB coupler
   !
     use iMOAB,  only       : iMOAB_SetDoubleTagStorage, iMOAB_WriteMesh
   !EOP
   !BOC
   type(ESMF_Clock) , intent(inout) :: EClock    ! Input synchronization clock from driver
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
         oceanCoupling,    &
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
         config_use_zaerosols,              &
         config_use_column_shortwave,       &
         config_use_data_icebergs,          &
         config_couple_biogeochemistry_fields

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
         frazilMassAdjust,            &
         bergFreshwaterFlux,          &
         bergLatentHeatFlux,          &
         oceanNitrateFlux,            &
         oceanSilicateFlux,           &
         oceanAmmoniumFlux,           &
         oceanDMSFlux,                &
         oceanDMSPpFlux,              &
         oceanDMSPdFlux,              &
         oceanHumicsFlux,             &
         oceanDustIronFlux,           &
         carbonToNitrogenRatioAlgae,  &
         carbonToNitrogenRatioDON

      real (kind=RKIND), dimension(:,:), pointer :: &
         oceanAlgaeFlux,              &
         oceanDOCFlux,                &
         oceanDICFlux,                &
         oceanDONFlux,                &
         oceanParticulateIronFlux,    &
         oceanDissolvedIronFlux

      integer :: ent_type, ierr, cur_ice_stepno
      character(len=32), parameter :: sub = 'ice_export_moab'

      character(len=100) :: outfile, wopts, localmeshfile, lnum
      character(CXX) :: tagname
   !-----------------------------------------------------------------------
      call shr_file_setLogUnit (iceLogUnit)
      n = 0
      i2x_im(: ,:) = 0.0_RKIND
      block_ptr => domain % blocklist
      do while(associated(block_ptr))

         configs => block_ptr % configs
         call MPAS_pool_get_config(configs, "config_rotate_cartesian_grid", config_rotate_cartesian_grid)
         call MPAS_pool_get_config(configs, "config_use_topo_meltponds", config_use_topo_meltponds)
         call MPAS_pool_get_config(configs, "config_use_zaerosols", config_use_zaerosols)
         call MPAS_pool_get_config(configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)
         call mpas_pool_get_config(configs, "config_couple_biogeochemistry_fields", config_couple_biogeochemistry_fields)
         call MPAS_pool_get_config(configs, "config_use_column_shortwave", config_use_column_shortwave)
         call MPAS_pool_get_config(configs, "config_use_data_icebergs", config_use_data_icebergs)

         call MPAS_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
         call MPAS_pool_get_subpool(block_ptr % structs, "tracers_aggregate", tracersAggregate)
         call MPAS_pool_get_subpool(block_ptr % structs, "velocity_solver", velocitySolver)
         call MPAS_pool_get_subpool(block_ptr % structs, "shortwave", shortwave)
         call MPAS_pool_get_subpool(block_ptr % structs, 'atmos_coupling', atmosCoupling)
         call MPAS_pool_get_subpool(block_ptr % structs, 'ocean_coupling', oceanCoupling)
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

         call MPAS_pool_get_array(oceanCoupling, 'frazilMassAdjust', frazilMassAdjust)

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

         if (config_use_column_biogeochemistry .and. config_couple_biogeochemistry_fields) then
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
         if (config_use_zaerosols .and. config_couple_biogeochemistry_fields) &
            call mpas_pool_get_array(biogeochemistry, 'oceanDustIronFlux', oceanDustIronFlux)

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
            i2x_im(n, index_i2x_Si_ifrac)    = ailohi

            if (config_use_data_icebergs) then
               i2x_im(n, index_i2x_Fioi_bergw) = bergFreshwaterFlux(i)
               i2x_im(n, index_i2x_Fioi_bergh) = bergLatentHeatFlux(i)
            endif

            if ( ailohi > 0.0_RKIND ) then

               !-------states--------------------
               i2x_im(n, index_i2x_Si_t)  = Tsrf
               i2x_im(n, index_i2x_Si_bpress)  = basalPressure
               i2x_im(n, index_i2x_Si_u10)  = atmosReferenceSpeed10m(i)
               i2x_im(n, index_i2x_Si_tref)  = atmosReferenceTemperature2m(i)
               i2x_im(n, index_i2x_Si_qref)  = atmosReferenceHumidity2m(i)
               i2x_im(n, index_i2x_Si_snowh)  = snowVolumeCell(i) / ailohi

               !--- a/i fluxes computed by ice
               i2x_im(n, index_i2x_Faii_taux) = tauxa
               i2x_im(n, index_i2x_Faii_tauy) = tauya
               i2x_im(n, index_i2x_Faii_lat ) = latentHeatFlux(i)
               i2x_im(n, index_i2x_Faii_sen ) = sensibleHeatFlux(i)
               i2x_im(n, index_i2x_Faii_lwup) = longwaveUp(i)
               i2x_im(n, index_i2x_Faii_evap) = evaporativeWaterFlux(i)
               i2x_im(n, index_i2x_Faii_swnet) = absorbedShortwaveFlux(i)
               i2x_im(n, index_i2x_Faii_evap) = evaporativeWaterFlux(i)

               if (config_use_column_shortwave) then
                  i2x_im(n, index_i2x_Si_avsdr)  = albedoVisibleDirectCell(i)
                  i2x_im(n, index_i2x_Si_anidr)  = albedoIRDirectCell(i)
                  i2x_im(n, index_i2x_Si_avsdf)  = albedoVisibleDiffuseCell(i)
                  i2x_im(n, index_i2x_Si_anidf)  = albedoIRDiffuseCell(i)

                  i2x_im(n, index_i2x_Faii_swnet) = absorbedShortwaveFlux(i)
               endif

               ! i/o fluxes computed by ice, as well as additional freshwater and salt calculated at the last
               ! coupling import and needed to grow sea ice from frazil passed from the ocean model in the
               ! field frazilMassAdjust.
               i2x_im(n, index_i2x_Fioi_melth) = oceanHeatFlux(i)
               i2x_im(n, index_i2x_Fioi_swpen) = oceanShortwaveFlux(i)
               i2x_im(n, index_i2x_Fioi_meltw) = oceanFreshWaterFlux(i) + frazilMassAdjust(i)/ailohi
               i2x_im(n, index_i2x_Fioi_salt ) = oceanSaltFlux(i) + seaiceReferenceSalinity*0.001_RKIND*frazilMassAdjust(i)/ailohi
               i2x_im(n, index_i2x_Fioi_taux ) = tauxo
               i2x_im(n, index_i2x_Fioi_tauy ) = tauyo

               ! export biogeochemistry fields, if configured
               if (config_use_column_biogeochemistry .and. config_couple_biogeochemistry_fields) then
                  ! convert from mmol N/m^3 to mmol C/m^3
                  i2x_im(n, index_i2x_Fioi_algae1) = oceanAlgaeFlux(1,i) * carbonToNitrogenRatioAlgae(1)
                  i2x_im(n, index_i2x_Fioi_algae2) = oceanAlgaeFlux(2,i) * carbonToNitrogenRatioAlgae(2)
                  i2x_im(n, index_i2x_Fioi_algae3) = oceanAlgaeFlux(3,i) * carbonToNitrogenRatioAlgae(3)
                  i2x_im(n, index_i2x_Fioi_doc1  ) = oceanDOCFlux(1,i)
                  i2x_im(n, index_i2x_Fioi_doc2  ) = oceanDOCFlux(2,i)
                  i2x_im(n, index_i2x_Fioi_doc3  ) = oceanDOCFlux(3,i) !JW set to 0?
                  i2x_im(n, index_i2x_Fioi_dic1  ) = oceanDICFlux(1,i)
                  i2x_im(n, index_i2x_Fioi_don1  ) = oceanDONFlux(1,i)
                  i2x_im(n, index_i2x_Fioi_no3   ) = oceanNitrateFlux(i)
                  i2x_im(n, index_i2x_Fioi_sio3  ) = oceanSilicateFlux(i)
                  i2x_im(n, index_i2x_Fioi_nh4   ) = oceanAmmoniumFlux(i)
                  i2x_im(n, index_i2x_Fioi_dms   ) = oceanDMSFlux(i)
                  i2x_im(n, index_i2x_Fioi_dmspp ) = oceanDMSPpFlux(i)
                  i2x_im(n, index_i2x_Fioi_dmspd ) = oceanDMSPdFlux(i)
                  i2x_im(n, index_i2x_Fioi_docr  ) = oceanHumicsFlux(i)
                  ! convert from umol Fe/m^3 to mmol Fe/m^3
                  i2x_im(n, index_i2x_Fioi_fep1  ) = oceanParticulateIronFlux(1,i) / 1000._RKIND
                  i2x_im(n, index_i2x_Fioi_fep2  ) = oceanParticulateIronFlux(2,i) / 1000._RKIND
                  i2x_im(n, index_i2x_Fioi_fed1  ) = oceanDissolvedIronFlux(1,i) / 1000._RKIND
                  i2x_im(n, index_i2x_Fioi_fed2  ) = oceanDissolvedIronFlux(2,i) / 1000._RKIND
               endif
               ! export dust, kg/m2/s
               if (config_use_zaerosols .and. config_couple_biogeochemistry_fields) &
                  i2x_im(n, index_i2x_Fioi_dust1 ) = oceanDustIronFlux(i)
            endif
         enddo

         block_ptr => block_ptr % next
      enddo

      ent_type = 1 ! cells
      !  set all tags in one method
      tagname = trim(seq_flds_i2x_fields)//C_NULL_CHAR
      ierr = iMOAB_SetDoubleTagStorage ( MPSIID, tagname, totalmbls , ent_type, i2x_im(1, 1) )
      if ( ierr /= 0 ) then
         write(iceLogUnit,*) 'Fail to set MOAB fields '
      endif

      call seq_timemgr_EClockGetData( EClock, stepno=cur_ice_stepno )
#ifdef MOABDEBUG
      write(lnum,"(I0.2)")cur_ice_stepno
      outfile = 'ice_export_'//trim(lnum)//'.h5m'//C_NULL_CHAR
      wopts   = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
      ierr = iMOAB_WriteMesh(MPSIID, outfile, wopts)
#endif
   end subroutine ice_export_moab
!***********************************************************************
!BOP
! !IROUTINE: ice_import_moab
! !INTERFACE:

 subroutine ice_import_moab(Eclock)!{{{

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
!    o  q        -- ocn frazil heat flux(+) / melt potential(-)
!    o  frazil   -- ocn frazil mass flux
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
     use iMOAB,  only       : iMOAB_GetDoubleTagStorage, iMOAB_WriteMesh

    type(ESMF_Clock)            , intent(inout) :: EClock

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
      frazilMassFlux, &
      frazilMassFluxRev

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
      frazilMassAdjustField,            &
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
      atmosDustFluxField,               &
      atmosWetDustFluxField,            &
      atmosDryDustFluxField

   real (kind=RKIND), dimension(:), pointer :: &
      seaSurfaceTemperature,       &
      seaSurfaceSalinity,          &
      seaFreezingTemperature,      &
      freezingMeltingPotential,    &
      frazilMassAdjust,            &
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
      atmosDustFlux,               &
      atmosWetDustFlux,            &
      atmosDryDustFlux
   character(CXX) :: tagname
   integer        :: ierr, ent_type
   integer :: cur_ice_stepno
#ifdef MOABDEBUG
    character*100 outfile, wopts, lnum
#endif

    call seq_timemgr_EClockGetData( EClock, stepno=cur_ice_stepno )
#ifdef MOABDEBUG
    write(lnum,"(I0.2)")cur_ice_stepno
    outfile = 'ice_import_'//trim(lnum)//'.h5m'//C_NULL_CHAR
    wopts   = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
    ierr = iMOAB_WriteMesh(MPSIID, outfile, wopts)
    if (ierr > 0 )  then
       write(iceLogUnit,*) 'Fail to write ice state '
    endif
#endif
!-----------------------------------------------------------------------
!
!  zero out padded cells
!
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------

   n = 0
   x2i_im(: ,:) = 0.0_RKIND

   ent_type = 1 ! cells
  !  set all tags in one method
   tagname = trim(seq_flds_x2i_fields)//C_NULL_CHAR
   ierr = iMOAB_GetDoubleTagStorage ( MPSIID, tagname, totalmblr , ent_type, x2i_im(1, 1) )
   if ( ierr /= 0 ) then
      write(iceLogUnit,*) 'Fail to get seq_flds_x2i_fields '
   endif
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
      call mpas_pool_get_array(oceanCoupling, 'frazilMassAdjust', frazilMassAdjust)
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
         call mpas_pool_get_array(biogeochemistry, 'carbonToNitrogenRatioAlgae', carbonToNitrogenRatioAlgae)
         call mpas_pool_get_array(biogeochemistry, 'carbonToNitrogenRatioDON', carbonToNitrogenRatioDON)
      endif

      call mpas_pool_get_config(configs, "config_use_zaerosols", config_use_zaerosols)
         if (config_use_zaerosols) then
            call mpas_pool_get_array(biogeochemistry, "atmosBlackCarbonFlux", atmosBlackCarbonFlux)
            call mpas_pool_get_array(biogeochemistry, "atmosDustFlux", atmosDustFlux)
            call mpas_pool_get_array(biogeochemistry, "atmosWetDustFlux", atmosWetDustFlux)
            call mpas_pool_get_array(biogeochemistry, "atmosDryDustFlux", atmosDryDustFlux)
         endif
      if (config_use_column_biogeochemistry .or. config_use_zaerosols) &
         call mpas_pool_get_array(biogeochemistry, 'oceanZAerosolConc', oceanZAerosolConc)

      do i = 1, nCellsSolve
        n = n + 1
        seaSurfaceTemperature(i)       = x2i_im(n,index_x2i_So_t)
        seaSurfaceSalinity(i)          = x2i_im(n,index_x2i_So_s)

        seaFreezingTemperature(i) = colpkg_sea_freezing_temperature(seaSurfaceSalinity(i))

        uOceanVelocity(i)              = x2i_im(n,index_x2i_So_u)
        vOceanVelocity(i)              = x2i_im(n,index_x2i_So_v)
        seaSurfaceTiltU(i)             = x2i_im(n,index_x2i_So_dhdx)
        seaSurfaceTiltV(i)             = x2i_im(n,index_x2i_So_dhdy)

        if (trim(config_ocean_surface_type) == "free") then ! free surface (MPAS-O)

           ! freezingMeltingPotential(i) is the ocean energy associated with frazil formation
           ! when it is positive and frazilMassFlux is positive. Conversely, freezingMeltingPotential(i)
           ! is negative when there is the melting potential in which case frazilMassFlux is zero.

           freezingMeltingPotential(i) = x2i_im(n,index_x2i_Fioo_q)

           frazilMassFlux              = x2i_im(n,index_x2i_Fioo_frazil)

           ! Now determine the sea ice mass associated with the frazil heat flux given when
           ! freezingMeltingPotential(i) is positive. This produces a revised mass flux, given
           ! in frazilMassFluxRev for the given sea surface salinity. The resulting difference
           ! is assigned to frazilMassAdjust(i) which is exported to the ocean in the subsequent
           ! coupling step as a freshwater and salt flux. This step is required to balance mass
           ! and heat with the ocean.

           call frazil_mass(freezingMeltingPotential(i), frazilMassFluxRev, seaSurfaceSalinity(i))

           frazilMassAdjust(i) = frazilMassFlux-frazilMassFluxRev

        else ! non-free surface (SOM)

           freezingMeltingPotential(i) = x2i_im(n,index_x2i_Fioo_q)

        endif

        airLevelHeight(i)              = x2i_im(n,index_x2i_Sa_z)
        airPotentialTemperature(i)     = x2i_im(n,index_x2i_Sa_ptem)
        airTemperature(i)              = x2i_im(n,index_x2i_Sa_tbot)
        airSpecificHumidity(i)         = x2i_im(n,index_x2i_Sa_shum)
        airDensity(i)                  = x2i_im(n,index_x2i_Sa_dens)
        shortwaveVisibleDirectDown(i)  = x2i_im(n,index_x2i_Faxa_swvdr)
        shortwaveVisibleDiffuseDown(i) = x2i_im(n,index_x2i_Faxa_swvdf)
        shortwaveIRDirectDown(i)       = x2i_im(n,index_x2i_Faxa_swndr)
        shortwaveIRDiffuseDown(i)      = x2i_im(n,index_x2i_Faxa_swndf)
        longwaveDown(i)                = x2i_im(n,index_x2i_Faxa_lwdn)
        rainfallRate(i)                = x2i_im(n,index_x2i_Faxa_rain)
        snowfallRate(i)                = x2i_im(n,index_x2i_Faxa_snow)
        uAirVelocity(i)                = x2i_im(n,index_x2i_Sa_u)
        vAirVelocity(i)                = x2i_im(n,index_x2i_Sa_v)

        ! set aerosols, if configured
        if (config_use_aerosols) then
           if (config_use_modal_aerosols) then
              atmosAerosolFlux(1,i) = x2i_im(n,index_x2i_Faxa_bcphodry) &
                                    + x2i_im(n,index_x2i_Faxa_bcphidry)
              atmosAerosolFlux(2,i) = x2i_im(n,index_x2i_Faxa_bcphiwet)
              ! combine all the dust into one category
              atmosAerosolFlux(3,i) = x2i_im(n,index_x2i_Faxa_dstwet1) &
                                    + x2i_im(n,index_x2i_Faxa_dstwet2) &
                                    + x2i_im(n,index_x2i_Faxa_dstwet3) &
                                    + x2i_im(n,index_x2i_Faxa_dstwet4) &
                                    + x2i_im(n,index_x2i_Faxa_dstdry1) &
                                    + x2i_im(n,index_x2i_Faxa_dstdry2) &
                                    + x2i_im(n,index_x2i_Faxa_dstdry3) &
                                    + x2i_im(n,index_x2i_Faxa_dstdry4)
           else
              atmosAerosolFlux(1,i) = x2i_im(n,index_x2i_Faxa_bcphodry)
              atmosAerosolFlux(2,i) = x2i_im(n,index_x2i_Faxa_bcphidry) &
                                    + x2i_im(n,index_x2i_Faxa_bcphiwet)
              ! combine all the dust into one category
              atmosAerosolFlux(3,i) = x2i_im(n,index_x2i_Faxa_dstwet1) &
                                    + x2i_im(n,index_x2i_Faxa_dstwet2) &
                                    + x2i_im(n,index_x2i_Faxa_dstwet3) &
                                    + x2i_im(n,index_x2i_Faxa_dstwet4) &
                                    + x2i_im(n,index_x2i_Faxa_dstdry1) &
                                    + x2i_im(n,index_x2i_Faxa_dstdry2) &
                                    + x2i_im(n,index_x2i_Faxa_dstdry3) &
                                    + x2i_im(n,index_x2i_Faxa_dstdry4)
           endif
        endif

        ! import biogeochemistry fields, if configured
        if (config_use_column_biogeochemistry) then
           oceanAlgaeConc(1,i)           = x2i_im(n,index_x2i_So_algae1)
           oceanAlgaeConc(2,i)           = x2i_im(n,index_x2i_So_algae2)
           oceanAlgaeConc(3,i)           = x2i_im(n,index_x2i_So_algae3)
           oceanDOCConc(1,i)             = x2i_im(n,index_x2i_So_doc1)
           oceanDOCConc(2,i)             = x2i_im(n,index_x2i_So_doc2)
           oceanDOCConc(3,i)             = 0.0_RKIND
           oceanDICConc(1,i)             = x2i_im(n,index_x2i_So_dic1) !JW not used, set to 0?
           oceanDONConc(1,i)             = x2i_im(n,index_x2i_So_don1)
           oceanNitrateConc(i)           = x2i_im(n,index_x2i_So_no3)
           oceanSilicateConc(i)          = x2i_im(n,index_x2i_So_sio3)
           oceanAmmoniumConc(i)          = x2i_im(n,index_x2i_So_nh4)
           oceanDMSConc(i)               = x2i_im(n,index_x2i_So_dms)
           oceanDMSPConc(i)              = x2i_im(n,index_x2i_So_dmsp)
           oceanHumicsConc(i)            = x2i_im(n,index_x2i_So_docr)
           oceanParticulateIronConc(1,i) = x2i_im(n,index_x2i_So_fep1)
           oceanParticulateIronConc(2,i) = x2i_im(n,index_x2i_So_fep2)
           oceanDissolvedIronConc(1,i)   = x2i_im(n,index_x2i_So_fed1)
           oceanDissolvedIronConc(2,i)   = x2i_im(n,index_x2i_So_fed2)
        endif
        if (config_use_zaerosols .or. config_use_column_biogeochemistry) then
           oceanZAerosolConc(1,i)        = 0.0_RKIND
           oceanZAerosolConc(2,i)        = 0.0_RKIND
           oceanZAerosolConc(3,i)        = 0.0_RKIND
           oceanZAerosolConc(4,i)        = 0.0_RKIND
           oceanZAerosolConc(5,i)        = 0.0_RKIND
           oceanZAerosolConc(6,i)        = 0.0_RKIND
        end if
        ! set zaerosols, if configured
        if (config_use_zaerosols) then
           if (config_use_modal_aerosols) then
              atmosBlackCarbonFlux(1,i) = x2i_im(n,index_x2i_Faxa_bcphodry) &
                                        + x2i_im(n,index_x2i_Faxa_bcphidry)
              atmosBlackCarbonFlux(2,i) = x2i_im(n,index_x2i_Faxa_bcphiwet)
              ! combine wet and dry dust
              atmosDustFlux(1,i) = x2i_im(n,index_x2i_Faxa_dstwet1) &
                                 + x2i_im(n,index_x2i_Faxa_dstdry1)
              atmosDustFlux(2,i) = x2i_im(n,index_x2i_Faxa_dstwet2) &
                                 + x2i_im(n,index_x2i_Faxa_dstdry2)
              atmosDustFlux(3,i) = x2i_im(n,index_x2i_Faxa_dstwet3) &
                                 + x2i_im(n,index_x2i_Faxa_dstdry3)
              atmosDustFlux(4,i) = x2i_im(n,index_x2i_Faxa_dstwet4) &
                                 + x2i_im(n,index_x2i_Faxa_dstdry4)

              ! wet dust
              atmosWetDustFlux(1,i) = x2i_im(n,index_x2i_Faxa_dstwet1)
              atmosWetDustFlux(2,i) = x2i_im(n,index_x2i_Faxa_dstwet2)
              atmosWetDustFlux(3,i) = x2i_im(n,index_x2i_Faxa_dstwet3)
              atmosWetDustFlux(4,i) = x2i_im(n,index_x2i_Faxa_dstwet4)

              ! dry dust
              atmosDryDustFlux(1,i) = x2i_im(n,index_x2i_Faxa_dstdry1)
              atmosDryDustFlux(2,i) = x2i_im(n,index_x2i_Faxa_dstdry2)
              atmosDryDustFlux(3,i) = x2i_im(n,index_x2i_Faxa_dstdry3)
              atmosDryDustFlux(4,i) = x2i_im(n,index_x2i_Faxa_dstdry4)
           else
              atmosBlackCarbonFlux(1,i) = x2i_im(n,index_x2i_Faxa_bcphodry)
              atmosBlackCarbonFlux(2,i) = x2i_im(n,index_x2i_Faxa_bcphidry) &
                                        + x2i_im(n,index_x2i_Faxa_bcphiwet)
              ! combine wet and dry dust
              atmosDustFlux(1,i) = x2i_im(n,index_x2i_Faxa_dstwet1) &
                                 + x2i_im(n,index_x2i_Faxa_dstdry1)
              atmosDustFlux(2,i) = x2i_im(n,index_x2i_Faxa_dstwet2) &
                                 + x2i_im(n,index_x2i_Faxa_dstdry2)
              atmosDustFlux(3,i) = x2i_im(n,index_x2i_Faxa_dstwet3) &
                                 + x2i_im(n,index_x2i_Faxa_dstdry3)
              atmosDustFlux(4,i) = x2i_im(n,index_x2i_Faxa_dstwet4) &
                                 + x2i_im(n,index_x2i_Faxa_dstdry4)

              ! wet dust
              atmosWetDustFlux(1,i) = x2i_im(n,index_x2i_Faxa_dstwet1)
              atmosWetDustFlux(2,i) = x2i_im(n,index_x2i_Faxa_dstwet2)
              atmosWetDustFlux(3,i) = x2i_im(n,index_x2i_Faxa_dstwet3)
              atmosWetDustFlux(4,i) = x2i_im(n,index_x2i_Faxa_dstwet4)

              ! dry dust
              atmosDryDustFlux(1,i) = x2i_im(n,index_x2i_Faxa_dstdry1)
              atmosDryDustFlux(2,i) = x2i_im(n,index_x2i_Faxa_dstdry2)
              atmosDryDustFlux(3,i) = x2i_im(n,index_x2i_Faxa_dstdry3)
              atmosDryDustFlux(4,i) = x2i_im(n,index_x2i_Faxa_dstdry4)
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
   call mpas_pool_get_field(oceanCoupling, 'frazilMassAdjust', frazilMassAdjustField)
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
   endif
   if (config_use_zaerosols .or. config_use_column_biogeochemistry) &
      call mpas_pool_get_field(biogeochemistry, 'oceanZAerosolConc', oceanZAerosolConcField)
   if (config_use_zaerosols) then
      call mpas_pool_get_field(biogeochemistry, "atmosBlackCarbonFlux", atmosBlackCarbonFluxField)
      call mpas_pool_get_field(biogeochemistry, "atmosDustFlux", atmosDustFluxField)
      call mpas_pool_get_field(biogeochemistry, "atmosWetDustFlux", atmosWetDustFluxField)
      call mpas_pool_get_field(biogeochemistry, "atmosDryDustFlux", atmosDryDustFluxField)
   endif
   call mpas_dmpar_exch_halo_field(seaSurfaceTemperatureField)
   call mpas_dmpar_exch_halo_field(seaSurfaceSalinityField)
   call mpas_dmpar_exch_halo_field(seaFreezingTemperatureField)
   call mpas_dmpar_exch_halo_field(freezingMeltingPotentialField)
   call mpas_dmpar_exch_halo_field(frazilMassAdjustField)
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
   endif
   if (config_use_zaerosols .or. config_use_column_biogeochemistry) &
      call mpas_dmpar_exch_halo_field(oceanZAerosolConcField)
   if (config_use_zaerosols) then
         call mpas_dmpar_exch_halo_field(atmosBlackCarbonFluxField)
         call mpas_dmpar_exch_halo_field(atmosDustFluxField)
         call mpas_dmpar_exch_halo_field(atmosWetDustFluxField)
         call mpas_dmpar_exch_halo_field(atmosDryDustFluxField)
   endif

! REVISION HISTORY:
!-----------------------------------------------------------------------
!EOC

 end subroutine ice_import_moab!}}}
#endif
end module ice_comp_mct

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
