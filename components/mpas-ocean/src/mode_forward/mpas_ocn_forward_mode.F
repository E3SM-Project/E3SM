! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_forward_mode
!
!> \brief Main driver for MPAS ocean core
!> \author Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date   September 2011
!> \details
!>  This module contains initialization and timestep drivers for
!>  the MPAS ocean core.
!
!-----------------------------------------------------------------------

module ocn_forward_mode

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_stream_manager
   use mpas_timekeeping
   use mpas_dmpar
   use mpas_timer
   use mpas_log
   use mpas_decomp
   use mpas_log

   use ocn_analysis_driver
   use ocn_init_routines

   use ocn_time_integration
   use ocn_time_integration_split
   use ocn_time_integration_si
   use ocn_tendency
   use ocn_diagnostics
   use ocn_diagnostics_variables
   use ocn_test
   use ocn_mesh

   use ocn_eddy_parameterization_helpers
   use ocn_thick_hadv
   use ocn_thick_vadv
   use ocn_thick_ale
   use ocn_thick_surface_flux

   use ocn_vel_pressure_grad
   use ocn_vel_vadv
   use ocn_vel_hmix
   use ocn_vel_forcing
   use ocn_vel_hadv_coriolis
   use ocn_vel_forcing_surface_stress
   use ocn_surface_bulk_forcing
   use ocn_surface_land_ice_fluxes
   use ocn_frazil_forcing
   use ocn_tidal_forcing
   use ocn_vel_tidal_potential
   use ocn_vel_self_attraction_loading

   use ocn_tracer_hmix
   use ocn_tracer_hmix_redi
   use ocn_tracer_surface_flux_to_tend
   use ocn_tracer_short_wave_absorption
   use ocn_tracer_short_wave_absorption_variable
   use ocn_tracer_nonlocalflux
   use ocn_tracer_advection
   use ocn_tracer_ecosys
   use ocn_tracer_DMS
   use ocn_tracer_MacroMolecules
   use ocn_tracer_ideal_age
   use ocn_tracer_CFC
   use ocn_tracer_surface_restoring
   use ocn_gm
   use ocn_submesoscale_eddies
   use ocn_stokes_drift 

   use ocn_high_freq_thickness_hmix_del2

   use ocn_equation_of_state

   use ocn_vmix
   use ocn_vmix_gotm, only: ocn_vmix_gotm_finalize

   use ocn_vertical_advection
   use ocn_vertical_remap, only: ocn_vertical_remap_init

   use ocn_forcing
   use ocn_time_varying_forcing

   use ocn_constants
   use ocn_config

   implicit none
   private

   public :: ocn_forward_mode_init, &
             ocn_forward_mode_run, &
             ocn_forward_mode_finalize, &
             ocn_forward_mode_setup_clock

!***********************************************************************

   contains

!***********************************************************************
!
!  function ocn_forward_mode_init
!
!> \brief   Initialize MPAS-Ocean core
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details
!>  This function calls all initializations required to begin a
!>  simulation with MPAS-Ocean
!
!-----------------------------------------------------------------------

   function ocn_forward_mode_init(domain, startTimeStamp) result(ierr)!{{{

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: &
         domain  !< [inout] structure containing most ocean variables

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer :: ierr    !< [out] returned error code

      character(len=*), intent(out) :: &
         startTimeStamp  !< [out] start time in time stamp string

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: err_tmp    ! local error codes from subroutine calls

      real (kind=RKIND) :: &
         dt,               &! time step in seconds
         maxDensity,       &! local max of mesh density
         maxDensity_global  ! global max of mesh density

      real (kind=RKIND), dimension(:), pointer :: &
         meshDensity        ! pointer to mesh density array

      type (block_type), pointer :: &
         block              ! structure with all data for a given block

      type (mpas_pool_type), pointer :: &
         diagnosticsPool    ! structure containing diagnostics data

      type (MPAS_Time_type) :: &
         xtime_timeType,               &! current time
         simulationStartTime_timeType, &! start time whole simulation
         referenceTime_timeType,       &! reference time for Time coordinate
         startTime,                    &! start time this run
         alarmTime                      ! alarm for salinity restoring

      type (field0DReal), pointer :: TimeField    ! to add units

      type (MPAS_TimeInterval_type) :: &
         timeStep,      &! time step in time manager form
         alarmTimeStep   ! alarm interval in time manager form

      character (len=StrKIND) :: units

      ! End preamble
      !-------------
      ! Begin code

      !-----------------------------------------------------------------
      ! Initialize variables needed for all
      !-----------------------------------------------------------------
      ierr = 0

      ! Setup ocean config pool
      call ocn_config_init(domain % configs)

      ! Set up ocean constants
      call ocn_constants_init(domain % configs, domain % packages)

      ! Set startTimeStamp based on start time of simulation clock
      startTime = mpas_get_clock_time(domain%clock, MPAS_START_TIME, &
                                      err_tmp)
      if (err_tmp /= 0) then
         call mpas_log_write('Error getting start time in init', &
              MPAS_LOG_ERR, masterOnly=.true., flushNow=.true.)
         ierr = ior(ierr, err_tmp)
      endif
      call mpas_get_time(startTime, dateTimeString=startTimeStamp)

      ! Extract block structure from domain
      ! Currently only support one block per task
      block => domain % blocklist

      !-----------------------------------------------------------------
      ! Initialize data from external file input
      !-----------------------------------------------------------------

      ! Read mesh data for model
      call mpas_timer_start('io_read')
      call MPAS_stream_mgr_read(domain%streamManager, streamID='mesh', &
                             whence=MPAS_STREAM_NEAREST, ierr=err_tmp)
      if (err_tmp /= 0) then
         call mpas_log_write('Error reading mesh in init', &
              MPAS_LOG_ERR, masterOnly=.true., flushNow=.true.)
         ierr = ior(ierr, err_tmp)
      endif

      ! Read initial state from either restart or init condition file
      if ( config_do_restart ) then
         if ( config_read_nearest_restart ) then
            call MPAS_stream_mgr_read(domain%streamManager, &
                                      streamID='restart', &
                                      whence=MPAS_STREAM_NEAREST, &
                                      ierr=err_tmp)
         else
            call MPAS_stream_mgr_read(domain%streamManager, &
                                      streamID='restart', ierr=err_tmp)
         end if
      else
         call MPAS_stream_mgr_read(domain % streamManager, &
                                   streamID='input', ierr=err_tmp)
      end if
      if (err_tmp /= 0) then
         call mpas_log_write('Error reading initial state in init', &
              MPAS_LOG_ERR, masterOnly=.true., flushNow=.true.)
         ierr = ior(ierr, err_tmp)
      endif

      ! Read any analysis initialization
      call ocn_analysis_bootstrap(domain, err=err_tmp)

      call mpas_timer_stop('io_read')
      if (ierr /= 0) return

      call mpas_timer_start('reset_io_alarms')
      call mpas_stream_mgr_reset_alarms(domain%streamManager, &
                                        streamID='input', ierr=err_tmp)
      call mpas_stream_mgr_reset_alarms(domain%streamManager, &
                                      streamID='restart', ierr=err_tmp)
      call mpas_stream_mgr_reset_alarms(domain%streamManager, &
                            direction=MPAS_STREAM_OUTPUT, ierr=err_tmp)
      if (err_tmp /= 0) then
         call mpas_log_write('Error resetting io alarms in init', &
              MPAS_LOG_ERR, masterOnly=.true., flushNow=.true.)
         ierr = ior(ierr, err_tmp)
      endif
      call mpas_timer_stop('reset_io_alarms')

      ! Read the remaining input streams
      call mpas_timer_start('io_read')
      call mpas_stream_mgr_read(domain%streamManager, ierr=err_tmp)
      ierr = ior(ierr, err_tmp)
      call mpas_timer_stop('io_read')
      if (ierr /= 0) &
         call mpas_log_write('Error reading remaining io streams in init', &
              MPAS_LOG_ERR, masterOnly=.true., flushNow=.true.)

      call mpas_timer_start('reset_io_alarms')
      call mpas_stream_mgr_reset_alarms(domain%streamManager, &
                       direction=MPAS_STREAM_INPUT, ierr=err_tmp)
      ierr = ior(ierr, err_tmp)
      call mpas_timer_stop('reset_io_alarms')
      if (ierr /= 0) then
         call mpas_log_write('Error resetting io alarms(2) in init', &
              MPAS_LOG_ERR, masterOnly=.true., flushNow=.true.)
         return
      endif

      !-----------------------------------------------------------------
      ! Initialize mesh before anything else since many depend on it
      !-----------------------------------------------------------------

      ! fill the ocean mesh structure
      call ocn_meshCreate(domain)

      !-----------------------------------------------------------------
      ! Initialize submodules and fields that are required by others
      ! or that have dependencies that require specific ordering.
      !-----------------------------------------------------------------

      ! Needed for later call to compute land ice fluxes
      call ocn_surface_land_ice_fluxes_init(err_tmp)
      ierr = ior(ierr, err_tmp)

      ! Needed for later diagnostic_solve computations
      call ocn_equation_of_state_init(domain, err_tmp)
      ierr = ior(ierr, err_tmp)
      call ocn_diagnostics_init(domain, err_tmp)
      ierr = ior(ierr,err_tmp)

      if(ierr.eq.1) then
          call mpas_log_write('An error was encountered while initializing the MPAS-Ocean forward mode', MPAS_LOG_CRIT)
      endif

      !
      ! Initialize core
      !

      timeStep = mpas_get_clock_timestep(domain % clock, ierr=err_tmp)
      call mpas_get_timeInterval(timeStep, dt=dt)

      call ocn_init_routines_block(block, dt, ierr)

      if(ierr.eq.1) then
          call mpas_log_write('An error was encountered in ocn_init_routines_block', MPAS_LOG_CRIT)
      endif

      xtime = startTimeStamp

      ! Set simulationStartTime only if that variable is not read from the restart file.
      if (trim(simulationStartTime)=="no_date_available") then
         simulationStartTime = startTimeStamp
      end if

      ! compute time since start of simulation, in days
      call mpas_set_time(xtime_timeType, dateTimeString=xtime)
      call mpas_set_time(simulationStartTime_timeType, dateTimeString=simulationStartTime)
      call mpas_get_timeInterval(xtime_timeType - simulationStartTime_timeType, dt=daysSinceStartOfSim)
      daysSinceStartOfSim = daysSinceStartOfSim*days_per_second

      call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

      ! add units of "days since <ref_date>" to Time
      write(units, '(a10,a20)') 'days since', &
           config_output_reference_time(1:10)//" "//config_output_reference_time(12:19)
      call mpas_pool_get_Field(diagnosticsPool, 'Time', TimeField)
      call mpas_add_att(TimeField % attLists(1) % attList, 'units', units)

      ! Add calendar attribute to Time
      call mpas_add_att(TimeField % attLists(1) % attList, 'calendar', config_calendar_type)

      ! compute Time as days since config_output_reference_time
      call mpas_set_time(referenceTime_timeType, dateTimeString=config_output_reference_time)
      call mpas_get_timeInterval(xtime_timeType - referenceTime_timeType, dt=Time)
      Time = Time*days_per_second
      ! bounds on instantaneous time are just the time itself
      Time_bnds(1) = Time
      Time_bnds(2) = Time

      ! Update any mesh fields that may have been modified during init
      call ocn_meshUpdateFields(domain)

      if (config_conduct_tests) then
         call mpas_timer_start("test suite")
         call ocn_test_suite(domain,ierr)
         call mpas_timer_stop("test suite")
      endif

      call ocn_analysis_init(domain, err_tmp)
      ierr = ior(ierr, err_tmp)
      if(ierr.eq.1) then
          call mpas_log_write('An error was encountered while initializing ' &
               // 'the analysis members in the MPAS-Ocean forward mode', MPAS_LOG_CRIT)
      endif

      if (config_use_surface_salinity_monthly_restoring) then
      !initialize the alarm for reading salinity data for restoring
      if ( config_surface_salinity_monthly_restoring_compute_interval == 'dt') then
         alarmTimeStep = mpas_get_clock_timestep(domain % clock, err_tmp)
         call mpas_get_timeInterval(alarmTimeStep, timeString= &
                      config_surface_salinity_monthly_restoring_compute_interval, ierr=err_tmp)
      endif
      alarmTime = mpas_get_clock_time(domain % clock, MPAS_START_TIME, ierr=err_tmp)

      call mpas_set_timeInterval(alarmTimeStep, timeString= &
              config_surface_salinity_monthly_restoring_compute_interval, ierr=err_tmp)

      call mpas_add_clock_alarm(domain % clock, 'salinityDataReadAlarm', alarmTime, &
                alarmTimeInterval=alarmTimeStep, ierr=err_tmp)
      call mpas_reset_clock_alarm(domain % clock, 'salinityDataReadAlarm', ierr=err_tmp)

      endif

      ! Initialize self-attraction and loading
      call ocn_vel_self_attraction_loading_init(domain,err_tmp)
      ierr = ior(ierr, err_tmp)

      !-----------------------------------------------------------------
      ! Complete all remaining module initializations.
      !-----------------------------------------------------------------

      ! add some mesh metadata to tau performance tools if enabled
      call ocn_init_metadata(domain)

      call ocn_timestep_init(domain, dt, ierr)
      if (ierr /= 0) then
         call mpas_log_write( &
              'Error initializing ocean time integration ', &
              MPAS_LOG_ERR, masterOnly=.true., flushNow=.true.)
         return
      endif

      ! Initialize thickness tendency modules
      call ocn_thick_hadv_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      call ocn_thick_vadv_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      call ocn_thick_surface_flux_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      call ocn_thick_ale_init(err_tmp)
      ierr = ior(ierr,err_tmp)
      if (ierr /= 0) then
         call mpas_log_write( &
              'Error initializing thick tendency modules', &
              MPAS_LOG_ERR, masterOnly=.true., flushNow=.true.)
         return
      endif

      ! Initialize velocity tendency modules
      call ocn_vel_hadv_coriolis_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      call ocn_vel_hmix_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      call ocn_vel_pressure_grad_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      call ocn_vel_vadv_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      call ocn_vel_forcing_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      call ocn_vel_forcing_surface_stress_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      if (ierr /= 0) then
         call mpas_log_write( &
              'Error initializing vel tendency modules', &
              MPAS_LOG_ERR, masterOnly=.true., flushNow=.true.)
         return
      endif

      ! Initialize forcing modules
      call ocn_surface_bulk_forcing_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      call ocn_frazil_forcing_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      call ocn_tidal_forcing_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      call ocn_vel_tidal_potential_init(domain,err_tmp)
      ierr = ior(ierr, err_tmp)
      call ocn_forcing_init(err_tmp)
      ierr = ior(ierr,err_tmp)
      if (ierr /= 0) then
         call mpas_log_write( &
              'Error initializing forcing modules', &
              MPAS_LOG_ERR, masterOnly=.true., flushNow=.true.)
         return
      endif

      ! Initialize tracer tendencies
      call ocn_tracer_hmix_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      call ocn_tracer_hmix_redi_init(domain,err_tmp)
      ierr = ior(ierr, err_tmp)
      call ocn_tracer_surface_flux_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      call ocn_tracer_advection_init(err_tmp)
      ierr = ior(ierr,err_tmp)
      call ocn_tracer_short_wave_absorption_init(domain,err_tmp)
      ierr = ior(ierr,err_tmp)
      call ocn_gm_init(domain, err_tmp)
      ierr = ior(ierr,err_tmp)
      call ocn_submesoscale_init(err_tmp)
      ierr = ior(ierr,err_tmp)
      call ocn_tracer_nonlocalflux_init(err_tmp)
      ierr = ior(ierr,err_tmp)
      call ocn_tracer_ecosys_init(domain, err_tmp)
      ierr = ior(ierr,err_tmp)
      call ocn_tracer_DMS_init(domain, err_tmp)
      ierr = ior(ierr,err_tmp)
      call ocn_tracer_MacroMolecules_init(domain, err_tmp)
      ierr = ior(ierr,err_tmp)
      call ocn_tracer_ideal_age_init(domain, err_tmp)
      ierr = ior(ierr,err_tmp)
      call ocn_tracer_CFC_init(domain, err_tmp)
      ierr = ior(ierr,err_tmp)
      if (ierr /= 0) then
         call mpas_log_write( &
              'Error initializing tracer tendency modules', &
              MPAS_LOG_ERR, masterOnly=.true., flushNow=.true.)
         return
      endif

      ! Initialize other modules
      call ocn_vmix_init(domain, err_tmp)
      ierr = ior(ierr, err_tmp)
      call ocn_tendency_init(err_tmp)
      ierr = ior(ierr,err_tmp)
      call ocn_vertical_advection_init(err_tmp)
      if ( configVertAdvMethod == vertAdvRemap ) &
         call ocn_vertical_remap_init(err_tmp)
      call ocn_stokes_drift_init(err_tmp)
      ierr = ior(ierr,err_tmp)

      if (ierr /= 0) then
         call mpas_log_write( &
              'Error initializing vmix,EOS,tend modules', &
              MPAS_LOG_ERR, masterOnly=.true., flushNow=.true.)
         return
      endif

   !--------------------------------------------------------------------

   end function ocn_forward_mode_init!}}}

!***********************************************************************
!
!  function ocn_forward_mode_setup_clock
!
!> \brief   Initialize timer variables
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details
!>  This routine initializes all timer variables
!
!-----------------------------------------------------------------------

   function ocn_forward_mode_setup_clock(core_clock, configs) result(ierr)!{{{

      implicit none

      type (MPAS_Clock_type), intent(inout) :: core_clock
      type (mpas_pool_type), intent(inout) :: configs
      integer :: ierr

      type (MPAS_Time_Type) :: startTime, stopTime, alarmStartTime
      type (MPAS_TimeInterval_type) :: runDuration, timeStep, alarmTimeStep
      character(len=StrKIND) :: restartTimeStamp
      character(len=StrKIND), pointer :: config_start_time, config_stop_time, config_run_duration
      character(len=StrKIND), pointer :: config_dt, config_restart_timestamp_name
      integer :: err_tmp

      ierr = 0

      call mpas_pool_get_config(configs, 'config_dt', config_dt)
      call mpas_pool_get_config(configs, 'config_start_time', config_start_time)
      call mpas_pool_get_config(configs, 'config_stop_time', config_stop_time)
      call mpas_pool_get_config(configs, 'config_run_duration', config_run_duration)
      call mpas_pool_get_config(configs, 'config_restart_timestamp_name', config_restart_timestamp_name)

      if ( trim(config_start_time) == "file" ) then
         open(22,file=config_restart_timestamp_name,form='formatted',status='old')
         read(22,*) restartTimeStamp
         close(22)
         call mpas_set_time(curr_time=startTime, dateTimeString=restartTimeStamp, ierr=ierr)
      else
         call mpas_set_time(curr_time=startTime, dateTimeString=config_start_time, ierr=err_tmp)
      end if

      call mpas_set_timeInterval(timeStep, timeString=config_dt, ierr=err_tmp)
      if (trim(config_run_duration) /= "none") then
         call mpas_set_timeInterval(runDuration, timeString=config_run_duration, ierr=err_tmp)
         call mpas_create_clock(core_clock, startTime=startTime, timeStep=timeStep, runDuration=runDuration, ierr=err_tmp)

         if (trim(config_stop_time) /= "none") then
            call mpas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, ierr=err_tmp)
            if(startTime + runduration /= stopTime) then
               call mpas_log_write('Warning: config_run_duration and config_stop_time are inconsitent: using config_run_duration.')
            end if
         end if
      else if (trim(config_stop_time) /= "none") then
         call mpas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, ierr=err_tmp)
         call mpas_create_clock(core_clock, startTime=startTime, timeStep=timeStep, stopTime=stopTime, ierr=err_tmp)
      else
          call mpas_log_write('Error: Neither config_run_duration nor config_stop_time were specified.')
          ierr = 1
      end if

   end function ocn_forward_mode_setup_clock!}}}

!***********************************************************************
!
!  function ocn_forward_mode_run
!
!> \brief   Main driver for MPAS-Ocean time-stepping
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details
!>  This function includes the time-stepping loop, and calls
!>  routines to write output and restart files.
!
!-----------------------------------------------------------------------

   function ocn_forward_mode_run(domain) result(ierr)!{{{

      type (domain_type), intent(inout) :: domain

      integer :: itimestep, err
      real (kind=RKIND) :: dt
      type (block_type), pointer :: block_ptr

      type (MPAS_Time_Type) :: currTime
      character(len=StrKIND) :: timeStamp
      integer :: ierr, err_tmp

      type (mpas_pool_type), pointer :: averagePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: scratchPool

      type (MPAS_timeInterval_type) :: timeStep

      ierr = 0

      ! Eventually, dt should be domain specific
      timeStep = mpas_get_clock_timestep(domain % clock, ierr=ierr)
      call mpas_get_timeInterval(timeStep, dt=dt)

      currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
      call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
#ifdef MPAS_DEBUG
      call mpas_log_write('Initial time ' // trim(timeStamp))
#endif
      call ocn_analysis_compute_startup(domain, err)

      if (config_write_output_on_startup) then
          call mpas_timer_start('io_write')
          call mpas_stream_mgr_write(domain % streamManager, 'output', forceWriteNow=.true., ierr=ierr)
          call mpas_timer_stop('io_write')
      endif

      ! read initial data required for variable shortwave
      call mpas_timer_start('io_shortwave',.false.)
      call ocn_get_shortWaveData(domain % streamManager, domain, domain % clock, .true.)
      call mpas_timer_stop('io_shortwave')

      ! read initial data required for ecosys forcing
      if (config_use_ecosysTracers) then
        call mpas_timer_start('io_ecosys',.false.)
        call ocn_get_ecosysData(domain % streamManager, domain, domain % clock, .true.)
        call mpas_timer_stop('io_ecosys')
      endif

      ! read initial data required for CFC forcing
      if (config_use_CFCTracers) then
        call mpas_timer_start('io_CFC',.false.)
        call ocn_get_CFCData(domain % streamManager, domain, domain % clock, .true.)
        call mpas_timer_stop('io_CFC')
      endif

      ! read initial data required for monthly surface salinity restoring
      ! always execute this call as initial data needed regardless of alarm
      if (config_use_activeTracers_surface_restoring .and. config_use_surface_salinity_monthly_restoring) then
        call mpas_timer_start('io_monthly_surface_salinity',.false.)
        call ocn_get_surfaceSalinityData(domain % streamManager, domain, domain % clock, .true.)
        call mpas_timer_stop('io_monthly_surface_salinity')
      endif

      ! initialize time-varying forcing
      call ocn_time_varying_forcing_init(domain)
      call ocn_time_varying_forcing_get(domain % streamManager, domain, domain % clock)

      ! During integration, time level 1 stores the model state at the beginning of the
      !   time step, and time level 2 stores the state advanced dt in time by timestep(...)
      itimestep = 0

      do while (.not. mpas_is_clock_stop_time(domain % clock))
         call mpas_timer_start('io_read')
         call mpas_stream_mgr_read(domain % streamManager, ierr=ierr)
         call mpas_timer_stop('io_read')
         call mpas_timer_start('reset_io_alarms')
         call mpas_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_INPUT, ierr=ierr)

         ! Also restart all block_* streams, if any are defined.
         if ( mpas_stream_mgr_stream_exists(domain % streamManager, 'block_.*') ) then
            call mpas_stream_mgr_reset_alarms(domain % streamManager, streamID='block_.*', ierr=ierr)
         end if
         call mpas_timer_stop('reset_io_alarms')

         itimestep = itimestep + 1
         call mpas_advance_clock(domain % clock)

         currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
         call mpas_log_write('Doing timestep ' // trim(timeStamp))

         ! pre-timestep analysis computation
         call ocn_analysis_precompute(domain, ierr)

#ifdef MPAS_DEBUG
         call mpas_log_write( '   Computing surface flux arrays')
#endif
         block_ptr => domain % blocklist
         do while(associated(block_ptr))
           call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
           call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
           call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
           call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)
           call ocn_forcing_build_fraction_absorbed_array(meshPool, statePool, forcingPool, ierr, 1)
           call mpas_timer_start("land_ice_build_arrays")
           call ocn_surface_land_ice_fluxes_build_arrays(meshPool, &
                                                         forcingPool, scratchPool, statePool, dt, err)
           call mpas_timer_stop("land_ice_build_arrays")

           call ocn_frazil_forcing_build_arrays(domain, meshPool, forcingPool, statePool, err)
           call ocn_tidal_forcing_build_array(domain, meshPool, forcingPool, statePool, err)

           ! Compute normalGMBolusVelocity, relativeSlope and RediDiffVertCoef if respective flags are turned on

           call ocn_eddy_compute_mixed_layer_depth(statePool, forcingPool)
           if (config_use_GM .or. config_submesoscale_enable) then
               call mpas_timer_start("submesoscale eddy velocity compute", .false.)
               call ocn_eddy_compute_buoyancy_gradient()
               call mpas_timer_stop("submesoscale eddy velocity compute")
           end if

           if (config_use_Redi.or.config_use_GM) then
               call ocn_gm_compute_Bolus_velocity(statePool, meshPool, scratchPool, timeLevelIn=1)
           end if

           block_ptr => block_ptr % next
         end do

         call mpas_timer_start("time integration")

#ifdef MPAS_DEBUG
         call mpas_log_write( '   Computing forward time step')
#endif
         call ocn_timestep(domain, dt, timeStamp)

         call mpas_timer_stop("time integration")

         ! Move time level 2 fields back into time level 1 for next time step
#ifdef MPAS_DEBUG
         call mpas_log_write( '   Shifting time levels')
#endif
         call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
         call mpas_pool_shift_time_levels(statePool)

#ifdef MPAS_DEBUG
         call mpas_log_write( '   Handling analysis members')
#endif
         call ocn_analysis_compute(domain, err)
         call ocn_analysis_restart(domain, err)
         call ocn_analysis_write(domain, err)

#ifdef MPAS_DEBUG
         call mpas_log_write( '   Performing I/O')
#endif
         call mpas_timer_start('io_write')
         call mpas_stream_mgr_write(domain % streamManager, streamID='output', ierr=ierr)
         call mpas_timer_stop('io_write')
         call mpas_timer_start('reset_io_alarms')
         call mpas_stream_mgr_reset_alarms(domain % streamManager, streamID='output', ierr=ierr)
         call mpas_timer_stop('reset_io_alarms')

         if ( mpas_stream_mgr_ringing_alarms(domain % streamManager, streamID='restart', direction=MPAS_STREAM_OUTPUT, &
                                             ierr=ierr) ) then
#ifdef MPAS_DEBUG
            call mpas_log_write( '   Writing restart timestamp file')
#endif
            if ( domain % dminfo % my_proc_id == 0 ) then
               open(22, file=config_restart_timestamp_name, form='formatted', status='replace')
               write(22, *) trim(timeStamp)
               close(22)
            end if

            ! write forcing restart if any type(s) time varying forcing is used
            if (trim(config_sw_absorption_type)=='ohlmann00' .or. &
                config_use_ecosysTracers .or. &
                (config_use_activeTracers_surface_restoring .and. config_use_surface_salinity_monthly_restoring) .or. &
                config_use_time_varying_atmospheric_forcing .or. &
                config_use_time_varying_land_ice_forcing) then

                call ocn_time_varying_forcing_write_restart_times(domain)
            end if
         end if

         call mpas_timer_start('io_write')
         call mpas_stream_mgr_write(domain % streamManager, streamID='restart', ierr=ierr)
         call mpas_timer_stop('io_write')


         call mpas_timer_start('reset_io_alarms')
         call mpas_stream_mgr_reset_alarms(domain % streamManager, streamID='restart', ierr=ierr)
         call mpas_timer_stop('reset_io_alarms')

         call mpas_timer_start('io_write')
         call mpas_stream_mgr_write(domain % streamManager, ierr=ierr)
         call mpas_timer_stop('io_write')
         call mpas_timer_start('reset_io_alarms')
         call mpas_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)
         call mpas_timer_stop('reset_io_alarms')

         ! read next time level data required for variable shortwave
         call mpas_timer_start('io_shortwave',.false.)
         call ocn_get_shortWaveData(domain % streamManager, domain, domain % clock, .false.)
         call mpas_timer_stop('io_shortwave')

         ! read next time level data required for ecosys forcing
         if (config_use_ecosysTracers) then
           call mpas_timer_start('io_ecosys',.false.)
           call ocn_get_ecosysData(domain % streamManager, domain, domain % clock, .false.)
           call mpas_timer_stop('io_ecosys')
         endif

         ! read next time level data required for CFC forcing
         if (config_use_CFCTracers) then
           call mpas_timer_start('io_CFC',.false.)
           call ocn_get_CFCData(domain % streamManager, domain, domain % clock, .false.)
           call mpas_timer_stop('io_CFC')
         endif

         ! read next time level data required for monthly surface salinity restoring
         if (config_use_activeTracers_surface_restoring .and. config_use_surface_salinity_monthly_restoring) then
           if ( mpas_is_alarm_ringing(domain % clock, 'salinityDataReadAlarm', ierr=err_tmp) ) then
              call mpas_reset_clock_alarm(domain % clock, 'salinityDataReadAlarm', ierr=err_tmp)
              call mpas_timer_start('io_monthly_surface_salinity',.false.)
              call ocn_get_surfaceSalinityData(domain % streamManager, domain, domain % clock, .false.)
              call mpas_timer_stop('io_monthly_surface_salinity')
           endif
         endif

         ! read in next time level data required for time-varying forcing
         call ocn_time_varying_forcing_get(domain % streamManager, domain, domain % clock)

         ! Validate that the state is OK to run with for the next timestep.
         call ocn_validate_state(domain, timeLevel=1)

      end do

   end function ocn_forward_mode_run!}}}

!***********************************************************************
!
!  function ocn_forward_mode_finalize
!
!> \brief   Finalize MPAS-Ocean Forward Mode
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details
!>  This function finalizes the MPAS-Ocean core in forward mode.
!
!-----------------------------------------------------------------------

   function ocn_forward_mode_finalize(domain) result(iErr)!{{{

      type (domain_type), intent(inout) :: domain

      integer :: ierr

      ! destroy the ocean mesh structure
      call ocn_meshDestroy(ierr)

      if ( config_time_integrator == 'split_explicit') then
         call mpas_dmpar_exch_group_destroy_reusable_buffers(domain, 'subcycleFields')
      else if ( config_time_integrator == 'split_implicit') then
         call ocn_time_integrator_si_variable_destroy()
      end if

      call ocn_analysis_finalize(domain, ierr)

      call ocn_vmix_gotm_finalize(domain, ierr)

      call ocn_submesoscale_finalize()

      call mpas_destroy_clock(domain % clock, ierr)

      call mpas_decomp_destroy_decomp_list(domain % decompositions)

   end function ocn_forward_mode_finalize!}}}

end module ocn_forward_mode

! vim: foldmethod=marker
