! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
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

   use mpas_framework
   use mpas_timekeeping
   use mpas_dmpar
   use mpas_timer
   use mpas_io_units

   use ocn_analysis_driver
   use ocn_init_routines

   use ocn_global_diagnostics
   use ocn_time_integration
   use ocn_tendency
   use ocn_diagnostics
   use ocn_test

   use ocn_thick_hadv
   use ocn_thick_vadv
   use ocn_thick_ale
   use ocn_thick_surface_flux

   use ocn_vel_pressure_grad
   use ocn_vel_vadv
   use ocn_vel_hmix
   use ocn_vel_forcing
   use ocn_vel_coriolis

   use ocn_tracer_hmix
   use ocn_tracer_surface_flux
   use ocn_tracer_short_wave_absorption
   use ocn_tracer_nonlocalflux
   use ocn_gm

   use ocn_high_freq_thickness_hmix_del2

   use ocn_equation_of_state

   use ocn_vmix

   use ocn_time_average

   use ocn_forcing
   use ocn_sea_ice

   use ocn_constants

   type (MPAS_Clock_type), pointer :: clock

   type (timer_node), pointer :: globalDiagTimer, timeIntTimer, testSuiteTimer

   contains

!***********************************************************************
!
!  routine ocn_forward_mode_init
!
!> \brief   Initialize MPAS-Ocean core
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details 
!>  This routine calls all initializations required to begin a
!>  simulation with MPAS-Ocean
!
!-----------------------------------------------------------------------

   subroutine ocn_forward_mode_init(domain, stream_manager, startTimeStamp)!{{{

      use mpas_grid_types
      use mpas_stream_manager
      use ocn_tracer_advection
      use ocn_time_integration_split

      implicit none

      type (domain_type), intent(inout) :: domain
      type (MPAS_streamManager_type), intent(inout) :: stream_manager
      character(len=*), intent(out) :: startTimeStamp

      real (kind=RKIND) :: dt
      type (block_type), pointer :: block
      type (dm_info) :: dminfo

      integer :: err, err_tmp
      integer, pointer :: nVertLevels
      real (kind=RKIND) :: maxDensity, maxDensity_global
      real (kind=RKIND), dimension(:), pointer :: meshDensity
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: diagnosticsPool

      character (len=StrKIND), pointer :: xtime
      type (MPAS_Time_Type) :: startTime
      type (MPAS_TimeInterval_type) :: timeStep

      logical, pointer :: config_do_restart, config_filter_btr_mode, config_conduct_tests
      logical, pointer :: config_write_stats_on_startup
      character (len=StrKIND), pointer :: config_vert_coord_movement, config_pressure_gradient_type
      real (kind=RKIND), pointer :: config_maxMeshDensity

      err = 0

      !
      ! Set "local" clock to point to the clock contained in the domain type
      !
      clock => domain % clock

      !
      ! Set startTimeStamp based on the start time of the simulation clock
      !
      startTime = mpas_get_clock_time(clock, MPAS_START_TIME, err_tmp)
      call mpas_get_time(startTime, dateTimeString=startTimeStamp) 
      err = ior(err, err_tmp)

      ! Setup ocean config pool
      call ocn_constants_init(domain % configs, domain % packages)

      dminfo = domain % dminfo

      call mpas_pool_get_config(domain % configs, 'config_do_restart', config_do_restart)
      call mpas_pool_get_config(domain % configs, 'config_vert_coord_movement', config_vert_coord_movement)
      call mpas_pool_get_config(domain % configs, 'config_pressure_gradient_type', config_pressure_gradient_type)
      call mpas_pool_get_config(domain % configs, 'config_filter_btr_mode', config_filter_btr_mode)
      call mpas_pool_get_config(domain % configs, 'config_maxMeshDensity', config_maxMeshDensity)
      call mpas_pool_get_config(domain % configs, 'config_conduct_tests', config_conduct_tests)
      call mpas_pool_get_config(domain % configs, 'config_write_stats_on_startup', config_write_stats_on_startup)

      !
      ! Read input data for model
      !
      if ( config_do_restart ) then
         call mpas_timer_start('io_read', .false.)
         call MPAS_stream_mgr_read(stream_manager, streamID='restart', ierr=err)
         call mpas_timer_stop('io_read')
      else
         call mpas_timer_start('io_read', .false.)
         call MPAS_stream_mgr_read(stream_manager, streamID='input', ierr=err)
         call mpas_timer_stop('io_read')
      end if
      call mpas_timer_start('reset_io_alarms', .false.)
      call mpas_stream_mgr_reset_alarms(stream_manager, streamID='input', ierr=err)
      call mpas_stream_mgr_reset_alarms(stream_manager, streamID='restart', ierr=err)
      call mpas_stream_mgr_reset_alarms(stream_manager, direction=MPAS_STREAM_OUTPUT, ierr=err)
      call mpas_timer_stop('reset_io_alarms')

      ! Initialize submodules before initializing blocks.
      call ocn_timestep_init(err)

      call ocn_thick_hadv_init(err_tmp)
      err = ior(err, err_tmp)
      call ocn_thick_vadv_init(err_tmp)
      err = ior(err, err_tmp)
      call ocn_thick_surface_flux_init(err_tmp)
      err = ior(err, err_tmp)
      call ocn_thick_ale_init(err_tmp)
      err = ior(err,err_tmp)

      call ocn_vel_coriolis_init(err_tmp)
      err = ior(err, err_tmp)
      call ocn_vel_pressure_grad_init(err_tmp)
      err = ior(err, err_tmp)
      call ocn_vel_vadv_init(err_tmp)
      err = ior(err, err_tmp)
      call ocn_vel_hmix_init(err_tmp)
      err = ior(err, err_tmp)
      call ocn_vel_forcing_init(err_tmp)
      err = ior(err, err_tmp)

      call ocn_tracer_hmix_init(err_tmp)
      err = ior(err, err_tmp)
      call ocn_tracer_surface_flux_init(err_tmp)
      err = ior(err, err_tmp)
      call ocn_tracer_advection_init(err_tmp)
      err = ior(err,err_tmp)
      call ocn_tracer_short_wave_absorption_init(err_tmp)
      err = ior(err,err_tmp)
      call ocn_gm_init(err_tmp)
      err = ior(err,err_tmp)
      call ocn_tracer_nonlocalflux_init(err_tmp)
      err = ior(err,err_tmp)

      call ocn_vmix_init(domain, err_tmp)
      err = ior(err, err_tmp)

      call ocn_equation_of_state_init(err_tmp)
      err = ior(err, err_tmp)

      call ocn_tendency_init(err_tmp)
      err = ior(err,err_tmp)
      call ocn_diagnostics_init(err_tmp)
      err = ior(err,err_tmp)

      call ocn_forcing_init(err_tmp)
      err = ior(err,err_tmp)

      call ocn_high_freq_thickness_hmix_del2_init(err_tmp)
      err = ior(err,err_tmp)

      call ocn_global_diagnostics_init(dminfo,err_tmp)
      err = ior(err, err_tmp)

      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nVertLevels', nVertLevels)
      call ocn_sea_ice_init(nVertLevels, err_tmp)
      err = ior(err, err_tmp)

      call ocn_analysis_init(domain, err_tmp)
      err = ior(err, err_tmp)

      call mpas_timer_init(domain)

      if(err.eq.1) then
          call mpas_dmpar_abort(dminfo)
      endif

      call ocn_init_routines_vert_coord(domain)

      call ocn_init_routines_compute_max_level(domain)

      if (.not.config_do_restart) call ocn_time_integration_split_init(domain)

      write (stdoutUnit,'(a,a)') ' Vertical coordinate movement is: ',trim(config_vert_coord_movement)

      if (config_vert_coord_movement.ne.'fixed'.and. &
          config_vert_coord_movement.ne.'uniform_stretching'.and. &
          config_vert_coord_movement.ne.'impermeable_interfaces'.and. &
          config_vert_coord_movement.ne.'user_specified') then
         write (stderrUnit,*) ' Incorrect choice of config_vert_coord_movement.'
         call mpas_dmpar_abort(dminfo)
      endif

      if(config_vert_coord_movement .ne. 'impermeable_interfaces' .and. config_pressure_gradient_type .eq. 'MontgomeryPotential') then
         write (stderrUnit,*) ' Incorrect combination of config_vert_coord_movement and config_pressure_gradient_type'
         call mpas_dmpar_abort(dminfo)
      end if

      if (config_filter_btr_mode.and. &
          config_vert_coord_movement.ne.'fixed')then
         write (stderrUnit,*) 'filter_btr_mode has only been tested with'// &
            ' config_vert_coord_movement=fixed.'
         call mpas_dmpar_abort(dminfo)
      endif

      ! find the maximum value of the meshDensity
      if (config_maxMeshDensity < 0.0) then
        maxDensity=-1
        block => domain % blocklist
        do while (associated(block))
          call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
          call mpas_pool_get_array(meshPool, 'meshDensity', meshDensity)
          maxDensity = max(maxDensity, maxval(meshDensity))
          block => block % next
        end do
        call mpas_dmpar_max_real(domain % dminfo, maxDensity, maxDensity_global)
        config_maxMeshDensity = maxDensity_global
      endif

      !
      ! Initialize core
      !
      timeStep = mpas_get_clock_timestep(clock, ierr=err)
      call mpas_get_timeInterval(timeStep, dt=dt)

      block => domain % blocklist
      do while (associated(block))
         call ocn_init_routines_block(block, dt, err)
         if(err.eq.1) then
             call mpas_dmpar_abort(dminfo)
         endif

         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_array(diagnosticsPool, 'xtime', xtime)
         xtime = startTimeStamp 
         block => block % next
      end do

      if (config_conduct_tests) then
         call mpas_timer_start("test suite", .false., testSuiteTimer)
         call ocn_test_suite(domain,err)
         call mpas_timer_stop("test suite", testSuiteTimer)
      endif

      if (config_write_stats_on_startup) then
         call mpas_timer_start("global diagnostics", .false., globalDiagTimer)
         call ocn_compute_global_diagnostics(domain, 1 , 0, dt)
         call mpas_timer_stop("global diagnostics", globalDiagTimer)
      endif

   end subroutine ocn_forward_mode_init!}}}

!***********************************************************************
!
!  routine ocn_forward_mode_simulation_clock_init
!
!> \brief   Initialize timer variables
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details 
!>  This routine initializes all timer variables
!
!-----------------------------------------------------------------------

   subroutine ocn_forward_mode_simulation_clock_init(core_clock, configs, ierr)!{{{

      implicit none

      type (MPAS_Clock_type), intent(inout) :: core_clock
      type (mpas_pool_type), intent(inout) :: configs
      integer, intent(out) :: ierr

      type (MPAS_Time_Type) :: startTime, stopTime, alarmStartTime
      type (MPAS_TimeInterval_type) :: runDuration, timeStep, alarmTimeStep
      character(len=StrKIND) :: restartTimeStamp
      character(len=StrKIND), pointer :: config_start_time, config_stop_time, config_run_duration
      character(len=StrKIND), pointer :: config_stats_interval, config_dt, config_restart_timestamp_name
      integer :: err_tmp


      ierr = 0

      call mpas_pool_get_config(configs, 'config_dt', config_dt)
      call mpas_pool_get_config(configs, 'config_start_time', config_start_time)
      call mpas_pool_get_config(configs, 'config_stop_time', config_stop_time)
      call mpas_pool_get_config(configs, 'config_run_duration', config_run_duration)
      call mpas_pool_get_config(configs, 'config_stats_interval', config_stats_interval)
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
               write(stderrUnit,*) 'Warning: config_run_duration and config_stop_time are inconsitent: using config_run_duration.'
            end if
         end if
      else if (trim(config_stop_time) /= "none") then
         call mpas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, ierr=err_tmp)
         call mpas_create_clock(core_clock, startTime=startTime, timeStep=timeStep, stopTime=stopTime, ierr=err_tmp)
      else
          write(stderrUnit,*) 'Error: Neither config_run_duration nor config_stop_time were specified.'
          ierr = 1
      end if

      if (trim(config_stats_interval) /= "none") then      
         call mpas_set_timeInterval(alarmTimeStep, timeString=config_stats_interval, ierr=err_tmp)
         alarmStartTime = startTime + alarmTimeStep
         call mpas_add_clock_alarm(core_clock, statsAlarmID, alarmStartTime, alarmTimeStep, ierr=err_tmp)
      end if

   end subroutine ocn_forward_mode_simulation_clock_init!}}}

!***********************************************************************
!
!  routine ocn_forward_mode_run
!
!> \brief   Main driver for MPAS-Ocean time-stepping
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details 
!>  This routine includes the time-stepping loop, and calls timer
!>  routines to write output and restart files.
!
!-----------------------------------------------------------------------

   subroutine ocn_forward_mode_run(domain, stream_manager)!{{{
   
      use mpas_kind_types
      use mpas_grid_types
      use mpas_stream_manager
      use mpas_timer
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      type (MPAS_streamManager_type), intent(inout) :: stream_manager
   
      integer :: itimestep, err
      real (kind=RKIND) :: dt
      type (block_type), pointer :: block_ptr

      type (MPAS_Time_Type) :: currTime
      character(len=StrKIND) :: timeStamp
      integer :: ierr

      type (mpas_pool_type), pointer :: averagePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: forcingPool

      type (MPAS_timeInterval_type) :: timeStep
      character(len=StrKIND), pointer :: config_restart_timestamp_name
      logical, pointer :: config_write_output_on_startup
   
      call mpas_pool_get_config(domain % configs, 'config_write_output_on_startup', config_write_output_on_startup)
      call mpas_pool_get_config(domain % configs, 'config_restart_timestamp_name', config_restart_timestamp_name)

      ! Eventually, dt should be domain specific
      timeStep = mpas_get_clock_timestep(clock, ierr=ierr)
      call mpas_get_timeInterval(timeStep, dt=dt)

      currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
      call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
      write(stderrUnit,*) 'Initial time ', trim(timeStamp)

      call ocn_analysis_compute_startup(domain, stream_manager, err) 

      if (config_write_output_on_startup) then
          call mpas_timer_start('io_write', .false.)
          call mpas_stream_mgr_write(stream_manager, 'output', forceWriteNow=.true., ierr=ierr)
          call mpas_timer_stop('io_write')
      endif

      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'average', averagePool)
        call ocn_time_average_init(averagePool)
        block_ptr => block_ptr % next
      end do

      ! During integration, time level 1 stores the model state at the beginning of the
      !   time step, and time level 2 stores the state advanced dt in time by timestep(...)
      itimestep = 0
      do while (.not. mpas_is_clock_stop_time(clock))
         call mpas_timer_start('io_read', .false.)
         call mpas_stream_mgr_read(stream_manager, ierr=ierr)
         call mpas_timer_stop('io_read')
         call mpas_timer_start('reset_io_alarms', .false.)
         call mpas_stream_mgr_reset_alarms(stream_manager, direction=MPAS_STREAM_INPUT, ierr=ierr)
         call mpas_timer_stop('reset_io_alarms')

         itimestep = itimestep + 1
         call mpas_advance_clock(clock)

         currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
         write(stderrUnit,*) 'Doing timestep ', trim(timeStamp)
   
         block_ptr => domain % blocklist
         do while(associated(block_ptr))
           call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
           call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
           call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
           call ocn_forcing_build_arrays(meshPool, statePool, forcingPool, ierr, 1)
           call ocn_forcing_build_transmission_array(meshPool, statePool, forcingpool, ierr, 1)
           block_ptr => block_ptr % next
         end do

         call mpas_timer_start("time integration", .false., timeIntTimer)
         call ocn_forward_mode_timestep(domain, itimestep, dt, timeStamp)
         call mpas_timer_stop("time integration", timeIntTimer)
   
         ! Move time level 2 fields back into time level 1 for next time step
         block_ptr => domain % blocklist
         do while(associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
            call mpas_pool_shift_time_levels(statePool)
            block_ptr => block_ptr % next
         end do
      
         call ocn_analysis_compute_w_alarms(stream_manager, domain, err) 
         call ocn_analysis_write(stream_manager, err)

         call mpas_timer_start('io_write', .false.)
         call mpas_stream_mgr_write(stream_manager, streamID='output', ierr=ierr)
         call mpas_timer_stop('io_write')
         call mpas_timer_start('reset_io_alarms', .false.)
         call mpas_stream_mgr_reset_alarms(stream_manager, streamID='output', ierr=ierr)
         call mpas_timer_stop('reset_io_alarms')

         call mpas_timer_start('io_write', .false.)
         call mpas_stream_mgr_write(stream_manager, streamID='restart', ierr=ierr)
         call mpas_timer_stop('io_write')

         if ( mpas_stream_mgr_ringing_alarms(stream_manager, streamID='restart', direction=MPAS_STREAM_OUTPUT, ierr=ierr) ) then
            open(22, file=config_restart_timestamp_name, form='formatted', status='replace')
            write(22, *) trim(timeStamp)
            close(22)
         end if

         call mpas_timer_start('reset_io_alarms', .false.)
         call mpas_stream_mgr_reset_alarms(stream_manager, streamID='restart', ierr=ierr)
         call mpas_timer_stop('reset_io_alarms')

         call mpas_timer_start('io_write', .false.)
         call mpas_stream_mgr_write(stream_manager, ierr=ierr)
         call mpas_timer_stop('io_write')
         call mpas_timer_start('reset_io_alarms', .false.)
         call mpas_stream_mgr_reset_alarms(stream_manager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)
         call mpas_timer_stop('reset_io_alarms')
      end do

   end subroutine ocn_forward_mode_run!}}}
   
!***********************************************************************
!
!  routine ocn_forward_mode_timestep
!
!> \brief   Sub-driver for MPAS-Ocean time-stepping
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details 
!>  This routine calls the time integration routine within a time-stepping loop.
!
!-----------------------------------------------------------------------

   subroutine ocn_forward_mode_timestep(domain, itimestep, dt, timeStamp)!{{{
   
      use mpas_kind_types
      use mpas_grid_types
   
      implicit none
   
      type (domain_type), intent(inout) :: domain 
      integer, intent(in) :: itimestep
      real (kind=RKIND), intent(in) :: dt
      character(len=*), intent(in) :: timeStamp

      type (block_type), pointer :: block_ptr
      integer :: ierr
   
      call ocn_timestep(domain, dt, timeStamp)

      !if (config_stats_interval > 0) then
      !    if (mod(itimestep, config_stats_interval) == 0) then
      !       call mpas_timer_start("global diagnostics", .false., globalDiagTimer)
      !       call ocn_compute_global_diagnostics(domain, 2, itimestep, dt);
      !       call mpas_timer_stop("global diagnostics", globalDiagTimer)
      !    end if
      !end if

      !TODO: replace the above code block with this if we desire to convert config_stats_interval to use alarms
      if (mpas_is_alarm_ringing(clock, statsAlarmID, ierr=ierr)) then
         call mpas_reset_clock_alarm(clock, statsAlarmID, ierr=ierr)

!        block_ptr => domain % blocklist
!        if (associated(block_ptr % next)) then
!           write(stderrUnit,*) 'Error: computeGlobalDiagnostics assumes ',&
!                      'that there is only one block per processor.'
!        end if
   
         call mpas_timer_start("global diagnostics")
         call ocn_compute_global_diagnostics(domain, 2, itimestep, dt);
      !  call ocn_compute_global_diagnostics(domain % dminfo, &
      !           block_ptr % state % time_levs(2) % state, block_ptr % mesh, &
      !           timeStamp, dt)
         call mpas_timer_stop("global diagnostics")
      end if

   end subroutine ocn_forward_mode_timestep!}}}

   subroutine ocn_forward_mode_finalize(domain, stream_manager)!{{{
   
      use mpas_grid_types
      use mpas_decomp
      use mpas_stream_manager
   
      implicit none

      type (domain_type), intent(inout) :: domain 
      type (MPAS_streamManager_type), intent(inout) :: stream_manager
      integer :: ierr

      call ocn_analysis_finalize(domain, ierr) 

      call mpas_destroy_clock(clock, ierr)

      call mpas_decomp_destroy_decomp_list(decompositions)

   end subroutine ocn_forward_mode_finalize!}}}

!***********************************************************************
!
!  routine ocn_forward_mode_setup_clock
!
!> \brief   Pacakge setup routine
!> \author  Michael Duda
!> \date    6 August 2014
!> \details 
!>  The purpose of this routine is to allow the core to set up a simulation
!>  clock that will be used by the I/O subsystem for timing reads and writes
!>  of I/O streams.
!>  This routine is called from the superstructure after the framework 
!>  has been initialized but before any fields have been allocated and 
!>  initial fields have been read from input files. However, all namelist
!>  options are available.
!
!-----------------------------------------------------------------------
   subroutine ocn_forward_mode_setup_clock(core_clock, configs, ierr)!{{{

      implicit none

      type (MPAS_Clock_type), intent(inout) :: core_clock
      type (mpas_pool_type), intent(inout) :: configs
      integer, intent(out) :: ierr

      call ocn_forward_mode_simulation_clock_init(core_clock, configs, ierr)

   end subroutine ocn_forward_mode_setup_clock!}}}

end module ocn_forward_mode

! vim: foldmethod=marker
