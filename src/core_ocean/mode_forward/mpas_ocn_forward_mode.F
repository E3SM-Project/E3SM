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
   use ocn_tracer_advection
   use ocn_gm

   use ocn_high_freq_thickness_hmix_del2

   use ocn_equation_of_state

   use ocn_vmix

   use ocn_time_average

   use ocn_forcing
   use ocn_sea_ice

   use ocn_constants

   implicit none
   private

   public :: ocn_forward_mode_init, ocn_forward_mode_run, ocn_forward_mode_finalize
   public :: ocn_forward_mode_setup_clock

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

      type (domain_type), intent(inout) :: domain
      character(len=*), intent(out) :: startTimeStamp
      integer :: ierr

      real (kind=RKIND) :: dt
      type (block_type), pointer :: block

      integer :: err_tmp
      integer, pointer :: nVertLevels
      real (kind=RKIND) :: maxDensity, maxDensity_global
      real (kind=RKIND), dimension(:), pointer :: meshDensity
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: diagnosticsPool

      character (len=StrKIND), pointer :: xtime
      type (MPAS_Time_Type) :: startTime
      type (MPAS_TimeInterval_type) :: timeStep

      logical, pointer :: config_do_restart, config_filter_btr_mode, config_conduct_tests
      character (len=StrKIND), pointer :: config_vert_coord_movement, config_pressure_gradient_type
      real (kind=RKIND), pointer :: config_maxMeshDensity

      ierr = 0

      !
      ! Set startTimeStamp based on the start time of the simulation clock
      !
      startTime = mpas_get_clock_time(domain % clock, MPAS_START_TIME, err_tmp)
      call mpas_get_time(startTime, dateTimeString=startTimeStamp)
      ierr = ior(ierr, err_tmp)

      ! Setup ocean config pool
      call ocn_constants_init(domain % configs, domain % packages)

      call mpas_pool_get_config(domain % configs, 'config_do_restart', config_do_restart)
      call mpas_pool_get_config(domain % configs, 'config_vert_coord_movement', config_vert_coord_movement)
      call mpas_pool_get_config(domain % configs, 'config_pressure_gradient_type', config_pressure_gradient_type)
      call mpas_pool_get_config(domain % configs, 'config_filter_btr_mode', config_filter_btr_mode)
      call mpas_pool_get_config(domain % configs, 'config_maxMeshDensity', config_maxMeshDensity)
      call mpas_pool_get_config(domain % configs, 'config_conduct_tests', config_conduct_tests)

      !
      ! Read input data for model
      !
      call mpas_timer_start('io_read')
      call MPAS_stream_mgr_read(domain % streamManager, streamID='mesh', whence=MPAS_STREAM_NEAREST, ierr=err_tmp)
      if ( config_do_restart ) then
         call MPAS_stream_mgr_read(domain % streamManager, streamID='restart', ierr=err_tmp)
      else
         call MPAS_stream_mgr_read(domain % streamManager, streamID='input', ierr=err_tmp)
      end if
      call mpas_timer_stop('io_read')
      call mpas_timer_start('reset_io_alarms')
      call mpas_stream_mgr_reset_alarms(domain % streamManager, streamID='input', ierr=err_tmp)
      call mpas_stream_mgr_reset_alarms(domain % streamManager, streamID='restart', ierr=err_tmp)
      call mpas_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=err_tmp)
      call mpas_timer_stop('reset_io_alarms')

      ! Initialize submodules before initializing blocks.
      call ocn_timestep_init(ierr)

      call ocn_thick_hadv_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      call ocn_thick_vadv_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      call ocn_thick_surface_flux_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      call ocn_thick_ale_init(err_tmp)
      ierr = ior(ierr,err_tmp)

      call ocn_vel_coriolis_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      call ocn_vel_pressure_grad_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      call ocn_vel_vadv_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      call ocn_vel_hmix_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      call ocn_vel_forcing_init(err_tmp)
      ierr = ior(ierr, err_tmp)

      call ocn_tracer_hmix_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      call ocn_tracer_surface_flux_init(err_tmp)
      ierr = ior(ierr, err_tmp)
      call ocn_tracer_advection_init(err_tmp)
      ierr = ior(ierr,err_tmp)
      call ocn_tracer_short_wave_absorption_init(err_tmp)
      ierr = ior(ierr,err_tmp)
      call ocn_gm_init(err_tmp)
      ierr = ior(ierr,err_tmp)
      call ocn_tracer_nonlocalflux_init(err_tmp)
      ierr = ior(ierr,err_tmp)

      call ocn_vmix_init(domain, err_tmp)
      ierr = ior(ierr, err_tmp)

      call ocn_equation_of_state_init(err_tmp)
      ierr = ior(ierr, err_tmp)

      call ocn_tendency_init(err_tmp)
      ierr = ior(ierr,err_tmp)
      call ocn_diagnostics_init(err_tmp)
      ierr = ior(ierr,err_tmp)

      call ocn_forcing_init(err_tmp)
      ierr = ior(ierr,err_tmp)

      call ocn_high_freq_thickness_hmix_del2_init(err_tmp)
      ierr = ior(ierr,err_tmp)

      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nVertLevels', nVertLevels)
      call ocn_sea_ice_init(nVertLevels, err_tmp)
      ierr = ior(ierr, err_tmp)

      call ocn_analysis_init(domain, err_tmp)
      ierr = ior(ierr, err_tmp)

      if(ierr.eq.1) then
          call mpas_log_write('ERROR: An error was encountered while initializing the MPAS-Ocean forward mode')
      endif

      call ocn_init_routines_vert_coord(domain)

      call ocn_init_routines_compute_max_level(domain)

      if (.not.config_do_restart) call ocn_time_integration_split_init(domain)

      call mpas_log_write(' Vertical coordinate movement is: ' // trim(config_vert_coord_movement))

      if (config_vert_coord_movement.ne.'fixed'.and. &
          config_vert_coord_movement.ne.'uniform_stretching'.and. &
          config_vert_coord_movement.ne.'impermeable_interfaces'.and. &
          config_vert_coord_movement.ne.'user_specified') then
         call mpas_log_write(' Incorrect choice of config_vert_coord_movement.')
         call mpas_log_write('ERROR: Incorrect choice of config_vert_coord_movement.')
      endif

      if(config_vert_coord_movement .ne. 'impermeable_interfaces' .and. config_pressure_gradient_type .eq. 'MontgomeryPotential') then
         call mpas_log_write(' Incorrect combination of config_vert_coord_movement and config_pressure_gradient_type')
         call mpas_log_write('ERROR: Incorrect combination of config_vert_coord_movement and config_pressure_gradient_type')
      end if

      if (config_filter_btr_mode.and. &
          config_vert_coord_movement.ne.'fixed')then
         call mpas_log_write('ERROR: filter_btr_mode has only been tested with config_vert_coord_movement=fixed.')
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
      timeStep = mpas_get_clock_timestep(domain % clock, ierr=err_tmp)
      call mpas_get_timeInterval(timeStep, dt=dt)

      block => domain % blocklist
      do while (associated(block))
         call ocn_init_routines_block(block, dt, ierr)
         if(ierr.eq.1) then
             call mpas_log_write('ERROR: An error was encountered in ocn_init_routines_block')
         endif

         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_array(diagnosticsPool, 'xtime', xtime)
         xtime = startTimeStamp
         block => block % next
      end do

      if (config_conduct_tests) then
         call mpas_timer_start("test suite")
         call ocn_test_suite(domain,ierr)
         call mpas_timer_stop("test suite")
      endif

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
      integer :: ierr

      type (mpas_pool_type), pointer :: averagePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: forcingPool

      type (MPAS_timeInterval_type) :: timeStep
      character(len=StrKIND), pointer :: config_restart_timestamp_name
      logical, pointer :: config_write_output_on_startup

      ierr = 0

      call mpas_pool_get_config(domain % configs, 'config_write_output_on_startup', config_write_output_on_startup)
      call mpas_pool_get_config(domain % configs, 'config_restart_timestamp_name', config_restart_timestamp_name)

      ! Eventually, dt should be domain specific
      timeStep = mpas_get_clock_timestep(domain % clock, ierr=ierr)
      call mpas_get_timeInterval(timeStep, dt=dt)

      currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
      call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
      call mpas_log_write('Initial time ' // trim(timeStamp))

      call ocn_analysis_compute_startup(domain, err)

      if (config_write_output_on_startup) then
          call mpas_timer_start('io_write')
          call mpas_stream_mgr_write(domain % streamManager, 'output', forceWriteNow=.true., ierr=ierr)
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
      do while (.not. mpas_is_clock_stop_time(domain % clock))
         call mpas_timer_start('io_read')
         call mpas_stream_mgr_read(domain % streamManager, ierr=ierr)
         call mpas_timer_stop('io_read')
         call mpas_timer_start('reset_io_alarms')
         call mpas_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_INPUT, ierr=ierr)
         call mpas_timer_stop('reset_io_alarms')

         itimestep = itimestep + 1
         call mpas_advance_clock(domain % clock)

         currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr)
         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
         call mpas_log_write('Doing timestep ' // trim(timeStamp))

         block_ptr => domain % blocklist
         do while(associated(block_ptr))
           call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
           call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
           call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
           call ocn_forcing_build_arrays(meshPool, statePool, forcingPool, ierr, 1)
           call ocn_forcing_build_fraction_absorbed_array(meshPool, statePool, forcingpool, ierr, 1)
           block_ptr => block_ptr % next
         end do

         call mpas_timer_start("time integration")
         call ocn_timestep(domain, dt, timeStamp)
         call mpas_timer_stop("time integration")

         ! Move time level 2 fields back into time level 1 for next time step
         block_ptr => domain % blocklist
         do while(associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
            call mpas_pool_shift_time_levels(statePool)
            block_ptr => block_ptr % next
         end do

         call ocn_analysis_compute(domain, err) 
         call ocn_analysis_write(domain, err)

         call mpas_timer_start('io_write')
         call mpas_stream_mgr_write(domain % streamManager, streamID='output', ierr=ierr)
         call mpas_timer_stop('io_write')
         call mpas_timer_start('reset_io_alarms')
         call mpas_stream_mgr_reset_alarms(domain % streamManager, streamID='output', ierr=ierr)
         call mpas_timer_stop('reset_io_alarms')

         call mpas_timer_start('io_write')
         call mpas_stream_mgr_write(domain % streamManager, streamID='restart', ierr=ierr)
         call mpas_timer_stop('io_write')

         if ( mpas_stream_mgr_ringing_alarms(domain % streamManager, streamID='restart', direction=MPAS_STREAM_OUTPUT, ierr=ierr) ) then
            open(22, file=config_restart_timestamp_name, form='formatted', status='replace')
            write(22, *) trim(timeStamp)
            close(22)
         end if

         call mpas_timer_start('reset_io_alarms')
         call mpas_stream_mgr_reset_alarms(domain % streamManager, streamID='restart', ierr=ierr)
         call mpas_timer_stop('reset_io_alarms')

         call mpas_timer_start('io_write')
         call mpas_stream_mgr_write(domain % streamManager, ierr=ierr)
         call mpas_timer_stop('io_write')
         call mpas_timer_start('reset_io_alarms')
         call mpas_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)
         call mpas_timer_stop('reset_io_alarms')
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

      call ocn_analysis_finalize(domain, ierr)

      call mpas_destroy_clock(domain % clock, ierr)

      call mpas_decomp_destroy_decomp_list(domain % decompositions)

   end function ocn_forward_mode_finalize!}}}

end module ocn_forward_mode

! vim: foldmethod=marker
