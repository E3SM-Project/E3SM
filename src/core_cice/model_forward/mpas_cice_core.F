! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module cice_core

   use mpas_framework
   use mpas_timekeeping
   use cice_analysis_driver
   use cice_column

   private
   public :: &
        cice_core_init, &
        cice_core_run, &
        cice_core_finalize, &
        cice_simulation_clock_init

   type (MPAS_Clock_type), pointer :: clock

   contains

   !***********************************************************************
   !
   !  routine cice_core_init
   !
   !> \brief
   !> \author
   !> \date
   !> \details
   !>
   !
   !-----------------------------------------------------------------------

   function cice_core_init(domain, startTimeStamp) result(iErr)

      use mpas_derived_types
      use mpas_stream_manager

      use cice_diagnostics, only: &
           cice_initialize_time_diagnostics

      implicit none

      type (domain_type), intent(inout) :: domain
      character(len=*), intent(out) :: startTimeStamp

      real (kind=RKIND) :: dt

      character(len=strKIND), pointer :: config_initial_condition_type
      logical, pointer :: config_do_restart
      real (kind=RKIND), pointer :: config_dt
      type (MPAS_Time_Type) :: startTime
      integer :: ierr

      call mpas_pool_set_error_level(MPAS_POOL_WARN)

      iErr = 0

      ! write out namelist options
      write(stderrUnit,*) ; write(stderrUnit,*) " ----- namelist options used -----" ;  write(stderrUnit,*)
      call mpas_pool_print_summary(domain % configs, MPAS_POOL_CONFIG)
      write(stderrUnit,*) ; write(stderrUnit,*) " ----- done listing namelist options used -----" ; write(stderrUnit,*)

      call mpas_pool_get_config(domain % configs, 'config_initial_condition_type', config_initial_condition_type)
      call mpas_pool_get_config(domain % configs, 'config_do_restart', config_do_restart)
      call mpas_pool_get_config(domain % configs, 'config_dt', config_dt)

      !
      ! Set "local" clock to point to the clock contained in the domain type
      !
      clock => domain % clock


      !
      ! Set startTimeStamp based on the start time of the simulation clock
      !
      startTime = mpas_get_clock_time(clock, MPAS_START_TIME, ierr)
      call mpas_get_time(startTime, dateTimeString=startTimeStamp)

      !
      ! If this is a restart run, read the restart stream, else read the input stream.
      ! Regardless of which stream we read for initial conditions, reset the
      ! input alarms for both input and restart before reading any remaining input streams.
      !
      if (config_do_restart) then
         call MPAS_stream_mgr_read(domain % streamManager, streamID='restart', ierr=ierr)
      else
         if (trim(config_initial_condition_type) == "restart") then
            call MPAS_stream_mgr_read(domain % streamManager, streamID='restart_ic', ierr=ierr)
         else
            call MPAS_stream_mgr_read(domain % streamManager, streamID='input', ierr=ierr)
         end if
      endif
      call MPAS_stream_mgr_reset_alarms(domain % streamManager, streamID='input', direction=MPAS_STREAM_INPUT, ierr=ierr)
      call MPAS_stream_mgr_reset_alarms(domain % streamManager, streamID='restart', direction=MPAS_STREAM_INPUT, ierr=ierr)
      call MPAS_stream_mgr_reset_alarms(domain % streamManager, streamID='restart_ic', direction=MPAS_STREAM_INPUT, ierr=ierr)

      ! bootstrap analysis members
      call cice_analysis_bootstrap(domain, ierr)

      ! Read all other inputs
      call MPAS_stream_mgr_read(domain % streamManager, ierr=ierr)
      call MPAS_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_INPUT, ierr=ierr)

      !
      ! Initialize core
      !
      dt = config_dt

      ! initialize time diagnostics
      call cice_initialize_time_diagnostics(domain)

      ! initialize analysis members
      call cice_analysis_init(domain, ierr)

      ! bootstrap analysis has used config_do_restart: now modify for rest of code
      if (trim(config_initial_condition_type) == "restart") config_do_restart = .true.

      call mpas_init_block(domain, dt)

   end function cice_core_init

   !***********************************************************************
   !
   !  routine cice_simulation_clock_init
   !
   !> \brief
   !> \author
   !> \date
   !> \details
   !>
   !
   !-----------------------------------------------------------------------

   subroutine cice_simulation_clock_init(core_clock, configs, ierr)

      implicit none

      type (MPAS_Clock_type), intent(inout) :: core_clock
      type (mpas_pool_type), intent(inout) :: configs
      integer, intent(out) :: ierr

      type (MPAS_Time_Type) :: startTime, stopTime, alarmStartTime
      type (MPAS_TimeInterval_type) :: runDuration, timeStep, alarmTimeStep
      integer :: local_err

      character (len=StrKIND), pointer :: config_start_time, config_run_duration, config_stop_time, config_restart_timestamp_name
      character (len=StrKIND) :: restartTimeStamp
      real (kind=RKIND), pointer :: config_dt


      ierr = 0

      call mpas_pool_get_config(configs, 'config_dt', config_dt)
      call mpas_pool_get_config(configs, 'config_start_time', config_start_time)
      call mpas_pool_get_config(configs, 'config_run_duration', config_run_duration)
      call mpas_pool_get_config(configs, 'config_stop_time', config_stop_time)
      call mpas_pool_get_config(configs, 'config_restart_timestamp_name', config_restart_timestamp_name)


      ! Set time to the user-specified start time OR use a restart time from file
      if ( trim(config_start_time) == "file" ) then
         open(22, file=config_restart_timestamp_name, form='formatted', status='old')
         read(22,*) restartTimeStamp
         close(22)
         call mpas_set_time(curr_time=startTime, dateTimeString=restartTimeStamp, ierr=local_err)
      else
         call mpas_set_time(curr_time=startTime, dateTimeString=config_start_time, ierr=local_err)
      end if

      call mpas_set_timeInterval(timeStep, dt=config_dt, ierr=local_err)

      if (trim(config_run_duration) /= "none") then
         call mpas_set_timeInterval(runDuration, timeString=config_run_duration, ierr=local_err)
         call mpas_create_clock(core_clock, startTime=startTime, timeStep=timeStep, runDuration=runDuration, ierr=local_err)

         if (trim(config_stop_time) /= "none") then
            call mpas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, ierr=local_err)
            if(startTime + runduration /= stopTime) then
               write(0,*) 'Warning: config_run_duration and config_stop_time are inconsitent: using config_run_duration.'
            end if
         end if
      else if (trim(config_stop_time) /= "none") then
         call mpas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, ierr=local_err)
         call mpas_create_clock(core_clock, startTime=startTime, timeStep=timeStep, stopTime=stopTime, ierr=local_err)
      else
          write(0,*) 'Error: Neither config_run_duration nor config_stop_time were specified.'
          ierr = 1
      end if

   end subroutine cice_simulation_clock_init

   !***********************************************************************
   !
   !  routine mpas_init_block
   !
   !> \brief
   !> \author
   !> \date
   !> \details
   !>
   !
   !-----------------------------------------------------------------------

   subroutine mpas_init_block(domain, dt)

      use mpas_derived_types

      use cice_initialize, only: cice_init

      implicit none

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(in) :: dt

      ! initialize general stuff
      call cice_init(domain, clock, dt)

   end subroutine mpas_init_block

   !***********************************************************************
   !
   !  routine cice_core_run
   !
   !> \brief
   !> \author
   !> \date
   !> \details
   !>
   !
   !-----------------------------------------------------------------------

   function cice_core_run(domain) result(iErr)

      use mpas_derived_types
      use mpas_kind_types
      use mpas_stream_manager
      use mpas_timer
      use cice_time_integration, only: cice_timestep_finalize
      use cice_forcing, only: cice_forcing_get, cice_forcing_write_restart_times
      use cice_initialize, only: cice_init_post_clock_advance
      use cice_debug, only: istep

      implicit none

      type (domain_type), intent(inout) :: domain

      integer :: itimestep
      type (block_type), pointer :: block_ptr

      type (MPAS_Time_Type) :: currTime
      character(len=StrKIND) :: timeStamp
      integer :: ierr

      character(len=StrKIND), pointer :: config_restart_timestamp_name

      logical, pointer :: config_write_output_on_startup

      iErr = 0

      ! configurations
      call mpas_pool_get_config(domain % configs, 'config_restart_timestamp_name', config_restart_timestamp_name)

      currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
      call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
      write(0,*) 'Initial timestep ', trim(timeStamp)

      ! compute analysis at startup if required
      ! (note: might require some reordering of call here if cice_forcing_get data is used in analysis members)
      call cice_analysis_compute_startup(domain, ierr)

      ! Avoid writing a restart file at the initial time
      call MPAS_stream_mgr_reset_alarms(domain % streamManager, streamID='restart', direction=MPAS_STREAM_OUTPUT, ierr=ierr)

      call mpas_pool_get_config(domain % configs, 'config_write_output_on_startup', config_write_output_on_startup)
      if (config_write_output_on_startup) then
         call mpas_stream_mgr_write(domain % streamManager, ierr=ierr)
      endif
      call mpas_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)

      itimestep = 1
      call mpas_advance_clock(clock)

      ! initial forcing
      call cice_forcing_get(domain % streamManager, domain, clock, .true.)

      ! final initialization after clock advance
      call cice_init_post_clock_advance(domain, clock)

      istep = istep + 1

      timeLoop: do

         currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
         write(0,*) 'Doing timestep ', trim(timeStamp)

         ! pre-timestep analysis computation
         call cice_analysis_precompute(domain, ierr)

         ! do the time integration
         call mpas_timer_start("time integration")
         call mpas_timestep(domain, itimestep, timeStamp)
         call mpas_timer_stop("time integration")

         ! update analysis members
         call cice_analysis_compute(domain, ierr)
         call cice_analysis_restart(domain, ierr)
         call cice_analysis_write(domain, ierr)

         call mpas_stream_mgr_write(domain % streamManager, ierr=ierr)

         ! update the restart_timestamp file with the new time if needed
         if ( mpas_stream_mgr_ringing_alarms(domain % streamManager, &
              streamID='restart', direction=MPAS_STREAM_OUTPUT, ierr=ierr) ) then
            open(22, file=config_restart_timestamp_name, form='formatted', status='replace')
            write(22, *) timeStamp
            close(22)
            call cice_forcing_write_restart_times(domain)
         end if

         call mpas_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)

         ! finalize the timstep
         call cice_timestep_finalize(domain)

         ! check if the simulation is finished
         if (mpas_is_clock_stop_time(clock)) exit timeLoop

         ! advance the clock
         itimestep = itimestep + 1
         call mpas_advance_clock(clock)
         istep = istep + 1

         ! forcing
         call cice_forcing_get(domain % streamManager, domain, clock, .false.)

      end do timeLoop

   end function cice_core_run

   !***********************************************************************
   !
   !  routine mpas_timestep
   !
   !> \brief
   !> \author
   !> \date
   !> \details
   !>
   !
   !-----------------------------------------------------------------------

   subroutine mpas_timestep(domain, itimestep, timeStamp)

      use mpas_derived_types
      use cice_time_integration
      use mpas_timer

      implicit none

      type (domain_type), intent(inout) :: domain
      integer, intent(in) :: itimestep
      character(len=*), intent(in) :: timeStamp

      call cice_timestep(domain, clock, itimestep)

   end subroutine mpas_timestep

   !***********************************************************************
   !
   !  routine cice_core_finalize
   !
   !> \brief
   !> \author
   !> \date
   !> \details
   !>
   !
   !-----------------------------------------------------------------------

   function cice_core_finalize(domain) result(iErr)

      use mpas_derived_types
      use mpas_decomp
      use cice_column, only: &
           cice_column_finalize

      implicit none

      type (domain_type), intent(inout) :: domain
      integer :: ierr

      iErr = 0

      ! finalize column
      call cice_column_finalize(domain)

      call cice_analysis_finalize(domain, ierr)

      call mpas_destroy_clock(clock, ierr)

      call mpas_decomp_destroy_decomp_list(domain % decompositions)

   end function cice_core_finalize

   !-----------------------------------------------------------------------

end module cice_core
