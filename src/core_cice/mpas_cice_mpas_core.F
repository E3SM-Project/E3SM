module mpas_core

   use mpas_framework
   use mpas_timekeeping

   type (io_output_object), save :: restart_obj
   integer :: current_outfile_frames

   type (MPAS_Clock_type) :: clock

   integer, parameter :: outputAlarmID = 1
   integer, parameter :: restartAlarmID = 2
   !integer, parameter :: statsAlarmID = 3

   contains

!--------------------------------------------------------------------------

   subroutine mpas_core_init(domain, startTimeStamp)
   
      use mpas_configure
      use mpas_grid_types
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      character(len=*), intent(out) :: startTimeStamp
   
      real (kind=RKIND) :: dt
      type (block_type), pointer :: block

      !
      ! Initialize core
      !
      dt = config_dt

      call simulation_clock_init(domain, dt, startTimeStamp)

      block => domain % blocklist
      do while (associated(block))
         call mpas_init_block(block, block % mesh, dt)
         block % state % time_levs(1) % state % xtime % scalar = startTimeStamp
         block => block % next
      end do

      current_outfile_frames = 0

   end subroutine mpas_core_init

!--------------------------------------------------------------------------

   subroutine simulation_clock_init(domain, dt, startTimeStamp)

      implicit none

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(in) :: dt
      character(len=*), intent(out) :: startTimeStamp

      type (MPAS_Time_Type) :: startTime, stopTime, alarmStartTime
      type (MPAS_TimeInterval_type) :: runDuration, timeStep, alarmTimeStep
      integer :: ierr

      call mpas_set_time(curr_time=startTime, dateTimeString=config_start_time, ierr=ierr)
      call mpas_set_timeInterval(timeStep, dt=dt, ierr=ierr)

      if (trim(config_run_duration) /= "none") then
         call mpas_set_timeInterval(runDuration, timeString=config_run_duration, ierr=ierr)
         call mpas_create_clock(clock, startTime=startTime, timeStep=timeStep, runDuration=runDuration, ierr=ierr)

         if (trim(config_stop_time) /= "none") then
            call mpas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, ierr=ierr)
            if(startTime + runduration /= stopTime) then
               write(0,*) 'Warning: config_run_duration and config_stop_time are inconsitent: using config_run_duration.'
            end if
         end if
      else if (trim(config_stop_time) /= "none") then
         call mpas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, ierr=ierr)
         call mpas_create_clock(clock, startTime=startTime, timeStep=timeStep, stopTime=stopTime, ierr=ierr)
      else
          write(0,*) 'Error: Neither config_run_duration nor config_stop_time were specified.'
          call mpas_dmpar_abort(domain % dminfo)
      end if

      ! set output alarm
      call mpas_set_timeInterval(alarmTimeStep, timeString=config_output_interval, ierr=ierr)
      alarmStartTime = startTime + alarmTimeStep
      call mpas_add_clock_alarm(clock, outputAlarmID, alarmStartTime, alarmTimeStep, ierr=ierr)

      ! set restart alarm, if necessary
      if (trim(config_restart_interval) /= "none") then
         call mpas_set_timeInterval(alarmTimeStep, timeString=config_restart_interval, ierr=ierr)
         alarmStartTime = startTime + alarmTimeStep
         call mpas_add_clock_alarm(clock, restartAlarmID, alarmStartTime, alarmTimeStep, ierr=ierr)
      end if

      call mpas_get_time(curr_time=startTime, dateTimeString=startTimeStamp, ierr=ierr)

   end subroutine simulation_clock_init

!--------------------------------------------------------------------------

   subroutine mpas_init_block(block, mesh, dt)
   
      use mpas_grid_types
      use mpas_rbf_interpolation
      use mpas_vector_reconstruction

      use cice_advection, only: init_advection

      use cice_dynamics, only: init_dynamics

      use cice_testing, only: plot_pop_grid_cells, init_spherical_test_case

      use cice_forcing, only: init_forcing, get_forcing

      use cice_therm_zerolayer, only: init_semtner_forcing, init_semtner_test_case

      implicit none

      type (block_type), intent(inout) :: block
      type (mesh_type), intent(inout) :: mesh
      real (kind=RKIND), intent(in) :: dt

      call mpas_rbf_interp_initialize(mesh)
      call mpas_init_reconstruct(mesh)

      ! init advection
      call init_advection()

      ! init forcing
      call init_forcing(mesh, &
           block % atmos_forcing)

      call get_forcing(block % mesh, clock, block % atmos_forcing)

      call init_spherical_test_case(block % mesh, &
           block % icestate % tracers % array(index_iceAreaCell, 1, :), &
           block % icestate % tracers % array(index_iceVolumeCell, 1, :), &
           block % icestate % totalMassCell % array, &
           block % icestate % surfaceTemperature % array, &
           block % atmos_forcing % airTemperature % array)

      call init_dynamics(block, dt)

      call plot_pop_grid_cells(mesh, block % icestate % tracers % array(index_iceAreaCell, 1, :), "iceAreaCell_ic.txt")
      call plot_pop_grid_cells(mesh, block % icestate % tracers % array(index_iceVolumeCell, 1, :), "iceVolumeCell_ic.txt")

   end subroutine mpas_init_block
   
!--------------------------------------------------------------------------
   
   subroutine mpas_core_run(domain, output_obj, output_frame)
   
      use mpas_grid_types
      use mpas_kind_types
      use mpas_io_output
      use mpas_timer
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      type (io_output_object), intent(inout) :: output_obj
      integer, intent(inout) :: output_frame

      integer :: itimestep
      real (kind=RKIND) :: dt
      type (block_type), pointer :: block_ptr

      type (MPAS_Time_Type) :: currTime
      character(len=StrKIND) :: timeStamp
      integer :: ierr

      ! Eventually, dt should be domain specific
      dt = config_dt

      currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
      call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)         
      write(0,*) 'Initial timestep ', trim(timeStamp)

      call write_output_frame(output_obj, output_frame, domain)

      ! During integration, time level 1 stores the model state at the beginning of the
      !   time step, and time level 2 stores the state advanced dt in time by timestep(...)
      itimestep = 0
      do while (.not. mpas_is_clock_stop_time(clock))

         itimestep = itimestep + 1
         call mpas_advance_clock(clock)

         currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)         
         write(0,*) 'Doing timestep ', trim(timeStamp)

         call mpas_timer_start("time integration")
         call mpas_timestep(domain, itimestep, dt, timeStamp)
         call mpas_timer_stop("time integration")

         ! Move time level 2 fields back into time level 1 for next time step
         block_ptr => domain % blocklist
         do while(associated(block_ptr))
            call mpas_shift_time_levels_state(block_ptr % state)
            block_ptr => block_ptr % next
         end do

         !TODO: mpas_get_clock_ringing_alarms is probably faster than multiple mpas_is_alarm_ringing...

         if (mpas_is_alarm_ringing(clock, outputAlarmID, ierr=ierr)) then
            call mpas_reset_clock_alarm(clock, outputAlarmID, ierr=ierr)
            ! output_frame will always be > 1 here unless it was reset after the maximum number of frames per outfile was reached
            if(output_frame == 1) then
               call mpas_output_state_finalize(output_obj, domain % dminfo)
               call mpas_output_state_init(output_obj, domain, "OUTPUT", trim(timeStamp))
            end if
            call write_output_frame(output_obj, output_frame, domain)
         end if

         if (mpas_is_alarm_ringing(clock, restartAlarmID, ierr=ierr)) then
            call mpas_reset_clock_alarm(clock, restartAlarmID, ierr=ierr)

            ! Write one restart time per file
            call mpas_output_state_init(restart_obj, domain, "RESTART", trim(timeStamp))
            call mpas_output_state_for_domain(restart_obj, domain, 1)
            call mpas_output_state_finalize(restart_obj, domain % dminfo)
         end if

      end do

   end subroutine mpas_core_run
   
!--------------------------------------------------------------------------
   
   subroutine write_output_frame(output_obj, output_frame, domain)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute diagnostic fields for a domain and write model state to output file
   !
   ! Input/Output: domain - contains model state; diagnostic field are computed
   !                        before returning
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      use mpas_grid_types
      use mpas_io_output
   
      implicit none

      type (io_output_object), intent(inout) :: output_obj
      integer, intent(inout) :: output_frame
      type (domain_type), intent(inout) :: domain
   
      integer :: i, j, k
      integer :: eoe
      type (block_type), pointer :: block_ptr
   
      block_ptr => domain % blocklist
      do while (associated(block_ptr))
         call compute_output_diagnostics(block_ptr % state % time_levs(1) % state, block_ptr % mesh)
         block_ptr => block_ptr % next
      end do
   
      call mpas_output_state_for_domain(output_obj, domain, output_frame)
      output_frame = output_frame + 1

      ! reset frame if the maximum number of frames per outfile has been reached
      if (config_frames_per_outfile > 0) then
         current_outfile_frames = current_outfile_frames + 1            
         if(current_outfile_frames >= config_frames_per_outfile) then
            current_outfile_frames = 0
            output_frame = 1
         end if
      end if

   end subroutine write_output_frame
   
!--------------------------------------------------------------------------
   
   subroutine compute_output_diagnostics(state, grid)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute diagnostic fields for a domain
   !
   ! Input: state - contains model prognostic fields
   !        grid  - contains grid metadata
   !
   ! Output: state - upon returning, diagnostic fields will have be computed
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      use mpas_grid_types
   
      implicit none
   
      type (state_type), intent(inout) :: state
      type (mesh_type), intent(in) :: grid
   
      integer :: i, eoe
      integer :: iEdge, k
   
   end subroutine compute_output_diagnostics
   
!--------------------------------------------------------------------------
   
   subroutine mpas_timestep(domain, itimestep, dt, timeStamp)
   
      use mpas_grid_types
      use cice_time_integration
      use mpas_timer
   
      implicit none
   
      type (domain_type), intent(inout) :: domain 
      integer, intent(in) :: itimestep
      real (kind=RKIND), intent(in) :: dt
      character(len=*), intent(in) :: timeStamp

      type (block_type), pointer :: block

      block => domain % blocklist
      do while (associated(block))
         call cice_timestep(block, clock, dt)
         block => block % next
      end do

   end subroutine mpas_timestep
   
!--------------------------------------------------------------------------   

   subroutine mpas_core_finalize(domain)
   
      use mpas_grid_types
   
      implicit none

      type (domain_type), intent(inout) :: domain 
      integer :: ierr
 
     call mpas_destroy_clock(clock, ierr)

   end subroutine mpas_core_finalize

!--------------------------------------------------------------------------

   subroutine mpas_core_setup_packages(ierr)

      use mpas_packages

      implicit none

      integer, intent(out) :: ierr

      ierr = 0

   end subroutine mpas_core_setup_packages

!--------------------------------------------------------------------------

end module mpas_core
