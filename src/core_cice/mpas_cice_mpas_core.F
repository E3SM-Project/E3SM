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
   
      use mpas_grid_types
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      character(len=*), intent(out) :: startTimeStamp
   
      real (kind=RKIND), pointer :: dt
      type (block_type), pointer :: block

      type (MPAS_pool_type), pointer :: &
           configs

      configs => domain % blocklist % configs

      call MPAS_pool_get_config(configs, "config_dt", dt)

      call simulation_clock_init(domain, dt, startTimeStamp)

      block => domain % blocklist
      do while (associated(block))
         call mpas_init_block(block, dt)
         !block % state % time_levs(1) % state % xtime % scalar = startTimeStamp
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

      type (MPAS_pool_type), pointer :: &
           configs

      character(len=strKIND), pointer :: &
           config_start_time, &
           config_run_duration, &
           config_stop_time, &
           config_output_interval, &
           config_restart_interval

      configs => domain % blocklist % configs
      
      call MPAS_pool_get_config(configs, "config_start_time", config_start_time)
      call MPAS_pool_get_config(configs, "config_run_duration", config_run_duration)
      call MPAS_pool_get_config(configs, "config_stop_time", config_stop_time)
      call MPAS_pool_get_config(configs, "config_output_interval", config_output_interval)
      call MPAS_pool_get_config(configs, "config_restart_interval", config_restart_interval)

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

   subroutine mpas_init_block(block, dt)
   
      use mpas_grid_types
      use mpas_rbf_interpolation
      use mpas_vector_reconstruction
      
      use cice_initialize, only: cice_init

      implicit none

      type (block_type), intent(inout) :: block
      real (kind=RKIND), intent(in) :: dt

      type (MPAS_pool_type), pointer :: mesh

      integer, pointer :: nCells, vertexDegree

      call MPAS_pool_get_subpool(block % structs, "mesh", mesh)

      call mpas_rbf_interp_initialize(mesh)
      call mpas_init_reconstruct(mesh)

      ! initialize general stuff
      call cice_init(block, clock, dt)

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
      real (kind=RKIND), pointer :: dt
      type (block_type), pointer :: block_ptr

      type (MPAS_Time_Type) :: currTime
      character(len=StrKIND) :: timeStamp
      integer :: ierr

      type (MPAS_pool_type), pointer :: configs

      ! Eventually, dt should be domain specific
      configs => domain % blocklist % configs

      call MPAS_pool_get_config(configs, "config_dt", dt)

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
         !block_ptr => domain % blocklist
         !do while(associated(block_ptr))
         !   call mpas_shift_time_levels_state(block_ptr % state)
         !   block_ptr => block_ptr % next
         !end do

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

      type (MPAS_pool_type), pointer :: &
           configs

      integer, pointer :: config_frames_per_outfile

      configs => domain % blocklist % configs

      call MPAS_pool_get_config(configs, "config_frames_per_outfile", config_frames_per_outfile)
   
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

   subroutine mpas_core_setup_packages(configs, packages, ierr)

      use mpas_packages

      implicit none

      type (MPAS_pool_type), pointer :: &
           configs, &
           packages

      integer, intent(out) :: ierr

      logical, pointer :: &
           config_use_velocity_solver, &
           pkgWeakActive, &
           pkgVariationalActive, &
           pkgWachspressActive, &
           pkgPieceWiseLinearActive

      character(len=strKIND), pointer :: &
           config_stress_divergence_scheme, &
           config_variational_basis

      ierr = 0

      !pkgWeak
      !pkgVariational
      !pkgWachspress
      !pkgPieceWiseLinear

      call MPAS_pool_get_config(configs, "config_use_velocity_solver", config_use_velocity_solver)
      call MPAS_pool_get_config(configs, "config_stress_divergence_scheme", config_stress_divergence_scheme)
      call MPAS_pool_get_config(configs, "config_variational_basis", config_variational_basis)

      call MPAS_pool_get_package(packages, "pkgWeakActive", pkgWeakActive)
      call MPAS_pool_get_package(packages, "pkgVariationalActive", pkgVariationalActive)
      call MPAS_pool_get_package(packages, "pkgWachspressActive", pkgWachspressActive)
      call MPAS_pool_get_package(packages, "pkgPieceWiseLinearActive", pkgPieceWiseLinearActive)

      ! is the velocity solver on
      if (config_use_velocity_solver) then

         ! stress divergence scheme
         if (trim(config_stress_divergence_scheme) == "weak") then

            ! weak stress divergence scheme
            pkgWeakActive = .true.

         else if (trim(config_stress_divergence_scheme) == "variational") then

            ! variational stress divergence scheme
            pkgVariationalActive = .true.

            ! variational basis function scheme
            if (trim(config_variational_basis) == "wachspress") then

               ! wachspress basis functions
               pkgWachspressActive = .true.
               
            else if (trim(config_variational_basis) == "pwl") then

               ! piecewise linear basis functions
               pkgPieceWiseLinearActive = .true.

            else

               ! unknown value of config_variational_basis
               write(*,*) "Unknown variational scheme: ", trim(config_variational_basis)
               ierr = 1

            endif ! config_variational_basis

         else

            ! unknown value of config_stress_divergence_scheme
            write(*,*) "Unknown stress divergence scheme: ", trim(config_stress_divergence_scheme)
            ierr = 1

         endif ! config_stress_divergence_scheme

      endif ! config_use_velocity_solver

   end subroutine mpas_core_setup_packages

!--------------------------------------------------------------------------

end module mpas_core
