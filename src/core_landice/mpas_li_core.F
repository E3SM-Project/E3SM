! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module li_core

   use mpas_framework
   use mpas_timekeeping
   use mpas_abort, only : mpas_dmpar_global_abort

   implicit none
   private

   type (MPAS_Clock_type), pointer :: clock

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: li_core_init, &
             li_core_run, &
             li_core_finalize, &
             li_simulation_clock_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------


!***********************************************************************

   contains

!***********************************************************************



!***********************************************************************
!
!  function li_core_init
!
!> \brief   Initializes land ice core
!> \author  Matt Hoffman
!> \date    11 September 2013
!> \details
!>  This routine initializes the land ice core.
!
!-----------------------------------------------------------------------

   function li_core_init(domain, startTimeStamp) result(err)
   
      use mpas_derived_types
      use mpas_pool_routines
      use mpas_stream_manager
      use li_velocity
      use li_setup
!!!      use mpas_tracer_advection
!!!      use li_global_diagnostics
   
      implicit none
   
      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: domain  !< Input/output: Domain 

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------
      character(len=*), intent(out) :: startTimeStamp   !< Output: starting time stamp
   
      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      type (MPAS_Time_Type) :: startTime
      integer :: i, err, err_tmp, globalErr

      logical, pointer :: config_do_restart


      err = 0
      err_tmp = 0
      globalErr = 0

      call mpas_pool_get_config(domain % configs, 'config_do_restart', config_do_restart)
      !
      ! Initialize config option settings as needed
      !
      call li_setup_config_options( domain, err_tmp )
      err = ior(err, err_tmp)

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

      if (config_do_restart) then
         call mpas_stream_mgr_read(domain % streamManager, streamID='restart', ierr=err_tmp)
      else
         call mpas_stream_mgr_read(domain % streamManager, streamID='input', ierr=err_tmp)
      end if
      call MPAS_stream_mgr_reset_alarms(domain % streamManager, streamID='restart', ierr=err_tmp)
      err = ior(err, err_tmp)
      call MPAS_stream_mgr_reset_alarms(domain % streamManager, streamID='input', ierr=err_tmp)
      err = ior(err, err_tmp)
      call mpas_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=err_tmp)
      err = ior(err, err_tmp)


      ! ===
      ! === Initialize modules ===
      ! ===

      call li_velocity_init(domain, err_tmp)
      err = ior(err, err_tmp)

!!!      call mpas_tracer_advection_init(err_tmp)  ! Calling signature may be incorrect here.
!!!      err = ior(err,err_tmp)


      ! ===
      ! === Initialize blocks ===
      ! ===
      block => domain % blocklist
      do while (associated(block))
         call landice_init_block(block, startTimeStamp, domain % dminfo)

         block => block % next
      end do


      ! check for errors and exit
      call mpas_dmpar_max_int(domain % dminfo, err, globalErr)  ! Find out if any blocks got an error
      if (globalErr > 0) then
          call mpas_dmpar_global_abort("An error has occurred in li_core_init. Aborting...")
      endif

   !--------------------------------------------------------------------
   end function li_core_init



!***********************************************************************
!
!  function li_core_run
!
!> \brief   Runs the land ice core
!> \author  Matt Hoffman
!> \date    11 September 2013
!> \details
!>  This routine runs the land ice core.
!
!-----------------------------------------------------------------------

   function li_core_run(domain) result(err)
   
      use mpas_derived_types
      use mpas_pool_routines
      use mpas_kind_types
      use mpas_stream_manager
      use mpas_timer
      use li_diagnostic_vars
      use li_setup

      implicit none
   
      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: domain   !< Input/output: Domain 

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      integer :: itimestep
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool
      logical, pointer :: config_do_restart, config_write_output_on_startup
      character(len=StrKIND), pointer :: config_restart_timestamp_name

      type (MPAS_Time_Type) :: currTime
      character(len=StrKIND) :: timeStamp
      integer :: err, err_tmp, globalErr
      logical :: solveVelo
      type (MPAS_TimeInterval_type) :: timeStepInterval  !< time step as an interval
      real (kind=RKIND) :: dtSeconds  !< time step in seconds


      err = 0
      err_tmp = 0
      globalErr = 0

      ! Get Pool stuff that will be needed
      call mpas_pool_get_config(liConfigs, 'config_do_restart', config_do_restart)
      call mpas_pool_get_config(liConfigs, 'config_write_output_on_startup', config_write_output_on_startup)
      call mpas_pool_get_config(liConfigs, 'config_restart_timestamp_name', config_restart_timestamp_name)

      call mpas_timer_start("land ice core run")
      currTime = mpas_get_clock_time(clock, MPAS_NOW, err_tmp)
      err = ior(err, err_tmp)
      call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=err_tmp)
      err = ior(err, err_tmp)         
      write(0,*) 'Initial timestep ', trim(timeStamp)
      write(6,*) 'Initial timestep ', trim(timeStamp)


      ! ===
      ! === Calculate Initial state
      ! ===
      call mpas_timer_start("initial state calculation")

      ! On a restart, we already have the exact velocity field we need, 
      ! so don't do the expensive calculation again.
      if (config_do_restart) then
         solveVelo = .false.
      else
      ! Otherwise, we need to calculate velocity for the initial state
      !  (Note: even if the velocity is supplied, we should still calculate it
      !   to ensure it is consistent with the current geometry/B.C.  If the 
      !   velocity solver is iterative, the supplied field will be used as an 
      !   initial guess, so the solution should be quick.
         solveVelo = .true.
      endif

      call li_calculate_diagnostic_vars(domain, timeLevel=1, solveVelo=solveVelo, err=err_tmp)
      err = ior(err, err_tmp)

      call mpas_timer_stop("initial state calculation")


      ! ===
      ! === Write Initial Output
      ! ===
      call mpas_timer_start("write output")
      if (config_write_output_on_startup) then
         call mpas_stream_mgr_write(domain % streamManager, 'output', forceWriteNow=.true., ierr=err_tmp)
      endif
      call mpas_timer_stop("write output")

      ! === error check and exit
      call mpas_dmpar_max_int(domain % dminfo, err, globalErr)  ! Find out if any blocks got an error
      if (globalErr > 0) then
          call mpas_dmpar_global_abort("An error has occurred in li_core_run before time-stepping. Aborting...")
      endif



      ! During integration, time level 1 stores the model state at the beginning of the
      !   time step, and time level 2 stores the state advanced dt in time by timestep(...)
      itimestep = 0
      ! ===
      ! === Time step loop
      ! ===
      do while (.not. mpas_is_clock_stop_time(clock))

         itimestep = itimestep + 1

         ! Get the interval at this point in time - currently this does not change during the simulation, but re-calculating it explicitly for generality
         timeStepInterval = mpas_get_clock_timestep(clock, ierr=err_tmp)
         err = ior(err,err_tmp)
         ! Convert the clock's time interval into a dt in seconds to be used by the time stepper, using the currTime as the start time for this interval.
         ! (We want to do this conversion before advancing the clock because the dt in seconds may change
         !  as the base time changes, and we want the old time as the base time.
         !  For example, the number of seconds in a year will be longer in a leap year.)
         call mpas_get_timeInterval(timeStepInterval, StartTimeIn=currTime, dt=dtSeconds, ierr=err_tmp)
         err = ior(err,err_tmp)

         call mpas_advance_clock(clock)

         currTime = mpas_get_clock_time(clock, MPAS_NOW, err_tmp)
         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=err_tmp)
         err = ior(err, err_tmp)
         write(0,*) 'Doing timestep ', trim(timeStamp)
         write(6,*) 'Doing timestep ', trim(timeStamp)

         !write(6,*) '  dt (s) = ', dtSeconds


         ! ===
         ! === Perform Timestep
         ! ===
         call mpas_timer_start("time integration")

         call landice_timestep(domain, itimestep, dtSeconds, timeStamp, err_tmp)
         err = ior(err,err_tmp)

         ! Move time level 2 fields back into time level 1 for next time step
         block => domain % blocklist
         do while(associated(block))
            call mpas_pool_get_subpool(block % structs, 'state', statePool)
            call mpas_pool_shift_time_levels(statePool)
            block => block % next
         end do
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
         call mpas_stream_mgr_read(domain % streamManager, whence=MPAS_STREAM_LATEST_BEFORE, ierr=err_tmp)
         err = ior(err, err_tmp)
         call mpas_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_INPUT, ierr=err_tmp)
         err = ior(err, err_tmp)


         ! ===
         ! === Write Output and/or Restart, if needed
         ! ===
         call mpas_timer_start("write output")
         ! Update the restart_timestamp file with the new time, if needed.
         if ( mpas_stream_mgr_ringing_alarms(domain % streamManager, streamID='restart', direction=MPAS_STREAM_OUTPUT, ierr=err_tmp) ) then
            open(22, file=config_restart_timestamp_name, form='formatted', status='replace')
            write(22, *) timeStamp
            close(22)
         end if
         err = ior(err, err_tmp)
         ! These calls will handle ALL output streams that need writing.
         ! [Could add them individually, as the ocean does, if some other actions are needed when a
         ! specific alarm is ringing (e.g., global stats calculated only when output stream gets written)]
         call mpas_stream_mgr_write(domain % streamManager, ierr=err_tmp)
         err = ior(err, err_tmp)
         call mpas_stream_mgr_reset_alarms(domain % streamManager,   direction=MPAS_STREAM_OUTPUT, ierr=err_tmp)
         err = ior(err, err_tmp)
         call mpas_timer_stop("write output")


         ! === error check and exit
         call mpas_dmpar_max_int(domain % dminfo, err, globalErr)  ! Find out if any blocks got an error
         if (globalErr > 0) then
             call mpas_dmpar_global_abort("An error has occurred in li_core_run. Aborting...")
         endif

      end do
      call mpas_timer_stop("land ice core run")

   !--------------------------------------------------------------------
   end function li_core_run



!***********************************************************************
!
!  function li_core_finalize
!
!> \brief   Finalizes the land ice core
!> \author  Matt Hoffman
!> \date    11 September 2013
!> \details
!>  This routine finalizes the land ice core.
!
!-----------------------------------------------------------------------
   function li_core_finalize(domain) result(err)
   
      use mpas_derived_types
      use mpas_decomp
      use li_velocity, only: li_velocity_finalize
   
      implicit none

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: domain    !< Input/output: Domain 

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      integer :: err, err_tmp, globalErr
 
      err = 0
      err_tmp = 0
      globalErr = 0

      call li_velocity_finalize(domain, err_tmp)
      err = ior(err, err_tmp)

      call mpas_destroy_clock(clock, err_tmp)
      err = ior(err, err_tmp)

      call mpas_decomp_destroy_decomp_list(domain % decompositions)

      ! === error check and exit
      call mpas_dmpar_max_int(domain % dminfo, err, globalErr)  ! Find out if any blocks got an error
      if (globalErr > 0) then
          call mpas_dmpar_global_abort("An error has occurred in li_core_finalize. Aborting...")
      endif

   !--------------------------------------------------------------------
   end function li_core_finalize


!***********************************************************************
!***********************************************************************
! Private subroutines:
!***********************************************************************
!***********************************************************************



!***********************************************************************
!
!  routine landice_init_block
!
!> \brief   Initializes blocks for the land ice core
!> \author  Matt Hoffman
!> \date    11 September 2013
!> \details
!>  This routine initializes blocks for the land ice core.
!
!-----------------------------------------------------------------------
   subroutine landice_init_block(block, startTimeStamp, dminfo)

      use mpas_derived_types
      use mpas_pool_routines
      use mpas_rbf_interpolation
      use mpas_vector_reconstruction
      use li_setup
      use li_mask
      use li_velocity
   
      implicit none
   
      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      type (dm_info), intent(in) :: dminfo             !< Input: Domain info
      character(len=*), intent(in) :: startTimeStamp   !< Input: time stamp at start 

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      type (block_type), intent(inout) :: block         !< Input/output: Block object

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: statePool
      character (len=StrKIND), pointer :: xtime
      type (MPAS_Time_Type) :: currTime
      integer :: err, err_tmp
      integer :: iCell, iLevel, i


      err = 0
      err_tmp = 0

      ! Get pool stuff
      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

      ! ===
      ! === Call init routines ===
      ! ===
      call li_setup_vertical_grid(meshPool, err_tmp)
      err = ior(err, err_tmp)

      call li_setup_sign_and_index_fields(meshPool)


! This was needed to init FCT once.
!!!      ! Init for FCT tracer advection
!!!      mesh % maxLevelCell % array = mesh % nVertLevels ! Needed for FCT tracer advection
!!!      mesh % maxLevelEdgeTop % array = mesh % nVertLevels ! Needed for FCT tracer advection
!!!      mesh % maxLevelEdgeBot % array = mesh % nVertLevels ! Needed for FCT tracer advection
!!!      call ocn_initialize_advection_rk(mesh, err)
!!!      call mpas_ocn_tracer_advection_coefficients(mesh, err_tmp)
!!!      err = ior(err, err_tmp)

      call mpas_timer_start("initialize velocity")
      call li_velocity_block_init(block, err_tmp)
      err = ior(err, err_tmp)
      call mpas_timer_stop("initialize velocity")

      ! Init for reconstruction of velocity
      call mpas_rbf_interp_initialize(meshPool)
      call mpas_init_reconstruct(meshPool)

      ! Assign initial time stamp
      call mpas_pool_get_array(statePool, 'xtime', xtime, timeLevel=1)
      xtime = startTimeStamp

      ! Mask init identifies initial ice extent
      call li_calculate_mask_init(meshPool, statePool, timeLevel=1, err=err_tmp)
      err = ior(err, err_tmp)

      ! Make sure all time levels have a copy of the initial state
      call mpas_pool_initialize_time_levels(statePool)

      ! === error check
      if (err > 0) then
          write (0,*) "An error has occurred in init_block."
      endif

   !--------------------------------------------------------------------
   end subroutine landice_init_block



!***********************************************************************
!
!  routine landice_timestep
!
!> \brief   Performs a time step
!> \author  Matt Hoffman
!> \date    11 September 2013
!> \details
!>  This routine performs a time step for the land ice core.
!
!-----------------------------------------------------------------------
   subroutine landice_timestep(domain, itimestep, dt, timeStamp, err)
   
      use mpas_derived_types
      use li_time_integration
      use mpas_timer
!!!      use li_global_diagnostics
   
      implicit none
   
      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      integer, intent(in) :: itimestep              !< Input: time step number 
      real (kind=RKIND), intent(in) :: dt           !< Input: time step, in seconds
      character(len=*), intent(in) :: timeStamp     !< Input: time stamp of current time step

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: domain   !< Input/output: Domain 

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      type (block_type), pointer :: block_ptr
      integer :: err_tmp

   
      err = 0
      err_tmp = 0


      call li_timestep(domain, dt, timeStamp, err_tmp)
      err = ior(err,err_tmp)
   
!!!      if (config_stats_interval .gt. 0) then
!!!          if(mod(itimestep, config_stats_interval) == 0) then
!!!              block_ptr => domain % blocklist
!!!              if(associated(block_ptr % next)) then
!!!                  write(0,*) 'Error: computeGlobalDiagnostics assumes ',&
!!!                             'that there is only one block per processor.'
!!!              end if
!!!   
!!!              call mpas_timer_start("global_diagnostics")
!!!              call li_compute_global_diagnostics(domain % dminfo, &
!!!                       block_ptr % state % time_levs(2) % state, block_ptr % mesh, &
!!!                       itimestep, dt)
!!!              call mpas_timer_stop("global_diagnostics")
!!!          end if
!!!      end if

      !TODO: replace the above code block with this if we desire to convert config_stats_interval to use alarms
      !if (mpas_is_alarm_ringing(clock, statsAlarmID, ierr=ierr)) then
      !   call mpas_reset_clock_alarm(clock, statsAlarmID, ierr=ierr)

      !   block_ptr => domain % blocklist
      !   if(associated(block_ptr % next)) then
      !      write(0,*) 'Error: computeGlobalDiagnostics assumes ',&
      !                 'that there is only one block per processor.'
      !   end if

      !   call mpas_timer_start("global_diagnostics")
      !   call sw_compute_global_diagnostics(domain % dminfo, &
      !            block_ptr % state % time_levs(2) % state, block_ptr % mesh, &
      !            timeStamp, dt)
      !   call mpas_timer_stop("global_diagnostics")
      !end if

      ! === error check
      if (err > 0) then
          write (0,*) "An error has occurred in mpas_timestep."
      endif

   
   end subroutine landice_timestep
   
   

!***********************************************************************
!
!  routine li_simulation_clock_init
!
!> \brief   Initializes the simulation clock
!> \author  ??
!> \date    ??
!> \details
!>  This routine initializes the simulation clock.
!
!-----------------------------------------------------------------------

   subroutine li_simulation_clock_init(core_clock, configs, ierr)

      use mpas_derived_types
      use mpas_pool_routines

      implicit none

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      type (MPAS_Clock_type), intent(inout) :: core_clock  !< Input/output: core_clock
      type (mpas_pool_type), intent(inout) :: configs      !< Input/output: configs

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------
      integer, intent(out) :: ierr !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      type (MPAS_Time_Type) :: startTime, stopTime, alarmStartTime
      type (MPAS_TimeInterval_type) :: runDuration, timeStep, alarmTimeStep
      character (len=StrKIND), pointer :: config_start_time, config_run_duration, config_stop_time, config_output_interval, config_restart_interval ! MPAS standard configs
      character (len=StrKIND), pointer :: config_dt  ! MPAS LI-specific config option
      character (len=StrKIND), pointer :: config_restart_timestamp_name
      character (len=StrKIND) :: restartTimeStamp !< string to be read from file
      integer :: err_tmp


      ierr = 0
      err_tmp = 0

      call mpas_pool_get_config(configs, 'config_dt', config_dt)
      call mpas_pool_get_config(configs, 'config_start_time', config_start_time)
      call mpas_pool_get_config(configs, 'config_run_duration', config_run_duration)
      call mpas_pool_get_config(configs, 'config_stop_time', config_stop_time)
      call mpas_pool_get_config(configs, 'config_output_interval', config_output_interval)
      call mpas_pool_get_config(configs, 'config_restart_interval', config_restart_interval)
      call mpas_pool_get_config(configs, 'config_restart_timestamp_name', config_restart_timestamp_name)


      ! Set time to the user-specified start time OR use a restart time from file
      if ( trim(config_start_time) == "file" ) then
         open(22, file=config_restart_timestamp_name, form='formatted', status='old')
         read(22,*) restartTimeStamp
         close(22)
         call mpas_set_time(curr_time=startTime, dateTimeString=restartTimeStamp, ierr=err_tmp)
         ierr = ior(ierr,err_tmp)
      else
         call mpas_set_time(curr_time=startTime, dateTimeString=config_start_time, ierr=err_tmp)
         ierr = ior(ierr,err_tmp)
      end if

      ! Set interval to the user-specified time interval string
      call mpas_set_timeInterval(timeStep, timeString=config_dt, ierr=err_tmp)
      ierr = ior(ierr,err_tmp)


      if (trim(config_run_duration) /= "none") then
         call mpas_set_timeInterval(runDuration, timeString=config_run_duration, ierr=err_tmp)
         ierr = ior(ierr,err_tmp)
         call mpas_create_clock(core_clock, startTime=startTime, timeStep=timeStep, runDuration=runDuration, ierr=err_tmp)
         ierr = ior(ierr,err_tmp)

         if (trim(config_stop_time) /= "none") then
            call mpas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, ierr=err_tmp)
            ierr = ior(ierr,err_tmp)
            if(startTime + runduration /= stopTime) then
               write(0,*) 'Warning: config_run_duration and config_stop_time are inconsistent: using config_run_duration.'
            end if
         end if
      else if (trim(config_stop_time) /= "none") then
         call mpas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, ierr=err_tmp)
         ierr = ior(ierr,err_tmp)
         call mpas_create_clock(core_clock, startTime=startTime, timeStep=timeStep, stopTime=stopTime, ierr=err_tmp)
         ierr = ior(ierr,err_tmp)
      else
          write(0,*) 'Error: Neither config_run_duration nor config_stop_time were specified.'
          ierr = 1
      end if


      ! === error check
      if (ierr > 0) then
          write (0,*) "An error has occurred in li_simulation_clock_init."
      endif

   !--------------------------------------------------------------------
   end subroutine li_simulation_clock_init



end module li_core
