! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_core

   use mpas_framework
   use mpas_timekeeping

   implicit none
   private

   type (io_output_object), save :: restart_obj
   integer :: current_outfile_frames

   type (MPAS_Clock_type) :: clock

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: mpas_core_init, &
             mpas_core_run, &
             mpas_core_finalize

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   integer, parameter :: outputAlarmID = 1
   integer, parameter :: restartAlarmID = 2
   !integer, parameter :: statsAlarmID = 3

!***********************************************************************

   contains

!***********************************************************************



!***********************************************************************
!
!  routine mpas_core_init
!
!> \brief   Initializes land ice core
!> \author  Matt Hoffman
!> \date    11 September 2013
!> \version SVN:$Id$
!> \details
!>  This routine initializes the land ice core.
!
!-----------------------------------------------------------------------

   subroutine mpas_core_init(domain, startTimeStamp)
   
      use mpas_configure
      use mpas_grid_types
      use li_velocity
      use li_setup
!!!      use mpas_ocn_tracer_advection
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

      integer :: i, err, err_tmp, globalErr



      err = 0
      err_tmp = 0
      globalErr = 0

      !
      ! Initialize core
      !

      call li_setup_timestep(err_tmp)
      err = ior(err, err_tmp)

      ! ===
      ! === Initialize modules ===
      ! ===

      call simulation_clock_init(domain, li_deltat, startTimeStamp, err_tmp)
      err = ior(err, err_tmp)

      call mpas_timer_init(domain)

      call li_velocity_init(domain, err_tmp)
      err = ior(err, err_tmp)

!!!      call mpas_ocn_tracer_advection_init(err_tmp)
!!!      err = ior(err,err_tmp)

      ! ===
      ! === Initialize blocks ===
      ! ===
      block => domain % blocklist
      do while (associated(block))
         call mpas_init_block(block, block % mesh, startTimeStamp, domain % dminfo)

         block => block % next
      end do

      current_outfile_frames = 0

      ! check for errors and exit
      call mpas_dmpar_max_int(domain % dminfo, err, globalErr)  ! Find out if any blocks got an error
      if (globalErr > 0) then
          call mpas_dmpar_global_abort("An error has occurred in mpas_core_init. Aborting...")
      endif

   !--------------------------------------------------------------------
   end subroutine mpas_core_init



!***********************************************************************
!
!  routine mpas_core_run
!
!> \brief   Runs the land ice core
!> \author  Matt Hoffman
!> \date    11 September 2013
!> \version SVN:$Id$
!> \details
!>  This routine runs the land ice core.
!
!-----------------------------------------------------------------------

   subroutine mpas_core_run(domain, output_obj, output_frame)
   
      use mpas_grid_types
      use mpas_kind_types
      use mpas_io_output
      use mpas_timer
      use li_diagnostic_vars
      use li_setup
!!!      use li_annual_forcing
   
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
      type (io_output_object), intent(inout) :: output_obj   !< Input/output: i/o output object
      integer, intent(inout) :: output_frame   !< Input/output: output frame number

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
      type (block_type), pointer :: block_ptr
      type (MPAS_Time_Type) :: currTime
      character(len=StrKIND) :: timeStamp
      integer :: err, err_tmp, globalErr
      logical :: solveVelo
      ! Note: li_deltat is a module variable in li_setup   


      err = 0
      err_tmp = 0
      globalErr = 0


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
      call write_output_frame(output_obj, output_frame, domain)

!!! These lines were used previously instead of the one above - but we may need to reassess how we want output to work.
!!!      if (config_write_output_on_startup) then
!!!         call mpas_timer_start("write output frame")
!!!         call write_output_frame(output_obj, output_frame, domain)
!!!         call mpas_timer_stop("write output frame")
!!!      elseif (config_frames_per_outfile /= 0) then
!!!         ! if we are not writing on startup then close the already open file since it will have a confusing time stamp (the time stamp would be the inital time, but the intial time is not present in the file!)
!!!         call mpas_output_state_finalize(output_obj, domain % dminfo)
!!!      endif


      ! === error check and exit
      call mpas_dmpar_max_int(domain % dminfo, err, globalErr)  ! Find out if any blocks got an error
      if (globalErr > 0) then
          call mpas_dmpar_global_abort("An error has occurred in mpas_core_run before time-stepping. Aborting...")
      endif


      ! During integration, time level 1 stores the model state at the beginning of the
      !   time step, and time level 2 stores the state advanced dt in time by timestep(...)
      itimestep = 0
      ! ===
      ! === Time step loop
      ! ===
      do while (.not. mpas_is_clock_stop_time(clock))

         itimestep = itimestep + 1
         call mpas_advance_clock(clock)

         currTime = mpas_get_clock_time(clock, MPAS_NOW, err_tmp)
         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=err_tmp)         
         err = ior(err, err_tmp)
         write(0,*) 'Doing timestep ', trim(timeStamp)
         write(6,*) 'Doing timestep ', trim(timeStamp)

!!! These lines were used previously to assign annual forcing fields.  We may want to reassess how to do that.
!!!         call mpas_timer_start("assign forcing fields")
!!!         block_ptr => domain % blocklist
!!!         do while(associated(block_ptr))
!!!           call li_assign_annual_forcing(currTime, block_ptr % mesh, ierr)
!!!           block_ptr => block_ptr % next
!!!         end do
!!!         call mpas_timer_stop("assign forcing fields")


         ! ===
         ! === Perform Timestep
         ! ===
         call mpas_timer_start("time integration")
         call mpas_timestep(domain, itimestep, li_deltat, timeStamp, err_tmp)
         err = ior(err,err_tmp)

         ! Move time level 2 fields back into time level 1 for next time step
         block_ptr => domain % blocklist
         do while(associated(block_ptr))
            call mpas_shift_time_levels_state(block_ptr % state)
            block_ptr => block_ptr % next
         end do
         call mpas_timer_stop("time integration")


         !TODO: mpas_get_clock_ringing_alarms is probably faster than multiple mpas_is_alarm_ringing...

         ! ===
         ! === Write Output, if needed
         ! ===
         if (mpas_is_alarm_ringing(clock, outputAlarmID, ierr=err_tmp)) then
            call mpas_timer_start("write output")
            call mpas_reset_clock_alarm(clock, outputAlarmID, ierr=err_tmp)
            ! output_frame will always be > 1 here unless it was reset after the maximum number of frames per outfile was reached
            ! MJH: See old code for adjustments made here to deal with some obscure, specific use-cases
            if(output_frame == 1) then
               call mpas_output_state_finalize(output_obj, domain % dminfo)
               call mpas_output_state_init(output_obj, domain, "OUTPUT", trim(timeStamp))
            end if
            call write_output_frame(output_obj, output_frame, domain)
            call mpas_timer_stop("write output")
         end if
         err = ior(err, err_tmp)

         ! ===
         ! === Write Restart, if needed
         ! ===
         if (mpas_is_alarm_ringing(clock, restartAlarmID, ierr=err_tmp)) then
            call mpas_timer_start("write restart")
            call mpas_reset_clock_alarm(clock, restartAlarmID, ierr=err_tmp)

            ! Write one restart time per file
            call mpas_output_state_init(restart_obj, domain, "RESTART", trim(timeStamp))
            call mpas_output_state_for_domain(restart_obj, domain, 1)
            call mpas_output_state_finalize(restart_obj, domain % dminfo)
            call mpas_timer_stop("write restart")
         end if
         err = ior(err, err_tmp)

         ! === error check and exit
         call mpas_dmpar_max_int(domain % dminfo, err, globalErr)  ! Find out if any blocks got an error
         if (globalErr > 0) then
             call mpas_dmpar_global_abort("An error has occurred in mpas_core_run. Aborting...")
         endif

      end do
      call mpas_timer_stop("land ice core run")



   !--------------------------------------------------------------------
   end subroutine mpas_core_run



!***********************************************************************
!
!  routine mpas_core_finalize
!
!> \brief   Finalizes the land ice core
!> \author  Matt Hoffman
!> \date    11 September 2013
!> \version SVN:$Id$
!> \details
!>  This routine finalizes the land ice core.
!
!-----------------------------------------------------------------------
   subroutine mpas_core_finalize(domain)
   
      use mpas_grid_types
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

      ! === error check and exit
      call mpas_dmpar_max_int(domain % dminfo, err, globalErr)  ! Find out if any blocks got an error
      if (globalErr > 0) then
          call mpas_dmpar_global_abort("An error has occurred in mpas_core_finalize. Aborting...")
      endif

   !--------------------------------------------------------------------
   end subroutine mpas_core_finalize


!***********************************************************************
!***********************************************************************
! Private subroutines:
!***********************************************************************
!***********************************************************************



!***********************************************************************
!
!  routine mpas_init_block
!
!> \brief   Initializes blocks for the land ice core
!> \author  Matt Hoffman
!> \date    11 September 2013
!> \version SVN:$Id$
!> \details
!>  This routine initializes blocks for the land ice core.
!
!-----------------------------------------------------------------------
   subroutine mpas_init_block(block, mesh, startTimeStamp, dminfo)
   
      use mpas_grid_types
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
      type (dm_info), intent(in) :: dminfo    !< Input: Domain info
      character(len=*), intent(in) :: startTimeStamp   !< Input: time stamp at start 

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      type (block_type), intent(inout) :: block   !< Input/output: Block object
      type (mesh_type), intent(inout) :: mesh     !< Input/output: Mesh object

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
      type (MPAS_Time_Type) :: currTime
      integer :: err, err_tmp

      type (state_type), pointer :: state  

      integer :: iCell, iLevel, i


      err = 0
      err_tmp = 0

      state => block % state % time_levs(1) % state  ! initial state

      ! Make sure all time levels have a copy of the initial state
      do i=2,nTimeLevs
         call mpas_copy_state(block % state % time_levs(i) % state,   state)
      end do
   
      ! Assign initial time stamp
      state % xtime % scalar = startTimeStamp

      ! ===
      ! === Call init routines ===
      ! ===
      call li_setup_vertical_grid(mesh, err_tmp)
      err = ior(err, err_tmp)

      call li_setup_sign_and_index_fields(mesh)

! This was needed for annual forcing, which may or not be implemented in the same way.
!!!      currTime = mpas_get_clock_time(clock, MPAS_NOW, err_tmp)
!!!      err = ior(err, err_tmp)
!!!      call li_assign_annual_forcing(currTime, mesh, err_tmp)
!!!      err = ior(err, err_tmp)

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
      call mpas_rbf_interp_initialize(mesh)
      call mpas_init_reconstruct(mesh)

      ! Mask init identifies initial ice extent
      call li_calculate_mask_init(mesh, state, err)
      err = ior(err, err_tmp)


      ! === error check
      if (err > 0) then
          write (0,*) "An error has occurred in mpas_init_block."
      endif

   !--------------------------------------------------------------------
   end subroutine mpas_init_block



!***********************************************************************
!
!  routine mpas_timestep
!
!> \brief   Performs a time step
!> \author  Matt Hoffman
!> \date    11 September 2013
!> \version SVN:$Id$
!> \details
!>  This routine performs a time step for the land ice core.
!
!-----------------------------------------------------------------------
   subroutine mpas_timestep(domain, itimestep, dt, timeStamp, err)
   
      use mpas_grid_types
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
      real (kind=RKIND), intent(in) :: dt           !< Input: time step 
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

   
   end subroutine mpas_timestep
   
   

!***********************************************************************
!
!  routine simulation_clock_init
!
!> \brief   Initializes the simulation clock
!> \author  ??
!> \date    ??
!> \version SVN:$Id$
!> \details
!>  This routine initializes the simulation clock.
!
!-----------------------------------------------------------------------

   subroutine simulation_clock_init(domain, dt, startTimeStamp, ierr)

      implicit none

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      real (kind=RKIND), intent(in) :: dt   !< Input: time step

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
      character(len=*), intent(out) :: startTimeStamp   !< Output: time stamp at start 
      integer, intent(out) :: ierr !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      type (MPAS_Time_Type) :: startTime, stopTime, alarmStartTime
      type (MPAS_TimeInterval_type) :: runDuration, timeStep, alarmTimeStep
      integer :: err_tmp


      ierr = 0
      err_tmp = 0

      call mpas_set_time(curr_time=startTime, dateTimeString=config_start_time, ierr=err_tmp)
      ierr = ior(ierr,err_tmp)
      call mpas_set_timeInterval(timeStep, dt=dt, ierr=err_tmp)
      ierr = ior(ierr,err_tmp)

      if (trim(config_run_duration) /= "none") then
         call mpas_set_timeInterval(runDuration, timeString=config_run_duration, ierr=err_tmp)
         ierr = ior(ierr,err_tmp)
         call mpas_create_clock(clock, startTime=startTime, timeStep=timeStep, runDuration=runDuration, ierr=err_tmp)
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
         call mpas_create_clock(clock, startTime=startTime, timeStep=timeStep, stopTime=stopTime, ierr=err_tmp)
         ierr = ior(ierr,err_tmp)
      else
          write(0,*) 'Error: Neither config_run_duration nor config_stop_time were specified.'
          call mpas_dmpar_abort(domain % dminfo)
      end if

      ! set output alarm
      call mpas_set_timeInterval(alarmTimeStep, timeString=config_output_interval, ierr=err_tmp)
      ierr = ior(ierr,err_tmp)
      alarmStartTime = startTime + alarmTimeStep
      call mpas_add_clock_alarm(clock, outputAlarmID, alarmStartTime, alarmTimeStep, ierr=err_tmp)
      ierr = ior(ierr,err_tmp)

      ! set restart alarm, if necessary
      if (trim(config_restart_interval) /= "none") then
         call mpas_set_timeInterval(alarmTimeStep, timeString=config_restart_interval, ierr=err_tmp)
         ierr = ior(ierr,err_tmp)
         alarmStartTime = startTime + alarmTimeStep
         call mpas_add_clock_alarm(clock, restartAlarmID, alarmStartTime, alarmTimeStep, ierr=err_tmp)
         ierr = ior(ierr,err_tmp)
      end if

      !TODO: use this code if we desire to convert config_stats_interval to alarms 
      !(must also change config_stats_interval type to character) 
      ! set stats alarm, if necessary
      !if (trim(config_stats_interval) /= "none") then      
      !   call mpas_set_timeInterval(alarmTimeStep, timeString=config_stats_interval, ierr=err_tmp)
      !   alarmStartTime = startTime + alarmTimeStep
      !   call mpas_add_clock_alarm(clock, statsAlarmID, alarmStartTime, alarmTimeStep, ierr=err_tmp)
      !end if

      call mpas_get_time(curr_time=startTime, dateTimeString=startTimeStamp, ierr=err_tmp)
      ierr = ior(ierr,err_tmp)

      ! === error check
      if (ierr > 0) then
          write (0,*) "An error has occurred in simulation_clock_init."
      endif

   !--------------------------------------------------------------------
   end subroutine simulation_clock_init
   

   
!***********************************************************************
!
!  routine write_output_frame
!
!> \brief   Writes an output frame
!> \author  ??
!> \date    ??
!> \version SVN:$Id$
!> \details
!>  This routine writes an output frame.
!
!-----------------------------------------------------------------------
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
! MJH: My old version did things a bit differently here.
            current_outfile_frames = 0
            output_frame = 1
         end if
      end if

   end subroutine write_output_frame

   
   
!***********************************************************************
!
!  routine compute_output_diagnostics
!
!> \brief   Computes diagnostic fields only needed for output
!> \author  ??
!> \date    ??
!> \version SVN:$Id$
!> \details
!>  This routine is a placeholder for the calculation of any
!>  fields that are only needed for diagnostic output.
!
!-----------------------------------------------------------------------
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
   
   

end module mpas_core
