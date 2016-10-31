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
   use li_analysis_driver

   implicit none
   private


   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: li_core_init, &
             li_core_run, &
             li_core_finalize, &
             li_simulation_clock_init, &
             li_core_initial_solve

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
      use li_thermal
      use li_setup
      use li_constants
      use li_subglacial_hydro
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
      real (kind=RKIND), pointer :: deltat  ! variable in each block
      real (kind=RKIND) :: dtSeconds ! local variable
      type (MPAS_Pool_type), pointer :: meshPool
      type (MPAS_TimeInterval_type) :: timeStepInterval
      character (len=StrKIND), pointer :: xtime, simulationStartTime
      real (kind=RKIND), pointer :: daysSinceStart
      integer, dimension(:), pointer :: cellProcID
      type (MPAS_Time_type) :: xtime_timeType, simulationStartTime_timeType


      err = 0
      err_tmp = 0
      globalErr = 0

      call mpas_pool_get_config(domain % configs, 'config_do_restart', config_do_restart)

      !
      ! Initialize config option settings as needed
      !
      call li_setup_config_options( domain, err_tmp )
      err = ior(err, err_tmp)

      if (config_do_restart) then
         call mpas_stream_mgr_read(domain % streamManager, streamID='restart', ierr=err_tmp)
         err = ior(err, err_tmp)
      else
         call mpas_stream_mgr_read(domain % streamManager, streamID='input', ierr=err_tmp)
         err = ior(err, err_tmp)
      end if
      call MPAS_stream_mgr_reset_alarms(domain % streamManager, streamID='restart', ierr=err_tmp)
      err = ior(err, err_tmp)
      call MPAS_stream_mgr_reset_alarms(domain % streamManager, streamID='input', ierr=err_tmp)
      err = ior(err, err_tmp)
      call mpas_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_OUTPUT, ierr=err_tmp)
      err = ior(err, err_tmp)

      !
      ! Read the remaining input streams
      !
      call mpas_timer_start('io_read', .false.)
      call mpas_stream_mgr_read(domain % streamManager, ierr=err_tmp)
      err = ior(err, err_tmp)
      call mpas_timer_stop('io_read')
      call mpas_timer_start('reset_io_alarms', .false.)
      call mpas_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_INPUT, ierr=err_tmp)
      err = ior(err, err_tmp)
      call mpas_timer_stop('reset_io_alarms')

      ! ===
      ! Initialize some time stuff on each block
      ! ===
      ! Set startTimeStamp based on the start time of the simulation clock
      startTime = mpas_get_clock_time(domain % clock, MPAS_START_TIME, err_tmp)
      call mpas_get_time(startTime, dateTimeString=startTimeStamp) ! Get the start time as a time stamp
      err = ior(err, err_tmp)

      timeStepInterval = mpas_get_clock_timestep(domain % clock, ierr=err_tmp)  ! get timestep interval object
      err = ior(err,err_tmp)
      call mpas_get_timeInterval(timeStepInterval, StartTimeIn=startTime, dt=dtSeconds, ierr=err_tmp)  ! Get config dt in seconds
      err = ior(err,err_tmp)

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

         ! Assign initial time stamp
         call mpas_pool_get_array(meshPool, 'xtime', xtime)
         xtime = startTimeStamp

         ! Set simulationStartTime on a cold start only - it should have been read from restart otherwise
         call mpas_pool_get_array(meshPool, 'simulationStartTime', simulationStartTime)
         if (trim(simulationStartTime)=="no_date_available") then
            if (config_do_restart) then
               write(stderrUnit,*) "Warning: simulationStartTime is undefined but this is a restart, so setting its value to the current start time which may be incorrect."
            endif
            simulationStartTime = startTimeStamp
         elseif ((trim(simulationStartTime)/="no_date_available") .and. (.not. config_do_restart)) then
            write(stderrUnit,*) "Warning: This is a cold start but simulationStartTime is already defined, so its value may be incorrect."
         endif

         ! compute time since start of simulation, in days
         call mpas_pool_get_array(meshPool, 'daysSinceStart',daysSinceStart)
         call mpas_set_time(xtime_timeType, dateTimeString=xtime)
         call mpas_set_time(simulationStartTime_timeType, dateTimeString=simulationStartTime)
         call mpas_get_timeInterval(xtime_timeType - simulationStartTime_timeType,dt=daysSinceStart)
         daysSinceStart = daysSinceStart / seconds_per_day

         ! Initialize dt in seconds
         call mpas_pool_get_array(meshPool, 'deltat', deltat)
         deltat = dtSeconds

         block => block % next
      end do

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
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_array(meshPool, 'cellProcID', cellProcID)
         cellProcID(:) = domain % dminfo % my_proc_id

         call landice_init_block(block, domain % dminfo, err_tmp)
         err = ior(err, err_tmp)

         block => block % next
      end do

      ! initialize thermal solver
      ! Note: This subroutine includes a loop over blocks
      call li_thermal_init(domain, err_tmp)
      err = ior(err, err_tmp)

      ! Initialize subglacial hydrology solver
      call li_SGH_init(domain, err_tmp)
      err = ior(err, err_tmp)

      ! initialize analysis driver
      call li_analysis_init(domain, err_tmp)
      err = ior(err, err_tmp)

      ! halo update for reconstruction coefficients      
      ! Note: Results on multiple processors may be incorrect without this update

      call mpas_timer_start("halo updates")
      call mpas_dmpar_field_halo_exch(domain, 'coeffs_reconstruct')
      call mpas_timer_stop("halo updates")

      ! check for errors and exit

      call mpas_dmpar_max_int(domain % dminfo, err, globalErr)  ! Find out if any blocks got an error
      if (globalErr > 0) then
          call mpas_dmpar_global_abort("MPAS-landice: An error has occurred in li_core_init. Aborting...")
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
      use li_statistics
      use li_time_integration
      use mpas_io_streams, only: MPAS_STREAM_LATEST_BEFORE

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
      type (mpas_pool_type), pointer :: geometryPool
      integer, pointer :: config_stats_interval   !< interval (number of timesteps) for writing stats
      character(len=StrKIND), pointer :: config_restart_timestamp_name
      character(len=StrKIND), pointer :: config_velocity_solver
      ! Variables needed for printing timestamps
      type (MPAS_Time_Type) :: currTime
      character(len=StrKIND) :: timeStamp

      integer :: err, err_tmp, globalErr


      err = 0
      err_tmp = 0
      globalErr = 0

      ! Get Pool stuff that will be needed
      call mpas_pool_get_config(liConfigs, 'config_restart_timestamp_name', config_restart_timestamp_name)
      call mpas_pool_get_config(liConfigs, 'config_stats_interval', config_stats_interval)
      call mpas_pool_get_config(liConfigs, 'config_velocity_solver', config_velocity_solver)

      call mpas_timer_start("land ice core run")

      ! ===
      ! Solve initial state before beginning time stepping
      ! ===
      err_tmp = li_core_initial_solve(domain)
      err = ior(err,err_tmp)  ! li_core_finalize would abort if there was an error, but being safe.


      ! During integration, time level 1 stores the model state at the beginning of the
      !   time step, and time level 2 stores the state advanced dt in time by timestep(...)
      itimestep = 0
      ! ===
      ! === Time step loop
      ! ===
      do while (.not. mpas_is_clock_stop_time(domain % clock))

         itimestep = itimestep + 1
         write(stderrUnit,*) 'Starting timestep number ', iTimeStep
         write(stdoutUnit,*) 'Starting timestep number ', iTimeStep


         ! ===
         ! === Perform Timestep
         ! ===
         call mpas_timer_start("time integration")

         call li_timestep(domain, err_tmp)
         err = ior(err,err_tmp)

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
         call mpas_stream_mgr_read(domain % streamManager, whence=MPAS_STREAM_LATEST_BEFORE, ierr=err_tmp)
         err = ior(err, err_tmp)
         call mpas_stream_mgr_reset_alarms(domain % streamManager, direction=MPAS_STREAM_INPUT, ierr=err_tmp)
         err = ior(err, err_tmp)

         ! call analysis driver compute, restart, write subroutines (note: alarms and timers are handled by the analysis member code)
         call li_analysis_compute(domain, err_tmp)
         err = ior(err, err_tmp)
         call li_analysis_restart(domain, err_tmp)
         err = ior(err, err_tmp)
         call li_analysis_write(domain, err_tmp)
         err = ior(err, err_tmp)

         ! ===
         ! === Write Output and/or Restart, if needed
         ! ===
         call mpas_timer_start("write output")
         ! Update the restart_timestamp file with the new time, if needed.
         if ( mpas_stream_mgr_ringing_alarms(domain % streamManager, streamID='restart', direction=MPAS_STREAM_OUTPUT, ierr=err_tmp) ) then
            ! Need updated timestamp for writing restart timestamp
            currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, err_tmp)
            err = ior(err, err_tmp)
            call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=err_tmp)
            err = ior(err, err_tmp)
            ! write the timestamp to file
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
             call mpas_dmpar_global_abort("MPAS-landice: An error has occurred in li_core_run. Aborting...")
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
      use li_subglacial_hydro, only: li_SGH_finalize
   
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
 
      call mpas_timer_start("land ice finalize")

      err = 0
      err_tmp = 0
      globalErr = 0

      ! velocity finalize
      call li_velocity_finalize(domain, err_tmp)
      err = ior(err, err_tmp)

      ! subglacial hydro finalize
      call li_SGH_finalize(domain, err_tmp)
      err = ior(err, err_tmp)

      ! call finalize subroutine in analysis driver 
      call li_analysis_finalize(domain, err_tmp)
      err = ior(err, err_tmp)

      call mpas_destroy_clock(domain % clock, err_tmp)
      err = ior(err, err_tmp)

      call mpas_decomp_destroy_decomp_list(domain % decompositions)

      err = ior(err, err_tmp)

      ! === error check and exit
      call mpas_dmpar_max_int(domain % dminfo, err, globalErr)  ! Find out if any blocks got an error
      if (globalErr > 0) then
          call mpas_dmpar_global_abort("MPAS-landice: An error has occurred in li_core_finalize. Aborting...")
      endif

      call mpas_timer_stop("land ice finalize")

   end function li_core_finalize
   !--------------------------------------------------------------------



!***********************************************************************
!
!  function li_core_initial_solve
!
!> \brief   Performs the initial diagnostic solve for the LI core
!> \author  Matt Hoffman
!> \date    23 September 2015
!> \details
!>  This routine performs the initial diagnostic solve for the LI core.
!>  Rather than inlining these calculations, there are done in this
!>  routine to keep them modular.
!>  This has been made public so it can be called from an ESM.
!
!-----------------------------------------------------------------------
   function li_core_initial_solve(domain) result(err)

      use mpas_derived_types
      use mpas_pool_routines
      use mpas_kind_types
      use mpas_stream_manager
      use mpas_timer
      use li_diagnostic_vars
      use li_setup
      use li_statistics
      use mpas_io_streams, only: MPAS_STREAM_LATEST_BEFORE

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
      type (mpas_pool_type), pointer :: geometryPool
      logical, pointer :: config_do_restart, config_write_output_on_startup, config_write_stats_on_startup
      character(len=StrKIND), pointer :: config_velocity_solver
      ! Variables needed for printing timestamps
      type (MPAS_Time_Type) :: currTime
      character(len=StrKIND) :: timeStamp

      integer :: err, err_tmp, globalErr
      logical :: solveVelo

      integer, dimension(:), pointer :: vertexMask


      err = 0
      err_tmp = 0
      globalErr = 0

      ! Get Pool stuff that will be needed
      call mpas_pool_get_config(liConfigs, 'config_do_restart', config_do_restart)
      call mpas_pool_get_config(liConfigs, 'config_write_output_on_startup', config_write_output_on_startup)
      call mpas_pool_get_config(liConfigs, 'config_write_stats_on_startup', config_write_stats_on_startup)
      call mpas_pool_get_config(liConfigs, 'config_velocity_solver', config_velocity_solver)

      currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, err_tmp)
      err = ior(err, err_tmp)
      call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=err_tmp)
      err = ior(err, err_tmp)
      write(stderrUnit,*) 'Initial timestep ', trim(timeStamp)
      write(stdoutUnit,*) 'Initial timestep ', trim(timeStamp)


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

      call li_calculate_diagnostic_vars(domain, solveVelo=solveVelo, err=err_tmp)
      err = ior(err, err_tmp)

      call mpas_timer_stop("initial state calculation")

      if (config_write_stats_on_startup) then
         call mpas_timer_start("compute_statistics")
         call li_compute_statistics(domain, 0)     ! itimestep = 0
                                                   ! (itimestep is initialized below)
         call mpas_timer_stop("compute_statistics")
      endif

      ! compute analysis members on startup if option activated
      call mpas_timer_start("analysis member startup calculations")
          call li_analysis_compute_startup(domain, err_tmp)
          err = ior(err, err_tmp)
      call mpas_timer_stop("analysis member startup calculations")

      ! ===
      ! === Write Initial Output
      ! ===
      call mpas_timer_start("write output")

      if (config_write_output_on_startup) then
         call mpas_stream_mgr_write(domain % streamManager, 'output', forceWriteNow=.true., ierr=err_tmp)
      endif

      call mpas_timer_stop("write output")

      ! Move time level 1 fields (current values) into time level 2 (old values) for next time step
      ! (for those fields with multiple time levels)
      block => domain % blocklist
      do while(associated(block))
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_shift_time_levels(geometryPool)
         block => block % next
      end do

      if (config_do_restart .and. (trim(config_velocity_solver) /= 'sia')) then
         ! On a restart with the HO dycore, we need to make sure the FEM mesh will be rebuilt
         ! on the first time step.  Force this by setting the vertexMask at the end of the
         ! initial time to 0.  (Do this after writing output.)
         block => domain % blocklist
         do while(associated(block))
            call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
            call mpas_pool_get_array(geometryPool, 'vertexMask', vertexMask, timeLevel=2)  ! Get the old vertexMask
            vertexMask = 0
            block => block % next
         end do
      endif

      ! === error check and exit
      call mpas_dmpar_max_int(domain % dminfo, err, globalErr)  ! Find out if any blocks got an error
      if (globalErr > 0) then
          call mpas_dmpar_global_abort("MPAS-landice: An error has occurred in li_core_initial_solve. Aborting...")
      endif


   end function li_core_initial_solve
   !--------------------------------------------------------------------



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
   subroutine landice_init_block(block, dminfo, err)

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
      integer, intent(out) :: err  !< error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: geometryPool
      integer, dimension(:), pointer :: vertexMask
      character (len=StrKIND), pointer :: config_velocity_solver
      logical, pointer :: config_do_velocity_reconstruction_for_external_dycore
      logical, pointer :: config_adaptive_timestep_include_DCFL
      integer :: err_tmp

      err = 0
      err_tmp = 0

      ! Get pool stuff
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
      call mpas_pool_get_config(liConfigs, 'config_velocity_solver', config_velocity_solver)
      call mpas_pool_get_config(liConfigs, 'config_do_velocity_reconstruction_for_external_dycore', config_do_velocity_reconstruction_for_external_dycore)
      call mpas_pool_get_config(liConfigs, 'config_adaptive_timestep_include_DCFL', config_adaptive_timestep_include_DCFL)

      ! Copy data from first time level into all other time levels
      call mpas_pool_initialize_time_levels(geometryPool)

      ! Initialize vertexMask on time level 2 to 0, so diagnostic_solve_before_velocity in li_diagnostic_vars says that the vertexMask has changed (needed by external dycore)
      call mpas_pool_get_array(geometryPool, 'vertexMask', vertexMask, timeLevel = 2)
      vertexMask = 0

      ! ===
      ! === Call init routines ===
      ! ===
      call li_setup_vertical_grid(meshPool, geometryPool, err_tmp)
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

      ! Init for reconstruction of velocity
      !TODO - Initialize the reconstruction regardless of the velocity solver, because
      !        these coefficients will be needed in the future for IR transport.
      !       For now, we avoid calling these subroutines when using the external FO solver,
      !        because mpas_init_reconstruct fails with the MISMIP periodic mesh
      !        (since x_period and y_period are not set properly).
      if ( (trim(config_velocity_solver) == 'sia') .or.    &
           (trim(config_velocity_solver) == 'simple') .or.    &
           (trim(config_velocity_solver) == 'none') .or.    &
           config_do_velocity_reconstruction_for_external_dycore .or.   &
           config_adaptive_timestep_include_DCFL) then
         call mpas_rbf_interp_initialize(meshPool)
         call mpas_init_reconstruct(meshPool)
      endif

      ! Initialize velocity solver
      !WHL - This is now after the call to mpas_init_reconstruct, so that the reconstruction coefficients are available.
      call mpas_timer_start("initialize velocity")
      call li_velocity_block_init(block, err_tmp)
      err = ior(err, err_tmp)
      call mpas_timer_stop("initialize velocity")

      ! Mask init identifies initial ice extent
      call li_calculate_mask_init(geometryPool, err=err_tmp)
      err = ior(err, err_tmp)


      ! === error check
      if (err > 0) then
          write (stderrUnit,*) "An error has occurred in landice_init_block."
      endif

   !--------------------------------------------------------------------
   end subroutine landice_init_block



!***********************************************************************
!
!  routine simulation_clock_init
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
      type (MPAS_TimeInterval_type) :: adaptDtForceInterval
      character (len=StrKIND), pointer :: config_start_time, config_run_duration, config_stop_time, config_output_interval, config_restart_interval ! MPAS standard configs
      character (len=StrKIND), pointer :: config_dt  ! MPAS LI-specific config option
      character (len=StrKIND), pointer :: config_adaptive_timestep_force_interval  ! MPAS LI-specific config option
      character (len=StrKIND), pointer :: config_restart_timestamp_name
      character (len=StrKIND) :: restartTimeStamp !< string to be read from file
      integer, pointer :: config_year_digits
      integer :: err_tmp

      ierr = 0
      err_tmp = 0

      ! Adjust number of digits representing the year
      call mpas_pool_get_config(configs, 'config_year_digits', config_year_digits)
      call mpas_timekeeping_set_year_width(config_year_digits)

      ! Get configs needed for rest of routine
      call mpas_pool_get_config(configs, 'config_dt', config_dt)
      call mpas_pool_get_config(configs, 'config_start_time', config_start_time)
      call mpas_pool_get_config(configs, 'config_run_duration', config_run_duration)
      call mpas_pool_get_config(configs, 'config_stop_time', config_stop_time)
      call mpas_pool_get_config(configs, 'config_restart_timestamp_name', config_restart_timestamp_name)
      call mpas_pool_get_config(configs, 'config_adaptive_timestep_force_interval', config_adaptive_timestep_force_interval)


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

      ! Setup start/stop/duration times
      ! config_run_duration takes precedence over config_stop_time
      if (trim(config_run_duration) /= "none") then
         call mpas_set_timeInterval(runDuration, timeString=config_run_duration, ierr=err_tmp)
         ierr = ior(ierr,err_tmp)
         call mpas_create_clock(core_clock, startTime=startTime, timeStep=timeStep, runDuration=runDuration, ierr=err_tmp)
         ierr = ior(ierr,err_tmp)

         if (trim(config_stop_time) /= "none") then
            call mpas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, ierr=err_tmp)
            ierr = ior(ierr,err_tmp)
            if(startTime + runduration /= stopTime) then
               write(stderrUnit,*) 'Warning: config_run_duration and config_stop_time are inconsistent: using config_run_duration.'
            end if
         end if
      else if (trim(config_stop_time) /= "none") then
         call mpas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, ierr=err_tmp)
         ierr = ior(ierr,err_tmp)
         if (stopTime .lt. startTime) then
            write(stderrUnit,*) 'Error: config_stop_time is earlier than config_start_time!'
            ierr = 1
         endif
         call mpas_create_clock(core_clock, startTime=startTime, timeStep=timeStep, stopTime=stopTime, ierr=err_tmp)
         ierr = ior(ierr,err_tmp)
      else
         write(stderrUnit,*) 'Error: Neither config_run_duration nor config_stop_time were specified.'
         ierr = 1
      end if

      ! Set up the adaptiveTimestepForceInterval alarm.
      ! This is only needed if the adaptive time stepper is being used, but can be set up regardless.
      call mpas_set_timeInterval(adaptDtForceInterval, timeString=config_adaptive_timestep_force_interval, ierr=err_tmp)
      ierr = ior(ierr,err_tmp)
      call mpas_add_clock_alarm(core_clock, 'adaptiveTimestepForceInterval', alarmTime=startTime, alarmTimeInterval=adaptDtForceInterval, ierr=err_tmp)
      ierr = ior(ierr,err_tmp)
      ! Reset the alarm for checking for force setting of the adaptive timestep interval
      if (mpas_is_alarm_ringing(core_clock, 'adaptiveTimestepForceInterval', ierr=err_tmp)) then
         ierr = ior(ierr, err_tmp)
         call mpas_reset_clock_alarm(core_clock, 'adaptiveTimestepForceInterval', ierr=err_tmp)
         ierr = ior(ierr, err_tmp)
      endif
      ierr = ior(ierr, err_tmp)

      ! === error check
      if (ierr > 0) then
          write (stderrUnit,*) "An error has occurred in li_simulation_clock_init."
      endif

   !--------------------------------------------------------------------
   end subroutine li_simulation_clock_init



end module li_core
