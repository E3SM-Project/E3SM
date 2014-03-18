! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_core
!
!> \brief Main driver for MPAS ocean core
!> \author Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date   September 2011
!> \details
!>  This module contains initialization and timestep drivers for 
!>  the MPAS ocean core.
!
!-----------------------------------------------------------------------

module mpas_core

   use mpas_configure
   use mpas_framework
   use mpas_timekeeping
   use mpas_dmpar
   use mpas_timer
   use mpas_io_units

   use ocn_am_driver
   use ocn_init

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
   use ocn_gm

   use ocn_equation_of_state

   use ocn_vmix

   use ocn_time_average

   use ocn_forcing
   use ocn_sea_ice

   use ocn_constants

   type (io_output_object), save :: restart_obj

   integer :: current_outfile_frames

   type (MPAS_Clock_type) :: clock

   integer, parameter :: outputAlarmID = 1
   integer, parameter :: restartAlarmID = 2
   integer, parameter :: statsAlarmID = 3
   integer, parameter :: coupleAlarmID = 4

   type (timer_node), pointer :: globalDiagTimer, timeIntTimer, testSuiteTimer
   type (timer_node), pointer :: initDiagSolveTimer

   contains

!***********************************************************************
!
!  routine mpas_core_init
!
!> \brief   Initialize MPAS-Ocean core
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details 
!>  This routine calls all initializations required to begin a
!>  simulation with MPAS-Ocean
!
!-----------------------------------------------------------------------

   subroutine mpas_core_init(domain, startTimeStamp)!{{{

      use mpas_grid_types
      use ocn_tracer_advection

      implicit none

      type (domain_type), intent(inout) :: domain
      character(len=*), intent(out) :: startTimeStamp

      real (kind=RKIND) :: dt
      type (block_type), pointer :: block
      type (dm_info) :: dminfo

      integer :: err, err_tmp
      real (kind=RKIND) :: maxDensity, maxDensity_global

      dminfo = domain % dminfo

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

      call ocn_vmix_init(domain,err_tmp)
      err = ior(err, err_tmp)

      call ocn_equation_of_state_init(err_tmp)
      err = ior(err, err_tmp)

      call ocn_tendency_init(err_tmp)
      err = ior(err,err_tmp)
      call ocn_diagnostics_init(err_tmp)
      err = ior(err,err_tmp)

      call ocn_forcing_init(err_tmp)
      err = ior(err,err_tmp)


      call ocn_global_diagnostics_init(dminfo,err_tmp)
      err = ior(err, err_tmp)

      call ocn_sea_ice_init(domain % blocklist % mesh % nVertLevels, err_tmp)
      err = ior(err, err_tmp)

      call ocn_am_driver_init(domain, err_tmp)
      err = ior(err, err_tmp)

      call ocn_constants_init()

      call mpas_timer_init(domain)

      if(err.eq.1) then
          call mpas_dmpar_abort(dminfo)
      endif

      call ocn_init_vert_coord(domain)

      call ocn_compute_max_level(domain)

      if (.not.config_do_restart) call ocn_init_split_timestep(domain)

      write (stdoutUnit,'(a,a)') ' Vertical coordinate movement is: ',trim(config_vert_coord_movement)

      if (config_vert_coord_movement.ne.'fixed'.and. &
          config_vert_coord_movement.ne.'uniform_stretching'.and. &
          config_vert_coord_movement.ne.'impermeable_interfaces'.and. &
          config_vert_coord_movement.ne.'user_specified') then
         write (stderrUnit,*) ' Incorrect choice of config_vert_coord_movement.'
         call mpas_dmpar_abort(dminfo)
      endif

      write (stdoutUnit,'(a,a)') ' Pressure type is: ',trim(config_pressure_gradient_type)
      if (config_pressure_gradient_type.ne.'pressure_and_zmid'.and. &
          config_pressure_gradient_type.ne.'MontgomeryPotential') then
         write (stderrUnit,*) ' Incorrect choice of config_pressure_gradient_type.'
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
          maxDensity = max(maxDensity, maxval(block % mesh % meshDensity % array))
          block => block % next
        end do
        call mpas_dmpar_max_real(domain % dminfo, maxDensity, maxDensity_global)
        config_maxMeshDensity = maxDensity_global
      endif

      !
      ! Initialize core
      !
      dt = config_dt

      call ocn_simulation_clock_init(domain, dt, startTimeStamp)

      block => domain % blocklist
      do while (associated(block))
         call mpas_init_block(block, block % mesh, dt, err)
         if(err.eq.1) then
             call mpas_dmpar_abort(dminfo)
         endif
         block % diagnostics % xtime % scalar = startTimeStamp 
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

      current_outfile_frames = 0

   end subroutine mpas_core_init!}}}

!***********************************************************************
!
!  routine mpas_simulation_clock_init
!
!> \brief   Initialize timer variables
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details 
!>  This routine initializes all timer variables
!
!-----------------------------------------------------------------------

   subroutine ocn_simulation_clock_init(domain, dt, startTimeStamp)!{{{

      implicit none

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(in) :: dt
      character(len=*), intent(out) :: startTimeStamp

      type (MPAS_Time_Type) :: startTime, stopTime, alarmStartTime
      type (MPAS_TimeInterval_type) :: runDuration, timeStep, alarmTimeStep
      character(len=StrKIND) :: restartTimeStamp
      integer :: ierr

      if(config_start_time == 'file') then
        open(22,file=config_restart_timestamp_name,form='formatted',status='old')
        read(22,*) restartTimeStamp
        close(22)
        call mpas_set_time(curr_time=startTime, dateTimeString=restartTimeStamp, ierr=ierr)
      else
        call mpas_set_time(curr_time=startTime, dateTimeString=config_start_time, ierr=ierr)
      end if

      call mpas_set_timeInterval(timeStep, dt=dt, ierr=ierr)
      if (trim(config_run_duration) /= "none") then
         call mpas_set_timeInterval(runDuration, timeString=config_run_duration, ierr=ierr)
         call mpas_create_clock(clock, startTime=startTime, timeStep=timeStep, runDuration=runDuration, ierr=ierr)

         if (trim(config_stop_time) /= "none") then
            call mpas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, ierr=ierr)
            if(startTime + runduration /= stopTime) then
               write(stderrUnit,*) 'Warning: config_run_duration and config_stop_time are inconsitent: using config_run_duration.'
            end if
         end if
      else if (trim(config_stop_time) /= "none") then
         call mpas_set_time(curr_time=stopTime, dateTimeString=config_stop_time, ierr=ierr)
         call mpas_create_clock(clock, startTime=startTime, timeStep=timeStep, stopTime=stopTime, ierr=ierr)
      else
          write(stderrUnit,*) 'Error: Neither config_run_duration nor config_stop_time were specified.'
          call mpas_dmpar_finalize(domain % dminfo)
      end if

      ! set output alarm
      call mpas_set_timeInterval(alarmTimeStep, timeString=config_output_interval, ierr=ierr)
      alarmStartTime = startTime + alarmTimeStep
      call mpas_add_clock_alarm(clock, outputAlarmID, alarmStartTime, alarmTimeStep, ierr=ierr)
print *, 'output alarmTimeStep', alarmTimeStep, ierr

      ! set restart alarm, if necessary
      if (trim(config_restart_interval) /= "none") then
         call mpas_set_timeInterval(alarmTimeStep, timeString=config_restart_interval, ierr=ierr)
         alarmStartTime = startTime + alarmTimeStep
         call mpas_add_clock_alarm(clock, restartAlarmID, alarmStartTime, alarmTimeStep, ierr=ierr)
      end if

      !TODO: use this code if we desire to convert config_stats_interval to alarms 
      !(must also change config_stats_interval type to character) 
      ! set stats alarm, if necessary
      if (trim(config_stats_interval) /= "none") then      
         call mpas_set_timeInterval(alarmTimeStep, timeString=config_stats_interval, ierr=ierr)
         alarmStartTime = startTime + alarmTimeStep
         call mpas_add_clock_alarm(clock, statsAlarmID, alarmStartTime, alarmTimeStep, ierr=ierr)
      end if

      call mpas_get_time(curr_time=startTime, dateTimeString=startTimeStamp, ierr=ierr)

      call ocn_am_driver_init_alarms(clock,startTime,ierr)

   end subroutine ocn_simulation_clock_init!}}}

!***********************************************************************
!
!  routine mpas_init_block
!
!> \brief   Initialize blocks within MPAS-Ocean core
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details 
!>  This routine calls all block-level initializations required to begin a
!>  simulation with MPAS-Ocean
!
!-----------------------------------------------------------------------

   subroutine mpas_init_block(block, mesh, dt, err)!{{{
   
      use mpas_grid_types
      use mpas_rbf_interpolation
      use mpas_vector_operations
      use mpas_vector_reconstruction
      use mpas_tracer_advection_helpers
   
      implicit none
   
      type (block_type), intent(inout) :: block
      type (mesh_type), intent(inout) :: mesh
      real (kind=RKIND), intent(in) :: dt
      integer, intent(out) :: err
      integer :: i, iEdge, iCell, k
      integer :: err1
   
      call ocn_setup_sign_and_index_fields(mesh)
      call mpas_initialize_deriv_two(mesh, mesh % derivTwo % array, err)
      call mpas_tracer_advection_coefficients(mesh, &
          config_horiz_tracer_adv_order, mesh % derivTwo % array, mesh % advCoefs % array, &
          mesh % advCoefs3rd % array, mesh % nAdvCellsForEdge % array, mesh % advCellsForEdge % array, &
          err1, mesh % maxLevelCell % array, mesh % highOrderAdvectionMask % array, &
          mesh % boundaryCell % array)
      err = ior(err, err1)

      call ocn_time_average_init(block % average)
   
      call mpas_timer_start("diagnostic solve", .false., initDiagSolveTimer)
      call ocn_diagnostic_solve(dt,  block % state % time_levs(1) % state, block % forcing, mesh, block % diagnostics, block % scratch)
      call mpas_timer_stop("diagnostic solve", initDiagSolveTimer)

      ! Compute velocity transport, used in advection terms of layerThickness and tracer tendency
      block % diagnostics % uTransport % array(:,:) &
            = block % state % time_levs(1) % state % normalVelocity % array(:,:) &
            + block % diagnostics % uBolusGM % array(:,:)

      call ocn_compute_mesh_scaling(mesh)
 
      call mpas_rbf_interp_initialize(mesh)
      call mpas_initialize_tangent_vectors(mesh, mesh % edgeTangentVectors % array)

      call mpas_init_reconstruct(mesh)
      call mpas_reconstruct(mesh, block % state % time_levs(1) % state % normalVelocity % array,                  &
                       block % diagnostics % normalVelocityX % array,            &
                       block % diagnostics % normalVelocityY % array,            &
                       block % diagnostics % normalVelocityZ % array,            &
                       block % diagnostics % normalVelocityZonal % array,        &
                       block % diagnostics % normalVelocityMeridional % array    &
                      )

      ! initialize velocities and tracers on land to be zero.

      block % mesh % areaCell % array(block % mesh % nCells+1) = -1.0e34

      block % state % time_levs(1) % state % layerThickness % array(:,block % mesh % nCells+1) = 0.0

      do iEdge=1,block % mesh % nEdges
         block % state % time_levs(1) % state % normalVelocity % array( &
             block % mesh % maxLevelEdgeTop % array(iEdge)+1 &
            :block % mesh % maxLevelEdgeBot % array(iEdge), iEdge) = 0.0

         block % state % time_levs(1) % state % normalVelocity % array( &
             block % mesh % maxLevelEdgeBot % array(iEdge)+1: &
             block % mesh % nVertLevels,iEdge) = -1.0e34
      end do
      do iCell=1,block % mesh % nCells
         block % state % time_levs(1) % state % tracers % array( &
            :, block % mesh % maxLevelCell % array(iCell)+1 &
              :block % mesh % nVertLevels,iCell) =  -1.0e34
      end do

      do i=2,nTimeLevs
        call mpas_copy_state(block % state % time_levs(i) % state, &
                             block % state % time_levs(1) % state)
      end do

   end subroutine mpas_init_block!}}}
   
!***********************************************************************
!
!  routine mpas_core_run
!
!> \brief   Main driver for MPAS-Ocean time-stepping
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details 
!>  This routine includes the time-stepping loop, and calls timer
!>  routines to write output and restart files.
!
!-----------------------------------------------------------------------

   subroutine mpas_core_run(domain, output_obj, output_frame)!{{{
   
      use mpas_kind_types
      use mpas_grid_types
      use mpas_io_output
      use mpas_timer
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      type (io_output_object), intent(inout) :: output_obj
      integer, intent(inout) :: output_frame
   
      integer :: itimestep, err
      real (kind=RKIND) :: dt
      type (block_type), pointer :: block_ptr

      type (MPAS_Time_Type) :: currTime
      character(len=StrKIND) :: timeStamp
      integer :: ierr
   
      ! Eventually, dt should be domain specific
      dt = config_dt

      currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
      call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
      write(stderrUnit,*) 'Initial time ', trim(timeStamp)

      call ocn_am_driver_compute_startup(domain, err) 

      if (config_write_output_on_startup) then
         call ocn_write_output_frame(output_obj, output_frame, domain)
      endif

      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call ocn_time_average_init(block_ptr % average)
        block_ptr => block_ptr % next
      end do

      ! During integration, time level 1 stores the model state at the beginning of the
      !   time step, and time level 2 stores the state advanced dt in time by timestep(...)
      itimestep = 0
      do while (.not. mpas_is_clock_stop_time(clock))

         itimestep = itimestep + 1
         call mpas_advance_clock(clock)

         currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
         write(stderrUnit,*) 'Doing timestep ', trim(timeStamp)
   
         block_ptr => domain % blocklist
         do while(associated(block_ptr))
           call ocn_forcing_build_arrays(block_ptr % mesh, block_ptr % state % time_levs(1) % state, block_ptr % forcing, ierr)
           call ocn_forcing_build_transmission_array(block_ptr % mesh, block_ptr % state % time_levs(1) % state, block_ptr % forcing, ierr)
           block_ptr => block_ptr % next
         end do

         call mpas_timer_start("time integration", .false., timeIntTimer)
         call mpas_timestep(domain, itimestep, dt, timeStamp)
         call mpas_timer_stop("time integration", timeIntTimer)
   
         ! Move time level 2 fields back into time level 1 for next time step
         block_ptr => domain % blocklist
         do while(associated(block_ptr))
            call mpas_shift_time_levels_state(block_ptr % state)
            block_ptr => block_ptr % next
         end do
      
         call ocn_am_driver_compute_w_alarms(clock, domain, err) 

         if (mpas_is_alarm_ringing(clock, outputAlarmID, ierr=ierr)) then
            call mpas_reset_clock_alarm(clock, outputAlarmID, ierr=ierr)
            ! output_frame will always be > 1 here unless it was reset after the 
            ! maximum number of frames per outfile was reached.
            if(output_frame == 1) then
               call mpas_output_state_finalize(output_obj, domain % dminfo)
               call mpas_output_state_init(output_obj, domain, "OUTPUT", trim(timeStamp))
            end if

            block_ptr => domain % blocklist
            do while (associated(block_ptr))
                call ocn_time_average_normalize(block_ptr % average)
                block_ptr => block_ptr % next
            end do

            call ocn_write_output_frame(output_obj, output_frame, domain)

            block_ptr => domain % blocklist
            do while (associated(block_ptr))
                call ocn_time_average_init(block_ptr % average)
                block_ptr => block_ptr % next
            end do
         end if

         if (mpas_is_alarm_ringing(clock, restartAlarmID, ierr=ierr)) then
            call mpas_reset_clock_alarm(clock, restartAlarmID, ierr=ierr)

            ! Write one restart time per file
            call mpas_output_state_init(restart_obj, domain, "RESTART", trim(timeStamp))
            call mpas_output_state_for_domain(restart_obj, domain, 1)
            call mpas_output_state_finalize(restart_obj, domain % dminfo)

            call ocn_am_driver_restart(domain, err)
         end if

      end do

   end subroutine mpas_core_run!}}}
   
   subroutine ocn_write_output_frame(output_obj, output_frame, domain)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute diagnostic fields for a domain and write model state to output file
   !
   ! Input/Output: domain - contains model state; diagnostic field are computed
   !                        before returning
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      use mpas_grid_types
      use mpas_io_output
      use mpas_vector_reconstruction
   
      implicit none
   
      integer, intent(inout) :: output_frame
      type (domain_type), intent(inout) :: domain
      type (io_output_object), intent(inout) :: output_obj
   
      integer :: i, j, k
      integer :: eoe
      type (block_type), pointer :: block_ptr
   
      ! Compute output diagnostics
      block_ptr => domain % blocklist
      do while (associated(block_ptr))
         call ocn_compute_output_diagnostics(block_ptr % state % time_levs(1) % state, block_ptr % mesh)

         call mpas_allocate_scratch_field(block_ptr % scratch % windStressFull, .true.)
         call mpas_allocate_scratch_field(block_ptr % scratch % windStressX, .true.)
         call mpas_allocate_scratch_field(block_ptr % scratch % windStressY, .true.)
         call mpas_allocate_scratch_field(block_ptr % scratch % windStressZ, .true.)
         call mpas_allocate_scratch_field(block_ptr % scratch % windStressZonal, .true.)
         call mpas_allocate_scratch_field(block_ptr % scratch % windStressMeridional, .true.)

         block_ptr % scratch % windStressFull % array = 0.0_RKIND
         block_ptr % scratch % windStressFull % array(1,:) = block_ptr % forcing % surfaceWindStress % array

         call mpas_reconstruct(block_ptr % mesh, block_ptr % scratch % windStressFull % array, &
                               block_ptr % scratch % windStressX % array, &
                               block_ptr % scratch % windStressY % array, &
                               block_ptr % scratch % windStressZ % array, &
                               block_ptr % scratch % windStressZonal % array, &
                               block_ptr % scratch % windStressMeridional % array)

         block_ptr % diagnostics % windStressZonalDiag % array(:) = block_ptr % scratch % windStressZonal % array(1,:)
         block_ptr % diagnostics % windStressMeridionalDiag % array(:) = block_ptr % scratch % windStressMeridional % array(1,:)

         call mpas_deallocate_scratch_field(block_ptr % scratch % windStressFull, .true.)
         call mpas_deallocate_scratch_field(block_ptr % scratch % windStressX, .true.)
         call mpas_deallocate_scratch_field(block_ptr % scratch % windStressY, .true.)
         call mpas_deallocate_scratch_field(block_ptr % scratch % windStressZ, .true.)
         call mpas_deallocate_scratch_field(block_ptr % scratch % windStressZonal, .true.)
         call mpas_deallocate_scratch_field(block_ptr % scratch % windStressMeridional, .true.)

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
   
   end subroutine ocn_write_output_frame!}}}
   
   subroutine ocn_compute_output_diagnostics(state, mesh)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute diagnostic fields for a domain
   !
   ! Input: state - contains model prognostic fields
   !        mesh  - contains mesh metadata
   !
   ! Output: state - upon returning, diagnostic fields will have be computed
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      use mpas_grid_types
   
      implicit none
   
      type (state_type), intent(inout) :: state
      type (mesh_type), intent(in) :: mesh
   
      integer :: i, eoe
      integer :: iEdge, k
   
   end subroutine ocn_compute_output_diagnostics!}}}
   
!***********************************************************************
!
!  routine mpas_core_run
!
!> \brief   Sub-driver for MPAS-Ocean time-stepping
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details 
!>  This routine calls the time integration routine within a time-stepping loop.
!
!-----------------------------------------------------------------------

   subroutine mpas_timestep(domain, itimestep, dt, timeStamp)!{{{
   
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

   end subroutine mpas_timestep!}}}

   subroutine ocn_init_split_timestep(domain)!{{{
   ! Initialize splitting variables

      use mpas_grid_types
      use mpas_configure

      implicit none

      type (domain_type), intent(inout) :: domain

      integer :: i, iCell, iEdge, iVertex, k
      type (block_type), pointer :: block

      integer :: iTracer, cell, cell1, cell2
      real (kind=RKIND) :: normalThicknessFluxSum, layerThicknessSum, layerThicknessEdge1
      real (kind=RKIND), dimension(:), pointer :: refBottomDepth
         
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness
      integer :: nVertLevels

      ! Initialize z-level mesh variables from h, read in from input file.
      block => domain % blocklist
      do while (associated(block))

         layerThickness          => block % state % time_levs(1) % state % layerThickness % array
         refBottomDepth => block % mesh % refBottomDepth % array
         nVertLevels = block % mesh % nVertLevels

         ! Compute barotropic velocity at first timestep
         ! This is only done upon start-up.
         if (trim(config_time_integrator) == 'unsplit_explicit') then
            block % state % time_levs(1) % state % normalBarotropicVelocity % array(:) = 0.0

              block % state % time_levs(1) % state % normalBaroclinicVelocity % array(:,:) &
            = block % state % time_levs(1) % state % normalVelocity % array(:,:) 

         elseif (trim(config_time_integrator) == 'split_explicit') then

            if (config_filter_btr_mode) then
               do iCell=1,block % mesh % nCells
                  block % state % time_levs(1) % state % layerThickness % array(1,iCell) & 
                = block % mesh % refBottomDepth % array(1)
               enddo
            endif 

            do iEdge=1,block % mesh % nEdges
               cell1 = block % mesh % cellsOnEdge % array(1,iEdge)
               cell2 = block % mesh % cellsOnEdge % array(2,iEdge)

               ! normalBarotropicVelocity = sum(u)/sum(h) on each column
               ! ocn_diagnostic_solve has not yet been called, so compute hEdge 
               ! just for this edge.

               ! thicknessSum is initialized outside the loop because on land boundaries 
               ! maxLevelEdgeTop=0, but I want to initialize thicknessSum with a 
               ! nonzero value to avoid a NaN.
               layerThicknessEdge1 = 0.5*( &
                   block % state % time_levs(1) % state % layerThickness % array(1,cell1) & 
                 + block % state % time_levs(1) % state % layerThickness % array(1,cell2) ) 
               normalThicknessFluxSum = layerThicknessEdge1*block % state % time_levs(1) % state % normalVelocity % array(1,iEdge)
               layerThicknessSum = layerThicknessEdge1

               do k=2,block % mesh % maxLevelEdgeTop % array(iEdge)
                  ! ocn_diagnostic_solve has not yet been called, so compute hEdge 
                  ! just for this edge.
                  layerThicknessEdge1 = 0.5*( &
                      block % state % time_levs(1) % state % layerThickness % array(k,cell1) & 
                    + block % state % time_levs(1) % state % layerThickness % array(k,cell2) ) 

                  normalThicknessFluxSum = normalThicknessFluxSum &
                     + layerThicknessEdge1*block % state % time_levs(1) % state % normalVelocity % array(k,iEdge)
                  layerThicknessSum = layerThicknessSum + layerThicknessEdge1

               enddo
               block % state % time_levs(1) % state % normalBarotropicVelocity % array(iEdge) = normalThicknessFluxSum/layerThicknessSum

               ! normalBaroclinicVelocity(k,iEdge) = normalVelocity(k,iEdge) - normalBarotropicVelocity(iEdge)
               do k=1,block % mesh % maxLevelEdgeTop % array(iEdge)
                 block % state % time_levs(1) % state % normalBaroclinicVelocity % array(k,iEdge) &
                 = block % state % time_levs(1) % state % normalVelocity % array(k,iEdge) &
                 - block % state % time_levs(1) % state % normalBarotropicVelocity % array(iEdge)
               enddo

               ! normalBaroclinicVelocity=0, normalVelocity=0 on land cells
               do k=block % mesh % maxLevelEdgeTop % array(iEdge)+1, block % mesh % nVertLevels
                 block % state % time_levs(1) % state % normalBaroclinicVelocity % array(k,iEdge) = 0.0
                 block % state % time_levs(1) % state % normalVelocity % array(k,iEdge) = 0.0
               enddo
            enddo

            if (config_filter_btr_mode) then
               ! filter normalBarotropicVelocity out of initial condition
                block % state % time_levs(1) % state % normalVelocity % array(:,:) &
              = block % state % time_levs(1) % state % normalBaroclinicVelocity % array(:,:)

               block % state % time_levs(1) % state % normalBarotropicVelocity % array(:) = 0.0
            endif 

         endif

      block => block % next
      end do

   end subroutine ocn_init_split_timestep!}}}
 
   subroutine mpas_core_finalize(domain)!{{{
   
      use mpas_grid_types
   
      implicit none

      type (domain_type), intent(inout) :: domain 
      integer :: ierr

      call ocn_am_driver_finalize(domain, ierr) 

      call mpas_destroy_clock(clock, ierr)

   end subroutine mpas_core_finalize!}}}

!***********************************************************************
!
!  routine mpas_core_setup_packages
!
!> \brief   Package setup routine
!> \author  Doug Jacobsen
!> \date    September 2011
!> \details 
!>  This routine is intended to correctly configure the packages for this MPAS
!>   core. It can use any Fortran logic to properly configure packages, and it
!>   can also make use of any namelist options. All variables in the model are
!>   *not* allocated until after this routine is called.
!
!-----------------------------------------------------------------------
   subroutine mpas_core_setup_packages(ierr)!{{{

      use mpas_packages

      implicit none

      integer, intent(out) :: ierr

      integer :: err_tmp

      ierr = 0

      if (config_use_freq_filtered_thickness) then
         thicknessFilterActive = .true.
      end if

      if (config_time_integrator == trim('split_explicit') &
         .or. config_time_integrator == trim('unsplit_explicit') ) then

         splitTimeIntegratorActive = .true.
      end if

      if (config_forcing_type == trim('restoring')) then
         surfaceRestoringActive = .true.
      else if (config_forcing_type == trim('bulk')) then
         bulkForcingActive = .true.
      end if

      if (config_frazil_ice_formation) then
         frazilIceActive = .true.
      end if

      call ocn_am_driver_setup_packages(err_tmp)
      ierr = ior(ierr, err_tmp)

   end subroutine mpas_core_setup_packages!}}}

end module mpas_core

! vim: foldmethod=marker
