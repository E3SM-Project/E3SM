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
!> \brief Main driver for MPAS ocean analysis core in post-processing mode
!> \author Mark Petersen
!> \date   November 2013
!> \details
!>  This module contains the drivers for the MPAS ocean analysis core in 
!>  post-processing mode.  None of these routines are used in run-time mode
!>  by the ocean core.
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
   use ocn_diagnostics
   use ocn_equation_of_state

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
!> \brief   Initialize MPAS-Ocean analysis post-processing core
!> \author  Mark Petersen
!> \date    November 2013
!> \details 
!>  This routine calls all initializations required to begin a
!>  simulation with MPAS-Ocean
!
!-----------------------------------------------------------------------

   subroutine mpas_core_init(domain, startTimeStamp)!{{{

      use mpas_grid_types

      implicit none

      type (domain_type), intent(inout) :: domain
      character(len=*), intent(out) :: startTimeStamp

      type (block_type), pointer :: block
      type (dm_info) :: dminfo

      integer :: err, err_tmp

      ! remove dt later
      real (kind=RKIND) :: dt

      dminfo = domain % dminfo
      err = 0

      ! Initialize submodules before initializing blocks.
      call ocn_equation_of_state_init(err_tmp)
      err = ior(err, err_tmp)

      call ocn_am_driver_init(domain, err_tmp)
      err = ior(err, err_tmp)

      call ocn_init_vert_coord(domain)

      call ocn_compute_max_level(domain)

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

      current_outfile_frames = 0

   end subroutine mpas_core_init!}}}

!***********************************************************************
!
!  routine mpas_simulation_clock_init
!
!> \brief   Initialize timer variables
!> \author  Mark Petersen
!> \date    November 2013
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

   end subroutine ocn_simulation_clock_init!}}}

!***********************************************************************
!
!  routine mpas_init_block
!
!> \brief   Initialize blocks within MPAS-Ocean analysis post-processing core
!> \author  Mark Petersen
!> \date    November 2013
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

!      call mpas_timer_start("diagnostic solve", .false., initDiagSolveTimer)
!      call ocn_diagnostic_solve(dt,  block % state % time_levs(1) % state, block % forcing, mesh, block % diagnostics, block % scratch)
!      call mpas_timer_stop("diagnostic solve", initDiagSolveTimer)

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
!> \author  Mark Petersen
!> \date    November 2013
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
   
      integer :: itimestep
      real (kind=RKIND) :: dt
      type (block_type), pointer :: block_ptr

      type (MPAS_Time_Type) :: currTime
      character(len=StrKIND) :: timeStamp
      integer :: err, ierr
   
      dt = config_dt

      currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
      call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
      write(stderrUnit,*) 'Initial time ', trim(timeStamp)

      ! fill in diagnostics variables
      call mpas_timer_start("diagnostic solve", .false., initDiagSolveTimer)
      block_ptr => domain % blocklist
      do while(associated(block_ptr))
         call ocn_diagnostic_solve(dt,  block_ptr % state % time_levs(1) % state, block_ptr % forcing, block_ptr % mesh, block_ptr % diagnostics, block_ptr % scratch)
         block_ptr => block_ptr % next
      end do
      call mpas_timer_stop("diagnostic solve", initDiagSolveTimer)

      if (config_write_output_on_startup) then
         call ocn_am_driver_compute(domain, err) 
         call ocn_write_output_frame(output_obj, output_frame, domain)
      endif

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

   subroutine mpas_core_finalize(domain)!{{{
   
      use mpas_grid_types
   
      implicit none

      type (domain_type), intent(inout) :: domain 
      integer :: ierr

      call mpas_destroy_clock(clock, ierr)

   end subroutine mpas_core_finalize!}}}

!***********************************************************************
!
!  routine mpas_core_setup_packages
!
!> \brief   Package setup routine
!> \author  Doug Jacobsen
!> \date    November 2013
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

      call ocn_am_driver_setup_packages(err_tmp)
      ierr = ior(ierr, err_tmp)

   end subroutine mpas_core_setup_packages!}}}

end module mpas_core

! vim: foldmethod=marker
