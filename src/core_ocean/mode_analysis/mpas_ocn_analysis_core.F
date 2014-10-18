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

   use mpas_framework
   use mpas_timekeeping
   use mpas_dmpar
   use mpas_timer
   use mpas_io_units

   use ocn_analysis_driver
   use ocn_init
   use ocn_diagnostics
   use ocn_equation_of_state
   use ocn_constants
   use ocn_time_average

   type (MPAS_Clock_type), pointer :: clock

   character(len=*), parameter :: statsAlarmID = 'stats'
   character(len=*), parameter :: coupleAlarmID = 'couple'

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

   subroutine mpas_core_init(domain, stream_manager, startTimeStamp)!{{{

      use mpas_grid_types
      use mpas_stream_manager

      implicit none

      type (domain_type), intent(inout) :: domain
      type (MPAS_streamManager_type), intent(inout) :: stream_manager
      character(len=*), intent(out) :: startTimeStamp

      type (block_type), pointer :: block
      type (dm_info) :: dminfo
      type (mpas_pool_type), pointer :: diagnosticsPool

      integer :: err, err_tmp

      ! remove dt later
      real (kind=RKIND) :: dt
      real (kind=RKIND), pointer :: config_dt
      character (len=StrKIND), pointer :: xtime
      type (MPAS_Time_Type) :: startTime


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
      call mpas_pool_get_config(domain % configs, 'config_dt', config_dt)

      !
      ! Read input data for model
      !
      call MPAS_stream_mgr_read(stream_manager, streamID='input', ierr=err)
      call MPAS_stream_mgr_reset_alarms(stream_manager, streamID='input', ierr=err)

      dminfo = domain % dminfo
      err = 0

      ! Initialize submodules before initializing blocks.
      call ocn_equation_of_state_init(err_tmp)
      err = ior(err, err_tmp)

      call ocn_analysis_init(domain, err_tmp)
      err = ior(err, err_tmp)

      call ocn_init_vert_coord(domain)

      call ocn_compute_max_level(domain)

      dt = config_dt

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_array(diagnosticsPool, 'xtime', xtime)

         call mpas_init_block(block, dt, err)
         if(err.eq.1) then
             call mpas_dmpar_abort(dminfo)
         endif

         xtime = startTimeStamp 
         block => block % next
      end do

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

   subroutine ocn_simulation_clock_init(core_clock, configs, ierr)!{{{

      implicit none

      type (MPAS_Clock_type), intent(inout) :: core_clock
      type (mpas_pool_type), intent(inout) :: configs
      integer, intent(out) :: ierr

      type (MPAS_Time_Type) :: startTime, stopTime, alarmStartTime
      type (MPAS_TimeInterval_type) :: runDuration, timeStep, alarmTimeStep
      character(len=StrKIND) :: restartTimeStamp
      integer :: err_tmp
      character (len=StrKIND), pointer :: config_start_time, config_stop_time
      character (len=StrKIND), pointer :: config_run_duration
      character (len=StrKIND), pointer :: config_stats_interval
      real(kind=RKIND), pointer :: config_dt


      ierr = 0

      call mpas_pool_get_config(configs, 'config_dt', config_dt)
      call mpas_pool_get_config(configs, 'config_start_time', config_start_time)
      call mpas_pool_get_config(configs, 'config_stop_time', config_stop_time)
      call mpas_pool_get_config(configs, 'config_run_duration', config_run_duration)
      call mpas_pool_get_config(configs, 'config_stats_interval', config_stats_interval)

      call mpas_set_time(curr_time=startTime, dateTimeString=config_start_time, ierr=err_tmp)

      call mpas_set_timeInterval(timeStep, dt=config_dt, ierr=err_tmp)
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
   subroutine mpas_init_block(block, dt, err)!{{{
   
      use mpas_grid_types
      use mpas_rbf_interpolation
      use mpas_vector_operations
      use mpas_vector_reconstruction
      use mpas_tracer_advection_helpers
   
      implicit none
   
      type (block_type), intent(inout) :: block
      real (kind=RKIND), intent(in) :: dt
      integer, intent(out) :: err

      type (mpas_pool_type), pointer :: meshPool, averagePool, statePool
      type (mpas_pool_type), pointer :: forcingPool, diagnosticsPool, scratchPool
      integer :: i, iEdge, iCell, k
      integer :: err1

      integer, dimension(:), pointer :: nAdvCellsForEdge, maxLevelCell
      integer, dimension(:), pointer :: maxLevelEdgeBot, maxLevelEdgeTop
      integer, dimension(:,:), pointer :: advCellsForEdge, highOrderAdvectionMask, boundaryCell
      real (kind=RKIND), dimension(:), pointer :: areaCell
      real (kind=RKIND), dimension(:,:), pointer :: advCoefs, advCoefs3rd, uTransport
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity, uBolusGM, edgeTangentVectors
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocityX, normalVelocityY, normalVelocityZ
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocityZonal, normalVelocityMeridional
      real (kind=RKIND), dimension(:,:,:), pointer :: derivTwo

      real (kind=RKIND), dimension(:,:,:), pointer :: tracers

      integer, pointer :: nCells, nEdges, nVertices, nVertLevels
      integer, pointer :: config_horiz_tracer_adv_order
      logical, pointer :: config_hmix_scaleWithMesh
      real (kind=RKIND), pointer :: config_maxMeshDensity

      call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
      call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
      call mpas_pool_get_dimension(block % dimensions, 'nVertices', nVertices)
      call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)

      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
      call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
      call mpas_pool_get_subpool(block % structs, 'average', averagePool)

      call mpas_pool_get_array(meshPool, 'derivTwo', derivTwo)
      call mpas_pool_get_array(meshPool, 'advCoefs', advCoefs)
      call mpas_pool_get_array(meshPool, 'advCoefs3rd', advCoefs3rd)
      call mpas_pool_get_array(meshPool, 'nAdvCellsForEdge', nAdvCellsForEdge)
      call mpas_pool_get_array(meshPool, 'advCellsForEdge', advCellsForEdge)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'highOrderAdvectionMask', highOrderAdvectionMask)
      call mpas_pool_get_array(meshPool, 'boundaryCell', boundaryCell)
      call mpas_pool_get_array(meshPool, 'edgeTangentVectors', edgeTangentVectors)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'boundaryCell', boundaryCell)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeBot', maxLevelEdgeBot)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)

      call mpas_pool_get_array(diagnosticsPool, 'uTransport', uTransport)
      call mpas_pool_get_array(diagnosticsPool, 'uBolusGM', uBolusGM)
      call mpas_pool_get_array(diagnosticsPool, 'normalVelocityX', normalVelocityX)
      call mpas_pool_get_array(diagnosticsPool, 'normalVelocityY', normalVelocityY)
      call mpas_pool_get_array(diagnosticsPool, 'normalVelocityZ', normalVelocityZ)
      call mpas_pool_get_array(diagnosticsPool, 'normalVelocityZonal', normalVelocityZonal)
      call mpas_pool_get_array(diagnosticsPool, 'normalVelocityMeridional', normalVelocityMeridional)

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, 1)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
      call mpas_pool_get_array(statePool, 'tracers', tracers, 1)

      call mpas_pool_get_config(block % configs, 'config_horiz_tracer_adv_order', config_horiz_tracer_adv_order)
      call mpas_pool_get_config(block % configs, 'config_hmix_scaleWithMesh', config_hmix_scaleWithMesh)
      call mpas_pool_get_config(block % configs, 'config_maxMeshDensity', config_maxMeshDensity)
   
      call ocn_setup_sign_and_index_fields(meshPool)
      call mpas_initialize_deriv_two(meshPool, derivTwo, err)
      call mpas_tracer_advection_coefficients(meshPool, &
          config_horiz_tracer_adv_order, derivTwo, advCoefs, &
          advCoefs3rd, nAdvCellsForEdge, advCellsForEdge, &
          err1, maxLevelCell, highOrderAdvectionMask, &
          boundaryCell)
      err = ior(err, err1)

      call ocn_time_average_init(averagePool)
   
      call mpas_timer_start("diagnostic solve", .false., initDiagSolveTimer)
      call ocn_diagnostic_solve(dt,  statePool, forcingPool, meshPool, diagnosticsPool, scratchPool)
      call mpas_timer_stop("diagnostic solve", initDiagSolveTimer)

      ! Compute velocity transport, used in advection terms of layerThickness and tracer tendency
      uTransport(:,:) = normalVelocity(:,:) + uBolusGM(:,:)

      call ocn_compute_mesh_scaling(meshPool, config_hmix_scaleWithMesh, config_maxMeshDensity)
 
      call mpas_rbf_interp_initialize(meshPool)
      call mpas_initialize_tangent_vectors(meshPool, edgeTangentVectors)

      call mpas_init_reconstruct(meshPool)
      call mpas_reconstruct(meshPool, normalVelocity,        &
                       normalVelocityX,            &
                       normalVelocityY,            &
                       normalVelocityZ,            &
                       normalVelocityZonal,        &
                       normalVelocityMeridional    &
                      )

      ! initialize velocities and tracers on land to be zero.
      areaCell(nCells+1) = -1.0e34

      layerThickness(:, nCells+1) = 0.0

      do iEdge=1, nEdges
         normalVelocity(maxLevelEdgeTop(iEdge)+1:maxLevelEdgeBot(iEdge), iEdge) = 0.0

         normalVelocity(maxLevelEdgeBot(iEdge)+1:nVertLevels,iEdge) = -1.0e34
      end do

      do iCell=1,nCells
         tracers(:, maxLevelCell(iCell)+1:nVertLevels,iCell) =  -1.0e34
      end do

      call mpas_pool_initialize_time_levels(statePool)

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

   subroutine mpas_core_run(domain, stream_manager)!{{{
   
      use mpas_kind_types
      use mpas_stream_manager
      use mpas_grid_types
      use mpas_timer
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      type (MPAS_streamManager_type), intent(inout) :: stream_manager
   
      integer :: itimestep
      real (kind=RKIND) :: dt
      type (block_type), pointer :: block_ptr

      type (MPAS_Time_Type) :: currTime
      character(len=StrKIND) :: timeStamp
      integer :: err, ierr

      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: scratchPool

      real (kind=RKIND), pointer :: config_dt
      logical, pointer :: config_write_output_on_startup

      call mpas_pool_get_config(ocnConfigs, 'config_dt', config_dt)
      call mpas_pool_get_config(ocnConfigs, 'config_write_output_on_startup', config_write_output_on_startup)
   
      dt = config_dt

      currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
      call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
      write(stderrUnit,*) 'Initial time ', trim(timeStamp)

      ! fill in diagnostics variables
      call mpas_timer_start("diagnostic solve", .false., initDiagSolveTimer)
      block_ptr => domain % blocklist
      do while(associated(block_ptr))
         call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
         call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)

         call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, 1)
         block_ptr => block_ptr % next
      end do
      call mpas_timer_stop("diagnostic solve", initDiagSolveTimer)

      if (config_write_output_on_startup) then
         call ocn_analysis_compute(domain, err) 
         call mpas_stream_mgr_write(stream_manager, ierr=ierr)
      endif

   end subroutine mpas_core_run!}}}
   
   subroutine mpas_core_finalize(domain, stream_manager)!{{{
   
      use mpas_grid_types
      use mpas_stream_manager
   
      implicit none

      type (domain_type), intent(inout) :: domain 
      type (MPAS_streamManager_type), intent(inout) :: stream_manager
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
   subroutine mpas_core_setup_packages(configPool, packagePool, ierr)!{{{

      use mpas_packages

      implicit none

      type (mpas_pool_type), intent(in) :: configPool
      type (mpas_pool_type), intent(in) :: packagePool

      integer, intent(out) :: ierr

      integer :: err_tmp

      ierr = 0

      call ocn_analysis_setup_packages(configPool, packagePool, err_tmp)
      ierr = ior(ierr, err_tmp)

   end subroutine mpas_core_setup_packages!}}}

   !***********************************************************************
   !
   !  routine mpas_core_setup_clock
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
   subroutine mpas_core_setup_clock(core_clock, configs, ierr)!{{{

      implicit none

      type (MPAS_Clock_type), intent(inout) :: core_clock
      type (mpas_pool_type), intent(inout) :: configs
      integer, intent(out) :: ierr

      call ocn_simulation_clock_init(core_clock, configs, ierr)

   end subroutine mpas_core_setup_clock!}}}

   !***********************************************************************
   !
   !  routine mpas_core_get_mesh_stream
   !
   !> \brief   Returns the name of the stream containing mesh information
   !> \author  Michael Duda
   !> \date    8 August 2014
   !> \details 
   !>  This routine returns the name of the I/O stream containing dimensions,
   !>  attributes, and mesh fields needed by the framework bootstrapping 
   !>  routine. At the time this routine is called, only namelist options 
   !>  are available.
   !
   !-----------------------------------------------------------------------
   subroutine mpas_core_get_mesh_stream(configs, stream, ierr)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: configs
      character(len=*), intent(out) :: stream
      integer, intent(out) :: ierr

      logical, pointer :: config_do_restart

      ierr = 0

      call mpas_pool_get_config(configs, 'config_do_restart', config_do_restart)

      if (.not. associated(config_do_restart)) then
         ierr = 1
      else if (config_do_restart) then
         write(stream,'(a)') 'restart'
      else
         write(stream,'(a)') 'input'
      end if

   end subroutine mpas_core_get_mesh_stream!}}}

end module mpas_core

! vim: foldmethod=marker
