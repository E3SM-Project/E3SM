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

   use mpas_framework
   use mpas_timekeeping
   use mpas_dmpar
   use mpas_timer
   use mpas_io_units

   use ocn_analysis_driver
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
   use ocn_tracer_nonlocalflux
   use ocn_gm

   use ocn_high_freq_thickness_hmix_del2

   use ocn_equation_of_state

   use ocn_vmix

   use ocn_time_average

   use ocn_forcing
   use ocn_sea_ice

   use ocn_constants

   type (MPAS_Clock_type), pointer :: clock

   character (len=*), parameter :: statsAlarmID = 'stats'
   character (len=*), parameter :: coupleAlarmID = 'coupling'

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

   subroutine mpas_core_init(domain, stream_manager, startTimeStamp)!{{{

      use mpas_grid_types
      use mpas_stream_manager
      use ocn_tracer_advection

      implicit none

      type (domain_type), intent(inout) :: domain
      type (MPAS_streamManager_type), intent(inout) :: stream_manager
      character(len=*), intent(out) :: startTimeStamp

      real (kind=RKIND) :: dt
      type (block_type), pointer :: block
      type (dm_info) :: dminfo

      integer :: err, err_tmp
      integer, pointer :: nVertLevels
      real (kind=RKIND) :: maxDensity, maxDensity_global
      real (kind=RKIND), dimension(:), pointer :: meshDensity
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: diagnosticsPool

      character (len=StrKIND), pointer :: xtime
      type (MPAS_Time_Type) :: startTime
      type (MPAS_TimeInterval_type) :: timeStep

      logical, pointer :: config_do_restart, config_filter_btr_mode, config_conduct_tests
      logical, pointer :: config_write_stats_on_startup
      character (len=StrKIND), pointer :: config_vert_coord_movement, config_pressure_gradient_type
      real (kind=RKIND), pointer :: config_maxMeshDensity

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

      dminfo = domain % dminfo

      call mpas_pool_get_config(domain % configs, 'config_do_restart', config_do_restart)
      call mpas_pool_get_config(domain % configs, 'config_vert_coord_movement', config_vert_coord_movement)
      call mpas_pool_get_config(domain % configs, 'config_pressure_gradient_type', config_pressure_gradient_type)
      call mpas_pool_get_config(domain % configs, 'config_filter_btr_mode', config_filter_btr_mode)
      call mpas_pool_get_config(domain % configs, 'config_maxMeshDensity', config_maxMeshDensity)
      call mpas_pool_get_config(domain % configs, 'config_conduct_tests', config_conduct_tests)
      call mpas_pool_get_config(domain % configs, 'config_write_stats_on_startup', config_write_stats_on_startup)

      !
      ! Read input data for model
      !
      if ( config_do_restart ) then
         call MPAS_stream_mgr_read(stream_manager, streamID='modelRestart', ierr=err)
      else
         call MPAS_stream_mgr_read(stream_manager, streamID='input', ierr=err)
      end if
      call mpas_stream_mgr_reset_alarms(stream_manager, streamID='input', ierr=err)
      call MPAS_stream_mgr_reset_alarms(stream_manager, streamID='modelRestart', ierr=err)

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
      call ocn_gm_init(err_tmp)
      err = ior(err,err_tmp)
      call ocn_tracer_nonlocalflux_init(err_tmp)
      err = ior(err,err_tmp)

      call ocn_vmix_init(domain, err_tmp)
      err = ior(err, err_tmp)

      call ocn_equation_of_state_init(err_tmp)
      err = ior(err, err_tmp)

      call ocn_tendency_init(err_tmp)
      err = ior(err,err_tmp)
      call ocn_diagnostics_init(err_tmp)
      err = ior(err,err_tmp)

      call ocn_forcing_init(err_tmp)
      err = ior(err,err_tmp)

      call ocn_high_freq_thickness_hmix_del2_init(err_tmp)
      err = ior(err,err_tmp)

      call ocn_global_diagnostics_init(dminfo,err_tmp)
      err = ior(err, err_tmp)

      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nVertLevels', nVertLevels)
      call ocn_sea_ice_init(nVertLevels, err_tmp)
      err = ior(err, err_tmp)

      call ocn_analysis_init(domain, err_tmp)
      err = ior(err, err_tmp)

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
      timeStep = mpas_get_clock_timestep(clock, ierr=err)
      call mpas_get_timeInterval(timeStep, dt=dt)

      block => domain % blocklist
      do while (associated(block))
         call mpas_init_block(block, dt, err)
         if(err.eq.1) then
             call mpas_dmpar_abort(dminfo)
         endif

         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_array(diagnosticsPool, 'xtime', xtime)
         xtime = startTimeStamp 
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

   subroutine ocn_simulation_clock_init(core_clock, configs, ierr)!{{{

      implicit none

      type (MPAS_Clock_type), intent(inout) :: core_clock
      type (mpas_pool_type), intent(inout) :: configs
      integer, intent(out) :: ierr

      type (MPAS_Time_Type) :: startTime, stopTime, alarmStartTime
      type (MPAS_TimeInterval_type) :: runDuration, timeStep, alarmTimeStep
      character(len=StrKIND) :: restartTimeStamp
      character(len=StrKIND), pointer :: config_start_time, config_stop_time, config_run_duration
      character(len=StrKIND), pointer :: config_stats_interval, config_dt
      integer :: err_tmp


      ierr = 0

      call mpas_pool_get_config(configs, 'config_dt', config_dt)
      call mpas_pool_get_config(configs, 'config_start_time', config_start_time)
      call mpas_pool_get_config(configs, 'config_stop_time', config_stop_time)
      call mpas_pool_get_config(configs, 'config_run_duration', config_run_duration)
      call mpas_pool_get_config(configs, 'config_stats_interval', config_stats_interval)

      call mpas_set_time(curr_time=startTime, dateTimeString=config_start_time, ierr=err_tmp)

      call mpas_set_timeInterval(timeStep, timeString=config_dt, ierr=err_tmp)
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
!> \brief   Initialize blocks within MPAS-Ocean core
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    September 2011
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
      real (kind=RKIND), dimension(:), pointer :: areaCell, boundaryLayerDepth
      real (kind=RKIND), dimension(:,:), pointer :: advCoefs, advCoefs3rd, normalTransportVelocity
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity, normalGMBolusVelocity, edgeTangentVectors
      real (kind=RKIND), dimension(:,:), pointer :: velocityX, velocityY, velocityZ
      real (kind=RKIND), dimension(:,:), pointer :: velocityZonal, velocityMeridional
      real (kind=RKIND), dimension(:,:,:), pointer :: derivTwo

      real (kind=RKIND), dimension(:,:,:), pointer :: tracers

      integer, pointer :: nCells, nEdges, nVertices, nVertLevels
      integer, pointer :: config_horiz_tracer_adv_order
      logical, pointer :: config_hmix_scaleWithMesh, config_do_restart
      logical, pointer :: config_use_standardGM
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

      call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)
      call mpas_pool_get_array(diagnosticsPool, 'normalGMBolusVelocity', normalGMBolusVelocity)
      call mpas_pool_get_array(diagnosticsPool, 'velocityX', velocityX)
      call mpas_pool_get_array(diagnosticsPool, 'velocityY', velocityY)
      call mpas_pool_get_array(diagnosticsPool, 'velocityZ', velocityZ)
      call mpas_pool_get_array(diagnosticsPool, 'velocityZonal', velocityZonal)
      call mpas_pool_get_array(diagnosticsPool, 'velocityMeridional', velocityMeridional)
      call mpas_pool_get_array(diagnosticsPool, 'boundaryLayerDepth', boundaryLayerDepth)

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, 1)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
      call mpas_pool_get_array(statePool, 'tracers', tracers, 1)

      call mpas_pool_get_config(block % configs, 'config_horiz_tracer_adv_order', config_horiz_tracer_adv_order)
      call mpas_pool_get_config(block % configs, 'config_hmix_scaleWithMesh', config_hmix_scaleWithMesh)
      call mpas_pool_get_config(block % configs, 'config_maxMeshDensity', config_maxMeshDensity)
      call mpas_pool_get_config(block % configs, 'config_use_standardGM', config_use_standardGM)
      call mpas_pool_get_config(block % configs, 'config_do_restart', config_do_restart)
   
      call ocn_setup_sign_and_index_fields(meshPool)
      call mpas_initialize_deriv_two(meshPool, derivTwo, err)
      call mpas_tracer_advection_coefficients(meshPool, &
          config_horiz_tracer_adv_order, derivTwo, advCoefs, &
          advCoefs3rd, nAdvCellsForEdge, advCellsForEdge, &
          err1, maxLevelCell, highOrderAdvectionMask, &
          boundaryCell)
      err = ior(err, err1)

      call ocn_time_average_init(averagePool)

      if (.not. config_do_restart) then
         do iCell=1,nCells
            boundaryLayerDepth(iCell) = layerThickness(1, iCell) * 0.5
         end do
      end if

      call mpas_timer_start("diagnostic solve", .false., initDiagSolveTimer)
      call ocn_diagnostic_solve(dt,  statePool, forcingPool, meshPool, diagnosticsPool, scratchPool)
      call mpas_timer_stop("diagnostic solve", initDiagSolveTimer)

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

      ! ------------------------------------------------------------------
      ! Accumulating various parametrizations of the transport velocity
      ! ------------------------------------------------------------------
      normalTransportVelocity(:,:) = normalVelocity(:,:)

      ! Compute normalGMBolusVelocity, relativeSlope and RediDiffVertCoef if respective flags are turned on
      if (config_use_standardGM) then
          call ocn_gm_compute_Bolus_velocity(diagnosticsPool, meshPool, scratchPool)
      end if
 
      if (config_use_standardGM) then
         normalTransportVelocity(:,:) = normalTransportVelocity(:,:) + normalGMBolusVelocity(:,:)
       end if
      ! ------------------------------------------------------------------
      ! End: Accumulating various parametrizations of the transport velocity
      ! ------------------------------------------------------------------

      call ocn_compute_mesh_scaling(meshPool, config_hmix_scaleWithMesh, config_maxMeshDensity)
 
      call mpas_rbf_interp_initialize(meshPool)
      call mpas_initialize_tangent_vectors(meshPool, edgeTangentVectors)

      call mpas_init_reconstruct(meshPool)

      call mpas_reconstruct(meshPool, normalVelocity,        &
                       velocityX,            &
                       velocityY,            &
                       velocityZ,            &
                       velocityZonal,        &
                       velocityMeridional    &
                      )

      if (config_use_standardGM) then
         call ocn_reconstruct_gm_vectors(diagnosticsPool, meshPool)
      end if

      call mpas_pool_initialize_time_levels(statePool)

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

   subroutine mpas_core_run(domain, stream_manager)!{{{
   
      use mpas_kind_types
      use mpas_grid_types
      use mpas_stream_manager
      use mpas_timer
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      type (MPAS_streamManager_type), intent(inout) :: stream_manager
   
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
      logical, pointer :: config_write_output_on_startup
   
      call mpas_pool_get_config(domain % configs, 'config_write_output_on_startup', config_write_output_on_startup)

      ! Eventually, dt should be domain specific
      timeStep = mpas_get_clock_timestep(clock, ierr=ierr)
      call mpas_get_timeInterval(timeStep, dt=dt)

      currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
      call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
      write(stderrUnit,*) 'Initial time ', trim(timeStamp)

      call ocn_analysis_compute_startup(domain, stream_manager, err) 

      if (config_write_output_on_startup) then
          call mpas_stream_mgr_write(stream_manager, 'output', forceWriteNow=.true., ierr=ierr)
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
      do while (.not. mpas_is_clock_stop_time(clock))
         call mpas_stream_mgr_read(stream_manager, ierr=ierr)
         call mpas_stream_mgr_reset_alarms(stream_manager, direction=MPAS_STREAM_INPUT, ierr=ierr)

         itimestep = itimestep + 1
         call mpas_advance_clock(clock)

         currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)
         call mpas_get_time(curr_time=currTime, dateTimeString=timeStamp, ierr=ierr)
         write(stderrUnit,*) 'Doing timestep ', trim(timeStamp)
   
         block_ptr => domain % blocklist
         do while(associated(block_ptr))
           call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
           call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
           call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
           call ocn_forcing_build_arrays(meshPool, statePool, forcingPool, ierr, 1)
           call ocn_forcing_build_transmission_array(meshPool, statePool, forcingpool, ierr, 1)
           block_ptr => block_ptr % next
         end do

         call mpas_timer_start("time integration", .false., timeIntTimer)
         call mpas_timestep(domain, itimestep, dt, timeStamp)
         call mpas_timer_stop("time integration", timeIntTimer)
   
         ! Move time level 2 fields back into time level 1 for next time step
         block_ptr => domain % blocklist
         do while(associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
            call mpas_pool_shift_time_levels(statePool)
            block_ptr => block_ptr % next
         end do
      
         call ocn_analysis_compute_w_alarms(stream_manager, domain, err) 
         call ocn_analysis_write(stream_manager, err)

         call mpas_stream_mgr_write(stream_manager, streamID='output', ierr=ierr)
         call mpas_stream_mgr_reset_alarms(stream_manager, streamID='output', ierr=ierr)

         call mpas_stream_mgr_write(stream_manager, streamID='modelCheckpoint', ierr=ierr)
         call mpas_stream_mgr_reset_alarms(stream_manager, streamID='modelCheckpoint', ierr=ierr)

         call mpas_stream_mgr_write(stream_manager, ierr=ierr)
         call mpas_stream_mgr_reset_alarms(stream_manager, direction=MPAS_STREAM_OUTPUT, ierr=ierr)
      end do

   end subroutine mpas_core_run!}}}
   
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

      type (mpas_pool_type), pointer :: statePool, meshPool

      integer :: iTracer, cell, cell1, cell2
      integer, dimension(:), pointer :: maxLevelEdgeTop
      integer, dimension(:,:), pointer :: cellsOnEdge
      real (kind=RKIND) :: normalThicknessFluxSum, layerThicknessSum, layerThicknessEdge1
      real (kind=RKIND), dimension(:), pointer :: refBottomDepth, normalBarotropicVelocity
         
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness
      real (kind=RKIND), dimension(:,:), pointer :: normalBaroclinicVelocity, normalVelocity
      integer, pointer :: nVertLevels, nCells, nEdges
      character (len=StrKIND), pointer :: config_time_integrator
      logical, pointer :: config_filter_btr_mode

      ! Initialize z-level mesh variables from h, read in from input file.
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_config(block % configs, 'config_time_integrator', config_time_integrator)
         call mpas_pool_get_config(block % configs, 'config_filter_btr_mode', config_filter_btr_mode)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
         call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)

         call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, 1)
         call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocity, 1)
         call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', normalBaroclinicVelocity, 1)
         
         call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
         call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
         call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)

         ! Compute barotropic velocity at first timestep
         ! This is only done upon start-up.
         if (trim(config_time_integrator) == 'unsplit_explicit') then
            call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', normalBarotropicVelocity)
            normalBarotropicVelocity(:) = 0.0

            normalBaroclinicVelocity(:,:) = normalVelocity(:,:) 

         elseif (trim(config_time_integrator) == 'split_explicit') then

            if (config_filter_btr_mode) then
               do iCell = 1, nCells
                  layerThickness(1,iCell) = refBottomDepth(1)
               enddo
            endif 

            do iEdge = 1, nEdges
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)

               ! normalBarotropicVelocity = sum(u)/sum(h) on each column
               ! ocn_diagnostic_solve has not yet been called, so compute hEdge 
               ! just for this edge.

               ! thicknessSum is initialized outside the loop because on land boundaries 
               ! maxLevelEdgeTop=0, but I want to initialize thicknessSum with a 
               ! nonzero value to avoid a NaN.
               layerThicknessEdge1 = 0.5*( layerThickness(1,cell1) + layerThickness(1,cell2) ) 
               normalThicknessFluxSum = layerThicknessEdge1 * normalVelocity(1,iEdge)
               layerThicknessSum = layerThicknessEdge1

               do k=2, maxLevelEdgeTop(iEdge)
                  ! ocn_diagnostic_solve has not yet been called, so compute hEdge 
                  ! just for this edge.
                  layerThicknessEdge1 = 0.5*( layerThickness(k,cell1) + layerThickness(k,cell2) ) 

                  normalThicknessFluxSum = normalThicknessFluxSum &
                     + layerThicknessEdge1 * normalVelocity(k,iEdge)
                  layerThicknessSum = layerThicknessSum + layerThicknessEdge1

               enddo
               normalBarotropicVelocity(iEdge) = normalThicknessFluxSum / layerThicknessSum

               ! normalBaroclinicVelocity(k,iEdge) = normalVelocity(k,iEdge) - normalBarotropicVelocity(iEdge)
               do k = 1, maxLevelEdgeTop(iEdge)
                  normalBaroclinicVelocity(k,iEdge) = normalVelocity(k,iEdge) - normalBarotropicVelocity(iEdge)
               enddo

               ! normalBaroclinicVelocity=0, normalVelocity=0 on land cells
               do k = maxLevelEdgeTop(iEdge)+1, nVertLevels
                  normalBaroclinicVelocity(k,iEdge) = 0.0
                  normalVelocity(k,iEdge) = 0.0
               enddo
            enddo

            if (config_filter_btr_mode) then
               ! filter normalBarotropicVelocity out of initial condition
                normalVelocity(:,:) = normalBaroclinicVelocity(:,:)

                normalBarotropicVelocity(:) = 0.0
            endif 

         endif

      block => block % next
      end do

   end subroutine ocn_init_split_timestep!}}}
 
   subroutine mpas_core_finalize(domain, stream_manager)!{{{
   
      use mpas_grid_types
      use mpas_stream_manager
   
      implicit none

      type (domain_type), intent(inout) :: domain 
      type (MPAS_streamManager_type), intent(inout) :: stream_manager
      integer :: ierr

      call ocn_analysis_finalize(domain, ierr) 

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
   subroutine mpas_core_setup_packages(configPool, packagePool, ierr)!{{{

      implicit none

      type (mpas_pool_type), intent(in) :: configPool
      type (mpas_pool_type), intent(in) :: packagePool

      integer, intent(out) :: ierr

      integer :: err_tmp

      logical, pointer :: thicknessFilterActive
      logical, pointer :: splitTimeIntegratorActive
      logical, pointer :: surfaceRestoringActive
      logical, pointer :: bulkForcingActive
      logical, pointer :: frazilIceActive
      logical, pointer :: inSituEOSActive

      logical, pointer :: config_use_freq_filtered_thickness
      logical, pointer :: config_frazil_ice_formation
      character (len=StrKIND), pointer :: config_time_integrator, config_forcing_type, config_pressure_gradient_type

      call mpas_pool_get_package(packagePool, 'thicknessFilterActive', thicknessFilterActive)
      call mpas_pool_get_package(packagePool, 'splitTimeIntegratorActive', splitTimeIntegratorActive)
      call mpas_pool_get_package(packagePool, 'surfaceRestoringActive', surfaceRestoringActive)
      call mpas_pool_get_package(packagePool, 'bulkForcingActive', bulkForcingActive)
      call mpas_pool_get_package(packagePool, 'frazilIceActive', frazilIceActive)
      call mpas_pool_get_package(packagePool, 'inSituEOSActive', inSituEOSActive)

      call mpas_pool_get_config(configPool, 'config_use_freq_filtered_thickness', config_use_freq_filtered_thickness)
      call mpas_pool_get_config(configPool, 'config_time_integrator', config_time_integrator)
      call mpas_pool_get_config(configPool, 'config_forcing_type', config_forcing_type)
      call mpas_pool_get_config(configPool, 'config_frazil_ice_formation', config_frazil_ice_formation)
      call mpas_pool_get_config(configPool, 'config_pressure_gradient_type', config_pressure_gradient_type)

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

      if (config_pressure_gradient_type.eq.'common_level_eos') then
         inSituEOSActive = .true.
      end if

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
         write(stream,'(a)') 'modelRestart'
      else
         write(stream,'(a)') 'input'
      end if

   end subroutine mpas_core_get_mesh_stream!}}}

end module mpas_core

! vim: foldmethod=marker
