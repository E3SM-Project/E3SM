! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_ziso
!
!> \brief MPAS ocean initialize case -- Zonally periodic Idealized Southern Ocean (ZISO)
!> \author Phillip J. Wolfram, Luke Van Roekel, Todd Ringler
!> \date   09/14/2015
!> \details
!>  This module contains the routines for initializing the
!>  ZISO initial condition.
!
!-----------------------------------------------------------------------

module ocn_init_ziso

   use mpas_kind_types
   use mpas_io_units
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_stream_manager
   use mpas_dmpar

   use ocn_constants
   use ocn_init_vertical_grids
   use ocn_init_cell_markers

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_init_setup_ziso, &
             ocn_init_validate_ziso

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_ziso
!
!> \brief   Setup for this initial condition
!> \author  Phillip J. Wolfram, Luke Van Roekel, Todd Ringler
!> \date    09/14/2015
!> \details
!>  This routine sets up the initial conditions for the ZISO configuration.
!
!-----------------------------------------------------------------------

  subroutine ocn_init_setup_ziso(domain, iErr)!{{{

    !--------------------------------------------------------------------

    type (domain_type), intent(inout) :: domain
    integer, intent(out) :: iErr

    ! local work variables
    type (block_type), pointer :: block_ptr
    type (mpas_pool_type), pointer :: meshPool, verticalMeshPool, statePool, forcingPool, tracersPool
    type (mpas_pool_type), pointer :: tracersSurfaceRestoringFieldsPool, tracersInteriorRestoringFieldsPool

    integer :: iCell, iEdge, iVertex, k, idx
    real (kind=RKIND), dimension(:), pointer :: interfaceLocations

    ! Define config variable pointers
    character (len=StrKIND), pointer :: config_init_configuration, config_vertical_grid
    logical, pointer :: config_write_cull_cell_mask, config_use_debugTracers

    ! ZISO test case run-time configuration parameters
    logical, pointer :: config_ziso_use_slopping_bathymetry
    real (kind=RKIND), pointer :: config_ziso_meridional_extent
    real (kind=RKIND), pointer :: config_ziso_bottom_depth
    real (kind=RKIND), pointer :: config_ziso_wind_stress_max
    real (kind=RKIND), pointer :: config_ziso_reference_coriolis
    real (kind=RKIND), pointer :: config_ziso_coriolis_gradient
    real (kind=RKIND), pointer :: config_ziso_shelf_depth
    real (kind=RKIND), pointer :: config_ziso_slope_center_position
    real (kind=RKIND), pointer :: config_ziso_slope_half_width
    real (kind=RKIND), pointer :: config_ziso_initial_temp_t1
    real (kind=RKIND), pointer :: config_ziso_initial_temp_t2
    real (kind=RKIND), pointer :: config_ziso_initial_temp_h1
    real (kind=RKIND), pointer :: config_ziso_initial_temp_mt
    real (kind=RKIND), pointer :: config_ziso_mean_restoring_temp
    real (kind=RKIND), pointer :: config_ziso_restoring_temp_dev_ta
    real (kind=RKIND), pointer :: config_ziso_restoring_temp_dev_tb
    real (kind=RKIND), pointer :: config_ziso_restoring_temp_piston_vel
    real (kind=RKIND), pointer :: config_ziso_restoring_sponge_l
    real (kind=RKIND), pointer :: config_ziso_restoring_temp_tau
    real (kind=RKIND), pointer :: config_ziso_restoring_temp_ts
    real (kind=RKIND), pointer :: config_ziso_restoring_temp_ze
    real (kind=RKIND), pointer :: config_ziso_wind_transition_position
    real (kind=RKIND), pointer :: config_ziso_antarctic_shelf_front_width
    real (kind=RKIND), pointer :: config_ziso_wind_stress_shelf_front_max
    logical, pointer :: config_ziso_add_easterly_wind_stress_ASF

    ! configure settings related to frazil
    logical, pointer :: config_ziso_frazil_enable
    real (kind=RKIND), pointer :: config_ziso_frazil_temperature_anomaly

    integer, pointer :: config_ziso_vert_levels

    ! Define dimension pointers
    integer, pointer :: nVertLevels, nCellsSolve, nEdgesSolve, nVerticesSolve, nVertLevelsP1
    integer, pointer :: index_temperature, index_salinity, index_tracer1

    ! Define variable pointers
    logical, pointer :: on_a_sphere
    integer, dimension(:), pointer :: maxLevelCell
    real (kind=RKIND), dimension(:), pointer :: xCell, yCell, xEdge, yEdge, xVertex, yVertex, refBottomDepth, refZMid, &
         vertCoordMovementWeights, bottomDepth, &
         fCell, fEdge, fVertex, dcEdge
    real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness
    real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers, debugTracers
    real (kind=RKIND), dimension(:, :), pointer ::    activeTracersPistonVelocity, activeTracersSurfaceRestoringValue
    real (kind=RKIND), dimension(:, :, :), pointer :: activeTracersInteriorRestoringValue, activeTracersInteriorRestoringRate
    real (kind=RKIND), dimension(:), pointer :: windStressZonal, windStressMeridional

    real (kind=RKIND) :: yMin, yMax, xMin, xMax, dcEdgeMin, dcEdgeMinGlobal
    real (kind=RKIND) :: yMinGlobal, yMaxGlobal, yMidGlobal, xMinGlobal, xMaxGlobal
    real(kind=RKIND), pointer :: y_period
    character (len=StrKIND) :: streamID
    integer :: directionProperty

    ! Local variable related to frazil
    real (kind=RKIND) :: distanceX, distanceY, distance, frazil_temperature, scaleFactor

    ! assume no error
    iErr = 0

    ! test if ZISO is the desired configuration
    call mpas_pool_get_config(ocnConfigs, 'config_init_configuration', config_init_configuration)
    if(config_init_configuration .ne. trim('ziso')) return

    write(stderrUnit,*) 'Starting initialization of Zonally periodic Idealized Southern Ocean (ZISO)'

    ! get config variables !{{{
    call mpas_pool_get_config(domain % configs, 'config_use_debugTracers', config_use_debugTracers)
    call mpas_pool_get_config(domain % configs, 'config_write_cull_cell_mask', config_write_cull_cell_mask)
    call mpas_pool_get_config(domain % configs, 'config_ziso_use_slopping_bathymetry', config_ziso_use_slopping_bathymetry)
    call mpas_pool_get_config(domain % configs, 'config_ziso_bottom_depth', config_ziso_bottom_depth)
    call mpas_pool_get_config(domain % configs, 'config_ziso_meridional_extent', config_ziso_meridional_extent)
    call mpas_pool_get_config(domain % configs, 'config_ziso_reference_coriolis', config_ziso_reference_coriolis)
    call mpas_pool_get_config(domain % configs, 'config_ziso_coriolis_gradient', config_ziso_coriolis_gradient)
    call mpas_pool_get_config(domain % configs, 'config_ziso_vert_levels', config_ziso_vert_levels)
    call mpas_pool_get_config(domain % configs, 'config_ziso_wind_stress_max', config_ziso_wind_stress_max)
    call mpas_pool_get_config(domain % configs, 'config_ziso_slope_half_width', config_ziso_slope_half_width)
    call mpas_pool_get_config(domain % configs, 'config_ziso_shelf_depth', config_ziso_shelf_depth)
    call mpas_pool_get_config(domain % configs, 'config_ziso_slope_center_position', config_ziso_slope_center_position)
    call mpas_pool_get_config(domain % configs, 'config_ziso_initial_temp_t1', config_ziso_initial_temp_t1)
    call mpas_pool_get_config(domain % configs, 'config_ziso_initial_temp_t2', config_ziso_initial_temp_t2)
    call mpas_pool_get_config(domain % configs, 'config_ziso_initial_temp_h1', config_ziso_initial_temp_h1)
    call mpas_pool_get_config(domain % configs, 'config_ziso_initial_temp_mt', config_ziso_initial_temp_mt)
    call mpas_pool_get_config(domain % configs, 'config_ziso_mean_restoring_temp', config_ziso_mean_restoring_temp)
    call mpas_pool_get_config(domain % configs, 'config_ziso_restoring_temp_dev_ta', config_ziso_restoring_temp_dev_ta)
    call mpas_pool_get_config(domain % configs, 'config_ziso_restoring_temp_dev_tb', config_ziso_restoring_temp_dev_tb)
    call mpas_pool_get_config(domain % configs, 'config_ziso_restoring_temp_piston_vel', config_ziso_restoring_temp_piston_vel)
    call mpas_pool_get_config(domain % configs, 'config_ziso_restoring_sponge_l', config_ziso_restoring_sponge_l)
    call mpas_pool_get_config(domain % configs, 'config_ziso_restoring_temp_tau', config_ziso_restoring_temp_tau)
    call mpas_pool_get_config(domain % configs, 'config_ziso_restoring_temp_ts', config_ziso_restoring_temp_ts)
    call mpas_pool_get_config(domain % configs, 'config_ziso_restoring_temp_ze', config_ziso_restoring_temp_ze)
    call mpas_pool_get_config(domain % configs, 'config_vertical_grid', config_vertical_grid)
    call mpas_pool_get_config(domain % configs, 'config_ziso_add_easterly_wind_stress_ASF', &
                              config_ziso_add_easterly_wind_stress_ASF)
    call mpas_pool_get_config(domain % configs, 'config_ziso_wind_transition_position', config_ziso_wind_transition_position)
    call mpas_pool_get_config(domain % configs, 'config_ziso_antarctic_shelf_front_width', config_ziso_antarctic_shelf_front_width)
    call mpas_pool_get_config(domain % configs, 'config_ziso_wind_stress_shelf_front_max', config_ziso_wind_stress_shelf_front_max)

    ! frazil configures
    call mpas_pool_get_config(domain % configs, 'config_ziso_frazil_enable', config_ziso_frazil_enable)
    call mpas_pool_get_config(domain % configs, 'config_ziso_frazil_temperature_anomaly', config_ziso_frazil_temperature_anomaly)
    !}}}

    ! Determine vertical grid for configuration
    call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
    call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
    call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
    call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

    ! test if configure settings are invalid
    if ( on_a_sphere ) call mpas_dmpar_global_abort('MPAS-ocean: ERROR: The ZISO configuration can only be ' &
                                                 // 'applied to a planar mesh. Exiting...')

    ! Define interface locations
    allocate(interfaceLocations(nVertLevelsP1))
    call ocn_generate_vertical_grid( config_vertical_grid, interfaceLocations )

    ! assign config variables
    nVertLevels  = config_ziso_vert_levels
    nVertLevelsP1 = nVertLevels + 1

    ! keep all cells on planar, periodic mesh (no culling)

    !--------------------------------------------------------------------
    ! Use this section to find min/max of grid to allow culling
    !--------------------------------------------------------------------

    ! Initalize min/max values to large positive and negative values
    yMin = 1.0E10_RKIND
    yMax = -1.0E10_RKIND
    xMin = 1.0E10_RKIND
    xMax = -1.0E10_RKIND
    dcEdgeMin = 1.0E10_RKIND

    ! Determine local min and max values.
    block_ptr => domain % blocklist
    do while(associated(block_ptr))
      call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)

      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)

      yMin = min( yMin, minval(yCell(1:nCellsSolve)))
      yMax = max( yMax, maxval(yCell(1:nCellsSolve)))
      xMin = min( xMin, minval(xCell(1:nCellsSolve)))
      xMax = max( xMax, maxval(xCell(1:nCellsSolve)))
      dcEdgeMin = min( dcEdgeMin, minval(dcEdge(1:nEdgesSolve)))

      block_ptr => block_ptr % next
    end do ! do while(associated(block_ptr))


    !--------------------------------------------------------------------
    ! Use this section to set initial values
    !--------------------------------------------------------------------

    block_ptr => domain % blocklist
    do while(associated(block_ptr))
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
       call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)
       call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)

       call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
       call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
       call mpas_pool_get_dimension(meshPool, 'nVerticesSolve', nVerticesSolve)

       call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
       call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)
       call mpas_pool_get_dimension(tracersPool, 'index_tracer1', index_tracer1)

       call mpas_pool_get_array(meshPool, 'xCell', xCell)
       call mpas_pool_get_array(meshPool, 'yCell', yCell)
       call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
       call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
       call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
       call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
       call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
       call mpas_pool_get_array(meshPool, 'vertCoordMovementWeights', vertCoordMovementWeights)
       call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
       call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
       call mpas_pool_get_array(meshPool, 'fCell', fCell)
       call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
       call mpas_pool_get_array(meshPool, 'fVertex', fVertex)

       call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
       if(config_use_debugTracers) call mpas_pool_get_array(tracersPool, 'debugTracers', debugTracers, 1)
       call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)

       call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)
       call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)

       call mpas_pool_get_array(forcingPool, 'windStressZonal', windStressZonal)
       call mpas_pool_get_array(forcingPool, 'windStressMeridional', windStressMeridional)
       ! tests to make sure these are allocated
       if (.not. associated(windStressZonal) .or. .not. associated(windStressMeridional)) then
         call mpas_dmpar_global_abort("MPAS-ocean: windStressZonal and / or windStressMeridional are not allocated")
       end if

       call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceRestoringFields', tracersSurfaceRestoringFieldsPool)
       if (.not. associated(tracersSurfaceRestoringFieldsPool)) then
         call mpas_dmpar_global_abort("MPAS-ocean: tracersSurfaceRestoringFieldsPool not allocated.")
       end if
       call mpas_pool_get_subpool(forcingPool, 'tracersInteriorRestoringFields', tracersInteriorRestoringFieldsPool)
       call mpas_pool_get_array(tracersSurfaceRestoringFieldsPool, 'activeTracersPistonVelocity', activeTracersPistonVelocity, 1)
       if (.not. associated(activeTracersPistonVelocity)) then
         call mpas_dmpar_global_abort("MPAS-ocean: activeTracersPistonVelocity not allocated.")
       end if
       call mpas_pool_get_array(tracersSurfaceRestoringFieldsPool, &
                 'activeTracersSurfaceRestoringValue', activeTracersSurfaceRestoringValue, 1)
       call mpas_pool_get_array(tracersInteriorRestoringFieldsPool, &
                 'activeTracersInteriorRestoringRate', activeTracersInteriorRestoringRate, 1)
       call mpas_pool_get_array(tracersInteriorRestoringFieldsPool, &
                 'activeTracersInteriorRestoringValue', activeTracersInteriorRestoringValue, 1)

       ! Determine global min and max values.
       call mpas_dmpar_min_real(domain % dminfo, yMin, yMinGlobal)
       call mpas_dmpar_max_real(domain % dminfo, yMax, yMaxGlobal)
       call mpas_dmpar_min_real(domain % dminfo, xMin, xMinGlobal)
       call mpas_dmpar_max_real(domain % dminfo, xMax, xMaxGlobal)
       call mpas_dmpar_min_real(domain % dminfo, dcEdgeMin, dcEdgeMinGlobal)

       ! mark north / south boundaries
       if(config_write_cull_cell_mask) then
         call ocn_mark_north_boundary(meshPool, yMaxGlobal, dcEdgeMinGlobal, iErr)
         call ocn_mark_south_boundary(meshPool, yMinGlobal, dcEdgeMinGlobal, iErr)
         call mpas_pool_get_config(meshPool, 'y_period', y_period)
         y_period = 0.0_RKIND
       endif
       call mpas_stream_mgr_begin_iteration(domain % streamManager)
       do while (mpas_stream_mgr_get_next_stream(domain % streamManager, streamID, directionProperty))
         if ( directionProperty == MPAS_STREAM_OUTPUT .or. directionProperty == MPAS_STREAM_INPUT_OUTPUT ) then
           call mpas_stream_mgr_add_att(domain % streamManager, 'y_period', 0.0_RKIND, streamID)
         end if
       end do

       activeTracersInteriorRestoringRate(:,:,:) = 0.0_RKIND
       activeTracersInteriorRestoringValue(:,:,:) = 0.0_RKIND
       activeTracersPistonVelocity(:,:)  = 0.0_RKIND
       activeTracersSurfaceRestoringValue(:,:) = 0.0_RKIND

       ! Set refBottomDepth and refZMid
       do k = 1, nVertLevels
          refBottomDepth(k) = config_ziso_bottom_depth * interfaceLocations(k+1)
          refZMid(k) = - 0.5_RKIND * (interfaceLocations(k+1) + interfaceLocations(k)) * config_ziso_bottom_depth
       end do

       ! set bottomDepth and maxLevelCell !{{{{
       bottomDepth(:) = 0.0_RKIND
       do iCell = 1, nCellsSolve

         if (config_ziso_use_slopping_bathymetry) then
            ! bottom depth function to be applied
            bottomDepth(iCell)  = config_ziso_shelf_depth + &
                                  0.5_RKIND*(config_ziso_bottom_depth - config_ziso_shelf_depth) * &
                                  (1.0_RKIND + tanh((yCell(iCell) - config_ziso_slope_center_position) / &
                                  config_ziso_slope_half_width))
          else
            bottomDepth(iCell)  = config_ziso_bottom_depth
          end if

         ! Determine maxLevelCell based on bottomDepth and refBottomDepth
         ! Also set botomDepth based on refBottomDepth, since
         ! above bottomDepth was set with continuous analytical functions,
         ! and needs to be discrete
         maxLevelCell(iCell) = nVertLevels
         if (nVertLevels > 1) then
           do k = 1, nVertLevels
             if (bottomDepth(iCell) < refBottomDepth(k)) then
               maxLevelCell(iCell) = k-1
               bottomDepth(iCell) = refBottomDepth(k-1)
               exit
             end if
           end do
         end if

       enddo ! Looping through with iCell !}}}

       ! Set vertCoordMovementWeights
       vertCoordMovementWeights(:) = 1.0_RKIND

       do iCell = 1, nCellsSolve

          ! Set initial temperature
          idx = index_temperature
          do k = 1, nVertLevels
             activeTracers(idx, k, iCell) = config_ziso_initial_temp_t1 + &
               config_ziso_initial_temp_t2*tanh(refZMid(k)/config_ziso_initial_temp_h1) + config_ziso_initial_temp_mt*refZMid(k)
          end do

          ! Set initial salinity
          idx = index_salinity
          do k = 1, nVertLevels
             activeTracers(idx, k, iCell) = 34.0_RKIND
          end do

          ! Set layerThickness and restingThickness
          ! Uniform layer thickness
          do k = 1, nVertLevels
            layerThickness(k, iCell) = config_ziso_bottom_depth * ( interfaceLocations(k+1) - interfaceLocations(k) )
            restingThickness(k, iCell) = layerThickness(k, iCell)
          end do

          ! set a passive tracer
          if(config_use_debugTracers) then
             idx = index_tracer1
             do k = 1, nVertLevels
                debugTracers(idx, k, iCell) = 1.0_RKIND + &
                   100000.0_RKIND*exp(-(refZMid(k)+1250.0_RKIND)**2/100.0_RKIND**2) * &
                   exp(-(yCell(iCell)-1250.0_RKIND*1000.0_RKIND)**2/(50.0_RKIND*1000.0_RKIND)**2)
             enddo
          endif

          ! set windstress
          if (config_ziso_add_easterly_wind_stress_ASF) then
            if(yCell(iCell) .ge. config_ziso_wind_transition_position) then
                   windStressZonal(iCell) = config_ziso_wind_stress_max*sin((pii*(yCell(iCell) - &
                                            config_ziso_wind_transition_position) / &
                                            (config_ziso_meridional_extent - config_ziso_wind_transition_position)))**2
            elseif(yCell(iCell) .ge. config_ziso_wind_transition_position - config_ziso_antarctic_shelf_front_width) then
              windStressZonal(iCell) = 0.0_RKIND
              if(yCell(iCell) .lt. config_ziso_wind_transition_position) then
                   windStressZonal(iCell) = config_ziso_wind_stress_shelf_front_max * &
                                            sin((pii*(config_ziso_wind_transition_position &
                                                 - yCell(iCell)))/config_ziso_antarctic_shelf_front_width)**2
              endif
            endif
          else
            windStressZonal(iCell) = config_ziso_wind_stress_max * exp(-((yCell(iCell) - &
                                     config_ziso_meridional_extent/2.0_RKIND) / &
                                     (config_ziso_meridional_extent/2.0_RKIND))**2.0_RKIND) * cos(pii/2.0_RKIND*(yCell(iCell) - &
                                      config_ziso_meridional_extent/2.0_RKIND)/(config_ziso_meridional_extent/2.0_RKIND))
          endif
          windStressMeridional(iCell) = 0.0_RKIND

          ! surface restoring
          idx = index_temperature
          activeTracersSurfaceRestoringValue(idx,iCell) = config_ziso_mean_restoring_temp &
            + config_ziso_restoring_temp_dev_ta * &
              tanh(2.0_RKIND*(yCell(iCell)-config_ziso_meridional_extent/2.0_RKIND)/(config_ziso_meridional_extent/2.0_RKIND)) &
            + config_ziso_restoring_temp_dev_tb * &
              (yCell(iCell)-config_ziso_meridional_extent/2.0_RKIND)/(config_ziso_meridional_extent/2.0_RKIND)
          activeTracersPistonVelocity(idx,iCell) = config_ziso_restoring_temp_piston_vel
          idx = index_salinity
          activeTracersSurfaceRestoringValue(idx,iCell) = 34.0_RKIND
          activeTracersPistonVelocity(idx,iCell) = 0.0_RKIND

          ! set restoring at equatorward (north) boundary
          do k = 1, nVertLevels
            !Interior restoring along northern wall
            if(config_ziso_meridional_extent-yCell(iCell) <= 1.5_RKIND*config_ziso_restoring_sponge_l) then
              idx = index_temperature
              activeTracersInteriorRestoringValue(idx, k, iCell) = activeTracersSurfaceRestoringValue(idx,iCell) &
                                * exp(refZMid(k)/config_ziso_restoring_temp_ze)
              activeTracersInteriorRestoringRate(idx, k, iCell) = &
                                  exp(-(config_ziso_meridional_extent-yCell(iCell))/config_ziso_restoring_sponge_l) &
                                * ( 1.0_RKIND / (config_ziso_restoring_temp_tau*86400.0_RKIND))
              idx = index_salinity
              activeTracersInteriorRestoringValue(idx, k, iCell) = 34.0_RKIND
              activeTracersInteriorRestoringRate(idx, k, iCell) = 0.0_RKIND
            end if
          end do


          ! set restoring at poleward (south) boundary
          do k = 1, nVertLevels
            !Interior restoring along southern wall
            if(yCell(iCell) <= 2.0_RKIND*config_ziso_restoring_sponge_l) then
              idx = index_temperature
              activeTracersInteriorRestoringValue(idx, k, iCell) = activeTracersSurfaceRestoringValue(idx,iCell)
              activeTracersInteriorRestoringRate(idx, k, iCell) =   exp(-yCell(iCell)/config_ziso_restoring_sponge_l) &
                                                                  * ( 1.0_RKIND / (config_ziso_restoring_temp_tau*86400.0_RKIND))
              idx = index_salinity
              activeTracersInteriorRestoringValue(idx, k, iCell) = 34.0_RKIND
              activeTracersInteriorRestoringRate(idx, k, iCell) = 0.0_RKIND
            end if
          enddo

!************************************************************************************************************************
! this test case is overloaded with the ability to evaluate the frazil algorithm
! if config_ziso_enable_frazil is true, some of the configure options are over written to make the test useful for frazil
!************************************************************************************************************************

          if(config_ziso_frazil_enable) then
             config_ziso_initial_temp_t1 = 0.0_RKIND
             config_ziso_initial_temp_t2 = -1.0_RKIND
             config_ziso_initial_temp_h1 = 300.0_RKIND
             config_ziso_initial_temp_mt = 0.0_RKIND

             ! recompute initial temperature with altered parameters
             idx = index_temperature
             do k = 1, nVertLevels
                activeTracers(idx, k, iCell) = config_ziso_initial_temp_t1 + &
                  config_ziso_initial_temp_t2*tanh(refZMid(k)/config_ziso_initial_temp_h1) + config_ziso_initial_temp_mt*refZMid(k)
             end do

             distanceX = config_ziso_meridional_extent/4.0_RKIND-xCell(iCell)
             distanceY = config_ziso_meridional_extent/2.0_RKIND-yCell(iCell)
             distance = sqrt(distanceY**2+distanceX**2)
             scaleFactor = exp(-distance/config_ziso_meridional_extent*20.0_RKIND)
             if (scaleFactor.gt.0.9_RKIND) write(stderrUnit,*) ' frazil production likely at this cell: ', iCell
             do k = 1, nVertLevels
                frazil_temperature = config_ziso_frazil_temperature_anomaly &
                    + config_ziso_initial_temp_t2 * tanh(refZMid(k) / config_ziso_initial_temp_h1) &
                    + config_ziso_initial_temp_mt * refZMid(k)
                if (refZMid(k).gt.-50.0) frazil_temperature = frazil_temperature + 1.0_RKIND*cos( refZMid(k) / 50.0_RKIND &
                                                            * pii / 2.0_RKIND)
                activeTracers(idx, k, iCell) = (1.0_RKIND-scaleFactor)* activeTracers(idx, k, iCell) + scaleFactor &
                                             * frazil_temperature
             end do
          end if

!*********************************************************************************
! end frazil overload
!*********************************************************************************

       end do  ! do iCell

       ! write warning to stderrUnit
       if(config_ziso_frazil_enable) then
         write(stderrUnit,*)
         write(stderrUnit,*) ' this test case is configured for the testing of the frazil algorithm'
         write(stderrUnit,*)
       endif

       ! Set Coriolis parameters, if other than zero
       do iCell = 1, nCellsSolve
          fCell(iCell) = config_ziso_reference_coriolis + yCell(iCell) * config_ziso_coriolis_gradient
       end do
       do iEdge = 1, nEdgesSolve
          fEdge(iEdge) = config_ziso_reference_coriolis + yEdge(iEdge) * config_ziso_coriolis_gradient
       end do
       do iVertex = 1, nVerticesSolve
          fVertex(iVertex) = config_ziso_reference_coriolis + yVertex(iVertex) * config_ziso_reference_coriolis
       end do

       block_ptr => block_ptr % next
    end do  ! do while(associated(block_ptr))

    write(stderrUnit,*) 'Finishing initialization of Zonally periodic Idealized Southern Ocean (ZISO)'
    !--------------------------------------------------------------------

  end subroutine ocn_init_setup_ziso!}}}

!***********************************************************************
!
!  routine ocn_init_validate_ziso
!
!> \brief   Validation for this initial condition
!> \author  Phillip J. Wolfram, Luke Van Roekel, Todd Ringler
!> \date    09/14/2015
!> \details
!>  This routine validates the configuration options for this case.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_validate_ziso(configPool, packagePool, iErr)!{{{

   !--------------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: configPool, packagePool
      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration
      integer, pointer :: config_vert_levels, config_ziso_vert_levels

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)
      if(config_init_configuration .ne. trim('ziso')) return

      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      call mpas_pool_get_config(configPool, 'config_ziso_vert_levels', config_ziso_vert_levels)

      if(config_vert_levels <= 0 .and. config_ziso_vert_levels > 0) then
         config_vert_levels = config_ziso_vert_levels
      else if (config_vert_levels <= 0) then
         write(stderrUnit,*) 'ERROR: Validation failed for ziso. Not given a usable value for vertical levels.'
         iErr = 1
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_ziso!}}}


!***********************************************************************

end module ocn_init_ziso

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
