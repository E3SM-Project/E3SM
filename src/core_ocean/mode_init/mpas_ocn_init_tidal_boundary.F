! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_tidal_boundary
!
!> \brief MPAS ocean initialize case -- Tests tidal boundary
!> \author  Phillip Wolfram
!> \date   04/05/2019
!> \details
!>  This module contains the routines for initializing the
!>  the tidal_boundary test case
!
!-----------------------------------------------------------------------

module ocn_init_tidal_boundary

   use mpas_kind_types
   use mpas_io_units
   use mpas_derived_types
   use mpas_pool_routines
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

   public :: ocn_init_setup_tidal_boundary, &
             ocn_init_validate_tidal_boundary

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_tidal_boundary
!
!> \brief   Setup for Tidal Boundary on drying slope
!> \author  Phillip Wolfram
!> \date    04/05/2019
!> \details
!>  This routine sets up the initial conditions for the tidal_boundary test case.
!
!-----------------------------------------------------------------------
  subroutine ocn_init_setup_tidal_boundary(domain, iErr)!{{{

    !--------------------------------------------------------------------

    implicit none

    type (domain_type), intent(inout) :: domain
    integer, intent(out) :: iErr
    real (kind=RKIND) :: yMin, yMax, dcEdgeMin
    real (kind=RKIND) :: yMinGlobal, yMaxGlobal, dcEdgeMinGlobal
    real (kind=RKIND) :: plugWidth
    real (kind=RKIND) :: slopeCenter, slopeWidth, depCenter, depWidth

    type (block_type), pointer :: block_ptr

    type (mpas_pool_type), pointer :: meshPool
    type (mpas_pool_type), pointer :: forcingPool
    type (mpas_pool_type), pointer :: statePool
    type (mpas_pool_type), pointer :: diagnosticsPool
    type (mpas_pool_type), pointer :: verticalMeshPool
    type (mpas_pool_type), pointer :: tracersPool

    integer :: iCell, k

    ! Define dimensions
    integer, pointer :: nCellsSolve, nEdgesSolve, nVertLevels, nVertLevelsP1
    integer, pointer :: index_temperature, index_salinity, index_tracer1

    ! Define arrays
    integer, dimension(:), pointer :: maxLevelCell
    real (kind=RKIND), dimension(:), pointer :: yCell, refBottomDepth, bottomDepth, vertCoordMovementWeights, dcEdge
    real (kind=RKIND), dimension(:), pointer :: tidalInputMask
    real (kind=RKIND), dimension(:), pointer :: ssh
    real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness, zMid
    real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers, debugTracers

    ! Define configs
    character (len=StrKIND), pointer :: config_init_configuration, &
                                        config_vertical_grid, &
                                        config_tidal_boundary_layer_type
    logical, pointer :: config_tidal_boundary_use_distances, config_alter_ICs_for_pbcs
    logical, pointer :: config_use_wetting_drying, config_tidal_start_dry
    real (kind=RKIND), pointer :: config_tidal_boundary_right_bottom_depth, &
                                  config_tidal_boundary_left_bottom_depth, &
                                  config_tidal_boundary_plug_width_dist, config_tidal_boundary_plug_width_frac, &
                                  config_tidal_boundary_plug_temperature, config_tidal_boundary_domain_temperature, &
                                  config_tidal_boundary_salinity, config_tidal_boundary_isopycnal_min_thickness, &
                                  config_tidal_boundary_left_value, &
                                  config_tidal_boundary_right_value, &
                                  config_drying_min_cell_height
    integer, pointer :: config_tidal_boundary_min_vert_levels


    real (kind=RKIND), dimension(:), pointer :: interfaceLocations
    real (kind=RKIND),  parameter :: eps=1.0e-12

    iErr = 0


    call mpas_pool_get_config(ocnConfigs, 'config_init_configuration', config_init_configuration)
    call mpas_pool_get_config(ocnConfigs, 'config_tidal_boundary_layer_type', config_tidal_boundary_layer_type)

    call mpas_log_write('here')
    if(config_init_configuration .ne. trim('tidal_boundary')) return

    call mpas_pool_get_config(ocnConfigs, 'config_vertical_grid', config_vertical_grid)

    call mpas_pool_get_config(ocnConfigs, 'config_tidal_boundary_min_vert_levels', config_tidal_boundary_min_vert_levels)
    call mpas_pool_get_config(ocnConfigs, 'config_tidal_boundary_right_value', config_tidal_boundary_right_value)
    call mpas_pool_get_config(ocnConfigs, 'config_tidal_boundary_left_value', config_tidal_boundary_left_value)
    call mpas_pool_get_config(ocnConfigs, 'config_tidal_boundary_right_bottom_depth', config_tidal_boundary_right_bottom_depth)
    call mpas_pool_get_config(ocnConfigs, 'config_tidal_boundary_left_bottom_depth', config_tidal_boundary_left_bottom_depth)
    call mpas_pool_get_config(ocnConfigs, 'config_tidal_boundary_plug_width_frac', config_tidal_boundary_plug_width_frac)
    call mpas_pool_get_config(ocnConfigs, 'config_tidal_boundary_plug_temperature', config_tidal_boundary_plug_temperature)
    call mpas_pool_get_config(ocnConfigs, 'config_tidal_boundary_domain_temperature', config_tidal_boundary_domain_temperature)
    call mpas_pool_get_config(ocnConfigs, 'config_tidal_boundary_salinity', config_tidal_boundary_salinity)
    call mpas_pool_get_config(ocnConfigs, 'config_alter_ICs_for_pbcs', config_alter_ICs_for_pbcs)

    call mpas_pool_get_config(ocnConfigs, 'config_tidal_start_dry', config_tidal_start_dry)
    call mpas_pool_get_config(ocnConfigs, 'config_use_wetting_drying', config_use_wetting_drying)
    call mpas_pool_get_config(ocnConfigs, 'config_drying_min_cell_height', config_drying_min_cell_height)

    call mpas_pool_get_config(ocnConfigs, 'config_tidal_boundary_use_distances', config_tidal_boundary_use_distances)

    call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
    call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)

    allocate(interfaceLocations(nVertLevelsP1))
    call ocn_generate_vertical_grid(config_vertical_grid, interfaceLocations)

    ! Initalize y values to large positive and negative values
    yMin = 1.0E10_RKIND
    yMax = -1.0E10_RKIND
    dcEdgeMin = 1.0E10_RKIND

    ! Determine local min and max y value.
    block_ptr => domain % blocklist
    do while(associated(block_ptr))
      call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)

      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)

      yMin = min( yMin, minval(yCell(1:nCellssolve)))
      yMax = max( yMax, maxval(yCell(1:nCellssolve)))
      dcEdgeMin = min( dcEdgeMin, minval(dcEdge(1:nEdgesSolve)))

      block_ptr => block_ptr % next
    end do

    ! Determine global min and max y value. This is so the domain
    ! can be split into right and left.
    call mpas_dmpar_min_real(domain % dminfo, yMin, yMinGlobal)
    call mpas_dmpar_max_real(domain % dminfo, yMax, yMaxGlobal)
    call mpas_dmpar_min_real(domain % dminfo, dcEdgeMin, dcEdgeMinGlobal)

    plugWidth = (yMaxGlobal - yMinGlobal) * config_tidal_boundary_plug_width_frac

    block_ptr => domain % blocklist
    do while(associated(block_ptr))
      call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
      call mpas_pool_get_array(forcingPool, 'tidalInputMask', tidalInputMask)

      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
      call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
      call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)
      call mpas_pool_get_dimension(tracersPool, 'index_tracer1', index_tracer1)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
      call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
      call mpas_pool_get_array(meshPool, 'vertCoordMovementWeights', vertCoordMovementWeights)

      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
      call mpas_pool_get_array(tracersPool, 'debugTracers', debugTracers, 1)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
      call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
      call mpas_pool_get_array(statePool, 'ssh', ssh, 1)

      call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)


      if (config_tidal_boundary_right_bottom_depth < config_tidal_boundary_left_bottom_depth) then
        call mpas_log_write('Right boundary must be deeper than left boundary!', MPAS_LOG_CRIT)
      end if

      ! Set refBottomDepth, bottomDepth, and maxLevelCell
      do k = 1, nVertLevels
        refBottomDepth(k) = config_tidal_boundary_right_bottom_depth * interfaceLocations(k+1)
      end do

      if (config_tidal_boundary_use_distances) then
        yMin = config_tidal_boundary_left_value
        yMax = config_tidal_boundary_right_value
      end if

      do iCell = 1, nCellsSolve
        bottomDepth(iCell) = config_tidal_boundary_left_bottom_depth &
          + (yCell(iCell) - yMin) / (yMax - yMin)  * &
          (config_tidal_boundary_right_bottom_depth - config_tidal_boundary_left_bottom_depth)
      end do

      if (config_use_wetting_drying .and. config_tidal_start_dry .and. &
          trim(config_tidal_boundary_layer_type) == 'zstar') then
        do iCell = 1, nCellsSolve
          ! make sure bottomDepth not smaller than 3 layers)
          bottomDepth(iCell) = max(refBottomDepth(2), bottomDepth(iCell))
        end do
      end if

      if (.not. config_alter_ICs_for_pbcs .and. trim(config_tidal_boundary_layer_type) == 'zstar') then
        do iCell = 1, nCellsSolve
          do k = 1,nVertLevels
            if (refBottomDepth(k) > bottomDepth(iCell)) then
              bottomDepth(iCell) = refBottomDepth(k)
              exit
            end if
          end do
        end do
      end if

      if(trim(config_tidal_boundary_layer_type) == 'zstar') then

        ! Determine maxLevelCell based on bottomDepth and refBottomDepth
        ! Also set botomDepth based on refBottomDepth, since
        ! above bottomDepth was set with continuous analytical functions,
        ! and needs to be discrete
        do iCell = 1, nCellsSolve
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

          ssh(iCell) = 0.0_RKIND
        end do

      else if (trim(config_tidal_boundary_layer_type) == 'sigma') then

        do iCell = 1, nCellsSolve
          maxLevelCell(iCell) = nVertLevels
        end do

      else if (trim(config_tidal_boundary_layer_type) == 'hybrid') then

        do iCell = 1, nCellsSolve
          maxLevelCell(iCell) = nVertLevels
          if (nVertLevels > 1) then
            do k = 1, nVertLevels
              if (bottomDepth(iCell) < refBottomDepth(k)) then
                if (k < config_tidal_boundary_min_vert_levels) then
                  ! Use sigma
                  maxLevelCell(iCell) = config_tidal_boundary_min_vert_levels
                  exit
                else
                  ! Use zstar
                  maxLevelCell(iCell) = k-1
                  bottomDepth(iCell) = refBottomDepth(k-1)
                  exit
                end if
              end if
            end do
          end if

          ssh(iCell) = 0.0_RKIND
        end do

      end if


      if (config_use_wetting_drying .and. config_tidal_start_dry) then
        do iCell = 1, nCellsSolve
          ssh(iCell) = -bottomDepth(iCell) + config_drying_min_cell_height*maxLevelCell(iCell)
          ! also computes zMid
          do k = 1, maxLevelCell(iCell)
            layerThickness(k,iCell) = (ssh(iCell) + bottomDepth(iCell))/maxLevelCell(iCell)
            restingThickness(k, iCell) = bottomDepth(iCell)/maxLevelCell(iCell)
          end do
        end do
      else
        do iCell = 1, nCellsSolve
          ! Set layerThickness and restingThickness
          call ocn_compute_layerThickness_zMid_from_bottomDepth(layerThickness(:,iCell), zMid(:,iCell), &
            refBottomDepth, bottomDepth(iCell), maxLevelCell(iCell), nVertLevels, iErr, &
            restingThickness=restingThickness(:,iCell), &
            ssh=ssh(iCell))
          ! also computes zMid
          do k = 1, maxLevelCell(iCell)
            restingThickness(k, iCell) = layerThickness(k, iCell)
          end do
        end do
      end if

      ! Set vertCoordMovementWeights
      vertCoordMovementWeights(:) = 1.0_RKIND

      ! Set tidal boundary mask
      do iCell = 1, nCellsSolve
        tidalInputMask(iCell) = 0.0_RKIND
        if (yCell(iCell) > (25.0e3 - dcEdgeMinGlobal/2.0_RKIND)) then
          tidalInputMask(iCell) = 1.0_RKIND
        ! spread it over multiple cells
        !if (yCell(iCell) > (25.0e3 - 3*dcEdgeMinGlobal)) then
        !  tidalInputMask(iCell) = exp(-((yCell(iCell)-25.0e3)/dcEdgeMinGlobal)**2.0)
        end if
      end do

      ! Set salinity
      if ( associated(activeTracers) ) then
        do iCell = 1, nCellsSolve
          activeTracers(index_salinity, :, iCell) = config_tidal_boundary_salinity
        end do
      end if

      ! Set temperature
      if ( associated(activeTracers) ) then
        do iCell = 1, nCellsSolve
          do k = 1, maxLevelCell(iCell)
            if(yCell(iCell) < yMinGlobal + plugWidth) then
              activeTracers(index_temperature, k, iCell) = config_tidal_boundary_plug_temperature
            else
              activeTracers(index_temperature, k, iCell) = config_tidal_boundary_domain_temperature
            end if
          end do
        end do
      end if

      ! Set debug tracer to match plug onshore
      if ( associated(debugTracers) ) then
        do iCell = 1, nCellsSolve
          do k = 1, maxLevelCell(iCell)
            if(yCell(iCell) < (yMinGlobal + plugWidth)) then
              debugTracers(index_tracer1, k, iCell) = 1.0_RKIND
            else
              debugTracers(index_tracer1, k, iCell) = 0.0_RKIND
            end if
          end do
        end do
      end if

      block_ptr => block_ptr % next
    end do

    deallocate(interfaceLocations)

    !--------------------------------------------------------------------

  end subroutine ocn_init_setup_tidal_boundary!}}}


!***********************************************************************
!
!  routine ocn_init_validate_tidal_boundary
!
!> \brief   Validation for Wetting and drying by tidal boundary test case
!> \author  Phillip Wolfram
!> \date    04/05/2019
!> \details
!>  This routine validates the configuration options for the tidal_boundary test case.
!
!-----------------------------------------------------------------------
   subroutine ocn_init_validate_tidal_boundary(configPool, packagePool, iocontext, iErr)!{{{

   !--------------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: configPool, packagePool
      type (mpas_io_context_type), intent(inout) :: iocontext

      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration
      integer, pointer :: config_tidal_boundary_vert_levels, config_vert_levels

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)

      if(config_init_configuration .ne. trim('tidal_boundary')) return

      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      call mpas_pool_get_config(configPool, 'config_tidal_boundary_vert_levels', config_tidal_boundary_vert_levels)

      if(config_vert_levels <= 0 .and. config_tidal_boundary_vert_levels > 0) then
         config_vert_levels = config_tidal_boundary_vert_levels
         !call mpas_log_write( 'Using value of $i', intArgs=(/ config_vert_levels /))
      else if(config_vert_levels <= 0) then
         call mpas_log_write( 'Validation failed for tidal_boundary test case. ' &
           // 'Not given a usable value for vertical levels.', MPAS_LOG_CRIT)
         iErr = 1
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_tidal_boundary!}}}


!***********************************************************************

end module ocn_init_tidal_boundary

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
