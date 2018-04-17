! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_overflow
!
!> \brief MPAS ocean initialize case -- Overflow
!> \author Doug Jacobsen
!> \date   02/18/2014
!> \details
!>  This module contains the routines for initializing the
!>  the overflow test case
!
!-----------------------------------------------------------------------

module ocn_init_overflow

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

   public :: ocn_init_setup_overflow, &
             ocn_init_validate_overflow

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_overflow
!
!> \brief   Setup for overflow test case
!> \author  Doug Jacobsen
!> \date    02/18/2014
!> \details
!>  This routine sets up the initial conditions for the overflow test case.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_setup_overflow(domain, iErr)!{{{

   !--------------------------------------------------------------------

      type (domain_type), intent(inout) :: domain
      integer, intent(out) :: iErr
      real (kind=RKIND) :: yMin, yMax, dcEdgeMin
      real (kind=RKIND) :: yMinGlobal, yMaxGlobal, dcEdgeMinGlobal
      real (kind=RKIND) :: plugWidth
      real (kind=RKIND) :: slopeCenter, slopeWidth

      type (block_type), pointer :: block_ptr

      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: verticalMeshPool
      type (mpas_pool_type), pointer :: tracersPool

      integer :: iCell, k

      ! Define dimensions
      integer, pointer :: nCellsSolve, nEdgesSolve, nVertLevels, nVertLevelsP1
      integer, pointer :: index_temperature, index_salinity, index_tracer1

      ! Define arrays
      integer, dimension(:), pointer :: maxLevelCell
      real (kind=RKIND), dimension(:), pointer :: yCell, refBottomDepth, bottomDepth, vertCoordMovementWeights, dcEdge
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers, debugTracers

      ! Define configs
      character (len=StrKIND), pointer :: config_init_configuration, config_vertical_grid, config_overflow_layer_type
      logical, pointer :: config_overflow_use_distances
      real (kind=RKIND), pointer :: config_overflow_plug_width_dist, config_overflow_slope_center_dist, &
                                    config_overflow_slope_width_dist, config_overflow_plug_width_frac, &
                                    config_overflow_slope_center_frac, config_overflow_slope_width_frac, &
                                    config_overflow_bottom_depth, config_overflow_ridge_depth, &
                                    config_overflow_plug_temperature, config_overflow_domain_temperature, &
                                    config_overflow_salinity, config_overflow_isopycnal_min_thickness


     real (kind=RKIND), dimension(:), pointer :: interfaceLocations

      iErr = 0


      call mpas_pool_get_config(ocnConfigs, 'config_init_configuration', config_init_configuration)

      if(config_init_configuration .ne. trim('overflow')) return

      call mpas_pool_get_config(ocnConfigs, 'config_vertical_grid', config_vertical_grid)

      call mpas_pool_get_config(ocnConfigs, 'config_overflow_use_distances', config_overflow_use_distances)
      call mpas_pool_get_config(ocnConfigs, 'config_overflow_plug_width_dist', config_overflow_plug_width_dist)
      call mpas_pool_get_config(ocnConfigs, 'config_overflow_slope_center_dist', config_overflow_slope_center_dist)
      call mpas_pool_get_config(ocnConfigs, 'config_overflow_slope_width_dist', config_overflow_slope_width_dist)
      call mpas_pool_get_config(ocnConfigs, 'config_overflow_plug_width_frac', config_overflow_plug_width_frac)
      call mpas_pool_get_config(ocnConfigs, 'config_overflow_slope_center_frac', config_overflow_slope_center_frac)
      call mpas_pool_get_config(ocnConfigs, 'config_overflow_slope_width_frac', config_overflow_slope_width_frac)
      call mpas_pool_get_config(ocnConfigs, 'config_overflow_bottom_depth', config_overflow_bottom_depth)
      call mpas_pool_get_config(ocnConfigs, 'config_overflow_ridge_depth', config_overflow_ridge_depth)
      call mpas_pool_get_config(ocnConfigs, 'config_overflow_plug_temperature', config_overflow_plug_temperature)
      call mpas_pool_get_config(ocnConfigs, 'config_overflow_domain_temperature', config_overflow_domain_temperature)
      call mpas_pool_get_config(ocnConfigs, 'config_overflow_salinity', config_overflow_salinity)
      call mpas_pool_get_config(ocnConfigs, 'config_overflow_layer_type', config_overflow_layer_type)
      call mpas_pool_get_config(ocnConfigs, 'config_overflow_isopycnal_min_thickness', config_overflow_isopycnal_min_thickness)

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
        dcEdgeMin = min( dcEdgeMin, minval(dcEdge(1:nEdgessolve)))

        block_ptr => block_ptr % next
      end do

      ! Determine global min and max y value. This is so the domain
      ! can be split into north and south.
      call mpas_dmpar_min_real(domain % dminfo, yMin, yMinGlobal)
      call mpas_dmpar_max_real(domain % dminfo, yMax, yMaxGlobal)
      call mpas_dmpar_min_real(domain % dminfo, dcEdgeMin, dcEdgeMinGlobal)

      if ( config_overflow_use_distances ) then
         plugWidth = config_overflow_plug_width_dist
         slopeCenter = yMinGlobal + config_overflow_slope_center_dist
         slopeWidth = config_overflow_slope_width_dist
      else
         plugWidth = (yMaxGlobal - yMinGlobal) * config_overflow_plug_width_frac
         slopeCenter = yMinGlobal + (yMaxGlobal - yMinGlobal) * config_overflow_slope_center_frac
         slopeWidth = (yMaxGlobal - yMinGlobal) * config_overflow_slope_width_frac
      end if

      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
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

        call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)

        call ocn_mark_north_boundary(meshPool, yMaxGlobal, dcEdgeMinGlobal, iErr)
        call ocn_mark_south_boundary(meshPool, yMinGlobal, dcEdgeMinGlobal, iErr)

        ! Set refBottomDepth, bottomDepth, and maxLevelCell
        do k = 1, nVertLevels
           refBottomDepth(k) = config_overflow_bottom_depth * interfaceLocations(k+1)
        end do

        do iCell = 1, nCellsSolve
           ! From Mehmet Ilicak:
           ! depth=2000
           ! val1 = 500 is top of ridge
           ! h(i,j) = val1 + 0.5*(depth-val1) * (1.0+TANH((lon(i,j)-40000.0)/7000.0))
           bottomDepth(iCell) = config_overflow_ridge_depth &
                     + 0.5_RKIND*(config_overflow_bottom_depth - config_overflow_ridge_depth) &
                     * (1.0_RKIND+tanh((yCell(iCell) - slopeCenter)/slopeWidth))

           if ( trim(config_overflow_layer_type) == 'sigma' .or. trim(config_overflow_layer_type) == 'isopycnal' ) then
              maxLevelCell(iCell) = nVertLevels
           else if ( trim(config_overflow_layer_type) == 'z-level' ) then
              maxLevelCell(iCell) = -1
              do k = 1, nVertLevels
                 if (bottomDepth(iCell) .le. refBottomDepth(k) .and. &
                     maxLevelCell(iCell) == -1) then

                     maxLevelCell(iCell) = k
                 end if
              end do
           end if
        end do

        do iCell = 1, nCellsSolve
           ! Set temperature
           if ( associated(activeTracers) ) then
              if ( trim(config_overflow_layer_type) == 'sigma' .or. trim(config_overflow_layer_type) == 'z-level' ) then
                 do k = 1, maxLevelCell(iCell)
                    if(yCell(iCell) < yMinGlobal + plugWidth) then
                       activeTracers(index_temperature, k, iCell) = config_overflow_plug_temperature
                    else
                       activeTracers(index_temperature, k, iCell) = config_overflow_domain_temperature
                    end if
                 end do
              else if ( trim(config_overflow_layer_type) == 'isopycnal' ) then
                 activeTracers(index_temperature, 1, :) = config_overflow_domain_temperature
                 activeTracers(index_temperature, 2:nVertLevels, :) = config_overflow_plug_temperature
              end if
           end if

           ! Set layerThickness and restingThickness
           if ( trim(config_overflow_layer_type) == 'z-level' ) then
              do k = 1, maxLevelCell(iCell)
                 layerThickness(k, iCell) = config_overflow_bottom_depth * (interfaceLocations(k+1) - interfaceLocations(k))
                 restingThickness(k, iCell) = layerThickness(k, iCell)
              end do
           else if ( trim(config_overflow_layer_type) == 'sigma' ) then
              do k = 1, nVertLevels
                 layerThickness(k, iCell) = bottomDepth(iCell) / nVertLevels
                 restingThickness(k, iCell) = layerThickness(k, iCell)
              end do
           else if ( trim(config_overflow_layer_type) == 'isopycnal' ) then
              ! Set layerThickness.  Normally isopycnal overflow has only two layers.
              if ( yCell(iCell) < yMinGlobal + plugWidth) then
                 layerThickness(1, iCell) = config_overflow_isopycnal_min_thickness
                 layerThickness(2:nVertLevels, iCell) = bottomDepth(iCell) - config_overflow_isopycnal_min_thickness
                 restingThickness(:, iCell) = layerThickness(:, iCell)
              else
                 layerThickness(1, iCell) = bottomDepth(iCell) - config_overflow_isopycnal_min_thickness
                 layerThickness(2:nVertLevels, iCell) = config_overflow_isopycnal_min_thickness
                 restingThickness(:, iCell) = layerThickness(:, iCell)
              end if
           end if

           ! Set salinity
           if ( associated(activeTracers) ) then
              activeTracers(index_salinity, :, iCell) = config_overflow_salinity
           end if

           ! Set debug tracer
           if ( associated(debugTracers) ) then
              do k = 1, nVertLevels
                debugTracers(index_tracer1, k, iCell) = 1.0_RKIND
              end do
           end if

        end do

        ! Set vertCoordMovementWeights
        vertCoordMovementWeights(:) = 1.0_RKIND

        block_ptr => block_ptr % next
      end do

      deallocate(interfaceLocations)

   !--------------------------------------------------------------------

   end subroutine ocn_init_setup_overflow!}}}

!***********************************************************************
!
!  routine ocn_init_validate_overflow
!
!> \brief   Validation for overflow test case
!> \author  Doug Jacobsen
!> \date    02/20/2014
!> \details
!>  This routine validates the configuration options for the overflow test case.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_validate_overflow(configPool, packagePool, iocontext, iErr)!{{{

   !--------------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: configPool, packagePool
      type (mpas_io_context_type), intent(inout) :: iocontext

      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration
      integer, pointer :: config_overflow_vert_levels, config_vert_levels

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)

      if(config_init_configuration .ne. trim('overflow')) return

      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      call mpas_pool_get_config(configPool, 'config_overflow_vert_levels', config_overflow_vert_levels)

      if(config_vert_levels <= 0 .and. config_overflow_vert_levels > 0) then
         config_vert_levels = config_overflow_vert_levels
      else if(config_vert_levels <= 0) then
         call mpas_log_write( 'Validation failed for overflow test case. Not given a usable value for vertical levels.', MPAS_LOG_CRIT)
         iErr = 1
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_overflow!}}}

!***********************************************************************

end module ocn_init_overflow

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
