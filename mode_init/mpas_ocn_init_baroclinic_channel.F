! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_baroclinic_channel
!
!> \brief MPAS ocean initialize case -- Baroclinic Channel
!> \author Doug Jacobsen
!> \date   02/18/2014
!> \details
!>  This module contains the routines for initializing the
!>  the baroclinic channel test case
!
!-----------------------------------------------------------------------

module ocn_init_baroclinic_channel

   use mpas_kind_types
   use mpas_io_units
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
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

   public :: ocn_init_setup_baroclinic_channel, &
             ocn_init_validate_baroclinic_channel

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_baroclinic_channel
!
!> \brief   Setup for baroclinic channel test case
!> \author  Doug Jacobsen
!> \date    02/19/2014
!> \details
!>  This routine sets up the initial conditions for the baroclinic channel test case.
!>  It should also ensure the mesh that was input is valid for the configuration.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_setup_baroclinic_channel(domain, iErr)!{{{

   !--------------------------------------------------------------------

      type (domain_type), intent(inout) :: domain
      integer, intent(out) :: iErr
      real (kind=RKIND) :: yMin, yMax, xMin, xMax, dcEdgeMin, dcEdgeMinGlobal
      real (kind=RKIND) :: yMinGlobal, yMaxGlobal, yMidGlobal, xMinGlobal, xMaxGlobal
      real (kind=RKIND) :: temperature, yOffset, xPerturbationMin, xPerturbationMax
      real (kind=RKIND) :: perturbationWidth

      type (block_type), pointer :: block_ptr

      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: verticalMeshPool

      integer :: iCell, k, idx

      ! Define config variable pointers
      character (len=StrKIND), pointer :: config_init_configuration, config_vertical_grid
      logical, pointer :: config_baroclinic_channel_use_distances
      real (kind=RKIND), pointer :: config_baroclinic_channel_gradient_width_dist, &
                                    config_baroclinic_channel_gradient_width_frac, &
                                    config_baroclinic_channel_bottom_depth, &
                                    config_baroclinic_channel_surface_temperature, &
                                    config_baroclinic_channel_bottom_temperature, &
                                    config_baroclinic_channel_temperature_difference, &
                                    config_baroclinic_channel_salinity, &
                                    config_baroclinic_channel_coriolis_parameter

      ! Define dimension pointers
      integer, pointer :: nCellsSolve, nEdgesSolve, nVertLevels, nVertLevelsP1
      integer, pointer :: index_temperature, index_salinity, index_tracer1

      ! Define variable pointers
      integer, dimension(:), pointer :: maxLevelCell
      real (kind=RKIND), dimension(:), pointer :: xCell, yCell,refBottomDepth, refZMid, &
                                                  vertCoordMovementWeights, bottomDepth, &
                                                  fCell, fEdge, fVertex, dcEdge
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers, debugTracers

      ! Define local interfaceLocations variable
      real (kind=RKIND), dimension(:), pointer :: interfaceLocations

      logical, pointer :: on_a_sphere

      iErr = 0

      call mpas_pool_get_config(ocnConfigs, 'config_init_configuration', config_init_configuration)

      if(config_init_configuration .ne. trim('baroclinic_channel')) return

      call mpas_pool_get_config(ocnConfigs, 'config_vertical_grid', config_vertical_grid)

      call mpas_pool_get_config(ocnConfigs, 'config_baroclinic_channel_use_distances', config_baroclinic_channel_use_distances)
      call mpas_pool_get_config(ocnConfigs, 'config_baroclinic_channel_gradient_width_dist', &
                                config_baroclinic_channel_gradient_width_dist)
      call mpas_pool_get_config(ocnConfigs, 'config_baroclinic_channel_gradient_width_frac', &
                                config_baroclinic_channel_gradient_width_frac)
      call mpas_pool_get_config(ocnConfigs, 'config_baroclinic_channel_bottom_depth', config_baroclinic_channel_bottom_depth)
      call mpas_pool_get_config(ocnConfigs, 'config_baroclinic_channel_surface_temperature', &
                                config_baroclinic_channel_surface_temperature)
      call mpas_pool_get_config(ocnConfigs, 'config_baroclinic_channel_bottom_temperature', &
                                config_baroclinic_channel_bottom_temperature)
      call mpas_pool_get_config(ocnConfigs, 'config_baroclinic_channel_temperature_difference', &
                                config_baroclinic_channel_temperature_difference)
      call mpas_pool_get_config(ocnConfigs, 'config_baroclinic_channel_salinity', config_baroclinic_channel_salinity)
      call mpas_pool_get_config(ocnConfigs, 'config_baroclinic_channel_coriolis_parameter', &
                                config_baroclinic_channel_coriolis_parameter)

      ! Determine vertical grid for configuration
      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
      call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
      call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

      if ( on_a_sphere ) call mpas_log_write('The baroclinic channel configuration ' &
              // 'can only be applied to a planar mesh. Exiting...', MPAS_LOG_CRIT)

      allocate(interfaceLocations(nVertLevelsP1))
      call ocn_generate_vertical_grid( config_vertical_grid, interfaceLocations )

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
      end do

      ! Determine global min and max values.
      call mpas_dmpar_min_real(domain % dminfo, yMin, yMinGlobal)
      call mpas_dmpar_max_real(domain % dminfo, yMax, yMaxGlobal)
      call mpas_dmpar_min_real(domain % dminfo, xMin, xMinGlobal)
      call mpas_dmpar_max_real(domain % dminfo, xMax, xMaxGlobal)
      call mpas_dmpar_min_real(domain % dminfo, dcEdgeMin, dcEdgeMinGlobal)

      yMidGlobal = (yMinGlobal + yMaxGlobal) * 0.5_RKIND
      xPerturbationMin = xMinGlobal + 4.0_RKIND * (xMaxGlobal - xMinGlobal) / 6.0_RKIND
      xPerturbationMax = xMinGlobal + 5.0_RKIND * (xMaxGlobal - xMinGlobal) / 6.0_RKIND
      if(config_baroclinic_channel_use_distances) then
         perturbationWidth = config_baroclinic_channel_gradient_width_dist
      else
         perturbationWidth = (yMaxGlobal - yMinGlobal) * config_baroclinic_channel_gradient_width_frac
      end if

      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
        call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)
        call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

        call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
        call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

        call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
        call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)
        call mpas_pool_get_dimension(tracersPool, 'index_tracer1', index_tracer1)

        call mpas_pool_get_array(meshPool, 'xCell', xCell)
        call mpas_pool_get_array(meshPool, 'yCell', yCell)
        call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
        call mpas_pool_get_array(meshPool, 'vertCoordMovementWeights', vertCoordMovementWeights)
        call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
        call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
        call mpas_pool_get_array(meshPool, 'fCell', fCell)
        call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
        call mpas_pool_get_array(meshPool, 'fVertex', fVertex)

        call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
        call mpas_pool_get_array(tracersPool, 'debugTracers', debugTracers, 1)
        call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)

        call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)
        call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)

        call ocn_mark_north_boundary(meshPool, yMaxGlobal, dcEdgeMinGlobal, iErr)
        call ocn_mark_south_boundary(meshPool, yMinGlobal, dcEdgeMinGlobal, iErr)

        ! Set refBottomDepth and refZMid
        do k = 1, nVertLevels
            refBottomDepth(k) = config_baroclinic_channel_bottom_depth * interfaceLocations(k+1)
            refZMid(k) = - 0.5_RKIND * (interfaceLocations(k+1) + interfaceLocations(k)) * config_baroclinic_channel_bottom_depth
        end do

        ! Set vertCoordMovementWeights
        vertCoordMovementWeights(:) = 1.0_RKIND

        do iCell = 1, nCellsSolve
           ! Determine cutoff location for large sin wave
           yOffset = perturbationWidth * sin (6.0_RKIND * pii * (xCell(iCell) - xMinGlobal) / (xMaxGlobal - xMinGlobal))

           ! Set debug tracer
           if ( associated(debugTracers) ) then
              idx = index_tracer1
              do k = 1, nVertLevels
                 debugTracers(idx, k, iCell) = 1.0_RKIND
              enddo
           end if

           ! Set stratification based on northern half of domain temperature
           if ( associated(activeTracers) ) then
              idx = index_temperature
              do k = nVertLevels, 1, -1
                 temperature = config_baroclinic_channel_bottom_temperature &
                             + (config_baroclinic_channel_surface_temperature - config_baroclinic_channel_bottom_temperature) &
                             * ( (refZMid(k) + refBottomDepth(nVertLevels)) / refBottomDepth(nVertLevels) )
                 activeTracers(idx, k, iCell) = temperature
              end do

              if(yCell(iCell) < yMidGlobal - yOffset) then
                   ! If cell is in the southern half, outside the sin width, subtract temperature difference
                   activeTracers(idx, :, iCell) = activeTracers(idx, :, iCell) - config_baroclinic_channel_temperature_difference
              else if(yCell(iCell) >= yMidGlobal - yOffset .and. &
                      yCell(iCell) < yMidGlobal - yOffset + perturbationWidth) then
                   activeTracers(idx, :, iCell) = activeTracers(idx, :, iCell) - config_baroclinic_channel_temperature_difference &
                                                * ( 1.0_RKIND - ( yCell(iCell) - ((yMaxGlobal + yMinGlobal) * 0.5_RKIND &
                                                - yOffset)) / perturbationWidth)
              end if

              ! Determine yOffset for 3rd crest in sin wave.
              yOffset = 0.5_RKIND * perturbationWidth * sin(pii * (xCell(iCell) - xPerturbationMin) &
                      / (xPerturbationMax - xPerturbationMin))

              if ( yCell(iCell) >= yMidGlobal - yOffset - 0.5_RKIND * perturbationWidth .and. &
                   yCell(iCell) <= yMidGlobal - yOffset + 0.5_RKIND * perturbationWidth .and. &
                   xCell(iCell) >= xPerturbationMin .and. &
                   xCell(iCell) <= xPerturbationMax) then


                   do k = 1, nVertLevels
                      activeTracers(idx, k, iCell) = activeTracers(idx, k, iCell) &
                                + 0.3_RKIND * ( 1.0_RKIND - ( ( yCell(iCell) - (yMidGlobal - yOffset)) &
                                / (0.5_RKIND * perturbationWidth)))
                   end do
              end if

              ! Set salinity
              idx = index_salinity
              activeTracers(idx, :, iCell) = config_baroclinic_channel_salinity
           end if

           ! Set layerThickness and restingThickness
           do k = 1, nVertLevels
              layerThickness(k, iCell) = config_baroclinic_channel_bottom_depth * ( interfaceLocations(k+1) &
                                       - interfaceLocations(k) )
              restingThickness(k, iCell) = config_baroclinic_channel_bottom_depth * ( interfaceLocations(k+1) &
                                         - interfaceLocations(k) )

           end do

           ! Set bottomDepth
           bottomDepth(iCell) = config_baroclinic_channel_bottom_depth

           ! Set maxLevelCell
           maxLevelCell(iCell) = nVertLevels
        end do

        ! Set Coriolis parameters
        fCell(:) = config_baroclinic_channel_coriolis_parameter
        fEdge(:) = config_baroclinic_channel_coriolis_parameter
        fVertex(:) = config_baroclinic_channel_coriolis_parameter

        block_ptr => block_ptr % next
      end do

      deallocate(interfaceLocations)

   !--------------------------------------------------------------------

   end subroutine ocn_init_setup_baroclinic_channel!}}}

!***********************************************************************
!
!  routine ocn_init_validate_baroclinic_channel
!
!> \brief   Validation for baroclinic channel test case
!> \author  Doug Jacobsen
!> \date    02/20/2014
!> \details
!>  This routine validates the configuration options for the baroclinic channel test case.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_validate_baroclinic_channel(configPool, packagePool, iocontext, iErr)!{{{

   !--------------------------------------------------------------------
      type (mpas_pool_type), intent(inout) :: configPool, packagePool
      type (mpas_io_context_type), intent(inout) :: iocontext

      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration
      integer, pointer :: config_vert_levels, config_baroclinic_channel_vert_levels

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)

      if(config_init_configuration .ne. trim('baroclinic_channel')) return

      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      call mpas_pool_get_config(configPool, 'config_baroclinic_channel_vert_levels', config_baroclinic_channel_vert_levels)

      if(config_vert_levels <= 0 .and. config_baroclinic_channel_vert_levels > 0) then
         config_vert_levels = config_baroclinic_channel_vert_levels
      else if (config_vert_levels <= 0) then
         call mpas_log_write( 'Validation failed for baroclinic channel. Not given a usable value for vertical levels.', MPAS_LOG_CRIT)
         iErr = 1
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_baroclinic_channel!}}}

!***********************************************************************

end module ocn_init_baroclinic_channel

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
