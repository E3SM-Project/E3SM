! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_lock_exchange
!
!> \brief MPAS ocean initialize case -- Lock Exchange
!> \author Doug Jacobsen
!> \date   02/18/2014
!> \details
!>  This module contains the routines for initializing the
!>  the lock exchange test case
!
!-----------------------------------------------------------------------

module ocn_init_lock_exchange

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

   public :: ocn_init_setup_lock_exchange, &
             ocn_init_validate_lock_exchange

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_lock_exchange
!
!> \brief   Setup for lock exchange test case
!> \author  Doug Jacobsen
!> \date    02/18/2014
!> \details
!>  This routine sets up the initial conditions for the lock exchange test case.
!>  It is setup in the y direction, such that everything in the southern half of
!>  the domain has a temperature of 5.0C and the northern half has a value of
!>  30.0C. Salinity is setup as a constant 35PSU.
!>  No windstress is specified, and layerThickness is constant depending on the input parameter
!>  config_lock_exchange_bottom_depth.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_setup_lock_exchange(domain, iErr)!{{{

   !--------------------------------------------------------------------

      type (domain_type), intent(inout) :: domain
      integer, intent(out) :: iErr
      real (kind=RKIND) :: yMin, yMax, xMin, xMax, dcEdgeMin
      real (kind=RKIND) :: yMinGlobal, yMaxGlobal, xMinGlobal, xMaxGlobal, dcEdgeMinGlobal

      type (block_type), pointer :: block_ptr

      type (mpas_pool_type), pointer :: meshPool, statePool, verticalMeshPool, tracersPool

      character (len=StrKIND), pointer :: config_init_configuration, config_vertical_grid, config_lock_exchange_layer_type, &
            config_lock_exchange_direction

      integer, pointer :: nCellsSolve, nEdgesSolve, nVertLevels, nVertLevelsP1, index_temperature, index_salinity, index_tracer1

      integer, dimension(:), pointer :: maxLevelCell

      real (kind=RKIND), pointer :: config_lock_exchange_cold_temperature, config_lock_exchange_warm_temperature, &
                                    config_lock_exchange_salinity, config_lock_exchange_bottom_depth, &
                                    config_lock_exchange_isopycnal_min_thickness

      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, bottomDepth, refBottomDepthTopOfCell, refBottomDepth, &
                                                  vertCoordMovementWeights, dcEdge
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers, debugTracers

      real (kind=RKIND), dimension(:), pointer :: interfaceLocations

      logical, pointer :: on_a_sphere

      integer :: iCell, k

      iErr = 0

      call mpas_pool_get_config(ocnConfigs, 'config_init_configuration', config_init_configuration)

      if (config_init_configuration .ne. trim('lock_exchange')) return

      call mpas_pool_get_config(ocnConfigs, 'config_vertical_grid', config_vertical_grid)

      call mpas_pool_get_config(ocnConfigs, 'config_lock_exchange_cold_temperature', config_lock_exchange_cold_temperature)
      call mpas_pool_get_config(ocnConfigs, 'config_lock_exchange_warm_temperature', config_lock_exchange_warm_temperature)
      call mpas_pool_get_config(ocnConfigs, 'config_lock_exchange_salinity', config_lock_exchange_salinity)
      call mpas_pool_get_config(ocnConfigs, 'config_lock_exchange_bottom_depth', config_lock_exchange_bottom_depth)
      call mpas_pool_get_config(ocnConfigs, 'config_lock_exchange_direction', config_lock_exchange_direction)
      call mpas_pool_get_config(ocnConfigs, 'config_lock_exchange_isopycnal_min_thickness', &
                                config_lock_exchange_isopycnal_min_thickness)
      call mpas_pool_get_config(ocnConfigs, 'config_lock_exchange_layer_type', config_lock_exchange_layer_type)

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
      call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
      call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

      if ( on_a_sphere ) call mpas_log_write('The lock exchange configuration can not be ' &
             // 'applied to spherical meshes', MPAS_LOG_CRIT)

      ! Define interface locations
      allocate( interfaceLocations( nVertLevelsP1 ) )
      call ocn_generate_vertical_grid( config_vertical_grid, interfaceLocations )

      ! Initalize y values to large positive and negative values
      yMin = 1.0E10_RKIND
      yMax = -1.0E10_RKIND
      xMin = 1.0E10_RKIND
      xMax = -1.0E10_RKIND
      dcEdgeMin = 1.0E10_RKIND

      ! Determine local min and max y value.
      block_ptr => domain % blocklist
      do while(associated(block_ptr))
         call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

         call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)

         call mpas_pool_get_array(meshPool, 'xCell', xCell)
         call mpas_pool_get_array(meshPool, 'yCell', yCell)
         call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)

         xMin = min( xMin, minval(xCell(1:nCellssolve)))
         xMax = max( xMax, maxval(xCell(1:nCellssolve)))
         yMin = min( yMin, minval(yCell(1:nCellssolve)))
         yMax = max( yMax, maxval(yCell(1:nCellssolve)))
         dcEdgeMin = min( dcEdgeMin, minval(dcEdge(1:nEdgessolve)))

         block_ptr => block_ptr % next
      end do

      ! Determine global min and max y value. This is so the domain
      ! can be split into north and south.
      call mpas_dmpar_min_real(domain % dminfo, xMin, xMinGlobal)
      call mpas_dmpar_max_real(domain % dminfo, xMax, xMaxGlobal)
      call mpas_dmpar_min_real(domain % dminfo, yMin, yMinGlobal)
      call mpas_dmpar_max_real(domain % dminfo, yMax, yMaxGlobal)
      call mpas_dmpar_min_real(domain % dminfo, dcEdgeMin, dcEdgeMinGlobal)

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

         call mpas_pool_get_array(meshPool, 'yCell', yCell)
         call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
         call mpas_pool_get_array(meshPool, 'refBottomDepthTopOfCell', refBottomDepthTopOfCell)
         call mpas_pool_get_array(meshPool, 'vertCoordMovementWeights', vertCoordMovementWeights)

         call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
         call mpas_pool_get_array(tracersPool, 'debugTracers', debugTracers, 1)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)

         call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)

         call ocn_mark_north_boundary(meshPool, yMaxGlobal, dcEdgeMinGlobal, iErr)
         call ocn_mark_south_boundary(meshPool, yMinGlobal, dcEdgeMinGlobal, iErr)

         do iCell = 1, nCellsSolve
            ! Set layerThickness, and restingThickness

            if ( trim(config_lock_exchange_layer_type) == 'z-level' ) then
               ! Set layerThickness and restingThickness
               do k = 1, nVertLevels
                  layerThickness(k, iCell) = config_lock_exchange_bottom_depth * ( interfaceLocations(k+1) &
                       - interfaceLocations(k) )
                  restingThickness(k, iCell) = layerThickness(k, iCell)
               end do

               ! Set temperature
               if ( associated(activeTracers) ) then
                  if ( trim(config_lock_exchange_direction) == 'x' ) then
                     if(xCell(iCell) < xMinGlobal + (xMaxGlobal - xMinGlobal) * 0.5_RKIND) then
                        activeTracers(index_temperature, :, iCell) = config_lock_exchange_cold_temperature
                     else
                        activeTracers(index_temperature, :, iCell) = config_lock_exchange_warm_temperature
                     end if

                  elseif ( trim(config_lock_exchange_direction) == 'y' ) then
                     if(yCell(iCell) < yMinGlobal + (yMaxGlobal - yMinGlobal) * 0.5_RKIND) then
                        activeTracers(index_temperature, :, iCell) = config_lock_exchange_cold_temperature
                     else
                        activeTracers(index_temperature, :, iCell) = config_lock_exchange_warm_temperature
                     end if

                  elseif ( trim(config_lock_exchange_direction) == 'z' ) then
                        activeTracers(index_temperature, 1:nVertLevels/2, iCell) = config_lock_exchange_warm_temperature
                        activeTracers(index_temperature, nVertLevels/2+1:nVertLevels, iCell) = &
                           config_lock_exchange_cold_temperature
                  else
                     call mpas_log_write('MPAS-ocean: Error: wrong choice of config_lock_exchange_direction')
                  end if
               end if

            else if ( trim(config_lock_exchange_layer_type) == 'isopycnal' ) then
               if ( associated(activeTracers) ) then
                  activeTracers(index_temperature, 1, iCell) = config_lock_exchange_warm_temperature
                  activeTracers(index_temperature, 2:nVertLevels, iCell) = config_lock_exchange_cold_temperature
               end if

               if(yCell(iCell) < (yMaxGlobal - yMinGlobal) * 0.5_RKIND) then
                  layerThickness(1, iCell) = config_lock_exchange_isopycnal_min_thickness
                  layerThickness(2:nVertLevels, iCell) = config_lock_exchange_bottom_depth &
                       - config_lock_exchange_isopycnal_min_thickness
               else
                  layerThickness(1, iCell) = config_lock_exchange_bottom_depth - config_lock_exchange_isopycnal_min_thickness
                  layerThickness(2:nVertLevels, iCell) = config_lock_exchange_isopycnal_min_thickness
               end if
            else
               call mpas_log_write('MPAS-ocean: Error: wrong choice of config_lock_exchange_layer_type')
            end if

            ! Set salinity
            if ( associated(activeTracers) ) then
               activeTracers(index_salinity, :, iCell) = config_lock_exchange_salinity
            end if

            ! Set debugging tracer
            if ( associated(debugTracers) ) then
               do k = 1, nVertLevels
                  debugTracers(index_tracer1, k, iCell) = 1.0_RKIND
               enddo
            end if

            ! Set bottomDepth
            bottomDepth(iCell) = config_lock_exchange_bottom_depth

            ! Set maxLevelCell
            maxLevelCell(iCell) = nVertLevels
         end do

         ! Set refBottomDepth and refBottomDepthTopOfCell
         do k = 1, nVertLevels
            refBottomDepth(k) = config_lock_exchange_bottom_depth * interfaceLocations(k+1)
            refBottomDepthTopOfCell(k) = config_lock_exchange_bottom_depth * interfaceLocations(k)
         end do

         refBottomDepthTopOfCell(nVertLevels+1) = interfaceLocations(nVertLevels+1) * config_lock_exchange_bottom_depth

         ! Set vertCoordMovementWeights
         vertCoordMovementWeights(:) = 1.0_RKIND

         block_ptr => block_ptr % next
      end do

      deallocate(interfaceLocations)

   !--------------------------------------------------------------------

   end subroutine ocn_init_setup_lock_exchange!}}}

!***********************************************************************
!
!  routine ocn_init_validate_lock_exchange
!
!> \brief   Validation for lock exchange test case
!> \author  Doug Jacobsen
!> \date    02/20/2014
!> \details
!>  This routine validates the configuration options for the lock exchange test case.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_validate_lock_exchange(configPool, packagePool, iocontext, iErr)!{{{

   !--------------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: configPool
      type (mpas_pool_type), intent(inout) :: packagePool
      type (mpas_io_context_type), intent(inout) :: iocontext

      integer, intent(out) :: iErr

      integer, pointer :: config_vert_levels, config_lock_exchange_vert_levels
      character (len=StrKIND), pointer :: config_init_configuration

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)

      if(config_init_configuration .ne. trim('lock_exchange')) return


      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      call mpas_pool_get_config(configPool, 'config_lock_exchange_vert_levels', config_lock_exchange_vert_levels)

      if(config_vert_levels <= 0 .and. config_lock_exchange_vert_levels > 0) then
         config_vert_levels = config_lock_exchange_vert_levels
      else if(config_vert_levels <= 0) then
         call mpas_log_write( 'Validation failed for lock exchange test case. Not given a usable value for vertical levels.', MPAS_LOG_CRIT)
         iErr = 1
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_lock_exchange!}}}

!***********************************************************************

end module ocn_init_lock_exchange

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
