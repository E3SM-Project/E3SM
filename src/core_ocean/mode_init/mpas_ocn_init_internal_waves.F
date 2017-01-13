! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_internal_waves
!
!> \brief MPAS ocean initialize case -- Internal waves
!> \author Doug Jacobsen
!> \date   02/18/2014
!> \details
!>  This module contains the routines for initializing the
!>  the internal waves test case
!
!-----------------------------------------------------------------------

module ocn_init_internal_waves

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar

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

   public :: ocn_init_setup_internal_waves, &
             ocn_init_validate_internal_waves

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_internal_waves
!
!> \brief   Setup for internal waves test case
!> \author  Doug Jacobsen
!> \date    02/19/2014
!> \details
!>  This routine sets up the initial conditions for the internal waves test case.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_setup_internal_waves(domain, iErr)!{{{

   !--------------------------------------------------------------------

      type (domain_type), intent(inout) :: domain
      integer, intent(out) :: iErr

      ! Define pool pointers
      type (mpas_pool_type), pointer :: meshPool, verticalMeshPool, statePool, tracersPool

      ! Define dimension pointers
      integer, pointer :: nVertLevels, nVertLevelsP1, nCells, nEdges, nVertices
      integer, pointer :: nCellsSolve, nEdgesSolve, index_temperature, index_salinity, index_tracer1

      ! Define array pointers
      integer, dimension(:), pointer :: maxLevelCell

      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, bottomDepth, dcEdge
      real (kind=RKIND), dimension(:), pointer :: refBottomDepth, refZMid, vertCoordMovementWeights
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers, debugTracers

      real (kind=RKIND) :: yMin, yMax, xMin, xMax, dcEdgeMin
      real (kind=RKIND) :: yMinGlobal, yMaxGlobal, yMidGlobal, xMinGlobal, xMaxGlobal, dcEdgeMinGlobal
      real (kind=RKIND) :: temperature, yOffset, perturbationWidth

      ! Define config pointers
      character (len=StrKIND), pointer :: config_init_configuration, config_vertical_grid, config_internal_waves_layer_type
      logical, pointer :: config_internal_waves_use_distances
      real (kind=RKIND), pointer :: config_internal_waves_amplitude_width_frac, config_internal_waves_amplitude_width_dist
      real (kind=RKIND), pointer :: config_internal_waves_bottom_depth, config_internal_waves_bottom_temperature
      real (kind=RKIND), pointer :: config_internal_waves_surface_temperature, config_internal_waves_temperature_difference
      real (kind=RKIND), pointer :: config_internal_waves_salinity, config_internal_waves_isopycnal_displacement

      type (block_type), pointer :: block_ptr

      integer :: iCell, k

      real (kind=RKIND) :: deltaTemperature
      real (kind=RKIND), dimension(:), pointer :: zTop, refTemperature, refTemperatureTop, refZTop
      real (kind=RKIND), dimension(:), pointer :: interfaceLocations

      iErr = 0

      call mpas_pool_get_config(domain % configs, 'config_init_configuration', config_init_configuration)

      if(config_init_configuration .ne. trim('internal_waves')) return

      ! Initalize min/max values to large positive and negative values
      yMin = 1.0E10_RKIND
      yMax = -1.0E10_RKIND
      xMin = 1.0E10_RKIND
      xMax = -1.0E10_RKIND
      dcEdgEMin = 1.0E10_RKIND

      ! Define locations of layer interfaces
      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
      call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
      call mpas_pool_get_config(domain % configs, 'config_vertical_grid', config_vertical_grid)
      allocate( interfaceLocations( nVertLevelsP1 ) )

      call ocn_generate_vertical_grid( config_vertical_grid, interfaceLocations )

      call mpas_pool_get_config(domain % configs, 'config_internal_waves_use_distances', config_internal_waves_use_distances)
      call mpas_pool_get_config(domain % configs, 'config_internal_waves_amplitude_width_frac', &
                                config_internal_waves_amplitude_width_frac)
      call mpas_pool_get_config(domain % configs, 'config_internal_waves_amplitude_width_dist', &
                                config_internal_waves_amplitude_width_dist)
      call mpas_pool_get_config(domain % configs, 'config_internal_waves_bottom_depth', config_internal_waves_bottom_depth)
      call mpas_pool_get_config(domain % configs, 'config_internal_waves_bottom_temperature', &
                                config_internal_waves_bottom_temperature)
      call mpas_pool_get_config(domain % configs, 'config_internal_waves_surface_temperature', &
                                config_internal_waves_surface_temperature)
      call mpas_pool_get_config(domain % configs, 'config_internal_waves_temperature_difference', &
                                config_internal_waves_temperature_difference)
      call mpas_pool_get_config(domain % configs, 'config_internal_waves_salinity', config_internal_waves_salinity)
      call mpas_pool_get_config(domain % configs, 'config_internal_waves_isopycnal_displacement', &
                                config_internal_waves_isopycnal_displacement)
      call mpas_pool_get_config(domain % configs, 'config_internal_waves_layer_type', config_internal_waves_layer_type)

      ! Determine local min and max values.
      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

        call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
        call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)

        call mpas_pool_get_array(meshPool, 'xCell', xCell)
        call mpas_pool_get_array(meshPool, 'yCell', yCell)
        call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)

        xMin = min( xMin, minval(xCell(1:nCellsSolve)))
        xMax = max( xMax, maxval(xCell(1:nCellsSolve)))
        yMin = min( yMin, minval(yCell(1:nCellsSolve)))
        yMax = max( yMax, maxval(yCell(1:nCellsSolve)))
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
      if(config_internal_waves_use_distances) then
         perturbationWidth = config_internal_waves_amplitude_width_dist
      else
         perturbationWidth = (yMaxGlobal - yMinGlobal) * config_internal_waves_amplitude_width_frac
      end if

      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
        call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

        call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
        call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

        call mpas_pool_get_array(meshPool, 'yCell', yCell)
        call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
        call mpas_pool_get_array(meshPool, 'vertCoordMovementWeights', vertCoordMovementWeights)
        call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
        call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)

        call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)
        call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)

        call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
        call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)
        call mpas_pool_get_dimension(tracersPool, 'index_tracer1', index_tracer1)

        call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
        call mpas_pool_get_array(tracersPool, 'debugTracers', debugTracers, 1)
        call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)

        call ocn_mark_north_boundary(meshPool, yMaxGlobal, dcEdgeMinGlobal,iErr)
        call ocn_mark_south_boundary(meshPool, yMinGlobal, dcEdgeMinGlobal,iErr)

        allocate(zTop(nVertLevels+1), refTemperature(nVertLevels), refTemperatureTop(nVertLevels+1), refZTop(nVertLevels+1))

        ! Set refBottomDepth and refBottomDepthTopOfCell
        do k = 1, nVertLevels
           refBottomDepth(k) = config_internal_waves_bottom_depth * interfaceLocations(k+1)
           refZMid(k) = -0.5_RKIND * config_internal_waves_bottom_depth * (interfaceLocations(k) + interfaceLocations(k+1))
        end do

        if ( trim(config_internal_waves_layer_type) == 'isopycnal' ) then

           refTemperatureTop(1) = config_internal_waves_surface_temperature
           refTemperatureTop(nVertLevels+1) = config_internal_waves_bottom_temperature
           deltaTemperature = (config_internal_waves_surface_temperature - config_internal_waves_bottom_temperature)/nVertLevels
           refTemperature(1) = config_internal_waves_surface_temperature - deltaTemperature/2.0_RKIND
           refZTop(1) = 0.0_RKIND
           do k = 2, nVertLevels
              refTemperatureTop(k) = refTemperatureTop(1) - (k-1)*deltaTemperature
              refTemperature(k) = refTemperature(1) - (k-1)*deltaTemperature
              refZTop(k) = refZTop(k-1) - config_internal_waves_bottom_depth / nVertLevels
           end do

        endif

        ! Set vertCoordMovementWeights
        vertCoordMovementWeights(:) = 1.0_RKIND

        do iCell = 1, nCellsSolve

           ! Set debug tracer
           if ( associated(debugTracers) ) then
              do k = 1, nVertLevels
                debugTracers(index_tracer1, k, iCell) = 1.0_RKIND
              enddo
           end if

           if ( trim(config_internal_waves_layer_type) == 'z-level' ) then

              ! Set stratified temperature
              if ( associated(activeTracers) ) then
                 do k = nVertLevels, 1, -1
                    temperature = config_internal_waves_bottom_temperature &
                         + (config_internal_waves_surface_temperature - config_internal_waves_bottom_temperature) &
                         * ( (refZMid(k) - refZMid(nVertLevels)) / (-refZMid(nVertLevels) ))
                    activeTracers(index_temperature, k, iCell) = temperature
                 end do

                 if ( abs(yCell(iCell) - yMidGlobal) < perturbationWidth ) then
                    ! If cell is in the southern half, outside the sin width, subtract temperature difference
                    do k = 2, nVertLevels
                       temperature = -config_internal_waves_temperature_difference * cos(0.5_RKIND * pii * (yCell(iCell) &
                                   - yMidGlobal) / perturbationWidth) * sin ( pii * refBottomDepth(k-1) &
                                   / refBottomDepth(nVertLevels-1) )

                       activeTracers(index_temperature, k, iCell) = activeTracers(index_temperature, k, iCell) + temperature
                    end do
                 end if
              end if

              ! Set layerThickness and restingThickness
              do k = 1, nVertLevels
                 layerThickness(k, iCell) = config_internal_waves_bottom_depth * ( interfaceLocations(k+1) &
                      - interfaceLocations(k) )
                 restingThickness(k, iCell) = config_internal_waves_bottom_depth * ( interfaceLocations(k+1) &
                      - interfaceLocations(k) )
              end do 

           else if ( trim(config_internal_waves_layer_type) == 'isopycnal' ) then

              ! Set stratified temperature
              if ( associated(activeTracers) ) then
                 activeTracers(index_temperature, :, iCell) =  refTemperature(:)
              end if

              ! Set layerThickness
              if ( abs(yCell(iCell) - yMidGlobal) < perturbationWidth) then
                 ! If cell is in the southern half, outside the sin width, subtract temperature difference
                 zTop(1) = 0.0_RKIND
                 do k = 2, nVertLevels
                    zTop(k) =  refZTop(k) + &
                          config_internal_waves_isopycnal_displacement * sin(pii * (k-1) / (nVertLevels+4)) &
                          * cos(0.5_RKIND * pii * (yCell(iCell) - yMidGlobal) / perturbationWidth)
                 end do
                 zTop(nVertLevels+1) = -config_internal_waves_bottom_depth

                 do k = 1, nVertLevels
                    layerThickness(k, iCell) = zTop(k) - zTop(k+1)
                    restingThickness(k, iCell) = layerThickness(k, iCell)
                 end do
              else
                 layerThickness(:, iCell) = config_internal_waves_bottom_depth / nVertLevels
                 restingThickness(:, iCell) = layerThickness(:, iCell)
              end if
           else
              call mpas_dmpar_global_abort('MPAS-ocean: Error: wrong choice of config_internal_waves_layer_type')
           endif

           ! Set salinity
           if ( associated(activeTracers) ) then
              activeTracers(index_salinity, :, iCell) = config_internal_waves_salinity
           end if

           ! Set bottomDepth
           bottomDepth(iCell) = config_internal_waves_bottom_depth

           ! Set maxLevelCell
           maxLevelCell(iCell) = nVertLevels
        end do

        deallocate(zTop, refTemperature, refTemperatureTop, refZTop)

        block_ptr => block_ptr % next
      end do



      deallocate(interfaceLocations)

   !--------------------------------------------------------------------

   end subroutine ocn_init_setup_internal_waves!}}}

!***********************************************************************
!
!  routine ocn_init_validate_internal_waves
!
!> \brief   Validation for internal waves test case
!> \author  Doug Jacobsen
!> \date    02/20/2014
!> \details
!>  This routine validates the configuration options for the internal waves test case.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_validate_internal_waves(configPool, packagePool, iocontext, iErr)!{{{

   !--------------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: configPool
      type (mpas_pool_type), intent(inout) :: packagePool
      type (mpas_io_context_type), intent(inout) :: iocontext

      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration
      integer, pointer :: config_vert_levels, config_internal_waves_vert_levels

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)

      if(config_init_configuration .ne. trim('internal_waves')) return

      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      call mpas_pool_get_config(configPool, 'config_internal_waves_vert_levels', config_internal_waves_vert_levels)

      if(config_vert_levels <= 0 .and. config_internal_waves_vert_levels > 0) then
         config_vert_levels = config_internal_waves_vert_levels
      else if(config_vert_levels <= 0) then
         write(stderrUnit,*) 'ERROR: Validation failed for internal waves. Not given a usable value for vertical levels.'
         iErr = 1
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_internal_waves!}}}

!***********************************************************************

end module ocn_init_internal_waves

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
