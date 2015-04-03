! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_cvmix_shear_unit_test
!
!> \brief MPAS ocean initialize case -- CVMix shear Mixing Unit Test
!> \author Doug Jacobsen
!> \date   04/01/2015
!> \details
!>  This module contains the routines for initializing the 
!>  the cvmix shear mixing unit test configuration
!
!-----------------------------------------------------------------------

module ocn_init_cvmix_shear_unit_test

   use mpas_kind_types
   use mpas_io_units
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants

   use ocn_init_cell_markers
   use ocn_init_vertical_grids

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

   public :: ocn_init_setup_cvmix_shear_unit_test, &
             ocn_init_validate_cvmix_shear_unit_test

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_cvmix_shear_unit_test
!
!> \brief   Setup for cvmix shear mixing unit test configuration
!> \author  Doug Jacobsen
!> \date    04/01/2015
!> \details 
!>  This routine sets up the initial conditions for the cvmix shear mixing unit test configuration.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_setup_cvmix_shear_unit_test(domain, iErr)!{{{

   !--------------------------------------------------------------------

      type (domain_type), intent(inout) :: domain
      integer, intent(out) :: iErr
      real (kind=RKIND) :: maxMidDepth, temperature

      type (block_type), pointer :: block_ptr

      type (mpas_pool_type), pointer :: meshPool, verticalMeshPool, statePool
      type (mpas_pool_type), pointer :: diagnosticsPool, forcingPool

      integer, pointer :: nVertLevels, nVertLevelsP1, nCellsSolve, nEdgesSolve
      integer, pointer :: index_temperature, index_salinity

      integer, dimension(:), pointer :: maxLevelCell
      real (kind=RKIND), dimension(:), pointer :: refBottomDepth, refZMid, vertCoordMovementWeights
      real (kind=RKIND), dimension(:), pointer :: surfaceWindStress, boundaryLayerDepth, temperatureRestore
      real (kind=RKIND), dimension(:), pointer :: salinityRestore, bottomDepth, angleEdge
      real (kind=RKIND), dimension(:, :), pointer :: layerThickness, restingThickness
      real (kind=RKIND), dimension(:, :, :), pointer :: tracers

      real (kind=RKIND), dimension(:), pointer :: interfaceLocations

      integer :: iCell, iEdge, k, idx

      character (len=StrKIND), pointer :: config_init_configuration, config_vertical_grid
      real (kind=RKIND), pointer :: config_cvmix_shear_unit_test_bottom_depth, config_cvmix_shear_unit_test_bottom_temperature, &
                                    config_cvmix_shear_unit_test_surface_temperature, config_cvmix_shear_unit_test_salinity, &
                                    config_cvmix_shear_unit_test_max_windstress

      iErr = 0

      call mpas_pool_get_config(domain % configs, 'config_init_configuration', config_init_configuration)

      if(config_init_configuration .ne. trim('cvmix_shear_unit_test')) return


      call mpas_pool_get_config(domain % configs, 'config_cvmix_shear_unit_test_bottom_depth', config_cvmix_shear_unit_test_bottom_depth)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_shear_unit_test_bottom_temperature', config_cvmix_shear_unit_test_bottom_temperature)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_shear_unit_test_surface_temperature', config_cvmix_shear_unit_test_surface_temperature)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_shear_unit_test_salinity', config_cvmix_shear_unit_test_salinity)
      call mpas_pool_get_config(domain % configs, 'config_cvmix_shear_unit_test_max_windstress', config_cvmix_shear_unit_test_max_windstress)

      call mpas_pool_get_config(domain % configs, 'config_vertical_grid', config_vertical_grid)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nVertLevelsP1', nVertLevelsP1)
      allocate(interfaceLocations(nVertLevelsP1))
      call ocn_generate_vertical_grid(config_vertical_grid, interfaceLocations)

      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
        call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)

        call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
        call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
        call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)

        call mpas_pool_get_dimension(statePool, 'index_temperature', index_temperature)
        call mpas_pool_get_dimension(statePool, 'index_salinity', index_salinity)

        call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
        call mpas_pool_get_array(meshPool, 'vertCoordMovementWeights', vertCoordMovementWeights)
        call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
        call mpas_pool_get_array(meshPool, 'temperatureRestore', temperatureRestore)
        call mpas_pool_get_array(meshPool, 'salinityRestore', salinityRestore)
        call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
        call mpas_pool_get_array(meshPool, 'angleEdge', angleEdge)

        call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)
        call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)

        call mpas_pool_get_array(statePool, 'tracers', tracers, 1)
        call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)

        call mpas_pool_get_array(diagnosticsPool, 'boundaryLayerDepth', boundaryLayerDepth)

        call mpas_pool_get_array(forcingPool, 'surfaceWindStress', surfaceWindStress)

        ! Set refBottomDepth and refBottomDepthTopOfCell
        do k = 1, nVertLevels
           refBottomDepth(k) = config_cvmix_shear_unit_test_bottom_depth * interfaceLocations(k+1)
           refZMid(k) = - 0.5_RKIND * config_cvmix_shear_unit_test_bottom_depth * (interfaceLocations(k) + interfaceLocations(k+1))
        end do

        maxMidDepth = -minval(refZMid(:))

        ! Set vertCoordMovementWeights
        vertCoordMovementWeights(:) = 1.0_RKIND

        do iCell = 1, nCellsSolve
           ! Set stratified temperature
           do k = nVertLevels, 1, -1
              temperature = config_cvmix_shear_unit_test_bottom_temperature &
                          + (config_cvmix_shear_unit_test_surface_temperature - config_cvmix_shear_unit_test_bottom_temperature) &
                          * ( (refZMid(k) - refZMid(nVertLevels)) / ( - refZMid(nVertLevels) ))
              tracers(index_temperature, k, iCell) = temperature
           end do

           ! Set salinity
           tracers(index_salinity, :, iCell) = config_cvmix_shear_unit_test_salinity

           ! Set layerThickness
           do k = 1, nVertLevels
              layerThickness(k, iCell) = config_cvmix_shear_unit_test_bottom_depth * (interfaceLocations(k+1) - interfaceLocations(k))
              restingThickness(k, iCell) = layerThickness(k, iCell)
           end do

           ! Set temperatureRestore
           temperatureRestore(iCell) = config_cvmix_shear_unit_test_surface_temperature + 10.0_RKIND

           ! Set salinityRestore
           salinityRestore(iCell) = config_cvmix_shear_unit_test_salinity

           ! Set boundary layer depth
           boundaryLayerDepth(iCell) = 2.0_RKIND * (config_cvmix_shear_unit_test_bottom_depth / nVertLevels) - 1.0-4_RKIND

           ! Set bottomDepth
           bottomDepth(iCell) = config_cvmix_shear_unit_test_bottom_depth

           ! Set maxLevelCell
           maxLevelCell(iCell) = nVertLevels
        end do

        do iEdge = 1, nEdgesSolve
           surfaceWindStress(iEdge) = config_cvmix_shear_unit_test_max_windstress * cos(angleEdge(iEdge))
        end do

        block_ptr => block_ptr % next
      end do

      deallocate(interfaceLocations)

   !--------------------------------------------------------------------

   end subroutine ocn_init_setup_cvmix_shear_unit_test!}}}

!***********************************************************************
!
!  routine ocn_init_validate_cvmix_shear_unit_test
!
!> \brief   Validation for CVMix shear mixing unit test case
!> \author  Doug Jacobsen
!> \date    04/01/2015
!> \details 
!>  This routine validates the configuration options for the CVMix shear mixing unit test configuration.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_validate_cvmix_shear_unit_test(configPool, packagePool, iErr)!{{{

   !--------------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: configPool
      type (mpas_pool_type), intent(in) :: packagePool
      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration
      integer, pointer :: config_vert_levels, config_cvmix_shear_unit_test_vert_levels

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)

      if(config_init_configuration .ne. trim('cvmix_shear_unit_test')) return

      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      call mpas_pool_get_config(configPool, 'config_cvmix_shear_unit_test_vert_levels', config_cvmix_shear_unit_test_vert_levels)

      if(config_vert_levels <= 0 .and. config_cvmix_shear_unit_test_vert_levels > 0) then
         config_vert_levels = config_cvmix_shear_unit_test_vert_levels
      else if(config_vert_levels <= 0) then
         write(stderrUnit,*) 'ERROR: Validation failed for CVMix shear mixing unit test case. Not given a usable value for vertical levels.'
         iErr = 1
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_cvmix_shear_unit_test!}}}

!***********************************************************************

end module ocn_init_cvmix_shear_unit_test

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
