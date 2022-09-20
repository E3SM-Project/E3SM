! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_mixed_layer_eddy
!
!> \brief MPAS ocean initialize case -- Mixed Layer Eddy
!> \author Qing Li
!> \date   01/17/2020
!> \details
!>  This module contains the routines for initializing the
!>  the mixed layer eddy test case
!
!-----------------------------------------------------------------------

module ocn_init_mixed_layer_eddy

   use mpas_kind_types
   use mpas_io_units
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar

   use ocn_constants
   use ocn_config
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

   public :: ocn_init_setup_mixed_layer_eddy, &
             ocn_init_validate_mixed_layer_eddy

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_mixed_layer_eddy
!
!> \brief   Setup for mixed layer eddy test case
!> \author  Qing Li
!> \date    01/17/2020
!> \details
!>  This routine sets up the initial conditions for the mixed layer eddy test case.
!>  It sets up initial fields with one or two fronts.
!>  It should also ensure the mesh that was input is valid for the configuration.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_setup_mixed_layer_eddy(domain, iErr)!{{{

   !--------------------------------------------------------------------

      type (domain_type), intent(inout) :: domain
      integer, intent(out) :: iErr
      real (kind=RKIND) :: yMin, yMax, xMin, xMax, dcEdgeMin, dcEdgeMinGlobal
      real (kind=RKIND) :: yMinGlobal, yMaxGlobal, yMidGlobal, xMinGlobal, xMaxGlobal, xMidGlobal
      real (kind=RKIND) :: y1Global, y2Global
      real (kind=RKIND) :: temperature, randomPerturbation, restoringRate

      type (block_type), pointer :: block_ptr

      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: verticalMeshPool
      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: tracersInteriorRestoringFieldsPool

      integer :: iCell, j, k, idx, idx_mixed_layer_base

      ! Define dimension pointers
      integer, pointer :: nCellsSolve, nEdgesSolve, nVertLevels, nVertLevelsP1
      integer, pointer :: index_temperature, index_salinity, index_tracer1

      ! Define variable pointers
      integer, dimension(:), pointer :: maxLevelCell, minLevelCell
      real (kind=RKIND), dimension(:), pointer :: xCell, yCell,refBottomDepth, refZMid, &
                                                  vertCoordMovementWeights, bottomDepth, &
                                                  fCell, fEdge, fVertex, dcEdge
      real (kind=RKIND), dimension(:), pointer :: windStressZonal, windStressMeridional, &
                                                  sensibleHeatFlux, evaporationFlux
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers, debugTracers
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracersInteriorRestoringValue, &
                                                      activeTracersInteriorRestoringRate

      ! Define local interfaceLocations variable
      real (kind=RKIND), dimension(:), pointer :: interfaceLocations

      logical, pointer :: on_a_sphere

      iErr = 0

      if(config_init_configuration .ne. trim('mixed_layer_eddy')) return

     ! Determine vertical grid for configuration
      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
      call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
      call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

      if ( on_a_sphere ) call mpas_log_write('The mixed layer eddy configuration ' &
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

      xMidGlobal = (xMinGlobal + xMaxGlobal) * 0.5_RKIND
      yMidGlobal = (yMinGlobal + yMaxGlobal) * 0.5_RKIND

      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
        call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)
        call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
        call mpas_pool_get_subpool(forcingPool, 'tracersInteriorRestoringFields', tracersInteriorRestoringFieldsPool)

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
        call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
        call mpas_pool_get_array(meshPool, 'fCell', fCell)
        call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
        call mpas_pool_get_array(meshPool, 'fVertex', fVertex)

        call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
        call mpas_pool_get_array(tracersPool, 'debugTracers', debugTracers, 1)
        call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
        call mpas_pool_get_array(forcingPool, 'sensibleHeatFlux', sensibleHeatFlux)
        call mpas_pool_get_array(forcingPool, 'evaporationFlux', evaporationFlux)
        call mpas_pool_get_array(forcingPool, 'windStressZonal', windStressZonal, 1)
        call mpas_pool_get_array(forcingPool, 'windStressMeridional', windStressMeridional, 1)

        call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)
        call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)

        call mpas_pool_get_array(tracersInteriorRestoringFieldsPool, &
             'activeTracersInteriorRestoringRate', activeTracersInteriorRestoringRate, 1)
        call mpas_pool_get_array(tracersInteriorRestoringFieldsPool, &
             'activeTracersInteriorRestoringValue', activeTracersInteriorRestoringValue, 1)

        call ocn_mark_north_boundary(meshPool, yMaxGlobal, dcEdgeMinGlobal, iErr)
        call ocn_mark_south_boundary(meshPool, yMinGlobal, dcEdgeMinGlobal, iErr)

        activeTracersInteriorRestoringRate(:,:,:) = 0.0_RKIND
        activeTracersInteriorRestoringValue(:,:,:) = 0.0_RKIND

        ! Set refBottomDepth and refZMid
        do k = 1, nVertLevels
           refBottomDepth(k) = config_mixed_layer_eddy_bottom_depth * interfaceLocations(k+1)
           refZMid(k) = - 0.5_RKIND * (interfaceLocations(k+1) + interfaceLocations(k)) * config_mixed_layer_eddy_bottom_depth
        end do

        ! Find the index of mixed layer base
        do k = 1, nVertLevels
           if ( refZMid(k) + config_mixed_layer_eddy_mixed_layer_depth <= 0.0_RKIND ) then
              idx_mixed_layer_base = k
              exit
           end if
        end do

        ! Set vertCoordMovementWeights
        vertCoordMovementWeights(:) = 1.0_RKIND

        ! Loop over cells
        do iCell = 1, nCellsSolve
           ! Set surface forcing
           sensibleHeatFlux(iCell) = config_mixed_layer_eddy_heat_flux
           evaporationFlux(iCell) = config_mixed_layer_eddy_evaporation_flux
           windStressZonal(iCell) = config_mixed_layer_eddy_wind_stress_zonal
           windStressMeridional(iCell) = config_mixed_layer_eddy_wind_stress_meridional

           ! Set layerThickness and restingThickness
           do k = 1, nVertLevels
              layerThickness(k, iCell) = config_mixed_layer_eddy_bottom_depth * ( interfaceLocations(k+1) &
                                       - interfaceLocations(k) )
              restingThickness(k, iCell) = config_mixed_layer_eddy_bottom_depth * ( interfaceLocations(k+1) &
                                         - interfaceLocations(k) )
           end do

           ! Set bottomDepth
           bottomDepth(iCell) = config_mixed_layer_eddy_bottom_depth

           ! Set maxLevelCell
           maxLevelCell(iCell) = nVertLevels

           ! Set minLevelCell
           minLevelCell(iCell) = 1

           ! Set debug tracer
           if ( associated(debugTracers) ) then
              idx = index_tracer1
              do k = 1, nVertLevels
                 debugTracers(idx, k, iCell) = 1.0_RKIND
              enddo
           end if
        end do

        ! Set initial condition for temperature and salinity
        if ( associated(activeTracers) ) then
           if ( config_mixed_layer_eddy_two_fronts ) then
              y1Global = ( yMinGlobal + yMidGlobal ) * 0.5_RKIND
              y2Global = ( yMaxGlobal + yMidGlobal ) * 0.5_RKIND
              ! Loop over cells
              do iCell = 1, nCellsSolve
                 ! Set temperature
                 idx = index_temperature
                 do k = nVertLevels, idx_mixed_layer_base, -1
                    temperature = config_mixed_layer_eddy_base_temperature &
                                + config_mixed_layer_eddy_temperature_stratification_interior &
                                * ( refZMid(k) + config_mixed_layer_eddy_mixed_layer_depth ) &
                                + 0.5_RKIND * config_mixed_layer_eddy_temperature_front_width &
                                * config_mixed_layer_eddy_temperature_horizontal_gradient &
                                * ( tanh( 2.0_RKIND * (yCell(iCell) - y1Global) / &
                                        config_mixed_layer_eddy_temperature_front_width ) &
                                -   tanh( 2.0_RKIND * (yCell(iCell) - y2Global) / &
                                        config_mixed_layer_eddy_temperature_front_width ) + 1.0_RKIND )
                    activeTracers(idx, k, iCell) = temperature
                 end do
                 do k = idx_mixed_layer_base-1, 1, -1
                    temperature = config_mixed_layer_eddy_base_temperature &
                                + config_mixed_layer_eddy_temperature_stratification_mixed_layer &
                                * ( refZMid(k) + config_mixed_layer_eddy_mixed_layer_depth ) &
                                + 0.5_RKIND * config_mixed_layer_eddy_temperature_front_width &
                                * config_mixed_layer_eddy_temperature_horizontal_gradient &
                                * ( tanh( 2.0_RKIND * (yCell(iCell) - y1Global) / &
                                        config_mixed_layer_eddy_temperature_front_width ) &
                                -   tanh( 2.0_RKIND * (yCell(iCell) - y2Global) / &
                                        config_mixed_layer_eddy_temperature_front_width ) + 1.0_RKIND )
                    activeTracers(idx, k, iCell) = temperature
                 end do

                 ! Set salinity
                 idx = index_salinity
                 activeTracers(idx, :, iCell) = config_mixed_layer_eddy_salinity
              end do
           else ! single front
              ! Loop over cells
              do iCell = 1, nCellsSolve
                 ! Set temperature
                 idx = index_temperature
                 do k = nVertLevels, idx_mixed_layer_base, -1
                    temperature = config_mixed_layer_eddy_base_temperature &
                                + config_mixed_layer_eddy_temperature_stratification_interior &
                                * ( refZMid(k) + config_mixed_layer_eddy_mixed_layer_depth ) &
                                + 0.5_RKIND * config_mixed_layer_eddy_temperature_front_width &
                                * config_mixed_layer_eddy_temperature_horizontal_gradient &
                                * tanh( 2.0_RKIND * (yCell(iCell) - yMidGlobal) / &
                                        config_mixed_layer_eddy_temperature_front_width )
                    activeTracers(idx, k, iCell) = temperature
                 end do
                 do k = idx_mixed_layer_base-1, 1, -1
                    temperature = config_mixed_layer_eddy_base_temperature &
                                + config_mixed_layer_eddy_temperature_stratification_mixed_layer &
                                * ( refZMid(k) + config_mixed_layer_eddy_mixed_layer_depth ) &
                                + 0.5_RKIND * config_mixed_layer_eddy_temperature_front_width &
                                * config_mixed_layer_eddy_temperature_horizontal_gradient &
                                * tanh( 2.0_RKIND * (yCell(iCell) - yMidGlobal) / &
                                        config_mixed_layer_eddy_temperature_front_width )
                    activeTracers(idx, k, iCell) = temperature
                 end do

                 ! Set salinity
                 idx = index_salinity
                 activeTracers(idx, :, iCell) = config_mixed_layer_eddy_salinity

                 ! Set restoring at boundaries
                 if (yMaxGlobal-yCell(iCell) <= 3.0_RKIND*config_mixed_layer_eddy_restoring_width) then
                    do k = 1, nVertLevels
                       restoringRate = exp( (yCell(iCell) - yMaxGlobal) / &
                                       config_mixed_layer_eddy_restoring_width ) * &
                                       ( 1.0_RKIND / (config_mixed_layer_eddy_restoring_tau*86400.0_RKIND))
                       idx = index_temperature
                       activeTracersInteriorRestoringValue(idx, k, iCell) = activeTracers(idx, k, iCell)
                       activeTracersInteriorRestoringRate(idx, k, iCell) = restoringRate
                       idx = index_salinity
                       activeTracersInteriorRestoringValue(idx, k, iCell) = activeTracers(idx, k, iCell)
                       activeTracersInteriorRestoringRate(idx, k, iCell) = restoringRate
                    end do
                 end if

                 if (yCell(iCell)-yMinGlobal <= 3.0_RKIND*config_mixed_layer_eddy_restoring_width) then
                    do k = 1, nVertLevels
                       restoringRate = exp( (yMinGlobal - yCell(iCell)) / &
                                       config_mixed_layer_eddy_restoring_width ) * &
                                       ( 1.0_RKIND / (config_mixed_layer_eddy_restoring_tau*86400.0_RKIND))
                       idx = index_temperature
                       activeTracersInteriorRestoringValue(idx, k, iCell) = activeTracers(idx, k, iCell)
                       activeTracersInteriorRestoringRate(idx, k, iCell) = restoringRate
                       idx = index_salinity
                       activeTracersInteriorRestoringValue(idx, k, iCell) = activeTracers(idx, k, iCell)
                       activeTracersInteriorRestoringRate(idx, k, iCell) = restoringRate
                    end do
                 end if
              end do
           end if

           ! Add random perturbation to initial temperature
           idx = index_temperature

           ! Initialize random number generator
           call random_seed()

           ! Loop over cells
           do iCell = 1, nCellsSolve
              ! Random perturbation [0,1)
              call random_number(randomPerturbation)
              ! [0,1) to [-1,1)
              randomPerturbation = 2.0_RKIND * randomPerturbation - 1.0_RKIND
              ! Add random perturbation, constant in z
              activeTracers(idx, :, iCell) = activeTracers(idx, :, iCell) &
                             + config_mixed_layer_eddy_temperature_perturbation_magnitude &
                             * randomPerturbation
           end do
        end if

        ! Set Coriolis parameters
        fCell(:) = config_mixed_layer_eddy_coriolis_parameter
        fEdge(:) = config_mixed_layer_eddy_coriolis_parameter
        fVertex(:) = config_mixed_layer_eddy_coriolis_parameter

        block_ptr => block_ptr % next
      end do

      deallocate(interfaceLocations)

   !--------------------------------------------------------------------

   end subroutine ocn_init_setup_mixed_layer_eddy!}}}

!***********************************************************************
!
!  routine ocn_init_validate_mixed_layer_eddy
!
!> \brief   Validation for mixed layer eddy test case
!> \author  Qing Li
!> \date    01/17/2020
!> \details
!>  This routine validates the configuration options for the mixed layer eddy test case.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_validate_mixed_layer_eddy(configPool, packagePool, iocontext, iErr)!{{{

   !--------------------------------------------------------------------
      type (mpas_pool_type), intent(inout) :: configPool, packagePool
      type (mpas_io_context_type), intent(inout) :: iocontext

      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration
      integer, pointer :: config_vert_levels, config_mixed_layer_eddy_vert_levels

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)

      if(config_init_configuration .ne. trim('mixed_layer_eddy')) return

      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      call mpas_pool_get_config(configPool, 'config_mixed_layer_eddy_vert_levels', config_mixed_layer_eddy_vert_levels)

      if(config_vert_levels <= 0 .and. config_mixed_layer_eddy_vert_levels > 0) then
         config_vert_levels = config_mixed_layer_eddy_vert_levels
      else if (config_vert_levels <= 0) then
         call mpas_log_write( 'Validation failed for mixed layer eddy. Not given a usable value for vertical levels.', MPAS_LOG_CRIT)
         iErr = 1
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_mixed_layer_eddy!}}}

!***********************************************************************

end module ocn_init_mixed_layer_eddy

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
