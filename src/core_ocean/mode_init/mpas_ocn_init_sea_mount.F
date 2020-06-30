! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_sea_mount
!
!> \brief MPAS ocean initialize case -- Sea Mount
!> \author Mark Petersen
!> \date   08/10/15
!> \details
!>  This module contains the routines for initializing the
!>  the sea mount test case
!
!-----------------------------------------------------------------------

module ocn_init_sea_mount

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

   public :: ocn_init_setup_sea_mount, &
             ocn_init_validate_sea_mount

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_sea_mount
!
!> \brief   Setup for sea mount test case
!> \author  Mark Petersen
!> \date    08/10/15
!> \details
!>  This routine sets up the initial conditions for the sea mount test case.
!>  It should also ensure the mesh that was input is valid for the configuration.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_setup_sea_mount(domain, iErr)!{{{

   !--------------------------------------------------------------------

      type (domain_type), intent(inout) :: domain
      integer, intent(out) :: iErr
      real (kind=RKIND) :: yMin, yMax, xMin, xMax, dcEdgeMin, dcEdgeMinGlobal
      real (kind=RKIND) :: yMinGlobal, yMaxGlobal, yMidGlobal, xMinGlobal, xMaxGlobal, xMidGlobal
      real (kind=RKIND) :: densityCell, z, radius

      type (block_type), pointer :: block_ptr

      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: verticalMeshPool
      type (mpas_pool_type), pointer :: tracersPool

      integer :: iCell, k, idx

      ! Define config variable pointers
      character (len=StrKIND), pointer :: config_init_configuration, config_vertical_grid, config_sea_mount_layer_type, &
                                          config_sea_mount_stratification_type
      real (kind=RKIND), pointer :: config_sea_mount_width, &
                                    config_sea_mount_bottom_depth, config_sea_mount_height,config_sea_mount_radius, &
                                    config_sea_mount_density_coef_linear, config_sea_mount_density_coef_exp, &
                                    config_sea_mount_density_gradient_linear, config_sea_mount_density_gradient_exp, &
                                    config_sea_mount_density_depth_linear, config_sea_mount_density_depth_exp, &
                                    config_sea_mount_density_ref, config_sea_mount_density_alpha, config_sea_mount_density_Tref, &
                                    config_sea_mount_salinity, config_sea_mount_coriolis_parameter

      ! Define dimension pointers
      integer, pointer :: nCellsSolve, nEdgesSolve, nVertLevels, nVertLevelsP1
      integer, pointer :: index_temperature, index_salinity

      ! Define variable pointers
      integer, dimension(:), pointer :: maxLevelCell
      real (kind=RKIND), dimension(:), pointer :: xCell, yCell,refBottomDepth, refZMid, refLayerThickness, &
                                                  vertCoordMovementWeights, bottomDepth, &
                                                  fCell, fEdge, fVertex, dcEdge
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers

      ! Define local interfaceLocations variable
      real (kind=RKIND), dimension(:), pointer :: interfaceLocations

      logical, pointer :: on_a_sphere

      iErr = 0

      call mpas_pool_get_config(ocnConfigs, 'config_init_configuration', config_init_configuration)

      if(config_init_configuration .ne. trim('sea_mount')) return

      call mpas_pool_get_config(ocnConfigs, 'config_vertical_grid', config_vertical_grid)

      call mpas_pool_get_config(ocnConfigs, 'config_sea_mount_layer_type', config_sea_mount_layer_type)
      call mpas_pool_get_config(ocnConfigs, 'config_sea_mount_stratification_type', config_sea_mount_stratification_type)
      call mpas_pool_get_config(ocnConfigs, 'config_sea_mount_width', config_sea_mount_width)
      call mpas_pool_get_config(ocnConfigs, 'config_sea_mount_bottom_depth', config_sea_mount_bottom_depth)
      call mpas_pool_get_config(ocnConfigs, 'config_sea_mount_height', config_sea_mount_height)
      call mpas_pool_get_config(ocnConfigs, 'config_sea_mount_radius', config_sea_mount_radius)
      call mpas_pool_get_config(ocnConfigs, 'config_sea_mount_density_coef_linear', config_sea_mount_density_coef_linear)
      call mpas_pool_get_config(ocnConfigs, 'config_sea_mount_density_gradient_linear', config_sea_mount_density_gradient_linear)
      call mpas_pool_get_config(ocnConfigs, 'config_sea_mount_density_depth_linear', config_sea_mount_density_depth_linear)
      call mpas_pool_get_config(ocnConfigs, 'config_sea_mount_density_coef_exp', config_sea_mount_density_coef_exp)
      call mpas_pool_get_config(ocnConfigs, 'config_sea_mount_density_gradient_exp', config_sea_mount_density_gradient_exp)
      call mpas_pool_get_config(ocnConfigs, 'config_sea_mount_density_depth_exp', config_sea_mount_density_depth_exp)
      call mpas_pool_get_config(ocnConfigs, 'config_sea_mount_density_ref', config_sea_mount_density_ref)
      call mpas_pool_get_config(ocnConfigs, 'config_sea_mount_density_alpha', config_sea_mount_density_alpha)
      call mpas_pool_get_config(ocnConfigs, 'config_sea_mount_density_Tref', config_sea_mount_density_Tref)
      call mpas_pool_get_config(ocnConfigs, 'config_sea_mount_salinity', config_sea_mount_salinity)
      call mpas_pool_get_config(ocnConfigs, 'config_sea_mount_coriolis_parameter', config_sea_mount_coriolis_parameter)

      ! Determine vertical grid for configuration
      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
      call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
      call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

      if ( on_a_sphere ) call mpas_log_write('The sea mount configuration can ' &
           // 'only be applied to a planar mesh. Exiting...', MPAS_LOG_CRIT)

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
        call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)

        call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
        call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

        call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
        call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
        call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)
        call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)

        call mpas_pool_get_array(meshPool, 'xCell', xCell)
        call mpas_pool_get_array(meshPool, 'yCell', yCell)
        call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
        call mpas_pool_get_array(meshPool, 'vertCoordMovementWeights', vertCoordMovementWeights)
        call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
        call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
        call mpas_pool_get_array(meshPool, 'fCell', fCell)
        call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
        call mpas_pool_get_array(meshPool, 'fVertex', fVertex)

        call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)

        call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)
        call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)
        call mpas_pool_get_array(verticalMeshPool, 'refLayerThickness', refLayerThickness)

        ! mrp note: doubly non-periodic.  Can delete these later.
        call ocn_mark_north_boundary(meshPool, yMaxGlobal, dcEdgeMinGlobal, iErr)
        call ocn_mark_south_boundary(meshPool, yMinGlobal, dcEdgeMinGlobal, iErr)
        call ocn_mark_east_boundary(meshPool, xMaxGlobal, dcEdgeMinGlobal, iErr)
        call ocn_mark_west_boundary(meshPool, xMinGlobal, dcEdgeMinGlobal, iErr)

        ! Set refBottomDepth and refZMid
        do k = 1, nVertLevels
            refBottomDepth(k) = config_sea_mount_bottom_depth * interfaceLocations(k+1)
            refZMid(k) = - 0.5_RKIND * (interfaceLocations(k+1) + interfaceLocations(k)) * config_sea_mount_bottom_depth
        end do

        refLayerThickness(1) = refBottomDepth(1)
        do k = 2, nVertLevels
            refLayerThickness(k) = refBottomDepth(k) - refBottomDepth(k-1)
        end do

        ! Set vertCoordMovementWeights
        vertCoordMovementWeights(:) = 1.0_RKIND

        do iCell = 1, nCellsSolve

           radius =  sqrt( (xCell(iCell)-xMidGlobal)**2 + (yCell(iCell)-yMidGlobal)**2 )

           ! Set bottomDepth.  See Beckmann and Haidvogel 1993 eqn 12, Shchepetkin 2003 eqn 4.2
           bottomDepth(iCell) = config_sea_mount_bottom_depth - config_sea_mount_height &
                              * exp(-(max(radius-config_sea_mount_radius, 0.0_RKIND))**2 / config_sea_mount_width**2)

           ! Set maxLevelCell and layerThickness
           if ( trim(config_sea_mount_layer_type) == 'z-level' ) then
              maxLevelCell(iCell) = -1
              do k = 1, nVertLevels
                 if (bottomDepth(iCell) .le. refBottomDepth(k)) then
                    maxLevelCell(iCell) = k
                    ! make full cell only:
                    bottomDepth(iCell) = refBottomDepth(k)
                    exit
                 end if
              end do
              do k = 1, maxLevelCell(iCell)
                 layerThickness(k, iCell) = refLayerThickness(k)
              end do
           else if ( trim(config_sea_mount_layer_type) == 'sigma') then
              maxLevelCell(iCell) = nVertLevels
              do k = 1, nVertLevels
                 layerThickness(k, iCell) = bottomDepth(iCell) / nVertLevels
              end do
           end if

           ! Set restingThickness
           do k = 1, maxLevelCell(iCell)
              restingThickness(k, iCell) = layerThickness(k, iCell)
           end do

           ! Set stratification using temperature.  See Beckmann and Haidvogel 1993 eqn 15-16.
           idx = index_temperature
           z = 0.0_RKIND
           do k = 1, maxLevelCell(iCell)

              z = z - 0.5_RKIND * layerThickness(k, iCell)

              if ( trim(config_sea_mount_stratification_type) == 'linear' ) then
                 densityCell =  config_sea_mount_density_coef_linear - config_sea_mount_density_gradient_linear * z &
                             / config_sea_mount_density_depth_linear
              elseif ( trim(config_sea_mount_stratification_type) == 'exponential' ) then
                 densityCell = config_sea_mount_density_coef_exp - config_sea_mount_density_gradient_exp * exp( z &
                             / config_sea_mount_density_depth_exp)
              else
                 call mpas_log_write('MPAS-ocean: Error: Incorrect config_sea_mount_stratification_type: ' &
                                              // config_sea_mount_stratification_type)
              end if

              ! Back-solve linear EOS for temperature, with S=S_ref
              ! T = T_ref - (rho - rho_ref)/alpha
              activeTracers(idx, k, iCell) = config_sea_mount_density_Tref - (densityCell - config_sea_mount_density_ref) &
                                           / config_sea_mount_density_alpha
              z = z - 0.5_RKIND * layerThickness(k, iCell)

           end do

           ! Set salinity
           idx = index_salinity
           activeTracers(idx, :, iCell) = config_sea_mount_salinity

        end do

        ! Set Coriolis parameters
        fCell(:) = config_sea_mount_coriolis_parameter
        fEdge(:) = config_sea_mount_coriolis_parameter
        fVertex(:) = config_sea_mount_coriolis_parameter

        block_ptr => block_ptr % next
      end do

      deallocate(interfaceLocations)

   !--------------------------------------------------------------------

   end subroutine ocn_init_setup_sea_mount!}}}

!***********************************************************************
!
!  routine ocn_init_validate_sea_mount
!
!> \brief   Validation for sea mount test case
!> \author  Mark Petersen
!> \date    08/10/15
!> \details
!>  This routine validates the configuration options for the sea mount test case.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_validate_sea_mount(configPool, packagePool, iocontext, iErr)!{{{

   !--------------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: configPool, packagePool
      type (mpas_io_context_type), intent(inout), target :: iocontext

      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration
      integer, pointer :: config_vert_levels, config_sea_mount_vert_levels

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)

      if(config_init_configuration .ne. trim('sea_mount')) return

      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      call mpas_pool_get_config(configPool, 'config_sea_mount_vert_levels', config_sea_mount_vert_levels)

      if(config_vert_levels <= 0 .and. config_sea_mount_vert_levels > 0) then
         config_vert_levels = config_sea_mount_vert_levels
      else if (config_vert_levels <= 0) then
         call mpas_log_write( 'Validation failed for sea mount. Not given a usable value for vertical levels.', MPAS_LOG_CRIT)
         iErr = 1
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_sea_mount!}}}

!***********************************************************************

end module ocn_init_sea_mount

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
