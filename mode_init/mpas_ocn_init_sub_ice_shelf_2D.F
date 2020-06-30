! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_sub_ice_shelf_2D
!
!> \brief MPAS ocean initialize case -- sub_ice_shelf_2D
!> \author Mark Petersen
!> \date   9/2/2015
!
!-----------------------------------------------------------------------

module ocn_init_sub_ice_shelf_2D

   use mpas_kind_types
   use mpas_io_units
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar

   use ocn_constants
   use ocn_init_vertical_grids
   use ocn_init_cell_markers

   use ocn_init_ssh_and_landIcePressure

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

   public :: ocn_init_setup_sub_ice_shelf_2D, &
             ocn_init_validate_sub_ice_shelf_2D

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

 contains

   !***********************************************************************
   !
   !  routine ocn_init_setup_sub_ice_shelf_2D
   !
   !> \brief   Setup for this initial condition
   !> \author  Mark Petersen
   !> \date    9/2/2015
   !> \details
   !>  This routine sets up the initial conditions for this case.
   !
   !-----------------------------------------------------------------------

   subroutine ocn_init_setup_sub_ice_shelf_2D(domain, iErr)!{{{
     !--------------------------------------------------------------------

     type (domain_type), intent(inout) :: domain
     integer, intent(out) :: iErr
     real (kind=RKIND) :: yMin, yMax, xMin, xMax, dcEdgeMin, dcEdgeMinGlobal, maxDepth
     real (kind=RKIND) :: yMinGlobal, yMaxGlobal, yMidGlobal, xMinGlobal, xMaxGlobal
     real (kind=RKIND) :: totalSubIceThickness, y1,y2,y3, d1,d2,d3, surfaceDepression, surfaceSalinity, bottomSalinity

     type (block_type), pointer :: block_ptr

     type (mpas_pool_type), pointer :: meshPool
     type (mpas_pool_type), pointer :: statePool
     type (mpas_pool_type), pointer :: tracersPool
     type (mpas_pool_type), pointer :: forcingPool
     type (mpas_pool_type), pointer :: verticalMeshPool

     integer :: iCell, k, idx

     ! Define config variable pointers
     character (len=StrKIND), pointer :: config_init_configuration, config_vertical_grid
     real (kind=RKIND), pointer :: config_sub_ice_shelf_2D_bottom_depth, &
          config_sub_ice_shelf_2D_cavity_thickness, config_sub_ice_shelf_2D_edge_width, config_sub_ice_shelf_2D_temperature, &
          config_sub_ice_shelf_2D_surface_salinity, config_sub_ice_shelf_2D_bottom_salinity, &
          config_sub_ice_shelf_2D_y1, config_sub_ice_shelf_2D_y2,config_sub_ice_shelf_2D_slope_height

     ! Define dimension pointers
     integer, pointer :: nCellsSolve, nEdgesSolve, nVertLevels, nVertLevelsP1, nCells
     integer, pointer :: index_temperature, index_salinity

     ! Define variable pointers
     integer, dimension(:), pointer :: maxLevelCell, modifySSHMask, landIceMask
     real (kind=RKIND), dimension(:), pointer :: xCell, yCell,refBottomDepth, refZMid, &
          vertCoordMovementWeights, bottomDepth, &
          fCell, fEdge, fVertex, dcEdge, refLayerThickness
     real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers

     ! Define local interfaceLocations variable
     real (kind=RKIND), dimension(:), pointer :: interfaceLocations

     logical, pointer :: on_a_sphere

     type (mpas_pool_type), pointer :: diagnosticsPool
     real(kind=RKIND), dimension(:), pointer :: landIceFraction, ssh
     real (kind=RKIND), dimension(:,:), pointer :: zMid

     iErr = 0

     call mpas_pool_get_config(ocnConfigs, 'config_init_configuration', config_init_configuration)
     if(config_init_configuration .ne. trim('sub_ice_shelf_2D')) return

     call mpas_pool_get_config(ocnConfigs, 'config_vertical_grid', config_vertical_grid)
     call mpas_pool_get_config(ocnConfigs, 'config_sub_ice_shelf_2D_bottom_depth',config_sub_ice_shelf_2D_bottom_depth)
     call mpas_pool_get_config(ocnConfigs, 'config_sub_ice_shelf_2D_cavity_thickness', config_sub_ice_shelf_2D_cavity_thickness)
     call mpas_pool_get_config(ocnConfigs, 'config_sub_ice_shelf_2D_edge_width', config_sub_ice_shelf_2D_edge_width)
     call mpas_pool_get_config(ocnConfigs, 'config_sub_ice_shelf_2D_temperature', config_sub_ice_shelf_2D_temperature)
     call mpas_pool_get_config(ocnConfigs, 'config_sub_ice_shelf_2D_surface_salinity', config_sub_ice_shelf_2D_surface_salinity)
     call mpas_pool_get_config(ocnConfigs, 'config_sub_ice_shelf_2D_bottom_salinity', config_sub_ice_shelf_2D_bottom_salinity)
     call mpas_pool_get_config(ocnConfigs, 'config_sub_ice_shelf_2D_y1', config_sub_ice_shelf_2D_y1)
     call mpas_pool_get_config(ocnConfigs, 'config_sub_ice_shelf_2D_y2', config_sub_ice_shelf_2D_y2)
     call mpas_pool_get_config(ocnConfigs, 'config_sub_ice_shelf_2D_slope_height', config_sub_ice_shelf_2D_slope_height)

     ! points 1 and 2 are where angles on ice shelf are located.
     ! point 3 is at the surface.
     ! d variables are total water thickness below ice shelf.
     y1=config_sub_ice_shelf_2D_y1
     y2=config_sub_ice_shelf_2D_y2
     y3=config_sub_ice_shelf_2D_y2 + config_sub_ice_shelf_2D_edge_width
     d1=config_sub_ice_shelf_2D_cavity_thickness
     d2=config_sub_ice_shelf_2D_cavity_thickness+config_sub_ice_shelf_2D_slope_height
     d3=config_sub_ice_shelf_2D_bottom_depth

     bottomSalinity = config_sub_ice_shelf_2D_bottom_salinity
     surfaceSalinity = config_sub_ice_shelf_2D_surface_salinity

     ! Determine vertical grid for configuration
     call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
     call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
     call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
     call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

     if ( on_a_sphere ) call mpas_log_write('The sub_ice_shelf_2D configuration can ' &
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

     block_ptr => domain % blocklist
     do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
        call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)

        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

        call mpas_pool_get_array(meshPool, 'xCell', xCell)
        call mpas_pool_get_array(meshPool, 'yCell', yCell)
        call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
        call mpas_pool_get_array(meshPool, 'vertCoordMovementWeights', vertCoordMovementWeights)
        call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
        call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

        call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)
        call mpas_pool_get_array(verticalMeshPool, 'refLayerThickness', refLayerThickness)

        call mpas_pool_get_array(forcingPool, 'landIceFraction', landIceFraction)
        call mpas_pool_get_array(forcingPool, 'landIceMask', landIceMask)

        call mpas_pool_get_array(statePool, 'ssh', ssh, 1)
        call mpas_pool_get_array(diagnosticsPool, 'modifySSHMask', modifySSHMask)
        call ocn_mark_north_boundary(meshPool, yMaxGlobal, dcEdgeMinGlobal, iErr)
        call ocn_mark_south_boundary(meshPool, yMinGlobal, dcEdgeMinGlobal, iErr)

        ! Set refBottomDepth and refZMid
        do k = 1, nVertLevels
           refBottomDepth(k) = config_sub_ice_shelf_2D_bottom_depth * interfaceLocations(k+1)
        end do

        ! Compute refLayerThickness and refZMid
        call ocn_compute_layerThickness_zMid_from_bottomDepth(refLayerThickness,refZMid, &
             refBottomDepth,refBottomDepth(nVertLevels), &
             nVertLevels,nVertLevels,iErr)

        ! Set vertCoordMovementWeights
        vertCoordMovementWeights(:) = 1.0_RKIND

        maxDepth = refBottomDepth(nVertLevels)

        if(associated(landIceFraction)) &
          landIceFraction(:) = 0.0_RKIND
        modifySSHMask(:) = 0
        if(associated(landIceMask)) &
          landIceMask(:) = 0

        do iCell = 1, nCells
           ! set up sub ice shelf thicknesses
           if (yCell(iCell) < y1 ) then
              totalSubIceThickness = d1
           elseif (yCell(iCell) < y2 ) then
              totalSubIceThickness = d1 + (d2-d1)*(yCell(iCell)-y1)/(y2-y1)
           elseif (yCell(iCell) < y3 ) then
              totalSubIceThickness = d2 + (d3-d2)*(yCell(iCell)-y2)/(y3-y2)
           else
              totalSubIceThickness = d3
           endif
           ssh(iCell) = -config_sub_ice_shelf_2D_bottom_depth + totalSubIceThickness

           if(ssh(iCell) < 0.0_RKIND) then
              modifySSHMask(iCell) = 1
           end if
        end do

        do iCell = 1, nCells
           if (yCell(iCell) < y3 ) then
              if(associated(landIceFraction)) &
                landIceFraction(iCell) = 1.0_RKIND
              if(associated(landIceMask)) &
                landIceMask(iCell) = 1
           end if

           ! Set bottomDepth
           bottomDepth(iCell) = refBottomDepth(nVertLevels)

           ! Set maxLevelCell
           maxLevelCell(iCell) = nVertLevels
        end do

        block_ptr => block_ptr % next
     end do

     ! compute the vertical grid (layerThickness, restingThickness, maxLevelCell, zMid) based on ssh,
     ! bottomDepth and refBottomDepth
     call ocn_init_ssh_and_landIcePressure_vertical_grid(domain, iErr)

     if(iErr .ne. 0) then
        call mpas_log_write( 'ocn_init_ssh_and_landIcePressure_vertical_grid failed.', MPAS_LOG_CRIT)
        call mpas_dmpar_finalize(domain % dminfo)
     end if

     block_ptr => domain % blocklist
     do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
        call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
        call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

        call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
        call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)

        call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)

        call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)

        ! compute active tracer fields
        ! If we are constructing an initial guess (rather than reading ssh in from a stream), these are reference activeTracers
        ! on a vertical grid that has not been displaced by the ssh
        do iCell = 1, nCells
           ! Set temperature
           idx = index_temperature
           do k = 1, nVertLevels
              activeTracers(idx, k, iCell) = config_sub_ice_shelf_2D_temperature
           end do

           ! Set up salinity stratification
           idx = index_salinity
           do k = 1, nVertLevels
              activeTracers(idx, k, iCell) = surfaceSalinity + (bottomSalinity - surfaceSalinity) &
                                                             * (zMid(k,iCell)/(-config_sub_ice_shelf_2D_bottom_depth))
           end do
        end do

        block_ptr => block_ptr % next
     end do

     ! compute or update the land-ice pressure (or possibly SSH), also computing density along the way
     ! If this is the initial guess, the vertical grid and activeTracers may also be recomputed based on SSH
     call ocn_init_ssh_and_landIcePressure_balance(domain, iErr)

     if(iErr .ne. 0) then
        call mpas_log_write( 'ocn_init_ssh_and_landIcePressure_balance failed.', MPAS_LOG_CRIT)
        call mpas_dmpar_finalize(domain % dminfo)
     end if

     call ocn_compute_Haney_number(domain, iErr)

     if(iErr .ne. 0) then
        call mpas_log_write( 'ocn_compute_Haney_number failed.', MPAS_LOG_CRIT)
        call mpas_dmpar_finalize(domain % dminfo)
     end if

     !--------------------------------------------------------------------

   end subroutine ocn_init_setup_sub_ice_shelf_2D!}}}

   !***********************************************************************
   !
   !  routine ocn_init_validate_sub_ice_shelf_2D
   !
   !> \brief   Validation for this initial condition
   !> \author  Mark Petersen
   !> \date    9/2/2015
   !> \details
   !>  This routine validates the configuration options for this case.
   !
   !-----------------------------------------------------------------------

   subroutine ocn_init_validate_sub_ice_shelf_2D(configPool, packagePool, iocontext, iErr)!{{{

     !--------------------------------------------------------------------
     type (mpas_pool_type), intent(in) :: configPool, packagePool
      type (mpas_io_context_type), intent(inout), target :: iocontext

     integer, intent(out) :: iErr

     character (len=StrKIND), pointer :: config_init_configuration
     integer, pointer :: config_vert_levels, config_sub_ice_shelf_2D_vert_levels

     type (mpas_io_context_type), pointer :: iocontext_ptr

     iocontext_ptr => iocontext

     iErr = 0

     call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)

     if(config_init_configuration .ne. trim('sub_ice_shelf_2D')) return

     call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
     call mpas_pool_get_config(configPool, 'config_sub_ice_shelf_2D_vert_levels', config_sub_ice_shelf_2D_vert_levels)

     if(config_vert_levels <= 0 .and. config_sub_ice_shelf_2D_vert_levels > 0) then
        config_vert_levels = config_sub_ice_shelf_2D_vert_levels
     else if (config_vert_levels <= 0) then
        call mpas_log_write( 'Validation failed for sub_ice_shelf_2D. Not given a usable value for vertical levels.', MPAS_LOG_CRIT)
        iErr = 1
     end if

     !--------------------------------------------------------------------

   end subroutine ocn_init_validate_sub_ice_shelf_2D!}}}


   !***********************************************************************

 end module ocn_init_sub_ice_shelf_2D

 !|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 ! vim: foldmethod=marker
