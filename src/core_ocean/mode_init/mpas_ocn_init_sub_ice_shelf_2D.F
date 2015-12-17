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
     real (kind=RKIND) :: totalSubIceThickness, y1,y2,y3, d1,d2,d3, surfaceDepression
     real (kind=RKIND), dimension(:), allocatable :: midDepth

     type (block_type), pointer :: block_ptr

     type (mpas_pool_type), pointer :: meshPool
     type (mpas_pool_type), pointer :: statePool
     type (mpas_pool_type), pointer :: tracersPool
     type (mpas_pool_type), pointer :: forcingPool
     type (mpas_pool_type), pointer :: verticalMeshPool

     integer :: iCell, k, idx

     ! Define config variable pointers
     character (len=StrKIND), pointer :: config_init_configuration, config_vertical_grid
     real (kind=RKIND), pointer :: config_sub_ice_shelf_2D_bottom_depth, config_sub_ice_shelf_2D_seaSurfacePressure_factor, &
          config_sub_ice_shelf_2D_cavity_thickness, config_sub_ice_shelf_2D_edge_width, config_sub_ice_shelf_2D_temperature, &
          config_sub_ice_shelf_2D_surface_salinity, config_sub_ice_shelf_2D_bottom_salinity, &
          config_sub_ice_shelf_2D_y1, config_sub_ice_shelf_2D_y2,config_sub_ice_shelf_2D_slope_height

     real (kind=RKIND), pointer :: config_eos_linear_alpha, config_eos_linear_beta, config_eos_linear_Tref, config_eos_linear_Sref, config_eos_linear_densityref

     ! Define dimension pointers
     integer, pointer :: nCellsSolve, nEdgesSolve, nVertLevels, nVertLevelsP1
     integer, pointer :: index_temperature, index_salinity

     ! Define variable pointers
     integer, dimension(:), pointer :: maxLevelCell
     real (kind=RKIND), dimension(:), pointer :: xCell, yCell,refBottomDepth, refZMid, &
          vertCoordMovementWeights, bottomDepth, seaSurfacePressure, &
          fCell, fEdge, fVertex, dcEdge
     real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness
     real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers

     ! Define local interfaceLocations variable
     real (kind=RKIND), dimension(:), pointer :: interfaceLocations

     logical, pointer :: on_a_sphere

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
     call mpas_pool_get_config(ocnConfigs, 'config_sub_ice_shelf_2D_seaSurfacePressure_factor',config_sub_ice_shelf_2D_seaSurfacePressure_factor)
     call mpas_pool_get_config(ocnConfigs, 'config_sub_ice_shelf_2D_y1', config_sub_ice_shelf_2D_y1)
     call mpas_pool_get_config(ocnConfigs, 'config_sub_ice_shelf_2D_y2', config_sub_ice_shelf_2D_y2)
     call mpas_pool_get_config(ocnConfigs, 'config_sub_ice_shelf_2D_slope_height', config_sub_ice_shelf_2D_slope_height)

     call mpas_pool_get_config(ocnConfigs, 'config_eos_linear_alpha', config_eos_linear_alpha)
     call mpas_pool_get_config(ocnConfigs, 'config_eos_linear_beta', config_eos_linear_beta)
     call mpas_pool_get_config(ocnConfigs, 'config_eos_linear_Tref', config_eos_linear_Tref)
     call mpas_pool_get_config(ocnConfigs, 'config_eos_linear_Sref', config_eos_linear_Sref)
     call mpas_pool_get_config(ocnConfigs, 'config_eos_linear_densityref', config_eos_linear_densityref)

     ! points 1 and 2 are where angles on ice shelf are located.
     ! point 3 is at the surface.
     ! d variables are total water thickness below ice shelf.
     y1=config_sub_ice_shelf_2D_y1
     y2=config_sub_ice_shelf_2D_y2
     y3=config_sub_ice_shelf_2D_y2 + config_sub_ice_shelf_2D_edge_width
     d1=config_sub_ice_shelf_2D_cavity_thickness
     d2=config_sub_ice_shelf_2D_cavity_thickness+config_sub_ice_shelf_2D_slope_height
     d3=config_sub_ice_shelf_2D_bottom_depth

     config_sub_ice_shelf_2D_bottom_salinity = config_sub_ice_shelf_2D_bottom_salinity
     config_sub_ice_shelf_2D_surface_salinity = config_sub_ice_shelf_2D_surface_salinity

     ! Determine vertical grid for configuration
     call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
     call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
     call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
     call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

     allocate(midDepth(nVertLevels))

     ! you may restrict your case geometry as follows:
     if ( on_a_sphere ) call mpas_dmpar_global_abort('ERROR: The sub_ice_shelf_2D configuration can only be applied to a planar mesh. Exiting...')

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
        call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

        call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

        call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
        call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)

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
        call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)

        call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)
        call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)

        call mpas_pool_get_array(forcingPool, 'seaSurfacePressure',seaSurfacePressure)

        call ocn_mark_north_boundary(meshPool, yMaxGlobal, dcEdgeMinGlobal, iErr)
        call ocn_mark_south_boundary(meshPool, yMinGlobal, dcEdgeMinGlobal, iErr)

        ! Set refBottomDepth and refZMid
        do k = 1, nVertLevels
           refBottomDepth(k) = config_sub_ice_shelf_2D_bottom_depth * interfaceLocations(k+1)
           refZMid(k) = - 0.5_RKIND * (interfaceLocations(k+1) + interfaceLocations(k)) * config_sub_ice_shelf_2D_bottom_depth
        end do

        ! Set vertCoordMovementWeights
        vertCoordMovementWeights(:) = 1.0_RKIND

        maxDepth = refBottomDepth(nVertLevels)

        do iCell = 1, nCellsSolve

           ! Set layerThickness and restingThickness
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

           do k = 1, nVertLevels
              layerThickness(k, iCell) = totalSubIceThickness/nVertLevels
              restingThickness(k, iCell) = config_sub_ice_shelf_2D_bottom_depth/nVertLevels
           end do

           ! Set bottomDepth
           bottomDepth(iCell) = refBottomDepth(nVertLevels)

           ! Set maxLevelCell
           maxLevelCell(iCell) = nVertLevels

           ! Set temperature
           idx = index_temperature
           do k = 1, nVertLevels
              activeTracers(idx, k, iCell) = config_sub_ice_shelf_2D_temperature
           end do

           ! Set up salinity stratification
           midDepth(nVertLevels) =  config_sub_ice_shelf_2D_bottom_depth - 0.5*layerThickness(nVertLevels,iCell)
           do k=nVertLevels-1,1,-1
              midDepth(k) = midDepth(k+1) - 0.5*(layerThickness(k+1,iCell)+layerThickness(k,iCell))
           enddo
           idx = index_salinity
           do k = nVertLevels, 1, -1
              activeTracers(idx, k, iCell) = config_sub_ice_shelf_2D_surface_salinity + (config_sub_ice_shelf_2D_bottom_salinity - config_sub_ice_shelf_2D_surface_salinity) * (midDepth(k)/config_sub_ice_shelf_2D_bottom_depth)
           end do

           ! Set sea surface pressure
           surfaceDepression = bottomDepth(iCell) - sum(layerThickness(1:maxLevelCell(iCell),iCell))

           ! Compute sea surface pressure as  
           !    p = g int_{z=-eta}}{0} rho(z) dz
           ! that is, integrate density from the bottom of the sea surface (bottom of ice) to z=0
           ! with some assumed density profile.  Here we use a linear equation of state for the density
           ! and integrate analytically:
           !   p = g c int_{z=-eta}}{0} ( rho_0 - alpha*(T-Tref) + beta*(S-Sref) ) dz
           ! Note this formulation is for a constant T and linear stratification is S.
           seaSurfacePressure(iCell) = gravity * surfaceDepression * config_sub_ice_shelf_2D_seaSurfacePressure_factor*(&
                config_eos_linear_densityref &
                - config_eos_linear_alpha * ( config_sub_ice_shelf_2D_temperature - config_eos_linear_Tref) &
                + config_eos_linear_beta  * ( -config_eos_linear_Sref + config_sub_ice_shelf_2D_surface_salinity  &
                + (config_sub_ice_shelf_2D_bottom_salinity-config_sub_ice_shelf_2D_surface_salinity)/maxDepth/2.0* surfaceDepression ))

        end do

        block_ptr => block_ptr % next
     end do

     deallocate(midDepth)
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
        write(stderrUnit,*) 'IERROR: Validation failed for sub_ice_shelf_2D. Not given a usable value for vertical levels.'
        iErr = 1
     end if

     !--------------------------------------------------------------------

   end subroutine ocn_init_validate_sub_ice_shelf_2D!}}}


   !***********************************************************************

 end module ocn_init_sub_ice_shelf_2D

 !|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 ! vim: foldmethod=marker
