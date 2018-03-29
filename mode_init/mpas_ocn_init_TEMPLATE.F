! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_TEMPLATE
!
!> \brief MPAS ocean initialize case -- TEMPLATE
!> \author FILL_IN_AUTHOR
!> \date   FILL_IN_DATE
!> \details
!>  This module contains the routines for initializing the
!>  TEMPLATE initial condition
!>
!>  In order to add a new analysis member, do the following:
!>  1. In src/core_ocean/mode_init, copy these to your new analysis member name:
!>     cp mpas_ocn_init_TEMPLATE.F mpas_ocn_init_your_new_name.F
!>     cp Registry_TEMPLATE.xml Registry_ocn_your_new_name.xml
!>
!>  2. In those two new files, replace the following text:
!>     TEMPLATE, FILL_IN_AUTHOR, FILL_IN_DATE
!>     TEMPLATE uses underscores (subroutine names), like your_new_name.
!>
!>  3. Add a #include line for your registry to
!>     src/core_ocean/mode_init/Registry.xml
!>
!>  4. Copy and change TEMPLATE lines in src/core_ocean/mode_init/mpas_ocn_init_mode.F
!>
!>  5. Add these dependency lines by following TEMPLATE examples in:
!>     in src/core_ocean/mode_init/Makefile
!
!-----------------------------------------------------------------------

module ocn_init_TEMPLATE

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

   public :: ocn_init_setup_TEMPLATE, &
             ocn_init_validate_TEMPLATE

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_TEMPLATE
!
!> \brief   Setup for this initial condition
!> \author  FILL_IN_AUTHOR
!> \date    FILL_IN_DATE
!> \details
!>  This routine sets up the initial conditions for this case.
!
!-----------------------------------------------------------------------

  subroutine ocn_init_setup_TEMPLATE(domain, iErr)!{{{

    !--------------------------------------------------------------------

    type (domain_type), intent(inout) :: domain
    integer, intent(out) :: iErr

    type (block_type), pointer :: block_ptr
    type (mpas_pool_type), pointer :: meshPool
    type (mpas_pool_type), pointer :: statePool
    type (mpas_pool_type), pointer :: tracersPool
    type (mpas_pool_type), pointer :: verticalMeshPool

    ! local variables
    integer :: iCell, k, idx
    real (kind=RKIND) :: yMin, yMax, xMin, xMax, dcEdgeMin, dcEdgeMinGlobal
    real (kind=RKIND) :: yMinGlobal, yMaxGlobal, yMidGlobal, xMinGlobal, xMaxGlobal
    real (kind=RKIND) :: localVar1, localVar2
    real (kind=RKIND), dimension(:), pointer :: interfaceLocations

    ! Define config variable pointers
    character (len=StrKIND), pointer :: config_init_configuration, config_vertical_grid
    logical, pointer :: config_TEMPLATE_example_flag1
    real (kind=RKIND), pointer :: config_TEMPLATE_example_flag2

    ! Define dimension pointers
    integer, pointer :: nCellsSolve, nEdgesSolve, nVertLevels, nVertLevelsP1
    integer, pointer :: index_temperature, index_salinity

    ! Define variable pointers
    logical, pointer :: on_a_sphere
    integer, dimension(:), pointer :: maxLevelCell
    real (kind=RKIND), dimension(:), pointer :: xCell, yCell,refBottomDepth, refZMid, &
         vertCoordMovementWeights, bottomDepth, &
         fCell, fEdge, fVertex, dcEdge
    real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness
    real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers

    iErr = 0

    call mpas_pool_get_config(ocnConfigs, 'config_init_configuration', config_init_configuration)

    if(config_init_configuration .ne. trim('TEMPLATE')) return

    ! Get config flag settings

    call mpas_pool_get_config(ocnConfigs, 'config_vertical_grid', config_vertical_grid)

    call mpas_pool_get_config(ocnConfigs, 'config_TEMPLATE_example_flag1', config_TEMPLATE_example_flag1)
    call mpas_pool_get_config(ocnConfigs, 'config_TEMPLATE_example_flag2', config_TEMPLATE_example_flag2)

    ! Determine vertical grid for configuration
    call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
    call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
    call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

    ! you may restrict your case geometry as follows:
    ! if ( on_a_sphere ) call mpas_log_write('The TEMPLATE configuration can only be applied ' &
    !     // 'to a planar mesh. Exiting...', MPAS_LOG_CRIT)

    allocate(interfaceLocations(nVertLevelsP1))
    call ocn_generate_vertical_grid( config_vertical_grid, interfaceLocations )

    !--------------------------------------------------------------------
    ! Use this section to make boundaries non-periodic
    !--------------------------------------------------------------------

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

       call ocn_mark_north_boundary(meshPool, yMaxGlobal, dcEdgeMinGlobal, iErr)
       call ocn_mark_south_boundary(meshPool, yMinGlobal, dcEdgeMinGlobal, iErr)

       block_ptr => block_ptr % next
    end do

    !--------------------------------------------------------------------
    ! Use this section to set initial values
    !--------------------------------------------------------------------

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


       ! ! Set refBottomDepth and refZMid
       do k = 1, nVertLevels
          refBottomDepth(k) = config_TEMPLATE_bottom_depth * interfaceLocations(k+1)
          refZMid(k) = - 0.5_RKIND * (interfaceLocations(k+1) + interfaceLocations(k)) * config_TEMPLATE_bottom_depth
       end do

       ! Set vertCoordMovementWeights
       vertCoordMovementWeights(:) = 1.0_RKIND

       do iCell = 1, nCellsSolve

          ! Set temperature
          idx = index_temperature
          do k = 1, nVertLevels
             ! activeTracers(idx, k, iCell) =
          end do

          ! Set salinity
          idx = index_salinity
          do k = 1, nVertLevels
             ! activeTracers(idx, k, iCell) =
          end do

          ! Set layerThickness and restingThickness
          do k = 1, nVertLevels
             ! layerThickness(k, iCell) =
             ! restingThickness(k, iCell) =
          end do

          ! Set bottomDepth
          ! bottomDepth(iCell) =

          ! Set maxLevelCell
          ! maxLevelCell(iCell) =

          ! Set Coriolis parameters, if other than zero
          fCell(iCell) = config_TEMPLATE_coriolis_parameter
          fEdge(iCell) = config_TEMPLATE_coriolis_parameter
          fVertex(iCell) = config_TEMPLATE_coriolis_parameter

       end do

       block_ptr => block_ptr % next
    end do

    deallocate(interfaceLocations)
    !--------------------------------------------------------------------

  end subroutine ocn_init_setup_TEMPLATE!}}}

!***********************************************************************
!
!  routine ocn_init_validate_TEMPLATE
!
!> \brief   Validation for this initial condition
!> \author  FILL_IN_AUTHOR
!> \date    FILL_IN_DATE
!> \details
!>  This routine validates the configuration options for this case.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_validate_TEMPLATE(configPool, packagePool, iocontext, iErr)!{{{

   !--------------------------------------------------------------------
      type (mpas_pool_type), intent(inout) :: configPool, packagePool
      type (mpas_io_context_type), intent(inout) :: iocontext

      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration
      integer, pointer :: config_vert_levels, config_TEMPLATE_vert_levels

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)

      if(config_init_configuration .ne. trim('TEMPLATE')) return

      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      call mpas_pool_get_config(configPool, 'config_TEMPLATE_vert_levels', config_TEMPLATE_vert_levels)

      if(config_vert_levels <= 0 .and. config_TEMPLATE_vert_levels > 0) then
         config_vert_levels = config_TEMPLATE_vert_levels
      else if (config_vert_levels <= 0) then
         call mpas_log_write( 'Validation failed for TEMPLATE. Not given a usable value for vertical levels.', MPAS_LOG_CRIT)
         iErr = 1
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_TEMPLATE!}}}


!***********************************************************************

end module ocn_init_TEMPLATE

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
