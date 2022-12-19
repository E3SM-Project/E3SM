! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_dam_break
!
!> \brief MPAS ocean initialize case -- Zonally periodic Idealized Southern Ocean (dam_break)
!> \author Phillip J. Wolfram, Luke Van Roekel, Todd Ringler
!> \date   09/14/2015
!> \details
!>  This module contains the routines for initializing the
!>  dam_break initial condition.
!
!-----------------------------------------------------------------------

module ocn_init_dam_break

   use mpas_kind_types
   use mpas_io_units
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_stream_manager
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

   public :: ocn_init_setup_dam_break, &
             ocn_init_validate_dam_break
   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains
!***********************************************************************
!
!  routine ocn_init_setup_dam_break
!
!> \brief   Setup for this initial condition
!> \author  Zhendong Cao
!> \date    01/17/2019
!> \details
!>  This routine sets up the initial conditions for the dam_break configuration.
!
!-----------------------------------------------------------------------

  subroutine ocn_init_setup_dam_break(domain, iErr)!{{{

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
    real (kind=RKIND) :: Tscale, dx, dy

    ! Define config variable pointers
    character (len=StrKIND), pointer :: config_init_configuration, config_vertical_grid
    real (kind=RKIND), pointer :: config_dam_break_eta0
    real (kind=RKIND), pointer :: config_dam_break_R0
    real (kind=RKIND), pointer :: config_drying_min_cell_height
    real (kind=RKIND), pointer :: config_dam_break_Xl,config_dam_break_Yl,config_dam_break_Inlet
    real (kind=RKIND), pointer :: config_dam_break_dc
    integer, pointer :: config_dam_break_vert_levels
    logical, pointer :: config_write_cull_cell_mask
    logical, pointer :: config_use_wetting_drying
    real (kind=RKIND), parameter :: eps=1.0e-12

    ! Define dimension pointers
    integer, pointer :: nCellsSolve, nEdgesSolve, nVertLevels, nVertLevelsP1
    integer, pointer :: index_temperature, index_salinity

    ! Define variable pointers
    logical, pointer :: on_a_sphere
    integer, dimension(:), pointer :: maxLevelCell
    real (kind=RKIND), dimension(:), pointer :: ssh
    real (kind=RKIND), dimension(:), pointer :: xCell, yCell,refBottomDepth, refZMid, &
         vertCoordMovementWeights, bottomDepth, &
         fCell, fEdge, fVertex, dcEdge
    real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness
    real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers
    integer, dimension(:), pointer :: cullCell

    iErr = 0

    call mpas_pool_get_config(ocnConfigs, 'config_init_configuration', config_init_configuration)

    if(config_init_configuration .ne. trim('dam_break')) return

    ! Get config flag settings
    call mpas_pool_get_config(domain % configs, 'config_write_cull_cell_mask', config_write_cull_cell_mask)
    call mpas_pool_get_config(ocnConfigs, 'config_vertical_grid', config_vertical_grid)
    call mpas_pool_get_config(ocnConfigs, 'config_dam_break_eta0', config_dam_break_eta0)
    call mpas_pool_get_config(ocnConfigs, 'config_dam_break_R0', config_dam_break_R0)
    call mpas_pool_get_config(ocnConfigs, 'config_drying_min_cell_height', config_drying_min_cell_height)
    call mpas_pool_get_config(ocnConfigs, 'config_dam_break_Xl', config_dam_break_Xl)
    call mpas_pool_get_config(ocnConfigs, 'config_dam_break_Yl', config_dam_break_Yl)
    call mpas_pool_get_config(ocnConfigs, 'config_dam_break_Inlet', config_dam_break_Inlet)
    call mpas_pool_get_config(ocnConfigs, 'config_dam_break_dc', config_dam_break_dc)
    call mpas_pool_get_config(ocnConfigs, 'config_dam_break_vert_levels', config_dam_break_vert_levels)

    ! Determine vertical grid for configuration
    call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
    call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
    call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
    call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

    ! you may restrict your case geometry as follows:
    if ( on_a_sphere ) call mpas_log_write('The dam_break configuration can only be applied ' &
        // 'to a planar mesh. Exiting...', MPAS_LOG_CRIT)

    nVertLevels = config_dam_break_vert_levels
    nVertLevelsP1 = nVertLevels + 1
    allocate(interfaceLocations(nVertLevelsP1))
    call ocn_generate_vertical_grid( config_vertical_grid, interfaceLocations, ocnConfigs)

    ! Intialize the min/max values to large postive and negative values
    yMin = 1.0E10_RKIND
    yMax = -1.0E10_RKIND
    xMin = 1.0E10_RKIND
    xMax = -1.0E10_RKIND
    dcEdgeMin = 1.0E10_RKIND

    ! Determine local min and max values.
    block_ptr => domain % blocklist
    do while (associated(block_ptr))
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
    end do ! do while(associated(block_ptr))

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
       call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
       call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)
       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
       call mpas_pool_get_array(meshPool, 'xCell', xCell)
       call mpas_pool_get_array(meshPool, 'yCell', yCell)
       call mpas_pool_get_array(meshPool,'bottomDepth',bottomDepth)
       call mpas_pool_get_array(meshPool,'vertCoordMovementWeights',vertCoordMovementWeights)
       call mpas_pool_get_array(meshPool,'refBottomDepth',refBottomDepth)
       call mpas_pool_get_array(meshPool,'maxLevelCell',maxLevelCell)
       call mpas_pool_get_array(meshPool,'cullCell',cullCell)

       call mpas_pool_get_array(statePool,'layerThickness',layerThickness,1)
       call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
       call mpas_pool_get_array(statePool,'ssh',ssh,1)
       call mpas_pool_get_array(verticalMeshPool,'restingThickness',restingThickness)
       call mpas_pool_get_array(verticalMeshPool,'refZMid',refZMid)

       !set refBottomDepth and refZMid
       do k = 1, nVertLevels
         refBottomDepth(k) = config_dam_break_eta0 * interfaceLocations(k+1)
         refZMid(k) = -0.5_RKIND * (interfaceLocations(k+1) + interfaceLocations(k)) * config_dam_break_eta0
       enddo

       !set vertCoordMovementWeights
       vertCoordMovementWeights(:) = 1.0_RKIND

       !cullCell dam
       do iCell = 1, nCellsSolve
         cullCell(iCell)=0
         if (xCell(iCell).ge.config_dam_break_R0-2.0_RKIND*config_dam_break_dc) then
           if (yCell(iCell).le.(config_dam_break_R0-0.5_RKIND*config_dam_break_Yl-3.0_RKIND*config_dam_break_dc) .or. &
             yCell(iCell).ge.(config_dam_break_R0+0.5_RKIND*config_dam_break_Yl+3.0_RKIND*config_dam_break_dc)) then
             cullCell(iCell)=1
           endif
         endif
         !cullCell dam mouth
         if (xCell(iCell).ge.(config_dam_break_R0-11.0_RKIND/4.0_RKIND*config_dam_break_dc) .AND.   &
           xCell(iCell).le.(config_dam_break_R0-5.0_RKIND/4.0_RKIND*config_dam_break_dc)) then
           if ( yCell(iCell) .le.(config_dam_break_R0-0.5_RKIND*config_dam_break_Inlet) .OR.  &
             yCell(iCell) .ge. (config_dam_break_R0+0.5_RKIND*config_dam_break_Inlet+config_dam_break_dc)) then
             cullCell(iCell)=1
           endif
         endif

         !initial bathymetry
         bottomDepth(iCell) = config_dam_break_eta0
         if (xCell(iCell).ge.(config_dam_break_R0-2.0_RKIND*config_dam_break_dc)) then
           ssh(iCell) = config_dam_break_eta0
         else
           ssh(iCell) = config_dam_break_vert_levels*config_drying_min_cell_height + eps
         endif

         ! reorient for coordinate system in vertical
         ssh(iCell) = -bottomDepth(iCell) + ssh(iCell)

         ! set maxLevelCell
         maxLevelCell(iCell) = config_dam_break_vert_levels

         !layerThickness
         do k=1, maxLevelCell(iCell)
           layerThickness(k,iCell) = max(config_drying_min_cell_height + eps, &
             1.0_RKIND/float(maxLevelCell(iCell))*(ssh(iCell)+bottomDepth(iCell)))
           restingThickness(k,iCell) = layerThickness(k,iCell)
         enddo
       enddo
       ! Set salinity
       if ( associated(activeTracers) ) then
         do iCell = 1, nCellsSolve
           activeTracers(index_salinity, :, iCell) = 30.0_RKIND
         end do
       end if

       ! Set temperature
       if ( associated(activeTracers) ) then
         do iCell = 1, nCellsSolve
           activeTracers(index_temperature, :, iCell) = 10.0_RKIND
         end do
       end if

       ! Determine global min and max values.
       call mpas_dmpar_min_real(domain % dminfo, yMin, yMinGlobal)
       call mpas_dmpar_max_real(domain % dminfo, yMax, yMaxGlobal)
       call mpas_dmpar_min_real(domain % dminfo, xMin, xMinGlobal)
       call mpas_dmpar_max_real(domain % dminfo, xMax, xMaxGlobal)
       call mpas_dmpar_min_real(domain % dminfo, dcEdgeMin, dcEdgeMinGlobal)

       !mark periodic boundaries
       if(config_write_cull_cell_mask) then
         call ocn_mark_north_boundary(meshPool, yMaxGlobal, dcEdgeMinGlobal, iErr)
         call ocn_mark_south_boundary(meshPool, yMinGlobal, dcEdgeMinGlobal, iErr)
         call ocn_mark_east_boundary(meshPool, xMaxGlobal, dcEdgeMinGlobal, iErr)
         call ocn_mark_west_boundary(meshPool, xMinGlobal, dcEdgeMinGlobal, iErr)
       endif
       block_ptr => block_ptr % next
     end do !!!(do while(associated(block_ptr)))

    deallocate(interfaceLocations)

  end subroutine ocn_init_setup_dam_break!}}}


!***********************************************************************
!
!  routine ocn_init_validate_dam_break
!
!> \brief   Validation for this initial condition
!> \author  Phillip J. Wolfram, Luke Van Roekel, Todd Ringler
!> \date    09/14/2015
!> \details
!>  This routine validates the configuration options for this case.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_validate_dam_break(configPool, packagePool, iErr)!{{{

   !--------------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: configPool, packagePool
      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration
      integer, pointer :: config_vert_levels, config_dam_break_vert_levels

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)
      if(config_init_configuration .ne. trim('dam_break')) return

      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      call mpas_pool_get_config(configPool, 'config_dam_break_vert_levels', config_dam_break_vert_levels)

      if(config_vert_levels <= 0 .and. config_dam_break_vert_levels > 0) then
         config_vert_levels = config_dam_break_vert_levels
      else if (config_vert_levels <= 0) then
         call mpas_log_write( 'Validation failed for dam_break. Not given a usable value for vertical levels.', MPAS_LOG_CRIT)
         iErr = 1
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_dam_break!}}}


!***********************************************************************

end module ocn_init_dam_break

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
