! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_deep_convection
!
!> \brief MPAS ocean initialize case -- Deep Convection
!> \author Luke Van Roekel
!> \date   April 2021
!> \details
!>  This module initializes a deep convection test case used for  
!>  simple testing of the nonhydrostatic mode for MPAS
!>  It is a simple doubly periodic box with uniform T&S a linear EOS
!>  and a constant cooling rate
!>
!
!-----------------------------------------------------------------------

module ocn_init_deep_convection

   use mpas_kind_types
   use mpas_io_units
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants

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

   public :: ocn_init_setup_deep_convection, &
             ocn_init_validate_deep_convection

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_deep_convection
!
!> \brief   Setup for this initial condition
!> \author  Luke Van Roekel
!> \date    April 2021
!> \details
!>  This routine sets up the initial conditions for this case.
!
!-----------------------------------------------------------------------

  subroutine ocn_init_setup_deep_convection(domain, iErr)!{{{

    !--------------------------------------------------------------------

    type (domain_type), intent(inout) :: domain
    integer, intent(out) :: iErr

    type (block_type), pointer :: block_ptr
    type (mpas_pool_type), pointer :: meshPool
    type (mpas_pool_type), pointer :: statePool
    type (mpas_pool_type), pointer :: tracersPool
    type (mpas_pool_type), pointer :: verticalMeshPool
    type (mpas_pool_type), pointer :: forcingPool

    ! local variables
    integer :: iCell, k, idx, edgeCount, iNeighbor, nEdges, iEdge
    real (kind=RKIND) :: yMin, yMax, xMin, xMax, dcEdgeMin, dcEdgeMinGlobal
    real (kind=RKIND) :: yMinGlobal, yMaxGlobal, yMidGlobal, xMinGlobal, xMaxGlobal
    real (kind=RKIND) :: localVar1, localVar2, areaSum
    real (kind=RKIND), dimension(:), pointer :: interfaceLocations

    ! Define dimension pointers
    integer, pointer :: nCellsSolve, nEdgesSolve, nVertLevels, nVertLevelsP1
    integer, pointer :: index_temperature, index_salinity

    ! Define variable pointers
    logical, pointer :: on_a_sphere
    integer, dimension(:), pointer :: maxLevelCell, nEdgesOnCell
    integer, dimension(:,:), pointer :: cellsOnCell, cellMask
    real (kind=RKIND), dimension(:), pointer :: refBottomDepth, refZMid, &
         vertCoordMovementWeights, bottomDepth, &
         fCell, fEdge, fVertex, latentHeatFlux, areaCell
    real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness
    real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers

    real (kind=RKIND), dimension(:), allocatable :: sstSmooth
    real (kind=RKIND) :: botDepth, randVal

    iErr = 0

    if(config_init_configuration .ne. trim('deep_convection')) return

    ! Determine vertical grid for configuration
    call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
    call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
    call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

    ! Enforce this to be planar only 
    if ( on_a_sphere ) call mpas_log_write('The deep_convection configuration can only be applied ' &
         // 'to a planar mesh. Exiting...', MPAS_LOG_CRIT)

    allocate(interfaceLocations(nVertLevelsP1))
    call ocn_generate_vertical_grid( config_vertical_grid, interfaceLocations )

    !--------------------------------------------------------------------
    ! Use this section to set initial values
    !--------------------------------------------------------------------

    block_ptr => domain % blocklist
    do while(associated(block_ptr))
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
       call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
       call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

       call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

       call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
       call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)

       call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
       call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
       call mpas_pool_get_array(meshPool, 'cellMask', cellMask)
       call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
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

       call mpas_pool_get_array(forcingPool, 'latentHeatFlux', latentHeatFlux)

       ! ! Set refBottomDepth and refZMid
       botDepth = config_deep_convection_vert_levels*config_deep_convection_vertical_resolution
       do k = 1, nVertLevels
          refBottomDepth(k) = botDepth * interfaceLocations(k+1)
          refZMid(k) = - 0.5_RKIND * (interfaceLocations(k+1) + interfaceLocations(k)) * botDepth
       end do

       allocate(sstSmooth(nCellsSolve))
       ! Set vertCoordMovementWeights
       vertCoordMovementWeights(:) = 1.0_RKIND

       call random_seed()
       do iCell = 1, nCellsSolve

          ! Set temperature
          idx = index_temperature
          do k = 1, nVertLevels
             activeTracers(idx, k, iCell) = config_deep_convection_surface_temperature + refZmid(k)* &
                      config_deep_convection_temperature_stratification
          end do

          ! add a random perturbation to the surface
          k = 1
          call random_number(randVal)
          activeTracers(idx, k, iCell) = activeTracers(idx, k, iCell) + config_deep_convection_random_perturbation_strength*(2.0_RKIND*randVal - 1.0_RKIND)

          ! Set salinity
          idx = index_salinity
          do k = 1, nVertLevels
             activeTracers(idx, k, iCell) = config_deep_convection_salinity
          end do

          ! Set layerThickness and restingThickness
          do k = 1, nVertLevels
             layerThickness(k, iCell) = config_deep_convection_vertical_resolution
             restingThickness(k, iCell) = layerThickness(k, iCell)
          end do

          ! Set bottomDepth
          bottomDepth(iCell) = botDepth

          ! Set maxLevelCell
          maxLevelCell(iCell) = config_deep_convection_vert_levels

          ! Set Coriolis parameters, if other than zero
          fCell(iCell) = config_deep_convection_coriolis_parameter
          fEdge(iCell) = config_deep_convection_coriolis_parameter
          fVertex(iCell) = config_deep_convection_coriolis_parameter

          ! Set latent heat flux
          latentHeatFlux(iCell) = config_deep_convection_surface_cooling_strength

       end do

       !Smooth the initial noise
       do k=1,config_deep_convection_smooth_passes
          idx = index_temperature
          do iCell=1,nCellsSolve
             sstSmooth(iCell) = 0.0_RKIND
             nEdges = nEdgesOnCell(iCell)
             areaSum = 0.0_RKIND
             edgeCount = 0

             do iEdge = 1,nEdges
                iNeighbor = cellsOnCell(iEdge,iCell)
                sstSmooth(iCell) = sstSmooth(iCell) + 2.0_RKIND * areaCell(iNeighbor) * &
                                   activeTracers(idx,1,iNeighbor)
                areaSum = areaSum + 2.0_RKIND * areaCell(iNeighbor)
                edgeCount = edgeCount + 1 
             end do
             areaSum = areaSum + edgeCount*areaCell(iCell)
             sstSmooth(iCell) = sstSmooth(iCell) + activeTracers(idx, 1, iCell) * edgeCount * areaCell(iCell)
             sstSmooth(iCell) = sstSmooth(iCell) / areaSum
          end do

          do iCell = 1,nCellsSolve
             activeTracers(idx,1,iCell) = sstSmooth(iCell)
          end do
       end do

       block_ptr => block_ptr % next
    end do

    deallocate(interfaceLocations)
    !--------------------------------------------------------------------

  end subroutine ocn_init_setup_deep_convection!}}}

!***********************************************************************
!
!  routine ocn_init_validate_deep_convection
!
!> \brief   Validation for this initial condition
!> \author  Luke Van Roekel
!> \date    April 2021
!> \details
!>  This routine validates the configuration options for this case.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_validate_deep_convection(configPool, packagePool, iocontext, iErr)!{{{

   !--------------------------------------------------------------------
      type (mpas_pool_type), intent(inout) :: configPool, packagePool
      type (mpas_io_context_type), intent(inout) :: iocontext

      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration
      integer, pointer :: config_vert_levels, config_deep_convection_vert_levels

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)

      if(config_init_configuration .ne. trim('deep_convection')) return

      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      call mpas_pool_get_config(configPool, 'config_deep_convection_vert_levels', config_deep_convection_vert_levels)

      if(config_vert_levels <= 0 .and. config_deep_convection_vert_levels > 0) then
         config_vert_levels = config_deep_convection_vert_levels
      else if (config_vert_levels <= 0) then
         call mpas_log_write( 'Validation failed for deep_convection. Not given a usable value for vertical levels.', MPAS_LOG_CRIT)
         iErr = 1
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_deep_convection!}}}


!***********************************************************************

end module ocn_init_deep_convection

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
