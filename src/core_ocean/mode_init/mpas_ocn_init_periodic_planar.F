! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_periodic_planar
!
!> \brief MPAS ocean initialize case -- periodic_planar
!> \author Phillip J. Wolfram
!> \date   10/14/2015
!> \details
!>  This module contains the routines for initializing the
!>  periodic_planar initial condition, which is a constant
!>  velocity in a periodic domain.
!
!-----------------------------------------------------------------------

module ocn_init_periodic_planar

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

   public :: ocn_init_setup_periodic_planar, &
             ocn_init_validate_periodic_planar

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_periodic_planar
!
!> \brief   Setup for this initial condition
!> \author  Phillip J. Wolfram
!> \date    10/14/2015
!> \details
!>  This routine sets up the initial conditions for this case.
!
!-----------------------------------------------------------------------

  subroutine ocn_init_setup_periodic_planar(domain, iErr)!{{{

    !--------------------------------------------------------------------

    type (domain_type), intent(inout) :: domain
    integer, intent(out) :: iErr

    ! local work variables
    type (block_type), pointer :: block_ptr
    type (mpas_pool_type), pointer :: meshPool, verticalMeshPool, statePool, forcingPool, tracersPool, scratchPool

    integer :: iCell, iEdge, iVertex, k, idx
    real (kind=RKIND), dimension(:), pointer :: interfaceLocations

    ! Define config variable pointers
    character (len=StrKIND), pointer :: config_init_configuration, config_vertical_grid
    logical, pointer :: config_write_cull_cell_mask

    ! periodic_planar test case run-time configuration parameters
    real (kind=RKIND), pointer :: config_periodic_planar_bottom_depth, config_periodic_planar_velocity_strength

    integer, pointer :: config_periodic_planar_vert_levels


    ! Define dimension pointers
    integer, pointer :: nVertLevels, nCellsSolve, nEdgesSolve, nVerticesSolve, nVertLevelsP1
    integer, pointer :: index_temperature, index_salinity

    ! Define variable pointers
    logical, pointer :: on_a_sphere
    integer, dimension(:), pointer :: maxLevelCell
    integer, dimension(:,:), pointer :: verticesOnEdge
    real (kind=RKIND), dimension(:), pointer :: xCell, yCell, xEdge, yEdge, xVertex, yVertex, refBottomDepth, refZMid, &
         vertCoordMovementWeights, bottomDepth, &
         fCell, fEdge, fVertex, dcEdge, dvEdge
    real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness, normalVelocity
    real (kind=RKIND), dimension(:), pointer :: psiVertex
    type (field1DReal), pointer :: psiVertexField
    real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers

    real (kind=RKIND) :: yMin, yMax, xMin, xMax, dcEdgeMin, dcEdgeMinGlobal
    real (kind=RKIND) :: yMinGlobal, yMaxGlobal, yMidGlobal, xMinGlobal, xMaxGlobal
    real(kind=RKIND), pointer :: y_period
    character (len=StrKIND) :: streamID
    integer :: directionProperty

    ! assume no error
    iErr = 0


    ! test if periodic_planar is the desired configuration
    call mpas_pool_get_config(ocnConfigs, 'config_init_configuration', config_init_configuration)
    if(config_init_configuration .ne. trim('periodic_planar')) return

    call mpas_log_write( 'Starting initialization of planar periodic grid')

    ! get config variables !{{{
    call mpas_pool_get_config(domain % configs, 'config_write_cull_cell_mask', config_write_cull_cell_mask)
    call mpas_pool_get_config(domain % configs, 'config_periodic_planar_bottom_depth', config_periodic_planar_bottom_depth)
    call mpas_pool_get_config(domain % configs, 'config_periodic_planar_vert_levels', config_periodic_planar_vert_levels)
    call mpas_pool_get_config(domain % configs, 'config_periodic_planar_velocity_strength', &
                              config_periodic_planar_velocity_strength)
    call mpas_pool_get_config(domain % configs, 'config_vertical_grid', config_vertical_grid)
    !}}}

    ! Determine vertical grid for configuration
    call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
    call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
    call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
    call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

    ! test if configure settings are invalid
    if ( on_a_sphere ) call mpas_log_write('The planar periodic configuration can ' &
           // 'only be applied to a planar mesh. Exiting...', MPAS_LOG_CRIT)

    ! Define interface locations
    allocate(interfaceLocations(nVertLevelsP1))
    call ocn_generate_vertical_grid( config_vertical_grid, interfaceLocations )

    ! assign config variables
    nVertLevels  = config_periodic_planar_vert_levels
    nVertLevelsP1 = nVertLevels + 1

    ! keep all cells on planar, periodic mesh (no culling)

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
       call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)

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
    call mpas_pool_get_subpool(domain % blocklist % structs, 'scratch', scratchPool)
    call mpas_pool_get_field(scratchPool, 'psiVertex', psiVertexField)
    call mpas_allocate_scratch_field(psiVertexField, .false.)

    block_ptr => domain % blocklist
    do while(associated(block_ptr))
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
       call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)
       call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)

       call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
       call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
       call mpas_pool_get_dimension(meshPool, 'nVerticesSolve', nVerticesSolve)

       call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
       call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)

       call mpas_pool_get_array(meshPool, 'xCell', xCell)
       call mpas_pool_get_array(meshPool, 'yCell', yCell)
       call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
       call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
       call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
       call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
       call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
       call mpas_pool_get_array(meshPool, 'vertCoordMovementWeights', vertCoordMovementWeights)
       call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
       call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
       call mpas_pool_get_array(meshPool, 'fCell', fCell)
       call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
       call mpas_pool_get_array(meshPool, 'fVertex', fVertex)

       call mpas_pool_get_array(scratchPool, 'psiVertex', psiVertex)
       call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
       call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
       call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel=1)

       call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
       call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)

       call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)
       call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)

       ! Determine global min and max values.
       call mpas_dmpar_min_real(domain % dminfo, yMin, yMinGlobal)
       call mpas_dmpar_max_real(domain % dminfo, yMax, yMaxGlobal)
       call mpas_dmpar_min_real(domain % dminfo, xMin, xMinGlobal)
       call mpas_dmpar_max_real(domain % dminfo, xMax, xMaxGlobal)
       call mpas_dmpar_min_real(domain % dminfo, dcEdgeMin, dcEdgeMinGlobal)

       ! mark north / south boundaries
       if(config_write_cull_cell_mask) then
         call ocn_mark_north_boundary(meshPool, yMaxGlobal, dcEdgeMinGlobal, iErr)
         call ocn_mark_south_boundary(meshPool, yMinGlobal, dcEdgeMinGlobal, iErr)
         call mpas_pool_get_config(meshPool, 'y_period', y_period)
         y_period = 0.0_RKIND
       endif
       call mpas_stream_mgr_begin_iteration(domain % streamManager)
       do while (mpas_stream_mgr_get_next_stream(domain % streamManager, streamID, directionProperty))
         if ( directionProperty == MPAS_STREAM_OUTPUT .or. directionProperty == MPAS_STREAM_INPUT_OUTPUT ) then
           call mpas_stream_mgr_add_att(domain % streamManager, 'y_period', 0.0_RKIND, streamID)
         end if
       end do

       ! Set refBottomDepth and refZMid
       do k = 1, nVertLevels
          refBottomDepth(k) = config_periodic_planar_bottom_depth * interfaceLocations(k+1)
          refZMid(k) = - 0.5_RKIND * (interfaceLocations(k+1) + interfaceLocations(k)) * config_periodic_planar_bottom_depth
       end do

       ! set bottomDepth and maxLevelCell !{{{{
       bottomDepth(:) = 0.0_RKIND
       do iCell = 1, nCellsSolve

         bottomDepth(iCell)  = config_periodic_planar_bottom_depth

         ! Determine maxLevelCell based on bottomDepth and refBottomDepth
         ! Also set botomDepth based on refBottomDepth, since
         ! above bottomDepth was set with continuous analytical functions,
         ! and needs to be discrete
         maxLevelCell(iCell) = nVertLevels
         if (nVertLevels > 1) then
       do k = 1, nVertLevels
             if (bottomDepth(iCell) < refBottomDepth(k)) then
               maxLevelCell(iCell) = k-1
               bottomDepth(iCell) = refBottomDepth(k-1)
               exit
             end if
       end do
         end if

       enddo ! Looping through with iCell !}}}

       ! Set vertCoordMovementWeights
       vertCoordMovementWeights(:) = 1.0_RKIND

       do iCell = 1, nCellsSolve

          ! Set initial temperature
          idx = index_temperature
          do k = 1, nVertLevels
             activeTracers(idx, k, iCell) = 0.0_RKIND
          end do

          ! Set initial salinity
          idx = index_salinity
          do k = 1, nVertLevels
             activeTracers(idx, k, iCell) = 0.0_RKIND
          end do

          ! Set layerThickness and restingThickness
          ! Uniform layer thickness
          do k = 1, nVertLevels
            layerThickness(k, iCell) = config_periodic_planar_bottom_depth * ( interfaceLocations(k+1) - interfaceLocations(k) )
            restingThickness(k, iCell) = layerThickness(k, iCell)
          end do

          ! Set bottomDepth (above)

          ! Set maxLevelCell (above)

       end do  ! do iCell

          ! Set Coriolis parameters, if other than zero
       do iCell = 1, nCellsSolve
          fCell(iCell) = 0.0_RKIND
       end do
       do iEdge = 1, nEdgesSolve
          fEdge(iEdge) = 0.0_RKIND
       end do
       do iVertex = 1, nVerticesSolve
          fVertex(iVertex) = 0.0_RKIND
       end do

       ! Setup stream function for velocity
       do iVertex = 1, nVerticesSolve ! need to loop over all vertices to ensure correct value for edges
         psiVertex(iVertex) = yVertex(iVertex)*config_periodic_planar_velocity_strength
       end do

       !boundaryVertex => block_ptr % mesh % boundaryVertex % array(1,:)
       !block_ptr % scratch % psiVertex % array = &
       !  boundaryVertex * &
       !  sum(boundaryVertex * block_ptr % scratch % psiVertex % array) &
       !  /sum(boundaryVertex) &
       !  + (1-boundaryVertex) * block_ptr % scratch % psiVertex % array

       ! Define normalVelocity as (grad psiVertex)
       do iEdge = 1, nEdgesSolve
         normalVelocity(:,iEdge) = -1.0_RKIND * (psiVertex(verticesOnEdge(1, iEdge)) &
                                 - psiVertex(verticesOnEdge(2, iEdge)))/dvEdge(iEdge)
       end do

       block_ptr => block_ptr % next
    end do  ! do while(associated(block_ptr))
    call mpas_deallocate_scratch_field(psiVertexField, .false.)

    call mpas_log_write( 'Finishing initialization of periodic_planar')
    !--------------------------------------------------------------------

  end subroutine ocn_init_setup_periodic_planar!}}}

!***********************************************************************
!
!  routine ocn_init_validate_periodic_planar
!
!> \brief   Validation for this initial condition
!> \author  Phillip J. Wolfram
!> \date    10/14/2015
!> \details
!>  This routine validates the configuration options for this case.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_validate_periodic_planar(configPool, packagePool, iocontext, iErr)!{{{

   !--------------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: configPool, packagePool
      type (mpas_io_context_type), intent(inout) :: iocontext
      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration
      integer, pointer :: config_vert_levels, config_periodic_planar_vert_levels

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)
      if(config_init_configuration .ne. trim('periodic_planar')) return

      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      call mpas_pool_get_config(configPool, 'config_periodic_planar_vert_levels', config_periodic_planar_vert_levels)

      if(config_vert_levels <= 0 .and. config_periodic_planar_vert_levels > 0) then
         config_vert_levels = config_periodic_planar_vert_levels
      else if (config_vert_levels <= 0) then
         call mpas_log_write( 'Validation failed for periodic_planar. Not given a usable value for vertical levels.', MPAS_LOG_CRIT)
         iErr = 1
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_periodic_planar!}}}


!***********************************************************************

end module ocn_init_periodic_planar

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
