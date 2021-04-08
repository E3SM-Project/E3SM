! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_init
!
!> \brief MPAS-ocean initialization routines.
!> \author Mark Petersen
!> \date   December 2013
!> \details
!>  This module contains routines to initialize variables at the
!>    beginning of an MPAS-Ocean simulation, or when starting the
!>    ocean analysis core.
!
!-----------------------------------------------------------------------

module ocn_init_routines

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer
   use mpas_dmpar
   use mpas_constants

   use mpas_rbf_interpolation
   use mpas_vector_operations
   use mpas_vector_reconstruction
   use mpas_tracer_advection_helpers

   use ocn_diagnostics
   use ocn_diagnostics_variables
   use ocn_gm
   use ocn_constants
   use ocn_config
   use ocn_mesh

   use ocn_surface_land_ice_fluxes
   use ocn_forcing

   interface ocn_init_add_tau_metadata
      module procedure ocn_init_add_tau_metadata_real
      module procedure ocn_init_add_tau_metadata_int
      module procedure ocn_init_add_tau_metadata_logical
      module procedure ocn_init_add_tau_metadata_character
   end interface

   private

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

   public :: &
      ocn_init_routines_compute_min_max_level, &
      ocn_init_routines_compute_mesh_scaling, &
      ocn_init_routines_setup_sign_and_index_fields, &
      ocn_init_routines_vert_coord, &
      ocn_init_routines_block, &
      ocn_init_metadata

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_routines_compute_min_max_level
!
!> \brief  initialize max level and boundary mask variables
!> \author Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date   September 2011
!> \details
!>  This routine initializes max level and boundary mask variables
!
!-----------------------------------------------------------------------

   subroutine ocn_init_routines_compute_min_max_level(domain)!{{{

      !-----------------------------------------------------------------
      ! Input/output variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: &
         domain           !< [inout] structure with all data for ocean

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: &
         iCell, iEdge, iVertex, &! loop indices for cells,edges,vertices
         i, k                    ! loop indices for neighbors, vertical

      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: meshPool

      integer, pointer :: &
         nCells, nEdges, nVertices, &! number of cells, edges, vertices
         nVertLevels,               &! number of vertical levels
         vertexDegree                ! max cell nbrs around vertex

      integer, dimension(:), pointer :: &
         minLevelCell,      &! min vert level in each cell column
         minLevelEdgeTop,   &! min vert level at top of edge
         minLevelEdgeBot,   &! min vert level at bot of edge
         minLevelVertexTop, &! min vert level at top of vertex
         minLevelVertexBot, &! min vert level at bot of vertex
         maxLevelCell,      &! max vert level in each cell column
         maxLevelEdgeTop,   &! max vert level at top of edge
         maxLevelEdgeBot,   &! max vert level at bot of edge
         maxLevelVertexTop, &! max vert level at top of vertex
         maxLevelVertexBot   ! max vert level at bot of vertex

      integer, dimension(:,:), pointer :: &
         cellsOnEdge,    &! indices for cells on an edge
         cellsOnVertex,  &! indices for cells on a vertex
         verticesOnEdge, &! indices for vertices on an edge
         boundaryEdge,   &! mask for edges    on boundaries
         boundaryCell,   &! mask for cells    on boundaries
         boundaryVertex, &! mask for vertices on boundaries
         edgeMask,       &! mask for active edge   points
         cellMask,       &! mask for active cell   points
         vertexMask       ! mask for active vertex points

      ! End preamble
      !-------------
      ! Begin code

      ! Retrieve data from block, pools
      block => domain % blocklist
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeBot', maxLevelEdgeBot)
      call mpas_pool_get_array(meshPool, 'minLevelEdgeTop', minLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'minLevelEdgeBot', minLevelEdgeBot)
      call mpas_pool_get_array(meshPool, 'maxLevelVertexTop', maxLevelVertexTop)
      call mpas_pool_get_array(meshPool, 'maxLevelVertexBot', maxLevelVertexBot)
      call mpas_pool_get_array(meshPool, 'minLevelVertexTop', minLevelVertexTop)
      call mpas_pool_get_array(meshPool, 'minLevelVertexBot', minLevelVertexBot)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'boundaryEdge', boundaryEdge)
      call mpas_pool_get_array(meshPool, 'boundaryCell', boundaryCell)
      call mpas_pool_get_array(meshPool, 'boundaryVertex', boundaryVertex)
      call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)
      call mpas_pool_get_array(meshPool, 'cellMask', cellMask)
      call mpas_pool_get_array(meshPool, 'vertexMask', vertexMask)

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nVertices ', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)

      ! minLevelEdgeTop is the minimum (shallowest) of the surrounding cells
      ! if the water column is dry, set minLevelCell outside of bounds
      ! to prevent dry cells from determining minimum
      do iCell = 1, nCells+1
         if (maxLevelCell(iCell) == 0) minLevelCell(iCell) = nVertLevels+1 
      end do
      do iEdge = 1, nEdges
         minLevelEdgeTop(iEdge) = &
            min( minLevelCell(cellsOnEdge(1,iEdge)), &
                 minLevelCell(cellsOnEdge(2,iEdge)) )
      end do
      minLevelEdgeTop(nEdges+1) = 0

      ! minLevelEdgeBot is the maximum (deepest) of surrounding cells
      ! prevent dry cells from determining maximum
      do iCell = 1, nCells+1
         if (maxLevelCell(iCell) == 0) minLevelCell(iCell) = 1 
      end do
      do iEdge = 1, nEdges
         minLevelEdgeBot(iEdge) = &
            max( minLevelCell(cellsOnEdge(1,iEdge)), &
                 minLevelCell(cellsOnEdge(2,iEdge)) )
      end do
      minLevelEdgeBot(nEdges+1) = 0

      ! minLevelVertexBot is the maximum (deepest) of surrounding cells
      do iVertex = 1,nVertices
         minLevelVertexBot(iVertex) = &
              minLevelCell(cellsOnVertex(1,iVertex))
         do i = 2, vertexDegree
            minLevelVertexBot(iVertex) = &
               max( minLevelVertexBot(iVertex), &
                    minLevelCell(cellsOnVertex(i,iVertex)))
         end do
      end do
      minLevelVertexBot(nVertices+1) = 0

      ! minLevelVertexTop is the minimum (shallowest) of surrounding cells
      ! if the water column is dry, set minLevelCell outside of bounds
      ! to prevent dry cells from determining minimum
      do iCell = 1, nCells+1
         if (maxLevelCell(iCell) == 0) minLevelCell(iCell) = nVertLevels+1 
      end do
      do iVertex = 1,nVertices
         minLevelVertexTop(iVertex) = minLevelCell(cellsOnVertex(1,iVertex))
         do i = 2, vertexDegree
            minLevelVertexTop(iVertex) = &
               min( minLevelVertexTop(iVertex), &
                    minLevelCell(cellsOnVertex(i,iVertex)))
         end do
      end do
      minLevelVertexTop(nVertices+1) = 0

      ! if the water column is dry, set minLevelCell = 1
      do iCell = 1, nCells+1
         if (maxLevelCell(iCell) == 0) minLevelCell(iCell) = 1
      end do
      
      ! maxLevelEdgeTop is the minimum (shallowest) of surrounding cells
      do iEdge = 1, nEdges
         maxLevelEdgeTop(iEdge) = &
            min( maxLevelCell(cellsOnEdge(1,iEdge)), &
                 maxLevelCell(cellsOnEdge(2,iEdge)) )
      end do
      maxLevelEdgeTop(nEdges+1) = 0

      ! maxLevelEdgeBot is the maximum (deepest) of surrounding cells
      do iEdge = 1, nEdges
         maxLevelEdgeBot(iEdge) = &
            max( maxLevelCell(cellsOnEdge(1,iEdge)), &
                 maxLevelCell(cellsOnEdge(2,iEdge)) )
      end do
      maxLevelEdgeBot(nEdges+1) = 0

      ! maxLevelVertexBot is the maximum (deepest) of surrounding cells
      do iVertex = 1,nVertices
         maxLevelVertexBot(iVertex) = maxLevelCell(cellsOnVertex(1,iVertex))
         do i = 2, vertexDegree
            maxLevelVertexBot(iVertex) = &
               max( maxLevelVertexBot(iVertex), &
                    maxLevelCell(cellsOnVertex(i,iVertex)))
         end do
      end do
      maxLevelVertexBot(nVertices+1) = 0

      ! maxLevelVertexTop is the minimum (shallowest) of surrounding cells
      do iVertex = 1,nVertices
         maxLevelVertexTop(iVertex) = maxLevelCell(cellsOnVertex(1,iVertex))
         do i = 2, vertexDegree
            maxLevelVertexTop(iVertex) = &
               min( maxLevelVertexTop(iVertex), &
                    maxLevelCell(cellsOnVertex(i,iVertex)))
         end do
      end do
      maxLevelVertexTop(nVertices+1) = 0

      ! set boundary edge
      boundaryEdge(:,1:nEdges+1)=1
      edgeMask(:,1:nEdges+1)=0

      do iEdge = 1, nEdges
      do k= minLevelEdgeBot(iEdge),maxLevelEdgeTop(iEdge)
         boundaryEdge(k,iEdge)=0
         edgeMask(k,iEdge)=1
      end do
      end do

      ! Find cells and vertices that have an edge on the boundary
      boundaryCell(:,1:nCells+1) = 0
      cellMask(:,1:nCells+1) = 0
      boundaryVertex(:,1:nVertices+1) = 0
      vertexMask(:,1:nVertices+1) = 0

      do iEdge = 1, nEdges
         do k = 1, nVertLevels
            if (boundaryEdge(k,iEdge) == 1) then
               boundaryCell(k,cellsOnEdge(1,iEdge)) = 1
               boundaryCell(k,cellsOnEdge(2,iEdge)) = 1
               boundaryVertex(k,verticesOnEdge(1,iEdge)) = 1
               boundaryVertex(k,verticesOnEdge(2,iEdge)) = 1
            endif
         end do
      end do

      do iCell = 1, nCells
         do k = 1, nVertLevels
            if ( minLevelCell(iCell) <= k .and. &
                 maxLevelCell(iCell) >= k ) then
               cellMask(k, iCell) = 1
            end if
         end do
      end do

      do iVertex = 1, nVertices
         do k = 1, nVertLevels
            if ( minLevelVertexTop(iVertex) <= k .and. &
                 maxLevelVertexBot(iVertex) >= k ) then
               vertexMask(k, iVertex) = 1
            end if
         end do
      end do


      ! Note: We do not update halos on maxLevel* variables.  We want
      ! the outside edge of a halo to be zero on each processor.

   !--------------------------------------------------------------------

end subroutine ocn_init_routines_compute_min_max_level!}}}

!***********************************************************************
!
!  routine ocn_init_routines_setup_sign_and_index_fields
!
!> \brief   set up sign and index fields
!> \author Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date   September 2011
!> \details
!>  This routine initializes edgeSignOnCell, edgeSignOnVertex, and
!>   kiteIndexOnCell.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_routines_setup_sign_and_index_fields()!{{{

       integer :: iCell, iEdge, iVertex, i, j, k

       edgeSignOnCell = 0.0_RKIND
       edgeSignOnVertex = 0.0_RKIND
       kiteIndexOnCell = 0.0_RKIND

       do iCell = 1, nCellsAll
         do i = 1, nEdgesOnCell(iCell)
           iEdge = edgesOnCell(i, iCell)
           iVertex = verticesOnCell(i, iCell)

           ! Vector points from cell 1 to cell 2
           if(iCell == cellsOnEdge(1, iEdge)) then
             edgeSignOnCell(i, iCell) = -1
           else
             edgeSignOnCell(i, iCell) =  1
           end if

           do j = 1, vertexDegree
             if(cellsOnVertex(j, iVertex) == iCell) then
               kiteIndexOnCell(i, iCell) = j
             end if
           end do
         end do
       end do

       do iVertex = 1, nVerticesAll
         do i = 1, vertexDegree
           iEdge = edgesOnVertex(i, iVertex)

           ! Vector points from vertex 1 to vertex 2
           if(iVertex == verticesOnEdge(1, iEdge)) then
             edgeSignOnVertex(i, iVertex) = -1
           else
             edgeSignOnVertex(i, iVertex) =  1
           end if
         end do
       end do

   end subroutine ocn_init_routines_setup_sign_and_index_fields!}}}

!***********************************************************************
!
!  routine ocn_init_routines_area_weights
!
!> \brief  set up area weighting
!> \author Mark Petersen
!> \date   June 2020
!> \details
!>  This routine initializes edgeAreaFractionOfCell
!
!-----------------------------------------------------------------------

   subroutine ocn_init_routines_area_weights(meshPool, edgeAreaFractionOfCell)!{{{

      type(mpas_pool_type), intent(inout) :: meshPool
      real (kind=RKIND), dimension(:,:), intent(inout) :: edgeAreaFractionOfCell

      integer :: iCell, iEdge, i, j, k

      do iCell = 1, nCellsAll
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            ! edgeAreaFractionOfCell is the fractional area of cell iCell
            ! encompassed by the triangle with edge iEdge connected to the cell
            ! center of cell iCell. On a perfect planar hex mesh this is always
            ! 1/6. This weighting is used to interpolate scalars from edges to
            ! cell centers. Use 2*edgeAreaFractionOfCell to interpolate normal
            ! vectors at edges to vector norms at cell centers."
            edgeAreaFractionOfCell(i, iCell) = &
               0.25_RKIND * dcEdge(iEdge) * dvEdge(iEdge) / areaCell(iCell)
         end do
      end do

   end subroutine ocn_init_routines_area_weights!}}}

!***********************************************************************
!
!  routine ocn_init_routines_compute_mesh_scaling
!
!> \brief   set up mesh scaling variables
!> \author Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date   September 2011
!> \details
!>  This routine initializes meshScalingDel2 and meshScalingDel4
!
!-----------------------------------------------------------------------
   subroutine ocn_init_routines_compute_mesh_scaling(meshPool, & !{{{
                 config_hmix_scaleWithMesh, config_maxMeshDensity, &
                 config_hmix_use_ref_cell_width, config_hmix_ref_cell_width)

      type (mpas_pool_type), intent(inout) :: meshPool
      logical, intent(in) :: config_hmix_scaleWithMesh
      real (kind=RKIND), intent(in) :: config_maxMeshDensity
      logical, intent(in) :: config_hmix_use_ref_cell_width
      real (kind=RKIND), intent(in) :: config_hmix_ref_cell_width

      integer :: iEdge, cell1, cell2
      integer, pointer :: nEdges
      integer, dimension(:,:), pointer :: cellsOnEdge
      real (kind=RKIND) :: cellWidth
      real (kind=RKIND), dimension(:), pointer :: meshDensity, meshScalingDel2, meshScalingDel4, meshScaling
      real (kind=RKIND), dimension(:), pointer :: areaCell

      call mpas_pool_get_array(meshPool, 'meshDensity', meshDensity)
      call mpas_pool_get_array(meshPool, 'meshScalingDel2', meshScalingDel2)
      call mpas_pool_get_array(meshPool, 'meshScalingDel4', meshScalingDel4)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)

      if (config_hmix_scaleWithMesh) then
         if (config_hmix_use_ref_cell_width) then
            ! Mesh scaling is set by areaCell and config_hmix_ref_cell_width
            ! See description of config_hmix_ref_cell_width in Registry.xml for
            ! more detail.
            do iEdge = 1, nEdges
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               ! Effective cell width at edge iEdge, assuming neighboring cells
               ! are circles for this calculation.
               cellWidth = 2.0_RKIND * sqrt((areaCell(cell1) + areaCell(cell2) ) / 2.0_RKIND / pii)
               meshScalingDel2(iEdge) = cellWidth / config_hmix_ref_cell_width
               meshScalingDel4(iEdge) = ( cellWidth / config_hmix_ref_cell_width )**3
            end do
         else
            ! Mesh scaling is set by meshDensity. This is both confusing and
            ! inconvenient, as the flags like config_mom_del2 need to be reset
            ! for every resolution. It is kept for backwards compatibility, but
            ! should become defunct.
            do iEdge = 1, nEdges
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               meshScalingDel2(iEdge) = 1.0_RKIND / ( ((meshDensity(cell1) + meshDensity(cell2) ) / 2.0_RKIND) &
                                      / config_maxMeshDensity)**(1.0_RKIND / 4.0_RKIND)  ! goes as dc**1
               meshScalingDel4(iEdge) = 1.0_RKIND / ( ((meshDensity(cell1) + meshDensity(cell2) ) / 2.0_RKIND) &
                                      / config_maxMeshDensity)**(3.0_RKIND / 4.0_RKIND)  ! goes as dc**3
            end do
         end if
      else
         ! If config_hmix_scaleWithMesh is false, hmix coefficients do not vary with
         ! cell size but remain constant across the domain.
         meshScalingDel2(:) = 1.0_RKIND
         meshScalingDel4(:) = 1.0_RKIND
      end if

   end subroutine ocn_init_routines_compute_mesh_scaling!}}}

!***********************************************************************
!
!  routine ocn_init_routines_vert_coord
!
!> \brief  initialize vertical coordinate variables
!> \author Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date   September 2011
!> \details
!>  This routine initializes vertical coordinate variables, including
!>  zlevel-type variables and adjusted initial conditions for partial
!>  bottom cells.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_routines_vert_coord(domain)!{{{

      !-----------------------------------------------------------------
      ! Input/output variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain !< [inout] ocean state

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      type (block_type), pointer :: block ! all ocn data in a block

      type (mpas_pool_type), pointer :: &
         meshPool,         &! mesh data
         verticalMeshPool, &! vertical mesh data
         statePool,        &! ocean state (layer thickness)
         forcingPool        ! forcing data (land ice mask)

      integer :: &
         iCell, k,         &! loop indices for cell, vertical loops
         kmin, kmax         ! min, max active levels in column

      integer, pointer ::  &
         nVertLevels,      &! number of vertical levels
         nCells             ! total number of cells

      logical :: &
         consistentSSH,    &! flag if SSH consistent with depth
         landIceFlag,      &! flag if land ice mask exists
         checkSSH           ! flag for checking SSH in cell column

      real (kind=RKIND) &
         thickDepth         ! depth based on sum of layer thicknesses

      integer, dimension(:), pointer :: &
         minLevelCell,     &! min vertical level in each cell column
         maxLevelCell,     &! max vertical level in each cell column
         landIceMask        ! land ice mask

      real (kind=RKIND), dimension(:), pointer :: &
         refBottomDepth,           &! reference bottom depth
         refBottomDepthTopOfCell,  &! reference depth at top of cell
         vertCoordMovementWeights, &! weights for distributing hgt change
         bottomDepth,              &! bottom depth
         refZMid,                  &! reference depth at mid-layer
         refLayerThickness          ! reference profile layer thickness

      real (kind=RKIND), dimension(:,:), pointer :: &
         layerThickness     ! current layer thickness

      ! End preamble
      !-------------
      ! Begin code

      ! Extract pools from domain
      block => domain % blocklist
      call mpas_pool_get_subpool(block%structs, 'state', statePool)
      call mpas_pool_get_subpool(block%structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block%structs, 'verticalMesh', &
                                                 verticalMeshPool)
      call mpas_pool_get_subpool(block%structs, 'forcing', forcingPool)

      ! Extract needed variables from pools
      call mpas_pool_get_array(statePool, 'layerThickness', &
                                           layerThickness, 1)
      call mpas_pool_get_array(meshPool,  'refBottomDepth', &
                                           refBottomDepth)
      call mpas_pool_get_array(meshPool,  'refBottomDepthTopOfCell', &
                                           refBottomDepthTopOfCell)
      call mpas_pool_get_array(meshPool,  'bottomDepth', bottomDepth)
      call mpas_pool_get_array(meshPool,  'vertCoordMovementWeights', &
                                           vertCoordMovementWeights)
      call mpas_pool_get_array(meshPool,  'minLevelCell', minLevelCell)
      call mpas_pool_get_array(meshPool,  'maxLevelCell', maxLevelCell)

      call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)
      call mpas_pool_get_array(verticalMeshPool, 'refLayerThickness', &
                                                  refLayerThickness)

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(forcingPool, 'landIceMask', landIceMask)
      if (associated(landIceMask)) then
         landIceFlag = .true.  ! ice shelves
      else
         landIceFlag = .false. ! ice shelves
      endif

      ! Initialize reference z coordinates based on refBottomDepth
      refBottomDepthTopOfCell(1) = 0.0_RKIND
      do k = 1, nVertLevels
         refBottomDepthTopOfCell(k+1) = refBottomDepth(k)
         refLayerThickness(k) = refBottomDepth(k) - &
                                refBottomDepthTopOfCell(k)
         refZMid(k) = - refBottomDepthTopOfCell(k) - &
                        refLayerThickness(k)/2.0_RKIND
      end do

      ! Initialization of vertCoordMovementWeights. 
      ! This determines how SSH perturbations are distributed throughout
      ! the column. Check for invalid choice of vert coord movement.

      call mpas_log_write(' Vertical coordinate movement is: ' // &
                          trim(config_vert_coord_movement))


      select case (trim(config_vert_coord_movement))
      case ('fixed')

         vertCoordMovementWeights = 0.0_RKIND
         vertCoordMovementWeights(1) = 1.0_RKIND

      case ('uniform_stretching')

         vertCoordMovementWeights = 1.0_RKIND

      case ('tapered')

         ! Set weight tapering:
         !  1.0 shallower than config_vert_taper_weight_depth_1
         !  linear in between
         !  0.0 deeper than config_vert_taper_weight_depth_2
         do k = 1, nVertLevels
            vertCoordMovementWeights(k) = 1.0_RKIND + &
                 (refZMid(k) + config_vert_taper_weight_depth_1) / &
                 (config_vert_taper_weight_depth_2 - &
                  config_vert_taper_weight_depth_1)
            ! limit range to (0,1)
            vertCoordMovementWeights(k) = max( 0.0_RKIND, &
                      min( 1.0_RKIND, vertCoordMovementWeights(k) ) )
         end do

      case ('impermeable_interfaces')

         ! weights are never used, no init needed

      case ('user_specified')

         ! weights should be input? probably should check

      case default

         call mpas_log_write( &
            'Incorrect choice of config_vert_coord_movement.', &
             MPAS_LOG_CRIT)

      end select

      !-----------------------------------------------------------------
      ! Check consistency and validity of options
      !-----------------------------------------------------------------

      ! SSH consistency
      if (config_check_ssh_consistency) then

         ! Check if abs(ssh)>20m.  If so, print warning.
         consistentSSH = .true.
         do iCell = 1,nCells

            ! Only check for open ocean when land ice is being used
            checkSSH = .true.
            if (landIceFlag) then
               if (landIceMask(iCell) /= 0) checkSSH = .false.
            endif

            ! Check whether thickness is consistent with bottom depth
            if (checkSSH) then
               kmin = minLevelCell(iCell)
               kmax = maxLevelCell(iCell)

               thickDepth = sum(layerThickness(kmin:kmax,iCell))

               if (abs(thickDepth - bottomDepth(iCell)) > 20.0_RKIND) then
                  consistentSSH = .false.

                  call mpas_log_write( &
                     ' Warning: Sea surface height outside of acceptable range (20m)', &
                     MPAS_LOG_ERR)
                  call mpas_log_write( &
                     ' iCell: $i, minLevelCell(iCell), maxLevelCell(iCell): $i, $i, bottomDepth(iCell): $r, sum(h): $r', &
                     intArgs=(/iCell, minLevelCell(iCell), maxLevelCell(iCell) /), &
                     realArgs=(/ bottomDepth(iCell), thickDepth /) )
               endif ! depth diff check
            endif ! checkSSH
         end do ! cell loop

         if (.not. consistentSSH) call mpas_log_write( &
            'Warning: SSH not consistent - ' // &
            'initial layerThickness does not match bottomDepth.')

      endif ! config_check_ssh_consistency

      ! Z-level consistency
      if (config_check_zlevel_consistency) then
         do iCell = 1,nCells
            ! Check that bottomDepth and maxLevelCell match.  Some
            ! older meshs do not have the bottomDepth variable.
            kmax = maxLevelCell(iCell)
            if (bottomDepth(iCell) > refBottomDepth(kmax) .or. &
                bottomDepth(iCell) < refBottomDepthTopOfCell(kmax)) then
               call mpas_log_write( &
                 ' fatal error: bottomDepth and maxLevelCell do not match:',  &
                 MPAS_LOG_ERR)
               call mpas_log_write( &
                 ' iCell: $i, maxLevelCell(iCell): $i, bottomDepth(iCell): $r', &
                 intArgs=(/iCell, kmax /), &
                 realArgs=(/ bottomDepth(iCell) /) )
            endif
         enddo ! cell loop
      endif ! check_zlevel_consistency

      ! pressure gradient compatibility
      if (config_vert_coord_movement /= 'impermeable_interfaces' .and. &
          config_pressure_gradient_type == 'MontgomeryPotential') then
         call mpas_log_write( 'Incompatible choice of ' // &
            'impermeable vertical interfaces and Montgomery Potential',&
            MPAS_LOG_CRIT)
      end if

      ! barotropic filter mode compatibility
      if (config_filter_btr_mode .and. &
          config_vert_coord_movement /= 'fixed')then
         call mpas_log_write( &
          'filter_btr_mode has only been tested with ' // &
          'config_vert_coord_movement=fixed.', &
          MPAS_LOG_CRIT)
      endif

   !--------------------------------------------------------------------

   end subroutine ocn_init_routines_vert_coord!}}}

!***********************************************************************
!
!  routine ocn_init_routines_block
!
!> \brief   Initialize blocks within MPAS-Ocean core
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details
!>  This routine calls all block-level initializations required to begin a
!>  simulation with MPAS-Ocean
!
!-----------------------------------------------------------------------

   subroutine ocn_init_routines_block(block, dt, err)!{{{

      type (block_type), intent(inout) :: block
      real (kind=RKIND), intent(in) :: dt
      integer, intent(out) :: err

      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: statePool, tracersPool
      type (mpas_pool_type), pointer :: forcingPool, scratchPool
      integer :: i, iEdge, iCell, k
      integer :: err1

      integer, dimension(:,:), pointer :: highOrderAdvectionMaskTmp, boundaryCellTmp
      integer, dimension(:,:), pointer :: edgeSignOnCellTmp, edgeSignOnVertexTmp
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity

      real (kind=RKIND), dimension(:,:,:), pointer :: tracersGroup

      integer, pointer :: nCells, nEdges, nVertices, nVertLevels

      type (mpas_pool_iterator_type) :: groupItr

      call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
      call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
      call mpas_pool_get_dimension(block % dimensions, 'nVertices', nVertices)
      call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)

      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      call mpas_pool_get_array(meshPool, 'highOrderAdvectionMask', highOrderAdvectionMaskTmp)
      call mpas_pool_get_array(meshPool, 'boundaryCell', boundaryCellTmp)

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, 1)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)

      call ocn_init_routines_setup_sign_and_index_fields()

      ! Update edgeSignOnCell and edgeSignOnVertex pionters in ocn_mesh.  The ocn_mesh allocatables
      ! are updated in ocn_init_routines_setup_sign_and_index_fields.  This updates the pointers to
      ! Registry to match the allocatable arrays
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCellTmp)
      call mpas_pool_get_array(meshPool, 'edgeSignOnVertex', edgeSignOnVertexTmp)
      edgeSignOnCellTmp = int(edgeSignOnCell)
      edgeSignOnVertexTmp = int(edgeSignOnVertex)

      call ocn_init_routines_area_weights(meshPool, edgeAreaFractionOfCell)
      call mpas_initialize_deriv_two(meshPool, derivTwo, err)
      call mpas_tracer_advection_coefficients(meshPool, &
          config_horiz_tracer_adv_order, derivTwo, advCoefs, &
          advCoefs3rd, nAdvCellsForEdge, advCellsForEdge, &
          err1, maxLevelCell, highOrderAdvectionMaskTmp, &
          boundaryCellTmp)
      err = ior(err, err1)

      ! Update ocn_mesh variables for highOrderAdvectionMask and boundaryCell
      highOrderAdvectionMask = real(highOrderAdvectionMaskTmp, kind=RKIND)
      boundaryCell = real(boundaryCellTmp, kind=RKIND)

      if (.not. config_do_restart) then
         do iCell=1,nCells
            boundaryLayerDepth(iCell) = layerThickness(minLevelCell(iCell), iCell) * 0.5_RKIND
            indMLD(iCell) = minLevelCell(iCell)
         end do
      end if

      call ocn_diagnostic_solve(dt,  statePool, forcingPool, meshPool, scratchPool, tracersPool)

      ! initialize velocities and active tracers on land to be zero.
      areaCell(nCells+1) = -1.0e34_RKIND

      layerThickness(:, nCells+1) = 0.0_RKIND


      do iEdge=1, nEdges
         normalVelocity(maxLevelEdgeTop(iEdge)+1:maxLevelEdgeBot(iEdge), iEdge) = 0.0_RKIND

         normalVelocity(maxLevelEdgeBot(iEdge)+1:nVertLevels,iEdge) = -1.0e34_RKIND
         
         normalVelocity(1:minLevelEdgeTop(iEdge)-1,iEdge) = -1.0e34_RKIND
      end do

      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
         if ( groupItr % memberType == MPAS_POOL_FIELD ) then
            call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroup, 1)
            if ( associated(tracersGroup) ) then
               do iCell=1,nCells
                  tracersGroup(:, maxLevelCell(iCell)+1:nVertLevels,iCell) =  -1.0e34_RKIND
                  tracersGroup(:, 1:minLevelCell(iCell)-1,iCell) =  -1.0e34_RKIND
               end do
            end if
         end if
      end do

      ! ------------------------------------------------------------------
      ! Accumulating various parametrizations of the transport velocity
      ! ------------------------------------------------------------------
      do iEdge = 1, nEdges
         normalTransportVelocity(:, iEdge) = normalVelocity(:, iEdge)
      end do

      ! ------------------------------------------------------------------
      ! End: Accumulating various parametrizations of the transport velocity
      ! ------------------------------------------------------------------

      call ocn_init_routines_compute_mesh_scaling(meshPool, config_hmix_scaleWithMesh, config_maxMeshDensity, &
              config_hmix_use_ref_cell_width, config_hmix_ref_cell_width)

      call mpas_rbf_interp_initialize(meshPool)
      call mpas_initialize_tangent_vectors(meshPool, edgeTangentVectors)

      call mpas_init_reconstruct(meshPool, includeHalos=.true.)

      call mpas_reconstruct(meshPool, normalVelocity,        &
                       velocityX,            &
                       velocityY,            &
                       velocityZ,            &
                       velocityZonal,        &
                       velocityMeridional    &
                      )

      call mpas_pool_initialize_time_levels(statePool)

      ! compute land-ice fluxes for potential output at startup
      call ocn_forcing_build_fraction_absorbed_array(meshPool, statePool, forcingPool, err1, 1)
      err = ior(err, err1)
      call ocn_surface_land_ice_fluxes_build_arrays(meshPool, &
                                                    forcingPool, scratchPool, statePool, dt, err1)
      err = ior(err, err1)



   end subroutine ocn_init_routines_block!}}}


!***********************************************************************
!
!  routine ocn_init_metadata
!
!> \brief   Initialize any metadata for this processor
!> \author  Doug Jacobsen
!> \date    08/05/2016
!> \details
!>  This routine sets up any metadata for this MPI task and it's associated threads.
!>  The meta data could be related to performance data, or information about
!>  all the blocks on this processor.
!
!-----------------------------------------------------------------------
   subroutine ocn_init_metadata(domain)!{{{
      type (domain_type), intent(inout) :: domain

      type (block_type), pointer :: block

      character (len=StrKIND) :: metaDataName

      integer :: iHalo
      integer :: numBlocks
      integer, dimension(:), pointer :: nCellsArray, nEdgesArray, nVerticesArray

      numBlocks = 0

      block => domain % blocklist
      do while ( associated(block) )
         numBlocks = numBlocks + 1

         call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
         call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
         call mpas_pool_get_dimension(block % dimensions, 'nVerticesArray', nVerticesArray)

         call ocn_init_add_tau_metadata( 'numCellHalos', size(nCellsArray) )
         call ocn_init_add_tau_metadata( 'numEdgeHalos', size(nEdgesArray) )
         call ocn_init_add_tau_metadata( 'numVertexHalos', size(nVerticesArray) )

         do iHalo = 1, size(nCellsArray)
            write(metaDataName, '(a8, i2)' ) 'cellHalo', iHalo
            call ocn_init_add_tau_metadata( metaDataName, nCellsArray(iHalo) )
         end do

         do iHalo = 1, size(nEdgesArray)
            write(metaDataName, '(a8, i2)' ) 'edgeHalo', iHalo
            call ocn_init_add_tau_metadata( metaDataName, nEdgesArray(iHalo) )
         end do

         do iHalo = 1, size(nVerticesArray)
            write(metaDataName, '(a10, i2)' ) 'vertexHalo', iHalo
            call ocn_init_add_tau_metadata( metaDataName, nVerticesArray(iHalo) )
         end do

         block => block % next
      end do

      call ocn_init_add_tau_metadata( 'numBlocks', numBlocks )

   end subroutine ocn_init_metadata!}}}


!***********************************************************************
!
!  routine ocn_init_add_tau_metadata_real
!
!> \brief   Add tau real metadata
!> \author  Doug Jacobsen
!> \date    08/05/2016
!> \details
!>  This routine adds a real value metadata for TAU to this task.
!
!-----------------------------------------------------------------------
   subroutine ocn_init_add_tau_metadata_real( dataName, dataValue )!{{{
      character (len=*), intent(in) :: dataName
      real (kind=RKIND), intent(in) :: dataValue

      character (len=StrKIND) :: dataString

#ifdef MPAS_TAU
      dataString = ''
      write( dataString, * ) dataValue
      call tau_metadata( trim(dataName), trim(dataString) )
#endif

   end subroutine ocn_init_add_tau_metadata_real!}}}

!***********************************************************************
!
!  routine ocn_init_add_tau_metadata_int
!
!> \brief   Add tau integer metadata
!> \author  Doug Jacobsen
!> \date    08/05/2016
!> \details
!>  This routine adds an integer value metadata for TAU to this task.
!
!-----------------------------------------------------------------------
   subroutine ocn_init_add_tau_metadata_int( dataName, dataValue )!{{{
      character (len=*), intent(in) :: dataName
      integer, intent(in) :: dataValue

      character (len=StrKIND) :: dataString

#ifdef MPAS_TAU
      dataString = ''
      write( dataString, * ) dataValue
      call tau_metadata( trim(dataName), trim(dataString) )
#endif

   end subroutine ocn_init_add_tau_metadata_int!}}}

!***********************************************************************
!
!  routine ocn_init_add_tau_metadata_logical
!
!> \brief   Add tau logical metadata
!> \author  Doug Jacobsen
!> \date    08/05/2016
!> \details
!>  This routine adds a logical value metadata for TAU to this task.
!
!-----------------------------------------------------------------------
   subroutine ocn_init_add_tau_metadata_logical( dataName, dataValue )!{{{
      character (len=*), intent(in) :: dataName
      logical, intent(in) :: dataValue

      character (len=StrKIND) :: dataString

#ifdef MPAS_TAU
      dataString = ''
      write( dataString, * ) dataValue
      call tau_metadata( trim(dataName), trim(dataString) )
#endif

   end subroutine ocn_init_add_tau_metadata_logical!}}}

!***********************************************************************
!
!  routine ocn_init_add_tau_metadata_character
!
!> \brief   Add tau character metadata
!> \author  Doug Jacobsen
!> \date    08/05/2016
!> \details
!>  This routine adds a character value metadata for TAU to this task.
!
!-----------------------------------------------------------------------
   subroutine ocn_init_add_tau_metadata_character( dataName, dataValue )!{{{
      character (len=*), intent(in) :: dataName
      character (len=*), intent(in) :: dataValue

#ifdef MPAS_TAU
      call tau_metadata( trim(dataName), trim(dataValue) )
#endif

   end subroutine ocn_init_add_tau_metadata_character!}}}



end module ocn_init_routines

! vim: foldmethod=marker
