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
   use ocn_gm
   use ocn_constants

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
      ocn_init_routines_compute_max_level, &
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
!  routine ocn_init_routines_compute_max_level
!
!> \brief  initialize max level and boundary mask variables
!> \author Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date   September 2011
!> \details
!>  This routine initializes max level and boundary mask variables
!
!-----------------------------------------------------------------------
subroutine ocn_init_routines_compute_max_level(domain)!{{{
! Initialize maxLevel and boundary mesh variables.

   type (domain_type), intent(inout) :: domain
   type (mpas_pool_type), pointer :: meshPool

   integer :: i, iCell, iEdge, iVertex, k
   type (block_type), pointer :: block

   integer, pointer :: nCells, nEdges, nVertices, nVertLevels, vertexDegree

   integer, dimension(:), pointer :: &
      maxLevelCell, maxLevelEdgeTop, maxLevelEdgeBot, &
      maxLevelVertexTop, maxLevelVertexBot
   integer, dimension(:,:), pointer :: &
      cellsOnEdge, cellsOnVertex, boundaryEdge, boundaryCell, &
      boundaryVertex, verticesOnEdge, edgeMask, cellMask, vertexMask

   ! Initialize z-level mesh variables from h, read in from input file.
   block => domain % blocklist
   do while (associated(block))
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeBot', maxLevelEdgeBot)
      call mpas_pool_get_array(meshPool, 'maxLevelVertexTop', maxLevelVertexTop)
      call mpas_pool_get_array(meshPool, 'maxLevelVertexBot', maxLevelVertexBot)
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

      ! maxLevelEdgeTop is the minimum (shallowest) of the surrounding cells
      do iEdge = 1, nEdges
         maxLevelEdgeTop(iEdge) = &
            min( maxLevelCell(cellsOnEdge(1,iEdge)), &
                 maxLevelCell(cellsOnEdge(2,iEdge)) )
      end do

      maxLevelEdgeTop(nEdges+1) = 0

      ! maxLevelEdgeBot is the maximum (deepest) of the surrounding cells
      do iEdge = 1, nEdges
         maxLevelEdgeBot(iEdge) = &
            max( maxLevelCell(cellsOnEdge(1,iEdge)), &
                 maxLevelCell(cellsOnEdge(2,iEdge)) )
      end do

      maxLevelEdgeBot(nEdges+1) = 0

      ! maxLevelVertexBot is the maximum (deepest) of the surrounding cells
      do iVertex = 1,nVertices
         maxLevelVertexBot(iVertex) = maxLevelCell(cellsOnVertex(1,iVertex))
         do i = 2, vertexDegree
            maxLevelVertexBot(iVertex) = &
               max( maxLevelVertexBot(iVertex), &
                    maxLevelCell(cellsOnVertex(i,iVertex)))
         end do
      end do

      maxLevelVertexBot(nVertices+1) = 0

      ! maxLevelVertexTop is the minimum (shallowest) of the surrounding cells
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
         boundaryEdge(1:maxLevelEdgeTop(iEdge),iEdge)=0
         edgeMask(1:maxLevelEdgeTop(iEdge),iEdge)=1
      end do

      !
      ! Find cells and vertices that have an edge on the boundary
      !
      boundaryCell(:,1:nCells+1) = 0
      cellMask(:,1:nCells+1) = 0
      boundaryVertex(:,1:nVertices+1) = 0
      vertexMask(:,1:nVertices+1) = 0


      do iEdge = 1, nEdges
         do k = 1, nVertLevels
            if (boundaryEdge(k,iEdge).eq.1) then
               boundaryCell(k,cellsOnEdge(1,iEdge)) = 1
               boundaryCell(k,cellsOnEdge(2,iEdge)) = 1
               boundaryVertex(k,verticesOnEdge(1,iEdge)) = 1
               boundaryVertex(k,verticesOnEdge(2,iEdge)) = 1
            endif
         end do
      end do

      do iCell = 1, nCells
         do k = 1, nVertLevels
            if ( maxLevelCell(iCell) >= k ) then
               cellMask(k, iCell) = 1
            end if
         end do
      end do

      do iVertex = 1, nVertices
         do k = 1, nVertLevels
            if ( maxLevelVertexBot(iVertex) >= k ) then
               vertexMask(k, iVertex) = 1
            end if
         end do
      end do

      block => block % next
   end do

   ! Note: We do not update halos on maxLevel* variables.  I want the
   ! outside edge of a halo to be zero on each processor.

end subroutine ocn_init_routines_compute_max_level!}}}

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
   subroutine ocn_init_routines_setup_sign_and_index_fields(meshPool)!{{{

       type (mpas_pool_type), intent(inout) :: meshPool

       integer, dimension(:), pointer :: nEdgesOnCell
       integer, dimension(:,:), pointer :: edgesOnCell, edgesOnVertex, cellsOnVertex, cellsOnEdge, verticesOnCell, verticesOnEdge
       integer, dimension(:,:), pointer :: edgeSignOnCell, edgeSignOnVertex, kiteIndexOnCell

       integer, pointer :: nCells, nEdges, nVertices, vertexDegree
       integer :: iCell, iEdge, iVertex, i, j, k

       call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
       call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
       call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
       call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)

       call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
       call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
       call mpas_pool_get_array(meshPool, 'edgesOnVertex', edgesOnVertex)
       call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)
       call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
       call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
       call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
       call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
       call mpas_pool_get_array(meshPool, 'edgeSignOnVertex', edgeSignOnVertex)
       call mpas_pool_get_array(meshPool, 'kiteIndexOnCell', kiteIndexOnCell)

       edgeSignOnCell = 0.0_RKIND
       edgeSignOnVertex = 0.0_RKIND
       kiteIndexOnCell = 0.0_RKIND

       do iCell = 1, nCells
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

       do iVertex = 1, nVertices
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
!  routine ocn_init_routines_compute_mesh_scaling
!
!> \brief   set up mesh scaling variables
!> \author Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date   September 2011
!> \details
!>  This routine initializes meshScaling, meshScalingDel2, and
!>   meshScalingDel4
!
!-----------------------------------------------------------------------
   subroutine ocn_init_routines_compute_mesh_scaling(meshPool, scaleHmixWithMesh, maxMeshDensity)!{{{

      type (mpas_pool_type), intent(inout) :: meshPool
      logical, intent(in) :: scaleHmixWithMesh
      real (kind=RKIND), intent(in) :: maxMeshDensity

      integer :: iEdge, cell1, cell2
      integer, pointer :: nEdges
      integer, dimension(:,:), pointer :: cellsOnEdge
      real (kind=RKIND), dimension(:), pointer :: meshDensity, meshScalingDel2, meshScalingDel4, meshScaling

      call mpas_pool_get_array(meshPool, 'meshDensity', meshDensity)
      call mpas_pool_get_array(meshPool, 'meshScalingDel2', meshScalingDel2)
      call mpas_pool_get_array(meshPool, 'meshScalingDel4', meshScalingDel4)
      call mpas_pool_get_array(meshPool, 'meshScaling', meshScaling)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)

      !
      ! Compute the scaling factors to be used in the del2 and del4 dissipation
      !
      ! Typical use cases have the minval(meshScaling)==1.
      ! meshScaling values of approximately 1 indicate the highest resolution of the domain.

      meshScalingDel2(:) = 1.0_RKIND
      meshScalingDel4(:) = 1.0_RKIND
      meshScaling(:)     = 1.0_RKIND

      if (scaleHmixWithMesh) then
         do iEdge = 1, nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            meshScalingDel2(iEdge) = 1.0_RKIND / ( ((meshDensity(cell1) + meshDensity(cell2) ) / 2.0_RKIND) &
                                   / maxMeshDensity)**(3.0_RKIND / 4.0_RKIND)  ! goes as dc**3
            meshScalingDel4(iEdge) = 1.0_RKIND / ( ((meshDensity(cell1) + meshDensity(cell2) ) / 2.0_RKIND) &
                                   / maxMeshDensity)**(3.0_RKIND / 4.0_RKIND)  ! goes as dc**3
            meshScaling(iEdge)     = 1.0_RKIND / ( ((meshDensity(cell1) + meshDensity(cell2) ) / 2.0_RKIND) &
                                   / maxMeshDensity)**(1.0_RKIND / 4.0_RKIND)
         end do
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
!>  This routine initializes vertical coordinate variables
!
!-----------------------------------------------------------------------
   subroutine ocn_init_routines_vert_coord(domain)!{{{
   ! Initialize zlevel-type variables and adjust initial conditions for
   ! partial bottom cells.

      type (domain_type), intent(inout) :: domain

      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: verticalMeshPool
      type (dm_info) :: dminfo

      integer :: i, iCell, iEdge, iVertex, k, km1
      type (block_type), pointer :: block

      integer :: iTracer, cell, cell1, cell2
      real (kind=RKIND) :: normalThicknessFluxSum, thicknessSum, hEdge1, zMidPBC

      integer, dimension(:), pointer :: maxLevelCell
      real (kind=RKIND), dimension(:), pointer :: refBottomDepth, &
         refBottomDepthTopOfCell, vertCoordMovementWeights, bottomDepth, refZMid, refLayerThickness
      real (kind=RKIND), dimension(:), allocatable :: minBottomDepth, minBottomDepthMid, zMidZLevel

      real (kind=RKIND), dimension(:,:), pointer :: layerThickness
      real (kind=RKIND), dimension(:,:,:), pointer :: tracersGroup
      integer, pointer :: nVertLevels, nCells
      logical :: consistentSSH

      logical, pointer :: config_do_restart, config_check_ssh_consistency
      logical, pointer :: config_check_zlevel_consistency
      character (len=StrKIND), pointer :: config_vert_coord_movement

      type (mpas_pool_iterator_type) :: groupItr

      call mpas_pool_get_config(domain % configs, 'config_vert_coord_movement', config_vert_coord_movement)
      call mpas_pool_get_config(domain % configs, 'config_do_restart', config_do_restart)
      call mpas_pool_get_config(domain % configs, 'config_check_ssh_consistency', config_check_ssh_consistency)
      call mpas_pool_get_config(domain % configs, 'config_check_zlevel_consistency', config_check_zlevel_consistency)

      ! Initialize z-level mesh variables from h, read in from input file.
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

         call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)

         call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
         call mpas_pool_get_array(meshPool, 'refBottomDepthTopOfCell', refBottomDepthTopOfCell)
         call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
         call mpas_pool_get_array(meshPool, 'vertCoordMovementWeights', vertCoordMovementWeights)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

         call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)
         call mpas_pool_get_array(verticalMeshPool, 'refLayerThickness', refLayerThickness)

         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

         ! TopOfCell needed where zero depth for the very top may be referenced.
         refBottomDepthTopOfCell(1) = 0.0_RKIND
         do k = 1, nVertLevels
            refBottomDepthTopOfCell(k+1) = refBottomDepth(k)
            refLayerThickness(k) = refBottomDepth(k) - refBottomDepthTopOfCell(k)
            refZMid(k) = - refBottomDepthTopOfCell(k) - refLayerThickness(k)/2.0_RKIND
         end do

         ! Initialization of vertCoordMovementWeights. This determines how SSH perturbations
         ! are distributed throughout the column.
         if (config_vert_coord_movement.eq.'fixed') then

           vertCoordMovementWeights = 0.0_RKIND
           vertCoordMovementWeights(1) = 1.0_RKIND

         elseif (config_vert_coord_movement.eq.'uniform_stretching') then

            vertCoordMovementWeights = 1.0_RKIND

         endif

         if (config_check_ssh_consistency) then
            consistentSSH = .true.
            do iCell = 1,nCells
               ! Check if abs(ssh)>2m.  If so, print warning.
               if (abs(sum(layerThickness(1:maxLevelCell(iCell),iCell))-bottomDepth(iCell))>2.0_RKIND) then
                  consistentSSH = .false.
                  write (stderrUnit,'(a)') ' Warning: abs(sum(h)-bottomDepth)>2m.  Most likely, initial layerThickness ' &
                                        // 'does not match bottomDepth.'
                  write (stderrUnit,*) ' iCell, K=maxLevelCell(iCell), bottomDepth(iCell),sum(h),bottomDepth: ', &
                                iCell, maxLevelCell(iCell), bottomDepth(iCell), sum( &
                                layerThickness(1:maxLevelCell(iCell), iCell) ), bottomDepth(iCell), &
                                layerThickness(maxLevelCell(iCell),iCell)
               endif
            enddo

            if (.not. consistentSSH) then
               write(stderrUnit,*) 'Warning: SSH is not consistent. Most likely, initial layerThickness ' &
                                // 'does not match bottomDepth.'
            end if

         endif ! config_check_ssh_consistency

         if (config_check_zlevel_consistency) then
            do iCell = 1,nCells
               ! Check that bottomDepth and maxLevelCell match.  Some older meshs do not have the bottomDepth variable.
               if (bottomDepth(iCell) > refBottomDepth(maxLevelCell(iCell)).or. &
                   bottomDepth(iCell) < refBottomDepthTopOfCell(maxLevelCell(iCell))) then
                  write (stderrUnit,'(a)') ' fatal error: bottomDepth and maxLevelCell do not match:'
                  write (stderrUnit,'(a,2i5,10f10.2)') ' iCell, maxLevelCell(iCell), bottomDepth(iCell): ', &
                                iCell, maxLevelCell(iCell), bottomDepth(iCell)
                  write (stderrUnit,'(a,10f10.2)') ' refBottomDepth(maxLevelCell(iCell)), ' &
                                                   // 'refBottomDepthTopOfCell(maxLevelCell(iCell)): ', &
                                refBottomDepth(maxLevelCell(iCell)), refBottomDepthTopOfCell(maxLevelCell(iCell))
                  call mpas_dmpar_global_abort('MPAS-ocean: Abort: bottomDepth and maxLevelCell do not match')
               endif

            enddo
         endif

      block => block % next
      end do

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

      type (mpas_pool_type), pointer :: meshPool, statePool, tracersPool
      type (mpas_pool_type), pointer :: forcingPool, diagnosticsPool, scratchPool
      integer :: i, iEdge, iCell, k
      integer :: err1

      integer, dimension(:), pointer :: nAdvCellsForEdge, maxLevelCell
      integer, dimension(:), pointer :: maxLevelEdgeBot, maxLevelEdgeTop
      integer, dimension(:,:), pointer :: advCellsForEdge, highOrderAdvectionMask, boundaryCell
      real (kind=RKIND), dimension(:), pointer :: areaCell, boundaryLayerDepth
      real (kind=RKIND), dimension(:,:), pointer :: advCoefs, advCoefs3rd, normalTransportVelocity
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity, normalGMBolusVelocity, edgeTangentVectors
      real (kind=RKIND), dimension(:,:), pointer :: velocityX, velocityY, velocityZ
      real (kind=RKIND), dimension(:,:), pointer :: velocityZonal, velocityMeridional
      real (kind=RKIND), dimension(:,:,:), pointer :: derivTwo

      real (kind=RKIND), dimension(:,:,:), pointer :: tracersGroup

      integer, pointer :: nCells, nEdges, nVertices, nVertLevels
      integer, pointer :: config_horiz_tracer_adv_order
      logical, pointer :: config_hmix_scaleWithMesh, config_do_restart
      logical, pointer :: config_use_standardGM
      real (kind=RKIND), pointer :: config_maxMeshDensity

      type (mpas_pool_iterator_type) :: groupItr

      call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
      call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
      call mpas_pool_get_dimension(block % dimensions, 'nVertices', nVertices)
      call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)

      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
      call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      call mpas_pool_get_array(meshPool, 'derivTwo', derivTwo)
      call mpas_pool_get_array(meshPool, 'advCoefs', advCoefs)
      call mpas_pool_get_array(meshPool, 'advCoefs3rd', advCoefs3rd)
      call mpas_pool_get_array(meshPool, 'nAdvCellsForEdge', nAdvCellsForEdge)
      call mpas_pool_get_array(meshPool, 'advCellsForEdge', advCellsForEdge)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'highOrderAdvectionMask', highOrderAdvectionMask)
      call mpas_pool_get_array(meshPool, 'boundaryCell', boundaryCell)
      call mpas_pool_get_array(meshPool, 'edgeTangentVectors', edgeTangentVectors)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'boundaryCell', boundaryCell)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeBot', maxLevelEdgeBot)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)

      call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)
      call mpas_pool_get_array(diagnosticsPool, 'normalGMBolusVelocity', normalGMBolusVelocity)
      call mpas_pool_get_array(diagnosticsPool, 'velocityX', velocityX)
      call mpas_pool_get_array(diagnosticsPool, 'velocityY', velocityY)
      call mpas_pool_get_array(diagnosticsPool, 'velocityZ', velocityZ)
      call mpas_pool_get_array(diagnosticsPool, 'velocityZonal', velocityZonal)
      call mpas_pool_get_array(diagnosticsPool, 'velocityMeridional', velocityMeridional)
      call mpas_pool_get_array(diagnosticsPool, 'boundaryLayerDepth', boundaryLayerDepth)

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, 1)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)

      call mpas_pool_get_config(block % configs, 'config_horiz_tracer_adv_order', config_horiz_tracer_adv_order)
      call mpas_pool_get_config(block % configs, 'config_hmix_scaleWithMesh', config_hmix_scaleWithMesh)
      call mpas_pool_get_config(block % configs, 'config_maxMeshDensity', config_maxMeshDensity)
      call mpas_pool_get_config(block % configs, 'config_use_standardGM', config_use_standardGM)
      call mpas_pool_get_config(block % configs, 'config_do_restart', config_do_restart)

      call ocn_init_routines_setup_sign_and_index_fields(meshPool)
      call mpas_initialize_deriv_two(meshPool, derivTwo, err)
      call mpas_tracer_advection_coefficients(meshPool, &
          config_horiz_tracer_adv_order, derivTwo, advCoefs, &
          advCoefs3rd, nAdvCellsForEdge, advCellsForEdge, &
          err1, maxLevelCell, highOrderAdvectionMask, &
          boundaryCell)
      err = ior(err, err1)

      if (.not. config_do_restart) then
         do iCell=1,nCells
            boundaryLayerDepth(iCell) = layerThickness(1, iCell) * 0.5_RKIND
         end do
      end if

      call ocn_diagnostic_solve(dt,  statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, tracersPool)

      ! initialize velocities and active tracers on land to be zero.
      areaCell(nCells+1) = -1.0e34_RKIND

      layerThickness(:, nCells+1) = 0.0_RKIND


      do iEdge=1, nEdges
         normalVelocity(maxLevelEdgeTop(iEdge)+1:maxLevelEdgeBot(iEdge), iEdge) = 0.0_RKIND

         normalVelocity(maxLevelEdgeBot(iEdge)+1:nVertLevels,iEdge) = -1.0e34_RKIND
      end do

      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
         if ( groupItr % memberType == MPAS_POOL_FIELD ) then
            call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroup, 1)
            if ( associated(tracersGroup) ) then
               do iCell=1,nCells
                  tracersGroup(:, maxLevelCell(iCell)+1:nVertLevels,iCell) =  -1.0e34_RKIND
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


      ! Compute normalGMBolusVelocity, relativeSlope and RediDiffVertCoef if respective flags are turned on
      if (config_use_standardGM) then
          call ocn_gm_compute_Bolus_velocity(diagnosticsPool, meshPool, scratchPool)
      end if

      if (config_use_standardGM) then
         do iEdge = 1, nEdges
            normalTransportVelocity(:, iEdge) = normalTransportVelocity(:, iEdge) + normalGMBolusVelocity(:, iEdge)
         end do
       end if

      ! ------------------------------------------------------------------
      ! End: Accumulating various parametrizations of the transport velocity
      ! ------------------------------------------------------------------

      call ocn_init_routines_compute_mesh_scaling(meshPool, config_hmix_scaleWithMesh, config_maxMeshDensity)

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

      if (config_use_standardGM) then
         call ocn_reconstruct_gm_vectors(diagnosticsPool, meshPool)
      end if

      call mpas_pool_initialize_time_levels(statePool)

      ! compute land-ice fluxes for potential output at startup
      call ocn_forcing_build_fraction_absorbed_array(meshPool, statePool, diagnosticsPool, forcingPool, err1, 1)
      err = ior(err, err1)
      call ocn_surface_land_ice_fluxes_build_arrays(meshPool, diagnosticsPool, &
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
