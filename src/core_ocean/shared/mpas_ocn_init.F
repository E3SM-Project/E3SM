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

module ocn_init

   use mpas_configure
   use mpas_framework
   use mpas_dmpar

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
      ocn_compute_max_level, &
      ocn_compute_mesh_scaling, &
      ocn_setup_sign_and_index_fields, &
      ocn_init_vert_coord

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_compute_max_level
!
!> \brief  initialize max level and boundary mask variables
!> \author Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date   September 2011
!> \details 
!>  This routine initializes max level and boundary mask variables
!
!-----------------------------------------------------------------------
subroutine ocn_compute_max_level(domain)!{{{
! Initialize maxLevel and boundary mesh variables.

   use mpas_grid_types
   use mpas_configure
   use mpas_constants

   implicit none

   type (domain_type), intent(inout) :: domain

   integer :: i, iCell, iEdge, iVertex, k
   type (block_type), pointer :: block

   integer :: nCells, nEdges, nVertices, nVertLevels, vertexDegree

   integer, dimension(:), pointer :: &
      maxLevelCell, maxLevelEdgeTop, maxLevelEdgeBot, &
      maxLevelVertexTop, maxLevelVertexBot
   integer, dimension(:,:), pointer :: &
      cellsOnEdge, cellsOnVertex, boundaryEdge, boundaryCell, &
      boundaryVertex, verticesOnEdge, edgeMask, cellMask, vertexMask

   ! Initialize z-level mesh variables from h, read in from input file.
   block => domain % blocklist
   do while (associated(block))

      maxLevelCell => block % mesh % maxLevelCell % array
      maxLevelEdgeTop => block % mesh % maxLevelEdgeTop % array
      maxLevelEdgeBot => block % mesh % maxLevelEdgeBot % array
      maxLevelVertexTop => block % mesh % maxLevelVertexTop % array
      maxLevelVertexBot => block % mesh % maxLevelVertexBot % array
      cellsOnEdge    => block % mesh % cellsOnEdge % array
      cellsOnVertex  => block % mesh % cellsOnVertex % array
      verticesOnEdge => block % mesh % verticesOnEdge % array
      boundaryEdge   => block % mesh % boundaryEdge % array
      boundaryCell   => block % mesh % boundaryCell % array
      boundaryVertex => block % mesh % boundaryVertex % array
      edgeMask => block % mesh % edgeMask % array
      cellMask => block % mesh % cellMask % array
      vertexMask => block % mesh % vertexMask % array

      nCells      = block % mesh % nCells
      nEdges      = block % mesh % nEdges
      nVertices   = block % mesh % nVertices
      nVertLevels = block % mesh % nVertLevels
      vertexDegree = block % mesh % vertexDegree

      ! maxLevelEdgeTop is the minimum (shallowest) of the surrounding cells
      do iEdge=1,nEdges
         maxLevelEdgeTop(iEdge) = &
            min( maxLevelCell(cellsOnEdge(1,iEdge)), &
                 maxLevelCell(cellsOnEdge(2,iEdge)) )
      end do 
      maxLevelEdgeTop(nEdges+1) = 0

      ! maxLevelEdgeBot is the maximum (deepest) of the surrounding cells
      do iEdge=1,nEdges
         maxLevelEdgeBot(iEdge) = &
            max( maxLevelCell(cellsOnEdge(1,iEdge)), &
                 maxLevelCell(cellsOnEdge(2,iEdge)) )
      end do 
      maxLevelEdgeBot(nEdges+1) = 0

      ! maxLevelVertexBot is the maximum (deepest) of the surrounding cells
      do iVertex = 1,nVertices
         maxLevelVertexBot(iVertex) = maxLevelCell(cellsOnVertex(1,iVertex))
         do i=2,vertexDegree
            maxLevelVertexBot(iVertex) = &
               max( maxLevelVertexBot(iVertex), &
                    maxLevelCell(cellsOnVertex(i,iVertex)))
         end do
      end do 
      maxLevelVertexBot(nVertices+1) = 0

      ! maxLevelVertexTop is the minimum (shallowest) of the surrounding cells
      do iVertex = 1,nVertices
         maxLevelVertexTop(iVertex) = maxLevelCell(cellsOnVertex(1,iVertex))
         do i=2,vertexDegree
            maxLevelVertexTop(iVertex) = &
               min( maxLevelVertexTop(iVertex), &
                    maxLevelCell(cellsOnVertex(i,iVertex)))
         end do
      end do 
      maxLevelVertexTop(nVertices+1) = 0

      ! set boundary edge
      boundaryEdge(:,1:nEdges+1)=1
      edgeMask(:,1:nEdges+1)=0
      do iEdge=1,nEdges
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
      do iEdge=1,nEdges
         do k=1,nVertLevels
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

end subroutine ocn_compute_max_level!}}}

!***********************************************************************
!
!  routine ocn_setup_sign_and_index_fields
!
!> \brief   set up sign and index fields
!> \author Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date   September 2011
!> \details 
!>  This routine initializes edgeSignOnCell, edgeSignOnVertex, and 
!>   kiteIndexOnCell.
!
!-----------------------------------------------------------------------
   subroutine ocn_setup_sign_and_index_fields(mesh)!{{{

       type (mesh_type), intent(inout) :: mesh

       integer, dimension(:), pointer :: nEdgesOnCell
       integer, dimension(:,:), pointer :: edgesOnCell, edgesOnVertex, cellsOnVertex, cellsOnEdge, verticesOnCell, verticesOnEdge
       integer, dimension(:,:), pointer :: edgeSignOnCell, edgeSignOnVertex, kiteIndexOnCell

       integer :: nCells, nEdges, nVertices, vertexDegree
       integer :: iCell, iEdge, iVertex, i, j, k

       nCells = mesh % nCells
       nEdges = mesh % nEdges
       nVertices = mesh % nVertices
       vertexDegree = mesh % vertexDegree

       nEdgesOnCell => mesh % nEdgesOnCell % array
       edgesOnCell => mesh % edgeSOnCell % array
       edgesOnVertex => mesh % edgesOnVertex % array
       cellsOnVertex => mesh % cellsOnVertex % array
       cellsOnEdge => mesh % cellsOnEdge % array
       verticesOnCell => mesh % verticesOnCell % array
       verticesOnEdge => mesh % verticesOnEdge % array
       edgeSignOnCell => mesh % edgeSignOnCell % array
       edgeSignOnVertex => mesh % edgeSignOnVertex % array
       kiteIndexOnCell => mesh % kiteIndexOnCell % array

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

   end subroutine ocn_setup_sign_and_index_fields!}}}

!***********************************************************************
!
!  routine ocn_compute_mesh_scaling
!
!> \brief   set up mesh scaling variables
!> \author Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date   September 2011
!> \details 
!>  This routine initializes meshScaling, meshScalingDel2, and
!>   meshScalingDel4
!
!-----------------------------------------------------------------------
   subroutine ocn_compute_mesh_scaling(mesh)!{{{

      use mpas_grid_types
      use mpas_configure

      implicit none

      type (mesh_type), intent(inout) :: mesh

      integer :: iEdge, cell1, cell2
      real (kind=RKIND), dimension(:), pointer :: meshDensity, meshScalingDel2, meshScalingDel4, meshScaling

      meshDensity => mesh % meshDensity % array
      meshScalingDel2 => mesh % meshScalingDel2 % array
      meshScalingDel4 => mesh % meshScalingDel4 % array
      meshScaling     => mesh % meshScaling     % array

      !
      ! Compute the scaling factors to be used in the del2 and del4 dissipation
      !
      meshScalingDel2(:) = 1.0
      meshScalingDel4(:) = 1.0
      meshScaling(:)     = 1.0
      if (config_hmix_ScaleWithMesh) then
         do iEdge=1,mesh%nEdges
            cell1 = mesh % cellsOnEdge % array(1,iEdge)
            cell2 = mesh % cellsOnEdge % array(2,iEdge)
            meshScalingDel2(iEdge) = 1.0 / ( ((meshDensity(cell1) + meshDensity(cell2) )/2.0)/config_maxMeshDensity)**(3.0/4.0)  ! goes as dc**3
            meshScalingDel4(iEdge) = 1.0 / ( ((meshDensity(cell1) + meshDensity(cell2) )/2.0)/config_maxMeshDensity)**(3.0/4.0)  ! goes as dc**3
            meshScaling(iEdge)     = 1.0 / ( ((meshDensity(cell1) + meshDensity(cell2) )/2.0)/config_maxMeshDensity)**(1.0/4.0)
         end do
      end if

   end subroutine ocn_compute_mesh_scaling!}}}

!***********************************************************************
!
!  routine ocn_init_vert_coord
!
!> \brief  initialize vertical coordinate variables
!> \author Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date   September 2011
!> \details 
!>  This routine initializes vertical coordinate variables
!
!-----------------------------------------------------------------------
   subroutine ocn_init_vert_coord(domain)!{{{
   ! Initialize zlevel-type variables and adjust initial conditions for
   ! partial bottom cells.

      use mpas_grid_types
      use mpas_configure

      implicit none

      type (domain_type), intent(inout) :: domain
      type (dm_info) :: dminfo

      integer :: i, iCell, iEdge, iVertex, k, km1, nCells, num_tracers
      type (block_type), pointer :: block

      integer :: iTracer, cell, cell1, cell2
      real (kind=RKIND) :: normalThicknessFluxSum, thicknessSum, hEdge1, zMidPBC

      integer, dimension(:), pointer :: maxLevelCell
      real (kind=RKIND), dimension(:), pointer :: refBottomDepth, &
         refBottomDepthTopOfCell, vertCoordMovementWeights, bottomDepth, refZMid, refLayerThickness
      real (kind=RKIND), dimension(:), allocatable :: minBottomDepth, minBottomDepthMid, zMidZLevel
         
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers
      integer :: nVertLevels
      logical :: consistentSSH

      ! Initialize z-level mesh variables from h, read in from input file.
      block => domain % blocklist
      do while (associated(block))

         layerThickness => block % state % time_levs(1) % state % layerThickness % array
         tracers => block % state % time_levs(1) % state % tracers % array
         refBottomDepth => block % mesh % refBottomDepth % array
         refBottomDepthTopOfCell => block % mesh % refBottomDepthTopOfCell % array
         bottomDepth => block % mesh % bottomDepth % array
         vertCoordMovementWeights => block % mesh % vertCoordMovementWeights % array
         maxLevelCell => block % mesh % maxLevelCell % array

         restingThickness => block % verticalMesh % restingThickness % array
         refZMid => block % verticalMesh % refZMid % array
         refLayerThickness => block % verticalMesh % refLayerThickness % array

         nCells      = block % mesh % nCells
         nVertLevels = block % mesh % nVertLevels
         num_tracers = size(tracers, dim=1)

         ! TopOfCell needed where zero depth for the very top may be referenced.
         refBottomDepthTopOfCell(1) = 0.0
         do k = 1,nVertLevels
            refBottomDepthTopOfCell(k+1) = refBottomDepth(k)
            refLayerThickness(k) = refBottomDepth(k) - refBottomDepthTopOfCell(k)
            refZMid(k) = - refBottomDepthTopOfCell(1) - refLayerThickness(k)/2.0
         end do

         ! Initialization of vertCoordMovementWeights. This determines how SSH perturbations
         ! are distributed throughout the column.
         if (config_vert_coord_movement.eq.'fixed') then

           vertCoordMovementWeights = 0.0
           vertCoordMovementWeights(1) = 1.0

         elseif (config_vert_coord_movement.eq.'uniform_stretching') then

            vertCoordMovementWeights = 1.0

         endif

         ! Initial condition files (ocean.nc, produced by basin) include a realistic
         ! bottomDepth variable and h,T,S variables for full thickness cells.
         ! If running with pbcs, set config_alter_ICs_for_pbc='zlevel_pbcs_on'. Then thin pbc cells
         !    will be changed, and h,T,S will be altered to match the pbcs.
         ! If running without pbcs, set config_alter_ICs_for_pbc='zlevel_pbcs_off'. Then 
         !    bottomDepth will be altered so it is full cells everywhere.
         !    If your input file does not include bottomDepth, the false option will
         !    initialize bottomDepth correctly for a non-pbc run.

         if (.not.config_do_restart.and.config_alter_ICs_for_pbcs) then

            if (config_pbc_alteration_type.eq.'partial_cell') then

               write (stdoutUnit,'(a)') ' Altering bottomDepth to avoid very thin cells.'
               write (stdoutUnit,'(a)') ' Altering layerThickness and tracer initial conditions to conform with partial bottom cells.'

               allocate(minBottomDepth(nVertLevels),minBottomDepthMid(nVertLevels),zMidZLevel(nVertLevels))

               ! min_pbc_fraction restricts pbcs from being too small.
               ! A typical value is 10%, so pbcs must occupy at least 10% of the cell thickness.
               ! If min_pbc_fraction = 0.0, bottomDepth gives the actual depth for that cell.
               ! If min_pbc_fraction = 1.0, bottomDepth reverts to discrete z-level depths, same 
               !    as partial_bottom_cells = .false.

               minBottomDepth(1) = (1.0-config_min_pbc_fraction)*refBottomDepth(1)
               minBottomDepthMid(1) = 0.5*(minBottomDepth(1) + refBottomDepthTopOfCell(1))
               zMidZLevel(1) = - 0.5*(refBottomDepth(1) + refBottomDepthTopOfCell(1))
               do k=2,nVertLevels
                  minBottomDepth(k) = refBottomDepth(k) - (1.0-config_min_pbc_fraction)*(refBottomDepth(k) - refBottomDepth(k-1))
                  minBottomDepthMid(k) = 0.5*(minBottomDepth(k) + refBottomDepthTopOfCell(k))
                  zMidZLevel(k) = - 0.5*(refBottomDepth(k) + refBottomDepthTopOfCell(k))
               enddo

               do iCell=1,nCells

                  ! Change value of maxLevelCell for partial bottom cells
                  k = maxLevelCell(iCell)
                  if (bottomDepth(iCell).lt.minBottomDepthMid(k)) then
                     ! Round up to cell above
                     maxLevelCell(iCell) = maxLevelCell(iCell) - 1
                     bottomDepth(iCell) = refBottomDepth(maxLevelCell(iCell))
                  elseif (bottomDepth(iCell).lt.minBottomDepth(k)) then
                     ! Round down cell to the min_pbc_fraction.
                     bottomDepth(iCell) = minBottomDepth(k)
                  endif
                  ! reset k to new value of maxLevelCell
                  k = maxLevelCell(iCell)

                  ! Alter thickness of bottom level to account for PBC
                  layerThickness(k,iCell) = bottomDepth(iCell) - refBottomDepthTopOfCell(k) 

                  ! Linearly interpolate the initial T&S for new location of bottom cell for PBCs
                  zMidPBC = -0.5*(bottomDepth(iCell) + refBottomDepthTopOfCell(k))
                  km1 = max(k-1,1)
                  do iTracer=1,num_tracers
                     tracers(iTracer,k,iCell) = tracers(iTracer,k,iCell) &
                        + (tracers(iTracer,km1,iCell) - tracers(iTracer,k,iCell)) &
                         /(zMidZLevel(km1)-zMidZLevel(k)+1.0e-16) &
                         *(zMidPBC - zMidZLevel(k))
                  enddo

               enddo  

               deallocate(minBottomDepth,zMidZLevel)

            elseif (config_pbc_alteration_type.eq.'full_cell') then

               do iCell = 1,nCells
                  bottomDepth(iCell) = refBottomDepth(maxLevelCell(iCell))
               enddo

            else

                write (stderrUnit,*) ' Incorrect choice of config_pbc_alteration_type.'
                call mpas_dmpar_abort(dminfo)

            endif

         endif ! .not.config_do_restart

         if (.not.config_do_restart) then

            ! Layer thickness when the ocean is at rest, i.e. without SSH or internal perturbations.
            ! This is applied only from the initial condition
            if (config_set_restingThickness_to_IC) then
                restingThickness = layerThickness
            endif

         endif ! .not.config_do_restart.and.config_alter_ICs_for_pbcs

         if (config_check_ssh_consistency) then
            consistentSSH = .true.
            do iCell = 1,nCells
               ! Check if abs(ssh)>2m.  If so, print warning.
               if (abs(sum(layerThickness(1:maxLevelCell(iCell),iCell))-bottomDepth(iCell))>2.0) then
                  consistentSSH = .false.
                  write (stderrUnit,'(a)') ' Warning: abs(sum(h)-bottomDepth)>2m.  Most likely, initial layerThickness does not match bottomDepth.'
                  write (stderrUnit,*) ' iCell, K=maxLevelCell(iCell), bottomDepth(iCell),sum(h),bottomDepth: ', &
                                iCell, maxLevelCell(iCell), bottomDepth(iCell),sum(layerThickness(1:maxLevelCell(iCell),iCell)),bottomDepth(iCell), &
                                layerThickness(maxLevelCell(iCell),iCell)
               endif
            enddo

            if (.not. consistentSSH) then
               write(stderrUnit,*) 'Warning: SSH is not consistent. Most likely, initial layerThickness does not match bottomDepth.'
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
                  write (stderrUnit,'(a,10f10.2)') ' refBottomDepth(maxLevelCell(iCell)), refBottomDepthTopOfCell(maxLevelCell(iCell)): ', &
                                refBottomDepth(maxLevelCell(iCell)), refBottomDepthTopOfCell(maxLevelCell(iCell))
                  call mpas_dmpar_abort(dminfo)
               endif

            enddo
         endif

      block => block % next
      end do

   end subroutine ocn_init_vert_coord!}}}

end module ocn_init

! vim: foldmethod=marker
