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
!> \brief Main driver for MPAS ocean analysis core in post-processing mode
!> \author Mark Petersen
!> \date   November 2013
!> \details
!>  This module contains the drivers for the MPAS ocean analysis core in 
!>  post-processing mode.  None of these routines are used in run-time mode
!>  by the ocean core.
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
      ocn_setup_sign_and_index_fields
!      ocn_init_vert_coord, &
!      ocn_compute_max_level, &
!      ocn_compute_mesh_scaling, &
      

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

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


end module ocn_init

! vim: foldmethod=marker
