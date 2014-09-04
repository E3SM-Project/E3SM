! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_diagnostics_routines
!
!> \brief MPAS ocean diagnostics driver
!> \author Mark Petersen
!> \date   November 2013
!> \details
!>  This module contains the routines for computing individual diagnostic variables
!
!-----------------------------------------------------------------------

module ocn_diagnostics_routines

   use mpas_grid_types
   use mpas_constants
   use mpas_timer

   implicit none
   private
   save

   type (timer_node), pointer :: diagEOSTimer

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

   public :: ocn_relativeVorticity_circulation

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_relativeVorticity_circulation
!
!> \brief   Computes relative vorticity and circulation
!> \author  Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date    November 2013
!> \details 
!>  Computes relative vorticity and circulation
!
!-----------------------------------------------------------------------

   subroutine ocn_relativeVorticity_circulation(relativeVorticity, circulation, meshPool, normalVelocity, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         relativeVorticity 

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         circulation       

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iVertex, iEdge, i, k
      integer, pointer :: nEdges, nVertices, vertexDegree
      integer, dimension(:), pointer :: maxLevelVertexBot
      integer, dimension(:,:), pointer :: edgesOnVertex, edgeSignOnVertex

      real (kind=RKIND) :: invAreaTri1, r_tmp
      real (kind=RKIND), dimension(:), pointer :: &
              dcEdge, areaTriangle

      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)

      call mpas_pool_get_array(meshPool, 'maxLevelVertexBot', maxLevelVertexBot)
      call mpas_pool_get_array(meshPool, 'edgesOnVertex', edgesOnVertex)
      call mpas_pool_get_array(meshPool, 'edgeSignOnVertex', edgeSignOnVertex)

      call mpas_pool_get_array(meshPool, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)

      err = 0

      circulation(:,:) = 0.0
      relativeVorticity(:,:) = 0.0
      do iVertex = 1, nVertices
         invAreaTri1 = 1.0 / areaTriangle(iVertex)
         do i = 1, vertexDegree
            iEdge = edgesOnVertex(i, iVertex)
            do k = 1, maxLevelVertexBot(iVertex)
              r_tmp = dcEdge(iEdge) * normalVelocity(k, iEdge)

              circulation(k, iVertex) = circulation(k, iVertex) + edgeSignOnVertex(i, iVertex) * r_tmp 
              relativeVorticity(k, iVertex) = relativeVorticity(k, iVertex) + edgeSignOnVertex(i, iVertex) * r_tmp * invAreaTri1
            end do
         end do
      end do


   !--------------------------------------------------------------------

   end subroutine ocn_relativeVorticity_circulation!}}}

end module ocn_diagnostics_routines

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
