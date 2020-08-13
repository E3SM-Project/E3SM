! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_hmix_del4
!
!> \brief Ocean horizontal mixing - biharmonic parameterization
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This module contains routines and variables for computing
!>  horizontal mixing tendencies using a biharmonic formulation.
!
!-----------------------------------------------------------------------

module ocn_vel_hmix_del4

   use mpas_timer
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_threading
   use mpas_vector_operations
   use mpas_matrix_operations
   use ocn_constants
   use ocn_config

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

   public :: ocn_vel_hmix_del4_tend, &
             ocn_vel_hmix_del4_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: hmixDel4On       !< local flag to determine whether del4 chosen

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_hmix_del4_tend
!
!> \brief   Computes tendency term for biharmonic horizontal momentum mixing
!> \author  Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date    September 2011
!> \details
!>  This routine computes the horizontal mixing tendency for momentum
!>  based on a biharmonic form for the mixing.  This mixing tendency
!>  takes the form  \f$-\nu_4 \nabla^4 u\f$
!>  but is computed as
!>  \f$\nabla^2 u = \nabla divergence + k \times \nabla relativeVorticity\f$
!>  applied recursively.
!>  This formulation is only valid for constant \f$\nu_4\f$ .
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_hmix_del4_tend(meshPool, divergence, relativeVorticity, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         divergence      !< Input: velocity divergence

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         relativeVorticity       !< Input: relative vorticity

      type (mpas_pool_type), intent(in) :: &
         meshPool           !< Input: mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend       !< Input/Output: velocity tendency

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, cell1, cell2, vertex1, vertex2, k, i
      integer :: iCell, iVertex, nEdges, nCells, nVertices
      integer, pointer :: nVertLevels, vertexDegree
      integer, dimension(:), pointer :: nEdgesArray, nCellsArray, nVerticesArray

      integer, dimension(:), pointer :: maxLevelEdgeTop, maxLevelVertexTop, &
            maxLevelCell, nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, verticesOnEdge, edgeMask, edgesOnVertex, edgesOnCell, edgeSignOnVertex, &
                                          edgeSignOnCell


      real (kind=RKIND) :: u_diffusion, invAreaCell1, invAreaCell2, invAreaTri1, &
            invAreaTri2, invDcEdge, invDvEdge, r_tmp
      real (kind=RKIND), dimension(:), pointer :: dcEdge, dvEdge, areaTriangle, &
            meshScalingDel4, areaCell

      ! Scratch Arrays
      real (kind=RKIND), dimension(:,:), allocatable :: delsq_divergence
      real (kind=RKIND), dimension(:,:), allocatable :: delsq_relativeVorticity
      real (kind=RKIND), dimension(:,:), allocatable :: delsq_u

      err = 0

      if(.not.hmixDel4On) return

      call mpas_timer_start("vel del4")

      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
      call mpas_pool_get_dimension(meshPool, 'nVerticesArray', nVerticesArray)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)

      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'maxLevelVertexTop', maxLevelVertexTop)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'meshScalingDel4', meshScalingDel4)
      call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnVertex', edgesOnVertex)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgeSignOnVertex', edgeSignOnVertex)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)

      nEdges = nEdgesArray(size(nEdgesArray))
      nCells = nCellsArray(size(nCellsArray))
      nVertices = nVerticesArray(size(nVerticesArray))
      allocate(delsq_u(nVertLevels, nEdges + 1), &
               delsq_divergence(nVertLevels, nCells + 1), &
               delsq_relativeVorticity(nVertLevels, nVertices + 1))

      nEdges = nEdgesArray( 3 )

      !Compute delsq_u
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(cell1, cell2, vertex1, vertex2, invDcEdge, invDvEdge, k)
      do iEdge = 1, nEdges
         delsq_u(:, iEdge) = 0.0_RKIND
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         vertex1 = verticesOnEdge(1,iEdge)
         vertex2 = verticesOnEdge(2,iEdge)

         invDcEdge = 1.0_RKIND / dcEdge(iEdge)
         invDvEdge = 1.0_RKIND / max(dvEdge(iEdge), 0.25_RKIND*dcEdge(iEdge))

         do k=1,maxLevelEdgeTop(iEdge)
            ! Compute \nabla^2 u = \nabla divergence + k \times \nabla relativeVorticity
            delsq_u(k, iEdge) = ( divergence(k,cell2)  - divergence(k,cell1) ) * invDcEdge  &
                               -( relativeVorticity(k,vertex2) - relativeVorticity(k,vertex1)) * invDvEdge
         end do
      end do
      !$omp end do
      !$omp end parallel

      nVertices = nVerticesArray( 2 )

      ! Compute delsq_relativeVorticity
      !$omp parallel
      !$omp do schedule(runtime) private(invAreaTri1, i, iEdge, k)
      do iVertex = 1, nVertices
         delsq_relativeVorticity(:, iVertex) = 0.0_RKIND
         invAreaTri1 = 1.0_RKIND / areaTriangle(iVertex)
         do i = 1, vertexDegree
            iEdge = edgesOnVertex(i, iVertex)
            do k = 1, maxLevelVertexTop(iVertex)
               delsq_relativeVorticity(k, iVertex) = delsq_relativeVorticity(k, iVertex) + edgeSignOnVertex(i, iVertex) &
                                                   * dcEdge(iEdge) * delsq_u(k, iEdge) * invAreaTri1
            end do
         end do
      end do
      !$omp end do
      !$omp end parallel

      nCells = nCellsArray( 2 )

      ! Compute delsq_divergence
      !$omp parallel
      !$omp do schedule(runtime) private(invAreaCell1, i, iEdge, k)
      do iCell = 1, nCells
         delsq_divergence(:, iCell) = 0.0_RKIND
         invAreaCell1 = 1.0_RKIND / areaCell(iCell)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            do k = 1, maxLevelCell(iCell)
               delsq_divergence(k, iCell) = delsq_divergence(k, iCell) - edgeSignOnCell(i, iCell) * dvEdge(iEdge) &
                                          * delsq_u(k, iEdge) * invAreaCell1
            end do
         end do
      end do
      !$omp end do
      !$omp end parallel

      nEdges = nEdgesArray( 1 )

      ! Compute - \kappa \nabla^4 u
      ! as  \nabla div(\nabla^2 u) + k \times \nabla ( k \cross curl(\nabla^2 u) )
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(cell1, cell2, vertex1, vertex2, invDcEdge, invDvEdge, r_tmp, k, u_diffusion)
      do iEdge = 1, nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         vertex1 = verticesOnEdge(1,iEdge)
         vertex2 = verticesOnEdge(2,iEdge)

         invDcEdge = 1.0_RKIND / dcEdge(iEdge)
         invDvEdge = 1.0_RKIND / dvEdge(iEdge)
         r_tmp = config_mom_del4 * meshScalingDel4(iEdge)

         do k=1,maxLevelEdgeTop(iEdge)
            u_diffusion = config_mom_del4_div_factor*(delsq_divergence(k,cell2) - delsq_divergence(k,cell1)) * invDcEdge  &
                        - (delsq_relativeVorticity(k,vertex2) - delsq_relativeVorticity(k,vertex1) ) * invDvEdge

            tend(k,iEdge) = tend(k,iEdge) - edgeMask(k, iEdge) * u_diffusion * r_tmp
         end do
      end do
      !$omp end do
      !$omp end parallel

      deallocate(delsq_u, &
                 delsq_divergence, &
                 delsq_relativeVorticity)

      call mpas_timer_stop("vel del4")

   !--------------------------------------------------------------------

   end subroutine ocn_vel_hmix_del4_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_hmix_del4_init
!
!> \brief   Initializes ocean momentum biharmonic horizontal mixing
!> \author  Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date    September 2011
!> \details
!>  This routine initializes a variety of quantities related to
!>  biharmonic horizontal tracer mixing in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_hmix_del4_init(err)!{{{

   integer, intent(out) :: err !< Output: error flag

   !--------------------------------------------------------------------
   !
   ! set some local module variables based on input config choices
   !
   !--------------------------------------------------------------------

   err = 0

   hmixDel4On = .false.

   if ( config_mom_del4 > 0.0_RKIND ) then
      hmixDel4On = .true.
   endif

   if(.not.config_use_mom_del4) hmixDel4On = .false.

   !--------------------------------------------------------------------

   end subroutine ocn_vel_hmix_del4_init!}}}

!***********************************************************************

end module ocn_vel_hmix_del4

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
