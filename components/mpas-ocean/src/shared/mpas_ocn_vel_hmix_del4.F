! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
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
   use mpas_log
   use ocn_constants
   use ocn_config
   use ocn_mesh

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

   logical :: &
      hmixDel4Off      !< on/off flag to determine whether del4 chosen

   real (kind=RKIND) :: &
      viscDel4,         &!< biharmonic viscosity coefficient
      divFactor          !< factor for divergence term

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_hmix_del4_tend
!
!> \brief   Computes tendency for biharmonic horizontal momentum mixing
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

   subroutine ocn_vel_hmix_del4_tend(div, relVort, tend, err)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         div,           &!< [in] velocity divergence
         relVort         !< [in] relative vorticity

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend            !< [inout] accumulated velocity tendency

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< [out] error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer ::                      &
         iEdge, iCell, iVertex, k, i, &! loop counters 
         cell1, cell2,     &! neighbor cell addresses across edge
         vertex1, vertex2, &! neighbor vertex addresses along edge
         nEdges, nCells, nVertices ! num edges,cells,vertices incl halos

      real (kind=RKIND) :: &
         uDiff,            &! diffusion operator temporary
         areaCellInv,      &! 1/area of cell
         areaTriInv,       &! 1/area of triangle
         dcEdgeInv,        &! 1/dcEdge
         dvEdgeInv,        &! 1/dvEdge
         visc4              ! scaled biharmonic viscosity coeff

      ! Scratch Arrays
      real (kind=RKIND), dimension(:,:), allocatable :: &
         del2div,     &!
         del2relVort, &!
         del2u

      integer :: kmin, kmax

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** initialize error code and return if not selected
      !*** otherwise start timer

      err = 0
      if (hmixDel4Off) return
      call mpas_timer_start("vel del4")

      !*** allocate temporaries

      allocate(del2u(nVertLevels, nEdgesAll + 1), &
               del2div(nVertLevels, nCellsAll + 1), &
               del2relVort(nVertLevels, nVerticesAll + 1))
      !$acc enter data create(del2u, del2div, del2relVort)

      nEdges = nEdgesHalo(2)

      !Compute Laplacian of velocity
#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(cellsOnEdge, verticesOnEdge, minLevelEdgeBot, maxLevelEdgeTop, &
      !$acc            dvEdge, dcEdge, del2u, div, relVort) &
      !$acc    private(k, cell1, cell2, vertex1, vertex2, &
      !$acc            dcEdgeInv, dvEdgeInv, kmin, kmax)
#else
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp    private(k, cell1, cell2, vertex1, vertex2, &
      !$omp            dcEdgeInv, dvEdgeInv, kmin, kmax)
#endif
      do iEdge = 1, nEdges
         del2u(:, iEdge) = 0.0_RKIND
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         kmin = minLevelEdgeBot(iEdge)
         kmax = maxLevelEdgeTop(iEdge)

         vertex1 = verticesOnEdge(1,iEdge)
         vertex2 = verticesOnEdge(2,iEdge)

         dcEdgeInv = 1.0_RKIND / dcEdge(iEdge)
         dvEdgeInv = 1.0_RKIND / max(dvEdge(iEdge), 0.25_RKIND*dcEdge(iEdge))

         do k=kmin, kmax
            ! Compute \nabla^2 u = \nabla divergence + k \times \nabla relativeVorticity
            del2u(k,iEdge) = (div(k,cell2) - div(k,cell1))*dcEdgeInv &
                            -(relVort(k,vertex2) - relVort(k,vertex1))*&
                                                           dvEdgeInv
         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      nVertices = nVerticesHalo(1)

      ! Compute del2 of relative vorticity
#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(edgesOnVertex, minLevelVertexBot, maxLevelVertexTop, dcEdge, &
      !$acc            areaTriangle, edgeSignOnVertex, del2u, &
      !$acc            del2relVort) &
      !$acc    private(i, k, iEdge, areaTriInv, kmin, kmax)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(i, k, iEdge, areaTriInv, kmin, kmax)
#endif
      do iVertex = 1, nVertices
         kmin = minLevelVertexBot(iVertex)
         kmax = maxLevelVertexTop(iVertex)
         del2relVort(:, iVertex) = 0.0_RKIND
         areaTriInv = 1.0_RKIND/areaTriangle(iVertex)
         do i = 1, vertexDegree
            iEdge = edgesOnVertex(i, iVertex)
            do k = kmin, kmax
               del2relVort(k,iVertex) = del2relVort(k,iVertex) &
                                      + edgeSignOnVertex(i,iVertex) &
                                       *dcEdge(iEdge)*del2u(k,iEdge)*areaTriInv
            end do
         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      nCells = nCellsHalo(1)

      ! Compute del2 of divergence
#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(nEdgesOnCell, edgesOnCell, minLevelCell, maxLevelCell, dvEdge,&
      !$acc            edgeSignOnCell, areaCell, del2u, del2div) &
      !$acc    private(i, k, iEdge, areaCellInv, kmin, kmax)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(i, k, iEdge, areaCellInv, kmin, kmax)
#endif
      do iCell = 1, nCells
         kmin = minLevelCell(iCell)
         kmax = maxLevelCell(iCell)
         del2div(:, iCell) = 0.0_RKIND
         areaCellInv = 1.0_RKIND / areaCell(iCell)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            do k = kmin, kmax
               del2div(k,iCell) = del2div(k,iCell) - &
                                  edgeSignOnCell(i,iCell)*dvEdge(iEdge) &
                                 *del2u(k,iEdge)*areaCellInv
            end do
         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      ! Compute final tendency \kappa \nabla^4 u
      ! as  \nabla div(\nabla^2 u) + k \times 
      !     \nabla ( k \cross curl(\nabla^2 u) )

#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(cellsOnEdge, verticesOnEdge, minLevelEdgeBot, maxLevelEdgeTop, &
      !$acc            dcEdge, dvEdge, meshScalingDel4, edgeMask, &
      !$acc            del2div, del2relVort, tend) &
      !$acc private(k, cell1, cell2, vertex1, vertex2, &
      !$acc         dcEdgeInv, dvEdgeInv, visc4, uDiff, kmin, kmax)
#else
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(k, cell1, cell2, vertex1, vertex2, &
      !$omp         dcEdgeInv, dvEdgeInv, visc4, uDiff, kmin, kmax)
#endif
      do iEdge = 1, nEdgesOwned
         kmin = minLevelEdgeBot(iEdge)
         kmax = maxLevelEdgeTop(iEdge)

         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         vertex1 = verticesOnEdge(1,iEdge)
         vertex2 = verticesOnEdge(2,iEdge)

         dcEdgeInv = 1.0_RKIND / dcEdge(iEdge)
         dvEdgeInv = 1.0_RKIND / dvEdge(iEdge)
         visc4 = viscDel4*meshScalingDel4(iEdge)

         do k = kmin, kmax
            uDiff = divFactor*(del2div(k,cell2) - del2div(k,cell1))* &
                                                          dcEdgeInv  &
                  - (del2relVort(k,vertex2) - del2relVort(k,vertex1))* &
                                                          dvEdgeInv

            tend(k,iEdge) = tend(k,iEdge) - &
                            edgeMask(k,iEdge)*uDiff*visc4
         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      !$acc exit data delete(del2u, del2div, del2relVort)
      deallocate(del2u, &
                 del2div, &
                 del2relVort)

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

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< [out] error flag

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** initialize return error code and set module variable defaults

      err = 0
      hmixDel4Off = .true.
      viscDel4    = 0.0_RKIND
      divFactor   = 1.0_RKIND

      !*** reset module variables based on input configuration

      if (config_use_mom_del4) then
         hmixDel4Off = .false.
         viscDel4    = config_mom_del4
         divFactor   = config_mom_del4_div_factor
         if ( config_mom_del4 <= 0.0_RKIND ) then
            call mpas_log_write( &
               'vel_hmix_del4_init: config_mom_del4 must be > 0 ', &
               MPAS_LOG_ERR)
            err = -1
         endif
      endif

   !--------------------------------------------------------------------

   end subroutine ocn_vel_hmix_del4_init!}}}

!***********************************************************************

end module ocn_vel_hmix_del4

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
