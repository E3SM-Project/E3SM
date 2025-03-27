! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_hmix_leith
!
!> \brief Ocean horizontal mixing - Leith parameterization
!> \author Mark Petersen
!> \date   22 October 2012
!> \details
!>  This module contains routines for computing horizontal mixing
!>  tendencies using the Leith parameterization.
!
!-----------------------------------------------------------------------

module ocn_vel_hmix_leith

   use mpas_timer
   use mpas_constants
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

   public :: ocn_vel_hmix_leith_del2_tend, &
             ocn_vel_hmix_leith_del4_tend, &
             ocn_vel_hmix_leith_init

   !-------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: &
      hmixLeith2Off, & !< on/off switch to determine whether leith del2 chosen
      hmixLeith4Off    !< on/off switch to determine whether leith del4 chosen

   real (kind=RKIND) :: &
      leithParam2,     &!< Leith del2 parameter
      leithParam4,     &!< Leith del4 parameter
      leithDivScaler2, &!< Leith del4 divergence scale factor
      leithDivScaler4, &!< Leith del4 divergence scale factor
      dxLeith,         &!< Leith length scale
      visc2MaxLeith,   &!< maximum del2 viscosity
      visc4MaxLeith,   &!< maximum del4 viscosity
      sqrt3fact         !< sqrt(3)

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_hmix_leith_tend
!
!> \brief  Computes horizontal momentum mixing with Leith formulation
!> \author Mark Petersen, Todd Ringler
!> \date   22 October 2012
!> \details
!> This routine computes the horizontal mixing tendency for momentum
!> based on the Leith closure.  The Leith closure is the
!> enstrophy-cascade analogy to the Smagorinsky (1963) energy-cascade
!> closure, i.e. Leith (1996) assumes an inertial range of enstrophy
!> flux moving toward the mesh scale. The assumption of an enstrophy
!> cascade and dimensional analysis produces right-hand-side
!> dissipation, $\bf{D}$, of velocity of the form
!> $ {\bf D} = \nabla \cdot \left( \nu_\ast \nabla {\bf u} \right)
!>    = \nabla \cdot \left( \gamma \left| \nabla \omega  \right|
!>      \left( \Delta x \right)^3 \nabla \bf{u} \right)
!> where $\omega$ is the relative vorticity and $\gamma$ is a
!> non-dimensional, $O(1)$ parameter. We set $\gamma=1$.
!
!-----------------------------------------------------------------------

!***********************************************************************
!  routine ocn_vel_hmix_leith_del2_tend
!***********************************************************************

   subroutine ocn_vel_hmix_leith_del2_tend(div, relVort, tend, err)!{{{

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
         tend             !< [inout] accumulated velocity tendency

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< [out] error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer ::           &
         iEdge, k,         &! loop indices for edge, vertical loops
         cell1, cell2,     &! neighbor cell addresses across edge
         vertex1, vertex2   ! neighbor vertex addresses along edge

      real (kind=RKIND) :: &!
         uDiff,            &! velocity diffusion operator
         dcEdgeInv,        &! 1/dcEdge
         dvEdgeInv,        &! 1/dvEdge
         visc2tmp,         &! common factor for visc2
         visc2              ! scaled viscosity coeff

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** initialize return code and exit if Leith not chosen
      !*** start timer if chosen

      err = 0
      if (hmixLeith2Off) return
      call mpas_timer_start("vel leith del2")

#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(cellsOnEdge, verticesOnEdge, minLevelEdgeBot, maxLevelEdgeTop, &
      !$acc            dcEdge, dvEdge, meshScalingDel2, div, relVort, &
      !$acc            tend, edgeMask) &
      !$acc    private(k, cell1, cell2, vertex1, vertex2, dcEdgeInv, &
      !$acc            dvEdgeInv, uDiff, visc2, visc2tmp)
#else
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp    private(k, cell1, cell2, vertex1, vertex2, dcEdgeInv, &
      !$omp            dvEdgeInv, uDiff, visc2, visc2tmp)
#endif
      do iEdge = 1, nEdgesOwned
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         vertex1 = verticesOnEdge(1,iEdge)
         vertex2 = verticesOnEdge(2,iEdge)

         dcEdgeInv = 1.0_RKIND / dcEdge(iEdge)
         dvEdgeInv = 1.0_RKIND / dvEdge(iEdge)

         visc2tmp = (leithParam2*dxLeith*meshScalingDel2(iEdge)/pi)**3

         do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)

            ! Here -( relativeVorticity(k,vertex2) -
            !         relativeVorticity(k,vertex1) ) / dvEdge(iEdge)
            ! is - \nabla relativeVorticity pointing from vertex 2 to
            ! vertex 1, or equivalently
            ! + k \times \nabla relativeVorticity pointing from cell1
            ! to cell2.

            uDiff = (div(k,cell2)  - div(k,cell1))*dcEdgeInv &
                   -(relVort(k,vertex2) - relVort(k,vertex1))*dvEdgeInv

            ! Here the first line is (\delta x)^3
            ! the second line is |\nabla \omega|
            ! and u_diffusion is \nabla^2 u (see formula for $\bf{D}$ above).
            visc2 = visc2tmp &
                   *sqrt( (abs(relVort(k,vertex2) - relVort(k,vertex1))* &
                           dcEdgeInv*sqrt3fact)**2.0 + &
                           leithDivScaler2 * &
                           abs((div(k,cell2)  - div(k,cell1))*dcEdgeInv)**2.0 )
            !visc2 = min(visc2, viscMaxLeith)
            visc2 = min(visc2, visc2MaxLeith)

            tend(k,iEdge) = tend(k,iEdge) + &
                            edgeMask(k,iEdge)*visc2*uDiff

         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      call mpas_timer_stop("vel leith del2")

   !--------------------------------------------------------------------
   end subroutine ocn_vel_hmix_leith_del2_tend!}}}


!***********************************************************************
!  routine ocn_vel_hmix_leith_del4_tend
!***********************************************************************

   !--------------------------------------------------------------------
   subroutine ocn_vel_hmix_leith_del4_tend(div, relVort, tend, err)!{{{

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
         tend             !< [inout] accumulated velocity tendency

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
         del2relVortEdge,  &! del2relVort on Edges
         del2divEdge,      &! del2div on Edges
         visc4tmp,         &! common factor for visc4
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

      !*** initialize return code and exit if Leith not chosen
      !*** start timer if chosen

      err = 0
      if (hmixLeith4Off) return
      call mpas_timer_start("vel leith del4")


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
      !$acc private(k, cell1, cell2, vertex1, vertex2, del2relVortEdge, del2divEdge, &
      !$acc         dcEdgeInv, dvEdgeInv, visc4, visc4tmp uDiff, kmin, kmax)
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
         visc4tmp = (leithParam4*dxLeith*meshScalingDel4(iEdge)/pi)**6

         do k = kmin, kmax
            uDiff = (del2div(k,cell2) - del2div(k,cell1))* &
                                                dcEdgeInv  &
                  - (del2relVort(k,vertex2) - del2relVort(k,vertex1))* &
                                                          dvEdgeInv

            del2relVortEdge = 0.5_RKIND * &
                             (del2relVort(k,vertex1)+del2relVort(k,vertex2))

            del2divEdge = 0.5_RKIND * &
                         (del2div(k,cell1)+del2div(k,cell2))

            visc4 = visc4tmp &
                   *sqrt((abs(del2relVortEdge)**2.0 + &
                          leithDivScaler4*abs(del2divEdge)**2.0))

            visc4 = min(visc4, visc4MaxLeith)

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

      call mpas_timer_stop("vel leith del4")

   !--------------------------------------------------------------------

   end subroutine ocn_vel_hmix_leith_del4_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_hmix_leith_init
!
!> \brief  Initializes ocean momentum horiz mixing in Leith formulation
!> \author Mark Petersen
!> \date   22 October 2012
!> \details
!>  This routine initializes a variety of quantities related to
!>  Leith parameterization for horizontal momentum mixing in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_hmix_leith_init(err)!{{{

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< [out] error flag

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** initialize return flag and set default values

      err = 0
      hmixLeith2Off   = .true.
      hmixLeith4Off   = .true.
      leithParam2     = 0.0_RKIND
      leithParam4     = 0.0_RKIND
      leithDivScaler2 = 0.0_RKIND
      leithDivScaler4 = 0.0_RKIND
      dxLeith         = 0.0_RKIND
      visc2MaxLeith   = 0.0_RKIND
      visc4MaxLeith   = 0.0_RKIND
      sqrt3fact       = sqrt(3.0_RKIND)

      !*** reset values based on input configuration

      if (config_use_leith_del2) then
         hmixLeith2Off   = .false.
         leithParam2     = config_leith_del2_parameter
         leithDivScaler2 = config_leith_del2_div_scaler ** 2.0
         dxLeith         = config_leith_dx
         visc2MaxLeith   = config_leith_visc2_max
      endif

      if (config_use_leith_del4) then
         hmixLeith4Off   = .false.
         leithParam4     = config_leith_del4_parameter
         leithDivScaler4 = config_leith_del4_div_scaler ** 2.0
         dxLeith         = config_leith_dx
         visc4MaxLeith   = config_leith_visc4_max
      endif

      !-----------------------------------------------------------------

   end subroutine ocn_vel_hmix_leith_init!}}}

!***********************************************************************

end module ocn_vel_hmix_leith

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
