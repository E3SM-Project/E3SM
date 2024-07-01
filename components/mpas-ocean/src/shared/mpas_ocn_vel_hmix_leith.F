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

   public :: ocn_vel_hmix_leith_tend, &
             ocn_vel_hmix_leith_init

   !-------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: &
      hmixLeithOff  !< on/off switch to determine whether leith chosen

   real (kind=RKIND) :: &
      leithParam,       &!< Leith parameter
      dxLeith,          &!< Leith length scale
      viscMaxLeith,     &!< maximum viscosity
      sqrt3fact          !< sqrt(3)

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

   subroutine ocn_vel_hmix_leith_tend(div, relVort, tend, err)!{{{

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
      if (hmixLeithOff) return
      call mpas_timer_start("vel leith")

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

         visc2tmp = (leithParam*dxLeith*meshScalingDel2(iEdge)/pi)**3

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
                   *abs(relVort(k,vertex2) - relVort(k,vertex1))* &
                    dcEdgeInv*sqrt3fact
            visc2 = min(visc2, viscMaxLeith)

            tend(k,iEdge) = tend(k,iEdge) + &
                            edgeMask(k,iEdge)*visc2*uDiff

         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      call mpas_timer_stop("vel leith")

   !--------------------------------------------------------------------

   end subroutine ocn_vel_hmix_leith_tend!}}}

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
      hmixLeithOff = .true.
      leithParam   = 0.0_RKIND
      dxLeith      = 0.0_RKIND
      viscMaxLeith = 0.0_RKIND
      sqrt3fact    = sqrt(3.0_RKIND)

      !*** reset values based on input configuration

      if (config_use_leith_del2) then
         hmixLeithOff = .false.
         leithParam   = config_leith_parameter
         dxLeith      = config_leith_dx
         viscMaxLeith = config_leith_visc2_max
      endif

      !-----------------------------------------------------------------

   end subroutine ocn_vel_hmix_leith_init!}}}

!***********************************************************************

end module ocn_vel_hmix_leith

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
