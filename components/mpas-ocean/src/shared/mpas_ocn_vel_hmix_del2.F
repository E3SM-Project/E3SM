! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_hmix_del2
!
!> \brief Ocean horizontal mixing - Laplacian parameterization
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This module contains routines for computing horizontal mixing
!>  tendencies using a Laplacian formulation.
!
!-----------------------------------------------------------------------

module ocn_vel_hmix_del2

   use mpas_timer
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

   public :: ocn_vel_hmix_del2_tend, &
             ocn_vel_hmix_del2_init

   !-------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical ::  &
      hmixDel2Off  !< on/off flag to determine whether del2 chosen
         
   real (kind=RKIND) :: &
      viscDel2     !< viscosity coefficient for del2 horz mixing

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_hmix_del2_tend
!
!> \brief   Computes tendency for Laplacian horizontal momentum mixing
!> \author  Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date    22 August 2011
!> \details
!>  This routine computes the horizontal mixing tendency for momentum
!>  based on a Laplacian form for the mixing, \f$\nu_2 \nabla^2 u\f$
!>  This tendency takes the form
!>  \f$\nu( \nabla divergence + k \times \nabla relativeVorticity )\f$,
!>  where \f$\nu\f$ is a viscosity and \f$k\f$ is the vertical unit 
!>  vector. This form is strictly only valid for constant \f$\nu\f$ .
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_hmix_del2_tend(div, relVort, tend, err)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         div,           &!< [in] velocity divergence
         relVort         !< [in] relative vorticity

      !-----------------------------------------------------------------
      ! input /output variables
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

      integer ::           &
         iEdge, k,         &! loop indices for edge, vertical loops
         cell1, cell2,     &! neighbor cell addresses across edge
         vertex1, vertex2   ! neighbor vertex addresses along edge

      real (kind=RKIND) :: &
         uDiff,            &! velocity diffusion
         dcEdgeInv,        &! 1/dcEdge
         dvEdgeInv,        &! 1/dvEdge
         visc2              ! scaled viscosity coefficient

      !-----------------------------------------------------------------
      !
      ! exit if this mixing is not selected
      !
      !-----------------------------------------------------------------

      err = 0

      if (hmixDel2Off) return

      call mpas_timer_start("vel del2")

#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(cellsOnEdge, maxLevelEdgeTop, minLevelEdgeBot, &
      !$acc            verticesOnEdge, edgeMask, dcEdge, dvEdge, &
      !$acc            meshScalingDel2, div, relVort, tend) &
      !$acc    private(k, cell1, cell2, uDiff, vertex1, vertex2, &
      !$acc            dcEdgeInv, dvEdgeInv, visc2)
#else
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp    private(k, cell1, cell2, uDiff, vertex1, vertex2, &
      !$omp            dcEdgeInv, dvEdgeInv, visc2)
#endif
      do iEdge = 1, nEdgesOwned
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         vertex1 = verticesOnEdge(1,iEdge)
         vertex2 = verticesOnEdge(2,iEdge)

         dcEdgeInv = 1.0_RKIND / dcEdge(iEdge)
         dvEdgeInv = 1.0_RKIND / dvEdge(iEdge)

         visc2 =  viscDel2*meshScalingDel2(iEdge)

         do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)

            ! Here -( relativeVorticity(k,vertex2) - 
            !         relativeVorticity(k,vertex1) ) / dvEdge(iEdge)
            ! is - \nabla relativeVorticity pointing from vertex 2 
            ! to vertex 1, or equivalently
            ! + k \times \nabla relativeVorticity pointing from cell1 
            ! to cell2.

            uDiff = (div(k,cell2) - div(k,cell1))*dcEdgeInv &
                   -(relVort(k,vertex2)-relVort(k,vertex1))*dvEdgeInv

            tend(k,iEdge) = tend(k,iEdge) + &
                            edgeMask(k,iEdge)*visc2*uDiff

         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      call mpas_timer_stop("vel del2")

   !--------------------------------------------------------------------

   end subroutine ocn_vel_hmix_del2_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_hmix_del2_init
!
!> \brief   Initializes ocean momentum Laplacian horizontal mixing
!> \author  Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date    September 2011
!> \details
!>  This routine initializes a variety of quantities related to
!>  Laplacian horizontal momentum mixing in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_hmix_del2_init(err)!{{{

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< [out] error flag

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** Initialize return code and default variables
      err = 0
      hmixDel2Off = .true.
      viscDel2    = 0.0_RKIND

      !*** reset values based on input configuration

      if (config_use_mom_del2) then 
         hmixDel2Off = .false.
         viscDel2    = config_mom_del2

         if (viscDel2 <= 0.0_RKIND) then
            call mpas_log_write( &
                'vel_hmix_del2_init: viscosity must be > 0', &
                 MPAS_LOG_ERR)
            err = -1
         endif
      endif

      !-----------------------------------------------------------------

   end subroutine ocn_vel_hmix_del2_init!}}}

!***********************************************************************

end module ocn_vel_hmix_del2

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
