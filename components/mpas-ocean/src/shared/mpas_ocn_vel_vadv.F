! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_vadv
!
!> \brief MPAS ocean vertical advection
!> \author Mark Petersen
!> \date   September 2011
!> \details
!>  This module contains the routine for computing
!>  tendencies for vertical advection.
!>
!
!-----------------------------------------------------------------------

module ocn_vel_vadv

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

   public :: ocn_vel_vadv_tend, &
             ocn_vel_vadv_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: &
      velVadvOff  ! on/off switch for vertical advection


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_vadv_tend
!
!> \brief   Computes tendency term for vertical advection
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  This routine computes the vertical advection tendency for momentum
!>  based on current state.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_vadv_tend(normalVelocity, layerThickEdgeFlux, &
                                vertAleTransportTop, tend, err)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity,     &!< [in] Horizontal velocity
         layerThickEdgeFlux,     &!< [in] Layer thickness at edge
         vertAleTransportTop  !< [in] Vertical velocity on top layer

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend          !< [inout] accumulated velocity tendency

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< [out] error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer ::       &
         iEdge, k,     &! loop indices for edge, vertical loops
         kmin, kmax,   &! shallowest and deepest active layer on edge
         cell1, cell2   ! neighbor cell indices across edge

      real (kind=RKIND) :: &
         wAvg           ! ALE transport velocity across top edge

      real (kind=RKIND), dimension(:,:), allocatable :: &
         w_dudzTopEdge  ! w*du/dz at top of edge

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** Set return error code and return early if not turned on
      !*** Start relevant timer

      err = 0
      if (velVadvOff) return
      call mpas_timer_start("vel vadv")


      allocate(w_dudzTopEdge(nVertLevels+1,nEdgesAll))
      !$acc enter data create(w_dudzTopEdge)

#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(cellsOnEdge, minLevelEdgeBot, maxLevelEdgeTop, w_dudzTopEdge, &
      !$acc            vertAleTransportTop, normalVelocity, &
      !$acc            layerThickEdgeFlux) &
      !$acc    private(cell1, cell2, k, kmin, kmax, wAvg)
#else
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp    private(cell1, cell2, k, kmin, kmax, wAvg)
#endif
      do iEdge = 1, nEdgesOwned
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         kmin  = minLevelEdgeBot(iEdge)
         kmax  = maxLevelEdgeTop(iEdge)

         w_dudzTopEdge(kmin,iEdge) = 0.0_RKIND ! flux is zero at top
         do k = kmin+1,kmax

            ! Average w from cell center to edge
            wAvg = 0.5_RKIND*(vertAleTransportTop(k,cell1) + &
                              vertAleTransportTop(k,cell2))

            ! compute dudz at vertical interface with first order derivative.
            w_dudzTopEdge(k,iEdge) = wAvg* &
                                     (normalVelocity(k-1,iEdge) - &
                                      normalVelocity(k,  iEdge))/ &
                         (0.5_RKIND*(layerThickEdgeFlux(k-1,iEdge) + &
                                     layerThickEdgeFlux(k  ,iEdge)))
         end do
         ! transport is zero at bottom
         w_dudzTopEdge(kmax+1,iEdge) = 0.0_RKIND
      end do
#ifndef MPAS_OPENACC
      !$omp end do
#endif

      ! Average w*du/dz from vertical interface to vertical middle of cell
#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(minLevelEdgeBot, maxLevelEdgeTop, tend, edgeMask, w_dudzTopEdge)&
      !$acc    private(k, kmin, kmax)
#else
      !$omp do schedule(runtime) &
      !$omp    private(k, kmin, kmax)
#endif
      do iEdge = 1, nEdgesOwned
      kmin = minLevelEdgeBot(iEdge)
      kmax = maxLevelEdgeTop(iEdge)
      do k = kmin, kmax
         tend(k,iEdge) = tend(k,iEdge) - edgeMask(k,iEdge)* &
                         0.5_RKIND*(w_dudzTopEdge(k  ,iEdge) + &
                                    w_dudzTopEdge(k+1,iEdge))
      enddo
      enddo
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      !$acc exit data delete(w_dudzTopEdge)
      deallocate(w_dudzTopEdge)

      call mpas_timer_stop("vel vadv")

   !--------------------------------------------------------------------

   end subroutine ocn_vel_vadv_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_vadv_init
!
!> \brief   Initializes ocean momentum vertical advection
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  This routine initializes a variety of quantities related to
!>  vertical velocity advection in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_vadv_init(err)!{{{

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< [out] error flag

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** initialize error code and default values

      err = 0
      velVadvOff = .true.

      !*** Set values based on input configuration

      if (config_vert_coord_movement == 'impermeable_interfaces') then
         velVadvOff = .true.
      else
         velVadvOff = .false.
      end if

      if (config_disable_vel_vadv ) velVadvOff = .true.

   !--------------------------------------------------------------------

   end subroutine ocn_vel_vadv_init!}}}

!***********************************************************************

end module ocn_vel_vadv

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
