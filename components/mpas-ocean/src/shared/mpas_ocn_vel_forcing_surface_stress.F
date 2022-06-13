! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_forcing_surface_stress
!
!> \brief MPAS ocean surface stress
!> \author Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date   September 2011
!> \details
!>  This module contains the routine for computing
!>  tendencies from surface stress.
!
!-----------------------------------------------------------------------

module ocn_vel_forcing_surface_stress

   use mpas_timer

   use ocn_constants
   use ocn_config
   use ocn_forcing
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

   public :: ocn_vel_forcing_surface_stress_tend, &
             ocn_vel_forcing_surface_stress_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: &
      surfaceStressOff  !< on/off switch for surface stress

   real (kind=RKIND), parameter :: &
      maxAttDepth = -100.0_RKIND ! max attenuation depth to prevent
                                 !  underflow in exponential

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_forcing_surface_stress_tend
!
!> \brief   Computes tendency term from surface stress
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    15 September 2011
!> \details
!>  This routine computes the surface stress tendency for momentum
!>  based on current state.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_forcing_surface_stress_tend(sfcFlxAttCoeff, &
                               sfcStress, layerThickEdgeMean, tend, err)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(in) :: &
         sfcStress,    & !< [in] Wind stress at surface
         sfcFlxAttCoeff  !< [in] attenuation coefficient for sfc fluxes

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickEdgeMean  !< [in] mean thickness at edge

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

      integer ::      &
         iEdge, k,    &! loop indices for edge, vertical loops
         kmax,        &! index of deepest active edge
         cell1, cell2  ! neighbor cell addresses across edge

      real (kind=RKIND) :: &
         transCoeffTop,    &! transmission coefficent at top of edge
         transCoeffBot,    &! transmission coefficent at bottom of edge
         zTop, zBot,       &! depths at top and bottom of edge
         remainingStress,  &! stress remaining at bottom of edge
         absorb,           &! absorption fraction within layer
         attDepth,         &! depth scaled by attenuation length
         attCoeff           ! attenuation coefficient on edge

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** Initialize error code and return if not on
      !*** If turned on, start timer.

      err = 0
      if (surfaceStressOff) return
      call mpas_timer_start('vel surface stress')

#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(cellsOnEdge, minLevelEdgeBot, maxLevelEdgeTop, edgeMask, &
      !$acc            sfcFlxAttCoeff, sfcStress, layerThickEdgeMean, &
      !$acc            tend) &
      !$acc    private(k, kmax, cell1, cell2, zBot, zTop, &
      !$acc            attCoeff, attDepth, remainingStress, &
      !$acc            transCoeffTop, transCoeffBot, absorb)
#else
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp    private(k, kmax, cell1, cell2, zBot, zTop, &
      !$omp            attCoeff, attDepth, remainingStress, &
      !$omp            transCoeffTop, transCoeffBot, absorb)
#endif
      do iEdge = 1, nEdgesOwned
        zTop     = 0.0_RKIND
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)
        attCoeff = 0.5_RKIND*(sfcFlxAttCoeff(cell1) &
                            + sfcFlxAttCoeff(cell2))
        attDepth = 0.0_RKIND
        transCoeffTop = 1.0_RKIND
        remainingStress = 1.0_RKIND
        kmax = maxLevelEdgeTop(iEdge)
        do k = minLevelEdgeBot(iEdge), kmax
           zBot = zTop - layerThickEdgeMean(k,iEdge)
           attDepth = max(zBot/attCoeff, maxAttDepth)

           transCoeffBot = exp(attDepth)

           absorb = transCoeffTop - transCoeffBot
           remainingStress = remainingStress - absorb

           tend(k,iEdge) =  tend(k,iEdge) + &
                            edgeMask(k,iEdge)*sfcStress(iEdge) * &
                            absorb/rho_sw/layerThickEdgeMean(k,iEdge)

           zTop = zBot
           transCoeffTop = transCoeffBot
        enddo

        !*** if there is any remaining stress at the bottom, add it
        !*** to the bottom layer

        if ( kmax > 0 .and. remainingStress > 0.0_RKIND) then
           tend(kmax,iEdge) = tend(kmax,iEdge) &
                            + edgeMask(kmax,iEdge)*sfcStress(iEdge)* &
                              remainingStress/rho_sw/ &
                              layerThickEdgeMean(kmax,iEdge)
        end if
      enddo
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      call mpas_timer_stop('vel surface stress')

   !--------------------------------------------------------------------

   end subroutine ocn_vel_forcing_surface_stress_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_forcing_surface_stress_init
!
!> \brief   Initializes ocean surface stress forcing
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details
!>  This routine initializes quantities related to surface stress
!>  in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_forcing_surface_stress_init(err)!{{{

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< [out] error flag

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** Initialize return error code and set module defaults

      err = 0
      surfaceStressOff = .false.

      !*** Reset based on input model configuration

      if (config_disable_vel_surface_stress) surfaceStressOff = .true.

   !--------------------------------------------------------------------

   end subroutine ocn_vel_forcing_surface_stress_init!}}}

!***********************************************************************

end module ocn_vel_forcing_surface_stress

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
