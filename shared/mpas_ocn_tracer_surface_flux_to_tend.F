! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_surface_flux_to_tend
!
!> \brief MPAS ocean tracer surface flux
!> \author Doug Jacobsen
!> \date   12/17/12
!> \details
!>  This module contains the routine for computing
!>  surface flux tendencies.
!
!-----------------------------------------------------------------------

module ocn_tracer_surface_flux_to_tend

   use mpas_timer
   use mpas_derived_types
   use mpas_pool_routines

   use ocn_constants
   use ocn_forcing

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

   public :: ocn_tracer_surface_flux_tend, &
             ocn_tracer_surface_flux_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: surfaceTracerFluxOn

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_surface_flux_tend
!
!> \brief   Computes tendency term for surface fluxes
!> \author  Doug Jacobsen
!> \date    12/17/12
!> \details
!>  This routine computes the tendency for tracers based on surface fluxes.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_surface_flux_tend(meshPool, fractionAbsorbed, fractionAbsorbedRunoff, layerThickness,  &
      surfaceTracerFlux, surfaceTracerFluxRunoff, tend, err)!{{{
      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
        layerThickness !< Input: Layer thickness

      real (kind=RKIND), dimension(:,:), intent(in) :: &
        surfaceTracerFlux !< Input: surface tracer fluxes

      real (kind=RKIND), dimension(:,:), intent(in) :: &
        fractionAbsorbed !< Input: Coefficients for the application of surface fluxes

      real (kind=RKIND), dimension(:,:), intent(in), pointer :: &
        surfaceTracerFluxRunoff !< Input: surface tracer fluxes from river runoff

      real (kind=RKIND), dimension(:,:), intent(in) :: &
        fractionAbsorbedRunoff !< Input: Coefficients for the application of surface fluxes due to river runoff

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         tend          !< Input/Output: velocity tendency

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

      integer :: iCell, k, iTracer, nTracers, nCells
      integer, pointer :: nVertLevels
      integer, dimension(:), pointer :: nCellsArray
      integer, dimension(:), pointer :: maxLevelCell

      real (kind=RKIND) :: remainingFlux

      err = 0

      if (.not. surfaceTracerFluxOn) return

      call mpas_timer_start("surface_tracer_flux")

      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      nTracers = size(tend, dim=1)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

      nCells = nCellsArray( 1 )

      !$omp do schedule(runtime) private(remainingFlux, k, iTracer)
      do iCell = 1, nCells
        remainingFlux = 1.0_RKIND
        do k = 1, maxLevelCell(iCell)
          remainingFlux = remainingFlux - fractionAbsorbed(k, iCell)

          do iTracer = 1, nTracers
            tend(iTracer, k, iCell) = tend(iTracer, k, iCell) + surfaceTracerFlux(iTracer, iCell) * fractionAbsorbed(k, iCell)
          end do
        end do

        if(maxLevelCell(iCell) > 0 .and. remainingFlux > 0.0_RKIND) then
          do iTracer = 1, nTracers
            tend(iTracer, maxLevelCell(iCell), iCell) = tend(iTracer, maxLevelCell(iCell), iCell) &
                                                      + surfaceTracerFlux(iTracer, iCell) * remainingFlux
          end do
        end if
      end do
      !$omp end do

      call mpas_timer_stop("surface_tracer_flux")

      ! now do runoff component

      if (associated(surfaceTracerFluxRunoff)) then
        call mpas_timer_start("surface_tracer_runoff_flux")

        !$omp do schedule(runtime) private(remainingFlux, k, iTracer)
        do iCell = 1, nCells
          remainingFlux = 1.0_RKIND
          do k = 1, maxLevelCell(iCell)
            remainingFlux = remainingFlux - fractionAbsorbedRunoff(k, iCell)

            do iTracer = 1, nTracers
              tend(iTracer, k, iCell) = tend(iTracer, k, iCell) +  &
                 surfaceTracerFluxRunoff(iTracer, iCell) * fractionAbsorbedRunoff(k, iCell)
            end do
          end do

          if(maxLevelCell(iCell) > 0 .and. remainingFlux > 0.0_RKIND) then
            do iTracer = 1, nTracers
              tend(iTracer, maxLevelCell(iCell), iCell) = tend(iTracer, maxLevelCell(iCell), iCell) &
                                                        + surfaceTracerFluxRunoff(iTracer, iCell) * remainingFlux
            end do
          end if
        end do
        !$omp end do

        call mpas_timer_stop("surface_tracer_runoff_flux")
      end if


   !--------------------------------------------------------------------

   end subroutine ocn_tracer_surface_flux_tend!}}}

!***********************************************************************
!
!  routine ocn_tracer_surface_flux_init
!
!> \brief   Initializes ocean tracer surface flux quantities
!> \author  Doug Jacobsen
!> \date    12/17/12
!> \details
!>  This routine initializes quantities related to surface fluxes in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_surface_flux_init(err)!{{{

   !--------------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      logical, pointer :: config_disable_tr_sflux

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_disable_tr_sflux', config_disable_tr_sflux)

      surfaceTracerFluxOn = .true.

      if (config_disable_tr_sflux) then
         surfaceTracerFluxOn = .false.
      end if

   end subroutine ocn_tracer_surface_flux_init!}}}

!***********************************************************************

end module ocn_tracer_surface_flux_to_tend

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
