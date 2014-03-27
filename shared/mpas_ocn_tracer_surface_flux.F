! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_surface_flux
!
!> \brief MPAS ocean tracer surface flux
!> \author Doug Jacobsen
!> \date   12/17/12
!> \details
!>  This module contains the routine for computing 
!>  surface flux tendencies.  
!
!-----------------------------------------------------------------------

module ocn_tracer_surface_flux

   use mpas_grid_types

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

   subroutine ocn_tracer_surface_flux_tend(meshPool, transmissionCoefficients, layerThickness, surfaceTracerFlux, tend, err)!{{{
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
        transmissionCoefficients !< Input: Coefficients for the application of surface fluxes

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

      integer :: iCell, k, iTracer, nTracers
      integer, pointer :: nCells, nVertLevels
      integer, dimension(:), pointer :: maxLevelCell
      integer, dimension(:,:), pointer :: cellMask

      real (kind=RKIND) :: remainingFlux

      err = 0

      if (.not. surfaceTracerFluxOn) return

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      nTracers = size(tend, dim=1)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'cellMask', cellMask)

      do iCell = 1, nCells
        remainingFlux = 1.0_RKIND
        do k = 1, maxLevelCell(iCell)
          remainingFlux = remainingFlux - transmissionCoefficients(k, iCell)

          do iTracer = 1, nTracers
            tend(iTracer, k, iCell) = tend(iTracer, k, iCell) + cellMask(k, icell) * surfaceTracerFlux(iTracer, iCell) * transmissionCoefficients(k, iCell)
          end do
        end do

        if(maxLevelCell(iCell) > 0 .and. remainingFlux > 0.0_RKIND) then
          do iTracer = 1, nTracers
            tend(iTracer, maxLevelCell(iCell), iCell) = tend(iTracer, maxLevelCell(iCell), iCell) + cellMask(k, iCell) * surfaceTracerFlux(iTracer, iCell) * remainingFlux
          end do
        end if
      end do

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
      character (len=StrKIND), pointer :: config_forcing_type

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_disable_tr_sflux', config_disable_tr_sflux)
      call mpas_pool_get_config(ocnConfigs, 'config_forcing_type', config_forcing_type)

      surfaceTracerFluxOn = .true.

      if (config_disable_tr_sflux) then
         surfaceTracerFluxOn = .false.
      end if

      if (config_forcing_type == trim('off')) then
         surfaceTracerFluxOn = .false.
      end if

   end subroutine ocn_tracer_surface_flux_init!}}}

!***********************************************************************

end module ocn_tracer_surface_flux

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
