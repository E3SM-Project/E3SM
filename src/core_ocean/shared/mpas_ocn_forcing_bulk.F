! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_forcing_bulk
!
!> \brief MPAS ocean bulk forcing
!> \author Doug Jacobsen
!> \date   04/25/12
!> \details
!>  This module contains routines for building the forcing arrays,
!>  if bulk forcing is used.
!
!-----------------------------------------------------------------------

module ocn_forcing_bulk

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timekeeping
   use ocn_constants

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

   public :: ocn_forcing_bulk_build_arrays, &
             ocn_forcing_bulk_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   real (kind=RKIND) :: refDensity

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_forcing_bulk_build_arrays
!
!> \brief   Determines the forcing array used for the bulk forcing.
!> \author  Doug Jacobsen
!> \date    04/25/12
!> \details 
!>  This routine computes the forcing arrays used later in MPAS.
!
!-----------------------------------------------------------------------

   subroutine ocn_forcing_bulk_build_arrays(meshPool, forcingPool, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: meshPool !< Input: mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(inout) :: forcingPool !< Input: Forcing information

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: Error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, cell1, cell2
      integer :: iCell, k
      integer, pointer :: index_temperature_flux, index_salinity_flux
      integer, pointer :: nCells, nEdges

      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND) :: meridionalAverage, zonalAverage
      real (kind=RKIND), dimension(:), pointer :: angleEdge
      real (kind=RKIND), dimension(:), pointer :: windStressZonal, windStressMeridional
      real (kind=RKIND), dimension(:), pointer :: latentHeatFlux, sensibleHeatFlux, longWaveHeatFluxUp, longWaveHeatFluxDown, evaporationFlux, seaIceHeatFlux, snowFlux
      real (kind=RKIND), dimension(:), pointer :: seaIceFreshWaterFlux, seaIceSalinityFlux, riverRunoffFlux, iceRunoffFlux
      real (kind=RKIND), dimension(:), pointer :: shortWaveHeatFlux, penetrativeTemperatureFlux

      real (kind=RKIND), dimension(:), pointer :: rainFlux
      real (kind=RKIND), dimension(:), pointer :: seaSurfacePressure, iceFraction

      real (kind=RKIND), dimension(:), pointer :: surfaceThicknessFlux, surfaceWindStress, surfaceWindStressMagnitude
      real (kind=RKIND), dimension(:,:), pointer :: surfaceTracerFlux

      err = 0

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)

      call mpas_pool_get_array(meshPool, 'angleEdge', angleEdge)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

      call mpas_pool_get_dimension(forcingPool, 'index_surfaceTemperatureFlux', index_temperature_flux)
      call mpas_pool_get_dimension(forcingPool, 'index_surfaceSalinityFlux', index_salinity_flux)

      call mpas_pool_get_array(forcingPool, 'surfaceWindStress', surfaceWindStress)
      call mpas_pool_get_array(forcingPool, 'surfaceWindStressMagnitude', surfaceWindStressMagnitude)
      call mpas_pool_get_array(forcingPool, 'windStressZonal', windStressZonal)
      call mpas_pool_get_array(forcingPool, 'windStressMeridional', windStressMeridional)
      call mpas_pool_get_array(forcingPool, 'latentHeatFlux', latentHeatFlux)
      call mpas_pool_get_array(forcingPool, 'sensibleHeatFlux', sensibleHeatFlux)
      call mpas_pool_get_array(forcingPool, 'longWaveHeatFluxUp', longWaveHeatFluxUp)
      call mpas_pool_get_array(forcingPool, 'longWaveHeatFluxDown', longWaveHeatFluxDown)
      call mpas_pool_get_array(forcingPool, 'evaporationFlux', evaporationFlux)
      call mpas_pool_get_array(forcingPool, 'seaIceHeatFlux', seaIceHeatFlux)
      call mpas_pool_get_array(forcingPool, 'snowFlux', snowFlux)
      call mpas_pool_get_array(forcingPool, 'shortWaveHeatFlux', shortWaveHeatFlux)

      call mpas_pool_get_array(forcingPool, 'seaIceFreshWaterFlux', seaIceFreshWaterFlux)
      call mpas_pool_get_array(forcingPool, 'seaIceSalinityFlux', seaIceSalinityFlux)
      call mpas_pool_get_array(forcingPool, 'riverRunoffFlux', riverRunoffFlux)
      call mpas_pool_get_array(forcingPool, 'iceRunoffFlux', iceRunoffFlux)

      call mpas_pool_get_array(forcingPool, 'rainFlux', rainFlux)

      call mpas_pool_get_array(forcingPool, 'seaSurfacePressure', seaSurfacePressure)
      call mpas_pool_get_array(forcingPool, 'iceFraction', iceFraction)

      call mpas_pool_get_array(forcingPool, 'surfaceThicknessFlux', surfaceThicknessFlux)
      call mpas_pool_get_array(forcingPool, 'surfaceTracerFlux', surfaceTracerFlux)
      call mpas_pool_get_array(forcingPool, 'penetrativeTemperatureFlux', penetrativeTemperatureFlux)

      ! Convert CESM wind stress to MPAS-O windstress
      do iEdge = 1, nEdges
        cell1 = cellsOnEdge(1, iEdge)
        cell2 = cellsOnEdge(2, iEdge)

        zonalAverage = 0.5 * (windStressZonal(cell1) + windStressZonal(cell2))
        meridionalAverage = 0.5 * (windStressMeridional(cell1) + windStressMeridional(cell2))

        surfaceWindStress(iEdge) = cos(angleEdge(iEdge)) * zonalAverage + sin(angleEdge(iEdge)) * meridionalAverage
      end do


      ! Build surface fluxes at cell centers
      do iCell = 1, nCells
        surfaceWindStressMagnitude(iCell) = sqrt(windStressZonal(iCell)**2 + windStressMeridional(iCell)**2)
        surfaceTracerFlux(index_temperature_flux, iCell) = (latentHeatFlux(iCell) + sensibleHeatFlux(iCell) + longWaveHeatFluxUp(iCell) + longWaveHeatFluxDown(iCell) &
                                                           + seaIceHeatFlux(iCell) - (snowFlux(iCell) + iceRunoffFlux(iCell)) * latent_heat_fusion_mks) * hflux_factor

        surfaceTracerFlux(index_salinity_flux, iCell) = seaIceSalinityFlux(iCell) * sflux_factor

        surfaceThicknessFlux(iCell) = ( snowFlux(iCell) + rainFlux(iCell) + evaporationFlux(iCell) + seaIceFreshWaterFlux(iCell) + iceRunoffFlux(iCell) + riverRunoffFlux(iCell) ) / refDensity
      end do

      penetrativeTemperatureFlux = shortWaveHeatFlux * hflux_factor

   end subroutine ocn_forcing_bulk_build_arrays!}}}

!***********************************************************************
!
!  routine ocn_forcing_bulk_init
!
!> \brief   Initializes bulk forcing module
!> \author  Doug Jacobsen
!> \date    04/25/12
!> \details 
!>  This routine initializes the bulk forcing module.
!
!-----------------------------------------------------------------------

   subroutine ocn_forcing_bulk_init(err)!{{{

      integer, intent(out) :: err !< Output: error flag

      real (kind=RKIND), pointer :: config_density0

      err = 0


      call mpas_pool_get_config(ocnConfigs, 'config_density0', config_density0)

      refDensity = config_density0

   end subroutine ocn_forcing_bulk_init!}}}

!***********************************************************************

end module ocn_forcing_bulk


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
