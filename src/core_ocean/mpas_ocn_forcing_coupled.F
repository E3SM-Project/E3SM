! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_forcing_coupled
!
!> \brief MPAS ocean coupled forcing
!> \author Doug Jacobsen
!> \date   04/25/12
!> \version SVN:$Id:$
!> \details
!>  This module contains routines for building the forcing arrays,
!>  if coupled forcing is used.
!
!-----------------------------------------------------------------------

module ocn_forcing_coupled

   use mpas_kind_types
   use mpas_grid_types
   use mpas_configure
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

   public :: ocn_forcing_coupled_build_arrays, &
             ocn_forcing_coupled_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: coupledForcingOn !< Flag to turn on/off coupled forcing

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_build_forcing_arrays
!
!> \brief   Determines the forcing array used for the coupled forcing.
!> \author  Doug Jacobsen
!> \date    04/25/12
!> \version SVN:$Id$
!> \details 
!>  This routine computes the forcing arrays used later in MPAS.
!
!-----------------------------------------------------------------------

   subroutine ocn_forcing_coupled_build_arrays(timeStamp, mesh, forcing, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type(MPAS_Time_type), intent(in) :: timeStamp

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: mesh !< Input: mesh information

      type (forcing_type), intent(inout) :: forcing !< Input: Forcing information


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
      integer :: index_temperature_flux, index_salinity_flux

      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND) :: meridionalAverage, zonalAverage
      real (kind=RKIND), dimension(:), pointer :: angleEdge
      real (kind=RKIND), dimension(:), pointer :: windStressZonal, windStressMeridional
      real (kind=RKIND), dimension(:), pointer :: latentHeatFlux, sensibleHeatFlux, longWaveHeatFluxUp, longWaveHeatFluxDown, evaporationFlux, seaIceHeatFlux, snowFlux
      real (kind=RKIND), dimension(:), pointer :: seaIceFreshWaterFlux, seaIceSalinityFlux, riverRunoffFlux, iceRunoffFlux
      real (kind=RKIND), dimension(:), pointer :: shortWaveHeatFlux, penetrativeTemperatureFlux

      real (kind=RKIND), dimension(:), pointer :: rainFlux
      real (kind=RKIND), dimension(:), pointer :: seaSurfacePressure, iceFraction

      real (kind=RKIND), dimension(:), pointer :: surfaceMassFlux, surfaceWindStress
      real (kind=RKIND), dimension(:,:), pointer :: surfaceTracerFlux

      if(.not.coupledForcingOn) return

      angleEdge => mesh % angleEdge % array
      cellsOnEdge => mesh % cellsOnEdge % array

      index_temperature_flux = forcing % index_surfaceTemperatureflux
      index_salinity_flux = forcing % index_surfaceSalinityFlux

      index_temperature_flux = forcing % index_surfaceTemperatureFlux
      index_salinity_flux = forcing % index_surfaceSalinityFlux

      surfaceWindStress    => forcing % surfaceWindStress % array
      windStressZonal      => forcing % windStressZonal % array
      windStressMeridional => forcing % windStressMeridional % array
      latentHeatFlux       => forcing % latentHeatFlux % array
      sensibleHeatFlux     => forcing % sensibleHeatFlux % array
      longWaveHeatFluxUp   => forcing % longWaveHeatFluxUp % array
      longWaveHeatFluxDown => forcing % longWaveHeatFluxDown % array
      evaporationFlux      => forcing % evaporationFlux % array
      seaIceHeatFlux       => forcing % seaIceHeatFlux % array
      snowFlux             => forcing % snowFlux % array
      shortWaveHeatFlux    => forcing % shortWaveHeatFlux % array

      seaIceFreshWaterFlux => forcing % seaIceFreshWaterFlux % array
      seaIceSalinityFlux   => forcing % seaIceSalinityFlux % array
      riverRunoffFlux      => forcing % riverRunoffFlux % array
      iceRunoffFlux        => forcing % iceRunoffFlux % array

      rainFlux             => forcing % rainFlux % array

      seaSurfacePressure   => forcing % seaSurfacePressure % array
      iceFraction          => forcing % iceFraction % array

      surfaceMassFlux      => forcing % surfaceMassFlux % array
      surfaceTracerFlux    => forcing % surfaceTracerFlux % array
      penetrativeTemperatureFlux  => forcing % penetrativeTemperatureFlux % array

      ! Convert CESM wind stress to MPAS-O windstress
      do iEdge = 1, mesh % nEdges
        cell1 = cellsOnEdge(1, iEdge)
        cell2 = cellsOnEdge(2, iEdge)

        zonalAverage = 0.5 * (windStressZonal(cell1) + windStressZonal(cell2))
        meridionalAverage = 0.5 * (windStressMeridional(cell1) + windStressMeridional(cell2))

        surfaceWindStress(iEdge) = cos(angleEdge(iEdge)) * zonalAverage + sin(angleEdge(iEdge)) * meridionalAverage
      end do


      ! Build surface fluxes at cell centers
      do iCell = 1, mesh % nCells
        surfaceTracerFlux(index_temperature_flux, iCell) = (latentHeatFlux(iCell) + sensibleHeatFlux(iCell) + longWaveHeatFluxUp(iCell) + longWaveHeatFluxDown(iCell) &
                                                           + seaIceHeatFlux(iCell) - (snowFlux(iCell) + iceRunoffFlux(iCell)) * latent_heat_fusion_mks) * hflux_factor

        surfaceTracerFlux(index_salinity_flux, iCell) = seaIceSalinityFlux(iCell) * sflux_factor

        surfaceMassFlux(iCell) = snowFlux(iCell) + rainFlux(iCell) + evaporationFlux(iCell) + seaIceFreshWaterFlux(iCell) + iceRunoffFlux(iCell) + riverRunoffFlux(iCell)
      end do

      penetrativeTemperatureFlux = shortWaveHeatFlux * hflux_factor

   end subroutine ocn_forcing_coupled_build_arrays!}}}

!***********************************************************************
!
!  routine ocn_forcing_coupled_init
!
!> \brief   Initializes coupled forcing module
!> \author  Doug Jacobsen
!> \date    04/25/12
!> \version SVN:$Id$
!> \details 
!>  This routine initializes the coupled forcing module.
!
!-----------------------------------------------------------------------

   subroutine ocn_forcing_coupled_init(err)!{{{

      integer, intent(out) :: err !< Output: error flag

      err = 0

      coupledForcingOn = .false.

      if(config_use_coupled_forcing) then
        coupledForcingOn = .true.
      end if

   end subroutine ocn_forcing_coupled_init!}}}

!***********************************************************************

end module ocn_forcing_coupled


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
