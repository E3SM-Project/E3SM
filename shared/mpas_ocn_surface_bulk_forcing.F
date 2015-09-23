! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_surface_bulk_forcing
!
!> \brief MPAS ocean bulk forcing
!> \author Doug Jacobsen
!> \date   04/25/12
!> \details
!>  This module contains routines for building the forcing arrays,
!>  if bulk forcing is used.
!
!-----------------------------------------------------------------------

module ocn_surface_bulk_forcing

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

   public :: ocn_surface_bulk_forcing_tracers, &
             ocn_surface_bulk_forcing_vel, &
             ocn_surface_bulk_forcing_thick, &
             ocn_surface_bulk_forcing_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: bulkWindStressOn, bulkThicknessFluxOn

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_surface_bulk_forcing_tracers
!
!> \brief   Determines the tracers forcing array used for the bulk forcing.
!> \author  Doug Jacobsen
!> \date    04/25/12
!> \details 
!>  This routine computes the tracers forcing arrays used later in MPAS.
!
!-----------------------------------------------------------------------

! TRACER-CLEAN-UP
! Currently, penetrativeTemperatureFlux is built into bulk forcing.. what should we do about that?
   subroutine ocn_surface_bulk_forcing_tracers(meshPool, groupName, forcingPool, tracersSurfaceFlux, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: meshPool !< Input: mesh information
      character (len=*) :: groupName !< Input: Name of tracer group

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(inout) :: forcingPool !< Input: Forcing information
      real (kind=RKIND), dimension(:,:), intent(inout) :: tracersSurfaceFlux !< Input/Output: Surface flux for tracer group

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

      err = 0

      if ( trim(groupName) == 'activeTracers' ) then
         call ocn_surface_bulk_forcing_active_tracers(meshPool, forcingPool, tracersSurfaceFlux, err)
      end if

   end subroutine ocn_surface_bulk_forcing_tracers!}}}

!***********************************************************************
!
!  routine ocn_surface_bulk_forcing_vel
!
!> \brief   Determines the velocity forcing array used for the bulk forcing.
!> \author  Doug Jacobsen
!> \date    04/25/12
!> \details 
!>  This routine computes the velocity forcing arrays used later in MPAS.
!
!-----------------------------------------------------------------------

   subroutine ocn_surface_bulk_forcing_vel(meshPool, forcingPool, surfaceStress, surfaceStressMagnitude, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: meshPool !< Input: mesh information
      type (mpas_pool_type), intent(in) :: forcingPool !< Input: Forcing information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      real (kind=RKIND), dimension(:), intent(inout) :: surfaceStress, & !< Input/Output: Array for surface stress
                                                        surfaceStressMagnitude !< Input/Output: Array for magnitude of surface stress

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

      integer :: iEdge, cell1, cell2, iCell
      integer, pointer :: nCells, nEdges

      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND) :: meridionalAverage, zonalAverage
      real (kind=RKIND), dimension(:), pointer :: angleEdge
      real (kind=RKIND), dimension(:), pointer :: windStressZonal, windStressMeridional

      err = 0

      if ( .not. bulkWindStressOn ) return

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)

      call mpas_pool_get_array(meshPool, 'angleEdge', angleEdge)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

      call mpas_pool_get_array(forcingPool, 'windStressZonal', windStressZonal)
      call mpas_pool_get_array(forcingPool, 'windStressMeridional', windStressMeridional)

      ! Convert CESM wind stress to MPAS-O wind stress
      do iEdge = 1, nEdges
        cell1 = cellsOnEdge(1, iEdge)
        cell2 = cellsOnEdge(2, iEdge)

        zonalAverage = 0.5 * (windStressZonal(cell1) + windStressZonal(cell2))
        meridionalAverage = 0.5 * (windStressMeridional(cell1) + windStressMeridional(cell2))

        surfaceStress(iEdge) = surfaceStress(iEdge) + cos(angleEdge(iEdge)) * zonalAverage + sin(angleEdge(iEdge)) * meridionalAverage
      end do


      ! Build surface fluxes at cell centers
      do iCell = 1, nCells
        surfaceStressMagnitude(iCell) = surfaceStressMagnitude(iCell) + sqrt(windStressZonal(iCell)**2 + windStressMeridional(iCell)**2)
      end do

   end subroutine ocn_surface_bulk_forcing_vel!}}}

!***********************************************************************
!
!  routine ocn_surface_bulk_forcing_thick
!
!> \brief   Determines the thickness forcing array used for the bulk forcing.
!> \author  Doug Jacobsen
!> \date    04/25/12
!> \details 
!>  This routine computes the thickness forcing arrays used later in MPAS.
!
!-----------------------------------------------------------------------

   subroutine ocn_surface_bulk_forcing_thick(meshPool, forcingPool, surfaceThicknessFlux, err)!{{{

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
      real (kind=RKIND), dimension(:), intent(inout) :: surfaceThicknessFlux !< Input/Output: Array for surface thickness flux

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

      integer :: iCell
      integer, pointer :: nCells, nEdges

      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND), dimension(:), pointer :: evaporationFlux, snowFlux
      real (kind=RKIND), dimension(:), pointer :: seaIceFreshWaterFlux, riverRunoffFlux, iceRunoffFlux
      real (kind=RKIND), dimension(:), pointer :: rainFlux

      err = 0

      if ( .not. bulkThicknessFluxOn ) return

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      call mpas_pool_get_array(forcingPool, 'evaporationFlux', evaporationFlux)
      call mpas_pool_get_array(forcingPool, 'snowFlux', snowFlux)
      call mpas_pool_get_array(forcingPool, 'seaIceFreshWaterFlux', seaIceFreshWaterFlux)
      call mpas_pool_get_array(forcingPool, 'riverRunoffFlux', riverRunoffFlux)
      call mpas_pool_get_array(forcingPool, 'iceRunoffFlux', iceRunoffFlux)
      call mpas_pool_get_array(forcingPool, 'rainFlux', rainFlux)


      ! Build surface fluxes at cell centers
      do iCell = 1, nCells
        surfaceThicknessFlux(iCell) = surfaceThicknessFlux(iCell) + ( snowFlux(iCell) + rainFlux(iCell) + evaporationFlux(iCell) &
                                    + seaIceFreshWaterFlux(iCell) + iceRunoffFlux(iCell) &
                                    + riverRunoffFlux(iCell) ) / rho_sw
      end do

   end subroutine ocn_surface_bulk_forcing_thick!}}}

!***********************************************************************
!
!  routine ocn_surface_bulk_forcing_init
!
!> \brief   Initializes bulk forcing module
!> \author  Doug Jacobsen
!> \date    04/25/12
!> \details 
!>  This routine initializes the bulk forcing module.
!
!-----------------------------------------------------------------------

   subroutine ocn_surface_bulk_forcing_init(err)!{{{

      integer, intent(out) :: err !< Output: error flag

      logical, pointer :: config_use_bulk_wind_stress, config_use_bulk_thickness_flux

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_use_bulk_wind_stress', config_use_bulk_wind_stress)
      call mpas_pool_get_config(ocnConfigs, 'config_use_bulk_thickness_flux', config_use_bulk_thickness_flux)

      bulkWindStressOn = config_use_bulk_wind_stress
      bulkThicknessFluxOn = config_use_bulk_thickness_flux

   end subroutine ocn_surface_bulk_forcing_init!}}}

!***********************************************************************
!
! Private module subroutines
!
!***********************************************************************


!***********************************************************************
!
!  routine ocn_surface_bulk_forcing_active_tracers
!
!> \brief   Determines the active tracers forcing array used for the bulk forcing.
!> \author  Doug Jacobsen
!> \date    04/25/12
!> \details 
!>  This routine computes the active tracers forcing arrays used later in MPAS.
!
!-----------------------------------------------------------------------

   subroutine ocn_surface_bulk_forcing_active_tracers(meshPool, forcingPool, tracersSurfaceFlux, err)!{{{

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
      real (kind=RKIND), dimension(:,:), intent(inout) :: tracersSurfaceFlux

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

      integer :: iCell
      integer, pointer :: index_temperature_flux, index_salinity_flux
      integer, pointer :: nCells

      real (kind=RKIND), dimension(:), pointer :: latentHeatFlux, sensibleHeatFlux, longWaveHeatFluxUp, longWaveHeatFluxDown, seaIceHeatFlux, snowFlux
      real (kind=RKIND), dimension(:), pointer :: seaIceFreshWaterFlux, seaIceSalinityFlux, iceRunoffFlux
      real (kind=RKIND), dimension(:), pointer :: shortWaveHeatFlux, penetrativeTemperatureFlux

      err = 0

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

 ! CLEANUP
   !  call mpas_pool_get_dimension(forcingPool, 'index_temperatureSurfaceFlux', index_temperature_flux)
   !  call mpas_pool_get_dimension(forcingPool, 'index_salinitySurfaceFlux', index_salinity_flux)

      call mpas_pool_get_array(forcingPool, 'latentHeatFlux', latentHeatFlux)
      call mpas_pool_get_array(forcingPool, 'sensibleHeatFlux', sensibleHeatFlux)
      call mpas_pool_get_array(forcingPool, 'longWaveHeatFluxUp', longWaveHeatFluxUp)
      call mpas_pool_get_array(forcingPool, 'longWaveHeatFluxDown', longWaveHeatFluxDown)
      call mpas_pool_get_array(forcingPool, 'seaIceHeatFlux', seaIceHeatFlux)
      call mpas_pool_get_array(forcingPool, 'snowFlux', snowFlux)
      call mpas_pool_get_array(forcingPool, 'shortWaveHeatFlux', shortWaveHeatFlux)

      call mpas_pool_get_array(forcingPool, 'seaIceSalinityFlux', seaIceSalinityFlux)
      call mpas_pool_get_array(forcingPool, 'iceRunoffFlux', iceRunoffFlux)

      call mpas_pool_get_array(forcingPool, 'penetrativeTemperatureFlux', penetrativeTemperatureFlux)

      ! Build surface fluxes at cell centers
      ! CLEANUP
      do iCell = 1, nCells
        tracersSurfaceFlux(1, iCell) = tracersSurfaceFlux(1, iCell) &
                                                           + (latentHeatFlux(iCell) + sensibleHeatFlux(iCell) + longWaveHeatFluxUp(iCell) + longWaveHeatFluxDown(iCell) &
                                                           + seaIceHeatFlux(iCell) - (snowFlux(iCell) + iceRunoffFlux(iCell)) * latent_heat_fusion_mks) * hflux_factor

        tracersSurfaceFlux(2, iCell) = tracersSurfaceFlux(2, iCell)  &
                                                        + seaIceSalinityFlux(iCell) * sflux_factor
      end do

      ! TRACER-CLEAN-UP
      ! Do we want this here still?
      penetrativeTemperatureFlux = shortWaveHeatFlux * hflux_factor

   end subroutine ocn_surface_bulk_forcing_active_tracers!}}}

end module ocn_surface_bulk_forcing


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
