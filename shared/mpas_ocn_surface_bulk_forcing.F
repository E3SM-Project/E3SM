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

   use mpas_timer
   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timekeeping
   use ocn_constants
   use ocn_equation_of_state

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

   subroutine ocn_surface_bulk_forcing_tracers(meshPool, groupName, forcingPool, tracerGroup,   &
      tracersSurfaceFlux, tracersSurfaceFluxRunoff, tracersSurfaceFluxRemoved, dt, layerThickness, err)!{{{

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
      real (kind=RKIND), dimension(:,:,:), intent(inout) :: tracerGroup
      real (kind=RKIND), dimension(:,:), intent(inout) :: tracersSurfaceFlux !< Input/Output: Surface flux for tracer group
      real (kind=RKIND), dimension(:,:), intent(inout) ::   &
         tracersSurfaceFluxRunoff !< Input/Output: Surface flux for tracer group due to river runoff
      real (kind=RKIND), dimension(:,:), intent(inout) ::   &
         tracersSurfaceFluxRemoved !< Input/Output: Accumulator for ignored Surface flux for tracer group
      real (kind=RKIND), dimension(:,:), intent(in) :: layerThickness
      real (kind=RKIND), intent(in) :: dt

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

      call mpas_timer_start("bulk_" // trim(groupName))
      if ( trim(groupName) == 'activeTracers' ) then
         call ocn_surface_bulk_forcing_active_tracers(meshPool, forcingPool, tracerGroup,  &
            tracersSurfaceFlux, tracersSurfaceFluxRunoff, tracersSurfaceFluxRemoved, layerThickness, dt, err)
      end if
      call mpas_timer_stop("bulk_" // trim(groupName))

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

      integer :: iEdge, cell1, cell2, iCell, nCells, nEdges
      integer, dimension(:), pointer :: nCellsArray, nEdgesArray

      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND) :: meridionalAverage, zonalAverage
      real (kind=RKIND), dimension(:), pointer :: angleEdge
      real (kind=RKIND), dimension(:), pointer :: windStressZonal, windStressMeridional

      err = 0

      if ( .not. bulkWindStressOn ) return

      call mpas_timer_start("bulk_ws", .false.)

      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)

      call mpas_pool_get_array(meshPool, 'angleEdge', angleEdge)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

      call mpas_pool_get_array(forcingPool, 'windStressZonal', windStressZonal)
      call mpas_pool_get_array(forcingPool, 'windStressMeridional', windStressMeridional)

      nEdges = nEdgesArray( 3 )
      nCells = nCellsArray( 2 )

      ! Convert CESM wind stress to MPAS-O wind stress
      !$omp do schedule(runtime) private(cell1, cell2, zonalAverage, meridionalAverage)
      do iEdge = 1, nEdges
        cell1 = cellsOnEdge(1, iEdge)
        cell2 = cellsOnEdge(2, iEdge)

        zonalAverage = 0.5_RKIND * (windStressZonal(cell1) + windStressZonal(cell2))
        meridionalAverage = 0.5_RKIND * (windStressMeridional(cell1) + windStressMeridional(cell2))

        surfaceStress(iEdge) = surfaceStress(iEdge) + cos(angleEdge(iEdge)) * zonalAverage + sin(angleEdge(iEdge)) &
                             * meridionalAverage
      end do
      !$omp end do

      ! Build surface fluxes at cell centers
      !$omp do schedule(runtime)
      do iCell = 1, nCells
        surfaceStressMagnitude(iCell) = surfaceStressMagnitude(iCell) + sqrt( windStressZonal(iCell)**2 &
                                      + windStressMeridional(iCell)**2 )
      end do
      !$omp end do

      call mpas_timer_stop("bulk_ws")

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

   subroutine ocn_surface_bulk_forcing_thick(meshPool, forcingPool, surfaceThicknessFlux, surfaceThicknessFluxRunoff, err)!{{{

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
      real (kind=RKIND), dimension(:), intent(inout) ::  &
         surfaceThicknessFluxRunoff !< Input/Output: Array for surface thickness flux due to river runoff

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

      integer :: iCell, nCells, nEdges
      integer, dimension(:), pointer :: nCellsArray

      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND), dimension(:), pointer :: evaporationFlux, snowFlux
      real (kind=RKIND), dimension(:), pointer :: seaIceFreshWaterFlux, riverRunoffFlux, iceRunoffFlux
      real (kind=RKIND), dimension(:), pointer :: rainFlux

      err = 0

      if ( .not. bulkThicknessFluxOn ) return

      call mpas_timer_start("bulk_thick", .false.)

      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)

      call mpas_pool_get_array(forcingPool, 'evaporationFlux', evaporationFlux)
      call mpas_pool_get_array(forcingPool, 'snowFlux', snowFlux)
      call mpas_pool_get_array(forcingPool, 'seaIceFreshWaterFlux', seaIceFreshWaterFlux)
      call mpas_pool_get_array(forcingPool, 'riverRunoffFlux', riverRunoffFlux)
      call mpas_pool_get_array(forcingPool, 'iceRunoffFlux', iceRunoffFlux)
      call mpas_pool_get_array(forcingPool, 'rainFlux', rainFlux)

      nCells = nCellsArray( 2 )

      ! Build surface fluxes at cell centers
      !$omp do schedule(runtime)
      do iCell = 1, nCells
        surfaceThicknessFlux(iCell) = surfaceThicknessFlux(iCell) + ( snowFlux(iCell) + rainFlux(iCell) + evaporationFlux(iCell) &
                                    + seaIceFreshWaterFlux(iCell) + iceRunoffFlux(iCell) ) / rho_sw
        surfaceThicknessFluxRunoff(iCell) = riverRunoffFlux(iCell) / rho_sw
      end do
      !$omp end do

      call mpas_timer_stop("bulk_thick")

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

   subroutine ocn_surface_bulk_forcing_active_tracers(meshPool, forcingPool, tracerGroup,  &
      tracersSurfaceFlux, tracersSurfaceFluxRunoff, tracersSurfaceFluxRemoved, layerThickness, dt, err)!{{{

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
      real (kind=RKIND), dimension(:,:), intent(inout) :: tracersSurfaceFluxRunoff
      real (kind=RKIND), dimension(:,:), intent(inout) :: tracersSurfaceFluxRemoved
      real (kind=RKIND), dimension(:,:,:), intent(inout) :: tracerGroup
      real (kind=RKIND), dimension(:,:), intent(in) :: layerThickness
      real (kind=RKIND), intent(in) :: dt

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

      integer :: iCell, nCells
      integer, pointer :: index_temperature_flux, index_salinity_flux
      integer, dimension(:), pointer :: nCellsArray

      type(mpas_pool_type),pointer :: tracersSurfaceFluxPool

      real (kind=RKIND), dimension(:), pointer :: latentHeatFlux, sensibleHeatFlux, longWaveHeatFluxUp, longWaveHeatFluxDown, &
                                                  seaIceHeatFlux, evaporationFlux, riverRunoffFlux
      real (kind=RKIND), dimension(:), pointer :: seaIceFreshWaterFlux, seaIceSalinityFlux, iceRunoffFlux
      real (kind=RKIND), dimension(:), pointer :: shortWaveHeatFlux, penetrativeTemperatureFlux
      real (kind=RKIND), dimension(:), pointer :: snowFlux, rainFlux
      real (kind=RKIND) :: requiredSalt, allowedSalt

      err = 0

      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)

      call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceFlux',tracersSurfaceFluxPool)

      call mpas_pool_get_dimension(tracersSurfaceFluxPool, 'index_temperatureSurfaceFlux', index_temperature_flux)
      call mpas_pool_get_dimension(tracersSurfaceFluxPool, 'index_salinitySurfaceFlux', index_salinity_flux)

      call mpas_pool_get_array(forcingPool, 'latentHeatFlux', latentHeatFlux)
      call mpas_pool_get_array(forcingPool, 'sensibleHeatFlux', sensibleHeatFlux)
      call mpas_pool_get_array(forcingPool, 'longWaveHeatFluxUp', longWaveHeatFluxUp)
      call mpas_pool_get_array(forcingPool, 'longWaveHeatFluxDown', longWaveHeatFluxDown)
      call mpas_pool_get_array(forcingPool, 'seaIceHeatFlux', seaIceHeatFlux)
      call mpas_pool_get_array(forcingPool, 'rainFlux', rainFlux)
      call mpas_pool_get_array(forcingPool, 'snowFlux', snowFlux)
      call mpas_pool_get_array(forcingPool, 'shortWaveHeatFlux', shortWaveHeatFlux)
      call mpas_pool_get_array(forcingPool, 'evaporationFlux', evaporationFlux)

      call mpas_pool_get_array(forcingPool, 'seaIceFreshWaterFlux', seaIceFreshWaterFlux)
      call mpas_pool_get_array(forcingPool, 'seaIceSalinityFlux', seaIceSalinityFlux)
      call mpas_pool_get_array(forcingPool, 'iceRunoffFlux', iceRunoffFlux)
      call mpas_pool_get_array(forcingPool, 'riverRunoffFlux', riverRunoffFlux)
      call mpas_pool_get_array(forcingPool, 'penetrativeTemperatureFlux', penetrativeTemperatureFlux)

      nCells = nCellsArray( 2 )

      ! Build surface fluxes at cell centers
      !$omp do schedule(runtime) private(allowedSalt, requiredSalt)
      do iCell = 1, nCells
        tracersSurfaceFlux(index_temperature_flux, iCell) = tracersSurfaceFlux(index_temperature_flux, iCell) &
                                                           + (latentHeatFlux(iCell) + sensibleHeatFlux(iCell) &
                                                           + longWaveHeatFluxUp(iCell) + longWaveHeatFluxDown(iCell) &
                                                           + seaIceHeatFlux(iCell) - (snowFlux(iCell) + iceRunoffFlux(iCell)) &
                                                           * latent_heat_fusion_mks) * hflux_factor

        ! Negative seaIceSalinityFlux is an extraction of salt from the ocean
        ! So, we negate seaIceSalinityFlux when determining how much salt this flux needs.
        requiredSalt = - seaIceSalinityFlux(iCell) * sflux_factor * dt / layerThickness(1, iCell)
        allowedSalt = min( 4.0_RKIND, tracerGroup(index_salinity_flux, 1, iCell) )

        if ( allowedSalt < requiredSalt ) then
           tracersSurfaceFluxRemoved(index_salinity_flux, iCell) = tracersSurfaceFluxRemoved(index_salinity_flux, iCell)  &
                                                                 + ( 1 - allowedSalt / requiredSalt ) * seaIceSalinityFlux(iCell) &
                                                                 * sflux_factor

           tracersSurfaceFlux(index_salinity_flux, iCell) = tracersSurfaceFlux(index_salinity_flux, iCell)  &
                                                          + ( allowedSalt / requiredSalt ) * seaIceSalinityFlux(iCell) &
                                                          * sflux_factor
        else
           tracersSurfaceFlux(index_salinity_flux, iCell) = tracersSurfaceFlux(index_salinity_flux, iCell)  &
                                                          + seaIceSalinityFlux(iCell) * sflux_factor
        end if
      end do
      !$omp end do
      ! assume that snow comes in at 0 C

      ! Surface fluxes of water have an associated heat content, but the coupled system does not account for this
      ! Assume surface fluxes of water have a temperature dependent on the incoming mass flux.
      ! Assume surface fluxes of water have zero salinity. So the RHS forcing is zero for salinity.
      ! Only include this heat forcing when bulk thickness is turned on
      ! indices on tracerGroup are (iTracer, iLevel, iCell)
      if (bulkThicknessFluxOn) then
         !$omp do schedule(runtime)
         do iCell = 1, nCells

           ! Accumulate fluxes that use the surface temperature
           tracersSurfaceFlux(index_temperature_flux, iCell) = tracersSurfaceFlux(index_temperature_flux, iCell) &
                      + (rainFlux(iCell) + evaporationFlux(iCell)) * tracerGroup(index_temperature_flux,1,iCell) / rho_sw

           ! Runoff can only have a minimum temperature of 0.0C, since it is fresh water.
           tracersSurfaceFluxRunoff(index_temperature_flux,iCell) = riverRunoffFlux(iCell) &
                      * max(tracerGroup(index_temperature_flux,1,iCell), 0.0_RKIND) / rho_sw

           ! Accumulate fluxes that use the freezing point
           tracersSurfaceFlux(index_temperature_flux, iCell) = tracersSurfaceFlux(index_temperature_flux, iCell) &
               + seaIceFreshWaterFlux(iCell) * ocn_freezing_temperature( tracerGroup(index_salinity_flux, 1, iCell) , &
                                                                         pressure=0.0_RKIND, &
                                                                         inLandIceCavity=.false.) / rho_sw

           ! Fields with zero temperature are not accumulated. These include:
           !    snowFlux
           !    iceRunoffFlux

         end do
         !$omp end do
      endif ! bulkThicknessFluxOn

      ! convert short wave heat flux to a temperature flux
      !$omp do schedule(runtime)
      do iCell = 1, nCells
         penetrativeTemperatureFlux(iCell) = shortWaveHeatFlux(iCell) * hflux_factor
      end do
      !$omp end do

   end subroutine ocn_surface_bulk_forcing_active_tracers!}}}

end module ocn_surface_bulk_forcing


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
