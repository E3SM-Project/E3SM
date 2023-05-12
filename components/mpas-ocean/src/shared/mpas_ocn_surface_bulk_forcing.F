! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
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
   use ocn_config
   use ocn_mesh
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

   logical ::              &
      bulkWindStressOff,   &! on/off switch for bulk wind stress
      bulkThicknessFluxOff  ! on/off switch for bulk thickness flux

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
!>  This routine computes the surface stress for velocity forcing
!>  based on bulk flux formulae.
!
!-----------------------------------------------------------------------

   subroutine ocn_surface_bulk_forcing_vel( &
                               forcingPool, sfcStress, &
                               sfcStressMag, err)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         forcingPool         !< [in] Forcing structure w/ sfc stresses

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(inout) :: &
         sfcStress, & !< [inout] Array for surface stress
         sfcStressMag !< [inout] Array for magnitude of surface stress

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< [out] Error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer ::         &
         iEdge, iCell,   &! loop indices for edge, cell loops
         nCells, nEdges, &! number of cells, edges participating
         cell1, cell2     ! neighbor cell indices across edge

      real (kind=RKIND) :: &
         meridWSEdge,      &! meridional wind stress avg to edge
         zonalWSEdge        ! zonal wind stress averaged to edge

      real (kind=RKIND), dimension(:), pointer :: &
         windStressZonal, &! zonal wind stress from forcing
         windStressMerid   ! meridional wind stress from forcing

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** Set error code and exit early if bulk wind stress not on
      !*** Otherwise start timer

      err = 0
      if (bulkWindStressOff) return

      call mpas_timer_start("bulk_ws", .false.)

      call mpas_pool_get_array(forcingPool, 'windStressZonal', windStressZonal)
      call mpas_pool_get_array(forcingPool, 'windStressMeridional', windStressMerid)

      !*** Transfer data to device
      !$acc enter data &
      !$acc    copyin(windStressZonal, windStressMerid)

      !*** Forcing needed over halo regions

      nEdges = nEdgesHalo(3)
      nCells = nCellsHalo(2)

      !*** Convert coupled climate model wind stress to internal
      !*** ocean model wind stress by averaging cell-centered wind
      !*** stress to edges and then rotating vectors

#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(cellsOnEdge, angleEdge, sfcStress, &
      !$acc            windStressZonal, windStressMerid) &
      !$acc    private(cell1, cell2, zonalWSEdge, meridWSEdge)
#else
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp    private(cell1, cell2, zonalWSEdge, meridWSEdge)
#endif
      do iEdge = 1, nEdges
        cell1 = cellsOnEdge(1, iEdge)
        cell2 = cellsOnEdge(2, iEdge)

        zonalWSEdge = 0.5_RKIND * (windStressZonal(cell1) + &
                                   windStressZonal(cell2))
        meridWSEdge = 0.5_RKIND * (windStressMerid(cell1) + &
                                   windStressMerid(cell2))

        sfcStress(iEdge) = sfcStress(iEdge) + &
                           cos(angleEdge(iEdge))*zonalWSEdge + &
                           sin(angleEdge(iEdge))*meridWSEdge
      end do
#ifndef MPAS_OPENACC
      !$omp end do
#endif

      ! Calculate surface flux magnitude at cell centers

#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(sfcStressMag, windStressZonal, &
      !$acc                          windStressMerid)
#else
      !$omp do schedule(runtime)
#endif
      do iCell = 1, nCells
        sfcStressMag(iCell) = sfcStressMag(iCell) + &
                              sqrt( windStressZonal(iCell)**2 &
                                  + windStressMerid(iCell)**2 )
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      !*** Transfer any needed data back to host and delete inputs
      !$acc exit data &
      !$acc    delete(windStressZonal, windStressMerid)

      call mpas_timer_stop("bulk_ws")

      !-----------------------------------------------------------------

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

   subroutine ocn_surface_bulk_forcing_thick(forcingPool, surfaceThicknessFlux, surfaceThicknessFluxRunoff, err)!{{{

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

      integer :: iCell, nCells

      real (kind=RKIND), dimension(:), pointer :: evaporationFlux, snowFlux
      real (kind=RKIND), dimension(:), pointer :: seaIceFreshWaterFlux, icebergFreshWaterFlux, riverRunoffFlux, iceRunoffFlux
      real (kind=RKIND), dimension(:), pointer :: rainFlux

      err = 0

      if (bulkThicknessFluxOff) return

      call mpas_timer_start("bulk_thick", .false.)

      call mpas_pool_get_array(forcingPool, 'evaporationFlux', evaporationFlux)
      call mpas_pool_get_array(forcingPool, 'snowFlux', snowFlux)
      call mpas_pool_get_array(forcingPool, 'seaIceFreshWaterFlux', seaIceFreshWaterFlux)
      call mpas_pool_get_array(forcingPool, 'icebergFreshWaterFlux', icebergFreshWaterFlux)
      call mpas_pool_get_array(forcingPool, 'riverRunoffFlux', riverRunoffFlux)
      call mpas_pool_get_array(forcingPool, 'iceRunoffFlux', iceRunoffFlux)
      call mpas_pool_get_array(forcingPool, 'rainFlux', rainFlux)

      nCells = nCellsHalo( 2 )

      ! Build surface fluxes at cell centers
#ifdef MPAS_OPENACC
      !$acc enter data copyin(evaporationFlux, snowFlux, seaIceFreshWaterFlux, icebergFreshWaterFlux, &
      !$acc     riverRunoffFlux, iceRunoffFlux, rainFlux)

      !$acc parallel loop &
      !$acc     present(surfaceThicknessFlux, surfaceThicknessFluxRunoff, evaporationFlux, snowFlux, &
      !$acc     seaIceFreshWaterFlux, icebergFreshWaterFlux, riverRunoffFlux, iceRunoffFlux, rainFlux)
#else
      !$omp parallel
      !$omp do schedule(runtime)
#endif
      do iCell = 1, nCells
        surfaceThicknessFlux(iCell) = surfaceThicknessFlux(iCell) + ( snowFlux(iCell) + rainFlux(iCell) + evaporationFlux(iCell) &
                                    + seaIceFreshWaterFlux(iCell) + icebergFreshWaterFlux(iCell) + iceRunoffFlux(iCell) ) / rho_sw
        surfaceThicknessFluxRunoff(iCell) = riverRunoffFlux(iCell) / rho_sw
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

#ifdef MPAS_OPENACC
      !$acc exit data delete(evaporationFlux, snowFlux, seaIceFreshWaterFlux, icebergFreshWaterFlux, &
      !$acc     riverRunoffFlux, iceRunoffFlux, rainFlux)
#endif

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

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< [out] error flag

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** Initialize error flag and module defaults

      err = 0
      bulkWindStressOff    = .true.
      bulkThicknessFluxOff = .true.

      !*** Set on/off flags based in user input configuration

      if (config_use_bulk_wind_stress) &
         bulkWindStressOff    = .false.
      if (config_use_bulk_thickness_flux) &
         bulkThicknessFluxOff = .false.

      !-----------------------------------------------------------------

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
      integer, pointer :: index_temperature_fluxPtr, index_salinity_fluxPtr
      integer :: index_temperature_flux, index_salinity_flux
      integer, dimension(:), pointer :: nCellsArray

      type(mpas_pool_type),pointer :: tracersSurfaceFluxPool

      real (kind=RKIND), dimension(:), pointer :: latentHeatFlux, sensibleHeatFlux, longWaveHeatFluxUp, longWaveHeatFluxDown, &
                                                  seaIceHeatFlux, icebergHeatFlux, evaporationFlux, riverRunoffFlux
      real (kind=RKIND), dimension(:), pointer :: seaIceFreshWaterFlux, icebergFreshWaterFlux, seaIceSalinityFlux, iceRunoffFlux
      real (kind=RKIND), dimension(:), pointer :: shortWaveHeatFlux, penetrativeTemperatureFlux
      real (kind=RKIND), dimension(:), pointer :: snowFlux, rainFlux
      real (kind=RKIND), dimension(:), pointer :: rainTemperatureFlux, evapTemperatureFlux, &
                                                  seaIceTemperatureFlux, icebergTemperatureFlux, &
                                                  totalFreshWaterTemperatureFlux
      real (kind=RKIND) :: requiredSalt, allowedSalt, surfaceTemperatureFluxWithoutRunoff

      err = 0

      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)

      call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceFlux',tracersSurfaceFluxPool)

      call mpas_pool_get_dimension(tracersSurfaceFluxPool, &
                                  'index_temperatureSurfaceFlux', &
                                   index_temperature_fluxPtr)
      call mpas_pool_get_dimension(tracersSurfaceFluxPool, &
                                  'index_salinitySurfaceFlux', &
                                   index_salinity_fluxPtr)
      index_temperature_flux = index_temperature_fluxPtr
      index_salinity_flux    = index_salinity_fluxPtr

      call mpas_pool_get_array(forcingPool, 'latentHeatFlux', latentHeatFlux)
      call mpas_pool_get_array(forcingPool, 'sensibleHeatFlux', sensibleHeatFlux)
      call mpas_pool_get_array(forcingPool, 'longWaveHeatFluxUp', longWaveHeatFluxUp)
      call mpas_pool_get_array(forcingPool, 'longWaveHeatFluxDown', longWaveHeatFluxDown)
      call mpas_pool_get_array(forcingPool, 'seaIceHeatFlux', seaIceHeatFlux)
      call mpas_pool_get_array(forcingPool, 'icebergHeatFlux', icebergHeatFlux)
      call mpas_pool_get_array(forcingPool, 'rainFlux', rainFlux)
      call mpas_pool_get_array(forcingPool, 'snowFlux', snowFlux)
      call mpas_pool_get_array(forcingPool, 'shortWaveHeatFlux', shortWaveHeatFlux)
      call mpas_pool_get_array(forcingPool, 'evaporationFlux', evaporationFlux)

      call mpas_pool_get_array(forcingPool, 'seaIceFreshWaterFlux', seaIceFreshWaterFlux)
      call mpas_pool_get_array(forcingPool, 'icebergFreshWaterFlux', icebergFreshWaterFlux)
      call mpas_pool_get_array(forcingPool, 'seaIceSalinityFlux', seaIceSalinityFlux)
      call mpas_pool_get_array(forcingPool, 'iceRunoffFlux', iceRunoffFlux)
      call mpas_pool_get_array(forcingPool, 'riverRunoffFlux', riverRunoffFlux)
      call mpas_pool_get_array(forcingPool, 'penetrativeTemperatureFlux', penetrativeTemperatureFlux)

      call mpas_pool_get_array(forcingPool, 'rainTemperatureFlux', rainTemperatureFlux)
      call mpas_pool_get_array(forcingPool, 'evapTemperatureFlux', evapTemperatureFlux)
      call mpas_pool_get_array(forcingPool, 'seaIceTemperatureFlux', seaIceTemperatureFlux)
      call mpas_pool_get_array(forcingPool, 'icebergTemperatureFlux', icebergTemperatureFlux)
      call mpas_pool_get_array(forcingPool, 'totalFreshWaterTemperatureFlux', totalFreshWaterTemperatureFlux)

      nCells = nCellsArray( 3 )

      ! Build surface fluxes at cell centers
      !$omp parallel
      !$omp do schedule(runtime) private(allowedSalt, requiredSalt)
      do iCell = 1, nCells
        tracersSurfaceFlux(index_temperature_flux, iCell) = tracersSurfaceFlux(index_temperature_flux, iCell) &
                                                           + (latentHeatFlux(iCell) + sensibleHeatFlux(iCell) &
                                                           + longWaveHeatFluxUp(iCell) + longWaveHeatFluxDown(iCell) &
                                                           + seaIceHeatFlux(iCell) + icebergHeatFlux(iCell) &
                                                           - (snowFlux(iCell) + iceRunoffFlux(iCell)) &
                                                           * latent_heat_fusion_mks) * hflux_factor

        ! Negative seaIceSalinityFlux is an extraction of salt from the ocean
        ! So, we negate seaIceSalinityFlux when determining how much salt this flux needs.
        requiredSalt = - seaIceSalinityFlux(iCell) * sflux_factor * dt / layerThickness(minLevelCell(iCell), iCell)
        allowedSalt = min( 4.0_RKIND, tracerGroup(index_salinity_flux, minLevelCell(iCell), iCell) )

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
      !$omp end parallel
      ! assume that snow comes in at 0 C

      ! Surface fluxes of water have an associated heat content, but the coupled system does not account for this
      ! Assume surface fluxes of water have a temperature dependent on the incoming mass flux.
      ! Assume surface fluxes of water have zero salinity. So the RHS forcing is zero for salinity.
      ! Only include this heat forcing when bulk thickness is turned on
      ! indices on tracerGroup are (iTracer, iLevel, iCell)
      if (.not. bulkThicknessFluxOff) then
         !$omp parallel
         !$omp do schedule(runtime) private(surfaceTemperatureFluxWithoutRunoff)
         do iCell = 1, nCells

           ! Accumulate fluxes that use the surface temperature
           rainTemperatureFlux(iCell) = rainFlux(iCell) * tracerGroup(index_temperature_flux,minLevelCell(iCell),iCell) / rho_sw
           evapTemperatureFlux(iCell) = evaporationFlux(iCell) * tracerGroup(index_temperature_flux,minLevelCell(iCell),iCell) / rho_sw

           ! Runoff can only have a minimum temperature of 0.0C, since it is fresh water.
           tracersSurfaceFluxRunoff(index_temperature_flux,iCell) = riverRunoffFlux(iCell) &
                      * max(tracerGroup(index_temperature_flux,minLevelCell(iCell),iCell), 0.0_RKIND) / rho_sw

           ! Accumulate fluxes that use the freezing point
! mrp performance note: should call ocn_freezing_temperature just once here 
           seaIceTemperatureFlux(iCell) = seaIceFreshWaterFlux(iCell) * &
               ocn_freezing_temperature( tracerGroup(index_salinity_flux, minLevelCell(iCell), iCell), pressure=0.0_RKIND, &
                                         inLandIceCavity=.false.) / rho_sw
           icebergTemperatureFlux(iCell) = icebergFreshWaterFlux(iCell) * &
               ocn_freezing_temperature( tracerGroup(index_salinity_flux, minLevelCell(iCell), iCell), pressure=0.0_RKIND, &
                                         inLandIceCavity=.false.) / rho_sw

           surfaceTemperatureFluxWithoutRunoff = rainTemperatureFlux(iCell) + evapTemperatureFlux(iCell) + &
                                                 seaIceTemperatureFlux(iCell) + icebergTemperatureFlux(iCell)

           tracersSurfaceFlux(index_temperature_flux, iCell) = tracersSurfaceFlux(index_temperature_flux, iCell) &
              + surfaceTemperatureFluxWithoutRunoff

           ! add runoff contribution for sending through coupler
           totalFreshWaterTemperatureFlux(iCell) = surfaceTemperatureFluxWithoutRunoff &
              + tracersSurfaceFluxRunoff(index_temperature_flux,iCell)

           ! Fields with zero temperature are not accumulated. These include:
           !    snowFlux
           !    iceRunoffFlux

         end do
         !$omp end do
         !$omp end parallel
      endif ! bulkThicknessFluxOn

      ! convert short wave heat flux to a temperature flux
      !$omp parallel
      !$omp do schedule(runtime)
      do iCell = 1, nCells
         penetrativeTemperatureFlux(iCell) = shortWaveHeatFlux(iCell) * hflux_factor
      end do
      !$omp end do
      !$omp end parallel

   end subroutine ocn_surface_bulk_forcing_active_tracers!}}}

end module ocn_surface_bulk_forcing


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
