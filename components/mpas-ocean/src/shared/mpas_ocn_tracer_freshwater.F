! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_freshwater
!
!> \brief MPAS freshwater tracer routines
!> \author Mathew Maltrud
!> \date   08/19/2021
!> \details
!>  This module contains routines for computing the tracer tendency due to freshwater fluxes
!
!-----------------------------------------------------------------------

module ocn_tracer_freshwater

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use ocn_constants
   use ocn_config

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

   public :: ocn_tracer_freshwater_surface_flux_compute, &
             ocn_tracer_freshwater_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_freshwater_surface_flux_compute
!
!> \brief   computes a tracer tendency due to CFC surface fluxes
!> \author  Mathew Maltrud
!> \date    12/22/2016
!> \details
!>  This routine computes a tracer tendency due to CFC surface fluxes
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_freshwater_surface_flux_compute(tracersPool, forcingPool, nTracers,   &
      nCellsSolve, freshwaterSurfaceFlux, freshwaterSurfaceFluxRunoff, err)!{{{


      use ocn_config

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      ! scalars
      integer, intent(in) :: nTracers, nCellsSolve

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: tracersPool, forcingPool

      ! two dimensional arrays
      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         freshwaterSurfaceFlux, freshwaterSurfaceFluxRunoff

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
      real (kind=RKIND), dimension(:), pointer :: &
         rainFlux, snowFlux, riverRunoffFlux, iceRunoffFlux, seaIceFreshWaterFlux, icebergFreshWaterFlux, landIceFreshwaterFlux

      integer :: iCell
      integer, pointer :: indexRain, indexSnow, indexRiverRunoff, indexIceRunoff, indexSeaIceFreshWater, indexIcebergFreshWater, indexLandIceFreshwater

      err = 0

      call mpas_pool_get_dimension(tracersPool, 'index_rainConcentration', indexRain)
      call mpas_pool_get_dimension(tracersPool, 'index_snowConcentration', indexSnow)
      call mpas_pool_get_dimension(tracersPool, 'index_riverRunoffConcentration', indexRiverRunoff)
      call mpas_pool_get_dimension(tracersPool, 'index_iceRunoffConcentration', indexIceRunoff)
      call mpas_pool_get_dimension(tracersPool, 'index_seaIceFreshWaterConcentration', indexSeaIceFreshWater)
      call mpas_pool_get_dimension(tracersPool, 'index_icebergFreshWaterConcentration', indexIcebergFreshWater)
      call mpas_pool_get_dimension(tracersPool, 'index_landIceFreshwaterConcentration', indexLandIceFreshwater)

      call mpas_pool_get_array(forcingPool, 'rainFlux', rainFlux)
      call mpas_pool_get_array(forcingPool, 'snowFlux', snowFlux)
      call mpas_pool_get_array(forcingPool, 'riverRunoffFlux', riverRunoffFlux)
      call mpas_pool_get_array(forcingPool, 'iceRunoffFlux', iceRunoffFlux)
      call mpas_pool_get_array(forcingPool, 'seaIceFreshWaterFlux', seaIceFreshWaterFlux)
      call mpas_pool_get_array(forcingPool, 'icebergFreshWaterFlux', icebergFreshWaterFlux)
      call mpas_pool_get_array(forcingPool, 'landIceFreshWaterFlux', landIceFreshwaterFlux)

      do iCell=1,nCellsSolve
          freshwaterSurfaceFlux(indexRain, iCell) = rainFlux(iCell) / rho_sw
          freshwaterSurfaceFlux(indexSnow, iCell) = snowFlux(iCell) / rho_sw
          freshwaterSurfaceFlux(indexIceRunoff, iCell) = iceRunoffFlux(iCell) / rho_sw
          freshwaterSurfaceFlux(indexSeaIceFreshWater, iCell) = seaIceFreshWaterFlux(iCell) / rho_sw
          freshwaterSurfaceFlux(indexIcebergFreshWater, iCell) = icebergFreshWaterFlux(iCell) / rho_sw
          freshwaterSurfaceFlux(indexLandIceFreshwater, iCell) = landIceFreshwaterFlux(iCell) / rho_sw
          freshwaterSurfaceFluxRunoff(indexRiverRunoff, iCell) = riverRunoffFlux(iCell) / rho_sw
      enddo  !  iCell

   end subroutine ocn_tracer_freshwater_surface_flux_compute!}}}

!***********************************************************************
!
!  routine ocn_tracer_freshwater_init
!
!> \brief   Initializes ocean ideal age
!> \author  Mathew Maltrud
!> \date    08/19/2021
!> \details
!>  This routine initializes fields required for tracer ideal age
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_freshwater_init(domain,err)!{{{

      type (domain_type), intent(inout) :: domain !< Input/Output: domain information

      integer, intent(out) :: err !< Output: error flag

      err = 0

      if (.not. config_use_freshwaterTracers) return

   end subroutine ocn_tracer_freshwater_init!}}}

!***********************************************************************

end module ocn_tracer_freshwater

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
