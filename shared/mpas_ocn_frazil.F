! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_frazil
!
!> \brief MPAS ocean frazil formation module
!> \author Todd Ringler
!> \date   10/19/2015
!> \details
!>  This module contains routines for the formation of frazil ice.
!
!-----------------------------------------------------------------------

module ocn_frazil

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

   public :: ocn_frazil_formation, &
             ocn_frazil_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   integer :: verticalLevelCap
   logical :: frazilFormationOn

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_frazil_formation
!
!> \brief   Performs the formation of frazil within the ocean.
!> \author  Todd Ringler
!> \date    10/19/2015
!> \details 
!>   ocn_frazil_formation compute the tendencies to layer thickness, temperature and salinity
!>	due to the creation and possible melting of frazil ice
!>
!>   these tendencies can be retrieved at any point by calling into ocn_frazil_*_tendency routines
!>	where * is layerThickness, temperature or salinity
!>
!>   the pressure exerted by the frazil on the ocean "surface" can be retrieved by calling into
!> 	ocn_frazil_surface_pressure
!
!-----------------------------------------------------------------------

   subroutine ocn_frazil_formation(meshPool, indexTemperature, indexSalinity, layerThickness, tracers, seaIceEnergy, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information

      integer :: indexTemperature !< Input: Index in tracers array for temperature
      integer :: indexSalinity !< Input: Index in tracers array for salinity

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(inout) :: seaIceEnergy !< Input/Output: Accumulated energy for sea ice formation
      real (kind=RKIND), dimension(:,:,:), intent(inout) :: tracers !< Input/Output: Array of tracers
      real (kind=RKIND), dimension(:,:), intent(inout) :: layerThickness !< Input/Output: Thickness of each layer
      integer, intent(inout) :: err !< Error flag

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: maxLevel, nTracers
      integer :: iCell, k, iTracer
      integer, pointer :: nCells, nVertLevels, nCellsSolve

      integer, dimension(:), pointer :: maxLevelCell

      real (kind=RKIND) :: temperatureTendency, thicknessTendency, salinityTendency


      real (kind=RKIND) :: netEnergyChange, availableEnergyChange, energyChange
      real (kind=RKIND) :: temperatureTendency, thicknessTendency, salinityTendency

      real (kind=RKIND) :: referenceSalinity, iceSalinity
      real (kind=RKIND) :: freezingTemp, density_ice
      real (kind=RKIND), dimension(:), allocatable :: iceTracer

      if(.not. frazilFormationOn) return

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)


      do iCell = 1, nCellsSolve

	columnTemperatureMin = min(tracers(indexTemperature,:,iCell))
	freezingTemp = ocn_freezing_temperature(tracers(indexSalinity, 1, iCell))

	if (columnTemperatureMin < freezingTemp) then

           do k=1,nVertLevels
             deltaTemperature = freezingTemp - tracers(indexTemperature,k,iCell)
	     energyPotential = cp_sw * density(k,iCell) * deltaTemperature 
             fractionalFreezePotential(k) = max(0.0_RKIND, energyPotential / (energyPotential + density_ice * latent_heat_fusion_mks))
             fractionalMeltPotential(k) = max(0.0_RKIND, -energyPotential / energyPotential + density_ice * latent_heat_fusion_mks)
           enddo
     
           ! set all accumulators to zero
           frazilThickness = 0.0_RKIND
           frazilMass = 0.0_RKIND

           do k=maxLevelCell(iCell),1,-1
             
             ! test to see if frazil is created in this layer
             if (fractionalFreezePotential(k).gt.0.0_RKIND) then
                thicknessFreeze = min(fractionalFreezePotential(k), fractionalFrazilLimiter)*layerThickness(k,iCell)
                thickessLiquid = frazilReferenceSalinity * thicknessFreeze / (tracers(indexSalinity,k,iCell)-frazilReferenceSalinity)
                frazilThicknessTend(k,iCell) = -(thicknessFreeze+thickessLiquid) / dt
                frazilSalinityTend(k,iCell) = -thickessLiquid*tracers(indexSalinity,k,iCell) / dt
                frazilTemperatureTend(k,iCell) = latent_heat_fusion_mks * thicknessFreeze / cp_sw / dt
                frazilThickness = frazilThickness + thicknessFreeze + thickessLiquid
                frazilMass = frazilMass + (thicknessFreeze + thickessLiquid) * density(k,iCell)
             else 
                if (frazilThickness.gt.0.0_RKIND) then
                   thicknessMelt = min(frazilThickness, fractionalMeltPotential(k)*layerThickness(k,iCell), fractionalFrazilLimiter*layerThickness(k,iCell))
                   frazilThicknessTend(k,iCell) = thicknessMelt / dt
                   frazilSalinityTend(k,iCell) = thicknessMelt * frazilReferenceSalinity / dt
                   frazilTemperatureTend(k,iCell) = -latent_heat_fusion_mks * thicknessFreeze / cp_sw / dt
                   frazilThickness = frazilThickness - thicknessMelt
                   frazilMass = frazilMass - thicknessMelt * density_ice
                endif
             endif  ! (fractionalFreezePotential(k).gt.0.0_RKIND)

          enddo  ! do k=maxLevelCell(iCell),1,-1

          frazilSurfacePressure(newTime) = frazilSurfacePressure(oldTime) + frazilMass*gravity

       endif

    enddo

           



	endif  ! (columnTemperatureMin < freezingTemp)

	



         maxLevel = min(maxLevelCell(iCell), verticalLevelCap)
         netEnergyChange = 0.0_RKIND

         ! Loop over vertical levels, starting from the bottom of a column
         do k = maxLevel, 1, -1
            freezingTemp = ocn_freezing_temperature(tracers(indexSalinity, k, iCell))
            ! availableEnergyChange is:
            !     positive when frazil ice is formed
            !     negative when frazil ice can be melted
            availableEnergyChange = rho_sw * cp_sw * layerThickness(k, iCell) &
                                      * (freezingTemp - tracers(indexTemperature, k, iCell))

            ! energyChange is capped when negative.
            !      melting energy can't be greater than the amount of energy
            !      available in formed ice.
            energyChange = max(availableEnergyChange, -netEnergyChange)

            ! Compute temperature change in ocean cell due to energy change
            temperatureChange = energyChange / ( rho_sw * cp_sw * layerThickness(k, iCell) )
            ! Compute thickness change in ocean cell due to energy change
            thicknessChange = energyChange / ( rho_sw * latent_heat_fusion_mks )
            ! Compute thickness change in sea ice due to energy change
            iceThicknessChange = energyChange / ( density_ice * latent_heat_fusion_mks )

            ! Update all tracers based on the thickness change
            do iTracer = 1, nTracers
               if(iTracer /= indexTemperature) then
                  ! computed as:
                  !    \rho_{ocn} h_{ocn}^{pre} \theta_{ocn}^{pre} = 
                  !               \rho_{ocn}^{new} h_{ocn}^{new} \theta_{ocn}^{new}  =  \rho_{si} h_{si} \theta_{si}
                  tracers(iTracer, k, iCell) = ( rho_sw * layerThickness(k,iCell) * tracers(iTracer, k, iCell) &
                                               - density_ice * iceThicknessChange * iceTracer(iTracer)) / &
                                               (rho_sw * (layerThickness(k,iCell) + thicknessChange))
               end if
            end do

            ! Adjust Temperature
            tracers(indexTemperature, k, iCell) = tracers(indexTemperature, k, iCell) + temperatureChange
            ! Adjust Thickness
            layerThickness(k,iCell) = layerThickness(k,iCell) + thicknessChange

            ! Add energyChange to netEnergyChange.
            ! netEnergyChange should always be >= 0.0
            netEnergyChange = netEnergychange + energyChange
         end do

         ! Add netEnergyChange to the cell's energy.
         ! seaIceEnergy should always be >= 0.0
         seaIceEnergy(iCell) = seaIceEnergy(iCell) + netEnergyChange

         ! Adjust top layer one more time, based on energy availabe in seaIceEnergy(iCell)
         ! This really only allows melting of previously formed ice to occur.
         if(maxLevelCell(iCell) >= 1 .and. seaIceEnergy(iCell) > 0.0_RKIND) then
            k = 1

            netEnergychange = 0.0_RKIND
            freezingTemp = ocn_freezing_temperature(tracers(indexSalinity, k, iCell))
            ! availableEnergyChange is:
            !     positive when frazil ice is formed
            !     negative when frazil ice can be melted
            availableEnergyChange = rho_sw * cp_sw * layerThickness(k, iCell) &
                                      * (freezingTemp - tracers(indexTemperature, k, iCell))

            ! energyChange is capped when negative.
            !      melting energy can't be greater than the amount of energy
            !      available in formed ice.
            !      compared with seaIceEnergy in this case, rather than netEnergyChange
            energyChange = max(availableEnergyChange, -seaIceEnergy(iCell))

            ! Compute temperature change in ocean cell due to energy change
            temperatureChange = energyChange / ( rho_sw * cp_sw * layerThickness(k, iCell) )
            ! Compute thickness change in ocean cell due to energy change
            thicknessChange = energyChange / ( rho_sw * latent_heat_fusion_mks )
            ! Compute thickness change in sea ice due to energy change
            iceThicknessChange = energyChange / ( density_ice * latent_heat_fusion_mks )

            ! Update all tracers based on the thickness change
            do iTracer = 1, nTracers
               if(iTracer /= indexTemperature) then
                  ! computed as:
                  !    \rho_{ocn} h_{ocn}^{pre} \theta_{ocn}^{pre} = 
                  !               \rho_{ocn}^{new} h_{ocn}^{new} \theta_{ocn}^{new}  =  \rho_{si} h_{si} \theta_{si}
                  tracers(iTracer, k, iCell) = ( rho_sw * layerThickness(k,iCell) * tracers(iTracer, k, iCell) &
                                               - density_ice * iceThicknessChange * iceTracer(iTracer)) / &
                                               (rho_sw * (layerThickness(k,iCell) + thicknessChange))
               end if
            end do

            ! Adjust Temperature
            tracers(indexTemperature, k, iCell) = tracers(indexTemperature, k, iCell) + temperatureChange
            ! Adjust Thickness
            layerThickness(k,iCell) = layerThickness(k,iCell) + thicknessChange

            ! Add energyChange to netEnergyChange.
            ! netEnergyChange should always be >= 0.0
            seaIceEnergy(iCell) = seaIceEnergy(iCell) + energyChange
         end if
      end do

      deallocate(iceTracer)

   end subroutine ocn_frazil_formation!}}}

!***********************************************************************
!
!  function ocn_freezing_temperature
!
!> \brief   Computes the freezing temperature of the ocean.
!> \author  Todd Ringler
!> \date    10/29/2015
!> \details 
!>  This routine computes the freezing temperature of the ocean at a given
!>  salinity value.
!
!-----------------------------------------------------------------------
    real (kind=RKIND) function ocn_freezing_temperature(salinity)!{{{
      real (kind=RKIND) :: salinity !< Input: Salinity value of water for freezing temperature
      ocn_freezing_temperature = -1.8
    end function ocn_freezing_temperature!}}}


!***********************************************************************
!
!  routine ocn_frazil_init
!
!> \brief   Initializes ocean frazil ice module.
!> \author  Todd Ringler
!> \date    10/19/2015
!> \details 
!>  This routine initializes the ocean frazil ice module and variables..
!
!-----------------------------------------------------------------------

   subroutine ocn_frazil_init(nVertLevels, err)!{{{

      integer, intent(in) :: nVertLevels !< Input: Number of vertical levels suggested for level cap
      integer, intent(out) :: err !< Output: error flag

      logical, pointer :: config_frazil_ice_formation, config_monotonic

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_frazil_ice_formation', config_frazil_ice_formation)
      call mpas_pool_get_config(ocnConfigs, 'config_monotonic', config_monotonic)

      frazilFormationOn = .false.

      if(config_frazil_ice_formation) then
        frazilFormationOn = .true.
      end if

      if(.not. config_monotonic) then
        verticalLevelCap = 1
      else
        verticalLevelCap = nVertLevels
      end if

   end subroutine ocn_frazil_init!}}}

!***********************************************************************

end module ocn_sea_ice

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
