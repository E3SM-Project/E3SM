! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_sea_ice
!
!> \brief MPAS ocean sea ice formation module
!> \author Doug Jacobsen
!> \date   08/19/2013
!> \version SVN:$Id:$
!> \details
!>  This module contains routines for the formation of sea ice.
!
!-----------------------------------------------------------------------

module ocn_sea_ice

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

   public :: ocn_sea_ice_formation, &
             ocn_sea_ice_init

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
!  routine ocn_sea_ice_formation
!
!> \brief   Performs the formation of Sea Ice within the ocean.
!> \author  Doug Jacobsen
!> \date    08/19/2013
!> \version SVN:$Id$
!> \details 
!>   ocn_sea_ice_formation performs the adjustment of tracer values
!>   and layerThickness based on the formation of frazil ice within the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_sea_ice_formation(grid, indexTemperature, indexSalinity, layerThickness, tracers, seaIceEnergy, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: grid !< Input: Grid/Mesh information

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

      integer :: nCells, nVertLevels, maxLevel, nTracers, nCellsSolve
      integer :: iCell, k, iTracer

      integer, dimension(:), pointer :: maxLevelCell

      real (kind=RKIND) :: netEnergyChange, availableEnergyChange, energyChange
      real (kind=RKIND) :: temperatureChange, thicknessChange, iceThicknessChange
      real (kind=RKIND) :: referenceSalinity, iceSalinity
      real (kind=RKIND) :: freezingTemp, density_ice
      real (kind=RKIND), dimension(:), allocatable :: iceTracer

      if(.not. frazilFormationOn) return

      nCells = grid % nCells
      nCellsSolve = grid % nCellsSolve
      nVertLevels = grid % nVertLevels
      nTracers = size(tracers, dim=1)

      maxLevelCell => grid % maxLevelCell % array

      allocate(iceTracer(nTracers))
      iceTracer = 0.0_RKIND
      iceTracer(indexSalinity) = sea_ice_salinity * ppt_to_salt
      density_ice = rho_ice

      do iCell = 1, nCellsSolve ! Check performance of these two loop definitions
!     do iCell = nCellsSolve, 1, -1
         maxLevel = min(maxLevelCell(iCell), verticalLevelCap)
         netEnergyChange = 0.0_RKIND

         ! Loop over vertical levels, starting from the bottom of a column
         do k = maxLevel, 1, -1
            freezingTemp = ocn_freezing_temperature(tracers(indexSalinity, k, iCell))
            ! availableEnergyChange is:
            !     positive when frazil ice is formed
            !     negative when frazil ice can be melted
            availableEnergyChange = config_density0 * cp_sw * layerThickness(k, iCell) &
                                      * (freezingTemp - tracers(indexTemperature, k, iCell))

            ! energyChange is capped when negative.
            !      melting energy can't be greater than the amount of energy
            !      available in formed ice.
            energyChange = max(availableEnergyChange, -netEnergyChange)

            ! Compute temperature change in ocean cell due to energy change
            temperatureChange = energyChange / ( config_density0 * cp_sw * layerThickness(k, iCell) )
            ! Compute thickness change in ocean cell due to energy change
            thicknessChange = energyChange / ( config_density0 * latent_heat_fusion_mks )
            ! Compute thickness change in sea ice due to energy change
            iceThicknessChange = energyChange / ( density_ice * latent_heat_fusion_mks )

            ! Update all tracers based on the thickness change
            do iTracer = 1, nTracers
               if(iTracer /= indexTemperature) then
                  ! computed as:
                  !    \rho_{ocn} h_{ocn}^{pre} \theta_{ocn}^{pre} = 
                  !               \rho_{ocn}^{new} h_{ocn}^{new} \theta_{ocn}^{new}  =  \rho_{si} h_{si} \theta_{si}
                  tracers(iTracer, k, iCell) = ( config_density0 * layerThickness(k,iCell) * tracers(iTracer, k, iCell) &
                                               - density_ice * iceThicknessChange * iceTracer(iTracer)) / &
                                               (config_density0 * (layerThickness(k,iCell) + thicknessChange))
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
            availableEnergyChange = config_density0 * cp_sw * layerThickness(k, iCell) &
                                      * (freezingTemp - tracers(indexTemperature, k, iCell))

            ! energyChange is capped when negative.
            !      melting energy can't be greater than the amount of energy
            !      available in formed ice.
            !      compared with seaIceEnergy in this case, rather than netEnergyChange
            energyChange = max(availableEnergyChange, -seaIceEnergy(iCell))

            ! Compute temperature change in ocean cell due to energy change
            temperatureChange = energyChange / ( config_density0 * cp_sw * layerThickness(k, iCell) )
            ! Compute thickness change in ocean cell due to energy change
            thicknessChange = energyChange / ( config_density0 * latent_heat_fusion_mks )
            ! Compute thickness change in sea ice due to energy change
            iceThicknessChange = energyChange / ( density_ice * latent_heat_fusion_mks )

            ! Update all tracers based on the thickness change
            do iTracer = 1, nTracers
               if(iTracer /= indexTemperature) then
                  ! computed as:
                  !    \rho_{ocn} h_{ocn}^{pre} \theta_{ocn}^{pre} = 
                  !               \rho_{ocn}^{new} h_{ocn}^{new} \theta_{ocn}^{new}  =  \rho_{si} h_{si} \theta_{si}
                  tracers(iTracer, k, iCell) = ( config_density0 * layerThickness(k,iCell) * tracers(iTracer, k, iCell) &
                                               - density_ice * iceThicknessChange * iceTracer(iTracer)) / &
                                               (config_density0 * (layerThickness(k,iCell) + thicknessChange))
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

   end subroutine ocn_sea_ice_formation!}}}

!***********************************************************************
!
!  function ocn_freezing_temperature
!
!> \brief   Computes the freezing temperature of the ocean.
!> \author  Doug Jacobsen
!> \date    08/29/2013
!> \version SVN:$Id$
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
!  routine ocn_sea_ice_init
!
!> \brief   Initializes ocean sea ice module.
!> \author  Doug Jacobsen
!> \date    08/19/2013
!> \version SVN:$Id$
!> \details 
!>  This routine initializes the ocean sea ice module and variables..
!
!-----------------------------------------------------------------------

   subroutine ocn_sea_ice_init(nVertLevels, err)!{{{

      integer, intent(in) :: nVertLevels !< Input: Number of vertical levels suggested for level cap
      integer, intent(out) :: err !< Output: error flag

      err = 0

      frazilFormationOn = .false.

      if(config_frazil_ice_formation) then
        frazilFormationOn = .true.
      end if

      if(.not. config_monotonic) then
        verticalLevelCap = 1
      else
        verticalLevelCap = nVertLevels
      end if

   end subroutine ocn_sea_ice_init!}}}

!***********************************************************************

end module ocn_sea_ice

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
