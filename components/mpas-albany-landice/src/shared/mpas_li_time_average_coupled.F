!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! \file mpas_li_time_average_coupled.F
!
! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_time_average_coupled
!
!> \brief MALI time averager for coupling 
!> \author Stephen Price and Matthew Hoffman, modified 
!>         after similarly named subroutine for MPAS Ocean 
!>         written by Doug Jacobsen.
!> \date   04 March 2025
!> \details
!> This module contains subroutines for time averaging of MALI fluxes 
!> for use in coupling to E3SM.
!
!-------------------------------------------------------------------------------

module li_time_average_coupled

    use mpas_kind_types
    use mpas_derived_types
    use mpas_pool_routines

    implicit none
    save
    public

    contains

!***********************************************************************
!
!  routine li_time_average_coupled_init
!
!> \brief   Coupled time averager initialization
!> \author  Stephen Price 
!> \date    04 March 2025
!> \details
!>  This routine initializes the coupled time averaging fields
!
!-----------------------------------------------------------------------
!    subroutine li_time_average_coupled_init(meshPool,timeAveragingPool)
    subroutine li_time_average_coupled_init()
            
!        type (mpas_pool_type), intent(inout) :: timeAveragingPool
!
!        real (kind=RKIND), dimension(:), pointer :: avgBareIceAblation, avgCalvingThickness !, &
!
!        integer :: iCell
!        integer, pointer :: nAccumulatedCoupled, nCells
!
!        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
!
!        call mpas_pool_get_array(timeAveragingPool, 'nAccumulatedCoupled', nAccumulatedCoupled)
!        call mpas_pool_get_array(timeAveragingPool, 'avgBareIceAblation', avgBareIceAbltion)
!        call mpas_pool_get_array(timeAveragingPool, 'avgCalvingThickness', avgCalvingThickness)
!
!        do iCell = 1, nCells
!           avgBareIceAblation(iCell) = 0.0_RKIND
!           avgCalvingThickness(iCell) = 0.0_RKIND
!        end do
!
!        nAccumulatedCoupled = 0

    end subroutine li_time_average_coupled_init

!***********************************************************************
!
!  routine li_time_average_coupled_accumulate
!
!> \brief   Coupled time averager accumulation
!> \author  Stephen Price
!> \date    04 March 2025
!> \details
!>  This routine accumulated the coupled time averaging fields
!
!-----------------------------------------------------------------------
!    subroutine li_time_average_coupled_accumulate(statePool, forcingPool, timeLevel)
!    subroutine li_time_average_coupled_accumulated(timeAveragingPool, geometryPool, meshPool)
    subroutine li_time_average_coupled_accumulate()
!        use li_constants, only: &
!             constants-list-here   

!        type (mpas_pool_type), intent(in) :: statePool
!        type (mpas_pool_type), intent(inout) :: forcingPool
!        integer, intent(in) :: timeLevel

!        real (kind=RKIND), dimension(:,:), pointer :: avgSurfaceVelocity
!        real (kind=RKIND), dimension(:,:), pointer :: avgTracersSurfaceValue
!        real (kind=RKIND), dimension(:,:), pointer :: avgSSHGradient
!        integer :: iCell
!        integer, pointer :: index_temperaturePtr, index_SSHzonalPtr, &
!                 index_SSHmeridionalPtr, nAccumulatedCoupled, nCells, nVertLevels
!        integer :: index_temperature, index_SSHzonal, index_SSHmeridional
!        real (kind=RKIND), dimension(:,:), pointer :: &
!                                                      avgLandIceBoundaryLayerTracers, avgLandIceTracerTransferVelocities
!        real (kind=RKIND), dimension(:), pointer :: effectiveDensityInLandIce, avgEffectiveDensityInLandIce, &
!                                                    totalFreshWaterTemperatureFlux, avgTotalFreshWaterTemperatureFlux, &
!                                                    landIceFreshwaterFlux, avgLandIceFreshwaterFlux, &
!                                                    landIceHeatFlux, avgLandIceHeatFlux, &
!                                                    removedRiverRunoffFlux, avgRemovedRiverRunoffFlux, &
!                                                    removedIceRunoffFlux, avgRemovedIceRunoffFlux, &
!                                                    avgRemovedIceRunoffHeatFlux, &
!                                                    avgThermalForcingAtCritDepth
!
!        type (mpas_pool_type), pointer :: tracersPool
!
!        real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers
!        integer, pointer :: indexTemperature, indexSalinity
!        integer :: iLevelCritDepth, iLevel
!        real (kind=RKIND) :: freezingTemp
!
!        real (kind=RKIND), dimension(:,:,:), pointer :: &
!           ecosysTracers,  &
!           DMSTracers,     &
!           MacroMoleculesTracers
!
!        type (mpas_pool_type), pointer :: ecosysSeaIceCoupling, &
!                                          DMSSeaIceCoupling,    &
!                                          MacroMoleculesSeaIceCoupling
!
!!        real (kind=RKIND), dimension(:), pointer :: CO2_gas_flux, avgCO2_gas_flux
!
!!        call mpas_pool_get_array(forcingPool, 'avgTracersSurfaceValue', avgTracersSurfaceValue)
!!        call mpas_pool_get_array(forcingPool, 'avgSurfaceVelocity', avgSurfaceVelocity)
!!        call mpas_pool_get_array(forcingPool, 'avgSSHGradient', avgSSHGradient)
!
!        call mpas_pool_get_array(forcingPool, 'totalFreshWaterTemperatureFlux', totalFreshWaterTemperatureFlux)
!        call mpas_pool_get_array(forcingPool, 'avgTotalFreshWaterTemperatureFlux', avgTotalFreshWaterTemperatureFlux)
!
!        call mpas_pool_get_dimension(forcingPool, 'nCells', nCells)
!        call mpas_pool_get_dimension(forcingPool, 'nVertLevels', nVertLevels)
!        call mpas_pool_get_dimension(forcingPool, &
!                                     'index_avgTemperatureSurfaceValue', &
!                                      index_temperaturePtr)
!        call mpas_pool_get_dimension(forcingPool, &
!                                     'index_avgSSHGradientZonal', &
!                                      index_SSHzonalPtr)
!        call mpas_pool_get_dimension(forcingPool, &
!                                     'index_avgSSHGradientMeridional', &
!                                      index_SSHmeridionalPtr)
!        index_temperature   = index_temperaturePtr
!        index_SSHzonal      = index_SSHzonalPtr
!        index_SSHmeridional = index_SSHmeridionalPtr
!
!        call mpas_pool_get_array(forcingPool, 'nAccumulatedCoupled', nAccumulatedCoupled)
!
!        do iCell = 1, nCells
!
!!           avgTracersSurfaceValue(:, iCell) = avgTracersSurfaceValue(:, iCell) * nAccumulatedCoupled &
!!                                            + tracersSurfaceValue(:, iCell)
!!           avgTracersSurfaceValue(index_temperature, iCell) = avgTracersSurfaceValue(index_temperature, iCell) + T0_Kelvin
!!           avgTracersSurfaceValue(:, iCell) = avgTracersSurfaceValue(:, iCell) / ( nAccumulatedCoupled + 1 )
!!
!!           avgSSHGradient(index_SSHzonal, iCell) = ( avgSSHGradient(index_SSHzonal, iCell) * nAccumulatedCoupled &
!!                                                 + gradSSHZonal(iCell) ) / ( nAccumulatedCoupled + 1 )
!!           avgSSHGradient(index_SSHmeridional, iCell) = ( avgSSHGradient(index_SSHmeridional, iCell) * nAccumulatedCoupled &
!!                                                 + gradSSHMeridional(iCell) ) / ( nAccumulatedCoupled + 1 )
!!           avgSurfaceVelocity(:, iCell) = ( avgSurfaceVelocity(:, iCell) * nAccumulatedCoupled + surfaceVelocity(:, iCell) ) &
!!                                        / ( nAccumulatedCoupled + 1 )
!!           avgTotalFreshWaterTemperatureFlux(iCell) = ( avgTotalFreshWaterTemperatureFlux(iCell) * nAccumulatedCoupled &
!                                               + totalFreshWaterTemperatureFlux(iCell) ) / ( nAccumulatedCoupled + 1 )
!
!        end do
!
!!           do iCell = 1, nCells
!!              avgLandIceFreshwaterFlux(iCell) = ( avgLandIceFreshwaterFlux(iCell) * nAccumulatedCoupled &
!!                                              + landIceFreshwaterFlux(iCell) ) / ( nAccumulatedCoupled + 1)
!!              avgLandIceHeatFlux(iCell) = ( avgLandIceHeatFlux(iCell) * nAccumulatedCoupled &
!!                                        + landIceHeatFlux(iCell) ) / ( nAccumulatedCoupled + 1)
!!           end do
!
!        nAccumulatedCoupled = nAccumulatedCoupled + 1

    end subroutine li_time_average_coupled_accumulate

end module li_time_average_coupled
