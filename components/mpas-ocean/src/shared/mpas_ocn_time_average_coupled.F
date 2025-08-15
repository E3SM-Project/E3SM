! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_average_coupled
!
!> \brief MPAS ocean coupled time averager
!> \author Doug Jacobsen
!> \date   06/08/2013
!> \details
!>  This module contains the routines for time averaging
!>  coupling fields for the ocean core.
!
!-----------------------------------------------------------------------

module ocn_time_average_coupled

    use mpas_kind_types
    use mpas_derived_types
    use mpas_pool_routines
    use ocn_constants
    use ocn_config
    use ocn_equation_of_state
    use ocn_mesh
    use ocn_tracer_ecosys
    use ocn_tracer_DMS
    use ocn_tracer_MacroMolecules
    use ocn_diagnostics_variables

    implicit none
    save
    public

    contains

!***********************************************************************
!
!  routine ocn_time_average_coupled_init
!
!> \brief   Coupled time averager initialization
!> \author  Doug Jacobsen
!> \date    06/08/2013
!> \details
!>  This routine initializes the coupled time averaging fields
!
!-----------------------------------------------------------------------
    subroutine ocn_time_average_coupled_init(forcingPool)!{{{
        type (mpas_pool_type), intent(inout) :: forcingPool

        real (kind=RKIND), dimension(:,:), pointer :: avgTracersSurfaceValue, avgSurfaceVelocity, avgSSHGradient, &
                                                      avgLandIceBoundaryLayerTracers, avgLandIceTracerTransferVelocities

        real (kind=RKIND), dimension(:), pointer :: avgEffectiveDensityInLandIce, avgTotalFreshWaterTemperatureFlux, &
                                                    avgLandIceFreshwaterFlux, &
                                                    avgRemovedRiverRunoffFlux, avgRemovedIceRunoffFlux, &
                                                    avgLandIceHeatFlux, avgRemovedIceRunoffHeatFlux

        real (kind=RKIND), dimension(:,:), pointer :: avgThermalForcingAtZLevels
        real (kind=RKIND), dimension(:,:), pointer :: avgThermalForcingAtZLevelsMask

        integer :: iCell
        integer, pointer :: nAccumulatedCoupled, nCells

        type (mpas_pool_type), pointer :: ecosysAuxiliary

        type (mpas_pool_type), pointer :: ecosysSeaIceCoupling, &
                                          DMSSeaIceCoupling,    &
                                          MacroMoleculesSeaIceCoupling

        real (kind=RKIND), dimension(:), pointer :: avgCO2_gas_flux

        real (kind=RKIND), dimension(:), pointer :: avgOceanSurfaceDIC, &
                                                    avgOceanSurfaceDON, &
                                                    avgOceanSurfaceNO3, &
                                                    avgOceanSurfaceSiO3, &
                                                    avgOceanSurfaceNH4, &
                                                    avgOceanSurfaceDMS, &
                                                    avgOceanSurfaceDMSP, &
                                                    avgOceanSurfaceDOCr, &
                                                    avgOceanSurfaceDOCSemiLabile, &
                                                    avgOceanSurfaceFeParticulate, &
                                                    avgOceanSurfaceFeDissolved

        real (kind=RKIND), dimension(:,:), pointer :: avgOceanSurfacePhytoC, &
                                                      avgOceanSurfaceDOC

        call mpas_pool_get_dimension(forcingPool, 'nCells', nCells)

        call mpas_pool_get_array(forcingPool, 'avgTracersSurfaceValue', avgTracersSurfaceValue)
        call mpas_pool_get_array(forcingPool, 'avgSurfaceVelocity', avgSurfaceVelocity)
        call mpas_pool_get_array(forcingPool, 'avgSSHGradient', avgSSHGradient)
        call mpas_pool_get_array(forcingPool, 'nAccumulatedCoupled', nAccumulatedCoupled)
        call mpas_pool_get_array(forcingPool, 'avgTotalFreshWaterTemperatureFlux', avgTotalFreshWaterTemperatureFlux)

        !$omp parallel
        !$omp do schedule(runtime)
        do iCell = 1, nCells
           avgSurfaceVelocity(:, iCell) = 0.0_RKIND
           avgTracersSurfaceValue(:, iCell) = 0.0_RKIND
           avgSSHGradient(:, iCell) = 0.0_RKIND
           avgTotalFreshWaterTemperatureFlux(iCell) = 0.0_RKIND
        end do
        !$omp end do
        !$omp end parallel

        if(trim(config_land_ice_flux_mode) == 'coupled') then
           call mpas_pool_get_array(forcingPool, 'avgLandIceBoundaryLayerTracers', avgLandIceBoundaryLayerTracers)
           call mpas_pool_get_array(forcingPool, 'avgLandIceTracerTransferVelocities', avgLandIceTracerTransferVelocities)
           call mpas_pool_get_array(forcingPool, 'avgEffectiveDensityInLandIce', avgEffectiveDensityInLandIce)

           !$omp parallel
           !$omp do schedule(runtime)
           do iCell = 1, nCells
              avgLandIceBoundaryLayerTracers(:, iCell) = 0.0_RKIND
              avgLandIceTracerTransferVelocities(:, iCell) = 0.0_RKIND
              avgEffectiveDensityInLandIce(iCell) = 0.0_RKIND
           end do
           !$omp end do
           !$omp end parallel
        end if

        ! Set up polar fields if necessary
        if(trim(config_land_ice_flux_mode)=='standalone' .or. trim(config_land_ice_flux_mode) == 'data') then
           call mpas_pool_get_array(forcingPool, 'avgLandIceFreshwaterFlux', avgLandIceFreshwaterFlux)
           call mpas_pool_get_array(forcingPool, 'avgLandIceHeatFlux', avgLandIceHeatFlux)

           !$omp parallel
           !$omp do schedule(runtime)
           do iCell = 1, nCells
              avgLandIceFreshwaterFlux(iCell) = 0.0_RKIND
              avgLandIceHeatFlux(iCell) = 0.0_RKIND
           end do
           !$omp end do
           !$omp end parallel
        end if

        if(config_remove_ais_river_runoff) then
           call mpas_pool_get_array(forcingPool, 'avgRemovedRiverRunoffFlux', avgRemovedRiverRunoffFlux)

           !$omp parallel
           !$omp do schedule(runtime)
           do iCell = 1, nCells
              avgRemovedRiverRunoffFlux(iCell) = 0.0_RKIND
           end do
           !$omp end do
           !$omp end parallel
        end if

        if(config_remove_ais_ice_runoff) then
           call mpas_pool_get_array(forcingPool, 'avgRemovedIceRunoffFlux', avgRemovedIceRunoffFlux)
           call mpas_pool_get_array(forcingPool, 'avgRemovedIceRunoffHeatFlux', avgRemovedIceRunoffHeatFlux)

           !$omp parallel
           !$omp do schedule(runtime)
           do iCell = 1, nCells
              avgRemovedIceRunoffFlux(iCell) = 0.0_RKIND
              avgRemovedIceRunoffHeatFlux(iCell) = 0.0_RKIND
           end do
           !$omp end do
           !$omp end parallel
        end if

        if (trim(config_glc_thermal_forcing_coupling_mode) == '3d') then
           call mpas_pool_get_array(forcingPool, 'avgThermalForcingAtZLevels', avgThermalForcingAtZLevels)
           call mpas_pool_get_array(forcingPool, 'avgThermalForcingAtZLevelsMask', avgThermalForcingAtZLevelsMask)

           !$omp parallel
           !$omp do schedule(runtime)
           do iCell = 1, nCells
              avgThermalForcingAtZLevels(:, iCell) = 0.0_RKIND
              avgThermalForcingAtZLevelsMask(:, iCell) = 0.0_RKIND
           end do
           !$omp end do
           !$omp end parallel
        endif

        !  set up BGC coupling fields if necessary
        if (config_use_ecosysTracers) then

           call mpas_pool_get_subpool(forcingPool, 'ecosysAuxiliary', ecosysAuxiliary)
           call mpas_pool_get_array(ecosysAuxiliary, 'avgCO2_gas_flux', avgCO2_gas_flux)

           call mpas_pool_get_subpool(forcingPool, 'ecosysSeaIceCoupling', ecosysSeaIceCoupling)

           call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfacePhytoC', avgOceanSurfacePhytoC)
           call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceDIC', avgOceanSurfaceDIC)
           call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceNO3', avgOceanSurfaceNO3)
           call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceSiO3', avgOceanSurfaceSiO3)
           call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceNH4', avgOceanSurfaceNH4)
           call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceDOCr', avgOceanSurfaceDOCr)
           call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceDOCSemiLabile', avgOceanSurfaceDOCSemiLabile)
           call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceFeParticulate', avgOceanSurfaceFeParticulate)
           call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceFeDissolved', avgOceanSurfaceFeDissolved)

           !$omp parallel
           !$omp do schedule(runtime)
           do iCell = 1, nCells
              avgCO2_gas_flux(iCell) = 0.0_RKIND

              avgOceanSurfacePhytoC(:,iCell) = 0.0_RKIND

              avgOceanSurfaceDIC(iCell) = 0.0_RKIND
              avgOceanSurfaceDOCSemiLabile(iCell) = 0.0_RKIND
              avgOceanSurfaceNO3(iCell) = 0.0_RKIND
              avgOceanSurfaceSiO3(iCell) = 0.0_RKIND
              avgOceanSurfaceNH4(iCell) = 0.0_RKIND
              avgOceanSurfaceDOCr(iCell) = 0.0_RKIND
              avgOceanSurfaceFeParticulate(iCell) = 0.0_RKIND
              avgOceanSurfaceFeDissolved(iCell) = 0.0_RKIND
           end do
           !$omp end do
           !$omp end parallel
        end if

        if (config_use_DMSTracers) then
           call mpas_pool_get_subpool(forcingPool, 'DMSSeaIceCoupling', DMSSeaIceCoupling)

           call mpas_pool_get_array(DMSSeaIceCoupling, 'avgOceanSurfaceDMS', avgOceanSurfaceDMS)
           call mpas_pool_get_array(DMSSeaIceCoupling, 'avgOceanSurfaceDMSP', avgOceanSurfaceDMSP)

           !$omp parallel
           !$omp do schedule(runtime)
           do iCell = 1, nCells
              avgOceanSurfaceDMS(iCell) = 0.0_RKIND
              avgOceanSurfaceDMSP(iCell) = 0.0_RKIND
           end do
           !$omp end do
           !$omp end parallel
        endif
        if (config_use_MacroMoleculesTracers) then
           call mpas_pool_get_subpool(forcingPool, 'MacroMoleculesSeaIceCoupling', MacroMoleculesSeaIceCoupling)

           call mpas_pool_get_array(MacroMoleculesSeaIceCoupling, 'avgOceanSurfaceDOC', avgOceanSurfaceDOC)
           call mpas_pool_get_array(MacroMoleculesSeaIceCoupling, 'avgOceanSurfaceDON', avgOceanSurfaceDON)

           !$omp parallel
           !$omp do schedule(runtime)
           do iCell = 1, nCells
              avgOceanSurfaceDOC(:,iCell) = 0.0_RKIND
              avgOceanSurfaceDON(iCell) = 0.0_RKIND
           end do
           !$omp end do
           !$omp end parallel
        endif

        nAccumulatedCoupled = 0

    end subroutine ocn_time_average_coupled_init!}}}

!***********************************************************************
!
!  routine ocn_time_average_coupled_accumulate
!
!> \brief   Coupled time averager accumulation
!> \author  Doug Jacobsen
!> \date    06/08/2013
!> \details
!>  This routine accumulated the coupled time averaging fields
!
!-----------------------------------------------------------------------
    subroutine ocn_time_average_coupled_accumulate(statePool, forcingPool, timeLevel)!{{{
        use ocn_constants, only: &
             latent_heat_fusion_mks

        type (mpas_pool_type), intent(in) :: statePool
        type (mpas_pool_type), intent(inout) :: forcingPool
        integer, intent(in) :: timeLevel

        real (kind=RKIND), dimension(:,:), pointer :: avgSurfaceVelocity
        real (kind=RKIND), dimension(:,:), pointer :: avgTracersSurfaceValue
        real (kind=RKIND), dimension(:,:), pointer :: avgSSHGradient
        integer :: iCell
        integer, pointer :: index_temperaturePtr, index_SSHzonalPtr, &
                 index_SSHmeridionalPtr, nAccumulatedCoupled, nCells
        integer :: index_temperature, index_SSHzonal, index_SSHmeridional
        real (kind=RKIND), dimension(:,:), pointer :: &
                                                      avgLandIceBoundaryLayerTracers, avgLandIceTracerTransferVelocities
        real (kind=RKIND), dimension(:), pointer :: effectiveDensityInLandIce, avgEffectiveDensityInLandIce, &
                                                    totalFreshWaterTemperatureFlux, avgTotalFreshWaterTemperatureFlux, &
                                                    landIceFreshwaterFlux, avgLandIceFreshwaterFlux, &
                                                    landIceHeatFlux, avgLandIceHeatFlux, &
                                                    removedRiverRunoffFlux, avgRemovedRiverRunoffFlux, &
                                                    removedIceRunoffFlux, avgRemovedIceRunoffFlux, &
                                                    avgRemovedIceRunoffHeatFlux

        type (mpas_pool_type), pointer :: tracersPool

        real (kind=RKIND), dimension(:,:,:), pointer :: &
           ecosysTracers,  &
           DMSTracers,     &
           MacroMoleculesTracers

        type (mpas_pool_type), pointer :: ecosysSeaIceCoupling, &
                                          DMSSeaIceCoupling,    &
                                          MacroMoleculesSeaIceCoupling

        type (mpas_pool_type), pointer :: ecosysAuxiliary
        real (kind=RKIND), dimension(:), pointer :: CO2_gas_flux, avgCO2_gas_flux

        real (kind=RKIND), dimension(:), pointer :: avgOceanSurfaceDIC, &
                                                    avgOceanSurfaceDON, &
                                                    avgOceanSurfaceDOCSemiLabile, &
                                                    avgOceanSurfaceNO3, &
                                                    avgOceanSurfaceSiO3, &
                                                    avgOceanSurfaceNH4, &
                                                    avgOceanSurfaceDMS, &
                                                    avgOceanSurfaceDMSP, &
                                                    avgOceanSurfaceDOCr, &
                                                    avgOceanSurfaceFeParticulate, &
                                                    avgOceanSurfaceFeDissolved

        real (kind=RKIND), dimension(:,:), pointer :: avgOceanSurfacePhytoC, &
                                                      avgOceanSurfaceDOC

        call mpas_pool_get_array(forcingPool, 'avgTracersSurfaceValue', avgTracersSurfaceValue)
        call mpas_pool_get_array(forcingPool, 'avgSurfaceVelocity', avgSurfaceVelocity)
        call mpas_pool_get_array(forcingPool, 'avgSSHGradient', avgSSHGradient)

        call mpas_pool_get_array(forcingPool, 'totalFreshWaterTemperatureFlux', totalFreshWaterTemperatureFlux)
        call mpas_pool_get_array(forcingPool, 'avgTotalFreshWaterTemperatureFlux', avgTotalFreshWaterTemperatureFlux)

        call mpas_pool_get_dimension(forcingPool, 'nCells', nCells)
        call mpas_pool_get_dimension(forcingPool, &
                                     'index_avgTemperatureSurfaceValue', &
                                      index_temperaturePtr)
        call mpas_pool_get_dimension(forcingPool, &
                                     'index_avgSSHGradientZonal', &
                                      index_SSHzonalPtr)
        call mpas_pool_get_dimension(forcingPool, &
                                     'index_avgSSHGradientMeridional', &
                                      index_SSHmeridionalPtr)
        index_temperature   = index_temperaturePtr
        index_SSHzonal      = index_SSHzonalPtr
        index_SSHmeridional = index_SSHmeridionalPtr

        call mpas_pool_get_array(forcingPool, 'nAccumulatedCoupled', nAccumulatedCoupled)

        !$omp parallel
        !$omp do schedule(runtime)
        do iCell = 1, nCells
           avgTracersSurfaceValue(:, iCell) = avgTracersSurfaceValue(:, iCell) * nAccumulatedCoupled &
                                            + tracersSurfaceValue(:, iCell)
           avgTracersSurfaceValue(index_temperature, iCell) = avgTracersSurfaceValue(index_temperature, iCell) + T0_Kelvin
           avgTracersSurfaceValue(:, iCell) = avgTracersSurfaceValue(:, iCell) / ( nAccumulatedCoupled + 1 )

           avgSSHGradient(index_SSHzonal, iCell) = ( avgSSHGradient(index_SSHzonal, iCell) * nAccumulatedCoupled &
                                                 + gradSSHZonal(iCell) ) / ( nAccumulatedCoupled + 1 )
           avgSSHGradient(index_SSHmeridional, iCell) = ( avgSSHGradient(index_SSHmeridional, iCell) * nAccumulatedCoupled &
                                                 + gradSSHMeridional(iCell) ) / ( nAccumulatedCoupled + 1 )
           avgSurfaceVelocity(:, iCell) = ( avgSurfaceVelocity(:, iCell) * nAccumulatedCoupled + surfaceVelocity(:, iCell) ) &
                                        / ( nAccumulatedCoupled + 1 )
           avgTotalFreshWaterTemperatureFlux(iCell) = ( avgTotalFreshWaterTemperatureFlux(iCell) * nAccumulatedCoupled &
                                               + totalFreshWaterTemperatureFlux(iCell) ) / ( nAccumulatedCoupled + 1 )

        end do
        !$omp end do
        !$omp end parallel

        if(trim(config_land_ice_flux_mode) == 'coupled') then
           call mpas_pool_get_array(statePool, 'effectiveDensityInLandIce', effectiveDensityInLandIce, timeLevel)

           call mpas_pool_get_array(forcingPool, 'avgLandIceBoundaryLayerTracers', avgLandIceBoundaryLayerTracers)
           call mpas_pool_get_array(forcingPool, 'avgLandIceTracerTransferVelocities', avgLandIceTracerTransferVelocities)
           call mpas_pool_get_array(forcingPool, 'avgEffectiveDensityInLandIce', avgEffectiveDensityInLandIce)


           !$omp parallel
           !$omp do schedule(runtime)
           do iCell = 1, nCells
              avgLandIceBoundaryLayerTracers(:, iCell) = ( avgLandIceBoundaryLayerTracers(:, iCell) * nAccumulatedCoupled &
                 + landIceBoundaryLayerTracers(:, iCell) ) / ( nAccumulatedCoupled + 1 )
              avgLandIceTracerTransferVelocities(:, iCell) = ( avgLandIceTracerTransferVelocities(:, iCell) * nAccumulatedCoupled &
                 + landIceTracerTransferVelocities(:, iCell) ) / ( nAccumulatedCoupled + 1)
              avgEffectiveDensityInLandIce(iCell) = ( avgEffectiveDensityInLandIce(iCell) * nAccumulatedCoupled &
                 + effectiveDensityInLandIce(iCell) ) / ( nAccumulatedCoupled + 1)
           end do
           !$omp end do
           !$omp end parallel
        end if

        ! Accumulate polar fields if necessary
        if(trim(config_land_ice_flux_mode) == 'standalone' .or. trim(config_land_ice_flux_mode) == 'data') then
           call mpas_pool_get_array(forcingPool, 'avgLandIceFreshwaterFlux', avgLandIceFreshwaterFlux)
           call mpas_pool_get_array(forcingPool, 'landIceFreshwaterFlux', landIceFreshwaterFlux)
           call mpas_pool_get_array(forcingPool, 'avgLandIceHeatFlux', avgLandIceHeatFlux)
           call mpas_pool_get_array(forcingPool, 'landIceHeatFlux', landIceHeatFlux)

           !$omp parallel
           !$omp do schedule(runtime)
           do iCell = 1, nCells
              avgLandIceFreshwaterFlux(iCell) = ( avgLandIceFreshwaterFlux(iCell) * nAccumulatedCoupled &
                                              + landIceFreshwaterFlux(iCell) ) / ( nAccumulatedCoupled + 1)
              avgLandIceHeatFlux(iCell) = ( avgLandIceHeatFlux(iCell) * nAccumulatedCoupled &
                                        + landIceHeatFlux(iCell) ) / ( nAccumulatedCoupled + 1)
           end do
           !$omp end do
           !$omp end parallel

        end if

        if (config_remove_ais_river_runoff) then
           call mpas_pool_get_array(forcingPool, 'avgRemovedRiverRunoffFlux', avgRemovedRiverRunoffFlux)
           call mpas_pool_get_array(forcingPool, 'removedRiverRunoffFlux', removedRiverRunoffFlux)

           !$omp parallel
           !$omp do schedule(runtime)
           do iCell = 1, nCells
              avgRemovedRiverRunoffFlux(iCell) = ( avgRemovedRiverRunoffFlux(iCell) * nAccumulatedCoupled &
                                               + removedRiverRunoffFlux(iCell) ) / ( nAccumulatedCoupled + 1)
           end do
           !$omp end do
           !$omp end parallel

        end if

        if (config_remove_ais_ice_runoff) then
           call mpas_pool_get_array(forcingPool, 'avgRemovedIceRunoffFlux', avgRemovedIceRunoffFlux)
           call mpas_pool_get_array(forcingPool, 'avgRemovedIceRunoffHeatFlux', avgRemovedIceRunoffHeatFlux)
           call mpas_pool_get_array(forcingPool, 'removedIceRunoffFlux', removedIceRunoffFlux)

           !$omp parallel
           !$omp do schedule(runtime)
           do iCell = 1, nCells
              avgRemovedIceRunoffFlux(iCell) = ( avgRemovedIceRunoffFlux(iCell) * nAccumulatedCoupled &
                                             + removedIceRunoffFlux(iCell) ) / ( nAccumulatedCoupled + 1)
              avgRemovedIceRunoffHeatFlux(iCell) = ( avgRemovedIceRunoffHeatFlux(iCell) * nAccumulatedCoupled &
                                             + removedIceRunoffFlux(iCell)*latent_heat_fusion_mks ) / ( nAccumulatedCoupled + 1)
           end do
           !$omp end do
           !$omp end parallel

        end if

        if (trim(config_glc_thermal_forcing_coupling_mode) == '3d') then
           call accumulate_glc_3d_thermal_forcing(statePool, forcingPool, timeLevel)
        endif

        !  accumulate BGC coupling fields if necessary
        if (config_use_ecosysTracers) then

         call mpas_pool_get_subpool(forcingPool, 'ecosysAuxiliary', ecosysAuxiliary)
         call mpas_pool_get_array(ecosysAuxiliary, 'CO2_gas_flux', CO2_gas_flux)
         call mpas_pool_get_array(ecosysAuxiliary, 'avgCO2_gas_flux', avgCO2_gas_flux)

         call mpas_pool_get_subpool(forcingPool, 'ecosysSeaIceCoupling', ecosysSeaIceCoupling)

         call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfacePhytoC', avgOceanSurfacePhytoC)
         call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceDIC', avgOceanSurfaceDIC)
         call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceDOCSemiLabile', avgOceanSurfaceDOCSemiLabile)
         call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceNO3', avgOceanSurfaceNO3)
         call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceSiO3', avgOceanSurfaceSiO3)
         call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceNH4', avgOceanSurfaceNH4)
         call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceDOCr', avgOceanSurfaceDOCr)
         call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceFeParticulate', avgOceanSurfaceFeParticulate)
         call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceFeDissolved', avgOceanSurfaceFeDissolved)

         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_array(tracersPool, 'ecosysTracers', ecosysTracers, 1)

           !$omp parallel
           !$omp do schedule(runtime)
           do iCell = 1, nCells

           avgCO2_gas_flux(iCell) = ( avgCO2_gas_flux(iCell) * nAccumulatedCoupled &
              + CO2_gas_flux(iCell) ) / ( nAccumulatedCoupled + 1)

           avgOceanSurfacePhytoC(1,iCell) = ( avgOceanSurfacePhytoC(1,iCell) * nAccumulatedCoupled &
              + ecosysTracers(diatC_ind_MPAS,1,iCell) ) / ( nAccumulatedCoupled + 1)
           avgOceanSurfacePhytoC(2,iCell) = ( avgOceanSurfacePhytoC(2,iCell) * nAccumulatedCoupled &
              + ecosysTracers(spC_ind_MPAS,1,iCell) ) / ( nAccumulatedCoupled + 1)
! Currently no Phaeo in MARBL
!          avgOceanSurfacePhytoC(3,iCell) = ( avgOceanSurfacePhytoC(3,iCell) * nAccumulatedCoupled &
!             + ecosysTracers(phaeoC_ind_MPAS,1,iCell) ) / ( nAccumulatedCoupled + 1)
           avgOceanSurfacePhytoC(3,iCell) = 0.0_RKIND

           avgOceanSurfaceDIC(iCell) = ( avgOceanSurfaceDIC(iCell) * nAccumulatedCoupled &
              + ecosysTracers(dic_ind_MPAS,1,iCell) ) / ( nAccumulatedCoupled + 1)
           avgOceanSurfaceDOCSemiLabile(iCell) = ( avgOceanSurfaceDOCSemiLabile(iCell) * nAccumulatedCoupled &
              + ecosysTracers(doc_ind_MPAS,1,iCell) ) / ( nAccumulatedCoupled + 1)
           avgOceanSurfaceSiO3(iCell) = ( avgOceanSurfaceSiO3(iCell) * nAccumulatedCoupled &
              + ecosysTracers(sio3_ind_MPAS,1,iCell) ) / ( nAccumulatedCoupled + 1)
           avgOceanSurfaceNO3(iCell) = ( avgOceanSurfaceNO3(iCell) * nAccumulatedCoupled &
              + ecosysTracers(no3_ind_MPAS,1,iCell) ) / ( nAccumulatedCoupled + 1)
           avgOceanSurfaceNH4(iCell) = ( avgOceanSurfaceNH4(iCell) * nAccumulatedCoupled &
              + ecosysTracers(nh4_ind_MPAS,1,iCell) ) / ( nAccumulatedCoupled + 1)
           avgOceanSurfaceFeDissolved(iCell) = ( avgOceanSurfaceFeDissolved(iCell) * nAccumulatedCoupled &
              + ecosysTracers(fe_ind_MPAS,1,iCell) ) / ( nAccumulatedCoupled + 1)
           avgOceanSurfaceDOCr(iCell) = ( avgOceanSurfaceDOCr(iCell) * nAccumulatedCoupled &
              + ecosysTracers(docr_ind_MPAS,1,iCell) ) / ( nAccumulatedCoupled + 1)

           avgOceanSurfaceFeParticulate(iCell) = 0.0_RKIND

           end do
           !$omp end do
           !$omp end parallel
        end if

        if (config_use_DMSTracers) then
           call mpas_pool_get_subpool(forcingPool, 'DMSSeaIceCoupling', DMSSeaIceCoupling)

           call mpas_pool_get_array(DMSSeaIceCoupling, 'avgOceanSurfaceDMS', avgOceanSurfaceDMS)
           call mpas_pool_get_array(DMSSeaIceCoupling, 'avgOceanSurfaceDMSP', avgOceanSurfaceDMSP)

           call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
           call mpas_pool_get_array(tracersPool, 'DMSTracers', DMSTracers, 1)

           !$omp parallel
           !$omp do schedule(runtime)
           do iCell = 1, nCells
           avgOceanSurfaceDMS(iCell) = ( avgOceanSurfaceDMS(iCell) * nAccumulatedCoupled &
              + DMSTracers(dmsIndices%dms_ind,1,iCell) ) / ( nAccumulatedCoupled + 1)
           avgOceanSurfaceDMSP(iCell) = ( avgOceanSurfaceDMSP(iCell) * nAccumulatedCoupled &
              + DMSTracers(dmsIndices%dmsp_ind,1,iCell) ) / ( nAccumulatedCoupled + 1)
           end do
           !$omp end do
           !$omp end parallel
        endif

        if (config_use_MacroMoleculesTracers) then
           call mpas_pool_get_subpool(forcingPool, 'MacroMoleculesSeaIceCoupling', MacroMoleculesSeaIceCoupling)

           call mpas_pool_get_array(MacroMoleculesSeaIceCoupling, 'avgOceanSurfaceDOC', avgOceanSurfaceDOC)
           call mpas_pool_get_array(MacroMoleculesSeaIceCoupling, 'avgOceanSurfaceDON', avgOceanSurfaceDON)

           call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
           call mpas_pool_get_array(tracersPool, 'MacroMoleculesTracers', MacroMoleculesTracers, 1)

           !$omp parallel
           !$omp do schedule(runtime)
           do iCell = 1, nCells
           avgOceanSurfaceDOC(1,iCell) = ( avgOceanSurfaceDOC(1,iCell) * nAccumulatedCoupled &
              + MacroMoleculesTracers(macrosIndices%poly_ind,1,iCell) ) / ( nAccumulatedCoupled + 1)
           avgOceanSurfaceDOC(2,iCell) = ( avgOceanSurfaceDOC(2,iCell) * nAccumulatedCoupled &
              + MacroMoleculesTracers(macrosIndices%lip_ind,1,iCell) ) / ( nAccumulatedCoupled + 1)
!maltrud need to renormalize
           avgOceanSurfaceDON(iCell) = ( avgOceanSurfaceDON(iCell) * nAccumulatedCoupled &
              + MacroMoleculesTracers(macrosIndices%prot_ind,1,iCell) ) / ( nAccumulatedCoupled + 1)
           end do
           !$omp end do
           !$omp end parallel
        endif

        nAccumulatedCoupled = nAccumulatedCoupled + 1

    end subroutine ocn_time_average_coupled_accumulate!}}}


!***********************************************************************
!
!  routine accumulate_glc_3d_thermal_forcing
!
!> \brief   Coupled time averager accumulation for 3D thermal forcing
!> \author  Xylar Asay-Davis
!> \date    03/05/2025
!> \details
!>  This routine accumulated the coupled time average of the 3D thermal
!>  forcing
!
!-----------------------------------------------------------------------
    subroutine accumulate_glc_3d_thermal_forcing(statePool, forcingPool, timeLevel)!{{{
        use ocn_constants, only: &
             latent_heat_fusion_mks

        type (mpas_pool_type), intent(in) :: statePool
        type (mpas_pool_type), intent(inout) :: forcingPool
        integer, intent(in) :: timeLevel

        integer :: iCell
        real (kind=RKIND), dimension(:,:), pointer :: avgThermalForcingAtZLevels
        real (kind=RKIND), dimension(:,:), pointer :: avgThermalForcingAtZLevelsMask
        real (kind=RKIND), dimension(:), pointer :: glcZLevels
        real (kind=RKIND), dimension(:), pointer :: ssh

        type (mpas_pool_type), pointer :: tracersPool

        real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers
        integer, pointer :: indexTemperature, indexSalinity, nGlcZLevels, &
                            nCells, nAccumulatedCoupled
        integer :: iLevel, iLevelGlc, iLevelAbove, iLevelBelow

        real (kind=RKIND) :: freezingTemp, glcz, zMin, zMax, zAbove, zBelow, &
                             alpha, salin, temp, press

        call mpas_pool_get_dimension(forcingPool, 'nCells', nCells)
        call mpas_pool_get_dimension(forcingPool, 'nGlcZLevels', nGlcZLevels)

        call mpas_pool_get_array(forcingPool, 'nAccumulatedCoupled', nAccumulatedCoupled)

        call mpas_pool_get_array(forcingPool, 'glcZLevels', glcZLevels)
        call mpas_pool_get_array(statePool, 'ssh', ssh, timeLevel)

        call mpas_pool_get_array(forcingPool, 'avgThermalForcingAtZLevels', avgThermalForcingAtZLevels)
        call mpas_pool_get_array(forcingPool, 'avgThermalForcingAtZLevelsMask', avgThermalForcingAtZLevelsMask)
        call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
        call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, timeLevel)
        call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperature)
        call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)

        !$omp parallel
        !$omp do collapse(2) schedule(runtime) &
        !$omp    private(iLevelGlc, iCell, glcz, zMin, zMax, iLevelAbove, &
        !$omp            iLevelBelow, zAbove, zBelow, alpha, salin, &
        !$omp            temp, press, freezingTemp)
        do iLevelGlc = 1, nGlcZLevels
            do iCell = 1, nCells
                ! calculate thermal forcing at identified level for each cell
                glcz = glcZLevels(iLevelGlc)

                if (glcz > ssh(iCell) .or. glcz < -bottomDepth(iCell) .or. &
                    minLevelCell(iCell) > maxLevelCell(iCell)) then
                    ! mask out glc levels that are too shallow or too deep
                    avgThermalForcingAtZLevels(iLevelGlc, iCell) = &
                        ( avgThermalForcingAtZLevels(iLevelGlc, iCell) * nAccumulatedCoupled &
                          + 0.0_RKIND) / (nAccumulatedCoupled + 1)
                    avgThermalForcingAtZLevelsMask(iLevelGlc, iCell) = &
                        ( avgThermalForcingAtZLevelsMask(iLevelGlc, iCell) * nAccumulatedCoupled &
                          + 0.0_RKIND) / ( nAccumulatedCoupled + 1)
                    cycle
                end if

                ! find which zMid values glcz lies between and the linear
                ! interpolation factor alpha between them

                zMin = zMid(minLevelCell(iCell), iCell)
                zMax = zMid(maxLevelCell(iCell), iCell)
                if (glcz > zMin .or. minLevelCell(iCell) == maxLevelCell(iCell)) then
                    ! above the middle of the top level so no interpolation
                    iLevelAbove = minLevelCell(iCell)
                    iLevelBelow = minLevelCell(iCell)
                    alpha = 0.0_RKIND
                else if (glcz < zMax) then
                    ! below the middle of the bottom level so no interpolation
                    iLevelAbove = maxLevelCell(iCell)
                    iLevelBelow = maxLevelCell(iCell)
                    alpha = 0.0_RKIND
                else
                    do iLevel = minLevelCell(iCell), maxLevelCell(iCell) - 1
                        iLevelAbove = iLevel
                        iLevelBelow = iLevel + 1
                        zAbove = zMid(iLevelAbove, iCell)
                        zBelow = zMid(iLevelBelow, iCell)
                        if (glcz >= zBelow) then
                            exit
                        end if
                    end do
                    alpha = (glcz - zAbove) / (zBelow - zAbove)
                end if

                ! linearly interpolate salinity, temperature, and pressure
                salin = activeTracers(indexSalinity, iLevelAbove, iCell) * (1.0_RKIND - alpha) &
                      + activeTracers(indexSalinity, iLevelBelow, iCell) * alpha

                temp = activeTracers(indexTemperature, iLevelAbove, iCell) * (1.0_RKIND - alpha) &
                     + activeTracers(indexTemperature, iLevelBelow, iCell) * alpha

                press = pressure(iLevelAbove, iCell) * (1.0_RKIND - alpha) &
                      + pressure(iLevelBelow, iCell) * alpha

                freezingTemp = ocn_freezing_temperature(salinity=salin, pressure=press, &
                                                        inLandIceCavity=.true.)

                avgThermalForcingAtZLevels(iLevelGlc, iCell) = &
                    ( avgThermalForcingAtZLevels(iLevelGlc, iCell) * nAccumulatedCoupled &
                      + temp - freezingTemp ) / ( nAccumulatedCoupled + 1)
                avgThermalForcingAtZLevelsMask(iLevelGlc, iCell) = &
                    ( avgThermalForcingAtZLevelsMask(iLevelGlc, iCell) * nAccumulatedCoupled &
                     + 1.0_RKIND) / ( nAccumulatedCoupled + 1)
            end do
        end do
        !$omp end do
        !$omp end parallel

    end subroutine accumulate_glc_3d_thermal_forcing!}}}

end module ocn_time_average_coupled
