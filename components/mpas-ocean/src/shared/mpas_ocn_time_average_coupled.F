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

        real (kind=RKIND), dimension(:), pointer :: avgEffectiveDensityInLandIce, avgTotalFreshWaterTemperatureFlux

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
                                                    totalFreshWaterTemperatureFlux, avgTotalFreshWaterTemperatureFlux

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

end module ocn_time_average_coupled
