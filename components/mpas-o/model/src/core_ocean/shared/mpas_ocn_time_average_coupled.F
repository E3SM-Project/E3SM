! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
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
    use ocn_tracer_ecosys
    use ocn_tracer_DMS
    use ocn_tracer_MacroMolecules

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

        real (kind=RKIND), dimension(:), pointer :: avgEffectiveDensityInLandIce
        character (len=StrKIND), pointer :: config_land_ice_flux_mode

        integer :: iCell
        integer, pointer :: nAccumulatedCoupled, nCells

        logical, pointer :: config_use_ecosysTracers, &
                            config_use_DMSTracers,    &
                            config_use_MacroMoleculesTracers

        type (mpas_pool_type), pointer :: ecosysSeaIceCoupling, &
                                          DMSSeaIceCoupling,    &
                                          MacroMoleculesSeaIceCoupling

        real (kind=RKIND), dimension(:), pointer :: avgOceanSurfaceDIC, &
                                                    avgOceanSurfaceDON, &
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

        call mpas_pool_get_dimension(forcingPool, 'nCells', nCells)

        call mpas_pool_get_array(forcingPool, 'avgTracersSurfaceValue', avgTracersSurfaceValue)
        call mpas_pool_get_array(forcingPool, 'avgSurfaceVelocity', avgSurfaceVelocity)
        call mpas_pool_get_array(forcingPool, 'avgSSHGradient', avgSSHGradient)
        call mpas_pool_get_array(forcingPool, 'nAccumulatedCoupled', nAccumulatedCoupled)

        !$omp do schedule(runtime)
        do iCell = 1, nCells
           avgSurfaceVelocity(:, iCell) = 0.0_RKIND
           avgTracersSurfaceValue(:, iCell) = 0.0_RKIND
           avgSSHGradient(:, iCell) = 0.0_RKIND
        end do
        !$omp end do

        call mpas_pool_get_config(ocnConfigs, 'config_land_ice_flux_mode', config_land_ice_flux_mode)
        if(trim(config_land_ice_flux_mode) == 'coupled') then
           call mpas_pool_get_array(forcingPool, 'avgLandIceBoundaryLayerTracers', avgLandIceBoundaryLayerTracers)
           call mpas_pool_get_array(forcingPool, 'avgLandIceTracerTransferVelocities', avgLandIceTracerTransferVelocities)
           call mpas_pool_get_array(forcingPool, 'avgEffectiveDensityInLandIce', avgEffectiveDensityInLandIce)

           !$omp do schedule(runtime)
           do iCell = 1, nCells
              avgLandIceBoundaryLayerTracers(:, iCell) = 0.0_RKIND
              avgLandIceTracerTransferVelocities(:, iCell) = 0.0_RKIND
              avgEffectiveDensityInLandIce(iCell) = 0.0_RKIND
           end do
           !$omp end do
        end if

        !  set up BGC coupling fields if necessary
        call mpas_pool_get_config(ocnConfigs, 'config_use_ecosysTracers', config_use_ecosysTracers)
        call mpas_pool_get_config(ocnConfigs, 'config_use_DMSTracers', config_use_DMSTracers)
        call mpas_pool_get_config(ocnConfigs, 'config_use_MacroMoleculesTracers', config_use_MacroMoleculesTracers)

        if (config_use_ecosysTracers) then

           call mpas_pool_get_subpool(forcingPool, 'ecosysSeaIceCoupling', ecosysSeaIceCoupling)

           call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfacePhytoC', avgOceanSurfacePhytoC)
           call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceDIC', avgOceanSurfaceDIC)
           call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceNO3', avgOceanSurfaceNO3)
           call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceSiO3', avgOceanSurfaceSiO3)
           call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceNH4', avgOceanSurfaceNH4)
           call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceDOCr', avgOceanSurfaceDOCr)
           call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceFeParticulate', avgOceanSurfaceFeParticulate)
           call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceFeDissolved', avgOceanSurfaceFeDissolved)

           !$omp do schedule(runtime)
           do iCell = 1, nCells
              avgOceanSurfacePhytoC(:,iCell) = 0.0_RKIND

              avgOceanSurfaceDIC(iCell) = 0.0_RKIND
              avgOceanSurfaceNO3(iCell) = 0.0_RKIND
              avgOceanSurfaceSiO3(iCell) = 0.0_RKIND
              avgOceanSurfaceNH4(iCell) = 0.0_RKIND
              avgOceanSurfaceDOCr(iCell) = 0.0_RKIND
              avgOceanSurfaceFeParticulate(iCell) = 0.0_RKIND
              avgOceanSurfaceFeDissolved(iCell) = 0.0_RKIND
           end do
           !$omp end do
        end if

        if (config_use_DMSTracers) then
           call mpas_pool_get_subpool(forcingPool, 'DMSSeaIceCoupling', DMSSeaIceCoupling)

           call mpas_pool_get_array(DMSSeaIceCoupling, 'avgOceanSurfaceDMS', avgOceanSurfaceDMS)
           call mpas_pool_get_array(DMSSeaIceCoupling, 'avgOceanSurfaceDMSP', avgOceanSurfaceDMSP)

           !$omp do schedule(runtime)
           do iCell = 1, nCells
              avgOceanSurfaceDMS(iCell) = 0.0_RKIND
              avgOceanSurfaceDMSP(iCell) = 0.0_RKIND
           end do
           !$omp end do
        endif
        if (config_use_MacroMoleculesTracers) then
           call mpas_pool_get_subpool(forcingPool, 'MacroMoleculesSeaIceCoupling', MacroMoleculesSeaIceCoupling)

           call mpas_pool_get_array(MacroMoleculesSeaIceCoupling, 'avgOceanSurfaceDOC', avgOceanSurfaceDOC)
           call mpas_pool_get_array(MacroMoleculesSeaIceCoupling, 'avgOceanSurfaceDON', avgOceanSurfaceDON)

           !$omp do schedule(runtime)
           do iCell = 1, nCells
              avgOceanSurfaceDOC(:,iCell) = 0.0_RKIND
              avgOceanSurfaceDON(iCell) = 0.0_RKIND
           end do
           !$omp end do
        endif

        call mpas_threading_barrier()
        !$omp master

        nAccumulatedCoupled = 0

        !$omp end master
        call mpas_threading_barrier()

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
    subroutine ocn_time_average_coupled_accumulate(diagnosticsPool, statePool, forcingPool, timeLevel)!{{{
        type (mpas_pool_type), intent(in) :: diagnosticsPool
        type (mpas_pool_type), intent(in) :: statePool
        type (mpas_pool_type), intent(inout) :: forcingPool
        integer, intent(in) :: timeLevel

        real (kind=RKIND), dimension(:,:), pointer :: surfaceVelocity, avgSurfaceVelocity
        real (kind=RKIND), dimension(:,:), pointer :: tracersSurfaceValue, avgTracersSurfaceValue
        real (kind=RKIND), dimension(:,:), pointer :: avgSSHGradient
        real (kind=RKIND), dimension(:), pointer :: gradSSHZonal, gradSSHMeridional
        integer :: iCell
        integer, pointer :: index_temperature, index_SSHzonal, index_SSHmeridional, nAccumulatedCoupled, nCells
        real (kind=RKIND), dimension(:,:), pointer :: landIceBoundaryLayerTracers, landIceTracerTransferVelocities, &
                                                      avgLandIceBoundaryLayerTracers, avgLandIceTracerTransferVelocities
        real (kind=RKIND), dimension(:), pointer :: effectiveDensityInLandIce, avgEffectiveDensityInLandIce
        character (len=StrKIND), pointer :: config_land_ice_flux_mode

        logical, pointer :: config_use_ecosysTracers, &
                            config_use_DMSTracers,    &
                            config_use_MacroMoleculesTracers

        type (mpas_pool_type), pointer :: tracersPool

        real (kind=RKIND), dimension(:,:,:), pointer :: &
           ecosysTracers,  &
           DMSTracers,     &
           MacroMoleculesTracers

        type (mpas_pool_type), pointer :: ecosysSeaIceCoupling, &
                                          DMSSeaIceCoupling,    &
                                          MacroMoleculesSeaIceCoupling

        real (kind=RKIND), dimension(:), pointer :: avgOceanSurfaceDIC, &
                                                    avgOceanSurfaceDON, &
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

        call mpas_pool_get_array(diagnosticsPool, 'tracersSurfaceValue', tracersSurfaceValue)
        call mpas_pool_get_array(diagnosticsPool, 'surfaceVelocity', surfaceVelocity)
        call mpas_pool_get_array(diagnosticsPool, 'gradSSHZonal', gradSSHZonal)
        call mpas_pool_get_array(diagnosticsPool, 'gradSSHMeridional', gradSSHMeridional)

        call mpas_pool_get_array(forcingPool, 'avgTracersSurfaceValue', avgTracersSurfaceValue)
        call mpas_pool_get_array(forcingPool, 'avgSurfaceVelocity', avgSurfaceVelocity)
        call mpas_pool_get_array(forcingPool, 'avgSSHGradient', avgSSHGradient)

        call mpas_pool_get_dimension(forcingPool, 'nCells', nCells)
        call mpas_pool_get_dimension(forcingPool, 'index_avgTemperatureSurfaceValue', index_temperature)
        call mpas_pool_get_dimension(forcingPool, 'index_avgSSHGradientZonal', index_SSHzonal)
        call mpas_pool_get_dimension(forcingPool, 'index_avgSSHGradientMeridional', index_SSHmeridional)

        call mpas_pool_get_array(forcingPool, 'nAccumulatedCoupled', nAccumulatedCoupled)

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
        end do
        !$omp end do

        call mpas_pool_get_config(ocnConfigs, 'config_land_ice_flux_mode', config_land_ice_flux_mode)
        if(trim(config_land_ice_flux_mode) == 'coupled') then
           call mpas_pool_get_array(diagnosticsPool, 'landIceBoundaryLayerTracers', landIceBoundaryLayerTracers)
           call mpas_pool_get_array(diagnosticsPool, 'landIceTracerTransferVelocities', landIceTracerTransferVelocities)
           call mpas_pool_get_array(statePool, 'effectiveDensityInLandIce', effectiveDensityInLandIce, timeLevel)

           call mpas_pool_get_array(forcingPool, 'avgLandIceBoundaryLayerTracers', avgLandIceBoundaryLayerTracers)
           call mpas_pool_get_array(forcingPool, 'avgLandIceTracerTransferVelocities', avgLandIceTracerTransferVelocities)
           call mpas_pool_get_array(forcingPool, 'avgEffectiveDensityInLandIce', avgEffectiveDensityInLandIce)

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
        end if

        !  accumulate BGC coupling fields if necessary
        call mpas_pool_get_config(ocnConfigs, 'config_use_ecosysTracers', config_use_ecosysTracers)
        call mpas_pool_get_config(ocnConfigs, 'config_use_DMSTracers', config_use_DMSTracers)
        call mpas_pool_get_config(ocnConfigs, 'config_use_MacroMoleculesTracers', config_use_MacroMoleculesTracers)

        if (config_use_ecosysTracers) then

         call mpas_pool_get_subpool(forcingPool, 'ecosysSeaIceCoupling', ecosysSeaIceCoupling)

         call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfacePhytoC', avgOceanSurfacePhytoC)
         call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceDIC', avgOceanSurfaceDIC)
         call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceNO3', avgOceanSurfaceNO3)
         call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceSiO3', avgOceanSurfaceSiO3)
         call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceNH4', avgOceanSurfaceNH4)
         call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceDOCr', avgOceanSurfaceDOCr)
         call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceFeParticulate', avgOceanSurfaceFeParticulate)
         call mpas_pool_get_array(ecosysSeaIceCoupling, 'avgOceanSurfaceFeDissolved', avgOceanSurfaceFeDissolved)

         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_array(tracersPool, 'ecosysTracers', ecosysTracers, 1)

           !$omp do schedule(runtime)
           do iCell = 1, nCells

           avgOceanSurfacePhytoC(1,iCell) = ( avgOceanSurfacePhytoC(1,iCell) * nAccumulatedCoupled &
              + ecosysTracers(ecosysIndices%diatC_ind,1,iCell) ) / ( nAccumulatedCoupled + 1)
           avgOceanSurfacePhytoC(2,iCell) = ( avgOceanSurfacePhytoC(2,iCell) * nAccumulatedCoupled &
              + ecosysTracers(ecosysIndices%spC_ind,1,iCell) ) / ( nAccumulatedCoupled + 1)
           avgOceanSurfacePhytoC(3,iCell) = ( avgOceanSurfacePhytoC(3,iCell) * nAccumulatedCoupled &
              + ecosysTracers(ecosysIndices%phaeoC_ind,1,iCell) ) / ( nAccumulatedCoupled + 1)

           avgOceanSurfaceDIC(iCell) = ( avgOceanSurfaceDIC(iCell) * nAccumulatedCoupled &
              + ecosysTracers(ecosysIndices%dic_ind,1,iCell) ) / ( nAccumulatedCoupled + 1)
           avgOceanSurfaceSiO3(iCell) = ( avgOceanSurfaceSiO3(iCell) * nAccumulatedCoupled &
              + ecosysTracers(ecosysIndices%sio3_ind,1,iCell) ) / ( nAccumulatedCoupled + 1)
           avgOceanSurfaceNO3(iCell) = ( avgOceanSurfaceNO3(iCell) * nAccumulatedCoupled &
              + ecosysTracers(ecosysIndices%no3_ind,1,iCell) ) / ( nAccumulatedCoupled + 1)
           avgOceanSurfaceNH4(iCell) = ( avgOceanSurfaceNH4(iCell) * nAccumulatedCoupled &
              + ecosysTracers(ecosysIndices%nh4_ind,1,iCell) ) / ( nAccumulatedCoupled + 1)
           avgOceanSurfaceFeDissolved(iCell) = ( avgOceanSurfaceFeDissolved(iCell) * nAccumulatedCoupled &
              + ecosysTracers(ecosysIndices%fe_ind,1,iCell) ) / ( nAccumulatedCoupled + 1)
!maltrud need to renormalize
           avgOceanSurfaceDOCr(iCell) = ( avgOceanSurfaceDOCr(iCell) * nAccumulatedCoupled &
              + ecosysTracers(ecosysIndices%donr_ind,1,iCell) ) / ( nAccumulatedCoupled + 1)

           avgOceanSurfaceFeParticulate(iCell) = 0.0_RKIND

           end do
           !$omp end do
        end if

        if (config_use_DMSTracers) then
           call mpas_pool_get_subpool(forcingPool, 'DMSSeaIceCoupling', DMSSeaIceCoupling)

           call mpas_pool_get_array(DMSSeaIceCoupling, 'avgOceanSurfaceDMS', avgOceanSurfaceDMS)
           call mpas_pool_get_array(DMSSeaIceCoupling, 'avgOceanSurfaceDMSP', avgOceanSurfaceDMSP)

           call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
           call mpas_pool_get_array(tracersPool, 'DMSTracers', DMSTracers, 1)

           !$omp do schedule(runtime)
           do iCell = 1, nCells
           avgOceanSurfaceDMS(iCell) = ( avgOceanSurfaceDMS(iCell) * nAccumulatedCoupled &
              + DMSTracers(dmsIndices%dms_ind,1,iCell) ) / ( nAccumulatedCoupled + 1)
           avgOceanSurfaceDMSP(iCell) = ( avgOceanSurfaceDMSP(iCell) * nAccumulatedCoupled &
              + DMSTracers(dmsIndices%dmsp_ind,1,iCell) ) / ( nAccumulatedCoupled + 1)
           end do
           !$omp end do
        endif

        if (config_use_MacroMoleculesTracers) then
           call mpas_pool_get_subpool(forcingPool, 'MacroMoleculesSeaIceCoupling', MacroMoleculesSeaIceCoupling)

           call mpas_pool_get_array(MacroMoleculesSeaIceCoupling, 'avgOceanSurfaceDOC', avgOceanSurfaceDOC)
           call mpas_pool_get_array(MacroMoleculesSeaIceCoupling, 'avgOceanSurfaceDON', avgOceanSurfaceDON)

           call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
           call mpas_pool_get_array(tracersPool, 'MacroMoleculesTracers', MacroMoleculesTracers, 1)

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
        endif

        call mpas_threading_barrier()
        !$omp master

        nAccumulatedCoupled = nAccumulatedCoupled + 1

        !$omp end master
        call mpas_threading_barrier()

    end subroutine ocn_time_average_coupled_accumulate!}}}

end module ocn_time_average_coupled
