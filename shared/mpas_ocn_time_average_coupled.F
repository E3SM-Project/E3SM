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

        nAccumulatedCoupled = nAccumulatedCoupled + 1

    end subroutine ocn_time_average_coupled_accumulate!}}}

end module ocn_time_average_coupled
