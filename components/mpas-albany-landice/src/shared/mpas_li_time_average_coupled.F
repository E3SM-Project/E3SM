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
    subroutine li_time_average_coupled_init(meshPool,timeAveragingPool)

        type (mpas_pool_type), intent(in) :: meshPool          
        type (mpas_pool_type), intent(inout) :: timeAveragingPool

        real (kind=RKIND), dimension(:), pointer :: avgBareIceAblationApplied, avgCalvingFlux, &
                avgFaceMeltFlux, avgFloatingBMBFlux
        ! fields for ISMIP reporting not included in the above
        real (kind=RKIND), dimension(:), pointer :: avgSMBFlux, avgGroundedBMBFlux, avgGroundingLineFlux, avgDhdt

        integer, pointer :: nCells

        real (kind=RKIND), pointer :: timeAccumulatedCoupled 

        integer :: iCell

        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

        call mpas_pool_get_array(timeAveragingPool, 'timeAccumulatedCoupled', timeAccumulatedCoupled)
        ! fields for E3SM coupling
        call mpas_pool_get_array(timeAveragingPool, 'avgBareIceAblationApplied', avgBareIceAblationApplied)
        call mpas_pool_get_array(timeAveragingPool, 'avgCalvingFlux', avgCalvingFlux)
        call mpas_pool_get_array(timeAveragingPool, 'avgFaceMeltFlux', avgFaceMeltFlux)
        call mpas_pool_get_array(timeAveragingPool, 'avgFloatingBMBFlux', avgFloatingBMBFlux)
        ! fields for ISMIP reporting not included in the above
        call mpas_pool_get_array(timeAveragingPool, 'avgSMBFlux', avgSMBFlux)
        call mpas_pool_get_array(timeAveragingPool, 'avgGroundedBMBFlux', avgGroundedBMBFlux)
        call mpas_pool_get_array(timeAveragingPool, 'avgGroundingLineFlux', avgGroundingLineFlux)
        call mpas_pool_get_array(timeAveragingPool, 'avgDhdt', avgDhdt)

        do iCell = 1, nCells
           avgBareIceAblationApplied(iCell) = 0.0_RKIND
           avgCalvingFlux(iCell) = 0.0_RKIND
           avgFaceMeltFlux(iCell) = 0.0_RKIND
           avgFloatingBMBFlux(iCell) = 0.0_RKIND
           ! fields for ISMIP reporting not included in the above
           avgSMBFlux(iCell) = 0.0_RKIND
           avgGroundedBMBFlux(iCell) = 0.0_RKIND
           avgGroundingLineFlux(iCell) = 0.0_RKIND
           avgDhdt(iCell) = 0.0_RKIND
        end do

        timeAccumulatedCoupled = 0.0_RKIND

    end subroutine li_time_average_coupled_init

!***********************************************************************
!
!  routine li_time_average_coupled_accumulate
!
!> \brief   Coupled time averager accumulation
!> \author  Stephen Price
!> \date    04 March 2025
!> \details
!>  This routine accumulates the coupled time averaging fields
!
!-----------------------------------------------------------------------
    subroutine li_time_average_coupled_accumulate(timeAveragingPool, geometryPool, meshPool, velocityPool)
            
        use li_setup 

        type (mpas_pool_type), intent(inout) :: timeAveragingPool
        type (mpas_pool_type), intent(in) :: geometryPool
        type (mpas_pool_type), intent(in) :: meshPool
        type (mpas_pool_type), intent(in) :: velocityPool

        real (kind=RKIND), dimension(:), pointer :: bareIceAblationApplied, avgBareIceAblationApplied, &
                                                    calvingThickness, avgCalvingFlux, &
                                                    faceMeltingThickness, avgFaceMeltFlux, avgFloatingBMBFlux, &
                                                    floatingBasalMassBalApplied
        ! fields for ISMIP reporting not included in the above
        real (kind=RKIND), dimension(:), pointer :: avgSMBFlux, avgGroundedBMBFlux, avgGroundingLineFlux, avgDhdt
        real (kind=RKIND), dimension(:), pointer :: sfcMassBalApplied, groundedBasalMassBalApplied, &
                                                    fluxAcrossGroundingLineOnCells, dHdt
        integer, pointer :: nCells

        real (kind=RKIND), pointer :: timeAccumulatedCoupled, rhoi, deltat

        integer :: iCell

        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
        call mpas_pool_get_config(liConfigs, 'config_ice_density', rhoi)
        call mpas_pool_get_array(meshPool, 'deltat', deltat)
        call mpas_pool_get_array(timeAveragingPool, 'timeAccumulatedCoupled', timeAccumulatedCoupled)
        call mpas_pool_get_array(geometryPool, 'bareIceAblationApplied', bareIceAblationApplied)
        call mpas_pool_get_array(geometryPool, 'calvingThickness', calvingThickness)
        call mpas_pool_get_array(geometryPool, 'faceMeltingThickness', faceMeltingThickness)
        call mpas_pool_get_array(geometryPool, 'floatingBasalMassBalApplied', floatingBasalMassBalApplied)
        call mpas_pool_get_array(timeAveragingPool, 'avgBareIceAblationApplied', avgBareIceAblationApplied)
        call mpas_pool_get_array(timeAveragingPool, 'avgCalvingFlux', avgCalvingFlux)
        call mpas_pool_get_array(timeAveragingPool, 'avgFaceMeltFlux', avgFaceMeltFlux)
        call mpas_pool_get_array(timeAveragingPool, 'avgFloatingBMBFlux', avgFloatingBMBFlux)
        ! fields for ISMIP reporting not included in the above
        call mpas_pool_get_array(timeAveragingPool, 'avgSMBFlux', avgSMBFlux)
        call mpas_pool_get_array(timeAveragingPool, 'avgGroundedBMBFlux', avgGroundedBMBFlux)
        call mpas_pool_get_array(timeAveragingPool, 'avgGroundingLineFlux', avgGroundingLineFlux)
        call mpas_pool_get_array(timeAveragingPool, 'avgDhdt', avgDhdt)
        call mpas_pool_get_array(geometryPool, 'sfcMassBalApplied', sfcMassBalApplied)
        call mpas_pool_get_array(geometryPool, 'groundedBasalMassBalApplied', groundedBasalMassBalApplied)
        call mpas_pool_get_array(velocityPool, 'fluxAcrossGroundingLineOnCells', fluxAcrossGroundingLineOnCells)
        call mpas_pool_get_array(geometryPool, 'dHdt', dHdt)

        do iCell = 1, nCells

           avgBareIceAblationApplied(iCell) = ( avgBareIceAblationApplied(iCell) * timeAccumulatedCoupled &
                                               + bareIceAblationApplied(iCell) * deltat ) / ( timeAccumulatedCoupled + deltat )

           avgCalvingFlux(iCell) = ( avgCalvingFlux(iCell) * timeAccumulatedCoupled &
                                                + calvingThickness(iCell) * deltat * (rhoi / deltat) ) / ( timeAccumulatedCoupled + deltat )

           avgFaceMeltFlux(iCell) = ( avgFaceMeltFlux(iCell) * timeAccumulatedCoupled &
                                                + faceMeltingThickness(iCell) * deltat * (rhoi / deltat) ) / ( timeAccumulatedCoupled + deltat )
                                        
           avgFloatingBMBFlux(iCell) = ( avgFloatingBMBFlux(iCell) * timeAccumulatedCoupled &
                                                + floatingBasalMassBalApplied(iCell) * deltat) / ( timeAccumulatedCoupled + deltat )

           ! fields for ISMIP reporting not included in the above
           avgSMBFlux(iCell) = ( avgSMBFlux(iCell) * timeAccumulatedCoupled &
                                                + sfcMassBalApplied(iCell) * deltat) / ( timeAccumulatedCoupled + deltat )

           avgGroundedBMBFlux(iCell) = ( avgGroundedBMBFlux(iCell) * timeAccumulatedCoupled &
                                                + groundedBasalMassBalApplied(iCell) * deltat) / ( timeAccumulatedCoupled + deltat )

           avgGroundingLineFlux(iCell) = ( avgGroundingLineFlux(iCell) * timeAccumulatedCoupled &
                                                + fluxAcrossGroundingLineOnCells(iCell) * deltat) / ( timeAccumulatedCoupled + deltat )

           avgDhdt(iCell) = ( avgDhdt(iCell) * timeAccumulatedCoupled &
                                                + dHdt(iCell) * deltat) / ( timeAccumulatedCoupled + deltat )

        end do

        timeAccumulatedCoupled = timeAccumulatedCoupled + deltat

    end subroutine li_time_average_coupled_accumulate

end module li_time_average_coupled
