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

        real (kind=RKIND), dimension(:), pointer :: avgBareIceAblation, avgCalvingThickness !, &

        integer :: iCell
        integer, pointer :: nAccumulatedCoupled, nCells

        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

        call mpas_pool_get_array(timeAveragingPool, 'nAccumulatedCoupled', nAccumulatedCoupled)
        call mpas_pool_get_array(timeAveragingPool, 'avgBareIceAblation', avgBareIceAblation)
        call mpas_pool_get_array(timeAveragingPool, 'avgCalvingThickness', avgCalvingThickness)

        do iCell = 1, nCells
           avgBareIceAblation(iCell) = 0.0_RKIND
           avgCalvingThickness(iCell) = 0.0_RKIND
        end do

        nAccumulatedCoupled = 0

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
    subroutine li_time_average_coupled_accumulate(timeAveragingPool, geometryPool, meshPool)
            
!        use li_constants, only: &
!             constants-list-here   

        type (mpas_pool_type), intent(inout) :: timeAveragingPool
        type (mpas_pool_type), intent(in) :: geometryPool
        type (mpas_pool_type), intent(in) :: meshPool

        real (kind=RKIND), dimension(:), pointer :: avgBareIceAblation, avgCalvingThickness, &
                                                    bareIceAblation, calvingThickness

        integer :: iCell
        integer, pointer :: nAccumulatedCoupled, nCells

        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

        call mpas_pool_get_array(timeAveragingPool, 'bareIceAblation', bareIceAblation)
        call mpas_pool_get_array(timeAveragingPool, 'calvingThickness', calvingThickness)
        call mpas_pool_get_array(timeAveragingPool, 'avgBareIceAblation', avgBareIceAblation)
        call mpas_pool_get_array(timeAveragingPool, 'avgCalvingThickness', avgCalvingThickness)
        call mpas_pool_get_array(timeAveragingPool, 'nAccumulatedCoupled', nAccumulatedCoupled)


        do iCell = 1, nCells

           avgBareIceAblation(iCell) = ( avgBareIceAblation(iCell) * nAccumulatedCoupled &
                                               + bareIceAblation(iCell) ) / ( nAccumulatedCoupled + 1 )

        end do

        nAccumulatedCoupled = nAccumulatedCoupled + 1

    end subroutine li_time_average_coupled_accumulate

end module li_time_average_coupled
