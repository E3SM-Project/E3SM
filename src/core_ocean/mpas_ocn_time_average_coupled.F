! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
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

    use mpas_grid_types
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
    subroutine ocn_time_average_coupled_init(forcing)!{{{
        type (forcing_type), intent(inout) :: forcing

        real (kind=RKIND), dimension(:,:), pointer :: avgTracersSurfaceValue, avgSurfaceVelocity

        avgTracersSurfaceValue => forcing % avgTracersSurfaceValue % array
        avgSurfaceVelocity => forcing % avgSurfaceVelocity % array

        avgTracersSurfaceValue(:,:) = 0.0_RKIND
        avgSurfaceVelocity(:,:) = 0.0_RKIND
        forcing % nAccumulatedCoupled % scalar = 0

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
    subroutine ocn_time_average_coupled_accumulate(diagnostics, forcing)!{{{
        type (diagnostics_type), intent(in) :: diagnostics
        type (forcing_type), intent(inout) :: forcing

        real (kind=RKIND), dimension(:,:), pointer :: surfaceVelocity, avgSurfaceVelocity
        real (kind=RKIND), dimension(:,:), pointer :: tracersSurfaceValue
        real (kind=RKIND), dimension(:,:), pointer :: avgTracersSurfaceValue
        integer :: index_temperature, index_salinity

        tracersSurfaceValue => diagnostics % tracersSurfaceValue % array
        surfaceVelocity => diagnostics % surfaceVelocity % array

        avgTracersSurfaceValue => forcing % avgTracersSurfaceValue % array
        avgSurfaceVelocity => forcing % avgSurfaceVelocity % array

        index_temperature = forcing % index_temperatureSurfaceValue

        avgTracersSurfaceValue(:,:) = avgTracersSurfaceValue(:,:) + tracersSurfaceValue(:,:)
        avgSurfaceVelocity(:,:) = avgSurfaceVelocity(:,:) + surfaceVelocity(:,:)

        avgTracersSurfaceValue(index_temperature,:) = avgTracersSurfaceValue(index_temperature,:) + tracersSurfaceValue(index_temperature,:) + T0_Kelvin

        forcing % nAccumulatedCoupled % scalar = forcing % nAccumulatedCoupled % scalar + 1

    end subroutine ocn_time_average_coupled_accumulate!}}}

!***********************************************************************
!
!  routine ocn_time_average_coupled_normalize
!
!> \brief   Coupled time averager normalization
!> \author  Doug Jacobsen
!> \date    06/08/2013
!> \details 
!>  This routine normalizes the coupled time averaging fields
!
!-----------------------------------------------------------------------
    subroutine ocn_time_average_coupled_normalize(forcing)!{{{

        type (forcing_type), intent(inout) :: forcing

        real (kind=RKIND), dimension(:,:), pointer :: avgTracersSurfaceValue, avgSurfaceVelocity

        avgTracersSurfaceValue => forcing % avgTracersSurfaceValue % array
        avgSurfaceVelocity => forcing % avgSurfaceVelocity % array

        if(forcing % nAccumulatedCoupled % scalar > 0) then
           avgTracersSurfaceValue = avgTracersSurfaceValue / forcing % nAccumulatedCoupled % scalar
           avgSurfaceVelocity = avgSurfaceVelocity / forcing % nAccumulatedCoupled % scalar
           forcing % nAccumulatedCoupled % scalar = 0
        end if

    end subroutine ocn_time_average_coupled_normalize!}}}

end module ocn_time_average_coupled
