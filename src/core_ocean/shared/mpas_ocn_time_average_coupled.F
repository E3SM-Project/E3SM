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
    subroutine ocn_time_average_coupled_init(forcingPool)!{{{
        type (mpas_pool_type), intent(inout) :: forcingPool

        real (kind=RKIND), dimension(:,:), pointer :: avgTracersSurfaceValue, avgSurfaceVelocity, avgSSHGradient

        integer, pointer :: nAccumulatedCoupled

        call mpas_pool_get_array(forcingPool, 'avgTracersSurfaceValue', avgTracersSurfaceValue)
        call mpas_pool_get_array(forcingPool, 'avgSurfaceVelocity', avgSurfaceVelocity)
        call mpas_pool_get_array(forcingPool, 'avgSSHGradient', avgSSHGradient)
        call mpas_pool_get_array(forcingPool, 'nAccumulatedCoupled', nAccumulatedCoupled)

        avgTracersSurfaceValue(:,:) = 0.0_RKIND
        avgSurfaceVelocity(:,:) = 0.0_RKIND
        avgSSHGradient(:,:) = 0.0_RKIND

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
    subroutine ocn_time_average_coupled_accumulate(diagnosticsPool, forcingPool)!{{{
        type (mpas_pool_type), intent(in) :: diagnosticsPool
        type (mpas_pool_type), intent(inout) :: forcingPool

        real (kind=RKIND), dimension(:,:), pointer :: surfaceVelocity, avgSurfaceVelocity
        real (kind=RKIND), dimension(:,:), pointer :: tracersSurfaceValue, avgTracersSurfaceValue
        real (kind=RKIND), dimension(:,:), pointer :: avgSSHGradient
        real (kind=RKIND), dimension(:,:), pointer :: gradSSHZonal, gradSSHMeridional
        integer, pointer :: index_temperature, index_SSHzonal, index_SSHmeridional, nAccumulatedCoupled

        call mpas_pool_get_array(diagnosticsPool, 'tracersSurfaceValue', tracersSurfaceValue)
        call mpas_pool_get_array(diagnosticsPool, 'surfaceVelocity', surfaceVelocity)
        call mpas_pool_get_array(diagnosticsPool, 'gradSSHZonal', gradSSHZonal)
        call mpas_pool_get_array(diagnosticsPool, 'gradSSHMeridional', gradSSHMeridional)

        call mpas_pool_get_array(forcingPool, 'avgTracersSurfaceValue', avgTracersSurfaceValue)
        call mpas_pool_get_array(forcingPool, 'avgSurfaceVelocity', avgSurfaceVelocity)
        call mpas_pool_get_array(forcingPool, 'avgSSHGradient', avgSSHGradient)

        call mpas_pool_get_dimension(forcingPool, 'index_temperatureSurfaceValue', index_temperature)
        call mpas_pool_get_dimension(forcingPool, 'index_avgSSHGradientZonal', index_SSHzonal)
        call mpas_pool_get_dimension(forcingPool, 'index_avgSSHGradientMeridional', index_SSHmeridional)

        call mpas_pool_get_array(forcingPool, 'nAccumulatedCoupled', nAccumulatedCoupled)

        avgTracersSurfaceValue(:,:) = avgTracersSurfaceValue(:,:) * nAccumulatedCoupled + tracersSurfaceValue(:,:)
        avgTracersSurfaceValue(index_temperature,:) = avgTracersSurfaceValue(index_temperature,:) + T0_Kelvin
        avgTracersSurfaceValue(:,:) = avgTracersSurfaceValue(:,:) / ( nAccumulatedCoupled + 1 )

        avgSurfaceVelocity(:,:)     = ( avgSurfaceVelocity(:,:)     * nAccumulatedCoupled + surfaceVelocity(:,:)     ) / ( nAccumulatedCoupled + 1 )

        avgSSHGradient(index_SSHzonal,:)      = ( avgSSHGradient(index_SSHzonal,:)      * nAccumulatedCoupled + gradSSHZonal(1,:) ) / ( nAccumulatedCoupled + 1 )
        avgSSHGradient(index_SSHmeridional,:) = ( avgSSHGradient(index_SSHmeridional,:) * nAccumulatedCoupled + gradSSHMeridional(1,:) ) / ( nAccumulatedCoupled + 1 )

        nAccumulatedCoupled = nAccumulatedCoupled + 1

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
    subroutine ocn_time_average_coupled_normalize(forcingPool)!{{{

        type (mpas_pool_type), intent(inout) :: forcingPool

!       real (kind=RKIND), dimension(:,:), pointer :: avgTracersSurfaceValue, avgSurfaceVelocity, avgSSHGradient

!       avgTracersSurfaceValue => forcing % avgTracersSurfaceValue % array
!       avgSurfaceVelocity => forcing % avgSurfaceVelocity % array
!       avgSSHGradient => forcing % avgSSHGradient % array

!       if(forcing % nAccumulatedCoupled % scalar > 0) then
!          avgTracersSurfaceValue = avgTracersSurfaceValue / forcing % nAccumulatedCoupled % scalar
!          avgSurfaceVelocity = avgSurfaceVelocity / forcing % nAccumulatedCoupled % scalar
!          avgSSHGradient = avgSSHGradient / forcing % nAccumulatedCoupled % scalar
!          forcing % nAccumulatedCoupled % scalar = 0
!       end if

    end subroutine ocn_time_average_coupled_normalize!}}}

end module ocn_time_average_coupled
