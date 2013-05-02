! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module ocn_time_average

    use mpas_grid_types

    implicit none
    save
    public

    contains 

    subroutine ocn_time_average_init(averagePool)!{{{
        type (mpas_pool_type), intent(inout) :: averagePool

        real (kind=RKIND), pointer :: nAverage

        real (kind=RKIND), dimension(:), pointer :: avgSSH, varSSH
        real (kind=RKIND), dimension(:,:), pointer :: &
           avgNormalVelocity, avgZonalVelocity, avgMeridionalVelocity, avgVertVelocityTop, &
           varNormalVelocity, varZonalVelocity, varMeridionalVelocity, &
           avgNormalTransportVelocity, avgZonalTransportVelocity, avgMeridionalTransportVelocity, avgVertTransportVelocityTop, &
           avgNormalGMBolusVelocity, avgZonalGMBolusVelocity, avgMeridionalGMBolusVelocity, avgVertGMBolusVelocityTop

        call mpas_pool_get_array(averagePool, 'nAverage', nAverage)
        call mpas_pool_get_array(averagePool, 'avgSsh', avgSSH)
        call mpas_pool_get_array(averagePool, 'varSsh', varSSH)
        call mpas_pool_get_array(averagePool, 'avgNormalVelocity', avgNormalVelocity)
        call mpas_pool_get_array(averagePool, 'avgZonalVelocity', avgZonalVelocity)
        call mpas_pool_get_array(averagePool, 'avgMeridionalVelocity', avgMeridionalVelocity)
        call mpas_pool_get_array(averagePool, 'avgVertVelocityTop', avgVertVelocityTop)
        call mpas_pool_get_array(averagePool, 'varNormalVelocity', varNormalVelocity)
        call mpas_pool_get_array(averagePool, 'varZonalVelocity', varZonalVelocity)
        call mpas_pool_get_array(averagePool, 'varMeridionalVelocity', varMeridionalVelocity)
        call mpas_pool_get_array(averagePool, 'avgNormalTransportVelocity', avgNormalTransportVelocity)
        call mpas_pool_get_array(averagePool, 'avgZonalTransportVelocity', avgZonalTransportVelocity)
        call mpas_pool_get_array(averagePool, 'avgMeridionalTransportVelocity', avgMeridionalTransportVelocity)
        call mpas_pool_get_array(averagePool, 'avgVertTransportVelocityTop', avgVertTransportVelocityTop)
        call mpas_pool_get_array(averagePool, 'avgNormalGMBolusVelocity', avgNormalGMBolusVelocity)
        call mpas_pool_get_array(averagePool, 'avgZonalGMBolusVelocity', avgZonalGMBolusVelocity)
        call mpas_pool_get_array(averagePool, 'avgMeridionalGMBolusVelocity', avgMeridionalGMBolusVelocity)
        call mpas_pool_get_array(averagePool, 'avgVertGMBolusVelocityTop', avgVertGMBolusVelocityTop)

        nAverage = 0

        avgSSH = 0.0
        varSSH = 0.0
        avgNormalVelocity = 0.0
        avgZonalVelocity = 0.0
        avgMeridionalVelocity = 0.0
        avgVertVelocityTop = 0.0
        varNormalVelocity = 0.0
        varZonalVelocity = 0.0
        varMeridionalVelocity = 0.0
        avgNormalTransportVelocity = 0.0
        avgZonalTransportVelocity = 0.0
        avgMeridionalTransportVelocity = 0.0
        avgVertTransportVelocityTop = 0.0
        avgNormalGMBolusVelocity = 0.0
        avgZonalGMBolusVelocity = 0.0
        avgMeridionalGMBolusVelocity = 0.0
        avgVertGMBolusVelocityTop = 0.0

    end subroutine ocn_time_average_init!}}}

    subroutine ocn_time_average_accumulate(averagePool, statePool, diagnosticsPool, timeLevelIn)!{{{
        type (mpas_pool_type), intent(inout) :: averagePool
        type (mpas_pool_type), intent(in) :: statePool
        type (mpas_pool_type), intent(in) :: diagnosticsPool
        integer, intent(in), optional :: timeLevelIn

        real (kind=RKIND), pointer :: nAverage, old_nAverage

        real (kind=RKIND), dimension(:), pointer :: ssh
        real (kind=RKIND), dimension(:,:), pointer :: &
           normalVelocityZonal, normalVelocityMeridional, normalVelocity, vertVelocityTop, &
           zonalTransportVelocity, meridionalTransportVelocity, normalTransportVelocity, vertTransportVelocityTop, &
           zonalGMBolusVelocity, meridionalGMBolusVelocity, normalGMBolusVelocity, vertGMBolusVelocityTop

        real (kind=RKIND), dimension(:), pointer :: avgSSH, varSSH
        real (kind=RKIND), dimension(:,:), pointer :: &
           avgNormalVelocity, avgZonalVelocity, avgMeridionalVelocity, avgVertVelocityTop, &
           varNormalVelocity, varZonalVelocity, varMeridionalVelocity, &
           avgNormalTransportVelocity, avgZonalTransportVelocity, avgMeridionalTransportVelocity, avgVertTransportVelocityTop, &
           avgNormalGMBolusVelocity, avgZonalGMBolusVelocity, avgMeridionalGMBolusVelocity, avgVertGMBolusVelocityTop

        real (kind=RKIND), dimension(:), pointer :: old_avgSSH, old_varSSH
        real (kind=RKIND), dimension(:,:), pointer :: &
           old_avgNormalVelocity, old_avgZonalVelocity, old_avgMeridionalVelocity, old_avgVertVelocityTop, &
           old_varNormalVelocity, old_varZonalVelocity, old_varMeridionalVelocity, &
           old_avgNormalTransportVelocity, old_avgZonalTransportVelocity, old_avgMeridionalTransportVelocity, old_avgVertTransportVelocityTop, &
           old_avgNormalGMBolusVelocity, old_avgZonalGMBolusVelocity, old_avgMeridionalGMBolusVelocity, old_avgVertGMBolusVelocityTop

        integer :: timeLevel

        if (present(timeLevelIn)) then
           timeLevel = timeLevelIn
        else
           timeLevel = 1
        end if

        call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)
        call mpas_pool_get_array(statePool, 'ssh', ssh, timeLevel)

        call mpas_pool_get_array(diagnosticsPool, 'normalVelocityZonal', normalVelocityZonal)
        call mpas_pool_get_array(diagnosticsPool, 'normalVelocityMeridional', normalVelocityMeridional)
        call mpas_pool_get_array(diagnosticsPool, 'vertVelocityTop', vertVelocityTop)
        call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity ', normalTransportVelocity)
        call mpas_pool_get_array(diagnosticsPool, 'zonalTransportVelocity', zonalTransportVelocity)
        call mpas_pool_get_array(diagnosticsPool, 'meridionalTransportVelocity', meridionalTransportVelocity)
        call mpas_pool_get_array(diagnosticsPool, 'vertTransportVelocityTop', vertTransportVelocityTop)
        call mpas_pool_get_array(diagnosticsPool, 'normalGMBolusVelocity', normalGMBolusVelocity)
        call mpas_pool_get_array(diagnosticsPool, 'zonalGMBolusVelocity', zonalGMBolusVelocity)
        call mpas_pool_get_array(diagnosticsPool, 'meridionalGMBolusVelocity', meridionalGMBolusVelocity)
        call mpas_pool_get_array(diagnosticsPool, 'vertGMBolusVelocityTop', vertGMBolusVelocityTop)

        call mpas_pool_get_array(averagePool, 'nAverage', nAverage)
        call mpas_pool_get_array(averagePool, 'avgSsh', avgSSH)
        call mpas_pool_get_array(averagePool, 'varSsh', varSSH)
        call mpas_pool_get_array(averagePool, 'avgNormalVelocity', avgNormalVelocity)
        call mpas_pool_get_array(averagePool, 'avgZonalVelocity', avgZonalVelocity)
        call mpas_pool_get_array(averagePool, 'avgMeridionalVelocity', avgMeridionalVelocity)
        call mpas_pool_get_array(averagePool, 'avgVertVelocityTop', avgVertVelocityTop)
        call mpas_pool_get_array(averagePool, 'varNormalVelocity', varNormalVelocity)
        call mpas_pool_get_array(averagePool, 'varZonalVelocity', varZonalVelocity)
        call mpas_pool_get_array(averagePool, 'varMeridionalVelocity', varMeridionalVelocity)
        call mpas_pool_get_array(averagePool, 'avgNormalTransportVelocity', avgNormalTransportVelocity)
        call mpas_pool_get_array(averagePool, 'avgZonalTransportVelocity', avgZonalTransportVelocity)
        call mpas_pool_get_array(averagePool, 'avgMeridionalTransportVelocity', avgMeridionalTransportVelocity)
        call mpas_pool_get_array(averagePool, 'avgVertTransportVelocityTop', avgVertTransportVelocityTop)
        call mpas_pool_get_array(averagePool, 'avgNormalGMBolusVelocity', avgNormalGMBolusVelocity)
        call mpas_pool_get_array(averagePool, 'avgZonalGMBolusVelocity', avgZonalGMBolusVelocity)
        call mpas_pool_get_array(averagePool, 'avgMeridionalGMBolusVelocity', avgMeridionalGMBolusVelocity)
        call mpas_pool_get_array(averagePool, 'avgVertGMBolusVelocityTop', avgVertGMBolusVelocityTop)

        avgSSH = avgSSH + ssh
        varSSH = varSSH + ssh**2
        avgNormalVelocity = avgNormalVelocity + normalVelocity
        avgZonalVelocity = avgZonalVelocity + normalVelocityZonal
        avgMeridionalVelocity = avgMeridionalVelocity + normalVelocityMeridional
        avgVertVelocityTop = avgVertVelocityTop + vertVelocityTop
        varNormalVelocity = varNormalVelocity + normalVelocity**2
        varZonalVelocity = varZonalVelocity + normalVelocityZonal**2
        varMeridionalVelocity = varMeridionalVelocity + normalVelocityMeridional**2
        avgNormalTransportVelocity = avgNormalTransportVelocity + normalTransportVelocity
        avgZonalTransportVelocity = avgZonalTransportVelocity + zonalTransportVelocity
        avgMeridionalTransportVelocity = avgMeridionalTransportVelocity + meridionalTransportVelocity
        avgVertTransportVelocityTop = avgVertTransportVelocityTop + vertTransportVelocityTop
        avgNormalGMBolusVelocity = avgNormalGMBolusVelocity + normalGMBolusVelocity
        avgZonalGMBolusVelocity = avgZonalGMBolusVelocity + zonalGMBolusVelocity
        avgMeridionalGMBolusVelocity = avgMeridionalGMBolusVelocity + meridionalGMBolusVelocity
        avgVertGMBolusVelocityTop = avgVertGMBolusVelocityTop + vertGMBolusVelocityTop

        nAverage = nAverage + 1
    end subroutine ocn_time_average_accumulate!}}}

    subroutine ocn_time_average_normalize(averagePool)!{{{
        type (mpas_pool_type), intent(inout) :: averagePool

        real (kind=RKIND), pointer :: nAverage

        real (kind=RKIND), dimension(:), pointer :: avgSSH, varSSH
        real (kind=RKIND), dimension(:,:), pointer :: &
           avgNormalVelocity, avgZonalVelocity, avgMeridionalVelocity, avgVertVelocityTop, &
           varNormalVelocity, varZonalVelocity, varMeridionalVelocity, &
           avgNormalTransportVelocity, avgZonalTransportVelocity, avgMeridionalTransportVelocity, avgVertTransportVelocityTop, &
           avgNormalGMBolusVelocity, avgZonalGMBolusVelocity, avgMeridionalGMBolusVelocity, avgVertGMBolusVelocityTop

        call mpas_pool_get_array(averagePool, 'nAverage', nAverage)
        call mpas_pool_get_array(averagePool, 'avgSsh', avgSSH)
        call mpas_pool_get_array(averagePool, 'varSsh', varSSH)
        call mpas_pool_get_array(averagePool, 'avgNormalVelocity', avgNormalVelocity)
        call mpas_pool_get_array(averagePool, 'avgZonalVelocity', avgZonalVelocity)
        call mpas_pool_get_array(averagePool, 'avgMeridionalVelocity', avgMeridionalVelocity)
        call mpas_pool_get_array(averagePool, 'avgVertVelocityTop', avgVertVelocityTop)
        call mpas_pool_get_array(averagePool, 'varNormalVelocity', varNormalVelocity)
        call mpas_pool_get_array(averagePool, 'varZonalVelocity', varZonalVelocity)
        call mpas_pool_get_array(averagePool, 'varMeridionalVelocity', varMeridionalVelocity)
        call mpas_pool_get_array(averagePool, 'avgNormalTransportVelocity', avgNormalTransportVelocity)
        call mpas_pool_get_array(averagePool, 'avgZonalTransportVelocity', avgZonalTransportVelocity)
        call mpas_pool_get_array(averagePool, 'avgMeridionalTransportVelocity', avgMeridionalTransportVelocity)
        call mpas_pool_get_array(averagePool, 'avgVertTransportVelocityTop', avgVertTransportVelocityTop)
        call mpas_pool_get_array(averagePool, 'avgNormalGMBolusVelocity', avgNormalGMBolusVelocity)
        call mpas_pool_get_array(averagePool, 'avgZonalGMBolusVelocity', avgZonalGMBolusVelocity)
        call mpas_pool_get_array(averagePool, 'avgMeridionalGMBolusVelocity', avgMeridionalGMBolusVelocity)
        call mpas_pool_get_array(averagePool, 'avgVertGMBolusVelocityTop', avgVertGMBolusVelocityTop)

        if(nAverage > 0) then
          avgSSH = avgSSH / nAverage
          varSSH = varSSH / nAverage
          avgNormalVelocity = avgNormalVelocity / nAverage
          avgZonalVelocity = avgZonalVelocity / nAverage
          avgMeridionalVelocity = avgMeridionalVelocity / nAverage
          avgVertVelocityTop = avgVertVelocityTop / nAverage
          varNormalVelocity = varNormalVelocity / nAverage
          varZonalVelocity = varZonalVelocity / nAverage
          varMeridionalVelocity = varMeridionalVelocity / nAverage
          avgNormalTransportVelocity = avgNormalTransportVelocity / nAverage
          avgZonalTransportVelocity = avgZonalTransportVelocity / nAverage
          avgMeridionalTransportVelocity = avgMeridionalTransportVelocity / nAverage
          avgVertTransportVelocityTop = avgVertTransportVelocityTop / nAverage
          avgNormalGMBolusVelocity = avgNormalGMBolusVelocity / nAverage
          avgZonalGMBolusVelocity = avgZonalGMBolusVelocity / nAverage
          avgMeridionalGMBolusVelocity = avgMeridionalGMBolusVelocity / nAverage
          avgVertGMBolusVelocityTop = avgVertGMBolusVelocityTop / nAverage
        end if
    end subroutine ocn_time_average_normalize!}}}

end module ocn_time_average
