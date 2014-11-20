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
           avgNormalVelocity, avgVelocityZonal, avgVelocityMeridional, avgVertVelocityTop, &
           varNormalVelocity, varVelocityZonal, varVelocityMeridional, &
           avgNormalTransportVelocity, avgTransportVelocityZonal, avgTransportVelocityMeridional, avgVertTransportVelocityTop, &
           avgNormalGMBolusVelocity, avgGMBolusVelocityZonal, avgGMBolusVelocityMeridional, avgVertGMBolusVelocityTop

        call mpas_pool_get_array(averagePool, 'nAverage', nAverage)
        call mpas_pool_get_array(averagePool, 'avgSSH', avgSSH)
        call mpas_pool_get_array(averagePool, 'varSSH', varSSH)
        call mpas_pool_get_array(averagePool, 'avgNormalVelocity', avgNormalVelocity)
        call mpas_pool_get_array(averagePool, 'avgVelocityZonal', avgVelocityZonal)
        call mpas_pool_get_array(averagePool, 'avgVelocityMeridional', avgVelocityMeridional)
        call mpas_pool_get_array(averagePool, 'avgVertVelocityTop', avgVertVelocityTop)
        call mpas_pool_get_array(averagePool, 'varNormalVelocity', varNormalVelocity)
        call mpas_pool_get_array(averagePool, 'varVelocityZonal', varVelocityZonal)
        call mpas_pool_get_array(averagePool, 'varVelocityMeridional', varVelocityMeridional)
        call mpas_pool_get_array(averagePool, 'avgNormalTransportVelocity', avgNormalTransportVelocity)
        call mpas_pool_get_array(averagePool, 'avgTransportVelocityZonal', avgTransportVelocityZonal)
        call mpas_pool_get_array(averagePool, 'avgTransportVelocityMeridional', avgTransportVelocityMeridional)
        call mpas_pool_get_array(averagePool, 'avgVertTransportVelocityTop', avgVertTransportVelocityTop)
        call mpas_pool_get_array(averagePool, 'avgNormalGMBolusVelocity', avgNormalGMBolusVelocity)
        call mpas_pool_get_array(averagePool, 'avgGMBolusVelocityZonal', avgGMBolusVelocityZonal)
        call mpas_pool_get_array(averagePool, 'avgGMBolusVelocityMeridional', avgGMBolusVelocityMeridional)
        call mpas_pool_get_array(averagePool, 'avgVertGMBolusVelocityTop', avgVertGMBolusVelocityTop)

        nAverage = 0

        avgSSH = 0.0
        varSSH = 0.0
        avgNormalVelocity = 0.0
        avgVelocityZonal = 0.0
        avgVelocityMeridional = 0.0
        avgVertVelocityTop = 0.0
        varNormalVelocity = 0.0
        varVelocityZonal = 0.0
        varVelocityMeridional = 0.0
        avgNormalTransportVelocity = 0.0
        avgTransportVelocityZonal = 0.0
        avgTransportVelocityMeridional = 0.0
        avgVertTransportVelocityTop = 0.0
        avgNormalGMBolusVelocity = 0.0
        avgGMBolusVelocityZonal = 0.0
        avgGMBolusVelocityMeridional = 0.0
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
           velocityZonal, velocityMeridional, normalVelocity, vertVelocityTop, &
           transportVelocityZonal, transportVelocityMeridional, normalTransportVelocity, vertTransportVelocityTop, &
           GMBolusVelocityZonal, GMBolusVelocityMeridional, normalGMBolusVelocity, vertGMBolusVelocityTop

        real (kind=RKIND), dimension(:), pointer :: avgSSH, varSSH
        real (kind=RKIND), dimension(:,:), pointer :: &
           avgNormalVelocity, avgVelocityZonal, avgVelocityMeridional, avgVertVelocityTop, &
           varNormalVelocity, varVelocityZonal, varVelocityMeridional, &
           avgNormalTransportVelocity, avgTransportVelocityZonal, avgTransportVelocityMeridional, avgVertTransportVelocityTop, &
           avgNormalGMBolusVelocity, avgGMBolusVelocityZonal, avgGMBolusVelocityMeridional, avgVertGMBolusVelocityTop

        real (kind=RKIND), dimension(:), pointer :: old_avgSSH, old_varSSH
        real (kind=RKIND), dimension(:,:), pointer :: &
           old_avgNormalVelocity, old_avgVelocityZonal, old_avgVelocityMeridional, old_avgVertVelocityTop, &
           old_varNormalVelocity, old_varVelocityZonal, old_varVelocityMeridional, &
           old_avgNormalTransportVelocity, old_avgTransportVelocityZonal, old_avgTransportVelocityMeridional, old_avgVertTransportVelocityTop, &
           old_avgNormalGMBolusVelocity, old_avgGMBolusVelocityZonal, old_avgGMBolusVelocityMeridional, old_avgVertGMBolusVelocityTop

        integer :: timeLevel

        if (present(timeLevelIn)) then
           timeLevel = timeLevelIn
        else
           timeLevel = 1
        end if

        call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)
        call mpas_pool_get_array(statePool, 'ssh', ssh, timeLevel)

        call mpas_pool_get_array(diagnosticsPool, 'velocityZonal', velocityZonal)
        call mpas_pool_get_array(diagnosticsPool, 'velocityMeridional', velocityMeridional)
        call mpas_pool_get_array(diagnosticsPool, 'vertVelocityTop', vertVelocityTop)
        call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity ', normalTransportVelocity)
        call mpas_pool_get_array(diagnosticsPool, 'transportVelocityZonal', transportVelocityZonal)
        call mpas_pool_get_array(diagnosticsPool, 'transportVelocityMeridional', transportVelocityMeridional)
        call mpas_pool_get_array(diagnosticsPool, 'vertTransportVelocityTop', vertTransportVelocityTop)
        call mpas_pool_get_array(diagnosticsPool, 'normalGMBolusVelocity', normalGMBolusVelocity)
        call mpas_pool_get_array(diagnosticsPool, 'GMBolusVelocityZonal', GMBolusVelocityZonal)
        call mpas_pool_get_array(diagnosticsPool, 'GMBolusVelocityMeridional', GMBolusVelocityMeridional)
        call mpas_pool_get_array(diagnosticsPool, 'vertGMBolusVelocityTop', vertGMBolusVelocityTop)

        call mpas_pool_get_array(averagePool, 'nAverage', nAverage)
        call mpas_pool_get_array(averagePool, 'avgSSH', avgSSH)
        call mpas_pool_get_array(averagePool, 'varSSH', varSSH)
        call mpas_pool_get_array(averagePool, 'avgNormalVelocity', avgNormalVelocity)
        call mpas_pool_get_array(averagePool, 'avgVelocityZonal', avgVelocityZonal)
        call mpas_pool_get_array(averagePool, 'avgVelocityMeridional', avgVelocityMeridional)
        call mpas_pool_get_array(averagePool, 'avgVertVelocityTop', avgVertVelocityTop)
        call mpas_pool_get_array(averagePool, 'varNormalVelocity', varNormalVelocity)
        call mpas_pool_get_array(averagePool, 'varVelocityZonal', varVelocityZonal)
        call mpas_pool_get_array(averagePool, 'varVelocityMeridional', varVelocityMeridional)
        call mpas_pool_get_array(averagePool, 'avgNormalTransportVelocity', avgNormalTransportVelocity)
        call mpas_pool_get_array(averagePool, 'avgTransportVelocityZonal', avgTransportVelocityZonal)
        call mpas_pool_get_array(averagePool, 'avgTransportVelocityMeridional', avgTransportVelocityMeridional)
        call mpas_pool_get_array(averagePool, 'avgVertTransportVelocityTop', avgVertTransportVelocityTop)
        call mpas_pool_get_array(averagePool, 'avgNormalGMBolusVelocity', avgNormalGMBolusVelocity)
        call mpas_pool_get_array(averagePool, 'avgGMBolusVelocityZonal', avgGMBolusVelocityZonal)
        call mpas_pool_get_array(averagePool, 'avgGMBolusVelocityMeridional', avgGMBolusVelocityMeridional)
        call mpas_pool_get_array(averagePool, 'avgVertGMBolusVelocityTop', avgVertGMBolusVelocityTop)

        avgSSH = avgSSH + ssh
        varSSH = varSSH + ssh**2
        avgNormalVelocity = avgNormalVelocity + normalVelocity
        avgVelocityZonal = avgVelocityZonal + velocityZonal
        avgVelocityMeridional = avgVelocityMeridional + velocityMeridional
        avgVertVelocityTop = avgVertVelocityTop + vertVelocityTop
        varNormalVelocity = varNormalVelocity + normalVelocity**2
        varVelocityZonal = varVelocityZonal + velocityZonal**2
        varVelocityMeridional = varVelocityMeridional + velocityMeridional**2
        avgNormalTransportVelocity = avgNormalTransportVelocity + normalTransportVelocity
        avgTransportVelocityZonal = avgTransportVelocityZonal + transportVelocityZonal
        avgTransportVelocityMeridional = avgTransportVelocityMeridional + transportVelocityMeridional
        avgVertTransportVelocityTop = avgVertTransportVelocityTop + vertTransportVelocityTop
        avgNormalGMBolusVelocity = avgNormalGMBolusVelocity + normalGMBolusVelocity
        avgGMBolusVelocityZonal = avgGMBolusVelocityZonal + GMBolusVelocityZonal
        avgGMBolusVelocityMeridional = avgGMBolusVelocityMeridional + GMBolusVelocityMeridional
        avgVertGMBolusVelocityTop = avgVertGMBolusVelocityTop + vertGMBolusVelocityTop

        nAverage = nAverage + 1
    end subroutine ocn_time_average_accumulate!}}}

    subroutine ocn_time_average_normalize(averagePool)!{{{
        type (mpas_pool_type), intent(inout) :: averagePool

        real (kind=RKIND), pointer :: nAverage

        real (kind=RKIND), dimension(:), pointer :: avgSSH, varSSH
        real (kind=RKIND), dimension(:,:), pointer :: &
           avgNormalVelocity, avgVelocityZonal, avgVelocityMeridional, avgVertVelocityTop, &
           varNormalVelocity, varVelocityZonal, varVelocityMeridional, &
           avgNormalTransportVelocity, avgTransportVelocityZonal, avgTransportVelocityMeridional, avgVertTransportVelocityTop, &
           avgNormalGMBolusVelocity, avgGMBolusVelocityZonal, avgGMBolusVelocityMeridional, avgVertGMBolusVelocityTop

        call mpas_pool_get_array(averagePool, 'nAverage', nAverage)
        call mpas_pool_get_array(averagePool, 'avgSSH', avgSSH)
        call mpas_pool_get_array(averagePool, 'varSSH', varSSH)
        call mpas_pool_get_array(averagePool, 'avgNormalVelocity', avgNormalVelocity)
        call mpas_pool_get_array(averagePool, 'avgVelocityZonal', avgVelocityZonal)
        call mpas_pool_get_array(averagePool, 'avgVelocityMeridional', avgVelocityMeridional)
        call mpas_pool_get_array(averagePool, 'avgVertVelocityTop', avgVertVelocityTop)
        call mpas_pool_get_array(averagePool, 'varNormalVelocity', varNormalVelocity)
        call mpas_pool_get_array(averagePool, 'varVelocityZonal', varVelocityZonal)
        call mpas_pool_get_array(averagePool, 'varVelocityMeridional', varVelocityMeridional)
        call mpas_pool_get_array(averagePool, 'avgNormalTransportVelocity', avgNormalTransportVelocity)
        call mpas_pool_get_array(averagePool, 'avgTransportVelocityZonal', avgTransportVelocityZonal)
        call mpas_pool_get_array(averagePool, 'avgTransportVelocityMeridional', avgTransportVelocityMeridional)
        call mpas_pool_get_array(averagePool, 'avgVertTransportVelocityTop', avgVertTransportVelocityTop)
        call mpas_pool_get_array(averagePool, 'avgNormalGMBolusVelocity', avgNormalGMBolusVelocity)
        call mpas_pool_get_array(averagePool, 'avgGMBolusVelocityZonal', avgGMBolusVelocityZonal)
        call mpas_pool_get_array(averagePool, 'avgGMBolusVelocityMeridional', avgGMBolusVelocityMeridional)
        call mpas_pool_get_array(averagePool, 'avgVertGMBolusVelocityTop', avgVertGMBolusVelocityTop)

        if(nAverage > 0) then
          avgSSH = avgSSH / nAverage
          varSSH = varSSH / nAverage
          avgNormalVelocity = avgNormalVelocity / nAverage
          avgVelocityZonal = avgVelocityZonal / nAverage
          avgVelocityMeridional = avgVelocityMeridional / nAverage
          avgVertVelocityTop = avgVertVelocityTop / nAverage
          varNormalVelocity = varNormalVelocity / nAverage
          varVelocityZonal = varVelocityZonal / nAverage
          varVelocityMeridional = varVelocityMeridional / nAverage
          avgNormalTransportVelocity = avgNormalTransportVelocity / nAverage
          avgTransportVelocityZonal = avgTransportVelocityZonal / nAverage
          avgTransportVelocityMeridional = avgTransportVelocityMeridional / nAverage
          avgVertTransportVelocityTop = avgVertTransportVelocityTop / nAverage
          avgNormalGMBolusVelocity = avgNormalGMBolusVelocity / nAverage
          avgGMBolusVelocityZonal = avgGMBolusVelocityZonal / nAverage
          avgGMBolusVelocityMeridional = avgGMBolusVelocityMeridional / nAverage
          avgVertGMBolusVelocityTop = avgVertGMBolusVelocityTop / nAverage
        end if
    end subroutine ocn_time_average_normalize!}}}

end module ocn_time_average
