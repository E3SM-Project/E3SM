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
        real (kind=RKIND), dimension(:,:), pointer :: avgVelocityZonal, avgVelocityMeridional, varVelocityZonal, varVelocityMeridional
        real (kind=RKIND), dimension(:,:), pointer :: avgNormalVelocity, varNormalVelocity, avgVertVelocityTop

        call mpas_pool_get_array(averagePool, 'nAverage', nAverage)

        call mpas_pool_get_array(averagePool, 'avgSsh', avgSSH)
        call mpas_pool_get_array(averagePool, 'varSsh', varSSH)
        call mpas_pool_get_array(averagePool, 'avgVelocityZonal', avgVelocityZonal)
        call mpas_pool_get_array(averagePool, 'avgVelocityMeridional', avgVelocityMeridional)
        call mpas_pool_get_array(averagePool, 'varVelocityZonal', varVelocityZonal)
        call mpas_pool_get_array(averagePool, 'varVelocityMeridional', varVelocityMeridional)
        call mpas_pool_get_array(averagePool, 'avgNormalVelocity', avgNormalVelocity)
        call mpas_pool_get_array(averagePool, 'varNormalVelocity', varNormalVelocity)
        call mpas_pool_get_array(averagePool, 'avgVertVelocityTop', avgVertVelocityTop)

        nAverage = 0

        avgSSH = 0.0
        varSSH = 0.0
        avgVelocityZonal = 0.0
        avgVelocityMeridional = 0.0
        varVelocityZonal = 0.0
        varVelocityMeridional = 0.0
        avgNormalVelocity = 0.0
        varNormalVelocity = 0.0
        avgVertVelocityTop = 0.0

    end subroutine ocn_time_average_init!}}}

    subroutine ocn_time_average_accumulate(averagePool, statePool, diagnosticsPool, timeLevelIn)!{{{
        type (mpas_pool_type), intent(inout) :: averagePool
        type (mpas_pool_type), intent(in) :: statePool
        type (mpas_pool_type), intent(in) :: diagnosticsPool
        integer, intent(in), optional :: timeLevelIn

        real (kind=RKIND), pointer :: nAverage, old_nAverage

        real (kind=RKIND), dimension(:), pointer :: ssh
        real (kind=RKIND), dimension(:,:), pointer :: normalVelocityZonal, normalVelocityMeridional, normalVelocity, vertVelocityTop

        real (kind=RKIND), dimension(:,:), pointer :: avgNormalVelocity, varNormalVelocity, avgVertVelocityTop
        real (kind=RKIND), dimension(:,:), pointer :: avgVelocityZonal, avgVelocityMeridional, varVelocityZonal, varVelocityMeridional
        real (kind=RKIND), dimension(:), pointer :: avgSSH, varSSH

        real (kind=RKIND), dimension(:,:), pointer :: old_avgNormalVelocity, old_varNormalVelocity, old_avgVertVelocityTop
        real (kind=RKIND), dimension(:,:), pointer :: old_avgVelocityZonal, old_avgVelocityMeridional, old_varVelocityZonal, old_varVelocityMeridional
        real (kind=RKIND), dimension(:), pointer :: old_avgSSH, old_varSSH

        integer :: timeLevel

        if (present(timeLevelIn)) then
           timeLevel = timeLevelIn
        else
           timeLevel = 1
        end if

        call mpas_pool_get_array(averagePool, 'nAverage ', nAverage)

        call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)
        call mpas_pool_get_array(statePool, 'ssh', ssh, timeLevel)

        call mpas_pool_get_array(diagnosticsPool, 'normalVelocityZonal', normalVelocityZonal)
        call mpas_pool_get_array(diagnosticsPool, 'normalVelocityMeridional', normalVelocityMeridional)
        call mpas_pool_get_array(diagnosticsPool, 'vertVelocityTop', vertVelocityTop)

        call mpas_pool_get_array(averagePool, 'avgSsh', avgSSH)
        call mpas_pool_get_array(averagePool, 'varSsh', varSSH)
        call mpas_pool_get_array(averagePool, 'avgVelocityZonal', avgVelocityZonal)
        call mpas_pool_get_array(averagePool, 'avgVelocityMeridional', avgVelocityMeridional)
        call mpas_pool_get_array(averagePool, 'varVelocityZonal', varVelocityZonal)
        call mpas_pool_get_array(averagePool, 'varVelocityMeridional', varVelocityMeridional)
        call mpas_pool_get_array(averagePool, 'avgNormalVelocity', avgNormalVelocity)
        call mpas_pool_get_array(averagePool, 'varNormalVelocity', varNormalVelocity)
        call mpas_pool_get_array(averagePool, 'avgVertVelocityTop', avgVertVelocityTop)

        avgSSH = avgSSH + ssh
        varSSH = varSSH + ssh**2
        avgVelocityZonal = avgVelocityZonal + normalVelocityZonal
        avgVelocityMeridional = avgVelocityMeridional + normalVelocityMeridional
        varVelocityZonal = varVelocityZonal + normalVelocityZonal**2
        varVelocityMeridional = varVelocityMeridional + normalVelocityMeridional**2
        avgNormalVelocity = avgNormalVelocity + normalVelocity
        varNormalVelocity = varNormalVelocity + normalVelocity**2
        avgVertVelocityTop = avgVertVelocityTop + vertVelocityTop

        nAverage = nAverage + 1
    end subroutine ocn_time_average_accumulate!}}}

    subroutine ocn_time_average_normalize(averagePool)!{{{
        type (mpas_pool_type), intent(inout) :: averagePool

        real (kind=RKIND), pointer :: nAverage

        real (kind=RKIND), dimension(:), pointer :: avgSSH, varSSH
        real (kind=RKIND), dimension(:,:), pointer :: avgVelocityZonal, avgVelocityMeridional, varVelocityZonal, varVelocityMeridional
        real (kind=RKIND), dimension(:,:), pointer :: avgNormalVelocity, varNormalVelocity, avgVertVelocityTop

        call mpas_pool_get_array(averagePool, 'nAverage ', nAverage)

        call mpas_pool_get_array(averagePool, 'avgSsh', avgSSH)
        call mpas_pool_get_array(averagePool, 'varSsh', varSSH)
        call mpas_pool_get_array(averagePool, 'avgVelocityZonal', avgVelocityZonal)
        call mpas_pool_get_array(averagePool, 'avgVelocityMeridional', avgVelocityMeridional)
        call mpas_pool_get_array(averagePool, 'varVelocityZonal', varVelocityZonal)
        call mpas_pool_get_array(averagePool, 'varVelocityMeridional', varVelocityMeridional)
        call mpas_pool_get_array(averagePool, 'avgNormalVelocity', avgNormalVelocity)
        call mpas_pool_get_array(averagePool, 'varNormalVelocity', varNormalVelocity)
        call mpas_pool_get_array(averagePool, 'avgVertVelocityTop', avgVertVelocityTop)

        if(nAverage > 0) then
          avgSSH = avgSSH / nAverage
          varSSH = varSSH / nAverage
          avgVelocityZonal = avgVelocityZonal / nAverage
          avgVelocityMeridional = avgVelocityMeridional / nAverage
          varVelocityZonal = varVelocityZonal / nAverage
          varVelocityMeridional = varVelocityMeridional / nAverage
          avgNormalVelocity = avgNormalVelocity / nAverage
          varNormalVelocity = varNormalVelocity / nAverage
          avgVertVelocityTop = avgVertVelocityTop / nAverage

          nAverage = 0
        end if
    end subroutine ocn_time_average_normalize!}}}

end module ocn_time_average
