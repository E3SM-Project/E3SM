! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
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

    subroutine ocn_time_average_init(state)!{{{
        type (state_type), intent(inout) :: state

        real (kind=RKIND), pointer :: nAverage

        real (kind=RKIND), dimension(:), pointer :: avgSSH, varSSH
        real (kind=RKIND), dimension(:,:), pointer :: avgnormalVelocityZonal, avgnormalVelocityMeridional, varnormalVelocityZonal, varnormalVelocityMeridional
        real (kind=RKIND), dimension(:,:), pointer :: avgNormalVelocity, varNormalVelocity, avgVertVelocityTop

        nAverage => state % nAverage % scalar

        avgSSH => state % avgSSH % array
        varSSH => state % varSSH % array
        avgnormalVelocityZonal => state % avgnormalVelocityZonal % array
        avgnormalVelocityMeridional => state % avgnormalVelocityMeridional % array
        varnormalVelocityZonal => state % varnormalVelocityZonal % array
        varnormalVelocityMeridional => state % varnormalVelocityMeridional % array
        avgNormalVelocity => state % avgNormalVelocity % array
        varNormalVelocity => state % varNormalVelocity % array
        avgVertVelocityTop => state % avgVertVelocityTop % array

        nAverage = 0

        avgSSH = 0.0
        varSSH = 0.0
        avgnormalVelocityZonal = 0.0
        avgnormalVelocityMeridional = 0.0
        varnormalVelocityZonal = 0.0
        varnormalVelocityMeridional = 0.0
        avgNormalVelocity = 0.0
        varNormalVelocity = 0.0
        avgVertVelocityTop = 0.0

    end subroutine ocn_time_average_init!}}}

    subroutine ocn_time_average_accumulate(state, old_state)!{{{
        type (state_type), intent(inout) :: state
        type (state_type), intent(in) :: old_state

        real (kind=RKIND), pointer :: nAverage, old_nAverage

        real (kind=RKIND), dimension(:), pointer :: ssh
        real (kind=RKIND), dimension(:,:), pointer :: normalVelocityZonal, normalVelocityMeridional, normalVelocity, vertVelocityTop

        real (kind=RKIND), dimension(:,:), pointer :: avgNormalVelocity, varNormalVelocity, avgVertVelocityTop
        real (kind=RKIND), dimension(:,:), pointer :: avgnormalVelocityZonal, avgnormalVelocityMeridional, varnormalVelocityZonal, varnormalVelocityMeridional
        real (kind=RKIND), dimension(:), pointer :: avgSSH, varSSH

        real (kind=RKIND), dimension(:,:), pointer :: old_avgNormalVelocity, old_varNormalVelocity, old_avgVertVelocityTop
        real (kind=RKIND), dimension(:,:), pointer :: old_avgnormalVelocityZonal, old_avgnormalVelocityMeridional, old_varnormalVelocityZonal, old_varnormalVelocityMeridional
        real (kind=RKIND), dimension(:), pointer :: old_avgSSH, old_varSSH

        old_nAverage => old_state % nAverage  % scalar
        nAverage => state % nAverage  % scalar

        ssh => state % ssh % array
        normalVelocityZonal => state % normalVelocityZonal % array
        normalVelocityMeridional => state % normalVelocityMeridional % array
        normalVelocity => state % normalVelocity % array
        vertVelocityTop => state % vertVelocityTop % array

        avgSSH => state % avgSSH % array
        varSSH => state % varSSH % array
        avgnormalVelocityZonal => state % avgnormalVelocityZonal % array
        avgnormalVelocityMeridional => state % avgnormalVelocityMeridional % array
        varnormalVelocityZonal => state % varnormalVelocityZonal % array
        varnormalVelocityMeridional => state % varnormalVelocityMeridional % array
        avgNormalVelocity => state % avgNormalVelocity % array
        varNormalVelocity => state % varNormalVelocity % array
        avgVertVelocityTop => state % avgVertVelocityTop % array

        old_avgSSH => old_state % avgSSH % array
        old_varSSH => old_state % varSSH % array
        old_avgnormalVelocityZonal => old_state % avgnormalVelocityZonal % array
        old_avgnormalVelocityMeridional => old_state % avgnormalVelocityMeridional % array
        old_varnormalVelocityZonal => old_state % varnormalVelocityZonal % array
        old_varnormalVelocityMeridional => old_state % varnormalVelocityMeridional % array
        old_avgNormalVelocity => old_state % avgNormalVelocity % array
        old_varNormalVelocity => old_state % varNormalVelocity % array
        old_avgVertVelocityTop => old_state % avgVertVelocityTop % array

        avgSSH = old_avgSSH + ssh
        varSSH = old_varSSH + ssh**2
        avgnormalVelocityZonal = old_avgnormalVelocityZonal + normalVelocityZonal
        avgnormalVelocityMeridional = old_avgnormalVelocityMeridional + normalVelocityMeridional
        varnormalVelocityZonal = old_varnormalVelocityZonal + normalVelocityZonal**2
        varnormalVelocityMeridional = old_varnormalVelocityMeridional + normalVelocityMeridional**2
        avgNormalVelocity = old_avgNormalVelocity + normalVelocity
        varNormalVelocity = old_varNormalVelocity + normalVelocity**2
        avgVertVelocityTop = old_avgVertVelocityTop + vertVelocityTop

        nAverage = old_nAverage + 1
    end subroutine ocn_time_average_accumulate!}}}

    subroutine ocn_time_average_normalize(state)!{{{
        type (state_type), intent(inout) :: state

        real (kind=RKIND), pointer :: nAverage

        real (kind=RKIND), dimension(:), pointer :: avgSSH, varSSH
        real (kind=RKIND), dimension(:,:), pointer :: avgnormalVelocityZonal, avgnormalVelocityMeridional, varnormalVelocityZonal, varnormalVelocityMeridional
        real (kind=RKIND), dimension(:,:), pointer :: avgNormalVelocity, varNormalVelocity, avgVertVelocityTop

        nAverage => state % nAverage  % scalar

        avgSSH => state % avgSSH % array
        varSSH => state % varSSH % array
        avgnormalVelocityZonal => state % avgnormalVelocityZonal % array
        avgnormalVelocityMeridional => state % avgnormalVelocityMeridional % array
        varnormalVelocityZonal => state % varnormalVelocityZonal % array
        varnormalVelocityMeridional => state % varnormalVelocityMeridional % array
        avgNormalVelocity => state % avgNormalVelocity % array
        varNormalVelocity => state % varNormalVelocity % array
        avgVertVelocityTop => state % avgVertVelocityTop % array

        if(nAverage > 0) then
          avgSSH = avgSSH / nAverage
          varSSH = varSSH / nAverage
          avgnormalVelocityZonal = avgnormalVelocityZonal / nAverage
          avgnormalVelocityMeridional = avgnormalVelocityMeridional / nAverage
          varnormalVelocityZonal = varnormalVelocityZonal / nAverage
          varnormalVelocityMeridional = varnormalVelocityMeridional / nAverage
          avgNormalVelocity = avgNormalVelocity / nAverage
          varNormalVelocity = varNormalVelocity / nAverage
          avgVertVelocityTop = avgVertVelocityTop / nAverage
        end if
    end subroutine ocn_time_average_normalize!}}}

end module ocn_time_average
