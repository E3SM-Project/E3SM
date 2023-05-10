!
MODULE MOSART_reservoir_mod
! Description: core code of MOSART-reservoir. Can be incoporated within any river model via a interface module
! 
! Developed by Hongyi Li, 11/2016.
! REVISION HISTORY:
!-----------------------------------------------------------------------

! !USES:
    use shr_kind_mod  , only : r8 => shr_kind_r8, SHR_KIND_CL
    use shr_const_mod , only : denh2o => SHR_CONST_RHOFW, denice => SHR_CONST_RHOICE, grav => SHR_CONST_G, SHR_CONST_REARTH, SHR_CONST_PI
    !use clm_varcon , only : denh2o, grav !!density of liquid water [kg/m3], gravity constant [m/s2]
    use RunoffMod, only : Tctl, TUnit, TRunoff, TPara
    use RunoffMod, only : rtmCTL
    use rof_cpl_indices, only : nt_rtm, rtm_tracers, nt_nliq, nt_nice, nt_nmud, nt_nsan
    use MOSART_RES_type, only : Tres_para, Tres
    implicit none
    real(r8), parameter :: TINYVALUE_r = 1.0e-14_r8  ! double precision variable has a significance of about 16 decimal digits

! !PUBLIC MEMBER FUNCTIONS:
    contains


    subroutine res_trapping(iunit, nt)
    ! !DESCRIPTION: trapping of particulate fluxes in the reservoirs on the main channels (e.g., sediment)
    ! !assume the trapping occurs after the channel routing, i.e., won't affect the current channel mass balance, but only change the outflow to the downstream
        implicit none

        integer, intent(in) :: iunit, nt
        !real(r8), intent(in) :: theDeltaT

        real(r8) :: trapping_eff  ! trapping efficiency, proportion
        trapping_eff = Tres_para%Eff_trapping(iunit)

        Tres%eres_in(iunit,nt) = -TRunoff%erout(iunit,nt)  * trapping_eff
        Tres%eres_out(iunit,nt) = 0._r8  ! currently assuming no sediment release from the dam bottom (this might be OK in U.S., but wrong in other plances such as China)
        Tres%dwres(iunit,nt) =  Tres%eres_in(iunit,nt) + Tres%eres_out(iunit,nt)

        TRunoff%erout(iunit,nt) = TRunoff%erout(iunit,nt) * (1._r8 - trapping_eff)

    end subroutine res_trapping

    subroutine res_trapping_t(iunit, nt)
    ! !DESCRIPTION: trapping of particulate fluxes in the reservoirs on the sub-network channels (e.g., sediment)
    ! !assume the trapping occurs after the channel routing, i.e., won't affect the current channel mass balance, but only change the outflow to the downstream
        implicit none

        integer, intent(in) :: iunit, nt
        !real(r8), intent(in) :: theDeltaT

        real(r8) :: trapping_eff_t  ! trapping efficiency, propoation
        trapping_eff_t = Tres_para%Eff_trapping_t(iunit)

        Tres%eres_in_t(iunit,nt) = -TRunoff%erlateral(iunit,nt)  * trapping_eff_t
        Tres%eres_out_t(iunit,nt) = 0._r8  ! currently assuming no sediment release from the dam bottom (this might be OK in U.S., but wrong in other plances such as China)
        Tres%dwres_t(iunit,nt) =  Tres%eres_in_t(iunit,nt) + Tres%eres_out_t(iunit,nt)

        TRunoff%erlateral(iunit,nt) = TRunoff%erlateral(iunit,nt) * (1._r8 - trapping_eff_t)

    end subroutine res_trapping_t

    subroutine res_trapping_r(iunit, nt)
    ! !DESCRIPTION: trapping of particulate fluxes in the reservoirs on the main channels (e.g., sediment)
    ! !assume the trapping occurs after the channel routing, i.e., won't affect the current channel mass balance, but only change the outflow to the downstream
        implicit none

        integer, intent(in) :: iunit, nt
        !real(r8), intent(in) :: theDeltaT

        real(r8) :: trapping_eff_r  ! trapping efficiency, proportion
        trapping_eff_r = Tres_para%Eff_trapping_r(iunit)

        Tres%eres_in(iunit,nt) = -TRunoff%erout(iunit,nt)  * trapping_eff_r
        Tres%eres_out(iunit,nt) = 0._r8  ! currently assuming no sediment release from the dam bottom (this might be OK in U.S., but wrong in other plances such as China)
        Tres%dwres(iunit,nt) =  Tres%eres_in(iunit,nt) + Tres%eres_out(iunit,nt)

        TRunoff%erout(iunit,nt) = TRunoff%erout(iunit,nt) * (1._r8 - trapping_eff_r)

    end subroutine res_trapping_r

end MODULE MOSART_reservoir_mod 