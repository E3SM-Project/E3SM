!
MODULE MOSART_heat_mod
! Description: core code of MOSART-heat. 
! 
! Developed by Hongyi Li, 12/29/2015. 
! REVISION HISTORY:
! April 2019, migrated to E3SM
!-----------------------------------------------------------------------
    
! !USES:
    use shr_kind_mod  , only : r8 => shr_kind_r8, SHR_KIND_CL
    use shr_const_mod , only : SHR_CONST_REARTH, SHR_CONST_PI
    use shr_const_mod , only : cpliq => SHR_CONST_CPFW, denh2o => SHR_CONST_RHOFW, denice => SHR_CONST_RHOICE, sb => SHR_CONST_STEBOL
    use RunoffMod , only : Tctl, TUnit, TRunoff, THeat, TPara
    use RunoffMod , only : rtmCTL
    use rof_cpl_indices, only : nt_rtm, rtm_tracers, nt_nliq, nt_nice
    implicit none
    real(r8), parameter :: TINYVALUE1 = 1.0e-14_r8  ! double precision variable has a significance of about 16 decimal digits
    real(r8), parameter :: Hthreshold = 0.05_r8  ! threshold value of channel water depth, if lower than this, headwater temp is calculated using the simplified formula, instead of heat balance equation
    real(r8), parameter :: WaterAreaRatio = 0.8_r8  ! dimensionless coefficient for effective surface area
    real(r8), parameter :: WindSheltering = 0.8_r8  ! dimensionless coefficient for the wind shltering by riparian vegetation

! !PUBLIC MEMBER FUNCTIONS:
    contains
    
    subroutine hillslopeHeat(iunit, theDeltaT)
    ! !DESCRIPTION: calculate the temperature of surface and subsurface runoff.
        implicit none
        
        integer, intent(in) :: iunit
        real(r8), intent(in) :: theDeltaT        
            THeat%Tqsur(iunit) = THeat%Tqsur(iunit)
            ! adjust surface runoff temperature estimated based on the top-layer soil temperature, i.e., no less than freezing point
            if(THeat%Tqsur(iunit) < 273.15_r8-TINYVALUE1) then
                THeat%Tqsur(iunit) = 273.15_r8
            end if
            
            THeat%Tqsub(iunit) = THeat%Tqsub(iunit)
            if(THeat%Tqsub(iunit) < 273.15_r8-TINYVALUE1) then
                THeat%Tqsub(iunit) = 273.15_r8
            end if
        !end if
    end subroutine hillslopeHeat
    
    subroutine subnetworkHeat(iunit, theDeltaT)
    ! !DESCRIPTION: calculate the net heat balance of subnetwork channel.
    use shr_sys_mod , only : shr_sys_flush
        implicit none
        integer, intent(in) :: iunit
        real(r8), intent(in) :: theDeltaT    
        
        real(r8) :: Qsur, Qsub ! flow rate of surface and subsurface runoff separately
        !if(TUnit%fdir(iunit) >= 0 .and. TUnit%areaTotal(iunit) > TINYVALUE1 .and. TUnit%tlen(iunit) >= TINYVALUE1) then
                TRunoff%tarea(iunit,nt_nliq) = WaterAreaRatio*TUnit%twidth(iunit) * TUnit%tlen(iunit)
                THeat%Hs_t(iunit) = cr_swrad(THeat%forc_solar(iunit), TRunoff%tarea(iunit,nt_nliq))
                THeat%Hl_t(iunit) = cr_lwrad(THeat%forc_lwrad(iunit), THeat%Tt(iunit), TRunoff%tarea(iunit,nt_nliq))
                THeat%He_t(iunit) = cr_latentheat(THeat%forc_t(iunit), THeat%forc_pbot(iunit), THeat%forc_vp(iunit), THeat%forc_wind(iunit), WindSheltering, THeat%Tt(iunit), TRunoff%tarea(iunit,nt_nliq))
                THeat%Hh_t(iunit) = cr_sensibleheat(THeat%forc_t(iunit), THeat%forc_pbot(iunit), THeat%forc_wind(iunit), WindSheltering, THeat%Tt(iunit), TRunoff%tarea(iunit,nt_nliq))
                THeat%Hc_t(iunit) = cr_condheat(THeat%Hs_t(iunit),THeat%Hl_t(iunit),TRunoff%tarea(iunit,nt_nliq))

                Qsur = (-TRunoff%ehout(iunit,nt_nliq)-TRunoff%ehout(iunit,nt_nice)) * TUnit%area(iunit) * TUnit%frac(iunit)
                Qsub = TRunoff%qsub(iunit,nt_nliq) * TUnit%area(iunit) * TUnit%frac(iunit)
                THeat%Ha_h2t(iunit) = cr_advectheat(Qsur, THeat%Tqsur(iunit)) + cr_advectheat(Qsub, THeat%Tqsub(iunit))
                THeat%Ha_t2r(iunit) = -cr_advectheat(abs(TRunoff%etout(iunit,nt_nliq)+TRunoff%etout(iunit,nt_nice)), THeat%Tt(iunit))
            ! change of energy due to heat exchange with the environment
            THeat%deltaH_t(iunit) = theDeltaT * (THeat%Hs_t(iunit) + THeat%Hl_t(iunit) + THeat%He_t(iunit) + THeat%Hc_t(iunit) + THeat%Hh_t(iunit))
            ! change of energy due to advective heat flux
            THeat%deltaM_t(iunit) = theDeltaT * (THeat%Ha_h2t(iunit)-cr_advectheat(Qsur + Qsub, THeat%Tt(iunit)))
    end subroutine subnetworkHeat

    subroutine subnetworkHeat_simple(iunit, theDeltaT)
    ! !DESCRIPTION: calculate the net heat balance of subnetwork channel.
    use shr_sys_mod , only : shr_sys_flush
        implicit none
        integer, intent(in) :: iunit
        real(r8), intent(in) :: theDeltaT    
        
        THeat%Hs_t(iunit) = 0._r8
        THeat%Hl_t(iunit) = 0._r8
        THeat%He_t(iunit) = 0._r8
        THeat%Hh_t(iunit) = 0._r8
        THeat%Hc_t(iunit) = 0._r8

        THeat%Ha_h2t(iunit) = 0._r8
        THeat%Ha_t2r(iunit) = -cr_advectheat(abs(TRunoff%etout(iunit,nt_nliq)+TRunoff%etout(iunit,nt_nice)), THeat%Tt(iunit))
        ! change of energy due to heat exchange with the environment
        THeat%deltaH_t(iunit) = 0._r8 
        ! change of energy due to advective heat flux
        THeat%deltaM_t(iunit) = 0._r8 

    end subroutine subnetworkHeat_simple

    
    subroutine mainchannelHeat(iunit, theDeltaT)
    ! !DESCRIPTION: calculate the net heat balance of main channel.
        use shr_sys_mod , only : shr_sys_flush
        implicit none
        integer, intent(in) :: iunit
        real(r8), intent(in) :: theDeltaT    
        integer    :: k
        real(r8) :: Ha_temp, Ha_temp1, Ha_temp2
        THeat%Ha_rin(iunit) = 0._r8
        !if(TUnit%fdir(iunit) >= 0 .and. TUnit%areaTotal(iunit) > TINYVALUE1 .and. TUnit%rlen(iunit) >= TINYVALUE1) then
                if(TRunoff%yr(iunit,nt_nliq) > TUnit%rdepth(iunit)) then ! TODO: the cases w/o flooding should be treated differently, here treated the same for now
                    !TRunoff%rarea(iunit,nt_nliq) = WaterAreaRatio*TUnit%rwidth0(iunit) * TUnit%rlen(iunit)
                    TRunoff%rarea(iunit,nt_nliq) = WaterAreaRatio*TUnit%rwidth(iunit) * TUnit%rlen(iunit)
                else
                    TRunoff%rarea(iunit,nt_nliq) = WaterAreaRatio*TUnit%rwidth(iunit) * TUnit%rlen(iunit)
                end if
                THeat%Hs_r(iunit) = cr_swrad(THeat%forc_solar(iunit), TRunoff%rarea(iunit,nt_nliq))
                THeat%Hl_r(iunit) = cr_lwrad(THeat%forc_lwrad(iunit), THeat%Tr(iunit), TRunoff%rarea(iunit,nt_nliq))
                THeat%He_r(iunit) = cr_latentheat(THeat%forc_t(iunit), THeat%forc_pbot(iunit), THeat%forc_vp(iunit), THeat%forc_wind(iunit), WindSheltering, THeat%Tr(iunit), TRunoff%rarea(iunit,nt_nliq))
                THeat%Hh_r(iunit) = cr_sensibleheat(THeat%forc_t(iunit), THeat%forc_pbot(iunit), THeat%forc_wind(iunit), WindSheltering, THeat%Tr(iunit), TRunoff%rarea(iunit,nt_nliq))
                THeat%Hc_r(iunit) = cr_condheat(THeat%Hs_r(iunit),THeat%Hl_r(iunit),TRunoff%rarea(iunit,nt_nliq))
                
                !do k=1,TUnit%nUp(iunit)
                !    THeat%Ha_rin(iunit) = THeat%Ha_rin(iunit) - THeat%Ha_rout(TUnit%iUp(iunit,k))
                !end do
                
                THeat%Ha_rin(iunit) = THeat%Ha_rin(iunit) - THeat%Ha_eroutUp(iunit) 
                !if(TUnit%indexDown(iunit) > 0) then
                    THeat%Ha_rout(iunit) = -cr_advectheat(abs(TRunoff%erout(iunit,nt_nliq)+TRunoff%erout(iunit,nt_nice)), THeat%Tr(iunit))
                !else
                !    THeat%Ha_rout(iunit) = 0._r8
                !end if
 
            ! change of energy due to heat exchange with the environment
            THeat%deltaH_r(iunit) = theDeltaT * (THeat%Hs_r(iunit) + THeat%Hl_r(iunit) + THeat%He_r(iunit) + THeat%Hc_r(iunit) + THeat%Hh_r(iunit))
            ! change of energy due to advective heat fluxes. Note here the advective heat flux is calculated differently from that by Huan Wu or van Vliet et al.
            ! Their routing model is based on source-to-sink, while our model is explicitly tracing inflow from each upstream channel.
            Ha_temp = cr_advectheat(TRunoff%erin(iunit,nt_nliq)+TRunoff%erin(iunit,nt_nice)+TRunoff%erlateral(iunit,nt_nliq)+TRunoff%erlateral(iunit,nt_nice),THeat%Tr(iunit))
            THeat%deltaM_r(iunit) = theDeltaT * (THeat%Ha_lateral(iunit) + THeat%Ha_rin(iunit) - Ha_temp)
            
        !end if
    end subroutine mainchannelHeat

    subroutine mainchannelHeat_simple(iunit, theDeltaT)
    ! !DESCRIPTION: calculate the net heat balance of main channel.
        use shr_sys_mod , only : shr_sys_flush
        implicit none
        integer, intent(in) :: iunit
        real(r8), intent(in) :: theDeltaT    
        integer    :: k
        real(r8) :: Ha_temp, Ha_temp1, Ha_temp2
        THeat%Ha_rin(iunit) = 0._r8
        !if(TUnit%fdir(iunit) >= 0 .and. TUnit%areaTotal(iunit) > TINYVALUE1) then
                THeat%Hs_r(iunit) = 0._r8
                THeat%Hl_r(iunit) = 0._r8
                THeat%He_r(iunit) = 0._r8
                THeat%Hh_r(iunit) = 0._r8
                THeat%Hc_r(iunit) = 0._r8
                THeat%Ha_rin(iunit) = 0._r8
                !if(TUnit%indexDown(iunit) > 0) then
                    THeat%Ha_rout(iunit) = -cr_advectheat(abs(TRunoff%erout(iunit,nt_nliq)+TRunoff%erout(iunit,nt_nice)), THeat%Tr(iunit))
                !else
                    !THeat%Ha_rout(iunit) = 0._r8
                !end if                

            ! change of energy due to heat exchange with the environment
            THeat%deltaH_r(iunit) = theDeltaT * (THeat%Hs_r(iunit) + THeat%Hl_r(iunit) + THeat%He_r(iunit) + THeat%Hc_r(iunit) + THeat%Hh_r(iunit))
            ! change of energy due to advective heat fluxes. Note here the advective heat flux is calculated differently from that by Huan Wu or van Vliet et al.
            ! Their routing model is based on source-to-sink, while our model is explicitly tracing inflow from each upstream channel.
            Ha_temp = cr_advectheat(TRunoff%erin(iunit,nt_nliq)+TRunoff%erin(iunit,nt_nice)+TRunoff%erlateral(iunit,nt_nliq)+TRunoff%erlateral(iunit,nt_nice),THeat%Tr(iunit))
            THeat%deltaM_r(iunit) = theDeltaT * (THeat%ha_lateral(iunit) + THeat%Ha_rin(iunit) - Ha_temp)
        !end if
    end subroutine mainchannelHeat_simple
    
    subroutine subnetworkTemp(iunit)
    ! !DESCRIPTION: calculate the water temperature of subnetwork channel.
        implicit none
        integer, intent(in) :: iunit
        
        real(r8) :: Mt  !mass of water (Kg)
        real(r8) :: Ttmp1, Ttmp2  !
        
        if((TRunoff%wt(iunit,nt_nliq)+TRunoff%wt(iunit,nt_nice)) > TINYVALUE1  .and. THeat%forc_t(iunit) > 200._r8) then
            Mt = TRunoff%wt(iunit,nt_nliq) * denh2o + TRunoff%wt(iunit,nt_nice) * denice
            THeat%Tt(iunit) = THeat%Tt(iunit) + (THeat%deltaH_t(iunit)+THeat%deltaM_t(iunit)) / (Mt * cpliq)
        else
            if(TRunoff%qsur(iunit,nt_nliq)+TRunoff%qsur(iunit,nt_nice) > TINYVALUE1) then
                THeat%Tt(iunit) = THeat%Tqsur(iunit) * (TRunoff%qsur(iunit,nt_nliq)+TRunoff%qsur(iunit,nt_nice)) + THeat%Tqsub(iunit) * (TRunoff%qsub(iunit,nt_nliq)+TRunoff%qsub(iunit,nt_nice))
                THeat%Tt(iunit) = THeat%Tt(iunit)/(TRunoff%qsur(iunit,nt_nliq)+TRunoff%qsur(iunit,nt_nice)+TRunoff%qsub(iunit,nt_nliq)+TRunoff%qsub(iunit,nt_nice))
            else
                THeat%Tt(iunit) = THeat%Tqsur(iunit)
            end if
        end if

        if(THeat%Tt(iunit) < 273.15_r8) then
            THeat%Tt(iunit) = 273.15_r8
        end if
        
    end subroutine subnetworkTemp

    subroutine subnetworkTemp_simple(iunit)
    ! !DESCRIPTION: calculate the water temperature of subnetwork channel.
        implicit none
        integer, intent(in) :: iunit
        
        real(r8) :: Mt  !mass of water (Kg)
        real(r8) :: Ttmp1, Ttmp2  !
        
        if(TRunoff%qsur(iunit,nt_nliq)+TRunoff%qsub(iunit,nt_nliq) > TINYVALUE1) then
            THeat%Tt(iunit) = THeat%Tqsur(iunit) * (TRunoff%qsur(iunit,nt_nliq)+TRunoff%qsur(iunit,nt_nice)) + THeat%Tqsub(iunit) * (TRunoff%qsub(iunit,nt_nliq)+TRunoff%qsub(iunit,nt_nice))
            THeat%Tt(iunit) = THeat%Tt(iunit)/(TRunoff%qsur(iunit,nt_nliq)+TRunoff%qsur(iunit,nt_nice)+TRunoff%qsub(iunit,nt_nliq)+TRunoff%qsub(iunit,nt_nice))
        else
            THeat%Tt(iunit) = THeat%Tqsur(iunit)
        end if
    end subroutine subnetworkTemp_simple
    
    
    subroutine mainchannelTemp(iunit)
    ! !DESCRIPTION: calculate the water temperature of subnetwork channel.
        use shr_sys_mod , only : shr_sys_flush
        implicit none
        integer, intent(in) :: iunit
        
        real(r8) :: Mr  !mass of water (Kg)
        real(r8) :: Ttmp1, Ttmp2  !
        
        if((TRunoff%wr(iunit,nt_nliq)+TRunoff%wr(iunit,nt_nice)) > TINYVALUE1 .and. THeat%forc_t(iunit) > 200._r8) then
            Mr = TRunoff%wr(iunit,nt_nliq) * denh2o + TRunoff%wr(iunit,nt_nice) * denice
            THeat%Tr(iunit) = THeat%Tr(iunit) + (THeat%deltaH_r(iunit)+THeat%deltaM_r(iunit)) / (Mr * cpliq)
        else
            THeat%Tr(iunit) = THeat%Tt(iunit)
        end if
        if(THeat%Tr(iunit) < 273.15_r8) then
            THeat%Tr(iunit) = 273.15_r8
        end if

    end subroutine mainchannelTemp

    subroutine mainchannelTemp_simple(iunit)
    ! !DESCRIPTION: calculate the water temperature of subnetwork channel.
    use shr_sys_mod , only : shr_sys_flush
        implicit none
        integer, intent(in) :: iunit
        
        real(r8) :: Mr  !mass of water (Kg)
        real(r8) :: Ttmp1, Ttmp2  !
        
        THeat%Tr(iunit) = THeat%Tt(iunit)
    end subroutine mainchannelTemp_simple

    subroutine reservoirHeat(iunit, theDeltaT)
    ! !DESCRIPTION: calculate the net heat balance of reservoir. invoked after the regulation subroutine
    ! simplified version as of 09/2014, now assuming the flow regulation won't modify the release water temperature directly
    ! to be extended later, e.g., incorporating the stratification processes
        implicit none
        integer, intent(in) :: iunit
        real(r8), intent(in) :: theDeltaT
        !if(TUnit%indexDown(iunit) > 0) then
            THeat%Ha_rout(iunit) = -cr_advectheat(abs(TRunoff%erout(iunit,nt_nliq)+TRunoff%erout(iunit,nt_nice)), THeat%Tr(iunit))
            !THeat%Ha_rout(iunit) = THeat%Ha_rout(iunit) - cr_advectheat(abs(TRunoff%erout(iunit,nt_nliq)+TRunoff%erout(iunit,nt_nice)), THeat%Tr(iunit))
        !else
        !    THeat%Ha_rout(iunit) = 0._r8
        !end if
    end subroutine reservoirHeat    

    subroutine reservoirTemp(iunit)
    ! !DESCRIPTION: calculate the water temperature of reservoir.
    ! simplified version as of 09/2014, to be extended later
        implicit none
        integer, intent(in) :: iunit
        
                
    end subroutine reservoirTemp        
    
    function cr_swrad(Hswin_, Aw_) result(Hsw_)
    ! closure relationship for net short-wave solar radiation
        implicit none
        real(r8), intent(in) :: Hswin_, Aw_  ! incoming short-wave radiation, surface area
        real(r8) :: Hsw_ ! [J/s]
        
        real(r8) :: Hswout_  ! short-wave radiation reflected by the water
        Hswout_ = 0.03_r8 * Hswin_  ! 3% of the incoming short-wave radiation [Wu et al., 2012; Herbert et al., 2011]
        Hsw_ = Hswin_ - Hswout_
        Hsw_ = Hsw_ * Aw_
        
        return
    end function cr_swrad 

    function cr_lwrad(Hlwin_, Tw_, Aw_) result(Hlw_)
    ! closure relationship for net long-wave radiation
        implicit none
        real(r8), intent(in) :: Hlwin_, Tw_, Aw_  ! incoming atmos. long-wave radiation, water temperature (kelvin), surface area (m2)
        real(r8) :: Hlw_ ! [J/s]
        
        real(r8) :: Hlwout_  ! long-wave radiation emitted by the water
        Hlwout_ = 0.97_r8 * sb * Tw_**4
        Hlw_ = Hlwin_ - Hlwout_
        Hlw_ = Hlw_ * Aw_
        
        return
    end function cr_lwrad 

    function cr_latentheat(Ta_, Pbot_, e_, U_, F_, Tw_, Aw_) result(He_)
    ! closure relationship for latent heat flux, [Wu et al., 2012]
        implicit none
        real(r8), intent(in) :: Ta_, Pbot_  ! air temperature (k), surface atmospheric pressure (pa)
        real(r8), intent(in) :: e_, U_, Tw_ ! atmos. vapor pressure (pa), wind speed (m/s), 
        real(r8), intent(in) :: F_, Aw_ ! dimensionless coefficient for the wind shltering by riparian vegetation, , surface area (m2)
        real(r8) :: He_ ! [J/s]
        
        real(r8) :: esat_  ! atmospheric saturated vapor pressure at certain temperature
        real(r8) :: esdT, qs, qsdT  ! d(es)/d(T), humidity (kg/kg), d(qs)/d(T) 
        real(r8) :: Kl_    ! empirical coefficient for the turbulent exchange of water vapor (mm/d/hpa)
        real(r8) :: Le_    ! latent heat of vaporization (J/Kg)
        real(r8) :: Evap_     ! evaporation rate (mm/d)
        
        call QSat(Tw_, Pbot_, esat_, esdT, qs, qsdT)        
        
        Kl_ = 0.211_r8 + 0.103_r8 * U_ * F_
        Le_ = (2.495_r8 - 2.36_r8 * 1.e-3 * (Tw_-273.15_r8)) * 1.e6    ! S. L. Dingman (2009), Fluivial Hydraulics
        Evap_  = Kl_ * (esat_ - e_)/100._r8  ! 100 here is for conversion from Pa to hPa
        He_ = -denh2o * Evap_ * Le_ / (86.4e6)
        He_ = He_ * Aw_
        
        return
        
    end function cr_latentheat
    
    function cr_sensibleheat(Ta_, Pbot_, U_, F_, Tw_, Aw_) result(Hh_)
    ! closure relationship for sensible heat flux, [Wu et al., 2012]
        implicit none
        real(r8), intent(in) :: Ta_, Pbot_  ! air temperature (k), surface atmospheric pressure (pa)
        real(r8), intent(in) :: U_, Tw_ ! wind speed (m/s), water temperature (K)
        real(r8), intent(in) :: F_ , Aw_ ! dimensionless coefficient for the wind shltering by riparian vegetation, surface area (m2)
        real(r8) :: Hh_ ! [J/s]
        
        real(r8) :: Kl_    ! empirical coefficient for the turbulent exchange of water vapor (mm/d/hpa)
        real(r8) :: Le_    ! latent heat of vaporization (J/Kg)
        real(r8) :: gamma = 0.655_r8  ! psychrometric constant at normal pressure, 0.655 hPa/Celcus
        real(r8) :: Pbot0 = 1013.25_r8 ! normal atmosphere pressure (hpa)
        
        Kl_ = 0.211_r8 + 0.103_r8 * U_ * F_
        Le_ = (2.495_r8 - 2.36_r8 * 1.e-3 * (Tw_-273.15_r8)) * 1.e6    ! S. L. Dingman (2009), Fluivial Hydraulics
        Hh_ = -gamma * (Pbot_/100._r8/Pbot0) * Kl_ * Le_ * (Tw_ - Ta_) * denh2o/(86.4e6)
        Hh_ = Hh_ * Aw_
        
        return
    end function cr_sensibleheat
    
    function cr_condheat(Hsw_, Hlw_, Aw_) result(Hc_)
    ! closure relationship for conductive heat flux
        implicit none
        real(r8), intent(in) :: Hsw_,Hlw_, Aw_  ! net short-wave radiation, surface area (m2)
        real(r8) :: Hc_ ! [J/s]
        
        Hc_ = 0.05_r8 * (Hsw_ + Hlw_)  !  [Wu et al., 2012]
        Hc_ = 0._r8    ! TODO: will be better represented along with the groundwater-river water interactions
        
        return
    end function cr_condheat 
    
    function cr_advectheat(Qin_, Twin_) result(Ha_)
    ! closure relationship for advective heat flux, assuming reference temperature (to define internal energy) is zero
        implicit none
        real(r8), intent(in) :: Qin_, Twin_ ! rate of water inflow (m3/s), temperature of water  inflow (K),
        real(r8) :: Ha_        ! conductive heat flux (J/s)
        
        Ha_ = denh2o * cpliq * Qin_ * Twin_
        
        return
    end function cr_advectheat 
    
    function cr_linear_Tw_Ta(iunit_,Ta_) result(Tw_)
    ! closure relationship to calculate water temperature based on the linear relationship proposed by Stefan and Preudhomme (1993)
        implicit none
        integer, intent(in) :: iunit_ ! 
        real(r8), intent(in) :: Ta_ ! temperature of air (Kelvin),
        real(r8) :: Tw_        ! temperature of water (Kelvin)

        Tw_ = 5.0_r8 + 0.75_r8 * (Ta_ - 273.15_r8) + 273.15_r8
    
        return
    end function cr_linear_Tw_Ta

    function cr_S_curve(iunit_,Ta_) result(Tw_)
    ! closure relationship to calculate water temperature based on the S-curve 
        implicit none
        integer, intent(in) :: iunit_ ! 
        real(r8), intent(in) :: Ta_ ! temperature of air (Kelvin),
        real(r8) :: Tw_        ! temperature of water (Kelvin)

        real(r8) :: Ttmp1, Ttmp2  !
        Ttmp1 = TPara%t_alpha(iunit_) - TPara%t_mu(iunit_)
        Ttmp2 = 1._r8 + exp(TPara%t_gamma(iunit_)*(TPara%t_beta(iunit_) - (Ta_-273.15_r8)))
        Tw_ = TPara%t_mu(iunit_) + Ttmp1/Ttmp2 + 273.15_r8
        if(Tw_ < 273.15_r8) then
            Tw_ = 273.15_r8
        end if
    
        return
    end function cr_S_curve

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: QSat
!
! !INTERFACE:
  subroutine QSat (T, p, es, esdT, qs, qsdT)
!
! !DESCRIPTION:
! Computes saturation mixing ratio and the change in saturation
! mixing ratio with respect to temperature.
! Reference:  Polynomial approximations from:
!             Piotr J. Flatau, et al.,1992:  Polynomial fits to saturation
!             vapor pressure.  Journal of Applied Meteorology, 31, 1507-1513.
!
! !USES:
    use shr_kind_mod , only: r8 => shr_kind_r8
    use shr_const_mod, only: SHR_CONST_TKFRZ
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: T        ! temperature (K)
    real(r8), intent(in)  :: p        ! surface atmospheric pressure (pa)
    real(r8), intent(out) :: es       ! vapor pressure (pa)
    real(r8), intent(out) :: esdT     ! d(es)/d(T)
    real(r8), intent(out) :: qs       ! humidity (kg/kg)
    real(r8), intent(out) :: qsdT     ! d(qs)/d(T)
!
! !CALLED FROM:
! subroutine Biogeophysics1 in module Biogeophysics1Mod
! subroutine BiogeophysicsLake in module BiogeophysicsLakeMod
! subroutine CanopyFluxesMod CanopyFluxesMod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
!
!
! !LOCAL VARIABLES:
!EOP
!
    real(r8) :: T_limit
    real(r8) :: td,vp,vp1,vp2
!
! For water vapor (temperature range 0C-100C)
!
    real(r8), parameter :: a0 =  6.11213476_r8
    real(r8), parameter :: a1 =  0.444007856_r8
    real(r8), parameter :: a2 =  0.143064234e-01_r8
    real(r8), parameter :: a3 =  0.264461437e-03_r8
    real(r8), parameter :: a4 =  0.305903558e-05_r8
    real(r8), parameter :: a5 =  0.196237241e-07_r8
    real(r8), parameter :: a6 =  0.892344772e-10_r8
    real(r8), parameter :: a7 = -0.373208410e-12_r8
    real(r8), parameter :: a8 =  0.209339997e-15_r8
!
! For derivative:water vapor
!
    real(r8), parameter :: b0 =  0.444017302_r8
    real(r8), parameter :: b1 =  0.286064092e-01_r8
    real(r8), parameter :: b2 =  0.794683137e-03_r8
    real(r8), parameter :: b3 =  0.121211669e-04_r8
    real(r8), parameter :: b4 =  0.103354611e-06_r8
    real(r8), parameter :: b5 =  0.404125005e-09_r8
    real(r8), parameter :: b6 = -0.788037859e-12_r8
    real(r8), parameter :: b7 = -0.114596802e-13_r8
    real(r8), parameter :: b8 =  0.381294516e-16_r8
!
! For ice (temperature range -75C-0C)
!
    real(r8), parameter :: c0 =  6.11123516_r8
    real(r8), parameter :: c1 =  0.503109514_r8
    real(r8), parameter :: c2 =  0.188369801e-01_r8
    real(r8), parameter :: c3 =  0.420547422e-03_r8
    real(r8), parameter :: c4 =  0.614396778e-05_r8
    real(r8), parameter :: c5 =  0.602780717e-07_r8
    real(r8), parameter :: c6 =  0.387940929e-09_r8
    real(r8), parameter :: c7 =  0.149436277e-11_r8
    real(r8), parameter :: c8 =  0.262655803e-14_r8
!
! For derivative:ice
!
    real(r8), parameter :: d0 =  0.503277922_r8
    real(r8), parameter :: d1 =  0.377289173e-01_r8
    real(r8), parameter :: d2 =  0.126801703e-02_r8
    real(r8), parameter :: d3 =  0.249468427e-04_r8
    real(r8), parameter :: d4 =  0.313703411e-06_r8
    real(r8), parameter :: d5 =  0.257180651e-08_r8
    real(r8), parameter :: d6 =  0.133268878e-10_r8
    real(r8), parameter :: d7 =  0.394116744e-13_r8
    real(r8), parameter :: d8 =  0.498070196e-16_r8
!-----------------------------------------------------------------------

    T_limit = T - SHR_CONST_TKFRZ
    if (T_limit > 100.0_r8) T_limit=100.0_r8
    if (T_limit < -75.0_r8) T_limit=-75.0_r8

    td       = T_limit
    if (td >= 0.0_r8) then
       es   = a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 &
            + td*(a5 + td*(a6 + td*(a7 + td*a8)))))))
       esdT = b0 + td*(b1 + td*(b2 + td*(b3 + td*(b4 &
            + td*(b5 + td*(b6 + td*(b7 + td*b8)))))))
    else
       es   = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 &
            + td*(c5 + td*(c6 + td*(c7 + td*c8)))))))
       esdT = d0 + td*(d1 + td*(d2 + td*(d3 + td*(d4 &
            + td*(d5 + td*(d6 + td*(d7 + td*d8)))))))
    endif

    es    = es    * 100._r8            ! pa
    esdT  = esdT  * 100._r8            ! pa/K

    vp    = 1.0_r8   / (p - 0.378_r8*es)
    vp1   = 0.622_r8 * vp
    vp2   = vp1   * vp

    qs    = es    * vp1             ! kg/kg
    qsdT  = esdT  * vp2 * p         ! 1 / K

  end subroutine QSat
    
    
  subroutine printTest1(nio)
      ! !DESCRIPTION: output the simulation results into external files
      implicit none
      integer, intent(in) :: nio        ! unit of the file to print
      
      !integer :: IDlist(1:5) = (/9958,7232,7241,6850,6852/)
      integer :: IDlist(1:4) = (/1233,1244,1138,1583/)
      integer :: ios,ii                    ! flag of io status
            

      write(unit=nio,fmt="(12(e20.11))") rtmCTL%runofflnd_nt1(IDlist(1)), rtmCTL%templand_Tchanr_nt1(IDlist(1)), rtmCTL%templand_Ttrib_nt1(IDlist(1)), &
                                         rtmCTL%runofflnd_nt1(IDlist(2)), rtmCTL%templand_Tchanr_nt1(IDlist(2)), rtmCTL%templand_Ttrib_nt1(IDlist(2)), &
                                         rtmCTL%runofflnd_nt1(IDlist(3)), rtmCTL%templand_Tchanr_nt1(IDlist(3)), rtmCTL%templand_Ttrib_nt1(IDlist(3)), &
                                         rtmCTL%runofflnd_nt1(IDlist(4)), rtmCTL%templand_Tchanr_nt1(IDlist(4)), rtmCTL%templand_Ttrib_nt1(IDlist(4))
      
      !write(unit=nio,fmt="(16(e20.11))") TRunoff%erout(IDlist(1),1), TRunoff%wr(IDlist(1),1), THeat%Tt(IDlist(1)), THeat%Tr(IDlist(1)), &
      !                                   TRunoff%erout(IDlist(2),1), TRunoff%wr(IDlist(2),1), THeat%Tt(IDlist(2)), THeat%Tr(IDlist(2)), &
      !                                     TRunoff%erout(IDlist(3),1), TRunoff%wr(IDlist(3),1), THeat%Tt(IDlist(3)), THeat%Tr(IDlist(3)), &
      !                                     TRunoff%erout(IDlist(4),1), TRunoff%wr(IDlist(4),1), THeat%Tt(IDlist(4)), THeat%Tr(IDlist(4))
      !write(unit=nio,fmt="(32(e20.11))") THeat%Hs_r(IDlist(1)), THeat%Hl_r(IDlist(1)), THeat%He_r(IDlist(1)), THeat%Hh_r(IDlist(1)), THeat%Hc_r(IDlist(1)), THeat%Ha_lateral(IDlist(1)), THeat%deltaH_r(IDlist(1)), THeat%deltaM_r(IDlist(1)), &
      !                                   THeat%Hs_r(IDlist(2)), THeat%Hl_r(IDlist(2)), THeat%He_r(IDlist(2)), THeat%Hh_r(IDlist(2)), THeat%Hc_r(IDlist(2)), THeat%Ha_lateral(IDlist(2)), THeat%deltaH_r(IDlist(2)), THeat%deltaM_r(IDlist(2)), &
      !                                     THeat%Hs_r(IDlist(3)), THeat%Hl_r(IDlist(3)), THeat%He_r(IDlist(3)), THeat%Hh_r(IDlist(3)), THeat%Hc_r(IDlist(3)), THeat%Ha_lateral(IDlist(3)), THeat%deltaH_r(IDlist(3)), THeat%deltaM_r(IDlist(3)), &
      !                                     THeat%Hs_r(IDlist(4)), THeat%Hl_r(IDlist(4)), THeat%He_r(IDlist(4)), THeat%Hh_r(IDlist(4)), THeat%Hc_r(IDlist(4)), THeat%Ha_lateral(IDlist(4)), THeat%deltaH_r(IDlist(4)), THeat%deltaM_r(IDlist(4))
      !write(unit=nio,fmt="((a10),(e20.11))") theTime, liqWater%flow(ii)
      !write(unit=nio,fmt="((a10),6(e20.11))") theTime, liqWater%qsur(ii), liqWater%qsub(ii), liqWater%etin(ii)/(TUnit%area(ii)*TUnit%frac(ii)), liqWater%erlateral(ii)/(TUnit%area(ii)*TUnit%frac(ii)), liqWater%erin(ii), liqWater%flow(ii)
      !if(liqWater%yr(ii) > 0._r8) then
      !    write(unit=nio,fmt="((a10),6(e20.11))") theTime, liqWater%mr(ii)/liqWater%yr(ii),liqWater%yr(ii), liqWater%vr(ii), liqWater%erin(ii), liqWater%erout(ii)/(TUnit%area(ii)*TUnit%frac(ii)), liqWater%flow(ii)
      !else
      !    write(unit=nio,fmt="((a10),6(e20.11))") theTime, liqWater%mr(ii)-liqWater%mr(ii),liqWater%yr(ii), liqWater%vr(ii), liqWater%erin(ii), liqWater%erout(ii)/(TUnit%area(ii)*TUnit%frac(ii)), liqWater%flow(ii)
      !end if
      !write(unit=nio,fmt="((a10),7(e20.11))") theTime, liqWater%erlateral(ii)/(TUnit%area(ii)*TUnit%frac(ii)), liqWater%wr(ii),liqWater%mr(ii), liqWater%yr(ii), liqWater%pr(ii), liqWater%rr(ii), liqWater%flow(ii)
      !write(unit=nio,fmt="((a10),7(e20.11))") theTime, liqWater%yh(ii), liqWater%dwh(ii),liqWater%etin(ii), liqWater%vr(ii), liqWater%erin(ii), liqWater%erout(ii)/(TUnit%area(ii)*TUnit%frac(ii)), liqWater%flow(ii)
  
  end subroutine printTest1
    
end MODULE MOSART_heat_mod