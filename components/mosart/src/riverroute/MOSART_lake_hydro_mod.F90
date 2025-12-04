!
MODULE MOSART_lake_hydro_mod
! Description: core code of lake hydrology module in MOSART framework
! A simple wier formula (Kindsvater and Carter (1959))is used to estimate outflow from lakes 
! Developed by Wondie Yigzaw, July 2020. 
! REVISION HISTORY:
! Revised by Hongyi Li, Jan. 2023
!-----------------------------------------------------------------------
    
! !USES:
    use shr_kind_mod  , only : r8 => shr_kind_r8
    use shr_const_mod , only : SHR_CONST_G
    use RunoffMod , only : Tctl, TUnit, TRunoff, THeat, TUnit_lake_r, TLake_r, TUnit_lake_t, TLake_t, TPara, rtmCTL
    use rof_cpl_indices, only : nt_nliq, nt_nice
    use RtmVar         , only : iulog, ngeom, nlayers, rstraflag, lakeflag
    use RtmTimeManager

    implicit none
    private

    real(r8), parameter :: TINYVALUE = 1.0e-14_r8  ! double precision variable has a significance of about 16 decimal digits
    
    public mosart_lake_hydro_sub_channel
    public mosart_lake_hydro_main_channel
    !public mosart_lake_routing
    
! !PUBLIC MEMBER FUNCTIONS:
    contains  
        
    subroutine mosart_lake_hydro_sub_channel(iunit,nt,localDeltaT)
    ! !DESCRIPTION: calculate the lake water balance.
    
        use shr_sys_mod , only : shr_sys_flush
    
        implicit none
        integer,  intent(in) :: iunit, nt
        real(r8), intent(in) :: localDeltaT
        character(len=*),parameter :: subname = '(mosart_lake_hydro)'
        
        integer :: j,i                   ! indices
        real(r8) :: delta_V              ! the excess volume (m3)
        real(r8) :: d_min                ! minimum lake depth to allow outflow (m)
        real(r8) :: delta_h              ! the hydraulic head difference (m)
        real(r8) :: dv_tmp               ! temperary storage change (m3)
        real(r8) :: inflow, outflow      ! Lake inflow, outflow (m3/s)
        
!**************************************************************************************************************************************************
        TLake_t%lake_outflow(iunit) = 0._r8
        if (TUnit_lake_t%lake_flg(iunit) >=1) then    ! Lake module active if there is natural lake
            !do i=2,nlayers
            !   if (TLake_t%v_zt(iunit,i)<TLake_t%v_zt(iunit,i-1) .or. TLake_t%d_v(iunit,i)<0.0_r8)return
            !end do
            !do j = 2, ngeom+1     
            !   if (TLake_t%a_di(iunit,j)<0._r8 .or. TLake_t%v_zti(iunit,j)<0._r8) return
            !end do
            
            !TLake_t%lake_inflow(iunit) = TRunoff%erlateral(iunit,nt) * TUnit_lake_t%F_local(iunit)
            TLake_t%lake_outflow(iunit) = 0._r8

            !! we don't have bankfull depth value for sub-network channel, so assume we can use the main channel bankfull depth
            !delta_h = min(TUnit%rdepth(iunit), TLake_t%d_lake(iunit) - TUnit_lake_t%h_min(iunit))
            !delta_h = min(2._r8-(TUnit_lake_t%h_lake(iunit)-TLake_t%d_lake(iunit)), TLake_t%d_lake(iunit) - TUnit_lake_t%h_min(iunit))
            if(TLake_t%d_lake(iunit) <= 0.1_r8) then ! no outflow when lake water depth is very shallow
                delta_h = 0._r8
            else
                d_min = max(TUnit_lake_t%h_lake(iunit) - 2._r8, 0._r8) !! TODO: here set the bankfull channel depth for t-zone as 2 meters
                delta_h = TLake_t%d_lake(iunit) - d_min  ! only allow outflow when lake water level exceeds the channel bed elevation
                if(delta_h < TINYVALUE) then
                    delta_h = 0._r8
                end if
                                        
                TLake_t%J_Min(iunit) = TLake_t%d_ns(iunit)
                if(delta_h >= TLake_t%d_lake(iunit)) then
                    delta_h = 0.95_r8 * TLake_t%d_lake(iunit)
                end if
            end if
            
            !if(iunit == 229628) then
            !    write(unit=1118,fmt="(i4,i4, i4, 4(e14.6))") 1, TLake_t%d_ns(iunit), TLake_t%J_Min(iunit), TLake_t%dd_z(iunit,TLake_t%d_ns(iunit)), TLake_t%d_v(iunit,TLake_t%d_ns(iunit)), TLake_t%lake_outflow(iunit) * localDeltaT, TLake_t%lake_outflow(iunit)
            !    write(unit=1119,fmt="(i4, i10, 4(e14.6))") 1, iunit, delta_h, d_min, TLake_t%dd_z(iunit,TLake_t%d_ns(iunit)), TLake_t%dd_z(iunit,TLake_t%d_ns(iunit)-1)
            !end if
            if (delta_h > TINYVALUE) then
                TLake_t%lake_outflow(iunit) = - CR_lake_Bernoulli(TUnit%twidth(iunit), delta_h)
                do i=TLake_t%d_ns(iunit), 1, -1
                    if(delta_h <= sum(TLake_t%dd_z(iunit,TLake_t%d_ns(iunit):i))) then
                        TLake_t%J_Min(iunit) = i
                        exit
                    end if
                end do
                
                !if(TLake_t%lake_outflow(iunit)*localDeltaT + delta_V < -TINYVALUE) then
                !    TLake_t%lake_outflow(iunit) = -delta_V / localDeltaT
                !end if
                
                ! make sure the lake storage is not depleted in a single time step
                !delta_V = TLake_t%V_str(iunit) - TUnit_lake_t%v_min(iunit)
                delta_V = sum(TLake_t%d_v(iunit,TLake_t%J_Min(iunit):TLake_t%d_ns(iunit)))
                if(delta_V > TINYVALUE) then
                    if(TLake_t%lake_outflow(iunit)*localDeltaT + 0.5_r8*delta_V < -TINYVALUE) then
                        TLake_t%lake_outflow(iunit) = -0.5*delta_V / localDeltaT
                    end if
                else
                    TLake_t%lake_outflow(iunit) = 0._r8
                end if
            !if(iunit == 229628) then
            !    write(unit=1118,fmt="(i4,i4, i4, 4(e14.6))") 2, TLake_t%d_ns(iunit), TLake_t%J_Min(iunit), TLake_t%dd_z(iunit,TLake_t%d_ns(iunit)), TLake_t%d_v(iunit,TLake_t%d_ns(iunit)), TLake_t%lake_outflow(iunit) * localDeltaT, TLake_t%lake_outflow(iunit)
            !    write(unit=1119,fmt="(i4,i10, 4(e14.6))") 2, iunit, delta_h, d_min, TLake_t%dd_z(iunit,TLake_t%d_ns(iunit)), TLake_t%dd_z(iunit,TLake_t%d_ns(iunit)-1)
            !end if
                
            !if(iunit == 229628) then
            !    write(unit=2118,fmt="(i10, i4, 3(e14.6))") TLake_t%d_ns(iunit), TLake_t%J_Min(iunit), TLake_t%d_v(iunit,TLake_t%d_ns(iunit)), TLake_t%lake_outflow(iunit) * localDeltaT, TLake_t%lake_outflow(iunit)
            !    write(unit=2119,fmt="(i10, 4(e14.6))") iunit, delta_h, d_min, TLake_t%dd_z(iunit,TLake_t%d_ns(iunit)), TLake_t%dd_z(iunit,TLake_t%d_ns(iunit)-1)
            !end if

            else
                TLake_t%lake_outflow(iunit) =  0._r8
            end if
            
            ! increase outflow if reservoir is too full already
            dv_tmp = sum(TLake_t%d_v(iunit,:)) - (TLake_t%v_zti(iunit,ngeom+1) - TLake_t%v_zti(iunit,1))*0.999_r8
            if((dv_tmp + TLake_t%lake_outflow(iunit)*localDeltaT) > TINYVALUE) then
                TLake_t%lake_outflow(iunit) = - dv_tmp / localDeltaT
                if(TLake_t%d_ns(iunit) > 1) then ! when flooded, allowing all layers to outflow
                    TLake_t%J_Min(iunit) = 1
                end if
            end if
            !if(iunit == 229628) then
            !    write(unit=1118,fmt="(i4,i4, i4, 4(e14.6))") 4, TLake_t%d_ns(iunit), TLake_t%J_Min(iunit), TLake_t%dd_z(iunit,TLake_t%d_ns(iunit)), TLake_t%d_v(iunit,TLake_t%d_ns(iunit)), TLake_t%lake_outflow(iunit) * localDeltaT, TLake_t%lake_outflow(iunit)
            !    write(unit=1119,fmt="(i4,i10, 4(e14.6))") 4, iunit, delta_h, d_min, TLake_t%dd_z(iunit,TLake_t%d_ns(iunit)), TLake_t%dd_z(iunit,TLake_t%d_ns(iunit)-1)
            !end if
            !TRunoff%erout(iunit,nt) = TLake_t%lake_outflow(iunit) ! now the lake outflow is discharged to downstream as the outflow from the current grid
            !if(TLake_t%lake_outflow(iunit)>TINYVALUE) then
            !   write(unit=18000,fmt="(i10, 5(e14.6))") iunit, TLake_t%lake_outflow(iunit), - CR_lake_Bernoulli(TUnit%twidth(iunit), delta_h), delta_h, delta_V, dv_tmp
            !end if        

        
        end if ! lake hydro
        
    end subroutine mosart_lake_hydro_sub_channel

    subroutine mosart_lake_hydro_main_channel(iunit,nt,localDeltaT)
    ! !DESCRIPTION: calculate the lake water balance.
    
        use shr_sys_mod , only : shr_sys_flush
    
        implicit none
        integer,  intent(in) :: iunit, nt
        real(r8), intent(in) :: localDeltaT
        character(len=*),parameter :: subname = '(mosart_lake_hydro)'
        
        integer :: j,i                   ! indices
        real(r8) :: delta_V              ! the excess volume (m3)
        real(r8) :: d_min                ! minimum lake depth to allow outflow (m)
        real(r8) :: delta_h              ! the hydraulic head difference (m)
        real(r8) :: dv_tmp               ! temperary storage change (m3)
        real(r8) :: inflow, outflow      ! Lake inflow, outflow (m3/s)
        
!**************************************************************************************************************************************************
        TLake_r%lake_outflow(iunit) = 0._r8
        if (TUnit_lake_r%lake_flg(iunit) >=1) then    ! Lake module active if there is natural lake
            
            !TLake_r%lake_inflow(iunit) = -TRunoff%erout(iunit,nt)
            !TLake_r%lake_outflow(iunit) = 0._r8
            
            if(TLake_r%d_lake(iunit) <= 0.1_r8) then ! no outflow when lake water depth is very shallow
                delta_h = 0._r8
            else
                d_min = max(TUnit_lake_r%h_lake(iunit) - TUnit%rdepth(iunit), 0._r8)
                delta_h = TLake_r%d_lake(iunit) - d_min  ! only allow outflow when lake water level exceeds the channel bed elevation
                if(delta_h < TINYVALUE) then
                    delta_h = 0._r8
                end if
            
                TLake_r%J_Min(iunit) = TLake_r%d_ns(iunit)
                if(delta_h >= TLake_r%d_lake(iunit)) then
                    delta_h = 0.95_r8 * TLake_r%d_lake(iunit)
                end if
            end if
            
            if (delta_h > TINYVALUE) then
                outflow = - CR_lake_Bernoulli(TUnit%twidth(iunit), delta_h)
                do i=TLake_r%d_ns(iunit), 1, -1
                    if(delta_h <= sum(TLake_r%dd_z(iunit,TLake_r%d_ns(iunit):i))) then
                        TLake_r%J_Min(iunit) = i
                        exit
                    end if
                end do
                
                !if(TLake_r%lake_outflow(iunit)*localDeltaT + delta_V < -TINYVALUE) then
                !    TLake_r%lake_outflow(iunit) = -delta_V / localDeltaT
                !end if
                
                ! make sure the lake storage is not depleted in a single time step
                !delta_V = TLake_r%V_str(iunit) - TUnit_lake_r%v_min(iunit)
                delta_V = sum(TLake_r%d_v(iunit,TLake_r%J_Min(iunit):TLake_r%d_ns(iunit)))
                if(delta_V > TINYVALUE) then
                    if(outflow*localDeltaT + 0.5_r8*delta_V < -TINYVALUE) then
                        outflow = -0.5*delta_V / localDeltaT
                    end if
                else
                    outflow = 0._r8
                end if
                
            else
                outflow = 0._r8
            end if
            
            ! increase outflow if reservoir is too full already
            dv_tmp = sum(TLake_r%d_v(iunit,:)) - (TLake_r%v_zti(iunit,ngeom+1) - TLake_r%v_zti(iunit,1))*0.999_r8
            if((dv_tmp + outflow*localDeltaT) > TINYVALUE) then
                outflow = - dv_tmp / localDeltaT
                if(TLake_r%d_ns(iunit) > 1) then ! when flooded, allowing all layers to outflow
                    TLake_r%J_Min(iunit) = 1
                end if
            end if
            
			TLake_r%lake_outflow(iunit) = TLake_r%lake_outflow(iunit) + outflow
            !if(iunit == 229628) then
            !    write(unit=1118,fmt="(i4, i10, 4(e14.6))") 1, TLake_r%d_ns(iunit), TUnit_lake_r%V_max(iunit)*1e6, TLake_r%d_v(iunit,TLake_r%d_ns(iunit)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), TLake_r%lake_outflow(iunit) * localDeltaT
            !    write(unit=1119,fmt="(i10, 4(e14.6))") iunit, delta_h, d_min, TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)-1)
            !end if
            
            !if(iunit == 229628) then
            !    write(unit=1118,fmt="(i4, i10, 4(e14.6))") 2, TLake_r%d_ns(iunit), TUnit_lake_r%V_max(iunit)*1e6, TLake_r%d_v(iunit,TLake_r%d_ns(iunit)), TLake_r%v_zt(iunit, TLake_r%d_ns(iunit) + 1), TLake_r%lake_outflow(iunit) * localDeltaT
            !    write(unit=1119,fmt="(i10, 4(e14.6))") iunit, delta_h, d_min, TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)), TLake_r%dd_z(iunit,TLake_r%d_ns(iunit)-1)
            !end if
            !TRunoff%erout(iunit,nt) = TLake_r%lake_outflow(iunit) ! now the lake outflow is discharged to downstream as the outflow from the current grid
           
        end if ! lake hydro
        
    end subroutine mosart_lake_hydro_main_channel


    function mosart_lake_routing(iunit,nt,localDeltaT) result(outflow_)
    ! !DESCRIPTION: calculate the lake outflow.
    
        use shr_sys_mod , only : shr_sys_flush
    
        implicit none
        integer,  intent(in) :: iunit, nt
        real(r8), intent(in) :: localDeltaT
        real(r8) :: outflow_    ! Lake outflow (m3/s)
        character(len=*),parameter :: subname = '(mosart_lake_routing)'
        
        integer  :: j,i                   ! indices
        real(r8) :: delta_V              ! the excess volume (m3)
        real(r8) :: delta_h              ! the hydraulic head difference (m)
        
!**************************************************************************************************************************************************
        outflow_ = 0._r8
            
        delta_V = TLake_r%V_str(iunit)-0.75*TLake_r%v_zti(iunit,ngeom+1)
        
        if (delta_V > TINYVALUE) then
            delta_h = delta_V/(TUnit_lake_r%A_max(iunit)*1e6) ! here unit for V_max is mcm, and for A_max is km2
            outflow_ = CR_lake_Bernoulli(TUnit%rwidth(iunit), delta_h)
            ! a treatment for water balance to avoid negative lake storage when evap. and outflow both large
            if(outflow_*localDeltaT - delta_V > TINYVALUE) then
                outflow_ = delta_V / localDeltaT
            end if      
        else
            outflow_ = 0._r8
        end if
                        
        return
        
    end function mosart_lake_routing

  
!-----------------------------------------------------------------------

    function CR_lake_Bernoulli(L_, Delta_h_) result(elake_out_)
    ! Function for calculating the outflow from a natural lake using the Bernoulli's equation
      implicit none
      real(r8), intent(in) :: L_! weir crest width [m]
      real(r8), intent(in) :: Delta_h_ ! difference between lake surface and weir crest [m]
      real(r8)             :: elake_out_            ! outflow from the lake, [m3/s]
      
      character(len=*),parameter :: subname = '(CR_lake_Bernoulli)'
      
      real(r8):: C_! weir coefficient [-]
      real(r8):: g_! gravity [m/s^2]
      
      C_ = 0.6_r8 ! Here we take a generic value for all sharp and broad weirs
      g_ = SHR_CONST_G
      
      if(Delta_h_ <= TINYVALUE) then
          elake_out_ = 0._r8
      else 
          elake_out_ = 2._r8/3._r8*C_*sqrt(2_r8*g_)*L_*(Delta_h_**(3._r8/2._r8))
      end if

        return
                  
    end function CR_lake_Bernoulli
  
    function cr_pet(Ta_, Pbot_, e_, U_, F_, Tw_) result(evap_)
    ! closure relationship for water surface evaporation, [Wu et al., 2012]
        implicit none
        real(r8), intent(in) :: Ta_, Pbot_  ! air temperature (k), surface atmospheric pressure (pa)
        real(r8), intent(in) :: e_, U_, Tw_ ! atmos. vapor pressure (pa), wind speed (m/s), 
        real(r8), intent(in) :: F_ ! dimensionless coefficient for the wind shltering by riparian vegetation, , surface area (m2)
        real(r8) :: Evap_     ! evaporation rate (mm/d)
        
        real(r8) :: esat_  ! atmospheric saturated vapor pressure at certain temperature
        real(r8) :: esdT, qs, qsdT  ! d(es)/d(T), humidity (kg/kg), d(qs)/d(T) 
        real(r8) :: Kl_    ! empirical coefficient for the turbulent exchange of water vapor (mm/d/hpa)
        real(r8) :: Le_    ! latent heat of vaporization (J/Kg)
        !real(r8) :: Evap_     ! evaporation rate (mm/d)
        
        call QSat0(Tw_, Pbot_, esat_, esdT, qs, qsdT)        
        
        Kl_ = 0.211_r8 + 0.103_r8 * U_ * F_
        Le_ = (2.495_r8 - 2.36_r8 * 1.e-3 * (Tw_-273.15_r8)) * 1.e6    ! S. L. Dingman (2009), Fluivial Hydraulics
        Evap_  = Kl_ * (esat_ - e_)/100._r8  ! 100 here is for conversion from Pa to hPa
        
        Evap_ = Evap_ /86.4e6_r8  ! mm/d --> m/s

        return
        
    end function cr_pet

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: QSat
!
! !INTERFACE:
  subroutine QSat0(T, p, es, esdT, qs, qsdT)
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

  end subroutine QSat0
     
end MODULE MOSART_lake_hydro_mod