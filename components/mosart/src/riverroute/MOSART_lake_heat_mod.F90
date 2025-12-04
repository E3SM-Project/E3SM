!
MODULE mosart_lake_heat_mod
! Description: MOSART-lake heat subroutines
! Developed by Hongyi Li, Jan. 2023. 
! REVISION HISTORY:
! 
!-----------------------------------------------------------------------
    
! !USES:
    use shr_kind_mod  , only : r8 => shr_kind_r8
    use shr_const_mod   , only : SHR_CONST_PI
    use shr_sys_mod   , only : shr_sys_abort
    use RtmTimeManager
    use RunoffMod , only : Tctl, TUnit, TRunoff, THeat, TUnit_lake_r, TLake_r, TPara, rtmCTL, TUnit_lake_t, TLake_t
    use rof_cpl_indices, only : nt_nliq, nt_nice
    use RtmVar         , only : iulog

    implicit none
    
    public temp_shallow_lake
    public den
    public QSat_str
    public QSat_huang
    public doy
    
! !PUBLIC MEMBER FUNCTIONS:
    contains  
    
    function den(t_z) result(rho)
    ! calculate density from temperature 
        implicit none
        real(r8), intent(in) :: t_z! Temperature (k) 
        real(r8) :: rho! ! density (kg/m3)
        
        ! varying density with temperature is causing some numerical instability
        rho = 1000._r8*( 1._r8 - 1.9549e-05*(abs(t_z-277._r8))**1.68_r8) ! modified from Subin et al, 2011 with lake ice fraction = 0
        
    end function den 
 
    function doy() result(doy_)
    ! calculate the day of a year 
        implicit none
        integer  :: doy_

        integer  :: yr, mon, day, ymd, tod      ! time information
		integer  :: days_mon(12)  ! number of days in each month
        integer  :: im
        
		days_mon = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
		call get_curr_date(yr, mon, day, tod)
		
		doy_ = day
		do im=1,mon-1
		    doy_ = doy_ + days_mon(im)
		end do
        
    end function doy 
 
    
    function QSat_huang (Ta) result(es)
!
! !DESCRIPTION:
! Computes saturation mixing ratio and the change in saturation
! mixing ratio with respect to temperature.
! Reference:  Huang, J. (2018). A simple accurate formula for calculating saturation vapor pressure of water and ice. 
!             Journal of Applied Meteorology and Climatology, 57(6), 1265-1272.

! !USES:
        use shr_kind_mod , only: r8 => shr_kind_r8
        use shr_const_mod, only: SHR_CONST_TKFRZ
    !
    ! !ARGUMENTS:
        implicit none
        real(r8), intent(in)  :: Ta       ! temperature (Celcus)
        real(r8) :: es       ! vapor pressure (pa)
        real(r8) :: tmp1, tmp2 ! 
		
!! !LOCAL VARIABLES:
!EOP
! 
        if(Ta > 0._r8) then ! water
		    tmp1 = exp(34.494 - 4924.99/(Ta+237.1))
			tmp2 = (Ta + 105)**1.57
		    es = tmp1 / tmp2
		else ! ice
		    tmp1 = exp(43.494 - 6545.8/(Ta+278))
			tmp2 = (Ta + 868)**2
			es = tmp1/tmp2
		end if

if(0>1) then  ! the improved Magnus formula
        if(Ta > 0._r8) then ! water
		    es = 610.94_r8 * exp((17.625_r8* Ta)/(Ta + 243.04))
		else ! ice
		    es = 611.21_r8 * exp((22.587_r8* Ta)/(Ta + 273.86))
		end if
end if
    
    end function QSat_huang

    subroutine QSat_str (T, p, es)
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
!! !LOCAL VARIABLES:
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
    !-----------------------------------------------------------------------

        T_limit = T - SHR_CONST_TKFRZ
        if (T_limit > 100.0_r8) T_limit=100.0_r8
        if (T_limit < -75.0_r8) T_limit=-75.0_r8

        td       = T_limit
        if (td >= 0.0_r8) then
           es   = a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 &
                + td*(a5 + td*(a6 + td*(a7 + td*a8)))))))
        else
           es   = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 &
                + td*(c5 + td*(c6 + td*(c7 + td*c8)))))))
        endif

        es    = es    * 100._r8            ! pa
    
    end subroutine QSat_str
  
    function temp_shallow_lake(myarray) result(Tw)
!*******************************************************************************************************
!     calculate the lake water column (including surface) temperature following Zhao and Gao, 2019
!*******************************************************************************************************    
        use RtmVar          , only : iulog
        
        implicit none
        character(len=*),parameter :: subname = '(temp_shallow_lake)'
        real(r8),dimension(14), intent(in) :: myarray!  input array
        real(r8) :: Tw     !  waterbody temperature (celcus)            

        real(r8) :: es     !  saturated vapor pressure (kPa)            
        real(r8) :: ea     !  vapor pressure (kPa)            
        real(r8) :: Ta     !  air temperature (celcus)            
        real(r8) :: u2     !  wind speed at the height of 2m (m/s)        
        real(r8) :: pbot   !  atmospheric pressure (kPa)        
        real(r8) :: Lf     !  fetch length (m)        
        real(r8) :: h_     !  lake depth (m)        
        real(r8) :: phi    !  location latitude (radius)      
        real(r8) :: solar  !  downward shortwave radiation radiation (MJ·m-2·d-1)
        real(r8) :: elev   !  lake surface elevation (m)
        real(r8) :: alb_s  !  surface albedo shortwave radiation (-)
        real(r8) :: Tw0    !  temperature of water body from previous time step (celcus)
        real(r8) :: doy    !  day of the year (day)
        real(r8) :: deltaT_!  time step (day)

        integer  :: yr, mon, day, ymd, tod      ! time information
        real(r8) :: pi = 3.1415916  ! Water density  (kg/m3)
        real(r8) :: Gsc = 4.92      ! solar constant (4.92 MJ·m-2·h-1)
        real(r8) :: rho_w = 1.e3    ! Water density  (kg/m3)
        real(r8) :: c_w = 0.0042    ! specific heat of water (MJ·kg−1·°C−1);
        real(r8) :: sigma = 4.9e-9  ! Stefan-Boltzman constant ( MJ·m−2·K−4·d−1)
        real(r8) :: k = 0.46        ! constant (MJ·m−2·d−1·°C−1)
        real(r8) :: b = 28.38       ! constant (28.38 MJ·m−2·d−1)
        real(r8) :: emi_w = 0.97    ! emisivty of water (-)
        real(r8) :: emi_a           ! emisivty of air with cloudiness factor (-)
        real(r8) :: Te              ! dew-point temperature (celcus)
        real(r8) :: Td              ! dew-point temperature (celcus)
        real(r8) :: Twb             ! wet-bulb temperature (celcus)
        real(r8) :: Swb             ! the slope of the saturation vapor pressure curve at Twb (kpa/celcus)
        real(r8) :: f_u2            ! the wind function (MJ·m−2·d−1·kPa−1)
        real(r8) :: lambda_v        ! latent heat of vaporization (MJ·kg−1)
        real(r8) :: gamma           ! psychrometric constant(kPa·°C-1)
        real(r8) :: tau             ! lab time(day)
        real(r8) :: J_              ! the of year(-)
        real(r8) :: delta           ! solar declination (radians)
        real(r8) :: omega_s         ! sunset hour angle (radians)
        real(r8) :: dr              ! inverse relative earth-sun distance factor (-)
        real(r8) :: Ket             ! extraterrestrial radiation (MJ·m-2·d-1)
        real(r8) :: Kso             ! clear-sky shortwave radiation (MJ·m-2·d-1)
        real(r8) :: Kr              ! ratio of incoming shortwave radiation to clear-sky shortwave radiation (-)
        real(r8) :: fcd             ! cloudiness factor (-)
        real(r8) :: Ses             ! the slope of the saturation vapor pressure curve at air temperature (kpa/celcus)
        
        real(r8) :: temp1, temp2 ! local variables
        
		es      = myarray(1)
		ea      = myarray(2)
		Ta      = myarray(3)
		u2      = myarray(4)
		pbot    = myarray(5)
	    Lf      = myarray(6)
	    h_      = myarray(7)
		phi     = myarray(8)
		solar   = myarray(9)
		elev    = myarray(10)
		alb_s   = myarray(11)
		Tw0     = myarray(12)
        doy     = myarray(13)
        deltaT_ = myarray(14)		
	 
        Td = (116.9_r8 + 237.3_r8 * log(ea)) / (16.78 - log(ea))        
        temp1 = 0.00066_r8 * 100_r8
        temp2 = 4098._r8 * ea / ((Td + 237.3_r8)*(Td + 237.3_r8))
        Twb = (temp1 * Ta + temp2 * Td) / (temp1 + temp2)
        
        temp1 = (Twb + 237.3_r8)
        temp2 = exp(17.27_r8 * Twb / temp1)
        Swb = 4098_r8 * 0.6108_r8 * temp2 / (temp1 * temp1)
        
        lambda_v = 2.501_r8 - 2.361 * 1.e-3 * Ta
        f_u2 = lambda_v * (2.33 + 1.65 * u2) * (Lf**(-0.1))
        !p = p * 0.001_r8  ! pa to kpa
        gamma = 0.00163_r8 * pbot / lambda_v
        tau = rho_w * c_w * h_ /(4.0*sigma *((Twb+273.15)**3) + f_u2 * (Swb + gamma))
        
        !! now calculate equilibrium temperature
        !Couliness factor
        J_ = doy !76.2 !doy
        delta = 0.409 * sin(2*pi*J_/365 - 1.39)
        temp1 = -tan(phi) * tan(delta)
        temp2 = sqrt(1 - ((tan(phi))**2)*((tan(delta))**2))
        omega_s = pi/2 - atan(temp1/temp2)
        dr = 1 + 0.033*cos(2*pi*J_/365)
        Ket = (24/pi)*Gsc*dr*(omega_s*sin(phi)*sin(delta) + cos(phi)*cos(delta)*sin(omega_s))
		Kso = (0.75 + 2*1e-5*elev)*Ket
        !solar = solar * 1e-6 * 86400 ! from watt·m-2 (J·m-2·s-1) to MJ·m-2·d-1
		Kr = min(solar / Kso, 1._r8)
		fcd = 1 - Kr
		! emissivity of air
		temp1 = 1 + 0.22 * (fcd**2.75)
		temp2 = 1.08 *(1-exp(-(10*ea)**((273.15+Ta)/2016)))
		emi_a = temp1 * temp2
		! now calculate Te
        temp1 = (Ta + 237.3_r8)
        temp2 = exp(17.27_r8 * Ta / temp1)
        Ses = 4098_r8 * 0.6108_r8 * temp2 / (temp1 * temp1)
		
		temp1 = (k*emi_a + f_u2*(Ses+gamma)) * Ta + (1-alb_s) * solar - b*(emi_w - emi_a) - f_u2 * (es - ea)
		temp2 = k*emi_w + f_u2 * (Ses + gamma)
		Te = temp1/temp2
		
		! now calculate updated lake water temperature
		!deltaT = deltaT / 86400  ! seconds to day
		Tw = Te + (Tw0 - Te) * exp(-deltaT_/tau)
        
    end function temp_shallow_lake

    function temp_shallow_lake_bkp(myarray) result(Tw)
!*******************************************************************************************************
!     calculate the lake water column (including surface) temperature following Zhao and Gao, 2019
!*******************************************************************************************************    
        use RtmVar          , only : iulog
        
        implicit none
        character(len=*),parameter :: subname = '(temp_shallow_lake)'
        real(r8),dimension(14), intent(in) :: myarray!  input array
        real(r8) :: Tw     !  waterbody temperature (celcus)            

        real(r8) :: es     !  saturated vapor pressure (kPa)            
        real(r8) :: ea     !  vapor pressure (kPa)            
        real(r8) :: Ta     !  air temperature (celcus)            
        real(r8) :: u2     !  wind speed at the height of 2m (m/s)        
        real(r8) :: P      !  atmospheric pressure (kPa)        
        real(r8) :: Lf     !  fetch length (m)        
        real(r8) :: h_     !  lake depth (m)        
        real(r8) :: phi    !  location latitude (radius)      
        real(r8) :: solar  !  downward shortwave radiation radiation (MJ·m-2·d-1)
        real(r8) :: elev   !  lake surface elevation (m)
        real(r8) :: alb_s  !  surface albedo shortwave radiation (-)
        real(r8) :: Tw0    !  temperature of water body from previous time step (celcus)
        real(r8) :: doy    !  day of the year (day)
        real(r8) :: deltaT_!  time step (day)

        integer  :: yr, mon, day, ymd, tod      ! time information
        real(r8) :: pi = 3.1415916  ! Water density  (kg/m3)
        real(r8) :: Gsc = 4.92      ! solar constant (4.92 MJ·m-2·h-1)
        real(r8) :: rho_w = 1.e3    ! Water density  (kg/m3)
        real(r8) :: c_w = 0.0042    ! specific heat of water (MJ·kg−1·°C−1);
        real(r8) :: sigma = 4.9e-9  ! Stefan-Boltzman constant ( MJ·m−2·K−4·d−1)
        real(r8) :: k = 0.46        ! constant (MJ·m−2·d−1·°C−1)
        real(r8) :: b = 28.38       ! constant (28.38 MJ·m−2·d−1)
        real(r8) :: emi_w = 0.97    ! emisivty of water (-)
        real(r8) :: emi_a           ! emisivty of air with cloudiness factor (-)
        real(r8) :: Te              ! dew-point temperature (celcus)
        real(r8) :: Td              ! dew-point temperature (celcus)
        real(r8) :: Twb             ! wet-bulb temperature (celcus)
        real(r8) :: Swb             ! the slope of the saturation vapor pressure curve at Twb (kpa/celcus)
        real(r8) :: f_u2            ! the wind function (MJ·m−2·d−1·kPa−1)
        real(r8) :: lambda_v        ! latent heat of vaporization (MJ·kg−1)
        real(r8) :: gamma           ! psychrometric constant(kPa·°C-1)
        real(r8) :: tau             ! lab time(day)
        real(r8) :: J_              ! the of year(-)
        real(r8) :: delta           ! solar declination (radians)
        real(r8) :: omega_s         ! sunset hour angle (radians)
        real(r8) :: dr              ! inverse relative earth-sun distance factor (-)
        real(r8) :: Ket             ! extraterrestrial radiation (MJ·m-2·d-1)
        real(r8) :: Kso             ! clear-sky shortwave radiation (MJ·m-2·d-1)
        real(r8) :: Kr              ! ratio of incoming shortwave radiation to clear-sky shortwave radiation (-)
        real(r8) :: fcd             ! cloudiness factor (-)
        real(r8) :: Ses             ! the slope of the saturation vapor pressure curve at air temperature (kpa/celcus)
        
        real(r8) :: temp1, temp2 ! local variables
        
		es      = myarray(1)
		ea      = myarray(2)
		Ta      = myarray(3)
		u2      = myarray(4)
		P       = myarray(5)
	    Lf      = myarray(6)
	    h_      = myarray(7)
		phi     = myarray(8)
		solar   = myarray(9)
		elev    = myarray(10)
		alb_s   = myarray(11)
		Tw0     = myarray(12)
        doy     = myarray(13)
        deltaT_ = myarray(14)		
	 
        Td = (116.9_r8 + 237.3_r8 * log(ea)) / (16.78 - log(ea))
		write(iulog,*) 'lake evap function result - h_ ', h_
		write(iulog,*) 'lake evap function result - Ta ', Ta
		write(iulog,*) 'lake evap function result - es ', es
		write(iulog,*) 'lake evap function result - ea ', ea
		write(iulog,*) 'lake evap function result - Td ', Td
        !Ta = Ta - 273.15_r8 ! kelvin to celcus
        
        temp1 = 0.00066_r8 * 100_r8
        temp2 = 4098._r8 * ea / ((Td + 237.3_r8)*(Td + 237.3_r8))
        Twb = (temp1 * Ta + temp2 * Td) / (temp1 + temp2)
		write(iulog,*) 'lake evap function result - Twb ', Twb
        
        temp1 = (Twb + 237.3_r8)
        temp2 = exp(17.27_r8 * Twb / temp1)
        Swb = 4098_r8 * 0.6108_r8 * temp2 / (temp1 * temp1)
		write(iulog,*) 'lake evap function result - Swb ', Swb
        
        lambda_v = 2.501_r8 - 2.361 * 1.e-3 * Ta
        f_u2 = lambda_v * (2.33 + 1.65 * u2) * (Lf**(-0.1))
		write(iulog,*) 'lake evap function result - f_u2 ', f_u2
        !p = p * 0.001_r8  ! pa to kpa
        gamma = 0.00163_r8 * P / lambda_v
        tau = rho_w * c_w * h_ /(4.0*sigma *((Twb+273.15)**3) + f_u2 * (Swb + gamma))
		write(iulog,*) 'lake evap function result - tau ', tau
        
        !! now calculate equilibrium temperature
        !Couliness factor
        J_ = doy !76.2 !doy
        delta = 0.409 * sin(2*pi*J_/365 - 1.39)
		write(iulog,*) 'lake evap function result - J_, delta ', J_, delta
        temp1 = -tan(phi) * tan(delta)
        temp2 = sqrt(1 - ((tan(phi))**2)*((tan(delta))**2))
        omega_s = pi/2 - atan(temp1/temp2)
		write(iulog,*) 'lake evap function result - omega_s ', omega_s
        dr = 1 + 0.033*cos(2*pi*J_/365)
        Ket = (24/pi)*Gsc*dr*(omega_s*sin(phi)*sin(delta) + cos(phi)*cos(delta)*sin(omega_s))
		Kso = (0.75 + 2*1e-5*elev)*Ket
		write(iulog,*) 'lake evap function result - Ket ', Ket
		write(iulog,*) 'lake evap function result - Kso ', Kso
        !solar = solar * 1e-6 * 86400 ! from watt·m-2 (J·m-2·s-1) to MJ·m-2·d-1
		Kr = min(solar / Kso, 1._r8)
		write(iulog,*) 'lake evap function result - solar ', solar
		write(iulog,*) 'lake evap function result - Kr ', Kr
		fcd = 1 - Kr
		write(iulog,*) 'lake evap function result - fcd ', fcd
		! emissivity of air
		temp1 = 1 + 0.22 * (fcd**2.75)
		temp2 = 1.08 *(1-exp(-(10*ea)**((273.15+Ta)/2016)))
		emi_a = temp1 * temp2
		write(iulog,*) 'lake evap function result - emi_a ', emi_a
		! now calculate Te
        temp1 = (Ta + 237.3_r8)
        temp2 = exp(17.27_r8 * Ta / temp1)
        Ses = 4098_r8 * 0.6108_r8 * temp2 / (temp1 * temp1)
		write(iulog,*) 'lake evap function result - Ses ', Ses
		
		temp1 = (k*emi_a + f_u2*(Ses+gamma)) * Ta + (1-alb_s) * solar - b*(emi_w - emi_a) - f_u2 * (es - ea)
		temp2 = k*emi_w + f_u2 * (Ses + gamma)
		Te = temp1/temp2
		write(iulog,*) 'lake evap function result - Te ', Te
		
		! now calculate updated lake water temperature
		!deltaT = deltaT / 86400  ! seconds to day
		write(iulog,*) 'lake evap function result - Tw0 ', Tw0
		write(iulog,*) 'lake evap function result - deltaT ', deltaT_
		Tw = Te + (Tw0 - Te) * exp(-deltaT_/tau)
        
    end function temp_shallow_lake_bkp
  
end MODULE mosart_lake_heat_mod