module atm2lndType

  use shr_kind_mod  , only : r8 => shr_kind_r8
  use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
  use clm_varpar    , only : numrad
  use decompMod     , only : bounds_type
  use clm_varcon    , only : rair, grav, cpair, hfus, tfrz, spval
  use clm_varctl    , only : iulog, use_c13, use_cn, use_lch4, use_cndv, use_fates

  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save

  type, public :: atm2lnd_type
     ! atm->lnd not downscaled
     real(r8), pointer :: forc_u_grc                    (:)   => null() ! atm wind speed, east direction (m/s)
     real(r8), pointer :: forc_v_grc                    (:)   => null() ! atm wind speed, north direction (m/s)
     real(r8), pointer :: forc_wind_grc                 (:)   => null() ! atmospheric wind speed
     real(r8), pointer :: forc_hgt_grc                  (:)   => null() ! atmospheric reference height (m)
     real(r8), pointer :: forc_hgt_u_grc                (:)   => null() ! obs height of wind [m] (new)
     real(r8), pointer :: forc_hgt_t_grc                (:)   => null() ! obs height of temperature [m] (new)
     real(r8), pointer :: forc_hgt_q_grc                (:)   => null() ! obs height of humidity [m] (new)
     real(r8), pointer :: forc_vp_grc                   (:)   => null() ! atmospheric vapor pressure (Pa)
     real(r8), pointer :: forc_rh_grc                   (:)   => null() ! atmospheric relative humidity (%)
     real(r8), pointer :: forc_pco2_grc                 (:)   => null() ! CO2 partial pressure (Pa)
     real(r8), pointer :: forc_solad_grc                (:,:) => null() ! direct beam radiation (numrad) (vis=forc_sols , nir=forc_soll )
     real(r8), pointer :: forc_solai_grc                (:,:) => null() ! diffuse radiation (numrad) (vis=forc_solsd, nir=forc_solld)
     real(r8), pointer :: forc_solar_grc                (:)   => null() ! incident solar radiation
     real(r8), pointer :: forc_ndep_grc                 (:)   => null() ! nitrogen deposition rate (gN/m2/s)
     real(r8), pointer :: forc_pdep_grc                 (:)   => null() ! phosphorus deposition rate (gP/m2/s)
     real(r8), pointer :: forc_pc13o2_grc               (:)   => null() ! C13O2 partial pressure (Pa)
     real(r8), pointer :: forc_po2_grc                  (:)   => null() ! O2 partial pressure (Pa)
     real(r8), pointer :: forc_aer_grc                  (:,:) => null() ! aerosol deposition array
     real(r8), pointer :: forc_pch4_grc                 (:)   => null() ! CH4 partial pressure (Pa)

     real(r8), pointer :: forc_t_not_downscaled_grc     (:)   => null() ! not downscaled atm temperature (Kelvin)       
     real(r8), pointer :: forc_th_not_downscaled_grc    (:)   => null() ! not downscaled atm potential temperature (Kelvin)    
     real(r8), pointer :: forc_q_not_downscaled_grc     (:)   => null() ! not downscaled atm specific humidity (kg/kg)  
     real(r8), pointer :: forc_pbot_not_downscaled_grc  (:)   => null() ! not downscaled atm pressure (Pa)              
     real(r8), pointer :: forc_rho_not_downscaled_grc   (:)   => null() ! not downscaled atm density (kg/m**3)                      
     real(r8), pointer :: forc_rain_not_downscaled_grc  (:)   => null() ! not downscaled atm rain rate [mm/s]                       
     real(r8), pointer :: forc_snow_not_downscaled_grc  (:)   => null() ! not downscaled atm snow rate [mm/s]                       
     real(r8), pointer :: forc_lwrad_not_downscaled_grc (:)   => null() ! not downscaled atm downwrd IR longwave radiation (W/m**2) 

     ! atm->lnd downscaled
     real(r8), pointer :: forc_t_downscaled_col         (:)   => null() ! downscaled atm temperature (Kelvin)
     real(r8), pointer :: forc_th_downscaled_col        (:)   => null() ! downscaled atm potential temperature (Kelvin)
     real(r8), pointer :: forc_q_downscaled_col         (:)   => null() ! downscaled atm specific humidity (kg/kg)
     real(r8), pointer :: forc_pbot_downscaled_col      (:)   => null() ! downscaled atm pressure (Pa)
     real(r8), pointer :: forc_rho_downscaled_col       (:)   => null() ! downscaled atm density (kg/m**3)
     real(r8), pointer :: forc_rain_downscaled_col      (:)   => null() ! downscaled atm rain rate [mm/s]
     real(r8), pointer :: forc_snow_downscaled_col      (:)   => null() ! downscaled atm snow rate [mm/s]
     real(r8), pointer :: forc_lwrad_downscaled_col     (:)   => null() ! downscaled atm downwrd IR longwave radiation (W/m**2)

     !  rof->lnd
     real(r8), pointer :: forc_flood_grc                (:)   => null() ! rof flood (mm/s)
     real(r8), pointer :: volr_grc                      (:)   => null() ! rof volr total volume (m3)
     real(r8), pointer :: volrmch_grc                   (:)   => null() ! rof volr main channel (m3)

     ! anomaly forcing
     real(r8), pointer :: af_precip_grc                 (:)   => null() ! anomaly forcing 
     real(r8), pointer :: af_uwind_grc                  (:)   => null() ! anomaly forcing 
     real(r8), pointer :: af_vwind_grc                  (:)   => null() ! anomaly forcing 
     real(r8), pointer :: af_tbot_grc                   (:)   => null() ! anomaly forcing 
     real(r8), pointer :: af_pbot_grc                   (:)   => null() ! anomaly forcing 
     real(r8), pointer :: af_shum_grc                   (:)   => null() ! anomaly forcing 
     real(r8), pointer :: af_swdn_grc                   (:)   => null() ! anomaly forcing 
     real(r8), pointer :: af_lwdn_grc                   (:)   => null() ! anomaly forcing 
     real(r8), pointer :: bc_precip_grc                 (:)   => null() ! anomaly forcing - add bias correction

     ! time averaged quantities
     real(r8) , pointer :: fsd24_patch                  (:)   => null() ! patch 24hr average of direct beam radiation 
     real(r8) , pointer :: fsd240_patch                 (:)   => null() ! patch 240hr average of direct beam radiation 
     real(r8) , pointer :: fsi24_patch                  (:)   => null() ! patch 24hr average of diffuse beam radiation 
     real(r8) , pointer :: fsi240_patch                 (:)   => null() ! patch 240hr average of diffuse beam radiation 
     real(r8) , pointer :: prec365_patch                (:)   => null() ! patch 365-day running mean of tot. precipitation
     real(r8) , pointer :: prec60_patch                 (:)   => null() ! patch 60-day running mean of tot. precipitation (mm/s) 
     real(r8) , pointer :: prec10_patch                 (:)   => null() ! patch 10-day running mean of tot. precipitation (mm/s) 
     real(r8) , pointer :: prec24_patch                 (:)   => null() ! patch 24-hour running mean of tot. precipitation (mm/s) 
     real(r8) , pointer :: rh24_patch                   (:)   => null() ! patch 24-hour running mean of relative humidity
     real(r8) , pointer :: wind24_patch                 (:)   => null() ! patch 24-hour running mean of wind
     real(r8) , pointer :: t_mo_patch                   (:)   => null() ! patch 30-day average temperature (Kelvin)
     real(r8) , pointer :: t_mo_min_patch               (:)   => null() ! patch annual min of t_mo (Kelvin)

  contains
     procedure, public  :: Init
     procedure, private :: InitAllocate 
  end type atm2lnd_type

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(atm2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize atm2lnd derived type
    !
    ! !ARGUMENTS:
    class(atm2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ival  = 0.0_r8  ! initial value
    integer  :: begg, endg
    integer  :: begc, endc
    integer  :: begp, endp

    begg = bounds%begg; endg= bounds%endg
    begc = bounds%begc; endc= bounds%endc
    begp = bounds%begp; endp= bounds%endp

    allocate(this%forc_u_grc                    (begg:endg))        ; this%forc_u_grc                    (:)   = ival
    allocate(this%forc_v_grc                    (begg:endg))        ; this%forc_v_grc                    (:)   = ival
    allocate(this%forc_wind_grc                 (begg:endg))        ; this%forc_wind_grc                 (:)   = ival
    allocate(this%forc_rh_grc                   (begg:endg))        ; this%forc_rh_grc                   (:)   = ival
    allocate(this%forc_hgt_grc                  (begg:endg))        ; this%forc_hgt_grc                  (:)   = ival
    allocate(this%forc_hgt_u_grc                (begg:endg))        ; this%forc_hgt_u_grc                (:)   = ival
    allocate(this%forc_hgt_t_grc                (begg:endg))        ; this%forc_hgt_t_grc                (:)   = ival
    allocate(this%forc_hgt_q_grc                (begg:endg))        ; this%forc_hgt_q_grc                (:)   = ival
    allocate(this%forc_vp_grc                   (begg:endg))        ; this%forc_vp_grc                   (:)   = ival
    allocate(this%forc_pco2_grc                 (begg:endg))        ; this%forc_pco2_grc                 (:)   = ival
    allocate(this%forc_solad_grc                (begg:endg,numrad)) ; this%forc_solad_grc                (:,:) = ival
    allocate(this%forc_solai_grc                (begg:endg,numrad)) ; this%forc_solai_grc                (:,:) = ival
    allocate(this%forc_solar_grc                (begg:endg))        ; this%forc_solar_grc                (:)   = ival
    allocate(this%forc_ndep_grc                 (begg:endg))        ; this%forc_ndep_grc                 (:)   = ival
    allocate(this%forc_pdep_grc                 (begg:endg))        ; this%forc_pdep_grc                 (:)   = ival
    allocate(this%forc_pc13o2_grc               (begg:endg))        ; this%forc_pc13o2_grc               (:)   = ival
    allocate(this%forc_po2_grc                  (begg:endg))        ; this%forc_po2_grc                  (:)   = ival
    allocate(this%forc_aer_grc                  (begg:endg,14))     ; this%forc_aer_grc                  (:,:) = ival
    allocate(this%forc_pch4_grc                 (begg:endg))        ; this%forc_pch4_grc                 (:)   = ival

    ! atm->lnd not downscaled
    allocate(this%forc_t_not_downscaled_grc     (begg:endg))        ; this%forc_t_not_downscaled_grc     (:)   = ival
    allocate(this%forc_q_not_downscaled_grc     (begg:endg))        ; this%forc_q_not_downscaled_grc     (:)   = ival
    allocate(this%forc_pbot_not_downscaled_grc  (begg:endg))        ; this%forc_pbot_not_downscaled_grc  (:)   = ival
    allocate(this%forc_th_not_downscaled_grc    (begg:endg))        ; this%forc_th_not_downscaled_grc    (:)   = ival
    allocate(this%forc_rho_not_downscaled_grc   (begg:endg))        ; this%forc_rho_not_downscaled_grc   (:)   = ival
    allocate(this%forc_lwrad_not_downscaled_grc (begg:endg))        ; this%forc_lwrad_not_downscaled_grc (:)   = ival
    allocate(this%forc_rain_not_downscaled_grc  (begg:endg))        ; this%forc_rain_not_downscaled_grc  (:)   = ival
    allocate(this%forc_snow_not_downscaled_grc  (begg:endg))        ; this%forc_snow_not_downscaled_grc  (:)   = ival
    
    ! atm->lnd downscaled
    allocate(this%forc_t_downscaled_col         (begc:endc))        ; this%forc_t_downscaled_col         (:)   = ival
    allocate(this%forc_q_downscaled_col         (begc:endc))        ; this%forc_q_downscaled_col         (:)   = ival
    allocate(this%forc_pbot_downscaled_col      (begc:endc))        ; this%forc_pbot_downscaled_col      (:)   = ival
    allocate(this%forc_th_downscaled_col        (begc:endc))        ; this%forc_th_downscaled_col        (:)   = ival
    allocate(this%forc_rho_downscaled_col       (begc:endc))        ; this%forc_rho_downscaled_col       (:)   = ival
    allocate(this%forc_lwrad_downscaled_col     (begc:endc))        ; this%forc_lwrad_downscaled_col     (:)   = ival
    allocate(this%forc_rain_downscaled_col      (begc:endc))        ; this%forc_rain_downscaled_col      (:)   = ival
    allocate(this%forc_snow_downscaled_col      (begc:endc))        ; this%forc_snow_downscaled_col      (:)   = ival

    ! rof->lnd
    allocate(this%forc_flood_grc                (begg:endg))        ; this%forc_flood_grc                (:)   = ival
    allocate(this%volr_grc                      (begg:endg))        ; this%volr_grc                      (:)   = ival
    allocate(this%volrmch_grc                   (begg:endg))        ; this%volrmch_grc                   (:)   = ival

    ! anomaly forcing
    allocate(this%bc_precip_grc                 (begg:endg))        ; this%bc_precip_grc                 (:)   = ival
    allocate(this%af_precip_grc                 (begg:endg))        ; this%af_precip_grc                 (:)   = ival
    allocate(this%af_uwind_grc                  (begg:endg))        ; this%af_uwind_grc                  (:)   = ival
    allocate(this%af_vwind_grc                  (begg:endg))        ; this%af_vwind_grc                  (:)   = ival
    allocate(this%af_tbot_grc                   (begg:endg))        ; this%af_tbot_grc                   (:)   = ival
    allocate(this%af_pbot_grc                   (begg:endg))        ; this%af_pbot_grc                   (:)   = ival
    allocate(this%af_shum_grc                   (begg:endg))        ; this%af_shum_grc                   (:)   = ival
    allocate(this%af_swdn_grc                   (begg:endg))        ; this%af_swdn_grc                   (:)   = ival
    allocate(this%af_lwdn_grc                   (begg:endg))        ; this%af_lwdn_grc                   (:)   = ival

    allocate(this%fsd24_patch                   (begp:endp))        ; this%fsd24_patch                   (:)   = nan
    allocate(this%fsd240_patch                  (begp:endp))        ; this%fsd240_patch                  (:)   = nan
    allocate(this%fsi24_patch                   (begp:endp))        ; this%fsi24_patch                   (:)   = nan
    allocate(this%fsi240_patch                  (begp:endp))        ; this%fsi240_patch                  (:)   = nan
    allocate(this%prec10_patch                  (begp:endp))        ; this%prec10_patch                  (:)   = nan
    allocate(this%prec60_patch                  (begp:endp))        ; this%prec60_patch                  (:)   = nan
    allocate(this%prec365_patch                 (begp:endp))        ; this%prec365_patch                 (:)   = nan
    if (use_fates) then
       allocate(this%prec24_patch               (begp:endp))        ; this%prec24_patch                  (:)   = nan
       allocate(this%rh24_patch                 (begp:endp))        ; this%rh24_patch                    (:)   = nan
       allocate(this%wind24_patch               (begp:endp))        ; this%wind24_patch                  (:)   = nan
    end if
    allocate(this%t_mo_patch                    (begp:endp))        ; this%t_mo_patch               (:)   = nan
    allocate(this%t_mo_min_patch                (begp:endp))        ; this%t_mo_min_patch           (:)   = spval ! TODO - initialize this elsewhere

  end subroutine InitAllocate

end module atm2lndType
