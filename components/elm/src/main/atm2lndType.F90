module atm2lndType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle atm2lnd, lnd2atm mapping
  !
  ! !USES:
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod   , only : errMsg => shr_log_errMsg
  use shr_megan_mod , only : shr_megan_mechcomps_n
  use elm_varpar    , only : numrad, ndst, nlevgrnd !ndst = number of dust bins.
  use elm_varcon    , only : rair, grav, cpair, hfus, tfrz, spval
  use clm_varctl    , only : iulog, use_c13, use_cn, use_lch4, use_fates
  use seq_drydep_mod, only : n_drydep, drydep_method, DD_XLND
  use decompMod     , only : bounds_type
  use abortutils    , only : endrun
  use VegetationType     , only : veg_pp
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save
  !
  ! !PUBLIC DATA TYPES:
  !----------------------------------------------------
  ! atmosphere -> land variables structure
  !
  ! NOTE:
  ! IF there are forcing variables that are downscaled - then the
  ! non-downscaled versions SHOULD NOT be used in the code. Currently
  ! the non-downscaled versions are only used n a handful of places in
  ! the code (and needs to be used in lnd_import_export and the
  ! downscaling routines), but in general should NOT be used in new
  ! code. Instead use the datatype variables that have a _col suffix
  ! which gives the downscaled versions of these fields.
  !----------------------------------------------------
  type, public :: atm2lnd_type
      !DMR additions for CPL_BYPASS option
#ifdef CPL_BYPASS
      integer*2, pointer :: atm_input                (:,:,:,:) => null()  !Single-site meteorological input
      integer, pointer  :: loaded_bypassdata                   => null()
      real(r8), pointer :: add_offsets                     (:) => null()  !offsets for compressed met drivers
      real(r8), pointer :: scale_factors                   (:) => null()  !scale factors for compressed met drivers      
      integer(r8), pointer :: startyear_met                    => null()  !staring driver met year
      integer(r8), pointer :: endyear_met_spinup               => null()  !end driver met year for spinup cycle
      integer(r8), pointer :: endyear_met_trans                => null()  !end driver met year for transient simulation
      real(r8), pointer :: timeres                         (:) => null()  !time resolution of driver met (hours)
      real(r8), pointer :: var_offset                  (:,:,:) => null()  !correction offset for grid->site driver met (monthly)
      real(r8), pointer :: var_mult                    (:,:,:) => null()  !correction factor for grid->site driver met (monthly)
      integer,  pointer :: timelen                         (:) => null()  !length of input meteorology
      integer,  pointer :: timelen_spinup                  (:) => null()  !length of spinup meteorology
      integer,  pointer :: tindex                      (:,:,:) => null()  !current index for meteorolgoical data
      integer,  pointer :: metsource                           => null()  !Meteorogical source (0=Qian, 1=cruncep)
      real(r8), pointer :: npf                            (:)  => null()  !number of model timesteps per forcing timestep
      real(r8), pointer :: co2_input                   (:,:,:) => null()  !annual CO2 input data
      real(r8), pointer :: c13o2_input                 (:,:,:) => null()  !annual C13O2 input data
      integer, pointer :: ndepind                        (:,:) => null()  !annual nitrogen deposition data
      integer, pointer :: hdmind                         (:,:) => null()  !popluation density
      real(r8), pointer :: forc_hdm                      (:)   => null() 
      real(r8), pointer :: forc_lnfm                     (:)   => null()
      real(r8), pointer ::  hdm1                       (:,:,:) => null() 
      real(r8), pointer ::  hdm2                       (:,:,:) => null()
      real(r8), pointer ::  lnfm_all                   (:,:,:) => null()
      real(r8), pointer ::  lnfm                         (:,:) => null()
      real(r8), pointer ::  ndep1                      (:,:,:) => null()
      real(r8), pointer ::  ndep2                      (:,:,:) => null()
      real(r8), pointer ::  aerodata                 (:,:,:,:) => null()
#endif
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
     real(r8), pointer :: supply_grc                    (:)   => null() ! rof volr supply (mm/s)
     real(r8), pointer :: deficit_grc                   (:)   => null() ! rof volr deficit (mm/s)
	 
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
     procedure, private :: InitHistory  
     procedure, public  :: InitAccBuffer
     procedure, public  :: InitAccVars
     procedure, public  :: UpdateAccVars
     procedure, public  :: Restart

  end type atm2lnd_type
  !----------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(atm2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    
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
    real     :: ival_float = 0.0
    integer  :: ival_int = 0
    integer*2 :: ival_short = 0
    integer  :: begg, endg
    integer  :: begc, endc
    integer  :: begp, endp
    !------------------------------------------------------------------------

    begg = bounds%begg; endg= bounds%endg
    begc = bounds%begc; endc= bounds%endc
    begp = bounds%begp; endp= bounds%endp

    ! atm->lnd

    !DMR - variables added for CPL_BYPASS option 
#ifdef CPL_BYPASS
    allocate(this%timelen                            (1:14))        ; this%timelen                       (:)   = ival_int
    allocate(this%timelen_spinup                     (1:14))        ; this%timelen_spinup                (:)   = ival_int
    allocate(this%tindex               (begg:endg,1:14,1:2))        ; this%tindex                    (:,:,:)   = ival_int
    allocate(this%metsource                                )        ; this%metsource                           = ival_int   
    allocate(this%npf                                (1:14))        ; this%npf                           (:)   = ival
    !allocate(this%atm_input       (14,begg:endg,1,1:600000))        ; this%atm_input               (:,:,:,:)   = ival_short
    allocate(this%loaded_bypassdata                        )        ; this%loaded_bypassdata                   = 0
    allocate(this%add_offsets                        (1:14))        ; this%add_offsets                   (:)   = ival_float 
    allocate(this%scale_factors                      (1:14))        ; this%scale_factors                 (:)   = ival_float
    allocate(this%startyear_met                            )        ; this%startyear_met                       = ival_int
    allocate(this%endyear_met_spinup                       )        ; this%endyear_met_spinup                  = ival_int
    allocate(this%endyear_met_trans                        )        ; this%endyear_met_trans                   = ival_int
    allocate(this%timeres                            (1:14))        ; this%timeres                       (:)   = ival
    allocate(this%var_offset              (14,begg:endg,12))        ; this%var_offset                (:,:,:)   = ival
    allocate(this%var_mult                (14,begg:endg,12))        ; this%var_mult                  (:,:,:)   = ival
    allocate(this%co2_input                      (1,1,3000))        ; this%co2_input                 (:,:,:)   = ival    
    allocate(this%c13o2_input                    (1,1,3000))        ; this%c13o2_input               (:,:,:)   = ival
    allocate(this%ndepind                     (begg:endg,2))        ; this%ndepind                     (:,:)   = ival_int
    allocate(this%hdmind                      (begg:endg,2))        ; this%hdmind                      (:,:)   = ival_int
    allocate(this%forc_hdm                      (begg:endg))        ; this%forc_hdm                      (:)   = ival
    allocate(this%forc_lnfm                     (begg:endg))        ; this%forc_lnfm                     (:)   = ival
    allocate(this%hdm1                          (720,360,1))        ; this%hdm1                      (:,:,:)   = ival
    allocate(this%hdm2                          (720,360,1))        ; this%hdm2                      (:,:,:)   = ival
    allocate(this%lnfm                     (begg:endg,2920))        ; this%lnfm                        (:,:)   = ival 
    allocate(this%ndep1                          (144,96,1))        ; this%ndep1                     (:,:,:)   = ival
    allocate(this%ndep2                          (144,96,1))        ; this%ndep2                     (:,:,:)   = ival
    allocate(this%aerodata                   (14,144,96,14))        ; this%aerodata                (:,:,:,:)   = ival
    !END DMR
#endif
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
    allocate(this%supply_grc                    (begg:endg))        ; this%supply_grc                    (:)   = ival
    allocate(this%deficit_grc                   (begg:endg))        ; this%deficit_grc                   (:)   = ival

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

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !USES:
    use histFileMod, only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(atm2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer  :: begg, endg
    integer  :: begp, endp
    !---------------------------------------------------------------------

    begg = bounds%begg; endg= bounds%endg
    begp = bounds%begp; endp= bounds%endp

    this%forc_flood_grc(begg:endg) = spval
    call hist_addfld1d (fname='QFLOOD',  units='mm/s',  &
         avgflag='A', long_name='runoff from river flooding', &
         ptr_lnd=this%forc_flood_grc)

    this%volr_grc(begg:endg) = spval
    call hist_addfld1d (fname='VOLR',  units='m3',  &
         avgflag='A', long_name='river channel total water storage', &
         ptr_lnd=this%volr_grc)

    this%volrmch_grc(begg:endg) = spval
    call hist_addfld1d (fname='VOLRMCH',  units='m3',  &
         avgflag='A', long_name='river channel main channel water storage', &
         ptr_lnd=this%volrmch_grc)
		 
    this%supply_grc(begg:endg) = spval
    call hist_addfld1d (fname='SUPPLY',  units='mm/s',  &
         avgflag='A', long_name='runoff supply for land use', &
         ptr_lnd=this%supply_grc)
         
    this%deficit_grc(begg:endg) = spval
    call hist_addfld1d (fname='DEFICIT',  units='mm/s',  &
         avgflag='A', long_name='runoff supply deficit', &
         ptr_lnd=this%deficit_grc)

!    this%forc_wind_grc(begg:endg) = spval
!    call hist_addfld1d (fname='WIND', units='m/s',  &
!         avgflag='A', long_name='atmospheric wind velocity magnitude', &
!         ptr_lnd=this%forc_wind_grc)
    ! Rename of WIND for Urban intercomparision project
    call hist_addfld1d (fname='Wind', units='m/s',  &
         avgflag='A', long_name='atmospheric wind velocity magnitude', &
         ptr_gcell=this%forc_wind_grc, default = 'inactive')

!    this%forc_hgt_grc(begg:endg) = spval
!    call hist_addfld1d (fname='ZBOT', units='m',  &
!         avgflag='A', long_name='atmospheric reference height', &
!         ptr_lnd=this%forc_hgt_grc)

!    this%forc_solar_grc(begg:endg) = spval
!    call hist_addfld1d (fname='FSDS', units='W/m^2',  &
!         avgflag='A', long_name='atmospheric incident solar radiation', &
!         ptr_lnd=this%forc_solar_grc)

!    this%forc_pco2_grc(begg:endg) = spval
!    call hist_addfld1d (fname='PCO2', units='Pa',  &
!         avgflag='A', long_name='atmospheric partial pressure of CO2', &
!         ptr_lnd=this%forc_pco2_grc)

    this%forc_solar_grc(begg:endg) = spval
    call hist_addfld1d (fname='SWdown', units='W/m^2',  &
         avgflag='A', long_name='atmospheric incident solar radiation', &
         ptr_gcell=this%forc_solar_grc, default='inactive')

!    this%forc_rh_grc(begg:endg) = spval
!    call hist_addfld1d (fname='RH', units='%',  &
!         avgflag='A', long_name='atmospheric relative humidity', &
!         ptr_gcell=this%forc_rh_grc, default='inactive')

!    if (use_lch4) then
!       this%forc_pch4_grc(begg:endg) = spval
!       call hist_addfld1d (fname='PCH4', units='Pa',  &
!            avgflag='A', long_name='atmospheric partial pressure of CH4', &
!            ptr_lnd=this%forc_pch4_grc)
!    end if

    this%forc_t_not_downscaled_grc(begg:endg) = spval
    call hist_addfld1d (fname='Tair', units='K',  &
         avgflag='A', long_name='atmospheric air temperature', &
         ptr_gcell=this%forc_t_not_downscaled_grc, default='inactive')

    this%forc_pbot_not_downscaled_grc(begg:endg) = spval
    call hist_addfld1d (fname='PSurf', units='Pa',  &
         avgflag='A', long_name='surface pressure', &
         ptr_gcell=this%forc_pbot_not_downscaled_grc, default='inactive')

    this%forc_rain_not_downscaled_grc(begg:endg) = spval
    call hist_addfld1d (fname='Rainf', units='mm/s',  &
         avgflag='A', long_name='atmospheric rain', &
         ptr_gcell=this%forc_rain_not_downscaled_grc, default='inactive')

    this%forc_lwrad_not_downscaled_grc(begg:endg) = spval
    call hist_addfld1d (fname='LWdown', units='W/m^2',  &
         avgflag='A', long_name='atmospheric longwave radiation', &
         ptr_gcell=this%forc_lwrad_not_downscaled_grc, default='inactive')

!    this%forc_rain_not_downscaled_grc(begg:endg) = spval
!    call hist_addfld1d (fname='RAIN', units='mm/s',  &
!         avgflag='A', long_name='atmospheric rain', &
!         ptr_lnd=this%forc_rain_not_downscaled_grc)

!    this%forc_snow_not_downscaled_grc(begg:endg) = spval
!    call hist_addfld1d (fname='SNOW', units='mm/s',  &
!         avgflag='A', long_name='atmospheric snow', &
!         ptr_lnd=this%forc_snow_not_downscaled_grc)

!   this%forc_t_not_downscaled_grc(begg:endg) = spval
!   call hist_addfld1d (fname='TBOT', units='K',  &
!        avgflag='A', long_name='atmospheric air temperature', &
!        ptr_lnd=this%forc_t_not_downscaled_grc)

!    this%forc_th_not_downscaled_grc(begg:endg) = spval
!    call hist_addfld1d (fname='THBOT', units='K',  &
!         avgflag='A', long_name='atmospheric air potential temperature', &
!         ptr_lnd=this%forc_th_not_downscaled_grc)

!    this%forc_q_not_downscaled_grc(begg:endg) = spval
!    call hist_addfld1d (fname='QBOT', units='kg/kg',  &
!         avgflag='A', long_name='atmospheric specific humidity', &
!         ptr_lnd=this%forc_q_not_downscaled_grc)
    ! Rename of QBOT for Urban intercomparision project
    this%forc_q_not_downscaled_grc(begg:endg) = spval
    call hist_addfld1d (fname='Qair', units='kg/kg',  &
         avgflag='A', long_name='atmospheric specific humidity', &
         ptr_lnd=this%forc_q_not_downscaled_grc, default='inactive')

!    this%forc_lwrad_not_downscaled_grc(begg:endg) = spval
!    call hist_addfld1d (fname='FLDS', units='W/m^2',  &
!         avgflag='A', long_name='atmospheric longwave radiation', &
!         ptr_lnd=this%forc_lwrad_not_downscaled_grc)

!    this%forc_pbot_not_downscaled_grc(begg:endg) = spval
!    call hist_addfld1d (fname='PBOT', units='Pa',  &
!         avgflag='A', long_name='atmospheric pressure', &
!         ptr_lnd=this%forc_pbot_not_downscaled_grc)

#ifdef CPL_BYPASS
   call hist_addfld1d (fname='HDM', units='counts/km^2',      &
         avgflag='A', long_name='human population density',   &
         ptr_lnd=this%forc_hdm, default='inactive')
   
    call hist_addfld1d (fname='LNFM', units='counts/km^2/hr',  &
         avgflag='A', long_name='Lightning frequency',        &
         ptr_lnd=this%forc_lnfm, default='inactive')
#endif

    ! Time averaged quantities
    this%fsi24_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSI24', units='K',  &
         avgflag='A', long_name='indirect radiation (last 24hrs)', &
         ptr_patch=this%fsi24_patch, default='inactive')

    this%fsi240_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSI240', units='K',  &
         avgflag='A', long_name='indirect radiation (last 240hrs)', &
         ptr_patch=this%fsi240_patch, default='inactive')

    this%fsd24_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSD24', units='K',  &
         avgflag='A', long_name='direct radiation (last 24hrs)', &
         ptr_patch=this%fsd24_patch, default='inactive')

    this%fsd240_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSD240', units='K',  &
         avgflag='A', long_name='direct radiation (last 240hrs)', &
         ptr_patch=this%fsd240_patch, default='inactive')

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitAccBuffer (this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize accumulation buffer for all required module accumulated fields
    ! This routine set defaults values that are then overwritten by the
    ! restart file for restart or branch runs
    !
    ! !USES 
    use elm_varcon  , only : spval
    use accumulMod  , only : init_accum_field
    !
    ! !ARGUMENTS:
    class(atm2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !---------------------------------------------------------------------

    this%fsd24_patch(bounds%begp:bounds%endp) = spval
    call init_accum_field (name='FSD24', units='W/m2',                                             &
         desc='24hr average of direct solar radiation',  accum_type='runmean', accum_period=-1,    &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    this%fsd240_patch(bounds%begp:bounds%endp) = spval
    call init_accum_field (name='FSD240', units='W/m2',                                            &
         desc='240hr average of direct solar radiation',  accum_type='runmean', accum_period=-10,  &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    this%fsi24_patch(bounds%begp:bounds%endp) = spval
    call init_accum_field (name='FSI24', units='W/m2',                                             &
         desc='24hr average of diffuse solar radiation',  accum_type='runmean', accum_period=-1,   &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    this%fsi240_patch(bounds%begp:bounds%endp) = spval
    call init_accum_field (name='FSI240', units='W/m2',                                            &
         desc='240hr average of diffuse solar radiation',  accum_type='runmean', accum_period=-10, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    if (use_cn) then
       call init_accum_field (name='PREC10', units='MM H2O/S', &
            desc='10-day running mean of total precipitation', accum_type='runmean', accum_period=-10, &
            subgrid_type='pft', numlev=1, init_value=0._r8)

       call init_accum_field (name='PREC60', units='MM H2O/S', &
            desc='60-day running mean of total precipitation', accum_type='runmean', accum_period=-60, &
            subgrid_type='pft', numlev=1, init_value=0._r8)
    end if

    if ( use_fates ) then
       call init_accum_field (name='PREC24', units='m', &
            desc='24hr sum of precipitation', accum_type='runmean', accum_period=-1, &
            subgrid_type='pft', numlev=1, init_value=0._r8)

       ! Fudge - this neds to be initialized from the restat file eventually. 
       call init_accum_field (name='RH24', units='m', &
            desc='24hr average of RH', accum_type='runmean', accum_period=-1, &
            subgrid_type='pft', numlev=1, init_value=100._r8) 

       call init_accum_field (name='WIND24', units='m', &
            desc='24hr average of wind', accum_type='runmean', accum_period=-1, &
            subgrid_type='pft', numlev=1, init_value=0._r8)
    end if

  end subroutine InitAccBuffer

  !-----------------------------------------------------------------------
  subroutine InitAccVars(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module variables that are associated with
    ! time accumulated fields. This routine is called for both an initial run
    ! and a restart run (and must therefore must be called after the restart file 
    ! is read in and the accumulation buffer is obtained)
    !
    ! !USES 
    use accumulMod       , only : extract_accum_field
    use clm_time_manager , only : get_nstep
    !
    ! !ARGUMENTS:
    class(atm2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer  :: begp, endp
    integer  :: nstep
    integer  :: ier
    real(r8), pointer :: rbufslp(:)  ! temporary
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp

    ! Allocate needed dynamic memory for single level pft field
    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(iulog,*)' in '
       call endrun(msg="extract_accum_hist allocation error for rbufslp"//&
            errMsg(__FILE__, __LINE__))
    endif

    ! Determine time step
    nstep = get_nstep()

    call extract_accum_field ('FSD24', rbufslp, nstep)
    this%fsd24_patch(begp:endp) = rbufslp(begp:endp)

    call extract_accum_field ('FSD240', rbufslp, nstep)
    this%fsd240_patch(begp:endp) = rbufslp(begp:endp)

    call extract_accum_field ('FSI24', rbufslp, nstep)
    this%fsi24_patch(begp:endp) = rbufslp(begp:endp)

    call extract_accum_field ('FSI240', rbufslp, nstep)
    this%fsi240_patch(begp:endp) = rbufslp(begp:endp)

    if (use_cn) then
       call extract_accum_field ('PREC10', rbufslp, nstep)
       this%prec10_patch(begp:endp) = rbufslp(begp:endp)

       call extract_accum_field ('PREC60', rbufslp, nstep)
       this%prec60_patch(begp:endp) = rbufslp(begp:endp)
    end if

    if (use_fates) then
       call extract_accum_field ('PREC24', rbufslp, nstep)
       this%prec24_patch(begp:endp) = rbufslp(begp:endp)

       call extract_accum_field ('RH24', rbufslp, nstep)
       this%rh24_patch(begp:endp) = rbufslp(begp:endp)

       call extract_accum_field ('WIND24', rbufslp, nstep)
       this%wind24_patch(begp:endp) = rbufslp(begp:endp)
    end if

    deallocate(rbufslp)

  end subroutine InitAccVars

  !-----------------------------------------------------------------------
  subroutine UpdateAccVars (this, bounds)
    !
    ! USES
    use clm_time_manager, only : get_nstep
    use accumulMod      , only : update_accum_field, extract_accum_field
    !
    ! !ARGUMENTS:
    class(atm2lnd_type)                 :: this
    type(bounds_type)      , intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: g,c,p                     ! indices
    integer :: dtime                     ! timestep size [seconds]
    integer :: nstep                     ! timestep number
    integer :: ier                       ! error status
    integer :: begp, endp
    real(r8), pointer :: rbufslp(:)      ! temporary single level - pft level
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp

    nstep = get_nstep()

    ! Allocate needed dynamic memory for single level pft field

    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(iulog,*)'update_accum_hist allocation error for rbuf1dp'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    ! Accumulate and extract forc_solad24 & forc_solad240 
    do p = begp,endp
       g = veg_pp%gridcell(p)
       rbufslp(p) = this%forc_solad_grc(g,1)
    end do
    call update_accum_field  ('FSD240', rbufslp               , nstep)
    call extract_accum_field ('FSD240', this%fsd240_patch     , nstep)
    call update_accum_field  ('FSD24' , rbufslp               , nstep)
    call extract_accum_field ('FSD24' , this%fsd24_patch      , nstep)

    ! Accumulate and extract forc_solai24 & forc_solai240 
    do p = begp,endp
       g = veg_pp%gridcell(p)
       rbufslp(p) = this%forc_solai_grc(g,1)
    end do
    call update_accum_field  ('FSI24' , rbufslp               , nstep)
    call extract_accum_field ('FSI24' , this%fsi24_patch      , nstep)
    call update_accum_field  ('FSI240', rbufslp               , nstep)
    call extract_accum_field ('FSI240', this%fsi240_patch     , nstep)

    do p = begp,endp
       c = veg_pp%column(p)
       rbufslp(p) = this%forc_rain_downscaled_col(c) + this%forc_snow_downscaled_col(c)
    end do

    if (use_cn) then
       ! Accumulate and extract PREC60 (accumulates total precipitation as 60-day running mean)
       call update_accum_field  ('PREC60', rbufslp, nstep)
       call extract_accum_field ('PREC60', this%prec60_patch, nstep)

       ! Accumulate and extract PREC10 (accumulates total precipitation as 10-day running mean)
       call update_accum_field  ('PREC10', rbufslp, nstep)
       call extract_accum_field ('PREC10', this%prec10_patch, nstep)
    end if

    if (use_fates) then
       call update_accum_field  ('PREC24', rbufslp, nstep)
       call extract_accum_field ('PREC24', this%prec24_patch, nstep)

       do p = bounds%begp,bounds%endp
          c = veg_pp%column(p)
          g = veg_pp%gridcell(p)
          rbufslp(p) = this%forc_wind_grc(g) 
       end do
       call update_accum_field  ('WIND24', rbufslp, nstep)
       call extract_accum_field ('WIND24', this%wind24_patch, nstep)

       do p = bounds%begp,bounds%endp
          c = veg_pp%column(p)
          g = veg_pp%gridcell(p)
          rbufslp(p) = this%forc_rh_grc(g) 
       end do
       call update_accum_field  ('RH24', rbufslp, nstep)
       call extract_accum_field ('RH24', this%rh24_patch, nstep)
    end if

    deallocate(rbufslp)

  end subroutine UpdateAccVars

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !USES:
    use restUtilMod
    use ncdio_pio
    !
    ! !ARGUMENTS:
    class(atm2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds  
    type(file_desc_t), intent(inout) :: ncid   
    character(len=*) , intent(in)    :: flag   
    !
    ! !LOCAL VARIABLES:
    logical            :: readvar 
    !------------------------------------------------------------------------

    call restartvar(ncid=ncid, flag=flag, varname='qflx_floodg', xtype=ncd_double, &
         dim1name='gridcell', &
         long_name='flood water flux', units='mm/s', &
         interpinic_flag='skip', readvar=readvar, data=this%forc_flood_grc)
    if (flag == 'read' .and. .not. readvar) then
       ! initial run, readvar=readvar, not restart: initialize flood to zero
       this%forc_flood_grc = 0._r8
    endif

  end subroutine Restart

end module atm2lndType
