module UrbanFluxesMod

  !----------------------------------------------------------------------- 
  ! !DESCRIPTION: 
  ! Calculate solar and longwave radiation, and turbulent fluxes for urban landunit
  !
  ! !USES:
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use shr_sys_mod          , only : shr_sys_flush 
  use shr_log_mod          , only : errMsg => shr_log_errMsg
  use decompMod            , only : bounds_type
  use clm_varpar           , only : numrad
  use clm_varcon           , only : isecspday, degpsec, namel
  use clm_varctl           , only : iulog
  use abortutils           , only : endrun  
  use UrbanParamsType      , only : urbanparams_type
  use UrbanParamsType      , only : urban_wasteheat_on, urban_hac_on, urban_hac 
  use atm2lndType          , only : atm2lnd_type
  use SoilStateType        , only : soilstate_type
  use TemperatureType      , only : temperature_type
  use WaterstateType       , only : waterstate_type
  use FrictionVelocityType , only : frictionvel_type
  use EnergyFluxType       , only : energyflux_type
  use WaterfluxType        , only : waterflux_type
  use GridcellType         , only : grc                
  use LandunitType         , only : lun                
  use ColumnType           , only : col                
  use PatchType            , only : pft                
  use SurfaceResistanceMod , only : do_soilevap_beta
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: UrbanFluxes       ! Urban physics - turbulent fluxes
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine UrbanFluxes (bounds, num_nourbanl, filter_nourbanl,                        &
       num_urbanl, filter_urbanl, num_urbanc, filter_urbanc, num_urbanp, filter_urbanp, &
       atm2lnd_vars, urbanparams_vars, soilstate_vars, temperature_vars,                &
       waterstate_vars, frictionvel_vars, energyflux_vars, waterflux_vars) 
    !
    ! !DESCRIPTION: 
    ! Turbulent and momentum fluxes from urban canyon (consisting of roof, sunwall, 
    ! shadewall, pervious and impervious road).

    ! !USES:
    use clm_varcon          , only : cpair, vkc, spval, grav, pondmx_urban, rpi, rgas
    use clm_varcon          , only : ht_wasteheat_factor, ac_wasteheat_factor, wasteheat_limit
    use column_varcon       , only : icol_shadewall, icol_road_perv, icol_road_imperv
    use column_varcon       , only : icol_roof, icol_sunwall
    use filterMod           , only : filter
    use FrictionVelocityMod , only : FrictionVelocity, MoninObukIni
    use QSatMod             , only : QSat
    use clm_varpar          , only : maxpatch_urb, nlevurb, nlevgrnd
    use clm_time_manager    , only : get_curr_date, get_step_size, get_nstep
    use clm_varctl          , only : use_vsfm
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds    
    integer                , intent(in)    :: num_nourbanl       ! number of non-urban landunits in clump
    integer                , intent(in)    :: filter_nourbanl(:) ! non-urban landunit filter
    integer                , intent(in)    :: num_urbanl         ! number of urban landunits in clump
    integer                , intent(in)    :: filter_urbanl(:)   ! urban landunit filter
    integer                , intent(in)    :: num_urbanc         ! number of urban columns in clump
    integer                , intent(in)    :: filter_urbanc(:)   ! urban column filter
    integer                , intent(in)    :: num_urbanp         ! number of urban patches in clump
    integer                , intent(in)    :: filter_urbanp(:)   ! urban pft filter
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_vars
    type(urbanparams_type) , intent(in)    :: urbanparams_vars
    type(soilstate_type)   , intent(inout) :: soilstate_vars
    type(temperature_type) , intent(inout) :: temperature_vars
    type(waterstate_type)  , intent(inout) :: waterstate_vars
    type(frictionvel_type) , intent(inout) :: frictionvel_vars
    type(waterflux_type)   , intent(inout) :: waterflux_vars
    type(energyflux_type)  , intent(inout) :: energyflux_vars
    !
    ! !LOCAL VARIABLES:
    character(len=*), parameter :: sub="UrbanFluxes"
    integer  :: fp,fc,fl,f,p,c,l,g,j,pi,i     ! indices

    real(r8) :: canyontop_wind(bounds%begl:bounds%endl)              ! wind at canyon top (m/s) 
    real(r8) :: canyon_u_wind(bounds%begl:bounds%endl)               ! u-component of wind speed inside canyon (m/s)
    real(r8) :: canyon_wind(bounds%begl:bounds%endl)                 ! net wind speed inside canyon (m/s)
    real(r8) :: canyon_resistance(bounds%begl:bounds%endl)           ! resistance to heat and moisture transfer from canyon road/walls to canyon air (s/m)

    real(r8) :: ur(bounds%begl:bounds%endl)                          ! wind speed at reference height (m/s)
    real(r8) :: ustar(bounds%begl:bounds%endl)                       ! friction velocity (m/s)
    real(r8) :: ramu(bounds%begl:bounds%endl)                        ! aerodynamic resistance (s/m)
    real(r8) :: rahu(bounds%begl:bounds%endl)                        ! thermal resistance (s/m)
    real(r8) :: rawu(bounds%begl:bounds%endl)                        ! moisture resistance (s/m)
    real(r8) :: temp1(bounds%begl:bounds%endl)                       ! relation for potential temperature profile
    real(r8) :: temp12m(bounds%begl:bounds%endl)                     ! relation for potential temperature profile applied at 2-m
    real(r8) :: temp2(bounds%begl:bounds%endl)                       ! relation for specific humidity profile
    real(r8) :: temp22m(bounds%begl:bounds%endl)                     ! relation for specific humidity profile applied at 2-m
    real(r8) :: thm_g(bounds%begl:bounds%endl)                       ! intermediate variable (forc_t+0.0098*forc_hgt_t)
    real(r8) :: thv_g(bounds%begl:bounds%endl)                       ! virtual potential temperature (K)
    real(r8) :: dth(bounds%begl:bounds%endl)                         ! diff of virtual temp. between ref. height and surface
    real(r8) :: dqh(bounds%begl:bounds%endl)                         ! diff of humidity between ref. height and surface
    real(r8) :: zldis(bounds%begl:bounds%endl)                       ! reference height "minus" zero displacement height (m)
    real(r8) :: um(bounds%begl:bounds%endl)                          ! wind speed including the stablity effect (m/s)
    real(r8) :: obu(bounds%begl:bounds%endl)                         ! Monin-Obukhov length (m)
    real(r8) :: taf_numer(bounds%begl:bounds%endl)                   ! numerator of taf equation (K m/s)
    real(r8) :: taf_denom(bounds%begl:bounds%endl)                   ! denominator of taf equation (m/s)
    real(r8) :: qaf_numer(bounds%begl:bounds%endl)                   ! numerator of qaf equation (kg m/kg s)
    real(r8) :: qaf_denom(bounds%begl:bounds%endl)                   ! denominator of qaf equation (m/s)
    real(r8) :: wtas(bounds%begl:bounds%endl)                        ! sensible heat conductance for urban air to atmospheric air (m/s)
    real(r8) :: wtaq(bounds%begl:bounds%endl)                        ! latent heat conductance for urban air to atmospheric air (m/s)
    real(r8) :: wts_sum(bounds%begl:bounds%endl)                     ! sum of wtas, wtus_roof, wtus_road_perv, wtus_road_imperv, wtus_sunwall, wtus_shadewall
    real(r8) :: wtq_sum(bounds%begl:bounds%endl)                     ! sum of wtaq, wtuq_roof, wtuq_road_perv, wtuq_road_imperv, wtuq_sunwall, wtuq_shadewall
    real(r8) :: beta(bounds%begl:bounds%endl)                        ! coefficient of convective velocity
    real(r8) :: zii(bounds%begl:bounds%endl)                         ! convective boundary layer height (m)
    real(r8) :: fm(bounds%begl:bounds%endl)                          ! needed for BGC only to diagnose 10m wind speed
    real(r8) :: wtus(bounds%begc:bounds%endc)                        ! sensible heat conductance for urban columns (scaled) (m/s)
    real(r8) :: wtuq(bounds%begc:bounds%endc)                        ! latent heat conductance for urban columns (scaled) (m/s)
    integer  :: iter                                                 ! iteration index
    real(r8) :: dthv                                                 ! diff of vir. poten. temp. between ref. height and surface
    real(r8) :: tstar                                                ! temperature scaling parameter
    real(r8) :: qstar                                                ! moisture scaling parameter
    real(r8) :: thvstar                                              ! virtual potential temperature scaling parameter
    real(r8) :: wtus_roof(bounds%begl:bounds%endl)                   ! sensible heat conductance for roof (scaled) (m/s)
    real(r8) :: wtuq_roof(bounds%begl:bounds%endl)                   ! latent heat conductance for roof (scaled) (m/s)
    real(r8) :: wtus_road_perv(bounds%begl:bounds%endl)              ! sensible heat conductance for pervious road (scaled) (m/s)
    real(r8) :: wtuq_road_perv(bounds%begl:bounds%endl)              ! latent heat conductance for pervious road (scaled) (m/s)
    real(r8) :: wtus_road_imperv(bounds%begl:bounds%endl)            ! sensible heat conductance for impervious road (scaled) (m/s)
    real(r8) :: wtuq_road_imperv(bounds%begl:bounds%endl)            ! latent heat conductance for impervious road (scaled) (m/s)
    real(r8) :: wtus_sunwall(bounds%begl:bounds%endl)                ! sensible heat conductance for sunwall (scaled) (m/s)
    real(r8) :: wtuq_sunwall(bounds%begl:bounds%endl)                ! latent heat conductance for sunwall (scaled) (m/s)
    real(r8) :: wtus_shadewall(bounds%begl:bounds%endl)              ! sensible heat conductance for shadewall (scaled) (m/s)
    real(r8) :: wtuq_shadewall(bounds%begl:bounds%endl)              ! latent heat conductance for shadewall (scaled) (m/s)
    real(r8) :: wtus_roof_unscl(bounds%begl:bounds%endl)             ! sensible heat conductance for roof (not scaled) (m/s)
    real(r8) :: wtuq_roof_unscl(bounds%begl:bounds%endl)             ! latent heat conductance for roof (not scaled) (m/s)
    real(r8) :: wtus_road_perv_unscl(bounds%begl:bounds%endl)        ! sensible heat conductance for pervious road (not scaled) (m/s)
    real(r8) :: wtuq_road_perv_unscl(bounds%begl:bounds%endl)        ! latent heat conductance for pervious road (not scaled) (m/s)
    real(r8) :: wtus_road_imperv_unscl(bounds%begl:bounds%endl)      ! sensible heat conductance for impervious road (not scaled) (m/s)
    real(r8) :: wtuq_road_imperv_unscl(bounds%begl:bounds%endl)      ! latent heat conductance for impervious road (not scaled) (m/s)
    real(r8) :: wtus_sunwall_unscl(bounds%begl:bounds%endl)          ! sensible heat conductance for sunwall (not scaled) (m/s)
    real(r8) :: wtuq_sunwall_unscl(bounds%begl:bounds%endl)          ! latent heat conductance for sunwall (not scaled) (m/s)
    real(r8) :: wtus_shadewall_unscl(bounds%begl:bounds%endl)        ! sensible heat conductance for shadewall (not scaled) (m/s)
    real(r8) :: wtuq_shadewall_unscl(bounds%begl:bounds%endl)        ! latent heat conductance for shadewall (not scaled) (m/s)
    real(r8) :: t_sunwall_innerl(bounds%begl:bounds%endl)            ! temperature of inner layer of sunwall (K)
    real(r8) :: t_shadewall_innerl(bounds%begl:bounds%endl)          ! temperature of inner layer of shadewall (K)
    real(r8) :: t_roof_innerl(bounds%begl:bounds%endl)               ! temperature of inner layer of roof (K)
    real(r8) :: lngth_roof                                           ! length of roof (m)
    real(r8) :: wc                                                   ! convective velocity (m/s)
    real(r8) :: zeta                                                 ! dimensionless height used in Monin-Obukhov theory 
    real(r8) :: eflx_sh_grnd_scale(bounds%begp:bounds%endp)          ! scaled sensible heat flux from ground (W/m**2) [+ to atm] 
    real(r8) :: qflx_evap_soi_scale(bounds%begp:bounds%endp)         ! scaled soil evaporation (mm H2O/s) (+ = to atm) 
    real(r8) :: eflx_wasteheat_roof(bounds%begl:bounds%endl)         ! sensible heat flux from urban heating/cooling sources of waste heat for roof (W/m**2)
    real(r8) :: eflx_wasteheat_sunwall(bounds%begl:bounds%endl)      ! sensible heat flux from urban heating/cooling sources of waste heat for sunwall (W/m**2)
    real(r8) :: eflx_wasteheat_shadewall(bounds%begl:bounds%endl)    ! sensible heat flux from urban heating/cooling sources of waste heat for shadewall (W/m**2)
    real(r8) :: eflx_heat_from_ac_roof(bounds%begl:bounds%endl)      ! sensible heat flux put back into canyon due to heat removal by AC for roof (W/m**2)
    real(r8) :: eflx_heat_from_ac_sunwall(bounds%begl:bounds%endl)   ! sensible heat flux put back into canyon due to heat removal by AC for sunwall (W/m**2)
    real(r8) :: eflx_heat_from_ac_shadewall(bounds%begl:bounds%endl) ! sensible heat flux put back into canyon due to heat removal by AC for shadewall (W/m**2)
    real(r8) :: eflx(bounds%begl:bounds%endl)                        ! total sensible heat flux for error check (W/m**2)
    real(r8) :: qflx(bounds%begl:bounds%endl)                        ! total water vapor flux for error check (kg/m**2/s)
    real(r8) :: eflx_scale(bounds%begl:bounds%endl)                  ! sum of scaled sensible heat fluxes for urban columns for error check (W/m**2)
    real(r8) :: qflx_scale(bounds%begl:bounds%endl)                  ! sum of scaled water vapor fluxes for urban columns for error check (kg/m**2/s)
    real(r8) :: eflx_err(bounds%begl:bounds%endl)                    ! sensible heat flux error (W/m**2)
    real(r8) :: qflx_err(bounds%begl:bounds%endl)                    ! water vapor flux error (kg/m**2/s)
    real(r8) :: fwet_roof                                            ! fraction of roof surface that is wet (-)
    real(r8) :: fwet_road_imperv                                     ! fraction of impervious road surface that is wet (-)
    integer  :: local_secp1(bounds%begl:bounds%endl)                 ! seconds into current date in local time (sec)
    real(r8) :: dtime                                                ! land model time step (sec)
    integer  :: year,month,day,secs                                  ! calendar info for current time step
    logical  :: found                                                ! flag in search loop
    integer  :: indexl                                               ! index of first found in search loop
    integer  :: nstep                                                ! time step number
    real(r8) :: e_ref2m                                              ! 2 m height surface saturated vapor pressure [Pa]
    real(r8) :: de2mdT                                               ! derivative of 2 m height surface saturated vapor pressure on t_ref2m
    real(r8) :: qsat_ref2m                                           ! 2 m height surface saturated specific humidity [kg/kg]
    real(r8) :: dqsat2mdT                                            ! derivative of 2 m height surface saturated specific humidity on t_ref2m
    !
    real(r8), parameter :: lapse_rate = 0.0098_r8 ! Dry adiabatic lapse rate (K/m)
    integer , parameter  :: niters = 3            ! maximum number of iterations for surface temperature
    !-----------------------------------------------------------------------

    associate(                                                                & 
         snl                 =>   col%snl                                   , & ! Input:  [integer  (:)   ]  number of snow layers                              
         ctype               =>   col%itype                                 , & ! Input:  [integer  (:)   ]  column type                                        
         z_0_town            =>   lun%z_0_town                              , & ! Input:  [real(r8) (:)   ]  momentum roughness length of urban landunit (m)   
         z_d_town            =>   lun%z_d_town                              , & ! Input:  [real(r8) (:)   ]  displacement height of urban landunit (m)         
         ht_roof             =>   lun%ht_roof                               , & ! Input:  [real(r8) (:)   ]  height of urban roof (m)                          
         wtlunit_roof        =>   lun%wtlunit_roof                          , & ! Input:  [real(r8) (:)   ]  weight of roof with respect to landunit           
         canyon_hwr          =>   lun%canyon_hwr                            , & ! Input:  [real(r8) (:)   ]  ratio of building height to street width          
         wtroad_perv         =>   lun%wtroad_perv                           , & ! Input:  [real(r8) (:)   ]  weight of pervious road wrt total road            

         forc_t              =>   atm2lnd_vars%forc_t_not_downscaled_grc    , & ! Input:  [real(r8) (:)   ]  atmospheric temperature (K)                       
         forc_th             =>   atm2lnd_vars%forc_th_not_downscaled_grc   , & ! Input:  [real(r8) (:)   ]  atmospheric potential temperature (K)             
         forc_rho            =>   atm2lnd_vars%forc_rho_not_downscaled_grc  , & ! Input:  [real(r8) (:)   ]  density (kg/m**3)                                 
         forc_q              =>   atm2lnd_vars%forc_q_not_downscaled_grc    , & ! Input:  [real(r8) (:)   ]  atmospheric specific humidity (kg/kg)             
         forc_pbot           =>   atm2lnd_vars%forc_pbot_not_downscaled_grc , & ! Input:  [real(r8) (:)   ]  atmospheric pressure (Pa)                         
         forc_u              =>   atm2lnd_vars%forc_u_grc                   , & ! Input:  [real(r8) (:)   ]  atmospheric wind speed in east direction (m/s)    
         forc_v              =>   atm2lnd_vars%forc_v_grc                   , & ! Input:  [real(r8) (:)   ]  atmospheric wind speed in north direction (m/s)   

         wind_hgt_canyon     =>   urbanparams_vars%wind_hgt_canyon          , & ! Input:  [real(r8) (:)   ]  height above road at which wind in canyon is to be computed (m)
         eflx_traffic_factor =>   urbanparams_vars%eflx_traffic_factor      , & ! Input:  [real(r8) (:)   ]  multiplicative urban traffic factor for sensible heat flux

         rootr_road_perv     =>   soilstate_vars%rootr_road_perv_col        , & ! Input:  [real(r8) (:,:) ]  effective fraction of roots in each soil layer for urban pervious road
         soilalpha_u         =>   soilstate_vars%soilalpha_u_col            , & ! Input:  [real(r8) (:)   ]  Urban factor that reduces ground saturated specific humidity (-)
         soilbeta            =>   soilstate_vars%soilbeta_col               , & ! Input:  [real(r8) (:)   ]  soil wetness relative to field capacity
         rootr               =>   soilstate_vars%rootr_patch                , & ! Output: [real(r8) (:,:) ]  effective fraction of roots in each soil layer  

         t_grnd              =>   temperature_vars%t_grnd_col               , & ! Input:  [real(r8) (:)   ]  ground surface temperature (K)                    
         t_soisno            =>   temperature_vars%t_soisno_col             , & ! Input:  [real(r8) (:,:) ]  soil temperature (K)                            
         t_ref2m             =>   temperature_vars%t_ref2m_patch            , & ! Output: [real(r8) (:)   ]  2 m height surface air temperature (K)            
         t_ref2m_u           =>   temperature_vars%t_ref2m_u_patch          , & ! Output: [real(r8) (:)   ]  Urban 2 m height surface air temperature (K)     
         t_veg               =>   temperature_vars%t_veg_patch              , & ! Output: [real(r8) (:)   ]  vegetation temperature (K)                        
         t_building          =>   temperature_vars%t_building_lun           , & ! Output: [real(r8) (:)   ]  internal building temperature (K)                 
         taf                 =>   temperature_vars%taf_lun                  , & ! Output: [real(r8) (:)   ]  urban canopy air temperature (K)                  

         frac_sno            =>   waterstate_vars%frac_sno_col              , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)       
         snow_depth          =>   waterstate_vars%snow_depth_col            , & ! Input:  [real(r8) (:)   ]  snow height (m)                                   
         dqgdT               =>   waterstate_vars%dqgdT_col                 , & ! Input:  [real(r8) (:)   ]  temperature derivative of "qg"                    
         qg                  =>   waterstate_vars%qg_col                    , & ! Input:  [real(r8) (:)   ]  specific humidity at ground surface (kg/kg)       
         h2osoi_ice          =>   waterstate_vars%h2osoi_ice_col            , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)                                
         h2osoi_liq          =>   waterstate_vars%h2osoi_liq_col            , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)                            
         h2osno              =>   waterstate_vars%h2osno_col                , & ! Input:  [real(r8) (:)   ]  snow water (mm H2O)                               
         qaf                 =>   waterstate_vars%qaf_lun                   , & ! Output: [real(r8) (:)   ]  urban canopy air specific humidity (kg/kg)        
         q_ref2m             =>   waterstate_vars%q_ref2m_patch             , & ! Output: [real(r8) (:)   ]  2 m height surface specific humidity (kg/kg)      
         rh_ref2m            =>   waterstate_vars%rh_ref2m_patch            , & ! Output: [real(r8) (:)   ]  2 m height surface relative humidity (%)          
         rh_ref2m_u          =>   waterstate_vars%rh_ref2m_u_patch          , & ! Output: [real(r8) (:)   ]  2 m height surface relative humidity (%)          

         forc_hgt_u_patch    =>   frictionvel_vars%forc_hgt_u_patch         , & ! Input:  [real(r8) (:)   ]  observational height of wind at pft-level (m)     
         forc_hgt_t_patch    =>   frictionvel_vars%forc_hgt_t_patch         , & ! Input:  [real(r8) (:)   ]  observational height of temperature at pft-level (m)
         ram1                =>   frictionvel_vars%ram1_patch               , & ! Output: [real(r8) (:)   ]  aerodynamical resistance (s/m)                    

         htvp                =>   energyflux_vars%htvp_col                  , & ! Input:  [real(r8) (:)   ]  latent heat of evaporation (/sublimation) (J/kg)  
         eflx_urban_ac       =>   energyflux_vars%eflx_urban_ac_col         , & ! Input:  [real(r8) (:)   ]  urban air conditioning flux (W/m**2)              
         eflx_urban_heat     =>   energyflux_vars%eflx_urban_heat_col       , & ! Input:  [real(r8) (:)   ]  urban heating flux (W/m**2)                       
         dlrad               =>   energyflux_vars%dlrad_patch               , & ! Output: [real(r8) (:)   ]  downward longwave radiation below the canopy (W/m**2)
         ulrad               =>   energyflux_vars%ulrad_patch               , & ! Output: [real(r8) (:)   ]  upward longwave radiation above the canopy (W/m**2)
         cgrnds              =>   energyflux_vars%cgrnds_patch              , & ! Output: [real(r8) (:)   ]  deriv, of soil sensible heat flux wrt soil temp (W/m**2/K)
         cgrndl              =>   energyflux_vars%cgrndl_patch              , & ! Output: [real(r8) (:)   ]  deriv of soil latent heat flux wrt soil temp (W/m**2/K)
         cgrnd               =>   energyflux_vars%cgrnd_patch               , & ! Output: [real(r8) (:)   ]  deriv. of soil energy flux wrt to soil temp (W/m**2/K)
         eflx_sh_grnd        =>   energyflux_vars%eflx_sh_grnd_patch        , & ! Output: [real(r8) (:)   ]  sensible heat flux from ground (W/m**2) [+ to atm]
         eflx_sh_tot         =>   energyflux_vars%eflx_sh_tot_patch         , & ! Output: [real(r8) (:)   ]  total sensible heat flux (W/m**2) [+ to atm]      
         eflx_sh_tot_u       =>   energyflux_vars%eflx_sh_tot_u_patch       , & ! Output: [real(r8) (:)   ]  urban total sensible heat flux (W/m**2) [+ to atm]
         eflx_sh_snow        =>   energyflux_vars%eflx_sh_snow_patch        , & ! Output: [real(r8) (:)   ]  sensible heat flux from snow (W/m**2) [+ to atm]  
         eflx_sh_soil        =>   energyflux_vars%eflx_sh_soil_patch        , & ! Output: [real(r8) (:)   ]  sensible heat flux from soil (W/m**2) [+ to atm]  
         eflx_sh_h2osfc      =>   energyflux_vars%eflx_sh_h2osfc_patch      , & ! Output: [real(r8) (:)   ]  sensible heat flux from soil (W/m**2) [+ to atm]  
         eflx_traffic        =>   energyflux_vars%eflx_traffic_lun          , & ! Output: [real(r8) (:)   ]  traffic sensible heat flux (W/m**2)               
         eflx_wasteheat      =>   energyflux_vars%eflx_wasteheat_lun        , & ! Output: [real(r8) (:)   ]  sensible heat flux from urban heating/cooling sources of waste heat (W/m**2)
         eflx_heat_from_ac   =>   energyflux_vars%eflx_heat_from_ac_lun     , & ! Output: [real(r8) (:)   ]  sensible heat flux put back into canyon due to removal by AC (W/m**2)
         taux                =>   energyflux_vars%taux_patch                , & ! Output: [real(r8) (:)   ]  wind (shear) stress: e-w (kg/m/s**2)              
         tauy                =>   energyflux_vars%tauy_patch                , & ! Output: [real(r8) (:)   ]  wind (shear) stress: n-s (kg/m/s**2)               

         qflx_evap_soi       =>   waterflux_vars%qflx_evap_soi_patch        , & ! Output: [real(r8) (:)   ]  soil evaporation (mm H2O/s) (+ = to atm)          
         qflx_tran_veg       =>   waterflux_vars%qflx_tran_veg_patch        , & ! Output: [real(r8) (:)   ]  vegetation transpiration (mm H2O/s) (+ = to atm)  
         qflx_evap_veg       =>   waterflux_vars%qflx_evap_veg_patch        , & ! Output: [real(r8) (:)   ]  vegetation evaporation (mm H2O/s) (+ = to atm)    
         qflx_evap_tot       =>   waterflux_vars%qflx_evap_tot_patch        , & ! Output: [real(r8) (:)   ]  qflx_evap_soi + qflx_evap_can + qflx_tran_veg     

         begl                =>   bounds%begl                               , &
         endl                =>   bounds%endl                                 &
         )

      ! Define fields that appear on the restart file for non-urban landunits 
      
      do fl = 1,num_nourbanl
         l = filter_nourbanl(fl)
         taf(l) = spval
         qaf(l) = spval
      end do

      ! Get time step
      nstep = get_nstep()

      ! Set constants (same as in Biogeophysics1Mod)
      beta(begl:endl) = 1._r8             ! Should be set to the same values as in Biogeophysics1Mod
      zii(begl:endl)  = 1000._r8          ! Should be set to the same values as in Biogeophysics1Mod

      ! Get current date
      dtime = get_step_size()
      call get_curr_date (year, month, day, secs)

      ! Compute canyontop wind using Masson (2000)

      do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         g = lun%gridcell(l)

         local_secp1(l)        = secs + nint((grc%londeg(g)/degpsec)/dtime)*dtime
         local_secp1(l)        = mod(local_secp1(l),isecspday)

         ! Error checks

         if (ht_roof(l) - z_d_town(l) <= z_0_town(l)) then
            write (iulog,*) 'aerodynamic parameter error in UrbanFluxes'
            write (iulog,*) 'h_r - z_d <= z_0'
            write (iulog,*) 'ht_roof, z_d_town, z_0_town: ', ht_roof(l), z_d_town(l), &
                 z_0_town(l)
            write (iulog,*) 'clm model is stopping'
            call endrun(decomp_index=l, clmlevel=namel, msg=errmsg(__FILE__, __LINE__))
         end if
         if (forc_hgt_u_patch(lun%pfti(l)) - z_d_town(l) <= z_0_town(l)) then
            write (iulog,*) 'aerodynamic parameter error in UrbanFluxes'
            write (iulog,*) 'h_u - z_d <= z_0'
            write (iulog,*) 'forc_hgt_u_patch, z_d_town, z_0_town: ', forc_hgt_u_patch(lun%pfti(l)), z_d_town(l), &
                 z_0_town(l)
            write (iulog,*) 'clm model is stopping'
            call endrun(decomp_index=l, clmlevel=namel, msg=errmsg(__FILE__, __LINE__))
         end if

         ! Magnitude of atmospheric wind

         ur(l) = max(1.0_r8,sqrt(forc_u(g)*forc_u(g)+forc_v(g)*forc_v(g)))

         ! Canyon top wind

         canyontop_wind(l) = ur(l) * &
              log( (ht_roof(l)-z_d_town(l)) / z_0_town(l) ) / &
              log( (forc_hgt_u_patch(lun%pfti(l))-z_d_town(l)) / z_0_town(l) )

         ! U component of canyon wind 

         if (canyon_hwr(l) < 0.5_r8) then  ! isolated roughness flow
            canyon_u_wind(l) = canyontop_wind(l) * exp( -0.5_r8*canyon_hwr(l)* &
                 (1._r8-(wind_hgt_canyon(l)/ht_roof(l))) )
         else if (canyon_hwr(l) < 1.0_r8) then ! wake interference flow
            canyon_u_wind(l) = canyontop_wind(l) * (1._r8+2._r8*(2._r8/rpi - 1._r8)* &
                 (ht_roof(l)/(ht_roof(l)/canyon_hwr(l)) - 0.5_r8)) * &
                 exp(-0.5_r8*canyon_hwr(l)*(1._r8-(wind_hgt_canyon(l)/ht_roof(l))))
         else  ! skimming flow
            canyon_u_wind(l) = canyontop_wind(l) * (2._r8/rpi) * &
                 exp(-0.5_r8*canyon_hwr(l)*(1._r8-(wind_hgt_canyon(l)/ht_roof(l))))
         end if

      end do

      ! Compute fluxes - Follows CLM approach for bare soils (Oleson et al 2004)

      do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         g = lun%gridcell(l)

         thm_g(l) = forc_t(g) + lapse_rate*forc_hgt_t_patch(lun%pfti(l))
         thv_g(l) = forc_th(g)*(1._r8+0.61_r8*forc_q(g))
         dth(l)   = thm_g(l)-taf(l)
         dqh(l)   = forc_q(g)-qaf(l)
         dthv     = dth(l)*(1._r8+0.61_r8*forc_q(g))+0.61_r8*forc_th(g)*dqh(l)
         zldis(l) = forc_hgt_u_patch(lun%pfti(l)) - z_d_town(l)

         ! Initialize Monin-Obukhov length and wind speed including convective velocity

         call MoninObukIni(ur(l), thv_g(l), dthv, zldis(l), z_0_town(l), um(l), obu(l))

      end do

      ! Initialize conductances
      wtus_roof(begl:endl)        = 0._r8
      wtus_road_perv(begl:endl)   = 0._r8
      wtus_road_imperv(begl:endl) = 0._r8
      wtus_sunwall(begl:endl)     = 0._r8
      wtus_shadewall(begl:endl)   = 0._r8
      wtuq_roof(begl:endl)        = 0._r8
      wtuq_road_perv(begl:endl)   = 0._r8
      wtuq_road_imperv(begl:endl) = 0._r8
      wtuq_sunwall(begl:endl)     = 0._r8
      wtuq_shadewall(begl:endl)   = 0._r8
      wtus_roof_unscl(begl:endl)        = 0._r8
      wtus_road_perv_unscl(begl:endl)   = 0._r8
      wtus_road_imperv_unscl(begl:endl) = 0._r8
      wtus_sunwall_unscl(begl:endl)     = 0._r8
      wtus_shadewall_unscl(begl:endl)   = 0._r8
      wtuq_roof_unscl(begl:endl)        = 0._r8
      wtuq_road_perv_unscl(begl:endl)   = 0._r8
      wtuq_road_imperv_unscl(begl:endl) = 0._r8
      wtuq_sunwall_unscl(begl:endl)     = 0._r8
      wtuq_shadewall_unscl(begl:endl)   = 0._r8

      ! Start stability iteration

      do iter = 1,niters

         ! Get friction velocity, relation for potential
         ! temperature and humidity profiles of surface boundary layer.

         if (num_urbanl > 0) then
            call FrictionVelocity(begl, endl, &
                 num_urbanl, filter_urbanl, &
                 z_d_town(begl:endl), z_0_town(begl:endl), z_0_town(begl:endl), z_0_town(begl:endl), &
                 obu(begl:endl), iter, ur(begl:endl), um(begl:endl), ustar(begl:endl), &
                 temp1(begl:endl), temp2(begl:endl), temp12m(begl:endl), temp22m(begl:endl), fm(begl:endl), &
                 frictionvel_vars, landunit_index=.true.)
         end if

         do fl = 1, num_urbanl
            l = filter_urbanl(fl)
            g = lun%gridcell(l)

            ! Determine aerodynamic resistance to fluxes from urban canopy air to
            ! atmosphere

            ramu(l) = 1._r8/(ustar(l)*ustar(l)/um(l))
            rahu(l) = 1._r8/(temp1(l)*ustar(l))
            rawu(l) = 1._r8/(temp2(l)*ustar(l))

            ! Determine magnitude of canyon wind by using horizontal wind determined
            ! previously and vertical wind from friction velocity (Masson 2000)

            canyon_wind(l) = sqrt(canyon_u_wind(l)**2._r8 + ustar(l)**2._r8)

            ! Determine canyon_resistance (currently this single resistance determines the
            ! resistance from urban surfaces (roof, pervious and impervious road, sunlit and
            ! shaded walls) to urban canopy air, since it is only dependent on wind speed
            ! Also from Masson 2000.

            canyon_resistance(l) = cpair * forc_rho(g) / (11.8_r8 + 4.2_r8*canyon_wind(l))

         end do

         ! This is the first term in the equation solutions for urban canopy air temperature
         ! and specific humidity (numerator) and is a landunit quantity
         do fl = 1, num_urbanl
            l = filter_urbanl(fl)
            g = lun%gridcell(l)

            taf_numer(l) = thm_g(l)/rahu(l)
            taf_denom(l) = 1._r8/rahu(l)
            qaf_numer(l) = forc_q(g)/rawu(l)
            qaf_denom(l) = 1._r8/rawu(l)

            ! First term needed for derivative of heat fluxes
            wtas(l) = 1._r8/rahu(l)
            wtaq(l) = 1._r8/rawu(l)

         end do


         ! Gather other terms for other urban columns for numerator and denominator of
         ! equations for urban canopy air temperature and specific humidity

         do fc = 1,num_urbanc
            c = filter_urbanc(fc)
            l = col%landunit(c)

            if (ctype(c) == icol_roof) then

               ! scaled sensible heat conductance
               wtus(c) = wtlunit_roof(l)/canyon_resistance(l)
               wtus_roof(l) = wtus(c)
               ! unscaled sensible heat conductance
               wtus_roof_unscl(l) = 1._r8/canyon_resistance(l)

               if (snow_depth(c) > 0._r8) then
                  fwet_roof = min(snow_depth(c)/0.05_r8, 1._r8)
               else
                  fwet_roof = (max(0._r8, h2osoi_liq(c,1)+h2osoi_ice(c,1))/pondmx_urban)**0.666666666666_r8
                  fwet_roof = min(fwet_roof,1._r8)
               end if
               if (qaf(l) > qg(c)) then 
                  fwet_roof = 1._r8
               end if
               ! scaled latent heat conductance
               wtuq(c) = fwet_roof*(wtlunit_roof(l)/canyon_resistance(l))
               wtuq_roof(l) = wtuq(c)
               ! unscaled latent heat conductance
               wtuq_roof_unscl(l) = fwet_roof*(1._r8/canyon_resistance(l))

               ! wasteheat from heating/cooling
               if (trim(urban_hac) == urban_wasteheat_on) then
                  eflx_wasteheat_roof(l) = ac_wasteheat_factor * eflx_urban_ac(c) + &
                       ht_wasteheat_factor * eflx_urban_heat(c)
               else
                  eflx_wasteheat_roof(l) = 0._r8
               end if

               ! If air conditioning on, always replace heat removed with heat into canyon
               if (trim(urban_hac) == urban_hac_on .or. trim(urban_hac) == urban_wasteheat_on) then
                  eflx_heat_from_ac_roof(l) = abs(eflx_urban_ac(c))
               else
                  eflx_heat_from_ac_roof(l) = 0._r8
               end if

            else if (ctype(c) == icol_road_perv) then

               ! scaled sensible heat conductance
               wtus(c) = wtroad_perv(l)*(1._r8-wtlunit_roof(l))/canyon_resistance(l)
               wtus_road_perv(l) = wtus(c)
               ! unscaled sensible heat conductance
               wtus_road_perv_unscl(l) = 1._r8/canyon_resistance(l)

               ! scaled latent heat conductance
               wtuq(c) = wtroad_perv(l)*(1._r8-wtlunit_roof(l))/canyon_resistance(l)
               wtuq_road_perv(l) = wtuq(c)
               ! unscaled latent heat conductance
               wtuq_road_perv_unscl(l) = 1._r8/canyon_resistance(l)

               if (use_vsfm) then
                  if (qaf(l) < qg(c)) then
                     if (do_soilevap_beta()) then
                        wtuq_road_perv(l)       = soilbeta(c)*wtuq_road_perv(l)
                        wtuq_road_perv_unscl(l) = soilbeta(c)*wtuq_road_perv_unscl(l)
                     endif
                  endif
               endif
            else if (ctype(c) == icol_road_imperv) then

               ! scaled sensible heat conductance
               wtus(c) = (1._r8-wtroad_perv(l))*(1._r8-wtlunit_roof(l))/canyon_resistance(l)
               wtus_road_imperv(l) = wtus(c)
               ! unscaled sensible heat conductance
               wtus_road_imperv_unscl(l) = 1._r8/canyon_resistance(l)

               if (snow_depth(c) > 0._r8) then
                  fwet_road_imperv = min(snow_depth(c)/0.05_r8, 1._r8)
               else
                  fwet_road_imperv = (max(0._r8, h2osoi_liq(c,1)+h2osoi_ice(c,1))/pondmx_urban)**0.666666666666_r8
                  fwet_road_imperv = min(fwet_road_imperv,1._r8)
               end if
               if (qaf(l) > qg(c)) then 
                  fwet_road_imperv = 1._r8
               end if
               ! scaled latent heat conductance
               wtuq(c) = fwet_road_imperv*(1._r8-wtroad_perv(l))*(1._r8-wtlunit_roof(l))/canyon_resistance(l)
               wtuq_road_imperv(l) = wtuq(c)
               ! unscaled latent heat conductance
               wtuq_road_imperv_unscl(l) = fwet_road_imperv*(1._r8/canyon_resistance(l))

            else if (ctype(c) == icol_sunwall) then

               ! scaled sensible heat conductance
               wtus(c) = canyon_hwr(l)*(1._r8-wtlunit_roof(l))/canyon_resistance(l)
               wtus_sunwall(l) = wtus(c)
               ! unscaled sensible heat conductance
               wtus_sunwall_unscl(l) = 1._r8/canyon_resistance(l)

               ! scaled latent heat conductance
               wtuq(c) = 0._r8
               wtuq_sunwall(l) = wtuq(c)
               ! unscaled latent heat conductance
               wtuq_sunwall_unscl(l) = 0._r8

               ! wasteheat from heating/cooling
               if (trim(urban_hac) == urban_wasteheat_on) then
                  eflx_wasteheat_sunwall(l) = ac_wasteheat_factor * eflx_urban_ac(c) + &
                       ht_wasteheat_factor * eflx_urban_heat(c)
               else
                  eflx_wasteheat_sunwall(l) = 0._r8
               end if

               ! If air conditioning on, always replace heat removed with heat into canyon
               if (trim(urban_hac) == urban_hac_on .or. trim(urban_hac) == urban_wasteheat_on) then
                  eflx_heat_from_ac_sunwall(l) = abs(eflx_urban_ac(c))
               else
                  eflx_heat_from_ac_sunwall(l) = 0._r8
               end if

            else if (ctype(c) == icol_shadewall) then

               ! scaled sensible heat conductance
               wtus(c) = canyon_hwr(l)*(1._r8-wtlunit_roof(l))/canyon_resistance(l)
               wtus_shadewall(l) = wtus(c)
               ! unscaled sensible heat conductance
               wtus_shadewall_unscl(l) = 1._r8/canyon_resistance(l)

               ! scaled latent heat conductance
               wtuq(c) = 0._r8
               wtuq_shadewall(l) = wtuq(c)
               ! unscaled latent heat conductance
               wtuq_shadewall_unscl(l) = 0._r8

               ! wasteheat from heating/cooling
               if (trim(urban_hac) == urban_wasteheat_on) then
                  eflx_wasteheat_shadewall(l) = ac_wasteheat_factor * eflx_urban_ac(c) + &
                       ht_wasteheat_factor * eflx_urban_heat(c)
               else
                  eflx_wasteheat_shadewall(l) = 0._r8
               end if

               ! If air conditioning on, always replace heat removed with heat into canyon
               if (trim(urban_hac) == urban_hac_on .or. trim(urban_hac) == urban_wasteheat_on) then
                  eflx_heat_from_ac_shadewall(l) = abs(eflx_urban_ac(c))
               else
                  eflx_heat_from_ac_shadewall(l) = 0._r8
               end if
            else
               write(iulog,*) 'c, ctype, pi = ', c, ctype(c), pi
               write(iulog,*) 'Column indices for: shadewall, sunwall, road_imperv, road_perv, roof: '
               write(iulog,*) icol_shadewall, icol_sunwall, icol_road_imperv, icol_road_perv, icol_roof
               call endrun(decomp_index=l, clmlevel=namel, msg="ERROR, ctype out of range"//errmsg(__FILE__, __LINE__))
            end if

            taf_numer(l) = taf_numer(l) + t_grnd(c)*wtus(c)
            taf_denom(l) = taf_denom(l) + wtus(c)
            qaf_numer(l) = qaf_numer(l) + qg(c)*wtuq(c)
            qaf_denom(l) = qaf_denom(l) + wtuq(c)

         end do

         ! Calculate new urban canopy air temperature and specific humidity

         do fl = 1, num_urbanl
            l = filter_urbanl(fl)
            g = lun%gridcell(l)

            ! Total waste heat and heat from AC is sum of heat for walls and roofs
            ! accounting for different surface areas
            eflx_wasteheat(l) = wtlunit_roof(l)*eflx_wasteheat_roof(l) + &
                 (1._r8-wtlunit_roof(l))*(canyon_hwr(l)*(eflx_wasteheat_sunwall(l) + &
                 eflx_wasteheat_shadewall(l)))

            ! Limit wasteheat to ensure that we don't get any unrealistically strong 
            ! positive feedbacks due to AC in a warmer climate
            eflx_wasteheat(l) = min(eflx_wasteheat(l),wasteheat_limit)

            eflx_heat_from_ac(l) = wtlunit_roof(l)*eflx_heat_from_ac_roof(l) + &
                 (1._r8-wtlunit_roof(l))*(canyon_hwr(l)*(eflx_heat_from_ac_sunwall(l) + &
                 eflx_heat_from_ac_shadewall(l)))

            ! Calculate traffic heat flux
            ! Only comes from impervious road
            eflx_traffic(l) = (1._r8-wtlunit_roof(l))*(1._r8-wtroad_perv(l))* &
                 eflx_traffic_factor(l)

            taf(l) = taf_numer(l)/taf_denom(l)
            qaf(l) = qaf_numer(l)/qaf_denom(l)

            wts_sum(l) = wtas(l) + wtus_roof(l) + wtus_road_perv(l) + &
                 wtus_road_imperv(l) + wtus_sunwall(l) + wtus_shadewall(l)

            wtq_sum(l) = wtaq(l) + wtuq_roof(l) + wtuq_road_perv(l) + &
                 wtuq_road_imperv(l) + wtuq_sunwall(l) + wtuq_shadewall(l)

         end do

         ! This section of code is not required if niters = 1
         ! Determine stability using new taf and qaf
         ! TODO: Some of these constants replicate what is in FrictionVelocity and BareGround fluxes should consildate. EBK
         do fl = 1, num_urbanl
            l = filter_urbanl(fl)
            g = lun%gridcell(l)

            dth(l) = thm_g(l)-taf(l)
            dqh(l) = forc_q(g)-qaf(l)
            tstar = temp1(l)*dth(l)
            qstar = temp2(l)*dqh(l)
            thvstar = tstar*(1._r8+0.61_r8*forc_q(g)) + 0.61_r8*forc_th(g)*qstar
            zeta = zldis(l)*vkc*grav*thvstar/(ustar(l)**2*thv_g(l))

            if (zeta >= 0._r8) then                   !stable
               zeta = min(2._r8,max(zeta,0.01_r8))
               um(l) = max(ur(l),0.1_r8)
            else                                      !unstable
               zeta = max(-100._r8,min(zeta,-0.01_r8))
               wc = beta(l)*(-grav*ustar(l)*thvstar*zii(l)/thv_g(l))**0.333_r8
               um(l) = sqrt(ur(l)*ur(l) + wc*wc)
            end if

            obu(l) = zldis(l)/zeta
         end do

      end do   ! end iteration

      ! Determine fluxes from canyon surfaces

      ! the following initializations are needed to ensure that the values are 0 over non-
      ! active urban Patches
      eflx_sh_grnd_scale(bounds%begp : bounds%endp) = 0._r8
      qflx_evap_soi_scale(bounds%begp : bounds%endp) = 0._r8

      do f = 1, num_urbanp

         p = filter_urbanp(f)
         c = pft%column(p)
         g = pft%gridcell(p)
         l = pft%landunit(p)

         ram1(p) = ramu(l)  !pass value to global variable

         ! Upward and downward canopy longwave are zero

         ulrad(p)  = 0._r8
         dlrad(p)  = 0._r8

         ! Derivative of sensible and latent heat fluxes with respect to 
         ! ground temperature

         if (ctype(c) == icol_roof) then
            cgrnds(p) = forc_rho(g) * cpair * (wtas(l) + wtus_road_perv(l) +  &
                 wtus_road_imperv(l) + wtus_sunwall(l) + wtus_shadewall(l)) * &
                 (wtus_roof_unscl(l)/wts_sum(l))
            cgrndl(p) = forc_rho(g) * (wtaq(l) + wtuq_road_perv(l) +  &
                 wtuq_road_imperv(l) + wtuq_sunwall(l) + wtuq_shadewall(l)) * &
                 (wtuq_roof_unscl(l)/wtq_sum(l))*dqgdT(c)
         else if (ctype(c) == icol_road_perv) then
            cgrnds(p) = forc_rho(g) * cpair * (wtas(l) + wtus_roof(l) +  &
                 wtus_road_imperv(l) + wtus_sunwall(l) + wtus_shadewall(l)) * &
                 (wtus_road_perv_unscl(l)/wts_sum(l))
            cgrndl(p) = forc_rho(g) * (wtaq(l) + wtuq_roof(l) +  &
                 wtuq_road_imperv(l) + wtuq_sunwall(l) + wtuq_shadewall(l)) * &
                 (wtuq_road_perv_unscl(l)/wtq_sum(l))*dqgdT(c)
         else if (ctype(c) == icol_road_imperv) then
            cgrnds(p) = forc_rho(g) * cpair * (wtas(l) + wtus_roof(l) +  &
                 wtus_road_perv(l) + wtus_sunwall(l) + wtus_shadewall(l)) * &
                 (wtus_road_imperv_unscl(l)/wts_sum(l))
            cgrndl(p) = forc_rho(g) * (wtaq(l) + wtuq_roof(l) +  &
                 wtuq_road_perv(l) + wtuq_sunwall(l) + wtuq_shadewall(l)) * &
                 (wtuq_road_imperv_unscl(l)/wtq_sum(l))*dqgdT(c)
         else if (ctype(c) == icol_sunwall) then
            cgrnds(p) = forc_rho(g) * cpair * (wtas(l) + wtus_roof(l) +  &
                 wtus_road_perv(l) + wtus_road_imperv(l) + wtus_shadewall(l)) * &
                 (wtus_sunwall_unscl(l)/wts_sum(l))
            cgrndl(p) = 0._r8
         else if (ctype(c) == icol_shadewall) then
            cgrnds(p) = forc_rho(g) * cpair * (wtas(l) + wtus_roof(l) +  &
                 wtus_road_perv(l) + wtus_road_imperv(l) + wtus_sunwall(l)) * &
                 (wtus_shadewall_unscl(l)/wts_sum(l))
            cgrndl(p) = 0._r8
         end if
         cgrnd(p)  = cgrnds(p) + cgrndl(p)*htvp(c)

         ! Surface fluxes of momentum, sensible and latent heat

         taux(p)          = -forc_rho(g)*forc_u(g)/ramu(l)
         tauy(p)          = -forc_rho(g)*forc_v(g)/ramu(l)

         ! Use new canopy air temperature
         dth(l) = taf(l) - t_grnd(c)

         if (ctype(c) == icol_roof) then
            eflx_sh_grnd(p)  = -forc_rho(g)*cpair*wtus_roof_unscl(l)*dth(l)
            eflx_sh_snow(p)  = 0._r8
            eflx_sh_soil(p)  = 0._r8
            eflx_sh_h2osfc(p)= 0._r8
         else if (ctype(c) == icol_road_perv) then
            eflx_sh_grnd(p)  = -forc_rho(g)*cpair*wtus_road_perv_unscl(l)*dth(l)
            eflx_sh_snow(p)  = 0._r8
            eflx_sh_soil(p)  = 0._r8
            eflx_sh_h2osfc(p)= 0._r8
         else if (ctype(c) == icol_road_imperv) then
            eflx_sh_grnd(p)  = -forc_rho(g)*cpair*wtus_road_imperv_unscl(l)*dth(l)
            eflx_sh_snow(p)  = 0._r8
            eflx_sh_soil(p)  = 0._r8
            eflx_sh_h2osfc(p)= 0._r8
         else if (ctype(c) == icol_sunwall) then
            eflx_sh_grnd(p)  = -forc_rho(g)*cpair*wtus_sunwall_unscl(l)*dth(l)
            eflx_sh_snow(p)  = 0._r8
            eflx_sh_soil(p)  = 0._r8
            eflx_sh_h2osfc(p)= 0._r8
         else if (ctype(c) == icol_shadewall) then
            eflx_sh_grnd(p)  = -forc_rho(g)*cpair*wtus_shadewall_unscl(l)*dth(l)
            eflx_sh_snow(p)  = 0._r8
            eflx_sh_soil(p)  = 0._r8
            eflx_sh_h2osfc(p)= 0._r8
         end if

         eflx_sh_tot(p)   = eflx_sh_grnd(p)
         eflx_sh_tot_u(p) = eflx_sh_tot(p)

         dqh(l) = qaf(l) - qg(c)

         if (ctype(c) == icol_roof) then
            qflx_evap_soi(p) = -forc_rho(g)*wtuq_roof_unscl(l)*dqh(l)
         else if (ctype(c) == icol_road_perv) then
            ! Evaporation assigned to soil term if dew or snow
            ! or if no liquid water available in soil column
            if (dqh(l) > 0._r8 .or. frac_sno(c) > 0._r8 .or. soilalpha_u(c) <= 0._r8) then
               qflx_evap_soi(p) = -forc_rho(g)*wtuq_road_perv_unscl(l)*dqh(l)
               qflx_tran_veg(p) = 0._r8
               ! Otherwise, evaporation assigned to transpiration term
            else
               qflx_evap_soi(p) = 0._r8
               qflx_tran_veg(p) = -forc_rho(g)*wtuq_road_perv_unscl(l)*dqh(l)
            end if
            qflx_evap_veg(p) = qflx_tran_veg(p)
         else if (ctype(c) == icol_road_imperv) then
            qflx_evap_soi(p) = -forc_rho(g)*wtuq_road_imperv_unscl(l)*dqh(l)
         else if (ctype(c) == icol_sunwall) then
            qflx_evap_soi(p) = 0._r8
         else if (ctype(c) == icol_shadewall) then
            qflx_evap_soi(p) = 0._r8
         end if

         ! SCALED sensible and latent heat flux for error check
         eflx_sh_grnd_scale(p)  = -forc_rho(g)*cpair*wtus(c)*dth(l)
         qflx_evap_soi_scale(p) = -forc_rho(g)*wtuq(c)*dqh(l)

      end do

      ! Check to see that total sensible and latent heat equal the sum of
      ! the scaled heat fluxes above
      do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         g = lun%gridcell(l)
         eflx(l)       = -(forc_rho(g)*cpair/rahu(l))*(thm_g(l) - taf(l))
         qflx(l)       = -(forc_rho(g)/rawu(l))*(forc_q(g) - qaf(l))
         eflx_scale(l) = sum(eflx_sh_grnd_scale(lun%pfti(l):lun%pftf(l)))
         qflx_scale(l) = sum(qflx_evap_soi_scale(lun%pfti(l):lun%pftf(l)))
         eflx_err(l)   = eflx_scale(l) - eflx(l)
         qflx_err(l)   = qflx_scale(l) - qflx(l)
      end do

      found = .false.
      do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         if (abs(eflx_err(l)) > 0.01_r8) then
            found = .true.
            indexl = l
            exit
         end if
      end do
      if ( found ) then
         write(iulog,*)'WARNING:  Total sensible heat does not equal sum of scaled heat fluxes for urban columns ',&
              ' nstep = ',nstep,' indexl= ',indexl,' eflx_err= ',eflx_err(indexl)
         if (abs(eflx_err(indexl)) > .01_r8) then
            write(iulog,*)'clm model is stopping - error is greater than .01 W/m**2'
            write(iulog,*)'eflx_scale    = ',eflx_scale(indexl)
            write(iulog,*)'eflx_sh_grnd_scale: ',eflx_sh_grnd_scale(lun%pfti(indexl):lun%pftf(indexl))
            write(iulog,*)'eflx          = ',eflx(indexl)
            call endrun(decomp_index=indexl, clmlevel=namel, msg=errmsg(__FILE__, __LINE__))
         end if
      end if

      found = .false.
      do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         ! 4.e-9 kg/m**2/s = 0.01 W/m**2
         if (abs(qflx_err(l)) > 4.e-9_r8) then
            found = .true.
            indexl = l
            exit
         end if
      end do
      if ( found ) then
         write(iulog,*)'WARNING:  Total water vapor flux does not equal sum of scaled water vapor fluxes for urban columns ',&
              ' nstep = ',nstep,' indexl= ',indexl,' qflx_err= ',qflx_err(indexl)
         if (abs(qflx_err(indexl)) > 4.e-9_r8) then
            write(iulog,*)'clm model is stopping - error is greater than 4.e-9 kg/m**2/s'
            write(iulog,*)'qflx_scale    = ',qflx_scale(indexl)
            write(iulog,*)'qflx          = ',qflx(indexl)
            call endrun(decomp_index=indexl, clmlevel=namel, msg=errmsg(__FILE__, __LINE__))
         end if
      end if

      ! Gather terms required to determine internal building temperature

      do fc = 1,num_urbanc
         c = filter_urbanc(fc)
         l = col%landunit(c)

         if (ctype(c) == icol_roof) then
            t_roof_innerl(l) = t_soisno(c,nlevurb)
         else if (ctype(c) == icol_sunwall) then
            t_sunwall_innerl(l) = t_soisno(c,nlevurb)
         else if (ctype(c) == icol_shadewall) then
            t_shadewall_innerl(l) = t_soisno(c,nlevurb)
         end if

      end do

      ! Calculate internal building temperature
      do fl = 1, num_urbanl
         l = filter_urbanl(fl)

         lngth_roof = (ht_roof(l)/canyon_hwr(l))*wtlunit_roof(l)/(1._r8-wtlunit_roof(l))
         t_building(l) = (ht_roof(l)*(t_shadewall_innerl(l) + t_sunwall_innerl(l)) &
              +lngth_roof*t_roof_innerl(l))/(2._r8*ht_roof(l)+lngth_roof)
      end do

      ! No roots for urban except for pervious road

      do j = 1, nlevgrnd
         do f = 1, num_urbanp
            p = filter_urbanp(f)
            c = pft%column(p)
            if (ctype(c) == icol_road_perv) then
               rootr(p,j) = rootr_road_perv(c,j)
            else
               rootr(p,j) = 0._r8
            end if
         end do
      end do

      do f = 1, num_urbanp

         p = filter_urbanp(f)
         c = pft%column(p)
         g = pft%gridcell(p)
         l = pft%landunit(p)

         ! Use urban canopy air temperature and specific humidity to represent 
         ! 2-m temperature and humidity

         t_ref2m(p) = taf(l)
         q_ref2m(p) = qaf(l)
         t_ref2m_u(p) = taf(l)

         ! 2 m height relative humidity

         call QSat(t_ref2m(p), forc_pbot(g), e_ref2m, de2mdT, qsat_ref2m, dqsat2mdT)
         rh_ref2m(p) = min(100._r8, q_ref2m(p) / qsat_ref2m * 100._r8)
         rh_ref2m_u(p) = rh_ref2m(p)

         ! Variables needed by history tape

         t_veg(p) = forc_t(g)

      end do

    end associate

  end subroutine UrbanFluxes

end module UrbanFluxesMod
