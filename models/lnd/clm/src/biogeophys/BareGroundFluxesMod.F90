module BareGroundFluxesMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Compute sensible and latent fluxes and their derivatives with respect
  ! to ground temperature using ground temperatures from previous time step.
  !
  ! !USES:
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use decompMod            , only : bounds_type
  use ch4Mod               , only : ch4_type
  use atm2lndType          , only : atm2lnd_type
  use EnergyFluxType       , only : energyflux_type
  use FrictionVelocityMod  , only : frictionvel_type
  use SoilStateType        , only : soilstate_type
  use TemperatureType      , only : temperature_type
  use PhotosynthesisMod    , only : photosyns_type
  use WaterfluxType        , only : waterflux_type
  use WaterstateType       , only : waterstate_type
  use HumanIndexMod        , only : humanindex_type
  use LandunitType         , only : lun                
  use ColumnType           , only : col                
  use PatchType            , only : patch                
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: BareGroundFluxes   ! Calculate sensible and latent heat fluxes
  !------------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  subroutine BareGroundFluxes(bounds, num_noexposedvegp, filter_noexposedvegp, &
       atm2lnd_inst, soilstate_inst, &
       frictionvel_inst, ch4_inst, energyflux_inst, temperature_inst, &
       waterflux_inst, waterstate_inst, photosyns_inst, humanindex_inst)
    !
    ! !DESCRIPTION:
    ! Compute sensible and latent fluxes and their derivatives with respect
    ! to ground temperature using ground temperatures from previous time step.
    !
    ! !USES:
    use shr_const_mod        , only : SHR_CONST_RGAS
    use clm_varpar           , only : nlevgrnd
    use clm_varcon           , only : cpair, vkc, grav, denice, denh2o
    use clm_varctl           , only : use_lch4
    use landunit_varcon      , only : istsoil, istcrop
    use FrictionVelocityMod  , only : FrictionVelocity, MoninObukIni
    use QSatMod              , only : QSat
    use SurfaceResistanceMod , only : do_soilevap_beta
    use HumanIndexMod        , only : calc_human_stress_indices, Wet_Bulb, Wet_BulbS, HeatIndex, AppTemp, &
                                      swbgt, hmdex, dis_coi, dis_coiS, THIndex, &
                                      SwampCoolEff, KtoC, VaporPres
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds  
    integer                , intent(in)    :: num_noexposedvegp       ! number of points in filter_noexposedvegp
    integer                , intent(in)    :: filter_noexposedvegp(:) ! patch filter where frac_veg_nosno is 0 
                                                                      ! (but does NOT include lake or urban)
    type(atm2lnd_type)     , intent(in)    :: atm2lnd_inst
    type(soilstate_type)   , intent(inout) :: soilstate_inst
    type(frictionvel_type) , intent(inout) :: frictionvel_inst
    type(ch4_type)         , intent(inout) :: ch4_inst
    type(energyflux_type)  , intent(inout) :: energyflux_inst
    type(temperature_type) , intent(inout) :: temperature_inst
    type(waterflux_type)   , intent(inout) :: waterflux_inst
    type(waterstate_type)  , intent(inout) :: waterstate_inst
    type(photosyns_type)   , intent(inout) :: photosyns_inst
    type(humanindex_type)  , intent(inout) :: humanindex_inst
    !
    ! !LOCAL VARIABLES:
    integer, parameter  :: niters = 3            ! maximum number of iterations for surface temperature
    integer  :: p,c,g,f,j,l                      ! indices
    integer  :: iter                             ! iteration index
    real(r8) :: zldis(bounds%begp:bounds%endp)   ! reference height "minus" zero displacement height [m]
    real(r8) :: displa(bounds%begp:bounds%endp)  ! displacement height [m]
    real(r8) :: zeta                             ! dimensionless height used in Monin-Obukhov theory
    real(r8) :: wc                               ! convective velocity [m/s]
    real(r8) :: dth(bounds%begp:bounds%endp)     ! diff of virtual temp. between ref. height and surface
    real(r8) :: dthv                             ! diff of vir. poten. temp. between ref. height and surface
    real(r8) :: dqh(bounds%begp:bounds%endp)     ! diff of humidity between ref. height and surface
    real(r8) :: obu(bounds%begp:bounds%endp)     ! Monin-Obukhov length (m)
    real(r8) :: ur(bounds%begp:bounds%endp)      ! wind speed at reference height [m/s]
    real(r8) :: um(bounds%begp:bounds%endp)      ! wind speed including the stablity effect [m/s]
    real(r8) :: temp1(bounds%begp:bounds%endp)   ! relation for potential temperature profile
    real(r8) :: temp12m(bounds%begp:bounds%endp) ! relation for potential temperature profile applied at 2-m
    real(r8) :: temp2(bounds%begp:bounds%endp)   ! relation for specific humidity profile
    real(r8) :: temp22m(bounds%begp:bounds%endp) ! relation for specific humidity profile applied at 2-m
    real(r8) :: ustar(bounds%begp:bounds%endp)   ! friction velocity [m/s]
    real(r8) :: tstar                            ! temperature scaling parameter
    real(r8) :: qstar                            ! moisture scaling parameter
    real(r8) :: thvstar                          ! virtual potential temperature scaling parameter
    real(r8) :: cf_bare                          ! heat transfer coefficient from bare ground [-]
    real(r8) :: ram                              ! aerodynamical resistance [s/m]
    real(r8) :: rah                              ! thermal resistance [s/m]
    real(r8) :: raw                              ! moisture resistance [s/m]
    real(r8) :: raih                             ! temporary variable [kg/m2/s]
    real(r8) :: raiw                             ! temporary variable [kg/m2/s]
    real(r8) :: fm(bounds%begp:bounds%endp)      ! needed for BGC only to diagnose 10m wind speed
    real(r8) :: z0mg_patch(bounds%begp:bounds%endp)
    real(r8) :: z0hg_patch(bounds%begp:bounds%endp)
    real(r8) :: z0qg_patch(bounds%begp:bounds%endp)
    real(r8) :: e_ref2m                          ! 2 m height surface saturated vapor pressure [Pa]
    real(r8) :: de2mdT                           ! derivative of 2 m height surface saturated vapor pressure on t_ref2m
    real(r8) :: qsat_ref2m                       ! 2 m height surface saturated specific humidity [kg/kg]
    real(r8) :: dqsat2mdT                        ! derivative of 2 m height surface saturated specific humidity on t_ref2m 
    real(r8) :: www                              ! surface soil wetness [-]
    !------------------------------------------------------------------------------

    associate(                                                                       & 
         snl                    => col%snl                                      , & ! Input:  [integer  (:)   ]  number of snow layers                                                  
         dz                     => col%dz                                       , & ! Input:  [real(r8) (:,:) ]  layer depth (m)                                                     
         zii                    => col%zii                                      , & ! Input:  [real(r8) (:)   ]  convective boundary height [m]                                        

         tc_ref2m               => humanindex_inst%tc_ref2m_patch               , & ! Output: [real(r8) (:)   ]  2 m height surface air temperature (C)
         vap_ref2m              => humanindex_inst%vap_ref2m_patch              , & ! Output: [real(r8) (:)   ]  2 m height vapor pressure (Pa)
         appar_temp_ref2m       => humanindex_inst%appar_temp_ref2m_patch       , & ! Output: [real(r8) (:)   ]  2 m apparent temperature (C)
         appar_temp_ref2m_r     => humanindex_inst%appar_temp_ref2m_r_patch     , & ! Output: [real(r8) (:)   ]  Rural 2 m apparent temperature (C)
         swbgt_ref2m            => humanindex_inst%swbgt_ref2m_patch            , & ! Output: [real(r8) (:)   ]  2 m Simplified Wetbulb Globe temperature (C)
         swbgt_ref2m_r          => humanindex_inst%swbgt_ref2m_r_patch          , & ! Output: [real(r8) (:)   ]  Rural 2 m Simplified Wetbulb Globe temperature (C)
         humidex_ref2m          => humanindex_inst%humidex_ref2m_patch          , & ! Output: [real(r8) (:)   ]  2 m Humidex (C)
         humidex_ref2m_r        => humanindex_inst%humidex_ref2m_r_patch        , & ! Output: [real(r8) (:)   ]  Rural 2 m Humidex (C)
         wbt_ref2m              => humanindex_inst%wbt_ref2m_patch              , & ! Output: [real(r8) (:)   ]  2 m Stull Wet Bulb temperature (C)
         wbt_ref2m_r            => humanindex_inst%wbt_ref2m_r_patch            , & ! Output: [real(r8) (:)   ]  Rural 2 m Stull Wet Bulb temperature (C)
         wb_ref2m               => humanindex_inst%wb_ref2m_patch               , & ! Output: [real(r8) (:)   ]  2 m Wet Bulb temperature (C)
         wb_ref2m_r             => humanindex_inst%wb_ref2m_r_patch             , & ! Output: [real(r8) (:)   ]  Rural 2 m Wet Bulb temperature (C)
         teq_ref2m              => humanindex_inst%teq_ref2m_patch              , & ! Output: [real(r8) (:)   ]  2 m height Equivalent temperature (K)
         teq_ref2m_r            => humanindex_inst%teq_ref2m_r_patch            , & ! Output: [real(r8) (:)   ]  Rural 2 m Equivalent temperature (K)
         ept_ref2m              => humanindex_inst%ept_ref2m_patch              , & ! Output: [real(r8) (:)   ]  2 m height Equivalent Potential temperature (K)
         ept_ref2m_r            => humanindex_inst%ept_ref2m_r_patch            , & ! Output: [real(r8) (:)   ]  Rural 2 m height Equivalent Potential temperature (K)
         discomf_index_ref2m    => humanindex_inst%discomf_index_ref2m_patch    , & ! Output: [real(r8) (:)   ]  2 m Discomfort Index temperature (C)
         discomf_index_ref2m_r  => humanindex_inst%discomf_index_ref2m_r_patch  , & ! Output: [real(r8) (:)   ]  Rural 2 m Discomfort Index temperature (C)
         discomf_index_ref2mS   => humanindex_inst%discomf_index_ref2mS_patch   , & ! Output: [real(r8) (:)   ]  2 m height Discomfort Index Stull temperature (C)
         discomf_index_ref2mS_r => humanindex_inst%discomf_index_ref2mS_r_patch , & ! Output: [real(r8) (:)   ]  Rural 2 m Discomfort Index Stull temperature (K)
         nws_hi_ref2m           => humanindex_inst%nws_hi_ref2m_patch           , & ! Output: [real(r8) (:)   ]  2 m NWS Heat Index (C)
         nws_hi_ref2m_r         => humanindex_inst%nws_hi_ref2m_r_patch         , & ! Output: [real(r8) (:)   ]  Rural 2 m NWS Heat Index (C)
         thip_ref2m             => humanindex_inst%thip_ref2m_patch             , & ! Output: [real(r8) (:)   ]  2 m Temperature Humidity Index Physiology (C)
         thip_ref2m_r           => humanindex_inst%thip_ref2m_r_patch           , & ! Output: [real(r8) (:)   ]  Rural 2 m Temperature Humidity Index Physiology (C)
         thic_ref2m             => humanindex_inst%thic_ref2m_patch             , & ! Output: [real(r8) (:)   ]  2 m Temperature Humidity Index Comfort (C)
         thic_ref2m_r           => humanindex_inst%thic_ref2m_r_patch           , & ! Output: [real(r8) (:)   ]  Rural 2 m Temperature Humidity Index Comfort (C)
         swmp65_ref2m           => humanindex_inst%swmp65_ref2m_patch           , & ! Output: [real(r8) (:)   ]  2 m Swamp Cooler temperature 65% effi (C)
         swmp65_ref2m_r         => humanindex_inst%swmp65_ref2m_r_patch         , & ! Output: [real(r8) (:)   ]  Rural 2 m Swamp Cooler temperature 65% effi (C)
         swmp80_ref2m           => humanindex_inst%swmp80_ref2m_patch           , & ! Output: [real(r8) (:)   ]  2 m Swamp Cooler temperature 80% effi (C)
         swmp80_ref2m_r         => humanindex_inst%swmp80_ref2m_r_patch         , & ! Output: [real(r8) (:)   ]  Rural 2 m Swamp Cooler temperature 80% effi (C)

         forc_u                 => atm2lnd_inst%forc_u_grc                      , & ! Input:  [real(r8) (:)   ]  atmospheric wind speed in east direction (m/s)                        
         forc_v                 => atm2lnd_inst%forc_v_grc                      , & ! Input:  [real(r8) (:)   ]  atmospheric wind speed in north direction (m/s)                       
         forc_th                => atm2lnd_inst%forc_th_downscaled_col          , & ! Input:  [real(r8) (:)   ]  atmospheric potential temperature (Kelvin)                            
         forc_t                 => atm2lnd_inst%forc_t_downscaled_col           , & ! Input:  [real(r8) (:)   ]  atmospheric temperature (Kelvin) 
         forc_pbot              => atm2lnd_inst%forc_pbot_downscaled_col        , & ! Input:  [real(r8) (:)   ]  atmospheric pressure (Pa)                                             
         forc_rho               => atm2lnd_inst%forc_rho_downscaled_col         , & ! Input:  [real(r8) (:)   ]  density (kg/m**3)                                                     
         forc_q                 => atm2lnd_inst%forc_q_downscaled_col           , & ! Input:  [real(r8) (:)   ]  atmospheric specific humidity (kg/kg)                                 

         watsat                 => soilstate_inst%watsat_col                    , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)                      
         soilbeta               => soilstate_inst%soilbeta_col                  , & ! Input:  [real(r8) (:)   ]  soil wetness relative to field capacity                               
         rootr                  => soilstate_inst%rootr_patch                   , & ! Output: [real(r8) (:,:) ]  effective fraction of roots in each soil layer                      
         t_soisno               => temperature_inst%t_soisno_col                , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)                                           
         t_grnd                 => temperature_inst%t_grnd_col                  , & ! Input:  [real(r8) (:)   ]  ground surface temperature [K]                                        
         thv                    => temperature_inst%thv_col                     , & ! Input:  [real(r8) (:)   ]  virtual potential temperature (kelvin)                                
         thm                    => temperature_inst%thm_patch                   , & ! Input:  [real(r8) (:)   ]  intermediate variable (forc_t+0.0098*forc_hgt_t_patch)                  
         t_h2osfc               => temperature_inst%t_h2osfc_col                , & ! Input:  [real(r8) (:)   ]  surface water temperature                                             
         beta                   => temperature_inst%beta_col                    , & ! Input:  [real(r8) (:)   ]  coefficient of conective velocity [-]                                 

         frac_sno               => waterstate_inst%frac_sno_col                 , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)                           
         qg_snow                => waterstate_inst%qg_snow_col                  , & ! Input:  [real(r8) (:)   ]  specific humidity at snow surface [kg/kg]                             
         qg_soil                => waterstate_inst%qg_soil_col                  , & ! Input:  [real(r8) (:)   ]  specific humidity at soil surface [kg/kg]                             
         qg_h2osfc              => waterstate_inst%qg_h2osfc_col                , & ! Input:  [real(r8) (:)   ]  specific humidity at h2osfc surface [kg/kg]                           
         qg                     => waterstate_inst%qg_col                       , & ! Input:  [real(r8) (:)   ]  specific humidity at ground surface [kg/kg]                           
         dqgdT                  => waterstate_inst%dqgdT_col                    , & ! Input:  [real(r8) (:)   ]  temperature derivative of "qg"                                        
         h2osoi_ice             => waterstate_inst%h2osoi_ice_col               , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)                                                    
         h2osoi_liq             => waterstate_inst%h2osoi_liq_col               , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)                                                
         grnd_ch4_cond          => ch4_inst%grnd_ch4_cond_patch                 , & ! Output: [real(r8) (:)   ]  tracer conductance for boundary layer [m/s]

         eflx_sh_snow           => energyflux_inst%eflx_sh_snow_patch           , & ! Output: [real(r8) (:)   ]  sensible heat flux from snow (W/m**2) [+ to atm]                      
         eflx_sh_soil           => energyflux_inst%eflx_sh_soil_patch           , & ! Output: [real(r8) (:)   ]  sensible heat flux from soil (W/m**2) [+ to atm]                      
         eflx_sh_h2osfc         => energyflux_inst%eflx_sh_h2osfc_patch         , & ! Output: [real(r8) (:)   ]  sensible heat flux from soil (W/m**2) [+ to atm]                      
         eflx_sh_grnd           => energyflux_inst%eflx_sh_grnd_patch           , & ! Output: [real(r8) (:)   ]  sensible heat flux from ground (W/m**2) [+ to atm]                    
         eflx_sh_tot            => energyflux_inst%eflx_sh_tot_patch            , & ! Output: [real(r8) (:)   ]  total sensible heat flux (W/m**2) [+ to atm]                          
         taux                   => energyflux_inst%taux_patch                   , & ! Output: [real(r8) (:)   ]  wind (shear) stress: e-w (kg/m/s**2)                                  
         tauy                   => energyflux_inst%tauy_patch                   , & ! Output: [real(r8) (:)   ]  wind (shear) stress: n-s (kg/m/s**2)                                  
         dlrad                  => energyflux_inst%dlrad_patch                  , & ! Output: [real(r8) (:)   ]  downward longwave radiation below the canopy [W/m2]                   
         ulrad                  => energyflux_inst%ulrad_patch                  , & ! Output: [real(r8) (:)   ]  upward longwave radiation above the canopy [W/m2]                     
         cgrnds                 => energyflux_inst%cgrnds_patch                 , & ! Output: [real(r8) (:)   ]  deriv, of soil sensible heat flux wrt soil temp [w/m2/k]              
         cgrndl                 => energyflux_inst%cgrndl_patch                 , & ! Output: [real(r8) (:)   ]  deriv of soil latent heat flux wrt soil temp [w/m**2/k]               
         cgrnd                  => energyflux_inst%cgrnd_patch                  , & ! Output: [real(r8) (:)   ]  deriv. of soil energy flux wrt to soil temp [w/m2/k]                  
         btran                  => energyflux_inst%btran_patch                  , & ! Output: [real(r8) (:)   ]  transpiration wetness factor (0 to 1)                                 
         rresis                 => energyflux_inst%rresis_patch                 , & ! Output: [real(r8) (:,:) ]  root resistance by layer (0-1)  (nlevgrnd)                          

         t_ref2m                => temperature_inst%t_ref2m_patch               , & ! Output: [real(r8) (:)   ]  2 m height surface air temperature (Kelvin)                           
         t_ref2m_r              => temperature_inst%t_ref2m_r_patch             , & ! Output: [real(r8) (:)   ]  Rural 2 m height surface air temperature (Kelvin)                     
         t_veg                  => temperature_inst%t_veg_patch                 , & ! Output: [real(r8) (:)   ]  vegetation temperature (Kelvin)                                       

         q_ref2m                => waterstate_inst%q_ref2m_patch                , & ! Output: [real(r8) (:)   ]  2 m height surface specific humidity (kg/kg)                          
         rh_ref2m_r             => waterstate_inst%rh_ref2m_r_patch             , & ! Output: [real(r8) (:)   ]  Rural 2 m height surface relative humidity (%)                        
         rh_ref2m               => waterstate_inst%rh_ref2m_patch               , & ! Output: [real(r8) (:)   ]  2 m height surface relative humidity (%)                              

         forc_hgt_u_patch       => frictionvel_inst%forc_hgt_u_patch            , & ! Input:
         u10_clm                => frictionvel_inst%u10_clm_patch               , & ! Input:  [real(r8) (:)   ]  10 m height winds (m/s)
         z0mg_col               => frictionvel_inst%z0mg_col                    , & ! Output: [real(r8) (:)   ]  roughness length, momentum [m]                                        
         z0hg_col               => frictionvel_inst%z0hg_col                    , & ! Output: [real(r8) (:)   ]  roughness length, sensible heat [m]                                   
         z0qg_col               => frictionvel_inst%z0qg_col                    , & ! Output: [real(r8) (:)   ]  roughness length, latent heat [m]                                     
         ram1                   => frictionvel_inst%ram1_patch                  , & ! Output: [real(r8) (:)   ]  aerodynamical resistance (s/m)                                        

         htvp                   => energyflux_inst%htvp_col                     , & ! Input:  [real(r8) (:)   ]  latent heat of evaporation (/sublimation) [J/kg]                      
         qflx_ev_snow           => waterflux_inst%qflx_ev_snow_patch            , & ! Output: [real(r8) (:)   ]  evaporation flux from snow (W/m**2) [+ to atm]                        
         qflx_ev_soil           => waterflux_inst%qflx_ev_soil_patch            , & ! Output: [real(r8) (:)   ]  evaporation flux from soil (W/m**2) [+ to atm]                        
         qflx_ev_h2osfc         => waterflux_inst%qflx_ev_h2osfc_patch          , & ! Output: [real(r8) (:)   ]  evaporation flux from h2osfc (W/m**2) [+ to atm]                      
         qflx_evap_soi          => waterflux_inst%qflx_evap_soi_patch           , & ! Output: [real(r8) (:)   ]  soil evaporation (mm H2O/s) (+ = to atm)                              
         qflx_evap_tot          => waterflux_inst%qflx_evap_tot_patch           , & ! Output: [real(r8) (:)   ]  qflx_evap_soi + qflx_evap_can + qflx_tran_veg                         

         rssun                  => photosyns_inst%rssun_patch                   , & ! Output: [real(r8) (:)   ]  leaf sunlit stomatal resistance (s/m) (output from Photosynthesis)
         rssha                  => photosyns_inst%rssha_patch                   , & ! Output: [real(r8) (:)   ]  leaf shaded stomatal resistance (s/m) (output from Photosynthesis)

         begp                   => bounds%begp                                  , &
         endp                   => bounds%endp                                    &
         )

      ! First do some simple settings of values over points where frac vegetation covered
      ! by snow is zero

      do f = 1, num_noexposedvegp
         p = filter_noexposedvegp(f)
         c = patch%column(p)
         btran(p) = 0._r8     
         t_veg(p) = forc_t(c) 
         cf_bare  = forc_pbot(c)/(SHR_CONST_RGAS*0.001_r8*thm(p))*1.e06_r8
         rssun(p) = 1._r8/1.e15_r8 * cf_bare
         rssha(p) = 1._r8/1.e15_r8 * cf_bare
         do j = 1, nlevgrnd
            rootr(p,j)  = 0._r8
            rresis(p,j) = 0._r8
         end do
      end do

      ! Compute sensible and latent fluxes and their derivatives with respect
      ! to ground temperature using ground temperatures from previous time step

      do f = 1, num_noexposedvegp
         p = filter_noexposedvegp(f)
         c = patch%column(p)
         g = patch%gridcell(p)

         ! Initialization variables

         displa(p) = 0._r8
         dlrad(p)  = 0._r8
         ulrad(p)  = 0._r8

         ur(p)    = max(1.0_r8,sqrt(forc_u(g)*forc_u(g)+forc_v(g)*forc_v(g)))
         dth(p)   = thm(p)-t_grnd(c)
         dqh(p)   = forc_q(c) - qg(c)
         dthv     = dth(p)*(1._r8+0.61_r8*forc_q(c))+0.61_r8*forc_th(c)*dqh(p)
         zldis(p) = forc_hgt_u_patch(p)

         ! Copy column roughness to local patch-level arrays

         z0mg_patch(p) = z0mg_col(c)
         z0hg_patch(p) = z0hg_col(c)
         z0qg_patch(p) = z0qg_col(c)

         ! Initialize Monin-Obukhov length and wind speed

         call MoninObukIni(ur(p), thv(c), dthv, zldis(p), z0mg_patch(p), um(p), obu(p))

      end do

      ! Perform stability iteration
      ! Determine friction velocity, and potential temperature and humidity
      ! profiles of the surface boundary layer

      do iter = 1, niters

         call FrictionVelocity(begp, endp, num_noexposedvegp, filter_noexposedvegp, &
              displa(begp:endp), z0mg_patch(begp:endp), z0hg_patch(begp:endp), z0qg_patch(begp:endp), &
              obu(begp:endp), iter, ur(begp:endp), um(begp:endp), ustar(begp:endp), &
              temp1(begp:endp), temp2(begp:endp), temp12m(begp:endp), temp22m(begp:endp), fm(begp:endp), &
              frictionvel_inst)

         do f = 1, num_noexposedvegp
            p = filter_noexposedvegp(f)
            c = patch%column(p)
            g = patch%gridcell(p)

            tstar = temp1(p)*dth(p)
            qstar = temp2(p)*dqh(p)
            z0hg_patch(p) = z0mg_patch(p)/exp(0.13_r8 * (ustar(p)*z0mg_patch(p)/1.5e-5_r8)**0.45_r8)
            z0qg_patch(p) = z0hg_patch(p)
            thvstar = tstar*(1._r8+0.61_r8*forc_q(c)) + 0.61_r8*forc_th(c)*qstar
            zeta = zldis(p)*vkc*grav*thvstar/(ustar(p)**2*thv(c))

            if (zeta >= 0._r8) then                   !stable
               zeta = min(2._r8,max(zeta,0.01_r8))
               um(p) = max(ur(p),0.1_r8)
            else                                      !unstable
               zeta = max(-100._r8,min(zeta,-0.01_r8))
               wc = beta(c)*(-grav*ustar(p)*thvstar*zii(c)/thv(c))**0.333_r8
               um(p) = sqrt(ur(p)*ur(p) + wc*wc)
            end if
            obu(p) = zldis(p)/zeta
         end do

      end do ! end stability iteration

      do f = 1, num_noexposedvegp
         p = filter_noexposedvegp(f)
         c = patch%column(p)
         g = patch%gridcell(p)
         l = patch%landunit(p)

         ! Determine aerodynamic resistances

         ram  = 1._r8/(ustar(p)*ustar(p)/um(p))
         rah  = 1._r8/(temp1(p)*ustar(p))
         raw  = 1._r8/(temp2(p)*ustar(p))
         raih = forc_rho(c)*cpair/rah
         if (use_lch4) then
            grnd_ch4_cond(p) = 1._r8/raw
         end if

         ! Soil evaporation resistance
         www = (h2osoi_liq(c,1)/denh2o+h2osoi_ice(c,1)/denice)/dz(c,1)/watsat(c,1)
         www = min(max(www,0.0_r8),1._r8)

         !changed by K.Sakaguchi. Soilbeta is used for evaporation
         if (dqh(p) > 0._r8) then  !dew  (beta is not applied, just like rsoil used to be) 
            raiw = forc_rho(c)/(raw)
         else
            if(do_soilevap_beta())then
               ! Lee and Pielke 1992 beta is applied
               raiw    = soilbeta(c)*forc_rho(c)/(raw)
            endif
         end if

         ram1(p) = ram  !pass value to global variable

         ! Output to patch-level data structures
         ! Derivative of fluxes with respect to ground temperature
         cgrnds(p) = raih
         cgrndl(p) = raiw*dqgdT(c)
         cgrnd(p)  = cgrnds(p) + htvp(c)*cgrndl(p)


         ! Variables needed by history tape

         ! Surface fluxes of momentum, sensible and latent heat
         ! using ground temperatures from previous time step
         taux(p)          = -forc_rho(c)*forc_u(g)/ram
         tauy(p)          = -forc_rho(c)*forc_v(g)/ram
         eflx_sh_grnd(p)  = -raih*dth(p)
         eflx_sh_tot(p)   = eflx_sh_grnd(p)

         ! compute sensible heat fluxes individually
         eflx_sh_snow(p)   = -raih*(thm(p)-t_soisno(c,snl(c)+1))
         eflx_sh_soil(p)   = -raih*(thm(p)-t_soisno(c,1))
         eflx_sh_h2osfc(p) = -raih*(thm(p)-t_h2osfc(c))

         ! water fluxes from soil
         qflx_evap_soi(p)  = -raiw*dqh(p)
         qflx_evap_tot(p)  = qflx_evap_soi(p)

         ! compute latent heat fluxes individually
         qflx_ev_snow(p)   = -raiw*(forc_q(c) - qg_snow(c))
         qflx_ev_soil(p)   = -raiw*(forc_q(c) - qg_soil(c))
         qflx_ev_h2osfc(p) = -raiw*(forc_q(c) - qg_h2osfc(c))

         ! 2 m height air temperature
         t_ref2m(p) = thm(p) + temp1(p)*dth(p)*(1._r8/temp12m(p) - 1._r8/temp1(p))

         ! 2 m height specific humidity
         q_ref2m(p) = forc_q(c) + temp2(p)*dqh(p)*(1._r8/temp22m(p) - 1._r8/temp2(p))

         ! 2 m height relative humidity
         call QSat(t_ref2m(p), forc_pbot(c), e_ref2m, de2mdT, qsat_ref2m, dqsat2mdT)

         rh_ref2m(p) = min(100._r8, q_ref2m(p) / qsat_ref2m * 100._r8)

         if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
            rh_ref2m_r(p) = rh_ref2m(p)
            t_ref2m_r(p) = t_ref2m(p)
         end if

         ! Human Heat Stress
         if ( calc_human_stress_indices )then
            call KtoC(t_ref2m(p), tc_ref2m(p))
            call VaporPres(rh_ref2m(p), e_ref2m, vap_ref2m(p))
            call Wet_Bulb(t_ref2m(p), vap_ref2m(p), forc_pbot(c), rh_ref2m(p), q_ref2m(p), &
                            teq_ref2m(p), ept_ref2m(p), wb_ref2m(p))
            call Wet_BulbS(tc_ref2m(p),rh_ref2m(p), wbt_ref2m(p))
            call HeatIndex(tc_ref2m(p), rh_ref2m(p), nws_hi_ref2m(p))
            call AppTemp(tc_ref2m(p), vap_ref2m(p), u10_clm(p), appar_temp_ref2m(p))
            call swbgt(tc_ref2m(p), vap_ref2m(p), swbgt_ref2m(p))
            call hmdex(tc_ref2m(p), vap_ref2m(p), humidex_ref2m(p))
            call dis_coi(tc_ref2m(p), wb_ref2m(p), discomf_index_ref2m(p))
            call dis_coiS(tc_ref2m(p), rh_ref2m(p), wbt_ref2m(p), discomf_index_ref2mS(p))
            call THIndex(tc_ref2m(p), wb_ref2m(p), thic_ref2m(p), thip_ref2m(p))
            call SwampCoolEff(tc_ref2m(p), wb_ref2m(p), swmp80_ref2m(p), swmp65_ref2m(p))
  
            if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
              teq_ref2m_r(p)            = teq_ref2m(p)
              ept_ref2m_r(p)            = ept_ref2m(p)
              wb_ref2m_r(p)             = wb_ref2m(p)
              wbt_ref2m_r(p)            = wbt_ref2m(p)
              nws_hi_ref2m_r(p)         = nws_hi_ref2m(p)
              appar_temp_ref2m_r(p)     = appar_temp_ref2m(p)
              swbgt_ref2m_r(p)          = swbgt_ref2m(p)
              humidex_ref2m_r(p)        = humidex_ref2m(p)
              discomf_index_ref2m_r(p)  = discomf_index_ref2m(p)
              discomf_index_ref2mS_r(p) = discomf_index_ref2mS(p)
              thic_ref2m_r(p)           = thic_ref2m(p)
              thip_ref2m_r(p)           = thip_ref2m(p)
              swmp80_ref2m_r(p)         = swmp80_ref2m(p)
              swmp65_ref2m_r(p)         = swmp65_ref2m(p)
            end if
         end if
      end do

    end associate

  end subroutine BareGroundFluxes

end module BareGroundFluxesMod
