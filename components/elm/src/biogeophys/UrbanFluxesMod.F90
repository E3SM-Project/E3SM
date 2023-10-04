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
  use elm_varpar           , only : numrad
  use elm_varcon           , only : isecspday, degpsec, namel
  use elm_varctl           , only : iulog
  use abortutils           , only : endrun
  use UrbanParamsType      , only : urbanparams_type
  use UrbanParamsType      , only : urban_wasteheat_on, urban_hac_on, urban_hac
  use atm2lndType          , only : atm2lnd_type
  use SoilStateType        , only : soilstate_type
  use FrictionVelocityType , only : frictionvel_type
  use EnergyFluxType       , only : energyflux_type
  use SurfaceResistanceMod , only : do_soilevap_beta
  use GridcellType         , only : grc_pp
  use TopounitDataType     , only : top_as
  use LandunitType         , only : lun_pp
  use LandunitDataType     , only : lun_es, lun_ef, lun_ws
  use ColumnType           , only : col_pp
  use ColumnDataType       , only : col_es, col_ef, col_ws
  use VegetationType       , only : veg_pp
  use VegetationDataType   , only : veg_es, veg_ef, veg_ws, veg_wf
  use clm_time_manager    , only : get_curr_date, get_step_size, get_nstep

  use timeinfoMod  , only : nstep_mod, year_curr, mon_curr, day_curr, secs_curr
  use timeinfoMod  , only : dtime_mod

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
       urbanparams_vars, soilstate_vars,    &
      frictionvel_vars )
    !
    ! !DESCRIPTION:
    ! Turbulent and momentum fluxes from urban canyon (consisting of roof, sunwall,
    ! shadewall, pervious and impervious road).

    ! !USES:
    use shr_flux_mod         , only : shr_flux_update_stress
    use elm_varcon          , only : cpair, vkc, spval, grav, pondmx_urban, rpi, rgas
    use elm_varcon          , only : ht_wasteheat_factor, ac_wasteheat_factor, wasteheat_limit
    use column_varcon       , only : icol_shadewall, icol_road_perv, icol_road_imperv
    use column_varcon       , only : icol_roof, icol_sunwall
    use filterMod           , only : filter
    use FrictionVelocityMod , only : FrictionVelocity_loops, MoninObukIni, implicit_stress
    use QSatMod             , only : QSat
    use elm_varpar          , only : maxpatch_urb, nlevurb, nlevgrnd
    use elm_varctl          , only : use_vsfm
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
    type(urbanparams_type) , intent(in)    :: urbanparams_vars
    type(soilstate_type)   , intent(inout) :: soilstate_vars
    type(frictionvel_type) , intent(inout) :: frictionvel_vars
    real(r8) :: dtime                                                ! land model time step (sec)
    integer  :: year,month,day,secs
    !
    ! !LOCAL VARIABLES:
    integer  :: fp,fc,fl,f,p,c,l,t,g,j,pi,i     ! indices

   !  integer  :: filter_copyl(num_urbanl)                ! iteration copy of filter_urbanl
   !  integer  :: filter_copyc(num_urbanc)                ! iteration copy of filter_urbanc
    integer  :: num_copyl                                ! iteration num_urbanl
    integer  :: num_copyl_old                            ! previous iteration num_copyl
    integer  :: num_copyc                                ! iteration num_urbanc
    integer  :: num_copyc_old                            ! previous iteration num_copyc

    real(r8) :: canyontop_wind(1:num_urbanl)              ! wind at canyon top (m/s) 
    real(r8) :: canyon_u_wind(1:num_urbanl)               ! u-component of wind speed inside canyon (m/s)
    real(r8) :: canyon_wind(1:num_urbanl)                 ! net wind speed inside canyon (m/s)
    real(r8) :: canyon_resistance(1:num_urbanl)           ! resistance to heat and moisture transfer from canyon road/walls to canyon air (s/m)

    real(r8) :: ur(1:num_urbanl)                          ! wind speed at reference height (m/s)
    real(r8) :: ustar(1:num_urbanl)                       ! friction velocity (m/s)
    real(r8) :: ramu(1:num_urbanl)                        ! aerodynamic resistance (s/m)
    real(r8) :: rahu(1:num_urbanl)                        ! thermal resistance (s/m)
    real(r8) :: rawu(1:num_urbanl)                        ! moisture resistance (s/m)
    real(r8) :: temp1(1:num_urbanl)                       ! relation for potential temperature profile
    real(r8) :: temp12m(1:num_urbanl)                     ! relation for potential temperature profile applied at 2-m
    real(r8) :: temp2(1:num_urbanl)                       ! relation for specific humidity profile
    real(r8) :: temp22m(1:num_urbanl)                     ! relation for specific humidity profile applied at 2-m
    real(r8) :: thm_g(1:num_urbanl)                       ! intermediate variable (forc_t+0.0098*forc_hgt_t)
    real(r8) :: thv_g(1:num_urbanl)                       ! virtual potential temperature (K)
    real(r8) :: dth(1:num_urbanl)                         ! diff of virtual temp. between ref. height and surface
    real(r8) :: dqh(1:num_urbanl)                         ! diff of humidity between ref. height and surface
    real(r8) :: zldis(1:num_urbanl)                       ! reference height "minus" zero displacement height (m)
    real(r8) :: um(1:num_urbanl)                          ! wind speed including the stablity effect (m/s)
    real(r8) :: obu(1:num_urbanl)                         ! Monin-Obukhov length (m)
    real(r8) :: taf_numer(1:num_urbanl)                   ! numerator of taf equation (K m/s)
    real(r8) :: taf_denom(1:num_urbanl)                   ! denominator of taf equation (m/s)
    real(r8) :: taf_numer_test(1:num_urbanl)                   ! numerator of taf equation (K m/s)
    real(r8) :: taf_denom_test(1:num_urbanl)                   ! denominator of taf equation (m/s)
    real(r8) :: qaf_numer(1:num_urbanl)                   ! numerator of qaf equation (kg m/kg s)
    real(r8) :: qaf_denom(1:num_urbanl)                   ! denominator of qaf equation (m/s)
    real(r8) :: wtas(1:num_urbanl)                        ! sensible heat conductance for urban air to atmospheric air (m/s)
    real(r8) :: wtaq(1:num_urbanl)                        ! latent heat conductance for urban air to atmospheric air (m/s)
    real(r8) :: wts_sum(1:num_urbanl)                     ! sum of wtas, wtus_roof, wtus_road_perv, wtus_road_imperv, wtus_sunwall, wtus_shadewall
    real(r8) :: wtq_sum(1:num_urbanl)                     ! sum of wtaq, wtuq_roof, wtuq_road_perv, wtuq_road_imperv, wtuq_sunwall, wtuq_shadewall
    real(r8) :: fm(1:num_urbanl)                          ! needed for BGC only to diagnose 10m wind speed
    real(r8) :: wtus(1:num_urbanc)                        ! sensible heat conductance for urban columns (scaled) (m/s)
    real(r8) :: wtuq(1:num_urbanc)                        ! latent heat conductance for urban columns (scaled) (m/s)
    integer  :: iter                                                 ! iteration index
    integer  :: iter_final                                           ! number of iterations used
    real(r8) :: dthv                                                 ! diff of vir. poten. temp. between ref. height and surface
    real(r8) :: tstar                                                ! temperature scaling parameter
    real(r8) :: qstar                                                ! moisture scaling parameter
    real(r8) :: thvstar                                              ! virtual potential temperature scaling parameter
    real(r8) :: wtus_roof(1:num_urbanl)                   ! sensible heat conductance for roof (scaled) (m/s)
    real(r8) :: wtuq_roof(1:num_urbanl)                   ! latent heat conductance for roof (scaled) (m/s)
    real(r8) :: wtus_road_perv(1:num_urbanl)              ! sensible heat conductance for pervious road (scaled) (m/s)
    real(r8) :: wtuq_road_perv(1:num_urbanl)              ! latent heat conductance for pervious road (scaled) (m/s)
    real(r8) :: wtus_road_imperv(1:num_urbanl)            ! sensible heat conductance for impervious road (scaled) (m/s)
    real(r8) :: wtuq_road_imperv(1:num_urbanl)            ! latent heat conductance for impervious road (scaled) (m/s)
    real(r8) :: wtus_sunwall(1:num_urbanl)                ! sensible heat conductance for sunwall (scaled) (m/s)
    real(r8) :: wtuq_sunwall(1:num_urbanl)                ! latent heat conductance for sunwall (scaled) (m/s)
    real(r8) :: wtus_shadewall(1:num_urbanl)              ! sensible heat conductance for shadewall (scaled) (m/s)
    real(r8) :: wtuq_shadewall(1:num_urbanl)              ! latent heat conductance for shadewall (scaled) (m/s)
    real(r8) :: wtus_roof_unscl(1:num_urbanl)             ! sensible heat conductance for roof (not scaled) (m/s)
    real(r8) :: wtuq_roof_unscl(1:num_urbanl)             ! latent heat conductance for roof (not scaled) (m/s)
    real(r8) :: wtus_road_perv_unscl(1:num_urbanl)        ! sensible heat conductance for pervious road (not scaled) (m/s)
    real(r8) :: wtuq_road_perv_unscl(1:num_urbanl)        ! latent heat conductance for pervious road (not scaled) (m/s)
    real(r8) :: wtus_road_imperv_unscl(1:num_urbanl)      ! sensible heat conductance for impervious road (not scaled) (m/s)
    real(r8) :: wtuq_road_imperv_unscl(1:num_urbanl)      ! latent heat conductance for impervious road (not scaled) (m/s)
    real(r8) :: wtus_sunwall_unscl(1:num_urbanl)          ! sensible heat conductance for sunwall (not scaled) (m/s)
    real(r8) :: wtuq_sunwall_unscl(1:num_urbanl)          ! latent heat conductance for sunwall (not scaled) (m/s)
    real(r8) :: wtus_shadewall_unscl(1:num_urbanl)        ! sensible heat conductance for shadewall (not scaled) (m/s)
    real(r8) :: wtuq_shadewall_unscl(1:num_urbanl)        ! latent heat conductance for shadewall (not scaled) (m/s)
    real(r8) :: t_sunwall_innerl(1:num_urbanl)            ! temperature of inner layer of sunwall (K)
    real(r8) :: t_shadewall_innerl(1:num_urbanl)          ! temperature of inner layer of shadewall (K)
    real(r8) :: t_roof_innerl(1:num_urbanl)               ! temperature of inner layer of roof (K)
    real(r8) :: lngth_roof                                           ! length of roof (m)
    real(r8) :: wc                                                   ! convective velocity (m/s)
    real(r8) :: zeta                                                 ! dimensionless height used in Monin-Obukhov theory
    real(r8) :: eflx_sh_grnd_scale(bounds%begp:bounds%endp)          ! scaled sensible heat flux from ground (W/m**2) [+ to atm]
    real(r8) :: qflx_evap_soi_scale(bounds%begp:bounds%endp)         ! scaled soil evaporation (mm H2O/s) (+ = to atm)
    real(r8) :: eflx_wasteheat_roof(1:num_urbanl)         ! sensible heat flux from urban heating/cooling sources of waste heat for roof (W/m**2)
    real(r8) :: eflx_wasteheat_sunwall(1:num_urbanl)      ! sensible heat flux from urban heating/cooling sources of waste heat for sunwall (W/m**2)
    real(r8) :: eflx_wasteheat_shadewall(1:num_urbanl)    ! sensible heat flux from urban heating/cooling sources of waste heat for shadewall (W/m**2)
    real(r8) :: eflx_heat_from_ac_roof(1:num_urbanl)      ! sensible heat flux put back into canyon due to heat removal by AC for roof (W/m**2)
    real(r8) :: eflx_heat_from_ac_sunwall(1:num_urbanl)   ! sensible heat flux put back into canyon due to heat removal by AC for sunwall (W/m**2)
    real(r8) :: eflx_heat_from_ac_shadewall(1:num_urbanl) ! sensible heat flux put back into canyon due to heat removal by AC for shadewall (W/m**2)
    real(r8) :: eflx(1:num_urbanl)                        ! total sensible heat flux for error check (W/m**2)
    real(r8) :: qflx(1:num_urbanl)                        ! total water vapor flux for error check (kg/m**2/s)
    real(r8) :: eflx_scale(1:num_urbanl)                  ! sum of scaled sensible heat fluxes for urban columns for error check (W/m**2)
    real(r8) :: qflx_scale(1:num_urbanl)                  ! sum of scaled water vapor fluxes for urban columns for error check (kg/m**2/s)
    real(r8) :: eflx_err(1:num_urbanl)                    ! sensible heat flux error (W/m**2)
    real(r8) :: qflx_err(1:num_urbanl)                    ! water vapor flux error (kg/m**2/s)
    real(r8) :: fwet_roof                                            ! fraction of roof surface that is wet (-)
    real(r8) :: fwet_road_imperv                                     ! fraction of impervious road surface that is wet (-)
    integer  :: local_secp1(1:num_urbanl)                 ! seconds into current date in local time (sec)
                                     ! calendar info for current time step
    logical  :: found                                                ! flag in search loop
    integer  :: indexl                                               ! index of first found in search loop
    integer  :: nstep                                                ! time step number
    real(r8) :: e_ref2m                                              ! 2 m height surface saturated vapor pressure [Pa]
    real(r8) :: de2mdT                                               ! derivative of 2 m height surface saturated vapor pressure on t_ref2m
    real(r8) :: qsat_ref2m                                           ! 2 m height surface saturated specific humidity [kg/kg]
    real(r8) :: dqsat2mdT                                            ! derivative of 2 m height surface saturated specific humidity on t_ref2m
    !
    real(r8), parameter :: lapse_rate = 0.0098_r8 ! Dry adiabatic lapse rate (K/m)
    real(r8), parameter :: dtaumin = 0.01_r8      ! max limit for stress convergence [Pa]
    integer, parameter  :: itmin = 3              ! minimum number of iterations
    integer, parameter  :: itmax = 30             ! maximum number of iterations
    integer  :: loopmax                           ! bound for iteration loop
    real(r8) :: wind_speed0(1:num_urbanl) ! Wind speed from atmosphere at start of iteration
    real(r8) :: wind_speed_adj(1:num_urbanl) ! Adjusted wind speed for iteration
    real(r8) :: tau(1:num_urbanl)      ! Stress used in iteration
    real(r8) :: tau_diff(1:num_urbanl) ! Difference from previous iteration tau
    real(r8) :: prev_tau(1:num_urbanl) ! Previous iteration tau
    real(r8) :: prev_tau_diff(bounds%begl:bounds%endl) ! Previous difference in iteration tau
    real(r8), parameter :: beta = 1._r8           ! coefficient of convective velocity
    real(r8), parameter :: zii  = 1000._r8        ! convective boundary layer height (m)
    integer :: lnd_to_urban_filter(bounds%begl:bounds%endl) !
    integer :: col_to_urban_filter(bounds%begc:bounds%endc)
    logical :: converged_landunits(bounds%begl:bounds%endl)
    integer :: begl, endl, begc,endc 
    real(r8) :: sum_denom, sum_numer
    !-----------------------------------------------------------------------

    associate(                                                                & 
         snl                 =>   col_pp%snl                                   , & ! Input:  [integer  (:)   ]  number of snow layers                              
         ctype               =>   col_pp%itype                                 , & ! Input:  [integer  (:)   ]  column type                                        
         z_0_town            =>   lun_pp%z_0_town                              , & ! Input:  [real(r8) (:)   ]  momentum roughness length of urban landunit (m)   
         z_d_town            =>   lun_pp%z_d_town                              , & ! Input:  [real(r8) (:)   ]  displacement height of urban landunit (m)         
         ht_roof             =>   lun_pp%ht_roof                               , & ! Input:  [real(r8) (:)   ]  height of urban roof (m)                          
         wtlunit_roof        =>   lun_pp%wtlunit_roof                          , & ! Input:  [real(r8) (:)   ]  weight of roof with respect to landunit           
         canyon_hwr          =>   lun_pp%canyon_hwr                            , & ! Input:  [real(r8) (:)   ]  ratio of building height to street width          
         wtroad_perv         =>   lun_pp%wtroad_perv                           , & ! Input:  [real(r8) (:)   ]  weight of pervious road wrt total road            

         forc_t              =>   top_as%tbot                               , & ! Input:  [real(r8) (:)   ]  atmospheric temperature (K)                       
         forc_th             =>   top_as%thbot                              , & ! Input:  [real(r8) (:)   ]  atmospheric potential temperature (K)             
         forc_rho            =>   top_as%rhobot                             , & ! Input:  [real(r8) (:)   ]  air density (kg/m**3)                                 
         forc_q              =>   top_as%qbot                               , & ! Input:  [real(r8) (:)   ]  atmospheric specific humidity (kg/kg)             
         forc_pbot           =>   top_as%pbot                               , & ! Input:  [real(r8) (:)   ]  atmospheric pressure (Pa)                         
         forc_u              =>   top_as%ubot                               , & ! Input:  [real(r8) (:)   ]  atmospheric wind speed in east direction (m/s)    
         forc_v              =>   top_as%vbot                               , & ! Input:  [real(r8) (:)   ]  atmospheric wind speed in north direction (m/s)   
         wsresp              =>   top_as%wsresp                             , & ! Input:  [real(r8) (:)   ]  response of wind to surface stress (m/s/Pa)
         tau_est             =>   top_as%tau_est                            , & ! Input:  [real(r8) (:)   ]  approximate atmosphere change to zonal wind (m/s)
         ugust               =>   top_as%ugust                              , & ! Input:  [real(r8) (:)   ]  gustiness from atmosphere (m/s)

         wind_hgt_canyon     =>   urbanparams_vars%wind_hgt_canyon          , & ! Input:  [real(r8) (:)   ]  height above road at which wind in canyon is to be computed (m)
         eflx_traffic_factor =>   urbanparams_vars%eflx_traffic_factor      , & ! Input:  [real(r8) (:)   ]  multiplicative urban traffic factor for sensible heat flux

         rootr_road_perv     =>   soilstate_vars%rootr_road_perv_col        , & ! Input:  [real(r8) (:,:) ]  effective fraction of roots in each soil layer for urban pervious road
         soilalpha_u         =>   soilstate_vars%soilalpha_u_col            , & ! Input:  [real(r8) (:)   ]  Urban factor that reduces ground saturated specific humidity (-)
         soilbeta            =>   soilstate_vars%soilbeta_col               , & ! Input:  [real(r8) (:)   ]  soil wetness relative to field capacity
         rootr               =>   soilstate_vars%rootr_patch                , & ! Output: [real(r8) (:,:) ]  effective fraction of roots in each soil layer

         t_grnd              =>   col_es%t_grnd               , & ! Input:  [real(r8) (:)   ]  ground surface temperature (K)
         t_soisno            =>   col_es%t_soisno             , & ! Input:  [real(r8) (:,:) ]  soil temperature (K)
         t_ref2m             =>   veg_es%t_ref2m            , & ! Output: [real(r8) (:)   ]  2 m height surface air temperature (K)
         t_ref2m_u           =>   veg_es%t_ref2m_u          , & ! Output: [real(r8) (:)   ]  Urban 2 m height surface air temperature (K)
         t_veg               =>   veg_es%t_veg                , & ! Output: [real(r8) (:)   ]  vegetation temperature (K)
         t_building          =>   lun_es%t_building           , & ! Output: [real(r8) (:)   ]  internal building temperature (K)
         taf                 =>   lun_es%taf                  , & ! Output: [real(r8) (:)   ]  urban canopy air temperature (K)

         frac_sno            =>   col_ws%frac_sno              , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)
         snow_depth          =>   col_ws%snow_depth            , & ! Input:  [real(r8) (:)   ]  snow height (m)
         dqgdT               =>   col_ws%dqgdT                 , & ! Input:  [real(r8) (:)   ]  temperature derivative of "qg"
         qg                  =>   col_ws%qg                    , & ! Input:  [real(r8) (:)   ]  specific humidity at ground surface (kg/kg)
         h2osoi_ice          =>   col_ws%h2osoi_ice            , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)
         h2osoi_liq          =>   col_ws%h2osoi_liq            , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)
         h2osno              =>   col_ws%h2osno                , & ! Input:  [real(r8) (:)   ]  snow water (mm H2O)
         qaf                 =>   lun_ws%qaf                   , & ! Output: [real(r8) (:)   ]  urban canopy air specific humidity (kg/kg)
         q_ref2m             =>   veg_ws%q_ref2m             , & ! Output: [real(r8) (:)   ]  2 m height surface specific humidity (kg/kg)
         rh_ref2m            =>   veg_ws%rh_ref2m            , & ! Output: [real(r8) (:)   ]  2 m height surface relative humidity (%)
         rh_ref2m_u          =>   veg_ws%rh_ref2m_u          , & ! Output: [real(r8) (:)   ]  2 m height surface relative humidity (%)

         forc_hgt_u_patch    =>   frictionvel_vars%forc_hgt_u_patch         , & ! Input:  [real(r8) (:)   ]  observational height of wind at pft-level (m)
         forc_hgt_t_patch    =>   frictionvel_vars%forc_hgt_t_patch         , & ! Input:  [real(r8) (:)   ]  observational height of temperature at pft-level (m)
         ram1                =>   frictionvel_vars%ram1_patch               , & ! Output: [real(r8) (:)   ]  aerodynamical resistance (s/m)

         htvp                =>   col_ef%htvp                  , & ! Input:  [real(r8) (:)   ]  latent heat of evaporation (/sublimation) (J/kg)
         eflx_urban_ac       =>   col_ef%eflx_urban_ac         , & ! Input:  [real(r8) (:)   ]  urban air conditioning flux (W/m**2)
         eflx_urban_heat     =>   col_ef%eflx_urban_heat       , & ! Input:  [real(r8) (:)   ]  urban heating flux (W/m**2)
         dlrad               =>   veg_ef%dlrad               , & ! Output: [real(r8) (:)   ]  downward longwave radiation below the canopy (W/m**2)
         ulrad               =>   veg_ef%ulrad               , & ! Output: [real(r8) (:)   ]  upward longwave radiation above the canopy (W/m**2)
         cgrnds              =>   veg_ef%cgrnds              , & ! Output: [real(r8) (:)   ]  deriv, of soil sensible heat flux wrt soil temp (W/m**2/K)
         cgrndl              =>   veg_ef%cgrndl              , & ! Output: [real(r8) (:)   ]  deriv of soil latent heat flux wrt soil temp (W/m**2/K)
         cgrnd               =>   veg_ef%cgrnd               , & ! Output: [real(r8) (:)   ]  deriv. of soil energy flux wrt to soil temp (W/m**2/K)
         eflx_sh_grnd        =>   veg_ef%eflx_sh_grnd        , & ! Output: [real(r8) (:)   ]  sensible heat flux from ground (W/m**2) [+ to atm]
         eflx_sh_tot         =>   veg_ef%eflx_sh_tot         , & ! Output: [real(r8) (:)   ]  total sensible heat flux (W/m**2) [+ to atm]
         eflx_sh_tot_u       =>   veg_ef%eflx_sh_tot_u       , & ! Output: [real(r8) (:)   ]  urban total sensible heat flux (W/m**2) [+ to atm]
         eflx_sh_snow        =>   veg_ef%eflx_sh_snow        , & ! Output: [real(r8) (:)   ]  sensible heat flux from snow (W/m**2) [+ to atm]
         eflx_sh_soil        =>   veg_ef%eflx_sh_soil        , & ! Output: [real(r8) (:)   ]  sensible heat flux from soil (W/m**2) [+ to atm]
         eflx_sh_h2osfc      =>   veg_ef%eflx_sh_h2osfc      , & ! Output: [real(r8) (:)   ]  sensible heat flux from soil (W/m**2) [+ to atm]
         eflx_traffic        =>   lun_ef%eflx_traffic          , & ! Output: [real(r8) (:)   ]  traffic sensible heat flux (W/m**2)
         eflx_wasteheat      =>   lun_ef%eflx_wasteheat        , & ! Output: [real(r8) (:)   ]  sensible heat flux from urban heating/cooling sources of waste heat (W/m**2)
         eflx_heat_from_ac   =>   lun_ef%eflx_heat_from_ac     , & ! Output: [real(r8) (:)   ]  sensible heat flux put back into canyon due to removal by AC (W/m**2)
         taux                =>   veg_ef%taux                , & ! Output: [real(r8) (:)   ]  wind (shear) stress: e-w (kg/m/s**2)
         tauy                =>   veg_ef%tauy                , & ! Output: [real(r8) (:)   ]  wind (shear) stress: n-s (kg/m/s**2)

         qflx_evap_soi       =>   veg_wf%qflx_evap_soi        , & ! Output: [real(r8) (:)   ]  soil evaporation (mm H2O/s) (+ = to atm)
         qflx_tran_veg       =>   veg_wf%qflx_tran_veg        , & ! Output: [real(r8) (:)   ]  vegetation transpiration (mm H2O/s) (+ = to atm)
         qflx_evap_veg       =>   veg_wf%qflx_evap_veg        , & ! Output: [real(r8) (:)   ]  vegetation evaporation (mm H2O/s) (+ = to atm)
         qflx_evap_tot       =>   veg_wf%qflx_evap_tot         & ! Output: [real(r8) (:)   ]  qflx_evap_soi + qflx_evap_can + qflx_tran_veg

        
         )
    !$acc enter data create(&
    !$acc canyontop_wind(:), &
    !$acc canyon_u_wind(:), &
    !$acc canyon_wind(:), &
    !$acc canyon_resistance(:), &
    !$acc ur(:), &
    !$acc ustar(:), &
    !$acc ramu(:), &
    !$acc rahu(:), &
    !$acc rawu(:), &
    !$acc temp1(:), &
    !$acc temp12m(:), &
    !$acc temp2(:), &
    !$acc temp22m(:), &
    !$acc thm_g(:), &
    !$acc thv_g(:), &
    !$acc dth(:), &
    !$acc dqh(:), &
    !$acc zldis(:), &
    !$acc um(:), &
    !$acc obu(:), &
    !$acc taf_numer(:), &
    !$acc taf_denom(:), &
    !$acc qaf_numer(:), &
    !$acc qaf_denom(:), &
    !$acc wtas(:), &
    !$acc wtaq(:), &
    !$acc wts_sum(:), &
    !$acc wtq_sum(:), &
    !$acc fm(:), &
    !$acc wtus(:), &
    !$acc wtuq(:), &
    !$acc wtus_roof(:), &
    !$acc wtuq_roof(:), &
    !$acc wtus_road_perv(:), &
    !$acc wtuq_road_perv(:), &
    !$acc wtus_road_imperv(:), &
    !$acc wtuq_road_imperv(:), &
    !$acc wtus_sunwall(:), &
    !$acc wtuq_sunwall(:), &
    !$acc wtus_shadewall(:), &
    !$acc wtuq_shadewall(:), &
    !$acc wtus_roof_unscl(:), &
    !$acc wtuq_roof_unscl(:), &
    !$acc wtus_road_perv_unscl(:), &
    !$acc wtuq_road_perv_unscl(:), &
    !$acc wtus_road_imperv_unscl(:), &
    !$acc wtuq_road_imperv_unscl(:), &
    !$acc wtus_sunwall_unscl(:), &
    !$acc wtuq_sunwall_unscl(:), &
    !$acc wtus_shadewall_unscl(:), &
    !$acc wtuq_shadewall_unscl(:), &
    !$acc t_sunwall_innerl(:), &
    !$acc t_shadewall_innerl(:), &
    !$acc t_roof_innerl(:), &
    !$acc eflx_sh_grnd_scale(:), &
    !$acc qflx_evap_soi_scale(:), &
    !$acc eflx_wasteheat_roof(:), &
    !$acc eflx_wasteheat_sunwall(:), &
    !$acc eflx_wasteheat_shadewall(:), &
    !$acc eflx_heat_from_ac_roof(:), &
    !$acc eflx_heat_from_ac_sunwall(:), &
    !$acc eflx_heat_from_ac_shadewall(:), &
    !$acc eflx(:), &
    !$acc qflx(:), &
    !$acc eflx_scale(:), &
    !$acc qflx_scale(:), &
    !$acc eflx_err(:), &
    !$acc qflx_err(:), &
    !$acc local_secp1(:), &
    !$acc wind_speed0(:), &
    !$acc wind_speed_adj(:), &
    !$acc tau(:), &
    !$acc tau_diff(:), &
    !$acc prev_tau(:), &
    !$acc prev_tau_diff(:), &
    !$acc lnd_to_urban_filter(:), &
    !$acc converged_landunits(:), &
    !$acc num_copyl, &
    !$acc zeta, &
    !$acc fwet_roof, &
    !$acc fwet_road_imperv)

         
       begl =   bounds%begl 
       endl =   bounds%endl
       begc =   bounds%begc
       endc =   bounds%endc 
       ! Define fields that appear on the restart file for non-urban landunits
       !$acc parallel loop independent gang vector default(present) 
       do fl = 1,num_nourbanl
         l = filter_nourbanl(fl)
         taf(l) = spval
         qaf(l) = spval
       end do

      ! Get time step

      nstep = nstep_mod
      dtime = dtime_mod
      year = year_curr
      month = mon_curr
      day = day_curr
      secs = secs_curr
      lnd_to_urban_filter(:) = -9999
      ! Compute canyontop wind using Masson (2000)
      !$acc parallel loop independent gang vector default(present) 
      do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         lnd_to_urban_filter(l) = fl 
         g = lun_pp%gridcell(l)
         t = lun_pp%topounit(l)

         local_secp1(fl)        = secs + nint((grc_pp%londeg(g)/degpsec)/dtime)*dtime
         local_secp1(fl)        = mod(local_secp1(fl),isecspday)

         ! Error checks

         if (ht_roof(l) - z_d_town(l) <= z_0_town(l)) then
            write (iulog,*) 'aerodynamic parameter error in UrbanFluxes'
            write (iulog,*) 'h_r - z_d <= z_0'
            write (iulog,*) 'ht_roof, z_d_town, z_0_town: ', ht_roof(l), z_d_town(l), &
                 z_0_town(l)
            write (iulog,*) 'elm model is stopping'
            call endrun(decomp_index=l, elmlevel=namel, msg=errmsg(__FILE__, __LINE__))
         end if
         if (forc_hgt_u_patch(lun_pp%pfti(l)) - z_d_town(l) <= z_0_town(l)) then
            write (iulog,*) 'aerodynamic parameter error in UrbanFluxes'
            write (iulog,*) 'h_u - z_d <= z_0'
            write (iulog,*) 'forc_hgt_u_patch, z_d_town, z_0_town: ', forc_hgt_u_patch(lun_pp%pfti(l)), z_d_town(l), &
                 z_0_town(l)
            write (iulog,*) 'elm model is stopping'
            call endrun(decomp_index=l, elmlevel=namel, msg=errmsg(__FILE__, __LINE__))
         end if
         ! Initialize winds for iteration.
         if (implicit_stress) then
            wind_speed0(fl) = max(0.01_r8, hypot(forc_u(t), forc_v(t)))
            wind_speed_adj(fl) = wind_speed0(fl)
            ur(fl) = max(1.0_r8, wind_speed_adj(fl) + ugust(t))

            prev_tau(fl) = tau_est(t)
         else
            ur(fl) = max(1.0_r8,sqrt(forc_u(t)*forc_u(t)+forc_v(t)*forc_v(t)) + ugust(t))
         end if
         tau_diff(fl) = 1.e100_r8

      end do

      ! Compute fluxes - Follows elm approach for bare soils (Oleson et al 2004)

      !$acc parallel loop independent gang vector default(present)
      do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         t = lun_pp%topounit(l)
         g = lun_pp%gridcell(l)

         thm_g(fl) = forc_t(t) + lapse_rate*forc_hgt_t_patch(lun_pp%pfti(l))
         thv_g(fl) = forc_th(t)*(1._r8+0.61_r8*forc_q(t))
         dth(fl)   = thm_g(fl)-taf(l)
         dqh(fl)   = forc_q(t)-qaf(l)
         dthv     = dth(fl)*(1._r8+0.61_r8*forc_q(t))+0.61_r8*forc_th(t)*dqh(fl)
         zldis(fl) = forc_hgt_u_patch(lun_pp%pfti(l)) - z_d_town(l)

         ! Initialize Monin-Obukhov length and wind speed including convective velocity

         call MoninObukIni(ur(fl), thv_g(fl), dthv, zldis(fl), z_0_town(l), um(fl), obu(fl))
         ! Initialize conductances
         wtus_roof(fl)              = 0._r8
         wtus_road_perv(fl)         = 0._r8
         wtus_road_imperv(fl)       = 0._r8
         wtus_sunwall(fl)           = 0._r8
         wtus_shadewall(fl)         = 0._r8
         wtuq_roof(fl)              = 0._r8
         wtuq_road_perv(fl)         = 0._r8
         wtuq_road_imperv(fl)       = 0._r8
         wtuq_sunwall(fl)           = 0._r8
         wtuq_shadewall(fl)         = 0._r8
         wtus_roof_unscl(fl)        = 0._r8
         wtus_road_perv_unscl(fl)   = 0._r8
         wtus_road_imperv_unscl(fl) = 0._r8
         wtus_sunwall_unscl(fl)     = 0._r8
         wtus_shadewall_unscl(fl)   = 0._r8
         wtuq_roof_unscl(fl)        = 0._r8
         wtuq_road_perv_unscl(fl)   = 0._r8
         wtuq_road_imperv_unscl(fl) = 0._r8
         wtuq_sunwall_unscl(fl)     = 0._r8
         wtuq_shadewall_unscl(fl)   = 0._r8
      end do

      

      ! Start stability iteration
      num_copyl = num_urbanl
      num_copyc = num_urbanc
      ! filter_copyl(1:num_urbanl) = filter_urbanl(1:num_urbanl)
      ! filter_copyc(1:num_urbanc) = filter_urbanc(1:num_urbanc)

      if (implicit_stress) then
         loopmax = itmax
      else
         loopmax = itmin
      end if
      ! converged_cols(begc:endc) = .false.
      converged_landunits(begl:endl) = .false.
      col_to_urban_filter(:) = -9999

      ITERATION: do iter = 1, loopmax

         ! Get friction velocity, relation for potential
         ! temperature and humidity profiles of surface boundary layer.

         if (num_copyl > 0) then
            call FrictionVelocity_loops(begl, endl, &
                 num_urbanl, filter_urbanl, &
                 z_d_town(begl:endl), z_0_town(begl:endl), z_0_town(begl:endl), z_0_town(begl:endl), &
                 obu(1:num_urbanl), iter, ur(1:num_urbanl), um(1:num_urbanl), ustar(1:num_urbanl), &
                 temp1(1:num_urbanl), temp2(1:num_urbanl), temp12m(1:num_urbanl), temp22m(1:num_urbanl), fm(1:num_urbanl), &
                 frictionvel_vars, converged_landunits(begl:endl),landunit_index=.true.)
         end if

         !$acc parallel loop independent gang vector default(present)
         do fl = 1, num_urbanl
            l = filter_urbanl(fl)
            t = lun_pp%topounit(l)
            g = lun_pp%gridcell(l)
            if(converged_landunits(l)) cycle

            ! Determine aerodynamic resistance to fluxes from urban canopy air to
            ! atmosphere
            ramu(fl) = 1._r8/(ustar(fl)*ustar(fl)/um(fl))
            rahu(fl) = 1._r8/(temp1(fl)*ustar(fl))
            rawu(fl) = 1._r8/(temp2(fl)*ustar(fl))

            ! Calculate magnitude of stress and update wind speed.
            if (implicit_stress) then
               tau(fl) = forc_rho(t)*wind_speed_adj(fl)/ramu(fl)
               call shr_flux_update_stress(wind_speed0(l), wsresp(t), tau_est(t), &
                    tau(l), prev_tau(l), tau_diff(l), prev_tau_diff(l), &
                    wind_speed_adj(l))
               ur(fl) = max(1.0_r8, wind_speed_adj(fl) + ugust(t))
            end if

            ! Canyon top wind
            ! If the wind does not change in this loop (explicit stress), then
            ! we only need to calculate this on the first iteration.
            if (implicit_stress .or. iter == 1) then
               canyontop_wind(fl) = ur(fl) * &
                    log( (ht_roof(l)-z_d_town(l)) / z_0_town(l) ) / &
                    log( (forc_hgt_u_patch(lun_pp%pfti(l))-z_d_town(l)) / z_0_town(l) )

               ! U component of canyon wind 

               if (canyon_hwr(l) < 0.5_r8) then  ! isolated roughness flow
                  canyon_u_wind(fl) = canyontop_wind(fl) * exp( -0.5_r8*canyon_hwr(l)* &
                       (1._r8-(wind_hgt_canyon(l)/ht_roof(l))) )
               else if (canyon_hwr(l) < 1.0_r8) then ! wake interference flow
                  canyon_u_wind(fl) = canyontop_wind(fl) * (1._r8+2._r8*(2._r8/rpi - 1._r8)* &
                       (ht_roof(l)/(ht_roof(l)/canyon_hwr(l)) - 0.5_r8)) * &
                       exp(-0.5_r8*canyon_hwr(l)*(1._r8-(wind_hgt_canyon(l)/ht_roof(l))))
               else  ! skimming flow
                  canyon_u_wind(fl) = canyontop_wind(fl) * (2._r8/rpi) * &
                       exp(-0.5_r8*canyon_hwr(l)*(1._r8-(wind_hgt_canyon(l)/ht_roof(l))))
               end if
            end if

            ! Determine magnitude of canyon wind by using horizontal wind determined
            ! previously and vertical wind from friction velocity (Masson 2000)

            canyon_wind(fl) = sqrt(canyon_u_wind(fl)**2._r8 + ustar(fl)**2._r8)

            ! Determine canyon_resistance (currently this single resistance determines the
            ! resistance from urban surfaces (roof, pervious and impervious road, sunlit and
            ! shaded walls) to urban canopy air, since it is only dependent on wind speed
            ! Also from Masson 2000.

            canyon_resistance(fl) = cpair * forc_rho(t) / (11.8_r8 + 4.2_r8*canyon_wind(fl))

         end do

         ! This is the first term in the equation solutions for urban canopy air temperature
         ! and specific humidity (numerator) and is a landunit quantity
         !$acc parallel loop independent gang vector default(present)
         do fl = 1, num_urbanl
            l = filter_urbanl(fl)
            t = lun_pp%topounit(l)
            g = lun_pp%gridcell(l)
            if(converged_landunits(l)) cycle

            taf_numer(fl) = thm_g(fl)/rahu(fl)
            taf_denom(fl) = 1._r8/rahu(fl)
            taf_numer_test(fl) = thm_g(fl)/rahu(fl)
            taf_denom_test(fl) = 1._r8/rahu(fl)
            qaf_numer(fl) = forc_q(t)/rawu(fl)
            qaf_denom(fl) = 1._r8/rawu(fl)

            ! First term needed for derivative of heat fluxes
            wtas(fl) = 1._r8/rahu(fl)
            wtaq(fl) = 1._r8/rawu(fl)

         end do


         ! Gather other terms for other urban columns for numerator and denominator of
         ! equations for urban canopy air temperature and specific humidity
         do fc = 1, num_urbanc
            c = filter_urbanc(fc)
            l = col_pp%landunit(c)
            fl = lnd_to_urban_filter(l)
            col_to_urban_filter(c) = fc 
            if(converged_landunits(l)) cycle

            if (ctype(c) == icol_roof) then

               ! scaled sensible heat conductance
               wtus(fc) = wtlunit_roof(l)/canyon_resistance(fl)
               wtus_roof(fl) = wtus(fc)
               ! unscaled sensible heat conductance
               wtus_roof_unscl(fl) = 1._r8/canyon_resistance(fl)

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
               wtuq(fc) = fwet_roof*(wtlunit_roof(l)/canyon_resistance(fl))
               wtuq_roof(fl) = wtuq(fc)
               ! unscaled latent heat conductance
               wtuq_roof_unscl(fl) = fwet_roof*(1._r8/canyon_resistance(fl))

               ! wasteheat from heating/cooling
               if (trim(urban_hac) == urban_wasteheat_on) then
                  eflx_wasteheat_roof(fl) = ac_wasteheat_factor * eflx_urban_ac(c) + &
                       ht_wasteheat_factor * eflx_urban_heat(c)
               else
                  eflx_wasteheat_roof(fl) = 0._r8
               end if

               ! If air conditioning on, always replace heat removed with heat into canyon
               if (trim(urban_hac) == urban_hac_on .or. trim(urban_hac) == urban_wasteheat_on) then
                  eflx_heat_from_ac_roof(fl) = abs(eflx_urban_ac(c))
               else
                  eflx_heat_from_ac_roof(fl) = 0._r8
               end if

            else if (ctype(c) == icol_road_perv) then

               ! scaled sensible heat conductance
               wtus(fc) = wtroad_perv(l)*(1._r8-wtlunit_roof(l))/canyon_resistance(fl)
               wtus_road_perv(fl) = wtus(fc)
               ! unscaled sensible heat conductance
               wtus_road_perv_unscl(fl) = 1._r8/canyon_resistance(fl)

               ! scaled latent heat conductance
               wtuq(fc) = wtroad_perv(l)*(1._r8-wtlunit_roof(l))/canyon_resistance(fl)
               wtuq_road_perv(fl) = wtuq(fc)
               ! unscaled latent heat conductance
               wtuq_road_perv_unscl(fl) = 1._r8/canyon_resistance(fl)

               if (use_vsfm) then
                  if (qaf(l) < qg(c)) then
                     if (do_soilevap_beta()) then
                        wtuq_road_perv(fl)       = soilbeta(c)*wtuq_road_perv(fl)
                        wtuq_road_perv_unscl(fl) = soilbeta(c)*wtuq_road_perv_unscl(fl)
                     endif
                  endif
               endif
            else if (ctype(c) == icol_road_imperv) then

               ! scaled sensible heat conductance
               wtus(fc) = (1._r8-wtroad_perv(l))*(1._r8-wtlunit_roof(l))/canyon_resistance(fl)
               wtus_road_imperv(fl) = wtus(fc)
               ! unscaled sensible heat conductance
               wtus_road_imperv_unscl(fl) = 1._r8/canyon_resistance(fl)

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
               wtuq(fc) = fwet_road_imperv*(1._r8-wtroad_perv(l))*(1._r8-wtlunit_roof(l))/canyon_resistance(fl)
               wtuq_road_imperv(fl) = wtuq(fc)
               ! unscaled latent heat conductance
               wtuq_road_imperv_unscl(fl) = fwet_road_imperv*(1._r8/canyon_resistance(fl))

            else if (ctype(c) == icol_sunwall) then

               ! scaled sensible heat conductance
               wtus(fc) = canyon_hwr(l)*(1._r8-wtlunit_roof(l))/canyon_resistance(fl)
               wtus_sunwall(fl) = wtus(fc)
               ! unscaled sensible heat conductance
               wtus_sunwall_unscl(fl) = 1._r8/canyon_resistance(fl)

               ! scaled latent heat conductance
               wtuq(fc) = 0._r8
               wtuq_sunwall(fl) = wtuq(fc)
               ! unscaled latent heat conductance
               wtuq_sunwall_unscl(fl) = 0._r8

               ! wasteheat from heating/cooling
               if (trim(urban_hac) == urban_wasteheat_on) then
                  eflx_wasteheat_sunwall(fl) = ac_wasteheat_factor * eflx_urban_ac(c) + &
                       ht_wasteheat_factor * eflx_urban_heat(c)
               else
                  eflx_wasteheat_sunwall(fl) = 0._r8
               end if

               ! If air conditioning on, always replace heat removed with heat into canyon
               if (trim(urban_hac) == urban_hac_on .or. trim(urban_hac) == urban_wasteheat_on) then
                  eflx_heat_from_ac_sunwall(fl) = abs(eflx_urban_ac(c))
               else
                  eflx_heat_from_ac_sunwall(fl) = 0._r8
               end if

            else if (ctype(c) == icol_shadewall) then

               ! scaled sensible heat conductance
               wtus(fc) = canyon_hwr(l)*(1._r8-wtlunit_roof(l))/canyon_resistance(fl)
               wtus_shadewall(fl) = wtus(fc)
               ! unscaled sensible heat conductance
               wtus_shadewall_unscl(fl) = 1._r8/canyon_resistance(fl)

               ! scaled latent heat conductance
               wtuq(fc) = 0._r8
               wtuq_shadewall(fl) = wtuq(fc)
               ! unscaled latent heat conductance
               wtuq_shadewall_unscl(fl) = 0._r8

               ! wasteheat from heating/cooling
               if (trim(urban_hac) == urban_wasteheat_on) then
                  eflx_wasteheat_shadewall(fl) = ac_wasteheat_factor * eflx_urban_ac(c) + &
                       ht_wasteheat_factor * eflx_urban_heat(c)
               else
                  eflx_wasteheat_shadewall(fl) = 0._r8
               end if

               ! If air conditioning on, always replace heat removed with heat into canyon
               if (trim(urban_hac) == urban_hac_on .or. trim(urban_hac) == urban_wasteheat_on) then
                  eflx_heat_from_ac_shadewall(fl) = abs(eflx_urban_ac(c))
               else
                  eflx_heat_from_ac_shadewall(fl) = 0._r8
               end if
            else
               write(iulog,*) 'c, ctype, pi = ', c, ctype(c), pi
               write(iulog,*) 'Column indices for: shadewall, sunwall, road_imperv, road_perv, roof: '
               write(iulog,*) icol_shadewall, icol_sunwall, icol_road_imperv, icol_road_perv, icol_roof
               call endrun(decomp_index=l, elmlevel=namel, msg="ERROR, ctype out of range"//errmsg(__FILE__, __LINE__))
            end if

            taf_numer(fl) = taf_numer(fl) + t_grnd(c)*wtus(fc)
            taf_denom(fl) = taf_denom(fl) + wtus(fc)
            qaf_numer(fl) = qaf_numer(fl) + qg(c)*wtuq(fc)
            qaf_denom(fl) = qaf_denom(fl) + wtuq(fc)

         end do
        ! !$acc parallel loop independent gang worker default(present) private(sum_denom,sum_numer)
        ! do fl = 1, num_urbanl
        !    sum_denom = 0._r8 
        !    sum_numer = 0._r8 
        !    l = filter_urbanl(fl) 
        !    if(converged_landunits(l)) cycle 
        !    !$acc loop vector reduction(+:sum_denom,sum_numer)
        !    do c = lun_pp%coli(l), lun_pp%colf(l) 
        !       if(col_pp%active(c)) then 
        !          fc = col_to_urban_filter(c)
        !          sum_denom = sum_denom + wtus(fc)
        !          sum_numer = sum_numer + t_grnd(c)*wtus(fc)
        !       end if 
        !    end do 
        !    taf_denom_test(fl) = taf_denom_test(fl) + sum_denom
        !    taf_numer(fl) = taf_numer(fl) + sum_numer  
        ! end do 
         
         ! !$acc parallel loop independent gang worker default(present) private(sum_denom,sum_numer)
         ! do fl = 1, num_urbanl
         !    sum_denom = 0._r8 
         !    sum_numer = 0._r8 
         !    l = filter_urbanl(fl) 
         !    !$acc loop vector reduction(+:sum_denom,sum_numer)
         !    do c = lun_pp%coli(l), lun_pp%colf(l) 
         !       if(col_pp%active(c)) then 
         !          fc = col_to_urban_filter(c)
         !          sum_denom = sum_denom + wtuq(fc)
         !          sum_numer = sum_numer + qg(c)*wtuq(fc)
         !       end if 
         !    end do 
         !    qaf_denom(fl) = qaf_denom(fl) + sum_denom
         !    qaf_numer(fl) = qaf_numer(fl) + sum_numer 
         ! end do 

         ! Calculate new urban canopy air temperature and specific humidity

         !$acc parallel loop independent gang vector default(present)
         do fl = 1, num_urbanl
            l = filter_urbanl(fl)
            g = lun_pp%gridcell(l)
            if(converged_landunits(l)) cycle
            ! Total waste heat and heat from AC is sum of heat for walls and roofs
            ! accounting for different surface areas
            eflx_wasteheat(l) = wtlunit_roof(l)*eflx_wasteheat_roof(fl) + &
                 (1._r8-wtlunit_roof(l))*(canyon_hwr(l)*(eflx_wasteheat_sunwall(fl) + &
                 eflx_wasteheat_shadewall(fl)))

            ! Limit wasteheat to ensure that we don't get any unrealistically strong
            ! positive feedbacks due to AC in a warmer climate
            eflx_wasteheat(l) = min(eflx_wasteheat(l),wasteheat_limit)

            eflx_heat_from_ac(l) = wtlunit_roof(l)*eflx_heat_from_ac_roof(fl) + &
                 (1._r8-wtlunit_roof(l))*(canyon_hwr(l)*(eflx_heat_from_ac_sunwall(fl) + &
                 eflx_heat_from_ac_shadewall(fl)))

            ! Calculate traffic heat flux
            ! Only comes from impervious road
            eflx_traffic(l) = (1._r8-wtlunit_roof(l))*(1._r8-wtroad_perv(l))* &
                 eflx_traffic_factor(l)

            taf(l) = taf_numer(fl)/taf_denom(fl)
            qaf(l) = qaf_numer(fl)/qaf_denom(fl)

            wts_sum(fl) = wtas(fl) + wtus_roof(fl) + wtus_road_perv(fl) + &
                 wtus_road_imperv(fl) + wtus_sunwall(fl) + wtus_shadewall(fl)

            wtq_sum(fl) = wtaq(fl) + wtuq_roof(fl) + wtuq_road_perv(fl) + &
                 wtuq_road_imperv(fl) + wtuq_sunwall(fl) + wtuq_shadewall(fl)

         end do

         ! This section of code is not required if niters = 1
         ! Determine stability using new taf and qaf
         ! TODO: Some of these constants replicate what is in FrictionVelocity and BareGround fluxes should consildate. EBK
         do fl = 1, num_urbanl
            l = filter_urbanl(fl)
            t = lun_pp%topounit(l)
            g = lun_pp%gridcell(l)
            if(converged_landunits(l)) cycle

            dth(fl) = thm_g(fl)-taf(l)
            dqh(fl) = forc_q(t)-qaf(l)
            tstar = temp1(fl)*dth(fl)
            qstar = temp2(fl)*dqh(fl)
            thvstar = tstar*(1._r8+0.61_r8*forc_q(t)) + 0.61_r8*forc_th(t)*qstar
            zeta = zldis(fl)*vkc*grav*thvstar/(ustar(fl)**2*thv_g(fl))

            if (zeta >= 0._r8) then                   !stable
               zeta = min(2._r8,max(zeta,0.01_r8))
               um(fl) = max(ur(fl),0.1_r8)
            else                                      !unstable
               zeta = max(-100._r8,min(zeta,-0.01_r8))
               wc = beta*(-grav*ustar(fl)*thvstar*zii/thv_g(fl))**0.333_r8
               um(fl) = sqrt(ur(fl)*ur(fl) + wc*wc)
            end if

            obu(fl) = zldis(fl)/zeta
         end do

         ! Test for convergence
         iter_final = iter
         if (iter >= itmin) then
            num_copyl_old = num_copyl
            num_copyl = 0
            do fl = 1, num_urbanl
               l = filter_urbanl(fl)
               if (.not. (abs(tau_diff(fl)) < dtaumin)) then
                  num_copyl = num_copyl + 1
                  ! filter_copyl(num_copyl) = l
               else
                  converged_landunits(l) = .true. 
               end if
            end do
            if (num_copyl == 0) then
               exit ITERATION
            end if
         end if

      end do ITERATION ! end iteration

      ! Determine fluxes from canyon surfaces

      ! the following initializations are needed to ensure that the values are 0 over non-
      ! active urban Patches
      eflx_sh_grnd_scale(bounds%begp : bounds%endp) = 0._r8
      qflx_evap_soi_scale(bounds%begp : bounds%endp) = 0._r8

      do f = 1, num_urbanp

         p = filter_urbanp(f)
         c = veg_pp%column(p)
         g = veg_pp%gridcell(p)
         t = veg_pp%topounit(p)
         l = veg_pp%landunit(p)
         fl = lnd_to_urban_filter(l)
         fc = col_to_urban_filter(c)
         ram1(p) = ramu(fl)  !pass value to global variable

         ! Upward and downward canopy longwave are zero

         ulrad(p)  = 0._r8
         dlrad(p)  = 0._r8

         ! Derivative of sensible and latent heat fluxes with respect to
         ! ground temperature

         if (ctype(c) == icol_roof) then
            cgrnds(p) = forc_rho(t) * cpair * (wtas(fl) + wtus_road_perv(fl) +  &
                 wtus_road_imperv(fl) + wtus_sunwall(fl) + wtus_shadewall(fl)) * &
                 (wtus_roof_unscl(fl)/wts_sum(fl))
            cgrndl(p) = forc_rho(t) * (wtaq(fl) + wtuq_road_perv(fl) +  &
                 wtuq_road_imperv(fl) + wtuq_sunwall(fl) + wtuq_shadewall(fl)) * &
                 (wtuq_roof_unscl(fl)/wtq_sum(fl))*dqgdT(c)
         else if (ctype(c) == icol_road_perv) then
            cgrnds(p) = forc_rho(t) * cpair * (wtas(fl) + wtus_roof(fl) +  &
                 wtus_road_imperv(fl) + wtus_sunwall(fl) + wtus_shadewall(fl)) * &
                 (wtus_road_perv_unscl(fl)/wts_sum(fl))
            cgrndl(p) = forc_rho(t) * (wtaq(fl) + wtuq_roof(fl) +  &
                 wtuq_road_imperv(fl) + wtuq_sunwall(fl) + wtuq_shadewall(fl)) * &
                 (wtuq_road_perv_unscl(fl)/wtq_sum(fl))*dqgdT(c)
         else if (ctype(c) == icol_road_imperv) then
            cgrnds(p) = forc_rho(t) * cpair * (wtas(fl) + wtus_roof(fl) +  &
                 wtus_road_perv(fl) + wtus_sunwall(fl) + wtus_shadewall(fl)) * &
                 (wtus_road_imperv_unscl(fl)/wts_sum(fl))
            cgrndl(p) = forc_rho(t) * (wtaq(fl) + wtuq_roof(fl) +  &
                 wtuq_road_perv(fl) + wtuq_sunwall(fl) + wtuq_shadewall(fl)) * &
                 (wtuq_road_imperv_unscl(fl)/wtq_sum(fl))*dqgdT(c)
         else if (ctype(c) == icol_sunwall) then
            cgrnds(p) = forc_rho(t) * cpair * (wtas(fl) + wtus_roof(fl) +  &
                 wtus_road_perv(fl) + wtus_road_imperv(fl) + wtus_shadewall(fl)) * &
                 (wtus_sunwall_unscl(fl)/wts_sum(fl))
            cgrndl(p) = 0._r8
         else if (ctype(c) == icol_shadewall) then
            cgrnds(p) = forc_rho(t) * cpair * (wtas(fl) + wtus_roof(fl) +  &
                 wtus_road_perv(fl) + wtus_road_imperv(fl) + wtus_sunwall(fl)) * &
                 (wtus_shadewall_unscl(fl)/wts_sum(fl))
            cgrndl(p) = 0._r8
         end if
         cgrnd(p)  = cgrnds(p) + cgrndl(p)*htvp(c)

         ! Surface fluxes of momentum, sensible and latent heat

         taux(p) = -forc_rho(t)*forc_u(t)/ramu(fl)
         tauy(p) = -forc_rho(t)*forc_v(t)/ramu(fl)
         if (implicit_stress) then
            taux(p) = taux(p) * (wind_speed_adj(fl) / wind_speed0(fl))
            tauy(p) = tauy(p) * (wind_speed_adj(fl) / wind_speed0(fl))
         end if

         ! Use new canopy air temperature
         dth(fl) = taf(l) - t_grnd(c)

         if (ctype(c) == icol_roof) then
            eflx_sh_grnd(p)  = -forc_rho(t)*cpair*wtus_roof_unscl(fl)*dth(fl)
            eflx_sh_snow(p)  = 0._r8
            eflx_sh_soil(p)  = 0._r8
            eflx_sh_h2osfc(p)= 0._r8
         else if (ctype(c) == icol_road_perv) then
            eflx_sh_grnd(p)  = -forc_rho(t)*cpair*wtus_road_perv_unscl(fl)*dth(fl)
            eflx_sh_snow(p)  = 0._r8
            eflx_sh_soil(p)  = 0._r8
            eflx_sh_h2osfc(p)= 0._r8
         else if (ctype(c) == icol_road_imperv) then
            eflx_sh_grnd(p)  = -forc_rho(t)*cpair*wtus_road_imperv_unscl(fl)*dth(fl)
            eflx_sh_snow(p)  = 0._r8
            eflx_sh_soil(p)  = 0._r8
            eflx_sh_h2osfc(p)= 0._r8
         else if (ctype(c) == icol_sunwall) then
            eflx_sh_grnd(p)  = -forc_rho(t)*cpair*wtus_sunwall_unscl(fl)*dth(fl)
            eflx_sh_snow(p)  = 0._r8
            eflx_sh_soil(p)  = 0._r8
            eflx_sh_h2osfc(p)= 0._r8
         else if (ctype(c) == icol_shadewall) then
            eflx_sh_grnd(p)  = -forc_rho(t)*cpair*wtus_shadewall_unscl(fl)*dth(fl)
            eflx_sh_snow(p)  = 0._r8
            eflx_sh_soil(p)  = 0._r8
            eflx_sh_h2osfc(p)= 0._r8
         end if

         eflx_sh_tot(p)   = eflx_sh_grnd(p)
         eflx_sh_tot_u(p) = eflx_sh_tot(p)

         dqh(fl) = qaf(l) - qg(c)

         if (ctype(c) == icol_roof) then
            qflx_evap_soi(p) = -forc_rho(t)*wtuq_roof_unscl(fl)*dqh(fl)
         else if (ctype(c) == icol_road_perv) then
            ! Evaporation assigned to soil term if dew or snow
            ! or if no liquid water available in soil column
            if (dqh(fl) > 0._r8 .or. frac_sno(c) > 0._r8 .or. soilalpha_u(c) <= 0._r8) then
               qflx_evap_soi(p) = -forc_rho(t)*wtuq_road_perv_unscl(fl)*dqh(fl)
               qflx_tran_veg(p) = 0._r8
               ! Otherwise, evaporation assigned to transpiration term
            else
               qflx_evap_soi(p) = 0._r8
               qflx_tran_veg(p) = -forc_rho(t)*wtuq_road_perv_unscl(fl)*dqh(fl)
            end if
            qflx_evap_veg(p) = qflx_tran_veg(p)
         else if (ctype(c) == icol_road_imperv) then
            qflx_evap_soi(p) = -forc_rho(t)*wtuq_road_imperv_unscl(fl)*dqh(fl)
         else if (ctype(c) == icol_sunwall) then
            qflx_evap_soi(p) = 0._r8
         else if (ctype(c) == icol_shadewall) then
            qflx_evap_soi(p) = 0._r8
         end if

         ! SCALED sensible and latent heat flux for error check
         eflx_sh_grnd_scale(p)  = -forc_rho(t)*cpair*wtus(fc)*dth(fl)
         qflx_evap_soi_scale(p) = -forc_rho(t)*wtuq(fc)*dqh(fl)

      end do

      ! Check to see that total sensible and latent heat equal the sum of
      ! the scaled heat fluxes above
      !$acc parallel loop independent gang vector default(present)
      do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         t = lun_pp%topounit(l)
         g = lun_pp%gridcell(l)
         eflx(fl)       = -(forc_rho(t)*cpair/rahu(fl))*(thm_g(fl) - taf(l))
         qflx(fl)       = -(forc_rho(t)/rawu(fl))*(forc_q(t) - qaf(l))
         eflx_scale(fl) = sum(eflx_sh_grnd_scale(lun_pp%pfti(l):lun_pp%pftf(l)))
         qflx_scale(fl) = sum(qflx_evap_soi_scale(lun_pp%pfti(l):lun_pp%pftf(l)))
         eflx_err(fl)   = eflx_scale(fl) - eflx(fl)
         qflx_err(fl)   = qflx_scale(fl) - qflx(fl)
      end do

      found = .false.
      !$acc parallel loop independent gang vector default(present)
      do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         if (abs(eflx_err(fl)) > 0.01_r8) then
            found = .true.
            indexl = fl
            exit
         end if
      end do

      if ( found ) then
         write(iulog,*)'WARNING:  Total sensible heat does not equal sum of scaled heat fluxes for urban columns ',&
              ' nstep = ',nstep,' indexl= ',indexl,' eflx_err= ',eflx_err(indexl)
         if (abs(eflx_err(indexl)) > .01_r8) then
            l = filter_urbanl(indexl)
            write(iulog,*)'elm model is stopping - error is greater than .01 W/m**2'
            write(iulog,*)'eflx_scale    = ',eflx_scale(indexl)
            write(iulog,*)'eflx_sh_grnd_scale: ',eflx_sh_grnd_scale(lun_pp%pfti(l):lun_pp%pftf(l))
            write(iulog,*)'eflx          = ',eflx(indexl)
            ! test code, PET
            write(iulog,*)'tbot          = ',forc_t(lun_pp%topounit(l))
            call endrun(decomp_index=indexl, elmlevel=namel, msg=errmsg(__FILE__, __LINE__))
         end if
      end if

      found = .false.
      !$acc parallel loop independent gang vector default(present)
      do fl = 1, num_urbanl
         l = filter_urbanl(fl)
         ! 4.e-9 kg/m**2/s = 0.01 W/m**2
         if (abs(qflx_err(fl)) > 4.e-9_r8) then
            found = .true.
            indexl = l
            exit
         end if
      end do

      if ( found ) then
         write(iulog,*)'WARNING:  Total water vapor flux does not equal sum of scaled water vapor fluxes for urban columns ',&
              ' nstep = ',nstep,' indexl= ',indexl,' qflx_err= ',qflx_err(indexl)
         if (abs(qflx_err(indexl)) > 4.e-9_r8) then
            write(iulog,*)'elm model is stopping - error is greater than 4.e-9 kg/m**2/s'
            write(iulog,*)'qflx_scale    = ',qflx_scale(indexl)
            write(iulog,*)'qflx          = ',qflx(indexl)
            call endrun(decomp_index=indexl, elmlevel=namel, msg=errmsg(__FILE__, __LINE__))
         end if
      end if

      ! Check for convergence of stress.
      if (implicit_stress) then
         !$acc parallel loop independent gang vector default(present)
         do fl = 1, num_urbanl
            l = filter_urbanl(fl)
            if (abs(tau_diff(fl)) > dtaumin) then
               if (nstep > 0) then ! Suppress common warnings on the first time step.
                  write(iulog,*)'WARNING: Stress did not converge for urban columns ',&
                       ' nstep = ',nstep,' indexl= ',l,' prev_tau_diff= ',prev_tau_diff(l),&
                       ' tau_diff= ',tau_diff(fl),' tau= ',tau(fl),&
                       ' wind_speed_adj= ',wind_speed_adj(fl),' iter_final= ',iter_final
               end if
            end if
         end do
      end if

      ! Gather terms required to determine internal building temperature

      !$acc parallel loop independent gang vector default(present)
      do fc = 1,num_urbanc
         c = filter_urbanc(fc)
         l = col_pp%landunit(c)
         fl = lnd_to_urban_filter(l)
         if (ctype(c) == icol_roof) then
            t_roof_innerl(fl) = t_soisno(c,nlevurb)
         else if (ctype(c) == icol_sunwall) then
            t_sunwall_innerl(fl) = t_soisno(c,nlevurb)
         else if (ctype(c) == icol_shadewall) then
            t_shadewall_innerl(fl) = t_soisno(c,nlevurb)
         end if

      end do

      ! Calculate internal building temperature
      !$acc parallel loop independent gang vector default(present)
      do fl = 1, num_urbanl
         l = filter_urbanl(fl)

         lngth_roof = (ht_roof(l)/canyon_hwr(l))*wtlunit_roof(l)/(1._r8-wtlunit_roof(l))
         t_building(l) = (ht_roof(l)*(t_shadewall_innerl(fl) + t_sunwall_innerl(fl)) &
              +lngth_roof*t_roof_innerl(fl))/(2._r8*ht_roof(l)+lngth_roof)
      end do

      ! No roots for urban except for pervious road

      !$acc parallel loop independent gang vector default(present) collapse(2) 
      do j = 1, nlevgrnd
         do f = 1, num_urbanp
            p = filter_urbanp(f)
            c = veg_pp%column(p)
            if (ctype(c) == icol_road_perv) then
               rootr(p,j) = rootr_road_perv(c,j)
            else

               rootr(p,j) = 0._r8
            end if
         end do
      end do

      !$acc parallel loop independent gang vector default(present)
      do f = 1, num_urbanp

         p = filter_urbanp(f)
         c = veg_pp%column(p)
         g = veg_pp%gridcell(p)
         t = veg_pp%topounit(p)
         l = veg_pp%landunit(p)

         ! Use urban canopy air temperature and specific humidity to represent
         ! 2-m temperature and humidity

         t_ref2m(p) = taf(l)
         q_ref2m(p) = qaf(l)
         t_ref2m_u(p) = taf(l)

         ! 2 m height relative humidity

         call QSat(t_ref2m(p), forc_pbot(t), e_ref2m, de2mdT, qsat_ref2m, dqsat2mdT)
         rh_ref2m(p) = min(100._r8, q_ref2m(p) / qsat_ref2m * 100._r8)
         rh_ref2m_u(p) = rh_ref2m(p)

         ! Variables needed by history tape
         t_veg(p) = forc_t(t)

      end do

    !$acc exit data delete(&
    !$acc canyontop_wind(:), &
    !$acc canyon_u_wind(:), &
    !$acc canyon_wind(:), &
    !$acc canyon_resistance(:), &
    !$acc ur(:), &
    !$acc ustar(:), &
    !$acc ramu(:), &
    !$acc rahu(:), &
    !$acc rawu(:), &
    !$acc temp1(:), &
    !$acc temp12m(:), &
    !$acc temp2(:), &
    !$acc temp22m(:), &
    !$acc thm_g(:), &
    !$acc thv_g(:), &
    !$acc dth(:), &
    !$acc dqh(:), &
    !$acc zldis(:), &
    !$acc um(:), &
    !$acc obu(:), &
    !$acc taf_numer(:), &
    !$acc taf_denom(:), &
    !$acc qaf_numer(:), &
    !$acc qaf_denom(:), &
    !$acc wtas(:), &
    !$acc wtaq(:), &
    !$acc wts_sum(:), &
    !$acc wtq_sum(:), &
    !$acc fm(:), &
    !$acc wtus(:), &
    !$acc wtuq(:), &
    !$acc wtus_roof(:), &
    !$acc wtuq_roof(:), &
    !$acc wtus_road_perv(:), &
    !$acc wtuq_road_perv(:), &
    !$acc wtus_road_imperv(:), &
    !$acc wtuq_road_imperv(:), &
    !$acc wtus_sunwall(:), &
    !$acc wtuq_sunwall(:), &
    !$acc wtus_shadewall(:), &
    !$acc wtuq_shadewall(:), &
    !$acc wtus_roof_unscl(:), &
    !$acc wtuq_roof_unscl(:), &
    !$acc wtus_road_perv_unscl(:), &
    !$acc wtuq_road_perv_unscl(:), &
    !$acc wtus_road_imperv_unscl(:), &
    !$acc wtuq_road_imperv_unscl(:), &
    !$acc wtus_sunwall_unscl(:), &
    !$acc wtuq_sunwall_unscl(:), &
    !$acc wtus_shadewall_unscl(:), &
    !$acc wtuq_shadewall_unscl(:), &
    !$acc t_sunwall_innerl(:), &
    !$acc t_shadewall_innerl(:), &
    !$acc t_roof_innerl(:), &
    !$acc eflx_sh_grnd_scale(:), &
    !$acc qflx_evap_soi_scale(:), &
    !$acc eflx_wasteheat_roof(:), &
    !$acc eflx_wasteheat_sunwall(:), &
    !$acc eflx_wasteheat_shadewall(:), &
    !$acc eflx_heat_from_ac_roof(:), &
    !$acc eflx_heat_from_ac_sunwall(:), &
    !$acc eflx_heat_from_ac_shadewall(:), &
    !$acc eflx(:), &
    !$acc qflx(:), &
    !$acc eflx_scale(:), &
    !$acc qflx_scale(:), &
    !$acc eflx_err(:), &
    !$acc qflx_err(:), &
    !$acc local_secp1(:), &
    !$acc wind_speed0(:), &
    !$acc wind_speed_adj(:), &
    !$acc tau(:), &
    !$acc tau_diff(:), &
    !$acc prev_tau(:), &
    !$acc prev_tau_diff(:), &
    !$acc lnd_to_urban_filter(:), &
    !$acc converged_landunits(:), &
    !$acc num_copyl, &
    !$acc zeta, &
    !$acc fwet_roof, &
    !$acc fwet_road_imperv)

    end associate

  end subroutine UrbanFluxes

end module UrbanFluxesMod
