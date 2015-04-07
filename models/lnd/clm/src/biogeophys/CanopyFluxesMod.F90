module CanopyFluxesMod

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Performs calculation of leaf temperature and surface fluxes.
  ! SoilFluxes then determines soil/snow and ground temperatures and updates the surface 
  ! fluxes for the new ground temperature.
  !
  ! !USES:
  use shr_sys_mod           , only : shr_sys_flush
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use abortutils            , only : endrun
  use clm_varctl            , only : iulog, use_cn, use_lch4, use_c13, use_c14, use_cndv, use_ed
  use clm_varpar            , only : nlevgrnd, nlevsno
  use clm_varcon            , only : namep 
  use pftconMod             , only : nbrdlf_dcd_tmp_shrub, nsoybean , nsoybeanirrig, pftcon
  use decompMod             , only : bounds_type
  use PhotosynthesisMod     , only : Photosynthesis, PhotosynthesisTotal, Fractionation
  use EDPhotosynthesisMod   , only : Photosynthesis_ED
  use EDAccumulateFluxesMod , only : AccumulateFluxes_ED
  use EDBtranMod            , only : Btran_ED
  use SoilMoistStressMod    , only : calc_effective_soilporosity, calc_volumetric_h2oliq
  use SoilMoistStressMod    , only : calc_root_moist_stress, set_perchroot_opt
  use SimpleMathMod         , only : array_div_vector
  use SurfaceResistanceMod  , only : do_soilevap_beta
  use atm2lndType           , only : atm2lnd_type
  use CanopyStateType       , only : canopystate_type
  use CNVegStateType        , only : cnveg_state_type
  use EnergyFluxType        , only : energyflux_type
  use FrictionvelocityMod   , only : frictionvel_type
  use SoilStateType         , only : soilstate_type
  use SolarAbsorbedType     , only : solarabs_type
  use SurfaceAlbedoType     , only : surfalb_type
  use TemperatureType       , only : temperature_type
  use WaterfluxType         , only : waterflux_type
  use WaterstateType        , only : waterstate_type
  use HumanIndexMod         , only : humanindex_type
  use ch4Mod                , only : ch4_type
  use PhotosynthesisMod     , only : photosyns_type
  use GridcellType          , only : grc                
  use ColumnType            , only : col                
  use PatchType             , only : patch                
  use EDTypesMod            , only : ed_site_type
  use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CanopyFluxes
  !
  ! !PUBLIC DATA MEMBERS:
  ! true => btran is based only on unfrozen soil levels
  logical,  public :: perchroot     = .false.  

  ! true  => btran is based on active layer (defined over two years); 
  ! false => btran is based on currently unfrozen levels
  logical,  public :: perchroot_alt = .false.  
  !------------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  subroutine CanopyFluxes(bounds,  num_exposedvegp, filter_exposedvegp, &
       ed_allsites_inst,  atm2lnd_inst, canopystate_inst, cnveg_state_inst,            &
       energyflux_inst, frictionvel_inst, soilstate_inst, solarabs_inst, surfalb_inst, &
       temperature_inst, waterflux_inst, waterstate_inst, ch4_inst, photosyns_inst,    &
       humanindex_inst, soil_water_retention_curve) 
    !
    ! !DESCRIPTION:
    ! 1. Calculates the leaf temperature:
    ! 2. Calculates the leaf fluxes, transpiration, photosynthesis and
    !    updates the dew accumulation due to evaporation.
    !
    ! Method:
    ! Use the Newton-Raphson iteration to solve for the foliage
    ! temperature that balances the surface energy budget:
    !
    ! f(t_veg) = Net radiation - Sensible - Latent = 0
    ! f(t_veg) + d(f)/d(t_veg) * dt_veg = 0     (*)
    !
    ! Note:
    ! (1) In solving for t_veg, t_grnd is given from the previous timestep.
    ! (2) The partial derivatives of aerodynamical resistances, which cannot
    !     be determined analytically, are ignored for d(H)/dT and d(LE)/dT
    ! (3) The weighted stomatal resistance of sunlit and shaded foliage is used
    ! (4) Canopy air temperature and humidity are derived from => Hc + Hg = Ha
    !                                                          => Ec + Eg = Ea
    ! (5) Energy loss is due to: numerical truncation of energy budget equation
    !     (*); and "ecidif" (see the code) which is dropped into the sensible
    !     heat
    ! (6) The convergence criteria: the difference, del = t_veg(n+1)-t_veg(n)
    !     and del2 = t_veg(n)-t_veg(n-1) less than 0.01 K, and the difference
    !     of water flux from the leaf between the iteration step (n+1) and (n)
    !     less than 0.1 W/m2; or the iterative steps over 40.
    !
    ! !USES:
    use shr_const_mod      , only : SHR_CONST_RGAS
    use clm_time_manager   , only : get_step_size, get_prev_date
    use clm_varcon         , only : sb, cpair, hvap, vkc, grav, denice
    use clm_varcon         , only : denh2o, tfrz, csoilc, tlsai_crit, alpha_aero
    use clm_varcon         , only : c14ratio
    use perf_mod           , only : t_startf, t_stopf
    use QSatMod            , only : QSat
    use FrictionVelocityMod, only : FrictionVelocity, MoninObukIni
    use HumanIndexMod      , only : calc_human_stress_indices, Wet_Bulb, Wet_BulbS, HeatIndex, AppTemp, &
                                    swbgt, hmdex, dis_coi, dis_coiS, THIndex, &
                                    SwampCoolEff, KtoC, VaporPres
    use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
    !
    ! !ARGUMENTS:
    type(bounds_type)                      , intent(in)            :: bounds 
    integer                                , intent(in)            :: num_exposedvegp        ! number of points in filter_exposedvegp
    integer                                , intent(in)            :: filter_exposedvegp(:)  ! patch filter for non-snow-covered veg
    type(ed_site_type)                     , intent(inout), target :: ed_allsites_inst( bounds%begg: )
    type(atm2lnd_type)                     , intent(in)            :: atm2lnd_inst
    type(canopystate_type)                 , intent(inout)         :: canopystate_inst
    type(cnveg_state_type)                 , intent(in)            :: cnveg_state_inst
    type(energyflux_type)                  , intent(inout)         :: energyflux_inst
    type(frictionvel_type)                 , intent(inout)         :: frictionvel_inst
    type(solarabs_type)                    , intent(in)            :: solarabs_inst
    type(surfalb_type)                     , intent(in)            :: surfalb_inst
    type(soilstate_type)                   , intent(inout)         :: soilstate_inst
    type(temperature_type)                 , intent(inout)         :: temperature_inst
    type(waterstate_type)                  , intent(inout)         :: waterstate_inst
    type(waterflux_type)                   , intent(inout)         :: waterflux_inst
    type(ch4_type)                         , intent(inout)         :: ch4_inst
    type(photosyns_type)                   , intent(inout)         :: photosyns_inst
    type(humanindex_type)                  , intent(inout)         :: humanindex_inst
    class(soil_water_retention_curve_type) , intent(in)            :: soil_water_retention_curve
    !
    ! !LOCAL VARIABLES:
    real(r8), parameter :: btran0 = 0.0_r8  ! initial value
    real(r8), parameter :: zii = 1000.0_r8  ! convective boundary layer height [m]
    real(r8), parameter :: beta = 1.0_r8    ! coefficient of conective velocity [-]
    real(r8), parameter :: delmax = 1.0_r8  ! maxchange in  leaf temperature [K]
    real(r8), parameter :: dlemin = 0.1_r8  ! max limit for energy flux convergence [w/m2]
    real(r8), parameter :: dtmin = 0.01_r8  ! max limit for temperature convergence [K]
    integer , parameter :: itmax = 40       ! maximum number of iteration [-]
    integer , parameter :: itmin = 2        ! minimum number of iteration [-]

    !added by K.Sakaguchi for litter resistance
    real(r8), parameter :: lai_dl = 0.5_r8           ! placeholder for (dry) plant litter area index (m2/m2)
    real(r8), parameter :: z_dl = 0.05_r8            ! placeholder for (dry) litter layer thickness (m)

    !added by K.Sakaguchi for stability formulation
    real(r8), parameter :: ria  = 0.5_r8             ! free parameter for stable formulation (currently = 0.5, "gamma" in Sakaguchi&Zeng,2008)

    real(r8) :: dtime                                ! land model time step (sec)
    real(r8) :: zldis(bounds%begp:bounds%endp)       ! reference height "minus" zero displacement height [m]
    real(r8) :: zeta                                 ! dimensionless height used in Monin-Obukhov theory
    real(r8) :: wc                                   ! convective velocity [m/s]
    real(r8) :: dth(bounds%begp:bounds%endp)         ! diff of virtual temp. between ref. height and surface
    real(r8) :: dthv(bounds%begp:bounds%endp)        ! diff of vir. poten. temp. between ref. height and surface
    real(r8) :: dqh(bounds%begp:bounds%endp)         ! diff of humidity between ref. height and surface
    real(r8) :: obu(bounds%begp:bounds%endp)         ! Monin-Obukhov length (m)
    real(r8) :: um(bounds%begp:bounds%endp)          ! wind speed including the stablity effect [m/s]
    real(r8) :: ur(bounds%begp:bounds%endp)          ! wind speed at reference height [m/s]
    real(r8) :: uaf(bounds%begp:bounds%endp)         ! velocity of air within foliage [m/s]
    real(r8) :: temp1(bounds%begp:bounds%endp)       ! relation for potential temperature profile
    real(r8) :: temp12m(bounds%begp:bounds%endp)     ! relation for potential temperature profile applied at 2-m
    real(r8) :: temp2(bounds%begp:bounds%endp)       ! relation for specific humidity profile
    real(r8) :: temp22m(bounds%begp:bounds%endp)     ! relation for specific humidity profile applied at 2-m
    real(r8) :: ustar(bounds%begp:bounds%endp)       ! friction velocity [m/s]
    real(r8) :: tstar                                ! temperature scaling parameter
    real(r8) :: qstar                                ! moisture scaling parameter
    real(r8) :: thvstar                              ! virtual potential temperature scaling parameter
    real(r8) :: taf(bounds%begp:bounds%endp)         ! air temperature within canopy space [K]
    real(r8) :: qaf(bounds%begp:bounds%endp)         ! humidity of canopy air [kg/kg]
    real(r8) :: rpp                                  ! fraction of potential evaporation from leaf [-]
    real(r8) :: rppdry                               ! fraction of potential evaporation through transp [-]
    real(r8) :: cf                                   ! heat transfer coefficient from leaves [-]
    real(r8) :: rb(bounds%begp:bounds%endp)          ! leaf boundary layer resistance [s/m]
    real(r8) :: rah(bounds%begp:bounds%endp,2)       ! thermal resistance [s/m]
    real(r8) :: raw(bounds%begp:bounds%endp,2)       ! moisture resistance [s/m]
    real(r8) :: wta                                  ! heat conductance for air [m/s]
    real(r8) :: wtg(bounds%begp:bounds%endp)         ! heat conductance for ground [m/s]
    real(r8) :: wtl                                  ! heat conductance for leaf [m/s]
    real(r8) :: wta0(bounds%begp:bounds%endp)        ! normalized heat conductance for air [-]
    real(r8) :: wtl0(bounds%begp:bounds%endp)        ! normalized heat conductance for leaf [-]
    real(r8) :: wtg0                                 ! normalized heat conductance for ground [-]
    real(r8) :: wtal(bounds%begp:bounds%endp)        ! normalized heat conductance for air and leaf [-]
    real(r8) :: wtga                                 ! normalized heat cond. for air and ground  [-]
    real(r8) :: wtaq                                 ! latent heat conductance for air [m/s]
    real(r8) :: wtlq                                 ! latent heat conductance for leaf [m/s]
    real(r8) :: wtgq(bounds%begp:bounds%endp)        ! latent heat conductance for ground [m/s]
    real(r8) :: wtaq0(bounds%begp:bounds%endp)       ! normalized latent heat conductance for air [-]
    real(r8) :: wtlq0(bounds%begp:bounds%endp)       ! normalized latent heat conductance for leaf [-]
    real(r8) :: wtgq0                                ! normalized heat conductance for ground [-]
    real(r8) :: wtalq(bounds%begp:bounds%endp)       ! normalized latent heat cond. for air and leaf [-]
    real(r8) :: wtgaq                                ! normalized latent heat cond. for air and ground [-]
    real(r8) :: el(bounds%begp:bounds%endp)          ! vapor pressure on leaf surface [pa]
    real(r8) :: deldT                                ! derivative of "el" on "t_veg" [pa/K]
    real(r8) :: qsatl(bounds%begp:bounds%endp)       ! leaf specific humidity [kg/kg]
    real(r8) :: qsatldT(bounds%begp:bounds%endp)     ! derivative of "qsatl" on "t_veg"
    real(r8) :: e_ref2m                              ! 2 m height surface saturated vapor pressure [Pa]
    real(r8) :: de2mdT                               ! derivative of 2 m height surface saturated vapor pressure on t_ref2m
    real(r8) :: qsat_ref2m                           ! 2 m height surface saturated specific humidity [kg/kg]
    real(r8) :: dqsat2mdT                            ! derivative of 2 m height surface saturated specific humidity on t_ref2m
    real(r8) :: air(bounds%begp:bounds%endp)         ! atmos. radiation temporay set
    real(r8) :: bir(bounds%begp:bounds%endp)         ! atmos. radiation temporay set
    real(r8) :: cir(bounds%begp:bounds%endp)         ! atmos. radiation temporay set
    real(r8) :: dc1,dc2                              ! derivative of energy flux [W/m2/K]
    real(r8) :: delt                                 ! temporary
    real(r8) :: delq(bounds%begp:bounds%endp)        ! temporary
    real(r8) :: del(bounds%begp:bounds%endp)         ! absolute change in leaf temp in current iteration [K]
    real(r8) :: del2(bounds%begp:bounds%endp)        ! change in leaf temperature in previous iteration [K]
    real(r8) :: dele(bounds%begp:bounds%endp)        ! change in latent heat flux from leaf [K]
    real(r8) :: dels                                 ! change in leaf temperature in current iteration [K]
    real(r8) :: det(bounds%begp:bounds%endp)         ! maximum leaf temp. change in two consecutive iter [K]
    real(r8) :: efeb(bounds%begp:bounds%endp)        ! latent heat flux from leaf (previous iter) [mm/s]
    real(r8) :: efeold                               ! latent heat flux from leaf (previous iter) [mm/s]
    real(r8) :: efpot                                ! potential latent energy flux [kg/m2/s]
    real(r8) :: efe(bounds%begp:bounds%endp)         ! water flux from leaf [mm/s]
    real(r8) :: efsh                                 ! sensible heat from leaf [mm/s]
    real(r8) :: obuold(bounds%begp:bounds%endp)      ! monin-obukhov length from previous iteration
    real(r8) :: tlbef(bounds%begp:bounds%endp)       ! leaf temperature from previous iteration [K]
    real(r8) :: ecidif                               ! excess energies [W/m2]
    real(r8) :: err(bounds%begp:bounds%endp)         ! balance error
    real(r8) :: erre                                 ! balance error
    real(r8) :: co2(bounds%begp:bounds%endp)         ! atmospheric co2 partial pressure (pa)
    real(r8) :: c13o2(bounds%begp:bounds%endp)       ! atmospheric c13o2 partial pressure (pa)
    real(r8) :: o2(bounds%begp:bounds%endp)          ! atmospheric o2 partial pressure (pa)
    real(r8) :: svpts(bounds%begp:bounds%endp)       ! saturation vapor pressure at t_veg (pa)
    real(r8) :: eah(bounds%begp:bounds%endp)         ! canopy air vapor pressure (pa)
    real(r8) :: s_node                               ! vol_liq/eff_porosity
    real(r8) :: smp_node                             ! matrix potential
    real(r8) :: smp_node_lf                          ! F. Li and S. Levis
    real(r8) :: vol_liq                              ! partial volume of liquid water in layer
    integer  :: itlef                                ! counter for leaf temperature iteration [-]
    integer  :: nmozsgn(bounds%begp:bounds%endp)     ! number of times stability changes sign
    real(r8) :: w                                    ! exp(-LSAI)
    real(r8) :: csoilcn                              ! interpolated csoilc for less than dense canopies
    real(r8) :: fm(bounds%begp:bounds%endp)          ! needed for BGC only to diagnose 10m wind speed
    real(r8) :: wtshi                                ! sensible heat resistance for air, grnd and leaf [-]
    real(r8) :: wtsqi                                ! latent heat resistance for air, grnd and leaf [-]
    integer  :: j                                    ! soil/snow level index
    integer  :: p                                    ! patch index
    integer  :: c                                    ! column index
    integer  :: l                                    ! landunit index
    integer  :: g                                    ! gridcell index
    integer  :: fn                                   ! number of values in vegetated patch filter
    integer  :: filterp(bounds%endp-bounds%begp+1)   ! vegetated patch filter
    integer  :: fnorig                               ! number of values in patch filter copy
    integer  :: fporig(bounds%endp-bounds%begp+1)    ! temporary filter
    integer  :: fnold                                ! temporary copy of patch count
    integer  :: f                                    ! filter index
    logical  :: found                                ! error flag for canopy above forcing hgt
    integer  :: index                                ! patch index for error
    real(r8) :: egvf                                 ! effective green vegetation fraction
    real(r8) :: lt                                   ! elai+esai
    real(r8) :: ri                                   ! stability parameter for under canopy air (unitless)
    real(r8) :: csoilb                               ! turbulent transfer coefficient over bare soil (unitless)
    real(r8) :: ricsoilc                             ! modified transfer coefficient under dense canopy (unitless)
    real(r8) :: snow_depth_c                         ! critical snow depth to cover plant litter (m)
    real(r8) :: rdl                                  ! dry litter layer resistance for water vapor  (s/m)
    real(r8) :: elai_dl                              ! exposed (dry) plant litter area index
    real(r8) :: fsno_dl                              ! effective snow cover over plant litter
    real(r8) :: dayl_factor(bounds%begp:bounds%endp) ! scalar (0-1) for daylength effect on Vcmax
    ! If no unfrozen layers, put all in the top layer.
    real(r8) :: rootsum(bounds%begp:bounds%endp)
    real(r8) :: delt_snow
    real(r8) :: delt_soil
    real(r8) :: delt_h2osfc
    real(r8) :: lw_grnd
    real(r8) :: delq_snow
    real(r8) :: delq_soil
    real(r8) :: delq_h2osfc
    real(r8) :: dt_veg(bounds%begp:bounds%endp)          ! change in t_veg, last iteration (Kelvin)                              
    integer  :: jtop(bounds%begc:bounds%endc)            ! lbning
    integer  :: filterc_tmp(bounds%endp-bounds%begp+1)   ! temporary variable
    integer  :: ft                                       ! plant functional type index
    real(r8) :: temprootr                 
    real(r8) :: dt_veg_temp(bounds%begp:bounds%endp)
    integer  :: iv
    !------------------------------------------------------------------------------

    associate(                                                                    & 
         snl                    => col%snl                                      , & ! Input:  [integer  (:)   ]  number of snow layers                                                  
         dayl                   => grc%dayl                                     , & ! Input:  [real(r8) (:)   ]  daylength (s)
         max_dayl               => grc%max_dayl                                 , & ! Input:  [real(r8) (:)   ]  maximum daylength for this grid cell (s)

         dleaf                  => pftcon%dleaf                                 , & ! Input:  characteristic leaf dimension (m)                                     

         forc_lwrad             => atm2lnd_inst%forc_lwrad_downscaled_col       , & ! Input:  [real(r8) (:)   ]  downward infrared (longwave) radiation (W/m**2)                       
         forc_q                 => atm2lnd_inst%forc_q_downscaled_col           , & ! Input:  [real(r8) (:)   ]  atmospheric specific humidity (kg/kg)                                 
         forc_pbot              => atm2lnd_inst%forc_pbot_downscaled_col        , & ! Input:  [real(r8) (:)   ]  atmospheric pressure (Pa)                                             
         forc_th                => atm2lnd_inst%forc_th_downscaled_col          , & ! Input:  [real(r8) (:)   ]  atmospheric potential temperature (Kelvin)                            
         forc_rho               => atm2lnd_inst%forc_rho_downscaled_col         , & ! Input:  [real(r8) (:)   ]  density (kg/m**3)                                                     
         forc_t                 => atm2lnd_inst%forc_t_downscaled_col           , & ! Input:  [real(r8) (:)   ]  atmospheric temperature (Kelvin)                                      
         forc_u                 => atm2lnd_inst%forc_u_grc                      , & ! Input:  [real(r8) (:)   ]  atmospheric wind speed in east direction (m/s)                        
         forc_v                 => atm2lnd_inst%forc_v_grc                      , & ! Input:  [real(r8) (:)   ]  atmospheric wind speed in north direction (m/s)                       
         forc_pco2              => atm2lnd_inst%forc_pco2_grc                   , & ! Input:  [real(r8) (:)   ]  partial pressure co2 (Pa)                                             
         forc_pc13o2            => atm2lnd_inst%forc_pc13o2_grc                 , & ! Input:  [real(r8) (:)   ]  partial pressure c13o2 (Pa)                                           
         forc_po2               => atm2lnd_inst%forc_po2_grc                    , & ! Input:  [real(r8) (:)   ]  partial pressure o2 (Pa)                                              

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

         sabv                   => solarabs_inst%sabv_patch                     , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by vegetation (W/m**2)                       

         frac_veg_nosno         => canopystate_inst%frac_veg_nosno_patch        , & ! Input:  [integer  (:)   ]  fraction of vegetation not covered by snow (0 OR 1) [-]
         elai                   => canopystate_inst%elai_patch                  , & ! Input:  [real(r8) (:)   ]  one-sided leaf area index with burying by snow                        
         esai                   => canopystate_inst%esai_patch                  , & ! Input:  [real(r8) (:)   ]  one-sided stem area index with burying by snow                        
         laisun                 => canopystate_inst%laisun_patch                , & ! Input:  [real(r8) (:)   ]  sunlit leaf area                                                      
         laisha                 => canopystate_inst%laisha_patch                , & ! Input:  [real(r8) (:)   ]  shaded leaf area                                                      
         displa                 => canopystate_inst%displa_patch                , & ! Input:  [real(r8) (:)   ]  displacement height (m)                                               
         htop                   => canopystate_inst%htop_patch                  , & ! Input:  [real(r8) (:)   ]  canopy top(m)                                                         
         altmax_lastyear_indx   => canopystate_inst%altmax_lastyear_indx_col    , & ! Input:  [integer  (:)   ]  prior year maximum annual depth of thaw                                
         altmax_indx            => canopystate_inst%altmax_indx_col             , & ! Input:  [integer  (:)   ]  maximum annual depth of thaw                                           
         rscanopy               => canopystate_inst%rscanopy_patch              , & ! Output: [real(r8) (:,:)]   canopy resistance s/m (ED)
         
         watsat                 => soilstate_inst%watsat_col                    , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)   (constant)                     
         watdry                 => soilstate_inst%watdry_col                    , & ! Input:  [real(r8) (:,:) ]  btran parameter for btran=0                      (constant)                                        
         watopt                 => soilstate_inst%watopt_col                    , & ! Input:  [real(r8) (:,:) ]  btran parameter for btran=1                      (constant)                                      
         eff_porosity           => soilstate_inst%eff_porosity_col              , & ! Output: [real(r8) (:,:) ]  effective soil porosity
         soilbeta               => soilstate_inst%soilbeta_col                  , & ! Input:  [real(r8) (:)   ]  soil wetness relative to field capacity                               
         rootr                  => soilstate_inst%rootr_patch                   , & ! Output: [real(r8) (:,:) ]  effective fraction of roots in each soil layer                      

         u10_clm                => frictionvel_inst%u10_clm_patch               , & ! Input:  [real(r8) (:)   ]  10 m height winds (m/s)
         forc_hgt_u_patch       => frictionvel_inst%forc_hgt_u_patch            , & ! Input:  [real(r8) (:)   ]  observational height of wind at patch level [m]                          
         z0mg                   => frictionvel_inst%z0mg_col                    , & ! Input:  [real(r8) (:)   ]  roughness length of ground, momentum [m]                              
         ram1                   => frictionvel_inst%ram1_patch                  , & ! Output: [real(r8) (:)   ]  aerodynamical resistance (s/m)                                        
         z0mv                   => frictionvel_inst%z0mv_patch                  , & ! Output: [real(r8) (:)   ]  roughness length over vegetation, momentum [m]                        
         z0hv                   => frictionvel_inst%z0hv_patch                  , & ! Output: [real(r8) (:)   ]  roughness length over vegetation, sensible heat [m]                   
         z0qv                   => frictionvel_inst%z0qv_patch                  , & ! Output: [real(r8) (:)   ]  roughness length over vegetation, latent heat [m]                     
         rb1                    => frictionvel_inst%rb1_patch                   , & ! Output: [real(r8) (:)   ]  boundary layer resistance (s/m)                                       

         t_h2osfc               => temperature_inst%t_h2osfc_col                , & ! Input:  [real(r8) (:)   ]  surface water temperature                                             
         t_soisno               => temperature_inst%t_soisno_col                , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)                                           
         t_grnd                 => temperature_inst%t_grnd_col                  , & ! Input:  [real(r8) (:)   ]  ground surface temperature [K]                                        
         thv                    => temperature_inst%thv_col                     , & ! Input:  [real(r8) (:)   ]  virtual potential temperature (kelvin)                                
         thm                    => temperature_inst%thm_patch                   , & ! Input:  [real(r8) (:)   ]  intermediate variable (forc_t+0.0098*forc_hgt_t_patch)                  
         emv                    => temperature_inst%emv_patch                   , & ! Input:  [real(r8) (:)   ]  vegetation emissivity                                                     
         emg                    => temperature_inst%emg_col                     , & ! Input:  [real(r8) (:)   ]  vegetation emissivity                                                 
         t_veg                  => temperature_inst%t_veg_patch                 , & ! Output: [real(r8) (:)   ]  vegetation temperature (Kelvin)                                       
         t_ref2m                => temperature_inst%t_ref2m_patch               , & ! Output: [real(r8) (:)   ]  2 m height surface air temperature (Kelvin)                           
         t_ref2m_r              => temperature_inst%t_ref2m_r_patch             , & ! Output: [real(r8) (:)   ]  Rural 2 m height surface air temperature (Kelvin)                     

         frac_h2osfc            => waterstate_inst%frac_h2osfc_col              , & ! Input:  [real(r8) (:)   ]  fraction of surface water                                             
         fwet                   => waterstate_inst%fwet_patch                   , & ! Input:  [real(r8) (:)   ]  fraction of canopy that is wet (0 to 1)                               
         fdry                   => waterstate_inst%fdry_patch                   , & ! Input:  [real(r8) (:)   ]  fraction of foliage that is green and dry [-]                         
         frac_sno               => waterstate_inst%frac_sno_eff_col             , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)                           
         snow_depth             => waterstate_inst%snow_depth_col               , & ! Input:  [real(r8) (:)   ]  snow height (m)                                                       
         qg_snow                => waterstate_inst%qg_snow_col                  , & ! Input:  [real(r8) (:)   ]  specific humidity at snow surface [kg/kg]                             
         qg_soil                => waterstate_inst%qg_soil_col                  , & ! Input:  [real(r8) (:)   ]  specific humidity at soil surface [kg/kg]                             
         qg_h2osfc              => waterstate_inst%qg_h2osfc_col                , & ! Input:  [real(r8) (:)   ]  specific humidity at h2osfc surface [kg/kg]                           
         qg                     => waterstate_inst%qg_col                       , & ! Input:  [real(r8) (:)   ]  specific humidity at ground surface [kg/kg]                           
         dqgdT                  => waterstate_inst%dqgdT_col                    , & ! Input:  [real(r8) (:)   ]  temperature derivative of "qg"                                        
         h2osoi_ice             => waterstate_inst%h2osoi_ice_col               , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)                                                    
         h2osoi_vol             => waterstate_inst%h2osoi_vol_col               , & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3] by F. Li and S. Levis
         h2osoi_liq             => waterstate_inst%h2osoi_liq_col               , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)                                                
         h2osoi_liqvol          => waterstate_inst%h2osoi_liqvol_col            , & ! Output: [real(r8) (:,:) ]  volumetric liquid water (v/v) 
         h2ocan                 => waterstate_inst%h2ocan_patch                 , & ! Output: [real(r8) (:)   ]  canopy water (mm H2O)                                                 
         q_ref2m                => waterstate_inst%q_ref2m_patch                , & ! Output: [real(r8) (:)   ]  2 m height surface specific humidity (kg/kg)                          
         rh_ref2m_r             => waterstate_inst%rh_ref2m_r_patch             , & ! Output: [real(r8) (:)   ]  Rural 2 m height surface relative humidity (%)                        
         rh_ref2m               => waterstate_inst%rh_ref2m_patch               , & ! Output: [real(r8) (:)   ]  2 m height surface relative humidity (%)                              
         rhaf                   => waterstate_inst%rh_af_patch                  , & ! Output: [real(r8) (:)   ]  fractional humidity of canopy air [dimensionless]                     

         qflx_tran_veg          => waterflux_inst%qflx_tran_veg_patch           , & ! Output: [real(r8) (:)   ]  vegetation transpiration (mm H2O/s) (+ = to atm)                      
         qflx_evap_veg          => waterflux_inst%qflx_evap_veg_patch           , & ! Output: [real(r8) (:)   ]  vegetation evaporation (mm H2O/s) (+ = to atm)                        
         qflx_evap_soi          => waterflux_inst%qflx_evap_soi_patch           , & ! Output: [real(r8) (:)   ]  soil evaporation (mm H2O/s) (+ = to atm)                              
         qflx_ev_snow           => waterflux_inst%qflx_ev_snow_patch            , & ! Output: [real(r8) (:)   ]  evaporation flux from snow (W/m**2) [+ to atm]                        
         qflx_ev_soil           => waterflux_inst%qflx_ev_soil_patch            , & ! Output: [real(r8) (:)   ]  evaporation flux from soil (W/m**2) [+ to atm]                        
         qflx_ev_h2osfc         => waterflux_inst%qflx_ev_h2osfc_patch          , & ! Output: [real(r8) (:)   ]  evaporation flux from h2osfc (W/m**2) [+ to atm]                      

         rssun                  => photosyns_inst%rssun_patch                   , & ! Output: [real(r8) (:)   ]  leaf sunlit stomatal resistance (s/m) (output from Photosynthesis)
         rssha                  => photosyns_inst%rssha_patch                   , & ! Output: [real(r8) (:)   ]  leaf shaded stomatal resistance (s/m) (output from Photosynthesis)

         grnd_ch4_cond          => ch4_inst%grnd_ch4_cond_patch                 , & ! Output: [real(r8) (:)   ]  tracer conductance for boundary layer [m/s] 

         htvp                   => energyflux_inst%htvp_col                     , & ! Input:  [real(r8) (:)   ]  latent heat of evaporation (/sublimation) [J/kg] (constant)                      
         btran2                 => energyflux_inst%btran2_patch                 , & ! Output: [real(r8) (:)   ]  F. Li and S. Levis                                                     
         btran                  => energyflux_inst%btran_patch                  , & ! Output: [real(r8) (:)   ]  transpiration wetness factor (0 to 1)                                 
         rresis                 => energyflux_inst%rresis_patch                 , & ! Output: [real(r8) (:,:) ]  root resistance by layer (0-1)  (nlevgrnd)                          
         taux                   => energyflux_inst%taux_patch                   , & ! Output: [real(r8) (:)   ]  wind (shear) stress: e-w (kg/m/s**2)                                  
         tauy                   => energyflux_inst%tauy_patch                   , & ! Output: [real(r8) (:)   ]  wind (shear) stress: n-s (kg/m/s**2)                                  
         canopy_cond            => energyflux_inst%canopy_cond_patch            , & ! Output: [real(r8) (:)   ]  tracer conductance for canopy [m/s] 
         cgrnds                 => energyflux_inst%cgrnds_patch                 , & ! Output: [real(r8) (:)   ]  deriv. of soil sensible heat flux wrt soil temp [w/m2/k]              
         cgrndl                 => energyflux_inst%cgrndl_patch                 , & ! Output: [real(r8) (:)   ]  deriv. of soil latent heat flux wrt soil temp [w/m**2/k]              
         dlrad                  => energyflux_inst%dlrad_patch                  , & ! Output: [real(r8) (:)   ]  downward longwave radiation below the canopy [W/m2]                   
         ulrad                  => energyflux_inst%ulrad_patch                  , & ! Output: [real(r8) (:)   ]  upward longwave radiation above the canopy [W/m2]                     
         cgrnd                  => energyflux_inst%cgrnd_patch                  , & ! Output: [real(r8) (:)   ]  deriv. of soil energy flux wrt to soil temp [w/m2/k]                  
         eflx_sh_snow           => energyflux_inst%eflx_sh_snow_patch           , & ! Output: [real(r8) (:)   ]  sensible heat flux from snow (W/m**2) [+ to atm]                      
         eflx_sh_h2osfc         => energyflux_inst%eflx_sh_h2osfc_patch         , & ! Output: [real(r8) (:)   ]  sensible heat flux from soil (W/m**2) [+ to atm]                      
         eflx_sh_soil           => energyflux_inst%eflx_sh_soil_patch           , & ! Output: [real(r8) (:)   ]  sensible heat flux from soil (W/m**2) [+ to atm]                      
         eflx_sh_veg            => energyflux_inst%eflx_sh_veg_patch            , & ! Output: [real(r8) (:)   ]  sensible heat flux from leaves (W/m**2) [+ to atm]                    
         eflx_sh_grnd           => energyflux_inst%eflx_sh_grnd_patch           , & ! Output: [real(r8) (:)   ]  sensible heat flux from ground (W/m**2) [+ to atm]                    

         begp                   => bounds%begp                                  , &
         endp                   => bounds%endp                                  , &
         begg                   => bounds%begg                                  , &
         endg                   => bounds%endg                                    &
         )

      ! Determine step size

      dtime = get_step_size()

      ! Make a local copy of the exposedvegp filter. With the current implementation,
      ! this is needed because the filter is modified in the iteration loop.
      !
      ! TODO(wjs, 2014-09-24) Determine if this is really needed. I suspect that we could
      ! do away with either this temporary fn/filterp, or the temporary fnorig/fporig,
      ! with one of these simply using the passed-in filter (num_exposedvegp /
      ! filter_exposedvegp)

      fn = num_exposedvegp
      filterp(1:fn) = filter_exposedvegp(1:fn)

      ! -----------------------------------------------------------------
      ! Time step initialization of photosynthesis variables
      ! -----------------------------------------------------------------

      call photosyns_inst%TimeStepInit(bounds)

      ! Initialize

      do f = 1, fn
         p = filterp(f)
         del(p)    = 0._r8  ! change in leaf temperature from previous iteration
         efeb(p)   = 0._r8  ! latent head flux from leaf for previous iteration
         wtlq0(p)  = 0._r8
         wtalq(p)  = 0._r8
         wtgq(p)   = 0._r8
         wtaq0(p)  = 0._r8
         obuold(p) = 0._r8
         btran(p)  = btran0
         btran2(p)  = btran0
      end do

      ! calculate daylength control for Vcmax
      do f = 1, fn
         p=filterp(f)
         g=patch%gridcell(p)
         ! calculate dayl_factor as the ratio of (current:max dayl)^2
         ! set a minimum of 0.01 (1%) for the dayl_factor
         dayl_factor(p)=min(1._r8,max(0.01_r8,(dayl(g)*dayl(g))/(max_dayl(g)*max_dayl(g))))
      end do

      rb1(begp:endp) = 0._r8

      ! FIX(FIX(SPM,032414),032414) refactor this...
      if ( use_ed ) then  

         do f = 1, fn
            p = filterp(f)
            call btran_ed(bounds, p, ed_allsites_inst(begg:endg), &
                 soilstate_inst, waterstate_inst, temperature_inst, energyflux_inst)
         enddo

      else

         !assign the temporary filter
         do f = 1, fn
            p = filterp(f)
            filterc_tmp(f)=patch%column(p)
         enddo

         !compute effective soil porosity
         call calc_effective_soilporosity(bounds,                          &
              ubj = nlevgrnd,                                              &
              numf = fn,                                                   &
              filter = filterc_tmp(1:fn),                                  &
              watsat = watsat(bounds%begc:bounds%endc, 1:nlevgrnd),        &
              h2osoi_ice = h2osoi_ice(bounds%begc:bounds%endc,1:nlevgrnd), &
              denice = denice,                                             &
              eff_por=eff_porosity(bounds%begc:bounds%endc, 1:nlevgrnd) )

         !compute volumetric liquid water content
         jtop(bounds%begc:bounds%endc) = 1

         call calc_volumetric_h2oliq(bounds,                                    &
              jtop = jtop(bounds%begc:bounds%endc),                             &
              lbj = 1,                                                          &
              ubj = nlevgrnd,                                                   &
              numf = fn,                                                        &
              filter = filterc_tmp(1:fn),                                       &
              eff_porosity = eff_porosity(bounds%begc:bounds%endc, 1:nlevgrnd), &
              h2osoi_liq = h2osoi_liq(bounds%begc:bounds%endc, 1:nlevgrnd),     &
              denh2o = denh2o,                                                  &
              vol_liq = h2osoi_liqvol(bounds%begc:bounds%endc, 1:nlevgrnd) )

         !set up perchroot options
         call set_perchroot_opt(perchroot, perchroot_alt)

         !calculate root moisture stress
         call calc_root_moist_stress(bounds,     &
              nlevgrnd = nlevgrnd,               &
              fn = fn,                           &
              filterp = filterp,                 &
              canopystate_inst=canopystate_inst, &
              energyflux_inst=energyflux_inst,   &
              soilstate_inst=soilstate_inst,     &
              temperature_inst=temperature_inst, &
              waterstate_inst=waterstate_inst,   &
              soil_water_retention_curve=soil_water_retention_curve)

      end if !use_ed

      ! Modify aerodynamic parameters for sparse/dense canopy (X. Zeng)
      do f = 1, fn
         p = filterp(f)
         c = patch%column(p)

         lt = min(elai(p)+esai(p), tlsai_crit)
         egvf =(1._r8 - alpha_aero * exp(-lt)) / (1._r8 - alpha_aero * exp(-tlsai_crit))
         displa(p) = egvf * displa(p)
         z0mv(p)   = exp(egvf * log(z0mv(p)) + (1._r8 - egvf) * log(z0mg(c)))
         z0hv(p)   = z0mv(p)
         z0qv(p)   = z0mv(p)
      end do

      found = .false.
      do f = 1, fn
         p = filterp(f)
         c = patch%column(p)
         g = patch%gridcell(p)

         ! Net absorbed longwave radiation by canopy and ground
         ! =air+bir*t_veg**4+cir*t_grnd(c)**4

         air(p) =   emv(p) * (1._r8+(1._r8-emv(p))*(1._r8-emg(c))) * forc_lwrad(c)
         bir(p) = - (2._r8-emv(p)*(1._r8-emg(c))) * emv(p) * sb
         cir(p) =   emv(p)*emg(c)*sb

         ! Saturated vapor pressure, specific humidity, and their derivatives
         ! at the leaf surface

         call QSat (t_veg(p), forc_pbot(c), el(p), deldT, qsatl(p), qsatldT(p))

         ! Determine atmospheric co2 and o2

         co2(p) = forc_pco2(g)
         o2(p)  = forc_po2(g)

         if ( use_c13 ) then
            c13o2(p) = forc_pc13o2(g)
         end if

         ! Initialize flux profile

         nmozsgn(p) = 0

         taf(p) = (t_grnd(c) + thm(p))/2._r8
         qaf(p) = (forc_q(c)+qg(c))/2._r8

         ur(p) = max(1.0_r8,sqrt(forc_u(g)*forc_u(g)+forc_v(g)*forc_v(g)))
         dth(p) = thm(p)-taf(p)
         dqh(p) = forc_q(c)-qaf(p)
         delq(p) = qg(c) - qaf(p)
         dthv(p) = dth(p)*(1._r8+0.61_r8*forc_q(c))+0.61_r8*forc_th(c)*dqh(p)
         zldis(p) = forc_hgt_u_patch(p) - displa(p)

         ! Check to see if the forcing height is below the canopy height
         if (zldis(p) < 0._r8) then
            found = .true.
            index = p
         end if

      end do

      if (found) then
         if ( .not. use_ed ) then
            write(iulog,*)'Error: Forcing height is below canopy height for patch index '
            call endrun(decomp_index=index, clmlevel=namep, msg=errmsg(__FILE__, __LINE__))
         end if
      end if

      do f = 1, fn
         p = filterp(f)
         c = patch%column(p)

         ! Initialize Monin-Obukhov length and wind speed

         call MoninObukIni(ur(p), thv(c), dthv(p), zldis(p), z0mv(p), um(p), obu(p))

      end do

      ! Set counter for leaf temperature iteration (itlef)

      itlef = 0    
      fnorig = fn
      fporig(1:fn) = filterp(1:fn)

      ! Begin stability iteration

      call t_startf('can_iter')
      ITERATION : do while (itlef <= itmax .and. fn > 0)

         ! Determine friction velocity, and potential temperature and humidity
         ! profiles of the surface boundary layer

         call FrictionVelocity (begp, endp, fn, filterp, &
              displa(begp:endp), z0mv(begp:endp), z0hv(begp:endp), z0qv(begp:endp), &
              obu(begp:endp), itlef+1, ur(begp:endp), um(begp:endp), ustar(begp:endp), &
              temp1(begp:endp), temp2(begp:endp), temp12m(begp:endp), temp22m(begp:endp), fm(begp:endp), &
              frictionvel_inst)

         do f = 1, fn
            p = filterp(f)
            c = patch%column(p)
            g = patch%gridcell(p)

            tlbef(p) = t_veg(p)
            del2(p) = del(p)

            ! Determine aerodynamic resistances

            ram1(p)  = 1._r8/(ustar(p)*ustar(p)/um(p))
            rah(p,1) = 1._r8/(temp1(p)*ustar(p))
            raw(p,1) = 1._r8/(temp2(p)*ustar(p))

            ! Bulk boundary layer resistance of leaves

            uaf(p) = um(p)*sqrt( 1._r8/(ram1(p)*um(p)) )
            cf  = 0.01_r8/(sqrt(uaf(p))*sqrt(dleaf(patch%itype(p))))
            rb(p)  = 1._r8/(cf*uaf(p))
            rb1(p) = rb(p)

            ! Parameterization for variation of csoilc with canopy density from
            ! X. Zeng, University of Arizona

            w = exp(-(elai(p)+esai(p)))

            ! changed by K.Sakaguchi from here
            ! transfer coefficient over bare soil is changed to a local variable
            ! just for readability of the code (from line 680)
            csoilb = (vkc/(0.13_r8*(z0mg(c)*uaf(p)/1.5e-5_r8)**0.45_r8))

            !compute the stability parameter for ricsoilc  ("S" in Sakaguchi&Zeng,2008)

            ri = ( grav*htop(p) * (taf(p) - t_grnd(c)) ) / (taf(p) * uaf(p) **2.00_r8)

            !! modify csoilc value (0.004) if the under-canopy is in stable condition

            if ( (taf(p) - t_grnd(c) ) > 0._r8) then
               ! decrease the value of csoilc by dividing it with (1+gamma*min(S, 10.0))
               ! ria ("gmanna" in Sakaguchi&Zeng, 2008) is a constant (=0.5)
               ricsoilc = csoilc / (1.00_r8 + ria*min( ri, 10.0_r8) )
               csoilcn = csoilb*w + ricsoilc*(1._r8-w)
            else
               csoilcn = csoilb*w + csoilc*(1._r8-w)
            end if

            !! Sakaguchi changes for stability formulation ends here

            rah(p,2) = 1._r8/(csoilcn*uaf(p))
            raw(p,2) = rah(p,2)
            if (use_lch4) then
               grnd_ch4_cond(p) = 1._r8/(raw(p,1)+raw(p,2))
            end if

            ! Stomatal resistances for sunlit and shaded fractions of canopy.
            ! Done each iteration to account for differences in eah, tv.

            svpts(p) = el(p)                         ! pa
            eah(p) = forc_pbot(c) * qaf(p) / 0.622_r8   ! pa
            rhaf(p) = eah(p)/svpts(p)
         end do

         ! Modification for shrubs proposed by X.D.Z 
         ! Equivalent modification for soy following AgroIBIS
         ! NOTE: the following block of code was moved out of Photosynthesis subroutine and 
         ! into here by M. Vertenstein on 4/6/2014 as part of making the photosynthesis
         ! routine a separate module. This move was also suggested by S. Levis in the previous
         ! version of the code. 
         ! BUG MV 4/7/2014 - is this the correct place to have it in the iteration? 
         ! THIS SHOULD BE MOVED OUT OF THE ITERATION but will change answers -
         
         do f = 1, fn
            p = filterp(f)
            c = patch%column(p)
            if (use_cndv) then
               if (patch%itype(p) == nbrdlf_dcd_tmp_shrub) then
                  btran(p) = min(1._r8, btran(p) * 3.33_r8)
               end if
            end if
            if (patch%itype(p) == nsoybean .or. patch%itype(p) == nsoybeanirrig) then
               btran(p) = min(1._r8, btran(p) * 1.25_r8)
            end if
         end do

         if ( use_ed ) then      

            call t_startf('edpsn')
            ! FIX(FIX(SPM,032414),032414) Photo*_ED will need refactoring
            call Photosynthesis_ED (bounds, fn, filterp, &
                 svpts(begp:endp), eah(begp:endp), o2(begp:endp), &
                 co2(begp:endp), rb(begp:endp), dayl_factor(begp:endp), &
                 ed_allsites_inst(begg:endg), atm2lnd_inst, temperature_inst, canopystate_inst, photosyns_inst)

            ! zero all of these things, not just the ones in the filter. 
            do p = bounds%begp,bounds%endp 
               photosyns_inst%rssun_patch(p)    = 0._r8
               photosyns_inst%rssha_patch(p)    = 0._r8
               photosyns_inst%psnsun_patch(p)   = 0._r8
               photosyns_inst%psnsha_patch(p)   = 0._r8
               photosyns_inst%fpsn_patch(p)     = 0._r8
               canopystate_inst%laisun_patch(p) = 0._r8
               canopystate_inst%laisha_patch(p) = 0._r8
            enddo

            call t_stopf('edpsn')

         else ! not use_ed

            call Photosynthesis (bounds, fn, filterp, &
                 svpts(begp:endp), eah(begp:endp), o2(begp:endp), co2(begp:endp), rb(begp:endp), btran(begp:endp), &
                 dayl_factor(begp:endp), atm2lnd_inst, temperature_inst, surfalb_inst, solarabs_inst, &
                 canopystate_inst, photosyns_inst, phase='sun')

            if ( use_cn .and. use_c13 ) then
               call Fractionation (bounds, fn, filterp, &
                    atm2lnd_inst, canopystate_inst, cnveg_state_inst, solarabs_inst, surfalb_inst, photosyns_inst, &
                    phase='sun')
            endif

            do f = 1, fn
               p = filterp(f)
               c = patch%column(p)
               if (use_cndv) then
                  if (patch%itype(p) == nbrdlf_dcd_tmp_shrub) then
                     btran(p) = min(1._r8, btran(p) * 3.33_r8)
                  end if
               end if
               if (patch%itype(p) == nsoybean .or. patch%itype(p) == nsoybeanirrig) then
                  btran(p) = min(1._r8, btran(p) * 1.25_r8)
               end if
            end do

            call Photosynthesis (bounds, fn, filterp, &
                 svpts(begp:endp), eah(begp:endp), o2(begp:endp), co2(begp:endp), rb(begp:endp), btran(begp:endp), &
                 dayl_factor(begp:endp), atm2lnd_inst, temperature_inst, surfalb_inst, solarabs_inst, &
                 canopystate_inst, photosyns_inst, phase='sha')

            if ( use_cn .and. use_c13 ) then
               call Fractionation (bounds, fn, filterp,  &
                    atm2lnd_inst, canopystate_inst, cnveg_state_inst, solarabs_inst, surfalb_inst, photosyns_inst, &
                    phase='sha')
            end if

         end if ! end of if use_ed

         do f = 1, fn
            p = filterp(f)
            c = patch%column(p)
            g = patch%gridcell(p)

            ! Sensible heat conductance for air, leaf and ground
            ! Moved the original subroutine in-line...

            wta    = 1._r8/rah(p,1)             ! air
            wtl    = (elai(p)+esai(p))/rb(p)    ! leaf
            wtg(p) = 1._r8/rah(p,2)             ! ground
            wtshi  = 1._r8/(wta+wtl+wtg(p))

            wtl0(p) = wtl*wtshi         ! leaf
            wtg0    = wtg(p)*wtshi      ! ground
            wta0(p) = wta*wtshi         ! air

            wtga    = wta0(p)+wtg0      ! ground + air
            wtal(p) = wta0(p)+wtl0(p)   ! air + leaf

            ! Fraction of potential evaporation from leaf

            if ( use_ed ) then
               
               if (fdry(p)  >  0._r8) then
                  rppdry  = fdry(p)*rb(p)/(rb(p)+rscanopy(p))
               else
                  rppdry = 0._r8
               end if
               if (use_lch4) then
                  ! Calculate canopy conductance for methane / oxygen (e.g. stomatal conductance & leaf bdy cond)
                  canopy_cond(p) = 1.0_r8/(rb(p)+rscanopy(p))
               end if

            else ! NOT use_ed

               if (fdry(p) > 0._r8) then
                  rppdry  = fdry(p)*rb(p)*(laisun(p)/(rb(p)+rssun(p)) + laisha(p)/(rb(p)+rssha(p)))/elai(p)
               else
                  rppdry = 0._r8
               end if

               ! Calculate canopy conductance for methane / oxygen (e.g. stomatal conductance & leaf bdy cond)
               if (use_lch4) then
                  canopy_cond(p) = (laisun(p)/(rb(p)+rssun(p)) + laisha(p)/(rb(p)+rssha(p)))/max(elai(p), 0.01_r8)
               end if

            end if ! end of if use_ed         

            efpot = forc_rho(c)*wtl*(qsatl(p)-qaf(p))

            if (efpot > 0._r8) then
               if (btran(p) > btran0) then
                  qflx_tran_veg(p) = efpot*rppdry
                  rpp = rppdry + fwet(p)
               else
                  !No transpiration if btran below 1.e-10
                  rpp = fwet(p)
                  qflx_tran_veg(p) = 0._r8
               end if
               !Check total evapotranspiration from leaves
               rpp = min(rpp, (qflx_tran_veg(p)+h2ocan(p)/dtime)/efpot)
            else
               !No transpiration if potential evaporation less than zero
               rpp = 1._r8
               qflx_tran_veg(p) = 0._r8
            end if

            ! Update conductances for changes in rpp
            ! Latent heat conductances for ground and leaf.
            ! Air has same conductance for both sensible and latent heat.
            ! Moved the original subroutine in-line...

            wtaq    = frac_veg_nosno(p)/raw(p,1)                        ! air
            wtlq    = frac_veg_nosno(p)*(elai(p)+esai(p))/rb(p) * rpp   ! leaf

            !Litter layer resistance. Added by K.Sakaguchi
            snow_depth_c = z_dl ! critical depth for 100% litter burial by snow (=litter thickness)
            fsno_dl = snow_depth(c)/snow_depth_c    ! effective snow cover for (dry)plant litter
            elai_dl = lai_dl*(1._r8 - min(fsno_dl,1._r8)) ! exposed (dry)litter area index
            rdl = ( 1._r8 - exp(-elai_dl) ) / ( 0.004_r8*uaf(p)) ! dry litter layer resistance

            ! add litter resistance and Lee and Pielke 1992 beta
            if (delq(p) < 0._r8) then  !dew. Do not apply beta for negative flux (follow old rsoil)
               wtgq(p) = frac_veg_nosno(p)/(raw(p,2)+rdl)
            else
               if (do_soilevap_beta()) then
                  wtgq(p) = soilbeta(c)*frac_veg_nosno(p)/(raw(p,2)+rdl)
               endif
            end if

            wtsqi   = 1._r8/(wtaq+wtlq+wtgq(p))

            wtgq0    = wtgq(p)*wtsqi      ! ground
            wtlq0(p) = wtlq*wtsqi         ! leaf
            wtaq0(p) = wtaq*wtsqi         ! air

            wtgaq    = wtaq0(p)+wtgq0     ! air + ground
            wtalq(p) = wtaq0(p)+wtlq0(p)  ! air + leaf

            dc1 = forc_rho(c)*cpair*wtl
            dc2 = hvap*forc_rho(c)*wtlq

            efsh   = dc1*(wtga*t_veg(p)-wtg0*t_grnd(c)-wta0(p)*thm(p))
            efe(p) = dc2*(wtgaq*qsatl(p)-wtgq0*qg(c)-wtaq0(p)*forc_q(c))

            ! Evaporation flux from foliage

            erre = 0._r8
            if (efe(p)*efeb(p) < 0._r8) then
               efeold = efe(p)
               efe(p)  = 0.1_r8*efeold
               erre = efe(p) - efeold
            end if
            ! fractionate ground emitted longwave
            lw_grnd=(frac_sno(c)*t_soisno(c,snl(c)+1)**4 &
                 +(1._r8-frac_sno(c)-frac_h2osfc(c))*t_soisno(c,1)**4 &
                 +frac_h2osfc(c)*t_h2osfc(c)**4)

            dt_veg(p) = (sabv(p) + air(p) + bir(p)*t_veg(p)**4 + &
                 cir(p)*lw_grnd - efsh - efe(p)) / &
                 (- 4._r8*bir(p)*t_veg(p)**3 +dc1*wtga +dc2*wtgaq*qsatldT(p))
            t_veg(p) = tlbef(p) + dt_veg(p)
            dels = dt_veg(p)
            del(p)  = abs(dels)
            err(p) = 0._r8
            if (del(p) > delmax) then
               dt_veg(p) = delmax*dels/del(p)
               t_veg(p) = tlbef(p) + dt_veg(p)
               err(p) = sabv(p) + air(p) + bir(p)*tlbef(p)**3*(tlbef(p) + &
                    4._r8*dt_veg(p)) + cir(p)*lw_grnd - &
                    (efsh + dc1*wtga*dt_veg(p)) - (efe(p) + &
                    dc2*wtgaq*qsatldT(p)*dt_veg(p))
            end if

            ! Fluxes from leaves to canopy space
            ! "efe" was limited as its sign changes frequently.  This limit may
            ! result in an imbalance in "hvap*qflx_evap_veg" and
            ! "efe + dc2*wtgaq*qsatdt_veg"

            efpot = forc_rho(c)*wtl*(wtgaq*(qsatl(p)+qsatldT(p)*dt_veg(p)) &
                 -wtgq0*qg(c)-wtaq0(p)*forc_q(c))
            qflx_evap_veg(p) = rpp*efpot

            ! Calculation of evaporative potentials (efpot) and
            ! interception losses; flux in kg m**-2 s-1.  ecidif
            ! holds the excess energy if all intercepted water is evaporated
            ! during the timestep.  This energy is later added to the
            ! sensible heat flux.

            ecidif = 0._r8
            if (efpot > 0._r8 .and. btran(p) > btran0) then
               qflx_tran_veg(p) = efpot*rppdry
            else
               qflx_tran_veg(p) = 0._r8
            end if
            ecidif = max(0._r8, qflx_evap_veg(p)-qflx_tran_veg(p)-h2ocan(p)/dtime)
            qflx_evap_veg(p) = min(qflx_evap_veg(p),qflx_tran_veg(p)+h2ocan(p)/dtime)

            ! The energy loss due to above two limits is added to
            ! the sensible heat flux.

            eflx_sh_veg(p) = efsh + dc1*wtga*dt_veg(p) + err(p) + erre + hvap*ecidif

            ! Re-calculate saturated vapor pressure, specific humidity, and their
            ! derivatives at the leaf surface

            call QSat(t_veg(p), forc_pbot(c), el(p), deldT, qsatl(p), qsatldT(p))

            ! Update vegetation/ground surface temperature, canopy air
            ! temperature, canopy vapor pressure, aerodynamic temperature, and
            ! Monin-Obukhov stability parameter for next iteration.

            taf(p) = wtg0*t_grnd(c) + wta0(p)*thm(p) + wtl0(p)*t_veg(p)
            qaf(p) = wtlq0(p)*qsatl(p) + wtgq0*qg(c) + forc_q(c)*wtaq0(p)

            ! Update Monin-Obukhov length and wind speed including the
            ! stability effect

            dth(p) = thm(p)-taf(p)
            dqh(p) = forc_q(c)-qaf(p)
            delq(p) = wtalq(p)*qg(c)-wtlq0(p)*qsatl(p)-wtaq0(p)*forc_q(c)

            tstar = temp1(p)*dth(p)
            qstar = temp2(p)*dqh(p)

            thvstar = tstar*(1._r8+0.61_r8*forc_q(c)) + 0.61_r8*forc_th(c)*qstar
            zeta = zldis(p)*vkc*grav*thvstar/(ustar(p)**2*thv(c))

            if (zeta >= 0._r8) then     !stable
               zeta = min(2._r8,max(zeta,0.01_r8))
               um(p) = max(ur(p),0.1_r8)
            else                     !unstable
               zeta = max(-100._r8,min(zeta,-0.01_r8))
               wc = beta*(-grav*ustar(p)*thvstar*zii/thv(c))**0.333_r8
               um(p) = sqrt(ur(p)*ur(p)+wc*wc)
            end if
            obu(p) = zldis(p)/zeta

            if (obuold(p)*obu(p) < 0._r8) nmozsgn(p) = nmozsgn(p)+1
            if (nmozsgn(p) >= 4) obu(p) = zldis(p)/(-0.01_r8)
            obuold(p) = obu(p)

         end do   ! end of filtered patch loop

         ! Test for convergence

         itlef = itlef+1
         if (itlef > itmin) then
            do f = 1, fn
               p = filterp(f)
               dele(p) = abs(efe(p)-efeb(p))
               efeb(p) = efe(p)
               det(p)  = max(del(p),del2(p))
            end do
            fnold = fn
            fn = 0
            do f = 1, fnold
               p = filterp(f)
               if (.not. (det(p) < dtmin .and. dele(p) < dlemin)) then
                  fn = fn + 1
                  filterp(fn) = p
               end if
            end do
         end if

      end do ITERATION     ! End stability iteration
      call t_stopf('can_iter')

      fn = fnorig
      filterp(1:fn) = fporig(1:fn)

      do f = 1, fn
         p = filterp(f)
         c = patch%column(p)
         g = patch%gridcell(p)

         ! Energy balance check in canopy

         lw_grnd=(frac_sno(c)*t_soisno(c,snl(c)+1)**4 &
              +(1._r8-frac_sno(c)-frac_h2osfc(c))*t_soisno(c,1)**4 &
              +frac_h2osfc(c)*t_h2osfc(c)**4)

         err(p) = sabv(p) + air(p) + bir(p)*tlbef(p)**3*(tlbef(p) + 4._r8*dt_veg(p)) &
                                !+ cir(p)*t_grnd(c)**4 - eflx_sh_veg(p) - hvap*qflx_evap_veg(p)
              + cir(p)*lw_grnd - eflx_sh_veg(p) - hvap*qflx_evap_veg(p)

         ! Fluxes from ground to canopy space

         delt    = wtal(p)*t_grnd(c)-wtl0(p)*t_veg(p)-wta0(p)*thm(p)
         taux(p) = -forc_rho(c)*forc_u(g)/ram1(p)
         tauy(p) = -forc_rho(c)*forc_v(g)/ram1(p)
         eflx_sh_grnd(p) = cpair*forc_rho(c)*wtg(p)*delt

         ! compute individual sensible heat fluxes
         delt_snow = wtal(p)*t_soisno(c,snl(c)+1)-wtl0(p)*t_veg(p)-wta0(p)*thm(p)
         eflx_sh_snow(p) = cpair*forc_rho(c)*wtg(p)*delt_snow

         delt_soil  = wtal(p)*t_soisno(c,1)-wtl0(p)*t_veg(p)-wta0(p)*thm(p)
         eflx_sh_soil(p) = cpair*forc_rho(c)*wtg(p)*delt_soil

         delt_h2osfc  = wtal(p)*t_h2osfc(c)-wtl0(p)*t_veg(p)-wta0(p)*thm(p)
         eflx_sh_h2osfc(p) = cpair*forc_rho(c)*wtg(p)*delt_h2osfc
         qflx_evap_soi(p) = forc_rho(c)*wtgq(p)*delq(p)

         ! compute individual latent heat fluxes
         delq_snow = wtalq(p)*qg_snow(c)-wtlq0(p)*qsatl(p)-wtaq0(p)*forc_q(c)
         qflx_ev_snow(p) = forc_rho(c)*wtgq(p)*delq_snow

         delq_soil = wtalq(p)*qg_soil(c)-wtlq0(p)*qsatl(p)-wtaq0(p)*forc_q(c)
         qflx_ev_soil(p) = forc_rho(c)*wtgq(p)*delq_soil

         delq_h2osfc = wtalq(p)*qg_h2osfc(c)-wtlq0(p)*qsatl(p)-wtaq0(p)*forc_q(c)
         qflx_ev_h2osfc(p) = forc_rho(c)*wtgq(p)*delq_h2osfc

         ! 2 m height air temperature

         t_ref2m(p) = thm(p) + temp1(p)*dth(p)*(1._r8/temp12m(p) - 1._r8/temp1(p))
         t_ref2m_r(p) = t_ref2m(p)

         ! 2 m height specific humidity

         q_ref2m(p) = forc_q(c) + temp2(p)*dqh(p)*(1._r8/temp22m(p) - 1._r8/temp2(p))

         ! 2 m height relative humidity

         call QSat(t_ref2m(p), forc_pbot(c), e_ref2m, de2mdT, qsat_ref2m, dqsat2mdT)
         rh_ref2m(p) = min(100._r8, q_ref2m(p) / qsat_ref2m * 100._r8)
         rh_ref2m_r(p) = rh_ref2m(p)

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

         ! Downward longwave radiation below the canopy

         dlrad(p) = (1._r8-emv(p))*emg(c)*forc_lwrad(c) + &
              emv(p)*emg(c)*sb*tlbef(p)**3*(tlbef(p) + 4._r8*dt_veg(p))

         ! Upward longwave radiation above the canopy

         ulrad(p) = ((1._r8-emg(c))*(1._r8-emv(p))*(1._r8-emv(p))*forc_lwrad(c) &
              + emv(p)*(1._r8+(1._r8-emg(c))*(1._r8-emv(p)))*sb*tlbef(p)**3*(tlbef(p) + &
              4._r8*dt_veg(p)) + emg(c)*(1._r8-emv(p))*sb*lw_grnd)

         ! Derivative of soil energy flux with respect to soil temperature

         cgrnds(p) = cgrnds(p) + cpair*forc_rho(c)*wtg(p)*wtal(p)
         cgrndl(p) = cgrndl(p) + forc_rho(c)*wtgq(p)*wtalq(p)*dqgdT(c)
         cgrnd(p)  = cgrnds(p) + cgrndl(p)*htvp(c)

         ! Update dew accumulation (kg/m2)

         h2ocan(p) = max(0._r8,h2ocan(p)+(qflx_tran_veg(p)-qflx_evap_veg(p))*dtime)

         if ( use_ed ) then
            call AccumulateFluxes_ED(bounds, p, ed_allsites_inst(begg:endg), photosyns_inst)
         end if

      end do

      ! Determine total photosynthesis

      call PhotosynthesisTotal(fn, filterp, &
           atm2lnd_inst, canopystate_inst, photosyns_inst)

      ! Filter out patches which have small energy balance errors; report others

      fnold = fn
      fn = 0
      do f = 1, fnold
         p = filterp(f)
         if (abs(err(p)) > 0.1_r8) then
            fn = fn + 1
            filterp(fn) = p
         end if
      end do

      do f = 1, fn
         p = filterp(f)
         write(iulog,*) 'energy balance in canopy ',p,', err=',err(p)
      end do

      if ( use_ed ) then      
         ! zero all of the array,  not just the ones in the filter. 
         do p = bounds%begp,bounds%endp 
            photosyns_inst%rssun_patch(p)    = 0._r8
            photosyns_inst%rssha_patch(p)    = 0._r8
            photosyns_inst%psnsun_patch(p)   = 0._r8
            photosyns_inst%psnsha_patch(p)   = 0._r8
            photosyns_inst%fpsn_patch(p)     = 0._r8
            canopystate_inst%laisun_patch(p) = 0._r8
            canopystate_inst%laisha_patch(p) = 0._r8
         enddo
      end if

    end associate


  end subroutine CanopyFluxes

end module CanopyFluxesMod

