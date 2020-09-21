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
  use elm_varctl            , only : iulog, use_cn, use_lch4, use_c13, use_c14, use_fates
  use elm_varctl            , only : use_hydrstress
  use elm_varpar            , only : nlevgrnd, nlevsno
  use elm_varcon            , only : namep 
  use pftvarcon             , only : nbrdlf_dcd_tmp_shrub, nsoybean , nsoybeanirrig
  use decompMod             , only : bounds_type
  use PhotosynthesisMod     , only : Photosynthesis, PhotosynthesisTotal, Fractionation, PhotoSynthesisHydraulicStress
  use SoilMoistStressMod    , only : calc_effective_soilporosity, calc_volumetric_h2oliq
  use SoilMoistStressMod    , only : calc_root_moist_stress, set_perchroot_opt
  use SimpleMathMod         , only : array_div_vector
  use SurfaceResistanceMod  , only : do_soilevap_beta
  use VegetationPropertiesType        , only : veg_vp
  use atm2lndType           , only : atm2lnd_type
  use CanopyStateType       , only : canopystate_type
  use CNStateType           , only : cnstate_type
  use EnergyFluxType        , only : energyflux_type
  use FrictionvelocityType  , only : frictionvel_type
  use SoilStateType         , only : soilstate_type
  use SolarAbsorbedType     , only : solarabs_type
  use SurfaceAlbedoType     , only : surfalb_type
  use TemperatureType       , only : temperature_type
  use WaterfluxType         , only : waterflux_type
  use WaterstateType        , only : waterstate_type
  use CH4Mod                , only : ch4_type
  use PhotosynthesisType    , only : photosyns_type
  use PhosphorusStateType   , only : phosphorusstate_type
  use CNNitrogenStateType   , only : nitrogenstate_type
  use CLMFatesInterfaceMod  , only : hlm_fates_interface_type
  use GridcellType          , only : grc_pp 
  use TopounitDataType      , only : top_as, top_af  
  use ColumnType            , only : col_pp
  use ColumnDataType        , only : col_es, col_ef, col_ws               
  use VegetationType        , only : veg_pp                
  use VegetationDataType    , only : veg_es, veg_ef, veg_ws, veg_wf  
  !
  ! !PUBLIC TYPES:
  implicit none
  save
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
  subroutine CanopyFluxes(bounds,  num_nolakeurbanp, filter_nolakeurbanp, &
       atm2lnd_vars, canopystate_vars, cnstate_vars, energyflux_vars, &
       frictionvel_vars, soilstate_vars, solarabs_vars, surfalb_vars, &
       temperature_vars, waterflux_vars, waterstate_vars, ch4_vars, photosyns_vars, &
       soil_water_retention_curve, nitrogenstate_vars, phosphorusstate_vars, &
       alm_fates) 
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
    use shr_const_mod      , only : SHR_CONST_TKFRZ, SHR_CONST_RGAS
    use clm_time_manager   , only : get_step_size, get_prev_date, get_nstep
    use elm_varcon         , only : sb, cpair, hvap, vkc, grav, denice
    use elm_varcon         , only : denh2o, tfrz, csoilc, tlsai_crit, alpha_aero
    use elm_varcon         , only : isecspday, degpsec
    use pftvarcon          , only : irrigated
    use elm_varcon         , only : c14ratio
    use perf_mod           , only : t_startf, t_stopf
    use domainMod          , only : ldomain
    use QSatMod            , only : QSat
    use FrictionVelocityMod, only : FrictionVelocity, MoninObukIni
    use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
    use SurfaceResistanceMod, only : getlblcef
    !
    ! !ARGUMENTS:
    type(bounds_type)         , intent(in)    :: bounds 
    integer                   , intent(in)    :: num_nolakeurbanp       ! number of column non-lake, non-urban points in pft filter
    integer                   , intent(in)    :: filter_nolakeurbanp(:) ! patch filter for non-lake, non-urban points
    type(atm2lnd_type)        , intent(in)    :: atm2lnd_vars
    type(canopystate_type)    , intent(inout) :: canopystate_vars
    type(cnstate_type)        , intent(in)    :: cnstate_vars
    type(energyflux_type)     , intent(inout) :: energyflux_vars
    type(frictionvel_type)    , intent(inout) :: frictionvel_vars
    type(solarabs_type)       , intent(in)    :: solarabs_vars
    type(surfalb_type)        , intent(in)    :: surfalb_vars
    type(soilstate_type)      , intent(inout) :: soilstate_vars
    type(temperature_type)    , intent(inout) :: temperature_vars
    type(waterstate_type)     , intent(inout) :: waterstate_vars
    type(waterflux_type)      , intent(inout) :: waterflux_vars
    type(ch4_type)            , intent(inout) :: ch4_vars
    type(photosyns_type)      , intent(inout) :: photosyns_vars
    class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve
    type(nitrogenstate_type)  , intent(inout) :: nitrogenstate_vars
    type(phosphorusstate_type), intent(inout) :: phosphorusstate_vars
    type(hlm_fates_interface_type) , intent(inout) :: alm_fates
    !
    ! !LOCAL VARIABLES:
    real(r8), pointer   :: bsun(:)          ! sunlit canopy transpiration wetness factor (0 to 1)
    real(r8), pointer   :: bsha(:)          ! shaded canopy transpiration wetness factor (0 to 1)
    real(r8), parameter :: btran0 = 0.0_r8  ! initial value
    real(r8), parameter :: zii = 1000.0_r8  ! convective boundary layer height [m]
    real(r8), parameter :: beta = 1.0_r8    ! coefficient of conective velocity [-]
    real(r8), parameter :: delmax = 1.0_r8  ! maxchange in  leaf temperature [K]
    real(r8), parameter :: dlemin = 0.1_r8  ! max limit for energy flux convergence [w/m2]
    real(r8), parameter :: dtmin = 0.01_r8  ! max limit for temperature convergence [K]
    integer , parameter :: itmax = 40       ! maximum number of iteration [-]
    integer , parameter :: itmin = 2        ! minimum number of iteration [-]
    real(r8), parameter :: irrig_min_lai = 0.0_r8           ! Minimum LAI for irrigation
    real(r8), parameter :: irrig_btran_thresh = 0.999999_r8 ! Irrigate when btran falls below 0.999999 rather than 1 to allow for round-off error
    integer , parameter :: irrig_start_time = isecspday/4   ! (6AM) Time of day to check whether we need irrigation, seconds (0 = midnight). 

    ! We start applying the irrigation in the time step FOLLOWING this time, 
    ! since we won't begin irrigating until the next call to CanopyHydrology
    ! Desired amount of time to irrigate per day (sec). Actual time may 
    ! differ if this is not a multiple of dtime. Irrigation won't work properly 
    ! if dtime > secsperday
    integer , parameter :: irrig_length = isecspday/6     ! 4 hours irrigation

    ! Determines target soil moisture level for irrigation. If h2osoi_liq_so 
    ! is the soil moisture level at which stomata are fully open and 
    ! h2osoi_liq_sat is the soil moisture level at saturation (eff_porosity), 
    ! then the target soil moisture level is 
    !     (h2osoi_liq_so + irrig_factor*(h2osoi_liq_sat - h2osoi_liq_so)). 
    ! A value of 0 means that the target soil moisture level is h2osoi_liq_so. 
    ! A value of 1 means that the target soil moisture level is h2osoi_liq_sat
    real(r8), parameter :: irrig_factor = 0.7_r8            

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
    real(r8) :: cf_bare                              ! heat transfer coefficient from bare ground [-]
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
    integer  :: t                                    ! topounit index
    integer  :: g                                    ! gridcell index
    integer  :: fp                                   ! lake filter pft index
    integer  :: fn_noveg                             ! number of values in bare ground pft filter
    integer  :: filterp_noveg(bounds%endp-bounds%begp+1) ! bare ground pft filter
    integer  :: fn                                   ! number of values in vegetated pft filter
    integer  :: filterp(bounds%endp-bounds%begp+1)   ! vegetated pft filter
    integer  :: fnorig                               ! number of values in pft filter copy
    integer  :: fporig(bounds%endp-bounds%begp+1)    ! temporary filter
    integer  :: fnold                                ! temporary copy of pft count
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
    integer  :: yr                                       ! year at start of time step
    integer  :: mon                                      ! month at start of time step
    integer  :: day                                      ! day at start of time step
    integer  :: time                                     ! time at start of time step (seconds after 0Z)
    integer  :: local_time                               ! local time at start of time step (seconds after solar midnight)
    integer  :: seconds_since_irrig_start_time
    integer  :: irrig_nsteps_per_day                     ! number of time steps per day in which we irrigate
    logical  :: check_for_irrig(bounds%begp:bounds%endp) ! where do we need to check soil moisture to see if we need to irrigate?
    logical  :: frozen_soil(bounds%begp:bounds%endp)     ! set to true if we have encountered a frozen soil layer
    real(r8) :: vol_liq_so                               ! partial volume of liquid water in layer for which smp_node = smpso
    real(r8) :: h2osoi_liq_so                            ! liquid water corresponding to vol_liq_so for this layer [kg/m2]
    real(r8) :: h2osoi_liq_sat                           ! liquid water corresponding to eff_porosity for this layer [kg/m2]
    real(r8) :: deficit                                  ! difference between desired soil moisture level for this layer and 
                                                         ! current soil moisture level [kg/m2]
    real(r8) :: dt_veg(bounds%begp:bounds%endp)          ! change in t_veg, last iteration (Kelvin)                              
    integer  :: jtop(bounds%begc:bounds%endc)            ! lbning
    integer  :: filterc_tmp(bounds%endp-bounds%begp+1)   ! temporary variable
    integer  :: ft                                       ! plant functional type index
    real(r8) :: temprootr                 
    real(r8) :: dt_veg_temp(bounds%begp:bounds%endp)
    integer  :: iv
    !------------------------------------------------------------------------------

    associate(                                                               & 
         snl                  => col_pp%snl                                   , & ! Input:  [integer  (:)   ]  number of snow layers                                                  
         dayl                 => grc_pp%dayl                                  , & ! Input:  [real(r8) (:)   ]  daylength (s)
         max_dayl             => grc_pp%max_dayl                              , & ! Input:  [real(r8) (:)   ]  maximum daylength for this grid cell (s)

         forc_lwrad           => top_af%lwrad                              , & ! Input:  [real(r8) (:)   ]  downward infrared (longwave) radiation (W/m**2)                       
         forc_q               => top_as%qbot                               , & ! Input:  [real(r8) (:)   ]  atmospheric specific humidity (kg/kg)                                 
         forc_pbot            => top_as%pbot                               , & ! Input:  [real(r8) (:)   ]  atmospheric pressure (Pa)                                             
         forc_th              => top_as%thbot                              , & ! Input:  [real(r8) (:)   ]  atmospheric potential temperature (Kelvin)                            
         forc_rho             => top_as%rhobot                             , & ! Input:  [real(r8) (:)   ]  air density (kg/m**3)                                                     
         forc_t               => top_as%tbot                               , & ! Input:  [real(r8) (:)   ]  atmospheric temperature (Kelvin)                                      
         forc_u               => top_as%ubot                               , & ! Input:  [real(r8) (:)   ]  atmospheric wind speed in east direction (m/s)                        
         forc_v               => top_as%vbot                               , & ! Input:  [real(r8) (:)   ]  atmospheric wind speed in north direction (m/s)                       
         forc_pco2            => top_as%pco2bot                            , & ! Input:  [real(r8) (:)   ]  partial pressure co2 (Pa)                                             
         forc_pc13o2          => top_as%pc13o2bot                          , & ! Input:  [real(r8) (:)   ]  partial pressure c13o2 (Pa)                                           
         forc_po2             => top_as%po2bot                             , & ! Input:  [real(r8) (:)   ]  partial pressure o2 (Pa)                                              

         dleaf                => veg_vp%dleaf                          , & ! Input:  [real(r8) (:)   ]  characteristic leaf dimension (m)                                     
         smpso                => veg_vp%smpso                          , & ! Input:  [real(r8) (:)   ]  soil water potential at full stomatal opening (mm)                    
         smpsc                => veg_vp%smpsc                          , & ! Input:  [real(r8) (:)   ]  soil water potential at full stomatal closure (mm)                    

         htvp                 => col_ef%htvp                  , & ! Input:  [real(r8) (:)   ]  latent heat of evaporation (/sublimation) [J/kg] (constant)                      

         sabv                 => solarabs_vars%sabv_patch                  , & ! Input:  [real(r8) (:)   ]  solar radiation absorbed by vegetation (W/m**2)                       

         lbl_rsc_h2o          => canopystate_vars%lbl_rsc_h2o_patch        , & ! Output: [real(r8) (:)   ] laminar boundary layer resistance for h2o
         frac_veg_nosno       => canopystate_vars%frac_veg_nosno_patch     , & ! Input:  [integer  (:)   ]  fraction of vegetation not covered by snow (0 OR 1) [-]
         elai                 => canopystate_vars%elai_patch               , & ! Input:  [real(r8) (:)   ]  one-sided leaf area index with burying by snow                        
         esai                 => canopystate_vars%esai_patch               , & ! Input:  [real(r8) (:)   ]  one-sided stem area index with burying by snow                        
         laisun               => canopystate_vars%laisun_patch             , & ! Input:  [real(r8) (:)   ]  sunlit leaf area                                                      
         laisha               => canopystate_vars%laisha_patch             , & ! Input:  [real(r8) (:)   ]  shaded leaf area                                                      
         displa               => canopystate_vars%displa_patch             , & ! Input:  [real(r8) (:)   ]  displacement height (m)                                               
         htop                 => canopystate_vars%htop_patch               , & ! Input:  [real(r8) (:)   ]  canopy top(m)                                                         
         altmax_lastyear_indx => canopystate_vars%altmax_lastyear_indx_col , & ! Input:  [integer  (:)   ]  prior year maximum annual depth of thaw                                
         altmax_indx          => canopystate_vars%altmax_indx_col          , & ! Input:  [integer  (:)   ]  maximum annual depth of thaw                                           

         dleaf_patch          => canopystate_vars%dleaf_patch                 , & ! Output: [real(r8) (:)   ]  mean leaf diameter for this patch/pft
         watsat               => soilstate_vars%watsat_col                 , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)   (constant)                     
         watdry               => soilstate_vars%watdry_col                 , & ! Input:  [real(r8) (:,:) ]  btran parameter for btran=0                      (constant)                                        
         watopt               => soilstate_vars%watopt_col                 , & ! Input:  [real(r8) (:,:) ]  btran parameter for btran=1                      (constant)                                      
         eff_porosity         => soilstate_vars%eff_porosity_col           , & ! Output: [real(r8) (:,:) ]  effective soil porosity

         sucsat               => soilstate_vars%sucsat_col                 , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)                        (constant)                                        
         bsw                  => soilstate_vars%bsw_col                    , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"                         (constant)                                        
         rootfr               => soilstate_vars%rootfr_patch               , & ! Input:  [real(r8) (:,:) ]  fraction of roots in each soil layer                                
         soilbeta             => soilstate_vars%soilbeta_col               , & ! Input:  [real(r8) (:)   ]  soil wetness relative to field capacity                               
         rootr                => soilstate_vars%rootr_patch                , & ! Output: [real(r8) (:,:) ]  effective fraction of roots in each soil layer                      

         forc_hgt_u_patch     => frictionvel_vars%forc_hgt_u_patch         , & ! Input:  [real(r8) (:)   ]  observational height of wind at pft level [m]                          
         z0mg                 => frictionvel_vars%z0mg_col                 , & ! Input:  [real(r8) (:)   ]  roughness length of ground, momentum [m]                              
         ram1                 => frictionvel_vars%ram1_patch               , & ! Output: [real(r8) (:)   ]  aerodynamical resistance (s/m)                                        
         z0mv                 => frictionvel_vars%z0mv_patch               , & ! Output: [real(r8) (:)   ]  roughness length over vegetation, momentum [m]                        
         z0hv                 => frictionvel_vars%z0hv_patch               , & ! Output: [real(r8) (:)   ]  roughness length over vegetation, sensible heat [m]                   
         z0qv                 => frictionvel_vars%z0qv_patch               , & ! Output: [real(r8) (:)   ]  roughness length over vegetation, latent heat [m]                     
         rb1                  => frictionvel_vars%rb1_patch                , & ! Output: [real(r8) (:)   ]  boundary layer resistance (s/m)                                       

         t_h2osfc             => col_es%t_h2osfc             , & ! Input:  [real(r8) (:)   ]  surface water temperature                                             
         t_soisno             => col_es%t_soisno             , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)                                           
         t_grnd               => col_es%t_grnd               , & ! Input:  [real(r8) (:)   ]  ground surface temperature [K]                                        
         thv                  => col_es%thv                  , & ! Input:  [real(r8) (:)   ]  virtual potential temperature (kelvin)                                
         thm                  => veg_es%thm                  , & ! Input:  [real(r8) (:)   ]  intermediate variable (forc_t+0.0098*forc_hgt_t_patch)                  
         emv                  => veg_es%emv                  , & ! Input:  [real(r8) (:)   ]  vegetation emissivity                                                     
         emg                  => col_es%emg                  , & ! Input:  [real(r8) (:)   ]  vegetation emissivity                                                 
         t_veg                => veg_es%t_veg                , & ! Output: [real(r8) (:)   ]  vegetation temperature (Kelvin)                                       
         t_ref2m              => veg_es%t_ref2m              , & ! Output: [real(r8) (:)   ]  2 m height surface air temperature (Kelvin)                           
         t_ref2m_r            => veg_es%t_ref2m_r            , & ! Output: [real(r8) (:)   ]  Rural 2 m height surface air temperature (Kelvin)                     

         frac_h2osfc          => col_ws%frac_h2osfc           , & ! Input:  [real(r8) (:)   ]  fraction of surface water                                             
         fwet                 => veg_ws%fwet                , & ! Input:  [real(r8) (:)   ]  fraction of canopy that is wet (0 to 1)                               
         fdry                 => veg_ws%fdry                , & ! Input:  [real(r8) (:)   ]  fraction of foliage that is green and dry [-]                         
         frac_sno             => col_ws%frac_sno_eff          , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)                           
         snow_depth           => col_ws%snow_depth            , & ! Input:  [real(r8) (:)   ]  snow height (m)                                                       
         qg_snow              => col_ws%qg_snow               , & ! Input:  [real(r8) (:)   ]  specific humidity at snow surface [kg/kg]                             
         qg_soil              => col_ws%qg_soil               , & ! Input:  [real(r8) (:)   ]  specific humidity at soil surface [kg/kg]                             
         qg_h2osfc            => col_ws%qg_h2osfc             , & ! Input:  [real(r8) (:)   ]  specific humidity at h2osfc surface [kg/kg]                           
         qg                   => col_ws%qg                    , & ! Input:  [real(r8) (:)   ]  specific humidity at ground surface [kg/kg]                           
         dqgdT                => col_ws%dqgdT                 , & ! Input:  [real(r8) (:)   ]  temperature derivative of "qg"                                        
         h2osoi_ice           => col_ws%h2osoi_ice            , & ! Input:  [real(r8) (:,:) ]  ice lens (kg/m2)                                                    
         h2osoi_vol           => col_ws%h2osoi_vol            , & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3] by F. Li and S. Levis
         h2osoi_liq           => col_ws%h2osoi_liq            , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)                                                
         h2osoi_liqvol        => col_ws%h2osoi_liqvol         , & ! Output: [real(r8) (:,:) ]  volumetric liquid water (v/v) 

         h2ocan               => veg_ws%h2ocan              , & ! Output: [real(r8) (:)   ]  canopy water (mm H2O)                                                 
         q_ref2m              => veg_ws%q_ref2m             , & ! Output: [real(r8) (:)   ]  2 m height surface specific humidity (kg/kg)                          
         rh_ref2m_r           => veg_ws%rh_ref2m_r          , & ! Output: [real(r8) (:)   ]  Rural 2 m height surface relative humidity (%)                        
         rh_ref2m             => veg_ws%rh_ref2m            , & ! Output: [real(r8) (:)   ]  2 m height surface relative humidity (%)                              
         rhaf                 => veg_ws%rh_af               , & ! Output: [real(r8) (:)   ]  fractional humidity of canopy air [dimensionless]                     

         pgwgt                => veg_pp%wtgcell              , & ! Input:  [integer  (:)   ]  pft's weight in gridcell
         n_irrig_steps_left   => veg_wf%n_irrig_steps_left   , & ! Output: [integer  (:)   ]  number of time steps for which we still need to irrigate today              
         irrig_rate           => veg_wf%irrig_rate           , & ! Output: [real(r8) (:)   ]  current irrigation rate [mm/s]                                        
         qflx_tran_veg        => veg_wf%qflx_tran_veg        , & ! Output: [real(r8) (:)   ]  vegetation transpiration (mm H2O/s) (+ = to atm)                      
         qflx_evap_veg        => veg_wf%qflx_evap_veg        , & ! Output: [real(r8) (:)   ]  vegetation evaporation (mm H2O/s) (+ = to atm)                        
         qflx_evap_soi        => veg_wf%qflx_evap_soi        , & ! Output: [real(r8) (:)   ]  soil evaporation (mm H2O/s) (+ = to atm)                              
         qflx_ev_snow         => veg_wf%qflx_ev_snow         , & ! Output: [real(r8) (:)   ]  evaporation flux from snow (W/m**2) [+ to atm]                        
         qflx_ev_soil         => veg_wf%qflx_ev_soil         , & ! Output: [real(r8) (:)   ]  evaporation flux from soil (W/m**2) [+ to atm]                        
         qflx_ev_h2osfc       => veg_wf%qflx_ev_h2osfc       , & ! Output: [real(r8) (:)   ]  evaporation flux from h2osfc (W/m**2) [+ to atm]                      

         rssun                => photosyns_vars%rssun_patch                , & ! Output: [real(r8) (:)   ]  leaf sunlit stomatal resistance (s/m) (output from Photosynthesis)
         rssha                => photosyns_vars%rssha_patch                , & ! Output: [real(r8) (:)   ]  leaf shaded stomatal resistance (s/m) (output from Photosynthesis)

         grnd_ch4_cond        => ch4_vars%grnd_ch4_cond_patch              , & ! Output: [real(r8) (:)   ]  tracer conductance for boundary layer [m/s] 

         btran2               => energyflux_vars%btran2_patch              , & ! Output: [real(r8) (:)   ]  F. Li and S. Levis                                                     
         btran                => energyflux_vars%btran_patch               , & ! Output: [real(r8) (:)   ]  transpiration wetness factor (0 to 1)                                 
         rresis               => energyflux_vars%rresis_patch              , & ! Output: [real(r8) (:,:) ]  root resistance by layer (0-1)  (nlevgrnd)                          
         taux                 => veg_ef%taux                , & ! Output: [real(r8) (:)   ]  wind (shear) stress: e-w (kg/m/s**2)                                  
         tauy                 => veg_ef%tauy                , & ! Output: [real(r8) (:)   ]  wind (shear) stress: n-s (kg/m/s**2)                                  
         canopy_cond          => energyflux_vars%canopy_cond_patch         , & ! Output: [real(r8) (:)   ]  tracer conductance for canopy [m/s] 
         cgrnds               => veg_ef%cgrnds              , & ! Output: [real(r8) (:)   ]  deriv. of soil sensible heat flux wrt soil temp [w/m2/k]              
         cgrndl               => veg_ef%cgrndl              , & ! Output: [real(r8) (:)   ]  deriv. of soil latent heat flux wrt soil temp [w/m**2/k]              
         dlrad                => veg_ef%dlrad               , & ! Output: [real(r8) (:)   ]  downward longwave radiation below the canopy [W/m2]                   
         ulrad                => veg_ef%ulrad               , & ! Output: [real(r8) (:)   ]  upward longwave radiation above the canopy [W/m2]                     
         cgrnd                => veg_ef%cgrnd               , & ! Output: [real(r8) (:)   ]  deriv. of soil energy flux wrt to soil temp [w/m2/k]                  
         eflx_sh_snow         => veg_ef%eflx_sh_snow        , & ! Output: [real(r8) (:)   ]  sensible heat flux from snow (W/m**2) [+ to atm]                      
         eflx_sh_h2osfc       => veg_ef%eflx_sh_h2osfc      , & ! Output: [real(r8) (:)   ]  sensible heat flux from soil (W/m**2) [+ to atm]                      
         eflx_sh_soil         => veg_ef%eflx_sh_soil        , & ! Output: [real(r8) (:)   ]  sensible heat flux from soil (W/m**2) [+ to atm]                      
         eflx_sh_veg          => veg_ef%eflx_sh_veg         , & ! Output: [real(r8) (:)   ]  sensible heat flux from leaves (W/m**2) [+ to atm]                    
         eflx_sh_grnd         => veg_ef%eflx_sh_grnd        , & ! Output: [real(r8) (:)   ]  sensible heat flux from ground (W/m**2) [+ to atm]
         begp                 => bounds%begp                               , &
         endp                 => bounds%endp                                 &
         )
      if (use_hydrstress) then
        bsun                    => energyflux_vars%bsun_patch ! Output:[real(r8) (:)   ]  sunlit canopy transpiration wetness factor (0 to 1)
        bsha                    => energyflux_vars%bsha_patch ! Output:[real(r8) (:)   ]  sunlit canopy transpiration wetness factor (0 to 1)
      end if
      ! Determine step size

      dtime = get_step_size()
      irrig_nsteps_per_day = ((irrig_length + (dtime - 1))/dtime)  ! round up
      ! irrig_nsteps_per_day = 12 !6 hrs
      ! First - set the following values over points where frac vegetation covered by snow is zero
      ! (e.g. btran, t_veg, rootr, rresis)

      do fp = 1,num_nolakeurbanp
         p = filter_nolakeurbanp(fp)
         c = veg_pp%column(p)
         t = veg_pp%topounit(p)
         if (frac_veg_nosno(p) == 0) then
            btran(p) = 0._r8     
            t_veg(p) = forc_t(t) 
            cf_bare  = forc_pbot(t)/(SHR_CONST_RGAS*0.001_r8*thm(p))*1.e06_r8
            rssun(p) = 1._r8/1.e15_r8 * cf_bare
            rssha(p) = 1._r8/1.e15_r8 * cf_bare
            lbl_rsc_h2o(p)=0._r8
            do j = 1, nlevgrnd
               rootr(p,j)  = 0._r8
               rresis(p,j) = 0._r8
            end do
         end if
      end do

      ! -----------------------------------------------------------------
      ! Time step initialization of photosynthesis variables
      ! -----------------------------------------------------------------

      call photosyns_vars%TimeStepInit(bounds)


      ! -----------------------------------------------------------------
      ! Filter patches where frac_veg_nosno IS NON-ZERO
      ! -----------------------------------------------------------------

      fn = 0
      do fp = 1,num_nolakeurbanp
         p = filter_nolakeurbanp(fp)
         if (frac_veg_nosno(p) /= 0) then
            fn = fn + 1
            filterp(fn) = p
         end if
      end do


      if (use_fates) then
         call alm_fates%prep_canopyfluxes( bounds )
      end if


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
         g=veg_pp%gridcell(p)
         ! calculate dayl_factor as the ratio of (current:max dayl)^2
         ! set a minimum of 0.01 (1%) for the dayl_factor
         dayl_factor(p)=min(1._r8,max(0.01_r8,(dayl(g)*dayl(g))/(max_dayl(g)*max_dayl(g))))
      end do

      rb1(begp:endp) = 0._r8

      !assign the temporary filter
      do f = 1, fn
         p = filterp(f)
         filterc_tmp(f)=veg_pp%column(p)
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
   
      ! --------------------------------------------------------------------------
      ! if this is a FATES simulation
      ! ask fates to calculate btran functions and distribution of uptake
      ! this will require boundary conditions from CLM, boundary conditions which
      ! may only be available from a smaller subset of patches that meet the
      ! exposed veg.  
      ! calc_root_moist_stress already calculated root soil water stress 'rresis'
      ! this is the input boundary condition to calculate the transpiration
      ! wetness factor btran and the root weighting factors for FATES.  These
      ! values require knowledge of the belowground root structure.
      ! --------------------------------------------------------------------------
      
      if(use_fates)then
         call alm_fates%wrap_btran(bounds, fn, filterc_tmp(1:fn), soilstate_vars, &
               temperature_vars, energyflux_vars, soil_water_retention_curve)
         
      else
         !calculate root moisture stress
         call calc_root_moist_stress(bounds,     &
              nlevgrnd = nlevgrnd,               &
              fn = fn,                           &
              filterp = filterp,                 &
              canopystate_vars=canopystate_vars, &
              energyflux_vars=energyflux_vars,   &
              soilstate_vars=soilstate_vars,     &
              temperature_vars=temperature_vars, &
              waterstate_vars=waterstate_vars,   &
              soil_water_retention_curve=soil_water_retention_curve)
         
      end if !use_fates

      ! Determine if irrigation is needed (over irrigated soil columns)

      ! First, determine in what grid cells we need to bother 'measuring' soil water, to see if we need irrigation
      ! Also set n_irrig_steps_left for these grid cells
      ! n_irrig_steps_left(p) > 0 is ok even if irrig_rate(p) ends up = 0
      ! in this case, we'll irrigate by 0 for the given number of time steps

      call get_prev_date(yr, mon, day, time)  ! get time as of beginning of time step

      do f = 1, fn
         p = filterp(f)
         c = veg_pp%column(p)
         g = veg_pp%gridcell(p)
         if ( .not.veg_pp%is_fates(p)             .and. &
              irrigated(veg_pp%itype(p)) == 1._r8 .and. &
              elai(p) > irrig_min_lai          .and. &
              btran(p) < irrig_btran_thresh ) then
            
            ! see if it's the right time of day to start irrigating:
            local_time = modulo(time + nint(grc_pp%londeg(g)/degpsec), isecspday)
            ! local_time = time ! if not using local time
            seconds_since_irrig_start_time = modulo(local_time - irrig_start_time, isecspday)
            if (seconds_since_irrig_start_time < dtime) then
               ! it's time to start irrigating
               check_for_irrig(p)    = .true.
               n_irrig_steps_left(p) = irrig_nsteps_per_day
               irrig_rate(p)         = 0._r8  ! reset; we'll add to this later
            else
               check_for_irrig(p)    = .false.
            end if
         else  ! non-irrig pft or elai<=irrig_min_lai or btran>irrig_btran_thresh
            check_for_irrig(p)       = .false.
         end if
         
      end do


      ! Now 'measure' soil water for the grid cells identified above and see if the 
      ! soil is dry enough to warrant irrigation
      ! (Note: frozen_soil could probably be a column-level variable, but that would be
      ! slightly less robust to potential future modifications)
      ! This should not be operating on FATES patches (see is_fates filter above, pushes
      ! check_for_irrig = false
      frozen_soil(bounds%begp : bounds%endp) = .false.
      do j = 1,nlevgrnd
         do f = 1, fn
            p = filterp(f)
            c = veg_pp%column(p)
            g = veg_pp%gridcell(p)
            if (check_for_irrig(p) .and. .not. frozen_soil(p)) then
               ! if level L was frozen, then we don't look at any levels below L
               if (t_soisno(c,j) <= SHR_CONST_TKFRZ) then
                  frozen_soil(p) = .true.
               else if (rootfr(p,j) > 0._r8) then
                  ! determine soil water deficit in this layer:

                  ! Calculate vol_liq_so - i.e., vol_liq at which smp_node = smpso - by inverting the above equations 
                  ! for the root resistance factors
                  vol_liq_so   = eff_porosity(c,j) * (-smpso(veg_pp%itype(p))/sucsat(c,j))**(-1/bsw(c,j))

                  ! Translate vol_liq_so and eff_porosity into h2osoi_liq_so and h2osoi_liq_sat and calculate deficit
                  h2osoi_liq_so  = vol_liq_so * denh2o * col_pp%dz(c,j)
                  h2osoi_liq_sat = eff_porosity(c,j) * denh2o * col_pp%dz(c,j)
                  deficit        = max((h2osoi_liq_so + ldomain%firrig(g)*(h2osoi_liq_sat - h2osoi_liq_so)) - h2osoi_liq(c,j), 0._r8)

                  ! Add deficit to irrig_rate, converting units from mm to mm/sec
                  irrig_rate(p)  = irrig_rate(p) + deficit/(dtime*irrig_nsteps_per_day)

               end if  ! else if (rootfr(p,j) > 0)
            end if     ! if (check_for_irrig(p) .and. .not. frozen_soil(p))
         end do        ! do f
      end do           ! do j

      ! Modify aerodynamic parameters for sparse/dense canopy (X. Zeng)
      do f = 1, fn
         p = filterp(f)
         c = veg_pp%column(p)

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
         c = veg_pp%column(p)
         t = veg_pp%topounit(p)
         g = veg_pp%gridcell(p)

         ! Net absorbed longwave radiation by canopy and ground
         ! =air+bir*t_veg**4+cir*t_grnd(c)**4

         air(p) =   emv(p) * (1._r8+(1._r8-emv(p))*(1._r8-emg(c))) * forc_lwrad(t)
         bir(p) = - (2._r8-emv(p)*(1._r8-emg(c))) * emv(p) * sb
         cir(p) =   emv(p)*emg(c)*sb

         ! Saturated vapor pressure, specific humidity, and their derivatives
         ! at the leaf surface

         call QSat (t_veg(p), forc_pbot(t), el(p), deldT, qsatl(p), qsatldT(p))

         ! Determine atmospheric co2 and o2

         co2(p) = forc_pco2(t)
         o2(p)  = forc_po2(t)

         if ( use_c13 ) then
            c13o2(p) = forc_pc13o2(t)
         end if

         ! Initialize flux profile

         nmozsgn(p) = 0

         taf(p) = (t_grnd(c) + thm(p))/2._r8
         qaf(p) = (forc_q(t)+qg(c))/2._r8

         ur(p) = max(1.0_r8,sqrt(forc_u(t)*forc_u(t)+forc_v(t)*forc_v(t)))
         dth(p) = thm(p)-taf(p)
         dqh(p) = forc_q(t)-qaf(p)
         delq(p) = qg(c) - qaf(p)
         dthv(p) = dth(p)*(1._r8+0.61_r8*forc_q(t))+0.61_r8*forc_th(t)*dqh(p)
         zldis(p) = forc_hgt_u_patch(p) - displa(p)

         ! Check to see if the forcing height is below the canopy height
         if (zldis(p) < 0._r8) then
            found = .true.
            index = p
         end if

      end do

      if (found) then
         if ( .not. use_fates ) then
            write(iulog,*)'Error: Forcing height is below canopy height for pft index '
            call endrun(decomp_index=index, clmlevel=namep, msg=errmsg(__FILE__, __LINE__))
         end if
      end if

      do f = 1, fn
         p = filterp(f)
         c = veg_pp%column(p)

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
              frictionvel_vars)

         do f = 1, fn
            p = filterp(f)
            c = veg_pp%column(p)
            t = veg_pp%topounit(p)
            g = veg_pp%gridcell(p)

            tlbef(p) = t_veg(p)
            del2(p) = del(p)

            ! Determine aerodynamic resistances

            ram1(p)  = 1._r8/(ustar(p)*ustar(p)/um(p))
            rah(p,1) = 1._r8/(temp1(p)*ustar(p))
            raw(p,1) = 1._r8/(temp2(p)*ustar(p))

            ! Bulk boundary layer resistance of leaves

            uaf(p) = um(p)*sqrt( 1._r8/(ram1(p)*um(p)) )

            ! Use pft parameter for leaf characteristic width
            ! dleaf_patch if this is not an ed patch.
            ! Otherwise, the value has already been loaded
            ! during the FATES dynamics and/or initialization call
            if(.not.veg_pp%is_fates(p)) then  
               dleaf_patch(p) = dleaf(veg_pp%itype(p))
            end if
            
            
            cf  = 0.01_r8/(sqrt(uaf(p))*sqrt( dleaf_patch(p) ))
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
            eah(p) = forc_pbot(t) * qaf(p) / 0.622_r8   ! pa
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
            c = veg_pp%column(p)
            if(.not.veg_pp%is_fates(p)) then
               if (veg_pp%itype(p) == nsoybean .or. veg_pp%itype(p) == nsoybeanirrig) then
                  btran(p) = min(1._r8, btran(p) * 1.25_r8)
               end if
            end if
         end do

         if ( use_fates ) then      

            call alm_fates%wrap_photosynthesis(bounds, fn, filterp(1:fn), &
                  svpts(begp:endp), eah(begp:endp), o2(begp:endp), &
                  co2(begp:endp), rb(begp:endp), dayl_factor(begp:endp), &
                  atm2lnd_vars, temperature_vars, canopystate_vars, photosyns_vars)

         else ! not use_fates

            if ( use_hydrstress ) then
               call PhotosynthesisHydraulicStress (bounds, fn, filterp, &
                    svpts(begp:endp), eah(begp:endp), o2(begp:endp), co2(begp:endp), rb(begp:endp), bsun(begp:endp), &
                    bsha(begp:endp), btran(begp:endp), dayl_factor(begp:endp), &
                    qsatl(begp:endp), qaf(begp:endp),     &
                    atm2lnd_vars, temperature_vars, soilstate_vars, waterstate_vars, surfalb_vars, solarabs_vars,    &
                    canopystate_vars, photosyns_vars, waterflux_vars, &
                    nitrogenstate_vars, phosphorusstate_vars)
            else
               call Photosynthesis (bounds, fn, filterp, &
                 svpts(begp:endp), eah(begp:endp), o2(begp:endp), co2(begp:endp), rb(begp:endp), btran(begp:endp), &
                 dayl_factor(begp:endp), atm2lnd_vars, temperature_vars, surfalb_vars, solarabs_vars, &
                 canopystate_vars, photosyns_vars, nitrogenstate_vars, phosphorusstate_vars, phase='sun')
            end if

            if ( use_c13 ) then
               call Fractionation (bounds, fn, filterp, &
                    atm2lnd_vars, canopystate_vars, cnstate_vars, solarabs_vars, surfalb_vars, photosyns_vars, &
                    phase='sun')
            endif

            do f = 1, fn
               p = filterp(f)
               c = veg_pp%column(p)
               if (veg_pp%itype(p) == nsoybean .or. veg_pp%itype(p) == nsoybeanirrig) then
                  btran(p) = min(1._r8, btran(p) * 1.25_r8)
               end if
            end do
            if ( .not. use_hydrstress ) then
              call Photosynthesis (bounds, fn, filterp, &
                 svpts(begp:endp), eah(begp:endp), o2(begp:endp), co2(begp:endp), rb(begp:endp), btran(begp:endp), &
                 dayl_factor(begp:endp), atm2lnd_vars, temperature_vars, surfalb_vars, solarabs_vars, &
                 canopystate_vars, photosyns_vars, nitrogenstate_vars, phosphorusstate_vars, phase='sha')
            end if

            if ( use_c13 ) then
               call Fractionation (bounds, fn, filterp,  &
                    atm2lnd_vars, canopystate_vars, cnstate_vars, solarabs_vars, surfalb_vars, photosyns_vars, &
                    phase='sha')
            end if

         end if ! end of if use_fates

         do f = 1, fn
            p = filterp(f)
            c = veg_pp%column(p)
            t = veg_pp%topounit(p)
            g = veg_pp%gridcell(p)

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

            if (fdry(p) > 0._r8) then
               rppdry  = fdry(p)*rb(p)*(laisun(p)/(rb(p)+rssun(p)) + &
                    laisha(p)/(rb(p)+rssha(p)))/elai(p)
            else
               rppdry = 0._r8
            end if
               
            ! Calculate canopy conductance for methane / oxygen (e.g. stomatal conductance & leaf bdy cond)
            if (use_lch4) then
               canopy_cond(p) = (laisun(p)/(rb(p)+rssun(p)) + laisha(p)/(rb(p)+rssha(p)))/max(elai(p), 0.01_r8)
            end if

            efpot = forc_rho(t)*wtl*(qsatl(p)-qaf(p))
            ! When the hydraulic stress parameterization is active calculate rpp
            ! but not transpiration
            if ( use_hydrstress ) then
              if (efpot > 0._r8) then
                 if (btran(p) > btran0) then
                   rpp = rppdry + fwet(p)
                 else
                   rpp = fwet(p)
                 end if
                 !Check total evapotranspiration from leaves
                 rpp = min(rpp, (qflx_tran_veg(p)+h2ocan(p)/dtime)/efpot)
              else
                 rpp = 1._r8
              end if
            else

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

            dc1 = forc_rho(t)*cpair*wtl
            dc2 = hvap*forc_rho(t)*wtlq

            efsh   = dc1*(wtga*t_veg(p)-wtg0*t_grnd(c)-wta0(p)*thm(p))
            efe(p) = dc2*(wtgaq*qsatl(p)-wtgq0*qg(c)-wtaq0(p)*forc_q(t))

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

            efpot = forc_rho(t)*wtl*(wtgaq*(qsatl(p)+qsatldT(p)*dt_veg(p)) &
                 -wtgq0*qg(c)-wtaq0(p)*forc_q(t))
            qflx_evap_veg(p) = rpp*efpot

            ! Calculation of evaporative potentials (efpot) and
            ! interception losses; flux in kg m**-2 s-1.  ecidif
            ! holds the excess energy if all intercepted water is evaporated
            ! during the timestep.  This energy is later added to the
            ! sensible heat flux.
            if ( use_hydrstress ) then
               ecidif = max(0._r8,qflx_evap_veg(p)-qflx_tran_veg(p)-h2ocan(p)/dtime)
               qflx_evap_veg(p) = min(qflx_evap_veg(p),qflx_tran_veg(p)+h2ocan(p)/dtime)
            else

              ecidif = 0._r8
              if (efpot > 0._r8 .and. btran(p) > btran0) then
               qflx_tran_veg(p) = efpot*rppdry
              else
               qflx_tran_veg(p) = 0._r8
              end if
              ecidif = max(0._r8, qflx_evap_veg(p)-qflx_tran_veg(p)-h2ocan(p)/dtime)
              qflx_evap_veg(p) = min(qflx_evap_veg(p),qflx_tran_veg(p)+h2ocan(p)/dtime)
            end if

            ! The energy loss due to above two limits is added to
            ! the sensible heat flux.

            eflx_sh_veg(p) = efsh + dc1*wtga*dt_veg(p) + err(p) + erre + hvap*ecidif

            ! Re-calculate saturated vapor pressure, specific humidity, and their
            ! derivatives at the leaf surface

            call QSat(t_veg(p), forc_pbot(t), el(p), deldT, qsatl(p), qsatldT(p))

            ! Update vegetation/ground surface temperature, canopy air
            ! temperature, canopy vapor pressure, aerodynamic temperature, and
            ! Monin-Obukhov stability parameter for next iteration.

            taf(p) = wtg0*t_grnd(c) + wta0(p)*thm(p) + wtl0(p)*t_veg(p)
            qaf(p) = wtlq0(p)*qsatl(p) + wtgq0*qg(c) + forc_q(t)*wtaq0(p)

            ! Update Monin-Obukhov length and wind speed including the
            ! stability effect

            dth(p) = thm(p)-taf(p)
            dqh(p) = forc_q(t)-qaf(p)
            delq(p) = wtalq(p)*qg(c)-wtlq0(p)*qsatl(p)-wtaq0(p)*forc_q(t)

            tstar = temp1(p)*dth(p)
            qstar = temp2(p)*dqh(p)

            thvstar = tstar*(1._r8+0.61_r8*forc_q(t)) + 0.61_r8*forc_th(t)*qstar
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

         end do   ! end of filtered pft loop

         do f = 1, fn
           p = filterp(f)
           t = veg_pp%topounit(p)
           lbl_rsc_h2o(p) = getlblcef(forc_rho(t),t_veg(p))*uaf(p)/(uaf(p)**2._r8+1.e-10_r8)   !laminar boundary resistance for h2o over leaf, should I make this consistent for latent heat calculation?
         enddo
            
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
         c = veg_pp%column(p)
         t = veg_pp%topounit(p)
         g = veg_pp%gridcell(p)

         ! Energy balance check in canopy

         lw_grnd=(frac_sno(c)*t_soisno(c,snl(c)+1)**4 &
              +(1._r8-frac_sno(c)-frac_h2osfc(c))*t_soisno(c,1)**4 &
              +frac_h2osfc(c)*t_h2osfc(c)**4)

         err(p) = sabv(p) + air(p) + bir(p)*tlbef(p)**3*(tlbef(p) + 4._r8*dt_veg(p)) &
                                !+ cir(p)*t_grnd(c)**4 - eflx_sh_veg(p) - hvap*qflx_evap_veg(p)
              + cir(p)*lw_grnd - eflx_sh_veg(p) - hvap*qflx_evap_veg(p)

         ! Fluxes from ground to canopy space

         delt    = wtal(p)*t_grnd(c)-wtl0(p)*t_veg(p)-wta0(p)*thm(p)
         taux(p) = -forc_rho(t)*forc_u(t)/ram1(p)
         tauy(p) = -forc_rho(t)*forc_v(t)/ram1(p)
         eflx_sh_grnd(p) = cpair*forc_rho(t)*wtg(p)*delt

         ! compute individual sensible heat fluxes
         delt_snow = wtal(p)*t_soisno(c,snl(c)+1)-wtl0(p)*t_veg(p)-wta0(p)*thm(p)
         eflx_sh_snow(p) = cpair*forc_rho(t)*wtg(p)*delt_snow

         delt_soil  = wtal(p)*t_soisno(c,1)-wtl0(p)*t_veg(p)-wta0(p)*thm(p)
         eflx_sh_soil(p) = cpair*forc_rho(t)*wtg(p)*delt_soil

         delt_h2osfc  = wtal(p)*t_h2osfc(c)-wtl0(p)*t_veg(p)-wta0(p)*thm(p)
         eflx_sh_h2osfc(p) = cpair*forc_rho(t)*wtg(p)*delt_h2osfc
         qflx_evap_soi(p) = forc_rho(t)*wtgq(p)*delq(p)

         ! compute individual latent heat fluxes
         delq_snow = wtalq(p)*qg_snow(c)-wtlq0(p)*qsatl(p)-wtaq0(p)*forc_q(t)
         qflx_ev_snow(p) = forc_rho(t)*wtgq(p)*delq_snow

         delq_soil = wtalq(p)*qg_soil(c)-wtlq0(p)*qsatl(p)-wtaq0(p)*forc_q(t)
         qflx_ev_soil(p) = forc_rho(t)*wtgq(p)*delq_soil

         delq_h2osfc = wtalq(p)*qg_h2osfc(c)-wtlq0(p)*qsatl(p)-wtaq0(p)*forc_q(t)
         qflx_ev_h2osfc(p) = forc_rho(t)*wtgq(p)*delq_h2osfc

         ! 2 m height air temperature

         t_ref2m(p) = thm(p) + temp1(p)*dth(p)*(1._r8/temp12m(p) - 1._r8/temp1(p))
         t_ref2m_r(p) = t_ref2m(p)

         ! 2 m height specific humidity

         q_ref2m(p) = forc_q(t) + temp2(p)*dqh(p)*(1._r8/temp22m(p) - 1._r8/temp2(p))

         ! 2 m height relative humidity

         call QSat(t_ref2m(p), forc_pbot(t), e_ref2m, de2mdT, qsat_ref2m, dqsat2mdT)
         rh_ref2m(p) = min(100._r8, q_ref2m(p) / qsat_ref2m * 100._r8)
         rh_ref2m_r(p) = rh_ref2m(p)

         ! Downward longwave radiation below the canopy

         dlrad(p) = (1._r8-emv(p))*emg(c)*forc_lwrad(t) + &
              emv(p)*emg(c)*sb*tlbef(p)**3*(tlbef(p) + 4._r8*dt_veg(p))

         ! Upward longwave radiation above the canopy

         ulrad(p) = ((1._r8-emg(c))*(1._r8-emv(p))*(1._r8-emv(p))*forc_lwrad(t) &
              + emv(p)*(1._r8+(1._r8-emg(c))*(1._r8-emv(p)))*sb*tlbef(p)**3*(tlbef(p) + &
              4._r8*dt_veg(p)) + emg(c)*(1._r8-emv(p))*sb*lw_grnd)

         ! Derivative of soil energy flux with respect to soil temperature

         cgrnds(p) = cgrnds(p) + cpair*forc_rho(t)*wtg(p)*wtal(p)
         cgrndl(p) = cgrndl(p) + forc_rho(t)*wtgq(p)*wtalq(p)*dqgdT(c)
         cgrnd(p)  = cgrnds(p) + cgrndl(p)*htvp(c)

         ! Update dew accumulation (kg/m2)

         h2ocan(p) = max(0._r8,h2ocan(p)+(qflx_tran_veg(p)-qflx_evap_veg(p))*dtime)

      end do

      if ( use_fates ) then
         call alm_fates%wrap_accumulatefluxes(bounds,fn,filterp(1:fn))
         call alm_fates%wrap_hydraulics_drive(bounds,fn,filterp(1:fn),soilstate_vars, &
              solarabs_vars,energyflux_vars)

      else

         ! Determine total photosynthesis
         
         call PhotosynthesisTotal(fn, filterp, &
              atm2lnd_vars, cnstate_vars, canopystate_vars, photosyns_vars)
         
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
         
      end if

    end associate


  end subroutine CanopyFluxes

end module CanopyFluxesMod

