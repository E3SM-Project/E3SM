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
  use elm_varctl            , only : use_hydrstress, use_finetop_rad
  use elm_varpar            , only : nlevgrnd, nlevsno
  use elm_varcon            , only : namep
  use elm_varcon            , only : mm_epsilon
  use elm_varcon            , only : pa_to_kpa
  use pftvarcon             , only : crop, nfixer
  use decompMod             , only : bounds_type
  use PhotosynthesisMod     , only : Photosynthesis, PhotosynthesisTotal, Fractionation, PhotoSynthesisHydraulicStress
  use SoilMoistStressMod    , only : calc_effective_soilporosity, calc_volumetric_h2oliq
  use SoilMoistStressMod    , only : calc_root_moist_stress, set_perchroot_opt
  use SurfaceResistanceMod  , only : do_soilevap_beta
  use VegetationPropertiesType        , only : veg_vp
  use atm2lndType           , only : atm2lnd_type
  use CanopyStateType       , only : canopystate_type, perchroot, perchroot_alt
  use CNStateType           , only : cnstate_type
  use EnergyFluxType        , only : energyflux_type
  use FrictionvelocityType  , only : frictionvel_type
  use SoilStateType         , only : soilstate_type
  use SolarAbsorbedType     , only : solarabs_type
  use SurfaceAlbedoType     , only : surfalb_type
  use CH4Mod                , only : ch4_type
  use PhotosynthesisType    , only : photosyns_type
  use GridcellType          , only : grc_pp
  use TopounitDataType      , only : top_as, top_af
  use ColumnType            , only : col_pp
  use ColumnDataType        , only : col_es, col_ef, col_ws
  use VegetationType        , only : veg_pp
  use VegetationDataType    , only : veg_es, veg_ef, veg_ws, veg_wf
  ! using elm_instMod messes with the compilation order
  use elm_instMod           , only : alm_fates, soil_water_retention_curve, atm2lnd_vars
  use perf_mod, only: t_startf, t_stopf
  use timeinfoMod
  use spmdmod          , only: masterproc
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CanopyFluxes

contains

  !------------------------------------------------------------------------------
  subroutine CanopyFluxes(bounds,  num_nolu_barep, filter_nolu_barep, &
       num_nolu_vegp, filter_nolu_vegp , &
       canopystate_vars, cnstate_vars, energyflux_vars, &
       frictionvel_vars, soilstate_vars, solarabs_vars, surfalb_vars, &
       ch4_vars, photosyns_vars)

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
    use shr_flux_mod       , only : shr_flux_update_stress
    use elm_varcon         , only : sb, cpair, hvap, vkc, grav, denice
    use elm_varcon         , only : denh2o, tfrz, csoilc, tlsai_crit, alpha_aero
    use elm_varcon         , only : isecspday, degpsec
    use pftvarcon          , only : irrigated
    use elm_varcon         , only : c14ratio
    use shr_const_mod      , only : SHR_CONST_PI
    use elm_varsur         , only : firrig
    use TopounitType       , only : top_pp
    use QSatMod            , only : QSat
    use FrictionVelocityMod, only : FrictionVelocity, MoninObukIni, &
         implicit_stress, atm_gustiness, force_land_gustiness
    use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
    use SurfaceResistanceMod, only : getlblcef
    use PhotosynthesisType, only : photosyns_vars_TimeStepInit
    !
    ! !ARGUMENTS:
    type(bounds_type)         , intent(in)    :: bounds
    integer                   , intent(in)    :: num_nolu_barep
    integer                   , intent(in)    :: filter_nolu_barep(:)
    integer                   , intent(in)    :: num_nolu_vegp
    integer                   , intent(in)    :: filter_nolu_vegp(:)
    type(canopystate_type)    , intent(inout) :: canopystate_vars
    type(cnstate_type)        , intent(inout) :: cnstate_vars  
    type(energyflux_type)     , intent(inout) :: energyflux_vars
    type(frictionvel_type)    , intent(inout) :: frictionvel_vars
    type(solarabs_type)       , intent(inout) :: solarabs_vars
    type(surfalb_type)        , intent(inout) :: surfalb_vars
    type(soilstate_type)      , intent(inout) :: soilstate_vars
    type(ch4_type)            , intent(inout) :: ch4_vars
    type(photosyns_type)      , intent(inout) :: photosyns_vars
    integer  ::  time
    !
    ! !LOCAL VARIABLES:
    real(r8), pointer   :: bsun(:)          ! sunlit canopy transpiration wetness factor (0 to 1)
    real(r8), pointer   :: bsha(:)          ! shaded canopy transpiration wetness factor (0 to 1)
    real(r8), parameter :: btran0 = 0.0_r8  ! initial value
    real(r8), parameter :: zii = 1000.0_r8  ! convective boundary layer height [m]
    real(r8), parameter :: beta = 1.0_r8    ! coefficient of convective velocity [-]
    real(r8), parameter :: delmax = 1.0_r8  ! maxchange in  leaf temperature [K]
    real(r8), parameter :: dlemin = 0.1_r8  ! max limit for energy flux convergence [w/m2]
    real(r8), parameter :: dtmin = 0.01_r8  ! max limit for temperature convergence [K]
    real(r8), parameter :: dtaumin = 0.01_r8! max limit for stress convergence [Pa]
    integer , parameter :: itmax = 41       ! maximum number of iteration [-]
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

    real(r8) :: zldis(num_nolu_vegp)   ! reference height "minus" zero displacement height [m]
    real(r8) :: ugust_total(num_nolu_vegp) ! gustiness including convective velocity [m/s]
    integer  :: filterc_tmp(num_nolu_vegp) ! column of each filter_nolu_vegp patch (FATES btran)
    real(r8) :: wc                     ! convective velocity [m/s]
    real(r8) :: dth(num_nolu_vegp)     ! diff of virtual temp. between ref. height and surface
    real(r8) :: dthv(num_nolu_vegp)    ! diff of vir. poten. temp. between ref. height and surface
    real(r8) :: dqh(num_nolu_vegp)     ! diff of humidity between ref. height and surface
    real(r8) :: obu(num_nolu_vegp)     ! Monin-Obukhov length (m)
    real(r8) :: um (num_nolu_vegp)     ! wind speed including the stablity effect [m/s]
    real(r8) :: ur (num_nolu_vegp)     ! wind speed at reference height [m/s]
    real(r8) :: uaf(num_nolu_vegp)     ! velocity of air within foliage [m/s]
    real(r8) :: temp1(num_nolu_vegp)   ! relation for potential temperature profile
    real(r8) :: temp12m(num_nolu_vegp) ! relation for potential temperature profile applied at 2-m
    real(r8) :: temp2  (num_nolu_vegp) ! relation for specific humidity profile
    real(r8) :: temp22m(num_nolu_vegp) ! relation for specific humidity profile applied at 2-m
    real(r8) :: ustar  (num_nolu_vegp) ! friction velocity [m/s]
    real(r8) :: tstar                  ! temperature scaling parameter
    real(r8) :: qstar                  ! moisture scaling parameter
    real(r8) :: thvstar                ! virtual potential temperature scaling parameter
    real(r8) :: taf(num_nolu_vegp)     ! air temperature within canopy space [K]
    real(r8) :: qaf(num_nolu_vegp)     ! humidity of canopy air [kg/kg]
    real(r8) :: rpp                    ! fraction of potential evaporation from leaf [-]
    real(r8) :: rppdry                 ! fraction of potential evaporation through transp [-]
    real(r8) :: cf                     ! heat transfer coefficient from leaves [-]
    real(r8) :: cf_bare                ! heat transfer coefficient from bare ground [-]
    real(r8) :: rb(num_nolu_vegp)      ! leaf boundary layer resistance [s/m]
    real(r8) :: rah(num_nolu_vegp,2)   ! thermal resistance [s/m]
    real(r8) :: raw(num_nolu_vegp,2)   ! moisture resistance [s/m]
    real(r8) :: wta                    ! heat conductance for air [m/s]
    real(r8) :: wtg(num_nolu_vegp)     ! heat conductance for ground [m/s]
    real(r8) :: wtl                    ! heat conductance for leaf [m/s]
    real(r8) :: wta0(num_nolu_vegp)    ! normalized heat conductance for air [-]
    real(r8) :: wtl0(num_nolu_vegp)    ! normalized heat conductance for leaf [-]
    real(r8) :: wtg0                   ! normalized heat conductance for ground [-]
    real(r8) :: wtal(num_nolu_vegp)    ! normalized heat conductance for air and leaf [-]
    real(r8) :: wtga                   ! normalized heat cond. for air and ground  [-]
    real(r8) :: wtaq                   ! latent heat conductance for air [m/s]
    real(r8) :: wtlq                   ! latent heat conductance for leaf [m/s]
    real(r8) :: wtgq(num_nolu_vegp)    ! latent heat conductance for ground [m/s]
    real(r8) :: wtaq0(num_nolu_vegp)   ! normalized latent heat conductance for air [-]
    real(r8) :: wtlq0(num_nolu_vegp)   ! normalized latent heat conductance for leaf [-]
    real(r8) :: wtgq0                  ! normalized heat conductance for ground [-]
    real(r8) :: wtalq(num_nolu_vegp)   ! normalized latent heat cond. for air and leaf [-]
    real(r8) :: wtgaq                  ! normalized latent heat cond. for air and ground [-]
    real(r8) :: el(num_nolu_vegp)      ! vapor pressure on leaf surface [pa]
    real(r8) :: deldT                  ! derivative of "el" on "t_veg" [pa/K]
    real(r8) :: qsatl(num_nolu_vegp)   ! leaf specific humidity [kg/kg]
    real(r8) :: qsatldT(num_nolu_vegp) ! derivative of "qsatl" on "t_veg"
    real(r8) :: e_ref2m                ! 2 m height surface saturated vapor pressure [Pa]
    real(r8) :: de2mdT                 ! derivative of 2 m height surface saturated vapor pressure on t_ref2m
    real(r8) :: qsat_ref2m             ! 2 m height surface saturated specific humidity [kg/kg]
    real(r8) :: dqsat2mdT              ! derivative of 2 m height surface saturated specific humidity on t_ref2m
    real(r8) :: air(num_nolu_vegp)     ! atmos. radiation temporay set
    real(r8) :: bir(num_nolu_vegp)     ! atmos. radiation temporay set
    real(r8) :: cir(num_nolu_vegp)     ! atmos. radiation temporay set
    real(r8) :: dc1,dc2                ! derivative of energy flux [W/m2/K]
    real(r8) :: delt                   ! temporary
    real(r8) :: delq(num_nolu_vegp)    ! temporary
    real(r8) :: del(num_nolu_vegp)     ! absolute change in leaf temp in current iteration [K]
    real(r8) :: del2(num_nolu_vegp)    ! change in leaf temperature in previous iteration [K]
    real(r8) :: dele(num_nolu_vegp)    ! change in latent heat flux from leaf [K]
    real(r8) :: dels                   ! change in leaf temperature in current iteration [K]
    real(r8) :: det(num_nolu_vegp)     ! maximum leaf temp. change in two consecutive iter [K]
    real(r8) :: efeb(num_nolu_vegp)    ! latent heat flux from leaf (previous iter) [mm/s]
    real(r8) :: efeold                 ! latent heat flux from leaf (previous iter) [mm/s]
    real(r8) :: efpot                  ! potential latent energy flux [kg/m2/s]
    real(r8) :: efe(num_nolu_vegp)     ! water flux from leaf [mm/s]
    real(r8) :: efsh                   ! sensible heat from leaf [mm/s]
    real(r8) :: obuold(num_nolu_vegp)  ! monin-obukhov length from previous iteration
    real(r8) :: tlbef(num_nolu_vegp)   ! leaf temperature from previous iteration [K]
    real(r8) :: ecidif                 ! excess energies [W/m2]
    real(r8) :: err(num_nolu_vegp)     ! balance error
    real(r8) :: erre                   ! balance error
    real(r8) :: co2(num_nolu_vegp)     ! atmospheric co2 partial pressure (pa)
    real(r8) :: o2(num_nolu_vegp)      ! atmospheric o2 partial pressure (pa)
    real(r8) :: svpts(num_nolu_vegp)   ! saturation vapor pressure at t_veg (pa)
    real(r8) :: eah(num_nolu_vegp)     ! canopy air vapor pressure (pa)
    real(r8) :: s_node                 ! vol_liq/eff_porosity
    real(r8) :: smp_node               ! matrix potential
    real(r8) :: smp_node_lf            ! F. Li and S. Levis
    real(r8) :: vol_liq                ! partial volume of liquid water in layer
    integer  :: itlef                  ! counter for leaf temperature iteration [-]
    integer  :: nmozsgn(num_nolu_vegp) ! number of times stability changes sign
    real(r8) :: w                      ! exp(-LSAI)
    real(r8) :: csoilcn                ! interpolated csoilc for less than dense canopies
    real(r8) :: fm(num_nolu_vegp)    ! needed for BGC only to diagnose 10m wind speed
    real(r8) :: wtshi                  ! sensible heat resistance for air, grnd and leaf [-]
    real(r8) :: wtsqi                  ! latent heat resistance for air, grnd and leaf [-]
    integer  :: j                      ! soil/snow level index
    integer  :: p                      ! patch index
    integer  :: c                      ! column index
    integer  :: l                      ! landunit index
    integer  :: t                      ! topounit index
    integer  :: g                      ! gridcell index
    integer  :: tpu_ind                ! index of topounit to grid
    integer  :: fp                     ! lake filter pft index
    integer  :: fn                     ! number of values in vegetated pft filter
    integer  :: fnorig                 ! number of values in pft filter copy
    integer  :: fporig(num_nolu_vegp)  ! temporary filter
    integer  :: fnold                       ! temporary copy of pft count
    integer  :: filter_index                ! filter index
    logical  :: found                       ! error flag for canopy above forcing hgt
    integer  :: index                       ! patch index for error
    real(r8) :: egvf                        ! effective green vegetation fraction
    real(r8) :: lt                          ! elai+esai
    real(r8) :: ri                          ! stability parameter for under canopy air (unitless)
    real(r8) :: csoilb                      ! turbulent transfer coefficient over bare soil (unitless)
    real(r8) :: ricsoilc                    ! modified transfer coefficient under dense canopy (unitless)
    real(r8) :: snow_depth_c                ! critical snow depth to cover plant litter (m)
    real(r8) :: rdl                         ! dry litter layer resistance for water vapor  (s/m)
    real(r8) :: elai_dl                     ! exposed (dry) plant litter area index
    real(r8) :: fsno_dl                     ! effective snow cover over plant litter
    real(r8) :: dayl_factor(num_nolu_vegp) ! scalar (0-1) for daylength effect on Vcmax

    ! If no unfrozen layers, put all in the top layer.
    real(r8) :: delt_snow
    real(r8) :: delt_soil
    real(r8) :: delt_h2osfc
    real(r8) :: lw_grnd
    real(r8) :: delq_snow
    real(r8) :: delq_soil
    real(r8) :: delq_h2osfc
    integer  :: local_time                     ! local time at start of time step (seconds after solar midnight)
    integer  :: seconds_since_irrig_start_time
    integer  :: irrig_nsteps_per_day           ! number of time steps per day in which we irrigate
    logical  :: check_for_irrig(num_nolu_vegp) ! where do we need to check soil moisture to see if we need to irrigate?
    logical  :: frozen_soil(num_nolu_vegp)     ! set to true if we have encountered a frozen soil layer
    real(r8) :: vol_liq_so                     ! partial volume of liquid water in layer for which smp_node = smpso
    real(r8) :: h2osoi_liq_so                  ! liquid water corresponding to vol_liq_so for this layer [kg/m2]
    real(r8) :: h2osoi_liq_sat                 ! liquid water corresponding to eff_porosity for this layer [kg/m2]
    real(r8) :: deficit                        ! difference between desired soil moisture level for this layer and
                                               ! current soil moisture level [kg/m2]
    real(r8) :: dt_veg(num_nolu_vegp)          ! change in t_veg, last iteration (Kelvin)
    integer  :: ft                             ! plant functional type index
    real(r8) :: temprootr,sum1
    integer  :: iv
    real(r8) :: wind_speed0(bounds%begp:bounds%endp) ! Wind speed from atmosphere at start of iteration
    real(r8) :: wind_speed_adj(bounds%begp:bounds%endp) ! Adjusted wind speed for iteration
    real(r8) :: tau(bounds%begp:bounds%endp)      ! Stress used in iteration
    real(r8) :: tau_diff(bounds%begp:bounds%endp) ! Difference from previous iteration tau
    real(r8) :: prev_tau(bounds%begp:bounds%endp) ! Previous iteration tau
    real(r8) :: prev_tau_diff(bounds%begp:bounds%endp) ! Previous difference in iteration tau
    real(r8) :: slope_rad, deg2rad
    ! iter_filterp(1:fn) holds the patch index "p" for patches still iterating (unconverged).
    ! iter_filter_map(1:fn) holds, for each entry, the ORIGINAL compressed filter_index
    ! (i.e. its position in filter_nolu_vegp) which is how all the num_nolu_vegp-sized
    ! local arrays below (del, rb, obu, etc.) are permanently indexed. As patches converge,
    ! iter_filterp/iter_filter_map are compacted (fn shrinks) but the local arrays are not
    ! moved, so iter_filter_map must be used to address them.
    integer :: iter_filterp(num_nolu_vegp), iter_filter_map(num_nolu_vegp)  ! filter for iteration loop
    integer :: active_index                                                ! position within the compacted iteration filter

    ! Indices for raw and rah
    integer, parameter :: above_canopy = 1         ! Above canopy
    integer, parameter :: below_canopy = 2         ! Below canopy

    ! Lower bound for VPD (based on CLM)
    real(r8), parameter :: vpd_min = 50._r8
    !------------------------------------------------------------------------------

    associate(                                                               &
         snl                  => col_pp%snl                                   , & ! Input:  [integer  (:)   ]  number of snow layers
         dayl                 => grc_pp%dayl                                  , & ! Input:  [real(r8) (:)   ]  daylength (s)
         max_dayl             => grc_pp%max_dayl                              , & ! Input:  [real(r8) (:)   ]  maximum daylength for this grid cell (s)
         slope_deg            => grc_pp%slope_deg                             , &

         forc_lwrad           => top_af%lwrad                              , & ! Input:  [real(r8) (:)   ]  downward infrared (longwave) radiation (W/m**2)
         forc_q               => top_as%qbot                               , & ! Input:  [real(r8) (:)   ]  atmospheric specific humidity (kg/kg)
         forc_pbot            => top_as%pbot                               , & ! Input:  [real(r8) (:)   ]  atmospheric pressure (Pa)
         forc_th              => top_as%thbot                              , & ! Input:  [real(r8) (:)   ]  atmospheric potential temperature (Kelvin)
         forc_rho             => top_as%rhobot                             , & ! Input:  [real(r8) (:)   ]  air density (kg/m**3)
         forc_t               => top_as%tbot                               , & ! Input:  [real(r8) (:)   ]  atmospheric temperature (Kelvin)
         forc_u               => top_as%ubot                               , & ! Input:  [real(r8) (:)   ]  atmospheric wind speed in east direction (m/s)
         forc_v               => top_as%vbot                               , & ! Input:  [real(r8) (:)   ]  atmospheric wind speed in north direction (m/s)
         wsresp               => top_as%wsresp                             , & ! Input:  [real(r8) (:)   ]  response of wind to surface stress (m/s/Pa)
         tau_est              => top_as%tau_est                            , & ! Input:  [real(r8) (:)   ]  approximate atmosphere change to zonal wind (m/s)
         ugust                => top_as%ugust                              , & ! Input:  [real(r8) (:)   ]  gustiness from atmosphere (m/s)
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
         num_iter             => frictionvel_vars%num_iter_patch           , & ! Output: number of iterations required
         forc_hgt_t_patch => frictionvel_vars%forc_hgt_t_patch , & ! Input:  [real(r8) (:) ] observational height of temperature at pft level [m]
         forc_hgt_q_patch => frictionvel_vars%forc_hgt_q_patch , & ! Input:  [real(r8) (:) ] observational height of specific humidity at pft level [m]
         vds              => frictionvel_vars%vds_patch        , & ! Output: [real(r8) (:) ] dry deposition velocity term (m/s) (for SO4 NH4NO3)
         u10              => frictionvel_vars%u10_patch        , & ! Output: [real(r8) (:) ] 10-m wind (m/s) (for dust model)
         u10_elm          => frictionvel_vars%u10_elm_patch    , & ! Output: [real(r8) (:) ] 10-m wind (m/s)
         va               => frictionvel_vars%va_patch         , & ! Output: [real(r8) (:) ] atmospheric wind speed plus convective velocity (m/s)
         fv               => frictionvel_vars%fv_patch         ,  & ! Output: [real(r8) (:) ] friction velocity (m/s) (for dust model)

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

         !pgwgt                => veg_pp%wtgcell              , & ! Input:  [integer  (:)   ]  pft's weight in gridcell
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
         rah_above            => frictionvel_vars%rah_above_patch , & ! Output: [real(r8) (:)   ]  above-canopy sensible heat flux resistance [s/m]
         rah_below            => frictionvel_vars%rah_above_patch , & ! Output: [real(r8) (:)   ]  below-canopy sensible heat flux resistance [s/m]
         raw_above            => frictionvel_vars%raw_below_patch , & ! Output: [real(r8) (:)   ]  above-canopy water vapour flux resistance [s/m]
         raw_below            => frictionvel_vars%raw_below_patch , & ! Output: [real(r8) (:)   ]  below-canopy water vapour flux resistance [s/m]
         ustar                => frictionvel_vars%ustar_patch     , & ! Output: [real(r8) (:)   ]  friction velocity [m/s]
         um                   => frictionvel_vars%um_patch        , & ! Output: [real(r8) (:)   ]  wind speed including the stablity effect [m/s]
         uaf                  => frictionvel_vars%uaf_patch       , & ! Output: [real(r8) (:)   ]  canopy air wind speed [m/s]
         taf                  => frictionvel_vars%taf_patch       , & ! Output: [real(r8) (:)   ]  canopy air temperature [K]
         qaf                  => frictionvel_vars%qaf_patch       , & ! Output: [real(r8) (:)   ]  canopy air specific humidity [kg/kg]
         obu                  => frictionvel_vars%obu_patch       , & ! Output: [real(r8) (:)   ]  Obukhov length scale [m]
         zeta                 => frictionvel_vars%zeta_patch      , & ! Output: [real(r8) (:)   ]  dimensionless stability parameter 
         vpd                  => frictionvel_vars%vpd_patch       , & ! Output: [real(r8) (:)   ]  vapour pressure deficit [kPa]
         begp                 => bounds%begp                               , &
         endp                 => bounds%endp                                 &
         )
      if (use_hydrstress) then
        bsun                    => energyflux_vars%bsun_patch ! Output:[real(r8) (:)   ]  sunlit canopy transpiration wetness factor (0 to 1)
        bsha                    => energyflux_vars%bsha_patch ! Output:[real(r8) (:)   ]  sunlit canopy transpiration wetness factor (0 to 1)
      end if

     
      fn = num_nolu_vegp
      ! First - set the following values over points where frac vegetation covered by snow is zero
      ! (e.g. btran, t_veg, rootr, rresis)
      if(num_nolu_barep > 0) then
         !$acc parallel loop independent gang vector private(p,c,t) default(present) &
         !$acc   present(t_veg(:), btran(:), rssun(:), rssha(:), lbl_rsc_h2o(:), thm(:) ) 
         do fp = 1,num_nolu_barep
            p = filter_nolu_barep(fp)
            c = veg_pp%column(p)
            t = veg_pp%topounit(p)
            btran(p) = 0._r8
            t_veg(p) = forc_t(t)
            cf_bare  = forc_pbot(t)/(SHR_CONST_RGAS*0.001_r8*thm(p))*1.e06_r8
            rssun(p) = 1._r8/1.e15_r8 * cf_bare
            rssha(p) = 1._r8/1.e15_r8 * cf_bare
            lbl_rsc_h2o(p)=0._r8
         end do
         !$acc parallel loop independent gang default(present)
         do j = 1, nlevgrnd
            !$acc loop vector private(p)
            do fp = 1,num_nolu_barep
               p = filter_nolu_barep(fp)
               rootr(p,j)  = 0._r8
               rresis(p,j) = 0._r8
            end do
         end do
         
            
      end if

      deg2rad = SHR_CONST_PI/180._r8
      if(num_nolu_vegp == 0) return
      time = secs_curr;
      irrig_nsteps_per_day = ((irrig_length + (dtime_mod - 1))/dtime_mod)  ! round up

      !$acc enter data copyin(time,irrig_nsteps_per_day) 

      !$acc enter data create(del(:), efeb(:), wtlq0(:),wtalq(:), &
      !$acc    wtgq(:), wtaq0(:), obuold(:),dayl_factor(:),check_for_irrig(:),zldis(:) ) 
      !$acc enter data create(iter_filterp(:), iter_filter_map(:))

      ! Initialize
      !$acc parallel loop independent gang vector private(p) default(present) present(btran(:), btran2(:))
      do filter_index = 1, fn
         iter_filterp(filter_index) = filter_nolu_vegp(filter_index) 
         iter_filter_map(filter_index) = filter_index
         p = iter_filterp(filter_index)
         del(filter_index)    = 0._r8  ! change in leaf temperature from previous iteration
         efeb(filter_index)   = 0._r8  ! latent head flux from leaf for previous iteration
         wtlq0(filter_index)  = 0._r8
         wtalq(filter_index)  = 0._r8
         wtgq(filter_index)   = 0._r8
         wtaq0(filter_index)  = 0._r8
         obuold(filter_index) = 0._r8
         btran(p)  = btran0
         btran2(p)  = btran0
      end do

      ! calculate daylength control for Vcmax
      !$acc parallel loop independent gang vector private(p,g) default(present)
      do filter_index = 1, fn
         p=filter_nolu_vegp(filter_index)
         g=veg_pp%gridcell(p)
         ! calculate dayl_factor as the ratio of (current:max dayl)^2
         ! set a minimum of 0.01 (1%) for the dayl_factor
         dayl_factor(filter_index)=min(1._r8,max(0.01_r8,(dayl(g)*dayl(g))/(max_dayl(g)*max_dayl(g))))
      end do
      ! -----------------------------------------------------------------
      ! Time step initialization of photosynthesis variables
      ! -----------------------------------------------------------------
      !NOTE: This likely shouldn't init based on bounds but on a filter !!
      call photosyns_vars_TimeStepInit(photosyns_vars,bounds)

      if (use_fates) then
         call alm_fates%prep_canopyfluxes( bounds )
      end if
      
      rb1(begp:endp) = 0._r8


      !NOTE: this filter set up means doing the same calculations for a column
      !      redundantly (eg. computing the same column variable 4 times)
      !      It would be best to make nolu_vegc filter.
      !assign the temporary filter
      
      ! compute effective soil porosity
      call calc_effective_soilporosity(bounds,                          &
           ubj = nlevgrnd,                                              &
           numf = fn,                                                   &
           filter = filter_nolu_vegp(1:fn),                             &
           watsat = watsat(bounds%begc:bounds%endc, 1:nlevgrnd),        &
           h2osoi_ice = h2osoi_ice(bounds%begc:bounds%endc,1:nlevgrnd), &
           denice = denice,                                             &
           eff_por=eff_porosity(bounds%begc:bounds%endc, 1:nlevgrnd) )

      !compute volumetric liquid water content

      call calc_volumetric_h2oliq(bounds,                                    &
           lbj = 1,                                                          &
           ubj = nlevgrnd,                                                   &
           numf = fn,                                                        &
           filter = filter_nolu_vegp(1:fn),                                       &
           eff_porosity = eff_porosity(bounds%begc:bounds%endc, 1:nlevgrnd), &
           h2osoi_liq = h2osoi_liq(bounds%begc:bounds%endc, 1:nlevgrnd),     &
           denh2o = denh2o,                                                  &
           vol_liq = h2osoi_liqvol(bounds%begc:bounds%endc, 1:nlevgrnd) )

      ! set up perchroot options
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
         do fp = 1, fn
            filterc_tmp(fp) = veg_pp%column(filter_nolu_vegp(fp))
         end do
         call alm_fates%wrap_btran(bounds, fn, filterc_tmp(1:fn), soilstate_vars, &
               energyflux_vars, soil_water_retention_curve)
      else
         !calculate root moisture stress
         call calc_root_moist_stress(bounds,     &
              nlevgrnd = nlevgrnd,               &
              fn = fn,                           &
              filterp = filter_nolu_vegp,                 &
              canopystate_vars=canopystate_vars, &
              energyflux_vars=energyflux_vars,   &
              soilstate_vars=soilstate_vars      &
              )
      end if !use_fates

      ! Determine if irrigation is needed (over irrigated soil columns)
      ! First, determine in what grid cells we need to bother 'measuring' soil water, to see if we need irrigation
      ! Also set n_irrig_steps_left for these grid cells
      ! n_irrig_steps_left(p) > 0 is ok even if irrig_rate(p) ends up = 0
      ! in this case, we'll irrigate by 0 for the given number of time steps

      !$acc parallel loop independent gang vector default(present) present(btran(:),elai(:),n_irrig_steps_left(:), irrig_rate(:)) private(p,g,local_time,seconds_since_irrig_start_time)
      do filter_index = 1, fn
         p = filter_nolu_vegp(filter_index)
         g = veg_pp%gridcell(p)
         if ( .not.veg_pp%is_fates(p)             .and. &
              irrigated(veg_pp%itype(p)) == 1._r8 .and. &
              elai(p) > irrig_min_lai  .and. btran(p) < irrig_btran_thresh ) then

            ! see if it's the right time of day to start irrigating:
            local_time = modulo(secs_curr + nint(grc_pp%londeg(g)/degpsec), isecspday)
            seconds_since_irrig_start_time = modulo(local_time - irrig_start_time, isecspday)
            if (seconds_since_irrig_start_time < dtime_mod) then
               ! it's time to start irrigating
               check_for_irrig(filter_index)    = .true.
               n_irrig_steps_left(p) = irrig_nsteps_per_day
               irrig_rate(p)         = 0._r8  ! reset; we'll add to this later
            else
               check_for_irrig(filter_index)    = .false.
            end if
         else  ! non-irrig pft or elai<=irrig_min_lai or btran>irrig_btran_thresh
            check_for_irrig(filter_index)       = .false.
         end if

      end do

      ! Now 'measure' soil water for the grid cells identified above and see if the
      ! soil is dry enough to warrant irrigation
      ! (Note: frozen_soil could probably be a column-level variable, but that would be
      ! slightly less robust to potential future modifications)
      ! This should not be operating on FATES patches (see is_fates filter above, pushes
      ! check_for_irrig = false
      ! frozen_soil(1:fn) = .false.
      !$acc parallel loop independent gang worker default(present) private(p,c,g)
      do filter_index = 1, fn
         p = filter_nolu_vegp(filter_index)
         c = veg_pp%column(p)
         g = veg_pp%gridcell(p)
         tpu_ind = top_pp%topo_grc_ind(t)  !Get topounit index on the grid
         if (check_for_irrig(filter_index)) then
            !$acc loop vector reduction(+:sum1) private(vol_liq_so,h2osoi_liq_so,h2osoi_liq_sat,deficit)
            do j = 1,nlevgrnd
               ! if level L was frozen, then we don't look at any levels below L
               if (t_soisno(c,j) > SHR_CONST_TKFRZ .and. rootfr(p,j) > 0._r8) then
                  ! determine soil water deficit in this layer:
                  ! Calculate vol_liq_so - i.e., vol_liq at which smp_node = smpso - by inverting the above equations
                  ! for the root resistance factors
                  vol_liq_so   = eff_porosity(c,j) * (-smpso(veg_pp%itype(p))/sucsat(c,j))**(-1/bsw(c,j))

                  ! Translate vol_liq_so and eff_porosity into h2osoi_liq_so and h2osoi_liq_sat and calculate deficit
                  h2osoi_liq_so  = vol_liq_so * denh2o * col_pp%dz(c,j)
                  h2osoi_liq_sat = eff_porosity(c,j) * denh2o * col_pp%dz(c,j)
                  deficit        = max((h2osoi_liq_so + firrig(g,tpu_ind)*(h2osoi_liq_sat - h2osoi_liq_so)) - h2osoi_liq(c,j), 0._r8)

                  ! Add deficit to irrig_rate, converting units from mm to mm/sec
                  sum1  = sum1 + deficit/(dtime_mod*irrig_nsteps_per_day)

               end if  ! else if (rootfr(p,j) > 0)
            end do     ! do j
            irrig_rate(p) = sum1
         end if        ! if (check_for_irrig(f) .and. .not. frozen_soil(f))
      end do           ! do f

      found = .false.

      ! Modify aerodynamic parameters for sparse/dense canopy (X. Zeng)
      !$acc parallel loop independent gang vector default(present) private(p,c,egvf,lt) &
      !$acc present(z0qv(:),forc_hgt_u_patch(:),elai(:),displa(:),esai(:),z0hv(:),z0mg(:),z0mv(:))
      do filter_index = 1, fn
         p = filter_nolu_vegp(filter_index)
         c = veg_pp%column(p)

         lt = min(elai(p)+esai(p), tlsai_crit)
         egvf =(1._r8 - alpha_aero * exp(-lt)) / (1._r8 - alpha_aero * exp(-tlsai_crit))
         displa(p) = egvf * displa(p)
         z0mv(p)   = exp(egvf * log(z0mv(p)) + (1._r8 - egvf) * log(z0mg(c)))
         z0hv(p)   = z0mv(p)
         z0qv(p)   = z0mv(p)

         !!Moved this here to allow async compute/data create w/ loop below
         zldis(filter_index) = forc_hgt_u_patch(p) - displa(p)

      end do

      !$acc enter data create(air(:),bir(:), cir(:), co2(:),o2(:),&
      !$acc nmozsgn(:), taf(:),qaf(:), ur(:),dth(:),dqh(:),delq(:), &
      !$acc  dthv(:), obu(:),el(:),qsatl(:),qsatldT(:), um(:), &
      !$acc  wta0(:), err(:), det(:))

      !$acc parallel loop independent gang vector default(present) private(p,c,t,g,deldT) present(thm(:),emv(:))
      do filter_index = 1, fn
         p = filter_nolu_vegp(filter_index)
         c = veg_pp%column(p)
         t = veg_pp%topounit(p)
         g = veg_pp%gridcell(p)

         ! Net absorbed longwave radiation by canopy and ground
         ! =air+bir*t_veg**4+cir*t_grnd(c)**4
         air(filter_index) =   emv(p) * (1._r8+(1._r8-emv(p))*(1._r8-emg(c))) * forc_lwrad(t)
         bir(filter_index) = - (2._r8-emv(p)*(1._r8-emg(c))) * emv(p) * sb
         cir(filter_index) =   emv(p)*emg(c)*sb

         if (use_finetop_rad) then
            slope_rad = slope_deg(g) * deg2rad
            bir(p) = bir(p) / cos(slope_rad)
            cir(p) = cir(p) / cos(slope_rad)
         endif
         ! Saturated vapor pressure, specific humidity, and their derivatives
         ! at the leaf surface

         call QSat (t_veg(p), forc_pbot(t), el(filter_index), deldT, qsatl(filter_index), qsatldT(filter_index))

         ! Determine atmospheric co2 and o2

         co2(filter_index) = forc_pco2(t)
         o2(filter_index)  = forc_po2(t)

         ! Initialize flux profile
         nmozsgn(filter_index) = 0

         taf(filter_index) = (t_grnd(c) + thm(p))/2._r8
         qaf(filter_index) = (forc_q(t)+qg(c))/2._r8

         ugust_total(filter_index) = ugust(t)

         ur(filter_index)    = max(1.0_r8,sqrt(forc_u(t)*forc_u(t)+forc_v(t)*forc_v(t)))
         dth(filter_index)   = thm(p)-taf(filter_index)
         dqh(filter_index)   = forc_q(t)-qaf(filter_index)
         delq(filter_index)  = qg(c) - qaf(filter_index)
         dthv(filter_index)  = dth(filter_index)*(1._r8+0.61_r8*forc_q(t))+0.61_r8*forc_th(t)*dqh(filter_index)

      end do

      if (found) then
         if ( .not. use_fates ) then
            write(*,*)'Error: Forcing height is below canopy height for pft index '
            call endrun(decomp_index=index, elmlevel=namep, msg=errmsg(__FILE__, __LINE__))
         end if
      end if
      
      !$acc parallel loop independent gang vector default(present) private(p,c)
      do filter_index = 1, fn
         p = filter_nolu_vegp(filter_index)
         c = veg_pp%column(p)

         ! Initialize Monin-Obukhov length and wind speed
         num_iter(p) = 0._r8
         call MoninObukIni(ur(filter_index), thv(c), dthv(filter_index), zldis(filter_index), z0mv(p), um(filter_index), obu(filter_index))

      end do

      ! Set counter for leaf temperature iteration (itlef)
      itlef = 0
      !$acc enter data copyin(itlef) create(temp1(:), temp2(:),temp12m(:),&
      !$acc    temp22m(:),ustar(:),rah(:,:),raw(:,:), uaf(:),rb(:), &
      !$acc     tlbef(:), del(:),del2(:),svpts(:),eah(:), dt_veg(:),wtg(:), &
      !$acc     wtl0(:), wtal(:), efe(:), dele(:),fm(:))
      
      ! Begin stability iteration
      call t_startf('can_iter')
      ITERATION : do while (itlef <= itmax .and. fn > 0)
        !$acc update device(itlef)  
         call FrictionVelocity (begp, endp, fn, iter_filterp, iter_filter_map, num_nolu_vegp, &
              displa(begp:endp), z0mv(begp:endp), z0hv(begp:endp), z0qv(begp:endp), &
              obu(1:num_nolu_vegp), itlef, ur(1:num_nolu_vegp), um(1:num_nolu_vegp), ugust_total(1:num_nolu_vegp), ustar(1:num_nolu_vegp), &
              temp1(1:num_nolu_vegp), temp2(1:num_nolu_vegp), temp12m(1:num_nolu_vegp), temp22m(1:num_nolu_vegp), fm(1:num_nolu_vegp), &
              frictionvel_vars)

        !$acc parallel loop independent gang vector default(present) private(p,filter_index,c,t,g,&
        !$acc  cf, w,csoilb,ri, ricsoilc, csoilcn) present(ram1(:), rb1(:), rhaf(:),grnd_ch4_cond(:),t_veg(:),elai(:),btran(:),&
        !$acc  esai(:), temp2(:), htop(:), dleaf_patch(:), rah(:,:))
        do active_index = 1, fn
           p = iter_filterp(active_index)
           filter_index = iter_filter_map(active_index)

           c = veg_pp%column(p)
           t = veg_pp%topounit(p)
           g = veg_pp%gridcell(p)

            tlbef(filter_index) = t_veg(p) !not used right now?
            del2(filter_index) = del(filter_index)   ! also not used in this loop

            ! Determine aerodynamic resistances
            ram1(p)  = 1._r8/(ustar(filter_index)*ustar(filter_index)/um(filter_index))
            rah(filter_index,above_canopy) = 1._r8/(temp1(filter_index)*ustar(filter_index))
            raw(filter_index,above_canopy) = 1._r8/(temp2(filter_index)*ustar(filter_index))

            ! Forbid removing more than 99% of wind speed in a time step.
            ! This is mainly to avoid convergence issues since this is such a
            ! basic form of iteration in this loop...
            if (implicit_stress) then
               tau(p) = forc_rho(t)*wind_speed_adj(p)/ram1(p)
               call shr_flux_update_stress(wind_speed0(p), wsresp(t), tau_est(t), &
                    tau(p), prev_tau(p), tau_diff(p), prev_tau_diff(p), &
                    wind_speed_adj(p))
               ur(p) = max(1.0_r8, sqrt(wind_speed_adj(p)**2 + ugust(t)**2))
            end if

            ! Bulk boundary layer resistance of leaves
            uaf(filter_index) = um(filter_index)*sqrt( 1._r8/(ram1(p)*um(filter_index)) )

            ! Use pft parameter for leaf characteristic width
            ! dleaf_patch if this is not an ed patch.
            ! Otherwise, the value has already been loaded
            ! during the FATES dynamics and/or initialization call
            if(.not.veg_pp%is_fates(p)) then
               dleaf_patch(p) = dleaf(veg_pp%itype(p))
            end if


            cf  = 0.01_r8/(sqrt(uaf(filter_index))*sqrt( dleaf_patch(p) ))
            rb(filter_index)  = 1._r8/(cf*uaf(filter_index))
            rb1(p) = rb(filter_index) !NOTE: this doesn't need to be updated every iteration

            ! Parameterization for variation of csoilc with canopy density from
            ! X. Zeng, University of Arizona

            w = exp(-(elai(p)+esai(p)))

            ! changed by K.Sakaguchi from here
            ! transfer coefficient over bare soil is changed to a local variable
            ! just for readability of the code (from line 680)
            csoilb = (vkc/(0.13_r8*(z0mg(c)*uaf(filter_index)/1.5e-5_r8)**0.45_r8))

            !compute the stability parameter for ricsoilc  ("S" in Sakaguchi&Zeng,2008)

            ri = ( grav*htop(p) * (taf(filter_index) - t_grnd(c)) ) / (taf(filter_index) * uaf(filter_index) **2.00_r8)

            !! modify csoilc value (0.004) if the under-canopy is in stable condition

            if ( (taf(filter_index) - t_grnd(c) ) > 0._r8) then
               ! decrease the value of csoilc by dividing it with (1+gamma*min(S, 10.0))
               ! ria ("gmanna" in Sakaguchi&Zeng, 2008) is a constant (=0.5)
               ricsoilc = csoilc / (1.00_r8 + ria*min( ri, 10.0_r8) )
               csoilcn = csoilb*w + ricsoilc*(1._r8-w)
            else
               csoilcn = csoilb*w + csoilc*(1._r8-w)
            end if

            !! Sakaguchi changes for stability formulation ends here

            rah(filter_index,below_canopy) = 1._r8/(csoilcn*uaf(filter_index))
            raw(filter_index,below_canopy) = rah(filter_index,below_canopy)
            if (use_lch4) then
               grnd_ch4_cond(p) = 1._r8/(raw(filter_index,above_canopy)+raw(filter_index,below_canopy))
            end if

            ! Stomatal resistances for sunlit and shaded fractions of canopy.
            ! Done each iteration to account for differences in eah, tv.

            svpts(filter_index) = el(filter_index)                         ! pa
            eah(filter_index) = forc_pbot(t) * qaf(filter_index) / 0.622_r8   ! pa
            rhaf(p) = eah(filter_index)/svpts(filter_index)

         ! Modification for shrubs proposed by X.D.Z
         ! Equivalent modification for soy following AgroIBIS
         ! NOTE: the following block of code was moved out of Photosynthesis subroutine and
         ! into here by M. Vertenstein on 4/6/2014 as part of making the photosynthesis
         ! routine a separate module. This move was also suggested by S. Levis in the previous
         ! version of the code.
         ! BUG MV 4/7/2014 - is this the correct place to have it in the iteration?
         ! THIS SHOULD BE MOVED OUT OF THE ITERATION but will change answers -

            if(.not.veg_pp%is_fates(p)) then
               if (crop(veg_pp%itype(p)) >= 1 .and. nfixer(veg_pp%itype(p)) == 1) then
                  btran(p) = min(1._r8, btran(p) * 1.25_r8)
               end if
            end if

        end do

         if ( use_fates ) then
            call alm_fates%wrap_photosynthesis(bounds, fn, iter_filterp(1:fn), &
                  svpts(begp:endp), eah(begp:endp), o2(begp:endp), &
                  co2(begp:endp), rb(begp:endp), dayl_factor(begp:endp), &
                  atm2lnd_vars, canopystate_vars, photosyns_vars)
         else ! not use_fates

            if ( use_hydrstress ) then
               call PhotosynthesisHydraulicStress (bounds, fn, iter_filterp, &
                    svpts(begp:endp), eah(begp:endp), o2(begp:endp), co2(begp:endp), rb(begp:endp), bsun(begp:endp), &
                    bsha(begp:endp), btran(begp:endp), dayl_factor(begp:endp), &
                    qsatl(begp:endp), qaf(begp:endp),     &
                     soilstate_vars, surfalb_vars, solarabs_vars,    &
                    canopystate_vars, photosyns_vars)
            else
              call Photosynthesis(bounds,fn,iter_filterp,iter_filter_map,num_nolu_vegp,&
                        svpts(1:num_nolu_vegp), eah(1:num_nolu_vegp),o2(1:num_nolu_vegp),&
                        co2(1:num_nolu_vegp), rb(1:num_nolu_vegp), btran(begp:endp), dayl_factor(1:num_nolu_vegp),&
                        surfalb_vars, solarabs_vars, canopystate_vars, photosyns_vars, 'sun', &
                        solarabs_vars%parsun_z_patch(begp:endp,:),  canopystate_vars%laisun_z_patch(begp:endp,:), &
                        surfalb_vars%vcmaxcintsun_patch(begp:endp),  photosyns_vars%alphapsnsun_patch(begp:endp), &
                        photosyns_vars%cisun_z_patch(begp:endp,:), photosyns_vars%rssun_patch(begp:endp), &
                        photosyns_vars%rssun_z_patch(begp:endp,:), photosyns_vars%lmrsun_patch(begp:endp), &
                        photosyns_vars%lmrsun_z_patch(begp:endp,:), photosyns_vars%psnsun_patch(begp:endp), &
                        photosyns_vars%psnsun_z_patch(begp:endp,:),photosyns_vars%psnsun_wc_patch(begp:endp), &
                        photosyns_vars%psnsun_wj_patch(begp:endp),photosyns_vars%psnsun_wp_patch(begp:endp)   )

            end if

            if ( use_c13 ) then
               call Fractionation (bounds, fn, iter_filterp, &
                     cnstate_vars, solarabs_vars, surfalb_vars, photosyns_vars, 1)
            endif

            !$acc parallel loop independent gang vector default(present) private(p,c)
            do active_index = 1, fn
               p = iter_filterp(active_index)
               c = veg_pp%column(p)
               ! soybean (crop with N fixation)
               if (crop(veg_pp%itype(p)) >= 1 .and. nfixer(veg_pp%itype(p)) == 1) then
                  btran(p) = min(1._r8, btran(p) * 1.25_r8)
               end if
            end do

            if ( .not. use_hydrstress ) then
               call Photosynthesis(bounds,fn,iter_filterp,iter_filter_map,num_nolu_vegp, &
                        svpts(1:num_nolu_vegp), eah(1:num_nolu_vegp),o2(1:num_nolu_vegp),&
                        co2(1:num_nolu_vegp),rb(1:num_nolu_vegp), btran(begp:endp), dayl_factor(1:num_nolu_vegp),&
                        surfalb_vars, solarabs_vars, canopystate_vars, photosyns_vars, 'sha', &
                        solarabs_vars%parsha_z_patch(begp:endp,:), canopystate_vars%laisha_z_patch(begp:endp,:), &
                        surfalb_vars%vcmaxcintsha_patch(begp:endp), photosyns_vars%alphapsnsha_patch(begp:endp), &
                        photosyns_vars%cisha_z_patch(begp:endp,:),photosyns_vars%rssha_patch(begp:endp), &
                        photosyns_vars%rssha_z_patch(begp:endp,:),photosyns_vars%lmrsha_patch(begp:endp), &
                        photosyns_vars%lmrsha_z_patch(begp:endp,:),photosyns_vars%psnsha_patch(begp:endp),&
                        photosyns_vars%psnsha_z_patch(begp:endp,:),photosyns_vars%psnsha_wc_patch(begp:endp),&
                        photosyns_vars%psnsha_wj_patch(begp:endp),photosyns_vars%psnsha_wp_patch(begp:endp)   )

            end if

            if ( use_c13 ) then
               call Fractionation (bounds, fn, iter_filterp,  &
                     cnstate_vars, solarabs_vars, surfalb_vars, photosyns_vars, 0)
            end if

         end if ! end of if use_fates

         !$acc parallel loop independent gang vector default(present) private(p,filter_index) present(laisun(:),&
         !$acc  thm(:), canopy_cond(:),temp2(:), frac_veg_nosno(:), esai(:), fdry(:), wta0(:), h2ocan(:), &
         !$acc  laisha(:),rssha(:),btran(:), fwet(:), qflx_evap_veg(:), qflx_tran_veg(:),sabv(:), eflx_sh_veg(:) )
         do active_index = 1, fn
            p = iter_filterp(active_index)
            filter_index = iter_filter_map(active_index)
            c = veg_pp%column(p)
            t = veg_pp%topounit(p)
            g = veg_pp%gridcell(p)

            ! Sensible heat conductance for air, leaf and ground
            ! Moved the original subroutine in-line...

            wta    = 1._r8/rah(filter_index,above_canopy)  ! air
            wtl    = (elai(p)+esai(p))/rb(filter_index)    ! leaf
            wtg(filter_index) = 1._r8/rah(filter_index,below_canopy)  ! ground
            wtshi  = 1._r8/(wta+wtl+wtg(filter_index))
            wtl0(filter_index) = wtl*wtshi         ! leaf
            wtg0    = wtg(filter_index)*wtshi      ! ground
            wta0(filter_index) = wta*wtshi         ! air

            wtga    = wta0(filter_index)+wtg0      ! ground + air
            wtal(filter_index) = wta0(filter_index)+wtl0(filter_index)   ! air + leaf

            ! Fraction of potential evaporation from leaf

            if (fdry(p) > 0._r8) then
               rppdry  = fdry(p)*rb(filter_index)*(laisun(p)/(rb(filter_index)+rssun(p)) + &
                    laisha(p)/(rb(filter_index)+rssha(p)))/elai(p)
            else
               rppdry = 0._r8
            end if

            ! Calculate canopy conductance for methane / oxygen (e.g. stomatal conductance & leaf bdy cond)
            if (use_lch4) then
               canopy_cond(p) = (laisun(p)/(rb(filter_index)+rssun(p)) + laisha(p)/(rb(filter_index)+rssha(p)))/max(elai(p), 0.01_r8)
            end if

            efpot = forc_rho(t)*wtl*(qsatl(filter_index)-qaf(filter_index))
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
                 rpp = min(rpp, (qflx_tran_veg(p)+h2ocan(p)/dtime_mod)/efpot)
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
               rpp = min(rpp, (qflx_tran_veg(p)+h2ocan(p)/dtime_mod)/efpot)
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

            wtaq    = frac_veg_nosno(p)/raw(filter_index,above_canopy)             ! air
            wtlq    = frac_veg_nosno(p)*(elai(p)+esai(p))/rb(filter_index) * rpp   ! leaf

            !Litter layer resistance. Added by K.Sakaguchi
            snow_depth_c = z_dl ! critical depth for 100% litter burial by snow (=litter thickness)
            fsno_dl = snow_depth(c)/snow_depth_c    ! effective snow cover for (dry)plant litter
            elai_dl = lai_dl*(1._r8 - min(fsno_dl,1._r8)) ! exposed (dry)litter area index
            rdl = ( 1._r8 - exp(-elai_dl) ) / ( 0.004_r8*uaf(filter_index)) ! dry litter layer resistance

            ! add litter resistance and Lee and Pielke 1992 beta
            if (delq(filter_index) < 0._r8) then  !dew. Do not apply beta for negative flux (follow old rsoil)
               wtgq(filter_index) = frac_veg_nosno(p)/(raw(filter_index,below_canopy)+rdl)
            else
               if (do_soilevap_beta()) then
                  wtgq(filter_index) = soilbeta(c)*frac_veg_nosno(p)/(raw(filter_index,below_canopy)+rdl)
               endif
            end if

            wtsqi   = 1._r8/(wtaq+wtlq+wtgq(filter_index))

            wtgq0    = wtgq(filter_index)*wtsqi      ! ground
            wtlq0(filter_index) = wtlq*wtsqi         ! leaf
            wtaq0(filter_index) = wtaq*wtsqi         ! air

            wtgaq    = wtaq0(filter_index)+wtgq0     ! air + ground
            wtalq(filter_index) = wtaq0(filter_index)+wtlq0(filter_index)  ! air + leaf

            dc1 = forc_rho(t)*cpair*wtl
            dc2 = hvap*forc_rho(t)*wtlq

            efsh   = dc1*(wtga*t_veg(p)-wtg0*t_grnd(c)-wta0(filter_index)*thm(p))
            efe(filter_index) = dc2*(wtgaq*qsatl(filter_index)-wtgq0*qg(c)-wtaq0(filter_index)*forc_q(t))

            ! Evaporation flux from foliage

            erre = 0._r8
            if (efe(filter_index)*efeb(filter_index) < 0._r8) then
               efeold = efe(filter_index)
               efe(filter_index)  = 0.1_r8*efeold
               erre = efe(filter_index) - efeold
            end if
            ! fractionate ground emitted longwave
            lw_grnd=(frac_sno(c)*t_soisno(c,snl(c)+1)**4 &
                 +(1._r8-frac_sno(c)-frac_h2osfc(c))*t_soisno(c,1)**4 &
                 +frac_h2osfc(c)*t_h2osfc(c)**4)

            dt_veg(filter_index) = (sabv(p) + air(filter_index) + bir(filter_index)*t_veg(p)**4 + &
                 cir(filter_index)*lw_grnd - efsh - efe(filter_index)) / &
                 (- 4._r8*bir(filter_index)*t_veg(p)**3 +dc1*wtga +dc2*wtgaq*qsatldT(filter_index))
            t_veg(p) = tlbef(filter_index) + dt_veg(filter_index)
            dels = dt_veg(filter_index)
            del(filter_index)  = abs(dels)
            err(filter_index) = 0._r8
            if (del(filter_index) > delmax) then
               dt_veg(filter_index) = delmax*dels/del(filter_index)
               t_veg(p) = tlbef(filter_index) + dt_veg(filter_index)
               err(filter_index) = sabv(p) + air(filter_index) + bir(filter_index)*tlbef(filter_index)**3*(tlbef(filter_index) + &
                    4._r8*dt_veg(filter_index)) + cir(filter_index)*lw_grnd - &
                    (efsh + dc1*wtga*dt_veg(filter_index)) - (efe(filter_index) + &
                    dc2*wtgaq*qsatldT(filter_index)*dt_veg(filter_index))
            end if

            ! Fluxes from leaves to canopy space
            ! "efe" was limited as its sign changes frequently.  This limit may
            ! result in an imbalance in "hvap*qflx_evap_veg" and
            ! "efe + dc2*wtgaq*qsatdt_veg"

            efpot = forc_rho(t)*wtl*(wtgaq*(qsatl(filter_index)+qsatldT(filter_index)*dt_veg(filter_index)) &
                 -wtgq0*qg(c)-wtaq0(filter_index)*forc_q(t))
            qflx_evap_veg(p) = rpp*efpot

            ! Calculation of evaporative potentials (efpot) and
            ! interception losses; flux in kg m**-2 s-1.  ecidif
            ! holds the excess energy if all intercepted water is evaporated
            ! during the timestep.  This energy is later added to the
            ! sensible heat flux.
            if ( use_hydrstress ) then
               ecidif = max(0._r8,qflx_evap_veg(p)-qflx_tran_veg(p)-h2ocan(p)/dtime_mod)
               qflx_evap_veg(p) = min(qflx_evap_veg(p),qflx_tran_veg(p)+h2ocan(p)/dtime_mod)
            else

              ecidif = 0._r8
              if (efpot > 0._r8 .and. btran(p) > btran0) then
               qflx_tran_veg(p) = efpot*rppdry
              else
               qflx_tran_veg(p) = 0._r8
              end if
              ecidif = max(0._r8, qflx_evap_veg(p)-qflx_tran_veg(p)-h2ocan(p)/dtime_mod)
              qflx_evap_veg(p) = min(qflx_evap_veg(p),qflx_tran_veg(p)+h2ocan(p)/dtime_mod)
            end if

            ! The energy loss due to above two limits is added to
            ! the sensible heat flux.
            eflx_sh_veg(p) = efsh + dc1*wtga*dt_veg(filter_index) + err(filter_index) + erre + hvap*ecidif

            ! Re-calculate saturated vapor pressure, specific humidity, and their
            ! derivatives at the leaf surface

            call QSat(t_veg(p), forc_pbot(t), el(filter_index), deldT, qsatl(filter_index), qsatldT(filter_index))

            ! Update vegetation/ground surface temperature, canopy air
            ! temperature, canopy vapor pressure, aerodynamic temperature, and
            ! Monin-Obukhov stability parameter for next iteration.

            taf(filter_index) = wtg0*t_grnd(c) + wta0(filter_index)*thm(p) + wtl0(filter_index)*t_veg(p)
            qaf(filter_index) = wtlq0(filter_index)*qsatl(filter_index) + wtgq0*qg(c) + forc_q(t)*wtaq0(filter_index)

            ! Update Obukhov length scale and wind speed including the
            ! stability effect

            dth(filter_index) = thm(p)-taf(filter_index)
            dqh(filter_index) = forc_q(t)-qaf(filter_index)
            delq(filter_index) = wtalq(filter_index)*qg(c)-wtlq0(filter_index)*qsatl(filter_index)-wtaq0(filter_index)*forc_q(t)

            tstar = temp1(filter_index)*dth(filter_index)
            qstar = temp2(filter_index)*dqh(filter_index)

            thvstar = tstar*(1._r8+0.61_r8*forc_q(t)) + 0.61_r8*forc_th(t)*qstar

            zeta(p) = zldis(filter_index)*vkc*grav*thvstar/(ustar(filter_index)**2*thv(c))
            if (zeta(p) >= 0._r8) then     !stable
               zeta(p) = min(2._r8,max(zeta(p),0.01_r8))
               um(filter_index) = max(ur(filter_index),0.1_r8)
            else                     !unstable
               zeta(p) = max(-100._r8,min(zeta(p),-0.01_r8))
               if ((.not. atm_gustiness) .or. force_land_gustiness) then
                  wc = beta*(-grav*ustar(filter_index)*thvstar*zii/thv(c))**0.333_r8
                  ugust_total(filter_index) = sqrt(ugust(t)**2 + wc**2)
                  um(filter_index) = sqrt(ur(filter_index)*ur(filter_index)+wc*wc)
               else
                  um(filter_index) = max(ur(filter_index),0.1_r8)
               end if
            end if
            obu(filter_index) = zldis(filter_index)/zeta(p)

            if (obuold(filter_index)*obu(filter_index) < 0._r8) nmozsgn(filter_index) = nmozsgn(filter_index)+1
            if (nmozsgn(filter_index) >= 4) obu(filter_index) = zldis(filter_index)/(-0.01_r8)
            obuold(filter_index) = obu(filter_index)

         end do   ! end of filtered pft loop

         !$acc parallel loop independent gang vector default(present) private(p,filter_index,t)
         do active_index = 1, fn
           p = iter_filterp(active_index)
           filter_index = iter_filter_map(active_index)
           t = veg_pp%topounit(p)
           !laminar boundary resistance for h2o over leaf, should I make this consistent for latent heat calculation?
           lbl_rsc_h2o(p) = getlblcef(forc_rho(t),t_veg(p))*uaf(filter_index)/(uaf(filter_index)**2._r8+1.e-10_r8)   
         enddo

         ! Test for convergence.
         ! Compact iter_filterp/iter_filter_map in place, keeping only the patches that have
         ! NOT yet converged; fn shrinks to the new active count. This loop has a sequential
         ! write-index dependency (fn_new), so it must run in order.
         itlef = itlef+1
         if (itlef > itmin) then
            fnold = fn
            fn = 0
            !$acc parallel loop seq default(present) private(p,filter_index) present(det(1:fnold), dele(1:fnold))
            do active_index = 1, fnold
               p = iter_filterp(active_index)
               filter_index = iter_filter_map(active_index)
               num_iter(p) = real(itlef,r8)
               dele(filter_index) = abs(efe(filter_index) - efeb(filter_index))
               efeb(filter_index) = efe(filter_index)
               det(filter_index)  = max(del(filter_index),del2(filter_index))
               if (.not. (det(filter_index) < dtmin .and. dele(filter_index) < dlemin)) then
                  ! still unconverged: keep it in the compacted filter for the next iteration
                  fn = fn + 1
                  iter_filterp(fn)     = p
                  iter_filter_map(fn)  = filter_index
               end if
            end do
         end if
      end do ITERATION     ! End stability iteration

      call t_stopf('can_iter')
      
      !$acc parallel loop independent gang vector default(present)
      do filter_index = 1, num_nolu_vegp
         p = filter_nolu_vegp(filter_index)
         c = veg_pp%column(p)
         t = veg_pp%topounit(p)
         g = veg_pp%gridcell(p)

         ! Energy balance check in canopy

         lw_grnd=(frac_sno(c)*t_soisno(c,snl(c)+1)**4 &
              +(1._r8-frac_sno(c)-frac_h2osfc(c))*t_soisno(c,1)**4 &
              +frac_h2osfc(c)*t_h2osfc(c)**4)

         err(filter_index) = sabv(p) + air(filter_index) + bir(filter_index)*tlbef(filter_index)**3*(tlbef(filter_index) + 4._r8*dt_veg(filter_index)) &
              + cir(filter_index)*lw_grnd - eflx_sh_veg(p) - hvap*qflx_evap_veg(p)

         ! Fluxes from ground to canopy space

         delt    = wtal(filter_index)*t_grnd(c)-wtl0(filter_index)*t_veg(p)-wta0(filter_index)*thm(p)
         taux(p) = -forc_rho(t)*forc_u(t)/ram1(p)
         tauy(p) = -forc_rho(t)*forc_v(t)/ram1(p)
         eflx_sh_grnd(p) = cpair*forc_rho(t)*wtg(filter_index)*delt

         ! compute individual sensible heat fluxes
         delt_snow = wtal(filter_index)*t_soisno(c,snl(c)+1)-wtl0(filter_index)*t_veg(p)-wta0(filter_index)*thm(p)
         eflx_sh_snow(p) = cpair*forc_rho(t)*wtg(filter_index)*delt_snow

         delt_soil  = wtal(filter_index)*t_soisno(c,1)-wtl0(filter_index)*t_veg(p)-wta0(filter_index)*thm(p)
         eflx_sh_soil(p) = cpair*forc_rho(t)*wtg(filter_index)*delt_soil

         delt_h2osfc  = wtal(filter_index)*t_h2osfc(c)-wtl0(filter_index)*t_veg(p)-wta0(filter_index)*thm(p)
         eflx_sh_h2osfc(p) = cpair*forc_rho(t)*wtg(filter_index)*delt_h2osfc
         qflx_evap_soi(p) = forc_rho(t)*wtgq(filter_index)*delq(filter_index)

         ! compute individual latent heat fluxes
         delq_snow = wtalq(filter_index)*qg_snow(c)-wtlq0(filter_index)*qsatl(filter_index)-wtaq0(filter_index)*forc_q(t)
         qflx_ev_snow(p) = forc_rho(t)*wtgq(filter_index)*delq_snow

         delq_soil = wtalq(filter_index)*qg_soil(c)-wtlq0(filter_index)*qsatl(filter_index)-wtaq0(filter_index)*forc_q(t)
         qflx_ev_soil(p) = forc_rho(t)*wtgq(filter_index)*delq_soil

         delq_h2osfc = wtalq(filter_index)*qg_h2osfc(c)-wtlq0(filter_index)*qsatl(filter_index)-wtaq0(filter_index)*forc_q(t)
         qflx_ev_h2osfc(p) = forc_rho(t)*wtgq(filter_index)*delq_h2osfc

         ! 2 m height air temperature

         t_ref2m(p) = thm(p) + temp1(filter_index)*dth(filter_index)*(1._r8/temp12m(filter_index) - 1._r8/temp1(filter_index))
         t_ref2m_r(p) = t_ref2m(p)

         ! 2 m height specific humidity

         q_ref2m(p) = forc_q(t) + temp2(filter_index)*dqh(filter_index)*(1._r8/temp22m(filter_index) - 1._r8/temp2(filter_index))

         ! 2 m height relative humidity

         call QSat(t_ref2m(p), forc_pbot(t), e_ref2m, de2mdT, qsat_ref2m, dqsat2mdT)
         rh_ref2m(p) = min(100._r8, q_ref2m(p) / qsat_ref2m * 100._r8)
         rh_ref2m_r(p) = rh_ref2m(p)


         if (use_finetop_rad) then
            slope_rad = slope_deg(g) * deg2rad

            ! Downward longwave radiation below the canopy
            dlrad(p) = (1._r8-emv(p))*emg(c)*forc_lwrad(t) + &
                  emv(p)*emg(c)*sb*tlbef(p)**3*(tlbef(filter_index) + 4._r8*dt_veg(filter_index))/cos(slope_rad)

            ! Upward longwave radiation above the canopy
            ulrad(p) = ((1._r8-emg(c))*(1._r8-emv(p))*(1._r8-emv(p))*forc_lwrad(t) &
                + emv(p)*(1._r8+(1._r8-emg(c))*(1._r8-emv(p)))*sb*tlbef(filter_index)**3*(tlbef(filter_index) + &
                4._r8*dt_veg(filter_index))/cos(slope_rad) + emg(c)*(1._r8-emv(p))*sb*lw_grnd/cos(slope_rad))
         else
            dlrad(p) = (1._r8-emv(p))*emg(c)*forc_lwrad(t) + &
                  emv(p)*emg(c)*sb*tlbef(p)**3*(tlbef(p) + 4._r8*dt_veg(p))

            ulrad(p) = ((1._r8-emg(c))*(1._r8-emv(p))*(1._r8-emv(p))*forc_lwrad(t) &
                + emv(p)*(1._r8+(1._r8-emg(c))*(1._r8-emv(p)))*sb*tlbef(filter_index)**3*(tlbef(filter_index) + &
                4._r8*dt_veg(filter_index)) + emg(c)*(1._r8-emv(p))*sb*lw_grnd)
         endif

         ! Derivative of soil energy flux with respect to soil temperature

         cgrnds(p) = cgrnds(p) + cpair*forc_rho(t)*wtg(filter_index)*wtal(filter_index)
         cgrndl(p) = cgrndl(p) + forc_rho(t)*wtgq(filter_index)*wtalq(filter_index)*dqgdT(c)
         cgrnd(p)  = cgrnds(p) + cgrndl(p)*htvp(c)

         ! Update dew accumulation (kg/m2)

         h2ocan(p) = max(0._r8,h2ocan(p)+(qflx_tran_veg(p)-qflx_evap_veg(p))*dtime_mod)

         ! Check for convergence of stress.
         if (implicit_stress .and. abs(tau_diff(p)) > dtaumin) then
            if (nstep_mod > 0) then ! Suppress common warnings on the first time step.
               write(iulog,*)'WARNING: Stress did not converge for canopy ',&
                    ' nstep = ',nstep_mod,' p= ',p,' prev_tau_diff= ',prev_tau_diff(p),&
                    ' tau_diff= ',tau_diff(p),' tau= ',tau(p),&
                    ' wind_speed_adj= ',wind_speed_adj(p),' iter_final= ',itlef
            end if
         end if

      end do
            ! variables for history fields
            rah_above(p)  = rah(p,above_canopy)
            raw_above(p)  = raw(p,above_canopy)
            rah_below(p)  = rah(p,below_canopy)
            raw_below(p)  = raw(p,below_canopy)
            vpd(p)        = max((svpts(p) - eah(p)), vpd_min) * pa_to_kpa ! kPa

      if ( use_fates ) then

        call alm_fates%wrap_accumulatefluxes(bounds,num_nolu_vegp,filter_nolu_vegp(1:num_nolu_vegp))
        call alm_fates%wrap_hydraulics_drive(bounds,num_nolu_vegp,filter_nolu_vegp(1:num_nolu_vegp),soilstate_vars, &
                                            solarabs_vars,energyflux_vars)
      else

         ! Determine total photosynthesis
         call PhotosynthesisTotal(num_nolu_vegp, filter_nolu_vegp, &
               cnstate_vars, canopystate_vars, photosyns_vars)
         ! Filter out patches which have small energy balance errors; report others
         ! NOTE: filter out patches for what? This is the end of the subroutine.
         fnold = num_nolu_vegp
         fn = 0
         !$acc parallel loop independent gang vector default(present) 
         do filter_index = 1, fnold
            p = filter_nolu_vegp(filter_index)
            if (abs(err(filter_index)) > 0.1_r8) then
               fn = fn + 1
               iter_filterp(fn) = p
               write(iulog,*) 'energy balance in canopy ',p,', err=',err(filter_index)
               write(iulog,*) "sabv  :", sabv(p) 
               write(iulog,*) "air   :",air(p)
               write(iulog,*) "bir   :" ,bir(p)
               write(iulog,*) "cir   :" ,cir(p)
               write(iulog,*) "tlbef :",tlbef(p)
               write(iulog,*) "dt_veg:",dt_veg(p) 
               write(iulog,*) "eflx_sh_veg:",eflx_sh_veg(p) 
               write(iulog,*) "qflx_evap_veg:",qflx_evap_veg(p)
            end if
         end do

      end if
      !$acc exit data delete(del(:), efeb(:), wtlq0(:),wtalq(:), &
      !$acc  wtgq(:), wtaq0(:), obuold(:),dayl_factor(:) , &
      !$acc  check_for_irrig(:), filterp(:),zldis(:), &
      !$acc  air(:),bir(:), cir(:), co2(:),o2(:),&
      !$acc  nmozsgn(:), taf(:),qaf(:), ur(:),dth(:),dqh(:),delq(:), &
      !$acc  dthv(:), obu(:),el(:),qsatl(:),qsatldT(:), &
      !$acc  temp1(:), temp2(:),temp12m(:),&
      !$acc  temp22m(:),ustar(:), um(:),rah(:,:),raw(:,:), uaf(:),rb(:), &
      !$acc  tlbef(:), del(:),del2(:),svpts(:),eah(:),wta0(:), err(:), dt_veg(:) ,wtg(:), &
      !$acc  wtal(:), wtl0(:), efe(:), det(:), dele(:), fm(:)  )
      !$acc exit data delete(time,irrig_nsteps_per_day, itlef) 
    end associate

  end subroutine CanopyFluxes

end module CanopyFluxesMod
