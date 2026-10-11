module  PhotosynthesisMod

#include "shr_assert.h"
  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Leaf photosynthesis and stomatal conductance calculation as described by
  ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 and extended to
  ! a multi-layer canopy.
  !
  ! The sunlit and shaded canopy fractions are solved separately: CanopyFluxes
  ! calls Photosynthesis once per phase (sunlit_phase / shaded_phase).
  !
  ! The plant hydraulic stress (PHS) variant, which solves sunlit and shaded
  ! leaves simultaneously, lives in PhotosynthesisHydraulicStressMod.
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use abortutils          , only : endrun
  use elm_varctl          , only : iulog, use_c13, use_c14, use_cn, use_fates
  use elm_varctl          , only : carbon_only, carbonnitrogen_only, carbonphosphorus_only
  use elm_varpar          , only : nlevcan
  use elm_varcon          , only : namep, rgas, tfrz
  use decompMod           , only : bounds_type
  use QuadraticMod        , only : quadratic
  use VegetationPropertiesType      , only : veg_vp
  use CNStateType         , only : cnstate_type
  use CanopyStateType     , only : canopystate_type
  use SolarAbsorbedType   , only : solarabs_type
  use SurfaceAlbedoType   , only : surfalb_type
  use PhotosynthesisType  , only : photosyns_type
  use VegetationType      , only : veg_pp
  use AllocationMod       , only : nu_com_leaf_physiology
  use pftvarcon           , only : noveg
  use pftvarcon           , only : vcmax_np1, vcmax_np2, vcmax_np3, vcmax_np4, jmax_np1, jmax_np2, jmax_np3
  use SharedParamsMod     , only : ParamsShareInst
  use TopounitDataType    , only : top_as
  use VegetationDataType  , only : veg_es, veg_ns, veg_ps
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: Photosynthesis       ! Leaf stomatal resistance and leaf photosynthesis (one canopy phase)
  public :: PhotosynthesisTotal  ! Determine of total photosynthesis
  public :: Fractionation        ! C13 fractionation during photosynthesis
  public :: ft                   ! photosynthesis temperature response
  public :: fth                  ! photosynthesis temperature inhibition
  public :: fth25                ! scaling factor for photosynthesis temperature inhibition

  ! !PUBLIC DATA:
  ! Canopy phase selectors for Photosynthesis and Fractionation
  integer, parameter, public :: sunlit_phase = 1
  integer, parameter, public :: shaded_phase = 0

  ! !PRIVATE MEMBER FUNCTIONS:
  private :: PhaseLeafPhotosynthesis ! Photosynthesis for one canopy phase given its phase-specific arrays
  private :: canopy_top_capacity     ! vcmax25/jmax25 and leaf N at the canopy top
  private :: temperature_scaled      ! Apply Arrhenius response with high temperature inhibition to a 25C rate
  private :: hybrid                  ! hybrid solver for ci
  private :: ci_func                 ! ci function
  private :: brent                   ! brent solver for root of a single variable function

  ! !PRIVATE DATA:
  ! Leaf biochemistry constants, Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
  real(r8), parameter :: fnps       = 0.15_r8  ! fraction of light absorbed by non-photosynthetic pigments
  real(r8), parameter :: theta_psii = 0.7_r8   ! empirical curvature parameter for electron transport rate
  real(r8), parameter :: theta_ip   = 0.95_r8  ! empirical curvature parameter for ap photosynthesis co-limitation
  real(r8), parameter :: rsmax0     = 2.e4_r8  ! maximum stomatal resistance [s/m]
  real(r8), parameter :: par_to_photons = 4.6_r8 ! converts PAR from W/m2 to umol photons/m**2/s

  ! Leaf environment held fixed while the ci / stomatal conductance solvers iterate
  type :: leaf_env_type
     real(r8) :: pbot     ! atmospheric pressure (Pa)
     real(r8) :: gb_mol   ! leaf boundary layer conductance (umol H2O/m**2/s)
     real(r8) :: je       ! electron transport rate (umol electrons/m**2/s)
     real(r8) :: cair     ! atmospheric CO2 partial pressure (Pa)
     real(r8) :: oair     ! atmospheric O2 partial pressure (Pa)
     real(r8) :: lmr      ! canopy layer leaf maintenance respiration rate (umol CO2/m**2/s)
     real(r8) :: par      ! par absorbed per unit lai for canopy layer (w/m**2)
     real(r8) :: rh_can   ! canopy air relative humidity
  end type leaf_env_type

contains

  !------------------------------------------------------------------------------
  subroutine Photosynthesis ( bounds, fn, filterp, filter_map, norig, &
       esat_tv, eair, oair, cair, rb, btran, dayl_factor, &
       surfalb_vars, solarabs_vars, canopystate_vars, photosyns_vars, phase)
    !
    ! !DESCRIPTION:
    ! Leaf photosynthesis and stomatal conductance for the sunlit or shaded
    ! fraction of the canopy. Selects the phase-specific arrays and passes them
    ! explicitly to PhaseLeafPhotosynthesis (avoids deep copying derived-type
    ! pointers to the GPU).
    !
    ! Note: This subroutine is not called via FATES (RGK)
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: fn                  ! number of still-active patches
    integer                , intent(in)    :: filterp(1:fn)       ! patch filter (active, unconverged patches)
    integer                , intent(in)    :: filter_map(1:norig) ! active-filter position -> index into the norig-sized arrays
    integer                , intent(in)    :: norig               ! size of the original (uncompacted) per-patch arrays below
    real(r8)               , intent(in)    :: esat_tv(1:norig)    ! saturation vapor pressure at t_veg (Pa)
    real(r8)               , intent(in)    :: eair(1:norig)       ! vapor pressure of canopy air (Pa)
    real(r8)               , intent(in)    :: oair(1:norig)       ! atmospheric O2 partial pressure (Pa)
    real(r8)               , intent(in)    :: cair(1:norig)       ! atmospheric CO2 partial pressure (Pa)
    real(r8)               , intent(in)    :: rb(1:norig)         ! boundary layer resistance (s/m)
    real(r8)               , intent(in)    :: btran(bounds%begp:) ! transpiration wetness factor (0 to 1) [patch]
    real(r8)               , intent(in)    :: dayl_factor(1:norig)! scalar (0-1) for daylength
    type(surfalb_type)     , intent(inout) :: surfalb_vars
    type(solarabs_type)    , intent(inout) :: solarabs_vars
    type(canopystate_type) , intent(inout) :: canopystate_vars
    type(photosyns_type)   , intent(inout) :: photosyns_vars
    integer                , intent(in)    :: phase               ! sunlit_phase or shaded_phase
    !------------------------------------------------------------------------------

    associate(begp => bounds%begp, endp => bounds%endp)
      if (phase == sunlit_phase) then
         call PhaseLeafPhotosynthesis(bounds, fn, filterp, filter_map, norig, &
              esat_tv, eair, oair, cair, rb, btran, dayl_factor, surfalb_vars, photosyns_vars, &
              par_z     = solarabs_vars%parsun_z_patch(begp:endp,:), &
              lai_z     = canopystate_vars%laisun_z_patch(begp:endp,:), &
              vcmaxcint = surfalb_vars%vcmaxcintsun_patch(begp:endp), &
              alphapsn  = photosyns_vars%alphapsnsun_patch(begp:endp), &
              ci_z      = photosyns_vars%cisun_z_patch(begp:endp,:), &
              rs        = photosyns_vars%rssun_patch(begp:endp), &
              rs_z      = photosyns_vars%rssun_z_patch(begp:endp,:), &
              lmr       = photosyns_vars%lmrsun_patch(begp:endp), &
              lmr_z     = photosyns_vars%lmrsun_z_patch(begp:endp,:), &
              psn       = photosyns_vars%psnsun_patch(begp:endp), &
              psn_z     = photosyns_vars%psnsun_z_patch(begp:endp,:), &
              psn_wc    = photosyns_vars%psnsun_wc_patch(begp:endp), &
              psn_wj    = photosyns_vars%psnsun_wj_patch(begp:endp), &
              psn_wp    = photosyns_vars%psnsun_wp_patch(begp:endp))
      else
         call PhaseLeafPhotosynthesis(bounds, fn, filterp, filter_map, norig, &
              esat_tv, eair, oair, cair, rb, btran, dayl_factor, surfalb_vars, photosyns_vars, &
              par_z     = solarabs_vars%parsha_z_patch(begp:endp,:), &
              lai_z     = canopystate_vars%laisha_z_patch(begp:endp,:), &
              vcmaxcint = surfalb_vars%vcmaxcintsha_patch(begp:endp), &
              alphapsn  = photosyns_vars%alphapsnsha_patch(begp:endp), &
              ci_z      = photosyns_vars%cisha_z_patch(begp:endp,:), &
              rs        = photosyns_vars%rssha_patch(begp:endp), &
              rs_z      = photosyns_vars%rssha_z_patch(begp:endp,:), &
              lmr       = photosyns_vars%lmrsha_patch(begp:endp), &
              lmr_z     = photosyns_vars%lmrsha_z_patch(begp:endp,:), &
              psn       = photosyns_vars%psnsha_patch(begp:endp), &
              psn_z     = photosyns_vars%psnsha_z_patch(begp:endp,:), &
              psn_wc    = photosyns_vars%psnsha_wc_patch(begp:endp), &
              psn_wj    = photosyns_vars%psnsha_wj_patch(begp:endp), &
              psn_wp    = photosyns_vars%psnsha_wp_patch(begp:endp))
      end if
    end associate

  end subroutine Photosynthesis

  !------------------------------------------------------------------------------
  subroutine PhaseLeafPhotosynthesis ( bounds, fn, filterp, filter_map, norig, &
       esat_tv, eair, oair, cair, rb, btran, dayl_factor, surfalb_vars, photosyns_vars, &
       par_z, lai_z, vcmaxcint, alphapsn, ci_z, rs, rs_z, lmr, lmr_z, &
       psn, psn_z, psn_wc, psn_wj, psn_wp)
    !
    ! !DESCRIPTION:
    ! Photosynthesis for one canopy phase. Four steps, each one OpenACC loop:
    !   1. Patch kinetic parameters (C3/C4, Ball-Berry, kc/ko/cp at leaf temperature)
    !   2. Canopy-layer capacities (vcmax/jmax/tpu/kp/lmr) scaled by the leaf N profile
    !   3. Leaf-level ci / stomatal conductance solution per canopy layer
    !   4. Canopy sums -> effective leaf-level fluxes for this phase
    !
    ! !ARGUMENTS:
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: fn
    integer                , intent(in)    :: filterp(1:fn)
    integer                , intent(in)    :: filter_map(1:norig)
    integer                , intent(in)    :: norig
    real(r8)               , intent(in)    :: esat_tv(1:norig)
    real(r8)               , intent(in)    :: eair(1:norig)
    real(r8)               , intent(in)    :: oair(1:norig)
    real(r8)               , intent(in)    :: cair(1:norig)
    real(r8)               , intent(in)    :: rb(1:norig)
    real(r8)               , intent(in)    :: btran(bounds%begp:)
    real(r8)               , intent(in)    :: dayl_factor(1:norig)
    type(surfalb_type)     , intent(in)    :: surfalb_vars
    type(photosyns_type)   , intent(inout) :: photosyns_vars
    ! Phase-specific (sunlit or shaded) arrays:
    real(r8), intent(in)    :: par_z    (bounds%begp:,:) ! par absorbed per unit lai for canopy layer (w/m**2)
    real(r8), intent(in)    :: lai_z    (bounds%begp:,:) ! leaf area index for canopy layer
    real(r8), intent(in)    :: vcmaxcint(bounds%begp:)   ! leaf to canopy scaling coefficient
    real(r8), intent(inout) :: alphapsn (bounds%begp:)   ! 13C fractionation factor for PSN ()
    real(r8), intent(inout) :: ci_z     (bounds%begp:,:) ! intracellular leaf CO2 (Pa)
    real(r8), intent(inout) :: rs       (bounds%begp:)   ! leaf stomatal resistance (s/m)
    real(r8), intent(inout) :: rs_z     (bounds%begp:,:) ! canopy layer: leaf stomatal resistance (s/m)
    real(r8), intent(inout) :: lmr      (bounds%begp:)   ! leaf maintenance respiration rate (umol CO2/m**2/s)
    real(r8), intent(inout) :: lmr_z    (bounds%begp:,:) ! canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
    real(r8), intent(inout) :: psn      (bounds%begp:)   ! foliage photosynthesis (umol co2 /m**2/ s) [always +]
    real(r8), intent(inout) :: psn_z    (bounds%begp:,:) ! canopy layer: foliage photosynthesis (umol co2 /m**2/ s) [always +]
    real(r8), intent(inout) :: psn_wc   (bounds%begp:)   ! Rubisco-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]
    real(r8), intent(inout) :: psn_wj   (bounds%begp:)   ! RuBP-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]
    real(r8), intent(inout) :: psn_wp   (bounds%begp:)   ! product-limited foliage photosynthesis (umol co2 /m**2/ s) [always +]
    !
    ! !LOCAL VARIABLES:
    ! Work arrays carried between the steps (indexed by active-filter position f)
    real(r8) :: jmax_z(1:fn,nlevcan)   ! maximum electron transport rate (umol electrons/m**2/s)
    real(r8) :: kn(1:fn)               ! leaf nitrogen decay coefficient
    real(r8) :: psn_wc_z(1:fn,nlevcan) ! Rubisco-limited contribution to psn_z (umol CO2/m**2/s)
    real(r8) :: psn_wj_z(1:fn,nlevcan) ! RuBP-limited contribution to psn_z (umol CO2/m**2/s)
    real(r8) :: psn_wp_z(1:fn,nlevcan) ! product-limited contribution to psn_z (umol CO2/m**2/s)

    ! Step 1: kinetic parameters
    real(r8) :: kc25           ! Michaelis-Menten constant for CO2 at 25C (Pa)
    real(r8) :: ko25           ! Michaelis-Menten constant for O2 at 25C (Pa)
    real(r8) :: cp25           ! CO2 compensation point at 25C (Pa)
    real(r8) :: sco            ! relative specificity of rubisco

    ! Step 2: canopy-top and layer capacities
    real(r8) :: lnc            ! leaf N concentration at the canopy top (gN leaf/m^2)
    real(r8) :: vcmax25top     ! canopy top: maximum rate of carboxylation at 25C (umol CO2/m**2/s)
    real(r8) :: jmax25top      ! canopy top: maximum electron transport rate at 25C (umol electrons/m**2/s)
    real(r8) :: tpu25top       ! canopy top: triose phosphate utilization rate at 25C (umol CO2/m**2/s)
    real(r8) :: lmr25top       ! canopy top: leaf maintenance respiration rate at 25C (umol CO2/m**2/s)
    real(r8) :: kp25top        ! canopy top: initial slope of CO2 response curve (C4 plants) at 25C
    real(r8) :: vcmax25        ! leaf layer: maximum rate of carboxylation at 25C (umol CO2/m**2/s)
    real(r8) :: jmax25         ! leaf layer: maximum electron transport rate at 25C (umol electrons/m**2/s)
    real(r8) :: tpu25          ! leaf layer: triose phosphate utilization rate at 25C (umol CO2/m**2/s)
    real(r8) :: lmr25          ! leaf layer: leaf maintenance respiration rate at 25C (umol CO2/m**2/s)
    real(r8) :: kp25           ! leaf layer: initial slope of CO2 response curve (C4 plants) at 25C
    real(r8) :: vcmaxse        ! entropy term for vcmax (J/mol/K), acclimated to t10
    real(r8) :: jmaxse         ! entropy term for jmax (J/mol/K), acclimated to t10
    real(r8) :: nscaler        ! leaf nitrogen scaling coefficient
    real(r8) :: laican         ! cumulative lai at the middle of the canopy layer
    real(r8) :: total_lai      ! total canopy lai over the radiative layers
    real(r8) :: sum_nscaler    ! sum of the leaf nitrogen profile scaling over the radiative layers

    ! Step 3: leaf-level solution
    type(leaf_env_type) :: leaf ! leaf environment passed to the ci solver
    real(r8) :: cf             ! s m**2/umol -> s/m
    real(r8) :: cs             ! CO2 partial pressure at leaf surface (Pa)
    real(r8) :: gs             ! leaf stomatal conductance (m/s)
    real(r8) :: hs             ! fractional humidity at leaf surface (dimensionless)
    real(r8) :: ciold          ! initial guess for ci
    real(r8) :: gs_mol_err     ! gs_mol for error check
    real(r8) :: qabs           ! PAR absorbed by PS II (umol photons/m**2/s)
    real(r8) :: r1,r2          ! roots of quadratic equation
    real(r8) :: ceair          ! vapor pressure of air, constrained (Pa)
    integer  :: niter          ! solver iteration count (diagnostic only)

    ! Step 4: canopy sums
    real(r8) :: psncan, psncan_wc, psncan_wj, psncan_wp ! canopy sums of psn_z and its limiting contributions
    real(r8) :: lmrcan         ! canopy sum of lmr_z
    real(r8) :: gscan          ! canopy sum of leaf conductance

    integer  :: f,p,t,iv       ! indices
    integer  :: og             ! index into the caller's norig-sized arrays
    integer  :: i_type         ! patch vegetation type
    !------------------------------------------------------------------------------

    associate(                                                  &
         c3psn         => veg_vp%c3psn     , & ! Input:  [real(r8) (:)   ]  photosynthetic pathway: 0. = c4, 1. = c3

         forc_pbot     => top_as%pbot      , & ! Input:  [real(r8) (:)   ]  atmospheric pressure (Pa)

         t_veg         => veg_es%t_veg             , & ! Input:  [real(r8) (:)   ]  vegetation temperature (Kelvin)
         t10           => veg_es%t_a10             , & ! Input:  [real(r8) (:)   ]  10-day running mean of the 2 m temperature (K)
         tgcm          => veg_es%thm               , & ! Input:  [real(r8) (:)   ]  air temperature at agcm reference height (kelvin)

         nrad          => surfalb_vars%nrad_patch   , & ! Input:  [integer  (:)   ]  pft number of canopy layers, above snow for radiative transfer
         tlai_z        => surfalb_vars%tlai_z_patch , & ! Input:  [real(r8) (:,:) ]  pft total leaf area index for canopy layer

         c3flag        => photosyns_vars%c3flag_patch              , & ! Output: [logical  (:)   ]  true if C3 and false if C4
         ac            => photosyns_vars%ac_patch                  , & ! Output: [real(r8) (:,:) ]  Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
         aj            => photosyns_vars%aj_patch                  , & ! Output: [real(r8) (:,:) ]  RuBP-limited gross photosynthesis (umol CO2/m**2/s)
         ap            => photosyns_vars%ap_patch                  , & ! Output: [real(r8) (:,:) ]  product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
         ag            => photosyns_vars%ag_patch                  , & ! Output: [real(r8) (:,:) ]  co-limited gross leaf photosynthesis (umol CO2/m**2/s)
         an            => photosyns_vars%an_patch                  , & ! Output: [real(r8) (:,:) ]  net leaf photosynthesis (umol CO2/m**2/s)
         gb_mol        => photosyns_vars%gb_mol_patch              , & ! Output: [real(r8) (:)   ]  leaf boundary layer conductance (umol H2O/m**2/s)
         gs_mol        => photosyns_vars%gs_mol_patch              , & ! Output: [real(r8) (:,:) ]  leaf stomatal conductance (umol H2O/m**2/s)
         vcmax_z       => photosyns_vars%vcmax_z_patch             , & ! Output: [real(r8) (:,:) ]  maximum rate of carboxylation (umol co2/m**2/s)
         vcmax25_top   => photosyns_vars%vcmax25_top_patch         , & ! Output: [real(r8) (:)   ]  maximum rate of carboxylation at top canopy at 25oC (umol co2/m**2/s)
         cp            => photosyns_vars%cp_patch                  , & ! Output: [real(r8) (:)   ]  CO2 compensation point (Pa)
         kc            => photosyns_vars%kc_patch                  , & ! Output: [real(r8) (:)   ]  Michaelis-Menten constant for CO2 (Pa)
         ko            => photosyns_vars%ko_patch                  , & ! Output: [real(r8) (:)   ]  Michaelis-Menten constant for O2 (Pa)
         qe            => photosyns_vars%qe_patch                  , & ! Output: [real(r8) (:)   ]  quantum efficiency, used only for C4 (mol CO2 / mol photons)
         tpu_z         => photosyns_vars%tpu_z_patch               , & ! Output: [real(r8) (:,:) ]  triose phosphate utilization rate (umol CO2/m**2/s)
         kp_z          => photosyns_vars%kp_z_patch                , & ! Output: [real(r8) (:,:) ]  initial slope of CO2 response curve (C4 plants)
         theta_cj      => photosyns_vars%theta_cj_patch            , & ! Output: [real(r8) (:)   ]  empirical curvature parameter for ac, aj photosynthesis co-limitation
         bbb           => photosyns_vars%bbb_patch                 , & ! Output: [real(r8) (:)   ]  Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
         mbb           => photosyns_vars%mbb_patch                 , & ! Output: [real(r8) (:)   ]  Ball-Berry slope of conductance-photosynthesis relationship
         rh_leaf       => photosyns_vars%rh_leaf_patch               & ! Output: [real(r8) (:)   ]  fractional humidity at leaf surface (dimensionless)
         )

      !$acc enter data create(jmax_z(:,:), kn(:), psn_wc_z(:,:), psn_wj_z(:,:), psn_wp_z(:,:))

      !==============================================================================!
      ! Step 1: patch kinetic parameters at leaf temperature
      !==============================================================================!
      !$acc parallel loop independent gang vector default(present) private(p,og,t,i_type,kc25,ko25,sco,cp25)
      do f = 1, fn
         p = filterp(f)
         og = filter_map(f)
         t = veg_pp%topounit(p)
         i_type = veg_pp%itype(p)
         ! C3 or C4 photosynthesis logical variable
         if (nint(c3psn(i_type)) == 1) then
            c3flag(p) = .true.
         else if (nint(c3psn(i_type)) == 0) then
            c3flag(p) = .false.
         end if
         ! Soil water stress applied to Ball-Berry parameters
         theta_cj(p) = veg_vp%theta_cj(i_type)
         qe(p)  = veg_vp%qe(i_type)
         bbb(p) = max (veg_vp%bbbopt(i_type) * btran(p), 1._r8)
         mbb(p) = veg_vp%mbbopt(i_type)
         ! kc, ko, cp, from: Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
         !
         !       kc25 = 404.9 umol/mol
         !       ko25 = 278.4 mmol/mol
         !       cp25 = 42.75 umol/mol
         !
         ! Derive sco from cp and O2 using present-day O2 (0.209 mol/mol) and re-calculate
         ! cp to account for variation in O2 using cp = 0.5 O2 / sco
         kc25 = (404.9_r8 / 1.e06_r8) * forc_pbot(t)
         ko25 = (278.4_r8 / 1.e03_r8) * forc_pbot(t)
         sco  = 0.5_r8 * 0.209_r8 / (42.75_r8 / 1.e06_r8)
         cp25 = 0.5_r8 * oair(og) / sco

         kc(p) = kc25 * ft(t_veg(p), veg_vp%kcha(i_type))
         ko(p) = ko25 * ft(t_veg(p), veg_vp%koha(i_type))
         cp(p) = cp25 * ft(t_veg(p), veg_vp%cpha(i_type))
      end do

      !==============================================================================!
      ! Step 2: canopy-layer capacities scaled by the leaf nitrogen profile
      ! (cumulative lai at the midpoint of each layer). Respiration is computed
      ! every timestep; the other capacities only for daytime layers.
      !==============================================================================!
      !$acc parallel loop independent gang vector default(present) &
      !$acc private(p,og,i_type,iv,lnc,total_lai,sum_nscaler,vcmax25top,jmax25top,tpu25top,kp25top,lmr25top,&
      !$acc laican,nscaler,lmr25,vcmax25,jmax25,tpu25,kp25,vcmaxse,jmaxse) present(&
      !$acc btran(:),kp_z(:,:),lmr_z(:,:),par_z(:,:),nrad(:),vcmaxcint(:),&
      !$acc vcmax_z(:,:),t10(:),tpu_z(:,:),tlai_z(:,:),alphapsn(:),c3flag(:))
      do f = 1, fn
         p = filterp(f)
         og = filter_map(f)
         i_type = veg_pp%itype(p)

         ! Leaf N profile summed over the radiative layers; needed only when leaf N
         ! (and P) pools set the canopy-top capacity.
         ! NOTE: kn(f) is read here before it is set below for this call; this
         ! preserves the original behavior.
         total_lai   = 0.0_r8
         sum_nscaler = 0.0_r8
         if (nu_com_leaf_physiology .and. .not. (carbon_only .or. carbonphosphorus_only)) then
            laican = 0.0_r8
            !$acc loop seq
            do iv = 1, nrad(p)
               if (iv == 1) then
                  laican = 0.5_r8 * tlai_z(p,iv)
               else
                  laican = laican + 0.5_r8 * (tlai_z(p,iv-1)+tlai_z(p,iv))
               end if
               total_lai = total_lai + tlai_z(p,iv)
               ! If multi-layer code, use explicit profile. If sun/shade big leaf code,
               ! use canopy integrated factor.
               if (nlevcan == 1) then
                  nscaler = 1.0_r8
               else if (nlevcan > 1) then
                  nscaler = exp(-kn(f) * laican)
               end if
               sum_nscaler = sum_nscaler + nscaler
            end do
         end if

         call canopy_top_capacity(p, i_type, dayl_factor(og), total_lai, sum_nscaler, &
              lnc, vcmax25top, jmax25top)

         tpu25top  = 0.167_r8 * vcmax25top
         kp25top   = 20000._r8 * vcmax25top

         ! Nitrogen scaling factor. Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 used
         ! kn = 0.11. Here, derive kn from vcmax25 as in Lloyd et al (2010) Biogeosciences, 7, 1833-1859
         ! Remove daylength factor from vcmax25 so that kn is based on maximum vcmax25
         ! But not used as defined here if using sun/shade big leaf code. Instead,
         ! will use canopy integrated scaling factors from SurfaceAlbedo.
         if (dayl_factor(og) .eq. 0._r8) then
            kn(f) =  0._r8
         else
            kn(f) = exp(0.00963_r8 * vcmax25top/dayl_factor(og) - 2.43_r8)
         end if

         if (use_cn) then
            ! Leaf maintenance respiration to match the base rate used in CN
            ! but with the new temperature functions for C3 and C4 plants.
            !
            ! Base rate for maintenance respiration is from:
            ! M. Ryan, 1991. Effects of climate change on plant respiration.
            ! Ecological Applications, 1(2), 157-167.
            ! Original expression is br = 0.0106 molC/(molN h)
            ! Conversion by molecular weights of C and N gives 2.525e-6 gC/(gN s)
            !
            ! Base rate is at 20C. Adjust to 25C using the CN Q10 = 1.5
            !
            ! CN respiration has units:  g C / g N [leaf] / s. This needs to be
            ! converted from g C / g N [leaf] / s to umol CO2 / m**2 [leaf] / s
            !
            ! Then scale this value at the top of the canopy for canopy depth
            lmr25top = 2.525e-6_r8 * (ParamsShareInst%Q10_mr ** ((25._r8 - 20._r8)/10._r8))
            lmr25top = lmr25top * lnc / 12.e-06_r8
         else
            ! Leaf maintenance respiration in proportion to vcmax25top
            if (c3flag(p)) then
               lmr25top = vcmax25top * 0.015_r8
            else
               lmr25top = vcmax25top * 0.025_r8
            end if
         end if

         ! Loop through canopy layers (above snow).
         laican = 0._r8
         !$acc loop seq
         do iv = 1, nrad(p)
            ! Cumulative lai at middle of layer
            if (iv == 1) then
               laican = 0.5_r8 * tlai_z(p,iv)
            else
               laican = laican + 0.5_r8 * (tlai_z(p,iv-1)+tlai_z(p,iv))
            end if
            ! Scale for leaf nitrogen profile. If multi-layer code, use explicit
            ! profile. If sun/shade big leaf code, use canopy integrated factor.
            if (nlevcan == 1) then
               nscaler = vcmaxcint(p)
               if (nu_com_leaf_physiology) nscaler = 1
            else if (nlevcan > 1) then
               nscaler = exp(-kn(f) * laican)
            end if

            ! Maintenance respiration
            lmr25 = lmr25top * nscaler
            if (c3flag(p)) then
               lmr_z(p,iv) = temperature_scaled(lmr25, t_veg(p), veg_vp%lmrha(i_type), &
                    veg_vp%lmrhd(i_type), veg_vp%lmrse(i_type))
            else
               lmr_z(p,iv) = lmr25 * 2._r8**((t_veg(p)-(tfrz+25._r8))/10._r8)
               lmr_z(p,iv) = lmr_z(p,iv) / (1._r8 + exp( 1.3_r8*(t_veg(p)-(tfrz+55._r8)) ))
            end if

            if (par_z(p,iv) <= 0._r8) then ! night time
               vcmax_z(p,iv) = 0._r8
               jmax_z(f,iv)  = 0._r8
               tpu_z(p,iv)   = 0._r8
               kp_z(p,iv)    = 0._r8
               if ( use_c13 ) then
                  alphapsn(p) = 1._r8
               end if
            else                           ! day time
               vcmax25 = vcmax25top * nscaler
               jmax25 = jmax25top * nscaler
               tpu25 = tpu25top * nscaler
               kp25 = kp25top * nscaler

               ! Adjust for temperature. Activation energies from Bernacchi et al (2001, 2003),
               ! TPU from Harley et al (1992); high temperature deactivation from Leuning (2002)
               ! with entropy terms acclimated to the 10-day mean temperature.
               vcmaxse = 668.39_r8 - 1.07_r8 * min(max((t10(p)-tfrz),11._r8),35._r8)
               jmaxse  = 659.70_r8 - 0.75_r8 * min(max((t10(p)-tfrz),11._r8),35._r8)
               vcmax_z(p,iv) = temperature_scaled(vcmax25, t_veg(p), veg_vp%vcmaxha(i_type), veg_vp%vcmaxhd(i_type), vcmaxse)
               jmax_z(f,iv)  = temperature_scaled(jmax25,  t_veg(p), veg_vp%jmaxha(i_type),  veg_vp%jmaxhd(i_type),  jmaxse)
               tpu_z(p,iv)   = temperature_scaled(tpu25,   t_veg(p), veg_vp%tpuha(i_type),   veg_vp%tpuhd(i_type),   vcmaxse)

               if (.not. c3flag(p)) then
                  vcmax_z(p,iv) = vcmax25 * 2._r8**((t_veg(p)-(tfrz+25._r8))/10._r8)
                  vcmax_z(p,iv) = vcmax_z(p,iv) / (1._r8 + exp( 0.2_r8*((tfrz+15._r8)-t_veg(p)) ))
                  vcmax_z(p,iv) = vcmax_z(p,iv) / (1._r8 + exp( 0.3_r8*(t_veg(p)-(tfrz+40._r8)) ))
               end if

               kp_z(p,iv) = kp25 * 2._r8**((t_veg(p)-(tfrz+25._r8))/10._r8)
            end if
            ! Adjust for soil water
            vcmax_z(p,iv) = vcmax_z(p,iv) * btran(p)
            lmr_z(p,iv) = lmr_z(p,iv) * btran(p)

            ! output variable
            vcmax25_top(p) = vcmax25top
         end do       ! canopy layer loop
      end do     ! patch loop

      !==============================================================================!
      ! Step 3: leaf-level photosynthesis and stomatal conductance
      !==============================================================================!
      !$acc parallel loop independent gang vector default(present) &
      !$acc private(p,og,iv,t,cf,ceair,qabs,r1,r2,niter,ciold,cs,gs,hs,gs_mol_err,leaf) present(&
      !$acc c3flag(:),nrad(:),bbb(:),rh_leaf(:),ap(:,:),ac(:,:),ag(:,:),psn_z(:,:),&
      !$acc rs_z(:,:),an(:,:),aj(:,:),ci_z(:,:))
      do f = 1, fn
         p = filterp(f)
         og = filter_map(f)
         t = veg_pp%topounit(p)

         !$acc loop seq
         do iv = 1, nrad(p)
            ! Leaf boundary layer conductance, umol/m**2/s
            cf = forc_pbot(t)/(rgas*1.e-3_r8*tgcm(p))*1.e06_r8
            gb_mol(p) = 1._r8/rb(og) * cf

            if (par_z(p,iv) <= 0._r8) then           ! night time
               ac(p,iv) = 0._r8
               aj(p,iv) = 0._r8
               ap(p,iv) = 0._r8
               ag(p,iv) = 0._r8
               an(p,iv) = ag(p,iv) - lmr_z(p,iv)
               psn_z(p,iv) = 0._r8
               psn_wc_z(f,iv) = 0._r8
               psn_wj_z(f,iv) = 0._r8
               psn_wp_z(f,iv) = 0._r8
               rs_z(p,iv) = min(rsmax0, 1._r8/bbb(p) * cf)
               ci_z(p,iv) = 0._r8
               rh_leaf(p) = 0._r8

            else                                     ! day time
               ceair = min( eair(og),  esat_tv(og) )

               leaf%pbot   = forc_pbot(t)
               leaf%gb_mol = gb_mol(p)
               leaf%cair   = cair(og)
               leaf%oair   = oair(og)
               leaf%lmr    = lmr_z(p,iv)
               leaf%par    = par_z(p,iv)
               leaf%rh_can = ceair / esat_tv(og)

               ! Electron transport rate for C3 plants.
               qabs  = 0.5_r8 * (1._r8 - fnps) * par_z(p,iv) * par_to_photons
               call quadratic (theta_psii, -(qabs + jmax_z(f,iv)), qabs * jmax_z(f,iv), r1, r2)
               leaf%je = min(r1,r2)

               ! Initial guess for ci, then solve for ci and stomatal conductance
               if (c3flag(p)) then
                  ciold = 0.7_r8 * cair(og)
               else
                  ciold = 0.4_r8 * cair(og)
               end if
               call hybrid(ciold, p, iv, leaf, gs_mol(p,iv), niter, photosyns_vars)

               ! Check for an < 0, in which case gs_mol = bbb
               if (an(p,iv) < 0._r8) gs_mol(p,iv) = bbb(p)

               ! Final estimates for cs and ci (needed for early exit of ci iteration when an < 0)
               cs = cair(og) - 1.4_r8/gb_mol(p) * an(p,iv) * forc_pbot(t)
               cs = max(cs,1.e-06_r8)
               ci_z(p,iv) = cair(og) - an(p,iv) * forc_pbot(t) * (1.4_r8*gs_mol(p,iv)+1.6_r8*gb_mol(p)) / (gb_mol(p)*gs_mol(p,iv))

               ! Convert gs_mol (umol H2O/m**2/s) to gs (m/s) and then to rs (s/m)
               gs = gs_mol(p,iv) / cf
               rs_z(p,iv) = min(1._r8/gs, rsmax0)

               ! Photosynthesis. Save rate-limiting photosynthesis
               psn_z(p,iv) = ag(p,iv)
               psn_wc_z(f,iv) = 0._r8
               psn_wj_z(f,iv) = 0._r8
               psn_wp_z(f,iv) = 0._r8
               if (ac(p,iv) <= aj(p,iv) .and. ac(p,iv) <= ap(p,iv)) then
                  psn_wc_z(f,iv) =  psn_z(p,iv)
               else if (aj(p,iv) < ac(p,iv) .and. aj(p,iv) <= ap(p,iv)) then
                  psn_wj_z(f,iv) =  psn_z(p,iv)
               else if (ap(p,iv) < ac(p,iv) .and. ap(p,iv) < aj(p,iv)) then
                  psn_wp_z(f,iv) =  psn_z(p,iv)
               end if

               ! Make sure iterative solution is correct
               if (gs_mol(p,iv) < 0._r8) then
                  write(iulog,*) 'Negative stomatal conductance:'
                  write (iulog,*)'p,iv,gs_mol= ',p,iv,gs_mol(p,iv)
                  call endrun(decomp_index=p, elmlevel=namep, msg=errmsg(__FILE__, __LINE__))
               end if

               ! Compare with Ball-Berry model: gs_mol = m * an * hs/cs p + b
               hs = (gb_mol(p)*ceair + gs_mol(p,iv)*esat_tv(og)) / ((gb_mol(p)+gs_mol(p,iv))*esat_tv(og))
               rh_leaf(p) = hs
               gs_mol_err = mbb(p)*max(an(p,iv), 0._r8)*hs/cs*forc_pbot(t) + bbb(p)
               if (abs(gs_mol(p,iv)-gs_mol_err) > 1.e-01_r8) then
                  print *, 'Ball-Berry error check - stomatal conductance error:'
                  write (iulog,*) gs_mol(p,iv), gs_mol_err
               end if
            end if    ! night or day if branch
         end do       ! canopy layer loop
      end do          ! patch loop

      !==============================================================================!
      ! Step 4: canopy photosynthesis and stomatal conductance
      ! Sum canopy layer fluxes and then derive effective leaf-level fluxes (per
      ! unit leaf area), which are used in other parts of the model. Here, laican
      ! sums to either laisun or laisha.
      !==============================================================================!
      !$acc parallel loop gang worker independent default(present) &
      !$acc private(p,og,iv,psncan,psncan_wc,psncan_wj,psncan_wp,lmrcan,gscan,laican) present(&
      !$acc lai_z(:,:),nrad(:),rs_z(:,:),rs(:),psn_z(:,:),psn_wp(:),psn_wc(:),&
      !$acc psn_wj(:),lmr(:),lmr_z(:,:),psn(:))
      do f = 1, fn
         p = filterp(f)
         og = filter_map(f)
         psncan = 0._r8
         psncan_wc = 0._r8
         psncan_wj = 0._r8
         psncan_wp = 0._r8
         lmrcan = 0._r8
         gscan = 0._r8
         laican = 0._r8
         !$acc loop seq
         do iv = 1, nrad(p)
            psncan = psncan + psn_z(p,iv) * lai_z(p,iv)
            psncan_wc = psncan_wc + psn_wc_z(f,iv) * lai_z(p,iv)
            psncan_wj = psncan_wj + psn_wj_z(f,iv) * lai_z(p,iv)
            psncan_wp = psncan_wp + psn_wp_z(f,iv) * lai_z(p,iv)
            lmrcan = lmrcan + lmr_z(p,iv) * lai_z(p,iv)
            gscan = gscan + lai_z(p,iv) / (rb(og)+rs_z(p,iv))
            laican = laican + lai_z(p,iv)
         end do
         if (laican > 0._r8) then
            psn(p) = psncan / laican
            psn_wc(p) = psncan_wc / laican
            psn_wj(p) = psncan_wj / laican
            psn_wp(p) = psncan_wp / laican
            lmr(p) = lmrcan / laican
            rs(p) = laican / gscan - rb(og)
         else
            psn(p) =  0._r8
            psn_wc(p) =  0._r8
            psn_wj(p) =  0._r8
            psn_wp(p) =  0._r8
            lmr(p) = 0._r8
            rs(p) = 0._r8
         end if
      end do

      !$acc exit data delete(jmax_z(:,:), kn(:), psn_wc_z(:,:), psn_wj_z(:,:), psn_wp_z(:,:))

    end associate

  end subroutine PhaseLeafPhotosynthesis

  !------------------------------------------------------------------------------
  subroutine canopy_top_capacity(p, i_type, dayl_factor, total_lai, sum_nscaler, &
       lnc, vcmax25top, jmax25top)
    !$acc routine seq
    !
    ! !DESCRIPTION:
    ! Leaf nitrogen concentration and vcmax25/jmax25 at the canopy top, either from
    ! the PFT leaf C:N (default) or from the leaf N (and P) pools when
    ! nu_com_leaf_physiology.
    !
    ! !ARGUMENTS:
    integer , intent(in)  :: p            ! patch index
    integer , intent(in)  :: i_type       ! patch vegetation type
    real(r8), intent(in)  :: dayl_factor  ! scalar (0-1) for daylength
    real(r8), intent(in)  :: total_lai    ! total canopy lai over the radiative layers
    real(r8), intent(in)  :: sum_nscaler  ! sum of the leaf nitrogen profile scaling over the radiative layers
    real(r8), intent(out) :: lnc          ! leaf N concentration at the canopy top (gN leaf/m^2)
    real(r8), intent(out) :: vcmax25top   ! canopy top: maximum rate of carboxylation at 25C (umol CO2/m**2/s)
    real(r8), intent(out) :: jmax25top    ! canopy top: maximum electron transport rate at 25C (umol electrons/m**2/s)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: fnr          ! (gRubisco/gN in Rubisco)
    real(r8) :: act25        ! Rubisco activity at 25 C (umol/gRubisco/s)
    real(r8) :: lpc          ! leaf P concentration (gP leaf/m^2)
    !------------------------------------------------------------------------------

    associate(                                &
         slatop => veg_vp%slatop            , & ! Input:  specific leaf area at top of canopy, projected area basis [m^2/gC]
         leafcn => veg_vp%leafcn            , & ! Input:  leaf C:N (gC/gN)
         flnr   => veg_vp%flnr              , & ! Input:  fraction of leaf N in the Rubisco enzyme (gN Rubisco / gN leaf)
         fnitr  => veg_vp%fnitr             , & ! Input:  foliage nitrogen limitation factor (-)
         t10    => veg_es%t_a10             , & ! Input:  10-day running mean of the 2 m temperature (K)
         leafn  => veg_ns%leafn             , & ! Input:  leaf N (gN/m2)
         leafp  => veg_ps%leafp               & ! Input:  leaf P (gP/m2)
         )

      ! vcmax25 parameters, from CN
      fnr   = veg_vp%fnr(i_type)
      ! Convert rubisco activity units from umol/mgRubisco/min -> umol/gRubisco/s
      act25 = veg_vp%act25(i_type) * 1000.0_r8 / 60.0_r8

      if ( .not. nu_com_leaf_physiology) then
         ! Leaf nitrogen concentration at the top of the canopy (g N leaf / m**2 leaf)
         lnc = 1._r8 / (slatop(i_type) * leafcn(i_type))

         ! vcmax25 at canopy top, as in CN but using lnc at top of the canopy
         vcmax25top = lnc * flnr(i_type) * fnr * act25 * dayl_factor
         if (.not. use_cn) then
            vcmax25top = vcmax25top * fnitr(i_type)
         else if ( carbon_only ) then
            vcmax25top = vcmax25top * fnitr(i_type)
         end if
         ! Parameters derived from vcmax25top. Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
         ! used jmax25 = 1.97 vcmax25, from Wullschleger (1993) Journal of Experimental Botany 44:907-920.
         jmax25top = (2.59_r8 - 0.035_r8*min(max((t10(p)-tfrz),11._r8),35._r8)) * vcmax25top

      else if ( carbon_only  .or.  carbonphosphorus_only ) then
         ! leaf level nutrient control on photosynthesis rate added by Q. Zhu Aug 2015
         lnc = 1._r8 / (slatop(i_type) * leafcn(i_type))
         vcmax25top = lnc * flnr(i_type) * fnr * act25 * dayl_factor
         vcmax25top = vcmax25top * fnitr(i_type)
         jmax25top = (2.59_r8 - 0.035_r8*min(max((t10(p)-tfrz),11._r8),35._r8)) * vcmax25top

      else if ( carbonnitrogen_only ) then
         ! only N control, from Kattge 2009 Global Change Biology 15 (4), 976-991
         if (total_lai > 0.0_r8 .and. sum_nscaler > 0.0_r8) then
            ! dividing by LAI to convert total leaf nitrogen from m2 ground to m2 leaf;
            ! dividing by sum_nscaler to convert total leaf N to leaf N at canopy top
            lnc = leafn(p) / (total_lai * sum_nscaler)
            lnc = min(max(lnc,0.25_r8),3.0_r8) ! based on doi: 10.1002/ece3.1173
         else
            lnc = 0.0_r8
         end if

         vcmax25top = (veg_vp%i_vc(i_type) + veg_vp%s_vc(i_type) * lnc) * dayl_factor
         jmax25top = (2.59_r8 - 0.035_r8*min(max((t10(p)-tfrz),11._r8),35._r8)) * vcmax25top
         vcmax25top = min(max(vcmax25top, 10.0_r8), 150.0_r8)
         jmax25top = min(max(jmax25top, 10.0_r8), 250.0_r8)

      else if (veg_pp%active(p) .and. (i_type .ne. noveg)) then
         ! vcmax25, jmax25 derived from leafn, leafp concentration
         ! Anthony Walker 2014 DOI: 10.1002/ece3.1173
         if (total_lai > 0.0_r8 .and. sum_nscaler > 0.0_r8) then
            lnc = leafn(p) / (total_lai * sum_nscaler)
            lpc = leafp(p) / (total_lai * sum_nscaler)
            lnc = min(max(lnc,0.25_r8),3.0_r8)  ! based on doi: 10.1002/ece3.1173
            lpc = min(max(lpc,0.014_r8),0.85_r8) ! based on doi: 10.1002/ece3.1173
            vcmax25top = exp(vcmax_np1(i_type) + vcmax_np2(i_type)*log(lnc) + &
                 vcmax_np3(i_type)*log(lpc) + vcmax_np4(i_type)*log(lnc)*log(lpc ))&
                 * dayl_factor
            jmax25top = exp(jmax_np1 + jmax_np2*log(vcmax25top) + jmax_np3*log(lpc )) * dayl_factor
            vcmax25top = min(max(vcmax25top, 10.0_r8), 150.0_r8)
            jmax25top = min(max(jmax25top, 10.0_r8), 250.0_r8)
         else
            lnc = 0.0_r8
            vcmax25top = 0.0_r8
            jmax25top = 0.0_r8
         end if
      else
         lnc        = 0.0_r8
         vcmax25top = 0.0_r8
         jmax25top  = 0.0_r8
      end if

    end associate

  end subroutine canopy_top_capacity


  !------------------------------------------------------------------------------
  subroutine PhotosynthesisTotal (fn, filterp, &
       cnstate_vars, canopystate_vars, photosyns_vars)

    ! Note: This subroutine is not called via FATES (RGK)
    !
    ! Determine total photosynthesis
    !
    ! !ARGUMENTS:
    integer                , intent(in)    :: fn                             ! size of pft filter
    integer                , intent(in)    :: filterp(fn)                    ! patch filter
    type(cnstate_type)     , intent(in)    :: cnstate_vars
    type(canopystate_type) , intent(in)    :: canopystate_vars
    type(photosyns_type)   , intent(inout) :: photosyns_vars
    !
    ! !LOCAL VARIABLES:
    integer :: f,p,t,g                    ! indices
    !-----------------------------------------------------------------------

    associate(                                             &
         forc_pco2   => top_as%pco2bot                   , & ! Input:  [real(r8) (:) ]  partial pressure co2 (Pa)
         forc_pc13o2 => top_as%pc13o2bot                 , & ! Input:  [real(r8) (:) ]  partial pressure c13o2 (Pa)
         forc_po2    => top_as%po2bot                    , & ! Input:  [real(r8) (:) ]  partial pressure o2 (Pa)

         rc14_atm    => cnstate_vars%rc14_atm_patch      , & ! Input : [real(r8) (:) ]  C14O2/C12O2 in atmosphere

         laisun      => canopystate_vars%laisun_patch    , & ! Input:  [real(r8) (:) ]  sunlit leaf area
         laisha      => canopystate_vars%laisha_patch    , & ! Input:  [real(r8) (:) ]  shaded leaf area

         psnsun      => photosyns_vars%psnsun_patch      , & ! Input:  [real(r8) (:) ]  sunlit leaf photosynthesis (umol CO2 /m**2/ s)
         psnsha      => photosyns_vars%psnsha_patch      , & ! Input:  [real(r8) (:) ]  shaded leaf photosynthesis (umol CO2 /m**2/ s)
         rc13_canair => photosyns_vars%rc13_canair_patch , & ! Output: [real(r8) (:) ]  C13O2/C12O2 in canopy air
         rc13_psnsun => photosyns_vars%rc13_psnsun_patch , & ! Output: [real(r8) (:) ]  C13O2/C12O2 in sunlit canopy psn flux
         rc13_psnsha => photosyns_vars%rc13_psnsha_patch , & ! Output: [real(r8) (:) ]  C13O2/C12O2 in shaded canopy psn flux
         alphapsnsun => photosyns_vars%alphapsnsun_patch , & ! Output: [real(r8) (:) ]  fractionation factor in sunlit canopy psn flux
         alphapsnsha => photosyns_vars%alphapsnsha_patch , & ! Output: [real(r8) (:) ]  fractionation factor in shaded canopy psn flux
         psnsun_wc   => photosyns_vars%psnsun_wc_patch   , & ! Output: [real(r8) (:) ]  Rubsico-limited sunlit leaf photosynthesis (umol CO2 /m**2/ s)
         psnsun_wj   => photosyns_vars%psnsun_wj_patch   , & ! Output: [real(r8) (:) ]  RuBP-limited sunlit leaf photosynthesis (umol CO2 /m**2/ s)
         psnsun_wp   => photosyns_vars%psnsun_wp_patch   , & ! Output: [real(r8) (:) ]  product-limited sunlit leaf photosynthesis (umol CO2 /m**2/ s)
         psnsha_wc   => photosyns_vars%psnsha_wc_patch   , & ! Output: [real(r8) (:) ]  Rubsico-limited shaded leaf photosynthesis (umol CO2 /m**2/ s)
         psnsha_wj   => photosyns_vars%psnsha_wj_patch   , & ! Output: [real(r8) (:) ]  RuBP-limited shaded leaf photosynthesis (umol CO2 /m**2/ s)
         psnsha_wp   => photosyns_vars%psnsha_wp_patch   , & ! Output: [real(r8) (:) ]  product-limited shaded leaf photosynthesis (umol CO2 /m**2/ s)
         c13_psnsun  => photosyns_vars%c13_psnsun_patch  , & ! Output: [real(r8) (:) ]  sunlit leaf photosynthesis (umol 13CO2 /m**2/ s)
         c13_psnsha  => photosyns_vars%c13_psnsha_patch  , & ! Output: [real(r8) (:) ]  shaded leaf photosynthesis (umol 13CO2 /m**2/ s)
         c14_psnsun  => photosyns_vars%c14_psnsun_patch  , & ! Output: [real(r8) (:) ]  sunlit leaf photosynthesis (umol 14CO2 /m**2/ s)
         c14_psnsha  => photosyns_vars%c14_psnsha_patch  , & ! Output: [real(r8) (:) ]  shaded leaf photosynthesis (umol 14CO2 /m**2/ s)
         fpsn        => photosyns_vars%fpsn_patch        , & ! Output: [real(r8) (:) ]  photosynthesis (umol CO2 /m**2 /s)
         fpsn_wc     => photosyns_vars%fpsn_wc_patch     , & ! Output: [real(r8) (:) ]  Rubisco-limited photosynthesis (umol CO2 /m**2 /s)
         fpsn_wj     => photosyns_vars%fpsn_wj_patch     , & ! Output: [real(r8) (:) ]  RuBP-limited photosynthesis (umol CO2 /m**2 /s)
         fpsn_wp     => photosyns_vars%fpsn_wp_patch       & ! Output: [real(r8) (:) ]  product-limited photosynthesis (umol CO2 /m**2 /s)
         )

      !$acc parallel loop independent gang vector default(present)
      do f = 1, fn
         p = filterp(f)
         g = veg_pp%gridcell(p)
         t = veg_pp%topounit(p)

         if (.not.use_fates) then
            fpsn(p)    = psnsun(p)   *laisun(p) + psnsha(p)   *laisha(p)
            fpsn_wc(p) = psnsun_wc(p)*laisun(p) + psnsha_wc(p)*laisha(p)
            fpsn_wj(p) = psnsun_wj(p)*laisun(p) + psnsha_wj(p)*laisha(p)
            fpsn_wp(p) = psnsun_wp(p)*laisun(p) + psnsha_wp(p)*laisha(p)
         end if

         if (use_cn) then
            if ( use_c13 ) then
               rc13_canair(p) = forc_pc13o2(t)/(forc_pco2(t) - forc_pc13o2(t))
               rc13_psnsun(p) = rc13_canair(p)/alphapsnsun(p)
               rc13_psnsha(p) = rc13_canair(p)/alphapsnsha(p)
               c13_psnsun(p)  = psnsun(p) * (rc13_psnsun(p)/(1._r8+rc13_psnsun(p)))
               c13_psnsha(p)  = psnsha(p) * (rc13_psnsha(p)/(1._r8+rc13_psnsha(p)))

               ! use fixed c13 ratio with del13C of -25 to test the overall c13 structure
               ! c13_psnsun(p) = 0.01095627 * psnsun(p)
               ! c13_psnsha(p) = 0.01095627 * psnsha(p)
            endif
            if ( use_c14 ) then
               c14_psnsun(p) = rc14_atm(p) * psnsun(p)
               c14_psnsha(p) = rc14_atm(p) * psnsha(p)
            endif
         end if

      end do

    end associate

  end subroutine PhotosynthesisTotal

  !------------------------------------------------------------------------------
  subroutine Fractionation(bounds, fn, filterp, &
       cnstate_vars, solarabs_vars, surfalb_vars, photosyns_vars, &
       phase)
    !
    ! !DESCRIPTION:
    ! C13 fractionation during photosynthesis is calculated here after the nitrogen
    ! limitation is taken into account in the CNAllocation module.
    !
    ! !ARGUMENTS:
    type(bounds_type)     , intent(in   )    :: bounds
    integer               , intent(in   )    :: fn                   ! size of pft filter
    integer               , intent(in   )    :: filterp(fn)          ! patch filter
    type(cnstate_type)    , intent(in)    :: cnstate_vars
    type(solarabs_type)   , intent(in)    :: solarabs_vars
    type(surfalb_type)    , intent(in)    :: surfalb_vars
    type(photosyns_type)  , intent(inout)    :: photosyns_vars
    integer      , intent(in   )    :: phase                ! sunlit_phase or shaded_phase
    !
    ! !LOCAL VARIABLES:
    integer  :: f,p,c,t,g,iv            ! indices
    real(r8) :: co2  ! atmospheric co2 partial pressure (pa)
    real(r8) :: ci
    !------------------------------------------------------------------------------

    associate(                                                  &
         forc_pbot   => top_as%pbot                           , & ! Input:  [real(r8) (:)   ]  atmospheric pressure (Pa)
         forc_pco2   => top_as%pco2bot                        , & ! Input:  [real(r8) (:)   ]  partial pressure co2 (Pa)

         c3psn       => veg_vp%c3psn                          , & ! Input:  [real(r8) (:)   ]  photosynthetic pathway: 0. = c4, 1. = c3

         nrad        => surfalb_vars%nrad_patch               , & ! Input:  [integer  (:)   ]  number of canopy layers, above snow for radiative transfer

         downreg     => cnstate_vars%downreg_patch            , & ! Input:  [real(r8) (:)   ]  fractional reduction in GPP due to N limitation (DIM)

         an          => photosyns_vars%an_patch               , & ! Input:  [real(r8) (:,:) ]  net leaf photosynthesis (umol CO2/m**2/s)
         gb_mol      => photosyns_vars%gb_mol_patch           , & ! Input:  [real(r8) (:)   ]  leaf boundary layer conductance (umol H2O/m**2/s)
         gs_mol      => photosyns_vars%gs_mol_patch           , & ! Input:  [real(r8) (:,:) ]  leaf stomatal conductance (umol H2O/m**2/s)
         par_z_sun    =>    solarabs_vars%parsun_z_patch      , & ! Input :  [real(r8) (:,:)]  par absorbed per unit lai for canopy layer (w/m**2)
         alphapsn_sun =>    photosyns_vars%alphapsnsun_patch  , & ! Output:  [real(r8) (:)]
         par_z_sha    =>    solarabs_vars%parsha_z_patch      , & ! Input :  [real(r8) (:,:)]  par absorbed per unit lai for canopy layer (w/m**2)
         alphapsn_sha =>    photosyns_vars%alphapsnsha_patch   &  ! Output:  [real(r8) (:)]
         )

      if (phase == sunlit_phase) then
         !$acc parallel loop independent gang vector collapse(2) default(present)
         do iv = 1,nrad(p)
            do f = 1, fn
               p = filterp(f)
               c = veg_pp%column(p)
               t = veg_pp%topounit(p)
               g = veg_pp%gridcell(p)

               co2 = forc_pco2(t)
               if (par_z_sun(p,iv) <= 0._r8) then           ! night time
                  alphapsn_sun(p) = 1._r8
               else                                     ! day time
                  ci = co2 - ((an(p,iv) * (1._r8-downreg(p)) ) * &
                        forc_pbot(t) * &
                        (1.4_r8*gs_mol(p,iv)+1.6_r8*gb_mol(p)) / (gb_mol(p)*gs_mol(p,iv)))
                   alphapsn_sun(p) = 1._r8 + (((c3psn(veg_pp%itype(p)) * &
                               (4.4_r8 + (22.6_r8*(ci/co2)))) + &
                               ((1._r8 - c3psn(veg_pp%itype(p))) * 4.4_r8))/1000._r8)
               end if
            end do
         end do
      else if (phase == shaded_phase) then
         !$acc parallel loop independent gang vector collapse(2) default(present)
         do iv = 1,nrad(p)
            do f = 1, fn
               p = filterp(f)
               c = veg_pp%column(p)
               t = veg_pp%topounit(p)
               g = veg_pp%gridcell(p)

               co2 = forc_pco2(t)
               if (par_z_sha(p,iv) <= 0._r8) then           ! night time
                  alphapsn_sha(p) = 1._r8
               else                                     ! day time
                  ci = co2 - ((an(p,iv) * (1._r8-downreg(p)) ) * &
                        forc_pbot(t) * &
                        (1.4_r8*gs_mol(p,iv)+1.6_r8*gb_mol(p)) / (gb_mol(p)*gs_mol(p,iv)))

                   alphapsn_sha(p) = 1._r8 + (((c3psn(veg_pp%itype(p)) * &
                               (4.4_r8 + (22.6_r8*(ci/co2)))) + &
                               ((1._r8 - c3psn(veg_pp%itype(p))) * 4.4_r8))/1000._r8)
               end if
            end do
         end do
      end if

    end associate

  end subroutine Fractionation

  !------------------------------------------------------------------------------
  real(r8) function temperature_scaled(rate25, tl, ha, hd, se) result(ans)
    !$acc routine seq
    !
    ! !DESCRIPTION:
    ! Scale a 25C rate to leaf temperature: Arrhenius activation (ft) times high
    ! temperature deactivation (fth), normalized to 1 at 25C.
    !   Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
    !   Bernacchi et al (2003) Plant, Cell and Environment 26:1419-1430
    !   Harley et al (1992) Plant, Cell and Environment 15:271-282 (TPU)
    !   Leuning (2002) Plant, Cell and Environment 25:1205-1210 (deactivation)
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: rate25 ! rate at 25C
    real(r8), intent(in) :: tl     ! leaf temperature (K)
    real(r8), intent(in) :: ha     ! activation energy (J/mol)
    real(r8), intent(in) :: hd     ! deactivation energy (J/mol)
    real(r8), intent(in) :: se     ! entropy term (J/mol/K)
    !------------------------------------------------------------------------------

    ans = rate25 * ft(tl, ha) * fth(tl, hd, se, fth25(hd, se))

  end function temperature_scaled

  !-------------------------------------------------------------------------------
  subroutine hybrid(x0, p, iv, env, gs_mol, iter, photosyns_vars)
    !
    ! DESCRIPTION:
    ! use a hybrid solver to find the root of equation
    ! f(x) = x- h(x),
    !we want to find x, s.t. f(x) = 0.
    !the hybrid approach combines the strength of the newton secant approach (find the solution domain)
    !and the bisection approach implemented with the Brent's method to guarrantee convergence.
    !
    ! REVISION HISTORY:
    !Dec 14/2012: created by Jinyun Tang
    !
    ! ARGUMENTS:
    !$acc routine seq
    implicit none
    real(r8), intent(inout) :: x0              !initial guess and final value of the solution
    integer,  intent(in) :: p, iv              ! patch and canopy layer index
    type(leaf_env_type), intent(in) :: env     ! leaf environment
    real(r8), intent(out) :: gs_mol            ! leaf stomatal conductance (umol H2O/m**2/s)
    integer,  intent(out) :: iter              !number of iterations used, for record only
    type(photosyns_type), intent(inout) :: photosyns_vars
    !
    !! LOCAL VARIABLES
    real(r8) :: x1, f0, f1
    real(r8) :: x, dx
    real(r8), parameter :: eps = 1.e-2_r8      !relative accuracy
    real(r8), parameter :: eps1= 1.e-4_r8
    integer,  parameter :: itmax = 40          !maximum number of iterations
    real(r8) :: tol,minx,minf

    call ci_func(x0, f0, p, iv, env, gs_mol, photosyns_vars)

    if(f0 == 0._r8)return

    minx=x0
    minf=f0
    x1 = x0 * 0.99_r8

    call ci_func(x1, f1, p, iv, env, gs_mol, photosyns_vars)

    if(f1==0._r8)then
       x0 = x1
       return
    endif
    if(f1<minf)then
       minx=x1
       minf=f1
    endif

    !first use the secant approach, then use the brent approach as a backup
    iter = 0
    do
       iter = iter + 1
       dx = - f1 * (x1-x0)/(f1-f0)
       x = x1 + dx
       tol = abs(x) * eps
       if(abs(dx)<tol)then
          x0 = x
          exit
       endif
       x0 = x1
       f0 = f1
       x1 = x

       call ci_func(x1, f1, p, iv, env, gs_mol, photosyns_vars)

       if(f1<minf)then
          minx=x1
          minf=f1
       endif
       if(abs(f1)<=eps1)then
          x0 = x1
          exit
       endif

       !if a root zone is found, use the brent method for a robust backup strategy
       if(f1 * f0 < 0._r8)then
          call brent(x, x0, x1, f0, f1, tol, p, iv, env, gs_mol, photosyns_vars)
          x0=x
          exit
       endif
       if(iter>itmax)then
          !in case of failing to converge within itmax iterations
          !stop at the minimum function
          !this happens because of some other issues besides the stomatal conductance calculation
          !and it happens usually in very dry places and more likely with c4 plants.
          call ci_func(minx, f1, p, iv, env, gs_mol, photosyns_vars)
          exit
       endif
    enddo

  end subroutine hybrid

  !------------------------------------------------------------------------------
  subroutine brent(x, x1, x2, f1, f2, tol, ip, iv, env, gs_mol, photosyns_vars)
    !$acc routine seq
    !
    !!DESCRIPTION:
    !Use Brent's method to find the root of a single variable function ci_func, which is known to exist between x1 and x2.
    !The found root will be updated until its accuracy is tol.
    !
    !!REVISION HISTORY:
    !Dec 14/2012: Jinyun Tang, modified from numerical recipes in F90 by press et al. 1188-1189
    !
    !!ARGUMENTS:
    real(r8), intent(out) :: x                ! indepedent variable of the single value function ci_func(x)
    real(r8), intent(in) :: x1, x2, f1, f2    ! minimum and maximum of the variable domain to search for the solution ci_func(x1) = f1, ci_func(x2)=f2
    real(r8), intent(in) :: tol               ! the error tolerance
    integer,  intent(in) :: ip, iv            ! patch and canopy layer index
    type(leaf_env_type), intent(in) :: env    ! leaf environment
    real(r8), intent(out) :: gs_mol           ! leaf stomatal conductance (umol H2O/m**2/s)
    type(photosyns_type), intent(inout) :: photosyns_vars
    !
    !!LOCAL VARIABLES:
    integer, parameter :: ITMAX=20            !maximum number of iterations
    real(r8), parameter :: EPS=1.e-2_r8       !relative error tolerance
    integer :: iter
    real(r8)  :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
    !------------------------------------------------------------------------------

    a=x1
    b=x2
    fa=f1
    fb=f2
    if((fa > 0._r8 .and. fb > 0._r8).or.(fa < 0._r8 .and. fb < 0._r8))then
       write(iulog,*) 'root must be bracketed for brent'
       call endrun(msg=errmsg(__FILE__, __LINE__))
    endif
    c=b
    fc=fb
    iter = 0
    do
       if(iter==ITMAX)exit
       iter=iter+1
       if((fb > 0._r8 .and. fc > 0._r8) .or. (fb < 0._r8 .and. fc < 0._r8))then
          c=a   !Rename a, b, c and adjust bounding interval d.
          fc=fa
          d=b-a
          e=d
       endif
       if( abs(fc) < abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
       endif
       tol1=2._r8*EPS*abs(b)+0.5_r8*tol  !Convergence check.
       xm=0.5_r8*(c-b)
       if(abs(xm) <= tol1 .or. fb == 0.)then
          x=b
          return
       endif
       if(abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
          s=fb/fa !Attempt inverse quadratic interpolation.
          if(a == c) then
             p=2._r8*xm*s
             q=1._r8-s
          else
             q=fa/fc
             r=fb/fc
             p=s*(2._r8*xm*q*(q-r)-(b-a)*(r-1._r8))
             q=(q-1._r8)*(r-1._r8)*(s-1._r8)
          endif
          if(p > 0._r8) q=-q !Check whether in bounds.
          p=abs(p)
          if(2._r8*p < min(3._r8*xm*q-abs(tol1*q),abs(e*q))) then
             e=d !Accept interpolation.
             d=p/q
          else
             d=xm  !Interpolation failed, use bisection.
             e=d
          endif
       else !Bounds decreasing too slowly, use bisection.
          d=xm
          e=d
       endif
       a=b !Move last best guess to a.
       fa=fb
       if(abs(d) > tol1) then !Evaluate new trial root.
          b=b+d
       else
          b=b+sign(tol1,xm)
       endif

       call ci_func(b, fb, ip, iv, env, gs_mol, photosyns_vars)

       if(fb==0._r8)exit

    enddo

    !if(iter==ITMAX)write(iulog,*) 'brent exceeding maximum iterations', b, fb
    x=b

    return
  end subroutine brent

  !------------------------------------------------------------------------------
  subroutine ci_func(ci, fval, p, iv, env, gs_mol, photosyns_vars)
    !$acc routine seq
    ! DESCRIPTION:
    ! evaluate the function
    ! f(ci)=ci - (ca - (1.37rb+1.65rs))*patm*an
    !
    ! remark:  I am attempting to maintain the original code structure, also
    ! considering one may be interested to output relevant variables for the
    ! photosynthesis model, I have decided to add these relevant variables to
    ! the relevant data types.
    !
    ! ARGUMENTS:
    real(r8)             , intent(in)    :: ci       ! intracellular leaf CO2 (Pa)
    integer              , intent(in)    :: p,iv     ! patch and canopy layer index
    type(leaf_env_type)  , intent(in)    :: env      ! leaf environment
    real(r8)             , intent(out)   :: fval     ! return function of the value f(ci)
    real(r8)             , intent(out)   :: gs_mol   ! leaf stomatal conductance (umol H2O/m**2/s)
    type(photosyns_type) , intent(inout) :: photosyns_vars
    !
    !local variables
    real(r8) :: ai                  ! intermediate co-limited photosynthesis (umol CO2/m**2/s)
    real(r8) :: cs                  ! CO2 partial pressure at leaf surface (Pa)
    real(r8) :: r1, r2              ! roots of quadratic equation
    !------------------------------------------------------------------------------

    associate(&
         c3flag     => photosyns_vars%c3flag_patch             , & ! Input:  [logical  (:)   ]  true if C3 and false if C4
         ac         => photosyns_vars%ac_patch                 , & ! Output: [real(r8) (:,:) ]  Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
         aj         => photosyns_vars%aj_patch                 , & ! Output: [real(r8) (:,:) ]  RuBP-limited gross photosynthesis (umol CO2/m**2/s)
         ap         => photosyns_vars%ap_patch                 , & ! Output: [real(r8) (:,:) ]  product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
         ag         => photosyns_vars%ag_patch                 , & ! Output: [real(r8) (:,:) ]  co-limited gross leaf photosynthesis (umol CO2/m**2/s)
         an         => photosyns_vars%an_patch                 , & ! Output: [real(r8) (:,:) ]  net leaf photosynthesis (umol CO2/m**2/s)
         vcmax_z    => photosyns_vars%vcmax_z_patch            , & ! Input:  [real(r8) (:,:) ]  maximum rate of carboxylation (umol co2/m**2/s)
         cp         => photosyns_vars%cp_patch                 , & ! Input:  [real(r8) (:)   ]  CO2 compensation point (Pa)
         kc         => photosyns_vars%kc_patch                 , & ! Input:  [real(r8) (:)   ]  Michaelis-Menten constant for CO2 (Pa)
         ko         => photosyns_vars%ko_patch                 , & ! Input:  [real(r8) (:)   ]  Michaelis-Menten constant for O2 (Pa)
         qe         => photosyns_vars%qe_patch                 , & ! Input:  [real(r8) (:)   ]  quantum efficiency, used only for C4 (mol CO2 / mol photons)
         tpu_z      => photosyns_vars%tpu_z_patch              , & ! Input:  [real(r8) (:,:) ]  triose phosphate utilization rate (umol CO2/m**2/s)
         kp_z       => photosyns_vars%kp_z_patch               , & ! Input:  [real(r8) (:,:) ]  initial slope of CO2 response curve (C4 plants)
         theta_cj   => photosyns_vars%theta_cj_patch           , & ! Input:  [real(r8) (:)   ]  empirical curvature parameter for ac, aj photosynthesis co-limitation
         bbb        => photosyns_vars%bbb_patch                , & ! Input:  [real(r8) (:)   ]  Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
         mbb        => photosyns_vars%mbb_patch                  & ! Input:  [real(r8) (:)   ]  Ball-Berry slope of conductance-photosynthesis relationship
         )

      if (c3flag(p)) then
         ! C3: Rubisco-limited photosynthesis
         ac(p,iv) = vcmax_z(p,iv) * max(ci-cp(p), 0._r8) / (ci+kc(p)*(1._r8+env%oair/ko(p)))

         ! C3: RuBP-limited photosynthesis
         aj(p,iv) = env%je * max(ci-cp(p), 0._r8) / (4._r8*ci+8._r8*cp(p))

         ! C3: Product-limited photosynthesis
         ap(p,iv) = 3._r8 * tpu_z(p,iv)
      else
         ! C4: Rubisco-limited photosynthesis
         ac(p,iv) = vcmax_z(p,iv)

         ! C4: RuBP-limited photosynthesis
         aj(p,iv) = qe(p) * env%par * par_to_photons

         ! C4: PEP carboxylase-limited (CO2-limited)
         ap(p,iv) = kp_z(p,iv) * max(ci, 0._r8) / env%pbot
      end if

      ! Gross photosynthesis. First co-limit ac and aj. Then co-limit ap
      call quadratic (theta_cj(p), -(ac(p,iv) + aj(p,iv)), ac(p,iv) * aj(p,iv), r1, r2)
      ai = min(r1,r2)

      call quadratic (theta_ip, -(ai + ap(p,iv)), ai * ap(p,iv), r1, r2)
      ag(p,iv) = min(r1,r2)

      ! Net photosynthesis. Exit iteration if an < 0
      an(p,iv) = ag(p,iv) - env%lmr
      if (an(p,iv) < 0._r8) then
         fval = 0._r8
         return
      endif

      ! Quadratic gs_mol calculation with an known. Valid for an >= 0.
      ! With an <= 0, then gs_mol = bbb
      cs = env%cair - 1.4_r8/env%gb_mol * an(p,iv) * env%pbot
      cs = max(cs,1.e-06_r8)
      call quadratic (cs, &
           cs*(env%gb_mol - bbb(p)) - mbb(p)*an(p,iv)*env%pbot, &
           -env%gb_mol*(cs*bbb(p) + mbb(p)*an(p,iv)*env%pbot*env%rh_can), r1, r2)
      gs_mol = max(r1,r2)

      ! Derive new estimate for ci
      fval = ci - env%cair + an(p,iv) * env%pbot * (1.4_r8*gs_mol+1.6_r8*env%gb_mol) / (env%gb_mol*gs_mol)

    end associate

  end subroutine ci_func

  !-------------------------------------------------------------------------------
  function ft(tl, ha) result(ans)
    !$acc routine seq
    !!DESCRIPTION:
    ! photosynthesis temperature response
    !
    ! !REVISION HISTORY
    ! Jinyun Tang separated it out from Photosynthesis, Feb. 07/2013
    !
    !!USES
    use elm_varcon  , only : rgas, tfrz
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: tl  ! leaf temperature in photosynthesis temperature function (K)
    real(r8), intent(in) :: ha  ! activation energy in photosynthesis temperature function (J/mol)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ans
    !-------------------------------------------------------------------------------

    ans = exp( ha / (rgas*1.e-3_r8*(tfrz+25._r8)) * (1._r8 - (tfrz+25._r8)/tl) )

    return
  end function ft

  !-------------------------------------------------------------------------------
  function fth(tl,hd,se,scaleFactor) result(ans)
    !$acc routine seq
    !!DESCRIPTION:
    !photosynthesis temperature inhibition
    !
    ! !REVISION HISTORY
    ! Jinyun Tang separated it out from Photosynthesis, Feb. 07/2013
    !
    use elm_varcon  , only : rgas, tfrz
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: tl  ! leaf temperature in photosynthesis temperature function (K)
    real(r8), intent(in) :: hd  ! deactivation energy in photosynthesis temperature function (J/mol)
    real(r8), intent(in) :: se  ! entropy term in photosynthesis temperature function (J/mol/K)
    real(r8), intent(in) :: scaleFactor  ! scaling factor for high temperature inhibition (25 C = 1.0)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ans
    !-------------------------------------------------------------------------------

    ans = scaleFactor / ( 1._r8 + exp( (-hd+se*tl) / (rgas*1.e-3_r8*tl) ) )

    return
  end function fth

  !-------------------------------------------------------------------------------
  function fth25(hd,se)result(ans)
    !$acc routine seq
    !!DESCRIPTION:
    ! scaling factor for photosynthesis temperature inhibition
    !
    ! !REVISION HISTORY:
    ! Jinyun Tang separated it out from Photosynthesis, Feb. 07/2013
    !
    !!USES
    use elm_varcon  , only : rgas, tfrz
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: hd    ! deactivation energy in photosynthesis temperature function (J/mol)
    real(r8), intent(in) :: se    ! entropy term in photosynthesis temperature function (J/mol/K)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ans
    !-------------------------------------------------------------------------------

    ans = 1._r8 + exp( (-hd+se*(tfrz+25._r8)) / (rgas*1.e-3_r8*(tfrz+25._r8)) )

    return
  end function fth25

end module PhotosynthesisMod
