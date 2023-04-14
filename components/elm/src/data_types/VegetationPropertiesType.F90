module VegetationPropertiesType

  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use elm_varpar     , only : nlevdecomp
  use elm_varpar     , only : nsoilorder
  use elm_varctl     , only : nu_com
  use elm_varcon     , only : spval, ispval
  !
  implicit none
  save
  public
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  !DW    public :: VegetationPropertiesInit
  !
  ! !PUBLIC TYPES:
  type, public :: vegetation_properties_type
     integer , pointer :: noveg         (:)   => null() ! value for not vegetated
     integer , pointer :: tree          (:)   => null() ! tree or not?
     real(r8), pointer :: smpso         (:)   => null() ! soil water potential at full stomatal opening (mm)
     real(r8), pointer :: smpsc         (:)   => null() ! soil water potential at full stomatal closure (mm)
     real(r8), pointer :: fnitr         (:)   => null() ! foliage nitrogen limitation factor (-)
     real(r8), pointer :: foln          (:)   => null() ! foliage nitrogen (%)
     real(r8), pointer :: dleaf         (:)   => null() ! characteristic leaf dimension (m)
     real(r8), pointer :: c3psn         (:)   => null() ! photosynthetic pathway: 0. = c4, 1. = c3
     real(r8), pointer :: xl            (:)   => null() ! leaf/stem orientation index
     real(r8), pointer :: rhol          (:,:) => null() ! leaf reflectance: 1=vis, 2=nir   (numrad)
     real(r8), pointer :: rhos          (:,:) => null() ! stem reflectance: 1=vis, 2=nir   (numrad)
     real(r8), pointer :: taul          (:,:) => null() ! leaf transmittance: 1=vis, 2=nir (numrad)
     real(r8), pointer :: taus          (:,:) => null() ! stem transmittance: 1=vis, 2=nir (numrad)
     real(r8), pointer :: z0mr          (:)   => null() ! ratio of momentum roughness length to canopy top height (-)
     real(r8), pointer :: displar       (:)   => null() ! ratio of displacement height to canopy top height (-)
     real(r8), pointer :: roota_par     (:)   => null() ! rooting distribution parameter [1/m]
     real(r8), pointer :: rootb_par     (:)   => null() ! rooting distribution parameter [1/m]
     real(r8), pointer :: rootprof_beta (:)   => null() ! rooting distribution parameter for C and N inputs [unitless]
     real(r8), pointer :: dwood         (:)   => null() ! wood density (gC/m3)
     real(r8), pointer :: slatop        (:)   => null() ! specific leaf area at top of canopy, projected area basis [m^2/gC]
! for plant hydraulics
     real(r8), pointer :: root_radius   (:) => null()   ! root radius (m)
     real(r8), pointer :: root_density  (:) => null()   ! root density (gC/m3)
!
     real(r8), pointer :: dsladlai      (:) => null()  ! dSLA/dLAI, projected area basis [m^2/gC]
     real(r8), pointer :: leafcn        (:) => null()  ! leaf C:N (gC/gN)
     real(r8), pointer :: flnr          (:) => null()  ! fraction of leaf N in the Rubisco enzyme (gN Rubisco / gN leaf)
     real(r8), pointer :: woody         (:) => null()  ! binary flag for woody lifeform (1=woody, 0=not woody)
     real(r8), pointer :: lflitcn       (:) => null()  ! leaf litter C:N (gC/gN)
     real(r8), pointer :: frootcn       (:) => null()  ! fine root C:N (gC/gN)
     real(r8), pointer :: livewdcn      (:) => null()  ! live wood (phloem and ray parenchyma) C:N (gC/gN)
     real(r8), pointer :: deadwdcn      (:) => null()  ! dead wood (xylem and heartwood) C:N (gC/gN)
     real(r8), pointer :: graincn       (:) => null()  ! grain C:N (gC/gN) for prognostic crop model
     real(r8), pointer :: froot_leaf    (:) => null()  ! allocation parameter: new fine root C per new leaf C (gC/gC)
     real(r8), pointer :: stem_leaf     (:) => null()  ! allocation parameter: new stem c per new leaf C (gC/gC)
     real(r8), pointer :: croot_stem    (:) => null()  ! allocation parameter: new coarse root C per new stem C (gC/gC)
     real(r8), pointer :: flivewd       (:) => null()  ! allocation parameter: fraction of new wood that is live
                                                       ! (phloem and ray parenchyma) (no units)
     real(r8), pointer :: fcur          (:) => null()  ! allocation parameter: fraction of allocation that goes
                                                       ! to currently displayed growth, remainder to storage
     real(r8), pointer :: lf_flab       (:) => null()  ! leaf litter labile fraction
     real(r8), pointer :: lf_fcel       (:) => null()  ! leaf litter cellulose fraction
     real(r8), pointer :: lf_flig       (:) => null()  ! leaf litter lignin fraction
     real(r8), pointer :: fr_flab       (:) => null()  ! fine root litter labile fraction
     real(r8), pointer :: fr_fcel       (:) => null()  ! fine root litter cellulose fraction
     real(r8), pointer :: fr_flig       (:) => null()  ! fine root litter lignin fraction
     real(r8), pointer :: leaf_long     (:) => null()  ! leaf longevity (yrs)
     real(r8), pointer :: froot_long    (:) => null()  ! fine root longevity (yrs)
     real(r8), pointer :: evergreen     (:) => null()  ! binary flag for evergreen leaf habit (0 or 1)
     real(r8), pointer :: stress_decid  (:) => null()  ! binary flag for stress-deciduous leaf habit (0 or 1)
     real(r8), pointer :: season_decid  (:) => null()  ! binary flag for seasonal-deciduous leaf habit (0 or 1)
     real(r8), pointer :: cc_leaf       (:) => null()  ! combustion completeness factor for leaf (0 to 1)
     real(r8), pointer :: cc_lstem      (:) => null()  ! combustion completeness factor for live stem (0 to 1)
     real(r8), pointer :: cc_dstem      (:) => null()  ! combustion completeness factor for dead stem (0 to 1)
     real(r8), pointer :: cc_other      (:) => null()  ! combustion completeness factor for other plant tissues (0 to 1)
     real(r8), pointer :: fm_leaf       (:) => null()  ! fire-related mortality factor for leaf (0 to 1)
     real(r8), pointer :: fm_lstem      (:) => null()  ! fire-related mortality factor for live stem (0 to 1)
     real(r8), pointer :: fm_dstem      (:) => null()  ! fire-related mortality factor for dead stem (0 to 1)
     real(r8), pointer :: fm_other      (:) => null()  ! fire-related mortality factor for other plant tissues (0 to 1)
     real(r8), pointer :: fm_root       (:) => null()  ! fire-related mortality factor for fine roots (0 to 1)
     real(r8), pointer :: fm_lroot      (:) => null()  ! fire-related mortality factor for live roots (0 to 1)
     real(r8), pointer :: fm_droot      (:) => null()  ! fire-related mortality factor for dead roots (0 to 1)
     real(r8), pointer :: fertnitro     (:) => null()  ! fertilizer applied (crop)
     real(r8), pointer :: fleafcn       (:) => null()  ! C:N during grain fill; leaf (crop)
     real(r8), pointer :: ffrootcn      (:) => null()  ! C:N during grain fill; froot (crop)
     real(r8), pointer :: fstemcn       (:) => null()  ! C:N during grain fill; stem (crop)
     real(r8), pointer :: presharv      (:) => null()  ! porportion of residue harvested (crop)
     real(r8), pointer :: convfact      (:) => null()  ! converstion factor to bu/acre (crop)
     real(r8), pointer :: fyield        (:) => null()  ! fraction of grain that is actually harvested (crop)

     real(r8), pointer :: leafcp        (:) => null()  ! leaf C:P (gC/gP)
     real(r8), pointer :: lflitcp       (:) => null()  ! leaf litter C:P (gC/gP)
     real(r8), pointer :: frootcp       (:) => null()  ! fine root C:P (gC/gP)
     real(r8), pointer :: livewdcp      (:) => null()  ! live wood (phloem and ray parenchyma) C:P (gC/gP)
     real(r8), pointer :: deadwdcp      (:) => null()  ! dead wood (xylem and heartwood) C:P (gC/gP)
     real(r8), pointer :: graincp       (:) => null()  ! grain C:P (gC/gP) for prognostic crop model

     ! pft dependent parameters for phosphorus for nutrient competition
     real(r8), pointer :: vmax_plant_nh4(:)      => null()   ! vmax for plant nh4 uptake
     real(r8), pointer :: vmax_plant_no3(:)      => null()   ! vmax for plant no3 uptake
     real(r8), pointer :: vmax_plant_p(:)        => null()   ! vmax for plant p uptake
     real(r8), pointer :: vmax_minsurf_p_vr(:,:) => null()   ! vmax for p adsorption
     real(r8), pointer :: km_plant_nh4(:)        => null()   ! km for plant nh4 uptake
     real(r8), pointer :: km_plant_no3(:)        => null()   ! km for plant no3 uptake
     real(r8), pointer :: km_plant_p(:)          => null()   ! km for plant p uptake
     real(r8), pointer :: km_minsurf_p_vr(:,:)   => null()   ! km for p adsorption
     real(r8), pointer :: km_decomp_nh4          => null()  ! km for microbial decomposer nh4 uptake
     real(r8), pointer :: km_decomp_no3          => null()  ! km for microbial decomposer no3 uptake
     real(r8), pointer :: km_decomp_p            => null()  ! km for microbial decomposer p uptake
     real(r8), pointer :: km_nit                 => null()  ! km for nitrifier nh4 uptake
     real(r8), pointer :: km_den                 => null()  ! km for denitrifier no3 uptake
     real(r8), pointer :: decompmicc_patch_vr(:,:) => null()! microbial decomposer biomass gc/m3
     real(r8), pointer :: vmax_nfix(:)             => null()! vmax of symbiotic n2 fixation
     real(r8), pointer :: km_nfix(:)               => null()! km of symbiotic n2 fixation
     real(r8), pointer :: vmax_ptase(:)            => null()! vmax of biochemical p production
     real(r8), pointer :: km_ptase                 => null()! km of biochemical p production
     real(r8), pointer :: lamda_ptase              => null()! critical value that incur biochemical production
     real(r8), pointer :: i_vc(:)          => null()        ! intercept of photosynthesis vcmax ~ leaf n content regression model
     real(r8), pointer :: s_vc(:)          => null()        ! slope of photosynthesis vcmax ~ leaf n content regression model
     real(r8), pointer :: nsc_rtime(:)     => null()        ! non-structural carbon residence time 
     real(r8), pointer :: pinit_beta1(:)   => null()        ! shaping parameter for P initialization
     real(r8), pointer :: pinit_beta2(:)   => null()        ! shaping parameter for P initialization
     real(r8), pointer :: alpha_nfix(:)    => null()        ! fraction of fixed N goes directly to plant
     real(r8), pointer :: alpha_ptase(:)   => null()        ! fraction of phosphatase produced P goes directly to plant
     real(r8), pointer :: ccost_nfix(:)    => null()        ! plant C cost per unit N produced by N2 fixation
     real(r8), pointer :: ccost_ptase(:)   => null()        ! plant C cost per unit P produced by phosphatase
     real(r8), pointer :: fnr(:)           => null()   !fraction of nitrogen in RuBisCO
     real(r8), pointer :: act25(:)         => null()   !
     real(r8), pointer :: kcha(:)          => null()   !Activation energy for kc
     real(r8), pointer :: koha(:)          => null()   !Activation energy for ko
     real(r8), pointer :: cpha(:)          => null()   !Activation energy for cp
     real(r8), pointer :: vcmaxha(:)       => null()   !Activation energy for vcmax
     real(r8), pointer :: jmaxha(:)        => null()   !Activation energy for jmax
     real(r8), pointer :: tpuha(:)         => null()   !Activation energy for tpu
     real(r8), pointer :: lmrha(:)         => null()   !Acitivation energy for lmr
     real(r8), pointer :: vcmaxhd(:)       => null()   !Deactivation energy for vcmax
     real(r8), pointer :: jmaxhd(:)        => null()   !Deactivation energy for jmax
     real(r8), pointer :: tpuhd(:)         => null()   !Deactivation energy for tpu
     real(r8), pointer :: lmrhd(:)         => null()   !Deacitivation energy for lmr
     real(r8), pointer :: lmrse(:)         => null()   !SE for lmr
     real(r8), pointer :: qe(:)            => null()   !Quantum efficiency
     real(r8), pointer :: theta_cj(:)      => null()   !
     real(r8), pointer :: bbbopt(:)        => null()   !Ball-Berry stomatal conductance intercept
     real(r8), pointer :: mbbopt(:)        => null()   !Ball-Berry stomatal conductance slope
     real(r8), pointer :: nstor(:)         => null()   !Nitrogen storage pool timescale
     real(r8), pointer :: br_xr(:)         => null()   !Base rate for excess respiration
     real(r8), pointer :: tc_stress        => null()   !Critial temperature for moisture stress


   contains
   procedure, public :: Init => veg_vp_init

   end type vegetation_properties_type

  type(vegetation_properties_type), public :: veg_vp ! patch ecophysiological constants structure
  !$acc declare create(veg_vp)
contains

  !-----------------------------------------------------------------------
  subroutine veg_vp_init(this)
    !
    ! !USES:
    use elm_varpar, only : numrad, numpft
    use pftvarcon , only : ntree, smpso, smpsc, fnitr
    use pftvarcon , only : z0mr, displar, dleaf, rhol, rhos, taul, taus, xl
    use pftvarcon , only : c3psn, slatop, dsladlai, leafcn, flnr, woody
    use pftvarcon , only : lflitcn, frootcn, livewdcn, deadwdcn, froot_leaf, stem_leaf, croot_stem
    use pftvarcon , only : flivewd, fcur, lf_flab, lf_fcel, lf_flig, fr_flab, fr_fcel, fr_flig
    use pftvarcon , only : leaf_long, froot_long, evergreen, stress_decid, season_decid
    use pftvarcon , only : fertnitro, graincn, fleafcn, ffrootcn, fstemcn, dwood
    use pftvarcon , only : presharv, convfact, fyield
    use pftvarcon , only : leafcp,lflitcp, frootcp, livewdcp, deadwdcp,graincp
    use pftvarcon , only : vmax_plant_nh4, vmax_plant_no3, vmax_plant_p, vmax_minsurf_p_vr
    use pftvarcon , only : km_plant_nh4, km_plant_no3, km_plant_p, km_minsurf_p_vr
    use pftvarcon , only : km_decomp_nh4, km_decomp_no3, km_decomp_p, km_nit, km_den
    use pftvarcon , only : decompmicc_patch_vr
    use pftvarcon , only : vmax_nfix, km_nfix
    use pftvarcon , only : alpha_nfix, alpha_ptase,ccost_nfix,ccost_ptase
    use pftvarcon , only : vmax_ptase, km_ptase, lamda_ptase
    use pftvarcon , only : i_vc, s_vc, nsc_rtime, pinit_beta1, pinit_beta2
    use pftvarcon , only : leafcn_obs, frootcn_obs, livewdcn_obs, deadwdcn_obs
    use pftvarcon , only : leafcp_obs, frootcp_obs, livewdcp_obs, deadwdcp_obs
    use pftvarcon , only : fnr, act25, kcha, koha, cpha, vcmaxha, jmaxha, tpuha
    use pftvarcon , only : lmrha, vcmaxhd, jmaxhd, tpuhd, lmrse, qe, theta_cj
    use pftvarcon , only : bbbopt, mbbopt, nstor, br_xr, tc_stress, lmrhd
    !

    class (vegetation_properties_type) :: this

    !LOCAL VARIABLES:
    integer :: m, ib, j
    !------------------------------------------------------------------------

    allocate(this%noveg         (0:numpft))        ; this%noveg        (:)   =huge(1)
    allocate(this%tree          (0:numpft))        ; this%tree         (:)   =huge(1)
    allocate(this%smpso         (0:numpft))        ; this%smpso        (:)   =spval
    allocate(this%smpsc         (0:numpft))        ; this%smpsc        (:)   =spval
    allocate(this%fnitr         (0:numpft))        ; this%fnitr        (:)   =spval
    allocate(this%foln          (0:numpft))        ; this%foln         (:)   =spval
    allocate(this%dleaf         (0:numpft))        ; this%dleaf        (:)   =spval
    allocate(this%c3psn         (0:numpft))        ; this%c3psn        (:)   =spval
    allocate(this%xl            (0:numpft))        ; this%xl           (:)   =spval
    allocate(this%rhol          (0:numpft,numrad)) ; this%rhol         (:,:) =spval
    allocate(this%rhos          (0:numpft,numrad)) ; this%rhos         (:,:) =spval
    allocate(this%taul          (0:numpft,numrad)) ; this%taul         (:,:) =spval
    allocate(this%taus          (0:numpft,numrad)) ; this%taus         (:,:) =spval
    allocate(this%z0mr          (0:numpft))        ; this%z0mr         (:)   =spval
    allocate(this%displar       (0:numpft))        ; this%displar      (:)   =spval
    allocate(this%roota_par     (0:numpft))        ; this%roota_par    (:)   =spval
    allocate(this%rootb_par     (0:numpft))        ; this%rootb_par    (:)   =spval
    allocate(this%slatop        (0:numpft))        ; this%slatop       (:)   =spval
    allocate(this%dsladlai      (0:numpft))        ; this%dsladlai     (:)   =spval
    allocate(this%leafcn        (0:numpft))        ; this%leafcn       (:)   =spval
    allocate(this%flnr          (0:numpft))        ; this%flnr         (:)   =spval
    allocate(this%woody         (0:numpft))        ; this%woody        (:)   =spval
    allocate(this%lflitcn       (0:numpft))        ; this%lflitcn      (:)   =spval
    allocate(this%frootcn       (0:numpft))        ; this%frootcn      (:)   =spval
    allocate(this%livewdcn      (0:numpft))        ; this%livewdcn     (:)   =spval
    allocate(this%deadwdcn      (0:numpft))        ; this%deadwdcn     (:)   =spval
    allocate(this%graincn       (0:numpft))        ; this%graincn      (:)   =spval
    allocate(this%froot_leaf    (0:numpft))        ; this%froot_leaf   (:)   =spval
    allocate(this%stem_leaf     (0:numpft))        ; this%stem_leaf    (:)   =spval
    allocate(this%croot_stem    (0:numpft))        ; this%croot_stem   (:)   =spval
    allocate(this%flivewd       (0:numpft))        ; this%flivewd      (:)   =spval
    allocate(this%fcur          (0:numpft))        ; this%fcur         (:)   =spval
    allocate(this%lf_flab       (0:numpft))        ; this%lf_flab      (:)   =spval
    allocate(this%lf_fcel       (0:numpft))        ; this%lf_fcel      (:)   =spval
    allocate(this%lf_flig       (0:numpft))        ; this%lf_flig      (:)   =spval
    allocate(this%fr_flab       (0:numpft))        ; this%fr_flab      (:)   =spval
    allocate(this%fr_fcel       (0:numpft))        ; this%fr_fcel      (:)   =spval
    allocate(this%fr_flig       (0:numpft))        ; this%fr_flig      (:)   =spval
    allocate(this%leaf_long     (0:numpft))        ; this%leaf_long    (:)   =spval
    allocate(this%froot_long    (0:numpft))        ; this%froot_long   (:)   =spval
    allocate(this%evergreen     (0:numpft))        ; this%evergreen    (:)   =spval
    allocate(this%stress_decid  (0:numpft))        ; this%stress_decid (:)   =spval
    allocate(this%season_decid  (0:numpft))        ; this%season_decid (:)   =spval
    allocate(this%dwood         (0:numpft))        ; this%dwood        (:)   =spval
    allocate(this%root_radius   (0:numpft))        ; this%root_radius  (:)   =spval
    allocate(this%root_density  (0:numpft))        ; this%root_density (:)   =spval
    allocate(this%rootprof_beta (0:numpft))        ; this%rootprof_beta(:)   =spval
    allocate(this%fertnitro     (0:numpft))        ; this%fertnitro    (:)   =spval
    allocate(this%fleafcn       (0:numpft))        ; this%fleafcn      (:)   =spval
    allocate(this%ffrootcn      (0:numpft))        ; this%ffrootcn     (:)   =spval
    allocate(this%fstemcn       (0:numpft))        ; this%fstemcn      (:)   =spval
    allocate(this%presharv      (0:numpft))        ; this%presharv     (:)   =spval
    allocate(this%convfact      (0:numpft))        ; this%convfact     (:)   =spval
    allocate(this%fyield        (0:numpft))        ; this%fyield       (:)   =spval


    allocate(this%leafcp        (0:numpft))        ; this%leafcp       (:)   =spval
    allocate(this%lflitcp       (0:numpft))        ; this%lflitcp      (:)   =spval
    allocate(this%frootcp       (0:numpft))        ; this%frootcp      (:)   =spval
    allocate(this%livewdcp      (0:numpft))        ; this%livewdcp     (:)   =spval
    allocate(this%deadwdcp      (0:numpft))        ; this%deadwdcp     (:)   =spval
    allocate(this%graincp       (0:numpft))        ; this%graincp      (:)   =spval

    allocate( this%alpha_nfix    (0:numpft))                     ; this%alpha_nfix    (:)        =spval
    allocate( this%alpha_ptase   (0:numpft))                     ; this%alpha_ptase   (:)        =spval
    allocate( this%ccost_nfix    (0:numpft))                     ; this%ccost_nfix    (:)        =spval
    allocate( this%ccost_ptase   (0:numpft))                     ; this%ccost_ptase   (:)        =spval
    allocate( this%vmax_plant_nh4(0:numpft))                     ; this%vmax_plant_nh4(:)        =spval
    allocate( this%vmax_plant_no3(0:numpft))                     ; this%vmax_plant_no3(:)        =spval
    allocate( this%vmax_plant_p(0:numpft))                       ; this%vmax_plant_p(:)          =spval
    allocate( this%vmax_minsurf_p_vr(0:nsoilorder,1:nlevdecomp)) ; this%vmax_minsurf_p_vr(:,:)   =spval
    allocate( this%km_plant_nh4(0:numpft))                       ; this%km_plant_nh4(:)          =spval
    allocate( this%km_plant_no3(0:numpft))                       ; this%km_plant_no3(:)          =spval
    allocate( this%km_plant_p(0:numpft))                         ; this%km_plant_p(:)            =spval
    allocate( this%km_minsurf_p_vr(0:nsoilorder,1:nlevdecomp))   ; this%km_minsurf_p_vr(:,:)     =spval
    allocate( this%decompmicc_patch_vr(0:numpft,1:nlevdecomp))   ; this%decompmicc_patch_vr(:,:) =spval
    allocate( this%vmax_ptase(0:numpft))                         ; this%vmax_ptase(:)            =spval
    allocate( this%i_vc(0:numpft))                               ; this%i_vc(:)                  =spval
    allocate( this%s_vc(0:numpft))                               ; this%s_vc(:)                  =spval
    allocate( this%nsc_rtime(0:numpft))                          ; this%nsc_rtime(:)             =spval
    allocate( this%pinit_beta1(0:nsoilorder))                    ; this%pinit_beta1(:)           =spval
    allocate( this%pinit_beta2(0:nsoilorder))                    ; this%pinit_beta2(:)           =spval
    allocate( this%vmax_nfix(0:numpft))                          ; this%vmax_nfix(:)             =spval
    allocate( this%km_nfix(0:numpft))                            ; this%km_nfix(:)               =spval
    allocate( this%fnr(0:numpft))                                ; this%fnr(:)                   =spval
    allocate( this%act25(0:numpft))                              ; this%act25(:)                 =spval
    allocate( this%kcha(0:numpft))                               ; this%kcha(:)                  =spval
    allocate( this%koha(0:numpft))                               ; this%koha(:)                  =spval
    allocate( this%cpha(0:numpft))                               ; this%cpha(:)                  =spval
    allocate( this%vcmaxha(0:numpft))                            ; this%vcmaxha(:)               =spval
    allocate( this%jmaxha(0:numpft))                             ; this%jmaxha(:)                =spval
    allocate( this%tpuha(0:numpft))                              ; this%tpuha(:)                 =spval
    allocate( this%lmrha(0:numpft))                              ; this%lmrha(:)                 =spval
    allocate( this%vcmaxhd(0:numpft))                            ; this%vcmaxhd(:)               =spval
    allocate( this%jmaxhd(0:numpft))                             ; this%jmaxhd(:)                =spval
    allocate( this%tpuhd(0:numpft))                              ; this%tpuhd(:)                 =spval
    allocate( this%lmrhd(0:numpft))                              ; this%lmrhd(:)                 =spval
    allocate( this%lmrse(0:numpft))                              ; this%lmrse(:)                 =spval
    allocate( this%qe(0:numpft))                                 ; this%qe(:)                    =spval
    allocate( this%theta_cj(0:numpft))                           ; this%theta_cj(:)              =spval
    allocate( this%bbbopt(0:numpft))                             ; this%bbbopt(:)                =spval
    allocate( this%mbbopt(0:numpft))                             ; this%mbbopt(:)                =spval
    allocate( this%nstor(0:numpft))                              ; this%nstor(:)                 =spval
    allocate( this%br_xr(0:numpft))                              ; this%br_xr(:)                 =spval

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    allocate(this%km_decomp_nh4)
    allocate(this%km_decomp_no3)
    allocate(this%km_decomp_p  )
    allocate(this%km_nit       )
    allocate(this%km_den       )
    allocate(this%km_ptase     )
    allocate(this%lamda_ptase  )
    allocate(this%tc_stress    )
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    do m = 0,numpft

       if (m <= ntree) then
          this%tree(m) = 1
       else
          this%tree(m) = 0
       end if

       do ib = 1,numrad
          this%rhol(m,ib)   = rhol(m,ib)
          this%rhos(m,ib)   = rhos(m,ib)
          this%taul(m,ib)   = taul(m,ib)
          this%taus(m,ib)   = taus(m,ib)
       end do

       this%z0mr(m)         = z0mr(m)
       this%displar(m)      = displar(m)
       this%dleaf(m)        = dleaf(m)
       this%xl(m)           = xl(m)
       this%c3psn(m)        = c3psn(m)
       this%slatop(m)       = slatop(m)
       this%dsladlai(m)     = dsladlai(m)
       this%leafcn(m)       = leafcn(m)
       this%flnr(m)         = flnr(m)
       this%smpso(m)        = smpso(m)
       this%smpsc(m)        = smpsc(m)
       this%fnitr(m)        = fnitr(m)
       this%woody(m)        = woody(m)
       this%lflitcn(m)      = lflitcn(m)
       this%frootcn(m)      = frootcn(m)
       this%livewdcn(m)     = livewdcn(m)
       this%deadwdcn(m)     = deadwdcn(m)
       this%graincn(m)      = graincn(m)
       this%froot_leaf(m)   = froot_leaf(m)
       this%stem_leaf(m)    = stem_leaf(m)
       this%croot_stem(m)   = croot_stem(m)
       this%flivewd(m)      = flivewd(m)
       this%fcur(m)         = fcur(m)
       this%lf_flab(m)      = lf_flab(m)
       this%lf_fcel(m)      = lf_fcel(m)
       this%lf_flig(m)      = lf_flig(m)
       this%fr_flab(m)      = fr_flab(m)
       this%fr_fcel(m)      = fr_fcel(m)
       this%fr_flig(m)      = fr_flig(m)
       this%leaf_long(m)    = leaf_long(m)
       this%froot_long(m)   = froot_long(m)
       this%evergreen(m)    = evergreen(m)
       this%stress_decid(m) = stress_decid(m)
       this%season_decid(m) = season_decid(m)
       this%dwood(m)        = dwood
       this%root_radius(m)  = 0.29e-03_r8 !(m)
       this%root_density(m) = 0.31e06_r8 !(g biomass / m3 root)
       this%fertnitro(m)    = fertnitro(m)
       this%fleafcn(m)      = fleafcn(m)
       this%ffrootcn(m)     = ffrootcn(m)
       this%fstemcn(m)      = fstemcn(m)
       this%presharv(m)     = presharv(m)
       this%convfact(m)     = convfact(m)
       this%fyield(m)       = fyield(m)

       this%leafcp(m)       = leafcp(m)
       this%lflitcp(m)      = lflitcp(m)
       this%frootcp(m)      = frootcp(m)
       this%livewdcp(m)     = livewdcp(m)
       this%deadwdcp(m)     = deadwdcp(m)
       this%graincp(m)      = graincp(m)
       this%fnr(m)          = fnr(m)
       this%act25(m)        = act25(m)
       this%kcha(m)         = kcha(m)
       this%koha(m)         = koha(m)
       this%cpha(m)         = cpha(m)
       this%vcmaxha(m)      = vcmaxha(m)
       this%jmaxha(m)       = jmaxha(m)
       this%tpuha(m)        = tpuha(m)
       this%lmrha(m)        = lmrha(m)
       this%vcmaxhd(m)      = vcmaxhd(m)
       this%jmaxhd(m)       = jmaxhd(m)
       this%tpuhd(m)        = tpuhd(m)
       this%lmrhd(m)        = lmrhd(m)
       this%lmrse(m)        = lmrse(m)
       this%qe(m)           = qe(m)
       this%theta_cj(m)     = theta_cj(m)
       this%bbbopt(m)       = bbbopt(m)
       this%mbbopt(m)       = mbbopt(m)
       this%nstor(m)        = nstor(m)
       this%br_xr(m)        = br_xr(m)

    end do

    do m = 0,numpft
        this%alpha_nfix(m)     = alpha_nfix(m)
        this%alpha_ptase(m)    = alpha_ptase(m)
        this%ccost_nfix(m)     = ccost_nfix(m)
        this%ccost_ptase(m)    = ccost_ptase(m)
        this%vmax_plant_nh4(m) = vmax_plant_nh4(m)
        this%vmax_plant_no3(m) = vmax_plant_no3(m)
        this%vmax_plant_p(m)   = vmax_plant_p(m)
        this%km_plant_nh4(m)   = km_plant_nh4(m)
        this%km_plant_no3(m)   = km_plant_no3(m)
        this%km_plant_p(m)     = km_plant_p(m)
        this%i_vc(m)           = i_vc(m)
        this%s_vc(m)           = s_vc(m)
        this%nsc_rtime(m)      = nsc_rtime(m)
        this%vmax_nfix(m)      = vmax_nfix(m)
        this%km_nfix(m)        = km_nfix(m)
        this%vmax_ptase(m)     = vmax_ptase(m)

        do j = 1 , nlevdecomp
           this%decompmicc_patch_vr(m,j) = decompmicc_patch_vr(j,m)
        end do

        if (nu_com .ne. 'RD') then ! use new stoichiometry for eca and mic competition
           this%leafcn(m)     = leafcn_obs(m)
           this%frootcn(m)    = frootcn_obs(m)
           this%livewdcn(m)   = livewdcn_obs(m)
           this%deadwdcn(m)   = deadwdcn_obs(m)
           this%leafcp(m)     = leafcp_obs(m)
           this%frootcp(m)    = frootcp_obs(m)
           this%livewdcp(m)   = livewdcp_obs(m)
           this%deadwdcp(m)   = deadwdcp_obs(m)
        end if
    end do

    do m = 0, nsoilorder
       this%pinit_beta1(m) = pinit_beta1(m)
       this%pinit_beta2(m) = pinit_beta2(m)
       do j = 1 , nlevdecomp
          this%vmax_minsurf_p_vr(m,j) = vmax_minsurf_p_vr(j,m)
          this%km_minsurf_p_vr(m,j) = km_minsurf_p_vr(j,m)
       end do
    end do

    this%km_decomp_nh4 = km_decomp_nh4
    this%km_decomp_no3 = km_decomp_no3
    this%km_decomp_p   = km_decomp_p
    this%km_nit        = km_nit
    this%km_den        = km_den
    this%km_ptase      = km_ptase
    this%lamda_ptase   = lamda_ptase
    this%tc_stress     = tc_stress

  end subroutine veg_vp_init

end module VegetationPropertiesType
