module VegetationPropertiesType

  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use abortutils     , only : endrun
  use elm_varpar     , only : nlevdecomp
  use elm_varpar     , only : nsoilorder
  use elm_varctl     , only : nu_com
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
     integer , allocatable :: noveg         (:)   ! value for not vegetated
     integer , allocatable :: tree          (:)   ! tree or not?
     real(r8), allocatable :: smpso         (:)   ! soil water potential at full stomatal opening (mm)
     real(r8), allocatable :: smpsc         (:)   ! soil water potential at full stomatal closure (mm)
     real(r8), allocatable :: fnitr         (:)   ! foliage nitrogen limitation factor (-)
     real(r8), allocatable :: foln          (:)   ! foliage nitrogen (%)
     real(r8), allocatable :: dleaf         (:)   ! characteristic leaf dimension (m)
     real(r8), allocatable :: c3psn         (:)   ! photosynthetic pathway: 0. = c4, 1. = c3
     real(r8), allocatable :: xl            (:)   ! leaf/stem orientation index
     real(r8), allocatable :: rhol          (:,:) ! leaf reflectance: 1=vis, 2=nir   (numrad)
     real(r8), allocatable :: rhos          (:,:) ! stem reflectance: 1=vis, 2=nir   (numrad)
     real(r8), allocatable :: taul          (:,:) ! leaf transmittance: 1=vis, 2=nir (numrad)
     real(r8), allocatable :: taus          (:,:) ! stem transmittance: 1=vis, 2=nir (numrad)
     real(r8), allocatable :: z0mr          (:)   ! ratio of momentum roughness length to canopy top height (-)
     real(r8), allocatable :: displar       (:)   ! ratio of displacement height to canopy top height (-)
     real(r8), allocatable :: roota_par     (:)   ! rooting distribution parameter [1/m]
     real(r8), allocatable :: rootb_par     (:)   ! rooting distribution parameter [1/m]
     real(r8), allocatable :: rootprof_beta (:)   ! rooting distribution parameter for C and N inputs [unitless]
     real(r8), allocatable :: dwood         (:)   ! wood density (gC/m3)
     real(r8), allocatable :: slatop        (:)   ! specific leaf area at top of canopy, projected area basis [m^2/gC]
! for plant hydraulics
     real(r8), allocatable :: root_radius   (:)   ! root radius (m)
     real(r8), allocatable :: root_density  (:)   ! root density (gC/m3)
!
     real(r8), allocatable :: dsladlai      (:)   ! dSLA/dLAI, projected area basis [m^2/gC]
     real(r8), allocatable :: leafcn        (:)   ! leaf C:N (gC/gN)
     real(r8), allocatable :: flnr          (:)   ! fraction of leaf N in the Rubisco enzyme (gN Rubisco / gN leaf)
     real(r8), allocatable :: woody         (:)   ! binary flag for woody lifeform (1=woody, 0=not woody)
     real(r8), allocatable :: lflitcn       (:)   ! leaf litter C:N (gC/gN)
     real(r8), allocatable :: frootcn       (:)   ! fine root C:N (gC/gN)
     real(r8), allocatable :: livewdcn      (:)   ! live wood (phloem and ray parenchyma) C:N (gC/gN)
     real(r8), allocatable :: deadwdcn      (:)   ! dead wood (xylem and heartwood) C:N (gC/gN)
     real(r8), allocatable :: graincn       (:)   ! grain C:N (gC/gN) for prognostic crop model
     real(r8), allocatable :: froot_leaf    (:)   ! allocation parameter: new fine root C per new leaf C (gC/gC)
     real(r8), allocatable :: stem_leaf     (:)   ! allocation parameter: new stem c per new leaf C (gC/gC)
     real(r8), allocatable :: croot_stem    (:)   ! allocation parameter: new coarse root C per new stem C (gC/gC)
     real(r8), allocatable :: flivewd       (:)   ! allocation parameter: fraction of new wood that is live 
                                                  ! (phloem and ray parenchyma) (no units)
     real(r8), allocatable :: fcur          (:)   ! allocation parameter: fraction of allocation that goes 
                                                  ! to currently displayed growth, remainder to storage
     real(r8), allocatable :: lf_flab       (:)   ! leaf litter labile fraction
     real(r8), allocatable :: lf_fcel       (:)   ! leaf litter cellulose fraction
     real(r8), allocatable :: lf_flig       (:)   ! leaf litter lignin fraction
     real(r8), allocatable :: fr_flab       (:)   ! fine root litter labile fraction
     real(r8), allocatable :: fr_fcel       (:)   ! fine root litter cellulose fraction
     real(r8), allocatable :: fr_flig       (:)   ! fine root litter lignin fraction
     real(r8), allocatable :: leaf_long     (:)   ! leaf longevity (yrs)
     real(r8), allocatable :: froot_long    (:)   ! fine root longevity (yrs)
     real(r8), allocatable :: evergreen     (:)   ! binary flag for evergreen leaf habit (0 or 1)
     real(r8), allocatable :: stress_decid  (:)   ! binary flag for stress-deciduous leaf habit (0 or 1)
     real(r8), allocatable :: season_decid  (:)   ! binary flag for seasonal-deciduous leaf habit (0 or 1)
     real(r8), allocatable :: cc_leaf       (:)   ! combustion completeness factor for leaf (0 to 1)
     real(r8), allocatable :: cc_lstem      (:)   ! combustion completeness factor for live stem (0 to 1)
     real(r8), allocatable :: cc_dstem      (:)   ! combustion completeness factor for dead stem (0 to 1)
     real(r8), allocatable :: cc_other      (:)   ! combustion completeness factor for other plant tissues (0 to 1)
     real(r8), allocatable :: fm_leaf       (:)   ! fire-related mortality factor for leaf (0 to 1)
     real(r8), allocatable :: fm_lstem      (:)   ! fire-related mortality factor for live stem (0 to 1)
     real(r8), allocatable :: fm_dstem      (:)   ! fire-related mortality factor for dead stem (0 to 1)
     real(r8), allocatable :: fm_other      (:)   ! fire-related mortality factor for other plant tissues (0 to 1)
     real(r8), allocatable :: fm_root       (:)   ! fire-related mortality factor for fine roots (0 to 1)
     real(r8), allocatable :: fm_lroot      (:)   ! fire-related mortality factor for live roots (0 to 1)
     real(r8), allocatable :: fm_droot      (:)   ! fire-related mortality factor for dead roots (0 to 1)
     real(r8), allocatable :: fertnitro     (:)   ! fertilizer applied (crop)
     real(r8), allocatable :: fleafcn       (:)   ! C:N during grain fill; leaf (crop)
     real(r8), allocatable :: ffrootcn      (:)   ! C:N during grain fill; froot (crop)
     real(r8), allocatable :: fstemcn       (:)   ! C:N during grain fill; stem (crop)
     real(r8), allocatable :: presharv      (:)   ! porportion of residue harvested (crop)
     real(r8), allocatable :: convfact      (:)   ! converstion factor to bu/acre (crop)
     real(r8), allocatable :: fyield        (:)   ! fraction of grain that is actually harvested (crop)

     real(r8), allocatable :: leafcp        (:)   ! leaf C:P (gC/gP)
     real(r8), allocatable :: lflitcp       (:)   ! leaf litter C:P (gC/gP)
     real(r8), allocatable :: frootcp       (:)   ! fine root C:P (gC/gP)
     real(r8), allocatable :: livewdcp      (:)   ! live wood (phloem and ray parenchyma) C:P (gC/gP)
     real(r8), allocatable :: deadwdcp      (:)   ! dead wood (xylem and heartwood) C:P (gC/gP)
     real(r8), allocatable :: graincp       (:)   ! grain C:P (gC/gP) for prognostic crop model
     
     ! pft dependent parameters for phosphorus for nutrient competition
     real(r8), allocatable :: vmax_plant_nh4(:)        ! vmax for plant nh4 uptake
     real(r8), allocatable :: vmax_plant_no3(:)        ! vmax for plant no3 uptake
     real(r8), allocatable :: vmax_plant_p(:)          ! vmax for plant p uptake
     real(r8), allocatable :: vmax_minsurf_p_vr(:,:)   ! vmax for p adsorption
     real(r8), allocatable :: km_plant_nh4(:)          ! km for plant nh4 uptake
     real(r8), allocatable :: km_plant_no3(:)          ! km for plant no3 uptake
     real(r8), allocatable :: km_plant_p(:)            ! km for plant p uptake
     real(r8), allocatable :: km_minsurf_p_vr(:,:)     ! km for p adsorption
     real(r8)              :: km_decomp_nh4            ! km for microbial decomposer nh4 uptake
     real(r8)              :: km_decomp_no3            ! km for microbial decomposer no3 uptake
     real(r8)              :: km_decomp_p              ! km for microbial decomposer p uptake
     real(r8)              :: km_nit                   ! km for nitrifier nh4 uptake
     real(r8)              :: km_den                   ! km for denitrifier no3 uptake
     real(r8), allocatable :: decompmicc_patch_vr(:,:) ! microbial decomposer biomass gc/m3
     real(r8), allocatable :: vmax_nfix(:)             ! vmax of symbiotic n2 fixation
     real(r8), allocatable :: km_nfix(:)               ! km of symbiotic n2 fixation
     real(r8), allocatable :: vmax_ptase(:)            ! vmax of biochemical p production
     real(r8)              :: km_ptase                 ! km of biochemical p production
     real(r8)              :: lamda_ptase              ! critical value that incur biochemical production
     real(r8), allocatable :: i_vc(:)                  ! intercept of photosynthesis vcmax ~ leaf n content regression model
     real(r8), allocatable :: s_vc(:)                  ! slope of photosynthesis vcmax ~ leaf n content regression model

     real(r8), allocatable :: alpha_nfix(:)            ! fraction of fixed N goes directly to plant
     real(r8), allocatable :: alpha_ptase(:)           ! fraction of phosphatase produced P goes directly to plant
     real(r8), allocatable :: ccost_nfix(:)            ! plant C cost per unit N produced by N2 fixation
     real(r8), allocatable :: ccost_ptase(:)           ! plant C cost per unit P produced by phosphatase
     real(r8), allocatable :: fnr(:)              !fraction of nitrogen in RuBisCO
     real(r8), allocatable :: act25(:)
     real(r8), allocatable :: kcha(:)             !Activation energy for kc
     real(r8), allocatable :: koha(:)             !Activation energy for ko
     real(r8), allocatable :: cpha(:)             !Activation energy for cp
     real(r8), allocatable :: vcmaxha(:)          !Activation energy for vcmax
     real(r8), allocatable :: jmaxha(:)           !Activation energy for jmax
     real(r8), allocatable :: tpuha(:)            !Activation energy for tpu
     real(r8), allocatable :: lmrha(:)            !Acitivation energy for lmr
     real(r8), allocatable :: vcmaxhd(:)          !Deactivation energy for vcmax
     real(r8), allocatable :: jmaxhd(:)           !Deactivation energy for jmax
     real(r8), allocatable :: tpuhd(:)            !Deactivation energy for tpu
     real(r8), allocatable :: lmrhd(:)            !Deacitivation energy for lmr
     real(r8), allocatable :: lmrse(:)            !SE for lmr
     real(r8), allocatable :: qe(:)               !Quantum efficiency
     real(r8), allocatable :: theta_cj(:)         !
     real(r8), allocatable :: bbbopt(:)           !Ball-Berry stomatal conductance intercept
     real(r8), allocatable :: mbbopt(:)           !Ball-Berry stomatal conductance slope
     real(r8), allocatable :: nstor(:)            !Nitrogen storage pool timescale 
     real(r8), allocatable :: br_xr(:)            !Base rate for excess respiration
     real(r8)              :: tc_stress           !Critial temperature for moisture stress


   contains
   procedure, public :: Init => veg_vp_init
  
   end type vegetation_properties_type

  type(vegetation_properties_type), public :: veg_vp ! patch ecophysiological constants structure

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
    use pftvarcon , only : i_vc, s_vc
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
    allocate(this%smpso         (0:numpft))        ; this%smpso        (:)   =nan
    allocate(this%smpsc         (0:numpft))        ; this%smpsc        (:)   =nan
    allocate(this%fnitr         (0:numpft))        ; this%fnitr        (:)   =nan
    allocate(this%foln          (0:numpft))        ; this%foln         (:)   =nan
    allocate(this%dleaf         (0:numpft))        ; this%dleaf        (:)   =nan
    allocate(this%c3psn         (0:numpft))        ; this%c3psn        (:)   =nan
    allocate(this%xl            (0:numpft))        ; this%xl           (:)   =nan
    allocate(this%rhol          (0:numpft,numrad)) ; this%rhol         (:,:) =nan
    allocate(this%rhos          (0:numpft,numrad)) ; this%rhos         (:,:) =nan
    allocate(this%taul          (0:numpft,numrad)) ; this%taul         (:,:) =nan
    allocate(this%taus          (0:numpft,numrad)) ; this%taus         (:,:) =nan
    allocate(this%z0mr          (0:numpft))        ; this%z0mr         (:)   =nan
    allocate(this%displar       (0:numpft))        ; this%displar      (:)   =nan
    allocate(this%roota_par     (0:numpft))        ; this%roota_par    (:)   =nan
    allocate(this%rootb_par     (0:numpft))        ; this%rootb_par    (:)   =nan
    allocate(this%slatop        (0:numpft))        ; this%slatop       (:)   =nan
    allocate(this%dsladlai      (0:numpft))        ; this%dsladlai     (:)   =nan
    allocate(this%leafcn        (0:numpft))        ; this%leafcn       (:)   =nan
    allocate(this%flnr          (0:numpft))        ; this%flnr         (:)   =nan
    allocate(this%woody         (0:numpft))        ; this%woody        (:)   =nan
    allocate(this%lflitcn       (0:numpft))        ; this%lflitcn      (:)   =nan
    allocate(this%frootcn       (0:numpft))        ; this%frootcn      (:)   =nan
    allocate(this%livewdcn      (0:numpft))        ; this%livewdcn     (:)   =nan
    allocate(this%deadwdcn      (0:numpft))        ; this%deadwdcn     (:)   =nan
    allocate(this%graincn       (0:numpft))        ; this%graincn      (:)   =nan
    allocate(this%froot_leaf    (0:numpft))        ; this%froot_leaf   (:)   =nan
    allocate(this%stem_leaf     (0:numpft))        ; this%stem_leaf    (:)   =nan
    allocate(this%croot_stem    (0:numpft))        ; this%croot_stem   (:)   =nan
    allocate(this%flivewd       (0:numpft))        ; this%flivewd      (:)   =nan
    allocate(this%fcur          (0:numpft))        ; this%fcur         (:)   =nan
    allocate(this%lf_flab       (0:numpft))        ; this%lf_flab      (:)   =nan
    allocate(this%lf_fcel       (0:numpft))        ; this%lf_fcel      (:)   =nan
    allocate(this%lf_flig       (0:numpft))        ; this%lf_flig      (:)   =nan
    allocate(this%fr_flab       (0:numpft))        ; this%fr_flab      (:)   =nan
    allocate(this%fr_fcel       (0:numpft))        ; this%fr_fcel      (:)   =nan
    allocate(this%fr_flig       (0:numpft))        ; this%fr_flig      (:)   =nan
    allocate(this%leaf_long     (0:numpft))        ; this%leaf_long    (:)   =nan
    allocate(this%froot_long    (0:numpft))        ; this%froot_long   (:)   =nan
    allocate(this%evergreen     (0:numpft))        ; this%evergreen    (:)   =nan
    allocate(this%stress_decid  (0:numpft))        ; this%stress_decid (:)   =nan
    allocate(this%season_decid  (0:numpft))        ; this%season_decid (:)   =nan
    allocate(this%dwood         (0:numpft))        ; this%dwood        (:)   =nan
    allocate(this%root_radius   (0:numpft))        ; this%root_radius  (:)   =nan
    allocate(this%root_density  (0:numpft))        ; this%root_density (:)   =nan
    allocate(this%rootprof_beta (0:numpft))        ; this%rootprof_beta(:)   =nan
    allocate(this%fertnitro     (0:numpft))        ; this%fertnitro    (:)   =nan
    allocate(this%fleafcn       (0:numpft))        ; this%fleafcn      (:)   =nan
    allocate(this%ffrootcn      (0:numpft))        ; this%ffrootcn     (:)   =nan
    allocate(this%fstemcn       (0:numpft))        ; this%fstemcn      (:)   =nan
    allocate(this%presharv      (0:numpft))        ; this%presharv     (:)   =nan
    allocate(this%convfact      (0:numpft))        ; this%convfact     (:)   =nan
    allocate(this%fyield        (0:numpft))        ; this%fyield       (:)   =nan


    allocate(this%leafcp        (0:numpft))        ; this%leafcp       (:)   =nan
    allocate(this%lflitcp       (0:numpft))        ; this%lflitcp      (:)   =nan
    allocate(this%frootcp       (0:numpft))        ; this%frootcp      (:)   =nan
    allocate(this%livewdcp      (0:numpft))        ; this%livewdcp     (:)   =nan
    allocate(this%deadwdcp      (0:numpft))        ; this%deadwdcp     (:)   =nan
    allocate(this%graincp       (0:numpft))        ; this%graincp      (:)   =nan
    
    allocate( this%alpha_nfix    (0:numpft))                     ; this%alpha_nfix    (:)        =nan
    allocate( this%alpha_ptase   (0:numpft))                     ; this%alpha_ptase   (:)        =nan
    allocate( this%ccost_nfix    (0:numpft))                     ; this%ccost_nfix    (:)        =nan
    allocate( this%ccost_ptase   (0:numpft))                     ; this%ccost_ptase   (:)        =nan
    allocate( this%vmax_plant_nh4(0:numpft))                     ; this%vmax_plant_nh4(:)        =nan
    allocate( this%vmax_plant_no3(0:numpft))                     ; this%vmax_plant_no3(:)        =nan
    allocate( this%vmax_plant_p(0:numpft))                       ; this%vmax_plant_p(:)          =nan
    allocate( this%vmax_minsurf_p_vr(0:nsoilorder,1:nlevdecomp)) ; this%vmax_minsurf_p_vr(:,:)   =nan
    allocate( this%km_plant_nh4(0:numpft))                       ; this%km_plant_nh4(:)          =nan
    allocate( this%km_plant_no3(0:numpft))                       ; this%km_plant_no3(:)          =nan
    allocate( this%km_plant_p(0:numpft))                         ; this%km_plant_p(:)            =nan
    allocate( this%km_minsurf_p_vr(0:nsoilorder,1:nlevdecomp))   ; this%km_minsurf_p_vr(:,:)     =nan
    allocate( this%decompmicc_patch_vr(0:numpft,1:nlevdecomp))   ; this%decompmicc_patch_vr(:,:) =nan
    allocate( this%vmax_ptase(0:numpft))                         ; this%vmax_ptase(:)            =nan
    allocate( this%i_vc(0:numpft))                               ; this%i_vc(:)                  =nan
    allocate( this%s_vc(0:numpft))                               ; this%s_vc(:)                  =nan
    allocate( this%vmax_nfix(0:numpft))                          ; this%vmax_nfix(:)             =nan
    allocate( this%km_nfix(0:numpft))                            ; this%km_nfix(:)               =nan
    allocate( this%fnr(0:numpft))                                ; this%fnr(:)                   =nan
    allocate( this%act25(0:numpft))                              ; this%act25(:)                 =nan
    allocate( this%kcha(0:numpft))                               ; this%kcha(:)                  =nan
    allocate( this%koha(0:numpft))                               ; this%koha(:)                  =nan
    allocate( this%cpha(0:numpft))                               ; this%cpha(:)                  =nan
    allocate( this%vcmaxha(0:numpft))                            ; this%vcmaxha(:)               =nan
    allocate( this%jmaxha(0:numpft))                             ; this%jmaxha(:)                =nan
    allocate( this%tpuha(0:numpft))                              ; this%tpuha(:)                 =nan
    allocate( this%lmrha(0:numpft))                              ; this%lmrha(:)                 =nan
    allocate( this%vcmaxhd(0:numpft))                            ; this%vcmaxhd(:)               =nan
    allocate( this%jmaxhd(0:numpft))                             ; this%jmaxhd(:)                =nan
    allocate( this%tpuhd(0:numpft))                              ; this%tpuhd(:)                 =nan
    allocate( this%lmrhd(0:numpft))                              ; this%lmrhd(:)                 =nan
    allocate( this%lmrse(0:numpft))                              ; this%lmrse(:)                 =nan
    allocate( this%qe(0:numpft))                                 ; this%qe(:)                    =nan
    allocate( this%theta_cj(0:numpft))                           ; this%theta_cj(:)              =nan
    allocate( this%bbbopt(0:numpft))                             ; this%bbbopt(:)                =nan
    allocate( this%mbbopt(0:numpft))                             ; this%mbbopt(:)                =nan
    allocate( this%nstor(0:numpft))                              ; this%nstor(:)                 =nan
    allocate( this%br_xr(0:numpft))                              ; this%br_xr(:)                 =nan

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
