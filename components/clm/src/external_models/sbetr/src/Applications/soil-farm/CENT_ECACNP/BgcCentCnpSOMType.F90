module BgcCentSOMType
!
!DESCRIPTION
!module defines the century decomposition
  use bshr_kind_mod       , only : r8 => shr_kind_r8
  use bshr_log_mod  , only : errMsg => shr_log_errMsg
implicit none
  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  !totally 7 om pools
  !met, cell, lig, cwd, som1, som2, som3
  !during decomposition, the cnp ratios of som1, som2, and som3 are all fixed
  !the litter pools have their cnp ratios varying with nutrient status
  !We consider som1 and som2 as DOM, som3 as humus pool
  !
  integer :: ncentpools
  type, public :: CentSom_type
    real(r8), pointer :: cn_ratios(:)  => null()
    real(r8), pointer :: cp_ratios(:) => null()
    real(r8), pointer :: cc14_ratios(:)=> null()
    real(r8), pointer :: cc13_ratios(:)=> null()
    real(r8), pointer :: def_cn(:)=> null()
    real(r8), pointer :: def_cp(:)=> null()
    real(r8), pointer :: def_cc13(:)=> null()
    real(r8), pointer :: def_cc14(:)=> null()

    !private parameters
    real(r8) :: rf_l1s1_bgc    !co2 production when metabolic carbon is decomposed
    real(r8) :: rf_l2s1_bgc    !co2 production when cellulose is decomposed
    real(r8) :: rf_l3s2_bgc    !co2 production when lignin is decomposed
    real(r8) :: rf_s2s1_bgc
    real(r8) :: rf_s3s1_bgc
    real(r8) :: cwd_fcel
    real(r8) :: cwd_flig
    real(r8) :: lwd_fcel
    real(r8) :: lwd_flig
    real(r8) :: fwd_fcel
    real(r8) :: fwd_flig
    real(r8) :: lit_flig
    real(r8) :: k_decay_lit1
    real(r8) :: k_decay_lit2
    real(r8) :: k_decay_lit3
    real(r8) :: k_decay_som1
    real(r8) :: k_decay_som2
    real(r8) :: k_decay_som3
    real(r8) :: k_decay_cwd  !coarse root
    real(r8) :: k_decay_lwd  !large wood
    real(r8) :: k_decay_fwd  !fine branch wood

    logical  :: use_c13
    logical  :: use_c14
  contains
    procedure, public :: Init
    procedure, public :: run_decomp
    procedure, public :: calc_pot_min_np_flx
    procedure, private :: calc_cascade_matrix
    procedure, private :: calc_som_decay_k
    procedure, private :: calc_som_decay_r
    procedure, private :: calc_cnp_ratios
    procedure, private :: InitAllocate
    procedure, private :: InitPar
    procedure, private :: calc_potential_aerobic_hr
  end type CentSom_type
contains

  subroutine Init(this, centurybgc_index, biogeo_con, bstatus)

  use BgcCentCnpIndexType , only : centurybgc_index_type
  use BiogeoConType       , only : BiogeoCon_type
  use BetrStatusType      , only : betr_status_type
  implicit none
  class(CentSom_type)         , intent(inout) :: this
  type(centurybgc_index_type) , intent(in)    :: centurybgc_index
  type(BiogeoCon_type)        , intent(in)    :: biogeo_con
  type(betr_status_type)      , intent(out)   :: bstatus

  call bstatus%reset()
  ncentpools = centurybgc_index%nom_pools

  call this%InitAllocate()

  call this%InitPar(centurybgc_index, biogeo_con)
  end subroutine Init
!------------------------------------------
  subroutine InitAllocate (this)
  !
  !
  implicit none
  class(CentSom_type), intent(inout) :: this

  allocate(this%cn_ratios(ncentpools));
  allocate(this%cp_ratios(ncentpools));
  allocate(this%cc14_ratios(ncentpools)); this%cc14_ratios(:) = 0._r8
  allocate(this%cc13_ratios(ncentpools)); this%cc13_ratios(:) = 0._r8
  allocate(this%def_cn(ncentpools));
  allocate(this%def_cp(ncentpools));
  allocate(this%def_cc13(ncentpools));this%def_cc13(:) = 0._r8
  allocate(this%def_cc14(ncentpools));this%def_cc14(:) = 0._r8
  end subroutine InitAllocate
!------------------------------------------
  subroutine InitPar(this,centurybgc_index,  biogeo_con)
  !
  ! intialize model parameters
  use BiogeoConType , only : BiogeoCon_type
  use BgcCentCnpIndexType , only : centurybgc_index_type
  use tracer_varcon    , only : catomw, natomw, patomw
  implicit none
  class(CentSom_type)  , intent(inout) :: this
  type(centurybgc_index_type) , intent(in) :: centurybgc_index
  type(BiogeoCon_type) , intent(in)    :: biogeo_con

  this%rf_l1s1_bgc    = biogeo_con%rf_l1s1_bgc
  this%rf_l2s1_bgc    = biogeo_con%rf_l2s1_bgc
  this%rf_l3s2_bgc    = biogeo_con%rf_l3s2_bgc
  this%rf_s2s1_bgc    = biogeo_con%rf_s2s1_bgc
  this%rf_s3s1_bgc    = biogeo_con%rf_s3s1_bgc
  this%cwd_fcel   = biogeo_con%cwd_fcel_bgc
  this%cwd_flig   = biogeo_con%cwd_flig_bgc
  this%lwd_fcel   = biogeo_con%lwd_fcel_bgc
  this%lwd_flig   = biogeo_con%lwd_flig_bgc
  this%fwd_fcel   = biogeo_con%fwd_fcel_bgc
  this%fwd_flig   = biogeo_con%fwd_flig_bgc

  this%k_decay_lit1   =  biogeo_con%k_decay_lit1
  this%k_decay_lit2   =  biogeo_con%k_decay_lit2
  this%k_decay_lit3   =  biogeo_con%k_decay_lit3
  this%k_decay_som1   =  biogeo_con%k_decay_som1
  this%k_decay_som2   =  biogeo_con%k_decay_som2
  this%k_decay_som3   =  biogeo_con%k_decay_som3
  this%k_decay_cwd    =  biogeo_con%k_decay_cwd
  this%k_decay_lwd    =  biogeo_con%k_decay_lwd
  this%k_decay_fwd    =  biogeo_con%k_decay_fwd


  this%def_cn(centurybgc_index%lit1) = biogeo_con%init_cn_met * natomw/catomw
  this%def_cn(centurybgc_index%lit2) = biogeo_con%init_cn_cel * natomw/catomw
  this%def_cn(centurybgc_index%lit3) = biogeo_con%init_cn_lig * natomw/catomw
  this%def_cn(centurybgc_index%cwd)  = biogeo_con%init_cn_cwd * natomw/catomw
  this%def_cn(centurybgc_index%lwd)  = biogeo_con%init_cn_lwd * natomw/catomw
  this%def_cn(centurybgc_index%fwd)  = biogeo_con%init_cn_fwd * natomw/catomw

  this%def_cn(centurybgc_index%som1) = biogeo_con%init_cn_som1 * natomw/catomw
  this%def_cn(centurybgc_index%som2) = biogeo_con%init_cn_som2 * natomw/catomw
  this%def_cn(centurybgc_index%som3) = biogeo_con%init_cn_som3 * natomw/catomw

  this%def_cp(centurybgc_index%lit1) = biogeo_con%init_cp_met * patomw/catomw
  this%def_cp(centurybgc_index%lit2) = biogeo_con%init_cp_cel * patomw/catomw
  this%def_cp(centurybgc_index%lit3) = biogeo_con%init_cp_lig * patomw/catomw
  this%def_cp(centurybgc_index%cwd)  = biogeo_con%init_cp_cwd * patomw/catomw
  this%def_cp(centurybgc_index%lwd)  = biogeo_con%init_cp_lwd * patomw/catomw
  this%def_cp(centurybgc_index%fwd)  = biogeo_con%init_cp_fwd * patomw/catomw

  this%def_cp(centurybgc_index%som1) = biogeo_con%init_cp_som1 * patomw/catomw
  this%def_cp(centurybgc_index%som2) = biogeo_con%init_cp_som2 * patomw/catomw
  this%def_cp(centurybgc_index%som3) = biogeo_con%init_cp_som3 * patomw/catomw

  this%use_c13=biogeo_con%use_c13
  this%use_c14=biogeo_con%use_c14

  if(this%use_c13)then
    this%def_cc13(centurybgc_index%lit1) = biogeo_con%init_cc13_met
    this%def_cc13(centurybgc_index%lit2) = biogeo_con%init_cc13_cel
    this%def_cc13(centurybgc_index%lit3) = biogeo_con%init_cc13_lig
    this%def_cc13(centurybgc_index%cwd)  = biogeo_con%init_cc13_cwd
    this%def_cc13(centurybgc_index%lwd)  = biogeo_con%init_cc13_lwd
    this%def_cc13(centurybgc_index%fwd)  = biogeo_con%init_cc13_fwd
    this%def_cc13(centurybgc_index%som1) = biogeo_con%init_cc13_som1
    this%def_cc13(centurybgc_index%som2) = biogeo_con%init_cc13_som2
    this%def_cc13(centurybgc_index%som3) = biogeo_con%init_cc13_som3
  endif

  if(this%use_c14)then
    this%def_cc14(centurybgc_index%lit1) = biogeo_con%init_cc14_met
    this%def_cc14(centurybgc_index%lit2) = biogeo_con%init_cc14_cel
    this%def_cc14(centurybgc_index%lit3) = biogeo_con%init_cc14_lig
    this%def_cc14(centurybgc_index%cwd)  = biogeo_con%init_cc14_cwd
    this%def_cc14(centurybgc_index%lwd)  = biogeo_con%init_cc14_lwd
    this%def_cc14(centurybgc_index%fwd)  = biogeo_con%init_cc14_fwd
    this%def_cc14(centurybgc_index%som1) = biogeo_con%init_cc14_som1
    this%def_cc14(centurybgc_index%som2) = biogeo_con%init_cc14_som2
    this%def_cc14(centurybgc_index%som3) = biogeo_con%init_cc14_som3
  endif

  end subroutine InitPar
!------------------------------------------

  subroutine run_decomp(this, is_surf, centurybgc_index, dtime, ystates,&
      decompkf_eca, pct_sand, pct_clay, alpha_n, alpha_p, cascade_matrix, &
      k_decay, pot_co2_hr, bstatus)
  !
  !DESCRIPTION
  !
  use BgcCentCnpIndexType , only : centurybgc_index_type
  use BgcCentCnpDecompType      , only : DecompCent_type
  use BetrStatusType      , only : betr_status_type
  implicit none
  class(CentSom_type)         , intent(inout) :: this
  type(centurybgc_index_type) , intent(in) :: centurybgc_index
  real(r8)                    , intent(in) :: dtime
  real(r8)                    , intent(in) :: ystates(1:centurybgc_index%nom_tot_elms)
  type(DecompCent_type)       , intent(in) :: decompkf_eca
  logical                     , intent(in) :: is_surf
  real(r8)                    , intent(in) :: pct_sand
  real(r8)                    , intent(in) :: pct_clay
  real(r8)                    , intent(inout) :: cascade_matrix(centurybgc_index%nstvars, centurybgc_index%nreactions)
  real(r8)                    , intent(out) :: k_decay(1:ncentpools)
  real(r8)                    , intent(out) :: pot_co2_hr
  real(r8)                    , intent(out) :: alpha_n(1:ncentpools)
  real(r8)                    , intent(out) :: alpha_p(1:ncentpools)
  type(betr_status_type)      , intent(out) :: bstatus

  !local variables
  real(r8) :: pot_om_decay_rates(1:ncentpools)
  integer :: kc, jj

  associate(                                      &
    nelms => centurybgc_index%nelms,              &
    nom_tot_elms=> centurybgc_index%nom_tot_elms, &
    c_loc => centurybgc_index%c_loc               &
  )
  call bstatus%reset()

  call this%calc_cnp_ratios(centurybgc_index, ystates)

  !calculate potential decay coefficient (1/s)
  call this%calc_som_decay_k(centurybgc_index, decompkf_eca, k_decay)

  !calculate potential decay rates (mol C / s)
  call this%calc_som_decay_r(centurybgc_index, dtime, k_decay(1:ncentpools), &
      ystates(1:nom_tot_elms), pot_om_decay_rates)

  do jj = 1, ncentpools
    kc = (jj-1) * nelms + c_loc
    !the following avoids over-estimation of potential hr which is used for nitri-denit estimation
    pot_om_decay_rates(jj) = min(pot_om_decay_rates(jj), ystates(kc)/dtime)
  enddo

  call this%calc_cascade_matrix(is_surf, centurybgc_index, pct_sand, pct_clay, alpha_n, alpha_p, cascade_matrix)

  !calculate potential respiration rates by summarizing all om decomposition pathways
  call this%calc_potential_aerobic_hr(centurybgc_index, pot_om_decay_rates, &
    cascade_matrix, pot_co2_hr, bstatus)
  end associate
  end subroutine run_decomp
!------------------------------------------

  subroutine calc_cascade_matrix(this, is_surf, centurybgc_index, pct_sand, pct_clay, alpha_n, alpha_p, cascade_matrix)
  !
  ! DESCRIPTION
  ! calculate cascade matrix for decomposition
  ! in all the reactions, the nominal carbon oxidation status is assumed as zero, which is apparently not correct.
  ! It is also assumed the recycling of nitrogen and phosphorus during decomposition is 100%, which is likely
  ! not quite right as well.
  use BgcCentCnpIndexType , only : centurybgc_index_type
  use MathfuncMod         , only : safe_div, fpmax
  implicit none
  class(CentSom_type),           intent(inout) :: this
  type(centurybgc_index_type)  , intent(in)    :: centurybgc_index
  logical                      , intent(in)    :: is_surf
  real(r8)                     , intent(in)    :: pct_sand
  real(r8)                     , intent(in)    :: pct_clay
  real(r8)                     , intent(out)   :: alpha_n(ncentpools) !indicating factor for nitrogen limitation
  real(r8)                     , intent(out)   :: alpha_p(ncentpools) !indicating factor for phosphorus limitation
  real(r8)                     , intent(inout) :: cascade_matrix(centurybgc_index%nstvars, centurybgc_index%nreactions)

  integer  :: reac,jj
  real(r8) :: f1, f2, rf_s1

  associate(                                                   &
    lit1      => centurybgc_index%lit1                       , & !
    lit2      => centurybgc_index%lit2                       , & !
    lit3      => centurybgc_index%lit3                       , & !
    som1      => centurybgc_index%som1                       , & !
    som2      => centurybgc_index%som2                       , & !
    som3      => centurybgc_index%som3                       , & !
    cwd       => centurybgc_index%cwd                        , & !
    lwd       => centurybgc_index%lwd                        , & !
    fwd       => centurybgc_index%fwd                        , & !
    c_loc     => centurybgc_index%c_loc                      , & !
    n_loc     => centurybgc_index%n_loc                      , & !
    p_loc     => centurybgc_index%p_loc                      , & !
    c13_loc   => centurybgc_index%c13_loc                    , & !
    c14_loc   => centurybgc_index%c14_loc                    , & !
    nelms     => centurybgc_index%nelms                      , & !
    lid_o2    => centurybgc_index%lid_o2                     , & !
    lid_co2   => centurybgc_index%lid_co2                    , & !
    lid_nh4   => centurybgc_index%lid_nh4                    , & !
    lid_c14_co2=> centurybgc_index%lid_c14_co2               , & !
    lid_c13_co2=> centurybgc_index%lid_c13_co2               , & !
    lid_co2_hr => centurybgc_index%lid_co2_hr                , &
    lid_minn_nh4_immob=> centurybgc_index%lid_minn_nh4_immob , &
    lid_minp_immob => centurybgc_index%lid_minp_immob        , &
    lid_minp_soluble=> centurybgc_index%lid_minp_soluble     , &
    lit1_dek_reac => centurybgc_index%lit1_dek_reac          , &
    lit2_dek_reac => centurybgc_index%lit2_dek_reac          , &
    lit3_dek_reac => centurybgc_index%lit3_dek_reac          , &
    som1_dek_reac => centurybgc_index%som1_dek_reac          , &
    som2_dek_reac => centurybgc_index%som2_dek_reac          , &
    som3_dek_reac => centurybgc_index%som3_dek_reac          , &
    cwd_dek_reac => centurybgc_index%cwd_dek_reac            , &
    lwd_dek_reac => centurybgc_index%lwd_dek_reac            , &
    fwd_dek_reac => centurybgc_index%fwd_dek_reac            , &
    cwd_fcel     => this%cwd_fcel                            , &
    cwd_flig     => this%cwd_flig                            , &
    lwd_fcel     => this%lwd_fcel                            , &
    lwd_flig     => this%lwd_flig                            , &
    fwd_fcel     => this%fwd_fcel                            , &
    fwd_flig     => this%fwd_flig                            , &
    rf_l2s1_bgc  => this%rf_l2s1_bgc                         , &
    rf_l3s2_bgc  => this%rf_l3s2_bgc                         , &
    rf_s2s1_bgc  => this%rf_s2s1_bgc                         , &
    rf_s3s1_bgc  => this%rf_s3s1_bgc                         , &
    rf_l1s1_bgc  => this%rf_l1s1_bgc                         , &
    debug        => centurybgc_index%debug                     &
  )
    alpha_n = 0._r8; alpha_p = 0._r8
    !---------------------------------------------------------------------------------
    !reaction1, lit1 -> s1
    reac=lit1_dek_reac
    !lit1 + 0.55*o2 -> 0.45 som1 + 0.55co2 + (1/cn_ratios(lit1) - 0.45/cn_ratios(som1))min_n+ (1/cp_ratios(lit1)-0.45/cp_ratios(som1))min_p
    cascade_matrix((lit1-1)*nelms+c_loc   ,reac)  = -1._r8
    cascade_matrix((lit1-1)*nelms+n_loc   ,reac)  = -safe_div(1._r8,this%cn_ratios(lit1))
    cascade_matrix((lit1-1)*nelms+p_loc   ,reac)  = -safe_div(1._r8,this%cp_ratios(lit1))

    cascade_matrix((som1-1)*nelms+c_loc   ,reac)  = 1._r8-rf_l1s1_bgc
    cascade_matrix((som1-1)*nelms+n_loc   ,reac)  = safe_div(1._r8-rf_l1s1_bgc,this%cn_ratios(som1))
    cascade_matrix((som1-1)*nelms+p_loc   ,reac)  = safe_div(1._r8-rf_l1s1_bgc,this%cp_ratios(som1))

    cascade_matrix(lid_co2                ,reac)  = -cascade_matrix((lit1-1)*nelms+c_loc   ,reac)- &
                                                     cascade_matrix((som1-1)*nelms+c_loc   ,reac)

    cascade_matrix(lid_o2                 ,reac)  = -cascade_matrix(lid_co2                ,reac)
    cascade_matrix(lid_nh4                ,reac)  = -cascade_matrix((lit1-1)*nelms+n_loc   ,reac) - &
        cascade_matrix((som1-1)*nelms+n_loc   ,reac)
    cascade_matrix(lid_minp_soluble         ,reac)  = -cascade_matrix((lit1-1)*nelms+p_loc   ,reac) - &
        cascade_matrix((som1-1)*nelms+p_loc   ,reac)

    cascade_matrix(lid_minn_nh4_immob     ,reac)  = -cascade_matrix(lid_nh4         ,reac)
    cascade_matrix(lid_co2_hr             ,reac) = cascade_matrix(lid_co2           ,reac)
    cascade_matrix(lid_minp_immob         ,reac)  = -cascade_matrix(lid_minp_soluble  ,reac)

    if(this%use_c14)then
      cascade_matrix((lit1-1)*nelms+c14_loc   , reac) = -safe_div(1._r8,this%cc14_ratios(lit1))
      cascade_matrix(lid_c14_co2              , reac) = safe_div(rf_l1s1_bgc,this%cc14_ratios(lit1))
      cascade_matrix((som1-1)*nelms+c14_loc   , reac) = safe_div(1._r8-rf_l1s1_bgc,this%cc14_ratios(lit1))
    endif

    if(this%use_c13)then
      cascade_matrix((lit1-1)*nelms+c13_loc   , reac) = -safe_div(1._r8,this%cc13_ratios(lit1))
      cascade_matrix(lid_c13_co2              , reac) = safe_div(rf_l1s1_bgc,this%cc13_ratios(lit1))
      cascade_matrix((som1-1)*nelms+c13_loc   , reac) = safe_div(1._r8-rf_l1s1_bgc,this%cc13_ratios(lit1))
    endif

    if (cascade_matrix(lid_nh4, reac) < 0._r8)alpha_n(reac)=1._r8
    if (cascade_matrix(lid_minp_soluble,reac) < 0._r8)alpha_p(reac)=1._r8
    if(debug)then
      write(*,*)'lit1 carbon',cascade_matrix((lit1-1)*nelms+c_loc   ,reac)+cascade_matrix((som1-1)*nelms+c_loc   ,reac)+&
        cascade_matrix(lid_co2             ,reac)
      write(*,*)'lit1 nitrogen',cascade_matrix((lit1-1)*nelms+n_loc   ,reac)+cascade_matrix((som1-1)*nelms+n_loc   ,reac)+&
       cascade_matrix(lid_nh4         ,reac)
      write(*,*)'lit2 phosp',cascade_matrix((lit1-1)*nelms+p_loc   ,reac)+cascade_matrix((som1-1)*nelms+p_loc   ,reac)+&
       cascade_matrix(lid_minp_soluble        ,reac)
    endif
    !---------------------------------------------------------------------------------
    !reaction 2, lit2 -> s1
    reac = lit2_dek_reac
    !lit2 + 0.5 o2  -> 0.5 som1 + 0.5 co2 + (1/cn_ratios(lit2)-0.5/cn_ratios(som1))min_n +(1/cp_ratios(lit2)-0.5/cp_ratios(som1))min_p
    cascade_matrix((lit2-1)*nelms+c_loc   ,reac)   = -1._r8
    cascade_matrix((lit2-1)*nelms+n_loc   ,reac)   = -safe_div(1._r8,this%cn_ratios(lit2))
    cascade_matrix((lit2-1)*nelms+p_loc   ,reac)   = -safe_div(1._r8,this%cp_ratios(lit2))

    cascade_matrix((som1-1)*nelms+c_loc   ,reac)   =  1._r8-rf_l2s1_bgc
    cascade_matrix((som1-1)*nelms+n_loc   ,reac)   =  safe_div(1._r8-rf_l2s1_bgc,this%cn_ratios(som1))
    cascade_matrix((som1-1)*nelms+p_loc   ,reac)   =  safe_div(1._r8-rf_l2s1_bgc,this%cp_ratios(som1))

    cascade_matrix(lid_co2                ,reac)   =  -cascade_matrix((lit2-1)*nelms+c_loc   ,reac) - &
                                                       cascade_matrix((som1-1)*nelms+c_loc   ,reac)
    cascade_matrix(lid_o2                 ,reac)   = -cascade_matrix(lid_co2   ,reac)
    cascade_matrix(lid_nh4                ,reac)   = -cascade_matrix((lit2-1)*nelms+n_loc   ,reac) - &
                                                      cascade_matrix((som1-1)*nelms+n_loc   ,reac)

    cascade_matrix(lid_minp_soluble         ,reac)   = -cascade_matrix((lit2-1)*nelms+p_loc   ,reac) - &
                                                       cascade_matrix((som1-1)*nelms+p_loc   ,reac)

    cascade_matrix(lid_minn_nh4_immob     ,reac)   = -cascade_matrix(lid_nh4         ,reac)
    cascade_matrix(lid_minp_immob         ,reac)   = -cascade_matrix(lid_minp_soluble  ,reac)
    cascade_matrix(lid_co2_hr             ,reac)   = cascade_matrix(lid_co2        ,reac)


    if(cascade_matrix(lid_nh4, reac) < 0._r8)alpha_n(reac)=1._r8
    if(cascade_matrix(lid_minp_soluble,reac) < 0._r8)alpha_p(reac)=1._r8

    if(this%use_c14)then
      cascade_matrix((lit2-1)*nelms+c14_loc   , reac) = -safe_div(1._r8,this%cc14_ratios(lit2))
      cascade_matrix(lid_c14_co2              , reac) = safe_div(rf_l2s1_bgc,this%cc14_ratios(lit2))
      cascade_matrix((som1-1)*nelms+c14_loc   , reac) = safe_div(1._r8-rf_l2s1_bgc,this%cc14_ratios(lit2))
    endif
    if(this%use_c13)then
      cascade_matrix((lit2-1)*nelms+c13_loc   , reac) = -safe_div(1._r8,this%cc13_ratios(lit2))
      cascade_matrix(lid_c13_co2              , reac) =  safe_div(rf_l2s1_bgc,this%cc13_ratios(lit2))
      cascade_matrix((som1-1)*nelms+c13_loc   , reac) =  safe_div(1._r8-rf_l2s1_bgc,this%cc13_ratios(lit2))
    endif
    if(debug)then
      write(*,*)'lit2 carbon',cascade_matrix((lit2-1)*nelms+c_loc   ,reac) + cascade_matrix((som1-1)*nelms+c_loc   ,reac)+&
         cascade_matrix(lid_co2                ,reac)
      write(*,*)'lit2 nitrogen',cascade_matrix((lit2-1)*nelms+n_loc   ,reac) + cascade_matrix((som1-1)*nelms+n_loc   ,reac)+&
        cascade_matrix(lid_nh4                ,reac)
      write(*,*)'lit2 phosp',cascade_matrix((lit2-1)*nelms+p_loc   ,reac) + cascade_matrix((som1-1)*nelms+p_loc   ,reac)+&
        cascade_matrix(lid_minp_soluble                ,reac)
    endif
    !---------------------------------------------------------------------------------
    !reaction 3, lit3->s2
    reac = lit3_dek_reac
    !lit3 + 0.5 o2 -> 0.5 som2 + 0.5 co2 + (1/cn_ratios(lit3) - 0.5/cn_ratios(som2))min_n + (1/cp_ratios(lit3)-0.5_r8/cp_ratios(som2))minp
    cascade_matrix((lit3-1)*nelms+c_loc   ,reac) = -1._r8
    cascade_matrix((lit3-1)*nelms+n_loc   ,reac) = -safe_div(1._r8,this%cn_ratios(lit3))
    cascade_matrix((lit3-1)*nelms+p_loc   ,reac) = -safe_div(1._r8,this%cp_ratios(lit3))

    cascade_matrix((som2-1)*nelms+c_loc   ,reac) =  1._r8-rf_l3s2_bgc
    cascade_matrix((som2-1)*nelms+n_loc   ,reac) =  safe_div(1._r8-rf_l3s2_bgc,this%cn_ratios(som2))
    cascade_matrix((som2-1)*nelms+p_loc   ,reac) =  safe_div(1._r8-rf_l3s2_bgc,this%cp_ratios(som2))

    cascade_matrix(lid_co2                ,reac) = -cascade_matrix((lit3-1)*nelms+c_loc   ,reac) - &
                                                    cascade_matrix((som2-1)*nelms+c_loc   ,reac)
    cascade_matrix(lid_o2                 ,reac) = -cascade_matrix(lid_co2   ,reac)
    cascade_matrix(lid_nh4                ,reac) = -cascade_matrix((lit3-1)*nelms+n_loc   ,reac) - &
                                                    cascade_matrix((som2-1)*nelms+n_loc   ,reac)
    cascade_matrix(lid_minp_soluble      ,reac) = -cascade_matrix((lit3-1)*nelms+p_loc   ,reac) - &
                                                   cascade_matrix((som2-1)*nelms+p_loc   ,reac)

    cascade_matrix(lid_minn_nh4_immob     ,reac) = -cascade_matrix(lid_nh4         ,reac)
    cascade_matrix(lid_minp_immob         ,reac) = -cascade_matrix(lid_minp_soluble  ,reac)
    cascade_matrix(lid_co2_hr             ,reac) = cascade_matrix(lid_co2        ,reac)


    if (cascade_matrix(lid_nh4, reac) < 0._r8)alpha_n(reac)=1._r8
    if (cascade_matrix(lid_minp_soluble,reac) < 0._r8)alpha_p(reac)=1._r8

    if(this%use_c14)then
      cascade_matrix((lit3-1)*nelms+c14_loc   , reac) = -safe_div(1._r8,this%cc14_ratios(lit3))
      cascade_matrix(lid_c14_co2              , reac) =  safe_div(rf_l3s2_bgc,this%cc14_ratios(lit3))
      cascade_matrix((som2-1)*nelms+c14_loc   , reac) =  safe_div(1._r8-rf_l3s2_bgc,this%cc14_ratios(lit3))
    endif

    if(this%use_c13)then
      cascade_matrix((lit3-1)*nelms+c13_loc   , reac) = -safe_div(1._r8,this%cc13_ratios(lit3))
      cascade_matrix(lid_c13_co2              , reac) =  safe_div(rf_l3s2_bgc,this%cc13_ratios(lit3))
      cascade_matrix((som2-1)*nelms+c13_loc   , reac) =  safe_div(1._r8-rf_l3s2_bgc,this%cc13_ratios(lit3))
    endif
    if(debug)then
      write(*,*)'lit3 carbon',cascade_matrix((lit3-1)*nelms+c_loc   ,reac) + cascade_matrix((som2-1)*nelms+c_loc   ,reac) + &
         cascade_matrix(lid_co2                ,reac)
      write(*,*)'lit3 nitrogen',cascade_matrix((lit3-1)*nelms+n_loc   ,reac) + cascade_matrix((som2-1)*nelms+n_loc   ,reac) + &
         cascade_matrix(lid_nh4                ,reac)
      write(*,*)'lit3 phosp',cascade_matrix((lit3-1)*nelms+p_loc   ,reac) + cascade_matrix((som2-1)*nelms+p_loc   ,reac) + &
         cascade_matrix(lid_minp_soluble                ,reac)
    endif
    !---------------------------------------------------------------------------------
    !double check those stoichiometry parameters
    !reaction 4, the partition into som2 and som3 is soil texture dependent
    !SOM1 -> f1*SOM2 + f2*SOm3 + rf_s1*CO2 + (1/cn_ratios(SOM1)-f1/cn_ratios(SOM2)-f2/cn_ratios(SOM3))*min_n
    ! +(1/cp_ratios(SOM1)-f1/cp_ratios(SOM2)-f2/cp_ratios(SOM3))*min_p
    reac = som1_dek_reac
    f2=0.003_r8 + 0.00032_r8*pct_clay
    if(is_surf)then
      rf_s1 = 0.6_r8
    else
      rf_s1 = 0.17_r8 + 0.0068_r8*pct_sand
    endif

    f1 = 1._r8 - rf_s1 - f2

    cascade_matrix((som1-1)*nelms+c_loc   ,reac)  = -1._r8
    cascade_matrix((som1-1)*nelms+n_loc   ,reac)  = -safe_div(1._r8,this%cn_ratios(som1))
    cascade_matrix((som1-1)*nelms+p_loc   ,reac)  = -safe_div(1._r8,this%cp_ratios(som1))

    cascade_matrix((som3-1)*nelms+c_loc   ,reac)  = f2
    cascade_matrix((som3-1)*nelms+n_loc   ,reac)  = safe_div(f2,this%cn_ratios(som3))
    cascade_matrix((som3-1)*nelms+p_loc   ,reac)  = safe_div(f2,this%cp_ratios(som3))

    cascade_matrix((som2-1)*nelms+c_loc   ,reac) = f1
    cascade_matrix((som2-1)*nelms+n_loc   ,reac) = safe_div(f1,this%cn_ratios(som2))
    cascade_matrix((som2-1)*nelms+p_loc   ,reac) = safe_div(f1,this%cp_ratios(som2))

    cascade_matrix(lid_co2, reac)     = -cascade_matrix((som1-1)*nelms+c_loc   ,reac) - &
                                         cascade_matrix((som2-1)*nelms+c_loc   ,reac) - &
                                         cascade_matrix((som3-1)*nelms+c_loc   ,reac)

    cascade_matrix(lid_o2                 ,reac) = -cascade_matrix(lid_co2, reac)
    cascade_matrix(lid_nh4                ,reac) = -cascade_matrix((som1-1)*nelms+n_loc   ,reac) - &
           cascade_matrix((som2-1)*nelms+n_loc   ,reac)- &
           cascade_matrix((som3-1)*nelms+n_loc   ,reac)

    cascade_matrix(lid_minp_soluble         ,reac) = -cascade_matrix((som1-1)*nelms+p_loc   ,reac)-&
           cascade_matrix((som2-1)*nelms+p_loc   ,reac)- &
           cascade_matrix((som3-1)*nelms+p_loc   ,reac)

    cascade_matrix(lid_minn_nh4_immob     ,reac) = -cascade_matrix(lid_nh4         ,reac)
    cascade_matrix(lid_minp_immob         ,reac) = -cascade_matrix(lid_minp_soluble  ,reac)
    cascade_matrix(lid_co2_hr             ,reac) = cascade_matrix(lid_co2        ,reac)


    if (cascade_matrix(lid_nh4, reac) < 0._r8)alpha_n(reac)=1._r8
    if (cascade_matrix(lid_minp_soluble,reac) < 0._r8)alpha_p(reac)=1._r8

    if(this%use_c14)then
      cascade_matrix((som1-1)*nelms+c14_loc   , reac) = -safe_div(1._r8,this%cc14_ratios(som1))
      cascade_matrix(lid_c14_co2              , reac) =  safe_div(rf_s1,this%cc14_ratios(som1))
      cascade_matrix((som2-1)*nelms+c14_loc   , reac) =  safe_div(f1,this%cc14_ratios(som1))
      cascade_matrix((som3-1)*nelms+c14_loc   , reac) =  safe_div(f2,this%cc14_ratios(som1))
    endif

    if(this%use_c13)then
      cascade_matrix((som1-1)*nelms+c13_loc   , reac) = -safe_div(1._r8,this%cc13_ratios(som1))
      cascade_matrix(lid_c13_co2              , reac) = safe_div(rf_s1,this%cc13_ratios(som1))
      cascade_matrix((som2-1)*nelms+c13_loc   , reac) = safe_div(f1,this%cc13_ratios(som1))
      cascade_matrix((som3-1)*nelms+c13_loc   , reac) = safe_div(f2,this%cc13_ratios(som1))
    endif

    if(debug)then
      write(*,*)'som1 carbon',cascade_matrix((som1-1)*nelms+c_loc   ,reac) +cascade_matrix((som2-1)*nelms+c_loc   ,reac)+&
         cascade_matrix((som3-1)*nelms+c_loc   ,reac) + cascade_matrix(lid_co2                ,reac)
      write(*,*)'som1 nitrogen',cascade_matrix((som1-1)*nelms+n_loc   ,reac) +cascade_matrix((som2-1)*nelms+n_loc   ,reac)+&
         cascade_matrix((som3-1)*nelms+n_loc   ,reac) + cascade_matrix(lid_nh4                ,reac)
      write(*,*)'som1 phosp',cascade_matrix((som1-1)*nelms+p_loc   ,reac) +cascade_matrix((som2-1)*nelms+p_loc   ,reac)+&
         cascade_matrix((som3-1)*nelms+p_loc   ,reac) + cascade_matrix(lid_minp_soluble       ,reac)
    endif
    !---------------------------------------------------------------------------------
    !reaction 5, som2->som1, som3
    reac = som2_dek_reac
    !som2 + 0.55 o2 -> (0.45-f1) som1 + f1*som3 + 0.55co2 + (1/cn_ratios(som2)-0.42/cn_ratios(som1)-0.03/cn_ratios(som3)) + (1/cp_raitos(som2)-0.42/cp_ratios(som1)-0.03/cp_ratios(som3))
    f1 = 0.003_r8+0.00009_r8*pct_clay
    cascade_matrix((som2-1)*nelms+c_loc   ,reac)   = -1._r8
    cascade_matrix((som2-1)*nelms+n_loc   ,reac)   = -safe_div(1._r8,this%cn_ratios(som2))
    cascade_matrix((som2-1)*nelms+p_loc   ,reac)   = -safe_div(1._r8,this%cp_ratios(som2))


    cascade_matrix((som1-1)*nelms+c_loc   ,reac)   =  1._r8-rf_s2s1_bgc-f1
    cascade_matrix((som1-1)*nelms+n_loc   ,reac)   =  safe_div(1._r8-rf_s2s1_bgc-f1,this%cn_ratios(som1))
    cascade_matrix((som1-1)*nelms+p_loc   ,reac)   =  safe_div(1._r8-rf_s2s1_bgc-f1,this%cp_ratios(som1))

    cascade_matrix((som3-1)*nelms+c_loc   ,reac)   =  f1
    cascade_matrix((som3-1)*nelms+n_loc   ,reac)   =  safe_div(f1,this%cn_ratios(som3))
    cascade_matrix((som3-1)*nelms+p_loc   ,reac)   =  safe_div(f1,this%cp_ratios(som3))

    cascade_matrix(lid_co2                ,reac)   = -cascade_matrix((som2-1)*nelms+c_loc   ,reac) - &
                                                      cascade_matrix((som1-1)*nelms+c_loc   ,reac) - &
                                                      cascade_matrix((som3-1)*nelms+c_loc   ,reac)
    cascade_matrix(lid_o2                 ,reac)   = -cascade_matrix(lid_co2                ,reac)
    cascade_matrix(lid_nh4                ,reac)   =  -cascade_matrix((som2-1)*nelms+n_loc   ,reac)-&
                                             cascade_matrix((som1-1)*nelms+n_loc   ,reac) -&
                                             cascade_matrix((som3-1)*nelms+n_loc   ,reac)
    cascade_matrix(lid_minp_soluble         ,reac) = -cascade_matrix((som2-1)*nelms+p_loc   ,reac) -&
                                            cascade_matrix((som1-1)*nelms+p_loc   ,reac) - &
                                            cascade_matrix((som3-1)*nelms+p_loc   ,reac)

    cascade_matrix(lid_minn_nh4_immob     ,reac)   = -cascade_matrix(lid_nh4         ,reac)
    cascade_matrix(lid_minp_immob         ,reac)   = -cascade_matrix(lid_minp_soluble  ,reac)
    cascade_matrix(lid_co2_hr             ,reac)   = cascade_matrix(lid_co2        ,reac)

    if (cascade_matrix(lid_nh4, reac) < 0._r8)alpha_n(reac)=1._r8
    if (cascade_matrix(lid_minp_soluble,reac) < 0._r8)alpha_p(reac)=1._r8

    if(this%use_c14)then
      cascade_matrix((som2-1)*nelms+c14_loc   , reac) = -safe_div(1._r8,this%cc14_ratios(som2))
      cascade_matrix(lid_c14_co2              , reac) = safe_div(rf_s2s1_bgc,this%cc14_ratios(som2))
      cascade_matrix((som1-1)*nelms+c14_loc   , reac) = safe_div(1._r8-rf_s2s1_bgc-f1,this%cc14_ratios(som2))
      cascade_matrix((som3-1)*nelms+c14_loc   , reac) = safe_div(f1,this%cc14_ratios(som2))
    endif

    if(this%use_c13)then
      cascade_matrix((som2-1)*nelms+c13_loc   , reac) = -safe_div(1._r8,this%cc13_ratios(som2))
      cascade_matrix(lid_c13_co2              , reac) = safe_div(rf_s2s1_bgc,this%cc13_ratios(som2))
      cascade_matrix((som1-1)*nelms+c13_loc   , reac) = safe_div(1._r8-rf_s2s1_bgc-f1,this%cc13_ratios(som2))
      cascade_matrix((som3-1)*nelms+c13_loc   , reac) = safe_div(f1,this%cc13_ratios(som2))
    endif
    if(debug)then
      write(*,*)'som2 carbon',cascade_matrix((som2-1)*nelms+c_loc   ,reac) + cascade_matrix((som1-1)*nelms+c_loc   ,reac) + &
        cascade_matrix((som3-1)*nelms+c_loc   ,reac) + cascade_matrix(lid_co2                ,reac)
      write(*,*)'som2 nitrogen',cascade_matrix((som2-1)*nelms+n_loc   ,reac) + cascade_matrix((som1-1)*nelms+n_loc   ,reac) + &
        cascade_matrix((som3-1)*nelms+n_loc   ,reac) + cascade_matrix(lid_nh4                ,reac)
      write(*,*)'som2 phosp',cascade_matrix((som2-1)*nelms+p_loc   ,reac) + cascade_matrix((som1-1)*nelms+p_loc   ,reac) + &
        cascade_matrix((som3-1)*nelms+p_loc   ,reac) + cascade_matrix(lid_minp_soluble                ,reac)
    endif
    !---------------------------------------------------------------------------------
    !reaction 6, s3-> s1
    reac = som3_dek_reac
    !som3 + 0.55 o2 -> 0.45*som1 + 0.55co2 + (1/cn_ratios(som3)-0.45/cn_ratios(som1)) + (1/cp_ratios(som3)-0.45/cp_ratios(som1))
    cascade_matrix((som3-1)*nelms+c_loc   ,reac) = -1._r8
    cascade_matrix((som3-1)*nelms+n_loc   ,reac) = -safe_div(1._r8,this%cn_ratios(som3))
    cascade_matrix((som3-1)*nelms+p_loc   ,reac) = -safe_div(1._r8,this%cp_ratios(som3))

    cascade_matrix((som1-1)*nelms+c_loc   ,reac) = 1._r8-rf_s3s1_bgc
    cascade_matrix((som1-1)*nelms+n_loc   ,reac) = safe_div(1._r8-rf_s3s1_bgc,this%cn_ratios(som1))
    cascade_matrix((som1-1)*nelms+p_loc   ,reac) = safe_div(1._r8-rf_s3s1_bgc,this%cp_ratios(som1))

    cascade_matrix(lid_co2                ,reac) = -cascade_matrix((som3-1)*nelms+c_loc   ,reac) - &
                                                    cascade_matrix((som1-1)*nelms+c_loc   ,reac)
    cascade_matrix(lid_o2                 ,reac) = -cascade_matrix(lid_co2                ,reac)
    cascade_matrix(lid_nh4                ,reac) = -cascade_matrix((som3-1)*nelms+n_loc   ,reac)  - &
                                                    cascade_matrix((som1-1)*nelms+n_loc   ,reac)

    cascade_matrix(lid_minp_soluble       ,reac) = -cascade_matrix((som3-1)*nelms+p_loc   ,reac) - &
                                                    cascade_matrix((som1-1)*nelms+p_loc   ,reac)

    cascade_matrix(lid_minn_nh4_immob     ,reac) = -cascade_matrix(lid_nh4         ,reac)
    cascade_matrix(lid_minp_immob         ,reac) = -cascade_matrix(lid_minp_soluble  ,reac)
    cascade_matrix(lid_co2_hr             ,reac) = cascade_matrix(lid_co2        ,reac)


    if (cascade_matrix(lid_nh4, reac) < 0._r8)alpha_n(reac)=1._r8
    if (cascade_matrix(lid_minp_soluble,reac) < 0._r8)alpha_p(reac)=1._r8

    if(this%use_c14)then
      cascade_matrix((som3-1)*nelms+c14_loc   , reac) = -safe_div(1._r8,this%cc14_ratios(som3))
      cascade_matrix(lid_c14_co2              , reac) =  safe_div(rf_s3s1_bgc, this%cc14_ratios(som3))
      cascade_matrix((som1-1)*nelms+c14_loc   , reac) =  safe_div(1._r8-rf_s3s1_bgc,this%cc14_ratios(som3))
    endif

    if(this%use_c13)then
      cascade_matrix((som3-1)*nelms+c13_loc   , reac) = -safe_div(1._r8,this%cc13_ratios(som3))
      cascade_matrix(lid_c13_co2              , reac) =  safe_div(rf_s3s1_bgc,this%cc13_ratios(som3))
      cascade_matrix((som1-1)*nelms+c13_loc   , reac) =  safe_div(1._r8-rf_s3s1_bgc,this%cc13_ratios(som3))
    endif
    if(debug)then
      write(*,*)'som3 carbon',cascade_matrix((som3-1)*nelms+c_loc   ,reac)  + cascade_matrix((som1-1)*nelms+c_loc   ,reac) + &
         cascade_matrix(lid_co2                ,reac)
      write(*,*)'som3 nitrogen',cascade_matrix((som3-1)*nelms+n_loc   ,reac)  + cascade_matrix((som1-1)*nelms+n_loc   ,reac) + &
         cascade_matrix(lid_nh4                ,reac)
      write(*,*)'som3 phosp',cascade_matrix((som3-1)*nelms+p_loc   ,reac)  + cascade_matrix((som1-1)*nelms+p_loc   ,reac) + &
         cascade_matrix(lid_minp_soluble                ,reac)
    !print*,'cp_ratio',maxval(this%cp_ratios),minval(this%cp_ratios)
    !print*,'cn_ratio',maxval(this%cn_ratios),minval(this%cn_ratios)

     endif

    !---------------------------------------------------------------------------------
    !reaction 7, the partition cwd into som1 and som2
    reac = cwd_dek_reac
    !cwd + o2 -> (1-flig)((1-rf_l2s1_bgc)*SOM1+rf_l2s1_bgc*CO2) + flig*((1-rf_l3s2_bgc)*SOM2+rf_l3s2_bgc*CO2)
    !    + (1/cn_ratios(cwd)-f1/cn_ratios(som1)-f2/cn_ratios(som2))
    !    + (1/cp_ratios(cwd)-f1/cp_ratios(som1)-f2/cp_ratios(som2))
    f1 = cwd_fcel*(1-rf_l2s1_bgc)
    f2 = (1._r8-cwd_fcel)*(1-rf_l3s2_bgc)

    cascade_matrix((cwd-1)*nelms+c_loc    ,reac) = -1._r8
    cascade_matrix((cwd-1)*nelms+n_loc    ,reac) = -safe_div(1._r8,this%cn_ratios(cwd))
    cascade_matrix((cwd-1)*nelms+p_loc    ,reac) = -safe_div(1._r8,this%cp_ratios(cwd))

    cascade_matrix((som1-1)*nelms+c_loc   ,reac) = f1
    cascade_matrix((som1-1)*nelms+n_loc   ,reac) = safe_div(f1,this%cn_ratios(som1))
    cascade_matrix((som1-1)*nelms+p_loc   ,reac) = safe_div(f1,this%cp_ratios(som1))

    cascade_matrix((som2-1)*nelms+c_loc   ,reac) = f2
    cascade_matrix((som2-1)*nelms+n_loc   ,reac) = safe_div(f2,this%cn_ratios(lit3))
    cascade_matrix((som2-1)*nelms+p_loc   ,reac) = safe_div(f2,this%cp_ratios(lit3))

    cascade_matrix(lid_co2                ,reac) = - cascade_matrix((cwd-1)*nelms+c_loc    ,reac) - &
                                                     cascade_matrix((som1-1)*nelms+c_loc   ,reac) - &
                                                     cascade_matrix((som2-1)*nelms+c_loc   ,reac)
    cascade_matrix(lid_o2                 ,reac) = -cascade_matrix(lid_co2                ,reac)
    cascade_matrix(lid_nh4                ,reac) = -cascade_matrix((cwd-1)*nelms+n_loc    ,reac) - &
       cascade_matrix((som1-1)*nelms+n_loc   ,reac) - &
       cascade_matrix((som2-1)*nelms+n_loc   ,reac)

    cascade_matrix(lid_minp_soluble         ,reac) = -cascade_matrix((cwd-1)*nelms+p_loc    ,reac) - &
       cascade_matrix((som1-1)*nelms+p_loc   ,reac) - &
       cascade_matrix((som2-1)*nelms+p_loc   ,reac)

    cascade_matrix(lid_minn_nh4_immob     ,reac) = -cascade_matrix(lid_nh4         ,reac)
    cascade_matrix(lid_minp_immob         ,reac) = -cascade_matrix(lid_minp_soluble  ,reac)
    cascade_matrix(lid_co2_hr             ,reac) = cascade_matrix(lid_co2        ,reac)

    if (cascade_matrix(lid_nh4, reac) < 0._r8)alpha_n(reac)=1._r8
    if (cascade_matrix(lid_minp_soluble,reac) < 0._r8)alpha_p(reac)=1._r8

    if(this%use_c14)then
      cascade_matrix((cwd-1)*nelms+c14_loc   , reac) = -safe_div(1._r8,this%cc14_ratios(cwd))
      cascade_matrix((som1-1)*nelms+c14_loc  , reac) =  safe_div(f1,this%cc14_ratios(cwd))
      cascade_matrix((som2-1)*nelms+c14_loc  , reac) =  safe_div(f2,this%cc14_ratios(cwd))
    endif

    if(this%use_c14)then
      cascade_matrix((cwd-1)*nelms+c13_loc   , reac) = -safe_div(1._r8,this%cc13_ratios(cwd))
      cascade_matrix((som1-1)*nelms+c13_loc  , reac) =  safe_div(f1,this%cc13_ratios(cwd))
      cascade_matrix((som2-1)*nelms+c13_loc  , reac) =  safe_div(f2,this%cc13_ratios(cwd))
    endif
    if(debug)then
      !write(*,*)'cwd f1 f2',f1,f2
      write(*,*)'cwd carbon', cascade_matrix((cwd-1)*nelms+c_loc    ,reac) + cascade_matrix((som1-1)*nelms+c_loc   ,reac) + &
        cascade_matrix((som2-1)*nelms+c_loc   ,reac) + cascade_matrix(lid_co2                ,reac)
      write(*,*)'cwd nitrogen', cascade_matrix((cwd-1)*nelms+n_loc    ,reac) + cascade_matrix((som1-1)*nelms+n_loc   ,reac) + &
        cascade_matrix((som2-1)*nelms+n_loc   ,reac) +cascade_matrix(lid_nh4         ,reac)
      write(*,*)'cwd phosp', cascade_matrix((cwd-1)*nelms+p_loc    ,reac) + cascade_matrix((som1-1)*nelms+p_loc   ,reac) + &
        cascade_matrix((som2-1)*nelms+p_loc   ,reac) +cascade_matrix(lid_minp_soluble         ,reac)
    endif
    !---------------------------------------------------------------------------------
    !reaction 8, the partition lwd into som1 and som2
    reac = lwd_dek_reac
    !lwd + o2 -> (1-flig)((1-rf_l2s1_bgc)*SOM1+rf_l2s1_bgc*CO2) + flig*((1-rf_l3s2_bgc)*SOM2+rf_l3s2_bgc*CO2)
    !    + (1/cn_ratios(cwd)-f1/cn_ratios(som1)-f2/cn_ratios(som2))
    !    + (1/cp_ratios(cwd)-f1/cp_ratios(som1)-f2/cp_ratios(som2))
    f1 = lwd_fcel*(1-rf_l2s1_bgc)
    f2 = (1._r8-lwd_fcel)*(1-rf_l3s2_bgc)

    cascade_matrix((lwd-1)*nelms+c_loc    ,reac) = -1._r8
    cascade_matrix((lwd-1)*nelms+n_loc    ,reac) = -safe_div(1._r8,this%cn_ratios(lwd))
    cascade_matrix((lwd-1)*nelms+p_loc    ,reac) = -safe_div(1._r8,this%cp_ratios(lwd))

    cascade_matrix((som1-1)*nelms+c_loc   ,reac) = f1
    cascade_matrix((som1-1)*nelms+n_loc   ,reac) = safe_div(f1,this%cn_ratios(som1))
    cascade_matrix((som1-1)*nelms+p_loc   ,reac) = safe_div(f1,this%cp_ratios(som1))

    cascade_matrix((som2-1)*nelms+c_loc   ,reac) = f2
    cascade_matrix((som2-1)*nelms+n_loc   ,reac) = safe_div(f2,this%cn_ratios(som2))
    cascade_matrix((som2-1)*nelms+p_loc   ,reac) = safe_div(f2,this%cp_ratios(som2))

    cascade_matrix(lid_co2                ,reac) = -cascade_matrix((lwd-1)*nelms+c_loc    ,reac) - &
                                                    cascade_matrix((som1-1)*nelms+c_loc   ,reac) - &
                                                    cascade_matrix((som2-1)*nelms+c_loc   ,reac)
    cascade_matrix(lid_o2                 ,reac) = -cascade_matrix(lid_co2  ,reac)

    cascade_matrix(lid_nh4                ,reac) = -cascade_matrix((lwd-1)*nelms+n_loc    ,reac) - &
       cascade_matrix((som1-1)*nelms+n_loc   ,reac) - &
       cascade_matrix((som2-1)*nelms+n_loc   ,reac)

    cascade_matrix(lid_minp_soluble         ,reac) = -cascade_matrix((lwd-1)*nelms+p_loc    ,reac) - &
       cascade_matrix((som1-1)*nelms+p_loc   ,reac) - &
       cascade_matrix((som2-1)*nelms+p_loc   ,reac)

    cascade_matrix(lid_minn_nh4_immob     ,reac) = -cascade_matrix(lid_nh4         ,reac)
    cascade_matrix(lid_minp_immob         ,reac) = -cascade_matrix(lid_minp_soluble  ,reac)
    cascade_matrix(lid_co2_hr             ,reac) = cascade_matrix(lid_co2        ,reac)

    if (cascade_matrix(lid_nh4, reac) < 0._r8)alpha_n(reac)=1._r8
    if (cascade_matrix(lid_minp_soluble,reac) < 0._r8)alpha_p(reac)=1._r8

    if(this%use_c14)then
      cascade_matrix((lwd-1)*nelms+c14_loc   , reac) = -safe_div(1._r8,this%cc14_ratios(lwd))
      cascade_matrix((som1-1)*nelms+c14_loc  , reac) =  safe_div(f1,this%cc14_ratios(lwd))
      cascade_matrix((som2-1)*nelms+c14_loc  , reac) =  safe_div(f2,this%cc14_ratios(lwd))
    endif

    if(this%use_c14)then
      cascade_matrix((lwd-1)*nelms+c13_loc   , reac) = -safe_div(1._r8,this%cc13_ratios(lwd))
      cascade_matrix((som1-1)*nelms+c13_loc  , reac) =  safe_div(f1,this%cc13_ratios(lwd))
      cascade_matrix((som2-1)*nelms+c13_loc  , reac) =  safe_div(f2,this%cc13_ratios(lwd))
    endif
    if(debug)then
      !write(*,*)'lwd f1 f2', f1, f2
      write(*,*)'lwd carbon', cascade_matrix((lwd-1)*nelms+c_loc    ,reac) + cascade_matrix((som1-1)*nelms+c_loc   ,reac) + &
        cascade_matrix((som2-1)*nelms+c_loc   ,reac) + cascade_matrix(lid_co2                ,reac)
      write(*,*)'lwd nitrogen', cascade_matrix((lwd-1)*nelms+n_loc    ,reac) + cascade_matrix((som1-1)*nelms+n_loc   ,reac) + &
        cascade_matrix((som2-1)*nelms+n_loc   ,reac) +cascade_matrix(lid_nh4         ,reac)
      write(*,*)'lwd phosp', cascade_matrix((lwd-1)*nelms+p_loc    ,reac) + cascade_matrix((som1-1)*nelms+p_loc   ,reac) + &
        cascade_matrix((som2-1)*nelms+p_loc   ,reac) +cascade_matrix(lid_minp_soluble         ,reac)
    endif
    !---------------------------------------------------------------------------------
    !reaction 9, the partition fwd into som1 and som2
    reac = fwd_dek_reac
    !fwd + o2 -> (1-flig)((1-rf_l2s1_bgc)*SOM1+rf_l2s1_bgc*CO2) + flig*((1-rf_l3s2_bgc)*SOM2+rf_l3s2_bgc*CO2)
    !    + (1/cn_ratios(cwd)-f1/cn_ratios(som1)-f2/cn_ratios(som2))
    !    + (1/cp_ratios(cwd)-f1/cp_ratios(som1)-f2/cp_ratios(som2))
    f1 = fwd_fcel*(1-rf_l2s1_bgc)
    f2 = (1._r8-fwd_fcel)*(1-rf_l3s2_bgc)

    cascade_matrix((fwd-1)*nelms+c_loc    ,reac) = -1._r8
    cascade_matrix((fwd-1)*nelms+n_loc    ,reac) = -safe_div(1._r8,this%cn_ratios(fwd))
    cascade_matrix((fwd-1)*nelms+p_loc    ,reac) = -safe_div(1._r8,this%cp_ratios(fwd))

    cascade_matrix((som1-1)*nelms+c_loc   ,reac) = f1
    cascade_matrix((som1-1)*nelms+n_loc   ,reac) = safe_div(f1,this%cn_ratios(som1))
    cascade_matrix((som1-1)*nelms+p_loc   ,reac) = safe_div(f1,this%cp_ratios(som1))

    cascade_matrix((som2-1)*nelms+c_loc   ,reac) = f2
    cascade_matrix((som2-1)*nelms+n_loc   ,reac) = safe_div(f2,this%cn_ratios(som2))
    cascade_matrix((som2-1)*nelms+p_loc   ,reac) = safe_div(f2,this%cp_ratios(som2))

    cascade_matrix(lid_co2                ,reac) = -cascade_matrix((fwd-1)*nelms+c_loc    ,reac) - &
                                                    cascade_matrix((som1-1)*nelms+c_loc   ,reac) - &
                                                    cascade_matrix((som2-1)*nelms+c_loc   ,reac)
    cascade_matrix(lid_o2                 ,reac) = -cascade_matrix(lid_co2  ,reac)

    cascade_matrix(lid_nh4                ,reac) = -cascade_matrix((fwd-1)*nelms+n_loc    ,reac)  - &
       cascade_matrix((som1-1)*nelms+n_loc   ,reac) - &
       cascade_matrix((som2-1)*nelms+n_loc   ,reac)

    cascade_matrix(lid_minp_soluble         ,reac) = -cascade_matrix((fwd-1)*nelms+p_loc    ,reac) - &
       cascade_matrix((som1-1)*nelms+p_loc   ,reac) - &
       cascade_matrix((som2-1)*nelms+p_loc   ,reac)

    cascade_matrix(lid_minn_nh4_immob     ,reac) = -cascade_matrix(lid_nh4         ,reac)
    cascade_matrix(lid_minp_immob         ,reac) = -cascade_matrix(lid_minp_soluble  ,reac)
    cascade_matrix(lid_co2_hr             ,reac) = cascade_matrix(lid_co2        ,reac)

    if (cascade_matrix(lid_nh4, reac) < 0._r8)alpha_n(reac)=1._r8
    if (cascade_matrix(lid_minp_soluble,reac) < 0._r8)alpha_p(reac)=1._r8

    if(this%use_c14)then
      cascade_matrix((fwd-1)*nelms+c14_loc   , reac) = -safe_div(1._r8,this%cc14_ratios(fwd))
      cascade_matrix((som1-1)*nelms+c14_loc  , reac) =  safe_div(f1,this%cc14_ratios(fwd))
      cascade_matrix((som2-1)*nelms+c14_loc  , reac) =  safe_div(f2,this%cc14_ratios(fwd))
    endif

    if(this%use_c14)then
      cascade_matrix((fwd-1)*nelms+c13_loc   , reac) = -safe_div(1._r8,this%cc13_ratios(fwd))
      cascade_matrix((som1-1)*nelms+c13_loc  , reac) =  safe_div(f1,this%cc13_ratios(fwd))
      cascade_matrix((som2-1)*nelms+c13_loc  , reac) =  safe_div(f2,this%cc13_ratios(fwd))
    endif
    if(debug)then
      !write(*,*)'fwd f1 f2', f1, f2
      write(*,*)'fwd carbon', cascade_matrix((fwd-1)*nelms+c_loc    ,reac) + cascade_matrix((som1-1)*nelms+c_loc   ,reac) + &
        cascade_matrix((som2-1)*nelms+c_loc   ,reac) + cascade_matrix(lid_co2                ,reac)
      write(*,*)'fwd nitrogen', cascade_matrix((fwd-1)*nelms+n_loc    ,reac) + cascade_matrix((som1-1)*nelms+n_loc   ,reac) + &
        cascade_matrix((som2-1)*nelms+n_loc   ,reac) +cascade_matrix(lid_nh4         ,reac)
      write(*,*)'fwd phosp', cascade_matrix((fwd-1)*nelms+p_loc    ,reac) + cascade_matrix((som1-1)*nelms+p_loc   ,reac) + &
        cascade_matrix((som2-1)*nelms+p_loc   ,reac) +cascade_matrix(lid_minp_soluble ,reac)
    endif
  end associate
  end subroutine calc_cascade_matrix

  !-----------------------------------------------------------------------
  subroutine calc_potential_aerobic_hr(this, centurybgc_index, pot_decay_rates, &
    cascade_matrix, pot_co2_hr, bstatus)
    !
    ! DESCRIPTION:
    ! calculate potential aerobic heteorotrophic respiration, and potential oxygen consumption based on cascade_matrix
    ! !USES:
    use MathfuncMod         , only : dot_sum
    use MathfuncMod         , only : safe_div
    use BgcCentCnpIndexType       , only : centurybgc_index_type
    use BetrStatusType, only : betr_status_type
    implicit none
    ! !ARGUMENTS:
    class(CentSom_type), intent(inout) :: this
    type(centurybgc_index_type)   , intent(in) :: centurybgc_index
    real(r8)                , intent(in) :: pot_decay_rates(ncentpools)
    real(r8)                , intent(in) :: cascade_matrix(centurybgc_index%nstvars, centurybgc_index%nreactions)
    real(r8)                , intent(out):: pot_co2_hr
    type(betr_status_type)  , intent(out) :: bstatus
    ! !LOCAL VARIABLES:
    real(r8) :: cascade_matrix_hr(ncentpools)
    integer  :: reac

    associate(                                           & !
         nom_pools => centurybgc_index%nom_pools        , & !
         lid_co2_hr=> centurybgc_index%lid_co2_hr       , & !
         lit1      => centurybgc_index%lit1             , & !
         lit2      => centurybgc_index%lit2             , & !
         lit3      => centurybgc_index%lit3             , & !
         som1      => centurybgc_index%som1             , & !
         som2      => centurybgc_index%som2             , & !
         som3      => centurybgc_index%som3             , & !
         cwd       => centurybgc_index%cwd              , & !
         lwd       => centurybgc_index%lwd              , & !
         fwd       => centurybgc_index%fwd              , & !
         lit1_dek_reac=> centurybgc_index%lit1_dek_reac , & !
         lit2_dek_reac=> centurybgc_index%lit2_dek_reac , & !
         lit3_dek_reac=> centurybgc_index%lit3_dek_reac , & !
         som1_dek_reac=> centurybgc_index%som1_dek_reac , & !
         som2_dek_reac=> centurybgc_index%som2_dek_reac , & !
         som3_dek_reac=> centurybgc_index%som3_dek_reac , & !
         cwd_dek_reac=> centurybgc_index%cwd_dek_reac   , & !
         lwd_dek_reac=> centurybgc_index%lwd_dek_reac   , & !
         fwd_dek_reac=> centurybgc_index%fwd_dek_reac     & !
         )

    cascade_matrix_hr = 0._r8
    reac=lit1_dek_reac; cascade_matrix_hr(lit1)=cascade_matrix(lid_co2_hr,reac)
    reac=lit2_dek_reac; cascade_matrix_hr(lit2)=cascade_matrix(lid_co2_hr,reac)
    reac=lit3_dek_reac; cascade_matrix_hr(lit3)=cascade_matrix(lid_co2_hr,reac)
    reac=cwd_dek_reac ; cascade_matrix_hr(cwd) =cascade_matrix(lid_co2_hr,reac)
    reac=lwd_dek_reac ; cascade_matrix_hr(lwd) =cascade_matrix(lid_co2_hr,reac)
    reac=fwd_dek_reac ; cascade_matrix_hr(fwd) =cascade_matrix(lid_co2_hr,reac)
    reac=som1_dek_reac; cascade_matrix_hr(som1)=cascade_matrix(lid_co2_hr,reac)
    reac=som2_dek_reac; cascade_matrix_hr(som2)=cascade_matrix(lid_co2_hr,reac)
    reac=som3_dek_reac; cascade_matrix_hr(som3)=cascade_matrix(lid_co2_hr,reac)

    pot_co2_hr = dot_sum(cascade_matrix_hr, pot_decay_rates, bstatus)  !mol CO2/m3/s
    end associate
  end subroutine calc_potential_aerobic_hr

  !-----------------------------------------------------------------------
  subroutine calc_cnp_ratios(this, centurybgc_index, ystates)
  !
  ! DESCRIPTION
  ! compute the cnp ratios for the om pools
  use MathfuncMod         , only : safe_div
  use BgcCentCnpIndexType       , only : centurybgc_index_type
  implicit none
  class(CentSom_type), intent(inout) :: this
  type(centurybgc_index_type)   , intent(in) :: centurybgc_index
  real(r8), intent(in) :: ystates(centurybgc_index%nstvars)
  integer :: jj
  integer :: kc, kn, kp, kc1, kc2

  associate(                         &
    nelms => centurybgc_index%nelms, &
    c_loc => centurybgc_index%c_loc, &
    n_loc => centurybgc_index%n_loc, &
    p_loc => centurybgc_index%p_loc, &
    lit2  => centurybgc_index%lit2 , &
    lit3  => centurybgc_index%lit3   &
  )

  !for om pools
  do jj = 1, ncentpools
    kc = (jj-1) * nelms + c_loc
    kn = (jj-1) * nelms + n_loc
    kp = (jj-1) * nelms + p_loc

    if(ystates(kc)==0._r8)then
      this%cn_ratios(jj) = this%def_cn(jj)
      this%cp_ratios(jj) = this%def_cp(jj)
      if(this%use_c13)then
      endif
      if(this%use_c14)then
        this%cc14_ratios(jj) = this%def_cc14(jj)
      endif
      if(this%use_c13)then
        this%cc13_ratios(jj) = this%def_cc13(jj)
      endif
    else
      this%cn_ratios(jj) = safe_div(ystates(kc),ystates(kn))
      this%cp_ratios(jj) = safe_div(ystates(kc),ystates(kp))
    endif
  enddo
  kc1 = (lit2-1)*nelms+c_loc
  kc2 = (lit3-1)*nelms+c_loc
  !lignin fraction of the structural carbon
  this%lit_flig = safe_div(ystates(kc2),ystates(kc1)+ystates(kc2))

  end associate
  end subroutine calc_cnp_ratios

  !-------------------------------------------------------------------------------
  subroutine calc_som_decay_r(this, centurybgc_index, dtime, om_k_decay, om_pools, om_decay_rates)
    !
    ! !DESCRIPTION:
    ! calculate degradation for all different pools
    !
    ! !USES:
    use BgcCentCnpIndexType       , only : centurybgc_index_type
   implicit none
   class(CentSom_type)     , intent(inout) :: this
   type(centurybgc_index_type) , intent(in)    :: centurybgc_index
    real(r8)  , intent(in)    :: dtime
    real(r8)  , intent(in)    :: om_k_decay(ncentpools)
    real(r8)  , intent(in)    :: om_pools(centurybgc_index%nom_tot_elms)
    real(r8)  , intent(out)   :: om_decay_rates(ncentpools)

    ! !LOCAL VARIABLES:
    integer :: jj, fc, c, j
    integer :: kc, kn
    associate(                                        &
         nelms => centurybgc_index%nelms            , &
         nom_pools => centurybgc_index%nom_pools    , &
         c_loc => centurybgc_index%c_loc              &
    )

    !for om pools
    do jj = 1, nom_pools
      kc = (jj-1) * nelms + c_loc
      om_decay_rates(jj) = om_pools(kc) * om_k_decay(jj)
    enddo
    end associate
  end subroutine calc_som_decay_r

  !-------------------------------------------------------------------------------
  subroutine calc_som_decay_k(this, centurybgc_index, decompkf_eca, k_decay)

  use BgcCentCnpIndexType       , only : centurybgc_index_type
  use BgcCentCnpDecompType      , only : DecompCent_type
  implicit none
  class(CentSom_type)     , intent(inout) :: this
  type(DecompCent_type), intent(in) :: decompkf_eca
  type(centurybgc_index_type)   , intent(in)    :: centurybgc_index
  real(r8)                      , intent(out)   :: k_decay(ncentpools)
  integer :: jj

  associate(   &
   t_scalar       => decompkf_eca%t_scalar        , & ! Intput: [real(r8) (:,:)   ]  soil temperature scalar for decomp
   w_scalar       => decompkf_eca%w_scalar        , & ! Intput: [real(r8) (:,:)   ]  soil water scalar for decomp
   o_scalar       => decompkf_eca%o_scalar        , & ! Intput: [real(r8) (:,:)   ]  fraction by which decomposition is limited by anoxia
   depth_scalar   => decompkf_eca%depth_scalar   , & ! Intput: [real(r8) (:,:)   ]  rate constant for decomposition (1./sec)
   lit1           => centurybgc_index%lit1               , & !
   lit2           => centurybgc_index%lit2               , & !
   lit3           => centurybgc_index%lit3               , & !
   som1           => centurybgc_index%som1               , & !
   som2           => centurybgc_index%som2               , & !
   som3           => centurybgc_index%som3               , & !
   cwd            => centurybgc_index%cwd                , & !
   lwd            => centurybgc_index%lwd                , & !
   fwd            => centurybgc_index%fwd                  & !
  )
  k_decay(lit1) = this%k_decay_lit1 * t_scalar * w_scalar * o_scalar * depth_scalar
  k_decay(lit2) = this%k_decay_lit2 * t_scalar * w_scalar * o_scalar * depth_scalar
  k_decay(lit3) = this%k_decay_lit3 * t_scalar * w_scalar * o_scalar * depth_scalar
  k_decay(som1) = this%k_decay_som1 * t_scalar * w_scalar * o_scalar * depth_scalar
  k_decay(som2) = this%k_decay_som2 * t_scalar * w_scalar * o_scalar * depth_scalar
  k_decay(som3) = this%k_decay_som3 * t_scalar * w_scalar * o_scalar * depth_scalar
  k_decay(cwd)  = this%k_decay_cwd  * t_scalar * w_scalar * o_scalar * depth_scalar
  k_decay(lwd)  = this%k_decay_lwd  * t_scalar * w_scalar * o_scalar * depth_scalar
  k_decay(fwd)  = this%k_decay_fwd  * t_scalar * w_scalar * o_scalar * depth_scalar
  !impose the ligin effect
  k_decay(cwd)  = k_decay(cwd) * exp(-3._r8*this%cwd_flig)
  k_decay(lwd)  = k_decay(lwd) * exp(-3._r8*this%lwd_flig)
  k_decay(fwd)  = k_decay(fwd) * exp(-3._r8*this%fwd_flig)
  k_decay(lit2) = k_decay(lit2)* exp(-3._r8*this%lit_flig)
  k_decay(lit3) = k_decay(lit3)* exp(-3._r8*this%lit_flig)
  end associate
  end subroutine calc_som_decay_k
  !-------------------------------------------------------------------------------
  subroutine calc_pot_min_np_flx(this, dtime, centurybgc_index, ystates, k_decay, cascade_matrix, &
    alpha_n, alpha_p, pot_decomp, pot_nn_flx, pot_np_flx)
  use BgcCentCnpIndexType       , only : centurybgc_index_type
  implicit none
  class(CentSom_type)         , intent(inout) :: this
  real(r8)                    , intent(in) :: dtime
  type(centurybgc_index_type) , intent(in) :: centurybgc_index
  real(r8)                    , intent(in) :: ystates(1:centurybgc_index%nom_tot_elms)
  real(r8)                    , intent(in) :: k_decay(1:ncentpools)
  real(r8)                    , intent(in) :: cascade_matrix(centurybgc_index%nstvars, centurybgc_index%nreactions)
  real(r8)                    , intent(in) :: alpha_n(ncentpools)
  real(r8)                    , intent(in) :: alpha_p(ncentpools)
  real(r8)                    , intent(out) :: pot_decomp(ncentpools)
  real(r8)                    , intent(out):: pot_nn_flx
  real(r8)                    , intent(out):: pot_np_flx

  integer :: reac
  integer :: reacs(ncentpools)

  associate(                                                    & !
       nom_pools => centurybgc_index%nom_pools                , & !
       nom_tot_elms=> centurybgc_index%nom_tot_elms           , & !
       lid_nh4   => centurybgc_index%lid_nh4                  , & !
       lid_minp_soluble  => centurybgc_index%lid_minp_soluble , & !
       lit1      => centurybgc_index%lit1                     , & !
       lit2      => centurybgc_index%lit2                     , & !
       lit3      => centurybgc_index%lit3                     , & !
       som1      => centurybgc_index%som1                     , & !
       som2      => centurybgc_index%som2                     , & !
       som3      => centurybgc_index%som3                     , & !
       cwd       => centurybgc_index%cwd                      , & !
       lit1_dek_reac=> centurybgc_index%lit1_dek_reac         , & !
       lit2_dek_reac=> centurybgc_index%lit2_dek_reac         , & !
       lit3_dek_reac=> centurybgc_index%lit3_dek_reac         , & !
       som1_dek_reac=> centurybgc_index%som1_dek_reac         , & !
       som2_dek_reac=> centurybgc_index%som2_dek_reac         , & !
       som3_dek_reac=> centurybgc_index%som3_dek_reac         , & !
       cwd_dek_reac=> centurybgc_index%cwd_dek_reac           , & !
       lwd_dek_reac=> centurybgc_index%lwd_dek_reac           , & !
       fwd_dek_reac=> centurybgc_index%fwd_dek_reac             & !
   )

  !calculate potential decay rates (mol C / s)
  call this%calc_som_decay_r(centurybgc_index, dtime, k_decay(1:nom_pools), &
      ystates(1:nom_tot_elms), pot_decomp)

  pot_nn_flx = 0._r8; pot_np_flx = 0._r8

  reacs=(/lit1_dek_reac, lit2_dek_reac, lit3_dek_reac, &
    cwd_dek_reac, lwd_dek_reac, fwd_dek_reac, &
    som1_dek_reac, som2_dek_reac, som3_dek_reac/)

  do reac = 1, nom_pools
    if(alpha_n(reac)>0._r8)then
      pot_nn_flx = pot_nn_flx - cascade_matrix(lid_nh4, reacs(reac)) * pot_decomp(reac)
    endif
    if(alpha_p(reac)>0._r8)then
      pot_np_flx = pot_np_flx - cascade_matrix(lid_minp_soluble, reacs(reac)) * pot_decomp(reac)
    endif
  enddo
  end associate
  end subroutine calc_pot_min_np_flx

end module BgcCentSOMType
