module BgcCentCnpType
#include "bshr_assert.h"
  !
  ! !DESCRIPTION:
  ! subroutines for stoichiometric configuration of the century bgc
  ! !History, created by Jinyun Tang, Dec, 2014.
  ! !USES:
  use bshr_kind_mod       , only : r8 => shr_kind_r8
  use bshr_log_mod        , only : errMsg => shr_log_errMsg
  use betr_varcon         , only : spval => bspval, spinup_state => bspinup_state
  use BgcCentCnpDecompType      , only : DecompCent_type
  use BgcCentCnpIndexType , only : centurybgc_index_type
  use BgcCentCnpNitDenType   , only : century_nitden_type
  use gbetrType           , only : gbetr_type
  use BgcCentSOMType         , only : CentSom_type
  use BgcCentCnpCompetType       , only : Compet_ECA_type
  use BiogeoConType       , only : BiogeoCon_type
  use BetrStatusType      , only : betr_status_type
  implicit none
  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  !Note:
  !Keeping centurybgc_index as a private member is a workaround to call the ode solver
  !it increase the memory for each instance of centurybgceca_type, but enables
  !the ode function to be called by the ode solver

  type, extends(gbetr_type), public :: centurybgceca_type
    type(DecompCent_type),private        :: decompkf_eca
    type(century_nitden_type), private   :: nitden
    type(CentSom_type), private          :: censom
    type(Compet_ECA_type), public        :: competECA
    type(centurybgc_index_type), private :: centurybgc_index
    real(r8), pointer                    :: ystates0(:)
    real(r8), pointer                    :: ystates1(:)
    real(r8), pointer                    :: k_decay(:)
    real(r8), pointer                    :: cascade_matrix(:,:)
    real(r8), pointer                    :: cascade_matrixd(:,:)
    real(r8), pointer                    :: cascade_matrixp(:,:)
    real(r8), pointer                    :: alpha_n(:)
    real(r8), pointer                    :: alpha_p(:)
    real(r8)                             :: pot_f_nit
    real(r8)                             :: pot_f_denit
    real(r8)                             :: rt_ar
    real(r8)                             :: frac_p_sec_to_sol
    real(r8), pointer                    :: minp_secondary_decay(:)
    real(r8), pointer                    :: mumax_minp_soluble_to_secondary(:)
    integer                              :: plant_ntypes

    real(r8), pointer                    :: scal_f(:)
    real(r8), pointer                    :: conv_f(:)
    real(r8), pointer                    :: conc_f(:)
    integer                              :: soilorder
    real(r8)                             :: msurf_nh4
    real(r8)                             :: msurf_minp
    logical, private                     :: use_c13
    logical, private                     :: use_c14
  contains
    procedure, public  :: init
    procedure, public  :: runbgc
    procedure, private :: calc_cascade_matrix
    procedure, private :: init_states
    procedure, private :: add_ext_input
    procedure, private :: InitAllocate
    procedure, private :: arenchyma_gas_transport
    procedure, private :: sumup_cnp_mass
    procedure, private :: sumup_cnp_msflx
    procedure, private :: bgc_integrate
  end type centurybgceca_type
  logical, public :: ldebug_bgc =.false.

  public :: create_centuryeca_type
contains

  function create_centuryeca_type()
  ! DESCRIPTION
  ! constructor
    implicit none
    type(centurybgceca_type), pointer :: create_centuryeca_type
    type(centurybgceca_type), pointer :: bgc

    allocate(bgc)
    create_centuryeca_type => bgc

  end function create_centuryeca_type
  !-------------------------------------------------------------------------------
  subroutine init(this,  biogeo_con,  bstatus)
  use betr_varcon         , only : betr_maxpatch_pft
  implicit none
  class(centurybgceca_type), intent(inout) :: this
  type(BiogeoCon_type),intent(in) :: biogeo_con
  type(betr_status_type), intent(out) :: bstatus

  call bstatus%reset()
  call this%centurybgc_index%Init(biogeo_con%use_c13, biogeo_con%use_c14, betr_maxpatch_pft)

  call this%InitAllocate(this%centurybgc_index)

  call this%censom%Init(this%centurybgc_index, biogeo_con, bstatus)

  if(bstatus%check_status())return

  call this%decompkf_eca%Init(biogeo_con)

  call this%nitden%Init(biogeo_con)

  call this%competECA%Init()

  this%frac_p_sec_to_sol = biogeo_con%frac_p_sec_to_sol

  this%minp_secondary_decay = biogeo_con%minp_secondary_decay

  this%mumax_minp_soluble_to_secondary = biogeo_con%vmax_minp_soluble_to_secondary

  this%use_c13 = biogeo_con%use_c13

  this%use_c14 = biogeo_con%use_c14

  end subroutine init
  !-------------------------------------------------------------------------------

  subroutine InitAllocate(this, centurybgc_index)
  use BgcCentCnpIndexType , only : centurybgc_index_type
  use betr_varcon         , only : betr_maxpatch_pft, betr_max_soilorder
  implicit none
  class(centurybgceca_type)   , intent(inout) :: this
  type(centurybgc_index_type) , intent(in):: centurybgc_index

  associate(                                   &
    nom_pools=> centurybgc_index%nom_pools,    &
    nstvars => centurybgc_index%nstvars  ,     &
    nreactions => centurybgc_index%nreactions, &
    nprimvars => centurybgc_index%nprimvars    &
  )

  allocate(this%ystates0(nstvars)); this%ystates0(:) = 0._r8
  allocate(this%ystates1(nstvars)); this%ystates1(:) = 0._r8
  allocate(this%k_decay(nom_pools)); this%k_decay(:) = 0._r8
  allocate(this%scal_f(nprimvars));  this%scal_f(:) = 0._r8
  allocate(this%conv_f(nprimvars));  this%conv_f(:) = 0._r8
  allocate(this%conc_f(nprimvars));  this%conc_f(:) = 0._r8

  allocate(this%cascade_matrix(nstvars,  nreactions)); this%cascade_matrix(:,:) = 0._r8
  allocate(this%cascade_matrixd(1:nprimvars, 1:nreactions)); this%cascade_matrixd(:,:) = 0._r8
  allocate(this%cascade_matrixp(1:nprimvars, 1:nreactions)); this%cascade_matrixp(:,:) = 0._r8


  allocate(this%alpha_n(nom_pools)); this%alpha_n(:) = 0._r8
  allocate(this%alpha_p(nom_pools)); this%alpha_p(:) = 0._r8
  allocate(this%minp_secondary_decay(betr_max_soilorder)); this%minp_secondary_decay(:) = 0._r8
  allocate(this%mumax_minp_soluble_to_secondary(betr_max_soilorder)); this%mumax_minp_soluble_to_secondary(:) = 0._r8
  end associate
  end subroutine InitAllocate

  !-------------------------------------------------------------------------------
  subroutine runbgc(this,  is_surf, dtime, centuryeca_forc, nstates, ystates0, ystatesf, bstatus)

  !DESCRIPTION
  !do bgc model integration for one step
  use BgcCentCnpForcType        , only : centuryeca_forc_type
  use MathfuncMod               , only : pd_decomp
  use BetrStatusType            , only : betr_status_type
  use MathfuncMod               , only : safe_div
  use tracer_varcon             , only : catomw, natomw, patomw
  implicit none
  class(centurybgceca_type)  , intent(inout) :: this
  logical                    , intent(in) :: is_surf
  real(r8)                   , intent(in) :: dtime
  type(centuryeca_forc_type) , intent(in) :: centuryeca_forc
  integer                    , intent(in) :: nstates
  real(r8)                   , intent(out):: ystates0(nstates)
  real(r8)                   , intent(out):: ystatesf(nstates)
  type(betr_status_type)     , intent(out) :: bstatus


  !local variables
  real(r8) :: pot_om_decay_rates(this%centurybgc_index%nom_pools)
  real(r8) :: pot_co2_hr
  real(r8) :: pot_f_nit_mol_per_sec
  real(r8) :: n2_n2o_ratio_denit
  real(r8) :: yf(this%centurybgc_index%nstvars)
  real(r8) :: o2_decomp_depth
  real(r8) :: time
  real(r8) :: frc_c13, frc_c14
  real(r8) :: c_mass1, n_mass1, p_mass1
  real(r8) :: c_mass2, n_mass2, p_mass2
  real(r8) :: c_flx, n_flx, p_flx
  real(r8) :: c_inf, n_inf, p_inf
  integer :: jj
  character(len=*),parameter :: subname = 'runbgc'
  associate(                                            &
    pctsand   => centuryeca_forc%pct_sand             , &  !sand in %
    rt_ar     => centuryeca_forc%rt_ar                , &  !root autotrophic respiration
    rt_ar_c13 => centuryeca_forc%rt_ar_c13            , &  !root autotrophic respiration
    rt_ar_c14 => centuryeca_forc%rt_ar_c14            , &  !root autotrophic respiration
    lid_nh4   => this%centurybgc_index%lid_nh4        , &  !position id of nh4
    lid_no3   => this%centurybgc_index%lid_no3        , &  !
    lid_o2    => this%centurybgc_index%lid_o2         , &  !
    nom_pools => this%centurybgc_index%nom_pools      , &  !number om pools
    nom_tot_elms=> this%centurybgc_index%nom_tot_elms , &
    nstvars   => this%centurybgc_index%nstvars        , &
    nprimvars => this%centurybgc_index%nprimvars      , &
    nreactions=> this%centurybgc_index%nreactions     , &
    lid_plant_minn_nh4  => this%centurybgc_index%lid_plant_minn_nh4       , &
    lid_plant_minn_no3  => this%centurybgc_index%lid_plant_minn_no3       , &
    lid_n2o_nit => this%centurybgc_index%lid_n2o_nit,&
    lid_no3_den => this%centurybgc_index%lid_no3_den,&
    cascade_matrix=> this%cascade_matrix              , &
    cascade_matrixp=> this%cascade_matrixp            , &
    cascade_matrixd=> this%cascade_matrixd            , &
    ystates1 => this%ystates1                           &
  )
  !this%centurybgc_index%debug = centuryeca_forc%debug
  this%rt_ar = rt_ar
  frc_c13 = safe_div(rt_ar_c13,rt_ar); frc_c14 = safe_div(rt_ar_c14,rt_ar)
  call bstatus%reset()

  !initialize state variables
  call this%init_states(this%centurybgc_index, centuryeca_forc)
  ystates0(:) = this%ystates0(:)

  call this%sumup_cnp_msflx(ystates0, c_mass1,n_mass1,p_mass1)

  !add all external input
  call this%add_ext_input(dtime, this%centurybgc_index, centuryeca_forc, c_inf, n_inf, p_inf)
  c_inf = c_inf * catomw; n_inf=n_inf * natomw; p_inf=p_inf * patomw
!  call this%sumup_cnp_mass('afext input')

  !initialize decomposition scaling factors
  call this%decompkf_eca%set_decompk_scalar(ystates1(lid_o2), centuryeca_forc)

  !initialize all entries to zero
  cascade_matrix(:,:) = 0._r8

  !calculate default stoichiometry entries
  call this%calc_cascade_matrix(this%centurybgc_index, cascade_matrix, frc_c13, frc_c14)

  !run century decomposition, return decay rates, cascade matrix, potential hr
  call this%censom%run_decomp(is_surf, this%centurybgc_index, dtime, ystates1(1:nom_tot_elms),&
      this%decompkf_eca, centuryeca_forc%pct_sand, centuryeca_forc%pct_clay,this%alpha_n, this%alpha_p, &
      cascade_matrix, this%k_decay(1:nom_pools), pot_co2_hr, bstatus)

!  call this%sumup_cnp_mass('af run_decomp')
  if(bstatus%check_status())return

  call this%nitden%calc_pot_nitr(ystates1(lid_nh4), centuryeca_forc, this%decompkf_eca, pot_f_nit_mol_per_sec)

  !calculate potential o2 consumption
  o2_decomp_depth = pot_co2_hr + rt_ar + pot_f_nit_mol_per_sec * this%nitden%get_nit_o2_scef()

  !take a minimum > 0 to avoid singularity in calculating anaerobic fractions
  o2_decomp_depth = max(o2_decomp_depth,1.e-40_r8)

  !run nitrification-denitrification, returns cascade_matrix, decay rates
  call this%nitden%run_nitden(this%centurybgc_index, centuryeca_forc, this%decompkf_eca, &
    ystates1(lid_nh4), ystates1(lid_no3), ystates1(lid_o2), o2_decomp_depth, &
    pot_f_nit_mol_per_sec, pot_co2_hr, this%pot_f_nit, this%pot_f_denit, cascade_matrix)

  !---------------------
  !turn off nitrification and denitrification
  !this%pot_f_denit = 0._r8
  !this%pot_f_nit = 0._r8
  !---------------------
  !do integration, in each integration, the stoichiometric matrix is kept as constant
  !so the reaction rate is a function of state variables. Further, for simplicity,
  !the nitrification and denitrification rates have been assumed as linear function
  !nh4 and no3 in soil.

  call this%arenchyma_gas_transport(this%centurybgc_index, dtime)

  !do the stoichiometric matrix separation
  call pd_decomp(nprimvars, nreactions, cascade_matrix(1:nprimvars, 1:nreactions), &
     cascade_matrixp, cascade_matrixd, bstatus)
  if(bstatus%check_status())return

  time = 0._r8
  yf(:) = ystates1(:)

!  call this%sumup_cnp_mass('bfdecomp',c_mass1,n_mass1,p_mass1)

  call ode_adapt_ebbks1(this, yf, nprimvars, nstvars, time, dtime, ystates1)

!  call this%sumup_cnp_mass('afdecomp',c_mass2,n_mass2,p_mass2)

!  write(*,'(A,3(X,E20.10))')'cnp mass diff',c_mass2-c_mass1,n_mass2-n_mass1,p_mass2-p_mass1
  !print*,'sz',size(ystatesf),size(ystates1)

  ystatesf(:) = ystates1(:)

!  print*,'szx',maxval(ystatesf),maxval(ystates1)
  call this%sumup_cnp_msflx(ystatesf, c_mass2,n_mass2,p_mass2, c_flx, n_flx, p_flx)

  if(centuryeca_forc%debug)then
     write(*,'(A,10(1X,E20.10))')'cnp bal', &
     c_mass2*centuryeca_forc%dzsoi, &
     n_mass2*centuryeca_forc%dzsoi, &
     p_mass2*centuryeca_forc%dzsoi, &
     c_mass2 - c_mass1-c_inf+c_flx, &
     n_mass2 - n_mass1-n_inf+n_flx, &
     p_mass2 - p_mass1-p_inf+p_flx, &
     ystatesf(lid_plant_minn_nh4), ystatesf(lid_plant_minn_no3), &
     ystatesf(lid_n2o_nit), ystatesf(lid_no3_den)
  endif

  end associate
  end subroutine runbgc
  !-------------------------------------------------------------------------------
  subroutine calc_cascade_matrix(this,centurybgc_index, cascade_matrix, frc_c13, frc_c14)
    !
    ! !DESCRIPTION:
    ! calculate cascade matrix for the decomposition model
    !
    ! !USES:
    use MathfuncMod               , only : safe_div
    use BgcCentCnpIndexType       , only : centurybgc_index_type
    implicit none
    ! !ARGUMENTS:
    class(centurybgceca_type)     , intent(in) :: this
    type(centurybgc_index_type)   , intent(in) :: centurybgc_index
    real(r8)                      , intent(inout)   :: cascade_matrix(centurybgc_index%nstvars, centurybgc_index%nreactions)
    real(r8)                      , intent(in) :: frc_c13, frc_c14
    ! !LOCAL VARIABLES:
    real(r8) :: ftxt, f1, f2
    integer :: k, reac

    associate(                                                             & !
         lid_autr_rt => centurybgc_index%lid_autr_rt                      , & !
         lid_o2    => centurybgc_index%lid_o2                             , & !
         lid_co2   => centurybgc_index%lid_co2                            , & !
         lid_c13_co2=> centurybgc_index%lid_c13_co2                       , & !
         lid_c14_co2=> centurybgc_index%lid_c14_co2                       , & !
         lid_nh4   => centurybgc_index%lid_nh4                            , & !
         lid_ch4   => centurybgc_index%lid_ch4                            , & !
         lid_ar    => centurybgc_index%lid_ar                             , & !
         lid_no3   => centurybgc_index%lid_no3                            , & !
         lid_n2o   => centurybgc_index%lid_n2o                            , & !
         lid_n2    => centurybgc_index%lid_n2                             , & !
         lid_co2_hr=> centurybgc_index%lid_co2_hr                         , & !
         lid_minn_nh4_immob => centurybgc_index%lid_minn_nh4_immob        , & !
         lid_minn_no3_immob => centurybgc_index%lid_minn_no3_immob        , & !
         lid_n2_paere => centurybgc_index%lid_n2_paere                    , & !
         lid_ch4_paere => centurybgc_index%lid_ch4_paere                  , & !
         lid_n2o_paere => centurybgc_index%lid_n2o_paere                  , & !
         lid_o2_paere => centurybgc_index%lid_o2_paere                    , & !
         lid_ar_paere => centurybgc_index%lid_ar_paere                    , & !
         lid_co2_paere => centurybgc_index%lid_co2_paere                  , & !
         lid_c13_co2_paere => centurybgc_index%lid_c13_co2_paere          , & !
         lid_c14_co2_paere => centurybgc_index%lid_c14_co2_paere          , & !
         lid_minp_soluble => centurybgc_index%lid_minp_soluble            , & !
         lid_minp_secondary => centurybgc_index%lid_minp_secondary        , & !
         lid_minp_occlude =>  centurybgc_index%lid_minp_occlude           , & !
         lid_plant_minp => centurybgc_index%lid_plant_minp                , & !
         lid_minp_immob => centurybgc_index%lid_minp_immob                , & !

         lid_autr_rt_reac=> centurybgc_index%lid_autr_rt_reac                 , & !
         lid_no3_den  => centurybgc_index%lid_no3_den                     , & !
         lid_plant_minn_nh4_up_reac=> centurybgc_index%lid_plant_minn_nh4_up_reac , & !
         lid_plant_minn_no3_up_reac=> centurybgc_index%lid_plant_minn_no3_up_reac , & !
         lid_plant_minn_nh4  => centurybgc_index%lid_plant_minn_nh4       , &
         lid_plant_minn_no3  => centurybgc_index%lid_plant_minn_no3       , &
         lid_minp_secondary_to_sol_occ_reac => centurybgc_index%lid_minp_secondary_to_sol_occ_reac    , & !
         lid_minp_soluble_to_secp_reac => centurybgc_index%lid_minp_soluble_to_secp_reac      , & !
         lid_plant_minp_up_reac => centurybgc_index%lid_plant_minp_up_reac, & !

         lid_n2_aren_reac => centurybgc_index%lid_n2_aren_reac            , & !
         lid_ch4_aren_reac=> centurybgc_index%lid_ch4_aren_reac           , & !
         lid_n2o_aren_reac=> centurybgc_index%lid_n2o_aren_reac           , & !
         lid_o2_aren_reac => centurybgc_index%lid_o2_aren_reac            , & !
         lid_ar_aren_reac => centurybgc_index%lid_ar_aren_reac            , & !
         lid_co2_aren_reac=> centurybgc_index%lid_co2_aren_reac           , & !
         lid_c13_co2_aren_reac=> centurybgc_index%lid_c13_co2_aren_reac   , & !
         lid_c14_co2_aren_reac=> centurybgc_index%lid_c14_co2_aren_reac     & !
         )

    !higher [nh4] makes lower [no3] competitiveness
    !note all reactions are in the form products - substrates = 0, therefore
    !mass balance is automatically ensured.
    !set up first order reactions

    !---------------------------------------------------------------------------------
    !reaction 10, inorganic P non-equilibrium adsorption
    !P_soluble -> p_secondary
    reac = lid_minp_soluble_to_secp_reac
    cascade_matrix(lid_minp_soluble,  reac) = -1._r8
    cascade_matrix(lid_minp_secondary, reac) = 1._r8

    !----------------------------------------------------------------------
    !reaction 11, inorganic P non-equilibrium desorption
    ! p_secondary -> P_soluble + P_occlude
    reac = lid_minp_secondary_to_sol_occ_reac
    cascade_matrix(lid_minp_soluble,  reac) = this%frac_p_sec_to_sol
    cascade_matrix(lid_minp_occlude  ,  reac) = 1._r8 - this%frac_p_sec_to_sol
    cascade_matrix(lid_minp_secondary, reac) = -1._r8

    !----------------------------------------------------------------------
    !reaction 12, plant mineral nitrogen nh4 uptake
    reac = lid_plant_minn_nh4_up_reac
    !  nh4  -> plant_nitrogen
    cascade_matrix(lid_nh4, reac)        = -1._r8
    cascade_matrix(lid_plant_minn_nh4, reac) = 1._r8

    !----------------------------------------------------------------------
    !reaction 13, plant mineral nitrogen no3 uptake
    reac = lid_plant_minn_no3_up_reac
    !  no3  -> plant_nitrogen
    cascade_matrix(lid_no3, reac)        = -1._r8
    cascade_matrix(lid_plant_minn_no3, reac) = 1._r8

    !----------------------------------------------------------------------
    !reaction 14, plant mineral phosphorus uptake
    reac = lid_plant_minp_up_reac
    ! p_solution -> plant_p
    cascade_matrix(lid_minp_soluble, reac) = -1._r8
    cascade_matrix(lid_plant_minp, reac) = 1._r8

    !----------------------------------------------------------------------
    !reaction 15, ar + o2 -> co2
    reac = lid_autr_rt_reac
    cascade_matrix(lid_co2, reac) =  1._r8
    cascade_matrix(lid_o2,  reac) = -1._r8

    if(this%use_c13)then
      cascade_matrix(lid_c13_co2, reac) =  1._r8 * frc_c13
    endif

    if(this%use_c14)then
      cascade_matrix(lid_c14_co2, reac) =  1._r8 * frc_c14
    endif
    !--------------------------------------------------------------------
    !arenchyma transport
    !second primary variables
    reac                               = lid_o2_aren_reac
    cascade_matrix(lid_o2, reac)       = -1._r8
    cascade_matrix(lid_o2_paere, reac) = 1._r8

    if ( spinup_state /= 1 ) then
       reac                                = lid_ch4_aren_reac
       cascade_matrix(lid_ch4, reac)       = -1._r8
       cascade_matrix(lid_ch4_paere, reac) = 1._r8

       reac                                = lid_ar_aren_reac
       cascade_matrix(lid_ar, reac)        = -1._r8
       cascade_matrix(lid_ar_paere, reac)  = 1._r8

       reac                                = lid_co2_aren_reac
       cascade_matrix(lid_co2, reac)       = -1._r8
       cascade_matrix(lid_co2_paere, reac) = 1._r8

       if(this%use_c13)then
         reac                                = lid_c13_co2_aren_reac
         cascade_matrix(lid_c13_co2, reac)       = -1._r8
         cascade_matrix(lid_c13_co2_paere, reac) = 1._r8
       endif

       if(this%use_c14)then
         reac                                = lid_c14_co2_aren_reac
         cascade_matrix(lid_c14_co2, reac)       = -1._r8
         cascade_matrix(lid_c14_co2_paere, reac) = 1._r8
       endif

       reac                                = lid_n2o_aren_reac
       cascade_matrix(lid_n2o, reac)       = -1._r8
       cascade_matrix(lid_n2o_paere, reac) = 1._r8

       reac                                = lid_n2_aren_reac
       cascade_matrix(lid_n2, reac)        = -1._r8
       cascade_matrix(lid_n2_paere, reac)  = 1._r8
    endif

  end associate
  end subroutine calc_cascade_matrix
  !--------------------------------------------------------------------
  subroutine init_states(this, centurybgc_index, centuryeca_forc)

  use BgcCentCnpIndexType       , only : centurybgc_index_type
  use BgcCentCnpForcType        , only : centuryeca_forc_type
  implicit none
  class(centurybgceca_type)     , intent(inout) :: this
  type(centurybgc_index_type)  , intent(in) :: centurybgc_index
  type(centuryeca_forc_type)  , intent(in) :: centuryeca_forc

  integer :: j
  associate(                           &
    lid_n2 => centurybgc_index%lid_n2, &
    lid_o2 => centurybgc_index%lid_o2, &
    lid_co2 => centurybgc_index%lid_co2, &
    lid_c13_co2 => centurybgc_index%lid_c13_co2, &
    lid_c14_co2 => centurybgc_index%lid_c14_co2, &
    lid_n2o => centurybgc_index%lid_n2o, &
    lid_ar => centurybgc_index%lid_ar, &
    lid_ch4 => centurybgc_index%lid_ch4  &
  )
  this%ystates0(:) = centuryeca_forc%ystates(:)
  this%ystates1(:) = this%ystates0(:)

  !set conversion parameters for arenchyma transport
  this%scal_f(lid_n2) = centuryeca_forc%aren_cond_n2
  this%conc_f(lid_n2) = centuryeca_forc%conc_atm_n2
  this%conv_f(lid_n2) = 1._r8/centuryeca_forc%n2_g2b

  this%scal_f(lid_o2) = centuryeca_forc%aren_cond_o2
  this%conc_f(lid_o2) = centuryeca_forc%conc_atm_o2
  this%conv_f(lid_o2) = 1._r8/centuryeca_forc%o2_g2b

  this%scal_f(lid_ar) = centuryeca_forc%aren_cond_ar
  this%conc_f(lid_ar) = centuryeca_forc%conc_atm_ar
  this%conv_f(lid_ar) = 1._r8/centuryeca_forc%ar_g2b

  this%scal_f(lid_co2) = centuryeca_forc%aren_cond_co2
  this%conc_f(lid_co2) = centuryeca_forc%conc_atm_co2
  this%conv_f(lid_co2) = 1._r8/centuryeca_forc%co2_g2b

  if(this%use_c13)then
    this%scal_f(lid_c13_co2) = centuryeca_forc%aren_cond_co2_c13
    this%conc_f(lid_c13_co2) = centuryeca_forc%conc_atm_co2_c13
    this%conv_f(lid_c13_co2) = 1._r8/centuryeca_forc%co2_g2b
  endif

  if(this%use_c14)then
    this%scal_f(lid_c14_co2) = centuryeca_forc%aren_cond_co2_c14
    this%conc_f(lid_c14_co2) = centuryeca_forc%conc_atm_co2_c14
  endif

  this%scal_f(lid_ch4) = centuryeca_forc%aren_cond_ch4
  this%conc_f(lid_ch4) = centuryeca_forc%conc_atm_ch4
  this%conv_f(lid_ch4) = 1._r8/centuryeca_forc%ch4_g2b

  this%scal_f(lid_n2o) = centuryeca_forc%aren_cond_n2o
  this%conc_f(lid_n2o) = centuryeca_forc%conc_atm_n2o
  this%conv_f(lid_n2o) = 1._r8/centuryeca_forc%n2o_g2b

  this%plant_ntypes = centuryeca_forc%plant_ntypes
  this%soilorder = centuryeca_forc%soilorder

  this%msurf_nh4 = centuryeca_forc%msurf_nh4
  this%msurf_minp = centuryeca_forc%msurf_minp

  end associate
  end subroutine init_states
  !--------------------------------------------------------------------
  subroutine add_ext_input(this, dtime, centurybgc_index, centuryeca_forc, c_inf, n_inf, p_inf)
  use BgcCentCnpIndexType       , only : centurybgc_index_type
  use BgcCentCnpForcType        , only : centuryeca_forc_type
  use tracer_varcon             , only : catomw, natomw, patomw,c13atomw,c14atomw
  implicit none
  class(centurybgceca_type)     , intent(inout) :: this
  real(r8), intent(in) :: dtime
  type(centurybgc_index_type)  , intent(in) :: centurybgc_index
  type(centuryeca_forc_type)  , intent(in) :: centuryeca_forc
  real(r8), optional, intent(out) :: c_inf, n_inf, p_inf
  integer :: kc, kn, kp,kc13,kc14
  integer :: jj

  associate(                        &
    lit1 =>  centurybgc_index%lit1, &
    lit2 =>  centurybgc_index%lit2, &
    lit3 =>  centurybgc_index%lit3, &
    cwd =>   centurybgc_index%cwd, &
    fwd =>   centurybgc_index%fwd, &
    lwd =>   centurybgc_index%lwd, &
    c_loc=>  centurybgc_index%c_loc,&
    c13_loc=>  centurybgc_index%c13_loc,&
    c14_loc=>  centurybgc_index%c14_loc,&
    n_loc=>  centurybgc_index%n_loc,&
    p_loc=>  centurybgc_index%p_loc,&
    som1 =>  centurybgc_index%som1, &
    som2 =>  centurybgc_index%som2, &
    som3 =>  centurybgc_index%som3, &
    nelms => centurybgc_index%nelms, &
    lid_nh4=> centurybgc_index%lid_nh4, &
    lid_minp_soluble =>  centurybgc_index%lid_minp_soluble  &
  )
  jj=lit1;kc = (jj-1)*nelms+c_loc;kn=(jj-1)*nelms+n_loc;kp=(jj-1)*nelms+p_loc
  this%ystates1(kc) =this%ystates0(kc) + centuryeca_forc%cflx_input_litr_met*dtime/catomw
  this%ystates1(kn) =this%ystates0(kn) + centuryeca_forc%nflx_input_litr_met*dtime/natomw
  this%ystates1(kp) =this%ystates0(kp) + centuryeca_forc%pflx_input_litr_met*dtime/patomw

  if(present(c_inf))then
    c_inf = this%ystates1(kc) - this%ystates0(kc)
  endif
  if(present(n_inf))then
    n_inf = this%ystates1(kn) - this%ystates0(kn)
  endif
  if(present(p_inf))then
    p_inf = this%ystates1(kp) - this%ystates0(kp)
  endif

  if(this%use_c13)then
    kc13=(jj-1)*nelms+c13_loc
    this%ystates1(kc13) =this%ystates0(kc13) + centuryeca_forc%cflx_input_litr_met_c13*dtime/c13atomw
  endif
  if(this%use_c14)then
    kc14=(jj-1)*nelms+c14_loc
    this%ystates1(kc14) =this%ystates0(kc14) + centuryeca_forc%cflx_input_litr_met_c14*dtime/c14atomw
  endif

  jj=lit2;kc = (jj-1)*nelms+c_loc;kn=(jj-1)*nelms+n_loc;kp=(jj-1)*nelms+p_loc
  this%ystates1(kc) =this%ystates0(kc) + centuryeca_forc%cflx_input_litr_cel*dtime/catomw
  this%ystates1(kn) =this%ystates0(kn) + centuryeca_forc%nflx_input_litr_cel*dtime/natomw
  this%ystates1(kp) =this%ystates0(kp) + centuryeca_forc%pflx_input_litr_cel*dtime/patomw

  if(present(c_inf))then
    c_inf = c_inf + this%ystates1(kc) - this%ystates0(kc)
  endif
  if(present(n_inf))then
    n_inf = n_inf + this%ystates1(kn) - this%ystates0(kn)
  endif
  if(present(p_inf))then
    p_inf = p_inf + this%ystates1(kp) - this%ystates0(kp)
  endif
  if(this%use_c13)then
    kc13=(jj-1)*nelms+c13_loc
    this%ystates1(kc13) =this%ystates0(kc13) + centuryeca_forc%cflx_input_litr_cel_c13*dtime/c13atomw
  endif
  if(this%use_c14)then
    kc14=(jj-1)*nelms+c14_loc
    this%ystates1(kc14) =this%ystates0(kc14) + centuryeca_forc%cflx_input_litr_cel_c14*dtime/c14atomw
  endif

  jj=lit3;kc = (jj-1)*nelms+c_loc;kn=(jj-1)*nelms+n_loc;kp=(jj-1)*nelms+p_loc
  this%ystates1(kc) =this%ystates0(kc) + centuryeca_forc%cflx_input_litr_lig*dtime/catomw
  this%ystates1(kn) =this%ystates0(kn) + centuryeca_forc%nflx_input_litr_lig*dtime/natomw
  this%ystates1(kp) =this%ystates0(kp) + centuryeca_forc%pflx_input_litr_lig*dtime/patomw

  if(present(c_inf))then
    c_inf = c_inf + this%ystates1(kc) - this%ystates0(kc)
  endif
  if(present(n_inf))then
    n_inf = n_inf + this%ystates1(kn) - this%ystates0(kn)
  endif
  if(present(p_inf))then
    p_inf = p_inf + this%ystates1(kp) - this%ystates0(kp)
  endif
  if(this%use_c13)then
    kc13=(jj-1)*nelms+c13_loc
    this%ystates1(kc13) =this%ystates0(kc13) + centuryeca_forc%cflx_input_litr_lig_c13*dtime/c13atomw
  endif
  if(this%use_c14)then
    kc14=(jj-1)*nelms+c14_loc
    this%ystates1(kc14) =this%ystates0(kc14) + centuryeca_forc%cflx_input_litr_lig_c14*dtime/c14atomw
  endif

  jj=cwd;kc = (jj-1)*nelms+c_loc;kn=(jj-1)*nelms+n_loc;kp=(jj-1)*nelms+p_loc
  this%ystates1(kc) =this%ystates0(kc) + centuryeca_forc%cflx_input_litr_cwd*dtime/catomw
  this%ystates1(kn) =this%ystates0(kn) + centuryeca_forc%nflx_input_litr_cwd*dtime/natomw
  this%ystates1(kp) =this%ystates0(kp) + centuryeca_forc%pflx_input_litr_cwd*dtime/patomw

  if(present(c_inf))then
    c_inf = c_inf + this%ystates1(kc) - this%ystates0(kc)
  endif
  if(present(n_inf))then
    n_inf = n_inf + this%ystates1(kn) - this%ystates0(kn)
  endif
  if(present(p_inf))then
    p_inf = p_inf + this%ystates1(kp) - this%ystates0(kp)
  endif
  if(this%use_c13)then
    kc13=(jj-1)*nelms+c13_loc
    this%ystates1(kc13) =this%ystates0(kc13) + centuryeca_forc%cflx_input_litr_cwd_c13*dtime/c13atomw
  endif
  if(this%use_c14)then
    kc14=(jj-1)*nelms+c14_loc
    this%ystates1(kc14) =this%ystates0(kc14) + centuryeca_forc%cflx_input_litr_cwd_c14*dtime/c14atomw
  endif

  jj=fwd;kc = (jj-1)*nelms+c_loc;kn=(jj-1)*nelms+n_loc;kp=(jj-1)*nelms+p_loc
  this%ystates1(kc) =this%ystates0(kc) + centuryeca_forc%cflx_input_litr_fwd*dtime/catomw
  this%ystates1(kn) =this%ystates0(kn) + centuryeca_forc%nflx_input_litr_fwd*dtime/natomw
  this%ystates1(kp) =this%ystates0(kp) + centuryeca_forc%pflx_input_litr_fwd*dtime/patomw

  if(present(c_inf))then
    c_inf = c_inf + this%ystates1(kc) - this%ystates0(kc)
  endif
  if(present(n_inf))then
    n_inf = n_inf + this%ystates1(kn) - this%ystates0(kn)
  endif
  if(present(p_inf))then
    p_inf = p_inf + this%ystates1(kp) - this%ystates0(kp)
  endif
  if(this%use_c13)then
    kc13=(jj-1)*nelms+c13_loc
    this%ystates1(kc13) =this%ystates0(kc13) + centuryeca_forc%cflx_input_litr_fwd_c13*dtime/c13atomw
  endif
  if(this%use_c14)then
    kc14=(jj-1)*nelms+c14_loc
    this%ystates1(kc14) =this%ystates0(kc14) + centuryeca_forc%cflx_input_litr_fwd_c14*dtime/c14atomw
  endif

  jj=lwd;kc = (jj-1)*nelms+c_loc;kn=(jj-1)*nelms+n_loc;kp=(jj-1)*nelms+p_loc
  this%ystates1(kc) =this%ystates0(kc) + centuryeca_forc%cflx_input_litr_lwd*dtime/catomw
  this%ystates1(kn) =this%ystates0(kn) + centuryeca_forc%nflx_input_litr_lwd*dtime/natomw
  this%ystates1(kp) =this%ystates0(kp) + centuryeca_forc%pflx_input_litr_lwd*dtime/patomw

  if(present(c_inf))then
    c_inf = c_inf + this%ystates1(kc) - this%ystates0(kc)
  endif
  if(present(n_inf))then
    n_inf = n_inf + this%ystates1(kn) - this%ystates0(kn)
  endif
  if(present(p_inf))then
    p_inf = p_inf + this%ystates1(kp) - this%ystates0(kp)
  endif
  if(this%use_c13)then
    kc13=(jj-1)*nelms+c13_loc
    this%ystates1(kc13) =this%ystates0(kc13) + centuryeca_forc%cflx_input_litr_lwd_c13*dtime/c13atomw
  endif
  if(this%use_c14)then
    kc14=(jj-1)*nelms+c14_loc
    this%ystates1(kc14) =this%ystates0(kc14) + centuryeca_forc%cflx_input_litr_lwd_c14*dtime/c14atomw
  endif

  this%ystates1(lid_nh4) =this%ystates0(lid_nh4) + dtime * &
      (centuryeca_forc%sflx_minn_input_nh4 + &
        centuryeca_forc%sflx_minn_nh4_fix_nomic)/natomw

  this%ystates1(lid_minp_soluble) =this%ystates0(lid_minp_soluble) + dtime * &
      (centuryeca_forc%sflx_minp_input_po4 + &
        centuryeca_forc%sflx_minp_weathering_po4)/patomw

  if(present(n_inf))then
    n_inf = n_inf + this%ystates1(lid_nh4) - this%ystates0(lid_nh4)
  endif

  if(present(p_inf))then
    p_inf = p_inf + this%ystates1(lid_minp_soluble) - this%ystates0(lid_minp_soluble)
  endif
  end associate
  end subroutine add_ext_input


  !--------------------------------------------------------------------
  subroutine bgc_integrate(this, ystate, dtime, time, nprimvars, nstvars, dydt)
  !
  !DESCRIPTION
  !calculate the reaction rates
  !In current implementation, no active adsorption of NH4 is involved.
  !The inorganic phosphorus does the transition from soluble->secondary->occlude
  use SOMStateVarUpdateMod , only : calc_dtrend_som_bgc
  use MathfuncMod          , only : lom_type, safe_div
  implicit none
  class(centurybgceca_type) , intent(inout) :: this
  integer                   , intent(in) :: nstvars
  integer                   , intent(in) :: nprimvars
  real(r8)                  , intent(in) :: dtime
  real(r8)                  , intent(in) :: time
  real(r8)                  , intent(in) :: ystate(nstvars)
  real(r8)                  , intent(out) :: dydt(nstvars)

  !local variables
  real(r8) :: mic_pot_nn_flx  !potential nitrogen uptake to support decomposition
  real(r8) :: mic_pot_np_flx  !potential phosphorus uptake to support decomposition
  real(r8) :: pot_decomp(this%centurybgc_index%nom_pools)
  real(r8) :: rrates(this%centurybgc_index%nreactions)
  real(r8) :: p_dt(1:this%centurybgc_index%nprimvars)
  real(r8) :: d_dt(1:this%centurybgc_index%nprimvars)
  real(r8) :: dydt1(nstvars)
  real(r8) :: pscal(1:nprimvars)
  real(r8) :: rscal(1:this%centurybgc_index%nreactions)
  real(r8) :: ECA_flx_nh4_plants(this%plant_ntypes)
  real(r8) :: ECA_flx_no3_plants(this%plant_ntypes)
  real(r8) :: ECA_factor_msurf_nh4
  real(r8) :: ECA_flx_phosphorus_plants(this%plant_ntypes)
  real(r8) :: ECA_factor_minp_msurf
  real(r8) :: ECA_factor_phosphorus_mic
  real(r8) :: ECA_factor_nh4_mic
  real(r8) :: ECA_factor_no3_mic
  real(r8) :: ECA_factor_nitrogen_mic
  real(r8) :: ECA_factor_den
  real(r8) :: ECA_factor_nit
  integer  :: jj, it
  integer, parameter  :: itmax = 100
  type(lom_type) :: lom
  type(betr_status_type) :: bstatus
  logical :: lneg
  real(r8) :: scal

  associate(                                                                                      &
    nreactions => this%centurybgc_index%nreactions                                              , &
    nstvars  => this%centurybgc_index%nstvars                                                   , &
    nom_pools => this%centurybgc_index%nom_pools                                                , &
    lid_nh4 => this%centurybgc_index%lid_nh4                                                    , &
    lid_no3 => this%centurybgc_index%lid_no3                                                    , &
    lid_plant_minn_no3_pft=> this%centurybgc_index%lid_plant_minn_no3_pft                       , &
    lid_plant_minn_nh4_pft=> this%centurybgc_index%lid_plant_minn_nh4_pft                       , &
    lid_plant_minp_pft=> this%centurybgc_index%lid_plant_minp_pft                               , &
    lid_plant_minp    => this%centurybgc_index%lid_plant_minp                                   , &
    lid_plant_minn_nh4 => this%centurybgc_index%lid_plant_minn_nh4                              , &
    lid_plant_minn_no3 => this%centurybgc_index%lid_plant_minn_no3                              , &
    lid_minp_soluble => this%centurybgc_index%lid_minp_soluble                                  , &
    lid_minp_secondary => this%centurybgc_index%lid_minp_secondary                              , &
    lid_minp_soluble_to_secp_reac=> this%centurybgc_index%lid_minp_soluble_to_secp_reac         , &
    lid_autr_rt_reac => this%centurybgc_index%lid_autr_rt_reac                                  , &
    lid_nh4_nit_reac => this%centurybgc_index%lid_nh4_nit_reac                                  , &
    lid_no3_den_reac => this%centurybgc_index%lid_no3_den_reac                                  , &
    lid_plant_minn_nh4_up_reac => this%centurybgc_index%lid_plant_minn_nh4_up_reac              , &
    lid_plant_minn_no3_up_reac => this%centurybgc_index%lid_plant_minn_no3_up_reac              , &
    lid_plant_minp_up_reac => this%centurybgc_index%lid_plant_minp_up_reac                      , &
    lid_minp_secondary_to_sol_occ_reac=> this%centurybgc_index%lid_minp_secondary_to_sol_occ_reac &
  )

  dydt(:) = 0._r8
  rrates(:) = 0._r8
  !calculate reaction rates, because arenchyma transport is
  !done, now only calculate for those that are actively
  !reacting. These include: OM pools, plant nutrient uptake
  !microbial nutrient uptake

  call this%censom%calc_pot_min_np_flx(dtime, this%centurybgc_index,  ystate, this%k_decay,&
    this%cascade_matrix, this%alpha_n, this%alpha_p, pot_decomp, mic_pot_nn_flx, mic_pot_np_flx)

  !do ECA nutrient scaling
  !
  call this%competECA%run_compet_nitrogen(ystate(lid_nh4),ystate(lid_no3),mic_pot_nn_flx,&
     this%pot_f_nit, this%pot_f_denit, this%plant_ntypes, &
     this%msurf_nh4, ECA_factor_nit, &
     ECA_factor_den, ECA_factor_nh4_mic, ECA_factor_no3_mic, &
     ECA_flx_nh4_plants,ECA_flx_no3_plants, ECA_factor_msurf_nh4)

  ECA_factor_nitrogen_mic = ECA_factor_nh4_mic + ECA_factor_no3_mic
  call this%competECA%run_compet_phosphorus(ystate(lid_minp_soluble), mic_pot_np_flx, &
     this%plant_ntypes, this%msurf_minp, &
     ECA_factor_phosphorus_mic, ECA_factor_minp_msurf, ECA_flx_phosphorus_plants)

  !apply ECA factor to obtain actual reaction rate, decomposition
  !plant, nit, den nutrient uptake,
  do jj = 1, nom_pools
    scal = 1._r8
    if(this%alpha_n(jj)>0._r8)then
      scal = min(scal, ECA_factor_nitrogen_mic)
      this%cascade_matrixd(lid_no3,jj) = this%cascade_matrix(lid_nh4,jj) * &
          safe_div(ECA_factor_no3_mic,ECA_factor_nitrogen_mic)
      this%cascade_matrixd(lid_nh4,jj) = this%cascade_matrix(lid_nh4,jj)-this%cascade_matrixd(lid_no3,jj)
    endif
    if(this%alpha_p(jj)>0._r8)scal = min(scal, ECA_factor_phosphorus_mic)
    if(scal /= 1._r8)pot_decomp(jj)=pot_decomp(jj)*scal
    rrates(jj) = pot_decomp(jj)
  enddo

  rrates(lid_nh4_nit_reac) = this%pot_f_nit*ECA_factor_nit
  rrates(lid_no3_den_reac) = this%pot_f_denit*ECA_factor_den
  rrates(lid_minp_soluble_to_secp_reac) =  ECA_factor_minp_msurf * this%msurf_minp &
       * this%mumax_minp_soluble_to_secondary(this%soilorder) !calculate from eca competition
  rrates(lid_autr_rt_reac) = this%rt_ar                            !authotrophic respiration
  rrates(lid_plant_minn_no3_up_reac) = sum(ECA_flx_no3_plants)     !calculate by ECA competition
  rrates(lid_plant_minn_nh4_up_reac) = sum(ECA_flx_nh4_plants)     !calculate by ECA competition
  rrates(lid_plant_minp_up_reac) =     sum(ECA_flx_phosphorus_plants) !calculate by ECA competition
  rrates(lid_minp_secondary_to_sol_occ_reac)= ystate(lid_minp_secondary) * this%minp_secondary_decay(this%soilorder)

  !the following treatment is to ensure mass balance
  if(this%plant_ntypes==1)then
    do jj = 1, this%plant_ntypes
      this%cascade_matrixd(lid_plant_minn_no3_pft(jj),lid_plant_minn_no3_up_reac) = 1._r8
      this%cascade_matrixd(lid_plant_minn_nh4_pft(jj),lid_plant_minn_nh4_up_reac) = 1._r8
      this%cascade_matrixd(lid_plant_minp_pft(jj),lid_plant_minp_up_reac) = 1._r8
    enddo
  elseif(this%plant_ntypes>=2)then
    do jj = 1, this%plant_ntypes-1
      this%cascade_matrixd(lid_plant_minn_no3_pft(jj),lid_plant_minn_no3_up_reac) = &
           safe_div(ECA_flx_no3_plants(jj),rrates(lid_plant_minn_no3_up_reac))
      this%cascade_matrixd(lid_plant_minn_nh4_pft(jj),lid_plant_minn_nh4_up_reac) = &
           safe_div(ECA_flx_nh4_plants(jj),rrates(lid_plant_minn_nh4_up_reac))
      this%cascade_matrixd(lid_plant_minp_pft(jj),lid_plant_minp_up_reac) = &
           safe_div(ECA_flx_phosphorus_plants(jj),rrates(lid_plant_minp_up_reac))
    enddo
    jj = this%plant_ntypes
    this%cascade_matrixd(lid_plant_minn_no3_pft(jj),lid_plant_minn_no3_up_reac) = &
      1._r8 - sum(this%cascade_matrixd(lid_plant_minn_no3_pft(1:jj-1),lid_plant_minn_no3_up_reac))
    this%cascade_matrixd(lid_plant_minn_nh4_pft(jj),lid_plant_minn_nh4_up_reac) = &
      1._r8 - sum(this%cascade_matrixd(lid_plant_minn_nh4_pft(1:jj-1),lid_plant_minn_nh4_up_reac))
    this%cascade_matrixd(lid_plant_minp_pft(jj),lid_plant_minp_up_reac) = &
      1._r8 - sum(this%cascade_matrixd(lid_plant_minp_pft(1:jj-1),lid_plant_minp_up_reac))
  endif

  it=0
  do
    call calc_dtrend_som_bgc(nprimvars, nreactions, this%cascade_matrixp(1:nprimvars, 1:nreactions), rrates, p_dt)

    call calc_dtrend_som_bgc(nprimvars, nreactions, this%cascade_matrixd(1:nprimvars, 1:nreactions), rrates, d_dt)

    !update the state variables
    call lom%calc_state_pscal(nprimvars, dtime, ystate(1:nprimvars), p_dt(1:nprimvars),  d_dt(1:nprimvars), &
        pscal(1:nprimvars), lneg, bstatus)

    if(lneg .and. it<=itmax)then
      call lom%calc_reaction_rscal(nprimvars, nreactions,  pscal(1:nprimvars), &
        this%cascade_matrixd(1:nprimvars, 1:nreactions),rscal, bstatus)

      call lom%apply_reaction_rscal(nreactions, rscal(1:nreactions), rrates(1:nreactions))
    else
      call calc_dtrend_som_bgc(nstvars, nreactions, this%cascade_matrix(1:nstvars, 1:nreactions), &
         rrates(1:nreactions), dydt)
      exit
    endif
    it = it + 1
  enddo

  end associate
  end subroutine bgc_integrate
  !--------------------------------------------------------------------
  subroutine arenchyma_gas_transport(this, centurybgc_index, dtime)
  use BgcCentCnpIndexType       , only : centurybgc_index_type
  implicit none
  class(centurybgceca_type)     , intent(inout) :: this
  type(centurybgc_index_type)  , intent(in) :: centurybgc_index
  real(r8), intent(in) :: dtime

  integer :: j
  real(r8) :: y0
  associate(                             &
    lid_n2 => centurybgc_index%lid_n2,   &
    lid_o2 => centurybgc_index%lid_o2,   &
    lid_co2 => centurybgc_index%lid_co2, &
    lid_c13_co2 => centurybgc_index%lid_c13_co2, &
    lid_c14_co2 => centurybgc_index%lid_c14_co2, &
    lid_n2o => centurybgc_index%lid_n2o, &
    lid_ar => centurybgc_index%lid_ar,   &
    lid_ch4 => centurybgc_index%lid_ch4  &
  )
  j = lid_n2; y0=this%ystates1(j)
  call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
  this%ystates1(centurybgc_index%lid_n2_paere) = this%ystates1(j)-y0

  j = lid_o2; y0=this%ystates1(j)
  call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
  this%ystates1(centurybgc_index%lid_o2_paere) = this%ystates1(j)-y0

  j = lid_ar; y0=this%ystates1(j)
  call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
  this%ystates1(centurybgc_index%lid_ar_paere) = this%ystates1(j)-y0

  j = lid_ch4; y0=this%ystates1(j)
  call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
  this%ystates1(centurybgc_index%lid_ch4_paere) = this%ystates1(j)-y0

  j = lid_co2; y0=this%ystates1(j)
  call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
  this%ystates1(centurybgc_index%lid_co2_paere) = this%ystates1(j)-y0

  if(this%use_c13)then
    j = lid_c13_co2; y0=this%ystates1(j)
    call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
    this%ystates1(centurybgc_index%lid_c13_co2_paere) = this%ystates1(j)-y0

  endif

  if(this%use_c14)then
    j = lid_c14_co2; y0=this%ystates1(j)
    call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
    this%ystates1(centurybgc_index%lid_c14_co2_paere) = this%ystates1(j)-y0

  endif
  j = lid_n2o; y0=this%ystates1(j)
  call exp_ode_int(dtime, this%scal_f(j), this%conv_f(j), this%conc_f(j), this%ystates1(j))
  this%ystates1(centurybgc_index%lid_n2o_paere) = this%ystates1(j)-y0

  end associate
  contains
    subroutine exp_ode_int(dt, c1, c2, c3, y0)
    !
    ! DESCRIPTION
    ! solve dy/dt=-c1*(c2*y-c3) using analytic solution
    implicit none
    real(r8), intent(in) :: dt
    real(r8), intent(in) :: c1
    real(r8), intent(in) :: c2
    real(r8), intent(in) :: c3
    real(r8), intent(inout) :: y0

    if(c1>0._r8)then
      y0 = c3/c2+(y0-c3/c2)*exp(-c1/c2*dtime)
    endif

    end subroutine exp_ode_int
  end subroutine arenchyma_gas_transport

  !--------------------------------------------------------------------
  subroutine sumup_cnp_msflx(this, ystates1, c_mass, n_mass, p_mass,c_flx,n_flx,p_flx)
  use tracer_varcon, only : catomw, natomw, patomw
  implicit none
  class(centurybgceca_type)     , intent(in) :: this
  real(r8), intent(in)  :: ystates1(this%centurybgc_index%nstvars)
  real(r8), intent(out) :: c_mass, n_mass, p_mass
  real(r8), optional, intent(out) :: c_flx, n_flx, p_flx
  !local variables

  integer :: kc, kn, kp, jj
  associate(                        &
    c_loc=>  this%centurybgc_index%c_loc,&
    n_loc=>  this%centurybgc_index%n_loc,&
    p_loc=>  this%centurybgc_index%p_loc,&
    lid_n2o_nit => this%centurybgc_index%lid_n2o_nit,&
    lid_no3_den => this%centurybgc_index%lid_no3_den,&
    lit1 =>  this%centurybgc_index%lit1, &
    lit2 =>  this%centurybgc_index%lit2, &
    lit3 =>  this%centurybgc_index%lit3, &
    cwd =>   this%centurybgc_index%cwd, &
    lwd =>   this%centurybgc_index%lwd, &
    fwd =>   this%centurybgc_index%fwd, &
    som1 =>  this%centurybgc_index%som1, &
    som2 =>  this%centurybgc_index%som2, &
    som3 =>  this%centurybgc_index%som3, &
    nelms => this%centurybgc_index%nelms, &
    lid_nh4=> this%centurybgc_index%lid_nh4, &
    lid_no3=> this%centurybgc_index%lid_no3, &
    lid_plant_minn_nh4 => this%centurybgc_index%lid_plant_minn_nh4, &
    lid_plant_minn_no3 => this%centurybgc_index%lid_plant_minn_no3, &
    lid_co2_hr => this%centurybgc_index%lid_co2_hr, &
    lid_minp_soluble =>  this%centurybgc_index%lid_minp_soluble,  &
    lid_minp_secondary => this%centurybgc_index%lid_minp_secondary, &
    lid_minp_occlude => this%centurybgc_index%lid_minp_occlude, &
    lid_plant_minp => this%centurybgc_index%lid_plant_minp  &
  )

  c_mass = 0._r8; n_mass = 0._r8; p_mass = 0._r8;
  if(present(c_flx))c_flx=0._r8
  if(present(n_flx))n_flx=0._r8
  if(present(p_flx))p_flx=0._r8

  jj=lit1;kc = (jj-1)*nelms + c_loc;kn = (jj-1)*nelms + n_loc;kp = (jj-1)*nelms + p_loc
  c_mass = c_mass + ystates1(kc);n_mass=n_mass + ystates1(kn); p_mass = p_mass + ystates1(kp)

  jj=lit2;kc = (jj-1)*nelms + c_loc;kn = (jj-1)*nelms + n_loc;kp = (jj-1)*nelms + p_loc
  c_mass = c_mass + ystates1(kc);n_mass=n_mass + ystates1(kn); p_mass = p_mass + ystates1(kp)

  jj=lit3;kc = (jj-1)*nelms + c_loc;kn = (jj-1)*nelms + n_loc;kp = (jj-1)*nelms + p_loc
  c_mass = c_mass + ystates1(kc);n_mass=n_mass + ystates1(kn); p_mass = p_mass + ystates1(kp)

  jj=cwd;kc = (jj-1)*nelms + c_loc;kn = (jj-1)*nelms + n_loc;kp = (jj-1)*nelms + p_loc
  c_mass = c_mass + ystates1(kc);n_mass=n_mass + ystates1(kn); p_mass = p_mass + ystates1(kp)

  jj=lwd;kc = (jj-1)*nelms + c_loc;kn = (jj-1)*nelms + n_loc;kp = (jj-1)*nelms + p_loc
  c_mass = c_mass + ystates1(kc);n_mass=n_mass + ystates1(kn); p_mass = p_mass + ystates1(kp)

  jj=fwd;kc = (jj-1)*nelms + c_loc;kn = (jj-1)*nelms + n_loc;kp = (jj-1)*nelms + p_loc
  c_mass = c_mass + ystates1(kc);n_mass=n_mass + ystates1(kn); p_mass = p_mass + ystates1(kp)

  jj=som1;kc = (jj-1)*nelms + c_loc;kn = (jj-1)*nelms + n_loc;kp = (jj-1)*nelms + p_loc
  c_mass = c_mass + ystates1(kc);n_mass=n_mass + ystates1(kn); p_mass = p_mass + ystates1(kp)

  jj=som2;kc = (jj-1)*nelms + c_loc;kn = (jj-1)*nelms + n_loc;kp = (jj-1)*nelms + p_loc
  c_mass = c_mass + ystates1(kc);n_mass=n_mass + ystates1(kn); p_mass = p_mass + ystates1(kp)

  jj=som3;kc = (jj-1)*nelms + c_loc;kn = (jj-1)*nelms + n_loc;kp = (jj-1)*nelms + p_loc
  c_mass = c_mass + ystates1(kc);n_mass=n_mass + ystates1(kn); p_mass = p_mass + ystates1(kp)


  n_mass = n_mass + ystates1(lid_nh4) + ystates1(lid_no3)

  p_mass = p_mass + ystates1(lid_minp_soluble) + ystates1(lid_minp_secondary) + ystates1(lid_minp_occlude)

  c_mass = c_mass * catomw
  n_mass = n_mass * natomw
  p_mass = p_mass * patomw
  if(present(c_flx))then
    c_flx = c_flx + ystates1(lid_co2_hr)
    c_flx = c_flx * catomw
  endif
  if(present(n_flx))then
    n_flx=n_flx + ystates1(lid_plant_minn_nh4) + ystates1(lid_plant_minn_no3) &
     + ystates1(lid_n2o_nit) + ystates1(lid_no3_den)
    n_flx = n_flx * natomw
  endif
  if(present(p_flx))then
    p_flx=p_flx + ystates1(lid_plant_minp)
    p_flx = p_flx * patomw
  endif

  end associate
  end subroutine sumup_cnp_msflx

  !--------------------------------------------------------------------
  subroutine sumup_cnp_mass(this, header, c_mass, n_mass, p_mass)
  implicit none
  class(centurybgceca_type)     , intent(in) :: this
  character(len=*), intent(in) :: header
  real(r8), intent(out) :: c_mass, n_mass, p_mass
  !local variables

  integer :: kc, kn, kp, jj
  associate(                        &
    c_loc=>  this%centurybgc_index%c_loc,&
    n_loc=>  this%centurybgc_index%n_loc,&
    p_loc=>  this%centurybgc_index%p_loc,&
    lid_n2o => this%centurybgc_index%lid_n2o,&
    lid_n2 => this%centurybgc_index%lid_n2,&
    lit1 =>  this%centurybgc_index%lit1, &
    lit2 =>  this%centurybgc_index%lit2, &
    lit3 =>  this%centurybgc_index%lit3, &
    cwd =>   this%centurybgc_index%cwd, &
    lwd =>   this%centurybgc_index%lwd, &
    fwd =>   this%centurybgc_index%fwd, &
    som1 =>  this%centurybgc_index%som1, &
    som2 =>  this%centurybgc_index%som2, &
    som3 =>  this%centurybgc_index%som3, &
    nelms => this%centurybgc_index%nelms, &
    lid_nh4=> this%centurybgc_index%lid_nh4, &
    lid_no3=> this%centurybgc_index%lid_no3, &
    lid_plant_minn_nh4 => this%centurybgc_index%lid_plant_minn_nh4, &
    lid_plant_minn_no3 => this%centurybgc_index%lid_plant_minn_no3, &
    lid_co2 => this%centurybgc_index%lid_co2, &
    lid_minp_soluble =>  this%centurybgc_index%lid_minp_soluble,  &
    lid_minp_secondary => this%centurybgc_index%lid_minp_secondary, &
    lid_minp_occlude => this%centurybgc_index%lid_minp_occlude, &
    lid_plant_minp => this%centurybgc_index%lid_plant_minp, &
    ystates1 => this%ystates1  &
  )
  print*,trim(header)
  c_mass = 0._r8; n_mass = 0._r8; p_mass = 0._r8;

  jj=lit1;kc = (jj-1)*nelms + c_loc;kn = (jj-1)*nelms + n_loc;kp = (jj-1)*nelms + p_loc
  c_mass = c_mass + ystates1(kc);n_mass=n_mass + ystates1(kn); p_mass = p_mass + ystates1(kp)
!  write(*,'(A,3(X,E20.10))')'lit1 c n p',ystates1(kc),ystates1(kn),ystates1(kp)
  jj=lit2;kc = (jj-1)*nelms + c_loc;kn = (jj-1)*nelms + n_loc;kp = (jj-1)*nelms + p_loc
  c_mass = c_mass + ystates1(kc);n_mass=n_mass + ystates1(kn); p_mass = p_mass + ystates1(kp)

  jj=lit3;kc = (jj-1)*nelms + c_loc;kn = (jj-1)*nelms + n_loc;kp = (jj-1)*nelms + p_loc
  c_mass = c_mass + ystates1(kc);n_mass=n_mass + ystates1(kn); p_mass = p_mass + ystates1(kp)

  jj=cwd;kc = (jj-1)*nelms + c_loc;kn = (jj-1)*nelms + n_loc;kp = (jj-1)*nelms + p_loc
  c_mass = c_mass + ystates1(kc);n_mass=n_mass + ystates1(kn); p_mass = p_mass + ystates1(kp)

  jj=lwd;kc = (jj-1)*nelms + c_loc;kn = (jj-1)*nelms + n_loc;kp = (jj-1)*nelms + p_loc
  c_mass = c_mass + ystates1(kc);n_mass=n_mass + ystates1(kn); p_mass = p_mass + ystates1(kp)

  jj=fwd;kc = (jj-1)*nelms + c_loc;kn = (jj-1)*nelms + n_loc;kp = (jj-1)*nelms + p_loc
  c_mass = c_mass + ystates1(kc);n_mass=n_mass + ystates1(kn); p_mass = p_mass + ystates1(kp)

  jj=som1;kc = (jj-1)*nelms + c_loc;kn = (jj-1)*nelms + n_loc;kp = (jj-1)*nelms + p_loc
  c_mass = c_mass + ystates1(kc);n_mass=n_mass + ystates1(kn); p_mass = p_mass + ystates1(kp)

  jj=som2;kc = (jj-1)*nelms + c_loc;kn = (jj-1)*nelms + n_loc;kp = (jj-1)*nelms + p_loc
  c_mass = c_mass + ystates1(kc);n_mass=n_mass + ystates1(kn); p_mass = p_mass + ystates1(kp)

  jj=som3;kc = (jj-1)*nelms + c_loc;kn = (jj-1)*nelms + n_loc;kp = (jj-1)*nelms + p_loc
  c_mass = c_mass + ystates1(kc);n_mass=n_mass + ystates1(kn); p_mass = p_mass + ystates1(kp)


!x  print*,'total om',c_mass,n_mass,p_mass
  c_mass = c_mass + ystates1(lid_co2)
!x  print*,'hr',ystates1(lid_co2_hr)
!x  print*,'no3, n2o, no3_den',ystates1(lid_nh4),ystates1(lid_no3), ystates1(lid_n2o_nit), ystates1(lid_no3_den)
!  write(*,'(A,3(X,E20.10))')'minp',ystates1(lid_minp_soluble), ystates1(lid_minp_secondary), ystates1(lid_minp_occlude)
  n_mass = n_mass + ystates1(lid_nh4) + ystates1(lid_no3)+ 2._r8*ystates1(lid_n2) +2._r8*ystates1(lid_n2o) + &
           ystates1(lid_plant_minn_nh4) + ystates1(lid_plant_minn_no3)

  p_mass = p_mass + ystates1(lid_minp_soluble) + ystates1(lid_minp_secondary) + ystates1(lid_minp_occlude) + &
           ystates1(lid_plant_minp)

!  write(*,'(A,6(X,E20.10))')'c,n,p mass',c_mass,n_mass,p_mass,ystates1(lid_plant_minn_nh4), &
!     ystates1(lid_plant_minn_no3), ystates1(lid_plant_minp)
  if(p_mass>1.e10_r8)stop
  end associate
  end subroutine sumup_cnp_mass

  !-------------------------------------------------------------------------------
  subroutine ode_adapt_ebbks1(me, y0, nprimeq, neq, t, dt, y)
    ! !DESCRIPTION:
    !first order implicit bkks ode integration with the adaptive time stepping
    !This could be used as an example for the implementation of time-adaptive
    !mbbks1.
    ! !NOTE:
    ! this code should only be used for mass positive ODE integration
    use ODEMOD, only : ebbks, get_rerr, get_tscal
    implicit none
    ! !ARGUMENTS:
    class(centurybgceca_type),  intent(inout)  :: me
    integer,  intent(in)  :: neq      ! number of equations
    real(r8), intent(in)  :: y0(neq)  ! state variable at previous time step
    real(r8), intent(in)  :: t        ! time stamp
    real(r8), intent(in)  :: dt       ! time stepping
    integer,  intent(in)  :: nprimeq  !
    real(r8), intent(out) :: y(neq)   ! updated state variable

    ! !LOCAL VARIABLES:
    real(r8) :: yc(neq)    !coarse time stepping solution
    real(r8) :: yf(neq)    !fine time stepping solution
    real(r8) :: ycp(neq)   !temporary variable
    real(r8) :: f(neq)   ! derivative
    real(r8) :: dt2
    real(r8) :: dtr
    real(r8) :: dt05
    real(r8) :: dtmin
    real(r8) :: tt,tt2     !temporary variables
    logical  :: acc
    real(r8) :: rerr, dt_scal, pscal
    integer  :: n, nJ

    dt2=dt
    dtmin=dt/64._r8
    dtr=dt
    tt=0._r8
    !make a copy of the solution at the current time step
    y=y0
    do
       if(dt2<=dtmin)then
          call me%bgc_integrate(y, dt2, tt, nprimeq, neq, f)
          call ebbks(y, f, nprimeq, neq, dt2, yc, pscal)
          dtr=dtr-dt2
          tt=tt+dt2
          y=yc
       else
          !get coarse grid solution
          call me%bgc_integrate(y, dt2, tt, nprimeq, neq, f)
          call ebbks(y, f, nprimeq, neq, dt2, yc, pscal)

          !get fine grid solution
          dt05=dt2*0.5_r8
          call ebbks(y,f,nprimeq, neq,dt05, yf, pscal)
          tt2=tt+dt05
          ycp=yf
          call me%bgc_integrate(ycp, dt05, tt, nprimeq, neq, f)
          call ebbks(ycp,f,nprimeq, neq,dt05,yf,pscal)

          !determine the relative error
          rerr=get_rerr(yc,yf, neq)*exp(1._r8-1._r8/pscal)

          !determine time scalar factor
          call get_tscal(rerr,dt_scal,acc)

          if(acc)then
             dtr=dtr-dt2
             tt=tt+dt2
             y=yf
          endif
          dt2=dt2*dt_scal
          dt2=min(dt2,dtr)
       endif
       if(abs(dtr/dt)<1.e-4_r8)exit
    enddo
  end subroutine ode_adapt_ebbks1
end module BgcCentCnpType
