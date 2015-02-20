module BGCCenturySubMod
#include "shr_assert.h"
!
! DESCRIPTION
! module contains subroutine for the century bgc implementation
  use shr_kind_mod       , only : r8 => shr_kind_r8
  use shr_log_mod        , only : errMsg => shr_log_errMsg
  use decompMod          , only : bounds_type
  use clm_varcon         , only : spval  
  implicit none
  save
  
  
  type, public :: centurybgc_type
  
  integer :: nom_pools   !not include coarse wood debris
  integer :: lit1  
  integer :: lit2  
  integer :: lit3  
  integer :: som1  
  integer :: som2  
  integer :: som3    
  integer :: cwd   
  integer :: c_loc
  integer :: n_loc
  integer :: nelms                 !number of chemical elements in an om pool
  real(r8) :: k_decay_lit1
  real(r8) :: k_decay_lit2 
  real(r8) :: k_decay_lit3 
  real(r8) :: k_decay_som1 
  real(r8) :: k_decay_som2 
  real(r8) :: k_decay_som3 
  real(r8) :: k_decay_cwd  
  
  
  integer :: lid_o2          !local position of o2 in the state variable vector
  integer :: lid_co2         !local position of co2 in the state variable vector
  integer :: lid_nh4         !local position of nh4 in the state variable vector
  integer :: lid_no3         !local position of no3 in the state variable vector
  integer :: lid_n2
  integer :: lid_n2o
  integer :: lid_plant_minn  !local position of plant consumption of mineral nitrogen in the state variable vector
  integer :: lid_n2o_nit     !n2o production from nitrification, used to for mass balance book keeping
  integer :: lid_co2_hr      !co2 production from heterotrophic respiration
  
  
  integer :: nstvars            !number of equations for the state variabile vector
  integer :: nreactions      !seven decomposition pathways plus nitrification, denitrification and plant immobilization

  real(r8), pointer :: t_scalar_col(:,:)
  real(r8), pointer :: w_scalar_col(:,:)
  real(r8), pointer :: o_scalar_col(:,:)
  real(r8), pointer :: depth_scalar_col(:,:)
   
  contains
    procedure, public  :: Init
    procedure, private :: Init_pars
    procedure, private :: InitAllocate
  end type centurybgc_type
  
    
  
 contains
 
  subroutine Init(this, bounds, lbj, ubj)
  !
  ! DESCRIPTION
  ! Initialize centurybgc type
  
  use ncdio_pio    , only: file_desc_t  

  class(centurybgc_type) :: this
  type(bounds_type),                intent(in) :: bounds
  integer,                          intent(in) :: lbj, ubj                           ! lower and upper bounds, make sure they are > 0

  type(file_desc_t) :: ncid
  
  ncid%fh=10

  call this%init_pars()
  
  call this%InitAllocate(bounds, lbj, ubj)
  

  
  end subroutine Init
!-------------------------------------------------------------------------------

  subroutine Init_pars(this)

  !
  !DESCRIPTION
  !describe the layout of the stoichiometric matrix for the reactions
  !           r{1} r{2} r{3} r{4} ... r{n}
  ! s{1}
  ! s{2}
  ! s{3}
  ! s{4}
  ! ...
  ! s{n}
  ! s{n+1}
  ! ...
  ! s{m}
  ! each reaction is associated with a primary species, the secondary species follows after primary species
  ! for the century model, the primary species are seven om pools and nh4, no3 and plant nitrogen 
  ! 

  class(centurybgc_type) :: this
  
  this%nom_pools = 7   !not include coarse wood debris
  this%lit1 = 1
  this%lit2 = 2
  this%lit3 = 3
  this%som1 = 4
  this%som2 = 5
  this%som3 = 6  
  this%cwd  = 7
  
  this%nelms = 2   !carbon and nitrogen
  this%c_loc = 1
  this%n_loc = 2
  this%lid_nh4 = this%nom_pools*this%nelms        + 1   !this is also used to indicate the nitrification reaction
  this%lid_no3 = this%nom_pools*this%nelms        + 2   !this is also used to indicate the denitrification reaction
  this%lid_plant_minn = this%nom_pools*this%nelms + 3  
  this%lid_o2  = this%nom_pools*this%nelms        + 4
  this%lid_co2 = this%nom_pools*this%nelms        + 5
  this%lid_n2o = this%nom_pools*this%nelms        + 6
  this%lid_n2  = this%nom_pools*this%nelms        + 7
  this%lid_n2o_nit = this%nom_pools*this%nelms    + 8
  this%lid_co2_hr  = this%nom_pools*this%nelms    + 9
  
  this%nstvars = this%nom_pools*this%nelms + 9        !totally 16 state variables
  
  this%nreactions = this%nom_pools + 3     !seven decomposition pathways plus nitrification, denitrification and plant immobilization  
  end subroutine Init_pars
!-------------------------------------------------------------------------------

  subroutine InitAllocate(this, bounds, lbj, ubj)
  

  class(centurybgc_type)                       :: this
  type(bounds_type)               , intent(in) :: bounds
  integer                         , intent(in) :: lbj, ubj                           ! lower and upper bounds, make sure they are > 0
  
  
  allocate(this%t_scalar_col(bounds%begc:bounds%endc,     lbj:ubj))
  allocate(this%w_scalar_col(bounds%begc:bounds%endc,     lbj:ubj))
  allocate(this%o_scalar_col(bounds%begc:bounds%endc,     lbj:ubj))
  allocate(this%depth_scalar_col(bounds%begc:bounds%endc, lbj:ubj))
  
  end subroutine InitAllocate
  


!-------------------------------------------------------------------------------
  subroutine init_state_vector(bounds, lbj, ubj, numf, filter, jtops, neq, tracerstate_vars, betrtracer_vars, centurybgc_vars, y0)
  !
  ! number of equations, total number of carbon pools + o2 + co2
  use tracerstatetype       , only : tracerstate_type
  use BeTRTracerType        , only : betrtracer_type
  

  type(bounds_type)            , intent(in) :: bounds
  integer                      , intent(in) :: lbj, ubj
  integer                      , intent(in) :: jtops(bounds%begc:bounds%endc)        ! top label of each column
  integer                      , intent(in) :: numf
  integer                      , intent(in) :: filter(:)
  integer                      , intent(in) :: neq
  type(betrtracer_type)        , intent(in) :: betrtracer_vars                    ! betr configuration information
  type(tracerstate_type)       , intent(in) :: tracerstate_vars
  type(centurybgc_type)        , intent(in) :: centurybgc_vars
  real(r8)                     , intent(out):: y0(neq, bounds%begc:bounds%endc, lbj:ubj)

  integer :: fc, c, j
  ! all organic matter pools are distributed into solid passive tracers
  do fc = 1, numf
    c = filter(fc)
    do j = jtops(c), ubj
      y0(1:centurybgc_vars%nom_pools*centurybgc_vars%nelms, j, c)    = tracerstate_vars%tracer_conc_solid_passive_col(c, j, :)

      y0(centurybgc_vars%lid_o2,  j, c)        = tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_o2) 

      y0(centurybgc_vars%lid_co2, j, c)        = tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_co2x)
  
      y0(centurybgc_vars%lid_nh4,j, c)         = tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_nh3x)
      
      y0(centurybgc_vars%lid_no3,j, c)         = tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_no3x)
      
      y0(centurybgc_vars%lid_plant_minn, j, c) = 0._r8  !initialize plant nitrogen uptake during the time step to zero
    enddo
  enddo
  end subroutine init_state_vector



!-------------------------------------------------------------------------------
  subroutine calc_som_deacyK(bounds, lbj, ubj, numf, filter, jtops, nom_pools, tracercoeff_vars, tracerstate_vars, &
    betrtracer_vars, centurybgc_vars,  k_decay)
  !
  ! DESCRIPTION
  ! calculate decay coefficients for different pools
  use tracercoeffType       , only : tracercoeff_type 
  use BetrTracerType        , only : betrtracer_type
  use tracerstatetype       , only : tracerstate_type
  use BeTRTracerType        , only : betrtracer_type
  use BGCCenturyParMod      , only : CNDecompBgcParamsInst 

  integer,                     intent(in) :: nom_pools
  type(bounds_type),           intent(in) :: bounds
  integer,                     intent(in) :: lbj, ubj
  integer,                     intent(in) :: jtops(bounds%begc:bounds%endc)        ! top label of each column
  integer,                     intent(in) :: numf
  integer,                     intent(in) :: filter(:)
  type(betrtracer_type),       intent(in) :: betrtracer_vars                    ! betr configuration information
  type(centurybgc_type),       intent(in) :: centurybgc_vars    
  type(tracercoeff_type),      intent(in) :: tracercoeff_vars
  type(tracerstate_type),   intent(inout) :: tracerstate_vars
  real(r8),                   intent(out) :: k_decay(nom_pools, bounds%begc:bounds%endc, lbj:ubj)
  
  !local variables
  integer :: fc, c, j

  associate(                                              & 
    t_scalar       => centurybgc_vars%t_scalar_col      , & ! Output: [real(r8) (:,:)   ]  soil temperature scalar for decomp                     
    w_scalar       => centurybgc_vars%w_scalar_col      , & ! Output: [real(r8) (:,:)   ]  soil water scalar for decomp                           
    o_scalar       => centurybgc_vars%o_scalar_col      , & ! Output: [real(r8) (:,:)   ]  fraction by which decomposition is limited by anoxia   
    depth_scalar   => centurybgc_vars%depth_scalar_col  , & ! Output: [real(r8) (:,:,:) ]  rate constant for decomposition (1./sec)
    lit1           => centurybgc_vars%lit1              , & !
    lit2           => centurybgc_vars%lit2              , & !
    lit3           => centurybgc_vars%lit3              , & !
    som1           => centurybgc_vars%som1              , & !
    som2           => centurybgc_vars%som2              , & !
    som3           => centurybgc_vars%som3              , & !
    cwd            => centurybgc_vars%cwd               , & !
    k_decay_lit1   => CNDecompBgcParamsInst%k_decay_lit1      , & !  
    k_decay_lit2   => CNDecompBgcParamsInst%k_decay_lit2      , & !  
    k_decay_lit3   => CNDecompBgcParamsInst%k_decay_lit3      , & !  
    k_decay_som1   => CNDecompBgcParamsInst%k_decay_som1      , & !  
    k_decay_som2   => CNDecompBgcParamsInst%k_decay_som2      , & !  
    k_decay_som3   => CNDecompBgcParamsInst%k_decay_som3      , & !  
    k_decay_cwd    => CNDecompBgcParamsInst%k_decay_cwd         & !  
  ) 
  
  k_decay(:, :, :) = spval
  do fc = 1, numf
    c = filter(fc)
    do j = jtops(c), ubj
      k_decay(lit1, j, c) = k_decay_lit1 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) * depth_scalar(c,j)
      k_decay(lit2, j, c) = k_decay_lit2 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) * depth_scalar(c,j)
      k_decay(lit3, j, c) = k_decay_lit3 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) * depth_scalar(c,j)
      k_decay(som1, j, c) = k_decay_som1 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) * depth_scalar(c,j)
      k_decay(som2, j, c) = k_decay_som2 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) * depth_scalar(c,j)
      k_decay(som3, j, c) = k_decay_som3 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) * depth_scalar(c,j)
      k_decay(cwd,  j, c) = k_decay_cwd  * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) * depth_scalar(c,j)
    enddo  
  enddo
  end associate
  end subroutine calc_som_deacyK  

!-------------------------------------------------------------------------------
  subroutine calc_cascade_matrix(nstvars, nreactions, cn_ratios, cp_ratios, n2_n2o_ratio_denit, nh4_no3_ratio, pct_sand, centurybgc_vars, cascade_matrix)
  !
  ! DESCRIPTION
  ! calculate cascade matrix for the decomposition model
  !
  use clm_varcon                , only: nitrif_n2o_loss_frac

  integer                       , intent(in) :: nstvars
  integer                       , intent(in) :: nreactions
  type(centurybgc_type)         , intent(in) :: centurybgc_vars
  real(r8)                      , intent(in) :: cn_ratios(centurybgc_vars%nom_pools)
  real(r8)                      , intent(in) :: cp_ratios(centurybgc_vars%nom_pools)
  real(r8)                      , intent(in) :: n2_n2o_ratio_denit                   !ratio of n2 to n2o during denitrification
  real(r8)                      , intent(in) :: nh4_no3_ratio                        !ratio of nh4 to no3 at current time step
  real(r8)                      , intent(in) :: pct_sand
  
  real(r8),intent(out) :: cascade_matrix(nstvars, nreactions)
  
  real(r8) :: ftxt, f1, f2
  integer :: k
  
  
  associate(                                            & !
    lit1      => centurybgc_vars%lit1                 , & !
    lit2      => centurybgc_vars%lit2                 , & !
    lit3      => centurybgc_vars%lit3                 , & !
    som1      => centurybgc_vars%som1                 , & !
    som2      => centurybgc_vars%som2                 , & !
    som3      => centurybgc_vars%som3                 , & !
    cwd       => centurybgc_vars%cwd                  , & !
    c_loc     => centurybgc_vars%c_loc                , & !
    n_loc     => centurybgc_vars%n_loc                , & !
    nelms     => centurybgc_vars%nelms                , & !
    lid_o2    => centurybgc_vars%lid_o2               , & !
    lid_co2   => centurybgc_vars%lid_co2              , & !
    lid_nh4   => centurybgc_vars%lid_nh4              , & !
    lid_no3   => centurybgc_vars%lid_no3              , & !
    lid_n2o   => centurybgc_vars%lid_n2o              , & !
    lid_n2    => centurybgc_vars%lid_n2               , & !  
    lid_plant_minn => centurybgc_vars%lid_plant_minn    & !
  )
  !initialize all entries to zero
  cascade_matrix = 0._r8
  
  !note all reactions are in the form products - substrates = 0, therefore
  !mass balance is automatically ensured.
  
  
  !reaction1, lit1 -> s1
  !lit1 + 0.55*o2 -> 0.45 som1 + 0.55co2 + (1/cn_ratios(lit1) - 0.45/cn_ratios(som1))min_n+ (1/cp_ratios(lit1)-0.45/cp_ratios(som1))min_p
  cascade_matrix((lit1-1)*nelms+c_loc   ,1)  = -1._r8
  cascade_matrix((lit1-1)*nelms+n_loc   ,1)  = -1._r8
  
  cascade_matrix(lid_o2 ,1)  = -CNDecompBgcParamsInst%rf_l1s1_bgc
  cascade_matrix((som1-1)*nelms+c_loc   ,1)  = 1._r8-CNDecompBgcParamsInst%rf_l1s1_bgc
  cascade_matrix((som1-1)*nelms+n_loc   ,1)  = 1._r8-CNDecompBgcParamsInst%rf_l1s1_bgc
  cascade_matrix(lid_co2,1)  =  CNDecompBgcParamsInst%rf_l1s1_bgc
  cascade_matrix(lid_nh4,1)  =  1._r8/cn_ratios(lit1) - 1._r8/cn_ratios(som1)


  !reaction 2, lit2 -> s1
  !lit2 + 0.5 o2  -> 0.5 som1 + 0.5 co2 + (1/cn_ratios(lit2)-0.5/cn_ratios(som1))min_n +(1/cp_ratios(lit2)-0.5/cp_ratios(som1))min_p
  cascade_matrix((lit2-1)*nelms+c_loc   ,2)   = -1._r8
  cascade_matrix((lit2-1)*nelms+n_loc   ,2)   = -1._r8
  
  cascade_matrix(lid_o2 ,2)   = -CNDecompBgcParamsInst%rf_l2s1_bgc
  cascade_matrix((som1-1)*nelms+c_loc   ,2)   =  1._r8-CNDecompBgcParamsInst%rf_l2s1_bgc
  cascade_matrix((som1-1)*nelms+n_loc   ,2)   =  1._r8-CNDecompBgcParamsInst%rf_l2s1_bgc
  
  cascade_matrix(lid_co2,2)   =  CNDecompBgcParamsInst%rf_l2s1_bgc
  cascade_matrix(lid_nh4,2)   = 1._r8/cn_ratios(lit2) - 0.5_r8/cn_ratios(som1)
  
  !reaction 3, lit3->s2
  !lit3 + 0.5 o2 -> 0.5 som2 + 0.5 co2 + (1/cn_ratios(lit3) - 0.5/cn_ratios(som2))min_n + (1/cp_ratios(lit3)-0.5_r8/cp_ratios(som2))minp
  cascade_matrix((lit3-1)*nelms+c_loc   ,3) = -1._r8
  cascade_matrix((lit3-1)*nelms+n_loc   ,3) = -1._r8
  
  cascade_matrix(lid_o2 ,3) = -CNDecompBgcParamsInst%rf_l3s2_bgc
  cascade_matrix((som2-1)*nelms+c_loc   ,3) =  1._r8-CNDecompBgcParamsInst%rf_l3s2_bgc
  cascade_matrix((som2-1)*nelms+n_loc   ,3) =  1._r8-CNDecompBgcParamsInst%rf_l3s2_bgc
  
  cascade_matrix(lid_co2,3) = CNDecompBgcParamsInst%rf_l3s2_bgc
  cascade_matrix(lid_nh4,3) = 1._r8/cn_ratios(lit3) - 0.5_r8/cn_ratios(som2)
  
  !double check those stoichiometry parameters
  !reaction 4, the partition into som2 and som3 is soil texture dependent
  !som1 + f(txt) o2 -> f1*som2 + f2*som3 + f(txt) co2 + (1/cn_ratios(som1)-f1/cn_ratios(som2)-f2/cn_ratios(som3)) +(1/cp_ratios(som1)-f1/cp_ratios(som2)-f2/cp_ratios(som3))
  !f(txt) = 0.85_r8 - 0.68_r8 * 0.01_r8 * (100._r8 - sand)
  !f1+f2+f(txt)=1._r8
  ftxt = 0.85_r8 - 0.68_r8 * 0.01_r8 * (100._r8 - pct_sand)
  f1 = 0.996*(1._r8-ftxt)
  f2 = 0.004*(1._r8-ftxt)
  cascade_matrix((som1-1)*nelms+c_loc   ,4)  = -1._r8
  cascade_matrix((som1-1)*nelms+n_loc   ,4)  = -1._r8
  
  cascade_matrix(lid_o2 ,4) = -ftxt
  cascade_matrix((som3-1)*nelms+c_loc   ,4)  = 0.004*(1._r8-ftxt)
  cascade_matrix((som3-1)*nelms+n_loc   ,4)  = 0.004*(1._r8-ftxt)
  
  cascade_matrix((som2-1)*nelms+c_loc   ,4) = 0.996*(1._r8-ftxt)
  cascade_matrix((som2-1)*nelms+n_loc   ,4) = 0.996*(1._r8-ftxt)
  
  cascade_matrix(lid_co2,4) = ftxt
  cascade_matrix(lid_nh4,4) = 1._r8/cn_ratios(som1)-f1/cn_ratios(som2)-f2/cn_ratios(som3)
  
  !reaction 5, som2->som1, som3
  !som2 + 0.55 o2 -> 0.42 som1 + 0.03som3 + 0.55co2 + (1/cn_ratios(som2)-0.42/cn_ratios(som1)-0.03/cn_ratios(som3)) + (1/cp_raitos(som2)-0.42/cp_ratios(som1)-0.03/cp_ratios(som3))
  cascade_matrix((som2-1)*nelms+c_loc   ,5)   = -1._r8
  cascade_matrix((som2-1)*nelms+n_loc   ,5)   = -1._r8
  
  cascade_matrix(lid_o2 ,5)   = -0.55_r8
  cascade_matrix((som1-1)*nelms+c_loc   ,5)   =  0.42_r8
  cascade_matrix((som1-1)*nelms+n_loc   ,5)   =  0.42_r8
  
  cascade_matrix((som3-1)*nelms+c_loc   ,5)   =  0.03_r8
  cascade_matrix((som3-1)*nelms+n_loc   ,5)   =  0.03_r8
  
  cascade_matrix(lid_co2,5)   =  0.55_r8
  cascade_matrix(lid_nh4,5)   = 1._r8/cn_ratios(som2)-0.42_r8/cn_ratios(som1)-0.03_r8/cn_ratios(som3)
  
  !reaction 6, s3-> s1
  !som3 + 0.55 o2 -> 0.45*som1 + 0.55co2 + (1/cn_ratios(som3)-0.45/cn_ratios(som1)) + (1/cp_ratios(som3)-0.45/cp_ratios(som1))
  cascade_matrix((som3-1)*nelms+c_loc   ,6) = -1._r8
  cascade_matrix((som3-1)*nelms+n_loc   ,6) = -1._r8
  
  cascade_matrix(lid_o2 ,6) = -CNDecompBgcParamsInst%rf_s3s1_bgc
  cascade_matrix((som1-1)*nelms+c_loc   ,6) = 1._r8-CNDecompBgcParamsInst%rf_s3s1_bgc
  cascade_matrix((som1-1)*nelms+n_loc   ,6) = 1._r8-CNDecompBgcParamsInst%rf_s3s1_bgc
  
  cascade_matrix(lid_co2,6) = CNDecompBgcParamsInst%rf_s3s1_bgc
  cascade_matrix(lid_nh4,6) = 1._r8/cn_ratios(som3) - 0.45_r8/cn_ratios(som1)
  
  !reaction 7, the partition into lit1 and lit2 is nutrient dependent, respires co2?
  !cwd + o2 -> 0.76lit2 + 0.24*lit3 + (1/cn_ratios(cwd)-0.76/cn_ratios(lit2)-0.24/cn_ratios(lit3)) + (1/cp_ratios(cwd)-0.76/cp_ratios(lit2)-0.24/cp_ratios(lit3))
  cascade_matrix((cwd-1)*nelms+c_loc    ,7) = -1._r8
  cascade_matrix((cwd-1)*nelms+n_loc    ,7) = -1._r8
  
  cascade_matrix((lit2-1)*nelms+c_loc   ,7) = CNDecompBgcParamsInst%cwd_fcel_bgc
  cascade_matrix((lit2-1)*nelms+n_loc   ,7) = CNDecompBgcParamsInst%cwd_fcel_bgc
  
  cascade_matrix((lit3-1)*nelms+c_loc   ,7) = CNDecompBgcParamsInst%cwd_flig_bgc
  cascade_matrix((lit3-1)*nelms+n_loc   ,7) = CNDecompBgcParamsInst%cwd_flig_bgc
  
  cascade_matrix(lid_nh4,7) = 1._r8/cn_ratios(cwd) - 0.76_r8/cn_ratios(lit2) - 0.24_r8/cn_ratios(lit3)
    
  do k = 1, 7
    !Note: Jinyun Tang, Dec 26, 2014
    !When a reaction needs mineral nitrogen to balance the elements, it takes mineral nitrogen proportionally from nh4 and no3.
    !This formulation assumes that the nitrogen mineralized from om decomposition is equally accessible to plants and decomposers. Such
    !a formulation is different from the century BGC in CLM4.5. Rather, CLM4.5 bgc assumes that the nitrogen mineralized from nitrogen releasing
    !decomposition pathways is first used to meet the nitrogen demand from nitrogen immobilizing decomposition pathways. In the later case, the stoichiometry becomes
    !rate dependent. 
    if(cascade_matrix(lid_nh4, k)<0._r8)then
      !it requires nitrogen uptake
      cascade_matrix(lid_no3, k) = cascade_matrix(lid_nh4, k) / (1._r8 + nh4_no3_ratio)
      cascade_matrix(lid_nh4, k) = cascade_matrix(lid_nh4, k) - cascade_matrix(lid_no3, k)
    endif
  enddo
  !reaction 8, nitrification
  !NH4(+) + (2-f)O2 + (2-f)OH(-)-> (1-f)NO3(-) + (f/2)N2O + (3-f/2) H2O 
  cascade_matrix(lid_nh4 ,8) = -1._r8
  cascade_matrix(lid_o2  ,8) = -(2._r8 - nitrif_n2o_loss_frac)
  cascade_matrix(lid_no3 ,8) = 1._r8 - nitrif_n2o_loss_frac
  cascade_matrix(lid_n2o, 8) = 0.5_r8 * nitrif_n2o_loss_frac
  
  !reaction 9, denitrification
  !NO3(-) -> 0.5*f N2O + 0.5* (1-f) N2, where f is a function determined from the century denitrification model
  cascade_matrix(lid_no3 ,9) = -1._r8
  cascade_matrix(lid_n2o ,9) = 0.5_r8 * 1._r8/(1._r8+n2_n2o_ratio_denit)
  cascade_matrix(lid_n2  ,9) = 0.5_r8 * n2_n2o_ratio_denit/(1._r8+n2_n2o_ratio_denit)

  !reaction 10, plant mineral nitrogen uptake
  ! f nh4 + (1-f) no3 -> plant_nitrogen
  cascade_matrix(lid_nh4, 10)        = -nh4_no3_ratio/(1._r8+nh4_no3_ratio)
  cascade_matrix(lid_no3, 10)        = -1._r8/(1._r8 + nh4_no3_ratio)
  cascade_matrix(lid_plant_minn, 10) = 1._r8
  
  end associate
  end subroutine calc_cascade_matrix
  
!-------------------------------------------------------------------------------
  subroutine calc_sompool_decay(bounds, lbj, ubj, numf, filter, jtops, nom_pools, k_decay, om_pools, decay_rates)
  !
  ! DESCRIPTION
  ! calculate degradation for all different pools

  integer,  intent(in) :: nom_pools
  type(bounds_type),           intent(in) :: bounds
  integer,                     intent(in) :: lbj, ubj
  integer,                     intent(in) :: jtops(bounds%begc:bounds%endc)        ! top label of each column
  integer,                     intent(in) :: numf
  integer,                     intent(in) :: filter(:)  
  real(r8), intent(in) :: k_decay(nom_pools, bounds%begc:bounds%endc, lbj:ubj)
  real(r8), intent(in) :: om_pools(nom_pools,bounds%begc:bounds%endc, lbj:ubj)
  real(r8), intent(out):: decay_rates(nom_pools, bounds%begc:bounds%endc, lbj:ubj)
  
  integer :: jj, fc, c, j
  
  do fc = 1, numf
    c = filter(fc)
    do j = jtops(c), ubj
      !for om pools
      do jj = 1, nom_pools
        decay_rates(jj, c, j) = om_pools(jj, c, j) * k_decay(jj, c, j)
      enddo
    enddo
  enddo
  
  !for nitrification and denitrification
  
  end subroutine calc_sompool_decay

  
!-------------------------------------------------------------------------------  
  subroutine retrieve_nutrient_flux(bounds, lbj, ubj, numf, filter, jtops, neq, dtime, yf, y0,  tracerflux_vars, centurybgc_vars, plantsoilnutrientflux_vars)
  !
  ! DESCRIPTION
  ! retrieve the fluxes
  use tracerfluxType           , only : tracerflux_type
  use PlantSoilnutrientFluxType, only : plantsoilnutrientflux_type
  

  type(bounds_type),      intent(in) :: bounds
  integer,                intent(in) :: lbj, ubj
  integer,                intent(in) :: jtops(bounds%begc:bounds%endc)        ! top label of each column
  integer,                intent(in) :: numf
  integer,                intent(in) :: filter(:)
  
  integer                         , intent(in) :: neq
  real(r8)                        , intent(in) :: dtime
  real(r8)                        , intent(in) :: yf(neq, bounds%begc:bounds%endc, lbj:ubj) !
  real(r8)                        , intent(in) :: y0(neq, bounds%begc:bounds%endc, lbj:ubj) !
  type(centurybgc_type)           , intent(in) :: centurybgc_vars 
  type(tracerflux_type)           , intent(inout) :: tracerflux_vars
  class(plantsoilnutrientflux_type), intent(inout) :: plantsoilnutrientflux_vars
  
  integer :: fc, c, j
  
  do fc = 1, numf
    c = filter(fc)    
    do j = jtops(c), ubj
      plantsoilnutrientflux_vars%plant_minn_active_yield_flx_vr_col(c,j) = (yf(centurybgc_vars%lid_plant_minn, c, j) - y0(centurybgc_vars%lid_plant_minn, c, j))/dtime
    enddo
  enddo
  
  
  end subroutine retrieve_nutrient_flux
!-------------------------------------------------------------------------------

  subroutine retrieve_state_vector(bounds, lbj, ubj, numf, filter, jtops, neq, yf, centurybgc_vars, betrtracer_vars, tracerstate_vars)
  !
  ! number of equations, total number of carbon pools + o2 + co2
  use tracerstatetype       , only : tracerstate_type
  use BeTRTracerType        , only : betrtracer_type
  

  type(bounds_type)       , intent(in) :: bounds
  integer                 , intent(in) :: lbj, ubj
  integer                 , intent(in) :: jtops(bounds%begc:bounds%endc)        ! top label of each column
  integer                 , intent(in) :: numf
  integer                 , intent(in) :: filter(:)
  
  integer                 , intent(in) :: neq
  type(betrtracer_type)   , intent(in) :: betrtracer_vars                          ! betr configuration information
  type(centurybgc_type)   , intent(in) :: centurybgc_vars 
  real(r8)                , intent(in) :: yf(neq, bounds%begc:bounds%endc, lbj:ubj) !
  type(tracerstate_type)  , intent(inout) :: tracerstate_vars
  
  
  integer :: fc, c, j
  ! all organic matter pools are distributed into solid passive tracers
  
  do fc = 1, numf
    c = filter(fc)
    do j = jtops(c), ubj   !currently, om is added only to soil layers       
  
      tracerstate_vars%tracer_conc_solid_passive_col(c, j, :) = yf(1:centurybgc_vars%nom_pools, c, j)

      tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_o2) = yf(centurybgc_vars%lid_o2, c, j)

      tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_co2x) = yf(centurybgc_vars%lid_co2, c, j)
      
      tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_nh3x) = yf(centurybgc_vars%lid_nh4, c, j)
      
      tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_no3x) = yf(centurybgc_vars%lid_no3, c, j)      
    enddo
  enddo
  
  end subroutine retrieve_state_vector
  
!-------------------------------------------------------------------------------  

  subroutine calc_nitrif_denitrif_rate(bounds, lbj, ubj, numf, filter, jtops, dz, t_soisno, pH, &
    pot_co2_hr, anaerobic_frac, smin_nh4_vr, smin_no3_vr, soilstate_vars, waterstate_vars, centurybgc_vars, &
    n2_n2o_ratio_denit, nh4_no3_ratio, decay_nh4, decay_no3)
  !
  ! calculate nitrification denitrification rate
  ! the actual nitrification rate will be f_nitr * [nh4]
  ! and the actual denitri rate will be of f_denit * [no3]
  use clm_varcon          , only : rpi, secspday, catomw, natomw
  use SoilStatetype       , only : soilstate_type
  use WaterStateType      , only : waterstate_type
  use MathfuncMod         , only : safe_div
  use shr_const_mod       , only : SHR_CONST_TKFRZ
  
  type(bounds_type),      intent(in) :: bounds
  integer,                intent(in) :: lbj, ubj
  integer,                intent(in) :: jtops(bounds%begc: )        ! top label of each column
  integer,                intent(in) :: numf
  integer,                intent(in) :: filter(:)
  real(r8),               intent(in) :: dz(bounds%begc: , lbj: )
  real(r8),               intent(in) :: pH(bounds%begc: , lbj: )                  !pH of soil
  real(r8),               intent(in) :: t_soisno(bounds%begc: , lbj: )            !soil temperature
  real(r8),               intent(in) :: pot_co2_hr(bounds%begc: , lbj: )          !potential aerobic heteotrophic respiration, mol CO2/m3/s
  real(r8),               intent(in) :: anaerobic_frac(bounds%begc: , lbj: )      !fraction of anaerobic soil
  real(r8),               intent(in) :: smin_nh4_vr(bounds%begc: , lbj: )
  real(r8),               intent(in) :: smin_no3_vr(bounds%begc: , lbj: )   !soil no3 concentration [mol N/m3]
  type(waterstate_type) , intent(in) :: waterstate_vars  
  type(soilstate_type)  , intent(in) :: soilstate_vars  
  type(centurybgc_type) , intent(in) :: centurybgc_vars
  real(r8)             , intent(out) :: n2_n2o_ratio_denit(bounds%begc: , lbj: )  !ratio of n2 to n2o in denitrification
  real(r8)             , intent(out) :: nh4_no3_ratio(bounds%begc: , lbj: )       !ratio of soil nh4 and no3
  real(r8)             , intent(out) :: decay_nh4(bounds%begc: ,lbj: )            !1/s, decay rate of nh4
  real(r8)             , intent(out) :: decay_no3(bounds%begc: ,lbj: )            !1/s, decay rate of no3
  

  logical, parameter :: no_frozen_nitrif_denitrif = .false.                     !this is a testing parameter, just to make the model run
  
  !local variables
  real(r8) :: soil_hr_vr(bounds%begc:bounds%endc,   lbj:ubj)
  real(r8) :: k_nitr_t_vr(bounds%begc:bounds%endc,  lbj:ubj)
  real(r8) :: k_nitr_ph_vr(bounds%begc:bounds%endc, lbj:ubj)
  real(r8) :: k_nitr_h2o_vr(bounds%begc:bounds%endc,lbj:ubj)
  real(r8) :: k_nitr_vr(bounds%begc:bounds%endc,    lbj:ubj)
  real(r8) :: pot_f_nit_vr(bounds%begc:bounds%endc, lbj:ubj)
  real(r8) :: soil_bulkdensity(bounds%begc:bounds%endc, lbj:ubj)
  real(r8) :: smin_no3_massdens_vr(bounds%begc:bounds%endc,lbj:ubj)
  real(r8) :: soil_co2_prod(bounds%begc:bounds%endc, lbj:ubj)
  real(r8) :: fmax_denit_carbonsubstrate_vr(bounds%begc:bounds%endc, lbj:ubj)
  real(r8) :: fmax_denit_nitrate_vr(bounds%begc:bounds%endc, lbj:ubj)
  real(r8) :: f_denit_base_vr(bounds%begc:bounds%endc, lbj:ubj)
  real(r8) :: ratio_k1(bounds%begc:bounds%endc, lbj:ubj)
  real(r8) :: diffus(bounds%begc:bounds%endc, lbj:ubj)
  real(r8) :: ratio_no3_co2(bounds%begc:bounds%endc, lbj:ubj)
  real(r8) :: wfps_vr(bounds%begc:bounds%endc, lbj:ubj)
  real(r8) :: fr_WFPS(bounds%begc:bounds%endc, lbj:ubj)
  real(r8) :: pot_f_denit_vr(bounds%begc:bounds%endc, lbj:ubj)
  real(r8) :: co2diff_con(2) = (/0.1325_r8, 0.0009_r8/)
  real(r8) :: g_per_m3__to__ug_per_gsoil
  real(r8) :: g_per_m3_sec__to__ug_per_gsoil_day
  real(r8) :: k_nitr_max
  integer  :: fc, c, j

  SHR_ASSERT_ALL((ubound(jtops) == (/bounds%endc/)), errMsg(__FILE__,__LINE__))
  SHR_ASSERT_ALL((ubound(pH) == (/bounds%endc, ubj/)), errMsg(__FILE__,__LINE__))
  SHR_ASSERT_ALL((ubound(t_soisno) == (/bounds%endc, ubj/)), errMsg(__FILE__,__LINE__))
  SHR_ASSERT_ALL((ubound(pot_co2_hr) == (/bounds%endc, ubj/)), errMsg(__FILE__,__LINE__))
  SHR_ASSERT_ALL((ubound(anaerobic_frac) == (/bounds%endc, ubj/)), errMsg(__FILE__,__LINE__))
  SHR_ASSERT_ALL((ubound(smin_nh4_vr) == (/bounds%endc, ubj/)), errMsg(__FILE__,__LINE__))  
  SHR_ASSERT_ALL((ubound(smin_no3_vr) == (/bounds%endc, ubj/)), errMsg(__FILE__,__LINE__))
  SHR_ASSERT_ALL((ubound(n2_n2o_ratio_denit) == (/bounds%endc, ubj/)), errMsg(__FILE__,__LINE__))
  SHR_ASSERT_ALL((ubound(nh4_no3_ratio) == (/bounds%endc, ubj/)), errMsg(__FILE__,__LINE__))
  SHR_ASSERT_ALL((ubound(decay_nh4) == (/bounds%endc, ubj/)), errMsg(__FILE__,__LINE__))
  SHR_ASSERT_ALL((ubound(decay_no3) == (/bounds%endc, ubj/)), errMsg(__FILE__,__LINE__))
  
  associate(                                              &
    bd                            =>    soilstate_vars%bd_col          , & !
    watsat                        =>    soilstate_vars%watsat_col      , & !
    h2osoi_vol                    =>    waterstate_vars%h2osoi_vol_col , & !
    h2osoi_liq                    =>    waterstate_vars%h2osoi_liq_col , & !
    finundated                    =>    waterstate_vars%finundated_col , & ! Input: [real(r8) (:)]
    t_scalar                      =>    centurybgc_vars%t_scalar_col   , & ! Input: [real(r8) (:,:)   ]  soil temperature scalar for decomp                     
    w_scalar                      =>    centurybgc_vars%w_scalar_col     & ! Input: [real(r8) (:,:)   ]  soil water scalar for decomp                           
  )
  ! Set maximum nitrification rate constant 
  k_nitr_max =  0.1_r8 / secspday   ! [1/sec] 10%/day  Parton et al., 2001 


  decay_nh4(:, :) = spval
  decay_no3(:, :) = spval
  do j = lbj, ubj
    do fc = 1,numf
      c = filter(fc)
      if(j<jtops(c))cycle
      nh4_no3_ratio(c,j) = safe_div(smin_nh4_vr(c,j), smin_no3_vr(c,j))
      !---------------- nitrification
      ! follows CENTURY nitrification scheme (Parton et al., (2001, 1996))

      ! assume nitrification temp function equal to the HR scalar
      k_nitr_t_vr(c,j) = min(t_scalar(c,j), 1._r8)

      ! ph function from Parton et al., (2001, 1996)
      k_nitr_ph_vr(c,j) = 0.56 + atan(rpi * 0.45 * (-5.+ pH(c,j)))/rpi

      ! moisture function-- assume the same moisture function as limits heterotrophic respiration
      ! Parton et al. base their nitrification- soil moisture rate constants based on heterotrophic rates-- can we do the same?
      k_nitr_h2o_vr(c,j) = w_scalar(c,j)

      ! nitrification constant is a set scalar * temp, moisture, and ph scalars
      k_nitr_vr(c,j) = k_nitr_max * k_nitr_t_vr(c,j) * k_nitr_h2o_vr(c,j) * k_nitr_ph_vr(c,j)

      ! first-order decay of ammonium pool with scalar defined above
      pot_f_nit_vr(c,j) = max(k_nitr_vr(c,j), 0._r8)

      ! limit to oxic fraction of soils
      pot_f_nit_vr(c,j)  = pot_f_nit_vr(c,j) * (1._r8 - anaerobic_frac(c,j))   ![1/s]

      ! limit to non-frozen soil layers
      if ( t_soisno(c,j) <= SHR_CONST_TKFRZ .and. no_frozen_nitrif_denitrif) then
        pot_f_nit_vr(c,j) = 0._r8
      endif  
  
      !pot_f_nit_vr is the rate parameter to be returned
      decay_nh4(j,c) = pot_f_nit_vr(c,j) 
      !---------------- denitrification
      ! first some input variables an unit conversions
      soil_hr_vr(c,j) = pot_co2_hr(c,j) * catomw

      ! CENTURY papers give denitrification in units of per gram soil; need to convert from volumetric to mass-based units here
      soil_bulkdensity(c,j) = bd(c,j) + h2osoi_liq(c,j)/dz(c,j)         

      g_per_m3__to__ug_per_gsoil = 1.e3_r8 / soil_bulkdensity(c,j)

      g_per_m3_sec__to__ug_per_gsoil_day = g_per_m3__to__ug_per_gsoil * secspday

      smin_no3_massdens_vr(c,j) = max(smin_no3_vr(c,j), 0._r8) * g_per_m3__to__ug_per_gsoil * natomw

      soil_co2_prod(c,j) = (soil_hr_vr(c,j) * (g_per_m3_sec__to__ug_per_gsoil_day))

      !! maximum potential denitrification rates based on heterotrophic respiration rates or nitrate concentrations, 
      !! from (del Grosso et al., 2000)
      fmax_denit_carbonsubstrate_vr(c,j) = (0.1_r8 * (soil_co2_prod(c,j)**1.3_r8)) &
                 / g_per_m3_sec__to__ug_per_gsoil_day
      !  
      fmax_denit_nitrate_vr(c,j) = (1.15_r8 * smin_no3_massdens_vr(c,j)**0.57_r8)  &
                 / g_per_m3_sec__to__ug_per_gsoil_day

      ! find limiting denitrification rate
      f_denit_base_vr(c,j) = max(min(fmax_denit_carbonsubstrate_vr(c,j), fmax_denit_nitrate_vr(c,j)),0._r8) 

      ! limit to non-frozen soil layers
      if ( t_soisno(c,j) <= SHR_CONST_TKFRZ .and. no_frozen_nitrif_denitrif ) then
        f_denit_base_vr(c,j) = 0._r8
      endif

      ! limit to anoxic fraction of soils
      pot_f_denit_vr(c,j) = f_denit_base_vr(c,j) * anaerobic_frac(c,j)
      
      decay_no3(j,c) = safe_div(pot_f_denit_vr(c,j), smin_no3_vr(c,j)) ![1/s]

      ! now calculate the ratio of N2O to N2 from denitrifictaion, following Del Grosso et al., 2000
      ! diffusivity constant (figure 6b)
      ratio_k1(c,j) = max(1.7_r8, 38.4_r8 - 350._r8 * diffus(c,j))

      ! ratio function (figure 7c)
      if ( soil_co2_prod(c,j) > 0 ) then
        ratio_no3_co2(c,j) = smin_no3_massdens_vr(c,j) / soil_co2_prod(c,j)
      else
        ! fucntion saturates at large no3/co2 ratios, so set as some nominally large number
        ratio_no3_co2(c,j) = 100._r8
      endif

      ! total water limitation function (Del Grosso et al., 2000, figure 7a)
      wfps_vr(c,j) = max(min(h2osoi_vol(c,j)/watsat(c, j), 1._r8), 0._r8) * 100._r8
      fr_WFPS(c,j) = max(0.1_r8, 0.015_r8 * wfps_vr(c,j) - 0.32_r8)
      
      ! final ratio expression 
      n2_n2o_ratio_denit(c,j) = max(0.16*ratio_k1(c,j), ratio_k1(c,j)*exp(-0.8 * ratio_no3_co2(c,j))) * fr_WFPS(c,j)
    enddo
  enddo
  
  
  end associate
  end subroutine calc_nitrif_denitrif_rate
  
!----------------------------------------------------------------------------------------------------  
  subroutine calc_anaerobic_frac(bounds, lbj, ubj, numf, filter, jtops, t_soisno, soilstate_vars, &
    h2osoi_vol, o2_decomp_depth_unsat, conc_o2_unsat, anaerobic_frac)
  !
  ! DESCRIPTIONS
  ! calculate soil anoxia state for doing nitrification and denitrification
  ! Rewritten based on Charlie Koven's code by Jinyun Tang
  use CNSharedParamsMod   , only : CNParamsShareInst
  use clm_varcon          , only : d_con_g, grav, d_con_w
  use SoilStatetype       , only : soilstate_type
  use BGCCenturyParMod    , only : CNNitrifDenitrifParamsInst

  type(bounds_type),      intent(in) :: bounds
  integer,                intent(in) :: lbj, ubj
  integer,                intent(in) :: jtops(bounds%begc: )            ! top label of each column
  integer,                intent(in) :: numf
  integer,                intent(in) :: filter(:)                       !indices
  type(soilstate_type),   intent(in) :: soilstate_vars
  real(r8),               intent(in) :: t_soisno(bounds%begc: , lbj: )
  real(r8),               intent(in) :: h2osoi_vol(bounds%begc:, lbj: )
  real(r8),               intent(in) :: o2_decomp_depth_unsat(bounds%begc: ,lbj: )    !potential o2 consumption, as deduced from aerobic heteorotrophic decomposition, mol o2/m3/s
  real(r8),               intent(in) :: conc_o2_unsat(bounds%begc: , lbj: )           !bulk soil o2 concentration, mol/m3
  real(r8),              intent(out) :: anaerobic_frac(bounds%begc: , lbj: )          !fraction of aerobic soil
  
  real(r8), parameter :: rho_w  = 1.e3_r8        ![kg/m3]
  real(r8) :: f_a
  real(r8) :: eps
  real(r8) :: om_frac
  real(r8) :: diffus(bounds%begc:bounds%endc, lbj:ubj)
  real(r8) :: organic_max
  real(r8) :: rij_kro_a
  real(r8) :: rij_kro_alpha
  real(r8) :: rij_kro_beta
  real(r8) :: rij_kro_gamma
  real(r8) :: rij_kro_delta
  real(r8) :: surface_tension_water
  real(r8) :: r_min_sat
  real(r8) :: r_psi_sat
  real(r8) :: r_max(bounds%begc:bounds%endc, lbj:ubj)
  real(r8) :: r_min(bounds%begc:bounds%endc, lbj:ubj)
  real(r8) :: r_psi(bounds%begc:bounds%endc, lbj:ubj)  
  real(r8) :: anaerobic_frac_sat
  real(r8) :: ratio_diffusivity_water_gas(bounds%begc:bounds%endc, lbj:ubj)
  
  integer  :: fc, c, j
  
  SHR_ASSERT_ALL((ubound(jtops) == (/bounds%endc/)), errMsg(__FILE__,__LINE__))
  SHR_ASSERT_ALL((ubound(t_soisno) == (/bounds%endc, ubj/)), errMsg(__FILE__,__LINE__))
  SHR_ASSERT_ALL((ubound(o2_decomp_depth_unsat) == (/bounds%endc, ubj/)), errMsg(__FILE__,__LINE__))
  SHR_ASSERT_ALL((ubound(conc_o2_unsat) == (/bounds%endc, ubj/)), errMsg(__FILE__,__LINE__))
  SHR_ASSERT_ALL((ubound(anaerobic_frac) == (/bounds%endc, ubj/)), errMsg(__FILE__,__LINE__))
  SHR_ASSERT_ALL((ubound(h2osoi_vol) == (/bounds%endc, ubj/)), errMsg(__FILE__,__LINE__))
  
  associate(                                                         &
    watsat                        =>    soilstate_vars%watsat_col  , &
    watfc                         =>    soilstate_vars%watfc_col   , &
    bsw                           =>    soilstate_vars%bsw_col     , &
    sucsat                        =>    soilstate_vars%sucsat_col  , &
    soilpsi                       =>    soilstate_vars%soilpsi_col , &
    cellorg                       =>    soilstate_vars%cellorg_col   &
  )
  ! Todo:  FIX(SPM,032414) - the explicit divide gives different results than when that
  ! value is placed in the parameters netcdf file.  To get bfb, keep the 
  ! divide in source.
  !k_nitr_max = CNNitrifDenitrifParamsInst%k_nitr_max

  surface_tension_water = CNNitrifDenitrifParamsInst%surface_tension_water

  ! Set parameters from simple-structure model to calculate anoxic fratction (Arah and Vinten 1995)
  rij_kro_a     = CNNitrifDenitrifParamsInst%rij_kro_a
  rij_kro_alpha = CNNitrifDenitrifParamsInst%rij_kro_alpha
  rij_kro_beta  = CNNitrifDenitrifParamsInst%rij_kro_beta
  rij_kro_gamma = CNNitrifDenitrifParamsInst%rij_kro_gamma
  rij_kro_delta = CNNitrifDenitrifParamsInst%rij_kro_delta

  organic_max = CNParamsShareInst%organic_max
  

  do j = lbj, ubj
    do fc = 1,numf
      c = filter(fc)
      if(j<jtops(c))cycle
      ! calculate gas diffusivity of soil at field capacity here
      ! use expression from methane code, but neglect OM for now
      f_a = 1._r8 - watfc(c,j) / watsat(c,j)
      eps =  watsat(c,j)-watfc(c,j) ! Air-filled fraction of total soil volume
            
      if (organic_max > 0._r8) then
        om_frac = min(cellorg(c,j)/organic_max, 1._r8)
        ! Use first power, not square as in iniTimeConst
      else
        om_frac = 1._r8
      end if
      diffus (c,j) = (d_con_g(2,1) + d_con_g(2,2)*t_soisno(c,j)) * 1.e-4_r8 * &
                    (om_frac * f_a**(10._r8/3._r8) / watsat(c,j)**2 + &
                    (1._r8-om_frac) * eps**2 * f_a**(3._r8 / bsw(c,j)) ) 

      ! calculate anoxic fraction of soils
      ! use rijtema and kroess model after Riley et al., 2000
      ! caclulated r_psi as a function of psi
      r_min(c,j) = 2 * surface_tension_water / (rho_w * grav * abs(soilpsi(c,j)))
      r_max(c,j) = 2 * surface_tension_water / (rho_w * grav * 0.1_r8)
      r_psi(c,j) = sqrt(r_min(c,j) * r_max(c,j))
      ratio_diffusivity_water_gas(c,j) = (d_con_g(2,1) + d_con_g(2,2)*t_soisno(c,j) ) * 1.e-4_r8 / &
                    ((d_con_w(2,1) + d_con_w(2,2)*t_soisno(c,j) + d_con_w(2,3)*t_soisno(c,j)**2) * 1.e-9_r8)

      if (o2_decomp_depth_unsat(c,j) /= spval .and. conc_o2_unsat(c,j) /= spval .and.  & 
        o2_decomp_depth_unsat(c,j) > 0._r8) then
        anaerobic_frac(c,j) = exp(-rij_kro_a * r_psi(c,j)**(-rij_kro_alpha) * &
                       o2_decomp_depth_unsat(c,j)**(-rij_kro_beta) * &
                       conc_o2_unsat(c,j)**rij_kro_gamma * (h2osoi_vol(c,j) + ratio_diffusivity_water_gas(c,j) * &
                       watsat(c,j))**rij_kro_delta)
      else
        anaerobic_frac(c,j) = 0._r8
      endif

    enddo  
  enddo
  end associate 
  end subroutine calc_anaerobic_frac

!----------------------------------------------------------------------------------------------------
  subroutine calc_potential_aerobic_hr(bounds, lbj, ubj, numf, filter, jtops, nom_pools, pot_decay_rates, pct_sand, pot_co2_hr)
  !
  ! DESCRIPTION
  ! calculate potential aerobic heteorotrophic respiration, and potential oxygen consumption based on cascade_matrix
  use MathfuncMod         , only : dot_sum  
  use CNSharedParamsMod   , only : CNParamsShareInst
  use CNSharedParamsMod   , only : CNParamsShareInst
  
  type(bounds_type),      intent(in) :: bounds
  integer,                intent(in) :: lbj, ubj
  integer,                intent(in) :: jtops(bounds%begc:bounds%endc)        ! top label of each column
  integer,                intent(in) :: numf
  integer,                intent(in) :: filter(:)
  integer,                intent(in) :: nom_pools
  real(r8),               intent(in) :: pot_decay_rates(nom_pools, bounds%begc:bounds%endc, lbj:ubj)  
  real(r8),               intent(in) :: pct_sand(bounds%begc:bounds%endc, lbj:ubj)
  real(r8),               intent(out):: pot_co2_hr(bounds%begc:bounds%endc, lbj:ubj)
  
  real(r8) :: ftxt
  real(r8) :: cascade_matrix_hr(nom_pools)
  integer  :: fc, c, j
  
  cascade_matrix_hr = 0._r8
  !reaction1, lit1 -> som1
  !lit1 + 0.55*o2 -> 0.45 som1 + 0.55co2 + (1/cn_ratios(lit1) - 0.45/cn_ratios(som1))+ (1/cp_ratios(lit1)-0.45/cp_ratios(som1))

  cascade_matrix_hr(1)=  CNDecompBgcParamsInst%rf_l1s1_bgc 
  
  !reaction 2, lit2 -> som1
  !lit2 + 0.5 o2  -> 0.5 som1 + 0.5 co2 + (1/cn_ratios(lit2)-0.5/cn_ratios(som1)) +(1/cp_ratios(lit2)-0.5/cp_ratios(som1))
  cascade_matrix_hr(2) =  CNDecompBgcParamsInst%rf_l2s1_bgc
  
  !reaction 3, lit3 -> som2
  !lit3 + 0.5 o2 -> 0.5 som2 + 0.5 co2
  cascade_matrix_hr(3)=  CNDecompBgcParamsInst%rf_l3s2_bgc
  
  
  !reaction 5, som2 -> som1
  !som2 + 0.55 o2 -> 0.42 som1 + 0.03som3 + 0.55co2 + (1/cn_ratios(som2)-0.42/cn_ratios(som1)-0.03/cn_ratios(som3)) + (1/cp_raitos(som2)-0.42/cp_ratios(som1)-0.03/cp_ratios(som3))
  cascade_matrix_hr(5)   =  CNDecompBgcParamsInst%rf_s2s1_bgc
  
  !reaction 6, s3 -> s1
  !som3 + 0.55 o2 -> 0.45*som1 + 0.55co2 + (1/cn_ratios(som3)-0.45/cn_ratios(som1)) + (1/cp_ratios(som3)-0.45/cp_ratios(som1))
  cascade_matrix_hr(6) = CNDecompBgcParamsInst%rf_s3s1_bgc
  
  !obtain the potential respiration

  do fc = 1, numf
    c = filter(fc)
    do j = jtops(c), ubj
      !reaction 4, the partition into som2 and som3 is soil texture dependent, som1->som2, som3
      !som1 + f(txt) o2 -> f1*som2 + f2*som3 + f(txt) co2 + (1/cn_ratios(som1)-f1/cn_ratios(som2)-f2/cn_ratios(som3)) +(1/cp_ratios(som1)-f1/cp_ratios(som2)-f2/cp_ratios(som3))
      !f(txt) = 0.85_r8 - 0.68_r8 * 0.01_r8 * (100._r8 - sand), assuming sand = 30%
      !f1+f2+f(txt)=1._r8
      ftxt = 0.85_r8 - 0.68_r8 * 0.01_r8 * (100._r8 - pct_sand(c,j))  
      cascade_matrix_hr(4) = ftxt
      
      pot_co2_hr(c,j) = dot_sum(cascade_matrix_hr, pot_decay_rates(:, c, j))  !mol CO2/m3/s
    enddo  
  enddo
  end subroutine calc_potential_aerobic_hr
!----------------------------------------------------------------------------------------------------  
  subroutine calc_decompK_multiply_scalar(bounds, lbj, ubj, numf, filter, jtops, finundated, zsoi, &
    t_soisno, o2_bulk, o2_aqu2bulkcef, soilstate_vars, centurybgc_vars)
  !
  ! DESCRIPTION
  ! compute scalar multipliers for aerobic om decomposition
  ! because temperature and moisture scalars will not be independent from each other for microbe explicit models,
  ! I decide to put temp_scalar and moist_scalar as private variables
  
  ! uses
  !
  use CNSharedParamsMod   , only : CNParamsShareInst
  use shr_const_mod       , only : SHR_CONST_TKFRZ, SHR_CONST_PI
  use SoilStatetype       , only : soilstate_type

  type(bounds_type),         intent(in) :: bounds
  integer,                   intent(in) :: lbj, ubj
  integer,                   intent(in) :: jtops(bounds%begc:bounds%endc)        ! top label of each column
  integer,                   intent(in) :: numf
  integer,                   intent(in) :: filter(:)
  real(r8),                  intent(in) :: t_soisno(bounds%begc:bounds%endc, lbj:ubj)
  real(r8),                  intent(in) :: o2_bulk(bounds%begc:bounds%endc,lbj:ubj)
  real(r8),                  intent(in) :: zsoi(bounds%begc:bounds%endc, lbj:ubj)
  real(r8),                  intent(in) :: finundated(bounds%begc:bounds%endc)  
  real(r8),                  intent(in) :: o2_aqu2bulkcef(bounds%begc:bounds%endc, lbj:ubj)
  type(soilstate_type),      intent(in) :: soilstate_vars
  type(centurybgc_type),  intent(inout) :: centurybgc_vars  
  
  integer :: fc, c, j  !indices
  real(r8), parameter :: normalization_tref = 15._r8        ! reference temperature for normalizaion (degrees C)
  
  real(r8) :: decomp_depth_efolding      ! [m] a testing parameter, which will be replaced,  
  real(r8) :: Q10                                           ! a number taken from CLM4.5bgc
  real(r8) :: froz_q10
  
  !local variables
  real(r8) :: minpsi
  real(r8) :: maxpsi
  real(r8) :: normalization_factor
  real(r8) :: catanf_30
  real(r8) :: catanf
  real(r8) :: t1
  real(r8) :: o2w
  real(r8) :: psi
  !----- CENTURY T response function
  catanf(t1) = 11.75_r8 +(29.7_r8 / SHR_CONST_PI) * atan( SHR_CONST_PI * 0.031_r8  * ( t1 - 15.4_r8 ))
    
  associate(                                              &
    sucsat         => soilstate_vars%sucsat_col         , & ! Input:  [real(r8) (:,:)] minimum soil suction [mm]
    soilpsi        => soilstate_vars%soilpsi_col        , & ! Input:  [real(r8) (:,:)] soilwater pontential in each soil layer [MPa]
    t_scalar       => centurybgc_vars%t_scalar_col      , & ! Output: [real(r8) (:,:)   ]  soil temperature scalar for decomp                     
    w_scalar       => centurybgc_vars%w_scalar_col      , & ! Output: [real(r8) (:,:)   ]  soil water scalar for decomp                           
    o_scalar       => centurybgc_vars%o_scalar_col      , & ! Output: [real(r8) (:,:)   ]  fraction by which decomposition is limited by anoxia   
    depth_scalar   => centurybgc_vars%depth_scalar_col    & ! Output: [real(r8) (:,:,:) ]  rate constant for decomposition (1./sec)             
  )
  catanf_30 = catanf(30._r8)
  
  ! set "Q10" parameter
  Q10 = CNParamsShareInst%Q10

  ! set "froz_q10" parameter
  froz_q10  = CNParamsShareInst%froz_q10

  ! Set "decomp_depth_efolding" parameter
  decomp_depth_efolding = CNParamsShareInst%decomp_depth_efolding
      
  do fc = 1, numf
    c = filter(fc)
    do j = jtops(c), ubj
      !temperature scalar
      t_scalar(c,j)     = 1._r8
      !use Charlie's Q10 based temperature scalar
      if (t_soisno(c,j) >= SHR_CONST_TKFRZ) then
        t_scalar(c,j)= (Q10**((t_soisno(c,j)-(SHR_CONST_TKFRZ+25._r8))/10._r8))
      else
        t_scalar(c,j)= (Q10**(-25._r8/10._r8))*(froz_q10**((t_soisno(c,j)-SHR_CONST_TKFRZ)/10._r8))
      endif
      ! scale all decomposition rates by a constant to compensate for offset between original CENTURY temp func and Q10
      normalization_factor = (catanf(normalization_tref)/catanf_30) / (q10**((normalization_tref-25._r8)/10._r8))
      t_scalar(c,j) = t_scalar(c,j) * normalization_factor
      !moisture scalar, also follows what Charlie has done
      minpsi = -10.0_r8
      w_scalar(c,j)     = 1._r8
      
      maxpsi = sucsat(c,j) * (-9.8e-6_r8)
      psi = min(soilpsi(c,j),maxpsi)
      ! decomp only if soilpsi is higher than minpsi
      if (psi > minpsi) then
        w_scalar(c,j) = (log(minpsi/psi)/log(minpsi/maxpsi))
      else
        w_scalar(c,j) = 0._r8
      end if
            
      !oxygen scalar, this is different from what CLM4.5bgc does, I use a M-M formulation to indicate O2 stress
      !and the O2 budget is done on the fly
      o2w = o2_bulk(c,j) / o2_aqu2bulkcef(c, j)   
      o_scalar(c,j)     = o2w/(o2w++0.22_r8)   !the value 0.22 mol O3/m3 is from Arah and Kirk, 2000
      
      !depth scalar, according to Koven et al. (2013), BG, the depth scalar is needed to resolve the radiocarbon profile
      depth_scalar(c,j) = exp(-zsoi(c,j)/decomp_depth_efolding)
    enddo
  enddo
  
  end associate
  end subroutine calc_decompK_multiply_scalar


  
  !-----------------------------------------------------------------------  
  subroutine calc_nuptake_prof(bounds, nlevdecomp, num_soilc, filter_soilc, sminn_nh4_vr, sminn_no3_vr, &
     dzsoi, nfixation_prof, nuptake_prof)
  !
  !DESCRIPTION
  ! calculate the nitrogen uptake profile
  !

  type(bounds_type)        , intent(in)   :: bounds  
  integer                  , intent(in)   :: nlevdecomp                        ! number of vertical layers
  integer                  , intent(in)   :: num_soilc                         ! number of soil columns in filter
  integer                  , intent(in)   :: filter_soilc(:)                   ! filter for soil columns
  real(r8)                 , intent(in)   :: sminn_nh4_vr(bounds%begc: , 1: )                        ! soil mineral nitrogen profile
  real(r8)                 , intent(in)   :: sminn_no3_vr(bounds%begc: , 1: )                        ! soil mineral nitrogen profile  
  real(r8)                 , intent(in)   :: dzsoi(bounds%begc: , 1: )                                   ! layer thickness
  real(r8)                 , intent(in)   :: nfixation_prof(bounds%begc: , 1: )                  ! nitrogen fixation profile
  real(r8)                 , intent(inout):: nuptake_prof(bounds%begc: , 1: ) ! nitrogen uptake profile
  
  !local variables
  integer :: fc, j, c      ! indices
  real(r8):: sminn_tot(bounds%begc:bounds%endc)  !vertically integrated mineral nitrogen
  
 
  SHR_ASSERT_ALL((ubound(dzsoi_decomp)     == (/nlevdecomp/)), errMsg(__FILE__, __LINE__))   
  SHR_ASSERT_ALL((ubound(sminn_nh4_vr)     == (/bounds%endc, nlevdecomp/)), errMsg(__FILE__, __LINE__))
  SHR_ASSERT_ALL((ubound(sminn_no3_vr)     == (/bounds%endc, nlevdecomp/)), errMsg(__FILE__, __LINE__))
  SHR_ASSERT_ALL((ubound(nfixation_prof)     == (/bounds%endc, nlevdecomp/)), errMsg(__FILE__, __LINE__))
  SHR_ASSERT_ALL((ubound(dzsoi)     == (/bounds%endc, nlevdecomp/)), errMsg(__FILE__, __LINE__))  
  SHR_ASSERT_ALL((ubound(nuptake_prof)     == (/bounds%endc, nlevdecomp/)), errMsg(__FILE__, __LINE__))
  
  ! init sminn_tot
  do fc=1,num_soilc
    c = filter_soilc(fc)
    sminn_tot(c) = 0.
  end do

  do j = 1, nlevdecomp
    do fc=1,num_soilc
      c = filter_soilc(fc)
      sminn_tot(c) = sminn_tot(c) + (sminn_nh4_vr(c,j)+sminn_no3_vr(c,j)) * dzsoi(c,j)
    end do
  end do

  do j = 1, nlevdecomp
    do fc=1,num_soilc
      c = filter_soilc(fc)      
      if (sminn_tot(c)  >  0.) then
        nuptake_prof(c,j) = (sminn_nh4_vr(c,j)+sminn_no3_vr(c,j)) / sminn_tot(c)
      else
        nuptake_prof(c,j) = nfixation_prof(c,j)
      endif

   end do
  end do
         
  
  end subroutine calc_nuptake_prof
  
 
  !-----------------------------------------------------------------------  
  subroutine calc_plant_nitrogen_uptake_prof(bounds, nlevdecomp, num_soilc, filter_soilc, &
    dzsoi, plant_totn_demand_flx_col, nuptake_prof, plant_demand_vr)
  
  !caluate depth specific demand
  use clm_varcon           , only : natomw

  type(bounds_type)        , intent(in)    :: bounds  
  integer                  , intent(in)    :: nlevdecomp                        ! number of vertical layers
  integer                  , intent(in)    :: num_soilc                         ! number of soil columns in filter
  integer                  , intent(in)    :: filter_soilc(:)                   ! filter for soil columns
  real(r8)                 , intent(in)    :: dzsoi(bounds%begc:bounds%endc,1:nlevdecomp)                                   ! layer thickness  
  real(r8)                 , intent(in)    :: plant_totn_demand_flx_col(bounds%begc:bounds%endc)
  real(r8)                 , intent(in)    :: nuptake_prof(bounds%begc:bounds%endc, 1:nlevdecomp)
  real(r8)                 , intent(inout) :: plant_demand_vr(1,bounds%begc:bounds%endc, 1:nlevdecomp)      !mol N/m3/s
  
  integer :: fc, c, j
  
  do j = 1, nlevdecomp
    do fc = 1, num_soilc  
      c = filter_soilc(fc)
      plant_demand_vr(1,c,j) = plant_totn_demand_flx_col(c) * nuptake_prof(c,j) / dzsoi(c,j) /natomw   
    enddo
  enddo  
  
  
  end subroutine calc_plant_nitrogen_uptake_prof
  
 
  !-----------------------------------------------------------------------  
  
  subroutine calc_extneral_bgc_input(bounds, lbj, ubj, num_soilc, filter_soilc, carbonflux_vars, nitrogenflux_vars, &
    centurybgc_vars, betrtracer_vars, tracerstate_vars, cn_ratios, cp_ratios)


  use MathfuncMod         , only :  safe_div

  type(bounds_type)                  , intent(in) :: bounds                             ! bounds
  integer                            , intent(in) :: num_soilc                               ! number of columns in column filter
  integer                            , intent(in) :: filter_soilc(:)                          ! column filter
  type(carbonflux_type)              , intent(in) :: carbonflux_vars
  type(nitrogenflux_type)            , intent(in) :: nitrogenflux_vars
  type(centurybgc_type)              , intent(in) :: centurybgc_vars  
  type(tracerstate_type)             , intent(inout) :: tracerstate_vars
  real(r8)                           , intent(inout) :: cn_ratios(centurybgc_vars%nom_pools, bounds%begc:bounds%endc, lbj:ubj)
  real(r8)                           , intent(inout) :: cp_ratios(centurybgc_vars%nom_pools, bounds%begc:bounds%endc, lbj:ubj)
  

  integer :: k, fc, c, j
 
  associate(                                                                         & !  
    lit1           => centurybgc_vars%lit1                                         , & !
    lit2           => centurybgc_vars%lit2                                         , & !
    lit3           => centurybgc_vars%lit3                                         , & !
    som1           => centurybgc_vars%som1                                         , & !
    som2           => centurybgc_vars%som2                                         , & !
    som3           => centurybgc_vars%som3                                         , & !
    cwd            => centurybgc_vars%cwd                                          , & !
    nelm           => centurybgc_vars%nelms                                        , & !
    id_trc_nh3x    => betrtracer_vars%id_trc_nh3x                                  , & !
    id_trc_no3x    => betrtracer_vars%id_trc_no3x                                  , & !   
    ngwmobile_tracers => betrtracer_vars%ngwmobile_tracers                         , & !
    tracer_conc_solid_passive => tracerstate_vars%tracer_conc_solid_passive_col    , & !
    tracer_conc_mobile => tracerstate_vars%tracer_conc_mobile_col                  , & !
    bgc_cpool_inputs_vr => carbonflux_vars%bgc_cpool_inputs_vr_col                 , & !
    bgc_npool_inputs_vr => nitrogenflux_vars%bgc_npool_inputs_vr_col               , & !
    sminn_nh4_input_vr  => nitrogenflux_vars%sminn_nh4_input_vr_col                , & !
    sminn_no3_input_vr  => nitrogenflux_vars%sminn_no3_input_vr_col                  & 
  )
  

  do fc = 1, num_soilc
    c = filter_soilc(fc)
    
    do j = lbj, ubj
      k = lit1
      tracer_conc_solid_passive(c,j,k*nelm-1) = tracer_conc_solid_passive(c,j,k*nelm-1) + bgc_cpool_inputs_vr(c,n,k)/catomw      
      tracer_conc_solid_passive(c,j,k*nelm)   = tracer_conc_solid_passive(c,j,k*nelm) +   bgc_npool_inputs_vr(c,n,k)/natomw
      
      cn_ratios(k, c,j) = safe_div(tracer_conc_solid_passive(c,j,k*nelm-1), tracer_conc_solid_passive(c,j,k*nelm))

      k = lit2
      tracer_conc_solid_passive(c,j,k*nelm-1) = tracer_conc_solid_passive(c,j,k*nelm-1) + bgc_cpool_inputs_vr(c,n,k)/catomw      
      tracer_conc_solid_passive(c,j,k*nelm)   = tracer_conc_solid_passive(c,j,k*nelm) +   bgc_npool_inputs_vr(c,n,k)/natomw
      
      cn_ratios(k, c,j) = safe_div(tracer_conc_solid_passive(c,j,k*nelm-1), tracer_conc_solid_passive(c,j,k*nelm))
      
      k = lit3
      tracer_conc_solid_passive(c,j,k*nelm-1) = tracer_conc_solid_passive(c,j,k*nelm-1) + bgc_cpool_inputs_vr(c,n,k)/catomw      
      tracer_conc_solid_passive(c,j,k*nelm)   = tracer_conc_solid_passive(c,j,k*nelm) +   bgc_npool_inputs_vr(c,n,k)/natomw
      
      cn_ratios(k,c,j) = safe_div(tracer_conc_solid_passive(c,j,k*nelm-1), tracer_conc_solid_passive(c,j,k*nelm))

      k = cwd
      tracer_conc_solid_passive(c,j,k*nelm-1) = tracer_conc_solid_passive(c,j,k*nelm-1) + bgc_cpool_inputs_vr(c,n,k)/catomw      
      tracer_conc_solid_passive(c,j,k*nelm)   = tracer_conc_solid_passive(c,j,k*nelm) +   bgc_npool_inputs_vr(c,n,k)/natomw
      
      cn_ratios(k, c,j) = safe_div(tracer_conc_solid_passive(c,j,k*nelm-1), tracer_conc_solid_passive(c,j,k*nelm))
      
      k = som1
      cn_ratios(k, c,j) = safe_div(tracer_conc_solid_passive(c,j,k*nelm-1), tracer_conc_solid_passive(c,j,k*nelm))

      k = som2
      cn_ratios(k, c,j) = safe_div(tracer_conc_solid_passive(c,j,k*nelm-1), tracer_conc_solid_passive(c,j,k*nelm))

      k = som3
      cn_ratios(k, c,j) = safe_div(tracer_conc_solid_passive(c,j,k*nelm-1), tracer_conc_solid_passive(c,j,k*nelm))

      tracer_conc_mobile(c, j, id_trc_nh3x) = tracer_conc_mobile(c, j, id_trc_nh3x) + sminn_nh4_input_vr_col(c,j)/natomw
      tracer_conc_mobile(c, j, id_trc_no3x) = tracer_conc_mobile(c, j, id_trc_no3x) + sminn_no3_input_vr_col(c,j)/natomw
      
    enddo
  enddo
  
  end associate
  end subroutine calc_extneral_bgc_input
  

end module BGCCenturySubMod
