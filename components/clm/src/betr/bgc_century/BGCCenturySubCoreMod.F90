module BGCCenturySubCoreMod
#include "shr_assert.h"
  !
  ! !DESCRIPTION:
  ! module contains subroutine for the century bgc implementation
  ! !USES:
  use shr_kind_mod       , only : r8 => shr_kind_r8
  use shr_log_mod        , only : errMsg => shr_log_errMsg
  use decompMod          , only : bounds_type
  use clm_varcon         , only : spval
  use clm_varcon         , only : catomw, natomw
  use clm_varpar         , only : ndecomp_pools
  use ColumnType         , only : col
  use clm_varctl         , only : spinup_state,iulog
  use clm_varctl         , only : CNAllocate_Carbon_only
  implicit none
  save
  private

  public :: apply_plant_root_respiration_prof
  public :: apply_plant_root_nuptake_prof
  public :: calc_som_deacyK
  public :: calc_sompool_decay
  public :: init_state_vector
  public :: retrieve_flux_vars
  public :: retrieve_state_vars
  public :: calc_nitrif_denitrif_rate
  public :: calc_anaerobic_frac
  public :: calc_potential_aerobic_hr
  public :: calc_decompK_multiply_scalar
  public :: calc_nuptake_prof
  public :: calc_plant_nitrogen_uptake_prof
  public :: bgcstate_ext_update_bfdecomp
  public :: bgcstate_ext_update_afdecomp
  public :: set_reaction_order
  public :: calc_nutrient_compet_rescal
  public :: assign_nitrogen_hydroloss
  public :: assign_OM_CNpools

  type, public :: centurybgc_type

     integer           :: nom_pools                              !not include coarse wood debris
     integer           :: nom_totelms
     integer           :: lit1, lit1_dek_reac
     integer           :: lit2, lit2_dek_reac
     integer           :: lit3, lit3_dek_reac
     integer           :: som1, som1_dek_reac
     integer           :: som2, som2_dek_reac
     integer           :: som3, som3_dek_reac
     integer           :: cwd,  cwd_dek_reac

     integer           :: c_loc
     integer           :: n_loc
     integer           :: nelms                                  !number of chemical elements in an om pool

                                                                 !reactive primary variables
     real(r8)          :: k_decay_lit1
     real(r8)          :: k_decay_lit2
     real(r8)          :: k_decay_lit3
     real(r8)          :: k_decay_som1
     real(r8)          :: k_decay_som2
     real(r8)          :: k_decay_som3
     real(r8)          :: k_decay_cwd
     integer           :: lid_nh4, lid_nh4_nit_reac              !local position of nh4 in the state variable vector
     integer           :: lid_no3, lid_no3_den_reac              !local position of no3 in the state variable vector
     integer           :: lid_plant_minn, lid_plant_minn_up_reac !local position of plant consumption of mineral nitrogen in the state variable vector
     integer           :: lid_at_rt, lid_at_rt_reac              !root autotrophic respiration

                                                                 !non reactive primary variables
     integer           :: lid_ar, lid_ar_aere_reac               !local position of ar in the state variable vector
     integer           :: lid_ch4, lid_ch4_aere_reac             !nonreactive primary variables

                                                                 !secondary variables
     integer           :: lid_o2,  lid_o2_aere_reac              !local position of o2 in the state variable vector
     integer           :: lid_co2, lid_co2_aere_reac             !local position of co2 in the state variable vector
     integer           :: lid_n2,  lid_n2_aere_reac
     integer           :: lid_n2o, lid_n2o_aere_reac
                                                                 !diagnostic variables
     integer           :: lid_n2o_nit                            !n2o production from nitrification, used to for mass balance book keeping
     integer           :: lid_co2_hr                             !co2 production from heterotrophic respiration
     integer           :: lid_no3_den                            !no3 consumption due to denitrification
     integer           :: lid_minn_nh4_immob                     !net mineral N immobilization for decomposition
     integer           :: lid_minn_no3_immob
     integer           :: lid_minn_nh4_plant
     integer           :: lid_minn_no3_plant
     integer           :: lid_nh4_nit
                                                                 !aerechyma transport, diagnostic efflux

     integer           :: lid_ar_paere
     integer           :: lid_n2_paere
     integer           :: lid_o2_paere
     integer           :: lid_co2_paere
     integer           :: lid_ch4_paere
     integer           :: lid_n2o_paere

     integer           :: lid_nh4_mpcbuf                         !nh4 buffer providing the tight to loose coupling between immobilizing and mineralizing decomposers
     integer           :: lid_nh4_supp
     integer           :: nstvars                                !number of equations for the state variabile vector
     integer           :: nprimvars                              !total number of primary variables
     integer           :: nreactions                             !seven decomposition pathways plus nitrification, denitrification and plant immobilization
     integer           :: ncompets                               !decomposers, + nitrifiers, + denitrifiers, + plants, + adsorption surface

     integer           :: lid_lit1_compet
     integer           :: lid_lit2_compet
     integer           :: lid_lit3_compet
     integer           :: lid_cwd_compet
     integer           :: lid_som1_compet
     integer           :: lid_som2_compet
     integer           :: lid_som3_compet
     integer           :: lid_plant_compet
     integer           :: lid_nitri_compet
     integer           :: lid_denit_compet
     integer           :: lid_clay_compet

     real(r8), pointer :: t_scalar_col(:,:)
     real(r8), pointer :: w_scalar_col(:,:)
     real(r8), pointer :: o_scalar_col(:,:)
     real(r8), pointer :: depth_scalar_col(:,:)                  !
     integer , pointer :: primvarid(:)
     logical , pointer :: is_aerobic_reac(:)

   contains
     procedure, public  :: Init
     procedure, private :: Init_pars
     procedure, private :: InitAllocate
  end type centurybgc_type

contains

  subroutine Init(this, bounds, lbj, ubj, do_mpc)
    !
    ! DESCRIPTION:
    ! Initialize centurybgc type
    ! !USES:
    use ncdio_pio    , only: file_desc_t

    ! !ARGUMENTS:
    class(centurybgc_type)            :: this
    type(bounds_type),     intent(in) :: bounds
    integer,               intent(in) :: lbj, ubj                           ! lower and upper bounds, make sure they are > 0

    ! !LOCAL VARIABLES:
    logical, optional :: do_mpc
    logical :: do_mpc_loc

    if(present(do_mpc))then
       do_mpc_loc = do_mpc
    else
       do_mpc_loc = .false.
    endif
    call this%init_pars(do_mpc_loc)

    call this%InitAllocate(bounds, lbj, ubj)

  end subroutine Init
  !-------------------------------------------------------------------------------

  subroutine Init_pars(this, do_mpc)
    !
    ! !DESCRIPTION:
    !  describe the layout of the stoichiometric matrix for the reactions
    !           r{1} r{2} r{3} r{4} ... r{n}
    ! s{1}
    ! s{2}
    ! s{3}
    ! s{4}
    ! ...
    ! s{n}
    ! s{n+1}  nonreactive primary variables
    ! s{n+2}
    ! ...
    ! s{m}
    ! s{m+1} diagnostic variables
    ! s{p}
    ! each reaction is associated with a primary species, the secondary species follows after primary species
    ! for the century model, the primary species are seven om pools and nh4, no3 and plant nitrogen
    !
    ! !USES:
    use MathfuncMod            , only : addone
    class(centurybgc_type) :: this
    logical   , intent(in) :: do_mpc

    ! !LOCAL VARIABLES:
    integer :: itemp
    integer :: ireac   !counter of reactions
    integer :: itemp1

    itemp = 0
    ireac = 0
    this%nom_pools = 7   !not include coarse wood debris
    this%lit1 = addone(itemp); this%lit1_dek_reac = addone(ireac)
    this%lit2 = addone(itemp); this%lit2_dek_reac = addone(ireac)
    this%lit3 = addone(itemp); this%lit3_dek_reac = addone(ireac)
    this%cwd  = addone(itemp); this%cwd_dek_reac  = addone(ireac)
    this%som1 = addone(itemp); this%som1_dek_reac = addone(ireac)
    this%som2 = addone(itemp); this%som2_dek_reac = addone(ireac)
    this%som3 = addone(itemp); this%som3_dek_reac = addone(ireac)

    this%nelms = 2   !carbon and nitrogen
    this%c_loc = 1
    this%n_loc = 2

    itemp               = this%nom_pools*this%nelms
    this%nom_totelms    = itemp
    this%lid_nh4        = addone(itemp); this%lid_nh4_nit_reac = addone(ireac)       !this is also used to indicate the nitrification reaction
    this%lid_no3        = addone(itemp); this%lid_no3_den_reac = addone(ireac)       !this is also used to indicate the denitrification reaction
    this%lid_plant_minn = addone(itemp); this%lid_plant_minn_up_reac = addone(ireac) !this is used to indicate plant mineral nitrogen uptake
    this%lid_at_rt      = addone(itemp); this%lid_at_rt_reac = addone(ireac)         !this is used to indicate plant autotrophic root respiration

    !non-reactive primary variables
    this%lid_ch4        = addone(itemp);
    this%lid_ar         = addone(itemp);


    !second primary variables
    this%lid_o2         = addone(itemp);
    this%lid_co2        = addone(itemp);
    this%lid_n2o        = addone(itemp);
    this%lid_n2         = addone(itemp);

    this%lid_o2_aere_reac  = addone(ireac)

    if(do_mpc)then
       this%lid_nh4_mpcbuf = addone(itemp)
    endif
    this%nprimvars      = itemp     !primary state variables 14 + 6

    !diagnostic variables
    this%lid_n2o_nit        = addone(itemp)
    this%lid_co2_hr         = addone(itemp)
    this%lid_no3_den        = addone(itemp)
    this%lid_minn_nh4_immob = addone(itemp)
    this%lid_minn_no3_immob = addone(itemp)
    this%lid_minn_nh4_plant = addone(itemp)
    this%lid_minn_no3_plant = addone(itemp)
    this%lid_nh4_nit        = addone(itemp)

    if(CNAllocate_Carbon_only())then
       this%lid_nh4_supp = addone(itemp)
    endif
    !aerechyma transport
    this%lid_o2_paere   = addone(itemp)   !
    if ( spinup_state /= 1 ) then
       this%lid_ar_paere   = addone(itemp);  this%lid_ar_aere_reac  = addone(ireac)   !
       this%lid_n2_paere   = addone(itemp);  this%lid_n2_aere_reac  = addone(ireac)   !
       this%lid_co2_paere  = addone(itemp);  this%lid_co2_aere_reac = addone(ireac)   !
       this%lid_ch4_paere  = addone(itemp);  this%lid_ch4_aere_reac = addone(ireac)   !
       this%lid_n2o_paere  = addone(itemp);  this%lid_n2o_aere_reac = addone(ireac)   !
    endif
    this%nstvars          = itemp          !totally 14+32 state variables

    this%nreactions = ireac            !seven decomposition pathways plus root auto respiration, nitrification, denitrification and plant immobilization
    allocate(this%primvarid(ireac)); this%primvarid(:) = -1
    allocate(this%is_aerobic_reac(ireac)); this%is_aerobic_reac(:)=.false.

    !decomposers, + nitrifiers, + denitrifiers, + plants, + adsorption surface
    itemp1=0
    this%lid_lit1_compet      = addone(itemp1)
    this%lid_lit2_compet      = addone(itemp1)
    this%lid_lit3_compet      = addone(itemp1)
    this%lid_cwd_compet       = addone(itemp1)
    this%lid_som1_compet      = addone(itemp1)
    this%lid_som2_compet      = addone(itemp1)
    this%lid_som3_compet      = addone(itemp1)
    this%lid_plant_compet     = addone(itemp1)
    this%lid_nitri_compet     = addone(itemp1)
    this%lid_denit_compet     = addone(itemp1)
    this%lid_clay_compet      = addone(itemp1)

    this%ncompets  = itemp1

  end subroutine Init_pars
  !-------------------------------------------------------------------------------

  subroutine InitAllocate(this, bounds, lbj, ubj)
    !
    ! !DESCRIPTION:
    ! memory allocation for the data type specified by this
    !
    ! !ARGUMENTS:
    class(centurybgc_type)         :: this
    type(bounds_type) , intent(in) :: bounds
    integer           , intent(in) :: lbj, ubj                           ! lower and upper bounds, make sure they are > 0


    allocate(this%t_scalar_col(bounds%begc:bounds%endc,     lbj:ubj))
    allocate(this%w_scalar_col(bounds%begc:bounds%endc,     lbj:ubj))
    allocate(this%o_scalar_col(bounds%begc:bounds%endc,     lbj:ubj))
    allocate(this%depth_scalar_col(bounds%begc:bounds%endc, lbj:ubj))

  end subroutine InitAllocate

  !-------------------------------------------------------------------------------
  subroutine init_state_vector(bounds, lbj, ubj, numf, filter, jtops, neq, &
       tracerstate_vars, betrtracer_vars, centurybgc_vars, y0)
    !
    ! !DESCRIPTION:
    ! number of equations, total number of carbon pools + o2 + co2
    !
    ! !USES:
    use tracerstatetype       , only : tracerstate_type
    use BeTRTracerType        , only : betrtracer_type
    ! !ARGUMENTS:
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

    ! !LOCAL VARIABLES:
    integer :: fc, c, j

    ! all organic matter pools are distributed into solid passive tracers
    do j = lbj, ubj
       do fc = 1, numf
          c = filter(fc)
          if(j>=jtops(c))then
             !zero out everything
             y0(:, c, j ) = 0._r8

             !set up nonzero variables
             y0(1:centurybgc_vars%nom_pools*centurybgc_vars%nelms, c, j)    = tracerstate_vars%tracer_conc_solid_passive_col(c, j, :)

             y0(centurybgc_vars%lid_n2,  c, j)        = max(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_n2)  ,0._r8)

             y0(centurybgc_vars%lid_o2,  c, j)        = max(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_o2)  ,0._r8)

             y0(centurybgc_vars%lid_ar,  c, j)        = max(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_ar)  ,0._r8)

             y0(centurybgc_vars%lid_co2, c, j)        = max(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_co2x),0._r8)

             y0(centurybgc_vars%lid_ch4, c, j)        = max(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_ch4) ,0._r8)

             y0(centurybgc_vars%lid_nh4, c, j)        = max(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_nh3x),0._r8)

             y0(centurybgc_vars%lid_no3, c, j)        = max(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_no3x),0._r8)

             y0(centurybgc_vars%lid_n2o, c, j)        = max(tracerstate_vars%tracer_conc_mobile_col(c, j, betrtracer_vars%id_trc_n2o), 0._r8)
          endif
       enddo
    enddo
  end subroutine init_state_vector

  !-------------------------------------------------------------------------------
  subroutine calc_som_deacyK(bounds, lbj, ubj, numf, filter, jtops, nom_pools, &
       tracercoeff_vars, tracerstate_vars, betrtracer_vars, centurybgc_vars, &
       carbonflux_vars, dtime, k_decay)
    !
    ! !DESCRIPTION:
    ! calculate decay coefficients for different pools
    !
    ! !USES:
    use tracercoeffType    , only          : tracercoeff_type
    use BetrTracerType     , only          : betrtracer_type
    use tracerstatetype    , only          : tracerstate_type
    use BeTRTracerType     , only          : betrtracer_type
    use BGCCenturyParMod   , only          : CNDecompBgcParamsInst
    use CNCarbonFluxType   , only          : carbonflux_type
    integer                , intent(in)    :: nom_pools
    type(bounds_type)      , intent(in)    :: bounds
    integer                , intent(in)    :: lbj, ubj
    integer                , intent(in)    :: jtops(bounds%begc:bounds%endc)        ! top label of each column
    integer                , intent(in)    :: numf
    integer                , intent(in)    :: filter(:)
    real(r8)               , intent(in)    :: dtime
    type(betrtracer_type)  , intent(in)    :: betrtracer_vars                    ! betr configuration information
    type(centurybgc_type)  , intent(in)    :: centurybgc_vars
    type(carbonflux_type)  , intent(in)    :: carbonflux_vars
    type(tracercoeff_type) , intent(in)    :: tracercoeff_vars
    type(tracerstate_type) , intent(inout) :: tracerstate_vars
    real(r8)               , intent(out)   :: k_decay(nom_pools, bounds%begc:bounds%endc, lbj:ubj)
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c, j
    real(r8):: dtimei
    associate(                                                  &
         t_scalar       => carbonflux_vars%t_scalar_col       , & ! Output: [real(r8) (:,:)   ]  soil temperature scalar for decomp
         w_scalar       => carbonflux_vars%w_scalar_col       , & ! Output: [real(r8) (:,:)   ]  soil water scalar for decomp
         o_scalar       => carbonflux_vars%o_scalar_col       , & ! Output: [real(r8) (:,:)   ]  fraction by which decomposition is limited by anoxia
         depth_scalar   => centurybgc_vars%depth_scalar_col   , & ! Output: [real(r8) (:,:,:) ]  rate constant for decomposition (1./sec)
         lit1           => centurybgc_vars%lit1               , & !
         lit2           => centurybgc_vars%lit2               , & !
         lit3           => centurybgc_vars%lit3               , & !
         som1           => centurybgc_vars%som1               , & !
         som2           => centurybgc_vars%som2               , & !
         som3           => centurybgc_vars%som3               , & !
         cwd            => centurybgc_vars%cwd                , & !
         k_decay_lit1   => CNDecompBgcParamsInst%k_decay_lit1 , & !
         k_decay_lit2   => CNDecompBgcParamsInst%k_decay_lit2 , & !
         k_decay_lit3   => CNDecompBgcParamsInst%k_decay_lit3 , & !
         k_decay_som1   => CNDecompBgcParamsInst%k_decay_som1 , & !
         k_decay_som2   => CNDecompBgcParamsInst%k_decay_som2 , & !
         k_decay_som3   => CNDecompBgcParamsInst%k_decay_som3 , & !
         k_decay_cwd    => CNDecompBgcParamsInst%k_decay_cwd    & !
         )

      dtimei=1._r8/dtime
      k_decay(:, :, :) = spval
      do j = lbj, ubj
         do fc = 1, numf
            c = filter(fc)
            if(j>=jtops(c))then
               k_decay(lit1, c, j) = min(k_decay_lit1 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) * depth_scalar(c,j),dtimei)
               k_decay(lit2, c, j) = min(k_decay_lit2 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) * depth_scalar(c,j),dtimei)
               k_decay(lit3, c, j) = min(k_decay_lit3 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) * depth_scalar(c,j),dtimei)
               k_decay(som1, c, j) = min(k_decay_som1 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) * depth_scalar(c,j),dtimei)
               k_decay(som2, c, j) = min(k_decay_som2 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) * depth_scalar(c,j),dtimei)
               k_decay(som3, c, j) = min(k_decay_som3 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) * depth_scalar(c,j),dtimei)
               k_decay(cwd,  c, j) = min(k_decay_cwd  * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) * depth_scalar(c,j),dtimei)
            endif
         enddo
      enddo
    end associate
  end subroutine calc_som_deacyK


  !-------------------------------------------------------------------------------
  subroutine calc_sompool_decay(bounds, lbj, ubj, numf, filter, jtops, &
       centurybgc_vars, k_decay, om_pools, decay_rates)
    !
    ! !DESCRIPTION:
    ! calculate degradation for all different pools
    !
    ! !USES:
    type(centurybgc_type) , intent(in)    :: centurybgc_vars
    type(bounds_type)     , intent(in)    :: bounds
    integer               , intent(in)    :: lbj, ubj
    integer               , intent(in)    :: jtops(bounds%begc:bounds%endc)        ! top label of each column
    integer               , intent(in)    :: numf
    integer               , intent(in)    :: filter(:)
    real(r8)              , intent(inout) :: k_decay(centurybgc_vars%nom_pools, bounds%begc:bounds%endc, lbj:ubj)
    real(r8)              , intent(in)    :: om_pools(centurybgc_vars%nom_totelms,bounds%begc:bounds%endc, lbj:ubj)
    real(r8)              , intent(out)   :: decay_rates(centurybgc_vars%nom_pools, bounds%begc:bounds%endc, lbj:ubj)

    ! !LOCAL VARIABLES:
    integer :: jj, fc, c, j
    integer :: kc, kn
    associate(                                       &
         nelms => centurybgc_vars%nelms            , &
         nom_pools => centurybgc_vars%nom_pools    , &
         nom_totelms => centurybgc_vars%nom_totelms, &
         c_loc => centurybgc_vars%c_loc            , &
         n_loc => centurybgc_vars%n_loc              &
         )

      do j = lbj, ubj
         do fc = 1, numf
            c = filter(fc)
            if(j>=jtops(c))then
               !for om pools
               do jj = 1, nom_pools
                  kc = (jj-1) * nelms + c_loc
                  kn = (jj-1) * nelms + n_loc
                  if(min(om_pools(kc, c, j),om_pools(kn, c, j))<1.e-10_r8)then
                     k_decay(jj,c,j) = 0._r8
                  endif
                  decay_rates(jj, c, j) = om_pools(kc, c, j) * k_decay(jj, c, j)
               enddo
            endif
         enddo
      enddo

      !for nitrification and denitrification
    end associate
  end subroutine calc_sompool_decay

  !-------------------------------------------------------------------------------
  
  subroutine retrieve_flux_vars(bounds, lbj, ubj, numf, filter, jtops, neq, dtime, yf, y0,    &
       centurybgc_vars, betrtracer_vars, tracerflux_vars, carbonflux_vars, nitrogenflux_vars, &
       plantsoilnutrientflux_vars)
    !
    ! !DESCRIPTION:
    !
    ! retrieve the fluxes
    ! !USES:
    use tracerfluxType           , only : tracerflux_type
    use PlantSoilnutrientFluxType, only : plantsoilnutrientflux_type
    use CNCarbonFluxType         , only : carbonflux_type
    use CNNitrogenFluxType       , only : nitrogenflux_type
    use BeTRTracerType           , only : betrtracer_type

    ! !ARGUMENTS:
    type(bounds_type),                  intent(in)    :: bounds
    integer,                            intent(in)    :: lbj, ubj
    integer,                            intent(in)    :: jtops(bounds%begc:bounds%endc)        ! top label of each column
    integer,                            intent(in)    :: numf
    integer,                            intent(in)    :: filter(:)
    integer                          ,  intent(in)    :: neq
    real(r8)                         ,  intent(in)    :: dtime
    real(r8)                         ,  intent(in)    :: yf(neq, bounds%begc:bounds%endc, lbj:ubj) !
    real(r8)                         ,  intent(in)    :: y0(neq, bounds%begc:bounds%endc, lbj:ubj) !
    type(centurybgc_type)            ,  intent(in)    :: centurybgc_vars
    type(betrtracer_type)            ,  intent(in)    :: betrtracer_vars
    type(carbonflux_type)            ,  intent(inout) :: carbonflux_vars
    type(nitrogenflux_type)          ,  intent(inout) :: nitrogenflux_vars
    type(tracerflux_type)            ,  intent(inout) :: tracerflux_vars
    type(plantsoilnutrientflux_type),   intent(inout) :: plantsoilnutrientflux_vars

    ! !LOCAL VARIABLES:
    real(r8) :: deltac, fnit
    real(r8) :: delta_nh4, delta_no3
    real(r8) :: delta_nh4_m,delta_no3_m
    real(r8) :: sminn_plant, sminn_plant2
    real(r8) :: err,hr, immob
    real(r8) :: f_nit_n2o, f_den
    integer :: fc, c, j, k

    associate(                                                                   & !
         nom_pools             => centurybgc_vars%nom_pools                    , & !
         nelms                 => centurybgc_vars%nelms                        , & !
         c_loc                 => centurybgc_vars%c_loc                        , & !
         n_loc                 => centurybgc_vars%n_loc                        , & !
         f_n2o_nit_vr          => nitrogenflux_vars%f_n2o_nit_vr_col           , & !
         f_denit_vr            => nitrogenflux_vars%f_denit_vr_col             , & !
         f_nit_vr              => nitrogenflux_vars%f_nit_vr_col               , & !
         actual_immob_no3_vr   => nitrogenflux_vars%actual_immob_no3_vr_col    , & !
         actual_immob_nh4_vr   => nitrogenflux_vars%actual_immob_nh4_vr_col    , & !
         smin_no3_to_plant_vr  => nitrogenflux_vars%smin_no3_to_plant_vr_col   , & !
         smin_nh4_to_plant_vr  => nitrogenflux_vars%smin_nh4_to_plant_vr_col   , & !
         supplement_to_sminn_vr=> nitrogenflux_vars%supplement_to_sminn_vr_col , & !
         hr_vr                 => carbonflux_vars%hr_vr_col                    , & !
         volatileid            => betrtracer_vars%volatileid                   , & !
         ngwmobile_tracers     => betrtracer_vars%ngwmobile_tracers            , & !
         tracer_flx_netpro_vr  => tracerflux_vars%tracer_flx_netpro_vr_col     , & !
         tracer_flx_parchm_vr  => tracerflux_vars%tracer_flx_parchm_vr_col       & !
         )

      if(CNAllocate_Carbon_only())then
         do j = lbj, ubj
            do fc = 1, numf
               c = filter(fc)
               supplement_to_sminn_vr(c,j) = (y0(centurybgc_vars%lid_nh4_supp, c, j) - yf(centurybgc_vars%lid_nh4_supp, c, j))*natomw/dtime
            enddo
         enddo
      endif

      do j = lbj, ubj
         do fc = 1, numf
            c = filter(fc)

            if(j>=jtops(c))then
               plantsoilnutrientflux_vars%plant_minn_active_yield_flx_vr_col(c,j) = (yf(centurybgc_vars%lid_plant_minn, c, j) - y0(centurybgc_vars%lid_plant_minn, c, j))*natomw

               smin_no3_to_plant_vr(c,j) = (yf(centurybgc_vars%lid_minn_no3_plant, c, j) - y0(centurybgc_vars%lid_minn_no3_plant, c, j))*natomw/dtime
               smin_nh4_to_plant_vr(c,j) = (yf(centurybgc_vars%lid_minn_nh4_plant, c, j) - y0(centurybgc_vars%lid_minn_nh4_plant, c, j))*natomw/dtime

               hr_vr       (c,j)  = (yf(centurybgc_vars%lid_co2_hr, c, j) - y0(centurybgc_vars%lid_co2_hr, c, j))*catomw/dtime
               f_nit_vr    (c,j)  = (yf(centurybgc_vars%lid_nh4_nit,c, j) - y0(centurybgc_vars%lid_nh4_nit,c, j))*natomw/dtime
               f_n2o_nit_vr(c,j)  = (yf(centurybgc_vars%lid_n2o_nit,c, j) - y0(centurybgc_vars%lid_n2o_nit,c, j))*natomw/dtime
               f_denit_vr  (c,j)  = (yf(centurybgc_vars%lid_no3_den,c, j) - y0(centurybgc_vars%lid_no3_den,c, j))*natomw/dtime

               actual_immob_no3_vr(c,j) = (yf(centurybgc_vars%lid_minn_no3_immob,c, j) - y0(centurybgc_vars%lid_minn_no3_immob,c, j))*natomw/dtime
               actual_immob_nh4_vr(c,j) = (yf(centurybgc_vars%lid_minn_nh4_immob,c, j) - y0(centurybgc_vars%lid_minn_nh4_immob,c, j))*natomw/dtime

               !the temporal averaging for fluxes below will be done later

               tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_o2)  ) = yf(centurybgc_vars%lid_o2_paere  ,c, j)  - y0(centurybgc_vars%lid_o2_paere , c, j)

               if ( spinup_state /= 1 ) then
                  tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_n2)  ) = yf(centurybgc_vars%lid_n2_paere  ,c, j)  - y0(centurybgc_vars%lid_n2_paere , c, j)
                  tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_ar)  ) = yf(centurybgc_vars%lid_ar_paere  ,c, j)  - y0(centurybgc_vars%lid_ar_paere , c, j)
                  tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_co2x)) = yf(centurybgc_vars%lid_co2_paere ,c, j)  - y0(centurybgc_vars%lid_co2_paere, c, j)
                  tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_ch4) ) = yf(centurybgc_vars%lid_ch4_paere ,c, j)  - y0(centurybgc_vars%lid_ch4_paere, c, j)
                  tracer_flx_parchm_vr(c,j,volatileid(betrtracer_vars%id_trc_n2o) ) = yf(centurybgc_vars%lid_n2o_paere ,c, j)  - y0(centurybgc_vars%lid_n2o_paere, c, j)
               endif

               tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_nh3x) =      &
                    tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_nh3x) + &
                    yf(centurybgc_vars%lid_nh4,c,j) - y0(centurybgc_vars%lid_nh4,c,j)                    

               tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_no3x)  =      &
                    tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_no3x)  + &
                    yf(centurybgc_vars%lid_no3,c,j) - y0(centurybgc_vars%lid_no3,c,j)                    

               tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_n2) = &
                    yf(centurybgc_vars%lid_n2,c,j) - y0(centurybgc_vars%lid_n2,c,j)

               tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_co2x ) = &
                    yf(centurybgc_vars%lid_co2,c,j) - y0(centurybgc_vars%lid_co2,c,j)

               tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_n2o  ) = &
                    yf(centurybgc_vars%lid_n2o,c,j) - y0(centurybgc_vars%lid_n2o,c,j)                    

               tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_o2   ) = &
                    yf(centurybgc_vars%lid_o2,c,j) - y0(centurybgc_vars%lid_o2,c,j)                   

               tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_ch4  ) = &
                    yf(centurybgc_vars%lid_ch4,c,j) - y0(centurybgc_vars%lid_ch4,c,j)                   

               tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_ar) = &
                    yf(centurybgc_vars%lid_ar,c,j) - y0(centurybgc_vars%lid_ar,c,j)

               !get net production for om pools
               deltac=0._r8
               do k = 1, nom_pools
                  tracer_flx_netpro_vr(c,j,ngwmobile_tracers+(k-1)*nelms+c_loc) =      &
                       tracer_flx_netpro_vr(c,j,ngwmobile_tracers+(k-1)*nelms+c_loc) + &
                       yf((k-1)*nelms+c_loc, c, j) - y0((k-1)*nelms+c_loc, c, j)
                       
                  tracer_flx_netpro_vr(c,j,ngwmobile_tracers+(k-1)*nelms+n_loc) =      &
                       tracer_flx_netpro_vr(c,j,ngwmobile_tracers+(k-1)*nelms+n_loc) + &
                       yf((k-1)*nelms+n_loc, c, j) - y0((k-1)*nelms+n_loc, c, j)                       
               enddo
            endif
         enddo

      enddo

    end associate
  end subroutine retrieve_flux_vars
  !-------------------------------------------------------------------------------

  subroutine retrieve_state_vars(bounds, lbj, ubj, numf, filter, jtops, neq, yf, &
       centurybgc_vars, betrtracer_vars, tracerstate_vars)
    !
    ! !DESCRIPTION:
    ! retrieve state variables
    ! number of equations, total number of carbon pools + o2 + co2
    ! !USES:
    use tracerstatetype       , only : tracerstate_type
    use BeTRTracerType        , only : betrtracer_type
    use MathfuncMod           , only : addone

    ! !ARGUMENTS:
    type(bounds_type)       , intent(in)    :: bounds
    integer                 , intent(in)    :: lbj, ubj
    integer                 , intent(in)    :: jtops(bounds%begc:bounds%endc)            ! top label of each column
    integer                 , intent(in)    :: numf
    integer                 , intent(in)    :: filter(:)
    integer                 , intent(in)    :: neq
    type(betrtracer_type)   , intent(in)    :: betrtracer_vars                           ! betr configuration information
    type(centurybgc_type)   , intent(in)    :: centurybgc_vars
    real(r8)                , intent(in)    :: yf(neq, bounds%begc:bounds%endc, lbj:ubj) !
    type(tracerstate_type)  , intent(inout) :: tracerstate_vars

    ! !LOCAL VARIABLES:
    integer :: fc, c, j, k1, k2, k, ll, l
    real(r8):: totsomc
    ! all organic matter pools are distributed into solid passive tracers
    associate(   &
         ngwtracers             => betrtracer_vars%ngwmobile_tracers       , & !
         tracer_conc_mobile_col => tracerstate_vars%tracer_conc_mobile_col   &
         )

      !only retrieve non-mobile tracers
      
      do k = 1, ndecomp_pools
         do l = 1, centurybgc_vars%nelms
            ll = (k-1)*centurybgc_vars%nelms + l
            do j = lbj, ubj   !currently, om is added only to soil layers
               do fc = 1, numf
                  c = filter(fc)
                  if(j>=jtops(c))then
                     tracerstate_vars%tracer_conc_solid_passive_col(c, j, ll) = yf(ll, c, j)
                  endif
               enddo
            enddo
         enddo
      enddo

      k1 = betrtracer_vars%id_trc_o2   ; k2 = centurybgc_vars%lid_o2  ;
      call assign_A2B(bounds, lbj, ubj, neq, ngwtracers, numf, filter, jtops, k1, k2, yf(:,bounds%begc:bounds%endc, lbj:ubj), tracer_conc_mobile_col(bounds%begc:bounds%endc,lbj:ubj,:))
      k1 = betrtracer_vars%id_trc_co2x ; k2 = centurybgc_vars%lid_co2 ;
      call assign_A2B(bounds, lbj, ubj, neq, ngwtracers, numf, filter, jtops, k1, k2, yf(:,bounds%begc:bounds%endc, lbj:ubj), tracer_conc_mobile_col(bounds%begc:bounds%endc,lbj:ubj,:))
      k1 = betrtracer_vars%id_trc_nh3x ; k2 = centurybgc_vars%lid_nh4 ;
      call assign_A2B(bounds, lbj, ubj, neq, ngwtracers, numf, filter, jtops, k1, k2, yf(:,bounds%begc:bounds%endc, lbj:ubj), tracer_conc_mobile_col(bounds%begc:bounds%endc,lbj:ubj,:))
      k1 = betrtracer_vars%id_trc_no3x ; k2 = centurybgc_vars%lid_no3 ;
      call assign_A2B(bounds, lbj, ubj, neq, ngwtracers, numf, filter, jtops, k1, k2, yf(:,bounds%begc:bounds%endc, lbj:ubj), tracer_conc_mobile_col(bounds%begc:bounds%endc,lbj:ubj,:))
      k1 = betrtracer_vars%id_trc_n2   ; k2 = centurybgc_vars%lid_n2  ;
      call assign_A2B(bounds, lbj, ubj, neq, ngwtracers, numf, filter, jtops, k1, k2, yf(:,bounds%begc:bounds%endc, lbj:ubj), tracer_conc_mobile_col(bounds%begc:bounds%endc,lbj:ubj,:))
      k1 = betrtracer_vars%id_trc_ar   ; k2 = centurybgc_vars%lid_ar  ;
      call assign_A2B(bounds, lbj, ubj, neq, ngwtracers, numf, filter, jtops, k1, k2, yf(:,bounds%begc:bounds%endc, lbj:ubj), tracer_conc_mobile_col(bounds%begc:bounds%endc,lbj:ubj,:))
      k1 = betrtracer_vars%id_trc_ch4  ; k2 = centurybgc_vars%lid_ch4 ;
      call assign_A2B(bounds, lbj, ubj, neq, ngwtracers, numf, filter, jtops, k1, k2, yf(:,bounds%begc:bounds%endc, lbj:ubj), tracer_conc_mobile_col(bounds%begc:bounds%endc,lbj:ubj,:))
      k1 =  betrtracer_vars%id_trc_n2o ; k2 = centurybgc_vars%lid_n2o ;
      call assign_A2B(bounds, lbj, ubj, neq, ngwtracers, numf, filter, jtops, k1, k2, yf(:,bounds%begc:bounds%endc, lbj:ubj), tracer_conc_mobile_col(bounds%begc:bounds%endc,lbj:ubj,:))

    end associate
  end subroutine retrieve_state_vars
  !-------------------------------------------------------------------------------

  subroutine assign_A2B(bounds, lbj, ubj, neq, ngwtracers,  numf, filter, jtops, &
       k1, k2, yf, tracer_conc_mobile_col)

    !
    ! !DESCRIPTION:
    ! assign state variables
    !
    implicit none
    ! !ARGUMENTS:
    type(bounds_type),      intent(in)    :: bounds
    integer,                intent(in)    :: lbj, ubj
    integer,                intent(in)    :: jtops(bounds%begc:bounds%endc)        ! top label of each column
    integer,                intent(in)    :: numf
    integer,                intent(in)    :: filter(:)
    integer,                intent(in)    :: k1, k2, neq, ngwtracers
    real(r8),               intent(in)    :: yf(neq, bounds%begc:bounds%endc, lbj:ubj)
    real(r8),               intent(inout) :: tracer_conc_mobile_col(bounds%begc:bounds%endc,lbj:ubj,ngwtracers)

    ! !LOCAL VARIABLES:
    integer :: j, fc, c

    do j = lbj, ubj   !currently, om is added only to soil layers
       do fc = 1, numf
          c = filter(fc)
          if(j>=jtops(c))then
             tracer_conc_mobile_col(c, j, k1) = yf(k2, c, j)
          endif
       enddo
    enddo
  end subroutine assign_A2B

  !-------------------------------------------------------------------------------
  subroutine calc_nitrif_denitrif_rate(bounds, lbj, ubj, numf, filter, jtops, dz, t_soisno, pH, &
       pot_co2_hr, anaerobic_frac, smin_nh4_vr, smin_no3_vr, soilstate_vars, waterstate_vars,   &
       carbonflux_vars, n2_n2o_ratio_denit, nh4_no3_ratio, decay_nh4, decay_no3)
    !
    ! !DESCRIPTION:
    ! calculate nitrification denitrification rate
    ! the actual nitrification rate will be f_nitr * [nh4]
    ! and the actual denitri rate will be of f_denit * [no3]
    !
    ! !USES:
    use clm_varcon          , only : rpi, secspday
    use SoilStatetype       , only : soilstate_type
    use WaterStateType      , only : waterstate_type
    use MathfuncMod         , only : safe_div
    use shr_const_mod       , only : SHR_CONST_TKFRZ
    use CNCarbonFluxType    , only : carbonflux_type

    ! !ARGUMENTS:
    type(bounds_type),      intent(in)  :: bounds
    integer,                intent(in)  :: lbj, ubj
    integer,                intent(in)  :: jtops(bounds%begc: )                                            ! top label of each column
    integer,                intent(in)  :: numf
    integer,                intent(in)  :: filter(:)
    real(r8),               intent(in)  :: dz(bounds%begc: , lbj: )
    real(r8),               intent(in)  :: pH(bounds%begc: , lbj: )                                        !pH of soil
    real(r8),               intent(in)  :: t_soisno(bounds%begc: , lbj: )                                  !soil temperature
    real(r8),               intent(in)  :: pot_co2_hr(bounds%begc: , lbj: )                                !potential aerobic heteotrophic respiration, mol CO2/m3/s
    real(r8),               intent(in)  :: anaerobic_frac(bounds%begc: , lbj: )                            !fraction of anaerobic soil
    real(r8),               intent(in)  :: smin_nh4_vr(bounds%begc: , lbj: )
    real(r8),               intent(in)  :: smin_no3_vr(bounds%begc: , lbj: )                               !soil no3 concentration [mol N/m3]
    type(waterstate_type) , intent(in)  :: waterstate_vars
    type(soilstate_type)  , intent(in)  :: soilstate_vars
    type(carbonflux_type) , intent(in)  :: carbonflux_vars
    real(r8)             ,  intent(out) :: n2_n2o_ratio_denit(bounds%begc: , lbj: )                        !ratio of n2 to n2o in denitrification
    real(r8)             ,  intent(out) :: nh4_no3_ratio(bounds%begc: , lbj: )                             !ratio of soil nh4 and no3
    real(r8)             ,  intent(out) :: decay_nh4(bounds%begc: ,lbj: )                                  !1/s, decay rate of nh4
    real(r8)             ,  intent(out) :: decay_no3(bounds%begc: ,lbj: )                                  !1/s, decay rate of no3

                                                                                                           ! !LOCAL VARIABLES:
    logical, parameter                  :: no_frozen_nitrif_denitrif = .false.                             !this is a testing parameter, just to make the model run
    real(r8)                            :: soil_hr_vr(bounds%begc:bounds%endc,   lbj:ubj)
    real(r8)                            :: k_nitr_t_vr(bounds%begc:bounds%endc,  lbj:ubj)
    real(r8)                            :: k_nitr_ph_vr(bounds%begc:bounds%endc, lbj:ubj)
    real(r8)                            :: k_nitr_h2o_vr(bounds%begc:bounds%endc,lbj:ubj)
    real(r8)                            :: k_nitr_vr(bounds%begc:bounds%endc,    lbj:ubj)
    real(r8)                            :: pot_f_nit_vr(bounds%begc:bounds%endc, lbj:ubj)
    real(r8)                            :: soil_bulkdensity(bounds%begc:bounds%endc, lbj:ubj)
    real(r8)                            :: smin_no3_massdens_vr(bounds%begc:bounds%endc,lbj:ubj)
    real(r8)                            :: soil_co2_prod(bounds%begc:bounds%endc, lbj:ubj)
    real(r8)                            :: fmax_denit_carbonsubstrate_vr(bounds%begc:bounds%endc, lbj:ubj) !
    real(r8)                            :: fmax_denit_nitrate_vr(bounds%begc:bounds%endc, lbj:ubj)
    real(r8)                            :: f_denit_base_vr(bounds%begc:bounds%endc, lbj:ubj)
    real(r8)                            :: ratio_k1(bounds%begc:bounds%endc, lbj:ubj)
    real(r8)                            :: diffus(bounds%begc:bounds%endc, lbj:ubj)
    real(r8)                            :: ratio_no3_co2(bounds%begc:bounds%endc, lbj:ubj)
    real(r8)                            :: wfps_vr(bounds%begc:bounds%endc, lbj:ubj)
    real(r8)                            :: fr_WFPS(bounds%begc:bounds%endc, lbj:ubj)
    real(r8)                            :: pot_f_denit_vr(bounds%begc:bounds%endc, lbj:ubj)
    real(r8)                            :: co2diff_con(2) = (/0.1325_r8, 0.0009_r8/)
    real(r8)                            :: g_per_m3__to__ug_per_gsoil
    real(r8)                            :: g_per_m3_sec__to__ug_per_gsoil_day
    real(r8)                            :: k_nitr_max
    integer                             :: fc, c, j

    SHR_ASSERT_ALL((ubound(jtops)              == (/bounds%endc/))      , errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(pH)                 == (/bounds%endc, ubj/)) , errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(t_soisno)           == (/bounds%endc, ubj/)) , errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(pot_co2_hr)         == (/bounds%endc, ubj/)) , errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(anaerobic_frac)     == (/bounds%endc, ubj/)) , errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(smin_nh4_vr)        == (/bounds%endc, ubj/)) , errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(smin_no3_vr)        == (/bounds%endc, ubj/)) , errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(n2_n2o_ratio_denit) == (/bounds%endc, ubj/)) , errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(nh4_no3_ratio)      == (/bounds%endc, ubj/)) , errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(decay_nh4)          == (/bounds%endc, ubj/)) , errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(decay_no3)          == (/bounds%endc, ubj/)) , errMsg(__FILE__,__LINE__))

    associate(                                                                &
         bd                            =>    soilstate_vars%bd_col          , & !
         watsat                        =>    soilstate_vars%watsat_col      , & !
         h2osoi_vol                    =>    waterstate_vars%h2osoi_vol_col , & !
         h2osoi_liq                    =>    waterstate_vars%h2osoi_liq_col , & !
         finundated                    =>    waterstate_vars%finundated_col , & ! Input: [real(r8) (:)]
         t_scalar                      =>    carbonflux_vars%t_scalar_col   , & ! Input: [real(r8) (:,:)   ]  soil temperature scalar for decomp
         w_scalar                      =>    carbonflux_vars%w_scalar_col     & ! Input: [real(r8) (:,:)   ]  soil water scalar for decomp
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
            decay_nh4(c,j) = pot_f_nit_vr(c,j)
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

            decay_no3(c,j) = safe_div(pot_f_denit_vr(c,j), smin_no3_vr(c,j)) ![1/s]
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
    ! !DESCRIPTION:
    !
    ! calculate soil anoxia state for doing nitrification and denitrification
    ! Rewritten based on Charlie Koven's code by Jinyun Tang
    ! !USES:
    use CNSharedParamsMod   , only : CNParamsShareInst
    use clm_varcon          , only : d_con_g, grav, d_con_w
    use SoilStatetype       , only : soilstate_type
    use BGCCenturyParMod    , only : CNNitrifDenitrifParamsInst

    ! !ARGUMENTS:
    type(bounds_type),      intent(in) :: bounds
    integer,                intent(in) :: lbj, ubj
    integer,                intent(in) :: jtops(bounds%begc: )                       ! top label of each column
    integer,                intent(in) :: numf
    integer,                intent(in) :: filter(:)                                  !indices
    type(soilstate_type),   intent(in) :: soilstate_vars
    real(r8),               intent(in) :: t_soisno(bounds%begc: , lbj: )
    real(r8),               intent(in) :: h2osoi_vol(bounds%begc:, lbj: )
    real(r8),               intent(in) :: o2_decomp_depth_unsat(bounds%begc: ,lbj: ) !potential o2 consumption, as deduced from aerobic heteorotrophic decomposition, mol o2/m3/s
    real(r8),               intent(in) :: conc_o2_unsat(bounds%begc: , lbj: )        !bulk soil o2 concentration, mol/m3
    real(r8),              intent(out) :: anaerobic_frac(bounds%begc: , lbj: )       !fraction of aerobic soil
    ! !LOCAL VARIABLES:
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

    SHR_ASSERT_ALL((ubound(jtops)                 == (/bounds%endc/))      , errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(t_soisno)              == (/bounds%endc, ubj/)) , errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(o2_decomp_depth_unsat) == (/bounds%endc, ubj/)) , errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(conc_o2_unsat)         == (/bounds%endc, ubj/)) , errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(anaerobic_frac)        == (/bounds%endc, ubj/)) , errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(h2osoi_vol)            == (/bounds%endc, ubj/)) , errMsg(__FILE__,__LINE__))

    associate(                                                            &
         watsat                        =>    soilstate_vars%watsat_col  , &
         watfc                         =>    soilstate_vars%watfc_col   , &
         bsw                           =>    soilstate_vars%bsw_col     , &
         sucsat                        =>    soilstate_vars%sucsat_col  , &
         soilpsi                       =>    soilstate_vars%soilpsi_col , &
         cellorg                       =>    soilstate_vars%cellorg_col   &
         )

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

  subroutine calc_potential_aerobic_hr(bounds, lbj, ubj, numf, filter, jtops, cn_ratios, cp_ratios, &
       centurybgc_vars, pot_decay_rates, pct_sand, pot_co2_hr, pot_nh3_immob)
    !
    ! DESCRIPTION:
    ! calculate potential aerobic heteorotrophic respiration, and potential oxygen consumption based on cascade_matrix
    ! !USES:
    use MathfuncMod         , only : dot_sum
    use CNSharedParamsMod   , only : CNParamsShareInst
    use CNSharedParamsMod   , only : CNParamsShareInst
    use BGCCenturyParMod    , only : CNDecompBgcParamsInst
    use MathfuncMod         , only : safe_div

    ! !ARGUMENTS:
    type(bounds_type)       , intent(in) :: bounds
    integer                 , intent(in) :: lbj, ubj
    integer                 , intent(in) :: jtops(bounds%begc:bounds%endc)        ! top label of each column
    integer                 , intent(in) :: numf
    integer                 , intent(in) :: filter(:)
    type(centurybgc_type)   , intent(in) :: centurybgc_vars
    real(r8)                , intent(in) :: cn_ratios(centurybgc_vars%nom_pools, bounds%begc:bounds%endc, lbj:ubj)
    real(r8)                , intent(in) :: cp_ratios(centurybgc_vars%nom_pools, bounds%begc:bounds%endc, lbj:ubj)
    real(r8)                , intent(in) :: pot_decay_rates(centurybgc_vars%nom_pools, bounds%begc:bounds%endc, lbj:ubj)
    real(r8)                , intent(in) :: pct_sand(bounds%begc:bounds%endc, lbj:ubj)
    real(r8)                , intent(out):: pot_co2_hr(bounds%begc:bounds%endc, lbj:ubj)
    real(r8)                , intent(out):: pot_nh3_immob(bounds%begc:bounds%endc, lbj:ubj)

    ! !LOCAL VARIABLES:
    real(r8) :: ftxt, f1, f2
    real(r8) :: cascade_matrix_hr(centurybgc_vars%nom_pools)
    real(r8) :: cascade_matrix_nh3(centurybgc_vars%nom_pools)

    integer  :: fc, c, j, reac

    associate(                                           & !
         nom_pools => centurybgc_vars%nom_pools        , & !
         lit1      => centurybgc_vars%lit1             , & !
         lit2      => centurybgc_vars%lit2             , & !
         lit3      => centurybgc_vars%lit3             , & !
         som1      => centurybgc_vars%som1             , & !
         som2      => centurybgc_vars%som2             , & !
         som3      => centurybgc_vars%som3             , & !
         cwd       => centurybgc_vars%cwd              , & !
         lit1_dek_reac=> centurybgc_vars%lit1_dek_reac , & !
         lit2_dek_reac=> centurybgc_vars%lit2_dek_reac , & !
         lit3_dek_reac=> centurybgc_vars%lit3_dek_reac , & !
         som1_dek_reac=> centurybgc_vars%som1_dek_reac , & !
         som2_dek_reac=> centurybgc_vars%som2_dek_reac , & !
         som3_dek_reac=> centurybgc_vars%som3_dek_reac , & !
         cwd_dek_reac=> centurybgc_vars%cwd_dek_reac     & !
         )

      do j = lbj, ubj
         do fc = 1, numf
            c = filter(fc)
            if(j<jtops(c))cycle
            cascade_matrix_hr = 0._r8
            !reaction1, lit1 -> som1
            !lit1 + 0.55*o2 -> 0.45 som1 + 0.55co2 + (1/cn_ratios(lit1) - 0.45/cn_ratios(som1))+ (1/cp_ratios(lit1)-0.45/cp_ratios(som1))
            reac=lit1_dek_reac
            cascade_matrix_hr(reac)  =  CNDecompBgcParamsInst%rf_l1s1_bgc
            cascade_matrix_nh3(reac) =  safe_div(1._r8,cn_ratios(lit1,c,j)) - safe_div(1._r8-CNDecompBgcParamsInst%rf_l1s1_bgc,cn_ratios(som1,c,j))

            cascade_matrix_nh3(reac) = min(cascade_matrix_nh3(reac) , 0._r8)
            !reaction 2, lit2 -> som1
            !lit2 + 0.5 o2  -> 0.5 som1 + 0.5 co2 + (1/cn_ratios(lit2)-0.5/cn_ratios(som1)) +(1/cp_ratios(lit2)-0.5/cp_ratios(som1))
            reac = lit2_dek_reac
            cascade_matrix_hr(reac)  =  CNDecompBgcParamsInst%rf_l2s1_bgc
            cascade_matrix_nh3(reac) = safe_div(1._r8,cn_ratios(lit2,c,j)) - safe_div(1._r8-CNDecompBgcParamsInst%rf_l2s1_bgc,cn_ratios(som1,c,j))

            cascade_matrix_nh3(reac) = min(cascade_matrix_nh3(reac) , 0._r8)

            !reaction 3, lit3 -> som2
            !lit3 + 0.5 o2 -> 0.5 som2 + 0.5 co2
            reac = lit3_dek_reac
            cascade_matrix_hr(reac)  =  CNDecompBgcParamsInst%rf_l3s2_bgc
            cascade_matrix_nh3(reac) = safe_div(1._r8,cn_ratios(lit3,c,j)) - safe_div(1._r8-CNDecompBgcParamsInst%rf_l3s2_bgc,cn_ratios(som2,c,j))

            cascade_matrix_nh3(reac) = min(cascade_matrix_nh3(reac) , 0._r8)

            !reaction 4, the partition into som2 and som3 is soil texture dependent, som1->som2, som3
            !som1 + f(txt) o2 -> f1*som2 + f2*som3 + f(txt) co2 + (1/cn_ratios(som1)-f1/cn_ratios(som2)-f2/cn_ratios(som3)) +(1/cp_ratios(som1)-f1/cp_ratios(som2)-f2/cp_ratios(som3))
            reac = som1_dek_reac
            !f(txt) = 0.85_r8 - 0.68_r8 * 0.01_r8 * (100._r8 - sand), assuming sand = 30%
            !f1+f2+f(txt)=1._r8
            ftxt = 0.85_r8 - 0.68_r8 * 0.01_r8 * (100._r8 - pct_sand(c,j))
            f1 = 0.996*(1._r8-ftxt)
            f2 = 0.004*(1._r8-ftxt)

            cascade_matrix_hr(reac)  = ftxt
            cascade_matrix_nh3(reac) = safe_div(1._r8,cn_ratios(som1,c,j))-safe_div(f1,cn_ratios(som2,c,j))-safe_div(f2,cn_ratios(som3,c,j))

            cascade_matrix_nh3(reac) = min(cascade_matrix_nh3(reac) , 0._r8)

            !reaction 5, som2 -> som1, som3
            !som2 + 0.55 o2 -> 0.42 som1 + 0.03som3 + 0.55co2 + (1/cn_ratios(som2)-0.42/cn_ratios(som1)-0.03/cn_ratios(som3)) + (1/cp_raitos(som2)-0.42/cp_ratios(som1)-0.03/cp_ratios(som3))
            reac = som2_dek_reac
            cascade_matrix_hr(reac)   =  CNDecompBgcParamsInst%rf_s2s1_bgc
            cascade_matrix_nh3(reac)  = safe_div(1._r8,cn_ratios(som2,c,j))-0.93_r8*safe_div(1._r8-CNDecompBgcParamsInst%rf_s2s1_bgc,cn_ratios(som1,c,j)) &
                 -0.07_r8*safe_div(1._r8-CNDecompBgcParamsInst%rf_s2s1_bgc,cn_ratios(som3,c,j))
            cascade_matrix_nh3(reac) = min(cascade_matrix_nh3(reac) , 0._r8)

            !reaction 6, s3 -> s1
            !som3 + 0.55 o2 -> 0.45*som1 + 0.55co2 + (1/cn_ratios(som3)-0.45/cn_ratios(som1)) + (1/cp_ratios(som3)-0.45/cp_ratios(som1))
            reac = som3_dek_reac
            cascade_matrix_hr(reac) = CNDecompBgcParamsInst%rf_s3s1_bgc
            cascade_matrix_nh3(reac)= safe_div(1._r8,cn_ratios(som3,c,j)) - safe_div(1._r8-CNDecompBgcParamsInst%rf_s3s1_bgc,cn_ratios(som1,c,j))

            cascade_matrix_nh3(reac) = min(cascade_matrix_nh3(reac) , 0._r8)

            !cwd -> lit1, lit2
            reac = cwd_dek_reac
            cascade_matrix_nh3(reac) = safe_div(1._r8,cn_ratios(cwd,c,j)) - safe_div(CNDecompBgcParamsInst%cwd_fcel_bgc,cn_ratios(lit2,c,j)) - &
                 safe_div(CNDecompBgcParamsInst%cwd_flig_bgc,cn_ratios(lit3,c,j))

            cascade_matrix_nh3(reac) = min(cascade_matrix_nh3(reac) , 0._r8)

            !obtain the potential respiration
            pot_co2_hr(c,j) = dot_sum(cascade_matrix_hr, pot_decay_rates(:, c, j))  !mol CO2/m3/s
            pot_nh3_immob(c,j) = dot_sum(cascade_matrix_nh3,pot_decay_rates(:,c,j)) !mol NH3/m3/s, this does not include mineralization
         enddo
      enddo
    end associate
  end subroutine calc_potential_aerobic_hr

  !----------------------------------------------------------------------------------------------------
  subroutine calc_decompK_multiply_scalar(bounds, lbj, ubj, numf, filter, jtops, finundated, zsoi, &
       t_soisno, o2_bulk, o2_aqu2bulkcef, soilstate_vars, centurybgc_vars, carbonflux_vars)
    !
    ! !DESCRIPTION:
    ! compute scalar multipliers for aerobic om decomposition
    ! because temperature and moisture scalars will not be independent from each other for microbe explicit models,
    ! temp_scalar and moist_scalar are used as private variables

    ! !USES:
    !
    use CNSharedParamsMod   , only : CNParamsShareInst
    use shr_const_mod       , only : SHR_CONST_TKFRZ, SHR_CONST_PI
    use SoilStatetype       , only : soilstate_type
    use CNCarbonFluxType    , only : carbonflux_type
    !
    ! !ARGUMENTS:
    type(bounds_type),         intent(in)  :: bounds
    integer,                   intent(in)  :: lbj, ubj
    integer,                   intent(in)  :: jtops(bounds%begc:bounds%endc)        ! top label of each column
    integer,                   intent(in)  :: numf
    integer,                   intent(in)  :: filter(:)
    real(r8),                  intent(in)  :: t_soisno(bounds%begc:bounds%endc, lbj:ubj)
    real(r8),                  intent(in)  :: o2_bulk(bounds%begc:bounds%endc,lbj:ubj)
    real(r8),                  intent(in)  :: zsoi(bounds%begc:bounds%endc, lbj:ubj)
    real(r8),                  intent(in)  :: finundated(bounds%begc:bounds%endc)
    real(r8),                  intent(in)  :: o2_aqu2bulkcef(bounds%begc:bounds%endc, lbj:ubj)
    type(soilstate_type),      intent(in)  :: soilstate_vars
    type(centurybgc_type)  , intent(inout) :: centurybgc_vars
    type(carbonflux_type)  , intent(inout) :: carbonflux_vars

    ! !LOCAL VARIABLES:
    integer                                :: fc, c, j                    !indices
    real(r8), parameter                    :: normalization_tref = 15._r8 ! reference temperature for normalizaion (degrees C)
    real(r8)                               :: decomp_depth_efolding       ! [m] a testing parameter, which will be replaced,
    real(r8)                               :: Q10                         ! a number taken from CLM4.5bgc
    real(r8)                               :: froz_q10
    real(r8)                               :: minpsi
    real(r8)                               :: maxpsi
    real(r8)                               :: normalization_factor
    real(r8)                               :: catanf_30
    real(r8)                               :: catanf
    real(r8)                               :: t1
    real(r8)                               :: o2w
    real(r8)                               :: psi

    !----- CENTURY T response function
    catanf(t1) = 11.75_r8 +(29.7_r8 / SHR_CONST_PI) * atan( SHR_CONST_PI * 0.031_r8  * ( t1 - 15.4_r8 ))

    associate(                                                 &
         sucsat         => soilstate_vars%sucsat_col         , & ! Input:  [real(r8) (:,:)] minimum soil suction [mm]
         soilpsi        => soilstate_vars%soilpsi_col        , & ! Input:  [real(r8) (:,:)] soilwater pontential in each soil layer [MPa]
         t_scalar       => carbonflux_vars%t_scalar_col      , & ! Output: [real(r8) (:,:)   ]  soil temperature scalar for decomp
         w_scalar       => carbonflux_vars%w_scalar_col      , & ! Output: [real(r8) (:,:)   ]  soil water scalar for decomp
         o_scalar       => carbonflux_vars%o_scalar_col      , & ! Output: [real(r8) (:,:)   ]  fraction by which decomposition is limited by anoxia
         depth_scalar   => centurybgc_vars%depth_scalar_col    & ! Output: [real(r8) (:,:,:) ]  rate constant for decomposition (1./sec)
         )

      catanf_30 = catanf(30._r8)
      
      ! set "Q10" parameter
      Q10 = CNParamsShareInst%Q10

      ! set "froz_q10" parameter
      froz_q10  = CNParamsShareInst%froz_q10

      ! Set "decomp_depth_efolding" parameter
      decomp_depth_efolding = CNParamsShareInst%decomp_depth_efolding

      do j = lbj, ubj
         do fc = 1, numf
            c = filter(fc)
            if(j<jtops(c))cycle
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
            o_scalar(c,j)     = o2w/(o2w+0.02_r8)   !the value 0.22 mol O3/m3 is from Arah and Kirk, 2000

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
    ! !DESCRIPTION:
    ! calculate the nitrogen uptake profile
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)   :: bounds
    integer                  , intent(in)   :: nlevdecomp                         ! number of vertical layers
    integer                  , intent(in)   :: num_soilc                          ! number of soil columns in filter
    integer                  , intent(in)   :: filter_soilc(:)                    ! filter for soil columns
    real(r8)                 , intent(in)   :: sminn_nh4_vr(bounds%begc: , 1: )   ! soil mineral nitrogen profile
    real(r8)                 , intent(in)   :: sminn_no3_vr(bounds%begc: , 1: )   ! soil mineral nitrogen profile
    real(r8)                 , intent(in)   :: dzsoi(bounds%begc: , 1: )          ! layer thickness
    real(r8)                 , intent(in)   :: nfixation_prof(bounds%begc: , 1: ) ! nitrogen fixation profile
    real(r8)                 , intent(inout):: nuptake_prof(bounds%begc: , 1: )   ! nitrogen uptake profile

    ! !LOCAL VARIABLES:
    integer :: fc, j, c      ! indices
    real(r8):: sminn_tot(bounds%begc:bounds%endc)  !vertically integrated mineral nitrogen


    SHR_ASSERT_ALL((ubound(dzsoi)          == (/bounds%endc, nlevdecomp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sminn_nh4_vr)   == (/bounds%endc, nlevdecomp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(sminn_no3_vr)   == (/bounds%endc, nlevdecomp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(nfixation_prof) == (/bounds%endc, nlevdecomp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dzsoi)          == (/bounds%endc, nlevdecomp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(nuptake_prof)   == (/bounds%endc, nlevdecomp/)), errMsg(__FILE__, __LINE__))

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
    !
    ! !DESCRIPTION:
    !
    !caluate depth specific demand
    ! !USES:
    use clm_varcon    , only : natomw
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)    :: bounds
    integer           , intent(in)    :: nlevdecomp                                               ! number of vertical layers
    integer           , intent(in)    :: num_soilc                                                ! number of soil columns in filter
    integer           , intent(in)    :: filter_soilc(:)                                          ! filter for soil columns
    real(r8)          , intent(in)    :: dzsoi(bounds%begc:bounds%endc,1:nlevdecomp)              ! layer thickness
    real(r8)          , intent(in)    :: plant_totn_demand_flx_col(bounds%begc:bounds%endc)
    real(r8)          , intent(in)    :: nuptake_prof(bounds%begc:bounds%endc, 1:nlevdecomp)
    real(r8)          , intent(inout) :: plant_demand_vr(1,bounds%begc:bounds%endc, 1:nlevdecomp) !mol N/m3/s
    ! !LOCAL VARIABLES:
    integer :: fc, c, j


    do j = 1, nlevdecomp
       do fc = 1, num_soilc
          c = filter_soilc(fc)
          plant_demand_vr(1,c,j) = plant_totn_demand_flx_col(c) * nuptake_prof(c,j) / dzsoi(c,j) /natomw
       enddo
    enddo


  end subroutine calc_plant_nitrogen_uptake_prof


  !-----------------------------------------------------------------------

  subroutine bgcstate_ext_update_bfdecomp(bounds, lbj, ubj, num_soilc, filter_soilc, &
       carbonflux_vars, nitrogenflux_vars, centurybgc_vars, betrtracer_vars, &
       tracerflux_vars, y0, cn_ratios, cp_ratios)
    ! !DESCRIPTION:
    ! update om pools with external input before doing decomposition
    !
    ! !USES:
    use MathfuncMod              , only :  safe_div
    use CNCarbonFluxType         , only : carbonflux_type
    use CNNitrogenFluxType       , only : nitrogenflux_type
    use BetrTracerType           , only : betrtracer_type
    use tracerstatetype          , only : tracerstate_type
    use tracerfluxType           , only : tracerflux_type
    use CNDecompCascadeConType , only : decomp_cascade_con
    !
    ! !ARGUMENTS:
    type(bounds_type)       , intent(in) :: bounds                                                                    ! bounds
    integer                 , intent(in) :: num_soilc                                                                 ! number of columns in column filter
    integer                 , intent(in) :: filter_soilc(:)                                                           ! column filter
    integer                 , intent(in) :: lbj, ubj
    type(carbonflux_type)   , intent(in) :: carbonflux_vars
    type(nitrogenflux_type) , intent(in) :: nitrogenflux_vars
    type(betrtracer_type)   , intent(in) :: betrtracer_vars                                                           ! betr configuration information
    type(centurybgc_type)   , intent(in) :: centurybgc_vars
    type(tracerflux_type)   , intent(inout) :: tracerflux_vars
    real(r8)                , intent(inout) :: y0(centurybgc_vars%nstvars, bounds%begc:bounds%endc, lbj:ubj)
    real(r8)                , intent(inout) :: cn_ratios(centurybgc_vars%nom_pools, bounds%begc:bounds%endc, lbj:ubj) !
    real(r8)                , intent(inout) :: cp_ratios(centurybgc_vars%nom_pools, bounds%begc:bounds%endc, lbj:ubj)
    ! !LOCAL VARIABLES:
    real(r8) :: delta_no3, delta_nh4
    real(r8) :: delta_somn
    integer :: k, fc, c, j

    associate(                                                                                 & !
         lid_nh4                        => centurybgc_vars%lid_nh4                           , & !
         lid_no3                        => centurybgc_vars%lid_no3                           , & !
         nelm                           => centurybgc_vars%nelms                             , & !
         c_loc                          => centurybgc_vars%c_loc                             , & !
         n_loc                          => centurybgc_vars%n_loc                             , & !
         id_trc_nh3x                    => betrtracer_vars%id_trc_nh3x                       , & !
         id_trc_no3x                    => betrtracer_vars%id_trc_no3x                       , & !
         ngwmobile_tracers              => betrtracer_vars%ngwmobile_tracers                 , & !
         tracer_flx_netpro_vr           => tracerflux_vars%tracer_flx_netpro_vr_col          , & !
         bgc_cpool_ext_loss_vr          => carbonflux_vars%bgc_cpool_ext_loss_vr_col         , & !
         bgc_npool_ext_loss_vr          => nitrogenflux_vars%bgc_npool_ext_loss_vr_col       , & !
         sminn_nh4_input_vr             => nitrogenflux_vars%sminn_nh4_input_vr_col          , & !
         sminn_no3_input_vr             => nitrogenflux_vars%sminn_no3_input_vr_col          , & !
         initial_cn_ratio               => decomp_cascade_con%initial_cn_ratio               , & ! Output: [real(r8)          (:)     ]  c:n ratio for initialization of pools
         floating_cn_ratio_decomp_pools => decomp_cascade_con%floating_cn_ratio_decomp_pools   & ! Output: [logical           (:)     ]  TRUE => pool has fixed C:N ratio         
         )
    
      do k = 1, ndecomp_pools
         do j = 1, ubj
            do fc = 1, num_soilc
               c = filter_soilc(fc)
               y0((k-1)*nelm+c_loc,c,j) = y0((k-1)*nelm+c_loc,c,j) - bgc_cpool_ext_loss_vr(c,j,k)/catomw
               y0((k-1)*nelm+n_loc,c,j) = y0((k-1)*nelm+n_loc,c,j) - bgc_npool_ext_loss_vr(c,j,k)/natomw
               if(floating_cn_ratio_decomp_pools(k))then
                  cn_ratios(k, c,j) = safe_div(y0((k-1)*nelm+c_loc,c,j), y0((k-1)*nelm+n_loc,c,j))
               else
                  cn_ratios(k,c,j) = initial_cn_ratio(k)*natomw/catomw

               endif

               tracer_flx_netpro_vr(c,j,ngwmobile_tracers+(k-1)*nelm+c_loc) =      &
                    tracer_flx_netpro_vr(c,j,ngwmobile_tracers+(k-1)*nelm+c_loc) - &
                    bgc_cpool_ext_loss_vr(c,j,k)/catomw

               tracer_flx_netpro_vr(c,j,ngwmobile_tracers+(k-1)*nelm+n_loc) =      &
                    tracer_flx_netpro_vr(c,j,ngwmobile_tracers+(k-1)*nelm+n_loc) - &
                    bgc_npool_ext_loss_vr(c,j,k)/natomw

            enddo
         enddo
      enddo

      do j = 1, ubj
         do fc = 1, num_soilc
            c = filter_soilc(fc)
            y0(lid_nh4, c, j) = y0(lid_nh4, c, j) + sminn_nh4_input_vr(c,j)/natomw
            y0(lid_no3, c, j) = y0(lid_no3, c, j) + sminn_no3_input_vr(c,j)/natomw

            tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_no3x   ) =      &
                 tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_no3x   ) + &
                 sminn_no3_input_vr(c,j)/natomw

            tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_nh3x   ) =      &
                 tracer_flx_netpro_vr(c,j,betrtracer_vars%id_trc_nh3x   ) + &
                 sminn_nh4_input_vr(c,j)/natomw
         enddo

      enddo

    end associate
  end subroutine bgcstate_ext_update_bfdecomp
  !-----------------------------------------------------------------------

  subroutine bgcstate_ext_update_afdecomp(bounds, lbj, ubj, num_soilc, filter_soilc, &
       carbonflux_vars, nitrogenflux_vars,  centurybgc_vars, betrtracer_vars, tracerflux_vars, yf)

    ! !DESCRIPTION:
    ! update om state variables after doing decomposition.
    ! !USES:
    use MathfuncMod         , only :  safe_div
    use CNCarbonFluxType    , only : carbonflux_type
    use CNNitrogenFluxType  , only : nitrogenflux_type
    use BetrTracerType      , only : betrtracer_type
    use tracerstatetype     , only : tracerstate_type
    use tracerfluxType      , only : tracerflux_type
    !
    ! !ARGUMENTS:
    type(bounds_type)       , intent(in)    :: bounds                                                        ! bounds
    integer                 , intent(in)    :: num_soilc                                                     ! number of columns in column filter
    integer                 , intent(in)    :: filter_soilc(:)                                               ! column filter
    integer                 , intent(in)    :: lbj, ubj
    type(carbonflux_type)   , intent(in)    :: carbonflux_vars
    type(nitrogenflux_type) , intent(in)    :: nitrogenflux_vars
    type(betrtracer_type)   , intent(in)    :: betrtracer_vars                                               ! betr configuration information
    type(centurybgc_type)   , intent(in)    :: centurybgc_vars
    type(tracerflux_type)   , intent(inout) :: tracerflux_vars
    real(r8)                , intent(inout) :: yf(centurybgc_vars%nstvars, bounds%begc:bounds%endc, lbj:ubj) !

    ! !LOCAL VARIABLES:
    real(r8) :: delta_no3, delta_nh4
    real(r8) :: delta_somn
    integer :: k, fc, c, j

    associate(                                                                    & !
         nelm                    => centurybgc_vars%nelms                       , & !
         c_loc                   => centurybgc_vars%c_loc                       , & !
         n_loc                   => centurybgc_vars%n_loc                       , & !
         tracer_flx_netpro_vr    => tracerflux_vars%tracer_flx_netpro_vr_col    , & !
         bgc_cpool_ext_inputs_vr => carbonflux_vars%bgc_cpool_ext_inputs_vr_col , & !
         bgc_npool_ext_inputs_vr => nitrogenflux_vars%bgc_npool_ext_inputs_vr_col & !
         )    

      do k = 1, ndecomp_pools
         do j = 1, ubj
            do fc = 1, num_soilc
               c = filter_soilc(fc)

               yf((k-1)*nelm+c_loc,c,j) = yf((k-1)*nelm+c_loc,c,j) + bgc_cpool_ext_inputs_vr(c,j,k)/catomw
               yf((k-1)*nelm+n_loc,c,j) = yf((k-1)*nelm+n_loc,c,j) + bgc_npool_ext_inputs_vr(c,j,k)/natomw

            enddo
         enddo
      enddo

    end associate
  end subroutine bgcstate_ext_update_afdecomp
  !-----------------------------------------------------------------------
  subroutine apply_plant_root_respiration_prof(bounds, ubj, num_soilc, filter_soilc, &
       rr_col, root_prof_col, rr_col_vr)
    !
    ! !DESCRIPTION:
    ! obtain root respiration profile
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in) :: bounds                                        ! bounds
    integer           , intent(in) :: ubj
    integer           , intent(in) :: num_soilc                                     ! number of columns in column filter
    integer           , intent(in) :: filter_soilc(:)                               ! column filter
    real(r8)          , intent(in) :: rr_col(bounds%begc:bounds%endc)
    real(r8)          , intent(in) :: root_prof_col(bounds%begc:bounds%endc, 1:ubj)
    real(r8)          , intent(inout):: rr_col_vr(1,bounds%begc:bounds%endc, 1:ubj) !
    ! !LOCAL VARIABLES:
    integer :: fc, c, j

    do j = 1, ubj
       do fc = 1, num_soilc
          c = filter_soilc(fc)
          rr_col_vr(1,c,j) = rr_col(c) * root_prof_col(c,j)
       enddo
    enddo

  end subroutine apply_plant_root_respiration_prof

  !-----------------------------------------------------------------------
  subroutine set_reaction_order( nreact, centurybgc_vars, is_zero_order)
    !
    ! !DESCRIPTION:
    ! set order of the reactions, 0 or 1
    !
    ! !ARGUMENTS:
    integer                      , intent(in)  :: nreact
    type(centurybgc_type)        , intent(in)  :: centurybgc_vars
    logical                      , intent(out) :: is_zero_order(nreact)



    is_zero_order(:) = .false.
    is_zero_order(centurybgc_vars%lid_o2_aere_reac)  = .true.
    if(spinup_state /= 1)then
       is_zero_order(centurybgc_vars%lid_n2o_aere_reac) = .true.
       is_zero_order(centurybgc_vars%lid_ar_aere_reac)  = .true.
       is_zero_order(centurybgc_vars%lid_ch4_aere_reac) = .true.
       is_zero_order(centurybgc_vars%lid_o2_aere_reac)  = .true.

       is_zero_order(centurybgc_vars%lid_co2_aere_reac) = .true.
       is_zero_order(centurybgc_vars%lid_n2_aere_reac)  = .true.
    endif

    is_zero_order(centurybgc_vars%lid_plant_minn_up_reac) = .true.
    is_zero_order(centurybgc_vars%lid_at_rt_reac)         = .true.

  end subroutine set_reaction_order


  !-----------------------------------------------------------------------
  subroutine calc_nutrient_compet_rescal(bounds, ubj, num_soilc, filter_soilc, &
       dtime, centurybgc_vars,  k_nit, decomp_nh4_immob, plant_ndemand, smin_nh4_vr, nh4_compet)

    !
    ! !DESCRIPTION:
    ! scaling factor for nitrogen competition
    ! !USES:
    use MathfuncMod       , only : safe_div
                                                                                      ! !ARGUMENTS:
    type(bounds_type)     , intent(in) :: bounds                                      ! bounds
    integer               , intent(in) :: ubj
    integer               , intent(in) :: num_soilc                                   ! number of columns in column filter
    integer               , intent(in) :: filter_soilc(:)                             ! column filter
    type(centurybgc_type) , intent(in) :: centurybgc_vars
    real(r8)              , intent(in) :: dtime
    real(r8)              , intent(in) :: k_nit(bounds%begc: , 1: )
    real(r8)              , intent(in) :: decomp_nh4_immob(bounds%begc: , 1: )
    real(r8)              , intent(in) :: plant_ndemand(bounds%begc: , 1: )
    real(r8)              , intent(in) :: smin_nh4_vr(bounds%begc: , 1: )
    real(r8)              , intent(inout):: nh4_compet(bounds%begc:bounds%endc,1:ubj) !
    ! !LOCAL VARIABLES:
    integer :: j, fc, c
    real(r8):: tot_demand

    SHR_ASSERT_ALL((ubound(k_nit)            == (/bounds%endc, ubj/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(decomp_nh4_immob) == (/bounds%endc, ubj/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(plant_ndemand)    == (/bounds%endc, ubj/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(smin_nh4_vr)      == (/bounds%endc, ubj/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(nh4_compet)       == (/bounds%endc, ubj/)), errMsg(__FILE__, __LINE__))

    do j = 1, ubj
       do fc = 1, num_soilc
          c = filter_soilc(fc)

          tot_demand= (k_nit(c,j) * smin_nh4_vr(c,j) + decomp_nh4_immob(c,j) + plant_ndemand(c,j))*dtime
          if(tot_demand<=smin_nh4_vr(c,j))then
             nh4_compet(c,j)=1._r8
          else
             nh4_compet(c,j) = smin_nh4_vr(c,j)/tot_demand
          endif
       enddo
    enddo
  end subroutine calc_nutrient_compet_rescal

  !-----------------------------------------------------------------------
  subroutine assign_OM_CNpools(bounds, num_soilc, filter_soilc,  carbonstate_vars, &
       nitrogenstate_vars, tracerstate_vars, betrtracer_vars, centurybgc_vars)

    ! !DESCRIPTION:
    ! update OM pools
    ! !USES:
    use clm_varpar           , only : i_cwd, i_met_lit, i_cel_lit, i_lig_lit
    use CNCarbonStateType    , only : carbonstate_type
    use CNNitrogenStateType  , only : nitrogenstate_type
    use tracerstatetype      , only : tracerstate_type
    use BetrTracerType       , only : betrtracer_type
    use clm_varpar           , only : nlevtrc_soil

    ! !ARGUMENTS:
    type(bounds_type)        , intent(in) :: bounds                ! bounds
    integer                  , intent(in) :: num_soilc             ! number of columns in column filter
    integer                  , intent(in) :: filter_soilc(:)       ! column filter
    type(tracerstate_type)   , intent(in) :: tracerstate_vars
    type(centurybgc_type)    , intent(in) :: centurybgc_vars
    type(betrtracer_type)    , intent(in) :: betrtracer_vars       ! betr configuration information
    type(carbonstate_type)   , intent(inout) :: carbonstate_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars !

    ! !LOCAL VARIABLES:
    integer, parameter :: i_soil1 = 5
    integer, parameter :: i_soil2 = 6
    integer, parameter :: i_soil3 = 7

    integer :: fc, c, j, k

    associate(                                                                         &
         id_trc_no3x               => betrtracer_vars%id_trc_no3x                    , &
         id_trc_nh3x               => betrtracer_vars%id_trc_nh3x                    , &
         decomp_cpools_vr          => carbonstate_vars%decomp_cpools_vr_col          , &
         decomp_npools_vr          => nitrogenstate_vars%decomp_npools_vr_col        , &
         smin_no3_vr_col           => nitrogenstate_vars%smin_no3_vr_col             , &
         smin_nh4_vr_col           => nitrogenstate_vars%smin_nh4_vr_col             , &
         sminn_vr_col              => nitrogenstate_vars%sminn_vr_col                , &
         tracer_conc_mobile        => tracerstate_vars%tracer_conc_mobile_col        , &
         tracer_conc_solid_passive => tracerstate_vars%tracer_conc_solid_passive_col , &
         c_loc                     => centurybgc_vars%c_loc                          , &
         n_loc                     => centurybgc_vars%n_loc                          , &
         lit1                      => centurybgc_vars%lit1                           , &
         lit2                      => centurybgc_vars%lit2                           , &
         lit3                      => centurybgc_vars%lit3                           , &
         som1                      => centurybgc_vars%som1                           , &
         som2                      => centurybgc_vars%som2                           , &
         som3                      => centurybgc_vars%som3                           , &
         cwd                       => centurybgc_vars%cwd                            , &
         nelms                     => centurybgc_vars%nelms                            &
         )

      do j = 1, nlevtrc_soil
         do fc = 1, num_soilc
            c = filter_soilc(fc)

            smin_no3_vr_col(c,j) = tracer_conc_mobile(c,j,id_trc_no3x)*natomw
            smin_nh4_vr_col(c,j) = tracer_conc_mobile(c,j,id_trc_nh3x)*natomw
            sminn_vr_col   (c,j) = smin_no3_vr_col(c,j) + smin_nh4_vr_col(c,j)

            k = lit1; decomp_cpools_vr(c,j,i_met_lit) = tracer_conc_solid_passive(c,j,(k-1)*nelms+c_loc) * catomw
            k = lit2; decomp_cpools_vr(c,j,i_cel_lit) = tracer_conc_solid_passive(c,j,(k-1)*nelms+c_loc) * catomw
            k = lit3; decomp_cpools_vr(c,j,i_lig_lit) = tracer_conc_solid_passive(c,j,(k-1)*nelms+c_loc) * catomw
            k = cwd ; decomp_cpools_vr(c,j,i_cwd    ) = tracer_conc_solid_passive(c,j,(k-1)*nelms+c_loc) * catomw
            k = som1; decomp_cpools_vr(c,j,i_soil1  ) = tracer_conc_solid_passive(c,j,(k-1)*nelms+c_loc) * catomw
            k = som2; decomp_cpools_vr(c,j,i_soil2  ) = tracer_conc_solid_passive(c,j,(k-1)*nelms+c_loc) * catomw
            k = som3; decomp_cpools_vr(c,j,i_soil3  ) = tracer_conc_solid_passive(c,j,(k-1)*nelms+c_loc) * catomw

            k = lit1; decomp_npools_vr(c,j,i_met_lit) = tracer_conc_solid_passive(c,j,(k-1)*nelms+n_loc) * natomw
            k = lit2; decomp_npools_vr(c,j,i_cel_lit) = tracer_conc_solid_passive(c,j,(k-1)*nelms+n_loc) * natomw
            k = lit3; decomp_npools_vr(c,j,i_lig_lit) = tracer_conc_solid_passive(c,j,(k-1)*nelms+n_loc) * natomw
            k = cwd ; decomp_npools_vr(c,j,i_cwd    ) = tracer_conc_solid_passive(c,j,(k-1)*nelms+n_loc) * natomw
            k = som1; decomp_npools_vr(c,j,i_soil1  ) = tracer_conc_solid_passive(c,j,(k-1)*nelms+n_loc) * natomw
            k = som2; decomp_npools_vr(c,j,i_soil2  ) = tracer_conc_solid_passive(c,j,(k-1)*nelms+n_loc) * natomw
            k = som3; decomp_npools_vr(c,j,i_soil3  ) = tracer_conc_solid_passive(c,j,(k-1)*nelms+n_loc) * natomw
         enddo
      enddo


    end associate
  end subroutine assign_OM_CNpools
  !-------------------------------------------------------------------------------
  subroutine assign_nitrogen_hydroloss(bounds, num_soilc, filter_soilc, &
       tracerflux_vars, nitrogenflux_vars, betrtracer_vars)

    !
    ! !DESCRIPTION:
    ! feedback the nitrogen hydrological fluxes, this comes after tracer mass balance, so the flux is with the unit of gN/m2/s
    !
    ! !USES:
    use tracerfluxType      , only : tracerflux_type
    use BetrTracerType      , only : betrtracer_type
    use CNNitrogenFluxType  , only : nitrogenflux_type
    use clm_varcon          , only : natomw

    ! !ARGUMENTS:
    type(bounds_type)       , intent(in)    :: bounds            ! bounds
    integer                 , intent(in)    :: num_soilc         ! number of columns in column filter
    integer                 , intent(in)    :: filter_soilc(:)   ! column filter
    type(tracerflux_type)   , intent(in)    :: tracerflux_vars
    type(betrtracer_type)   , intent(in)    :: betrtracer_vars   ! betr configuration information
    type(nitrogenflux_type) , intent(inout) :: nitrogenflux_vars !

    ! !LOCAL VARIABLES:
    integer :: fc, c
    !get nitrogen leaching, and loss through surface runoff

    associate(                                      &
         id_trc_no3x => betrtracer_vars%id_trc_no3x &
         )

      do fc = 1, num_soilc
         c = filter_soilc(fc)
         nitrogenflux_vars%smin_no3_leached_col(c) = tracerflux_vars%tracer_flx_totleached_col(c,id_trc_no3x)*natomw
         nitrogenflux_vars%smin_no3_runoff_col(c)  = tracerflux_vars%tracer_flx_surfrun_col(c,id_trc_no3x)*natomw
      enddo

    end associate
  end subroutine assign_nitrogen_hydroloss

  !-------------------------------------------------------------------------------
  subroutine apply_plant_root_nuptake_prof(bounds, ubj, num_soilc, filter_soilc, &
       root_prof_col, plantsoilnutrientflux_vars)
    !
    ! !DESCRIPTION:
    ! nitroge uptake profile
    ! !USES:
    use PlantSoilnutrientFluxType, only : plantsoilnutrientflux_type
    implicit none
    ! !ARGUMENTS:
    type(bounds_type)                , intent(in)    :: bounds                                        ! bounds
    integer                          , intent(in)    :: ubj
    integer                          , intent(in)    :: num_soilc                                     ! number of columns in column filter
    integer                          , intent(in)    :: filter_soilc(:)                               ! column filter
    real(r8)                         , intent(in)    :: root_prof_col(bounds%begc:bounds%endc, 1:ubj) !
    type(plantsoilnutrientflux_type) , intent(inout) :: plantsoilnutrientflux_vars

    ! !LOCAL VARIABLES:
    integer :: fc, c, j

    associate(                                                                &
         plant_frootsc_col =>  plantsoilnutrientflux_vars%plant_frootsc_col , &
         plant_frootsc_vr  => plantsoilnutrientflux_vars%plant_frootsc_vr_col &
         )

      do j = 1, ubj
         do fc = 1, num_soilc
            c = filter_soilc(fc)
            plant_frootsc_vr(c, j) = root_prof_col(c,j)  * plant_frootsc_col(c)
         enddo
      enddo
    end associate

  end subroutine apply_plant_root_nuptake_prof

end module BGCCenturySubCoreMod
