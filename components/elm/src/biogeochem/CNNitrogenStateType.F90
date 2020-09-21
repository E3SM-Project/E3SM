module CNNitrogenStateType

#include "shr_assert.h"

  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_infnan_mod         , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use elm_varpar             , only : ndecomp_cascade_transitions, ndecomp_pools, nlevcan
  use elm_varpar             , only : nlevdecomp_full, nlevdecomp, crop_prog
  use elm_varcon             , only : spval, ispval, dzsoi_decomp, zisoi
  use landunit_varcon        , only : istcrop, istsoil 
  use elm_varctl             , only : use_nitrif_denitrif, use_vertsoilc, use_century_decomp
  use elm_varctl             , only : iulog, override_bgc_restart_mismatch_dump, spinup_state
  use decompMod              , only : bounds_type
  use pftvarcon              , only : npcropmin, nstor
  use CNDecompCascadeConType , only : decomp_cascade_con
  use VegetationPropertiesType         , only : veg_vp
  use abortutils             , only : endrun
  use spmdMod                , only : masterproc 
  use LandunitType           , only : lun_pp                
  use ColumnType             , only : col_pp                
  use VegetationType         , only : veg_pp
  use elm_varctl             , only : use_pflotran, pf_cmode
  use elm_varctl             , only : nu_com, use_crop
  use dynPatchStateUpdaterMod, only : patch_state_updater_type               
  use SpeciesMod           , only : CN_SPECIES_N
  ! 
  ! !PUBLIC TYPES:
  implicit none
  save
  private

  real(r8) , parameter :: npool_seed_param     = 0.1_r8

  type, public :: nitrogenstate_type

     real(r8), pointer :: grainn_patch                 (:)     ! patch (gN/m2) grain N (crop)
     real(r8), pointer :: grainn_storage_patch         (:)     ! patch (gN/m2) grain N storage (crop)
     real(r8), pointer :: grainn_xfer_patch            (:)     ! patch (gN/m2) grain N transfer (crop)
     real(r8), pointer :: leafn_patch                  (:)     ! patch (gN/m2) leaf N 
     real(r8), pointer :: leafn_storage_patch          (:)     ! patch (gN/m2) leaf N storage
     real(r8), pointer :: leafn_xfer_patch             (:)     ! patch (gN/m2) leaf N transfer
     real(r8), pointer :: frootn_patch                 (:)     ! patch (gN/m2) fine root N
     real(r8), pointer :: frootn_storage_patch         (:)     ! patch (gN/m2) fine root N storage
     real(r8), pointer :: frootn_xfer_patch            (:)     ! patch (gN/m2) fine root N transfer
     real(r8), pointer :: livestemn_patch              (:)     ! patch (gN/m2) live stem N
     real(r8), pointer :: livestemn_storage_patch      (:)     ! patch (gN/m2) live stem N storage
     real(r8), pointer :: livestemn_xfer_patch         (:)     ! patch (gN/m2) live stem N transfer
     real(r8), pointer :: deadstemn_patch              (:)     ! patch (gN/m2) dead stem N
     real(r8), pointer :: deadstemn_storage_patch      (:)     ! patch (gN/m2) dead stem N storage
     real(r8), pointer :: deadstemn_xfer_patch         (:)     ! patch (gN/m2) dead stem N transfer
     real(r8), pointer :: livecrootn_patch             (:)     ! patch (gN/m2) live coarse root N
     real(r8), pointer :: livecrootn_storage_patch     (:)     ! patch (gN/m2) live coarse root N storage
     real(r8), pointer :: livecrootn_xfer_patch        (:)     ! patch (gN/m2) live coarse root N transfer
     real(r8), pointer :: deadcrootn_patch             (:)     ! patch (gN/m2) dead coarse root N
     real(r8), pointer :: deadcrootn_storage_patch     (:)     ! patch (gN/m2) dead coarse root N storage
     real(r8), pointer :: deadcrootn_xfer_patch        (:)     ! patch (gN/m2) dead coarse root N transfer
     real(r8), pointer :: retransn_patch               (:)     ! patch (gN/m2) plant pool of retranslocated N
     real(r8), pointer :: npool_patch                  (:)     ! patch (gN/m2) temporary plant N pool
     real(r8), pointer :: ntrunc_patch                 (:)     ! patch (gN/m2) pft-level sink for N truncation
     real(r8), pointer :: plant_n_buffer_patch         (:)     ! patch (gN/m2) pft-level abstract N storage
     real(r8), pointer :: plant_n_buffer_col           (:)     ! patch (gN/m2) col-level abstract N storage
     real(r8), pointer :: decomp_npools_vr_col         (:,:,:) ! col (gN/m3) vertically-resolved decomposing (litter, cwd, soil) N pools
     real(r8), pointer :: sminn_vr_col                 (:,:)   ! col (gN/m3) vertically-resolved soil mineral N
     real(r8), pointer :: ntrunc_vr_col                (:,:)   ! col (gN/m3) vertically-resolved column-level sink for N truncation

     ! NITRIF_DENITRIF
     real(r8), pointer :: smin_no3_vr_col              (:,:)   ! col (gN/m3) vertically-resolved soil mineral NO3
     real(r8), pointer :: smin_no3_col                 (:)     ! col (gN/m2) soil mineral NO3 pool
     real(r8), pointer :: smin_nh4_vr_col              (:,:)   ! col (gN/m3) vertically-resolved soil mineral NH4
     real(r8), pointer :: smin_nh4_col                 (:)     ! col (gN/m2) soil mineral NH4 pool

     ! wood product pools, for dynamic landcover
     real(r8), pointer :: cropseedn_deficit_patch      (:)     ! (gN/m2) pool for seeding new crop growth; this is a NEGATIVE term, indicating the amount of seed usage that needs to be repaid     
     real(r8), pointer :: seedn_grc                    (:)     ! (gN/m2) gridcell-level pool for seeding new PFTs via dynamic landcover
     real(r8), pointer :: seedn_col                    (:)     ! col (gN/m2) column-level pool for seeding new Patches
     real(r8), pointer :: prod1n_col                   (:)     ! col (gN/m2) crop product N pool, 1-year lifespan
     real(r8), pointer :: prod10n_col                  (:)     ! col (gN/m2) wood product N pool, 10-year lifespan
     real(r8), pointer :: prod100n_col                 (:)     ! col (gN/m2) wood product N pool, 100-year lifespan
     real(r8), pointer :: totprodn_col                 (:)     ! col (gN/m2) total wood product N
     real(r8), pointer :: dyn_nbal_adjustments_col     (:)     ! (gN/m2) adjustments to each column made in this timestep via dynamic column area adjustments

     ! summary (diagnostic) state variables, not involved in mass balance
     real(r8), pointer :: dispvegn_patch               (:)     ! patch (gN/m2) displayed veg nitrogen, excluding storage
     real(r8), pointer :: storvegn_patch               (:)     ! patch (gN/m2) stored vegetation nitrogen
     real(r8), pointer :: totvegn_patch                (:)     ! patch (gN/m2) total vegetation nitrogen
     real(r8), pointer :: totpftn_patch                (:)     ! patch (gN/m2) total pft-level nitrogen
     real(r8), pointer :: decomp_npools_col            (:,:)   ! col (gN/m2)  decomposing (litter, cwd, soil) N pools
     real(r8), pointer :: decomp_npools_1m_col         (:,:)   ! col (gN/m2)  diagnostic: decomposing (litter, cwd, soil) N pools to 1 meter
     real(r8), pointer :: sminn_col                    (:)     ! col (gN/m2) soil mineral N
     real(r8), pointer :: ntrunc_col                   (:)     ! col (gN/m2) column-level sink for N truncation
     real(r8), pointer :: cwdn_col                     (:)     ! col (gN/m2) Diagnostic: coarse woody debris N
     real(r8), pointer :: totlitn_col                  (:)     ! col (gN/m2) total litter nitrogen
     real(r8), pointer :: totsomn_col                  (:)     ! col (gN/m2) total soil organic matter nitrogen
     real(r8), pointer :: totlitn_1m_col               (:)     ! col (gN/m2) total litter nitrogen to 1 meter
     real(r8), pointer :: totsomn_1m_col               (:)     ! col (gN/m2) total soil organic matter nitrogen to 1 meter
     real(r8), pointer :: totecosysn_col               (:)     ! col (gN/m2) total ecosystem nitrogen, incl veg 
     real(r8), pointer :: totcoln_col                  (:)     ! col (gN/m2) total column nitrogen, incl veg
     real(r8), pointer :: totabgn_col                  (:)     ! col (gN/m2)
     real(r8), pointer :: totblgn_col                  (:)     ! col (gN/m2) total below ground nitrogen
     ! patch averaged to column variables 
     real(r8), pointer :: totvegn_col                  (:)     ! col (gN/m2) total vegetation nitrogen (p2c)
     real(r8), pointer :: totpftn_col                  (:)     ! col (gN/m2) total pft-level nitrogen (p2c)

     ! col balance checks
     real(r8), pointer :: begnb_patch                  (:)     ! patch nitrogen mass, beginning of time step (gN/m**2)
     real(r8), pointer :: endnb_patch                  (:)     ! patch nitrogen mass, end of time step (gN/m**2)
     real(r8), pointer :: errnb_patch                  (:)     ! patch nitrogen balance error for the timestep (gN/m**2)
     real(r8), pointer :: begnb_col                    (:)     ! col nitrogen mass, beginning of time step (gN/m**2)
     real(r8), pointer :: endnb_col                    (:)     ! col nitrogen mass, end of time step (gN/m**2)
     real(r8), pointer :: errnb_col                    (:)     ! colnitrogen balance error for the timestep (gN/m**2)
     real(r8), pointer :: begnb_grc                    (:)     ! grid cell nitrogen mass, beginning of time step (gN/m**2)
     real(r8), pointer :: endnb_grc                    (:)     ! grid cell nitrogen mass, end of time step (gN/m**2)
     real(r8), pointer :: errnb_grc                    (:)     ! grid cell nitrogen balance error for the timestep (gN/m**2)

     ! for newly-added coupled codes with pflotran (it should be included in total 'sminn' defined above when doing summation)
     real(r8), pointer :: smin_nh4sorb_vr_col          (:,:)   ! col (gN/m3) vertically-resolved soil mineral NH4 absorbed
     real(r8), pointer :: smin_nh4sorb_col             (:)     ! col (gN/m2) soil mineral NH4 pool absorbed

     real(r8), pointer :: plant_nbuffer_col            (:)     ! col plant nitrogen buffer, (gN/m2), used to exchange info with betr 

     real(r8), pointer :: totpftn_beg_col              (:)
     real(r8), pointer :: cwdn_beg_col                 (:)
     real(r8), pointer :: totlitn_beg_col              (:)
     real(r8), pointer :: totsomn_beg_col              (:)
     real(r8), pointer :: sminn_beg_col                (:)
     real(r8), pointer :: smin_no3_beg_col             (:)
     real(r8), pointer :: smin_nh4_beg_col             (:)
     real(r8), pointer :: totprodn_beg_col             (:)
     real(r8), pointer :: seedn_beg_col                (:)
     real(r8), pointer :: ntrunc_beg_col               (:)

     
     real(r8), pointer :: totpftn_end_col              (:)
     real(r8), pointer :: cwdn_end_col                 (:)
     real(r8), pointer :: totlitn_end_col              (:)
     real(r8), pointer :: totsomn_end_col              (:)
     real(r8), pointer :: sminn_end_col                (:)
     real(r8), pointer :: smin_no3_end_col             (:)
     real(r8), pointer :: smin_nh4_end_col             (:)
     real(r8), pointer :: totprodn_end_col             (:)
     real(r8), pointer :: seedn_end_col                (:)
     real(r8), pointer :: ntrunc_end_col               (:)

     ! for dynamic C/N/P allocation cost-benefit analysis
     real(r8), pointer :: npimbalance_patch                         (:)
     real(r8), pointer :: pnup_pfrootc_patch                        (:)
     real(r8), pointer :: ppup_pfrootc_patch                        (:)
     real(r8), pointer :: ptlai_pleafc_patch                        (:)
     
     real(r8), pointer :: ppsnsun_ptlai_patch                       (:)
     real(r8), pointer :: ppsnsun_pleafn_patch                      (:)
     real(r8), pointer :: ppsnsun_pleafp_patch                      (:)
     
     real(r8), pointer :: plmrsun_ptlai_patch                       (:)
     real(r8), pointer :: plmrsun_pleafn_patch                      (:)
     real(r8), pointer :: plaisun_ptlai_patch                       (:)
     
     real(r8), pointer :: ppsnsha_ptlai_patch                       (:)
     real(r8), pointer :: ppsnsha_pleafn_patch                      (:)
     real(r8), pointer :: ppsnsha_pleafp_patch                      (:)
     
     real(r8), pointer :: plmrsha_ptlai_patch                       (:)
     real(r8), pointer :: plmrsha_pleafn_patch                      (:)
     real(r8), pointer :: plaisha_ptlai_patch                       (:)
     
     real(r8), pointer :: benefit_pgpp_pleafc_patch                 (:)     ! partial gpp / partial leaf carbon (used by symbiotic n2 fixation and dynamic allocation)
     real(r8), pointer :: benefit_pgpp_pleafn_patch                 (:)     ! partial gpp / partial leaf nitrogen (used by phosphatase activity and dynamic allocation)
     real(r8), pointer :: benefit_pgpp_pleafp_patch                 (:)     ! partial gpp / partial leaf phosphorus (used by phosphatase activity and dynamic allocation)
     real(r8), pointer :: cost_pgpp_pfrootc_patch                   (:)     ! partial gpp /  partial fine root carbon (used by dynamic allocation)
     real(r8), pointer :: cost_plmr_pleafc_patch                    (:)     ! partial maintenance respiration /  partial leaf carbon (used by dynamic allocation)
     real(r8), pointer :: cost_plmr_pleafn_patch                    (:)     ! partial maintenance respiration /  partial leaf nitrogen (used by dynamic allocation)
     
     real(r8), pointer :: ppsn_ptlai_z                              (:,:)
     real(r8), pointer :: ppsn_pleafn_z                             (:,:)
     real(r8), pointer :: ppsn_pleafp_z                             (:,:)
     
     real(r8), pointer :: ppsn_ptlai_z_vcmax                        (:,:)
     real(r8), pointer :: ppsn_pleafn_z_vcmax                       (:,:)
     real(r8), pointer :: ppsn_pleafp_z_vcmax                       (:,:)
     
     real(r8), pointer :: ppsn_ptlai_z_jmax                         (:,:)
     real(r8), pointer :: ppsn_pleafn_z_jmax                        (:,:)
     real(r8), pointer :: ppsn_pleafp_z_jmax                        (:,:)
     
     real(r8), pointer :: ppsn_ptlai_z_tpu                          (:,:)
     real(r8), pointer :: ppsn_pleafn_z_tpu                         (:,:)
     real(r8), pointer :: ppsn_pleafp_z_tpu                         (:,:)
    
     real(r8), pointer :: plmr_ptlai_z                              (:,:)
     real(r8), pointer :: plmr_pleafn_z                             (:,:)

   contains

     procedure , public  :: Init   
     procedure , public  :: Restart
     procedure , public  :: SetValues
     procedure , public  :: ZeroDWT
     procedure , public  :: Summary
     procedure , private :: InitAllocate 
     procedure , private :: InitHistory  
     procedure , private :: InitCold     

  end type nitrogenstate_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds,                           &
       leafc_patch, leafc_storage_patch, frootc_patch, frootc_storage_patch, &
       deadstemc_patch, decomp_cpools_vr_col, decomp_cpools_col, decomp_cpools_1m_col)

    class(nitrogenstate_type)         :: this
    type(bounds_type) , intent(in)    :: bounds  
    real(r8)          , intent(in)    :: leafc_patch          (bounds%begp:)
    real(r8)          , intent(in)    :: leafc_storage_patch  (bounds%begp:)
    real(r8)          , intent(in)    :: frootc_patch         (bounds%begp:)
    real(r8)          , intent(in)    :: frootc_storage_patch (bounds%begp:)
    real(r8)          , intent(in)    :: deadstemc_patch      (bounds%begp:)
    real(r8)          , intent(in)    :: decomp_cpools_vr_col (bounds%begc:, 1:, 1:)
    real(r8)          , intent(in)    :: decomp_cpools_col    (bounds%begc:, 1:)
    real(r8)          , intent(in)    :: decomp_cpools_1m_col (bounds%begc:, 1:)

    call this%InitAllocate (bounds )

    call this%InitHistory (bounds)

    call this%InitCold ( bounds, leafc_patch, leafc_storage_patch, &
         frootc_patch, frootc_storage_patch, deadstemc_patch, &
         decomp_cpools_vr_col, decomp_cpools_col, decomp_cpools_1m_col)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !ARGUMENTS:
    class (nitrogenstate_type) :: this
    type(bounds_type) , intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begp,endp
    integer           :: begc,endc
    integer           :: begg,endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc
    begg = bounds%begg; endg = bounds%endg

    allocate(this%grainn_patch             (begp:endp))                   ; this%grainn_patch             (:)   = nan     
    allocate(this%grainn_storage_patch     (begp:endp))                   ; this%grainn_storage_patch     (:)   = nan
    allocate(this%grainn_xfer_patch        (begp:endp))                   ; this%grainn_xfer_patch        (:)   = nan     
    allocate(this%leafn_patch              (begp:endp))                   ; this%leafn_patch              (:)   = nan
    allocate(this%leafn_storage_patch      (begp:endp))                   ; this%leafn_storage_patch      (:)   = nan     
    allocate(this%leafn_xfer_patch         (begp:endp))                   ; this%leafn_xfer_patch         (:)   = nan     
    allocate(this%frootn_patch             (begp:endp))                   ; this%frootn_patch             (:)   = nan
    allocate(this%frootn_storage_patch     (begp:endp))                   ; this%frootn_storage_patch     (:)   = nan     
    allocate(this%frootn_xfer_patch        (begp:endp))                   ; this%frootn_xfer_patch        (:)   = nan     
    allocate(this%livestemn_patch          (begp:endp))                   ; this%livestemn_patch          (:)   = nan
    allocate(this%livestemn_storage_patch  (begp:endp))                   ; this%livestemn_storage_patch  (:)   = nan
    allocate(this%livestemn_xfer_patch     (begp:endp))                   ; this%livestemn_xfer_patch     (:)   = nan
    allocate(this%deadstemn_patch          (begp:endp))                   ; this%deadstemn_patch          (:)   = nan
    allocate(this%deadstemn_storage_patch  (begp:endp))                   ; this%deadstemn_storage_patch  (:)   = nan
    allocate(this%deadstemn_xfer_patch     (begp:endp))                   ; this%deadstemn_xfer_patch     (:)   = nan
    allocate(this%livecrootn_patch         (begp:endp))                   ; this%livecrootn_patch         (:)   = nan
    allocate(this%livecrootn_storage_patch (begp:endp))                   ; this%livecrootn_storage_patch (:)   = nan
    allocate(this%livecrootn_xfer_patch    (begp:endp))                   ; this%livecrootn_xfer_patch    (:)   = nan
    allocate(this%deadcrootn_patch         (begp:endp))                   ; this%deadcrootn_patch         (:)   = nan
    allocate(this%deadcrootn_storage_patch (begp:endp))                   ; this%deadcrootn_storage_patch (:)   = nan
    allocate(this%deadcrootn_xfer_patch    (begp:endp))                   ; this%deadcrootn_xfer_patch    (:)   = nan
    allocate(this%retransn_patch           (begp:endp))                   ; this%retransn_patch           (:)   = nan
    allocate(this%npool_patch              (begp:endp))                   ; this%npool_patch              (:)   = nan
    allocate(this%ntrunc_patch             (begp:endp))                   ; this%ntrunc_patch             (:)   = nan
    allocate(this%dispvegn_patch           (begp:endp))                   ; this%dispvegn_patch           (:)   = nan
    allocate(this%storvegn_patch           (begp:endp))                   ; this%storvegn_patch           (:)   = nan
    allocate(this%totvegn_patch            (begp:endp))                   ; this%totvegn_patch            (:)   = nan
    allocate(this%totpftn_patch            (begp:endp))                   ; this%totpftn_patch            (:)   = nan
    allocate(this%plant_n_buffer_patch    (begp:endp))                    ; this%plant_n_buffer_patch     (:)   = nan
    allocate(this%plant_n_buffer_col    (begc:endc))                      ; this%plant_n_buffer_col       (:)   = nan
    allocate(this%sminn_vr_col             (begc:endc,1:nlevdecomp_full)) ; this%sminn_vr_col             (:,:) = nan
    allocate(this%ntrunc_vr_col            (begc:endc,1:nlevdecomp_full)) ; this%ntrunc_vr_col            (:,:) = nan
    allocate(this%smin_no3_vr_col          (begc:endc,1:nlevdecomp_full)) ; this%smin_no3_vr_col          (:,:) = nan
    allocate(this%smin_nh4_vr_col          (begc:endc,1:nlevdecomp_full)) ; this%smin_nh4_vr_col          (:,:) = nan
    allocate(this%smin_no3_col             (begc:endc))                   ; this%smin_no3_col             (:)   = nan
    allocate(this%smin_nh4_col             (begc:endc))                   ; this%smin_nh4_col             (:)   = nan
    allocate(this%cwdn_col                 (begc:endc))                   ; this%cwdn_col                 (:)   = nan
    allocate(this%sminn_col                (begc:endc))                   ; this%sminn_col                (:)   = nan
    allocate(this%ntrunc_col               (begc:endc))                   ; this%ntrunc_col               (:)   = nan

    allocate(this%cropseedn_deficit_patch  (begp:endp))                   ; this%cropseedn_deficit_patch  (:)   = nan
    allocate(this%seedn_grc                (begg:endg))                   ; this%seedn_grc                (:)   = nan
    allocate(this%seedn_col                (begc:endc))                   ; this%seedn_col                (:)   = nan
    allocate(this%prod1n_col               (begc:endc))                   ; this%prod1n_col               (:)   = nan
    allocate(this%prod10n_col              (begc:endc))                   ; this%prod10n_col              (:)   = nan
    allocate(this%prod100n_col             (begc:endc))                   ; this%prod100n_col             (:)   = nan
    allocate(this%totprodn_col             (begc:endc))                   ; this%totprodn_col             (:)   = nan
    allocate(this%dyn_nbal_adjustments_col (begc:endc))                   ; this%dyn_nbal_adjustments_col (:)   = nan
    allocate(this%totlitn_col              (begc:endc))                   ; this%totlitn_col              (:)   = nan
    allocate(this%totsomn_col              (begc:endc))                   ; this%totsomn_col              (:)   = nan
    allocate(this%totlitn_1m_col           (begc:endc))                   ; this%totlitn_1m_col           (:)   = nan
    allocate(this%totsomn_1m_col           (begc:endc))                   ; this%totsomn_1m_col           (:)   = nan
    allocate(this%totecosysn_col           (begc:endc))                   ; this%totecosysn_col           (:)   = nan
    allocate(this%totcoln_col              (begc:endc))                   ; this%totcoln_col              (:)   = nan
    allocate(this%decomp_npools_col        (begc:endc,1:ndecomp_pools))   ; this%decomp_npools_col        (:,:) = nan
    allocate(this%decomp_npools_1m_col     (begc:endc,1:ndecomp_pools))   ; this%decomp_npools_1m_col     (:,:) = nan
    allocate(this%totpftn_col              (begc:endc))                   ; this%totpftn_col              (:)   = nan
    allocate(this%totvegn_col              (begc:endc))                   ; this%totvegn_col              (:)   = nan
    allocate(this%totabgn_col              (begc:endc))                   ; this%totabgn_col              (:)   = nan
    allocate(this%totblgn_col              (begc:endc))                   ; this%totblgn_col              (:)   = nan
    allocate(this%decomp_npools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools));
    this%decomp_npools_vr_col(:,:,:)= nan

    allocate(this%begnb_patch (begp:endp));     this%begnb_patch (:) =nan
    allocate(this%begnb_col   (begc:endc));     this%begnb_col   (:) =nan
    allocate(this%endnb_patch (begp:endp));     this%endnb_patch (:) =nan
    allocate(this%endnb_col   (begc:endc));     this%endnb_col   (:) =nan
    allocate(this%errnb_patch (begp:endp));     this%errnb_patch (:) =nan
    allocate(this%errnb_col   (begc:endc));     this%errnb_col   (:) =nan 

    allocate(this%begnb_grc   (begg:endg));     this%begnb_grc   (:) =nan
    allocate(this%endnb_grc   (begg:endg));     this%endnb_grc   (:) =nan
    allocate(this%errnb_grc   (begg:endg));     this%errnb_grc   (:) =nan
    
    allocate(this%totpftn_beg_col     (begc:endc))   ; this%totpftn_beg_col     (:) = nan
    allocate(this%cwdn_beg_col        (begc:endc))   ; this%cwdn_beg_col        (:) = nan
    allocate(this%totlitn_beg_col     (begc:endc))   ; this%totlitn_beg_col     (:) = nan
    allocate(this%totsomn_beg_col     (begc:endc))   ; this%totsomn_beg_col     (:) = nan
    allocate(this%sminn_beg_col       (begc:endc))   ; this%sminn_beg_col       (:) = nan
    allocate(this%smin_no3_beg_col    (begc:endc))   ; this%smin_no3_beg_col    (:) = nan
    allocate(this%smin_nh4_beg_col    (begc:endc))   ; this%smin_nh4_beg_col    (:) = nan
    allocate(this%totprodn_beg_col    (begc:endc))   ; this%totprodn_beg_col    (:) = nan
    allocate(this%seedn_beg_col       (begc:endc))   ; this%seedn_beg_col       (:) = nan
    allocate(this%ntrunc_beg_col      (begc:endc))   ; this%ntrunc_beg_col      (:) = nan
    
    allocate(this%totpftn_end_col     (begc:endc))   ; this%totpftn_end_col     (:) = nan
    allocate(this%cwdn_end_col        (begc:endc))   ; this%cwdn_end_col        (:) = nan
    allocate(this%totlitn_end_col     (begc:endc))   ; this%totlitn_end_col     (:) = nan
    allocate(this%totsomn_end_col     (begc:endc))   ; this%totsomn_end_col     (:) = nan
    allocate(this%sminn_end_col       (begc:endc))   ; this%sminn_end_col       (:) = nan
    allocate(this%smin_no3_end_col    (begc:endc))   ; this%smin_no3_end_col    (:) = nan
    allocate(this%smin_nh4_end_col    (begc:endc))   ; this%smin_nh4_end_col    (:) = nan
    allocate(this%totprodn_end_col    (begc:endc))   ; this%totprodn_end_col    (:) = nan
    allocate(this%seedn_end_col       (begc:endc))   ; this%seedn_end_col       (:) = nan
    allocate(this%ntrunc_end_col      (begc:endc))   ; this%ntrunc_end_col      (:) = nan

    ! for dynamic C/N/P allocation
    allocate(this%npimbalance_patch           (begp:endp)) ;             this%npimbalance_patch           (:) = nan
    allocate(this%pnup_pfrootc_patch          (begp:endp)) ;             this%pnup_pfrootc_patch          (:) = nan
    allocate(this%ppup_pfrootc_patch          (begp:endp)) ;             this%ppup_pfrootc_patch          (:) = nan
    allocate(this%ptlai_pleafc_patch          (begp:endp)) ;             this%ptlai_pleafc_patch          (:) = nan
    allocate(this%ppsnsun_ptlai_patch         (begp:endp)) ;             this%ppsnsun_ptlai_patch         (:) = nan
    allocate(this%ppsnsun_pleafn_patch        (begp:endp)) ;             this%ppsnsun_pleafn_patch        (:) = nan
    allocate(this%ppsnsun_pleafp_patch        (begp:endp)) ;             this%ppsnsun_pleafp_patch        (:) = nan
    allocate(this%plmrsun_ptlai_patch         (begp:endp)) ;             this%plmrsun_ptlai_patch         (:) = nan
    allocate(this%plmrsun_pleafn_patch        (begp:endp)) ;             this%plmrsun_pleafn_patch        (:) = nan
    allocate(this%plaisun_ptlai_patch         (begp:endp)) ;             this%plaisun_ptlai_patch         (:) = nan
    allocate(this%ppsnsha_ptlai_patch         (begp:endp)) ;             this%ppsnsha_ptlai_patch         (:) = nan
    allocate(this%ppsnsha_pleafn_patch        (begp:endp)) ;             this%ppsnsha_pleafn_patch        (:) = nan
    allocate(this%ppsnsha_pleafp_patch        (begp:endp)) ;             this%ppsnsha_pleafp_patch        (:) = nan
    allocate(this%plmrsha_ptlai_patch         (begp:endp)) ;             this%plmrsha_ptlai_patch         (:) = nan
    allocate(this%plmrsha_pleafn_patch        (begp:endp)) ;             this%plmrsha_pleafn_patch        (:) = nan
    allocate(this%plaisha_ptlai_patch         (begp:endp)) ;             this%plaisha_ptlai_patch         (:) = nan
    allocate(this%benefit_pgpp_pleafc_patch   (begp:endp)) ;             this%benefit_pgpp_pleafc_patch   (:) = nan
    allocate(this%benefit_pgpp_pleafn_patch   (begp:endp)) ;             this%benefit_pgpp_pleafn_patch   (:) = nan
    allocate(this%benefit_pgpp_pleafp_patch   (begp:endp)) ;             this%benefit_pgpp_pleafp_patch   (:) = nan
    allocate(this%cost_pgpp_pfrootc_patch     (begp:endp)) ;             this%cost_pgpp_pfrootc_patch     (:) = nan
    allocate(this%cost_plmr_pleafc_patch      (begp:endp)) ;             this%cost_plmr_pleafc_patch      (:) = nan
    allocate(this%cost_plmr_pleafn_patch      (begp:endp)) ;             this%cost_plmr_pleafn_patch      (:) = nan
    allocate(this%ppsn_ptlai_z                (begp:endp,1:nlevcan)) ;   this%ppsn_ptlai_z                (:,:) = nan
    allocate(this%ppsn_pleafn_z               (begp:endp,1:nlevcan)) ;   this%ppsn_pleafn_z               (:,:) = nan
    allocate(this%ppsn_pleafp_z               (begp:endp,1:nlevcan)) ;   this%ppsn_pleafp_z               (:,:) = nan
    allocate(this%ppsn_ptlai_z_vcmax          (begp:endp,1:nlevcan)) ;   this%ppsn_ptlai_z_vcmax          (:,:) = nan
    allocate(this%ppsn_pleafn_z_vcmax         (begp:endp,1:nlevcan)) ;   this%ppsn_pleafn_z_vcmax         (:,:) = nan
    allocate(this%ppsn_pleafp_z_vcmax         (begp:endp,1:nlevcan)) ;   this%ppsn_pleafp_z_vcmax         (:,:) = nan
    allocate(this%ppsn_ptlai_z_jmax           (begp:endp,1:nlevcan)) ;   this%ppsn_ptlai_z_jmax           (:,:) = nan
    allocate(this%ppsn_pleafn_z_jmax          (begp:endp,1:nlevcan)) ;   this%ppsn_pleafn_z_jmax          (:,:) = nan
    allocate(this%ppsn_pleafp_z_jmax          (begp:endp,1:nlevcan)) ;   this%ppsn_pleafp_z_jmax          (:,:) = nan
    allocate(this%ppsn_ptlai_z_tpu            (begp:endp,1:nlevcan)) ;   this%ppsn_ptlai_z_tpu            (:,:) = nan
    allocate(this%ppsn_pleafn_z_tpu           (begp:endp,1:nlevcan)) ;   this%ppsn_pleafn_z_tpu           (:,:) = nan
    allocate(this%ppsn_pleafp_z_tpu           (begp:endp,1:nlevcan)) ;   this%ppsn_pleafp_z_tpu           (:,:) = nan
    allocate(this%plmr_ptlai_z                (begp:endp,1:nlevcan)) ;   this%plmr_ptlai_z                (:,:) = nan
    allocate(this%plmr_pleafn_z               (begp:endp,1:nlevcan)) ;   this%plmr_pleafn_z               (:,:) = nan

    allocate(this%smin_nh4sorb_vr_col      (begc:endc,1:nlevdecomp_full)) ; this%smin_nh4sorb_vr_col      (:,:) = nan
    allocate(this%smin_nh4sorb_col         (begc:endc))                   ; this%smin_nh4sorb_col         (:)   = nan

    allocate(this%plant_nbuffer_col(begc:endc));this%plant_nbuffer_col(:) = nan

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! add history fields for all CN variables, always set as default='inactive'
    !
    ! !USES:
    use elm_varpar , only : ndecomp_cascade_transitions, ndecomp_pools
    use elm_varpar , only : nlevdecomp, nlevdecomp_full,crop_prog, nlevgrnd
    use histFileMod, only : hist_addfld1d, hist_addfld2d, hist_addfld_decomp 
    use decompMod  , only : bounds_type
    !
    ! !ARGUMENTS:
    class(nitrogenstate_type) :: this
    type(bounds_type)         , intent(in) :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer           :: k,l,ii,jj 
    character(10)     :: active
    character(8)      :: vr_suffix
    integer           :: begp,endp
    integer           :: begc,endc
    integer           :: begg,endg 
    character(24)     :: fieldname
    character(100)    :: longname
    real(r8), pointer :: data1dptr(:)   ! temp. pointer for slicing larger arrays
    real(r8), pointer :: data2dptr(:,:) ! temp. pointer for slicing larger arrays
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc
    begg = bounds%begg; endg = bounds%endg

    !-------------------------------
    ! N state variables - native to PFT
    !-------------------------------
    
    
    !-------------------------------
    ! N state variables - native to column
    !-------------------------------


  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, &
       leafc_patch, leafc_storage_patch, frootc_patch, frootc_storage_patch, &
       deadstemc_patch, decomp_cpools_vr_col, decomp_cpools_col, decomp_cpools_1m_col)
    !
    ! !DESCRIPTION:
    ! Initializes time varying variables used only in coupled carbon-nitrogen mode (CN):
    !
    ! !USES:
    use elm_varpar     , only : crop_prog
    use decompMod      , only : bounds_type
    use pftvarcon      , only : noveg, npcropmin
    !
    ! !ARGUMENTS:
    class(nitrogenstate_type)      :: this
    type(bounds_type) , intent(in) :: bounds  
    real(r8)          , intent(in) :: leafc_patch(bounds%begp:)
    real(r8)          , intent(in) :: leafc_storage_patch(bounds%begp:)
    real(r8)          , intent(in) :: frootc_patch(bounds%begp:)
    real(r8)          , intent(in) :: frootc_storage_patch(bounds%begp:)
    real(r8)          , intent(in) :: deadstemc_patch(bounds%begp:)
    real(r8)          , intent(in) :: decomp_cpools_vr_col(bounds%begc:,:,:)
    real(r8)          , intent(in) :: decomp_cpools_col(bounds%begc:,:)
    real(r8)          , intent(in) :: decomp_cpools_1m_col(bounds%begc:,:)
    !
    ! !LOCAL VARIABLES:
    integer :: fc,fp,g,l,c,p,j,k                       ! indices
    integer :: num_special_col                         ! number of good values in special_col filter
    integer :: num_special_patch                         ! number of good values in special_patch filter
    integer :: special_col   (bounds%endc-bounds%begc+1) ! special landunit filter - columns
    integer :: special_patch (bounds%endp-bounds%begp+1) ! special landunit filter - patches
    !------------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(leafc_patch)          == (/bounds%endp/)),                               errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(leafc_storage_patch)  == (/bounds%endp/)),                               errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(frootc_patch)         == (/bounds%endp/)),                               errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(frootc_storage_patch) == (/bounds%endp/)),                               errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(deadstemc_patch)      == (/bounds%endp/)),                               errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(decomp_cpools_col)    == (/bounds%endc,ndecomp_pools/)),                 errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(decomp_cpools_1m_col) == (/bounds%endc,ndecomp_pools/)),                 errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(decomp_cpools_vr_col) == (/bounds%endc,nlevdecomp_full,ndecomp_pools/)), errMsg(__FILE__, __LINE__))

    ! Set column filters

    num_special_patch = 0
    do p = bounds%begp,bounds%endp
       l = veg_pp%landunit(p)
       if (lun_pp%ifspecial(l)) then
          num_special_patch = num_special_patch + 1
          special_patch(num_special_patch) = p
       end if
    end do

    ! Set patch filters

    num_special_col = 0
    do c = bounds%begc, bounds%endc
       l = col_pp%landunit(c)
       if (lun_pp%ifspecial(l)) then
          num_special_col = num_special_col + 1
          special_col(num_special_col) = c
       end if
    end do

    !-------------------------------------------
    ! initialize pft-level variables
    !-------------------------------------------
    
    do p = bounds%begp,bounds%endp

       l = veg_pp%landunit(p)
       if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then       
          if (veg_pp%itype(p) == noveg) then
             this%leafn_patch(p) = 0._r8
             this%leafn_storage_patch(p) = 0._r8
          else
             this%leafn_patch(p)         = leafc_patch(p)         / veg_vp%leafcn(veg_pp%itype(p))
             this%leafn_storage_patch(p) = leafc_storage_patch(p) / veg_vp%leafcn(veg_pp%itype(p))
          end if

          this%leafn_xfer_patch(p)        = 0._r8
          if ( crop_prog )then
             this%grainn_patch(p)            = 0._r8
             this%grainn_storage_patch(p)    = 0._r8
             this%grainn_xfer_patch(p)       = 0._r8
             this%cropseedn_deficit_patch(p) = 0._r8
          end if
          this%frootn_patch(p)            = 0._r8
          this%frootn_storage_patch(p)    = 0._r8
          this%frootn_xfer_patch(p)       = 0._r8
          this%livestemn_patch(p)         = 0._r8
          this%livestemn_storage_patch(p) = 0._r8
          this%livestemn_xfer_patch(p)    = 0._r8

          ! tree types need to be initialized with some stem mass so that
          ! roughness length is not zero in canopy flux calculation

          if (veg_vp%woody(veg_pp%itype(p)) == 1._r8) then
             this%deadstemn_patch(p) = deadstemc_patch(p) / veg_vp%deadwdcn(veg_pp%itype(p))
          else
             this%deadstemn_patch(p) = 0._r8
          end if
          
          if (nu_com .ne. 'RD') then
              ! ECA competition calculate root NP uptake as a function of fine root biomass
              ! better to initialize root CNP pools with a non-zero value
              if (veg_pp%itype(p) .ne. noveg) then
                 this%frootn_patch(p) = frootc_patch(p) / veg_vp%frootcn(veg_pp%itype(p))
                 this%frootn_storage_patch(p) = frootc_storage_patch(p) / veg_vp%frootcn(veg_pp%itype(p))
              end if
          end if

          this%deadstemn_storage_patch(p)  = 0._r8
          this%deadstemn_xfer_patch(p)     = 0._r8
          this%livecrootn_patch(p)         = 0._r8
          this%livecrootn_storage_patch(p) = 0._r8
          this%livecrootn_xfer_patch(p)    = 0._r8
          this%deadcrootn_patch(p)         = 0._r8
          this%deadcrootn_storage_patch(p) = 0._r8
          this%deadcrootn_xfer_patch(p)    = 0._r8
          this%retransn_patch(p)           = 0._r8
          this%npool_patch(p)              = 0._r8
          if (nstor(veg_pp%itype(p)) .gt. 1e-6_r8) then 
              this%npool_patch(p)          = 10.0_r8
          end if
          this%ntrunc_patch(p)             = 0._r8
          this%dispvegn_patch(p)           = 0._r8
          this%storvegn_patch(p)           = 0._r8
          this%totvegn_patch(p)            = 0._r8
          this%totpftn_patch(p)            = 0._r8          
          this%plant_n_buffer_patch(p)     = 1._r8
       end if

       this%npimbalance_patch(p) = 0.0_r8
       this%pnup_pfrootc_patch(p) = 0.0_r8 
       this%benefit_pgpp_pleafc_patch(p) = 0.0_r8   
    end do

    !-------------------------------------------
    ! initialize column-level variables
    !-------------------------------------------

    do c = bounds%begc, bounds%endc
       l = col_pp%landunit(c)
       if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then

          ! column nitrogen state variables
          this%ntrunc_col(c) = 0._r8
          this%sminn_col(c) = 0._r8
          do j = 1, nlevdecomp
             do k = 1, ndecomp_pools
                this%decomp_npools_vr_col(c,j,k) = decomp_cpools_vr_col(c,j,k) / decomp_cascade_con%initial_cn_ratio(k)
             end do
             this%sminn_vr_col(c,j) = 0._r8
             this%ntrunc_vr_col(c,j) = 0._r8
          end do
          if ( nlevdecomp > 1 ) then
             do j = nlevdecomp+1, nlevdecomp_full
                do k = 1, ndecomp_pools
                   this%decomp_npools_vr_col(c,j,k) = 0._r8
                end do
                this%sminn_vr_col(c,j) = 0._r8
                this%ntrunc_vr_col(c,j) = 0._r8
             end do
          end if
          do k = 1, ndecomp_pools
             this%decomp_npools_col(c,k)    = decomp_cpools_col(c,k)    / decomp_cascade_con%initial_cn_ratio(k)
             this%decomp_npools_1m_col(c,k) = decomp_cpools_1m_col(c,k) / decomp_cascade_con%initial_cn_ratio(k)
          end do
          if (use_nitrif_denitrif .or. (use_pflotran .and. pf_cmode)) then
             do j = 1, nlevdecomp_full
                this%smin_nh4_vr_col(c,j) = 0._r8
                this%smin_no3_vr_col(c,j) = 0._r8
                if(use_pflotran .and. pf_cmode) then
                    this%smin_nh4sorb_vr_col(c,j) = 0._r8
                end if
             end do
             this%smin_nh4_col(c) = 0._r8
             this%smin_no3_col(c) = 0._r8
             if(use_pflotran .and. pf_cmode) then
                this%smin_nh4sorb_col(c) = 0._r8
             end if
          end if
          this%totlitn_col(c)    = 0._r8
          this%totsomn_col(c)    = 0._r8
          this%totlitn_1m_col(c) = 0._r8
          this%totsomn_1m_col(c) = 0._r8
          this%totecosysn_col(c) = 0._r8
          this%totcoln_col(c)    = 0._r8
          this%cwdn_col(c)       = 0._r8

          ! dynamic landcover state variables
          this%seedn_col(c)         = 0._r8
          this%prod1n_col(c)        = 0._r8
          this%prod10n_col(c)       = 0._r8
          this%prod100n_col(c)      = 0._r8
          this%totprodn_col(c)      = 0._r8
       end if
    end do

    ! now loop through special filters and explicitly set the variables that
    ! have to be in place for biogeophysics

    do fc = 1,num_special_col
       c = special_col(fc)

       this%seedn_col(c)    = 0._r8
       this%prod1n_col(c)   = 0._r8
       this%prod10n_col(c)  = 0._r8	  
       this%prod100n_col(c) = 0._r8	  
       this%totprodn_col(c) = 0._r8	  
    end do

    do g = bounds%begg, bounds%endg
       this%seedn_grc(g) = 0._r8
    end do

    ! initialize fields for special filters

    call this%SetValues (&
         num_patch=num_special_patch, filter_patch=special_patch, value_patch=0._r8, &
         num_column=num_special_col, filter_column=special_col, value_column=0._r8)

  end subroutine InitCold

  !-----------------------------------------------------------------------

  subroutine Restart ( this,  bounds, ncid, flag, cnstate_vars )
    !
    ! !DESCRIPTION: 
    ! Read/write CN restart data for carbon state
    !
    ! !USES:
    use shr_infnan_mod      , only : isnan => shr_infnan_isnan, nan => shr_infnan_nan, assignment(=)
    use clm_time_manager    , only : is_restart, get_nstep
    use elm_varctl          , only : spinup_mortality_factor
    use CNStateType         , only : cnstate_type
    use restUtilMod
    use ncdio_pio
    !
    ! !ARGUMENTS:
    class (nitrogenstate_type) :: this
    type(bounds_type)          , intent(in)    :: bounds 
    type(file_desc_t)          , intent(inout) :: ncid   
    type(cnstate_type)         , intent(in)    :: cnstate_vars
    character(len=*)           , intent(in)    :: flag   !'read' or 'write' or 'define'
    !
    ! !LOCAL VARIABLES:
    integer            :: i,j,k,l,c
    logical            :: readvar
    integer            :: idata
    logical            :: exit_spinup = .false.
    logical            :: enter_spinup = .false.
    real(r8)           :: m, m_veg          ! multiplier for the exit_spinup code
    real(r8), pointer  :: ptr2d(:,:) ! temp. pointers for slicing larger arrays
    real(r8), pointer  :: ptr1d(:)   ! temp. pointers for slicing larger arrays
    character(len=128) :: varname    ! temporary
    integer            :: itemp      ! temporary 
    integer , pointer  :: iptemp(:)  ! pointer to memory to be allocated
    ! spinup state as read from restart file, for determining whether to enter or exit spinup mode.
    integer            :: restart_file_spinup_state 
    ! flags for comparing the model and restart decomposition cascades
    integer            :: decomp_cascade_state, restart_file_decomp_cascade_state 
    !------------------------------------------------------------------------

    !--------------------------------
    ! patch nitrogen state variables
    !--------------------------------

 
    !--------------------------------
    ! column nitrogen state variables
    !--------------------------------

    !--------------------------------
    ! Spinup state
    !--------------------------------

    ! Do nothing for write
    ! Note that the call to write spinup_state out was done in CNCarbonStateType and
    ! cannot be called again because it will try to define the variable twice
    ! when the flag below is set to define
    ! now compare the model and restart file spinup states, and either take the 
    ! model into spinup mode or out of it if they are not identical
    ! taking model out of spinup mode requires multiplying each decomposing pool 
    ! by the associated AD factor.
    ! putting model into spinup mode requires dividing each decomposing pool 
    ! by the associated AD factor.
    ! only allow this to occur on first timestep of model run.


  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine SetValues ( this, &
       num_patch, filter_patch, value_patch, &
       num_column, filter_column, value_column)
    !
    ! !DESCRIPTION:
    ! Set nitrogen state variables
    !
    ! !ARGUMENTS:
    class (nitrogenstate_type) :: this
    integer , intent(in) :: num_patch
    integer , intent(in) :: filter_patch(:)
    real(r8), intent(in) :: value_patch
    integer , intent(in) :: num_column
    integer , intent(in) :: filter_column(:)
    real(r8), intent(in) :: value_column
    !
    ! !LOCAL VARIABLES:
    integer :: fi,i     ! loop index
    integer :: j,k      ! indices
    !------------------------------------------------------------------------

    do fi = 1,num_patch
       i = filter_patch(fi)

       this%leafn_patch(i)              = value_patch
       this%leafn_storage_patch(i)      = value_patch
       this%leafn_xfer_patch(i)         = value_patch
       this%frootn_patch(i)             = value_patch
       this%frootn_storage_patch(i)     = value_patch
       this%frootn_xfer_patch(i)        = value_patch
       this%livestemn_patch(i)          = value_patch
       this%livestemn_storage_patch(i)  = value_patch
       this%livestemn_xfer_patch(i)     = value_patch
       this%deadstemn_patch(i)          = value_patch
       this%deadstemn_storage_patch(i)  = value_patch
       this%deadstemn_xfer_patch(i)     = value_patch
       this%livecrootn_patch(i)         = value_patch
       this%livecrootn_storage_patch(i) = value_patch
       this%livecrootn_xfer_patch(i)    = value_patch
       this%deadcrootn_patch(i)         = value_patch
       this%deadcrootn_storage_patch(i) = value_patch
       this%deadcrootn_xfer_patch(i)    = value_patch
       this%retransn_patch(i)           = value_patch
       this%npool_patch(i)              = value_patch
       this%ntrunc_patch(i)             = value_patch
       this%dispvegn_patch(i)           = value_patch
       this%storvegn_patch(i)           = value_patch
       this%totvegn_patch(i)            = value_patch
       this%totpftn_patch(i)            = value_patch
    end do

    if ( crop_prog )then
       do fi = 1,num_patch
          i = filter_patch(fi)
          this%grainn_patch(i)            = value_patch
          this%grainn_storage_patch(i)    = value_patch
          this%grainn_xfer_patch(i)       = value_patch
          this%cropseedn_deficit_patch(i) = value_patch
       end do
    end if

    do fi = 1,num_column
       i = filter_column(fi)

       this%sminn_col(i)       = value_column
       this%ntrunc_col(i)      = value_column
       this%cwdn_col(i)        = value_column
       if (use_nitrif_denitrif .or. (use_pflotran .and. pf_cmode)) then
          this%smin_no3_col(i) = value_column
          this%smin_nh4_col(i) = value_column
          if(use_pflotran .and. pf_cmode) then
             this%smin_nh4sorb_col(i) = value_column
          end if
       end if
       this%totlitn_col(i)     = value_column
       this%totsomn_col(i)     = value_column
       this%totecosysn_col(i)  = value_column
       this%totcoln_col(i)     = value_column
       this%totsomn_1m_col(i)  = value_column
       this%totlitn_1m_col(i)  = value_column
    end do

    do j = 1,nlevdecomp_full
       do fi = 1,num_column
          i = filter_column(fi)
          this%sminn_vr_col(i,j)       = value_column
          this%ntrunc_vr_col(i,j)      = value_column
          if (use_nitrif_denitrif  .or. (use_pflotran .and. pf_cmode)) then
             this%smin_no3_vr_col(i,j) = value_column
             this%smin_nh4_vr_col(i,j) = value_column
             if(use_pflotran .and. pf_cmode) then
               this%smin_nh4sorb_vr_col(i,j) = value_column
             end if
          end if
       end do
    end do

    ! column and decomp_pools
    do k = 1, ndecomp_pools
       do fi = 1,num_column
          i = filter_column(fi)
          this%decomp_npools_col(i,k)    = value_column
          this%decomp_npools_1m_col(i,k) = value_column
       end do
    end do

    ! column levdecomp, and decomp_pools
    do j = 1,nlevdecomp_full
       do k = 1, ndecomp_pools
          do fi = 1,num_column
             i = filter_column(fi)
             this%decomp_npools_vr_col(i,j,k) = value_column
          end do
       end do
    end do

  end subroutine SetValues

  !-----------------------------------------------------------------------
  subroutine ZeroDwt( this, bounds )
    !
    ! !DESCRIPTION
    ! Initialize variables needed for dynamic land use.
    !
    ! !ARGUMENTS:
    class(nitrogenstate_type) :: this
    type(bounds_type), intent(in)  :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer  :: p          ! indices
    !-----------------------------------------------------------------------

    do p = bounds%begp,bounds%endp
       this%dispvegn_patch(p) = 0._r8
       this%storvegn_patch(p) = 0._r8
       this%totvegn_patch(p)  = 0._r8
       this%totpftn_patch(p)  = 0._r8
    end do

  end subroutine ZeroDwt

  !-----------------------------------------------------------------------
  subroutine Summary(this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
    !
    ! !USES:
    use elm_varpar    , only: nlevdecomp,ndecomp_cascade_transitions,ndecomp_pools
    use elm_varctl    , only: use_nitrif_denitrif
    use subgridAveMod , only: p2c
    use elm_varpar    , only: nlevdecomp_full
    !
    ! !ARGUMENTS:
    class (nitrogenstate_type) :: this
    type(bounds_type) , intent(in) :: bounds  
    integer           , intent(in) :: num_soilc       ! number of soil columns in filter
    integer           , intent(in) :: filter_soilc(:) ! filter for soil columns
    integer           , intent(in) :: num_soilp       ! number of soil patches in filter
    integer           , intent(in) :: filter_soilp(:) ! filter for soil patches
    !
    ! !LOCAL VARIABLES:
    integer  :: c,p,j,k,l   ! indices
    integer  :: fp,fc       ! lake filter indices
    real(r8) :: maxdepth    ! depth to integrate soil variables
    integer  :: nlev
    real(r8) :: cropseedn_deficit_col(bounds%begc:bounds%endc)
    !-----------------------------------------------------------------------

    do fp = 1,num_soilp
       p = filter_soilp(fp)

       ! displayed vegetation nitrogen, excluding storage (DISPVEGN)
       this%dispvegn_patch(p) = &
            this%leafn_patch(p)      + &
            this%frootn_patch(p)     + &
            this%livestemn_patch(p)  + &
            this%deadstemn_patch(p)  + &
            this%livecrootn_patch(p) + &
            this%deadcrootn_patch(p)
       
      ! stored vegetation nitrogen, including retranslocated N pool (STORVEGN)
      this%storvegn_patch(p) = &
           this%leafn_storage_patch(p)      + &
           this%frootn_storage_patch(p)     + &
           this%livestemn_storage_patch(p)  + &
           this%deadstemn_storage_patch(p)  + &
           this%livecrootn_storage_patch(p) + &
           this%deadcrootn_storage_patch(p) + &
           this%leafn_xfer_patch(p)         + &
           this%frootn_xfer_patch(p)        + &
           this%livestemn_xfer_patch(p)     + &
           this%deadstemn_xfer_patch(p)     + &
           this%livecrootn_xfer_patch(p)    + &
           this%deadcrootn_xfer_patch(p)    + &
           this%npool_patch(p)              + &
           this%retransn_patch(p)

      if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
         this%dispvegn_patch(p) = &
              this%dispvegn_patch(p) + &
              this%grainn_patch(p)

         this%storvegn_patch(p) = &
              this%storvegn_patch(p) + &
              this%grainn_storage_patch(p)     + &
              this%grainn_xfer_patch(p)
      end if

      ! total vegetation nitrogen (TOTVEGN)
      this%totvegn_patch(p) = &
           this%dispvegn_patch(p) + &
           this%storvegn_patch(p)

      ! total pft-level carbon (add pft_ntrunc)
      this%totpftn_patch(p) = &
           this%totvegn_patch(p) + &
           this%ntrunc_patch(p)

   end do

   call p2c(bounds, num_soilc, filter_soilc, &
        this%plant_n_buffer_patch(bounds%begp:bounds%endp), &
        this%plant_n_buffer_col(bounds%begc:bounds%endc))

   call p2c(bounds, num_soilc, filter_soilc, &
        this%totvegn_patch(bounds%begp:bounds%endp), &
        this%totvegn_col(bounds%begc:bounds%endc))

   call p2c(bounds, num_soilc, filter_soilc, &
        this%totpftn_patch(bounds%begp:bounds%endp), &
        this%totpftn_col(bounds%begc:bounds%endc))

   call p2c(bounds, num_soilc, filter_soilc, &
        this%cropseedn_deficit_patch(bounds%begp:bounds%endp), &
        cropseedn_deficit_col(bounds%begc:bounds%endc))

   ! vertically integrate NO3 NH4 N2O pools
   nlev = nlevdecomp
   if (use_pflotran .and. pf_cmode) nlev = nlevdecomp_full

   if (use_nitrif_denitrif .or. (use_pflotran .and. pf_cmode)) then
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         this%smin_no3_col(c) = 0._r8
         this%smin_nh4_col(c) = 0._r8
         if(use_pflotran .and. pf_cmode) then
            this%smin_nh4sorb_col(c) = 0._r8
         end if
      end do
      do j = 1, nlev
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            this%smin_no3_col(c) = &
                 this%smin_no3_col(c) + &
                 this%smin_no3_vr_col(c,j) * dzsoi_decomp(j)
            
            this%smin_nh4_col(c) = &
                 this%smin_nh4_col(c) + &
                 this%smin_nh4_vr_col(c,j) * dzsoi_decomp(j)
            if(use_pflotran .and. pf_cmode) then
               this%smin_nh4sorb_col(c) = &
                 this%smin_nh4sorb_col(c) + &
                 this%smin_nh4sorb_vr_col(c,j) * dzsoi_decomp(j)
            end if
          end do 
       end do

    end if

   ! vertically integrate each of the decomposing N pools
   do l = 1, ndecomp_pools
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         this%decomp_npools_col(c,l) = 0._r8
      end do
      do j = 1, nlev
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            this%decomp_npools_col(c,l) = &
                 this%decomp_npools_col(c,l) + &
                 this%decomp_npools_vr_col(c,j,l) * dzsoi_decomp(j)
         end do
      end do
   end do

   ! for vertically-resolved soil biogeochemistry, calculate some diagnostics of carbon pools to a given depth
   if ( nlevdecomp > 1) then

      do l = 1, ndecomp_pools
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            this%decomp_npools_1m_col(c,l) = 0._r8
         end do
      end do

      ! vertically integrate each of the decomposing n pools to 1 meter
      maxdepth = 1._r8
      do l = 1, ndecomp_pools
         do j = 1, nlevdecomp
            if ( zisoi(j) <= maxdepth ) then
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  this%decomp_npools_1m_col(c,l) = &
                       this%decomp_npools_1m_col(c,l) + &
                       this%decomp_npools_vr_col(c,j,l) * dzsoi_decomp(j)
               end do
            elseif ( zisoi(j-1) < maxdepth ) then
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  this%decomp_npools_1m_col(c,l) = &
                       this%decomp_npools_1m_col(c,l) + &
                       this%decomp_npools_vr_col(c,j,l) * (maxdepth - zisoi(j-1))
               end do
            endif
         end do
      end do
      
      ! total litter nitrogen to 1 meter (TOTLITN_1m)
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         this%totlitn_1m_col(c) = 0._r8
      end do
      do l = 1, ndecomp_pools
         if ( decomp_cascade_con%is_litter(l) ) then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               this%totlitn_1m_col(c) = &
                    this%totlitn_1m_col(c) + &
                    this%decomp_npools_1m_col(c,l)
            end do
         end if
      end do
      
      ! total soil organic matter nitrogen to 1 meter (TOTSOMN_1m)
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         this%totsomn_1m_col(c) = 0._r8
      end do
      do l = 1, ndecomp_pools
         if ( decomp_cascade_con%is_soil(l) ) then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               this%totsomn_1m_col(c) = &
                    this%totsomn_1m_col(c) + &
                    this%decomp_npools_1m_col(c,l)
            end do
         end if
      end do
      
   endif
   
   ! total litter nitrogen (TOTLITN)
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      this%totlitn_col(c)    = 0._r8
   end do
   do l = 1, ndecomp_pools
      if ( decomp_cascade_con%is_litter(l) ) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            this%totlitn_col(c) = &
                 this%totlitn_col(c) + &
                 this%decomp_npools_col(c,l)
         end do
      end if
   end do
   
   ! total soil organic matter nitrogen (TOTSOMN)
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      this%totsomn_col(c)    = 0._r8
   end do
   do l = 1, ndecomp_pools
      if ( decomp_cascade_con%is_soil(l) ) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            this%totsomn_col(c) = &
                 this%totsomn_col(c) + &
                 this%decomp_npools_col(c,l)
         end do
      end if
   end do
   
   ! total cwdn
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      this%cwdn_col(c) = 0._r8
   end do
   do l = 1, ndecomp_pools
      if ( decomp_cascade_con%is_cwd(l) ) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            this%cwdn_col(c) = &
                 this%cwdn_col(c) + &
                 this%decomp_npools_col(c,l)
         end do
      end if
   end do

   ! total sminn
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      this%sminn_col(c)      = 0._r8
   end do
   do j = 1, nlev
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         this%sminn_col(c) = &
              this%sminn_col(c) + &
              this%sminn_vr_col(c,j) * dzsoi_decomp(j)
      end do
   end do

   ! total col_ntrunc
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      this%ntrunc_col(c) = 0._r8
   end do
   do j = 1, nlev
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         this%ntrunc_col(c) = &
              this%ntrunc_col(c) + &
              this%ntrunc_vr_col(c,j) * dzsoi_decomp(j)
      end do
   end do

   do fc = 1,num_soilc
      c = filter_soilc(fc)

      ! total wood product nitrogen
      this%totprodn_col(c) = &
           this%prod1n_col(c) + &
           this%prod10n_col(c) + &
           this%prod100n_col(c)	 

      ! total ecosystem nitrogen, including veg (TOTECOSYSN)
      this%totecosysn_col(c) = &
           this%cwdn_col(c) + &
           this%totlitn_col(c) + &
           this%totsomn_col(c) + &
           this%sminn_col(c) + &
           this%totprodn_col(c) + &
           this%totvegn_col(c)

      ! total column nitrogen, including pft (TOTCOLN)
      this%totcoln_col(c) = &
           this%totpftn_col(c) + &
           this%cwdn_col(c) + &
           this%totlitn_col(c) + &
           this%totsomn_col(c) + &
           this%sminn_col(c) + &
           this%prod1n_col(c) + &
           this%ntrunc_col(c)+ &
           this%plant_n_buffer_col(c) + &
           cropseedn_deficit_col(c)
           
      this%totabgn_col (c) =  &
           this%totpftn_col(c) + &
           this%totprodn_col(c) + &
           this%seedn_col(c) + &
           this%ntrunc_col(c)+ &
           this%plant_n_buffer_col(c) 

      this%totblgn_col(c) = &
           this%cwdn_col(c) + &
           this%totlitn_col(c) + &
           this%totsomn_col(c) + &
           this%sminn_col(c) 
           
   end do

 end subroutine Summary
 

end module CNNitrogenStateType
