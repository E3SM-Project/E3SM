module CNNitrogenStateType

#include "shr_assert.h"

  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_infnan_mod         , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use clm_varpar             , only : ndecomp_cascade_transitions, ndecomp_pools, nlevcan
  use clm_varpar             , only : nlevdecomp_full, nlevdecomp, crop_prog
  use clm_varcon             , only : spval, ispval, dzsoi_decomp, zisoi
  use landunit_varcon        , only : istcrop, istsoil 
  use clm_varctl             , only : use_nitrif_denitrif, use_vertsoilc, use_century_decomp
  use clm_varctl             , only : iulog, override_bgc_restart_mismatch_dump, spinup_state
  use decompMod              , only : bounds_type
  use pftvarcon              , only : npcropmin, nstor
  use CNDecompCascadeConType , only : decomp_cascade_con
  use VegetationPropertiesType         , only : veg_vp
  use abortutils             , only : endrun
  use spmdMod                , only : masterproc 
  use LandunitType           , only : lun_pp                
  use ColumnType             , only : col_pp                
  use VegetationType         , only : veg_pp
  use clm_varctl             , only : use_pflotran, pf_cmode
  use clm_varctl             , only : nu_com, use_crop
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
     procedure , public  :: DynamicPatchAdjustments
     procedure , public  :: DynamicColumnAdjustments
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
    use clm_varpar , only : ndecomp_cascade_transitions, ndecomp_pools
    use clm_varpar , only : nlevdecomp, nlevdecomp_full,crop_prog, nlevgrnd
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
    
    if (crop_prog) then
       this%grainn_patch(begp:endp) = spval
       call hist_addfld1d (fname='GRAINN', units='gN/m^2', &
            avgflag='A', long_name='grain N', &
            ptr_patch=this%grainn_patch, default='inactive')

       this%cropseedn_deficit_patch(begp:endp) = spval
       call hist_addfld1d (fname='CROPSEEDN_DEFICIT', units='gN/m^2', &
            avgflag='A', long_name='N used for crop seed that needs to be repaid', &
            ptr_patch=this%cropseedn_deficit_patch)
    end if

    this%leafn_patch(begp:endp) = spval
    call hist_addfld1d (fname='LEAFN', units='gN/m^2', &
         avgflag='A', long_name='leaf N', &
         ptr_patch=this%leafn_patch)

    this%leafn_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='LEAFN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='leaf N storage', &
         ptr_patch=this%leafn_storage_patch, default='inactive')

    this%leafn_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='LEAFN_XFER', units='gN/m^2', &
         avgflag='A', long_name='leaf N transfer', &
         ptr_patch=this%leafn_xfer_patch, default='inactive')

    this%frootn_patch(begp:endp) = spval
    call hist_addfld1d (fname='FROOTN', units='gN/m^2', &
         avgflag='A', long_name='fine root N', &
         ptr_patch=this%frootn_patch)

    this%frootn_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='FROOTN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='fine root N storage', &
         ptr_patch=this%frootn_storage_patch, default='inactive')

    this%frootn_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='FROOTN_XFER', units='gN/m^2', &
         avgflag='A', long_name='fine root N transfer', &
         ptr_patch=this%frootn_xfer_patch, default='inactive')

    this%livestemn_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMN', units='gN/m^2', &
         avgflag='A', long_name='live stem N', &
         ptr_patch=this%livestemn_patch)

    this%livestemn_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='live stem N storage', &
         ptr_patch=this%livestemn_storage_patch, default='inactive')

    this%livestemn_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMN_XFER', units='gN/m^2', &
         avgflag='A', long_name='live stem N transfer', &
         ptr_patch=this%livestemn_xfer_patch, default='inactive')

    this%deadstemn_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADSTEMN', units='gN/m^2', &
         avgflag='A', long_name='dead stem N', &
         ptr_patch=this%deadstemn_patch)

    this%deadstemn_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADSTEMN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='dead stem N storage', &
         ptr_patch=this%deadstemn_storage_patch, default='inactive')

    this%deadstemn_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADSTEMN_XFER', units='gN/m^2', &
         avgflag='A', long_name='dead stem N transfer', &
         ptr_patch=this%deadstemn_xfer_patch, default='inactive')

    this%livecrootn_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTN', units='gN/m^2', &
         avgflag='A', long_name='live coarse root N', &
         ptr_patch=this%livecrootn_patch)

    this%livecrootn_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='live coarse root N storage', &
         ptr_patch=this%livecrootn_storage_patch, default='inactive')

    this%livecrootn_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTN_XFER', units='gN/m^2', &
         avgflag='A', long_name='live coarse root N transfer', &
         ptr_patch=this%livecrootn_xfer_patch, default='inactive')

    this%deadcrootn_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADCROOTN', units='gN/m^2', &
         avgflag='A', long_name='dead coarse root N', &
         ptr_patch=this%deadcrootn_patch)

    this%deadcrootn_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADCROOTN_STORAGE', units='gN/m^2', &
         avgflag='A', long_name='dead coarse root N storage', &
         ptr_patch=this%deadcrootn_storage_patch, default='inactive')

    this%deadcrootn_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADCROOTN_XFER', units='gN/m^2', &
         avgflag='A', long_name='dead coarse root N transfer', &
         ptr_patch=this%deadcrootn_xfer_patch, default='inactive')

    this%retransn_patch(begp:endp) = spval
    call hist_addfld1d (fname='RETRANSN', units='gN/m^2', &
         avgflag='A', long_name='plant pool of retranslocated N', &
         ptr_patch=this%retransn_patch)

    this%npool_patch(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL', units='gN/m^2', &
         avgflag='A', long_name='temporary plant N pool', &
         ptr_patch=this%npool_patch, default='inactive')

    this%ntrunc_patch(begp:endp) = spval
    call hist_addfld1d (fname='PFT_NTRUNC', units='gN/m^2', &
         avgflag='A', long_name='pft-level sink for N truncation', &
         ptr_patch=this%ntrunc_patch, default='inactive')

    this%dispvegn_patch(begp:endp) = spval
    call hist_addfld1d (fname='DISPVEGN', units='gN/m^2', &
         avgflag='A', long_name='displayed vegetation nitrogen', &
         ptr_patch=this%dispvegn_patch)

    this%storvegn_patch(begp:endp) = spval
    call hist_addfld1d (fname='STORVEGN', units='gN/m^2', &
         avgflag='A', long_name='stored vegetation nitrogen', &
         ptr_patch=this%storvegn_patch)

    this%totvegn_patch(begp:endp) = spval
    call hist_addfld1d (fname='TOTVEGN', units='gN/m^2', &
         avgflag='A', long_name='total vegetation nitrogen', &
         ptr_patch=this%totvegn_patch)

    this%totpftn_patch(begp:endp) = spval
    call hist_addfld1d (fname='TOTPFTN', units='gN/m^2', &
         avgflag='A', long_name='total PFT-level nitrogen', &
         ptr_patch=this%totpftn_patch)

    this%npimbalance_patch(begp:endp) = spval
    call hist_addfld1d (fname='leaf_npimbalance', units='gN/gP', &
         avgflag='A', long_name='leaf np imbalance partial C partial P/partial C partial N', &
         ptr_patch=this%npimbalance_patch)
     
    !-------------------------------
    ! N state variables - native to column
    !-------------------------------

    if ( nlevdecomp_full > 1 ) then
       this%decomp_npools_vr_col(begc:endc,:,:) = spval
       this%decomp_npools_1m_col(begc:endc,:) = spval
    end if
    this%decomp_npools_col(begc:endc,:) = spval
    do l  = 1, ndecomp_pools
       if ( nlevdecomp_full > 1 ) then
          data2dptr => this%decomp_npools_vr_col(:,:,l)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'N_vr'
          longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' N (vertically resolved)'
          call hist_addfld2d (fname=fieldname, units='gN/m^3',  type2d='levdcmp', &
               avgflag='A', long_name=longname, &
               ptr_col=data2dptr)
       endif

       data1dptr => this%decomp_npools_col(:,l)
       fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'N'
       longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' N'
       call hist_addfld1d (fname=fieldname, units='gN/m^2', &
            avgflag='A', long_name=longname, &
            ptr_col=data1dptr)

       if ( nlevdecomp_full > 1 ) then
          data1dptr => this%decomp_npools_1m_col(:,l)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'N_1m'
          longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' N to 1 meter'
          call hist_addfld1d (fname=fieldname, units='gN/m^2', &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr, default = 'inactive')
       endif
    end do


    if ( nlevdecomp_full > 1 ) then

       this%sminn_col(begc:endc) = spval
       call hist_addfld1d (fname='SMINN', units='gN/m^2', &
            avgflag='A', long_name='soil mineral N', &
            ptr_col=this%sminn_col)

       this%totlitn_1m_col(begc:endc) = spval
       call hist_addfld1d (fname='TOTLITN_1m', units='gN/m^2', &
            avgflag='A', long_name='total litter N to 1 meter', &
            ptr_col=this%totlitn_1m_col, default='inactive')

       this%totsomn_1m_col(begc:endc) = spval
       call hist_addfld1d (fname='TOTSOMN_1m', units='gN/m^2', &
            avgflag='A', long_name='total soil organic matter N to 1 meter', &
            ptr_col=this%totsomn_1m_col, default='inactive')
    endif

    this%plant_n_buffer_patch(begp:endp) = spval
    call hist_addfld1d (fname='PLANTN_BUFFER', units='gN/m^2', &
            avgflag='A', long_name='plant nitrogen stored as buffer', &
            ptr_col=this%plant_n_buffer_patch,default='inactive')
    
    this%ntrunc_col(begc:endc) = spval
    call hist_addfld1d (fname='COL_NTRUNC', units='gN/m^2',  &
         avgflag='A', long_name='column-level sink for N truncation', &
         ptr_col=this%ntrunc_col, default='inactive')

    ! add suffix if number of soil decomposition depths is greater than 1
    if (nlevdecomp > 1) then
       vr_suffix = "_vr"
    else 
       vr_suffix = ""
    endif

    if (use_nitrif_denitrif .or. (use_pflotran .and. pf_cmode)) then
       this%smin_no3_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SMIN_NO3'//trim(vr_suffix), units='gN/m^3',  type2d='levdcmp', &
            avgflag='A', long_name='soil mineral NO3 (vert. res.)', &
            ptr_col=this%smin_no3_vr_col)

       this%smin_nh4_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SMIN_NH4'//trim(vr_suffix), units='gN/m^3',  type2d='levdcmp', &
            avgflag='A', long_name='soil mineral NH4 (vert. res.)', &
            ptr_col=this%smin_nh4_vr_col)

       ! pflotran
       if(use_pflotran .and. pf_cmode) then
          this%smin_nh4sorb_vr_col(begc:endc,:) = spval
          call hist_addfld_decomp (fname='SMIN_NH4SORB'//trim(vr_suffix), units='gN/m^3',  type2d='levdcmp', &
            avgflag='A', long_name='soil mineral NH4 absorbed (vert. res.)', &
            ptr_col=this%smin_nh4sorb_vr_col)
       end if

       if ( nlevdecomp_full > 1 ) then
          this%smin_no3_col(begc:endc) = spval
          call hist_addfld1d (fname='SMIN_NO3', units='gN/m^2', &
               avgflag='A', long_name='soil mineral NO3', &
               ptr_col=this%smin_no3_col)

          this%smin_nh4_col(begc:endc) = spval
          call hist_addfld1d (fname='SMIN_NH4', units='gN/m^2', &
               avgflag='A', long_name='soil mineral NH4', &
               ptr_col=this%smin_nh4_col)

          ! pflotran
          if(use_pflotran .and. pf_cmode) then
            this%smin_nh4sorb_col(begc:endc) = spval
            call hist_addfld1d (fname='SMIN_NH4SORB', units='gN/m^2', &
               avgflag='A', long_name='soil mineral NH4 absorbed', &
               ptr_col=this%smin_nh4sorb_col)
          end if
       end if

       this%sminn_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SMINN'//trim(vr_suffix), units='gN/m^3',  type2d='levdcmp', &
            avgflag='A', long_name='soil mineral N', &
            ptr_col=this%sminn_vr_col, default = 'inactive')
    else
       this%sminn_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SMINN'//trim(vr_suffix), units='gN/m^3',  type2d='levdcmp', &
            avgflag='A', long_name='soil mineral N', &
            ptr_col=this%sminn_vr_col)
    end if

    this%totlitn_col(begc:endc) = spval
    call hist_addfld1d (fname='TOTLITN', units='gN/m^2', &
         avgflag='A', long_name='total litter N', &
         ptr_col=this%totlitn_col)

    this%totsomn_col(begc:endc) = spval
    call hist_addfld1d (fname='TOTSOMN', units='gN/m^2', &
         avgflag='A', long_name='total soil organic matter N', &
         ptr_col=this%totsomn_col)

    this%totecosysn_col(begc:endc) = spval
    call hist_addfld1d (fname='TOTECOSYSN', units='gN/m^2', &
         avgflag='A', long_name='total ecosystem N but excl product pools', &
         ptr_col=this%totecosysn_col)

    this%totcoln_col(begc:endc) = spval
    call hist_addfld1d (fname='TOTCOLN', units='gN/m^2', &
         avgflag='A', long_name='total column-level N but excl product pools', &
         ptr_col=this%totcoln_col)

    this%seedn_grc(begg:endg) = spval
    call hist_addfld1d (fname='SEEDN_GRC', units='gN/m^2', &
         avgflag='A', long_name='pool for seeding new PFTs ', &
         ptr_gcell=this%seedn_grc, default='inactive')

    this%seedn_col(begc:endc) = spval
    call hist_addfld1d (fname='SEEDN', units='gN/m^2', &
         avgflag='A', long_name='pool for seeding new PFTs ', &
         ptr_col=this%seedn_col, default='inactive')

    this%prod10n_col(begc:endc) = spval
    call hist_addfld1d (fname='PROD10N', units='gN/m^2', &
         avgflag='A', long_name='10-yr wood product N', &
         ptr_col=this%prod10n_col, default='inactive')

    this%prod100n_col(begc:endc) = spval
    call hist_addfld1d (fname='PROD100N', units='gN/m^2', &
         avgflag='A', long_name='100-yr wood product N', &
         ptr_col=this%prod100n_col, default='inactive')

    this%prod1n_col(begc:endc) = spval
    call hist_addfld1d (fname='PROD1N', units='gN/m^2', &
         avgflag='A', long_name='1-yr crop product N', &
         ptr_col=this%prod1n_col, default='inactive')

    this%totprodn_col(begc:endc) = spval
    call hist_addfld1d (fname='TOTPRODN', units='gN/m^2', &
         avgflag='A', long_name='total wood product N', &
         ptr_col=this%totprodn_col, default='inactive')

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
    use clm_varpar     , only : crop_prog
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
    use clm_varctl          , only : spinup_mortality_factor
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

    call restartvar(ncid=ncid, flag=flag, varname='leafn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%leafn_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='leafn_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%leafn_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='leafn_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%leafn_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='frootn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%frootn_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='frootn_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%frootn_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='frootn_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%frootn_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livestemn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livestemn_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livestemn_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livestemn_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livestemn_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livestemn_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadstemn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadstemn_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadstemn_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadstemn_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadstemn_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadstemn_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livecrootn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livecrootn_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livecrootn_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livecrootn_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='livecrootn_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%livecrootn_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadcrootn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadcrootn_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadcrootn_storage', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadcrootn_storage_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='deadcrootn_xfer', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%deadcrootn_xfer_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='retransn', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%retransn_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='npool', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%npool_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='pft_ntrunc', xtype=ncd_double,  &
         dim1name='pft', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%ntrunc_patch) 

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag,  varname='grainn', xtype=ncd_double,  &
            dim1name='pft',    long_name='grain N', units='gN/m2', &
            interpinic_flag='interp', readvar=readvar, data=this%grainn_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='grainn_storage', xtype=ncd_double,  &
            dim1name='pft',    long_name='grain N storage', units='gN/m2', &
            interpinic_flag='interp', readvar=readvar, data=this%grainn_storage_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='grainn_xfer', xtype=ncd_double,  &
            dim1name='pft',    long_name='grain N transfer', units='gN/m2', &
            interpinic_flag='interp', readvar=readvar, data=this%grainn_xfer_patch)
    end if
    
    call restartvar(ncid=ncid, flag=flag,  varname='npimbalance_patch', xtype=ncd_double,  &
        dim1name='pft',    long_name='npimbalance_patch', units='-', &
        interpinic_flag='interp', readvar=readvar, data=this%npimbalance_patch)
     call restartvar(ncid=ncid, flag=flag,  varname='pnup_pfrootc_patch', xtype=ncd_double,  &
        dim1name='pft',    long_name='pnup_pfrootc_patch', units='-', &
        interpinic_flag='interp', readvar=readvar, data=this%pnup_pfrootc_patch)
    call restartvar(ncid=ncid, flag=flag,  varname='benefit_pgpp_pleafc_patch', xtype=ncd_double,  &
        dim1name='pft',    long_name='benefit_pgpp_pleafc_patch', units='-', &
        interpinic_flag='interp', readvar=readvar, data=this%benefit_pgpp_pleafc_patch)
 
    !--------------------------------
    ! column nitrogen state variables
    !--------------------------------

    ! sminn
    if (use_vertsoilc) then
       ptr2d => this%sminn_vr_col
       call restartvar(ncid=ncid, flag=flag, varname="sminn_vr", xtype=ncd_double,  &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)
    else
       ptr1d => this%sminn_vr_col(:,1)
       call restartvar(ncid=ncid, flag=flag, varname="sminn", xtype=ncd_double,  &
            dim1name='column', &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp' , readvar=readvar, data=ptr1d)
    end if
    if (flag=='read' .and. .not. readvar) then
       call endrun(msg='ERROR::'//trim(varname)//' is required on an initialization dataset'//&
            errMsg(__FILE__, __LINE__))
    end if

    ! decomposing N pools
    do k = 1, ndecomp_pools
       varname=trim(decomp_cascade_con%decomp_pool_name_restart(k))//'n'
       if (use_vertsoilc) then
          ptr2d => this%decomp_npools_vr_col(:,:,k)
          call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_vr", xtype=ncd_double, &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=ptr2d) 
       else
          ptr1d => this%decomp_npools_vr_col(:,1,k)
          call restartvar(ncid=ncid, flag=flag, varname=varname, xtype=ncd_double,  &
               dim1name='column', &
               long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp' , readvar=readvar, data=ptr1d)
       end if
       if (flag=='read' .and. .not. readvar) then
          call endrun(msg='ERROR:: '//trim(varname)//' is required on an initialization dataset'//&
               errMsg(__FILE__, __LINE__))
       end if
    end do

    if (use_vertsoilc) then
       ptr2d => this%ntrunc_vr_col
       call restartvar(ncid=ncid, flag=flag, varname="col_ntrunc_vr", xtype=ncd_double,  &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)
    else
       ptr1d => this%ntrunc_vr_col(:,1)
       call restartvar(ncid=ncid, flag=flag, varname="col_ntrunc", xtype=ncd_double,  &
            dim1name='column', &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp' , readvar=readvar, data=ptr1d)
    end if

    if (use_nitrif_denitrif .or. (use_pflotran .and. pf_cmode)) then
       ! smin_no3_vr
       if (use_vertsoilc) then
          ptr2d => this%smin_no3_vr_col(:,:)
          call restartvar(ncid=ncid, flag=flag, varname='smin_no3_vr', xtype=ncd_double, &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
       else
          ptr1d => this%smin_no3_vr_col(:,1)
          call restartvar(ncid=ncid, flag=flag, varname='smin_no3', xtype=ncd_double, &
               dim1name='column', &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=ptr1d)
       end if
       if (flag=='read' .and. .not. readvar) then
          call endrun(msg= 'ERROR:: smin_no3_vr'//' is required on an initialization dataset' )
       end if
    end if

    if (use_nitrif_denitrif .or. (use_pflotran .and. pf_cmode)) then
       ! smin_nh4
       if (use_vertsoilc) then
          ptr2d => this%smin_nh4_vr_col(:,:)
          call restartvar(ncid=ncid, flag=flag, varname='smin_nh4_vr', xtype=ncd_double, &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=ptr2d) 
       else
          ptr1d => this%smin_nh4_vr_col(:,1)
          call restartvar(ncid=ncid, flag=flag, varname='smin_nh4', xtype=ncd_double, &
               dim1name='column', &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=ptr1d)
       end if
       if (flag=='read' .and. .not. readvar) then
          call endrun(msg= 'ERROR:: smin_nh4_vr'//' is required on an initialization dataset' )
       end if
    end if

    ! pflotran: smin_nh4sorb
    if (use_pflotran .and. pf_cmode) then
       if (use_vertsoilc) then
          ptr2d => this%smin_nh4sorb_vr_col(:,:)
          call restartvar(ncid=ncid, flag=flag, varname='smin_nh4sorb_vr', xtype=ncd_double, &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
        else
          ptr1d => this%smin_nh4sorb_vr_col(:,1)
          call restartvar(ncid=ncid, flag=flag, varname='smin_nh4sorb', xtype=ncd_double, &
               dim1name='column', &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=ptr1d)
       end if
       if (flag=='read' .and. .not. readvar) then
          call endrun(msg= 'ERROR:: smin_nh4sorb_vr'//' is required on an initialization dataset' )
       end if
    end if

    ! Set the integrated sminn based on sminn_vr, as is done in CNSummaryMod (this may
    ! not be the most appropriate method or place to do this)

    this%sminn_col(bounds%begc:bounds%endc) = 0._r8
    do j = 1, nlevdecomp
       do c = bounds%begc, bounds%endc
          this%sminn_col(c) = &
               this%sminn_col(c) + &
               this%sminn_vr_col(c,j) * dzsoi_decomp(j)
       end do
    end do

    call restartvar(ncid=ncid, flag=flag, varname='totcoln', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%totcoln_col) 

    call restartvar(ncid=ncid, flag=flag, varname='seedn', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%seedn_col) 

    call restartvar(ncid=ncid, flag=flag, varname='prod10n', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%prod10n_col) 

    call restartvar(ncid=ncid, flag=flag, varname='prod100n', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%prod100n_col) 

    call restartvar(ncid=ncid, flag=flag, varname='prod1n', xtype=ncd_double,  &
         dim1name='column', long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%prod1n_col)

    ! decomp_cascade_state - the purpose of this is to check to make sure the bgc used 
    ! matches what the restart file was generated with.  
    ! add info about the SOM decomposition cascade

    if (use_century_decomp) then
       decomp_cascade_state = 1
    else
       decomp_cascade_state = 0
    end if
    ! add info about the nitrification / denitrification state
    if (use_nitrif_denitrif .or. (use_pflotran .and. pf_cmode)) then
       decomp_cascade_state = decomp_cascade_state + 10
    end if
    if (flag == 'write') itemp = decomp_cascade_state    
    call restartvar(ncid=ncid, flag=flag, varname='decomp_cascade_state', xtype=ncd_int,  &
         long_name='BGC of the model that wrote this restart file:' &
         // '  1s column: 0 = CLM-CN cascade, 1 = Century cascade;' &
         // ' 10s column: 0 = CLM-CN denitrification, 10 = Century denitrification', units='', &
         interpinic_flag='skip', readvar=readvar, data=itemp)
    if (flag=='read') then
       if (.not. readvar) then
          ! assume, for sake of backwards compatibility, that if decomp_cascade_state 
          ! is not in the restart file, then the current model state is the same as 
          ! the prior model state
          restart_file_decomp_cascade_state = decomp_cascade_state
          if ( masterproc ) write(iulog,*) ' CNRest: WARNING!  Restart file does not ' &
               // ' contain info on decomp_cascade_state used to generate the restart file.  '
          if ( masterproc ) write(iulog,*) '   Assuming the same as current setting: ', decomp_cascade_state
       else
          restart_file_decomp_cascade_state = itemp  
          if (decomp_cascade_state /= restart_file_decomp_cascade_state ) then
             if ( masterproc ) then
                write(iulog,*) 'CNRest: ERROR--the decomposition cascade differs between the current ' &
                     // ' model state and the model that wrote the restart file. '
                write(iulog,*) 'The model will be horribly out of equilibrium until after a lengthy spinup. '
                write(iulog,*) 'Stopping here since this is probably an error in configuring the run. '
                write(iulog,*) 'If you really wish to proceed, then override by setting '
                write(iulog,*) 'override_bgc_restart_mismatch_dump to .true. in the namelist'
                if ( .not. override_bgc_restart_mismatch_dump ) then
                   call endrun(msg= ' CNRest: Stopping. Decomposition cascade mismatch error.'//&
                        errMsg(__FILE__, __LINE__))
                endif
             endif
          endif
       end if
    end if

    !--------------------------------
    ! Spinup state
    !--------------------------------

    ! Do nothing for write
    ! Note that the call to write spinup_state out was done in CNCarbonStateType and
    ! cannot be called again because it will try to define the variable twice
    ! when the flag below is set to define
    if (flag == 'read') then
       call restartvar(ncid=ncid, flag=flag, varname='spinup_state', xtype=ncd_int,  &
            long_name='Spinup state of the model that wrote this restart file: ' &
            // ' 0 = normal model mode, 1 = AD spinup', units='', &
            interpinic_flag='copy', readvar=readvar,  data=idata)
       if (readvar) then
          restart_file_spinup_state = idata
       else
          ! assume, for sake of backwards compatibility, that if spinup_state is not in 
          ! the restart file then current model state is the same as prior model state
          restart_file_spinup_state = spinup_state
          if ( masterproc ) then
             write(iulog,*) ' WARNING!  Restart file does not contain info ' &
                  // ' on spinup state used to generate the restart file. '
             write(iulog,*) '   Assuming the same as current setting: ', spinup_state
          end if
       end if
    end if

    ! now compare the model and restart file spinup states, and either take the 
    ! model into spinup mode or out of it if they are not identical
    ! taking model out of spinup mode requires multiplying each decomposing pool 
    ! by the associated AD factor.
    ! putting model into spinup mode requires dividing each decomposing pool 
    ! by the associated AD factor.
    ! only allow this to occur on first timestep of model run.

    if (flag == 'read' .and. spinup_state /= restart_file_spinup_state ) then
       if (spinup_state == 0 .and. restart_file_spinup_state == 1 ) then
          if ( masterproc ) write(iulog,*) ' NitrogenStateType Restart: taking SOM pools out of AD spinup mode'
          exit_spinup = .true.
       else if (spinup_state == 1 .and. restart_file_spinup_state == 0 ) then
          if ( masterproc ) write(iulog,*) ' NitrogenStateType Restart: taking SOM pools into AD spinup mode'
          enter_spinup = .true.
       else
          call endrun(msg=' Error in entering/exiting spinup.  spinup_state ' &
               // ' != restart_file_spinup_state, but do not know what to do'//&
               errMsg(__FILE__, __LINE__))
       end if
       if (get_nstep() >= 2) then
          call endrun(msg=' Error in entering/exiting spinup - should occur only when nstep = 1'//&
               errMsg(__FILE__, __LINE__))
       endif
       do k = 1, ndecomp_pools
          do c = bounds%begc, bounds%endc
             do j = 1, nlevdecomp
	       if ( exit_spinup ) then
		 m = decomp_cascade_con%spinup_factor(k)
                 if (decomp_cascade_con%spinup_factor(k) > 1) m = m / cnstate_vars%scalaravg_col(c,j)
               else if ( enter_spinup ) then 
                 m = 1. / decomp_cascade_con%spinup_factor(k)
		 if (decomp_cascade_con%spinup_factor(k) > 1) m = m * cnstate_vars%scalaravg_col(c,j)
               end if 
               this%decomp_npools_vr_col(c,j,k) = this%decomp_npools_vr_col(c,j,k) * m
             end do
          end do
       end do
       do i = bounds%begp, bounds%endp
         if (exit_spinup) then
            m_veg = spinup_mortality_factor
         else if (enter_spinup) then
            m_veg = 1._r8 / spinup_mortality_factor
         end if
         this%deadstemn_patch(i)  = this%deadstemn_patch(i) * m_veg
         this%deadcrootn_patch(i) = this%deadcrootn_patch(i) * m_veg
       end do
    end if

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
    use clm_varpar    , only: nlevdecomp,ndecomp_cascade_transitions,ndecomp_pools
    use clm_varctl    , only: use_nitrif_denitrif
    use subgridAveMod , only: p2c
    use clm_varpar    , only: nlevdecomp_full
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
 
  !-----------------------------------------------------------------------
  subroutine DynamicPatchAdjustments( this, &
       bounds,                              &
       num_filterp_with_inactive,           &
       filterp_with_inactive,               &
       prior_weights,                       &
       patch_state_updater,                 &
       dwt_leafn_seed,                      &
       dwt_deadstemn_seed,                  &
       dwt_npool_seed,                      &
       conv_nflux,                          &
       dwt_frootn_to_litter,                &
       dwt_livecrootn_to_litter,            &
       dwt_deadcrootn_to_litter,            &
       prod10_nflux,                        &
       prod100_nflux,                       &
       crop_product_nflux                   &
       )
    !
    ! !DESCRIPTION:
    ! Adjust state variables and compute associated fluxes when patch areas change due to
    ! dynamic landuse
    !
    ! !USES:
    use pftvarcon          , only : pconv, pprod10, pprod100
    use dynPriorWeightsMod , only : prior_weights_type
    use landunit_varcon    , only : istsoil, istcrop
    use ComputeSeedMod   , only : ComputeSeedAmounts
    !
    ! !ARGUMENTS:
    class(nitrogenstate_type)      , intent(inout) :: this
    type(bounds_type)              , intent(in)    :: bounds
    integer                        , intent(in)    :: num_filterp_with_inactive
    integer                        , intent(in)    :: filterp_with_inactive(:)
    type(prior_weights_type)       , intent(in)    :: prior_weights
    type(patch_state_updater_type) , intent(in)    :: patch_state_updater
    real(r8)                       , intent(inout) :: dwt_leafn_seed           (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_deadstemn_seed       (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_npool_seed           (bounds%begp:)
    real(r8)                       , intent(inout) :: conv_nflux               (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_frootn_to_litter     (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_livecrootn_to_litter (bounds%begp:)
    real(r8)                       , intent(inout) :: dwt_deadcrootn_to_litter (bounds%begp:)
    real(r8)                       , intent(inout) :: prod10_nflux             (bounds%begp:)
    real(r8)                       , intent(inout) :: prod100_nflux            (bounds%begp:)
    real(r8)                       , intent(inout) :: crop_product_nflux       (bounds%begp:)
    !
    ! !LOCAL VARIABLES:
    integer                     :: begp, endp
    integer                     :: l, c, p
    logical                     :: old_weight_was_zero      (bounds%begp:bounds%endp)
    logical                     :: patch_grew               (bounds%begp:bounds%endp)

    ! The following are only set for growing patches:
    real(r8)                    :: seed_leafn_patch         (bounds%begp:bounds%endp)
    real(r8)                    :: seed_leafn_storage_patch (bounds%begp:bounds%endp)
    real(r8)                    :: seed_leafn_xfer_patch    (bounds%begp:bounds%endp)
    real(r8)                    :: seed_deadstemn_patch     (bounds%begp:bounds%endp)
    real(r8)                    :: seed_npool_patch         (bounds%begp:bounds%endp)

    real(r8)                    :: wood_product_nflux       (bounds%begp:bounds%endp)
    real(r8)                    :: deadstemn_patch_temp     (bounds%begp:bounds%endp)

    character(len=*), parameter :: subname = 'NStateDynamicPatchAdjustments'
    !-----------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    SHR_ASSERT_ALL((ubound(dwt_leafn_seed           ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_deadstemn_seed       ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_npool_seed           ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(conv_nflux               ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_frootn_to_litter     ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_livecrootn_to_litter ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(dwt_deadcrootn_to_litter ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(prod10_nflux             ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(prod100_nflux            ) == (/endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(crop_product_nflux       ) == (/endp/)), errMsg(__FILE__, __LINE__))
   
    old_weight_was_zero = patch_state_updater%old_weight_was_zero(bounds)
    patch_grew          = patch_state_updater%patch_grew(bounds)

    call ComputeSeedAmounts(bounds                                        , &
         species                    = CN_SPECIES_N                        , &
         leaf_patch                 = this%leafn_patch(begp:endp)           , &
         leaf_storage_patch         = this%leafn_storage_patch(begp:endp)   , &
         leaf_xfer_patch            = this%leafn_xfer_patch(begp:endp)      , &

         ! Calculations only needed for patches that grew:
         compute_here_patch         = patch_grew(begp:endp)               , &

         ! For patches that previously had zero area, ignore the current state for the
         ! sake of computing leaf proportions:
         ignore_current_state_patch = old_weight_was_zero(begp:endp)      , &

         seed_leaf_patch            = seed_leafn_patch(begp:endp)         , &
         seed_leaf_storage_patch    = seed_leafn_storage_patch(begp:endp) , &
         seed_leaf_xfer_patch       = seed_leafn_xfer_patch(begp:endp)    , &
         seed_deadstem_patch        = seed_deadstemn_patch(begp:endp)     , &
         pool_seed_param            = npool_seed_param                    , &
         pool_seed_patch            = seed_npool_patch(begp:endp))

    ! 1) LEAFN_PATCH
    call patch_state_updater%update_patch_state(            &
         bounds                                           , &
         num_filterp_with_inactive                        , &
         filterp_with_inactive                            , &
         var               = this%leafn_patch   (begp:endp) , &
         flux_out_grc_area = conv_nflux       (begp:endp) , &
         seed              = seed_leafn_patch (begp:endp) , &
         seed_addition     = dwt_leafn_seed   (begp:endp))

    ! 2) LEAFN_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                    &
         bounds                                                   , &
         num_filterp_with_inactive                                , &
         filterp_with_inactive                                    , &
         var               = this%leafn_storage_patch   (begp:endp) , &
         flux_out_grc_area = conv_nflux               (begp:endp) , &
         seed              = seed_leafn_storage_patch (begp:endp) , &
         seed_addition     = dwt_leafn_seed           (begp:endp))

    ! 3) LEAFN_XFER_PATCH
    call patch_state_updater%update_patch_state( &
         bounds                                                        , &
         num_filterp_with_inactive                                     , &
         filterp_with_inactive                                         , &
         var               = this%leafn_xfer_patch   (begp:endp), &
         flux_out_grc_area = conv_nflux            (begp:endp), &
         seed              = seed_leafn_xfer_patch (begp:endp), &
         seed_addition     = dwt_leafn_seed        (begp:endp))

    ! 4) FROOTN_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%frootn_patch(begp:endp)             , &
         flux_out_col_area = dwt_frootn_to_litter(begp:endp))

    ! 5) FROOTN_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%frootn_storage_patch(begp:endp)     , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 6) FROOTN_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%frootn_xfer_patch(begp:endp)        , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 7) LIVESTEMN_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livestemn_patch(begp:endp)          , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 8) LIVESTEMN_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livestemn_storage_patch(begp:endp)  , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 9) LIVESTEMN_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livestemn_xfer_patch(begp:endp)     , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 10) PROD10_NFLUX
    wood_product_nflux(begp:endp)      = 0._r8
    deadstemn_patch_temp(begp:endp)    = this%deadstemn_patch(begp:endp)
    call patch_state_updater%update_patch_state_partition_flux_by_type(     &
         bounds                                                           , &
         num_filterp_with_inactive                                        , &
         filterp_with_inactive                                            , &
         flux1_fraction_by_pft_type = pprod10                             , &
         var                        = deadstemn_patch_temp    (begp:endp) , &
         flux1_out                  = prod10_nflux            (begp:endp) , &
         flux2_out                  = wood_product_nflux      (begp:endp) , &
         seed                       = seed_deadstemn_patch    (begp:endp) )

    ! 11) PROD100_NFLUX
    wood_product_nflux(begp:endp)      = 0._r8
    deadstemn_patch_temp(begp:endp)    = this%deadstemn_patch(begp:endp)
    call patch_state_updater%update_patch_state_partition_flux_by_type(     &
         bounds                                                           , &
         num_filterp_with_inactive                                        , &
         filterp_with_inactive                                            , &
         flux1_fraction_by_pft_type = pprod100                            , &
         var                        = deadstemn_patch_temp    (begp:endp) , &
         flux1_out                  = prod100_nflux           (begp:endp) , &
         flux2_out                  = wood_product_nflux      (begp:endp) , &
         seed                       = seed_deadstemn_patch    (begp:endp))

    ! 12) DEADSTEMN_PATCH
    wood_product_nflux(begp:endp)      = 0._r8
    call patch_state_updater%update_patch_state_partition_flux_by_type(     &
         bounds                                                           , &
         num_filterp_with_inactive                                        , &
         filterp_with_inactive                                            , &
         flux1_fraction_by_pft_type = pconv                               , &
         var                        = this%deadstemn_patch   (begp:endp)    , &
         flux1_out                  = conv_nflux           (begp:endp)    , &
         flux2_out                  = wood_product_nflux   (begp:endp)    , &
         seed                       = seed_deadstemn_patch (begp:endp)    , &
         seed_addition              = dwt_deadstemn_seed   (begp:endp))

    ! 13) DEADSTEMN_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%deadstemn_storage_patch(begp:endp)  , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 14) DEADSTEMN_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%deadstemn_xfer_patch(begp:endp)     , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 15) LIVECROOTN_PATCH 
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livecrootn_patch(begp:endp)         , &
         flux_out_col_area = dwt_livecrootn_to_litter(begp:endp))

    ! 16) LIVECROOTN_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livecrootn_storage_patch(begp:endp) , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 17) LIVECROOTN_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%livecrootn_xfer_patch(begp:endp)    , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 18) DEADCROOTN_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%deadcrootn_patch(begp:endp)         , &
         flux_out_col_area = dwt_deadcrootn_to_litter(begp:endp))

    ! 19) DEADCROOTN_STORAGE_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%deadcrootn_storage_patch(begp:endp) , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 20) DEADCROOT_XFER_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%deadcrootn_xfer_patch(begp:endp)    , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 21) RETRANSN_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%retransn_patch(begp:endp)           , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 22) NTRUNC_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%ntrunc_patch(begp:endp)             , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 23) CPOOL_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%npool_patch(begp:endp)              , &
         flux_out_grc_area = conv_nflux(begp:endp))

    ! 24) DISPVEGN_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%dispvegn_patch(begp:endp))

    ! 25) STORVEGN_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%storvegn_patch(begp:endp))

    ! 26) TOTVEGN_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%totvegn_patch(begp:endp))

    ! 27) TOTPFTN_PATCH
    call patch_state_updater%update_patch_state(                      &
         bounds                                                     , &
         num_filterp_with_inactive                                  , &
         filterp_with_inactive                                      , &
         var               = this%totpftn_patch(begp:endp))

    if (use_crop) then
       ! This is a negative pool. So any deficit that we haven't repaid gets sucked out
       ! of the atmosphere.
       call patch_state_updater%update_patch_state(         &
            bounds                                        , &
            num_filterp_with_inactive                     , &
            filterp_with_inactive                         , &
            var = this%cropseedn_deficit_patch(begp:endp) , &
            flux_out_grc_area = conv_nflux(begp:endp))
    end if

    ! These fluxes are computed as negative quantities, but are expected to be positive,
    ! so flip the signs
    do p = begp,endp
       dwt_frootn_to_litter(p)     = -1._r8 * dwt_frootn_to_litter(p)
       dwt_livecrootn_to_litter(p) = -1._r8 * dwt_livecrootn_to_litter(p)
       dwt_deadcrootn_to_litter(p) = -1._r8 * dwt_deadcrootn_to_litter(p)
    end do

  end subroutine DynamicPatchAdjustments 

  !-----------------------------------------------------------------------
  subroutine DynamicColumnAdjustments( this, &
       bounds, clump_index, column_state_updater)
    !
    ! !DESCRIPTION:
    ! Adjust state variables and compute associated fluxes when patch areas change due to
    ! dynamic landuse
    !
    ! !USES:
    use dynPriorWeightsMod       , only : prior_weights_type
    use landunit_varcon          , only : istsoil, istcrop
    use dynColumnStateUpdaterMod , only : column_state_updater_type
    !
    ! !ARGUMENTS:
    class(nitrogenstate_type)       , intent(inout) :: this
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: clump_index
    type(column_state_updater_type) , intent(in)    :: column_state_updater
    !
    ! !LOCAL VARIABLES:
    integer                     :: l, j
    integer                     :: begc, endc
    real(r8)                    :: adjustment_one_level(bounds%begc:bounds%endc)

    character(len=*), parameter :: subname = 'NStateDynamicColumnAdjustments'
    !-----------------------------------------------------------------------

    begc = bounds%begc
    endc = bounds%endc

    !this%dyn_nbal_adjustments_col(begc:endc) = 0._r8

    do l = 1, ndecomp_pools
       do j = 1, nlevdecomp
          call column_state_updater%update_column_state_no_special_handling( &
               bounds      = bounds,                                         &
               clump_index = clump_index,                                    &
               var         = this%decomp_npools_vr_col(begc:endc, j, l),     &
               adjustment  = adjustment_one_level(begc:endc))

          this%dyn_nbal_adjustments_col(begc:endc) = &
               this%dyn_nbal_adjustments_col(begc:endc) + &
               adjustment_one_level(begc:endc) * dzsoi_decomp(j)

       end do
    end do

    do j = 1, nlevdecomp
       call column_state_updater%update_column_state_no_special_handling( &
            bounds      = bounds,                                         &
            clump_index = clump_index,                                    &
            var         = this%ntrunc_vr_col(begc:endc,j),     &
            adjustment  = adjustment_one_level(begc:endc))

       this%dyn_nbal_adjustments_col(begc:endc) = &
            this%dyn_nbal_adjustments_col(begc:endc) + &
            adjustment_one_level(begc:endc) * dzsoi_decomp(j)

       call column_state_updater%update_column_state_no_special_handling( &
           bounds      = bounds                          , &
           clump_index = clump_index                     , &
           var         = this%sminn_vr_col(begc:endc, j) , &
           adjustment  = adjustment_one_level(begc:endc))

       this%dyn_nbal_adjustments_col(begc:endc) = &
           this%dyn_nbal_adjustments_col(begc:endc) + &
           adjustment_one_level(begc:endc) * dzsoi_decomp(j)
    end do

  end subroutine DynamicColumnAdjustments

end module CNNitrogenStateType
