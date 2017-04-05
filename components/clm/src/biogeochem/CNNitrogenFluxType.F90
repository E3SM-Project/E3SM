module CNNitrogenFluxType

  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_infnan_mod         , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use clm_varpar             , only : ndecomp_cascade_transitions, ndecomp_pools
  use clm_varpar             , only : nlevdecomp_full, nlevdecomp, crop_prog
  use clm_varcon             , only : spval, ispval, dzsoi_decomp
  use decompMod              , only : bounds_type
  use clm_varctl             , only : use_nitrif_denitrif, use_vertsoilc
  use CNDecompCascadeConType , only : decomp_cascade_con
  use abortutils             , only : endrun
  use LandunitType           , only : lun
  use ColumnType             , only : col
  use PatchType              , only : pft
  !! bgc interface & pflotran:
  use clm_varctl             , only : use_bgc_interface, use_pflotran, pf_cmode, pf_hmode, use_vertsoilc
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  type, public :: nitrogenflux_type

     ! gap mortality fluxes
     real(r8), pointer :: m_leafn_to_litter_patch                   (:)     ! patch leaf N mortality (gN/m2/s)
     real(r8), pointer :: m_frootn_to_litter_patch                  (:)     ! patch fine root N mortality (gN/m2/s)
     real(r8), pointer :: m_leafn_storage_to_litter_patch           (:)     ! patch leaf N storage mortality (gN/m2/s)
     real(r8), pointer :: m_frootn_storage_to_litter_patch          (:)     ! patch fine root N storage mortality (gN/m2/s)
     real(r8), pointer :: m_livestemn_storage_to_litter_patch       (:)     ! patch live stem N storage mortality (gN/m2/s)
     real(r8), pointer :: m_deadstemn_storage_to_litter_patch       (:)     ! patch dead stem N storage mortality (gN/m2/s)
     real(r8), pointer :: m_livecrootn_storage_to_litter_patch      (:)     ! patch live coarse root N storage mortality (gN/m2/s)
     real(r8), pointer :: m_deadcrootn_storage_to_litter_patch      (:)     ! patch dead coarse root N storage mortality (gN/m2/s)
     real(r8), pointer :: m_leafn_xfer_to_litter_patch              (:)     ! patch leaf N transfer mortality (gN/m2/s)
     real(r8), pointer :: m_frootn_xfer_to_litter_patch             (:)     ! patch fine root N transfer mortality (gN/m2/s)
     real(r8), pointer :: m_livestemn_xfer_to_litter_patch          (:)     ! patch live stem N transfer mortality (gN/m2/s)
     real(r8), pointer :: m_deadstemn_xfer_to_litter_patch          (:)     ! patch dead stem N transfer mortality (gN/m2/s)
     real(r8), pointer :: m_livecrootn_xfer_to_litter_patch         (:)     ! patch live coarse root N transfer mortality (gN/m2/s)
     real(r8), pointer :: m_deadcrootn_xfer_to_litter_patch         (:)     ! patch dead coarse root N transfer mortality (gN/m2/s)
     real(r8), pointer :: m_livestemn_to_litter_patch               (:)     ! patch live stem N mortality (gN/m2/s)
     real(r8), pointer :: m_deadstemn_to_litter_patch               (:)     ! patch dead stem N mortality (gN/m2/s)
     real(r8), pointer :: m_livecrootn_to_litter_patch              (:)     ! patch live coarse root N mortality (gN/m2/s)
     real(r8), pointer :: m_deadcrootn_to_litter_patch              (:)     ! patch dead coarse root N mortality (gN/m2/s)
     real(r8), pointer :: m_retransn_to_litter_patch                (:)     ! patch retranslocated N pool mortality (gN/m2/s)

     ! harvest fluxes
     real(r8), pointer :: hrv_leafn_to_litter_patch                 (:)     ! patch leaf N harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_frootn_to_litter_patch                (:)     ! patch fine root N harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_leafn_storage_to_litter_patch         (:)     ! patch leaf N storage harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_frootn_storage_to_litter_patch        (:)     ! patch fine root N storage harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_livestemn_storage_to_litter_patch     (:)     ! patch live stem N storage harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_deadstemn_storage_to_litter_patch     (:)     ! patch dead stem N storage harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_livecrootn_storage_to_litter_patch    (:)     ! patch live coarse root N storage harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_deadcrootn_storage_to_litter_patch    (:)     ! patch dead coarse root N storage harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_leafn_xfer_to_litter_patch            (:)     ! patch leaf N transfer harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_frootn_xfer_to_litter_patch           (:)     ! patch fine root N transfer harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_livestemn_xfer_to_litter_patch        (:)     ! patch live stem N transfer harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_deadstemn_xfer_to_litter_patch        (:)     ! patch dead stem N transfer harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_livecrootn_xfer_to_litter_patch       (:)     ! patch live coarse root N transfer harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_deadcrootn_xfer_to_litter_patch       (:)     ! patch dead coarse root N transfer harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_livestemn_to_litter_patch             (:)     ! patch live stem N harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_deadstemn_to_prod10n_patch            (:)     ! patch dead stem N harvest to 10-year product pool (gN/m2/s)
     real(r8), pointer :: hrv_deadstemn_to_prod100n_patch           (:)     ! patch dead stem N harvest to 100-year product pool (gN/m2/s)
     real(r8), pointer :: hrv_livecrootn_to_litter_patch            (:)     ! patch live coarse root N harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_deadcrootn_to_litter_patch            (:)     ! patch dead coarse root N harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_retransn_to_litter_patch              (:)     ! patch retranslocated N pool harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_deadstemn_to_prod10n_col              (:)     ! col dead stem N harvest mortality to 10-year product pool (gN/m2/s)
     real(r8), pointer :: hrv_deadstemn_to_prod100n_col             (:)     ! col dead stem N harvest mortality to 100-year product pool (gN/m2/s)
     real(r8), pointer :: m_n_to_litr_met_fire_col                  (:,:)   ! col N from leaf, froot, xfer and storage N to litter labile N by fire (gN/m3/s)
     real(r8), pointer :: m_n_to_litr_cel_fire_col                  (:,:)   ! col N from leaf, froot, xfer and storage N to litter cellulose N by fire (gN/m3/s)
     real(r8), pointer :: m_n_to_litr_lig_fire_col                  (:,:)   ! col N from leaf, froot, xfer and storage N to litter lignin N by fire (gN/m3/s)
     real(r8), pointer :: harvest_n_to_litr_met_n_col               (:,:)   ! col N fluxes associated with harvest to litter metabolic pool (gN/m3/s)
     real(r8), pointer :: harvest_n_to_litr_cel_n_col               (:,:)   ! col N fluxes associated with harvest to litter cellulose pool (gN/m3/s)
     real(r8), pointer :: harvest_n_to_litr_lig_n_col               (:,:)   ! col N fluxes associated with harvest to litter lignin pool (gN/m3/s)
     real(r8), pointer :: harvest_n_to_cwdn_col                     (:,:)   ! col N fluxes associated with harvest to CWD pool (gN/m3/s)
     ! crop harvest
     real(r8), pointer :: hrv_leafn_to_prod1n_patch                 (:)     ! crop leaf N harvested (gN/m2/s)
     real(r8), pointer :: hrv_livestemn_to_prod1n_patch             (:)     ! crop stem N harvested (gN/m2/s)
     real(r8), pointer :: hrv_grainn_to_prod1n_patch                (:)     ! crop grain N harvested (gN/m2/s)
     real(r8), pointer :: hrv_cropn_to_prod1n_patch                 (:)     ! total amount of crop N harvested (gN/m2/s)
     real(r8), pointer :: hrv_cropn_to_prod1n_col                   (:)     ! crop N harvest mortality to 1-yr product pool (gN/m2/s)

     ! fire N fluxes
     real(r8), pointer :: m_decomp_npools_to_fire_vr_col            (:,:,:) ! col vertically-resolved decomposing N fire loss (gN/m3/s)
     real(r8), pointer :: m_decomp_npools_to_fire_col               (:,:)   ! col vertically-integrated (diagnostic) decomposing N fire loss (gN/m2/s)
     real(r8), pointer :: m_leafn_to_fire_patch                     (:)     ! patch (gN/m2/s) fire N emissions from leafn
     real(r8), pointer :: m_leafn_storage_to_fire_patch             (:)     ! patch (gN/m2/s) fire N emissions from leafn_storage
     real(r8), pointer :: m_leafn_xfer_to_fire_patch                (:)     ! patch (gN/m2/s) fire N emissions from leafn_xfer
     real(r8), pointer :: m_livestemn_to_fire_patch                 (:)     ! patch (gN/m2/s) fire N emissions from livestemn
     real(r8), pointer :: m_livestemn_storage_to_fire_patch         (:)     ! patch (gN/m2/s) fire N emissions from livestemn_storage
     real(r8), pointer :: m_livestemn_xfer_to_fire_patch            (:)     ! patch (gN/m2/s) fire N emissions from livestemn_xfer
     real(r8), pointer :: m_deadstemn_to_fire_patch                 (:)     ! patch (gN/m2/s) fire N emissions from deadstemn
     real(r8), pointer :: m_deadstemn_storage_to_fire_patch         (:)     ! patch (gN/m2/s) fire N emissions from deadstemn_storage
     real(r8), pointer :: m_deadstemn_xfer_to_fire_patch            (:)     ! patch (gN/m2/s) fire N emissions from deadstemn_xfer
     real(r8), pointer :: m_frootn_to_fire_patch                    (:)     ! patch (gN/m2/s) fire N emissions from frootn
     real(r8), pointer :: m_frootn_storage_to_fire_patch            (:)     ! patch (gN/m2/s) fire N emissions from frootn_storage
     real(r8), pointer :: m_frootn_xfer_to_fire_patch               (:)     ! patch (gN/m2/s) fire N emissions from frootn_xfer
     real(r8), pointer :: m_livecrootn_to_fire_patch                (:)     ! patch (gN/m2/s) fire N emissions from m_livecrootn_to_fire
     real(r8), pointer :: m_livecrootn_storage_to_fire_patch        (:)     ! patch (gN/m2/s) fire N emissions from livecrootn_storage
     real(r8), pointer :: m_livecrootn_xfer_to_fire_patch           (:)     ! patch (gN/m2/s) fire N emissions from livecrootn_xfer
     real(r8), pointer :: m_deadcrootn_to_fire_patch                (:)     ! patch (gN/m2/s) fire N emissions from deadcrootn
     real(r8), pointer :: m_deadcrootn_storage_to_fire_patch        (:)     ! patch (gN/m2/s) fire N emissions from deadcrootn_storage
     real(r8), pointer :: m_deadcrootn_xfer_to_fire_patch           (:)     ! patch (gN/m2/s) fire N emissions from deadcrootn_xfer
     real(r8), pointer :: m_retransn_to_fire_patch                  (:)     ! patch (gN/m2/s) fire N emissions from retransn
     real(r8), pointer :: m_leafn_to_litter_fire_patch              (:)     ! patch (gN/m2/s) from leafn to litter N  due to fire
     real(r8), pointer :: m_leafn_storage_to_litter_fire_patch      (:)     ! patch (gN/m2/s) from leafn_storage to litter N  due to fire
     real(r8), pointer :: m_leafn_xfer_to_litter_fire_patch         (:)     ! patch (gN/m2/s) from leafn_xfer to litter N  due to fire
     real(r8), pointer :: m_livestemn_to_litter_fire_patch          (:)     ! patch (gN/m2/s) from livestemn to litter N  due to fire
     real(r8), pointer :: m_livestemn_storage_to_litter_fire_patch  (:)     ! patch (gN/m2/s) from livestemn_storage to litter N  due to fire
     real(r8), pointer :: m_livestemn_xfer_to_litter_fire_patch     (:)     ! patch (gN/m2/s) from livestemn_xfer to litter N  due to fire
     real(r8), pointer :: m_livestemn_to_deadstemn_fire_patch       (:)     ! patch (gN/m2/s) from livestemn to deadstemn N  due to fire
     real(r8), pointer :: m_deadstemn_to_litter_fire_patch          (:)     ! patch (gN/m2/s) from deadstemn to litter N  due to fire
     real(r8), pointer :: m_deadstemn_storage_to_litter_fire_patch  (:)     ! patch (gN/m2/s) from deadstemn_storage to litter N  due to fire
     real(r8), pointer :: m_deadstemn_xfer_to_litter_fire_patch     (:)     ! patch (gN/m2/s) from deadstemn_xfer to litter N  due to fire
     real(r8), pointer :: m_frootn_to_litter_fire_patch             (:)     ! patch (gN/m2/s) from frootn to litter N  due to fire
     real(r8), pointer :: m_frootn_storage_to_litter_fire_patch     (:)     ! patch (gN/m2/s) from frootn_storage to litter N  due to fire
     real(r8), pointer :: m_frootn_xfer_to_litter_fire_patch        (:)     ! patch (gN/m2/s) from frootn_xfer to litter N  due to fire
     real(r8), pointer :: m_livecrootn_to_litter_fire_patch         (:)     ! patch (gN/m2/s) from livecrootn to litter N  due to fire
     real(r8), pointer :: m_livecrootn_storage_to_litter_fire_patch (:)     ! patch (gN/m2/s) from livecrootn_storage to litter N  due to fire
     real(r8), pointer :: m_livecrootn_xfer_to_litter_fire_patch    (:)     ! patch (gN/m2/s) from livecrootn_xfer to litter N  due to fire
     real(r8), pointer :: m_livecrootn_to_deadcrootn_fire_patch     (:)     ! patch (gN/m2/s) from livecrootn_xfer to deadcrootn due to fire
     real(r8), pointer :: m_deadcrootn_to_litter_fire_patch         (:)     ! patch (gN/m2/s) from deadcrootn to deadcrootn due to fire
     real(r8), pointer :: m_deadcrootn_storage_to_litter_fire_patch (:)     ! patch (gN/m2/s) from deadcrootn_storage to deadcrootn due to fire
     real(r8), pointer :: m_deadcrootn_xfer_to_litter_fire_patch    (:)     ! patch (gN/m2/s) from deadcrootn_xfer to deadcrootn due to fire
     real(r8), pointer :: m_retransn_to_litter_fire_patch           (:)     ! patch (gN/m2/s) from retransn to deadcrootn due to fire
     real(r8), pointer :: fire_nloss_patch                          (:)     ! patch total pft-level fire N loss (gN/m2/s)
     real(r8), pointer :: fire_nloss_col                            (:)     ! col total column-level fire N loss (gN/m2/s)
     real(r8), pointer :: fire_decomp_nloss_col                     (:)     ! col fire N loss from decomposable pools (gN/m2/s)
     real(r8), pointer :: fire_nloss_p2c_col                        (:)     ! col patch2col column-level fire N loss (gN/m2/s) (p2c)
     real(r8), pointer :: fire_mortality_n_to_cwdn_col              (:,:)   ! col N fluxes associated with fire mortality to CWD pool (gN/m3/s)

     ! phenology fluxes from transfer pool
     real(r8), pointer :: grainn_xfer_to_grainn_patch               (:)     ! patch grain N growth from storage for prognostic crop model (gN/m2/s)
     real(r8), pointer :: leafn_xfer_to_leafn_patch                 (:)     ! patch leaf N growth from storage (gN/m2/s)
     real(r8), pointer :: frootn_xfer_to_frootn_patch               (:)     ! patch fine root N growth from storage (gN/m2/s)
     real(r8), pointer :: livestemn_xfer_to_livestemn_patch         (:)     ! patch live stem N growth from storage (gN/m2/s)
     real(r8), pointer :: deadstemn_xfer_to_deadstemn_patch         (:)     ! patch dead stem N growth from storage (gN/m2/s)
     real(r8), pointer :: livecrootn_xfer_to_livecrootn_patch       (:)     ! patch live coarse root N growth from storage (gN/m2/s)
     real(r8), pointer :: deadcrootn_xfer_to_deadcrootn_patch       (:)     ! patch dead coarse root N growth from storage (gN/m2/s)

     ! litterfall fluxes
     real(r8), pointer :: livestemn_to_litter_patch                 (:)     ! patch livestem N to litter (gN/m2/s)
     real(r8), pointer :: grainn_to_food_patch                      (:)     ! patch grain N to food for prognostic crop (gN/m2/s)
     real(r8), pointer :: leafn_to_litter_patch                     (:)     ! patch leaf N litterfall (gN/m2/s)
     real(r8), pointer :: leafn_to_retransn_patch                   (:)     ! patch leaf N to retranslocated N pool (gN/m2/s)
     real(r8), pointer :: frootn_to_retransn_patch                  (:)     ! patch fine root N to retranslocated N pool (gN/m2/s)
     real(r8), pointer :: frootn_to_litter_patch                    (:)     ! patch fine root N litterfall (gN/m2/s)

     ! allocation fluxes
     real(r8), pointer :: retransn_to_npool_patch                   (:)     ! patch deployment of retranslocated N (gN/m2/s)
     real(r8), pointer :: sminn_to_npool_patch                      (:)     ! patch deployment of soil mineral N uptake (gN/m2/s)
     real(r8), pointer :: npool_to_grainn_patch                     (:)     ! patch allocation to grain N for prognostic crop (gN/m2/s)
     real(r8), pointer :: npool_to_grainn_storage_patch             (:)     ! patch allocation to grain N storage for prognostic crop (gN/m2/s)
     real(r8), pointer :: npool_to_leafn_patch                      (:)     ! patch allocation to leaf N (gN/m2/s)
     real(r8), pointer :: npool_to_leafn_storage_patch              (:)     ! patch allocation to leaf N storage (gN/m2/s)
     real(r8), pointer :: npool_to_frootn_patch                     (:)     ! patch allocation to fine root N (gN/m2/s)
     real(r8), pointer :: npool_to_frootn_storage_patch             (:)     ! patch allocation to fine root N storage (gN/m2/s)
     real(r8), pointer :: npool_to_livestemn_patch                  (:)     ! patch allocation to live stem N (gN/m2/s)
     real(r8), pointer :: npool_to_livestemn_storage_patch          (:)     ! patch allocation to live stem N storage (gN/m2/s)
     real(r8), pointer :: npool_to_deadstemn_patch                  (:)     ! patch allocation to dead stem N (gN/m2/s)
     real(r8), pointer :: npool_to_deadstemn_storage_patch          (:)     ! patch allocation to dead stem N storage (gN/m2/s)
     real(r8), pointer :: npool_to_livecrootn_patch                 (:)     ! patch allocation to live coarse root N (gN/m2/s)
     real(r8), pointer :: npool_to_livecrootn_storage_patch         (:)     ! patch allocation to live coarse root N storage (gN/m2/s)
     real(r8), pointer :: npool_to_deadcrootn_patch                 (:)     ! patch allocation to dead coarse root N (gN/m2/s)
     real(r8), pointer :: npool_to_deadcrootn_storage_patch         (:)     ! patch allocation to dead coarse root N storage (gN/m2/s)

     ! annual turnover of storage to transfer pools
     real(r8), pointer :: grainn_storage_to_xfer_patch              (:)     ! patch grain N shift storage to transfer for prognostic crop (gN/m2/s)
     real(r8), pointer :: leafn_storage_to_xfer_patch               (:)     ! patch leaf N shift storage to transfer (gN/m2/s)
     real(r8), pointer :: frootn_storage_to_xfer_patch              (:)     ! patch fine root N shift storage to transfer (gN/m2/s)
     real(r8), pointer :: livestemn_storage_to_xfer_patch           (:)     ! patch live stem N shift storage to transfer (gN/m2/s)
     real(r8), pointer :: deadstemn_storage_to_xfer_patch           (:)     ! patch dead stem N shift storage to transfer (gN/m2/s)
     real(r8), pointer :: livecrootn_storage_to_xfer_patch          (:)     ! patch live coarse root N shift storage to transfer (gN/m2/s)
     real(r8), pointer :: deadcrootn_storage_to_xfer_patch          (:)     ! patch dead coarse root N shift storage to transfer (gN/m2/s)
     real(r8), pointer :: fert_patch                                (:)     ! patch applied fertilizer (gN/m2/s)
     real(r8), pointer :: fert_counter_patch                        (:)     ! patch >0 fertilize; <=0 not
     real(r8), pointer :: soyfixn_patch                             (:)     ! patch soybean fixed N (gN/m2/s)

     ! turnover of livewood to deadwood, with retranslocation
     real(r8), pointer :: livestemn_to_deadstemn_patch              (:)     ! patch live stem N turnover (gN/m2/s)
     real(r8), pointer :: livestemn_to_retransn_patch               (:)     ! patch live stem N to retranslocated N pool (gN/m2/s)
     real(r8), pointer :: livecrootn_to_deadcrootn_patch            (:)     ! patch live coarse root N turnover (gN/m2/s)
     real(r8), pointer :: livecrootn_to_retransn_patch              (:)     ! patch live coarse root N to retranslocated N pool (gN/m2/s)

     ! summary (diagnostic) flux variables, not involved in mass balance
     real(r8), pointer :: ndeploy_patch                             (:)     ! patch total N deployed to growth and storage (gN/m2/s)
     real(r8), pointer :: ninputs_patch                             (:)     ! patch total N inputs to pft-level (gN/m2/s)
     real(r8), pointer :: noutputs_patch                            (:)     ! patch total N outputs from pft-level (gN/m2/s)
     real(r8), pointer :: wood_harvestn_patch                       (:)     ! patch total N losses to wood product pools (gN/m2/s)
     real(r8), pointer :: wood_harvestn_col                         (:)     ! col total N losses to wood product pools (gN/m2/s) (p2c)

     ! deposition fluxes
     real(r8), pointer :: ndep_to_sminn_col                         (:)     ! col atmospheric N deposition to soil mineral N (gN/m2/s)
     real(r8), pointer :: nfix_to_sminn_col                         (:)     ! col symbiotic/asymbiotic N fixation to soil mineral N (gN/m2/s)
     real(r8), pointer :: fert_to_sminn_col                         (:)     ! col fertilizer N to soil mineral N (gN/m2/s)
     real(r8), pointer :: soyfixn_to_sminn_col                      (:)     ! col soybean fixation to soil mineral N (gN/m2/s)

     ! phenology: litterfall and crop fluxes
     real(r8), pointer :: phenology_n_to_litr_met_n_col             (:,:)   ! col N fluxes associated with phenology (litterfall and crop) to litter metabolic pool (gN/m3/s)
     real(r8), pointer :: phenology_n_to_litr_cel_n_col             (:,:)   ! col N fluxes associated with phenology (litterfall and crop) to litter cellulose pool (gN/m3/s)
     real(r8), pointer :: phenology_n_to_litr_lig_n_col             (:,:)   ! col N fluxes associated with phenology (litterfall and crop) to litter lignin pool (gN/m3/s)

     ! gap mortality fluxes
     real(r8), pointer :: gap_mortality_n_to_litr_met_n_col         (:,:)   ! col N fluxes associated with gap mortality to litter metabolic pool (gN/m3/s)
     real(r8), pointer :: gap_mortality_n_to_litr_cel_n_col         (:,:)   ! col N fluxes associated with gap mortality to litter cellulose pool (gN/m3/s)
     real(r8), pointer :: gap_mortality_n_to_litr_lig_n_col         (:,:)   ! col N fluxes associated with gap mortality to litter lignin pool (gN/m3/s)
     real(r8), pointer :: gap_mortality_n_to_cwdn_col               (:,:)   ! col N fluxes associated with gap mortality to CWD pool (gN/m3/s)

     ! decomposition fluxes
     real(r8), pointer :: decomp_cascade_ntransfer_vr_col           (:,:,:) ! col vert-res transfer of N from donor to receiver pool along decomp. cascade (gN/m3/s)
     real(r8), pointer :: decomp_cascade_ntransfer_col              (:,:)   ! col vert-int (diagnostic) transfer of N from donor to receiver pool along decomp. cascade (gN/m2/s)
     real(r8), pointer :: decomp_cascade_sminn_flux_vr_col          (:,:,:) ! col vert-res mineral N flux for transition along decomposition cascade (gN/m3/s)
     real(r8), pointer :: decomp_cascade_sminn_flux_col             (:,:)   ! col vert-int (diagnostic) mineral N flux for transition along decomposition cascade (gN/m2/s)
                                                                            ! Used to update concentrations concurrently with vertical transport
     ! vertically-resolved immobilization fluxes
     real(r8), pointer :: potential_immob_vr_col                    (:,:)   ! col vertically-resolved potential N immobilization (gN/m3/s) at each level
     real(r8), pointer :: potential_immob_col                       (:)     ! col vert-int (diagnostic) potential N immobilization (gN/m2/s)
     real(r8), pointer :: actual_immob_vr_col                       (:,:)   ! col vertically-resolved actual N immobilization (gN/m3/s) at each level
     real(r8), pointer :: actual_immob_col                          (:)     ! col vert-int (diagnostic) actual N immobilization (gN/m2/s)
     real(r8), pointer :: sminn_to_plant_vr_col                     (:,:)   ! col vertically-resolved plant uptake of soil mineral N (gN/m3/s)
     real(r8), pointer :: sminn_to_plant_col                        (:)     ! col vert-int (diagnostic) plant uptake of soil mineral N (gN/m2/s)
     real(r8), pointer :: supplement_to_sminn_vr_col                (:,:)   ! col vertically-resolved supplemental N supply (gN/m3/s)
     real(r8), pointer :: supplement_to_sminn_col                   (:)     ! col vert-int (diagnostic) supplemental N supply (gN/m2/s)
     real(r8), pointer :: gross_nmin_vr_col                         (:,:)   ! col vertically-resolved gross rate of N mineralization (gN/m3/s)
     real(r8), pointer :: gross_nmin_col                            (:)     ! col vert-int (diagnostic) gross rate of N mineralization (gN/m2/s)
     real(r8), pointer :: net_nmin_vr_col                           (:,:)   ! col vertically-resolved net rate of N mineralization (gN/m3/s)
     real(r8), pointer :: net_nmin_col                              (:)     ! col vert-int (diagnostic) net rate of N mineralization (gN/m2/s)

     real(r8), pointer :: sminn_no3_input_vr_col                    (:,:)   !col no3 input, gN/m3/time step
     real(r8), pointer :: sminn_nh4_input_vr_col                    (:,:)   !col nh4 input, gN/m3/time step
     real(r8), pointer :: sminn_no3_input_col                       (:)     !col no3 input, gN/m2
     real(r8), pointer :: sminn_nh4_input_col                       (:)     !col nh4 input, gN/m2
     real(r8), pointer :: sminn_input_col                           (:)     !col minn input, gN/m2

     ! ---------- NITRIF_DENITRIF  ---------------------

     ! nitrification / denitrification fluxes
     real(r8), pointer :: f_nit_vr_col                              (:,:)   ! col (gN/m3/s) soil nitrification flux
     real(r8), pointer :: f_denit_vr_col                            (:,:)   ! col (gN/m3/s) soil denitrification flux
     real(r8), pointer :: f_nit_col                                 (:)     ! col (gN/m2/s) soil nitrification flux
     real(r8), pointer :: f_denit_col                               (:)     ! col (gN/m2/s) soil denitrification flux

     real(r8), pointer :: pot_f_nit_vr_col                          (:,:)   ! col (gN/m3/s) potential soil nitrification flux
     real(r8), pointer :: pot_f_denit_vr_col                        (:,:)   ! col (gN/m3/s) potential soil denitrification flux
     real(r8), pointer :: pot_f_nit_col                             (:)     ! col (gN/m2/s) potential soil nitrification flux
     real(r8), pointer :: pot_f_denit_col                           (:)     ! col (gN/m2/s) potential soil denitrification flux
     real(r8), pointer :: n2_n2o_ratio_denit_vr_col                 (:,:)   ! col ratio of N2 to N2O production by denitrification [gN/gN]
     real(r8), pointer :: f_n2o_denit_vr_col                        (:,:)   ! col flux of N2o from denitrification [gN/m^3/s]
     real(r8), pointer :: f_n2o_denit_col                           (:)     ! col flux of N2o from denitrification [gN/m^2/s]
     real(r8), pointer :: f_n2o_nit_vr_col                          (:,:)   ! col flux of N2o from nitrification [gN/m^3/s]
     real(r8), pointer :: f_n2o_nit_col                             (:)     ! col flux of N2o from nitrification [gN/m^2/s]

     ! immobilization / uptake fluxes
     real(r8), pointer :: actual_immob_no3_vr_col                   (:,:)   ! col vertically-resolved actual immobilization of NO3 (gN/m3/s)
     real(r8), pointer :: actual_immob_nh4_vr_col                   (:,:)   ! col vertically-resolved actual immobilization of NH4 (gN/m3/s)
     real(r8), pointer :: smin_no3_to_plant_vr_col                  (:,:)   ! col vertically-resolved plant uptake of soil NO3 (gN/m3/s)
     real(r8), pointer :: smin_nh4_to_plant_vr_col                  (:,:)   ! col vertically-resolved plant uptake of soil NH4 (gN/m3/s)
     real(r8), pointer :: actual_immob_no3_col                      (:)     ! col actual immobilization of NO3 (gN/m2/s)
     real(r8), pointer :: actual_immob_nh4_col                      (:)     ! col actual immobilization of NH4 (gN/m2/s)
     real(r8), pointer :: smin_no3_to_plant_col                     (:)     ! col plant uptake of soil NO3 (gN/m2/s)
     real(r8), pointer :: smin_nh4_to_plant_col                     (:)     ! col plant uptake of soil Nh4 (gN/m2/s)

     ! leaching fluxes
     real(r8), pointer :: smin_no3_leached_vr_col                   (:,:)   ! col vertically-resolved soil mineral NO3 loss to leaching (gN/m3/s)
     real(r8), pointer :: smin_no3_leached_col                      (:)     ! col soil mineral NO3 pool loss to leaching (gN/m2/s)
     real(r8), pointer :: smin_no3_runoff_vr_col                    (:,:)   ! col vertically-resolved rate of mineral NO3 loss with runoff (gN/m3/s)
     real(r8), pointer :: smin_no3_runoff_col                       (:)     ! col soil mineral NO3 pool loss to runoff (gN/m2/s)

     ! nitrification /denitrification diagnostic quantities
     real(r8), pointer :: smin_no3_massdens_vr_col                  (:,:)   ! col (ugN / g soil) soil nitrate concentration
     real(r8), pointer :: soil_bulkdensity_col                      (:,:)   ! col (kg soil / m3) bulk density of soil
     real(r8), pointer :: k_nitr_t_vr_col                           (:,:)
     real(r8), pointer :: k_nitr_ph_vr_col                          (:,:)
     real(r8), pointer :: k_nitr_h2o_vr_col                         (:,:)
     real(r8), pointer :: k_nitr_vr_col                             (:,:)
     real(r8), pointer :: wfps_vr_col                               (:,:)
     real(r8), pointer :: fmax_denit_carbonsubstrate_vr_col         (:,:)
     real(r8), pointer :: fmax_denit_nitrate_vr_col                 (:,:)
     real(r8), pointer :: f_denit_base_vr_col                       (:,:)   ! col nitrification and denitrification fluxes
     real(r8), pointer :: diffus_col                                (:,:)   ! col diffusivity (m2/s)
     real(r8), pointer :: ratio_k1_col                              (:,:)
     real(r8), pointer :: ratio_no3_co2_col                         (:,:)
     real(r8), pointer :: soil_co2_prod_col                         (:,:)
     real(r8), pointer :: fr_WFPS_col                               (:,:)

     real(r8), pointer :: r_psi_col                                 (:,:)
     real(r8), pointer :: anaerobic_frac_col                        (:,:)

     !----------- no NITRIF_DENITRIF--------------

     ! denitrification fluxes
     real(r8), pointer :: sminn_to_denit_decomp_cascade_vr_col      (:,:,:) ! col vertically-resolved denitrification along decomp cascade (gN/m3/s)
     real(r8), pointer :: sminn_to_denit_decomp_cascade_col         (:,:)   ! col vertically-integrated (diagnostic) denitrification along decomp cascade (gN/m2/s)
     real(r8), pointer :: sminn_to_denit_excess_vr_col              (:,:)   ! col vertically-resolved denitrification from excess mineral N pool (gN/m3/s)
     real(r8), pointer :: sminn_to_denit_excess_col                 (:)     ! col vertically-integrated (diagnostic) denitrification from excess mineral N pool (gN/m2/s)

     ! leaching fluxes
     real(r8), pointer :: sminn_leached_vr_col                      (:,:)   ! col vertically-resolved soil mineral N pool loss to leaching (gN/m3/s)
     real(r8), pointer :: sminn_leached_col                         (:)     ! col soil mineral N pool loss to leaching (gN/m2/s)

     ! dynamic landcover fluxes
     real(r8), pointer :: dwt_seedn_to_leaf_col                     (:)     ! col (gN/m2/s) seed source to PFT-level
     real(r8), pointer :: dwt_seedn_to_deadstem_col                 (:)     ! col (gN/m2/s) seed source to PFT-level
     real(r8), pointer :: dwt_conv_nflux_col                        (:)     ! col (gN/m2/s) conversion N flux (immediate loss to atm)
     real(r8), pointer :: dwt_prod10n_gain_col                      (:)     ! col (gN/m2/s) addition to 10-yr wood product pool
     real(r8), pointer :: dwt_prod100n_gain_col                     (:)     ! col (gN/m2/s) addition to 100-yr wood product pool
     real(r8), pointer :: dwt_frootn_to_litr_met_n_col              (:,:)   ! col (gN/m3/s) fine root to litter due to landcover change
     real(r8), pointer :: dwt_frootn_to_litr_cel_n_col              (:,:)   ! col (gN/m3/s) fine root to litter due to landcover change
     real(r8), pointer :: dwt_frootn_to_litr_lig_n_col              (:,:)   ! col (gN/m3/s) fine root to litter due to landcover change
     real(r8), pointer :: dwt_livecrootn_to_cwdn_col                (:,:)   ! col (gN/m3/s) live coarse root to CWD due to landcover change
     real(r8), pointer :: dwt_deadcrootn_to_cwdn_col                (:,:)   ! col (gN/m3/s) dead coarse root to CWD due to landcover change
     real(r8), pointer :: dwt_nloss_col                             (:)     ! col (gN/m2/s) total nitrogen loss from product pools and conversion

     ! wood product pool loss fluxes
     real(r8), pointer :: prod1n_loss_col                           (:)     ! col (gN/m2/s) decomposition loss from 1-yr crop product pool
     real(r8), pointer :: prod10n_loss_col                          (:)     ! col (gN/m2/s) decomposition loss from 10-yr wood product pool
     real(r8), pointer :: prod100n_loss_col                         (:)     ! col (gN/m2/s) decomposition loss from 100-yr wood product pool
     real(r8), pointer :: product_nloss_col                         (:)     ! col (gN/m2/s) total wood product nitrogen loss


     ! summary (diagnostic) flux variables, not involved in mass balance
     real(r8), pointer :: denit_col                                 (:)     ! col total rate of denitrification (gN/m2/s)
     real(r8), pointer :: ninputs_col                               (:)     ! col column-level N inputs (gN/m2/s)
     real(r8), pointer :: noutputs_col                              (:)     ! col column-level N outputs (gN/m2/s)
     real(r8), pointer :: som_n_leached_col                         (:)     ! col total SOM N loss from vertical transport (gN/m^2/s)
     real(r8), pointer :: decomp_npools_leached_col                 (:,:)   ! col N loss from vertical transport from each decomposing N pool (gN/m^2/s)
     real(r8), pointer :: decomp_npools_transport_tendency_col      (:,:,:) ! col N tendency due to vertical transport in decomposing N pools (gN/m^3/s)

     ! all n pools involved in decomposition
     real(r8), pointer :: decomp_npools_sourcesink_col              (:,:,:) ! col (gN/m3) change in decomposing n pools
                                                                            !     (sum of all additions and subtractions from stateupdate1).

     ! Misc
     real(r8), pointer :: plant_ndemand_patch                       (:)     ! N flux required to support initial GPP (gN/m2/s)
     real(r8), pointer :: avail_retransn_patch                      (:)     ! N flux available from retranslocation pool (gN/m2/s)
     real(r8), pointer :: plant_nalloc_patch                        (:)     ! total allocated N flux (gN/m2/s)

     ! bgc interfacepflotran
     !------------------------------------------------------------------------
     real(r8), pointer :: plant_ndemand_col                         (:)     ! col N flux required to support initial GPP (gN/m2/s)
     !! pflotran
     real(r8), pointer :: plant_ndemand_vr_col                      (:,:)   ! col vertically-resolved N flux required to support initial GPP (gN/m3/s)

     real(r8), pointer :: f_ngas_decomp_vr_col                      (:,:)   ! col vertically-resolved N emission from excess mineral N pool due to mineralization (gN/m3/s)
     real(r8), pointer :: f_ngas_decomp_col                         (:)     ! col N emission from excess mineral N pool due to mineralization (gN/m2/s)
     real(r8), pointer :: f_ngas_nitri_vr_col                       (:,:)   ! col vertically-resolved N emission from nitrification (gN/m3/s)
     real(r8), pointer :: f_ngas_nitri_col                          (:)     ! col vertically-resolved N emission from nitrification (gN/m2/s)
     real(r8), pointer :: f_ngas_denit_vr_col                       (:,:)   ! col vertically-resolved N emission from denitrification (gN/m3/s)
     real(r8), pointer :: f_ngas_denit_col                          (:)     ! col vertically-resolved N emission from denitrification (gN/m2/s)
    ! (from PF bgc disgassing-solving)
     real(r8), pointer :: f_n2o_soil_vr_col                         (:,:)   ! col flux of N2o from soil-N processes [gN/m^3/s]
     real(r8), pointer :: f_n2o_soil_col                            (:)     ! col flux of N2o from soil-N processes [gN/m^2/s]
     real(r8), pointer :: f_n2_soil_vr_col                          (:,:)   ! col flux of N2 from soil-N processes [gN/m^3/s]
     real(r8), pointer :: f_n2_soil_col                             (:)     ! col flux of N2 from soil-N processes [gN/m^2/s]

      ! for PF-bgc mass-balance error checking
     real(r8), pointer :: externaln_to_decomp_npools_col            (:,:,:) ! col net N fluxes associated with litter/som-adding/removal to decomp pools (gN/m3/s)
                                                                            ! (sum of all external N additions and removals, excluding decomposition/hr).
     real(r8), pointer :: externaln_to_decomp_delta_col             (:)     ! col summarized net N i/o changes associated with litter/som-adding/removal to decomp pools  btw time-step (gN/m2)
     real(r8), pointer :: no3_net_transport_vr_col                  (:,:)   ! col net NO3 transport associated with runoff/leaching (gN/m3/s)
     real(r8), pointer :: no3_net_transport_delta_col               (:)     ! col summarized net change of column-level N leaching to NO3 bwtn time-step (for balance checking) (gN/m2)
    !------------------------------------------------------------------------

     real(r8), pointer :: smin_no3_to_plant_patch                   (:)     ! pft level plant uptake of soil NO3 (gN/m2/s) BGC mode
     real(r8), pointer :: smin_nh4_to_plant_patch                   (:)     ! pft level plant uptake of soil Nh4 (gN/m2/s) BGC mode
     real(r8), pointer :: sminn_to_plant_patch                      (:)     ! pft level plant uptake of soil N (gN/m2/s) CN mode
     real(r8), pointer :: col_plant_ndemand_vr                      (:,:)   ! column-level plant N demand
     real(r8), pointer :: col_plant_nh4demand_vr                    (:,:)   ! column-level plant NH4 demand
     real(r8), pointer :: col_plant_no3demand_vr                    (:,:)   ! column-level plant NO3 demand
     real(r8), pointer :: col_plant_pdemand_vr                      (:,:)   ! column-level plant P demand
     real(r8), pointer :: plant_nh4demand_vr_patch                  (:,:)   ! pft-level plant NH4 demand BGC mode
     real(r8), pointer :: plant_no3demand_vr_patch                  (:,:)   ! pft-level plant NO3 demand BGC mode
     real(r8), pointer :: plant_ndemand_vr_patch                    (:,:)   ! pft-level plant N demand CN mode
     real(r8), pointer :: prev_leafn_to_litter_patch                (:)     ! previous timestep leaf N litterfall flux (gN/m2/s)
     real(r8), pointer :: prev_frootn_to_litter_patch               (:)     ! previous timestep froot N litterfall flux (gN/m2/s)
     real(r8), pointer :: pmnf_decomp_cascade                       (:,:,:) !potential mineral N flux, from one pool to another

     real(r8), pointer :: plant_n_uptake_flux                       (:)     ! for the purpose of mass balance check
     real(r8), pointer :: soil_n_immob_flux                         (:)     ! for the purpose of mass balance check
     real(r8), pointer :: soil_n_immob_flux_vr                      (:,:)   ! for the purpose of mass balance check
     real(r8), pointer :: soil_n_grossmin_flux                      (:)     ! for the purpose of mass balance check
     real(r8), pointer :: plant_to_litter_nflux                     (:)     ! for the purpose of mass balance check
     real(r8), pointer :: plant_to_cwd_nflux                        (:)     ! for the purpose of mass balance check

   contains

     procedure , public  :: Init
     procedure , public  :: Restart
     procedure , public  :: SetValues
     procedure , public  :: ZeroDWT
     procedure , public  :: Summary
     procedure , private :: InitAllocate
     procedure , private :: InitHistory
     procedure , private :: InitCold
     procedure , private :: Summary_betr
     procedure , private :: NSummary_interface

  end type nitrogenflux_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(nitrogenflux_type) :: this
    type(bounds_type), intent(in) :: bounds

    call this%InitAllocate (bounds)
    call this%InitHistory (bounds)
    call this%InitCold (bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize pft nitrogen flux
    !
    ! !ARGUMENTS:
    class (nitrogenflux_type) :: this
    type(bounds_type) , intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer           :: begp,endp
    integer           :: begc,endc
    !------------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc

    allocate(this%m_leafn_to_litter_patch                   (begp:endp)) ; this%m_leafn_to_litter_patch                   (:) = nan
    allocate(this%m_frootn_to_litter_patch                  (begp:endp)) ; this%m_frootn_to_litter_patch                  (:) = nan
    allocate(this%m_leafn_storage_to_litter_patch           (begp:endp)) ; this%m_leafn_storage_to_litter_patch           (:) = nan
    allocate(this%m_frootn_storage_to_litter_patch          (begp:endp)) ; this%m_frootn_storage_to_litter_patch          (:) = nan
    allocate(this%m_livestemn_storage_to_litter_patch       (begp:endp)) ; this%m_livestemn_storage_to_litter_patch       (:) = nan
    allocate(this%m_deadstemn_storage_to_litter_patch       (begp:endp)) ; this%m_deadstemn_storage_to_litter_patch       (:) = nan
    allocate(this%m_livecrootn_storage_to_litter_patch      (begp:endp)) ; this%m_livecrootn_storage_to_litter_patch      (:) = nan
    allocate(this%m_deadcrootn_storage_to_litter_patch      (begp:endp)) ; this%m_deadcrootn_storage_to_litter_patch      (:) = nan
    allocate(this%m_leafn_xfer_to_litter_patch              (begp:endp)) ; this%m_leafn_xfer_to_litter_patch              (:) = nan
    allocate(this%m_frootn_xfer_to_litter_patch             (begp:endp)) ; this%m_frootn_xfer_to_litter_patch             (:) = nan
    allocate(this%m_livestemn_xfer_to_litter_patch          (begp:endp)) ; this%m_livestemn_xfer_to_litter_patch          (:) = nan
    allocate(this%m_deadstemn_xfer_to_litter_patch          (begp:endp)) ; this%m_deadstemn_xfer_to_litter_patch          (:) = nan
    allocate(this%m_livecrootn_xfer_to_litter_patch         (begp:endp)) ; this%m_livecrootn_xfer_to_litter_patch         (:) = nan
    allocate(this%m_deadcrootn_xfer_to_litter_patch         (begp:endp)) ; this%m_deadcrootn_xfer_to_litter_patch         (:) = nan
    allocate(this%m_livestemn_to_litter_patch               (begp:endp)) ; this%m_livestemn_to_litter_patch               (:) = nan
    allocate(this%m_deadstemn_to_litter_patch               (begp:endp)) ; this%m_deadstemn_to_litter_patch               (:) = nan
    allocate(this%m_livecrootn_to_litter_patch              (begp:endp)) ; this%m_livecrootn_to_litter_patch              (:) = nan
    allocate(this%m_deadcrootn_to_litter_patch              (begp:endp)) ; this%m_deadcrootn_to_litter_patch              (:) = nan
    allocate(this%m_retransn_to_litter_patch                (begp:endp)) ; this%m_retransn_to_litter_patch                (:) = nan
    allocate(this%hrv_leafn_to_litter_patch                 (begp:endp)) ; this%hrv_leafn_to_litter_patch                 (:) = nan
    allocate(this%hrv_frootn_to_litter_patch                (begp:endp)) ; this%hrv_frootn_to_litter_patch                (:) = nan
    allocate(this%hrv_leafn_storage_to_litter_patch         (begp:endp)) ; this%hrv_leafn_storage_to_litter_patch         (:) = nan
    allocate(this%hrv_frootn_storage_to_litter_patch        (begp:endp)) ; this%hrv_frootn_storage_to_litter_patch        (:) = nan
    allocate(this%hrv_livestemn_storage_to_litter_patch     (begp:endp)) ; this%hrv_livestemn_storage_to_litter_patch     (:) = nan
    allocate(this%hrv_deadstemn_storage_to_litter_patch     (begp:endp)) ; this%hrv_deadstemn_storage_to_litter_patch     (:) = nan
    allocate(this%hrv_livecrootn_storage_to_litter_patch    (begp:endp)) ; this%hrv_livecrootn_storage_to_litter_patch    (:) = nan
    allocate(this%hrv_deadcrootn_storage_to_litter_patch    (begp:endp)) ; this%hrv_deadcrootn_storage_to_litter_patch    (:) = nan
    allocate(this%hrv_leafn_xfer_to_litter_patch            (begp:endp)) ; this%hrv_leafn_xfer_to_litter_patch            (:) = nan
    allocate(this%hrv_frootn_xfer_to_litter_patch           (begp:endp)) ; this%hrv_frootn_xfer_to_litter_patch           (:) = nan
    allocate(this%hrv_livestemn_xfer_to_litter_patch        (begp:endp)) ; this%hrv_livestemn_xfer_to_litter_patch        (:) = nan
    allocate(this%hrv_deadstemn_xfer_to_litter_patch        (begp:endp)) ; this%hrv_deadstemn_xfer_to_litter_patch        (:) = nan
    allocate(this%hrv_livecrootn_xfer_to_litter_patch       (begp:endp)) ; this%hrv_livecrootn_xfer_to_litter_patch       (:) = nan
    allocate(this%hrv_deadcrootn_xfer_to_litter_patch       (begp:endp)) ; this%hrv_deadcrootn_xfer_to_litter_patch       (:) = nan
    allocate(this%hrv_livestemn_to_litter_patch             (begp:endp)) ; this%hrv_livestemn_to_litter_patch             (:) = nan
    allocate(this%hrv_deadstemn_to_prod10n_patch            (begp:endp)) ; this%hrv_deadstemn_to_prod10n_patch            (:) = nan
    allocate(this%hrv_deadstemn_to_prod100n_patch           (begp:endp)) ; this%hrv_deadstemn_to_prod100n_patch           (:) = nan
    allocate(this%hrv_leafn_to_prod1n_patch                 (begp:endp)) ; this%hrv_leafn_to_prod1n_patch                 (:) = nan
    allocate(this%hrv_livestemn_to_prod1n_patch             (begp:endp)) ; this%hrv_livestemn_to_prod1n_patch             (:) = nan
    allocate(this%hrv_grainn_to_prod1n_patch                (begp:endp)) ; this%hrv_grainn_to_prod1n_patch                (:) = nan
    allocate(this%hrv_cropn_to_prod1n_patch                 (begp:endp)) ; this%hrv_cropn_to_prod1n_patch                 (:) = nan
    allocate(this%hrv_livecrootn_to_litter_patch            (begp:endp)) ; this%hrv_livecrootn_to_litter_patch            (:) = nan
    allocate(this%hrv_deadcrootn_to_litter_patch            (begp:endp)) ; this%hrv_deadcrootn_to_litter_patch            (:) = nan
    allocate(this%hrv_retransn_to_litter_patch              (begp:endp)) ; this%hrv_retransn_to_litter_patch              (:) = nan

    allocate(this%m_leafn_to_fire_patch                     (begp:endp)) ; this%m_leafn_to_fire_patch                     (:) = nan
    allocate(this%m_leafn_storage_to_fire_patch             (begp:endp)) ; this%m_leafn_storage_to_fire_patch             (:) = nan
    allocate(this%m_leafn_xfer_to_fire_patch                (begp:endp)) ; this%m_leafn_xfer_to_fire_patch                (:) = nan
    allocate(this%m_livestemn_to_fire_patch                 (begp:endp)) ; this%m_livestemn_to_fire_patch                 (:) = nan
    allocate(this%m_livestemn_storage_to_fire_patch         (begp:endp)) ; this%m_livestemn_storage_to_fire_patch         (:) = nan
    allocate(this%m_livestemn_xfer_to_fire_patch            (begp:endp)) ; this%m_livestemn_xfer_to_fire_patch            (:) = nan
    allocate(this%m_deadstemn_to_fire_patch                 (begp:endp)) ; this%m_deadstemn_to_fire_patch                 (:) = nan
    allocate(this%m_deadstemn_storage_to_fire_patch         (begp:endp)) ; this%m_deadstemn_storage_to_fire_patch         (:) = nan
    allocate(this%m_deadstemn_xfer_to_fire_patch            (begp:endp)) ; this%m_deadstemn_xfer_to_fire_patch            (:) = nan
    allocate(this%m_frootn_to_fire_patch                    (begp:endp)) ; this%m_frootn_to_fire_patch                    (:) = nan
    allocate(this%m_frootn_storage_to_fire_patch            (begp:endp)) ; this%m_frootn_storage_to_fire_patch            (:) = nan
    allocate(this%m_frootn_xfer_to_fire_patch               (begp:endp)) ; this%m_frootn_xfer_to_fire_patch               (:) = nan
    allocate(this%m_livecrootn_to_fire_patch                (begp:endp)) ;
    allocate(this%m_livecrootn_storage_to_fire_patch        (begp:endp)) ; this%m_livecrootn_storage_to_fire_patch        (:) = nan
    allocate(this%m_livecrootn_xfer_to_fire_patch           (begp:endp)) ; this%m_livecrootn_xfer_to_fire_patch           (:) = nan
    allocate(this%m_deadcrootn_to_fire_patch                (begp:endp)) ; this%m_deadcrootn_to_fire_patch                (:) = nan
    allocate(this%m_deadcrootn_storage_to_fire_patch        (begp:endp)) ; this%m_deadcrootn_storage_to_fire_patch        (:) = nan
    allocate(this%m_deadcrootn_xfer_to_fire_patch           (begp:endp)) ; this%m_deadcrootn_xfer_to_fire_patch           (:) = nan
    allocate(this%m_retransn_to_fire_patch                  (begp:endp)) ; this%m_retransn_to_fire_patch                  (:) = nan

    allocate(this%m_leafn_to_litter_fire_patch              (begp:endp)) ; this%m_leafn_to_litter_fire_patch              (:) = nan
    allocate(this%m_leafn_storage_to_litter_fire_patch      (begp:endp)) ; this%m_leafn_storage_to_litter_fire_patch      (:) = nan
    allocate(this%m_leafn_xfer_to_litter_fire_patch         (begp:endp)) ; this%m_leafn_xfer_to_litter_fire_patch         (:) = nan
    allocate(this%m_livestemn_to_litter_fire_patch          (begp:endp)) ; this%m_livestemn_to_litter_fire_patch          (:) = nan
    allocate(this%m_livestemn_storage_to_litter_fire_patch  (begp:endp)) ; this%m_livestemn_storage_to_litter_fire_patch  (:) = nan
    allocate(this%m_livestemn_xfer_to_litter_fire_patch     (begp:endp)) ; this%m_livestemn_xfer_to_litter_fire_patch     (:) = nan
    allocate(this%m_livestemn_to_deadstemn_fire_patch       (begp:endp)) ; this%m_livestemn_to_deadstemn_fire_patch       (:) = nan
    allocate(this%m_deadstemn_to_litter_fire_patch          (begp:endp)) ; this%m_deadstemn_to_litter_fire_patch          (:) = nan
    allocate(this%m_deadstemn_storage_to_litter_fire_patch  (begp:endp)) ; this%m_deadstemn_storage_to_litter_fire_patch  (:) = nan
    allocate(this%m_deadstemn_xfer_to_litter_fire_patch     (begp:endp)) ; this%m_deadstemn_xfer_to_litter_fire_patch     (:) = nan
    allocate(this%m_frootn_to_litter_fire_patch             (begp:endp)) ; this%m_frootn_to_litter_fire_patch             (:) = nan
    allocate(this%m_frootn_storage_to_litter_fire_patch     (begp:endp)) ; this%m_frootn_storage_to_litter_fire_patch     (:) = nan
    allocate(this%m_frootn_xfer_to_litter_fire_patch        (begp:endp)) ; this%m_frootn_xfer_to_litter_fire_patch        (:) = nan
    allocate(this%m_livecrootn_to_litter_fire_patch         (begp:endp)) ; this%m_livecrootn_to_litter_fire_patch         (:) = nan
    allocate(this%m_livecrootn_storage_to_litter_fire_patch (begp:endp)) ; this%m_livecrootn_storage_to_litter_fire_patch (:) = nan
    allocate(this%m_livecrootn_xfer_to_litter_fire_patch    (begp:endp)) ; this%m_livecrootn_xfer_to_litter_fire_patch    (:) = nan
    allocate(this%m_livecrootn_to_deadcrootn_fire_patch     (begp:endp)) ; this%m_livecrootn_to_deadcrootn_fire_patch     (:) = nan
    allocate(this%m_deadcrootn_to_litter_fire_patch         (begp:endp)) ; this%m_deadcrootn_to_litter_fire_patch         (:) = nan
    allocate(this%m_deadcrootn_storage_to_litter_fire_patch (begp:endp)) ; this%m_deadcrootn_storage_to_litter_fire_patch (:) = nan
    allocate(this%m_deadcrootn_xfer_to_litter_fire_patch    (begp:endp)) ; this%m_deadcrootn_xfer_to_litter_fire_patch    (:) = nan
    allocate(this%m_retransn_to_litter_fire_patch           (begp:endp)) ; this%m_retransn_to_litter_fire_patch           (:) = nan

    allocate(this%leafn_xfer_to_leafn_patch                 (begp:endp)) ; this%leafn_xfer_to_leafn_patch                 (:) = nan
    allocate(this%frootn_xfer_to_frootn_patch               (begp:endp)) ; this%frootn_xfer_to_frootn_patch               (:) = nan
    allocate(this%livestemn_xfer_to_livestemn_patch         (begp:endp)) ; this%livestemn_xfer_to_livestemn_patch         (:) = nan
    allocate(this%deadstemn_xfer_to_deadstemn_patch         (begp:endp)) ; this%deadstemn_xfer_to_deadstemn_patch         (:) = nan
    allocate(this%livecrootn_xfer_to_livecrootn_patch       (begp:endp)) ; this%livecrootn_xfer_to_livecrootn_patch       (:) = nan
    allocate(this%deadcrootn_xfer_to_deadcrootn_patch       (begp:endp)) ; this%deadcrootn_xfer_to_deadcrootn_patch       (:) = nan
    allocate(this%leafn_to_litter_patch                     (begp:endp)) ; this%leafn_to_litter_patch                     (:) = nan
    allocate(this%leafn_to_retransn_patch                   (begp:endp)) ; this%leafn_to_retransn_patch                   (:) = nan
    allocate(this%frootn_to_retransn_patch                  (begp:endp)) ; this%frootn_to_retransn_patch                  (:) = nan
    allocate(this%frootn_to_litter_patch                    (begp:endp)) ; this%frootn_to_litter_patch                    (:) = nan
    allocate(this%retransn_to_npool_patch                   (begp:endp)) ; this%retransn_to_npool_patch                   (:) = nan
    allocate(this%sminn_to_npool_patch                      (begp:endp)) ; this%sminn_to_npool_patch                      (:) = nan

    allocate(this%npool_to_leafn_patch              (begp:endp)) ; this%npool_to_leafn_patch              (:) = nan
    allocate(this%npool_to_leafn_storage_patch      (begp:endp)) ; this%npool_to_leafn_storage_patch      (:) = nan
    allocate(this%npool_to_frootn_patch             (begp:endp)) ; this%npool_to_frootn_patch             (:) = nan
    allocate(this%npool_to_frootn_storage_patch     (begp:endp)) ; this%npool_to_frootn_storage_patch     (:) = nan
    allocate(this%npool_to_livestemn_patch          (begp:endp)) ; this%npool_to_livestemn_patch          (:) = nan
    allocate(this%npool_to_livestemn_storage_patch  (begp:endp)) ; this%npool_to_livestemn_storage_patch  (:) = nan
    allocate(this%npool_to_deadstemn_patch          (begp:endp)) ; this%npool_to_deadstemn_patch          (:) = nan
    allocate(this%npool_to_deadstemn_storage_patch  (begp:endp)) ; this%npool_to_deadstemn_storage_patch  (:) = nan
    allocate(this%npool_to_livecrootn_patch         (begp:endp)) ; this%npool_to_livecrootn_patch         (:) = nan
    allocate(this%npool_to_livecrootn_storage_patch (begp:endp)) ; this%npool_to_livecrootn_storage_patch (:) = nan
    allocate(this%npool_to_deadcrootn_patch         (begp:endp)) ; this%npool_to_deadcrootn_patch         (:) = nan
    allocate(this%npool_to_deadcrootn_storage_patch (begp:endp)) ; this%npool_to_deadcrootn_storage_patch (:) = nan
    allocate(this%leafn_storage_to_xfer_patch       (begp:endp)) ; this%leafn_storage_to_xfer_patch       (:) = nan
    allocate(this%frootn_storage_to_xfer_patch      (begp:endp)) ; this%frootn_storage_to_xfer_patch      (:) = nan
    allocate(this%livestemn_storage_to_xfer_patch   (begp:endp)) ; this%livestemn_storage_to_xfer_patch   (:) = nan
    allocate(this%deadstemn_storage_to_xfer_patch   (begp:endp)) ; this%deadstemn_storage_to_xfer_patch   (:) = nan
    allocate(this%livecrootn_storage_to_xfer_patch  (begp:endp)) ; this%livecrootn_storage_to_xfer_patch  (:) = nan
    allocate(this%deadcrootn_storage_to_xfer_patch  (begp:endp)) ; this%deadcrootn_storage_to_xfer_patch  (:) = nan
    allocate(this%livestemn_to_deadstemn_patch      (begp:endp)) ; this%livestemn_to_deadstemn_patch      (:) = nan
    allocate(this%livestemn_to_retransn_patch       (begp:endp)) ; this%livestemn_to_retransn_patch       (:) = nan
    allocate(this%livecrootn_to_deadcrootn_patch    (begp:endp)) ; this%livecrootn_to_deadcrootn_patch    (:) = nan
    allocate(this%livecrootn_to_retransn_patch      (begp:endp)) ; this%livecrootn_to_retransn_patch      (:) = nan
    allocate(this%ndeploy_patch                     (begp:endp)) ; this%ndeploy_patch                     (:) = nan
    allocate(this%ninputs_patch                     (begp:endp)) ; this%ninputs_patch                     (:) = nan
    allocate(this%noutputs_patch                    (begp:endp)) ; this%noutputs_patch                    (:) = nan
    allocate(this%wood_harvestn_patch               (begp:endp)) ; this%wood_harvestn_patch               (:) = nan
    allocate(this%fire_nloss_patch                  (begp:endp)) ; this%fire_nloss_patch                  (:) = nan
    allocate(this%npool_to_grainn_patch             (begp:endp)) ; this%npool_to_grainn_patch             (:) = nan
    allocate(this%npool_to_grainn_storage_patch     (begp:endp)) ; this%npool_to_grainn_storage_patch     (:) = nan
    allocate(this%livestemn_to_litter_patch         (begp:endp)) ; this%livestemn_to_litter_patch         (:) = nan
    allocate(this%grainn_to_food_patch              (begp:endp)) ; this%grainn_to_food_patch              (:) = nan
    allocate(this%grainn_xfer_to_grainn_patch       (begp:endp)) ; this%grainn_xfer_to_grainn_patch       (:) = nan
    allocate(this%grainn_storage_to_xfer_patch      (begp:endp)) ; this%grainn_storage_to_xfer_patch      (:) = nan
    allocate(this%fert_patch                        (begp:endp)) ; this%fert_patch                        (:) = nan
    allocate(this%fert_counter_patch                (begp:endp)) ; this%fert_counter_patch                (:) = nan
    allocate(this%soyfixn_patch                     (begp:endp)) ; this%soyfixn_patch                     (:) = nan

    allocate(this%ndep_to_sminn_col             (begc:endc))    ; this%ndep_to_sminn_col	     (:) = nan
    allocate(this%nfix_to_sminn_col             (begc:endc))    ; this%nfix_to_sminn_col	     (:) = nan
    allocate(this%fert_to_sminn_col             (begc:endc))    ; this%fert_to_sminn_col	     (:) = nan
    allocate(this%soyfixn_to_sminn_col          (begc:endc))    ; this%soyfixn_to_sminn_col          (:) = nan
    allocate(this%hrv_deadstemn_to_prod10n_col  (begc:endc))    ; this%hrv_deadstemn_to_prod10n_col  (:) = nan
    allocate(this%hrv_deadstemn_to_prod100n_col (begc:endc))    ; this%hrv_deadstemn_to_prod100n_col (:) = nan
    allocate(this%hrv_cropn_to_prod1n_col       (begc:endc))    ; this%hrv_cropn_to_prod1n_col       (:) = nan
    allocate(this%sminn_to_plant_col            (begc:endc))    ; this%sminn_to_plant_col	     (:) = nan
    allocate(this%potential_immob_col           (begc:endc))    ; this%potential_immob_col           (:) = nan
    allocate(this%actual_immob_col              (begc:endc))    ; this%actual_immob_col              (:) = nan
    allocate(this%gross_nmin_col                (begc:endc))    ; this%gross_nmin_col                (:) = nan
    allocate(this%net_nmin_col                  (begc:endc))    ; this%net_nmin_col                  (:) = nan
    allocate(this%denit_col                     (begc:endc))    ; this%denit_col		     (:) = nan
    allocate(this%supplement_to_sminn_col       (begc:endc))    ; this%supplement_to_sminn_col       (:) = nan
    allocate(this%prod1n_loss_col               (begc:endc))    ; this%prod1n_loss_col               (:) = nan
    allocate(this%prod10n_loss_col              (begc:endc))    ; this%prod10n_loss_col              (:) = nan
    allocate(this%prod100n_loss_col             (begc:endc))    ; this%prod100n_loss_col	     (:) = nan
    allocate(this%product_nloss_col             (begc:endc))    ; this%product_nloss_col	     (:) = nan
    allocate(this%ninputs_col                   (begc:endc))    ; this%ninputs_col                   (:) = nan
    allocate(this%noutputs_col                  (begc:endc))    ; this%noutputs_col                  (:) = nan
    allocate(this%fire_nloss_col                (begc:endc))    ; this%fire_nloss_col                (:) = nan
    allocate(this%fire_decomp_nloss_col         (begc:endc))    ; this%fire_decomp_nloss_col         (:) = nan
    allocate(this%fire_nloss_p2c_col            (begc:endc))    ; this%fire_nloss_p2c_col            (:) = nan
    allocate(this%som_n_leached_col             (begc:endc))    ; this%som_n_leached_col	     (:) = nan

    allocate(this%m_n_to_litr_met_fire_col   (begc:endc,1:nlevdecomp_full)) ; this%m_n_to_litr_met_fire_col   (:,:) = nan
    allocate(this%m_n_to_litr_cel_fire_col   (begc:endc,1:nlevdecomp_full)) ; this%m_n_to_litr_cel_fire_col   (:,:) = nan
    allocate(this%m_n_to_litr_lig_fire_col   (begc:endc,1:nlevdecomp_full)) ; this%m_n_to_litr_lig_fire_col   (:,:) = nan
    allocate(this%r_psi_col                  (begc:endc,1:nlevdecomp_full)) ; this%r_psi_col                  (:,:) = spval
    allocate(this%anaerobic_frac_col         (begc:endc,1:nlevdecomp_full)) ; this%anaerobic_frac_col         (:,:) = spval
    allocate(this%potential_immob_vr_col     (begc:endc,1:nlevdecomp_full)) ; this%potential_immob_vr_col     (:,:) = nan
    allocate(this%actual_immob_vr_col        (begc:endc,1:nlevdecomp_full)) ; this%actual_immob_vr_col        (:,:) = nan
    allocate(this%sminn_to_plant_vr_col      (begc:endc,1:nlevdecomp_full)) ; this%sminn_to_plant_vr_col      (:,:) = nan
    allocate(this%supplement_to_sminn_vr_col (begc:endc,1:nlevdecomp_full)) ; this%supplement_to_sminn_vr_col (:,:) = nan
    allocate(this%gross_nmin_vr_col          (begc:endc,1:nlevdecomp_full)) ; this%gross_nmin_vr_col          (:,:) = nan
    allocate(this%net_nmin_vr_col            (begc:endc,1:nlevdecomp_full)) ; this%net_nmin_vr_col            (:,:) = nan

    allocate(this%dwt_seedn_to_leaf_col      (begc:endc))                   ; this%dwt_seedn_to_leaf_col      (:)   = nan
    allocate(this%dwt_seedn_to_deadstem_col  (begc:endc))                   ; this%dwt_seedn_to_deadstem_col  (:)   = nan
    allocate(this%dwt_conv_nflux_col         (begc:endc))                   ; this%dwt_conv_nflux_col         (:)   = nan
    allocate(this%dwt_prod10n_gain_col       (begc:endc))                   ; this%dwt_prod10n_gain_col       (:)   = nan
    allocate(this%dwt_prod100n_gain_col      (begc:endc))                   ; this%dwt_prod100n_gain_col      (:)   = nan
    allocate(this%dwt_nloss_col              (begc:endc))                   ; this%dwt_nloss_col              (:)   = nan
    allocate(this%wood_harvestn_col          (begc:endc))                   ; this%wood_harvestn_col          (:)   = nan

    allocate(this%dwt_frootn_to_litr_met_n_col(begc:endc,1:nlevdecomp_full)) ; this%dwt_frootn_to_litr_met_n_col     (:,:) = nan
    allocate(this%dwt_frootn_to_litr_cel_n_col(begc:endc,1:nlevdecomp_full)) ; this%dwt_frootn_to_litr_cel_n_col     (:,:) = nan
    allocate(this%dwt_frootn_to_litr_lig_n_col(begc:endc,1:nlevdecomp_full)) ; this%dwt_frootn_to_litr_lig_n_col     (:,:) = nan
    allocate(this%dwt_livecrootn_to_cwdn_col  (begc:endc,1:nlevdecomp_full)) ; this%dwt_livecrootn_to_cwdn_col       (:,:) = nan
    allocate(this%dwt_deadcrootn_to_cwdn_col  (begc:endc,1:nlevdecomp_full)) ; this%dwt_deadcrootn_to_cwdn_col       (:,:) = nan
    allocate(this%f_nit_vr_col                (begc:endc,1:nlevdecomp_full)) ; this%f_nit_vr_col                     (:,:) = nan
    allocate(this%f_denit_vr_col              (begc:endc,1:nlevdecomp_full)) ; this%f_denit_vr_col                   (:,:) = nan
    allocate(this%smin_no3_leached_vr_col     (begc:endc,1:nlevdecomp_full)) ; this%smin_no3_leached_vr_col          (:,:) = nan
    allocate(this%smin_no3_leached_col        (begc:endc))                   ; this%smin_no3_leached_col             (:)   = nan
    allocate(this%smin_no3_runoff_vr_col      (begc:endc,1:nlevdecomp_full)) ; this%smin_no3_runoff_vr_col           (:,:) = nan
    allocate(this%smin_no3_runoff_col         (begc:endc))                   ; this%smin_no3_runoff_col              (:)   = nan
    allocate(this%pot_f_nit_vr_col            (begc:endc,1:nlevdecomp_full)) ; this%pot_f_nit_vr_col                 (:,:) = nan
    allocate(this%pot_f_nit_col               (begc:endc))                   ; this%pot_f_nit_col                    (:)   = nan
    allocate(this%pot_f_denit_vr_col          (begc:endc,1:nlevdecomp_full)) ; this%pot_f_denit_vr_col               (:,:) = nan
    allocate(this%pot_f_denit_col             (begc:endc))                   ; this%pot_f_denit_col                  (:)   = nan
    allocate(this%actual_immob_no3_vr_col     (begc:endc,1:nlevdecomp_full)) ; this%actual_immob_no3_vr_col          (:,:) = nan
    allocate(this%actual_immob_nh4_vr_col     (begc:endc,1:nlevdecomp_full)) ; this%actual_immob_nh4_vr_col          (:,:) = nan
    allocate(this%smin_no3_to_plant_vr_col    (begc:endc,1:nlevdecomp_full)) ; this%smin_no3_to_plant_vr_col         (:,:) = nan
    allocate(this%smin_nh4_to_plant_vr_col    (begc:endc,1:nlevdecomp_full)) ; this%smin_nh4_to_plant_vr_col         (:,:) = nan
    allocate(this%smin_no3_to_plant_col       (begc:endc))                   ; this%smin_no3_to_plant_col            (:) = nan
    allocate(this%smin_nh4_to_plant_col       (begc:endc))                   ; this%smin_nh4_to_plant_col            (:) = nan

    allocate(this%f_nit_col                   (begc:endc))                   ; this%f_nit_col                        (:)   = nan
    allocate(this%f_denit_col                 (begc:endc))                   ; this%f_denit_col                      (:)   = nan
    allocate(this%n2_n2o_ratio_denit_vr_col   (begc:endc,1:nlevdecomp_full)) ; this%n2_n2o_ratio_denit_vr_col        (:,:) = nan
    allocate(this%f_n2o_denit_col             (begc:endc))                   ; this%f_n2o_denit_col                  (:)   = nan
    allocate(this%f_n2o_denit_vr_col          (begc:endc,1:nlevdecomp_full)) ; this%f_n2o_denit_vr_col               (:,:) = nan
    allocate(this%f_n2o_nit_col               (begc:endc))                   ; this%f_n2o_nit_col                    (:)   = nan
    allocate(this%f_n2o_nit_vr_col            (begc:endc,1:nlevdecomp_full)) ; this%f_n2o_nit_vr_col                 (:,:) = nan

    allocate(this%sminn_no3_input_vr_col      (begc:endc,1:nlevdecomp_full)) ; this%sminn_no3_input_vr_col           (:,:) = nan
    allocate(this%sminn_nh4_input_vr_col      (begc:endc,1:nlevdecomp_full)) ; this%sminn_nh4_input_vr_col           (:,:) = nan
    allocate(this%sminn_nh4_input_col         (begc:endc))                   ; this%sminn_nh4_input_col              (:)   = nan
    allocate(this%sminn_no3_input_col         (begc:endc))                   ; this%sminn_no3_input_col              (:)   = nan
    allocate(this%sminn_input_col             (begc:endc))                   ; this%sminn_input_col                  (:)   = nan

    allocate(this%smin_no3_massdens_vr_col    (begc:endc,1:nlevdecomp_full)) ; this%smin_no3_massdens_vr_col         (:,:) = nan
    allocate(this%soil_bulkdensity_col        (begc:endc,1:nlevdecomp_full)) ; this%soil_bulkdensity_col             (:,:) = nan
    allocate(this%k_nitr_t_vr_col             (begc:endc,1:nlevdecomp_full)) ; this%k_nitr_t_vr_col                  (:,:) = nan
    allocate(this%k_nitr_ph_vr_col            (begc:endc,1:nlevdecomp_full)) ; this%k_nitr_ph_vr_col                 (:,:) = nan
    allocate(this%k_nitr_h2o_vr_col           (begc:endc,1:nlevdecomp_full)) ; this%k_nitr_h2o_vr_col                (:,:) = nan
    allocate(this%k_nitr_vr_col               (begc:endc,1:nlevdecomp_full)) ; this%k_nitr_vr_col                    (:,:) = nan
    allocate(this%wfps_vr_col                 (begc:endc,1:nlevdecomp_full)) ; this%wfps_vr_col                      (:,:) = nan
    allocate(this%f_denit_base_vr_col         (begc:endc,1:nlevdecomp_full)) ; this%f_denit_base_vr_col              (:,:) = nan
    allocate(this%diffus_col                  (begc:endc,1:nlevdecomp_full)) ; this%diffus_col                       (:,:) = spval
    allocate(this%ratio_k1_col                (begc:endc,1:nlevdecomp_full)) ; this%ratio_k1_col                     (:,:) = nan
    allocate(this%ratio_no3_co2_col           (begc:endc,1:nlevdecomp_full)) ; this%ratio_no3_co2_col                (:,:) = spval
    allocate(this%soil_co2_prod_col           (begc:endc,1:nlevdecomp_full)) ; this%soil_co2_prod_col                (:,:) = nan
    allocate(this%fr_WFPS_col                 (begc:endc,1:nlevdecomp_full)) ; this%fr_WFPS_col                      (:,:) = spval

    allocate(this%fmax_denit_carbonsubstrate_vr_col(begc:endc,1:nlevdecomp_full)); this%fmax_denit_carbonsubstrate_vr_col(:,:) = nan
    allocate(this%fmax_denit_nitrate_vr_col        (begc:endc,1:nlevdecomp_full)); this%fmax_denit_nitrate_vr_col        (:,:) = nan

    allocate(this%decomp_cascade_ntransfer_vr_col   (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions ))
    allocate(this%decomp_cascade_sminn_flux_vr_col  (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions ))
    allocate(this%m_decomp_npools_to_fire_vr_col    (begc:endc,1:nlevdecomp_full,1:ndecomp_pools               ))
    allocate(this%decomp_cascade_ntransfer_col      (begc:endc,1:ndecomp_cascade_transitions                   ))
    allocate(this%decomp_cascade_sminn_flux_col     (begc:endc,1:ndecomp_cascade_transitions                   ))
    allocate(this%m_decomp_npools_to_fire_col       (begc:endc,1:ndecomp_pools                                 ))

    this%decomp_cascade_ntransfer_vr_col  (:,:,:) = nan
    this%decomp_cascade_sminn_flux_vr_col (:,:,:) = nan
    this%m_decomp_npools_to_fire_vr_col   (:,:,:) = nan
    this%m_decomp_npools_to_fire_col      (:,:)   = nan
    this%decomp_cascade_ntransfer_col     (:,:)   = nan
    this%decomp_cascade_sminn_flux_col    (:,:)   = nan

    allocate(this%phenology_n_to_litr_met_n_col     (begc:endc, 1:nlevdecomp_full))
    allocate(this%phenology_n_to_litr_cel_n_col     (begc:endc, 1:nlevdecomp_full))
    allocate(this%phenology_n_to_litr_lig_n_col     (begc:endc, 1:nlevdecomp_full))
    allocate(this%gap_mortality_n_to_litr_met_n_col (begc:endc, 1:nlevdecomp_full))
    allocate(this%gap_mortality_n_to_litr_cel_n_col (begc:endc, 1:nlevdecomp_full))
    allocate(this%gap_mortality_n_to_litr_lig_n_col (begc:endc, 1:nlevdecomp_full))
    allocate(this%gap_mortality_n_to_cwdn_col       (begc:endc, 1:nlevdecomp_full))
    allocate(this%fire_mortality_n_to_cwdn_col      (begc:endc, 1:nlevdecomp_full))
    allocate(this%harvest_n_to_litr_met_n_col       (begc:endc, 1:nlevdecomp_full))
    allocate(this%harvest_n_to_litr_cel_n_col       (begc:endc, 1:nlevdecomp_full))
    allocate(this%harvest_n_to_litr_lig_n_col       (begc:endc, 1:nlevdecomp_full))
    allocate(this%harvest_n_to_cwdn_col             (begc:endc, 1:nlevdecomp_full))

    this%phenology_n_to_litr_met_n_col     (:,:) = nan
    this%phenology_n_to_litr_cel_n_col     (:,:) = nan
    this%phenology_n_to_litr_lig_n_col     (:,:) = nan
    this%gap_mortality_n_to_litr_met_n_col (:,:) = nan
    this%gap_mortality_n_to_litr_cel_n_col (:,:) = nan
    this%gap_mortality_n_to_litr_lig_n_col (:,:) = nan
    this%gap_mortality_n_to_cwdn_col       (:,:) = nan
    this%fire_mortality_n_to_cwdn_col      (:,:) = nan
    this%harvest_n_to_litr_met_n_col       (:,:) = nan
    this%harvest_n_to_litr_cel_n_col       (:,:) = nan
    this%harvest_n_to_litr_lig_n_col       (:,:) = nan
    this%harvest_n_to_cwdn_col             (:,:) = nan

    allocate(this%sminn_to_denit_decomp_cascade_vr_col (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions ))
    allocate(this%sminn_to_denit_decomp_cascade_col    (begc:endc,1:ndecomp_cascade_transitions                   ))
    allocate(this%sminn_to_denit_excess_vr_col         (begc:endc,1:nlevdecomp_full                               ))
    allocate(this%sminn_to_denit_excess_col            (begc:endc                                                 ))
    allocate(this%sminn_leached_vr_col                 (begc:endc,1:nlevdecomp_full                               ))
    allocate(this%sminn_leached_col                    (begc:endc                                                 ))
    allocate(this%decomp_npools_leached_col            (begc:endc,1:ndecomp_pools                                 ))
    allocate(this%decomp_npools_transport_tendency_col (begc:endc,1:nlevdecomp_full,1:ndecomp_pools               ))

    this%sminn_to_denit_decomp_cascade_vr_col (:,:,:) = nan
    this%sminn_to_denit_decomp_cascade_col    (:,:)   = nan
    this%sminn_to_denit_excess_vr_col         (:,:)   = nan
    this%sminn_to_denit_excess_col            (:)     = nan
    this%sminn_leached_vr_col                 (:,:)   = nan
    this%sminn_leached_col                    (:)     = nan
    this%decomp_npools_leached_col            (:,:)   = nan
    this%decomp_npools_transport_tendency_col (:,:,:) = nan

    allocate(this%decomp_npools_sourcesink_col (begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
    this%decomp_npools_sourcesink_col (:,:,:) = nan

    allocate(this%plant_ndemand_patch         (begp:endp)) ;    this%plant_ndemand_patch         (:) = nan
    allocate(this%avail_retransn_patch        (begp:endp)) ;    this%avail_retransn_patch        (:) = nan
    allocate(this%plant_nalloc_patch          (begp:endp)) ;    this%plant_nalloc_patch          (:) = nan

    ! bgc interface & pflotran
    !------------------------------------------------------------------------
    allocate(this%plant_ndemand_col           (begc:endc))                   ; this%plant_ndemand_col                (:)    = nan
    allocate(this%plant_ndemand_vr_col        (begc:endc,1:nlevdecomp_full)) ; this%plant_ndemand_vr_col             (:,:)  = nan

    allocate(this%f_ngas_decomp_vr_col        (begc:endc,1:nlevdecomp_full)) ; this%f_ngas_decomp_vr_col             (:,:)  = nan
    allocate(this%f_ngas_nitri_vr_col         (begc:endc,1:nlevdecomp_full)) ; this%f_ngas_nitri_vr_col              (:,:)  = nan
    allocate(this%f_ngas_denit_vr_col         (begc:endc,1:nlevdecomp_full)) ; this%f_ngas_denit_vr_col              (:,:)  = nan
    allocate(this%f_n2o_soil_vr_col           (begc:endc,1:nlevdecomp_full)) ; this%f_n2o_soil_vr_col                (:,:)  = nan
    allocate(this%f_n2_soil_vr_col            (begc:endc,1:nlevdecomp_full)) ; this%f_n2_soil_vr_col                 (:,:)  = nan

    allocate(this%f_ngas_decomp_col           (begc:endc                  )) ; this%f_ngas_decomp_col                (:)    = nan
    allocate(this%f_ngas_nitri_col            (begc:endc                  )) ; this%f_ngas_nitri_col                 (:)    = nan
    allocate(this%f_ngas_denit_col            (begc:endc                  )) ; this%f_ngas_denit_col                 (:)    = nan
    allocate(this%f_n2o_soil_col              (begc:endc                  )) ; this%f_n2o_soil_col                   (:)    = nan
    allocate(this%f_n2_soil_col               (begc:endc                  )) ; this%f_n2_soil_col                    (:)    = nan

    allocate(this%externaln_to_decomp_npools_col    (begc:endc, 1:nlevdecomp_full, 1:ndecomp_pools)); this%externaln_to_decomp_npools_col    (:,:,:) = spval
    allocate(this%externaln_to_decomp_delta_col     (begc:endc))                                    ; this%externaln_to_decomp_delta_col     (:)     = spval
    allocate(this%no3_net_transport_vr_col          (begc:endc, 1:nlevdecomp_full))                 ; this%no3_net_transport_vr_col          (:,:)   = spval
    allocate(this%no3_net_transport_delta_col       (begc:endc))                                    ; this%no3_net_transport_delta_col       (:)     = spval
    !------------------------------------------------------------------------

    allocate(this%smin_no3_to_plant_patch     (begp:endp)) ;             this%smin_no3_to_plant_patch     (:) = nan
    allocate(this%smin_nh4_to_plant_patch     (begp:endp)) ;             this%smin_nh4_to_plant_patch     (:) = nan
    allocate(this%sminn_to_plant_patch        (begp:endp)) ;             this%sminn_to_plant_patch        (:) = nan
    allocate(this%col_plant_ndemand_vr        (begc:endc,1:nlevdecomp)); this%col_plant_ndemand_vr        (:,:) = nan
    allocate(this%col_plant_nh4demand_vr      (begc:endc,1:nlevdecomp)); this%col_plant_nh4demand_vr      (:,:) = nan
    allocate(this%col_plant_no3demand_vr      (begc:endc,1:nlevdecomp)); this%col_plant_no3demand_vr      (:,:) = nan
    allocate(this%col_plant_pdemand_vr        (begc:endc,1:nlevdecomp)); this%col_plant_pdemand_vr        (:,:) = nan
    allocate(this%plant_nh4demand_vr_patch    (begp:endp,1:nlevdecomp)); this%plant_nh4demand_vr_patch    (:,:) = nan
    allocate(this%plant_no3demand_vr_patch    (begp:endp,1:nlevdecomp)); this%plant_no3demand_vr_patch    (:,:) = nan
    allocate(this%plant_ndemand_vr_patch      (begp:endp,1:nlevdecomp)); this%plant_ndemand_vr_patch      (:,:) = nan
    allocate(this%prev_leafn_to_litter_patch  (begp:endp)) ;             this%prev_leafn_to_litter_patch  (:) = nan
    allocate(this%prev_frootn_to_litter_patch (begp:endp)) ;             this%prev_frootn_to_litter_patch (:) = nan
    allocate(this%pmnf_decomp_cascade(begc:endc,1:nlevdecomp,1:ndecomp_cascade_transitions)); this%pmnf_decomp_cascade(:,:,:) = nan

    allocate(this%plant_n_uptake_flux         (begc:endc)) ;	         this%plant_n_uptake_flux   (:)   = nan
    allocate(this%soil_n_immob_flux           (begc:endc)) ;	         this%soil_n_immob_flux	    (:)   = nan
    allocate(this%soil_n_immob_flux_vr        (begc:endc,1:nlevdecomp)); this%soil_n_immob_flux_vr  (:,:) = nan
    allocate(this%soil_n_grossmin_flux        (begc:endc)) ;	         this%soil_n_grossmin_flux  (:)   = nan
    allocate(this%actual_immob_no3_col        (begc:endc)) ;             this%actual_immob_no3_col  (:)   = nan
    allocate(this%actual_immob_nh4_col        (begc:endc)) ;             this%actual_immob_nh4_col  (:)   = nan
    allocate(this%smin_no3_to_plant_col       (begc:endc)) ;             this%smin_no3_to_plant_col (:)   = nan
    allocate(this%smin_nh4_to_plant_col       (begc:endc)) ;             this%smin_nh4_to_plant_col (:)   = nan
    allocate(this%plant_to_litter_nflux       (begc:endc)) ;             this%plant_to_litter_nflux (:)   = nan
    allocate(this%plant_to_cwd_nflux          (begc:endc)) ;             this%plant_to_cwd_nflux    (:)   = nan

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use clm_varpar     , only : nlevsno, nlevgrnd, crop_prog
    use histFileMod    , only : hist_addfld1d, hist_addfld2d, hist_addfld_decomp
    !
    ! !ARGUMENTS:
    class(nitrogenflux_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer        :: k,l
    integer        :: begp, endp
    integer        :: begc, endc
    character(10)  :: active
    character(24)  :: fieldname
    character(100) :: longname
    character(8)   :: vr_suffix
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    ! add suffix if number of soil decomposition depths is greater than 1
    if (nlevdecomp > 1) then
       vr_suffix = "_vr"
    else
       vr_suffix = ""
    endif

    this%m_leafn_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LEAFN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N mortality', &
         ptr_patch=this%m_leafn_to_litter_patch, default='inactive')

    this%m_frootn_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_FROOTN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N mortality', &
         ptr_patch=this%m_frootn_to_litter_patch, default='inactive')

    this%m_leafn_storage_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LEAFN_STORAGE_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N storage mortality', &
         ptr_patch=this%m_leafn_storage_to_litter_patch, default='inactive')

    this%m_frootn_storage_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_FROOTN_STORAGE_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N storage mortality', &
         ptr_patch=this%m_frootn_storage_to_litter_patch, default='inactive')

    this%m_livestemn_storage_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVESTEMN_STORAGE_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N storage mortality', &
         ptr_patch=this%m_livestemn_storage_to_litter_patch, default='inactive')

    this%m_deadstemn_storage_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMN_STORAGE_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N storage mortality', &
         ptr_patch=this%m_deadstemn_storage_to_litter_patch, default='inactive')

    this%m_livecrootn_storage_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVECROOTN_STORAGE_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N storage mortality', &
         ptr_patch=this%m_livecrootn_storage_to_litter_patch, default='inactive')

    this%m_deadcrootn_storage_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTN_STORAGE_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N storage mortality', &
         ptr_patch=this%m_deadcrootn_storage_to_litter_patch, default='inactive')

    this%m_leafn_xfer_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LEAFN_XFER_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N transfer mortality', &
         ptr_patch=this%m_leafn_xfer_to_litter_patch, default='inactive')

    this%m_frootn_xfer_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_FROOTN_XFER_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N transfer mortality', &
         ptr_patch=this%m_frootn_xfer_to_litter_patch, default='inactive')

    this%m_livestemn_xfer_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVESTEMN_XFER_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N transfer mortality', &
         ptr_patch=this%m_livestemn_xfer_to_litter_patch, default='inactive')

    this%m_deadstemn_xfer_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMN_XFER_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N transfer mortality', &
         ptr_patch=this%m_deadstemn_xfer_to_litter_patch, default='inactive')

    this%m_livecrootn_xfer_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVECROOTN_XFER_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N transfer mortality', &
         ptr_patch=this%m_livecrootn_xfer_to_litter_patch, default='inactive')

    this%m_deadcrootn_xfer_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTN_XFER_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N transfer mortality', &
         ptr_patch=this%m_deadcrootn_xfer_to_litter_patch, default='inactive')

    this%m_livestemn_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVESTEMN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N mortality', &
         ptr_patch=this%m_livestemn_to_litter_patch, default='inactive')

    this%m_deadstemn_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N mortality', &
         ptr_patch=this%m_deadstemn_to_litter_patch, default='inactive')

    this%m_livecrootn_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVECROOTN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N mortality', &
         ptr_patch=this%m_livecrootn_to_litter_patch, default='inactive')

    this%m_deadcrootn_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N mortality', &
         ptr_patch=this%m_deadcrootn_to_litter_patch, default='inactive')

    this%m_retransn_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_RETRANSN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='retranslocated N pool mortality', &
         ptr_patch=this%m_retransn_to_litter_patch, default='inactive')

    this%m_leafn_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LEAFN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N fire loss', &
         ptr_patch=this%m_leafn_to_fire_patch, default='inactive')

    this%m_frootn_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_FROOTN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N fire loss ', &
         ptr_patch=this%m_frootn_to_fire_patch, default='inactive')

    this%m_leafn_storage_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LEAFN_STORAGE_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N storage fire loss', &
         ptr_patch=this%m_leafn_storage_to_fire_patch, default='inactive')

    this%m_frootn_storage_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_FROOTN_STORAGE_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N storage fire loss', &
         ptr_patch=this%m_frootn_storage_to_fire_patch, default='inactive')

    this%m_livestemn_storage_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVESTEMN_STORAGE_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N storage fire loss', &
         ptr_patch=this%m_livestemn_storage_to_fire_patch, default='inactive')

    this%m_deadstemn_storage_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMN_STORAGE_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N storage fire loss', &
         ptr_patch=this%m_deadstemn_storage_to_fire_patch, default='inactive')

    this%m_livecrootn_storage_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVECROOTN_STORAGE_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N storage fire loss', &
         ptr_patch=this%m_livecrootn_storage_to_fire_patch, default='inactive')

    this%m_deadcrootn_storage_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTN_STORAGE_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N storage fire loss', &
         ptr_patch=this%m_deadcrootn_storage_to_fire_patch, default='inactive')

    this%m_leafn_xfer_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LEAFN_XFER_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N transfer fire loss', &
         ptr_patch=this%m_leafn_xfer_to_fire_patch, default='inactive')

    this%m_frootn_xfer_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_FROOTN_XFER_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N transfer fire loss', &
         ptr_patch=this%m_frootn_xfer_to_fire_patch, default='inactive')

    this%m_livestemn_xfer_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVESTEMN_XFER_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N transfer fire loss', &
         ptr_patch=this%m_livestemn_xfer_to_fire_patch, default='inactive')

    this%m_deadstemn_xfer_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMN_XFER_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N transfer fire loss', &
         ptr_patch=this%m_deadstemn_xfer_to_fire_patch, default='inactive')

    this%m_livecrootn_xfer_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVECROOTN_XFER_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N transfer fire loss', &
         ptr_patch=this%m_livecrootn_xfer_to_fire_patch, default='inactive')

    this%m_deadcrootn_xfer_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTN_XFER_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N transfer fire loss', &
         ptr_patch=this%m_deadcrootn_xfer_to_fire_patch, default='inactive')

    this%m_livestemn_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVESTEMN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N fire loss', &
         ptr_patch=this%m_livestemn_to_fire_patch, default='inactive')

    this%m_deadstemn_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N fire loss', &
         ptr_patch=this%m_deadstemn_to_fire_patch, default='inactive')

    this%m_deadstemn_to_litter_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMN_TO_LITTER_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N fire mortality to litter', &
         ptr_patch=this%m_deadstemn_to_litter_fire_patch, default='inactive')

    this%m_livecrootn_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVECROOTN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N fire loss', &
         ptr_patch=this%m_livecrootn_to_fire_patch, default='inactive')

    this%m_deadcrootn_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N fire loss', &
         ptr_patch=this%m_deadcrootn_to_fire_patch, default='inactive')

    this%m_deadcrootn_to_litter_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTN_TO_LITTER_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N fire mortality to litter', &
         ptr_patch=this%m_deadcrootn_to_litter_fire_patch, default='inactive')

    this%m_retransn_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_RETRANSN_TO_FIRE', units='gN/m^2/s', &
         avgflag='A', long_name='retranslocated N pool fire loss', &
         ptr_patch=this%m_retransn_to_fire_patch, default='inactive')

    this%leafn_xfer_to_leafn_patch(begp:endp) = spval
    call hist_addfld1d (fname='LEAFN_XFER_TO_LEAFN', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N growth from storage', &
         ptr_patch=this%leafn_xfer_to_leafn_patch, default='inactive')

    this%frootn_xfer_to_frootn_patch(begp:endp) = spval
    call hist_addfld1d (fname='FROOTN_XFER_TO_FROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N growth from storage', &
         ptr_patch=this%frootn_xfer_to_frootn_patch, default='inactive')

    this%livestemn_xfer_to_livestemn_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMN_XFER_TO_LIVESTEMN', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N growth from storage', &
         ptr_patch=this%livestemn_xfer_to_livestemn_patch, default='inactive')

    this%deadstemn_xfer_to_deadstemn_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADSTEMN_XFER_TO_DEADSTEMN', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N growth from storage', &
         ptr_patch=this%deadstemn_xfer_to_deadstemn_patch, default='inactive')

    this%livecrootn_xfer_to_livecrootn_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTN_XFER_TO_LIVECROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N growth from storage', &
         ptr_patch=this%livecrootn_xfer_to_livecrootn_patch, default='inactive')

    this%deadcrootn_xfer_to_deadcrootn_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADCROOTN_XFER_TO_DEADCROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N growth from storage', &
         ptr_patch=this%deadcrootn_xfer_to_deadcrootn_patch, default='inactive')

    this%leafn_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='LEAFN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N litterfall', &
         ptr_patch=this%leafn_to_litter_patch, default='inactive')

    this%leafn_to_retransn_patch(begp:endp) = spval
    call hist_addfld1d (fname='LEAFN_TO_RETRANSN', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N to retranslocated N pool', &
         ptr_patch=this%leafn_to_retransn_patch, default='inactive')

    this%frootn_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='FROOTN_TO_LITTER', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N litterfall', &
         ptr_patch=this%frootn_to_litter_patch, default='inactive')

    this%retransn_to_npool_patch(begp:endp) = spval
    call hist_addfld1d (fname='RETRANSN_TO_NPOOL', units='gN/m^2/s', &
         avgflag='A', long_name='deployment of retranslocated N', &
         ptr_patch=this%retransn_to_npool_patch)

    this%sminn_to_npool_patch(begp:endp) = spval
    call hist_addfld1d (fname='SMINN_TO_NPOOL', units='gN/m^2/s', &
         avgflag='A', long_name='deployment of soil mineral N uptake', &
         ptr_patch=this%sminn_to_npool_patch)

    this%npool_to_leafn_patch(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_LEAFN', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to leaf N', &
         ptr_patch=this%npool_to_leafn_patch, default='inactive')

    this%npool_to_leafn_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_LEAFN_STORAGE', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to leaf N storage', &
         ptr_patch=this%npool_to_leafn_storage_patch, default='inactive')

    this%npool_to_frootn_patch(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_FROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to fine root N', &
         ptr_patch=this%npool_to_frootn_patch, default='inactive')

    this%npool_to_frootn_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_FROOTN_STORAGE', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to fine root N storage', &
         ptr_patch=this%npool_to_frootn_storage_patch, default='inactive')

    this%npool_to_livestemn_patch(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_LIVESTEMN', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to live stem N', &
         ptr_patch=this%npool_to_livestemn_patch, default='inactive')

    this%npool_to_livestemn_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_LIVESTEMN_STORAGE', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to live stem N storage', &
         ptr_patch=this%npool_to_livestemn_storage_patch, default='inactive')

    this%npool_to_deadstemn_patch(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_DEADSTEMN', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to dead stem N', &
         ptr_patch=this%npool_to_deadstemn_patch, default='inactive')

    this%npool_to_deadstemn_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_DEADSTEMN_STORAGE', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to dead stem N storage', &
         ptr_patch=this%npool_to_deadstemn_storage_patch, default='inactive')

    this%npool_to_livecrootn_patch(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_LIVECROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to live coarse root N', &
         ptr_patch=this%npool_to_livecrootn_patch, default='inactive')

    this%npool_to_livecrootn_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_LIVECROOTN_STORAGE', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to live coarse root N storage', &
         ptr_patch=this%npool_to_livecrootn_storage_patch, default='inactive')

    this%npool_to_deadcrootn_patch(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_DEADCROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to dead coarse root N', &
         ptr_patch=this%npool_to_deadcrootn_patch, default='inactive')

    this%npool_to_deadcrootn_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='NPOOL_TO_DEADCROOTN_STORAGE', units='gN/m^2/s', &
         avgflag='A', long_name='allocation to dead coarse root N storage', &
         ptr_patch=this%npool_to_deadcrootn_storage_patch, default='inactive')

    this%leafn_storage_to_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='LEAFN_STORAGE_TO_XFER', units='gN/m^2/s', &
         avgflag='A', long_name='leaf N shift storage to transfer', &
         ptr_patch=this%leafn_storage_to_xfer_patch, default='inactive')

    this%frootn_storage_to_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='FROOTN_STORAGE_TO_XFER', units='gN/m^2/s', &
         avgflag='A', long_name='fine root N shift storage to transfer', &
         ptr_patch=this%frootn_storage_to_xfer_patch, default='inactive')

    this%livestemn_storage_to_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMN_STORAGE_TO_XFER', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N shift storage to transfer', &
         ptr_patch=this%livestemn_storage_to_xfer_patch, default='inactive')

    this%deadstemn_storage_to_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADSTEMN_STORAGE_TO_XFER', units='gN/m^2/s', &
         avgflag='A', long_name='dead stem N shift storage to transfer', &
         ptr_patch=this%deadstemn_storage_to_xfer_patch, default='inactive')

    this%livecrootn_storage_to_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTN_STORAGE_TO_XFER', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N shift storage to transfer', &
         ptr_patch=this%livecrootn_storage_to_xfer_patch, default='inactive')

    this%deadcrootn_storage_to_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADCROOTN_STORAGE_TO_XFER', units='gN/m^2/s', &
         avgflag='A', long_name='dead coarse root N shift storage to transfer', &
         ptr_patch=this%deadcrootn_storage_to_xfer_patch, default='inactive')

    this%livestemn_to_deadstemn_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMN_TO_DEADSTEMN', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N turnover', &
         ptr_patch=this%livestemn_to_deadstemn_patch, default='inactive')

    this%livestemn_to_retransn_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMN_TO_RETRANSN', units='gN/m^2/s', &
         avgflag='A', long_name='live stem N to retranslocated N pool', &
         ptr_patch=this%livestemn_to_retransn_patch, default='inactive')

    this%livecrootn_to_deadcrootn_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTN_TO_DEADCROOTN', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N turnover', &
         ptr_patch=this%livecrootn_to_deadcrootn_patch, default='inactive')

    this%livecrootn_to_retransn_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTN_TO_RETRANSN', units='gN/m^2/s', &
         avgflag='A', long_name='live coarse root N to retranslocated N pool', &
         ptr_patch=this%livecrootn_to_retransn_patch, default='inactive')

    this%ndeploy_patch(begp:endp) = spval
    call hist_addfld1d (fname='NDEPLOY', units='gN/m^2/s', &
         avgflag='A', long_name='total N deployed in new growth', &
         ptr_patch=this%ndeploy_patch)

    this%wood_harvestn_patch(begp:endp) = spval
    call hist_addfld1d (fname='WOOD_HARVESTN', units='gN/m^2/s', &
         avgflag='A', long_name='wood harvest N (to product pools)', &
         ptr_patch=this%wood_harvestn_patch)

    this%fire_nloss_patch(begp:endp) = spval
    call hist_addfld1d (fname='PFT_FIRE_NLOSS', units='gN/m^2/s', &
         avgflag='A', long_name='total pft-level fire N loss', &
         ptr_patch=this%fire_nloss_patch)

    if (crop_prog) then
       this%fert_patch(begp:endp) = spval
       call hist_addfld1d (fname='FERT', units='gN/m^2/s', &
            avgflag='A', long_name='fertilizer added', &
            ptr_patch=this%fert_patch)
    end if

    if (crop_prog) then
       this%soyfixn_patch(begp:endp) = spval
       call hist_addfld1d (fname='SOYFIXN', units='gN/m^2/s', &
            avgflag='A', long_name='soybean fixation', &
            ptr_patch=this%soyfixn_patch)
    end if

    if (crop_prog) then
       this%fert_counter_patch(begp:endp) = spval
       call hist_addfld1d (fname='FERT_COUNTER', units='seconds', &
            avgflag='A', long_name='time left to fertilize', &
            ptr_patch=this%fert_counter_patch)
    end if

    !-------------------------------
    ! N flux variables - native to column
    !-------------------------------

    this%ndep_to_sminn_col(begc:endc) = spval
    call hist_addfld1d (fname='NDEP_TO_SMINN', units='gN/m^2/s', &
         avgflag='A', long_name='atmospheric N deposition to soil mineral N', &
         ptr_col=this%ndep_to_sminn_col)

    this%nfix_to_sminn_col(begc:endc) = spval
    call hist_addfld1d (fname='NFIX_TO_SMINN', units='gN/m^2/s', &
         avgflag='A', long_name='symbiotic/asymbiotic N fixation to soil mineral N', &
         ptr_col=this%nfix_to_sminn_col)

    do k = 1, ndecomp_pools
       if ( decomp_cascade_con%is_litter(k) .or. decomp_cascade_con%is_cwd(k) ) then
          this%m_decomp_npools_to_fire_col(begc:endc,k) = spval
          data1dptr => this%m_decomp_npools_to_fire_col(:,k)
          fieldname = 'M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'N_TO_FIRE'
          longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' N fire loss'
          call hist_addfld1d (fname=fieldname, units='gN/m^2',  &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr, default='inactive')

          if ( nlevdecomp_full > 1 ) then
             this%m_decomp_npools_to_fire_vr_col(begc:endc,:,k) = spval
             data2dptr => this%m_decomp_npools_to_fire_vr_col(:,:,k)
             fieldname = 'M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'N_TO_FIRE'//trim(vr_suffix)
             longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' N fire loss'
             call hist_addfld_decomp (fname=fieldname, units='gN/m^3',  type2d='levdcmp', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data2dptr, default='inactive')
          endif
       endif
    end do

       do l = 1, ndecomp_cascade_transitions
          ! vertically integrated fluxes
          !-- mineralization/immobilization fluxes (none from CWD)
          if ( .not. decomp_cascade_con%is_cwd(decomp_cascade_con%cascade_donor_pool(l)) ) then
             this%decomp_cascade_sminn_flux_col(begc:endc,l) = spval
             data1dptr => this%decomp_cascade_sminn_flux_col(:,l)
             if ( decomp_cascade_con%cascade_receiver_pool(l) /= 0 ) then
                fieldname = 'SMINN_TO_'//&
                     trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))//'N_'//&
                     trim(decomp_cascade_con%decomp_pool_name_short(decomp_cascade_con%cascade_donor_pool(l)))
                longname =  'mineral N flux for decomp. of '&
                     //trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//&
                     'to '//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))
             else
                fieldname = trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))&
                     //'N_TO_SMINN'
                longname =  'mineral N flux for decomp. of '&
                     //trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))
             endif
             call hist_addfld1d (fname=fieldname, units='gN/m^2', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data1dptr)
          end if

          !-- transfer fluxes (none from terminal pool, if present)
          if ( decomp_cascade_con%cascade_receiver_pool(l) /= 0 ) then
             this%decomp_cascade_ntransfer_col(begc:endc,l) = spval
             data1dptr => this%decomp_cascade_ntransfer_col(:,l)
             fieldname = trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'N_TO_'//&
                  trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))//'N'
             longname =  'decomp. of '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))//&
                  ' N to '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_receiver_pool(l)))//' N'
             call hist_addfld1d (fname=fieldname, units='gN/m^2',  &
                  avgflag='A', long_name=longname, &
                  ptr_col=data1dptr)
          end if

          ! vertically resolved fluxes
          if ( nlevdecomp_full > 1 ) then
             !-- mineralization/immobilization fluxes (none from CWD)
             if ( .not. decomp_cascade_con%is_cwd(decomp_cascade_con%cascade_donor_pool(l)) ) then
                this%decomp_cascade_sminn_flux_vr_col(begc:endc,:,l) = spval
                data2dptr => this%decomp_cascade_sminn_flux_vr_col(:,:,l)
                if ( decomp_cascade_con%cascade_receiver_pool(l) /= 0 ) then
                   fieldname = 'SMINN_TO_'&
                        //trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))//'N_'//&
                        trim(decomp_cascade_con%decomp_pool_name_short(decomp_cascade_con%cascade_donor_pool(l)))//trim(vr_suffix)
                   longname =  'mineral N flux for decomp. of '&
                        //trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//&
                        'to '//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))
                else
                   fieldname = trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))&
                        //'N_TO_SMINN'//trim(vr_suffix)
                   longname =  'mineral N flux for decomp. of '&
                        //trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))
                endif
                call hist_addfld_decomp (fname=fieldname, units='gN/m^3',  type2d='levdcmp', &
                     avgflag='A', long_name=longname, &
                     ptr_col=data2dptr, default='inactive')
             endif

             !-- transfer fluxes (none from terminal pool, if present)
             if ( decomp_cascade_con%cascade_receiver_pool(l) /= 0 ) then
                this%decomp_cascade_ntransfer_vr_col(begc:endc,:,l) = spval
                data2dptr => this%decomp_cascade_ntransfer_vr_col(:,:,l)
                fieldname = trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'N_TO_'//&
                     trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))&
                     //'N'//trim(vr_suffix)
                longname =  'decomp. of '&
                     //trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))//&
                     ' N to '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_receiver_pool(l)))//' N'
                call hist_addfld_decomp (fname=fieldname, units='gN/m^3',  type2d='levdcmp', &
                     avgflag='A', long_name=longname, &
                     ptr_col=data2dptr, default='inactive')
             endif

          endif
       end do


    this%sminn_no3_input_vr_col(begc:endc,:) = spval
    data2dptr => this%sminn_no3_input_vr_col(:,:)
    fieldname='SMINN_NO3_INPUT_vr'
    call hist_addfld_decomp (fname=fieldname, units='gN/m^3/s',  type2d='levdcmp', &
        avgflag='A', long_name=longname, &
        ptr_col=data2dptr, default='inactive')

    this%sminn_nh4_input_vr_col(begc:endc,:)  = spval
    data2dptr => this%sminn_nh4_input_vr_col(:,:)
    fieldname='SMINN_NH4_INPUT_vr'
    call hist_addfld_decomp (fname=fieldname, units='gN/m^3/s',  type2d='levdcmp', &
        avgflag='A', long_name=longname, &
        ptr_col=data2dptr, default='inactive')

    this%denit_col(begc:endc) = spval
    call hist_addfld1d (fname='DENIT', units='gN/m^2/s', &
         avgflag='A', long_name='total rate of denitrification', &
         ptr_col=this%denit_col)

    this%som_n_leached_col(begc:endc) = spval
    call hist_addfld1d (fname='SOM_N_LEACHED', units='gN/m^2/s', &
         avgflag='A', long_name='total flux of N from SOM pools due to leaching', &
         ptr_col=this%som_n_leached_col, default='inactive')

    do k = 1, ndecomp_pools
       if ( .not. decomp_cascade_con%is_cwd(k) ) then
          this%decomp_npools_leached_col(begc:endc,k) = spval
          data1dptr => this%decomp_npools_leached_col(:,k)
          fieldname = 'M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'N_TO_LEACHING'
          longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' N leaching loss'
          call hist_addfld1d (fname=fieldname, units='gN/m^2/s', &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr, default='inactive')

          this%decomp_npools_transport_tendency_col(begc:endc,:,k) = spval
          data2dptr => this%decomp_npools_transport_tendency_col(:,:,k)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(k))//'N_TNDNCY_VERT_TRANSPORT'
          longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' N tendency due to vertical transport'
          call hist_addfld_decomp (fname=fieldname, units='gN/m^3/s',  type2d='levdcmp', &
               avgflag='A', long_name=longname, &
               ptr_col=data2dptr)
       end if
    end do


    if (.not. use_nitrif_denitrif) then
       do l = 1, ndecomp_cascade_transitions
          !-- denitrification fluxes (none from CWD)
          if ( .not. decomp_cascade_con%is_cwd(decomp_cascade_con%cascade_donor_pool(l)) ) then
             this%sminn_to_denit_decomp_cascade_col(begc:endc,l) = spval
             data1dptr => this%sminn_to_denit_decomp_cascade_col(:,l)
             fieldname = 'SMINN_TO_DENIT_'//trim(decomp_cascade_con%cascade_step_name(l))
             longname =  'denitrification for decomp. of '&
                  //trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))//&
                  'to '//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))
             call hist_addfld1d (fname=fieldname, units='gN/m^2',  &
                  avgflag='A', long_name=longname, &
                  ptr_col=data1dptr)
          endif

          if ( nlevdecomp_full > 1 ) then
             !-- denitrification fluxes (none from CWD)
             if ( .not. decomp_cascade_con%is_cwd(decomp_cascade_con%cascade_donor_pool(l)) ) then
                this%sminn_to_denit_decomp_cascade_vr_col(begc:endc,:,l) = spval
                data2dptr => this%sminn_to_denit_decomp_cascade_vr_col(:,:,l)
                fieldname = 'SMINN_TO_DENIT_'//trim(decomp_cascade_con%cascade_step_name(l))//trim(vr_suffix)
                longname =  'denitrification for decomp. of '&
                     //trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))//&
                     'to '//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))
                call hist_addfld_decomp (fname=fieldname, units='gN/m^3',  type2d='levdcmp', &
                     avgflag='A', long_name=longname, &
                     ptr_col=data2dptr, default='inactive')
             endif
          endif
       end do
    end if

    if (.not. use_nitrif_denitrif) then
       this%sminn_to_denit_excess_col(begc:endc) = spval
       call hist_addfld1d (fname='SMINN_TO_DENIT_EXCESS', units='gN/m^2/s',  &
            avgflag='A', long_name='denitrification from excess mineral N pool', &
            ptr_col=this%sminn_to_denit_excess_col, default='inactive')
    end if

    if (.not. use_nitrif_denitrif) then
       this%sminn_leached_col(begc:endc) = spval
       call hist_addfld1d (fname='SMINN_LEACHED', units='gN/m^2/s',   &
            avgflag='A', long_name='soil mineral N pool loss to leaching', &
            ptr_col=this%sminn_leached_col)
    end if

    if (.not. use_nitrif_denitrif) then
       if ( nlevdecomp_full > 1 ) then
          this%sminn_to_denit_excess_vr_col(begc:endc,:) = spval
          call hist_addfld_decomp (fname='SMINN_TO_DENIT_EXCESS'//trim(vr_suffix), units='gN/m^3/s',  type2d='levdcmp', &
               avgflag='A', long_name='denitrification from excess mineral N pool', &
               ptr_col=this%sminn_to_denit_excess_vr_col, default='inactive')

          this%sminn_leached_vr_col(begc:endc,:) = spval
          call hist_addfld_decomp (fname='SMINN_LEACHED'//trim(vr_suffix), units='gN/m^3/s',  type2d='levdcmp', &
               avgflag='A', long_name='soil mineral N pool loss to leaching', &
               ptr_col=this%sminn_leached_vr_col, default='inactive')
       endif
    end if

    if (use_nitrif_denitrif) then
       this%f_nit_col(begc:endc) = spval
       call hist_addfld1d (fname='F_NIT', units='gN/m^2/s',  &
            avgflag='A', long_name='nitrification flux', &
            ptr_col=this%f_nit_col)
    end if

    if (use_nitrif_denitrif) then
       this%f_denit_col(begc:endc) = spval
       call hist_addfld1d (fname='F_DENIT', units='gN/m^2/s', &
            avgflag='A', long_name='denitrification flux', &
            ptr_col=this%f_denit_col)
    end if

    if (use_nitrif_denitrif) then
       this%pot_f_nit_col(begc:endc) = spval
       call hist_addfld1d (fname='POT_F_NIT', units='gN/m^2/s', &
            avgflag='A', long_name='potential nitrification flux', &
            ptr_col=this%pot_f_nit_col)
    end if

    if (use_nitrif_denitrif) then
       this%pot_f_denit_col(begc:endc) = spval
       call hist_addfld1d (fname='POT_F_DENIT', units='gN/m^2/s', &
            avgflag='A', long_name='potential denitrification flux', &
            ptr_col=this%pot_f_denit_col)
    end if

    if (use_nitrif_denitrif) then
       this%smin_no3_leached_col(begc:endc) = spval
       call hist_addfld1d (fname='SMIN_NO3_LEACHED', units='gN/m^2/s', &
            avgflag='A', long_name='soil NO3 pool loss to leaching', &
            ptr_col=this%smin_no3_leached_col)
    end if

    if (use_nitrif_denitrif) then
       this%smin_no3_runoff_col(begc:endc) = spval
       call hist_addfld1d (fname='SMIN_NO3_RUNOFF', units='gN/m^2/s', &
            avgflag='A', long_name='soil NO3 pool loss to runoff', &
            ptr_col=this%smin_no3_runoff_col)
    end if

    if (use_nitrif_denitrif .and.  nlevdecomp_full > 1 ) then
       this%f_nit_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='F_NIT'//trim(vr_suffix), units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='nitrification flux', &
            ptr_col=this%f_nit_vr_col)
    end if

    if (use_nitrif_denitrif .and.  nlevdecomp_full > 1 ) then
       this%f_denit_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='F_DENIT'//trim(vr_suffix), units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='denitrification flux', &
            ptr_col=this%f_denit_vr_col)
    end if

    if (use_nitrif_denitrif .and.  nlevdecomp_full > 1 ) then
       this%pot_f_nit_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='POT_F_NIT'//trim(vr_suffix), units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='potential nitrification flux', &
            ptr_col=this%pot_f_nit_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif .and.  nlevdecomp_full > 1 ) then
       this%pot_f_denit_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='POT_F_DENIT'//trim(vr_suffix), units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='potential denitrification flux', &
            ptr_col=this%pot_f_denit_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif .and.  nlevdecomp_full > 1 ) then
       this%smin_no3_leached_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SMIN_NO3_LEACHED'//trim(vr_suffix), units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='soil NO3 pool loss to leaching', &
            ptr_col=this%smin_no3_leached_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif .and.  nlevdecomp_full > 1 ) then
       this%smin_no3_runoff_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SMIN_NO3_RUNOFF'//trim(vr_suffix), units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='soil NO3 pool loss to runoff', &
            ptr_col=this%smin_no3_runoff_vr_col, default='inactive')
    endif

    if (use_nitrif_denitrif) then
       this%n2_n2o_ratio_denit_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='n2_n2o_ratio_denit', units='gN/gN', type2d='levdcmp', &
            avgflag='A', long_name='n2_n2o_ratio_denit', &
            ptr_col=this%n2_n2o_ratio_denit_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%actual_immob_no3_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='ACTUAL_IMMOB_NO3', units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='immobilization of NO3', &
            ptr_col=this%actual_immob_no3_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%actual_immob_nh4_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='ACTUAL_IMMOB_NH4', units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='immobilization of NH4', &
            ptr_col=this%actual_immob_nh4_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%smin_no3_to_plant_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SMIN_NO3_TO_PLANT', units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='plant uptake of NO3', &
            ptr_col=this%smin_no3_to_plant_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%smin_nh4_to_plant_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SMIN_NH4_TO_PLANT', units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='plant uptake of NH4', &
            ptr_col=this%smin_nh4_to_plant_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%smin_no3_massdens_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SMIN_NO3_MASSDENS', units='ugN/cm^3 soil', type2d='levdcmp', &
            avgflag='A', long_name='SMIN_NO3_MASSDENS', &
            ptr_col=this%smin_no3_massdens_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%k_nitr_t_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='K_NITR_T', units='unitless', type2d='levdcmp', &
            avgflag='A', long_name='K_NITR_T', &
            ptr_col=this%k_nitr_t_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%k_nitr_ph_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='K_NITR_PH', units='unitless', type2d='levdcmp', &
            avgflag='A', long_name='K_NITR_PH', &
            ptr_col=this%k_nitr_ph_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%k_nitr_h2o_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='K_NITR_H2O', units='unitless', type2d='levdcmp', &
            avgflag='A', long_name='K_NITR_H2O', &
            ptr_col=this%k_nitr_h2o_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%k_nitr_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='K_NITR', units='1/s', type2d='levdcmp', &
            avgflag='A', long_name='K_NITR', &
            ptr_col=this%k_nitr_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%wfps_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='WFPS', units='percent', type2d='levdcmp', &
            avgflag='A', long_name='WFPS', &
            ptr_col=this%wfps_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%fmax_denit_carbonsubstrate_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='FMAX_DENIT_CARBONSUBSTRATE', units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='FMAX_DENIT_CARBONSUBSTRATE', &
            ptr_col=this%fmax_denit_carbonsubstrate_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%fmax_denit_nitrate_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='FMAX_DENIT_NITRATE', units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='FMAX_DENIT_NITRATE', &
            ptr_col=this%fmax_denit_nitrate_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%f_denit_base_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='F_DENIT_BASE', units='gN/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='F_DENIT_BASE', &
            ptr_col=this%f_denit_base_vr_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%diffus_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='diffus', units='m^2/s', type2d='levdcmp', &
            avgflag='A', long_name='diffusivity', &
            ptr_col=this%diffus_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%ratio_k1_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='ratio_k1', units='none', type2d='levdcmp', &
            avgflag='A', long_name='ratio_k1', &
            ptr_col=this%ratio_k1_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%ratio_no3_co2_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='ratio_no3_co2', units='ratio', type2d='levdcmp', &
            avgflag='A', long_name='ratio_no3_co2', &
            ptr_col=this%ratio_no3_co2_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%soil_co2_prod_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='soil_co2_prod', units='ug C / g soil / day', type2d='levdcmp', &
            avgflag='A', long_name='soil_co2_prod', &
            ptr_col=this%soil_co2_prod_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%fr_WFPS_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='fr_WFPS', units='fraction', type2d='levdcmp', &
            avgflag='A', long_name='fr_WFPS', &
            ptr_col=this%fr_WFPS_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%soil_bulkdensity_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='soil_bulkdensity', units='kg/m3', type2d='levdcmp', &
            avgflag='A', long_name='soil_bulkdensity', &
            ptr_col=this%soil_bulkdensity_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%anaerobic_frac_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='anaerobic_frac', units='m3/m3', type2d='levdcmp', &
            avgflag='A', long_name='anaerobic_frac', &
            ptr_col=this%anaerobic_frac_col, default='inactive')
    end if

    if (use_nitrif_denitrif) then
       this%r_psi_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='r_psi', units='m', type2d='levdcmp', &
            avgflag='A', long_name='r_psi', &
            ptr_col=this%r_psi_col, default='inactive')
    end if


    if ( use_nitrif_denitrif .and. nlevdecomp_full > 1 ) then
       this%potential_immob_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='POTENTIAL_IMMOB'//trim(vr_suffix), units='gN/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='potential N immobilization', &
            ptr_col=this%potential_immob_vr_col, default='inactive')
    end if

    if ( use_nitrif_denitrif .and. nlevdecomp_full > 1 ) then
       this%actual_immob_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='ACTUAL_IMMOB'//trim(vr_suffix), units='gN/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='actual N immobilization', &
            ptr_col=this%actual_immob_vr_col, default='inactive')
    end if

    if ( use_nitrif_denitrif .and. nlevdecomp_full > 1 ) then
       this%sminn_to_plant_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SMINN_TO_PLANT'//trim(vr_suffix), units='gN/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='plant uptake of soil mineral N', &
            ptr_col=this%sminn_to_plant_vr_col, default='inactive')
    end if


    if ( use_nitrif_denitrif .and. nlevdecomp_full > 1 ) then
       this%supplement_to_sminn_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SUPPLEMENT_TO_SMINN'//trim(vr_suffix), units='gN/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='supplemental N supply', &
            ptr_col=this%supplement_to_sminn_vr_col, default='inactive')
    end if

    if ( use_nitrif_denitrif .and. nlevdecomp_full > 1 ) then
       this%gross_nmin_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='GROSS_NMIN'//trim(vr_suffix), units='gN/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='gross rate of N mineralization', &
            ptr_col=this%gross_nmin_vr_col, default='inactive')
    end if

    if ( use_nitrif_denitrif .and. nlevdecomp_full > 1 ) then
       this%net_nmin_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='NET_NMIN'//trim(vr_suffix), units='gN/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='net rate of N mineralization', &
            ptr_col=this%net_nmin_vr_col, default='inactive')
    end if

    this%potential_immob_col(begc:endc) = spval
    call hist_addfld1d (fname='POTENTIAL_IMMOB', units='gN/m^2/s', &
         avgflag='A', long_name='potential N immobilization', &
         ptr_col=this%potential_immob_col)

    this%actual_immob_col(begc:endc) = spval
    call hist_addfld1d (fname='ACTUAL_IMMOB', units='gN/m^2/s', &
         avgflag='A', long_name='actual N immobilization', &
         ptr_col=this%actual_immob_col)

    this%sminn_to_plant_col(begc:endc) = spval
    call hist_addfld1d (fname='SMINN_TO_PLANT', units='gN/m^2/s', &
         avgflag='A', long_name='plant uptake of soil mineral N', &
         ptr_col=this%sminn_to_plant_col)

    this%supplement_to_sminn_col(begc:endc) = spval
    call hist_addfld1d (fname='SUPPLEMENT_TO_SMINN', units='gN/m^2/s', &
         avgflag='A', long_name='supplemental N supply', &
         ptr_col=this%supplement_to_sminn_col)

    this%gross_nmin_col(begc:endc) = spval
    call hist_addfld1d (fname='GROSS_NMIN', units='gN/m^2/s', &
         avgflag='A', long_name='gross rate of N mineralization', &
         ptr_col=this%gross_nmin_col)

    this%net_nmin_col(begc:endc) = spval
    call hist_addfld1d (fname='NET_NMIN', units='gN/m^2/s', &
         avgflag='A', long_name='net rate of N mineralization', &
         ptr_col=this%net_nmin_col)

    if (use_nitrif_denitrif) then
       this%f_n2o_nit_col(begc:endc) = spval
       call hist_addfld1d (fname='F_N2O_NIT', units='gN/m^2/s', &
            avgflag='A', long_name='nitrification N2O flux', &
            ptr_col=this%f_n2o_nit_col)

       this%f_n2o_denit_col(begc:endc) = spval
       call hist_addfld1d (fname='F_N2O_DENIT', units='gN/m^2/s', &
            avgflag='A', long_name='denitrification N2O flux', &
            ptr_col=this%f_n2o_denit_col)
    end if

    this%fire_nloss_col(begc:endc) = spval
    call hist_addfld1d (fname='COL_FIRE_NLOSS', units='gN/m^2/s', &
         avgflag='A', long_name='total column-level fire N loss', &
         ptr_col=this%fire_nloss_col, default='inactive')

    this%fire_decomp_nloss_col(begc:endc) = spval
    call hist_addfld1d (fname='DECOMP_FIRE_NLOSS', units='gN/m^2/s', &
         avgflag='A', long_name='fire N loss from decomposable pools', &
         ptr_col=this%fire_decomp_nloss_col, default='inactive')

    this%dwt_seedn_to_leaf_col(begc:endc) = spval
    call hist_addfld1d (fname='DWT_SEEDN_TO_LEAF', units='gN/m^2/s', &
         avgflag='A', long_name='seed source to PFT-level leaf', &
         ptr_col=this%dwt_seedn_to_leaf_col, default='inactive')

    this%dwt_seedn_to_deadstem_col(begc:endc) = spval
    call hist_addfld1d (fname='DWT_SEEDN_TO_DEADSTEM', units='gN/m^2/s', &
         avgflag='A', long_name='seed source to PFT-level deadstem', &
         ptr_col=this%dwt_seedn_to_deadstem_col, default='inactive')

    this%dwt_conv_nflux_col(begc:endc) = spval
    call hist_addfld1d (fname='DWT_CONV_NFLUX', units='gN/m^2/s', &
         avgflag='A', long_name='conversion N flux (immediate loss to atm)', &
         ptr_col=this%dwt_conv_nflux_col, default='inactive')

    this%dwt_prod10n_gain_col(begc:endc) = spval
    call hist_addfld1d (fname='DWT_PROD10N_GAIN', units='gN/m^2/s', &
         avgflag='A', long_name='addition to 10-yr wood product pool', &
         ptr_col=this%dwt_prod10n_gain_col, default='inactive')

    this%prod10n_loss_col(begc:endc) = spval
    call hist_addfld1d (fname='PROD10N_LOSS', units='gN/m^2/s', &
         avgflag='A', long_name='loss from 10-yr wood product pool', &
         ptr_col=this%prod10n_loss_col, default='inactive')

    this%dwt_prod100n_gain_col(begc:endc) = spval
    call hist_addfld1d (fname='DWT_PROD100N_GAIN', units='gN/m^2/s', &
         avgflag='A', long_name='addition to 100-yr wood product pool', &
         ptr_col=this%dwt_prod100n_gain_col, default='inactive')

    this%prod100n_loss_col(begc:endc) = spval
    call hist_addfld1d (fname='PROD100N_LOSS', units='gN/m^2/s', &
         avgflag='A', long_name='loss from 100-yr wood product pool', &
         ptr_col=this%prod100n_loss_col, default='inactive')

    this%prod1n_loss_col(begc:endc) = spval
    call hist_addfld1d (fname='PROD1N_LOSS', units='gN/m^2/s', &
         avgflag='A', long_name='loss from 1-yr crop product pool', &
         ptr_col=this%prod1n_loss_col, default='inactive')

    this%product_nloss_col(begc:endc) = spval
    call hist_addfld1d (fname='PRODUCT_NLOSS', units='gN/m^2/s', &
         avgflag='A', long_name='total N loss from wood product pools', &
         ptr_col=this%product_nloss_col, default='inactive')

    this%dwt_frootn_to_litr_met_n_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='DWT_FROOTN_TO_LITR_MET_N', units='gN/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='fine root to litter due to landcover change', &
         ptr_col=this%dwt_frootn_to_litr_met_n_col, default='inactive')

    this%dwt_frootn_to_litr_cel_n_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='DWT_FROOTN_TO_LITR_CEL_N', units='gN/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='fine root to litter due to landcover change', &
         ptr_col=this%dwt_frootn_to_litr_cel_n_col, default='inactive')

    this%dwt_frootn_to_litr_lig_n_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='DWT_FROOTN_TO_LITR_LIG_N', units='gN/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='fine root to litter due to landcover change', &
         ptr_col=this%dwt_frootn_to_litr_lig_n_col, default='inactive')

    this%dwt_livecrootn_to_cwdn_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='DWT_LIVECROOTN_TO_CWDN', units='gN/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='live coarse root to CWD due to landcover change', &
         ptr_col=this%dwt_livecrootn_to_cwdn_col, default='inactive')

    this%dwt_deadcrootn_to_cwdn_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='DWT_DEADCROOTN_TO_CWDN', units='gN/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='dead coarse root to CWD due to landcover change', &
         ptr_col=this%dwt_deadcrootn_to_cwdn_col, default='inactive')

    this%dwt_nloss_col(begc:endc) = spval
    call hist_addfld1d (fname='DWT_NLOSS', units='gN/m^2/s', &
         avgflag='A', long_name='total nitrogen loss from landcover conversion', &
         ptr_col=this%dwt_nloss_col, default='inactive')

    if (crop_prog) then
       this%fert_to_sminn_col(begc:endc) = spval
       call hist_addfld1d (fname='FERT_TO_SMINN', units='gN/m^2/s', &
            avgflag='A', long_name='fertilizer to soil mineral N', &
            ptr_col=this%fert_to_sminn_col)
    end if

    if (crop_prog) then
       this%soyfixn_to_sminn_col(begc:endc) = spval
       call hist_addfld1d (fname='SOYFIXN_TO_SMINN', units='gN/m^2/s', &
            avgflag='A', long_name='Soybean fixation to soil mineral N', &
            ptr_col=this%soyfixn_to_sminn_col)
    end if

    this%plant_ndemand_patch(begp:endp) = spval
    call hist_addfld1d (fname='PLANT_NDEMAND', units='gN/m^2/s', &
         avgflag='A', long_name='N flux required to support initial GPP', &
         ptr_patch=this%plant_ndemand_patch)

    this%avail_retransn_patch(begp:endp) = spval
    call hist_addfld1d (fname='AVAIL_RETRANSN', units='gN/m^2/s', &
         avgflag='A', long_name='N flux available from retranslocation pool', &
         ptr_patch=this%avail_retransn_patch, default='inactive')

    this%plant_nalloc_patch(begp:endp) = spval
    call hist_addfld1d (fname='PLANT_NALLOC', units='gN/m^2/s', &
         avgflag='A', long_name='total allocated N flux', &
         ptr_patch=this%plant_nalloc_patch, default='inactive')

    !!-----------------------------------------------------------
    !! bgc interface & pflotran
    this%plant_ndemand_col(begc:endc) = spval
       call hist_addfld1d (fname='PLANT_NDEMAND_COL', units='gN/m^2/s', &
            avgflag='A', long_name='N flux required to support initial GPP', &
            ptr_col=this%plant_ndemand_col)

    if (use_pflotran.and.pf_cmode) then
          this%f_ngas_decomp_col(begc:endc) = spval
          call hist_addfld1d (fname='F_NGAS_DECOMP', units='gN/m^2/s',  &
                avgflag='A', long_name='N gas emission from excess mineral N pool due to mineralization', &
                ptr_col=this%f_ngas_decomp_col, default='inactive')

          this%f_ngas_nitri_col(begc:endc) = spval
          call hist_addfld1d (fname='F_NGAS_NITRI', units='gN/m^2/s',  &
                avgflag='A', long_name='N gas emission from nitrification', &
                ptr_col=this%f_ngas_nitri_col, default='inactive')

          this%f_ngas_denit_col(begc:endc) = spval
          call hist_addfld1d (fname='F_NGAS_DENIT', units='gN/m^2/s',  &
                avgflag='A', long_name='N gas emission from denitrification', &
                ptr_col=this%f_ngas_denit_col, default='inactive')

          this%f_n2o_soil_col(begc:endc) = spval
          call hist_addfld1d (fname='F_N2O_SOIL', units='gN/m^2/s',  &
                avgflag='A', long_name='soil n2o exchange flux', &
                ptr_col=this%f_n2o_soil_col)

          this%f_n2_soil_col(begc:endc) = spval
          call hist_addfld1d (fname='F_N2_SOIL', units='gN/m^2/s',  &
                avgflag='A', long_name='soil n2 exchange flux', &
                ptr_col=this%f_n2_soil_col)

          this%smin_nh4_to_plant_col(begc:endc) = spval
          call hist_addfld1d (fname='SMIN_NH4_TO_PLANT', units='gN/m^2/s', &
               avgflag='A', long_name='plant uptake of NH4', &
               ptr_col=this%smin_nh4_to_plant_col, default='inactive')

          this%smin_no3_to_plant_col(begc:endc) = spval
          call hist_addfld1d (fname='SMIN_NO3_TO_PLANT', units='gN/m^2/s', &
               avgflag='A', long_name='plant uptake of NO3', &
               ptr_col=this%smin_no3_to_plant_col, default='inactive')
          !!---------------------------------------------------------------
          this%f_ngas_decomp_vr_col(begc:endc,:) = spval
            call hist_addfld_decomp (fname='F_NGAS_DECOMP'//trim(vr_suffix), units='gN/m^3/s',  type2d='levdcmp', &
               avgflag='A', long_name='n gas emission from excess mineral N pool due to mineralization', &
               ptr_col=this%f_ngas_decomp_vr_col, default='inactive')

            this%f_ngas_nitri_vr_col(begc:endc,:) = spval
            call hist_addfld_decomp (fname='F_NGAS_NITRI'//trim(vr_suffix), units='gN/m^3/s',  type2d='levdcmp', &
               avgflag='A', long_name='n gas emission in nitrification', &
               ptr_col=this%f_ngas_nitri_vr_col, default='inactive')

            this%f_ngas_denit_vr_col(begc:endc,:) = spval
            call hist_addfld_decomp (fname='F_NGAS_DENIT'//trim(vr_suffix), units='gN/m^3/s',  type2d='levdcmp', &
               avgflag='A', long_name='n gas emission in denitrification', &
               ptr_col=this%f_ngas_denit_vr_col, default='inactive')

            this%f_n2o_soil_vr_col(begc:endc,:) = spval
            call hist_addfld_decomp (fname='F_N2O_SOIL'//trim(vr_suffix), units='gN/m^3/s', type2d='levdcmp', &
               avgflag='A', long_name='soil N2O exchange flux', &
               ptr_col=this%f_n2o_soil_vr_col)

            this%f_n2_soil_vr_col(begc:endc,:) = spval
            call hist_addfld_decomp (fname='F_N2_SOIL'//trim(vr_suffix), units='gN/m^3/s', type2d='levdcmp', &
               avgflag='A', long_name='soil N2 exchange flux', &
               ptr_col=this%f_n2_soil_vr_col)

            this%plant_ndemand_vr_col(begc:endc,:) = spval
            call hist_addfld_decomp (fname='PLANT_NDEMAND'//trim(vr_suffix), units='gN/m^3/s', type2d='levdcmp', &
               avgflag='A', long_name='plant N demand distribution via roots', &
               ptr_col=this%plant_ndemand_vr_col, default='inactive')
    end if !! if (use_pflotran.and.pf_cmode)
    !!-----------------------------------------------------------
  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initializes time varying variables used only in coupled carbon-nitrogen mode (CN):
    !
    ! !USES:
    use clm_varpar      , only : crop_prog
    use landunit_varcon , only : istsoil, istcrop
    !
    ! !ARGUMENTS:
    class(nitrogenflux_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: p,c,l
    integer :: fp, fc                                    ! filter indices
    integer :: num_special_col                           ! number of good values in special_col filter
    integer :: num_special_patch                         ! number of good values in special_patch filter
    integer :: special_col(bounds%endc-bounds%begc+1)    ! special landunit filter - columns
    integer :: special_patch(bounds%endp-bounds%begp+1)  ! special landunit filter - patches
    !---------------------------------------------------------------------

    ! Set column filters

    num_special_col = 0
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%ifspecial(l)) then
          num_special_col = num_special_col + 1
          special_col(num_special_col) = c
       end if
    end do

    ! Set patch filters

    num_special_patch = 0
    do p = bounds%begp,bounds%endp
       l = pft%landunit(p)
       if (lun%ifspecial(l)) then
          num_special_patch = num_special_patch + 1
          special_patch(num_special_patch) = p
       end if
    end do

    !-----------------------------------------------
    ! initialize nitrogen flux variables
    !-----------------------------------------------

    do p = bounds%begp,bounds%endp
       l = pft%landunit(p)

       this%prev_leafn_to_litter_patch(p)  = 0._r8
       this%prev_frootn_to_litter_patch(p) = 0._r8

       if ( crop_prog )then
          this%fert_counter_patch(p)  = spval
          this%fert_patch(p)          = 0._r8
          this%soyfixn_patch(p)       = 0._r8
       end if

       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          this%fert_counter_patch(p)  = 0._r8
       end if

       if (lun%ifspecial(l)) then
          this%plant_ndemand_patch(p)  = spval
          this%avail_retransn_patch(p) = spval
          this%plant_nalloc_patch(p)   = spval
       end if
    end do

    ! initialize fields for special filters

    do fc = 1,num_special_col
       c = special_col(fc)
       this%dwt_nloss_col(c) = 0._r8
    end do

    call this%SetValues (&
         num_patch=num_special_patch, filter_patch=special_patch, value_patch=0._r8, &
         num_column=num_special_col, filter_column=special_col, value_column=0._r8)

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine Restart (this,  bounds, ncid, flag )
    !
    ! !DESCRIPTION:
    ! Read/write CN restart data for carbon state
    !
    ! !USES:
    use clm_varpar, only : crop_prog
    use restUtilMod
    use ncdio_pio
    ! pflotran
!    use clm_varctl, only : use_pflotran, pf_cmode, pf_hmode
    !
    ! !ARGUMENTS:
    class (nitrogenflux_type) :: this
    type(bounds_type) , intent(in)    :: bounds
    type(file_desc_t) , intent(inout) :: ncid   ! netcdf id
    character(len=*)  , intent(in)    :: flag   !'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer :: j,c ! indices
    logical :: readvar      ! determine if variable is on initial file
    real(r8), pointer :: ptr2d(:,:) ! temp. pointers for slicing larger arrays
    real(r8), pointer :: ptr1d(:)   ! temp. pointers for slicing larger arrays

    ! pflotran
    integer :: k
    character(len=128) :: varname   ! temporary
    !------------------------------------------------------------------------

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag, varname='fert_counter', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%fert_counter_patch)

       call restartvar(ncid=ncid, flag=flag, varname='fert', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%fert_patch)
    end if

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag,  varname='grainn_xfer_to_grainn', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='grain N growth from storage', units='gN/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%grainn_xfer_to_grainn_patch)
    end if

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag,  varname='livestemn_to_litter', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='livestem N to litter', units='gN/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemn_to_litter_patch)
    end if

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag,  varname='grainn_to_food', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='grain N to food', units='gN/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%grainn_to_food_patch)
    end if

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag,  varname='npool_to_grainn', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='allocation to grain N', units='gN/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%npool_to_grainn_patch)
    end if

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag,  varname='npool_to_grainn_storage', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='allocation to grain N storage', units='gN/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%npool_to_grainn_storage_patch)
    end if

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag, varname='grainn_storage_to_xfer', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='grain N shift storage to transfer', units='gN/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%grainn_storage_to_xfer_patch)
    end if

    if (use_nitrif_denitrif) then
       ! pot_f_nit_vr
       if (use_vertsoilc) then
          ptr2d => this%pot_f_nit_vr_col(:,:)
          call restartvar(ncid=ncid, flag=flag, varname='pot_f_nit_vr', xtype=ncd_double, &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='potential soil nitrification flux', units='gN/m3/s', &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
       else
          ptr1d => this%pot_f_nit_vr_col(:,1)
          call restartvar(ncid=ncid, flag=flag, varname='pot_f_nit_vr', xtype=ncd_double, &
               dim1name='column', &
               long_name='soil nitrification flux', units='gN/m3/s', &
               interpinic_flag='interp', readvar=readvar, data=ptr1d)
       end if
       if (flag=='read' .and. .not. readvar) then
          call endrun(msg= 'ERROR:: pot_f_nit_vr'//' is required on an initialization dataset' )
       end if
    end if

    if (use_nitrif_denitrif) then
       ! f_nit_vr
       if (use_vertsoilc) then
          ptr2d => this%f_nit_vr_col(:,:)
          call restartvar(ncid=ncid, flag=flag, varname='f_nit_vr', xtype=ncd_double, &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='soil nitrification flux', units='gN/m3/s', &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
       else
          ptr1d => this%f_nit_vr_col(:,1)
          call restartvar(ncid=ncid, flag=flag, varname='f_nit_vr', xtype=ncd_double, &
               dim1name='column', &
               long_name='soil nitrification flux', units='gN/m3/s', &
               interpinic_flag='interp', readvar=readvar, data=ptr1d)
       end if
       if (flag=='read' .and. .not. readvar) then
          call endrun(msg='ERROR:: f_nit_vr'//' is required on an initialization dataset'//&
               errMsg(__FILE__, __LINE__))
       end if
    end if

    call restartvar(ncid=ncid, flag=flag, varname='plant_ndemand', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%plant_ndemand_patch)

    call restartvar(ncid=ncid, flag=flag, varname='avail_retransn', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%avail_retransn_patch)

    call restartvar(ncid=ncid, flag=flag, varname='plant_nalloc', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%plant_nalloc_patch)

    ! pflotran
    !------------------------------------------------------------------------
    if (use_pflotran .and. pf_cmode) then
       ! externaln_to_decomp_npools_col
       do k = 1, ndecomp_pools
          varname=trim(decomp_cascade_con%decomp_pool_name_restart(k))//'external_n'
          if (use_vertsoilc) then
             ptr2d => this%externaln_to_decomp_npools_col(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_vr", xtype=ncd_double,  &
                  dim1name='column', dim2name='levgrnd', switchdim=.true., &
                  long_name='net organic N adding/removal/transport to soil', units='gN/m3/s', fill_value=spval, &
                  interpinic_flag='interp', readvar=readvar, data=ptr2d)
          else
             ptr1d => this%externaln_to_decomp_npools_col(:,1,k) ! nlevdecomp = 1; so treat as 1D variable
             call restartvar(ncid=ncid, flag=flag, varname=varname, xtype=ncd_double,  &
                  dim1name='column', &
                  long_name='net organic N adding/removal/transport to soil', units='gN/m3/s', fill_value=spval, &
                  interpinic_flag='interp' , readvar=readvar, data=ptr1d)
          end if
          if (flag=='read' .and. .not. readvar) then
          !   call endrun(msg='ERROR:: '//trim(varname)//' is required on an initialization dataset'//&
          !        errMsg(__FILE__, __LINE__))
             this%externaln_to_decomp_npools_col(:,:,k) = 0._r8
          end if
       end do

       !no3_net_transport_vr
       if (.not.pf_hmode) then
          if (use_vertsoilc) then
             ptr2d => this%no3_net_transport_vr_col(:,:)
             call restartvar(ncid=ncid, flag=flag, varname='no3_net_transport_vr', xtype=ncd_double, &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='net soil NO3-N transport', units='gN/m3/s', &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
          else
             ptr1d => this%no3_net_transport_vr_col(:,1)
             call restartvar(ncid=ncid, flag=flag, varname='no3_net_transport_vr', xtype=ncd_double, &
               dim1name='column', &
               long_name='net soil  NO3-N transport', units='gN/m3/s', &
               interpinic_flag='interp', readvar=readvar, data=ptr1d)
          end if
          if (flag=='read' .and. .not. readvar) then
          !   call endrun(msg='ERROR:: no3_net_transport_vr'//' is required on an initialization dataset'//&
          !     errMsg(__FILE__, __LINE__))
             this%no3_net_transport_vr_col(:,:) = 0._r8
          end if
       end if

    end if !! if (use_pflotran .and. pf_cmode)
    !------------------------------------------------------------------------

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine SetValues ( this, &
       num_patch, filter_patch, value_patch, &
       num_column, filter_column, value_column)
    !
    ! !DESCRIPTION:
    ! Set nitrogen flux variables
    !
    ! !ARGUMENTS:
    ! !ARGUMENTS:
    class (nitrogenflux_type) :: this
    integer , intent(in) :: num_patch
    integer , intent(in) :: filter_patch(:)
    real(r8), intent(in) :: value_patch
    integer , intent(in) :: num_column
    integer , intent(in) :: filter_column(:)
    real(r8), intent(in) :: value_column
    !
    ! !LOCAL VARIABLES:
    integer :: fi,i,j,k,l     ! loop index
    !------------------------------------------------------------------------

    do fi = 1,num_patch
       i=filter_patch(fi)

       this%m_leafn_to_litter_patch(i)                   = value_patch
       this%m_frootn_to_litter_patch(i)                  = value_patch
       this%m_leafn_storage_to_litter_patch(i)           = value_patch
       this%m_frootn_storage_to_litter_patch(i)          = value_patch
       this%m_livestemn_storage_to_litter_patch(i)       = value_patch
       this%m_deadstemn_storage_to_litter_patch(i)       = value_patch
       this%m_livecrootn_storage_to_litter_patch(i)      = value_patch
       this%m_deadcrootn_storage_to_litter_patch(i)      = value_patch
       this%m_leafn_xfer_to_litter_patch(i)              = value_patch
       this%m_frootn_xfer_to_litter_patch(i)             = value_patch
       this%m_livestemn_xfer_to_litter_patch(i)          = value_patch
       this%m_deadstemn_xfer_to_litter_patch(i)          = value_patch
       this%m_livecrootn_xfer_to_litter_patch(i)         = value_patch
       this%m_deadcrootn_xfer_to_litter_patch(i)         = value_patch
       this%m_livestemn_to_litter_patch(i)               = value_patch
       this%m_deadstemn_to_litter_patch(i)               = value_patch
       this%m_livecrootn_to_litter_patch(i)              = value_patch
       this%m_deadcrootn_to_litter_patch(i)              = value_patch
       this%m_retransn_to_litter_patch(i)                = value_patch
       this%hrv_leafn_to_litter_patch(i)                 = value_patch
       this%hrv_frootn_to_litter_patch(i)                = value_patch
       this%hrv_leafn_storage_to_litter_patch(i)         = value_patch
       this%hrv_frootn_storage_to_litter_patch(i)        = value_patch
       this%hrv_livestemn_storage_to_litter_patch(i)     = value_patch
       this%hrv_deadstemn_storage_to_litter_patch(i)     = value_patch
       this%hrv_livecrootn_storage_to_litter_patch(i)    = value_patch
       this%hrv_deadcrootn_storage_to_litter_patch(i)    = value_patch
       this%hrv_leafn_xfer_to_litter_patch(i)            = value_patch
       this%hrv_frootn_xfer_to_litter_patch(i)           = value_patch
       this%hrv_livestemn_xfer_to_litter_patch(i)        = value_patch
       this%hrv_deadstemn_xfer_to_litter_patch(i)        = value_patch
       this%hrv_livecrootn_xfer_to_litter_patch(i)       = value_patch
       this%hrv_deadcrootn_xfer_to_litter_patch(i)       = value_patch
       this%hrv_livestemn_to_litter_patch(i)             = value_patch
       this%hrv_deadstemn_to_prod10n_patch(i)            = value_patch
       this%hrv_deadstemn_to_prod100n_patch(i)           = value_patch
       this%hrv_leafn_to_prod1n_patch(i)                 = value_patch
       this%hrv_livestemn_to_prod1n_patch(i)             = value_patch
       this%hrv_grainn_to_prod1n_patch(i)                = value_patch
       this%hrv_cropn_to_prod1n_patch(i)                 = value_patch
       this%hrv_livecrootn_to_litter_patch(i)            = value_patch
       this%hrv_deadcrootn_to_litter_patch(i)            = value_patch
       this%hrv_retransn_to_litter_patch(i)              = value_patch

       this%m_leafn_to_fire_patch(i)                     = value_patch
       this%m_leafn_storage_to_fire_patch(i)             = value_patch
       this%m_leafn_xfer_to_fire_patch(i)                = value_patch
       this%m_livestemn_to_fire_patch(i)                 = value_patch
       this%m_livestemn_storage_to_fire_patch(i)         = value_patch
       this%m_livestemn_xfer_to_fire_patch(i)            = value_patch
       this%m_deadstemn_to_fire_patch(i)                 = value_patch
       this%m_deadstemn_storage_to_fire_patch(i)         = value_patch
       this%m_deadstemn_xfer_to_fire_patch(i)            = value_patch
       this%m_frootn_to_fire_patch(i)                    = value_patch
       this%m_frootn_storage_to_fire_patch(i)            = value_patch
       this%m_frootn_xfer_to_fire_patch(i)               = value_patch
       this%m_livecrootn_to_fire_patch(i)                = value_patch
       this%m_livecrootn_storage_to_fire_patch(i)        = value_patch
       this%m_livecrootn_xfer_to_fire_patch(i)           = value_patch
       this%m_deadcrootn_to_fire_patch(i)                = value_patch
       this%m_deadcrootn_storage_to_fire_patch(i)        = value_patch
       this%m_deadcrootn_xfer_to_fire_patch(i)           = value_patch
       this%m_retransn_to_fire_patch(i)                  = value_patch


       this%m_leafn_to_litter_fire_patch(i)              = value_patch
       this%m_leafn_storage_to_litter_fire_patch(i)      = value_patch
       this%m_leafn_xfer_to_litter_fire_patch(i)         = value_patch
       this%m_livestemn_to_litter_fire_patch(i)          = value_patch
       this%m_livestemn_storage_to_litter_fire_patch(i)  = value_patch
       this%m_livestemn_xfer_to_litter_fire_patch(i)     = value_patch
       this%m_livestemn_to_deadstemn_fire_patch(i)       = value_patch
       this%m_deadstemn_to_litter_fire_patch(i)          = value_patch
       this%m_deadstemn_storage_to_litter_fire_patch(i)  = value_patch
       this%m_deadstemn_xfer_to_litter_fire_patch(i)     = value_patch
       this%m_frootn_to_litter_fire_patch(i)             = value_patch
       this%m_frootn_storage_to_litter_fire_patch(i)     = value_patch
       this%m_frootn_xfer_to_litter_fire_patch(i)        = value_patch
       this%m_livecrootn_to_litter_fire_patch(i)         = value_patch
       this%m_livecrootn_storage_to_litter_fire_patch(i) = value_patch
       this%m_livecrootn_xfer_to_litter_fire_patch(i)    = value_patch
       this%m_livecrootn_to_deadcrootn_fire_patch(i)     = value_patch
       this%m_deadcrootn_to_litter_fire_patch(i)         = value_patch
       this%m_deadcrootn_storage_to_litter_fire_patch(i) = value_patch
       this%m_deadcrootn_xfer_to_litter_fire_patch(i)    = value_patch
       this%m_retransn_to_litter_fire_patch(i)           = value_patch

       this%leafn_xfer_to_leafn_patch(i)                 = value_patch
       this%frootn_xfer_to_frootn_patch(i)               = value_patch
       this%livestemn_xfer_to_livestemn_patch(i)         = value_patch
       this%deadstemn_xfer_to_deadstemn_patch(i)         = value_patch
       this%livecrootn_xfer_to_livecrootn_patch(i)       = value_patch
       this%deadcrootn_xfer_to_deadcrootn_patch(i)       = value_patch
       this%leafn_to_litter_patch(i)                     = value_patch
       this%leafn_to_retransn_patch(i)                   = value_patch
       this%frootn_to_litter_patch(i)                    = value_patch
       this%retransn_to_npool_patch(i)                   = value_patch
       this%sminn_to_npool_patch(i)                      = value_patch
       this%npool_to_leafn_patch(i)                      = value_patch
       this%npool_to_leafn_storage_patch(i)              = value_patch
       this%npool_to_frootn_patch(i)                     = value_patch
       this%npool_to_frootn_storage_patch(i)             = value_patch
       this%npool_to_livestemn_patch(i)                  = value_patch
       this%npool_to_livestemn_storage_patch(i)          = value_patch
       this%npool_to_deadstemn_patch(i)                  = value_patch
       this%npool_to_deadstemn_storage_patch(i)          = value_patch
       this%npool_to_livecrootn_patch(i)                 = value_patch
       this%npool_to_livecrootn_storage_patch(i)         = value_patch
       this%npool_to_deadcrootn_patch(i)                 = value_patch
       this%npool_to_deadcrootn_storage_patch(i)         = value_patch
       this%leafn_storage_to_xfer_patch(i)               = value_patch
       this%frootn_storage_to_xfer_patch(i)              = value_patch
       this%livestemn_storage_to_xfer_patch(i)           = value_patch
       this%deadstemn_storage_to_xfer_patch(i)           = value_patch
       this%livecrootn_storage_to_xfer_patch(i)          = value_patch
       this%deadcrootn_storage_to_xfer_patch(i)          = value_patch
       this%livestemn_to_deadstemn_patch(i)              = value_patch
       this%livestemn_to_retransn_patch(i)               = value_patch
       this%livecrootn_to_deadcrootn_patch(i)            = value_patch
       this%livecrootn_to_retransn_patch(i)              = value_patch
       this%ndeploy_patch(i)                             = value_patch
       this%ninputs_patch(i)                             = value_patch
       this%noutputs_patch(i)                            = value_patch
       this%wood_harvestn_patch(i)                       = value_patch
       this%fire_nloss_patch(i)                          = value_patch
    end do

    if ( crop_prog )then
       do fi = 1,num_patch
          i = filter_patch(fi)
          this%livestemn_to_litter_patch(i)              = value_patch
          this%grainn_to_food_patch(i)                   = value_patch
          this%grainn_xfer_to_grainn_patch(i)            = value_patch
          this%npool_to_grainn_patch(i)                  = value_patch
          this%npool_to_grainn_storage_patch(i)          = value_patch
          this%grainn_storage_to_xfer_patch(i)           = value_patch
          this%soyfixn_patch(i)                          = value_patch
          this%frootn_to_retransn_patch(i)               = value_patch
       end do
    end if

    do j = 1, nlevdecomp_full
       do fi = 1,num_column
          i = filter_column(fi)

          ! phenology: litterfall and crop fluxes associated wit
          this%phenology_n_to_litr_met_n_col(i,j)        = value_column
          this%phenology_n_to_litr_cel_n_col(i,j)        = value_column
          this%phenology_n_to_litr_lig_n_col(i,j)        = value_column

          ! gap mortality
          this%gap_mortality_n_to_litr_met_n_col(i,j)    = value_column
          this%gap_mortality_n_to_litr_cel_n_col(i,j)    = value_column
          this%gap_mortality_n_to_litr_lig_n_col(i,j)    = value_column
          this%gap_mortality_n_to_cwdn_col(i,j)          = value_column

          ! fire
          this%fire_mortality_n_to_cwdn_col(i,j)         = value_column
          this%m_n_to_litr_met_fire_col(i,j)             = value_column
          this%m_n_to_litr_cel_fire_col(i,j)             = value_column
          this%m_n_to_litr_lig_fire_col(i,j)             = value_column

          ! harvest
          this%harvest_n_to_litr_met_n_col(i,j)          = value_column
          this%harvest_n_to_litr_cel_n_col(i,j)          = value_column
          this%harvest_n_to_litr_lig_n_col(i,j)          = value_column
          this%harvest_n_to_cwdn_col(i,j)                = value_column

          if (.not. use_nitrif_denitrif) then
             this%sminn_to_denit_excess_vr_col(i,j)      = value_column
             this%sminn_leached_vr_col(i,j)              = value_column
          else
             this%f_nit_vr_col(i,j)                      = value_column
             this%f_denit_vr_col(i,j)                    = value_column
             this%smin_no3_leached_vr_col(i,j)           = value_column
             this%smin_no3_runoff_vr_col(i,j)            = value_column
             this%n2_n2o_ratio_denit_vr_col(i,j)         = value_column
             this%pot_f_nit_vr_col(i,j)                  = value_column
             this%pot_f_denit_vr_col(i,j)                = value_column
             this%actual_immob_no3_vr_col(i,j)           = value_column
             this%actual_immob_nh4_vr_col(i,j)           = value_column
             this%smin_no3_to_plant_vr_col(i,j)          = value_column
             this%smin_nh4_to_plant_vr_col(i,j)          = value_column
             this%f_n2o_denit_vr_col(i,j)                = value_column
             this%f_n2o_nit_vr_col(i,j)                  = value_column

             this%smin_no3_massdens_vr_col(i,j)          = value_column
             this%k_nitr_t_vr_col(i,j)                   = value_column
             this%k_nitr_ph_vr_col(i,j)                  = value_column
             this%k_nitr_h2o_vr_col(i,j)                 = value_column
             this%k_nitr_vr_col(i,j)                     = value_column
             this%wfps_vr_col(i,j)                       = value_column
             this%fmax_denit_carbonsubstrate_vr_col(i,j) = value_column
             this%fmax_denit_nitrate_vr_col(i,j)         = value_column
             this%f_denit_base_vr_col(i,j)               = value_column

             this%diffus_col(i,j)                        = value_column
             this%ratio_k1_col(i,j)                      = value_column
             this%ratio_no3_co2_col(i,j)                 = value_column
             this%soil_co2_prod_col(i,j)                 = value_column
             this%fr_WFPS_col(i,j)                       = value_column
             this%soil_bulkdensity_col(i,j)              = value_column

             this%r_psi_col(i,j)                         = value_column
             this%anaerobic_frac_col(i,j)                = value_column

             ! pflotran
             this%plant_ndemand_vr_col(i,j)              = value_column
             this%f_ngas_decomp_vr_col(i,j)              = value_column
             this%f_ngas_nitri_vr_col(i,j)               = value_column
             this%f_ngas_denit_vr_col(i,j)               = value_column
             this%f_n2o_soil_vr_col(i,j)                 = value_column
             this%f_n2_soil_vr_col(i,j)                  = value_column
          end if
          this%potential_immob_vr_col(i,j)               = value_column
          this%actual_immob_vr_col(i,j)                  = value_column
          this%sminn_to_plant_vr_col(i,j)                = value_column
          this%supplement_to_sminn_vr_col(i,j)           = value_column
          this%gross_nmin_vr_col(i,j)                    = value_column
          this%net_nmin_vr_col(i,j)                      = value_column
          this%sminn_nh4_input_vr_col(i,j)               = value_column
          this%sminn_no3_input_vr_col(i,j)               = value_column
       end do
    end do

    do fi = 1,num_column
       i = filter_column(fi)

       this%ndep_to_sminn_col(i)             = value_column
       this%nfix_to_sminn_col(i)             = value_column
       this%fert_to_sminn_col(i)             = value_column
       this%soyfixn_to_sminn_col(i)          = value_column
       this%hrv_deadstemn_to_prod10n_col(i)  = value_column
       this%hrv_deadstemn_to_prod100n_col(i) = value_column
       this%hrv_cropn_to_prod1n_col(i)       = value_column
       this%prod10n_loss_col(i)              = value_column
       this%prod100n_loss_col(i)             = value_column
       this%prod1n_loss_col(i)               = value_column
       this%product_nloss_col(i)             = value_column
       this%potential_immob_col(i)           = value_column
       this%actual_immob_col(i)              = value_column
       this%sminn_to_plant_col(i)            = value_column
       this%supplement_to_sminn_col(i)       = value_column
       this%gross_nmin_col(i)                = value_column
       this%net_nmin_col(i)                  = value_column
       this%denit_col(i)                     = value_column
       if (use_nitrif_denitrif) then
          this%f_nit_col(i)                  = value_column
          this%pot_f_nit_col(i)              = value_column
          this%f_denit_col(i)                = value_column
          this%pot_f_denit_col(i)            = value_column
          this%f_n2o_denit_col(i)            = value_column
          this%f_n2o_nit_col(i)              = value_column
          this%smin_no3_leached_col(i)       = value_column
          this%smin_no3_runoff_col(i)        = value_column

          this%f_ngas_decomp_col(i)         = value_column
          this%f_ngas_nitri_col(i)          = value_column
          this%f_ngas_denit_col(i)          = value_column
          this%f_n2o_soil_col(i)            = value_column
          this%f_n2_soil_col(i)             = value_column

          this%smin_nh4_to_plant_col(i)      = value_column
          this%smin_no3_to_plant_col(i)      = value_column
       else
          this%sminn_to_denit_excess_col(i)  = value_column
          this%sminn_leached_col(i)          = value_column
       end if
       this%ninputs_col(i)                   = value_column
       this%noutputs_col(i)                  = value_column
       this%fire_nloss_col(i)                = value_column
       this%som_n_leached_col(i)             = value_column
       this%sminn_input_col(i)               = value_column
       this%sminn_nh4_input_col(i)           = value_column
       this%sminn_no3_input_col(i)           = value_column
       ! Zero p2c column fluxes
       this%fire_nloss_col(i) = value_column
       this%wood_harvestn_col(i) = value_column

       !! bgc-interface
       this%plant_ndemand_col(i) = value_column
    end do

    do k = 1, ndecomp_pools
       do fi = 1,num_column
          i = filter_column(fi)
          this%decomp_npools_leached_col(i,k) = value_column
          this%m_decomp_npools_to_fire_col(i,k) = value_column
       end do
    end do

    do k = 1, ndecomp_pools
       do j = 1, nlevdecomp_full
          do fi = 1,num_column
             i = filter_column(fi)
             this%m_decomp_npools_to_fire_vr_col(i,j,k) = value_column
             this%decomp_npools_transport_tendency_col(i,j,k) = value_column
          end do
       end do
    end do


       do l = 1, ndecomp_cascade_transitions
          do fi = 1,num_column
             i = filter_column(fi)
             this%decomp_cascade_ntransfer_col(i,l) = value_column
             this%decomp_cascade_sminn_flux_col(i,l) = value_column
             if (.not. use_nitrif_denitrif) then
                this%sminn_to_denit_decomp_cascade_col(i,l) = value_column
             end if
          end do
       end do

       do l = 1, ndecomp_cascade_transitions
          do j = 1, nlevdecomp_full
             do fi = 1,num_column
                i = filter_column(fi)
                this%decomp_cascade_ntransfer_vr_col(i,j,l) = value_column
                this%decomp_cascade_sminn_flux_vr_col(i,j,l) = value_column
                if (.not. use_nitrif_denitrif) then
                   this%sminn_to_denit_decomp_cascade_vr_col(i,j,l) = value_column
                end if
             end do
          end do
       end do


    do k = 1, ndecomp_pools
       do j = 1, nlevdecomp_full
          do fi = 1,num_column
             i = filter_column(fi)
             this%decomp_npools_sourcesink_col(i,j,k) = value_column
          end do
       end do
    end do

    ! pflotran
    !------------------------------------------------------------------------
    do k = 1, ndecomp_pools
       do j = 1, nlevdecomp_full
          do fi = 1,num_column
             i = filter_column(fi)
             ! only initializing in the first time-step
             if ( this%externaln_to_decomp_npools_col(i,j,k) == spval ) then
                this%externaln_to_decomp_npools_col(i,j,k) = value_column
             end if
          end do
       end do
    end do

    do j = 1, nlevdecomp_full
       do fi = 1,num_column
          i = filter_column(fi)
          ! only initializing in the first time-step
          if ( this%no3_net_transport_vr_col(i,j) == spval ) then
             this%no3_net_transport_vr_col(i,j) = value_column
          end if
       end do
    end do

    do fi = 1,num_column
       i = filter_column(fi)
       ! only initializing in the first time-step
       if ( this%externaln_to_decomp_delta_col(i) == spval ) then
          this%externaln_to_decomp_delta_col(i) = value_column
       end if
       if ( this%no3_net_transport_delta_col(i) == spval ) then
          this%no3_net_transport_delta_col(i)   = value_column
       end if
    end do
    !------------------------------------------------------------------------

  end subroutine SetValues

  !-----------------------------------------------------------------------
  subroutine ZeroDwt( this, bounds )
    !
    ! !DESCRIPTION
    ! Initialize flux variables needed for dynamic land use.
    !
    ! !ARGUMENTS:
    class(nitrogenflux_type) :: this
    type(bounds_type), intent(in)  :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: c, j          ! indices
    !-----------------------------------------------------------------------

    do c = bounds%begc,bounds%endc
       this%dwt_seedn_to_leaf_col(c)     = 0._r8
       this%dwt_seedn_to_deadstem_col(c) = 0._r8
       this%dwt_conv_nflux_col(c)        = 0._r8
       this%dwt_prod10n_gain_col(c)      = 0._r8
       this%dwt_prod100n_gain_col(c)     = 0._r8
    end do

    do j = 1, nlevdecomp_full
       do c = bounds%begc,bounds%endc
          this%dwt_frootn_to_litr_met_n_col(c,j) = 0._r8
          this%dwt_frootn_to_litr_cel_n_col(c,j) = 0._r8
          this%dwt_frootn_to_litr_lig_n_col(c,j) = 0._r8
          this%dwt_livecrootn_to_cwdn_col(c,j)   = 0._r8
          this%dwt_deadcrootn_to_cwdn_col(c,j)   = 0._r8
       end do
    end do

  end subroutine ZeroDwt

  !-----------------------------------------------------------------------
   subroutine Summary_betr(this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
     !
     ! !USES:
     use clm_varpar    , only: nlevdecomp,ndecomp_cascade_transitions,ndecomp_pools
     use clm_varctl    , only: use_nitrif_denitrif
     use subgridAveMod , only: p2c
     use pftvarcon     , only : npcropmin
     !
     ! !ARGUMENTS:
     class (nitrogenflux_type) :: this
     type(bounds_type) , intent(in) :: bounds
     integer           , intent(in) :: num_soilc       ! number of soil columns in filter
     integer           , intent(in) :: filter_soilc(:) ! filter for soil columns
     integer           , intent(in) :: num_soilp       ! number of soil patches in filter
     integer           , intent(in) :: filter_soilp(:) ! filter for soil patches

     integer :: fc, c, j

     ! total column-level fire N losses
     do fc = 1,num_soilc
        c = filter_soilc(fc)
        this%fire_nloss_col(c) = this%fire_nloss_p2c_col(c) + this%fire_decomp_nloss_col(c)
     end do

     ! supplementary N supplement_to_sminn
     do j = 1, nlevdecomp
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           this%supplement_to_sminn_col(c) = &
                this%supplement_to_sminn_col(c) + &
                this%supplement_to_sminn_vr_col(c,j) * dzsoi_decomp(j)

           this%sminn_input_col(c) = &
                this%sminn_input_col(c) + &
                (this%sminn_nh4_input_vr_col(c,j)+this%sminn_no3_input_vr_col(c,j))*dzsoi_decomp(j)

           this%sminn_nh4_input_col(c) = &
                this%sminn_nh4_input_col(c) + &
                this%sminn_nh4_input_vr_col(c,j)*dzsoi_decomp(j)

           this%sminn_no3_input_col(c) = &
                this%sminn_no3_input_col(c) + &
                this%sminn_no3_input_vr_col(c,j)*dzsoi_decomp(j)
        end do
     end do

     do fc = 1,num_soilc
        c = filter_soilc(fc)

        ! column-level N losses due to landcover change
        this%dwt_nloss_col(c) = &
             this%dwt_conv_nflux_col(c)

        ! total wood product N loss
        this%product_nloss_col(c) = &
             this%prod10n_loss_col(c) + &
             this%prod100n_loss_col(c)+ &
             this%prod1n_loss_col(c)
     end do

   end subroutine Summary_betr
 !-----------------------------------------------------------------------
  subroutine Summary(this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
    !
    ! !USES:
    use clm_varpar    , only: nlevdecomp,ndecomp_cascade_transitions,ndecomp_pools
    use clm_varctl    , only: use_nitrif_denitrif
    use subgridAveMod , only: p2c
    use pftvarcon     , only : npcropmin
    use tracer_varcon , only: is_active_betr_bgc
    !
    ! !ARGUMENTS:
    class (nitrogenflux_type) :: this
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
    !-----------------------------------------------------------------------

    do fp = 1,num_soilp
       p = filter_soilp(fp)

       ! total N deployment (from sminn and retranslocated N pool) (NDEPLOY)
       this%ndeploy_patch(p) = &
            this%sminn_to_npool_patch(p) + &
            this%retransn_to_npool_patch(p)

       ! pft-level wood harvest
       this%wood_harvestn_patch(p) = &
            this%hrv_deadstemn_to_prod10n_patch(p) + &
            this%hrv_deadstemn_to_prod100n_patch(p)
       if ( crop_prog .and. pft%itype(p) >= npcropmin )then
            this%wood_harvestn_patch(p) = &
            this%wood_harvestn_patch(p) + &
            this%hrv_cropn_to_prod1n_patch(p)
       end if

       ! total pft-level fire N losses
       this%fire_nloss_patch(p) = &
            this%m_leafn_to_fire_patch(p)               + &
            this%m_leafn_storage_to_fire_patch(p)       + &
            this%m_leafn_xfer_to_fire_patch(p)          + &
            this%m_frootn_to_fire_patch(p)              + &
            this%m_frootn_storage_to_fire_patch(p)      + &
            this%m_frootn_xfer_to_fire_patch(p)         + &
            this%m_livestemn_to_fire_patch(p)           + &
            this%m_livestemn_storage_to_fire_patch(p)   + &
            this%m_livestemn_xfer_to_fire_patch(p)      + &
            this%m_deadstemn_to_fire_patch(p)           + &
            this%m_deadstemn_storage_to_fire_patch(p)   + &
            this%m_deadstemn_xfer_to_fire_patch(p)      + &
            this%m_livecrootn_to_fire_patch(p)          + &
            this%m_livecrootn_storage_to_fire_patch(p)  + &
            this%m_livecrootn_xfer_to_fire_patch(p)     + &
            this%m_deadcrootn_to_fire_patch(p)          + &
            this%m_deadcrootn_storage_to_fire_patch(p)  + &
            this%m_deadcrootn_xfer_to_fire_patch(p)     + &
            this%m_retransn_to_fire_patch(p)

    end do

    call p2c(bounds, num_soilc, filter_soilc, &
         this%fire_nloss_patch(bounds%begp:bounds%endp), &
         this%fire_nloss_p2c_col(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%wood_harvestn_patch(bounds%begp:bounds%endp), &
         this%wood_harvestn_col(bounds%begc:bounds%endc))
    if(is_active_betr_bgc)then
      call this%Summary_betr(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
      return
    endif

    do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%denit_col(c) = 0._r8
       this%supplement_to_sminn_col(c) = 0._r8
       this%som_n_leached_col(c)       = 0._r8
    end do


    if ( (.not. (use_pflotran .and. pf_cmode)) ) then

       ! BeTR is off AND PFLOTRAN's pf_cmode is false

       ! vertically integrate decomposing N cascade fluxes and
       !soil mineral N fluxes associated with decomposition cascade

       do k = 1, ndecomp_cascade_transitions
          do j = 1,nlevdecomp
             do fc = 1,num_soilc
                c = filter_soilc(fc)

                this%decomp_cascade_ntransfer_col(c,k) = &
                     this%decomp_cascade_ntransfer_col(c,k) + &
                     this%decomp_cascade_ntransfer_vr_col(c,j,k) * dzsoi_decomp(j)

                this%decomp_cascade_sminn_flux_col(c,k) = &
                     this%decomp_cascade_sminn_flux_col(c,k) + &
                     this%decomp_cascade_sminn_flux_vr_col(c,j,k) * dzsoi_decomp(j)
             end do
          end do
       end do

       if (.not. use_nitrif_denitrif) then
          ! vertically integrate each denitrification flux
          do l = 1, ndecomp_cascade_transitions
             do j = 1, nlevdecomp
                do fc = 1,num_soilc
                   c = filter_soilc(fc)
                   this%sminn_to_denit_decomp_cascade_col(c,l) = &
                        this%sminn_to_denit_decomp_cascade_col(c,l) + &
                        this%sminn_to_denit_decomp_cascade_vr_col(c,j,l) * dzsoi_decomp(j)
                end do
             end do
          end do

          ! vertically integrate bulk denitrification and  leaching flux
          do j = 1, nlevdecomp
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                this%sminn_to_denit_excess_col(c) = &
                     this%sminn_to_denit_excess_col(c) + &
                     this%sminn_to_denit_excess_vr_col(c,j) * dzsoi_decomp(j)

                this%sminn_leached_col(c) = &
                     this%sminn_leached_col(c) + &
                     this%sminn_leached_vr_col(c,j) * dzsoi_decomp(j)
             end do
          end do

          ! total N denitrification (DENIT)
          do l = 1, ndecomp_cascade_transitions
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                this%denit_col(c) = &
                     this%denit_col(c) + &
                     this%sminn_to_denit_decomp_cascade_col(c,l)
             end do
          end do

          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%denit_col(c) =  &
                  this%denit_col(c) + &
                  this%sminn_to_denit_excess_col(c)
          end do

       else

          ! vertically integrate NO3 NH4 N2O fluxes and pools
          do j = 1, nlevdecomp
             do fc = 1,num_soilc
                c = filter_soilc(fc)

                ! nitrification and denitrification fluxes
                this%f_nit_col(c) = &
                     this%f_nit_col(c) + &
                     this%f_nit_vr_col(c,j) * dzsoi_decomp(j)

                this%f_denit_col(c) = &
                     this%f_denit_col(c) + &
                     this%f_denit_vr_col(c,j) * dzsoi_decomp(j)

                this%pot_f_nit_col(c) = &
                     this%pot_f_nit_col(c) + &
                     this%pot_f_nit_vr_col(c,j) * dzsoi_decomp(j)

                this%pot_f_denit_col(c) = &
                     this%pot_f_denit_col(c) + &
                     this%pot_f_denit_vr_col(c,j) * dzsoi_decomp(j)

                this%f_n2o_nit_col(c) = &
                     this%f_n2o_nit_col(c) + &
                     this%f_n2o_nit_vr_col(c,j) * dzsoi_decomp(j)

                this%f_n2o_denit_col(c) = &
                     this%f_n2o_denit_col(c) + &
                     this%f_n2o_denit_vr_col(c,j) * dzsoi_decomp(j)

                ! leaching/runoff flux
                this%smin_no3_leached_col(c) = &
                     this%smin_no3_leached_col(c) + &
                     this%smin_no3_leached_vr_col(c,j) * dzsoi_decomp(j)

                this%smin_no3_runoff_col(c) = &
                     this%smin_no3_runoff_col(c) + &
                     this%smin_no3_runoff_vr_col(c,j) * dzsoi_decomp(j)
             end do
          end do

          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%denit_col(c) = this%f_denit_col(c)
          end do

       end if

    endif
    ! vertically integrate column-level fire N losses
    do k = 1, ndecomp_pools
       do j = 1, nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%m_decomp_npools_to_fire_col(c,k) = &
                  this%m_decomp_npools_to_fire_col(c,k) + &
                  this%m_decomp_npools_to_fire_vr_col(c,j,k) * dzsoi_decomp(j)
          end do
       end do
    end do

    ! total column-level fire N losses
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%fire_nloss_col(c) = this%fire_nloss_p2c_col(c)
    end do
    do k = 1, ndecomp_pools
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%fire_nloss_col(c) = &
               this%fire_nloss_col(c) + &
               this%m_decomp_npools_to_fire_col(c,k)
       end do
    end do

    ! supplementary N supplement_to_sminn
    do j = 1, nlevdecomp
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%supplement_to_sminn_col(c) = &
               this%supplement_to_sminn_col(c) + &
               this%supplement_to_sminn_vr_col(c,j) * dzsoi_decomp(j)

          this%sminn_input_col(c) = &
               this%sminn_input_col(c) + &
               (this%sminn_nh4_input_vr_col(c,j)+this%sminn_no3_input_vr_col(c,j))*dzsoi_decomp(j)

          this%sminn_nh4_input_col(c) = &
               this%sminn_nh4_input_col(c) + &
               this%sminn_nh4_input_vr_col(c,j)*dzsoi_decomp(j)

          this%sminn_no3_input_col(c) = &
               this%sminn_no3_input_col(c) + &
               this%sminn_no3_input_vr_col(c,j)*dzsoi_decomp(j)
       end do
    end do

    do fc = 1,num_soilc
       c = filter_soilc(fc)

       ! column-level N losses due to landcover change
       this%dwt_nloss_col(c) = &
            this%dwt_conv_nflux_col(c)

       ! total wood product N loss
       this%product_nloss_col(c) = &
            this%prod10n_loss_col(c) + &
            this%prod100n_loss_col(c)+ &
            this%prod1n_loss_col(c)
    end do

      ! add up all vertical transport tendency terms and calculate total som leaching loss as the sum of these
      do l = 1, ndecomp_pools
        do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%decomp_npools_leached_col(c,l) = 0._r8
        end do

        do j = 1, nlevdecomp
           do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%decomp_npools_leached_col(c,l) = &
                  this%decomp_npools_leached_col(c,l) + &
                  this%decomp_npools_transport_tendency_col(c,j,l) * dzsoi_decomp(j)

           end do
        end do

        do fc = 1,num_soilc
           c = filter_soilc(fc)
           this%som_n_leached_col(c) = &
               this%som_n_leached_col(c) + &
               this%decomp_npools_leached_col(c,l)
        end do
      end do

      do fc = 1,num_soilc
        c = filter_soilc(fc)
        this%smin_no3_to_plant_col(c) = 0._r8
        this%smin_nh4_to_plant_col(c) = 0._r8
        this%plant_to_litter_nflux(c) = 0._r8
        this%plant_to_cwd_nflux(c) = 0._r8
        do j = 1, nlevdecomp
          this%plant_to_litter_nflux(c) = &
               this%plant_to_litter_nflux(c)  + &
               this%phenology_n_to_litr_met_n_col(c,j)* dzsoi_decomp(j) + &
               this%phenology_n_to_litr_cel_n_col(c,j)* dzsoi_decomp(j) + &
               this%phenology_n_to_litr_lig_n_col(c,j)* dzsoi_decomp(j) + &
               this%gap_mortality_n_to_litr_met_n_col(c,j)* dzsoi_decomp(j) + &
               this%gap_mortality_n_to_litr_cel_n_col(c,j)* dzsoi_decomp(j) + &
               this%gap_mortality_n_to_litr_lig_n_col(c,j)* dzsoi_decomp(j) + &
               this%m_n_to_litr_met_fire_col(c,j)* dzsoi_decomp(j) + &
               this%m_n_to_litr_cel_fire_col(c,j)* dzsoi_decomp(j) + &
               this%m_n_to_litr_lig_fire_col(c,j)* dzsoi_decomp(j)
          this%plant_to_cwd_nflux(c) = &
               this%plant_to_cwd_nflux(c) + &
               this%gap_mortality_n_to_cwdn_col(c,j)* dzsoi_decomp(j) + &
               this%fire_mortality_n_to_cwdn_col(c,j)* dzsoi_decomp(j)
        end do
      end do

      if (use_nitrif_denitrif) then
        do fc = 1,num_soilc
          c = filter_soilc(fc)
          do j = 1, nlevdecomp
             this%smin_no3_to_plant_col(c)= this%smin_no3_to_plant_col(c) + &
                  this%smin_no3_to_plant_vr_col(c,j) * dzsoi_decomp(j)
             this%smin_nh4_to_plant_col(c)= this%smin_nh4_to_plant_col(c) + &
                  this%smin_nh4_to_plant_vr_col(c,j) * dzsoi_decomp(j)
          enddo
        enddo
      endif

    ! bgc interface & pflotran
    !----------------------------------------------------------------
    if (use_bgc_interface) then
        call NSummary_interface(this, bounds, num_soilc, filter_soilc)
    end if
    !----------------------------------------------------------------
  end subroutine Summary

!!-------------------------------------------------------------------------------------------------
! !INTERFACE:
subroutine NSummary_interface(this,bounds,num_soilc, filter_soilc)
!
! !DESCRIPTION:
!! bgc interface & pflotran:
! On the radiation time step, perform olumn-level nitrogen
! summary calculations, which mainly from PFLOTRAN bgc coupling
!
! !USES:
   use clm_varpar  , only: nlevdecomp, ndecomp_pools
   use clm_varpar  , only: i_met_lit, i_cel_lit, i_lig_lit, i_cwd
   use clm_time_manager    , only : get_step_size

!   use clm_varctl    , only: pf_hmode
!
! !ARGUMENTS:
   implicit none
   class (nitrogenflux_type)       :: this
   type(bounds_type) ,  intent(in) :: bounds
   integer,             intent(in) :: num_soilc       ! number of soil columns in filter
   integer,             intent(in) :: filter_soilc(:) ! filter for soil columns
!
! !CALLED FROM:
! subroutine NSummary (if pflotran coupled)
!
! !REVISION HISTORY:
!!06/17/2015: modified by Gangsheng Wang
!
! !LOCAL VARIABLES:
   integer :: c,j, l      ! indices
   integer :: fc          ! column filter indices
   real(r8):: dtime             ! radiation time step (seconds)

   ! set time steps
    dtime = real( get_step_size(), r8 )

    if (use_pflotran .and. pf_cmode) then
! nitrification-denitrification rates (not yet passing out from PF, but will)
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         this%f_nit_col(c)   = 0._r8
         this%f_denit_col(c) = 0._r8
         do j = 1, nlevdecomp
            this%f_nit_vr_col(c,j) = 0._r8
            this%f_nit_col(c)  = this%f_nit_col(c) + &
                                 this%f_nit_vr_col(c,j)*dzsoi_decomp(j)

            this%f_denit_vr_col(c,j) = 0._r8
            this%f_denit_col(c) = this%f_denit_col(c) + &
                                 this%f_denit_vr_col(c,j)*dzsoi_decomp(j)

         end do
         this%denit_col(c)      = this%f_denit_col(c)

       end do

       ! the following are from pflotran bgc
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%f_n2_soil_col(c)    = 0._r8
          this%f_n2o_soil_col(c)   = 0._r8
          this%f_ngas_decomp_col(c)= 0._r8
          this%f_ngas_nitri_col(c) = 0._r8
          this%f_ngas_denit_col(c) = 0._r8
          this%smin_no3_leached_col(c) = 0._r8
          this%smin_no3_runoff_col(c)  = 0._r8

          do j = 1, nlevdecomp

            ! all N2/N2O gas exchange between atm. and soil (i.e., dissolving - degassing)
            this%f_n2_soil_col(c)  = this%f_n2_soil_col(c) + &
                                  this%f_n2_soil_vr_col(c,j)*dzsoi_decomp(j)
            this%f_n2o_soil_col(c) = this%f_n2o_soil_col(c) + &
                                  this%f_n2o_soil_vr_col(c,j)*dzsoi_decomp(j)

            ! all N2/N2O production from soil bgc N processes (mineralization-nitrification-denitrification)
            ! note: those are directly dissolved into aq. gas species, which would be exchanging with atm.
            this%f_ngas_decomp_col(c) = this%f_ngas_decomp_col(c) + &
                                     this%f_ngas_decomp_vr_col(c,j)*dzsoi_decomp(j)
            this%f_ngas_nitri_col(c)  = this%f_ngas_nitri_col(c) + &
                                     this%f_ngas_nitri_vr_col(c,j)*dzsoi_decomp(j)
            this%f_ngas_denit_col(c)  = this%f_ngas_denit_col(c) + &
                                     this%f_ngas_denit_vr_col(c,j)*dzsoi_decomp(j)

            ! leaching/runoff flux (if not hydroloy-coupled, from CLM-CN; otherwise from PF)
            this%smin_no3_leached_col(c) = this%smin_no3_leached_col(c) + &
                                        this%smin_no3_leached_vr_col(c,j) * dzsoi_decomp(j)
            this%smin_no3_runoff_col(c)  = this%smin_no3_runoff_col(c) + &
                                        this%smin_no3_runoff_vr_col(c,j) * dzsoi_decomp(j)

            ! assign all no3-N leaching/runoff to all mineral-N
            this%sminn_leached_vr_col(c,j) = this%smin_no3_leached_vr_col(c,j) + &
                                               this%smin_no3_runoff_vr_col(c,j)

          end do

          ! for balance-checking
          this%denit_col(c)     = this%f_ngas_denit_col(c)
          this%f_n2o_nit_col(c) = this%f_ngas_decomp_col(c) + this%f_ngas_nitri_col(c)

          ! assign all no3-N leaching/runoff to all mineral-N
          this%sminn_leached_col(c) = this%smin_no3_leached_col(c) + this%smin_no3_runoff_col(c)

      end do
    end if !! if (use_pflotran .and. pf_cmode)


       ! summarize at column-level vertically-resolved littering/removal for PFLOTRAN bgc input needs
       ! first it needs to save the total column-level N rate btw plant pool and decomposible pools at previous time step
       ! for adjusting difference when doing balance check

       do fc = 1,num_soilc
         c = filter_soilc(fc)
         this%externaln_to_decomp_delta_col(c) = 0._r8
         this%no3_net_transport_delta_col(c)   = 0._r8
         do j = 1, nlevdecomp
            do l = 1, ndecomp_pools
               this%externaln_to_decomp_delta_col(c) =    &
                  this%externaln_to_decomp_delta_col(c) + &
                    this%externaln_to_decomp_npools_col(c,j,l)*dzsoi_decomp(j)
            end do

            ! NO3 leaching/runoff at previous time-step, which may be as source by PFLOTRAN
            this%no3_net_transport_delta_col(c) = &
               this%no3_net_transport_delta_col(c) + &
                    this%no3_net_transport_vr_col(c,j)*dzsoi_decomp(j)

         end do
       end do

       ! do the initialization for the following 2 variables here.
       ! DON'T do so in the beginning of CLM-CN time-step (otherwise the above saved will not work)
       this%externaln_to_decomp_npools_col(:,:,:) = 0._r8
       this%no3_net_transport_vr_col(:,:) = 0._r8

       ! add up all vertically-resolved addition/removal rates (gC/m3/s) of decomp_pools
       do l = 1, ndecomp_pools
          do j = 1, nlevdecomp
             do fc = 1,num_soilc
                c = filter_soilc(fc)

                ! for litter C pools
                if (l==i_met_lit) then
                   this%externaln_to_decomp_npools_col(c,j,l) =              &
                       this%externaln_to_decomp_npools_col(c,j,l)            &
                        + this%phenology_n_to_litr_met_n_col(c,j)            &
                        + this%dwt_frootn_to_litr_met_n_col(c,j)             &
                        + this%gap_mortality_n_to_litr_met_n_col(c,j)        &
                        + this%harvest_n_to_litr_met_n_col(c,j)              !!&
!                        + this%m_n_to_litr_met_fire_col(c,j)                 &
!                        - this%m_decomp_npools_to_fire_vr_col(c,j,l)

                elseif (l==i_cel_lit) then
                   this%externaln_to_decomp_npools_col(c,j,l) =              &
                       this%externaln_to_decomp_npools_col(c,j,l)            &
                        + this%phenology_n_to_litr_cel_n_col(c,j)            &
                        + this%dwt_frootn_to_litr_cel_n_col(c,j)             &
                        + this%gap_mortality_n_to_litr_cel_n_col(c,j)        &
                        + this%harvest_n_to_litr_cel_n_col(c,j)              !!&
!                        + this%m_n_to_litr_cel_fire_col(c,j)                 &
!                        - this%m_decomp_npools_to_fire_vr_col(c,j,l)

                elseif (l==i_lig_lit) then
                   this%externaln_to_decomp_npools_col(c,j,l) =              &
                       this%externaln_to_decomp_npools_col(c,j,l)            &
                        + this%phenology_n_to_litr_lig_n_col(c,j)            &
                        + this%dwt_frootn_to_litr_lig_n_col(c,j)             &
                        + this%gap_mortality_n_to_litr_lig_n_col(c,j)        &
                        + this%harvest_n_to_litr_lig_n_col(c,j)              !!&
!                        + this%m_n_to_litr_lig_fire_col(c,j)                 &
!                        - this%m_decomp_npools_to_fire_vr_col(c,j,l)

                ! for cwd
                elseif (l==i_cwd) then
                   this%externaln_to_decomp_npools_col(c,j,l) =              &
                       this%externaln_to_decomp_npools_col(c,j,l)            &
                        + this%dwt_livecrootn_to_cwdn_col(c,j)               &
                        + this%dwt_deadcrootn_to_cwdn_col(c,j)               &
                        + this%gap_mortality_n_to_cwdn_col(c,j)              &
                        + this%harvest_n_to_cwdn_col(c,j)                    !!&
!                        + this%fire_mortality_n_to_cwdn_col(c,j)

             ! for som n
!                else
!                   this%externaln_to_decomp_npools_col(c,j,l) =              &
!                       this%externaln_to_decomp_npools_col(c,j,l)            &
!                        - this%m_decomp_npools_to_fire_vr_col(c,j,l)

                end if

             ! the following is the net changes of plant N to decompible N poools between time-step
             ! in pflotran, decomposible N pools increments ARE from previous time-step (saved above);
             ! while, in CLM-CN all plant N pools are updated with current N fluxes among plant and ground/soil.
             ! therefore, when do balance check it is needed to adjust the time-lag of changes.
                this%externaln_to_decomp_delta_col(c) =   &
                            this%externaln_to_decomp_delta_col(c) - &
                            this%externaln_to_decomp_npools_col(c,j,l)*dzsoi_decomp(j)

                if (abs(this%externaln_to_decomp_npools_col(c,j,l))<=1.e-21_r8) then
                    this%externaln_to_decomp_npools_col(c,j,l) = 0._r8
                end if

             end do
          end do
       end do

       ! if pflotran hydrology NOT coupled, need to do adjusting for NO3 leaching for balance error checking
       if (.not. pf_hmode) then
          do j = 1, nlevdecomp
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                !! wgs: EXCLUDE leaching from external input
                this%no3_net_transport_vr_col(c,j) = 0._r8
!                this%no3_net_transport_vr_col(c,j) = this%smin_no3_runoff_vr_col(c,j) + &
!                                               this%smin_no3_leached_vr_col(c,j)
                this%no3_net_transport_delta_col(c) = &
                            this%no3_net_transport_delta_col(c) - &
                            this%no3_net_transport_vr_col(c,j)*dzsoi_decomp(j)
             end do
          end do
       end if

       ! change the sign so that it is the increments from the previous time-step (unit: g/m2/s)
       do fc = 1, num_soilc
          c = filter_soilc(fc)
          this%externaln_to_decomp_delta_col(c) = -this%externaln_to_decomp_delta_col(c)
          this%no3_net_transport_delta_col(c)   = -this%no3_net_transport_delta_col(c)
       end do
end subroutine NSummary_interface
!!-------------------------------------------------------------------------------------------------

end module CNNitrogenFluxType
