module CNNitrogenFluxType

  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_infnan_mod         , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use elm_varpar             , only : ndecomp_cascade_transitions, ndecomp_pools
  use elm_varpar             , only : nlevdecomp_full, nlevdecomp, crop_prog
  use elm_varcon             , only : spval, ispval, dzsoi_decomp
  use decompMod              , only : bounds_type
  use clm_varctl             , only : use_nitrif_denitrif, use_vertsoilc
  use CNDecompCascadeConType , only : decomp_cascade_con
  use abortutils             , only : endrun
  use LandunitType           , only : lun_pp                
  use ColumnType             , only : col_pp                
  use VegetationType              , only : veg_pp
  ! bgc interface & pflotran:
  use clm_varctl             , only : use_clm_interface, use_pflotran, pf_cmode, pf_hmode, use_vertsoilc
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
     real(r8), pointer :: m_npool_to_litter_patch                   (:)     ! patch npool mortality (gN/m2/s)
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
     real(r8), pointer :: hrv_npool_to_litter_patch                 (:)     ! patch npool to harvest mortalty (gN/m2/s)
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
     real(r8), pointer :: m_npool_to_fire_patch                     (:)     ! patch (gN/m2/s) fire N emissions from npool
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
     real(r8), pointer :: m_npool_to_litter_fire_patch              (:)     ! patch (gN/m2/s) from npool to litter due to fire
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
     real(r8), pointer :: nfix_to_plantn_patch                      (:)     ! nitrogen fixation goes to plant
     real(r8), pointer :: nfix_to_ecosysn_col                       (:)     ! total nitrogen fixation
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
     real(r8), pointer :: bgc_npool_ext_inputs_vr_col               (:,:,:) !col organic nitrogen input, gN/m3/time step
     real(r8), pointer :: bgc_npool_ext_loss_vr_col                 (:,:,:) !col extneral organic nitrogen loss, gN/m3/time step
     
     real(r8), pointer :: bgc_npool_inputs_col                      (:,:)   !col organic N input, gN/m2/time step
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

     ! crop fluxes
     real(r8), pointer :: crop_seedn_to_leaf_patch                  (:)     ! patch (gN/m2/s) seed source to leaf, for crops

     ! dynamic landcover fluxes
     real(r8), pointer :: dwt_seedn_to_leaf_patch                   (:)     ! (gN/m2/s) seed source to patch-level; although this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_seedn_to_leaf_grc                     (:)     ! (gN/m2/s) dwt_seedn_to_leaf_patch summed to the gridcell-level
     real(r8), pointer :: dwt_seedn_to_deadstem_patch               (:)     ! (gN/m2/s) seed source to patch-level; although this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_seedn_to_deadstem_grc                 (:)     ! (gN/m2/s) dwt_seedn_to_deadstem_patch summed to the gridcell-level
     real(r8), pointer :: dwt_conv_nflux_patch                      (:)     ! (gN/m2/s) conversion N flux (immediate loss to atm); although this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_conv_nflux_grc                        (:)     ! (gN/m2/s) dwt_conv_nflux_patch summed to the gridcell-level
     real(r8), pointer :: dwt_prod10n_gain_patch                    (:)     ! patch (gN/m2/s) addition to 10-yr wood product pool; even though this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_prod100n_gain_patch                   (:)     ! patch (gN/m2/s) addition to 100-yr wood product pool; even though this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_crop_productn_gain_patch              (:)     ! patch (gN/m2/s) addition to crop product pool from landcover change; even though this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_slash_nflux_col                       (:)     ! (gN/m2/s) conversion slash flux due to landcover change

     real(r8), pointer :: dwt_conv_nflux_col                        (:)     ! col (gN/m2/s) conversion N flux (immediate loss to atm)
     real(r8), pointer :: dwt_prod10n_gain_col                      (:)     ! col (gN/m2/s) addition to 10-yr wood product pool
     real(r8), pointer :: dwt_prod100n_gain_col                     (:)     ! col (gN/m2/s) addition to 100-yr wood product pool
     real(r8), pointer :: dwt_frootn_to_litr_met_n_col              (:,:)   ! col (gN/m3/s) fine root to litter due to landcover change
     real(r8), pointer :: dwt_frootn_to_litr_cel_n_col              (:,:)   ! col (gN/m3/s) fine root to litter due to landcover change
     real(r8), pointer :: dwt_frootn_to_litr_lig_n_col              (:,:)   ! col (gN/m3/s) fine root to litter due to landcover change
     real(r8), pointer :: dwt_livecrootn_to_cwdn_col                (:,:)   ! col (gN/m3/s) live coarse root to CWD due to landcover change
     real(r8), pointer :: dwt_deadcrootn_to_cwdn_col                (:,:)   ! col (gN/m3/s) dead coarse root to CWD due to landcover change
     real(r8), pointer :: dwt_nloss_col                             (:)     ! col (gN/m2/s) total nitrogen loss from product pools and conversion

     real(r8), pointer :: dwt_seedn_to_npool_patch                  (:)     ! col (gN/m2/s) seed source to PFT level
     real(r8), pointer :: dwt_seedn_to_npool_grc                    (:)     ! col (gN/m2/s) seed source to PFT level
     real(r8), pointer :: dwt_prod10n_gain_grc                      (:)     ! col (gN/m2/s) addition to 10-yr wood product pool
     real(r8), pointer :: dwt_prod100n_gain_grc                     (:)     ! col (gN/m2/s) addition to 100-yr wood product pool

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

     ! bgc interface/pflotran
     !------------------------------------------------------------------------
     real(r8), pointer :: plant_ndemand_col                         (:)     ! col N flux required to support initial GPP (gN/m2/s)
     ! pflotran
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
     real(r8), pointer :: no3_net_transport_vr_col                  (:,:)   ! col net NO3 transport associated with runoff/leaching/diffusion (gN/m3/s)
     real(r8), pointer :: nh4_net_transport_vr_col                  (:,:)   ! col net NH4 transport associated with runoff/leaching/diffusion (gN/m3/s)
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
     real(r8), pointer :: supplement_to_plantn                      (:)     ! supplementary N flux for plant

     real(r8), pointer :: gap_nloss_litter                          (:)     ! total nloss from veg to litter pool due to gap mortality
     real(r8), pointer :: fire_nloss_litter                         (:)     ! total nloss from veg to litter pool due to fire
     real(r8), pointer :: hrv_nloss_litter                          (:)     ! total nloss from veg to litter pool due to harvest mortality
     real(r8), pointer :: sen_nloss_litter                          (:)     ! total nloss from veg to litter pool due to senescence

     ! C4MIP output variable
     real(r8), pointer :: plant_n_to_cwdn                           (:) ! sum of gap, fire, dynamic land use, and harvest mortality, plant nitrogen flux to CWD

   contains

     procedure , public  :: Init   
     procedure , public  :: Restart
     procedure , public  :: SetValues
     procedure , public  :: ZeroDWT
     procedure , public  :: Summary
     procedure , private :: InitAllocate
     procedure , private :: InitHistory
     procedure , private :: InitCold

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
    integer           :: begg,endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc
    begg = bounds%begg; endg = bounds%endg

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
    allocate(this%m_npool_to_litter_patch                   (begp:endp)) ; this%m_npool_to_litter_patch                   (:) = nan
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
    allocate(this%hrv_npool_to_litter_patch                 (begp:endp)) ; this%hrv_npool_to_litter_patch                 (:) = nan

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
    allocate(this%m_npool_to_fire_patch                     (begp:endp)) ; this%m_npool_to_fire_patch                     (:) = nan

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
    allocate(this%m_npool_to_litter_fire_patch              (begp:endp)) ; this%m_npool_to_litter_fire_patch              (:) = nan

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
    allocate(this%nfix_to_plantn_patch              (begp:endp)) ; this%nfix_to_plantn_patch              (:) = nan

    allocate(this%ndep_to_sminn_col             (begc:endc))    ; this%ndep_to_sminn_col	     (:) = nan
    allocate(this%nfix_to_sminn_col             (begc:endc))    ; this%nfix_to_sminn_col	     (:) = nan
    allocate(this%nfix_to_ecosysn_col           (begc:endc))    ; this%nfix_to_ecosysn_col           (:) = nan
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

    allocate(this%crop_seedn_to_leaf_patch          (begp:endp))                  ; this%crop_seedn_to_leaf_patch  (:)  =nan

    allocate(this%dwt_seedn_to_leaf_patch           (begp:endp))                  ; this%dwt_seedn_to_leaf_patch      (:) =nan
    allocate(this%dwt_seedn_to_leaf_grc             (begg:endg))                  ; this%dwt_seedn_to_leaf_grc        (:) =nan
    allocate(this%dwt_seedn_to_deadstem_patch       (begp:endp))                  ; this%dwt_seedn_to_deadstem_patch  (:) =nan
    allocate(this%dwt_seedn_to_deadstem_grc         (begg:endg))                  ; this%dwt_seedn_to_deadstem_grc    (:) =nan
    allocate(this%dwt_conv_nflux_patch              (begp:endp))                  ; this%dwt_conv_nflux_patch         (:) =nan
    allocate(this%dwt_conv_nflux_grc                (begg:endg))                  ; this%dwt_conv_nflux_grc           (:) =nan
    allocate(this%dwt_prod10n_gain_patch            (begp:endp))                  ; this%dwt_prod10n_gain_patch       (:) =nan
    allocate(this%dwt_prod100n_gain_patch           (begp:endp))                  ; this%dwt_prod100n_gain_patch      (:) =nan
    allocate(this%dwt_crop_productn_gain_patch      (begp:endp))                  ; this%dwt_crop_productn_gain_patch (:) =nan
    allocate(this%dwt_slash_nflux_col               (begc:endc))                  ; this%dwt_slash_nflux_col          (:) =nan

    allocate(this%dwt_conv_nflux_col         (begc:endc))                   ; this%dwt_conv_nflux_col         (:)   = nan
    allocate(this%dwt_prod10n_gain_col       (begc:endc))                   ; this%dwt_prod10n_gain_col       (:)   = nan
    allocate(this%dwt_prod100n_gain_col      (begc:endc))                   ; this%dwt_prod100n_gain_col      (:)   = nan
    allocate(this%dwt_nloss_col              (begc:endc))                   ; this%dwt_nloss_col              (:)   = nan
    allocate(this%wood_harvestn_col          (begc:endc))                   ; this%wood_harvestn_col          (:)   = nan

    allocate(this%dwt_seedn_to_npool_grc     (begg:endg))                   ; this%dwt_seedn_to_npool_grc     (:)   = nan
    allocate(this%dwt_seedn_to_npool_patch   (begp:endp))                   ; this%dwt_seedn_to_npool_patch   (:)   = nan

    allocate(this%dwt_prod10n_gain_grc      (begg:endg))                   ; this%dwt_prod10n_gain_grc       (:)   = nan
    allocate(this%dwt_prod100n_gain_grc     (begg:endg))                   ; this%dwt_prod100n_gain_grc      (:)   = nan

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
    allocate(this%bgc_npool_ext_inputs_vr_col (begc:endc,1:nlevdecomp_full,ndecomp_pools)) 
    this%bgc_npool_ext_inputs_vr_col    (:,:,:) = nan
    allocate(this%bgc_npool_ext_loss_vr_col   (begc:endc,1:nlevdecomp_full,ndecomp_pools)) 
    this%bgc_npool_ext_loss_vr_col      (:,:,:) = nan

    allocate(this%bgc_npool_inputs_col        (begc:endc,ndecomp_pools))     ;this%bgc_npool_inputs_col              (:,:) = nan
     
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

    allocate(this%externaln_to_decomp_npools_col    (begc:endc, 1:nlevdecomp_full, 1:ndecomp_pools))
    this%externaln_to_decomp_npools_col    (:,:,:) = spval
    allocate(this%externaln_to_decomp_delta_col     (begc:endc))
    this%externaln_to_decomp_delta_col     (:)     = spval
    allocate(this%no3_net_transport_vr_col          (begc:endc, 1:nlevdecomp_full))
    this%no3_net_transport_vr_col          (:,:)   = spval
    allocate(this%nh4_net_transport_vr_col          (begc:endc, 1:nlevdecomp_full))
    this%nh4_net_transport_vr_col          (:,:)   = spval
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
    allocate(this%supplement_to_plantn        (begp:endp)) ;             this%supplement_to_plantn  (:)   = 0.d0

    allocate(this%gap_nloss_litter            (begp:endp)) ; this%gap_nloss_litter                  (:) = nan
    allocate(this%fire_nloss_litter           (begp:endp)) ; this%fire_nloss_litter                 (:) = nan
    allocate(this%hrv_nloss_litter            (begp:endp)) ; this%hrv_nloss_litter                  (:) = nan
    allocate(this%sen_nloss_litter            (begp:endp)) ; this%sen_nloss_litter                  (:) = nan

    ! C4MIP output variable
    allocate(this%plant_n_to_cwdn             (begc:endc)) ; this%plant_n_to_cwdn                   (:)  =nan
 
  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use elm_varpar     , only : nlevsno, nlevgrnd, crop_prog 
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
    integer        :: begg, endg
    character(10)  :: active
    character(24)  :: fieldname
    character(100) :: longname
    character(8)   :: vr_suffix
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    ! add suffix if number of soil decomposition depths is greater than 1
    if (nlevdecomp > 1) then
       vr_suffix = "_vr"
    else 
       vr_suffix = ""
    endif

    !-------------------------------
    ! N flux variables - native to column
    !-------------------------------


    !-----------------------------------------------------------
  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initializes time varying variables used only in coupled carbon-nitrogen mode (CN):
    !
    ! !USES:
    use elm_varpar      , only : crop_prog
    use landunit_varcon , only : istsoil, istcrop
    !
    ! !ARGUMENTS:
    class(nitrogenflux_type) :: this 
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: p,c,l, g
    integer :: fp, fc                                    ! filter indices
    integer :: num_special_col                           ! number of good values in special_col filter
    integer :: num_special_patch                         ! number of good values in special_patch filter
    integer :: special_col(bounds%endc-bounds%begc+1)    ! special landunit filter - columns
    integer :: special_patch(bounds%endp-bounds%begp+1)  ! special landunit filter - patches
    !---------------------------------------------------------------------


    !-----------------------------------------------
    ! initialize nitrogen flux variables
    !-----------------------------------------------


    ! initialize fields for special filters

   end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine Restart (this,  bounds, ncid, flag )
    !
    ! !DESCRIPTION: 
    ! Read/write CN restart data for carbon state
    !
    ! !USES:
    use elm_varpar, only : crop_prog
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


    ! pflotran
    !------------------------------------------------------------------------
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
    integer  :: c, g, j          ! indices
    !-----------------------------------------------------------------------

  end subroutine ZeroDwt

 !-----------------------------------------------------------------------
  subroutine Summary(this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
    !
    ! !USES:
    use elm_varpar    , only: nlevdecomp,ndecomp_cascade_transitions,ndecomp_pools
    use clm_varctl    , only: use_nitrif_denitrif
    use subgridAveMod , only: p2c
    use pftvarcon     , only : npcropmin 
    use tracer_varcon , only: is_active_betr_bgc
    use elm_varpar    , only: nlevdecomp_full
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
    integer  :: nlev
    !-----------------------------------------------------------------------

  end subroutine Summary

!-------------------------------------------------------------------------------------------------
! !INTERFACE:
subroutine NSummary_interface(this,bounds,num_soilc, filter_soilc)
!
! !DESCRIPTION:
! bgc interface & pflotran:
! On the radiation time step, perform column-level nitrogen
! summary calculations, which mainly from PFLOTRAN bgc coupling
!
! !USES:
   use elm_varpar  , only: nlevdecomp_full, ndecomp_pools
   use elm_varpar  , only: i_met_lit, i_cel_lit, i_lig_lit, i_cwd
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
! subroutine NSummary (if pflotran coupled) vertically from 1 to 'nlevdecomp_full' (not 'nlevdecomp')
!
!
! !LOCAL VARIABLES:
   integer :: c,j, l      ! indices
   integer :: fc          ! column filter indices
   real(r8):: dtime             ! radiation time step (seconds)

   ! set time steps
    dtime = real( get_step_size(), r8 )
      ! nitrification-denitrification rates (not yet passing out from PF, but will)
      !------------------------------------------------
      ! NOT used currently
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         this%f_nit_col(c)   = 0._r8
         this%f_denit_col(c) = 0._r8
         do j = 1, nlevdecomp_full
            this%f_nit_vr_col(c,j) = 0._r8
            this%f_nit_col(c)  = this%f_nit_col(c) + &
                                 this%f_nit_vr_col(c,j)*dzsoi_decomp(j)

            this%f_denit_vr_col(c,j) = 0._r8
            this%f_denit_col(c) = this%f_denit_col(c) + &
                                 this%f_denit_vr_col(c,j)*dzsoi_decomp(j)

         end do
         this%denit_col(c)      = this%f_denit_col(c)

       end do
      !end------------------------------------------------

       ! the following are from pflotran bgc, and vertically down to 'nlevdecomp_full'
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%f_n2_soil_col(c)    = 0._r8
          this%f_n2o_soil_col(c)   = 0._r8
          this%f_ngas_decomp_col(c)= 0._r8
          this%f_ngas_nitri_col(c) = 0._r8
          this%f_ngas_denit_col(c) = 0._r8
          this%smin_no3_leached_col(c) = 0._r8
          this%smin_no3_runoff_col(c)  = 0._r8
          this%sminn_leached_col(c)    = 0._r8

          do j = 1, nlevdecomp_full

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

            ! leaching/runoff fluxes summed vertically
            ! (1) if not hydroloy-coupled, advection from CLM-CN, plus diffusion from PF
            ! (2) if hydrology-coupled, all from PF (i.e. 'no3_net_transport_vr_col');
            this%smin_no3_leached_col(c) = this%smin_no3_leached_col(c) + &
                                        this%no3_net_transport_vr_col(c,j) * dzsoi_decomp(j)

            if(.not. pf_hmode) then ! this is from CLM-CN's leaching subroutine
                this%smin_no3_leached_col(c) = this%smin_no3_leached_col(c) + &
                                        this%smin_no3_leached_vr_col(c,j) * dzsoi_decomp(j)
                this%smin_no3_runoff_col(c)  = this%smin_no3_runoff_col(c) + &
                                        this%smin_no3_runoff_vr_col(c,j) * dzsoi_decomp(j)
            endif

            ! assign all no3-N leaching/runof,including diffusion from PF, to all mineral-N
            this%sminn_leached_vr_col(c,j) = this%smin_no3_leached_vr_col(c,j) + &
                                             this%smin_no3_runoff_vr_col(c,j) +  &
                                        this%nh4_net_transport_vr_col(c,j) * dzsoi_decomp(j)

            this%sminn_leached_col(c) = this%sminn_leached_col(c) + &
                                        this%sminn_leached_vr_col(c,j)*dzsoi_decomp(j)

          end do !j = 1, nlevdecomp_full

          ! for balance-checking
          this%denit_col(c)     = this%f_ngas_denit_col(c)
          this%f_n2o_nit_col(c) = this%f_ngas_decomp_col(c) + this%f_ngas_nitri_col(c)

      end do !fc = 1,num_soilc


       ! summarize at column-level vertically-resolved littering/removal for PFLOTRAN bgc input needs
       ! first it needs to save the total column-level N rate btw plant pool and decomposible pools at previous time step
       ! for adjusting difference when doing balance check

       do fc = 1,num_soilc
         c = filter_soilc(fc)
         this%externaln_to_decomp_delta_col(c) = 0._r8
         do j = 1, nlevdecomp_full
            do l = 1, ndecomp_pools
               this%externaln_to_decomp_delta_col(c) =    &
                  this%externaln_to_decomp_delta_col(c) + &
                    this%externaln_to_decomp_npools_col(c,j,l)*dzsoi_decomp(j)
            end do

         end do
       end do

       ! do the initialization for the following variable here.
       ! DON'T do so in the beginning of CLM-CN time-step (otherwise the above saved will not work)
       do fc = 1,num_soilc
            c = filter_soilc(fc)
            this%externaln_to_decomp_npools_col(c, 1:nlevdecomp_full, 1:ndecomp_pools) = 0._r8
       end do

       ! add up all vertically-resolved addition/removal rates (gC/m3/s) of decomp_pools
       do fc = 1,num_soilc
            c = filter_soilc(fc)
            do j = 1, nlevdecomp_full
                do l = 1, ndecomp_pools
                ! for litter C pools
                if (l==i_met_lit) then
                   this%externaln_to_decomp_npools_col(c,j,l) =              &
                       this%externaln_to_decomp_npools_col(c,j,l)            &
                        + this%phenology_n_to_litr_met_n_col(c,j)            &
                        + this%dwt_frootn_to_litr_met_n_col(c,j)             &
                        + this%gap_mortality_n_to_litr_met_n_col(c,j)        &
                        + this%harvest_n_to_litr_met_n_col(c,j)              &
                        + this%m_n_to_litr_met_fire_col(c,j)                 

                elseif (l==i_cel_lit) then
                   this%externaln_to_decomp_npools_col(c,j,l) =              &
                       this%externaln_to_decomp_npools_col(c,j,l)            &
                        + this%phenology_n_to_litr_cel_n_col(c,j)            &
                        + this%dwt_frootn_to_litr_cel_n_col(c,j)             &
                        + this%gap_mortality_n_to_litr_cel_n_col(c,j)        &
                        + this%harvest_n_to_litr_cel_n_col(c,j)              &
                        + this%m_n_to_litr_cel_fire_col(c,j)                  

                elseif (l==i_lig_lit) then
                   this%externaln_to_decomp_npools_col(c,j,l) =              &
                       this%externaln_to_decomp_npools_col(c,j,l)            &
                        + this%phenology_n_to_litr_lig_n_col(c,j)            &
                        + this%dwt_frootn_to_litr_lig_n_col(c,j)             &
                        + this%gap_mortality_n_to_litr_lig_n_col(c,j)        &
                        + this%harvest_n_to_litr_lig_n_col(c,j)              &
                        + this%m_n_to_litr_lig_fire_col(c,j)                 

                ! for cwd
                elseif (l==i_cwd) then
                   this%externaln_to_decomp_npools_col(c,j,l) =              &
                       this%externaln_to_decomp_npools_col(c,j,l)            &
                        + this%dwt_livecrootn_to_cwdn_col(c,j)               &
                        + this%dwt_deadcrootn_to_cwdn_col(c,j)               &
                        + this%gap_mortality_n_to_cwdn_col(c,j)              &
                        + this%harvest_n_to_cwdn_col(c,j)                    &
                        + this%fire_mortality_n_to_cwdn_col(c,j)

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

             end do !l = 1, ndecomp_pools
          end do !j = 1, nlevdecomp_full
       end do !fc = 1,num_soilc


       ! if pflotran hydrology NOT coupled, need to do:
       ! saving for (next time-step) possible including of RT mass-transfer in PFLOTRAN bgc coupling.
       ! (NOT USED anymore - 04/26/2017)
       if (.not. pf_hmode) then
          do j = 1, nlevdecomp_full
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                this%no3_net_transport_vr_col(c,j) = this%smin_no3_runoff_vr_col(c,j) + &
                                               this%smin_no3_leached_vr_col(c,j)
             end do
          end do
       else
          do j = 1, nlevdecomp_full
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                this%no3_net_transport_vr_col(c,j) = 0._r8
             end do
          end do
       end if

       ! change the sign so that it is the increments from the previous time-step (unit: g/m2/s)
       do fc = 1, num_soilc
          c = filter_soilc(fc)
          this%externaln_to_decomp_delta_col(c) = -this%externaln_to_decomp_delta_col(c)
       end do

end subroutine NSummary_interface
!-------------------------------------------------------------------------------------------------

end module CNNitrogenFluxType

