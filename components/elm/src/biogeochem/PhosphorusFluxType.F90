module PhosphorusFluxType

  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_infnan_mod         , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use elm_varpar             , only : ndecomp_cascade_transitions, ndecomp_pools
  use elm_varpar             , only : nlevdecomp_full, nlevdecomp, crop_prog
  use elm_varcon             , only : spval, ispval, dzsoi_decomp
  use decompMod              , only : bounds_type
  use elm_varctl             , only : use_nitrif_denitrif, use_vertsoilc
  use CNDecompCascadeConType , only : decomp_cascade_con
  use abortutils             , only : endrun
  use LandunitType           , only : lun_pp                
  use ColumnType             , only : col_pp                
  use VegetationType              , only : veg_pp
  ! bgc interface & pflotran:
  use elm_varctl             , only : use_clm_interface, use_pflotran, pf_cmode, pf_hmode, use_vertsoilc
  ! 
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  type, public :: phosphorusflux_type

     ! gap mortality fluxes
     real(r8), pointer :: m_leafp_to_litter_patch                   (:)     ! patch leaf P mortality (gP/m2/s)
     real(r8), pointer :: m_frootp_to_litter_patch                  (:)     ! patch fine root P mortality (gP/m2/s)
     real(r8), pointer :: m_leafp_storage_to_litter_patch           (:)     ! patch leaf P storage mortality (gP/m2/s)
     real(r8), pointer :: m_frootp_storage_to_litter_patch          (:)     ! patch fine root P storage mortality (gP/m2/s)
     real(r8), pointer :: m_livestemp_storage_to_litter_patch       (:)     ! patch live stem P storage mortality (gP/m2/s)
     real(r8), pointer :: m_deadstemp_storage_to_litter_patch       (:)     ! patch dead stem P storage mortality (gP/m2/s)
     real(r8), pointer :: m_livecrootp_storage_to_litter_patch      (:)     ! patch live coarse root P storage mortality (gP/m2/s)
     real(r8), pointer :: m_deadcrootp_storage_to_litter_patch      (:)     ! patch dead coarse root P storage mortality (gP/m2/s)
     real(r8), pointer :: m_leafp_xfer_to_litter_patch              (:)     ! patch leaf P transfer mortality (gP/m2/s)
     real(r8), pointer :: m_frootp_xfer_to_litter_patch             (:)     ! patch fine root P transfer mortality (gP/m2/s)
     real(r8), pointer :: m_livestemp_xfer_to_litter_patch          (:)     ! patch live stem P transfer mortality (gP/m2/s)
     real(r8), pointer :: m_deadstemp_xfer_to_litter_patch          (:)     ! patch dead stem P transfer mortality (gP/m2/s)
     real(r8), pointer :: m_livecrootp_xfer_to_litter_patch         (:)     ! patch live coarse root P transfer mortality (gP/m2/s)
     real(r8), pointer :: m_deadcrootp_xfer_to_litter_patch         (:)     ! patch dead coarse root P transfer mortality (gP/m2/s)
     real(r8), pointer :: m_livestemp_to_litter_patch               (:)     ! patch live stem P mortality (gP/m2/s)
     real(r8), pointer :: m_deadstemp_to_litter_patch               (:)     ! patch dead stem P mortality (gP/m2/s)
     real(r8), pointer :: m_livecrootp_to_litter_patch              (:)     ! patch live coarse root P mortality (gP/m2/s)
     real(r8), pointer :: m_deadcrootp_to_litter_patch              (:)     ! patch dead coarse root P mortality (gP/m2/s)
     real(r8), pointer :: m_retransp_to_litter_patch                (:)     ! patch retranslocated P pool mortality (gP/m2/s)
     real(r8), pointer :: m_ppool_to_litter_patch                   (:)     ! patch storage P pool mortality (gP/m2/s)

     ! harvest fluxes
     real(r8), pointer :: hrv_leafp_to_litter_patch                 (:)     ! patch leaf P harvest mortality (gP/m2/s)
     real(r8), pointer :: hrv_frootp_to_litter_patch                (:)     ! patch fine root P harvest mortality (gP/m2/s)
     real(r8), pointer :: hrv_leafp_storage_to_litter_patch         (:)     ! patch leaf P storage harvest mortality (gP/m2/s)
     real(r8), pointer :: hrv_frootp_storage_to_litter_patch        (:)     ! patch fine root P storage harvest mortality (gP/m2/s)
     real(r8), pointer :: hrv_livestemp_storage_to_litter_patch     (:)     ! patch live stem P storage harvest mortality (gP/m2/s)
     real(r8), pointer :: hrv_deadstemp_storage_to_litter_patch     (:)     ! patch dead stem P storage harvest mortality (gP/m2/s)
     real(r8), pointer :: hrv_livecrootp_storage_to_litter_patch    (:)     ! patch live coarse root P storage harvest mortality (gP/m2/s)
     real(r8), pointer :: hrv_deadcrootp_storage_to_litter_patch    (:)     ! patch dead coarse root P storage harvest mortality (gP/m2/s)
     real(r8), pointer :: hrv_leafp_xfer_to_litter_patch            (:)     ! patch leaf P transfer harvest mortality (gP/m2/s)
     real(r8), pointer :: hrv_frootp_xfer_to_litter_patch           (:)     ! patch fine root P transfer harvest mortality (gP/m2/s)
     real(r8), pointer :: hrv_livestemp_xfer_to_litter_patch        (:)     ! patch live stem P transfer harvest mortality (gP/m2/s)
     real(r8), pointer :: hrv_deadstemp_xfer_to_litter_patch        (:)     ! patch dead stem P transfer harvest mortality (gP/m2/s)
     real(r8), pointer :: hrv_livecrootp_xfer_to_litter_patch       (:)     ! patch live coarse root P transfer harvest mortality (gP/m2/s)
     real(r8), pointer :: hrv_deadcrootp_xfer_to_litter_patch       (:)     ! patch dead coarse root P transfer harvest mortality (gP/m2/s)
     real(r8), pointer :: hrv_livestemp_to_litter_patch             (:)     ! patch live stem P harvest mortality (gP/m2/s)
     real(r8), pointer :: hrv_deadstemp_to_prod10p_patch            (:)     ! patch dead stem P harvest to 10-year product pool (gP/m2/s)
     real(r8), pointer :: hrv_deadstemp_to_prod100p_patch           (:)     ! patch dead stem P harvest to 100-year product pool (gP/m2/s)
     real(r8), pointer :: hrv_livecrootp_to_litter_patch            (:)     ! patch live coarse root P harvest mortality (gP/m2/s)
     real(r8), pointer :: hrv_deadcrootp_to_litter_patch            (:)     ! patch dead coarse root P harvest mortality (gP/m2/s)
     real(r8), pointer :: hrv_retransp_to_litter_patch              (:)     ! patch retranslocated P pool harvest mortality (gP/m2/s)
     real(r8), pointer :: hrv_ppool_to_litter_patch                 (:)     ! patch retranslocated P pool harvest mortality (gP/m2/s)
     real(r8), pointer :: hrv_deadstemp_to_prod10p_col              (:)     ! col dead stem P harvest mortality to 10-year product pool (gP/m2/s)
     real(r8), pointer :: hrv_deadstemp_to_prod100p_col             (:)     ! col dead stem P harvest mortality to 100-year product pool (gP/m2/s)
     real(r8), pointer :: m_p_to_litr_met_fire_col                  (:,:)   ! col P from leaf, froot, xfer and storage P to litter labile P by fire (gP/m3/s) 
     real(r8), pointer :: m_p_to_litr_cel_fire_col                  (:,:)   ! col P from leaf, froot, xfer and storage P to litter cellulose P by fire (gP/m3/s) 
     real(r8), pointer :: m_p_to_litr_lig_fire_col                  (:,:)   ! col P from leaf, froot, xfer and storage P to litter lignin P by fire (gP/m3/s) 
     real(r8), pointer :: harvest_p_to_litr_met_p_col               (:,:)   ! col P fluxes associated with harvest to litter metabolic pool (gP/m3/s)
     real(r8), pointer :: harvest_p_to_litr_cel_p_col               (:,:)   ! col P fluxes associated with harvest to litter cellulose pool (gP/m3/s)
     real(r8), pointer :: harvest_p_to_litr_lig_p_col               (:,:)   ! col P fluxes associated with harvest to litter lignin pool (gP/m3/s)
     real(r8), pointer :: harvest_p_to_cwdp_col                     (:,:)   ! col P fluxes associated with harvest to CWD pool (gP/m3/s)

     ! crop harvest
     real(r8), pointer :: hrv_leafp_to_prod1p_patch                 (:)     ! crop leafp harvested (gP/m2/s)
     real(r8), pointer :: hrv_livestemp_to_prod1p_patch             (:)     ! crop stemp harvested (gP/m2/s)
     real(r8), pointer :: hrv_grainp_to_prod1p_patch                (:)     ! crop grain harvested (gP/m2/s)
     real(r8), pointer :: hrv_cropp_to_prod1p_patch                 (:)     ! total amount of crop P harvested (gP/m2/s)
     real(r8), pointer :: hrv_cropp_to_prod1p_col                   (:)     ! crop P harvest mortality to 1-yr product pool (gP/m2/s)

     ! fire P fluxes 
     real(r8), pointer :: m_decomp_ppools_to_fire_vr_col            (:,:,:) ! col vertically-resolved decomposing P fire loss (gP/m3/s)
     real(r8), pointer :: m_decomp_ppools_to_fire_col               (:,:)   ! col vertically-integrated (diagnostic) decomposing P fire loss (gP/m2/s)
     real(r8), pointer :: m_leafp_to_fire_patch                     (:)     ! patch (gP/m2/s) fire P emissions from leafp 
     real(r8), pointer :: m_leafp_storage_to_fire_patch             (:)     ! patch (gP/m2/s) fire P emissions from leafp_storage            
     real(r8), pointer :: m_leafp_xfer_to_fire_patch                (:)     ! patch (gP/m2/s) fire P emissions from leafp_xfer     
     real(r8), pointer :: m_livestemp_to_fire_patch                 (:)     ! patch (gP/m2/s) fire P emissions from livestemp 
     real(r8), pointer :: m_livestemp_storage_to_fire_patch         (:)     ! patch (gP/m2/s) fire P emissions from livestemp_storage      
     real(r8), pointer :: m_livestemp_xfer_to_fire_patch            (:)     ! patch (gP/m2/s) fire P emissions from livestemp_xfer
     real(r8), pointer :: m_deadstemp_to_fire_patch                 (:)     ! patch (gP/m2/s) fire P emissions from deadstemp
     real(r8), pointer :: m_deadstemp_storage_to_fire_patch         (:)     ! patch (gP/m2/s) fire P emissions from deadstemp_storage         
     real(r8), pointer :: m_deadstemp_xfer_to_fire_patch            (:)     ! patch (gP/m2/s) fire P emissions from deadstemp_xfer
     real(r8), pointer :: m_frootp_to_fire_patch                    (:)     ! patch (gP/m2/s) fire P emissions from frootp
     real(r8), pointer :: m_frootp_storage_to_fire_patch            (:)     ! patch (gP/m2/s) fire P emissions from frootp_storage
     real(r8), pointer :: m_frootp_xfer_to_fire_patch               (:)     ! patch (gP/m2/s) fire P emissions from frootp_xfer
     real(r8), pointer :: m_livecrootp_to_fire_patch                (:)     ! patch (gP/m2/s) fire P emissions from m_livecrootp_to_fire
     real(r8), pointer :: m_livecrootp_storage_to_fire_patch        (:)     ! patch (gP/m2/s) fire P emissions from livecrootp_storage     
     real(r8), pointer :: m_livecrootp_xfer_to_fire_patch           (:)     ! patch (gP/m2/s) fire P emissions from livecrootp_xfer
     real(r8), pointer :: m_deadcrootp_to_fire_patch                (:)     ! patch (gP/m2/s) fire P emissions from deadcrootp
     real(r8), pointer :: m_deadcrootp_storage_to_fire_patch        (:)     ! patch (gP/m2/s) fire P emissions from deadcrootp_storage  
     real(r8), pointer :: m_deadcrootp_xfer_to_fire_patch           (:)     ! patch (gP/m2/s) fire P emissions from deadcrootp_xfer
     real(r8), pointer :: m_retransp_to_fire_patch                  (:)     ! patch (gP/m2/s) fire P emissions from retransp
     real(r8), pointer :: m_ppool_to_fire_patch                     (:)     ! patch (gP/m2/s) fire P emissions from ppool
     real(r8), pointer :: m_leafp_to_litter_fire_patch              (:)     ! patch (gP/m2/s) from leafp to litter P  due to fire               
     real(r8), pointer :: m_leafp_storage_to_litter_fire_patch      (:)     ! patch (gP/m2/s) from leafp_storage to litter P  due to fire                              
     real(r8), pointer :: m_leafp_xfer_to_litter_fire_patch         (:)     ! patch (gP/m2/s) from leafp_xfer to litter P  due to fire                              
     real(r8), pointer :: m_livestemp_to_litter_fire_patch          (:)     ! patch (gP/m2/s) from livestemp to litter P  due to fire                              
     real(r8), pointer :: m_livestemp_storage_to_litter_fire_patch  (:)     ! patch (gP/m2/s) from livestemp_storage to litter P  due to fire                                     
     real(r8), pointer :: m_livestemp_xfer_to_litter_fire_patch     (:)     ! patch (gP/m2/s) from livestemp_xfer to litter P  due to fire                                     
     real(r8), pointer :: m_livestemp_to_deadstemp_fire_patch       (:)     ! patch (gP/m2/s) from livestemp to deadstemp P  due to fire                                     
     real(r8), pointer :: m_deadstemp_to_litter_fire_patch          (:)     ! patch (gP/m2/s) from deadstemp to litter P  due to fire                                     
     real(r8), pointer :: m_deadstemp_storage_to_litter_fire_patch  (:)     ! patch (gP/m2/s) from deadstemp_storage to litter P  due to fire                                               
     real(r8), pointer :: m_deadstemp_xfer_to_litter_fire_patch     (:)     ! patch (gP/m2/s) from deadstemp_xfer to litter P  due to fire                                               
     real(r8), pointer :: m_frootp_to_litter_fire_patch             (:)     ! patch (gP/m2/s) from frootp to litter P  due to fire                                               
     real(r8), pointer :: m_frootp_storage_to_litter_fire_patch     (:)     ! patch (gP/m2/s) from frootp_storage to litter P  due to fire                                               
     real(r8), pointer :: m_frootp_xfer_to_litter_fire_patch        (:)     ! patch (gP/m2/s) from frootp_xfer to litter P  due to fire                                               
     real(r8), pointer :: m_livecrootp_to_litter_fire_patch         (:)     ! patch (gP/m2/s) from livecrootp to litter P  due to fire                                               
     real(r8), pointer :: m_livecrootp_storage_to_litter_fire_patch (:)     ! patch (gP/m2/s) from livecrootp_storage to litter P  due to fire                                                     
     real(r8), pointer :: m_livecrootp_xfer_to_litter_fire_patch    (:)     ! patch (gP/m2/s) from livecrootp_xfer to litter P  due to fire                                                     
     real(r8), pointer :: m_livecrootp_to_deadcrootp_fire_patch     (:)     ! patch (gP/m2/s) from livecrootp_xfer to deadcrootp due to fire                                                     
     real(r8), pointer :: m_deadcrootp_to_litter_fire_patch         (:)     ! patch (gP/m2/s) from deadcrootp to deadcrootp due to fire                                                       
     real(r8), pointer :: m_deadcrootp_storage_to_litter_fire_patch (:)     ! patch (gP/m2/s) from deadcrootp_storage to deadcrootp due to fire                                                        
     real(r8), pointer :: m_deadcrootp_xfer_to_litter_fire_patch    (:)     ! patch (gP/m2/s) from deadcrootp_xfer to deadcrootp due to fire 
     real(r8), pointer :: m_retransp_to_litter_fire_patch           (:)     ! patch (gP/m2/s) from retransp to deadcrootp due to fire                                                               
     real(r8), pointer :: m_ppool_to_litter_fire_patch              (:)     ! patch (gP/m2/s) from ppool to deadcrootp due to fire                                                         
     real(r8), pointer :: fire_ploss_patch                          (:)     ! patch total pft-level fire P loss (gP/m2/s) 
     real(r8), pointer :: fire_ploss_col                            (:)     ! col total column-level fire P loss (gP/m2/s)
     real(r8), pointer :: fire_decomp_ploss_col                     (:)     ! col fire p loss from decomposable pools (gP/m2/s)
     real(r8), pointer :: fire_ploss_p2c_col                        (:)     ! col patch2col column-level fire P loss (gP/m2/s) (p2c)
     real(r8), pointer :: fire_mortality_p_to_cwdp_col              (:,:)   ! col P fluxes associated with fire mortality to CWD pool (gP/m3/s)

     ! phenology fluxes from transfer pool
     real(r8), pointer :: grainp_xfer_to_grainp_patch               (:)     ! patch grain P growth from storage for prognostic crop model (gP/m2/s)
     real(r8), pointer :: leafp_xfer_to_leafp_patch                 (:)     ! patch leaf P growth from storage (gP/m2/s)
     real(r8), pointer :: frootp_xfer_to_frootp_patch               (:)     ! patch fine root P growth from storage (gP/m2/s)
     real(r8), pointer :: livestemp_xfer_to_livestemp_patch         (:)     ! patch live stem P growth from storage (gP/m2/s)
     real(r8), pointer :: deadstemp_xfer_to_deadstemp_patch         (:)     ! patch dead stem P growth from storage (gP/m2/s)
     real(r8), pointer :: livecrootp_xfer_to_livecrootp_patch       (:)     ! patch live coarse root P growth from storage (gP/m2/s)
     real(r8), pointer :: deadcrootp_xfer_to_deadcrootp_patch       (:)     ! patch dead coarse root P growth from storage (gP/m2/s)

     ! litterfall fluxes
     real(r8), pointer :: livestemp_to_litter_patch                 (:)     ! patch livestem P to litter (gP/m2/s)
     real(r8), pointer :: grainp_to_food_patch                      (:)     ! patch grain P to food for prognostic crop (gP/m2/s)
     real(r8), pointer :: leafp_to_litter_patch                     (:)     ! patch leaf P litterfall (gP/m2/s)
     real(r8), pointer :: leafp_to_retransp_patch                   (:)     ! patch leaf P to retranslocated P pool (gP/m2/s)
     real(r8), pointer :: frootp_to_retransp_patch                  (:)     ! patch fine root P to retranslocated P pool (gP/m2/s)
     real(r8), pointer :: frootp_to_litter_patch                    (:)     ! patch fine root P litterfall (gP/m2/s)

     ! allocation fluxes
     real(r8), pointer :: retransp_to_ppool_patch                   (:)     ! patch deployment of retranslocated P (gP/m2/s)       
     real(r8), pointer :: sminp_to_ppool_patch                      (:)     ! patch deployment of soil mineral P uptake (gP/m2/s)
     real(r8), pointer :: ppool_to_grainp_patch                     (:)     ! patch allocation to grain P for prognostic crop (gP/m2/s)
     real(r8), pointer :: ppool_to_grainp_storage_patch             (:)     ! patch allocation to grain P storage for prognostic crop (gP/m2/s)
     real(r8), pointer :: ppool_to_leafp_patch                      (:)     ! patch allocation to leaf P (gP/m2/s)
     real(r8), pointer :: ppool_to_leafp_storage_patch              (:)     ! patch allocation to leaf P storage (gP/m2/s)
     real(r8), pointer :: ppool_to_frootp_patch                     (:)     ! patch allocation to fine root P (gP/m2/s)
     real(r8), pointer :: ppool_to_frootp_storage_patch             (:)     ! patch allocation to fine root P storage (gP/m2/s)
     real(r8), pointer :: ppool_to_livestemp_patch                  (:)     ! patch allocation to live stem P (gP/m2/s)
     real(r8), pointer :: ppool_to_livestemp_storage_patch          (:)     ! patch allocation to live stem P storage (gP/m2/s)
     real(r8), pointer :: ppool_to_deadstemp_patch                  (:)     ! patch allocation to dead stem P (gP/m2/s)
     real(r8), pointer :: ppool_to_deadstemp_storage_patch          (:)     ! patch allocation to dead stem P storage (gP/m2/s)
     real(r8), pointer :: ppool_to_livecrootp_patch                 (:)     ! patch allocation to live coarse root P (gP/m2/s)
     real(r8), pointer :: ppool_to_livecrootp_storage_patch         (:)     ! patch allocation to live coarse root P storage (gP/m2/s)
     real(r8), pointer :: ppool_to_deadcrootp_patch                 (:)     ! patch allocation to dead coarse root P (gP/m2/s)
     real(r8), pointer :: ppool_to_deadcrootp_storage_patch         (:)     ! patch allocation to dead coarse root P storage (gP/m2/s)

     ! annual turnover of storage to transfer pools           
     real(r8), pointer :: grainp_storage_to_xfer_patch              (:)     ! patch grain P shift storage to transfer for prognostic crop (gP/m2/s)
     real(r8), pointer :: leafp_storage_to_xfer_patch               (:)     ! patch leaf P shift storage to transfer (gP/m2/s)
     real(r8), pointer :: frootp_storage_to_xfer_patch              (:)     ! patch fine root P shift storage to transfer (gP/m2/s)
     real(r8), pointer :: livestemp_storage_to_xfer_patch           (:)     ! patch live stem P shift storage to transfer (gP/m2/s)
     real(r8), pointer :: deadstemp_storage_to_xfer_patch           (:)     ! patch dead stem P shift storage to transfer (gP/m2/s)
     real(r8), pointer :: livecrootp_storage_to_xfer_patch          (:)     ! patch live coarse root P shift storage to transfer (gP/m2/s)
     real(r8), pointer :: deadcrootp_storage_to_xfer_patch          (:)     ! patch dead coarse root P shift storage to transfer (gP/m2/s)
     real(r8), pointer :: fert_p_patch                                (:)     ! patch applied fertilizer (gP/m2/s)
     real(r8), pointer :: fert_p_counter_patch                        (:)     ! patch >0 fertilize; <=0 not

     ! turnover of livewood to deadwood, with retranslocation 
     real(r8), pointer :: livestemp_to_deadstemp_patch              (:)     ! patch live stem P turnover (gP/m2/s)
     real(r8), pointer :: livestemp_to_retransp_patch               (:)     ! patch live stem P to retranslocated P pool (gP/m2/s)
     real(r8), pointer :: livecrootp_to_deadcrootp_patch            (:)     ! patch live coarse root P turnover (gP/m2/s)
     real(r8), pointer :: livecrootp_to_retransp_patch              (:)     ! patch live coarse root P to retranslocated P pool (gP/m2/s)

     ! summary (diagnostic) flux variables, not involved in mass balance
     real(r8), pointer :: pdeploy_patch                             (:)     ! patch total P deployed to growth and storage (gP/m2/s)
     real(r8), pointer :: pinputs_patch                             (:)     ! patch total P inputs to pft-level (gP/m2/s)
     real(r8), pointer :: poutputs_patch                            (:)     ! patch total P outputs from pft-level (gP/m2/s)
     real(r8), pointer :: wood_harvestp_patch                       (:)     ! patch total P losses to wood product pools (gP/m2/s)
     real(r8), pointer :: wood_harvestp_col                         (:)     ! col total P losses to wood product pools (gP/m2/s) (p2c)

     ! deposition fluxes
     real(r8), pointer :: pdep_to_sminp_col                         (:)     ! col atmospheric P deposition to soil mineral P (gP/m2/s)
     real(r8), pointer :: fert_p_to_sminp_col                            (:)     ! col fertilizer P to soil mineral P (gP/m2/s)

     ! phenology: litterfall and crop fluxes
     real(r8), pointer :: phenology_p_to_litr_met_p_col             (:,:)   ! col P fluxes associated with phenology (litterfall and crop) to litter metabolic pool (gP/m3/s)
     real(r8), pointer :: phenology_p_to_litr_cel_p_col             (:,:)   ! col P fluxes associated with phenology (litterfall and crop) to litter cellulose pool (gP/m3/s)
     real(r8), pointer :: phenology_p_to_litr_lig_p_col             (:,:)   ! col P fluxes associated with phenology (litterfall and crop) to litter lignin pool (gP/m3/s)

     ! gap mortality fluxes
     real(r8), pointer :: gap_mortality_p_to_litr_met_p_col         (:,:)   ! col P fluxes associated with gap mortality to litter metabolic pool (gP/m3/s)
     real(r8), pointer :: gap_mortality_p_to_litr_cel_p_col         (:,:)   ! col P fluxes associated with gap mortality to litter cellulose pool (gP/m3/s)
     real(r8), pointer :: gap_mortality_p_to_litr_lig_p_col         (:,:)   ! col P fluxes associated with gap mortality to litter lignin pool (gP/m3/s)
     real(r8), pointer :: gap_mortality_p_to_cwdp_col               (:,:)   ! col P fluxes associated with gap mortality to CWD pool (gP/m3/s)

     ! decomposition fluxes
     real(r8), pointer :: decomp_cascade_ptransfer_vr_col           (:,:,:) ! col vert-res transfer of P from donor to receiver pool along decomp. cascade (gP/m3/s)
     real(r8), pointer :: decomp_cascade_ptransfer_col              (:,:)   ! col vert-int (diagnostic) transfer of P from donor to receiver pool along decomp. cascade (gP/m2/s)
     real(r8), pointer :: decomp_cascade_sminp_flux_vr_col          (:,:,:) ! col vert-res mineral P flux for transition along decomposition cascade (gP/m3/s)
     real(r8), pointer :: decomp_cascade_sminp_flux_col             (:,:)   ! col vert-int (diagnostic) mineral P flux for transition along decomposition cascade (gP/m2/s)
                                                                            ! Used to update concentrations concurrently with vertical transport
     ! vertically-resolved immobilization fluxes
     real(r8), pointer :: potential_immob_p_vr_col                    (:,:)   ! col vertically-resolved potential P immobilization (gP/m3/s) at each level
     real(r8), pointer :: potential_immob_p_col                       (:)     ! col vert-int (diagnostic) potential P immobilization (gP/m2/s)
     real(r8), pointer :: actual_immob_p_vr_col                       (:,:)   ! col vertically-resolved actual P immobilization (gP/m3/s) at each level
     real(r8), pointer :: actual_immob_p_col                          (:)     ! col vert-int (diagnostic) actual P immobilization (gP/m2/s)
     real(r8), pointer :: sminp_to_plant_vr_col                     (:,:)   ! col vertically-resolved plant uptake of soil mineral P (gP/m3/s)
     real(r8), pointer :: sminp_to_plant_col                        (:)     ! col vert-int (diagnostic) plant uptake of soil mineral P (gP/m2/s)
     real(r8), pointer :: supplement_to_sminp_vr_col                (:,:)   ! col vertically-resolved supplemental P supply (gP/m3/s)
     real(r8), pointer :: supplement_to_sminp_col                   (:)     ! col vert-int (diagnostic) supplemental P supply (gP/m2/s)
     real(r8), pointer :: gross_pmin_vr_col                         (:,:)   ! col vertically-resolved gross rate of P mineralization (gP/m3/s)
     real(r8), pointer :: gross_pmin_col                            (:)     ! col vert-int (diagnostic) gross rate of P mineralization (gP/m2/s)
     real(r8), pointer :: net_pmin_vr_col                           (:,:)   ! col vertically-resolved net rate of P mineralization (gP/m3/s)
     real(r8), pointer :: net_pmin_col                              (:)     ! col vert-int (diagnostic) net rate of P mineralization (gP/m2/s)

     real(r8), pointer :: biochem_pmin_ppools_vr_col                (:,:,:) ! col vertically-resolved biochemical P mineralization for each soi pool (gP/m3/s)
     real(r8), pointer :: biochem_pmin_vr_col                       (:,:)   ! col vertically-resolved total biochemical P mineralization (gP/m3/s)
     real(r8), pointer :: biochem_pmin_to_plant_patch               (:)     ! biochemical P mineralization directly goes to plant (gP/m2/s)
     real(r8), pointer :: biochem_pmin_to_ecosysp_vr_col            (:,:)   ! biochemical P mineralization directly goes to soil (gP/m3/s)
     real(r8), pointer :: biochem_pmin_col                          (:)     ! col vert-int (diagnostic) total biochemical P mineralization (gP/m3/s)


     ! new variables for phosphorus code
     ! inorganic P transformation fluxes
     real(r8), pointer :: primp_to_labilep_vr_col                     (:,:)   ! col (gP/m3/s) flux of P from primary mineral to labile 
     real(r8), pointer :: primp_to_labilep_col                        (:)     ! col (gP/m3/s) flux of P from primary mineral to labile 
     real(r8), pointer :: labilep_to_secondp_vr_col                   (:,:)   ! col (gP/m3/s) flux of labile P to secondary mineral P 
     real(r8), pointer :: labilep_to_secondp_col                      (:)     ! col (gP/m3/s) flux of labile P to secondary mineral P 
     real(r8), pointer :: secondp_to_labilep_vr_col                   (:,:)   ! col (gP/m3/s) flux of the desorption of secondary mineral P to labile P
     real(r8), pointer :: secondp_to_labilep_col                      (:)     ! col (gP/m3/s) flux of the desorption of secondary mineral P to labile P
     real(r8), pointer :: secondp_to_occlp_vr_col                     (:,:)   ! col (gP/m3/s) flux of the occlusion of secondary P to occluded P
     real(r8), pointer :: secondp_to_occlp_col                        (:)     ! col (gP/m3/s) flux of the occlusion of secondary P to occluded P


     ! leaching fluxes
     real(r8), pointer :: sminp_leached_vr_col                      (:,:)   ! col vertically-resolved soil mineral P pool loss to leaching (gP/m3/s)
     real(r8), pointer :: sminp_leached_col                         (:)     ! col soil mineral P pool loss to leaching (gP/m2/s)

     ! crop fluxes
     real(r8), pointer :: crop_seedp_to_leaf_patch                  (:)     ! patch (gP/m2/s) seed source to leaf, for crops

     ! dynamic landcover fluxes
     real(r8), pointer :: dwt_seedp_to_leaf_patch                   (:)     ! (gP/m2/s) seed source to patch-level; although this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_seedp_to_leaf_grc                     (:)     ! (gP/m2/s) dwt_seedn_to_leaf_patch summed to the gridcell-level
     real(r8), pointer :: dwt_seedp_to_deadstem_patch               (:)     ! (gP/m2/s) seed source to patch-level; although this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_seedp_to_deadstem_grc                 (:)     ! (gP/m2/s) dwt_seedn_to_deadstem_patch summed to the gridcell-level
     real(r8), pointer :: dwt_conv_pflux_patch                      (:)     ! (gP/m2/s) conversion N flux (immediate loss to atm); although this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_conv_pflux_grc                        (:)     ! (gP/m2/s) dwt_conv_nflux_patch summed to the gridcell-level
     real(r8), pointer :: dwt_prod10p_gain_patch                    (:)     ! patch (gP/m2/s) addition to 10-yr wood product pool; even though this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_prod100p_gain_patch                   (:)     ! patch (gP/m2/s) addition to 100-yr wood product pool; even though this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_crop_productp_gain_patch              (:)     ! patch (gP/m2/s) addition to crop product pool from landcover change; even though this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_slash_pflux_col                       (:)     ! (gP/m2/s) conversion slash flux due to landcover change

     real(r8), pointer :: dwt_seedp_to_ppool_grc                    (:)     ! col (gP/m2/s) seed source to PFT-level
     real(r8), pointer :: dwt_seedp_to_ppool_patch                  (:)     ! col (gP/m2/s) seed source to PFT-level

     real(r8), pointer :: dwt_conv_pflux_col                        (:)     ! col (gP/m2/s) conversion P flux (immediate loss to atm)
     real(r8), pointer :: dwt_prod10p_gain_col                      (:)     ! col (gP/m2/s) addition to 10-yr wood product pool
     real(r8), pointer :: dwt_prod100p_gain_col                     (:)     ! col (gP/m2/s) addition to 100-yr wood product pool
     real(r8), pointer :: dwt_frootp_to_litr_met_p_col              (:,:)   ! col (gP/m3/s) fine root to litter due to landcover change
     real(r8), pointer :: dwt_frootp_to_litr_cel_p_col              (:,:)   ! col (gP/m3/s) fine root to litter due to landcover change
     real(r8), pointer :: dwt_frootp_to_litr_lig_p_col              (:,:)   ! col (gP/m3/s) fine root to litter due to landcover change
     real(r8), pointer :: dwt_livecrootp_to_cwdp_col                (:,:)   ! col (gP/m3/s) live coarse root to CWD due to landcover change
     real(r8), pointer :: dwt_deadcrootp_to_cwdp_col                (:,:)   ! col (gP/m3/s) dead coarse root to CWD due to landcover change
     real(r8), pointer :: dwt_ploss_col                             (:)     ! col (gP/m2/s) total phosphorus loss from product pools and conversion
     real(r8), pointer :: dwt_prod10p_gain_grc                      (:)     ! col (gP/m2/s) addition to 10-yr wood product pool
     real(r8), pointer :: dwt_prod100p_gain_grc                     (:)     ! col (gP/m2/s) addition to 100-yr wood product pool

     ! wood product pool loss fluxes
     real(r8), pointer :: prod1p_loss_col                           (:)     ! col (gP/m2/s) decomposition loss from 1-yr crop product pool
     real(r8), pointer :: prod10p_loss_col                          (:)     ! col (gP/m2/s) decomposition loss from 10-yr wood product pool
     real(r8), pointer :: prod100p_loss_col                         (:)     ! col (gP/m2/s) decomposition loss from 100-yr wood product pool
     real(r8), pointer :: product_ploss_col                         (:)     ! col (gP/m2/s) total wood product phosphorus loss


     ! summary (diagnostic) flux variables, not involved in mass balance
     real(r8), pointer :: pinputs_col                               (:)     ! col column-level P inputs (gP/m2/s)
     real(r8), pointer :: poutputs_col                              (:)     ! col column-level P outputs (gP/m2/s)
     real(r8), pointer :: som_p_leached_col                         (:)     ! col total SOM P loss from vertical transport (gP/m^2/s)
     real(r8), pointer :: decomp_ppools_leached_col                 (:,:)   ! col P loss from vertical transport from each decomposing P pool (gP/m^2/s)
     real(r8), pointer :: decomp_ppools_transport_tendency_col      (:,:,:) ! col P tendency due to vertical transport in decomposing P pools (gP/m^3/s)

     ! all n pools involved in decomposition
     real(r8), pointer :: decomp_ppools_sourcesink_col              (:,:,:) ! col (gP/m3) change in decomposing P pools
                                                                            !     (sum of all additions and subtractions from stateupdate1).  

     ! Misc
     real(r8), pointer :: plant_pdemand_patch                       (:)     ! P flux required to support initial GPP (gP/m2/s)
     real(r8), pointer :: avail_retransp_patch                      (:)     ! P flux available from retranslocation pool (gP/m2/s)
     real(r8), pointer :: plant_palloc_patch                        (:)     ! total allocated P flux (gP/m2/s)

     ! clm_interface & pflotran
     !------------------------------------------------------------------------
     real(r8), pointer :: plant_pdemand_col                         (:)     ! col P flux required to support initial GPP (gN/m2/s)
     real(r8), pointer :: plant_pdemand_vr_col                      (:,:)   ! col vertically-resolved P flux required to support initial GPP (gP/m3/s)
     ! for PF-bgc mass-balance error checking
     real(r8), pointer :: externalp_to_decomp_ppools_col            (:,:,:) ! col net N fluxes associated with litter/som-adding/removal to decomp pools (gP/m3/s)
                                                                            ! (sum of all external P additions and removals, excluding decomposition/hr).
     real(r8), pointer :: externalp_to_decomp_delta_col             (:)     ! col summarized net N i/o changes associated with litter/som-adding/removal to decomp pools  btw time-step (gP/m2)
     real(r8), pointer :: sminp_net_transport_vr_col                (:,:)   ! col net sminp transport associated with runoff/leaching (gP/m3/s)
     real(r8), pointer :: sminp_net_transport_delta_col             (:)     ! col summarized net change of column-level sminp leaching bwtn time-step (for balance checking) (gP/m2)
     !------------------------------------------------------------------------

     real(r8), pointer :: sminp_to_plant_patch                      (:)     ! pft-level plant p uptake (gP/m2/s)
     real(r8), pointer :: plant_pdemand_vr_patch                    (:,:)   ! pft-level plant P demand
     real(r8), pointer :: prev_leafp_to_litter_patch                (:)     ! previous timestep leaf P litterfall flux (gP/m2/s)
     real(r8), pointer :: prev_frootp_to_litter_patch               (:)     ! previous timestep froot P litterfall flux (gP/m2/s)
     real(r8), pointer :: adsorb_to_labilep_vr                      (:,:)
     real(r8), pointer :: desorb_to_solutionp_vr                    (:,:)
     real(r8), pointer :: adsorb_to_labilep_col                     (:)
     real(r8), pointer :: desorb_to_solutionp_col                   (:)
     real(r8), pointer :: pmpf_decomp_cascade                       (:,:,:)

     real(r8), pointer :: plant_p_uptake_flux                       (:)     ! for the purpose of mass balance check  
     real(r8), pointer :: soil_p_immob_flux                         (:)     ! for the purpose of mass balance check
     real(r8), pointer :: soil_p_immob_flux_vr                      (:,:)   ! for the purpose of mass balance check
     real(r8), pointer :: soil_p_grossmin_flux                      (:)     ! for the purpose of mass balance check
     real(r8), pointer :: smin_p_to_plant_col                       (:)     ! for the purpose of mass balance check
     real(r8), pointer :: plant_to_litter_pflux                     (:)     ! for the purpose of mass balance check
     real(r8), pointer :: plant_to_cwd_pflux                        (:)     ! for the purpose of mass balance check
     real(r8), pointer :: supplement_to_plantp                      (:)     ! supplementary P flux for plant 

     real(r8), pointer :: gap_ploss_litter                          (:)     ! total ploss from veg to litter pool due to gap mortality
     real(r8), pointer :: fire_ploss_litter                         (:)     ! total ploss from veg to litter pool due to fire
     real(r8), pointer :: hrv_ploss_litter                          (:)     ! total ploss from veg to litter pool due to harvest mortality
     real(r8), pointer :: sen_ploss_litter                          (:)     ! total ploss from veg to litter pool due to senescence

     ! C4MIP output variable
     real(r8), pointer :: plant_p_to_cwdp                           (:) ! sum of gap, fire, dynamic land use, and harvest mortality, plant phosphorus flux to CWD

   contains

     procedure , public  :: Init   
     procedure , public  :: Restart
     procedure , public  :: SetValues
     procedure , public  :: ZeroDWT
     procedure , public  :: Summary
     procedure , private :: InitAllocate
     procedure , private :: InitHistory
     procedure , private :: InitCold
     ! bgc & pflotran interface

  end type phosphorusflux_type
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(phosphorusflux_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate (bounds)
    call this%InitHistory (bounds)
    call this%InitCold (bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize pft phosphorus flux
    !
    ! !ARGUMENTS:
    class (phosphorusflux_type) :: this
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
    !------------------------------------------------------------------------
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
    class(phosphorusflux_type) :: this
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
    character(1)   :: aa 
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
    ! P flux variables - native to column
    !-------------------------------


    ! bgc interface
  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initializes time varying variables used only in coupled carbon-phosphorus mode (CN):
    !
    ! !USES:
    use elm_varpar      , only : crop_prog
    use landunit_varcon , only : istsoil, istcrop
    !
    ! !ARGUMENTS:
    class(phosphorusflux_type) :: this 
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: p,c,l,g
    integer :: fp, fc                                    ! filter indices
    integer :: num_special_col                           ! number of good values in special_col filter
    integer :: num_special_patch                         ! number of good values in special_patch filter
    integer :: special_col(bounds%endc-bounds%begc+1)    ! special landunit filter - columns
    integer :: special_patch(bounds%endp-bounds%begp+1)  ! special landunit filter - patches
    !---------------------------------------------------------------------

    ! Set column filters

    ! Set patch filters


    !-----------------------------------------------
    ! initialize phosphorus flux variables
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
    !
    ! !ARGUMENTS:
    class (phosphorusflux_type) :: this
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

    ! clm_interface & pflotran
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine SetValues ( this, &
       num_patch, filter_patch, value_patch, &
       num_column, filter_column, value_column)
    !
    ! !DESCRIPTION:
    ! Set phosphorus flux variables
    !
    ! !ARGUMENTS:
    ! !ARGUMENTS:
    class (phosphorusflux_type) :: this
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

  end subroutine SetValues

  !-----------------------------------------------------------------------
  subroutine ZeroDwt( this, bounds )
    !
    ! !DESCRIPTION
    ! Initialize flux variables needed for dynamic land use.
    !
    ! !ARGUMENTS:
    class(phosphorusflux_type) :: this
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
    use elm_varctl    , only: use_nitrif_denitrif
    use subgridAveMod , only: p2c
    use pftvarcon     , only: npcropmin
    ! pflotran
!    use elm_varctl    , only: use_pflotran, pf_cmode
    !
    ! !ARGUMENTS:
    class (phosphorusflux_type) :: this
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


  end subroutine Summary

!-------------------------------------------------------------------------------------------------

end module PhosphorusFluxType

