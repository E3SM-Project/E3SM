module PhosphorusFluxType

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
     real(r8), pointer :: biochem_pmin_col                          (:)   ! col vert-int (diagnostic) total biochemical P mineralization (gP/m3/s)


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

     ! dynamic landcover fluxes
     real(r8), pointer :: dwt_seedp_to_leaf_col                     (:)     ! col (gP/m2/s) seed source to PFT-level
     real(r8), pointer :: dwt_seedp_to_deadstem_col                 (:)     ! col (gP/m2/s) seed source to PFT-level
     real(r8), pointer :: dwt_conv_pflux_col                        (:)     ! col (gP/m2/s) conversion P flux (immediate loss to atm)
     real(r8), pointer :: dwt_prod10p_gain_col                      (:)     ! col (gP/m2/s) addition to 10-yr wood product pool
     real(r8), pointer :: dwt_prod100p_gain_col                     (:)     ! col (gP/m2/s) addition to 100-yr wood product pool
     real(r8), pointer :: dwt_frootp_to_litr_met_p_col              (:,:)   ! col (gP/m3/s) fine root to litter due to landcover change
     real(r8), pointer :: dwt_frootp_to_litr_cel_p_col              (:,:)   ! col (gP/m3/s) fine root to litter due to landcover change
     real(r8), pointer :: dwt_frootp_to_litr_lig_p_col              (:,:)   ! col (gP/m3/s) fine root to litter due to landcover change
     real(r8), pointer :: dwt_livecrootp_to_cwdp_col                (:,:)   ! col (gP/m3/s) live coarse root to CWD due to landcover change
     real(r8), pointer :: dwt_deadcrootp_to_cwdp_col                (:,:)   ! col (gP/m3/s) dead coarse root to CWD due to landcover change
     real(r8), pointer :: dwt_ploss_col                             (:)     ! col (gP/m2/s) total phosphorus loss from product pools and conversion

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

     ! clm_bgc_interface & pflotran
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

   contains

     procedure , public  :: Init
     procedure , public  :: Restart
     procedure , public  :: SetValues
     procedure , public  :: ZeroDWT
     procedure , public  :: Summary
     procedure , private :: InitAllocate
     procedure , private :: InitHistory
     procedure , private :: InitCold
     procedure , private  :: Summary_betr
     !! bgc & pflotran interface
     procedure , private :: PSummary_interface

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
    !------------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc

    allocate(this%m_leafp_to_litter_patch                   (begp:endp)) ; this%m_leafp_to_litter_patch                   (:) = nan
    allocate(this%m_frootp_to_litter_patch                  (begp:endp)) ; this%m_frootp_to_litter_patch                  (:) = nan
    allocate(this%m_leafp_storage_to_litter_patch           (begp:endp)) ; this%m_leafp_storage_to_litter_patch           (:) = nan
    allocate(this%m_frootp_storage_to_litter_patch          (begp:endp)) ; this%m_frootp_storage_to_litter_patch          (:) = nan
    allocate(this%m_livestemp_storage_to_litter_patch       (begp:endp)) ; this%m_livestemp_storage_to_litter_patch       (:) = nan
    allocate(this%m_deadstemp_storage_to_litter_patch       (begp:endp)) ; this%m_deadstemp_storage_to_litter_patch       (:) = nan
    allocate(this%m_livecrootp_storage_to_litter_patch      (begp:endp)) ; this%m_livecrootp_storage_to_litter_patch      (:) = nan
    allocate(this%m_deadcrootp_storage_to_litter_patch      (begp:endp)) ; this%m_deadcrootp_storage_to_litter_patch      (:) = nan
    allocate(this%m_leafp_xfer_to_litter_patch              (begp:endp)) ; this%m_leafp_xfer_to_litter_patch              (:) = nan
    allocate(this%m_frootp_xfer_to_litter_patch             (begp:endp)) ; this%m_frootp_xfer_to_litter_patch             (:) = nan
    allocate(this%m_livestemp_xfer_to_litter_patch          (begp:endp)) ; this%m_livestemp_xfer_to_litter_patch          (:) = nan
    allocate(this%m_deadstemp_xfer_to_litter_patch          (begp:endp)) ; this%m_deadstemp_xfer_to_litter_patch          (:) = nan
    allocate(this%m_livecrootp_xfer_to_litter_patch         (begp:endp)) ; this%m_livecrootp_xfer_to_litter_patch         (:) = nan
    allocate(this%m_deadcrootp_xfer_to_litter_patch         (begp:endp)) ; this%m_deadcrootp_xfer_to_litter_patch         (:) = nan
    allocate(this%m_livestemp_to_litter_patch               (begp:endp)) ; this%m_livestemp_to_litter_patch               (:) = nan
    allocate(this%m_deadstemp_to_litter_patch               (begp:endp)) ; this%m_deadstemp_to_litter_patch               (:) = nan
    allocate(this%m_livecrootp_to_litter_patch              (begp:endp)) ; this%m_livecrootp_to_litter_patch              (:) = nan
    allocate(this%m_deadcrootp_to_litter_patch              (begp:endp)) ; this%m_deadcrootp_to_litter_patch              (:) = nan
    allocate(this%m_retransp_to_litter_patch                (begp:endp)) ; this%m_retransp_to_litter_patch                (:) = nan
    allocate(this%hrv_leafp_to_litter_patch                 (begp:endp)) ; this%hrv_leafp_to_litter_patch                 (:) = nan
    allocate(this%hrv_frootp_to_litter_patch                (begp:endp)) ; this%hrv_frootp_to_litter_patch                (:) = nan
    allocate(this%hrv_leafp_storage_to_litter_patch         (begp:endp)) ; this%hrv_leafp_storage_to_litter_patch         (:) = nan
    allocate(this%hrv_frootp_storage_to_litter_patch        (begp:endp)) ; this%hrv_frootp_storage_to_litter_patch        (:) = nan
    allocate(this%hrv_livestemp_storage_to_litter_patch     (begp:endp)) ; this%hrv_livestemp_storage_to_litter_patch     (:) = nan
    allocate(this%hrv_deadstemp_storage_to_litter_patch     (begp:endp)) ; this%hrv_deadstemp_storage_to_litter_patch     (:) = nan
    allocate(this%hrv_livecrootp_storage_to_litter_patch    (begp:endp)) ; this%hrv_livecrootp_storage_to_litter_patch    (:) = nan
    allocate(this%hrv_deadcrootp_storage_to_litter_patch    (begp:endp)) ; this%hrv_deadcrootp_storage_to_litter_patch    (:) = nan
    allocate(this%hrv_leafp_xfer_to_litter_patch            (begp:endp)) ; this%hrv_leafp_xfer_to_litter_patch            (:) = nan
    allocate(this%hrv_frootp_xfer_to_litter_patch           (begp:endp)) ; this%hrv_frootp_xfer_to_litter_patch           (:) = nan
    allocate(this%hrv_livestemp_xfer_to_litter_patch        (begp:endp)) ; this%hrv_livestemp_xfer_to_litter_patch        (:) = nan
    allocate(this%hrv_deadstemp_xfer_to_litter_patch        (begp:endp)) ; this%hrv_deadstemp_xfer_to_litter_patch        (:) = nan
    allocate(this%hrv_livecrootp_xfer_to_litter_patch       (begp:endp)) ; this%hrv_livecrootp_xfer_to_litter_patch       (:) = nan
    allocate(this%hrv_deadcrootp_xfer_to_litter_patch       (begp:endp)) ; this%hrv_deadcrootp_xfer_to_litter_patch       (:) = nan
    allocate(this%hrv_livestemp_to_litter_patch             (begp:endp)) ; this%hrv_livestemp_to_litter_patch             (:) = nan
    allocate(this%hrv_deadstemp_to_prod10p_patch            (begp:endp)) ; this%hrv_deadstemp_to_prod10p_patch            (:) = nan
    allocate(this%hrv_deadstemp_to_prod100p_patch           (begp:endp)) ; this%hrv_deadstemp_to_prod100p_patch           (:) = nan
    allocate(this%hrv_leafp_to_prod1p_patch                 (begp:endp)) ; this%hrv_leafp_to_prod1p_patch                 (:) = nan
    allocate(this%hrv_livestemp_to_prod1p_patch             (begp:endp)) ; this%hrv_livestemp_to_prod1p_patch             (:) = nan
    allocate(this%hrv_grainp_to_prod1p_patch                (begp:endp)) ; this%hrv_grainp_to_prod1p_patch                (:) = nan
    allocate(this%hrv_cropp_to_prod1p_patch                 (begp:endp)) ; this%hrv_cropp_to_prod1p_patch                 (:) = nan
    allocate(this%hrv_livecrootp_to_litter_patch            (begp:endp)) ; this%hrv_livecrootp_to_litter_patch            (:) = nan
    allocate(this%hrv_deadcrootp_to_litter_patch            (begp:endp)) ; this%hrv_deadcrootp_to_litter_patch            (:) = nan
    allocate(this%hrv_retransp_to_litter_patch              (begp:endp)) ; this%hrv_retransp_to_litter_patch              (:) = nan

    allocate(this%m_leafp_to_fire_patch                     (begp:endp)) ; this%m_leafp_to_fire_patch                     (:) = nan
    allocate(this%m_leafp_storage_to_fire_patch             (begp:endp)) ; this%m_leafp_storage_to_fire_patch             (:) = nan
    allocate(this%m_leafp_xfer_to_fire_patch                (begp:endp)) ; this%m_leafp_xfer_to_fire_patch                (:) = nan
    allocate(this%m_livestemp_to_fire_patch                 (begp:endp)) ; this%m_livestemp_to_fire_patch                 (:) = nan
    allocate(this%m_livestemp_storage_to_fire_patch         (begp:endp)) ; this%m_livestemp_storage_to_fire_patch         (:) = nan
    allocate(this%m_livestemp_xfer_to_fire_patch            (begp:endp)) ; this%m_livestemp_xfer_to_fire_patch            (:) = nan
    allocate(this%m_deadstemp_to_fire_patch                 (begp:endp)) ; this%m_deadstemp_to_fire_patch                 (:) = nan
    allocate(this%m_deadstemp_storage_to_fire_patch         (begp:endp)) ; this%m_deadstemp_storage_to_fire_patch         (:) = nan
    allocate(this%m_deadstemp_xfer_to_fire_patch            (begp:endp)) ; this%m_deadstemp_xfer_to_fire_patch            (:) = nan
    allocate(this%m_frootp_to_fire_patch                    (begp:endp)) ; this%m_frootp_to_fire_patch                    (:) = nan
    allocate(this%m_frootp_storage_to_fire_patch            (begp:endp)) ; this%m_frootp_storage_to_fire_patch            (:) = nan
    allocate(this%m_frootp_xfer_to_fire_patch               (begp:endp)) ; this%m_frootp_xfer_to_fire_patch               (:) = nan
    allocate(this%m_livecrootp_to_fire_patch                (begp:endp)) ;
    allocate(this%m_livecrootp_storage_to_fire_patch        (begp:endp)) ; this%m_livecrootp_storage_to_fire_patch        (:) = nan
    allocate(this%m_livecrootp_xfer_to_fire_patch           (begp:endp)) ; this%m_livecrootp_xfer_to_fire_patch           (:) = nan
    allocate(this%m_deadcrootp_to_fire_patch                (begp:endp)) ; this%m_deadcrootp_to_fire_patch                (:) = nan
    allocate(this%m_deadcrootp_storage_to_fire_patch        (begp:endp)) ; this%m_deadcrootp_storage_to_fire_patch        (:) = nan
    allocate(this%m_deadcrootp_xfer_to_fire_patch           (begp:endp)) ; this%m_deadcrootp_xfer_to_fire_patch           (:) = nan
    allocate(this%m_retransp_to_fire_patch                  (begp:endp)) ; this%m_retransp_to_fire_patch                  (:) = nan

    allocate(this%m_leafp_to_litter_fire_patch              (begp:endp)) ; this%m_leafp_to_litter_fire_patch              (:) = nan
    allocate(this%m_leafp_storage_to_litter_fire_patch      (begp:endp)) ; this%m_leafp_storage_to_litter_fire_patch      (:) = nan
    allocate(this%m_leafp_xfer_to_litter_fire_patch         (begp:endp)) ; this%m_leafp_xfer_to_litter_fire_patch         (:) = nan
    allocate(this%m_livestemp_to_litter_fire_patch          (begp:endp)) ; this%m_livestemp_to_litter_fire_patch          (:) = nan
    allocate(this%m_livestemp_storage_to_litter_fire_patch  (begp:endp)) ; this%m_livestemp_storage_to_litter_fire_patch  (:) = nan
    allocate(this%m_livestemp_xfer_to_litter_fire_patch     (begp:endp)) ; this%m_livestemp_xfer_to_litter_fire_patch     (:) = nan
    allocate(this%m_livestemp_to_deadstemp_fire_patch       (begp:endp)) ; this%m_livestemp_to_deadstemp_fire_patch       (:) = nan
    allocate(this%m_deadstemp_to_litter_fire_patch          (begp:endp)) ; this%m_deadstemp_to_litter_fire_patch          (:) = nan
    allocate(this%m_deadstemp_storage_to_litter_fire_patch  (begp:endp)) ; this%m_deadstemp_storage_to_litter_fire_patch  (:) = nan
    allocate(this%m_deadstemp_xfer_to_litter_fire_patch     (begp:endp)) ; this%m_deadstemp_xfer_to_litter_fire_patch     (:) = nan
    allocate(this%m_frootp_to_litter_fire_patch             (begp:endp)) ; this%m_frootp_to_litter_fire_patch             (:) = nan
    allocate(this%m_frootp_storage_to_litter_fire_patch     (begp:endp)) ; this%m_frootp_storage_to_litter_fire_patch     (:) = nan
    allocate(this%m_frootp_xfer_to_litter_fire_patch        (begp:endp)) ; this%m_frootp_xfer_to_litter_fire_patch        (:) = nan
    allocate(this%m_livecrootp_to_litter_fire_patch         (begp:endp)) ; this%m_livecrootp_to_litter_fire_patch         (:) = nan
    allocate(this%m_livecrootp_storage_to_litter_fire_patch (begp:endp)) ; this%m_livecrootp_storage_to_litter_fire_patch (:) = nan
    allocate(this%m_livecrootp_xfer_to_litter_fire_patch    (begp:endp)) ; this%m_livecrootp_xfer_to_litter_fire_patch    (:) = nan
    allocate(this%m_livecrootp_to_deadcrootp_fire_patch     (begp:endp)) ; this%m_livecrootp_to_deadcrootp_fire_patch     (:) = nan
    allocate(this%m_deadcrootp_to_litter_fire_patch         (begp:endp)) ; this%m_deadcrootp_to_litter_fire_patch         (:) = nan
    allocate(this%m_deadcrootp_storage_to_litter_fire_patch (begp:endp)) ; this%m_deadcrootp_storage_to_litter_fire_patch (:) = nan
    allocate(this%m_deadcrootp_xfer_to_litter_fire_patch    (begp:endp)) ; this%m_deadcrootp_xfer_to_litter_fire_patch    (:) = nan
    allocate(this%m_retransp_to_litter_fire_patch           (begp:endp)) ; this%m_retransp_to_litter_fire_patch           (:) = nan

    allocate(this%leafp_xfer_to_leafp_patch                 (begp:endp)) ; this%leafp_xfer_to_leafp_patch                 (:) = nan
    allocate(this%frootp_xfer_to_frootp_patch               (begp:endp)) ; this%frootp_xfer_to_frootp_patch               (:) = nan
    allocate(this%livestemp_xfer_to_livestemp_patch         (begp:endp)) ; this%livestemp_xfer_to_livestemp_patch         (:) = nan
    allocate(this%deadstemp_xfer_to_deadstemp_patch         (begp:endp)) ; this%deadstemp_xfer_to_deadstemp_patch         (:) = nan
    allocate(this%livecrootp_xfer_to_livecrootp_patch       (begp:endp)) ; this%livecrootp_xfer_to_livecrootp_patch       (:) = nan
    allocate(this%deadcrootp_xfer_to_deadcrootp_patch       (begp:endp)) ; this%deadcrootp_xfer_to_deadcrootp_patch       (:) = nan
    allocate(this%leafp_to_litter_patch                     (begp:endp)) ; this%leafp_to_litter_patch                     (:) = nan
    allocate(this%leafp_to_retransp_patch                   (begp:endp)) ; this%leafp_to_retransp_patch                   (:) = nan
    allocate(this%frootp_to_retransp_patch                  (begp:endp)) ; this%frootp_to_retransp_patch                  (:) = nan
    allocate(this%frootp_to_litter_patch                    (begp:endp)) ; this%frootp_to_litter_patch                    (:) = nan
    allocate(this%retransp_to_ppool_patch                   (begp:endp)) ; this%retransp_to_ppool_patch                   (:) = nan
    allocate(this%sminp_to_ppool_patch                      (begp:endp)) ; this%sminp_to_ppool_patch                      (:) = nan

    allocate(this%ppool_to_leafp_patch              (begp:endp)) ; this%ppool_to_leafp_patch              (:) = nan
    allocate(this%ppool_to_leafp_storage_patch      (begp:endp)) ; this%ppool_to_leafp_storage_patch      (:) = nan
    allocate(this%ppool_to_frootp_patch             (begp:endp)) ; this%ppool_to_frootp_patch             (:) = nan
    allocate(this%ppool_to_frootp_storage_patch     (begp:endp)) ; this%ppool_to_frootp_storage_patch     (:) = nan
    allocate(this%ppool_to_livestemp_patch          (begp:endp)) ; this%ppool_to_livestemp_patch          (:) = nan
    allocate(this%ppool_to_livestemp_storage_patch  (begp:endp)) ; this%ppool_to_livestemp_storage_patch  (:) = nan
    allocate(this%ppool_to_deadstemp_patch          (begp:endp)) ; this%ppool_to_deadstemp_patch          (:) = nan
    allocate(this%ppool_to_deadstemp_storage_patch  (begp:endp)) ; this%ppool_to_deadstemp_storage_patch  (:) = nan
    allocate(this%ppool_to_livecrootp_patch         (begp:endp)) ; this%ppool_to_livecrootp_patch         (:) = nan
    allocate(this%ppool_to_livecrootp_storage_patch (begp:endp)) ; this%ppool_to_livecrootp_storage_patch (:) = nan
    allocate(this%ppool_to_deadcrootp_patch         (begp:endp)) ; this%ppool_to_deadcrootp_patch         (:) = nan
    allocate(this%ppool_to_deadcrootp_storage_patch (begp:endp)) ; this%ppool_to_deadcrootp_storage_patch (:) = nan
    allocate(this%leafp_storage_to_xfer_patch       (begp:endp)) ; this%leafp_storage_to_xfer_patch       (:) = nan
    allocate(this%frootp_storage_to_xfer_patch      (begp:endp)) ; this%frootp_storage_to_xfer_patch      (:) = nan
    allocate(this%livestemp_storage_to_xfer_patch   (begp:endp)) ; this%livestemp_storage_to_xfer_patch   (:) = nan
    allocate(this%deadstemp_storage_to_xfer_patch   (begp:endp)) ; this%deadstemp_storage_to_xfer_patch   (:) = nan
    allocate(this%livecrootp_storage_to_xfer_patch  (begp:endp)) ; this%livecrootp_storage_to_xfer_patch  (:) = nan
    allocate(this%deadcrootp_storage_to_xfer_patch  (begp:endp)) ; this%deadcrootp_storage_to_xfer_patch  (:) = nan
    allocate(this%livestemp_to_deadstemp_patch      (begp:endp)) ; this%livestemp_to_deadstemp_patch      (:) = nan
    allocate(this%livestemp_to_retransp_patch       (begp:endp)) ; this%livestemp_to_retransp_patch       (:) = nan
    allocate(this%livecrootp_to_deadcrootp_patch    (begp:endp)) ; this%livecrootp_to_deadcrootp_patch    (:) = nan
    allocate(this%livecrootp_to_retransp_patch      (begp:endp)) ; this%livecrootp_to_retransp_patch      (:) = nan
    allocate(this%pdeploy_patch                     (begp:endp)) ; this%pdeploy_patch                     (:) = nan
    allocate(this%pinputs_patch                     (begp:endp)) ; this%pinputs_patch                     (:) = nan
    allocate(this%poutputs_patch                    (begp:endp)) ; this%poutputs_patch                    (:) = nan
    allocate(this%wood_harvestp_patch               (begp:endp)) ; this%wood_harvestp_patch               (:) = nan
    allocate(this%fire_ploss_patch                  (begp:endp)) ; this%fire_ploss_patch                  (:) = nan
    allocate(this%ppool_to_grainp_patch             (begp:endp)) ; this%ppool_to_grainp_patch             (:) = nan
    allocate(this%ppool_to_grainp_storage_patch     (begp:endp)) ; this%ppool_to_grainp_storage_patch     (:) = nan
    allocate(this%livestemp_to_litter_patch         (begp:endp)) ; this%livestemp_to_litter_patch         (:) = nan
    allocate(this%grainp_to_food_patch              (begp:endp)) ; this%grainp_to_food_patch              (:) = nan
    allocate(this%grainp_xfer_to_grainp_patch       (begp:endp)) ; this%grainp_xfer_to_grainp_patch       (:) = nan
    allocate(this%grainp_storage_to_xfer_patch      (begp:endp)) ; this%grainp_storage_to_xfer_patch      (:) = nan
    allocate(this%fert_p_patch                      (begp:endp)) ; this%fert_p_patch                      (:) = nan
    allocate(this%fert_p_counter_patch                (begp:endp)) ; this%fert_p_counter_patch                (:) = nan

    allocate(this%pdep_to_sminp_col             (begc:endc))    ; this%pdep_to_sminp_col     (:) = nan
    allocate(this%fert_p_to_sminp_col             (begc:endc))    ; this%fert_p_to_sminp_col     (:) = nan
    allocate(this%hrv_deadstemp_to_prod10p_col  (begc:endc))    ; this%hrv_deadstemp_to_prod10p_col  (:) = nan
    allocate(this%hrv_deadstemp_to_prod100p_col (begc:endc))    ; this%hrv_deadstemp_to_prod100p_col (:) = nan
    allocate(this%hrv_cropp_to_prod1p_col       (begc:endc))    ; this%hrv_cropp_to_prod1p_col       (:) = nan
    allocate(this%sminp_to_plant_col            (begc:endc))    ; this%sminp_to_plant_col     (:) = nan
    allocate(this%potential_immob_p_col           (begc:endc))    ; this%potential_immob_p_col           (:) = nan
    allocate(this%actual_immob_p_col              (begc:endc))    ; this%actual_immob_p_col              (:) = nan
    allocate(this%gross_pmin_col                (begc:endc))    ; this%gross_pmin_col                (:) = nan
    allocate(this%net_pmin_col                  (begc:endc))    ; this%net_pmin_col                  (:) = nan
    allocate(this%supplement_to_sminp_col       (begc:endc))    ; this%supplement_to_sminp_col       (:) = nan
    allocate(this%prod1p_loss_col               (begc:endc))    ; this%prod1p_loss_col              (:) = nan
    allocate(this%prod10p_loss_col              (begc:endc))    ; this%prod10p_loss_col              (:) = nan
    allocate(this%prod100p_loss_col             (begc:endc))    ; this%prod100p_loss_col     (:) = nan
    allocate(this%product_ploss_col             (begc:endc))    ; this%product_ploss_col     (:) = nan
    allocate(this%pinputs_col                   (begc:endc))    ; this%pinputs_col                   (:) = nan
    allocate(this%poutputs_col                  (begc:endc))    ; this%poutputs_col                  (:) = nan
    allocate(this%fire_ploss_col                (begc:endc))    ; this%fire_ploss_col                (:) = nan
    allocate(this%fire_decomp_ploss_col         (begc:endc))    ; this%fire_decomp_ploss_col         (:) = nan
    allocate(this%fire_ploss_p2c_col            (begc:endc))    ; this%fire_ploss_p2c_col            (:) = nan
    allocate(this%som_p_leached_col             (begc:endc))    ; this%som_p_leached_col     (:) = nan

    allocate(this%m_p_to_litr_met_fire_col   (begc:endc,1:nlevdecomp_full)) ; this%m_p_to_litr_met_fire_col   (:,:) = nan
    allocate(this%m_p_to_litr_cel_fire_col   (begc:endc,1:nlevdecomp_full)) ; this%m_p_to_litr_cel_fire_col   (:,:) = nan
    allocate(this%m_p_to_litr_lig_fire_col   (begc:endc,1:nlevdecomp_full)) ; this%m_p_to_litr_lig_fire_col   (:,:) = nan
    allocate(this%potential_immob_p_vr_col   (begc:endc,1:nlevdecomp_full)) ; this%potential_immob_p_vr_col     (:,:) = nan
    allocate(this%actual_immob_p_vr_col      (begc:endc,1:nlevdecomp_full)) ; this%actual_immob_p_vr_col        (:,:) = nan
    allocate(this%sminp_to_plant_vr_col      (begc:endc,1:nlevdecomp_full)) ; this%sminp_to_plant_vr_col      (:,:) = nan
    allocate(this%supplement_to_sminp_vr_col (begc:endc,1:nlevdecomp_full)) ; this%supplement_to_sminp_vr_col (:,:) = nan
    allocate(this%gross_pmin_vr_col          (begc:endc,1:nlevdecomp_full)) ; this%gross_pmin_vr_col          (:,:) = nan
    allocate(this%net_pmin_vr_col            (begc:endc,1:nlevdecomp_full)) ; this%net_pmin_vr_col            (:,:) = nan

    allocate(this%biochem_pmin_ppools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
    allocate(this%biochem_pmin_vr_col       (begc:endc,1:nlevdecomp_full))
    allocate(this%biochem_pmin_col          (begc:endc))

    this%biochem_pmin_ppools_vr_col  (:,:,:) = nan
    this%biochem_pmin_vr_col         (:,:) = nan
    this%biochem_pmin_col            (:) = nan

    allocate(this%dwt_seedp_to_leaf_col      (begc:endc))                   ; this%dwt_seedp_to_leaf_col      (:)   = nan
    allocate(this%dwt_seedp_to_deadstem_col  (begc:endc))                   ; this%dwt_seedp_to_deadstem_col  (:)   = nan
    allocate(this%dwt_conv_pflux_col         (begc:endc))                   ; this%dwt_conv_pflux_col         (:)   = nan
    allocate(this%dwt_prod10p_gain_col       (begc:endc))                   ; this%dwt_prod10p_gain_col       (:)   = nan
    allocate(this%dwt_prod100p_gain_col      (begc:endc))                   ; this%dwt_prod100p_gain_col      (:)   = nan
    allocate(this%dwt_ploss_col              (begc:endc))                   ; this%dwt_ploss_col              (:)   = nan
    allocate(this%wood_harvestp_col          (begc:endc))                   ; this%wood_harvestp_col          (:)   = nan

    allocate(this%dwt_frootp_to_litr_met_p_col(begc:endc,1:nlevdecomp_full)) ; this%dwt_frootp_to_litr_met_p_col     (:,:) = nan
    allocate(this%dwt_frootp_to_litr_cel_p_col(begc:endc,1:nlevdecomp_full)) ; this%dwt_frootp_to_litr_cel_p_col     (:,:) = nan
    allocate(this%dwt_frootp_to_litr_lig_p_col(begc:endc,1:nlevdecomp_full)) ; this%dwt_frootp_to_litr_lig_p_col     (:,:) = nan
    allocate(this%dwt_livecrootp_to_cwdp_col  (begc:endc,1:nlevdecomp_full)) ; this%dwt_livecrootp_to_cwdp_col       (:,:) = nan
    allocate(this%dwt_deadcrootp_to_cwdp_col  (begc:endc,1:nlevdecomp_full)) ; this%dwt_deadcrootp_to_cwdp_col       (:,:) = nan


    allocate(this%decomp_cascade_ptransfer_vr_col   (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions ))
    allocate(this%decomp_cascade_sminp_flux_vr_col  (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions ))
    allocate(this%m_decomp_ppools_to_fire_vr_col    (begc:endc,1:nlevdecomp_full,1:ndecomp_pools               ))
    allocate(this%decomp_cascade_ptransfer_col      (begc:endc,1:ndecomp_cascade_transitions                   ))
    allocate(this%decomp_cascade_sminp_flux_col     (begc:endc,1:ndecomp_cascade_transitions                   ))
    allocate(this%m_decomp_ppools_to_fire_col       (begc:endc,1:ndecomp_pools                                 ))

    this%decomp_cascade_ptransfer_vr_col  (:,:,:) = nan
    this%decomp_cascade_sminp_flux_vr_col (:,:,:) = nan
    this%m_decomp_ppools_to_fire_vr_col   (:,:,:) = nan
    this%m_decomp_ppools_to_fire_col      (:,:)   = nan
    this%decomp_cascade_ptransfer_col     (:,:)   = nan
    this%decomp_cascade_sminp_flux_col    (:,:)   = nan

    allocate(this%phenology_p_to_litr_met_p_col     (begc:endc, 1:nlevdecomp_full))
    allocate(this%phenology_p_to_litr_cel_p_col     (begc:endc, 1:nlevdecomp_full))
    allocate(this%phenology_p_to_litr_lig_p_col     (begc:endc, 1:nlevdecomp_full))
    allocate(this%gap_mortality_p_to_litr_met_p_col (begc:endc, 1:nlevdecomp_full))
    allocate(this%gap_mortality_p_to_litr_cel_p_col (begc:endc, 1:nlevdecomp_full))
    allocate(this%gap_mortality_p_to_litr_lig_p_col (begc:endc, 1:nlevdecomp_full))
    allocate(this%gap_mortality_p_to_cwdp_col       (begc:endc, 1:nlevdecomp_full))
    allocate(this%fire_mortality_p_to_cwdp_col      (begc:endc, 1:nlevdecomp_full))
    allocate(this%harvest_p_to_litr_met_p_col       (begc:endc, 1:nlevdecomp_full))
    allocate(this%harvest_p_to_litr_cel_p_col       (begc:endc, 1:nlevdecomp_full))
    allocate(this%harvest_p_to_litr_lig_p_col       (begc:endc, 1:nlevdecomp_full))
    allocate(this%harvest_p_to_cwdp_col             (begc:endc, 1:nlevdecomp_full))

    this%phenology_p_to_litr_met_p_col     (:,:) = nan
    this%phenology_p_to_litr_cel_p_col     (:,:) = nan
    this%phenology_p_to_litr_lig_p_col     (:,:) = nan
    this%gap_mortality_p_to_litr_met_p_col (:,:) = nan
    this%gap_mortality_p_to_litr_cel_p_col (:,:) = nan
    this%gap_mortality_p_to_litr_lig_p_col (:,:) = nan
    this%gap_mortality_p_to_cwdp_col       (:,:) = nan
    this%fire_mortality_p_to_cwdp_col      (:,:) = nan
    this%harvest_p_to_litr_met_p_col       (:,:) = nan
    this%harvest_p_to_litr_cel_p_col       (:,:) = nan
    this%harvest_p_to_litr_lig_p_col       (:,:) = nan
    this%harvest_p_to_cwdp_col             (:,:) = nan



    allocate(this%primp_to_labilep_vr_col                 (begc:endc,1:nlevdecomp_full                               ))
    allocate(this%primp_to_labilep_col                    (begc:endc                                                 ))
    allocate(this%labilep_to_secondp_vr_col                 (begc:endc,1:nlevdecomp_full                               ))
    allocate(this%labilep_to_secondp_col                    (begc:endc                                                 ))
    allocate(this%secondp_to_labilep_vr_col                 (begc:endc,1:nlevdecomp_full                               ))
    allocate(this%secondp_to_labilep_col                    (begc:endc                                                 ))
    allocate(this%secondp_to_occlp_vr_col                 (begc:endc,1:nlevdecomp_full                               ))
    allocate(this%secondp_to_occlp_col                    (begc:endc                                                 ))

    this%primp_to_labilep_vr_col                 (:,:)   = nan
    this%primp_to_labilep_col                    (:)     = nan
    this%labilep_to_secondp_vr_col                 (:,:)   = nan
    this%labilep_to_secondp_col                    (:)     = nan
    this%secondp_to_labilep_vr_col                 (:,:)   = nan
    this%secondp_to_labilep_col                    (:)     = nan
    this%secondp_to_occlp_vr_col                 (:,:)   = nan
    this%secondp_to_occlp_col                    (:)     = nan

    allocate(this%sminp_leached_vr_col                 (begc:endc,1:nlevdecomp_full                               ))
    allocate(this%sminp_leached_col                    (begc:endc                                                 ))
    allocate(this%decomp_ppools_leached_col            (begc:endc,1:ndecomp_pools                                 ))
    allocate(this%decomp_ppools_transport_tendency_col (begc:endc,1:nlevdecomp_full,1:ndecomp_pools               ))

    this%sminp_leached_vr_col                 (:,:)   = nan
    this%sminp_leached_col                    (:)     = nan
    this%decomp_ppools_leached_col            (:,:)   = nan
    this%decomp_ppools_transport_tendency_col (:,:,:) = nan

    allocate(this%decomp_ppools_sourcesink_col (begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
    this%decomp_ppools_sourcesink_col (:,:,:) = nan

    allocate(this%plant_pdemand_patch         (begp:endp)) ;    this%plant_pdemand_patch         (:) = nan
    allocate(this%avail_retransp_patch        (begp:endp)) ;    this%avail_retransp_patch        (:) = nan
    allocate(this%plant_palloc_patch          (begp:endp)) ;    this%plant_palloc_patch          (:) = nan

    allocate(this%sminp_to_plant_patch        (begp:endp                   )) ; this%sminp_to_plant_patch        (:)   = nan
    allocate(this%plant_pdemand_vr_patch      (begp:endp,1:nlevdecomp_full )) ; this%plant_pdemand_vr_patch      (:,:) = nan
    allocate(this%prev_leafp_to_litter_patch  (begp:endp                   )) ; this%prev_leafp_to_litter_patch  (:)   = nan
    allocate(this%prev_frootp_to_litter_patch (begp:endp                   )) ; this%prev_frootp_to_litter_patch (:)   = nan
    allocate(this%adsorb_to_labilep_vr        (begc:endc,1:nlevdecomp_full )) ; this%adsorb_to_labilep_vr        (:,:) = nan
    allocate(this%desorb_to_solutionp_vr      (begc:endc,1:nlevdecomp_full )) ; this%desorb_to_solutionp_vr      (:,:) = nan
    allocate(this%adsorb_to_labilep_col       (begc:endc                   )) ; this%adsorb_to_labilep_col       (:)   = nan
    allocate(this%desorb_to_solutionp_col     (begc:endc                   )) ; this%desorb_to_solutionp_col     (:)   = nan
    allocate(this%pmpf_decomp_cascade(begc:endc,1:nlevdecomp,1:ndecomp_cascade_transitions)); this%pmpf_decomp_cascade(:,:,:) = nan

    allocate(this%plant_p_uptake_flux         (begc:endc                   )) ; this%plant_p_uptake_flux         (:)   = nan
    allocate(this%soil_p_immob_flux           (begc:endc                   )) ; this%soil_p_immob_flux           (:)   = nan
    allocate(this%soil_p_immob_flux_vr        (begc:endc,1:nlevdecomp_full )) ; this%soil_p_immob_flux_vr        (:,:) = nan
    allocate(this%soil_p_grossmin_flux        (begc:endc                   )) ; this%soil_p_grossmin_flux        (:)   = nan
    allocate(this%smin_p_to_plant_col         (begc:endc                   )) ; this%smin_p_to_plant_col         (:)   = nan
    allocate(this%plant_to_litter_pflux       (begc:endc                   )) ; this%plant_to_litter_pflux       (:)   = nan
    allocate(this%plant_to_cwd_pflux          (begc:endc                   )) ; this%plant_to_cwd_pflux          (:)   = nan

    ! clm_bgc_interface & pflotran
    !------------------------------------------------------------------------
    allocate(this%plant_pdemand_col                 (begc:endc))                                    ; this%plant_pdemand_col                 (:)     = nan
    allocate(this%plant_pdemand_vr_col              (begc:endc,1:nlevdecomp_full))                  ; this%plant_pdemand_vr_col (:,:) = nan
    allocate(this%externalp_to_decomp_ppools_col    (begc:endc, 1:nlevdecomp_full, 1:ndecomp_pools)); this%externalp_to_decomp_ppools_col    (:,:,:) = spval
    allocate(this%externalp_to_decomp_delta_col     (begc:endc))                                    ; this%externalp_to_decomp_delta_col     (:)     = spval
    allocate(this%sminp_net_transport_vr_col        (begc:endc, 1:nlevdecomp_full))                 ; this%sminp_net_transport_vr_col        (:,:)   = spval
    allocate(this%sminp_net_transport_delta_col     (begc:endc))                                    ; this%sminp_net_transport_delta_col     (:)     = spval
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
    use clm_varpar     , only : nlevsno, nlevgrnd, crop_prog
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
    character(10)  :: active
    character(24)  :: fieldname
    character(100) :: longname
    character(8)   :: vr_suffix
    character(1)   :: aa
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

    this%m_leafp_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LEAFP_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='leaf P mortality', &
         ptr_patch=this%m_leafp_to_litter_patch, default='inactive')

    this%m_frootp_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_FROOTP_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='fine root P mortality', &
         ptr_patch=this%m_frootp_to_litter_patch, default='inactive')

    this%m_leafp_storage_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LEAFP_STORAGE_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='leaf P storage mortality', &
         ptr_patch=this%m_leafp_storage_to_litter_patch, default='inactive')

    this%m_frootp_storage_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_FROOTP_STORAGE_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='fine root P storage mortality', &
         ptr_patch=this%m_frootp_storage_to_litter_patch, default='inactive')

    this%m_livestemp_storage_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVESTEMP_STORAGE_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='live stem P storage mortality', &
         ptr_patch=this%m_livestemp_storage_to_litter_patch, default='inactive')

    this%m_deadstemp_storage_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMP_STORAGE_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='dead stem P storage mortality', &
         ptr_patch=this%m_deadstemp_storage_to_litter_patch, default='inactive')

    this%m_livecrootp_storage_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVECROOTP_STORAGE_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='live coarse root P storage mortality', &
         ptr_patch=this%m_livecrootp_storage_to_litter_patch, default='inactive')

    this%m_deadcrootp_storage_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTP_STORAGE_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='dead coarse root P storage mortality', &
         ptr_patch=this%m_deadcrootp_storage_to_litter_patch, default='inactive')

    this%m_leafp_xfer_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LEAFP_XFER_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='leaf P transfer mortality', &
         ptr_patch=this%m_leafp_xfer_to_litter_patch, default='inactive')

    this%m_frootp_xfer_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_FROOTP_XFER_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='fine root P transfer mortality', &
         ptr_patch=this%m_frootp_xfer_to_litter_patch, default='inactive')

    this%m_livestemp_xfer_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVESTEMP_XFER_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='live stem P transfer mortality', &
         ptr_patch=this%m_livestemp_xfer_to_litter_patch, default='inactive')

    this%m_deadstemp_xfer_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMP_XFER_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='dead stem P transfer mortality', &
         ptr_patch=this%m_deadstemp_xfer_to_litter_patch, default='inactive')

    this%m_livecrootp_xfer_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVECROOTP_XFER_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='live coarse root P transfer mortality', &
         ptr_patch=this%m_livecrootp_xfer_to_litter_patch, default='inactive')

    this%m_deadcrootp_xfer_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTP_XFER_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='dead coarse root P transfer mortality', &
         ptr_patch=this%m_deadcrootp_xfer_to_litter_patch, default='inactive')

    this%m_livestemp_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVESTEMP_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='live stem P mortality', &
         ptr_patch=this%m_livestemp_to_litter_patch, default='inactive')

    this%m_deadstemp_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMP_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='dead stem P mortality', &
         ptr_patch=this%m_deadstemp_to_litter_patch, default='inactive')

    this%m_livecrootp_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVECROOTP_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='live coarse root P mortality', &
         ptr_patch=this%m_livecrootp_to_litter_patch, default='inactive')

    this%m_deadcrootp_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTP_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='dead coarse root P mortality', &
         ptr_patch=this%m_deadcrootp_to_litter_patch, default='inactive')

    this%m_retransp_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_RETRANSP_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='retranslocated P pool mortality', &
         ptr_patch=this%m_retransp_to_litter_patch, default='inactive')

    this%m_leafp_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LEAFP_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='leaf P fire loss', &
         ptr_patch=this%m_leafp_to_fire_patch, default='inactive')

    this%m_frootp_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_FROOTP_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='fine root P fire loss ', &
         ptr_patch=this%m_frootp_to_fire_patch, default='inactive')

    this%m_leafp_storage_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LEAFP_STORAGE_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='leaf P storage fire loss', &
         ptr_patch=this%m_leafp_storage_to_fire_patch, default='inactive')

    this%m_frootp_storage_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_FROOTP_STORAGE_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='fine root P storage fire loss', &
         ptr_patch=this%m_frootp_storage_to_fire_patch, default='inactive')

    this%m_livestemp_storage_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVESTEMP_STORAGE_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='live stem P storage fire loss', &
         ptr_patch=this%m_livestemp_storage_to_fire_patch, default='inactive')

    this%m_deadstemp_storage_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMP_STORAGE_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='dead stem P storage fire loss', &
         ptr_patch=this%m_deadstemp_storage_to_fire_patch, default='inactive')

    this%m_livecrootp_storage_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVECROOTP_STORAGE_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='live coarse root P storage fire loss', &
         ptr_patch=this%m_livecrootp_storage_to_fire_patch, default='inactive')

    this%m_deadcrootp_storage_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTP_STORAGE_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='dead coarse root P storage fire loss', &
         ptr_patch=this%m_deadcrootp_storage_to_fire_patch, default='inactive')

    this%m_leafp_xfer_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LEAFP_XFER_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='leaf P transfer fire loss', &
         ptr_patch=this%m_leafp_xfer_to_fire_patch, default='inactive')

    this%m_frootp_xfer_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_FROOTP_XFER_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='fine root P transfer fire loss', &
         ptr_patch=this%m_frootp_xfer_to_fire_patch, default='inactive')

    this%m_livestemp_xfer_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVESTEMP_XFER_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='live stem P transfer fire loss', &
         ptr_patch=this%m_livestemp_xfer_to_fire_patch, default='inactive')

    this%m_deadstemp_xfer_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMP_XFER_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='dead stem P transfer fire loss', &
         ptr_patch=this%m_deadstemp_xfer_to_fire_patch, default='inactive')

    this%m_livecrootp_xfer_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVECROOTP_XFER_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='live coarse root P transfer fire loss', &
         ptr_patch=this%m_livecrootp_xfer_to_fire_patch, default='inactive')

    this%m_deadcrootp_xfer_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTP_XFER_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='dead coarse root P transfer fire loss', &
         ptr_patch=this%m_deadcrootp_xfer_to_fire_patch, default='inactive')

    this%m_livestemp_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVESTEMP_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='live stem P fire loss', &
         ptr_patch=this%m_livestemp_to_fire_patch, default='inactive')

    this%m_deadstemp_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMP_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='dead stem P fire loss', &
         ptr_patch=this%m_deadstemp_to_fire_patch, default='inactive')

    this%m_deadstemp_to_litter_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADSTEMP_TO_LITTER_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='dead stem P fire mortality to litter', &
         ptr_patch=this%m_deadstemp_to_litter_fire_patch, default='inactive')

    this%m_livecrootp_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_LIVECROOTP_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='live coarse root P fire loss', &
         ptr_patch=this%m_livecrootp_to_fire_patch, default='inactive')

    this%m_deadcrootp_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTP_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='dead coarse root P fire loss', &
         ptr_patch=this%m_deadcrootp_to_fire_patch, default='inactive')

    this%m_deadcrootp_to_litter_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_DEADCROOTP_TO_LITTER_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='dead coarse root P fire mortality to litter', &
         ptr_patch=this%m_deadcrootp_to_litter_fire_patch, default='inactive')

    this%m_retransp_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname='M_RETRANSP_TO_FIRE', units='gP/m^2/s', &
         avgflag='A', long_name='retranslocated P pool fire loss', &
         ptr_patch=this%m_retransp_to_fire_patch, default='inactive')

    this%leafp_xfer_to_leafp_patch(begp:endp) = spval
    call hist_addfld1d (fname='LEAFP_XFER_TO_LEAFP', units='gP/m^2/s', &
         avgflag='A', long_name='leaf P growth from storage', &
         ptr_patch=this%leafp_xfer_to_leafp_patch, default='inactive')

    this%frootp_xfer_to_frootp_patch(begp:endp) = spval
    call hist_addfld1d (fname='FROOTP_XFER_TO_FROOTP', units='gP/m^2/s', &
         avgflag='A', long_name='fine root P growth from storage', &
         ptr_patch=this%frootp_xfer_to_frootp_patch, default='inactive')

    this%livestemp_xfer_to_livestemp_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMP_XFER_TO_LIVESTEMP', units='gP/m^2/s', &
         avgflag='A', long_name='live stem P growth from storage', &
         ptr_patch=this%livestemp_xfer_to_livestemp_patch, default='inactive')

    this%deadstemp_xfer_to_deadstemp_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADSTEMP_XFER_TO_DEADSTEMP', units='gP/m^2/s', &
         avgflag='A', long_name='dead stem P growth from storage', &
         ptr_patch=this%deadstemp_xfer_to_deadstemp_patch, default='inactive')

    this%livecrootp_xfer_to_livecrootp_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTP_XFER_TO_LIVECROOTP', units='gP/m^2/s', &
         avgflag='A', long_name='live coarse root P growth from storage', &
         ptr_patch=this%livecrootp_xfer_to_livecrootp_patch, default='inactive')

    this%deadcrootp_xfer_to_deadcrootp_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADCROOTP_XFER_TO_DEADCROOTP', units='gP/m^2/s', &
         avgflag='A', long_name='dead coarse root P growth from storage', &
         ptr_patch=this%deadcrootp_xfer_to_deadcrootp_patch, default='inactive')

    this%leafp_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='LEAFP_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='leaf P litterfall', &
         ptr_patch=this%leafp_to_litter_patch, default='inactive')

    this%leafp_to_retransp_patch(begp:endp) = spval
    call hist_addfld1d (fname='LEAFP_TO_RETRANSP', units='gP/m^2/s', &
         avgflag='A', long_name='leaf P to retranslocated P pool', &
         ptr_patch=this%leafp_to_retransp_patch, default='inactive')

    this%frootp_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname='FROOTP_TO_LITTER', units='gP/m^2/s', &
         avgflag='A', long_name='fine root P litterfall', &
         ptr_patch=this%frootp_to_litter_patch, default='inactive')

    this%retransp_to_ppool_patch(begp:endp) = spval
    call hist_addfld1d (fname='RETRANSP_TO_PPOOL', units='gP/m^2/s', &
         avgflag='A', long_name='deployment of retranslocated P', &
         ptr_patch=this%retransp_to_ppool_patch)

    this%sminp_to_ppool_patch(begp:endp) = spval
    call hist_addfld1d (fname='SMINP_TO_PPOOL', units='gP/m^2/s', &
         avgflag='A', long_name='deployment of soil mineral P uptake', &
         ptr_patch=this%sminp_to_ppool_patch)

    this%ppool_to_leafp_patch(begp:endp) = spval
    call hist_addfld1d (fname='PPOOL_TO_LEAFP', units='gP/m^2/s', &
         avgflag='A', long_name='allocation to leaf P', &
         ptr_patch=this%ppool_to_leafp_patch, default='inactive')

    this%ppool_to_leafp_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='PPOOL_TO_LEAFP_STORAGE', units='gP/m^2/s', &
         avgflag='A', long_name='allocation to leaf P storage', &
         ptr_patch=this%ppool_to_leafp_storage_patch, default='inactive')

    this%ppool_to_frootp_patch(begp:endp) = spval
    call hist_addfld1d (fname='PPOOL_TO_FROOTP', units='gP/m^2/s', &
         avgflag='A', long_name='allocation to fine root P', &
         ptr_patch=this%ppool_to_frootp_patch, default='inactive')

    this%ppool_to_frootp_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='PPOOL_TO_FROOTP_STORAGE', units='gP/m^2/s', &
         avgflag='A', long_name='allocation to fine root P storage', &
         ptr_patch=this%ppool_to_frootp_storage_patch, default='inactive')

    this%ppool_to_livestemp_patch(begp:endp) = spval
    call hist_addfld1d (fname='PPOOL_TO_LIVESTEMP', units='gP/m^2/s', &
         avgflag='A', long_name='allocation to live stem P', &
         ptr_patch=this%ppool_to_livestemp_patch, default='inactive')

    this%ppool_to_livestemp_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='PPOOL_TO_LIVESTEMP_STORAGE', units='gP/m^2/s', &
         avgflag='A', long_name='allocation to live stem P storage', &
         ptr_patch=this%ppool_to_livestemp_storage_patch, default='inactive')

    this%ppool_to_deadstemp_patch(begp:endp) = spval
    call hist_addfld1d (fname='PPOOL_TO_DEADSTEMP', units='gP/m^2/s', &
         avgflag='A', long_name='allocation to dead stem P', &
         ptr_patch=this%ppool_to_deadstemp_patch, default='inactive')

    this%ppool_to_deadstemp_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='PPOOL_TO_DEADSTEMP_STORAGE', units='gP/m^2/s', &
         avgflag='A', long_name='allocation to dead stem P storage', &
         ptr_patch=this%ppool_to_deadstemp_storage_patch, default='inactive')

    this%ppool_to_livecrootp_patch(begp:endp) = spval
    call hist_addfld1d (fname='PPOOL_TO_LIVECROOTP', units='gP/m^2/s', &
         avgflag='A', long_name='allocation to live coarse root P', &
         ptr_patch=this%ppool_to_livecrootp_patch, default='inactive')

    this%ppool_to_livecrootp_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='PPOOL_TO_LIVECROOTP_STORAGE', units='gP/m^2/s', &
         avgflag='A', long_name='allocation to live coarse root P storage', &
         ptr_patch=this%ppool_to_livecrootp_storage_patch, default='inactive')

    this%ppool_to_deadcrootp_patch(begp:endp) = spval
    call hist_addfld1d (fname='PPOOL_TO_DEADCROOTP', units='gP/m^2/s', &
         avgflag='A', long_name='allocation to dead coarse root P', &
         ptr_patch=this%ppool_to_deadcrootp_patch, default='inactive')

    this%ppool_to_deadcrootp_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname='PPOOL_TO_DEADCROOTP_STORAGE', units='gP/m^2/s', &
         avgflag='A', long_name='allocation to dead coarse root P storage', &
         ptr_patch=this%ppool_to_deadcrootp_storage_patch, default='inactive')

    this%leafp_storage_to_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='LEAFP_STORAGE_TO_XFER', units='gP/m^2/s', &
         avgflag='A', long_name='leaf P shift storage to transfer', &
         ptr_patch=this%leafp_storage_to_xfer_patch, default='inactive')

    this%frootp_storage_to_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='FROOTP_STORAGE_TO_XFER', units='gP/m^2/s', &
         avgflag='A', long_name='fine root P shift storage to transfer', &
         ptr_patch=this%frootp_storage_to_xfer_patch, default='inactive')

    this%livestemp_storage_to_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMP_STORAGE_TO_XFER', units='gP/m^2/s', &
         avgflag='A', long_name='live stem P shift storage to transfer', &
         ptr_patch=this%livestemp_storage_to_xfer_patch, default='inactive')

    this%deadstemp_storage_to_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADSTEMP_STORAGE_TO_XFER', units='gP/m^2/s', &
         avgflag='A', long_name='dead stem P shift storage to transfer', &
         ptr_patch=this%deadstemp_storage_to_xfer_patch, default='inactive')

    this%livecrootp_storage_to_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTP_STORAGE_TO_XFER', units='gP/m^2/s', &
         avgflag='A', long_name='live coarse root P shift storage to transfer', &
         ptr_patch=this%livecrootp_storage_to_xfer_patch, default='inactive')

    this%deadcrootp_storage_to_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname='DEADCROOTP_STORAGE_TO_XFER', units='gP/m^2/s', &
         avgflag='A', long_name='dead coarse root P shift storage to transfer', &
         ptr_patch=this%deadcrootp_storage_to_xfer_patch, default='inactive')

    this%livestemp_to_deadstemp_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMP_TO_DEADSTEMP', units='gP/m^2/s', &
         avgflag='A', long_name='live stem P turnover', &
         ptr_patch=this%livestemp_to_deadstemp_patch, default='inactive')

    this%livestemp_to_retransp_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVESTEMP_TO_RETRANSP', units='gP/m^2/s', &
         avgflag='A', long_name='live stem P to retranslocated P pool', &
         ptr_patch=this%livestemp_to_retransp_patch, default='inactive')

    this%livecrootp_to_deadcrootp_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTP_TO_DEADCROOTP', units='gP/m^2/s', &
         avgflag='A', long_name='live coarse root P turnover', &
         ptr_patch=this%livecrootp_to_deadcrootp_patch, default='inactive')

    this%livecrootp_to_retransp_patch(begp:endp) = spval
    call hist_addfld1d (fname='LIVECROOTP_TO_RETRANSP', units='gP/m^2/s', &
         avgflag='A', long_name='live coarse root P to retranslocated N pool', &
         ptr_patch=this%livecrootp_to_retransp_patch, default='inactive')

    this%pdeploy_patch(begp:endp) = spval
    call hist_addfld1d (fname='PDEPLOY', units='gP/m^2/s', &
         avgflag='A', long_name='total P deployed in new growth', &
         ptr_patch=this%pdeploy_patch)

    this%wood_harvestp_patch(begp:endp) = spval
    call hist_addfld1d (fname='WOOD_HARVESTP', units='gP/m^2/s', &
         avgflag='A', long_name='wood harvest P (to product pools)', &
         ptr_patch=this%wood_harvestp_patch, default='inactive')

    this%fire_ploss_patch(begp:endp) = spval
    call hist_addfld1d (fname='PFT_FIRE_PLOSS', units='gP/m^2/s', &
         avgflag='A', long_name='total pft-level fire P loss', &
         ptr_patch=this%fire_ploss_patch, default='inactive')

    if (crop_prog) then
       this%fert_p_patch(begp:endp) = spval
       call hist_addfld1d (fname='FERT_P', units='gP/m^2/s', &
            avgflag='A', long_name='fertilizer P added', &
            ptr_patch=this%fert_p_patch)
    end if


    if (crop_prog) then
       this%fert_p_counter_patch(begp:endp) = spval
       call hist_addfld1d (fname='FERT_P_COUNTER', units='seconds', &
            avgflag='A', long_name='time left to fertilize', &
            ptr_patch=this%fert_p_counter_patch)
    end if

    !-------------------------------
    ! P flux variables - native to column
    !-------------------------------

    this%pdep_to_sminp_col(begc:endc) = spval
    call hist_addfld1d (fname='PDEP_TO_SMINP', units='gP/m^2/s', &
         avgflag='A', long_name='atmospheric P deposition to soil mineral P', &
         ptr_col=this%pdep_to_sminp_col)


    do k = 1, ndecomp_pools
       if ( decomp_cascade_con%is_litter(k) .or. decomp_cascade_con%is_cwd(k) ) then
          this%m_decomp_ppools_to_fire_col(begc:endc,k) = spval
          data1dptr => this%m_decomp_ppools_to_fire_col(:,k)
          fieldname = 'M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'P_TO_FIRE'
          longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' P fire loss'
          call hist_addfld1d (fname=fieldname, units='gP/m^2',  &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr, default='inactive')

          if ( nlevdecomp_full > 1 ) then
             this%m_decomp_ppools_to_fire_vr_col(begc:endc,:,k) = spval
             data2dptr => this%m_decomp_ppools_to_fire_vr_col(:,:,k)
             fieldname = 'M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'P_TO_FIRE'//trim(vr_suffix)
             longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' P fire loss'
             call hist_addfld_decomp (fname=fieldname, units='gP/m^3',  type2d='levdcmp', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data2dptr, default='inactive')
          endif
       endif

       !! biochemical P mineralization for each soil pool at each soil level for
       !! each column
       if ( k >= 5 )then

          if ( nlevdecomp_full > 1 ) then
             this%biochem_pmin_ppools_vr_col(begc:endc,:,k) = spval
             data2dptr => this%biochem_pmin_ppools_vr_col(:,:,k)
             write(aa,'(i1)') k
             fieldname = 'BIOCHEM_PMIN_PPOOL'//aa//trim(vr_suffix)
             longname  = 'Biochemical mineralization of ppool'//aa
             call hist_addfld_decomp (fname=fieldname, units='gP/m^2/s',type2d='levdcmp', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data2dptr, default='inactive')
          endif
       endif
    end do
    do l = 1, ndecomp_cascade_transitions
       ! vertically integrated fluxes
       !-- mineralization/immobilization fluxes (none from CWD)
       if ( .not. decomp_cascade_con%is_cwd(decomp_cascade_con%cascade_donor_pool(l)) ) then
          this%decomp_cascade_sminp_flux_col(begc:endc,l) = spval
          data1dptr => this%decomp_cascade_sminp_flux_col(:,l)
          if ( decomp_cascade_con%cascade_receiver_pool(l) /= 0 ) then
             fieldname = 'SMINP_TO_'//&
                  trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))//'P_'//&
                  trim(decomp_cascade_con%decomp_pool_name_short(decomp_cascade_con%cascade_donor_pool(l)))
             longname =  'mineral P flux for decomp. of '&
                  //trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//&
                  'to '//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))
          else
             fieldname = trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))&
                  //'P_TO_SMINP'
             longname =  'mineral P flux for decomp. of '&
                  //trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))
          endif
          call hist_addfld1d (fname=fieldname, units='gP/m^2', &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr)
       end if

       !-- transfer fluxes (none from terminal pool, if present)
       if ( decomp_cascade_con%cascade_receiver_pool(l) /= 0 ) then
          this%decomp_cascade_ptransfer_col(begc:endc,l) = spval
          data1dptr => this%decomp_cascade_ptransfer_col(:,l)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'P_TO_'//&
               trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))//'P'
          longname =  'decomp. of '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))//&
               ' P to '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_receiver_pool(l)))//' N'
          call hist_addfld1d (fname=fieldname, units='gP/m^2',  &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr)
       end if

       ! vertically resolved fluxes
       if ( nlevdecomp_full > 1 ) then
          !-- mineralization/immobilization fluxes (none from CWD)
          if ( .not. decomp_cascade_con%is_cwd(decomp_cascade_con%cascade_donor_pool(l)) ) then
             this%decomp_cascade_sminp_flux_vr_col(begc:endc,:,l) = spval
             data2dptr => this%decomp_cascade_sminp_flux_vr_col(:,:,l)
             if ( decomp_cascade_con%cascade_receiver_pool(l) /= 0 ) then
                fieldname = 'SMINP_TO_'&
                     //trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))//'P_'//&
                     trim(decomp_cascade_con%decomp_pool_name_short(decomp_cascade_con%cascade_donor_pool(l)))//trim(vr_suffix)
                longname =  'mineral P flux for decomp. of '&
                     //trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//&
                     'to '//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))
             else
                fieldname = trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))&
                     //'P_TO_SMINP'//trim(vr_suffix)
                longname =  'mineral P flux for decomp. of '&
                     //trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))
             endif
             call hist_addfld_decomp (fname=fieldname, units='gP/m^3',  type2d='levdcmp', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data2dptr, default='inactive')
          endif

          !-- transfer fluxes (none from terminal pool, if present)
          if ( decomp_cascade_con%cascade_receiver_pool(l) /= 0 ) then
             this%decomp_cascade_ptransfer_vr_col(begc:endc,:,l) = spval
             data2dptr => this%decomp_cascade_ptransfer_vr_col(:,:,l)
             fieldname = trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'P_TO_'//&
                  trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))&
                  //'P'//trim(vr_suffix)
             longname =  'decomp. of '&
                  //trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))//&
                  ' P to '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_receiver_pool(l)))//' P'
             call hist_addfld_decomp (fname=fieldname, units='gP/m^3',  type2d='levdcmp', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data2dptr, default='inactive')
          endif

       endif
    end do

    this%som_p_leached_col(begc:endc) = spval
    call hist_addfld1d (fname='SOM_P_LEACHED', units='gP/m^2/s', &
         avgflag='A', long_name='total flux of P from SOM pools due to leaching', &
         ptr_col=this%som_p_leached_col, default='inactive')

    do k = 1, ndecomp_pools
       if ( .not. decomp_cascade_con%is_cwd(k) ) then
          this%decomp_ppools_leached_col(begc:endc,k) = spval
          data1dptr => this%decomp_ppools_leached_col(:,k)
          fieldname = 'M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'P_TO_LEACHING'
          longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' P leaching loss'
          call hist_addfld1d (fname=fieldname, units='gP/m^2/s', &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr, default='inactive')

          this%decomp_ppools_transport_tendency_col(begc:endc,:,k) = spval
          data2dptr => this%decomp_ppools_transport_tendency_col(:,:,k)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(k))//'P_TNDNCY_VERT_TRANSPORT'
          longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' P tendency due to vertical transport'
          call hist_addfld_decomp (fname=fieldname, units='gP/m^3/s',  type2d='levdcmp', &
               avgflag='A', long_name=longname, &
               ptr_col=data2dptr)
       end if
    end do

    this%primp_to_labilep_col(begc:endc) = spval
    call hist_addfld1d (fname='PRIMP_TO_LABILEP', units='gP/m^2/s',   &
         avgflag='A', long_name='PRIMARY MINERAL P TO LABILE P', &
         ptr_col=this%primp_to_labilep_col)

    if ( nlevdecomp_full > 1 ) then
       this%primp_to_labilep_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='PRIMP_TO_LABILEP'//trim(vr_suffix), units='gP/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='PRIMARY MINERAL P TO LABILE P', &
            ptr_col=this%primp_to_labilep_vr_col, default='inactive')
    endif

    this%labilep_to_secondp_col(begc:endc) = spval
    call hist_addfld1d (fname='LABILEP_TO_SECONDP', units='gP/m^2/s',   &
         avgflag='A', long_name='LABILE P TO SECONDARY MINERAL P', &
         ptr_col=this%labilep_to_secondp_col)

    if ( nlevdecomp_full > 1 ) then
       this%labilep_to_secondp_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='LABILEP_TO_SECONDP'//trim(vr_suffix), units='gP/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='LABILE P TO SECONDARY MINERAL P', &
            ptr_col=this%labilep_to_secondp_vr_col, default='inactive')
    endif


    this%secondp_to_labilep_col(begc:endc) = spval
    call hist_addfld1d (fname='SECONDP_TO_LABILEP', units='gP/m^2/s',   &
         avgflag='A', long_name='SECONDARY MINERAL P TO LABILE P', &
         ptr_col=this%secondp_to_labilep_col)

    if ( nlevdecomp_full > 1 ) then
       this%secondp_to_labilep_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SECONDP_TO_LABILEP'//trim(vr_suffix), units='gP/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='SECONDARY MINERAL P TO LABILE P', &
            ptr_col=this%secondp_to_labilep_vr_col, default='inactive')
    endif

    this%secondp_to_occlp_col(begc:endc) = spval
    call hist_addfld1d (fname='SECONDP_TO_OCCLP', units='gP/m^2/s',   &
         avgflag='A', long_name='SECONDARY MINERAL P TO OCCLUDED P', &
         ptr_col=this%secondp_to_occlp_col)

    if ( nlevdecomp_full > 1 ) then
       this%secondp_to_occlp_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SECONDP_TO_OCCLP'//trim(vr_suffix), units='gP/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='SECONDARY MINERAL P TO OCCLUDED P', &
            ptr_col=this%secondp_to_occlp_vr_col, default='inactive')
    endif

    this%sminp_leached_col(begc:endc) = spval
    call hist_addfld1d (fname='SMINP_LEACHED', units='gP/m^2/s',   &
         avgflag='A', long_name='soil mineral P pool loss to leaching', &
         ptr_col=this%sminp_leached_col)

    if ( nlevdecomp_full > 1 ) then
       this%sminp_leached_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SMINP_LEACHED'//trim(vr_suffix), units='gP/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='soil mineral P pool loss to leaching', &
            ptr_col=this%sminp_leached_vr_col, default='inactive')
    endif


    if ( nlevdecomp_full > 1 ) then
       this%potential_immob_p_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='POTENTIAL_IMMOB_P'//trim(vr_suffix), units='gP/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='potential P immobilization', &
            ptr_col=this%potential_immob_p_vr_col, default='inactive')
    end if

    if ( nlevdecomp_full > 1 ) then
       this%actual_immob_p_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='ACTUAL_IMMOB_P'//trim(vr_suffix), units='gP/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='actual P immobilization', &
            ptr_col=this%actual_immob_p_vr_col, default='inactive')
    end if

    if ( nlevdecomp_full > 1 ) then
       this%sminp_to_plant_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SMINP_TO_PLANT'//trim(vr_suffix), units='gP/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='plant uptake of soil mineral P', &
            ptr_col=this%sminp_to_plant_vr_col, default='inactive')
    end if

    if ( nlevdecomp_full > 1 ) then
       this%supplement_to_sminp_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SUPPLEMENT_TO_SMINP'//trim(vr_suffix), units='gP/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='supplemental P supply', &
            ptr_col=this%supplement_to_sminp_vr_col, default='inactive')
    end if

    if ( nlevdecomp_full > 1 ) then
       this%gross_pmin_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='GROSS_PMIN'//trim(vr_suffix), units='gP/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='gross rate of P mineralization', &
            ptr_col=this%gross_pmin_vr_col, default='inactive')
    end if

    if ( nlevdecomp_full > 1 ) then
       this%net_pmin_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='NET_PMIN'//trim(vr_suffix), units='gP/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='net rate of P mineralization', &
            ptr_col=this%net_pmin_vr_col, default='inactive')
    end if

    if ( nlevdecomp_full > 1 ) then
       this%biochem_pmin_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='BIOCHEM_PMIN'//trim(vr_suffix), units='gP/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='biochemical rate of P mineralization', &
            ptr_col=this%biochem_pmin_vr_col, default='inactive')
    end if

    this%potential_immob_p_col(begc:endc) = spval
    call hist_addfld1d (fname='POTENTIAL_IMMOB_P', units='gP/m^2/s', &
         avgflag='A', long_name='potential P immobilization', &
         ptr_col=this%potential_immob_p_col)

    this%actual_immob_p_col(begc:endc) = spval
    call hist_addfld1d (fname='ACTUAL_IMMOB_P', units='gP/m^2/s', &
         avgflag='A', long_name='actual P immobilization', &
         ptr_col=this%actual_immob_p_col)

    this%sminp_to_plant_col(begc:endc) = spval
    call hist_addfld1d (fname='SMINP_TO_PLANT', units='gP/m^2/s', &
         avgflag='A', long_name='plant uptake of soil mineral P', &
         ptr_col=this%sminp_to_plant_col)

    this%supplement_to_sminp_col(begc:endc) = spval
    call hist_addfld1d (fname='SUPPLEMENT_TO_SMINP', units='gP/m^2/s', &
         avgflag='A', long_name='supplemental P supply', &
         ptr_col=this%supplement_to_sminp_col)

    this%gross_pmin_col(begc:endc) = spval
    call hist_addfld1d (fname='GROSS_PMIN', units='gP/m^2/s', &
         avgflag='A', long_name='gross rate of P mineralization', &
         ptr_col=this%gross_pmin_col)

    this%net_pmin_col(begc:endc) = spval
    call hist_addfld1d (fname='NET_PMIN', units='gP/m^2/s', &
         avgflag='A', long_name='net rate of P mineralization', &
         ptr_col=this%net_pmin_col)

    this%biochem_pmin_col(begc:endc) = spval
    call hist_addfld1d (fname='BIOCHEM_PMIN', units='gP/m^2/s', &
         avgflag='A', long_name='biochemical rate of P mineralization', &
         ptr_col=this%biochem_pmin_col)

    this%fire_ploss_col(begc:endc) = spval
    call hist_addfld1d (fname='COL_FIRE_PLOSS', units='gP/m^2/s', &
         avgflag='A', long_name='total column-level fire P loss', &
         ptr_col=this%fire_ploss_col, default='inactive')

    this%fire_decomp_ploss_col(begc:endc) = spval
    call hist_addfld1d (fname='DECOMP_FIRE_PLOSS', units='gP/m^2/s', &
         avgflag='A', long_name='fire P loss of decomposable pools', &
         ptr_col=this%fire_decomp_ploss_col, default='inactive')

    this%dwt_seedp_to_leaf_col(begc:endc) = spval
    call hist_addfld1d (fname='DWT_SEEDP_TO_LEAF', units='gP/m^2/s', &
         avgflag='A', long_name='seed source to PFT-level leaf', &
         ptr_col=this%dwt_seedp_to_leaf_col, default='inactive')

    this%dwt_seedp_to_deadstem_col(begc:endc) = spval
    call hist_addfld1d (fname='DWT_SEEDP_TO_DEADSTEM', units='gP/m^2/s', &
         avgflag='A', long_name='seed source to PFT-level deadstem', &
         ptr_col=this%dwt_seedp_to_deadstem_col, default='inactive')

    this%dwt_conv_pflux_col(begc:endc) = spval
    call hist_addfld1d (fname='DWT_CONV_PFLUX', units='gP/m^2/s', &
         avgflag='A', long_name='conversion P flux (immediate loss to atm)', &
         ptr_col=this%dwt_conv_pflux_col, default='inactive')

    this%dwt_prod10p_gain_col(begc:endc) = spval
    call hist_addfld1d (fname='DWT_PROD10P_GAIN', units='gP/m^2/s', &
         avgflag='A', long_name='addition to 10-yr wood product pool', &
         ptr_col=this%dwt_prod10p_gain_col, default='inactive')

    this%prod10p_loss_col(begc:endc) = spval
    call hist_addfld1d (fname='PROD10P_LOSS', units='gP/m^2/s', &
         avgflag='A', long_name='loss from 10-yr wood product pool', &
         ptr_col=this%prod10p_loss_col, default='inactive')

    this%dwt_prod100p_gain_col(begc:endc) = spval
    call hist_addfld1d (fname='DWT_PROD100P_GAIN', units='gP/m^2/s', &
         avgflag='A', long_name='addition to 100-yr wood product pool', &
         ptr_col=this%dwt_prod100p_gain_col, default='inactive')

    this%prod100p_loss_col(begc:endc) = spval
    call hist_addfld1d (fname='PROD100P_LOSS', units='gP/m^2/s', &
         avgflag='A', long_name='loss from 100-yr wood product pool', &
         ptr_col=this%prod100p_loss_col, default='inactive')

    this%prod1p_loss_col(begc:endc) = spval
    call hist_addfld1d (fname='PROD1P_LOSS', units='gP/m^2/s', &
         avgflag='A', long_name='loss from 1-yr crop product pool', &
         ptr_col=this%prod1p_loss_col)

    this%product_ploss_col(begc:endc) = spval
    call hist_addfld1d (fname='PRODUCT_PLOSS', units='gP/m^2/s', &
         avgflag='A', long_name='total P loss from wood product pools', &
         ptr_col=this%product_ploss_col, default='inactive')

    this%dwt_frootp_to_litr_met_p_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='DWT_FROOTP_TO_LITR_MET_P', units='gP/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='fine root to litter due to landcover change', &
         ptr_col=this%dwt_frootp_to_litr_met_p_col, default='inactive')

    this%dwt_frootp_to_litr_cel_p_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='DWT_FROOTP_TO_LITR_CEL_P', units='gP/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='fine root to litter due to landcover change', &
         ptr_col=this%dwt_frootp_to_litr_cel_p_col, default='inactive')

    this%dwt_frootp_to_litr_lig_p_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='DWT_FROOTP_TO_LITR_LIG_P', units='gP/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='fine root to litter due to landcover change', &
         ptr_col=this%dwt_frootp_to_litr_lig_p_col, default='inactive')

    this%dwt_livecrootp_to_cwdp_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='DWT_LIVECROOTP_TO_CWDP', units='gP/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='live coarse root to CWD due to landcover change', &
         ptr_col=this%dwt_livecrootp_to_cwdp_col, default='inactive')

    this%dwt_deadcrootp_to_cwdp_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='DWT_DEADCROOTP_TO_CWDP', units='gP/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='dead coarse root to CWD due to landcover change', &
         ptr_col=this%dwt_deadcrootp_to_cwdp_col, default='inactive')

    this%dwt_ploss_col(begc:endc) = spval
    call hist_addfld1d (fname='DWT_PLOSS', units='gP/m^2/s', &
         avgflag='A', long_name='total phosphorus loss from landcover conversion', &
         ptr_col=this%dwt_ploss_col, default='inactive')

    if (crop_prog) then
       this%fert_p_to_sminp_col(begc:endc) = spval
       call hist_addfld1d (fname='FERT_TO_LABILEP', units='gP/m^2/s', &
            avgflag='A', long_name='fertilizer to soil mineral P', &
            ptr_col=this%fert_p_to_sminp_col)
    end if


    this%plant_pdemand_patch(begp:endp) = spval
    call hist_addfld1d (fname='PLANT_PDEMAND', units='gP/m^2/s', &
         avgflag='A', long_name='P flux required to support initial GPP', &
         ptr_patch=this%plant_pdemand_patch)

    this%avail_retransp_patch(begp:endp) = spval
    call hist_addfld1d (fname='AVAIL_RETRANSP', units='gP/m^2/s', &
         avgflag='A', long_name='P flux available from retranslocation pool', &
         ptr_patch=this%avail_retransp_patch, default='active')

    this%plant_palloc_patch(begp:endp) = spval
    call hist_addfld1d (fname='PLANT_PALLOC', units='gP/m^2/s', &
         avgflag='A', long_name='total allocated P flux', &
         ptr_patch=this%plant_palloc_patch, default='active')

    !! bgc interface
    this%plant_pdemand_col(begc:endc) = spval
    call hist_addfld1d (fname='PLANT_PDEMAND_COL', units='gN/m^2/s', &
        avgflag='A', long_name='P flux required to support initial GPP', &
        ptr_col=this%plant_pdemand_col)

    this%adsorb_to_labilep_col(begc:endc) = spval
    call hist_addfld1d (fname='ADSORBTION_P', units='gP/m^2/s', &
         avgflag='A', long_name='adsorb P flux', &
         ptr_col=this%adsorb_to_labilep_col, default='active')

    this%desorb_to_solutionp_col(begc:endc) = spval
    call hist_addfld1d (fname='DESORPTION_P', units='gP/m^2/s', &
         avgflag='A', long_name='desorp P flux', &
         ptr_col=this%desorb_to_solutionp_col, default='active')

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initializes time varying variables used only in coupled carbon-phosphorus mode (CN):
    !
    ! !USES:
    use clm_varpar      , only : crop_prog
    use landunit_varcon , only : istsoil, istcrop
    !
    ! !ARGUMENTS:
    class(phosphorusflux_type) :: this
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
    ! initialize phosphorus flux variables
    !-----------------------------------------------

    do p = bounds%begp,bounds%endp
       l = pft%landunit(p)

       this%prev_leafp_to_litter_patch (p)  = 0._r8
       this%prev_frootp_to_litter_patch(p)  = 0._r8

       if ( crop_prog )then
          this%fert_p_counter_patch(p)  = spval
          this%fert_p_patch(p)          = 0._r8
       end if

       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          this%fert_p_counter_patch(p)  = 0._r8
       end if

       if (lun%ifspecial(l)) then
          this%plant_pdemand_patch(p)  = spval
          this%avail_retransp_patch(p) = spval
          this%plant_palloc_patch(p)   = spval
       end if
    end do

    ! initialize fields for special filters

    do fc = 1,num_special_col
       c = special_col(fc)
       this%dwt_ploss_col(c) = 0._r8
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

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag, varname='fert_p_counter', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%fert_p_counter_patch)

       call restartvar(ncid=ncid, flag=flag, varname='fert_p', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%fert_p_patch)
    end if

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag,  varname='grainp_xfer_to_grainp', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='grain P growth from storage', units='gP/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%grainp_xfer_to_grainp_patch)
    end if

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag,  varname='livestemp_to_litter', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='livestem P to litter', units='gP/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemp_to_litter_patch)
    end if

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag,  varname='grainp_to_food', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='grain P to food', units='gP/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%grainp_to_food_patch)
    end if

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag,  varname='ppool_to_grainp', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='allocation to grain P', units='gP/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%ppool_to_grainp_patch)
    end if

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag,  varname='ppool_to_grainp_storage', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='allocation to grain P storage', units='gP/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%ppool_to_grainp_storage_patch)
    end if

    if (crop_prog) then
       call restartvar(ncid=ncid, flag=flag, varname='grainp_storage_to_xfer', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='grain P shift storage to transfer', units='gP/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%grainp_storage_to_xfer_patch)
    end if

    call restartvar(ncid=ncid, flag=flag, varname='plant_pdemand', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%plant_pdemand_patch)

    call restartvar(ncid=ncid, flag=flag, varname='avail_retransp', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%avail_retransp_patch)

    call restartvar(ncid=ncid, flag=flag, varname='plant_palloc', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%plant_palloc_patch)

    ! clm_bgc_interface & pflotran
    !------------------------------------------------------------------------
    if (use_pflotran .and. pf_cmode) then
       ! externalp_to_decomp_ppools_col
       do k = 1, ndecomp_pools
          varname=trim(decomp_cascade_con%decomp_pool_name_restart(k))//'external_p'
          if (use_vertsoilc) then
             ptr2d => this%externalp_to_decomp_ppools_col(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_vr", xtype=ncd_double,  &
                  dim1name='column', dim2name='levgrnd', switchdim=.true., &
                  long_name='net organic P adding/removal/transport to soil', units='gP/m3/s', fill_value=spval, &
                  interpinic_flag='interp', readvar=readvar, data=ptr2d)
          else
             ptr1d => this%externalp_to_decomp_ppools_col(:,1,k) ! nlevdecomp = 1; so treat as 1D variable
             call restartvar(ncid=ncid, flag=flag, varname=varname, xtype=ncd_double,  &
                  dim1name='column', &
                  long_name='net organic P adding/removal/transport to soil', units='gP/m3/s', fill_value=spval, &
                  interpinic_flag='interp' , readvar=readvar, data=ptr1d)
          end if
          if (flag=='read' .and. .not. readvar) then
          !   call endrun(msg='ERROR:: '//trim(varname)//' is required on an initialization dataset'//&
          !        errMsg(__FILE__, __LINE__))
             this%externalp_to_decomp_ppools_col(:,:,k) = 0._r8
          end if
       end do

       !sminp_net_transport_vr
       if (.not.pf_hmode) then
          if (use_vertsoilc) then
             ptr2d => this%sminp_net_transport_vr_col(:,:)
             call restartvar(ncid=ncid, flag=flag, varname='sminp_net_transport_vr', xtype=ncd_double, &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='net soil mineral-P transport', units='gP/m3/s', &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
          else
             ptr1d => this%sminp_net_transport_vr_col(:,1)
             call restartvar(ncid=ncid, flag=flag, varname='sminp_net_transport_vr', xtype=ncd_double, &
               dim1name='column', &
               long_name='net soil  mineral-P transport', units='gP/m3/s', &
               interpinic_flag='interp', readvar=readvar, data=ptr1d)
          end if
          if (flag=='read' .and. .not. readvar) then
          !   call endrun(msg='ERROR:: no3_net_transport_vr'//' is required on an initialization dataset'//&
          !     errMsg(__FILE__, __LINE__))
             this%sminp_net_transport_vr_col(:,:) = 0._r8
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

    do fi = 1,num_patch
       i=filter_patch(fi)

       this%m_leafp_to_litter_patch(i)                   = value_patch
       this%m_frootp_to_litter_patch(i)                  = value_patch
       this%m_leafp_storage_to_litter_patch(i)           = value_patch
       this%m_frootp_storage_to_litter_patch(i)          = value_patch
       this%m_livestemp_storage_to_litter_patch(i)       = value_patch
       this%m_deadstemp_storage_to_litter_patch(i)       = value_patch
       this%m_livecrootp_storage_to_litter_patch(i)      = value_patch
       this%m_deadcrootp_storage_to_litter_patch(i)      = value_patch
       this%m_leafp_xfer_to_litter_patch(i)              = value_patch
       this%m_frootp_xfer_to_litter_patch(i)             = value_patch
       this%m_livestemp_xfer_to_litter_patch(i)          = value_patch
       this%m_deadstemp_xfer_to_litter_patch(i)          = value_patch
       this%m_livecrootp_xfer_to_litter_patch(i)         = value_patch
       this%m_deadcrootp_xfer_to_litter_patch(i)         = value_patch
       this%m_livestemp_to_litter_patch(i)               = value_patch
       this%m_deadstemp_to_litter_patch(i)               = value_patch
       this%m_livecrootp_to_litter_patch(i)              = value_patch
       this%m_deadcrootp_to_litter_patch(i)              = value_patch
       this%m_retransp_to_litter_patch(i)                = value_patch
       this%hrv_leafp_to_litter_patch(i)                 = value_patch
       this%hrv_frootp_to_litter_patch(i)                = value_patch
       this%hrv_leafp_storage_to_litter_patch(i)         = value_patch
       this%hrv_frootp_storage_to_litter_patch(i)        = value_patch
       this%hrv_livestemp_storage_to_litter_patch(i)     = value_patch
       this%hrv_deadstemp_storage_to_litter_patch(i)     = value_patch
       this%hrv_livecrootp_storage_to_litter_patch(i)    = value_patch
       this%hrv_deadcrootp_storage_to_litter_patch(i)    = value_patch
       this%hrv_leafp_xfer_to_litter_patch(i)            = value_patch
       this%hrv_frootp_xfer_to_litter_patch(i)           = value_patch
       this%hrv_livestemp_xfer_to_litter_patch(i)        = value_patch
       this%hrv_deadstemp_xfer_to_litter_patch(i)        = value_patch
       this%hrv_livecrootp_xfer_to_litter_patch(i)       = value_patch
       this%hrv_deadcrootp_xfer_to_litter_patch(i)       = value_patch
       this%hrv_livestemp_to_litter_patch(i)             = value_patch
       this%hrv_deadstemp_to_prod10p_patch(i)            = value_patch
       this%hrv_deadstemp_to_prod100p_patch(i)           = value_patch
       this%hrv_leafp_to_prod1p_patch(i)                 = value_patch
       this%hrv_livestemp_to_prod1p_patch(i)             = value_patch
       this%hrv_grainp_to_prod1p_patch(i)                = value_patch
       this%hrv_cropp_to_prod1p_patch(i)                 = value_patch
       this%hrv_livecrootp_to_litter_patch(i)            = value_patch
       this%hrv_deadcrootp_to_litter_patch(i)            = value_patch
       this%hrv_retransp_to_litter_patch(i)              = value_patch

       this%m_leafp_to_fire_patch(i)                     = value_patch
       this%m_leafp_storage_to_fire_patch(i)             = value_patch
       this%m_leafp_xfer_to_fire_patch(i)                = value_patch
       this%m_livestemp_to_fire_patch(i)                 = value_patch
       this%m_livestemp_storage_to_fire_patch(i)         = value_patch
       this%m_livestemp_xfer_to_fire_patch(i)            = value_patch
       this%m_deadstemp_to_fire_patch(i)                 = value_patch
       this%m_deadstemp_storage_to_fire_patch(i)         = value_patch
       this%m_deadstemp_xfer_to_fire_patch(i)            = value_patch
       this%m_frootp_to_fire_patch(i)                    = value_patch
       this%m_frootp_storage_to_fire_patch(i)            = value_patch
       this%m_frootp_xfer_to_fire_patch(i)               = value_patch
       this%m_livecrootp_to_fire_patch(i)                = value_patch
       this%m_livecrootp_storage_to_fire_patch(i)        = value_patch
       this%m_livecrootp_xfer_to_fire_patch(i)           = value_patch
       this%m_deadcrootp_to_fire_patch(i)                = value_patch
       this%m_deadcrootp_storage_to_fire_patch(i)        = value_patch
       this%m_deadcrootp_xfer_to_fire_patch(i)           = value_patch
       this%m_retransp_to_fire_patch(i)                  = value_patch


       this%m_leafp_to_litter_fire_patch(i)              = value_patch
       this%m_leafp_storage_to_litter_fire_patch(i)      = value_patch
       this%m_leafp_xfer_to_litter_fire_patch(i)         = value_patch
       this%m_livestemp_to_litter_fire_patch(i)          = value_patch
       this%m_livestemp_storage_to_litter_fire_patch(i)  = value_patch
       this%m_livestemp_xfer_to_litter_fire_patch(i)     = value_patch
       this%m_livestemp_to_deadstemp_fire_patch(i)       = value_patch
       this%m_deadstemp_to_litter_fire_patch(i)          = value_patch
       this%m_deadstemp_storage_to_litter_fire_patch(i)  = value_patch
       this%m_deadstemp_xfer_to_litter_fire_patch(i)     = value_patch
       this%m_frootp_to_litter_fire_patch(i)             = value_patch
       this%m_frootp_storage_to_litter_fire_patch(i)     = value_patch
       this%m_frootp_xfer_to_litter_fire_patch(i)        = value_patch
       this%m_livecrootp_to_litter_fire_patch(i)         = value_patch
       this%m_livecrootp_storage_to_litter_fire_patch(i) = value_patch
       this%m_livecrootp_xfer_to_litter_fire_patch(i)    = value_patch
       this%m_livecrootp_to_deadcrootp_fire_patch(i)     = value_patch
       this%m_deadcrootp_to_litter_fire_patch(i)         = value_patch
       this%m_deadcrootp_storage_to_litter_fire_patch(i) = value_patch
       this%m_deadcrootp_xfer_to_litter_fire_patch(i)    = value_patch
       this%m_retransp_to_litter_fire_patch(i)           = value_patch

       this%leafp_xfer_to_leafp_patch(i)                 = value_patch
       this%frootp_xfer_to_frootp_patch(i)               = value_patch
       this%livestemp_xfer_to_livestemp_patch(i)         = value_patch
       this%deadstemp_xfer_to_deadstemp_patch(i)         = value_patch
       this%livecrootp_xfer_to_livecrootp_patch(i)       = value_patch
       this%deadcrootp_xfer_to_deadcrootp_patch(i)       = value_patch
       this%leafp_to_litter_patch(i)                     = value_patch
       this%leafp_to_retransp_patch(i)                   = value_patch
       this%frootp_to_litter_patch(i)                    = value_patch
       this%retransp_to_ppool_patch(i)                   = value_patch
       this%sminp_to_ppool_patch(i)                      = value_patch
       this%ppool_to_leafp_patch(i)                      = value_patch
       this%ppool_to_leafp_storage_patch(i)              = value_patch
       this%ppool_to_frootp_patch(i)                     = value_patch
       this%ppool_to_frootp_storage_patch(i)             = value_patch
       this%ppool_to_livestemp_patch(i)                  = value_patch
       this%ppool_to_livestemp_storage_patch(i)          = value_patch
       this%ppool_to_deadstemp_patch(i)                  = value_patch
       this%ppool_to_deadstemp_storage_patch(i)          = value_patch
       this%ppool_to_livecrootp_patch(i)                 = value_patch
       this%ppool_to_livecrootp_storage_patch(i)         = value_patch
       this%ppool_to_deadcrootp_patch(i)                 = value_patch
       this%ppool_to_deadcrootp_storage_patch(i)         = value_patch
       this%leafp_storage_to_xfer_patch(i)               = value_patch
       this%frootp_storage_to_xfer_patch(i)              = value_patch
       this%livestemp_storage_to_xfer_patch(i)           = value_patch
       this%deadstemp_storage_to_xfer_patch(i)           = value_patch
       this%livecrootp_storage_to_xfer_patch(i)          = value_patch
       this%deadcrootp_storage_to_xfer_patch(i)          = value_patch
       this%livestemp_to_deadstemp_patch(i)              = value_patch
       this%livestemp_to_retransp_patch(i)               = value_patch
       this%livecrootp_to_deadcrootp_patch(i)            = value_patch
       this%livecrootp_to_retransp_patch(i)              = value_patch
       this%pdeploy_patch(i)                             = value_patch
       this%pinputs_patch(i)                             = value_patch
       this%poutputs_patch(i)                            = value_patch
       this%wood_harvestp_patch(i)                       = value_patch
       this%fire_ploss_patch(i)                          = value_patch
    end do

    if ( crop_prog )then
       do fi = 1,num_patch
          i = filter_patch(fi)
          this%livestemp_to_litter_patch(i)              = value_patch
          this%grainp_to_food_patch(i)                   = value_patch
          this%grainp_xfer_to_grainp_patch(i)            = value_patch
          this%ppool_to_grainp_patch(i)                  = value_patch
          this%ppool_to_grainp_storage_patch(i)          = value_patch
          this%grainp_storage_to_xfer_patch(i)           = value_patch
          this%frootp_to_retransp_patch(i)               = value_patch
       end do
    end if

    do j = 1, nlevdecomp_full
       do fi = 1,num_column
          i = filter_column(fi)

          ! phenology: litterfall and crop fluxes associated wit
          this%phenology_p_to_litr_met_p_col(i,j)        = value_column
          this%phenology_p_to_litr_cel_p_col(i,j)        = value_column
          this%phenology_p_to_litr_lig_p_col(i,j)        = value_column

          ! gap mortality
          this%gap_mortality_p_to_litr_met_p_col(i,j)    = value_column
          this%gap_mortality_p_to_litr_cel_p_col(i,j)    = value_column
          this%gap_mortality_p_to_litr_lig_p_col(i,j)    = value_column
          this%gap_mortality_p_to_cwdp_col(i,j)          = value_column

          ! fire
          this%fire_mortality_p_to_cwdp_col(i,j)         = value_column
          this%m_p_to_litr_met_fire_col(i,j)             = value_column
          this%m_p_to_litr_cel_fire_col(i,j)             = value_column
          this%m_p_to_litr_lig_fire_col(i,j)             = value_column

          ! harvest
          this%harvest_p_to_litr_met_p_col(i,j)          = value_column
          this%harvest_p_to_litr_cel_p_col(i,j)          = value_column
          this%harvest_p_to_litr_lig_p_col(i,j)          = value_column
          this%harvest_p_to_cwdp_col(i,j)                = value_column

          this%primp_to_labilep_vr_col(i,j)              = value_column
          this%labilep_to_secondp_vr_col(i,j)            = value_column
          this%secondp_to_labilep_vr_col(i,j)            = value_column
          this%secondp_to_occlp_vr_col(i,j)              = value_column

          this%sminp_leached_vr_col(i,j)                 = value_column

          this%potential_immob_p_vr_col(i,j)             = value_column
          this%actual_immob_p_vr_col(i,j)                = value_column
          this%sminp_to_plant_vr_col(i,j)                = value_column
          this%supplement_to_sminp_vr_col(i,j)           = value_column
          this%gross_pmin_vr_col(i,j)                    = value_column
          this%net_pmin_vr_col(i,j)                      = value_column
          this%biochem_pmin_vr_col(i,j)                  = value_column

          ! bgc interface & pflotran
          this%plant_pdemand_vr_col(i,j)                 = value_column

          this%adsorb_to_labilep_vr(i,j)                 = value_column
          this%desorb_to_solutionp_vr(i,j)               = value_column

       end do
    end do

    do fi = 1,num_column
       i = filter_column(fi)

       this%pdep_to_sminp_col(i)             = value_column
       this%fert_p_to_sminp_col(i)           = value_column
       this%hrv_deadstemp_to_prod10p_col(i)  = value_column
       this%hrv_deadstemp_to_prod100p_col(i) = value_column
       this%hrv_cropp_to_prod1p_col(i)       = value_column
       this%prod10p_loss_col(i)              = value_column
       this%prod100p_loss_col(i)             = value_column
       this%product_ploss_col(i)             = value_column
       this%prod1p_loss_col(i)               = value_column
       this%potential_immob_p_col(i)         = value_column
       this%actual_immob_p_col(i)            = value_column
       this%sminp_to_plant_col(i)            = value_column
       this%supplement_to_sminp_col(i)       = value_column
       this%gross_pmin_col(i)                = value_column
       this%net_pmin_col(i)                  = value_column
       this%biochem_pmin_col(i)              = value_column
       this%primp_to_labilep_col(i)          = value_column
       this%labilep_to_secondp_col(i)        = value_column
       this%secondp_to_labilep_col(i)        = value_column
       this%secondp_to_occlp_col(i)          = value_column
       this%sminp_leached_col(i)             = value_column
       this%pinputs_col(i)                   = value_column
       this%poutputs_col(i)                  = value_column
       this%fire_ploss_col(i)                = value_column
       this%som_p_leached_col(i)             = value_column

       ! Zero p2c column fluxes
       this%fire_ploss_col(i) = value_column
       this%wood_harvestp_col(i) = value_column

       !! bgc-interface
       this%plant_pdemand_col(i) = value_column

       this%fire_ploss_col(i)                = value_column
       this%wood_harvestp_col(i)             = value_column

       this%adsorb_to_labilep_col(i)         = value_column
       this%desorb_to_solutionp_col(i)       = value_column

    end do

    do k = 1, ndecomp_pools
       do fi = 1,num_column
          i = filter_column(fi)
          this%decomp_ppools_leached_col(i,k) = value_column
          this%m_decomp_ppools_to_fire_col(i,k) = value_column
       end do
    end do

    do k = 1, ndecomp_pools
       do j = 1, nlevdecomp_full
          do fi = 1,num_column
             i = filter_column(fi)
             this%m_decomp_ppools_to_fire_vr_col(i,j,k) = value_column
             this%decomp_ppools_transport_tendency_col(i,j,k) = value_column
          end do
       end do
    end do

    do l = 1, ndecomp_cascade_transitions
       do fi = 1,num_column
          i = filter_column(fi)
          this%decomp_cascade_ptransfer_col(i,l) = value_column
          this%decomp_cascade_sminp_flux_col(i,l) = value_column
       end do
    end do

    do l = 1, ndecomp_cascade_transitions
       do j = 1, nlevdecomp_full
          do fi = 1,num_column
             i = filter_column(fi)
             this%decomp_cascade_ptransfer_vr_col(i,j,l) = value_column
             this%decomp_cascade_sminp_flux_vr_col(i,j,l) = value_column
          end do
       end do
    end do

    do k = 1, ndecomp_pools
       do j = 1, nlevdecomp_full
          do fi = 1,num_column
             i = filter_column(fi)
             this%decomp_ppools_sourcesink_col(i,j,k) = value_column
          end do
       end do
    end do

    ! pflotran
    do k = 1, ndecomp_pools
       do j = 1, nlevdecomp_full
          do fi = 1,num_column
             i = filter_column(fi)
             ! only initializing in the first time-step
             if ( this%externalp_to_decomp_ppools_col(i,j,k) == spval ) then
                this%externalp_to_decomp_ppools_col(i,j,k) = value_column
             end if
          end do
       end do
    end do

    do j = 1, nlevdecomp_full
       do fi = 1,num_column
          i = filter_column(fi)
          ! only initializing in the first time-step
          if ( this%sminp_net_transport_vr_col(i,j) == spval ) then
             this%sminp_net_transport_vr_col(i,j) = value_column
          end if
       end do
    end do

    do fi = 1,num_column
       i = filter_column(fi)
       ! only initializing in the first time-step
       if ( this%externalp_to_decomp_delta_col(i) == spval ) then
          this%externalp_to_decomp_delta_col(i) = value_column
       end if
       if ( this%sminp_net_transport_delta_col(i) == spval ) then
          this%sminp_net_transport_delta_col(i)   = value_column
       end if
    end do
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
    integer  :: c, j          ! indices
    !-----------------------------------------------------------------------

    do c = bounds%begc,bounds%endc
       this%dwt_seedp_to_leaf_col(c)     = 0._r8
       this%dwt_seedp_to_deadstem_col(c) = 0._r8
       this%dwt_conv_pflux_col(c)        = 0._r8
       this%dwt_prod10p_gain_col(c)      = 0._r8
       this%dwt_prod100p_gain_col(c)     = 0._r8
    end do

    do j = 1, nlevdecomp_full
       do c = bounds%begc,bounds%endc
          this%dwt_frootp_to_litr_met_p_col(c,j) = 0._r8
          this%dwt_frootp_to_litr_cel_p_col(c,j) = 0._r8
          this%dwt_frootp_to_litr_lig_p_col(c,j) = 0._r8
          this%dwt_livecrootp_to_cwdp_col(c,j)   = 0._r8
          this%dwt_deadcrootp_to_cwdp_col(c,j)   = 0._r8
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
     use pftvarcon     , only: npcropmin

     ! !ARGUMENTS:
     class (phosphorusflux_type) :: this
     type(bounds_type) , intent(in) :: bounds
     integer           , intent(in) :: num_soilc       ! number of soil columns in filter
     integer           , intent(in) :: filter_soilc(:) ! filter for soil columns
     integer           , intent(in) :: num_soilp       ! number of soil patches in filter
     integer           , intent(in) :: filter_soilp(:) ! filter for soil patches

     integer :: fc, c, j
     ! total column-level fire P losses
     do fc = 1,num_soilc
        c = filter_soilc(fc)
        this%fire_ploss_col(c) = this%fire_ploss_p2c_col(c)  + this%fire_decomp_ploss_col(c)
     end do
     ! vertically integrate inorganic P flux
     do j = 1, nlevdecomp
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           this%primp_to_labilep_col(c) = &
                this%primp_to_labilep_col(c) + &
                this%primp_to_labilep_vr_col(c,j) * dzsoi_decomp(j)
        end do
     end do
     
     do fc = 1,num_soilc
        c = filter_soilc(fc)

        ! column-level P losses due to landcover change
        this%dwt_ploss_col(c) = &
             this%dwt_conv_pflux_col(c)

        ! total wood product P loss
        this%product_ploss_col(c) = &
             this%prod1p_loss_col(c) + &
             this%prod10p_loss_col(c) + &
             this%prod100p_loss_col(c)
     end do

     ! supplementary P supplement_to_sminp
     do j = 1, nlevdecomp
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           this%supplement_to_sminp_col(c) = &
                this%supplement_to_sminp_col(c) + &
                this%supplement_to_sminp_vr_col(c,j) * dzsoi_decomp(j)
        end do
     end do

  end subroutine Summary_betr
 !-----------------------------------------------------------------------
  subroutine Summary(this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
    !
    ! !USES:
    use clm_varpar    , only: nlevdecomp,ndecomp_cascade_transitions,ndecomp_pools
    use clm_varctl    , only: use_nitrif_denitrif
    use subgridAveMod , only: p2c
    use pftvarcon     , only: npcropmin
    use tracer_varcon  , only : is_active_betr_bgc
    ! pflotran
!    use clm_varctl    , only: use_pflotran, pf_cmode
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

    do fp = 1,num_soilp
       p = filter_soilp(fp)

       ! total P deployment (from sminn and retranslocated P pool) (PDEPLOY)
       this%pdeploy_patch(p) = &
            this%sminp_to_ppool_patch(p) + &
            this%retransp_to_ppool_patch(p)

       ! pft-level wood harvest
       this%wood_harvestp_patch(p) = &
            this%hrv_deadstemp_to_prod10p_patch(p) + &
            this%hrv_deadstemp_to_prod100p_patch(p)
       if ( crop_prog .and. pft%itype(p) >= npcropmin )then
            this%wood_harvestp_patch(p) = &
            this%wood_harvestp_patch(p) + &
            this%hrv_cropp_to_prod1p_patch(p)
       end if

       ! total pft-level fire P losses
       this%fire_ploss_patch(p) = &
            this%m_leafp_to_fire_patch(p)               + &
            this%m_leafp_storage_to_fire_patch(p)       + &
            this%m_leafp_xfer_to_fire_patch(p)          + &
            this%m_frootp_to_fire_patch(p)              + &
            this%m_frootp_storage_to_fire_patch(p)      + &
            this%m_frootp_xfer_to_fire_patch(p)         + &
            this%m_livestemp_to_fire_patch(p)           + &
            this%m_livestemp_storage_to_fire_patch(p)   + &
            this%m_livestemp_xfer_to_fire_patch(p)      + &
            this%m_deadstemp_to_fire_patch(p)           + &
            this%m_deadstemp_storage_to_fire_patch(p)   + &
            this%m_deadstemp_xfer_to_fire_patch(p)      + &
            this%m_livecrootp_to_fire_patch(p)          + &
            this%m_livecrootp_storage_to_fire_patch(p)  + &
            this%m_livecrootp_xfer_to_fire_patch(p)     + &
            this%m_deadcrootp_to_fire_patch(p)          + &
            this%m_deadcrootp_storage_to_fire_patch(p)  + &
            this%m_deadcrootp_xfer_to_fire_patch(p)     + &
            this%m_retransp_to_fire_patch(p)

    end do


    call p2c(bounds, num_soilc, filter_soilc, &
         this%fire_ploss_patch(bounds%begp:bounds%endp), &
         this%fire_ploss_p2c_col(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%wood_harvestp_patch(bounds%begp:bounds%endp), &
         this%wood_harvestp_col(bounds%begc:bounds%endc))

    if(is_active_betr_bgc)then
      call this%Summary_betr(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
      return
    endif
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%supplement_to_sminp_col(c) = 0._r8
       this%som_p_leached_col(c)       = 0._r8
    end do

    ! pflotran
    !----------------------------------------------------------------
    if (.not.(use_pflotran .and. pf_cmode)) then
    ! vertically integrate decomposing P cascade fluxes and soil mineral P fluxes associated with decomposition cascade
    do k = 1, ndecomp_cascade_transitions
       do j = 1,nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)

             this%decomp_cascade_ptransfer_col(c,k) = &
                  this%decomp_cascade_ptransfer_col(c,k) + &
                  this%decomp_cascade_ptransfer_vr_col(c,j,k) * dzsoi_decomp(j)

             this%decomp_cascade_sminp_flux_col(c,k) = &
                  this%decomp_cascade_sminp_flux_col(c,k) + &
                  this%decomp_cascade_sminp_flux_vr_col(c,j,k) * dzsoi_decomp(j)
          end do
       end do
    end do
    end if !!if (.not.(use_pflotran .and. pf_cmode))
    !-----------------------------------------------------------------

    ! vertically integrate inorganic P flux
    do j = 1, nlevdecomp
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%primp_to_labilep_col(c) = &
               this%primp_to_labilep_col(c) + &
               this%primp_to_labilep_vr_col(c,j) * dzsoi_decomp(j)
       end do
    end do

    do j = 1, nlevdecomp
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%labilep_to_secondp_col(c) = &
               this%labilep_to_secondp_col(c) + &
               this%labilep_to_secondp_vr_col(c,j) * dzsoi_decomp(j)
       end do
    end do

    do j = 1, nlevdecomp
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%secondp_to_labilep_col(c) = &
               this%secondp_to_labilep_col(c) + &
               this%secondp_to_labilep_vr_col(c,j) * dzsoi_decomp(j)
       end do
    end do

    do j = 1, nlevdecomp
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%secondp_to_occlp_col(c) = &
               this%secondp_to_occlp_col(c) + &
               this%secondp_to_occlp_vr_col(c,j) * dzsoi_decomp(j)
       end do
    end do


    ! vertically integrate leaching flux
    do j = 1, nlevdecomp
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%sminp_leached_col(c) = &
               this%sminp_leached_col(c) + &
               this%sminp_leached_vr_col(c,j) * dzsoi_decomp(j)
       end do
    end do

    ! vertically integrate column-level fire P losses
    do k = 1, ndecomp_pools
       do j = 1, nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%m_decomp_ppools_to_fire_col(c,k) = &
                  this%m_decomp_ppools_to_fire_col(c,k) + &
                  this%m_decomp_ppools_to_fire_vr_col(c,j,k) * dzsoi_decomp(j)
          end do
       end do
    end do

    ! total column-level fire P losses
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%fire_ploss_col(c) = this%fire_ploss_p2c_col(c)
    end do
    do k = 1, ndecomp_pools
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%fire_ploss_col(c) = &
               this%fire_ploss_col(c) + &
               this%m_decomp_ppools_to_fire_col(c,k)
       end do
    end do

    ! supplementary P supplement_to_sminp
    do j = 1, nlevdecomp
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%supplement_to_sminp_col(c) = &
               this%supplement_to_sminp_col(c) + &
               this%supplement_to_sminp_vr_col(c,j) * dzsoi_decomp(j)
       end do
    end do

    do fc = 1,num_soilc
       c = filter_soilc(fc)

       ! column-level P losses due to landcover change
       this%dwt_ploss_col(c) = &
            this%dwt_conv_pflux_col(c)


       ! total wood product P loss
       this%product_ploss_col(c) = &
            this%prod10p_loss_col(c) + &
            this%prod100p_loss_col(c) + &
            this%prod1p_loss_col(c)
    end do

    ! add up all vertical transport tendency terms and calculate total som leaching loss as the sum of these
    do l = 1, ndecomp_pools
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%decomp_ppools_leached_col(c,l) = 0._r8
       end do

       do j = 1, nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%decomp_ppools_leached_col(c,l) = &
                  this%decomp_ppools_leached_col(c,l) + &
                  this%decomp_ppools_transport_tendency_col(c,j,l) * dzsoi_decomp(j)
          end do
       end do

       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%som_p_leached_col(c) = &
               this%som_p_leached_col(c) + &
               this%decomp_ppools_leached_col(c,l)
       end do
    end do

    do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%biochem_pmin_col(c) = 0.0_r8
    end do

    do j = 1, nlevdecomp
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%biochem_pmin_col(c) = this%biochem_pmin_col(c) + &
               this%biochem_pmin_vr_col(c,j)* dzsoi_decomp(j)
       end do
    end do

    do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%adsorb_to_labilep_col(c) = 0._r8
       this%desorb_to_solutionp_col(c) = 0._r8
    end do

    do j = 1, nlevdecomp
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%adsorb_to_labilep_col(c) = this%adsorb_to_labilep_col(c) + &
               this%adsorb_to_labilep_vr(c,j)* dzsoi_decomp(j)
          this%desorb_to_solutionp_col(c) = this%desorb_to_solutionp_col(c) + &
               this%desorb_to_solutionp_vr(c,j)* dzsoi_decomp(j)
       end do
    end do

    do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%actual_immob_p_col(c) = 0._r8
       this%smin_p_to_plant_col(c) = 0._r8
       this%plant_to_litter_pflux(c) = 0._r8
       this%plant_to_cwd_pflux(c) = 0._r8
       do j = 1, nlevdecomp
          this%actual_immob_p_col(c)= this%actual_immob_p_col(c) + &
               this%actual_immob_p_vr_col(c,j) * dzsoi_decomp(j)
          this%smin_p_to_plant_col(c)= this%smin_p_to_plant_col(c) + &
               this%sminp_to_plant_vr_col(c,j) * dzsoi_decomp(j)
          this%plant_to_litter_pflux(c) = &
               this%plant_to_litter_pflux(c)  + &
               this%phenology_p_to_litr_met_p_col(c,j)* dzsoi_decomp(j) + &
               this%phenology_p_to_litr_cel_p_col(c,j)* dzsoi_decomp(j) + &
               this%phenology_p_to_litr_lig_p_col(c,j)* dzsoi_decomp(j) + &
               this%gap_mortality_p_to_litr_met_p_col(c,j)* dzsoi_decomp(j) + &
               this%gap_mortality_p_to_litr_cel_p_col(c,j)* dzsoi_decomp(j) + &
               this%gap_mortality_p_to_litr_lig_p_col(c,j)* dzsoi_decomp(j) + &
               this%m_p_to_litr_met_fire_col(c,j)* dzsoi_decomp(j) + &
               this%m_p_to_litr_cel_fire_col(c,j)* dzsoi_decomp(j) + &
               this%m_p_to_litr_lig_fire_col(c,j)* dzsoi_decomp(j)
          this%plant_to_cwd_pflux(c) = &
               this%plant_to_cwd_pflux(c) + &
               this%gap_mortality_p_to_cwdp_col(c,j)* dzsoi_decomp(j) + &
               this%fire_mortality_p_to_cwdp_col(c,j)* dzsoi_decomp(j)
       end do
    end do

    !! bgc interface & pflotran:
    !----------------------------------------------------------------
    if (use_bgc_interface) then
        call PSummary_interface(this, bounds, num_soilc, filter_soilc)
    end if
    !----------------------------------------------------------------

  end subroutine Summary

!!-------------------------------------------------------------------------------------------------
! !INTERFACE:
subroutine PSummary_interface(this,bounds,num_soilc, filter_soilc)
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
   class (phosphorusflux_type)     :: this
   type(bounds_type) ,  intent(in) :: bounds
   integer,             intent(in) :: num_soilc       ! number of soil columns in filter
   integer,             intent(in) :: filter_soilc(:) ! filter for soil columns
!
! !REVISION HISTORY:
!!08/26/2015: created by Gangsheng Wang
!
! !LOCAL VARIABLES:
   integer :: c,j, l      ! indices
   integer :: fc          ! column filter indices
   real(r8):: dtime             ! radiation time step (seconds)

   ! set time steps
    dtime = real( get_step_size(), r8 )
    if (use_pflotran .and. pf_cmode) then
        !! TODO
    end if

       ! summarize at column-level vertically-resolved littering/removal for PFLOTRAN bgc input needs
       ! first it needs to save the total column-level N rate btw plant pool and decomposible pools at previous time step
       ! for adjusting difference when doing balance check

       do fc = 1,num_soilc
         c = filter_soilc(fc)
         this%externalp_to_decomp_delta_col(c) = 0._r8
         this%sminp_net_transport_delta_col(c)   = 0._r8
         do j = 1, nlevdecomp
            do l = 1, ndecomp_pools
               this%externalp_to_decomp_delta_col(c) =    &
                  this%externalp_to_decomp_delta_col(c) + &
                    this%externalp_to_decomp_ppools_col(c,j,l)*dzsoi_decomp(j)
            end do

            ! sminp leaching/runoff at previous time-step, which may be as source by PFLOTRAN
            this%sminp_net_transport_delta_col(c) = &
               this%sminp_net_transport_delta_col(c) + &
                    this%sminp_net_transport_vr_col(c,j)*dzsoi_decomp(j)

         end do
       end do

       ! do the initialization for the following 2 variables here.
       ! DON'T do so in the beginning of CLM-CN time-step (otherwise the above saved will not work)
       this%externalp_to_decomp_ppools_col(:,:,:) = 0._r8
       this%sminp_net_transport_vr_col(:,:) = 0._r8

       ! add up all vertically-resolved addition/removal rates (gC/m3/s) of decomp_pools
       do l = 1, ndecomp_pools
          do j = 1, nlevdecomp
             do fc = 1,num_soilc
                c = filter_soilc(fc)

                ! for litter C pools
                if (l==i_met_lit) then
                   this%externalp_to_decomp_ppools_col(c,j,l) =              &
                       this%externalp_to_decomp_ppools_col(c,j,l)            &
                        + this%phenology_p_to_litr_met_p_col(c,j)            &
                        + this%dwt_frootp_to_litr_met_p_col(c,j)             &
                        + this%gap_mortality_p_to_litr_met_p_col(c,j)        &
                        + this%harvest_p_to_litr_met_p_col(c,j)              !!&
!                        + this%m_p_to_litr_met_fire_col(c,j)                 &
!                        - this%m_decomp_ppools_to_fire_vr_col(c,j,l)

                elseif (l==i_cel_lit) then
                   this%externalp_to_decomp_ppools_col(c,j,l) =              &
                       this%externalp_to_decomp_ppools_col(c,j,l)            &
                        + this%phenology_p_to_litr_cel_p_col(c,j)            &
                        + this%dwt_frootp_to_litr_cel_p_col(c,j)             &
                        + this%gap_mortality_p_to_litr_cel_p_col(c,j)        &
                        + this%harvest_p_to_litr_cel_p_col(c,j)              !!&
!                        + this%m_p_to_litr_cel_fire_col(c,j)                 &
!                        - this%m_decomp_ppools_to_fire_vr_col(c,j,l)

                elseif (l==i_lig_lit) then
                   this%externalp_to_decomp_ppools_col(c,j,l) =              &
                       this%externalp_to_decomp_ppools_col(c,j,l)            &
                        + this%phenology_p_to_litr_lig_p_col(c,j)            &
                        + this%dwt_frootp_to_litr_lig_p_col(c,j)             &
                        + this%gap_mortality_p_to_litr_lig_p_col(c,j)        &
                        + this%harvest_p_to_litr_lig_p_col(c,j)              !!&
!                        + this%m_p_to_litr_lig_fire_col(c,j)                 &
!                        - this%m_decomp_ppools_to_fire_vr_col(c,j,l)

                ! for cwd
                elseif (l==i_cwd) then
                   this%externalp_to_decomp_ppools_col(c,j,l) =              &
                       this%externalp_to_decomp_ppools_col(c,j,l)            &
                        + this%dwt_livecrootp_to_cwdp_col(c,j)               &
                        + this%dwt_deadcrootp_to_cwdp_col(c,j)               &
                        + this%gap_mortality_p_to_cwdp_col(c,j)              &
                        + this%harvest_p_to_cwdp_col(c,j)                    !!&
!                        + this%fire_mortality_p_to_cwdp_col(c,j)
!
             ! for som n
!                else
!                   this%externalp_to_decomp_ppools_col(c,j,l) =              &
!                       this%externalp_to_decomp_ppools_col(c,j,l)            &
!                        - this%m_decomp_ppools_to_fire_vr_col(c,j,l)

                end if

             ! the following is the net changes of plant N to decompible N poools between time-step
             ! in pflotran, decomposible N pools increments ARE from previous time-step (saved above);
             ! while, in CLM-CN all plant N pools are updated with current N fluxes among plant and ground/soil.
             ! therefore, when do balance check it is needed to adjust the time-lag of changes.
                this%externalp_to_decomp_delta_col(c) =   &
                            this%externalp_to_decomp_delta_col(c) - &
                            this%externalp_to_decomp_ppools_col(c,j,l)*dzsoi_decomp(j)

                if (abs(this%externalp_to_decomp_ppools_col(c,j,l))<=1.e-21_r8) then
                    this%externalp_to_decomp_ppools_col(c,j,l) = 0._r8
                end if

             end do
          end do
       end do

       ! if pflotran hydrology NOT coupled, need to do adjusting for sminp leaching for balance error checking
       if (.not. pf_hmode) then
          do j = 1, nlevdecomp
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                this%sminp_net_transport_vr_col(c,j) = 0._r8
!                this%sminp_net_transport_vr_col(c,j) = this%sminp_leached_vr_col(c,j)
                this%sminp_net_transport_delta_col(c) = &
                            this%sminp_net_transport_delta_col(c) - &
                            this%sminp_net_transport_vr_col(c,j)*dzsoi_decomp(j)
             end do
          end do
       end if

       ! change the sign so that it is the increments from the previous time-step (unit: g/m2/s)
       do fc = 1, num_soilc
          c = filter_soilc(fc)
          this%externalp_to_decomp_delta_col(c)     = -this%externalp_to_decomp_delta_col(c)
          this%sminp_net_transport_delta_col(c)     = -this%sminp_net_transport_delta_col(c)
       end do
end subroutine PSummary_interface
!!-------------------------------------------------------------------------------------------------

end module PhosphorusFluxType
