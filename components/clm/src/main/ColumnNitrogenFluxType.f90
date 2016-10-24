module ColumnNitrogenFluxType


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
     real(r8), pointer :: hrv_cropn_to_prod1n_col                   (:)     ! crop N harvest mortality to 1-yr product pool (gN/m2/s)

     ! fire N fluxes 
     real(r8), pointer :: m_decomp_npools_to_fire_vr_col            (:,:,:) ! col vertically-resolved decomposing N fire loss (gN/m3/s)
     real(r8), pointer :: m_decomp_npools_to_fire_col               (:,:)   ! col vertically-integrated (diagnostic) decomposing N fire loss (gN/m2/s)
     real(r8), pointer :: fire_nloss_col                            (:)     ! col total column-level fire N loss (gN/m2/s)
     real(r8), pointer :: fire_nloss_p2c_col                        (:)     ! col patch2col column-level fire N loss (gN/m2/s) (p2c)
     real(r8), pointer :: fire_mortality_n_to_cwdn_col              (:,:)   ! col N fluxes associated with fire mortality to CWD pool (gN/m3/s)
     ! summary (diagnostic) flux variables, not involved in mass balance
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


    ! --- KB Are these column level variable??? -------!
	   ! real(r8), pointer :: col_plant_ndemand_vr                      (:,:)   ! column-level plant N demand
    !  real(r8), pointer :: col_plant_nh4demand_vr                    (:,:)   ! column-level plant NH4 demand
    !  real(r8), pointer :: col_plant_no3demand_vr                    (:,:)   ! column-level plant NO3 demand
    !  real(r8), pointer :: col_plant_pdemand_vr                      (:,:)   ! column-level plant P demand

    !  real(r8), pointer :: pmnf_decomp_cascade                       (:,:,:) !potential mineral N flux, from one pool to another

    !  real(r8), pointer :: plant_n_uptake_flux                       (:)     ! for the purpose of mass balance check  
    !  real(r8), pointer :: soil_n_immob_flux                         (:)     ! for the purpose of mass balance check
    !  real(r8), pointer :: soil_n_immob_flux_vr                      (:,:)   ! for the purpose of mass balance check
    !  real(r8), pointer :: soil_n_grossmin_flux                      (:)     ! for the purpose of mass balance check
    !  real(r8), pointer :: plant_to_litter_nflux                     (:)     ! for the purpose of mass balance check
    !  real(r8), pointer :: plant_to_cwd_nflux                        (:)     ! for the purpose of mass balance check
    ! -----------------------------------------------------------------------------

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

end module ColumnNitrogenFluxType