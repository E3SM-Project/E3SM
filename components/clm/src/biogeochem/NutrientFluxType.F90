module NutrientFluxType

#include "shr_assert.h"

  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_infnan_mod         , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use decompMod              , only : bounds_type
  use clm_varcon             , only : spval, ispval
  use abortutils             , only : endrun
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private

  type, public :: nutrientflux_type

     ! dynamic landcover fluxes
     real(r8), pointer :: dwt_seed_to_deadstem_patch               (:)     !
     real(r8), pointer :: dwt_seed_to_leaf_patch                   (:)     !
     real(r8), pointer :: dwt_seed_to_deadstem_grc                 (:)     !
     real(r8), pointer :: dwt_seed_to_leaf_grc                     (:)     !
     real(r8), pointer :: dwt_slash_flux_col                       (:)     !
     real(r8), pointer :: dwt_loss_col                             (:)     !
     real(r8), pointer :: dwt_deadcroot_to_cwd_col                 (:,:)   !
     real(r8), pointer :: dwt_froot_to_litr_cel_col                (:,:)   !
     real(r8), pointer :: dwt_froot_to_litr_lig_col                (:,:)   !
     real(r8), pointer :: dwt_froot_to_litr_met_col                (:,:)   !
     real(r8), pointer :: dwt_livecroot_to_cwd_col                 (:,:)   !

     real(r8), pointer :: dwt_conv_flux_patch                      (:)     !
     real(r8), pointer :: dwt_conv_flux_col                        (:)     !
     real(r8), pointer :: dwt_conv_flux_grc                        (:)     !

     real(r8), pointer :: dwt_prod10_gain_patch                    (:)     !
     real(r8), pointer :: dwt_prod100_gain_patch                   (:)     !
     real(r8), pointer :: dwt_prod10_gain_col                      (:)     !
     real(r8), pointer :: dwt_prod100_gain_col                     (:)     !
     real(r8), pointer :: dwt_prod10_gain_grc                      (:)     !
     real(r8), pointer :: dwt_prod100_gain_grc                     (:)     !
     real(r8), pointer :: dwt_crop_product_gain_patch              (:)     !

     !wood product pool loss fluxes
     real(r8), pointer :: product_loss_col                         (:)     !

     ! annual turnover of storage to transfer pools
     real(r8), pointer :: deadcroot_storage_to_xfer_patch          (:)     !
     real(r8), pointer :: deadstem_storage_to_xfer_patch           (:)     !
     real(r8), pointer :: froot_storage_to_xfer_patch              (:)     !
     real(r8), pointer :: grain_storage_to_xfer_patch              (:)     !
     real(r8), pointer :: leaf_storage_to_xfer_patch               (:)     !
     real(r8), pointer :: livecroot_storage_to_xfer_patch          (:)     !
     real(r8), pointer :: livestem_storage_to_xfer_patch           (:)     !

     ! transport in soil column
     real(r8), pointer :: decomp_pools_leached_col                 (:,:)   !
     real(r8), pointer :: decomp_pools_sourcesink_col              (:,:,:) !
     real(r8), pointer :: decomp_pools_transport_tendency_col      (:,:,:) !

     ! leaf and fine root litterfall fluxes                          
     real(r8), pointer :: leaf_to_litter_patch                     (:)     !
     real(r8), pointer :: froot_to_litter_patch                    (:)     !
     real(r8), pointer :: livestem_to_litter_patch                 (:)     !
     real(r8), pointer :: grain_to_food_patch                      (:)     !

     ! gap mortality fluxes
     real(r8), pointer :: m_leaf_to_litter_patch                   (:)     !
     real(r8), pointer :: m_leaf_storage_to_litter_patch           (:)     !
     real(r8), pointer :: m_leaf_xfer_to_litter_patch              (:)     !
     real(r8), pointer :: m_froot_to_litter_patch                  (:)     !
     real(r8), pointer :: m_froot_storage_to_litter_patch          (:)     !
     real(r8), pointer :: m_froot_xfer_to_litter_patch             (:)     !
     real(r8), pointer :: m_livestem_to_litter_patch               (:)     !
     real(r8), pointer :: m_livestem_storage_to_litter_patch       (:)     !
     real(r8), pointer :: m_livestem_xfer_to_litter_patch          (:)     !
     real(r8), pointer :: m_deadstem_storage_to_litter_patch       (:)     !
     real(r8), pointer :: m_deadcroot_xfer_to_litter_patch         (:)     !
     real(r8), pointer :: m_deadstem_to_litter_patch               (:)     !
     real(r8), pointer :: m_livecroot_storage_to_litter_patch      (:)     !
     real(r8), pointer :: m_livecroot_to_litter_patch              (:)     !
     real(r8), pointer :: m_livecroot_xfer_to_litter_patch         (:)     !
     real(r8), pointer :: m_deadcroot_to_litter_patch              (:)     !
     real(r8), pointer :: m_deadcroot_storage_to_litter_patch      (:)     !
     real(r8), pointer :: m_deadstem_xfer_to_litter_patch          (:)     !
     real(r8), pointer :: m_pool_to_litter_patch                   (:)     !
     ! gap mortality fluxes
     real(r8), pointer :: gap_mortality_to_litr_met_col            (:,:)   !
     real(r8), pointer :: gap_mortality_to_litr_cel_col            (:,:)   !
     real(r8), pointer :: gap_mortality_to_litr_lig_col            (:,:)   !
     real(r8), pointer :: gap_mortality_to_cwd_col                 (:,:)   !

     ! harvest mortality fluxes
     real(r8), pointer :: hrv_leaf_storage_to_litter_patch         (:)     !
     real(r8), pointer :: hrv_leaf_to_litter_patch                 (:)     !
     real(r8), pointer :: hrv_leaf_xfer_to_litter_patch            (:)     !
     real(r8), pointer :: hrv_froot_storage_to_litter_patch        (:)     !
     real(r8), pointer :: hrv_froot_to_litter_patch                (:)     !
     real(r8), pointer :: hrv_froot_xfer_to_litter_patch           (:)     !
     real(r8), pointer :: hrv_livestem_storage_to_litter_patch     (:)     !
     real(r8), pointer :: hrv_livestem_to_litter_patch             (:)     !
     real(r8), pointer :: hrv_livestem_xfer_to_litter_patch        (:)     !
     real(r8), pointer :: hrv_deadstem_to_prod10_patch             (:)     !
     real(r8), pointer :: hrv_deadstem_to_prod100_patch            (:)     !
     real(r8), pointer :: hrv_deadstem_storage_to_litter_patch     (:)     !
     real(r8), pointer :: hrv_deadstem_xfer_to_litter_patch        (:)     !
     real(r8), pointer :: hrv_deadcroot_storage_to_litter_patch    (:)     !
     real(r8), pointer :: hrv_deadcroot_to_litter_patch            (:)     !
     real(r8), pointer :: hrv_deadcroot_xfer_to_litter_patch       (:)     !
     real(r8), pointer :: hrv_livecroot_storage_to_litter_patch    (:)     !
     real(r8), pointer :: hrv_livecroot_to_litter_patch            (:)     !
     real(r8), pointer :: hrv_livecroot_xfer_to_litter_patch       (:)     !
     real(r8), pointer :: hrv_pool_to_litter_patch                 (:)     !
     real(r8), pointer :: hrv_leaf_to_prod1_patch                  (:)     !
     real(r8), pointer :: hrv_livestem_to_prod1_patch              (:)     !
     real(r8), pointer :: hrv_grain_to_prod1_patch                 (:)     !
     real(r8), pointer :: hrv_crop_to_prod1_patch                  (:)     !
     real(r8), pointer :: harvest_to_cwd_col                       (:,:)   !
     real(r8), pointer :: harvest_to_litr_cel_col                  (:,:)   !
     real(r8), pointer :: harvest_to_litr_lig_col                  (:,:)   !
     real(r8), pointer :: harvest_to_litr_met_col                  (:,:)   !

     ! fire fluxes      
     real(r8), pointer :: m_pool_to_fire_patch                     (:)     !
     real(r8), pointer :: m_deadcroot_storage_to_fire_patch        (:)     !
     real(r8), pointer :: m_deadcroot_to_fire_patch                (:)     !
     real(r8), pointer :: m_deadcroot_xfer_to_fire_patch           (:)     !
     real(r8), pointer :: m_deadstem_storage_to_fire_patch         (:)     !
     real(r8), pointer :: m_deadstem_to_fire_patch                 (:)     !
     real(r8), pointer :: m_deadstem_xfer_to_fire_patch            (:)     !
     real(r8), pointer :: m_froot_storage_to_fire_patch            (:)     !
     real(r8), pointer :: m_froot_to_fire_patch                    (:)     !
     real(r8), pointer :: m_froot_xfer_to_fire_patch               (:)     !
     real(r8), pointer :: m_leaf_storage_to_fire_patch             (:)     !
     real(r8), pointer :: m_leaf_to_fire_patch                     (:)     !
     real(r8), pointer :: m_leaf_xfer_to_fire_patch                (:)     !
     real(r8), pointer :: m_livecroot_storage_to_fire_patch        (:)     !
     real(r8), pointer :: m_livecroot_to_fire_patch                (:)     !
     real(r8), pointer :: m_livecroot_xfer_to_fire_patch           (:)     !
     real(r8), pointer :: m_livestem_storage_to_fire_patch         (:)     !
     real(r8), pointer :: m_livestem_to_fire_patch                 (:)     !
     real(r8), pointer :: m_livestem_xfer_to_fire_patch            (:)     !
     real(r8), pointer :: m_livecroot_to_deadcroot_fire_patch      (:)     !
     real(r8), pointer :: m_livestem_to_deadstem_fire_patch        (:)     !
     real(r8), pointer :: fire_mortality_to_cwd_col                (:,:)   !

     real(r8), pointer :: m_pool_to_litter_fire_patch              (:)     !
     real(r8), pointer :: m_deadcroot_storage_to_litter_fire_patch (:)     !
     real(r8), pointer :: m_deadcroot_to_litter_fire_patch         (:)     !
     real(r8), pointer :: m_deadcroot_xfer_to_litter_fire_patch    (:)     !
     real(r8), pointer :: m_deadstem_storage_to_litter_fire_patch  (:)     !
     real(r8), pointer :: m_deadstem_to_litter_fire_patch          (:)     !
     real(r8), pointer :: m_deadstem_xfer_to_litter_fire_patch     (:)     !
     real(r8), pointer :: m_froot_storage_to_litter_fire_patch     (:)     !
     real(r8), pointer :: m_froot_to_litter_fire_patch             (:)     !
     real(r8), pointer :: m_froot_xfer_to_litter_fire_patch        (:)     !
     real(r8), pointer :: m_leaf_storage_to_litter_fire_patch      (:)     !
     real(r8), pointer :: m_leaf_to_litter_fire_patch              (:)     !
     real(r8), pointer :: m_leaf_xfer_to_litter_fire_patch         (:)     !
     real(r8), pointer :: m_livecroot_storage_to_litter_fire_patch (:)     !
     real(r8), pointer :: m_livecroot_to_litter_fire_patch         (:)     !
     real(r8), pointer :: m_livecroot_xfer_to_litter_fire_patch    (:)     !
     real(r8), pointer :: m_livestem_storage_to_litter_fire_patch  (:)     !
     real(r8), pointer :: m_livestem_to_litter_fire_patch          (:)     !
     real(r8), pointer :: m_livestem_xfer_to_litter_fire_patch     (:)     !

     real(r8), pointer :: m_to_litr_cel_fire_col                   (:,:)   !
     real(r8), pointer :: m_to_litr_lig_fire_col                   (:,:)   !
     real(r8), pointer :: m_to_litr_met_fire_col                   (:,:)   !

     real(r8), pointer :: fire_loss_col                            (:)     !
     real(r8), pointer :: fire_loss_patch                          (:)     !

     ! wood product pool losses
     real(r8), pointer :: prod100_loss_col                         (:)     !
     real(r8), pointer :: prod10_loss_col                          (:)     !
     real(r8), pointer :: prod1_loss_col                           (:)     !

     ! diagnostic
     real(r8), pointer :: inputs_patch                             (:)     !
     real(r8), pointer :: outputs_patch                            (:)     !
     real(r8), pointer :: plant_alloc_patch                        (:)     !
     real(r8), pointer :: som_leached_col                          (:)     !
     real(r8), pointer :: wood_harvest_col                         (:)     !
     real(r8), pointer :: wood_harvest_patch                       (:)     !

     ! allocation fluxes, from current GPP                     
     real(r8), pointer :: pool_to_deadcroot_patch                  (:)     !
     real(r8), pointer :: pool_to_deadcroot_storage_patch          (:)     !
     real(r8), pointer :: pool_to_deadstem_patch                   (:)     !
     real(r8), pointer :: pool_to_deadstem_storage_patch           (:)     !
     real(r8), pointer :: pool_to_froot_patch                      (:)     !
     real(r8), pointer :: pool_to_froot_storage_patch              (:)     !
     real(r8), pointer :: pool_to_grain_patch                      (:)     !
     real(r8), pointer :: pool_to_grain_storage_patch              (:)     !
     real(r8), pointer :: pool_to_leaf_patch                       (:)     !
     real(r8), pointer :: pool_to_leaf_storage_patch               (:)     !
     real(r8), pointer :: pool_to_livecroot_patch                  (:)     !
     real(r8), pointer :: pool_to_livecroot_storage_patch          (:)     !
     real(r8), pointer :: pool_to_livestem_patch                   (:)     !
     real(r8), pointer :: pool_to_livestem_storage_patch           (:)     !

     ! phenology fluxes from transfer pools                     
     real(r8), pointer :: phenology_to_litr_met_col                (:,:)   !
     real(r8), pointer :: phenology_to_litr_cel_col                (:,:)   !
     real(r8), pointer :: phenology_to_litr_lig_col                (:,:)   !

     real(r8), pointer :: livecroot_to_deadcroot_patch             (:)     !
     real(r8), pointer :: deadcroot_xfer_to_deadcroot_patch        (:)     !
     real(r8), pointer :: deadstem_xfer_to_deadstem_patch          (:)     !
     real(r8), pointer :: livecroot_xfer_to_livecroot_patch        (:)     !
     real(r8), pointer :: froot_xfer_to_froot_patch                (:)     !
     real(r8), pointer :: grain_xfer_to_grain_patch                (:)     !
     real(r8), pointer :: leaf_xfer_to_leaf_patch                  (:)     !

     ! turnover of livewood to deadwood
     real(r8), pointer :: livestem_to_deadstem_patch               (:)     !
     real(r8), pointer :: livestem_xfer_to_livestem_patch          (:)     !

  end type nutrientflux_type

end module NutrientFluxType
