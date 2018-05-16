module NutrientFluxType

#include "shr_assert.h"

  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_infnan_mod         , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use decompMod              , only : bounds_type
  use clm_varcon             , only : spval, ispval
  use abortutils             , only : endrun
  use clm_varpar             , only : ndecomp_pools, nlevdecomp_full
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

     real(r8), pointer :: prev_leaf_to_litter_patch                (:)     !
     real(r8), pointer :: prev_froot_to_litter_patch               (:)     !
     
  end type nutrientflux_type

  public :: NutrientFluxInitAllocate

contains
  
  !------------------------------------------------------------------------
  subroutine NutrientFluxInitAllocate(nutrient_flux, bounds)
    !
    implicit none
    !
    class (nutrientflux_type)    :: nutrient_flux
    type(bounds_type), intent(in) :: bounds
    !
    integer           :: begp,endp
    integer           :: begc,endc
    integer           :: begg,endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc
    begg = bounds%begg; endg = bounds%endg

    allocate(nutrient_flux%dwt_seed_to_deadstem_patch               (begp:endp                                   )) 
    allocate(nutrient_flux%dwt_seed_to_leaf_patch                   (begp:endp                                   )) 
    allocate(nutrient_flux%dwt_seed_to_deadstem_grc                 (begg:endg                                   )) 
    allocate(nutrient_flux%dwt_seed_to_leaf_grc                     (begg:endg                                   )) 
    allocate(nutrient_flux%dwt_slash_flux_col                       (begc:endc                                   )) 
    allocate(nutrient_flux%dwt_loss_col                             (begc:endc                                   )) 
    allocate(nutrient_flux%dwt_deadcroot_to_cwd_col                 (begc:endc,1:nlevdecomp_full                 )) 
    allocate(nutrient_flux%dwt_froot_to_litr_cel_col                (begc:endc,1:nlevdecomp_full                 )) 
    allocate(nutrient_flux%dwt_froot_to_litr_lig_col                (begc:endc,1:nlevdecomp_full                 )) 
    allocate(nutrient_flux%dwt_froot_to_litr_met_col                (begc:endc,1:nlevdecomp_full                 )) 
    allocate(nutrient_flux%dwt_livecroot_to_cwd_col                 (begc:endc,1:nlevdecomp_full                 )) 
    allocate(nutrient_flux%dwt_conv_flux_patch                      (begp:endp                                   )) 
    allocate(nutrient_flux%dwt_conv_flux_col                        (begc:endc                                   )) 
    allocate(nutrient_flux%dwt_conv_flux_grc                        (begg:endg                                   )) 
    allocate(nutrient_flux%dwt_prod10_gain_patch                    (begp:endp                                   )) 
    allocate(nutrient_flux%dwt_prod100_gain_patch                   (begp:endp                                   )) 
    allocate(nutrient_flux%dwt_prod10_gain_col                      (begc:endc                                   )) 
    allocate(nutrient_flux%dwt_prod100_gain_col                     (begc:endc                                   )) 
    allocate(nutrient_flux%dwt_prod10_gain_grc                      (begg:endg                                   )) 
    allocate(nutrient_flux%dwt_prod100_gain_grc                     (begg:endg                                   )) 
    allocate(nutrient_flux%dwt_crop_product_gain_patch              (begp:endp                                   )) 
    allocate(nutrient_flux%product_loss_col                         (begc:endc                                   )) 
    allocate(nutrient_flux%deadcroot_storage_to_xfer_patch          (begp:endp                                   )) 
    allocate(nutrient_flux%deadstem_storage_to_xfer_patch           (begp:endp                                   )) 
    allocate(nutrient_flux%froot_storage_to_xfer_patch              (begp:endp                                   )) 
    allocate(nutrient_flux%grain_storage_to_xfer_patch              (begp:endp                                   )) 
    allocate(nutrient_flux%leaf_storage_to_xfer_patch               (begp:endp                                   )) 
    allocate(nutrient_flux%livecroot_storage_to_xfer_patch          (begp:endp                                   )) 
    allocate(nutrient_flux%livestem_storage_to_xfer_patch           (begp:endp                                   )) 
    allocate(nutrient_flux%decomp_pools_leached_col                 (begc:endc,1:ndecomp_pools                   )) 
    allocate(nutrient_flux%decomp_pools_sourcesink_col              (begc:endc,1:nlevdecomp_full,1:ndecomp_pools )) 
    allocate(nutrient_flux%decomp_pools_transport_tendency_col      (begc:endc,1:nlevdecomp_full,1:ndecomp_pools )) 
    allocate(nutrient_flux%leaf_to_litter_patch                     (begp:endp                                   )) 
    allocate(nutrient_flux%froot_to_litter_patch                    (begp:endp                                   )) 
    allocate(nutrient_flux%livestem_to_litter_patch                 (begp:endp                                   )) 
    allocate(nutrient_flux%grain_to_food_patch                      (begp:endp                                   )) 
    allocate(nutrient_flux%m_leaf_to_litter_patch                   (begp:endp                                   )) 
    allocate(nutrient_flux%m_leaf_storage_to_litter_patch           (begp:endp                                   )) 
    allocate(nutrient_flux%m_leaf_xfer_to_litter_patch              (begp:endp                                   )) 
    allocate(nutrient_flux%m_froot_to_litter_patch                  (begp:endp                                   )) 
    allocate(nutrient_flux%m_froot_storage_to_litter_patch          (begp:endp                                   )) 
    allocate(nutrient_flux%m_froot_xfer_to_litter_patch             (begp:endp                                   )) 
    allocate(nutrient_flux%m_livestem_to_litter_patch               (begp:endp                                   )) 
    allocate(nutrient_flux%m_livestem_storage_to_litter_patch       (begp:endp                                   )) 
    allocate(nutrient_flux%m_livestem_xfer_to_litter_patch          (begp:endp                                   )) 
    allocate(nutrient_flux%m_deadstem_storage_to_litter_patch       (begp:endp                                   )) 
    allocate(nutrient_flux%m_deadcroot_xfer_to_litter_patch         (begp:endp                                   )) 
    allocate(nutrient_flux%m_deadstem_to_litter_patch               (begp:endp                                   )) 
    allocate(nutrient_flux%m_livecroot_storage_to_litter_patch      (begp:endp                                   )) 
    allocate(nutrient_flux%m_livecroot_to_litter_patch              (begp:endp                                   )) 
    allocate(nutrient_flux%m_livecroot_xfer_to_litter_patch         (begp:endp                                   )) 
    allocate(nutrient_flux%m_deadcroot_to_litter_patch              (begp:endp                                   )) 
    allocate(nutrient_flux%m_deadcroot_storage_to_litter_patch      (begp:endp                                   )) 
    allocate(nutrient_flux%m_deadstem_xfer_to_litter_patch          (begp:endp                                   )) 
    allocate(nutrient_flux%m_pool_to_litter_patch                   (begp:endp                                   )) 
    allocate(nutrient_flux%gap_mortality_to_litr_met_col            (begc:endc,1:nlevdecomp_full                 )) 
    allocate(nutrient_flux%gap_mortality_to_litr_cel_col            (begc:endc,1:nlevdecomp_full                 )) 
    allocate(nutrient_flux%gap_mortality_to_litr_lig_col            (begc:endc,1:nlevdecomp_full                 )) 
    allocate(nutrient_flux%gap_mortality_to_cwd_col                 (begc:endc, 1:nlevdecomp_full                )) 
    allocate(nutrient_flux%hrv_leaf_storage_to_litter_patch         (begp:endp                                   )) 
    allocate(nutrient_flux%hrv_leaf_to_litter_patch                 (begp:endp                                   )) 
    allocate(nutrient_flux%hrv_leaf_xfer_to_litter_patch            (begp:endp                                   )) 
    allocate(nutrient_flux%hrv_froot_storage_to_litter_patch        (begp:endp                                   )) 
    allocate(nutrient_flux%hrv_froot_to_litter_patch                (begp:endp                                   )) 
    allocate(nutrient_flux%hrv_froot_xfer_to_litter_patch           (begp:endp                                   )) 
    allocate(nutrient_flux%hrv_livestem_storage_to_litter_patch     (begp:endp                                   )) 
    allocate(nutrient_flux%hrv_livestem_to_litter_patch             (begp:endp                                   )) 
    allocate(nutrient_flux%hrv_livestem_xfer_to_litter_patch        (begp:endp                                   )) 
    allocate(nutrient_flux%hrv_deadstem_to_prod10_patch             (begp:endp                                   )) 
    allocate(nutrient_flux%hrv_deadstem_to_prod100_patch            (begp:endp                                   )) 
    allocate(nutrient_flux%hrv_deadstem_storage_to_litter_patch     (begp:endp                                   )) 
    allocate(nutrient_flux%hrv_deadstem_xfer_to_litter_patch        (begp:endp                                   )) 
    allocate(nutrient_flux%hrv_deadcroot_storage_to_litter_patch    (begp:endp                                   )) 
    allocate(nutrient_flux%hrv_deadcroot_to_litter_patch            (begp:endp                                   )) 
    allocate(nutrient_flux%hrv_deadcroot_xfer_to_litter_patch       (begp:endp                                   )) 
    allocate(nutrient_flux%hrv_livecroot_storage_to_litter_patch    (begp:endp                                   )) 
    allocate(nutrient_flux%hrv_livecroot_to_litter_patch            (begp:endp                                   )) 
    allocate(nutrient_flux%hrv_livecroot_xfer_to_litter_patch       (begp:endp                                   )) 
    allocate(nutrient_flux%hrv_pool_to_litter_patch                 (begp:endp                                   )) 
    allocate(nutrient_flux%hrv_leaf_to_prod1_patch                  (begp:endp                                   )) 
    allocate(nutrient_flux%hrv_livestem_to_prod1_patch              (begp:endp                                   )) 
    allocate(nutrient_flux%hrv_grain_to_prod1_patch                 (begp:endp                                   )) 
    allocate(nutrient_flux%hrv_crop_to_prod1_patch                  (begp:endp                                   )) 
    allocate(nutrient_flux%harvest_to_cwd_col                       (begc:endc,1:nlevdecomp_full                 )) 
    allocate(nutrient_flux%harvest_to_litr_cel_col                  (begc:endc,1:nlevdecomp_full                 )) 
    allocate(nutrient_flux%harvest_to_litr_lig_col                  (begc:endc,1:nlevdecomp_full                 )) 
    allocate(nutrient_flux%harvest_to_litr_met_col                  (begc:endc,1:nlevdecomp_full                 )) 
    allocate(nutrient_flux%m_pool_to_fire_patch                     (begp:endp                                   )) 
    allocate(nutrient_flux%m_deadcroot_storage_to_fire_patch        (begp:endp                                   )) 
    allocate(nutrient_flux%m_deadcroot_to_fire_patch                (begp:endp                                   )) 
    allocate(nutrient_flux%m_deadcroot_xfer_to_fire_patch           (begp:endp                                   )) 
    allocate(nutrient_flux%m_deadstem_storage_to_fire_patch         (begp:endp                                   )) 
    allocate(nutrient_flux%m_deadstem_to_fire_patch                 (begp:endp                                   )) 
    allocate(nutrient_flux%m_deadstem_xfer_to_fire_patch            (begp:endp                                   )) 
    allocate(nutrient_flux%m_froot_storage_to_fire_patch            (begp:endp                                   )) 
    allocate(nutrient_flux%m_froot_to_fire_patch                    (begp:endp                                   )) 
    allocate(nutrient_flux%m_froot_xfer_to_fire_patch               (begp:endp                                   )) 
    allocate(nutrient_flux%m_leaf_storage_to_fire_patch             (begp:endp                                   )) 
    allocate(nutrient_flux%m_leaf_to_fire_patch                     (begp:endp                                   )) 
    allocate(nutrient_flux%m_leaf_xfer_to_fire_patch                (begp:endp                                   )) 
    allocate(nutrient_flux%m_livecroot_storage_to_fire_patch        (begp:endp                                   )) 
    allocate(nutrient_flux%m_livecroot_to_fire_patch                (begp:endp                                   )) 
    allocate(nutrient_flux%m_livecroot_xfer_to_fire_patch           (begp:endp                                   )) 
    allocate(nutrient_flux%m_livestem_storage_to_fire_patch         (begp:endp                                   )) 
    allocate(nutrient_flux%m_livestem_to_fire_patch                 (begp:endp                                   )) 
    allocate(nutrient_flux%m_livestem_xfer_to_fire_patch            (begp:endp                                   )) 
    allocate(nutrient_flux%m_livecroot_to_deadcroot_fire_patch      (begp:endp                                   )) 
    allocate(nutrient_flux%m_livestem_to_deadstem_fire_patch        (begp:endp                                   )) 
    allocate(nutrient_flux%fire_mortality_to_cwd_col                (begc:endc, 1:nlevdecomp_full                )) 
    allocate(nutrient_flux%m_pool_to_litter_fire_patch              (begp:endp                                   )) 
    allocate(nutrient_flux%m_deadcroot_storage_to_litter_fire_patch (begp:endp                                   )) 
    allocate(nutrient_flux%m_deadcroot_to_litter_fire_patch         (begp:endp                                   )) 
    allocate(nutrient_flux%m_deadcroot_xfer_to_litter_fire_patch    (begp:endp                                   )) 
    allocate(nutrient_flux%m_deadstem_storage_to_litter_fire_patch  (begp:endp                                   )) 
    allocate(nutrient_flux%m_deadstem_to_litter_fire_patch          (begp:endp                                   )) 
    allocate(nutrient_flux%m_deadstem_xfer_to_litter_fire_patch     (begp:endp                                   )) 
    allocate(nutrient_flux%m_froot_storage_to_litter_fire_patch     (begp:endp                                   )) 
    allocate(nutrient_flux%m_froot_to_litter_fire_patch             (begp:endp                                   )) 
    allocate(nutrient_flux%m_froot_xfer_to_litter_fire_patch        (begp:endp                                   )) 
    allocate(nutrient_flux%m_leaf_storage_to_litter_fire_patch      (begp:endp                                   )) 
    allocate(nutrient_flux%m_leaf_to_litter_fire_patch              (begp:endp                                   )) 
    allocate(nutrient_flux%m_leaf_xfer_to_litter_fire_patch         (begp:endp                                   )) 
    allocate(nutrient_flux%m_livecroot_storage_to_litter_fire_patch (begp:endp                                   )) 
    allocate(nutrient_flux%m_livecroot_to_litter_fire_patch         (begp:endp                                   )) 
    allocate(nutrient_flux%m_livecroot_xfer_to_litter_fire_patch    (begp:endp                                   )) 
    allocate(nutrient_flux%m_livestem_storage_to_litter_fire_patch  (begp:endp                                   )) 
    allocate(nutrient_flux%m_livestem_to_litter_fire_patch          (begp:endp                                   )) 
    allocate(nutrient_flux%m_livestem_xfer_to_litter_fire_patch     (begp:endp                                   )) 
    allocate(nutrient_flux%m_to_litr_cel_fire_col                   (begc:endc,1:nlevdecomp_full                 )) 
    allocate(nutrient_flux%m_to_litr_lig_fire_col                   (begc:endc,1:nlevdecomp_full                 )) 
    allocate(nutrient_flux%m_to_litr_met_fire_col                   (begc:endc,1:nlevdecomp_full                 )) 
    allocate(nutrient_flux%fire_loss_col                            (begc:endc                                   )) 
    allocate(nutrient_flux%fire_loss_patch                          (begp:endp                                   )) 
    allocate(nutrient_flux%prod100_loss_col                         (begc:endc                                   )) 
    allocate(nutrient_flux%prod10_loss_col                          (begc:endc                                   )) 
    allocate(nutrient_flux%prod1_loss_col                           (begc:endc                                   )) 
    allocate(nutrient_flux%inputs_patch                             (begp:endp                                   )) 
    allocate(nutrient_flux%outputs_patch                            (begp:endp                                   )) 
    allocate(nutrient_flux%plant_alloc_patch                        (begp:endp                                   )) 
    allocate(nutrient_flux%som_leached_col                          (begc:endc                                   )) 
    allocate(nutrient_flux%wood_harvest_col                         (begc:endc                                   )) 
    allocate(nutrient_flux%wood_harvest_patch                       (begp:endp                                   )) 
    allocate(nutrient_flux%pool_to_deadcroot_patch                  (begp:endp                                   )) 
    allocate(nutrient_flux%pool_to_deadcroot_storage_patch          (begp:endp                                   )) 
    allocate(nutrient_flux%pool_to_deadstem_patch                   (begp:endp                                   )) 
    allocate(nutrient_flux%pool_to_deadstem_storage_patch           (begp:endp                                   )) 
    allocate(nutrient_flux%pool_to_froot_patch                      (begp:endp                                   )) 
    allocate(nutrient_flux%pool_to_froot_storage_patch              (begp:endp                                   )) 
    allocate(nutrient_flux%pool_to_grain_patch                      (begp:endp                                   )) 
    allocate(nutrient_flux%pool_to_grain_storage_patch              (begp:endp                                   )) 
    allocate(nutrient_flux%pool_to_leaf_patch                       (begp:endp                                   )) 
    allocate(nutrient_flux%pool_to_leaf_storage_patch               (begp:endp                                   )) 
    allocate(nutrient_flux%pool_to_livecroot_patch                  (begp:endp                                   )) 
    allocate(nutrient_flux%pool_to_livecroot_storage_patch          (begp:endp                                   )) 
    allocate(nutrient_flux%pool_to_livestem_patch                   (begp:endp                                   )) 
    allocate(nutrient_flux%pool_to_livestem_storage_patch           (begp:endp                                   )) 
    allocate(nutrient_flux%phenology_to_litr_met_col                (begc:endc,1:nlevdecomp_full                 )) 
    allocate(nutrient_flux%phenology_to_litr_cel_col                (begc:endc,1:nlevdecomp_full                 )) 
    allocate(nutrient_flux%phenology_to_litr_lig_col                (begc:endc,1:nlevdecomp_full                 )) 
    allocate(nutrient_flux%livecroot_to_deadcroot_patch             (begp:endp                                   )) 
    allocate(nutrient_flux%deadcroot_xfer_to_deadcroot_patch        (begp:endp                                   )) 
    allocate(nutrient_flux%deadstem_xfer_to_deadstem_patch          (begp:endp                                   )) 
    allocate(nutrient_flux%livecroot_xfer_to_livecroot_patch        (begp:endp                                   )) 
    allocate(nutrient_flux%froot_xfer_to_froot_patch                (begp:endp                                   )) 
    allocate(nutrient_flux%grain_xfer_to_grain_patch                (begp:endp                                   )) 
    allocate(nutrient_flux%leaf_xfer_to_leaf_patch                  (begp:endp                                   )) 
    allocate(nutrient_flux%livestem_to_deadstem_patch               (begp:endp                                   )) 
    allocate(nutrient_flux%livestem_xfer_to_livestem_patch          (begp:endp                                   )) 
    allocate(nutrient_flux%prev_leaf_to_litter_patch                (begp:endp                                   )) 
    allocate(nutrient_flux%prev_froot_to_litter_patch               (begp:endp                                   )) 

    nutrient_flux%dwt_seed_to_deadstem_patch               (:)     = nan
    nutrient_flux%dwt_seed_to_leaf_patch                   (:)     = nan
    nutrient_flux%dwt_seed_to_deadstem_grc                 (:)     = nan
    nutrient_flux%dwt_seed_to_leaf_grc                     (:)     = nan
    nutrient_flux%dwt_slash_flux_col                       (:)     = nan
    nutrient_flux%dwt_loss_col                             (:)     = nan
    nutrient_flux%dwt_deadcroot_to_cwd_col                 (:,:)   = nan
    nutrient_flux%dwt_froot_to_litr_cel_col                (:,:)   = nan
    nutrient_flux%dwt_froot_to_litr_lig_col                (:,:)   = nan
    nutrient_flux%dwt_froot_to_litr_met_col                (:,:)   = nan
    nutrient_flux%dwt_livecroot_to_cwd_col                 (:,:)   = nan
    nutrient_flux%dwt_conv_flux_patch                      (:)     = nan
    nutrient_flux%dwt_conv_flux_col                        (:)     = nan
    nutrient_flux%dwt_conv_flux_grc                        (:)     = nan
    nutrient_flux%dwt_prod10_gain_patch                    (:)     = nan
    nutrient_flux%dwt_prod100_gain_patch                   (:)     = nan
    nutrient_flux%dwt_prod10_gain_col                      (:)     = nan
    nutrient_flux%dwt_prod100_gain_col                     (:)     = nan
    nutrient_flux%dwt_prod10_gain_grc                      (:)     = nan
    nutrient_flux%dwt_prod100_gain_grc                     (:)     = nan
    nutrient_flux%dwt_crop_product_gain_patch              (:)     = nan
    nutrient_flux%product_loss_col                         (:)     = nan
    nutrient_flux%deadcroot_storage_to_xfer_patch          (:)     = nan
    nutrient_flux%deadstem_storage_to_xfer_patch           (:)     = nan
    nutrient_flux%froot_storage_to_xfer_patch              (:)     = nan
    nutrient_flux%grain_storage_to_xfer_patch              (:)     = nan
    nutrient_flux%leaf_storage_to_xfer_patch               (:)     = nan
    nutrient_flux%livecroot_storage_to_xfer_patch          (:)     = nan
    nutrient_flux%livestem_storage_to_xfer_patch           (:)     = nan
    nutrient_flux%decomp_pools_leached_col                 (:,:)   = nan
    nutrient_flux%decomp_pools_sourcesink_col              (:,:,:) = nan
    nutrient_flux%decomp_pools_transport_tendency_col      (:,:,:) = nan
    nutrient_flux%leaf_to_litter_patch                     (:)     = nan
    nutrient_flux%froot_to_litter_patch                    (:)     = nan
    nutrient_flux%livestem_to_litter_patch                 (:)     = nan
    nutrient_flux%grain_to_food_patch                      (:)     = nan
    nutrient_flux%m_leaf_to_litter_patch                   (:)     = nan
    nutrient_flux%m_leaf_storage_to_litter_patch           (:)     = nan
    nutrient_flux%m_leaf_xfer_to_litter_patch              (:)     = nan
    nutrient_flux%m_froot_to_litter_patch                  (:)     = nan
    nutrient_flux%m_froot_storage_to_litter_patch          (:)     = nan
    nutrient_flux%m_froot_xfer_to_litter_patch             (:)     = nan
    nutrient_flux%m_livestem_to_litter_patch               (:)     = nan
    nutrient_flux%m_livestem_storage_to_litter_patch       (:)     = nan
    nutrient_flux%m_livestem_xfer_to_litter_patch          (:)     = nan
    nutrient_flux%m_deadstem_storage_to_litter_patch       (:)     = nan
    nutrient_flux%m_deadcroot_xfer_to_litter_patch         (:)     = nan
    nutrient_flux%m_deadstem_to_litter_patch               (:)     = nan
    nutrient_flux%m_livecroot_storage_to_litter_patch      (:)     = nan
    nutrient_flux%m_livecroot_to_litter_patch              (:)     = nan
    nutrient_flux%m_livecroot_xfer_to_litter_patch         (:)     = nan
    nutrient_flux%m_deadcroot_to_litter_patch              (:)     = nan
    nutrient_flux%m_deadcroot_storage_to_litter_patch      (:)     = nan
    nutrient_flux%m_deadstem_xfer_to_litter_patch          (:)     = nan
    nutrient_flux%m_pool_to_litter_patch                   (:)     = nan
    nutrient_flux%gap_mortality_to_litr_met_col            (:,:)   = nan
    nutrient_flux%gap_mortality_to_litr_cel_col            (:,:)   = nan
    nutrient_flux%gap_mortality_to_litr_lig_col            (:,:)   = nan
    nutrient_flux%gap_mortality_to_cwd_col                 (:,:)   = nan
    nutrient_flux%hrv_leaf_storage_to_litter_patch         (:)     = nan
    nutrient_flux%hrv_leaf_to_litter_patch                 (:)     = nan
    nutrient_flux%hrv_leaf_xfer_to_litter_patch            (:)     = nan
    nutrient_flux%hrv_froot_storage_to_litter_patch        (:)     = nan
    nutrient_flux%hrv_froot_to_litter_patch                (:)     = nan
    nutrient_flux%hrv_froot_xfer_to_litter_patch           (:)     = nan
    nutrient_flux%hrv_livestem_storage_to_litter_patch     (:)     = nan
    nutrient_flux%hrv_livestem_to_litter_patch             (:)     = nan
    nutrient_flux%hrv_livestem_xfer_to_litter_patch        (:)     = nan
    nutrient_flux%hrv_deadstem_to_prod10_patch             (:)     = nan
    nutrient_flux%hrv_deadstem_to_prod100_patch            (:)     = nan
    nutrient_flux%hrv_deadstem_storage_to_litter_patch     (:)     = nan
    nutrient_flux%hrv_deadstem_xfer_to_litter_patch        (:)     = nan
    nutrient_flux%hrv_deadcroot_storage_to_litter_patch    (:)     = nan
    nutrient_flux%hrv_deadcroot_to_litter_patch            (:)     = nan
    nutrient_flux%hrv_deadcroot_xfer_to_litter_patch       (:)     = nan
    nutrient_flux%hrv_livecroot_storage_to_litter_patch    (:)     = nan
    nutrient_flux%hrv_livecroot_to_litter_patch            (:)     = nan
    nutrient_flux%hrv_livecroot_xfer_to_litter_patch       (:)     = nan
    nutrient_flux%hrv_pool_to_litter_patch                 (:)     = nan
    nutrient_flux%hrv_leaf_to_prod1_patch                  (:)     = nan
    nutrient_flux%hrv_livestem_to_prod1_patch              (:)     = nan
    nutrient_flux%hrv_grain_to_prod1_patch                 (:)     = nan
    nutrient_flux%hrv_crop_to_prod1_patch                  (:)     = nan
    nutrient_flux%harvest_to_cwd_col                       (:,:)   = nan
    nutrient_flux%harvest_to_litr_cel_col                  (:,:)   = nan
    nutrient_flux%harvest_to_litr_lig_col                  (:,:)   = nan
    nutrient_flux%harvest_to_litr_met_col                  (:,:)   = nan
    nutrient_flux%m_pool_to_fire_patch                     (:)     = nan
    nutrient_flux%m_deadcroot_storage_to_fire_patch        (:)     = nan
    nutrient_flux%m_deadcroot_to_fire_patch                (:)     = nan
    nutrient_flux%m_deadcroot_xfer_to_fire_patch           (:)     = nan
    nutrient_flux%m_deadstem_storage_to_fire_patch         (:)     = nan
    nutrient_flux%m_deadstem_to_fire_patch                 (:)     = nan
    nutrient_flux%m_deadstem_xfer_to_fire_patch            (:)     = nan
    nutrient_flux%m_froot_storage_to_fire_patch            (:)     = nan
    nutrient_flux%m_froot_to_fire_patch                    (:)     = nan
    nutrient_flux%m_froot_xfer_to_fire_patch               (:)     = nan
    nutrient_flux%m_leaf_storage_to_fire_patch             (:)     = nan
    nutrient_flux%m_leaf_to_fire_patch                     (:)     = nan
    nutrient_flux%m_leaf_xfer_to_fire_patch                (:)     = nan
    nutrient_flux%m_livecroot_storage_to_fire_patch        (:)     = nan
    nutrient_flux%m_livecroot_to_fire_patch                (:)     = nan
    nutrient_flux%m_livecroot_xfer_to_fire_patch           (:)     = nan
    nutrient_flux%m_livestem_storage_to_fire_patch         (:)     = nan
    nutrient_flux%m_livestem_to_fire_patch                 (:)     = nan
    nutrient_flux%m_livestem_xfer_to_fire_patch            (:)     = nan
    nutrient_flux%m_livecroot_to_deadcroot_fire_patch      (:)     = nan
    nutrient_flux%m_livestem_to_deadstem_fire_patch        (:)     = nan
    nutrient_flux%fire_mortality_to_cwd_col                (:,:)   = nan
    nutrient_flux%m_pool_to_litter_fire_patch              (:)     = nan
    nutrient_flux%m_deadcroot_storage_to_litter_fire_patch (:)     = nan
    nutrient_flux%m_deadcroot_to_litter_fire_patch         (:)     = nan
    nutrient_flux%m_deadcroot_xfer_to_litter_fire_patch    (:)     = nan
    nutrient_flux%m_deadstem_storage_to_litter_fire_patch  (:)     = nan
    nutrient_flux%m_deadstem_to_litter_fire_patch          (:)     = nan
    nutrient_flux%m_deadstem_xfer_to_litter_fire_patch     (:)     = nan
    nutrient_flux%m_froot_storage_to_litter_fire_patch     (:)     = nan
    nutrient_flux%m_froot_to_litter_fire_patch             (:)     = nan
    nutrient_flux%m_froot_xfer_to_litter_fire_patch        (:)     = nan
    nutrient_flux%m_leaf_storage_to_litter_fire_patch      (:)     = nan
    nutrient_flux%m_leaf_to_litter_fire_patch              (:)     = nan
    nutrient_flux%m_leaf_xfer_to_litter_fire_patch         (:)     = nan
    nutrient_flux%m_livecroot_storage_to_litter_fire_patch (:)     = nan
    nutrient_flux%m_livecroot_to_litter_fire_patch         (:)     = nan
    nutrient_flux%m_livecroot_xfer_to_litter_fire_patch    (:)     = nan
    nutrient_flux%m_livestem_storage_to_litter_fire_patch  (:)     = nan
    nutrient_flux%m_livestem_to_litter_fire_patch          (:)     = nan
    nutrient_flux%m_livestem_xfer_to_litter_fire_patch     (:)     = nan
    nutrient_flux%m_to_litr_cel_fire_col                   (:,:)   = nan
    nutrient_flux%m_to_litr_lig_fire_col                   (:,:)   = nan
    nutrient_flux%m_to_litr_met_fire_col                   (:,:)   = nan
    nutrient_flux%fire_loss_col                            (:)     = nan
    nutrient_flux%fire_loss_patch                          (:)     = nan
    nutrient_flux%prod100_loss_col                         (:)     = nan
    nutrient_flux%prod10_loss_col                          (:)     = nan
    nutrient_flux%prod1_loss_col                           (:)     = nan
    nutrient_flux%inputs_patch                             (:)     = nan
    nutrient_flux%outputs_patch                            (:)     = nan
    nutrient_flux%plant_alloc_patch                        (:)     = nan
    nutrient_flux%som_leached_col                          (:)     = nan
    nutrient_flux%wood_harvest_col                         (:)     = nan
    nutrient_flux%wood_harvest_patch                       (:)     = nan
    nutrient_flux%pool_to_deadcroot_patch                  (:)     = nan
    nutrient_flux%pool_to_deadcroot_storage_patch          (:)     = nan
    nutrient_flux%pool_to_deadstem_patch                   (:)     = nan
    nutrient_flux%pool_to_deadstem_storage_patch           (:)     = nan
    nutrient_flux%pool_to_froot_patch                      (:)     = nan
    nutrient_flux%pool_to_froot_storage_patch              (:)     = nan
    nutrient_flux%pool_to_grain_patch                      (:)     = nan
    nutrient_flux%pool_to_grain_storage_patch              (:)     = nan
    nutrient_flux%pool_to_leaf_patch                       (:)     = nan
    nutrient_flux%pool_to_leaf_storage_patch               (:)     = nan
    nutrient_flux%pool_to_livecroot_patch                  (:)     = nan
    nutrient_flux%pool_to_livecroot_storage_patch          (:)     = nan
    nutrient_flux%pool_to_livestem_patch                   (:)     = nan
    nutrient_flux%pool_to_livestem_storage_patch           (:)     = nan
    nutrient_flux%phenology_to_litr_met_col                (:,:)   = nan
    nutrient_flux%phenology_to_litr_cel_col                (:,:)   = nan
    nutrient_flux%phenology_to_litr_lig_col                (:,:)   = nan
    nutrient_flux%livecroot_to_deadcroot_patch             (:)     = nan
    nutrient_flux%deadcroot_xfer_to_deadcroot_patch        (:)     = nan
    nutrient_flux%deadstem_xfer_to_deadstem_patch          (:)     = nan
    nutrient_flux%livecroot_xfer_to_livecroot_patch        (:)     = nan
    nutrient_flux%froot_xfer_to_froot_patch                (:)     = nan
    nutrient_flux%grain_xfer_to_grain_patch                (:)     = nan
    nutrient_flux%leaf_xfer_to_leaf_patch                  (:)     = nan
    nutrient_flux%livestem_to_deadstem_patch               (:)     = nan
    nutrient_flux%livestem_xfer_to_livestem_patch          (:)     = nan
    nutrient_flux%prev_leaf_to_litter_patch                (:)     = nan
    nutrient_flux%prev_froot_to_litter_patch               (:)     = nan

  end subroutine NutrientFluxInitAllocate

end module NutrientFluxType
