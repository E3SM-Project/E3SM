module NutrientFluxType

#include "shr_assert.h"

  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_infnan_mod         , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use decompMod              , only : bounds_type
  use clm_varcon             , only : spval, ispval
  use abortutils             , only : endrun
  use clm_varpar             , only : ndecomp_pools, nlevdecomp_full
  use clm_varpar             , only : crop_prog
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private

  type, public :: nutrientflux_type

     character(len=3)  :: name
     character(len=4)  :: history_name_prefix
     character(len=3)  :: restart_name

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

  public :: NutrientFluxInitAllocate, &
       NutrientFluxInitHistory

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

  !------------------------------------------------------------------------
  subroutine NutrientFluxInitHistory(nf, bounds)
    !
    use histFileMod, only : hist_addfld1d, hist_addfld2d, hist_addfld_decomp
    !
    implicit none
    !
    ! !ARGUMENTS:
    class (nutrientflux_type)    :: nf
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer           :: l
    integer           :: begp,endp
    integer           :: begc,endc
    integer           :: begg,endg
    character(len=10) :: unit
    character(len=3)  :: name
    character(len=4)  :: name_prefix
    character(24)     :: fieldname
    character(100)    :: longname
    real(r8), pointer :: data1dptr(:)   ! temp. pointer for slicing larger arrays
    real(r8), pointer :: data2dptr(:,:) ! temp. pointer for slicing larger arrays
    !------------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc
    begg = bounds%begg; endg = bounds%endg

    write(unit,*) 'g' // trim(nf%name) // '/m^2'
    unit = trim(unit)

    name = nf%name
    name_prefix = nf%history_name_prefix

    if (crop_prog) then
       nf%grain_to_food_patch(begp:endp) = spval
       call hist_addfld1d (fname=trim(name_prefix)//'GRAIN'//trim(name) //'_TO_FOOD', &
            units='g'//trim(nf%name)//'/m^2/s', &
            avgflag='A', long_name='grain ' // trim(nf%name) // ' to food', &
            ptr_patch=nf%grain_to_food_patch, default='inactive')
    end if

    nf%dwt_seed_to_deadstem_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'DWT_SEED'//trim(nf%name)//'_TO_DEADSTEM_PATCH', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', &
         long_name='patch-level seed source to patch-level deadstem ' // &
         '(per-area-gridcell; only makes sense with dov2xy=.false.)', &
         ptr_patch=nf%dwt_seed_to_deadstem_patch, default='inactive')

    nf%dwt_seed_to_leaf_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'DWT_SEED'//trim(nf%name)//'_TO_LEAF_PATCH', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', &
         long_name='patch-level seed source to patch-level leaf ' // &
         '(per-area-gridcell; only makes sense with dov2xy=.false.)', &
         ptr_patch=nf%dwt_seed_to_leaf_patch, default='inactive')

    nf%dwt_seed_to_deadstem_grc(begg:endg) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'DWT_SEED'//trim(nf%name)//'_TO_DEADSTEM_GRC', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='seed source to patch-level deadstem', &
         ptr_gcell=nf%dwt_seed_to_deadstem_grc, default='inactive')

    nf%dwt_seed_to_leaf_grc(begg:endg) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'DWT_SEED'//trim(nf%name)//'_TO_LEAF_GRC', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='seed source to patch-level leaf', &
         ptr_gcell=nf%dwt_seed_to_leaf_grc, default='inactive')

    nf%dwt_slash_flux_col(begc:endc) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'DWT_SLASH_'//trim(nf%name)//'FLUX', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='slash ' // trim(nf%name) // ' flux to litter and CWD due to land use', &
         ptr_col=nf%dwt_slash_flux_col)

    nf%dwt_loss_col(begc:endc) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'DWT_'//trim(nf%name)//'LOSS', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='total ' // trim(nf%name) // ' loss from land cover conversion', &
         ptr_col=nf%dwt_loss_col, default='inactive')

    nf%dwt_deadcroot_to_cwd_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname=trim(name_prefix)//'DWT_DEADCROOT'//trim(nf%name)//'_TO_CWD'//trim(nf%name), &
         units='g'//trim(nf%name)//'/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='dead coarse root to CWD due to landcover change', &
         ptr_col=nf%dwt_deadcroot_to_cwd_col, default='inactive')

    nf%dwt_froot_to_litr_cel_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname=trim(name_prefix)//'DWT_FROOT'//trim(nf%name)//'_TO_LITR_CEL_'//trim(nf%name), &
         units='g'//trim(nf%name)//'/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='fine root to litter due to landcover change', &
         ptr_col=nf%dwt_froot_to_litr_cel_col, default='inactive')

    nf%dwt_froot_to_litr_lig_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname=trim(name_prefix)//'DWT_FROOT'//trim(nf%name)//'_TO_LITR_LIG_'//trim(nf%name), &
         units='g'//trim(nf%name)//'/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='fine root to litter due to landcover change', &
         ptr_col=nf%dwt_froot_to_litr_lig_col, default='inactive')

    nf%dwt_froot_to_litr_met_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname=trim(name_prefix)//'DWT_FROOT'//trim(nf%name)//'_TO_LITR_MET_'//trim(nf%name), &
         units='g'//trim(nf%name)//'/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='fine root to litter due to landcover change', &
         ptr_col=nf%dwt_froot_to_litr_met_col, default='inactive')

    nf%dwt_livecroot_to_cwd_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname=trim(name_prefix)//'DWT_LIVECROOT'//trim(nf%name)//'_TO_CWD'//trim(nf%name), &
         units='g'//trim(nf%name)//'/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='live coarse root to CWD due to landcover change', &
         ptr_col=nf%dwt_livecroot_to_cwd_col, default='inactive')


    nf%dwt_conv_flux_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'DWT_CONV_'//trim(nf%name)//'FLUX_PATCH', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', &
         long_name='patch-level conversion ' // trim(nf%name) // ' flux (immediate loss to atm) ' // &
         '(0 at all times except first timestep of year) ' // &
         '(per-area-gridcell; only makes sense with dov2xy=.false.)', &
         ptr_patch=nf%dwt_conv_flux_patch, default='inactive')

    nf%dwt_conv_flux_col(begc:endc) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'DWT_CONV_'//trim(nf%name)//'FLUX', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='conversion ' // trim(nf%name) // ' flux (immediate loss to atm)', &
         ptr_col=nf%dwt_conv_flux_col, default='inactive')

    nf%dwt_conv_flux_grc(begg:endg) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'DWT_CONV_'//trim(nf%name)//'FLUX_GRC', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', &
         long_name='conversion ' // trim(nf%name) // ' flux (immediate loss to atm) (0 at all times except first timestep of year)', &
         ptr_gcell=nf%dwt_conv_flux_grc)

    nf%dwt_prod10_gain_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'DWT_PROD10' // trim(nf%name)//'_GAIN_PATCH', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='landcover change-driven addition to 10-yr wood product pool', &
         ptr_col=nf%dwt_prod10_gain_patch, default='inactive')
    
    nf%dwt_prod100_gain_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'DWT_PROD100' // trim(nf%name)//'_GAIN_PATCH', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='landcover change-driven addition to 100-yr wood product pool', &
         ptr_col=nf%dwt_prod100_gain_patch, default='inactive')

    nf%dwt_prod10_gain_col(begc:endc) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'DWT_PROD10' // trim(nf%name)//'_GAIN', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='landcover change-driven addition to 10-yr wood product pool', &
         ptr_col=nf%dwt_prod10_gain_col, default='inactive')

    nf%dwt_prod100_gain_col(begc:endc) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'DWT_PROD100' // trim(nf%name)//'_GAIN', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='landcover change-driven addition to 100-yr wood product pool', &
         ptr_col=nf%dwt_prod100_gain_col, default='inactive')

    nf%dwt_prod10_gain_grc(begg:endg) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'DWT_PROD10' // trim(nf%name)//'_GAIN_GRC', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='landcover change-driven addition to 10-yr wood product pool', &
         ptr_col=nf%dwt_prod10_gain_grc, default='inactive')

    nf%dwt_prod100_gain_grc(begg:endg) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'DWT_PROD100' // trim(nf%name)//'_GAIN_GRC', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='landcover change-driven addition to 100-yr wood product pool', &
         ptr_col=nf%dwt_prod100_gain_grc, default='inactive')

    nf%product_loss_col(begc:endc) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'PRODUCT_'//trim(nf%name)//'LOSS', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='total ' // trim(nf%name) // ' loss from wood product pools', &
         ptr_col=nf%product_loss_col, default='inactive')

    nf%deadcroot_storage_to_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'DEADCROOT'//trim(nf%name)//'_STORAGE_TO_XFER', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='dead coarse root ' // trim(nf%name) // ' shift storage to transfer', &
         ptr_patch=nf%deadcroot_storage_to_xfer_patch, default='inactive')
    nf%deadstem_storage_to_xfer_patch(begp:endp) = spval

    call hist_addfld1d (fname=trim(name_prefix)//'DEADSTEM' // trim(nf%name)//'_STORAGE_TO_XFER', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='dead stem ' // trim(nf%name) // ' shift storage to transfer', &
         ptr_patch=nf%deadstem_storage_to_xfer_patch, default='inactive')

    nf%froot_storage_to_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'FROOT'//trim(nf%name)//'_STORAGE_TO_XFER', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='fine root ' // trim(nf%name) // ' shift storage to transfer', &
         ptr_patch=nf%froot_storage_to_xfer_patch, default='inactive')

    nf%leaf_storage_to_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'LEAF' // trim(nf%name)//'_STORAGE_TO_XFER', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='leaf ' // trim(nf%name) // ' shift storage to transfer', &
         ptr_patch=nf%leaf_storage_to_xfer_patch, default='inactive')

    nf%livecroot_storage_to_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'LIVECROOT'//trim(nf%name)//'_STORAGE_TO_XFER', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='live coarse root ' // trim(nf%name) // ' shift storage to transfer', &
         ptr_patch=nf%livecroot_storage_to_xfer_patch, default='inactive')

    nf%livestem_storage_to_xfer_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'LIVESTEM' // trim(nf%name)//'_STORAGE_TO_XFER', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='live stem ' // trim(nf%name) // ' shift storage to transfer', &
         ptr_patch=nf%livestem_storage_to_xfer_patch, default='inactive')

    nf%froot_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'FROOT' // trim(nf%name)//'_TO_LITTER', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='fine root C litterfall', &
         ptr_patch=nf%froot_to_litter_patch, default='inactive')

    nf%m_leaf_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_LEAF' // trim(nf%name)//'_TO_LITTER', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='leaf ' // trim(nf%name) // ' mortality', &
         ptr_patch=nf%m_leaf_to_litter_patch, default='inactive')

    nf%m_froot_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_FROOT'//trim(nf%name)//'_TO_LITTER', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='fine root ' // trim(nf%name) // ' mortality', &
         ptr_patch=nf%m_froot_to_litter_patch, default='inactive')

    nf%m_livestem_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_LIVESTEM' // trim(nf%name)//'_TO_LITTER', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='live stem ' // trim(nf%name) // ' mortality', &
         ptr_patch=nf%m_livestem_to_litter_patch, default='inactive')

    nf%m_leaf_storage_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_LEAF' // trim(nf%name)//'_STORAGE_TO_LITTER', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='leaf ' // trim(nf%name) // ' storage mortality', &
         ptr_patch=nf%m_leaf_storage_to_litter_patch, default='inactive')

    nf%m_leaf_xfer_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_LEAF' // trim(nf%name)//'_XFER_TO_LITTER', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='leaf ' // trim(nf%name) // ' transfer mortality', &
         ptr_patch=nf%m_leaf_xfer_to_litter_patch, default='inactive')

    nf%m_froot_storage_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_FROOT'//trim(nf%name)//'_STORAGE_TO_LITTER', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='fine root ' // trim(nf%name) // ' storage mortality', &
         ptr_patch=nf%m_froot_storage_to_litter_patch, default='inactive')

    nf%m_froot_xfer_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_FROOT'//trim(nf%name)//'_XFER_TO_LITTER', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='fine root ' // trim(nf%name) // ' transfer mortality', &
         ptr_patch=nf%m_froot_xfer_to_litter_patch, default='inactive')

    nf%m_livestem_storage_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_LIVESTEM' // trim(nf%name)//'_STORAGE_TO_LITTER', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='live stem ' // trim(nf%name) // ' storage mortality', &
         ptr_patch=nf%m_livestem_storage_to_litter_patch, default='inactive')

    nf%m_livestem_xfer_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_LIVESTEM' // trim(nf%name)//'_XFER_TO_LITTER', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='live stem ' // trim(nf%name) // ' transfer mortality', &
         ptr_patch=nf%m_livestem_xfer_to_litter_patch, default='inactive')

    nf%m_deadstem_storage_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_DEADSTEM' // trim(nf%name)//'_STORAGE_TO_LITTER', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='dead stem ' // trim(nf%name) // ' storage mortality', &
         ptr_patch=nf%m_deadstem_storage_to_litter_patch, default='inactive')

    nf%m_deadcroot_xfer_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_DEADCROOT'//trim(nf%name)//'_XFER_TO_LITTER', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='dead coarse root ' // trim(nf%name) // ' transfer mortality', &
         ptr_patch=nf%m_deadcroot_xfer_to_litter_patch, default='inactive')
    
    nf%m_deadstem_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_DEADSTEM' // trim(nf%name)//'_TO_LITTER', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='dead stem ' // trim(nf%name) // ' mortality', &
         ptr_patch=nf%m_deadstem_to_litter_patch, default='inactive')

    nf%m_livecroot_storage_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_LIVECROOT'//trim(nf%name)//'_STORAGE_TO_LITTER', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='live coarse root ' // trim(nf%name) // ' storage mortality', &
         ptr_patch=nf%m_livecroot_storage_to_litter_patch, default='inactive')

    nf%m_livecroot_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_LIVECROOT'//trim(nf%name)//'_TO_LITTER', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='live coarse root ' // trim(nf%name) // ' mortality', &
         ptr_patch=nf%m_livecroot_to_litter_patch, default='inactive')

    nf%m_livecroot_xfer_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_LIVECROOT'//trim(nf%name)//'_XFER_TO_LITTER', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='live coarse root ' // trim(nf%name) // ' transfer mortality', &
         ptr_patch=nf%m_livecroot_xfer_to_litter_patch, default='inactive')

    nf%m_deadcroot_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_DEADCROOT'//trim(nf%name)//'_TO_LITTER', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='dead coarse root ' // trim(nf%name) // ' mortality', &
         ptr_patch=nf%m_deadcroot_to_litter_patch, default='inactive')
    
    nf%m_deadcroot_storage_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_DEADCROOT'//trim(nf%name)//'_STORAGE_TO_LITTER', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='dead coarse root ' // trim(nf%name) // ' storage mortality', &
         ptr_patch=nf%m_deadcroot_storage_to_litter_patch, default='inactive')

    nf%m_deadstem_xfer_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_DEADSTEM' // trim(nf%name)//'_XFER_TO_LITTER', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='dead stem ' // trim(nf%name) // ' transfer mortality', &
         ptr_patch=nf%m_deadstem_xfer_to_litter_patch, default='inactive')

    nf%m_pool_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_'//trim(nf%name)//'POOL_TO_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='cpool fire loss', &
         ptr_patch=nf%m_pool_to_fire_patch, default='inactive')

    nf%m_deadcroot_storage_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_DEADCROOT'//trim(name)//'_STORAGE_TO_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='dead root ' // trim(nf%name) // ' storage fire loss', &
         ptr_patch=nf%m_deadcroot_storage_to_fire_patch, default='inactive')

    nf%m_deadcroot_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_DEADCROOT'//trim(name)//'_TO_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='dead root ' // trim(nf%name) // ' fire loss', &
         ptr_patch=nf%m_deadcroot_to_fire_patch, default='inactive')

    nf%m_deadcroot_xfer_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_DEADCROOT'//trim(name)//'_XFER_TO_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='dead root ' // trim(nf%name) // ' transfer fire loss', &
         ptr_patch=nf%m_deadcroot_xfer_to_fire_patch, default='inactive')

    nf%m_deadstem_storage_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_DEADSTEM' // trim(nf%name)//'_STORAGE_TO_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='dead stem ' // trim(nf%name) // ' storage fire loss', &
         ptr_patch=nf%m_deadstem_storage_to_fire_patch, default='inactive')

    nf%m_deadstem_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_DEADSTEM' // trim(nf%name)//'_TO_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='dead stem ' // trim(nf%name) // ' fire loss', &
         ptr_patch=nf%m_deadstem_to_fire_patch, default='inactive')

    nf%m_deadstem_xfer_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_DEADSTEM' // trim(nf%name)//'_XFER_TO_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='dead stem ' // trim(nf%name) // ' transfer fire loss', &
         ptr_patch=nf%m_deadstem_xfer_to_fire_patch, default='inactive')

    nf%m_froot_storage_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_FROOT'//trim(nf%name)//'_STORAGE_TO_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='fine root ' // trim(nf%name) // ' storage fire loss', &
         ptr_patch=nf%m_froot_storage_to_fire_patch, default='inactive')

    nf%m_froot_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_FROOT'//trim(nf%name)//'_TO_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='fine root ' // trim(nf%name) // ' fire loss', &
         ptr_patch=nf%m_froot_to_fire_patch, default='inactive')

    nf%m_froot_xfer_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_FROOT'//trim(nf%name)//'_XFER_TO_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='fine root ' // trim(nf%name) // ' transfer fire loss', &
         ptr_patch=nf%m_froot_xfer_to_fire_patch, default='inactive')

    nf%m_leaf_storage_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_LEAF' // trim(nf%name)//'_STORAGE_TO_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='leaf ' // trim(nf%name) // ' storage fire loss', &
         ptr_patch=nf%m_leaf_storage_to_fire_patch, default='inactive')

    nf%m_leaf_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_LEAF' // trim(nf%name)//'_TO_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='leaf ' // trim(nf%name) // ' fire loss', &
         ptr_patch=nf%m_leaf_to_fire_patch, default='inactive')

    nf%m_leaf_xfer_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_LEAF' // trim(nf%name)//'_XFER_TO_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='leaf ' // trim(nf%name) // ' transfer fire loss', &
         ptr_patch=nf%m_leaf_xfer_to_fire_patch, default='inactive')

    nf%m_livecroot_storage_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_LIVECROOT'//trim(name)//'_STORAGE_TO_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='live root ' // trim(nf%name) // ' storage fire loss', &
         ptr_patch=nf%m_livecroot_storage_to_fire_patch, default='inactive')

    nf%m_livecroot_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_LIVECROOT'//trim(name)//'_TO_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='live root ' // trim(nf%name) // ' fire loss', &
         ptr_patch=nf%m_livecroot_to_fire_patch, default='inactive')

    nf%m_livecroot_xfer_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_LIVECROOT'//trim(name)//'_XFER_TO_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='live root ' // trim(nf%name) // ' transfer fire loss', &
         ptr_patch=nf%m_livecroot_xfer_to_fire_patch, default='inactive')

    nf%m_livestem_storage_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_LIVESTEM' // trim(nf%name)//'_STORAGE_TO_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='live stem ' // trim(nf%name) // ' storage fire loss', &
         ptr_patch=nf%m_livestem_storage_to_fire_patch, default='inactive')

    nf%m_livestem_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_LIVESTEM' // trim(nf%name)//'_TO_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='live stem ' // trim(nf%name) // ' fire loss', &
         ptr_patch=nf%m_livestem_to_fire_patch, default='inactive')

    nf%m_livestem_xfer_to_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_LIVESTEM' // trim(nf%name)//'_XFER_TO_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='live stem ' // trim(nf%name) // ' transfer fire loss', &
         ptr_patch=nf%m_livestem_xfer_to_fire_patch, default='inactive')

    nf%m_livecroot_to_deadcroot_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_LIVECROOT'//trim(name)//'_TO_DEADCROOT'//trim(name)//'_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='live root ' // trim(nf%name) // ' fire mortality to dead root C', &
         ptr_patch=nf%m_livecroot_to_deadcroot_fire_patch, default='inactive')

    nf%m_livestem_to_deadstem_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_LIVESTEM' // trim(nf%name)//'_TO_DEADSTEM' // trim(nf%name)//'_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='live stem ' // trim(nf%name) // ' fire mortality to dead stem C', &
         ptr_patch=nf%m_livestem_to_deadstem_fire_patch, default='inactive')

    nf%m_pool_to_litter_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_'//trim(nf%name)//'POOL_TO_LITTER_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='cpool fire mortality to litter', &
         ptr_patch=nf%m_pool_to_litter_fire_patch, default='inactive')

    nf%m_deadcroot_storage_to_litter_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_DEADCROOT'//trim(name)//'_STORAGE_TO_LITTER_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='dead root ' // trim(nf%name) // ' storage fire mortality to litter', &
         ptr_patch=nf%m_deadcroot_storage_to_litter_fire_patch, default='inactive')

    nf%m_deadcroot_to_litter_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_DEADCROOT'//trim(name)//'_TO_LITTER_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='dead root ' // trim(nf%name) // ' fire mortality to litter', &
         ptr_patch=nf%m_deadcroot_to_litter_fire_patch, default='inactive')

    nf%m_deadcroot_xfer_to_litter_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_DEADCROOT'//trim(name)//'_XFER_TO_LITTER_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='dead root ' // trim(nf%name) // ' transfer fire mortality to litter', &
         ptr_patch=nf%m_deadcroot_xfer_to_litter_fire_patch, default='inactive')

    nf%m_deadstem_storage_to_litter_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_DEADSTEM' // trim(nf%name)//'_STORAGE_TO_LITTER_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='dead stem ' // trim(nf%name) // ' storage fire mortality to litter', &
         ptr_patch=nf%m_deadstem_storage_to_litter_fire_patch, default='inactive')

    nf%m_deadstem_to_litter_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_DEADSTEM' // trim(nf%name)//'_TO_LITTER_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='dead stem ' // trim(nf%name) // ' fire mortality to litter', &
         ptr_patch=nf%m_deadstem_to_litter_fire_patch, default='inactive')

    nf%m_deadstem_xfer_to_litter_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_DEADSTEM' // trim(nf%name)//'_XFER_TO_LITTER_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='dead stem ' // trim(nf%name) // ' transfer fire mortality to litter', &
         ptr_patch=nf%m_deadstem_xfer_to_litter_fire_patch, default='inactive')

    nf%m_froot_storage_to_litter_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_FROOT'//trim(nf%name)//'_STORAGE_TO_LITTER_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='fine root ' // trim(nf%name) // ' storage fire mortality to litter', &
         ptr_patch=nf%m_froot_storage_to_litter_fire_patch, default='inactive')

    nf%m_froot_to_litter_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_FROOT'//trim(nf%name)//'_TO_LITTER_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='fine root ' // trim(nf%name) // ' fire mortality to litter', &
         ptr_patch=nf%m_froot_to_litter_fire_patch, default='inactive')

    nf%m_froot_xfer_to_litter_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_FROOT'//trim(nf%name)//'_XFER_TO_LITTER_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='fine root ' // trim(nf%name) // ' transfer fire mortality to litter', &
         ptr_patch=nf%m_froot_xfer_to_litter_fire_patch, default='inactive')
    nf%m_leaf_storage_to_litter_fire_patch(begp:endp) = spval

    call hist_addfld1d (fname=trim(name_prefix)//'M_LEAF' // trim(nf%name)//'_STORAGE_TO_LITTER_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='leaf ' // trim(nf%name) // ' fire mortality to litter', &
         ptr_patch=nf%m_leaf_storage_to_litter_fire_patch, default='inactive')

    nf%m_leaf_to_litter_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_LEAF' // trim(nf%name)//'_TO_LITTER_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='leaf ' // trim(nf%name) // ' fire mortality to litter', &
         ptr_patch=nf%m_leaf_to_litter_fire_patch, default='inactive')

    nf%m_leaf_xfer_to_litter_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_LEAF' // trim(nf%name)//'_XFER_TO_LITTER_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='leaf ' // trim(nf%name) // ' transfer fire mortality to litter', &
         ptr_patch=nf%m_leaf_xfer_to_litter_fire_patch, default='inactive')

    nf%m_livecroot_storage_to_litter_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_LIVECROOT'//trim(name)//'_STORAGE_TO_LITTER_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='live root ' // trim(nf%name) // ' storage fire mortality to litter', &
         ptr_patch=nf%m_livecroot_storage_to_litter_fire_patch, default='inactive')

    nf%m_livecroot_to_litter_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_LIVECROOT'//trim(name)//'_TO_LITTER_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='live root ' // trim(nf%name) // ' fire mortality to litter', &
         ptr_patch=nf%m_livecroot_to_litter_fire_patch, default='inactive')

    nf%m_livecroot_xfer_to_litter_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_LIVECROOT'//trim(name)//'_XFER_TO_LITTER_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='live root ' // trim(nf%name) // ' transfer fire mortality to litter', &
         ptr_patch=nf%m_livecroot_xfer_to_litter_fire_patch, default='inactive')

    nf%m_livestem_storage_to_litter_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_LIVESTEM' // trim(nf%name)//'_STORAGE_TO_LITTER_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='live stem ' // trim(nf%name) // ' storage fire mortality to litter', &
         ptr_patch=nf%m_livestem_storage_to_litter_fire_patch, default='inactive')

    nf%m_livestem_to_litter_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_LIVESTEM' // trim(nf%name)//'_TO_LITTER_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='live stem ' // trim(nf%name) // ' fire mortality to litter', &
         ptr_patch=nf%m_livestem_to_litter_fire_patch, default='inactive')

    nf%m_livestem_xfer_to_litter_fire_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'M_LIVESTEM' // trim(nf%name)//'_XFER_TO_LITTER_FIRE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='live stem ' // trim(nf%name) // ' transfer fire mortality to litter', &
         ptr_patch=nf%m_livestem_xfer_to_litter_fire_patch, default='inactive')

    nf%fire_loss_col(begc:endc) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'COL_FIRE_'//trim(nf%name)//'LOSS', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', &
         long_name='total column-level fire ' // trim(nf%name) // ' loss for non-peat fires outside land-type converted region', &
         ptr_col=nf%fire_loss_col, default='inactive')

    nf%fire_loss_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'PFT_FIRE_'//trim(nf%name)//'LOSS', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', &
         long_name='total patch-level fire ' // trim(nf%name) // ' loss for non-peat fires outside land-type converted region', &
         ptr_patch=nf%fire_loss_patch)

    nf%prod100_loss_col(begc:endc) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'PROD100' // trim(nf%name)//'_LOSS', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='loss from 100-yr wood product pool', &
         ptr_col=nf%prod100_loss_col, default='inactive')

    nf%prod10_loss_col(begc:endc) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'PROD10' // trim(nf%name)//'_LOSS', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='loss from 10-yr wood product pool', &
         ptr_col=nf%prod10_loss_col, default='inactive')

    nf%prod1_loss_col(begc:endc) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'PROD1' // trim(nf%name)//'_LOSS', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='loss from 1-yr crop product pool', &
         ptr_col=nf%prod1_loss_col, default='inactive')

    nf%plant_alloc_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'PLANT_'//trim(nf%name)//'ALLOC', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='total allocated ' // trim(nf%name) // ' flux', &
         ptr_patch=nf%plant_alloc_patch, default='active')

    nf%som_leached_col(begc:endc) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'SOM_' //trim(nf%name) // '_LEACHED', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='total flux of ' // trim(nf%name) // ' from SOM pools due to leaching', &
         ptr_col=nf%som_leached_col)!, default='inactive')

    nf%wood_harvest_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'WOOD_HARVEST' // trim(nf%name), &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='wood harvest ' // trim(nf%name) // ' (to product pools)', &
         ptr_patch=nf%wood_harvest_patch)

    nf%pool_to_deadcroot_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//trim(nf%name)//'POOL_TO_DEADCROOT'//trim(name), &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='allocation to dead coarse root C', &
         ptr_patch=nf%pool_to_deadcroot_patch, default='inactive')

    nf%pool_to_deadcroot_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//trim(nf%name)//'POOL_TO_DEADCROOT'//trim(nf%name)//'_STORAGE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='allocation to dead coarse root ' // trim(nf%name) // ' storage', &
         ptr_patch=nf%pool_to_deadcroot_storage_patch, default='inactive')

    nf%pool_to_deadstem_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//trim(nf%name)//'POOL_TO_DEADSTEM' // trim(nf%name), &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='allocation to dead stem C', &
         ptr_patch=nf%pool_to_deadstem_patch, default='inactive')

    nf%pool_to_deadstem_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//trim(nf%name)//'POOL_TO_DEADSTEM' // trim(nf%name)//'_STORAGE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='allocation to dead stem ' // trim(nf%name) // ' storage', &
         ptr_patch=nf%pool_to_deadstem_storage_patch, default='inactive')

    nf%pool_to_froot_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//trim(nf%name)//'POOL_TO_FROOT'//trim(nf%name), &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='allocation to fine root C', &
         ptr_patch=nf%pool_to_froot_patch, default='inactive')

    nf%pool_to_froot_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//trim(nf%name)//'POOL_TO_FROOT'//trim(nf%name)//'_STORAGE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='allocation to fine root ' // trim(nf%name) // ' storage', &
         ptr_patch=nf%pool_to_froot_storage_patch, default='inactive')

    nf%pool_to_leaf_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//trim(nf%name)//'POOL_TO_LEAF' // trim(nf%name), &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='allocation to leaf C', &
         ptr_patch=nf%pool_to_leaf_patch, default='inactive')

    nf%pool_to_leaf_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//trim(nf%name)//'POOL_TO_LEAF' // trim(nf%name)//'_STORAGE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='allocation to leaf ' // trim(nf%name) // ' storage', &
         ptr_patch=nf%pool_to_leaf_storage_patch, default='inactive')

    nf%pool_to_livecroot_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//trim(nf%name)//'POOL_TO_LIVECROOTC', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='allocation to live coarse root C', &
         ptr_patch=nf%pool_to_livecroot_patch, default='inactive')

    nf%pool_to_livecroot_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//trim(nf%name)//'POOL_TO_LIVECROOT'//trim(nf%name)//'_STORAGE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='allocation to live coarse root ' // trim(nf%name) // ' storage', &
         ptr_patch=nf%pool_to_livecroot_storage_patch, default='inactive')

    nf%pool_to_livestem_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//trim(nf%name)//'POOL_TO_LIVESTEM' // trim(nf%name), &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='allocation to live stem C', &
         ptr_patch=nf%pool_to_livestem_patch, default='inactive')

    nf%pool_to_livestem_storage_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//trim(nf%name)//'POOL_TO_LIVESTEM' // trim(nf%name)//'_STORAGE', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='allocation to live stem ' // trim(nf%name) // ' storage', &
         ptr_patch=nf%pool_to_livestem_storage_patch, default='inactive')

    nf%livecroot_to_deadcroot_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'LIVECROOT'//trim(nf%name)//'_TO_DEADCROOT'//trim(name), &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='live coarse root ' // trim(nf%name) // ' turnover', &
         ptr_patch=nf%livecroot_to_deadcroot_patch, default='inactive')

    nf%deadcroot_xfer_to_deadcroot_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'DEADCROOT'//trim(nf%name)//'_XFER_TO_DEADCROOT'//trim(name), &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='dead coarse root ' // trim(nf%name) // ' growth from storage', &
         ptr_patch=nf%deadcroot_xfer_to_deadcroot_patch, default='inactive')
    
    nf%deadstem_xfer_to_deadstem_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'DEADSTEM' // trim(nf%name)//'_XFER_TO_DEADSTEM' // trim(nf%name), &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='dead stem ' // trim(nf%name) // ' growth from storage', &
         ptr_patch=nf%deadstem_xfer_to_deadstem_patch, default='inactive')

    nf%livecroot_xfer_to_livecroot_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'LIVECROOT'//trim(nf%name)//'_XFER_TO_LIVECROOTC', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='live coarse root ' // trim(nf%name) // ' growth from storage', &
         ptr_patch=nf%livecroot_xfer_to_livecroot_patch, default='inactive')

    nf%froot_xfer_to_froot_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'FROOT'//trim(nf%name)//'_XFER_TO_FROOT'//trim(nf%name), &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='fine root ' // trim(nf%name) // ' growth from storage', &
         ptr_patch=nf%froot_xfer_to_froot_patch, default='inactive')

    nf%leaf_xfer_to_leaf_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'LEAF' // trim(nf%name)//'_XFER_TO_LEAF' // trim(nf%name), &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='leaf ' // trim(nf%name) // ' growth from storage', &
         ptr_patch=nf%leaf_xfer_to_leaf_patch, default='inactive')
    
    nf%livestem_to_deadstem_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'LIVESTEM' // trim(nf%name)//'_TO_DEADSTEM' // trim(nf%name), &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='live stem ' // trim(nf%name) // ' turnover', &
         ptr_patch=nf%livestem_to_deadstem_patch, default='inactive')

    nf%livestem_xfer_to_livestem_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'LIVESTEM' // trim(nf%name)//'_XFER_TO_LIVESTEM' // trim(nf%name), &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='live stem ' // trim(nf%name) // ' growth from storage', &
         ptr_patch=nf%livestem_xfer_to_livestem_patch, default='inactive')

    nf%prev_leaf_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'PREV_LEAF' // trim(nf%name)//'_TO_LITTER', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='previous timestep leaf ' // trim(nf%name) // ' litterfall flux', &
         ptr_patch=nf%prev_leaf_to_litter_patch, default='inactive')

    nf%prev_froot_to_litter_patch(begp:endp) = spval
    call hist_addfld1d (fname=trim(name_prefix)//'PREV_FROOT'//trim(nf%name)//'_TO_LITTER', &
         units='g'//trim(nf%name)//'/m^2/s', &
         avgflag='A', long_name='previous timestep froot ' // trim(nf%name) // ' litterfall flux', &
         ptr_patch=nf%prev_froot_to_litter_patch, default='inactive')

  end subroutine NutrientFluxInitHistory

end module NutrientFluxType
