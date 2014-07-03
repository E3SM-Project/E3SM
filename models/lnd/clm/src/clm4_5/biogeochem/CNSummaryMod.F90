module CNSummaryMod

#ifdef CN

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNSummaryMod
!
! !DESCRIPTION:
! Module for carbon and nitrogen summary calculations
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clm_varcon, only: dzsoi_decomp, zisoi
    use pftvarcon   , only: npcropmin
    use surfrdMod   , only: crop_prog
    use abortutils  , only: endrun
    implicit none
    save
    private
! !PUBLIC MEMBER FUNCTIONS:
    public :: CSummary
    public :: NSummary
!
! !REVISION HISTORY:
! 4/23/2004: Created by Peter Thornton
! F. Li and S. Levis (11/06/12)
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CSummary
!
! !INTERFACE:
subroutine CSummary(num_soilc, filter_soilc, num_soilp, filter_soilp, isotope)
!
! !DESCRIPTION:
! On the radiation time step, perform pft and column-level carbon
! summary calculations
!
! !USES:
   use clmtype
   use pft2colMod, only: p2c
   use clm_varctl, only: iulog
   use shr_sys_mod, only: shr_sys_flush
   use clm_varpar  , only: nlevdecomp,ndecomp_pools,ndecomp_cascade_transitions
   use CNNDynamicsMod, only: nfix_timeconst
   use clm_time_manager    , only : get_step_size
   use clm_varcon      , only: secspday, spval
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
   character(len=*), intent(in) :: isotope         ! 'bulk', 'c13' or 'c14'
!
! !CALLED FROM:
! subroutine CNEcosystemDyn
!
! !REVISION HISTORY:
! 12/9/03: Created by Peter Thornton
! 11/6/12: revised by F. Li and S. Levis
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
   integer , pointer :: ivt(:)                ! pft vegetation type
   real(r8), pointer :: col_fire_closs(:)     ! (gC/m2/s) total column-level fire C loss
   real(r8), pointer :: er(:)                 ! (gC/m2/s) total ecosystem respiration, autotrophic + heterotrophic
   real(r8), pointer :: hr(:)                 ! (gC/m2/s) total heterotrophic respiration
   real(r8), pointer :: litfire(:)            ! (gC/m2/s) litter fire losses
   real(r8), pointer :: lithr(:)              ! (gC/m2/s) litter heterotrophic respiration 
   real(r8), pointer :: cwdc(:)               ! (gC/m2) coarse woody debris C
   real(r8), pointer :: col_ctrunc(:)         ! (gC/m2) column-level sink for C truncation
   real(r8), pointer :: decomp_cascade_hr_vr(:,:,:)       
   real(r8), pointer :: decomp_cascade_hr(:,:)
   real(r8), pointer :: hr_vr(:,:)       ! total vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
   real(r8), pointer :: m_decomp_cpools_to_fire_vr(:,:,:)
   real(r8), pointer :: m_decomp_cpools_to_fire(:,:)
   
   real(r8), pointer :: decomp_cpools(:,:)             ! (gC/m2)  decomposing (litter, cwd, soil) c pools
   real(r8), pointer :: decomp_cpools_1m(:,:)          ! (gC/m2)  decomposing (litter, cwd, soil) c pools to 1 meter
   real(r8), pointer :: decomp_cpools_vr(:,:,:)        ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
   integer, pointer :: altmax_indx(:)                  ! maximum annual depth of thaw
   integer, pointer :: altmax_lastyear_indx(:)         ! prior year maximum annual depth of thaw
   real(r8), pointer :: col_ctrunc_vr(:,:)         ! (gC/m3) column-level sink for C truncation
   integer,  pointer :: cascade_donor_pool(:)             ! which pool is C taken from for a given decomposition step
   logical, pointer :: is_litter(:)                       ! TRUE => pool is a litter pool
   logical, pointer :: is_soil(:)                         ! TRUE => pool is a soil pool
   logical, pointer :: is_cwd(:)                          ! TRUE => pool is a cwd pool
   real(r8), pointer :: nee(:)                ! (gC/m2/s) net ecosystem exchange of carbon, includes fire, land-use, harvest, and hrv_xsmrpool flux, positive for source
   real(r8), pointer :: nep(:)                ! (gC/m2/s) net ecosystem production, excludes fire, land-use, and harvest flux, positive for sink
   real(r8), pointer :: nbp(:)                ! (gC/m2/s) net biome production, includes fire, land-use, and harvest flux, positive for sink
   real(r8), pointer :: col_ar(:)             ! (gC/m2/s) autotrophic respiration (MR + GR)
   real(r8), pointer :: col_gpp(:)            ! GPP flux before downregulation (gC/m2/s)
   real(r8), pointer :: col_npp(:)            ! (gC/m2/s) net primary production
   real(r8), pointer :: col_lag_npp(:)        ! (gC/m2/s) lagged net primary production
   real(r8), pointer :: col_pft_fire_closs(:) ! (gC/m2/s) total pft-level fire C loss 
   real(r8), pointer :: col_litfall(:)        ! (gC/m2/s) total pft-level litterfall C loss 
   real(r8), pointer :: col_rr(:)             ! (gC/m2/s) root respiration (fine root MR + total root GR)
   real(r8), pointer :: col_vegfire(:)        ! (gC/m2/s) pft-level fire loss (obsolete, mark for removal)
   real(r8), pointer :: col_wood_harvestc(:)
   real(r8), pointer :: somfire(:)            ! (gC/m2/s) soil organic matter fire losses
   real(r8), pointer :: somhr(:)              ! (gC/m2/s) soil organic matter heterotrophic respiration
   real(r8), pointer :: sr(:)                 ! (gC/m2/s) total soil respiration (HR + root resp)
   real(r8), pointer :: totfire(:)            ! (gC/m2/s) total ecosystem fire losses
   real(r8), pointer :: col_totpftc(:)        ! (gC/m2) total pft-level carbon, including cpool
   real(r8), pointer :: col_totvegc(:)        ! (gC/m2) total vegetation carbon, excluding cpool
   real(r8), pointer :: totcolc(:)            ! (gC/m2) total column carbon, incl veg and cpool
   real(r8), pointer :: totecosysc(:)         ! (gC/m2) total ecosystem carbon, incl veg but excl cpool
   real(r8), pointer :: totlitc(:)            ! (gC/m2) total litter carbon
   real(r8), pointer :: totlitc_1m(:)         ! (gC/m2) total litter carbon to 1 meter
   real(r8), pointer :: totsomc(:)            ! (gC/m2) total soil organic matter carbon
   real(r8), pointer :: totsomc_1m(:)         ! (gC/m2) total soil organic matter carbon to 1 meter
   real(r8), pointer :: agnpp(:)              ! (gC/m2/s) aboveground NPP
   real(r8), pointer :: ar(:)                 ! (gC/m2/s) autotrophic respiration (MR + GR)
   real(r8), pointer :: bgnpp(:)                  ! (gC/m2/s) belowground NPP
   real(r8), pointer :: xsmrpool_to_atm(:)        ! excess MR pool harvest mortality (gC/m2/s)
   real(r8), pointer :: cpool_grain_gr(:)         ! grain growth respiration (gC/m2/s)
   real(r8), pointer :: cpool_grain_storage_gr(:) ! grain growth respiration to storage (gC/m2/s)
   real(r8), pointer :: cpool_to_grainc(:)        ! allocation to grain C storage (gC/m2/s)
   real(r8), pointer :: grainc_xfer_to_grainc(:)  ! grain C growth from storage (gC/m2/s)
   real(r8), pointer :: transfer_grain_gr(:)      ! grain growth respiration from storage (gC/m2/s)
   real(r8), pointer :: grainc_to_food(:)         ! grain C to food (gC/m2/s)
   real(r8), pointer :: livestemc_to_litter(:)    ! live stem C litterfall (gC/m2/s)
   real(r8), pointer :: grainc(:)                 ! (gC/m2) grain C
   real(r8), pointer :: grainc_storage(:)         ! (gC/m2) grain C storage
   real(r8), pointer :: grainc_xfer(:)            ! (gC/m2) grain C transfer
   real(r8), pointer :: cpool_deadcroot_gr(:)     ! dead coarse root growth respiration (gC/m2/s)
   real(r8), pointer :: cpool_deadcroot_storage_gr(:) ! dead coarse root growth respiration to storage (gC/m2/s)
   real(r8), pointer :: cpool_deadstem_gr(:)          ! dead stem growth respiration (gC/m2/s)
   real(r8), pointer :: cpool_deadstem_storage_gr(:)  ! dead stem growth respiration to storage (gC/m2/s)
   real(r8), pointer :: cpool_froot_gr(:)             ! fine root growth respiration (gC/m2/s)
   real(r8), pointer :: cpool_froot_storage_gr(:)     ! fine root  growth respiration to storage (gC/m2/s)
   real(r8), pointer :: cpool_leaf_gr(:)              ! leaf growth respiration (gC/m2/s)
   real(r8), pointer :: cpool_leaf_storage_gr(:)      ! leaf growth respiration to storage (gC/m2/s)
   real(r8), pointer :: cpool_livecroot_gr(:)         ! live coarse root growth respiration (gC/m2/s)
   real(r8), pointer :: cpool_livecroot_storage_gr(:) ! live coarse root growth respiration to storage (gC/m2/s)
   real(r8), pointer :: cpool_livestem_gr(:)          ! live stem growth respiration (gC/m2/s)
   real(r8), pointer :: cpool_livestem_storage_gr(:)  ! live stem growth respiration to storage (gC/m2/s)
   real(r8), pointer :: cpool_to_deadcrootc(:)        ! allocation to dead coarse root C (gC/m2/s)
   real(r8), pointer :: cpool_to_deadstemc(:)         ! allocation to dead stem C (gC/m2/s)
   real(r8), pointer :: cpool_to_frootc(:)            ! allocation to fine root C (gC/m2/s)
   real(r8), pointer :: cpool_to_leafc(:)             ! allocation to leaf C (gC/m2/s)
   real(r8), pointer :: cpool_to_livecrootc(:)        ! allocation to live coarse root C (gC/m2/s)
   real(r8), pointer :: cpool_to_livestemc(:)         ! allocation to live stem C (gC/m2/s)
   real(r8), pointer :: current_gr(:)                 ! (gC/m2/s) growth resp for new growth displayed in this timestep
   real(r8), pointer :: deadcrootc_xfer_to_deadcrootc(:)
   real(r8), pointer :: deadstemc_xfer_to_deadstemc(:) 
   real(r8), pointer :: frootc_to_litter(:)
   real(r8), pointer :: frootc_xfer_to_frootc(:)       
   real(r8), pointer :: froot_mr(:)     
   real(r8), pointer :: grain_mr(:)
   real(r8), pointer :: froot_curmr(:)     
   real(r8), pointer :: froot_xsmr(:)     
   real(r8), pointer :: gpp(:)                !GPP flux before downregulation (gC/m2/s)
   real(r8), pointer :: gr(:)                 ! (gC/m2/s) total growth respiration
   real(r8), pointer :: leafc_to_litter(:)
   real(r8), pointer :: leafc_xfer_to_leafc(:)         
   real(r8), pointer :: leaf_mr(:)
   real(r8), pointer :: leaf_curmr(:)
   real(r8), pointer :: leaf_xsmr(:)
   real(r8), pointer :: litfall(:)            ! (gC/m2/s) litterfall (leaves and fine roots)
   real(r8), pointer :: livecrootc_xfer_to_livecrootc(:)
   real(r8), pointer :: livecroot_mr(:)
   real(r8), pointer :: livecroot_curmr(:)
   real(r8), pointer :: livecroot_xsmr(:)
   real(r8), pointer :: livestemc_xfer_to_livestemc(:) 
   real(r8), pointer :: livestem_mr(:)  
   real(r8), pointer :: livestem_curmr(:)  
   real(r8), pointer :: livestem_xsmr(:)  

! fire variables changed by F. Li and S. Levis
   real(r8), pointer :: m_leafc_to_fire(:)             
   real(r8), pointer :: m_leafc_storage_to_fire(:)                
   real(r8), pointer :: m_leafc_xfer_to_fire(:)        
   real(r8), pointer :: m_livestemc_to_fire(:)         
   real(r8), pointer :: m_livestemc_storage_to_fire(:)       
   real(r8), pointer :: m_livestemc_xfer_to_fire(:)    
   real(r8), pointer :: m_deadstemc_to_fire(:)         
   real(r8), pointer :: m_deadstemc_storage_to_fire(:)          
   real(r8), pointer :: m_deadstemc_xfer_to_fire(:)    
   real(r8), pointer :: m_frootc_to_fire(:)            
   real(r8), pointer :: m_frootc_storage_to_fire(:)    
   real(r8), pointer :: m_frootc_xfer_to_fire(:)       
   real(r8), pointer :: m_livecrootc_to_fire(:)        
   real(r8), pointer :: m_livecrootc_storage_to_fire(:)     
   real(r8), pointer :: m_livecrootc_xfer_to_fire(:)   
   real(r8), pointer :: m_deadcrootc_to_fire(:)        
   real(r8), pointer :: m_deadcrootc_storage_to_fire(:)
   real(r8), pointer :: m_deadcrootc_xfer_to_fire(:)  
   real(r8), pointer :: m_gresp_storage_to_fire(:)     
   real(r8), pointer :: m_gresp_xfer_to_fire(:)       
   real(r8), pointer :: m_leafc_to_litter_fire(:)   
   real(r8), pointer :: m_leafc_storage_to_litter_fire(:)                
   real(r8), pointer :: m_leafc_xfer_to_litter_fire(:)  
   real(r8), pointer :: m_livestemc_to_litter_fire(:)    
   real(r8), pointer :: m_livestemc_storage_to_litter_fire(:)        
   real(r8), pointer :: m_livestemc_xfer_to_litter_fire(:) 
   real(r8), pointer :: m_livestemc_to_deadstemc_fire(:)    
   real(r8), pointer :: m_deadstemc_to_litter_fire(:) 
   real(r8), pointer :: m_deadstemc_storage_to_litter_fire(:)           
   real(r8), pointer :: m_deadstemc_xfer_to_litter_fire(:) 
   real(r8), pointer :: m_frootc_to_litter_fire(:)        
   real(r8), pointer :: m_frootc_storage_to_litter_fire(:)  
   real(r8), pointer :: m_frootc_xfer_to_litter_fire(:)
   real(r8), pointer :: m_livecrootc_to_litter_fire(:)    
   real(r8), pointer :: m_livecrootc_storage_to_litter_fire(:)      
   real(r8), pointer :: m_livecrootc_xfer_to_litter_fire(:)
   real(r8), pointer :: m_livecrootc_to_deadcrootc_fire(:)    
   real(r8), pointer :: m_deadcrootc_to_litter_fire(:)        
   real(r8), pointer :: m_deadcrootc_storage_to_litter_fire(:)  
   real(r8), pointer :: m_deadcrootc_xfer_to_litter_fire(:)
   real(r8), pointer :: m_gresp_storage_to_litter_fire(:)      
   real(r8), pointer :: m_gresp_xfer_to_litter_fire(:)    

   real(r8), pointer :: m_deadcrootc_storage_to_litter(:)
   real(r8), pointer :: m_deadcrootc_to_litter(:)
   real(r8), pointer :: m_deadcrootc_xfer_to_litter(:)
   real(r8), pointer :: m_deadstemc_storage_to_litter(:)
   real(r8), pointer :: m_deadstemc_to_litter(:)
   real(r8), pointer :: m_deadstemc_xfer_to_litter(:)
   real(r8), pointer :: m_frootc_storage_to_litter(:)
   real(r8), pointer :: m_frootc_to_litter(:)
   real(r8), pointer :: m_frootc_xfer_to_litter(:)
   real(r8), pointer :: m_gresp_storage_to_litter(:)
   real(r8), pointer :: m_gresp_xfer_to_litter(:)
   real(r8), pointer :: m_leafc_storage_to_litter(:)
   real(r8), pointer :: m_leafc_to_litter(:)
   real(r8), pointer :: m_leafc_xfer_to_litter(:)
   real(r8), pointer :: m_livecrootc_storage_to_litter(:)
   real(r8), pointer :: m_livecrootc_to_litter(:)
   real(r8), pointer :: m_livecrootc_xfer_to_litter(:)
   real(r8), pointer :: m_livestemc_storage_to_litter(:)
   real(r8), pointer :: m_livestemc_to_litter(:)
   real(r8), pointer :: m_livestemc_xfer_to_litter(:)
   real(r8), pointer :: hrv_leafc_to_litter(:)              
   real(r8), pointer :: hrv_leafc_storage_to_litter(:)      
   real(r8), pointer :: hrv_leafc_xfer_to_litter(:)         
   real(r8), pointer :: hrv_frootc_to_litter(:)             
   real(r8), pointer :: hrv_frootc_storage_to_litter(:)     
   real(r8), pointer :: hrv_frootc_xfer_to_litter(:)        
   real(r8), pointer :: hrv_livestemc_to_litter(:)          
   real(r8), pointer :: hrv_livestemc_storage_to_litter(:)  
   real(r8), pointer :: hrv_livestemc_xfer_to_litter(:)     
   real(r8), pointer :: hrv_deadstemc_to_prod10c(:)         
   real(r8), pointer :: hrv_deadstemc_to_prod100c(:)        
   real(r8), pointer :: hrv_deadstemc_storage_to_litter(:)  
   real(r8), pointer :: hrv_deadstemc_xfer_to_litter(:)     
   real(r8), pointer :: hrv_livecrootc_to_litter(:)         
   real(r8), pointer :: hrv_livecrootc_storage_to_litter(:) 
   real(r8), pointer :: hrv_livecrootc_xfer_to_litter(:)    
   real(r8), pointer :: hrv_deadcrootc_to_litter(:)         
   real(r8), pointer :: hrv_deadcrootc_storage_to_litter(:) 
   real(r8), pointer :: hrv_deadcrootc_xfer_to_litter(:)    
   real(r8), pointer :: hrv_gresp_storage_to_litter(:)      
   real(r8), pointer :: hrv_gresp_xfer_to_litter(:)         
   real(r8), pointer :: hrv_xsmrpool_to_atm(:)              
   real(r8), pointer :: col_hrv_xsmrpool_to_atm(:)              
   real(r8), pointer :: mr(:)                 ! (gC/m2/s) maintenance respiration
   real(r8), pointer :: npp(:)                ! (gC/m2/s) net primary production
   real(r8), pointer :: pft_fire_closs(:)     ! (gC/m2/s) total pft-level fire C loss 
   real(r8), pointer :: psnshade_to_cpool(:)
   real(r8), pointer :: psnsun_to_cpool(:) 
   real(r8), pointer :: rr(:)                 ! (gC/m2/s) root respiration (fine root MR + total root GR)
   real(r8), pointer :: storage_gr(:)         ! (gC/m2/s) growth resp for growth sent to storage for later display
   real(r8), pointer :: transfer_deadcroot_gr(:)
   real(r8), pointer :: transfer_deadstem_gr(:)      
   real(r8), pointer :: transfer_froot_gr(:)         
   real(r8), pointer :: transfer_gr(:)        ! (gC/m2/s) growth resp for transfer growth displayed in this timestep
   real(r8), pointer :: transfer_leaf_gr(:)          
   real(r8), pointer :: transfer_livecroot_gr(:)     
   real(r8), pointer :: transfer_livestem_gr(:)      
   real(r8), pointer :: wood_harvestc(:)      ! (gC/m2/s) pft-level wood harvest (to product pools)
   real(r8), pointer :: vegfire(:)            ! (gC/m2/s) pft-level fire loss (obsolete, mark for removal)
   real(r8), pointer :: cpool(:)              ! (gC/m2) temporary photosynthate C pool
   real(r8), pointer :: xsmrpool(:)           ! (gC/m2) temporary photosynthate C pool
   real(r8), pointer :: pft_ctrunc(:)         ! (gC/m2) pft-level sink for C truncation
   real(r8), pointer :: deadcrootc(:)         ! (gC/m2) dead coarse root C
   real(r8), pointer :: deadcrootc_storage(:) ! (gC/m2) dead coarse root C storage
   real(r8), pointer :: deadcrootc_xfer(:)    ! (gC/m2) dead coarse root C transfer
   real(r8), pointer :: deadstemc(:)          ! (gC/m2) dead stem C
   real(r8), pointer :: deadstemc_storage(:)  ! (gC/m2) dead stem C storage
   real(r8), pointer :: deadstemc_xfer(:)     ! (gC/m2) dead stem C transfer
   real(r8), pointer :: dispvegc(:)           ! (gC/m2) displayed veg carbon, excluding storage and cpool
   real(r8), pointer :: frootc(:)             ! (gC/m2) fine root C
   real(r8), pointer :: frootc_storage(:)     ! (gC/m2) fine root C storage
   real(r8), pointer :: frootc_xfer(:)        ! (gC/m2) fine root C transfer
   real(r8), pointer :: gresp_storage(:)      ! (gC/m2) growth respiration storage
   real(r8), pointer :: gresp_xfer(:)         ! (gC/m2) growth respiration transfer
   real(r8), pointer :: leafc(:)              ! (gC/m2) leaf C
   real(r8), pointer :: leafc_storage(:)      ! (gC/m2) leaf C storage
   real(r8), pointer :: leafc_xfer(:)         ! (gC/m2) leaf C transfer
   real(r8), pointer :: livecrootc(:)         ! (gC/m2) live coarse root C
   real(r8), pointer :: livecrootc_storage(:) ! (gC/m2) live coarse root C storage
   real(r8), pointer :: livecrootc_xfer(:)    ! (gC/m2) live coarse root C transfer
   real(r8), pointer :: livestemc(:)          ! (gC/m2) live stem C
   real(r8), pointer :: livestemc_storage(:)  ! (gC/m2) live stem C storage
   real(r8), pointer :: livestemc_xfer(:)     ! (gC/m2) live stem C transfer
   real(r8), pointer :: storvegc(:)           ! (gC/m2) stored vegetation carbon, excluding cpool
   real(r8), pointer :: totpftc(:)            ! (gC/m2) total pft-level carbon, including cpool
   real(r8), pointer :: totvegc(:)            ! (gC/m2) total vegetation carbon, excluding cpool
   real(r8), pointer :: tempsum_npp(:)        ! temporary annual sum of NPP (gC/m2/yr)
#if (defined CNDV)
   real(r8), pointer :: tempsum_litfall(:)      !temporary annual sum of litfall (gC/m2/yr)
#endif
   ! for landcover change
   real(r8), pointer :: landuseflux(:)        ! (gC/m2/s) dwt_closs+product_closs
   real(r8), pointer :: landuptake(:)         ! (gC/m2/s) nee-landuseflux
   real(r8), pointer :: dwt_closs(:)          ! (gC/m2/s) total carbon loss from land cover conversion
   real(r8), pointer :: dwt_conv_cflux(:)     ! (gC/m2/s) conversion C flux (immediate loss to atm)
   real(r8), pointer :: prod10c_loss(:)       ! (gC/m2/s) loss from 10-yr wood product pool
   real(r8), pointer :: prod100c_loss(:)      ! (gC/m2/s) loss from 100-yr wood product pool
   real(r8), pointer :: product_closs(:)      ! (gC/m2/s) total wood product carbon loss
   real(r8), pointer :: seedc(:)              ! (gC/m2) column-level pool for seeding new PFTs
   real(r8), pointer :: prod10c(:)            ! (gC/m2) wood product C pool, 10-year lifespan
   real(r8), pointer :: prod100c(:)           ! (gC/m2) wood product C pool, 100-year lifespan
   real(r8), pointer :: totprodc(:)           ! (gC/m2) total wood product C

   real(r8), pointer :: frootc_alloc(:)       ! fine root C allocation (gC/m2/s)
   real(r8), pointer :: frootc_loss(:)        ! fine root C loss (gC/m2/s)
   real(r8), pointer :: leafc_alloc(:)        ! leaf C allocation (gC/m2/s)
   real(r8), pointer :: leafc_loss(:)         ! leaf C loss (gC/m2/s)
   real(r8), pointer :: woodc(:)              ! wood C (gC/m2)
   real(r8), pointer :: woodc_alloc(:)        ! wood C allocation (gC/m2/s)
   real(r8), pointer :: woodc_loss(:)         ! wood C loss (gC/m2/s)
   real(r8), pointer :: cwdc_hr(:)            ! coarse woody debris C heterotrophic respiration (gC/m2/s)
   real(r8), pointer :: cwdc_loss(:)          ! coarse woody debris C loss (gC/m2/s)
   real(r8), pointer :: litterc_loss(:)       ! litter C loss (gC/m2/s)
   real(r8), pointer :: decomp_cascade_ctransfer_vr(:,:,:)  
   real(r8), pointer :: decomp_cascade_ctransfer(:,:)
   real(r8), pointer :: som_c_leached(:)                           ! total SOM C loss from vertical transport (gC/m^2/s)
   real(r8), pointer :: decomp_cpools_leached(:,:)                 ! C loss from vertical transport from each decomposing C pool (gC/m^2/s)
   real(r8), pointer :: decomp_cpools_transport_tendency(:,:,:)    ! C tendency due to vertical transport in decomposing C pools (gC/m^3/s)
   real(r8) :: nfixlags, dtime                ! temp variables for making lagged npp

!
!
! local pointers to implicit in/out scalars
!
!
! local pointers to implicit out scalars
!
!
! !OTHER LOCAL VARIABLES:
   type(pft_cflux_type), pointer :: pcisof
   type(pft_cstate_type), pointer :: pcisos
   type(column_cflux_type), pointer :: ccisof
   type(column_cstate_type), pointer :: ccisos
   integer :: c,p,j,k,l        ! indices
   integer :: fp,fc        ! lake filter indices
   real(r8) :: maxdepth    ! depth to integrate soil variables

!EOP
!-----------------------------------------------------------------------
   ! select which isotope
   select case (isotope)
   case ('bulk')
      pcisof => clm3%g%l%c%p%pcf
      pcisos => clm3%g%l%c%p%pcs
      ccisof => clm3%g%l%c%ccf
      ccisos => clm3%g%l%c%ccs
   case ('c14')
      pcisof => clm3%g%l%c%p%pc14f
      pcisos => clm3%g%l%c%p%pc14s
      ccisof => clm3%g%l%c%cc14f
      ccisos => clm3%g%l%c%cc14s
   case ('c13')
      pcisof => clm3%g%l%c%p%pc13f
      pcisos => clm3%g%l%c%p%pc13s
      ccisof => clm3%g%l%c%cc13f
      ccisos => clm3%g%l%c%cc13s
   case default
      call endrun('CNCIsoSummaryMod: iso must be bulk, c13 or c14')
   end select

   ! assign local pointers
    ivt                            => clm3%g%l%c%p%itype
    col_fire_closs                 => ccisof%col_fire_closs
    er                             => ccisof%er
    hr                             => ccisof%hr
    litfire                        => ccisof%litfire
    lithr                          => ccisof%lithr
    col_totpftc                    => ccisos%pcs_a%totpftc
    col_totvegc                    => ccisos%pcs_a%totvegc
    cwdc                           => ccisos%cwdc
    col_ctrunc                     => ccisos%col_ctrunc
    decomp_cascade_hr_vr              => ccisof%decomp_cascade_hr_vr
    decomp_cascade_hr                 => ccisof%decomp_cascade_hr
    hr_vr                             => ccisof%hr_vr
    m_decomp_cpools_to_fire_vr        => ccisof%m_decomp_cpools_to_fire_vr
    m_decomp_cpools_to_fire           => ccisof%m_decomp_cpools_to_fire
    decomp_cascade_ctransfer_vr       => ccisof%decomp_cascade_ctransfer_vr
    decomp_cascade_ctransfer          => ccisof%decomp_cascade_ctransfer
    decomp_cpools_vr                  => ccisos%decomp_cpools_vr
    decomp_cpools                     => ccisos%decomp_cpools
    decomp_cpools_1m                  => ccisos%decomp_cpools_1m
    altmax_indx                       => clm3%g%l%c%cps%altmax_indx
    altmax_lastyear_indx              => clm3%g%l%c%cps%altmax_lastyear_indx
    col_ctrunc_vr                     => ccisos%col_ctrunc_vr
    cascade_donor_pool                => decomp_cascade_con%cascade_donor_pool
    is_litter                         => decomp_cascade_con%is_litter
    is_soil                           => decomp_cascade_con%is_soil
    is_cwd                            => decomp_cascade_con%is_cwd
    nee                            => ccisof%nee
    nep                            => ccisof%nep
    nbp                            => ccisof%nbp
    col_ar                         => ccisof%pcf_a%ar
    col_gpp                        => ccisof%pcf_a%gpp
    col_npp                        => ccisof%pcf_a%npp
    col_pft_fire_closs             => ccisof%pcf_a%pft_fire_closs
    col_litfall                    => ccisof%pcf_a%litfall
    col_rr                         => ccisof%pcf_a%rr
    col_vegfire                    => ccisof%pcf_a%vegfire
    col_wood_harvestc              => ccisof%pcf_a%wood_harvestc
    somfire                        => ccisof%somfire
    somhr                          => ccisof%somhr
    sr                             => ccisof%sr
    totfire                        => ccisof%totfire
    cwdc_hr                        => ccisof%cwdc_hr
    cwdc_loss                      => ccisof%cwdc_loss
    litterc_loss                   => ccisof%litterc_loss
    ! dynamic landcover pointers
    dwt_closs                      => ccisof%dwt_closs
    landuseflux                    => ccisof%landuseflux
    landuptake                     => ccisof%landuptake
    dwt_conv_cflux                 => ccisof%dwt_conv_cflux
    seedc                          => ccisos%seedc
    
    ! wood product pointers
    prod10c_loss                   => ccisof%prod10c_loss
    prod100c_loss                  => ccisof%prod100c_loss
    product_closs                  => ccisof%product_closs
    prod10c                        => ccisos%prod10c
    prod100c                       => ccisos%prod100c
    totprodc                       => ccisos%totprodc
    
    totcolc                        => ccisos%totcolc
    totecosysc                     => ccisos%totecosysc
    totlitc                        => ccisos%totlitc
    totlitc_1m                     => ccisos%totlitc_1m
    totsomc                        => ccisos%totsomc
    totsomc_1m                     => ccisos%totsomc_1m
    agnpp                          => pcisof%agnpp
    ar                             => pcisof%ar
    bgnpp                          => pcisof%bgnpp
    xsmrpool_to_atm                => pcisof%xsmrpool_to_atm
    cpool_grain_gr                 => pcisof%cpool_grain_gr
    cpool_grain_storage_gr         => pcisof%cpool_grain_storage_gr
    cpool_to_grainc                => pcisof%cpool_to_grainc
    grainc_xfer_to_grainc          => pcisof%grainc_xfer_to_grainc
    transfer_grain_gr              => pcisof%transfer_grain_gr
    grainc_to_food                 => pcisof%grainc_to_food
    livestemc_to_litter            => pcisof%livestemc_to_litter
    grainc                         => pcisos%grainc
    grainc_storage                 => pcisos%grainc_storage
    grainc_xfer                    => pcisos%grainc_xfer
    cpool_deadcroot_gr             => pcisof%cpool_deadcroot_gr
    cpool_deadcroot_storage_gr     => pcisof%cpool_deadcroot_storage_gr
    cpool_deadstem_gr              => pcisof%cpool_deadstem_gr
    cpool_deadstem_storage_gr      => pcisof%cpool_deadstem_storage_gr
    cpool_froot_gr                 => pcisof%cpool_froot_gr
    cpool_froot_storage_gr         => pcisof%cpool_froot_storage_gr
    cpool_leaf_gr                  => pcisof%cpool_leaf_gr
    cpool_leaf_storage_gr          => pcisof%cpool_leaf_storage_gr
    cpool_livecroot_gr             => pcisof%cpool_livecroot_gr
    cpool_livecroot_storage_gr     => pcisof%cpool_livecroot_storage_gr
    cpool_livestem_gr              => pcisof%cpool_livestem_gr
    cpool_livestem_storage_gr      => pcisof%cpool_livestem_storage_gr
    cpool_to_deadcrootc            => pcisof%cpool_to_deadcrootc
    cpool_to_deadstemc             => pcisof%cpool_to_deadstemc
    cpool_to_frootc                => pcisof%cpool_to_frootc
    cpool_to_leafc                 => pcisof%cpool_to_leafc
    cpool_to_livecrootc            => pcisof%cpool_to_livecrootc
    cpool_to_livestemc             => pcisof%cpool_to_livestemc
    current_gr                     => pcisof%current_gr
    deadcrootc_xfer_to_deadcrootc  => pcisof%deadcrootc_xfer_to_deadcrootc
    deadstemc_xfer_to_deadstemc    => pcisof%deadstemc_xfer_to_deadstemc
    frootc_to_litter               => pcisof%frootc_to_litter
    frootc_xfer_to_frootc          => pcisof%frootc_xfer_to_frootc
    froot_mr                       => pcisof%froot_mr
    froot_curmr                    => pcisof%froot_curmr
    froot_xsmr                     => pcisof%froot_xsmr
    grain_mr                       => pcisof%grain_mr
    gpp                            => pcisof%gpp
    gr                             => pcisof%gr
    leafc_to_litter                => pcisof%leafc_to_litter
    leafc_xfer_to_leafc            => pcisof%leafc_xfer_to_leafc
    leaf_mr                        => pcisof%leaf_mr
    leaf_curmr                     => pcisof%leaf_curmr
    leaf_xsmr                      => pcisof%leaf_xsmr
    litfall                        => pcisof%litfall
    livecrootc_xfer_to_livecrootc  => pcisof%livecrootc_xfer_to_livecrootc
    livecroot_mr                   => pcisof%livecroot_mr
    livecroot_curmr                => pcisof%livecroot_curmr
    livecroot_xsmr                 => pcisof%livecroot_xsmr
    livestemc_xfer_to_livestemc    => pcisof%livestemc_xfer_to_livestemc
    livestem_mr                    => pcisof%livestem_mr
    livestem_curmr                 => pcisof%livestem_curmr
    livestem_xsmr                  => pcisof%livestem_xsmr

! fire variables changed by F. Li and S. Levis
    m_leafc_to_fire                => pcisof%m_leafc_to_fire
    m_leafc_storage_to_fire        => pcisof%m_leafc_storage_to_fire
    m_leafc_xfer_to_fire           => pcisof%m_leafc_xfer_to_fire
    m_livestemc_to_fire            => pcisof%m_livestemc_to_fire
    m_livestemc_storage_to_fire    => pcisof%m_livestemc_storage_to_fire
    m_livestemc_xfer_to_fire       => pcisof%m_livestemc_xfer_to_fire
    m_deadstemc_to_fire            => pcisof%m_deadstemc_to_fire
    m_deadstemc_storage_to_fire    => pcisof%m_deadstemc_storage_to_fire
    m_deadstemc_xfer_to_fire       => pcisof%m_deadstemc_xfer_to_fire
    m_frootc_to_fire               => pcisof%m_frootc_to_fire
    m_frootc_storage_to_fire       => pcisof%m_frootc_storage_to_fire
    m_frootc_xfer_to_fire          => pcisof%m_frootc_xfer_to_fire
    m_livecrootc_to_fire           => pcisof%m_livecrootc_to_fire
    m_livecrootc_storage_to_fire   => pcisof%m_livecrootc_storage_to_fire
    m_livecrootc_xfer_to_fire      => pcisof%m_livecrootc_xfer_to_fire
    m_deadcrootc_to_fire           => pcisof%m_deadcrootc_to_fire
    m_deadcrootc_storage_to_fire   => pcisof%m_deadcrootc_storage_to_fire
    m_deadcrootc_xfer_to_fire      => pcisof%m_deadcrootc_xfer_to_fire
    m_gresp_storage_to_fire        => pcisof%m_gresp_storage_to_fire
    m_gresp_xfer_to_fire           => pcisof%m_gresp_xfer_to_fire
    m_leafc_to_litter_fire               => pcisof%m_leafc_to_litter_fire
    m_leafc_storage_to_litter_fire      => pcisof%m_leafc_storage_to_litter_fire
    m_leafc_xfer_to_litter_fire         => pcisof%m_leafc_xfer_to_litter_fire
    m_livestemc_to_litter_fire          => pcisof%m_livestemc_to_litter_fire
    m_livestemc_storage_to_litter_fire  => pcisof%m_livestemc_storage_to_litter_fire
    m_livestemc_xfer_to_litter_fire     => pcisof%m_livestemc_xfer_to_litter_fire
    m_livestemc_to_deadstemc_fire       => pcisof%m_livestemc_to_deadstemc_fire
    m_deadstemc_to_litter_fire          => pcisof%m_deadstemc_to_litter_fire
    m_deadstemc_storage_to_litter_fire  => pcisof%m_deadstemc_storage_to_litter_fire
    m_deadstemc_xfer_to_litter_fire     => pcisof%m_deadstemc_xfer_to_litter_fire
    m_frootc_to_litter_fire             => pcisof%m_frootc_to_litter_fire
    m_frootc_storage_to_litter_fire     => pcisof%m_frootc_storage_to_litter_fire
    m_frootc_xfer_to_litter_fire        => pcisof%m_frootc_xfer_to_litter_fire
    m_livecrootc_to_litter_fire         => pcisof%m_livecrootc_to_litter_fire
    m_livecrootc_storage_to_litter_fire => pcisof%m_livecrootc_storage_to_litter_fire
    m_livecrootc_xfer_to_litter_fire    => pcisof%m_livecrootc_xfer_to_litter_fire
    m_livecrootc_to_deadcrootc_fire     => pcisof%m_livecrootc_to_deadcrootc_fire
    m_deadcrootc_to_litter_fire         => pcisof%m_deadcrootc_to_litter_fire
    m_deadcrootc_storage_to_litter_fire => pcisof%m_deadcrootc_storage_to_litter_fire
    m_deadcrootc_xfer_to_litter_fire    => pcisof%m_deadcrootc_xfer_to_litter_fire
    m_gresp_storage_to_litter_fire      => pcisof%m_gresp_storage_to_litter_fire
    m_gresp_xfer_to_litter_fire         => pcisof%m_gresp_xfer_to_litter_fire

    m_deadcrootc_storage_to_litter => pcisof%m_deadcrootc_storage_to_litter
    m_deadcrootc_to_litter         => pcisof%m_deadcrootc_to_litter
    m_deadcrootc_xfer_to_litter    => pcisof%m_deadcrootc_xfer_to_litter
    m_deadstemc_storage_to_litter  => pcisof%m_deadstemc_storage_to_litter
    m_deadstemc_to_litter          => pcisof%m_deadstemc_to_litter
    m_deadstemc_xfer_to_litter     => pcisof%m_deadstemc_xfer_to_litter
    m_frootc_storage_to_litter     => pcisof%m_frootc_storage_to_litter
    m_frootc_to_litter             => pcisof%m_frootc_to_litter
    m_frootc_xfer_to_litter        => pcisof%m_frootc_xfer_to_litter
    m_gresp_storage_to_litter      => pcisof%m_gresp_storage_to_litter
    m_gresp_xfer_to_litter         => pcisof%m_gresp_xfer_to_litter
    m_leafc_storage_to_litter      => pcisof%m_leafc_storage_to_litter
    m_leafc_to_litter              => pcisof%m_leafc_to_litter
    m_leafc_xfer_to_litter         => pcisof%m_leafc_xfer_to_litter
    m_livecrootc_storage_to_litter => pcisof%m_livecrootc_storage_to_litter
    m_livecrootc_to_litter         => pcisof%m_livecrootc_to_litter
    m_livecrootc_xfer_to_litter    => pcisof%m_livecrootc_xfer_to_litter
    m_livestemc_storage_to_litter  => pcisof%m_livestemc_storage_to_litter
    m_livestemc_to_litter          => pcisof%m_livestemc_to_litter
    m_livestemc_xfer_to_litter     => pcisof%m_livestemc_xfer_to_litter
    hrv_leafc_to_litter               => pcisof%hrv_leafc_to_litter               
    hrv_leafc_storage_to_litter       => pcisof%hrv_leafc_storage_to_litter     
    hrv_leafc_xfer_to_litter          => pcisof%hrv_leafc_xfer_to_litter        
    hrv_frootc_to_litter              => pcisof%hrv_frootc_to_litter            
    hrv_frootc_storage_to_litter      => pcisof%hrv_frootc_storage_to_litter    
    hrv_frootc_xfer_to_litter         => pcisof%hrv_frootc_xfer_to_litter       
    hrv_livestemc_to_litter           => pcisof%hrv_livestemc_to_litter         
    hrv_livestemc_storage_to_litter   => pcisof%hrv_livestemc_storage_to_litter 
    hrv_livestemc_xfer_to_litter      => pcisof%hrv_livestemc_xfer_to_litter    
    hrv_deadstemc_to_prod10c          => pcisof%hrv_deadstemc_to_prod10c        
    hrv_deadstemc_to_prod100c         => pcisof%hrv_deadstemc_to_prod100c       
    hrv_deadstemc_storage_to_litter   => pcisof%hrv_deadstemc_storage_to_litter 
    hrv_deadstemc_xfer_to_litter      => pcisof%hrv_deadstemc_xfer_to_litter    
    hrv_livecrootc_to_litter          => pcisof%hrv_livecrootc_to_litter        
    hrv_livecrootc_storage_to_litter  => pcisof%hrv_livecrootc_storage_to_litter
    hrv_livecrootc_xfer_to_litter     => pcisof%hrv_livecrootc_xfer_to_litter   
    hrv_deadcrootc_to_litter          => pcisof%hrv_deadcrootc_to_litter        
    hrv_deadcrootc_storage_to_litter  => pcisof%hrv_deadcrootc_storage_to_litter
    hrv_deadcrootc_xfer_to_litter     => pcisof%hrv_deadcrootc_xfer_to_litter   
    hrv_gresp_storage_to_litter       => pcisof%hrv_gresp_storage_to_litter     
    hrv_gresp_xfer_to_litter          => pcisof%hrv_gresp_xfer_to_litter        
    hrv_xsmrpool_to_atm               => pcisof%hrv_xsmrpool_to_atm             
    col_hrv_xsmrpool_to_atm           => ccisof%pcf_a%hrv_xsmrpool_to_atm             
    mr                             => pcisof%mr
    npp                            => pcisof%npp
    pft_fire_closs                 => pcisof%pft_fire_closs
    psnshade_to_cpool              => pcisof%psnshade_to_cpool
    psnsun_to_cpool                => pcisof%psnsun_to_cpool
    rr                             => pcisof%rr
    storage_gr                     => pcisof%storage_gr
    transfer_deadcroot_gr          => pcisof%transfer_deadcroot_gr
    transfer_deadstem_gr           => pcisof%transfer_deadstem_gr
    transfer_froot_gr              => pcisof%transfer_froot_gr
    transfer_gr                    => pcisof%transfer_gr
    transfer_leaf_gr               => pcisof%transfer_leaf_gr
    transfer_livecroot_gr          => pcisof%transfer_livecroot_gr
    transfer_livestem_gr           => pcisof%transfer_livestem_gr
    vegfire                        => pcisof%vegfire
    wood_harvestc                  => pcisof%wood_harvestc
    frootc_alloc                   => pcisof%frootc_alloc
    frootc_loss                    => pcisof%frootc_loss
    leafc_alloc                    => pcisof%leafc_alloc
    leafc_loss                     => pcisof%leafc_loss
    woodc_alloc                    => pcisof%woodc_alloc
    woodc_loss                     => pcisof%woodc_loss
    cpool                          => pcisos%cpool
    xsmrpool                       => pcisos%xsmrpool
    pft_ctrunc                     => pcisos%pft_ctrunc
    deadcrootc                     => pcisos%deadcrootc
    deadcrootc_storage             => pcisos%deadcrootc_storage
    deadcrootc_xfer                => pcisos%deadcrootc_xfer
    deadstemc                      => pcisos%deadstemc
    deadstemc_storage              => pcisos%deadstemc_storage
    deadstemc_xfer                 => pcisos%deadstemc_xfer
    dispvegc                       => pcisos%dispvegc
    frootc                         => pcisos%frootc
    frootc_storage                 => pcisos%frootc_storage
    frootc_xfer                    => pcisos%frootc_xfer
    gresp_storage                  => pcisos%gresp_storage
    gresp_xfer                     => pcisos%gresp_xfer
    leafc                          => pcisos%leafc
    leafc_storage                  => pcisos%leafc_storage
    leafc_xfer                     => pcisos%leafc_xfer
    livecrootc                     => pcisos%livecrootc
    livecrootc_storage             => pcisos%livecrootc_storage
    livecrootc_xfer                => pcisos%livecrootc_xfer
    livestemc                      => pcisos%livestemc
    livestemc_storage              => pcisos%livestemc_storage
    livestemc_xfer                 => pcisos%livestemc_xfer
    storvegc                       => pcisos%storvegc
    totpftc                        => pcisos%totpftc
    totvegc                        => pcisos%totvegc
    woodc                          => pcisos%woodc
    tempsum_npp                    => clm3%g%l%c%p%pepv%tempsum_npp
#if (defined CNDV)
    tempsum_litfall                => clm3%g%l%c%p%pepv%tempsum_litfall
#endif
    som_c_leached                      => ccisof%som_c_leached
    decomp_cpools_leached              => ccisof%decomp_cpools_leached
    decomp_cpools_transport_tendency   => ccisof%decomp_cpools_transport_tendency

   ! pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)

      ! calculate pft-level summary carbon fluxes and states

      ! gross primary production (GPP)
      gpp(p) = &
         psnsun_to_cpool(p) + &
         psnshade_to_cpool(p)


      ! maintenance respiration (MR)
      if ( isotope .eq. 'c13' .or. isotope .eq. 'c14') then
         leaf_mr(p)      = leaf_curmr(p)      + leaf_xsmr(p)
         froot_mr(p)     = froot_curmr(p)     + froot_xsmr(p)
         livestem_mr(p)  = livestem_curmr(p)  + livestem_xsmr(p)
         livecroot_mr(p) = livecroot_curmr(p) + livecroot_xsmr(p)
      endif
      mr(p)  = &
         leaf_mr(p)     + &
         froot_mr(p)    + &
         livestem_mr(p) + &
         livecroot_mr(p)

      ! growth respiration (GR)
      ! current GR is respired this time step for new growth displayed in this timestep
      current_gr(p) = &
         cpool_leaf_gr(p)      + &
         cpool_froot_gr(p)     + &
         cpool_livestem_gr(p)  + &
         cpool_deadstem_gr(p)  + &
         cpool_livecroot_gr(p) + &
         cpool_deadcroot_gr(p)

      ! transfer GR is respired this time step for transfer growth displayed in this timestep
      transfer_gr(p) = &
         transfer_leaf_gr(p)      + &
         transfer_froot_gr(p)     + &
         transfer_livestem_gr(p)  + &
         transfer_deadstem_gr(p)  + &
         transfer_livecroot_gr(p) + &
         transfer_deadcroot_gr(p)

      ! storage GR is respired this time step for growth sent to storage for later display
      storage_gr(p) = &
         cpool_leaf_storage_gr(p)      + &
         cpool_froot_storage_gr(p)     + &
         cpool_livestem_storage_gr(p)  + &
         cpool_deadstem_storage_gr(p)  + &
         cpool_livecroot_storage_gr(p) + &
         cpool_deadcroot_storage_gr(p)

      if ( crop_prog .and. ivt(p) >= npcropmin )then
         mr(p) = mr(p) + &
            grain_mr(p)
         current_gr(p) = current_gr(p) + &
            cpool_grain_gr(p)
         transfer_gr(p) = transfer_gr(p) + &
            transfer_grain_gr(p)
         storage_gr(p) = storage_gr(p) + &
            cpool_grain_storage_gr(p)
      end if

      ! GR is the sum of current + transfer + storage GR
      gr(p) = &
         current_gr(p)  + &
         transfer_gr(p) + &
         storage_gr(p)

      ! autotrophic respiration (AR)
      if ( crop_prog .and. ivt(p) >= npcropmin )then
         ar(p) = mr(p) + gr(p) + xsmrpool_to_atm(p) ! xsmr... is -ve (slevis)
      else
         ar(p) = mr(p) + gr(p)
      end if

      ! root respiration (RR)
      rr(p) = &
         froot_mr(p) + &
         cpool_froot_gr(p) + &
         cpool_livecroot_gr(p) + &
         cpool_deadcroot_gr(p) + &
         transfer_froot_gr(p) + &
         transfer_livecroot_gr(p) + &
         transfer_deadcroot_gr(p) + &
         cpool_froot_storage_gr(p) + &
         cpool_livecroot_storage_gr(p) + &
         cpool_deadcroot_storage_gr(p)

      ! net primary production (NPP)
      npp(p) = gpp(p) - ar(p)

      ! update the annual NPP accumulator, for use in allocation code
      if (isotope == 'bulk') then      
         tempsum_npp(p) = tempsum_npp(p) + npp(p)
      end if

      ! aboveground NPP: leaf, live stem, dead stem (AGNPP)
      ! This is supposed to correspond as closely as possible to
      ! field measurements of AGNPP, so it ignores the storage pools
      ! and only treats the fluxes into displayed pools.
      agnpp(p) = &
         cpool_to_leafc(p)                  + &
         leafc_xfer_to_leafc(p)             + &
         cpool_to_livestemc(p)              + &
         livestemc_xfer_to_livestemc(p)     + &
         cpool_to_deadstemc(p)              + &
         deadstemc_xfer_to_deadstemc(p)

     ! belowground NPP: fine root, live coarse root, dead coarse root (BGNPP)
      ! This is supposed to correspond as closely as possible to
      ! field measurements of BGNPP, so it ignores the storage pools
      ! and only treats the fluxes into displayed pools.
      bgnpp(p) = &
         cpool_to_frootc(p)                   + &
         frootc_xfer_to_frootc(p)             + &
         cpool_to_livecrootc(p)               + &
         livecrootc_xfer_to_livecrootc(p)     + &
         cpool_to_deadcrootc(p)               + &
         deadcrootc_xfer_to_deadcrootc(p)

      ! litterfall (LITFALL)
      litfall(p) = &
         leafc_to_litter(p)                 + &
         frootc_to_litter(p)                + &
         m_leafc_to_litter(p)               + &
         m_leafc_storage_to_litter(p)       + &
         m_leafc_xfer_to_litter(p)          + &
         m_frootc_to_litter(p)              + &
         m_frootc_storage_to_litter(p)      + &
         m_frootc_xfer_to_litter(p)         + &
         m_livestemc_to_litter(p)           + &
         m_livestemc_storage_to_litter(p)   + &
         m_livestemc_xfer_to_litter(p)      + &
         m_deadstemc_to_litter(p)           + &
         m_deadstemc_storage_to_litter(p)   + &
         m_deadstemc_xfer_to_litter(p)      + &
         m_livecrootc_to_litter(p)          + &
         m_livecrootc_storage_to_litter(p)  + &
         m_livecrootc_xfer_to_litter(p)     + &
         m_deadcrootc_to_litter(p)          + &
         m_deadcrootc_storage_to_litter(p)  + &
         m_deadcrootc_xfer_to_litter(p)     + &
         m_gresp_storage_to_litter(p)       + &
         m_gresp_xfer_to_litter(p)          + &
! F. Li and S. Levis
         m_leafc_to_litter_fire(p)              + &
         m_leafc_storage_to_litter_fire(p)      + &
         m_leafc_xfer_to_litter_fire(p)         + &
         m_livestemc_to_litter_fire(p)          + &
         m_livestemc_storage_to_litter_fire(p)  + &
         m_livestemc_xfer_to_litter_fire(p)     + &
         m_deadstemc_to_litter_fire(p)          + &
         m_deadstemc_storage_to_litter_fire(p)  + &
         m_deadstemc_xfer_to_litter_fire(p)     + &
         m_frootc_to_litter_fire(p)             + &
         m_frootc_storage_to_litter_fire(p)     + &
         m_frootc_xfer_to_litter_fire(p)        + &
         m_livecrootc_to_litter_fire(p)         + &
         m_livecrootc_storage_to_litter_fire(p) + &
         m_livecrootc_xfer_to_litter_fire(p)    + &
         m_deadcrootc_to_litter_fire(p)         + &
         m_deadcrootc_storage_to_litter_fire(p) + &
         m_deadcrootc_xfer_to_litter_fire(p)    + &
         m_gresp_storage_to_litter_fire(p)      + &
         m_gresp_xfer_to_litter_fire(p)         + &

         hrv_leafc_to_litter(p)             + &
         hrv_leafc_storage_to_litter(p)     + &
         hrv_leafc_xfer_to_litter(p)        + &
         hrv_frootc_to_litter(p)            + &
         hrv_frootc_storage_to_litter(p)    + &
         hrv_frootc_xfer_to_litter(p)       + &
         hrv_livestemc_to_litter(p)         + &
         hrv_livestemc_storage_to_litter(p) + &
         hrv_livestemc_xfer_to_litter(p)    + &
         hrv_deadstemc_storage_to_litter(p) + &
         hrv_deadstemc_xfer_to_litter(p)    + &
         hrv_livecrootc_to_litter(p)        + &
         hrv_livecrootc_storage_to_litter(p)+ &
         hrv_livecrootc_xfer_to_litter(p)   + &
         hrv_deadcrootc_to_litter(p)        + &
         hrv_deadcrootc_storage_to_litter(p)+ &
         hrv_deadcrootc_xfer_to_litter(p)   + &
         hrv_gresp_storage_to_litter(p)     + &
         hrv_gresp_xfer_to_litter(p)
                 
#if (defined CNDV)
      ! update the annual litfall accumulator, for use in mortality code
      tempsum_litfall(p) = tempsum_litfall(p) + leafc_to_litter(p) + frootc_to_litter(p)
#endif

      ! pft-level fire losses (VEGFIRE)
      vegfire(p) = 0._r8
      
      ! pft-level wood harvest
      wood_harvestc(p) = &
         hrv_deadstemc_to_prod10c(p) + &
         hrv_deadstemc_to_prod100c(p)

      ! pft-level carbon losses to fire changed by F. Li and S. Levis
      pft_fire_closs(p) = &
         m_leafc_to_fire(p)                + &
         m_leafc_storage_to_fire(p)        + &
         m_leafc_xfer_to_fire(p)           + &
         m_frootc_to_fire(p)               + &
         m_frootc_storage_to_fire(p)       + &
         m_frootc_xfer_to_fire(p)          + &
         m_livestemc_to_fire(p)            + &
         m_livestemc_storage_to_fire(p)    + &
         m_livestemc_xfer_to_fire(p)       + &
         m_deadstemc_to_fire(p)            + &
         m_deadstemc_storage_to_fire(p)    + &
         m_deadstemc_xfer_to_fire(p)       + &
         m_livecrootc_to_fire(p)           + &
         m_livecrootc_storage_to_fire(p)   + &
         m_livecrootc_xfer_to_fire(p)      + &
         m_deadcrootc_to_fire(p)           + &
         m_deadcrootc_storage_to_fire(p)   + &
         m_deadcrootc_xfer_to_fire(p)      + &
         m_gresp_storage_to_fire(p)        + &
         m_gresp_xfer_to_fire(p)

      ! displayed vegetation carbon, excluding storage and cpool (DISPVEGC)
      dispvegc(p) = &
         leafc(p)      + &
         frootc(p)     + &
         livestemc(p)  + &
         deadstemc(p)  + &
         livecrootc(p) + &
         deadcrootc(p)

      ! stored vegetation carbon, excluding cpool (STORVEGC)
      storvegc(p) = &
      	cpool(p)              + &
         leafc_storage(p)      + &
         frootc_storage(p)     + &
         livestemc_storage(p)  + &
         deadstemc_storage(p)  + &
         livecrootc_storage(p) + &
         deadcrootc_storage(p) + &
         leafc_xfer(p)         + &
         frootc_xfer(p)        + &
         livestemc_xfer(p)     + &
         deadstemc_xfer(p)     + &
         livecrootc_xfer(p)    + &
         deadcrootc_xfer(p)    + &
         gresp_storage(p)      + &
         gresp_xfer(p)

      if ( crop_prog .and. ivt(p) >= npcropmin )then
         storvegc(p) = storvegc(p) + &
            grainc_storage(p)  + &
            grainc_xfer(p)
         agnpp(p) = agnpp(p) + &
            cpool_to_grainc(p)           + &
            grainc_xfer_to_grainc(p)
         litfall(p) = litfall(p) + &
            livestemc_to_litter(p)           + &
            grainc_to_food(p)
         dispvegc(p) = dispvegc(p) + &
            grainc(p)
      end if

      ! total vegetation carbon, excluding cpool (TOTVEGC)
      totvegc(p) = dispvegc(p) + storvegc(p)

      ! total pft-level carbon, including xsmrpool, ctrunc
      totpftc(p) = totvegc(p) + xsmrpool(p) + pft_ctrunc(p)
      
      ! new summary variables for CLAMP
      
      ! (FROOTC_ALLOC) - fine root C allocation
      frootc_alloc(p) = &
        frootc_xfer_to_frootc(p)    + &
        cpool_to_frootc(p)     
              
      ! (FROOTC_LOSS) - fine root C loss changed by F. Li and S. Levis
      frootc_loss(p) = &
        m_frootc_to_litter(p)       + &
        m_frootc_to_fire(p)         + &
        m_frootc_to_litter_fire(p)  + &
        hrv_frootc_to_litter(p)     + &
        frootc_to_litter(p)
      
      ! (LEAFC_ALLOC) - leaf C allocation
      leafc_alloc(p) = &
        leafc_xfer_to_leafc(p)    + &
        cpool_to_leafc(p)     

      ! (LEAFC_LOSS) - leaf C loss changed by F. Li and S. Levis
      leafc_loss(p) = &
        m_leafc_to_litter(p)      + &
        m_leafc_to_fire(p)        + &
        m_leafc_to_litter_fire(p) + &
        hrv_leafc_to_litter(p)    + &
        leafc_to_litter(p)
      
      ! (WOODC) - wood C
      woodc(p) = &
        deadstemc(p)    + &
        livestemc(p)    + &
        deadcrootc(p)   + &
        livecrootc(p)
      
      ! (WOODC_ALLOC) - wood C allocation
      woodc_alloc(p) = &
        livestemc_xfer_to_livestemc(p)  + &
        deadstemc_xfer_to_deadstemc(p)  + &
        livecrootc_xfer_to_livecrootc(p)    + &
        deadcrootc_xfer_to_deadcrootc(p)    + &
        cpool_to_livestemc(p)   + &
        cpool_to_deadstemc(p)   + &
        cpool_to_livecrootc(p)  + &
        cpool_to_deadcrootc(p)
      
      ! (WOODC_LOSS) - wood C loss
      woodc_loss(p) = &
        m_livestemc_to_litter(p)    + &
        m_deadstemc_to_litter(p)    + &
        m_livecrootc_to_litter(p)   + &
        m_deadcrootc_to_litter(p)   + &
        m_livestemc_to_fire(p)      + &
        m_deadstemc_to_fire(p)      + &
        m_livecrootc_to_fire(p)     + &
        m_deadcrootc_to_fire(p)     + &
        hrv_livestemc_to_litter(p)  + &
        hrv_livestemc_storage_to_litter(p) + &
        hrv_livestemc_xfer_to_litter(p)    + &
        hrv_deadstemc_to_prod10c(p)        + &
        hrv_deadstemc_to_prod100c(p)       + &
        hrv_deadstemc_storage_to_litter(p) + &
        hrv_deadstemc_xfer_to_litter(p)    + &
        hrv_livecrootc_to_litter(p)        + &
        hrv_livecrootc_storage_to_litter(p)+ &
        hrv_livecrootc_xfer_to_litter(p)   + &
        hrv_deadcrootc_to_litter(p)        + &
        hrv_deadcrootc_storage_to_litter(p)+ &
        hrv_deadcrootc_xfer_to_litter(p)   

   end do  ! end of pfts loop

   ! use p2c routine to get selected column-average pft-level fluxes and states
   call p2c(num_soilc, filter_soilc, gpp, col_gpp)
   call p2c(num_soilc, filter_soilc, ar, col_ar)
   call p2c(num_soilc, filter_soilc, rr, col_rr)
   call p2c(num_soilc, filter_soilc, npp, col_npp)
   call p2c(num_soilc, filter_soilc, vegfire, col_vegfire)
   call p2c(num_soilc, filter_soilc, wood_harvestc, col_wood_harvestc)
   call p2c(num_soilc, filter_soilc, totvegc, col_totvegc)
   call p2c(num_soilc, filter_soilc, totpftc, col_totpftc)
   call p2c(num_soilc, filter_soilc, pft_fire_closs, col_pft_fire_closs)
   call p2c(num_soilc, filter_soilc, litfall, col_litfall)
   call p2c(num_soilc, filter_soilc, hrv_xsmrpool_to_atm, col_hrv_xsmrpool_to_atm)

   if ( isotope .eq. 'bulk') then
      if (nfix_timeconst .gt. 0._r8 .and. nfix_timeconst .lt. 500._r8 ) then
         ! this code is to calculate an exponentially-relaxed npp value for use in NDynamics code
         col_lag_npp                    => clm3%g%l%c%cps%col_lag_npp
         dtime = get_step_size()
         nfixlags = nfix_timeconst * secspday
         ! column loop
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            if ( col_lag_npp(c) .ne. spval ) then
               col_lag_npp(c) = col_lag_npp(c) * exp(-dtime/nfixlags) &
                    + col_npp(c) * (1._r8 - exp(-dtime/nfixlags))
            else
               ! first timestep
               col_lag_npp(c) = col_npp(c)
            endif
         end do
      endif
   endif

   ! column loop
   do fc = 1,num_soilc
      c = filter_soilc(fc)

      ! some zeroing
      lithr(c) = 0._r8
      somhr(c) = 0._r8
      totlitc(c) = 0._r8
      totsomc(c) = 0._r8
      cwdc(c) = 0._r8
      col_ctrunc(c) = 0._r8
      cwdc_loss(c) = 0._r8
      som_c_leached(c) = 0._r8
      do l = 1, ndecomp_pools
         decomp_cpools(c,l) = 0._r8
      end do
      totlitc_1m(c) = 0._r8
      totsomc_1m(c) = 0._r8
   end do

   ! vertically integrate HR and decomposition cascade fluxes
   do k = 1, ndecomp_cascade_transitions
      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            decomp_cascade_hr(c,k) = decomp_cascade_hr(c,k) + decomp_cascade_hr_vr(c,j,k)*dzsoi_decomp(j) 
            decomp_cascade_ctransfer(c,k) = decomp_cascade_ctransfer(c,k) + decomp_cascade_ctransfer_vr(c,j,k) * dzsoi_decomp(j) 
         end do
      end do
   end do

   ! litter heterotrophic respiration (LITHR)
   do k = 1, ndecomp_cascade_transitions
      if ( is_litter(cascade_donor_pool(k)) ) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            lithr(c) = lithr(c) + decomp_cascade_hr(c,k)
         end do
      end if
   end do

   ! soil organic matter heterotrophic respiration (SOMHR)
   do k = 1, ndecomp_cascade_transitions
      if ( is_soil(cascade_donor_pool(k)) ) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            somhr(c) = somhr(c) + decomp_cascade_hr(c,k)
         end do
      end if
   end do
   
   ! total heterotrophic respiration (HR)
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      hr(c) = lithr(c) + somhr(c)
   end do
   
   ! total heterotrophic respiration, vertically resolved (HR)
   do j = 1,nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         hr_vr(c,j) = 0._r8
      end do
   end do
   do k = 1, ndecomp_cascade_transitions
      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            hr_vr(c,j) = hr_vr(c,j) + decomp_cascade_hr_vr(c,j,k)
         end do
      end do
   end do

   do fc = 1,num_soilc
      c = filter_soilc(fc)
      ! total soil respiration, heterotrophic + root respiration (SR)
      sr(c) = col_rr(c) + hr(c)

      ! total ecosystem respiration, autotrophic + heterotrophic (ER)
      er(c) = col_ar(c) + hr(c)

      ! litter fire losses (LITFIRE)
      litfire(c) = 0._r8
      
      ! total wood product loss
      product_closs(c) = &
         prod10c_loss(c) + &
         prod100c_loss(c) 

      ! soil organic matter fire losses (SOMFIRE)
      somfire(c) = 0._r8

      ! total ecosystem fire losses (TOTFIRE)
      totfire(c) = &
         litfire(c) + &
         somfire(c) + &
         col_vegfire(c)
   end do

   ! vertically integrate column-level carbon fire losses
   do l = 1, ndecomp_pools
      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            m_decomp_cpools_to_fire(c,l) = m_decomp_cpools_to_fire(c,l) + &
                 m_decomp_cpools_to_fire_vr(c,j,l)*dzsoi_decomp(j)
         end do
      end do
   end do

   ! column-level carbon losses to fire, including pft losses
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      col_fire_closs(c) = col_pft_fire_closs(c)
      do l = 1, ndecomp_pools
            col_fire_closs(c) = col_fire_closs(c) + m_decomp_cpools_to_fire(c,l)
      end do
         
      ! column-level carbon losses due to landcover change
      dwt_closs(c) = &
         dwt_conv_cflux(c)

      ! net ecosystem production, excludes fire flux, landcover change, and loss from wood products, positive for sink (NEP)
      nep(c) = col_gpp(c) - er(c)

      ! net biome production of carbon, includes depletion from: fire flux, landcover change flux, and loss
      ! from wood products pools, positive for sink (NBP)
      nbp(c) = nep(c) - col_fire_closs(c) - dwt_closs(c) - product_closs(c)

      ! net ecosystem exchange of carbon, includes fire flux, landcover change flux, loss
      ! from wood products pools, and hrv_xsmrpool flux, positive for source (NEE)
      nee(c) = -nep(c) + col_fire_closs(c) + dwt_closs(c) + product_closs(c) + col_hrv_xsmrpool_to_atm(c)
      ! land use flux and land uptake
      landuseflux(c) = dwt_closs(c) + product_closs(c)
      landuptake(c) = nee(c) - landuseflux(c)
   end do

   ! vertically integrate each of the decomposing C pools
   do l = 1, ndecomp_pools
      do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            decomp_cpools(c,l) = decomp_cpools(c,l) + &
                 decomp_cpools_vr(c,j,l) * dzsoi_decomp(j)
         end do
      end do
   end do

   ! for vertically-resolved soil biogeochemistry, calculate some diagnostics of carbon pools to a given depth
   if ( nlevdecomp .gt. 1) then
      
      ! zero some pools
      do l = 1, ndecomp_pools
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            decomp_cpools_1m(c,l) = 0._r8
         end do
      end do

      ! vertically integrate each of the decomposing C pools to 1 meter
      maxdepth = 1._r8
      do l = 1, ndecomp_pools
         do j = 1, nlevdecomp
            if ( zisoi(j) .le. maxdepth ) then
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  decomp_cpools_1m(c,l) = decomp_cpools_1m(c,l) + &
                       decomp_cpools_vr(c,j,l) * dzsoi_decomp(j)
               end do
            elseif ( zisoi(j-1) .lt. maxdepth ) then
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  decomp_cpools_1m(c,l) = decomp_cpools_1m(c,l) + &
                       decomp_cpools_vr(c,j,l) * (maxdepth - zisoi(j-1))
               end do
            endif
         end do
      end do
      
      ! total litter carbon in the top meter (TOTLITC_1m)
      do l = 1, ndecomp_pools
         if ( is_litter(l) ) then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               totlitc_1m(c) = totlitc_1m(c) + &
                    decomp_cpools_1m(c,l)
            end do
         endif
      end do
      
      ! total soil organic matter carbon in the top meter (TOTSOMC_1m)
      do l = 1, ndecomp_pools
         if ( is_soil(l) ) then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               totsomc_1m(c) = totsomc_1m(c) + &
                    decomp_cpools_1m(c,l)
            end do
         end if
      end do
      
   endif
   
   ! total litter carbon (TOTLITC)
   do l = 1, ndecomp_pools
      if ( is_litter(l) ) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            totlitc(c) = totlitc(c) + &
                 decomp_cpools(c,l)
         end do
      endif
   end do

   ! total soil organic matter carbon (TOTSOMC)
   do l = 1, ndecomp_pools
      if ( is_soil(l) ) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            totsomc(c) = totsomc(c) + &
                 decomp_cpools(c,l)
         end do
      end if
   end do

   ! coarse woody debris carbon
   do l = 1, ndecomp_pools
      if ( is_cwd(l) ) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            cwdc(c) = cwdc(c) + &
                 decomp_cpools(c,l)
         end do
      end if
   end do

   ! truncation carbon
   do j = 1, nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         col_ctrunc(c) = col_ctrunc(c) + &
              col_ctrunc_vr(c,j) * dzsoi_decomp(j)
      end do
   end do
      
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      ! total wood product carbon
      totprodc(c) = &
           prod10c(c) + &
           prod100c(c)	 
      
      ! total ecosystem carbon, including veg but excluding cpool (TOTECOSYSC)
      totecosysc(c) = &
           cwdc(c) + &
           totlitc(c) + &
           totsomc(c) + &
           totprodc(c) + &
           col_totvegc(c)
      
      ! total column carbon, including veg and cpool (TOTCOLC)
      ! adding col_ctrunc, seedc
      totcolc(c) = &
         col_totpftc(c) + &
         cwdc(c) + &
         totlitc(c) + &
         totsomc(c) + &
         totprodc(c) + &
         seedc(c) + &
         col_ctrunc(c)
      
      ! new summary variables for CLAMP
      
      ! (CWDC_HR) - coarse woody debris heterotrophic respiration
      cwdc_hr(c) = 0._r8
   end do

   ! (CWDC_LOSS) - coarse woody debris C loss
   do l = 1, ndecomp_pools
      if ( is_cwd(l) ) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            cwdc_loss(c) = cwdc_loss(c) + m_decomp_cpools_to_fire(c,l)
         end do
      end if
   end do

   do k = 1, ndecomp_cascade_transitions
      if ( is_cwd(cascade_donor_pool(k)) ) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            cwdc_loss(c) = cwdc_loss(c) + decomp_cascade_ctransfer(c,k)
         end do
      end if
   end do
   
   ! (LITTERC_LOSS) - litter C loss      
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      litterc_loss(c) = lithr(c)  
   end do
   do l = 1, ndecomp_pools
      if ( is_litter(l) ) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            litterc_loss(c) = litterc_loss(c) + m_decomp_cpools_to_fire(c,l)
         end do
      end if
   end do
   do k = 1, ndecomp_cascade_transitions
      if ( is_litter(cascade_donor_pool(k)) ) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            litterc_loss(c) = litterc_loss(c) + decomp_cascade_ctransfer(c,k)
         end do
      end if
   end do


   ! add up all vertical transport tendency terms and calculate total som leaching loss as the sum of these
   do l = 1, ndecomp_pools
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         decomp_cpools_leached(c,l) = 0._r8
      end do
      do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            decomp_cpools_leached(c,l) = decomp_cpools_leached(c,l) + decomp_cpools_transport_tendency(c,j,l) * dzsoi_decomp(j)
         end do
      end do
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         som_c_leached(c) = som_c_leached(c) + decomp_cpools_leached(c,l)
      end do
   end do



end subroutine CSummary
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: NSummary
!
! !INTERFACE:
subroutine NSummary(num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! On the radiation time step, perform pft and column-level nitrogen
! summary calculations
!
! !USES:
   use clmtype
   use pft2colMod, only: p2c
   use clm_varpar  , only: nlevdecomp,ndecomp_cascade_transitions,ndecomp_pools
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNEcosystemDyn
!
! !REVISION HISTORY:
! 6/28/04: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
   integer , pointer :: ivt(:)            ! pft vegetation type
   real(r8), pointer :: col_fire_nloss(:) ! (gN/m2/s) total column-level fire N loss
   real(r8), pointer :: col_wood_harvestn(:)
   real(r8), pointer :: denit(:)
   real(r8), pointer :: col_pft_fire_nloss(:) ! (gN/m2/s) total pft-level fire C loss 
   real(r8), pointer :: col_totpftn(:)        ! (gN/m2) total pft-level nitrogen
   real(r8), pointer :: col_totvegn(:)            ! (gN/m2) total vegetation nitrogen
   real(r8), pointer :: cwdn(:)               ! (gN/m2) coarse woody debris N
   real(r8), pointer :: col_ntrunc(:)         ! (gN/m2) column-level sink for N truncation
   real(r8), pointer :: sminn(:)              ! (gN/m2) soil mineral N
   real(r8), pointer :: m_decomp_npools_to_fire_vr(:,:,:)              
   real(r8), pointer :: m_decomp_npools_to_fire(:,:)              
   logical, pointer :: is_litter(:)                       ! TRUE => pool is a litter pool
   logical, pointer :: is_soil(:)                         ! TRUE => pool is a soil pool
   logical, pointer :: is_cwd(:)                          ! TRUE => pool is a cwd pool
#ifndef NITRIF_DENITRIF
   real(r8), pointer :: sminn_to_denit_excess_vr(:,:)
   real(r8), pointer :: sminn_to_denit_excess(:)
   real(r8), pointer :: sminn_leached_vr(:,:)
   real(r8), pointer :: sminn_leached(:)
   real(r8), pointer :: sminn_to_denit_decomp_cascade_vr(:,:,:)   ! vertically-resolved denitrification along decomp cascade (gN/m3/s) 
   real(r8), pointer :: sminn_to_denit_decomp_cascade(:,:)   ! vertically-integrated denitrification along decomp cascade (gN/m2/s) 
#else
   real(r8), pointer :: smin_no3(:)
   real(r8), pointer :: smin_nh4(:)
   real(r8), pointer :: smin_no3_vr(:,:)
   real(r8), pointer :: smin_nh4_vr(:,:)
   real(r8), pointer :: f_nit_vr(:,:)
   real(r8), pointer :: f_nit(:)
   real(r8), pointer :: f_denit_vr(:,:)
   real(r8), pointer :: f_denit(:)
   real(r8), pointer :: pot_f_nit_vr(:,:)
   real(r8), pointer :: pot_f_nit(:)
   real(r8), pointer :: pot_f_denit_vr(:,:)
   real(r8), pointer :: pot_f_denit(:)
   real(r8), pointer :: f_n2o_denit_vr(:,:)    ! flux of N2o from denitrification [gN/m3/s]
   real(r8), pointer :: f_n2o_denit(:)         ! flux of N2o from denitrification [gN/m2/s]
   real(r8), pointer :: f_n2o_nit_vr(:,:)      ! flux of N2o from nitrification [gN/m3/s]
   real(r8), pointer :: f_n2o_nit(:)           ! flux of N2o from nitrification [gN/m2/s]
   real(r8), pointer :: smin_no3_leached_vr(:,:)
   real(r8), pointer :: smin_no3_leached(:)
   real(r8), pointer :: smin_no3_runoff_vr(:,:)
   real(r8), pointer :: smin_no3_runoff(:)
#endif
   real(r8), pointer :: decomp_npools(:,:)         ! (gN/m2)  decomposing (litter, cwd, soil) N pools
   real(r8), pointer :: decomp_npools_vr(:,:,:)    ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
   real(r8), pointer :: decomp_npools_1m(:,:)           ! (gN/m2)  diagnostic: decomposing (litter, cwd, soil) N pools to 1 meter
   integer, pointer :: altmax_indx(:)                  ! maximum annual depth of thaw
   integer, pointer :: altmax_lastyear_indx(:)         ! prior year maximum annual depth of thaw
   real(r8), pointer :: sminn_vr(:,:)              ! (gN/m3) soil mineral N
   real(r8), pointer :: col_ntrunc_vr(:,:)         ! (gN/m3) column-level sink for N truncation
   real(r8), pointer :: supplement_to_sminn(:)
   real(r8), pointer :: supplement_to_sminn_vr(:,:)
   real(r8), pointer :: totcoln(:)            ! (gN/m2) total column nitrogen, incl veg
   real(r8), pointer :: totecosysn(:)         ! (gN/m2) total ecosystem nitrogen, incl veg 
   real(r8), pointer :: totlitn(:)            ! (gN/m2) total litter nitrogen
   real(r8), pointer :: totlitn_1m(:)         ! (gN/m2) total litter nitrogen to 1 meter
   real(r8), pointer :: totsomn(:)            ! (gN/m2) total soil organic matter nitrogen

! fire related variables changed by F. Li and S. Levis
   real(r8), pointer :: m_leafn_to_fire(:)   
   real(r8), pointer :: m_leafn_storage_to_fire(:)                
   real(r8), pointer :: m_leafn_xfer_to_fire(:)  
   real(r8), pointer :: m_livestemn_to_fire(:)    
   real(r8), pointer :: m_livestemn_storage_to_fire(:)        
   real(r8), pointer :: m_livestemn_xfer_to_fire(:) 
   real(r8), pointer :: m_deadstemn_to_fire(:) 
   real(r8), pointer :: m_deadstemn_storage_to_fire(:)           
   real(r8), pointer :: totsomn_1m(:)         ! (gN/m2) total soil organic matter nitrogen to 1 meter
   real(r8), pointer :: m_deadcrootn_storage_to_fire(:) 
   real(r8), pointer :: m_deadcrootn_to_fire(:)         
   real(r8), pointer :: m_deadcrootn_xfer_to_fire(:)
   real(r8), pointer :: m_deadstemn_xfer_to_fire(:) 
   real(r8), pointer :: m_frootn_to_fire(:)        
   real(r8), pointer :: m_frootn_storage_to_fire(:)  
   real(r8), pointer :: m_frootn_xfer_to_fire(:)
   real(r8), pointer :: m_livecrootn_to_fire(:)    
   real(r8), pointer :: m_livecrootn_storage_to_fire(:)      
   real(r8), pointer :: m_livecrootn_xfer_to_fire(:)
   real(r8), pointer :: m_retransn_to_fire(:)   
    
   real(r8), pointer :: hrv_deadstemn_to_prod10n(:)        
   real(r8), pointer :: hrv_deadstemn_to_prod100n(:)       
   real(r8), pointer :: ndeploy(:)
   real(r8), pointer :: pft_fire_nloss(:) ! (gN/m2/s) total pft-level fire C loss 
   real(r8), pointer :: retransn_to_npool(:)          
   real(r8), pointer :: sminn_to_npool(:)             
   real(r8), pointer :: deadcrootn(:)         ! (gN/m2) dead coarse root N
   real(r8), pointer :: deadcrootn_storage(:) ! (gN/m2) dead coarse root N storage
   real(r8), pointer :: deadcrootn_xfer(:)    ! (gN/m2) dead coarse root N transfer
   real(r8), pointer :: deadstemn(:)          ! (gN/m2) dead stem N
   real(r8), pointer :: deadstemn_storage(:)  ! (gN/m2) dead stem N storage
   real(r8), pointer :: deadstemn_xfer(:)     ! (gN/m2) dead stem N transfer
   real(r8), pointer :: dispvegn(:)           ! (gN/m2) displayed veg nitrogen, excluding storage
   real(r8), pointer :: frootn(:)             ! (gN/m2) fine root N
   real(r8), pointer :: frootn_storage(:)     ! (gN/m2) fine root N storage
   real(r8), pointer :: frootn_xfer(:)        ! (gN/m2) fine root N transfer
   real(r8), pointer :: leafn(:)              ! (gN/m2) leaf N 
   real(r8), pointer :: leafn_storage(:)      ! (gN/m2) leaf N storage
   real(r8), pointer :: leafn_xfer(:)         ! (gN/m2) leaf N transfer
   real(r8), pointer :: livecrootn(:)         ! (gN/m2) live coarse root N
   real(r8), pointer :: livecrootn_storage(:) ! (gN/m2) live coarse root N storage
   real(r8), pointer :: livecrootn_xfer(:)    ! (gN/m2) live coarse root N transfer
   real(r8), pointer :: grainn(:)             ! (gN/m2) grain N
   real(r8), pointer :: grainn_storage(:)     ! (gN/m2) grain N storage
   real(r8), pointer :: grainn_xfer(:)        ! (gN/m2) grain N transfer
   real(r8), pointer :: livestemn(:)          ! (gN/m2) live stem N
   real(r8), pointer :: livestemn_storage(:)  ! (gN/m2) live stem N storage
   real(r8), pointer :: livestemn_xfer(:)     ! (gN/m2) live stem N transfer
   real(r8), pointer :: retransn(:)           ! (gN/m2) plant pool of retranslocated N
   real(r8), pointer :: npool(:)              ! (gN/m2) temporary plant N pool
   real(r8), pointer :: pft_ntrunc(:)         ! (gN/m2) pft-level sink for N truncation
   real(r8), pointer :: storvegn(:)           ! (gN/m2) stored vegetation nitrogen
   real(r8), pointer :: totpftn(:)            ! (gN/m2) total pft-level nitrogen
   real(r8), pointer :: totvegn(:)            ! (gN/m2) total vegetation nitrogen
   ! for landcover change
   real(r8), pointer :: wood_harvestn(:)                    ! total N losses to wood product pools (gN/m2/s)
   real(r8), pointer :: dwt_nloss(:)          ! (gN/m2/s) total nitrogen loss from product pools and conversion
   real(r8), pointer :: dwt_conv_nflux(:)     ! (gN/m2/s) conversion N flux (immediate loss to atm)
   real(r8), pointer :: seedn(:)              ! (gN/m2) column-level pool for seeding new PFTs
   real(r8), pointer :: prod10n_loss(:)       ! (gN/m2/s) loss from 10-yr wood product pool
   real(r8), pointer :: prod100n_loss(:)      ! (gN/m2/s) loss from 100-yr wood product pool
   real(r8), pointer :: product_nloss(:)      ! (gN/m2/s) total wood product nitrogen loss
   real(r8), pointer :: prod10n(:)            ! (gN/m2) wood product N pool, 10-year lifespan
   real(r8), pointer :: prod100n(:)           ! (gN/m2) wood product N pool, 100-year lifespan
   real(r8), pointer :: totprodn(:)           ! (gN/m2) total wood product N

   real(r8), pointer :: decomp_cascade_ntransfer_vr(:,:,:)
   real(r8), pointer :: decomp_cascade_ntransfer(:,:)
   real(r8), pointer :: decomp_cascade_sminn_flux_vr(:,:,:)   ! vert-res mineral N flux for transition along decomposition cascade (gN/m3/s)
   real(r8), pointer :: decomp_cascade_sminn_flux(:,:)        ! vert-int (diagnostic) mineral N flux for transition along decomposition cascade (gN/m2/s)

   real(r8), pointer :: som_n_leached(:)                           ! total SOM N loss from vertical transport (gN/m^2/s)
   real(r8), pointer :: decomp_npools_leached(:,:)                 ! N loss from vertical transport from each decomposing N pool (gN/m^2/s)
   real(r8), pointer :: decomp_npools_transport_tendency(:,:,:)    ! N tendency due to vertical transport in decomposing N pools (gN/m^3/s)

!
! local pointers to implicit in/out scalars
!
! local pointers to implicit out scalars
!
! !OTHER LOCAL VARIABLES:
   integer :: c,p,j,k,l       ! indices
   integer :: fp,fc       ! lake filter indices
   real(r8) :: maxdepth    ! depth to integrate soil variables

!EOP
!-----------------------------------------------------------------------
    ! assign local pointers
    ivt                            => clm3%g%l%c%p%itype
    col_fire_nloss                 => clm3%g%l%c%cnf%col_fire_nloss
    denit                          => clm3%g%l%c%cnf%denit
    col_pft_fire_nloss             => clm3%g%l%c%cnf%pnf_a%pft_fire_nloss
    cwdn                           => clm3%g%l%c%cns%cwdn
    col_ntrunc                     => clm3%g%l%c%cns%col_ntrunc
    sminn                          => clm3%g%l%c%cns%sminn
    m_decomp_npools_to_fire_vr     => clm3%g%l%c%cnf%m_decomp_npools_to_fire_vr
    m_decomp_npools_to_fire        => clm3%g%l%c%cnf%m_decomp_npools_to_fire
    is_litter                               => decomp_cascade_con%is_litter
    is_soil                                 => decomp_cascade_con%is_soil
    is_cwd                                  => decomp_cascade_con%is_cwd
#ifndef NITRIF_DENITRIF
    sminn_to_denit_excess_vr          => clm3%g%l%c%cnf%sminn_to_denit_excess_vr
    sminn_to_denit_excess             => clm3%g%l%c%cnf%sminn_to_denit_excess
    sminn_to_denit_decomp_cascade_vr  => clm3%g%l%c%cnf%sminn_to_denit_decomp_cascade_vr
    sminn_to_denit_decomp_cascade     => clm3%g%l%c%cnf%sminn_to_denit_decomp_cascade
    sminn_leached_vr                  => clm3%g%l%c%cnf%sminn_leached_vr
    sminn_leached                     => clm3%g%l%c%cnf%sminn_leached
#else
    smin_no3                          => clm3%g%l%c%cns%smin_no3
    smin_nh4                          => clm3%g%l%c%cns%smin_nh4
    smin_no3_vr                       => clm3%g%l%c%cns%smin_no3_vr
    smin_nh4_vr                       => clm3%g%l%c%cns%smin_nh4_vr
    f_nit_vr                          => clm3%g%l%c%cnf%f_nit_vr
    f_nit                             => clm3%g%l%c%cnf%f_nit
    f_denit_vr                        => clm3%g%l%c%cnf%f_denit_vr
    f_denit                           => clm3%g%l%c%cnf%f_denit
    pot_f_nit_vr                      => clm3%g%l%c%cnf%pot_f_nit_vr
    pot_f_nit                         => clm3%g%l%c%cnf%pot_f_nit
    pot_f_denit_vr                    => clm3%g%l%c%cnf%pot_f_denit_vr
    pot_f_denit                       => clm3%g%l%c%cnf%pot_f_denit
    f_n2o_denit_vr                    => clm3%g%l%c%cnf%f_n2o_denit_vr
    f_n2o_nit_vr                      => clm3%g%l%c%cnf%f_n2o_nit_vr
    f_n2o_denit                       => clm3%g%l%c%cnf%f_n2o_denit
    f_n2o_nit                         => clm3%g%l%c%cnf%f_n2o_nit
    smin_no3_leached_vr               => clm3%g%l%c%cnf%smin_no3_leached_vr
    smin_no3_leached                  => clm3%g%l%c%cnf%smin_no3_leached
    smin_no3_runoff_vr                => clm3%g%l%c%cnf%smin_no3_runoff_vr
    smin_no3_runoff                   => clm3%g%l%c%cnf%smin_no3_runoff
#endif
    decomp_npools                     => clm3%g%l%c%cns%decomp_npools
    decomp_npools_vr                  => clm3%g%l%c%cns%decomp_npools_vr
    decomp_npools_1m                  => clm3%g%l%c%cns%decomp_npools_1m
    altmax_indx                       => clm3%g%l%c%cps%altmax_indx
    altmax_lastyear_indx              => clm3%g%l%c%cps%altmax_lastyear_indx
    sminn_vr                          => clm3%g%l%c%cns%sminn_vr
    col_ntrunc_vr                     => clm3%g%l%c%cns%col_ntrunc_vr
    supplement_to_sminn               => clm3%g%l%c%cnf%supplement_to_sminn
    supplement_to_sminn_vr            => clm3%g%l%c%cnf%supplement_to_sminn_vr
    col_totpftn                    => clm3%g%l%c%cns%pns_a%totpftn
    col_totvegn                    => clm3%g%l%c%cns%pns_a%totvegn
    totcoln                        => clm3%g%l%c%cns%totcoln
    totecosysn                     => clm3%g%l%c%cns%totecosysn
    totlitn                        => clm3%g%l%c%cns%totlitn
    totlitn_1m                     => clm3%g%l%c%cns%totlitn_1m
    totsomn                        => clm3%g%l%c%cns%totsomn
    m_leafn_to_fire                => clm3%g%l%c%p%pnf%m_leafn_to_fire
    m_leafn_storage_to_fire        => clm3%g%l%c%p%pnf%m_leafn_storage_to_fire
    m_leafn_xfer_to_fire           => clm3%g%l%c%p%pnf%m_leafn_xfer_to_fire
    m_livestemn_to_fire            => clm3%g%l%c%p%pnf%m_livestemn_to_fire
    m_livestemn_storage_to_fire    => clm3%g%l%c%p%pnf%m_livestemn_storage_to_fire
    m_livestemn_xfer_to_fire       => clm3%g%l%c%p%pnf%m_livestemn_xfer_to_fire
    m_deadstemn_to_fire            => clm3%g%l%c%p%pnf%m_deadstemn_to_fire
    totsomn_1m                     => clm3%g%l%c%cns%totsomn_1m
    m_deadcrootn_storage_to_fire   => clm3%g%l%c%p%pnf%m_deadcrootn_storage_to_fire
    m_deadcrootn_to_fire           => clm3%g%l%c%p%pnf%m_deadcrootn_to_fire
    m_deadcrootn_xfer_to_fire      => clm3%g%l%c%p%pnf%m_deadcrootn_xfer_to_fire
    m_deadstemn_storage_to_fire    => clm3%g%l%c%p%pnf%m_deadstemn_storage_to_fire
    m_deadstemn_xfer_to_fire       => clm3%g%l%c%p%pnf%m_deadstemn_xfer_to_fire
    m_frootn_to_fire               => clm3%g%l%c%p%pnf%m_frootn_to_fire
    m_frootn_storage_to_fire       => clm3%g%l%c%p%pnf%m_frootn_storage_to_fire
    m_frootn_xfer_to_fire          => clm3%g%l%c%p%pnf%m_frootn_xfer_to_fire
    m_livecrootn_to_fire           => clm3%g%l%c%p%pnf%m_livecrootn_to_fire
    m_livecrootn_storage_to_fire   => clm3%g%l%c%p%pnf%m_livecrootn_storage_to_fire
    m_livecrootn_xfer_to_fire      => clm3%g%l%c%p%pnf%m_livecrootn_xfer_to_fire
    m_deadcrootn_to_fire           => clm3%g%l%c%p%pnf%m_deadcrootn_to_fire
    m_deadcrootn_storage_to_fire   => clm3%g%l%c%p%pnf%m_deadcrootn_storage_to_fire
    m_deadcrootn_xfer_to_fire      => clm3%g%l%c%p%pnf%m_deadcrootn_xfer_to_fire
    m_retransn_to_fire             => clm3%g%l%c%p%pnf%m_retransn_to_fire
   

    hrv_deadstemn_to_prod10n         => clm3%g%l%c%p%pnf%hrv_deadstemn_to_prod10n        
    hrv_deadstemn_to_prod100n        => clm3%g%l%c%p%pnf%hrv_deadstemn_to_prod100n       
    ndeploy                        => clm3%g%l%c%p%pnf%ndeploy
    pft_fire_nloss                 => clm3%g%l%c%p%pnf%pft_fire_nloss
    retransn_to_npool              => clm3%g%l%c%p%pnf%retransn_to_npool
    sminn_to_npool                 => clm3%g%l%c%p%pnf%sminn_to_npool
    deadcrootn                     => clm3%g%l%c%p%pns%deadcrootn
    deadcrootn_storage             => clm3%g%l%c%p%pns%deadcrootn_storage
    deadcrootn_xfer                => clm3%g%l%c%p%pns%deadcrootn_xfer
    deadstemn                      => clm3%g%l%c%p%pns%deadstemn
    deadstemn_storage              => clm3%g%l%c%p%pns%deadstemn_storage
    deadstemn_xfer                 => clm3%g%l%c%p%pns%deadstemn_xfer
    dispvegn                       => clm3%g%l%c%p%pns%dispvegn
    frootn                         => clm3%g%l%c%p%pns%frootn
    frootn_storage                 => clm3%g%l%c%p%pns%frootn_storage
    frootn_xfer                    => clm3%g%l%c%p%pns%frootn_xfer
    leafn                          => clm3%g%l%c%p%pns%leafn
    leafn_storage                  => clm3%g%l%c%p%pns%leafn_storage
    leafn_xfer                     => clm3%g%l%c%p%pns%leafn_xfer
    livecrootn                     => clm3%g%l%c%p%pns%livecrootn
    livecrootn_storage             => clm3%g%l%c%p%pns%livecrootn_storage
    livecrootn_xfer                => clm3%g%l%c%p%pns%livecrootn_xfer
    grainn                         => clm3%g%l%c%p%pns%grainn
    grainn_storage                 => clm3%g%l%c%p%pns%grainn_storage
    grainn_xfer                    => clm3%g%l%c%p%pns%grainn_xfer
    livestemn                      => clm3%g%l%c%p%pns%livestemn
    livestemn_storage              => clm3%g%l%c%p%pns%livestemn_storage
    livestemn_xfer                 => clm3%g%l%c%p%pns%livestemn_xfer
    retransn                       => clm3%g%l%c%p%pns%retransn
    npool                          => clm3%g%l%c%p%pns%npool
    pft_ntrunc                     => clm3%g%l%c%p%pns%pft_ntrunc
    storvegn                       => clm3%g%l%c%p%pns%storvegn
    totpftn                        => clm3%g%l%c%p%pns%totpftn
    totvegn                        => clm3%g%l%c%p%pns%totvegn
    ! dynamic landcover pointers
    wood_harvestn                  => clm3%g%l%c%p%pnf%wood_harvestn
    col_wood_harvestn              => clm3%g%l%c%cnf%pnf_a%wood_harvestn 
    dwt_nloss                      => clm3%g%l%c%cnf%dwt_nloss
    dwt_conv_nflux                 => clm3%g%l%c%cnf%dwt_conv_nflux
    prod10n_loss                   => clm3%g%l%c%cnf%prod10n_loss
    prod100n_loss                  => clm3%g%l%c%cnf%prod100n_loss
    product_nloss                  => clm3%g%l%c%cnf%product_nloss
    seedn                          => clm3%g%l%c%cns%seedn
    prod10n                        => clm3%g%l%c%cns%prod10n
    prod100n                       => clm3%g%l%c%cns%prod100n
    totprodn                       => clm3%g%l%c%cns%totprodn
    som_n_leached                      => clm3%g%l%c%cnf%som_n_leached
    decomp_npools_leached              => clm3%g%l%c%cnf%decomp_npools_leached
    decomp_npools_transport_tendency   => clm3%g%l%c%cnf%decomp_npools_transport_tendency
    decomp_cascade_ntransfer_vr       => clm3%g%l%c%cnf%decomp_cascade_ntransfer_vr
    decomp_cascade_ntransfer          => clm3%g%l%c%cnf%decomp_cascade_ntransfer
    decomp_cascade_sminn_flux_vr      => clm3%g%l%c%cnf%decomp_cascade_sminn_flux_vr
    decomp_cascade_sminn_flux         => clm3%g%l%c%cnf%decomp_cascade_sminn_flux

   ! pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)

      ! calculate pft-level summary nitrogen fluxes and states

      ! total N deployment (from sminn and retranslocated N pool) (NDEPLOY)
      ndeploy(p) = &
         sminn_to_npool(p) + &
         retransn_to_npool(p)

      ! pft-level wood harvest
      wood_harvestn(p) = &
         hrv_deadstemn_to_prod10n(p) + &
         hrv_deadstemn_to_prod100n(p)

      ! total pft-level fire N losses
      pft_fire_nloss(p) = &
         m_leafn_to_fire(p)               + &
         m_leafn_storage_to_fire(p)       + &
         m_leafn_xfer_to_fire(p)          + &
         m_frootn_to_fire(p)              + &
         m_frootn_storage_to_fire(p)      + &
         m_frootn_xfer_to_fire(p)         + &
         m_livestemn_to_fire(p)           + &
         m_livestemn_storage_to_fire(p)   + &
         m_livestemn_xfer_to_fire(p)      + &
         m_deadstemn_to_fire(p)           + &
         m_deadstemn_storage_to_fire(p)   + &
         m_deadstemn_xfer_to_fire(p)      + &
         m_livecrootn_to_fire(p)          + &
         m_livecrootn_storage_to_fire(p)  + &
         m_livecrootn_xfer_to_fire(p)     + &
         m_deadcrootn_to_fire(p)          + &
         m_deadcrootn_storage_to_fire(p)  + &
         m_deadcrootn_xfer_to_fire(p)     + &
         m_retransn_to_fire(p)

      ! displayed vegetation nitrogen, excluding storage (DISPVEGN)
      dispvegn(p) = &
         leafn(p)      + &
         frootn(p)     + &
         livestemn(p)  + &
         deadstemn(p)  + &
         livecrootn(p) + &
         deadcrootn(p)

      ! stored vegetation nitrogen, including retranslocated N pool (STORVEGN)
      storvegn(p) = &
         leafn_storage(p)      + &
         frootn_storage(p)     + &
         livestemn_storage(p)  + &
         deadstemn_storage(p)  + &
         livecrootn_storage(p) + &
         deadcrootn_storage(p) + &
         leafn_xfer(p)         + &
         frootn_xfer(p)        + &
         livestemn_xfer(p)     + &
         deadstemn_xfer(p)     + &
         livecrootn_xfer(p)    + &
         deadcrootn_xfer(p)    + &
		 npool(p)              + &
         retransn(p)

      if ( crop_prog .and. ivt(p) >= npcropmin )then
         dispvegn(p) = dispvegn(p) + &
            grainn(p)
          storvegn(p) = storvegn(p) + &
             grainn_storage(p)     + &
             grainn_xfer(p)
      end if

      ! total vegetation nitrogen (TOTVEGN)
      totvegn(p) = dispvegn(p) + storvegn(p)

      ! total pft-level carbon (add pft_ntrunc)
      totpftn(p) = totvegn(p) + pft_ntrunc(p)


   end do  ! end of pfts loop

   ! use p2c routine to get selected column-average pft-level fluxes and states
   call p2c(num_soilc, filter_soilc, pft_fire_nloss, col_pft_fire_nloss)
   call p2c(num_soilc, filter_soilc, wood_harvestn, col_wood_harvestn)
   call p2c(num_soilc, filter_soilc, totvegn, col_totvegn)
   call p2c(num_soilc, filter_soilc, totpftn, col_totpftn)

   ! column loops
   do fc = 1,num_soilc
      c = filter_soilc(fc)

      ! some zeroing
      denit(c) = 0._r8
#ifdef NITRIF_DENITRIF
      smin_no3(c) = 0._r8
      smin_nh4(c) = 0._r8
#endif
      totlitn(c) = 0._r8
      totsomn(c) = 0._r8
      cwdn(c) = 0._r8
      sminn(c) = 0._r8
      col_ntrunc(c) = 0._r8
      supplement_to_sminn(c) = 0._r8
      som_n_leached(c) = 0._r8
      totlitn_1m(c) = 0._r8
      totsomn_1m(c) = 0._r8
   end do

   ! vertically integrate decomposing N cascade fluxes and soil mineral N fluxes associated with decomposition cascade
   do k = 1, ndecomp_cascade_transitions
      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            decomp_cascade_ntransfer(c,k) = decomp_cascade_ntransfer(c,k) + decomp_cascade_ntransfer_vr(c,j,k) * dzsoi_decomp(j) 
            decomp_cascade_sminn_flux(c,k) = decomp_cascade_sminn_flux(c,k) + decomp_cascade_sminn_flux_vr(c,j,k) * dzsoi_decomp(j) 
         end do
      end do
   end do
   
#ifndef NITRIF_DENITRIF
   ! vertically integrate each denitrification flux
   do l = 1, ndecomp_cascade_transitions
      do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            sminn_to_denit_decomp_cascade(c,l) = sminn_to_denit_decomp_cascade(c,l) + &
                 sminn_to_denit_decomp_cascade_vr(c,j,l) * dzsoi_decomp(j)
         end do
      end do
   end do

   ! vertically integrate bulk denitrification and  leaching flux
   do j = 1, nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         sminn_to_denit_excess(c) = sminn_to_denit_excess(c) + sminn_to_denit_excess_vr(c,j) * dzsoi_decomp(j)
         sminn_leached(c) = sminn_leached(c) + sminn_leached_vr(c,j) * dzsoi_decomp(j)
      end do
   end do

   ! total N denitrification (DENIT)
   do l = 1, ndecomp_cascade_transitions
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         denit(c) = denit(c) + sminn_to_denit_decomp_cascade(c,l)
      end do
   end do

   do fc = 1,num_soilc
      c = filter_soilc(fc)
      denit(c) =  denit(c) + sminn_to_denit_excess(c)
   end do

#else

   ! vertically integrate NO3 NH4 N2O fluxes and pools
   do j = 1, nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         ! nitrification and denitrification fluxes
         f_nit(c) = f_nit(c) + f_nit_vr(c,j) * dzsoi_decomp(j)
         f_denit(c) = f_denit(c) + f_denit_vr(c,j) * dzsoi_decomp(j)
         pot_f_nit(c) = pot_f_nit(c) + pot_f_nit_vr(c,j) * dzsoi_decomp(j)
         pot_f_denit(c) = pot_f_denit(c) + pot_f_denit_vr(c,j) * dzsoi_decomp(j)
         f_n2o_nit(c) = f_n2o_nit(c) + f_n2o_nit_vr(c,j) * dzsoi_decomp(j)
         f_n2o_denit(c) = f_n2o_denit(c) + f_n2o_denit_vr(c,j) * dzsoi_decomp(j)
         
         ! leaching/runoff flux
         smin_no3_leached(c) = smin_no3_leached(c) + smin_no3_leached_vr(c,j) * dzsoi_decomp(j)
         smin_no3_runoff(c) = smin_no3_runoff(c) + smin_no3_runoff_vr(c,j) * dzsoi_decomp(j)
         
         ! mineral N pools (must set to zero first since they are state rather than flux variables)
         smin_no3(c) = smin_no3(c) + smin_no3_vr(c,j) * dzsoi_decomp(j)
         smin_nh4(c) = smin_nh4(c) + smin_nh4_vr(c,j) * dzsoi_decomp(j)
      end do
   end do

   do fc = 1,num_soilc
      c = filter_soilc(fc)
      denit(c) = f_denit(c)
   end do

#endif

   ! vertically integrate column-level fire N losses
   do k = 1, ndecomp_pools
      do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            m_decomp_npools_to_fire(c,k) = m_decomp_npools_to_fire(c,k) + &
                 m_decomp_npools_to_fire_vr(c,j,k) * dzsoi_decomp(j)
         end do
      end do
   end do
   
   ! total column-level fire N losses
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      col_fire_nloss(c) = col_pft_fire_nloss(c)
   end do
   do k = 1, ndecomp_pools
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         col_fire_nloss(c) = col_fire_nloss(c) + &
              m_decomp_npools_to_fire(c,k)
      end do
   end do
   
   
   ! vertically integrate each of the decomposing N pools
   do l = 1, ndecomp_pools
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         decomp_npools(c,l) = 0._r8
      end do
      do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            decomp_npools(c,l) = decomp_npools(c,l) + &
                 decomp_npools_vr(c,j,l) * dzsoi_decomp(j)
         end do
      end do
   end do

   ! for vertically-resolved soil biogeochemistry, calculate some diagnostics of carbon pools to a given depth
   if ( nlevdecomp .gt. 1) then

      do l = 1, ndecomp_pools
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            decomp_npools_1m(c,l) = 0._r8
         end do
      end do

      ! vertically integrate each of the decomposing n pools to 1 meter
      maxdepth = 1._r8
      do l = 1, ndecomp_pools
         do j = 1, nlevdecomp
            if ( zisoi(j) .le. maxdepth ) then
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  decomp_npools_1m(c,l) = decomp_npools_1m(c,l) + &
                       decomp_npools_vr(c,j,l) * dzsoi_decomp(j)
               end do
            elseif ( zisoi(j-1) .lt. maxdepth ) then
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  decomp_npools_1m(c,l) = decomp_npools_1m(c,l) + &
                       decomp_npools_vr(c,j,l) * (maxdepth - zisoi(j-1))
               end do
            endif
         end do
      end do
      
      
      ! total litter nitrogen to 1 meter (TOTLITN_1m)
      do l = 1, ndecomp_pools
         if ( is_litter(l) ) then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               totlitn_1m(c) = totlitn_1m(c) + &
                    decomp_npools_1m(c,l)
            end do
         end if
      end do
      
      ! total soil organic matter nitrogen to 1 meter (TOTSOMN_1m)
      do l = 1, ndecomp_pools
         if ( is_soil(l) ) then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               totsomn_1m(c) = totsomn_1m(c) + &
                    decomp_npools_1m(c,l)
            end do
         end if
      end do
      
   endif
   
   ! total litter nitrogen (TOTLITN)
   do l = 1, ndecomp_pools
      if ( is_litter(l) ) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            totlitn(c) = totlitn(c) + &
                 decomp_npools(c,l)
         end do
      end if
   end do
   
   ! total soil organic matter nitrogen (TOTSOMN)
   do l = 1, ndecomp_pools
      if ( is_soil(l) ) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            totsomn(c) = totsomn(c) + &
                 decomp_npools(c,l)
         end do
      end if
   end do
   
   ! total cwdn
   do l = 1, ndecomp_pools
      if ( is_cwd(l) ) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            cwdn(c) = cwdn(c) + &
                 decomp_npools(c,l)
         end do
      end if
   end do
   
   ! total sminn
   do j = 1, nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         sminn(c) = sminn(c) + &
              sminn_vr(c,j) * dzsoi_decomp(j)
      end do
   end do

   ! total col_ntrunc
   do j = 1, nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         col_ntrunc(c) = col_ntrunc(c) + &
              col_ntrunc_vr(c,j) * dzsoi_decomp(j)
      end do
   end do

   ! supplementary N supplement_to_sminn
   do j = 1, nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         supplement_to_sminn(c) = supplement_to_sminn(c) + &
              supplement_to_sminn_vr(c,j) * dzsoi_decomp(j)
      end do
   end do

   do fc = 1,num_soilc
      c = filter_soilc(fc)
      
      ! column-level N losses due to landcover change
      dwt_nloss(c) = &
         dwt_conv_nflux(c)
         
      ! total wood product N loss
      product_nloss(c) = &
         prod10n_loss(c) + &
         prod100n_loss(c) 

      ! total wood product nitrogen
      totprodn(c) = &
         prod10n(c) + &
	     prod100n(c)	 

      ! total ecosystem nitrogen, including veg (TOTECOSYSN)
      totecosysn(c) = &
         cwdn(c) + &
         totlitn(c) + &
         totsomn(c) + &
         sminn(c) + &
	     totprodn(c) + &
         col_totvegn(c)

      ! total column nitrogen, including pft (TOTCOLN)
      totcoln(c) = &
         col_totpftn(c) + &
         cwdn(c) + &
         totlitn(c) + &
         totsomn(c) + &
         sminn(c) + &
	     totprodn(c) + &
		 seedn(c) + &
		 col_ntrunc(c)
   end do
   
   ! add up all vertical transport tendency terms and calculate total som leaching loss as the sum of these
   do l = 1, ndecomp_pools
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         decomp_npools_leached(c,l) = 0._r8
      end do
      do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            decomp_npools_leached(c,l) = decomp_npools_leached(c,l) + decomp_npools_transport_tendency(c,j,l) * dzsoi_decomp(j)
         end do
      end do
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         som_n_leached(c) = som_n_leached(c) + decomp_npools_leached(c,l)
      end do
   end do



end subroutine NSummary
!-----------------------------------------------------------------------

#endif

end module CNSummaryMod
