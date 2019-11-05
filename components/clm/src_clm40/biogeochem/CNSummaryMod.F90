module CNSummaryMod

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
    use pftvarcon   , only: npcropmin, nc3crop
    use surfrdMod   , only: crop_prog
    implicit none
    save
    private
! !PUBLIC MEMBER FUNCTIONS:
    public :: CSummary
    public :: NSummary
!
! !REVISION HISTORY:
! 4/23/2004: Created by Peter Thornton
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CSummary
!
! !INTERFACE:
subroutine CSummary(num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! On the radiation time step, perform pft and column-level carbon
! summary calculations
!
! !USES:
   use clmtype
   use pft2colMod, only: p2c
   use clm_varctl, only: iulog, use_cndv
   use shr_sys_mod, only: shr_sys_flush
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
! 12/9/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
   integer , pointer :: ivt(:)                ! pft vegetation type
   real(r8), pointer :: col_fire_closs(:)     ! (gC/m2/s) total column-level fire C loss
   real(r8), pointer :: er(:)                 ! (gC/m2/s) total ecosystem respiration, autotrophic + heterotrophic
   real(r8), pointer :: hr(:)                 ! (gC/m2/s) total heterotrophic respiration
   real(r8), pointer :: litfire(:)            ! (gC/m2/s) litter fire losses
   real(r8), pointer :: lithr(:)              ! (gC/m2/s) litter heterotrophic respiration 
   real(r8), pointer :: litr1_hr(:)       
   real(r8), pointer :: litr2_hr(:)        
   real(r8), pointer :: litr3_hr(:)        
   real(r8), pointer :: m_cwdc_to_fire(:)
   real(r8), pointer :: m_litr1c_to_fire(:)             
   real(r8), pointer :: m_litr2c_to_fire(:)             
   real(r8), pointer :: m_litr3c_to_fire(:)             
   real(r8), pointer :: nee(:)                ! (gC/m2/s) net ecosystem exchange of carbon, includes fire, land-use, harvest, and hrv_xsmrpool flux, positive for source
   real(r8), pointer :: nep(:)                ! (gC/m2/s) net ecosystem production, excludes fire, land-use, and harvest flux, positive for sink
   real(r8), pointer :: nbp(:)                ! (gC/m2/s) net biome production, includes fire, land-use, and harvest flux, positive for sink
   real(r8), pointer :: col_ar(:)             ! (gC/m2/s) autotrophic respiration (MR + GR)
   real(r8), pointer :: col_gpp(:)            ! GPP flux before downregulation (gC/m2/s)
   real(r8), pointer :: col_npp(:)            ! (gC/m2/s) net primary production
   real(r8), pointer :: col_pft_fire_closs(:) ! (gC/m2/s) total pft-level fire C loss 
   real(r8), pointer :: col_litfall(:)        ! (gC/m2/s) total pft-level litterfall C loss 
   real(r8), pointer :: col_rr(:)             ! (gC/m2/s) root respiration (fine root MR + total root GR)
   real(r8), pointer :: col_vegfire(:)        ! (gC/m2/s) pft-level fire loss (obsolete, mark for removal)
   real(r8), pointer :: col_wood_harvestc(:)
   real(r8), pointer :: soil1_hr(:)        
   real(r8), pointer :: soil2_hr(:)        
   real(r8), pointer :: soil3_hr(:) 
   real(r8), pointer :: soil4_hr(:) 
   real(r8), pointer :: somfire(:)            ! (gC/m2/s) soil organic matter fire losses
   real(r8), pointer :: somhr(:)              ! (gC/m2/s) soil organic matter heterotrophic respiration
   real(r8), pointer :: sr(:)                 ! (gC/m2/s) total soil respiration (HR + root resp)
   real(r8), pointer :: totfire(:)            ! (gC/m2/s) total ecosystem fire losses
   real(r8), pointer :: cwdc(:)               ! (gC/m2) coarse woody debris C
   real(r8), pointer :: litr1c(:)             ! (gC/m2) litter labile C
   real(r8), pointer :: litr2c(:)             ! (gC/m2) litter cellulose C
   real(r8), pointer :: litr3c(:)             ! (gC/m2) litter lignin C
   real(r8), pointer :: col_totpftc(:)        ! (gC/m2) total pft-level carbon, including cpool
   real(r8), pointer :: col_totvegc(:)        ! (gC/m2) total vegetation carbon, excluding cpool
   real(r8), pointer :: soil1c(:)             ! (gC/m2) soil organic matter C (fast pool)
   real(r8), pointer :: soil2c(:)             ! (gC/m2) soil organic matter C (medium pool)
   real(r8), pointer :: soil3c(:)             ! (gC/m2) soil organic matter C (slow pool)
   real(r8), pointer :: soil4c(:)             ! (gC/m2) soil organic matter C (slowest pool)
   real(r8), pointer :: col_ctrunc(:)         ! (gC/m2) column-level sink for C truncation
   real(r8), pointer :: totcolc(:)            ! (gC/m2) total column carbon, incl veg and cpool
   real(r8), pointer :: totecosysc(:)         ! (gC/m2) total ecosystem carbon, incl veg but excl cpool
   real(r8), pointer :: totlitc(:)            ! (gC/m2) total litter carbon
   real(r8), pointer :: totsomc(:)            ! (gC/m2) total soil organic matter carbon
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
   real(r8), pointer :: gpp(:)                !GPP flux before downregulation (gC/m2/s)
   real(r8), pointer :: gr(:)                 ! (gC/m2/s) total growth respiration
   real(r8), pointer :: leafc_to_litter(:)
   real(r8), pointer :: leafc_xfer_to_leafc(:)         
   real(r8), pointer :: leaf_mr(:)
   real(r8), pointer :: litfall(:)            ! (gC/m2/s) litterfall (leaves and fine roots)
   real(r8), pointer :: livecrootc_xfer_to_livecrootc(:)
   real(r8), pointer :: livecroot_mr(:)
   real(r8), pointer :: livestemc_xfer_to_livestemc(:) 
   real(r8), pointer :: livestem_mr(:)  
   real(r8), pointer :: m_deadcrootc_storage_to_fire(:) 
   real(r8), pointer :: m_deadcrootc_storage_to_litter(:) 
   real(r8), pointer :: m_deadcrootc_to_fire(:)         
   real(r8), pointer :: m_deadcrootc_to_litter(:)           
   real(r8), pointer :: m_deadcrootc_to_litter_fire(:)         
   real(r8), pointer :: m_deadcrootc_xfer_to_fire(:)
   real(r8), pointer :: m_deadcrootc_xfer_to_litter(:)
   real(r8), pointer :: m_deadstemc_storage_to_fire(:)  
   real(r8), pointer :: m_deadstemc_storage_to_litter(:)  
   real(r8), pointer :: m_deadstemc_to_fire(:)
   real(r8), pointer :: m_deadstemc_to_litter(:)            
   real(r8), pointer :: m_deadstemc_to_litter_fire(:)
   real(r8), pointer :: m_deadstemc_xfer_to_fire(:) 
   real(r8), pointer :: m_deadstemc_xfer_to_litter(:) 
   real(r8), pointer :: m_frootc_storage_to_fire(:)     
   real(r8), pointer :: m_frootc_storage_to_litter(:)     
   real(r8), pointer :: m_frootc_to_fire(:)             
   real(r8), pointer :: m_frootc_to_litter(:)             
   real(r8), pointer :: m_frootc_xfer_to_fire(:)    
   real(r8), pointer :: m_frootc_xfer_to_litter(:)    
   real(r8), pointer :: m_gresp_storage_to_fire(:)      
   real(r8), pointer :: m_gresp_storage_to_litter(:)      
   real(r8), pointer :: m_gresp_xfer_to_fire(:)    
   real(r8), pointer :: m_gresp_xfer_to_litter(:)
   real(r8), pointer :: m_leafc_storage_to_fire(:)      
   real(r8), pointer :: m_leafc_storage_to_litter(:)      
   real(r8), pointer :: m_leafc_to_fire(:)             
   real(r8), pointer :: m_leafc_to_litter(:)
   real(r8), pointer :: m_leafc_xfer_to_fire(:)     
   real(r8), pointer :: m_leafc_xfer_to_litter(:)    
   real(r8), pointer :: m_livecrootc_storage_to_fire(:) 
   real(r8), pointer :: m_livecrootc_storage_to_litter(:) 
   real(r8), pointer :: m_livecrootc_to_fire(:)         
   real(r8), pointer :: m_livecrootc_to_litter(:)           
   real(r8), pointer :: m_livecrootc_xfer_to_fire(:)
   real(r8), pointer :: m_livecrootc_xfer_to_litter(:)
   real(r8), pointer :: m_livestemc_storage_to_fire(:)  
   real(r8), pointer :: m_livestemc_storage_to_litter(:)  
   real(r8), pointer :: m_livestemc_to_fire(:)          
   real(r8), pointer :: m_livestemc_to_litter(:)            
   real(r8), pointer :: m_livestemc_xfer_to_fire(:) 
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
   real(r8), pointer :: tempsum_litfall(:)      !temporary annual sum of litfall (gC/m2/yr)

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
   real(r8), pointer :: litr1c_to_soil1c(:)   ! litter1 C loss to soil1 (gC/m2/s)
   real(r8), pointer :: litr2c_to_soil2c(:)   ! litter2 C loss to soil2 (gC/m2/s)
   real(r8), pointer :: litr3c_to_soil3c(:)   ! litter3 C loss to soil3 (gC/m2/s)
   real(r8), pointer :: cwdc_to_litr2c(:)     ! cwdc C to soil2 (gC/m2/s)
   real(r8), pointer :: cwdc_to_litr3c(:)     ! cwdc C to soil3 (gC/m2/s)
!
!
! local pointers to implicit in/out scalars
!
!
! local pointers to implicit out scalars
!
!
! !OTHER LOCAL VARIABLES:
   integer :: c,p          ! indices
   integer :: fp,fc        ! lake filter indices

!EOP
!-----------------------------------------------------------------------
   ! assign local pointers
    ivt                            => pft%itype
    col_fire_closs                 => ccf%col_fire_closs
    er                             => ccf%er
    hr                             => ccf%hr
    litfire                        => ccf%litfire
    lithr                          => ccf%lithr
    litr1_hr                       => ccf%litr1_hr
    litr2_hr                       => ccf%litr2_hr
    litr3_hr                       => ccf%litr3_hr
    m_cwdc_to_fire                 => ccf%m_cwdc_to_fire
    m_litr1c_to_fire               => ccf%m_litr1c_to_fire
    m_litr2c_to_fire               => ccf%m_litr2c_to_fire
    m_litr3c_to_fire               => ccf%m_litr3c_to_fire
    cwdc_to_litr2c                 => ccf%cwdc_to_litr2c
    cwdc_to_litr3c                 => ccf%cwdc_to_litr3c
    litr1c_to_soil1c               => ccf%litr1c_to_soil1c
    litr2c_to_soil2c               => ccf%litr2c_to_soil2c
    litr3c_to_soil3c               => ccf%litr3c_to_soil3c
    nee                            => ccf%nee
    nep                            => ccf%nep
    nbp                            => ccf%nbp
    col_ar                         => pcf_a%ar
    col_gpp                        => pcf_a%gpp
    col_npp                        => pcf_a%npp
    col_pft_fire_closs             => pcf_a%pft_fire_closs
    col_litfall                    => pcf_a%litfall
    col_rr                         => pcf_a%rr
    col_vegfire                    => pcf_a%vegfire
    col_wood_harvestc              => pcf_a%wood_harvestc
    soil1_hr                       => ccf%soil1_hr
    soil2_hr                       => ccf%soil2_hr
    soil3_hr                       => ccf%soil3_hr
    soil4_hr                       => ccf%soil4_hr
    somfire                        => ccf%somfire
    somhr                          => ccf%somhr
    sr                             => ccf%sr
    totfire                        => ccf%totfire
    cwdc_hr                        => ccf%cwdc_hr
    cwdc_loss                      => ccf%cwdc_loss
    litterc_loss                   => ccf%litterc_loss
    ! dynamic landcover pointers
    dwt_closs                      => ccf%dwt_closs
    landuseflux                    => ccf%landuseflux
    landuptake                     => ccf%landuptake
    dwt_conv_cflux                 => ccf%dwt_conv_cflux
    seedc                          => ccs%seedc
    
    ! wood product pointers
    prod10c_loss                   => ccf%prod10c_loss
    prod100c_loss                  => ccf%prod100c_loss
    product_closs                  => ccf%product_closs
    prod10c                        => ccs%prod10c
    prod100c                       => ccs%prod100c
    totprodc                       => ccs%totprodc
    
    cwdc                           => ccs%cwdc
    litr1c                         => ccs%litr1c
    litr2c                         => ccs%litr2c
    litr3c                         => ccs%litr3c
    col_totpftc                    => pcs_a%totpftc
    col_totvegc                    => pcs_a%totvegc
    soil1c                         => ccs%soil1c
    soil2c                         => ccs%soil2c
    soil3c                         => ccs%soil3c
    soil4c                         => ccs%soil4c
    col_ctrunc                     => ccs%col_ctrunc
    totcolc                        => ccs%totcolc
    totecosysc                     => ccs%totecosysc
    totlitc                        => ccs%totlitc
    totsomc                        => ccs%totsomc
    agnpp                          => pcf%agnpp
    ar                             => pcf%ar
    bgnpp                          => pcf%bgnpp
    xsmrpool_to_atm                => pcf%xsmrpool_to_atm
    cpool_grain_gr                 => pcf%cpool_grain_gr
    cpool_grain_storage_gr         => pcf%cpool_grain_storage_gr
    cpool_to_grainc                => pcf%cpool_to_grainc
    grainc_xfer_to_grainc          => pcf%grainc_xfer_to_grainc
    transfer_grain_gr              => pcf%transfer_grain_gr
    grainc_to_food                 => pcf%grainc_to_food
    livestemc_to_litter            => pcf%livestemc_to_litter
    grainc                         => pcs%grainc
    grainc_storage                 => pcs%grainc_storage
    grainc_xfer                    => pcs%grainc_xfer
    cpool_deadcroot_gr             => pcf%cpool_deadcroot_gr
    cpool_deadcroot_storage_gr     => pcf%cpool_deadcroot_storage_gr
    cpool_deadstem_gr              => pcf%cpool_deadstem_gr
    cpool_deadstem_storage_gr      => pcf%cpool_deadstem_storage_gr
    cpool_froot_gr                 => pcf%cpool_froot_gr
    cpool_froot_storage_gr         => pcf%cpool_froot_storage_gr
    cpool_leaf_gr                  => pcf%cpool_leaf_gr
    cpool_leaf_storage_gr          => pcf%cpool_leaf_storage_gr
    cpool_livecroot_gr             => pcf%cpool_livecroot_gr
    cpool_livecroot_storage_gr     => pcf%cpool_livecroot_storage_gr
    cpool_livestem_gr              => pcf%cpool_livestem_gr
    cpool_livestem_storage_gr      => pcf%cpool_livestem_storage_gr
    cpool_to_deadcrootc            => pcf%cpool_to_deadcrootc
    cpool_to_deadstemc             => pcf%cpool_to_deadstemc
    cpool_to_frootc                => pcf%cpool_to_frootc
    cpool_to_leafc                 => pcf%cpool_to_leafc
    cpool_to_livecrootc            => pcf%cpool_to_livecrootc
    cpool_to_livestemc             => pcf%cpool_to_livestemc
    current_gr                     => pcf%current_gr
    deadcrootc_xfer_to_deadcrootc  => pcf%deadcrootc_xfer_to_deadcrootc
    deadstemc_xfer_to_deadstemc    => pcf%deadstemc_xfer_to_deadstemc
    frootc_to_litter               => pcf%frootc_to_litter
    frootc_xfer_to_frootc          => pcf%frootc_xfer_to_frootc
    froot_mr                       => pcf%froot_mr
    gpp                            => pcf%gpp
    gr                             => pcf%gr
    leafc_to_litter                => pcf%leafc_to_litter
    leafc_xfer_to_leafc            => pcf%leafc_xfer_to_leafc
    leaf_mr                        => pcf%leaf_mr
    litfall                        => pcf%litfall
    livecrootc_xfer_to_livecrootc  => pcf%livecrootc_xfer_to_livecrootc
    livecroot_mr                   => pcf%livecroot_mr
    livestemc_xfer_to_livestemc    => pcf%livestemc_xfer_to_livestemc
    livestem_mr                    => pcf%livestem_mr
    m_deadcrootc_storage_to_fire   => pcf%m_deadcrootc_storage_to_fire
    m_deadcrootc_storage_to_litter => pcf%m_deadcrootc_storage_to_litter
    m_deadcrootc_to_fire           => pcf%m_deadcrootc_to_fire
    m_deadcrootc_to_litter         => pcf%m_deadcrootc_to_litter
    m_deadcrootc_to_litter_fire    => pcf%m_deadcrootc_to_litter_fire
    m_deadcrootc_xfer_to_fire      => pcf%m_deadcrootc_xfer_to_fire
    m_deadcrootc_xfer_to_litter    => pcf%m_deadcrootc_xfer_to_litter
    m_deadstemc_storage_to_fire    => pcf%m_deadstemc_storage_to_fire
    m_deadstemc_storage_to_litter  => pcf%m_deadstemc_storage_to_litter
    m_deadstemc_to_fire            => pcf%m_deadstemc_to_fire
    m_deadstemc_to_litter          => pcf%m_deadstemc_to_litter
    m_deadstemc_to_litter_fire     => pcf%m_deadstemc_to_litter_fire
    m_deadstemc_xfer_to_fire       => pcf%m_deadstemc_xfer_to_fire
    m_deadstemc_xfer_to_litter     => pcf%m_deadstemc_xfer_to_litter
    m_frootc_storage_to_fire       => pcf%m_frootc_storage_to_fire
    m_frootc_storage_to_litter     => pcf%m_frootc_storage_to_litter
    m_frootc_to_fire               => pcf%m_frootc_to_fire
    m_frootc_to_litter             => pcf%m_frootc_to_litter
    m_frootc_xfer_to_fire          => pcf%m_frootc_xfer_to_fire
    m_frootc_xfer_to_litter        => pcf%m_frootc_xfer_to_litter
    m_gresp_storage_to_fire        => pcf%m_gresp_storage_to_fire
    m_gresp_storage_to_litter      => pcf%m_gresp_storage_to_litter
    m_gresp_xfer_to_fire           => pcf%m_gresp_xfer_to_fire
    m_gresp_xfer_to_litter         => pcf%m_gresp_xfer_to_litter
    m_leafc_storage_to_fire        => pcf%m_leafc_storage_to_fire
    m_leafc_storage_to_litter      => pcf%m_leafc_storage_to_litter
    m_leafc_to_fire                => pcf%m_leafc_to_fire
    m_leafc_to_litter              => pcf%m_leafc_to_litter
    m_leafc_xfer_to_fire           => pcf%m_leafc_xfer_to_fire
    m_leafc_xfer_to_litter         => pcf%m_leafc_xfer_to_litter
    m_livecrootc_storage_to_fire   => pcf%m_livecrootc_storage_to_fire
    m_livecrootc_storage_to_litter => pcf%m_livecrootc_storage_to_litter
    m_livecrootc_to_fire           => pcf%m_livecrootc_to_fire
    m_livecrootc_to_litter         => pcf%m_livecrootc_to_litter
    m_livecrootc_xfer_to_fire      => pcf%m_livecrootc_xfer_to_fire
    m_livecrootc_xfer_to_litter    => pcf%m_livecrootc_xfer_to_litter
    m_livestemc_storage_to_fire    => pcf%m_livestemc_storage_to_fire
    m_livestemc_storage_to_litter  => pcf%m_livestemc_storage_to_litter
    m_livestemc_to_fire            => pcf%m_livestemc_to_fire
    m_livestemc_to_litter          => pcf%m_livestemc_to_litter
    m_livestemc_xfer_to_fire       => pcf%m_livestemc_xfer_to_fire
    m_livestemc_xfer_to_litter     => pcf%m_livestemc_xfer_to_litter
    hrv_leafc_to_litter               => pcf%hrv_leafc_to_litter               
    hrv_leafc_storage_to_litter       => pcf%hrv_leafc_storage_to_litter     
    hrv_leafc_xfer_to_litter          => pcf%hrv_leafc_xfer_to_litter        
    hrv_frootc_to_litter              => pcf%hrv_frootc_to_litter            
    hrv_frootc_storage_to_litter      => pcf%hrv_frootc_storage_to_litter    
    hrv_frootc_xfer_to_litter         => pcf%hrv_frootc_xfer_to_litter       
    hrv_livestemc_to_litter           => pcf%hrv_livestemc_to_litter         
    hrv_livestemc_storage_to_litter   => pcf%hrv_livestemc_storage_to_litter 
    hrv_livestemc_xfer_to_litter      => pcf%hrv_livestemc_xfer_to_litter    
    hrv_deadstemc_to_prod10c          => pcf%hrv_deadstemc_to_prod10c        
    hrv_deadstemc_to_prod100c         => pcf%hrv_deadstemc_to_prod100c       
    hrv_deadstemc_storage_to_litter   => pcf%hrv_deadstemc_storage_to_litter 
    hrv_deadstemc_xfer_to_litter      => pcf%hrv_deadstemc_xfer_to_litter    
    hrv_livecrootc_to_litter          => pcf%hrv_livecrootc_to_litter        
    hrv_livecrootc_storage_to_litter  => pcf%hrv_livecrootc_storage_to_litter
    hrv_livecrootc_xfer_to_litter     => pcf%hrv_livecrootc_xfer_to_litter   
    hrv_deadcrootc_to_litter          => pcf%hrv_deadcrootc_to_litter        
    hrv_deadcrootc_storage_to_litter  => pcf%hrv_deadcrootc_storage_to_litter
    hrv_deadcrootc_xfer_to_litter     => pcf%hrv_deadcrootc_xfer_to_litter   
    hrv_gresp_storage_to_litter       => pcf%hrv_gresp_storage_to_litter     
    hrv_gresp_xfer_to_litter          => pcf%hrv_gresp_xfer_to_litter        
    hrv_xsmrpool_to_atm               => pcf%hrv_xsmrpool_to_atm             
    col_hrv_xsmrpool_to_atm           => pcf_a%hrv_xsmrpool_to_atm             
    mr                             => pcf%mr
    npp                            => pcf%npp
    pft_fire_closs                 => pcf%pft_fire_closs
    psnshade_to_cpool              => pcf%psnshade_to_cpool
    psnsun_to_cpool                => pcf%psnsun_to_cpool
    rr                             => pcf%rr
    storage_gr                     => pcf%storage_gr
    transfer_deadcroot_gr          => pcf%transfer_deadcroot_gr
    transfer_deadstem_gr           => pcf%transfer_deadstem_gr
    transfer_froot_gr              => pcf%transfer_froot_gr
    transfer_gr                    => pcf%transfer_gr
    transfer_leaf_gr               => pcf%transfer_leaf_gr
    transfer_livecroot_gr          => pcf%transfer_livecroot_gr
    transfer_livestem_gr           => pcf%transfer_livestem_gr
    vegfire                        => pcf%vegfire
    wood_harvestc                  => pcf%wood_harvestc
    frootc_alloc                   => pcf%frootc_alloc
    frootc_loss                    => pcf%frootc_loss
    leafc_alloc                    => pcf%leafc_alloc
    leafc_loss                     => pcf%leafc_loss
    woodc_alloc                    => pcf%woodc_alloc
    woodc_loss                     => pcf%woodc_loss
    cpool                          => pcs%cpool
    xsmrpool                       => pcs%xsmrpool
	pft_ctrunc                     => pcs%pft_ctrunc
    deadcrootc                     => pcs%deadcrootc
    deadcrootc_storage             => pcs%deadcrootc_storage
    deadcrootc_xfer                => pcs%deadcrootc_xfer
    deadstemc                      => pcs%deadstemc
    deadstemc_storage              => pcs%deadstemc_storage
    deadstemc_xfer                 => pcs%deadstemc_xfer
    dispvegc                       => pcs%dispvegc
    frootc                         => pcs%frootc
    frootc_storage                 => pcs%frootc_storage
    frootc_xfer                    => pcs%frootc_xfer
    gresp_storage                  => pcs%gresp_storage
    gresp_xfer                     => pcs%gresp_xfer
    leafc                          => pcs%leafc
    leafc_storage                  => pcs%leafc_storage
    leafc_xfer                     => pcs%leafc_xfer
    livecrootc                     => pcs%livecrootc
    livecrootc_storage             => pcs%livecrootc_storage
    livecrootc_xfer                => pcs%livecrootc_xfer
    livestemc                      => pcs%livestemc
    livestemc_storage              => pcs%livestemc_storage
    livestemc_xfer                 => pcs%livestemc_xfer
    storvegc                       => pcs%storvegc
    totpftc                        => pcs%totpftc
    totvegc                        => pcs%totvegc
    woodc                          => pcs%woodc
    tempsum_npp                    => pepv%tempsum_npp
    tempsum_litfall                => pepv%tempsum_litfall

   ! pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)

      ! calculate pft-level summary carbon fluxes and states

      ! gross primary production (GPP)
      gpp(p) = &
         psnsun_to_cpool(p) + &
         psnshade_to_cpool(p)

      ! maintenance respiration (MR)
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

      ! GR is the sum of current + transfer + storage GR
      gr(p) = &
         current_gr(p)  + &
         transfer_gr(p) + &
         storage_gr(p)

      ! autotrophic respiration (AR)
      if ( ivt(p) >= npcropmin )then
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
      tempsum_npp(p) = tempsum_npp(p) + npp(p)

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
         m_deadstemc_to_litter_fire(p)      + &
         m_deadcrootc_to_litter_fire(p)     + &
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
                 
      if (use_cndv) then
         ! update the annual litfall accumulator, for use in mortality code
         tempsum_litfall(p) = tempsum_litfall(p) + leafc_to_litter(p) + frootc_to_litter(p)
      end if

      ! pft-level fire losses (VEGFIRE)
      vegfire(p) = 0._r8
      
      ! pft-level wood harvest
      wood_harvestc(p) = &
         hrv_deadstemc_to_prod10c(p) + &
         hrv_deadstemc_to_prod100c(p)

      ! pft-level carbon losses to fire
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

      if ( crop_prog .and. ivt(p) >= nc3crop )then
         current_gr(p) = current_gr(p) + &
            cpool_grain_gr(p)
         storvegc(p) = storvegc(p) + &
            grainc_storage(p)  + &
            grainc_xfer(p)
         transfer_gr(p) = transfer_gr(p) + &
            transfer_grain_gr(p)
         storage_gr(p) = storage_gr(p) + &
            cpool_grain_storage_gr(p)
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
              
      ! (FROOTC_LOSS) - fine root C loss
      frootc_loss(p) = &
        m_frootc_to_litter(p)   + &
        m_frootc_to_fire(p)     + &
        hrv_frootc_to_litter(p) + &
        frootc_to_litter(p)
      
      ! (LEAFC_ALLOC) - leaf C allocation
      leafc_alloc(p) = &
        leafc_xfer_to_leafc(p)    + &
        cpool_to_leafc(p)     

      ! (LEAFC_LOSS) - leaf C loss
      leafc_loss(p) = &
        m_leafc_to_litter(p)   + &
        m_leafc_to_fire(p)     + &
        hrv_leafc_to_litter(p) + &
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

   ! column loop
   do fc = 1,num_soilc
      c = filter_soilc(fc)

      ! litter heterotrophic respiration (LITHR)
      lithr(c) = &
         litr1_hr(c) + &
         litr2_hr(c) + &
         litr3_hr(c)

      ! soil organic matter heterotrophic respiration (SOMHR)
      somhr(c) = &
         soil1_hr(c) + &
         soil2_hr(c) + &
         soil3_hr(c) + &
         soil4_hr(c)

      ! total heterotrophic respiration (HR)
      hr(c) = lithr(c) + somhr(c)

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

      ! column-level carbon losses to fire, including pft losses
      col_fire_closs(c) = &
         m_litr1c_to_fire(c)  + &
         m_litr2c_to_fire(c)  + &
         m_litr3c_to_fire(c)  + &
         m_cwdc_to_fire(c)    + &
         col_pft_fire_closs(c)
      
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

      ! total litter carbon (TOTLITC)
      totlitc(c) = &
         litr1c(c) + &
         litr2c(c) + &
         litr3c(c)

      ! total soil organic matter carbon (TOTSOMC)
      totsomc(c) = &
         soil1c(c) + &
         soil2c(c) + &
         soil3c(c) + &
         soil4c(c)
      
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
      
      ! (CWDC_LOSS) - coarse woody debris C loss
      cwdc_loss(c) = & 
        m_cwdc_to_fire(c) + &
        cwdc_to_litr2c(c) + &
        cwdc_to_litr3c(c)
      
      ! (LITTERC_LOSS) - litter C loss
      litterc_loss(c) = &
        lithr(c)            + &
        m_litr1c_to_fire(c) + &
        m_litr2c_to_fire(c) + &
        m_litr3c_to_fire(c) + &
        litr1c_to_soil1c(c) + &
        litr2c_to_soil2c(c) + &
        litr3c_to_soil3c(c)
      
   end do ! end of columns loop


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
   real(r8), pointer :: m_cwdn_to_fire(:)              
   real(r8), pointer :: m_litr1n_to_fire(:)             
   real(r8), pointer :: m_litr2n_to_fire(:)             
   real(r8), pointer :: m_litr3n_to_fire(:)             
   real(r8), pointer :: col_pft_fire_nloss(:) ! (gN/m2/s) total pft-level fire C loss 
   real(r8), pointer :: sminn_to_denit_excess(:)
   real(r8), pointer :: sminn_to_denit_l1s1(:)
   real(r8), pointer :: sminn_to_denit_l2s2(:)
   real(r8), pointer :: sminn_to_denit_l3s3(:)
   real(r8), pointer :: sminn_to_denit_s1s2(:)
   real(r8), pointer :: sminn_to_denit_s2s3(:)
   real(r8), pointer :: sminn_to_denit_s3s4(:)
   real(r8), pointer :: sminn_to_denit_s4(:)  
   real(r8), pointer :: cwdn(:)               ! (gN/m2) coarse woody debris N
   real(r8), pointer :: litr1n(:)             ! (gN/m2) litter labile N
   real(r8), pointer :: litr2n(:)             ! (gN/m2) litter cellulose N
   real(r8), pointer :: litr3n(:)             ! (gN/m2) litter lignin N
   real(r8), pointer :: col_totpftn(:)        ! (gN/m2) total pft-level nitrogen
   real(r8), pointer :: col_totvegn(:)            ! (gN/m2) total vegetation nitrogen
   real(r8), pointer :: sminn(:)              ! (gN/m2) soil mineral N
   real(r8), pointer :: soil1n(:)             ! (gN/m2) soil organic matter N (fast pool)
   real(r8), pointer :: soil2n(:)             ! (gN/m2) soil organic matter N (medium pool)
   real(r8), pointer :: soil3n(:)             ! (gN/m2) soil orgainc matter N (slow pool)
   real(r8), pointer :: soil4n(:)             ! (gN/m2) soil orgainc matter N (slowest pool)
   real(r8), pointer :: col_ntrunc(:)         ! (gN/m2) column-level sink for N truncation
   real(r8), pointer :: totcoln(:)            ! (gN/m2) total column nitrogen, incl veg
   real(r8), pointer :: totecosysn(:)         ! (gN/m2) total ecosystem nitrogen, incl veg 
   real(r8), pointer :: totlitn(:)            ! (gN/m2) total litter nitrogen
   real(r8), pointer :: totsomn(:)            ! (gN/m2) total soil organic matter nitrogen
   real(r8), pointer :: m_deadcrootn_storage_to_fire(:) 
   real(r8), pointer :: m_deadcrootn_to_fire(:)         
   real(r8), pointer :: m_deadcrootn_xfer_to_fire(:)
   real(r8), pointer :: m_deadstemn_storage_to_fire(:)  
   real(r8), pointer :: m_deadstemn_to_fire(:)          
   real(r8), pointer :: m_deadstemn_xfer_to_fire(:) 
   real(r8), pointer :: m_frootn_storage_to_fire(:)     
   real(r8), pointer :: m_frootn_to_fire(:)             
   real(r8), pointer :: m_frootn_xfer_to_fire(:)    
   real(r8), pointer :: m_leafn_storage_to_fire(:)      
   real(r8), pointer :: m_leafn_to_fire(:)              
   real(r8), pointer :: m_leafn_xfer_to_fire(:)     
   real(r8), pointer :: m_livecrootn_storage_to_fire(:) 
   real(r8), pointer :: m_livecrootn_to_fire(:)         
   real(r8), pointer :: m_livecrootn_xfer_to_fire(:)
   real(r8), pointer :: m_livestemn_storage_to_fire(:)  
   real(r8), pointer :: m_livestemn_to_fire(:)          
   real(r8), pointer :: m_livestemn_xfer_to_fire(:) 
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
!
! local pointers to implicit in/out scalars
!
! local pointers to implicit out scalars
!
! !OTHER LOCAL VARIABLES:
   integer :: c,p         ! indices
   integer :: fp,fc       ! lake filter indices

!EOP
!-----------------------------------------------------------------------
    ! assign local pointers
    ivt                            => pft%itype
    col_fire_nloss                 => cnf%col_fire_nloss
    denit                          => cnf%denit
    m_cwdn_to_fire                 => cnf%m_cwdn_to_fire
    m_litr1n_to_fire               => cnf%m_litr1n_to_fire
    m_litr2n_to_fire               => cnf%m_litr2n_to_fire
    m_litr3n_to_fire               => cnf%m_litr3n_to_fire
    col_pft_fire_nloss             => pnf_a%pft_fire_nloss
    sminn_to_denit_excess          => cnf%sminn_to_denit_excess
    sminn_to_denit_l1s1            => cnf%sminn_to_denit_l1s1
    sminn_to_denit_l2s2            => cnf%sminn_to_denit_l2s2
    sminn_to_denit_l3s3            => cnf%sminn_to_denit_l3s3
    sminn_to_denit_s1s2            => cnf%sminn_to_denit_s1s2
    sminn_to_denit_s2s3            => cnf%sminn_to_denit_s2s3
    sminn_to_denit_s3s4            => cnf%sminn_to_denit_s3s4
    sminn_to_denit_s4              => cnf%sminn_to_denit_s4
    cwdn                           => cns%cwdn
    litr1n                         => cns%litr1n
    litr2n                         => cns%litr2n
    litr3n                         => cns%litr3n
    col_totpftn                    => pns_a%totpftn
    col_totvegn                    => pns_a%totvegn
    sminn                          => cns%sminn
    col_ntrunc                     => cns%col_ntrunc
    soil1n                         => cns%soil1n
    soil2n                         => cns%soil2n
    soil3n                         => cns%soil3n
    soil4n                         => cns%soil4n
    totcoln                        => cns%totcoln
    totecosysn                     => cns%totecosysn
    totlitn                        => cns%totlitn
    totsomn                        => cns%totsomn
    m_deadcrootn_storage_to_fire   => pnf%m_deadcrootn_storage_to_fire
    m_deadcrootn_to_fire           => pnf%m_deadcrootn_to_fire
    m_deadcrootn_xfer_to_fire      => pnf%m_deadcrootn_xfer_to_fire
    m_deadstemn_storage_to_fire    => pnf%m_deadstemn_storage_to_fire
    m_deadstemn_to_fire            => pnf%m_deadstemn_to_fire
    m_deadstemn_xfer_to_fire       => pnf%m_deadstemn_xfer_to_fire
    m_frootn_storage_to_fire       => pnf%m_frootn_storage_to_fire
    m_frootn_to_fire               => pnf%m_frootn_to_fire
    m_frootn_xfer_to_fire          => pnf%m_frootn_xfer_to_fire
    m_leafn_storage_to_fire        => pnf%m_leafn_storage_to_fire
    m_leafn_to_fire                => pnf%m_leafn_to_fire
    m_leafn_xfer_to_fire           => pnf%m_leafn_xfer_to_fire
    m_livecrootn_storage_to_fire   => pnf%m_livecrootn_storage_to_fire
    m_livecrootn_to_fire           => pnf%m_livecrootn_to_fire
    m_livecrootn_xfer_to_fire      => pnf%m_livecrootn_xfer_to_fire
    m_livestemn_storage_to_fire    => pnf%m_livestemn_storage_to_fire
    m_livestemn_to_fire            => pnf%m_livestemn_to_fire
    m_livestemn_xfer_to_fire       => pnf%m_livestemn_xfer_to_fire
    m_retransn_to_fire             => pnf%m_retransn_to_fire
    hrv_deadstemn_to_prod10n         => pnf%hrv_deadstemn_to_prod10n        
    hrv_deadstemn_to_prod100n        => pnf%hrv_deadstemn_to_prod100n       
    ndeploy                        => pnf%ndeploy
    pft_fire_nloss                 => pnf%pft_fire_nloss
    retransn_to_npool              => pnf%retransn_to_npool
    sminn_to_npool                 => pnf%sminn_to_npool
    deadcrootn                     => pns%deadcrootn
    deadcrootn_storage             => pns%deadcrootn_storage
    deadcrootn_xfer                => pns%deadcrootn_xfer
    deadstemn                      => pns%deadstemn
    deadstemn_storage              => pns%deadstemn_storage
    deadstemn_xfer                 => pns%deadstemn_xfer
    dispvegn                       => pns%dispvegn
    frootn                         => pns%frootn
    frootn_storage                 => pns%frootn_storage
    frootn_xfer                    => pns%frootn_xfer
    leafn                          => pns%leafn
    leafn_storage                  => pns%leafn_storage
    leafn_xfer                     => pns%leafn_xfer
    livecrootn                     => pns%livecrootn
    livecrootn_storage             => pns%livecrootn_storage
    livecrootn_xfer                => pns%livecrootn_xfer
    grainn                         => pns%grainn
    grainn_storage                 => pns%grainn_storage
    grainn_xfer                    => pns%grainn_xfer
    livestemn                      => pns%livestemn
    livestemn_storage              => pns%livestemn_storage
    livestemn_xfer                 => pns%livestemn_xfer
    retransn                       => pns%retransn
    npool                          => pns%npool
    pft_ntrunc                     => pns%pft_ntrunc
    storvegn                       => pns%storvegn
    totpftn                        => pns%totpftn
    totvegn                        => pns%totvegn
    ! dynamic landcover pointers
    wood_harvestn                  => pnf%wood_harvestn
    col_wood_harvestn              => pnf_a%wood_harvestn 
    dwt_nloss                      => cnf%dwt_nloss
    dwt_conv_nflux                 => cnf%dwt_conv_nflux
    prod10n_loss                   => cnf%prod10n_loss
    prod100n_loss                  => cnf%prod100n_loss
    product_nloss                  => cnf%product_nloss
    seedn                          => cns%seedn
    prod10n                        => cns%prod10n
    prod100n                       => cns%prod100n
    totprodn                       => cns%totprodn

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

      if ( crop_prog .and. ivt(p) >= nc3crop )then
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

   ! column loop
   do fc = 1,num_soilc
      c = filter_soilc(fc)

      ! total N denitrification (DENIT)
      denit(c) = &
         sminn_to_denit_l1s1(c) + &
         sminn_to_denit_l2s2(c) + &
         sminn_to_denit_l3s3(c) + &
         sminn_to_denit_s1s2(c) + &
         sminn_to_denit_s2s3(c) + &
         sminn_to_denit_s3s4(c) + &
         sminn_to_denit_s4(c) + &
         sminn_to_denit_excess(c)

      ! total column-level fire N losses
      col_fire_nloss(c) = &
         m_litr1n_to_fire(c) + &
         m_litr2n_to_fire(c) + &
         m_litr3n_to_fire(c) + &
         m_cwdn_to_fire(c)   + &
         col_pft_fire_nloss(c)

      ! column-level N losses due to landcover change
      dwt_nloss(c) = &
         dwt_conv_nflux(c)
         
      ! total wood product N loss
      product_nloss(c) = &
         prod10n_loss(c) + &
         prod100n_loss(c) 

      ! total litter nitrogen (TOTLITN)
      totlitn(c) = &
         litr1n(c) + &
         litr2n(c) + &
         litr3n(c)

      ! total soil organic matter nitrogen (TOTSOMN)
      totsomn(c) = &
         soil1n(c) + &
         soil2n(c) + &
         soil3n(c) + &
         soil4n(c)

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

   end do ! end of columns loop


end subroutine NSummary
!-----------------------------------------------------------------------

end module CNSummaryMod
