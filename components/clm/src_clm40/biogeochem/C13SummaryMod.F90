module C13SummaryMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: C13SummaryMod
!
! !DESCRIPTION:
! Module for isotope carbon summary calculations
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    implicit none
    save
    private
! !PUBLIC MEMBER FUNCTIONS:
    public :: C13Summary
!
! !REVISION HISTORY:
! 7/13/2005: Created by Peter Thornton
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: C13Summary
!
! !INTERFACE:
subroutine C13Summary(num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! On the radiation time step, perform pft and column-level carbon
! summary calculations
!
! !USES:
   use clmtype
   use pft2colMod, only: p2c
   use clm_varctl, only: iulog, use_c13
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
   real(r8), pointer :: col_fire_closs(:) ! (gC/m2/s) total column-level fire C loss
   real(r8), pointer :: er(:)            ! (gC/m2/s) total ecosystem respiration, autotrophic + heterotrophic
   real(r8), pointer :: hr(:)            ! (gC/m2/s) total heterotrophic respiration
   real(r8), pointer :: litfire(:)       ! (gC/m2/s) litter fire losses
   real(r8), pointer :: lithr(:)         ! (gC/m2/s) litter heterotrophic respiration 
   real(r8), pointer :: litr1_hr(:)       
   real(r8), pointer :: litr2_hr(:)        
   real(r8), pointer :: litr3_hr(:)        
   real(r8), pointer :: m_cwdc_to_fire(:)
   real(r8), pointer :: m_litr1c_to_fire(:)             
   real(r8), pointer :: m_litr2c_to_fire(:)             
   real(r8), pointer :: m_litr3c_to_fire(:)             
   real(r8), pointer :: nee(:)           ! (gC/m2/s) net ecosystem exchange of carbon, includes fire, land-use, and wood products flux, positive for source
   real(r8), pointer :: nbp(:)           ! (gC/m2/s) net biome production, includes fire, land-use, and wood products flux, positive for sink
   real(r8), pointer :: nep(:)           ! (gC/m2/s) net ecosystem production, excludes fire, land-use, and wood products flux, positive for sink
   real(r8), pointer :: col_ar(:)             ! (gC/m2/s) autotrophic respiration (MR + GR)
   real(r8), pointer :: col_gpp(:)                  !GPP flux before downregulation (gC/m2/s)
   real(r8), pointer :: col_npp(:)            ! (gC/m2/s) net primary production
   real(r8), pointer :: col_pft_fire_closs(:) ! (gC/m2/s) total pft-level fire C loss 
   real(r8), pointer :: col_rr(:)             ! (gC/m2/s) root respiration (fine root MR + total root GR)
   real(r8), pointer :: col_vegfire(:)        ! (gC/m2/s) pft-level fire loss (obsolete, mark for removal)
   real(r8), pointer :: col_wood_harvestc(:)
   real(r8), pointer :: soil1_hr(:)        
   real(r8), pointer :: soil2_hr(:)        
   real(r8), pointer :: soil3_hr(:) 
   real(r8), pointer :: soil4_hr(:) 
   real(r8), pointer :: somfire(:)       ! (gC/m2/s) soil organic matter fire losses
   real(r8), pointer :: somhr(:)         ! (gC/m2/s) soil organic matter heterotrophic respiration
   real(r8), pointer :: sr(:)            ! (gC/m2/s) total soil respiration (HR + root resp)
   real(r8), pointer :: totfire(:)       ! (gC/m2/s) total ecosystem fire losses
   real(r8), pointer :: cwdc(:)               ! (gC/m2) coarse woody debris C
   real(r8), pointer :: litr1c(:)             ! (gC/m2) litter labile C
   real(r8), pointer :: litr2c(:)             ! (gC/m2) litter cellulose C
   real(r8), pointer :: litr3c(:)             ! (gC/m2) litter lignin C
   real(r8), pointer :: col_totpftc(:)        ! (gC/m2) total pft-level carbon, including cpool
   real(r8), pointer :: col_totvegc(:)        ! (gC/m2) total vegetation carbon, excluding cpool
   real(r8), pointer :: soil1c(:)             ! (gC/m2) soil organic matter C (fast pool)
   real(r8), pointer :: soil2c(:)             ! (gC/m2) soil organic matter C (medium pool)
   real(r8), pointer :: soil3c(:)             ! (gC/m2) soil organic matter C (slow pool)
   real(r8), pointer :: soil4c(:)             ! (gC/m2) soil organic matter C (slow pool)
   real(r8), pointer :: totcolc(:)            ! (gC/m2) total column carbon, incl veg and cpool
   real(r8), pointer :: totecosysc(:)         ! (gC/m2) total ecosystem carbon, incl veg but excl cpool
   real(r8), pointer :: totlitc(:)            ! (gC/m2) total litter carbon
   real(r8), pointer :: totsomc(:)            ! (gC/m2) total soil organic matter carbon
   real(r8), pointer :: agnpp(:)          ! (gC/m2/s) aboveground NPP
   real(r8), pointer :: ar(:)             ! (gC/m2/s) autotrophic respiration (MR + GR)
   real(r8), pointer :: bgnpp(:)          ! (gC/m2/s) belowground NPP
   real(r8), pointer :: cpool_deadcroot_gr(:)        
   real(r8), pointer :: cpool_deadcroot_storage_gr(:)
   real(r8), pointer :: cpool_deadstem_gr(:)         
   real(r8), pointer :: cpool_deadstem_storage_gr(:) 
   real(r8), pointer :: cpool_froot_gr(:)            
   real(r8), pointer :: cpool_froot_storage_gr(:)    
   real(r8), pointer :: cpool_leaf_gr(:)             
   real(r8), pointer :: cpool_leaf_storage_gr(:)     
   real(r8), pointer :: cpool_livecroot_gr(:)        
   real(r8), pointer :: cpool_livecroot_storage_gr(:)
   real(r8), pointer :: cpool_livestem_gr(:)         
   real(r8), pointer :: cpool_livestem_storage_gr(:) 
   real(r8), pointer :: cpool_to_deadcrootc(:)        
   real(r8), pointer :: cpool_to_deadstemc(:)         
   real(r8), pointer :: cpool_to_frootc(:)            
   real(r8), pointer :: cpool_to_leafc(:)             
   real(r8), pointer :: cpool_to_livecrootc(:)        
   real(r8), pointer :: cpool_to_livestemc(:)         
   real(r8), pointer :: current_gr(:)     ! (gC/m2/s) growth resp for new growth displayed in this timestep
   real(r8), pointer :: deadcrootc_xfer_to_deadcrootc(:)
   real(r8), pointer :: deadstemc_xfer_to_deadstemc(:) 
   real(r8), pointer :: frootc_to_litter(:)
   real(r8), pointer :: frootc_xfer_to_frootc(:)       
   real(r8), pointer :: froot_mr(:)     
   real(r8), pointer :: froot_curmr(:)     
   real(r8), pointer :: froot_xsmr(:)     
   real(r8), pointer :: gpp(:)                  !GPP flux before downregulation (gC/m2/s)
   real(r8), pointer :: gr(:)             ! (gC/m2/s) total growth respiration
   real(r8), pointer :: leafc_to_litter(:)
   real(r8), pointer :: leafc_xfer_to_leafc(:)         
   real(r8), pointer :: leaf_mr(:)
   real(r8), pointer :: leaf_curmr(:)
   real(r8), pointer :: leaf_xsmr(:)
   real(r8), pointer :: litfall(:)        ! (gC/m2/s) litterfall (leaves and fine roots)
   real(r8), pointer :: livecrootc_xfer_to_livecrootc(:)
   real(r8), pointer :: livecroot_mr(:)
   real(r8), pointer :: livecroot_curmr(:)
   real(r8), pointer :: livecroot_xsmr(:)
   real(r8), pointer :: livestemc_xfer_to_livestemc(:) 
   real(r8), pointer :: livestem_mr(:)  
   real(r8), pointer :: livestem_curmr(:)  
   real(r8), pointer :: livestem_xsmr(:)  
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
   real(r8), pointer :: mr(:)             ! (gC/m2/s) maintenance respiration
   real(r8), pointer :: npp(:)            ! (gC/m2/s) net primary production
   real(r8), pointer :: pft_fire_closs(:) ! (gC/m2/s) total pft-level fire C loss 
   real(r8), pointer :: psnshade_to_cpool(:)
   real(r8), pointer :: psnsun_to_cpool(:) 
   real(r8), pointer :: rr(:)             ! (gC/m2/s) root respiration (fine root MR + total root GR)
   real(r8), pointer :: storage_gr(:)     ! (gC/m2/s) growth resp for growth sent to storage for later display
   real(r8), pointer :: transfer_deadcroot_gr(:)
   real(r8), pointer :: transfer_deadstem_gr(:)      
   real(r8), pointer :: transfer_froot_gr(:)         
   real(r8), pointer :: transfer_gr(:)    ! (gC/m2/s) growth resp for transfer growth displayed in this timestep
   real(r8), pointer :: transfer_leaf_gr(:)          
   real(r8), pointer :: transfer_livecroot_gr(:)     
   real(r8), pointer :: transfer_livestem_gr(:)      
   real(r8), pointer :: wood_harvestc(:)      ! (gC/m2/s) pft-level wood harvest (to product pools)
   real(r8), pointer :: vegfire(:)        ! (gC/m2/s) pft-level fire loss (obsolete, mark for removal)
   real(r8), pointer :: cpool(:)              ! (gC/m2) temporary photosynthate C pool
   real(r8), pointer :: xsmrpool(:)              ! (gC/m2) temporary photosynthate C pool
   real(r8), pointer :: deadcrootc(:)         ! (gC/m2) dead coarse root C
   real(r8), pointer :: deadcrootc_storage(:) ! (gC/m2) dead coarse root C storage
   real(r8), pointer :: deadcrootc_xfer(:)    !(gC/m2) dead coarse root C transfer
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
   real(r8), pointer :: livecrootc_xfer(:)    !(gC/m2) live coarse root C transfer
   real(r8), pointer :: livestemc(:)          ! (gC/m2) live stem C
   real(r8), pointer :: livestemc_storage(:)  ! (gC/m2) live stem C storage
   real(r8), pointer :: livestemc_xfer(:)     ! (gC/m2) live stem C transfer
   real(r8), pointer :: storvegc(:)           ! (gC/m2) stored vegetation carbon, excluding cpool
   real(r8), pointer :: totpftc(:)            ! (gC/m2) total pft-level carbon, including cpool
   real(r8), pointer :: totvegc(:)            ! (gC/m2) total vegetation carbon, excluding cpool
   ! for landcover change
   real(r8), pointer :: dwt_closs(:)          ! (gC/m2/s) total carbon loss from product pools and conversion
   real(r8), pointer :: dwt_conv_cflux(:)     ! (gC/m2/s) conversion C flux (immediate loss to atm)
   real(r8), pointer :: prod10c_loss(:)       ! (gC/m2/s) loss from 10-yr wood product pool
   real(r8), pointer :: prod100c_loss(:)      ! (gC/m2/s) loss from 100-yr wood product pool
   real(r8), pointer :: product_closs(:)      ! (gC/m2/s) total wood product carbon loss
   real(r8), pointer :: prod10c(:)            ! (gC/m2) wood product C pool, 10-year lifespan
   real(r8), pointer :: prod100c(:)           ! (gC/m2) wood product C pool, 100-year lifespan
   real(r8), pointer :: totprodc(:)           ! (gC/m2) total wood product C
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

   if (.not. use_c13) then
      RETURN
   end if

   ! assign local pointers
    col_fire_closs                 => cc13f%col_fire_closs
    er                             => cc13f%er
    hr                             => cc13f%hr
    litfire                        => cc13f%litfire
    lithr                          => cc13f%lithr
    litr1_hr                       => cc13f%litr1_hr
    litr2_hr                       => cc13f%litr2_hr
    litr3_hr                       => cc13f%litr3_hr
    m_cwdc_to_fire                 => cc13f%m_cwdc_to_fire
    m_litr1c_to_fire               => cc13f%m_litr1c_to_fire
    m_litr2c_to_fire               => cc13f%m_litr2c_to_fire
    m_litr3c_to_fire               => cc13f%m_litr3c_to_fire
    nee                            => cc13f%nee
    nep                            => cc13f%nep
    nbp                            => cc13f%nbp
    col_ar                         => pc13f_a%ar
    col_gpp                        => pc13f_a%gpp
    col_npp                        => pc13f_a%npp
    col_pft_fire_closs             => pc13f_a%pft_fire_closs
    col_rr                         => pc13f_a%rr
    col_vegfire                    => pc13f_a%vegfire
    col_wood_harvestc              => pc13f_a%wood_harvestc
    soil1_hr                       => cc13f%soil1_hr
    soil2_hr                       => cc13f%soil2_hr
    soil3_hr                       => cc13f%soil3_hr
    soil4_hr                       => cc13f%soil4_hr
    somfire                        => cc13f%somfire
    somhr                          => cc13f%somhr
    sr                             => cc13f%sr
    totfire                        => cc13f%totfire
    
    ! dynamic landcover pointers
    dwt_closs                      => cc13f%dwt_closs
    dwt_conv_cflux                 => cc13f%dwt_conv_cflux

    ! wood product pointers
    prod10c_loss                   => cc13f%prod10c_loss
    prod100c_loss                  => cc13f%prod100c_loss
    product_closs                  => cc13f%product_closs
    prod10c                        => cc13s%prod10c
    prod100c                       => cc13s%prod100c
    totprodc                       => cc13s%totprodc
    
    cwdc                           => cc13s%cwdc
    litr1c                         => cc13s%litr1c
    litr2c                         => cc13s%litr2c
    litr3c                         => cc13s%litr3c
    col_totpftc                    => pc13s_a%totpftc
    col_totvegc                    => pc13s_a%totvegc
    soil1c                         => cc13s%soil1c
    soil2c                         => cc13s%soil2c
    soil3c                         => cc13s%soil3c
    soil4c                         => cc13s%soil4c
    totcolc                        => cc13s%totcolc
    totecosysc                     => cc13s%totecosysc
    totlitc                        => cc13s%totlitc
    totsomc                        => cc13s%totsomc
    agnpp                          => pc13f%agnpp
    ar                             => pc13f%ar
    bgnpp                          => pc13f%bgnpp
    cpool_deadcroot_gr             => pc13f%cpool_deadcroot_gr
    cpool_deadcroot_storage_gr     => pc13f%cpool_deadcroot_storage_gr
    cpool_deadstem_gr              => pc13f%cpool_deadstem_gr
    cpool_deadstem_storage_gr      => pc13f%cpool_deadstem_storage_gr
    cpool_froot_gr                 => pc13f%cpool_froot_gr
    cpool_froot_storage_gr         => pc13f%cpool_froot_storage_gr
    cpool_leaf_gr                  => pc13f%cpool_leaf_gr
    cpool_leaf_storage_gr          => pc13f%cpool_leaf_storage_gr
    cpool_livecroot_gr             => pc13f%cpool_livecroot_gr
    cpool_livecroot_storage_gr     => pc13f%cpool_livecroot_storage_gr
    cpool_livestem_gr              => pc13f%cpool_livestem_gr
    cpool_livestem_storage_gr      => pc13f%cpool_livestem_storage_gr
    cpool_to_deadcrootc            => pc13f%cpool_to_deadcrootc
    cpool_to_deadstemc             => pc13f%cpool_to_deadstemc
    cpool_to_frootc                => pc13f%cpool_to_frootc
    cpool_to_leafc                 => pc13f%cpool_to_leafc
    cpool_to_livecrootc            => pc13f%cpool_to_livecrootc
    cpool_to_livestemc             => pc13f%cpool_to_livestemc
    current_gr                     => pc13f%current_gr
    deadcrootc_xfer_to_deadcrootc  => pc13f%deadcrootc_xfer_to_deadcrootc
    deadstemc_xfer_to_deadstemc    => pc13f%deadstemc_xfer_to_deadstemc
    frootc_to_litter               => pc13f%frootc_to_litter
    frootc_xfer_to_frootc          => pc13f%frootc_xfer_to_frootc
    froot_mr                       => pc13f%froot_mr
    froot_curmr                    => pc13f%froot_curmr
    froot_xsmr                     => pc13f%froot_xsmr
    gpp                            => pc13f%gpp
    gr                             => pc13f%gr
    leafc_to_litter                => pc13f%leafc_to_litter
    leafc_xfer_to_leafc            => pc13f%leafc_xfer_to_leafc
    leaf_mr                        => pc13f%leaf_mr
    leaf_curmr                     => pc13f%leaf_curmr
    leaf_xsmr                      => pc13f%leaf_xsmr
    litfall                        => pc13f%litfall
    livecrootc_xfer_to_livecrootc  => pc13f%livecrootc_xfer_to_livecrootc
    livecroot_mr                   => pc13f%livecroot_mr
    livecroot_curmr                => pc13f%livecroot_curmr
    livecroot_xsmr                 => pc13f%livecroot_xsmr
    livestemc_xfer_to_livestemc    => pc13f%livestemc_xfer_to_livestemc
    livestem_mr                    => pc13f%livestem_mr
    livestem_curmr                 => pc13f%livestem_curmr
    livestem_xsmr                  => pc13f%livestem_xsmr
    m_deadcrootc_storage_to_fire   => pc13f%m_deadcrootc_storage_to_fire
    m_deadcrootc_storage_to_litter => pc13f%m_deadcrootc_storage_to_litter
    m_deadcrootc_to_fire           => pc13f%m_deadcrootc_to_fire
    m_deadcrootc_to_litter         => pc13f%m_deadcrootc_to_litter
    m_deadcrootc_to_litter_fire    => pc13f%m_deadcrootc_to_litter_fire
    m_deadcrootc_xfer_to_fire      => pc13f%m_deadcrootc_xfer_to_fire
    m_deadcrootc_xfer_to_litter    => pc13f%m_deadcrootc_xfer_to_litter
    m_deadstemc_storage_to_fire    => pc13f%m_deadstemc_storage_to_fire
    m_deadstemc_storage_to_litter  => pc13f%m_deadstemc_storage_to_litter
    m_deadstemc_to_fire            => pc13f%m_deadstemc_to_fire
    m_deadstemc_to_litter          => pc13f%m_deadstemc_to_litter
    m_deadstemc_to_litter_fire     => pc13f%m_deadstemc_to_litter_fire
    m_deadstemc_xfer_to_fire       => pc13f%m_deadstemc_xfer_to_fire
    m_deadstemc_xfer_to_litter     => pc13f%m_deadstemc_xfer_to_litter
    m_frootc_storage_to_fire       => pc13f%m_frootc_storage_to_fire
    m_frootc_storage_to_litter     => pc13f%m_frootc_storage_to_litter
    m_frootc_to_fire               => pc13f%m_frootc_to_fire
    m_frootc_to_litter             => pc13f%m_frootc_to_litter
    m_frootc_xfer_to_fire          => pc13f%m_frootc_xfer_to_fire
    m_frootc_xfer_to_litter        => pc13f%m_frootc_xfer_to_litter
    m_gresp_storage_to_fire        => pc13f%m_gresp_storage_to_fire
    m_gresp_storage_to_litter      => pc13f%m_gresp_storage_to_litter
    m_gresp_xfer_to_fire           => pc13f%m_gresp_xfer_to_fire
    m_gresp_xfer_to_litter         => pc13f%m_gresp_xfer_to_litter
    m_leafc_storage_to_fire        => pc13f%m_leafc_storage_to_fire
    m_leafc_storage_to_litter      => pc13f%m_leafc_storage_to_litter
    m_leafc_to_fire                => pc13f%m_leafc_to_fire
    m_leafc_to_litter              => pc13f%m_leafc_to_litter
    m_leafc_xfer_to_fire           => pc13f%m_leafc_xfer_to_fire
    m_leafc_xfer_to_litter         => pc13f%m_leafc_xfer_to_litter
    m_livecrootc_storage_to_fire   => pc13f%m_livecrootc_storage_to_fire
    m_livecrootc_storage_to_litter => pc13f%m_livecrootc_storage_to_litter
    m_livecrootc_to_fire           => pc13f%m_livecrootc_to_fire
    m_livecrootc_to_litter         => pc13f%m_livecrootc_to_litter
    m_livecrootc_xfer_to_fire      => pc13f%m_livecrootc_xfer_to_fire
    m_livecrootc_xfer_to_litter    => pc13f%m_livecrootc_xfer_to_litter
    m_livestemc_storage_to_fire    => pc13f%m_livestemc_storage_to_fire
    m_livestemc_storage_to_litter  => pc13f%m_livestemc_storage_to_litter
    m_livestemc_to_fire            => pc13f%m_livestemc_to_fire
    m_livestemc_to_litter          => pc13f%m_livestemc_to_litter
    m_livestemc_xfer_to_fire       => pc13f%m_livestemc_xfer_to_fire
    m_livestemc_xfer_to_litter     => pc13f%m_livestemc_xfer_to_litter
    hrv_leafc_to_litter               => pc13f%hrv_leafc_to_litter               
    hrv_leafc_storage_to_litter       => pc13f%hrv_leafc_storage_to_litter     
    hrv_leafc_xfer_to_litter          => pc13f%hrv_leafc_xfer_to_litter        
    hrv_frootc_to_litter              => pc13f%hrv_frootc_to_litter            
    hrv_frootc_storage_to_litter      => pc13f%hrv_frootc_storage_to_litter    
    hrv_frootc_xfer_to_litter         => pc13f%hrv_frootc_xfer_to_litter       
    hrv_livestemc_to_litter           => pc13f%hrv_livestemc_to_litter         
    hrv_livestemc_storage_to_litter   => pc13f%hrv_livestemc_storage_to_litter 
    hrv_livestemc_xfer_to_litter      => pc13f%hrv_livestemc_xfer_to_litter    
    hrv_deadstemc_to_prod10c          => pc13f%hrv_deadstemc_to_prod10c        
    hrv_deadstemc_to_prod100c         => pc13f%hrv_deadstemc_to_prod100c       
    hrv_deadstemc_storage_to_litter   => pc13f%hrv_deadstemc_storage_to_litter 
    hrv_deadstemc_xfer_to_litter      => pc13f%hrv_deadstemc_xfer_to_litter    
    hrv_livecrootc_to_litter          => pc13f%hrv_livecrootc_to_litter        
    hrv_livecrootc_storage_to_litter  => pc13f%hrv_livecrootc_storage_to_litter
    hrv_livecrootc_xfer_to_litter     => pc13f%hrv_livecrootc_xfer_to_litter   
    hrv_deadcrootc_to_litter          => pc13f%hrv_deadcrootc_to_litter        
    hrv_deadcrootc_storage_to_litter  => pc13f%hrv_deadcrootc_storage_to_litter
    hrv_deadcrootc_xfer_to_litter     => pc13f%hrv_deadcrootc_xfer_to_litter   
    hrv_gresp_storage_to_litter       => pc13f%hrv_gresp_storage_to_litter     
    hrv_gresp_xfer_to_litter          => pc13f%hrv_gresp_xfer_to_litter        
    hrv_xsmrpool_to_atm               => pc13f%hrv_xsmrpool_to_atm             
    mr                             => pc13f%mr
    npp                            => pc13f%npp
    pft_fire_closs                 => pc13f%pft_fire_closs
    psnshade_to_cpool              => pc13f%psnshade_to_cpool
    psnsun_to_cpool                => pc13f%psnsun_to_cpool
    rr                             => pc13f%rr
    storage_gr                     => pc13f%storage_gr
    transfer_deadcroot_gr          => pc13f%transfer_deadcroot_gr
    transfer_deadstem_gr           => pc13f%transfer_deadstem_gr
    transfer_froot_gr              => pc13f%transfer_froot_gr
    transfer_gr                    => pc13f%transfer_gr
    transfer_leaf_gr               => pc13f%transfer_leaf_gr
    transfer_livecroot_gr          => pc13f%transfer_livecroot_gr
    transfer_livestem_gr           => pc13f%transfer_livestem_gr
    vegfire                        => pc13f%vegfire
    wood_harvestc                  => pc13f%wood_harvestc
    cpool                          => pc13s%cpool
    xsmrpool                       => pc13s%xsmrpool
    deadcrootc                     => pc13s%deadcrootc
    deadcrootc_storage             => pc13s%deadcrootc_storage
    deadcrootc_xfer                => pc13s%deadcrootc_xfer
    deadstemc                      => pc13s%deadstemc
    deadstemc_storage              => pc13s%deadstemc_storage
    deadstemc_xfer                 => pc13s%deadstemc_xfer
    dispvegc                       => pc13s%dispvegc
    frootc                         => pc13s%frootc
    frootc_storage                 => pc13s%frootc_storage
    frootc_xfer                    => pc13s%frootc_xfer
    gresp_storage                  => pc13s%gresp_storage
    gresp_xfer                     => pc13s%gresp_xfer
    leafc                          => pc13s%leafc
    leafc_storage                  => pc13s%leafc_storage
    leafc_xfer                     => pc13s%leafc_xfer
    livecrootc                     => pc13s%livecrootc
    livecrootc_storage             => pc13s%livecrootc_storage
    livecrootc_xfer                => pc13s%livecrootc_xfer
    livestemc                      => pc13s%livestemc
    livestemc_storage              => pc13s%livestemc_storage
    livestemc_xfer                 => pc13s%livestemc_xfer
    storvegc                       => pc13s%storvegc
    totpftc                        => pc13s%totpftc
    totvegc                        => pc13s%totvegc

   ! pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)

      ! calculate pft-level summary carbon fluxes and states

      ! gross primary production (GPP)
      gpp(p) = &
         psnsun_to_cpool(p) + &
         psnshade_to_cpool(p)

      ! maintenance respiration (MR)
      
      leaf_mr(p)      = leaf_curmr(p)      + leaf_xsmr(p)
      froot_mr(p)     = froot_curmr(p)     + froot_xsmr(p)
      livestem_mr(p)  = livestem_curmr(p)  + livestem_xsmr(p)
      livecroot_mr(p) = livecroot_curmr(p) + livecroot_xsmr(p)
      
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
      ar(p) = mr(p) + gr(p)

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

      ! total vegetation carbon, excluding cpool (TOTVEGC)
      totvegc(p) = dispvegc(p) + storvegc(p)

      ! total pft-level carbon, including cpool (TOTPFTC)
      totpftc(p) = totvegc(p) + xsmrpool(p)

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

      ! net ecosystem production, excludes fire flux, positive for sink (NEP)
      nep(c) = col_gpp(c) - er(c)

      ! net ecosystem exchange of carbon, includes fire flux, positive for source (NBP)
      nbp(c) = nep(c) - col_fire_closs(c) - dwt_closs(c) - product_closs(c)

      ! net ecosystem exchange of carbon, includes fire flux, positive for source (NEE)
      nee(c) = -nep(c) + col_fire_closs(c) + dwt_closs(c) + product_closs(c)

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
      totcolc(c) = &
         cwdc(c) + &
         totlitc(c) + &
         totsomc(c) + &
	      totprodc(c) + &
         col_totpftc(c)

   end do ! end of columns loop


end subroutine C13Summary
!-----------------------------------------------------------------------

end module C13SummaryMod
