module CNSummaryMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for carbon and nitrogen summary calculations
  !
  ! !USES:
  use shr_kind_mod , only: r8 => shr_kind_r8
  use clm_varcon   , only: dzsoi_decomp, zisoi
  use pftvarcon    , only: npcropmin
  use clm_varpar   , only: crop_prog
  use abortutils   , only: endrun
  use decompMod    , only: bounds_type
  use shr_log_mod  , only: errMsg => shr_log_errMsg
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CSummary
  public :: NSummary
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CSummary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, isotope)
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
  use clm_varpar  , only: nlevdecomp,ndecomp_pools,ndecomp_cascade_transitions
  use CNNDynamicsMod, only: nfix_timeconst
  use clm_time_manager    , only : get_step_size
  use clm_varcon      , only: secspday, spval
  !
  ! !ARGUMENTS:
  implicit none
  type(bounds_type), intent(in) :: bounds  ! bounds
  integer, intent(in) :: num_soilc       ! number of soil columns in filter
  integer, intent(in) :: filter_soilc(:) ! filter for soil columns
  integer, intent(in) :: num_soilp       ! number of soil pfts in filter
  integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
  character(len=*), intent(in) :: isotope         ! 'bulk', 'c13' or 'c14'
  !
  ! !LOCAL VARIABLES:
  real(r8), pointer :: woodc(:)               ! wood C (gC/m2)
  real(r8), pointer :: tempsum_npp(:)         ! temporary annual sum of NPP (gC/m2/yr)
  real(r8), pointer :: tempsum_litfall(:)    ! temporary annual sum of litfall (gC/m2/yr)
  real(r8), pointer :: col_lag_npp(:)        ! (gC/m2/s) lagged net primary production
  real(r8), pointer :: som_c_leached(:)                           ! total SOM C loss from vertical transport (gC/m^2/s)
  real(r8), pointer :: decomp_cpools_leached(:,:)                 ! C loss from vertical transport from each decomposing C pool (gC/m^2/s)
  real(r8), pointer :: decomp_cpools_transport_tendency(:,:,:)    ! C tendency due to vertical transport in decomposing C pools (gC/m^3/s)
  real(r8) :: nfixlags, dtime                ! temp variables for making lagged npp
  integer :: c,p,j,k,l        ! indices
  integer :: fp,fc        ! lake filter indices
  real(r8) :: maxdepth    ! depth to integrate soil variables
  type(pft_cflux_type)    , pointer :: pcisof
  type(pft_cstate_type)   , pointer :: pcisos
  type(column_cflux_type) , pointer :: ccisof
  type(column_cstate_type), pointer :: ccisos
  type(pft_cflux_type)    , pointer :: pcisof_a
  type(pft_cstate_type)   , pointer :: pcisos_a
!-----------------------------------------------------------------------

   ! select which isotope
   select case (isotope)
   case ('bulk')
      pcisof =>  pcf
      pcisos =>  pcs
      ccisof =>  ccf
      ccisos =>  ccs
      pcisof_a => pcf_a
      pcisos_a => pcs_a
   case ('c14')
      pcisof =>  pc14f
      pcisos =>  pc14s
      ccisof =>  cc14f
      ccisos =>  cc14s
      pcisof_a => pc14f_a
      pcisos_a => pc14s_a
   case ('c13')
      pcisof =>  pc13f
      pcisos =>  pc13s
      ccisof =>  cc13f
      ccisos =>  cc13s
      pcisof_a => pc13f_a
      pcisos_a => pc13s_a
   case default
      call endrun(msg='CNCIsoSummaryMod: iso must be bulk, c13 or c14'//&
           errMsg(__FILE__, __LINE__))
   end select

   associate(& 
   col_totpftc                     =>    pcisos_a%totpftc                            , & ! Input:  [real(r8) (:)]  (gC/m2) total pft-level carbon , including cpool                       
   col_totvegc                     =>    pcisos_a%totvegc                            , & ! Input:  [real(r8) (:)]  (gC/m2) total vegetation carbon , excluding cpool                      
   ivt                             =>    pft%itype                                   , & ! Input:  [integer (:)]  pft vegetation type                                
   col_fire_closs                  =>    ccisof%col_fire_closs                       , & ! Input:  [real(r8) (:)]  (gC/m2/s) total column-level fire C loss          
   er                              =>    ccisof%er                                   , & ! Input:  [real(r8) (:)]  (gC/m2/s) total ecosystem respiration, autotrophic + heterotrophic
   hr                              =>    ccisof%hr                                   , & ! Input:  [real(r8) (:)]  (gC/m2/s) total heterotrophic respiration         
   litfire                         =>    ccisof%litfire                              , & ! Input:  [real(r8) (:)]  (gC/m2/s) litter fire losses                      
   lithr                           =>    ccisof%lithr                                , & ! Input:  [real(r8) (:)]  (gC/m2/s) litter heterotrophic respiration        
   cwdc                            =>    ccisos%cwdc                                 , & ! Input:  [real(r8) (:)]  (gC/m2) coarse woody debris C                     
   col_ctrunc                      =>    ccisos%col_ctrunc                           , & ! Input:  [real(r8) (:)]  (gC/m2) column-level sink for C truncation        
   decomp_cascade_hr_vr            =>    ccisof%decomp_cascade_hr_vr                 , & ! Input:  [real(r8) (:,:,:)]                                                
   decomp_cascade_hr               =>    ccisof%decomp_cascade_hr                    , & ! Input:  [real(r8) (:,:)]                                                  
   hr_vr                           =>    ccisof%hr_vr                                , & ! Input:  [real(r8) (:,:)]  total vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
   m_decomp_cpools_to_fire_vr      =>    ccisof%m_decomp_cpools_to_fire_vr           , & ! Input:  [real(r8) (:,:,:)]                                                
   m_decomp_cpools_to_fire         =>    ccisof%m_decomp_cpools_to_fire              , & ! Input:  [real(r8) (:,:)]                                                  
   decomp_cascade_ctransfer_vr     =>    ccisof%decomp_cascade_ctransfer_vr          , & ! Input:  [real(r8) (:,:,:)]                                                
   decomp_cascade_ctransfer        =>    ccisof%decomp_cascade_ctransfer             , & ! Input:  [real(r8) (:,:)]                                                  
   decomp_cpools_vr                =>    ccisos%decomp_cpools_vr                     , & ! Input:  [real(r8) (:,:,:)]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
   decomp_cpools                   =>    ccisos%decomp_cpools                        , & ! Input:  [real(r8) (:,:)]  (gC/m2)  decomposing (litter, cwd, soil) c pools
   decomp_cpools_1m                =>    ccisos%decomp_cpools_1m                     , & ! Input:  [real(r8) (:,:)]  (gC/m2)  decomposing (litter, cwd, soil) c pools to 1 meter
   altmax_indx                     =>    cps%altmax_indx                             , & ! Input:  [integer (:)]  maximum annual depth of thaw                       
   altmax_lastyear_indx            =>    cps%altmax_lastyear_indx                    , & ! Input:  [integer (:)]  prior year maximum annual depth of thaw            
   col_ctrunc_vr                   =>    ccisos%col_ctrunc_vr                        , & ! Input:  [real(r8) (:,:)]  (gC/m3) column-level sink for C truncation      
   cascade_donor_pool              =>    decomp_cascade_con%cascade_donor_pool       , & ! Input:  [integer (:)]  which pool is C taken from for a given decomposition step
   is_litter                       =>    decomp_cascade_con%is_litter                , & ! Input:  [logical (:)]  TRUE => pool is a litter pool                      
   is_soil                         =>    decomp_cascade_con%is_soil                  , & ! Input:  [logical (:)]  TRUE => pool is a soil pool                        
   is_cwd                          =>    decomp_cascade_con%is_cwd                   , & ! Input:  [logical (:)]  TRUE => pool is a cwd pool                         
   nee                             =>    ccisof%nee                                  , & ! Input:  [real(r8) (:)]  (gC/m2/s) net ecosystem exchange of carbon, includes fire, land-use, harvest, and hrv_xsmrpool flux, positive for source
   nep                             =>    ccisof%nep                                  , & ! Input:  [real(r8) (:)]  (gC/m2/s) net ecosystem production, excludes fire, land-use, and harvest flux, positive for sink
   nbp                             =>    ccisof%nbp                                  , & ! Input:  [real(r8) (:)]  (gC/m2/s) net biome production, includes fire, land-use, and harvest flux, positive for sink
   col_ar                          =>    pcisof_a%ar                                 , & ! Input:  [real(r8) (:)]  (gC/m2/s) autotrophic respiration (MR + GR)       
   col_gpp                         =>    pcisof_a%gpp                                , & ! Input:  [real(r8) (:)]  GPP flux before downregulation (gC/m2/s)          
   col_npp                         =>    pcisof_a%npp                                , & ! Input:  [real(r8) (:)]  (gC/m2/s) net primary production                  
   col_pft_fire_closs              =>    pcisof_a%pft_fire_closs                     , & ! Input:  [real(r8) (:)]  (gC/m2/s) total pft-level fire C loss             
   col_litfall                     =>    pcisof_a%litfall                            , & ! Input:  [real(r8) (:)]  (gC/m2/s) total pft-level litterfall C loss       
   col_rr                          =>    pcisof_a%rr                                 , & ! Input:  [real(r8) (:)]  (gC/m2/s) root respiration (fine root MR + total root GR)
   col_vegfire                     =>    pcisof_a%vegfire                            , & ! Input:  [real(r8) (:)]  (gC/m2/s) pft-level fire loss (obsolete, mark for removal)
   col_wood_harvestc               =>    pcisof_a%wood_harvestc                      , & ! Input:  [real(r8) (:)]                                                    
   somfire                         =>    ccisof%somfire                              , & ! Input:  [real(r8) (:)]  (gC/m2/s) soil organic matter fire losses         
   somhr                           =>    ccisof%somhr                                , & ! Input:  [real(r8) (:)]  (gC/m2/s) soil organic matter heterotrophic respiration
   sr                              =>    ccisof%sr                                   , & ! Input:  [real(r8) (:)]  (gC/m2/s) total soil respiration (HR + root resp) 
   totfire                         =>    ccisof%totfire                              , & ! Input:  [real(r8) (:)]  (gC/m2/s) total ecosystem fire losses             
   cwdc_hr                         =>    ccisof%cwdc_hr                              , & ! Input:  [real(r8) (:)]  coarse woody debris C heterotrophic respiration (gC/m2/s)
   cwdc_loss                       =>    ccisof%cwdc_loss                            , & ! Input:  [real(r8) (:)]  coarse woody debris C loss (gC/m2/s)              
   litterc_loss                    =>    ccisof%litterc_loss                         , & ! Input:  [real(r8) (:)]  litter C loss (gC/m2/s)                           
   ! dynamic landcover pointers
   dwt_closs                       =>    ccisof%dwt_closs                            , & ! Input:  [real(r8) (:)]  (gC/m2/s) total carbon loss from land cover conversion
   landuseflux                     =>    ccisof%landuseflux                          , & ! Input:  [real(r8) (:)]  (gC/m2/s) dwt_closs+product_closs                 
   landuptake                      =>    ccisof%landuptake                           , & ! Input:  [real(r8) (:)]  (gC/m2/s) nee-landuseflux                         
   dwt_conv_cflux                  =>    ccisof%dwt_conv_cflux                       , & ! Input:  [real(r8) (:)]  (gC/m2/s) conversion C flux (immediate loss to atm)
   seedc                           =>    ccisos%seedc                                , & ! Input:  [real(r8) (:)]  (gC/m2) column-level pool for seeding new PFTs    
   ! wood product pointers
   prod10c_loss                    =>    ccisof%prod10c_loss                         , & ! Input:  [real(r8) (:)]  (gC/m2/s) loss from 10-yr wood product pool       
   prod100c_loss                   =>    ccisof%prod100c_loss                        , & ! Input:  [real(r8) (:)]  (gC/m2/s) loss from 100-yr wood product pool      
   product_closs                   =>    ccisof%product_closs                        , & ! Input:  [real(r8) (:)]  (gC/m2/s) total wood product carbon loss          
   prod10c                         =>    ccisos%prod10c                              , & ! Input:  [real(r8) (:)]  (gC/m2) wood product C pool, 10-year lifespan     
   prod100c                        =>    ccisos%prod100c                             , & ! Input:  [real(r8) (:)]  (gC/m2) wood product C pool, 100-year lifespan    
   totprodc                        =>    ccisos%totprodc                             , & ! Input:  [real(r8) (:)]  (gC/m2) total wood product C                      
   totcolc                         =>    ccisos%totcolc                              , & ! Input:  [real(r8) (:)]  (gC/m2) total column carbon, incl veg and cpool   
   totecosysc                      =>    ccisos%totecosysc                           , & ! Input:  [real(r8) (:)]  (gC/m2) total ecosystem carbon, incl veg but excl cpool
   totlitc                         =>    ccisos%totlitc                              , & ! Input:  [real(r8) (:)]  (gC/m2) total litter carbon                       
   totlitc_1m                      =>    ccisos%totlitc_1m                           , & ! Input:  [real(r8) (:)]  (gC/m2) total litter carbon to 1 meter            
   totsomc                         =>    ccisos%totsomc                              , & ! Input:  [real(r8) (:)]  (gC/m2) total soil organic matter carbon          
   totsomc_1m                      =>    ccisos%totsomc_1m                           , & ! Input:  [real(r8) (:)]  (gC/m2) total soil organic matter carbon to 1 meter
   agnpp                           =>    pcisof%agnpp                                , & ! Input:  [real(r8) (:)]  (gC/m2/s) aboveground NPP                         
   ar                              =>    pcisof%ar                                   , & ! Input:  [real(r8) (:)]  (gC/m2/s) autotrophic respiration (MR + GR)       
   bgnpp                           =>    pcisof%bgnpp                                , & ! Input:  [real(r8) (:)]  (gC/m2/s) belowground NPP                         
   xsmrpool_to_atm                 =>    pcisof%xsmrpool_to_atm                      , & ! Input:  [real(r8) (:)]  excess MR pool harvest mortality (gC/m2/s)        
   cpool_grain_gr                  =>    pcisof%cpool_grain_gr                       , & ! Input:  [real(r8) (:)]  grain growth respiration (gC/m2/s)                
   cpool_grain_storage_gr          =>    pcisof%cpool_grain_storage_gr               , & ! Input:  [real(r8) (:)]  grain growth respiration to storage (gC/m2/s)     
   cpool_to_grainc                 =>    pcisof%cpool_to_grainc                      , & ! Input:  [real(r8) (:)]  allocation to grain C storage (gC/m2/s)           
   grainc_xfer_to_grainc           =>    pcisof%grainc_xfer_to_grainc                , & ! Input:  [real(r8) (:)]  grain C growth from storage (gC/m2/s)             
   transfer_grain_gr               =>    pcisof%transfer_grain_gr                    , & ! Input:  [real(r8) (:)]  grain growth respiration from storage (gC/m2/s)   
   grainc_to_food                  =>    pcisof%grainc_to_food                       , & ! Input:  [real(r8) (:)]  grain C to food (gC/m2/s)                         
   livestemc_to_litter             =>    pcisof%livestemc_to_litter                  , & ! Input:  [real(r8) (:)]  live stem C litterfall (gC/m2/s)                  
   grainc                          =>    pcisos%grainc                               , & ! Input:  [real(r8) (:)]  (gC/m2) grain C                                   
   grainc_storage                  =>    pcisos%grainc_storage                       , & ! Input:  [real(r8) (:)]  (gC/m2) grain C storage                           
   grainc_xfer                     =>    pcisos%grainc_xfer                          , & ! Input:  [real(r8) (:)]  (gC/m2) grain C transfer                          
   cpool_deadcroot_gr              =>    pcisof%cpool_deadcroot_gr                   , & ! Input:  [real(r8) (:)]  dead coarse root growth respiration (gC/m2/s)     
   cpool_deadcroot_storage_gr      =>    pcisof%cpool_deadcroot_storage_gr           , & ! Input:  [real(r8) (:)]  dead coarse root growth respiration to storage (gC/m2/s)
   cpool_deadstem_gr               =>    pcisof%cpool_deadstem_gr                    , & ! Input:  [real(r8) (:)]  dead stem growth respiration (gC/m2/s)            
   cpool_deadstem_storage_gr       =>    pcisof%cpool_deadstem_storage_gr            , & ! Input:  [real(r8) (:)]  dead stem growth respiration to storage (gC/m2/s) 
   cpool_froot_gr                  =>    pcisof%cpool_froot_gr                       , & ! Input:  [real(r8) (:)]  fine root growth respiration (gC/m2/s)            
   cpool_froot_storage_gr          =>    pcisof%cpool_froot_storage_gr               , & ! Input:  [real(r8) (:)]  fine root  growth respiration to storage (gC/m2/s)
   cpool_leaf_gr                   =>    pcisof%cpool_leaf_gr                        , & ! Input:  [real(r8) (:)]  leaf growth respiration (gC/m2/s)                 
   cpool_leaf_storage_gr           =>    pcisof%cpool_leaf_storage_gr                , & ! Input:  [real(r8) (:)]  leaf growth respiration to storage (gC/m2/s)      
   cpool_livecroot_gr              =>    pcisof%cpool_livecroot_gr                   , & ! Input:  [real(r8) (:)]  live coarse root growth respiration (gC/m2/s)     
   cpool_livecroot_storage_gr      =>    pcisof%cpool_livecroot_storage_gr           , & ! Input:  [real(r8) (:)]  live coarse root growth respiration to storage (gC/m2/s)
   cpool_livestem_gr               =>    pcisof%cpool_livestem_gr                    , & ! Input:  [real(r8) (:)]  live stem growth respiration (gC/m2/s)            
   cpool_livestem_storage_gr       =>    pcisof%cpool_livestem_storage_gr            , & ! Input:  [real(r8) (:)]  live stem growth respiration to storage (gC/m2/s) 
   cpool_to_deadcrootc             =>    pcisof%cpool_to_deadcrootc                  , & ! Input:  [real(r8) (:)]  allocation to dead coarse root C (gC/m2/s)        
   cpool_to_deadstemc              =>    pcisof%cpool_to_deadstemc                   , & ! Input:  [real(r8) (:)]  allocation to dead stem C (gC/m2/s)               
   cpool_to_frootc                 =>    pcisof%cpool_to_frootc                      , & ! Input:  [real(r8) (:)]  allocation to fine root C (gC/m2/s)               
   cpool_to_leafc                  =>    pcisof%cpool_to_leafc                       , & ! Input:  [real(r8) (:)]  allocation to leaf C (gC/m2/s)                    
   cpool_to_livecrootc             =>    pcisof%cpool_to_livecrootc                  , & ! Input:  [real(r8) (:)]  allocation to live coarse root C (gC/m2/s)        
   cpool_to_livestemc              =>    pcisof%cpool_to_livestemc                   , & ! Input:  [real(r8) (:)]  allocation to live stem C (gC/m2/s)               
   current_gr                      =>    pcisof%current_gr                           , & ! Input:  [real(r8) (:)]  (gC/m2/s) growth resp for new growth displayed in this timestep
   deadcrootc_xfer_to_deadcrootc   =>    pcisof%deadcrootc_xfer_to_deadcrootc        , & ! Input:  [real(r8) (:)]                                                    
   deadstemc_xfer_to_deadstemc     =>    pcisof%deadstemc_xfer_to_deadstemc          , & ! Input:  [real(r8) (:)]                                                    
   frootc_to_litter                =>    pcisof%frootc_to_litter                     , & ! Input:  [real(r8) (:)]                                                    
   frootc_xfer_to_frootc           =>    pcisof%frootc_xfer_to_frootc                , & ! Input:  [real(r8) (:)]                                                    
   froot_mr                        =>    pcisof%froot_mr                             , & ! Input:  [real(r8) (:)]                                                    
   froot_curmr                     =>    pcisof%froot_curmr                          , & ! Input:  [real(r8) (:)]                                                    
   froot_xsmr                      =>    pcisof%froot_xsmr                           , & ! Input:  [real(r8) (:)]                                                    
   grain_mr                        =>    pcisof%grain_mr                             , & ! Input:  [real(r8) (:)]                                                    
   gpp                             =>    pcisof%gpp                                  , & ! Input:  [real(r8) (:)] GPP flux before downregulation (gC/m2/s)           
   gr                              =>    pcisof%gr                                   , & ! Input:  [real(r8) (:)]  (gC/m2/s) total growth respiration                
   leafc_to_litter                 =>    pcisof%leafc_to_litter                      , & ! Input:  [real(r8) (:)]                                                    
   leafc_xfer_to_leafc             =>    pcisof%leafc_xfer_to_leafc                  , & ! Input:  [real(r8) (:)]                                                    
   leaf_mr                         =>    pcisof%leaf_mr                              , & ! Input:  [real(r8) (:)]                                                    
   leaf_curmr                      =>    pcisof%leaf_curmr                           , & ! Input:  [real(r8) (:)]                                                    
   leaf_xsmr                       =>    pcisof%leaf_xsmr                            , & ! Input:  [real(r8) (:)]                                                    
   litfall                         =>    pcisof%litfall                              , & ! Input:  [real(r8) (:)]  (gC/m2/s) litterfall (leaves and fine roots)      
   livecrootc_xfer_to_livecrootc   =>    pcisof%livecrootc_xfer_to_livecrootc        , & ! Input:  [real(r8) (:)]                                                    
   livecroot_mr                    =>    pcisof%livecroot_mr                         , & ! Input:  [real(r8) (:)]                                                    
   livecroot_curmr                 =>    pcisof%livecroot_curmr                      , & ! Input:  [real(r8) (:)]                                                    
   livecroot_xsmr                  =>    pcisof%livecroot_xsmr                       , & ! Input:  [real(r8) (:)]                                                    
   livestemc_xfer_to_livestemc     =>    pcisof%livestemc_xfer_to_livestemc          , & ! Input:  [real(r8) (:)]                                                    
   livestem_mr                     =>    pcisof%livestem_mr                          , & ! Input:  [real(r8) (:)]                                                    
   livestem_curmr                  =>    pcisof%livestem_curmr                       , & ! Input:  [real(r8) (:)]                                                    
   livestem_xsmr                   =>    pcisof%livestem_xsmr                        , & ! Input:  [real(r8) (:)]                                                    
   m_leafc_to_fire                 =>    pcisof%m_leafc_to_fire                      , & ! Input:  [real(r8) (:)]                                                    
   m_leafc_storage_to_fire         =>    pcisof%m_leafc_storage_to_fire              , & ! Input:  [real(r8) (:)]                                                    
   m_leafc_xfer_to_fire            =>    pcisof%m_leafc_xfer_to_fire                 , & ! Input:  [real(r8) (:)]                                                    
   m_livestemc_to_fire             =>    pcisof%m_livestemc_to_fire                  , & ! Input:  [real(r8) (:)]                                                    
   m_livestemc_storage_to_fire     =>    pcisof%m_livestemc_storage_to_fire          , & ! Input:  [real(r8) (:)]                                                    
   m_livestemc_xfer_to_fire        =>    pcisof%m_livestemc_xfer_to_fire             , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemc_to_fire             =>    pcisof%m_deadstemc_to_fire                  , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemc_storage_to_fire     =>    pcisof%m_deadstemc_storage_to_fire          , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemc_xfer_to_fire        =>    pcisof%m_deadstemc_xfer_to_fire             , & ! Input:  [real(r8) (:)]                                                    
   m_frootc_to_fire                =>    pcisof%m_frootc_to_fire                     , & ! Input:  [real(r8) (:)]                                                    
   m_frootc_storage_to_fire        =>    pcisof%m_frootc_storage_to_fire             , & ! Input:  [real(r8) (:)]                                                    
   m_frootc_xfer_to_fire           =>    pcisof%m_frootc_xfer_to_fire                , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootc_to_fire            =>    pcisof%m_livecrootc_to_fire                 , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootc_storage_to_fire    =>    pcisof%m_livecrootc_storage_to_fire         , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootc_xfer_to_fire       =>    pcisof%m_livecrootc_xfer_to_fire            , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootc_to_fire            =>    pcisof%m_deadcrootc_to_fire                 , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootc_storage_to_fire    =>    pcisof%m_deadcrootc_storage_to_fire         , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootc_xfer_to_fire       =>    pcisof%m_deadcrootc_xfer_to_fire            , & ! Input:  [real(r8) (:)]                                                    
   m_gresp_storage_to_fire         =>    pcisof%m_gresp_storage_to_fire              , & ! Input:  [real(r8) (:)]                                                    
   m_gresp_xfer_to_fire            =>    pcisof%m_gresp_xfer_to_fire                 , & ! Input:  [real(r8) (:)]                                                    
   m_leafc_to_litter_fire          =>    pcisof%m_leafc_to_litter_fire               , & ! Input:  [real(r8) (:)]                                                    
   m_leafc_storage_to_litter_fire  =>    pcisof%m_leafc_storage_to_litter_fire       , & ! Input:  [real(r8) (:)]                                                    
   m_leafc_xfer_to_litter_fire     =>    pcisof%m_leafc_xfer_to_litter_fire          , & ! Input:  [real(r8) (:)]                                                    
   m_livestemc_to_litter_fire      =>    pcisof%m_livestemc_to_litter_fire           , & ! Input:  [real(r8) (:)]                                                    
   m_livestemc_storage_to_litter_fire  =>    pcisof%m_livestemc_storage_to_litter_fire, & ! Input:  [real(r8) (:)]                                                    
   m_livestemc_xfer_to_litter_fire =>    pcisof%m_livestemc_xfer_to_litter_fire      , & ! Input:  [real(r8) (:)]                                                    
   m_livestemc_to_deadstemc_fire   =>    pcisof%m_livestemc_to_deadstemc_fire        , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemc_to_litter_fire      =>    pcisof%m_deadstemc_to_litter_fire           , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemc_storage_to_litter_fire  =>    pcisof%m_deadstemc_storage_to_litter_fire, & ! Input:  [real(r8) (:)]                                                    
   m_deadstemc_xfer_to_litter_fire =>    pcisof%m_deadstemc_xfer_to_litter_fire      , & ! Input:  [real(r8) (:)]                                                    
   m_frootc_to_litter_fire         =>    pcisof%m_frootc_to_litter_fire              , & ! Input:  [real(r8) (:)]                                                    
   m_frootc_storage_to_litter_fire =>    pcisof%m_frootc_storage_to_litter_fire      , & ! Input:  [real(r8) (:)]                                                    
   m_frootc_xfer_to_litter_fire    =>    pcisof%m_frootc_xfer_to_litter_fire         , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootc_to_litter_fire     =>    pcisof%m_livecrootc_to_litter_fire          , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootc_storage_to_litter_fire  =>    pcisof%m_livecrootc_storage_to_litter_fire, & ! Input:  [real(r8) (:)]                                                    
   m_livecrootc_xfer_to_litter_fire=>    pcisof%m_livecrootc_xfer_to_litter_fire     , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootc_to_deadcrootc_fire =>    pcisof%m_livecrootc_to_deadcrootc_fire      , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootc_to_litter_fire     =>    pcisof%m_deadcrootc_to_litter_fire          , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootc_storage_to_litter_fire  =>    pcisof%m_deadcrootc_storage_to_litter_fire, & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootc_xfer_to_litter_fire=>    pcisof%m_deadcrootc_xfer_to_litter_fire     , & ! Input:  [real(r8) (:)]                                                    
   m_gresp_storage_to_litter_fire  =>    pcisof%m_gresp_storage_to_litter_fire       , & ! Input:  [real(r8) (:)]                                                    
   m_gresp_xfer_to_litter_fire     =>    pcisof%m_gresp_xfer_to_litter_fire          , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootc_storage_to_litter  =>    pcisof%m_deadcrootc_storage_to_litter       , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootc_to_litter          =>    pcisof%m_deadcrootc_to_litter               , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootc_xfer_to_litter     =>    pcisof%m_deadcrootc_xfer_to_litter          , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemc_storage_to_litter   =>    pcisof%m_deadstemc_storage_to_litter        , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemc_to_litter           =>    pcisof%m_deadstemc_to_litter                , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemc_xfer_to_litter      =>    pcisof%m_deadstemc_xfer_to_litter           , & ! Input:  [real(r8) (:)]                                                    
   m_frootc_storage_to_litter      =>    pcisof%m_frootc_storage_to_litter           , & ! Input:  [real(r8) (:)]                                                    
   m_frootc_to_litter              =>    pcisof%m_frootc_to_litter                   , & ! Input:  [real(r8) (:)]                                                    
   m_frootc_xfer_to_litter         =>    pcisof%m_frootc_xfer_to_litter              , & ! Input:  [real(r8) (:)]                                                    
   m_gresp_storage_to_litter       =>    pcisof%m_gresp_storage_to_litter            , & ! Input:  [real(r8) (:)]                                                    
   m_gresp_xfer_to_litter          =>    pcisof%m_gresp_xfer_to_litter               , & ! Input:  [real(r8) (:)]                                                    
   m_leafc_storage_to_litter       =>    pcisof%m_leafc_storage_to_litter            , & ! Input:  [real(r8) (:)]                                                    
   m_leafc_to_litter               =>    pcisof%m_leafc_to_litter                    , & ! Input:  [real(r8) (:)]                                                    
   m_leafc_xfer_to_litter          =>    pcisof%m_leafc_xfer_to_litter               , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootc_storage_to_litter  =>    pcisof%m_livecrootc_storage_to_litter       , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootc_to_litter          =>    pcisof%m_livecrootc_to_litter               , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootc_xfer_to_litter     =>    pcisof%m_livecrootc_xfer_to_litter          , & ! Input:  [real(r8) (:)]                                                    
   m_livestemc_storage_to_litter   =>    pcisof%m_livestemc_storage_to_litter        , & ! Input:  [real(r8) (:)]                                                    
   m_livestemc_to_litter           =>    pcisof%m_livestemc_to_litter                , & ! Input:  [real(r8) (:)]                                                    
   m_livestemc_xfer_to_litter      =>    pcisof%m_livestemc_xfer_to_litter           , & ! Input:  [real(r8) (:)]                                                    
   hrv_leafc_to_litter             =>    pcisof%hrv_leafc_to_litter                  , & ! Input:  [real(r8) (:)]                                                    
   hrv_leafc_storage_to_litter     =>    pcisof%hrv_leafc_storage_to_litter          , & ! Input:  [real(r8) (:)]                                                    
   hrv_leafc_xfer_to_litter        =>    pcisof%hrv_leafc_xfer_to_litter             , & ! Input:  [real(r8) (:)]                                                    
   hrv_frootc_to_litter            =>    pcisof%hrv_frootc_to_litter                 , & ! Input:  [real(r8) (:)]                                                    
   hrv_frootc_storage_to_litter    =>    pcisof%hrv_frootc_storage_to_litter         , & ! Input:  [real(r8) (:)]                                                    
   hrv_frootc_xfer_to_litter       =>    pcisof%hrv_frootc_xfer_to_litter            , & ! Input:  [real(r8) (:)]                                                    
   hrv_livestemc_to_litter         =>    pcisof%hrv_livestemc_to_litter              , & ! Input:  [real(r8) (:)]                                                    
   hrv_livestemc_storage_to_litter =>    pcisof%hrv_livestemc_storage_to_litter      , & ! Input:  [real(r8) (:)]                                                    
   hrv_livestemc_xfer_to_litter    =>    pcisof%hrv_livestemc_xfer_to_litter         , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadstemc_to_prod10c        =>    pcisof%hrv_deadstemc_to_prod10c             , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadstemc_to_prod100c       =>    pcisof%hrv_deadstemc_to_prod100c            , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadstemc_storage_to_litter =>    pcisof%hrv_deadstemc_storage_to_litter      , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadstemc_xfer_to_litter    =>    pcisof%hrv_deadstemc_xfer_to_litter         , & ! Input:  [real(r8) (:)]                                                    
   hrv_livecrootc_to_litter        =>    pcisof%hrv_livecrootc_to_litter             , & ! Input:  [real(r8) (:)]                                                    
   hrv_livecrootc_storage_to_litter=>    pcisof%hrv_livecrootc_storage_to_litter     , & ! Input:  [real(r8) (:)]                                                    
   hrv_livecrootc_xfer_to_litter   =>    pcisof%hrv_livecrootc_xfer_to_litter        , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadcrootc_to_litter        =>    pcisof%hrv_deadcrootc_to_litter             , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadcrootc_storage_to_litter=>    pcisof%hrv_deadcrootc_storage_to_litter     , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadcrootc_xfer_to_litter   =>    pcisof%hrv_deadcrootc_xfer_to_litter        , & ! Input:  [real(r8) (:)]                                                    
   hrv_gresp_storage_to_litter     =>    pcisof%hrv_gresp_storage_to_litter          , & ! Input:  [real(r8) (:)]                                                    
   hrv_gresp_xfer_to_litter        =>    pcisof%hrv_gresp_xfer_to_litter             , & ! Input:  [real(r8) (:)]                                                    
   hrv_xsmrpool_to_atm             =>    pcisof%hrv_xsmrpool_to_atm                  , & ! Input:  [real(r8) (:)]                                                    
   col_hrv_xsmrpool_to_atm         =>    pcisof_a%hrv_xsmrpool_to_atm                , & ! Input:  [real(r8) (:)]                                                    
   mr                              =>    pcisof%mr                                   , & ! Input:  [real(r8) (:)]  (gC/m2/s) maintenance respiration                 
   npp                             =>    pcisof%npp                                  , & ! Input:  [real(r8) (:)]  (gC/m2/s) net primary production                  
   pft_fire_closs                  =>    pcisof%pft_fire_closs                       , & ! Input:  [real(r8) (:)]  (gC/m2/s) total pft-level fire C loss             
   psnshade_to_cpool               =>    pcisof%psnshade_to_cpool                    , & ! Input:  [real(r8) (:)]                                                    
   psnsun_to_cpool                 =>    pcisof%psnsun_to_cpool                      , & ! Input:  [real(r8) (:)]                                                    
   rr                              =>    pcisof%rr                                   , & ! Input:  [real(r8) (:)]  (gC/m2/s) root respiration (fine root MR + total root GR)
   storage_gr                      =>    pcisof%storage_gr                           , & ! Input:  [real(r8) (:)]  (gC/m2/s) growth resp for growth sent to storage for later display
   transfer_deadcroot_gr           =>    pcisof%transfer_deadcroot_gr                , & ! Input:  [real(r8) (:)]                                                    
   transfer_deadstem_gr            =>    pcisof%transfer_deadstem_gr                 , & ! Input:  [real(r8) (:)]                                                    
   transfer_froot_gr               =>    pcisof%transfer_froot_gr                    , & ! Input:  [real(r8) (:)]                                                    
   transfer_gr                     =>    pcisof%transfer_gr                          , & ! Input:  [real(r8) (:)]  (gC/m2/s) growth resp for transfer growth displayed in this timestep
   transfer_leaf_gr                =>    pcisof%transfer_leaf_gr                     , & ! Input:  [real(r8) (:)]                                                    
   transfer_livecroot_gr           =>    pcisof%transfer_livecroot_gr                , & ! Input:  [real(r8) (:)]                                                    
   transfer_livestem_gr            =>    pcisof%transfer_livestem_gr                 , & ! Input:  [real(r8) (:)]                                                    
   vegfire                         =>    pcisof%vegfire                              , & ! Input:  [real(r8) (:)]  (gC/m2/s) pft-level fire loss (obsolete, mark for removal)
   wood_harvestc                   =>    pcisof%wood_harvestc                        , & ! Input:  [real(r8) (:)]  (gC/m2/s) pft-level wood harvest (to product pools)
   frootc_alloc                    =>    pcisof%frootc_alloc                         , & ! Input:  [real(r8) (:)]  fine root C allocation (gC/m2/s)                  
   frootc_loss                     =>    pcisof%frootc_loss                          , & ! Input:  [real(r8) (:)]  fine root C loss (gC/m2/s)                        
   leafc_alloc                     =>    pcisof%leafc_alloc                          , & ! Input:  [real(r8) (:)]  leaf C allocation (gC/m2/s)                       
   leafc_loss                      =>    pcisof%leafc_loss                           , & ! Input:  [real(r8) (:)]  leaf C loss (gC/m2/s)                             
   woodc_alloc                     =>    pcisof%woodc_alloc                          , & ! Input:  [real(r8) (:)]  wood C allocation (gC/m2/s)                       
   woodc_loss                      =>    pcisof%woodc_loss                           , & ! Input:  [real(r8) (:)]  wood C loss (gC/m2/s)                             
   cpool                           =>    pcisos%cpool                                , & ! Input:  [real(r8) (:)]  (gC/m2) temporary photosynthate C pool            
   xsmrpool                        =>    pcisos%xsmrpool                             , & ! Input:  [real(r8) (:)]  (gC/m2) temporary photosynthate C pool            
   pft_ctrunc                      =>    pcisos%pft_ctrunc                           , & ! Input:  [real(r8) (:)]  (gC/m2) pft-level sink for C truncation           
   deadcrootc                      =>    pcisos%deadcrootc                           , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C                        
   deadcrootc_storage              =>    pcisos%deadcrootc_storage                   , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C storage                
   deadcrootc_xfer                 =>    pcisos%deadcrootc_xfer                      , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C transfer               
   deadstemc                       =>    pcisos%deadstemc                            , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C                               
   deadstemc_storage               =>    pcisos%deadstemc_storage                    , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C storage                       
   deadstemc_xfer                  =>    pcisos%deadstemc_xfer                       , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C transfer                      
   dispvegc                        =>    pcisos%dispvegc                             , & ! Input:  [real(r8) (:)]  (gC/m2) displayed veg carbon, excluding storage and cpool
   frootc                          =>    pcisos%frootc                               , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C                               
   frootc_storage                  =>    pcisos%frootc_storage                       , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C storage                       
   frootc_xfer                     =>    pcisos%frootc_xfer                          , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C transfer                      
   gresp_storage                   =>    pcisos%gresp_storage                        , & ! Input:  [real(r8) (:)]  (gC/m2) growth respiration storage                
   gresp_xfer                      =>    pcisos%gresp_xfer                           , & ! Input:  [real(r8) (:)]  (gC/m2) growth respiration transfer               
   leafc                           =>    pcisos%leafc                                , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C                                    
   leafc_storage                   =>    pcisos%leafc_storage                        , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C storage                            
   leafc_xfer                      =>    pcisos%leafc_xfer                           , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C transfer                           
   livecrootc                      =>    pcisos%livecrootc                           , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C                        
   livecrootc_storage              =>    pcisos%livecrootc_storage                   , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C storage                
   livecrootc_xfer                 =>    pcisos%livecrootc_xfer                      , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C transfer               
   livestemc                       =>    pcisos%livestemc                            , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C                               
   livestemc_storage               =>    pcisos%livestemc_storage                    , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C storage                       
   livestemc_xfer                  =>    pcisos%livestemc_xfer                       , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C transfer                      
   storvegc                        =>    pcisos%storvegc                             , & ! Input:  [real(r8) (:)]  (gC/m2) stored vegetation carbon, excluding cpool 
   totpftc                         =>    pcisos%totpftc                              , & ! Input:  [real(r8) (:)]  (gC/m2) total pft-level carbon, including cpool   
   totvegc                         =>    pcisos%totvegc &
   )

   woodc                           =>    pcisos%woodc                                 ! Input:  [real(r8) (:)]  wood C (gC/m2)                                    
   tempsum_npp                     =>    pepv%tempsum_npp                             ! Input:  [real(r8) (:)]  temporary annual sum of NPP (gC/m2/yr)            
   tempsum_litfall                 =>    pepv%tempsum_litfall                         ! Input:  [real(r8) (:)] temporary annual sum of litfall (gC/m2/yr)         
   col_lag_npp                     =>    cps%col_lag_npp                              ! Input:  [real(r8) (:)]  (gC/m2/s) lagged net primary production           
   som_c_leached                   =>    ccisof%som_c_leached                         ! Input:  [real(r8) (:)]  total SOM C loss from vertical transport (gC/m^2/s)
   decomp_cpools_leached           =>    ccisof%decomp_cpools_leached                 ! Input:  [real(r8) (:,:)]  C loss from vertical transport from each decomposing C pool (gC/m^2/s)
   decomp_cpools_transport_tendency=>    ccisof%decomp_cpools_transport_tendency      ! Input:  [real(r8) (:,:,:)]  C tendency due to vertical transport in decomposing C pools (gC/m^3/s)

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
   call p2c(bounds, num_soilc, filter_soilc, &
        gpp(bounds%begp:bounds%endp), &
        col_gpp(bounds%begc:bounds%endc))
   call p2c(bounds, num_soilc, filter_soilc, &
        ar(bounds%begp:bounds%endp), &
        col_ar(bounds%begc:bounds%endc))
   call p2c(bounds, num_soilc, filter_soilc, &
        rr(bounds%begp:bounds%endp), &
        col_rr(bounds%begc:bounds%endc))
   call p2c(bounds, num_soilc, filter_soilc, &
        npp(bounds%begp:bounds%endp), &
        col_npp(bounds%begc:bounds%endc))
   call p2c(bounds, num_soilc, filter_soilc, &
        vegfire(bounds%begp:bounds%endp), &
        col_vegfire(bounds%begc:bounds%endc))
   call p2c(bounds, num_soilc, filter_soilc, &
        wood_harvestc(bounds%begp:bounds%endp), &
        col_wood_harvestc(bounds%begc:bounds%endc))
   call p2c(bounds, num_soilc, filter_soilc, &
        totvegc(bounds%begp:bounds%endp), &
        col_totvegc(bounds%begc:bounds%endc))
   call p2c(bounds, num_soilc, filter_soilc, &
        totpftc(bounds%begp:bounds%endp), &
        col_totpftc(bounds%begc:bounds%endc))
   call p2c(bounds, num_soilc, filter_soilc, &
        pft_fire_closs(bounds%begp:bounds%endp), &
        col_pft_fire_closs(bounds%begc:bounds%endc))
   call p2c(bounds, num_soilc, filter_soilc, &
        litfall(bounds%begp:bounds%endp), &
        col_litfall(bounds%begc:bounds%endc))
   call p2c(bounds, num_soilc, filter_soilc, &
        hrv_xsmrpool_to_atm(bounds%begp:bounds%endp), &
        col_hrv_xsmrpool_to_atm(bounds%begc:bounds%endc))

   if ( isotope .eq. 'bulk') then
      if (nfix_timeconst .gt. 0._r8 .and. nfix_timeconst .lt. 500._r8 ) then
         ! this code is to calculate an exponentially-relaxed npp value for use in NDynamics code
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

    end associate 
 end subroutine CSummary

 !-----------------------------------------------------------------------
 subroutine NSummary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
   !
   ! !DESCRIPTION:
   ! On the radiation time step, perform pft and column-level nitrogen
   ! summary calculations
   !
   ! !USES:
   use clmtype
   use pft2colMod, only: p2c
   use clm_varpar, only: nlevdecomp,ndecomp_cascade_transitions,ndecomp_pools
   use clm_varctl, only: use_nitrif_denitrif
   !
   ! !ARGUMENTS:
   implicit none
   type(bounds_type), intent(in) :: bounds  ! bounds
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
   !
   ! !LOCAL VARIABLES:
   integer :: c,p,j,k,l       ! indices
   integer :: fp,fc       ! lake filter indices
   real(r8) :: maxdepth    ! depth to integrate soil variables
   !-----------------------------------------------------------------------

   associate(& 
   ivt                             =>   pft%itype                                    , & ! Input:  [integer (:)]  pft vegetation type                                
   col_fire_nloss                  =>    cnf%col_fire_nloss                          , & ! Input:  [real(r8) (:)]  (gN/m2/s) total column-level fire N loss          
   denit                           =>    cnf%denit                                   , & ! Input:  [real(r8) (:)]                                                    
   col_pft_fire_nloss              =>    pnf_a%pft_fire_nloss                        , & ! Input:  [real(r8) (:)]  (gN/m2/s) total pft-level fire C loss             
   cwdn                            =>    cns%cwdn                                    , & ! Input:  [real(r8) (:)]  (gN/m2) coarse woody debris N                     
   col_ntrunc                      =>    cns%col_ntrunc                              , & ! Input:  [real(r8) (:)]  (gN/m2) column-level sink for N truncation        
   sminn                           =>    cns%sminn                                   , & ! Input:  [real(r8) (:)]  (gN/m2) soil mineral N                            
   m_decomp_npools_to_fire_vr      =>    cnf%m_decomp_npools_to_fire_vr              , & ! Input:  [real(r8) (:,:,:)]                                                
   m_decomp_npools_to_fire         =>    cnf%m_decomp_npools_to_fire                 , & ! Input:  [real(r8) (:,:)]                                                  
   is_litter                       =>    decomp_cascade_con%is_litter                , & ! Input:  [logical (:)]  TRUE => pool is a litter pool                      
   is_soil                         =>    decomp_cascade_con%is_soil                  , & ! Input:  [logical (:)]  TRUE => pool is a soil pool                        
   is_cwd                          =>    decomp_cascade_con%is_cwd                   , & ! Input:  [logical (:)]  TRUE => pool is a cwd pool                         
   sminn_to_denit_excess_vr        =>    cnf%sminn_to_denit_excess_vr                , & ! Input:  [real(r8) (:,:)]                                                  
   sminn_to_denit_excess           =>    cnf%sminn_to_denit_excess                   , & ! Input:  [real(r8) (:)]                                                    
   sminn_to_denit_decomp_cascade_vr=>    cnf%sminn_to_denit_decomp_cascade_vr        , & ! Input:  [real(r8) (:,:,:)]  vertically-resolved denitrification along decomp cascade (gN/m3/s)
   sminn_to_denit_decomp_cascade   =>    cnf%sminn_to_denit_decomp_cascade           , & ! Input:  [real(r8) (:,:)]  vertically-integrated denitrification along decomp cascade (gN/m2/s)
   sminn_leached_vr                =>    cnf%sminn_leached_vr                        , & ! Input:  [real(r8) (:,:)]                                                  
   sminn_leached                   =>    cnf%sminn_leached                           , & ! Input:  [real(r8) (:)]                                                    
   smin_no3                        =>    cns%smin_no3                                , & ! Input:  [real(r8) (:)]                                                    
   smin_nh4                        =>    cns%smin_nh4                                , & ! Input:  [real(r8) (:)]                                                    
   smin_no3_vr                     =>    cns%smin_no3_vr                             , & ! Input:  [real(r8) (:,:)]                                                  
   smin_nh4_vr                     =>    cns%smin_nh4_vr                             , & ! Input:  [real(r8) (:,:)]                                                  
   f_nit_vr                        =>    cnf%f_nit_vr                                , & ! Input:  [real(r8) (:,:)]                                                  
   f_nit                           =>    cnf%f_nit                                   , & ! Input:  [real(r8) (:)]                                                    
   f_denit_vr                      =>    cnf%f_denit_vr                              , & ! Input:  [real(r8) (:,:)]                                                  
   f_denit                         =>    cnf%f_denit                                 , & ! Input:  [real(r8) (:)]                                                    
   pot_f_nit_vr                    =>    cnf%pot_f_nit_vr                            , & ! Input:  [real(r8) (:,:)]                                                  
   pot_f_nit                       =>    cnf%pot_f_nit                               , & ! Input:  [real(r8) (:)]                                                    
   pot_f_denit_vr                  =>    cnf%pot_f_denit_vr                          , & ! Input:  [real(r8) (:,:)]                                                  
   pot_f_denit                     =>    cnf%pot_f_denit                             , & ! Input:  [real(r8) (:)]                                                    
   f_n2o_denit_vr                  =>    cnf%f_n2o_denit_vr                          , & ! Input:  [real(r8) (:,:)]  flux of N2o from denitrification [gN/m3/s]      
   f_n2o_nit_vr                    =>    cnf%f_n2o_nit_vr                            , & ! Input:  [real(r8) (:,:)]  flux of N2o from nitrification [gN/m3/s]        
   f_n2o_denit                     =>    cnf%f_n2o_denit                             , & ! Input:  [real(r8) (:)]  flux of N2o from denitrification [gN/m2/s]        
   f_n2o_nit                       =>    cnf%f_n2o_nit                               , & ! Input:  [real(r8) (:)]  flux of N2o from nitrification [gN/m2/s]          
   smin_no3_leached_vr             =>    cnf%smin_no3_leached_vr                     , & ! Input:  [real(r8) (:,:)]                                                  
   smin_no3_leached                =>    cnf%smin_no3_leached                        , & ! Input:  [real(r8) (:)]                                                    
   smin_no3_runoff_vr              =>    cnf%smin_no3_runoff_vr                      , & ! Input:  [real(r8) (:,:)]                                                  
   smin_no3_runoff                 =>    cnf%smin_no3_runoff                         , & ! Input:  [real(r8) (:)]                                                    
   decomp_npools                   =>    cns%decomp_npools                           , & ! Input:  [real(r8) (:,:)]  (gN/m2)  decomposing (litter, cwd, soil) N pools
   decomp_npools_vr                =>    cns%decomp_npools_vr                        , & ! Input:  [real(r8) (:,:,:)]  (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
   decomp_npools_1m                =>    cns%decomp_npools_1m                        , & ! Input:  [real(r8) (:,:)]  (gN/m2)  diagnostic: decomposing (litter, cwd, soil) N pools to 1 meter
   altmax_indx                     =>    cps%altmax_indx                             , & ! Input:  [integer (:)]  maximum annual depth of thaw                       
   altmax_lastyear_indx            =>    cps%altmax_lastyear_indx                    , & ! Input:  [integer (:)]  prior year maximum annual depth of thaw            
   sminn_vr                        =>    cns%sminn_vr                                , & ! Input:  [real(r8) (:,:)]  (gN/m3) soil mineral N                          
   col_ntrunc_vr                   =>    cns%col_ntrunc_vr                           , & ! Input:  [real(r8) (:,:)]  (gN/m3) column-level sink for N truncation      
   supplement_to_sminn             =>    cnf%supplement_to_sminn                     , & ! Input:  [real(r8) (:)]                                                    
   supplement_to_sminn_vr          =>    cnf%supplement_to_sminn_vr                  , & ! Input:  [real(r8) (:,:)]                                                  
   col_totpftn                     =>    pns_a%totpftn                               , & ! Input:  [real(r8) (:)]  (gN/m2) total pft-level nitrogen                  
   col_totvegn                     =>    pns_a%totvegn                               , & ! Input:  [real(r8) (:)]  (gN/m2) total vegetation nitrogen                 
   totcoln                         =>    cns%totcoln                                 , & ! Input:  [real(r8) (:)]  (gN/m2) total column nitrogen, incl veg           
   totecosysn                      =>    cns%totecosysn                              , & ! Input:  [real(r8) (:)]  (gN/m2) total ecosystem nitrogen, incl veg        
   totlitn                         =>    cns%totlitn                                 , & ! Input:  [real(r8) (:)]  (gN/m2) total litter nitrogen                     
   totlitn_1m                      =>    cns%totlitn_1m                              , & ! Input:  [real(r8) (:)]  (gN/m2) total litter nitrogen to 1 meter          
   totsomn                         =>    cns%totsomn                                 , & ! Input:  [real(r8) (:)]  (gN/m2) total soil organic matter nitrogen        
   m_leafn_to_fire                 =>    pnf%m_leafn_to_fire                         , & ! Input:  [real(r8) (:)]                                                    
   m_leafn_storage_to_fire         =>    pnf%m_leafn_storage_to_fire                 , & ! Input:  [real(r8) (:)]                                                    
   m_leafn_xfer_to_fire            =>    pnf%m_leafn_xfer_to_fire                    , & ! Input:  [real(r8) (:)]                                                    
   m_livestemn_to_fire             =>    pnf%m_livestemn_to_fire                     , & ! Input:  [real(r8) (:)]                                                    
   m_livestemn_storage_to_fire     =>    pnf%m_livestemn_storage_to_fire             , & ! Input:  [real(r8) (:)]                                                    
   m_livestemn_xfer_to_fire        =>    pnf%m_livestemn_xfer_to_fire                , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemn_to_fire             =>    pnf%m_deadstemn_to_fire                     , & ! Input:  [real(r8) (:)]                                                    
   totsomn_1m                      =>    cns%totsomn_1m                              , & ! Input:  [real(r8) (:)]  (gN/m2) total soil organic matter nitrogen to 1 meter
   m_deadcrootn_storage_to_fire    =>    pnf%m_deadcrootn_storage_to_fire            , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootn_to_fire            =>    pnf%m_deadcrootn_to_fire                    , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootn_xfer_to_fire       =>    pnf%m_deadcrootn_xfer_to_fire               , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemn_storage_to_fire     =>    pnf%m_deadstemn_storage_to_fire             , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemn_xfer_to_fire        =>    pnf%m_deadstemn_xfer_to_fire                , & ! Input:  [real(r8) (:)]                                                    
   m_frootn_to_fire                =>    pnf%m_frootn_to_fire                        , & ! Input:  [real(r8) (:)]                                                    
   m_frootn_storage_to_fire        =>    pnf%m_frootn_storage_to_fire                , & ! Input:  [real(r8) (:)]                                                    
   m_frootn_xfer_to_fire           =>    pnf%m_frootn_xfer_to_fire                   , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootn_to_fire            =>    pnf%m_livecrootn_to_fire                    , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootn_storage_to_fire    =>    pnf%m_livecrootn_storage_to_fire            , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootn_xfer_to_fire       =>    pnf%m_livecrootn_xfer_to_fire               , & ! Input:  [real(r8) (:)]                                                    
   m_retransn_to_fire              =>    pnf%m_retransn_to_fire                      , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadstemn_to_prod10n        =>    pnf%hrv_deadstemn_to_prod10n                , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadstemn_to_prod100n       =>    pnf%hrv_deadstemn_to_prod100n               , & ! Input:  [real(r8) (:)]                                                    
   ndeploy                         =>    pnf%ndeploy                                 , & ! Input:  [real(r8) (:)]                                                    
   pft_fire_nloss                  =>    pnf%pft_fire_nloss                          , & ! Input:  [real(r8) (:)]  (gN/m2/s) total pft-level fire C loss             
   retransn_to_npool               =>    pnf%retransn_to_npool                       , & ! Input:  [real(r8) (:)]                                                    
   sminn_to_npool                  =>    pnf%sminn_to_npool                          , & ! Input:  [real(r8) (:)]                                                    
   deadcrootn                      =>    pns%deadcrootn                              , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N                        
   deadcrootn_storage              =>    pns%deadcrootn_storage                      , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N storage                
   deadcrootn_xfer                 =>    pns%deadcrootn_xfer                         , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N transfer               
   deadstemn                       =>    pns%deadstemn                               , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N                               
   deadstemn_storage               =>    pns%deadstemn_storage                       , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N storage                       
   deadstemn_xfer                  =>    pns%deadstemn_xfer                          , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N transfer                      
   dispvegn                        =>    pns%dispvegn                                , & ! Input:  [real(r8) (:)]  (gN/m2) displayed veg nitrogen, excluding storage 
   frootn                          =>    pns%frootn                                  , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N                               
   frootn_storage                  =>    pns%frootn_storage                          , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N storage                       
   frootn_xfer                     =>    pns%frootn_xfer                             , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N transfer                      
   leafn                           =>    pns%leafn                                   , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N                                    
   leafn_storage                   =>    pns%leafn_storage                           , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N storage                            
   leafn_xfer                      =>    pns%leafn_xfer                              , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N transfer                           
   livecrootn                      =>    pns%livecrootn                              , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N                        
   livecrootn_storage              =>    pns%livecrootn_storage                      , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N storage                
   livecrootn_xfer                 =>    pns%livecrootn_xfer                         , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N transfer               
   grainn                          =>    pns%grainn                                  , & ! Input:  [real(r8) (:)]  (gN/m2) grain N                                   
   grainn_storage                  =>    pns%grainn_storage                          , & ! Input:  [real(r8) (:)]  (gN/m2) grain N storage                           
   grainn_xfer                     =>    pns%grainn_xfer                             , & ! Input:  [real(r8) (:)]  (gN/m2) grain N transfer                          
   livestemn                       =>    pns%livestemn                               , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N                               
   livestemn_storage               =>    pns%livestemn_storage                       , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N storage                       
   livestemn_xfer                  =>    pns%livestemn_xfer                          , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N transfer                      
   retransn                        =>    pns%retransn                                , & ! Input:  [real(r8) (:)]  (gN/m2) plant pool of retranslocated N            
   npool                           =>    pns%npool                                   , & ! Input:  [real(r8) (:)]  (gN/m2) temporary plant N pool                    
   pft_ntrunc                      =>    pns%pft_ntrunc                              , & ! Input:  [real(r8) (:)]  (gN/m2) pft-level sink for N truncation           
   storvegn                        =>    pns%storvegn                                , & ! Input:  [real(r8) (:)]  (gN/m2) stored vegetation nitrogen                
   totpftn                         =>    pns%totpftn                                 , & ! Input:  [real(r8) (:)]  (gN/m2) total pft-level nitrogen                  
   totvegn                         =>    pns%totvegn                                 , & ! Input:  [real(r8) (:)]  (gN/m2) total vegetation nitrogen                 
   wood_harvestn                   =>    pnf%wood_harvestn                           , & ! Input:  [real(r8) (:)]  total N losses to wood product pools (gN/m2/s)    
   col_wood_harvestn               =>    pnf_a%wood_harvestn                         , & ! Input:  [real(r8) (:)]                                                    
   dwt_nloss                       =>    cnf%dwt_nloss                               , & ! Input:  [real(r8) (:)]  (gN/m2/s) total nitrogen loss from product pools and conversion
   dwt_conv_nflux                  =>    cnf%dwt_conv_nflux                          , & ! Input:  [real(r8) (:)]  (gN/m2/s) conversion N flux (immediate loss to atm)
   prod10n_loss                    =>    cnf%prod10n_loss                            , & ! Input:  [real(r8) (:)]  (gN/m2/s) loss from 10-yr wood product pool       
   prod100n_loss                   =>    cnf%prod100n_loss                           , & ! Input:  [real(r8) (:)]  (gN/m2/s) loss from 100-yr wood product pool      
   product_nloss                   =>    cnf%product_nloss                           , & ! Input:  [real(r8) (:)]  (gN/m2/s) total wood product nitrogen loss        
   seedn                           =>    cns%seedn                                   , & ! Input:  [real(r8) (:)]  (gN/m2) column-level pool for seeding new PFTs    
   prod10n                         =>    cns%prod10n                                 , & ! Input:  [real(r8) (:)]  (gN/m2) wood product N pool, 10-year lifespan     
   prod100n                        =>    cns%prod100n                                , & ! Input:  [real(r8) (:)]  (gN/m2) wood product N pool, 100-year lifespan    
   totprodn                        =>    cns%totprodn                                , & ! Input:  [real(r8) (:)]  (gN/m2) total wood product N                      
   som_n_leached                   =>    cnf%som_n_leached                           , & ! Input:  [real(r8) (:)]  total SOM N loss from vertical transport (gN/m^2/s)
   decomp_npools_leached           =>    cnf%decomp_npools_leached                   , & ! Input:  [real(r8) (:,:)]  N loss from vertical transport from each decomposing N pool (gN/m^2/s)
   decomp_npools_transport_tendency=>    cnf%decomp_npools_transport_tendency        , & ! Input:  [real(r8) (:,:,:)]  N tendency due to vertical transport in decomposing N pools (gN/m^3/s)
   decomp_cascade_ntransfer_vr     =>    cnf%decomp_cascade_ntransfer_vr             , & ! Input:  [real(r8) (:,:,:)]                                                
   decomp_cascade_ntransfer        =>    cnf%decomp_cascade_ntransfer                , & ! Input:  [real(r8) (:,:)]                                                  
   decomp_cascade_sminn_flux_vr    =>    cnf%decomp_cascade_sminn_flux_vr            , & ! vert-res mineral N flux for transition along decomposition cascade (gN/m3/s)]
   decomp_cascade_sminn_flux       =>    cnf%decomp_cascade_sminn_flux                 & ! Input:  [real(r8) (:,:)]  vert-int (diagnostic) mineral N flux for transition along decomposition cascade (gN/m2/s)
   )

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
   call p2c(bounds, num_soilc, filter_soilc, &
        pft_fire_nloss(bounds%begp:bounds%endp), &
        col_pft_fire_nloss(bounds%begc:bounds%endc))
   call p2c(bounds, num_soilc, filter_soilc, &
        wood_harvestn(bounds%begp:bounds%endp), &
        col_wood_harvestn(bounds%begc:bounds%endc))
   call p2c(bounds, num_soilc, filter_soilc, &
        totvegn(bounds%begp:bounds%endp), &
        col_totvegn(bounds%begc:bounds%endc))
   call p2c(bounds, num_soilc, filter_soilc, &
        totpftn(bounds%begp:bounds%endp), &
        col_totpftn(bounds%begc:bounds%endc))

   ! column loops
   do fc = 1,num_soilc
      c = filter_soilc(fc)

      ! some zeroing
      denit(c) = 0._r8
      if (use_nitrif_denitrif) then
         smin_no3(c) = 0._r8
         smin_nh4(c) = 0._r8
      end if
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
            decomp_cascade_ntransfer(c,k) = decomp_cascade_ntransfer(c,k) &
               + decomp_cascade_ntransfer_vr(c,j,k) * dzsoi_decomp(j) 
            decomp_cascade_sminn_flux(c,k) = decomp_cascade_sminn_flux(c,k) &
               + decomp_cascade_sminn_flux_vr(c,j,k) * dzsoi_decomp(j) 
         end do
      end do
   end do
   
   if (.not. use_nitrif_denitrif) then
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

   else

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

   end if

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

    end associate 
 end subroutine NSummary

end module CNSummaryMod
