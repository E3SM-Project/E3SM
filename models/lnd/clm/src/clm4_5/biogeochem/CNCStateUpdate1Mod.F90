module CNCStateUpdate1Mod
  !-----------------------------------------------------------------------
  ! Module for carbon state variable update, non-mortality fluxes.
  !
  ! !USES:
  use shr_kind_mod , only: r8 => shr_kind_r8
  use clm_varpar   , only: ndecomp_cascade_transitions, nlevdecomp
  use abortutils   , only: endrun
  use shr_log_mod  , only: errMsg => shr_log_errMsg
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: CStateUpdate1
  public:: CStateUpdate0
  !-----------------------------------------------------------------------

contains


  !-----------------------------------------------------------------------
  subroutine CStateUpdate0(num_soilp, filter_soilp, isotope)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update cpool carbon state
    !
    ! !USES:
    use clmtype
    use clm_time_manager, only: get_step_size
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_soilp       ! number of soil pfts in filter
    integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
    character(len=*), intent(in) :: isotope         ! 'bulk', 'c13' or 'c14'
    !
    ! !LOCAL VARIABLES:
    type(pft_cflux_type), pointer :: pcisof
    type(pft_cstate_type), pointer :: pcisos
    integer :: p     ! indices
    integer :: fp   ! lake filter indices
    real(r8):: dt      ! radiation time step (seconds)
    !-----------------------------------------------------------------------

   ! select which isotope
   select case (isotope)
   case ('bulk')
      pcisof =>  pcf
      pcisos =>  pcs
   case ('c14')
      pcisof =>  pc14f
      pcisos =>  pc14s
   case ('c13')
      pcisof =>  pc13f
      pcisos =>  pc13s
   case default
      call endrun(msg='CNCIsoStateUpdate1Mod: iso must be bulk, c13 or c14'//&
           errMsg(__FILE__, __LINE__))
   end select

   associate(& 
   cpool                               =>    pcisos%cpool                                , & ! InOut:  [real(r8) (:)]  (gC/m2) temporary photosynthate C pool            
   psnshade_to_cpool                   =>    pcisof%psnshade_to_cpool                    , & ! Input:  [real(r8) (:)]                                                    
   psnsun_to_cpool                     =>    pcisof%psnsun_to_cpool                        & ! Input:  [real(r8) (:)]                                                    
   )

    ! set time steps
    dt = real( get_step_size(), r8 )

    ! pft loop
    do fp = 1,num_soilp
       p = filter_soilp(fp)
       ! gross photosynthesis fluxes
       cpool(p) = cpool(p) + psnsun_to_cpool(p)*dt
       cpool(p) = cpool(p) + psnshade_to_cpool(p)*dt
    end do

    end associate 
 end subroutine CStateUpdate0

 !-----------------------------------------------------------------------
 subroutine CStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, isotope)
   !
   ! !DESCRIPTION:
   ! On the radiation time step, update all the prognostic carbon state
   ! variables (except for gap-phase mortality and fire fluxes)
   !
   ! !USES:
   use clmtype
   use clm_time_manager, only: get_step_size
   use clm_varpar   , only: i_met_lit, i_cel_lit, i_lig_lit, i_cwd
   use pftvarcon , only: npcropmin, nc3crop
   !
   ! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
   character(len=*), intent(in) :: isotope         ! 'bulk', 'c13' or 'c14'
   !
   ! !LOCAL VARIABLES:
   type(pft_cflux_type), pointer :: pcisof
   type(pft_cstate_type), pointer :: pcisos
   type(column_cflux_type), pointer :: ccisof
   type(column_cstate_type), pointer :: ccisos
   integer :: c,p,j,k,l     ! indices
   integer :: fp,fc   ! lake filter indices
   real(r8):: dt      ! radiation time step (seconds)
   !-----------------------------------------------------------------------

   ! select which isotope
   select case (isotope)
   case ('bulk')
      pcisof =>  pcf
      pcisos =>  pcs
      ccisof =>  ccf
      ccisos =>  ccs
   case ('c14')
      pcisof =>  pc14f
      pcisos =>  pc14s
      ccisof =>  cc14f
      ccisos =>  cc14s
   case ('c13')
      pcisof =>  pc13f
      pcisos =>  pc13s
      ccisof =>  cc13f
      ccisos =>  cc13s
   case default
      call endrun(msg='CNCIsoStateUpdate1Mod: iso must be bulk, c13 or c14'//&
           errMsg(__FILE__, __LINE__))
   end select

   associate(& 
   woody                               =>    pftcon%woody                                , & ! Input:  [real(r8) (:)]  binary flag for woody lifeform (1=woody, 0=not woody)
   decomp_cpools_vr                    =>    ccisos%decomp_cpools_vr                     , & ! InOut:  [real(r8) (:,:,:)]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
   decomp_cpools_sourcesink            =>    ccisof%decomp_cpools_sourcesink             , & ! InOut:  [real(r8) (:,:,:)]  (gC/m3/timestep)  change in decomposing c pools.  Used to update concentrations concurrently with vertical transport equation
   decomp_cascade_hr_vr                =>    ccisof%decomp_cascade_hr_vr                 , & ! InOut:  [real(r8) (:,:,:)]  vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
   decomp_cascade_ctransfer_vr         =>    ccisof%decomp_cascade_ctransfer_vr          , & ! InOut:  [real(r8) (:,:,:)]  vertically-resolved C transferred along deomposition cascade (gC/m3/s)
   cascade_donor_pool                  =>    decomp_cascade_con%cascade_donor_pool       , & ! InOut:  [integer (:)]  which pool is C taken from for a given decomposition step
   cascade_receiver_pool               =>    decomp_cascade_con%cascade_receiver_pool    , & ! InOut:  [integer (:)]  which pool is C added to for a given decomposition step
   phenology_c_to_litr_met_c           =>    ccisof%phenology_c_to_litr_met_c            , & ! InOut:  [real(r8) (:,:)]  C fluxes associated with phenology (litterfall and crop) to litter metabolic pool (gC/m3/s)
   phenology_c_to_litr_cel_c           =>    ccisof%phenology_c_to_litr_cel_c            , & ! InOut:  [real(r8) (:,:)]  C fluxes associated with phenology (litterfall and crop) to litter cellulose pool (gC/m3/s)
   phenology_c_to_litr_lig_c           =>    ccisof%phenology_c_to_litr_lig_c            , & ! InOut:  [real(r8) (:,:)]  C fluxes associated with phenology (litterfall and crop) to litter lignin pool (gC/m3/s)
   dwt_seedc_to_leaf                   =>    ccisof%dwt_seedc_to_leaf                    , & ! InOut:  [real(r8) (:)]                                                    
   dwt_seedc_to_deadstem               =>    ccisof%dwt_seedc_to_deadstem                , & ! InOut:  [real(r8) (:)]                                                    
   dwt_frootc_to_litr_met_c            =>    ccisof%dwt_frootc_to_litr_met_c             , & ! InOut:  [real(r8) (:,:)]                                                  
   dwt_frootc_to_litr_cel_c            =>    ccisof%dwt_frootc_to_litr_cel_c             , & ! InOut:  [real(r8) (:,:)]                                                  
   dwt_frootc_to_litr_lig_c            =>    ccisof%dwt_frootc_to_litr_lig_c             , & ! InOut:  [real(r8) (:,:)]                                                  
   dwt_livecrootc_to_cwdc              =>    ccisof%dwt_livecrootc_to_cwdc               , & ! InOut:  [real(r8) (:,:)]  (gC/m2/s) live coarse root to CWD due to landcover change
   dwt_deadcrootc_to_cwdc              =>    ccisof%dwt_deadcrootc_to_cwdc               , & ! InOut:  [real(r8) (:,:)]  (gC/m2/s) dead coarse root to CWD due to landcover change
   seedc                               =>    ccisos%seedc                                , & ! InOut:  [real(r8) (:)]                                                    
   ivt                                 =>   pft%itype                                    , & ! Input:  [integer (:)]  pft vegetation type                                
   cpool_deadcroot_gr                  =>    pcisof%cpool_deadcroot_gr                   , & ! Input:  [real(r8) (:)]                                                    
   cpool_deadcroot_storage_gr          =>    pcisof%cpool_deadcroot_storage_gr           , & ! Input:  [real(r8) (:)]                                                    
   cpool_deadstem_gr                   =>    pcisof%cpool_deadstem_gr                    , & ! Input:  [real(r8) (:)]                                                    
   cpool_deadstem_storage_gr           =>    pcisof%cpool_deadstem_storage_gr            , & ! Input:  [real(r8) (:)]                                                    
   cpool_froot_gr                      =>    pcisof%cpool_froot_gr                       , & ! Input:  [real(r8) (:)]                                                    
   cpool_froot_storage_gr              =>    pcisof%cpool_froot_storage_gr               , & ! Input:  [real(r8) (:)]                                                    
   cpool_leaf_gr                       =>    pcisof%cpool_leaf_gr                        , & ! Input:  [real(r8) (:)]                                                    
   cpool_leaf_storage_gr               =>    pcisof%cpool_leaf_storage_gr                , & ! Input:  [real(r8) (:)]                                                    
   cpool_livecroot_gr                  =>    pcisof%cpool_livecroot_gr                   , & ! Input:  [real(r8) (:)]                                                    
   cpool_livecroot_storage_gr          =>    pcisof%cpool_livecroot_storage_gr           , & ! Input:  [real(r8) (:)]                                                    
   cpool_livestem_gr                   =>    pcisof%cpool_livestem_gr                    , & ! Input:  [real(r8) (:)]  live stem growth respiration (gC/m2/s)            
   cpool_livestem_storage_gr           =>    pcisof%cpool_livestem_storage_gr            , & ! Input:  [real(r8) (:)]  live stem growth respiration to storage (gC/m2/s) 
   cpool_to_xsmrpool                   =>    pcisof%cpool_to_xsmrpool                    , & ! Input:  [real(r8) (:)]                                                    
   cpool_to_deadcrootc                 =>    pcisof%cpool_to_deadcrootc                  , & ! Input:  [real(r8) (:)]                                                    
   cpool_to_deadcrootc_storage         =>    pcisof%cpool_to_deadcrootc_storage          , & ! Input:  [real(r8) (:)]                                                    
   cpool_to_deadstemc                  =>    pcisof%cpool_to_deadstemc                   , & ! Input:  [real(r8) (:)]                                                    
   cpool_to_deadstemc_storage          =>    pcisof%cpool_to_deadstemc_storage           , & ! Input:  [real(r8) (:)]                                                    
   cpool_to_frootc                     =>    pcisof%cpool_to_frootc                      , & ! Input:  [real(r8) (:)]                                                    
   cpool_to_frootc_storage             =>    pcisof%cpool_to_frootc_storage              , & ! Input:  [real(r8) (:)]                                                    
   cpool_to_gresp_storage              =>    pcisof%cpool_to_gresp_storage               , & ! Input:  [real(r8) (:)]                                                    
   cpool_to_leafc                      =>    pcisof%cpool_to_leafc                       , & ! Input:  [real(r8) (:)]                                                    
   cpool_to_leafc_storage              =>    pcisof%cpool_to_leafc_storage               , & ! Input:  [real(r8) (:)]                                                    
   cpool_to_livecrootc                 =>    pcisof%cpool_to_livecrootc                  , & ! Input:  [real(r8) (:)]                                                    
   cpool_to_livecrootc_storage         =>    pcisof%cpool_to_livecrootc_storage          , & ! Input:  [real(r8) (:)]                                                    
   cpool_to_livestemc                  =>    pcisof%cpool_to_livestemc                   , & ! Input:  [real(r8) (:)]                                                    
   cpool_to_livestemc_storage          =>    pcisof%cpool_to_livestemc_storage           , & ! Input:  [real(r8) (:)]                                                    
   deadcrootc_storage_to_xfer          =>    pcisof%deadcrootc_storage_to_xfer           , & ! Input:  [real(r8) (:)]                                                    
   deadcrootc_xfer_to_deadcrootc       =>    pcisof%deadcrootc_xfer_to_deadcrootc        , & ! Input:  [real(r8) (:)]                                                    
   deadstemc_storage_to_xfer           =>    pcisof%deadstemc_storage_to_xfer            , & ! Input:  [real(r8) (:)]                                                    
   deadstemc_xfer_to_deadstemc         =>    pcisof%deadstemc_xfer_to_deadstemc          , & ! Input:  [real(r8) (:)]                                                    
   froot_curmr                         =>    pcisof%froot_curmr                          , & ! Input:  [real(r8) (:)]                                                    
   froot_xsmr                          =>    pcisof%froot_xsmr                           , & ! Input:  [real(r8) (:)]                                                    
   frootc_storage_to_xfer              =>    pcisof%frootc_storage_to_xfer               , & ! Input:  [real(r8) (:)]                                                    
   frootc_to_litter                    =>    pcisof%frootc_to_litter                     , & ! Input:  [real(r8) (:)]                                                    
   frootc_xfer_to_frootc               =>    pcisof%frootc_xfer_to_frootc                , & ! Input:  [real(r8) (:)]                                                    
   gresp_storage_to_xfer               =>    pcisof%gresp_storage_to_xfer                , & ! Input:  [real(r8) (:)]                                                    
   leaf_curmr                          =>    pcisof%leaf_curmr                           , & ! Input:  [real(r8) (:)]                                                    
   leaf_xsmr                           =>    pcisof%leaf_xsmr                            , & ! Input:  [real(r8) (:)]                                                    
   leafc_storage_to_xfer               =>    pcisof%leafc_storage_to_xfer                , & ! Input:  [real(r8) (:)]                                                    
   leafc_to_litter                     =>    pcisof%leafc_to_litter                      , & ! Input:  [real(r8) (:)]                                                    
   leafc_xfer_to_leafc                 =>    pcisof%leafc_xfer_to_leafc                  , & ! Input:  [real(r8) (:)]                                                    
   livecroot_curmr                     =>    pcisof%livecroot_curmr                      , & ! Input:  [real(r8) (:)]                                                    
   livecroot_xsmr                      =>    pcisof%livecroot_xsmr                       , & ! Input:  [real(r8) (:)]                                                    
   livecrootc_storage_to_xfer          =>    pcisof%livecrootc_storage_to_xfer           , & ! Input:  [real(r8) (:)]                                                    
   livecrootc_to_deadcrootc            =>    pcisof%livecrootc_to_deadcrootc             , & ! Input:  [real(r8) (:)]                                                    
   livecrootc_xfer_to_livecrootc       =>    pcisof%livecrootc_xfer_to_livecrootc        , & ! Input:  [real(r8) (:)]                                                    
   livestem_curmr                      =>    pcisof%livestem_curmr                       , & ! Input:  [real(r8) (:)]                                                    
   livestem_xsmr                       =>    pcisof%livestem_xsmr                        , & ! Input:  [real(r8) (:)]                                                    
   livestemc_storage_to_xfer           =>    pcisof%livestemc_storage_to_xfer            , & ! Input:  [real(r8) (:)]                                                    
   livestemc_to_deadstemc              =>    pcisof%livestemc_to_deadstemc               , & ! Input:  [real(r8) (:)]                                                    
   livestemc_xfer_to_livestemc         =>    pcisof%livestemc_xfer_to_livestemc          , & ! Input:  [real(r8) (:)]                                                    
   transfer_deadcroot_gr               =>    pcisof%transfer_deadcroot_gr                , & ! Input:  [real(r8) (:)]  dead coarse root growth respiration from storage (gC/m2/s)
   transfer_deadstem_gr                =>    pcisof%transfer_deadstem_gr                 , & ! Input:  [real(r8) (:)]  dead stem growth respiration from storage (gC/m2/s)
   transfer_froot_gr                   =>    pcisof%transfer_froot_gr                    , & ! Input:  [real(r8) (:)]  fine root  growth respiration from storage (gC/m2/s)
   transfer_leaf_gr                    =>    pcisof%transfer_leaf_gr                     , & ! Input:  [real(r8) (:)]  leaf growth respiration from storage (gC/m2/s)    
   transfer_livecroot_gr               =>    pcisof%transfer_livecroot_gr                , & ! Input:  [real(r8) (:)]  live coarse root growth respiration from storage (gC/m2/s)
   transfer_livestem_gr                =>    pcisof%transfer_livestem_gr                 , & ! Input:  [real(r8) (:)]  live stem growth respiration from storage (gC/m2/s)
   harvdate                            =>    pps%harvdate                                , & ! Input:  [integer (:)]  harvest date                                       
   xsmrpool_to_atm                     =>    pcisof%xsmrpool_to_atm                      , & ! Input:  [real(r8) (:)]  excess MR pool harvest mortality (gC/m2/s)        
   cpool_grain_gr                      =>    pcisof%cpool_grain_gr                       , & ! Input:  [real(r8) (:)]  grain growth respiration (gC/m2/s)                
   cpool_grain_storage_gr              =>    pcisof%cpool_grain_storage_gr               , & ! Input:  [real(r8) (:)]  grain growth respiration to storage (gC/m2/s)     
   cpool_to_grainc                     =>    pcisof%cpool_to_grainc                      , & ! Input:  [real(r8) (:)]  allocation to grain C (gC/m2/s)                   
   cpool_to_grainc_storage             =>    pcisof%cpool_to_grainc_storage              , & ! Input:  [real(r8) (:)]  allocation to grain C storage (gC/m2/s)           
   livestemc_to_litter                 =>    pcisof%livestemc_to_litter                  , & ! Input:  [real(r8) (:)]  live stem C litterfall (gC/m2/s)                  
   grain_curmr                         =>    pcf%grain_curmr                             , & ! Input:  [real(r8) (:)]                                                    
   grain_xsmr                          =>    pcf%grain_xsmr                              , & ! Input:  [real(r8) (:)]                                                    
   grainc_storage_to_xfer              =>    pcisof%grainc_storage_to_xfer               , & ! Input:  [real(r8) (:)]  grain C shift storage to transfer (gC/m2/s)       
   grainc_to_food                      =>    pcisof%grainc_to_food                       , & ! Input:  [real(r8) (:)]  grain C to food (gC/m2/s)                         
   grainc_xfer_to_grainc               =>    pcisof%grainc_xfer_to_grainc                , & ! Input:  [real(r8) (:)]  grain C growth from storage (gC/m2/s)             
   transfer_grain_gr                   =>    pcisof%transfer_grain_gr                    , & ! Input:  [real(r8) (:)]  grain growth respiration from storage (gC/m2/s)   
   grainc                              =>    pcisos%grainc                               , & ! InOut:  [real(r8) (:)]  (gC/m2) grain C                                   
   grainc_storage                      =>    pcisos%grainc_storage                       , & ! InOut:  [real(r8) (:)]  (gC/m2) grain C storage                           
   grainc_xfer                         =>    pcisos%grainc_xfer                          , & ! InOut:  [real(r8) (:)]  (gC/m2) grain C transfer                          
   cpool                               =>    pcisos%cpool                                , & ! InOut:  [real(r8) (:)]  (gC/m2) temporary photosynthate C pool            
   xsmrpool                            =>    pcisos%xsmrpool                             , & ! InOut:  [real(r8) (:)]  (gC/m2) execss maint resp C pool                  
   deadcrootc                          =>    pcisos%deadcrootc                           , & ! InOut:  [real(r8) (:)]  (gC/m2) dead coarse root C                        
   deadcrootc_storage                  =>    pcisos%deadcrootc_storage                   , & ! InOut:  [real(r8) (:)]  (gC/m2) dead coarse root C storage                
   deadcrootc_xfer                     =>    pcisos%deadcrootc_xfer                      , & ! InOut:  [real(r8) (:)]  (gC/m2) dead coarse root C transfer               
   deadstemc                           =>    pcisos%deadstemc                            , & ! InOut:  [real(r8) (:)]  (gC/m2) dead stem C                               
   deadstemc_storage                   =>    pcisos%deadstemc_storage                    , & ! InOut:  [real(r8) (:)]  (gC/m2) dead stem C storage                       
   deadstemc_xfer                      =>    pcisos%deadstemc_xfer                       , & ! InOut:  [real(r8) (:)]  (gC/m2) dead stem C transfer                      
   frootc                              =>    pcisos%frootc                               , & ! InOut:  [real(r8) (:)]  (gC/m2) fine root C                               
   frootc_storage                      =>    pcisos%frootc_storage                       , & ! InOut:  [real(r8) (:)]  (gC/m2) fine root C storage                       
   frootc_xfer                         =>    pcisos%frootc_xfer                          , & ! InOut:  [real(r8) (:)]  (gC/m2) fine root C transfer                      
   gresp_storage                       =>    pcisos%gresp_storage                        , & ! InOut:  [real(r8) (:)]  (gC/m2) growth respiration storage                
   gresp_xfer                          =>    pcisos%gresp_xfer                           , & ! InOut:  [real(r8) (:)]  (gC/m2) growth respiration transfer               
   leafc                               =>    pcisos%leafc                                , & ! InOut:  [real(r8) (:)]  (gC/m2) leaf C                                    
   leafc_storage                       =>    pcisos%leafc_storage                        , & ! InOut:  [real(r8) (:)]  (gC/m2) leaf C storage                            
   leafc_xfer                          =>    pcisos%leafc_xfer                           , & ! InOut:  [real(r8) (:)]  (gC/m2) leaf C transfer                           
   livecrootc                          =>    pcisos%livecrootc                           , & ! InOut:  [real(r8) (:)]  (gC/m2) live coarse root C                        
   livecrootc_storage                  =>    pcisos%livecrootc_storage                   , & ! InOut:  [real(r8) (:)]  (gC/m2) live coarse root C storage                
   livecrootc_xfer                     =>    pcisos%livecrootc_xfer                      , & ! InOut:  [real(r8) (:)]  (gC/m2) live coarse root C transfer               
   livestemc                           =>    pcisos%livestemc                            , & ! InOut:  [real(r8) (:)]  (gC/m2) live stem C                               
   livestemc_storage                   =>    pcisos%livestemc_storage                    , & ! InOut:  [real(r8) (:)]  (gC/m2) live stem C storage                       
   livestemc_xfer                      =>    pcisos%livestemc_xfer                         & ! InOut:  [real(r8) (:)]  (gC/m2) live stem C transfer                      
   )

    ! set time steps
    dt = real( get_step_size(), r8 )

    ! column level fluxes
    
    ! column loop
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       ! seeding fluxes, from dynamic landcover
       seedc(c) = seedc(c) - dwt_seedc_to_leaf(c) * dt
       seedc(c) = seedc(c) - dwt_seedc_to_deadstem(c) * dt
    end do
    
    ! plant to litter fluxes
    do j = 1,nlevdecomp
       ! column loop
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          ! phenology and dynamic land cover fluxes
          decomp_cpools_sourcesink(c,j,i_met_lit) = ( phenology_c_to_litr_met_c(c,j) + dwt_frootc_to_litr_met_c(c,j) ) *dt
          decomp_cpools_sourcesink(c,j,i_cel_lit) = ( phenology_c_to_litr_cel_c(c,j) + dwt_frootc_to_litr_cel_c(c,j) ) *dt
          decomp_cpools_sourcesink(c,j,i_lig_lit) = ( phenology_c_to_litr_lig_c(c,j) + dwt_frootc_to_litr_lig_c(c,j) ) *dt
          decomp_cpools_sourcesink(c,j,i_cwd) = ( dwt_livecrootc_to_cwdc(c,j) + dwt_deadcrootc_to_cwdc(c,j) ) *dt
       end do
    end do
    
    ! litter and SOM HR fluxes
    do k = 1, ndecomp_cascade_transitions
       do j = 1,nlevdecomp
          ! column loop
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             decomp_cpools_sourcesink(c,j,cascade_donor_pool(k)) = &
                  decomp_cpools_sourcesink(c,j,cascade_donor_pool(k)) &
                    - ( decomp_cascade_hr_vr(c,j,k) + decomp_cascade_ctransfer_vr(c,j,k)) *dt
          end do
       end do
    end do
    do k = 1, ndecomp_cascade_transitions
       if ( cascade_receiver_pool(k) .ne. 0 ) then  ! skip terminal transitions
          do j = 1,nlevdecomp
             ! column loop
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                decomp_cpools_sourcesink(c,j,cascade_receiver_pool(k)) = &
                     decomp_cpools_sourcesink(c,j,cascade_receiver_pool(k)) &
                        + decomp_cascade_ctransfer_vr(c,j,k)*dt
             end do
          end do
       end if
    end do

    ! pft loop
    do fp = 1,num_soilp
       p = filter_soilp(fp)
 
       ! phenology: transfer growth fluxes
       leafc(p)           = leafc(p)       + leafc_xfer_to_leafc(p)*dt
       leafc_xfer(p)      = leafc_xfer(p)  - leafc_xfer_to_leafc(p)*dt
       frootc(p)          = frootc(p)      + frootc_xfer_to_frootc(p)*dt
       frootc_xfer(p)     = frootc_xfer(p) - frootc_xfer_to_frootc(p)*dt
       if (woody(ivt(p)) == 1._r8) then
           livestemc(p)       = livestemc(p)           + livestemc_xfer_to_livestemc(p)*dt
           livestemc_xfer(p)  = livestemc_xfer(p)  - livestemc_xfer_to_livestemc(p)*dt
           deadstemc(p)       = deadstemc(p)           + deadstemc_xfer_to_deadstemc(p)*dt
           deadstemc_xfer(p)  = deadstemc_xfer(p)  - deadstemc_xfer_to_deadstemc(p)*dt
           livecrootc(p)      = livecrootc(p)          + livecrootc_xfer_to_livecrootc(p)*dt
           livecrootc_xfer(p) = livecrootc_xfer(p) - livecrootc_xfer_to_livecrootc(p)*dt
           deadcrootc(p)      = deadcrootc(p)          + deadcrootc_xfer_to_deadcrootc(p)*dt
           deadcrootc_xfer(p) = deadcrootc_xfer(p) - deadcrootc_xfer_to_deadcrootc(p)*dt
       end if
       if (ivt(p) >= npcropmin) then ! skip 2 generic crops
           ! lines here for consistency; the transfer terms are zero
           livestemc(p)       = livestemc(p)      + livestemc_xfer_to_livestemc(p)*dt
           livestemc_xfer(p)  = livestemc_xfer(p) - livestemc_xfer_to_livestemc(p)*dt
           grainc(p)          = grainc(p)         + grainc_xfer_to_grainc(p)*dt
           grainc_xfer(p)     = grainc_xfer(p)    - grainc_xfer_to_grainc(p)*dt
       end if
 
       ! phenology: litterfall fluxes
       leafc(p) = leafc(p) - leafc_to_litter(p)*dt
       frootc(p) = frootc(p) - frootc_to_litter(p)*dt
 
       ! livewood turnover fluxes
       if (woody(ivt(p)) == 1._r8) then
           livestemc(p)  = livestemc(p)  - livestemc_to_deadstemc(p)*dt
           deadstemc(p)  = deadstemc(p)  + livestemc_to_deadstemc(p)*dt
           livecrootc(p) = livecrootc(p) - livecrootc_to_deadcrootc(p)*dt
           deadcrootc(p) = deadcrootc(p) + livecrootc_to_deadcrootc(p)*dt
       end if
       if (ivt(p) >= npcropmin) then ! skip 2 generic crops
           livestemc(p)  = livestemc(p)  - livestemc_to_litter(p)*dt
           grainc(p)     = grainc(p)     - grainc_to_food(p)*dt
       end if
 
       ! maintenance respiration fluxes from cpool
       cpool(p) = cpool(p) - cpool_to_xsmrpool(p)*dt
       cpool(p) = cpool(p) - leaf_curmr(p)*dt
       cpool(p) = cpool(p) - froot_curmr(p)*dt
       if (woody(ivt(p)) == 1._r8) then
           cpool(p) = cpool(p) - livestem_curmr(p)*dt
           cpool(p) = cpool(p) - livecroot_curmr(p)*dt
       end if
       if (ivt(p) >= npcropmin) then ! skip 2 generic crops
           cpool(p) = cpool(p) - livestem_curmr(p)*dt
           cpool(p) = cpool(p) - grain_curmr(p)*dt
       end if

       ! maintenance respiration fluxes from xsmrpool
       xsmrpool(p) = xsmrpool(p) + cpool_to_xsmrpool(p)*dt
       xsmrpool(p) = xsmrpool(p) - leaf_xsmr(p)*dt
       xsmrpool(p) = xsmrpool(p) - froot_xsmr(p)*dt
       if (woody(ivt(p)) == 1._r8) then
           xsmrpool(p) = xsmrpool(p) - livestem_xsmr(p)*dt
           xsmrpool(p) = xsmrpool(p) - livecroot_xsmr(p)*dt
       end if
       if (ivt(p) >= npcropmin) then ! skip 2 generic crops
           xsmrpool(p) = xsmrpool(p) - livestem_xsmr(p)*dt
           xsmrpool(p) = xsmrpool(p) - grain_xsmr(p)*dt
           if (harvdate(p) < 999) then ! beginning at harvest, send to atm
              xsmrpool_to_atm(p) = xsmrpool_to_atm(p) + xsmrpool(p)/dt
              xsmrpool(p) = xsmrpool(p) - xsmrpool_to_atm(p)*dt
           end if
       end if
 
       ! allocation fluxes
       cpool(p)           = cpool(p)          - cpool_to_leafc(p)*dt
       leafc(p)           = leafc(p)          + cpool_to_leafc(p)*dt
       cpool(p)           = cpool(p)          - cpool_to_leafc_storage(p)*dt
       leafc_storage(p)   = leafc_storage(p)  + cpool_to_leafc_storage(p)*dt
       cpool(p)           = cpool(p)          - cpool_to_frootc(p)*dt
       frootc(p)          = frootc(p)         + cpool_to_frootc(p)*dt
       cpool(p)           = cpool(p)          - cpool_to_frootc_storage(p)*dt
       frootc_storage(p)  = frootc_storage(p) + cpool_to_frootc_storage(p)*dt
       if (woody(ivt(p)) == 1._r8) then
           cpool(p)               = cpool(p)              - cpool_to_livestemc(p)*dt
           livestemc(p)           = livestemc(p)          + cpool_to_livestemc(p)*dt
           cpool(p)               = cpool(p)              - cpool_to_livestemc_storage(p)*dt
           livestemc_storage(p)   = livestemc_storage(p)  + cpool_to_livestemc_storage(p)*dt
           cpool(p)               = cpool(p)              - cpool_to_deadstemc(p)*dt
           deadstemc(p)           = deadstemc(p)          + cpool_to_deadstemc(p)*dt
           cpool(p)               = cpool(p)              - cpool_to_deadstemc_storage(p)*dt
           deadstemc_storage(p)   = deadstemc_storage(p)  + cpool_to_deadstemc_storage(p)*dt
           cpool(p)               = cpool(p)              - cpool_to_livecrootc(p)*dt
           livecrootc(p)          = livecrootc(p)         + cpool_to_livecrootc(p)*dt
           cpool(p)               = cpool(p)              - cpool_to_livecrootc_storage(p)*dt
           livecrootc_storage(p)  = livecrootc_storage(p) + cpool_to_livecrootc_storage(p)*dt
           cpool(p)               = cpool(p)              - cpool_to_deadcrootc(p)*dt
           deadcrootc(p)          = deadcrootc(p)         + cpool_to_deadcrootc(p)*dt
           cpool(p)               = cpool(p)              - cpool_to_deadcrootc_storage(p)*dt
           deadcrootc_storage(p)  = deadcrootc_storage(p) + cpool_to_deadcrootc_storage(p)*dt
       end if
       if (ivt(p) >= npcropmin) then ! skip 2 generic crops
           cpool(p)               = cpool(p)              - cpool_to_livestemc(p)*dt
           livestemc(p)           = livestemc(p)          + cpool_to_livestemc(p)*dt
           cpool(p)               = cpool(p)              - cpool_to_livestemc_storage(p)*dt
           livestemc_storage(p)   = livestemc_storage(p)  + cpool_to_livestemc_storage(p)*dt
           cpool(p)               = cpool(p)              - cpool_to_grainc(p)*dt
           grainc(p)              = grainc(p)             + cpool_to_grainc(p)*dt
           cpool(p)               = cpool(p)              - cpool_to_grainc_storage(p)*dt
           grainc_storage(p)      = grainc_storage(p)     + cpool_to_grainc_storage(p)*dt
       end if
 
       ! growth respiration fluxes for current growth
       cpool(p) = cpool(p) - cpool_leaf_gr(p)*dt
       cpool(p) = cpool(p) - cpool_froot_gr(p)*dt
       if (woody(ivt(p)) == 1._r8) then
           cpool(p) = cpool(p) - cpool_livestem_gr(p)*dt
           cpool(p) = cpool(p) - cpool_deadstem_gr(p)*dt
           cpool(p) = cpool(p) - cpool_livecroot_gr(p)*dt
           cpool(p) = cpool(p) - cpool_deadcroot_gr(p)*dt
       end if
       if (ivt(p) >= npcropmin) then ! skip 2 generic crops
           cpool(p) = cpool(p) - cpool_livestem_gr(p)*dt
           cpool(p) = cpool(p) - cpool_grain_gr(p)*dt
       end if
 
       ! growth respiration for transfer growth
       gresp_xfer(p) = gresp_xfer(p) - transfer_leaf_gr(p)*dt
       gresp_xfer(p) = gresp_xfer(p) - transfer_froot_gr(p)*dt
       if (woody(ivt(p)) == 1._r8) then
           gresp_xfer(p) = gresp_xfer(p) - transfer_livestem_gr(p)*dt
           gresp_xfer(p) = gresp_xfer(p) - transfer_deadstem_gr(p)*dt
           gresp_xfer(p) = gresp_xfer(p) - transfer_livecroot_gr(p)*dt
           gresp_xfer(p) = gresp_xfer(p) - transfer_deadcroot_gr(p)*dt
       end if
       if (ivt(p) >= npcropmin) then ! skip 2 generic crops
           gresp_xfer(p) = gresp_xfer(p) - transfer_livestem_gr(p)*dt
           gresp_xfer(p) = gresp_xfer(p) - transfer_grain_gr(p)*dt
       end if
 
       ! growth respiration at time of storage
       cpool(p) = cpool(p) - cpool_leaf_storage_gr(p)*dt
       cpool(p) = cpool(p) - cpool_froot_storage_gr(p)*dt
       if (woody(ivt(p)) == 1._r8) then
           cpool(p) = cpool(p) - cpool_livestem_storage_gr(p)*dt
           cpool(p) = cpool(p) - cpool_deadstem_storage_gr(p)*dt
           cpool(p) = cpool(p) - cpool_livecroot_storage_gr(p)*dt
           cpool(p) = cpool(p) - cpool_deadcroot_storage_gr(p)*dt
       end if
       if (ivt(p) >= npcropmin) then ! skip 2 generic crops
           cpool(p) = cpool(p) - cpool_livestem_storage_gr(p)*dt
           cpool(p) = cpool(p) - cpool_grain_storage_gr(p)*dt
       end if
 
       ! growth respiration stored for release during transfer growth
       cpool(p)         = cpool(p)         - cpool_to_gresp_storage(p)*dt
       gresp_storage(p) = gresp_storage(p) + cpool_to_gresp_storage(p)*dt
 
       ! move storage pools into transfer pools
       leafc_storage(p)   = leafc_storage(p)   - leafc_storage_to_xfer(p)*dt
       leafc_xfer(p)  = leafc_xfer(p)  + leafc_storage_to_xfer(p)*dt
       frootc_storage(p)  = frootc_storage(p)  - frootc_storage_to_xfer(p)*dt
       frootc_xfer(p) = frootc_xfer(p) + frootc_storage_to_xfer(p)*dt
       if (woody(ivt(p)) == 1._r8) then
           livestemc_storage(p)  = livestemc_storage(p)   - livestemc_storage_to_xfer(p)*dt
           livestemc_xfer(p)     = livestemc_xfer(p)  + livestemc_storage_to_xfer(p)*dt
           deadstemc_storage(p)  = deadstemc_storage(p)   - deadstemc_storage_to_xfer(p)*dt
           deadstemc_xfer(p)     = deadstemc_xfer(p)  + deadstemc_storage_to_xfer(p)*dt
           livecrootc_storage(p) = livecrootc_storage(p)  - livecrootc_storage_to_xfer(p)*dt
           livecrootc_xfer(p)    = livecrootc_xfer(p) + livecrootc_storage_to_xfer(p)*dt
           deadcrootc_storage(p) = deadcrootc_storage(p)  - deadcrootc_storage_to_xfer(p)*dt
           deadcrootc_xfer(p)    = deadcrootc_xfer(p) + deadcrootc_storage_to_xfer(p)*dt
           gresp_storage(p)      = gresp_storage(p)       - gresp_storage_to_xfer(p)*dt
           gresp_xfer(p)         = gresp_xfer(p)      + gresp_storage_to_xfer(p)*dt
       end if
       if (ivt(p) >= npcropmin) then ! skip 2 generic crops
           ! lines here for consistency; the transfer terms are zero
           livestemc_storage(p)  = livestemc_storage(p) - livestemc_storage_to_xfer(p)*dt
           livestemc_xfer(p)     = livestemc_xfer(p)    + livestemc_storage_to_xfer(p)*dt
           grainc_storage(p)     = grainc_storage(p)    - grainc_storage_to_xfer(p)*dt
           grainc_xfer(p)        = grainc_xfer(p)       + grainc_storage_to_xfer(p)*dt
       end if
 
    end do ! end of pft loop

    end associate 
 end subroutine CStateUpdate1

end module CNCStateUpdate1Mod
