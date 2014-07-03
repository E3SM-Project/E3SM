module CNNStateUpdate2Mod
!-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for nitrogen state variable update, mortality fluxes.
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varpar   , only: nlevsoi, nlevdecomp
  implicit none
  save
  private
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: NStateUpdate2
  public:: NStateUpdate2h
  !
  ! !REVISION HISTORY:
  ! 4/23/2004: Created by Peter Thornton
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
subroutine NStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp)
  !
  ! !DESCRIPTION:
  ! On the radiation time step, update all the prognostic nitrogen state
  ! variables affected by gap-phase mortality fluxes
  !
  ! !USES:
  use clmtype
  use clm_time_manager, only: get_step_size
  use clm_varctl  , only: iulog
  use clm_varpar   , only: i_met_lit, i_cel_lit, i_lig_lit, i_cwd
  !
  ! !ARGUMENTS:
  implicit none
  integer, intent(in) :: num_soilc       ! number of soil columns in filter
  integer, intent(in) :: filter_soilc(:) ! filter for soil columns
  integer, intent(in) :: num_soilp       ! number of soil pfts in filter
  integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
  !
  ! !LOCAL VARIABLES:
  integer :: c,p,j,l         ! indices
  integer :: fp,fc       ! lake filter indices
  real(r8):: dt          ! radiation time step (seconds)
!-----------------------------------------------------------------------

   associate(& 
   gap_mortality_n_to_litr_met_n       =>    cnf%gap_mortality_n_to_litr_met_n           , & ! Input:  [real(r8) (:,:)]  N fluxes associated with gap mortality to litter metabolic pool (gN/m3/s)
   gap_mortality_n_to_litr_cel_n       =>    cnf%gap_mortality_n_to_litr_cel_n           , & ! Input:  [real(r8) (:,:)]  N fluxes associated with gap mortality to litter cellulose pool (gN/m3/s)
   gap_mortality_n_to_litr_lig_n       =>    cnf%gap_mortality_n_to_litr_lig_n           , & ! Input:  [real(r8) (:,:)]  N fluxes associated with gap mortality to litter lignin pool (gN/m3/s)
   gap_mortality_n_to_cwdn             =>    cnf%gap_mortality_n_to_cwdn                 , & ! Input:  [real(r8) (:,:)]  N fluxes associated with gap mortality to CWD pool (gN/m3/s)
   decomp_npools_vr                    =>    cns%decomp_npools_vr                        , & ! InOut:  [real(r8) (:,:,:)]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
   m_deadcrootn_storage_to_litter      =>    pnf%m_deadcrootn_storage_to_litter          , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootn_to_litter              =>    pnf%m_deadcrootn_to_litter                  , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootn_xfer_to_litter         =>    pnf%m_deadcrootn_xfer_to_litter             , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemn_storage_to_litter       =>    pnf%m_deadstemn_storage_to_litter           , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemn_to_litter               =>    pnf%m_deadstemn_to_litter                   , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemn_xfer_to_litter          =>    pnf%m_deadstemn_xfer_to_litter              , & ! Input:  [real(r8) (:)]                                                    
   m_frootn_storage_to_litter          =>    pnf%m_frootn_storage_to_litter              , & ! Input:  [real(r8) (:)]                                                    
   m_frootn_to_litter                  =>    pnf%m_frootn_to_litter                      , & ! Input:  [real(r8) (:)]                                                    
   m_frootn_xfer_to_litter             =>    pnf%m_frootn_xfer_to_litter                 , & ! Input:  [real(r8) (:)]                                                    
   m_leafn_storage_to_litter           =>    pnf%m_leafn_storage_to_litter               , & ! Input:  [real(r8) (:)]                                                    
   m_leafn_to_litter                   =>    pnf%m_leafn_to_litter                       , & ! Input:  [real(r8) (:)]                                                    
   m_leafn_xfer_to_litter              =>    pnf%m_leafn_xfer_to_litter                  , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootn_storage_to_litter      =>    pnf%m_livecrootn_storage_to_litter          , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootn_to_litter              =>    pnf%m_livecrootn_to_litter                  , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootn_xfer_to_litter         =>    pnf%m_livecrootn_xfer_to_litter             , & ! Input:  [real(r8) (:)]                                                    
   m_livestemn_storage_to_litter       =>    pnf%m_livestemn_storage_to_litter           , & ! Input:  [real(r8) (:)]                                                    
   m_livestemn_to_litter               =>    pnf%m_livestemn_to_litter                   , & ! Input:  [real(r8) (:)]                                                    
   m_livestemn_xfer_to_litter          =>    pnf%m_livestemn_xfer_to_litter              , & ! Input:  [real(r8) (:)]                                                    
   m_retransn_to_litter                =>    pnf%m_retransn_to_litter                    , & ! Input:  [real(r8) (:)]                                                    
   deadcrootn                          =>    pns%deadcrootn                              , & ! InOut:  [real(r8) (:)]  (gN/m2) dead coarse root N                        
   deadcrootn_storage                  =>    pns%deadcrootn_storage                      , & ! InOut:  [real(r8) (:)]  (gN/m2) dead coarse root N storage                
   deadcrootn_xfer                     =>    pns%deadcrootn_xfer                         , & ! InOut:  [real(r8) (:)]  (gN/m2) dead coarse root N transfer               
   deadstemn                           =>    pns%deadstemn                               , & ! InOut:  [real(r8) (:)]  (gN/m2) dead stem N                               
   deadstemn_storage                   =>    pns%deadstemn_storage                       , & ! InOut:  [real(r8) (:)]  (gN/m2) dead stem N storage                       
   deadstemn_xfer                      =>    pns%deadstemn_xfer                          , & ! InOut:  [real(r8) (:)]  (gN/m2) dead stem N transfer                      
   frootn                              =>    pns%frootn                                  , & ! InOut:  [real(r8) (:)]  (gN/m2) fine root N                               
   frootn_storage                      =>    pns%frootn_storage                          , & ! InOut:  [real(r8) (:)]  (gN/m2) fine root N storage                       
   frootn_xfer                         =>    pns%frootn_xfer                             , & ! InOut:  [real(r8) (:)]  (gN/m2) fine root N transfer                      
   leafn                               =>    pns%leafn                                   , & ! InOut:  [real(r8) (:)]  (gN/m2) leaf N                                    
   leafn_storage                       =>    pns%leafn_storage                           , & ! InOut:  [real(r8) (:)]  (gN/m2) leaf N storage                            
   leafn_xfer                          =>    pns%leafn_xfer                              , & ! InOut:  [real(r8) (:)]  (gN/m2) leaf N transfer                           
   livecrootn                          =>    pns%livecrootn                              , & ! InOut:  [real(r8) (:)]  (gN/m2) live coarse root N                        
   livecrootn_storage                  =>    pns%livecrootn_storage                      , & ! InOut:  [real(r8) (:)]  (gN/m2) live coarse root N storage                
   livecrootn_xfer                     =>    pns%livecrootn_xfer                         , & ! InOut:  [real(r8) (:)]  (gN/m2) live coarse root N transfer               
   livestemn                           =>    pns%livestemn                               , & ! InOut:  [real(r8) (:)]  (gN/m2) live stem N                               
   livestemn_storage                   =>    pns%livestemn_storage                       , & ! InOut:  [real(r8) (:)]  (gN/m2) live stem N storage                       
   livestemn_xfer                      =>    pns%livestemn_xfer                          , & ! InOut:  [real(r8) (:)]  (gN/m2) live stem N transfer                      
   retransn                            =>    pns%retransn                                  & ! InOut:  [real(r8) (:)]  (gN/m2) plant pool of retranslocated N            
   )

   ! set time steps
   dt = real( get_step_size(), r8 )

   do j = 1, nlevdecomp
      ! column loop
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         
         ! column-level nitrogen fluxes from gap-phase mortality
         decomp_npools_vr(c,j,i_met_lit) = decomp_npools_vr(c,j,i_met_lit) + gap_mortality_n_to_litr_met_n(c,j) * dt
         decomp_npools_vr(c,j,i_cel_lit) = decomp_npools_vr(c,j,i_cel_lit) + gap_mortality_n_to_litr_cel_n(c,j) * dt
         decomp_npools_vr(c,j,i_lig_lit) = decomp_npools_vr(c,j,i_lig_lit) + gap_mortality_n_to_litr_lig_n(c,j) * dt
         decomp_npools_vr(c,j,i_cwd) = decomp_npools_vr(c,j,i_cwd) + gap_mortality_n_to_cwdn(c,j)  * dt

      end do ! end of column loop
   end do

   ! pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)

      ! pft-level nitrogen fluxes from gap-phase mortality
      ! displayed pools
      leafn(p)      = leafn(p)      - m_leafn_to_litter(p)      * dt
      frootn(p)     = frootn(p)     - m_frootn_to_litter(p)     * dt
      livestemn(p)  = livestemn(p)  - m_livestemn_to_litter(p)  * dt
      deadstemn(p)  = deadstemn(p)  - m_deadstemn_to_litter(p)  * dt
      livecrootn(p) = livecrootn(p) - m_livecrootn_to_litter(p) * dt
      deadcrootn(p) = deadcrootn(p) - m_deadcrootn_to_litter(p) * dt
      retransn(p)   = retransn(p)   - m_retransn_to_litter(p)   * dt

      ! storage pools
      leafn_storage(p)      = leafn_storage(p)      - m_leafn_storage_to_litter(p)      * dt
      frootn_storage(p)     = frootn_storage(p)     - m_frootn_storage_to_litter(p)     * dt
      livestemn_storage(p)  = livestemn_storage(p)  - m_livestemn_storage_to_litter(p)  * dt
      deadstemn_storage(p)  = deadstemn_storage(p)  - m_deadstemn_storage_to_litter(p)  * dt
      livecrootn_storage(p) = livecrootn_storage(p) - m_livecrootn_storage_to_litter(p) * dt
      deadcrootn_storage(p) = deadcrootn_storage(p) - m_deadcrootn_storage_to_litter(p) * dt

      ! transfer pools
      leafn_xfer(p)      = leafn_xfer(p)      - m_leafn_xfer_to_litter(p)      * dt
      frootn_xfer(p)     = frootn_xfer(p)     - m_frootn_xfer_to_litter(p)     * dt
      livestemn_xfer(p)  = livestemn_xfer(p)  - m_livestemn_xfer_to_litter(p)  * dt
      deadstemn_xfer(p)  = deadstemn_xfer(p)  - m_deadstemn_xfer_to_litter(p)  * dt
      livecrootn_xfer(p) = livecrootn_xfer(p) - m_livecrootn_xfer_to_litter(p) * dt
      deadcrootn_xfer(p) = deadcrootn_xfer(p) - m_deadcrootn_xfer_to_litter(p) * dt

   end do

    end associate 
 end subroutine NStateUpdate2
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
subroutine NStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp)
  !
  ! !DESCRIPTION:
  ! Update all the prognostic nitrogen state
  ! variables affected by harvest mortality fluxes
  !
  ! !USES:
  use clmtype
  use clm_time_manager, only: get_step_size
  use clm_varpar   , only: i_met_lit, i_cel_lit, i_lig_lit, i_cwd
  !
  ! !ARGUMENTS:
  implicit none
  integer, intent(in) :: num_soilc       ! number of soil columns in filter
  integer, intent(in) :: filter_soilc(:) ! filter for soil columns
  integer, intent(in) :: num_soilp       ! number of soil pfts in filter
  integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
  !
  ! !LOCAL VARIABLES:
  integer :: c,p,j,l         ! indices
  integer :: fp,fc       ! lake filter indices
  real(r8):: dt          ! radiation time step (seconds)
!-----------------------------------------------------------------------

   associate(& 
   harvest_n_to_litr_met_n             =>    cnf%harvest_n_to_litr_met_n                 , & ! Input:  [real(r8) (:,:)]  N fluxes associated with harvest to litter metabolic pool (gN/m3/s)
   harvest_n_to_litr_cel_n             =>    cnf%harvest_n_to_litr_cel_n                 , & ! Input:  [real(r8) (:,:)]  N fluxes associated with harvest to litter cellulose pool (gN/m3/s)
   harvest_n_to_litr_lig_n             =>    cnf%harvest_n_to_litr_lig_n                 , & ! Input:  [real(r8) (:,:)]  N fluxes associated with harvest to litter lignin pool (gN/m3/s)
   harvest_n_to_cwdn                   =>    cnf%harvest_n_to_cwdn                       , & ! Input:  [real(r8) (:,:)]  N fluxes associated with harvest to CWD pool (gN/m3/s)
   decomp_npools_vr                    =>    cns%decomp_npools_vr                        , & ! InOut:  [real(r8) (:,:,:)]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools

   hrv_deadcrootn_storage_to_litter    =>    pnf%hrv_deadcrootn_storage_to_litter        , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadcrootn_to_litter            =>    pnf%hrv_deadcrootn_to_litter                , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadcrootn_xfer_to_litter       =>    pnf%hrv_deadcrootn_xfer_to_litter           , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadstemn_storage_to_litter     =>    pnf%hrv_deadstemn_storage_to_litter         , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadstemn_to_prod10n            =>    pnf%hrv_deadstemn_to_prod10n                , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadstemn_to_prod100n           =>    pnf%hrv_deadstemn_to_prod100n               , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadstemn_xfer_to_litter        =>    pnf%hrv_deadstemn_xfer_to_litter            , & ! Input:  [real(r8) (:)]                                                    
   hrv_frootn_storage_to_litter        =>    pnf%hrv_frootn_storage_to_litter            , & ! Input:  [real(r8) (:)]                                                    
   hrv_frootn_to_litter                =>    pnf%hrv_frootn_to_litter                    , & ! Input:  [real(r8) (:)]                                                    
   hrv_frootn_xfer_to_litter           =>    pnf%hrv_frootn_xfer_to_litter               , & ! Input:  [real(r8) (:)]                                                    
   hrv_leafn_storage_to_litter         =>    pnf%hrv_leafn_storage_to_litter             , & ! Input:  [real(r8) (:)]                                                    
   hrv_leafn_to_litter                 =>    pnf%hrv_leafn_to_litter                     , & ! Input:  [real(r8) (:)]                                                    
   hrv_leafn_xfer_to_litter            =>    pnf%hrv_leafn_xfer_to_litter                , & ! Input:  [real(r8) (:)]                                                    
   hrv_livecrootn_storage_to_litter    =>    pnf%hrv_livecrootn_storage_to_litter        , & ! Input:  [real(r8) (:)]                                                    
   hrv_livecrootn_to_litter            =>    pnf%hrv_livecrootn_to_litter                , & ! Input:  [real(r8) (:)]                                                    
   hrv_livecrootn_xfer_to_litter       =>    pnf%hrv_livecrootn_xfer_to_litter           , & ! Input:  [real(r8) (:)]                                                    
   hrv_livestemn_storage_to_litter     =>    pnf%hrv_livestemn_storage_to_litter         , & ! Input:  [real(r8) (:)]                                                    
   hrv_livestemn_to_litter             =>    pnf%hrv_livestemn_to_litter                 , & ! Input:  [real(r8) (:)]                                                    
   hrv_livestemn_xfer_to_litter        =>    pnf%hrv_livestemn_xfer_to_litter            , & ! Input:  [real(r8) (:)]                                                    
   hrv_retransn_to_litter              =>    pnf%hrv_retransn_to_litter                  , & ! Input:  [real(r8) (:)]                                                    
   deadcrootn                          =>    pns%deadcrootn                              , & ! InOut:  [real(r8) (:)]  (gN/m2) dead coarse root N                        
   deadcrootn_storage                  =>    pns%deadcrootn_storage                      , & ! InOut:  [real(r8) (:)]  (gN/m2) dead coarse root N storage                
   deadcrootn_xfer                     =>    pns%deadcrootn_xfer                         , & ! InOut:  [real(r8) (:)]  (gN/m2) dead coarse root N transfer               
   deadstemn                           =>    pns%deadstemn                               , & ! InOut:  [real(r8) (:)]  (gN/m2) dead stem N                               
   deadstemn_storage                   =>    pns%deadstemn_storage                       , & ! InOut:  [real(r8) (:)]  (gN/m2) dead stem N storage                       
   deadstemn_xfer                      =>    pns%deadstemn_xfer                          , & ! InOut:  [real(r8) (:)]  (gN/m2) dead stem N transfer                      
   frootn                              =>    pns%frootn                                  , & ! InOut:  [real(r8) (:)]  (gN/m2) fine root N                               
   frootn_storage                      =>    pns%frootn_storage                          , & ! InOut:  [real(r8) (:)]  (gN/m2) fine root N storage                       
   frootn_xfer                         =>    pns%frootn_xfer                             , & ! InOut:  [real(r8) (:)]  (gN/m2) fine root N transfer                      
   leafn                               =>    pns%leafn                                   , & ! InOut:  [real(r8) (:)]  (gN/m2) leaf N                                    
   leafn_storage                       =>    pns%leafn_storage                           , & ! InOut:  [real(r8) (:)]  (gN/m2) leaf N storage                            
   leafn_xfer                          =>    pns%leafn_xfer                              , & ! InOut:  [real(r8) (:)]  (gN/m2) leaf N transfer                           
   livecrootn                          =>    pns%livecrootn                              , & ! InOut:  [real(r8) (:)]  (gN/m2) live coarse root N                        
   livecrootn_storage                  =>    pns%livecrootn_storage                      , & ! InOut:  [real(r8) (:)]  (gN/m2) live coarse root N storage                
   livecrootn_xfer                     =>    pns%livecrootn_xfer                         , & ! InOut:  [real(r8) (:)]  (gN/m2) live coarse root N transfer               
   livestemn                           =>    pns%livestemn                               , & ! InOut:  [real(r8) (:)]  (gN/m2) live stem N                               
   livestemn_storage                   =>    pns%livestemn_storage                       , & ! InOut:  [real(r8) (:)]  (gN/m2) live stem N storage                       
   livestemn_xfer                      =>    pns%livestemn_xfer                          , & ! InOut:  [real(r8) (:)]  (gN/m2) live stem N transfer                      
   retransn                            =>    pns%retransn                                  & ! InOut:  [real(r8) (:)]  (gN/m2) plant pool of retranslocated N            
   )

   ! set time steps
   dt = real( get_step_size(), r8 )

   do j = 1,nlevdecomp
      ! column loop
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         
         ! column-level nitrogen fluxes from harvest mortality
         decomp_npools_vr(c,j,i_met_lit) = decomp_npools_vr(c,j,i_met_lit) + harvest_n_to_litr_met_n(c,j) * dt
         decomp_npools_vr(c,j,i_cel_lit) = decomp_npools_vr(c,j,i_cel_lit) + harvest_n_to_litr_cel_n(c,j) * dt
         decomp_npools_vr(c,j,i_lig_lit) = decomp_npools_vr(c,j,i_lig_lit) + harvest_n_to_litr_lig_n(c,j) * dt
         decomp_npools_vr(c,j,i_cwd) = decomp_npools_vr(c,j,i_cwd) + harvest_n_to_cwdn(c,j)  * dt
         
      end do ! end of column loop
   end do

   ! pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)

      ! pft-level nitrogen fluxes from harvest mortality
      ! displayed pools
      leafn(p)      = leafn(p)      - hrv_leafn_to_litter(p)      * dt
      frootn(p)     = frootn(p)     - hrv_frootn_to_litter(p)     * dt
      livestemn(p)  = livestemn(p)  - hrv_livestemn_to_litter(p)  * dt
      deadstemn(p)  = deadstemn(p)  - hrv_deadstemn_to_prod10n(p) * dt
      deadstemn(p)  = deadstemn(p)  - hrv_deadstemn_to_prod100n(p)* dt
      livecrootn(p) = livecrootn(p) - hrv_livecrootn_to_litter(p) * dt
      deadcrootn(p) = deadcrootn(p) - hrv_deadcrootn_to_litter(p) * dt
      retransn(p)   = retransn(p)   - hrv_retransn_to_litter(p)   * dt

      ! storage pools
      leafn_storage(p)      = leafn_storage(p)      - hrv_leafn_storage_to_litter(p)      * dt
      frootn_storage(p)     = frootn_storage(p)     - hrv_frootn_storage_to_litter(p)     * dt
      livestemn_storage(p)  = livestemn_storage(p)  - hrv_livestemn_storage_to_litter(p)  * dt
      deadstemn_storage(p)  = deadstemn_storage(p)  - hrv_deadstemn_storage_to_litter(p)  * dt
      livecrootn_storage(p) = livecrootn_storage(p) - hrv_livecrootn_storage_to_litter(p) * dt
      deadcrootn_storage(p) = deadcrootn_storage(p) - hrv_deadcrootn_storage_to_litter(p) * dt

      ! transfer pools
      leafn_xfer(p)      = leafn_xfer(p)      - hrv_leafn_xfer_to_litter(p)      * dt
      frootn_xfer(p)     = frootn_xfer(p)     - hrv_frootn_xfer_to_litter(p)     * dt
      livestemn_xfer(p)  = livestemn_xfer(p)  - hrv_livestemn_xfer_to_litter(p)  * dt
      deadstemn_xfer(p)  = deadstemn_xfer(p)  - hrv_deadstemn_xfer_to_litter(p)  * dt
      livecrootn_xfer(p) = livecrootn_xfer(p) - hrv_livecrootn_xfer_to_litter(p) * dt
      deadcrootn_xfer(p) = deadcrootn_xfer(p) - hrv_deadcrootn_xfer_to_litter(p) * dt

   end do

    end associate 
 end subroutine NStateUpdate2h

end module CNNStateUpdate2Mod
