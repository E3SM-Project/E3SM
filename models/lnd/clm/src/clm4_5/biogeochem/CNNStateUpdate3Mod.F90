module CNNStateUpdate3Mod
!-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for nitrogen state variable update, mortality fluxes.
  ! Also, sminn leaching flux.
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varpar   , only: nlevdecomp, ndecomp_pools
  implicit none
  save
  private
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: NStateUpdate3
  !
  ! !REVISION HISTORY:
  ! 7/27/2004: Created by Peter Thornton
  ! F. Li and S. Levis (11/06/12)
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
  subroutine NStateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic nitrogen state
    ! variables affected by gap-phase mortality fluxes. Also the Sminn leaching flux.
    !
    ! !USES:
    use clmtype
    use clm_time_manager, only: get_step_size
    use clm_varctl      , only: iulog, use_nitrif_denitrif
    use clm_varpar      , only: i_cwd, i_met_lit, i_cel_lit, i_lig_lit
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_soilc       ! number of soil columns in filter
    integer, intent(in) :: filter_soilc(:) ! filter for soil columns
    integer, intent(in) :: num_soilp       ! number of soil pfts in filter
    integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,l,k        ! indices
    integer :: fp,fc      ! lake filter indices
    real(r8):: dt         ! radiation time step (seconds)
!-----------------------------------------------------------------------

   associate(& 
   fire_mortality_n_to_cwdn            =>    cnf%fire_mortality_n_to_cwdn                , & ! Input:  [real(r8) (:,:)]  N fluxes associated with fire mortality to CWD pool (gN/m3/s)
   sminn_leached_vr                    =>    cnf%sminn_leached_vr                        , & ! Input:  [real(r8) (:,:)]                                                  
   smin_no3_leached_vr                 =>    cnf%smin_no3_leached_vr                     , & ! Input:  [real(r8) (:,:)]                                                  
   smin_no3_runoff_vr                  =>    cnf%smin_no3_runoff_vr                      , & ! Input:  [real(r8) (:,:)]  vertically-resolved rate of mineral NO3 loss with runoff (gN/m3/s)
   smin_no3_vr                         =>    cns%smin_no3_vr                             , & ! Input:  [real(r8) (:,:)]                                                  
   smin_nh4_vr                         =>    cns%smin_nh4_vr                             , & ! Input:  [real(r8) (:,:)]                                                  
   m_decomp_npools_to_fire_vr          =>    cnf%m_decomp_npools_to_fire_vr              , & ! Input:  [real(r8) (:,:,:)]                                                
   m_n_to_litr_met_fire                =>    cnf%m_n_to_litr_met_fire                    , & ! Input:  [real(r8) (:,:)]                                                  
   m_n_to_litr_cel_fire                =>    cnf%m_n_to_litr_cel_fire                    , & ! Input:  [real(r8) (:,:)]                                                  
   m_n_to_litr_lig_fire                =>    cnf%m_n_to_litr_lig_fire                    , & ! Input:  [real(r8) (:,:)]                                                  
   decomp_npools_vr                    =>    cns%decomp_npools_vr                        , & ! InOut:  [real(r8) (:,:,:)]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
   sminn_vr                            =>    cns%sminn_vr                                , & ! InOut:  [real(r8) (:,:)]  (gN/m3) soil mineral N                          
   m_leafn_to_fire                     =>    pnf%m_leafn_to_fire                         , & ! Input:  [real(r8) (:)]                                                    
   m_leafn_storage_to_fire             =>    pnf%m_leafn_storage_to_fire                 , & ! Input:  [real(r8) (:)]                                                    
   m_leafn_xfer_to_fire                =>    pnf%m_leafn_xfer_to_fire                    , & ! Input:  [real(r8) (:)]                                                    
   m_livestemn_to_fire                 =>    pnf%m_livestemn_to_fire                     , & ! Input:  [real(r8) (:)]                                                    
   m_livestemn_storage_to_fire         =>    pnf%m_livestemn_storage_to_fire             , & ! Input:  [real(r8) (:)]                                                    
   m_livestemn_xfer_to_fire            =>    pnf%m_livestemn_xfer_to_fire                , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemn_to_fire                 =>    pnf%m_deadstemn_to_fire                     , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemn_storage_to_fire         =>    pnf%m_deadstemn_storage_to_fire             , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemn_xfer_to_fire            =>    pnf%m_deadstemn_xfer_to_fire                , & ! Input:  [real(r8) (:)]                                                    
   m_frootn_to_fire                    =>    pnf%m_frootn_to_fire                        , & ! Input:  [real(r8) (:)]                                                    
   m_frootn_storage_to_fire            =>    pnf%m_frootn_storage_to_fire                , & ! Input:  [real(r8) (:)]                                                    
   m_frootn_xfer_to_fire               =>    pnf%m_frootn_xfer_to_fire                   , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootn_to_fire                =>    pnf%m_livecrootn_to_fire                    , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootn_storage_to_fire        =>    pnf%m_livecrootn_storage_to_fire            , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootn_xfer_to_fire           =>    pnf%m_livecrootn_xfer_to_fire               , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootn_to_fire                =>    pnf%m_deadcrootn_to_fire                    , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootn_storage_to_fire        =>    pnf%m_deadcrootn_storage_to_fire            , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootn_xfer_to_fire           =>    pnf%m_deadcrootn_xfer_to_fire               , & ! Input:  [real(r8) (:)]                                                    
   m_retransn_to_fire                  =>    pnf%m_retransn_to_fire                      , & ! Input:  [real(r8) (:)]                                                    
   m_leafn_to_litter_fire              =>    pnf%m_leafn_to_litter_fire                  , & ! Input:  [real(r8) (:)]                                                    
   m_leafn_storage_to_litter_fire      =>    pnf%m_leafn_storage_to_litter_fire          , & ! Input:  [real(r8) (:)]                                                    
   m_leafn_xfer_to_litter_fire         =>    pnf%m_leafn_xfer_to_litter_fire             , & ! Input:  [real(r8) (:)]                                                    
   m_livestemn_to_litter_fire          =>    pnf%m_livestemn_to_litter_fire              , & ! Input:  [real(r8) (:)]                                                    
   m_livestemn_storage_to_litter_fire  =>    pnf%m_livestemn_storage_to_litter_fire      , & ! Input:  [real(r8) (:)]                                                    
   m_livestemn_xfer_to_litter_fire     =>    pnf%m_livestemn_xfer_to_litter_fire         , & ! Input:  [real(r8) (:)]                                                    
   m_livestemn_to_deadstemn_fire       =>    pnf%m_livestemn_to_deadstemn_fire           , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemn_to_litter_fire          =>    pnf%m_deadstemn_to_litter_fire              , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemn_storage_to_litter_fire  =>    pnf%m_deadstemn_storage_to_litter_fire      , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemn_xfer_to_litter_fire     =>   pnf%m_deadstemn_xfer_to_litter_fire          , & ! Input:  [real(r8) (:)]                                                    
   m_frootn_to_litter_fire             =>    pnf%m_frootn_to_litter_fire                 , & ! Input:  [real(r8) (:)]                                                    
   m_frootn_storage_to_litter_fire     =>    pnf%m_frootn_storage_to_litter_fire         , & ! Input:  [real(r8) (:)]                                                    
   m_frootn_xfer_to_litter_fire        =>    pnf%m_frootn_xfer_to_litter_fire            , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootn_to_litter_fire         =>    pnf%m_livecrootn_to_litter_fire             , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootn_storage_to_litter_fire  =>    pnf%m_livecrootn_storage_to_litter_fire     , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootn_xfer_to_litter_fire    =>    pnf%m_livecrootn_xfer_to_litter_fire        , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootn_to_deadcrootn_fire     =>    pnf%m_livecrootn_to_deadcrootn_fire         , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootn_to_litter_fire         =>    pnf%m_deadcrootn_to_litter_fire             , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootn_storage_to_litter_fire  =>    pnf%m_deadcrootn_storage_to_litter_fire     , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootn_xfer_to_litter_fire    =>    pnf%m_deadcrootn_xfer_to_litter_fire        , & ! Input:  [real(r8) (:)]                                                    
   m_retransn_to_litter_fire           =>    pnf%m_retransn_to_litter_fire               , & ! Input:  [real(r8) (:)]                                                    
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

         if (.not. use_nitrif_denitrif) then
            ! mineral N loss due to leaching
            sminn_vr(c,j) = sminn_vr(c,j) - sminn_leached_vr(c,j) * dt
         else
            ! mineral N loss due to leaching and runoff
            smin_no3_vr(c,j) = max(smin_no3_vr(c,j) - ( smin_no3_leached_vr(c,j) + smin_no3_runoff_vr(c,j) ) * dt, 0._r8)
            sminn_vr(c,j) = smin_no3_vr(c,j) + smin_nh4_vr(c,j)
         end if
         
         ! column level nitrogen fluxes from fire
         ! pft-level wood to column-level CWD (uncombusted wood)
         decomp_npools_vr(c,j,i_cwd) = decomp_npools_vr(c,j,i_cwd) + fire_mortality_n_to_cwdn(c,j) * dt
         
         ! pft-level wood to column-level litter (uncombusted wood)
         decomp_npools_vr(c,j,i_met_lit) = decomp_npools_vr(c,j,i_met_lit) + m_n_to_litr_met_fire(c,j)* dt
         decomp_npools_vr(c,j,i_cel_lit) = decomp_npools_vr(c,j,i_cel_lit) + m_n_to_litr_cel_fire(c,j)* dt
         decomp_npools_vr(c,j,i_lig_lit) = decomp_npools_vr(c,j,i_lig_lit) + m_n_to_litr_lig_fire(c,j)* dt
      end do ! end of column loop
   end do
   
   ! litter and CWD losses to fire
   do l = 1, ndecomp_pools
      do j = 1, nlevdecomp
         ! column loop
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            decomp_npools_vr(c,j,l) = decomp_npools_vr(c,j,l) - m_decomp_npools_to_fire_vr(c,j,l) * dt
         end do
      end do
   end do

   ! pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)

      ! pft-level nitrogen fluxes from fire
      ! displayed pools
      leafn(p)      = leafn(p)      - m_leafn_to_fire(p)             * dt
      frootn(p)     = frootn(p)     - m_frootn_to_fire(p)            * dt
      livestemn(p)  = livestemn(p)  - m_livestemn_to_fire(p)         * dt
      deadstemn(p)  = deadstemn(p)  - m_deadstemn_to_fire(p)         * dt
      livecrootn(p) = livecrootn(p) - m_livecrootn_to_fire(p)        * dt
      deadcrootn(p) = deadcrootn(p) - m_deadcrootn_to_fire(p)        * dt
      
      leafn(p)      = leafn(p)      - m_leafn_to_litter_fire(p)           * dt
      frootn(p)     = frootn(p)     - m_frootn_to_litter_fire(p)          * dt
      livestemn(p)  = livestemn(p)  - m_livestemn_to_litter_fire(p)       * dt
      deadstemn(p)  = deadstemn(p)  - m_deadstemn_to_litter_fire(p)       * dt
      livecrootn(p) = livecrootn(p) - m_livecrootn_to_litter_fire(p)      * dt
      deadcrootn(p) = deadcrootn(p) - m_deadcrootn_to_litter_fire(p)      * dt
     
      ! storage pools
      leafn_storage(p)      = leafn_storage(p)      - m_leafn_storage_to_fire(p)      * dt
      frootn_storage(p)     = frootn_storage(p)     - m_frootn_storage_to_fire(p)     * dt
      livestemn_storage(p)  = livestemn_storage(p)  - m_livestemn_storage_to_fire(p)  * dt
      deadstemn_storage(p)  = deadstemn_storage(p)  - m_deadstemn_storage_to_fire(p)  * dt
      livecrootn_storage(p) = livecrootn_storage(p) - m_livecrootn_storage_to_fire(p) * dt
      deadcrootn_storage(p) = deadcrootn_storage(p) - m_deadcrootn_storage_to_fire(p) * dt

      leafn_storage(p)      = leafn_storage(p)      - m_leafn_storage_to_litter_fire(p)      * dt
      frootn_storage(p)     = frootn_storage(p)     - m_frootn_storage_to_litter_fire(p)     * dt
      livestemn_storage(p)  = livestemn_storage(p)  - m_livestemn_storage_to_litter_fire(p)  * dt
      deadstemn_storage(p)  = deadstemn_storage(p)  - m_deadstemn_storage_to_litter_fire(p)  * dt
      livecrootn_storage(p) = livecrootn_storage(p) - m_livecrootn_storage_to_litter_fire(p) * dt
      deadcrootn_storage(p) = deadcrootn_storage(p) - m_deadcrootn_storage_to_litter_fire(p) * dt


      ! transfer pools
      leafn_xfer(p)      = leafn_xfer(p)      - m_leafn_xfer_to_fire(p)      * dt
      frootn_xfer(p)     = frootn_xfer(p)     - m_frootn_xfer_to_fire(p)     * dt
      livestemn_xfer(p)  = livestemn_xfer(p)  - m_livestemn_xfer_to_fire(p)  * dt
      deadstemn_xfer(p)  = deadstemn_xfer(p)  - m_deadstemn_xfer_to_fire(p)  * dt
      livecrootn_xfer(p) = livecrootn_xfer(p) - m_livecrootn_xfer_to_fire(p) * dt
      deadcrootn_xfer(p) = deadcrootn_xfer(p) - m_deadcrootn_xfer_to_fire(p) * dt
  
      leafn_xfer(p)      = leafn_xfer(p)      - m_leafn_xfer_to_litter_fire(p)      * dt
      frootn_xfer(p)     = frootn_xfer(p)     - m_frootn_xfer_to_litter_fire(p)     * dt
      livestemn_xfer(p)  = livestemn_xfer(p)  - m_livestemn_xfer_to_litter_fire(p)  * dt
      deadstemn_xfer(p)  = deadstemn_xfer(p)  - m_deadstemn_xfer_to_litter_fire(p)  * dt
      livecrootn_xfer(p) = livecrootn_xfer(p) - m_livecrootn_xfer_to_litter_fire(p) * dt
      deadcrootn_xfer(p) = deadcrootn_xfer(p) - m_deadcrootn_xfer_to_litter_fire(p) * dt

      ! retranslocated N pool
      retransn(p) = retransn(p) - m_retransn_to_fire(p)        * dt
      retransn(p) = retransn(p) - m_retransn_to_litter_fire(p) * dt

   end do

    end associate 
 end subroutine NStateUpdate3

end module CNNStateUpdate3Mod
