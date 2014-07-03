module CNNStateUpdate1Mod
!-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for nitrogen state variable updates, non-mortality fluxes.
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  save
  private
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: NStateUpdate1
  !
  ! !REVISION HISTORY:
  ! 4/23/2004: Created by Peter Thornton
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
subroutine NStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp)
  !
  ! !DESCRIPTION:
  ! On the radiation time step, update all the prognostic nitrogen state
  ! variables (except for gap-phase mortality and fire fluxes)
  !
  ! !USES:
  use clmtype
  use clm_time_manager, only: get_step_size
  use clm_varpar      , only: nlevdecomp, ndecomp_pools, ndecomp_cascade_transitions, &
                              crop_prog, i_met_lit, i_cel_lit, i_lig_lit, i_cwd
  use clm_varctl      , only: iulog, use_nitrif_denitrif
  use clm_varcon      , only: nitrif_n2o_loss_frac
  use pftvarcon       , only: npcropmin, nc3crop
  !
  ! !ARGUMENTS:
  implicit none
  integer, intent(in) :: num_soilc       ! number of soil columns in filter
  integer, intent(in) :: filter_soilc(:) ! filter for soil columns
  integer, intent(in) :: num_soilp       ! number of soil pfts in filter
  integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
  !
  ! !LOCAL VARIABLES:
  integer :: c,p,j,l,k      ! indices
  integer :: fp,fc    ! lake filter indices
  real(r8):: dt       ! radiation time step (seconds)
!-----------------------------------------------------------------------

   associate(& 
   woody                               =>    pftcon%woody                                , & ! Input:  [real(r8) (:)]  binary flag for woody lifeform (1=woody, 0=not woody)
   ndep_to_sminn                       =>    cnf%ndep_to_sminn                           , & ! Input:  [real(r8) (:)]                                                    
   nfix_to_sminn                       =>    cnf%nfix_to_sminn                           , & ! Input:  [real(r8) (:)]  symbiotic/asymbiotic N fixation to soil mineral N (gN/m2/s)
   fert_to_sminn                       =>    cnf%fert_to_sminn                           , & ! Input:  [real(r8) (:)]                                                    
   soyfixn_to_sminn                    =>    cnf%soyfixn_to_sminn                        , & ! Input:  [real(r8) (:)]                                                    
   sminn_to_denit_excess_vr            =>    cnf%sminn_to_denit_excess_vr                , & ! Input:  [real(r8) (:,:)]                                                  
   sminn_to_denit_decomp_cascade_vr    =>    cnf%sminn_to_denit_decomp_cascade_vr        , & ! Input:  [real(r8) (:,:,:)]  vertically-resolved denitrification along decomp cascade (gN/m3/s)
   smin_no3_vr                         =>    cns%smin_no3_vr                             , & ! InOut:  [real(r8) (:,:)]  (gN/m3) soil NO3                                
   smin_nh4_vr                         =>    cns%smin_nh4_vr                             , & ! InOut:  [real(r8) (:,:)]  (gN/m3) soil NH4                                
   f_nit_vr                            =>    cnf%f_nit_vr                                , & ! InOut:  [real(r8) (:,:)]  (gN/m3/s) soil nitrification flux               
   f_denit_vr                          =>    cnf%f_denit_vr                              , & ! InOut:  [real(r8) (:,:)]  (gN/m3/s) soil denitrification flux             
   actual_immob_no3_vr                 =>    cnf%actual_immob_no3_vr                     , & ! InOut:  [real(r8) (:,:)]  (gN/m3/s)                                       
   actual_immob_nh4_vr                 =>    cnf%actual_immob_nh4_vr                     , & ! InOut:  [real(r8) (:,:)]  (gN/m3/s)                                       
   smin_no3_to_plant_vr                =>    cnf%smin_no3_to_plant_vr                    , & ! InOut:  [real(r8) (:,:)]  (gN/m3/s)                                       
   smin_nh4_to_plant_vr                =>    cnf%smin_nh4_to_plant_vr                    , & ! InOut:  [real(r8) (:,:)]  (gN/m3/s)                                       
   gross_nmin_vr                       =>    cnf%gross_nmin_vr                           , & ! InOut:  [real(r8) (:,:)]  (gN/m3/s)                                       
   sminn_to_plant_vr                   =>    cnf%sminn_to_plant_vr                       , & ! Input:  [real(r8) (:,:)]                                                  
   decomp_cascade_sminn_flux_vr        =>    cnf%decomp_cascade_sminn_flux_vr            , & ! InOut:  [real(r8) (:,:,:)]  vert-res mineral N flux for transition along decomposition cascade (gN/m3/s)
   decomp_cascade_ntransfer_vr         =>    cnf%decomp_cascade_ntransfer_vr             , & ! InOut:  [real(r8) (:,:,:)]  vert-res transfer of N from donor to receiver pool along decomp. cascade (gN/m3/s)
   cascade_donor_pool                  =>    decomp_cascade_con%cascade_donor_pool       , & ! InOut:  [integer (:)]  which pool is C taken from for a given decomposition step
   cascade_receiver_pool               =>    decomp_cascade_con%cascade_receiver_pool    , & ! InOut:  [integer (:)]  which pool is C added to for a given decomposition step
   supplement_to_sminn_vr              =>    cnf%supplement_to_sminn_vr                  , & ! Input:  [real(r8) (:,:)]                                                  
   decomp_npools_vr                    =>    cns%decomp_npools_vr                        , & ! InOut:  [real(r8) (:,:,:)]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
   decomp_npools_sourcesink            =>    cnf%decomp_npools_sourcesink                , & ! InOut:  [real(r8) (:,:,:)]  (gC/m3)  change in decomposing N pools over a timestep.  Used to update concentrations concurrently with vertical transport
   sminn_vr                            =>    cns%sminn_vr                                , & ! InOut:  [real(r8) (:,:)]  (gN/m3) soil mineral N                          
   ndep_prof                           =>    cps%ndep_prof                               , & ! InOut:  [real(r8) (:,:)]  profile over which N deposition is distributed through column (1/m)
   nfixation_prof                      =>    cps%nfixation_prof                          , & ! InOut:  [real(r8) (:,:)]  profile over which N fixation is distributed through column (1/m)
   phenology_n_to_litr_met_n           =>    cnf%phenology_n_to_litr_met_n               , & ! InOut:  [real(r8) (:,:)]  N fluxes associated with phenology (litterfall and crop) to litter metabolic pool (gN/m3/s)
   phenology_n_to_litr_cel_n           =>    cnf%phenology_n_to_litr_cel_n               , & ! InOut:  [real(r8) (:,:)]  N fluxes associated with phenology (litterfall and crop) to litter cellulose pool (gN/m3/s)
   phenology_n_to_litr_lig_n           =>    cnf%phenology_n_to_litr_lig_n               , & ! InOut:  [real(r8) (:,:)]  N fluxes associated with phenology (litterfall and crop) to litter lignin pool (gN/m3/s)
   dwt_seedn_to_leaf                   =>    cnf%dwt_seedn_to_leaf                       , & ! InOut:  [real(r8) (:)]                                                    
   dwt_seedn_to_deadstem               =>    cnf%dwt_seedn_to_deadstem                   , & ! InOut:  [real(r8) (:)]                                                    
   dwt_frootn_to_litr_met_n            =>    cnf%dwt_frootn_to_litr_met_n                , & ! InOut:  [real(r8) (:,:)]                                                  
   dwt_frootn_to_litr_cel_n            =>    cnf%dwt_frootn_to_litr_cel_n                , & ! InOut:  [real(r8) (:,:)]                                                  
   dwt_frootn_to_litr_lig_n            =>    cnf%dwt_frootn_to_litr_lig_n                , & ! InOut:  [real(r8) (:,:)]                                                  
   dwt_livecrootn_to_cwdn              =>    cnf%dwt_livecrootn_to_cwdn                  , & ! InOut:  [real(r8) (:,:)]                                                  
   dwt_deadcrootn_to_cwdn              =>    cnf%dwt_deadcrootn_to_cwdn                  , & ! InOut:  [real(r8) (:,:)]                                                  
   seedn                               =>    cns%seedn                                   , & ! InOut:  [real(r8) (:)]                                                    
   ivt                                 =>    pft%itype                                   , & ! Input:  [integer (:)]  pft vegetation type                                
   deadcrootn_storage_to_xfer          =>    pnf%deadcrootn_storage_to_xfer              , & ! Input:  [real(r8) (:)]                                                    
   deadcrootn_xfer_to_deadcrootn       =>    pnf%deadcrootn_xfer_to_deadcrootn           , & ! Input:  [real(r8) (:)]                                                    
   deadstemn_storage_to_xfer           =>    pnf%deadstemn_storage_to_xfer               , & ! Input:  [real(r8) (:)]                                                    
   deadstemn_xfer_to_deadstemn         =>    pnf%deadstemn_xfer_to_deadstemn             , & ! Input:  [real(r8) (:)]                                                    
   frootn_storage_to_xfer              =>    pnf%frootn_storage_to_xfer                  , & ! Input:  [real(r8) (:)]                                                    
   frootn_to_litter                    =>    pnf%frootn_to_litter                        , & ! Input:  [real(r8) (:)]                                                    
   frootn_to_retransn                  =>    pnf%frootn_to_retransn                      , & ! Input:  [real(r8) (:)]                                                    
   frootn_xfer_to_frootn               =>    pnf%frootn_xfer_to_frootn                   , & ! Input:  [real(r8) (:)]                                                    
   leafn_storage_to_xfer               =>    pnf%leafn_storage_to_xfer                   , & ! Input:  [real(r8) (:)]                                                    
   leafn_to_litter                     =>    pnf%leafn_to_litter                         , & ! Input:  [real(r8) (:)]                                                    
   leafn_to_retransn                   =>    pnf%leafn_to_retransn                       , & ! Input:  [real(r8) (:)]                                                    
   leafn_xfer_to_leafn                 =>    pnf%leafn_xfer_to_leafn                     , & ! Input:  [real(r8) (:)]                                                    
   livecrootn_storage_to_xfer          =>    pnf%livecrootn_storage_to_xfer              , & ! Input:  [real(r8) (:)]                                                    
   livecrootn_to_deadcrootn            =>    pnf%livecrootn_to_deadcrootn                , & ! Input:  [real(r8) (:)]                                                    
   livecrootn_to_retransn              =>    pnf%livecrootn_to_retransn                  , & ! Input:  [real(r8) (:)]                                                    
   livecrootn_xfer_to_livecrootn       =>    pnf%livecrootn_xfer_to_livecrootn           , & ! Input:  [real(r8) (:)]                                                    
   livestemn_storage_to_xfer           =>    pnf%livestemn_storage_to_xfer               , & ! Input:  [real(r8) (:)]                                                    
   livestemn_to_deadstemn              =>    pnf%livestemn_to_deadstemn                  , & ! Input:  [real(r8) (:)]                                                    
   livestemn_to_retransn               =>    pnf%livestemn_to_retransn                   , & ! Input:  [real(r8) (:)]                                                    
   livestemn_xfer_to_livestemn         =>    pnf%livestemn_xfer_to_livestemn             , & ! Input:  [real(r8) (:)]                                                    
   npool_to_deadcrootn                 =>    pnf%npool_to_deadcrootn                     , & ! Input:  [real(r8) (:)]                                                    
   npool_to_deadcrootn_storage         =>    pnf%npool_to_deadcrootn_storage             , & ! Input:  [real(r8) (:)]                                                    
   npool_to_deadstemn                  =>    pnf%npool_to_deadstemn                      , & ! Input:  [real(r8) (:)]                                                    
   npool_to_deadstemn_storage          =>    pnf%npool_to_deadstemn_storage              , & ! Input:  [real(r8) (:)]                                                    
   npool_to_frootn                     =>    pnf%npool_to_frootn                         , & ! Input:  [real(r8) (:)]                                                    
   npool_to_frootn_storage             =>    pnf%npool_to_frootn_storage                 , & ! Input:  [real(r8) (:)]                                                    
   npool_to_leafn                      =>    pnf%npool_to_leafn                          , & ! Input:  [real(r8) (:)]                                                    
   npool_to_leafn_storage              =>    pnf%npool_to_leafn_storage                  , & ! Input:  [real(r8) (:)]                                                    
   npool_to_livecrootn                 =>    pnf%npool_to_livecrootn                     , & ! Input:  [real(r8) (:)]                                                    
   npool_to_livecrootn_storage         =>    pnf%npool_to_livecrootn_storage             , & ! Input:  [real(r8) (:)]                                                    
   npool_to_livestemn                  =>    pnf%npool_to_livestemn                      , & ! Input:  [real(r8) (:)]  allocation to live stem N (gN/m2/s)               
   npool_to_livestemn_storage          =>    pnf%npool_to_livestemn_storage              , & ! Input:  [real(r8) (:)]  allocation to live stem N storage (gN/m2/s)       
   retransn_to_npool                   =>    pnf%retransn_to_npool                       , & ! Input:  [real(r8) (:)]  deployment of retranslocated N (gN/m2/s)          
   sminn_to_npool                      =>    pnf%sminn_to_npool                          , & ! Input:  [real(r8) (:)]  deployment of soil mineral N uptake (gN/m2/s)     
   grainn_storage_to_xfer              =>    pnf%grainn_storage_to_xfer                  , & ! Input:  [real(r8) (:)]  grain N shift storage to transfer (gN/m2/s)       
   grainn_to_food                      =>    pnf%grainn_to_food                          , & ! Input:  [real(r8) (:)]  grain N to food (gN/m2/s)                         
   grainn_xfer_to_grainn               =>    pnf%grainn_xfer_to_grainn                   , & ! Input:  [real(r8) (:)]  grain N growth from storage (gN/m2/s)             
   livestemn_to_litter                 =>    pnf%livestemn_to_litter                     , & ! Input:  [real(r8) (:)]  livestem N to litter (gN/m2/s)                    
   npool_to_grainn                     =>    pnf%npool_to_grainn                         , & ! Input:  [real(r8) (:)]  allocation to grain N (gN/m2/s)                   
   npool_to_grainn_storage             =>    pnf%npool_to_grainn_storage                 , & ! Input:  [real(r8) (:)]  allocation to grain N storage (gN/m2/s)           
   grainn                              =>    pns%grainn                                  , & ! InOut:  [real(r8) (:)]  (gN/m2) grain N                                   
   grainn_storage                      =>    pns%grainn_storage                          , & ! InOut:  [real(r8) (:)]  (gN/m2) grain N storage                           
   grainn_xfer                         =>    pns%grainn_xfer                             , & ! InOut:  [real(r8) (:)]  (gN/m2) grain N transfer                          
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
   npool                               =>    pns%npool                                   , & ! InOut:  [real(r8) (:)]  (gN/m2) temporary plant N pool                    
   retransn                            =>    pns%retransn                                  & ! InOut:  [real(r8) (:)]  (gN/m2) plant pool of retranslocated N            
   )

   ! set time steps
   dt = real( get_step_size(), r8 )
   ! column-level fluxes
   
   ! column loop
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      ! seeding fluxes, from dynamic landcover
      seedn(c) = seedn(c) - dwt_seedn_to_leaf(c) * dt
      seedn(c) = seedn(c) - dwt_seedn_to_deadstem(c) * dt
   end do
   
   do j = 1, nlevdecomp
      ! column loop
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         
         if (.not. use_nitrif_denitrif) then
            ! N deposition and fixation
            sminn_vr(c,j) = sminn_vr(c,j) + ndep_to_sminn(c)*dt * ndep_prof(c,j)
            sminn_vr(c,j) = sminn_vr(c,j) + nfix_to_sminn(c)*dt * nfixation_prof(c,j)
         else
            ! N deposition and fixation (put all into NH4 pool)
            smin_nh4_vr(c,j) = smin_nh4_vr(c,j) + ndep_to_sminn(c)*dt * ndep_prof(c,j)
            smin_nh4_vr(c,j) = smin_nh4_vr(c,j) + nfix_to_sminn(c)*dt * nfixation_prof(c,j)
         end if

         ! plant to litter fluxes
         ! phenology and dynamic landcover fluxes
         decomp_npools_sourcesink(c,j,i_met_lit) = ( phenology_n_to_litr_met_n(c,j) + dwt_frootn_to_litr_met_n(c,j) ) *dt
         decomp_npools_sourcesink(c,j,i_cel_lit) = ( phenology_n_to_litr_cel_n(c,j) + dwt_frootn_to_litr_cel_n(c,j) ) *dt
         decomp_npools_sourcesink(c,j,i_lig_lit) = ( phenology_n_to_litr_lig_n(c,j) + dwt_frootn_to_litr_lig_n(c,j) ) *dt
         decomp_npools_sourcesink(c,j,i_cwd)	=  ( dwt_livecrootn_to_cwdn(c,j) + dwt_deadcrootn_to_cwdn(c,j) )*dt

      end do
   end do
   
   ! repeating N dep and fixation for crops
   if ( crop_prog )then
      do j = 1, nlevdecomp
         ! column loop
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            if (.not. use_nitrif_denitrif) then
               ! N deposition and fixation
               sminn_vr(c,j) = sminn_vr(c,j) + fert_to_sminn(c)*dt * ndep_prof(c,j)
               sminn_vr(c,j) = sminn_vr(c,j) + soyfixn_to_sminn(c)*dt * nfixation_prof(c,j)
            else
               ! N deposition and fixation (put all into NH4 pool)
               smin_nh4_vr(c,j) = smin_nh4_vr(c,j) + fert_to_sminn(c)*dt * ndep_prof(c,j)
               smin_nh4_vr(c,j) = smin_nh4_vr(c,j) + soyfixn_to_sminn(c)*dt * nfixation_prof(c,j)
            end if
         end do
      end do
   end if

   ! decomposition fluxes
   do k = 1, ndecomp_cascade_transitions
      do j = 1, nlevdecomp
         ! column loop
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            decomp_npools_sourcesink(c,j,cascade_donor_pool(k)) = decomp_npools_sourcesink(c,j,cascade_donor_pool(k)) - &
                 decomp_cascade_ntransfer_vr(c,j,k) * dt
         end do
      end do
   end do
   do k = 1, ndecomp_cascade_transitions
      if ( cascade_receiver_pool(k) .ne. 0 ) then  ! skip terminal transitions
         do j = 1, nlevdecomp
            ! column loop
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               decomp_npools_sourcesink(c,j,cascade_receiver_pool(k)) = decomp_npools_sourcesink(c,j,cascade_receiver_pool(k)) + &
                    (decomp_cascade_ntransfer_vr(c,j,k) + decomp_cascade_sminn_flux_vr(c,j,k)) * dt
            end do
         end do
      else  ! terminal transitions
         do j = 1, nlevdecomp
            ! column loop
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               decomp_npools_sourcesink(c,j,cascade_donor_pool(k)) = decomp_npools_sourcesink(c,j,cascade_donor_pool(k)) - &
                    decomp_cascade_sminn_flux_vr(c,j,k) * dt
            end do
         end do
      end if
   end do
   
   if (.not. use_nitrif_denitrif) then
      ! immobilization/mineralization in litter-to-SOM and SOM-to-SOM fluxes and denitrification fluxes
      do k = 1, ndecomp_cascade_transitions
         if ( cascade_receiver_pool(k) .ne. 0 ) then  ! skip terminal transitions
            do j = 1, nlevdecomp
               ! column loop
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  sminn_vr(c,j)  = sminn_vr(c,j) - &
                       (sminn_to_denit_decomp_cascade_vr(c,j,k) + decomp_cascade_sminn_flux_vr(c,j,k))* dt
               end do
            end do
         else
            do j = 1, nlevdecomp
               ! column loop
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  sminn_vr(c,j)  = sminn_vr(c,j) - sminn_to_denit_decomp_cascade_vr(c,j,k)* dt
                  sminn_vr(c,j)  = sminn_vr(c,j) + decomp_cascade_sminn_flux_vr(c,j,k)* dt

               end do
            end do
         endif
      end do

      do j = 1, nlevdecomp
         ! column loop
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            ! "bulk denitrification"
            sminn_vr(c,j) = sminn_vr(c,j) - sminn_to_denit_excess_vr(c,j) * dt

            ! total plant uptake from mineral N
            sminn_vr(c,j) = sminn_vr(c,j) - sminn_to_plant_vr(c,j)*dt

            ! flux that prevents N limitation (when Carbon_only is set)
            sminn_vr(c,j) = sminn_vr(c,j) + supplement_to_sminn_vr(c,j)*dt
         end do
      end do

   else   !-------------    NITRIF_DENITRIF --------------------

      do j = 1, nlevdecomp
         ! column loop
         do fc = 1,num_soilc
            c = filter_soilc(fc)

            ! mineralization fluxes (divert a fraction of this stream to nitrification flux, add the rest to NH4 pool)
            smin_nh4_vr(c,j) = smin_nh4_vr(c,j) + gross_nmin_vr(c,j)*dt

            ! immobilization fluxes
            smin_nh4_vr(c,j) = smin_nh4_vr(c,j) - actual_immob_nh4_vr(c,j)*dt
            smin_no3_vr(c,j) = smin_no3_vr(c,j) - actual_immob_no3_vr(c,j)*dt

            ! plant uptake fluxes
            smin_nh4_vr(c,j) = smin_nh4_vr(c,j) - smin_nh4_to_plant_vr(c,j)*dt
            smin_no3_vr(c,j) = smin_no3_vr(c,j) - smin_no3_to_plant_vr(c,j)*dt

            ! Account for nitrification fluxes
            smin_nh4_vr(c,j) = smin_nh4_vr(c,j) - f_nit_vr(c,j) * dt
            smin_no3_vr(c,j) = smin_no3_vr(c,j) + f_nit_vr(c,j) * dt * (1._r8 - nitrif_n2o_loss_frac)
            ! Account for denitrification fluxes
            smin_no3_vr(c,j) = smin_no3_vr(c,j) - f_denit_vr(c,j) * dt

            ! flux that prevents N limitation (when Carbon_only is set; put all into NH4)
            smin_nh4_vr(c,j) = smin_nh4_vr(c,j) + supplement_to_sminn_vr(c,j)*dt

            ! update diagnostic total
            sminn_vr(c,j) = smin_nh4_vr(c,j) + smin_no3_vr(c,j)

         end do ! end of column loop
      end do
   end if

   ! pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)

      ! phenology: transfer growth fluxes
      leafn(p)       = leafn(p)       + leafn_xfer_to_leafn(p)*dt
      leafn_xfer(p)  = leafn_xfer(p)  - leafn_xfer_to_leafn(p)*dt
      frootn(p)      = frootn(p)      + frootn_xfer_to_frootn(p)*dt
      frootn_xfer(p) = frootn_xfer(p) - frootn_xfer_to_frootn(p)*dt
      if (woody(ivt(p)) == 1.0_r8) then
          livestemn(p)       = livestemn(p)       + livestemn_xfer_to_livestemn(p)*dt
          livestemn_xfer(p)  = livestemn_xfer(p)  - livestemn_xfer_to_livestemn(p)*dt
          deadstemn(p)       = deadstemn(p)       + deadstemn_xfer_to_deadstemn(p)*dt
          deadstemn_xfer(p)  = deadstemn_xfer(p)  - deadstemn_xfer_to_deadstemn(p)*dt
          livecrootn(p)      = livecrootn(p)      + livecrootn_xfer_to_livecrootn(p)*dt
          livecrootn_xfer(p) = livecrootn_xfer(p) - livecrootn_xfer_to_livecrootn(p)*dt
          deadcrootn(p)      = deadcrootn(p)      + deadcrootn_xfer_to_deadcrootn(p)*dt
          deadcrootn_xfer(p) = deadcrootn_xfer(p) - deadcrootn_xfer_to_deadcrootn(p)*dt
      end if
      if (ivt(p) >= npcropmin) then ! skip 2 generic crops
          ! lines here for consistency; the transfer terms are zero
          livestemn(p)       = livestemn(p)      + livestemn_xfer_to_livestemn(p)*dt
          livestemn_xfer(p)  = livestemn_xfer(p) - livestemn_xfer_to_livestemn(p)*dt
          grainn(p)          = grainn(p)         + grainn_xfer_to_grainn(p)*dt
          grainn_xfer(p)     = grainn_xfer(p)    - grainn_xfer_to_grainn(p)*dt
      end if

      ! phenology: litterfall and retranslocation fluxes
      leafn(p)    = leafn(p)    - leafn_to_litter(p)*dt
      frootn(p)   = frootn(p)   - frootn_to_litter(p)*dt
      leafn(p)    = leafn(p)    - leafn_to_retransn(p)*dt
      retransn(p) = retransn(p) + leafn_to_retransn(p)*dt

      ! live wood turnover and retranslocation fluxes
      if (woody(ivt(p)) == 1._r8) then
          livestemn(p)  = livestemn(p)  - livestemn_to_deadstemn(p)*dt
          deadstemn(p)  = deadstemn(p)  + livestemn_to_deadstemn(p)*dt
          livestemn(p)  = livestemn(p)  - livestemn_to_retransn(p)*dt
          retransn(p)   = retransn(p)   + livestemn_to_retransn(p)*dt
          livecrootn(p) = livecrootn(p) - livecrootn_to_deadcrootn(p)*dt
          deadcrootn(p) = deadcrootn(p) + livecrootn_to_deadcrootn(p)*dt
          livecrootn(p) = livecrootn(p) - livecrootn_to_retransn(p)*dt
          retransn(p)   = retransn(p)   + livecrootn_to_retransn(p)*dt
      end if
      if (ivt(p) >= npcropmin) then ! Beth adds retrans from froot
          frootn(p)     = frootn(p)     - frootn_to_retransn(p)*dt
          retransn(p)   = retransn(p)   + frootn_to_retransn(p)*dt
          livestemn(p)  = livestemn(p)  - livestemn_to_litter(p)*dt
          livestemn(p)  = livestemn(p)  - livestemn_to_retransn(p)*dt
          retransn(p)   = retransn(p)   + livestemn_to_retransn(p)*dt
          grainn(p)     = grainn(p)     - grainn_to_food(p)*dt
      end if

      ! uptake from soil mineral N pool
      npool(p) = npool(p) + sminn_to_npool(p)*dt

      ! deployment from retranslocation pool
      npool(p)    = npool(p)    + retransn_to_npool(p)*dt
      retransn(p) = retransn(p) - retransn_to_npool(p)*dt

      ! allocation fluxes
      npool(p)           = npool(p)          - npool_to_leafn(p)*dt
      leafn(p)           = leafn(p)          + npool_to_leafn(p)*dt
      npool(p)           = npool(p)          - npool_to_leafn_storage(p)*dt
      leafn_storage(p)   = leafn_storage(p)  + npool_to_leafn_storage(p)*dt
      npool(p)           = npool(p)          - npool_to_frootn(p)*dt
      frootn(p)          = frootn(p)         + npool_to_frootn(p)*dt
      npool(p)           = npool(p)          - npool_to_frootn_storage(p)*dt
      frootn_storage(p)  = frootn_storage(p) + npool_to_frootn_storage(p)*dt
      if (woody(ivt(p)) == 1._r8) then
          npool(p)              = npool(p)              - npool_to_livestemn(p)*dt
          livestemn(p)          = livestemn(p)          + npool_to_livestemn(p)*dt
          npool(p)              = npool(p)              - npool_to_livestemn_storage(p)*dt
          livestemn_storage(p)  = livestemn_storage(p)  + npool_to_livestemn_storage(p)*dt
          npool(p)              = npool(p)              - npool_to_deadstemn(p)*dt
          deadstemn(p)          = deadstemn(p)          + npool_to_deadstemn(p)*dt
          npool(p)              = npool(p)              - npool_to_deadstemn_storage(p)*dt
          deadstemn_storage(p)  = deadstemn_storage(p)  + npool_to_deadstemn_storage(p)*dt
          npool(p)              = npool(p)              - npool_to_livecrootn(p)*dt
          livecrootn(p)         = livecrootn(p)         + npool_to_livecrootn(p)*dt
          npool(p)              = npool(p)              - npool_to_livecrootn_storage(p)*dt
          livecrootn_storage(p) = livecrootn_storage(p) + npool_to_livecrootn_storage(p)*dt
          npool(p)              = npool(p)              - npool_to_deadcrootn(p)*dt
          deadcrootn(p)         = deadcrootn(p)         + npool_to_deadcrootn(p)*dt
          npool(p)              = npool(p)              - npool_to_deadcrootn_storage(p)*dt
          deadcrootn_storage(p) = deadcrootn_storage(p) + npool_to_deadcrootn_storage(p)*dt
      end if
      if (ivt(p) >= npcropmin) then ! skip 2 generic crops
          npool(p)              = npool(p)              - npool_to_livestemn(p)*dt
          livestemn(p)          = livestemn(p)          + npool_to_livestemn(p)*dt
          npool(p)              = npool(p)              - npool_to_livestemn_storage(p)*dt
          livestemn_storage(p)  = livestemn_storage(p)  + npool_to_livestemn_storage(p)*dt
          npool(p)              = npool(p)              - npool_to_grainn(p)*dt
          grainn(p)             = grainn(p)             + npool_to_grainn(p)*dt
          npool(p)              = npool(p)              - npool_to_grainn_storage(p)*dt
          grainn_storage(p)     = grainn_storage(p)     + npool_to_grainn_storage(p)*dt
      end if

      ! move storage pools into transfer pools
      leafn_storage(p)  = leafn_storage(p)  - leafn_storage_to_xfer(p)*dt
      leafn_xfer(p)     = leafn_xfer(p)     + leafn_storage_to_xfer(p)*dt
      frootn_storage(p) = frootn_storage(p) - frootn_storage_to_xfer(p)*dt
      frootn_xfer(p)    = frootn_xfer(p)    + frootn_storage_to_xfer(p)*dt
      if (woody(ivt(p)) == 1._r8) then
          livestemn_storage(p)  = livestemn_storage(p)  - livestemn_storage_to_xfer(p)*dt
          livestemn_xfer(p)     = livestemn_xfer(p)     + livestemn_storage_to_xfer(p)*dt
          deadstemn_storage(p)  = deadstemn_storage(p)  - deadstemn_storage_to_xfer(p)*dt
          deadstemn_xfer(p)     = deadstemn_xfer(p)     + deadstemn_storage_to_xfer(p)*dt
          livecrootn_storage(p) = livecrootn_storage(p) - livecrootn_storage_to_xfer(p)*dt
          livecrootn_xfer(p)    = livecrootn_xfer(p)    + livecrootn_storage_to_xfer(p)*dt
          deadcrootn_storage(p) = deadcrootn_storage(p) - deadcrootn_storage_to_xfer(p)*dt
          deadcrootn_xfer(p)    = deadcrootn_xfer(p)    + deadcrootn_storage_to_xfer(p)*dt
      end if
      if (ivt(p) >= npcropmin) then ! skip 2 generic crops
          ! lines here for consistency; the transfer terms are zero
          livestemn_storage(p)  = livestemn_storage(p) - livestemn_storage_to_xfer(p)*dt
          livestemn_xfer(p)     = livestemn_xfer(p)    + livestemn_storage_to_xfer(p)*dt
          grainn_storage(p)     = grainn_storage(p)    - grainn_storage_to_xfer(p)*dt
          grainn_xfer(p)        = grainn_xfer(p)       + grainn_storage_to_xfer(p)*dt
      end if

   end do

 end associate
end subroutine NStateUpdate1

end module CNNStateUpdate1Mod
