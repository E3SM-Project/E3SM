module CNPrecisionControlMod
  !----------------------------------------------------------------------- 
  ! !MODULE: CNPrecisionControlMod
  ! 
  ! !DESCRIPTION:
  ! controls on very low values in critical state variables 
  ! 
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varpar  , only: ndecomp_pools
  implicit none
  save
  private
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: CNPrecisionControl
  !
  ! !REVISION HISTORY:
  ! 4/23/2004: Created by Peter Thornton
  !----------------------------------------------------------------------- 

contains

  !-----------------------------------------------------------------------
  subroutine CNPrecisionControl(num_soilc, filter_soilc, num_soilp, filter_soilp)
    !
    ! !DESCRIPTION: 
    ! On the radiation time step, force leaf and deadstem c and n to 0 if
    ! they get too small.
    !
    ! !USES:
    use clmtype
    use clm_varctl,   only: iulog, use_c13, use_c14, use_nitrif_denitrif
    use clm_varpar,   only: nlevdecomp, crop_prog
    use pftvarcon,    only: nc3crop
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_soilc       ! number of soil columns in filter
    integer, intent(in) :: filter_soilc(:) ! filter for soil columns
    integer, intent(in) :: num_soilp       ! number of soil pfts in filter
    integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
    !
    ! !REVISION HISTORY:
    ! 8/1/03: Created by Peter Thornton
    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,k  ! indices
    integer :: fp,fc    ! lake filter indices
    real(r8):: pc,pn    ! truncation terms for pft-level corrections
    real(r8):: cc,cn    ! truncation terms for column-level corrections
    real(r8):: pc13     ! truncation terms for pft-level corrections
    real(r8):: cc13     ! truncation terms for column-level corrections
    real(r8):: pc14     ! truncation terms for pft-level corrections
    real(r8):: cc14     ! truncation terms for column-level corrections
    real(r8):: ccrit    ! critical carbon state value for truncation
    real(r8):: ncrit    ! critical nitrogen state value for truncation
    !-----------------------------------------------------------------------

   associate(& 
   col_ctrunc_vr               =>  ccs%col_ctrunc_vr                   , & ! Input:  [real(r8) (:,:)]  (gC/m3) column-level sink for C truncation      
   decomp_cpools_vr            =>  ccs%decomp_cpools_vr                , & ! Input:  [real(r8) (:,:,:)]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
   col_ntrunc_vr               =>  cns%col_ntrunc_vr                   , & ! Input:  [real(r8) (:,:)]  (gN/m3) column-level sink for N truncation      
   decomp_npools_vr            =>  cns%decomp_npools_vr                , & ! Input:  [real(r8) (:,:,:)]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
   ivt                         =>  pft%itype                           , & ! Input:  [integer (:)]  pft vegetation type                                
   cpool                       =>  pcs%cpool                           , & ! Input:  [real(r8) (:)]  (gC/m2) temporary photosynthate C pool            
   deadcrootc                  =>  pcs%deadcrootc                      , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C                        
   deadcrootc_storage          =>  pcs%deadcrootc_storage              , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C storage                
   deadcrootc_xfer             =>  pcs%deadcrootc_xfer                 , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C transfer               
   deadstemc                   =>  pcs%deadstemc                       , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C                               
   deadstemc_storage           =>  pcs%deadstemc_storage               , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C storage                       
   deadstemc_xfer              =>  pcs%deadstemc_xfer                  , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C transfer                      
   frootc                      =>  pcs%frootc                          , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C                               
   frootc_storage              =>  pcs%frootc_storage                  , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C storage                       
   frootc_xfer                 =>  pcs%frootc_xfer                     , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C transfer                      
   gresp_storage               =>  pcs%gresp_storage                   , & ! Input:  [real(r8) (:)]  (gC/m2) growth respiration storage                
   gresp_xfer                  =>  pcs%gresp_xfer                      , & ! Input:  [real(r8) (:)]  (gC/m2) growth respiration transfer               
   leafc                       =>  pcs%leafc                           , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C                                    
   leafc_storage               =>  pcs%leafc_storage                   , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C storage                            
   leafc_xfer                  =>  pcs%leafc_xfer                      , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C transfer                           
   livecrootc                  =>  pcs%livecrootc                      , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C                        
   livecrootc_storage          =>  pcs%livecrootc_storage              , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C storage                
   livecrootc_xfer             =>  pcs%livecrootc_xfer                 , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C transfer               
   livestemc                   =>  pcs%livestemc                       , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C                               
   livestemc_storage           =>  pcs%livestemc_storage               , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C storage                       
   livestemc_xfer              =>  pcs%livestemc_xfer                  , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C transfer                      
   pft_ctrunc                  =>  pcs%pft_ctrunc                      , & ! Input:  [real(r8) (:)]  (gC/m2) pft-level sink for C truncation           
   xsmrpool                    =>  pcs%xsmrpool                        , & ! Input:  [real(r8) (:)]  (gC/m2) execss maint resp C pool                  
   grainc                      =>  pcs%grainc                          , & ! Input:  [real(r8) (:)]  (gC/m2) grain C                                   
   grainc_storage              =>  pcs%grainc_storage                  , & ! Input:  [real(r8) (:)]  (gC/m2) grain C storage                           
   grainc_xfer                 =>  pcs%grainc_xfer                     , & ! Input:  [real(r8) (:)]  (gC/m2) grain C transfer                          
   c13_col_ctrunc_vr           =>  cc13s%col_ctrunc_vr                 , & ! Input:  [real(r8) (:,:)]  (gC/m3) column-level sink for C truncation      
   decomp_c13pools_vr          =>  cc13s%decomp_cpools_vr              , & ! Input:  [real(r8) (:,:,:)]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
   c13_cpool                   =>  pc13s%cpool                         , & ! Input:  [real(r8) (:)]  (gC/m2) temporary photosynthate C pool            
   c13_deadcrootc              =>  pc13s%deadcrootc                    , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C                        
   c13_deadcrootc_storage      =>  pc13s%deadcrootc_storage            , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C storage                
   c13_deadcrootc_xfer         =>  pc13s%deadcrootc_xfer               , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C transfer               
   c13_deadstemc               =>  pc13s%deadstemc                     , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C                               
   c13_deadstemc_storage       =>  pc13s%deadstemc_storage             , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C storage                       
   c13_deadstemc_xfer          =>  pc13s%deadstemc_xfer                , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C transfer                      
   c13_frootc                  =>  pc13s%frootc                        , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C                               
   c13_frootc_storage          =>  pc13s%frootc_storage                , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C storage                       
   c13_frootc_xfer             =>  pc13s%frootc_xfer                   , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C transfer                      
   c13_gresp_storage           =>  pc13s%gresp_storage                 , & ! Input:  [real(r8) (:)]  (gC/m2) growth respiration storage                
   c13_gresp_xfer              =>  pc13s%gresp_xfer                    , & ! Input:  [real(r8) (:)]  (gC/m2) growth respiration transfer               
   c13_leafc                   =>  pc13s%leafc                         , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C                                    
   c13_leafc_storage           =>  pc13s%leafc_storage                 , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C storage                            
   c13_leafc_xfer              =>  pc13s%leafc_xfer                    , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C transfer                           
   c13_livecrootc              =>  pc13s%livecrootc                    , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C                        
   c13_livecrootc_storage      =>  pc13s%livecrootc_storage            , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C storage                
   c13_livecrootc_xfer         =>  pc13s%livecrootc_xfer               , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C transfer               
   c13_livestemc               =>  pc13s%livestemc                     , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C                               
   c13_livestemc_storage       =>  pc13s%livestemc_storage             , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C storage                       
   c13_livestemc_xfer          =>  pc13s%livestemc_xfer                , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C transfer                      
   c13_pft_ctrunc              =>  pc13s%pft_ctrunc                    , & ! Input:  [real(r8) (:)]  (gC/m2) pft-level sink for C truncation           
   c14_col_ctrunc_vr           =>  cc14s%col_ctrunc_vr                 , & ! Input:  [real(r8) (:,:)]  (gC/m3) column-level sink for C truncation      
   decomp_c14pools_vr          =>  cc14s%decomp_cpools_vr              , & ! Input:  [real(r8) (:,:,:)]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
   c14_cpool                   =>  pc14s%cpool                         , & ! Input:  [real(r8) (:)]  (gC/m2) temporary photosynthate C pool            
   c14_deadcrootc              =>  pc14s%deadcrootc                    , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C                        
   c14_deadcrootc_storage      =>  pc14s%deadcrootc_storage            , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C storage                
   c14_deadcrootc_xfer         =>  pc14s%deadcrootc_xfer               , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C transfer               
   c14_deadstemc               =>  pc14s%deadstemc                     , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C                               
   c14_deadstemc_storage       =>  pc14s%deadstemc_storage             , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C storage                       
   c14_deadstemc_xfer          =>  pc14s%deadstemc_xfer                , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C transfer                      
   c14_frootc                  =>  pc14s%frootc                        , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C                               
   c14_frootc_storage          =>  pc14s%frootc_storage                , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C storage                       
   c14_frootc_xfer             =>  pc14s%frootc_xfer                   , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C transfer                      
   c14_gresp_storage           =>  pc14s%gresp_storage                 , & ! Input:  [real(r8) (:)]  (gC/m2) growth respiration storage                
   c14_gresp_xfer              =>  pc14s%gresp_xfer                    , & ! Input:  [real(r8) (:)]  (gC/m2) growth respiration transfer               
   c14_leafc                   =>  pc14s%leafc                         , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C                                    
   c14_leafc_storage           =>  pc14s%leafc_storage                 , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C storage                            
   c14_leafc_xfer              =>  pc14s%leafc_xfer                    , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C transfer                           
   c14_livecrootc              =>  pc14s%livecrootc                    , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C                        
   c14_livecrootc_storage      =>  pc14s%livecrootc_storage            , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C storage                
   c14_livecrootc_xfer         =>  pc14s%livecrootc_xfer               , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C transfer               
   c14_livestemc               =>  pc14s%livestemc                     , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C                               
   c14_livestemc_storage       =>  pc14s%livestemc_storage             , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C storage                       
   c14_livestemc_xfer          =>  pc14s%livestemc_xfer                , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C transfer                      
   c14_pft_ctrunc              =>  pc14s%pft_ctrunc                    , & ! Input:  [real(r8) (:)]  (gC/m2) pft-level sink for C truncation           
   deadcrootn                  =>  pns%deadcrootn                      , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N                        
   deadcrootn_storage          =>  pns%deadcrootn_storage              , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N storage                
   deadcrootn_xfer             =>  pns%deadcrootn_xfer                 , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N transfer               
   deadstemn                   =>  pns%deadstemn                       , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N                               
   deadstemn_storage           =>  pns%deadstemn_storage               , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N storage                       
   deadstemn_xfer              =>  pns%deadstemn_xfer                  , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N transfer                      
   frootn                      =>  pns%frootn                          , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N                               
   frootn_storage              =>  pns%frootn_storage                  , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N storage                       
   frootn_xfer                 =>  pns%frootn_xfer                     , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N transfer                      
   leafn                       =>  pns%leafn                           , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N                                    
   leafn_storage               =>  pns%leafn_storage                   , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N storage                            
   leafn_xfer                  =>  pns%leafn_xfer                      , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N transfer                           
   livecrootn                  =>  pns%livecrootn                      , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N                        
   livecrootn_storage          =>  pns%livecrootn_storage              , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N storage                
   livecrootn_xfer             =>  pns%livecrootn_xfer                 , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N transfer               
   grainn                      =>  pns%grainn                          , & ! Input:  [real(r8) (:)]  (gC/m2) grain N                                   
   grainn_storage              =>  pns%grainn_storage                  , & ! Input:  [real(r8) (:)]  (gC/m2) grain N storage                           
   grainn_xfer                 =>  pns%grainn_xfer                     , & ! Input:  [real(r8) (:)]  (gC/m2) grain N transfer                          
   livestemn                   =>  pns%livestemn                       , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N                               
   livestemn_storage           =>  pns%livestemn_storage               , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N storage                       
   livestemn_xfer              =>  pns%livestemn_xfer                  , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N transfer                      
   npool                       =>  pns%npool                           , & ! Input:  [real(r8) (:)]  (gN/m2) temporary plant N pool                    
   pft_ntrunc                  =>  pns%pft_ntrunc                      , & ! Input:  [real(r8) (:)]  (gN/m2) pft-level sink for N truncation           
   smin_nh4_vr                 =>  cns%smin_nh4_vr                     , & ! Input:  [real(r8) (:,:)]  (gN/m3) soil mineral NH4                        
   smin_no3_vr                 =>  cns%smin_no3_vr                     , & ! Input:  [real(r8) (:,:)]  (gN/m3) soil mineral NO3                        
   retransn                    =>  pns%retransn                          & ! Input:  [real(r8) (:)]  (gN/m2) plant pool of retranslocated N            
   )

   ! set the critical carbon state value for truncation (gC/m2)
   ccrit = 1.e-8_r8
   ! set the critical nitrogen state value for truncation (gN/m2)
   ncrit = 1.e-8_r8
    
   ! pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)
      
      ! initialize the pft-level C and N truncation terms
      pc = 0._r8

      if ( use_c13 ) then
         pc13 = 0._r8
      end if
      if ( use_c14 ) then
         pc14 = 0._r8
      endif

      pn = 0._r8
      
      ! do tests on state variables for precision control
      ! for linked C-N state variables, perform precision test on
      ! the C component, but truncate C, C13, and N components
      
      ! leaf C and N
      if (abs(leafc(p)) < ccrit) then
          pc = pc + leafc(p)
          leafc(p) = 0._r8

          if ( use_c13 ) then
             pc13 = pc13 + c13_leafc(p)
             c13_leafc(p) = 0._r8
          endif

          if ( use_c14 ) then
             pc14 = pc14 + c14_leafc(p)
             c14_leafc(p) = 0._r8
          endif

          pn = pn + leafn(p)
          leafn(p) = 0._r8
      end if

      ! leaf storage C and N
      if (abs(leafc_storage(p)) < ccrit) then
          pc = pc + leafc_storage(p)
          leafc_storage(p) = 0._r8

          if ( use_c13 ) then
             pc13 = pc13 + c13_leafc_storage(p)
             c13_leafc_storage(p) = 0._r8
          endif

          if ( use_c14 ) then
             pc14 = pc14 + c14_leafc_storage(p)
             c14_leafc_storage(p) = 0._r8
          endif

          pn = pn + leafn_storage(p)
          leafn_storage(p) = 0._r8
      end if
          
      ! leaf transfer C and N
      if (abs(leafc_xfer(p)) < ccrit) then
          pc = pc + leafc_xfer(p)
          leafc_xfer(p) = 0._r8

          if ( use_c13 ) then
             pc13 = pc13 + c13_leafc_xfer(p)
             c13_leafc_xfer(p) = 0._r8
          endif
          if ( use_c14 ) then
             pc14 = pc14 + c14_leafc_xfer(p)
             c14_leafc_xfer(p) = 0._r8
          endif

          pn = pn + leafn_xfer(p)
          leafn_xfer(p) = 0._r8
      end if
          
      ! froot C and N
      if (abs(frootc(p)) < ccrit) then
          pc = pc + frootc(p)
          frootc(p) = 0._r8

          if ( use_c13 ) then
             pc13 = pc13 + c13_frootc(p)
             c13_frootc(p) = 0._r8
          endif
          if ( use_c14 ) then
             pc14 = pc14 + c14_frootc(p)
             c14_frootc(p) = 0._r8
          endif

          pn = pn + frootn(p)
          frootn(p) = 0._r8
      end if

      ! froot storage C and N
      if (abs(frootc_storage(p)) < ccrit) then
          pc = pc + frootc_storage(p)
          frootc_storage(p) = 0._r8

          if ( use_c13 ) then
             pc13 = pc13 + c13_frootc_storage(p)
             c13_frootc_storage(p) = 0._r8
          endif
          if ( use_c14 ) then
             pc14 = pc14 + c14_frootc_storage(p)
             c14_frootc_storage(p) = 0._r8
          endif

          pn = pn + frootn_storage(p)
          frootn_storage(p) = 0._r8
      end if
          
      ! froot transfer C and N
      if (abs(frootc_xfer(p)) < ccrit) then
          pc = pc + frootc_xfer(p)
          frootc_xfer(p) = 0._r8

          if ( use_c13 ) then
             pc13 = pc13 + c13_frootc_xfer(p)
             c13_frootc_xfer(p) = 0._r8
          endif
          if ( use_c14 ) then
             pc14 = pc14 + c14_frootc_xfer(p)
             c14_frootc_xfer(p) = 0._r8
          endif

          pn = pn + frootn_xfer(p)
          frootn_xfer(p) = 0._r8
      end if
          
      if ( crop_prog .and. ivt(p) >= nc3crop )then
         ! grain C and N
         if (abs(grainc(p)) < ccrit) then
             pc = pc + grainc(p)
             grainc(p) = 0._r8
             pn = pn + grainn(p)
             grainn(p) = 0._r8
         end if
   
         ! grain storage C and N
         if (abs(grainc_storage(p)) < ccrit) then
             pc = pc + grainc_storage(p)
             grainc_storage(p) = 0._r8
             pn = pn + grainn_storage(p)
             grainn_storage(p) = 0._r8
         end if
             
         ! grain transfer C and N
         if (abs(grainc_xfer(p)) < ccrit) then
             pc = pc + grainc_xfer(p)
             grainc_xfer(p) = 0._r8
             pn = pn + grainn_xfer(p)
             grainn_xfer(p) = 0._r8
         end if
      end if
          
      ! livestem C and N
      if (abs(livestemc(p)) < ccrit) then
          pc = pc + livestemc(p)
          livestemc(p) = 0._r8

          if ( use_c13 ) then
             pc13 = pc13 + c13_livestemc(p)
             c13_livestemc(p) = 0._r8
          endif
          if ( use_c14 ) then
             pc14 = pc14 + c14_livestemc(p)
             c14_livestemc(p) = 0._r8
          endif

          pn = pn + livestemn(p)
          livestemn(p) = 0._r8
      end if

      ! livestem storage C and N
      if (abs(livestemc_storage(p)) < ccrit) then
          pc = pc + livestemc_storage(p)
          livestemc_storage(p) = 0._r8
          if ( use_c13 ) then
             pc13 = pc13 + c13_livestemc_storage(p)
             c13_livestemc_storage(p) = 0._r8
          endif
          if ( use_c14 ) then
             pc14 = pc14 + c14_livestemc_storage(p)
             c14_livestemc_storage(p) = 0._r8
          endif
          pn = pn + livestemn_storage(p)
          livestemn_storage(p) = 0._r8
      end if
          
      ! livestem transfer C and N
      if (abs(livestemc_xfer(p)) < ccrit) then
          pc = pc + livestemc_xfer(p)
          livestemc_xfer(p) = 0._r8
          if ( use_c13 ) then
             pc13 = pc13 + c13_livestemc_xfer(p)
             c13_livestemc_xfer(p) = 0._r8
          endif
          if ( use_c14 ) then
             pc14 = pc14 + c14_livestemc_xfer(p)
             c14_livestemc_xfer(p) = 0._r8
          endif
          pn = pn + livestemn_xfer(p)
          livestemn_xfer(p) = 0._r8
      end if
          
      ! deadstem C and N
      if (abs(deadstemc(p)) < ccrit) then
          pc = pc + deadstemc(p)
          deadstemc(p) = 0._r8
          if ( use_c13 ) then
             pc13 = pc13 + c13_deadstemc(p)
             c13_deadstemc(p) = 0._r8
          endif
          if ( use_c14 ) then
             pc14 = pc14 + c14_deadstemc(p)
             c14_deadstemc(p) = 0._r8
          endif
          pn = pn + deadstemn(p)
          deadstemn(p) = 0._r8
      end if

      ! deadstem storage C and N
      if (abs(deadstemc_storage(p)) < ccrit) then
          pc = pc + deadstemc_storage(p)
          deadstemc_storage(p) = 0._r8
          if ( use_c13 ) then
             pc13 = pc13 + c13_deadstemc_storage(p)
             c13_deadstemc_storage(p) = 0._r8
          endif
          if ( use_c14 ) then
             pc14 = pc14 + c14_deadstemc_storage(p)
             c14_deadstemc_storage(p) = 0._r8
          endif
          pn = pn + deadstemn_storage(p)
          deadstemn_storage(p) = 0._r8
      end if
          
      ! deadstem transfer C and N
      if (abs(deadstemc_xfer(p)) < ccrit) then
          pc = pc + deadstemc_xfer(p)
          deadstemc_xfer(p) = 0._r8
          if ( use_c13 ) then
             pc13 = pc13 + c13_deadstemc_xfer(p)
             c13_deadstemc_xfer(p) = 0._r8
          endif
          if ( use_c14 ) then
             pc14 = pc14 + c14_deadstemc_xfer(p)
             c14_deadstemc_xfer(p) = 0._r8
          endif
          pn = pn + deadstemn_xfer(p)
          deadstemn_xfer(p) = 0._r8
      end if
          
      ! livecroot C and N
      if (abs(livecrootc(p)) < ccrit) then
          pc = pc + livecrootc(p)
          livecrootc(p) = 0._r8

          if ( use_c13 ) then
             pc13 = pc13 + c13_livecrootc(p)
             c13_livecrootc(p) = 0._r8
          endif
          if ( use_c14 ) then
             pc14 = pc14 + c14_livecrootc(p)
             c14_livecrootc(p) = 0._r8
          endif
          pn = pn + livecrootn(p)
          livecrootn(p) = 0._r8
      end if

      ! livecroot storage C and N
      if (abs(livecrootc_storage(p)) < ccrit) then
          pc = pc + livecrootc_storage(p)
          livecrootc_storage(p) = 0._r8

          if ( use_c13 ) then
             pc13 = pc13 + c13_livecrootc_storage(p)
             c13_livecrootc_storage(p) = 0._r8
          endif

          if ( use_c14 ) then
             pc14 = pc14 + c14_livecrootc_storage(p)
             c14_livecrootc_storage(p) = 0._r8
          endif

          pn = pn + livecrootn_storage(p)
          livecrootn_storage(p) = 0._r8
      end if
          
      ! livecroot transfer C and N
      if (abs(livecrootc_xfer(p)) < ccrit) then
          pc = pc + livecrootc_xfer(p)
          livecrootc_xfer(p) = 0._r8
          if ( use_c13 ) then
             pc13 = pc13 + c13_livecrootc_xfer(p)
             c13_livecrootc_xfer(p) = 0._r8
          endif
          if ( use_c14 ) then
             pc14 = pc14 + c14_livecrootc_xfer(p)
             c14_livecrootc_xfer(p) = 0._r8
          endif
          pn = pn + livecrootn_xfer(p)
          livecrootn_xfer(p) = 0._r8
      end if
          
      ! deadcroot C and N
      if (abs(deadcrootc(p)) < ccrit) then
          pc = pc + deadcrootc(p)
          deadcrootc(p) = 0._r8
          if ( use_c13 ) then
             pc13 = pc13 + c13_deadcrootc(p)
             c13_deadcrootc(p) = 0._r8
          endif
          if ( use_c14 ) then
             pc14 = pc14 + c14_deadcrootc(p)
             c14_deadcrootc(p) = 0._r8
          endif
          pn = pn + deadcrootn(p)
          deadcrootn(p) = 0._r8
      end if

      ! deadcroot storage C and N
      if (abs(deadcrootc_storage(p)) < ccrit) then
          pc = pc + deadcrootc_storage(p)
          deadcrootc_storage(p) = 0._r8
          if ( use_c13 ) then
             pc13 = pc13 + c13_deadcrootc_storage(p)
             c13_deadcrootc_storage(p) = 0._r8
          endif
          if ( use_c14 ) then
             pc14 = pc14 + c14_deadcrootc_storage(p)
             c14_deadcrootc_storage(p) = 0._r8
          endif
          pn = pn + deadcrootn_storage(p)
          deadcrootn_storage(p) = 0._r8
      end if
          
      ! deadcroot transfer C and N
      if (abs(deadcrootc_xfer(p)) < ccrit) then
          pc = pc + deadcrootc_xfer(p)
          deadcrootc_xfer(p) = 0._r8
          if ( use_c13 ) then
             pc13 = pc13 + c13_deadcrootc_xfer(p)
             c13_deadcrootc_xfer(p) = 0._r8
          endif
          if ( use_c14 ) then
             pc14 = pc14 + c14_deadcrootc_xfer(p)
             c14_deadcrootc_xfer(p) = 0._r8
          endif
          pn = pn + deadcrootn_xfer(p)
          deadcrootn_xfer(p) = 0._r8
      end if
          
      ! gresp_storage (C only)
      if (abs(gresp_storage(p)) < ccrit) then
          pc = pc + gresp_storage(p)
          gresp_storage(p) = 0._r8
          if ( use_c13 ) then
             pc13 = pc13 + c13_gresp_storage(p)
             c13_gresp_storage(p) = 0._r8
          endif
          if ( use_c14 ) then
             pc14 = pc14 + c14_gresp_storage(p)
             c14_gresp_storage(p) = 0._r8
          endif
      end if

      ! gresp_xfer (C only)
      if (abs(gresp_xfer(p)) < ccrit) then
          pc = pc + gresp_xfer(p)
          gresp_xfer(p) = 0._r8
          if ( use_c13 ) then
             pc13 = pc13 + c13_gresp_xfer(p)
             c13_gresp_xfer(p) = 0._r8
          endif
          if ( use_c14 ) then
             pc14 = pc14 + c14_gresp_xfer(p)
             c14_gresp_xfer(p) = 0._r8
          endif
      end if
          
      ! cpool (C only)
      if (abs(cpool(p)) < ccrit) then
          pc = pc + cpool(p)
          cpool(p) = 0._r8
          if ( use_c13 ) then
             pc13 = pc13 + c13_cpool(p)
             c13_cpool(p) = 0._r8
          endif
          if ( use_c14 ) then
             pc14 = pc14 + c14_cpool(p)
             c14_cpool(p) = 0._r8
          endif
      end if
          
      if ( crop_prog .and. ivt(p) >= nc3crop )then
         ! xsmrpool (C only)
         if (abs(xsmrpool(p)) < ccrit) then
             pc = pc + xsmrpool(p)
             xsmrpool(p) = 0._r8
         end if
      end if
          
      ! retransn (N only)
      if (abs(retransn(p)) < ncrit) then
          pn = pn + retransn(p)
          retransn(p) = 0._r8
      end if
          
      ! npool (N only)
      if (abs(npool(p)) < ncrit) then
          pn = pn + npool(p)
          npool(p) = 0._r8
      end if
      
      pft_ctrunc(p) = pft_ctrunc(p) + pc
      if ( use_c13 ) then
         c13_pft_ctrunc(p) = c13_pft_ctrunc(p) + pc13
      endif

      if ( use_c14 ) then
         c14_pft_ctrunc(p) = c14_pft_ctrunc(p) + pc14
      endif

      pft_ntrunc(p) = pft_ntrunc(p) + pn
          
   end do ! end of pft loop

   ! column loop
   do fc = 1,num_soilc
      c = filter_soilc(fc)

      do j = 1,nlevdecomp
         ! initialize the column-level C and N truncation terms
         cc = 0._r8
         if ( use_c13 ) then
            cc13 = 0._r8
         endif
         if ( use_c14 ) then
            cc14 = 0._r8
         endif
         cn = 0._r8
         
         ! do tests on state variables for precision control
         ! for linked C-N state variables, perform precision test on
         ! the C component, but truncate both C and N components
         
         ! all decomposing pools C and N
         do k = 1, ndecomp_pools

            if (abs(decomp_cpools_vr(c,j,k)) < ccrit) then
               cc = cc + decomp_cpools_vr(c,j,k)
               decomp_cpools_vr(c,j,k) = 0._r8
               if ( use_c13 ) then
                  cc13 = cc13 + decomp_c13pools_vr(c,j,k)
                  decomp_c13pools_vr(c,j,k) = 0._r8
               endif
               if ( use_c14 ) then
                  cc14 = cc14 + decomp_c14pools_vr(c,j,k)
                  decomp_c14pools_vr(c,j,k) = 0._r8
               endif
               cn = cn + decomp_npools_vr(c,j,k)
               decomp_npools_vr(c,j,k) = 0._r8
            end if

         end do
         
         ! not doing precision control on soil mineral N, since it will
         ! be getting the N truncation flux anyway.
         
         col_ctrunc_vr(c,j) = col_ctrunc_vr(c,j) + cc
         if ( use_c13 ) then
            c13_col_ctrunc_vr(c,j) = c13_col_ctrunc_vr(c,j) + cc13
         endif
         if ( use_c14 ) then
            c14_col_ctrunc_vr(c,j) = c14_col_ctrunc_vr(c,j) + cc14
         endif
         col_ntrunc_vr(c,j) = col_ntrunc_vr(c,j) + cn
      end do

   end do   ! end of column loop

   if (use_nitrif_denitrif) then
      ! remove small negative perturbations for stability purposes, if any should arise.
      
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         do j = 1,nlevdecomp
            if (abs(smin_no3_vr(c,j)) < ncrit/1e4_r8) then
               if ( smin_no3_vr(c,j) .lt. 0._r8 ) then
                  write(iulog, *) '-10^-12 < smin_no3 < 0. resetting to zero.'
                  write(iulog, *) 'smin_no3_vr(c,j), c, j: ', smin_no3_vr(c,j), c, j
                  smin_no3_vr(c,j) = 0._r8
               endif
            end if
            if (abs(smin_nh4_vr(c,j)) < ncrit/1e4_r8) then
               if ( smin_nh4_vr(c,j) .lt. 0._r8 ) then
                  write(iulog, *) '-10^-12 < smin_nh4 < 0. resetting to zero.'
                  write(iulog, *) 'smin_nh4_vr(c,j), c, j: ', smin_nh4_vr(c,j), c, j
                  smin_nh4_vr(c,j) = 0._r8
               endif
            end if
         end do
      end do
   endif
   
    end associate 
 end subroutine CNPrecisionControl

end module CNPrecisionControlMod
