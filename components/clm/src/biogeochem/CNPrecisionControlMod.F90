module CNPrecisionControlMod

  !----------------------------------------------------------------------- 
  ! !DESCRIPTION:
  ! controls on very low values in critical state variables 
  ! 
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use clm_varpar          , only : ndecomp_pools
  use CNCarbonStateType   , only : carbonstate_type
  use CNNitrogenStateType , only : nitrogenstate_type
  use PhosphorusStateType , only : phosphorusstate_type
  use VegetationType           , only : veg_pp
  use ColumnType          , only : col_pp
  use clm_varctl          , only : nu_com
  use abortutils          , only : endrun
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: CNPrecisionControl
  !----------------------------------------------------------------------- 

contains

  !-----------------------------------------------------------------------
  subroutine CNPrecisionControl(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, nitrogenstate_vars,&
       phosphorusstate_vars)
    !
    ! !DESCRIPTION: 
    ! On the radiation time step, force leaf and deadstem c and n to 0 if
    ! they get too small.
    !
    ! !USES:
    use clm_varctl , only : iulog, use_c13, use_c14, use_nitrif_denitrif, use_fates
    use clm_varpar , only : nlevdecomp_full, crop_prog
    use pftvarcon  , only : nc3crop
    use tracer_varcon          , only : is_active_betr_bgc    
    use CNDecompCascadeConType , only : decomp_cascade_con
    !
    ! !ARGUMENTS:
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                  , intent(in)    :: num_soilp       ! number of soil patchs in filter
    integer                  , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(carbonstate_type)   , intent(inout) :: carbonstate_vars
    type(carbonstate_type)   , intent(inout) :: c13_carbonstate_vars
    type(carbonstate_type)   , intent(inout) :: c14_carbonstate_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars
    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,k,l  ! indices
    integer :: fp,fc    ! lake filter indices
    real(r8):: pc,pn,pp    ! truncation terms for patch-level corrections
    real(r8):: cc,cn,cp    ! truncation terms for column-level corrections
    real(r8):: pc13     ! truncation terms for patch-level corrections
    real(r8):: cc13     ! truncation terms for column-level corrections
    real(r8):: pc14     ! truncation terms for patch-level corrections
    real(r8):: cc14     ! truncation terms for column-level corrections
    real(r8):: ccrit    ! critical carbon state value for truncation
    real(r8):: ncrit    ! critical nitrogen state value for truncation
    real(r8):: pcrit    ! critical phosphorus state value for truncation
    real(r8):: cc_eca
    real(r8):: cn_eca
    real(r8):: cp_eca
    !-----------------------------------------------------------------------

    ! carbonstate_vars%soil_trunc_vr_col                 Output:  [real(r8) (:,:)   ]  (gC/m3) column-level sink for C truncation      
    ! carbonstate_vars%decomp_pools_vr_col          Output:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
    ! carbonstate_vars%pool_patch                   Output:  [real(r8) (:)     ]  (gC/m2) temporary photosynthate C pool            
    ! carbonstate_vars%deadcroot_patch              Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C                        
    ! carbonstate_vars%deadcroot_storage_patch      Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C storage                
    ! carbonstate_vars%deadcroot_xfer_patch         Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C transfer               
    ! carbonstate_vars%deadstem_patch               Output:  [real(r8) (:)     ]  (gC/m2) dead stem C                               
    ! carbonstate_vars%deadstem_storage_patch       Output:  [real(r8) (:)     ]  (gC/m2) dead stem C storage                       
    ! carbonstate_vars%deadstem_xfer_patch          Output:  [real(r8) (:)     ]  (gC/m2) dead stem C transfer                      
    ! carbonstate_vars%froot_patch                  Output:  [real(r8) (:)     ]  (gC/m2) fine root C                               
    ! carbonstate_vars%froot_storage_patch          Output:  [real(r8) (:)     ]  (gC/m2) fine root C storage                       
    ! carbonstate_vars%froot_xfer_patch             Output:  [real(r8) (:)     ]  (gC/m2) fine root C transfer                      
    ! carbonstate_vars%gresp_storage_patch           Output:  [real(r8) (:)     ]  (gC/m2) growth respiration storage                
    ! carbonstate_vars%gresp_xfer_patch              Output:  [real(r8) (:)     ]  (gC/m2) growth respiration transfer               
    ! carbonstate_vars%leaf_patch                   Output:  [real(r8) (:)     ]  (gC/m2) leaf C                                    
    ! carbonstate_vars%leaf_storage_patch           Output:  [real(r8) (:)     ]  (gC/m2) leaf C storage                            
    ! carbonstate_vars%leaf_xfer_patch              Output:  [real(r8) (:)     ]  (gC/m2) leaf C transfer                           
    ! carbonstate_vars%livecroot_patch              Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C                        
    ! carbonstate_vars%livecroot_storage_patch      Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C storage                
    ! carbonstate_vars%livecroot_xfer_patch         Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C transfer               
    ! carbonstate_vars%livestem_patch               Output:  [real(r8) (:)     ]  (gC/m2) live stem C                               
    ! carbonstate_vars%livestem_storage_patch       Output:  [real(r8) (:)     ]  (gC/m2) live stem C storage                       
    ! carbonstate_vars%livestem_xfer_patch          Output:  [real(r8) (:)     ]  (gC/m2) live stem C transfer                      
    ! carbonstate_vars%veg_trunc_patch                  Output:  [real(r8) (:)     ]  (gC/m2) patch-level sink for C truncation           
    ! carbonstate_vars%xsmrpool_patch                Output:  [real(r8) (:)     ]  (gC/m2) execss maint resp C pool                  
    ! carbonstate_vars%grain_patch                  Output:  [real(r8) (:)     ]  (gC/m2) grain C                                   
    ! carbonstate_vars%grain_storage_patch          Output:  [real(r8) (:)     ]  (gC/m2) grain C storage                           
    ! carbonstate_vars%grain_xfer_patch             Output:  [real(r8) (:)     ]  (gC/m2) grain C transfer                          
    
    ! c13_carbonstate_vars%soil_trunc_vr_col             Output:  [real(r8) (:,:)   ]  (gC/m3) column-level sink for C truncation      
    ! c13_carbonstate_vars%decomp_pools_vr_col      Output:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
    ! c13_carbonstate_vars%pool_patch               Output:  [real(r8) (:)     ]  (gC/m2) temporary photosynthate C pool            
    ! c13_carbonstate_vars%deadcroot_patch          Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C                        
    ! c13_carbonstate_vars%deadcroot_storage_patch  Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C storage                
    ! c13_carbonstate_vars%deadcroot_xfer_patch     Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C transfer               
    ! c13_carbonstate_vars%deadstem_patch           Output:  [real(r8) (:)     ]  (gC/m2) dead stem C                               
    ! c13_carbonstate_vars%deadstem_storage_patch   Output:  [real(r8) (:)     ]  (gC/m2) dead stem C storage                       
    ! c13_carbonstate_vars%deadstem_xfer_patch      Output:  [real(r8) (:)     ]  (gC/m2) dead stem C transfer                      
    ! c13_carbonstate_vars%froot_patch              Output:  [real(r8) (:)     ]  (gC/m2) fine root C                               
    ! c13_carbonstate_vars%froot_storage_patch      Output:  [real(r8) (:)     ]  (gC/m2) fine root C storage                       
    ! c13_carbonstate_vars%froot_xfer_patch         Output:  [real(r8) (:)     ]  (gC/m2) fine root C transfer                      
    ! c13_carbonstate_vars%gresp_storage_patch       Output:  [real(r8) (:)     ]  (gC/m2) growth respiration storage                
    ! c13_carbonstate_vars%gresp_xfer_patch          Output:  [real(r8) (:)     ]  (gC/m2) growth respiration transfer               
    ! c13_carbonstate_vars%leaf_patch               Output:  [real(r8) (:)     ]  (gC/m2) leaf C                                    
    ! c13_carbonstate_vars%leaf_storage_patch       Output:  [real(r8) (:)     ]  (gC/m2) leaf C storage                            
    ! c13_carbonstate_vars%leaf_xfer_patch          Output:  [real(r8) (:)     ]  (gC/m2) leaf C transfer                           
    ! c13_carbonstate_vars%livecroot_patch          Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C                        
    ! c13_carbonstate_vars%livecroot_storage_patch  Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C storage                
    ! c13_carbonstate_vars%livecroot_xfer_patch     Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C transfer               
    ! c13_carbonstate_vars%livestem_patch           Output:  [real(r8) (:)     ]  (gC/m2) live stem C                               
    ! c13_carbonstate_vars%livestem_storage_patch   Output:  [real(r8) (:)     ]  (gC/m2) live stem C storage                       
    ! c13_carbonstate_vars%livestem_xfer_patch      Output:  [real(r8) (:)     ]  (gC/m2) live stem C transfer                      
    ! c13_carbonstate_vars%veg_trunc_patch              Output:  [real(r8) (:)     ]  (gC/m2) patch-level sink for C truncation           
    
    ! c14_carbonstate_vars%soil_trunc_vr_col             Output:  [real(r8) (:,:)   ]  (gC/m3) column-level sink for C truncation      
    ! c14_carbonstate_vars%decomp_pools_vr_col      Output:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
    ! c14_carbonstate_vars%pool_patch               Output:  [real(r8) (:)     ]  (gC/m2) temporary photosynthate C pool            
    ! c14_carbonstate_vars%deadcroot_patch          Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C                        
    ! c14_carbonstate_vars%deadcroot_storage_patch  Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C storage                
    ! c14_carbonstate_vars%deadcroot_xfer_patch     Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C transfer               
    ! c14_carbonstate_vars%deadstem_patch           Output:  [real(r8) (:)     ]  (gC/m2) dead stem C                               
    ! c14_carbonstate_vars%deadstem_storage_patch   Output:  [real(r8) (:)     ]  (gC/m2) dead stem C storage                       
    ! c14_carbonstate_vars%deadstem_xfer_patch      Output:  [real(r8) (:)     ]  (gC/m2) dead stem C transfer                      
    ! c14_carbonstate_vars%froot_patch              Output:  [real(r8) (:)     ]  (gC/m2) fine root C                               
    ! c14_carbonstate_vars%froot_storage_patch      Output:  [real(r8) (:)     ]  (gC/m2) fine root C storage                       
    ! c14_carbonstate_vars%froot_xfer_patch         Output:  [real(r8) (:)     ]  (gC/m2) fine root C transfer                      
    ! c14_carbonstate_vars%gresp_storage_patch       Output:  [real(r8) (:)     ]  (gC/m2) growth respiration storage                
    ! c14_carbonstate_vars%gresp_xfer_patch          Output:  [real(r8) (:)     ]  (gC/m2) growth respiration transfer               
    ! c14_carbonstate_vars%leaf_patch               Output:  [real(r8) (:)     ]  (gC/m2) leaf C                                    
    ! c14_carbonstate_vars%leaf_storage_patch       Output:  [real(r8) (:)     ]  (gC/m2) leaf C storage                            
    ! c14_carbonstate_vars%leaf_xfer_patch          Output:  [real(r8) (:)     ]  (gC/m2) leaf C transfer                           
    ! c14_carbonstate_vars%livecroot_patch          Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C                        
    ! c14_carbonstate_vars%livecroot_storage_patch  Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C storage                
    ! c14_carbonstate_vars%livecroot_xfer_patch     Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C transfer               
    ! c14_carbonstate_vars%livestem_patch           Output:  [real(r8) (:)     ]  (gC/m2) live stem C                               
    ! c14_carbonstate_vars%livestem_storage_patch   Output:  [real(r8) (:)     ]  (gC/m2) live stem C storage                       
    ! c14_carbonstate_vars%livestem_xfer_patch      Output:  [real(r8) (:)     ]  (gC/m2) live stem C transfer                      
    ! c14_carbonstate_vars%veg_trunc_patch              Output:  [real(r8) (:)     ]  (gC/m2) patch-level sink for C truncation           
    
    ! nitrogenstate_vars%soil_trunc_vr_col               Output:  [real(r8) (:,:)   ]  (gN/m3) column-level sink for N truncation      
    ! nitrogenstate_vars%decomp_pools_vr_col        Output:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
    ! nitrogenstate_vars%deadcroot_patch            Output:  [real(r8) (:)     ]  (gN/m2) dead coarse root N                        
    ! nitrogenstate_vars%deadcroot_storage_patch    Output:  [real(r8) (:)     ]  (gN/m2) dead coarse root N storage                
    ! nitrogenstate_vars%deadcroot_xfer_patch       Output:  [real(r8) (:)     ]  (gN/m2) dead coarse root N transfer               
    ! nitrogenstate_vars%deadstem_patch             Output:  [real(r8) (:)     ]  (gN/m2) dead stem N                               
    ! nitrogenstate_vars%deadstem_storage_patch     Output:  [real(r8) (:)     ]  (gN/m2) dead stem N storage                       
    ! nitrogenstate_vars%deadstem_xfer_patch        Output:  [real(r8) (:)     ]  (gN/m2) dead stem N transfer                      
    ! nitrogenstate_vars%froot_patch                Output:  [real(r8) (:)     ]  (gN/m2) fine root N                               
    ! nitrogenstate_vars%froot_storage_patch        Output:  [real(r8) (:)     ]  (gN/m2) fine root N storage                       
    ! nitrogenstate_vars%froot_xfer_patch           Output:  [real(r8) (:)     ]  (gN/m2) fine root N transfer                      
    ! nitrogenstate_vars%leaf_patch                 Output:  [real(r8) (:)     ]  (gN/m2) leaf N                                    
    ! nitrogenstate_vars%leaf_storage_patch         Output:  [real(r8) (:)     ]  (gN/m2) leaf N storage                            
    ! nitrogenstate_vars%leaf_xfer_patch            Output:  [real(r8) (:)     ]  (gN/m2) leaf N transfer                           
    ! nitrogenstate_vars%livecroot_patch            Output:  [real(r8) (:)     ]  (gN/m2) live coarse root N                        
    ! nitrogenstate_vars%livecroot_storage_patch    Output:  [real(r8) (:)     ]  (gN/m2) live coarse root N storage                
    ! nitrogenstate_vars%livecroot_xfer_patch       Output:  [real(r8) (:)     ]  (gN/m2) live coarse root N transfer               
    ! nitrogenstate_vars%grain_patch                Output:  [real(r8) (:)     ]  (gC/m2) grain N                                   
    ! nitrogenstate_vars%grain_storage_patch        Output:  [real(r8) (:)     ]  (gC/m2) grain N storage                           
    ! nitrogenstate_vars%grain_xfer_patch           Output:  [real(r8) (:)     ]  (gC/m2) grain N transfer                          
    ! nitrogenstate_vars%livestem_patch             Output:  [real(r8) (:)     ]  (gN/m2) live stem N                               
    ! nitrogenstate_vars%livestem_storage_patch     Output:  [real(r8) (:)     ]  (gN/m2) live stem N storage                       
    ! nitrogenstate_vars%livestem_xfer_patch        Output:  [real(r8) (:)     ]  (gN/m2) live stem N transfer                      
    ! nitrogenstate_vars%pool_patch                 Output:  [real(r8) (:)     ]  (gN/m2) temporary plant N pool                    
    ! nitrogenstate_vars%veg_trunc_patch                Output:  [real(r8) (:)     ]  (gN/m2) patch-level sink for N truncation           
    ! nitrogenstate_vars%retransn_patch              Output:  [real(r8) (:)     ]  (gN/m2) plant pool of retranslocated N            
    ! nitrogenstate_vars%smin_nh4_vr_col             Output:  [real(r8) (:,:)   ]  (gN/m3) soil mineral NH4                        
    ! nitrogenstate_vars%smin_no3_vr_col             Output:  [real(r8) (:,:)   ]  (gN/m3) soil mineral NO3                        
    
    associate(&
         cs    => carbonstate_vars     , &
         ns    => nitrogenstate_vars   , &
         ps    => phosphorusstate_vars , &
         c13cs => c13_carbonstate_vars , &
         c14cs => c14_carbonstate_vars , &
         floating_cn_ratio_decomp_pools   =>    decomp_cascade_con%floating_cn_ratio_decomp_pools , &
         floating_cp_ratio_decomp_pools   =>    decomp_cascade_con%floating_cp_ratio_decomp_pools , &
         initial_cn_ratio                 =>    decomp_cascade_con%initial_cn_ratio                 &
         )

      ! set the critical carbon state value for truncation (gC/m2)
      ccrit = 1.e-8_r8

      ! set the critical nitrogen state value for truncation (gN/m2)
      ncrit = 1.e-8_r8

      ! set the critical phosphorus state value for truncation (gN/m2)
      pcrit = 1.e-8_r8

      ! patch loop
      if (.not.use_fates) then
         do fp = 1,num_soilp
            p = filter_soilp(fp)

            ! initialize the patch-level C and N truncation terms
            pc = 0._r8
            pn = 0._r8
            pp = 0._r8
            if ( use_c13 ) pc13 = 0._r8
            if ( use_c14 ) pc14 = 0._r8

            ! do tests on state variables for precision control
            ! for linked C-N state variables, perform precision test on
            ! the C component, but truncate C, C13, and N components

            ! leaf C and N
            if (abs(cs%leaf_patch(p)) < ccrit) then
               pc = pc + cs%leaf_patch(p)
               cs%leaf_patch(p) = 0._r8
               pn = pn + ns%leaf_patch(p)
               ns%leaf_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%leaf_patch(p)
                  c13cs%leaf_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%leaf_patch(p)
                  c14cs%leaf_patch(p) = 0._r8
               endif

               pp = pp + ps%leaf_patch(p)
               ps%leaf_patch(p) = 0._r8
            end if

            ! leaf storage C and N
            if (abs(cs%leaf_storage_patch(p)) < ccrit) then
               pc = pc + cs%leaf_storage_patch(p)
               cs%leaf_storage_patch(p) = 0._r8
               pn = pn + ns%leaf_storage_patch(p)
               ns%leaf_storage_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%leaf_storage_patch(p)
                  c13cs%leaf_storage_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%leaf_storage_patch(p)
                  c14cs%leaf_storage_patch(p) = 0._r8
               endif

               pp = pp + ps%leaf_storage_patch(p)
               ps%leaf_storage_patch(p) = 0._r8
            end if

            ! leaf transfer C and N
            if (abs(cs%leaf_xfer_patch(p)) < ccrit) then
               pc = pc + cs%leaf_xfer_patch(p)
               cs%leaf_xfer_patch(p) = 0._r8
               pn = pn + ns%leaf_xfer_patch(p)
               ns%leaf_xfer_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%leaf_xfer_patch(p)
                  c13cs%leaf_xfer_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%leaf_xfer_patch(p)
                  c14cs%leaf_xfer_patch(p) = 0._r8
               endif

               pp = pp + ps%leaf_xfer_patch(p)
               ps%leaf_xfer_patch(p) = 0._r8
            end if

            ! froot C and N
            if (abs(cs%froot_patch(p)) < ccrit) then
               pc = pc + cs%froot_patch(p)
               cs%froot_patch(p) = 0._r8
               pn = pn + ns%froot_patch(p)
               ns%froot_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%froot_patch(p)
                  c13cs%froot_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%froot_patch(p)
                  c14cs%froot_patch(p) = 0._r8
               endif

               pp = pp + ps%froot_patch(p)
               ps%froot_patch(p) = 0._r8
            end if

            ! froot storage C and N
            if (abs(cs%froot_storage_patch(p)) < ccrit) then
               pc = pc + cs%froot_storage_patch(p)
               cs%froot_storage_patch(p) = 0._r8
               pn = pn + ns%froot_storage_patch(p)
               ns%froot_storage_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%froot_storage_patch(p)
                  c13cs%froot_storage_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%froot_storage_patch(p)
                  c14cs%froot_storage_patch(p) = 0._r8
               endif

               pp = pp + ps%froot_storage_patch(p)
               ps%froot_storage_patch(p) = 0._r8
            end if

            ! froot transfer C and N
            if (abs(cs%froot_xfer_patch(p)) < ccrit) then
               pc = pc + cs%froot_xfer_patch(p)
               cs%froot_xfer_patch(p) = 0._r8
               pn = pn + ns%froot_xfer_patch(p)
               ns%froot_xfer_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%froot_xfer_patch(p)
                  c13cs%froot_xfer_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%froot_xfer_patch(p)
                  c14cs%froot_xfer_patch(p) = 0._r8
               endif

               pp = pp + ps%froot_xfer_patch(p)
               ps%froot_xfer_patch(p) = 0._r8
            end if

            if ( crop_prog .and. veg_pp%itype(p) >= nc3crop )then
               ! grain C and N
               if (abs(cs%grain_patch(p)) < ccrit) then
                  pc = pc + cs%grain_patch(p)
                  cs%grain_patch(p) = 0._r8
                  pn = pn + ns%grain_patch(p)
                  ns%grain_patch(p) = 0._r8
                  pp = pp + ps%grain_patch(p)
                  ps%grain_patch(p) = 0._r8
               end if

               ! grain storage C and N
               if (abs(cs%grain_storage_patch(p)) < ccrit) then
                  pc = pc + cs%grain_storage_patch(p)
                  cs%grain_storage_patch(p) = 0._r8
                  pn = pn + ns%grain_storage_patch(p)
                  ns%grain_storage_patch(p) = 0._r8
                  pp = pp + ps%grain_storage_patch(p)
                  ps%grain_storage_patch(p) = 0._r8
               end if

               ! grain transfer C and N
               if (abs(cs%grain_xfer_patch(p)) < ccrit) then
                  pc = pc + cs%grain_xfer_patch(p)
                  cs%grain_xfer_patch(p) = 0._r8
                  pn = pn + ns%grain_xfer_patch(p)
                  ns%grain_xfer_patch(p) = 0._r8
                  pp = pp + ps%grain_xfer_patch(p)
                  ps%grain_xfer_patch(p) = 0._r8
               end if
            end if

            ! livestem C and N
            if (abs(cs%livestem_patch(p)) < ccrit) then
               pc = pc + cs%livestem_patch(p)
               cs%livestem_patch(p) = 0._r8
               pn = pn + ns%livestem_patch(p)
               ns%livestem_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%livestem_patch(p)
                  c13cs%livestem_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%livestem_patch(p)
                  c14cs%livestem_patch(p) = 0._r8
               endif

               pp = pp + ps%livestem_patch(p)
               ps%livestem_patch(p) = 0._r8
            end if

            ! livestem storage C and N
            if (abs(cs%livestem_storage_patch(p)) < ccrit) then
               pc = pc + cs%livestem_storage_patch(p)
               cs%livestem_storage_patch(p) = 0._r8
               pn = pn + ns%livestem_storage_patch(p)
               ns%livestem_storage_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%livestem_storage_patch(p)
                  c13cs%livestem_storage_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%livestem_storage_patch(p)
                  c14cs%livestem_storage_patch(p) = 0._r8
               endif

               pp = pp + ps%livestem_storage_patch(p)
               ps%livestem_storage_patch(p) = 0._r8
            end if

            ! livestem transfer C and N
            if (abs(cs%livestem_xfer_patch(p)) < ccrit) then
               pc = pc + cs%livestem_xfer_patch(p)
               cs%livestem_xfer_patch(p) = 0._r8
               pn = pn + ns%livestem_xfer_patch(p)
               ns%livestem_xfer_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%livestem_xfer_patch(p)
                  c13cs%livestem_xfer_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%livestem_xfer_patch(p)
                  c14cs%livestem_xfer_patch(p) = 0._r8
               endif

               pp = pp + ps%livestem_xfer_patch(p)
               ps%livestem_xfer_patch(p) = 0._r8
            end if

            ! deadstem C and N
            if (abs(cs%deadstem_patch(p)) < ccrit) then
               pc = pc + cs%deadstem_patch(p)
               cs%deadstem_patch(p) = 0._r8
               pn = pn + ns%deadstem_patch(p)
               ns%deadstem_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%deadstem_patch(p)
                  c13cs%deadstem_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%deadstem_patch(p)
                  c14cs%deadstem_patch(p) = 0._r8
               endif

               pp = pp + ps%deadstem_patch(p)
               ps%deadstem_patch(p) = 0._r8
            end if

            ! deadstem storage C and N
            if (abs(cs%deadstem_storage_patch(p)) < ccrit) then
               pc = pc + cs%deadstem_storage_patch(p)
               cs%deadstem_storage_patch(p) = 0._r8
               pn = pn + ns%deadstem_storage_patch(p)
               ns%deadstem_storage_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%deadstem_storage_patch(p)
                  c13cs%deadstem_storage_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%deadstem_storage_patch(p)
                  c14cs%deadstem_storage_patch(p) = 0._r8
               endif

               pp = pp + ps%deadstem_storage_patch(p)
               ps%deadstem_storage_patch(p) = 0._r8
            end if

            ! deadstem transfer C and N
            if (abs(cs%deadstem_xfer_patch(p)) < ccrit) then
               pc = pc + cs%deadstem_xfer_patch(p)
               cs%deadstem_xfer_patch(p) = 0._r8
               pn = pn + ns%deadstem_xfer_patch(p)
               ns%deadstem_xfer_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%deadstem_xfer_patch(p)
                  c13cs%deadstem_xfer_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%deadstem_xfer_patch(p)
                  c14cs%deadstem_xfer_patch(p) = 0._r8
               endif

               pp = pp + ps%deadstem_xfer_patch(p)
               ps%deadstem_xfer_patch(p) = 0._r8
            end if

            ! livecroot C and N
            if (abs(cs%livecroot_patch(p)) < ccrit) then
               pc = pc + cs%livecroot_patch(p)
               cs%livecroot_patch(p) = 0._r8
               pn = pn + ns%livecroot_patch(p)
               ns%livecroot_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%livecroot_patch(p)
                  c13cs%livecroot_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%livecroot_patch(p)
                  c14cs%livecroot_patch(p) = 0._r8
               endif

               pp = pp + ps%livecroot_patch(p)
               ps%livecroot_patch(p) = 0._r8
            end if

            ! livecroot storage C and N
            if (abs(cs%livecroot_storage_patch(p)) < ccrit) then
               pc = pc + cs%livecroot_storage_patch(p)
               cs%livecroot_storage_patch(p) = 0._r8
               pn = pn + ns%livecroot_storage_patch(p)
               ns%livecroot_storage_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%livecroot_storage_patch(p)
                  c13cs%livecroot_storage_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%livecroot_storage_patch(p)
                  c14cs%livecroot_storage_patch(p) = 0._r8
               endif

               pp = pp + ps%livecroot_storage_patch(p)
               ps%livecroot_storage_patch(p) = 0._r8
            end if

            ! livecroot transfer C and N
            if (abs(cs%livecroot_xfer_patch(p)) < ccrit) then
               pc = pc + cs%livecroot_xfer_patch(p)
               cs%livecroot_xfer_patch(p) = 0._r8
               pn = pn + ns%livecroot_xfer_patch(p)
               ns%livecroot_xfer_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%livecroot_xfer_patch(p)
                  c13cs%livecroot_xfer_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%livecroot_xfer_patch(p)
                  c14cs%livecroot_xfer_patch(p) = 0._r8
               endif

               pp = pp + ps%livecroot_xfer_patch(p)
               ps%livecroot_xfer_patch(p) = 0._r8
            end if

            ! deadcroot C and N
            if (abs(cs%deadcroot_patch(p)) < ccrit) then
               pc = pc + cs%deadcroot_patch(p)
               cs%deadcroot_patch(p) = 0._r8
               pn = pn + ns%deadcroot_patch(p)
               ns%deadcroot_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%deadcroot_patch(p)
                  c13cs%deadcroot_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%deadcroot_patch(p)
                  c14cs%deadcroot_patch(p) = 0._r8
               endif

               pp = pp + ps%deadcroot_patch(p)
               ps%deadcroot_patch(p) = 0._r8
            end if

            ! deadcroot storage C and N
            if (abs(cs%deadcroot_storage_patch(p)) < ccrit) then
               pc = pc + cs%deadcroot_storage_patch(p)
               cs%deadcroot_storage_patch(p) = 0._r8
               pn = pn + ns%deadcroot_storage_patch(p)
               ns%deadcroot_storage_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%deadcroot_storage_patch(p)
                  c13cs%deadcroot_storage_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%deadcroot_storage_patch(p)
                  c14cs%deadcroot_storage_patch(p) = 0._r8
               endif

               pp = pp + ps%deadcroot_storage_patch(p)
               ps%deadcroot_storage_patch(p) = 0._r8
            end if

            ! deadcroot transfer C and N
            if (abs(cs%deadcroot_xfer_patch(p)) < ccrit) then
               pc = pc + cs%deadcroot_xfer_patch(p)
               cs%deadcroot_xfer_patch(p) = 0._r8
               pn = pn + ns%deadcroot_xfer_patch(p)
               ns%deadcroot_xfer_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%deadcroot_xfer_patch(p)
                  c13cs%deadcroot_xfer_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%deadcroot_xfer_patch(p)
                  c14cs%deadcroot_xfer_patch(p) = 0._r8
               endif

               pp = pp + ps%deadcroot_xfer_patch(p)
               ps%deadcroot_xfer_patch(p) = 0._r8
            end if

            ! gresp_storage (C only)
            if (abs(cs%gresp_storage_patch(p)) < ccrit) then
               pc = pc + cs%gresp_storage_patch(p)
               cs%gresp_storage_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%gresp_storage_patch(p)
                  c13cs%gresp_storage_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%gresp_storage_patch(p)
                  c14cs%gresp_storage_patch(p) = 0._r8
               endif
            end if

            ! gresp_xfer(c only)
            if (abs(cs%gresp_xfer_patch(p)) < ccrit) then
               pc = pc + cs%gresp_xfer_patch(p)
               cs%gresp_xfer_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%gresp_xfer_patch(p)
                  c13cs%gresp_xfer_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%gresp_xfer_patch(p)
                  c14cs%gresp_xfer_patch(p) = 0._r8
               endif
            end if

            ! cpool (C only)
            if (abs(cs%pool_patch(p)) < ccrit) then
               pc = pc + cs%pool_patch(p)
               cs%pool_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%pool_patch(p)
                  c13cs%pool_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%pool_patch(p)
                  c14cs%pool_patch(p) = 0._r8
               endif
            end if

            if ( crop_prog .and. veg_pp%itype(p) >= nc3crop )then
               ! xsmrpool (C only)
               if (abs(cs%xsmrpool_patch(p)) < ccrit) then
                  pc = pc + cs%xsmrpool_patch(p)
                  cs%xsmrpool_patch(p) = 0._r8
               end if
            end if

            ! retransn (N only)
            if (abs(ns%retransn_patch(p)) < ncrit) then
               pn = pn + ns%retransn_patch(p)
               ns%retransn_patch(p) = 0._r8
            end if

            ! retransp (P only)
            if (abs(ps%retransp_patch(p)) < pcrit) then
               pp = pp + ps%retransp_patch(p)
               ps%retransp_patch(p) = 0._r8
            end if

            ! npool (N only)
            if (abs(ns%pool_patch(p)) < ncrit) then
               pn = pn + ns%pool_patch(p)
               ns%pool_patch(p) = 0._r8
            end if

            ! ppool (P only)
            if (abs(ps%pool_patch(p)) < pcrit) then
               pp = pp + ps%pool_patch(p)
               ps%pool_patch(p) = 0._r8
            end if

            cs%veg_trunc_patch(p) = cs%veg_trunc_patch(p) + pc
            ns%veg_trunc_patch(p) = ns%veg_trunc_patch(p) + pn
            ps%veg_trunc_patch(p) = ps%veg_trunc_patch(p) + pp

            if ( use_c13 ) then
               c13cs%veg_trunc_patch(p) = c13cs%veg_trunc_patch(p) + pc13
            endif
            if ( use_c14 ) then
               c14cs%veg_trunc_patch(p) = c14cs%veg_trunc_patch(p) + pc14
            endif

         end do ! end of pft loop
      end if ! end of if(not.use_fates)

      if (.not. is_active_betr_bgc) then

         ! column loop
         do fc = 1,num_soilc
            c = filter_soilc(fc)

            do j = 1,nlevdecomp_full
               ! initialize the column-level C and N truncation terms
               cc = 0._r8
               if ( use_c13 ) cc13 = 0._r8
               if ( use_c14 ) cc14 = 0._r8
               cn = 0._r8

               ! do tests on state variables for precision control
               ! for linked C-N state variables, perform precision test on
               ! the C component, but truncate both C and N components


               ! all decomposing pools C and N
               do k = 1, ndecomp_pools

                  if (abs(cs%decomp_pools_vr_col(c,j,k)) < ccrit) then
                     cc = cc + cs%decomp_pools_vr_col(c,j,k)
                     cs%decomp_pools_vr_col(c,j,k) = 0._r8
                     if (.not.use_fates) then
                        cn = cn + ns%decomp_pools_vr_col(c,j,k)
                        ns%decomp_pools_vr_col(c,j,k) = 0._r8
                     endif
                     if ( use_c13 ) then
                        cc13 = cc13 + c13cs%decomp_pools_vr_col(c,j,k)
                        c13cs%decomp_pools_vr_col(c,j,k) = 0._r8
                     endif
                     if ( use_c14 ) then
                        cc14 = cc14 + c14cs%decomp_pools_vr_col(c,j,k)
                        c14cs%decomp_pools_vr_col(c,j,k) = 0._r8
                     endif
                  end if

               end do

               ! not doing precision control on soil mineral N, since it will
               ! be getting the N truncation flux anyway.

               cs%soil_trunc_vr_col(c,j) = cs%soil_trunc_vr_col(c,j) + cc
               if (.not.use_fates) then
                  ns%soil_trunc_vr_col(c,j) = ns%soil_trunc_vr_col(c,j) + cn
               endif
               if ( use_c13 ) then
                  c13cs%soil_trunc_vr_col(c,j) = c13cs%soil_trunc_vr_col(c,j) + cc13
               endif
               if ( use_c14 ) then
                  c14cs%soil_trunc_vr_col(c,j) = c14cs%soil_trunc_vr_col(c,j) + cc14
               endif
            end do

         end do   ! end of column loop

         if (use_nitrif_denitrif) then
            ! remove small negative perturbations for stability purposes, if any should arise.

            do fc = 1,num_soilc
               c = filter_soilc(fc)
               do j = 1,nlevdecomp_full
                  if (abs(ns%smin_no3_vr_col(c,j)) < ncrit/1e4_r8) then
                     if ( ns%smin_no3_vr_col(c,j)  < 0._r8 ) then
                        write(iulog, *) '-10^-12 < smin_no3 < 0. resetting to zero.'
                        write(iulog, *) 'smin_no3_vr_col(c,j), c, j: ', ns%smin_no3_vr_col(c,j), c, j
                        ns%smin_no3_vr_col(c,j) = 0._r8
                     endif
                  end if
                  if (abs(ns%smin_nh4_vr_col(c,j)) < ncrit/1e4_r8) then
                     if ( ns%smin_nh4_vr_col(c,j)  < 0._r8 ) then
                        write(iulog, *) '-10^-12 < smin_nh4 < 0. resetting to zero.'
                        write(iulog, *) 'smin_nh4_vr_col(c,j), c, j: ', ns%smin_nh4_vr_col(c,j), c, j
                        ns%smin_nh4_vr_col(c,j) = 0._r8
                     endif
                  end if
               end do
            end do
         endif

         if (nu_com .eq. 'ECA') then
            ! decompose P pool adjust according to C pool
            !do fc = 1,num_soilc
            !   c = filter_soilc(fc)
            !   do j = 1,nlevdecomp_full
            !      cp_eca = 0.0_r8
            !      do l = 1,ndecomp_pools
            !         if (abs(cs%decomp_pools_vr_col(c,j,k)) < ccrit) then
            !            if (.not.use_fates) then
            !               cp_eca = cp_eca + ps%decomp_pools_vr_col(c,j,k)
            !               ps%decomp_pools_vr_col(c,j,k) = 0._r8
            !            endif
            !         endif
            !      end do
            !      ps%soil_trunc_vr_col(c,j) = ps%soil_trunc_vr_col(c,j) + cp_eca
            !   end do
            !end do

            ! fix soil CN ratio drift (normally < 0.01% drift)
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               do j = 1,nlevdecomp_full
                  cn_eca = 0.0_r8
                  do l = 1,ndecomp_pools
                     if ( cs%decomp_pools_vr_col(c,j,l) > 0.0_r8 .and.  &
                          abs(cs%decomp_pools_vr_col(c,j,l) / ns%decomp_pools_vr_col(c,j,l) - initial_cn_ratio(l) ) > 1.0e-3_r8 &
                          .and. (.not. floating_cn_ratio_decomp_pools(l)) ) then
                        cn_eca = cn_eca - ( cs%decomp_pools_vr_col(c,j,l) / initial_cn_ratio(l) - ns%decomp_pools_vr_col(c,j,l) )
                        ns%decomp_pools_vr_col(c,j,l) = cs%decomp_pools_vr_col(c,j,l) / initial_cn_ratio(l)
                     end if
                  end do
                  ns%soil_trunc_vr_col(c,j) = ns%soil_trunc_vr_col(c,j) + cn_eca
               end do
             end do

            ! remove small negative perturbations for stability purposes, if any should arise in N,P pools
            ! for floating CN, CP ratio pools
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               do j = 1,nlevdecomp_full

                  cn_eca = 0.0_r8
                  cp_eca = 0.0_r8
                  do l = 1,ndecomp_pools
                     if ( ns%decomp_pools_vr_col(c,j,l) < 0.0_r8 .and. floating_cn_ratio_decomp_pools(l) ) then
                        if ( abs(ns%decomp_pools_vr_col(c,j,l))  < ncrit ) then
                           cn_eca = cn_eca - ncrit + ns%decomp_pools_vr_col(c,j,l)
                           ns%decomp_pools_vr_col(c,j,l) = ncrit
                        else
                           write(iulog, "(A,2I8,E8.1)") 'error decomp_npools is negative: ',j,l,ns%decomp_pools_vr_col(c,j,l)
                           call endrun(msg=errMsg(__FILE__, __LINE__))
                        end if
                     end if
                     if ( ps%decomp_pools_vr_col(c,j,l)  < 0.0_r8 .and. floating_cp_ratio_decomp_pools(l) ) then
                        if ( abs(ps%decomp_pools_vr_col(c,j,l))  < ncrit/1e4_r8 ) then
                           cp_eca = cp_eca - ncrit/1e4_r8 + ps%decomp_pools_vr_col(c,j,l)
                           ps%decomp_pools_vr_col(c,j,l) = ncrit/1e4_r8
                         else 
                           write(iulog, "(A,2I8,E8.1)") 'error decomp_ppools is negative: ',j,l,ps%decomp_pools_vr_col(c,j,l)
                           call endrun(msg=errMsg(__FILE__, __LINE__))
                         end if
                     end if

                  end do

                  ns%soil_trunc_vr_col(c,j) = ns%soil_trunc_vr_col(c,j) + cn_eca
                  ps%soil_trunc_vr_col(c,j) = ps%soil_trunc_vr_col(c,j) + cp_eca

               end do
            end do

            do fp = 1,num_soilp
               p = filter_soilp(fp)
               if (ns%retransn_patch(p) < 0._r8) then
                  write(iulog, *) 'error retransn_patch is negative: ',p
                  write(iulog, *) 'retransn_patch: ', ns%retransn_patch(p)
                  call endrun(msg=errMsg(__FILE__, __LINE__))
               end if
               if (ns%pool_patch(p) < 0._r8) then
                  write(iulog, *) 'error npool_patch is negative: ',p
                  write(iulog, *) 'npool_patch: ', ns%pool_patch(p)
                  call endrun(msg=errMsg(__FILE__, __LINE__))
               end if
               if (ps%retransp_patch(p) < 0._r8) then
                  write(iulog, *) 'error retransp_patch is negative: ',p
                  write(iulog, *) 'retransp_patch: ', ps%retransp_patch(p)
                  call endrun(msg=errMsg(__FILE__, __LINE__))
               end if
               if (ps%pool_patch(p) < 0._r8) then
                  write(iulog, *) 'error ppool_patch is negative: ',p
                  write(iulog, *) 'ppool_patch: ', ps%pool_patch(p)
                  call endrun(msg=errMsg(__FILE__, __LINE__))
               end if
            end do

         endif

      endif ! if (.not. is_active_betr_bgc)

    end associate

 end subroutine CNPrecisionControl

end module CNPrecisionControlMod
