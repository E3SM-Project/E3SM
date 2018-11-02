module PrecisionControlMod

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
  public:: PrecisionControl
  !----------------------------------------------------------------------- 

contains

  !-----------------------------------------------------------------------
  subroutine PrecisionControl(num_soilc, filter_soilc, num_soilp, filter_soilp, &
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

    ! carbonstate_vars%ctrunc_vr_col                 Output:  [real(r8) (:,:)   ]  (gC/m3) column-level sink for C truncation      
    ! carbonstate_vars%decomp_cpools_vr_col          Output:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
    ! carbonstate_vars%cpool_patch                   Output:  [real(r8) (:)     ]  (gC/m2) temporary photosynthate C pool            
    ! carbonstate_vars%deadcrootc_patch              Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C                        
    ! carbonstate_vars%deadcrootc_storage_patch      Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C storage                
    ! carbonstate_vars%deadcrootc_xfer_patch         Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C transfer               
    ! carbonstate_vars%deadstemc_patch               Output:  [real(r8) (:)     ]  (gC/m2) dead stem C                               
    ! carbonstate_vars%deadstemc_storage_patch       Output:  [real(r8) (:)     ]  (gC/m2) dead stem C storage                       
    ! carbonstate_vars%deadstemc_xfer_patch          Output:  [real(r8) (:)     ]  (gC/m2) dead stem C transfer                      
    ! carbonstate_vars%frootc_patch                  Output:  [real(r8) (:)     ]  (gC/m2) fine root C                               
    ! carbonstate_vars%frootc_storage_patch          Output:  [real(r8) (:)     ]  (gC/m2) fine root C storage                       
    ! carbonstate_vars%frootc_xfer_patch             Output:  [real(r8) (:)     ]  (gC/m2) fine root C transfer                      
    ! carbonstate_vars%gresp_storage_patch           Output:  [real(r8) (:)     ]  (gC/m2) growth respiration storage                
    ! carbonstate_vars%gresp_xfer_patch              Output:  [real(r8) (:)     ]  (gC/m2) growth respiration transfer               
    ! carbonstate_vars%leafc_patch                   Output:  [real(r8) (:)     ]  (gC/m2) leaf C                                    
    ! carbonstate_vars%leafc_storage_patch           Output:  [real(r8) (:)     ]  (gC/m2) leaf C storage                            
    ! carbonstate_vars%leafc_xfer_patch              Output:  [real(r8) (:)     ]  (gC/m2) leaf C transfer                           
    ! carbonstate_vars%livecrootc_patch              Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C                        
    ! carbonstate_vars%livecrootc_storage_patch      Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C storage                
    ! carbonstate_vars%livecrootc_xfer_patch         Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C transfer               
    ! carbonstate_vars%livestemc_patch               Output:  [real(r8) (:)     ]  (gC/m2) live stem C                               
    ! carbonstate_vars%livestemc_storage_patch       Output:  [real(r8) (:)     ]  (gC/m2) live stem C storage                       
    ! carbonstate_vars%livestemc_xfer_patch          Output:  [real(r8) (:)     ]  (gC/m2) live stem C transfer                      
    ! carbonstate_vars%ctrunc_patch                  Output:  [real(r8) (:)     ]  (gC/m2) patch-level sink for C truncation           
    ! carbonstate_vars%xsmrpool_patch                Output:  [real(r8) (:)     ]  (gC/m2) execss maint resp C pool                  
    ! carbonstate_vars%grainc_patch                  Output:  [real(r8) (:)     ]  (gC/m2) grain C                                   
    ! carbonstate_vars%grainc_storage_patch          Output:  [real(r8) (:)     ]  (gC/m2) grain C storage                           
    ! carbonstate_vars%grainc_xfer_patch             Output:  [real(r8) (:)     ]  (gC/m2) grain C transfer                          
    
    ! c13_carbonstate_vars%ctrunc_vr_col             Output:  [real(r8) (:,:)   ]  (gC/m3) column-level sink for C truncation      
    ! c13_carbonstate_vars%decomp_cpools_vr_col      Output:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
    ! c13_carbonstate_vars%cpool_patch               Output:  [real(r8) (:)     ]  (gC/m2) temporary photosynthate C pool            
    ! c13_carbonstate_vars%deadcrootc_patch          Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C                        
    ! c13_carbonstate_vars%deadcrootc_storage_patch  Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C storage                
    ! c13_carbonstate_vars%deadcrootc_xfer_patch     Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C transfer               
    ! c13_carbonstate_vars%deadstemc_patch           Output:  [real(r8) (:)     ]  (gC/m2) dead stem C                               
    ! c13_carbonstate_vars%deadstemc_storage_patch   Output:  [real(r8) (:)     ]  (gC/m2) dead stem C storage                       
    ! c13_carbonstate_vars%deadstemc_xfer_patch      Output:  [real(r8) (:)     ]  (gC/m2) dead stem C transfer                      
    ! c13_carbonstate_vars%frootc_patch              Output:  [real(r8) (:)     ]  (gC/m2) fine root C                               
    ! c13_carbonstate_vars%frootc_storage_patch      Output:  [real(r8) (:)     ]  (gC/m2) fine root C storage                       
    ! c13_carbonstate_vars%frootc_xfer_patch         Output:  [real(r8) (:)     ]  (gC/m2) fine root C transfer                      
    ! c13_carbonstate_vars%gresp_storage_patch       Output:  [real(r8) (:)     ]  (gC/m2) growth respiration storage                
    ! c13_carbonstate_vars%gresp_xfer_patch          Output:  [real(r8) (:)     ]  (gC/m2) growth respiration transfer               
    ! c13_carbonstate_vars%leafc_patch               Output:  [real(r8) (:)     ]  (gC/m2) leaf C                                    
    ! c13_carbonstate_vars%leafc_storage_patch       Output:  [real(r8) (:)     ]  (gC/m2) leaf C storage                            
    ! c13_carbonstate_vars%leafc_xfer_patch          Output:  [real(r8) (:)     ]  (gC/m2) leaf C transfer                           
    ! c13_carbonstate_vars%livecrootc_patch          Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C                        
    ! c13_carbonstate_vars%livecrootc_storage_patch  Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C storage                
    ! c13_carbonstate_vars%livecrootc_xfer_patch     Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C transfer               
    ! c13_carbonstate_vars%livestemc_patch           Output:  [real(r8) (:)     ]  (gC/m2) live stem C                               
    ! c13_carbonstate_vars%livestemc_storage_patch   Output:  [real(r8) (:)     ]  (gC/m2) live stem C storage                       
    ! c13_carbonstate_vars%livestemc_xfer_patch      Output:  [real(r8) (:)     ]  (gC/m2) live stem C transfer                      
    ! c13_carbonstate_vars%ctrunc_patch              Output:  [real(r8) (:)     ]  (gC/m2) patch-level sink for C truncation           
    
    ! c14_carbonstate_vars%ctrunc_vr_col             Output:  [real(r8) (:,:)   ]  (gC/m3) column-level sink for C truncation      
    ! c14_carbonstate_vars%decomp_cpools_vr_col      Output:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
    ! c14_carbonstate_vars%cpool_patch               Output:  [real(r8) (:)     ]  (gC/m2) temporary photosynthate C pool            
    ! c14_carbonstate_vars%deadcrootc_patch          Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C                        
    ! c14_carbonstate_vars%deadcrootc_storage_patch  Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C storage                
    ! c14_carbonstate_vars%deadcrootc_xfer_patch     Output:  [real(r8) (:)     ]  (gC/m2) dead coarse root C transfer               
    ! c14_carbonstate_vars%deadstemc_patch           Output:  [real(r8) (:)     ]  (gC/m2) dead stem C                               
    ! c14_carbonstate_vars%deadstemc_storage_patch   Output:  [real(r8) (:)     ]  (gC/m2) dead stem C storage                       
    ! c14_carbonstate_vars%deadstemc_xfer_patch      Output:  [real(r8) (:)     ]  (gC/m2) dead stem C transfer                      
    ! c14_carbonstate_vars%frootc_patch              Output:  [real(r8) (:)     ]  (gC/m2) fine root C                               
    ! c14_carbonstate_vars%frootc_storage_patch      Output:  [real(r8) (:)     ]  (gC/m2) fine root C storage                       
    ! c14_carbonstate_vars%frootc_xfer_patch         Output:  [real(r8) (:)     ]  (gC/m2) fine root C transfer                      
    ! c14_carbonstate_vars%gresp_storage_patch       Output:  [real(r8) (:)     ]  (gC/m2) growth respiration storage                
    ! c14_carbonstate_vars%gresp_xfer_patch          Output:  [real(r8) (:)     ]  (gC/m2) growth respiration transfer               
    ! c14_carbonstate_vars%leafc_patch               Output:  [real(r8) (:)     ]  (gC/m2) leaf C                                    
    ! c14_carbonstate_vars%leafc_storage_patch       Output:  [real(r8) (:)     ]  (gC/m2) leaf C storage                            
    ! c14_carbonstate_vars%leafc_xfer_patch          Output:  [real(r8) (:)     ]  (gC/m2) leaf C transfer                           
    ! c14_carbonstate_vars%livecrootc_patch          Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C                        
    ! c14_carbonstate_vars%livecrootc_storage_patch  Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C storage                
    ! c14_carbonstate_vars%livecrootc_xfer_patch     Output:  [real(r8) (:)     ]  (gC/m2) live coarse root C transfer               
    ! c14_carbonstate_vars%livestemc_patch           Output:  [real(r8) (:)     ]  (gC/m2) live stem C                               
    ! c14_carbonstate_vars%livestemc_storage_patch   Output:  [real(r8) (:)     ]  (gC/m2) live stem C storage                       
    ! c14_carbonstate_vars%livestemc_xfer_patch      Output:  [real(r8) (:)     ]  (gC/m2) live stem C transfer                      
    ! c14_carbonstate_vars%ctrunc_patch              Output:  [real(r8) (:)     ]  (gC/m2) patch-level sink for C truncation           
    
    ! nitrogenstate_vars%ntrunc_vr_col               Output:  [real(r8) (:,:)   ]  (gN/m3) column-level sink for N truncation      
    ! nitrogenstate_vars%decomp_npools_vr_col        Output:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
    ! nitrogenstate_vars%deadcrootn_patch            Output:  [real(r8) (:)     ]  (gN/m2) dead coarse root N                        
    ! nitrogenstate_vars%deadcrootn_storage_patch    Output:  [real(r8) (:)     ]  (gN/m2) dead coarse root N storage                
    ! nitrogenstate_vars%deadcrootn_xfer_patch       Output:  [real(r8) (:)     ]  (gN/m2) dead coarse root N transfer               
    ! nitrogenstate_vars%deadstemn_patch             Output:  [real(r8) (:)     ]  (gN/m2) dead stem N                               
    ! nitrogenstate_vars%deadstemn_storage_patch     Output:  [real(r8) (:)     ]  (gN/m2) dead stem N storage                       
    ! nitrogenstate_vars%deadstemn_xfer_patch        Output:  [real(r8) (:)     ]  (gN/m2) dead stem N transfer                      
    ! nitrogenstate_vars%frootn_patch                Output:  [real(r8) (:)     ]  (gN/m2) fine root N                               
    ! nitrogenstate_vars%frootn_storage_patch        Output:  [real(r8) (:)     ]  (gN/m2) fine root N storage                       
    ! nitrogenstate_vars%frootn_xfer_patch           Output:  [real(r8) (:)     ]  (gN/m2) fine root N transfer                      
    ! nitrogenstate_vars%leafn_patch                 Output:  [real(r8) (:)     ]  (gN/m2) leaf N                                    
    ! nitrogenstate_vars%leafn_storage_patch         Output:  [real(r8) (:)     ]  (gN/m2) leaf N storage                            
    ! nitrogenstate_vars%leafn_xfer_patch            Output:  [real(r8) (:)     ]  (gN/m2) leaf N transfer                           
    ! nitrogenstate_vars%livecrootn_patch            Output:  [real(r8) (:)     ]  (gN/m2) live coarse root N                        
    ! nitrogenstate_vars%livecrootn_storage_patch    Output:  [real(r8) (:)     ]  (gN/m2) live coarse root N storage                
    ! nitrogenstate_vars%livecrootn_xfer_patch       Output:  [real(r8) (:)     ]  (gN/m2) live coarse root N transfer               
    ! nitrogenstate_vars%grainn_patch                Output:  [real(r8) (:)     ]  (gC/m2) grain N                                   
    ! nitrogenstate_vars%grainn_storage_patch        Output:  [real(r8) (:)     ]  (gC/m2) grain N storage                           
    ! nitrogenstate_vars%grainn_xfer_patch           Output:  [real(r8) (:)     ]  (gC/m2) grain N transfer                          
    ! nitrogenstate_vars%livestemn_patch             Output:  [real(r8) (:)     ]  (gN/m2) live stem N                               
    ! nitrogenstate_vars%livestemn_storage_patch     Output:  [real(r8) (:)     ]  (gN/m2) live stem N storage                       
    ! nitrogenstate_vars%livestemn_xfer_patch        Output:  [real(r8) (:)     ]  (gN/m2) live stem N transfer                      
    ! nitrogenstate_vars%npool_patch                 Output:  [real(r8) (:)     ]  (gN/m2) temporary plant N pool                    
    ! nitrogenstate_vars%ntrunc_patch                Output:  [real(r8) (:)     ]  (gN/m2) patch-level sink for N truncation           
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
            if (abs(cs%leafc_patch(p)) < ccrit) then
               pc = pc + cs%leafc_patch(p)
               cs%leafc_patch(p) = 0._r8
               pn = pn + ns%leafn_patch(p)
               ns%leafn_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%leafc_patch(p)
                  c13cs%leafc_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%leafc_patch(p)
                  c14cs%leafc_patch(p) = 0._r8
               endif

               pp = pp + ps%leafp_patch(p)
               ps%leafp_patch(p) = 0._r8
            end if

            ! leaf storage C and N
            if (abs(cs%leafc_storage_patch(p)) < ccrit) then
               pc = pc + cs%leafc_storage_patch(p)
               cs%leafc_storage_patch(p) = 0._r8
               pn = pn + ns%leafn_storage_patch(p)
               ns%leafn_storage_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%leafc_storage_patch(p)
                  c13cs%leafc_storage_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%leafc_storage_patch(p)
                  c14cs%leafc_storage_patch(p) = 0._r8
               endif

               pp = pp + ps%leafp_storage_patch(p)
               ps%leafp_storage_patch(p) = 0._r8
            end if

            ! leaf transfer C and N
            if (abs(cs%leafc_xfer_patch(p)) < ccrit) then
               pc = pc + cs%leafc_xfer_patch(p)
               cs%leafc_xfer_patch(p) = 0._r8
               pn = pn + ns%leafn_xfer_patch(p)
               ns%leafn_xfer_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%leafc_xfer_patch(p)
                  c13cs%leafc_xfer_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%leafc_xfer_patch(p)
                  c14cs%leafc_xfer_patch(p) = 0._r8
               endif

               pp = pp + ps%leafp_xfer_patch(p)
               ps%leafp_xfer_patch(p) = 0._r8
            end if

            ! froot C and N
            if (abs(cs%frootc_patch(p)) < ccrit) then
               pc = pc + cs%frootc_patch(p)
               cs%frootc_patch(p) = 0._r8
               pn = pn + ns%frootn_patch(p)
               ns%frootn_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%frootc_patch(p)
                  c13cs%frootc_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%frootc_patch(p)
                  c14cs%frootc_patch(p) = 0._r8
               endif

               pp = pp + ps%frootp_patch(p)
               ps%frootp_patch(p) = 0._r8
            end if

            ! froot storage C and N
            if (abs(cs%frootc_storage_patch(p)) < ccrit) then
               pc = pc + cs%frootc_storage_patch(p)
               cs%frootc_storage_patch(p) = 0._r8
               pn = pn + ns%frootn_storage_patch(p)
               ns%frootn_storage_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%frootc_storage_patch(p)
                  c13cs%frootc_storage_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%frootc_storage_patch(p)
                  c14cs%frootc_storage_patch(p) = 0._r8
               endif

               pp = pp + ps%frootp_storage_patch(p)
               ps%frootp_storage_patch(p) = 0._r8
            end if

            ! froot transfer C and N
            if (abs(cs%frootc_xfer_patch(p)) < ccrit) then
               pc = pc + cs%frootc_xfer_patch(p)
               cs%frootc_xfer_patch(p) = 0._r8
               pn = pn + ns%frootn_xfer_patch(p)
               ns%frootn_xfer_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%frootc_xfer_patch(p)
                  c13cs%frootc_xfer_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%frootc_xfer_patch(p)
                  c14cs%frootc_xfer_patch(p) = 0._r8
               endif

               pp = pp + ps%frootp_xfer_patch(p)
               ps%frootp_xfer_patch(p) = 0._r8
            end if

            if ( crop_prog .and. veg_pp%itype(p) >= nc3crop )then
               ! grain C and N
               if (abs(cs%grainc_patch(p)) < ccrit) then
                  pc = pc + cs%grainc_patch(p)
                  cs%grainc_patch(p) = 0._r8
                  pn = pn + ns%grainn_patch(p)
                  ns%grainn_patch(p) = 0._r8
                  pp = pp + ps%grainp_patch(p)
                  ps%grainp_patch(p) = 0._r8
               end if

               ! grain storage C and N
               if (abs(cs%grainc_storage_patch(p)) < ccrit) then
                  pc = pc + cs%grainc_storage_patch(p)
                  cs%grainc_storage_patch(p) = 0._r8
                  pn = pn + ns%grainn_storage_patch(p)
                  ns%grainn_storage_patch(p) = 0._r8
                  pp = pp + ps%grainp_storage_patch(p)
                  ps%grainp_storage_patch(p) = 0._r8
               end if

               ! grain transfer C and N
               if (abs(cs%grainc_xfer_patch(p)) < ccrit) then
                  pc = pc + cs%grainc_xfer_patch(p)
                  cs%grainc_xfer_patch(p) = 0._r8
                  pn = pn + ns%grainn_xfer_patch(p)
                  ns%grainn_xfer_patch(p) = 0._r8
                  pp = pp + ps%grainp_xfer_patch(p)
                  ps%grainp_xfer_patch(p) = 0._r8
               end if
            end if

            ! livestem C and N
            if (abs(cs%livestemc_patch(p)) < ccrit) then
               pc = pc + cs%livestemc_patch(p)
               cs%livestemc_patch(p) = 0._r8
               pn = pn + ns%livestemn_patch(p)
               ns%livestemn_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%livestemc_patch(p)
                  c13cs%livestemc_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%livestemc_patch(p)
                  c14cs%livestemc_patch(p) = 0._r8
               endif

               pp = pp + ps%livestemp_patch(p)
               ps%livestemp_patch(p) = 0._r8
            end if

            ! livestem storage C and N
            if (abs(cs%livestemc_storage_patch(p)) < ccrit) then
               pc = pc + cs%livestemc_storage_patch(p)
               cs%livestemc_storage_patch(p) = 0._r8
               pn = pn + ns%livestemn_storage_patch(p)
               ns%livestemn_storage_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%livestemc_storage_patch(p)
                  c13cs%livestemc_storage_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%livestemc_storage_patch(p)
                  c14cs%livestemc_storage_patch(p) = 0._r8
               endif

               pp = pp + ps%livestemp_storage_patch(p)
               ps%livestemp_storage_patch(p) = 0._r8
            end if

            ! livestem transfer C and N
            if (abs(cs%livestemc_xfer_patch(p)) < ccrit) then
               pc = pc + cs%livestemc_xfer_patch(p)
               cs%livestemc_xfer_patch(p) = 0._r8
               pn = pn + ns%livestemn_xfer_patch(p)
               ns%livestemn_xfer_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%livestemc_xfer_patch(p)
                  c13cs%livestemc_xfer_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%livestemc_xfer_patch(p)
                  c14cs%livestemc_xfer_patch(p) = 0._r8
               endif

               pp = pp + ps%livestemp_xfer_patch(p)
               ps%livestemp_xfer_patch(p) = 0._r8
            end if

            ! deadstem C and N
            if (abs(cs%deadstemc_patch(p)) < ccrit) then
               pc = pc + cs%deadstemc_patch(p)
               cs%deadstemc_patch(p) = 0._r8
               pn = pn + ns%deadstemn_patch(p)
               ns%deadstemn_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%deadstemc_patch(p)
                  c13cs%deadstemc_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%deadstemc_patch(p)
                  c14cs%deadstemc_patch(p) = 0._r8
               endif

               pp = pp + ps%deadstemp_patch(p)
               ps%deadstemp_patch(p) = 0._r8
            end if

            ! deadstem storage C and N
            if (abs(cs%deadstemc_storage_patch(p)) < ccrit) then
               pc = pc + cs%deadstemc_storage_patch(p)
               cs%deadstemc_storage_patch(p) = 0._r8
               pn = pn + ns%deadstemn_storage_patch(p)
               ns%deadstemn_storage_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%deadstemc_storage_patch(p)
                  c13cs%deadstemc_storage_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%deadstemc_storage_patch(p)
                  c14cs%deadstemc_storage_patch(p) = 0._r8
               endif

               pp = pp + ps%deadstemp_storage_patch(p)
               ps%deadstemp_storage_patch(p) = 0._r8
            end if

            ! deadstem transfer C and N
            if (abs(cs%deadstemc_xfer_patch(p)) < ccrit) then
               pc = pc + cs%deadstemc_xfer_patch(p)
               cs%deadstemc_xfer_patch(p) = 0._r8
               pn = pn + ns%deadstemn_xfer_patch(p)
               ns%deadstemn_xfer_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%deadstemc_xfer_patch(p)
                  c13cs%deadstemc_xfer_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%deadstemc_xfer_patch(p)
                  c14cs%deadstemc_xfer_patch(p) = 0._r8
               endif

               pp = pp + ps%deadstemp_xfer_patch(p)
               ps%deadstemp_xfer_patch(p) = 0._r8
            end if

            ! livecroot C and N
            if (abs(cs%livecrootc_patch(p)) < ccrit) then
               pc = pc + cs%livecrootc_patch(p)
               cs%livecrootc_patch(p) = 0._r8
               pn = pn + ns%livecrootn_patch(p)
               ns%livecrootn_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%livecrootc_patch(p)
                  c13cs%livecrootc_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%livecrootc_patch(p)
                  c14cs%livecrootc_patch(p) = 0._r8
               endif

               pp = pp + ps%livecrootp_patch(p)
               ps%livecrootp_patch(p) = 0._r8
            end if

            ! livecroot storage C and N
            if (abs(cs%livecrootc_storage_patch(p)) < ccrit) then
               pc = pc + cs%livecrootc_storage_patch(p)
               cs%livecrootc_storage_patch(p) = 0._r8
               pn = pn + ns%livecrootn_storage_patch(p)
               ns%livecrootn_storage_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%livecrootc_storage_patch(p)
                  c13cs%livecrootc_storage_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%livecrootc_storage_patch(p)
                  c14cs%livecrootc_storage_patch(p) = 0._r8
               endif

               pp = pp + ps%livecrootp_storage_patch(p)
               ps%livecrootp_storage_patch(p) = 0._r8
            end if

            ! livecroot transfer C and N
            if (abs(cs%livecrootc_xfer_patch(p)) < ccrit) then
               pc = pc + cs%livecrootc_xfer_patch(p)
               cs%livecrootc_xfer_patch(p) = 0._r8
               pn = pn + ns%livecrootn_xfer_patch(p)
               ns%livecrootn_xfer_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%livecrootc_xfer_patch(p)
                  c13cs%livecrootc_xfer_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%livecrootc_xfer_patch(p)
                  c14cs%livecrootc_xfer_patch(p) = 0._r8
               endif

               pp = pp + ps%livecrootp_xfer_patch(p)
               ps%livecrootp_xfer_patch(p) = 0._r8
            end if

            ! deadcroot C and N
            if (abs(cs%deadcrootc_patch(p)) < ccrit) then
               pc = pc + cs%deadcrootc_patch(p)
               cs%deadcrootc_patch(p) = 0._r8
               pn = pn + ns%deadcrootn_patch(p)
               ns%deadcrootn_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%deadcrootc_patch(p)
                  c13cs%deadcrootc_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%deadcrootc_patch(p)
                  c14cs%deadcrootc_patch(p) = 0._r8
               endif

               pp = pp + ps%deadcrootp_patch(p)
               ps%deadcrootp_patch(p) = 0._r8
            end if

            ! deadcroot storage C and N
            if (abs(cs%deadcrootc_storage_patch(p)) < ccrit) then
               pc = pc + cs%deadcrootc_storage_patch(p)
               cs%deadcrootc_storage_patch(p) = 0._r8
               pn = pn + ns%deadcrootn_storage_patch(p)
               ns%deadcrootn_storage_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%deadcrootc_storage_patch(p)
                  c13cs%deadcrootc_storage_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%deadcrootc_storage_patch(p)
                  c14cs%deadcrootc_storage_patch(p) = 0._r8
               endif

               pp = pp + ps%deadcrootp_storage_patch(p)
               ps%deadcrootp_storage_patch(p) = 0._r8
            end if

            ! deadcroot transfer C and N
            if (abs(cs%deadcrootc_xfer_patch(p)) < ccrit) then
               pc = pc + cs%deadcrootc_xfer_patch(p)
               cs%deadcrootc_xfer_patch(p) = 0._r8
               pn = pn + ns%deadcrootn_xfer_patch(p)
               ns%deadcrootn_xfer_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%deadcrootc_xfer_patch(p)
                  c13cs%deadcrootc_xfer_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%deadcrootc_xfer_patch(p)
                  c14cs%deadcrootc_xfer_patch(p) = 0._r8
               endif

               pp = pp + ps%deadcrootp_xfer_patch(p)
               ps%deadcrootp_xfer_patch(p) = 0._r8
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
            if (abs(cs%cpool_patch(p)) < ccrit) then
               pc = pc + cs%cpool_patch(p)
               cs%cpool_patch(p) = 0._r8
               if ( use_c13 ) then
                  pc13 = pc13 + c13cs%cpool_patch(p)
                  c13cs%cpool_patch(p) = 0._r8
               endif
               if ( use_c14 ) then
                  pc14 = pc14 + c14cs%cpool_patch(p)
                  c14cs%cpool_patch(p) = 0._r8
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
            if (abs(ns%npool_patch(p)) < ncrit) then
               pn = pn + ns%npool_patch(p)
               ns%npool_patch(p) = 0._r8
            end if

            ! ppool (P only)
            if (abs(ps%ppool_patch(p)) < pcrit) then
               pp = pp + ps%ppool_patch(p)
               ps%ppool_patch(p) = 0._r8
            end if

            cs%ctrunc_patch(p) = cs%ctrunc_patch(p) + pc
            ns%ntrunc_patch(p) = ns%ntrunc_patch(p) + pn
            ps%ptrunc_patch(p) = ps%ptrunc_patch(p) + pp

            if ( use_c13 ) then
               c13cs%ctrunc_patch(p) = c13cs%ctrunc_patch(p) + pc13
            endif
            if ( use_c14 ) then
               c14cs%ctrunc_patch(p) = c14cs%ctrunc_patch(p) + pc14
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

                  if (abs(cs%decomp_cpools_vr_col(c,j,k)) < ccrit) then
                     cc = cc + cs%decomp_cpools_vr_col(c,j,k)
                     cs%decomp_cpools_vr_col(c,j,k) = 0._r8
                     if (.not.use_fates) then
                        cn = cn + ns%decomp_npools_vr_col(c,j,k)
                        ns%decomp_npools_vr_col(c,j,k) = 0._r8
                     endif
                     if ( use_c13 ) then
                        cc13 = cc13 + c13cs%decomp_cpools_vr_col(c,j,k)
                        c13cs%decomp_cpools_vr_col(c,j,k) = 0._r8
                     endif
                     if ( use_c14 ) then
                        cc14 = cc14 + c14cs%decomp_cpools_vr_col(c,j,k)
                        c14cs%decomp_cpools_vr_col(c,j,k) = 0._r8
                     endif
                  end if

               end do

               ! not doing precision control on soil mineral N, since it will
               ! be getting the N truncation flux anyway.

               cs%ctrunc_vr_col(c,j) = cs%ctrunc_vr_col(c,j) + cc
               if (.not.use_fates) then
                  ns%ntrunc_vr_col(c,j) = ns%ntrunc_vr_col(c,j) + cn
               endif
               if ( use_c13 ) then
                  c13cs%ctrunc_vr_col(c,j) = c13cs%ctrunc_vr_col(c,j) + cc13
               endif
               if ( use_c14 ) then
                  c14cs%ctrunc_vr_col(c,j) = c14cs%ctrunc_vr_col(c,j) + cc14
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
            !         if (abs(cs%decomp_cpools_vr_col(c,j,k)) < ccrit) then
            !            if (.not.use_fates) then
            !               cp_eca = cp_eca + ps%decomp_ppools_vr_col(c,j,k)
            !               ps%decomp_ppools_vr_col(c,j,k) = 0._r8
            !            endif
            !         endif
            !      end do
            !      ps%ptrunc_vr_col(c,j) = ps%ptrunc_vr_col(c,j) + cp_eca
            !   end do
            !end do

            ! fix soil CN ratio drift (normally < 0.01% drift)
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               do j = 1,nlevdecomp_full
                  cn_eca = 0.0_r8
                  do l = 1,ndecomp_pools
                     if ( cs%decomp_cpools_vr_col(c,j,l) > 0.0_r8 .and.  &
                          abs(cs%decomp_cpools_vr_col(c,j,l) / ns%decomp_npools_vr_col(c,j,l) - initial_cn_ratio(l) ) > 1.0e-3_r8 &
                          .and. (.not. floating_cn_ratio_decomp_pools(l)) ) then
                        cn_eca = cn_eca - ( cs%decomp_cpools_vr_col(c,j,l) / initial_cn_ratio(l) - ns%decomp_npools_vr_col(c,j,l) )
                        ns%decomp_npools_vr_col(c,j,l) = cs%decomp_cpools_vr_col(c,j,l) / initial_cn_ratio(l)
                     end if
                  end do
                  ns%ntrunc_vr_col(c,j) = ns%ntrunc_vr_col(c,j) + cn_eca
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
                     if ( ns%decomp_npools_vr_col(c,j,l) < 0.0_r8 .and. floating_cn_ratio_decomp_pools(l) ) then
                        if ( abs(ns%decomp_npools_vr_col(c,j,l))  < ncrit ) then
                           cn_eca = cn_eca - ncrit + ns%decomp_npools_vr_col(c,j,l)
                           ns%decomp_npools_vr_col(c,j,l) = ncrit
                        else
                           write(iulog, "(A,2I8,E8.1)") 'error decomp_npools is negative: ',j,l,ns%decomp_npools_vr_col(c,j,l)
                           call endrun(msg=errMsg(__FILE__, __LINE__))
                        end if
                     end if
                     if ( ps%decomp_ppools_vr_col(c,j,l)  < 0.0_r8 .and. floating_cp_ratio_decomp_pools(l) ) then
                        if ( abs(ps%decomp_ppools_vr_col(c,j,l))  < ncrit/1e4_r8 ) then
                           cp_eca = cp_eca - ncrit/1e4_r8 + ps%decomp_ppools_vr_col(c,j,l)
                           ps%decomp_ppools_vr_col(c,j,l) = ncrit/1e4_r8
                         else 
                           write(iulog, "(A,2I8,E8.1)") 'error decomp_ppools is negative: ',j,l,ps%decomp_ppools_vr_col(c,j,l)
                           call endrun(msg=errMsg(__FILE__, __LINE__))
                         end if
                     end if

                  end do

                  ns%ntrunc_vr_col(c,j) = ns%ntrunc_vr_col(c,j) + cn_eca
                  ps%ptrunc_vr_col(c,j) = ps%ptrunc_vr_col(c,j) + cp_eca

               end do
            end do

            do fp = 1,num_soilp
               p = filter_soilp(fp)
               if (ns%retransn_patch(p) < 0._r8) then
                  write(iulog, *) 'error retransn_patch is negative: ',p
                  write(iulog, *) 'retransn_patch: ', ns%retransn_patch(p)
                  call endrun(msg=errMsg(__FILE__, __LINE__))
               end if
               if (ns%npool_patch(p) < 0._r8) then
                  write(iulog, *) 'error npool_patch is negative: ',p
                  write(iulog, *) 'npool_patch: ', ns%npool_patch(p)
                  call endrun(msg=errMsg(__FILE__, __LINE__))
               end if
               if (ps%retransp_patch(p) < 0._r8) then
                  write(iulog, *) 'error retransp_patch is negative: ',p
                  write(iulog, *) 'retransp_patch: ', ps%retransp_patch(p)
                  call endrun(msg=errMsg(__FILE__, __LINE__))
               end if
               if (ps%ppool_patch(p) < 0._r8) then
                  write(iulog, *) 'error ppool_patch is negative: ',p
                  write(iulog, *) 'ppool_patch: ', ps%ppool_patch(p)
                  call endrun(msg=errMsg(__FILE__, __LINE__))
               end if
            end do

         endif

      endif ! if (.not. is_active_betr_bgc)

    end associate

 end subroutine PrecisionControl

end module PrecisionControlMod
