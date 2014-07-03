module CNInitMod

  !-----------------------------------------------------------------------
  ! !MODULE: initCNMod
  !
  ! !DESCRIPTION:
  ! Contains cold start initial values and time constant (and flux / diagnostic vars) 
  ! for CN scheme.
  !
  ! !USES:
  use shr_kind_mod  , only: r8 => shr_kind_r8
  ! 
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: initColdCN
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine initColdCN(bounds)
    !
    ! !DESCRIPTION:
    ! Initializes time varying variables used only in coupled carbon-nitrogen mode (CN):
    !
    ! !USES:
    use clmtype       , only: pftcon, decomp_cascade_con
    use clmtype       , only: pft, pps, pcs, pcf, pns, pnf, pc13s, pc14s, pc13f, pc14f, pepv  
    use clmtype       , only: col, cps, ccs, ccf, cns, cnf, cc13s, cc14s, cc13f, cc14f, cws, cwf
    use clmtype       , only: pps_a, pcs_a, pcf_a, pc13f_a, pc14f_a, pns_a, pnf_a    
    use clmtype       , only: lun
    use clm_varpar    , only: nlevgrnd, nlevdecomp, ndecomp_pools, nlevdecomp_full, crop_prog
    use clm_varcon    , only: istsoil, zsoi, spval, denh2o
    use clm_varcon    , only: istcrop, c13ratio, c14ratio
    use clm_varctl    , only: use_c13, use_c14, use_nitrif_denitrif, use_cndv 
    use clm_varctl    , only: nsrest, nsrStartup
    use CNSetValueMod , only: cnsetpnf, cnsetpns, cnsetpcf, cnsetpcs, cnsetpps, cnsetpepv 
    use CNSetValueMod , only: cnsetcnf, cnsetcns, cnsetccf, cnsetccs, cnsetcps 
    use CNSetValueMod , only: cnzerofluxes_dwt
    use decompMod     , only: bounds_type
    use pftvarcon     , only: noveg, npcropmin
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer :: fc,fp,g,l,c,p,j,k                       ! indices
    integer :: num_special_col                         ! number of good values in special_col filter
    integer :: special_col(bounds%endc-bounds%begc+1)  ! special landunit filter - columns
    integer :: num_special_pft                         ! number of good values in special_pft filter
    integer :: special_pft(bounds%endp-bounds%begp+1)  ! special landunit filter - pfts
    integer :: num_soilcrop_col                        ! number of good values for soilcrop_col filter
    integer :: soilcrop_col(bounds%endc-bounds%begc+1) ! soil/crop landunit filter - columns
    integer :: num_soilcrop_pft                        ! number of good values for soilcrop_pft filter
    integer :: soilcrop_pft(bounds%endp-bounds%begp+1) ! soil/crop landunit filter - pfts
    integer :: num_all_col                             ! number of good values for all_col filter
    integer :: all_col(bounds%endc-bounds%begc+1)      ! all columns filter
    integer :: num_all_pft                             ! number of good values for all_pft filter
    integer :: all_pft(bounds%endp-bounds%begp+1)      ! all pfts filter
    real(r8):: vwc,psi                                 ! for calculating soilpsi
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    !  pftcon%evergreen                    Input:  [real(r8) (:)     ]  binary flag for evergreen leaf habit (0 or 1)     
    !  pftcon%woody                        Input:  [real(r8) (:)     ]  binary flag for woody lifeform (1=woody, 0=not woody)
    !  pftcon%leafcn                       Input:  [real(r8) (:)     ]  leaf C:N (gC/gN)                                  
    !  pftcon%deadwdcn                     Input:  [real(r8) (:)     ]  dead wood (xylem and heartwood) C:N (gC/gN)       

    !  cps%annsum_counter                  Output: [real(r8) (:)     ]  seconds since last annual accumulator turnover    
    !  cps%cannsum_npp                     Output: [real(r8) (:)     ]  annual sum of NPP, averaged from pft-level (gC/m2/yr)
    !  cps%cannavg_t2m                     Output: [real(r8) (:)     ]  annual average of 2m air temperature, averaged from pft-level (K)
    !  cps%wf                              Output: [real(r8) (:)     ]  soil moisture in top 0.05 m                       
    !  cps%wf2                             Output: [real(r8) (:)     ]                                                    
    !  cps%nfire                           Output: [real(r8) (:)     ]  fire counts/km2/timestep                          
    !  cps%baf_crop                        Output: [real(r8) (:)     ]  burned area fraction in crop                      
    !  cps%baf_peatf                       Output: [real(r8) (:)     ]  burned area fraction in peatland                  
    !  cps%fbac                            Output: [real(r8) (:)     ]                                                    
    !  cps%fbac1                           Output: [real(r8) (:)     ]                                                    
    !  cps%farea_burned                    Output: [real(r8) (:)     ]  timestep fractional area burned (proportion)      
    !  cps%nfixation_prof                  Output: [real(r8) (:,:)   ]  (1/m) profile for N fixation additions          
    !  cps%ndep_prof                       Output: [real(r8) (:,:)   ]  (1/m) profile for N fixation additions          
    !  cps%fpi_vr                          Output: [real(r8) (:,:)   ]                                                  
    !  cps%alt                             Output: [real(r8) (:)     ]                                                    
    !  cps%altmax                          Output: [real(r8) (:)     ]                                                    
    !  cps%altmax_lastyear                 Output: [real(r8) (:)     ]                                                    
    !  cps%som_adv_coef                    Output: [real(r8) (:,:)   ]                                                  
    !  cps%som_diffus_coef                 Output: [real(r8) (:,:)   ]                                                  
    !  cps%alt_indx                        Output: [integer (:)      ]                                                     
    !  cps%altmax_indx                     Output: [integer (:)      ]                                                     
    !  cps%altmax_lastyear_indx            Output: [integer (:)      ]                                                     

    !  cwf%qflx_drain                      Output: [real(r8) (:)     ]  sub-surface runoff (mm H2O /s)                    
    !  cwf%qflx_surf                       Output: [real(r8) (:)     ]  surface runoff (mm H2O /s)                        

    !  ccs%seedc                           Output: [real(r8) (:)     ]  (gC/m2) column-level pool for seeding new PFTs    
    !  ccs%prod10c                         Output: [real(r8) (:)     ]  (gC/m2) wood product C pool, 10-year lifespan     
    !  ccs%prod100c                        Output: [real(r8) (:)     ]  (gC/m2) wood product C pool, 100-year lifespan    
    !  ccs%totprodc                        Output: [real(r8) (:)     ]  (gC/m2) total wood product C                      
    !  ccs%totcolc                         Output: [real(r8) (:)     ]  (gC/m2) total column carbon, incl veg and cpool   
    !  ccs%col_ctrunc                      Output: [real(r8) (:)     ]  (gC/m2) column-level sink for C truncation (diagnostic)
    !  ccs%decomp_cpools                   Output: [real(r8) (:,:)   ]  (gC/m2)  decomposing (litter, cwd, soil) c pools
    !  ccs%decomp_cpools_1m                Output: [real(r8) (:,:)   ]  (gC/m2)  Diagnostic: decomposing (litter, cwd, soil) c pools to 1 meter
    !  ccs%decomp_cpools_vr                Output: [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
    !  ccs%cwdc                            Output: [real(r8) (:)     ]  (gC/m2) coarse woody debris C                     
    !  ccs%totecosysc                      Output: [real(r8) (:)     ]  (gC/m2) total ecosystem carbon, incl veg but excl cpool
    !  ccs%totlitc                         Output: [real(r8) (:)     ]  (gC/m2) total litter carbon                       
    !  ccs%totsomc                         Output: [real(r8) (:)     ]  (gC/m2) total soil organic matter carbon          
    !  ccs%totlitc_1m                      Output: [real(r8) (:)     ]  (gC/m2) total litter carbon to 1 meter            
    !  ccs%totsomc_1m                      Output: [real(r8) (:)     ]  (gC/m2) total soil organic matter carbon to 1 meter
    !  ccs%col_ctrunc_vr                   Output: [real(r8) (:,:)   ]  (gC/m3) column-level sink for C truncation (prognostic)

    !  cns%seedn                           Output: [real(r8) (:)     ]  (gN/m2) column-level pool for seeding new PFTs    
    !  cns%prod10n                         Output: [real(r8) (:)     ]  (gN/m2) wood product N pool, 10-year lifespan     
    !  cns%prod100n                        Output: [real(r8) (:)     ]  (gN/m2) wood product N pool, 100-year lifespan    
    !  cns%totprodn                        Output: [real(r8) (:)     ]  (gN/m2) total wood product N                      
    !  cns%totcoln                         Output: [real(r8) (:)     ]  (gN/m2) total column nitrogen, incl veg           
    !  cns%decomp_npools                   Output: [real(r8) (:,:)   ]  (gC/m2)  decomposing (litter, cwd, soil) N pools
    !  cns%decomp_npools_vr                Output: [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
    !  cns%decomp_npools_1m                Output: [real(r8) (:,:)   ]  (gN/m2)  diagnostic: decomposing (litter, cwd, soil) N pools to 1 meter
    !  cns%cwdn                            Output: [real(r8) (:)     ]  (gN/m2) coarse woody debris N                     
    !  cns%totecosysn                      Output: [real(r8) (:)     ]  (gN/m2) total ecosystem nitrogen, incl veg        
    !  cns%totlitn                         Output: [real(r8) (:)     ]  (gN/m2) total litter nitrogen                     
    !  cns%totsomn                         Output: [real(r8) (:)     ]  (gN/m2) total soil organic matter nitrogen        
    !  cns%totlitn_1m                      Output: [real(r8) (:)     ]  (gN/m2) total litter nitrogen to 1 meter          
    !  cns%totsomn_1m                      Output: [real(r8) (:)     ]  (gN/m2) total soil organic matter nitrogen to 1 meter
    !  cns%smin_nh4_vr                     Output: [real(r8) (:,:)   ]  (gN/m3) soil mineral NH4 pool                   
    !  cns%smin_no3_vr                     Output: [real(r8) (:,:)   ]  (gN/m3) soil mineral NO3 pool                   
    !  cns%smin_nh4                        Output: [real(r8) (:)     ]  (gN/m2) soil mineral NH4 pool                     
    !  cns%smin_no3                        Output: [real(r8) (:)     ]  (gN/m2) soil mineral NO3 pool                     
    !  cns%col_ntrunc_vr                   Output: [real(r8) (:,:)   ]  (gN/m3) column-level sink for N truncation      
    !  cns%sminn                           Output: [real(r8) (:)     ]  (gN/m2) soil mineral N                            
    !  cns%sminn_vr                        Output: [real(r8) (:,:)   ]  (gN/m3) soil mineral N                          

    !  cc13s%seedc                         Output: [real(r8) (:)     ]  (gC/m2) column-level pool for seeding new PFTs    
    !  cc13s%prod10c                       Output: [real(r8) (:)     ]  (gC/m2) wood product C13 pool, 10-year lifespan   
    !  cc13s%prod100c                      Output: [real(r8) (:)     ]  (gC/m2) wood product C13 pool, 100-year lifespan  
    !  cc13s%totprodc                      Output: [real(r8) (:)     ]  (gC/m2) total wood product C13                    
    !  cc13s%cwdc                          Output: [real(r8) (:)     ]  (gC/m2) coarse woody debris C                     
    !  cc13s%decomp_cpools                 Output: [real(r8) (:,:)   ]  (gC/m2)  decomposing (litter, cwd, soil) c pools
    !  cc13s%decomp_cpools_vr              Output: [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
    !  cc13s%col_ctrunc_vr                 Output: [real(r8) (:,:)   ]  (gC/m3) C truncation term                       
    !  cc13s%decomp_cpools_1m              Output: [real(r8) (:,:)   ]  (gC/m2)  Diagnostic: decomposing (litter, cwd, soil) c pools to 1 meter

    !  pps%forc_hgt_u_pft                  Output: [real(r8) (:)     ]  observational height of wind at pft-level [m]      
    !  pps%alphapsnsun                     Output: [real(r8) (:)     ] sunlit 13c fractionation ([])                      
    !  pps%alphapsnsha                     Output: [real(r8) (:)     ] shaded 13c fractionation ([])                      
    !  pps%laisha                          Output: [real(r8) (:)     ]  shaded projected leaf area index                  
    !  pps%laisun                          Output: [real(r8) (:)     ]  sunlit projected leaf area index                  

    !  pc13s%leafc                         Output: [real(r8) (:)     ]  (gC/m2) leaf C                                    
    !  pc13s%leafc_storage                 Output: [real(r8) (:)     ]  (gC/m2) leaf C storage                            
    !  pc13s%leafc_xfer                    Output: [real(r8) (:)     ]  (gC/m2) leaf C transfer                           
    !  pc13s%frootc                        Output: [real(r8) (:)     ]  (gC/m2) fine root C                               
    !  pc13s%frootc_storage                Output: [real(r8) (:)     ]  (gC/m2) fine root C storage                       
    !  pc13s%frootc_xfer                   Output: [real(r8) (:)     ]  (gC/m2) fine root C transfer                      
    !  pc13s%livestemc                     Output: [real(r8) (:)     ]  (gC/m2) live stem C                               
    !  pc13s%livestemc_storage             Output: [real(r8) (:)     ]  (gC/m2) live stem C storage                       
    !  pc13s%livestemc_xfer                Output: [real(r8) (:)     ]  (gC/m2) live stem C transfer                      
    !  pc13s%deadstemc                     Output: [real(r8) (:)     ]  (gC/m2) dead stem C                               
    !  pc13s%deadstemc_storage             Output: [real(r8) (:)     ]  (gC/m2) dead stem C storage                       
    !  pc13s%deadstemc_xfer                Output: [real(r8) (:)     ]  (gC/m2) dead stem C transfer                      
    !  pc13s%livecrootc                    Output: [real(r8) (:)     ]  (gC/m2) live coarse root C                        
    !  pc13s%livecrootc_storage            Output: [real(r8) (:)     ]  (gC/m2) live coarse root C storage                
    !  pc13s%livecrootc_xfer               Output: [real(r8) (:)     ]  (gC/m2) live coarse root C transfer               
    !  pc13s%deadcrootc                    Output: [real(r8) (:)     ]  (gC/m2) dead coarse root C                        
    !  pc13s%deadcrootc_storage            Output: [real(r8) (:)     ]  (gC/m2) dead coarse root C storage                
    !  pc13s%deadcrootc_xfer               Output: [real(r8) (:)     ]  (gC/m2) dead coarse root C transfer               
    !  pc13s%gresp_storage                 Output: [real(r8) (:)     ]  (gC/m2) growth respiration storage                
    !  pc13s%gresp_xfer                    Output: [real(r8) (:)     ]  (gC/m2) growth respiration transfer               
    !  pc13s%cpool                         Output: [real(r8) (:)     ]  (gC/m2) temporary photosynthate C pool            
    !  pc13s%xsmrpool                      Output: [real(r8) (:)     ]  (gC/m2) temporary photosynthate C pool            
    !  pc13s%pft_ctrunc                    Output: [real(r8) (:)     ]  (gC/m2) C truncation term                         
    !  pc13s%totvegc                       Output: [real(r8) (:)     ]  (gC/m2) total vegetation carbon, excluding cpool  
    !  pc13f%psnsun                        Output: [real(r8) (:)     ]  sunlit leaf photosynthesis (umol CO2 /m**2/ s)    
    !  pc13f%psnsha                        Output: [real(r8) (:)     ]  shaded leaf photosynthesis (umol CO2 /m**2/ s)    

    !  cc14s%seedc                         Output: [real(r8) (:)     ]  (gC/m2) column-level pool for seeding new PFTs    
    !  cc14s%prod10c                       Output: [real(r8) (:)     ]  (gC/m2) wood product C14 pool, 10-year lifespan   
    !  cc14s%prod100c                      Output: [real(r8) (:)     ]  (gC/m2) wood product C14 pool, 100-year lifespan  
    !  cc14s%totprodc                      Output: [real(r8) (:)     ]  (gC/m2) total wood product C14                    
    !  cc14s%cwdc                          Output: [real(r8) (:)     ]  (gC/m2) coarse woody debris C                     
    !  cc14s%decomp_cpools                 Output: [real(r8) (:,:)   ]  (gC/m2)  decomposing (litter, cwd, soil) c pools
    !  cc14s%decomp_cpools_vr              Output: [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
    !  cc14s%col_ctrunc_vr                 Output: [real(r8) (:,:)   ]  (gC/m3) C truncation term                       
    !  cc14s%decomp_cpools_1m              Output: [real(r8) (:,:)   ]  (gC/m2)  Diagnostic: decomposing (litter, cwd, soil) c pools to 1 meter

    !  pc14s%leafc                         Output: [real(r8) (:)     ]  (gC/m2) leaf C                                    
    !  pc14s%leafc_storage                 Output: [real(r8) (:)     ]  (gC/m2) leaf C storage                            
    !  pc14s%leafc_xfer                    Output: [real(r8) (:)     ]  (gC/m2) leaf C transfer                           
    !  pc14s%frootc                        Output: [real(r8) (:)     ]  (gC/m2) fine root C                               
    !  pc14s%frootc_storage                Output: [real(r8) (:)     ]  (gC/m2) fine root C storage                       
    !  pc14s%frootc_xfer                   Output: [real(r8) (:)     ]  (gC/m2) fine root C transfer                      
    !  pc14s%livestemc                     Output: [real(r8) (:)     ]  (gC/m2) live stem C                               
    !  pc14s%livestemc_storage             Output: [real(r8) (:)     ]  (gC/m2) live stem C storage                       
    !  pc14s%livestemc_xfer                Output: [real(r8) (:)     ]  (gC/m2) live stem C transfer                      
    !  pc14s%deadstemc                     Output: [real(r8) (:)     ]  (gC/m2) dead stem C                               
    !  pc14s%deadstemc_storage             Output: [real(r8) (:)     ]  (gC/m2) dead stem C storage                       
    !  pc14s%deadstemc_xfer                Output: [real(r8) (:)     ]  (gC/m2) dead stem C transfer                      
    !  pc14s%livecrootc                    Output: [real(r8) (:)     ]  (gC/m2) live coarse root C                        
    !  pc14s%livecrootc_storage            Output: [real(r8) (:)     ]  (gC/m2) live coarse root C storage                
    !  pc14s%livecrootc_xfer               Output: [real(r8) (:)     ]  (gC/m2) live coarse root C transfer               
    !  pc14s%deadcrootc                    Output: [real(r8) (:)     ]  (gC/m2) dead coarse root C                        
    !  pc14s%deadcrootc_storage            Output: [real(r8) (:)     ]  (gC/m2) dead coarse root C storage                
    !  pc14s%deadcrootc_xfer               Output: [real(r8) (:)     ]  (gC/m2) dead coarse root C transfer               
    !  pc14s%gresp_storage                 Output: [real(r8) (:)     ]  (gC/m2) growth respiration storage                
    !  pc14s%gresp_xfer                    Output: [real(r8) (:)     ]  (gC/m2) growth respiration transfer               
    !  pc14s%cpool                         Output: [real(r8) (:)     ]  (gC/m2) temporary photosynthate C pool            
    !  pc14s%xsmrpool                      Output: [real(r8) (:)     ]  (gC/m2) temporary photosynthate C pool            
    !  pc14s%pft_ctrunc                    Output: [real(r8) (:)     ]  (gC/m2) C truncation term                         
    !  pc14s%totvegc                       Output: [real(r8) (:)     ]  (gC/m2) total vegetation carbon, excluding cpool  
    !  pc14f%psnsun                        Output: [real(r8) (:)     ]  sunlit leaf photosynthesis (umol CO2 /m**2/ s)    
    !  pc14f%psnsha                        Output: [real(r8) (:)     ]  shaded leaf photosynthesis (umol CO2 /m**2/ s)    

    !  pcs%leafc                           Output: [real(r8) (:)     ]  (gC/m2) leaf C                                    
    !  pcs%leafc_storage                   Output: [real(r8) (:)     ]  (gC/m2) leaf C storage                            
    !  pcs%leafc_xfer                      Output: [real(r8) (:)     ]  (gC/m2) leaf C transfer                           
    !  pcs%grainc                          Output: [real(r8) (:)     ]  (gC/m2) grain C                                   
    !  pcs%grainc_storage                  Output: [real(r8) (:)     ]  (gC/m2) grain C storage                           
    !  pcs%grainc_xfer                     Output: [real(r8) (:)     ]  (gC/m2) grain C transfer                          
    !  pcs%frootc                          Output: [real(r8) (:)     ]  (gC/m2) fine root C                               
    !  pcs%frootc_storage                  Output: [real(r8) (:)     ]  (gC/m2) fine root C storage                       
    !  pcs%frootc_xfer                     Output: [real(r8) (:)     ]  (gC/m2) fine root C transfer                      
    !  pcs%livestemc                       Output: [real(r8) (:)     ]  (gC/m2) live stem C                               
    !  pcs%livestemc_storage               Output: [real(r8) (:)     ]  (gC/m2) live stem C storage                       
    !  pcs%livestemc_xfer                  Output: [real(r8) (:)     ]  (gC/m2) live stem C transfer                      
    !  pcs%deadstemc                       Output: [real(r8) (:)     ]  (gC/m2) dead stem C                               
    !  pcs%deadstemc_storage               Output: [real(r8) (:)     ]  (gC/m2) dead stem C storage                       
    !  pcs%deadstemc_xfer                  Output: [real(r8) (:)     ]  (gC/m2) dead stem C transfer                      
    !  pcs%livecrootc                      Output: [real(r8) (:)     ]  (gC/m2) live coarse root C                        
    !  pcs%livecrootc_storage              Output: [real(r8) (:)     ]  (gC/m2) live coarse root C storage                
    !  pcs%livecrootc_xfer                 Output: [real(r8) (:)     ]  (gC/m2) live coarse root C transfer               
    !  pcs%deadcrootc                      Output: [real(r8) (:)     ]  (gC/m2) dead coarse root C                        
    !  pcs%deadcrootc_storage              Output: [real(r8) (:)     ]  (gC/m2) dead coarse root C storage                
    !  pcs%deadcrootc_xfer                 Output: [real(r8) (:)     ]  (gC/m2) dead coarse root C transfer               
    !  pcs%gresp_storage                   Output: [real(r8) (:)     ]  (gC/m2) growth respiration storage                
    !  pcs%gresp_xfer                      Output: [real(r8) (:)     ]  (gC/m2) growth respiration transfer               
    !  pcs%cpool                           Output: [real(r8) (:)     ]  (gC/m2) temporary photosynthate C pool            
    !  pcs%xsmrpool                        Output: [real(r8) (:)     ]  (gC/m2) abstract C pool to meet excess MR demand  
    !  pcs%woodc                           Output: [real(r8) (:)     ]  (gC/m2) pft-level wood C                          
    !  pcs%dispvegc                        Output: [real(r8) (:)     ]  (gC/m2) displayed veg carbon, excluding storage and cpool
    !  pcs%pft_ctrunc                      Output: [real(r8) (:)     ]  (gC/m2) pft-level sink for C truncation           
    !  pcs%storvegc                        Output: [real(r8) (:)     ]  (gC/m2) stored vegetation carbon, excluding cpool 
    !  pcs%totpftc                         Output: [real(r8) (:)     ]  (gC/m2) total pft-level carbon, including cpool   
    !  pcs%totvegc                         Output: [real(r8) (:)     ]  (gC/m2) total vegetation carbon, excluding cpool  
    !  pcf%psnsun                          Output: [real(r8) (:)     ]  sunlit leaf photosynthesis (umol CO2 /m**2/ s)    
    !  pcf%psnsha                          Output: [real(r8) (:)     ]  shaded leaf photosynthesis (umol CO2 /m**2/ s)    

    !  pns%leafn                           Output: [real(r8) (:)     ]  (gN/m2) leaf N                                    
    !  pns%leafn_storage                   Output: [real(r8) (:)     ]  (gN/m2) leaf N storage                            
    !  pns%leafn_xfer                      Output: [real(r8) (:)     ]  (gN/m2) leaf N transfer                           
    !  pns%grainn                          Output: [real(r8) (:)     ]  (gN/m2) grain N                                   
    !  pns%grainn_storage                  Output: [real(r8) (:)     ]  (gN/m2) grain N storage                           
    !  pns%grainn_xfer                     Output: [real(r8) (:)     ]  (gN/m2) grain N transfer                          
    !  pns%frootn                          Output: [real(r8) (:)     ]  (gN/m2) fine root N                               
    !  pns%frootn_storage                  Output: [real(r8) (:)     ]  (gN/m2) fine root N storage                       
    !  pns%frootn_xfer                     Output: [real(r8) (:)     ]  (gN/m2) fine root N transfer                      
    !  pns%livestemn                       Output: [real(r8) (:)     ]  (gN/m2) live stem N                               
    !  pns%livestemn_storage               Output: [real(r8) (:)     ]  (gN/m2) live stem N storage                       
    !  pns%livestemn_xfer                  Output: [real(r8) (:)     ]  (gN/m2) live stem N transfer                      
    !  pns%deadstemn                       Output: [real(r8) (:)     ]  (gN/m2) dead stem N                               
    !  pns%deadstemn_storage               Output: [real(r8) (:)     ]  (gN/m2) dead stem N storage                       
    !  pns%deadstemn_xfer                  Output: [real(r8) (:)     ]  (gN/m2) dead stem N transfer                      
    !  pns%livecrootn                      Output: [real(r8) (:)     ]  (gN/m2) live coarse root N                        
    !  pns%livecrootn_storage              Output: [real(r8) (:)     ]  (gN/m2) live coarse root N storage                
    !  pns%livecrootn_xfer                 Output: [real(r8) (:)     ]  (gN/m2) live coarse root N transfer               
    !  pns%deadcrootn                      Output: [real(r8) (:)     ]  (gN/m2) dead coarse root N                        
    !  pns%deadcrootn_storage              Output: [real(r8) (:)     ]  (gN/m2) dead coarse root N storage                
    !  pns%deadcrootn_xfer                 Output: [real(r8) (:)     ]  (gN/m2) dead coarse root N transfer               
    !  pns%retransn                        Output: [real(r8) (:)     ]  (gN/m2) plant pool of retranslocated N            
    !  pns%npool                           Output: [real(r8) (:)     ]  (gN/m2) temporary plant N pool                    
    !  pns%dispvegn                        Output: [real(r8) (:)     ]  (gN/m2) displayed veg nitrogen, excluding storage 
    !  pns%pft_ntrunc                      Output: [real(r8) (:)     ]  (gN/m2) pft-level sink for N truncation           
    !  pns%storvegn                        Output: [real(r8) (:)     ]  (gN/m2) stored vegetation nitrogen                
    !  pns%totpftn                         Output: [real(r8) (:)     ]  (gN/m2) total pft-level nitrogen                  
    !  pns%totvegn                         Output: [real(r8) (:)     ]  (gN/m2) total vegetation nitrogen                 
    !  pnf%soyfixn                         Output: [real(r8) (:)     ]                                                    
    !  pnf%fert                            Output: [real(r8) (:)     ]                                                    

    !  pepv%xsmrpool_c13ratio              Output: [real(r8) (:)     ]  C flux assigned to recovery of negative cpool (gC/m2/s)
    !  pepv%rc13_canair                    Output: [real(r8) (:)     ]  C13O2/C12O2 in canopy air                          
    !  pepv%rc13_psnsun                    Output: [real(r8) (:)     ]  C13O2/C12O2 in sunlit canopy psn flux              
    !  pepv%rc13_psnsha                    Output: [real(r8) (:)     ]  C13O2/C12O2 in shaded canopy psn flux              
    !  pepv%rc14_atm                       Output: [real(r8) (:)     ] C14O2/C12O2 in atmosphere                          
    !  pepv%fert_counter                   Output: [real(r8) (:)     ]                                                    
    !  pepv%grain_flag                     Output: [real(r8) (:)     ]                                                    
    !  pepv%dormant_flag                   Output: [real(r8) (:)     ]  dormancy flag                                     
    !  pepv%days_active                    Output: [real(r8) (:)     ]  number of days since last dormancy                
    !  pepv%onset_flag                     Output: [real(r8) (:)     ]  onset flag                                        
    !  pepv%onset_counter                  Output: [real(r8) (:)     ]  onset days counter                                
    !  pepv%onset_gdf90dflag               Output: [real(r8) (:)     ]  onset flag for growing degree day sum             
    !  pepv%onset_fdd                      Output: [real(r8) (:)     ]  onset freezing degree days counter                
    !  pepv%onset_gdd                      Output: [real(r8) (:)     ]  onset growing degree days                         
    !  pepv%onset_swi                      Output: [real(r8) (:)     ]  onset soil water index                            
    !  pepv%offset_flag                    Output: [real(r8) (:)     ]  offset flag                                       
    !  pepv%offset_counter                 Output: [real(r8) (:)     ]  offset days counter                               
    !  pepv%offset_fdd                     Output: [real(r8) (:)     ]  offset freezing degree days counter               
    !  pepv%offset_swi                     Output: [real(r8) (:)     ]  offset soil water index                           
    !  pepv%lgsf                           Output: [real(r8) (:)     ]  long growing season factor [0-1]                  
    !  pepv%bglfr                          Output: [real(r8) (:)     ]  background litterfall rate (1/s)                  
    !  pepv%bgtr                           Output: [real(r8) (:)     ]  background transfer rate (1/s)                    
    !  pepv%annavg_t2m                     Output: [real(r8) (:)     ]  annual average 2m air temperature (K)             
    !  pepv%tempavg_t2m                    Output: [real(r8) (:)     ]  temporary average 2m air temperature (K)          
    !  pepv%gpp                            Output: [real(r8) (:)     ]  GPP flux before downregulation (gC/m2/s)          
    !  pepv%availc                         Output: [real(r8) (:)     ]  C flux available for allocation (gC/m2/s)         
    !  pepv%xsmrpool_recover               Output: [real(r8) (:)     ]  C flux assigned to recovery of negative cpool (gC/m2/s)
    !  pepv%alloc_pnow                     Output: [real(r8) (:)     ]  fraction of current allocation to display as new growth (DIM)
    !  pepv%c_allometry                    Output: [real(r8) (:)     ]  C allocation index (DIM)                          
    !  pepv%n_allometry                    Output: [real(r8) (:)     ]  N allocation index (DIM)                          
    !  pepv%plant_ndemand                  Output: [real(r8) (:)     ]  N flux required to support initial GPP (gN/m2/s)  
    !  pepv%tempsum_potential_gpp          Output: [real(r8) (:)     ]  temporary annual sum of plant_ndemand             
    !  pepv%annsum_potential_gpp           Output: [real(r8) (:)     ]  annual sum of plant_ndemand                       
    !  pepv%tempmax_retransn               Output: [real(r8) (:)     ]  temporary max of retranslocated N pool (gN/m2)    
    !  pepv%annmax_retransn                Output: [real(r8) (:)     ]  annual max of retranslocated N pool (gN/m2)       
    !  pepv%avail_retransn                 Output: [real(r8) (:)     ]  N flux available from retranslocation pool (gN/m2/s)
    !  pepv%plant_nalloc                   Output: [real(r8) (:)     ]  total allocated N flux (gN/m2/s)                  
    !  pepv%plant_calloc                   Output: [real(r8) (:)     ]  total allocated C flux (gC/m2/s)                  
    !  pepv%excess_cflux                   Output: [real(r8) (:)     ]  C flux not allocated due to downregulation (gC/m2/s)
    !  pepv%downreg                        Output: [real(r8) (:)     ]  fractional reduction in GPP due to N limitation (DIM)
    !  pepv%tempsum_npp                    Output: [real(r8) (:)     ]  temporary annual sum of NPP                       
    !  pepv%annsum_npp                     Output: [real(r8) (:)     ]  annual sum of NPP                                 
    !  pepv%tempsum_litfall                Output: [real(r8) (:)     ]  temporary annual sum of litfall                   
    !  pepv%annsum_litfall                 Output: [real(r8) (:)     ]  annual sum of litfall                             
    !  pepv%prev_frootc_to_litter          Output: [real(r8) (:)     ]  previous timestep froot C litterfall flux (gC/m2/s)
    !  pepv%prev_leafc_to_litter           Output: [real(r8) (:)     ]  previous timestep leaf C litterfall flux (gC/m2/s) 

    !  decomp_cascade_con%initial_cn_ratio Output: [real(r8) (:)     ]  c:n ratio for initialization of pools             
    !  decomp_cascade_con%initial_stock    Output: [real(r8) (:)     ]  initial concentration for seeding at spinup       
    !-----------------------------------------------------------------------

    ! Set column filters

    num_soilcrop_col = 0
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          num_soilcrop_col = num_soilcrop_col + 1
          soilcrop_col(num_soilcrop_col) = c
       end if
    end do

    num_special_col = 0
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%ifspecial(l)) then
          num_special_col = num_special_col + 1
          special_col(num_special_col) = c
       end if
    end do

    num_all_col = 0
    do c = bounds%begc, bounds%endc
       num_all_col = num_all_col + 1
       all_col(num_all_col) = c
    end do

    ! Set pft filters

    num_soilcrop_pft = 0
    do p = bounds%begp,bounds%endp
       l = pft%landunit(p)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          num_soilcrop_pft = num_soilcrop_pft + 1
          soilcrop_pft(num_soilcrop_pft) = p
       end if
    end do

    num_special_pft = 0
    do p = bounds%begp,bounds%endp
       l = pft%landunit(p)
       if (lun%ifspecial(l)) then
          num_special_pft = num_special_pft + 1
          special_pft(num_special_pft) = p
       end if
    end do

    num_all_pft = 0
    do p = bounds%begp, bounds%endp
       num_all_pft = num_all_pft + 1
       all_pft(num_all_pft) = p
    end do

    ! Added 5/4/04, PET: initialize forc_hgt_u (gridcell-level),
    ! since this is not initialized before first call to CNVegStructUpdate,
    ! and it is required to set the upper bound for canopy top height.
    ! Changed 3/21/08, KO: still needed but don't have sufficient information 
    ! to set this properly (e.g., pft-level displacement height and roughness 
    ! length). So leave at 30m.
    do p = bounds%begp,bounds%endp
       pps%forc_hgt_u_pft(p) = 30._r8
    end do

    ! initialize column-level variables
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then

          ! column physical state variables
          cps%annsum_counter(c) = 0._r8   !initialized to spval for all columns
          cps%cannsum_npp(c)    = 0._r8   !initialized to spval for all columns
          cps%cannavg_t2m(c)    = 280._r8 !initialized to nanr  for all columns

          ! fire related variables 
          ! cps%wf needs to be non zero so the first time step has no fires
          cps%wf(c)           = 1.0_r8   !initialized to nanr  for all columns
          cps%wf2(c)          = 1.0_r8   !initialized to nanr  for all columns
          cps%baf_crop(c)     = 0._r8    !initialized to nanr  for all columns
          cps%baf_peatf(c)    = 0._r8    !initialized to nanr  for all columns
          cps%fbac(c)         = 0._r8    !initialized to nanr  for all columns
          cps%fbac1(c)        = 0._r8    !initialized to nanr  for all columns
          cps%farea_burned(c) = 0._r8    !initialized to nanr  for all columns

          ! TODO: the following if-clause is there only for backwards compatibility
          ! without it answers are not the same between clm4_5_66 and clm4_5_65
          ! not clear which one is correct - note all of the other fire-related variables
          ! above are initialized to nanr in clmtineInitMod.F90
          if (nsrest == nsrStartup) then 
             cps%nfire(c) = 0._r8
          end if

          ! initialize fpi_vr so that levels below nlevsoi are not nans
          cps%fpi_vr(c,1:nlevdecomp_full)          = 0._r8 !initialized to nanr  for all columns
          cps%som_adv_coef(c,1:nlevdecomp_full)    = 0._r8 !initialized to spval for all columns
          cps%som_diffus_coef(c,1:nlevdecomp_full) = 0._r8 !initialized to spval for all columns

          ! initialize the profiles for converting to vertically resolved carbon pools
          cps%nfixation_prof(c,1:nlevdecomp_full)  = 0._r8 !iniitialized to spval for all columns
          cps%ndep_prof(c,1:nlevdecomp_full)       = 0._r8 !iniitialized to spval for all columns

          ! and define alt variables to be zero
          cps%alt(c)               = 0._r8 !iniitialized to spval for all columns
          cps%altmax(c)            = 0._r8 !iniitialized to spval for all columns
          cps%altmax_lastyear(c)   = 0._r8 !iniitialized to spval for all columns
          cps%alt_indx(c)          = 0     !initiialized to huge  for all columns
          cps%altmax_indx(c)       = 0     !initiialized to huge  for all columns
          cps%altmax_lastyear_indx = 0     !initiialized to huge  for all columns

          ! needed for CNNLeaching
          cwf%qflx_drain(c) = 0._r8
          cwf%qflx_surf(c)  = 0._r8

          ! column carbon state variable initialization
          do j = 1, nlevdecomp
             do k = 1, ndecomp_pools
                if (zsoi(j) .lt. 0.3 ) then  !! only initialize upper soil column
                   ccs%decomp_cpools_vr(c,j,k) = decomp_cascade_con%initial_stock(k)
                else
                   ccs%decomp_cpools_vr(c,j,k) = 0._r8
                endif
             end do
             ccs%col_ctrunc_vr(c,j) = 0._r8
          end do
          if ( nlevdecomp .gt. 1 ) then
             do j = nlevdecomp+1, nlevdecomp_full
                do k = 1, ndecomp_pools
                   ccs%decomp_cpools_vr(c,j,k) = 0._r8
                end do
                ccs%col_ctrunc_vr(c,j) = 0._r8
             end do
          end if
          ccs%decomp_cpools(c,1:ndecomp_pools)    = decomp_cascade_con%initial_stock(1:ndecomp_pools)
          ccs%decomp_cpools_1m(c,1:ndecomp_pools) = decomp_cascade_con%initial_stock(1:ndecomp_pools)

          ccs%cwdc(c)       = 0._r8
          ccs%col_ctrunc(c) = 0._r8
          ccs%totlitc(c)    = 0._r8
          ccs%totsomc(c)    = 0._r8
          ccs%totlitc_1m(c) = 0._r8
          ccs%totsomc_1m(c) = 0._r8
          ccs%totecosysc(c) = 0._r8
          ccs%totcolc(c)    = 0._r8

          if ( use_c13 ) then
             do j = 1, nlevdecomp
                do k = 1, ndecomp_pools
                   cc13s%decomp_cpools_vr(c,j,k) = ccs%decomp_cpools_vr(c,j,k) * c13ratio
                end do
                cc13s%col_ctrunc_vr(c,j) = ccs%col_ctrunc_vr(c,j) * c13ratio
             end do
             if ( nlevdecomp .gt. 1 ) then
                do j = nlevdecomp+1, nlevdecomp_full
                   do k = 1, ndecomp_pools
                      cc13s%decomp_cpools_vr(c,j,k) = 0._r8
                   end do
                   cc13s%col_ctrunc_vr(c,j) = 0._r8
                end do
             end if
             cc13s%cwdc(c) = ccs%cwdc(c) * c13ratio
             do k = 1, ndecomp_pools
                cc13s%decomp_cpools(c,k)    = ccs%decomp_cpools(c,k)    * c13ratio
                cc13s%decomp_cpools_1m(c,k) = ccs%decomp_cpools_1m(c,k) * c13ratio
             end do
          endif

          if ( use_c14 ) then
             do j = 1, nlevdecomp
                do k = 1, ndecomp_pools
                   cc14s%decomp_cpools_vr(c,j,k) = ccs%decomp_cpools_vr(c,j,k) * c14ratio
                end do
                cc14s%col_ctrunc_vr(c,j) = ccs%col_ctrunc_vr(c,j) * c14ratio
             end do
             if ( nlevdecomp .gt. 1 ) then
                do j = nlevdecomp+1, nlevdecomp_full
                   do k = 1, ndecomp_pools
                      cc14s%decomp_cpools_vr(c,j,k) = 0._r8
                   end do
                   cc14s%col_ctrunc_vr(c,j) = 0._r8
                end do
             end if
             cc14s%cwdc(c) = ccs%cwdc(c) * c14ratio
             do k = 1, ndecomp_pools
                cc14s%decomp_cpools(c,k)    = ccs%decomp_cpools(c,k)    * c14ratio
                cc14s%decomp_cpools_1m(c,k) = ccs%decomp_cpools_1m(c,k) * c14ratio
             end do
          endif

          ! column nitrogen state variables
          cns%sminn(c) = 0._r8
          do j = 1, nlevdecomp
             do k = 1, ndecomp_pools
                cns%decomp_npools_vr(c,j,k) = ccs%decomp_cpools_vr(c,j,k) / decomp_cascade_con%initial_cn_ratio(k)
             end do
             cns%sminn_vr(c,j) = 0._r8
             cns%col_ntrunc_vr(c,j) = 0._r8
          end do
          if ( nlevdecomp .gt. 1 ) then
             do j = nlevdecomp+1, nlevdecomp_full
                do k = 1, ndecomp_pools
                   cns%decomp_npools_vr(c,j,k) = 0._r8
                end do
                cns%sminn_vr(c,j) = 0._r8
                cns%col_ntrunc_vr(c,j) = 0._r8
             end do
          end if
          do k = 1, ndecomp_pools
             cns%decomp_npools(c,k)    = ccs%decomp_cpools(c,k)    / decomp_cascade_con%initial_cn_ratio(k)
             cns%decomp_npools_1m(c,k) = ccs%decomp_cpools_1m(c,k) / decomp_cascade_con%initial_cn_ratio(k)
          end do

          if (use_nitrif_denitrif) then
             do j = 1, nlevdecomp_full
                cns%smin_nh4_vr(c,j) = 0._r8
                cns%smin_no3_vr(c,j) = 0._r8
             end do
             cns%smin_nh4(c) = 0._r8
             cns%smin_no3(c) = 0._r8
          end if
          cns%totlitn(c)    = 0._r8
          cns%totsomn(c)    = 0._r8
          cns%totlitn_1m(c) = 0._r8
          cns%totsomn_1m(c) = 0._r8
          cns%totecosysn(c) = 0._r8
          cns%totcoln(c)    = 0._r8
          cns%cwdn(c)       = 0._r8

          ! dynamic landcover state variables
          ccs%seedc(c)         = 0._r8
          ccs%prod10c(c)       = 0._r8
          ccs%prod100c(c)      = 0._r8
          ccs%totprodc(c)      = 0._r8
          if ( use_c13 ) then
             cc13s%seedc(c)    = 0._r8
             cc13s%prod10c(c)  = 0._r8
             cc13s%prod100c(c) = 0._r8
             cc13s%totprodc(c) = 0._r8
          endif
          if ( use_c14 ) then
             cc14s%seedc(c)    = 0._r8
             cc14s%prod10c(c)   = 0._r8
             cc14s%prod100c(c)  = 0._r8
             cc14s%totprodc(c)  = 0._r8
          endif
          cns%seedn(c)         = 0._r8
          cns%prod10n(c)       = 0._r8
          cns%prod100n(c)      = 0._r8
          cns%totprodn(c)      = 0._r8

          ! also initialize dynamic landcover fluxes so that they have
          ! real values on first timestep, prior to calling pftdyn_cnbal
          ccf%dwt_seedc_to_leaf(c)     = 0._r8
          ccf%dwt_seedc_to_deadstem(c) = 0._r8
          ccf%dwt_conv_cflux(c)        = 0._r8
          ccf%lf_conv_cflux(c)         = 0._r8
          ccf%dwt_prod10c_gain(c)      = 0._r8
          ccf%prod10c_loss(c)          = 0._r8
          ccf%dwt_prod100c_gain(c)     = 0._r8
          ccf%prod100c_loss(c)         = 0._r8
          do j = 1, nlevdecomp_full
             ccf%dwt_frootc_to_litr_met_c(c,j) = 0._r8
             ccf%dwt_frootc_to_litr_cel_c(c,j) = 0._r8
             ccf%dwt_frootc_to_litr_lig_c(c,j) = 0._r8
             ccf%dwt_livecrootc_to_cwdc(c,j)   = 0._r8
             ccf%dwt_deadcrootc_to_cwdc(c,j)   = 0._r8
          end do
          ccf%dwt_closs(c)  = 0._r8

          if ( use_c13 ) then
             cc13f%dwt_seedc_to_leaf(c)     = 0._r8
             cc13f%dwt_seedc_to_deadstem(c) = 0._r8
             cc13f%dwt_conv_cflux(c)        = 0._r8
             cc13f%dwt_prod10c_gain(c)      = 0._r8
             cc13f%prod10c_loss(c)          = 0._r8
             cc13f%dwt_prod100c_gain(c)     = 0._r8
             cc13f%prod100c_loss(c)         = 0._r8
             do j = 1, nlevdecomp_full
                cc13f%dwt_frootc_to_litr_met_c(c,j) = 0._r8
                cc13f%dwt_frootc_to_litr_cel_c(c,j) = 0._r8
                cc13f%dwt_frootc_to_litr_lig_c(c,j) = 0._r8
                cc13f%dwt_livecrootc_to_cwdc(c,j)   = 0._r8
                cc13f%dwt_deadcrootc_to_cwdc(c,j)   = 0._r8
             end do
             cc13f%dwt_closs(c) = 0._r8
          endif

          if ( use_c14 ) then
             cc14f%dwt_seedc_to_leaf(c)     = 0._r8
             cc14f%dwt_seedc_to_deadstem(c) = 0._r8
             cc14f%dwt_conv_cflux(c)        = 0._r8
             cc14f%dwt_prod10c_gain(c)      = 0._r8
             cc14f%prod10c_loss(c)          = 0._r8
             cc14f%dwt_prod100c_gain(c)     = 0._r8
             cc14f%prod100c_loss(c)         = 0._r8
             do j = 1, nlevdecomp_full
                cc14f%dwt_frootc_to_litr_met_c(c,j) = 0._r8
                cc14f%dwt_frootc_to_litr_cel_c(c,j) = 0._r8
                cc14f%dwt_frootc_to_litr_lig_c(c,j) = 0._r8
                cc14f%dwt_livecrootc_to_cwdc(c,j)   = 0._r8
                cc14f%dwt_deadcrootc_to_cwdc(c,j)   = 0._r8
             end do
             cc14f%dwt_closs(c) = 0._r8
          endif

          cnf%dwt_seedn_to_leaf(c)     = 0._r8
          cnf%dwt_seedn_to_deadstem(c) = 0._r8
          cnf%dwt_conv_nflux(c)        = 0._r8
          cnf%dwt_prod10n_gain(c)      = 0._r8
          cnf%prod10n_loss(c)          = 0._r8
          cnf%dwt_prod100n_gain(c)     = 0._r8
          cnf%prod100n_loss(c)         = 0._r8
          do j = 1, nlevdecomp_full
             cnf%dwt_frootn_to_litr_met_n(c,j) = 0._r8
             cnf%dwt_frootn_to_litr_cel_n(c,j) = 0._r8
             cnf%dwt_frootn_to_litr_lig_n(c,j) = 0._r8
             cnf%dwt_livecrootn_to_cwdn(c,j)   = 0._r8
             cnf%dwt_deadcrootn_to_cwdn(c,j)   = 0._r8
          end do
          cnf%dwt_nloss(c) = 0._r8
       end if
    end do

    ! initialize pft-level variables
    do p = bounds%begp,bounds%endp
       l = pft%landunit(p)
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then

          ! carbon state variables
          if (pft%itype(p) == noveg) then
             pcs%leafc(p)         = 0._r8
             pcs%leafc_storage(p) = 0._r8
          else
             if (pftcon%evergreen(pft%itype(p)) == 1._r8) then
                pcs%leafc(p) = 1._r8
                pcs%leafc_storage(p) = 0._r8
             else if (pft%itype(p) >= npcropmin) then ! prognostic crop types
                pcs%leafc(p) = 0._r8
                pcs%leafc_storage(p) = 0._r8
             else
                pcs%leafc(p) = 0._r8
                pcs%leafc_storage(p) = 1._r8
             end if
          end if
          pcs%leafc_xfer(p) = 0._r8
          if ( crop_prog )then
             pcs%grainc(p)         = 0._r8
             pcs%grainc_storage(p) = 0._r8
             pcs%grainc_xfer(p)    = 0._r8
             pnf%fert(p)           = 0._r8
             pnf%soyfixn(p)        = 0._r8
          end if
          pcs%frootc(p)            = 0._r8
          pcs%frootc_storage(p)    = 0._r8
          pcs%frootc_xfer(p)       = 0._r8
          pcs%livestemc(p)         = 0._r8
          pcs%livestemc_storage(p) = 0._r8
          pcs%livestemc_xfer(p)    = 0._r8

          ! tree types need to be initialized with some stem mass so that
          ! roughness length is not zero in canopy flux calculation

          if (pftcon%woody(pft%itype(p)) == 1._r8) then
             pcs%deadstemc(p) = 0.1_r8
          else
             pcs%deadstemc(p) = 0._r8
          end if

          pcs%deadstemc_storage(p)  = 0._r8
          pcs%deadstemc_xfer(p)     = 0._r8
          pcs%livecrootc(p)         = 0._r8
          pcs%livecrootc_storage(p) = 0._r8
          pcs%livecrootc_xfer(p)    = 0._r8
          pcs%deadcrootc(p)         = 0._r8
          pcs%deadcrootc_storage(p) = 0._r8
          pcs%deadcrootc_xfer(p)    = 0._r8
          pcs%gresp_storage(p)      = 0._r8
          pcs%gresp_xfer(p)         = 0._r8
          pcs%cpool(p)              = 0._r8
          pcs%xsmrpool(p)           = 0._r8
          pcs%pft_ctrunc(p)         = 0._r8
          pcs%dispvegc(p)           = 0._r8
          pcs%storvegc(p)           = 0._r8
          pcs%totpftc(p)            = 0._r8
          pcs%woodc(p)              = 0._r8

          ! calculate totvegc explicitly so that it is available for the isotope 
          ! code on the first time step.

          pcs%totvegc(p)  =                                                                                          &
               pcs%leafc(p) + pcs%leafc_storage(p) + pcs%leafc_xfer(p) + pcs%frootc(p) +                             &
               pcs%frootc_storage(p) + pcs%frootc_xfer(p) + pcs%livestemc(p) + pcs%livestemc_storage(p) +            &
               pcs%livestemc_xfer(p) + pcs%deadstemc(p) + pcs%deadstemc_storage(p) + pcs%deadstemc_xfer(p) +         &
               pcs%livecrootc(p) + pcs%livecrootc_storage(p) + pcs%livecrootc_xfer(p) + pcs%deadcrootc(p) +          &
               pcs%deadcrootc_storage(p) + pcs%deadcrootc_xfer(p) + pcs%gresp_storage(p) +                           &
               pcs%gresp_xfer(p) + pcs%cpool(p)

          if ( crop_prog )then
             pcs%totvegc(p) = pcs%totvegc(p) + pcs%grainc(p) + pcs%grainc_storage(p) + pcs%grainc_xfer(p)
          end if

          if ( use_c13 ) then
             pc13s%leafc(p)              = pcs%leafc(p)               * c13ratio
             pc13s%leafc_storage(p)      = pcs%leafc_storage(p)       * c13ratio
             pc13s%leafc_xfer(p)         = pcs%leafc_xfer(p)          * c13ratio
             pc13s%frootc(p)             = pcs%frootc(p)              * c13ratio
             pc13s%frootc_storage(p)     = pcs%frootc_storage(p)      * c13ratio
             pc13s%frootc_xfer(p)        = pcs%frootc_xfer(p)         * c13ratio
             pc13s%livestemc(p)          = pcs%livestemc(p)           * c13ratio
             pc13s%livestemc_storage(p)  = pcs%livestemc_storage(p)   * c13ratio
             pc13s%livestemc_xfer(p)     = pcs%livestemc_xfer(p)      * c13ratio
             pc13s%deadstemc(p)          = pcs%deadstemc(p)           * c13ratio
             pc13s%deadstemc_storage(p)  = pcs%deadstemc_storage(p)   * c13ratio
             pc13s%deadstemc_xfer(p)     = pcs%deadstemc_xfer(p)      * c13ratio
             pc13s%livecrootc(p)         = pcs%livecrootc(p)          * c13ratio
             pc13s%livecrootc_storage(p) = pcs%livecrootc_storage(p)  * c13ratio
             pc13s%livecrootc_xfer(p)    = pcs%livecrootc_xfer(p)     * c13ratio
             pc13s%deadcrootc(p)         = pcs%deadcrootc(p)          * c13ratio
             pc13s%deadcrootc_storage(p) = pcs%deadcrootc_storage(p)  * c13ratio
             pc13s%deadcrootc_xfer(p)    = pcs%deadcrootc_xfer(p)     * c13ratio
             pc13s%gresp_storage(p)      = pcs%gresp_storage(p)       * c13ratio
             pc13s%gresp_xfer(p)         = pcs%gresp_xfer(p)          * c13ratio
             pc13s%cpool(p)              = pcs%cpool(p)               * c13ratio
             pc13s%xsmrpool(p)           = pcs%xsmrpool(p)            * c13ratio
             pc13s%pft_ctrunc(p)         = pcs%pft_ctrunc(p)          * c13ratio

             ! calculate totvegc explicitly so that it is available for the isotope 
             ! code on the first time step.

             pc13s%totvegc(p) = pc13s%leafc(p) + pc13s%leafc_storage(p) + pc13s%leafc_xfer(p) + pc13s%frootc(p) +       &
                  pc13s%frootc_storage(p) + pc13s%frootc_xfer(p) + pc13s%livestemc(p) + pc13s%livestemc_storage(p) +    &
                  pc13s%livestemc_xfer(p) + pc13s%deadstemc(p) + pc13s%deadstemc_storage(p) + pc13s%deadstemc_xfer(p) + &
                  pc13s%livecrootc(p) + pc13s%livecrootc_storage(p) + pc13s%livecrootc_xfer(p) + pc13s%deadcrootc(p) +  &
                  pc13s%deadcrootc_storage(p) + pc13s%deadcrootc_xfer(p) + pc13s%gresp_storage(p) +                     &
                  pc13s%gresp_xfer(p) + pc13s%cpool(p)
          endif

          if ( use_c14 ) then
             pc14s%leafc(p)              = pcs%leafc(p)               * c14ratio
             pc14s%leafc_storage(p)      = pcs%leafc_storage(p)       * c14ratio
             pc14s%leafc_xfer(p)         = pcs%leafc_xfer(p)          * c14ratio
             pc14s%frootc(p)             = pcs%frootc(p)              * c14ratio
             pc14s%frootc_storage(p)     = pcs%frootc_storage(p)      * c14ratio
             pc14s%frootc_xfer(p)        = pcs%frootc_xfer(p)         * c14ratio
             pc14s%livestemc(p)          = pcs%livestemc(p)           * c14ratio
             pc14s%livestemc_storage(p)  = pcs%livestemc_storage(p)   * c14ratio
             pc14s%livestemc_xfer(p)     = pcs%livestemc_xfer(p)      * c14ratio
             pc14s%deadstemc(p)          = pcs%deadstemc(p)           * c14ratio
             pc14s%deadstemc_storage(p)  = pcs%deadstemc_storage(p)   * c14ratio
             pc14s%deadstemc_xfer(p)     = pcs%deadstemc_xfer(p)      * c14ratio
             pc14s%livecrootc(p)         = pcs%livecrootc(p)          * c14ratio
             pc14s%livecrootc_storage(p) = pcs%livecrootc_storage(p)  * c14ratio
             pc14s%livecrootc_xfer(p)    = pcs%livecrootc_xfer(p)     * c14ratio
             pc14s%deadcrootc(p)         = pcs%deadcrootc(p)          * c14ratio
             pc14s%deadcrootc_storage(p) = pcs%deadcrootc_storage(p)  * c14ratio
             pc14s%deadcrootc_xfer(p)    = pcs%deadcrootc_xfer(p)     * c14ratio
             pc14s%gresp_storage(p)      = pcs%gresp_storage(p)       * c14ratio
             pc14s%gresp_xfer(p)         = pcs%gresp_xfer(p)          * c14ratio
             pc14s%cpool(p)              = pcs%cpool(p)               * c14ratio
             pc14s%xsmrpool(p)           = pcs%xsmrpool(p)            * c14ratio
             pc14s%pft_ctrunc(p)         = pcs%pft_ctrunc(p)          * c14ratio

             ! calculate totvegc explicitly so that it is available for the isotope 
             ! code on the first time step.

             pc14s%totvegc(p)  = pc14s%leafc(p) + pc14s%leafc_storage(p) + pc14s%leafc_xfer(p) + pc14s%frootc(p) +      &
                  pc14s%frootc_storage(p) + pc14s%frootc_xfer(p) + pc14s%livestemc(p) + pc14s%livestemc_storage(p) +    &
                  pc14s%livestemc_xfer(p) + pc14s%deadstemc(p) + pc14s%deadstemc_storage(p) + pc14s%deadstemc_xfer(p) + &
                  pc14s%livecrootc(p) + pc14s%livecrootc_storage(p) + pc14s%livecrootc_xfer(p) + pc14s%deadcrootc(p) +  &
                  pc14s%deadcrootc_storage(p) + pc14s%deadcrootc_xfer(p) + pc14s%gresp_storage(p) +                     &
                  pc14s%gresp_xfer(p) + pc14s%cpool(p)

             pepv%rc14_atm(p) = c14ratio
          endif

          ! nitrogen state variables
          if (pft%itype(p) == noveg) then
             pns%leafn(p) = 0._r8
             pns%leafn_storage(p) = 0._r8
          else
             pns%leafn(p) = pcs%leafc(p) / pftcon%leafcn(pft%itype(p))
             pns%leafn_storage(p) = pcs%leafc_storage(p) / pftcon%leafcn(pft%itype(p))
          end if

          pns%leafn_xfer(p) = 0._r8
          if ( crop_prog )then
             pns%grainn(p) = 0._r8
             pns%grainn_storage(p) = 0._r8
             pns%grainn_xfer(p) = 0._r8
          end if
          pns%frootn(p) = 0._r8
          pns%frootn_storage(p) = 0._r8
          pns%frootn_xfer(p) = 0._r8
          pns%livestemn(p) = 0._r8
          pns%livestemn_storage(p) = 0._r8
          pns%livestemn_xfer(p) = 0._r8

          ! tree types need to be initialized with some stem mass so that
          ! roughness length is not zero in canopy flux calculation

          if (pftcon%woody(pft%itype(p)) == 1._r8) then
             pns%deadstemn(p) = pcs%deadstemc(p) / pftcon%deadwdcn(pft%itype(p))
          else
             pns%deadstemn(p) = 0._r8
          end if

          pns%deadstemn_storage(p) = 0._r8
          pns%deadstemn_xfer(p) = 0._r8
          pns%livecrootn(p) = 0._r8
          pns%livecrootn_storage(p) = 0._r8
          pns%livecrootn_xfer(p) = 0._r8
          pns%deadcrootn(p) = 0._r8
          pns%deadcrootn_storage(p) = 0._r8
          pns%deadcrootn_xfer(p) = 0._r8
          pns%retransn(p) = 0._r8
          pns%npool(p) = 0._r8
          pns%pft_ntrunc(p) = 0._r8
          pns%dispvegn(p) = 0._r8
          pns%storvegn(p) = 0._r8
          pns%totvegn(p)  = 0._r8
          pns%totpftn(p)  = 0._r8

          ! initialization for psnsun and psnsha required for
          ! proper arbitrary initialization of allocation routine
          ! in initial ecosysdyn call

          pcf%psnsun(p) = 0._r8
          pcf%psnsha(p) = 0._r8
          if ( use_c13 ) then
             pc13f%psnsun(p) = 0._r8
             pc13f%psnsha(p) = 0._r8
          endif
          if ( use_c14 ) then
             pc14f%psnsun(p) = 0._r8
             pc14f%psnsha(p) = 0._r8
          endif

          pps%laisun(p) = 0._r8
          pps%laisha(p) = 0._r8

          ! ecophysiological variables
          ! phenology variables

          pepv%dormant_flag(p)   = 1._r8
          pepv%days_active(p)    = 0._r8
          pepv%onset_flag(p)     = 0._r8
          pepv%onset_counter(p)  = 0._r8
          pepv%onset_gddflag(p)  = 0._r8
          pepv%onset_fdd(p)      = 0._r8
          pepv%onset_gdd(p)      = 0._r8
          pepv%onset_swi(p)      = 0._r8
          pepv%offset_flag(p)    = 0._r8
          pepv%offset_counter(p) = 0._r8
          pepv%offset_fdd(p)     = 0._r8
          pepv%offset_swi(p)     = 0._r8
          pepv%lgsf(p)           = 0._r8
          pepv%bglfr(p)          = 0._r8
          pepv%bgtr(p)           = 0._r8
          pepv%annavg_t2m(p)     = 280._r8
          pepv%tempavg_t2m(p)    = 0._r8
          pepv%fert_counter(p)   = 0._r8
          pepv%grain_flag(p)     = 0._r8

          ! non-phenology variables
          pepv%gpp(p)                   = 0._r8
          pepv%availc(p)                = 0._r8
          pepv%xsmrpool_recover(p)      = 0._r8
          pepv%alloc_pnow(p)            = 1._r8
          pepv%c_allometry(p)           = 0._r8
          pepv%n_allometry(p)           = 0._r8
          pepv%plant_ndemand(p)         = 0._r8
          pepv%tempsum_potential_gpp(p) = 0._r8
          pepv%annsum_potential_gpp(p)  = 0._r8
          pepv%tempmax_retransn(p)      = 0._r8
          pepv%annmax_retransn(p)       = 0._r8
          pepv%avail_retransn(p)        = 0._r8
          pepv%plant_nalloc(p)          = 0._r8
          pepv%plant_calloc(p)          = 0._r8
          pepv%excess_cflux(p)          = 0._r8
          pepv%downreg(p)               = 0._r8
          pepv%prev_leafc_to_litter(p)  = 0._r8
          pepv%prev_frootc_to_litter(p) = 0._r8
          pepv%tempsum_npp(p)           = 0._r8
          pepv%annsum_npp(p)            = 0._r8
          if (use_cndv) then
             pepv%tempsum_litfall(p) = 0._r8
             pepv%annsum_litfall(p) = 0._r8
          end if
          if ( use_c13 ) then
             pepv%xsmrpool_c13ratio(p) = c13ratio
             pepv%rc13_canair(p) = 0._r8
             pepv%rc13_psnsun(p) = 0._r8
             pepv%rc13_psnsha(p) = 0._r8
             pps%alphapsnsun(p)  = 0._r8
             pps%alphapsnsha(p)  = 0._r8
          endif

       end if   ! end of if-istsoil block
    end do   ! end of loop over pfts  

    ! initialize column-level fields for special filters
    call CNSetCps(num_special_col, special_col, spval, cps)
    call CNSetCcs(num_special_col, special_col, 0._r8, ccs)
    call CNSetCns(num_special_col, special_col, 0._r8, cns)
    call CNSetCcf(num_special_col, special_col, 0._r8, ccf)
    call CNSetCnf(num_special_col, special_col, 0._r8, cnf)
    if ( use_c13 ) then
       call CNSetCcs(num_special_col, special_col, 0._r8, cc13s)
       call CNSetCcf(num_special_col, special_col, 0._r8, cc13f)
    endif
    if ( use_c14 ) then
       call CNSetCcs(num_special_col, special_col, 0._r8, cc14s)
       call CNSetCcf(num_special_col, special_col, 0._r8, cc14f)
    endif

    ! initialize column-average pft fields for special filters
    call CNSetPps (num_special_col, special_col, spval, pps_a)
    call CNSetPcs (num_special_col, special_col, 0._r8, pcs_a)
    call CNSetPns (num_special_col, special_col, 0._r8, pns_a)
    call CNSetPcf (num_special_col, special_col, 0._r8, pcf_a)
    call CNSetPnf (num_special_col, special_col, 0._r8, pnf_a)

    ! initialize pft-level fields for special filters ucall CNSetPepv
    call CNSetPepv (num_special_pft, special_pft, spval, pepv)
    call CNSetPps  (num_special_pft, special_pft, spval, pps)
    call CNSetPcs  (num_special_pft, special_pft, 0._r8, pcs)
    call CNSetPns  (num_special_pft, special_pft, 0._r8, pns)
    call CNSetPcf  (num_special_pft, special_pft, 0._r8, pcf)
    call CNSetPnf  (num_special_pft, special_pft, 0._r8, pnf)
    if ( use_c13 ) then
       call CNSetPcs(num_special_pft, special_pft, 0._r8, pc13s)
       call CNSetPcf(num_special_pft, special_pft, 0._r8, pc13f)
    endif
    if ( use_c14 ) then
       call CNSetPcs(num_special_pft, special_pft, 0._r8, pc14s)
       call CNSetPcf(num_special_pft, special_pft, 0._r8, pc14f)
    endif

    ! now loop through special filters and explicitly set the variables that
    ! have to be in place for biogeophysics
    ! also set pcf%psnsun and pcf%psnsha to 0 (not included in CNSetPcf())

    ! Loop over special pfts
    do fp = 1,num_special_pft
       p = special_pft(fp)

       pps%tlai(p)               = 0._r8
       pps%tsai(p)               = 0._r8
       pps%elai(p)               = 0._r8
       pps%esai(p)               = 0._r8
       pps%htop(p)               = 0._r8
       pps%hbot(p)               = 0._r8
       pps%fwet(p)               = 0._r8
       pps%fdry(p)               = 0._r8
       pps%frac_veg_nosno_alb(p) = 0._r8
       pps%frac_veg_nosno(p)     = 0._r8
       pcf%psnsun(p)             = 0._r8
       pcf%psnsha(p)             = 0._r8
       if ( use_c13 ) then
          pc13f%psnsun(p) = 0._r8
          pc13f%psnsha(p) = 0._r8
       endif
       if ( use_c14 ) then
          pc14f%psnsun(p) = 0._r8
          pc14f%psnsha(p) = 0._r8
       endif
       if (crop_prog) then
          pnf%fert(p) = 0._r8
       end if

    end do

    ! loop over special columns
    do fc = 1,num_special_col
       c = special_col(fc)

       pcf_a%psnsun(c) = 0._r8
       pcf_a%psnsha(c) = 0._r8
       if ( use_c13 ) then
          pc13f_a%psnsun(c) = 0._r8
          pc13f_a%psnsha(c) = 0._r8
       endif
       if ( use_c14 ) then
          pc14f_a%psnsun(c) = 0._r8
          pc14f_a%psnsha(c) = 0._r8
       endif
       
       ccs%seedc(c)         = 0._r8
       ccs%prod10c(c)       = 0._r8	  
       ccs%prod100c(c)      = 0._r8	  
       ccs%totprodc(c)      = 0._r8	  
       if ( use_c13 ) then
          cc13s%seedc(c)    = 0._r8
          cc13s%prod10c(c)  = 0._r8	  
          cc13s%prod100c(c) = 0._r8	  
          cc13s%totprodc(c) = 0._r8	  
       endif
       if ( use_c14 ) then
          cc14s%seedc(c)    = 0._r8
          cc14s%prod10c(c)  = 0._r8	  
          cc14s%prod100c(c) = 0._r8	  
          cc14s%totprodc(c) = 0._r8	  
       endif
       cns%seedn(c)                 = 0._r8
       cns%prod10n(c)               = 0._r8	  
       cns%prod100n(c)              = 0._r8	  
       cns%totprodn(c)              = 0._r8	  

       ccf%dwt_seedc_to_leaf(c)                          = 0._r8
       ccf%dwt_seedc_to_deadstem(c)                      = 0._r8
       ccf%dwt_conv_cflux(c)                             = 0._r8
       ccf%lf_conv_cflux(c)                              = 0._r8    
       ccf%dwt_prod10c_gain(c)                           = 0._r8
       ccf%prod10c_loss(c)                               = 0._r8
       ccf%dwt_prod100c_gain(c)                          = 0._r8
       ccf%prod100c_loss(c)                              = 0._r8
       ccf%dwt_closs(c)                                  = 0._r8
       ccf%landuseflux(c)                                = 0._r8
       ccf%landuptake(c)                                 = 0._r8
       ccf%dwt_frootc_to_litr_met_c(c,1:nlevdecomp_full) = 0._r8
       ccf%dwt_frootc_to_litr_cel_c(c,1:nlevdecomp_full) = 0._r8
       ccf%dwt_frootc_to_litr_lig_c(c,1:nlevdecomp_full) = 0._r8
       ccf%dwt_livecrootc_to_cwdc(c,1:nlevdecomp_full)   = 0._r8
       ccf%dwt_deadcrootc_to_cwdc(c,1:nlevdecomp_full)   = 0._r8
       if ( use_c13 ) then
          cc13f%dwt_seedc_to_leaf(c)                          = 0._r8
          cc13f%dwt_seedc_to_deadstem(c)                      = 0._r8
          cc13f%dwt_conv_cflux(c)                             = 0._r8
          cc13f%dwt_prod10c_gain(c)                           = 0._r8
          cc13f%prod10c_loss(c)                               = 0._r8
          cc13f%dwt_prod100c_gain(c)                          = 0._r8
          cc13f%prod100c_loss(c)                              = 0._r8
          cc13f%dwt_frootc_to_litr_met_c(c,1:nlevdecomp_full) = 0._r8
          cc13f%dwt_frootc_to_litr_cel_c(c,1:nlevdecomp_full) = 0._r8
          cc13f%dwt_frootc_to_litr_lig_c(c,1:nlevdecomp_full) = 0._r8
          cc13f%dwt_livecrootc_to_cwdc(c,1:nlevdecomp_full)   = 0._r8
          cc13f%dwt_deadcrootc_to_cwdc(c,1:nlevdecomp_full)   = 0._r8
          cc13f%dwt_closs(c)                                  = 0._r8
       endif
       if ( use_c14 ) then
          cc14f%dwt_seedc_to_leaf(c)                          = 0._r8
          cc14f%dwt_seedc_to_deadstem(c)                      = 0._r8
          cc14f%dwt_conv_cflux(c)                             = 0._r8
          cc14f%dwt_prod10c_gain(c)                           = 0._r8
          cc14f%prod10c_loss(c)                               = 0._r8
          cc14f%dwt_prod100c_gain(c)                          = 0._r8
          cc14f%prod100c_loss(c)                              = 0._r8
          cc14f%dwt_frootc_to_litr_met_c(c,1:nlevdecomp_full) = 0._r8
          cc14f%dwt_frootc_to_litr_cel_c(c,1:nlevdecomp_full) = 0._r8
          cc14f%dwt_frootc_to_litr_lig_c(c,1:nlevdecomp_full) = 0._r8
          cc14f%dwt_livecrootc_to_cwdc(c,1:nlevdecomp_full)   = 0._r8
          cc14f%dwt_deadcrootc_to_cwdc(c,1:nlevdecomp_full)   = 0._r8
          cc14f%dwt_closs(c)                                  = 0._r8
       endif
       cnf%dwt_seedn_to_leaf(c)                          = 0._r8
       cnf%dwt_seedn_to_deadstem(c)                      = 0._r8
       cnf%dwt_conv_nflux(c)                             = 0._r8
       cnf%dwt_prod10n_gain(c)                           = 0._r8
       cnf%prod10n_loss(c)                               = 0._r8
       cnf%dwt_prod100n_gain(c)                          = 0._r8
       cnf%prod100n_loss(c)                              = 0._r8
       cnf%dwt_frootn_to_litr_met_n(c,1:nlevdecomp_full) = 0._r8
       cnf%dwt_frootn_to_litr_cel_n(c,1:nlevdecomp_full) = 0._r8
       cnf%dwt_frootn_to_litr_lig_n(c,1:nlevdecomp_full) = 0._r8
       cnf%dwt_livecrootn_to_cwdn(c,1:nlevdecomp_full)   = 0._r8
       cnf%dwt_deadcrootn_to_cwdn(c,1:nlevdecomp_full)   = 0._r8
       cnf%dwt_nloss(c)                                  = 0._r8

    end do ! end loop over special columns

  end subroutine initColdCN

end module CNInitMod
