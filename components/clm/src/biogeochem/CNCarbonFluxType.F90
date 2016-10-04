module CNCarbonFluxType

  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_infnan_mod         , only : nan => shr_infnan_nan, assignment(=)
  use decompMod              , only : bounds_type
  use clm_varpar             , only : ndecomp_cascade_transitions, ndecomp_pools, nlevcan
  use clm_varpar             , only : crop_prog
  use clm_varpar             , only : nlevdecomp_full, nlevgrnd, nlevdecomp
  use clm_varcon             , only : spval, ispval, dzsoi_decomp
  use landunit_varcon        , only : istsoil, istcrop, istdlak 
  use clm_varctl             , only : use_cndv, use_c13, use_ed 
  use ch4varcon              , only : allowlakeprod
  use pftvarcon              , only : npcropmin
  use CNDecompCascadeConType , only : decomp_cascade_con
  use PatchType              , only : pft                
  use ColumnType             , only : col                
  use LandunitType           , only : lun
  use clm_varctl             , only : nu_com
  ! bgc interface & pflotran
  use clm_varctl             , only : use_bgc_interface, use_pflotran, pf_cmode, use_vertsoilc
  ! 
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! NOTE(bandre, 2013-10) according to Charlie Koven, nfix_timeconst
  ! is currently used as a flag and rate constant. Rate constant: time
  ! over which to exponentially relax the npp flux for N fixation term
  ! flag: (if  <=  0. or  >=  365; use old annual method). Default value is
  ! junk that should always be overwritten by the namelist or init function!
  !
  ! (days) time over which to exponentially relax the npp flux for N fixation term
  real(r8), public :: nfix_timeconst = -1.2345_r8 
  !
  type, public :: carbonflux_type

     ! gap mortality fluxes
     real(r8), pointer :: m_leafc_to_litter_patch                   (:)     ! leaf C mortality (gC/m2/s)
     real(r8), pointer :: m_leafc_storage_to_litter_patch           (:)     ! leaf C storage mortality (gC/m2/s)
     real(r8), pointer :: m_leafc_xfer_to_litter_patch              (:)     ! leaf C transfer mortality (gC/m2/s)
     real(r8), pointer :: m_frootc_to_litter_patch                  (:)     ! fine root C mortality (gC/m2/s)
     real(r8), pointer :: m_frootc_storage_to_litter_patch          (:)     ! fine root C storage mortality (gC/m2/s)
     real(r8), pointer :: m_frootc_xfer_to_litter_patch             (:)     ! fine root C transfer mortality (gC/m2/s)
     real(r8), pointer :: m_livestemc_to_litter_patch               (:)     ! live stem C mortality (gC/m2/s)
     real(r8), pointer :: m_livestemc_storage_to_litter_patch       (:)     ! live stem C storage mortality (gC/m2/s)
     real(r8), pointer :: m_livestemc_xfer_to_litter_patch          (:)     ! live stem C transfer mortality (gC/m2/s)
     real(r8), pointer :: m_deadstemc_to_litter_patch               (:)     ! dead stem C mortality (gC/m2/s)
     real(r8), pointer :: m_deadstemc_storage_to_litter_patch       (:)     ! dead stem C storage mortality (gC/m2/s)
     real(r8), pointer :: m_deadstemc_xfer_to_litter_patch          (:)     ! dead stem C transfer mortality (gC/m2/s)
     real(r8), pointer :: m_livecrootc_to_litter_patch              (:)     ! live coarse root C mortality (gC/m2/s)
     real(r8), pointer :: m_livecrootc_storage_to_litter_patch      (:)     ! live coarse root C storage mortality (gC/m2/s)
     real(r8), pointer :: m_livecrootc_xfer_to_litter_patch         (:)     ! live coarse root C transfer mortality (gC/m2/s)
     real(r8), pointer :: m_deadcrootc_to_litter_patch              (:)     ! dead coarse root C mortality (gC/m2/s)
     real(r8), pointer :: m_deadcrootc_storage_to_litter_patch      (:)     ! dead coarse root C storage mortality (gC/m2/s)
     real(r8), pointer :: m_deadcrootc_xfer_to_litter_patch         (:)     ! dead coarse root C transfer mortality (gC/m2/s)
     real(r8), pointer :: m_gresp_storage_to_litter_patch           (:)     ! growth respiration storage mortality (gC/m2/s)
     real(r8), pointer :: m_gresp_xfer_to_litter_patch              (:)     ! growth respiration transfer mortality (gC/m2/s)

     ! harvest mortality fluxes
     real(r8), pointer :: hrv_leafc_to_litter_patch                 (:)     ! leaf C harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_leafc_storage_to_litter_patch         (:)     ! leaf C storage harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_leafc_xfer_to_litter_patch            (:)     ! leaf C transfer harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_frootc_to_litter_patch                (:)     ! fine root C harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_frootc_storage_to_litter_patch        (:)     ! fine root C storage harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_frootc_xfer_to_litter_patch           (:)     ! fine root C transfer harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_livestemc_to_litter_patch             (:)     ! live stem C harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_livestemc_storage_to_litter_patch     (:)     ! live stem C storage harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_livestemc_xfer_to_litter_patch        (:)     ! live stem C transfer harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_deadstemc_to_prod10c_patch            (:)     ! dead stem C harvest to 10-year product pool (gC/m2/s)
     real(r8), pointer :: hrv_deadstemc_to_prod100c_patch           (:)     ! dead stem C harvest to 100-year product pool (gC/m2/s)
     real(r8), pointer :: hrv_deadstemc_storage_to_litter_patch     (:)     ! dead stem C storage harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_deadstemc_xfer_to_litter_patch        (:)     ! dead stem C transfer harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_livecrootc_to_litter_patch            (:)     ! live coarse root C harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_livecrootc_storage_to_litter_patch    (:)     ! live coarse root C storage harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_livecrootc_xfer_to_litter_patch       (:)     ! live coarse root C transfer harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_deadcrootc_to_litter_patch            (:)     ! dead coarse root C harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_deadcrootc_storage_to_litter_patch    (:)     ! dead coarse root C storage harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_deadcrootc_xfer_to_litter_patch       (:)     ! dead coarse root C transfer harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_gresp_storage_to_litter_patch         (:)     ! growth respiration storage harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_gresp_xfer_to_litter_patch            (:)     ! growth respiration transfer harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_xsmrpool_to_atm_patch                 (:)     ! excess MR pool harvest mortality (gC/m2/s)
     ! crop harvest
     real(r8), pointer :: hrv_leafc_to_prod1c_patch                 (:)     ! crop leafc harvested (gC/m2/s)
     real(r8), pointer :: hrv_livestemc_to_prod1c_patch             (:)     ! crop stemc harvested (gC/m2/s)
     real(r8), pointer :: hrv_grainc_to_prod1c_patch                (:)     ! crop grain harvested (gC/m2/s)
     real(r8), pointer :: hrv_cropc_to_prod1c_patch                 (:)     ! total amount of crop C harvested (gC/m2/s)

     ! fire C fluxes 
     real(r8), pointer :: m_leafc_to_fire_patch                     (:)     ! (gC/m2/s) fire C emissions from leafc 
     real(r8), pointer :: m_leafc_storage_to_fire_patch             (:)     ! (gC/m2/s) fire C emissions from leafc_storage             
     real(r8), pointer :: m_leafc_xfer_to_fire_patch                (:)     ! (gC/m2/s) fire C emissions from leafc_xfer
     real(r8), pointer :: m_livestemc_to_fire_patch                 (:)     ! (gC/m2/s) fire C emissions from livestemc
     real(r8), pointer :: m_livestemc_storage_to_fire_patch         (:)     ! (gC/m2/s) fire C emissions from livestemc_storage       
     real(r8), pointer :: m_livestemc_xfer_to_fire_patch            (:)     ! (gC/m2/s) fire C emissions from livestemc_xfer
     real(r8), pointer :: m_deadstemc_to_fire_patch                 (:)     ! (gC/m2/s) fire C emissions from deadstemc_xfer
     real(r8), pointer :: m_deadstemc_storage_to_fire_patch         (:)     ! (gC/m2/s) fire C emissions from deadstemc_storage         
     real(r8), pointer :: m_deadstemc_xfer_to_fire_patch            (:)     ! (gC/m2/s) fire C emissions from deadstemc_xfer
     real(r8), pointer :: m_frootc_to_fire_patch                    (:)     ! (gC/m2/s) fire C emissions from frootc
     real(r8), pointer :: m_frootc_storage_to_fire_patch            (:)     ! (gC/m2/s) fire C emissions from frootc_storage
     real(r8), pointer :: m_frootc_xfer_to_fire_patch               (:)     ! (gC/m2/s) fire C emissions from frootc_xfer
     real(r8), pointer :: m_livecrootc_to_fire_patch                (:)     ! (gC/m2/s) fire C emissions from livecrootc
     real(r8), pointer :: m_livecrootc_storage_to_fire_patch        (:)     ! (gC/m2/s) fire C emissions from livecrootc_storage     
     real(r8), pointer :: m_livecrootc_xfer_to_fire_patch           (:)     ! (gC/m2/s) fire C emissions from livecrootc_xfer
     real(r8), pointer :: m_deadcrootc_to_fire_patch                (:)     ! (gC/m2/s) fire C emissions from deadcrootc
     real(r8), pointer :: m_deadcrootc_storage_to_fire_patch        (:)     ! (gC/m2/s) fire C emissions from deadcrootc_storage 
     real(r8), pointer :: m_deadcrootc_xfer_to_fire_patch           (:)     ! (gC/m2/s) fire C emissions from deadcrootc_xfer
     real(r8), pointer :: m_gresp_storage_to_fire_patch             (:)     ! (gC/m2/s) fire C emissions from gresp_storage 
     real(r8), pointer :: m_gresp_xfer_to_fire_patch                (:)     ! (gC/m2/s) fire C emissions from gresp_xfer
     real(r8), pointer :: m_leafc_to_litter_fire_patch              (:)     ! (gC/m2/s) from leafc to litter c due to fire
     real(r8), pointer :: m_leafc_storage_to_litter_fire_patch      (:)     ! (gC/m2/s) from leafc_storage to litter C  due to fire               
     real(r8), pointer :: m_leafc_xfer_to_litter_fire_patch         (:)     ! (gC/m2/s) from leafc_xfer to litter C  due to fire               
     real(r8), pointer :: m_livestemc_to_litter_fire_patch          (:)     ! (gC/m2/s) from livestemc to litter C  due to fire               
     real(r8), pointer :: m_livestemc_storage_to_litter_fire_patch  (:)     ! (gC/m2/s) from livestemc_storage to litter C due to fire      
     real(r8), pointer :: m_livestemc_xfer_to_litter_fire_patch     (:)     ! (gC/m2/s) from livestemc_xfer to litter C due to fire      
     real(r8), pointer :: m_livestemc_to_deadstemc_fire_patch       (:)     ! (gC/m2/s) from livestemc to deadstemc due to fire       
     real(r8), pointer :: m_deadstemc_to_litter_fire_patch          (:)     ! (gC/m2/s) from deadstemc to litter C due to fire      
     real(r8), pointer :: m_deadstemc_storage_to_litter_fire_patch  (:)     ! (gC/m2/s) from deadstemc_storage to litter C due to fire               
     real(r8), pointer :: m_deadstemc_xfer_to_litter_fire_patch     (:)     ! (gC/m2/s) from deadstemc_xfer to litter C due to fire               
     real(r8), pointer :: m_frootc_to_litter_fire_patch             (:)     ! (gC/m2/s) from frootc to litter C due to fire               
     real(r8), pointer :: m_frootc_storage_to_litter_fire_patch     (:)     ! (gC/m2/s) from frootc_storage to litter C due to fire               
     real(r8), pointer :: m_frootc_xfer_to_litter_fire_patch        (:)     ! (gC/m2/s) from frootc_xfer to litter C due to fire               
     real(r8), pointer :: m_livecrootc_to_litter_fire_patch         (:)     ! (gC/m2/s) from livecrootc to litter C due to fire                     
     real(r8), pointer :: m_livecrootc_storage_to_litter_fire_patch (:)     ! (gC/m2/s) from livecrootc_storage to litter C due to fire                     
     real(r8), pointer :: m_livecrootc_xfer_to_litter_fire_patch    (:)     ! (gC/m2/s) from livecrootc_xfer to litter C due to fire                     
     real(r8), pointer :: m_livecrootc_to_deadcrootc_fire_patch     (:)     ! (gC/m2/s) from livecrootc to deadstemc due to fire        
     real(r8), pointer :: m_deadcrootc_to_litter_fire_patch         (:)     ! (gC/m2/s) from deadcrootc to litter C due to fire                       
     real(r8), pointer :: m_deadcrootc_storage_to_litter_fire_patch (:)     ! (gC/m2/s) from deadcrootc_storage to litter C due to fire                       
     real(r8), pointer :: m_deadcrootc_xfer_to_litter_fire_patch    (:)     ! (gC/m2/s) from deadcrootc_xfer to litter C due to fire                       
     real(r8), pointer :: m_gresp_storage_to_litter_fire_patch      (:)     ! (gC/m2/s) from gresp_storage to litter C due to fire                       
     real(r8), pointer :: m_gresp_xfer_to_litter_fire_patch         (:)     ! (gC/m2/s) from gresp_xfer to litter C due to fire                       

     ! phenology fluxes from transfer pools                     
     real(r8), pointer :: grainc_xfer_to_grainc_patch               (:)     ! grain C growth from storage for prognostic crop(gC/m2/s)
     real(r8), pointer :: leafc_xfer_to_leafc_patch                 (:)     ! leaf C growth from storage (gC/m2/s)
     real(r8), pointer :: frootc_xfer_to_frootc_patch               (:)     ! fine root C growth from storage (gC/m2/s)
     real(r8), pointer :: livestemc_xfer_to_livestemc_patch         (:)     ! live stem C growth from storage (gC/m2/s)
     real(r8), pointer :: deadstemc_xfer_to_deadstemc_patch         (:)     ! dead stem C growth from storage (gC/m2/s)
     real(r8), pointer :: livecrootc_xfer_to_livecrootc_patch       (:)     ! live coarse root C growth from storage (gC/m2/s)
     real(r8), pointer :: deadcrootc_xfer_to_deadcrootc_patch       (:)     ! dead coarse root C growth from storage (gC/m2/s)

     ! leaf and fine root litterfall fluxes                          
     real(r8), pointer :: leafc_to_litter_patch                     (:)     ! leaf C litterfall (gC/m2/s)
     real(r8), pointer :: frootc_to_litter_patch                    (:)     ! fine root C litterfall (gC/m2/s)
     real(r8), pointer :: livestemc_to_litter_patch                 (:)     ! live stem C litterfall (gC/m2/s)
     real(r8), pointer :: grainc_to_food_patch                      (:)     ! grain C to food for prognostic crop(gC/m2/s)

     ! maintenance respiration fluxes                          
     real(r8), pointer :: leaf_mr_patch                             (:)     ! leaf maintenance respiration (gC/m2/s)
     real(r8), pointer :: froot_mr_patch                            (:)     ! fine root maintenance respiration (gC/m2/s)
     real(r8), pointer :: livestem_mr_patch                         (:)     ! live stem maintenance respiration (gC/m2/s)
     real(r8), pointer :: livecroot_mr_patch                        (:)     ! live coarse root maintenance respiration (gC/m2/s)
     real(r8), pointer :: grain_mr_patch                            (:)     ! crop grain or organs maint. respiration (gC/m2/s)
     real(r8), pointer :: leaf_curmr_patch                          (:)     ! leaf maintenance respiration from current GPP (gC/m2/s)
     real(r8), pointer :: froot_curmr_patch                         (:)     ! fine root maintenance respiration from current GPP (gC/m2/s)
     real(r8), pointer :: livestem_curmr_patch                      (:)     ! live stem maintenance respiration from current GPP (gC/m2/s)
     real(r8), pointer :: livecroot_curmr_patch                     (:)     ! live coarse root maintenance respiration from current GPP (gC/m2/s)
     real(r8), pointer :: grain_curmr_patch                         (:)     ! crop grain or organs maint. respiration from current GPP (gC/m2/s)
     real(r8), pointer :: leaf_xsmr_patch                           (:)     ! leaf maintenance respiration from storage (gC/m2/s)
     real(r8), pointer :: froot_xsmr_patch                          (:)     ! fine root maintenance respiration from storage (gC/m2/s)
     real(r8), pointer :: livestem_xsmr_patch                       (:)     ! live stem maintenance respiration from storage (gC/m2/s)
     real(r8), pointer :: livecroot_xsmr_patch                      (:)     ! live coarse root maintenance respiration from storage (gC/m2/s)
     real(r8), pointer :: grain_xsmr_patch                          (:)     ! crop grain or organs maint. respiration from storage (gC/m2/s)

     ! photosynthesis fluxes                                   
     real(r8), pointer :: psnsun_to_cpool_patch                     (:)     ! C fixation from sunlit canopy (gC/m2/s)
     real(r8), pointer :: psnshade_to_cpool_patch                   (:)     ! C fixation from shaded canopy (gC/m2/s)

     ! allocation fluxes, from current GPP                     
     real(r8), pointer :: cpool_to_xsmrpool_patch                   (:)     ! allocation to maintenance respiration storage pool (gC/m2/s)
     real(r8), pointer :: cpool_to_grainc_patch                     (:)     ! allocation to grain C for prognostic crop(gC/m2/s)
     real(r8), pointer :: cpool_to_grainc_storage_patch             (:)     ! allocation to grain C storage for prognostic crop(gC/m2/s)
     real(r8), pointer :: cpool_to_leafc_patch                      (:)     ! allocation to leaf C (gC/m2/s)
     real(r8), pointer :: cpool_to_leafc_storage_patch              (:)     ! allocation to leaf C storage (gC/m2/s)
     real(r8), pointer :: cpool_to_frootc_patch                     (:)     ! allocation to fine root C (gC/m2/s)
     real(r8), pointer :: cpool_to_frootc_storage_patch             (:)     ! allocation to fine root C storage (gC/m2/s)
     real(r8), pointer :: cpool_to_livestemc_patch                  (:)     ! allocation to live stem C (gC/m2/s)
     real(r8), pointer :: cpool_to_livestemc_storage_patch          (:)     ! allocation to live stem C storage (gC/m2/s)
     real(r8), pointer :: cpool_to_deadstemc_patch                  (:)     ! allocation to dead stem C (gC/m2/s)
     real(r8), pointer :: cpool_to_deadstemc_storage_patch          (:)     ! allocation to dead stem C storage (gC/m2/s)
     real(r8), pointer :: cpool_to_livecrootc_patch                 (:)     ! allocation to live coarse root C (gC/m2/s)
     real(r8), pointer :: cpool_to_livecrootc_storage_patch         (:)     ! allocation to live coarse root C storage (gC/m2/s)
     real(r8), pointer :: cpool_to_deadcrootc_patch                 (:)     ! allocation to dead coarse root C (gC/m2/s)
     real(r8), pointer :: cpool_to_deadcrootc_storage_patch         (:)     ! allocation to dead coarse root C storage (gC/m2/s)
     real(r8), pointer :: cpool_to_gresp_storage_patch              (:)     ! allocation to growth respiration storage (gC/m2/s)

     ! growth respiration fluxes                               
     real(r8), pointer :: xsmrpool_to_atm_patch                     (:)     ! excess MR pool harvest mortality (gC/m2/s)
     real(r8), pointer :: cpool_leaf_gr_patch                       (:)     ! leaf growth respiration (gC/m2/s)
     real(r8), pointer :: cpool_leaf_storage_gr_patch               (:)     ! leaf growth respiration to storage (gC/m2/s)
     real(r8), pointer :: transfer_leaf_gr_patch                    (:)     ! leaf growth respiration from storage (gC/m2/s)
     real(r8), pointer :: cpool_froot_gr_patch                      (:)     ! fine root growth respiration (gC/m2/s)
     real(r8), pointer :: cpool_froot_storage_gr_patch              (:)     ! fine root  growth respiration to storage (gC/m2/s)
     real(r8), pointer :: transfer_froot_gr_patch                   (:)     ! fine root  growth respiration from storage (gC/m2/s)
     real(r8), pointer :: cpool_livestem_gr_patch                   (:)     ! live stem growth respiration (gC/m2/s)
     real(r8), pointer :: cpool_livestem_storage_gr_patch           (:)     ! live stem growth respiration to storage (gC/m2/s)
     real(r8), pointer :: transfer_livestem_gr_patch                (:)     ! live stem growth respiration from storage (gC/m2/s)
     real(r8), pointer :: cpool_deadstem_gr_patch                   (:)     ! dead stem growth respiration (gC/m2/s)
     real(r8), pointer :: cpool_deadstem_storage_gr_patch           (:)     ! dead stem growth respiration to storage (gC/m2/s)
     real(r8), pointer :: transfer_deadstem_gr_patch                (:)     ! dead stem growth respiration from storage (gC/m2/s)
     real(r8), pointer :: cpool_livecroot_gr_patch                  (:)     ! live coarse root growth respiration (gC/m2/s)
     real(r8), pointer :: cpool_livecroot_storage_gr_patch          (:)     ! live coarse root growth respiration to storage (gC/m2/s)
     real(r8), pointer :: transfer_livecroot_gr_patch               (:)     ! live coarse root growth respiration from storage (gC/m2/s)
     real(r8), pointer :: cpool_deadcroot_gr_patch                  (:)     ! dead coarse root growth respiration (gC/m2/s)
     real(r8), pointer :: cpool_deadcroot_storage_gr_patch          (:)     ! dead coarse root growth respiration to storage (gC/m2/s)
     real(r8), pointer :: transfer_deadcroot_gr_patch               (:)     ! dead coarse root growth respiration from storage (gC/m2/s)

     ! growth respiration for prognostic crop model
     real(r8), pointer :: cpool_grain_gr_patch                      (:)     ! grain growth respiration (gC/m2/s)
     real(r8), pointer :: cpool_grain_storage_gr_patch              (:)     ! grain growth respiration to storage (gC/m2/s)
     real(r8), pointer :: transfer_grain_gr_patch                   (:)     ! grain growth respiration from storage (gC/m2/s)

     ! annual turnover of storage to transfer pools            
     real(r8), pointer :: grainc_storage_to_xfer_patch              (:)     ! grain C shift storage to transfer for prognostic crop model (gC/m2/s)
     real(r8), pointer :: leafc_storage_to_xfer_patch               (:)     ! leaf C shift storage to transfer (gC/m2/s)
     real(r8), pointer :: frootc_storage_to_xfer_patch              (:)     ! fine root C shift storage to transfer (gC/m2/s)
     real(r8), pointer :: livestemc_storage_to_xfer_patch           (:)     ! live stem C shift storage to transfer (gC/m2/s)
     real(r8), pointer :: deadstemc_storage_to_xfer_patch           (:)     ! dead stem C shift storage to transfer (gC/m2/s)
     real(r8), pointer :: livecrootc_storage_to_xfer_patch          (:)     ! live coarse root C shift storage to transfer (gC/m2/s)
     real(r8), pointer :: deadcrootc_storage_to_xfer_patch          (:)     ! dead coarse root C shift storage to transfer (gC/m2/s)
     real(r8), pointer :: gresp_storage_to_xfer_patch               (:)     ! growth respiration shift storage to transfer (gC/m2/s)

     ! turnover of livewood to deadwood
     real(r8), pointer :: livestemc_to_deadstemc_patch              (:)     ! live stem C turnover (gC/m2/s)
     real(r8), pointer :: livecrootc_to_deadcrootc_patch            (:)     ! live coarse root C turnover (gC/m2/s)

     ! summary (diagnostic) flux variables, not involved in mass balance
     real(r8), pointer :: gpp_patch                                 (:)     ! (gC/m2/s) gross primary production 
     real(r8), pointer :: gpp_before_downreg_patch                  (:)     ! (gC/m2/s) gross primary production before down regulation
     real(r8), pointer :: mr_patch                                  (:)     ! (gC/m2/s) maintenance respiration
     real(r8), pointer :: current_gr_patch                          (:)     ! (gC/m2/s) growth resp for new growth displayed in this timestep
     real(r8), pointer :: transfer_gr_patch                         (:)     ! (gC/m2/s) growth resp for transfer growth displayed in this timestep
     real(r8), pointer :: storage_gr_patch                          (:)     ! (gC/m2/s) growth resp for growth sent to storage for later display
     real(r8), pointer :: gr_patch                                  (:)     ! (gC/m2/s) total growth respiration
     real(r8), pointer :: ar_patch                                  (:)     ! (gC/m2/s) autotrophic respiration (MR + GR)
     real(r8), pointer :: rr_patch                                  (:)     ! (gC/m2/s) root respiration (fine root MR + total root GR)
     real(r8), pointer :: npp_patch                                 (:)     ! (gC/m2/s) net primary production
     real(r8), pointer :: agnpp_patch                               (:)     ! (gC/m2/s) aboveground NPP
     real(r8), pointer :: bgnpp_patch                               (:)     ! (gC/m2/s) belowground NPP
     real(r8), pointer :: litfall_patch                             (:)     ! (gC/m2/s) litterfall (leaves and fine roots)
     real(r8), pointer :: vegfire_patch                             (:)     ! (gC/m2/s) patch-level fire loss (obsolete, mark for removal)
     real(r8), pointer :: wood_harvestc_patch                       (:)     ! (gC/m2/s) patch-level wood harvest (to product pools)
     real(r8), pointer :: cinputs_patch                             (:)     ! (gC/m2/s) patch-level carbon inputs (for balance checking)
     real(r8), pointer :: coutputs_patch                            (:)     ! (gC/m2/s) patch-level carbon outputs (for balance checking)

     real(r8), pointer :: plant_calloc_patch                        (:)     ! total allocated C flux (gC/m2/s)
     real(r8), pointer :: excess_cflux_patch                        (:)     ! C flux not allocated due to downregulation (gC/m2/s)
     real(r8), pointer :: prev_leafc_to_litter_patch                (:)     ! previous timestep leaf C litterfall flux (gC/m2/s)
     real(r8), pointer :: prev_frootc_to_litter_patch               (:)     ! previous timestep froot C litterfall flux (gC/m2/s)
     real(r8), pointer :: availc_patch                              (:)     ! C flux available for allocation (gC/m2/s)
     real(r8), pointer :: xsmrpool_recover_patch                    (:)     ! C flux assigned to recovery of negative cpool (gC/m2/s)
     real(r8), pointer :: xsmrpool_c13ratio_patch                   (:)     ! C13/C(12+13) ratio for xsmrpool (proportion)

     ! CN: CLAMP summary (diagnostic) variables, not involved in mass balance
     real(r8), pointer :: frootc_alloc_patch                        (:)     ! (gC/m2/s) patch-level fine root C alloc
     real(r8), pointer :: frootc_loss_patch                         (:)     ! (gC/m2/s) patch-level fine root C loss
     real(r8), pointer :: leafc_alloc_patch                         (:)     ! (gC/m2/s) patch-level leaf C alloc
     real(r8), pointer :: leafc_loss_patch                          (:)     ! (gC/m2/s) patch-level leaf C loss
     real(r8), pointer :: woodc_alloc_patch                         (:)     ! (gC/m2/s) patch-level wood C alloc
     real(r8), pointer :: woodc_loss_patch                          (:)     ! (gC/m2/s) patch-level wood C loss

     ! fire code
     real(r8), pointer :: fire_closs_patch                          (:)     ! (gC/m2/s) total patch-level fire C loss 

     ! For aerenchyma calculations in CH4 code
     real(r8), pointer :: annavg_agnpp_patch                        (:)     ! (gC/m2/s) annual average aboveground NPP
     real(r8), pointer :: annavg_bgnpp_patch                        (:)     ! (gC/m2/s) annual average belowground NPP
     real(r8), pointer :: tempavg_agnpp_patch                       (:)     ! (gC/m2/s) temp. average aboveground NPP
     real(r8), pointer :: tempavg_bgnpp_patch                       (:)     ! (gC/m2/s) temp. average belowground NPP

     !----------------------------------------------------
     ! column carbon flux variables  
     !----------------------------------------------------

     ! phenology: litterfall and crop fluxes
     real(r8), pointer :: phenology_c_to_litr_met_c_col             (:,:)   ! C fluxes associated with phenology (litterfall and crop) to litter metabolic pool (gC/m3/s)
     real(r8), pointer :: phenology_c_to_litr_cel_c_col             (:,:)   ! C fluxes associated with phenology (litterfall and crop) to litter cellulose pool (gC/m3/s)
     real(r8), pointer :: phenology_c_to_litr_lig_c_col             (:,:)   ! C fluxes associated with phenology (litterfall and crop) to litter lignin pool (gC/m3/s)

     ! gap mortality
     real(r8), pointer :: gap_mortality_c_to_litr_met_c_col         (:,:)   ! C fluxes associated with gap mortality to litter metabolic pool (gC/m3/s)
     real(r8), pointer :: gap_mortality_c_to_litr_cel_c_col         (:,:)   ! C fluxes associated with gap mortality to litter cellulose pool (gC/m3/s)
     real(r8), pointer :: gap_mortality_c_to_litr_lig_c_col         (:,:)   ! C fluxes associated with gap mortality to litter lignin pool (gC/m3/s)
     real(r8), pointer :: gap_mortality_c_to_cwdc_col               (:,:)   ! C fluxes associated with gap mortality to CWD pool (gC/m3/s)

     ! fire
     real(r8), pointer :: fire_mortality_c_to_cwdc_col              (:,:)   ! C fluxes associated with fire mortality to CWD pool (gC/m3/s)

     ! harvest
     real(r8), pointer :: harvest_c_to_litr_met_c_col               (:,:)   ! C fluxes associated with harvest to litter metabolic pool (gC/m3/s)
     real(r8), pointer :: harvest_c_to_litr_cel_c_col               (:,:)   ! C fluxes associated with harvest to litter cellulose pool (gC/m3/s)
     real(r8), pointer :: harvest_c_to_litr_lig_c_col               (:,:)   ! C fluxes associated with harvest to litter lignin pool (gC/m3/s)
     real(r8), pointer :: harvest_c_to_cwdc_col                     (:,:)   ! C fluxes associated with harvest to CWD pool (gC/m3/s)

     ! new variables for CN code
     real(r8), pointer :: hrv_deadstemc_to_prod10c_col              (:)     ! dead stem C harvest mortality to 10-year product pool (gC/m2/s)        
     real(r8), pointer :: hrv_deadstemc_to_prod100c_col             (:)     ! dead stem C harvest mortality to 100-year product pool (gC/m2/s)        
     real(r8), pointer :: hrv_cropc_to_prod1c_col                   (:)     ! crop C harvest mortality to 1-year product pool (gC/m2/s)

     ! column-level fire fluxes
     real(r8), pointer :: m_decomp_cpools_to_fire_vr_col            (:,:,:) ! vertically-resolved decomposing C fire loss (gC/m3/s)
     real(r8), pointer :: m_decomp_cpools_to_fire_col               (:,:)   ! vertically-integrated (diagnostic) decomposing C fire loss (gC/m2/s)
     real(r8), pointer :: m_c_to_litr_met_fire_col                  (:,:)   ! C from leaf, froot, xfer and storage C to litter labile C by fire (gC/m3/s) 
     real(r8), pointer :: m_c_to_litr_cel_fire_col                  (:,:)   ! C from leaf, froot, xfer and storage C to litter cellulose C by fire (gC/m3/s) 
     real(r8), pointer :: m_c_to_litr_lig_fire_col                  (:,:)   ! C from leaf, froot, xfer and storage C to litter lignin C by fire (gC/m3/s) 
     real(r8), pointer :: lf_conv_cflux_col                         (:)     ! (gC/m2/s) conversion C flux due to BET and BDT area decreasing (immediate loss to atm)
     real(r8), pointer :: somc_fire_col                             (:)     ! (gC/m2/s) carbon emissions due to peat burning

     ! decomposition fluxes
     real(r8), pointer :: decomp_cpools_sourcesink_col              (:,:,:) ! change in decomposing c pools. Used to update concentrations concurrently with vertical transport (gC/m3/timestep)  
     real(r8), pointer :: decomp_cascade_hr_vr_col                  (:,:,:) ! vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
     real(r8), pointer :: decomp_cascade_hr_col                     (:,:)   ! vertically-integrated (diagnostic) het. resp. from decomposing C pools (gC/m2/s)
     real(r8), pointer :: decomp_cascade_ctransfer_vr_col           (:,:,:) ! vertically-resolved C transferred along deomposition cascade (gC/m3/s)
     real(r8), pointer :: decomp_cascade_ctransfer_col              (:,:)   ! vertically-integrated (diagnostic) C transferred along deomposition cascade (gC/m2/s)
     real(r8), pointer :: decomp_k_col                              (:,:,:) ! rate constant for decomposition (1./sec)
     real(r8), pointer :: hr_vr_col                                 (:,:)   ! total vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
     real(r8), pointer :: o_scalar_col                              (:,:)   ! fraction by which decomposition is limited by anoxia
     real(r8), pointer :: w_scalar_col                              (:,:)   ! fraction by which decomposition is limited by moisture availability
     real(r8), pointer :: t_scalar_col                              (:,:)   ! fraction by which decomposition is limited by temperature
     real(r8), pointer :: som_c_leached_col                         (:)     ! total SOM C loss from vertical transport (gC/m^2/s)
     real(r8), pointer :: decomp_cpools_leached_col                 (:,:)   ! C loss from vertical transport from each decomposing C pool (gC/m^2/s)
     real(r8), pointer :: decomp_cpools_transport_tendency_col      (:,:,:) ! C tendency due to vertical transport in decomposing C pools (gC/m^3/s)

     ! nitrif_denitrif
     real(r8), pointer :: phr_vr_col                                (:,:)   ! potential hr (not N-limited) (gC/m3/s)
     real(r8), pointer :: fphr_col                                  (:,:)   ! fraction of potential heterotrophic respiration

     ! CN dynamic landcover fluxes
     real(r8), pointer :: dwt_seedc_to_leaf_col                     (:)     ! (gC/m2/s) seed source to patch-level
     real(r8), pointer :: dwt_seedc_to_deadstem_col                 (:)     ! (gC/m2/s) seed source to patch-level
     real(r8), pointer :: dwt_conv_cflux_col                        (:)     ! (gC/m2/s) conversion C flux (immediate loss to atm)
     real(r8), pointer :: dwt_prod10c_gain_col                      (:)     ! (gC/m2/s) addition to 10-yr wood product pool
     real(r8), pointer :: dwt_prod100c_gain_col                     (:)     ! (gC/m2/s) addition to 100-yr wood product pool
     real(r8), pointer :: dwt_frootc_to_litr_met_c_col              (:,:)   ! (gC/m3/s) fine root to litter due to landcover change
     real(r8), pointer :: dwt_frootc_to_litr_cel_c_col              (:,:)   ! (gC/m3/s) fine root to litter due to landcover change
     real(r8), pointer :: dwt_frootc_to_litr_lig_c_col              (:,:)   ! (gC/m3/s) fine root to litter due to landcover change
     real(r8), pointer :: dwt_livecrootc_to_cwdc_col                (:,:)   ! (gC/m3/s) live coarse root to CWD due to landcover change
     real(r8), pointer :: dwt_deadcrootc_to_cwdc_col                (:,:)   ! (gC/m3/s) dead coarse root to CWD due to landcover change
     real(r8), pointer :: dwt_closs_col                             (:)     ! (gC/m2/s) total carbon loss from product pools and conversion
     real(r8), pointer :: landuseflux_col                           (:)     ! (gC/m2/s) dwt_closs+product_closs
     real(r8), pointer :: landuptake_col                            (:)     ! (gC/m2/s) nee-landuseflux

     ! CN wood product pool loss fluxes
     real(r8), pointer :: prod1c_loss_col                           (:)     ! (gC/m2/s) decomposition loss from 1-year product pool
     real(r8), pointer :: prod10c_loss_col                          (:)     ! (gC/m2/s) decomposition loss from 10-yr wood product pool
     real(r8), pointer :: prod100c_loss_col                         (:)     ! (gC/m2/s) decomposition loss from 100-yr wood product pool
     real(r8), pointer :: product_closs_col                         (:)     ! (gC/m2/s) total wood product carbon loss

     ! summary (diagnostic) flux variables, not involved in mass balance
     real(r8), pointer :: lithr_col                                 (:)     ! (gC/m2/s) litter heterotrophic respiration 
     real(r8), pointer :: somhr_col                                 (:)     ! (gC/m2/s) soil organic matter heterotrophic respiration
     real(r8), pointer :: hr_col                                    (:)     ! (gC/m2/s) total heterotrophic respiration
     real(r8), pointer :: sr_col                                    (:)     ! (gC/m2/s) total soil respiration (HR + root resp)
     real(r8), pointer :: er_col                                    (:)     ! (gC/m2/s) total ecosystem respiration, autotrophic + heterotrophic
     real(r8), pointer :: litfire_col                               (:)     ! (gC/m2/s) litter fire losses
     real(r8), pointer :: somfire_col                               (:)     ! (gC/m2/s) soil organic matter fire losses
     real(r8), pointer :: totfire_col                               (:)     ! (gC/m2/s) total ecosystem fire losses
     real(r8), pointer :: nep_col                                   (:)     ! (gC/m2/s) net ecosystem production, excludes fire, landuse, and harvest flux, positive for sink
     real(r8), pointer :: nbp_col                                   (:)     ! (gC/m2/s) net biome production, includes fire, landuse, and harvest flux, positive for sink
     real(r8), pointer :: nee_col                                   (:)     ! (gC/m2/s) net ecosystem exchange of carbon, includes fire, landuse, harvest, and hrv_xsmrpool flux, positive for source

     ! CN CLAMP summary (diagnostic) flux variables, not involved in mass balance
     real(r8), pointer :: cwdc_hr_col                               (:)     ! (gC/m2/s) col-level coarse woody debris C heterotrophic respiration
     real(r8), pointer :: cwdc_loss_col                             (:)     ! (gC/m2/s) col-level coarse woody debris C loss
     real(r8), pointer :: litterc_loss_col                          (:)     ! (gC/m2/s) col-level litter C loss

     real(r8), pointer :: bgc_cpool_ext_inputs_vr_col               (:, :, :)  ! col-level extneral organic carbon input gC/m3 /time step
     real(r8), pointer :: bgc_cpool_ext_loss_vr_col                 (:, :, :)  ! col-level extneral organic carbon loss gC/m3 /time step
     ! patch averaged to column variables - to remove need for pcf_a instance
     real(r8), pointer :: rr_col                                    (:)     ! column (gC/m2/s) root respiration (fine root MR + total root GR) (p2c)
     real(r8), pointer :: ar_col                                    (:)     ! column (gC/m2/s) autotrophic respiration (MR + GR) (p2c)      
     real(r8), pointer :: gpp_col                                   (:)     ! column (gC/m2/s) GPP flux before downregulation  (p2c)         
     real(r8), pointer :: npp_col                                   (:)     ! column (gC/m2/s) net primary production (p2c)                  
     real(r8), pointer :: fire_closs_p2c_col                        (:)     ! column (gC/m2/s) patch2col averaged column-level fire C loss (p2c)
     real(r8), pointer :: fire_closs_col                            (:)     ! column (gC/m2/s) total patch-level fire C loss 
     real(r8), pointer :: litfall_col                               (:)     ! column (gC/m2/s) total patch-level litterfall C loss (p2c)       
     real(r8), pointer :: vegfire_col                               (:)     ! column (gC/m2/s) patch-level fire loss (obsolete, mark for removal) (p2c)
     real(r8), pointer :: wood_harvestc_col                         (:)     ! column (p2c)                                                  
     real(r8), pointer :: hrv_xsmrpool_to_atm_col                   (:)     ! column excess MR pool harvest mortality (gC/m2/s) (p2c)

     ! Temporary and annual sums
     real(r8), pointer :: tempsum_litfall_patch       (:) ! temporary annual sum of litfall (gC/m2/yr) (CNDV)
     real(r8), pointer :: tempsum_npp_patch           (:) ! patch temporary annual sum of NPP (gC/m2/yr)
     real(r8), pointer :: annsum_litfall_patch        (:) ! annual sum of litfall (gC/m2/yr) (CNDV)
     real(r8), pointer :: annsum_npp_patch            (:) ! patch annual sum of NPP (gC/m2/yr)
     real(r8), pointer :: annsum_npp_col              (:) ! col annual sum of NPP, averaged from pft-level (gC/m2/yr)
     real(r8), pointer :: lag_npp_col                 (:) ! col lagged net primary production (gC/m2/s)
     
     ! debug
     real(r8), pointer :: plant_to_litter_cflux		  (:) ! for the purpose of mass balance check
     real(r8), pointer :: plant_to_cwd_cflux		  (:) ! for the purpose of mass balance check
     real(r8), pointer :: allocation_leaf 		  (:) ! check allocation to leaf for dynamic allocation scheme
     real(r8), pointer :: allocation_stem 		  (:) ! check allocation to stem for dynamic allocation scheme
     real(r8), pointer :: allocation_froot 		  (:) ! check allocation to fine root for dynamic allocation scheme

     ! new variables for clm_bgc_interface & pflotran
     !------------------------------------------------------------------------
     real(r8), pointer :: externalc_to_decomp_cpools_col            (:,:,:) ! col (gC/m3/s) net C fluxes associated with litter/som-adding/removal to decomp pools
                                                                            ! (sum of all external C additions and removals, excluding decomposition/hr).
     real(r8), pointer :: externalc_to_decomp_delta_col             (:)     ! col (gC/m2) summarized net change of whole column C i/o to decomposing pool bwtn time-step
     real(r8), pointer :: f_co2_soil_vr_col                         (:,:)   ! total vertically-resolved soil-atm. CO2 exchange (gC/m3/s)
     real(r8), pointer :: f_co2_soil_col                            (:)     ! total soil-atm. CO2 exchange (gC/m2/s)
    !------------------------------------------------------------------------

   contains

     procedure , public  :: Init   
     procedure , public  :: Restart
     procedure , public  :: SetValues
     procedure , public  :: ZeroDWT
     procedure , public  :: Summary
     procedure , public  :: summary_cflux_for_ch4
     procedure , public  :: summary_rr
     procedure , private :: InitAllocate 
     procedure , private :: InitHistory
     procedure , private :: InitCold
     ! bgc & pflotran interface
     procedure , private :: CSummary_interface
  end type carbonflux_type
  !------------------------------------------------------------------------

contains
   
  !------------------------------------------------------------------------
  subroutine Init(this, bounds, carbon_type)

     class(carbonflux_type) :: this
     type(bounds_type), intent(in) :: bounds  
     character(len=3) , intent(in) :: carbon_type ! one of ['c12', c13','c14']

     call this%InitAllocate ( bounds)
     call this%InitHistory ( bounds, carbon_type )
     call this%InitCold (bounds )

   end subroutine Init

   !------------------------------------------------------------------------
   subroutine InitAllocate(this, bounds)
     !
     ! !ARGUMENTS:
     class (carbonflux_type) :: this 
     type(bounds_type), intent(in)    :: bounds 
     !
     ! !LOCAL VARIABLES:
     integer           :: begp,endp
     integer           :: begc,endc
     !------------------------------------------------------------------------

     begp = bounds%begp; endp = bounds%endp
     begc = bounds%begc; endc = bounds%endc

     allocate(this%m_leafc_to_litter_patch                   (begp:endp)) ; this%m_leafc_to_litter_patch                   (:) = nan
     allocate(this%m_frootc_to_litter_patch                  (begp:endp)) ; this%m_frootc_to_litter_patch                  (:) = nan
     allocate(this%m_leafc_storage_to_litter_patch           (begp:endp)) ; this%m_leafc_storage_to_litter_patch           (:) = nan
     allocate(this%m_frootc_storage_to_litter_patch          (begp:endp)) ; this%m_frootc_storage_to_litter_patch          (:) = nan
     allocate(this%m_livestemc_storage_to_litter_patch       (begp:endp)) ; this%m_livestemc_storage_to_litter_patch       (:) = nan
     allocate(this%m_deadstemc_storage_to_litter_patch       (begp:endp)) ; this%m_deadstemc_storage_to_litter_patch       (:) = nan
     allocate(this%m_livecrootc_storage_to_litter_patch      (begp:endp)) ; this%m_livecrootc_storage_to_litter_patch      (:) = nan
     allocate(this%m_deadcrootc_storage_to_litter_patch      (begp:endp)) ; this%m_deadcrootc_storage_to_litter_patch      (:) = nan
     allocate(this%m_leafc_xfer_to_litter_patch              (begp:endp)) ; this%m_leafc_xfer_to_litter_patch              (:) = nan
     allocate(this%m_frootc_xfer_to_litter_patch             (begp:endp)) ; this%m_frootc_xfer_to_litter_patch             (:) = nan
     allocate(this%m_livestemc_xfer_to_litter_patch          (begp:endp)) ; this%m_livestemc_xfer_to_litter_patch          (:) = nan
     allocate(this%m_deadstemc_xfer_to_litter_patch          (begp:endp)) ; this%m_deadstemc_xfer_to_litter_patch          (:) = nan
     allocate(this%m_livecrootc_xfer_to_litter_patch         (begp:endp)) ; this%m_livecrootc_xfer_to_litter_patch         (:) = nan
     allocate(this%m_deadcrootc_xfer_to_litter_patch         (begp:endp)) ; this%m_deadcrootc_xfer_to_litter_patch         (:) = nan
     allocate(this%m_livestemc_to_litter_patch               (begp:endp)) ; this%m_livestemc_to_litter_patch               (:) = nan
     allocate(this%m_deadstemc_to_litter_patch               (begp:endp)) ; this%m_deadstemc_to_litter_patch               (:) = nan
     allocate(this%m_livecrootc_to_litter_patch              (begp:endp)) ; this%m_livecrootc_to_litter_patch              (:) = nan
     allocate(this%m_deadcrootc_to_litter_patch              (begp:endp)) ; this%m_deadcrootc_to_litter_patch              (:) = nan
     allocate(this%m_gresp_storage_to_litter_patch           (begp:endp)) ; this%m_gresp_storage_to_litter_patch           (:) = nan
     allocate(this%m_gresp_xfer_to_litter_patch              (begp:endp)) ; this%m_gresp_xfer_to_litter_patch              (:) = nan
     allocate(this%hrv_leafc_to_litter_patch                 (begp:endp)) ; this%hrv_leafc_to_litter_patch                 (:) = nan
     allocate(this%hrv_leafc_storage_to_litter_patch         (begp:endp)) ; this%hrv_leafc_storage_to_litter_patch         (:) = nan
     allocate(this%hrv_leafc_xfer_to_litter_patch            (begp:endp)) ; this%hrv_leafc_xfer_to_litter_patch            (:) = nan
     allocate(this%hrv_frootc_to_litter_patch                (begp:endp)) ; this%hrv_frootc_to_litter_patch                (:) = nan
     allocate(this%hrv_frootc_storage_to_litter_patch        (begp:endp)) ; this%hrv_frootc_storage_to_litter_patch        (:) = nan
     allocate(this%hrv_frootc_xfer_to_litter_patch           (begp:endp)) ; this%hrv_frootc_xfer_to_litter_patch           (:) = nan
     allocate(this%hrv_livestemc_to_litter_patch             (begp:endp)) ; this%hrv_livestemc_to_litter_patch             (:) = nan
     allocate(this%hrv_livestemc_storage_to_litter_patch     (begp:endp)) ; this%hrv_livestemc_storage_to_litter_patch     (:) = nan
     allocate(this%hrv_livestemc_xfer_to_litter_patch        (begp:endp)) ; this%hrv_livestemc_xfer_to_litter_patch        (:) = nan
     allocate(this%hrv_deadstemc_to_prod10c_patch            (begp:endp)) ; this%hrv_deadstemc_to_prod10c_patch            (:) = nan
     allocate(this%hrv_deadstemc_to_prod100c_patch           (begp:endp)) ; this%hrv_deadstemc_to_prod100c_patch           (:) = nan
     allocate(this%hrv_leafc_to_prod1c_patch                 (begp:endp)) ; this%hrv_leafc_to_prod1c_patch                 (:) = nan
     allocate(this%hrv_livestemc_to_prod1c_patch             (begp:endp)) ; this%hrv_livestemc_to_prod1c_patch             (:) = nan
     allocate(this%hrv_grainc_to_prod1c_patch                (begp:endp)) ; this%hrv_grainc_to_prod1c_patch                (:) = nan
     allocate(this%hrv_cropc_to_prod1c_patch                 (begp:endp)) ; this%hrv_cropc_to_prod1c_patch                 (:) = nan
     allocate(this%hrv_deadstemc_storage_to_litter_patch     (begp:endp)) ; this%hrv_deadstemc_storage_to_litter_patch     (:) = nan
     allocate(this%hrv_deadstemc_xfer_to_litter_patch        (begp:endp)) ; this%hrv_deadstemc_xfer_to_litter_patch        (:) = nan
     allocate(this%hrv_livecrootc_to_litter_patch            (begp:endp)) ; this%hrv_livecrootc_to_litter_patch            (:) = nan
     allocate(this%hrv_livecrootc_storage_to_litter_patch    (begp:endp)) ; this%hrv_livecrootc_storage_to_litter_patch    (:) = nan
     allocate(this%hrv_livecrootc_xfer_to_litter_patch       (begp:endp)) ; this%hrv_livecrootc_xfer_to_litter_patch       (:) = nan
     allocate(this%hrv_deadcrootc_to_litter_patch            (begp:endp)) ; this%hrv_deadcrootc_to_litter_patch            (:) = nan
     allocate(this%hrv_deadcrootc_storage_to_litter_patch    (begp:endp)) ; this%hrv_deadcrootc_storage_to_litter_patch    (:) = nan
     allocate(this%hrv_deadcrootc_xfer_to_litter_patch       (begp:endp)) ; this%hrv_deadcrootc_xfer_to_litter_patch       (:) = nan
     allocate(this%hrv_gresp_storage_to_litter_patch         (begp:endp)) ; this%hrv_gresp_storage_to_litter_patch         (:) = nan
     allocate(this%hrv_gresp_xfer_to_litter_patch            (begp:endp)) ; this%hrv_gresp_xfer_to_litter_patch            (:) = nan
     allocate(this%hrv_xsmrpool_to_atm_patch                 (begp:endp)) ; this%hrv_xsmrpool_to_atm_patch                 (:) = nan
     allocate(this%m_leafc_to_fire_patch                     (begp:endp)) ; this%m_leafc_to_fire_patch                     (:) = nan
     allocate(this%m_leafc_storage_to_fire_patch             (begp:endp)) ; this%m_leafc_storage_to_fire_patch             (:) = nan
     allocate(this%m_leafc_xfer_to_fire_patch                (begp:endp)) ; this%m_leafc_xfer_to_fire_patch                (:) = nan
     allocate(this%m_livestemc_to_fire_patch                 (begp:endp)) ; this%m_livestemc_to_fire_patch                 (:) = nan
     allocate(this%m_livestemc_storage_to_fire_patch         (begp:endp)) ; this%m_livestemc_storage_to_fire_patch         (:) = nan
     allocate(this%m_livestemc_xfer_to_fire_patch            (begp:endp)) ; this%m_livestemc_xfer_to_fire_patch            (:) = nan
     allocate(this%m_deadstemc_to_fire_patch                 (begp:endp)) ; this%m_deadstemc_to_fire_patch                 (:) = nan
     allocate(this%m_deadstemc_storage_to_fire_patch         (begp:endp)) ; this%m_deadstemc_storage_to_fire_patch         (:) = nan
     allocate(this%m_deadstemc_xfer_to_fire_patch            (begp:endp)) ; this%m_deadstemc_xfer_to_fire_patch            (:) = nan
     allocate(this%m_frootc_to_fire_patch                    (begp:endp)) ; this%m_frootc_to_fire_patch                    (:) = nan
     allocate(this%m_frootc_storage_to_fire_patch            (begp:endp)) ; this%m_frootc_storage_to_fire_patch            (:) = nan
     allocate(this%m_frootc_xfer_to_fire_patch               (begp:endp)) ; this%m_frootc_xfer_to_fire_patch               (:) = nan
     allocate(this%m_livecrootc_to_fire_patch                (begp:endp)) ; this%m_livecrootc_to_fire_patch                (:) = nan
     allocate(this%m_livecrootc_storage_to_fire_patch        (begp:endp)) ; this%m_livecrootc_storage_to_fire_patch        (:) = nan
     allocate(this%m_livecrootc_xfer_to_fire_patch           (begp:endp)) ; this%m_livecrootc_xfer_to_fire_patch           (:) = nan
     allocate(this%m_deadcrootc_to_fire_patch                (begp:endp)) ; this%m_deadcrootc_to_fire_patch                (:) = nan
     allocate(this%m_deadcrootc_storage_to_fire_patch        (begp:endp)) ; this%m_deadcrootc_storage_to_fire_patch        (:) = nan
     allocate(this%m_deadcrootc_xfer_to_fire_patch           (begp:endp)) ; this%m_deadcrootc_xfer_to_fire_patch           (:) = nan
     allocate(this%m_gresp_storage_to_fire_patch             (begp:endp)) ; this%m_gresp_storage_to_fire_patch             (:) = nan
     allocate(this%m_gresp_xfer_to_fire_patch                (begp:endp)) ; this%m_gresp_xfer_to_fire_patch                (:) = nan
     allocate(this%m_leafc_to_litter_fire_patch              (begp:endp)) ; this%m_leafc_to_litter_fire_patch              (:) = nan
     allocate(this%m_leafc_storage_to_litter_fire_patch      (begp:endp)) ; this%m_leafc_storage_to_litter_fire_patch      (:) = nan
     allocate(this%m_leafc_xfer_to_litter_fire_patch         (begp:endp)) ; this%m_leafc_xfer_to_litter_fire_patch         (:) = nan
     allocate(this%m_livestemc_to_litter_fire_patch          (begp:endp)) ; this%m_livestemc_to_litter_fire_patch          (:) = nan
     allocate(this%m_livestemc_storage_to_litter_fire_patch  (begp:endp)) ; this%m_livestemc_storage_to_litter_fire_patch  (:) = nan
     allocate(this%m_livestemc_xfer_to_litter_fire_patch     (begp:endp)) ; this%m_livestemc_xfer_to_litter_fire_patch     (:) = nan
     allocate(this%m_livestemc_to_deadstemc_fire_patch       (begp:endp)) ; this%m_livestemc_to_deadstemc_fire_patch       (:) = nan
     allocate(this%m_deadstemc_to_litter_fire_patch          (begp:endp)) ; this%m_deadstemc_to_litter_fire_patch          (:) = nan
     allocate(this%m_deadstemc_storage_to_litter_fire_patch  (begp:endp)) ; this%m_deadstemc_storage_to_litter_fire_patch  (:) = nan
     allocate(this%m_deadstemc_xfer_to_litter_fire_patch     (begp:endp)) ; this%m_deadstemc_xfer_to_litter_fire_patch     (:) = nan
     allocate(this%m_frootc_to_litter_fire_patch             (begp:endp)) ; this%m_frootc_to_litter_fire_patch             (:) = nan
     allocate(this%m_frootc_storage_to_litter_fire_patch     (begp:endp)) ; this%m_frootc_storage_to_litter_fire_patch     (:) = nan
     allocate(this%m_frootc_xfer_to_litter_fire_patch        (begp:endp)) ; this%m_frootc_xfer_to_litter_fire_patch        (:) = nan
     allocate(this%m_livecrootc_to_litter_fire_patch         (begp:endp)) ; this%m_livecrootc_to_litter_fire_patch         (:) = nan
     allocate(this%m_livecrootc_storage_to_litter_fire_patch (begp:endp)) ; this%m_livecrootc_storage_to_litter_fire_patch (:) = nan
     allocate(this%m_livecrootc_xfer_to_litter_fire_patch    (begp:endp)) ; this%m_livecrootc_xfer_to_litter_fire_patch    (:) = nan
     allocate(this%m_livecrootc_to_deadcrootc_fire_patch     (begp:endp)) ; this%m_livecrootc_to_deadcrootc_fire_patch     (:) = nan
     allocate(this%m_deadcrootc_to_litter_fire_patch         (begp:endp)) ; this%m_deadcrootc_to_litter_fire_patch         (:) = nan
     allocate(this%m_deadcrootc_storage_to_litter_fire_patch (begp:endp)) ; this%m_deadcrootc_storage_to_litter_fire_patch (:) = nan
     allocate(this%m_deadcrootc_xfer_to_litter_fire_patch    (begp:endp)) ; this%m_deadcrootc_xfer_to_litter_fire_patch    (:) = nan
     allocate(this%m_gresp_storage_to_litter_fire_patch      (begp:endp)) ; this%m_gresp_storage_to_litter_fire_patch      (:) = nan
     allocate(this%m_gresp_xfer_to_litter_fire_patch         (begp:endp)) ; this%m_gresp_xfer_to_litter_fire_patch         (:) = nan
     allocate(this%leafc_xfer_to_leafc_patch                 (begp:endp)) ; this%leafc_xfer_to_leafc_patch                 (:) = nan
     allocate(this%frootc_xfer_to_frootc_patch               (begp:endp)) ; this%frootc_xfer_to_frootc_patch               (:) = nan
     allocate(this%livestemc_xfer_to_livestemc_patch         (begp:endp)) ; this%livestemc_xfer_to_livestemc_patch         (:) = nan
     allocate(this%deadstemc_xfer_to_deadstemc_patch         (begp:endp)) ; this%deadstemc_xfer_to_deadstemc_patch         (:) = nan
     allocate(this%livecrootc_xfer_to_livecrootc_patch       (begp:endp)) ; this%livecrootc_xfer_to_livecrootc_patch       (:) = nan
     allocate(this%deadcrootc_xfer_to_deadcrootc_patch       (begp:endp)) ; this%deadcrootc_xfer_to_deadcrootc_patch       (:) = nan
     allocate(this%leafc_to_litter_patch                     (begp:endp)) ; this%leafc_to_litter_patch                     (:) = nan
     allocate(this%frootc_to_litter_patch                    (begp:endp)) ; this%frootc_to_litter_patch                    (:) = nan
     allocate(this%leaf_mr_patch                             (begp:endp)) ; this%leaf_mr_patch                             (:) = nan
     allocate(this%froot_mr_patch                            (begp:endp)) ; this%froot_mr_patch                            (:) = nan
     allocate(this%livestem_mr_patch                         (begp:endp)) ; this%livestem_mr_patch                         (:) = nan
     allocate(this%livecroot_mr_patch                        (begp:endp)) ; this%livecroot_mr_patch                        (:) = nan
     allocate(this%grain_mr_patch                            (begp:endp)) ; this%grain_mr_patch                            (:) = nan
     allocate(this%leaf_curmr_patch                          (begp:endp)) ; this%leaf_curmr_patch                          (:) = nan
     allocate(this%froot_curmr_patch                         (begp:endp)) ; this%froot_curmr_patch                         (:) = nan
     allocate(this%livestem_curmr_patch                      (begp:endp)) ; this%livestem_curmr_patch                      (:) = nan
     allocate(this%livecroot_curmr_patch                     (begp:endp)) ; this%livecroot_curmr_patch                     (:) = nan
     allocate(this%grain_curmr_patch                         (begp:endp)) ; this%grain_curmr_patch                         (:) = nan
     allocate(this%leaf_xsmr_patch                           (begp:endp)) ; this%leaf_xsmr_patch                           (:) = nan
     allocate(this%froot_xsmr_patch                          (begp:endp)) ; this%froot_xsmr_patch                          (:) = nan
     allocate(this%livestem_xsmr_patch                       (begp:endp)) ; this%livestem_xsmr_patch                       (:) = nan
     allocate(this%livecroot_xsmr_patch                      (begp:endp)) ; this%livecroot_xsmr_patch                      (:) = nan
     allocate(this%grain_xsmr_patch                          (begp:endp)) ; this%grain_xsmr_patch                          (:) = nan
     allocate(this%psnsun_to_cpool_patch                     (begp:endp)) ; this%psnsun_to_cpool_patch                     (:) = nan
     allocate(this%psnshade_to_cpool_patch                   (begp:endp)) ; this%psnshade_to_cpool_patch                   (:) = nan
     allocate(this%cpool_to_xsmrpool_patch                   (begp:endp)) ; this%cpool_to_xsmrpool_patch                   (:) = nan
     allocate(this%cpool_to_leafc_patch                      (begp:endp)) ; this%cpool_to_leafc_patch                      (:) = nan
     allocate(this%cpool_to_leafc_storage_patch              (begp:endp)) ; this%cpool_to_leafc_storage_patch              (:) = nan
     allocate(this%cpool_to_frootc_patch                     (begp:endp)) ; this%cpool_to_frootc_patch                     (:) = nan
     allocate(this%cpool_to_frootc_storage_patch             (begp:endp)) ; this%cpool_to_frootc_storage_patch             (:) = nan
     allocate(this%cpool_to_livestemc_patch                  (begp:endp)) ; this%cpool_to_livestemc_patch                  (:) = nan
     allocate(this%cpool_to_livestemc_storage_patch          (begp:endp)) ; this%cpool_to_livestemc_storage_patch          (:) = nan
     allocate(this%cpool_to_deadstemc_patch                  (begp:endp)) ; this%cpool_to_deadstemc_patch                  (:) = nan
     allocate(this%cpool_to_deadstemc_storage_patch          (begp:endp)) ; this%cpool_to_deadstemc_storage_patch          (:) = nan
     allocate(this%cpool_to_livecrootc_patch                 (begp:endp)) ; this%cpool_to_livecrootc_patch                 (:) = nan
     allocate(this%cpool_to_livecrootc_storage_patch         (begp:endp)) ; this%cpool_to_livecrootc_storage_patch         (:) = nan
     allocate(this%cpool_to_deadcrootc_patch                 (begp:endp)) ; this%cpool_to_deadcrootc_patch                 (:) = nan
     allocate(this%cpool_to_deadcrootc_storage_patch         (begp:endp)) ; this%cpool_to_deadcrootc_storage_patch         (:) = nan
     allocate(this%cpool_to_gresp_storage_patch              (begp:endp)) ; this%cpool_to_gresp_storage_patch              (:) = nan
     allocate(this%cpool_leaf_gr_patch                       (begp:endp)) ; this%cpool_leaf_gr_patch                       (:) = nan
     allocate(this%cpool_leaf_storage_gr_patch               (begp:endp)) ; this%cpool_leaf_storage_gr_patch               (:) = nan
     allocate(this%transfer_leaf_gr_patch                    (begp:endp)) ; this%transfer_leaf_gr_patch                    (:) = nan
     allocate(this%cpool_froot_gr_patch                      (begp:endp)) ; this%cpool_froot_gr_patch                      (:) = nan
     allocate(this%cpool_froot_storage_gr_patch              (begp:endp)) ; this%cpool_froot_storage_gr_patch              (:) = nan
     allocate(this%transfer_froot_gr_patch                   (begp:endp)) ; this%transfer_froot_gr_patch                   (:) = nan
     allocate(this%cpool_livestem_gr_patch                   (begp:endp)) ; this%cpool_livestem_gr_patch                   (:) = nan
     allocate(this%cpool_livestem_storage_gr_patch           (begp:endp)) ; this%cpool_livestem_storage_gr_patch           (:) = nan
     allocate(this%transfer_livestem_gr_patch                (begp:endp)) ; this%transfer_livestem_gr_patch                (:) = nan
     allocate(this%cpool_deadstem_gr_patch                   (begp:endp)) ; this%cpool_deadstem_gr_patch                   (:) = nan
     allocate(this%cpool_deadstem_storage_gr_patch           (begp:endp)) ; this%cpool_deadstem_storage_gr_patch           (:) = nan
     allocate(this%transfer_deadstem_gr_patch                (begp:endp)) ; this%transfer_deadstem_gr_patch                (:) = nan
     allocate(this%cpool_livecroot_gr_patch                  (begp:endp)) ; this%cpool_livecroot_gr_patch                  (:) = nan
     allocate(this%cpool_livecroot_storage_gr_patch          (begp:endp)) ; this%cpool_livecroot_storage_gr_patch          (:) = nan
     allocate(this%transfer_livecroot_gr_patch               (begp:endp)) ; this%transfer_livecroot_gr_patch               (:) = nan
     allocate(this%cpool_deadcroot_gr_patch                  (begp:endp)) ; this%cpool_deadcroot_gr_patch                  (:) = nan
     allocate(this%cpool_deadcroot_storage_gr_patch          (begp:endp)) ; this%cpool_deadcroot_storage_gr_patch          (:) = nan
     allocate(this%transfer_deadcroot_gr_patch               (begp:endp)) ; this%transfer_deadcroot_gr_patch               (:) = nan
     allocate(this%leafc_storage_to_xfer_patch               (begp:endp)) ; this%leafc_storage_to_xfer_patch               (:) = nan
     allocate(this%frootc_storage_to_xfer_patch              (begp:endp)) ; this%frootc_storage_to_xfer_patch              (:) = nan
     allocate(this%livestemc_storage_to_xfer_patch           (begp:endp)) ; this%livestemc_storage_to_xfer_patch           (:) = nan
     allocate(this%deadstemc_storage_to_xfer_patch           (begp:endp)) ; this%deadstemc_storage_to_xfer_patch           (:) = nan
     allocate(this%livecrootc_storage_to_xfer_patch          (begp:endp)) ; this%livecrootc_storage_to_xfer_patch          (:) = nan
     allocate(this%deadcrootc_storage_to_xfer_patch          (begp:endp)) ; this%deadcrootc_storage_to_xfer_patch          (:) = nan
     allocate(this%gresp_storage_to_xfer_patch               (begp:endp)) ; this%gresp_storage_to_xfer_patch               (:) = nan
     allocate(this%livestemc_to_deadstemc_patch              (begp:endp)) ; this%livestemc_to_deadstemc_patch              (:) = nan
     allocate(this%livecrootc_to_deadcrootc_patch            (begp:endp)) ; this%livecrootc_to_deadcrootc_patch            (:) = nan
     allocate(this%mr_patch                                  (begp:endp)) ; this%mr_patch                                  (:) = nan
     allocate(this%current_gr_patch                          (begp:endp)) ; this%current_gr_patch                          (:) = nan
     allocate(this%transfer_gr_patch                         (begp:endp)) ; this%transfer_gr_patch                         (:) = nan
     allocate(this%storage_gr_patch                          (begp:endp)) ; this%storage_gr_patch                          (:) = nan
     allocate(this%gr_patch                                  (begp:endp)) ; this%gr_patch                                  (:) = nan
     allocate(this%ar_patch                                  (begp:endp)) ; this%ar_patch                                  (:) = nan
     allocate(this%rr_patch                                  (begp:endp)) ; this%rr_patch                                  (:) = nan
     allocate(this%npp_patch                                 (begp:endp)) ; this%npp_patch                                 (:) = nan
     allocate(this%agnpp_patch                               (begp:endp)) ; this%agnpp_patch                               (:) = nan
     allocate(this%bgnpp_patch                               (begp:endp)) ; this%bgnpp_patch                               (:) = nan
     allocate(this%litfall_patch                             (begp:endp)) ; this%litfall_patch                             (:) = nan
     allocate(this%vegfire_patch                             (begp:endp)) ; this%vegfire_patch                             (:) = nan
     allocate(this%wood_harvestc_patch                       (begp:endp)) ; this%wood_harvestc_patch                       (:) = nan
     allocate(this%cinputs_patch                             (begp:endp)) ; this%cinputs_patch                             (:) = nan
     allocate(this%coutputs_patch                            (begp:endp)) ; this%coutputs_patch                            (:) = nan

     allocate(this%plant_calloc_patch                        (begp:endp)) ; this%plant_calloc_patch                        (:) = nan
     allocate(this%excess_cflux_patch                        (begp:endp)) ; this%excess_cflux_patch                        (:) = nan
     allocate(this%prev_leafc_to_litter_patch                (begp:endp)) ; this%prev_leafc_to_litter_patch                (:) = nan
     allocate(this%prev_frootc_to_litter_patch               (begp:endp)) ; this%prev_frootc_to_litter_patch               (:) = nan
     allocate(this%gpp_patch                                 (begp:endp)) ; this%gpp_patch                                 (:) = nan
     allocate(this%gpp_before_downreg_patch                  (begp:endp)) ; this%gpp_before_downreg_patch                  (:) = nan
     allocate(this%availc_patch                              (begp:endp)) ; this%availc_patch                              (:) = nan
     allocate(this%xsmrpool_recover_patch                    (begp:endp)) ; this%xsmrpool_recover_patch                    (:) = nan
     allocate(this%xsmrpool_c13ratio_patch                   (begp:endp)) ; this%xsmrpool_c13ratio_patch                   (:) = nan

     allocate(this%fire_closs_patch                          (begp:endp)) ; this%fire_closs_patch                          (:) = nan
     allocate(this%cpool_to_grainc_patch                     (begp:endp)) ; this%cpool_to_grainc_patch                     (:) = nan
     allocate(this%cpool_to_grainc_storage_patch             (begp:endp)) ; this%cpool_to_grainc_storage_patch             (:) = nan
     allocate(this%livestemc_to_litter_patch                 (begp:endp)) ; this%livestemc_to_litter_patch                 (:) = nan
     allocate(this%grainc_to_food_patch                      (begp:endp)) ; this%grainc_to_food_patch                      (:) = nan
     allocate(this%grainc_xfer_to_grainc_patch               (begp:endp)) ; this%grainc_xfer_to_grainc_patch               (:) = nan
     allocate(this%cpool_grain_gr_patch                      (begp:endp)) ; this%cpool_grain_gr_patch                      (:) = nan
     allocate(this%cpool_grain_storage_gr_patch              (begp:endp)) ; this%cpool_grain_storage_gr_patch              (:) = nan
     allocate(this%transfer_grain_gr_patch                   (begp:endp)) ; this%transfer_grain_gr_patch                   (:) = nan
     allocate(this%xsmrpool_to_atm_patch                     (begp:endp)) ; this%xsmrpool_to_atm_patch                     (:) = nan
     allocate(this%grainc_storage_to_xfer_patch              (begp:endp)) ; this%grainc_storage_to_xfer_patch              (:) = nan
     allocate(this%frootc_alloc_patch                        (begp:endp)) ; this%frootc_alloc_patch                        (:) = nan
     allocate(this%frootc_loss_patch                         (begp:endp)) ; this%frootc_loss_patch                         (:) = nan
     allocate(this%leafc_alloc_patch                         (begp:endp)) ; this%leafc_alloc_patch                         (:) = nan
     allocate(this%leafc_loss_patch                          (begp:endp)) ; this%leafc_loss_patch                          (:) = nan
     allocate(this%woodc_alloc_patch                         (begp:endp)) ; this%woodc_alloc_patch                         (:) = nan
     allocate(this%woodc_loss_patch                          (begp:endp)) ; this%woodc_loss_patch                          (:) = nan          

     allocate(this%tempavg_agnpp_patch               (begp:endp))                  ; this%tempavg_agnpp_patch (:) = spval
     allocate(this%tempavg_bgnpp_patch               (begp:endp))                  ; this%tempavg_bgnpp_patch (:) = spval
     allocate(this%annavg_agnpp_patch                (begp:endp))                  ; this%annavg_agnpp_patch  (:) = spval ! To detect first year
     allocate(this%annavg_bgnpp_patch                (begp:endp))                  ; this%annavg_bgnpp_patch  (:) = spval ! To detect first year

     allocate(this%t_scalar_col                      (begc:endc,1:nlevdecomp_full)); this%t_scalar_col (:,:)=spval
     allocate(this%w_scalar_col                      (begc:endc,1:nlevdecomp_full)); this%w_scalar_col (:,:)=spval
     allocate(this%o_scalar_col                      (begc:endc,1:nlevdecomp_full)); this%o_scalar_col (:,:)=spval

     allocate(this%phenology_c_to_litr_met_c_col     (begc:endc,1:nlevdecomp_full)); this%phenology_c_to_litr_met_c_col (:,:)=nan
     allocate(this%phenology_c_to_litr_cel_c_col     (begc:endc,1:nlevdecomp_full)); this%phenology_c_to_litr_cel_c_col (:,:)=nan
     allocate(this%phenology_c_to_litr_lig_c_col     (begc:endc,1:nlevdecomp_full)); this%phenology_c_to_litr_lig_c_col (:,:)=nan

     allocate(this%gap_mortality_c_to_litr_met_c_col (begc:endc,1:nlevdecomp_full)); this%gap_mortality_c_to_litr_met_c_col(:,:)=nan
     allocate(this%gap_mortality_c_to_litr_cel_c_col (begc:endc,1:nlevdecomp_full)); this%gap_mortality_c_to_litr_cel_c_col(:,:)=nan
     allocate(this%gap_mortality_c_to_litr_lig_c_col (begc:endc,1:nlevdecomp_full)); this%gap_mortality_c_to_litr_lig_c_col(:,:)=nan

     allocate(this%gap_mortality_c_to_cwdc_col       (begc:endc,1:nlevdecomp_full)); this%gap_mortality_c_to_cwdc_col  (:,:)=nan
     allocate(this%fire_mortality_c_to_cwdc_col      (begc:endc,1:nlevdecomp_full)); this%fire_mortality_c_to_cwdc_col (:,:)=nan
     allocate(this%m_c_to_litr_met_fire_col          (begc:endc,1:nlevdecomp_full)); this%m_c_to_litr_met_fire_col     (:,:)=nan
     allocate(this%m_c_to_litr_cel_fire_col          (begc:endc,1:nlevdecomp_full)); this%m_c_to_litr_cel_fire_col     (:,:)=nan
     allocate(this%m_c_to_litr_lig_fire_col          (begc:endc,1:nlevdecomp_full)); this%m_c_to_litr_lig_fire_col     (:,:)=nan
     allocate(this%harvest_c_to_litr_met_c_col       (begc:endc,1:nlevdecomp_full)); this%harvest_c_to_litr_met_c_col  (:,:)=nan
     allocate(this%harvest_c_to_litr_cel_c_col       (begc:endc,1:nlevdecomp_full)); this%harvest_c_to_litr_cel_c_col  (:,:)=nan
     allocate(this%harvest_c_to_litr_lig_c_col       (begc:endc,1:nlevdecomp_full)); this%harvest_c_to_litr_lig_c_col  (:,:)=nan
     allocate(this%harvest_c_to_cwdc_col             (begc:endc,1:nlevdecomp_full)); this%harvest_c_to_cwdc_col        (:,:)=nan
     allocate(this%phr_vr_col                        (begc:endc,1:nlevdecomp_full)); this%phr_vr_col                   (:,:)=nan 
     allocate(this%fphr_col                          (begc:endc,1:nlevgrnd))       ; this%fphr_col                     (:,:)=nan 

     allocate(this%dwt_frootc_to_litr_met_c_col      (begc:endc,1:nlevdecomp_full)); this%dwt_frootc_to_litr_met_c_col (:,:)=nan
     allocate(this%dwt_frootc_to_litr_cel_c_col      (begc:endc,1:nlevdecomp_full)); this%dwt_frootc_to_litr_cel_c_col (:,:)=nan
     allocate(this%dwt_frootc_to_litr_lig_c_col      (begc:endc,1:nlevdecomp_full)); this%dwt_frootc_to_litr_lig_c_col (:,:)=nan
     allocate(this%dwt_livecrootc_to_cwdc_col        (begc:endc,1:nlevdecomp_full)); this%dwt_livecrootc_to_cwdc_col   (:,:)=nan
     allocate(this%dwt_deadcrootc_to_cwdc_col        (begc:endc,1:nlevdecomp_full)); this%dwt_deadcrootc_to_cwdc_col   (:,:)=nan

     allocate(this%dwt_closs_col                     (begc:endc))                  ; this%dwt_closs_col             (:)  =nan
     allocate(this%dwt_seedc_to_leaf_col             (begc:endc))                  ; this%dwt_seedc_to_leaf_col     (:)  =nan
     allocate(this%dwt_seedc_to_deadstem_col         (begc:endc))                  ; this%dwt_seedc_to_deadstem_col (:)  =nan
     allocate(this%dwt_conv_cflux_col                (begc:endc))                  ; this%dwt_conv_cflux_col        (:)  =nan
     allocate(this%dwt_prod10c_gain_col              (begc:endc))                  ; this%dwt_prod10c_gain_col      (:)  =nan
     allocate(this%dwt_prod100c_gain_col             (begc:endc))                  ; this%dwt_prod100c_gain_col     (:)  =nan
     allocate(this%som_c_leached_col                 (begc:endc))                  ; this%som_c_leached_col         (:)  =nan
     allocate(this%somc_fire_col                     (begc:endc))                  ; this%somc_fire_col             (:)  =nan
     allocate(this%landuseflux_col                   (begc:endc))                  ; this%landuseflux_col           (:)  =nan
     allocate(this%landuptake_col                    (begc:endc))                  ; this%landuptake_col            (:)  =nan
     allocate(this%prod1c_loss_col                   (begc:endc))                  ; this%prod1c_loss_col           (:)  =nan
     allocate(this%prod10c_loss_col                  (begc:endc))                  ; this%prod10c_loss_col          (:)  =nan
     allocate(this%prod100c_loss_col                 (begc:endc))                  ; this%prod100c_loss_col         (:)  =nan
     allocate(this%product_closs_col                 (begc:endc))                  ; this%product_closs_col         (:)  =nan

     allocate(this%bgc_cpool_ext_inputs_vr_col       (begc:endc, 1:nlevdecomp_full,ndecomp_pools));this%bgc_cpool_ext_inputs_vr_col (:,:,:) = nan
     allocate(this%bgc_cpool_ext_loss_vr_col         (begc:endc, 1:nlevdecomp_full,ndecomp_pools));this%bgc_cpool_ext_loss_vr_col   (:,:,:) = nan
     
     allocate(this%lf_conv_cflux_col                 (begc:endc))                  ; this%lf_conv_cflux_col         (:)  =nan
     allocate(this%lithr_col                         (begc:endc))                  ; this%lithr_col                 (:)  =nan
     allocate(this%somhr_col                         (begc:endc))                  ; this%somhr_col                 (:)  =nan
     allocate(this%hr_vr_col                         (begc:endc,1:nlevdecomp_full)); this%hr_vr_col                 (:,:)=nan
     allocate(this%hr_col                            (begc:endc))                  ; this%hr_col                    (:)  =nan
     allocate(this%sr_col                            (begc:endc))                  ; this%sr_col                    (:)  =nan
     allocate(this%er_col                            (begc:endc))                  ; this%er_col                    (:)  =nan
     allocate(this%litfire_col                       (begc:endc))                  ; this%litfire_col               (:)  =nan
     allocate(this%somfire_col                       (begc:endc))                  ; this%somfire_col               (:)  =nan
     allocate(this%totfire_col                       (begc:endc))                  ; this%totfire_col               (:)  =nan
     allocate(this%nep_col                           (begc:endc))                  ; this%nep_col                   (:)  =nan
     allocate(this%nbp_col                           (begc:endc))                  ; this%nbp_col                   (:)  =nan
     allocate(this%nee_col                           (begc:endc))                  ; this%nee_col                   (:)  =nan
     allocate(this%cwdc_hr_col                       (begc:endc))                  ; this%cwdc_hr_col               (:)  =nan
     allocate(this%cwdc_loss_col                     (begc:endc))                  ; this%cwdc_loss_col             (:)  =nan
     allocate(this%litterc_loss_col                  (begc:endc))                  ; this%litterc_loss_col          (:)  =nan
     allocate(this%rr_col                            (begc:endc))                  ; this%rr_col                    (:)  =nan
     allocate(this%ar_col                            (begc:endc))                  ; this%ar_col                    (:)  =nan
     allocate(this%gpp_col                           (begc:endc))                  ; this%gpp_col                   (:)  =nan
     allocate(this%npp_col                           (begc:endc))                  ; this%npp_col                   (:)  =nan
     allocate(this%fire_closs_p2c_col                (begc:endc))                  ; this%fire_closs_p2c_col        (:)  =nan
     allocate(this%fire_closs_col                    (begc:endc))                  ; this%fire_closs_col            (:)  =nan
     allocate(this%litfall_col                       (begc:endc))                  ; this%litfall_col               (:)  =nan
     allocate(this%vegfire_col                       (begc:endc))                  ; this%vegfire_col               (:)  =nan
     allocate(this%wood_harvestc_col                 (begc:endc))                  ; this%wood_harvestc_col         (:)  =nan
     allocate(this%hrv_xsmrpool_to_atm_col           (begc:endc))                  ; this%hrv_xsmrpool_to_atm_col   (:)  =nan 

     allocate(this%hrv_deadstemc_to_prod10c_col(begc:endc))                                                    
     this%hrv_deadstemc_to_prod10c_col(:)= nan

     allocate(this%hrv_deadstemc_to_prod100c_col(begc:endc))                                                   
     this%hrv_deadstemc_to_prod100c_col(:)= nan
 
     allocate(this%hrv_cropc_to_prod1c_col(begc:endc))
     this%hrv_cropc_to_prod1c_col(:) = nan

     allocate(this%m_decomp_cpools_to_fire_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))                
     this%m_decomp_cpools_to_fire_vr_col(:,:,:)= nan

     allocate(this%m_decomp_cpools_to_fire_col(begc:endc,1:ndecomp_pools))                                     
     this%m_decomp_cpools_to_fire_col(:,:)= nan

     allocate(this%decomp_cpools_sourcesink_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))                  
     this%decomp_cpools_sourcesink_col(:,:,:)= nan

     allocate(this%decomp_cascade_hr_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))        
     this%decomp_cascade_hr_vr_col(:,:,:)= spval

     allocate(this%decomp_cascade_hr_col(begc:endc,1:ndecomp_cascade_transitions))                             
     this%decomp_cascade_hr_col(:,:)= nan

     allocate(this%decomp_cascade_ctransfer_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)) 
     this%decomp_cascade_ctransfer_vr_col(:,:,:)= nan

     allocate(this%decomp_cascade_ctransfer_col(begc:endc,1:ndecomp_cascade_transitions))                      
     this%decomp_cascade_ctransfer_col(:,:)= nan

     allocate(this%decomp_k_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))                    
     this%decomp_k_col(:,:,:)= spval

     allocate(this%decomp_cpools_leached_col(begc:endc,1:ndecomp_pools))              
     this%decomp_cpools_leached_col(:,:)= nan

     allocate(this%decomp_cpools_transport_tendency_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))          
     this%decomp_cpools_transport_tendency_col(:,:,:)= nan


     allocate(this%tempsum_npp_patch     (begp:endp)) ; this%tempsum_npp_patch     (:) = nan
     allocate(this%annsum_npp_patch      (begp:endp)) ; this%annsum_npp_patch      (:) = nan
     allocate(this%tempsum_litfall_patch (begp:endp)) ; this%tempsum_litfall_patch (:) = nan
     allocate(this%annsum_litfall_patch  (begp:endp)) ; this%annsum_litfall_patch  (:) = nan
     allocate(this%annsum_npp_col        (begc:endc)) ; this%annsum_npp_col        (:) = nan
     allocate(this%lag_npp_col           (begc:endc)) ; this%lag_npp_col           (:) = spval
     
     ! debug
     allocate(this%plant_to_litter_cflux (begc:endc)) ;	this%plant_to_litter_cflux (:) = nan
     allocate(this%plant_to_cwd_cflux    (begc:endc)) ;	this%plant_to_cwd_cflux	   (:) = nan
     allocate(this%allocation_leaf       (begp:endp)) ; this%allocation_leaf       (:) = nan
     allocate(this%allocation_stem       (begp:endp)) ; this%allocation_stem       (:) = nan
     allocate(this%allocation_froot      (begp:endp)) ; this%allocation_froot      (:) = nan

     ! clm_bgc_interface & pflotran
     !------------------------------------------------------------------------
     allocate(this%externalc_to_decomp_cpools_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools)); this%externalc_to_decomp_cpools_col(:,:,:) = spval
     allocate(this%externalc_to_decomp_delta_col (begc:endc));                                   this%externalc_to_decomp_delta_col (:)     = spval
     allocate(this%f_co2_soil_vr_col             (begc:endc,1:nlevdecomp_full));                 this%f_co2_soil_vr_col             (:,:)   = nan
     allocate(this%f_co2_soil_col                (begc:endc))                  ;                 this%f_co2_soil_col                (:)     = nan
     !------------------------------------------------------------------------
   end subroutine InitAllocate; 

   !------------------------------------------------------------------------
   subroutine InitHistory(this, bounds, carbon_type)
     !
     ! !DESCRIPTION:
     ! add history fields for all CN variables, always set as default='inactive'
     !
     ! !USES:
     use clm_varpar , only : ndecomp_cascade_transitions, ndecomp_pools
     use clm_varpar , only : nlevdecomp, nlevdecomp_full, crop_prog, nlevgrnd
     use clm_varctl , only : hist_wrtch4diag
     use histFileMod, only : hist_addfld1d, hist_addfld2d, hist_addfld_decomp 
     use tracer_varcon    , only : is_active_betr_bgc
     use clm_varctl,  only : get_carbontag
    !
     ! !ARGUMENTS:
     class(carbonflux_type) :: this    
     type(bounds_type)         , intent(in) :: bounds 
     character(len=3)          , intent(in) :: carbon_type ! one of ['c12', c13','c14']
     !
     ! !LOCAL VARIABLES:
     integer           :: k,l,ii,jj 
     character(8)      :: vr_suffix
     character(10)     :: active
     integer           :: begp,endp
     integer           :: begc,endc
     character(24)     :: fieldname
     character(100)    :: longname
     real(r8), pointer :: data1dptr(:)   ! temp. pointer for slicing larger arrays
     real(r8), pointer :: data2dptr(:,:) ! temp. pointer for slicing larger arrays
     character(len=3)  :: ctag
     !---------------------------------------------------------------------

     begp = bounds%begp; endp = bounds%endp
     begc = bounds%begc; endc = bounds%endc

     if (nlevdecomp > 1) then
        vr_suffix = "_vr"
     else 
        vr_suffix = ""
     endif

     !-------------------------------
     ! C flux variables - native to PFT
     !-------------------------------

     ! add history fields for all CLAMP CN variables

     if (carbon_type == 'c12') then
        if (crop_prog) then
           this%grainc_to_food_patch(begp:endp) = spval
           call hist_addfld1d (fname='GRAINC_TO_FOOD', units='gC/m^2/s', &
                avgflag='A', long_name='grain C to food', &
                ptr_patch=this%grainc_to_food_patch, default='inactive')
        end if

        this%woodc_alloc_patch(begp:endp) = spval
        call hist_addfld1d (fname='WOODC_ALLOC', units='gC/m^2/s', &
             avgflag='A', long_name='wood C eallocation', &
             ptr_patch=this%woodc_alloc_patch)

        this%woodc_loss_patch(begp:endp) = spval
        call hist_addfld1d (fname='WOODC_LOSS', units='gC/m^2/s', &
             avgflag='A', long_name='wood C loss', &
             ptr_patch=this%woodc_loss_patch)

        this%leafc_loss_patch(begp:endp) = spval
        call hist_addfld1d (fname='LEAFC_LOSS', units='gC/m^2/s', &
             avgflag='A', long_name='leaf C loss', &
             ptr_patch=this%leafc_loss_patch)

        this%leafc_alloc_patch(begp:endp) = spval
        call hist_addfld1d (fname='LEAFC_ALLOC', units='gC/m^2/s', &
             avgflag='A', long_name='leaf C allocation', &
             ptr_patch=this%leafc_alloc_patch)

        this%frootc_loss_patch(begp:endp) = spval
        call hist_addfld1d (fname='FROOTC_LOSS', units='gC/m^2/s', &
             avgflag='A', long_name='fine root C loss', &
             ptr_patch=this%frootc_loss_patch)

        this%frootc_alloc_patch(begp:endp) = spval
        call hist_addfld1d (fname='FROOTC_ALLOC', units='gC/m^2/s', &
             avgflag='A', long_name='fine root C allocation', &
             ptr_patch=this%frootc_alloc_patch)

        this%m_leafc_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_LEAFC_TO_LITTER', units='gC/m^2/s', &
             avgflag='A', long_name='leaf C mortality', &
             ptr_patch=this%m_leafc_to_litter_patch, default='inactive')

        this%m_frootc_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_FROOTC_TO_LITTER', units='gC/m^2/s', &
             avgflag='A', long_name='fine root C mortality', &
             ptr_patch=this%m_frootc_to_litter_patch, default='inactive')

        this%m_leafc_storage_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_LEAFC_STORAGE_TO_LITTER', units='gC/m^2/s', &
             avgflag='A', long_name='leaf C storage mortality', &
             ptr_patch=this%m_leafc_storage_to_litter_patch, default='inactive')

        this%m_frootc_storage_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_FROOTC_STORAGE_TO_LITTER', units='gC/m^2/s', &
             avgflag='A', long_name='fine root C storage mortality', &
             ptr_patch=this%m_frootc_storage_to_litter_patch, default='inactive')

        this%m_livestemc_storage_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_LIVESTEMC_STORAGE_TO_LITTER', units='gC/m^2/s', &
             avgflag='A', long_name='live stem C storage mortality', &
             ptr_patch=this%m_livestemc_storage_to_litter_patch, default='inactive')

        this%m_deadstemc_storage_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_DEADSTEMC_STORAGE_TO_LITTER', units='gC/m^2/s', &
             avgflag='A', long_name='dead stem C storage mortality', &
             ptr_patch=this%m_deadstemc_storage_to_litter_patch, default='inactive')

        this%m_livecrootc_storage_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_LIVECROOTC_STORAGE_TO_LITTER', units='gC/m^2/s', &
             avgflag='A', long_name='live coarse root C storage mortality', &
             ptr_patch=this%m_livecrootc_storage_to_litter_patch, default='inactive')

        this%m_deadcrootc_storage_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_DEADCROOTC_STORAGE_TO_LITTER', units='gC/m^2/s', &
             avgflag='A', long_name='dead coarse root C storage mortality', &
             ptr_patch=this%m_deadcrootc_storage_to_litter_patch, default='inactive')

        this%m_leafc_xfer_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_LEAFC_XFER_TO_LITTER', units='gC/m^2/s', &
             avgflag='A', long_name='leaf C transfer mortality', &
             ptr_patch=this%m_leafc_xfer_to_litter_patch, default='inactive')

        this%m_frootc_xfer_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_FROOTC_XFER_TO_LITTER', units='gC/m^2/s', &
             avgflag='A', long_name='fine root C transfer mortality', &
             ptr_patch=this%m_frootc_xfer_to_litter_patch, default='inactive')

        this%m_livestemc_xfer_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_LIVESTEMC_XFER_TO_LITTER', units='gC/m^2/s', &
             avgflag='A', long_name='live stem C transfer mortality', &
             ptr_patch=this%m_livestemc_xfer_to_litter_patch, default='inactive')

        this%m_deadstemc_xfer_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_DEADSTEMC_XFER_TO_LITTER', units='gC/m^2/s', &
             avgflag='A', long_name='dead stem C transfer mortality', &
             ptr_patch=this%m_deadstemc_xfer_to_litter_patch, default='inactive')

        this%m_livecrootc_xfer_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_LIVECROOTC_XFER_TO_LITTER', units='gC/m^2/s', &
             avgflag='A', long_name='live coarse root C transfer mortality', &
             ptr_patch=this%m_livecrootc_xfer_to_litter_patch, default='inactive')

        this%m_deadcrootc_xfer_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_DEADCROOTC_XFER_TO_LITTER', units='gC/m^2/s', &
             avgflag='A', long_name='dead coarse root C transfer mortality', &
             ptr_patch=this%m_deadcrootc_xfer_to_litter_patch, default='inactive')

        this%m_livestemc_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_LIVESTEMC_TO_LITTER', units='gC/m^2/s', &
             avgflag='A', long_name='live stem C mortality', &
             ptr_patch=this%m_livestemc_to_litter_patch, default='inactive')

        this%m_deadstemc_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_DEADSTEMC_TO_LITTER', units='gC/m^2/s', &
             avgflag='A', long_name='dead stem C mortality', &
             ptr_patch=this%m_deadstemc_to_litter_patch, default='inactive')

        this%m_livecrootc_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_LIVECROOTC_TO_LITTER', units='gC/m^2/s', &
             avgflag='A', long_name='live coarse root C mortality', &
             ptr_patch=this%m_livecrootc_to_litter_patch, default='inactive')

        this%m_deadcrootc_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_DEADCROOTC_TO_LITTER', units='gC/m^2/s', &
             avgflag='A', long_name='dead coarse root C mortality', &
             ptr_patch=this%m_deadcrootc_to_litter_patch, default='inactive')

        this%m_gresp_storage_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_GRESP_STORAGE_TO_LITTER', units='gC/m^2/s', &
             avgflag='A', long_name='growth respiration storage mortality', &
             ptr_patch=this%m_gresp_storage_to_litter_patch, default='inactive')

        this%m_gresp_xfer_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_GRESP_XFER_TO_LITTER', units='gC/m^2/s', &
             avgflag='A', long_name='growth respiration transfer mortality', &
             ptr_patch=this%m_gresp_xfer_to_litter_patch, default='inactive')

        this%m_leafc_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_LEAFC_TO_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='leaf C fire loss', &
             ptr_patch=this%m_leafc_to_fire_patch, default='inactive')

        this%m_leafc_storage_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_LEAFC_STORAGE_TO_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='leaf C storage fire loss', &
             ptr_patch=this%m_leafc_storage_to_fire_patch, default='inactive')

        this%m_leafc_xfer_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_LEAFC_XFER_TO_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='leaf C transfer fire loss', &
             ptr_patch=this%m_leafc_xfer_to_fire_patch, default='inactive')

        this%m_livestemc_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_LIVESTEMC_TO_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='live stem C fire loss', &
             ptr_patch=this%m_livestemc_to_fire_patch, default='inactive')

        this%m_livestemc_storage_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_LIVESTEMC_STORAGE_TO_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='live stem C storage fire loss', &
             ptr_patch=this%m_livestemc_storage_to_fire_patch, default='inactive')

        this%m_livestemc_xfer_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_LIVESTEMC_XFER_TO_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='live stem C transfer fire loss', &
             ptr_patch=this%m_livestemc_xfer_to_fire_patch, default='inactive')

        this%m_deadstemc_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_DEADSTEMC_TO_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='dead stem C fire loss', &
             ptr_patch=this%m_deadstemc_to_fire_patch, default='inactive')

        this%m_deadstemc_storage_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_DEADSTEMC_STORAGE_TO_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='dead stem C storage fire loss', &
             ptr_patch=this%m_deadstemc_storage_to_fire_patch, default='inactive')

        this%m_deadstemc_xfer_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_DEADSTEMC_XFER_TO_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='dead stem C transfer fire loss', &
             ptr_patch=this%m_deadstemc_xfer_to_fire_patch, default='inactive')

        this%m_frootc_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_FROOTC_TO_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='fine root C fire loss', &
             ptr_patch=this%m_frootc_to_fire_patch, default='inactive')

        this%m_frootc_storage_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_FROOTC_STORAGE_TO_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='fine root C storage fire loss', &
             ptr_patch=this%m_frootc_storage_to_fire_patch, default='inactive')

        this%m_frootc_xfer_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_FROOTC_XFER_TO_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='fine root C transfer fire loss', &
             ptr_patch=this%m_frootc_xfer_to_fire_patch, default='inactive')

        this%m_livecrootc_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_LIVEROOTC_TO_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='live root C fire loss', &
             ptr_patch=this%m_livecrootc_to_fire_patch, default='inactive')

        this%m_livecrootc_storage_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_LIVEROOTC_STORAGE_TO_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='live root C storage fire loss', &
             ptr_patch=this%m_livecrootc_storage_to_fire_patch, default='inactive')

        this%m_livecrootc_xfer_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_LIVEROOTC_XFER_TO_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='live root C transfer fire loss', &
             ptr_patch=this%m_livecrootc_xfer_to_fire_patch, default='inactive')

        this%m_deadcrootc_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_DEADROOTC_TO_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='dead root C fire loss', &
             ptr_patch=this%m_deadcrootc_to_fire_patch, default='inactive')

        this%m_deadcrootc_storage_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_DEADROOTC_STORAGE_TO_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='dead root C storage fire loss', &
             ptr_patch=this%m_deadcrootc_storage_to_fire_patch, default='inactive')

        this%m_deadcrootc_xfer_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_DEADROOTC_XFER_TO_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='dead root C transfer fire loss', &
             ptr_patch=this%m_deadcrootc_xfer_to_fire_patch, default='inactive')

        this%m_gresp_storage_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_GRESP_STORAGE_TO_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='growth respiration storage fire loss', &
             ptr_patch=this%m_gresp_storage_to_fire_patch, default='inactive')

        this%m_gresp_xfer_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_GRESP_XFER_TO_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='growth respiration transfer fire loss', &
             ptr_patch=this%m_gresp_xfer_to_fire_patch, default='inactive')

        this%m_leafc_to_litter_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_LEAFC_TO_LITTER_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='leaf C fire mortality to litter', &
             ptr_patch=this%m_leafc_to_litter_fire_patch, default='inactive')

        ! add by F. Li and S. Levis
        this%m_leafc_storage_to_litter_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_LEAFC_STORAGE_TO_LITTER_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='leaf C fire mortality to litter', &
             ptr_patch=this%m_leafc_storage_to_litter_fire_patch, default='inactive')

        this%m_leafc_xfer_to_litter_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_LEAFC_XFER_TO_LITTER_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='leaf C transfer fire mortality to litter', &
             ptr_patch=this%m_leafc_xfer_to_litter_fire_patch, default='inactive')

        this%m_livestemc_to_litter_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_LIVESTEMC_TO_LITTER_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='live stem C fire mortality to litter', &
             ptr_patch=this%m_livestemc_to_litter_fire_patch, default='inactive')

        this%m_livestemc_storage_to_litter_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_LIVESTEMC_STORAGE_TO_LITTER_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='live stem C storage fire mortality to litter', &
             ptr_patch=this%m_livestemc_storage_to_litter_fire_patch, default='inactive')

        this%m_livestemc_xfer_to_litter_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_LIVESTEMC_XFER_TO_LITTER_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='live stem C transfer fire mortality to litter', &
             ptr_patch=this%m_livestemc_xfer_to_litter_fire_patch, default='inactive')

        this%m_livestemc_to_deadstemc_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_LIVESTEMC_TO_DEADSTEMC_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='live stem C fire mortality to dead stem C', &
             ptr_patch=this%m_livestemc_to_deadstemc_fire_patch, default='inactive')

        this%m_deadstemc_to_litter_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_DEADSTEMC_TO_LITTER_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='dead stem C fire mortality to litter', &
             ptr_patch=this%m_deadstemc_to_litter_fire_patch, default='inactive')

        this%m_deadstemc_storage_to_litter_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_DEADSTEMC_STORAGE_TO_LITTER_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='dead stem C storage fire mortality to litter', &
             ptr_patch=this%m_deadstemc_storage_to_litter_fire_patch, default='inactive')

        this%m_deadstemc_xfer_to_litter_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_DEADSTEMC_XFER_TO_LITTER_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='dead stem C transfer fire mortality to litter', &
             ptr_patch=this%m_deadstemc_xfer_to_litter_fire_patch, default='inactive')

        this%m_frootc_to_litter_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_FROOTC_TO_LITTER_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='fine root C fire mortality to litter', &
             ptr_patch=this%m_frootc_to_litter_fire_patch, default='inactive')

        this%m_frootc_storage_to_litter_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_FROOTC_STORAGE_TO_LITTER_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='fine root C storage fire mortality to litter', &
             ptr_patch=this%m_frootc_storage_to_litter_fire_patch, default='inactive')

        this%m_frootc_xfer_to_litter_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_FROOTC_XFER_TO_LITTER_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='fine root C transfer fire mortality to litter', &
             ptr_patch=this%m_frootc_xfer_to_litter_fire_patch, default='inactive')

        this%m_livecrootc_to_litter_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_LIVEROOTC_TO_LITTER_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='live root C fire mortality to litter', &
             ptr_patch=this%m_livecrootc_to_litter_fire_patch, default='inactive')

        this%m_livecrootc_storage_to_litter_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_LIVEROOTC_STORAGE_TO_LITTER_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='live root C storage fire mortality to litter', &
             ptr_patch=this%m_livecrootc_storage_to_litter_fire_patch, default='inactive')

        this%m_livecrootc_xfer_to_litter_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_LIVEROOTC_XFER_TO_LITTER_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='live root C transfer fire mortality to litter', &
             ptr_patch=this%m_livecrootc_xfer_to_litter_fire_patch, default='inactive')

        this%m_livecrootc_to_deadcrootc_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_LIVEROOTC_TO_DEADROOTC_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='live root C fire mortality to dead root C', &
             ptr_patch=this%m_livecrootc_to_deadcrootc_fire_patch, default='inactive')


        this%m_deadcrootc_to_litter_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_DEADROOTC_TO_LITTER_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='dead root C fire mortality to litter', &
             ptr_patch=this%m_deadcrootc_to_litter_fire_patch, default='inactive')

        this%m_deadcrootc_storage_to_litter_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_DEADROOTC_STORAGE_TO_LITTER_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='dead root C storage fire mortality to litter', &
             ptr_patch=this%m_deadcrootc_storage_to_litter_fire_patch, default='inactive')

        this%m_deadcrootc_xfer_to_litter_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_DEADROOTC_XFER_TO_LITTER_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='dead root C transfer fire mortality to litter', &
             ptr_patch=this%m_deadcrootc_xfer_to_litter_fire_patch, default='inactive')

        this%m_livecrootc_storage_to_litter_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_LIVECROOTC_STORAGE_TO_LITTER_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='live coarse root C fire mortality to litter', &
             ptr_patch=this%m_livecrootc_storage_to_litter_fire_patch, default='inactive')

        this%m_deadcrootc_storage_to_litter_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_DEADCROOTC_STORAGE_TO_LITTER_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='dead coarse root C storage fire mortality to litter', &
             ptr_patch=this%m_deadcrootc_storage_to_litter_fire_patch,  default='inactive')

        this%m_gresp_storage_to_litter_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_GRESP_STORAGE_TO_LITTER_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='growth respiration storage fire mortality to litter', &
             ptr_patch=this%m_gresp_storage_to_litter_fire_patch, default='inactive')

        this%m_gresp_xfer_to_litter_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='M_GRESP_XFER_TO_LITTER_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='growth respiration transfer fire mortality to litter', &
             ptr_patch=this%m_gresp_xfer_to_litter_fire_patch, default='inactive')   

        this%leafc_xfer_to_leafc_patch(begp:endp) = spval
        call hist_addfld1d (fname='LEAFC_XFER_TO_LEAFC', units='gC/m^2/s', &
             avgflag='A', long_name='leaf C growth from storage', &
             ptr_patch=this%leafc_xfer_to_leafc_patch, default='inactive')

        this%frootc_xfer_to_frootc_patch(begp:endp) = spval
        call hist_addfld1d (fname='FROOTC_XFER_TO_FROOTC', units='gC/m^2/s', &
             avgflag='A', long_name='fine root C growth from storage', &
             ptr_patch=this%frootc_xfer_to_frootc_patch, default='inactive')

        this%livestemc_xfer_to_livestemc_patch(begp:endp) = spval
        call hist_addfld1d (fname='LIVESTEMC_XFER_TO_LIVESTEMC', units='gC/m^2/s', &
             avgflag='A', long_name='live stem C growth from storage', &
             ptr_patch=this%livestemc_xfer_to_livestemc_patch, default='inactive')

        this%deadstemc_xfer_to_deadstemc_patch(begp:endp) = spval
        call hist_addfld1d (fname='DEADSTEMC_XFER_TO_DEADSTEMC', units='gC/m^2/s', &
             avgflag='A', long_name='dead stem C growth from storage', &
             ptr_patch=this%deadstemc_xfer_to_deadstemc_patch, default='inactive')

        this%livecrootc_xfer_to_livecrootc_patch(begp:endp) = spval
        call hist_addfld1d (fname='LIVECROOTC_XFER_TO_LIVECROOTC', units='gC/m^2/s', &
             avgflag='A', long_name='live coarse root C growth from storage', &
             ptr_patch=this%livecrootc_xfer_to_livecrootc_patch, default='inactive')

        this%deadcrootc_xfer_to_deadcrootc_patch(begp:endp) = spval
        call hist_addfld1d (fname='DEADCROOTC_XFER_TO_DEADCROOTC', units='gC/m^2/s', &
             avgflag='A', long_name='dead coarse root C growth from storage', &
             ptr_patch=this%deadcrootc_xfer_to_deadcrootc_patch, default='inactive')

        this%leafc_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='LEAFC_TO_LITTER', units='gC/m^2/s', &
             avgflag='A', long_name='leaf C litterfall', &
             ptr_patch=this%leafc_to_litter_patch, default='active')

        this%frootc_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='FROOTC_TO_LITTER', units='gC/m^2/s', &
             avgflag='A', long_name='fine root C litterfall', &
             ptr_patch=this%frootc_to_litter_patch, default='inactive')

        this%leaf_mr_patch(begp:endp) = spval
        call hist_addfld1d (fname='LEAF_MR', units='gC/m^2/s', &
             avgflag='A', long_name='leaf maintenance respiration', &
             ptr_patch=this%leaf_mr_patch)

        this%froot_mr_patch(begp:endp) = spval
        call hist_addfld1d (fname='FROOT_MR', units='gC/m^2/s', &
             avgflag='A', long_name='fine root maintenance respiration', &
             ptr_patch=this%froot_mr_patch, default='inactive')

        this%livestem_mr_patch(begp:endp) = spval
        call hist_addfld1d (fname='LIVESTEM_MR', units='gC/m^2/s', &
             avgflag='A', long_name='live stem maintenance respiration', &
             ptr_patch=this%livestem_mr_patch, default='inactive')

        this%livecroot_mr_patch(begp:endp) = spval
        call hist_addfld1d (fname='LIVECROOT_MR', units='gC/m^2/s', &
             avgflag='A', long_name='live coarse root maintenance respiration', &
             ptr_patch=this%livecroot_mr_patch, default='inactive')

        this%psnsun_to_cpool_patch(begp:endp) = spval
        call hist_addfld1d (fname='PSNSUN_TO_CPOOL', units='gC/m^2/s', &
             avgflag='A', long_name='C fixation from sunlit canopy', &
             ptr_patch=this%psnsun_to_cpool_patch)

        this%psnshade_to_cpool_patch(begp:endp) = spval
        call hist_addfld1d (fname='PSNSHADE_TO_CPOOL', units='gC/m^2/s', &
             avgflag='A', long_name='C fixation from shaded canopy', &
             ptr_patch=this%psnshade_to_cpool_patch)

        this%cpool_to_leafc_patch(begp:endp) = spval
        call hist_addfld1d (fname='CPOOL_TO_LEAFC', units='gC/m^2/s', &
             avgflag='A', long_name='allocation to leaf C', &
             ptr_patch=this%cpool_to_leafc_patch, default='inactive')

        this%cpool_to_leafc_storage_patch(begp:endp) = spval
        call hist_addfld1d (fname='CPOOL_TO_LEAFC_STORAGE', units='gC/m^2/s', &
             avgflag='A', long_name='allocation to leaf C storage', &
             ptr_patch=this%cpool_to_leafc_storage_patch, default='inactive')

        this%cpool_to_frootc_patch(begp:endp) = spval
        call hist_addfld1d (fname='CPOOL_TO_FROOTC', units='gC/m^2/s', &
             avgflag='A', long_name='allocation to fine root C', &
             ptr_patch=this%cpool_to_frootc_patch, default='inactive')

        this%cpool_to_frootc_storage_patch(begp:endp) = spval
        call hist_addfld1d (fname='CPOOL_TO_FROOTC_STORAGE', units='gC/m^2/s', &
             avgflag='A', long_name='allocation to fine root C storage', &
             ptr_patch=this%cpool_to_frootc_storage_patch, default='inactive')

        this%cpool_to_livestemc_patch(begp:endp) = spval
        call hist_addfld1d (fname='CPOOL_TO_LIVESTEMC', units='gC/m^2/s', &
             avgflag='A', long_name='allocation to live stem C', &
             ptr_patch=this%cpool_to_livestemc_patch, default='inactive')

        this%cpool_to_livestemc_storage_patch(begp:endp) = spval
        call hist_addfld1d (fname='CPOOL_TO_LIVESTEMC_STORAGE', units='gC/m^2/s', &
             avgflag='A', long_name='allocation to live stem C storage', &
             ptr_patch=this%cpool_to_livestemc_storage_patch, default='inactive')

        this%cpool_to_deadstemc_patch(begp:endp) = spval
        call hist_addfld1d (fname='CPOOL_TO_DEADSTEMC', units='gC/m^2/s', &
             avgflag='A', long_name='allocation to dead stem C', &
             ptr_patch=this%cpool_to_deadstemc_patch, default='inactive')

        this%cpool_to_deadstemc_storage_patch(begp:endp) = spval
        call hist_addfld1d (fname='CPOOL_TO_DEADSTEMC_STORAGE', units='gC/m^2/s', &
             avgflag='A', long_name='allocation to dead stem C storage', &
             ptr_patch=this%cpool_to_deadstemc_storage_patch, default='inactive')

        this%cpool_to_livecrootc_patch(begp:endp) = spval
        call hist_addfld1d (fname='CPOOL_TO_LIVECROOTC', units='gC/m^2/s', &
             avgflag='A', long_name='allocation to live coarse root C', &
             ptr_patch=this%cpool_to_livecrootc_patch, default='inactive')

        this%cpool_to_livecrootc_storage_patch(begp:endp) = spval
        call hist_addfld1d (fname='CPOOL_TO_LIVECROOTC_STORAGE', units='gC/m^2/s', &
             avgflag='A', long_name='allocation to live coarse root C storage', &
             ptr_patch=this%cpool_to_livecrootc_storage_patch, default='inactive')

        this%cpool_to_deadcrootc_patch(begp:endp) = spval
        call hist_addfld1d (fname='CPOOL_TO_DEADCROOTC', units='gC/m^2/s', &
             avgflag='A', long_name='allocation to dead coarse root C', &
             ptr_patch=this%cpool_to_deadcrootc_patch, default='inactive')

        this%cpool_to_deadcrootc_storage_patch(begp:endp) = spval
        call hist_addfld1d (fname='CPOOL_TO_DEADCROOTC_STORAGE', units='gC/m^2/s', &
             avgflag='A', long_name='allocation to dead coarse root C storage', &
             ptr_patch=this%cpool_to_deadcrootc_storage_patch, default='inactive')

        this%cpool_to_gresp_storage_patch(begp:endp) = spval
        call hist_addfld1d (fname='CPOOL_TO_GRESP_STORAGE', units='gC/m^2/s', &
             avgflag='A', long_name='allocation to growth respiration storage', &
             ptr_patch=this%cpool_to_gresp_storage_patch, default='inactive')

        this%cpool_leaf_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='CPOOL_LEAF_GR', units='gC/m^2/s', &
             avgflag='A', long_name='leaf growth respiration', &
             ptr_patch=this%cpool_leaf_gr_patch, default='inactive')

        this%cpool_leaf_storage_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='CPOOL_LEAF_STORAGE_GR', units='gC/m^2/s', &
             avgflag='A', long_name='leaf growth respiration to storage', &
             ptr_patch=this%cpool_leaf_storage_gr_patch, default='inactive')

        this%transfer_leaf_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='TRANSFER_LEAF_GR', units='gC/m^2/s', &
             avgflag='A', long_name='leaf growth respiration from storage', &
             ptr_patch=this%transfer_leaf_gr_patch, default='inactive')

        this%cpool_froot_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='CPOOL_FROOT_GR', units='gC/m^2/s', &
             avgflag='A', long_name='fine root growth respiration', &
             ptr_patch=this%cpool_froot_gr_patch, default='inactive')

        this%cpool_froot_storage_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='CPOOL_FROOT_STORAGE_GR', units='gC/m^2/s', &
             avgflag='A', long_name='fine root  growth respiration to storage', &
             ptr_patch=this%cpool_froot_storage_gr_patch, default='inactive')

        this%transfer_froot_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='TRANSFER_FROOT_GR', units='gC/m^2/s', &
             avgflag='A', long_name='fine root  growth respiration from storage', &
             ptr_patch=this%transfer_froot_gr_patch, default='inactive')

        this%cpool_livestem_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='CPOOL_LIVESTEM_GR', units='gC/m^2/s', &
             avgflag='A', long_name='live stem growth respiration', &
             ptr_patch=this%cpool_livestem_gr_patch, default='inactive')

        this%cpool_livestem_storage_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='CPOOL_LIVESTEM_STORAGE_GR', units='gC/m^2/s', &
             avgflag='A', long_name='live stem growth respiration to storage', &
             ptr_patch=this%cpool_livestem_storage_gr_patch, default='inactive')

        this%transfer_livestem_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='TRANSFER_LIVESTEM_GR', units='gC/m^2/s', &
             avgflag='A', long_name='live stem growth respiration from storage', &
             ptr_patch=this%transfer_livestem_gr_patch, default='inactive')

        this%cpool_deadstem_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='CPOOL_DEADSTEM_GR', units='gC/m^2/s', &
             avgflag='A', long_name='dead stem growth respiration', &
             ptr_patch=this%cpool_deadstem_gr_patch, default='inactive')

        this%cpool_deadstem_storage_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='CPOOL_DEADSTEM_STORAGE_GR', units='gC/m^2/s', &
             avgflag='A', long_name='dead stem growth respiration to storage', &
             ptr_patch=this%cpool_deadstem_storage_gr_patch, default='inactive')

        this%transfer_deadstem_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='TRANSFER_DEADSTEM_GR', units='gC/m^2/s', &
             avgflag='A', long_name='dead stem growth respiration from storage', &
             ptr_patch=this%transfer_deadstem_gr_patch, default='inactive')

        this%cpool_livecroot_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='CPOOL_LIVECROOT_GR', units='gC/m^2/s', &
             avgflag='A', long_name='live coarse root growth respiration', &
             ptr_patch=this%cpool_livecroot_gr_patch, default='inactive')

        this%cpool_livecroot_storage_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='CPOOL_LIVECROOT_STORAGE_GR', units='gC/m^2/s', &
             avgflag='A', long_name='live coarse root growth respiration to storage', &
             ptr_patch=this%cpool_livecroot_storage_gr_patch, default='inactive')

        this%transfer_livecroot_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='TRANSFER_LIVECROOT_GR', units='gC/m^2/s', &
             avgflag='A', long_name='live coarse root growth respiration from storage', &
             ptr_patch=this%transfer_livecroot_gr_patch, default='inactive')

        this%cpool_deadcroot_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='CPOOL_DEADCROOT_GR', units='gC/m^2/s', &
             avgflag='A', long_name='dead coarse root growth respiration', &
             ptr_patch=this%cpool_deadcroot_gr_patch, default='inactive')

        this%cpool_deadcroot_storage_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='CPOOL_DEADCROOT_STORAGE_GR', units='gC/m^2/s', &
             avgflag='A', long_name='dead coarse root growth respiration to storage', &
             ptr_patch=this%cpool_deadcroot_storage_gr_patch, default='inactive')

        this%transfer_deadcroot_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='TRANSFER_DEADCROOT_GR', units='gC/m^2/s', &
             avgflag='A', long_name='dead coarse root growth respiration from storage', &
             ptr_patch=this%transfer_deadcroot_gr_patch, default='inactive')

        this%leafc_storage_to_xfer_patch(begp:endp) = spval
        call hist_addfld1d (fname='LEAFC_STORAGE_TO_XFER', units='gC/m^2/s', &
             avgflag='A', long_name='leaf C shift storage to transfer', &
             ptr_patch=this%leafc_storage_to_xfer_patch, default='inactive')

        this%frootc_storage_to_xfer_patch(begp:endp) = spval
        call hist_addfld1d (fname='FROOTC_STORAGE_TO_XFER', units='gC/m^2/s', &
             avgflag='A', long_name='fine root C shift storage to transfer', &
             ptr_patch=this%frootc_storage_to_xfer_patch, default='inactive')

        this%livestemc_storage_to_xfer_patch(begp:endp) = spval
        call hist_addfld1d (fname='LIVESTEMC_STORAGE_TO_XFER', units='gC/m^2/s', &
             avgflag='A', long_name='live stem C shift storage to transfer', &
             ptr_patch=this%livestemc_storage_to_xfer_patch, default='inactive')

        this%deadstemc_storage_to_xfer_patch(begp:endp) = spval
        call hist_addfld1d (fname='DEADSTEMC_STORAGE_TO_XFER', units='gC/m^2/s', &
             avgflag='A', long_name='dead stem C shift storage to transfer', &
             ptr_patch=this%deadstemc_storage_to_xfer_patch, default='inactive')

        this%livecrootc_storage_to_xfer_patch(begp:endp) = spval
        call hist_addfld1d (fname='LIVECROOTC_STORAGE_TO_XFER', units='gC/m^2/s', &
             avgflag='A', long_name='live coarse root C shift storage to transfer', &
             ptr_patch=this%livecrootc_storage_to_xfer_patch, default='inactive')

        this%deadcrootc_storage_to_xfer_patch(begp:endp) = spval
        call hist_addfld1d (fname='DEADCROOTC_STORAGE_TO_XFER', units='gC/m^2/s', &
             avgflag='A', long_name='dead coarse root C shift storage to transfer', &
             ptr_patch=this%deadcrootc_storage_to_xfer_patch, default='inactive')

        this%gresp_storage_to_xfer_patch(begp:endp) = spval
        call hist_addfld1d (fname='GRESP_STORAGE_TO_XFER', units='gC/m^2/s', &
             avgflag='A', long_name='growth respiration shift storage to transfer', &
             ptr_patch=this%gresp_storage_to_xfer_patch, default='inactive')

        this%livestemc_to_deadstemc_patch(begp:endp) = spval
        call hist_addfld1d (fname='LIVESTEMC_TO_DEADSTEMC', units='gC/m^2/s', &
             avgflag='A', long_name='live stem C turnover', &
             ptr_patch=this%livestemc_to_deadstemc_patch, default='inactive')

        this%livecrootc_to_deadcrootc_patch(begp:endp) = spval
        call hist_addfld1d (fname='LIVECROOTC_TO_DEADCROOTC', units='gC/m^2/s', &
             avgflag='A', long_name='live coarse root C turnover', &
             ptr_patch=this%livecrootc_to_deadcrootc_patch, default='inactive')

        this%gpp_patch(begp:endp) = spval
        call hist_addfld1d (fname='GPP', units='gC/m^2/s', &
             avgflag='A', long_name='gross primary production', &
             ptr_patch=this%gpp_patch)

        this%gpp_before_downreg_patch(begp:endp) = spval
        call hist_addfld1d (fname='INIT_GPP', units='gC/m^2/s', &
             avgflag='A', long_name='GPP flux before downregulation', &
             ptr_patch=this%gpp_before_downreg_patch, default='inactive')

        this%mr_patch(begp:endp) = spval
        call hist_addfld1d (fname='MR', units='gC/m^2/s', &
             avgflag='A', long_name='maintenance respiration', &
             ptr_patch=this%mr_patch)

        this%current_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='CURRENT_GR', units='gC/m^2/s', &
             avgflag='A', long_name='growth resp for new growth displayed in this timestep', &
             ptr_patch=this%current_gr_patch, default='inactive')

        this%transfer_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='TRANSFER_GR', units='gC/m^2/s', &
             avgflag='A', long_name='growth resp for transfer growth displayed in this timestep', &
             ptr_patch=this%transfer_gr_patch, default='inactive')

        this%storage_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='STORAGE_GR', units='gC/m^2/s', &
             avgflag='A', long_name='growth resp for growth sent to storage for later display', &
             ptr_patch=this%storage_gr_patch, default='inactive')

        this%gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='GR', units='gC/m^2/s', &
             avgflag='A', long_name='total growth respiration', &
             ptr_patch=this%gr_patch)

        this%ar_patch(begp:endp) = spval
        call hist_addfld1d (fname='AR', units='gC/m^2/s', &
             avgflag='A', long_name='autotrophic respiration (MR + GR)', &
             ptr_patch=this%ar_patch)

        this%rr_patch(begp:endp) = spval
        call hist_addfld1d (fname='RR', units='gC/m^2/s', &
             avgflag='A', long_name='root respiration (fine root MR + total root GR)', &
             ptr_patch=this%rr_patch)

        this%npp_patch(begp:endp) = spval
        call hist_addfld1d (fname='NPP', units='gC/m^2/s', &
             avgflag='A', long_name='net primary production', &
             ptr_patch=this%npp_patch)

        this%agnpp_patch(begp:endp) = spval
        call hist_addfld1d (fname='AGNPP', units='gC/m^2/s', &
             avgflag='A', long_name='aboveground NPP', &
             ptr_patch=this%agnpp_patch)

        this%bgnpp_patch(begp:endp) = spval
        call hist_addfld1d (fname='BGNPP', units='gC/m^2/s', &
             avgflag='A', long_name='belowground NPP', &
             ptr_patch=this%bgnpp_patch)

        this%litfall_patch(begp:endp) = spval
        call hist_addfld1d (fname='LITFALL', units='gC/m^2/s', &
             avgflag='A', long_name='litterfall (leaves and fine roots)', &
             ptr_patch=this%litfall_patch)

        this%vegfire_patch(begp:endp) = spval
        call hist_addfld1d (fname='VEGFIRE', units='gC/m^2/s', &
             avgflag='A', long_name='patch-level fire loss', &
             ptr_patch=this%vegfire_patch, default='inactive')

        this%wood_harvestc_patch(begp:endp) = spval
        call hist_addfld1d (fname='WOOD_HARVESTC', units='gC/m^2/s', &
             avgflag='A', long_name='wood harvest carbon (to product pools)', &
             ptr_patch=this%wood_harvestc_patch)

        this%fire_closs_patch(begp:endp) = spval
        call hist_addfld1d (fname='PFT_FIRE_CLOSS', units='gC/m^2/s', &
             avgflag='A', long_name='total patch-level fire C loss for non-peat fires outside land-type converted region', &
             ptr_patch=this%fire_closs_patch)

        this%availc_patch(begp:endp) = spval
        call hist_addfld1d (fname='AVAILC', units='gC/m^2/s', &
             avgflag='A', long_name='C flux available for allocation', &
             ptr_patch=this%availc_patch, default='active')

        this%plant_calloc_patch(begp:endp) = spval
        call hist_addfld1d (fname='PLANT_CALLOC', units='gC/m^2/s', &
             avgflag='A', long_name='total allocated C flux', &
             ptr_patch=this%plant_calloc_patch, default='active')

        this%excess_cflux_patch(begp:endp) = spval
        call hist_addfld1d (fname='EXCESS_CFLUX', units='gC/m^2/s', &
             avgflag='A', long_name='C flux not allocated due to downregulation', &
             ptr_patch=this%excess_cflux_patch, default='inactive')

        this%prev_leafc_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='PREV_LEAFC_TO_LITTER', units='gC/m^2/s', &
             avgflag='A', long_name='previous timestep leaf C litterfall flux', &
             ptr_patch=this%prev_leafc_to_litter_patch, default='inactive')

        this%prev_frootc_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='PREV_FROOTC_TO_LITTER', units='gC/m^2/s', &
             avgflag='A', long_name='previous timestep froot C litterfall flux', &
             ptr_patch=this%prev_frootc_to_litter_patch, default='inactive')

        this%xsmrpool_recover_patch(begp:endp) = spval
        call hist_addfld1d (fname='XSMRPOOL_RECOVER', units='gC/m^2/s', &
             avgflag='A', long_name='C flux assigned to recovery of negative xsmrpool', &
             ptr_patch=this%xsmrpool_recover_patch, default='inactive')

        if (nu_com .ne. 'RD' ) then
            this%allocation_leaf(begp:endp) = spval
            call hist_addfld1d (fname='allocation_leaf', units='', &
               avgflag='A', long_name='fraction of availc allocated to leaf', &
               ptr_patch=this%allocation_leaf)
            this%allocation_stem(begp:endp) = spval
            call hist_addfld1d (fname='allocation_stem', units='', &
               avgflag='A', long_name='fraction of availc allocated to stem', &
               ptr_patch=this%allocation_stem)
            this%allocation_froot(begp:endp) = spval
            call hist_addfld1d (fname='allocation_froot', units='', &
               avgflag='A', long_name='fraction of availc allocated to fine root', &
               ptr_patch=this%allocation_froot)
        end if

     end if  ! end of if-c12

     !-------------------------------
     ! C13 flux variables - native to PFT
     !-------------------------------
     if ( carbon_type == 'c13') then

        this%m_leafc_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_LEAFC_TO_LITTER', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 leaf C mortality', &
             ptr_patch=this%m_leafc_to_litter_patch, default='inactive')

        this%m_frootc_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_FROOTC_TO_LITTER', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 fine root C mortality', &
             ptr_patch=this%m_frootc_to_litter_patch, default='inactive')

        this%m_leafc_storage_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_LEAFC_STORAGE_TO_LITTER', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 leaf C storage mortality', &
             ptr_patch=this%m_leafc_storage_to_litter_patch, default='inactive')

        this%m_frootc_storage_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_FROOTC_STORAGE_TO_LITTER', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 fine root C storage mortality', &
             ptr_patch=this%m_frootc_storage_to_litter_patch, default='inactive')

        this%m_livestemc_storage_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_LIVESTEMC_STORAGE_TO_LITTER', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 live stem C storage mortality', &
             ptr_patch=this%m_livestemc_storage_to_litter_patch, default='inactive')

        this%m_deadstemc_storage_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_DEADSTEMC_STORAGE_TO_LITTER', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 dead stem C storage mortality', &
             ptr_patch=this%m_deadstemc_storage_to_litter_patch, default='inactive')

        this%m_livecrootc_storage_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_LIVECROOTC_STORAGE_TO_LITTER', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 live coarse root C storage mortality', &
             ptr_patch=this%m_livecrootc_storage_to_litter_patch, default='inactive')

        this%m_deadcrootc_storage_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_DEADCROOTC_STORAGE_TO_LITTER', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 dead coarse root C storage mortality', &
             ptr_patch=this%m_deadcrootc_storage_to_litter_patch, default='inactive')

        this%m_leafc_xfer_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_LEAFC_XFER_TO_LITTER', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 leaf C transfer mortality', &
             ptr_patch=this%m_leafc_xfer_to_litter_patch, default='inactive')

        this%m_frootc_xfer_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_FROOTC_XFER_TO_LITTER', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 fine root C transfer mortality', &
             ptr_patch=this%m_frootc_xfer_to_litter_patch, default='inactive')

        this%m_livestemc_xfer_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_LIVESTEMC_XFER_TO_LITTER', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 live stem C transfer mortality', &
             ptr_patch=this%m_livestemc_xfer_to_litter_patch, default='inactive')

        this%m_deadstemc_xfer_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_DEADSTEMC_XFER_TO_LITTER', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 dead stem C transfer mortality', &
             ptr_patch=this%m_deadstemc_xfer_to_litter_patch, default='inactive')

        this%m_livecrootc_xfer_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_LIVECROOTC_XFER_TO_LITTER', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 live coarse root C transfer mortality', &
             ptr_patch=this%m_livecrootc_xfer_to_litter_patch, default='inactive')

        this%m_deadcrootc_xfer_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_DEADCROOTC_XFER_TO_LITTER', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 dead coarse root C transfer mortality', &
             ptr_patch=this%m_deadcrootc_xfer_to_litter_patch, default='inactive')

        this%m_livestemc_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_LIVESTEMC_TO_LITTER', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 live stem C mortality', &
             ptr_patch=this%m_livestemc_to_litter_patch, default='inactive')

        this%m_deadstemc_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_DEADSTEMC_TO_LITTER', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 dead stem C mortality', &
             ptr_patch=this%m_deadstemc_to_litter_patch, default='inactive')

        this%m_livecrootc_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_LIVECROOTC_TO_LITTER', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 live coarse root C mortality', &
             ptr_patch=this%m_livecrootc_to_litter_patch, default='inactive')

        this%m_deadcrootc_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_DEADCROOTC_TO_LITTER', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 dead coarse root C mortality', &
             ptr_patch=this%m_deadcrootc_to_litter_patch, default='inactive')

        this%m_gresp_storage_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_GRESP_STORAGE_TO_LITTER', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 growth respiration storage mortality', &
             ptr_patch=this%m_gresp_storage_to_litter_patch, default='inactive')

        this%m_gresp_xfer_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_GRESP_XFER_TO_LITTER', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 growth respiration transfer mortality', &
             ptr_patch=this%m_gresp_xfer_to_litter_patch, default='inactive')

        this%m_leafc_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_LEAFC_TO_FIRE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 leaf C fire loss', &
             ptr_patch=this%m_leafc_to_fire_patch, default='inactive')

        this%m_frootc_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_FROOTC_TO_FIRE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 fine root C fire loss', &
             ptr_patch=this%m_frootc_to_fire_patch, default='inactive')

        this%m_leafc_storage_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_LEAFC_STORAGE_TO_FIRE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 leaf C storage fire loss', &
             ptr_patch=this%m_leafc_storage_to_fire_patch, default='inactive')

        this%m_frootc_storage_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_FROOTC_STORAGE_TO_FIRE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 fine root C storage fire loss', &
             ptr_patch=this%m_frootc_storage_to_fire_patch, default='inactive')

        this%m_livestemc_storage_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_LIVESTEMC_STORAGE_TO_FIRE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 live stem C storage fire loss', &
             ptr_patch=this%m_livestemc_storage_to_fire_patch, default='inactive')

        this%m_deadstemc_storage_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_DEADSTEMC_STORAGE_TO_FIRE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 dead stem C storage fire loss', &
             ptr_patch=this%m_deadstemc_storage_to_fire_patch, default='inactive')

        this%m_livecrootc_storage_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_LIVECROOTC_STORAGE_TO_FIRE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 live coarse root C storage fire loss', &
             ptr_patch=this%m_livecrootc_storage_to_fire_patch, default='inactive')

        this%m_deadcrootc_storage_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_DEADCROOTC_STORAGE_TO_FIRE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 dead coarse root C storage fire loss', &
             ptr_patch=this%m_deadcrootc_storage_to_fire_patch,  default='inactive')

        this%m_leafc_xfer_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_LEAFC_XFER_TO_FIRE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 leaf C transfer fire loss', &
             ptr_patch=this%m_leafc_xfer_to_fire_patch, default='inactive')

        this%m_frootc_xfer_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_FROOTC_XFER_TO_FIRE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 fine root C transfer fire loss', &
             ptr_patch=this%m_frootc_xfer_to_fire_patch, default='inactive')

        this%m_livestemc_xfer_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_LIVESTEMC_XFER_TO_FIRE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 live stem C transfer fire loss', &
             ptr_patch=this%m_livestemc_xfer_to_fire_patch, default='inactive')

        this%m_deadstemc_xfer_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_DEADSTEMC_XFER_TO_FIRE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 dead stem C transfer fire loss', &
             ptr_patch=this%m_deadstemc_xfer_to_fire_patch, default='inactive')

        this%m_livecrootc_xfer_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_LIVECROOTC_XFER_TO_FIRE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 live coarse root C transfer fire loss', &
             ptr_patch=this%m_livecrootc_xfer_to_fire_patch, default='inactive')

        this%m_deadcrootc_xfer_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_DEADCROOTC_XFER_TO_FIRE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 dead coarse root C transfer fire loss', &
             ptr_patch=this%m_deadcrootc_xfer_to_fire_patch, default='inactive')

        this%m_livestemc_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_LIVESTEMC_TO_FIRE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 live stem C fire loss', &
             ptr_patch=this%m_livestemc_to_fire_patch, default='inactive')

        this%m_deadstemc_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_DEADSTEMC_TO_FIRE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 dead stem C fire loss', &
             ptr_patch=this%m_deadstemc_to_fire_patch, default='inactive')

        this%m_deadstemc_to_litter_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_DEADSTEMC_TO_LITTER_FIRE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 dead stem C fire mortality to litter', &
             ptr_patch=this%m_deadstemc_to_litter_fire_patch, default='inactive')

        this%m_livecrootc_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_LIVECROOTC_TO_FIRE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 live coarse root C fire loss', &
             ptr_patch=this%m_livecrootc_to_fire_patch, default='inactive')

        this%m_deadcrootc_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_DEADCROOTC_TO_FIRE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 dead coarse root C fire loss', &
             ptr_patch=this%m_deadcrootc_to_fire_patch, default='inactive')

        this%m_deadcrootc_to_litter_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_DEADCROOTC_TO_LITTER_FIRE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 dead coarse root C fire mortality to litter', &
             ptr_patch=this%m_deadcrootc_to_litter_fire_patch, default='inactive')

        this%m_gresp_storage_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_GRESP_STORAGE_TO_FIRE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 growth respiration storage fire loss', &
             ptr_patch=this%m_gresp_storage_to_fire_patch, default='inactive')

        this%m_gresp_xfer_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_M_GRESP_XFER_TO_FIRE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 growth respiration transfer fire loss', &
             ptr_patch=this%m_gresp_xfer_to_fire_patch, default='inactive')

        this%leafc_xfer_to_leafc_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_LEAFC_XFER_TO_LEAFC', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 leaf C growth from storage', &
             ptr_patch=this%leafc_xfer_to_leafc_patch, default='inactive')

        this%frootc_xfer_to_frootc_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_FROOTC_XFER_TO_FROOTC', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 fine root C growth from storage', &
             ptr_patch=this%frootc_xfer_to_frootc_patch, default='inactive')

        this%livestemc_xfer_to_livestemc_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_LIVESTEMC_XFER_TO_LIVESTEMC', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 live stem C growth from storage', &
             ptr_patch=this%livestemc_xfer_to_livestemc_patch, default='inactive')

        this%deadstemc_xfer_to_deadstemc_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_DEADSTEMC_XFER_TO_DEADSTEMC', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 dead stem C growth from storage', &
             ptr_patch=this%deadstemc_xfer_to_deadstemc_patch, default='inactive')

        this%livecrootc_xfer_to_livecrootc_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_LIVECROOTC_XFER_TO_LIVECROOTC', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 live coarse root C growth from storage', &
             ptr_patch=this%livecrootc_xfer_to_livecrootc_patch, default='inactive')

        this%deadcrootc_xfer_to_deadcrootc_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_DEADCROOTC_XFER_TO_DEADCROOTC', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 dead coarse root C growth from storage', &
             ptr_patch=this%deadcrootc_xfer_to_deadcrootc_patch, default='inactive')

        this%leafc_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_LEAFC_TO_LITTER', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 leaf C litterfall', &
             ptr_patch=this%leafc_to_litter_patch, default='inactive')

        this%frootc_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_FROOTC_TO_LITTER', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 fine root C litterfall', &
             ptr_patch=this%frootc_to_litter_patch, default='inactive')

        this%leaf_mr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_LEAF_MR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 leaf maintenance respiration', &
             ptr_patch=this%leaf_mr_patch, default='inactive')

        this%froot_mr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_FROOT_MR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 fine root maintenance respiration', &
             ptr_patch=this%froot_mr_patch, default='inactive')

        this%livestem_mr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_LIVESTEM_MR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 live stem maintenance respiration', &
             ptr_patch=this%livestem_mr_patch, default='inactive')

        this%livecroot_mr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_LIVECROOT_MR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 live coarse root maintenance respiration', &
             ptr_patch=this%livecroot_mr_patch, default='inactive')

        this%psnsun_to_cpool_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_PSNSUN_TO_CPOOL', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 C fixation from sunlit canopy', &
             ptr_patch=this%psnsun_to_cpool_patch)

        this%psnshade_to_cpool_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_PSNSHADE_TO_CPOOL', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 C fixation from shaded canopy', &
             ptr_patch=this%psnshade_to_cpool_patch)

        this%cpool_to_leafc_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_CPOOL_TO_LEAFC', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 allocation to leaf C', &
             ptr_patch=this%cpool_to_leafc_patch, default='inactive')

        this%cpool_to_leafc_storage_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_CPOOL_TO_LEAFC_STORAGE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 allocation to leaf C storage', &
             ptr_patch=this%cpool_to_leafc_storage_patch, default='inactive')

        this%cpool_to_frootc_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_CPOOL_TO_FROOTC', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 allocation to fine root C', &
             ptr_patch=this%cpool_to_frootc_patch, default='inactive')

        this%cpool_to_frootc_storage_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_CPOOL_TO_FROOTC_STORAGE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 allocation to fine root C storage', &
             ptr_patch=this%cpool_to_frootc_storage_patch, default='inactive')

        this%cpool_to_livestemc_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_CPOOL_TO_LIVESTEMC', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 allocation to live stem C', &
             ptr_patch=this%cpool_to_livestemc_patch, default='inactive')

        this%cpool_to_livestemc_storage_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_CPOOL_TO_LIVESTEMC_STORAGE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 allocation to live stem C storage', &
             ptr_patch=this%cpool_to_livestemc_storage_patch, default='inactive')

        this%cpool_to_deadstemc_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_CPOOL_TO_DEADSTEMC', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 allocation to dead stem C', &
             ptr_patch=this%cpool_to_deadstemc_patch, default='inactive')

        this%cpool_to_deadstemc_storage_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_CPOOL_TO_DEADSTEMC_STORAGE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 allocation to dead stem C storage', &
             ptr_patch=this%cpool_to_deadstemc_storage_patch, default='inactive')

        this%cpool_to_livecrootc_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_CPOOL_TO_LIVECROOTC', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 allocation to live coarse root C', &
             ptr_patch=this%cpool_to_livecrootc_patch, default='inactive')

        this%cpool_to_livecrootc_storage_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_CPOOL_TO_LIVECROOTC_STORAGE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 allocation to live coarse root C storage', &
             ptr_patch=this%cpool_to_livecrootc_storage_patch, default='inactive')

        this%cpool_to_deadcrootc_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_CPOOL_TO_DEADCROOTC', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 allocation to dead coarse root C', &
             ptr_patch=this%cpool_to_deadcrootc_patch, default='inactive')

        this%cpool_to_deadcrootc_storage_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_CPOOL_TO_DEADCROOTC_STORAGE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 allocation to dead coarse root C storage', &
             ptr_patch=this%cpool_to_deadcrootc_storage_patch, default='inactive')

        this%cpool_to_gresp_storage_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_CPOOL_TO_GRESP_STORAGE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 allocation to growth respiration storage', &
             ptr_patch=this%cpool_to_gresp_storage_patch, default='inactive')

        this%cpool_leaf_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_CPOOL_LEAF_GR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 leaf growth respiration', &
             ptr_patch=this%cpool_leaf_gr_patch, default='inactive')

        this%cpool_leaf_storage_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_CPOOL_LEAF_STORAGE_GR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 leaf growth respiration to storage', &
             ptr_patch=this%cpool_leaf_storage_gr_patch, default='inactive')

        this%transfer_leaf_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_TRANSFER_LEAF_GR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 leaf growth respiration from storage', &
             ptr_patch=this%transfer_leaf_gr_patch, default='inactive')

        this%cpool_froot_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_CPOOL_FROOT_GR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 fine root growth respiration', &
             ptr_patch=this%cpool_froot_gr_patch, default='inactive')

        this%cpool_froot_storage_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_CPOOL_FROOT_STORAGE_GR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 fine root  growth respiration to storage', &
             ptr_patch=this%cpool_froot_storage_gr_patch, default='inactive')

        this%transfer_froot_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_TRANSFER_FROOT_GR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 fine root  growth respiration from storage', &
             ptr_patch=this%transfer_froot_gr_patch, default='inactive')

        this%cpool_livestem_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_CPOOL_LIVESTEM_GR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 live stem growth respiration', &
             ptr_patch=this%cpool_livestem_gr_patch, default='inactive')

        this%cpool_livestem_storage_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_CPOOL_LIVESTEM_STORAGE_GR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 live stem growth respiration to storage', &
             ptr_patch=this%cpool_livestem_storage_gr_patch, default='inactive')

        this%transfer_livestem_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_TRANSFER_LIVESTEM_GR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 live stem growth respiration from storage', &
             ptr_patch=this%transfer_livestem_gr_patch, default='inactive')

        this%cpool_deadstem_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_CPOOL_DEADSTEM_GR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 dead stem growth respiration', &
             ptr_patch=this%cpool_deadstem_gr_patch, default='inactive')

        this%cpool_deadstem_storage_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_CPOOL_DEADSTEM_STORAGE_GR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 dead stem growth respiration to storage', &
             ptr_patch=this%cpool_deadstem_storage_gr_patch, default='inactive')

        this%transfer_deadstem_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_TRANSFER_DEADSTEM_GR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 dead stem growth respiration from storage', &
             ptr_patch=this%transfer_deadstem_gr_patch, default='inactive')

        this%cpool_livecroot_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_CPOOL_LIVECROOT_GR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 live coarse root growth respiration', &
             ptr_patch=this%cpool_livecroot_gr_patch, default='inactive')

        this%cpool_livecroot_storage_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_CPOOL_LIVECROOT_STORAGE_GR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 live coarse root growth respiration to storage', &
             ptr_patch=this%cpool_livecroot_storage_gr_patch, default='inactive')

        this%transfer_livecroot_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_TRANSFER_LIVECROOT_GR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 live coarse root growth respiration from storage', &
             ptr_patch=this%transfer_livecroot_gr_patch, default='inactive')

        this%cpool_deadcroot_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_CPOOL_DEADCROOT_GR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 dead coarse root growth respiration', &
             ptr_patch=this%cpool_deadcroot_gr_patch, default='inactive')

        this%cpool_deadcroot_storage_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_CPOOL_DEADCROOT_STORAGE_GR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 dead coarse root growth respiration to storage', &
             ptr_patch=this%cpool_deadcroot_storage_gr_patch, default='inactive')

        this%transfer_deadcroot_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_TRANSFER_DEADCROOT_GR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 dead coarse root growth respiration from storage', &
             ptr_patch=this%transfer_deadcroot_gr_patch, default='inactive')

        this%leafc_storage_to_xfer_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_LEAFC_STORAGE_TO_XFER', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 leaf C shift storage to transfer', &
             ptr_patch=this%leafc_storage_to_xfer_patch, default='inactive')

        this%frootc_storage_to_xfer_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_FROOTC_STORAGE_TO_XFER', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 fine root C shift storage to transfer', &
             ptr_patch=this%frootc_storage_to_xfer_patch, default='inactive')

        this%livestemc_storage_to_xfer_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_LIVESTEMC_STORAGE_TO_XFER', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 live stem C shift storage to transfer', &
             ptr_patch=this%livestemc_storage_to_xfer_patch, default='inactive')

        this%deadstemc_storage_to_xfer_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_DEADSTEMC_STORAGE_TO_XFER', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 dead stem C shift storage to transfer', &
             ptr_patch=this%deadstemc_storage_to_xfer_patch, default='inactive')

        this%livecrootc_storage_to_xfer_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_LIVECROOTC_STORAGE_TO_XFER', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 live coarse root C shift storage to transfer', &
             ptr_patch=this%livecrootc_storage_to_xfer_patch, default='inactive')

        this%deadcrootc_storage_to_xfer_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_DEADCROOTC_STORAGE_TO_XFER', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 dead coarse root C shift storage to transfer', &
             ptr_patch=this%deadcrootc_storage_to_xfer_patch, default='inactive')

        this%gresp_storage_to_xfer_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_GRESP_STORAGE_TO_XFER', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 growth respiration shift storage to transfer', &
             ptr_patch=this%gresp_storage_to_xfer_patch, default='inactive')

        this%livestemc_to_deadstemc_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_LIVESTEMC_TO_DEADSTEMC', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 live stem C turnover', &
             ptr_patch=this%livestemc_to_deadstemc_patch, default='inactive')

        this%livecrootc_to_deadcrootc_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_LIVECROOTC_TO_DEADCROOTC', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 live coarse root C turnover', &
             ptr_patch=this%livecrootc_to_deadcrootc_patch, default='inactive')

        this%gpp_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_GPP', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 gross primary production', &
             ptr_patch=this%gpp_patch)

        this%mr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_MR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 maintenance respiration', &
             ptr_patch=this%mr_patch)

        this%current_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_CURRENT_GR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 growth resp for new growth displayed in this timestep', &
             ptr_patch=this%current_gr_patch, default='inactive')

        this%transfer_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_TRANSFER_GR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 growth resp for transfer growth displayed in this timestep', &
             ptr_patch=this%transfer_gr_patch, default='inactive')

        this%storage_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_STORAGE_GR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 growth resp for growth sent to storage for later display', &
             ptr_patch=this%storage_gr_patch, default='inactive')

        this%gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_GR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 total growth respiration', &
             ptr_patch=this%gr_patch)

        this%ar_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_AR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 autotrophic respiration (MR + GR)', &
             ptr_patch=this%ar_patch)

        this%rr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_RR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 root respiration (fine root MR + total root GR)', &
             ptr_patch=this%rr_patch)

        this%npp_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_NPP', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 net primary production', &
             ptr_patch=this%npp_patch)

        this%agnpp_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_AGNPP', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 aboveground NPP', &
             ptr_patch=this%agnpp_patch)

        this%bgnpp_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_BGNPP', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 belowground NPP', &
             ptr_patch=this%bgnpp_patch)

        this%litfall_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_LITFALL', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 litterfall (leaves and fine roots)', &
             ptr_patch=this%litfall_patch, default='inactive')

        this%vegfire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_VEGFIRE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 patch-level fire loss', &
             ptr_patch=this%vegfire_patch, default='inactive')

        this%fire_closs_patch(begp:endp) = spval
        call hist_addfld1d (fname='C13_PFT_FIRE_CLOSS', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 total patch-level fire C loss', &
             ptr_patch=this%fire_closs_patch)
     endif

     !-------------------------------
     ! C14 flux variables - native to PFT
     !-------------------------------
     if ( carbon_type == 'c14' ) then

        this%m_leafc_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_LEAFC_TO_LITTER', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 leaf C mortality', &
             ptr_patch=this%m_leafc_to_litter_patch, default='inactive')

        this%m_frootc_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_FROOTC_TO_LITTER', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 fine root C mortality', &
             ptr_patch=this%m_frootc_to_litter_patch, default='inactive')

        this%m_leafc_storage_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_LEAFC_STORAGE_TO_LITTER', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 leaf C storage mortality', &
             ptr_patch=this%m_leafc_storage_to_litter_patch, default='inactive')

        this%m_frootc_storage_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_FROOTC_STORAGE_TO_LITTER', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 fine root C storage mortality', &
             ptr_patch=this%m_frootc_storage_to_litter_patch, default='inactive')

        this%m_livestemc_storage_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_LIVESTEMC_STORAGE_TO_LITTER', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 live stem C storage mortality', &
             ptr_patch=this%m_livestemc_storage_to_litter_patch, default='inactive')

        this%m_deadstemc_storage_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_DEADSTEMC_STORAGE_TO_LITTER', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 dead stem C storage mortality', &
             ptr_patch=this%m_deadstemc_storage_to_litter_patch, default='inactive')

        this%m_livecrootc_storage_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_LIVECROOTC_STORAGE_TO_LITTER', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 live coarse root C storage mortality', &
             ptr_patch=this%m_livecrootc_storage_to_litter_patch, default='inactive')

        this%m_deadcrootc_storage_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_DEADCROOTC_STORAGE_TO_LITTER', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 dead coarse root C storage mortality', &
             ptr_patch=this%m_deadcrootc_storage_to_litter_patch, default='inactive')

        this%m_leafc_xfer_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_LEAFC_XFER_TO_LITTER', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 leaf C transfer mortality', &
             ptr_patch=this%m_leafc_xfer_to_litter_patch, default='inactive')

        this%m_frootc_xfer_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_FROOTC_XFER_TO_LITTER', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 fine root C transfer mortality', &
             ptr_patch=this%m_frootc_xfer_to_litter_patch, default='inactive')

        this%m_livestemc_xfer_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_LIVESTEMC_XFER_TO_LITTER', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 live stem C transfer mortality', &
             ptr_patch=this%m_livestemc_xfer_to_litter_patch, default='inactive')

        this%m_deadstemc_xfer_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_DEADSTEMC_XFER_TO_LITTER', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 dead stem C transfer mortality', &
             ptr_patch=this%m_deadstemc_xfer_to_litter_patch, default='inactive')

        this%m_livecrootc_xfer_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_LIVECROOTC_XFER_TO_LITTER', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 live coarse root C transfer mortality', &
             ptr_patch=this%m_livecrootc_xfer_to_litter_patch, default='inactive')

        this%m_deadcrootc_xfer_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_DEADCROOTC_XFER_TO_LITTER', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 dead coarse root C transfer mortality', &
             ptr_patch=this%m_deadcrootc_xfer_to_litter_patch, default='inactive')

        this%m_livestemc_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_LIVESTEMC_TO_LITTER', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 live stem C mortality', &
             ptr_patch=this%m_livestemc_to_litter_patch, default='inactive')

        this%m_deadstemc_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_DEADSTEMC_TO_LITTER', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 dead stem C mortality', &
             ptr_patch=this%m_deadstemc_to_litter_patch, default='inactive')

        this%m_livecrootc_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_LIVECROOTC_TO_LITTER', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 live coarse root C mortality', &
             ptr_patch=this%m_livecrootc_to_litter_patch, default='inactive')

        this%m_deadcrootc_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_DEADCROOTC_TO_LITTER', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 dead coarse root C mortality', &
             ptr_patch=this%m_deadcrootc_to_litter_patch, default='inactive')

        this%m_gresp_storage_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_GRESP_STORAGE_TO_LITTER', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 growth respiration storage mortality', &
             ptr_patch=this%m_gresp_storage_to_litter_patch, default='inactive')

        this%m_gresp_xfer_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_GRESP_XFER_TO_LITTER', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 growth respiration transfer mortality', &
             ptr_patch=this%m_gresp_xfer_to_litter_patch, default='inactive')

        this%m_leafc_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_LEAFC_TO_FIRE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 leaf C fire loss', &
             ptr_patch=this%m_leafc_to_fire_patch, default='inactive')

        this%m_frootc_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_FROOTC_TO_FIRE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 fine root C fire loss', &
             ptr_patch=this%m_frootc_to_fire_patch, default='inactive')

        this%m_leafc_storage_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_LEAFC_STORAGE_TO_FIRE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 leaf C storage fire loss', &
             ptr_patch=this%m_leafc_storage_to_fire_patch, default='inactive')

        this%m_frootc_storage_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_FROOTC_STORAGE_TO_FIRE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 fine root C storage fire loss', &
             ptr_patch=this%m_frootc_storage_to_fire_patch, default='inactive')

        this%m_livestemc_storage_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_LIVESTEMC_STORAGE_TO_FIRE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 live stem C storage fire loss', &
             ptr_patch=this%m_livestemc_storage_to_fire_patch, default='inactive')

        this%m_deadstemc_storage_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_DEADSTEMC_STORAGE_TO_FIRE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 dead stem C storage fire loss', &
             ptr_patch=this%m_deadstemc_storage_to_fire_patch, default='inactive')

        this%m_livecrootc_storage_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_LIVECROOTC_STORAGE_TO_FIRE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 live coarse root C storage fire loss', &
             ptr_patch=this%m_livecrootc_storage_to_fire_patch, default='inactive')

        this%m_deadcrootc_storage_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_DEADCROOTC_STORAGE_TO_FIRE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 dead coarse root C storage fire loss', &
             ptr_patch=this%m_deadcrootc_storage_to_fire_patch,  default='inactive')

        this%m_leafc_xfer_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_LEAFC_XFER_TO_FIRE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 leaf C transfer fire loss', &
             ptr_patch=this%m_leafc_xfer_to_fire_patch, default='inactive')

        this%m_frootc_xfer_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_FROOTC_XFER_TO_FIRE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 fine root C transfer fire loss', &
             ptr_patch=this%m_frootc_xfer_to_fire_patch, default='inactive')

        this%m_livestemc_xfer_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_LIVESTEMC_XFER_TO_FIRE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 live stem C transfer fire loss', &
             ptr_patch=this%m_livestemc_xfer_to_fire_patch, default='inactive')

        this%m_deadstemc_xfer_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_DEADSTEMC_XFER_TO_FIRE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 dead stem C transfer fire loss', &
             ptr_patch=this%m_deadstemc_xfer_to_fire_patch, default='inactive')

        this%m_livecrootc_xfer_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_LIVECROOTC_XFER_TO_FIRE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 live coarse root C transfer fire loss', &
             ptr_patch=this%m_livecrootc_xfer_to_fire_patch, default='inactive')

        this%m_deadcrootc_xfer_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_DEADCROOTC_XFER_TO_FIRE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 dead coarse root C transfer fire loss', &
             ptr_patch=this%m_deadcrootc_xfer_to_fire_patch, default='inactive')

        this%m_livestemc_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_LIVESTEMC_TO_FIRE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 live stem C fire loss', &
             ptr_patch=this%m_livestemc_to_fire_patch, default='inactive')

        this%m_deadstemc_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_DEADSTEMC_TO_FIRE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 dead stem C fire loss', &
             ptr_patch=this%m_deadstemc_to_fire_patch, default='inactive')

        this%m_deadstemc_to_litter_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_DEADSTEMC_TO_LITTER_FIRE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 dead stem C fire mortality to litter', &
             ptr_patch=this%m_deadstemc_to_litter_fire_patch, default='inactive')

        this%m_livecrootc_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_LIVECROOTC_TO_FIRE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 live coarse root C fire loss', &
             ptr_patch=this%m_livecrootc_to_fire_patch, default='inactive')

        this%m_deadcrootc_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_DEADCROOTC_TO_FIRE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 dead coarse root C fire loss', &
             ptr_patch=this%m_deadcrootc_to_fire_patch, default='inactive')

        this%m_deadcrootc_to_litter_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_DEADCROOTC_TO_LITTER_FIRE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 dead coarse root C fire mortality to litter', &
             ptr_patch=this%m_deadcrootc_to_litter_fire_patch, default='inactive')

        this%m_gresp_storage_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_GRESP_STORAGE_TO_FIRE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 growth respiration storage fire loss', &
             ptr_patch=this%m_gresp_storage_to_fire_patch, default='inactive')

        this%m_gresp_xfer_to_fire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_M_GRESP_XFER_TO_FIRE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 growth respiration transfer fire loss', &
             ptr_patch=this%m_gresp_xfer_to_fire_patch, default='inactive')

        this%leafc_xfer_to_leafc_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_LEAFC_XFER_TO_LEAFC', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 leaf C growth from storage', &
             ptr_patch=this%leafc_xfer_to_leafc_patch, default='inactive')

        this%frootc_xfer_to_frootc_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_FROOTC_XFER_TO_FROOTC', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 fine root C growth from storage', &
             ptr_patch=this%frootc_xfer_to_frootc_patch, default='inactive')

        this%livestemc_xfer_to_livestemc_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_LIVESTEMC_XFER_TO_LIVESTEMC', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 live stem C growth from storage', &
             ptr_patch=this%livestemc_xfer_to_livestemc_patch, default='inactive')

        this%deadstemc_xfer_to_deadstemc_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_DEADSTEMC_XFER_TO_DEADSTEMC', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 dead stem C growth from storage', &
             ptr_patch=this%deadstemc_xfer_to_deadstemc_patch, default='inactive')

        this%livecrootc_xfer_to_livecrootc_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_LIVECROOTC_XFER_TO_LIVECROOTC', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 live coarse root C growth from storage', &
             ptr_patch=this%livecrootc_xfer_to_livecrootc_patch, default='inactive')

        this%deadcrootc_xfer_to_deadcrootc_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_DEADCROOTC_XFER_TO_DEADCROOTC', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 dead coarse root C growth from storage', &
             ptr_patch=this%deadcrootc_xfer_to_deadcrootc_patch, default='inactive')

        this%leafc_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_LEAFC_TO_LITTER', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 leaf C litterfall', &
             ptr_patch=this%leafc_to_litter_patch, default='inactive')

        this%frootc_to_litter_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_FROOTC_TO_LITTER', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 fine root C litterfall', &
             ptr_patch=this%frootc_to_litter_patch, default='inactive')

        this%leaf_mr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_LEAF_MR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 leaf maintenance respiration', &
             ptr_patch=this%leaf_mr_patch, default='inactive')

        this%froot_mr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_FROOT_MR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 fine root maintenance respiration', &
             ptr_patch=this%froot_mr_patch, default='inactive')

        this%livestem_mr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_LIVESTEM_MR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 live stem maintenance respiration', &
             ptr_patch=this%livestem_mr_patch, default='inactive')

        this%livecroot_mr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_LIVECROOT_MR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 live coarse root maintenance respiration', &
             ptr_patch=this%livecroot_mr_patch, default='inactive')

        this%psnsun_to_cpool_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_PSNSUN_TO_CPOOL', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 C fixation from sunlit canopy', &
             ptr_patch=this%psnsun_to_cpool_patch)

        this%psnshade_to_cpool_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_PSNSHADE_TO_CPOOL', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 C fixation from shaded canopy', &
             ptr_patch=this%psnshade_to_cpool_patch)

        this%cpool_to_leafc_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_CPOOL_TO_LEAFC', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 allocation to leaf C', &
             ptr_patch=this%cpool_to_leafc_patch, default='inactive')

        this%cpool_to_leafc_storage_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_CPOOL_TO_LEAFC_STORAGE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 allocation to leaf C storage', &
             ptr_patch=this%cpool_to_leafc_storage_patch, default='inactive')

        this%cpool_to_frootc_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_CPOOL_TO_FROOTC', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 allocation to fine root C', &
             ptr_patch=this%cpool_to_frootc_patch, default='inactive')

        this%cpool_to_frootc_storage_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_CPOOL_TO_FROOTC_STORAGE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 allocation to fine root C storage', &
             ptr_patch=this%cpool_to_frootc_storage_patch, default='inactive')

        this%cpool_to_livestemc_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_CPOOL_TO_LIVESTEMC', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 allocation to live stem C', &
             ptr_patch=this%cpool_to_livestemc_patch, default='inactive')

        this%cpool_to_livestemc_storage_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_CPOOL_TO_LIVESTEMC_STORAGE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 allocation to live stem C storage', &
             ptr_patch=this%cpool_to_livestemc_storage_patch, default='inactive')

        this%cpool_to_deadstemc_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_CPOOL_TO_DEADSTEMC', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 allocation to dead stem C', &
             ptr_patch=this%cpool_to_deadstemc_patch, default='inactive')

        this%cpool_to_deadstemc_storage_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_CPOOL_TO_DEADSTEMC_STORAGE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 allocation to dead stem C storage', &
             ptr_patch=this%cpool_to_deadstemc_storage_patch, default='inactive')

        this%cpool_to_livecrootc_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_CPOOL_TO_LIVECROOTC', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 allocation to live coarse root C', &
             ptr_patch=this%cpool_to_livecrootc_patch, default='inactive')

        this%cpool_to_livecrootc_storage_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_CPOOL_TO_LIVECROOTC_STORAGE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 allocation to live coarse root C storage', &
             ptr_patch=this%cpool_to_livecrootc_storage_patch, default='inactive')

        this%cpool_to_deadcrootc_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_CPOOL_TO_DEADCROOTC', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 allocation to dead coarse root C', &
             ptr_patch=this%cpool_to_deadcrootc_patch, default='inactive')

        this%cpool_to_deadcrootc_storage_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_CPOOL_TO_DEADCROOTC_STORAGE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 allocation to dead coarse root C storage', &
             ptr_patch=this%cpool_to_deadcrootc_storage_patch, default='inactive')

        this%cpool_to_gresp_storage_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_CPOOL_TO_GRESP_STORAGE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 allocation to growth respiration storage', &
             ptr_patch=this%cpool_to_gresp_storage_patch, default='inactive')

        this%cpool_leaf_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_CPOOL_LEAF_GR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 leaf growth respiration', &
             ptr_patch=this%cpool_leaf_gr_patch, default='inactive')

        this%cpool_leaf_storage_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_CPOOL_LEAF_STORAGE_GR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 leaf growth respiration to storage', &
             ptr_patch=this%cpool_leaf_storage_gr_patch, default='inactive')

        this%transfer_leaf_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_TRANSFER_LEAF_GR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 leaf growth respiration from storage', &
             ptr_patch=this%transfer_leaf_gr_patch, default='inactive')

        this%cpool_froot_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_CPOOL_FROOT_GR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 fine root growth respiration', &
             ptr_patch=this%cpool_froot_gr_patch, default='inactive')

        this%cpool_froot_storage_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_CPOOL_FROOT_STORAGE_GR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 fine root  growth respiration to storage', &
             ptr_patch=this%cpool_froot_storage_gr_patch, default='inactive')

        this%transfer_froot_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_TRANSFER_FROOT_GR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 fine root  growth respiration from storage', &
             ptr_patch=this%transfer_froot_gr_patch, default='inactive')

        this%cpool_livestem_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_CPOOL_LIVESTEM_GR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 live stem growth respiration', &
             ptr_patch=this%cpool_livestem_gr_patch, default='inactive')

        this%cpool_livestem_storage_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_CPOOL_LIVESTEM_STORAGE_GR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 live stem growth respiration to storage', &
             ptr_patch=this%cpool_livestem_storage_gr_patch, default='inactive')

        this%transfer_livestem_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_TRANSFER_LIVESTEM_GR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 live stem growth respiration from storage', &
             ptr_patch=this%transfer_livestem_gr_patch, default='inactive')

        this%cpool_deadstem_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_CPOOL_DEADSTEM_GR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 dead stem growth respiration', &
             ptr_patch=this%cpool_deadstem_gr_patch, default='inactive')

        this%cpool_deadstem_storage_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_CPOOL_DEADSTEM_STORAGE_GR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 dead stem growth respiration to storage', &
             ptr_patch=this%cpool_deadstem_storage_gr_patch, default='inactive')

        this%transfer_deadstem_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_TRANSFER_DEADSTEM_GR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 dead stem growth respiration from storage', &
             ptr_patch=this%transfer_deadstem_gr_patch, default='inactive')

        this%cpool_livecroot_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_CPOOL_LIVECROOT_GR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 live coarse root growth respiration', &
             ptr_patch=this%cpool_livecroot_gr_patch, default='inactive')

        this%cpool_livecroot_storage_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_CPOOL_LIVECROOT_STORAGE_GR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 live coarse root growth respiration to storage', &
             ptr_patch=this%cpool_livecroot_storage_gr_patch, default='inactive')

        this%transfer_livecroot_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_TRANSFER_LIVECROOT_GR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 live coarse root growth respiration from storage', &
             ptr_patch=this%transfer_livecroot_gr_patch, default='inactive')

        this%cpool_deadcroot_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_CPOOL_DEADCROOT_GR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 dead coarse root growth respiration', &
             ptr_patch=this%cpool_deadcroot_gr_patch, default='inactive')

        this%cpool_deadcroot_storage_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_CPOOL_DEADCROOT_STORAGE_GR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 dead coarse root growth respiration to storage', &
             ptr_patch=this%cpool_deadcroot_storage_gr_patch, default='inactive')

        this%transfer_deadcroot_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_TRANSFER_DEADCROOT_GR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 dead coarse root growth respiration from storage', &
             ptr_patch=this%transfer_deadcroot_gr_patch, default='inactive')

        this%leafc_storage_to_xfer_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_LEAFC_STORAGE_TO_XFER', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 leaf C shift storage to transfer', &
             ptr_patch=this%leafc_storage_to_xfer_patch, default='inactive')

        this%frootc_storage_to_xfer_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_FROOTC_STORAGE_TO_XFER', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 fine root C shift storage to transfer', &
             ptr_patch=this%frootc_storage_to_xfer_patch, default='inactive')

        this%livestemc_storage_to_xfer_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_LIVESTEMC_STORAGE_TO_XFER', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 live stem C shift storage to transfer', &
             ptr_patch=this%livestemc_storage_to_xfer_patch, default='inactive')

        this%deadstemc_storage_to_xfer_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_DEADSTEMC_STORAGE_TO_XFER', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 dead stem C shift storage to transfer', &
             ptr_patch=this%deadstemc_storage_to_xfer_patch, default='inactive')

        this%livecrootc_storage_to_xfer_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_LIVECROOTC_STORAGE_TO_XFER', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 live coarse root C shift storage to transfer', &
             ptr_patch=this%livecrootc_storage_to_xfer_patch, default='inactive')

        this%deadcrootc_storage_to_xfer_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_DEADCROOTC_STORAGE_TO_XFER', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 dead coarse root C shift storage to transfer', &
             ptr_patch=this%deadcrootc_storage_to_xfer_patch, default='inactive')

        this%gresp_storage_to_xfer_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_GRESP_STORAGE_TO_XFER', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 growth respiration shift storage to transfer', &
             ptr_patch=this%gresp_storage_to_xfer_patch, default='inactive')

        this%livestemc_to_deadstemc_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_LIVESTEMC_TO_DEADSTEMC', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 live stem C turnover', &
             ptr_patch=this%livestemc_to_deadstemc_patch, default='inactive')

        this%livecrootc_to_deadcrootc_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_LIVECROOTC_TO_DEADCROOTC', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 live coarse root C turnover', &
             ptr_patch=this%livecrootc_to_deadcrootc_patch, default='inactive')

        this%gpp_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_GPP', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 gross primary production', &
             ptr_patch=this%gpp_patch)

        this%mr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_MR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 maintenance respiration', &
             ptr_patch=this%mr_patch)

        this%current_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_CURRENT_GR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 growth resp for new growth displayed in this timestep', &
             ptr_patch=this%current_gr_patch, default='inactive')

        this%transfer_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_TRANSFER_GR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 growth resp for transfer growth displayed in this timestep', &
             ptr_patch=this%transfer_gr_patch, default='inactive')

        this%storage_gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_STORAGE_GR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 growth resp for growth sent to storage for later display', &
             ptr_patch=this%storage_gr_patch, default='inactive')

        this%gr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_GR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 total growth respiration', &
             ptr_patch=this%gr_patch)

        this%ar_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_AR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 autotrophic respiration (MR + GR)', &
             ptr_patch=this%ar_patch)

        this%rr_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_RR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 root respiration (fine root MR + total root GR)', &
             ptr_patch=this%rr_patch)

        this%npp_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_NPP', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 net primary production', &
             ptr_patch=this%npp_patch)

        this%agnpp_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_AGNPP', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 aboveground NPP', &
             ptr_patch=this%agnpp_patch)

        this%bgnpp_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_BGNPP', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 belowground NPP', &
             ptr_patch=this%bgnpp_patch)

        this%litfall_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_LITFALL', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 litterfall (leaves and fine roots)', &
             ptr_patch=this%litfall_patch, default='inactive')

        this%vegfire_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_VEGFIRE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 patch-level fire loss', &
             ptr_patch=this%vegfire_patch, default='inactive')

        this%fire_closs_patch(begp:endp) = spval
        call hist_addfld1d (fname='C14_PFT_FIRE_CLOSS', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 total patch-level fire C loss', &
             ptr_patch=this%fire_closs_patch)
     endif

     !-------------------------------
     ! C flux variables - native to column 
     !-------------------------------

     ! add history fields for all CLAMP CN variables

     if (carbon_type == 'c12') then

        if (hist_wrtch4diag) then
           this%fphr_col(begc:endc,1:nlevgrnd) = spval
           call hist_addfld_decomp (fname='FPHR'//trim(vr_suffix), units='unitless', type2d='levdcmp', &
                avgflag='A', long_name='fraction of potential HR due to N limitation', &
                ptr_col=this%fphr_col)
        end if

        this%cwdc_hr_col(begc:endc) = spval
        call hist_addfld1d (fname='CWDC_HR', units='gC/m^2/s', &
             avgflag='A', long_name='coarse woody debris C heterotrophic respiration', &
             ptr_col=this%cwdc_hr_col)

        this%cwdc_loss_col(begc:endc) = spval
        call hist_addfld1d (fname='CWDC_LOSS', units='gC/m^2/s', &
             avgflag='A', long_name='coarse woody debris C loss', &
             ptr_col=this%cwdc_loss_col)

        this%lithr_col(begc:endc) = spval
        call hist_addfld1d (fname='LITTERC_HR', units='gC/m^2/s', &
             avgflag='A', long_name='litter C heterotrophic respiration', &
             ptr_col=this%lithr_col)

        this%litterc_loss_col(begc:endc) = spval
        call hist_addfld1d (fname='LITTERC_LOSS', units='gC/m^2/s', &
             avgflag='A', long_name='litter C loss', &
             ptr_col=this%litterc_loss_col)

        this%somhr_col(begc:endc) = spval
        call hist_addfld1d (fname='SOILC_HR', units='gC/m^2/s', &
             avgflag='A', long_name='soil C heterotrophic respiration', &
             ptr_col=this%somhr_col)

        this%somhr_col(begc:endc) = spval
        call hist_addfld1d (fname='SOILC_LOSS', units='gC/m^2/s', &
             avgflag='A', long_name='soil C loss', &
             ptr_col=this%somhr_col)

        ! F. Li and S. Levis
        this%lf_conv_cflux_col(begc:endc) = spval
        call hist_addfld1d (fname='LF_CONV_CFLUX', units='gC/m^2/s', &
             avgflag='A', long_name='conversion carbon due to BET and BDT area decreasing', &
             ptr_col=this%lf_conv_cflux_col, default='inactive')   

        this%somc_fire_col(begc:endc) = spval
        call hist_addfld1d (fname='SOMC_FIRE', units='gC/m^2/s', &
             avgflag='A', long_name='C loss due to peat burning', &
             ptr_col=this%somc_fire_col, default='inactive')


        this%m_decomp_cpools_to_fire_col(begc:endc,:)      = spval
        this%m_decomp_cpools_to_fire_vr_col(begc:endc,:,:) = spval
        do k = 1, ndecomp_pools
           if ( decomp_cascade_con%is_litter(k) .or. decomp_cascade_con%is_cwd(k) ) then
              data1dptr => this%m_decomp_cpools_to_fire_col(:,k)
              fieldname = 'M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'C_TO_FIRE'
              longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' C fire loss'
              call hist_addfld1d (fname=fieldname, units='gC/m^2/s',  &
                   avgflag='A', long_name=longname, &
                   ptr_col=data1dptr, default='inactive')

              if ( nlevdecomp_full > 1 ) then
                 data2dptr => this%m_decomp_cpools_to_fire_vr_col(:,:,k)
                 fieldname = 'M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'C_TO_FIRE'//trim(vr_suffix)
                 longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' C fire loss'
                 call hist_addfld_decomp (fname=fieldname, units='gC/m^3/s', type2d='levdcmp', &
                      avgflag='A', long_name=longname, &
                      ptr_col=data2dptr, default='inactive')
              endif
           endif

           ! decomposition k
           data2dptr => this%decomp_k_col(:,:,k)
           fieldname = 'K_'//trim(decomp_cascade_con%decomp_pool_name_history(k))
           longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' potential loss coefficient'
           call hist_addfld_decomp (fname=fieldname, units='1/s',  type2d='levdcmp', &
                avgflag='A', long_name=longname, &
                ptr_col=data2dptr, default='inactive')
        end do

        if(.not. is_active_betr_bgc)then
          this%decomp_cascade_hr_col(begc:endc,:)             = spval
          this%decomp_cascade_hr_vr_col(begc:endc,:,:)        = spval
          this%decomp_cascade_ctransfer_col(begc:endc,:)      = spval
          this%decomp_cascade_ctransfer_vr_col(begc:endc,:,:) = spval
          do l = 1, ndecomp_cascade_transitions

           ! output the vertically integrated fluxes only as  default
           !-- HR fluxes (none from CWD)
           if ( .not. decomp_cascade_con%is_cwd(decomp_cascade_con%cascade_donor_pool(l)) ) then
              data1dptr => this%decomp_cascade_hr_col(:,l)
              ! check to see if there are multiple pathways that include respiration, and if so, note that in the history file
              ii = 0
              do jj = 1, ndecomp_cascade_transitions
                 if ( decomp_cascade_con%cascade_donor_pool(jj) == decomp_cascade_con%cascade_donor_pool(l) ) ii = ii+1
              end do
              if ( ii == 1 ) then
                 fieldname = trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'_HR'
              else
                 fieldname = trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'_HR_'//&
                      trim(decomp_cascade_con%decomp_pool_name_short(decomp_cascade_con%cascade_receiver_pool(l)))
              endif
              longname =  'Het. Resp. from '//&
                   trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))
              call hist_addfld1d (fname=fieldname, units='gC/m^2/s',  &
                   avgflag='A', long_name=longname, &
                   ptr_col=data1dptr)
           endif

           !-- transfer fluxes (none from terminal pool, if present)
           if ( decomp_cascade_con%cascade_receiver_pool(l) /= 0 ) then
              data1dptr => this%decomp_cascade_ctransfer_col(:,l)
              fieldname = trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'C_TO_'//&
                   trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))//'C'
              longname =  'decomp. of '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))//&
                   ' C to '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_receiver_pool(l)))//' C'
              call hist_addfld1d (fname=fieldname, units='gC/m^2/s', &
                   avgflag='A', long_name=longname, &
                   ptr_col=data1dptr)
           endif

           ! output the vertically resolved fluxes 
           if ( nlevdecomp_full > 1 ) then  
              !-- HR fluxes (none from CWD)
              if ( .not. decomp_cascade_con%is_cwd(decomp_cascade_con%cascade_donor_pool(l)) ) then
                 data2dptr => this%decomp_cascade_hr_vr_col(:,:,l)
                 ! check to see if there are multiple pathways that include respiration, and if so, note that in the history file
                 ii = 0
                 do jj = 1, ndecomp_cascade_transitions
                    if ( decomp_cascade_con%cascade_donor_pool(jj) == decomp_cascade_con%cascade_donor_pool(l) ) ii = ii+1
                 end do
                 if ( ii == 1 ) then
                    fieldname = &
                         trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))&
                         //'_HR'//trim(vr_suffix)
                 else
                    fieldname = &
                         trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'_HR_'//&
                         trim(decomp_cascade_con%decomp_pool_name_short(decomp_cascade_con%cascade_receiver_pool(l)))&
                         //trim(vr_suffix)
                 endif
                 longname =  'Het. Resp. from '//&
                      trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))
                 call hist_addfld_decomp (fname=fieldname, units='gC/m^3/s',  type2d='levdcmp', &
                      avgflag='A', long_name=longname, &
                      ptr_col=data2dptr, default='inactive')
              endif

              !-- transfer fluxes (none from terminal pool, if present)
              if ( decomp_cascade_con%cascade_receiver_pool(l) /= 0 ) then
                 data2dptr => this%decomp_cascade_ctransfer_vr_col(:,:,l)
                 fieldname = trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'C_TO_'//&
                      trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))&
                      //'C'//trim(vr_suffix)
                 longname =  'decomp. of '//&
                      trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))//&
                      ' C to '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_receiver_pool(l)))//' C'
                 call hist_addfld_decomp (fname=fieldname, units='gC/m^3/s',  type2d='levdcmp', &
                      avgflag='A', long_name=longname, &
                      ptr_col=data2dptr, default='inactive')
              endif
           end if

          end do
        endif
        
        this%t_scalar_col(begc:endc,:) = spval
        call hist_addfld_decomp (fname='T_SCALAR', units='unitless',  type2d='levdcmp', &
             avgflag='A', long_name='temperature inhibition of decomposition', &
             ptr_col=this%t_scalar_col)

        this%w_scalar_col(begc:endc,:) = spval
        call hist_addfld_decomp (fname='W_SCALAR', units='unitless',  type2d='levdcmp', &
             avgflag='A', long_name='Moisture (dryness) inhibition of decomposition', &
             ptr_col=this%w_scalar_col)

        this%o_scalar_col(begc:endc,:) = spval
        call hist_addfld_decomp (fname='O_SCALAR', units='unitless', type2d='levdcmp', &
             avgflag='A', long_name='fraction by which decomposition is reduced due to anoxia', &
             ptr_col=this%o_scalar_col)

        this%som_c_leached_col(begc:endc) = spval
        call hist_addfld1d (fname='SOM_C_LEACHED', units='gC/m^2/s', &
             avgflag='A', long_name='total flux of C from SOM pools due to leaching', &
             ptr_col=this%som_c_leached_col)!, default='inactive')

        if(.not. is_active_betr_bgc)then     
          this%decomp_cpools_leached_col(begc:endc,:) = spval
          this%decomp_cpools_transport_tendency_col(begc:endc,:,:) = spval
          do k = 1, ndecomp_pools
           if ( .not. decomp_cascade_con%is_cwd(k) ) then
              data1dptr => this%decomp_cpools_leached_col(:,k)
              fieldname = 'M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'C_TO_LEACHING'
              longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' C leaching loss'
              call hist_addfld1d (fname=fieldname, units='gC/m^2/s', &
                   avgflag='A', long_name=longname, &
                   ptr_col=data1dptr)!, default='inactive')

              data2dptr => this%decomp_cpools_transport_tendency_col(:,:,k)
              fieldname = trim(decomp_cascade_con%decomp_pool_name_history(k))//'C_TNDNCY_VERT_TRANSPORT'
              longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' C tendency due to vertical transport'
              call hist_addfld_decomp (fname=fieldname, units='gC/m^3/s',  type2d='levdcmp', &
                   avgflag='A', long_name=longname, &
                   ptr_col=data2dptr, default='inactive')
           endif
          end do
        endif
        this%lithr_col(begc:endc) = spval
        call hist_addfld1d (fname='LITHR', units='gC/m^2/s', &
             avgflag='A', long_name='litter heterotrophic respiration', &
             ptr_col=this%lithr_col)

        this%somhr_col(begc:endc) = spval
        call hist_addfld1d (fname='SOMHR', units='gC/m^2/s', &
             avgflag='A', long_name='soil organic matter heterotrophic respiration', &
             ptr_col=this%somhr_col)

        if ( nlevdecomp_full > 1 ) then
           this%hr_vr_col(begc:endc,:) = spval
           call hist_addfld2d (fname='HR_vr', units='gC/m^3/s', type2d='levdcmp', &
                avgflag='A', long_name='total vertically resolved heterotrophic respiration', &
                ptr_col=this%hr_vr_col)

           ! pflotran
           this%f_co2_soil_vr_col(begc:endc,:) = spval
           call hist_addfld2d (fname='F_CO2_SOIL_vr', units='gC/m^3/s', type2d='levdcmp', &
                avgflag='A', long_name='total vertically resolved soil-atm. CO2 exchange', &
                ptr_col=this%f_co2_soil_vr_col)
        endif

        this%hr_col(begc:endc) = spval
        call hist_addfld1d (fname='HR', units='gC/m^2/s', &
             avgflag='A', long_name='total heterotrophic respiration', &
             ptr_col=this%hr_col)

        !pflotran
        this%f_co2_soil_col(begc:endc) = spval
        call hist_addfld1d (fname='F_CO2_SOIL', units='gC/m^2/s', &
             avgflag='A', long_name='total soil-atm. CO2 exchange', &
             ptr_col=this%f_co2_soil_col)

        this%sr_col(begc:endc) = spval
        call hist_addfld1d (fname='SR', units='gC/m^2/s', &
             avgflag='A', long_name='total soil respiration (HR + root resp)', &
             ptr_col=this%sr_col)

        this%er_col(begc:endc) = spval
        call hist_addfld1d (fname='ER', units='gC/m^2/s', &
             avgflag='A', long_name='total ecosystem respiration, autotrophic + heterotrophic', &
             ptr_col=this%er_col)

        this%litfire_col(begc:endc) = spval
        call hist_addfld1d (fname='LITFIRE', units='gC/m^2/s', &
             avgflag='A', long_name='litter fire losses', &
             ptr_col=this%litfire_col, default='inactive')

        this%somfire_col(begc:endc) = spval
        call hist_addfld1d (fname='SOMFIRE', units='gC/m^2/s', &
             avgflag='A', long_name='soil organic matter fire losses', &
             ptr_col=this%somfire_col, default='inactive')

        this%totfire_col(begc:endc) = spval
        call hist_addfld1d (fname='TOTFIRE', units='gC/m^2/s', &
             avgflag='A', long_name='total ecosystem fire losses', &
             ptr_col=this%totfire_col, default='inactive')

        this%nep_col(begc:endc) = spval
        call hist_addfld1d (fname='NEP', units='gC/m^2/s', &
             avgflag='A', long_name='net ecosystem production, excludes fire, landuse, and harvest flux, positive for sink', &
             ptr_col=this%nep_col)

        this%nbp_col(begc:endc) = spval
        call hist_addfld1d (fname='NBP', units='gC/m^2/s', &
             avgflag='A', long_name='net biome production, includes fire, landuse, and harvest flux, positive for sink', &
             ptr_col=this%nbp_col)

        this%nee_col(begc:endc) = spval
        call hist_addfld1d (fname='NEE', units='gC/m^2/s', &
             avgflag='A', long_name='net ecosystem exchange of carbon, includes fire, landuse,'&
             //' harvest, and hrv_xsmrpool flux, positive for source', &
             ptr_col=this%nee_col)

        this%fire_closs_col(begc:endc) = spval
        call hist_addfld1d (fname='COL_FIRE_CLOSS', units='gC/m^2/s', &
             avgflag='A', long_name='total column-level fire C loss for non-peat fires outside land-type converted region', &
             ptr_col=this%fire_closs_col, default='inactive')

        this%dwt_seedc_to_leaf_col(begc:endc) = spval
        call hist_addfld1d (fname='DWT_SEEDC_TO_LEAF', units='gC/m^2/s', &
             avgflag='A', long_name='seed source to patch-level leaf', &
             ptr_col=this%dwt_seedc_to_leaf_col, default='inactive')

        this%dwt_seedc_to_deadstem_col(begc:endc) = spval
        call hist_addfld1d (fname='DWT_SEEDC_TO_DEADSTEM', units='gC/m^2/s', &
             avgflag='A', long_name='seed source to patch-level deadstem', &
             ptr_col=this%dwt_seedc_to_deadstem_col, default='inactive')

        this%dwt_conv_cflux_col(begc:endc) = spval
        call hist_addfld1d (fname='DWT_CONV_CFLUX', units='gC/m^2/s', &
             avgflag='A', long_name='conversion C flux (immediate loss to atm)', &
             ptr_col=this%dwt_conv_cflux_col, default='inactive')

        this%dwt_prod10c_gain_col(begc:endc) = spval
        call hist_addfld1d (fname='DWT_PROD10C_GAIN', units='gC/m^2/s', &
             avgflag='A', long_name='landcover change-driven addition to 10-yr wood product pool', &
             ptr_col=this%dwt_prod10c_gain_col, default='inactive')

        this%prod10c_loss_col(begc:endc) = spval
        call hist_addfld1d (fname='PROD10C_LOSS', units='gC/m^2/s', &
             avgflag='A', long_name='loss from 10-yr wood product pool', &
             ptr_col=this%prod10c_loss_col, default='inactive')

        this%dwt_prod100c_gain_col(begc:endc) = spval
        call hist_addfld1d (fname='DWT_PROD100C_GAIN', units='gC/m^2/s', &
             avgflag='A', long_name='landcover change-driven addition to 100-yr wood product pool', &
             ptr_col=this%dwt_prod100c_gain_col, default='inactive')

        this%prod100c_loss_col(begc:endc) = spval
        call hist_addfld1d (fname='PROD100C_LOSS', units='gC/m^2/s', &
             avgflag='A', long_name='loss from 100-yr wood product pool', &
             ptr_col=this%prod100c_loss_col, default='inactive')

        this%prod1c_loss_col(begc:endc) = spval
        call hist_addfld1d (fname='PROD1C_LOSS', units='gC/m^2/s', &
             avgflag='A', long_name='loss from 1-yr crop product pool', &
             ptr_col=this%prod1c_loss_col, default='inactive')

        this%dwt_frootc_to_litr_met_c_col(begc:endc,:) = spval
        call hist_addfld_decomp (fname='DWT_FROOTC_TO_LITR_MET_C', units='gC/m^2/s',  type2d='levdcmp', &
             avgflag='A', long_name='fine root to litter due to landcover change', &
             ptr_col=this%dwt_frootc_to_litr_met_c_col, default='inactive')

        this%dwt_frootc_to_litr_cel_c_col(begc:endc,:) = spval
        call hist_addfld_decomp (fname='DWT_FROOTC_TO_LITR_CEL_C', units='gC/m^2/s',  type2d='levdcmp', &
             avgflag='A', long_name='fine root to litter due to landcover change', &
             ptr_col=this%dwt_frootc_to_litr_cel_c_col, default='inactive')

        this%dwt_frootc_to_litr_lig_c_col(begc:endc,:) = spval
        call hist_addfld_decomp (fname='DWT_FROOTC_TO_LITR_LIG_C', units='gC/m^2/s',  type2d='levdcmp', &
             avgflag='A', long_name='fine root to litter due to landcover change', &
             ptr_col=this%dwt_frootc_to_litr_lig_c_col, default='inactive')

        this%dwt_livecrootc_to_cwdc_col(begc:endc,:) = spval
        call hist_addfld_decomp (fname='DWT_LIVECROOTC_TO_CWDC', units='gC/m^2/s',  type2d='levdcmp', &
             avgflag='A', long_name='live coarse root to CWD due to landcover change', &
             ptr_col=this%dwt_livecrootc_to_cwdc_col, default='inactive')

        this%dwt_deadcrootc_to_cwdc_col(begc:endc,:) = spval
        call hist_addfld_decomp (fname='DWT_DEADCROOTC_TO_CWDC', units='gC/m^2/s',  type2d='levdcmp', &
             avgflag='A', long_name='dead coarse root to CWD due to landcover change', &
             ptr_col=this%dwt_deadcrootc_to_cwdc_col, default='inactive')

        this%dwt_closs_col(begc:endc) = spval
        call hist_addfld1d (fname='DWT_CLOSS', units='gC/m^2/s', &
             avgflag='A', long_name='total carbon loss from land cover conversion', &
             ptr_col=this%dwt_closs_col, default='inactive')

        this%product_closs_col(begc:endc) = spval
        call hist_addfld1d (fname='PRODUCT_CLOSS', units='gC/m^2/s', &
             avgflag='A', long_name='total carbon loss from wood product pools', &
             ptr_col=this%product_closs_col, default='inactive')

        this%landuseflux_col(begc:endc) = spval
        call hist_addfld1d (fname='LAND_USE_FLUX', units='gC/m^2/s', &
             avgflag='A', long_name='total C emitted from land cover conversion and wood product pools', &
             ptr_col=this%landuseflux_col)

        this%landuptake_col(begc:endc) = spval
        call hist_addfld1d (fname='LAND_UPTAKE', units='gC/m^2/s', &
             avgflag='A', long_name='NEE minus LAND_USE_FLUX, negative for update', &
             ptr_col=this%landuptake_col)

        this%annsum_npp_patch(begp:endp) = spval
        call hist_addfld1d (fname='ANNSUM_NPP', units='gC/m^2/yr', &
             avgflag='A', long_name='annual sum of NPP', &
             ptr_patch=this%annsum_npp_patch, default='inactive')

        this%annsum_npp_col(begc:endc) = spval
        call hist_addfld1d (fname='CANNSUM_NPP', units='gC/m^2/s', &
             avgflag='A', long_name='annual sum of column-level NPP', &
             ptr_col=this%annsum_npp_col, default='inactive')

     end if

     ctag=get_carbontag(carbon_type)
     do k = 1, ndecomp_pools
       this%bgc_cpool_ext_inputs_vr_col(begc:endc, :, k) = spval    
       data2dptr => this%bgc_cpool_ext_inputs_vr_col(:,:,k)
       fieldname='BGC_'//trim(ctag)//'POOL_EINPUT_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'_vr'
       longname=trim(ctag)//' input to '//trim(decomp_cascade_con%decomp_pool_name_history(k))
       call hist_addfld_decomp (fname=fieldname, units='g'//ctag//'/m^3',  type2d='levdcmp', &
         avgflag='A', long_name=longname, &
         ptr_col=data2dptr, default='inactive')

       this%bgc_cpool_ext_loss_vr_col(begc:endc, :, k) = spval    
       data2dptr => this%bgc_cpool_ext_loss_vr_col(:,:,k)
       fieldname='BGC_'//trim(ctag)//'POOL_ELOSS_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'_vr'
       longname=trim(ctag)//' loss of '//trim(decomp_cascade_con%decomp_pool_name_history(k))
       call hist_addfld_decomp (fname=fieldname, units='g'//ctag//'/m^3',  type2d='levdcmp', &
         avgflag='A', long_name=longname, &
         ptr_col=data2dptr, default='inactive')
         
     enddo

     !-------------------------------
     ! C13 flux variables - native to column 
     !-------------------------------

     if ( carbon_type == 'c13' ) then

        this%m_decomp_cpools_to_fire_col(begc:endc,:) = spval
        this%m_decomp_cpools_to_fire_vr_col(begc:endc,:,:) = spval
        do k = 1, ndecomp_pools
           if ( decomp_cascade_con%is_litter(k) .or. decomp_cascade_con%is_cwd(k) ) then
              data1dptr => this%m_decomp_cpools_to_fire_col(:,k)
              fieldname = 'C13_M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'C_TO_FIRE'
              longname =  'C13 '//trim(decomp_cascade_con%decomp_pool_name_long(k))//' C fire loss'
              call hist_addfld1d (fname=fieldname, units='gC13/m^2',  &
                   avgflag='A', long_name=longname, &
                   ptr_col=data1dptr, default='inactive')

              if ( nlevdecomp_full > 1 ) then
                 data2dptr => this%m_decomp_cpools_to_fire_vr_col(:,:,k)
                 fieldname = 'C13_M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'C_TO_FIRE'//trim(vr_suffix)
                 longname =  'C13 '//trim(decomp_cascade_con%decomp_pool_name_long(k))//' C fire loss'
                 call hist_addfld_decomp (fname=fieldname, units='gC13/m^3',  type2d='levdcmp', &
                      avgflag='A', long_name=longname, &
                      ptr_col=data2dptr, default='inactive')
              end if
           endif
        end do
        if(.not. is_active_betr_bgc)then
          this%decomp_cascade_hr_col(begc:endc,:)             = spval
          this%decomp_cascade_hr_vr_col(begc:endc,:,:)        = spval
          this%decomp_cascade_ctransfer_col(begc:endc,:)      = spval
          this%decomp_cascade_ctransfer_vr_col(begc:endc,:,:) = spval
          do l = 1, ndecomp_cascade_transitions
           !-- HR fluxes (none from CWD)
           if ( .not. decomp_cascade_con%is_cwd(decomp_cascade_con%cascade_donor_pool(l)) ) then
              data2dptr => this%decomp_cascade_hr_vr_col(:,:,l)
              ! check to see if there are multiple pathways that include respiration, and if so, note that in the history file
              ii = 0
              do jj = 1, ndecomp_cascade_transitions
                 if ( decomp_cascade_con%cascade_donor_pool(jj) == decomp_cascade_con%cascade_donor_pool(l) ) ii = ii+1
              end do
              if ( ii == 1 ) then
                 fieldname = 'C13_'//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))&
                      //'_HR'//trim(vr_suffix)
              else
                 fieldname = 'C13_'//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))&
                      //'_HR_'//&
                      trim(decomp_cascade_con%decomp_pool_name_short(decomp_cascade_con%cascade_receiver_pool(l)))//&
                      trim(vr_suffix)
              endif
              longname =  'C13 Het. Resp. from '&
                   //trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))
              call hist_addfld_decomp (fname=fieldname, units='gC13/m^3',  type2d='levdcmp', &
                   avgflag='A', long_name=longname, &
                   ptr_col=data2dptr, default='inactive')
           endif
           !-- transfer fluxes (none from terminal pool, if present)
           if ( decomp_cascade_con%cascade_receiver_pool(l) /= 0 ) then
              data2dptr => this%decomp_cascade_ctransfer_vr_col(:,:,l)
              fieldname = 'C13_'//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))&
                   //'C_TO_'//&
                   trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))&
                   //'C'//trim(vr_suffix)
              longname =  'C13 decomp. of '&
                   //trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))&
                   //' C to '//&
                   trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_receiver_pool(l)))//' C'
              call hist_addfld_decomp (fname=fieldname, units='gC13/m^3',  type2d='levdcmp', &
                   avgflag='A', long_name=longname, &
                   ptr_col=data2dptr, default='inactive')
           endif
          end do
        endif
        
        this%lithr_col(begc:endc) = spval
        call hist_addfld1d (fname='C13_LITHR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 fine root C litterfall to litter 3 C', &
             ptr_col=this%lithr_col)

        this%somhr_col(begc:endc) = spval
        call hist_addfld1d (fname='C13_SOMHR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 soil organic matter heterotrophic respiration', &
             ptr_col=this%somhr_col)

        this%hr_col(begc:endc) = spval
        call hist_addfld1d (fname='C13_HR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 total heterotrophic respiration', &
             ptr_col=this%hr_col)

        this%sr_col(begc:endc) = spval
        call hist_addfld1d (fname='C13_SR', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 total soil respiration (HR + root resp)', &
             ptr_col=this%sr_col)

        this%er_col(begc:endc) = spval
        call hist_addfld1d (fname='C13_ER', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 total ecosystem respiration, autotrophic + heterotrophic', &
             ptr_col=this%er_col)

        this%litfire_col(begc:endc) = spval
        call hist_addfld1d (fname='C13_LITFIRE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 litter fire losses', &
             ptr_col=this%litfire_col, default='inactive')

        this%somfire_col(begc:endc) = spval
        call hist_addfld1d (fname='C13_SOMFIRE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 soil organic matter fire losses', &
             ptr_col=this%somfire_col, default='inactive')

        this%totfire_col(begc:endc) = spval
        call hist_addfld1d (fname='C13_TOTFIRE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 total ecosystem fire losses', &
             ptr_col=this%totfire_col, default='inactive')

        this%nep_col(begc:endc) = spval
        call hist_addfld1d (fname='C13_NEP', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 net ecosystem production, excludes fire flux, positive for sink', &
             ptr_col=this%nep_col)

        this%nee_col(begc:endc) = spval
        call hist_addfld1d (fname='C13_NEE', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 net ecosystem exchange of carbon, includes fire flux, positive for source', &
             ptr_col=this%nee_col)

        this%fire_closs_col(begc:endc) = spval
        call hist_addfld1d (fname='C13_COL_FIRE_CLOSS', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 total column-level fire C loss', &
             ptr_col=this%fire_closs_col)

        this%dwt_seedc_to_leaf_col(begc:endc) = spval
        call hist_addfld1d (fname='C13_DWT_SEEDC_TO_LEAF', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 seed source to patch-level leaf', &
             ptr_col=this%dwt_seedc_to_leaf_col)

        this%dwt_seedc_to_deadstem_col(begc:endc) = spval
        call hist_addfld1d (fname='C13_DWT_SEEDC_TO_DEADSTEM', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 seed source to patch-level deadstem', &
             ptr_col=this%dwt_seedc_to_deadstem_col)

        this%dwt_conv_cflux_col(begc:endc) = spval
        call hist_addfld1d (fname='C13_DWT_CONV_CFLUX', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 conversion C flux (immediate loss to atm)', &
             ptr_col=this%dwt_conv_cflux_col)

        this%dwt_prod10c_gain_col(begc:endc) = spval
        call hist_addfld1d (fname='C13_DWT_PROD10C_GAIN', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 addition to 10-yr wood product pool', &
             ptr_col=this%dwt_prod10c_gain_col)

        this%prod10c_loss_col(begc:endc) = spval
        call hist_addfld1d (fname='C13_PROD10C_LOSS', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 loss from 10-yr wood product pool', &
             ptr_col=this%prod10c_loss_col)

        this%dwt_prod100c_gain_col(begc:endc) = spval
        call hist_addfld1d (fname='C13_DWT_PROD100C_GAIN', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 addition to 100-yr wood product pool', &
             ptr_col=this%dwt_prod100c_gain_col)

        this%prod100c_loss_col(begc:endc) = spval
        call hist_addfld1d (fname='C13_PROD100C_LOSS', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 loss from 100-yr wood product pool', &
             ptr_col=this%prod100c_loss_col)

        this%prod1c_loss_col(begc:endc) = spval
        call hist_addfld1d (fname='C13_PROD1C_LOSS', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 loss from 1-yr crop product pool', &
             ptr_col=this%prod1c_loss_col)

        this%dwt_frootc_to_litr_met_c_col(begc:endc,:) = spval
        call hist_addfld_decomp (fname='C13_DWT_FROOTC_TO_LITR_MET_C', units='gC13/m^2/s',  type2d='levdcmp', &
             avgflag='A', long_name='C13 fine root to litter due to landcover change', &
             ptr_col=this%dwt_frootc_to_litr_met_c_col, default='inactive')

        this%dwt_frootc_to_litr_cel_c_col(begc:endc,:) = spval
        call hist_addfld_decomp (fname='C13_DWT_FROOTC_TO_LITR_CEL_C', units='gC13/m^2/s',  type2d='levdcmp', &
             avgflag='A', long_name='C13 fine root to litter due to landcover change', &
             ptr_col=this%dwt_frootc_to_litr_cel_c_col, default='inactive')

        this%dwt_frootc_to_litr_lig_c_col(begc:endc,:) = spval
        call hist_addfld_decomp (fname='C13_DWT_FROOTC_TO_LITR_LIG_C', units='gC13/m^2/s',  type2d='levdcmp', &
             avgflag='A', long_name='C13 fine root to litter due to landcover change', &
             ptr_col=this%dwt_frootc_to_litr_lig_c_col, default='inactive')

        this%dwt_livecrootc_to_cwdc_col(begc:endc,:) = spval
        call hist_addfld_decomp (fname='C13_DWT_LIVECROOTC_TO_CWDC', units='gC13/m^2/s',  type2d='levdcmp', &
             avgflag='A', long_name='C13 live coarse root to CWD due to landcover change', &
             ptr_col=this%dwt_livecrootc_to_cwdc_col, default='inactive')

        this%dwt_deadcrootc_to_cwdc_col(begc:endc,:) = spval
        call hist_addfld_decomp (fname='C13_DWT_DEADCROOTC_TO_CWDC', units='gC13/m^2/s',  type2d='levdcmp', &
             avgflag='A', long_name='C13 dead coarse root to CWD due to landcover change', &
             ptr_col=this%dwt_deadcrootc_to_cwdc_col, default='inactive')

        this%dwt_closs_col(begc:endc) = spval
        call hist_addfld1d (fname='C13_DWT_CLOSS', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 total carbon loss from land cover conversion', &
             ptr_col=this%dwt_closs_col)

        this%product_closs_col(begc:endc) = spval
        call hist_addfld1d (fname='C13_PRODUCT_CLOSS', units='gC13/m^2/s', &
             avgflag='A', long_name='C13 total carbon loss from wood product pools', &
             ptr_col=this%product_closs_col)
     endif

     !-------------------------------
     ! C14 flux variables - native to column 
     !-------------------------------

     if (carbon_type == 'c14') then

        this%m_decomp_cpools_to_fire_col(begc:endc,:)      = spval
        this%m_decomp_cpools_to_fire_vr_col(begc:endc,:,:) = spval
        do k = 1, ndecomp_pools
           if ( decomp_cascade_con%is_litter(k) .or. decomp_cascade_con%is_cwd(k) ) then
              data1dptr => this%m_decomp_cpools_to_fire_col(:,k)
              fieldname = 'C14_M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'C_TO_FIRE'
              longname =  'C14 '//trim(decomp_cascade_con%decomp_pool_name_long(k))//' C fire loss'
              call hist_addfld1d (fname=fieldname, units='gC14/m^2',  &
                   avgflag='A', long_name=longname, &
                   ptr_col=data1dptr, default='inactive')

              if ( nlevdecomp_full > 1 ) then
                 data2dptr => this%m_decomp_cpools_to_fire_vr_col(:,:,k)
                 fieldname = 'C14_M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'C_TO_FIRE'//trim(vr_suffix)
                 longname =  'C14 '//trim(decomp_cascade_con%decomp_pool_name_long(k))//' C fire loss'
                 call hist_addfld_decomp (fname=fieldname, units='gC14/m^3',  type2d='levdcmp', &
                      avgflag='A', long_name=longname, &
                      ptr_col=data2dptr, default='inactive')
              end if
           endif
        end do
        if(.not. is_active_betr_bgc)then
          this%decomp_cascade_hr_col(begc:endc,:)             = spval
          this%decomp_cascade_hr_vr_col(begc:endc,:,:)        = spval
          this%decomp_cascade_ctransfer_col(begc:endc,:)      = spval
          this%decomp_cascade_ctransfer_vr_col(begc:endc,:,:) = spval
          do l = 1, ndecomp_cascade_transitions
           !-- HR fluxes (none from CWD)
           if ( .not. decomp_cascade_con%is_cwd(decomp_cascade_con%cascade_donor_pool(l)) ) then
              data2dptr => this%decomp_cascade_hr_vr_col(:,:,l)
              ! check to see if there are multiple pathways that include respiration, and if so, note that in the history file
              ii = 0
              do jj = 1, ndecomp_cascade_transitions
                 if ( decomp_cascade_con%cascade_donor_pool(jj) == decomp_cascade_con%cascade_donor_pool(l) ) ii = ii+1
              end do
              if ( ii == 1 ) then
                 fieldname = 'C14_'//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))&
                      //'_HR'//trim(vr_suffix)
              else
                 fieldname = 'C14_'//&
                      trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))&
                      //'_HR_'//&
                      trim(decomp_cascade_con%decomp_pool_name_short(decomp_cascade_con%cascade_receiver_pool(l)))&
                      //trim(vr_suffix)
              endif
              longname =  'C14 Het. Resp. from '&
                   //trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))
              call hist_addfld_decomp (fname=fieldname, units='gC14/m^3',  type2d='levdcmp', &
                   avgflag='A', long_name=longname, &
                   ptr_col=data2dptr, default='inactive')
           endif
           !-- transfer fluxes (none from terminal pool, if present)
           if ( decomp_cascade_con%cascade_receiver_pool(l) /= 0 ) then
              data2dptr => this%decomp_cascade_ctransfer_vr_col(:,:,l)
              fieldname = 'C14_'//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))&
                   //'C_TO_'//&
                   trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))&
                   //'C'//trim(vr_suffix)
              longname =  'C14 decomp. of '&
                   //trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))//&
                   ' C to '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_receiver_pool(l)))//' C'
              call hist_addfld_decomp (fname=fieldname, units='gC14/m^3',  type2d='levdcmp', &
                   avgflag='A', long_name=longname, &
                   ptr_col=data2dptr, default='inactive')
           endif
          end do
        endif
        
        this%lithr_col(begc:endc) = spval
        call hist_addfld1d (fname='C14_LITHR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 fine root C litterfall to litter 3 C', &
             ptr_col=this%lithr_col)

        this%somhr_col(begc:endc) = spval
        call hist_addfld1d (fname='C14_SOMHR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 soil organic matter heterotrophic respiration', &
             ptr_col=this%somhr_col)

        this%hr_col(begc:endc) = spval
        call hist_addfld1d (fname='C14_HR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 total heterotrophic respiration', &
             ptr_col=this%hr_col)

        this%sr_col(begc:endc) = spval
        call hist_addfld1d (fname='C14_SR', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 total soil respiration (HR + root resp)', &
             ptr_col=this%sr_col)

        this%er_col(begc:endc) = spval
        call hist_addfld1d (fname='C14_ER', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 total ecosystem respiration, autotrophic + heterotrophic', &
             ptr_col=this%er_col)

        this%litfire_col(begc:endc) = spval
        call hist_addfld1d (fname='C14_LITFIRE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 litter fire losses', &
             ptr_col=this%litfire_col, default='inactive')

        this%somfire_col(begc:endc) = spval
        call hist_addfld1d (fname='C14_SOMFIRE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 soil organic matter fire losses', &
             ptr_col=this%somfire_col, default='inactive')

        this%totfire_col(begc:endc) = spval
        call hist_addfld1d (fname='C14_TOTFIRE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 total ecosystem fire losses', &
             ptr_col=this%totfire_col, default='inactive')

        this%nep_col(begc:endc) = spval
        call hist_addfld1d (fname='C14_NEP', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 net ecosystem production, excludes fire flux, positive for sink', &
             ptr_col=this%nep_col)

        this%nee_col(begc:endc) = spval
        call hist_addfld1d (fname='C14_NEE', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 net ecosystem exchange of carbon, includes fire flux, positive for source', &
             ptr_col=this%nee_col)

        this%fire_closs_col(begc:endc) = spval
        call hist_addfld1d (fname='C14_COL_FIRE_CLOSS', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 total column-level fire C loss', &
             ptr_col=this%fire_closs_col)

        this%dwt_seedc_to_leaf_col(begc:endc) = spval
        call hist_addfld1d (fname='C14_DWT_SEEDC_TO_LEAF', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 seed source to patch-level leaf', &
             ptr_col=this%dwt_seedc_to_leaf_col)

        this%dwt_seedc_to_deadstem_col(begc:endc) = spval
        call hist_addfld1d (fname='C14_DWT_SEEDC_TO_DEADSTEM', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 seed source to patch-level deadstem', &
             ptr_col=this%dwt_seedc_to_deadstem_col)

        this%dwt_conv_cflux_col(begc:endc) = spval
        call hist_addfld1d (fname='C14_DWT_CONV_CFLUX', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 conversion C flux (immediate loss to atm)', &
             ptr_col=this%dwt_conv_cflux_col)

        this%dwt_prod10c_gain_col(begc:endc) = spval
        call hist_addfld1d (fname='C14_DWT_PROD10C_GAIN', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 addition to 10-yr wood product pool', &
             ptr_col=this%dwt_prod10c_gain_col)

        this%prod10c_loss_col(begc:endc) = spval
        call hist_addfld1d (fname='C14_PROD10C_LOSS', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 loss from 10-yr wood product pool', &
             ptr_col=this%prod10c_loss_col)

        this%dwt_prod100c_gain_col(begc:endc) = spval
        call hist_addfld1d (fname='C14_DWT_PROD100C_GAIN', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 addition to 100-yr wood product pool', &
             ptr_col=this%dwt_prod100c_gain_col)

        this%prod100c_loss_col(begc:endc) = spval
        call hist_addfld1d (fname='C14_PROD100C_LOSS', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 loss from 100-yr wood product pool', &
             ptr_col=this%prod100c_loss_col)

        this%prod1c_loss_col(begc:endc) = spval
        call hist_addfld1d (fname='C14_PROD1C_LOSS', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 loss from 1-yr crop product pool', &
             ptr_col=this%prod1c_loss_col)

        this%dwt_frootc_to_litr_met_c_col(begc:endc,:) = spval
        call hist_addfld_decomp (fname='C14_DWT_FROOTC_TO_LITR_MET_C', units='gC14/m^2/s',  type2d='levdcmp', &
             avgflag='A', long_name='C14 fine root to litter due to landcover change', &
             ptr_col=this%dwt_frootc_to_litr_met_c_col, default='inactive')

        this%dwt_frootc_to_litr_cel_c_col(begc:endc,:) = spval
        call hist_addfld_decomp (fname='C14_DWT_FROOTC_TO_LITR_CEL_C', units='gC14/m^2/s',  type2d='levdcmp', &
             avgflag='A', long_name='C14 fine root to litter due to landcover change', &
             ptr_col=this%dwt_frootc_to_litr_cel_c_col, default='inactive')

        this%dwt_frootc_to_litr_lig_c_col(begc:endc,:) = spval
        call hist_addfld_decomp (fname='C14_DWT_FROOTC_TO_LITR_LIG_C', units='gC14/m^2/s',  type2d='levdcmp', &
             avgflag='A', long_name='C14 fine root to litter due to landcover change', &
             ptr_col=this%dwt_frootc_to_litr_lig_c_col, default='inactive')

        this%dwt_livecrootc_to_cwdc_col(begc:endc,:) = spval
        call hist_addfld_decomp (fname='C14_DWT_LIVECROOTC_TO_CWDC', units='gC14/m^2/s',  type2d='levdcmp', &
             avgflag='A', long_name='C14 live coarse root to CWD due to landcover change', &
             ptr_col=this%dwt_livecrootc_to_cwdc_col, default='inactive')

        this%dwt_deadcrootc_to_cwdc_col(begc:endc,:) = spval
        call hist_addfld_decomp (fname='C14_DWT_DEADCROOTC_TO_CWDC', units='gC14/m^2/s',  type2d='levdcmp', &
             avgflag='A', long_name='C14 dead coarse root to CWD due to landcover change', &
             ptr_col=this%dwt_deadcrootc_to_cwdc_col, default='inactive')

        this%dwt_closs_col(begc:endc) = spval
        call hist_addfld1d (fname='C14_DWT_CLOSS', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 total carbon loss from land cover conversion', &
             ptr_col=this%dwt_closs_col)

        this%product_closs_col(begc:endc) = spval
        call hist_addfld1d (fname='C14_PRODUCT_CLOSS', units='gC14/m^2/s', &
             avgflag='A', long_name='C14 total carbon loss from wood product pools', &
             ptr_col=this%product_closs_col)
     endif
     
     if (carbon_type == 'c13') then
        this%xsmrpool_c13ratio_patch(begp:endp) = spval
        call hist_addfld1d (fname='XSMRPOOL_C13RATIO', units='proportion', &
             avgflag='A', long_name='C13/C(12+13) ratio for xsmrpool', &
             ptr_patch=this%xsmrpool_c13ratio_patch, default='inactive')
     endif

   end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !ARGUMENTS:
    class(carbonflux_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: p, c, l, j
    integer :: fc                                        ! filter index
    integer :: num_special_col                           ! number of good values in special_col filter
    integer :: num_special_patch                         ! number of good values in special_patch filter
    integer :: special_col(bounds%endc-bounds%begc+1)    ! special landunit filter - columns
    integer :: special_patch(bounds%endp-bounds%begp+1)  ! special landunit filter - patches
    !-----------------------------------------------------------------------

    ! Set column filters

    num_special_col = 0
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%ifspecial(l)) then
          num_special_col = num_special_col + 1
          special_col(num_special_col) = c
       end if
    end do

    ! Set patch filters

    num_special_patch = 0
    do p = bounds%begp,bounds%endp
       l = pft%landunit(p)

       if (lun%ifspecial(l)) then
          num_special_patch = num_special_patch + 1
          special_patch(num_special_patch) = p
       end if
    end do

    do p = bounds%begp,bounds%endp
       l = pft%landunit(p)

       this%gpp_patch(p)                      = 0._r8
       this%gpp_before_downreg_patch(p)       = 0._r8

       if (lun%ifspecial(l)) then
          this%tempsum_npp_patch(p)           = spval
          this%annsum_npp_patch(p)            = spval
          this%availc_patch(p)                = spval
          this%xsmrpool_recover_patch(p)      = spval
          this%excess_cflux_patch(p)          = spval
          this%plant_calloc_patch(p)          = spval
          this%prev_leafc_to_litter_patch(p)  = spval
          this%prev_frootc_to_litter_patch(p) = spval
          if (use_cndv) then
             this%tempsum_litfall_patch(p)    = spval
             this%annsum_litfall_patch(p)     = spval
          end if
          if ( use_c13 ) then
             this%xsmrpool_c13ratio_patch(p)  = spval
          endif
       end if
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          this%tempsum_npp_patch(p)           = 0._r8
          this%annsum_npp_patch(p)            = 0._r8
          this%availc_patch(p)                = 0._r8
          this%xsmrpool_recover_patch(p)      = 0._r8
          this%excess_cflux_patch(p)          = 0._r8
          this%prev_leafc_to_litter_patch(p)  = 0._r8
          this%prev_frootc_to_litter_patch(p) = 0._r8
          if (use_cndv) then
             this%tempsum_litfall_patch(p)    = 0._r8
             this%annsum_litfall_patch(p)     = 0._r8
          end if
          this%plant_calloc_patch(p)          = 0._r8
       end if
    end do

    do c = bounds%begc, bounds%endc
       l = col%landunit(c)

       if (lun%ifspecial(l)) then
          this%annsum_npp_col(c) = spval
       end if

       this%fphr_col(c,nlevdecomp+1:nlevgrnd) = 0._r8 !used to be in ch4Mod
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          this%fphr_col(c,nlevdecomp+1:nlevgrnd) = 0._r8 
       else if (lun%itype(l) == istdlak .and. allowlakeprod) then
          this%fphr_col(c,:) = spval
       else  ! Inactive CH4 columns
          this%fphr_col(c,:) = spval
       end if

       ! also initialize dynamic landcover fluxes so that they have
       ! real values on first timestep, prior to calling pftdyn_cnbal
       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          this%lf_conv_cflux_col(c)         = 0._r8
          this%dwt_seedc_to_leaf_col(c)     = 0._r8
          this%dwt_seedc_to_deadstem_col(c) = 0._r8
          this%dwt_conv_cflux_col(c)        = 0._r8
          this%dwt_prod10c_gain_col(c)      = 0._r8
          this%dwt_prod100c_gain_col(c)     = 0._r8
          this%prod1c_loss_col(c)           = 0._r8
          this%prod10c_loss_col(c)          = 0._r8
          this%prod100c_loss_col(c)         = 0._r8
          do j = 1, nlevdecomp_full
             this%dwt_frootc_to_litr_met_c_col(c,j) = 0._r8
             this%dwt_frootc_to_litr_cel_c_col(c,j) = 0._r8
             this%dwt_frootc_to_litr_lig_c_col(c,j) = 0._r8
             this%dwt_livecrootc_to_cwdc_col(c,j)   = 0._r8
             this%dwt_deadcrootc_to_cwdc_col(c,j)   = 0._r8
          end do
          this%dwt_closs_col(c)  = 0._r8
          this%annsum_npp_col(c) = 0._r8   
       end if
    end do

    ! initialize fields for special filters

    do fc = 1,num_special_col
       c = special_col(fc)
       
       this%dwt_closs_col(c)   = 0._r8
       this%landuseflux_col(c) = 0._r8
       this%landuptake_col(c)  = 0._r8
    end do

    ! initialize fields for special filters

    call this%SetValues (&
         num_patch=num_special_patch, filter_patch=special_patch, value_patch=0._r8, &
         num_column=num_special_col, filter_column=special_col, value_column=0._r8)

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine Restart ( this, bounds, ncid, flag )
    !
    ! !DESCRIPTION: 
    ! Read/write CN restart data for carbon state
    !
    ! !USES:
    use shr_infnan_mod   , only : isnan => shr_infnan_isnan, nan => shr_infnan_nan, assignment(=)
    use clm_time_manager , only : is_restart
    use clm_varcon       , only : c13ratio, c14ratio
    use clm_varctl       , only : use_lch4, use_betr
    use restUtilMod
    use ncdio_pio

    ! pflotran
!    use clm_varctl       , only : use_pflotran, pf_cmode, use_vertsoilc
    !
    ! !ARGUMENTS:
    class (carbonflux_type) :: this
    type(bounds_type) , intent(in)    :: bounds 
    type(file_desc_t) , intent(inout) :: ncid   ! netcdf id
    character(len=*)  , intent(in)    :: flag   !'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer :: j,c ! indices
    logical :: readvar      ! determine if variable is on initial file

    ! pflotran
    integer :: k
    real(r8), pointer :: ptr2d(:,:) ! temp. pointers for slicing larger arrays
    real(r8), pointer :: ptr1d(:)   ! temp. pointers for slicing larger arrays
    character(len=128) :: varname   ! temporary
    !------------------------------------------------------------------------

    !-------------------------------
    ! Prognostic crop variables
    !-------------------------------

    if (crop_prog) then

       call restartvar(ncid=ncid, flag=flag,  varname='grainc_xfer_to_grainc', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='grain C growth from storage', units='gC/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%grainc_xfer_to_grainc_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='livestemc_to_litter', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='live stem C litterfall', units='gC/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%livestemc_to_litter_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='grainc_to_food', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='grain C to food', units='gC/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%grainc_to_food_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='cpool_to_grainc', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='allocation to grain C', units='gC/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%cpool_to_grainc_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='cpool_to_grainc_storage', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='allocation to grain C storage', units='gC/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%cpool_to_grainc_storage_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='cpool_grain_gr', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='grain growth respiration', units='gC/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%cpool_grain_gr_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='cpool_grain_storage_gr', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='grain growth respiration to storage', units='gC/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%cpool_grain_storage_gr_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='transfer_grain_gr', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='grain growth respiration from storage', units='gC/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%transfer_grain_gr_patch)

       call restartvar(ncid=ncid, flag=flag,  varname='grainc_storage_to_xfer', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='grain C shift storage to transfer', units='gC/m2/s', &
            interpinic_flag='interp', readvar=readvar, data=this%grainc_storage_to_xfer_patch)

    end if

    if (use_lch4 .or. use_betr) then
       call restartvar(ncid=ncid, flag=flag, varname='tempavg_agnpp', xtype=ncd_double,  &
            dim1name='pft',&
            long_name='Temp. Average AGNPP',units='gC/m^2/s', &
            readvar=readvar, interpinic_flag='interp', data=this%tempavg_agnpp_patch)
       
       call restartvar(ncid=ncid, flag=flag, varname='tempavg_bgnpp', xtype=ncd_double,  &
            dim1name='pft',&
            long_name='Temp. Average BGNPP',units='gC/m^2/s', &
            readvar=readvar, interpinic_flag='interp', data=this%tempavg_bgnpp_patch)
       
       call restartvar(ncid=ncid, flag=flag, varname='annavg_agnpp', xtype=ncd_double,  &
            dim1name='pft',&
            long_name='Ann. Average AGNPP',units='gC/m^2/s', &
            readvar=readvar, interpinic_flag='interp', data=this%annavg_agnpp_patch)
       
       call restartvar(ncid=ncid, flag=flag, varname='annavg_bgnpp', xtype=ncd_double,  &
            dim1name='pft',&
            long_name='Ann. Average BGNPP',units='gC/m^2/s', &
            readvar=readvar, interpinic_flag='interp', data=this%annavg_bgnpp_patch)
    end if

    call restartvar(ncid=ncid, flag=flag, varname='gpp_pepv', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%gpp_before_downreg_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='availc', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%availc_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='xsmrpool_recover', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%xsmrpool_recover_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='plant_calloc', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%plant_calloc_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='excess_cflux', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%excess_cflux_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='prev_leafc_to_litter', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%prev_leafc_to_litter_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='prev_frootc_to_litter', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%prev_frootc_to_litter_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='tempsum_npp', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%tempsum_npp_patch) 
 
    call restartvar(ncid=ncid, flag=flag, varname='annsum_npp', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%annsum_npp_patch) 

    call restartvar(ncid=ncid, flag=flag, varname='col_lag_npp', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%lag_npp_col) 

    call restartvar(ncid=ncid, flag=flag, varname='cannsum_npp', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%annsum_npp_col) 

    if (use_cndv) then
       call restartvar(ncid=ncid, flag=flag, varname='tempsum_litfall', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%tempsum_litfall_patch)

       call restartvar(ncid=ncid, flag=flag, varname='annsum_litfall', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=this%annsum_litfall_patch)
    end if

    ! clm_bgc_interface & pflotran
    !------------------------------------------------------------------------
    if (use_pflotran .and. pf_cmode) then
       ! externalc_to_decomp_npools_col
       do k = 1, ndecomp_pools
          varname=trim(decomp_cascade_con%decomp_pool_name_restart(k))//'external_c'
          if (use_vertsoilc) then
             ptr2d => this%externalc_to_decomp_cpools_col(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_vr",  &
                  xtype=ncd_double, dim1name='column', dim2name='levgrnd', switchdim=.true., &
                  long_name='net soil organic C adding/removal/transport', &
                  units='gC/m3/s', fill_value=spval, &
                  interpinic_flag='interp', readvar=readvar, data=ptr2d)
          else
             ptr1d => this%externalc_to_decomp_cpools_col(:,1,k) ! nlevdecomp = 1; so treat as 1D variable
             call restartvar(ncid=ncid, flag=flag, varname=varname, &
                  xtype=ncd_double, dim1name='column', &
                  long_name='net soil organic C adding/removal/transport', &
                  units='gC/m3/s', fill_value=spval, &
                  interpinic_flag='interp' , readvar=readvar, data=ptr1d)
          end if
          if (flag=='read' .and. .not. readvar) then
          !   call endrun(msg='ERROR:: '//trim(varname)//' is required on an initialization dataset'//&
          !        errMsg(__FILE__, __LINE__))
             this%externalc_to_decomp_cpools_col(:,:,k) = 0._r8
          end if
       end do
    end if
    !------------------------------------------------------------------------

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine SetValues ( this, &
       num_patch, filter_patch, value_patch, &
       num_column, filter_column, value_column)
    !
    ! !DESCRIPTION:
    ! Set carbon state fluxes
    !
    ! !ARGUMENTS:
    class (carbonflux_type) :: this
    integer , intent(in) :: num_patch
    integer , intent(in) :: filter_patch(:)
    real(r8), intent(in) :: value_patch
    integer , intent(in) :: num_column
    integer , intent(in) :: filter_column(:)
    real(r8), intent(in) :: value_column
    !
    ! !LOCAL VARIABLES:
    integer :: fi,i     ! loop index
    integer :: j,k,l    ! indices
    !------------------------------------------------------------------------

    do fi = 1,num_patch
       i = filter_patch(fi)

       this%m_leafc_to_litter_patch(i)                   = value_patch
       this%m_frootc_to_litter_patch(i)                  = value_patch
       this%m_leafc_storage_to_litter_patch(i)           = value_patch
       this%m_frootc_storage_to_litter_patch(i)          = value_patch
       this%m_livestemc_storage_to_litter_patch(i)       = value_patch
       this%m_deadstemc_storage_to_litter_patch(i)       = value_patch
       this%m_livecrootc_storage_to_litter_patch(i)      = value_patch
       this%m_deadcrootc_storage_to_litter_patch(i)      = value_patch
       this%m_leafc_xfer_to_litter_patch(i)              = value_patch
       this%m_frootc_xfer_to_litter_patch(i)             = value_patch
       this%m_livestemc_xfer_to_litter_patch(i)          = value_patch
       this%m_deadstemc_xfer_to_litter_patch(i)          = value_patch
       this%m_livecrootc_xfer_to_litter_patch(i)         = value_patch
       this%m_deadcrootc_xfer_to_litter_patch(i)         = value_patch
       this%m_livestemc_to_litter_patch(i)               = value_patch
       this%m_deadstemc_to_litter_patch(i)               = value_patch
       this%m_livecrootc_to_litter_patch(i)              = value_patch
       this%m_deadcrootc_to_litter_patch(i)              = value_patch
       this%m_gresp_storage_to_litter_patch(i)           = value_patch
       this%m_gresp_xfer_to_litter_patch(i)              = value_patch
       this%hrv_leafc_to_litter_patch(i)                 = value_patch             
       this%hrv_leafc_storage_to_litter_patch(i)         = value_patch     
       this%hrv_leafc_xfer_to_litter_patch(i)            = value_patch        
       this%hrv_frootc_to_litter_patch(i)                = value_patch            
       this%hrv_frootc_storage_to_litter_patch(i)        = value_patch    
       this%hrv_frootc_xfer_to_litter_patch(i)           = value_patch       
       this%hrv_livestemc_to_litter_patch(i)             = value_patch         
       this%hrv_livestemc_storage_to_litter_patch(i)     = value_patch 
       this%hrv_livestemc_xfer_to_litter_patch(i)        = value_patch    
       this%hrv_deadstemc_to_prod10c_patch(i)            = value_patch        
       this%hrv_deadstemc_to_prod100c_patch(i)           = value_patch       
       this%hrv_leafc_to_prod1c_patch(i)                 = value_patch
       this%hrv_livestemc_to_prod1c_patch(i)             = value_patch
       this%hrv_grainc_to_prod1c_patch(i)                = value_patch
       this%hrv_cropc_to_prod1c_patch(i)                 = value_patch
       this%hrv_deadstemc_storage_to_litter_patch(i)     = value_patch 
       this%hrv_deadstemc_xfer_to_litter_patch(i)        = value_patch    
       this%hrv_livecrootc_to_litter_patch(i)            = value_patch        
       this%hrv_livecrootc_storage_to_litter_patch(i)    = value_patch
       this%hrv_livecrootc_xfer_to_litter_patch(i)       = value_patch   
       this%hrv_deadcrootc_to_litter_patch(i)            = value_patch        
       this%hrv_deadcrootc_storage_to_litter_patch(i)    = value_patch
       this%hrv_deadcrootc_xfer_to_litter_patch(i)       = value_patch   
       this%hrv_gresp_storage_to_litter_patch(i)         = value_patch     
       this%hrv_gresp_xfer_to_litter_patch(i)            = value_patch        
       this%hrv_xsmrpool_to_atm_patch(i)                 = value_patch

       this%m_leafc_to_fire_patch(i)                     = value_patch
       this%m_leafc_storage_to_fire_patch(i)             = value_patch
       this%m_leafc_xfer_to_fire_patch(i)                = value_patch
       this%m_livestemc_to_fire_patch(i)                 = value_patch
       this%m_livestemc_storage_to_fire_patch(i)         = value_patch
       this%m_livestemc_xfer_to_fire_patch(i)            = value_patch
       this%m_deadstemc_to_fire_patch(i)                 = value_patch
       this%m_deadstemc_storage_to_fire_patch(i)         = value_patch
       this%m_deadstemc_xfer_to_fire_patch(i)            = value_patch
       this%m_frootc_to_fire_patch(i)                    = value_patch
       this%m_frootc_storage_to_fire_patch(i)            = value_patch
       this%m_frootc_xfer_to_fire_patch(i)               = value_patch
       this%m_livecrootc_to_fire_patch(i)                = value_patch
       this%m_livecrootc_storage_to_fire_patch(i)        = value_patch
       this%m_livecrootc_xfer_to_fire_patch(i)           = value_patch
       this%m_deadcrootc_to_fire_patch(i)                = value_patch
       this%m_deadcrootc_storage_to_fire_patch(i)        = value_patch
       this%m_deadcrootc_xfer_to_fire_patch(i)           = value_patch
       this%m_gresp_storage_to_fire_patch(i)             = value_patch
       this%m_gresp_xfer_to_fire_patch(i)                = value_patch

       this%m_leafc_to_litter_fire_patch(i)              = value_patch
       this%m_leafc_storage_to_litter_fire_patch(i)      = value_patch
       this%m_leafc_xfer_to_litter_fire_patch(i)         = value_patch
       this%m_livestemc_to_litter_fire_patch(i)          = value_patch
       this%m_livestemc_storage_to_litter_fire_patch(i)  = value_patch
       this%m_livestemc_xfer_to_litter_fire_patch(i)     = value_patch
       this%m_livestemc_to_deadstemc_fire_patch(i)       = value_patch
       this%m_deadstemc_to_litter_fire_patch(i)          = value_patch
       this%m_deadstemc_storage_to_litter_fire_patch(i)  = value_patch
       this%m_deadstemc_xfer_to_litter_fire_patch(i)     = value_patch
       this%m_frootc_to_litter_fire_patch(i)             = value_patch
       this%m_frootc_storage_to_litter_fire_patch(i)     = value_patch
       this%m_frootc_xfer_to_litter_fire_patch(i)        = value_patch
       this%m_livecrootc_to_litter_fire_patch(i)         = value_patch
       this%m_livecrootc_storage_to_litter_fire_patch(i) = value_patch
       this%m_livecrootc_xfer_to_litter_fire_patch(i)    = value_patch
       this%m_livecrootc_to_deadcrootc_fire_patch(i)     = value_patch
       this%m_deadcrootc_to_litter_fire_patch(i)         = value_patch
       this%m_deadcrootc_storage_to_litter_fire_patch(i) = value_patch
       this%m_deadcrootc_xfer_to_litter_fire_patch(i)    = value_patch
       this%m_gresp_storage_to_litter_fire_patch(i)      = value_patch
       this%m_gresp_xfer_to_litter_fire_patch(i)         = value_patch

       this%leafc_xfer_to_leafc_patch(i)                 = value_patch
       this%frootc_xfer_to_frootc_patch(i)               = value_patch
       this%livestemc_xfer_to_livestemc_patch(i)         = value_patch
       this%deadstemc_xfer_to_deadstemc_patch(i)         = value_patch
       this%livecrootc_xfer_to_livecrootc_patch(i)       = value_patch
       this%deadcrootc_xfer_to_deadcrootc_patch(i)       = value_patch
       this%leafc_to_litter_patch(i)                     = value_patch
       this%frootc_to_litter_patch(i)                    = value_patch
       this%leaf_mr_patch(i)                             = value_patch
       this%froot_mr_patch(i)                            = value_patch
       this%livestem_mr_patch(i)                         = value_patch
       this%livecroot_mr_patch(i)                        = value_patch
       this%grain_mr_patch(i)                            = value_patch
       this%leaf_curmr_patch(i)                          = value_patch
       this%froot_curmr_patch(i)                         = value_patch
       this%livestem_curmr_patch(i)                      = value_patch
       this%livecroot_curmr_patch(i)                     = value_patch
       this%grain_curmr_patch(i)                         = value_patch
       this%leaf_xsmr_patch(i)                           = value_patch
       this%froot_xsmr_patch(i)                          = value_patch
       this%livestem_xsmr_patch(i)                       = value_patch
       this%livecroot_xsmr_patch(i)                      = value_patch
       this%grain_xsmr_patch(i)                          = value_patch
       this%psnsun_to_cpool_patch(i)                     = value_patch
       this%psnshade_to_cpool_patch(i)                   = value_patch
       this%cpool_to_xsmrpool_patch(i)                   = value_patch
       this%cpool_to_leafc_patch(i)                      = value_patch
       this%cpool_to_leafc_storage_patch(i)              = value_patch
       this%cpool_to_frootc_patch(i)                     = value_patch
       this%cpool_to_frootc_storage_patch(i)             = value_patch
       this%cpool_to_livestemc_patch(i)                  = value_patch
       this%cpool_to_livestemc_storage_patch(i)          = value_patch
       this%cpool_to_deadstemc_patch(i)                  = value_patch
       this%cpool_to_deadstemc_storage_patch(i)          = value_patch
       this%cpool_to_livecrootc_patch(i)                 = value_patch
       this%cpool_to_livecrootc_storage_patch(i)         = value_patch
       this%cpool_to_deadcrootc_patch(i)                 = value_patch
       this%cpool_to_deadcrootc_storage_patch(i)         = value_patch
       this%cpool_to_gresp_storage_patch(i)              = value_patch
       this%cpool_leaf_gr_patch(i)                       = value_patch
       this%cpool_leaf_storage_gr_patch(i)               = value_patch
       this%transfer_leaf_gr_patch(i)                    = value_patch
       this%cpool_froot_gr_patch(i)                      = value_patch
       this%cpool_froot_storage_gr_patch(i)              = value_patch
       this%transfer_froot_gr_patch(i)                   = value_patch
       this%cpool_livestem_gr_patch(i)                   = value_patch
       this%cpool_livestem_storage_gr_patch(i)           = value_patch
       this%transfer_livestem_gr_patch(i)                = value_patch
       this%cpool_deadstem_gr_patch(i)                   = value_patch
       this%cpool_deadstem_storage_gr_patch(i)           = value_patch
       this%transfer_deadstem_gr_patch(i)                = value_patch
       this%cpool_livecroot_gr_patch(i)                  = value_patch
       this%cpool_livecroot_storage_gr_patch(i)          = value_patch
       this%transfer_livecroot_gr_patch(i)               = value_patch
       this%cpool_deadcroot_gr_patch(i)                  = value_patch
       this%cpool_deadcroot_storage_gr_patch(i)          = value_patch
       this%transfer_deadcroot_gr_patch(i)               = value_patch
       this%leafc_storage_to_xfer_patch(i)               = value_patch
       this%frootc_storage_to_xfer_patch(i)              = value_patch
       this%livestemc_storage_to_xfer_patch(i)           = value_patch
       this%deadstemc_storage_to_xfer_patch(i)           = value_patch
       this%livecrootc_storage_to_xfer_patch(i)          = value_patch
       this%deadcrootc_storage_to_xfer_patch(i)          = value_patch
       this%gresp_storage_to_xfer_patch(i)               = value_patch
       this%livestemc_to_deadstemc_patch(i)              = value_patch
       this%livecrootc_to_deadcrootc_patch(i)            = value_patch
       if (.not. use_ed) then
          !TODO - fix this when use_cn and use_ed are truly separate
          this%gpp_patch(i)                              = value_patch
          this%gpp_before_downreg_patch(i)               = value_patch
       end if
       this%mr_patch(i)                                  = value_patch
       this%current_gr_patch(i)                          = value_patch
       this%transfer_gr_patch(i)                         = value_patch
       this%storage_gr_patch(i)                          = value_patch
       this%gr_patch(i)                                  = value_patch
       this%ar_patch(i)                                  = value_patch
       this%rr_patch(i)                                  = value_patch
       if (.not. use_ed) then
          !TODO - fix this when use_cn and use_ed are truly separate
          this%npp_patch(i)                              = value_patch 
       end if
       this%agnpp_patch(i)                               = value_patch
       this%bgnpp_patch(i)                               = value_patch
       this%litfall_patch(i)                             = value_patch
       this%vegfire_patch(i)                             = value_patch
       this%wood_harvestc_patch(i)                       = value_patch
       this%cinputs_patch(i)                             = value_patch
       this%coutputs_patch(i)                            = value_patch
       this%fire_closs_patch(i)                          = value_patch
       this%frootc_alloc_patch(i)                        = value_patch
       this%frootc_loss_patch(i)                         = value_patch
       this%leafc_alloc_patch(i)                         = value_patch
       this%leafc_loss_patch(i)                          = value_patch
       this%woodc_alloc_patch(i)                         = value_patch
       this%woodc_loss_patch(i)                          = value_patch
    end do

    if ( crop_prog )then
       do fi = 1,num_patch
          i = filter_patch(fi)
          this%xsmrpool_to_atm_patch(i)         = value_patch
          this%livestemc_to_litter_patch(i)     = value_patch
          this%grainc_to_food_patch(i)          = value_patch
          this%grainc_xfer_to_grainc_patch(i)   = value_patch
          this%cpool_to_grainc_patch(i)         = value_patch
          this%cpool_to_grainc_storage_patch(i) = value_patch
          this%cpool_grain_gr_patch(i)          = value_patch
          this%cpool_grain_storage_gr_patch(i)  = value_patch
          this%transfer_grain_gr_patch(i)       = value_patch
          this%grainc_storage_to_xfer_patch(i)  = value_patch
       end do
    end if

    do j = 1, nlevdecomp_full
       do fi = 1,num_column
          i = filter_column(fi)

          this%phenology_c_to_litr_met_c_col(i,j)     = value_column
          this%phenology_c_to_litr_cel_c_col(i,j)     = value_column
          this%phenology_c_to_litr_lig_c_col(i,j)     = value_column

          this%gap_mortality_c_to_litr_met_c_col(i,j) = value_column
          this%gap_mortality_c_to_litr_cel_c_col(i,j) = value_column
          this%gap_mortality_c_to_litr_lig_c_col(i,j) = value_column
          this%gap_mortality_c_to_cwdc_col(i,j)       = value_column

          this%fire_mortality_c_to_cwdc_col(i,j)      = value_column
          this%m_c_to_litr_met_fire_col(i,j)          = value_column
          this%m_c_to_litr_cel_fire_col(i,j)          = value_column  
          this%m_c_to_litr_lig_fire_col(i,j)          = value_column

          this%harvest_c_to_litr_met_c_col(i,j)       = value_column             
          this%harvest_c_to_litr_cel_c_col(i,j)       = value_column             
          this%harvest_c_to_litr_lig_c_col(i,j)       = value_column             
          this%harvest_c_to_cwdc_col(i,j)             = value_column          

          this%hr_vr_col(i,j)                         = value_column
       end do
    end do

    do k = 1, ndecomp_pools
       do j = 1, nlevdecomp_full
          do fi = 1,num_column
             i = filter_column(fi)
             this%m_decomp_cpools_to_fire_vr_col(i,j,k) = value_column
             this%decomp_cpools_transport_tendency_col(i,j,k) = value_column
          end do
       end do
    end do

    do l = 1, ndecomp_cascade_transitions
       do fi = 1,num_column
          i = filter_column(fi)
          this%decomp_cascade_hr_col(i,l) = value_column
          this%decomp_cascade_ctransfer_col(i,l) = value_column
       end do
    end do

    do l = 1, ndecomp_cascade_transitions
       do j = 1, nlevdecomp_full
          do fi = 1,num_column
             i = filter_column(fi)
             this%decomp_cascade_hr_vr_col(i,j,l) = value_column
             this%decomp_cascade_ctransfer_vr_col(i,j,l) = value_column
             this%decomp_k_col(i,j,l) = value_column
          end do
       end do
    end do

    do k = 1, ndecomp_pools
       do fi = 1,num_column
          i = filter_column(fi)
          this%decomp_cpools_leached_col(i,k) = value_column
          this%m_decomp_cpools_to_fire_col(i,k) = value_column
          this%bgc_cpool_ext_inputs_vr_col(i,:, k) = value_column
          this%bgc_cpool_ext_loss_vr_col(i,:, k) = value_column
       end do
    end do

    do fi = 1,num_column
       i = filter_column(fi)

       this%hrv_deadstemc_to_prod10c_col(i)  = value_column        
       this%hrv_deadstemc_to_prod100c_col(i) = value_column  
       this%hrv_cropc_to_prod1c_col(i)       = value_column
       this%somc_fire_col(i)                 = value_column  
       this%prod1c_loss_col(i)               = value_column
       this%prod10c_loss_col(i)              = value_column
       this%prod100c_loss_col(i)             = value_column
       this%product_closs_col(i)             = value_column
       this%somhr_col(i)                     = value_column
       this%lithr_col(i)                     = value_column
       this%hr_col(i)                        = value_column
       this%sr_col(i)                        = value_column
       this%er_col(i)                        = value_column
       this%litfire_col(i)                   = value_column
       this%somfire_col(i)                   = value_column
       this%totfire_col(i)                   = value_column
       this%nep_col(i)                       = value_column
       this%nbp_col(i)                       = value_column
       this%nee_col(i)                       = value_column
       this%fire_closs_col(i)                = value_column
       this%cwdc_hr_col(i)                   = value_column
       this%cwdc_loss_col(i)                 = value_column
       this%litterc_loss_col(i)              = value_column
       this%som_c_leached_col(i)             = value_column

       ! Zero p2c column fluxes
       this%rr_col(i)                    = value_column  
       this%ar_col(i)                    = value_column  
       this%gpp_col(i)                   = value_column 
       this%npp_col(i)                   = value_column 
       this%fire_closs_col(i)            = value_column 
       this%litfall_col(i)               = value_column       
       this%vegfire_col(i)               = value_column       
       this%wood_harvestc_col(i)         = value_column 
       this%hrv_xsmrpool_to_atm_col(i)   = value_column
    end do

    do k = 1, ndecomp_pools
       do j = 1, nlevdecomp_full
          do fi = 1,num_column
             i = filter_column(fi)
             this%decomp_cpools_sourcesink_col(i,j,k) = value_column  
          end do
       end do
    end do

    ! pflotran
    do k = 1, ndecomp_pools
       do j = 1, nlevdecomp_full
          do fi = 1,num_column
             i = filter_column(fi)
             ! only initializing in the first time-step
             if ( this%externalc_to_decomp_cpools_col(i,j,k) == spval ) then
                this%externalc_to_decomp_cpools_col(i,j,k) = value_column
             end if
          end do
       end do
    end do

    do fi = 1,num_column
       i = filter_column(fi)
       this%f_co2_soil_col(i) = value_column
       ! only initializing in the first time-step
       if ( this%externalc_to_decomp_delta_col(i) == spval ) then
          this%externalc_to_decomp_delta_col(i) = value_column
       end if
    end do

    do j = 1, nlevdecomp_full
       do fi = 1,num_column
          i = filter_column(fi)
          this%f_co2_soil_vr_col(i,j) = value_column
       end do
    end do

  end subroutine SetValues

  !-----------------------------------------------------------------------
  subroutine ZeroDwt( this, bounds )
    !
    ! !DESCRIPTION
    ! Initialize flux variables needed for dynamic land use.
    !
    ! !ARGUMENTS:
    class(carbonflux_type) :: this
    type(bounds_type), intent(in)  :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer  :: c, j          ! indices
    !-----------------------------------------------------------------------

    ! set column-level conversion and product pool fluxes
    ! to 0 at the beginning of every timestep

    do c = bounds%begc,bounds%endc
       this%dwt_seedc_to_leaf_col(c)        = 0._r8
       this%dwt_seedc_to_deadstem_col(c)    = 0._r8
       this%dwt_conv_cflux_col(c)           = 0._r8
       this%lf_conv_cflux_col(c)            = 0._r8
       this%dwt_prod10c_gain_col(c)         = 0._r8
       this%dwt_prod100c_gain_col(c)        = 0._r8
    end do

    do j = 1, nlevdecomp_full
       do c = bounds%begc,bounds%endc
          this%dwt_frootc_to_litr_met_c_col(c,j)    = 0._r8
          this%dwt_frootc_to_litr_cel_c_col(c,j)    = 0._r8
          this%dwt_frootc_to_litr_lig_c_col(c,j)    = 0._r8
          this%dwt_livecrootc_to_cwdc_col(c,j)      = 0._r8
          this%dwt_deadcrootc_to_cwdc_col(c,j)      = 0._r8
       end do
    end do

  end subroutine ZeroDwt

  !-----------------------------------------------------------------------
  subroutine Summary(this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       isotope)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, perform patch and column-level carbon summary calculations
    !
    ! !USES:
    use clm_varctl       , only : iulog, use_cndv
    use clm_time_manager , only : get_step_size
    use clm_varcon       , only : secspday
    use clm_varpar       , only : nlevdecomp, ndecomp_pools, ndecomp_cascade_transitions
    use subgridAveMod    , only : p2c
    use tracer_varcon    , only : is_active_betr_bgc
    use MathfuncMod      , only : dot_sum
    !
    ! !ARGUMENTS:
    class(carbonflux_type)                 :: this
    type(bounds_type)      , intent(in)    :: bounds          
    integer                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    character(len=*)       , intent(in)    :: isotope   
    !
    ! !LOCAL VARIABLES:
    real(r8) :: nfixlags, dtime ! temp variables for making lagged npp
    integer  :: c,p,j,k,l       ! indices
    integer  :: fp,fc           ! lake filter indices
    real(r8) :: maxdepth        ! depth to integrate soil variables
    !-----------------------------------------------------------------------

   associate(& 
        is_litter =>    decomp_cascade_con%is_litter , & ! Input:  [logical (:) ]  TRUE => pool is a litter pool
        is_soil   =>    decomp_cascade_con%is_soil   , & ! Input:  [logical (:) ]  TRUE => pool is a soil pool  
        is_cwd    =>    decomp_cascade_con%is_cwd      & ! Input:  [logical (:) ]  TRUE => pool is a cwd pool   
        )

    ! patch loop
    do fp = 1,num_soilp
       p = filter_soilp(fp)

       ! calculate pft-level summary carbon fluxes and states

       ! gross primary production (GPP)
       this%gpp_patch(p) = &
            this%psnsun_to_cpool_patch(p) + &
            this%psnshade_to_cpool_patch(p)

       ! maintenance respiration (MR)
       if ( trim(isotope) == 'c13' .or. trim(isotope) == 'c14') then
          this%leaf_mr_patch(p)      = this%leaf_curmr_patch(p)      + this%leaf_xsmr_patch(p)
          this%froot_mr_patch(p)     = this%froot_curmr_patch(p)     + this%froot_xsmr_patch(p)
          this%livestem_mr_patch(p)  = this%livestem_curmr_patch(p)  + this%livestem_xsmr_patch(p)
          this%livecroot_mr_patch(p) = this%livecroot_curmr_patch(p) + this%livecroot_xsmr_patch(p)
       endif

       this%mr_patch(p)  = &
            this%leaf_mr_patch(p)     + &
            this%froot_mr_patch(p)    + &
            this%livestem_mr_patch(p) + &
            this%livecroot_mr_patch(p)

       ! growth respiration (GR)
       ! current GR is respired this time step for new growth displayed in this timestep
       this%current_gr_patch(p) = &
            this%cpool_leaf_gr_patch(p)      + &
            this%cpool_froot_gr_patch(p)     + &
            this%cpool_livestem_gr_patch(p)  + &
            this%cpool_deadstem_gr_patch(p)  + &
            this%cpool_livecroot_gr_patch(p) + &
            this%cpool_deadcroot_gr_patch(p)

       ! transfer GR is respired this time step for transfer growth displayed in this timestep
       this%transfer_gr_patch(p) = &
            this%transfer_leaf_gr_patch(p)      + &
            this%transfer_froot_gr_patch(p)     + &
            this%transfer_livestem_gr_patch(p)  + &
            this%transfer_deadstem_gr_patch(p)  + &
            this%transfer_livecroot_gr_patch(p) + &
            this%transfer_deadcroot_gr_patch(p)

       ! storage GR is respired this time step for growth sent to storage for later display
       this%storage_gr_patch(p) = &
            this%cpool_leaf_storage_gr_patch(p)      + &
            this%cpool_froot_storage_gr_patch(p)     + &
            this%cpool_livestem_storage_gr_patch(p)  + &
            this%cpool_deadstem_storage_gr_patch(p)  + &
            this%cpool_livecroot_storage_gr_patch(p) + &
            this%cpool_deadcroot_storage_gr_patch(p)

       if ( crop_prog .and. pft%itype(p) >= npcropmin )then
          this%mr_patch(p) = &
               this%mr_patch(p) + &
               this%grain_mr_patch(p)

          this%current_gr_patch(p) = &
               this%current_gr_patch(p) + &
               this%cpool_grain_gr_patch(p)

          this%transfer_gr_patch(p) = &
               this%transfer_gr_patch(p) + &
               this%transfer_grain_gr_patch(p)

          this%storage_gr_patch(p) = &
               this%storage_gr_patch(p) + &
               this%cpool_grain_storage_gr_patch(p)
       end if

       ! GR is the sum of current + transfer + storage GR
       this%gr_patch(p) = &
            this%current_gr_patch(p)  + &
            this%transfer_gr_patch(p) + &
            this%storage_gr_patch(p)

       ! autotrophic respiration (AR)
       if ( crop_prog .and. pft%itype(p) >= npcropmin )then
          this%ar_patch(p) = &
               this%mr_patch(p) + &
               this%gr_patch(p) + &
               this%xsmrpool_to_atm_patch(p) ! xsmr... is -ve (slevis)
       else
          this%ar_patch(p) = &
               this%mr_patch(p) + &
               this%gr_patch(p)
       end if



       ! net primary production (NPP)
       this%npp_patch(p) = &
            this%gpp_patch(p) - &
            this%ar_patch(p)

       ! update the annual NPP accumulator, for use in allocation code 
       if (trim(isotope) == 'bulk') then      
          this%tempsum_npp_patch(p) = &
               this%tempsum_npp_patch(p) + &
               this%npp_patch(p)
       end if

       ! litterfall (LITFALL)

       this%litfall_patch(p) = &
            this%leafc_to_litter_patch(p)                     + &
            this%frootc_to_litter_patch(p)                    + &
            this%m_leafc_to_litter_patch(p)                   + &
            this%m_leafc_storage_to_litter_patch(p)           + &
            this%m_leafc_xfer_to_litter_patch(p)              + &
            this%m_frootc_to_litter_patch(p)                  + &
            this%m_frootc_storage_to_litter_patch(p)          + &
            this%m_frootc_xfer_to_litter_patch(p)             + &
            this%m_livestemc_to_litter_patch(p)               + &
            this%m_livestemc_storage_to_litter_patch(p)       + &
            this%m_livestemc_xfer_to_litter_patch(p)          + &
            this%m_deadstemc_to_litter_patch(p)               + &
            this%m_deadstemc_storage_to_litter_patch(p)       + &
            this%m_deadstemc_xfer_to_litter_patch(p)          + &
            this%m_livecrootc_to_litter_patch(p)              + &
            this%m_livecrootc_storage_to_litter_patch(p)      + &
            this%m_livecrootc_xfer_to_litter_patch(p)         + &
            this%m_deadcrootc_to_litter_patch(p)              + &
            this%m_deadcrootc_storage_to_litter_patch(p)      + &
            this%m_deadcrootc_xfer_to_litter_patch(p)         + &
            this%m_gresp_storage_to_litter_patch(p)           + &
            this%m_gresp_xfer_to_litter_patch(p)              + &
            
            this%m_leafc_to_litter_fire_patch(p)              + &
            this%m_leafc_storage_to_litter_fire_patch(p)      + &
            this%m_leafc_xfer_to_litter_fire_patch(p)         + &
            this%m_livestemc_to_litter_fire_patch(p)          + &
            this%m_livestemc_storage_to_litter_fire_patch(p)  + &
            this%m_livestemc_xfer_to_litter_fire_patch(p)     + &
            this%m_deadstemc_to_litter_fire_patch(p)          + &
            this%m_deadstemc_storage_to_litter_fire_patch(p)  + &
            this%m_deadstemc_xfer_to_litter_fire_patch(p)     + &
            this%m_frootc_to_litter_fire_patch(p)             + &
            this%m_frootc_storage_to_litter_fire_patch(p)     + &
            this%m_frootc_xfer_to_litter_fire_patch(p)        + &
            this%m_livecrootc_to_litter_fire_patch(p)         + &
            this%m_livecrootc_storage_to_litter_fire_patch(p) + &
            this%m_livecrootc_xfer_to_litter_fire_patch(p)    + &
            this%m_deadcrootc_to_litter_fire_patch(p)         + &
            this%m_deadcrootc_storage_to_litter_fire_patch(p) + &
            this%m_deadcrootc_xfer_to_litter_fire_patch(p)    + &
            this%m_gresp_storage_to_litter_fire_patch(p)      + &
            this%m_gresp_xfer_to_litter_fire_patch(p)         + &
            
            this%hrv_leafc_to_litter_patch(p)                 + &
            this%hrv_leafc_storage_to_litter_patch(p)         + &
            this%hrv_leafc_xfer_to_litter_patch(p)            + &
            this%hrv_frootc_to_litter_patch(p)                + &
            this%hrv_frootc_storage_to_litter_patch(p)        + &
            this%hrv_frootc_xfer_to_litter_patch(p)           + &
            this%hrv_livestemc_to_litter_patch(p)             + &
            this%hrv_livestemc_storage_to_litter_patch(p)     + &
            this%hrv_livestemc_xfer_to_litter_patch(p)        + &
            this%hrv_deadstemc_storage_to_litter_patch(p)     + &
            this%hrv_deadstemc_xfer_to_litter_patch(p)        + &
            this%hrv_livecrootc_to_litter_patch(p)            + &
            this%hrv_livecrootc_storage_to_litter_patch(p)    + &
            this%hrv_livecrootc_xfer_to_litter_patch(p)       + &
            this%hrv_deadcrootc_to_litter_patch(p)            + &
            this%hrv_deadcrootc_storage_to_litter_patch(p)    + &
            this%hrv_deadcrootc_xfer_to_litter_patch(p)       + &
            this%hrv_gresp_storage_to_litter_patch(p)         + &
            this%hrv_gresp_xfer_to_litter_patch(p)

       ! update the annual litfall accumulator, for use in mortality code
       if (use_cndv) then
          this%tempsum_litfall_patch(p) = &
               this%tempsum_litfall_patch(p) + &
               this%leafc_to_litter_patch(p) + &
               this%frootc_to_litter_patch(p)
       end if

       ! patch-level fire losses (VEGFIRE)
       this%vegfire_patch(p) = 0._r8

       ! patch-level wood harvest
       this%wood_harvestc_patch(p) = &
            this%hrv_deadstemc_to_prod10c_patch(p) + &
            this%hrv_deadstemc_to_prod100c_patch(p)
       if ( crop_prog .and. pft%itype(p) >= npcropmin )then
            this%wood_harvestc_patch(p) = &
            this%wood_harvestc_patch(p) + &
            this%hrv_cropc_to_prod1c_patch(p)
      end if

       ! patch-level carbon losses to fire changed by F. Li and S. Levis
       this%fire_closs_patch(p) = &
            this%m_leafc_to_fire_patch(p)                + &
            this%m_leafc_storage_to_fire_patch(p)        + &
            this%m_leafc_xfer_to_fire_patch(p)           + &
            this%m_frootc_to_fire_patch(p)               + &
            this%m_frootc_storage_to_fire_patch(p)       + &
            this%m_frootc_xfer_to_fire_patch(p)          + &
            this%m_livestemc_to_fire_patch(p)            + &
            this%m_livestemc_storage_to_fire_patch(p)    + &
            this%m_livestemc_xfer_to_fire_patch(p)       + &
            this%m_deadstemc_to_fire_patch(p)            + &
            this%m_deadstemc_storage_to_fire_patch(p)    + &
            this%m_deadstemc_xfer_to_fire_patch(p)       + &
            this%m_livecrootc_to_fire_patch(p)           + &
            this%m_livecrootc_storage_to_fire_patch(p)   + &
            this%m_livecrootc_xfer_to_fire_patch(p)      + &
            this%m_deadcrootc_to_fire_patch(p)           + &
            this%m_deadcrootc_storage_to_fire_patch(p)   + &
            this%m_deadcrootc_xfer_to_fire_patch(p)      + &
            this%m_gresp_storage_to_fire_patch(p)        + &
            this%m_gresp_xfer_to_fire_patch(p)

       if ( crop_prog .and. pft%itype(p) >= npcropmin )then

          this%litfall_patch(p) =                  &
               this%litfall_patch(p)             + &
               this%livestemc_to_litter_patch(p) + &
               this%grainc_to_food_patch(p)
       end if

       ! new summary variables for CLAMP

       ! (FROOTC_ALLOC) - fine root C allocation
       this%frootc_alloc_patch(p) = &
            this%frootc_xfer_to_frootc_patch(p)    + &
            this%cpool_to_frootc_patch(p)     

       ! (FROOTC_LOSS) - fine root C loss changed by F. Li and S. Levis
       this%frootc_loss_patch(p) = &
            this%m_frootc_to_litter_patch(p)       + &
            this%m_frootc_to_fire_patch(p)         + &
            this%m_frootc_to_litter_fire_patch(p)  + &
            this%hrv_frootc_to_litter_patch(p)     + &
            this%frootc_to_litter_patch(p)
       
       ! (LEAFC_ALLOC) - leaf C allocation
       this%leafc_alloc_patch(p) = &
            this%leafc_xfer_to_leafc_patch(p)    + &
            this%cpool_to_leafc_patch(p)     

       ! (LEAFC_LOSS) - leaf C loss changed by F. Li and S. Levis
       this%leafc_loss_patch(p) = &
            this%m_leafc_to_litter_patch(p)      + &
            this%m_leafc_to_fire_patch(p)        + &
            this%m_leafc_to_litter_fire_patch(p) + &
            this%hrv_leafc_to_litter_patch(p)    + &
            this%leafc_to_litter_patch(p)

       if ( crop_prog .and. pft%itype(p) >= npcropmin )then
            this%leafc_loss_patch(p) = &
            this%leafc_loss_patch(p) + &
            this%hrv_leafc_to_prod1c_patch(p)
       end if


       ! (WOODC_ALLOC) - wood C allocation
       this%woodc_alloc_patch(p) = &
            this%livestemc_xfer_to_livestemc_patch(p)   + &
            this%deadstemc_xfer_to_deadstemc_patch(p)   + &
            this%livecrootc_xfer_to_livecrootc_patch(p) + &
            this%deadcrootc_xfer_to_deadcrootc_patch(p) + &
            this%cpool_to_livestemc_patch(p)            + &
            this%cpool_to_deadstemc_patch(p)            + &
            this%cpool_to_livecrootc_patch(p)           + &
            this%cpool_to_deadcrootc_patch(p)

       ! (WOODC_LOSS) - wood C loss
       this%woodc_loss_patch(p) = &
            this%m_livestemc_to_litter_patch(p)            + &
            this%m_deadstemc_to_litter_patch(p)            + &
            this%m_livecrootc_to_litter_patch(p)           + &
            this%m_deadcrootc_to_litter_patch(p)           + &
            this%m_livestemc_to_fire_patch(p)              + &
            this%m_deadstemc_to_fire_patch(p)              + &
            this%m_livecrootc_to_fire_patch(p)             + &
            this%m_deadcrootc_to_fire_patch(p)             + &
            this%hrv_livestemc_to_litter_patch(p)          + &
            this%hrv_livestemc_storage_to_litter_patch(p)  + &
            this%hrv_livestemc_xfer_to_litter_patch(p)     + &
            this%hrv_deadstemc_to_prod10c_patch(p)         + &
            this%hrv_deadstemc_to_prod100c_patch(p)        + &
            this%hrv_deadstemc_storage_to_litter_patch(p)  + &
            this%hrv_deadstemc_xfer_to_litter_patch(p)     + &
            this%hrv_livecrootc_to_litter_patch(p)         + &
            this%hrv_livecrootc_storage_to_litter_patch(p) + &
            this%hrv_livecrootc_xfer_to_litter_patch(p)    + &
            this%hrv_deadcrootc_to_litter_patch(p)         + &
            this%hrv_deadcrootc_storage_to_litter_patch(p) + &
            this%hrv_deadcrootc_xfer_to_litter_patch(p)   
        ! putting the harvested crop stem and grain in the wood loss bdrewniak
        if ( crop_prog .and. pft%itype(p) >= npcropmin )then
             this%woodc_loss_patch(p) = &
             this%woodc_loss_patch(p) + &
             this%hrv_grainc_to_prod1c_patch(p) + &
             this%hrv_livestemc_to_prod1c_patch(p)
        end if

    end do  ! end of patches loop

    ! use p2c routine to get selected column-average patch-level fluxes and states

    call p2c(bounds, num_soilc, filter_soilc, &
         this%gpp_patch(bounds%begp:bounds%endp), &
         this%gpp_col(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%ar_patch(bounds%begp:bounds%endp), &
         this%ar_col(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%npp_patch(bounds%begp:bounds%endp), &
         this%npp_col(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%vegfire_patch(bounds%begp:bounds%endp), &
         this%vegfire_col(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%wood_harvestc_patch(bounds%begp:bounds%endp), &
         this%wood_harvestc_col(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%fire_closs_patch(bounds%begp:bounds%endp), &
         this%fire_closs_p2c_col(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%litfall_patch(bounds%begp:bounds%endp), &
         this%litfall_col(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%hrv_xsmrpool_to_atm_patch(bounds%begp:bounds%endp), &
         this%hrv_xsmrpool_to_atm_col(bounds%begc:bounds%endc))

    if ( trim(isotope) == 'bulk') then
       if (nfix_timeconst > 0._r8 .and. nfix_timeconst < 500._r8 ) then

          ! this code is to calculate an exponentially-relaxed npp value for use in NDynamics code
          dtime = get_step_size()
          nfixlags = nfix_timeconst * secspday

          do fc = 1,num_soilc
             c = filter_soilc(fc)
             if ( this%lag_npp_col(c) /= spval ) then
                this%lag_npp_col(c) = &
                     this%lag_npp_col(c) * exp(-dtime/nfixlags) + &
                     this%npp_col(c) * (1._r8 - exp(-dtime/nfixlags))
             else
                ! first timestep
                this%lag_npp_col(c) = this%npp_col(c)
             endif
          end do
       endif
    endif

    ! column variables

    ! some zeroing
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%cwdc_loss_col(c)          = 0._r8
       this%som_c_leached_col(c)      = 0._r8
    end do
    
    if ( (.not. is_active_betr_bgc           ) .and. &
         (.not. (use_pflotran .and. pf_cmode))) then

      ! vertically integrate HR and decomposition cascade fluxes
      do k = 1, ndecomp_cascade_transitions

       do j = 1,nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)

             this%decomp_cascade_ctransfer_col(c,k) = &
                  this%decomp_cascade_ctransfer_col(c,k) + &
                  this%decomp_cascade_ctransfer_vr_col(c,j,k) * dzsoi_decomp(j) 
          end do
       end do
      end do


      ! total heterotrophic respiration (HR)
      do fc = 1,num_soilc
        c = filter_soilc(fc)
        this%hr_col(c) = &
            this%lithr_col(c) + &
            this%somhr_col(c)
      end do


    elseif (is_active_betr_bgc) then

       do fc = 1, num_soilc
        c = filter_soilc(fc)
        this%hr_col(c) = dot_sum(this%hr_vr_col(c,1:nlevdecomp),dzsoi_decomp(1:nlevdecomp)) 
      enddo
    endif
    

    ! bgc interface & pflotran:
    !----------------------------------------------------------------
    if (use_bgc_interface) then
        call CSummary_interface(this, bounds, num_soilc, filter_soilc)
    end if
    !! CSummary_interface: hr_col(c) will be used below
    !----------------------------------------------------------------
     
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       ! total soil respiration, heterotrophic + root respiration (SR)
       this%sr_col(c) = &
            this%rr_col(c) + &
            this%hr_col(c)

       ! total ecosystem respiration, autotrophic + heterotrophic (ER)
       this%er_col(c) = &
            this%ar_col(c) + &
            this%hr_col(c)

       ! litter fire losses (LITFIRE)
       this%litfire_col(c) = 0._r8

       ! total product loss
       this%product_closs_col(c) = &
            this%prod10c_loss_col(c)  + &
            this%prod100c_loss_col(c) + & 
            this%prod1c_loss_col(c)

       ! soil organic matter fire losses (SOMFIRE)
       this%somfire_col(c) = 0._r8

       ! total ecosystem fire losses (TOTFIRE)
       this%totfire_col(c) = &
            this%litfire_col(c) + &
            this%somfire_col(c) + &
            this%vegfire_col(c)
    end do

    ! vertically integrate column-level carbon fire losses
    do l = 1, ndecomp_pools
       do j = 1,nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%m_decomp_cpools_to_fire_col(c,l) = &
                  this%m_decomp_cpools_to_fire_col(c,l) + &
                  this%m_decomp_cpools_to_fire_vr_col(c,j,l)*dzsoi_decomp(j)
          end do
       end do
    end do

    ! column-level carbon losses to fire, including pft losses
    do fc = 1,num_soilc
       c = filter_soilc(fc)

       this%fire_closs_col(c) = this%fire_closs_p2c_col(c) 
       do l = 1, ndecomp_pools
          this%fire_closs_col(c) = &
               this%fire_closs_col(c) + &
               this%m_decomp_cpools_to_fire_col(c,l)
       end do

       ! column-level carbon losses due to landcover change
       this%dwt_closs_col(c) = &
            this%dwt_conv_cflux_col(c)

       ! net ecosystem production, excludes fire flux, landcover change, and loss from wood products, positive for sink (NEP)
       this%nep_col(c) = &
            this%gpp_col(c) - &
            this%er_col(c)

       ! net biome production of carbon, includes depletion from: fire flux, landcover change flux, and loss
       ! from wood products pools, positive for sink (NBP)
       this%nbp_col(c) =             &
            this%nep_col(c)        - &
            this%fire_closs_col(c) - &
            this%dwt_closs_col(c)  - &
            this%product_closs_col(c)

       ! net ecosystem exchange of carbon, includes fire flux, landcover change flux, loss
       ! from wood products pools, and hrv_xsmrpool flux, positive for source (NEE)
       this%nee_col(c) =                &
           -this%nep_col(c)           + &
            this%fire_closs_col(c)    + &
            this%dwt_closs_col(c)     + &
            this%product_closs_col(c) + &
            this%hrv_xsmrpool_to_atm_col(c)

       ! land use flux and land uptake
       this%landuseflux_col(c) = &
            this%dwt_closs_col(c) + &
            this%product_closs_col(c)

       this%landuptake_col(c) = &
            this%nee_col(c) - &
            this%landuseflux_col(c)
    end do

    ! for vertically-resolved soil biogeochemistry, calculate some diagnostics of carbon pools to a given depth

    if ( (.not. is_active_betr_bgc)           .and. &
         (.not.(use_pflotran .and. pf_cmode)) ) then

    ! _col(cWDC_HR) - coarse woody debris heterotrophic respiration
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%cwdc_hr_col(c) = 0._r8
    end do

    ! _col(cWDC_LOSS) - coarse woody debris C loss
    do l = 1, ndecomp_pools
       if ( is_cwd(l) ) then
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%cwdc_loss_col(c) = &
                  this%cwdc_loss_col(c) + &
                  this%m_decomp_cpools_to_fire_col(c,l)
          end do
       end if
    end do

    do k = 1, ndecomp_cascade_transitions
       if ( is_cwd(decomp_cascade_con%cascade_donor_pool(k)) ) then
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%cwdc_loss_col(c) = &
                  this%cwdc_loss_col(c) + &
                  this%decomp_cascade_ctransfer_col(c,k)
          end do
       end if
    end do

    ! (LITTERC_LOSS) - litter C loss      
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%litterc_loss_col(c) = this%lithr_col(c)  
    end do
    do l = 1, ndecomp_pools
       if ( is_litter(l) ) then
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%litterc_loss_col(c) = &
                  this%litterc_loss_col(c) + &
                  this%m_decomp_cpools_to_fire_col(c,l)
          end do
       end if
    end do
    
 
      do k = 1, ndecomp_cascade_transitions
       if ( is_litter(decomp_cascade_con%cascade_donor_pool(k)) ) then
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%litterc_loss_col(c) = &
                  this%litterc_loss_col(c) + &
                  this%decomp_cascade_ctransfer_col(c,k)
          end do
       end if
      end do

   else if ((use_pflotran .and. pf_cmode)) then

      ! add up all vertical transport tendency terms and calculate total som leaching loss as the sum of these
      do l = 1, ndecomp_pools
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%decomp_cpools_leached_col(c,l) = 0._r8
       end do
       do j = 1, nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%decomp_cpools_leached_col(c,l) = &
                  this%decomp_cpools_leached_col(c,l) + &
                  this%decomp_cpools_transport_tendency_col(c,j,l) * dzsoi_decomp(j)
          end do
       end do
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%som_c_leached_col(c) = &
               this%som_c_leached_col(c) + &
               this%decomp_cpools_leached_col(c,l)
       end do
      end do
    endif
    
    ! debug
    do fc = 1,num_soilc
        c = filter_soilc(fc)
        this%plant_to_litter_cflux(c) = 0._r8
        this%plant_to_cwd_cflux(c) = 0._r8
        do j = 1, nlevdecomp
            this%plant_to_litter_cflux(c) = &
                this%plant_to_litter_cflux(c)  + &
                this%phenology_c_to_litr_met_c_col(c,j)* dzsoi_decomp(j) + &
                this%phenology_c_to_litr_cel_c_col(c,j)* dzsoi_decomp(j) + &
                this%phenology_c_to_litr_lig_c_col(c,j)* dzsoi_decomp(j) + &
                this%gap_mortality_c_to_litr_met_c_col(c,j)* dzsoi_decomp(j) + &
                this%gap_mortality_c_to_litr_cel_c_col(c,j)* dzsoi_decomp(j) + &
                this%gap_mortality_c_to_litr_lig_c_col(c,j)* dzsoi_decomp(j) + &
                this%m_c_to_litr_met_fire_col(c,j)* dzsoi_decomp(j) + &
                this%m_c_to_litr_cel_fire_col(c,j)* dzsoi_decomp(j) + &
                this%m_c_to_litr_lig_fire_col(c,j)* dzsoi_decomp(j)
            this%plant_to_cwd_cflux(c) = &
                this%plant_to_cwd_cflux(c) + &
                this%gap_mortality_c_to_cwdc_col(c,j)* dzsoi_decomp(j) + &
                this%fire_mortality_c_to_cwdc_col(c,j)* dzsoi_decomp(j)
        end do
    end do
    
  end associate
  end subroutine Summary

!!-------------------------------------------------------------------------------------------------
! !INTERFACE:
subroutine CSummary_interface(this, bounds, num_soilc, filter_soilc)
!
! !DESCRIPTION:
!! bgc interface & pflotran:
! On the radiation time step, perform column-level carbon
! summary calculations, which mainly from PFLOTRAN bgc
!
! !USES:
   use shr_sys_mod, only: shr_sys_flush
   use clm_varpar , only: nlevdecomp,ndecomp_pools,ndecomp_cascade_transitions
   use clm_varpar , only: i_met_lit, i_cel_lit, i_lig_lit, i_cwd
   use clm_time_manager    , only : get_step_size
!
! !ARGUMENTS:
   implicit none
   class(carbonflux_type)          :: this
   type(bounds_type) ,  intent(in) :: bounds
   integer,             intent(in) :: num_soilc       ! number of soil columns in filter
   integer,             intent(in) :: filter_soilc(:) ! filter for soil columns
!
! !CALLED FROM:
! subroutine Summary (if plotran bgc coupled with CLM-CN
!
! !REVISION HISTORY:
!!06/17/2015: modified by Gangsheng Wang
! !
! !LOCAL VARIABLES:
   real(r8) :: dtime                ! time-step (s)
   integer :: c,j,l                 ! indices
   integer :: fc                    ! column filter indices

    associate(&
        is_litter =>    decomp_cascade_con%is_litter , & ! Input:  [logical (:) ]  TRUE => pool is a litter pool
        is_soil   =>    decomp_cascade_con%is_soil   , & ! Input:  [logical (:) ]  TRUE => pool is a soil pool
        is_cwd    =>    decomp_cascade_con%is_cwd      & ! Input:  [logical (:) ]  TRUE => pool is a cwd pool
        )

    dtime = get_step_size()
!!---------------------------------------------------------------------------------------------------
    if (use_pflotran.and.pf_cmode) then
     ! total heterotrophic respiration (HR)
       this%hr_col(:) = 0._r8
       do j = 1,nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%hr_col(c) = this%hr_col(c) + &
                this%hr_vr_col(c,j) * dzsoi_decomp(j)
          end do
       end do

       ! new variable to account for co2 exchange (not all HR goes to atm at current time-step)
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%f_co2_soil_col(c) = 0._r8
       end do
       do j = 1,nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%f_co2_soil_col(c) = this%f_co2_soil_col(c) + &
                this%f_co2_soil_vr_col(c,j) * dzsoi_decomp(j)
          end do
       end do


    ! ---------------------------------------------------------
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%cwdc_hr_col(c)      = 0._r8
          this%cwdc_loss_col(c)    = 0._r8
          this%litterc_loss_col(c) = 0._r8
       end do

       do l = 1, ndecomp_pools
          if ( is_cwd(l) ) then
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                do j = 1, nlevdecomp
                   this%cwdc_loss_col(c) = &
                      this%cwdc_loss_col(c) + &
                      this%decomp_cpools_sourcesink_col(c,j,l) / dtime
                end do
             end do
          end if

          if ( is_litter(l) ) then
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                do j = 1, nlevdecomp
                   this%litterc_loss_col(c) = &
                      this%litterc_loss_col(c) + &
                      this%decomp_cpools_sourcesink_col(c,j,l) / dtime
                end do
             end do
          end if

       end do
    end if !!if (use_pflotran.and.pf_cmode)

   ! add up all vertically-resolved addition/removal rates (gC/m3/s) of decomp_pools for PFLOTRAN-bgc
    ! (note: this can be for general purpose, although here added an 'if...endif' block for PF-bgc)
    ! first, need to save the total plant C adding/removing to decomposing pools at previous time-step
    ! for calculating the net changes, which are used to do balance check
    this%externalc_to_decomp_delta_col(:) = 0._r8
    do l = 1, ndecomp_pools
       do j = 1, nlevdecomp
          do fc = 1, num_soilc
             c = filter_soilc(fc)
             this%externalc_to_decomp_delta_col(c) = this%externalc_to_decomp_delta_col(c) + &
                                this%externalc_to_decomp_cpools_col(c,j,l)*dzsoi_decomp(j)
          end do
       end do
    end do
!write(*,'(A40,E14.6)')">>>DEBUG | externC[t-1]=",this%externalc_to_decomp_delta_col(1)*dtime
    !
    ! do the initialization for the following variable here.
    ! DON'T do so in the beginning of CLM-CN time-step (otherwise the above saved will not work)
    this%externalc_to_decomp_cpools_col(:,:,:) = 0._r8

    do l = 1, ndecomp_pools
       do j = 1, nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)

             ! for litter C pools
             if (l==i_met_lit) then
                this%externalc_to_decomp_cpools_col(c,j,l) =                 &
                    this%externalc_to_decomp_cpools_col(c,j,l)               &
                        + this%phenology_c_to_litr_met_c_col(c,j)            &
                        + this%dwt_frootc_to_litr_met_c_col(c,j)             &
                        + this%gap_mortality_c_to_litr_met_c_col(c,j)        &
                        + this%harvest_c_to_litr_met_c_col(c,j)              !!&
!                        + this%m_c_to_litr_met_fire_col(c,j)                 &
!                        + this%decomp_cpools_transport_tendency_col(c,j,l)   &
!                        - this%m_decomp_cpools_to_fire_vr_col(c,j,l)

             elseif (l==i_cel_lit) then
                this%externalc_to_decomp_cpools_col(c,j,l) =                 &
                    this%externalc_to_decomp_cpools_col(c,j,l)               &
                        + this%phenology_c_to_litr_cel_c_col(c,j)            &
                        + this%dwt_frootc_to_litr_cel_c_col(c,j)             &
                        + this%gap_mortality_c_to_litr_cel_c_col(c,j)        &
                        + this%harvest_c_to_litr_cel_c_col(c,j)              !!&
!                        + this%m_c_to_litr_cel_fire_col(c,j)                 &
!                        + this%decomp_cpools_transport_tendency_col(c,j,l)   &
!                        - this%m_decomp_cpools_to_fire_vr_col(c,j,l)

             elseif (l==i_lig_lit) then
                this%externalc_to_decomp_cpools_col(c,j,l) =                 &
                    this%externalc_to_decomp_cpools_col(c,j,l)               &
                        + this%phenology_c_to_litr_lig_c_col(c,j)            &
                        + this%dwt_frootc_to_litr_lig_c_col(c,j)             &
                        + this%gap_mortality_c_to_litr_lig_c_col(c,j)        &
                        + this%harvest_c_to_litr_lig_c_col(c,j)              !!&
!                        + this%m_c_to_litr_lig_fire_col(c,j)                 &
!                        + this%decomp_cpools_transport_tendency_col(c,j,l)   &
!                        - this%m_decomp_cpools_to_fire_vr_col(c,j,l)

             ! for cwd
             elseif (l==i_cwd) then
                this%externalc_to_decomp_cpools_col(c,j,l) =                 &
                    this%externalc_to_decomp_cpools_col(c,j,l)               &
                        + this%dwt_livecrootc_to_cwdc_col(c,j)               &
                        + this%dwt_deadcrootc_to_cwdc_col(c,j)               &
                        + this%gap_mortality_c_to_cwdc_col(c,j)              &
                        + this%harvest_c_to_cwdc_col(c,j)                    !!&
!                        + this%fire_mortality_c_to_cwdc_col(c,j)             &
!                        + this%decomp_cpools_transport_tendency_col(c,j,l)   &
!                        - this%m_decomp_cpools_to_fire_vr_col(c,j,l)

             ! for som
             ! no external input to som
!             else
!                this%externalc_to_decomp_cpools_col(c,j,l) =                 &
!                    this%externalc_to_decomp_cpools_col(c,j,l)               &
!                        + this%decomp_cpools_transport_tendency_col(c,j,l)   &
!                        - this%m_decomp_cpools_to_fire_vr_col(c,j,l)

             end if

             ! the following is the net changes of plant C to decompible C poools between time-step
             ! in pflotran, decomposible C pools increments ARE from previous time-step (saved above);
             ! while, in CLM-CN all plant C pools are updated with current C fluxes among plant and ground/soil.
             ! therefore, when do balance check it is needed to adjust the time-lag of changes.
             this%externalc_to_decomp_delta_col(c) = this%externalc_to_decomp_delta_col(c) - &
                                this%externalc_to_decomp_cpools_col(c,j,l)*dzsoi_decomp(j)

             if (abs(this%externalc_to_decomp_cpools_col(c,j,l))<=1.e-20_r8) then
                 this%externalc_to_decomp_cpools_col(c,j,l) = 0._r8
             end if

          end do
       end do
    end do

    ! change the sign so that it is the increments from the previous time-step (unit: from g/m2/s)
    do fc = 1, num_soilc
       c = filter_soilc(fc)
       this%externalc_to_decomp_delta_col(c) = -this%externalc_to_decomp_delta_col(c)
    end do

    end associate
end subroutine CSummary_interface
!!-------------------------------------------------------------------------------------------------

  !------------------------------------------------------------  
  subroutine summary_rr(this, bounds, num_soilp, filter_soilp, num_soilc, filter_soilc)
  !
  ! description
  ! summarize root respiration

  use subgridAveMod    , only: p2c
  class(carbonflux_type) :: this  

  type(bounds_type), intent(in) :: bounds  
  integer, intent(in) :: num_soilp
  integer, intent(in) :: filter_soilp(:)
  integer, intent(in) :: num_soilc
  integer, intent(in) :: filter_soilc(:)
  integer :: fp, p
   ! patch loop
  do fp = 1,num_soilp
    p = filter_soilp(fp)  
    ! root respiration (RR)
    this%rr_patch(p) = &
    this%froot_mr_patch(p) + &
    this%cpool_froot_gr_patch(p) + &
    this%cpool_livecroot_gr_patch(p) + &
    this%cpool_deadcroot_gr_patch(p) + &
    this%transfer_froot_gr_patch(p) + &
    this%transfer_livecroot_gr_patch(p) + &
    this%transfer_deadcroot_gr_patch(p) + &
    this%cpool_froot_storage_gr_patch(p) + &
    this%cpool_livecroot_storage_gr_patch(p) + &
    this%cpool_deadcroot_storage_gr_patch(p)
  enddo  
    call p2c(bounds, num_soilc, filter_soilc, &
         this%rr_patch(bounds%begp:bounds%endp), &
         this%rr_col(bounds%begc:bounds%endc))
         
  end subroutine summary_rr


    !-----------------------------------------------------------------------

  subroutine summary_cflux_for_ch4( this, bounds, num_soilp, filter_soilp, num_soilc, filter_soilc )

  !summarize heterotrophic respiration for methane calculation
  !
    use tracer_varcon    , only : is_active_betr_bgc
    use clm_varpar       , only : nlevdecomp, ndecomp_pools, ndecomp_cascade_transitions
  ! !ARGUMENTS:
    class(carbonflux_type) :: this
    type(bounds_type), intent(in)  :: bounds
    integer                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    integer  :: c,p,j,k,l       ! indices
    integer  :: fp,fc           ! lake filter indices

    associate(&
         is_litter =>    decomp_cascade_con%is_litter , & ! Input:  [logical (:) ]  TRUE => pool is a litter pool
         is_soil   =>    decomp_cascade_con%is_soil     & ! Input:  [logical (:) ]  TRUE => pool is a soil pool
    )

    ! patch loop
    do fp = 1,num_soilp
       p = filter_soilp(fp)   

       ! aboveground NPP: leaf, live stem, dead stem (AGNPP)
       ! This is supposed to correspond as closely as possible to
       ! field measurements of AGNPP, so it ignores the storage pools
       ! and only treats the fluxes into displayed pools.

       this%agnpp_patch(p) = &
            this%cpool_to_leafc_patch(p)                  + &
            this%leafc_xfer_to_leafc_patch(p)             + &
            this%cpool_to_livestemc_patch(p)              + &
            this%livestemc_xfer_to_livestemc_patch(p)     + &
            this%cpool_to_deadstemc_patch(p)              + &
            this%deadstemc_xfer_to_deadstemc_patch(p)

       if ( crop_prog .and. pft%itype(p) >= npcropmin )then
          this%agnpp_patch(p) =                    &
               this%agnpp_patch(p)               + &
               this%cpool_to_grainc_patch(p)     + &
               this%grainc_xfer_to_grainc_patch(p)
       endif

       ! belowground NPP: fine root, live coarse root, dead coarse root (BGNPP)
       ! This is supposed to correspond as closely as possible to
       ! field measurements of BGNPP, so it ignores the storage pools
       ! and only treats the fluxes into displayed pools.
       this%bgnpp_patch(p) = &
            this%cpool_to_frootc_patch(p)                   + &
            this%frootc_xfer_to_frootc_patch(p)             + &
            this%cpool_to_livecrootc_patch(p)               + &
            this%livecrootc_xfer_to_livecrootc_patch(p)     + &
            this%cpool_to_deadcrootc_patch(p)               + &
            this%deadcrootc_xfer_to_deadcrootc_patch(p)
    enddo
    ! some zeroing
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%somhr_col(c)              = 0._r8
       this%lithr_col(c)              = 0._r8
       this%decomp_cascade_hr_col(c,1:ndecomp_cascade_transitions)= 0._r8
       this%hr_vr_col(c,1:nlevdecomp) = 0._r8
    enddo

    if ( (.not. is_active_betr_bgc           ) .and. &
         (.not. (use_pflotran .and. pf_cmode))) then
      ! vertically integrate HR and decomposition cascade fluxes
      do k = 1, ndecomp_cascade_transitions

       do j = 1,nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)

             this%decomp_cascade_hr_col(c,k) = &
                this%decomp_cascade_hr_col(c,k) + &
                this%decomp_cascade_hr_vr_col(c,j,k) * dzsoi_decomp(j)

          end do
       end do
      end do

      ! litter heterotrophic respiration (LITHR)
      do k = 1, ndecomp_cascade_transitions
        if ( is_litter(decomp_cascade_con%cascade_donor_pool(k)) ) then
          do fc = 1,num_soilc
            c = filter_soilc(fc)
            this%lithr_col(c) = &
              this%lithr_col(c) + &
              this%decomp_cascade_hr_col(c,k)
          end do
        end if
      end do

      ! soil organic matter heterotrophic respiration (SOMHR)
      do k = 1, ndecomp_cascade_transitions
        if ( is_soil(decomp_cascade_con%cascade_donor_pool(k)) ) then
          do fc = 1,num_soilc
            c = filter_soilc(fc)
            this%somhr_col(c) = &
              this%somhr_col(c) + &
              this%decomp_cascade_hr_col(c,k)
          end do
        end if
      end do

      ! total heterotrophic respiration, vertically resolved (HR)

      do k = 1, ndecomp_cascade_transitions
        do j = 1,nlevdecomp
          do fc = 1,num_soilc
            c = filter_soilc(fc)
            this%hr_vr_col(c,j) = &
                this%hr_vr_col(c,j) + &
                this%decomp_cascade_hr_vr_col(c,j,k)
          end do
        end do
      end do
    endif

    end associate
  end subroutine summary_cflux_for_ch4
  
  
end module CNCarbonFluxType
