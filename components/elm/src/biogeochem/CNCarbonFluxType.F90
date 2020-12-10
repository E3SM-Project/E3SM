module CNCarbonFluxType

  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_infnan_mod         , only : nan => shr_infnan_nan, assignment(=)
  use decompMod              , only : bounds_type
  use elm_varpar             , only : ndecomp_cascade_transitions, ndecomp_pools, nlevcan
  use elm_varpar             , only : crop_prog
  use elm_varpar             , only : nlevdecomp_full, nlevgrnd, nlevdecomp
  use elm_varcon             , only : spval, ispval, dzsoi_decomp
  use landunit_varcon        , only : istsoil, istcrop, istdlak 
  use elm_varctl             , only : use_c13, use_fates 
  use CH4varcon              , only : allowlakeprod
  use pftvarcon              , only : npcropmin
  use CNDecompCascadeConType , only : decomp_cascade_con
  use VegetationType         , only : veg_pp
  use ColumnType             , only : col_pp                
  use LandunitType           , only : lun_pp
  use elm_varctl             , only : nu_com
  use elm_varctl             , only : use_elm_interface, use_pflotran, pf_cmode, use_vertsoilc
  use AnnualFluxDribbler     , only : annual_flux_dribbler_type, annual_flux_dribbler_gridcell
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
     real(r8), pointer :: m_cpool_to_litter_patch                   (:)     ! plant storage C pool to litter (gC/m2/s)

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
     real(r8), pointer :: hrv_cpool_to_litter_patch                 (:)     ! Harvest cpool to litter (gC/m2/s)
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
     real(r8), pointer :: m_cpool_to_fire_patch                     (:)     ! (gC/m2/s) fire C emissions from cpool

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
     real(r8), pointer :: m_cpool_to_litter_fire_patch              (:)     ! (gC/m2/s) from cpool to litter C due to fire              

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
     !turnover of excess carbon
     real(r8), pointer :: xr_patch                                  (:)     ! respiration from excess carbon cpool (gC/m2/s)

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
     real(r8), pointer :: xsmrpool_turnover_patch                   (:)     ! xsmrpool flux to atmosphere due to turnover

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

     ! For comparison with RAINFOR wood productivity data
     real(r8), pointer :: agwdnpp_patch                             (:)     !(gC/m2/s) aboveground NPP


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
     real(r8), pointer :: somc_fire_col                             (:)     ! (gC/m2/s) carbon emissions due to peat burning

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

     ! crop fluxes
     real(r8), pointer :: crop_seedc_to_leaf_patch                  (:)     ! (gC/m2/s) seed source to leaf, for crops

     ! CN dynamic landcover fluxes
     real(r8), pointer :: dwt_seedc_to_leaf_patch                   (:)     ! (gC/m2/s) seed source to patch-level; although this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_seedc_to_leaf_grc                     (:)     ! (gC/m2/s) dwt_seedc_to_leaf_patch summed to the gridcell-level
     real(r8), pointer :: dwt_seedc_to_deadstem_patch               (:)     ! (gC/m2/s) seed source to patch-level; although this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_seedc_to_deadstem_grc                 (:)     ! (gC/m2/s) dwt_seedc_to_leaf_patch summed to the gridcell-level
     real(r8), pointer :: dwt_conv_cflux_patch                      (:)     ! (gC/m2/s) conversion C flux (immediate loss to atm); although this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_conv_cflux_grc                        (:)     ! (gC/m2/s) dwt_conv_cflux_patch summed to the gridcell-level
     real(r8), pointer :: dwt_conv_cflux_dribbled_grc               (:)     ! (gC/m2/s) dwt_conv_cflux_grc dribbled evenly throughout the year
     real(r8), pointer :: dwt_prod10c_gain_patch                    (:)     ! (gC/m2/s) addition to 10-yr wood product pool; although this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_prod100c_gain_patch                   (:)     ! (gC/m2/s) addition to 100-yr wood product pool; although this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_crop_productc_gain_patch              (:)     ! (gC/m2/s) addition to crop product pools from landcover change; although this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_slash_cflux_col                       (:)     ! (gC/m2/s) conversion slash flux due to landcover change
     
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

     real(r8), pointer :: dwt_prod10c_gain_grc                      (:)     ! (gC/m2/s) dynamic landcover addition to 10-year wood product pool
     real(r8), pointer :: dwt_prod100c_gain_grc                     (:)     ! (gC/m2/s) dynamic landcover addition to 100-year wood product pool
     real(r8), pointer :: hrv_deadstemc_to_prod10c_grc              (:)     ! (gC/m2/s) dead stem harvest to 10-year wood product pool
     real(r8), pointer :: hrv_deadstemc_to_prod100c_grc             (:)     ! (gC/m2/s) dead stem harvest to 100-year wood product pool

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
     real(r8), pointer :: fire_decomp_closs_col                     (:)     ! column (gC/m2/s) carbon loss to fire for decomposable pools
     real(r8), pointer :: litfall_col                               (:)     ! column (gC/m2/s) total patch-level litterfall C loss (p2c)       
     real(r8), pointer :: vegfire_col                               (:)     ! column (gC/m2/s) patch-level fire loss (obsolete, mark for removal) (p2c)
     real(r8), pointer :: wood_harvestc_col                         (:)     ! column (p2c)                                                  
     real(r8), pointer :: hrv_xsmrpool_to_atm_col                   (:)     ! column excess MR pool harvest mortality (gC/m2/s) (p2c)
  
     ! Temporary and annual sums
     real(r8), pointer :: tempsum_npp_patch           (:) ! patch temporary annual sum of NPP (gC/m2/yr)
     real(r8), pointer :: annsum_npp_patch            (:) ! patch annual sum of NPP (gC/m2/yr)
     real(r8), pointer :: annsum_npp_col              (:) ! col annual sum of NPP, averaged from pft-level (gC/m2/yr)
     real(r8), pointer :: lag_npp_col                 (:) ! col lagged net primary production (gC/m2/s)
     
     ! debug
     real(r8), pointer :: plant_to_litter_cflux		  (:) ! for the purpose of mass balance check
     real(r8), pointer :: plant_to_cwd_cflux		  (:) ! for the purpose of mass balance check
     real(r8), pointer :: allocation_leaf 		  (:) ! check allocation to leaf for dynamic allocation scheme
     real(r8), pointer :: allocation_stem 		  (:) ! check allocation to stem for dynamic allocation scheme
     real(r8), pointer :: allocation_froot 		  (:) ! check allocation to fine root for dynamic allocation scheme

     ! C4MIP output variable
     real(r8), pointer :: plant_c_to_cwdc                 (:) ! sum of gap, fire, dynamic land use, and harvest mortality, plant carbon flux to CWD

     ! new variables for elm_interface_funcsMod & pflotran
     !------------------------------------------------------------------------
     real(r8), pointer :: externalc_to_decomp_cpools_col            (:,:,:) ! col (gC/m3/s) net C fluxes associated with litter/som-adding/removal to decomp pools
                                                                            ! (sum of all external C additions and removals, excluding decomposition/hr).
     real(r8), pointer :: externalc_to_decomp_delta_col             (:)     ! col (gC/m2) summarized net change of whole column C i/o to decomposing pool bwtn time-step
     real(r8), pointer :: f_co2_soil_vr_col                         (:,:)   ! total vertically-resolved soil-atm. CO2 exchange (gC/m3/s)
     real(r8), pointer :: f_co2_soil_col                            (:)     ! total soil-atm. CO2 exchange (gC/m2/s)
    !------------------------------------------------------------------------

     ! Objects that help convert once-per-year dynamic land cover changes into fluxes
     ! that are dribbled throughout the year
     type(annual_flux_dribbler_type) :: dwt_conv_cflux_dribbler
     type(annual_flux_dribbler_type) :: hrv_xsmrpool_to_atm_dribbler
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
     integer           :: begg,endg
     !------------------------------------------------------------------------

     begp = bounds%begp; endp = bounds%endp
     begc = bounds%begc; endc = bounds%endc
     begg = bounds%begg; endg = bounds%endg

     if (.not.use_fates) then
        allocate(this%m_leafc_to_litter_patch                  (begp:endp)) ; this%m_leafc_to_litter_patch                (:) = nan
        allocate(this%m_frootc_to_litter_patch                 (begp:endp)) ; this%m_frootc_to_litter_patch               (:) = nan
        allocate(this%m_leafc_storage_to_litter_patch          (begp:endp)) ; this%m_leafc_storage_to_litter_patch        (:) = nan
        allocate(this%m_frootc_storage_to_litter_patch         (begp:endp)) ; this%m_frootc_storage_to_litter_patch       (:) = nan
        allocate(this%m_livestemc_storage_to_litter_patch      (begp:endp)) ; this%m_livestemc_storage_to_litter_patch    (:) = nan
        allocate(this%m_deadstemc_storage_to_litter_patch      (begp:endp)) ; this%m_deadstemc_storage_to_litter_patch    (:) = nan
        allocate(this%m_livecrootc_storage_to_litter_patch     (begp:endp)) ; this%m_livecrootc_storage_to_litter_patch   (:) = nan
        allocate(this%m_deadcrootc_storage_to_litter_patch     (begp:endp)) ; this%m_deadcrootc_storage_to_litter_patch   (:) = nan
        allocate(this%m_leafc_xfer_to_litter_patch             (begp:endp)) ; this%m_leafc_xfer_to_litter_patch           (:) = nan
        allocate(this%m_frootc_xfer_to_litter_patch            (begp:endp)) ; this%m_frootc_xfer_to_litter_patch          (:) = nan
        allocate(this%m_livestemc_xfer_to_litter_patch         (begp:endp)) ; this%m_livestemc_xfer_to_litter_patch       (:) = nan
        allocate(this%m_deadstemc_xfer_to_litter_patch         (begp:endp)) ; this%m_deadstemc_xfer_to_litter_patch       (:) = nan
        allocate(this%m_livecrootc_xfer_to_litter_patch        (begp:endp)) ; this%m_livecrootc_xfer_to_litter_patch      (:) = nan
        allocate(this%m_deadcrootc_xfer_to_litter_patch        (begp:endp)) ; this%m_deadcrootc_xfer_to_litter_patch      (:) = nan
        allocate(this%m_livestemc_to_litter_patch              (begp:endp)) ; this%m_livestemc_to_litter_patch            (:) = nan
        allocate(this%m_deadstemc_to_litter_patch              (begp:endp)) ; this%m_deadstemc_to_litter_patch            (:) = nan
        allocate(this%m_livecrootc_to_litter_patch             (begp:endp)) ; this%m_livecrootc_to_litter_patch           (:) = nan
        allocate(this%m_deadcrootc_to_litter_patch             (begp:endp)) ; this%m_deadcrootc_to_litter_patch           (:) = nan
        allocate(this%m_gresp_storage_to_litter_patch          (begp:endp)) ; this%m_gresp_storage_to_litter_patch        (:) = nan
        allocate(this%m_gresp_xfer_to_litter_patch             (begp:endp)) ; this%m_gresp_xfer_to_litter_patch           (:) = nan
        allocate(this%m_cpool_to_litter_patch                  (begp:endp)) ; this%m_cpool_to_litter_patch                (:) = nan
        allocate(this%hrv_leafc_to_litter_patch                (begp:endp)) ; this%hrv_leafc_to_litter_patch              (:) = nan
        allocate(this%hrv_leafc_storage_to_litter_patch        (begp:endp)) ; this%hrv_leafc_storage_to_litter_patch      (:) = nan
        allocate(this%hrv_leafc_xfer_to_litter_patch           (begp:endp)) ; this%hrv_leafc_xfer_to_litter_patch         (:) = nan
        allocate(this%hrv_frootc_to_litter_patch               (begp:endp)) ; this%hrv_frootc_to_litter_patch             (:) = nan
        allocate(this%hrv_frootc_storage_to_litter_patch       (begp:endp)) ; this%hrv_frootc_storage_to_litter_patch     (:) = nan
        allocate(this%hrv_frootc_xfer_to_litter_patch          (begp:endp)) ; this%hrv_frootc_xfer_to_litter_patch        (:) = nan
        allocate(this%hrv_livestemc_to_litter_patch            (begp:endp)) ; this%hrv_livestemc_to_litter_patch          (:) = nan
        allocate(this%hrv_livestemc_storage_to_litter_patch    (begp:endp)) ; this%hrv_livestemc_storage_to_litter_patch  (:) = nan
        allocate(this%hrv_livestemc_xfer_to_litter_patch       (begp:endp)) ; this%hrv_livestemc_xfer_to_litter_patch     (:) = nan
        allocate(this%hrv_deadstemc_to_prod10c_patch           (begp:endp)) ; this%hrv_deadstemc_to_prod10c_patch         (:) = nan
        allocate(this%hrv_deadstemc_to_prod100c_patch          (begp:endp)) ; this%hrv_deadstemc_to_prod100c_patch        (:) = nan
        allocate(this%hrv_leafc_to_prod1c_patch                (begp:endp)) ; this%hrv_leafc_to_prod1c_patch              (:) = nan
        allocate(this%hrv_livestemc_to_prod1c_patch            (begp:endp)) ; this%hrv_livestemc_to_prod1c_patch          (:) = nan
        allocate(this%hrv_grainc_to_prod1c_patch               (begp:endp)) ; this%hrv_grainc_to_prod1c_patch             (:) = nan
        allocate(this%hrv_cropc_to_prod1c_patch                (begp:endp)) ; this%hrv_cropc_to_prod1c_patch              (:) = nan
        allocate(this%hrv_deadstemc_storage_to_litter_patch    (begp:endp)) ; this%hrv_deadstemc_storage_to_litter_patch  (:) = nan
        allocate(this%hrv_deadstemc_xfer_to_litter_patch       (begp:endp)) ; this%hrv_deadstemc_xfer_to_litter_patch     (:) = nan
        allocate(this%hrv_livecrootc_to_litter_patch           (begp:endp)) ; this%hrv_livecrootc_to_litter_patch         (:) = nan
        allocate(this%hrv_livecrootc_storage_to_litter_patch   (begp:endp)) ; this%hrv_livecrootc_storage_to_litter_patch (:) = nan
        allocate(this%hrv_livecrootc_xfer_to_litter_patch      (begp:endp)) ; this%hrv_livecrootc_xfer_to_litter_patch    (:) = nan
        allocate(this%hrv_deadcrootc_to_litter_patch           (begp:endp)) ; this%hrv_deadcrootc_to_litter_patch         (:) = nan
        allocate(this%hrv_deadcrootc_storage_to_litter_patch   (begp:endp)) ; this%hrv_deadcrootc_storage_to_litter_patch (:) = nan
        allocate(this%hrv_deadcrootc_xfer_to_litter_patch      (begp:endp)) ; this%hrv_deadcrootc_xfer_to_litter_patch    (:) = nan
        allocate(this%hrv_gresp_storage_to_litter_patch        (begp:endp)) ; this%hrv_gresp_storage_to_litter_patch      (:) = nan
        allocate(this%hrv_gresp_xfer_to_litter_patch           (begp:endp)) ; this%hrv_gresp_xfer_to_litter_patch         (:) = nan
        allocate(this%hrv_xsmrpool_to_atm_patch                (begp:endp)) ; this%hrv_xsmrpool_to_atm_patch              (:) = nan
        allocate(this%hrv_cpool_to_litter_patch                (begp:endp)) ; this%hrv_cpool_to_litter_patch              (:) = nan
        allocate(this%m_leafc_to_fire_patch                    (begp:endp)) ; this%m_leafc_to_fire_patch                  (:) = nan
        allocate(this%m_leafc_storage_to_fire_patch            (begp:endp)) ; this%m_leafc_storage_to_fire_patch          (:) = nan
        allocate(this%m_leafc_xfer_to_fire_patch               (begp:endp)) ; this%m_leafc_xfer_to_fire_patch             (:) = nan
        allocate(this%m_livestemc_to_fire_patch                (begp:endp)) ; this%m_livestemc_to_fire_patch              (:) = nan
        allocate(this%m_livestemc_storage_to_fire_patch        (begp:endp)) ; this%m_livestemc_storage_to_fire_patch      (:) = nan
        allocate(this%m_livestemc_xfer_to_fire_patch           (begp:endp)) ; this%m_livestemc_xfer_to_fire_patch         (:) = nan
        allocate(this%m_deadstemc_to_fire_patch                (begp:endp)) ; this%m_deadstemc_to_fire_patch              (:) = nan
        allocate(this%m_deadstemc_storage_to_fire_patch        (begp:endp)) ; this%m_deadstemc_storage_to_fire_patch      (:) = nan
        allocate(this%m_deadstemc_xfer_to_fire_patch           (begp:endp)) ; this%m_deadstemc_xfer_to_fire_patch         (:) = nan
        allocate(this%m_frootc_to_fire_patch                   (begp:endp)) ; this%m_frootc_to_fire_patch                 (:) = nan
        allocate(this%m_frootc_storage_to_fire_patch           (begp:endp)) ; this%m_frootc_storage_to_fire_patch         (:) = nan
        allocate(this%m_frootc_xfer_to_fire_patch              (begp:endp)) ; this%m_frootc_xfer_to_fire_patch            (:) = nan
        allocate(this%m_livecrootc_to_fire_patch               (begp:endp)) ; this%m_livecrootc_to_fire_patch             (:) = nan
        allocate(this%m_livecrootc_storage_to_fire_patch       (begp:endp)) ; this%m_livecrootc_storage_to_fire_patch     (:) = nan
        allocate(this%m_livecrootc_xfer_to_fire_patch          (begp:endp)) ; this%m_livecrootc_xfer_to_fire_patch        (:) = nan
        allocate(this%m_deadcrootc_to_fire_patch               (begp:endp)) ; this%m_deadcrootc_to_fire_patch             (:) = nan
        allocate(this%m_deadcrootc_storage_to_fire_patch       (begp:endp)) ; this%m_deadcrootc_storage_to_fire_patch     (:) = nan
        allocate(this%m_deadcrootc_xfer_to_fire_patch          (begp:endp)) ; this%m_deadcrootc_xfer_to_fire_patch        (:) = nan
        allocate(this%m_gresp_storage_to_fire_patch            (begp:endp)) ; this%m_gresp_storage_to_fire_patch          (:) = nan
        allocate(this%m_gresp_xfer_to_fire_patch               (begp:endp)) ; this%m_gresp_xfer_to_fire_patch             (:) = nan
        allocate(this%m_cpool_to_fire_patch                    (begp:endp)) ; this%m_cpool_to_fire_patch                  (:) = nan

        allocate(this%m_leafc_to_litter_fire_patch             (begp:endp)) ; this%m_leafc_to_litter_fire_patch           (:) = nan
        allocate(this%m_leafc_storage_to_litter_fire_patch     (begp:endp)) ; this%m_leafc_storage_to_litter_fire_patch   (:) = nan
        allocate(this%m_leafc_xfer_to_litter_fire_patch        (begp:endp)) ; this%m_leafc_xfer_to_litter_fire_patch      (:) = nan
        allocate(this%m_livestemc_to_litter_fire_patch         (begp:endp)) ; this%m_livestemc_to_litter_fire_patch       (:) = nan
        allocate(this%m_livestemc_storage_to_litter_fire_patch (begp:endp))
        this%m_livestemc_storage_to_litter_fire_patch(:) = nan
        allocate(this%m_livestemc_xfer_to_litter_fire_patch    (begp:endp)) ; this%m_livestemc_xfer_to_litter_fire_patch  (:) = nan
        allocate(this%m_livestemc_to_deadstemc_fire_patch      (begp:endp)) ; this%m_livestemc_to_deadstemc_fire_patch    (:) = nan
        allocate(this%m_deadstemc_to_litter_fire_patch         (begp:endp)) ; this%m_deadstemc_to_litter_fire_patch       (:) = nan
        allocate(this%m_deadstemc_storage_to_litter_fire_patch (begp:endp))
        this%m_deadstemc_storage_to_litter_fire_patch(:) = nan
        allocate(this%m_deadstemc_xfer_to_litter_fire_patch    (begp:endp)) ; this%m_deadstemc_xfer_to_litter_fire_patch  (:) = nan
        allocate(this%m_frootc_to_litter_fire_patch            (begp:endp)) ; this%m_frootc_to_litter_fire_patch          (:) = nan
        allocate(this%m_frootc_storage_to_litter_fire_patch    (begp:endp)) ; this%m_frootc_storage_to_litter_fire_patch  (:) = nan
        allocate(this%m_frootc_xfer_to_litter_fire_patch       (begp:endp)) ; this%m_frootc_xfer_to_litter_fire_patch     (:) = nan
        allocate(this%m_livecrootc_to_litter_fire_patch        (begp:endp)) ; this%m_livecrootc_to_litter_fire_patch      (:) = nan
        allocate(this%m_livecrootc_storage_to_litter_fire_patch(begp:endp)) 
        this%m_livecrootc_storage_to_litter_fire_patch(:) = nan
        allocate(this%m_livecrootc_xfer_to_litter_fire_patch   (begp:endp)) ; this%m_livecrootc_xfer_to_litter_fire_patch (:) = nan
        allocate(this%m_livecrootc_to_deadcrootc_fire_patch    (begp:endp)) ; this%m_livecrootc_to_deadcrootc_fire_patch  (:) = nan
        allocate(this%m_deadcrootc_to_litter_fire_patch        (begp:endp)) ; this%m_deadcrootc_to_litter_fire_patch      (:) = nan
        allocate(this%m_deadcrootc_storage_to_litter_fire_patch(begp:endp))
        this%m_deadcrootc_storage_to_litter_fire_patch (:) = nan
        allocate(this%m_deadcrootc_xfer_to_litter_fire_patch   (begp:endp)) ; this%m_deadcrootc_xfer_to_litter_fire_patch (:) = nan
        allocate(this%m_gresp_storage_to_litter_fire_patch     (begp:endp)) ; this%m_gresp_storage_to_litter_fire_patch   (:) = nan
        allocate(this%m_gresp_xfer_to_litter_fire_patch        (begp:endp)) ; this%m_gresp_xfer_to_litter_fire_patch      (:) = nan
        allocate(this%m_cpool_to_litter_fire_patch             (begp:endp)) ; this%m_cpool_to_litter_fire_patch           (:) = nan

        allocate(this%leafc_xfer_to_leafc_patch                (begp:endp)) ; this%leafc_xfer_to_leafc_patch              (:) = nan
        allocate(this%frootc_xfer_to_frootc_patch              (begp:endp)) ; this%frootc_xfer_to_frootc_patch            (:) = nan
        allocate(this%livestemc_xfer_to_livestemc_patch        (begp:endp)) ; this%livestemc_xfer_to_livestemc_patch      (:) = nan
        allocate(this%deadstemc_xfer_to_deadstemc_patch        (begp:endp)) ; this%deadstemc_xfer_to_deadstemc_patch      (:) = nan
        allocate(this%livecrootc_xfer_to_livecrootc_patch      (begp:endp)) ; this%livecrootc_xfer_to_livecrootc_patch    (:) = nan
        allocate(this%deadcrootc_xfer_to_deadcrootc_patch      (begp:endp)) ; this%deadcrootc_xfer_to_deadcrootc_patch    (:) = nan
        allocate(this%leafc_to_litter_patch                    (begp:endp)) ; this%leafc_to_litter_patch                  (:) = nan
        allocate(this%frootc_to_litter_patch                   (begp:endp)) ; this%frootc_to_litter_patch                 (:) = nan
        allocate(this%leaf_mr_patch                            (begp:endp)) ; this%leaf_mr_patch                          (:) = nan
        allocate(this%froot_mr_patch                           (begp:endp)) ; this%froot_mr_patch                         (:) = nan
        allocate(this%livestem_mr_patch                        (begp:endp)) ; this%livestem_mr_patch                      (:) = nan
        allocate(this%livecroot_mr_patch                       (begp:endp)) ; this%livecroot_mr_patch                     (:) = nan
        allocate(this%grain_mr_patch                           (begp:endp)) ; this%grain_mr_patch                         (:) = nan
        allocate(this%leaf_curmr_patch                         (begp:endp)) ; this%leaf_curmr_patch                       (:) = nan
        allocate(this%froot_curmr_patch                        (begp:endp)) ; this%froot_curmr_patch                      (:) = nan
        allocate(this%livestem_curmr_patch                     (begp:endp)) ; this%livestem_curmr_patch                   (:) = nan
        allocate(this%livecroot_curmr_patch                    (begp:endp)) ; this%livecroot_curmr_patch                  (:) = nan
        allocate(this%grain_curmr_patch                        (begp:endp)) ; this%grain_curmr_patch                      (:) = nan
        allocate(this%leaf_xsmr_patch                          (begp:endp)) ; this%leaf_xsmr_patch                        (:) = nan
        allocate(this%froot_xsmr_patch                         (begp:endp)) ; this%froot_xsmr_patch                       (:) = nan
        allocate(this%livestem_xsmr_patch                      (begp:endp)) ; this%livestem_xsmr_patch                    (:) = nan
        allocate(this%livecroot_xsmr_patch                     (begp:endp)) ; this%livecroot_xsmr_patch                   (:) = nan
        allocate(this%grain_xsmr_patch                         (begp:endp)) ; this%grain_xsmr_patch                       (:) = nan
        allocate(this%xr_patch                                 (begp:endp)) ; this%xr_patch                               (:) = nan
        allocate(this%psnsun_to_cpool_patch                    (begp:endp)) ; this%psnsun_to_cpool_patch                  (:) = nan
        allocate(this%psnshade_to_cpool_patch                  (begp:endp)) ; this%psnshade_to_cpool_patch                (:) = nan
        allocate(this%cpool_to_xsmrpool_patch                  (begp:endp)) ; this%cpool_to_xsmrpool_patch                (:) = nan
        allocate(this%cpool_to_leafc_patch                     (begp:endp)) ; this%cpool_to_leafc_patch                   (:) = nan
        allocate(this%cpool_to_leafc_storage_patch             (begp:endp)) ; this%cpool_to_leafc_storage_patch           (:) = nan
        allocate(this%cpool_to_frootc_patch                    (begp:endp)) ; this%cpool_to_frootc_patch                  (:) = nan
        allocate(this%cpool_to_frootc_storage_patch            (begp:endp)) ; this%cpool_to_frootc_storage_patch          (:) = nan
        allocate(this%cpool_to_livestemc_patch                 (begp:endp)) ; this%cpool_to_livestemc_patch               (:) = nan
        allocate(this%cpool_to_livestemc_storage_patch         (begp:endp)) ; this%cpool_to_livestemc_storage_patch       (:) = nan
        allocate(this%cpool_to_deadstemc_patch                 (begp:endp)) ; this%cpool_to_deadstemc_patch               (:) = nan
        allocate(this%cpool_to_deadstemc_storage_patch         (begp:endp)) ; this%cpool_to_deadstemc_storage_patch       (:) = nan
        allocate(this%cpool_to_livecrootc_patch                (begp:endp)) ; this%cpool_to_livecrootc_patch              (:) = nan
        allocate(this%cpool_to_livecrootc_storage_patch        (begp:endp)) ; this%cpool_to_livecrootc_storage_patch      (:) = nan
        allocate(this%cpool_to_deadcrootc_patch                (begp:endp)) ; this%cpool_to_deadcrootc_patch              (:) = nan
        allocate(this%cpool_to_deadcrootc_storage_patch        (begp:endp)) ; this%cpool_to_deadcrootc_storage_patch      (:) = nan
        allocate(this%cpool_to_gresp_storage_patch             (begp:endp)) ; this%cpool_to_gresp_storage_patch           (:) = nan
        allocate(this%cpool_leaf_gr_patch                      (begp:endp)) ; this%cpool_leaf_gr_patch                    (:) = nan
        allocate(this%cpool_leaf_storage_gr_patch              (begp:endp)) ; this%cpool_leaf_storage_gr_patch            (:) = nan
        allocate(this%transfer_leaf_gr_patch                   (begp:endp)) ; this%transfer_leaf_gr_patch                 (:) = nan
        allocate(this%cpool_froot_gr_patch                     (begp:endp)) ; this%cpool_froot_gr_patch                   (:) = nan
        allocate(this%cpool_froot_storage_gr_patch             (begp:endp)) ; this%cpool_froot_storage_gr_patch           (:) = nan
        allocate(this%transfer_froot_gr_patch                  (begp:endp)) ; this%transfer_froot_gr_patch                (:) = nan
        allocate(this%cpool_livestem_gr_patch                  (begp:endp)) ; this%cpool_livestem_gr_patch                (:) = nan
        allocate(this%cpool_livestem_storage_gr_patch          (begp:endp)) ; this%cpool_livestem_storage_gr_patch        (:) = nan
        allocate(this%transfer_livestem_gr_patch               (begp:endp)) ; this%transfer_livestem_gr_patch             (:) = nan
        allocate(this%cpool_deadstem_gr_patch                  (begp:endp)) ; this%cpool_deadstem_gr_patch                (:) = nan
        allocate(this%cpool_deadstem_storage_gr_patch          (begp:endp)) ; this%cpool_deadstem_storage_gr_patch        (:) = nan
        allocate(this%transfer_deadstem_gr_patch               (begp:endp)) ; this%transfer_deadstem_gr_patch             (:) = nan
        allocate(this%cpool_livecroot_gr_patch                 (begp:endp)) ; this%cpool_livecroot_gr_patch               (:) = nan
        allocate(this%cpool_livecroot_storage_gr_patch         (begp:endp)) ; this%cpool_livecroot_storage_gr_patch       (:) = nan
        allocate(this%transfer_livecroot_gr_patch              (begp:endp)) ; this%transfer_livecroot_gr_patch            (:) = nan
        allocate(this%cpool_deadcroot_gr_patch                 (begp:endp)) ; this%cpool_deadcroot_gr_patch               (:) = nan
        allocate(this%cpool_deadcroot_storage_gr_patch         (begp:endp)) ; this%cpool_deadcroot_storage_gr_patch       (:) = nan
        allocate(this%transfer_deadcroot_gr_patch              (begp:endp)) ; this%transfer_deadcroot_gr_patch            (:) = nan
        allocate(this%leafc_storage_to_xfer_patch              (begp:endp)) ; this%leafc_storage_to_xfer_patch            (:) = nan
        allocate(this%frootc_storage_to_xfer_patch             (begp:endp)) ; this%frootc_storage_to_xfer_patch           (:) = nan
        allocate(this%livestemc_storage_to_xfer_patch          (begp:endp)) ; this%livestemc_storage_to_xfer_patch        (:) = nan
        allocate(this%deadstemc_storage_to_xfer_patch          (begp:endp)) ; this%deadstemc_storage_to_xfer_patch        (:) = nan
        allocate(this%livecrootc_storage_to_xfer_patch         (begp:endp)) ; this%livecrootc_storage_to_xfer_patch       (:) = nan
        allocate(this%deadcrootc_storage_to_xfer_patch         (begp:endp)) ; this%deadcrootc_storage_to_xfer_patch       (:) = nan
        allocate(this%gresp_storage_to_xfer_patch              (begp:endp)) ; this%gresp_storage_to_xfer_patch            (:) = nan
        allocate(this%livestemc_to_deadstemc_patch             (begp:endp)) ; this%livestemc_to_deadstemc_patch           (:) = nan
        allocate(this%livecrootc_to_deadcrootc_patch           (begp:endp)) ; this%livecrootc_to_deadcrootc_patch         (:) = nan
        allocate(this%mr_patch                                 (begp:endp)) ; this%mr_patch                               (:) = nan
        allocate(this%current_gr_patch                         (begp:endp)) ; this%current_gr_patch                       (:) = nan
        allocate(this%transfer_gr_patch                        (begp:endp)) ; this%transfer_gr_patch                      (:) = nan
        allocate(this%storage_gr_patch                         (begp:endp)) ; this%storage_gr_patch                       (:) = nan
        allocate(this%gr_patch                                 (begp:endp)) ; this%gr_patch                               (:) = nan
        allocate(this%ar_patch                                 (begp:endp)) ; this%ar_patch                               (:) = nan
        allocate(this%rr_patch                                 (begp:endp)) ; this%rr_patch                               (:) = nan
        allocate(this%npp_patch                                (begp:endp)) ; this%npp_patch                              (:) = nan
        allocate(this%agnpp_patch                              (begp:endp)) ; this%agnpp_patch                            (:) = nan
        allocate(this%bgnpp_patch                              (begp:endp)) ; this%bgnpp_patch                            (:) = nan
        allocate(this%litfall_patch                            (begp:endp)) ; this%litfall_patch                          (:) = nan
        allocate(this%vegfire_patch                            (begp:endp)) ; this%vegfire_patch                          (:) = nan
        allocate(this%wood_harvestc_patch                      (begp:endp)) ; this%wood_harvestc_patch                    (:) = nan
        allocate(this%cinputs_patch                            (begp:endp)) ; this%cinputs_patch                          (:) = nan
        allocate(this%coutputs_patch                           (begp:endp)) ; this%coutputs_patch                         (:) = nan

        allocate(this%plant_calloc_patch                       (begp:endp)) ; this%plant_calloc_patch                     (:) = nan
        allocate(this%excess_cflux_patch                       (begp:endp)) ; this%excess_cflux_patch                     (:) = nan
        allocate(this%prev_leafc_to_litter_patch               (begp:endp)) ; this%prev_leafc_to_litter_patch             (:) = nan
        allocate(this%prev_frootc_to_litter_patch              (begp:endp)) ; this%prev_frootc_to_litter_patch            (:) = nan
        allocate(this%gpp_patch                                (begp:endp)) ; this%gpp_patch                              (:) = nan
        allocate(this%gpp_before_downreg_patch                 (begp:endp)) ; this%gpp_before_downreg_patch               (:) = nan
        allocate(this%availc_patch                             (begp:endp)) ; this%availc_patch                           (:) = nan
        allocate(this%xsmrpool_recover_patch                   (begp:endp)) ; this%xsmrpool_recover_patch                 (:) = nan
        allocate(this%xsmrpool_c13ratio_patch                  (begp:endp)) ; this%xsmrpool_c13ratio_patch                (:) = nan
        allocate(this%xsmrpool_turnover_patch                  (begp:endp)) ; this%xsmrpool_turnover_patch                (:) = nan

        allocate(this%fire_closs_patch                         (begp:endp)) ; this%fire_closs_patch                       (:) = nan
        allocate(this%cpool_to_grainc_patch                    (begp:endp)) ; this%cpool_to_grainc_patch                  (:) = nan
        allocate(this%cpool_to_grainc_storage_patch            (begp:endp)) ; this%cpool_to_grainc_storage_patch          (:) = nan
        allocate(this%livestemc_to_litter_patch                (begp:endp)) ; this%livestemc_to_litter_patch              (:) = nan
        allocate(this%grainc_to_food_patch                     (begp:endp)) ; this%grainc_to_food_patch                   (:) = nan
        allocate(this%grainc_xfer_to_grainc_patch              (begp:endp)) ; this%grainc_xfer_to_grainc_patch            (:) = nan
        allocate(this%cpool_grain_gr_patch                     (begp:endp)) ; this%cpool_grain_gr_patch                   (:) = nan
        allocate(this%cpool_grain_storage_gr_patch             (begp:endp)) ; this%cpool_grain_storage_gr_patch           (:) = nan
        allocate(this%transfer_grain_gr_patch                  (begp:endp)) ; this%transfer_grain_gr_patch                (:) = nan
        allocate(this%xsmrpool_to_atm_patch                    (begp:endp)) ; this%xsmrpool_to_atm_patch                  (:) = nan
        allocate(this%grainc_storage_to_xfer_patch             (begp:endp)) ; this%grainc_storage_to_xfer_patch           (:) = nan
        allocate(this%frootc_alloc_patch                       (begp:endp)) ; this%frootc_alloc_patch                     (:) = nan
        allocate(this%frootc_loss_patch                        (begp:endp)) ; this%frootc_loss_patch                      (:) = nan
        allocate(this%leafc_alloc_patch                        (begp:endp)) ; this%leafc_alloc_patch                      (:) = nan
        allocate(this%leafc_loss_patch                         (begp:endp)) ; this%leafc_loss_patch                       (:) = nan
        allocate(this%woodc_alloc_patch                        (begp:endp)) ; this%woodc_alloc_patch                      (:) = nan
        allocate(this%woodc_loss_patch                         (begp:endp)) ; this%woodc_loss_patch                       (:) = nan          

        allocate(this%tempavg_agnpp_patch               (begp:endp))                  ; this%tempavg_agnpp_patch (:) = spval
        allocate(this%tempavg_bgnpp_patch               (begp:endp))                  ; this%tempavg_bgnpp_patch (:) = spval
        allocate(this%annavg_agnpp_patch                (begp:endp))                  ; this%annavg_agnpp_patch  (:) = spval ! To detect first year
        allocate(this%annavg_bgnpp_patch                (begp:endp))                  ; this%annavg_bgnpp_patch  (:) = spval ! To detect first year

        allocate(this%agwdnpp_patch                             (begp:endp)) ; this%agwdnpp_patch                          (:) = nan


     end if ! if(.not.use_fates)

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
     allocate(this%crop_seedc_to_leaf_patch          (begp:endp))                  ; this%crop_seedc_to_leaf_patch  (:)  =nan

     allocate(this%dwt_seedc_to_leaf_patch           (begp:endp))                  ; this%dwt_seedc_to_leaf_patch      (:) =nan
     allocate(this%dwt_seedc_to_leaf_grc             (begg:endg))                  ; this%dwt_seedc_to_leaf_grc        (:) =nan
     allocate(this%dwt_seedc_to_deadstem_patch       (begp:endp))                  ; this%dwt_seedc_to_deadstem_patch  (:) =nan
     allocate(this%dwt_seedc_to_deadstem_grc         (begg:endg))                  ; this%dwt_seedc_to_deadstem_grc    (:) =nan
     allocate(this%dwt_conv_cflux_patch              (begp:endp))                  ; this%dwt_conv_cflux_patch         (:) =nan
     allocate(this%dwt_conv_cflux_grc                (begg:endg))                  ; this%dwt_conv_cflux_grc           (:) =nan
     allocate(this%dwt_conv_cflux_dribbled_grc       (begg:endg))                  ; this%dwt_conv_cflux_dribbled_grc  (:) =nan
     allocate(this%dwt_prod10c_gain_patch            (begp:endp))                  ; this%dwt_prod10c_gain_patch       (:) =nan
     allocate(this%dwt_prod100c_gain_patch           (begp:endp))                  ; this%dwt_prod100c_gain_patch      (:) =nan
     allocate(this%dwt_crop_productc_gain_patch      (begp:endp))                  ; this%dwt_crop_productc_gain_patch (:) =nan
     allocate(this%dwt_slash_cflux_col               (begc:endc))                  ; this%dwt_slash_cflux_col          (:) =nan

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

     allocate(this%dwt_prod10c_gain_grc              (begg:endg))                  ; this%dwt_prod10c_gain_grc      (:)  =nan
     allocate(this%dwt_prod100c_gain_grc             (begg:endg))                  ; this%dwt_prod100c_gain_grc     (:)  =nan
     allocate(this%hrv_deadstemc_to_prod10c_grc      (begg:endg))                  ; this%hrv_deadstemc_to_prod10c_grc (:) = nan
     allocate(this%hrv_deadstemc_to_prod100c_grc     (begg:endg))                  ; this%hrv_deadstemc_to_prod100c_grc(:) = nan

     allocate(this%bgc_cpool_ext_inputs_vr_col       (begc:endc, 1:nlevdecomp_full,ndecomp_pools))
     this%bgc_cpool_ext_inputs_vr_col(:,:,:) = nan
     allocate(this%bgc_cpool_ext_loss_vr_col         (begc:endc, 1:nlevdecomp_full,ndecomp_pools))
     this%bgc_cpool_ext_loss_vr_col(:,:,:) = nan

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
     allocate(this%fire_decomp_closs_col             (begc:endc))                  ; this%fire_decomp_closs_col     (:)  =nan
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
     allocate(this%annsum_npp_col        (begc:endc)) ; this%annsum_npp_col        (:) = nan
     allocate(this%lag_npp_col           (begc:endc)) ; this%lag_npp_col           (:) = spval

     ! debug
     allocate(this%plant_to_litter_cflux (begc:endc)) ;	this%plant_to_litter_cflux (:) = nan
     allocate(this%plant_to_cwd_cflux    (begc:endc)) ;	this%plant_to_cwd_cflux	   (:) = nan
     allocate(this%allocation_leaf       (begp:endp)) ; this%allocation_leaf       (:) = nan
     allocate(this%allocation_stem       (begp:endp)) ; this%allocation_stem       (:) = nan
     allocate(this%allocation_froot      (begp:endp)) ; this%allocation_froot      (:) = nan

     ! C4MIP output variable
     allocate(this%plant_c_to_cwdc       (begc:endc)) ; this%plant_c_to_cwdc       (:)  =nan

     ! clm_interface & pflotran
     !------------------------------------------------------------------------
     allocate(this%externalc_to_decomp_cpools_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
     this%externalc_to_decomp_cpools_col(:,:,:) = spval
     allocate(this%externalc_to_decomp_delta_col (begc:endc))
     this%externalc_to_decomp_delta_col (:)     = spval
     allocate(this%f_co2_soil_vr_col             (begc:endc,1:nlevdecomp_full))
     this%f_co2_soil_vr_col             (:,:)   = nan
     allocate(this%f_co2_soil_col                (begc:endc))
     this%f_co2_soil_col                (:)     = nan
     !------------------------------------------------------------------------
  end subroutine InitAllocate; 

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds, carbon_type)
    !
    ! !DESCRIPTION:
    ! add history fields for all CN variables, always set as default='inactive'
    !
    ! !USES:
    use elm_varpar , only : ndecomp_cascade_transitions, ndecomp_pools
    use elm_varpar , only : nlevdecomp, nlevdecomp_full, crop_prog, nlevgrnd
    use elm_varctl , only : hist_wrtch4diag
    use histFileMod, only : hist_addfld1d, hist_addfld2d, hist_addfld_decomp 
    use tracer_varcon    , only : is_active_betr_bgc
    use elm_varctl,  only : get_carbontag
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
    integer           :: begg,endg
    character(24)     :: fieldname
    character(100)    :: longname
    real(r8), pointer :: data1dptr(:)   ! temp. pointer for slicing larger arrays
    real(r8), pointer :: data2dptr(:,:) ! temp. pointer for slicing larger arrays
    character(len=3)  :: ctag
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc
    begg = bounds%begg; endg = bounds%endg

    if (nlevdecomp > 1) then
       vr_suffix = "_vr"
    else 
       vr_suffix = ""
    endif

    !-------------------------------
    ! C flux variables - native to PFT
    !-------------------------------

    ! add history fields for all CLAMP CN variables

    ! ------------------------------------------------------------------------------------
    ! History Diagnostics with FATES turned on is a very limited set, and only
    ! operates on C12 right now.
    ! ------------------------------------------------------------------------------------
    if (use_fates) then
       if (carbon_type == 'c12') then
       end if

       return

    end if



    if (carbon_type == 'c12') then
    end if

    !-------------------------------
    ! C13 flux variables - native to PFT
    !-------------------------------
    if ( carbon_type == 'c13') then

    end if
    !-------------------------------
    ! C14 flux variables - native to PFT
    !-------------------------------
    if ( carbon_type == 'c14' ) then

    end if

    !-------------------------------
    ! C flux variables - native to column 
    !-------------------------------

    ! add history fields for all CLAMP CN variables


    !-------------------------------
    ! C13 flux variables - native to column 
    !-------------------------------


    !-------------------------------
    ! C14 flux variables - native to column 
    !-------------------------------



  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !ARGUMENTS:
    class(carbonflux_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: g, p, c, l, j
    integer :: fc                                        ! filter index
    integer :: num_special_col                           ! number of good values in special_col filter
    integer :: num_special_patch                         ! number of good values in special_patch filter
    integer :: special_col(bounds%endc-bounds%begc+1)    ! special landunit filter - columns
    integer :: special_patch(bounds%endp-bounds%begp+1)  ! special landunit filter - patches
    !-----------------------------------------------------------------------

    ! Set column filters

    num_special_col = 0
    do c = bounds%begc, bounds%endc
       l = col_pp%landunit(c)
       if (lun_pp%ifspecial(l)) then
          num_special_col = num_special_col + 1
          special_col(num_special_col) = c
       end if
    end do

    ! Set patch filters

    num_special_patch = 0
    do p = bounds%begp,bounds%endp
       l = veg_pp%landunit(p)

       if (lun_pp%ifspecial(l)) then
          num_special_patch = num_special_patch + 1
          special_patch(num_special_patch) = p
       end if
    end do

    if (.not.use_fates) then
       
       do g = bounds%begg, bounds%endg
          this%dwt_prod10c_gain_grc(g)          = 0._r8
          this%dwt_prod100c_gain_grc(g)         = 0._r8
          this%hrv_deadstemc_to_prod10c_grc(g)  = 0._r8
          this%hrv_deadstemc_to_prod100c_grc(g) = 0._r8
       end do

       do p = bounds%begp,bounds%endp
          l = veg_pp%landunit(p)

          this%gpp_patch(p)                      = 0._r8
          this%gpp_before_downreg_patch(p)       = 0._r8

          if (lun_pp%ifspecial(l)) then
             this%tempsum_npp_patch(p)           = spval
             this%annsum_npp_patch(p)            = spval
             this%availc_patch(p)                = spval
             this%xsmrpool_recover_patch(p)      = spval
             this%excess_cflux_patch(p)          = spval
             this%plant_calloc_patch(p)          = spval
             this%prev_leafc_to_litter_patch(p)  = spval
             this%prev_frootc_to_litter_patch(p) = spval
             if ( use_c13 ) then
                this%xsmrpool_c13ratio_patch(p)  = spval
             endif
          end if
          if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then
             this%tempsum_npp_patch(p)           = 0._r8
             this%annsum_npp_patch(p)            = 0._r8
             this%availc_patch(p)                = 0._r8
             this%xsmrpool_recover_patch(p)      = 0._r8
             this%excess_cflux_patch(p)          = 0._r8
             this%prev_leafc_to_litter_patch(p)  = 0._r8
             this%prev_frootc_to_litter_patch(p) = 0._r8
             this%plant_calloc_patch(p)          = 0._r8
          end if
       end do

    end if !(.not.use_fates)

    do c = bounds%begc, bounds%endc
       l = col_pp%landunit(c)

       if (lun_pp%ifspecial(l)) then
          this%annsum_npp_col(c) = spval
       end if

       this%fphr_col(c,nlevdecomp+1:nlevgrnd) = 0._r8 !used to be in CH4Mod
       if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then
          this%fphr_col(c,nlevdecomp+1:nlevgrnd) = 0._r8 
       else if (lun_pp%itype(l) == istdlak .and. allowlakeprod) then
          this%fphr_col(c,:) = spval
       else  ! Inactive CH4 columns
          this%fphr_col(c,:) = spval
       end if

       ! also initialize dynamic landcover fluxes so that they have
       ! real values on first timestep, prior to calling pftdyn_cnbal
       if (lun_pp%itype(l) == istsoil .or. lun_pp%itype(l) == istcrop) then
          this%dwt_conv_cflux_col(c)        = 0._r8
          this%dwt_prod10c_gain_col(c)      = 0._r8
          this%dwt_prod100c_gain_col(c)     = 0._r8
          this%prod1c_loss_col(c)           = 0._r8
          this%prod10c_loss_col(c)          = 0._r8
          this%prod100c_loss_col(c)         = 0._r8
          this%dwt_slash_cflux_col(c)       = 0._r8
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
    use elm_varcon       , only : c13ratio, c14ratio
    use elm_varctl       , only : use_lch4, use_betr
    use restUtilMod
    use ncdio_pio

    ! pflotran
!    use elm_varctl       , only : use_pflotran, pf_cmode, use_vertsoilc
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

    
    ! -------------------------------------------
    ! None of these restarts are needed for FATES
    ! -------------------------------------------
    if (use_fates) return

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

    if(.not.use_fates) then
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
          this%m_cpool_to_litter_patch(i)                   = value_patch
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
          this%hrv_cpool_to_litter_patch(i)                 = value_patch

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
          this%m_cpool_to_fire_patch(i)                     = value_patch

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
          this%m_cpool_to_litter_fire_patch(i)              = value_patch

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
          this%xr_patch(i)                                  = value_patch
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
          this%gpp_patch(i)                                 = value_patch
          this%gpp_before_downreg_patch(i)                  = value_patch
          this%mr_patch(i)                                  = value_patch
          this%current_gr_patch(i)                          = value_patch
          this%transfer_gr_patch(i)                         = value_patch
          this%storage_gr_patch(i)                          = value_patch
          this%gr_patch(i)                                  = value_patch
          this%ar_patch(i)                                  = value_patch
          this%rr_patch(i)                                  = value_patch
          this%npp_patch(i)                                 = value_patch 
          this%agnpp_patch(i)                               = value_patch
          this%bgnpp_patch(i)                               = value_patch
          this%agwdnpp_patch(i)                               = value_patch
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
          this%xsmrpool_turnover_patch(i)                   = value_patch
       end do
    end if !(.not.use_fates)

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
          this%crop_seedc_to_leaf_patch(i)      = value_patch
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
    integer  :: g, c, j          ! indices
    !-----------------------------------------------------------------------

    ! set column-level conversion and product pool fluxes
    ! to 0 at the beginning of every timestep

    do g = bounds%begg, bounds%endg
       this%dwt_seedc_to_leaf_grc(g)         = 0._r8
       this%dwt_seedc_to_deadstem_grc(g)     = 0._r8
       this%dwt_conv_cflux_grc(g)            = 0._r8
       this%dwt_prod10c_gain_grc(g)          = 0._r8
       this%dwt_prod100c_gain_grc(g)         = 0._r8
       this%hrv_deadstemc_to_prod10c_grc(g)  = 0._r8
       this%hrv_deadstemc_to_prod100c_grc(g) = 0._r8
    end do
    
    do c = bounds%begc,bounds%endc
       this%dwt_conv_cflux_col(c)           = 0._r8
       this%dwt_prod10c_gain_col(c)         = 0._r8
       this%dwt_prod100c_gain_col(c)        = 0._r8
       this%dwt_slash_cflux_col(c)          = 0._r8
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
    use elm_varctl       , only : iulog
    use clm_time_manager , only : get_step_size
    use elm_varcon       , only : secspday
    use elm_varpar       , only : nlevdecomp, ndecomp_pools, ndecomp_cascade_transitions
    use subgridAveMod    , only : p2c
    use tracer_varcon    , only : is_active_betr_bgc
    use MathfuncMod      , only : dot_sum
    use elm_varpar       , only : nlevdecomp_full
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
    integer  :: nlev
    !-----------------------------------------------------------------------

    associate(& 
         is_litter =>    decomp_cascade_con%is_litter , & ! Input:  [logical (:) ]  TRUE => pool is a litter pool
         is_soil   =>    decomp_cascade_con%is_soil   , & ! Input:  [logical (:) ]  TRUE => pool is a soil pool  
         is_cwd    =>    decomp_cascade_con%is_cwd      & ! Input:  [logical (:) ]  TRUE => pool is a cwd pool   
         )

      ! Note that some of these variables and summary statistics are relevant to fates
      ! yet the great majority are not, and instead of riddling this subroutine
      ! with .not.use_fates filters, a wrapper will be created that selects the variables that should
      ! be used.
      if (use_fates) return


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

       if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
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
       if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
          this%ar_patch(p) = &
               this%mr_patch(p) + &
               this%gr_patch(p) + &
               this%xr_patch(p) + &
               this%xsmrpool_to_atm_patch(p) ! xsmr... is -ve (slevis)
          if (nu_com .ne. 'RD' ) then
             this%ar_patch(p) = this%ar_patch(p) + &
                  this%xsmrpool_turnover_patch(p)
          end if
       else
          this%ar_patch(p) = &
               this%mr_patch(p) + &
               this%gr_patch(p) + &
               this%xr_patch(p)
          if (nu_com .ne. 'RD' ) then
             this%ar_patch(p) = this%ar_patch(p) + &
                  this%xsmrpool_turnover_patch(p)
          end if
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
            this%hrv_gresp_xfer_to_litter_patch(p)            + &
            this%hrv_cpool_to_litter_patch(p)

       ! update the annual litfall accumulator, for use in mortality code

       ! patch-level fire losses (VEGFIRE)
       this%vegfire_patch(p) = 0._r8

       ! patch-level wood harvest
       this%wood_harvestc_patch(p) = &
            this%hrv_deadstemc_to_prod10c_patch(p) + &
            this%hrv_deadstemc_to_prod100c_patch(p)
       if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
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
            this%m_gresp_xfer_to_fire_patch(p)           + &
            this%m_cpool_to_fire_patch(p)

       if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then

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

       if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
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
       if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
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

    ! column soil variables
    ! column variables
    nlev = nlevdecomp
    if (use_pflotran .and. pf_cmode) nlev = nlevdecomp_full

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

       do j = 1,nlev
          do fc = 1,num_soilc
             c = filter_soilc(fc)

                this%decomp_cascade_ctransfer_col(c,k) = &
                     this%decomp_cascade_ctransfer_col(c,k) + &
                     this%decomp_cascade_ctransfer_vr_col(c,j,k) * dzsoi_decomp(j) 
             end do
          end do
       end do


    elseif (is_active_betr_bgc) then

       do fc = 1, num_soilc
          c = filter_soilc(fc)
          this%hr_col(c) = dot_sum(this%hr_vr_col(c,1:nlevdecomp),dzsoi_decomp(1:nlevdecomp)) 
       enddo
    endif
    
    ! some zeroing
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%somhr_col(c)              = 0._r8
       this%lithr_col(c)              = 0._r8
       this%decomp_cascade_hr_col(c,1:ndecomp_cascade_transitions)= 0._r8
       if (.not. (use_pflotran .and. pf_cmode)) then
       ! pflotran has returned 'hr_vr_col(begc:endc,1:nlevdecomp)' to ALM before this subroutine is called in EcosystemDynNoLeaching2
       ! thus 'hr_vr_col' should NOT be set to 0
            this%hr_vr_col(c,1:nlevdecomp) = 0._r8
       end if
    enddo

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
        if ( is_litter(decomp_cascade_con%cascade_donor_pool(k)) .or. is_cwd((decomp_cascade_con%cascade_donor_pool(k)))) then
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

    ! bgc interface & pflotran:
    !----------------------------------------------------------------
    if (use_elm_interface.and. (use_pflotran .and. pf_cmode)) then
        call CSummary_interface(this, bounds, num_soilc, filter_soilc)
    endif
    if(.not. (use_pflotran .and. pf_cmode))then
       ! total heterotrophic respiration (HR)
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%hr_col(c) = &
               this%lithr_col(c) + &
               this%somhr_col(c)
       end do

    end if
    ! CSummary_interface: hr_col(c) will be used below
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
       do j = 1,nlev
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

    ! C4MIP output variable, plant carbon flux to cwd (a part of fVegLitter)
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%plant_c_to_cwdc(c) = 0._r8
       do j = 1, nlevdecomp
          this%plant_c_to_cwdc(c) = this%plant_c_to_cwdc(c) + &
             this%gap_mortality_c_to_cwdc_col(c,j)* dzsoi_decomp(j) + &
             this%fire_mortality_c_to_cwdc_col(c,j)* dzsoi_decomp(j)+ &
             this%harvest_c_to_cwdc_col(c,j)* dzsoi_decomp(j)       + &
             this%dwt_livecrootc_to_cwdc_col(c,j)* dzsoi_decomp(j)  + &
             this%dwt_deadcrootc_to_cwdc_col(c,j)* dzsoi_decomp(j)
       end do
    end do

    if  (.not. is_active_betr_bgc) then

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

       if (.not.(use_pflotran .and. pf_cmode)) then
          ! (LITTERC_LOSS) - litter C loss
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%litterc_loss_col(c) = this%lithr_col(c)
          end do
       end if !(.not.(use_pflotran .and. pf_cmode))

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

       if (use_pflotran .and. pf_cmode) then
          ! note: the follwoing should be useful to non-pflotran-coupled, but seems cause 1 BFB test unmatching.
          ! add up all vertical transport tendency terms and calculate total som leaching loss as the sum of these
          do l = 1, ndecomp_pools
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                this%decomp_cpools_leached_col(c,l) = 0._r8
             end do
             do j = 1, nlev
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
       end if

    end if
    
    ! debug
    do fc = 1,num_soilc
        c = filter_soilc(fc)
        this%plant_to_litter_cflux(c) = 0._r8
        this%plant_to_cwd_cflux(c) = 0._r8
        do j = 1, nlev
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

!-------------------------------------------------------------------------------------------------
! !INTERFACE:
subroutine CSummary_interface(this, bounds, num_soilc, filter_soilc)
!
! !DESCRIPTION:
! bgc interface & pflotran:
! On the radiation time step, perform column-level carbon
! summary calculations, which mainly from PFLOTRAN bgc
!
! !USES:
   use shr_sys_mod, only: shr_sys_flush
   use elm_varpar , only: nlevdecomp_full,ndecomp_pools,ndecomp_cascade_transitions
   use elm_varpar , only: i_met_lit, i_cel_lit, i_lig_lit, i_cwd
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
! LOCAL VARIABLES:
   real(r8) :: dtime                ! time-step (s)
   integer :: c,j,l                 ! indices
   integer :: fc                    ! column filter indices

    associate(&
        is_litter =>    decomp_cascade_con%is_litter , & ! Input:  [logical (:) ]  TRUE => pool is a litter pool
        is_soil   =>    decomp_cascade_con%is_soil   , & ! Input:  [logical (:) ]  TRUE => pool is a soil pool
        is_cwd    =>    decomp_cascade_con%is_cwd      & ! Input:  [logical (:) ]  TRUE => pool is a cwd pool
        )

    dtime = get_step_size()
!---------------------------------------------------------------------------------------------------
   ! total heterotrophic respiration (HR)
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%hr_col(c) = 0._r8
          do j = 1,nlevdecomp_full
             this%hr_col(c) = this%hr_col(c) + &
                this%hr_vr_col(c,j) * dzsoi_decomp(j)
          end do
       end do

       ! new variable to account for co2 exchange (not all HR goes to atm at current time-step)
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%f_co2_soil_col(c) = 0._r8
       end do
       do j = 1,nlevdecomp_full
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
                do j = 1, nlevdecomp_full
                   this%cwdc_loss_col(c) = &
                      this%cwdc_loss_col(c) + &
                      this%decomp_cpools_sourcesink_col(c,j,l) / dtime
                end do
             end do
          end if

          if ( is_litter(l) ) then
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                do j = 1, nlevdecomp_full
                   this%litterc_loss_col(c) = &
                      this%litterc_loss_col(c) + &
                      this%decomp_cpools_sourcesink_col(c,j,l) / dtime
                end do
             end do
          end if

       end do

   ! add up all vertically-resolved addition/removal rates (gC/m3/s) of decomp_pools for PFLOTRAN-bgc
    ! (note: this can be for general purpose, although here added an 'if...endif' block for PF-bgc)
    ! first, need to save the total plant C adding/removing to decomposing pools at previous time-step
    ! for calculating the net changes, which are used to do balance check

    do fc = 1, num_soilc
        c = filter_soilc(fc)
        this%externalc_to_decomp_delta_col(c) = 0._r8
        do l = 1, ndecomp_pools
          do j = 1, nlevdecomp_full
            this%externalc_to_decomp_delta_col(c) = this%externalc_to_decomp_delta_col(c) + &
                                this%externalc_to_decomp_cpools_col(c,j,l)*dzsoi_decomp(j)
          end do
        end do
    end do
    !
    ! do the initialization for the following variable here.
    ! DON'T do so in the beginning of CLM-CN time-step (otherwise the above saved will not work)

    do fc = 1,num_soilc
        c = filter_soilc(fc)
        this%externalc_to_decomp_cpools_col(c, 1:nlevdecomp_full, 1:ndecomp_pools) = 0._r8
    end do

    do fc = 1,num_soilc
       c = filter_soilc(fc)
       do l = 1, ndecomp_pools
          do j = 1, nlevdecomp_full
             ! for litter C pools
             if (l==i_met_lit) then
                this%externalc_to_decomp_cpools_col(c,j,l) =                 &
                    this%externalc_to_decomp_cpools_col(c,j,l)               &
                        + this%phenology_c_to_litr_met_c_col(c,j)            &
                        + this%dwt_frootc_to_litr_met_c_col(c,j)             &
                        + this%gap_mortality_c_to_litr_met_c_col(c,j)        &
                        + this%harvest_c_to_litr_met_c_col(c,j)              &
                        + this%m_c_to_litr_met_fire_col(c,j)                 

             elseif (l==i_cel_lit) then
                this%externalc_to_decomp_cpools_col(c,j,l) =                 &
                    this%externalc_to_decomp_cpools_col(c,j,l)               &
                        + this%phenology_c_to_litr_cel_c_col(c,j)            &
                        + this%dwt_frootc_to_litr_cel_c_col(c,j)             &
                        + this%gap_mortality_c_to_litr_cel_c_col(c,j)        &
                        + this%harvest_c_to_litr_cel_c_col(c,j)              &
                        + this%m_c_to_litr_cel_fire_col(c,j)                 

             elseif (l==i_lig_lit) then
                this%externalc_to_decomp_cpools_col(c,j,l) =                 &
                    this%externalc_to_decomp_cpools_col(c,j,l)               &
                        + this%phenology_c_to_litr_lig_c_col(c,j)            &
                        + this%dwt_frootc_to_litr_lig_c_col(c,j)             &
                        + this%gap_mortality_c_to_litr_lig_c_col(c,j)        &
                        + this%harvest_c_to_litr_lig_c_col(c,j)              &
                        + this%m_c_to_litr_lig_fire_col(c,j)                 

             ! for cwd
             elseif (l==i_cwd) then
                this%externalc_to_decomp_cpools_col(c,j,l) =                 &
                    this%externalc_to_decomp_cpools_col(c,j,l)               &
                        + this%dwt_livecrootc_to_cwdc_col(c,j)               &
                        + this%dwt_deadcrootc_to_cwdc_col(c,j)               &
                        + this%gap_mortality_c_to_cwdc_col(c,j)              &
                        + this%harvest_c_to_cwdc_col(c,j)                    &
                        + this%fire_mortality_c_to_cwdc_col(c,j)             

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
!-------------------------------------------------------------------------------------------------

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
    use elm_varpar       , only : nlevdecomp, ndecomp_pools, ndecomp_cascade_transitions
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
         is_soil   =>    decomp_cascade_con%is_soil   , & ! Input:  [logical (:) ]  TRUE => pool is a soil pool
         is_cwd    =>    decomp_cascade_con%is_cwd      &
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

       if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
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

       this%agwdnpp_patch(p) = &
            this%cpool_to_livestemc_patch(p)              + &
            this%livestemc_xfer_to_livestemc_patch(p)     + &
            this%cpool_to_deadstemc_patch(p)              + &
            this%deadstemc_xfer_to_deadstemc_patch(p)

    enddo
    ! some zeroing
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%somhr_col(c)              = 0._r8
       this%lithr_col(c)              = 0._r8
       this%decomp_cascade_hr_col(c,1:ndecomp_cascade_transitions)= 0._r8
       if (.not. (use_pflotran .and. pf_cmode)) then
       ! pflotran has returned 'hr_vr_col(begc:endc,1:nlevdecomp)' to ALM before this subroutine is called in EcosystemDynNoLeaching2
       ! thus 'hr_vr_col' should NOT be set to 0
            this%hr_vr_col(c,1:nlevdecomp) = 0._r8
       end if
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
        if ( is_litter(decomp_cascade_con%cascade_donor_pool(k)) .or. is_cwd((decomp_cascade_con%cascade_donor_pool(k)))) then
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
