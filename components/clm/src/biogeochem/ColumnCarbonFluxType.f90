module ColumnCarbonFluxType
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
  use ColumnType             , only : col                
  use LandunitType           , only : lun
  use clm_varctl             , only : nu_com
  ! bgc interface & pflotran
  use clm_varctl             , only : use_bgc_interface, use_pflotran, pf_cmode, use_vertsoilc
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

	type, public :: soilcol_carbon_flux
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

  


  contains
      procedure, public  :: Init => init_col_cf
      procedure, public  :: Clean => clean_col_cf
      procedure, public  :: Restart => restart_col_cf
      procedure, private :: InitAllocate => initallocate_col_cf
      procedure, private :: InitHistory => inithistory_col_cf
      procedure, private :: InitCold => initcold_col_cf

  end type soilcol_carbon_flux

  ! Expose public interface 
  type(soilcol_carbon_flux), public, target :: col_cf

! Subroutines
  subroutine init_col_cf(this, bounds, carbon_type)

     class(soilcol_carbon_flux) :: this
     type(bounds_type), intent(in) :: bounds  
     character(len=3) , intent(in) :: carbon_type ! one of ['c12', c13','c14']

     call this%InitAllocate ( bounds)
     call this%InitHistory ( bounds, carbon_type )
     call this%InitCold (bounds )

   end subroutine init_cf_col

! Finished
	subroutine initallocate_col_cf(this, begc, endc)
    class(soilcol_carbon_flux) :: this
    integer, intent(in) :: begc   ! beginning soil column index
    integer, intent(in) :: endc   ! ending soil column index    

		allocate(this%t_scalar_col                      (begc:endc,1:nlevdecomp_full)); this%t_scalar_col (:,:)=spval
		allocate(this%w_scalar_col                      (begc:endc,1:nlevdecomp_full)); this%w_scalar_col (:,:)=spval
		allocate(this%o_scalar_col                      (begc:endc,1:nlevdecomp_full)); this%o_scalar_col (:,:)=spval
		allocate(this%phenology_c_to_litr_met_c_col     (begc:endc,1:nlevdecomp_full)); this%phenology_c_to_litr_met_c_col (:,:)=nan
		allocate(this%phenology_c_to_litr_cel_c_col     (begc:endc,1:nlevdecomp_full)); this%phenology_c_to_litr_cel_c_col (:,:)=nan
		allocate(this%phenology_c_to_litr_lig_c_col     (begc:endc,1:nlevdecomp_full)); this%phenology_c_to_litr_lig_c_col (:,:)=nan
		allocate(this%gap_mortality_c_to_litr_met_c_col (begc:endc,1:nlevdecomp_full)); this%gap_mortality_c_to_litrmet_c_col(:,:)=nan
		allocate(this%gap_mortality_c_to_litr_cel_c_col (begc:endc,1:nlevdecomp_full)); this%gap_mortality_c_to_litrcel_c_col(:,:)=nan
		allocate(this%gap_mortality_c_to_litr_lig_c_col (begc:endc,1:nlevdecomp_full)); this%gap_mortality_c_to_litrlig_c_col(:,:)=nan
		allocate(this%gap_mortality_c_to_cwdc_col       (begc:endc,1:nlevdecomp_full)); this%gap_mortality_c_to_cwdccol  (:,:)=nan
		allocate(this%fire_mortality_c_to_cwdc_col      (begc:endc,1:nlevdecomp_full)); this%fire_mortality_c_to_cwd_col (:,:)=nan
		allocate(this%m_c_to_litr_met_fire_col          (begc:endc,1:nlevdecomp_full)); this%m_c_to_litr_met_fire_co     (:,:)=nan
		allocate(this%m_c_to_litr_cel_fire_col          (begc:endc,1:nlevdecomp_full)); this%m_c_to_litr_cel_fire_co     (:,:)=nan
		allocate(this%m_c_to_litr_lig_fire_col          (begc:endc,1:nlevdecomp_full)); this%m_c_to_litr_lig_fire_co     (:,:)=nan
		allocate(this%harvest_c_to_litr_met_c_col       (begc:endc,1:nlevdecomp_full)); this%harvest_c_to_litr_met_ccol  (:,:)=nan
		allocate(this%harvest_c_to_litr_cel_c_col       (begc:endc,1:nlevdecomp_full)); this%harvest_c_to_litr_cel_ccol  (:,:)=nan
		allocate(this%harvest_c_to_litr_lig_c_col       (begc:endc,1:nlevdecomp_full)); this%harvest_c_to_litr_lig_ccol  (:,:)=nan
		allocate(this%harvest_c_to_cwdc_col             (begc:endc,1:nlevdecomp_full)); this%harvest_c_to_cwdc_col       (:,:)=nan
		allocate(this%phr_vr_col                        (begc:endc,1:nlevdecomp_full)); this%phr_vr_col                  (:,:)=nan 
		allocate(this%fphr_col                          (begc:endc,1:nlevgrnd))       ; this%fphr_col                    (:,:)=nan 
		allocate(this%dwt_frootc_to_litr_met_c_col      (begc:endc,1:nlevdecomp_full)); this%dwt_frootc_to_litr_met__col (:,:)=nan
		allocate(this%dwt_frootc_to_litr_cel_c_col      (begc:endc,1:nlevdecomp_full)); this%dwt_frootc_to_litr_cel__col (:,:)=nan
		allocate(this%dwt_frootc_to_litr_lig_c_col      (begc:endc,1:nlevdecomp_full)); this%dwt_frootc_to_litr_lig__col (:,:)=nan
		allocate(this%dwt_livecrootc_to_cwdc_col        (begc:endc,1:nlevdecomp_full)); this%dwt_livecrootc_to_cwdc_ol   (:,:)=nan
		allocate(this%dwt_deadcrootc_to_cwdc_col        (begc:endc,1:nlevdecomp_full)); this%dwt_deadcrootc_to_cwdc_ol   (:,:)=nan
		allocate(this%dwt_closs_col                     (begc:endc))                  ; this%dwt_closs_col            (:)  =nan
		allocate(this%dwt_seedc_to_leaf_col             (begc:endc))                  ; this%dwt_seedc_to_leaf_col    (:)  =nan
		allocate(this%dwt_seedc_to_deadstem_col         (begc:endc))                  ; this%dwt_seedc_to_deadstem_cl (:)  =nan
		allocate(this%dwt_conv_cflux_col                (begc:endc))                  ; this%dwt_conv_cflux_col       (:)  =nan
		allocate(this%dwt_prod10c_gain_col              (begc:endc))                  ; this%dwt_prod10c_gain_col     (:)  =nan
		allocate(this%dwt_prod100c_gain_col             (begc:endc))                  ; this%dwt_prod100c_gain_col    (:)  =nan
		allocate(this%som_c_leached_col                 (begc:endc))                  ; this%som_c_leached_col        (:)  =nan
		allocate(this%somc_fire_col                     (begc:endc))                  ; this%somc_fire_col            (:)  =nan
		allocate(this%landuseflux_col                   (begc:endc))                  ; this%landuseflux_col          (:)  =nan
		allocate(this%landuptake_col                    (begc:endc))                  ; this%landuptake_col           (:)  =nan
		allocate(this%prod1c_loss_col                   (begc:endc))                  ; this%prod1c_loss_col          (:)  =nan
		allocate(this%prod10c_loss_col                  (begc:endc))                  ; this%prod10c_loss_col         (:)  =nan
		allocate(this%prod100c_loss_col                 (begc:endc))                  ; this%prod100c_loss_col        (:)  =nan
		allocate(this%product_closs_col                 (begc:endc))                  ; this%product_closs_col        (:)  =nan
		allocate(this%bgc_cpool_ext_inputs_vr_col       (begc:endc, 1:nlevdecomp_full,ndecomp_pools));this%bgc_cpoolext_inputs_vr_col (:,:,:) = nan
		allocate(this%bgc_cpool_ext_loss_vr_col         (begc:endc, 1:nlevdecomp_full,ndecomp_pools));this%bgc_cpoolext_loss_vr_col   (:,:,:) = nan
		allocate(this%lf_conv_cflux_col                 (begc:endc))                  ; this%lf_conv_cflux_col        (:)  =nan
		allocate(this%lithr_col                         (begc:endc))                  ; this%lithr_col                (:)  =nan
		allocate(this%somhr_col                         (begc:endc))                  ; this%somhr_col                (:)  =nan
		allocate(this%hr_vr_col                         (begc:endc,1:nlevdecomp_full)); this%hr_vr_col                (:,:)=nan
		allocate(this%hr_col                            (begc:endc))                  ; this%hr_col                   (:)  =nan
		allocate(this%sr_col                            (begc:endc))                  ; this%sr_col                   (:)  =nan
		allocate(this%er_col                            (begc:endc))                  ; this%er_col                   (:)  =nan
		allocate(this%litfire_col                       (begc:endc))                  ; this%litfire_col              (:)  =nan
		allocate(this%somfire_col                       (begc:endc))                  ; this%somfire_col              (:)  =nan
		allocate(this%totfire_col                       (begc:endc))                  ; this%totfire_col              (:)  =nan
		allocate(this%nep_col                           (begc:endc))                  ; this%nep_col                  (:)  =nan
		allocate(this%nbp_col                           (begc:endc))                  ; this%nbp_col                  (:)  =nan
		allocate(this%nee_col                           (begc:endc))                  ; this%nee_col                  (:)  =nan
		allocate(this%cwdc_hr_col                       (begc:endc))                  ; this%cwdc_hr_col              (:)  =nan
		allocate(this%cwdc_loss_col                     (begc:endc))                  ; this%cwdc_loss_col            (:)  =nan
		allocate(this%litterc_loss_col                  (begc:endc))                  ; this%litterc_loss_col         (:)  =nan
		allocate(this%rr_col                            (begc:endc))                  ; this%rr_col                   (:)  =nan
		allocate(this%ar_col                            (begc:endc))                  ; this%ar_col                   (:)  =nan
		allocate(this%gpp_col                           (begc:endc))                  ; this%gpp_col                  (:)  =nan
		allocate(this%npp_col                           (begc:endc))                  ; this%npp_col                  (:)  =nan
		allocate(this%fire_closs_p2c_col                (begc:endc))                  ; this%fire_closs_p2c_col       (:)  =nan
		allocate(this%fire_closs_col                    (begc:endc))                  ; this%fire_closs_col           (:)  =nan
		allocate(this%litfall_col                       (begc:endc))                  ; this%litfall_col              (:)  =nan
		allocate(this%vegfire_col                       (begc:endc))                  ; this%vegfire_col              (:)  =nan
		allocate(this%wood_harvestc_col                 (begc:endc))                  ; this%wood_harvestc_col        (:)  =nan
		allocate(this%hrv_xsmrpool_to_atm_col           (begc:endc))                  ; this%hrv_xsmrpool_to_atm_col  (:)  =nan 
		allocate(this%hrv_deadstemc_to_prod10c_col      (begc:endc))                                                   
		this%hrv_deadstemc_to_prod10c_col(:)= nan


		allocate(this%hrv_deadstemc_to_prod100c_col(begc:endc))                                                   
		allocate(this%hrv_cropc_to_prod1c_col(begc:endc))
		allocate(this%m_decomp_cpools_to_fire_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))                
		allocate(this%m_decomp_cpools_to_fire_col(begc:endc,1:ndecomp_pools))                                     
		allocate(this%decomp_cpools_sourcesink_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))                  
		allocate(this%decomp_cascade_hr_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))        
		allocate(this%decomp_cascade_hr_col(begc:endc,1:ndecomp_cascade_transitions))                             
		allocate(this%decomp_cascade_ctransfer_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)) 
		allocate(this%decomp_cascade_ctransfer_col(begc:endc,1:ndecomp_cascade_transitions))                      
		allocate(this%decomp_k_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))                    
		allocate(this%decomp_cpools_leached_col(begc:endc,1:ndecomp_pools))              
		allocate(this%decomp_cpools_transport_tendency_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))          
		this%hrv_deadstemc_to_prod100c_col(:)= nan
		allocate(this%hrv_deadstemc_to_prod100c_col(begc:endc))                                                   
	 allocate(this%hrv_cropc_to_prod1c_col(begc:endc))
	 allocate(this%m_decomp_cpools_to_fire_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))                
	 allocate(this%m_decomp_cpools_to_fire_col(begc:endc,1:ndecomp_pools))                                     
	 allocate(this%decomp_cpools_sourcesink_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))                  
	 allocate(this%decomp_cascade_hr_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))        
	 allocate(this%decomp_cascade_hr_col(begc:endc,1:ndecomp_cascade_transitions))                             
	 allocate(this%decomp_cascade_ctransfer_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)) 
	 allocate(this%decomp_cascade_ctransfer_col(begc:endc,1:ndecomp_cascade_transitions))                      
	 allocate(this%decomp_k_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))                    
	 allocate(this%decomp_cpools_leached_col(begc:endc,1:ndecomp_pools))              
	 allocate(this%decomp_cpools_transport_tendency_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))          
	 this%hrv_cropc_to_prod1c_col(:) = nan
	 allocate(this%hrv_deadstemc_to_prod100c_col(begc:endc))                                                   
	 allocate(this%hrv_cropc_to_prod1c_col(begc:endc))
	 allocate(this%m_decomp_cpools_to_fire_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))                
	 allocate(this%m_decomp_cpools_to_fire_col(begc:endc,1:ndecomp_pools))                                     
	 allocate(this%decomp_cpools_sourcesink_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))                  
	 allocate(this%decomp_cascade_hr_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))        
	 allocate(this%decomp_cascade_hr_col(begc:endc,1:ndecomp_cascade_transitions))                             
	 allocate(this%decomp_cascade_ctransfer_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)) 
	 allocate(this%decomp_cascade_ctransfer_col(begc:endc,1:ndecomp_cascade_transitions))                      
	 allocate(this%decomp_k_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))                    
	 allocate(this%decomp_cpools_leached_col(begc:endc,1:ndecomp_pools))              
	 allocate(this%decomp_cpools_transport_tendency_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))          
	 this%m_decomp_cpools_to_fire_vr_col(:,:,:)= nan

	 allocate(this%hrv_deadstemc_to_prod100c_col(begc:endc))
	 allocate(this%hrv_cropc_to_prod1c_col(begc:endc))
	 allocate(this%m_decomp_cpools_to_fire_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
	 allocate(this%m_decomp_cpools_to_fire_col(begc:endc,1:ndecomp_pools))
	 allocate(this%decomp_cpools_sourcesink_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
	 allocate(this%decomp_cascade_hr_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
	 allocate(this%decomp_cascade_hr_col(begc:endc,1:ndecomp_cascade_transitions))
	 allocate(this%decomp_cascade_ctransfer_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
	 allocate(this%decomp_cascade_ctransfer_col(begc:endc,1:ndecomp_cascade_transitions))                     
	 allocate(this%decomp_k_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))                   
	 allocate(this%decomp_cpools_leached_col(begc:endc,1:ndecomp_pools))              
	 allocate(this%decomp_cpools_transport_tendency_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))          
	 this%m_decomp_cpools_to_fire_col(:,:)= nan

	 allocate(this%hrv_deadstemc_to_prod100c_col(begc:endc))
	 allocate(this%hrv_cropc_to_prod1c_col(begc:endc))
	 allocate(this%m_decomp_cpools_to_fire_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
	 allocate(this%m_decomp_cpools_to_fire_col(begc:endc,1:ndecomp_pools))
	 allocate(this%decomp_cpools_sourcesink_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
	 allocate(this%decomp_cascade_hr_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
	 allocate(this%decomp_cascade_hr_col(begc:endc,1:ndecomp_cascade_transitions))
	 allocate(this%decomp_cascade_ctransfer_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))
	 allocate(this%decomp_cascade_ctransfer_col(begc:endc,1:ndecomp_cascade_transitions))                     
	 allocate(this%decomp_k_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))                   
	 allocate(this%decomp_cpools_leached_col(begc:endc,1:ndecomp_pools))              
	 allocate(this%decomp_cpools_transport_tendency_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))          
	 this%decomp_cpools_sourcesink_col(:,:,:)= nan

	 allocate(this%hrv_deadstemc_to_prod100c_col(begc:endc))
	 allocate(this%hrv_cropc_to_prod1c_col(begc:endc))
	 allocate(this%m_decomp_cpools_to_fire_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
	 allocate(this%m_decomp_cpools_to_fire_col(begc:endc,1:ndecomp_pools))
	 allocate(this%decomp_cpools_sourcesink_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
	 allocate(this%decomp_cascade_hr_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))        
	 allocate(this%decomp_cascade_hr_col(begc:endc,1:ndecomp_cascade_transitions))                             
	 allocate(this%decomp_cascade_ctransfer_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)) 
	 allocate(this%decomp_cascade_ctransfer_col(begc:endc,1:ndecomp_cascade_transitions))                      
	 allocate(this%decomp_k_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))                    
	 allocate(this%decomp_cpools_leached_col(begc:endc,1:ndecomp_pools))              
	 allocate(this%decomp_cpools_transport_tendency_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
	 this%decomp_cascade_hr_vr_col(:,:,:)= spval

	 allocate(this%hrv_deadstemc_to_prod100c_col(begc:endc))
	 allocate(this%hrv_cropc_to_prod1c_col(begc:endc))
	 allocate(this%m_decomp_cpools_to_fire_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
	 allocate(this%m_decomp_cpools_to_fire_col(begc:endc,1:ndecomp_pools))                     
	 allocate(this%decomp_cpools_sourcesink_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))  
	 allocate(this%decomp_cascade_hr_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))        
	 allocate(this%decomp_cascade_hr_col(begc:endc,1:ndecomp_cascade_transitions))                             
	 allocate(this%decomp_cascade_ctransfer_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)) 
	 allocate(this%decomp_cascade_ctransfer_col(begc:endc,1:ndecomp_cascade_transitions))                      
	 allocate(this%decomp_k_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))                    
	 allocate(this%decomp_cpools_leached_col(begc:endc,1:ndecomp_pools))              
	 allocate(this%decomp_cpools_transport_tendency_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))          
	 this%decomp_cascade_hr_col(:,:)= nan

	 allocate(this%hrv_deadstemc_to_prod100c_col(begc:endc))
	 allocate(this%hrv_cropc_to_prod1c_col(begc:endc))
	 allocate(this%m_decomp_cpools_to_fire_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))                
	 allocate(this%m_decomp_cpools_to_fire_col(begc:endc,1:ndecomp_pools))                                     
	 allocate(this%decomp_cpools_sourcesink_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))                  
	 allocate(this%decomp_cascade_hr_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))        
	 allocate(this%decomp_cascade_hr_col(begc:endc,1:ndecomp_cascade_transitions))                             
	 allocate(this%decomp_cascade_ctransfer_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)) 
	 allocate(this%decomp_cascade_ctransfer_col(begc:endc,1:ndecomp_cascade_transitions))                      
	 allocate(this%decomp_k_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))                    
	 allocate(this%decomp_cpools_leached_col(begc:endc,1:ndecomp_pools))              
	 allocate(this%decomp_cpools_transport_tendency_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))          
	 this%decomp_cascade_ctransfer_vr_col(:,:,:)= nan
	 allocate(this%hrv_deadstemc_to_prod100c_col(begc:endc))                                                    
	 allocate(this%hrv_cropc_to_prod1c_col(begc:endc)) 
	 allocate(this%m_decomp_cpools_to_fire_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))                
	 allocate(this%m_decomp_cpools_to_fire_col(begc:endc,1:ndecomp_pools))                                     
	 allocate(this%decomp_cpools_sourcesink_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))                  
	 allocate(this%decomp_cascade_hr_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))        
	 allocate(this%decomp_cascade_hr_col(begc:endc,1:ndecomp_cascade_transitions))                             
	 allocate(this%decomp_cascade_ctransfer_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)) 
	 allocate(this%decomp_cascade_ctransfer_col(begc:endc,1:ndecomp_cascade_transitions))                      
	 allocate(this%decomp_k_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))                    
	 allocate(this%decomp_cpools_leached_col(begc:endc,1:ndecomp_pools))              
	 allocate(this%decomp_cpools_transport_tendency_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))          
	 this%decomp_cascade_ctransfer_col(:,:)= nan
	 allocate(this%hrv_deadstemc_to_prod100c_col(begc:endc))                                                   
	 allocate(this%hrv_cropc_to_prod1c_col(begc:endc))
	 allocate(this%m_decomp_cpools_to_fire_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))                
	 allocate(this%m_decomp_cpools_to_fire_col(begc:endc,1:ndecomp_pools))                                     
	 allocate(this%decomp_cpools_sourcesink_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))                  
	 allocate(this%decomp_cascade_hr_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))        
	 allocate(this%decomp_cascade_hr_col(begc:endc,1:ndecomp_cascade_transitions))                             
	 allocate(this%decomp_cascade_ctransfer_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)) 
	 allocate(this%decomp_cascade_ctransfer_col(begc:endc,1:ndecomp_cascade_transitions))                      
	 allocate(this%decomp_k_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))                    
	 allocate(this%decomp_cpools_leached_col(begc:endc,1:ndecomp_pools))              
	 allocate(this%decomp_cpools_transport_tendency_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))          
	 this%decomp_k_col(:,:,:)= spval
	 allocate(this%hrv_deadstemc_to_prod100c_col(begc:endc))                                                   
	 allocate(this%hrv_cropc_to_prod1c_col(begc:endc))
	 allocate(this%m_decomp_cpools_to_fire_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))                
	 allocate(this%m_decomp_cpools_to_fire_col(begc:endc,1:ndecomp_pools))                                     
	 allocate(this%decomp_cpools_sourcesink_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))                  
	 allocate(this%decomp_cascade_hr_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))        
	 allocate(this%decomp_cascade_hr_col(begc:endc,1:ndecomp_cascade_transitions))                             
	 allocate(this%decomp_cascade_ctransfer_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)) 
	 allocate(this%decomp_cascade_ctransfer_col(begc:endc,1:ndecomp_cascade_transitions))                      
	 allocate(this%decomp_k_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))                    
	 allocate(this%decomp_cpools_leached_col(begc:endc,1:ndecomp_pools))              
	 allocate(this%decomp_cpools_transport_tendency_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))          
	 this%decomp_cpools_leached_col(:,:)= nan
	 allocate(this%hrv_deadstemc_to_prod100c_col(begc:endc))                                                   
	 allocate(this%hrv_cropc_to_prod1c_col(begc:endc))
	 allocate(this%m_decomp_cpools_to_fire_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))                
	 allocate(this%m_decomp_cpools_to_fire_col(begc:endc,1:ndecomp_pools))                                     
	 allocate(this%decomp_cpools_sourcesink_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))                  
	 allocate(this%decomp_cascade_hr_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))        
	 allocate(this%decomp_cascade_hr_col(begc:endc,1:ndecomp_cascade_transitions))                             
	 allocate(this%decomp_cascade_ctransfer_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)) 
	 allocate(this%decomp_cascade_ctransfer_col(begc:endc,1:ndecomp_cascade_transitions))                      
	 allocate(this%decomp_k_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))                    
	 allocate(this%decomp_cpools_leached_col(begc:endc,1:ndecomp_pools))              
	 allocate(this%decomp_cpools_transport_tendency_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))          
	 this%decomp_cpools_transport_tendency_col(:,:,:)= nan

	  end subroutine init_col_cf

! Finished
  subroutine clean_col_cf(this)

		class(soilcol_carbon_flux) :: this

    deallocate(this%t_scalar_col)
    deallocate(this%w_scalar_col)
    deallocate(this%o_scalar_col)

    deallocate( this%phenology_c_to_litr_met_c_co)
    deallocate( this%phenology_c_to_litr_cel_c_co)
    deallocate( this%phenology_c_to_litr_lig_c_co)
    deallocate( this%gap_mortality_c_to_litr_met_)
    deallocate( this%gap_mortality_c_to_litr_cel_)
    deallocate( this%gap_mortality_c_to_litr_lig_)
    deallocate( this%gap_mortality_c_to_cwdc_col )
    deallocate( this%fire_mortality_c_to_cwdc_col)
    deallocate( this%m_c_to_litr_met_fire_col    )
    deallocate( this%m_c_to_litr_cel_fire_col    )
    deallocate( this%m_c_to_litr_lig_fire_col    )
    deallocate( this%harvest_c_to_litr_met_c_col )
    deallocate( this%harvest_c_to_litr_cel_c_col )
    deallocate( this%harvest_c_to_litr_lig_c_col )
    deallocate( this%harvest_c_to_cwdc_col       )
    deallocate( this%phr_vr_col                  )
    deallocate( this%fphr_col                    )
    deallocate( this%dwt_frootc_to_litr_met_c_col)
    deallocate( this%dwt_frootc_to_litr_cel_c_col)
    deallocate( this%dwt_frootc_to_litr_lig_c_col)
    deallocate( this%dwt_livecrootc_to_cwdc_col  )
    deallocate( this%dwt_deadcrootc_to_cwdc_col  )
    deallocate( this%dwt_closs_col               )
    deallocate( this%dwt_seedc_to_leaf_col       )
    deallocate( this%dwt_seedc_to_deadstem_col   )
    deallocate( this%dwt_conv_cflux_col          )
    deallocate( this%dwt_prod10c_gain_col        )
    deallocate( this%dwt_prod100c_gain_col       )
    deallocate( this%som_c_leached_col           )
    deallocate( this%somc_fire_col               )
    deallocate( this%landuseflux_col             )
    deallocate( this%landuptake_col              )
    deallocate( this%prod1c_loss_col             )
    deallocate( this%prod10c_loss_col            )
    deallocate( this%prod100c_loss_col           )
    deallocate( this%product_closs_col           )
    deallocate( this%bgc_cpool_ext_inputs_vr_col )
    deallocate( this%bgc_cpool_ext_loss_vr_col   )
    deallocate( this%lf_conv_cflux_col           )
    deallocate( this%lithr_col                   )
    deallocate( this%somhr_col                   )
    deallocate( this%hr_vr_col                   )
    deallocate( this%hr_col                      )
    deallocate( this%sr_col                      )
    deallocate( this%er_col                      )
    deallocate( this%litfire_col                 )
    deallocate( this%somfire_col                 )
    deallocate( this%totfire_col                 )
    deallocate( this%nep_col                     )
    deallocate( this%nbp_col                     )
    deallocate( this%nee_col                     )
    deallocate( this%cwdc_hr_col                 )
    deallocate( this%cwdc_loss_col               )
    deallocate( this%litterc_loss_col            )
    deallocate( this%rr_col                      )
    deallocate( this%ar_col                      )
    deallocate( this%gpp_col                     )
    deallocate( this%npp_col                     )
    deallocate( this%fire_closs_p2c_col          )
    deallocate( this%fire_closs_col              )
    deallocate( this%litfall_col                 )
    deallocate( this%vegfire_col                 )
    deallocate( this%wood_harvestc_col           )
    deallocate( this%hrv_xsmrpool_to_atm_col     )
    deallocate( this%hrv_deadstemc_to_prod10c_col)

    deallocate( this%hrv_deadstemc_to_prod100c_col)
    deallocate( this%hrv_cropc_to_prod1c_col)
    deallocate( this%m_decomp_cpools_to_fire_vr_col)
    deallocate( this%m_decomp_cpools_to_fire_col   )
    deallocate( this%decomp_cpools_sourcesink_col         )
    deallocate( this%decomp_cascade_hr_vr_col             )
    deallocate( this%decomp_cascade_hr_col                )
    deallocate( this%decomp_cascade_ctransfer_vr_col      )
    deallocate( this%decomp_cascade_ctransfer_col         )
    deallocate( this%decomp_k_col                         )
    deallocate( this%decomp_cpools_leached_col            )
    deallocate( this%decomp_cpools_transport_tendency_col )

! Finished
  end subroutine clean_col_cf

! Finished
  subroutine inithistory_col_cf(this, bounds, carbon_type)
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
    class(soilcol_carbon_flux) :: this    
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
    
    

  end subroutine inithistory_col_cf

! Finished

  subroutine Summary(this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, isotope)
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
    class(soilcol_carbon_flux)                 :: this
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


! Finished
  subroutine initcold_col_cf(this, bounds)
    !
    ! !ARGUMENTS:
    class(soilcol_carbon_flux) :: this
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

    !!!!!!!  KB - Do we need to keep this for the SetValues call below? !!!!!!
    ! Set patch filters

    ! num_special_patch = 0
    ! do p = bounds%begp,bounds%endp
    !    l = pft%landunit(p)

    !    if (lun%ifspecial(l)) then
    !       num_special_patch = num_special_patch + 1
    !       special_patch(num_special_patch) = p
    !    end if
    ! end do

    ! Set column filters

    num_special_col = 0
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%ifspecial(l)) then
          num_special_col = num_special_col + 1
          special_col(num_special_col) = c
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

  end subroutine initcold_cf_col

  ! Finished
  ! KB - Most of this seems to deal with patches. Keep?? !!!!!!!
  subroutine restart_col_cf ( this, bounds, ncid, flag )
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
    class (soilcol_carbon_flux) :: this
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
 

    call restartvar(ncid=ncid, flag=flag, varname='col_lag_npp', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%lag_npp_col) 

    call restartvar(ncid=ncid, flag=flag, varname='cannsum_npp', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%annsum_npp_col) 



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

  end subroutine restart_cf_col

! Finished
  subroutine SetValues ( this, num_patch, filter_patch, value_patch, num_column, filter_column, value_column)
    !
    ! !DESCRIPTION:
    ! Set carbon state fluxes
    !
    ! !ARGUMENTS:
    class (soilcol_carbon_flux) :: this
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
    class(soilcol_carbon_flux) :: this
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
	     class(soilcol_carbon_flux)          :: this
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


  !!!! KB - Keep this? It's all patches... !!!!
 subroutine summary_rr(this, bounds, num_soilp, filter_soilp, num_soilc, filter_soilc)
    !
    ! description
    ! summarize root respiration

    ! use subgridAveMod    , only: p2c
    ! class(soilcol_carbon_flux) :: this  

    ! type(bounds_type), intent(in) :: bounds  
    ! integer, intent(in) :: num_soilp
    ! integer, intent(in) :: filter_soilp(:)
    ! integer, intent(in) :: num_soilc
    ! integer, intent(in) :: filter_soilc(:)
    ! integer :: fp, p
    !  ! patch loop
    ! do fp = 1,num_soilp
    !   p = filter_soilp(fp)  
    !   ! root respiration (RR)
    !   this%rr_patch(p) = &
    !   this%froot_mr_patch(p) + &
    !   this%cpool_froot_gr_patch(p) + &
    !   this%cpool_livecroot_gr_patch(p) + &
    !   this%cpool_deadcroot_gr_patch(p) + &
    !   this%transfer_froot_gr_patch(p) + &
    !   this%transfer_livecroot_gr_patch(p) + &
    !   this%transfer_deadcroot_gr_patch(p) + &
    !   this%cpool_froot_storage_gr_patch(p) + &
    !   this%cpool_livecroot_storage_gr_patch(p) + &
    !   this%cpool_deadcroot_storage_gr_patch(p)
    ! enddo  
    !   call p2c(bounds, num_soilc, filter_soilc, &
    !        this%rr_patch(bounds%begp:bounds%endp), &
    !        this%rr_col(bounds%begc:bounds%endc))

  end subroutine summary_rr


    !-----------------------------------------------------------------------

  subroutine summary_cflux_for_ch4( this, bounds, num_soilp, filter_soilp, num_soilc, filter_soilc )

	  !summarize heterotrophic respiration for methane calculation
	  !
	    use tracer_varcon    , only : is_active_betr_bgc
	    use clm_varpar       , only : nlevdecomp, ndecomp_pools, ndecomp_cascade_transitions
	  ! !ARGUMENTS:
	    class(soilcol_carbon_flux) :: this
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





end module ColumnCarbonFluxType