module FatesHistoryInterfaceMod

  use FatesConstantsMod        , only : r8 => fates_r8
  use FatesConstantsMod        , only : fates_avg_flag_length, fates_short_string_length, fates_long_string_length
  use FatesConstantsMod        , only : itrue,ifalse
  use FatesGlobals             , only : fates_log
  use FatesGlobals             , only : endrun => fates_endrun

  use FatesIODimensionsMod     , only : fates_io_dimension_type
  use FatesIOVariableKindMod   , only : fates_io_variable_kind_type
  use FatesHistoryVariableType , only : fates_history_variable_type
  use FatesInterfaceMod        , only : hlm_hio_ignore_val
  use FatesInterfaceMod        , only : hlm_use_planthydro
  use FatesInterfaceMod        , only : hlm_use_ed_st3
  use FatesInterfaceMod        , only : numpft
  use EDParamsMod              , only : ED_val_comp_excln
  use FatesInterfaceMod        , only : nlevsclass, nlevage

  ! FIXME(bja, 2016-10) need to remove CLM dependancy 
  use EDPftvarcon              , only : EDPftvarcon_inst

  ! CIME Globals
  use shr_log_mod              , only : errMsg => shr_log_errMsg
  use shr_infnan_mod           , only : isnan => shr_infnan_isnan
  use FatesConstantsMod        , only : g_per_kg
  use FatesConstantsMod        , only : ha_per_m2
  use FatesConstantsMod        , only : days_per_sec
  use FatesConstantsMod        , only : sec_per_day
  use FatesConstantsMod        , only : days_per_year
  use FatesConstantsMod        , only : years_per_day

  implicit none

  ! These variables hold the index of the history output structure so we don't
  ! have to constantly do name lookup when we want to populate the dataset
  ! These indices are set during "define_history_vars()" call to "set_history_var()"
  ! during the initialize phase.  Definitions are not provide, for an explanation of
  ! the variable go to its registry.  (IH_ signifies "index history")
  
  ! Indices to 1D Patch variables

  integer, private :: ih_trimming_pa
  integer, private :: ih_area_plant_pa
  integer, private :: ih_area_treespread_pa
  integer, private :: ih_nesterov_fire_danger_pa
  integer, private :: ih_spitfire_ROS_pa
  integer, private :: ih_effect_wspeed_pa
  integer, private :: ih_TFC_ROS_pa
  integer, private :: ih_fire_intensity_pa
  integer, private :: ih_fire_area_pa
  integer, private :: ih_scorch_height_pa
  integer, private :: ih_fire_fuel_bulkd_pa
  integer, private :: ih_fire_fuel_eff_moist_pa
  integer, private :: ih_fire_fuel_sav_pa
  integer, private :: ih_fire_fuel_mef_pa
  integer, private :: ih_sum_fuel_pa
  integer, private :: ih_litter_in_si
  integer, private :: ih_litter_out_pa

  integer, private :: ih_efpot_pa        ! NA
  integer, private :: ih_rb_pa           ! NA

  integer, private :: ih_daily_temp
  integer, private :: ih_daily_rh
  integer, private :: ih_daily_prec
  integer, private :: ih_seed_bank_si
  integer, private :: ih_seeds_in_pa
  integer, private :: ih_seed_decay_pa
  integer, private :: ih_seed_germination_pa
  integer, private :: ih_bstore_pa
  integer, private :: ih_bdead_pa
  integer, private :: ih_balive_pa
  integer, private :: ih_bleaf_pa
  integer, private :: ih_btotal_pa
  integer, private :: ih_npp_pa
  integer, private :: ih_gpp_pa
  integer, private :: ih_aresp_pa
  integer, private :: ih_maint_resp_pa
  integer, private :: ih_growth_resp_pa
  integer, private :: ih_ar_canopy_pa
  integer, private :: ih_gpp_canopy_pa
  integer, private :: ih_ar_understory_pa
  integer, private :: ih_gpp_understory_pa
  integer, private :: ih_canopy_biomass_pa
  integer, private :: ih_understory_biomass_pa
  
  ! Indices to site by size-class by pft variables
  integer, private :: ih_nplant_si_scag
  integer, private :: ih_nplant_canopy_si_scag
  integer, private :: ih_nplant_understory_si_scag
  integer, private :: ih_ddbh_canopy_si_scag
  integer, private :: ih_ddbh_understory_si_scag
  integer, private :: ih_mortality_canopy_si_scag
  integer, private :: ih_mortality_understory_si_scag

  ! Indices to (site) variables
  integer, private :: ih_nep_si
  integer, private :: ih_nep_timeintegrated_si
  integer, private :: ih_npp_timeintegrated_si
  integer, private :: ih_hr_timeintegrated_si
  integer, private :: ih_nbp_si
  integer, private :: ih_npp_si
  integer, private :: ih_fire_c_to_atm_si
  integer, private :: ih_ed_to_bgc_this_edts_si
  integer, private :: ih_ed_to_bgc_last_edts_si
  integer, private :: ih_totecosysc_si
  integer, private :: ih_totecosysc_old_si
  integer, private :: ih_totedc_si
  integer, private :: ih_totedc_old_si
  integer, private :: ih_totbgcc_si
  integer, private :: ih_totbgcc_old_si
  integer, private :: ih_biomass_stock_si
  integer, private :: ih_litter_stock_si
  integer, private :: ih_cwd_stock_si
  integer, private :: ih_cbal_err_fates_si
  integer, private :: ih_cbal_err_bgc_si
  integer, private :: ih_cbal_err_tot_si
  integer, private :: ih_npatches_si
  integer, private :: ih_ncohorts_si
  integer, private :: ih_demotion_carbonflux_si
  integer, private :: ih_promotion_carbonflux_si
  integer, private :: ih_canopy_mortality_carbonflux_si
  integer, private :: ih_understory_mortality_carbonflux_si
  integer, private :: ih_canopy_spread_si
  
  ! Indices to (site x scpf) variables
  integer, private :: ih_nplant_si_scpf
  integer, private :: ih_gpp_si_scpf
  integer, private :: ih_npp_totl_si_scpf
  integer, private :: ih_npp_leaf_si_scpf
  integer, private :: ih_npp_seed_si_scpf
  integer, private :: ih_npp_fnrt_si_scpf
  integer, private :: ih_npp_bgsw_si_scpf
  integer, private :: ih_npp_bgdw_si_scpf
  integer, private :: ih_npp_agsw_si_scpf
  integer, private :: ih_npp_agdw_si_scpf
  integer, private :: ih_npp_stor_si_scpf

  integer, private :: ih_bstor_canopy_si_scpf
  integer, private :: ih_bstor_understory_si_scpf
  integer, private :: ih_bleaf_canopy_si_scpf
  integer, private :: ih_bleaf_understory_si_scpf
  integer, private :: ih_mortality_canopy_si_scpf
  integer, private :: ih_mortality_understory_si_scpf
  integer, private :: ih_nplant_canopy_si_scpf
  integer, private :: ih_nplant_understory_si_scpf
  integer, private :: ih_ddbh_canopy_si_scpf
  integer, private :: ih_ddbh_understory_si_scpf
  integer, private :: ih_gpp_canopy_si_scpf
  integer, private :: ih_gpp_understory_si_scpf
  integer, private :: ih_ar_canopy_si_scpf
  integer, private :: ih_ar_understory_si_scpf

  integer, private :: ih_ddbh_si_scpf
  integer, private :: ih_ba_si_scpf
  integer, private :: ih_m1_si_scpf
  integer, private :: ih_m2_si_scpf
  integer, private :: ih_m3_si_scpf
  integer, private :: ih_m4_si_scpf
  integer, private :: ih_m5_si_scpf
  integer, private :: ih_m6_si_scpf

   !LOGGING , make sure to add ih_m7_si_scpf and hio_m7_si_scpf
  integer, private :: ih_m7_si_scpf  

  integer, private :: ih_ar_si_scpf
  integer, private :: ih_ar_grow_si_scpf
  integer, private :: ih_ar_maint_si_scpf
  integer, private :: ih_ar_darkm_si_scpf
  integer, private :: ih_ar_agsapm_si_scpf
  integer, private :: ih_ar_crootm_si_scpf
  integer, private :: ih_ar_frootm_si_scpf

  ! indices to (site x scls) variables
  integer, private :: ih_ba_si_scls
  integer, private :: ih_nplant_canopy_si_scls
  integer, private :: ih_nplant_understory_si_scls
  integer, private :: ih_mortality_canopy_si_scls
  integer, private :: ih_mortality_understory_si_scls
  integer, private :: ih_demotion_rate_si_scls
  integer, private :: ih_promotion_rate_si_scls
  integer, private :: ih_trimming_canopy_si_scls
  integer, private :: ih_trimming_understory_si_scls
  integer, private :: ih_crown_area_canopy_si_scls
  integer, private :: ih_crown_area_understory_si_scls
  integer, private :: ih_ddbh_canopy_si_scls
  integer, private :: ih_ddbh_understory_si_scls

  ! lots of non-default diagnostics for understanding canopy versus understory carbon balances
  integer, private :: ih_rdark_canopy_si_scls
  integer, private :: ih_livestem_mr_canopy_si_scls
  integer, private :: ih_livecroot_mr_canopy_si_scls
  integer, private :: ih_froot_mr_canopy_si_scls
  integer, private :: ih_resp_g_canopy_si_scls
  integer, private :: ih_resp_m_canopy_si_scls
  integer, private :: ih_leaf_md_canopy_si_scls
  integer, private :: ih_root_md_canopy_si_scls
  integer, private :: ih_carbon_balance_canopy_si_scls
  integer, private :: ih_seed_prod_canopy_si_scls
  integer, private :: ih_dbalivedt_canopy_si_scls
  integer, private :: ih_dbdeaddt_canopy_si_scls
  integer, private :: ih_dbstoredt_canopy_si_scls
  integer, private :: ih_storage_flux_canopy_si_scls
  integer, private :: ih_npp_leaf_canopy_si_scls
  integer, private :: ih_npp_froot_canopy_si_scls
  integer, private :: ih_npp_bsw_canopy_si_scls
  integer, private :: ih_npp_bdead_canopy_si_scls
  integer, private :: ih_npp_bseed_canopy_si_scls
  integer, private :: ih_npp_store_canopy_si_scls

  integer, private :: ih_rdark_understory_si_scls
  integer, private :: ih_livestem_mr_understory_si_scls
  integer, private :: ih_livecroot_mr_understory_si_scls
  integer, private :: ih_froot_mr_understory_si_scls
  integer, private :: ih_resp_g_understory_si_scls
  integer, private :: ih_resp_m_understory_si_scls
  integer, private :: ih_leaf_md_understory_si_scls
  integer, private :: ih_root_md_understory_si_scls
  integer, private :: ih_carbon_balance_understory_si_scls
  integer, private :: ih_seed_prod_understory_si_scls
  integer, private :: ih_dbalivedt_understory_si_scls
  integer, private :: ih_dbdeaddt_understory_si_scls
  integer, private :: ih_dbstoredt_understory_si_scls
  integer, private :: ih_storage_flux_understory_si_scls
  integer, private :: ih_npp_leaf_understory_si_scls
  integer, private :: ih_npp_froot_understory_si_scls
  integer, private :: ih_npp_bsw_understory_si_scls
  integer, private :: ih_npp_bdead_understory_si_scls
  integer, private :: ih_npp_bseed_understory_si_scls
  integer, private :: ih_npp_store_understory_si_scls

  integer, private :: ih_yesterdaycanopylevel_canopy_si_scls
  integer, private :: ih_yesterdaycanopylevel_understory_si_scls

  ! indices to (site x pft) variables
  integer, private :: ih_biomass_si_pft
  integer, private :: ih_leafbiomass_si_pft
  integer, private :: ih_storebiomass_si_pft
  integer, private :: ih_nindivs_si_pft
  integer, private :: ih_recruitment_si_pft
  integer, private :: ih_mortality_si_pft


  ! indices to (site x patch-age) variables
  integer, private :: ih_area_si_age
  integer, private :: ih_lai_si_age
  integer, private :: ih_canopy_area_si_age
  integer, private :: ih_gpp_si_age
  integer, private :: ih_npp_si_age
  integer, private :: ih_ncl_si_age
  integer, private :: ih_npatches_si_age
  integer, private :: ih_zstar_si_age
  integer, private :: ih_biomass_si_age

  ! Indices to hydraulics variables
  
  integer, private :: ih_errh2o_scpf
  integer, private :: ih_tran_scpf
  integer, private :: ih_rootuptake_scpf
  integer, private :: ih_rootuptake01_scpf
  integer, private :: ih_rootuptake02_scpf
  integer, private :: ih_rootuptake03_scpf
  integer, private :: ih_rootuptake04_scpf
  integer, private :: ih_rootuptake05_scpf
  integer, private :: ih_rootuptake06_scpf
  integer, private :: ih_rootuptake07_scpf
  integer, private :: ih_rootuptake08_scpf
  integer, private :: ih_rootuptake09_scpf
  integer, private :: ih_rootuptake10_scpf
  integer, private :: ih_sapflow_scpf
  integer, private :: ih_iterh1_scpf          
  integer, private :: ih_iterh2_scpf           
  integer, private :: ih_supsub_scpf              
  integer, private :: ih_ath_scpf               
  integer, private :: ih_tth_scpf               
  integer, private :: ih_sth_scpf                     
  integer, private :: ih_lth_scpf                     
  integer, private :: ih_awp_scpf                     
  integer, private :: ih_twp_scpf  
  integer, private :: ih_swp_scpf                     
  integer, private :: ih_lwp_scpf                    
  integer, private :: ih_btran_scpf

  ! indices to (site x fuel class) variables
  integer, private :: ih_litter_moisture_si_fuel

  ! indices to (site x cwd size class) variables
  integer, private :: ih_cwd_ag_si_cwdsc
  integer, private :: ih_cwd_bg_si_cwdsc
  integer, private :: ih_cwd_ag_in_si_cwdsc
  integer, private :: ih_cwd_bg_in_si_cwdsc
  integer, private :: ih_cwd_ag_out_si_cwdsc
  integer, private :: ih_cwd_bg_out_si_cwdsc

  ! indices to (site x [canopy layer x leaf layer]) variables
  integer, private :: ih_parsun_z_si_cnlf
  integer, private :: ih_parsha_z_si_cnlf
  integer, private :: ih_laisun_z_si_cnlf
  integer, private :: ih_laisha_z_si_cnlf
  integer, private :: ih_fabd_sun_si_cnlf
  integer, private :: ih_fabd_sha_si_cnlf
  integer, private :: ih_fabi_sun_si_cnlf
  integer, private :: ih_fabi_sha_si_cnlf
  integer, private :: ih_ts_net_uptake_si_cnlf
  integer, private :: ih_year_net_uptake_si_cnlf
  integer, private :: ih_crownarea_si_cnlf


  ! indices to (site x [canopy layer x leaf layer x pft]) variables
  integer, private :: ih_parsun_z_si_cnlfpft
  integer, private :: ih_parsha_z_si_cnlfpft
  integer, private :: ih_laisun_z_si_cnlfpft
  integer, private :: ih_laisha_z_si_cnlfpft
  integer, private :: ih_fabd_sun_si_cnlfpft
  integer, private :: ih_fabd_sha_si_cnlfpft
  integer, private :: ih_fabi_sun_si_cnlfpft
  integer, private :: ih_fabi_sha_si_cnlfpft

  ! indices to (site x canopy layer) variables
  integer, private :: ih_parsun_top_si_can
  integer, private :: ih_parsha_top_si_can
  integer, private :: ih_laisun_top_si_can
  integer, private :: ih_laisha_top_si_can
  integer, private :: ih_fabd_sun_top_si_can
  integer, private :: ih_fabd_sha_top_si_can
  integer, private :: ih_fabi_sun_top_si_can
  integer, private :: ih_fabi_sha_top_si_can
  integer, private :: ih_crownarea_si_can

  ! The number of variable dim/kind types we have defined (static)
  integer, parameter :: fates_history_num_dimensions = 13
  integer, parameter :: fates_history_num_dim_kinds = 15
  

  
  ! This structure is allocated by thread, and must be calculated after the FATES
  ! sites are allocated, and their mapping to the HLM is identified.  This structure
  ! is not combined with iovar_bounds, because that one is multi-instanced.  This
  ! structure is used more during the update phase, wherease _bounds is used
  ! more for things like flushing
  type iovar_map_type
     integer, allocatable :: site_index(:)   ! maps site indexes to the HIO site position
     integer, allocatable :: patch1_index(:) ! maps site index to the HIO patch 1st position
  end type iovar_map_type


  type, public :: fates_history_interface_type
     
     ! Instance of the list of history output varialbes
     type(fates_history_variable_type), allocatable :: hvars(:)
     integer, private :: num_history_vars_
     
     ! Instanteat one registry of the different dimension/kinds (dk)
     ! All output variables will have a pointer to one of these dk's
     type(fates_io_variable_kind_type) :: dim_kinds(fates_history_num_dim_kinds)
     
     ! This is a structure that explains where FATES patch boundaries
     ! on each thread point to in the host IO array, this structure is
     ! allocated by number of threads. This could be dynamically
     ! allocated, but is unlikely to change...?
     type(fates_io_dimension_type) :: dim_bounds(fates_history_num_dimensions)
     
     type(iovar_map_type), pointer :: iovar_map(:)

     integer, private :: patch_index_, column_index_, levgrnd_index_, levscpf_index_
     integer, private :: levscls_index_, levpft_index_, levage_index_
     integer, private :: levfuel_index_, levcwdsc_index_, levscag_index_
     integer, private :: levcan_index_, levcnlf_index_, levcnlfpft_index_
   contains
     
     procedure, public :: Init
     procedure, public :: SetThreadBoundsEach
     procedure, public :: initialize_history_vars
     procedure, public :: assemble_history_output_types
     
     procedure, public :: update_history_dyn
     procedure, public :: update_history_prod
     procedure, public :: update_history_cbal
     procedure, public :: update_history_hydraulics

     ! 'get' methods used by external callers to access private read only data
     procedure, public :: num_history_vars
     procedure, public :: patch_index
     procedure, public :: column_index
     procedure, public :: levgrnd_index
     procedure, public :: levscpf_index
     procedure, public :: levscls_index
     procedure, public :: levpft_index
     procedure, public :: levage_index
     procedure, public :: levfuel_index
     procedure, public :: levcwdsc_index
     procedure, public :: levcan_index
     procedure, public :: levcnlf_index
     procedure, public :: levcnlfpft_index
     procedure, public :: levscag_index

     ! private work functions
     procedure, private :: define_history_vars
     procedure, private :: set_history_var
     procedure, private :: init_dim_kinds_maps
     procedure, private :: set_dim_indices
     procedure, private :: flush_hvars

     procedure, private :: set_patch_index
     procedure, private :: set_column_index
     procedure, private :: set_levgrnd_index
     procedure, private :: set_levscpf_index
     procedure, private :: set_levscls_index
     procedure, private :: set_levpft_index
     procedure, private :: set_levage_index
     procedure, private :: set_levfuel_index
     procedure, private :: set_levcwdsc_index
     procedure, private :: set_levcan_index
     procedure, private :: set_levcnlf_index
     procedure, private :: set_levcnlfpft_index
     procedure, private :: set_levscag_index

  end type fates_history_interface_type
   
  character(len=*), parameter, private :: sourcefile = &
         __FILE__

contains

  ! ======================================================================
  
  subroutine Init(this, num_threads, fates_bounds)

    use FatesIODimensionsMod, only : patch, column, levgrnd, levscpf
    use FatesIODimensionsMod, only : levscls, levpft, levage
    use FatesIODimensionsMod, only : levfuel, levcwdsc, levscag
    use FatesIODimensionsMod, only : levcan, levcnlf, levcnlfpft
    use FatesIODimensionsMod, only : fates_bounds_type

    implicit none

    class(fates_history_interface_type), intent(inout) :: this
    integer, intent(in) :: num_threads
    type(fates_bounds_type), intent(in) :: fates_bounds

    integer :: dim_count = 0

    dim_count = dim_count + 1
    call this%set_patch_index(dim_count)
    call this%dim_bounds(dim_count)%Init(patch, num_threads, &
         fates_bounds%patch_begin, fates_bounds%patch_end)

    dim_count = dim_count + 1
    call this%set_column_index(dim_count)
    call this%dim_bounds(dim_count)%Init(column, num_threads, &
         fates_bounds%column_begin, fates_bounds%column_end)

    dim_count = dim_count + 1
    call this%set_levgrnd_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levgrnd, num_threads, &
         fates_bounds%ground_begin, fates_bounds%ground_end)

    dim_count = dim_count + 1
    call this%set_levscpf_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levscpf, num_threads, &
         fates_bounds%sizepft_class_begin, fates_bounds%sizepft_class_end)

    dim_count = dim_count + 1
    call this%set_levscls_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levscls, num_threads, &
         fates_bounds%size_class_begin, fates_bounds%size_class_end)

    dim_count = dim_count + 1
    call this%set_levpft_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levpft, num_threads, &
         fates_bounds%pft_class_begin, fates_bounds%pft_class_end)

    dim_count = dim_count + 1
    call this%set_levage_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levage, num_threads, &
         fates_bounds%age_class_begin, fates_bounds%age_class_end)

    dim_count = dim_count + 1
    call this%set_levfuel_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levfuel, num_threads, &
         fates_bounds%fuel_begin, fates_bounds%fuel_end)

    dim_count = dim_count + 1
    call this%set_levcwdsc_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levcwdsc, num_threads, &
         fates_bounds%cwdsc_begin, fates_bounds%cwdsc_end)

    dim_count = dim_count + 1
    call this%set_levcan_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levcan, num_threads, &
         fates_bounds%can_begin, fates_bounds%can_end)

    dim_count = dim_count + 1
    call this%set_levcnlf_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levcnlf, num_threads, &
         fates_bounds%cnlf_begin, fates_bounds%cnlf_end)

    dim_count = dim_count + 1
    call this%set_levcnlfpft_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levcnlfpft, num_threads, &
         fates_bounds%cnlfpft_begin, fates_bounds%cnlfpft_end)

    dim_count = dim_count + 1
    call this%set_levscag_index(dim_count)
    call this%dim_bounds(dim_count)%Init(levscag, num_threads, &
         fates_bounds%sizeage_class_begin, fates_bounds%sizeage_class_end)
    

    ! FIXME(bja, 2016-10) assert(dim_count == FatesHistorydimensionmod::num_dimension_types)

    ! Allocate the mapping between FATES indices and the IO indices
    allocate(this%iovar_map(num_threads))
    
  end subroutine Init

  ! ======================================================================
  subroutine SetThreadBoundsEach(this, thread_index, thread_bounds)

    use FatesIODimensionsMod, only : fates_bounds_type

    implicit none

    class(fates_history_interface_type), intent(inout) :: this

    integer, intent(in) :: thread_index
    type(fates_bounds_type), intent(in) :: thread_bounds

    integer :: index
    
    index = this%patch_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%patch_begin, thread_bounds%patch_end)

    index = this%column_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%column_begin, thread_bounds%column_end)

    index = this%levgrnd_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%ground_begin, thread_bounds%ground_end)

    index = this%levscpf_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%sizepft_class_begin, thread_bounds%sizepft_class_end)

    index = this%levscls_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%size_class_begin, thread_bounds%size_class_end)

    index = this%levpft_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%pft_class_begin, thread_bounds%pft_class_end)
    
    index = this%levage_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%age_class_begin, thread_bounds%age_class_end)
    
    index = this%levfuel_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%fuel_begin, thread_bounds%fuel_end)
    
    index = this%levcwdsc_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%cwdsc_begin, thread_bounds%cwdsc_end)
    
    index = this%levcan_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%can_begin, thread_bounds%can_end)
    
    index = this%levcnlf_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
         thread_bounds%cnlf_begin, thread_bounds%cnlf_end)
    
    index = this%levcnlfpft_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
          thread_bounds%cnlfpft_begin, thread_bounds%cnlfpft_end)
    
    index = this%levscag_index()
    call this%dim_bounds(index)%SetThreadBounds(thread_index, &
          thread_bounds%sizeage_class_begin, thread_bounds%sizeage_class_end)
    
  end subroutine SetThreadBoundsEach
  
  ! ===================================================================================
  subroutine assemble_history_output_types(this)

    use FatesIOVariableKindMod, only : patch_r8, patch_ground_r8, patch_size_pft_r8
    use FatesIOVariableKindMod, only : site_r8, site_ground_r8, site_size_pft_r8
    use FatesIOVariableKindMod, only : site_size_r8, site_pft_r8, site_age_r8
    use FatesIOVariableKindMod, only : site_fuel_r8, site_cwdsc_r8, site_scag_r8
    use FatesIOVariableKindMod, only : site_can_r8, site_cnlf_r8, site_cnlfpft_r8

   implicit none

    class(fates_history_interface_type), intent(inout) :: this

    call this%init_dim_kinds_maps()

    call this%set_dim_indices(patch_r8, 1, this%patch_index())

    call this%set_dim_indices(site_r8, 1, this%column_index())

    call this%set_dim_indices(patch_ground_r8, 1, this%patch_index())
    call this%set_dim_indices(patch_ground_r8, 2, this%levgrnd_index())

    call this%set_dim_indices(site_ground_r8, 1, this%column_index())
    call this%set_dim_indices(site_ground_r8, 2, this%levgrnd_index())

    call this%set_dim_indices(patch_size_pft_r8, 1, this%patch_index())
    call this%set_dim_indices(patch_size_pft_r8, 2, this%levscpf_index())

    call this%set_dim_indices(site_size_pft_r8, 1, this%column_index())
    call this%set_dim_indices(site_size_pft_r8, 2, this%levscpf_index())

    call this%set_dim_indices(site_size_r8, 1, this%column_index())
    call this%set_dim_indices(site_size_r8, 2, this%levscls_index())

    call this%set_dim_indices(site_pft_r8, 1, this%column_index())
    call this%set_dim_indices(site_pft_r8, 2, this%levpft_index())

    call this%set_dim_indices(site_age_r8, 1, this%column_index())
    call this%set_dim_indices(site_age_r8, 2, this%levage_index())

    call this%set_dim_indices(site_fuel_r8, 1, this%column_index())
    call this%set_dim_indices(site_fuel_r8, 2, this%levfuel_index())

    call this%set_dim_indices(site_cwdsc_r8, 1, this%column_index())
    call this%set_dim_indices(site_cwdsc_r8, 2, this%levcwdsc_index())

    call this%set_dim_indices(site_can_r8, 1, this%column_index())
    call this%set_dim_indices(site_can_r8, 2, this%levcan_index())

    call this%set_dim_indices(site_cnlf_r8, 1, this%column_index())
    call this%set_dim_indices(site_cnlf_r8, 2, this%levcnlf_index())

    call this%set_dim_indices(site_cnlfpft_r8, 1, this%column_index())
    call this%set_dim_indices(site_cnlfpft_r8, 2, this%levcnlfpft_index())

    call this%set_dim_indices(site_scag_r8, 1, this%column_index())
    call this%set_dim_indices(site_scag_r8, 2, this%levscag_index())

  end subroutine assemble_history_output_types
  
  ! ===================================================================================
  
  subroutine set_dim_indices(this, dk_name, idim, dim_index)

    use FatesIOVariableKindMod , only : iotype_index

    implicit none

    ! arguments
    class(fates_history_interface_type), intent(inout) :: this
    character(len=*), intent(in)     :: dk_name
    integer, intent(in)              :: idim  ! dimension index
    integer, intent(in) :: dim_index


    ! local
    integer :: ityp

    ityp = iotype_index(trim(dk_name), fates_history_num_dim_kinds, this%dim_kinds)

    ! First check to see if the dimension is allocated
    if (this%dim_kinds(ityp)%ndims < idim) then
       write(fates_log(), *) 'Trying to define dimension size to a dim-type structure'
       write(fates_log(), *) 'but the dimension index does not exist'
       write(fates_log(), *) 'type: ',dk_name,' ndims: ',this%dim_kinds(ityp)%ndims,' input dim:',idim
       stop
       !end_run
    end if

    if (idim == 1) then
       this%dim_kinds(ityp)%dim1_index = dim_index
    else if (idim == 2) then
       this%dim_kinds(ityp)%dim2_index = dim_index
    end if

    ! With the map, we can set the dimension size
    this%dim_kinds(ityp)%dimsize(idim) = this%dim_bounds(dim_index)%upper_bound - &
         this%dim_bounds(dim_index)%lower_bound + 1

 end subroutine set_dim_indices
  
 ! =======================================================================
 subroutine set_patch_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%patch_index_ = index
 end subroutine set_patch_index

 integer function patch_index(this)
   implicit none
   class(fates_history_interface_type), intent(in) :: this
   patch_index = this%patch_index_
 end function patch_index

 ! =======================================================================
 subroutine set_column_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%column_index_ = index
 end subroutine set_column_index

 integer function column_index(this)
   implicit none
   class(fates_history_interface_type), intent(in) :: this
   column_index = this%column_index_
 end function column_index

 ! =======================================================================
 subroutine set_levgrnd_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%levgrnd_index_ = index
 end subroutine set_levgrnd_index

 integer function levgrnd_index(this)
   implicit none
   class(fates_history_interface_type), intent(in) :: this
   levgrnd_index = this%levgrnd_index_
 end function levgrnd_index

 ! =======================================================================
 subroutine set_levscpf_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%levscpf_index_ = index
 end subroutine set_levscpf_index

 integer function levscpf_index(this)
   implicit none
   class(fates_history_interface_type), intent(in) :: this
   levscpf_index = this%levscpf_index_
 end function levscpf_index

 ! =======================================================================
 subroutine set_levscls_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%levscls_index_ = index
 end subroutine set_levscls_index

 integer function levscls_index(this)
   implicit none
   class(fates_history_interface_type), intent(in) :: this
   levscls_index = this%levscls_index_
 end function levscls_index

 ! =======================================================================
 subroutine set_levpft_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%levpft_index_ = index
 end subroutine set_levpft_index

 integer function levpft_index(this)
   implicit none
   class(fates_history_interface_type), intent(in) :: this
   levpft_index = this%levpft_index_
 end function levpft_index

 ! =======================================================================
 subroutine set_levage_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%levage_index_ = index
 end subroutine set_levage_index

 integer function levage_index(this)
   implicit none
   class(fates_history_interface_type), intent(in) :: this
   levage_index = this%levage_index_
 end function levage_index

 ! =======================================================================
 subroutine set_levfuel_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%levfuel_index_ = index
 end subroutine set_levfuel_index

 integer function levfuel_index(this)
   implicit none
   class(fates_history_interface_type), intent(in) :: this
   levfuel_index = this%levfuel_index_
 end function levfuel_index

 ! =======================================================================
 subroutine set_levcwdsc_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%levcwdsc_index_ = index
 end subroutine set_levcwdsc_index

 integer function levcwdsc_index(this)
   implicit none
   class(fates_history_interface_type), intent(in) :: this
   levcwdsc_index = this%levcwdsc_index_
 end function levcwdsc_index

 ! =======================================================================
 subroutine set_levcan_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%levcan_index_ = index
 end subroutine set_levcan_index

 integer function levcan_index(this)
   implicit none
   class(fates_history_interface_type), intent(in) :: this
   levcan_index = this%levcan_index_
 end function levcan_index

 ! =======================================================================
 subroutine set_levcnlf_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%levcnlf_index_ = index
 end subroutine set_levcnlf_index

 integer function levcnlf_index(this)
   implicit none
   class(fates_history_interface_type), intent(in) :: this
   levcnlf_index = this%levcnlf_index_
 end function levcnlf_index

 ! =======================================================================
 subroutine set_levcnlfpft_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%levcnlfpft_index_ = index
 end subroutine set_levcnlfpft_index

 integer function levcnlfpft_index(this)
   implicit none
   class(fates_history_interface_type), intent(in) :: this
   levcnlfpft_index = this%levcnlfpft_index_
 end function levcnlfpft_index

 ! ======================================================================================
 subroutine set_levscag_index(this, index)
   implicit none
   class(fates_history_interface_type), intent(inout) :: this
   integer, intent(in) :: index
   this%levscag_index_ = index
 end subroutine set_levscag_index

 integer function levscag_index(this)
    implicit none
    class(fates_history_interface_type), intent(in) :: this
    levscag_index = this%levscag_index_
 end function levscag_index
 ! ======================================================================================


 subroutine flush_hvars(this,nc,upfreq_in)
 
   class(fates_history_interface_type)        :: this
   integer,intent(in)                     :: nc
   integer,intent(in)                     :: upfreq_in
   integer                      :: ivar
   integer                      :: lb1,ub1,lb2,ub2

   do ivar=1,ubound(this%hvars,1)
      if (this%hvars(ivar)%upfreq == upfreq_in) then ! Only flush variables with update on dynamics step
         call this%hvars(ivar)%flush(nc, this%dim_bounds, this%dim_kinds)
         
      end if
   end do
   
end subroutine flush_hvars

  
  ! =====================================================================================
   
  subroutine set_history_var(this, vname, units, long, use_default, avgflag, vtype, &
       hlms, flushval, upfreq, ivar, initialize, index)

    use FatesUtilsMod, only     : check_hlm_list
    use FatesInterfaceMod, only : hlm_name

    implicit none
    
    ! arguments
    class(fates_history_interface_type), intent(inout) :: this
    character(len=*), intent(in)  :: vname
    character(len=*), intent(in)  :: units
    character(len=*), intent(in)  :: long
    character(len=*), intent(in)  :: use_default
    character(len=*), intent(in)  :: avgflag
    character(len=*), intent(in)  :: vtype
    character(len=*), intent(in)  :: hlms
    real(r8), intent(in)          :: flushval ! IF THE TYPE IS AN INT WE WILL round with NINT
    integer, intent(in)           :: upfreq
    logical, intent(in) :: initialize
    integer, intent(inout)       :: ivar
    integer, intent(inout)       :: index  ! This is the index for the variable of
                                           ! interest that is associated with an
                                           ! explict name (for fast reference during update)
                                           ! A zero is passed back when the variable is
                                           ! not used

    ! locals
    integer :: ub1, lb1, ub2, lb2    ! Bounds for allocating the var
    integer :: ityp

    logical :: write_var

    write_var = check_hlm_list(trim(hlms), trim(hlm_name))
    if( write_var ) then
       ivar  = ivar+1
       index = ivar    
       
       if (initialize) then
          call this%hvars(ivar)%Init(vname, units, long, use_default, &
               vtype, avgflag, flushval, upfreq, fates_history_num_dim_kinds, this%dim_kinds, &
               this%dim_bounds)
       end if
    else
       index = 0
    end if
    
    return
  end subroutine set_history_var
  
  ! ====================================================================================
  
  subroutine init_dim_kinds_maps(this)
    
    ! ----------------------------------------------------------------------------------
    ! This subroutine simply initializes the structures that define the different
    ! array and type formats for different IO variables
    !
    ! PA_R8   : 1D patch scale 8-byte reals
    ! SI_R8   : 1D site scale 8-byte reals
    !
    ! The allocation on the structures is not dynamic and should only add up to the
    ! number of entries listed here.
    !
    ! ----------------------------------------------------------------------------------
    use FatesIOVariableKindMod, only : patch_r8, patch_ground_r8, patch_size_pft_r8
    use FatesIOVariableKindMod, only : site_r8, site_ground_r8, site_size_pft_r8
    use FatesIOVariableKindMod, only : site_size_r8, site_pft_r8, site_age_r8
    use FatesIOVariableKindMod, only : site_fuel_r8, site_cwdsc_r8, site_scag_r8
    use FatesIOVariableKindMod, only : site_can_r8, site_cnlf_r8, site_cnlfpft_r8
    
    implicit none
    
    ! Arguments
    class(fates_history_interface_type), intent(inout) :: this
       

    integer :: index

    ! 1d Patch
    index = 1
    call this%dim_kinds(index)%Init(patch_r8, 1)

    ! 1d Site
    index = index + 1
    call this%dim_kinds(index)%Init(site_r8, 1)

    ! patch x ground
    index = index + 1
    call this%dim_kinds(index)%Init(patch_ground_r8, 2)

    ! patch x size-class/pft
    index = index + 1
    call this%dim_kinds(index)%Init(patch_size_pft_r8, 2)

    ! site x ground
    index = index + 1
    call this%dim_kinds(index)%Init(site_ground_r8, 2)

    ! site x size-class/pft
    index = index + 1
    call this%dim_kinds(index)%Init(site_size_pft_r8, 2)

    ! site x size-class
    index = index + 1
    call this%dim_kinds(index)%Init(site_size_r8, 2)

    ! site x pft
    index = index + 1
    call this%dim_kinds(index)%Init(site_pft_r8, 2)

    ! site x patch-age class
    index = index + 1
    call this%dim_kinds(index)%Init(site_age_r8, 2)

    ! site x fuel size class
    index = index + 1
    call this%dim_kinds(index)%Init(site_fuel_r8, 2)

    ! site x cwd size class
    index = index + 1
    call this%dim_kinds(index)%Init(site_cwdsc_r8, 2)

    ! site x can class
    index = index + 1
    call this%dim_kinds(index)%Init(site_can_r8, 2)

    ! site x cnlf class
    index = index + 1
    call this%dim_kinds(index)%Init(site_cnlf_r8, 2)

    ! site x cnlfpft class
    index = index + 1
    call this%dim_kinds(index)%Init(site_cnlfpft_r8, 2)

    ! site x size-class x age class
    index = index + 1
    call this%dim_kinds(index)%Init(site_scag_r8, 2)

    ! FIXME(bja, 2016-10) assert(index == fates_history_num_dim_kinds)
  end subroutine init_dim_kinds_maps

 ! =======================================================================
 subroutine update_history_cbal(this,nc,nsites,sites)

     use EDtypesMod          , only : ed_site_type
     
     ! Arguments
     class(fates_history_interface_type)             :: this
     integer                 , intent(in)            :: nc   ! clump index
     integer                 , intent(in)            :: nsites
     type(ed_site_type)      , intent(inout), target :: sites(nsites)

     ! Locals
     integer  :: s        ! The local site index
     integer  :: io_si     ! The site index of the IO array
     
     
     associate( hio_nep_si            => this%hvars(ih_nep_si)%r81d, &
                 hio_nbp_si            => this%hvars(ih_nbp_si)%r81d, &
                 hio_fire_c_to_atm_si  => this%hvars(ih_fire_c_to_atm_si)%r81d, &
                 hio_totecosysc_si     => this%hvars(ih_totecosysc_si)%r81d, &
                 hio_cbal_err_fates_si => this%hvars(ih_cbal_err_fates_si)%r81d, &
                 hio_cbal_err_bgc_si   => this%hvars(ih_cbal_err_bgc_si)%r81d, &
                 hio_cbal_err_tot_si   => this%hvars(ih_cbal_err_tot_si)%r81d, &
                 hio_biomass_stock_si  => this%hvars(ih_biomass_stock_si)%r81d, &
                 hio_litter_stock_si   => this%hvars(ih_litter_stock_si)%r81d, &
                 hio_cwd_stock_si      => this%hvars(ih_cwd_stock_si)%r81d )

        ! ---------------------------------------------------------------------------------
        ! Flush arrays to values defined by %flushval (see registry entry in
        ! subroutine define_history_vars()
        ! ---------------------------------------------------------------------------------
        call this%flush_hvars(nc,upfreq_in=3)
        
        
        do s = 1,nsites
         
           io_si  = this%iovar_map(nc)%site_index(s)

           hio_nep_si(io_si) = sites(s)%nep
           hio_nbp_si(io_si) = sites(s)%nbp
           hio_fire_c_to_atm_si(io_si) = sites(s)%fire_c_to_atm
           hio_totecosysc_si(io_si) = sites(s)%totecosysc
           hio_cbal_err_fates_si(io_si) = sites(s)%cbal_err_fates
           hio_cbal_err_bgc_si(io_si) = sites(s)%cbal_err_bgc
           hio_cbal_err_tot_si(io_si) = sites(s)%cbal_err_tot
           hio_biomass_stock_si(io_si) = sites(s)%biomass_stock
           hio_litter_stock_si(io_si) = sites(s)%ed_litter_stock
           hio_cwd_stock_si(io_si) = sites(s)%cwd_stock

        end do

      end associate

   end subroutine update_history_cbal
   

  ! ====================================================================================
  
  subroutine update_history_dyn(this,nc,nsites,sites)
    
    ! ---------------------------------------------------------------------------------
    ! This is the call to update the history IO arrays that are expected to only change
    ! after Ecosystem Dynamics have been processed.
    ! ---------------------------------------------------------------------------------
    
    use EDtypesMod          , only : ed_site_type
    use EDtypesMod          , only : ed_cohort_type
    use EDtypesMod          , only : ed_patch_type
    use EDtypesMod          , only : AREA
    use EDtypesMod          , only : AREA_INV
    use EDtypesMod          , only : nfsc
    use EDtypesMod          , only : ncwd
    use EDtypesMod          , only : ican_upper
    use EDtypesMod          , only : ican_ustory
    use FatesSizeAgeTypeIndicesMod, only : get_sizeage_class_index
    use EDTypesMod        , only : nlevleaf

    ! Arguments
    class(fates_history_interface_type)             :: this
    integer                 , intent(in)            :: nc   ! clump index
    integer                 , intent(in)            :: nsites
    type(ed_site_type)      , intent(inout), target :: sites(nsites)
    
    ! Locals
    integer  :: s        ! The local site index
    integer  :: io_si     ! The site index of the IO array
    integer  :: ipa, ipa2 ! The local "I"ndex of "PA"tches 
    integer  :: io_pa    ! The patch index of the IO array
    integer  :: io_pa1   ! The first patch index in the IO array for each site
    integer  :: io_soipa 
    integer  :: lb1,ub1,lb2,ub2  ! IO array bounds for the calling thread
    integer  :: ivar             ! index of IO variable object vector
    integer  :: ft               ! functional type index
    integer  :: i_scpf,i_pft,i_scls     ! iterators for scpf, pft, and scls dims
    integer  :: i_cwd,i_fuel            ! iterators for cwd and fuel dims
    integer  :: iscag        ! size-class x age index
    integer  :: ican, ileaf, cnlf_indx  ! iterators for leaf and canopy level
    
    real(r8) :: n_density   ! individual of cohort per m2.
    real(r8) :: n_perm2     ! individuals per m2 for the whole column
    real(r8) :: patch_scaling_scalar ! ratio of canopy to patch area for counteracting patch scaling
    real(r8) :: dbh         ! diameter ("at breast height")

    type(ed_patch_type),pointer  :: cpatch
    type(ed_cohort_type),pointer :: ccohort

    real(r8), parameter :: tiny = 1.e-5_r8      ! some small number
    
    associate( hio_npatches_si         => this%hvars(ih_npatches_si)%r81d, &
               hio_ncohorts_si         => this%hvars(ih_ncohorts_si)%r81d, &
               hio_trimming_pa         => this%hvars(ih_trimming_pa)%r81d, &
               hio_area_plant_pa       => this%hvars(ih_area_plant_pa)%r81d, &
               hio_area_treespread_pa  => this%hvars(ih_area_treespread_pa)%r81d, & 
               hio_canopy_spread_si    => this%hvars(ih_canopy_spread_si)%r81d, &
               hio_biomass_si_pft      => this%hvars(ih_biomass_si_pft)%r82d, &
               hio_leafbiomass_si_pft  => this%hvars(ih_leafbiomass_si_pft)%r82d, &
               hio_storebiomass_si_pft => this%hvars(ih_storebiomass_si_pft)%r82d, &
               hio_nindivs_si_pft      => this%hvars(ih_nindivs_si_pft)%r82d, &
               hio_recruitment_si_pft  => this%hvars(ih_recruitment_si_pft)%r82d, &
               hio_mortality_si_pft  => this%hvars(ih_mortality_si_pft)%r82d, &
               hio_nesterov_fire_danger_pa => this%hvars(ih_nesterov_fire_danger_pa)%r81d, &
               hio_spitfire_ros_pa     => this%hvars(ih_spitfire_ROS_pa)%r81d, &
               hio_tfc_ros_pa          => this%hvars(ih_TFC_ROS_pa)%r81d, &
               hio_effect_wspeed_pa    => this%hvars(ih_effect_wspeed_pa)%r81d, &
               hio_fire_intensity_pa   => this%hvars(ih_fire_intensity_pa)%r81d, &
               hio_fire_area_pa        => this%hvars(ih_fire_area_pa)%r81d, &
               hio_scorch_height_pa    => this%hvars(ih_scorch_height_pa)%r81d, &
               hio_fire_fuel_bulkd_pa  => this%hvars(ih_fire_fuel_bulkd_pa)%r81d, &
               hio_fire_fuel_eff_moist_pa => this%hvars(ih_fire_fuel_eff_moist_pa)%r81d, &
               hio_fire_fuel_sav_pa    => this%hvars(ih_fire_fuel_sav_pa)%r81d, &
               hio_fire_fuel_mef_pa    => this%hvars(ih_fire_fuel_mef_pa)%r81d, &
               hio_sum_fuel_pa         => this%hvars(ih_sum_fuel_pa)%r81d,  &
               hio_litter_in_si        => this%hvars(ih_litter_in_si)%r81d, &
               hio_litter_out_pa       => this%hvars(ih_litter_out_pa)%r81d, &
               hio_seed_bank_si        => this%hvars(ih_seed_bank_si)%r81d, &
               hio_seeds_in_pa         => this%hvars(ih_seeds_in_pa)%r81d, &
               hio_seed_decay_pa       => this%hvars(ih_seed_decay_pa)%r81d, &
               hio_seed_germination_pa => this%hvars(ih_seed_germination_pa)%r81d, &
               hio_bstore_pa           => this%hvars(ih_bstore_pa)%r81d, &
               hio_bdead_pa            => this%hvars(ih_bdead_pa)%r81d, &
               hio_balive_pa           => this%hvars(ih_balive_pa)%r81d, &
               hio_bleaf_pa            => this%hvars(ih_bleaf_pa)%r81d, &
               hio_btotal_pa           => this%hvars(ih_btotal_pa)%r81d, &
               hio_canopy_biomass_pa   => this%hvars(ih_canopy_biomass_pa)%r81d, &
               hio_understory_biomass_pa   => this%hvars(ih_understory_biomass_pa)%r81d, &
               hio_gpp_si_scpf         => this%hvars(ih_gpp_si_scpf)%r82d, &
               hio_npp_totl_si_scpf    => this%hvars(ih_npp_totl_si_scpf)%r82d, &
               hio_npp_leaf_si_scpf    => this%hvars(ih_npp_leaf_si_scpf)%r82d, &
               hio_npp_seed_si_scpf    => this%hvars(ih_npp_seed_si_scpf)%r82d, &
               hio_npp_fnrt_si_scpf    => this%hvars(ih_npp_fnrt_si_scpf)%r82d, &
               hio_npp_bgsw_si_scpf    => this%hvars(ih_npp_bgsw_si_scpf)%r82d, &
               hio_npp_bgdw_si_scpf    => this%hvars(ih_npp_bgdw_si_scpf)%r82d, &
               hio_npp_agsw_si_scpf    => this%hvars(ih_npp_agsw_si_scpf)%r82d, &
               hio_npp_agdw_si_scpf    => this%hvars(ih_npp_agdw_si_scpf)%r82d, &
               hio_npp_stor_si_scpf    => this%hvars(ih_npp_stor_si_scpf)%r82d, &
               hio_bstor_canopy_si_scpf      => this%hvars(ih_bstor_canopy_si_scpf)%r82d, &
               hio_bstor_understory_si_scpf  => this%hvars(ih_bstor_understory_si_scpf)%r82d, &
               hio_bleaf_canopy_si_scpf      => this%hvars(ih_bleaf_canopy_si_scpf)%r82d, &
               hio_bleaf_understory_si_scpf  => this%hvars(ih_bleaf_understory_si_scpf)%r82d, &
               hio_mortality_canopy_si_scpf         => this%hvars(ih_mortality_canopy_si_scpf)%r82d, &
               hio_mortality_understory_si_scpf     => this%hvars(ih_mortality_understory_si_scpf)%r82d, &
               hio_nplant_canopy_si_scpf     => this%hvars(ih_nplant_canopy_si_scpf)%r82d, &
               hio_nplant_understory_si_scpf => this%hvars(ih_nplant_understory_si_scpf)%r82d, &
               hio_ddbh_canopy_si_scpf       => this%hvars(ih_ddbh_canopy_si_scpf)%r82d, &
               hio_ddbh_understory_si_scpf   => this%hvars(ih_ddbh_understory_si_scpf)%r82d, &
               hio_ddbh_canopy_si_scls       => this%hvars(ih_ddbh_canopy_si_scls)%r82d, &
               hio_ddbh_understory_si_scls   => this%hvars(ih_ddbh_understory_si_scls)%r82d, &
               hio_gpp_canopy_si_scpf        => this%hvars(ih_gpp_canopy_si_scpf)%r82d, &
               hio_gpp_understory_si_scpf    => this%hvars(ih_gpp_understory_si_scpf)%r82d, &
               hio_ar_canopy_si_scpf         => this%hvars(ih_ar_canopy_si_scpf)%r82d, &
               hio_ar_understory_si_scpf     => this%hvars(ih_ar_understory_si_scpf)%r82d, &
               hio_ddbh_si_scpf        => this%hvars(ih_ddbh_si_scpf)%r82d, &
               hio_ba_si_scpf          => this%hvars(ih_ba_si_scpf)%r82d, &
               hio_nplant_si_scpf      => this%hvars(ih_nplant_si_scpf)%r82d, &
               hio_m1_si_scpf          => this%hvars(ih_m1_si_scpf)%r82d, &
               hio_m2_si_scpf          => this%hvars(ih_m2_si_scpf)%r82d, &
               hio_m3_si_scpf          => this%hvars(ih_m3_si_scpf)%r82d, &
               hio_m4_si_scpf          => this%hvars(ih_m4_si_scpf)%r82d, &
               hio_m5_si_scpf          => this%hvars(ih_m5_si_scpf)%r82d, &
               hio_m6_si_scpf          => this%hvars(ih_m6_si_scpf)%r82d, &

               hio_m7_si_scpf          => this%hvars(ih_m7_si_scpf)%r82d, &                  

               hio_ba_si_scls          => this%hvars(ih_ba_si_scls)%r82d, &
               hio_nplant_canopy_si_scls         => this%hvars(ih_nplant_canopy_si_scls)%r82d, &
               hio_nplant_understory_si_scls     => this%hvars(ih_nplant_understory_si_scls)%r82d, &
               hio_mortality_canopy_si_scls      => this%hvars(ih_mortality_canopy_si_scls)%r82d, &
               hio_mortality_understory_si_scls  => this%hvars(ih_mortality_understory_si_scls)%r82d, &
               hio_demotion_rate_si_scls         => this%hvars(ih_demotion_rate_si_scls)%r82d, &
               hio_demotion_carbonflux_si        => this%hvars(ih_demotion_carbonflux_si)%r81d, &
               hio_promotion_rate_si_scls        => this%hvars(ih_promotion_rate_si_scls)%r82d, &
               hio_trimming_canopy_si_scls         => this%hvars(ih_trimming_canopy_si_scls)%r82d, &
               hio_trimming_understory_si_scls     => this%hvars(ih_trimming_understory_si_scls)%r82d, &
               hio_crown_area_canopy_si_scls         => this%hvars(ih_crown_area_canopy_si_scls)%r82d, &
               hio_crown_area_understory_si_scls     => this%hvars(ih_crown_area_understory_si_scls)%r82d, &
               hio_promotion_carbonflux_si       => this%hvars(ih_promotion_carbonflux_si)%r81d, &
               hio_canopy_mortality_carbonflux_si     => this%hvars(ih_canopy_mortality_carbonflux_si)%r81d, &
               hio_understory_mortality_carbonflux_si => this%hvars(ih_understory_mortality_carbonflux_si)%r81d, &
               hio_leaf_md_canopy_si_scls           => this%hvars(ih_leaf_md_canopy_si_scls)%r82d, &
               hio_root_md_canopy_si_scls           => this%hvars(ih_root_md_canopy_si_scls)%r82d, &
               hio_carbon_balance_canopy_si_scls    => this%hvars(ih_carbon_balance_canopy_si_scls)%r82d, &
               hio_seed_prod_canopy_si_scls         => this%hvars(ih_seed_prod_canopy_si_scls)%r82d, &
               hio_dbalivedt_canopy_si_scls         => this%hvars(ih_dbalivedt_canopy_si_scls)%r82d, &
               hio_dbdeaddt_canopy_si_scls          => this%hvars(ih_dbdeaddt_canopy_si_scls)%r82d, &
               hio_dbstoredt_canopy_si_scls         => this%hvars(ih_dbstoredt_canopy_si_scls)%r82d, &
               hio_storage_flux_canopy_si_scls      => this%hvars(ih_storage_flux_canopy_si_scls)%r82d, &
               hio_npp_leaf_canopy_si_scls          => this%hvars(ih_npp_leaf_canopy_si_scls)%r82d, &
               hio_npp_froot_canopy_si_scls         => this%hvars(ih_npp_froot_canopy_si_scls)%r82d, &
               hio_npp_bsw_canopy_si_scls           => this%hvars(ih_npp_bsw_canopy_si_scls)%r82d, &
               hio_npp_bdead_canopy_si_scls         => this%hvars(ih_npp_bdead_canopy_si_scls)%r82d, &
               hio_npp_bseed_canopy_si_scls         => this%hvars(ih_npp_bseed_canopy_si_scls)%r82d, &
               hio_npp_store_canopy_si_scls         => this%hvars(ih_npp_store_canopy_si_scls)%r82d, &
               hio_leaf_md_understory_si_scls       => this%hvars(ih_leaf_md_understory_si_scls)%r82d, &
               hio_root_md_understory_si_scls       => this%hvars(ih_root_md_understory_si_scls)%r82d, &
               hio_carbon_balance_understory_si_scls=> this%hvars(ih_carbon_balance_understory_si_scls)%r82d, &
               hio_seed_prod_understory_si_scls     => this%hvars(ih_seed_prod_understory_si_scls)%r82d, &
               hio_dbalivedt_understory_si_scls     => this%hvars(ih_dbalivedt_understory_si_scls)%r82d, &
               hio_dbdeaddt_understory_si_scls      => this%hvars(ih_dbdeaddt_understory_si_scls)%r82d, &
               hio_dbstoredt_understory_si_scls     => this%hvars(ih_dbstoredt_understory_si_scls)%r82d, &
               hio_storage_flux_understory_si_scls  => this%hvars(ih_storage_flux_understory_si_scls)%r82d, &
               hio_npp_leaf_understory_si_scls      => this%hvars(ih_npp_leaf_understory_si_scls)%r82d, &
               hio_npp_froot_understory_si_scls     => this%hvars(ih_npp_froot_understory_si_scls)%r82d, &
               hio_npp_bsw_understory_si_scls       => this%hvars(ih_npp_bsw_understory_si_scls)%r82d, &
               hio_npp_bdead_understory_si_scls     => this%hvars(ih_npp_bdead_understory_si_scls)%r82d, &
               hio_npp_bseed_understory_si_scls     => this%hvars(ih_npp_bseed_understory_si_scls)%r82d, &
               hio_npp_store_understory_si_scls     => this%hvars(ih_npp_store_understory_si_scls)%r82d, &
               hio_yesterdaycanopylevel_canopy_si_scls     => this%hvars(ih_yesterdaycanopylevel_canopy_si_scls)%r82d, &
               hio_yesterdaycanopylevel_understory_si_scls => this%hvars(ih_yesterdaycanopylevel_understory_si_scls)%r82d, &
               hio_area_si_age         => this%hvars(ih_area_si_age)%r82d, &
               hio_lai_si_age          => this%hvars(ih_lai_si_age)%r82d, &
               hio_canopy_area_si_age  => this%hvars(ih_canopy_area_si_age)%r82d, &
               hio_ncl_si_age          => this%hvars(ih_ncl_si_age)%r82d, &
               hio_npatches_si_age     => this%hvars(ih_npatches_si_age)%r82d, &
               hio_zstar_si_age        => this%hvars(ih_zstar_si_age)%r82d, &
               hio_biomass_si_age        => this%hvars(ih_biomass_si_age)%r82d, &
               hio_litter_moisture_si_fuel        => this%hvars(ih_litter_moisture_si_fuel)%r82d, &
               hio_cwd_ag_si_cwdsc                  => this%hvars(ih_cwd_ag_si_cwdsc)%r82d, &
               hio_cwd_bg_si_cwdsc                  => this%hvars(ih_cwd_bg_si_cwdsc)%r82d, &
               hio_cwd_ag_in_si_cwdsc               => this%hvars(ih_cwd_ag_in_si_cwdsc)%r82d, &
               hio_cwd_bg_in_si_cwdsc               => this%hvars(ih_cwd_bg_in_si_cwdsc)%r82d, &
               hio_cwd_ag_out_si_cwdsc              => this%hvars(ih_cwd_ag_out_si_cwdsc)%r82d, &
               hio_cwd_bg_out_si_cwdsc              => this%hvars(ih_cwd_bg_out_si_cwdsc)%r82d, &
               hio_crownarea_si_cnlf                => this%hvars(ih_crownarea_si_cnlf)%r82d, &
               hio_crownarea_si_can                 => this%hvars(ih_crownarea_si_can)%r82d, &
               hio_nplant_si_scag                   => this%hvars(ih_nplant_si_scag)%r82d, &
               hio_nplant_canopy_si_scag            => this%hvars(ih_nplant_canopy_si_scag)%r82d, &
               hio_nplant_understory_si_scag        => this%hvars(ih_nplant_understory_si_scag)%r82d, &
               hio_ddbh_canopy_si_scag              => this%hvars(ih_ddbh_canopy_si_scag)%r82d, &
               hio_ddbh_understory_si_scag          => this%hvars(ih_ddbh_understory_si_scag)%r82d, &
               hio_mortality_canopy_si_scag         => this%hvars(ih_mortality_canopy_si_scag)%r82d, &
               hio_mortality_understory_si_scag     => this%hvars(ih_mortality_understory_si_scag)%r82d)

               
      ! ---------------------------------------------------------------------------------
      ! Flush arrays to values defined by %flushval (see registry entry in
      ! subroutine define_history_vars()
      ! ---------------------------------------------------------------------------------
      call this%flush_hvars(nc,upfreq_in=1)


      ! If we don't have dynamics turned on, we just abort these diagnostics
      if (hlm_use_ed_st3.eq.itrue) return

      ! ---------------------------------------------------------------------------------
      ! Loop through the FATES scale hierarchy and fill the history IO arrays
      ! ---------------------------------------------------------------------------------
      
      do s = 1,nsites
         
         io_si  = this%iovar_map(nc)%site_index(s)
         io_pa1 = this%iovar_map(nc)%patch1_index(s)
         io_soipa = io_pa1-1
         
         ! Set trimming on the soil patch to 1.0
         hio_trimming_pa(io_soipa) = 1.0_r8

         ! The seed bank is a site level variable
         hio_seed_bank_si(io_si) = sum(sites(s)%seed_bank) * g_per_kg

         hio_canopy_spread_si(io_si)        = sites(s)%spread
            
         ipa = 0
         cpatch => sites(s)%oldest_patch
         do while(associated(cpatch))
            
            io_pa = io_pa1 + ipa

            ! Increment the number of patches per site
            hio_npatches_si(io_si) = hio_npatches_si(io_si) + 1._r8

            ! Increment the fractional area in each age class bin
            hio_area_si_age(io_si,cpatch%age_class) = hio_area_si_age(io_si,cpatch%age_class) &
                 + cpatch%area * AREA_INV

            ! Increment some patch-age-resolved diagnostics
            hio_lai_si_age(io_si,cpatch%age_class) = hio_lai_si_age(io_si,cpatch%age_class) &
                 + cpatch%lai * cpatch%area
            hio_ncl_si_age(io_si,cpatch%age_class) = hio_ncl_si_age(io_si,cpatch%age_class) &
                 + cpatch%ncl_p * cpatch%area
            hio_npatches_si_age(io_si,cpatch%age_class) = hio_npatches_si_age(io_si,cpatch%age_class) + 1._r8
            if ( ED_val_comp_excln .lt. 0._r8 ) then ! only valid when "strict ppa" enabled
               hio_zstar_si_age(io_si,cpatch%age_class) = hio_zstar_si_age(io_si,cpatch%age_class) &
                    + cpatch%zstar * cpatch%area * AREA_INV
            endif

            ccohort => cpatch%shortest
            do while(associated(ccohort))
               
               ft = ccohort%pft
               
               ! Increment the number of cohorts per site
               hio_ncohorts_si(io_si) = hio_ncohorts_si(io_si) + 1._r8
               
               if ((cpatch%area .gt. 0._r8) .and. (cpatch%total_canopy_area .gt. 0._r8)) then
                  
                  ! for quantities that are at the CLM patch level, because of the way 
                  ! that CLM patches are weighted for radiative purposes this # density needs 
                  ! to be over either ED patch canopy area or ED patch total area, whichever is less
                  n_density = ccohort%n/min(cpatch%area,cpatch%total_canopy_area) 
                  
                  ! for quantities that are natively at column level, calculate plant 
                  ! density using whole area
                  n_perm2   = ccohort%n * AREA_INV
                  
               else
                  n_density = 0.0_r8
                  n_perm2   = 0.0_r8
               endif
               
               if(associated(cpatch%tallest))then
                  hio_trimming_pa(io_pa) = cpatch%tallest%canopy_trim
               else
                  hio_trimming_pa(io_pa) = 0.0_r8
               endif
               
               hio_area_plant_pa(io_pa) = 1._r8
               
               if (min(cpatch%total_canopy_area,cpatch%area)>0.0_r8) then
                  hio_area_treespread_pa(io_pa) = cpatch%total_tree_area  &
                       / min(cpatch%total_canopy_area,cpatch%area)
               else
                  hio_area_treespread_pa(io_pa) = 0.0_r8
               end if
               
               hio_canopy_area_si_age(io_si,cpatch%age_class) = hio_canopy_area_si_age(io_si,cpatch%age_class) &
                    + ccohort%c_area * AREA_INV

               ! Update biomass components
               hio_bleaf_pa(io_pa)  = hio_bleaf_pa(io_pa)  + n_density * ccohort%bl       * g_per_kg
               hio_bstore_pa(io_pa) = hio_bstore_pa(io_pa) + n_density * ccohort%bstore   * g_per_kg
               hio_btotal_pa(io_pa) = hio_btotal_pa(io_pa) + n_density * ccohort%b        * g_per_kg
               hio_bdead_pa(io_pa)  = hio_bdead_pa(io_pa)  + n_density * ccohort%bdead    * g_per_kg
               hio_balive_pa(io_pa) = hio_balive_pa(io_pa) + n_density * ccohort%balive   * g_per_kg
               
               ! Update PFT partitioned biomass components
               hio_leafbiomass_si_pft(io_si,ft) = hio_leafbiomass_si_pft(io_si,ft) + &
                    (ccohort%n * AREA_INV) * ccohort%bl       * g_per_kg
             
               hio_storebiomass_si_pft(io_si,ft) = hio_storebiomass_si_pft(io_si,ft) + &
                    (ccohort%n * AREA_INV) * ccohort%bstore   * g_per_kg
               
               hio_nindivs_si_pft(io_si,ft) = hio_nindivs_si_pft(io_si,ft) + &
                    ccohort%n * AREA_INV

               hio_biomass_si_pft(io_si, ft) = hio_biomass_si_pft(io_si, ft) + &
                    (ccohort%n * AREA_INV) * ccohort%b * g_per_kg

               ! update total biomass per age bin
               hio_biomass_si_age(io_si,cpatch%age_class) = hio_biomass_si_age(io_si,cpatch%age_class) &
                    + ccohort%b * ccohort%n * AREA_INV


               ! Site by Size-Class x PFT (SCPF) 
               ! ------------------------------------------------------------------------

               dbh = ccohort%dbh !-0.5*(1./365.25)*ccohort%ddbhdt

               ! Flux Variables (cohorts must had experienced a day before any of these values
               ! have any meaning, otherwise they are just inialization values
               if( .not.(ccohort%isnew) ) then

                  associate( scpf => ccohort%size_by_pft_class, &
                             scls => ccohort%size_class )

                    hio_gpp_si_scpf(io_si,scpf)      = hio_gpp_si_scpf(io_si,scpf)      + &
                                                       n_perm2*ccohort%gpp_acc_hold  ! [kgC/m2/yr]
                    hio_npp_totl_si_scpf(io_si,scpf) = hio_npp_totl_si_scpf(io_si,scpf) + &
                                                       ccohort%npp_acc_hold *n_perm2
                    hio_npp_leaf_si_scpf(io_si,scpf) = hio_npp_leaf_si_scpf(io_si,scpf) + &
                                                       ccohort%npp_leaf*n_perm2
                    hio_npp_fnrt_si_scpf(io_si,scpf) = hio_npp_fnrt_si_scpf(io_si,scpf) + &
                                                       ccohort%npp_froot*n_perm2
                    hio_npp_bgsw_si_scpf(io_si,scpf) = hio_npp_bgsw_si_scpf(io_si,scpf) + &
                                                       ccohort%npp_bsw*n_perm2*           &
                                                       (1._r8-EDPftvarcon_inst%allom_agb_frac(ccohort%pft))
                    hio_npp_agsw_si_scpf(io_si,scpf) = hio_npp_agsw_si_scpf(io_si,scpf) + &
                                                       ccohort%npp_bsw*n_perm2*           &
                                                       EDPftvarcon_inst%allom_agb_frac(ccohort%pft)
                    hio_npp_bgdw_si_scpf(io_si,scpf) = hio_npp_bgdw_si_scpf(io_si,scpf) + &
                                                       ccohort%npp_bdead*n_perm2*         &
                                                       (1._r8-EDPftvarcon_inst%allom_agb_frac(ccohort%pft))
                    hio_npp_agdw_si_scpf(io_si,scpf) = hio_npp_agdw_si_scpf(io_si,scpf) + &
                                                       ccohort%npp_bdead*n_perm2*         &
                                                       EDPftvarcon_inst%allom_agb_frac(ccohort%pft)
                    hio_npp_seed_si_scpf(io_si,scpf) = hio_npp_seed_si_scpf(io_si,scpf) + &
                                                       ccohort%npp_bseed*n_perm2
                    hio_npp_stor_si_scpf(io_si,scpf) = hio_npp_stor_si_scpf(io_si,scpf) + &
                                                       ccohort%npp_store*n_perm2

                    if( abs(ccohort%npp_acc_hold-(ccohort%npp_leaf+ccohort%npp_froot+ &
                         ccohort%npp_bsw+ccohort%npp_bdead+ &
                         ccohort%npp_bseed+ccohort%npp_store))>1.e-9)  then
                       write(fates_log(),*) 'NPP Partitions are not balancing'
                       write(fates_log(),*) 'Fractional Error: ', &
                            abs(ccohort%npp_acc_hold-(ccohort%npp_leaf+ccohort%npp_froot+ &
                            ccohort%npp_bsw+ccohort%npp_bdead+ &
                            ccohort%npp_bseed+ccohort%npp_store))/ccohort%npp_acc_hold
                       write(fates_log(),*) 'Terms: ',ccohort%npp_acc_hold,ccohort%npp_leaf,ccohort%npp_froot, &
                            ccohort%npp_bsw,ccohort%npp_bdead, &
                            ccohort%npp_bseed,ccohort%npp_store
                       write(fates_log(),*) ' NPP components during FATES-HLM linking does not balance '
                       stop ! we need termination control for FATES!!!
                       ! call endrun(msg=errMsg(__FILE__, __LINE__))
                    end if
                  
                    ! Woody State Variables (basal area and number density and mortality)
                    if (EDPftvarcon_inst%woody(ft) == 1) then
                       
                       hio_m1_si_scpf(io_si,scpf) = hio_m1_si_scpf(io_si,scpf) + ccohort%bmort*ccohort%n
                       hio_m2_si_scpf(io_si,scpf) = hio_m2_si_scpf(io_si,scpf) + ccohort%hmort*ccohort%n
                       hio_m3_si_scpf(io_si,scpf) = hio_m3_si_scpf(io_si,scpf) + ccohort%cmort*ccohort%n
                       hio_m4_si_scpf(io_si,scpf) = hio_m4_si_scpf(io_si,scpf) + ccohort%imort*ccohort%n
                       hio_m5_si_scpf(io_si,scpf) = hio_m5_si_scpf(io_si,scpf) + ccohort%fmort*ccohort%n
                       

                      !Y.X. 
		       hio_m7_si_scpf(io_si,scpf) = hio_m7_si_scpf(io_si,scpf) + &
		       	    (ccohort%lmort_logging+ccohort%lmort_collateral+ccohort%lmort_infra) * ccohort%n


                       ! basal area  [m2/ha]
                       hio_ba_si_scpf(io_si,scpf) = hio_ba_si_scpf(io_si,scpf) + &
                            0.25_r8*3.14159_r8*((dbh/100.0_r8)**2.0_r8)*ccohort%n
                       ! also by size class only
                       hio_ba_si_scls(io_si,scls) = hio_ba_si_scls(io_si,scls) + &
                            0.25_r8*3.14159_r8*((dbh/100.0_r8)**2.0_r8)*ccohort%n
                       
                       ! number density [/ha]
                       hio_nplant_si_scpf(io_si,scpf) = hio_nplant_si_scpf(io_si,scpf) + ccohort%n
                       
                       ! growth increment
                       hio_ddbh_si_scpf(io_si,scpf) = hio_ddbh_si_scpf(io_si,scpf) + &
                            ccohort%ddbhdt*ccohort%n
                    end if

                    ! update size-class x patch-age related quantities

                    iscag = get_sizeage_class_index(ccohort%dbh,cpatch%age)
                    
                    hio_nplant_si_scag(io_si,iscag) = hio_nplant_si_scag(io_si,iscag) + ccohort%n

                    ! update SCPF/SCLS- and canopy/subcanopy- partitioned quantities
                    if (ccohort%canopy_layer .eq. 1) then
                       hio_nplant_canopy_si_scag(io_si,iscag) = hio_nplant_canopy_si_scag(io_si,iscag) + ccohort%n
                       hio_mortality_canopy_si_scag(io_si,iscag) = hio_mortality_canopy_si_scag(io_si,iscag) + &
                            (ccohort%bmort + ccohort%hmort + ccohort%cmort + ccohort%fmort) * ccohort%n
                       hio_ddbh_canopy_si_scag(io_si,iscag) = hio_ddbh_canopy_si_scag(io_si,iscag) + &
                            ccohort%ddbhdt*ccohort%n
                       hio_bstor_canopy_si_scpf(io_si,scpf) = hio_bstor_canopy_si_scpf(io_si,scpf) + &
                            ccohort%bstore * ccohort%n
                       hio_bleaf_canopy_si_scpf(io_si,scpf) = hio_bleaf_canopy_si_scpf(io_si,scpf) + &
                            ccohort%bl * ccohort%n
                       hio_canopy_biomass_pa(io_pa) = hio_canopy_biomass_pa(io_pa) + n_density * ccohort%b * g_per_kg

                       !hio_mortality_canopy_si_scpf(io_si,scpf) = hio_mortality_canopy_si_scpf(io_si,scpf)+ &
                       !    (ccohort%bmort + ccohort%hmort + ccohort%cmort + ccohort%imort + ccohort%fmort) * ccohort%n

                        hio_mortality_canopy_si_scpf(io_si,scpf) = hio_mortality_canopy_si_scpf(io_si,scpf)+ &
                            (ccohort%bmort + ccohort%hmort + ccohort%cmort + ccohort%imort + ccohort%fmort+ &
			    ccohort%lmort_logging + ccohort%lmort_collateral + ccohort%lmort_infra) * ccohort%n

                       hio_nplant_canopy_si_scpf(io_si,scpf) = hio_nplant_canopy_si_scpf(io_si,scpf) + ccohort%n
                       hio_nplant_canopy_si_scls(io_si,scls) = hio_nplant_canopy_si_scls(io_si,scls) + ccohort%n
                       hio_trimming_canopy_si_scls(io_si,scls) = hio_trimming_canopy_si_scls(io_si,scls) + &
                            ccohort%n * ccohort%canopy_trim
                       hio_crown_area_canopy_si_scls(io_si,scls) = hio_crown_area_canopy_si_scls(io_si,scls) + &
                            ccohort%c_area
                       hio_gpp_canopy_si_scpf(io_si,scpf)      = hio_gpp_canopy_si_scpf(io_si,scpf)      + &
                            n_perm2*ccohort%gpp_acc_hold
                       hio_ar_canopy_si_scpf(io_si,scpf)      = hio_ar_canopy_si_scpf(io_si,scpf)      + &
                            n_perm2*ccohort%resp_acc_hold
                       ! growth increment
                       hio_ddbh_canopy_si_scpf(io_si,scpf) = hio_ddbh_canopy_si_scpf(io_si,scpf) + &
                            ccohort%ddbhdt*ccohort%n
                       hio_ddbh_canopy_si_scls(io_si,scls) = hio_ddbh_canopy_si_scls(io_si,scls) + &
                            ccohort%ddbhdt*ccohort%n

                       ! sum of all mortality
                       hio_mortality_canopy_si_scls(io_si,scls) = hio_mortality_canopy_si_scls(io_si,scls) + &
                             (ccohort%bmort + ccohort%hmort + ccohort%cmort  + ccohort%fmort + &
                             ccohort%lmort_logging + ccohort%lmort_collateral + ccohort%lmort_infra) * ccohort%n

                       hio_canopy_mortality_carbonflux_si(io_si) = hio_canopy_mortality_carbonflux_si(io_si) + &
                            (ccohort%bmort + ccohort%hmort + ccohort%cmort + ccohort%fmort) * &
                            ccohort%b * ccohort%n * g_per_kg * days_per_sec * years_per_day * ha_per_m2 + &
                            (ccohort%lmort_logging + ccohort%lmort_collateral + ccohort%lmort_infra)* ccohort%b * &
                            ccohort%n * g_per_kg * ha_per_m2

                       hio_leaf_md_canopy_si_scls(io_si,scls) = hio_leaf_md_canopy_si_scls(io_si,scls) + &
                            ccohort%leaf_md * ccohort%n
                       hio_root_md_canopy_si_scls(io_si,scls) = hio_root_md_canopy_si_scls(io_si,scls) + &
                            ccohort%root_md * ccohort%n
                       hio_carbon_balance_canopy_si_scls(io_si,scls) = hio_carbon_balance_canopy_si_scls(io_si,scls) + &
                            ccohort%carbon_balance * ccohort%n
                       hio_seed_prod_canopy_si_scls(io_si,scls) = hio_seed_prod_canopy_si_scls(io_si,scls) + &
                            ccohort%seed_prod * ccohort%n
                       hio_dbalivedt_canopy_si_scls(io_si,scls) = hio_dbalivedt_canopy_si_scls(io_si,scls) + &
                            ccohort%dbalivedt * ccohort%n
                       hio_dbdeaddt_canopy_si_scls(io_si,scls) = hio_dbdeaddt_canopy_si_scls(io_si,scls) + &
                            ccohort%dbdeaddt * ccohort%n
                       hio_dbstoredt_canopy_si_scls(io_si,scls) = hio_dbstoredt_canopy_si_scls(io_si,scls) + &
                            ccohort%dbstoredt * ccohort%n
                       hio_storage_flux_canopy_si_scls(io_si,scls) = hio_storage_flux_canopy_si_scls(io_si,scls) + &
                            ccohort%storage_flux * ccohort%n
                       hio_npp_leaf_canopy_si_scls(io_si,scls) = hio_npp_leaf_canopy_si_scls(io_si,scls) + &
                            ccohort%npp_leaf * ccohort%n
                       hio_npp_froot_canopy_si_scls(io_si,scls) = hio_npp_froot_canopy_si_scls(io_si,scls) + &
                            ccohort%npp_froot * ccohort%n
                       hio_npp_bsw_canopy_si_scls(io_si,scls) = hio_npp_bsw_canopy_si_scls(io_si,scls) + &
                            ccohort%npp_bsw * ccohort%n
                       hio_npp_bdead_canopy_si_scls(io_si,scls) = hio_npp_bdead_canopy_si_scls(io_si,scls) + &
                            ccohort%npp_bdead * ccohort%n
                       hio_npp_bseed_canopy_si_scls(io_si,scls) = hio_npp_bseed_canopy_si_scls(io_si,scls) + &
                            ccohort%npp_bseed * ccohort%n
                       hio_npp_store_canopy_si_scls(io_si,scls) = hio_npp_store_canopy_si_scls(io_si,scls) + &
                            ccohort%npp_store * ccohort%n
                       hio_yesterdaycanopylevel_canopy_si_scls(io_si,scls) = &
                            hio_yesterdaycanopylevel_canopy_si_scls(io_si,scls) + &
                            ccohort%canopy_layer_yesterday * ccohort%n
                    else
                       hio_nplant_understory_si_scag(io_si,iscag) = hio_nplant_understory_si_scag(io_si,iscag) + ccohort%n
                       hio_mortality_understory_si_scag(io_si,iscag) = hio_mortality_understory_si_scag(io_si,iscag) + &
                            (ccohort%bmort + ccohort%hmort + ccohort%cmort + ccohort%fmort) * ccohort%n
                       hio_ddbh_understory_si_scag(io_si,iscag) = hio_ddbh_understory_si_scag(io_si,iscag) + &
                            ccohort%ddbhdt*ccohort%n
                       hio_bstor_understory_si_scpf(io_si,scpf) = hio_bstor_understory_si_scpf(io_si,scpf) + &
                            ccohort%bstore * ccohort%n
                       hio_bleaf_understory_si_scpf(io_si,scpf) = hio_bleaf_understory_si_scpf(io_si,scpf) + &
                            ccohort%bl * ccohort%n
                       hio_understory_biomass_pa(io_pa) = hio_understory_biomass_pa(io_pa) + n_density * ccohort%b * g_per_kg
                       !hio_mortality_understory_si_scpf(io_si,scpf) = hio_mortality_understory_si_scpf(io_si,scpf)+ &
                        !    (ccohort%bmort + ccohort%hmort + ccohort%cmort + ccohort%imort + ccohort%fmort) * ccohort%n

                       hio_mortality_understory_si_scpf(io_si,scpf) = hio_mortality_understory_si_scpf(io_si,scpf)+ &
                            (ccohort%bmort + ccohort%hmort + ccohort%cmort + ccohort%imort + ccohort%fmort + &
			    ccohort%lmort_logging + ccohort%lmort_collateral + ccohort%lmort_infra) * ccohort%n   

                       hio_nplant_understory_si_scpf(io_si,scpf) = hio_nplant_understory_si_scpf(io_si,scpf) + ccohort%n
                       hio_nplant_understory_si_scls(io_si,scls) = hio_nplant_understory_si_scls(io_si,scls) + ccohort%n
                       hio_trimming_understory_si_scls(io_si,scls) = hio_trimming_understory_si_scls(io_si,scls) + &
                            ccohort%n * ccohort%canopy_trim
                       hio_crown_area_understory_si_scls(io_si,scls) = hio_crown_area_understory_si_scls(io_si,scls) + &
                            ccohort%c_area
                       hio_gpp_understory_si_scpf(io_si,scpf)      = hio_gpp_understory_si_scpf(io_si,scpf)      + &
                            n_perm2*ccohort%gpp_acc_hold
                       hio_ar_understory_si_scpf(io_si,scpf)      = hio_ar_understory_si_scpf(io_si,scpf)      + &
                            n_perm2*ccohort%resp_acc_hold

                       ! growth increment
                       hio_ddbh_understory_si_scpf(io_si,scpf) = hio_ddbh_understory_si_scpf(io_si,scpf) + &
                            ccohort%ddbhdt*ccohort%n
                       hio_ddbh_understory_si_scls(io_si,scls) = hio_ddbh_understory_si_scls(io_si,scls) + &
                            ccohort%ddbhdt*ccohort%n

                       ! sum of all mortality
                       hio_mortality_understory_si_scls(io_si,scls) = hio_mortality_understory_si_scls(io_si,scls) + &
                             (ccohort%bmort + ccohort%hmort + ccohort%cmort + ccohort%imort + ccohort%fmort +&
                             ccohort%lmort_logging + ccohort%lmort_collateral + ccohort%lmort_infra) * ccohort%n
                       
                       hio_understory_mortality_carbonflux_si(io_si) = hio_understory_mortality_carbonflux_si(io_si) + &
                             (ccohort%bmort + ccohort%hmort + ccohort%cmort + ccohort%imort + ccohort%fmort) * &
                             ccohort%b * ccohort%n * g_per_kg * days_per_sec * years_per_day * ha_per_m2 + &
                             (ccohort%lmort_logging + ccohort%lmort_collateral + ccohort%lmort_infra) * ccohort%b * &
                             ccohort%n * g_per_kg * ha_per_m2

                       !
                       hio_leaf_md_understory_si_scls(io_si,scls) = hio_leaf_md_understory_si_scls(io_si,scls) + &
                            ccohort%leaf_md * ccohort%n
                       hio_root_md_understory_si_scls(io_si,scls) = hio_root_md_understory_si_scls(io_si,scls) + &
                            ccohort%root_md * ccohort%n
                       hio_carbon_balance_understory_si_scls(io_si,scls) = hio_carbon_balance_understory_si_scls(io_si,scls) + &
                            ccohort%carbon_balance * ccohort%n
                       hio_seed_prod_understory_si_scls(io_si,scls) = hio_seed_prod_understory_si_scls(io_si,scls) + &
                            ccohort%seed_prod * ccohort%n
                       hio_dbalivedt_understory_si_scls(io_si,scls) = hio_dbalivedt_understory_si_scls(io_si,scls) + &
                            ccohort%dbalivedt * ccohort%n
                       hio_dbdeaddt_understory_si_scls(io_si,scls) = hio_dbdeaddt_understory_si_scls(io_si,scls) + &
                            ccohort%dbdeaddt * ccohort%n
                       hio_dbstoredt_understory_si_scls(io_si,scls) = hio_dbstoredt_understory_si_scls(io_si,scls) + &
                            ccohort%dbstoredt * ccohort%n
                       hio_storage_flux_understory_si_scls(io_si,scls) = hio_storage_flux_understory_si_scls(io_si,scls) + &
                            ccohort%storage_flux * ccohort%n
                       hio_npp_leaf_understory_si_scls(io_si,scls) = hio_npp_leaf_understory_si_scls(io_si,scls) + &
                            ccohort%npp_leaf * ccohort%n
                       hio_npp_froot_understory_si_scls(io_si,scls) = hio_npp_froot_understory_si_scls(io_si,scls) + &
                            ccohort%npp_froot * ccohort%n
                       hio_npp_bsw_understory_si_scls(io_si,scls) = hio_npp_bsw_understory_si_scls(io_si,scls) + &
                            ccohort%npp_bsw * ccohort%n
                       hio_npp_bdead_understory_si_scls(io_si,scls) = hio_npp_bdead_understory_si_scls(io_si,scls) + &
                            ccohort%npp_bdead * ccohort%n
                       hio_npp_bseed_understory_si_scls(io_si,scls) = hio_npp_bseed_understory_si_scls(io_si,scls) + &
                            ccohort%npp_bseed * ccohort%n
                       hio_npp_store_understory_si_scls(io_si,scls) = hio_npp_store_understory_si_scls(io_si,scls) + &
                            ccohort%npp_store * ccohort%n
                       hio_yesterdaycanopylevel_understory_si_scls(io_si,scls) = &
                            hio_yesterdaycanopylevel_understory_si_scls(io_si,scls) + &
                            ccohort%canopy_layer_yesterday * ccohort%n
                    endif
                    !
                    ! consider imort as understory mortality even if it happens in 
                    ! cohorts that may have been promoted as part of the patch creation...
                    hio_mortality_understory_si_scpf(io_si,scpf) = hio_mortality_understory_si_scpf(io_si,scpf)+ &
                            (ccohort%imort) * ccohort%n
                    hio_mortality_understory_si_scls(io_si,scls) = hio_mortality_understory_si_scls(io_si,scls) + &
                         (ccohort%imort) * ccohort%n
                    hio_understory_mortality_carbonflux_si(io_si) = hio_understory_mortality_carbonflux_si(io_si) + &
                         (ccohort%imort) * &
                         ccohort%b * ccohort%n * g_per_kg * days_per_sec * years_per_day * ha_per_m2
                    hio_mortality_understory_si_scag(io_si,iscag) = hio_mortality_understory_si_scag(io_si,iscag) + &
                         (ccohort%imort) * ccohort%n
                    !
                    ccohort%canopy_layer_yesterday = real(ccohort%canopy_layer, r8)
                    
                  end associate
               end if

               ! resolve some canopy area profiles, both total and of occupied leaves
               ican = ccohort%canopy_layer
               !
               hio_crownarea_si_can(io_si, ican) = hio_crownarea_si_can(io_si, ican) + ccohort%c_area / AREA
               !
               do ileaf=1,ccohort%nv
                  cnlf_indx = ileaf + (ican-1) * nlevleaf
                  hio_crownarea_si_cnlf(io_si, cnlf_indx) = hio_crownarea_si_cnlf(io_si, cnlf_indx) + &
                       ccohort%c_area / AREA
               end do
               
               ccohort => ccohort%taller
            enddo ! cohort loop
            
            ! Patch specific variables that are already calculated
            ! These things are all duplicated. Should they all be converted to LL or array structures RF? 
            ! define scalar to counteract the patch albedo scaling logic for conserved quantities
            
            if (cpatch%area .gt. 0._r8 .and. cpatch%total_canopy_area .gt.0 ) then
               patch_scaling_scalar  = min(1._r8, cpatch%area / cpatch%total_canopy_area)
            else
               patch_scaling_scalar = 0._r8
            endif
            
            ! Update Fire Variables
            hio_nesterov_fire_danger_pa(io_pa) = sites(s)%acc_NI
            hio_spitfire_ros_pa(io_pa)         = cpatch%ROS_front 
            hio_effect_wspeed_pa(io_pa)        = cpatch%effect_wspeed
            hio_tfc_ros_pa(io_pa)              = cpatch%TFC_ROS
            hio_fire_intensity_pa(io_pa)       = cpatch%FI
            hio_fire_area_pa(io_pa)            = cpatch%frac_burnt
            hio_scorch_height_pa(io_pa)        = cpatch%SH
            hio_fire_fuel_bulkd_pa(io_pa)      = cpatch%fuel_bulkd
            hio_fire_fuel_eff_moist_pa(io_pa)  = cpatch%fuel_eff_moist
            hio_fire_fuel_sav_pa(io_pa)        = cpatch%fuel_sav
            hio_fire_fuel_mef_pa(io_pa)        = cpatch%fuel_mef
            hio_sum_fuel_pa(io_pa)             = cpatch%sum_fuel * g_per_kg * patch_scaling_scalar
            
            do i_fuel = 1,nfsc
               hio_litter_moisture_si_fuel(io_si, i_fuel) = hio_litter_moisture_si_fuel(io_si, i_fuel) + &
                    cpatch%litter_moisture(i_fuel) * cpatch%area * AREA_INV
            end do
            ! Update Litter Flux Variables

            ! put litter_in flux onto site level variable so as to be able to append site-level distubance-related input flux after patch loop
            hio_litter_in_si(io_si) = hio_litter_in_si(io_si) + &
                 (sum(cpatch%CWD_AG_in) +sum(cpatch%leaf_litter_in) + sum(cpatch%root_litter_in)) &
                 * g_per_kg * cpatch%area * AREA_INV * years_per_day * days_per_sec
            ! keep litter_out at patch level
            hio_litter_out_pa(io_pa)           = (sum(cpatch%CWD_AG_out)+sum(cpatch%leaf_litter_out) &
                 + sum(cpatch%root_litter_out)) &
                 * g_per_kg * patch_scaling_scalar * years_per_day * days_per_sec
            
            hio_seeds_in_pa(io_pa)             = sum(cpatch%seeds_in) * &
                 g_per_kg * patch_scaling_scalar * years_per_day * days_per_sec
            hio_seed_decay_pa(io_pa)           = sum(cpatch%seed_decay) * &
                 g_per_kg * patch_scaling_scalar * years_per_day * days_per_sec
            hio_seed_germination_pa(io_pa)     = sum(cpatch%seed_germination) * &
                 g_per_kg * patch_scaling_scalar * years_per_day * days_per_sec 

            
            do i_cwd = 1, ncwd
               hio_cwd_ag_si_cwdsc(io_si, i_cwd) = hio_cwd_ag_si_cwdsc(io_si, i_cwd) + &
                    cpatch%CWD_AG(i_cwd)*cpatch%area * AREA_INV * g_per_kg
               hio_cwd_bg_si_cwdsc(io_si, i_cwd) = hio_cwd_bg_si_cwdsc(io_si, i_cwd) + &
                    cpatch%CWD_BG(i_cwd)*cpatch%area * AREA_INV * g_per_kg
               hio_cwd_ag_in_si_cwdsc(io_si, i_cwd) = hio_cwd_ag_in_si_cwdsc(io_si, i_cwd) + &
                    cpatch%CWD_AG_IN(i_cwd)*cpatch%area * AREA_INV * g_per_kg
               hio_cwd_bg_in_si_cwdsc(io_si, i_cwd) = hio_cwd_bg_in_si_cwdsc(io_si, i_cwd) + &
                    cpatch%CWD_BG_IN(i_cwd)*cpatch%area * AREA_INV * g_per_kg
               hio_cwd_ag_out_si_cwdsc(io_si, i_cwd) = hio_cwd_ag_out_si_cwdsc(io_si, i_cwd) + &
                    cpatch%CWD_AG_OUT(i_cwd)*cpatch%area * AREA_INV * g_per_kg
               hio_cwd_bg_out_si_cwdsc(io_si, i_cwd) = hio_cwd_bg_out_si_cwdsc(io_si, i_cwd) + &
                    cpatch%CWD_BG_OUT(i_cwd)*cpatch%area * AREA_INV * g_per_kg
            end do

            ipa = ipa + 1
            cpatch => cpatch%younger
         end do !patch loop

         ! divide so-far-just-summed but to-be-averaged patch-age-class variables by patch-age-class area to get mean values
         do ipa2 = 1, nlevage
            if (hio_area_si_age(io_si, ipa2) .gt. tiny) then
               hio_lai_si_age(io_si, ipa2) = hio_lai_si_age(io_si, ipa2) / (hio_area_si_age(io_si, ipa2)*AREA)
               hio_ncl_si_age(io_si, ipa2) = hio_ncl_si_age(io_si, ipa2) / (hio_area_si_age(io_si, ipa2)*AREA)
            else
               hio_lai_si_age(io_si, ipa2) = 0._r8
               hio_ncl_si_age(io_si, ipa2) = 0._r8
            endif
         end do

         ! pass the cohort termination mortality as a flux to the history, and then reset the termination mortality buffer
         ! note there are various ways of reporting the total mortality, so pass to these as well
         do i_pft = 1, numpft
            do i_scls = 1,nlevsclass
               i_scpf = (i_pft-1)*nlevsclass + i_scls
               hio_m6_si_scpf(io_si,i_scpf) = (sites(s)%terminated_nindivs(i_scls,i_pft,1) + &
                    sites(s)%terminated_nindivs(i_scls,i_pft,2)) * days_per_year
               hio_mortality_canopy_si_scls(io_si,i_scls) = hio_mortality_canopy_si_scls(io_si,i_scls) + &
                    sites(s)%terminated_nindivs(i_scls,i_pft,1) * days_per_year
               hio_mortality_understory_si_scls(io_si,i_scls) = hio_mortality_understory_si_scls(io_si,i_scls) + &
                    sites(s)%terminated_nindivs(i_scls,i_pft,2) * days_per_year
               hio_mortality_canopy_si_scpf(io_si,i_scpf) = hio_mortality_canopy_si_scpf(io_si,i_scpf) + &
                    sites(s)%terminated_nindivs(i_scls,i_pft,1) * days_per_year
               hio_mortality_understory_si_scpf(io_si,i_scpf) = hio_mortality_understory_si_scpf(io_si,i_scpf) + &
                    sites(s)%terminated_nindivs(i_scls,i_pft,2) * days_per_year
            end do
         end do
         sites(s)%terminated_nindivs(:,:,:) = 0._r8

         ! pass the recruitment rate as a flux to the history, and then reset the recruitment buffer
         do i_pft = 1, numpft
            hio_recruitment_si_pft(io_si,i_pft) = sites(s)%recruitment_rate(i_pft) * days_per_year
         end do
         sites(s)%recruitment_rate(:) = 0._r8

         ! summarize all of the mortality fluxes by PFT
         do i_pft = 1, numpft
            do i_scls = 1,nlevsclass
               i_scpf = (i_pft-1)*nlevsclass + i_scls
               !hio_mortality_si_pft(io_si,i_pft) = hio_mortality_si_pft(io_si,i_pft) + &
               !     hio_m1_si_scpf(io_si,i_scpf) + &
               !     hio_m2_si_scpf(io_si,i_scpf) + &
               !     hio_m3_si_scpf(io_si,i_scpf) + &
               !     hio_m4_si_scpf(io_si,i_scpf) + &
               !     hio_m5_si_scpf(io_si,i_scpf) + &
               !     hio_m6_si_scpf(io_si,i_scpf) 
         
               hio_mortality_si_pft(io_si,i_pft) = hio_mortality_si_pft(io_si,i_pft) + &
                    hio_m1_si_scpf(io_si,i_scpf) + &
                    hio_m2_si_scpf(io_si,i_scpf) + &
                    hio_m3_si_scpf(io_si,i_scpf) + &
                    hio_m4_si_scpf(io_si,i_scpf) + &
                    hio_m5_si_scpf(io_si,i_scpf) + &
                    hio_m6_si_scpf(io_si,i_scpf) + &
		    hio_m7_si_scpf(io_si,i_scpf)



            end do
         end do

         ! pass demotion rates and associated carbon fluxes to history
         do i_scls = 1,nlevsclass
            hio_demotion_rate_si_scls(io_si,i_scls) = sites(s)%demotion_rate(i_scls) * days_per_year
            hio_promotion_rate_si_scls(io_si,i_scls) = sites(s)%promotion_rate(i_scls) * days_per_year
         end do
         !
         ! convert kg C / ha / day to gc / m2 / sec
         hio_demotion_carbonflux_si(io_si) = sites(s)%demotion_carbonflux * g_per_kg * ha_per_m2 * days_per_sec
         hio_promotion_carbonflux_si(io_si) = sites(s)%promotion_carbonflux * g_per_kg * ha_per_m2 * days_per_sec
         !
         ! mortality-associated carbon fluxes
         
         hio_canopy_mortality_carbonflux_si(io_si) = hio_canopy_mortality_carbonflux_si(io_si) + &
              sites(s)%termination_carbonflux(ican_upper) * g_per_kg * days_per_sec * ha_per_m2
         hio_understory_mortality_carbonflux_si(io_si) = hio_understory_mortality_carbonflux_si(io_si) + &
              sites(s)%termination_carbonflux(ican_ustory) * g_per_kg * days_per_sec * ha_per_m2
         ! and zero the site-level termination carbon flux variable
         sites(s)%termination_carbonflux(:) = 0._r8
         !
         ! add the site-level disturbance-associated cwd and litter input fluxes to thir respective flux fields
         do i_cwd = 1, ncwd
            hio_cwd_ag_in_si_cwdsc(io_si, i_cwd) = hio_cwd_ag_in_si_cwdsc(io_si, i_cwd) + &
                 sites(s)%CWD_AG_diagnostic_input_carbonflux(i_cwd) * g_per_kg
            hio_cwd_bg_in_si_cwdsc(io_si, i_cwd) = hio_cwd_bg_in_si_cwdsc(io_si, i_cwd) + &
                 sites(s)%CWD_BG_diagnostic_input_carbonflux(i_cwd) * g_per_kg
         end do
         hio_litter_in_si(io_si) = hio_litter_in_si(io_si) + &
              (sum(sites(s)%leaf_litter_diagnostic_input_carbonflux) + &
              sum(sites(s)%root_litter_diagnostic_input_carbonflux)) * g_per_kg * days_per_sec * years_per_day
         ! and reset the disturbance-related field buffers
         sites(s)%CWD_AG_diagnostic_input_carbonflux(:) = 0._r8
         sites(s)%CWD_BG_diagnostic_input_carbonflux(:) = 0._r8
         sites(s)%leaf_litter_diagnostic_input_carbonflux(:) = 0._r8
         sites(s)%root_litter_diagnostic_input_carbonflux(:) = 0._r8

      enddo ! site loop
      
    end associate

    return
  end subroutine update_history_dyn
 
 ! ======================================================================================

 subroutine update_history_prod(this,nc,nsites,sites,dt_tstep)

    ! ---------------------------------------------------------------------------------
    ! This is the call to update the history IO arrays that are expected to only change
    ! after rapid timescale productivity calculations (gpp and respiration).
    ! ---------------------------------------------------------------------------------
    
    use EDtypesMod          , only : ed_site_type,   &
                                     ed_cohort_type, &
                                     ed_patch_type,  &
                                     AREA,           &
                                     AREA_INV

    use EDTypesMod          , only : nclmax, nlevleaf
    !
    ! Arguments
    class(fates_history_interface_type)                 :: this
    integer                 , intent(in)            :: nc   ! clump index
    integer                 , intent(in)            :: nsites
    type(ed_site_type)      , intent(inout), target :: sites(nsites)
    real(r8)                , intent(in)            :: dt_tstep
    
    ! Locals
    integer  :: s        ! The local site index
    integer  :: io_si     ! The site index of the IO array
    integer  :: ipa      ! The local "I"ndex of "PA"tches 
    integer  :: io_pa    ! The patch index of the IO array
    integer  :: io_pa1   ! The first patch index in the IO array for each site
    integer  :: io_soipa 
    integer  :: lb1,ub1,lb2,ub2  ! IO array bounds for the calling thread
    integer  :: ivar             ! index of IO variable object vector
    integer  :: ft               ! functional type index
    real(r8) :: n_density   ! individual of cohort per m2.
    real(r8) :: n_perm2     ! individuals per m2 for the whole column
    real(r8) :: patch_area_by_age(nlevage) ! patch area in each bin for normalizing purposes
    real(r8), parameter :: tiny = 1.e-5_r8      ! some small number
    integer  :: ipa2     ! patch incrementer
    integer :: cnlfpft_indx, cnlf_indx, ipft, ican, ileaf ! more iterators and indices
    type(ed_patch_type),pointer  :: cpatch
    type(ed_cohort_type),pointer :: ccohort
    real(r8) :: per_dt_tstep          ! Time step in frequency units (/s)

    associate( hio_gpp_pa         => this%hvars(ih_gpp_pa)%r81d, &
               hio_npp_pa         => this%hvars(ih_npp_pa)%r81d, &
               hio_aresp_pa       => this%hvars(ih_aresp_pa)%r81d, &
               hio_maint_resp_pa  => this%hvars(ih_maint_resp_pa)%r81d, &
               hio_growth_resp_pa => this%hvars(ih_growth_resp_pa)%r81d, &
               hio_npp_si         => this%hvars(ih_npp_si)%r81d, &
               hio_ar_si_scpf     => this%hvars(ih_ar_si_scpf)%r82d, &
               hio_ar_grow_si_scpf   => this%hvars(ih_ar_grow_si_scpf)%r82d, &
               hio_ar_maint_si_scpf  => this%hvars(ih_ar_maint_si_scpf)%r82d, &
               hio_ar_agsapm_si_scpf => this%hvars(ih_ar_agsapm_si_scpf)%r82d, &
               hio_ar_darkm_si_scpf  => this%hvars(ih_ar_darkm_si_scpf)%r82d, &
               hio_ar_crootm_si_scpf => this%hvars(ih_ar_crootm_si_scpf)%r82d, &
               hio_ar_frootm_si_scpf => this%hvars(ih_ar_frootm_si_scpf)%r82d, &
               hio_gpp_canopy_pa     => this%hvars(ih_gpp_canopy_pa)%r81d, &
               hio_ar_canopy_pa      => this%hvars(ih_ar_canopy_pa)%r81d, &
               hio_gpp_understory_pa => this%hvars(ih_gpp_understory_pa)%r81d, &
               hio_ar_understory_pa  => this%hvars(ih_ar_understory_pa)%r81d, &
               hio_rdark_canopy_si_scls             => this%hvars(ih_rdark_canopy_si_scls)%r82d, &
               hio_livestem_mr_canopy_si_scls       => this%hvars(ih_livestem_mr_canopy_si_scls)%r82d, &
               hio_livecroot_mr_canopy_si_scls      => this%hvars(ih_livecroot_mr_canopy_si_scls)%r82d, &
               hio_froot_mr_canopy_si_scls          => this%hvars(ih_froot_mr_canopy_si_scls)%r82d, &
               hio_resp_g_canopy_si_scls            => this%hvars(ih_resp_g_canopy_si_scls)%r82d, &
               hio_resp_m_canopy_si_scls            => this%hvars(ih_resp_m_canopy_si_scls)%r82d, &
               hio_rdark_understory_si_scls         => this%hvars(ih_rdark_understory_si_scls)%r82d, &
               hio_livestem_mr_understory_si_scls   => this%hvars(ih_livestem_mr_understory_si_scls)%r82d, &
               hio_livecroot_mr_understory_si_scls  => this%hvars(ih_livecroot_mr_understory_si_scls)%r82d, &
               hio_froot_mr_understory_si_scls      => this%hvars(ih_froot_mr_understory_si_scls)%r82d, &
               hio_resp_g_understory_si_scls        => this%hvars(ih_resp_g_understory_si_scls)%r82d, &
               hio_resp_m_understory_si_scls        => this%hvars(ih_resp_m_understory_si_scls)%r82d, &
               hio_gpp_si_age         => this%hvars(ih_gpp_si_age)%r82d, &
               hio_npp_si_age         => this%hvars(ih_npp_si_age)%r82d, &
               hio_parsun_z_si_cnlf     => this%hvars(ih_parsun_z_si_cnlf)%r82d, &
               hio_parsha_z_si_cnlf     => this%hvars(ih_parsha_z_si_cnlf)%r82d, &
               hio_ts_net_uptake_si_cnlf     => this%hvars(ih_ts_net_uptake_si_cnlf)%r82d, &
               hio_year_net_uptake_si_cnlf     => this%hvars(ih_year_net_uptake_si_cnlf)%r82d, &
               hio_parsun_z_si_cnlfpft  => this%hvars(ih_parsun_z_si_cnlfpft)%r82d, &
               hio_parsha_z_si_cnlfpft  => this%hvars(ih_parsha_z_si_cnlfpft)%r82d, &
               hio_laisun_z_si_cnlf     => this%hvars(ih_laisun_z_si_cnlf)%r82d, &
               hio_laisha_z_si_cnlf     => this%hvars(ih_laisha_z_si_cnlf)%r82d, &
               hio_laisun_z_si_cnlfpft  => this%hvars(ih_laisun_z_si_cnlfpft)%r82d, &
               hio_laisha_z_si_cnlfpft  => this%hvars(ih_laisha_z_si_cnlfpft)%r82d, &
               hio_laisun_top_si_can     => this%hvars(ih_laisun_top_si_can)%r82d, &
               hio_laisha_top_si_can     => this%hvars(ih_laisha_top_si_can)%r82d, &
               hio_fabd_sun_si_cnlfpft  => this%hvars(ih_fabd_sun_si_cnlfpft)%r82d, &
               hio_fabd_sha_si_cnlfpft  => this%hvars(ih_fabd_sha_si_cnlfpft)%r82d, &
               hio_fabi_sun_si_cnlfpft  => this%hvars(ih_fabi_sun_si_cnlfpft)%r82d, &
               hio_fabi_sha_si_cnlfpft  => this%hvars(ih_fabi_sha_si_cnlfpft)%r82d, &
               hio_fabd_sun_si_cnlf  => this%hvars(ih_fabd_sun_si_cnlf)%r82d, &
               hio_fabd_sha_si_cnlf  => this%hvars(ih_fabd_sha_si_cnlf)%r82d, &
               hio_fabi_sun_si_cnlf  => this%hvars(ih_fabi_sun_si_cnlf)%r82d, &
               hio_fabi_sha_si_cnlf  => this%hvars(ih_fabi_sha_si_cnlf)%r82d, &
               hio_fabd_sun_top_si_can  => this%hvars(ih_fabd_sun_top_si_can)%r82d, &
               hio_fabd_sha_top_si_can  => this%hvars(ih_fabd_sha_top_si_can)%r82d, &
               hio_fabi_sun_top_si_can  => this%hvars(ih_fabi_sun_top_si_can)%r82d, &
               hio_fabi_sha_top_si_can  => this%hvars(ih_fabi_sha_top_si_can)%r82d, &
               hio_parsun_top_si_can     => this%hvars(ih_parsun_top_si_can)%r82d, &
               hio_parsha_top_si_can     => this%hvars(ih_parsha_top_si_can)%r82d &
 )


      ! Flush the relevant history variables 
      call this%flush_hvars(nc,upfreq_in=2)

      per_dt_tstep = 1.0_r8/dt_tstep

      do s = 1,nsites
         
         io_si  = this%iovar_map(nc)%site_index(s)
         io_pa1 = this%iovar_map(nc)%patch1_index(s)
         io_soipa = io_pa1-1
         
         ipa = 0
         cpatch => sites(s)%oldest_patch

         patch_area_by_age(:) = 0._r8

         do while(associated(cpatch))
            
            io_pa = io_pa1 + ipa

            patch_area_by_age(cpatch%age_class) = patch_area_by_age(cpatch%age_class) + cpatch%area

            ccohort => cpatch%shortest
            do while(associated(ccohort))
               
               ! TODO: we need a standardized logical function on this (used lots, RGK)
               if ((cpatch%area .gt. 0._r8) .and. (cpatch%total_canopy_area .gt. 0._r8)) then
                  n_density = ccohort%n/min(cpatch%area,cpatch%total_canopy_area) 
                  n_perm2   = ccohort%n * AREA_INV
               else
                  n_density = 0.0_r8
                  n_perm2   = 0.0_r8
               endif
               
               if ( .not. ccohort%isnew ) then

                  ! Calculate index for the scpf class
                  associate( scpf => ccohort%size_by_pft_class, &
                             scls => ccohort%size_class )

                  ! scale up cohort fluxes to their patches
                  hio_npp_pa(io_pa) = hio_npp_pa(io_pa) + &
                        ccohort%npp_tstep * g_per_kg * n_density * per_dt_tstep
                  hio_gpp_pa(io_pa) = hio_gpp_pa(io_pa) + &
                        ccohort%gpp_tstep * g_per_kg * n_density * per_dt_tstep
                  hio_aresp_pa(io_pa) = hio_aresp_pa(io_pa) + &
                        ccohort%resp_tstep * g_per_kg * n_density * per_dt_tstep
                  hio_growth_resp_pa(io_pa) = hio_growth_resp_pa(io_pa) + &
                        ccohort%resp_g * g_per_kg * n_density * per_dt_tstep
                  hio_maint_resp_pa(io_pa) = hio_maint_resp_pa(io_pa) + &
                        ccohort%resp_m * g_per_kg * n_density * per_dt_tstep
                  
                  ! map ed cohort-level npp fluxes to clm column fluxes
                  hio_npp_si(io_si) = hio_npp_si(io_si) + ccohort%npp_tstep * n_perm2 * g_per_kg * per_dt_tstep


                  ! Total AR (kgC/m2/yr) = (kgC/plant/step) / (s/step) * (plant/m2) * (s/yr)
                  hio_ar_si_scpf(io_si,scpf)    =   hio_ar_si_scpf(io_si,scpf) + &
                        (ccohort%resp_tstep/dt_tstep) * n_perm2 * sec_per_day * days_per_year

                  ! Growth AR (kgC/m2/yr)
                  hio_ar_grow_si_scpf(io_si,scpf) = hio_ar_grow_si_scpf(io_si,scpf) + &
                        (ccohort%resp_g/dt_tstep) * n_perm2 * sec_per_day * days_per_year

                  ! Maint AR (kgC/m2/yr)
                  hio_ar_maint_si_scpf(io_si,scpf) = hio_ar_maint_si_scpf(io_si,scpf) + &
                        (ccohort%resp_m/dt_tstep) * n_perm2 * sec_per_day * days_per_year
                  
                  ! Maintenance AR partition variables are stored as rates (kgC/plant/s)
                  ! (kgC/m2/yr) = (kgC/plant/s) * (plant/m2) * (s/yr)
                  hio_ar_agsapm_si_scpf(io_si,scpf) = hio_ar_agsapm_si_scpf(io_si,scpf) + &
                        ccohort%livestem_mr * n_perm2 * sec_per_day * days_per_year

                  ! (kgC/m2/yr) = (kgC/plant/s) * (plant/m2) * (s/yr)
                  hio_ar_darkm_si_scpf(io_si,scpf) = hio_ar_darkm_si_scpf(io_si,scpf) + &
                        ccohort%rdark * n_perm2 *  sec_per_day * days_per_year

                  ! (kgC/m2/yr) = (kgC/plant/s) * (plant/m2) * (s/yr)
                  hio_ar_crootm_si_scpf(io_si,scpf) = hio_ar_crootm_si_scpf(io_si,scpf) + &
                        ccohort%livecroot_mr * n_perm2 * sec_per_day * days_per_year

                  ! (kgC/m2/yr) = (kgC/plant/s) * (plant/m2) * (s/yr)
                  hio_ar_frootm_si_scpf(io_si,scpf) = hio_ar_frootm_si_scpf(io_si,scpf) + &
                        ccohort%froot_mr * n_perm2  * sec_per_day * days_per_year

                  ! accumulate fluxes per patch age bin
                  hio_gpp_si_age(io_si,cpatch%age_class) = hio_gpp_si_age(io_si,cpatch%age_class) &
                       + ccohort%gpp_tstep * ccohort%n * g_per_kg * per_dt_tstep
                  hio_npp_si_age(io_si,cpatch%age_class) = hio_npp_si_age(io_si,cpatch%age_class) &
                       + ccohort%npp_tstep * ccohort%n * g_per_kg * per_dt_tstep

                  ! accumulate fluxes on canopy- and understory- separated fluxes
                  if (ccohort%canopy_layer .eq. 1) then
                     !
                     ! bulk fluxes are in gC / m2 / s
                     hio_gpp_canopy_pa(io_pa) = hio_gpp_canopy_pa(io_pa) + &
                          ccohort%gpp_tstep * g_per_kg * n_density * per_dt_tstep                     
                     hio_ar_canopy_pa(io_pa) = hio_ar_canopy_pa(io_pa) + &
                          ccohort%resp_tstep * g_per_kg * n_density * per_dt_tstep                     
                     !
                     ! size-resolved respiration fluxes are in kg C / ha / yr
                     hio_rdark_canopy_si_scls(io_si,scls) = hio_rdark_canopy_si_scls(io_si,scls) + &
                          ccohort%rdark  * ccohort%n * sec_per_day * days_per_year
                     hio_livestem_mr_canopy_si_scls(io_si,scls) = hio_livestem_mr_canopy_si_scls(io_si,scls) + &
                          ccohort%livestem_mr  * ccohort%n * sec_per_day * days_per_year
                     hio_livecroot_mr_canopy_si_scls(io_si,scls) = hio_livecroot_mr_canopy_si_scls(io_si,scls) + &
                          ccohort%livecroot_mr  * ccohort%n * sec_per_day * days_per_year
                     hio_froot_mr_canopy_si_scls(io_si,scls) = hio_froot_mr_canopy_si_scls(io_si,scls) + &
                          ccohort%froot_mr  * ccohort%n * sec_per_day * days_per_year
                     hio_resp_g_canopy_si_scls(io_si,scls) = hio_resp_g_canopy_si_scls(io_si,scls) + &
                          ccohort%resp_g  * ccohort%n * sec_per_day * days_per_year * per_dt_tstep 
                     hio_resp_m_canopy_si_scls(io_si,scls) = hio_resp_m_canopy_si_scls(io_si,scls) + &
                          ccohort%resp_m  * ccohort%n * sec_per_day * days_per_year * per_dt_tstep 
                  else
                     !
                     ! bulk fluxes are in gC / m2 / s
                     hio_gpp_understory_pa(io_pa) = hio_gpp_understory_pa(io_pa) + &
                          ccohort%gpp_tstep * g_per_kg * n_density * per_dt_tstep                     
                     hio_ar_understory_pa(io_pa) = hio_ar_understory_pa(io_pa) + &
                          ccohort%resp_tstep * g_per_kg * n_density * per_dt_tstep                     
                     !
                     ! size-resolved respiration fluxes are in kg C / ha / yr
                     hio_rdark_understory_si_scls(io_si,scls) = hio_rdark_understory_si_scls(io_si,scls) + &
                          ccohort%rdark  * ccohort%n * sec_per_day * days_per_year
                     hio_livestem_mr_understory_si_scls(io_si,scls) = hio_livestem_mr_understory_si_scls(io_si,scls) + &
                          ccohort%livestem_mr  * ccohort%n * sec_per_day * days_per_year
                     hio_livecroot_mr_understory_si_scls(io_si,scls) = hio_livecroot_mr_understory_si_scls(io_si,scls) + &
                          ccohort%livecroot_mr  * ccohort%n * sec_per_day * days_per_year
                     hio_froot_mr_understory_si_scls(io_si,scls) = hio_froot_mr_understory_si_scls(io_si,scls) + &
                          ccohort%froot_mr  * ccohort%n * sec_per_day * days_per_year
                     hio_resp_g_understory_si_scls(io_si,scls) = hio_resp_g_understory_si_scls(io_si,scls) + &
                          ccohort%resp_g  * ccohort%n * sec_per_day * days_per_year * per_dt_tstep 
                     hio_resp_m_understory_si_scls(io_si,scls) = hio_resp_m_understory_si_scls(io_si,scls) + &
                          ccohort%resp_m  * ccohort%n * sec_per_day * days_per_year * per_dt_tstep 
                  endif
                end associate
               endif

               !!! resolve some canopy profile terms that are also on the cohort indices
               ican = ccohort%canopy_layer
               do ileaf=1,ccohort%nv
                  cnlf_indx = ileaf + (ican-1) * nlevleaf
                  hio_ts_net_uptake_si_cnlf(io_si, cnlf_indx) = hio_ts_net_uptake_si_cnlf(io_si, cnlf_indx) + &
                       ccohort%ts_net_uptake(ileaf) * ccohort%c_area / AREA
                  hio_year_net_uptake_si_cnlf(io_si, cnlf_indx) = hio_year_net_uptake_si_cnlf(io_si, cnlf_indx) + &
                       ccohort%year_net_uptake(ileaf) * ccohort%c_area / AREA
               end do

               ccohort => ccohort%taller
            enddo ! cohort loop

            ! summarize radiation profiles through the canopy
            do ipft=1,numpft
               do ican=1,nclmax
                  do ileaf=1,nlevleaf
                     ! calculate where we are on multiplexed dimensions
                     cnlfpft_indx = ileaf + (ican-1) * nlevleaf + (ipft-1) * nlevleaf * nclmax 
                     cnlf_indx = ileaf + (ican-1) * nlevleaf
                     !
                     ! first do all the canopy x leaf x pft calculations
                     hio_parsun_z_si_cnlfpft(io_si,cnlfpft_indx) = hio_parsun_z_si_cnlfpft(io_si,cnlfpft_indx) + &
                          cpatch%ed_parsun_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_parsha_z_si_cnlfpft(io_si,cnlfpft_indx) = hio_parsha_z_si_cnlfpft(io_si,cnlfpft_indx) + &
                          cpatch%ed_parsha_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     !
                     hio_laisun_z_si_cnlfpft(io_si,cnlfpft_indx) = hio_laisun_z_si_cnlfpft(io_si,cnlfpft_indx) + &
                          cpatch%ed_laisun_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_laisha_z_si_cnlfpft(io_si,cnlfpft_indx) = hio_laisha_z_si_cnlfpft(io_si,cnlfpft_indx) + &
                          cpatch%ed_laisha_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     !
                     hio_fabd_sun_si_cnlfpft(io_si,cnlfpft_indx) = hio_fabd_sun_si_cnlfpft(io_si,cnlfpft_indx) + &
                          cpatch%fabd_sun_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_fabd_sha_si_cnlfpft(io_si,cnlfpft_indx) = hio_fabd_sha_si_cnlfpft(io_si,cnlfpft_indx) + &
                          cpatch%fabd_sha_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_fabi_sun_si_cnlfpft(io_si,cnlfpft_indx) = hio_fabi_sun_si_cnlfpft(io_si,cnlfpft_indx) + &
                          cpatch%fabi_sun_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_fabi_sha_si_cnlfpft(io_si,cnlfpft_indx) = hio_fabi_sha_si_cnlfpft(io_si,cnlfpft_indx) + &
                          cpatch%fabi_sha_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     !
                     ! summarize across all PFTs
                     hio_parsun_z_si_cnlf(io_si,cnlf_indx) = hio_parsun_z_si_cnlf(io_si,cnlf_indx) + &
                          cpatch%ed_parsun_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_parsha_z_si_cnlf(io_si,cnlf_indx) = hio_parsha_z_si_cnlf(io_si,cnlf_indx) + &
                          cpatch%ed_parsha_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     !
                     hio_laisun_z_si_cnlf(io_si,cnlf_indx) = hio_laisun_z_si_cnlf(io_si,cnlf_indx) + &
                          cpatch%ed_laisun_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_laisha_z_si_cnlf(io_si,cnlf_indx) = hio_laisha_z_si_cnlf(io_si,cnlf_indx) + &
                          cpatch%ed_laisha_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     !
                     hio_fabd_sun_si_cnlf(io_si,cnlf_indx) = hio_fabd_sun_si_cnlf(io_si,cnlf_indx) + &
                          cpatch%fabd_sun_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_fabd_sha_si_cnlf(io_si,cnlf_indx) = hio_fabd_sha_si_cnlf(io_si,cnlf_indx) + &
                          cpatch%fabd_sha_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_fabi_sun_si_cnlf(io_si,cnlf_indx) = hio_fabi_sun_si_cnlf(io_si,cnlf_indx) + &
                          cpatch%fabi_sun_z(ican,ipft,ileaf) * cpatch%area * AREA_INV
                     hio_fabi_sha_si_cnlf(io_si,cnlf_indx) = hio_fabi_sha_si_cnlf(io_si,cnlf_indx) + &
                          cpatch%fabi_sha_z(ican,ipft,ileaf) * cpatch%area * AREA_INV

                  end do
                  !
                  ! summarize just the top leaf level across all PFTs, for each canopy level
                  hio_parsun_top_si_can(io_si,ican) = hio_parsun_top_si_can(io_si,ican) + &
                       cpatch%ed_parsun_z(ican,ipft,1) * cpatch%area * AREA_INV
                  hio_parsha_top_si_can(io_si,ican) = hio_parsha_top_si_can(io_si,ican) + &
                       cpatch%ed_parsha_z(ican,ipft,1) * cpatch%area * AREA_INV
                  !
                  hio_laisun_top_si_can(io_si,ican) = hio_laisun_top_si_can(io_si,ican) + &
                       cpatch%ed_laisun_z(ican,ipft,1) * cpatch%area * AREA_INV
                  hio_laisha_top_si_can(io_si,ican) = hio_laisha_top_si_can(io_si,ican) + &
                       cpatch%ed_laisha_z(ican,ipft,1) * cpatch%area * AREA_INV
                  !
                  hio_fabd_sun_top_si_can(io_si,ican) = hio_fabd_sun_top_si_can(io_si,ican) + &
                       cpatch%fabd_sun_z(ican,ipft,1) * cpatch%area * AREA_INV
                  hio_fabd_sha_top_si_can(io_si,ican) = hio_fabd_sha_top_si_can(io_si,ican) + &
                       cpatch%fabd_sha_z(ican,ipft,1) * cpatch%area * AREA_INV
                  hio_fabi_sun_top_si_can(io_si,ican) = hio_fabi_sun_top_si_can(io_si,ican) + &
                       cpatch%fabi_sun_z(ican,ipft,1) * cpatch%area * AREA_INV
                  hio_fabi_sha_top_si_can(io_si,ican) = hio_fabi_sha_top_si_can(io_si,ican) + &
                       cpatch%fabi_sha_z(ican,ipft,1) * cpatch%area * AREA_INV
                  !
               end do
            end do


            ipa = ipa + 1
            cpatch => cpatch%younger
         end do !patch loop

         do ipa2 = 1, nlevage
            if (patch_area_by_age(ipa2) .gt. tiny) then
               hio_gpp_si_age(io_si, ipa2) = hio_gpp_si_age(io_si, ipa2) / (patch_area_by_age(ipa2))
               hio_npp_si_age(io_si, ipa2) = hio_npp_si_age(io_si, ipa2) / (patch_area_by_age(ipa2))
            else
               hio_gpp_si_age(io_si, ipa2) = 0._r8
               hio_npp_si_age(io_si, ipa2) = 0._r8
            endif
         end do
         
      enddo ! site loop

    end associate
 
  end subroutine update_history_prod

  ! =====================================================================================

  subroutine update_history_hydraulics(this,nc,nsites,sites,dt_tstep)

    ! ---------------------------------------------------------------------------------
    ! This is the call to update the history IO arrays that are expected to only change
    ! after rapid timescale productivity calculations (gpp and respiration).
    ! ---------------------------------------------------------------------------------
    
    use EDtypesMod          , only : ed_site_type,   &
                                     ed_cohort_type, &
                                     ed_patch_type,  &
                                     AREA

    use FatesHydraulicsMemMod, only : ed_cohort_hydr_type
    use FatesHydraulicsMemMod, only : nlevsoi_hyd
    use EDTypesMod           , only : maxpft

    
    ! Arguments
    class(fates_history_interface_type)             :: this
    integer                 , intent(in)            :: nc   ! clump index
    integer                 , intent(in)            :: nsites
    type(ed_site_type)      , intent(inout), target :: sites(nsites)
    real(r8)                , intent(in)            :: dt_tstep
    
    ! Locals
    integer  :: s        ! The local site index
    integer  :: io_si     ! The site index of the IO array
    integer  :: ipa      ! The local "I"ndex of "PA"tches 
    integer  :: io_pa    ! The patch index of the IO array
    integer  :: io_pa1   ! The first patch index in the IO array for each site
    integer  :: ft               ! functional type index
    integer  :: scpf
    real(r8) :: n_density   ! individual of cohort per m2.
    real(r8) :: n_perm2     ! individuals per m2 for the whole column
    real(r8), parameter :: tiny = 1.e-5_r8      ! some small number
    real(r8) :: ncohort_scpf(nlevsclass*maxpft)  ! Bins to count up cohorts counts used in weighting
                                                   ! should be "hio_nplant_si_scpf"
    real(r8) :: number_fraction
    real(r8) :: number_fraction_rate
    integer  :: ipa2     ! patch incrementer
    integer  :: iscpf    ! index of the scpf group

    type(ed_patch_type),pointer  :: cpatch
    type(ed_cohort_type),pointer :: ccohort
    type(ed_cohort_hydr_type), pointer :: ccohort_hydr

    real(r8), parameter :: daysecs = 86400.0_r8 ! What modeler doesn't recognize 86400?
    real(r8), parameter :: yeardays = 365.0_r8  ! Should this be 365.25?
    
    if(hlm_use_planthydro.eq.ifalse) return

    associate( hio_errh2o_scpf  => this%hvars(ih_errh2o_scpf)%r82d, &
          hio_tran_scpf         => this%hvars(ih_tran_scpf)%r82d, &
          hio_rootuptake_scpf   => this%hvars(ih_rootuptake_scpf)%r82d, &
          hio_rootuptake01_scpf => this%hvars(ih_rootuptake01_scpf)%r82d, &
          hio_rootuptake02_scpf => this%hvars(ih_rootuptake02_scpf)%r82d, &
          hio_rootuptake03_scpf => this%hvars(ih_rootuptake03_scpf)%r82d, &
          hio_rootuptake04_scpf => this%hvars(ih_rootuptake04_scpf)%r82d, &
          hio_rootuptake05_scpf => this%hvars(ih_rootuptake05_scpf)%r82d, &
          hio_rootuptake06_scpf => this%hvars(ih_rootuptake06_scpf)%r82d, &
          hio_rootuptake07_scpf => this%hvars(ih_rootuptake07_scpf)%r82d, &
          hio_rootuptake08_scpf => this%hvars(ih_rootuptake08_scpf)%r82d, &
          hio_rootuptake09_scpf => this%hvars(ih_rootuptake09_scpf)%r82d, &
          hio_rootuptake10_scpf => this%hvars(ih_rootuptake10_scpf)%r82d, &
          hio_sapflow_scpf      => this%hvars(ih_sapflow_scpf)%r82d, &
          hio_iterh1_scpf       => this%hvars(ih_iterh1_scpf)%r82d, &          
          hio_iterh2_scpf       => this%hvars(ih_iterh2_scpf)%r82d, &           
          hio_ath_scpf          => this%hvars(ih_ath_scpf)%r82d, &               
          hio_tth_scpf          => this%hvars(ih_tth_scpf)%r82d, &               
          hio_sth_scpf          => this%hvars(ih_sth_scpf)%r82d, &                     
          hio_lth_scpf          => this%hvars(ih_lth_scpf)%r82d, &                     
          hio_awp_scpf          => this%hvars(ih_awp_scpf)%r82d, &                     
          hio_twp_scpf          => this%hvars(ih_twp_scpf)%r82d, &  
          hio_swp_scpf          => this%hvars(ih_swp_scpf)%r82d, &                     
          hio_lwp_scpf          => this%hvars(ih_lwp_scpf)%r82d, &                    
          hio_btran_scpf        => this%hvars(ih_btran_scpf)%r82d, &
          hio_nplant_si_scpf    => this%hvars(ih_nplant_si_scpf)%r82d )
      
      ! Flush the relevant history variables 
      call this%flush_hvars(nc,upfreq_in=4)

      do s = 1,nsites
         
         io_si  = this%iovar_map(nc)%site_index(s)
         io_pa1 = this%iovar_map(nc)%patch1_index(s)

         ncohort_scpf(:) = 0.0_r8  ! Counter for normalizing weighting 
                                   ! factors for cohort mean propoerties
                                   ! This is actually used as a check
                                   ! on hio_nplant_si_scpf

         cpatch => sites(s)%oldest_patch
         do while(associated(cpatch))
            ccohort => cpatch%shortest
            do while(associated(ccohort))
               if ( .not. ccohort%isnew ) then
                  ! Calculate index for the scpf class
                  iscpf = ccohort%size_by_pft_class
                  ncohort_scpf(iscpf) = ncohort_scpf(iscpf) + ccohort%n
               end if
               ccohort => ccohort%taller
            enddo ! cohort loop
            cpatch => cpatch%younger
         end do !patch loop
         

         ipa = 0
         cpatch => sites(s)%oldest_patch
         do while(associated(cpatch))
            
            io_pa = io_pa1 + ipa

            ccohort => cpatch%shortest
            do while(associated(ccohort))

               ccohort_hydr => ccohort%co_hydr
               
               ! TODO: we need a standardized logical function on this (used lots, RGK)
               if ((cpatch%area .gt. 0._r8) .and. (cpatch%total_canopy_area .gt. 0._r8)) then
                  n_density = ccohort%n/min(cpatch%area,cpatch%total_canopy_area) 
                  n_perm2   = ccohort%n/AREA   
               else
                  n_density = 0.0_r8
                  n_perm2   = 0.0_r8
               endif
               
               if ( .not. ccohort%isnew ) then

                  ! Calculate index for the scpf class
                  iscpf = ccohort%size_by_pft_class
                  
                  ! scale up cohort fluxes to their sites
                  number_fraction_rate = (ccohort%n / ncohort_scpf(iscpf))/dt_tstep
                  
                  ! scale cohorts to mean quantity
                  number_fraction = (ccohort%n / ncohort_scpf(iscpf))
                  
                  hio_errh2o_scpf(io_si,iscpf) = hio_errh2o_scpf(io_si,iscpf) + &
                        ccohort_hydr%errh2o * number_fraction_rate ! [kg/indiv/s]
                  
                  hio_tran_scpf(io_si,iscpf) = hio_tran_scpf(io_si,iscpf) + &
                        (ccohort_hydr%qtop_dt + ccohort_hydr%dqtopdth_dthdt) * number_fraction_rate ! [kg/indiv/s]
                  
                  hio_rootuptake_scpf(io_si,iscpf) = hio_rootuptake_scpf(io_si,iscpf) + &
                        ccohort_hydr%rootuptake * number_fraction_rate       ! [kg/indiv/s]
                  
                  if(nlevsoi_hyd == 10) then
                     hio_rootuptake01_scpf(io_si,iscpf) = hio_rootuptake01_scpf(io_si,iscpf) + &
                           ccohort_hydr%rootuptake01 * number_fraction_rate   ! [kg/indiv/s]

                     hio_rootuptake02_scpf(io_si,iscpf) = hio_rootuptake02_scpf(io_si,iscpf) + &
                           ccohort_hydr%rootuptake02 * number_fraction_rate     ! [kg/indiv/s]
                     
                     hio_rootuptake03_scpf(io_si,iscpf) = hio_rootuptake03_scpf(io_si,iscpf) + &
                           ccohort_hydr%rootuptake03 * number_fraction_rate     ! [kg/indiv/s]
                     
                     hio_rootuptake04_scpf(io_si,iscpf) = hio_rootuptake04_scpf(io_si,iscpf) + &
                           ccohort_hydr%rootuptake04 * number_fraction_rate     ! [kg/indiv/s]

                     hio_rootuptake05_scpf(io_si,iscpf) = hio_rootuptake05_scpf(io_si,iscpf) + &
                           ccohort_hydr%rootuptake05 * number_fraction_rate     ! [kg/indiv/s]
                     
                     hio_rootuptake06_scpf(io_si,iscpf) = hio_rootuptake06_scpf(io_si,iscpf) + &
                           ccohort_hydr%rootuptake06 * number_fraction_rate     ! [kg/indiv/s]
                     
                     hio_rootuptake07_scpf(io_si,iscpf) = hio_rootuptake07_scpf(io_si,iscpf) + &
                             ccohort_hydr%rootuptake07 * number_fraction_rate    ! [kg/indiv/s]
                     
                     hio_rootuptake08_scpf(io_si,iscpf) = hio_rootuptake08_scpf(io_si,iscpf) + &
                           ccohort_hydr%rootuptake08 * number_fraction_rate     ! [kg/indiv/s]
                     
                     hio_rootuptake09_scpf(io_si,iscpf) = hio_rootuptake09_scpf(io_si,iscpf) + &
                           ccohort_hydr%rootuptake09 * number_fraction_rate    ! [kg/indiv/s] 
                     
                     hio_rootuptake10_scpf(io_si,iscpf) = hio_rootuptake10_scpf(io_si,iscpf) + &
                           ccohort_hydr%rootuptake10 * number_fraction_rate     ! [kg/indiv/s]
                     
                  end if
                  
                  hio_sapflow_scpf(io_si,iscpf)         = hio_sapflow_scpf(io_si,iscpf)  + &
                        ccohort_hydr%sapflow * number_fraction_rate             ! [kg/indiv/s]
                  
                  hio_iterh1_scpf(io_si,iscpf)          = hio_iterh1_scpf(io_si,iscpf) + &
                        ccohort_hydr%iterh1  * number_fraction             ! [-]
                  
                  hio_iterh2_scpf(io_si,iscpf)          = hio_iterh2_scpf(io_si,iscpf) + &
                        ccohort_hydr%iterh2 * number_fraction             ! [-]
                  
                  hio_ath_scpf(io_si,iscpf)             = hio_ath_scpf(io_si,iscpf) + &
                        ccohort_hydr%th_aroot(1)   * number_fraction      ! [m3 m-3]
                  
                  hio_tth_scpf(io_si,iscpf)             = hio_tth_scpf(io_si,iscpf) + &
                        ccohort_hydr%th_bg(1)  * number_fraction         ! [m3 m-3]
                  
                  hio_sth_scpf(io_si,iscpf)             = hio_sth_scpf(io_si,iscpf) + &
                        ccohort_hydr%th_ag(2)  * number_fraction        ! [m3 m-3]
                  
                  hio_lth_scpf(io_si,iscpf)             =  hio_lth_scpf(io_si,iscpf) + &
                        ccohort_hydr%th_ag(1)  * number_fraction        ! [m3 m-3]
                  
                  hio_awp_scpf(io_si,iscpf)             = hio_awp_scpf(io_si,iscpf) + &
                        ccohort_hydr%psi_aroot(1)   * number_fraction     ! [MPa]
                  
                  hio_twp_scpf(io_si,iscpf)             = hio_twp_scpf(io_si,iscpf) + &
                        ccohort_hydr%psi_bg(1)  * number_fraction       ! [MPa]
                  
                  hio_swp_scpf(io_si,iscpf)             = hio_swp_scpf(io_si,iscpf) + &
                        ccohort_hydr%psi_ag(2)  * number_fraction       ! [MPa]
                  
                  hio_lwp_scpf(io_si,iscpf)             = hio_lwp_scpf(io_si,iscpf) + &
                        ccohort_hydr%psi_ag(1)  * number_fraction       ! [MPa]
                  
                  hio_btran_scpf(io_si,iscpf)           = hio_btran_scpf(io_si,iscpf) + &
                        ccohort_hydr%btran(1)  * number_fraction        ! [-]
                  
               endif

               ccohort => ccohort%taller
            enddo ! cohort loop
            ipa = ipa + 1
            cpatch => cpatch%younger
         end do !patch loop

         if(hlm_use_ed_st3.eq.ifalse) then
            do scpf=1,nlevsclass*numpft
               if( abs(hio_nplant_si_scpf(io_si, scpf)-ncohort_scpf(scpf)) > 1.0E-8_r8 ) then
                  write(fates_log(),*) 'nplant check on hio_nplant_si_scpf fails during hydraulics history updates'
                  call endrun(msg=errMsg(sourcefile, __LINE__))
               end if
            end do
         end if

      enddo ! site loop

    end associate
 
 end subroutine update_history_hydraulics

  ! ====================================================================================
  integer function num_history_vars(this)

    implicit none

    class(fates_history_interface_type), intent(in) :: this

    num_history_vars = this%num_history_vars_
    
  end function num_history_vars
  
  ! ====================================================================================
  
  subroutine initialize_history_vars(this)

    implicit none

    class(fates_history_interface_type), intent(inout) :: this

   ! Determine how many of the history IO variables registered in FATES
   ! are going to be allocated
   call this%define_history_vars(initialize_variables=.false.)

   ! Allocate the list of history output variable objects
   allocate(this%hvars(this%num_history_vars()))
   
   ! construct the object that defines all of the IO variables
   call this%define_history_vars(initialize_variables=.true.)
   
 end subroutine initialize_history_vars
  
  ! ====================================================================================
  
  subroutine define_history_vars(this, initialize_variables)
    
    ! ---------------------------------------------------------------------------------
    ! 
    !                    REGISTRY OF HISTORY OUTPUT VARIABLES
    !
    ! This subroutine is called in two contexts, either in count mode or inialize mode
    ! In count mode, we just walk through the list of registerred variables, compare
    ! if the variable of interest list the current host model and add it to the count
    ! if true.  This count is used just to allocate the variable space.  After this
    ! has been done, we go through the list a second time populating a memory structure.
    ! This phase is the "initialize" phase.  These two phases are differntiated by the
    ! string "callstep", which should be either "count" or "initialize".
    !
    ! Note 1 there are different ways you can flush or initialize the output fields.
    ! If you flush to a native type, (such as zero), the entire slab which covers
    ! indices which may not be relevant to FATES, are flushed to this value.  So
    ! in that case, lakes and crops that are not controlled by FATES will zero'd
    ! and when values are scaled up to the land-grid, the zero's for non FATES will
    ! be included.  This is good and correct if nothing is there.  
    !
    ! But, what if crops exist in the host model and occupy a fraction of the land-surface
    ! shared with natural vegetation? In that case, you want to flush your arrays
    ! with a value that the HLM treats as "do not average"
    ! 
    ! If your HLM makes use of, and you want, INTEGER OUTPUT, pass the flushval as
    ! a real.  The applied flush value will use the NINT() intrinsic function
    ! ---------------------------------------------------------------------------------

    use FatesIOVariableKindMod, only : patch_r8, patch_ground_r8, patch_size_pft_r8
    use FatesIOVariableKindMod, only : site_r8, site_ground_r8, site_size_pft_r8    
    use FatesIOVariableKindMod, only : site_size_r8, site_pft_r8, site_age_r8
    use FatesInterfaceMod     , only : hlm_use_planthydro
    
    use FatesIOVariableKindMod, only : site_fuel_r8, site_cwdsc_r8, site_scag_r8
    use FatesIOVariableKindMod, only : site_can_r8, site_cnlf_r8, site_cnlfpft_r8

    implicit none
    
    class(fates_history_interface_type), intent(inout) :: this
    logical, intent(in) :: initialize_variables  ! are we 'count'ing or 'initializ'ing?

    integer :: ivar
    character(len=10) :: tempstring 
    
    ivar=0
    
    ! Site level counting variables
    call this%set_history_var(vname='ED_NPATCHES', units='none',                &
         long='Total number of ED patches per site', use_default='active',      &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=1.0_r8, upfreq=1,    &
         ivar=ivar, initialize=initialize_variables, index = ih_npatches_si)

    call this%set_history_var(vname='ED_NCOHORTS', units='none',                &
         long='Total number of ED cohorts per site', use_default='active',      &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=1.0_r8, upfreq=1,    &
         ivar=ivar, initialize=initialize_variables, index = ih_ncohorts_si)
    
    ! Patch variables
    call this%set_history_var(vname='TRIMMING', units='none',                   &
         long='Degree to which canopy expansion is limited by leaf economics',  & 
         use_default='active', &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=1.0_r8, upfreq=1,    &
         ivar=ivar, initialize=initialize_variables, index = ih_trimming_pa)
    
    call this%set_history_var(vname='AREA_PLANT', units='m2',                   &
         long='area occupied by all plants', use_default='active',              &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,    &
         ivar=ivar, initialize=initialize_variables, index = ih_area_plant_pa)
    
    call this%set_history_var(vname='AREA_TREES', units='m2',                   &
         long='area occupied by woody plants', use_default='active',            &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,    &
         ivar=ivar, initialize=initialize_variables, index = ih_area_treespread_pa)

    call this%set_history_var(vname='CANOPY_SPREAD', units='0-1',               &
         long='Scaling factor between tree basal area and canopy area',         &
         use_default='active',                                                  &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,    &
         ivar=ivar, initialize=initialize_variables, index = ih_canopy_spread_si)

    call this%set_history_var(vname='PFTbiomass', units='gC/m2',                   &
         long='total PFT level biomass', use_default='active',                     &
         avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_biomass_si_pft )

    call this%set_history_var(vname='PFTleafbiomass', units='gC/m2',              &
         long='total PFT level leaf biomass', use_default='active',                &
         avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_leafbiomass_si_pft )

    call this%set_history_var(vname='PFTstorebiomass',  units='gC/m2',            &
         long='total PFT level stored biomass', use_default='active',              &
         avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_storebiomass_si_pft )
    
    call this%set_history_var(vname='PFTnindivs',  units='indiv / m2',            &
         long='total PFT level number of individuals', use_default='active',       &
         avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_nindivs_si_pft )

    call this%set_history_var(vname='RECRUITMENT',  units='indiv/ha/yr',            &
         long='Rate of recruitment by PFT', use_default='active',       &
         avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_recruitment_si_pft )

    call this%set_history_var(vname='MORTALITY',  units='indiv/ha/yr',            &
         long='Rate of total mortality by PFT', use_default='active',       &
         avgflag='A', vtype=site_pft_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_mortality_si_pft )

    ! patch age class variables
    call this%set_history_var(vname='PATCH_AREA_BY_AGE', units='m2/m2',             &
         long='patch area by age bin', use_default='active',                     &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_area_si_age )

    call this%set_history_var(vname='LAI_BY_AGE', units='m2/m2',                   &
         long='leaf area index by age bin', use_default='active',                     &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_lai_si_age )

    call this%set_history_var(vname='CANOPY_AREA_BY_AGE', units='m2/m2',             &
         long='canopy area by age bin', use_default='active',                     &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_canopy_area_si_age )
    
    call this%set_history_var(vname='NCL_BY_AGE', units='--',                   &
         long='number of canopy levels by age bin', use_default='inactive',             &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_ncl_si_age )

    call this%set_history_var(vname='NPATCH_BY_AGE', units='--',                   &
         long='number of patches by age bin', use_default='inactive',                     &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_npatches_si_age )

    if ( ED_val_comp_excln .lt. 0._r8 ) then ! only valid when "strict ppa" enabled
       tempstring = 'active'
    else
       tempstring = 'inactive'
    endif
    call this%set_history_var(vname='ZSTAR_BY_AGE', units='m',                   &
         long='product of zstar and patch area by age bin (divide by PATCH_AREA_BY_AGE to get mean zstar)', &
         use_default=trim(tempstring),                     &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_zstar_si_age )

    call this%set_history_var(vname='BIOMASS_BY_AGE', units='m',                   &
         long='Total Biomass within a given patch age bin (kg C)', &
         use_default='inactive',                     &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1, &
         ivar=ivar, initialize=initialize_variables, index = ih_biomass_si_age )

    ! Fire Variables

    call this%set_history_var(vname='FIRE_NESTEROV_INDEX', units='none',       &
         long='nesterov_fire_danger index', use_default='active',               &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_nesterov_fire_danger_pa)

    call this%set_history_var(vname='FIRE_ROS', units='m/min',                 &
         long='fire rate of spread m/min', use_default='active',                &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_spitfire_ROS_pa)

    call this%set_history_var(vname='EFFECT_WSPEED', units='none',             &
         long ='effective windspeed for fire spread', use_default='active',     &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_effect_wspeed_pa )

    call this%set_history_var(vname='FIRE_TFC_ROS', units='none',              &
         long ='total fuel consumed', use_default='active',                     &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_TFC_ROS_pa )

    call this%set_history_var(vname='FIRE_INTENSITY', units='kJ/m/s',          &
         long='spitfire fire intensity: kJ/m/s', use_default='active',          &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fire_intensity_pa )

    call this%set_history_var(vname='FIRE_AREA', units='fraction',             &
         long='spitfire fire area:m2', use_default='active',                    &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fire_area_pa )

    call this%set_history_var(vname='SCORCH_HEIGHT', units='m',                &
         long='spitfire fire area:m2', use_default='active',                    &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_scorch_height_pa )

    call this%set_history_var(vname='fire_fuel_mef', units='m',                &
         long='spitfire fuel moisture',  use_default='active',                  &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fire_fuel_mef_pa )

    call this%set_history_var(vname='fire_fuel_bulkd', units='m',              &
         long='spitfire fuel bulk density',  use_default='active',              &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fire_fuel_bulkd_pa )

    call this%set_history_var(vname='FIRE_FUEL_EFF_MOIST', units='m',          &
         long='spitfire fuel moisture', use_default='active',                   &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fire_fuel_eff_moist_pa )

    call this%set_history_var(vname='fire_fuel_sav', units='m',                &
         long='spitfire fuel surface/volume ',  use_default='active',           &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fire_fuel_sav_pa )

    call this%set_history_var(vname='SUM_FUEL', units='gC m-2',                &
         long='total ground fuel related to ros (omits 1000hr fuels)',          & 
         use_default='active',                                                  & 
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_sum_fuel_pa )

    call this%set_history_var(vname='FUEL_MOISTURE_NFSC', units='-',                &
         long='spitfire size-resolved fuel moisture', use_default='active',       &
         avgflag='A', vtype=site_fuel_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_litter_moisture_si_fuel )

    ! Litter Variables

    call this%set_history_var(vname='LITTER_IN', units='gC m-2 s-1',           &
         long='FATES litter flux in',  use_default='active',                   &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_litter_in_si )

    call this%set_history_var(vname='LITTER_OUT', units='gC m-2 s-1',          &
         long='FATES litter flux out',  use_default='active',                  & 
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_litter_out_pa )

    call this%set_history_var(vname='SEED_BANK', units='gC m-2',               &
         long='Total Seed Mass of all PFTs',  use_default='active',             &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_seed_bank_si )

    call this%set_history_var(vname='SEEDS_IN', units='gC m-2 s-1',            &
         long='Seed Production Rate',  use_default='active',                    &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_seeds_in_pa )

    call this%set_history_var(vname='SEED_GERMINATION', units='gC m-2 s-1',    &
         long='Seed mass converted into new cohorts',   use_default='active',   &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_seed_germination_pa )

    call this%set_history_var(vname='SEED_DECAY', units='gC m-2 s-1',          &
         long='Seed mass decay', use_default='active',                          &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_seed_decay_pa )
    
    call this%set_history_var(vname='ED_bstore', units='gC m-2',                  &
         long='Storage biomass', use_default='active',                          &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_bstore_pa )

    call this%set_history_var(vname='ED_bdead', units='gC m-2',                   &
         long='Dead (structural) biomass (live trees, not CWD)',                &
         use_default='active',                                                  &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_bdead_pa )

    call this%set_history_var(vname='ED_balive', units='gC m-2',                  &
         long='Live biomass', use_default='active',                             &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_balive_pa )

    call this%set_history_var(vname='ED_bleaf', units='gC m-2',                   &
         long='Leaf biomass',  use_default='active',                            &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_bleaf_pa )

    call this%set_history_var(vname='ED_biomass', units='gC m-2',                  &
         long='Total biomass',  use_default='active',                           &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_btotal_pa )

    call this%set_history_var(vname='BIOMASS_CANOPY', units='gC m-2',                   &
         long='Biomass of canopy plants',  use_default='active',                            &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_canopy_biomass_pa )

    call this%set_history_var(vname='BIOMASS_UNDERSTORY', units='gC m-2',                   &
         long='Biomass of understory plants',  use_default='active',                            &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_understory_biomass_pa )

    
    ! Ecosystem Carbon Fluxes (updated rapidly, upfreq=2)

    call this%set_history_var(vname='NPP_column', units='gC/m^2/s',                &
         long='net primary production on the site',  use_default='active',      &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_npp_si )

    call this%set_history_var(vname='GPP', units='gC/m^2/s',                   &
         long='gross primary production',  use_default='active',                &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_gpp_pa )

    call this%set_history_var(vname='NPP', units='gC/m^2/s',                   &
         long='net primary production', use_default='active',                   &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_npp_pa )

    call this%set_history_var(vname='AR', units='gC/m^2/s',                 &
         long='autotrophic respiration', use_default='active',                  &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_aresp_pa )

    call this%set_history_var(vname='GROWTH_RESP', units='gC/m^2/s',           &
         long='growth respiration', use_default='active',                       &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_growth_resp_pa )

    call this%set_history_var(vname='MAINT_RESP', units='gC/m^2/s',            &
         long='maintenance respiration', use_default='active',                  &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_maint_resp_pa )

    ! fast fluxes by age bin
    call this%set_history_var(vname='NPP_BY_AGE', units='gC/m^2/s',                   &
         long='net primary productivity by age bin', use_default='inactive',           &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2, &
         ivar=ivar, initialize=initialize_variables, index = ih_npp_si_age )

    call this%set_history_var(vname='GPP_BY_AGE', units='gC/m^2/s',                   &
         long='gross primary productivity by age bin', use_default='inactive',         &
         avgflag='A', vtype=site_age_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2, &
         ivar=ivar, initialize=initialize_variables, index = ih_gpp_si_age )

    ! fast fluxes separated canopy/understory
    call this%set_history_var(vname='GPP_CANOPY', units='gC/m^2/s',                   &
         long='gross primary production of canopy plants',  use_default='active',     &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_gpp_canopy_pa )

    call this%set_history_var(vname='AR_CANOPY', units='gC/m^2/s',                 &
         long='autotrophic respiration of canopy plants', use_default='active',       &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_ar_canopy_pa )

    call this%set_history_var(vname='GPP_UNDERSTORY', units='gC/m^2/s',                   &
         long='gross primary production of understory plants',  use_default='active',     &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_gpp_understory_pa )

    call this%set_history_var(vname='AR_UNDERSTORY', units='gC/m^2/s',                 &
         long='autotrophic respiration of understory plants', use_default='active',       &
         avgflag='A', vtype=patch_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_ar_understory_pa )


    ! fast radiative fluxes resolved through the canopy
    call this%set_history_var(vname='PARSUN_Z_CNLF', units='W/m2',                 &
         long='PAR absorbed in the sun by each canopy and leaf layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_parsun_z_si_cnlf )

    call this%set_history_var(vname='PARSHA_Z_CNLF', units='W/m2',                 &
         long='PAR absorbed in the shade by each canopy and leaf layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_parsha_z_si_cnlf )

    call this%set_history_var(vname='PARSUN_Z_CNLFPFT', units='W/m2',                 &
         long='PAR absorbed in the sun by each canopy, leaf, and PFT', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlfpft_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_parsun_z_si_cnlfpft )

    call this%set_history_var(vname='PARSHA_Z_CNLFPFT', units='W/m2',                 &
         long='PAR absorbed in the shade by each canopy, leaf, and PFT', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlfpft_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_parsha_z_si_cnlfpft )

    call this%set_history_var(vname='PARSUN_Z_CAN', units='W/m2',                 &
         long='PAR absorbed in the sun by top leaf layer in each canopy layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_can_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_parsun_top_si_can )

    call this%set_history_var(vname='PARSHA_Z_CAN', units='W/m2',                 &
         long='PAR absorbed in the shade by top leaf layer in each canopy layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_can_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_parsha_top_si_can )

    call this%set_history_var(vname='LAISUN_Z_CNLF', units='m2/m2',                 &
         long='LAI in the sun by each canopy and leaf layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_laisun_z_si_cnlf )

    call this%set_history_var(vname='LAISHA_Z_CNLF', units='m2/m2',                 &
         long='LAI in the shade by each canopy and leaf layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_laisha_z_si_cnlf )

    call this%set_history_var(vname='LAISUN_Z_CNLFPFT', units='m2/m2',                 &
         long='LAI in the sun by each canopy, leaf, and PFT', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlfpft_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_laisun_z_si_cnlfpft )

    call this%set_history_var(vname='LAISHA_Z_CNLFPFT', units='m2/m2',                 &
         long='LAI in the shade by each canopy, leaf, and PFT', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlfpft_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_laisha_z_si_cnlfpft )

    call this%set_history_var(vname='LAISUN_TOP_CAN', units='m2/m2',                 &
         long='LAI in the sun by the top leaf layer of each canopy layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_can_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_laisun_top_si_can )

    call this%set_history_var(vname='LAISHA_TOP_CAN', units='m2/m2',                 &
         long='LAI in the shade by the top leaf layer of each canopy layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_can_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_laisha_top_si_can )

    call this%set_history_var(vname='FABD_SUN_CNLFPFT', units='fraction',                 &
         long='sun fraction of direct light absorbed by each canopy, leaf, and PFT', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlfpft_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabd_sun_si_cnlfpft )

    call this%set_history_var(vname='FABD_SHA_CNLFPFT', units='fraction',                 &
         long='shade fraction of direct light absorbed by each canopy, leaf, and PFT', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlfpft_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabd_sha_si_cnlfpft )

    call this%set_history_var(vname='FABI_SUN_CNLFPFT', units='fraction',                 &
         long='sun fraction of indirect light absorbed by each canopy, leaf, and PFT', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlfpft_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabi_sun_si_cnlfpft )

    call this%set_history_var(vname='FABI_SHA_CNLFPFT', units='fraction',                 &
         long='shade fraction of indirect light absorbed by each canopy, leaf, and PFT', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlfpft_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabi_sha_si_cnlfpft )

    call this%set_history_var(vname='FABD_SUN_CNLF', units='fraction',                 &
         long='sun fraction of direct light absorbed by each canopy and leaf layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabd_sun_si_cnlf )

    call this%set_history_var(vname='FABD_SHA_CNLF', units='fraction',                 &
         long='shade fraction of direct light absorbed by each canopy and leaf layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabd_sha_si_cnlf )

    call this%set_history_var(vname='FABI_SUN_CNLF', units='fraction',                 &
         long='sun fraction of indirect light absorbed by each canopy and leaf layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabi_sun_si_cnlf )

    call this%set_history_var(vname='FABI_SHA_CNLF', units='fraction',                 &
         long='shade fraction of indirect light absorbed by each canopy and leaf layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabi_sha_si_cnlf )

    call this%set_history_var(vname='FABD_SUN_TOPLF_BYCANLAYER', units='fraction',                 &
         long='sun fraction of direct light absorbed by the top leaf layer of each canopy layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_can_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabd_sun_top_si_can )

    call this%set_history_var(vname='FABD_SHA_TOPLF_BYCANLAYER', units='fraction',                 &
         long='shade fraction of direct light absorbed by the top leaf layer of each canopy layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_can_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabd_sha_top_si_can )

    call this%set_history_var(vname='FABI_SUN_TOPLF_BYCANLAYER', units='fraction',                 &
         long='sun fraction of indirect light absorbed by the top leaf layer of each canopy layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_can_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabi_sun_top_si_can )

    call this%set_history_var(vname='FABI_SHA_TOPLF_BYCANLAYER', units='fraction',                 &
         long='shade fraction of indirect light absorbed by the top leaf layer of each canopy layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_can_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_fabi_sha_top_si_can )

    !!! canopy-resolved fluxes and structure
    call this%set_history_var(vname='TS_NET_UPTAKE_CNLF', units='kgC/m2/s',                 &
         long='net carbon uptake by each canopy and leaf layer er unit ground area (i.e. divide by CROWNAREA_CNLF)', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_ts_net_uptake_si_cnlf )

    call this%set_history_var(vname='YEAR_NET_UPTAKE_CNLF', units='kgC/m2/y',                 &
         long='yearly net carbon uptake by each canopy and leaf layer per unit ground area (i.e. divide by CROWNAREA_CNLF)', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=2,   &
         ivar=ivar, initialize=initialize_variables, index = ih_year_net_uptake_si_cnlf )

    call this%set_history_var(vname='CROWNAREA_CNLF', units='m2/m2',                 &
         long='total crown area that is occupied by leaves in each canopy and leaf layer', &
         use_default='inactive',       &
         avgflag='A', vtype=site_cnlf_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_crownarea_si_cnlf )

    call this%set_history_var(vname='CROWNAREA_CAN', units='m2/m2',                 &
         long='total crown area in each canopy layer', &
         use_default='active',       &
         avgflag='A', vtype=site_can_r8, hlms='CLM:ALM', flushval=0.0_r8, upfreq=1,   &
         ivar=ivar, initialize=initialize_variables, index = ih_crownarea_si_can )

    ! slow carbon fluxes associated with mortality from or transfer betweeen canopy and understory
    call this%set_history_var(vname='DEMOTION_CARBONFLUX', units = 'gC/m2/s',               &
          long='demotion-associated biomass carbon flux from canopy to understory', use_default='active',   &
          avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_demotion_carbonflux_si )

    call this%set_history_var(vname='PROMOTION_CARBONFLUX', units = 'gC/m2/s',               &
          long='promotion-associated biomass carbon flux from understory to canopy', use_default='active',   &
          avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_promotion_carbonflux_si )

    call this%set_history_var(vname='MORTALITY_CARBONFLUX_CANOPY', units = 'gC/m2/s',               &
          long='flux of biomass carbon from live to dead pools from mortality of canopy plants', use_default='inactive',   &
          avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_canopy_mortality_carbonflux_si )

    call this%set_history_var(vname='MORTALITY_CARBONFLUX_UNDERSTORY', units = 'gC/m2/s',               &
          long='flux of biomass carbon from live to dead pools from mortality of understory plants',use_default='inactive',&
          avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_understory_mortality_carbonflux_si )


    call this%set_history_var(vname='NPLANT_SCAG',units = 'plants/ha',               &
          long='number of plants per hectare in each size x age class', use_default='inactive',   &
          avgflag='A', vtype=site_scag_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_si_scag )

    call this%set_history_var(vname='NPLANT_CANOPY_SCAG',units = 'plants/ha',               &
          long='number of plants per hectare in canopy in each size x age class', use_default='inactive',   &
          avgflag='A', vtype=site_scag_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_canopy_si_scag )

    call this%set_history_var(vname='NPLANT_UNDERSTORY_SCAG',units = 'plants/ha',               &
          long='number of plants per hectare in understory in each size x age class', use_default='inactive',   &
          avgflag='A', vtype=site_scag_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_understory_si_scag )

    call this%set_history_var(vname='DDBH_CANOPY_SCAG',units = 'cm/yr/ha',               &
          long='growth rate of canopy plantsnumber of plants per hectare in canopy in each size x age class', use_default='inactive',   &
          avgflag='A', vtype=site_scag_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ddbh_canopy_si_scag )

    call this%set_history_var(vname='DDBH_UNDERSTORY_SCAG',units = 'cm/yr/ha',               &
          long='growth rate of understory plants in each size x age class', use_default='inactive',   &
          avgflag='A', vtype=site_scag_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ddbh_understory_si_scag )

    call this%set_history_var(vname='MORTALITY_CANOPY_SCAG',units = 'plants/ha/yr',               &
          long='mortality rate of canopy plants in each size x age class', use_default='inactive',   &
          avgflag='A', vtype=site_scag_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_mortality_canopy_si_scag )

    call this%set_history_var(vname='MORTALITY_UNDERSTORY_SCAG',units = 'plants/ha/yr',               &
          long='mortality rate of understory plantsin each size x age class', use_default='inactive',   &
          avgflag='A', vtype=site_scag_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_mortality_understory_si_scag )


    ! Carbon Flux (grid dimension x scpf) (THESE ARE DEFAULT INACTIVE!!!
    !                                     (BECAUSE THEY TAKE UP SPACE!!!
    ! ===================================================================================

    call this%set_history_var(vname='GPP_SCPF', units='kgC/m2/yr',            &
          long='gross primary production by pft/size', use_default='inactive',           &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_gpp_si_scpf )

    call this%set_history_var(vname='GPP_CANOPY_SCPF', units='kgC/m2/yr',            &
          long='gross primary production of canopy plants by pft/size ', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_gpp_canopy_si_scpf )

    call this%set_history_var(vname='AR_CANOPY_SCPF', units='kgC/m2/yr',            &
          long='autotrophic respiration of canopy plants by pft/size', use_default='inactive',           &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ar_canopy_si_scpf )

    call this%set_history_var(vname='GPP_UNDERSTORY_SCPF', units='kgC/m2/yr',            &
          long='gross primary production of understory plants by pft/size', use_default='inactive',           &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_gpp_understory_si_scpf )

    call this%set_history_var(vname='AR_UNDERSTORY_SCPF', units='kgC/m2/yr',            &
          long='autotrophic respiration of understory plants by pft/size', use_default='inactive',           &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ar_understory_si_scpf )

    call this%set_history_var(vname='NPP_SCPF', units='kgC/m2/yr',            &
          long='total net primary production by pft/size', use_default='inactive',       &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_totl_si_scpf )


    call this%set_history_var(vname='NPP_LEAF_SCPF', units='kgC/m2/yr',       &
          long='NPP flux into leaves by pft/size', use_default='inactive',               &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_leaf_si_scpf )

    call this%set_history_var(vname='NPP_SEED_SCPF', units='kgC/m2/yr',       &
          long='NPP flux into seeds by pft/size', use_default='inactive',                &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_seed_si_scpf )

    call this%set_history_var(vname='NPP_FNRT_SCPF', units='kgC/m2/yr',       &
          long='NPP flux into fine roots by pft/size', use_default='inactive',           &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_fnrt_si_scpf )

    call this%set_history_var(vname='NPP_BGSW_SCPF', units='kgC/m2/yr',       &
          long='NPP flux into below-ground sapwood by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_bgsw_si_scpf )

    call this%set_history_var(vname='NPP_BGDW_SCPF', units='kgC/m2/yr',       &
          long='NPP flux into below-ground deadwood by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_bgdw_si_scpf )

    call this%set_history_var(vname='NPP_AGSW_SCPF', units='kgC/m2/yr',       &
          long='NPP flux into above-ground sapwood by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_agsw_si_scpf )

    call this%set_history_var(vname = 'NPP_AGDW_SCPF', units='kgC/m2/yr',    &
          long='NPP flux into above-ground deadwood by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_agdw_si_scpf )

    call this%set_history_var(vname = 'NPP_STOR_SCPF', units='kgC/m2/yr',    &
          long='NPP flux into storage by pft/size', use_default='inactive',              &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_stor_si_scpf )

    call this%set_history_var(vname='DDBH_SCPF', units = 'cm/yr/ha',         &
          long='diameter growth increment by pft/size',use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ddbh_si_scpf )

    call this%set_history_var(vname='DDBH_CANOPY_SCPF', units = 'cm/yr/ha',         &
          long='diameter growth increment by pft/size',use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ddbh_canopy_si_scpf )

    call this%set_history_var(vname='DDBH_UNDERSTORY_SCPF', units = 'cm/yr/ha',         &
          long='diameter growth increment by pft/size',use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ddbh_understory_si_scpf )

    call this%set_history_var(vname='BA_SCPF', units = 'm2/ha',               &
          long='basal area by pft/size', use_default='inactive',   &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ba_si_scpf )

    call this%set_history_var(vname='NPLANT_SCPF', units = 'N/ha',         &
          long='stem number density by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_si_scpf )

    call this%set_history_var(vname='M1_SCPF', units = 'N/ha/yr',          &
          long='background mortality by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m1_si_scpf )
    
    call this%set_history_var(vname='M2_SCPF', units = 'N/ha/yr',          &
          long='hydraulic mortality by pft/size',use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m2_si_scpf )

    call this%set_history_var(vname='M3_SCPF', units = 'N/ha/yr',          &
          long='carbon starvation mortality by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m3_si_scpf )

    call this%set_history_var(vname='M4_SCPF', units = 'N/ha/yr',          &
          long='impact mortality by pft/size',use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m4_si_scpf )

    call this%set_history_var(vname='M5_SCPF', units = 'N/ha/yr',          &
          long='fire mortality by pft/size',use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m5_si_scpf )

    call this%set_history_var(vname='M6_SCPF', units = 'N/ha/yr',          &
          long='termination mortality by pft/size',use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m6_si_scpf )


    !Logging
    call this%set_history_var(vname='M7_SCPF', units = 'N/ha/event',               &
          long='logging mortalities by pft/size',use_default='inactive',           &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_m7_si_scpf )


    call this%set_history_var(vname='MORTALITY_CANOPY_SCPF', units = 'N/ha/yr',          &
          long='total mortality of canopy plants by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_mortality_canopy_si_scpf )

    call this%set_history_var(vname='BSTOR_CANOPY_SCPF', units = 'kgC/ha',          &
          long='biomass carbon in storage pools of canopy plants by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_bstor_canopy_si_scpf )

    call this%set_history_var(vname='BLEAF_CANOPY_SCPF', units = 'kgC/ha',          &
          long='biomass carbon in leaf of canopy plants by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_bleaf_canopy_si_scpf )

    call this%set_history_var(vname='NPLANT_CANOPY_SCPF', units = 'N/ha',         &
          long='stem number of canopy plants density by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_canopy_si_scpf )

    call this%set_history_var(vname='MORTALITY_UNDERSTORY_SCPF', units = 'N/ha/yr',          &
          long='total mortality of understory plants by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_mortality_understory_si_scpf )

    call this%set_history_var(vname='BSTOR_UNDERSTORY_SCPF', units = 'kgC/ha',          &
          long='biomass carbon in storage pools of understory plants by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_bstor_understory_si_scpf )

    call this%set_history_var(vname='BLEAF_UNDERSTORY_SCPF', units = 'kgC/ha',          &
          long='biomass carbon in leaf of understory plants by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_bleaf_understory_si_scpf )

    call this%set_history_var(vname='NPLANT_UNDERSTORY_SCPF', units = 'N/ha',         &
          long='stem number of understory plants density by pft/size', use_default='inactive', &
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_understory_si_scpf )

    call this%set_history_var(vname='CWD_AG_CWDSC', units='gC/m^2', &
          long='size-resolved AG CWD stocks', use_default='inactive', &
          avgflag='A', vtype=site_cwdsc_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_cwd_ag_si_cwdsc )

    call this%set_history_var(vname='CWD_BG_CWDSC', units='gC/m^2', &
          long='size-resolved BG CWD stocks', use_default='inactive', &
          avgflag='A', vtype=site_cwdsc_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_cwd_bg_si_cwdsc )

    call this%set_history_var(vname='CWD_AG_IN_CWDSC', units='gC/m^2/y', &
          long='size-resolved AG CWD input', use_default='inactive', &
          avgflag='A', vtype=site_cwdsc_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_cwd_ag_in_si_cwdsc )

    call this%set_history_var(vname='CWD_BG_IN_CWDSC', units='gC/m^2/y', &
          long='size-resolved BG CWD input', use_default='inactive', &
          avgflag='A', vtype=site_cwdsc_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_cwd_bg_in_si_cwdsc )

    call this%set_history_var(vname='CWD_AG_OUT_CWDSC', units='gC/m^2/y', &
          long='size-resolved AG CWD output', use_default='inactive', &
          avgflag='A', vtype=site_cwdsc_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_cwd_ag_out_si_cwdsc )

    call this%set_history_var(vname='CWD_BG_OUT_CWDSC', units='gC/m^2/y', &
          long='size-resolved BG CWD output', use_default='inactive', &
          avgflag='A', vtype=site_cwdsc_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_cwd_bg_out_si_cwdsc )

    ! Size structured diagnostics that require rapid updates (upfreq=2)

    call this%set_history_var(vname='AR_SCPF',units = 'kgC/m2/yr',          &
          long='total autotrophic respiration per m2 per year by pft/size',use_default='inactive',&
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_ar_si_scpf )
    
    call this%set_history_var(vname='AR_GROW_SCPF',units = 'kgC/m2/yr',          &
          long='growth autotrophic respiration per m2 per year by pft/size',use_default='inactive',&
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_ar_grow_si_scpf )

    call this%set_history_var(vname='AR_MAINT_SCPF',units = 'kgC/m2/yr',          &
          long='maintenance autotrophic respiration per m2 per year by pft/size',use_default='inactive',&
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_ar_maint_si_scpf )

    call this%set_history_var(vname='AR_DARKM_SCPF',units = 'kgC/m2/yr',          &
          long='dark portion of maintenance autotrophic respiration per m2 per year by pft/size',use_default='inactive',&
          avgflag='A', vtype=site_size_pft_r8,hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_ar_darkm_si_scpf )

    call this%set_history_var(vname='AR_AGSAPM_SCPF',units = 'kgC/m2/yr',          &
          long='above-ground sapwood maintenance autotrophic respiration per m2 per year by pft/size',use_default='inactive',&
          avgflag='A', vtype=site_size_pft_r8,hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_ar_agsapm_si_scpf )
    
    call this%set_history_var(vname='AR_CROOTM_SCPF',units = 'kgC/m2/yr',          &
          long='below-ground sapwood maintenance autotrophic respiration per m2 per year by pft/size',use_default='inactive',&
          avgflag='A', vtype=site_size_pft_r8,hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_ar_crootm_si_scpf )

    call this%set_history_var(vname='AR_FROOTM_SCPF',units = 'kgC/m2/yr',          &
          long='fine root maintenance autotrophic respiration per m2 per year by pft/size',use_default='inactive',&
          avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_ar_frootm_si_scpf )

    ! size-class only variables

    call this%set_history_var(vname='DDBH_CANOPY_SCLS', units = 'cm/yr/ha',         &
          long='diameter growth increment by pft/size',use_default='inactive', &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ddbh_canopy_si_scls )

    call this%set_history_var(vname='DDBH_UNDERSTORY_SCLS', units = 'cm/yr/ha',         &
          long='diameter growth increment by pft/size',use_default='inactive', &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ddbh_understory_si_scls )

    call this%set_history_var(vname='YESTERDAYCANLEV_CANOPY_SCLS', units = 'indiv/ha',               &
          long='Yesterdays canopy level for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_yesterdaycanopylevel_canopy_si_scls )

    call this%set_history_var(vname='YESTERDAYCANLEV_UNDERSTORY_SCLS', units = 'indiv/ha',               &
          long='Yesterdays canopy level for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_yesterdaycanopylevel_understory_si_scls )

    call this%set_history_var(vname='BA_SCLS', units = 'm2/ha',               &
          long='basal area by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_ba_si_scls )

    call this%set_history_var(vname='DEMOTION_RATE_SCLS', units = 'indiv/ha/yr',               &
          long='demotion rate from canopy to understory by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_demotion_rate_si_scls )

    call this%set_history_var(vname='PROMOTION_RATE_SCLS', units = 'indiv/ha/yr',               &
          long='promotion rate from understory to canopy by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_promotion_rate_si_scls )

    call this%set_history_var(vname='NPLANT_CANOPY_SCLS', units = 'indiv/ha',               &
          long='number of canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_canopy_si_scls )

    call this%set_history_var(vname='MORTALITY_CANOPY_SCLS', units = 'indiv/ha/yr',               &
          long='total mortality of canopy trees by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_mortality_canopy_si_scls )

    call this%set_history_var(vname='NPLANT_UNDERSTORY_SCLS', units = 'indiv/ha',               &
          long='number of understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_nplant_understory_si_scls )

    call this%set_history_var(vname='MORTALITY_UNDERSTORY_SCLS', units = 'indiv/ha/yr',               &
          long='total mortality of understory trees by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_mortality_understory_si_scls )

    call this%set_history_var(vname='TRIMMING_CANOPY_SCLS', units = 'indiv/ha',               &
          long='trimming term of canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_trimming_canopy_si_scls )

    call this%set_history_var(vname='TRIMMING_UNDERSTORY_SCLS', units = 'indiv/ha',               &
          long='trimming term of understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_trimming_understory_si_scls )

    call this%set_history_var(vname='CROWN_AREA_CANOPY_SCLS', units = 'm2/ha',               &
          long='total crown area of canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_crown_area_canopy_si_scls )

    call this%set_history_var(vname='CROWN_AREA_UNDERSTORY_SCLS', units = 'm2/ha',               &
          long='total crown area of understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_crown_area_understory_si_scls )

    call this%set_history_var(vname='LEAF_MD_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='LEAF_MD for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_leaf_md_canopy_si_scls )
    
    call this%set_history_var(vname='ROOT_MD_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='ROOT_MD for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_root_md_canopy_si_scls )
    
    call this%set_history_var(vname='CARBON_BALANCE_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='CARBON_BALANCE for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_carbon_balance_canopy_si_scls )
    
    call this%set_history_var(vname='SEED_PROD_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='SEED_PROD for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_seed_prod_canopy_si_scls )
    
    call this%set_history_var(vname='DBALIVEDT_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='DBALIVEDT for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_dbalivedt_canopy_si_scls )
    
    call this%set_history_var(vname='DBDEADDT_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='DBDEADDT for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_dbdeaddt_canopy_si_scls )
    
    call this%set_history_var(vname='DBSTOREDT_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='DBSTOREDT for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_dbstoredt_canopy_si_scls )
    
    call this%set_history_var(vname='STORAGE_FLUX_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='STORAGE_FLUX for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_storage_flux_canopy_si_scls )
    
    call this%set_history_var(vname='NPP_LEAF_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='NPP_LEAF for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_leaf_canopy_si_scls )
    
    call this%set_history_var(vname='NPP_FROOT_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='NPP_FROOT for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_froot_canopy_si_scls )
    
    call this%set_history_var(vname='NPP_BSW_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='NPP_BSW for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_bsw_canopy_si_scls )
    
    call this%set_history_var(vname='NPP_BDEAD_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='NPP_BDEAD for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_bdead_canopy_si_scls )
    
    call this%set_history_var(vname='NPP_BSEED_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='NPP_BSEED for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_bseed_canopy_si_scls )
    
    call this%set_history_var(vname='NPP_STORE_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='NPP_STORE for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_store_canopy_si_scls )
    
    call this%set_history_var(vname='RDARK_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='RDARK for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_rdark_canopy_si_scls )
    
    call this%set_history_var(vname='LIVESTEM_MR_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='LIVESTEM_MR for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_livestem_mr_canopy_si_scls )
    
    call this%set_history_var(vname='LIVECROOT_MR_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='LIVECROOT_MR for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_livecroot_mr_canopy_si_scls )
    
    call this%set_history_var(vname='FROOT_MR_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='FROOT_MR for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_froot_mr_canopy_si_scls )
    
    call this%set_history_var(vname='RESP_G_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='RESP_G for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_resp_g_canopy_si_scls )
    
    call this%set_history_var(vname='RESP_M_CANOPY_SCLS', units = 'kg C / ha / yr',               &
          long='RESP_M for canopy plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_resp_m_canopy_si_scls )

    call this%set_history_var(vname='LEAF_MD_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='LEAF_MD for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_leaf_md_understory_si_scls )
    
    call this%set_history_var(vname='ROOT_MD_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='ROOT_MD for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_root_md_understory_si_scls )
    
    call this%set_history_var(vname='CARBON_BALANCE_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='CARBON_BALANCE for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_carbon_balance_understory_si_scls )
    
    call this%set_history_var(vname='SEED_PROD_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='SEED_PROD for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_seed_prod_understory_si_scls )
    
    call this%set_history_var(vname='DBALIVEDT_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='DBALIVEDT for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_dbalivedt_understory_si_scls )
    
    call this%set_history_var(vname='DBDEADDT_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='DBDEADDT for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_dbdeaddt_understory_si_scls )
    
    call this%set_history_var(vname='DBSTOREDT_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='DBSTOREDT for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_dbstoredt_understory_si_scls )
    
    call this%set_history_var(vname='STORAGE_FLUX_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='STORAGE_FLUX for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_storage_flux_understory_si_scls )
    
    call this%set_history_var(vname='NPP_LEAF_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='NPP_LEAF for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_leaf_understory_si_scls )
    
    call this%set_history_var(vname='NPP_FROOT_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='NPP_FROOT for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_froot_understory_si_scls )
    
    call this%set_history_var(vname='NPP_BSW_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='NPP_BSW for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_bsw_understory_si_scls )
    
    call this%set_history_var(vname='NPP_BDEAD_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='NPP_BDEAD for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_bdead_understory_si_scls )
    
    call this%set_history_var(vname='NPP_BSEED_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='NPP_BSEED for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_bseed_understory_si_scls )
    
    call this%set_history_var(vname='NPP_STORE_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='NPP_STORE for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=1, ivar=ivar, initialize=initialize_variables, index = ih_npp_store_understory_si_scls )
    
    call this%set_history_var(vname='RDARK_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='RDARK for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_rdark_understory_si_scls )
    
    call this%set_history_var(vname='LIVESTEM_MR_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='LIVESTEM_MR for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_livestem_mr_understory_si_scls )
    
    call this%set_history_var(vname='LIVECROOT_MR_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='LIVECROOT_MR for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_livecroot_mr_understory_si_scls )
    
    call this%set_history_var(vname='FROOT_MR_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='FROOT_MR for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_froot_mr_understory_si_scls )
    
    call this%set_history_var(vname='RESP_G_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='RESP_G for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_resp_g_understory_si_scls )
    
    call this%set_history_var(vname='RESP_M_UNDERSTORY_SCLS', units = 'kg C / ha / yr',               &
          long='RESP_M for understory plants by size class', use_default='inactive',   &
          avgflag='A', vtype=site_size_r8, hlms='CLM:ALM', flushval=0.0_r8,    &
          upfreq=2, ivar=ivar, initialize=initialize_variables, index = ih_resp_m_understory_si_scls )


    ! CARBON BALANCE VARIABLES THAT DEPEND ON HLM BGC INPUTS

    call this%set_history_var(vname='NEP', units='gC/m^2/s', &
          long='net ecosystem production', use_default='active', &
          avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
          upfreq=3, ivar=ivar, initialize=initialize_variables, index = ih_nep_si )

    call this%set_history_var(vname='Fire_Closs', units='gC/m^2/s', &
          long='ED/SPitfire Carbon loss to atmosphere', use_default='active', &
          avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
          upfreq=3, ivar=ivar, initialize=initialize_variables, index = ih_fire_c_to_atm_si )
   
    call this%set_history_var(vname='NBP', units='gC/m^2/s', &
          long='net biosphere production', use_default='active', &
          avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
          upfreq=3, ivar=ivar, initialize=initialize_variables, index = ih_nbp_si )
   
    call this%set_history_var(vname='TOTECOSYSC', units='gC/m^2',  &
         long='total ecosystem carbon', use_default='active', &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
         upfreq=3, ivar=ivar, initialize=initialize_variables, index = ih_totecosysc_si )
    
    call this%set_history_var(vname='CBALANCE_ERROR_ED', units='gC/m^2/s',  &
         long='total carbon balance error on ED side', use_default='active', &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
         upfreq=3, ivar=ivar, initialize=initialize_variables, index = ih_cbal_err_fates_si )

    call this%set_history_var(vname='CBALANCE_ERROR_BGC', units='gC/m^2/s',  &
         long='total carbon balance error on HLMs BGC side', use_default='active', &
         avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
         upfreq=3, ivar=ivar, initialize=initialize_variables, index = ih_cbal_err_bgc_si )
    
    call this%set_history_var(vname='CBALANCE_ERROR_TOTAL', units='gC/m^2/s', &
          long='total carbon balance error total', use_default='active', &
          avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
          upfreq=3, ivar=ivar, initialize=initialize_variables, index = ih_cbal_err_tot_si )
    
    call this%set_history_var(vname='BIOMASS_STOCK_COL', units='gC/m^2',  &
          long='total ED biomass carbon at the column level', use_default='active', &
          avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
          upfreq=3, ivar=ivar, initialize=initialize_variables, index = ih_biomass_stock_si )
    
    call this%set_history_var(vname='ED_LITTER_STOCK_COL', units='gC/m^2', &
          long='total ED litter carbon at the column level', use_default='active', &
          avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
          upfreq=3, ivar=ivar, initialize=initialize_variables, index = ih_litter_stock_si )
    
    call this%set_history_var(vname='CWD_STOCK_COL', units='gC/m^2', &
          long='total CWD carbon at the column level', use_default='active', &
          avgflag='A', vtype=site_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
          upfreq=3, ivar=ivar, initialize=initialize_variables, index = ih_cwd_stock_si )
   

    ! PLANT HYDRAULICS

    if(hlm_use_planthydro.eq.itrue) then
       
       call this%set_history_var(vname='FATES_ERRH2O_SCPF', units='kg/indiv/s', &
             long='mean individual water balance error', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_errh2o_scpf )

       call this%set_history_var(vname='FATES_TRAN_SCPF', units='kg/indiv/s', &
             long='mean individual transpiration rate', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_tran_scpf )

       call this%set_history_var(vname='FATES_ROOTUPTAKE_SCPF', units='kg/indiv/s', &
             long='mean individual root uptake rate', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_rootuptake_scpf )

       call this%set_history_var(vname='FATES_ROOTUPTAKE01_SCPF', units='kg/indiv/s', &
             long='mean individual root uptake rate, layer 1', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_rootuptake01_scpf )
       
       call this%set_history_var(vname='FATES_ROOTUPTAKE02_SCPF', units='kg/indiv/s', &
             long='mean individual root uptake rate, layer 2', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_rootuptake02_scpf )
       
       call this%set_history_var(vname='FATES_ROOTUPTAKE03_SCPF', units='kg/indiv/s', &
             long='mean individual root uptake rate, layer 3', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_rootuptake03_scpf )
       
       call this%set_history_var(vname='FATES_ROOTUPTAKE04_SCPF', units='kg/indiv/s', &
             long='mean individual root uptake rate, layer 4', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_rootuptake04_scpf )
       
       call this%set_history_var(vname='FATES_ROOTUPTAKE05_SCPF', units='kg/indiv/s', &
             long='mean individual root uptake rate, layer 5', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_rootuptake05_scpf )
       
       call this%set_history_var(vname='FATES_ROOTUPTAKE06_SCPF', units='kg/indiv/s', &
             long='mean individual root uptake rate, layer 6', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_rootuptake06_scpf )
          
       call this%set_history_var(vname='FATES_ROOTUPTAKE07_SCPF', units='kg/indiv/s', &
             long='mean individual root uptake rate, layer 7', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_rootuptake07_scpf )
       
       call this%set_history_var(vname='FATES_ROOTUPTAKE08_SCPF', units='kg/indiv/s', &
             long='mean individual root uptake rate, layer 8', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_rootuptake08_scpf )
       
       call this%set_history_var(vname='FATES_ROOTUPTAKE09_SCPF', units='kg/indiv/s', &
             long='mean individual root uptake rate, layer 9', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_rootuptake09_scpf )
       
       call this%set_history_var(vname='FATES_ROOTUPTAKE10_SCPF', units='kg/indiv/s', &
             long='mean individual root uptake rate, layer 10', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_rootuptake10_scpf )

       call this%set_history_var(vname='FATES_SAPFLOW_COL_SCPF', units='kg/indiv/s', &
             long='individual sap flow rate', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_sapflow_scpf )
       
       call this%set_history_var(vname='FATES_ITERH1_COL_SCPF', units='count/indiv/step', &
             long='number of outer iterations required to achieve tolerable water balance error', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_iterh1_scpf )
       
       call this%set_history_var(vname='FATES_ITERH2_COL_SCPF', units='count/indiv/step', &
             long='number of inner iterations required to achieve tolerable water balance error', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_iterh2_scpf )
       
       call this%set_history_var(vname='FATES_ATH_COL_SCPF', units='m3 m-3', &
             long='absorbing root water content', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_ath_scpf )
       
       call this%set_history_var(vname='FATES_TTH_COL_SCPF', units='m3 m-3', &
             long='transporting root water content', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index =  ih_tth_scpf )
       
       call this%set_history_var(vname='FATES_STH_COL_SCPF', units='m3 m-3', &
             long='stem water contenet', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_sth_scpf )
       
       call this%set_history_var(vname='FATES_LTH_COL_SCPF', units='m3 m-3', &
             long='leaf water content', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_lth_scpf )

       call this%set_history_var(vname='FATES_AWP_COL_SCPF', units='MPa', &
             long='absorbing root water potential', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_awp_scpf )
       
       call this%set_history_var(vname='FATES_TWP_COL_SCPF', units='MPa', &
             long='transporting root water potential', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_twp_scpf )
       
       call this%set_history_var(vname='FATES_SWP_COL_SCPF', units='MPa', &
             long='stem water potential', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_swp_scpf )
       
       call this%set_history_var(vname='FATES_LWP_COL_SCPF', units='MPa', &
             long='leaf water potential', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_lwp_scpf )
       
       call this%set_history_var(vname='FATES_BTRAN_COL_SCPF', units='MPa', &
             long='mean individual level btran', use_default='active', &
             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_btran_scpf )
       
!       call this%set_history_var(vname='FATES_LAROOT_COL_SCPF', units='kg/indiv/s', &
!             long='Needs Description', use_default='active', &
!             avgflag='A', vtype=site_size_pft_r8, hlms='CLM:ALM', flushval=hlm_hio_ignore_val,    &
!             upfreq=4, ivar=ivar, initialize=initialize_variables, index = ih_laroot_scpf)

    end if

    ! Must be last thing before return
    this%num_history_vars_ = ivar
    
  end subroutine define_history_vars


   ! ====================================================================================
   ! DEPRECATED, TRANSITIONAL OR FUTURE CODE SECTION
   ! ====================================================================================

   !subroutine set_fates_hio_str(tag,iotype_name, iostr_val)

!       ! Arguments
!       character(len=*), intent(in)           :: tag
!       character(len=*), optional, intent(in) :: iotype_name
!       integer, optional, intent(in)         :: iostr_val

!       ! local variables
!       logical              :: all_set
!       integer,  parameter  :: unset_int = -999
!       real(r8), parameter  :: unset_double = -999.9
!       integer              :: ityp, idim

!       select case (trim(tag))
!       case('flush_to_unset')
!          write(*, *) ''
!          write(*, *) 'Flushing FATES IO types prior to transfer from host'
!          do ityp=1,ubound(iovar_str, 1)
!             iovar_str(ityp)%dimsize = unset_int
!             iovar_str(ityp)%active  = .false.
!          end do

!       case('check_allset')
!          do ityp=1,ubound(iovar_str, 1)
!             write(*, *) 'Checking to see if ',iovar_str(ityp)%name, ' IO communicators were sent to FATES'
!             if(iovar_str(ityp)%active)then
!                if(iovar_str(ityp)%offset .eq. unset_int) then
!                   write(*, *) 'FATES offset information of IO type:', iovar_str(ityp)%name
!                   write(*, *) 'was never set'
!                   ! end_run('MESSAGE')
!                end if
!                do idim=1, iovar_str(ityp)%ndims
!                   if(iovar_str(ityp)%dimsize(idim) .eq. unset_int) then
!                      write(*, *) 'FATES dimension information of IO type:', iovar_str(ityp)%name
!                      write(*, *) 'was never set'
!                      ! end_run('MESSAGE')
!                   end if
!                end do
!             end if
!          end do
!          write(*, *) 'Checked. All history IO specifications properly sent to FATES.'
!       case default

!          ! Must have two arguments if this is not a check or flush
!          if(present(iostr_val) .and. present(iotype_name))then
!
!             ! Tag in this case is dimsize or offset
!             select case (trim(tag))
!
!             case('offset')
!                ityp=iotype_index(trim(iotype_name))
!                iovar_str(ityp)%offset = iostr_val
!                write(*, *) 'Transfering offset for IOTYPE',iotype_name, ' to FATES'

!             case('dimsize1')
!                ityp=iotype_index(trim(iotype_name))
!                iovar_str(ityp)%dimsize(1) = iostr_val
!                write(*, *) 'Transfering 1st dimension size for IOTYPE',iotype_name, ' to FATES'

!             case('dimsize2')
!                ityp=iotype_index(trim(iotype_name))
!                if(ubound(iovar_str(ityp)%dimsize, 1)==1)then
!                   write(fates_log(), *) 'Transfering second dimensional bound to unallocated space'
!                   write(fates_log(), *) 'type:', iotype_name
!                   ! end_run
!                end if
!                iovar_str(ityp)%dimsize(2) = iostr_val
!                write(*, *) 'Transfering 2nd dimension size for IOTYPE',iotype_name, ' to FATES'

!             case('dimsize3')
!                ityp=iotype_index(trim(iotype_name))
!                if(ubound(iovar_str(ityp)%dimsize, 1)<3)then
!                   write(fates_log(), *) 'Transfering third dimensional bound to unallocated space'
!                   write(fates_log(), *) 'type:', iotype_name
!                   ! end_run
!                end if
!                iovar_str(ityp)%dimsize(3) = iostr_val
!                write(*, *) 'Transfering 3rd dimension size for IOTYPE',iotype_name, ' to FATES'

!             case default
!                write(*, *) 'IO parameter not recognized:', trim(tag)
!                ! end_run
!             end select
!          else
!             write(*, *) 'no value was provided for the tag'
!          end if
!
!       end select
!       return
!     end subroutine set_fates_hio_str



end module FatesHistoryInterfaceMod
