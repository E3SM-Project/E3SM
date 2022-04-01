module ColumnWorkRoutinesMod

   use shr_kind_mod   , only : r8 => shr_kind_r8
   use elm_varpar     , only : ndecomp_cascade_transitions, ndecomp_pools
   use elm_varpar     , only : nlevdecomp_full, nlevdecomp
   use elm_varcon     , only : spval
   use ColumnDataType , only : column_carbon_flux, column_carbon_state
   use ColumnDataType , only : column_nitrogen_flux, column_nitrogen_state
   use ColumnDataType , only : column_phosphorus_flux, column_phosphorus_state
   use decompMod      , only : bounds_type
   use CNDecompCascadeConType , only : decomp_cascade_con
   use elm_varcon      , only : watmin, bdsno, zsoi, zisoi, dzsoi_decomp
   use elm_varctl      , only : bound_h2osoi, use_cn, iulog, use_vertsoilc, spinup_state
   use elm_varpar      , only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
   use elm_varcon      , only : spval, ispval, zlnd, snw_rds_min, denice, denh2o, tfrz, pondmx
   use elm_varctl      , only : use_erosion
   use elm_varctl      , only : use_elm_interface, use_pflotran, pf_cmode
   use elm_varctl      , only : hist_wrtch4diag, use_century_decomp
   use elm_varctl      , only : get_carbontag, override_bgc_restart_mismatch_dump
   use elm_varctl      , only : pf_hmode, nu_com
   use ch4varcon       , only : allowlakeprod
   use pftvarcon       , only : VMAX_MINSURF_P_vr, KM_MINSURF_P_vr
   use soilorder_varcon, only : smax, ks_sorption


   implicit none
   public
   public :: col_cf_setvalues_acc 
   public :: col_nf_setvalues_acc
   public :: col_pf_setvalues_acc
   public :: col_cs_summary_acc
   public :: col_ns_summary_acc
   public :: col_ps_summary_acc



contains

   !-----------------------------------------------------------------------
   subroutine col_cf_setvalues_acc( this, num_column, filter_column)
     ! !DESCRIPTION:
     ! Set column-level carbon fluxes
     ! !ARGUMENTS:
     type (column_carbon_flux) :: this
     integer , intent(in) :: num_column
     integer , intent(in) :: filter_column(:)
     real(r8), parameter :: value_column=0.0_r8
     !
     ! !LOCAL VARIABLES:
     integer :: fi,i,j,k,l     ! loop index
     !------------------------------------------------------------------------

     !$acc parallel loop independent gang vector default(present) collapse(2) !async(1)
     do j = 1, nlevdecomp_full
        do fi = 1,num_column
           i = filter_column(fi)

           this%phenology_c_to_litr_met_c(i,j)     = value_column
           this%phenology_c_to_litr_cel_c(i,j)     = value_column
           this%phenology_c_to_litr_lig_c(i,j)     = value_column

           this%gap_mortality_c_to_litr_met_c(i,j) = value_column
           this%gap_mortality_c_to_litr_cel_c(i,j) = value_column
           this%gap_mortality_c_to_litr_lig_c(i,j) = value_column
           this%gap_mortality_c_to_cwdc(i,j)       = value_column

           this%fire_mortality_c_to_cwdc(i,j)      = value_column
           this%m_c_to_litr_met_fire(i,j)          = value_column
           this%m_c_to_litr_cel_fire(i,j)          = value_column
           this%m_c_to_litr_lig_fire(i,j)          = value_column

           this%harvest_c_to_litr_met_c(i,j)       = value_column
           this%harvest_c_to_litr_cel_c(i,j)       = value_column
           this%harvest_c_to_litr_lig_c(i,j)       = value_column
           this%harvest_c_to_cwdc(i,j)             = value_column

           this%hr_vr(i,j)                         = value_column
        end do
     end do

     !$acc parallel loop independent gang worker default(present) collapse(2) !async(2)
     do k = 1, ndecomp_pools
        do j = 1, nlevdecomp_full
           !$acc loop vector independent
           do fi = 1,num_column
              i = filter_column(fi)
              this%m_decomp_cpools_to_fire_vr(i,j,k) = value_column
              this%decomp_cpools_transport_tendency(i,j,k) = value_column
              this%decomp_cpools_yield_vr(i,j,k) = value_column
              this%bgc_cpool_ext_inputs_vr(i,j, k) = value_column
              this%bgc_cpool_ext_loss_vr(i,j, k) = value_column
              this%decomp_cpools_sourcesink(i,j,k) = value_column

           end do
        end do
     end do

     !$acc parallel loop independent gang vector default(present) collapse(2) !async(3)
     do l = 1, ndecomp_cascade_transitions
        do fi = 1,num_column
           i = filter_column(fi)
           this%decomp_cascade_hr(i,l) = value_column
           this%decomp_cascade_ctransfer(i,l) = value_column
        end do
     end do

     !$acc parallel loop independent gang worker default(present) collapse(2) !async(4)
     do l = 1, ndecomp_cascade_transitions
        do j = 1, nlevdecomp_full
           !$acc loop vector independent
           do fi = 1,num_column
              i = filter_column(fi)
              this%decomp_cascade_hr_vr(i,j,l) = value_column
              this%decomp_cascade_ctransfer_vr(i,j,l) = value_column
              this%decomp_k(i,j,l) = value_column
           end do
        end do
     end do

     !$acc parallel loop independent gang vector default(present) collapse(2) !async(5)
     do k = 1, ndecomp_pools
        do fi = 1,num_column
           i = filter_column(fi)
           this%decomp_cpools_leached(i,k) = value_column
           this%decomp_cpools_erode(i,k) = value_column
           this%decomp_cpools_deposit(i,k) = value_column
           this%decomp_cpools_yield(i,k) = value_column
           this%m_decomp_cpools_to_fire(i,k) = value_column
        end do
     end do

     !$acc parallel loop independent gang vector default(present) !async(6)
     do fi = 1,num_column
        i = filter_column(fi)

        this%hrv_deadstemc_to_prod10c(i)  = value_column
        this%hrv_deadstemc_to_prod100c(i) = value_column
        this%hrv_cropc_to_prod1c(i)       = value_column
        this%somc_fire(i)                 = value_column
        this%prod1c_loss(i)               = value_column
        this%prod10c_loss(i)              = value_column
        this%prod100c_loss(i)             = value_column
        this%product_closs(i)             = value_column
        this%somhr(i)                     = value_column
        this%lithr(i)                     = value_column
        this%hr(i)                        = value_column
        this%sr(i)                        = value_column
        this%er(i)                        = value_column
        this%litfire(i)                   = value_column
        this%somfire(i)                   = value_column
        this%totfire(i)                   = value_column
        this%nep(i)                       = value_column
        this%nbp(i)                       = value_column
        this%nee(i)                       = value_column
        this%cinputs(i)                   = value_column
        this%coutputs(i)                  = value_column
        this%fire_closs(i)                = value_column
        this%cwdc_hr(i)                   = value_column
        this%cwdc_loss(i)                 = value_column
        this%litterc_loss(i)              = value_column
        this%som_c_leached(i)             = value_column
        this%somc_erode(i)                = value_column
        this%somc_deposit(i)              = value_column
        this%somc_yield(i)                = value_column

        ! Zero p2c column fluxes
        this%rr(i)                    = value_column
        this%ar(i)                    = value_column
        this%gpp(i)                   = value_column
        this%npp(i)                   = value_column
        this%fire_closs(i)            = value_column
        this%litfall(i)               = value_column
        this%vegfire(i)               = value_column
        this%wood_harvestc(i)         = value_column
        this%hrv_xsmrpool_to_atm(i)   = value_column
     end do

     ! do k = 1, ndecomp_pools
     !    do j = 1, nlevdecomp_full
     !       do fi = 1,num_column
     !          i = filter_column(fi)
     !          this%decomp_cpools_sourcesink(i,j,k) = value_column
     !       end do
     !    end do
     ! end do

     !$acc parallel loop independent gang vector collapse(2) default(present) !async(7)
     do j = 1, nlevdecomp_full
        do fi = 1,num_column
           i = filter_column(fi)
           this%f_co2_soil_vr(i,j) = value_column
        end do
     end do
     ! !$acc wait
  end subroutine col_cf_setvalues_acc

   !-----------------------------------------------------------------------
   subroutine col_nf_setvalues_acc( this, num_column, filter_column)
     !
     ! !DESCRIPTION:
     ! Set column-level nitrogen fluxes
     ! !ARGUMENTS:
     type (column_nitrogen_flux)  :: this
     integer , intent(in)         :: num_column
     integer , intent(in)         :: filter_column(:)
     real(r8), parameter         :: value_column=0.0_r8
     !
     ! !LOCAL VARIABLES:
     integer :: fi,i,j,k,l     ! loop index
     !------------------------------------------------------------------------

     !$acc parallel loop independent gang vector collapse(2) default(present) !async(1)
     do j = 1, nlevdecomp_full
        do fi = 1,num_column
           i = filter_column(fi)

           ! phenology: litterfall and crop fluxes associated wit
           this%phenology_n_to_litr_met_n(i,j)        = value_column
           this%phenology_n_to_litr_cel_n(i,j)        = value_column
           this%phenology_n_to_litr_lig_n(i,j)        = value_column

           ! gap mortality
           this%gap_mortality_n_to_litr_met_n(i,j)    = value_column
           this%gap_mortality_n_to_litr_cel_n(i,j)    = value_column
           this%gap_mortality_n_to_litr_lig_n(i,j)    = value_column
           this%gap_mortality_n_to_cwdn(i,j)          = value_column

           ! fire
           this%fire_mortality_n_to_cwdn(i,j)         = value_column
           this%m_n_to_litr_met_fire(i,j)             = value_column
           this%m_n_to_litr_cel_fire(i,j)             = value_column
           this%m_n_to_litr_lig_fire(i,j)             = value_column

           ! harvest
           this%harvest_n_to_litr_met_n(i,j)          = value_column
           this%harvest_n_to_litr_cel_n(i,j)          = value_column
           this%harvest_n_to_litr_lig_n(i,j)          = value_column
           this%harvest_n_to_cwdn(i,j)                = value_column

           this%f_nit_vr(i,j)                      = value_column
           this%f_denit_vr(i,j)                    = value_column
           this%smin_no3_leached_vr(i,j)           = value_column
           this%smin_no3_runoff_vr(i,j)            = value_column
           this%n2_n2o_ratio_denit_vr(i,j)         = value_column
           this%pot_f_nit_vr(i,j)                  = value_column
           this%pot_f_denit_vr(i,j)                = value_column
           this%actual_immob_no3_vr(i,j)           = value_column
           this%actual_immob_nh4_vr(i,j)           = value_column
           this%smin_no3_to_plant_vr(i,j)          = value_column
           this%smin_nh4_to_plant_vr(i,j)          = value_column
           this%f_n2o_denit_vr(i,j)                = value_column
           this%f_n2o_nit_vr(i,j)                  = value_column

           this%smin_no3_massdens_vr(i,j)          = value_column
           this%k_nitr_t_vr(i,j)                   = value_column
           this%k_nitr_ph_vr(i,j)                  = value_column
           this%k_nitr_h2o_vr(i,j)                 = value_column
           this%k_nitr_vr(i,j)                     = value_column
           this%wfps_vr(i,j)                       = value_column
           this%fmax_denit_carbonsubstrate_vr(i,j) = value_column
           this%fmax_denit_nitrate_vr(i,j)         = value_column
           this%f_denit_base_vr(i,j)               = value_column

           this%diffus(i,j)                        = value_column
           this%ratio_k1(i,j)                      = value_column
           this%ratio_no3_co2(i,j)                 = value_column
           this%soil_co2_prod(i,j)                 = value_column
           this%fr_WFPS(i,j)                       = value_column
           this%soil_bulkdensity(i,j)              = value_column

           this%r_psi(i,j)                         = value_column
           this%anaerobic_frac(i,j)                = value_column

           ! pflotran
           this%plant_ndemand_vr(i,j)              = value_column
           this%f_ngas_decomp_vr(i,j)              = value_column
           this%f_ngas_nitri_vr(i,j)               = value_column
           this%f_ngas_denit_vr(i,j)               = value_column
           this%f_n2o_soil_vr(i,j)                 = value_column
           this%f_n2_soil_vr(i,j)                  = value_column

           this%potential_immob_vr(i,j)               = value_column
           this%actual_immob_vr(i,j)                  = value_column
           this%sminn_to_plant_vr(i,j)                = value_column
           this%supplement_to_sminn_vr(i,j)           = value_column
           this%gross_nmin_vr(i,j)                    = value_column
           this%net_nmin_vr(i,j)                      = value_column
           this%sminn_nh4_input_vr(i,j)               = value_column
           this%sminn_no3_input_vr(i,j)               = value_column
        end do
     end do

     !$acc parallel loop independent gang vector default(present) !async(2)
     do fi = 1,num_column
        i = filter_column(fi)

        this%ndep_to_sminn(i)             = value_column
        this%nfix_to_sminn(i)             = value_column
        this%nfix_to_ecosysn(i)           = value_column
        this%fert_to_sminn(i)             = value_column
        this%soyfixn_to_sminn(i)          = value_column
        this%hrv_deadstemn_to_prod10n(i)  = value_column
        this%hrv_deadstemn_to_prod100n(i) = value_column
        this%hrv_cropn_to_prod1n(i)       = value_column
        this%prod10n_loss(i)              = value_column
        this%prod100n_loss(i)             = value_column
        this%prod1n_loss(i)               = value_column
        this%product_nloss(i)             = value_column
        this%potential_immob(i)           = value_column
        this%actual_immob(i)              = value_column
        this%sminn_to_plant(i)            = value_column
        this%supplement_to_sminn(i)       = value_column
        this%gross_nmin(i)                = value_column
        this%net_nmin(i)                  = value_column
        this%denit(i)                     = value_column

        this%f_nit(i)                  = value_column
        this%pot_f_nit(i)              = value_column
        this%f_denit(i)                = value_column
        this%pot_f_denit(i)            = value_column
        this%f_n2o_denit(i)            = value_column
        this%f_n2o_nit(i)              = value_column
        this%smin_no3_leached(i)       = value_column
        this%smin_no3_runoff(i)        = value_column

        this%f_ngas_decomp(i)         = value_column
        this%f_ngas_nitri(i)          = value_column
        this%f_ngas_denit(i)          = value_column
        this%f_n2o_soil(i)            = value_column
        this%f_n2_soil(i)             = value_column

        this%smin_nh4_to_plant(i)      = value_column
        this%smin_no3_to_plant(i)      = value_column

        this%ninputs(i)                   = value_column
        this%noutputs(i)                  = value_column
        this%fire_nloss(i)                = value_column
        this%som_n_leached(i)             = value_column
        this%sminn_input(i)               = value_column
        this%sminn_nh4_input(i)           = value_column
        this%sminn_no3_input(i)           = value_column
        this%somn_erode(i)                = value_column
        this%somn_deposit(i)              = value_column
        this%somn_yield(i)                = value_column
        ! Zero p2c column fluxes
        this%fire_nloss(i) = value_column
        this%wood_harvestn(i) = value_column

        ! bgc-interface
        this%plant_ndemand(i) = value_column
     end do

     !$acc parallel loop independent gang vector collapse(2) default(present) !async(3)
     do k = 1, ndecomp_pools
        do fi = 1,num_column
           i = filter_column(fi)
           this%decomp_npools_leached(i,k) = value_column
           this%decomp_npools_erode(i,k) = value_column
           this%decomp_npools_deposit(i,k) = value_column
           this%decomp_npools_yield(i,k) = value_column
           this%m_decomp_npools_to_fire(i,k) = value_column
           this%bgc_npool_inputs(i,k) = value_column
        end do
     end do

     !$acc parallel loop independent gang worker collapse(2) default(present) !async(4)
     do k = 1, ndecomp_pools
        do j = 1, nlevdecomp_full
           !$acc loop vector independent
           do fi = 1,num_column
              i = filter_column(fi)
              this%m_decomp_npools_to_fire_vr(i,j,k) = value_column
              this%decomp_npools_transport_tendency(i,j,k) = value_column
              this%decomp_npools_yield_vr(i,j,k) = value_column
              this%bgc_npool_ext_inputs_vr(i,j,k) = value_column
              this%bgc_npool_ext_loss_vr(i,j,k) = value_column
              this%decomp_npools_sourcesink(i,j,k) = value_column
           end do
        end do
     end do

     !$acc parallel loop independent gang vector collapse(2) default(present) !async(5)
     do l = 1, ndecomp_cascade_transitions
      do fi = 1,num_column
         i = filter_column(fi)
         this%decomp_cascade_ntransfer(i,l) = value_column
         this%decomp_cascade_sminn_flux(i,l) = value_column
      end do
     end do

      !$acc parallel loop independent gang worker collapse(2) default(present) !async(6)
      do l = 1, ndecomp_cascade_transitions
       do j = 1, nlevdecomp_full
          !$acc loop vector independent
          do fi = 1,num_column
             i = filter_column(fi)
             this%decomp_cascade_ntransfer_vr(i,j,l) = value_column
             this%decomp_cascade_sminn_flux_vr(i,j,l) = value_column
          end do
       end do
      end do

     ! do k = 1, ndecomp_pools
     !    do j = 1, nlevdecomp_full
     !       do fi = 1,num_column
     !          i = filter_column(fi)
     !          this%decomp_npools_sourcesink(i,j,k) = value_column
     !       end do
     !    end do
     ! end do
     ! !$acc wait

  end subroutine col_nf_setvalues_acc

   !-----------------------------------------------------------------------
   subroutine col_pf_setvalues_acc( this, num_column, filter_column)
     !
     ! !DESCRIPTION:
     ! Set phosphorus flux variables
     ! !ARGUMENTS:
     type (column_phosphorus_flux) :: this
     integer , intent(in) :: num_column
     integer , intent(in) :: filter_column(:)
     real(r8), parameter :: value_column=0.0_r8
     !
     ! !LOCAL VARIABLES:
     integer :: fi,i,j,k,l     ! loop index
     !------------------------------------------------------------------------
     !$acc parallel loop independent gang vector collapse(2) default(present) !async(1)
     do j = 1, nlevdecomp_full
        do fi = 1,num_column
           i = filter_column(fi)

           ! phenology: litterfall and crop fluxes associated wit
           this%phenology_p_to_litr_met_p(i,j)        = value_column
           this%phenology_p_to_litr_cel_p(i,j)        = value_column
           this%phenology_p_to_litr_lig_p(i,j)        = value_column

           ! gap mortality
           this%gap_mortality_p_to_litr_met_p(i,j)    = value_column
           this%gap_mortality_p_to_litr_cel_p(i,j)    = value_column
           this%gap_mortality_p_to_litr_lig_p(i,j)    = value_column
           this%gap_mortality_p_to_cwdp(i,j)          = value_column

           ! fire
           this%fire_mortality_p_to_cwdp(i,j)         = value_column
           this%m_p_to_litr_met_fire(i,j)             = value_column
           this%m_p_to_litr_cel_fire(i,j)             = value_column
           this%m_p_to_litr_lig_fire(i,j)             = value_column

           ! harvest
           this%harvest_p_to_litr_met_p(i,j)          = value_column
           this%harvest_p_to_litr_cel_p(i,j)          = value_column
           this%harvest_p_to_litr_lig_p(i,j)          = value_column
           this%harvest_p_to_cwdp(i,j)                = value_column

           this%primp_to_labilep_vr(i,j)              = value_column
           this%labilep_to_secondp_vr(i,j)            = value_column
           this%secondp_to_labilep_vr(i,j)            = value_column
           this%secondp_to_occlp_vr(i,j)              = value_column

           this%sminp_leached_vr(i,j)                 = value_column

           this%labilep_yield_vr(i,j)                 = value_column
           this%secondp_yield_vr(i,j)                 = value_column
           this%occlp_yield_vr(i,j)                   = value_column
           this%primp_yield_vr(i,j)                   = value_column

           this%potential_immob_p_vr(i,j)             = value_column
           this%actual_immob_p_vr(i,j)                = value_column
           this%sminp_to_plant_vr(i,j)                = value_column
           this%supplement_to_sminp_vr(i,j)           = value_column
           this%gross_pmin_vr(i,j)                    = value_column
           this%net_pmin_vr(i,j)                      = value_column
           this%biochem_pmin_vr(i,j)                  = value_column
           this%biochem_pmin_to_ecosysp_vr(i,j)       = value_column

           ! bgc interface & pflotran
           this%plant_pdemand_vr(i,j)                 = value_column
           this%adsorb_to_labilep_vr(i,j)             = value_column
           this%desorb_to_solutionp_vr(i,j)           = value_column
        end do
     end do

     !$acc parallel loop independent gang vector default(present) !async(2)
     do fi = 1,num_column
        i = filter_column(fi)

        this%pdep_to_sminp(i)             = value_column
        this%fert_p_to_sminp(i)           = value_column
        this%hrv_deadstemp_to_prod10p(i)  = value_column
        this%hrv_deadstemp_to_prod100p(i) = value_column
        this%hrv_cropp_to_prod1p(i)       = value_column
        this%prod10p_loss(i)              = value_column
        this%prod100p_loss(i)             = value_column
        this%product_ploss(i)             = value_column
        this%prod1p_loss(i)               = value_column
        this%potential_immob_p(i)         = value_column
        this%actual_immob_p(i)            = value_column
        this%sminp_to_plant(i)            = value_column
        this%supplement_to_sminp(i)       = value_column
        this%gross_pmin(i)                = value_column
        this%net_pmin(i)                  = value_column
        this%biochem_pmin(i)              = value_column
        this%biochem_pmin_to_plant(i)     = value_column
        this%primp_to_labilep(i)          = value_column
        this%labilep_to_secondp(i)        = value_column
        this%secondp_to_labilep(i)        = value_column
        this%secondp_to_occlp(i)          = value_column
        this%sminp_leached(i)             = value_column
        this%fire_ploss(i)                = value_column
        this%pinputs(i)                   = value_column
        this%poutputs(i)                  = value_column
        this%som_p_leached(i)             = value_column
        this%somp_erode(i)                = value_column
        this%somp_deposit(i)              = value_column
        this%somp_yield(i)                = value_column
        this%labilep_erode(i)             = value_column
        this%labilep_deposit(i)           = value_column
        this%labilep_yield(i)             = value_column
        this%secondp_erode(i)             = value_column
        this%secondp_deposit(i)           = value_column
        this%secondp_yield(i)             = value_column
        this%occlp_erode(i)               = value_column
        this%occlp_deposit(i)             = value_column
        this%occlp_yield(i)               = value_column
        this%primp_erode(i)               = value_column
        this%primp_deposit(i)             = value_column
        this%primp_yield(i)               = value_column

        ! Zero p2c column fluxes
        this%fire_ploss(i)                = value_column
        this%wood_harvestp(i)             = value_column

        ! bgc-interface
        this%plant_pdemand(i)             = value_column

        this%fire_ploss(i)                = value_column
        this%wood_harvestp(i)             = value_column

        this%adsorb_to_labilep(i)         = value_column
        this%desorb_to_solutionp(i)       = value_column

     end do

     !$acc parallel loop independent gang vector collapse(2) default(present) !async(3)
     do k = 1, ndecomp_pools
        do fi = 1,num_column
           i = filter_column(fi)
           this%decomp_ppools_leached(i,k) = value_column
           this%decomp_ppools_erode(i,k) = value_column
           this%decomp_ppools_deposit(i,k) = value_column
           this%decomp_ppools_yield(i,k) = value_column
           this%m_decomp_ppools_to_fire(i,k) = value_column
        end do
     end do

     !$acc parallel loop independent gang worker collapse(2) default(present) !async(4)
     do k = 1, ndecomp_pools
        do j = 1, nlevdecomp_full
           !$acc loop vector independent
           do fi = 1,num_column
              i = filter_column(fi)
              this%m_decomp_ppools_to_fire_vr(i,j,k) = value_column
              this%decomp_ppools_transport_tendency(i,j,k) = value_column
              this%decomp_ppools_yield_vr(i,j,k) = value_column
              this%decomp_ppools_sourcesink(i,j,k) = value_column
              this%biochem_pmin_ppools_vr(i,j,k) = value_column

           end do
        end do
     end do

     !$acc parallel loop independent gang vector collapse(2) default(present) !async(5)
     do l = 1, ndecomp_cascade_transitions
        do fi = 1,num_column
           i = filter_column(fi)
           this%decomp_cascade_ptransfer(i,l) = value_column
           this%decomp_cascade_sminp_flux(i,l) = value_column
        end do
     end do

     !$acc parallel loop independent gang worker collapse(2) default(present) !async(6)
     do l = 1, ndecomp_cascade_transitions
        do j = 1, nlevdecomp_full
           !$acc loop vector independent
           do fi = 1,num_column
              i = filter_column(fi)
              this%decomp_cascade_ptransfer_vr(i,j,l) = value_column
              this%decomp_cascade_sminp_flux_vr(i,j,l) = value_column
           end do
        end do
     end do

     ! do k = 1, ndecomp_pools
     !    do j = 1, nlevdecomp_full
     !       do fi = 1,num_column
     !          i = filter_column(fi)
     !          this%decomp_ppools_sourcesink(i,j,k) = value_column
     !          this%biochem_pmin_ppools_vr(i,j,k) = value_column
     !       end do
     !    end do
     ! end do
     ! !$acc wait

  end subroutine col_pf_setvalues_acc

  !-----------------------------------------------------------------------
  subroutine col_cs_summary_acc(this, bounds, num_soilc, filter_soilc)
   !
   ! !DESCRIPTION:
   ! Column-level carbon state summary calculations
   ! !ARGUMENTS:
   type(column_carbon_state) :: this
   type(bounds_type)      , intent(in)    :: bounds
   integer                , intent(in)    :: num_soilc       ! number of soil columns in filter
   integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
   !
   ! !LOCAL VARIABLES:
   real(r8) :: nfixlags, dtime ! temp variables for making lagged npp
   integer  :: c,p,j,k,l       ! indices
   integer  :: fp,fc           ! lake filter indices
   real(r8), parameter :: maxdepth = 1._r8  ! depth to integrate soil variables
   integer  :: nlev
   real(r8) :: sum1, sum2, sum3
   !-----------------------------------------------------------------------

   nlev = nlevdecomp
   if (use_pflotran .and. pf_cmode) nlev = nlevdecomp_full
   !$acc enter data create(sum1, sum2, sum3) 

   ! vertically integrate each of the decomposing C pools
   !$acc parallel loop independent gang vector collapse(2) default(present)
   do l = 1, ndecomp_pools
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         this%decomp_cpools(c,l) = 0._r8
      end do
   end do
   
   !$acc parallel loop independent gang worker collapse(2) default(present) private(sum1) 
   do l = 1, ndecomp_pools
      do fc = 1,num_soilc
         sum1 = 0.0_r8 
         c = filter_soilc(fc)
         !$acc loop vector reduction(+:sum1)
         do j = 1, nlev
            sum1 = sum1 + this%decomp_cpools_vr(c,j,l) * dzsoi_decomp(j)
         end do
         this%decomp_cpools(c,l) = this%decomp_cpools(c,l) + sum1 
      end do
   end do

   if ( nlevdecomp > 1) then
      ! vertically integrate each of the decomposing C pools to 1 meter
      ! maxdepth = 1._r8
      ! !$acc parallel loop independent collapse(2) gang vector default(present)
      ! do l = 1, ndecomp_pools
      !    do fc = 1,num_soilc
      !       c = filter_soilc(fc)
      !       this%decomp_cpools_1m(c,l) = 0._r8
      !    end do
      ! end do

      !$acc parallel loop independent gang worker collapse(2) default(present) private(sum1)
      do l = 1, ndecomp_pools
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            sum1 = 0.0_r8 
            !$acc loop vector reduction(+:sum1) 
            do j = 1, nlevdecomp
               if ( zisoi(j) <= maxdepth ) then
                  sum1 = this%decomp_cpools_vr(c,j,l) * dzsoi_decomp(j)
               elseif ( zisoi(j-1) < maxdepth ) then
                  sum1 = sum1 + this%decomp_cpools_vr(c,j,l) * (maxdepth - zisoi(j-1))
               endif
            end do
            this%decomp_cpools_1m(c,l) = sum1 
         end do
      end do

      ! total litter carbon in the top meter (TOTLITC_1m)
      ! do fc = 1,num_soilc
      !    c = filter_soilc(fc)
      !    this%totlitc_1m(c)         = 0._r8
      ! end do
      !$acc parallel loop independent gang worker default(present) private(sum1,sum2)
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         sum1 = 0.0_r8 
         sum2 = 0.0_r8 
         !$acc loop vector reduction(+:sum1,sum2)
         do l = 1, ndecomp_pools
            ! total litter carbon in the top meter (TOTLITC_1m)
            if ( decomp_cascade_con%is_litter(l) ) then
               sum1 = sum1 + this%decomp_cpools_1m(c,l)
            endif
            ! total soil organic matter carbon in the top meter (TOTSOMC_1m)
            if(decomp_cascade_con%is_soil(l) ) then 
               sum2 = sum2 + this%decomp_cpools_1m(c,l)
            end if 

         end do
         this%totlitc_1m(c) = sum1
         this%totsomc_1m(c) = sum2 
      end do

      ! ! total soil organic matter carbon in the top meter (TOTSOMC_1m)
      ! do fc = 1,num_soilc
      !    c = filter_soilc(fc)
      !    this%totsomc_1m(c) = 0._r8
      ! end do
      ! do l = 1, ndecomp_pools
      !    if ( decomp_cascade_con%is_soil(l) ) then
      !       do fc = 1,num_soilc
      !          c = filter_soilc(fc)
      !          this%totsomc_1m(c) = &
      !               this%totsomc_1m(c) + &
      !               this%decomp_cpools_1m(c,l)
      !       end do
      !    end if
      ! end do


   end if ! nlevdecomp>1

   ! total litter carbon (TOTLITC)
   ! do fc = 1,num_soilc
   !    c = filter_soilc(fc)
   !    this%totlitc(c) = 0._r8
   ! end do
   !$acc parallel loop independent gang worker default(present) private(sum1,sum2,sum3)
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      sum1 = 0.0_r8 
      sum2 = 0.0_r8 
      sum3 = 0.0_r8
      !$acc loop vector reduction(+:sum1,sum2,sum3) 
      do l = 1, ndecomp_pools
         ! total litter carbon (TOTLITC)
         if ( decomp_cascade_con%is_litter(l) ) then
            sum1 = sum1 + this%decomp_cpools(c,l)
         endif
         ! total soil organic matter carbon (TOTSOMC)
         if ( decomp_cascade_con%is_soil(l) ) then
            sum2 = sum2 + this%decomp_cpools(c,l)
         endif
         ! coarse woody debris carbon
         if ( decomp_cascade_con%is_cwd(l) ) then
            sum3 = sum3 + this%decomp_cpools(c,l)
         end if 
      end do
      this%totlitc(c) = sum1
      this%totsomc(c) = sum2 
      this%cwdc(c)    = sum3
   end do

   ! ! total soil organic matter carbon (TOTSOMC)
   ! do fc = 1,num_soilc
   !    c = filter_soilc(fc)
   !    this%totsomc(c) = 0._r8
   ! end do
   ! do l = 1, ndecomp_pools
   !    if ( decomp_cascade_con%is_soil(l) ) then
   !       do fc = 1,num_soilc
   !          c = filter_soilc(fc)
   !          this%totsomc(c) = &
   !               this%totsomc(c) + &
   !               this%decomp_cpools(c,l)
   !       end do
   !    end if
   ! end do

   ! ! coarse woody debris carbon
   ! do fc = 1,num_soilc
   !    c = filter_soilc(fc)
   !    this%cwdc(c) = 0._r8
   ! end do
   ! do l = 1, ndecomp_pools
   !    if ( decomp_cascade_con%is_cwd(l) ) then
   !       do fc = 1,num_soilc
   !          c = filter_soilc(fc)
   !          this%cwdc(c) = &
   !               this%cwdc(c) + &
   !               this%decomp_cpools(c,l)
   !       end do
   !    end if
   ! end do

   ! truncation carbon
   ! do fc = 1,num_soilc
   !    c = filter_soilc(fc)
   !    this%ctrunc(c) = 0._r8
   ! end do
   !$acc parallel loop independent gang worker default(present) private(sum1)
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      sum1 = 0.0_r8 
      !$acc loop vector reduction(+:sum1) 
      do j = 1, nlev
         sum1 = sum1 + this%ctrunc_vr(c,j) * dzsoi_decomp(j)
      end do
      this%ctrunc(c) = sum1 
   end do

   !$acc parallel loop independent gang vector default(present) 
   do fc = 1,num_soilc
      c = filter_soilc(fc)

      ! total product carbon
      this%totprodc(c) =      &
           this%prod10c(c)  + &
           this%prod100c(c) + &
           this%prod1c(c)

      ! total ecosystem carbon, including veg but excluding cpool (TOTECOSYSC)
      this%totecosysc(c) =    &
           this%cwdc(c)     + &
           this%totlitc(c)  + &
           this%totsomc(c)  + &
           this%totprodc(c) + &
           this%totvegc(c)

      ! total column carbon, including veg and cpool (TOTCOLC)
      ! adding col_ctrunc, seedc
      this%totcolc(c) =       &
           this%totpftc(c)  + &
           this%cwdc(c)     + &
           this%totlitc(c)  + &
           this%totsomc(c)  + &
           this%totprodc(c) + &
           this%ctrunc(c)   + &
           this%cropseedc_deficit(c)

      this%totabgc(c) =       &
           this%totprodc(c) + &
           this%seedc(c)    + &
           this%ctrunc(c)   + &
           this%totpftc(c)
   end do
   !$acc exit data delete(sum1, sum2, sum3) 

 end subroutine col_cs_summary_acc

 !-----------------------------------------------------------------------
 subroutine col_ns_summary_acc(this, bounds, num_soilc, filter_soilc)
   ! !ARGUMENTS:
   type (column_nitrogen_state)  :: this
   type(bounds_type) , intent(in) :: bounds
   integer           , intent(in) :: num_soilc       ! number of soil columns in filter
   integer           , intent(in) :: filter_soilc(:) ! filter for soil columns
   !
   ! !LOCAL VARIABLES:
   integer  :: c,p,j,k,l   ! indices
   integer  :: fp,fc       ! lake filter indices
   real(r8), parameter :: maxdepth = 1._r8   ! depth to integrate soil variables
   integer  :: nlev
   real(r8) :: sum1,sum2,sum3 
   !-----------------------------------------------------------------------

   ! vertically integrate NO3 NH4 N2O pools
   nlev = nlevdecomp
   if (use_pflotran .and. pf_cmode) nlev = nlevdecomp_full
   !$acc enter data create(sum1,sum2,sum3)
   
   ! do fc = 1,num_soilc
   !    c = filter_soilc(fc)
   !    this%smin_no3(c) = 0._r8
   !    this%smin_nh4(c) = 0._r8
   !    if(use_pflotran .and. pf_cmode) then
   !       this%smin_nh4sorb(c) = 0._r8
   !    end if
   ! end do

   !$acc parallel loop gang worker independent default(present) private(sum1,sum2,sum3)
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      !$acc loop vector reduction(+:sum1,sum2,sum3)
      do j = 1, nlev
         sum1 = sum1 + this%smin_no3_vr(c,j) * dzsoi_decomp(j)

         sum2 = sum2 + this%smin_nh4_vr(c,j) * dzsoi_decomp(j)
         if(use_pflotran .and. pf_cmode) then
            sum3 = sum3 + this%smin_nh4sorb_vr(c,j) * dzsoi_decomp(j)
         end if
      end do
      this%smin_no3(c) = sum1 
      this%smin_nh4(c) = sum2 
      if(use_pflotran .and. pf_cmode) this%smin_nh4sorb(c) = sum3 
   end do

   ! vertically integrate each of the decomposing N pools
   ! do fc = 1,num_soilc
   !    c = filter_soilc(fc)
   !    this%decomp_npools(c,l) = 0._r8
   ! end do
   !$acc parallel loop independent gang worker collapse(2) default(present) private(sum1)
   do l = 1, ndecomp_pools
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         
         !$acc loop vector reduction(+:sum1)
         do j = 1, nlev
            sum1 = sum1 + this%decomp_npools_vr(c,j,l) * dzsoi_decomp(j)
         end do
         this%decomp_npools(c,l) = sum1 
      end do
   end do

   ! for vertically-resolved soil biogeochemistry, calculate some diagnostics of carbon pools to a given depth
   if ( nlevdecomp > 1) then

      ! do l = 1, ndecomp_pools
      !    do fc = 1,num_soilc
      !       c = filter_soilc(fc)
      !       this%decomp_npools_1m(c,l) = 0._r8
      !    end do
      ! end do

      ! vertically integrate each of the decomposing n pools to 1 meter
      !$acc parallel loop independent gang worker collapse(2) default(present) private(sum1)
      do l = 1, ndecomp_pools
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            !$acc loop vector reduction(+:sum1) 
            do j = 1, nlevdecomp
               if ( zisoi(j) <= maxdepth ) then
                  sum1 = sum1 + this%decomp_npools_vr(c,j,l) * dzsoi_decomp(j)
               elseif ( zisoi(j-1) < maxdepth ) then
                  sum1 = sum1 + this%decomp_npools_vr(c,j,l) * (maxdepth - zisoi(j-1))
               endif
            end do
            this%decomp_npools_1m(c,l) = sum1  
         end do
      end do

      ! do fc = 1,num_soilc
      !    c = filter_soilc(fc)
      !    this%totlitn_1m(c) = 0._r8
      ! end do
      
      !$acc parallel loop independent gang worker  default(present) private(sum1,sum2)
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         sum1 = 0.0_r8 
         sum2 = 0.0_r8
         !$acc loop vector reduction(+:sum1,sum2) 
         do l = 1, ndecomp_pools
            ! total litter nitrogen to 1 meter (TOTLITN_1m)
            if ( decomp_cascade_con%is_litter(l) ) then
               sum1 = sum1 + this%decomp_npools_1m(c,l)
            end if
            ! total soil organic matter nitrogen to 1 meter (TOTSOMN_1m)
            if ( decomp_cascade_con%is_soil(l) ) then
               sum2 = sum2 + this%decomp_npools_1m(c,l)
            end if 

         end do
         this%totlitn_1m(c) = sum1
         this%totsomn_1m(c) = sum2 

      end do

      ! ! total soil organic matter nitrogen to 1 meter (TOTSOMN_1m)
      ! do fc = 1,num_soilc
      !    c = filter_soilc(fc)
      !    this%totsomn_1m(c) = 0._r8
      ! end do
      ! do l = 1, ndecomp_pools
      !    if ( decomp_cascade_con%is_soil(l) ) then
      !       do fc = 1,num_soilc
      !          c = filter_soilc(fc)
      !          this%totsomn_1m(c) = &
      !               this%totsomn_1m(c) + &
      !               this%decomp_npools_1m(c,l)
      !       end do
      !    end if
      ! end do

   endif

   ! ! total litter nitrogen (TOTLITN)
   ! do fc = 1,num_soilc
   !    c = filter_soilc(fc)
   !    this%totlitn(c)    = 0._r8
   ! end do
   !$acc parallel loop independent gang worker default(present) private(sum1,sum2)
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      sum1 = 0.0_r8
      sum2 = 0.0_r8 
      sum3 = 0.0_r8
      !$acc loop vector reduction(+:sum1,sum2,sum3)
      do l = 1, ndecomp_pools
         if ( decomp_cascade_con%is_litter(l) ) then
            sum1 = sum1 + this%decomp_npools(c,l)
         end if
         if ( decomp_cascade_con%is_soil(l) ) then
            sum2 = sum2 + this%decomp_npools(c,l)
         end if 
         if ( decomp_cascade_con%is_cwd(l) ) then
            sum3 = sum3 + this%decomp_npools(c,l)
         endif 
      end do
      this%totlitn(c) = sum1 
      this%totsomn(c) = sum2 
      this%cwdn(c)    = sum3 
   end do

   ! ! total soil organic matter nitrogen (TOTSOMN)
   ! do fc = 1,num_soilc
   !    c = filter_soilc(fc)
   !    this%totsomn(c)    = 0._r8
   ! end do
   ! do l = 1, ndecomp_pools
   !    if ( decomp_cascade_con%is_soil(l) ) then
   !       do fc = 1,num_soilc
   !          c = filter_soilc(fc)
   !          this%totsomn(c) = &
   !               this%totsomn(c) + &
   !               this%decomp_npools(c,l)
   !       end do
   !    end if
   ! end do

   ! ! total cwdn
   ! do fc = 1,num_soilc
   !    c = filter_soilc(fc)
   !    this%cwdn(c) = 0._r8
   ! end do
   ! do l = 1, ndecomp_pools
   !    if ( decomp_cascade_con%is_cwd(l) ) then
   !       do fc = 1,num_soilc
   !          c = filter_soilc(fc)
   !          this%cwdn(c) = &
   !               this%cwdn(c) + &
   !               this%decomp_npools(c,l)
   !       end do
   !    end if
   ! end do

   ! total sminn
   ! do fc = 1,num_soilc
   !    c = filter_soilc(fc)
   !    this%sminn(c)      = 0._r8
   ! end do
   !$acc parallel loop independent gang worker default(present) private(sum1,sum2)
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      sum1 = 0.0_r8
      sum2 = 0.0_r8 
      do j = 1, nlev
         sum1 = sum1 + this%sminn_vr(c,j) * dzsoi_decomp(j)
         sum2 = sum2 + this%ntrunc_vr(c,j) * dzsoi_decomp(j)
      end do
      this%sminn(c) = sum1
      this%ntrunc(c) = sum2 

   end do

   ! ! total col_ntrunc
   ! do fc = 1,num_soilc
   !    c = filter_soilc(fc)
   !    this%ntrunc(c) = 0._r8
   ! end do
   ! do j = 1, nlev
   !    do fc = 1,num_soilc
   !       c = filter_soilc(fc)
   !       this%ntrunc(c) = &
   !            this%ntrunc(c) + &
   !            this%ntrunc_vr(c,j) * dzsoi_decomp(j)
   !    end do
   ! end do

   !$acc parallel loop independent gang vector default(present) 
   do fc = 1,num_soilc
      c = filter_soilc(fc)

      ! total wood product nitrogen
      this%totprodn(c) = &
           this%prod1n(c) + &
           this%prod10n(c) + &
           this%prod100n(c)

      ! total ecosystem nitrogen, including veg (TOTECOSYSN)
      this%totecosysn(c) = &
           this%cwdn(c) + &
           this%totlitn(c) + &
           this%totsomn(c) + &
           this%sminn(c) + &
           this%totprodn(c) + &
           this%totvegn(c)

      ! total column nitrogen, including pft (TOTCOLN)
      this%totcoln(c) = &
           this%totpftn(c) + &
           this%cwdn(c) + &
           this%totlitn(c) + &
           this%totsomn(c) + &
           this%sminn(c) + &
           this%totprodn(c) + &
           this%ntrunc(c)+ &
           this%plant_n_buffer(c) + &
           this%cropseedn_deficit(c)

      this%totabgn (c) =  &
           this%totpftn(c) + &
           this%totprodn(c) + &
           this%seedn(c) + &
           this%ntrunc(c)+ &
           this%plant_n_buffer(c)

      this%totblgn(c) = &
           this%cwdn(c) + &
           this%totlitn(c) + &
           this%totsomn(c) + &
           this%sminn(c)
   end do
   !$acc exit data delete(sum1,sum2,sum3)

 end subroutine col_ns_summary_acc

 !-----------------------------------------------------------------------
 subroutine col_ps_summary_acc(this, bounds, num_soilc, filter_soilc)
   ! !ARGUMENTS:
   type(column_phosphorus_state) :: this
   type(bounds_type) , intent(in)  :: bounds
   integer           , intent(in)  :: num_soilc       ! number of soil columns in filter
   integer           , intent(in)  :: filter_soilc(:) ! filter for soil columns
   !
   ! !LOCAL VARIABLES:
   integer  :: c,j,k,l  ! indices
   integer  :: fc       ! lake filter indices
   real(r8), parameter :: maxdepth=1.0_r8  ! depth to integrate soil variables
   real(r8) :: sum1, sum2, sum3, sum4, sum5 
   !-----------------------------------------------------------------------
   !   do fc = 1,num_soilc
   !      c = filter_soilc(fc)
   !      this%solutionp(c) =0._r8
   !      this%labilep(c)   =0._r8
   !      this%secondp(c)   =0._r8
   !      this%occlp(c)     =0._r8
   !      this%primp(c)     =0._r8
   !   end do
   !$acc enter data create( sum1, sum2, sum3, sum4, sum5)
  !$acc parallel loop independent gang vector default(present) private(sum1,sum2,sum3,sum4,sum5) 
  do fc = 1,num_soilc
   c = filter_soilc(fc)
   sum1 = 0.0_r8 
   sum2 = 0.0_r8 
   sum3 = 0.0_r8 
   sum4 = 0.0_r8 
   sum5 = 0.0_r8 
   !$acc loop vector reduction(+:sum1,sum2,sum3,sum4,sum5) 
   do j = 1, nlevdecomp
        sum1 = sum1 + this%solutionp_vr(c,j) * dzsoi_decomp(j)
        sum2 = sum2 + this%labilep_vr(c,j) * dzsoi_decomp(j)
        sum3 = sum3 + this%secondp_vr(c,j) * dzsoi_decomp(j)
        sum4 = sum4 + this%occlp_vr(c,j) * dzsoi_decomp(j)
        sum5 = sum5 + this%primp_vr(c,j) * dzsoi_decomp(j)
     end do
     this%solutionp(c) = sum1 
     this%labilep(c)   = sum2 
     this%secondp(c)   = sum3 
     this%occlp(c)     = sum4 
     this%primp(c)     = sum5 
  end do

  ! vertically integrate each of the decomposing P pools
  !   do fc = 1,num_soilc
  !      c = filter_soilc(fc)
  !      this%decomp_ppools(c,l) = 0._r8
  !   end do

  !$acc parallel loop independent collapse(2) gang worker default(present) private(sum1)
  do l = 1, ndecomp_pools
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         sum1 = 0.0_r8 
         !$acc loop vector reduction(+:sum1)
         do j = 1, nlevdecomp
           sum1 = sum1 + this%decomp_ppools_vr(c,j,l) * dzsoi_decomp(j)
        end do
        this%decomp_ppools(c,l) = sum1 
     end do
  end do

  ! for vertically-resolved soil biogeochemistry, calculate some diagnostics of carbon pools to a given depth
  if ( nlevdecomp > 1) then
   !   do l = 1, ndecomp_pools
   !      do fc = 1,num_soilc
   !         c = filter_soilc(fc)
   !         this%decomp_ppools_1m(c,l) = 0._r8
   !      end do
   !   end do

     ! vertically integrate each of the decomposing n pools to 1 meter
      !$acc parallel loop independent gang worker default(present) private(sum1)
      do l = 1, ndecomp_pools
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            sum1 = 0.0_r8 
            !$acc loop vector reduction(+:sum1) 
            do j = 1, nlevdecomp
               if ( zisoi(j) <= maxdepth ) then
                 sum1 = sum1 + this%decomp_ppools_vr(c,j,l) * dzsoi_decomp(j)
               elseif ( zisoi(j-1) < maxdepth ) then
                 sum1 = sum1 + this%decomp_ppools_vr(c,j,l) * (maxdepth - zisoi(j-1))
               endif
            end do
            this%decomp_ppools_1m(c,l) = sum1 
        end do
     end do

     ! total litter phosphorus to 1 meter (TOTLITN_1m)
   !   do fc = 1,num_soilc
   !      c = filter_soilc(fc)
   !      this%totlitp_1m(c) = 0._r8
   !   end do
      
      !$acc parallel loop independent gang worker default(present) private(sum1,sum2)
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         sum1 = 0.0_r8 
         sum2 = 0.0_r8 
         !$acc loop vector reduction(+:sum1) 
         do l = 1, ndecomp_pools
            if ( decomp_cascade_con%is_litter(l) ) then
              sum1 = sum1 + this%decomp_ppools_1m(c,l)
            end if
            ! total soil organic matter phosphorus to 1 meter (TOTSOMN_1m)
            if ( decomp_cascade_con%is_soil(l) ) then
               sum2 = sum2 + this%decomp_ppools_1m(c,l)
            end if 

         end do
         this%totlitp_1m(c) = sum1 
         this%totsomp_1m(c) = sum2 
      end do

   !   ! total soil organic matter phosphorus to 1 meter (TOTSOMN_1m)
   !   do fc = 1,num_soilc
   !      c = filter_soilc(fc)
   !      this%totsomp_1m(c) = 0._r8
   !   end do
   !   do l = 1, ndecomp_pools
   !      if ( decomp_cascade_con%is_soil(l) ) then
   !         do fc = 1,num_soilc
   !            c = filter_soilc(fc)
   !            this%totsomp_1m(c) = &
   !                 this%totsomp_1m(c) + &
   !                 this%decomp_ppools_1m(c,l)
   !         end do
   !      end if
   !   end do

  endif

  ! total litter phosphorus (TOTLITP)
   !   do fc = 1,num_soilc
   !      c = filter_soilc(fc)
   !      this%totlitp(c)    = 0._r8
   !   end do
  !$acc parallel loop independent gang worker default(present) private(sum1,sum2, sum3) 
  do fc = 1,num_soilc
      c = filter_soilc(fc)
      sum1 = 0.0_r8;
      sum2 = 0.0_r8;
      sum3 = 0.0_r8;

      !$acc loop vector reduction(+:sum1)
      do l = 1, ndecomp_pools
         ! total litter phosphorus (TOTLITP)
         if ( decomp_cascade_con%is_litter(l) ) then
           sum1 = sum1 + this%decomp_ppools(c,l)
         end if
         ! total soil organic matter phosphorus (TOTSOMP)
         if ( decomp_cascade_con%is_soil(l) ) then 
            sum2 = sum2 + this%decomp_ppools(c,l) 
         end if 

         ! total cwdn
         if ( decomp_cascade_con%is_cwd(l) ) then
            sum3 = sum3 + this%decomp_ppools(c,l)
         end if 
      end do
      this%totlitp(c) = sum1 
      this%totsomp(c) = sum2 
      this%cwdp(c)    = sum3 
  end do

!   ! total soil organic matter phosphorus (TOTSOMP)
!   do fc = 1,num_soilc
!      c = filter_soilc(fc)
!      this%totsomp(c)    = 0._r8
!   end do
!   do l = 1, ndecomp_pools
!      if ( decomp_cascade_con%is_soil(l) ) then
!         do fc = 1,num_soilc
!            c = filter_soilc(fc)
!            this%totsomp(c) = &
!                 this%totsomp(c) + &
!                 this%decomp_ppools(c,l)
!         end do
!      end if
!   end do

!   ! total cwdn
!   do fc = 1,num_soilc
!      c = filter_soilc(fc)
!      this%cwdp(c) = 0._r8
!   end do
!   do l = 1, ndecomp_pools
!      if ( decomp_cascade_con%is_cwd(l) ) then
!         do fc = 1,num_soilc
!            c = filter_soilc(fc)
!            this%cwdp(c) = &
!                 this%cwdp(c) + &
!                 this%decomp_ppools(c,l)
!         end do
!      end if
!   end do

  ! total sminp
!   do fc = 1,num_soilc
!      c = filter_soilc(fc)
!      this%sminp(c)      = 0._r8
!   end do
  !$acc parallel loop independent collapse(2) gang vector default(present)
  do j = 1, nlevdecomp
     do fc = 1,num_soilc
        c = filter_soilc(fc)
        this%sminp_vr(c,j) = this%solutionp_vr(c,j) + &
                                 this%labilep_vr(c,j) + &
                                 this%secondp_vr(c,j)
     end do
  end do

   !$acc parallel loop independent gang worker default(present) private(sum1,sum2)
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      sum1 = 0.0_r8
      sum2 = 0.0_r8 
      !$acc loop vector reduction(+:sum1,sum2) 
      do j = 1, nlevdecomp
        sum1 = sum1 + this%sminp_vr(c,j) * dzsoi_decomp(j)
        sum2 = sum2 + this%ptrunc_vr(c,j) * dzsoi_decomp(j)

      end do
      this%sminp(c) = sum1 
      this%ptrunc(c) = sum2 
   end do

  ! total col_ptrunc
!   do fc = 1,num_soilc
!      c = filter_soilc(fc)
!      this%ptrunc(c) = 0._r8
!   end do

!   do j = 1, nlevdecomp
!      do fc = 1,num_soilc
!         c = filter_soilc(fc)
!         this%ptrunc(c) = &
!              this%ptrunc(c) + &
!              this%ptrunc_vr(c,j) * dzsoi_decomp(j)
!      end do
!   end do

   !$acc parallel loop independent gang vector default(present) 
  do fc = 1,num_soilc
     c = filter_soilc(fc)

     ! total wood product phosphorus
     this%totprodp(c) = &
          this%prod1p(c) + &
          this%prod10p(c) + &
          this%prod100p(c)

     ! total ecosystem phosphorus, including veg (TOTECOSYSP)
     this%totecosysp(c) = &
          this%cwdp(c) + &
          this%totlitp(c) + &
          this%totsomp(c) + &
          this%solutionp(c) + &
          this%labilep(c) + &
          this%secondp(c) + &
          this%primp(c) + &
          this%occlp(c) + &
          this%totprodp(c) + &
          this%totvegp(c)

     ! total column phosphorus, including pft (TOTCOLP)
     this%totcolp(c) = &
          this%totpftp(c) + &
          this%cwdp(c) + &
          this%totlitp(c) + &
          this%totsomp(c) + &
          this%totprodp(c) + &
          this%solutionp(c) + &
          this%labilep(c) + &
          this%secondp(c) + &
          this%ptrunc(c) + &
          this%cropseedp_deficit(c)
  end do

  !$acc exit data delete( sum1, sum2, sum3, sum4, sum5)

 end subroutine col_ps_summary_acc

end module ColumnWorkRoutinesMod
