module Method_procs_acc
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_const_mod   , only : SHR_CONST_TKFRZ
  use shr_const_mod   , only : SHR_CONST_PDB
  use clm_varpar      , only : nlevsoi, nlevsno, nlevgrnd, nlevlak, nlevurb
  use clm_varpar      , only : ndecomp_cascade_transitions, ndecomp_pools, nlevcan
  use clm_varpar      , only : nlevdecomp_full, crop_prog, nlevdecomp
  use clm_varpar      , only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
  use clm_varcon      , only : spval, ispval, zlnd, snw_rds_min, denice, denh2o, tfrz, pondmx
  use clm_varcon      , only : watmin, bdsno, zsoi, zisoi, dzsoi_decomp
  use clm_varcon      , only : c13ratio, c14ratio, secspday
  use clm_varctl      , only : use_fates, use_fates_planthydro, create_glacier_mec_landunit
  use clm_varctl      , only : bound_h2osoi, use_cn, iulog, use_vertsoilc, spinup_state
  use clm_varctl      , only : use_clm_interface, use_pflotran, pf_cmode
  use clm_varctl      , only : hist_wrtch4diag, use_nitrif_denitrif, use_century_decomp
  use clm_varctl      , only : get_carbontag, override_bgc_restart_mismatch_dump
  use clm_varctl      , only : pf_hmode, nu_com
  use ch4varcon       , only : allowlakeprod
  use pftvarcon       , only : VMAX_MINSURF_P_vr, KM_MINSURF_P_vr
  use pftvarcon       , only : npcropmin, noveg, nstor
  use subgridAveMod   , only : p2c,p2c_1d_filter

  use soilorder_varcon, only : smax, ks_sorption
  use landunit_varcon , only : istice, istwet, istsoil, istdlak, istcrop, istice_mec
  use column_varcon   , only : icol_road_perv, icol_road_imperv, icol_roof, icol_sunwall, icol_shadewall
  use decompMod       , only : bounds_type
  use CNStateType     , only: cnstate_type
  use tracer_varcon   , only : is_active_betr_bgc
  use CNDecompCascadeConType , only : decomp_cascade_con
  use ColumnType      , only : col_pp
  use LandunitType    , only : lun_pp
  use VegetationType, only : veg_pp
  use ColumnDataType , only : column_nitrogen_flux, column_nitrogen_state, column_phosphorus_flux
  use ColumnDataType , only : column_phosphorus_state, column_carbon_flux, column_carbon_state
  use ColumnDataType , only : nfix_timeconst, column_water_flux
  use VegetationDataType, only: vegetation_nitrogen_flux, vegetation_nitrogen_state, vegetation_phosphorus_flux
  use VegetationDataType ,only : vegetation_phosphorus_state, vegetation_carbon_flux, vegetation_carbon_state
  !
  ! !PUBLIC TYPES:
  implicit none
  integer, parameter :: c13 = 0, c14 =1 , bulk = 2
  public

  contains

    !-----------------------------------------------------------------------
    subroutine colcf_setvalues_acc (colcf, num_column, filter_column, value_column)
      !$acc routine seq
      ! !DESCRIPTION:
      ! Set column-level carbon fluxes
      use ColumnDataType , only : column_carbon_flux
      !
      ! !ARGUMENTS:
      type(column_carbon_flux), intent(inout) :: colcf
      integer , intent(in) :: num_column
      integer , intent(in) :: filter_column(:)
      real(r8), intent(in) :: value_column
      !
      ! !LOCAL VARIABLES:
      integer :: fi,i,j,k,l     ! loop index
      !------------------------------------------------------------------------

      do j = 1, nlevdecomp_full
         do fi = 1,num_column
            i = filter_column(fi)

            colcf%phenology_c_to_litr_met_c(i,j)     = value_column
            colcf%phenology_c_to_litr_cel_c(i,j)     = value_column
            colcf%phenology_c_to_litr_lig_c(i,j)     = value_column

            colcf%gap_mortality_c_to_litr_met_c(i,j) = value_column
            colcf%gap_mortality_c_to_litr_cel_c(i,j) = value_column
            colcf%gap_mortality_c_to_litr_lig_c(i,j) = value_column
            colcf%gap_mortality_c_to_cwdc(i,j)       = value_column

            colcf%fire_mortality_c_to_cwdc(i,j)      = value_column
            colcf%m_c_to_litr_met_fire(i,j)          = value_column
            colcf%m_c_to_litr_cel_fire(i,j)          = value_column
            colcf%m_c_to_litr_lig_fire(i,j)          = value_column

            colcf%harvest_c_to_litr_met_c(i,j)       = value_column
            colcf%harvest_c_to_litr_cel_c(i,j)       = value_column
            colcf%harvest_c_to_litr_lig_c(i,j)       = value_column
            colcf%harvest_c_to_cwdc(i,j)             = value_column

            colcf%hr_vr(i,j)                         = value_column
         end do
      end do

      do k = 1, ndecomp_pools
         do j = 1, nlevdecomp_full
            do fi = 1,num_column
               i = filter_column(fi)
               colcf%m_decomp_cpools_to_fire_vr(i,j,k) = value_column
               colcf%decomp_cpools_transport_tendency(i,j,k) = value_column
            end do
         end do
      end do

      do l = 1, ndecomp_cascade_transitions
         do fi = 1,num_column
            i = filter_column(fi)
            colcf%decomp_cascade_hr(i,l) = value_column
            colcf%decomp_cascade_ctransfer(i,l) = value_column
         end do
      end do

      do l = 1, ndecomp_cascade_transitions
         do j = 1, nlevdecomp_full
            do fi = 1,num_column
               i = filter_column(fi)
               colcf%decomp_cascade_hr_vr(i,j,l) = value_column
               colcf%decomp_cascade_ctransfer_vr(i,j,l) = value_column
               colcf%decomp_k(i,j,l) = value_column
            end do
         end do
      end do

      do k = 1, ndecomp_pools
         do fi = 1,num_column
            i = filter_column(fi)
            colcf%decomp_cpools_leached(i,k) = value_column
            colcf%m_decomp_cpools_to_fire(i,k) = value_column
            colcf%bgc_cpool_ext_inputs_vr(i,:, k) = value_column
            colcf%bgc_cpool_ext_loss_vr(i,:, k) = value_column
         end do
      end do

      do fi = 1,num_column
         i = filter_column(fi)

         colcf%hrv_deadstemc_to_prod10c(i)  = value_column
         colcf%hrv_deadstemc_to_prod100c(i) = value_column
         colcf%hrv_cropc_to_prod1c(i)       = value_column
         colcf%somc_fire(i)                 = value_column
         colcf%prod1c_loss(i)               = value_column
         colcf%prod10c_loss(i)              = value_column
         colcf%prod100c_loss(i)             = value_column
         colcf%product_closs(i)             = value_column
         colcf%somhr(i)                     = value_column
         colcf%lithr(i)                     = value_column
         colcf%hr(i)                        = value_column
         colcf%sr(i)                        = value_column
         colcf%er(i)                        = value_column
         colcf%litfire(i)                   = value_column
         colcf%somfire(i)                   = value_column
         colcf%totfire(i)                   = value_column
         colcf%nep(i)                       = value_column
         colcf%nbp(i)                       = value_column
         colcf%nee(i)                       = value_column
         colcf%fire_closs(i)                = value_column
         colcf%cwdc_hr(i)                   = value_column
         colcf%cwdc_loss(i)                 = value_column
         colcf%litterc_loss(i)              = value_column
         colcf%som_c_leached(i)             = value_column

         ! Zero p2c column fluxes
         colcf%rr(i)                    = value_column
         colcf%ar(i)                    = value_column
         colcf%gpp(i)                   = value_column
         colcf%npp(i)                   = value_column
         colcf%fire_closs(i)            = value_column
         colcf%litfall(i)               = value_column
         colcf%vegfire(i)               = value_column
         colcf%wood_harvestc(i)         = value_column
         colcf%hrv_xsmrpool_to_atm(i)   = value_column
      end do

      do k = 1, ndecomp_pools
         do j = 1, nlevdecomp_full
            do fi = 1,num_column
               i = filter_column(fi)
               colcf%decomp_cpools_sourcesink(i,j,k) = value_column
            end do
         end do
      end do

      ! pflotran
      do k = 1, ndecomp_pools
         do j = 1, nlevdecomp_full
            do fi = 1,num_column
               i = filter_column(fi)
               ! only initializing in the first time-step
               if ( colcf%externalc_to_decomp_cpools(i,j,k) == spval ) then
                  colcf%externalc_to_decomp_cpools(i,j,k) = value_column
               end if
            end do
         end do
      end do

      do fi = 1,num_column
         i = filter_column(fi)
         colcf%f_co2_soil(i) = value_column
         ! only initializing in the first time-step
         if ( colcf%externalc_to_decomp_delta(i) == spval ) then
            colcf%externalc_to_decomp_delta(i) = value_column
         end if
      end do

      do j = 1, nlevdecomp_full
         do fi = 1,num_column
            i = filter_column(fi)
            colcf%f_co2_soil_vr(i,j) = value_column
         end do
      end do

    end subroutine colcf_setvalues_acc

    !-----------------------------------------------------------------------
    subroutine colnf_setvalues_acc( colnf, num_column, filter_column, value_column)
      !$acc routine seq
      ! !DESCRIPTION:
      ! Set column-level nitrogen fluxes
      !
      use ColumnDataType,  only : column_nitrogen_flux
      ! !ARGUMENTS:
      type(column_nitrogen_flux), intent(inout) :: colnf
      integer , intent(in)         :: num_column
      integer , intent(in)         :: filter_column(:)
      real(r8), intent(in)         :: value_column
      !
      ! !LOCAL VARIABLES:
      integer :: fi,i,j,k,l     ! loop index
      !------------------------------------------------------------------------
      do j = 1, nlevdecomp_full
         do fi = 1,num_column
            i = filter_column(fi)

            ! phenology: litterfall and crop fluxes associated wit
            colnf%phenology_n_to_litr_met_n(i,j)        = value_column
            colnf%phenology_n_to_litr_cel_n(i,j)        = value_column
            colnf%phenology_n_to_litr_lig_n(i,j)        = value_column

            ! gap mortality
            colnf%gap_mortality_n_to_litr_met_n(i,j)    = value_column
            colnf%gap_mortality_n_to_litr_cel_n(i,j)    = value_column
            colnf%gap_mortality_n_to_litr_lig_n(i,j)    = value_column
            colnf%gap_mortality_n_to_cwdn(i,j)          = value_column

            ! fire
            colnf%fire_mortality_n_to_cwdn(i,j)         = value_column
            colnf%m_n_to_litr_met_fire(i,j)             = value_column
            colnf%m_n_to_litr_cel_fire(i,j)             = value_column
            colnf%m_n_to_litr_lig_fire(i,j)             = value_column

            ! harvest
            colnf%harvest_n_to_litr_met_n(i,j)          = value_column
            colnf%harvest_n_to_litr_cel_n(i,j)          = value_column
            colnf%harvest_n_to_litr_lig_n(i,j)          = value_column
            colnf%harvest_n_to_cwdn(i,j)                = value_column

            if (.not. use_nitrif_denitrif) then
               colnf%sminn_to_denit_excess_vr(i,j)      = value_column
               colnf%sminn_leached_vr(i,j)              = value_column
            else
               colnf%f_nit_vr(i,j)                      = value_column
               colnf%f_denit_vr(i,j)                    = value_column
               colnf%smin_no3_leached_vr(i,j)           = value_column
               colnf%smin_no3_runoff_vr(i,j)            = value_column
               colnf%n2_n2o_ratio_denit_vr(i,j)         = value_column
               colnf%pot_f_nit_vr(i,j)                  = value_column
               colnf%pot_f_denit_vr(i,j)                = value_column
               colnf%actual_immob_no3_vr(i,j)           = value_column
               colnf%actual_immob_nh4_vr(i,j)           = value_column
               colnf%smin_no3_to_plant_vr(i,j)          = value_column
               colnf%smin_nh4_to_plant_vr(i,j)          = value_column
               colnf%f_n2o_denit_vr(i,j)                = value_column
               colnf%f_n2o_nit_vr(i,j)                  = value_column

               colnf%smin_no3_massdens_vr(i,j)          = value_column
               colnf%k_nitr_t_vr(i,j)                   = value_column
               colnf%k_nitr_ph_vr(i,j)                  = value_column
               colnf%k_nitr_h2o_vr(i,j)                 = value_column
               colnf%k_nitr_vr(i,j)                     = value_column
               colnf%wfps_vr(i,j)                       = value_column
               colnf%fmax_denit_carbonsubstrate_vr(i,j) = value_column
               colnf%fmax_denit_nitrate_vr(i,j)         = value_column
               colnf%f_denit_base_vr(i,j)               = value_column

               colnf%diffus(i,j)                        = value_column
               colnf%ratio_k1(i,j)                      = value_column
               colnf%ratio_no3_co2(i,j)                 = value_column
               colnf%soil_co2_prod(i,j)                 = value_column
               colnf%fr_WFPS(i,j)                       = value_column
               colnf%soil_bulkdensity(i,j)              = value_column

               colnf%r_psi(i,j)                         = value_column
               colnf%anaerobic_frac(i,j)                = value_column

               ! pflotran
               colnf%plant_ndemand_vr(i,j)              = value_column
               colnf%f_ngas_decomp_vr(i,j)              = value_column
               colnf%f_ngas_nitri_vr(i,j)               = value_column
               colnf%f_ngas_denit_vr(i,j)               = value_column
               colnf%f_n2o_soil_vr(i,j)                 = value_column
               colnf%f_n2_soil_vr(i,j)                  = value_column
            end if
            colnf%potential_immob_vr(i,j)               = value_column
            colnf%actual_immob_vr(i,j)                  = value_column
            colnf%sminn_to_plant_vr(i,j)                = value_column
            colnf%supplement_to_sminn_vr(i,j)           = value_column
            colnf%gross_nmin_vr(i,j)                    = value_column
            colnf%net_nmin_vr(i,j)                      = value_column
            colnf%sminn_nh4_input_vr(i,j)               = value_column
            colnf%sminn_no3_input_vr(i,j)               = value_column
         end do
      end do

      do fi = 1,num_column
         i = filter_column(fi)

         colnf%ndep_to_sminn(i)             = value_column
         colnf%nfix_to_sminn(i)             = value_column
         colnf%nfix_to_ecosysn(i)           = value_column
         colnf%fert_to_sminn(i)             = value_column
         colnf%soyfixn_to_sminn(i)          = value_column
         colnf%hrv_deadstemn_to_prod10n(i)  = value_column
         colnf%hrv_deadstemn_to_prod100n(i) = value_column
         colnf%hrv_cropn_to_prod1n(i)       = value_column
         colnf%prod10n_loss(i)              = value_column
         colnf%prod100n_loss(i)             = value_column
         colnf%prod1n_loss(i)               = value_column
         colnf%product_nloss(i)             = value_column
         colnf%potential_immob(i)           = value_column
         colnf%actual_immob(i)              = value_column
         colnf%sminn_to_plant(i)            = value_column
         colnf%supplement_to_sminn(i)       = value_column
         colnf%gross_nmin(i)                = value_column
         colnf%net_nmin(i)                  = value_column
         colnf%denit(i)                     = value_column
         if (use_nitrif_denitrif) then
            colnf%f_nit(i)                  = value_column
            colnf%pot_f_nit(i)              = value_column
            colnf%f_denit(i)                = value_column
            colnf%pot_f_denit(i)            = value_column
            colnf%f_n2o_denit(i)            = value_column
            colnf%f_n2o_nit(i)              = value_column
            colnf%smin_no3_leached(i)       = value_column
            colnf%smin_no3_runoff(i)        = value_column

            colnf%f_ngas_decomp(i)         = value_column
            colnf%f_ngas_nitri(i)          = value_column
            colnf%f_ngas_denit(i)          = value_column
            colnf%f_n2o_soil(i)            = value_column
            colnf%f_n2_soil(i)             = value_column

            colnf%smin_nh4_to_plant(i)      = value_column
            colnf%smin_no3_to_plant(i)      = value_column
         else
            colnf%sminn_to_denit_excess(i)  = value_column
            colnf%sminn_leached(i)          = value_column
         end if
         colnf%ninputs(i)                   = value_column
         colnf%noutputs(i)                  = value_column
         colnf%fire_nloss(i)                = value_column
         colnf%som_n_leached(i)             = value_column
         colnf%sminn_input(i)               = value_column
         colnf%sminn_nh4_input(i)           = value_column
         colnf%sminn_no3_input(i)           = value_column
         ! Zero p2c column fluxes
         colnf%fire_nloss(i) = value_column
         colnf%wood_harvestn(i) = value_column

         ! bgc-interface
         colnf%plant_ndemand(i) = value_column
      end do

      do k = 1, ndecomp_pools
         do fi = 1,num_column
            i = filter_column(fi)
            colnf%decomp_npools_leached(i,k) = value_column
            colnf%m_decomp_npools_to_fire(i,k) = value_column
            colnf%bgc_npool_ext_inputs_vr(i,:,k) = value_column
            colnf%bgc_npool_ext_loss_vr(i,:,k) = value_column
            colnf%bgc_npool_inputs(i,k) = value_column
         end do
      end do

      do k = 1, ndecomp_pools
         do j = 1, nlevdecomp_full
            do fi = 1,num_column
               i = filter_column(fi)
               colnf%m_decomp_npools_to_fire_vr(i,j,k) = value_column
               colnf%decomp_npools_transport_tendency(i,j,k) = value_column
            end do
         end do
      end do

         do l = 1, ndecomp_cascade_transitions
            do fi = 1,num_column
               i = filter_column(fi)
               colnf%decomp_cascade_ntransfer(i,l) = value_column
               colnf%decomp_cascade_sminn_flux(i,l) = value_column
               if (.not. use_nitrif_denitrif) then
                  colnf%sminn_to_denit_decomp_cascade(i,l) = value_column
               end if
            end do
         end do

         do l = 1, ndecomp_cascade_transitions
            do j = 1, nlevdecomp_full
               do fi = 1,num_column
                  i = filter_column(fi)
                  colnf%decomp_cascade_ntransfer_vr(i,j,l) = value_column
                  colnf%decomp_cascade_sminn_flux_vr(i,j,l) = value_column
                  if (.not. use_nitrif_denitrif) then
                     colnf%sminn_to_denit_decomp_cascade_vr(i,j,l) = value_column
                  end if
               end do
            end do
         end do

      do k = 1, ndecomp_pools
         do j = 1, nlevdecomp_full
            do fi = 1,num_column
               i = filter_column(fi)
               colnf%decomp_npools_sourcesink(i,j,k) = value_column
            end do
         end do
      end do

      ! pflotran
      !------------------------------------------------------------------------
      do k = 1, ndecomp_pools
         do j = 1, nlevdecomp_full
            do fi = 1,num_column
               i = filter_column(fi)
               ! only initializing in the first time-step
               if ( colnf%externaln_to_decomp_npools(i,j,k) == spval ) then
                  colnf%externaln_to_decomp_npools(i,j,k) = value_column
               end if
            end do
         end do
      end do

      do j = 1, nlevdecomp_full
         do fi = 1,num_column
            i = filter_column(fi)
            ! only initializing in the first time-step
            if ( colnf%no3_net_transport_vr(i,j) == spval ) then
               colnf%no3_net_transport_vr(i,j) = value_column
            end if
            if ( colnf%nh4_net_transport_vr(i,j) == spval ) then
               colnf%nh4_net_transport_vr(i,j) = value_column
            end if
         end do
      end do

      do fi = 1,num_column
         i = filter_column(fi)
         ! only initializing in the first time-step
         if ( colnf%externaln_to_decomp_delta(i) == spval ) then
            colnf%externaln_to_decomp_delta(i) = value_column
         end if
      end do

    end subroutine colnf_setvalues_acc

    !-----------------------------------------------------------------------
    subroutine colnf_summary_acc(colnf, bounds, num_soilc, filter_soilc, dtime)
      !$acc routine seq
      ! !DESCRIPTION:
      ! Column-level nitrogen summary calculations
      !
      ! !ARGUMENTS:
      type(column_nitrogen_flux)            :: colnf
      type(bounds_type)      , intent(in)    :: bounds
      integer                , intent(in)    :: num_soilc       ! number of soil columns in filter
      integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
      real(r8)               , intent(in)    :: dtime       ! radiation time step (seconds)

      !
      ! !LOCAL VARIABLES:
      integer  :: fc,c,j,k,l       ! indices
      integer  :: nlev
      !-----------------------------------------------------------------------

      nlev = nlevdecomp
      if (use_pflotran .and. pf_cmode) nlev = nlevdecomp_full

      do fc = 1,num_soilc
         c = filter_soilc(fc)
         colnf%denit(c) = 0._r8
         colnf%supplement_to_sminn(c) = 0._r8
         colnf%som_n_leached(c)       = 0._r8
      end do

      if (  (.not. (use_pflotran .and. pf_cmode)) ) then
         ! BeTR is off AND PFLOTRAN's pf_cmode is false
         ! vertically integrate decomposing N cascade fluxes and
         !soil mineral N fluxes associated with decomposition cascade
         do k = 1, ndecomp_cascade_transitions
            do j = 1,nlev
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  colnf%decomp_cascade_ntransfer(c,k) = &
                       colnf%decomp_cascade_ntransfer(c,k) + &
                       colnf%decomp_cascade_ntransfer_vr(c,j,k) * dzsoi_decomp(j)

                  colnf%decomp_cascade_sminn_flux(c,k) = &
                       colnf%decomp_cascade_sminn_flux(c,k) + &
                       colnf%decomp_cascade_sminn_flux_vr(c,j,k) * dzsoi_decomp(j)
               end do
            end do
         end do

         if (.not. use_nitrif_denitrif) then
            ! vertically integrate each denitrification flux
            do l = 1, ndecomp_cascade_transitions
               do j = 1, nlev
                  do fc = 1,num_soilc
                     c = filter_soilc(fc)
                     colnf%sminn_to_denit_decomp_cascade(c,l) = &
                          colnf%sminn_to_denit_decomp_cascade(c,l) + &
                          colnf%sminn_to_denit_decomp_cascade_vr(c,j,l) * dzsoi_decomp(j)
                  end do
               end do
            end do
            ! vertically integrate bulk denitrification and  leaching flux
            do j = 1, nlev
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  colnf%sminn_to_denit_excess(c) = &
                       colnf%sminn_to_denit_excess(c) + &
                       colnf%sminn_to_denit_excess_vr(c,j) * dzsoi_decomp(j)

                  colnf%sminn_leached(c) = &
                       colnf%sminn_leached(c) + &
                       colnf%sminn_leached_vr(c,j) * dzsoi_decomp(j)
               end do
            end do
            ! total N denitrification (DENIT)
            do l = 1, ndecomp_cascade_transitions
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  colnf%denit(c) = &
                       colnf%denit(c) + &
                       colnf%sminn_to_denit_decomp_cascade(c,l)
               end do
            end do
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               colnf%denit(c) =  &
                    colnf%denit(c) + &
                    colnf%sminn_to_denit_excess(c)
            end do
         else
            ! vertically integrate NO3 NH4 N2O fluxes and pools
            do j = 1, nlev
               do fc = 1,num_soilc
                  c = filter_soilc(fc)

                  ! nitrification and denitrification fluxes
                  colnf%f_nit(c) = &
                       colnf%f_nit(c) + &
                       colnf%f_nit_vr(c,j) * dzsoi_decomp(j)

                  colnf%f_denit(c) = &
                       colnf%f_denit(c) + &
                       colnf%f_denit_vr(c,j) * dzsoi_decomp(j)

                  colnf%pot_f_nit(c) = &
                       colnf%pot_f_nit(c) + &
                       colnf%pot_f_nit_vr(c,j) * dzsoi_decomp(j)

                  colnf%pot_f_denit(c) = &
                       colnf%pot_f_denit(c) + &
                       colnf%pot_f_denit_vr(c,j) * dzsoi_decomp(j)

                  colnf%f_n2o_nit(c) = &
                       colnf%f_n2o_nit(c) + &
                       colnf%f_n2o_nit_vr(c,j) * dzsoi_decomp(j)

                  colnf%f_n2o_denit(c) = &
                       colnf%f_n2o_denit(c) + &
                       colnf%f_n2o_denit_vr(c,j) * dzsoi_decomp(j)

                  ! leaching/runoff flux
                  colnf%smin_no3_leached(c) = &
                          colnf%smin_no3_leached(c) + &
                          colnf%smin_no3_leached_vr(c,j) * dzsoi_decomp(j)

                  colnf%smin_no3_runoff(c) = &
                          colnf%smin_no3_runoff(c) + &
                          colnf%smin_no3_runoff_vr(c,j) * dzsoi_decomp(j)
               end do
            end do
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               colnf%denit(c) = colnf%f_denit(c)
            end do
         end if
      end if

      ! vertically integrate column-level fire N losses
      do k = 1, ndecomp_pools
         do j = 1, nlev
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               colnf%m_decomp_npools_to_fire(c,k) = &
                    colnf%m_decomp_npools_to_fire(c,k) + &
                    colnf%m_decomp_npools_to_fire_vr(c,j,k) * dzsoi_decomp(j)
            end do
         end do
      end do

      ! total column-level fire N losses
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         colnf%fire_nloss(c) = colnf%fire_nloss_p2c(c)
      end do
      do k = 1, ndecomp_pools
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            colnf%fire_nloss(c) = &
                 colnf%fire_nloss(c) + &
                 colnf%m_decomp_npools_to_fire(c,k)
         end do
      end do

      ! supplementary N supplement_to_sminn
      do j = 1, nlev
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            colnf%supplement_to_sminn(c) = &
                 colnf%supplement_to_sminn(c) + &
                 colnf%supplement_to_sminn_vr(c,j) * dzsoi_decomp(j)

            colnf%sminn_input(c) = &
                 colnf%sminn_input(c) + &
                 (colnf%sminn_nh4_input_vr(c,j)+colnf%sminn_no3_input_vr(c,j))*dzsoi_decomp(j)

            colnf%sminn_nh4_input(c) = &
                 colnf%sminn_nh4_input(c) + &
                 colnf%sminn_nh4_input_vr(c,j)*dzsoi_decomp(j)

            colnf%sminn_no3_input(c) = &
                 colnf%sminn_no3_input(c) + &
                 colnf%sminn_no3_input_vr(c,j)*dzsoi_decomp(j)
         end do
      end do

      do fc = 1,num_soilc
         c = filter_soilc(fc)
         ! column-level N losses due to landcover change
         colnf%dwt_nloss(c) = &
              colnf%dwt_conv_nflux(c)
         ! total wood product N loss
         colnf%product_nloss(c) = &
              colnf%prod10n_loss(c) + &
              colnf%prod100n_loss(c)+ &
              colnf%prod1n_loss(c)
      end do

      ! add up all vertical transport tendency terms and calculate total som leaching loss as the sum of these
      do l = 1, ndecomp_pools
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            colnf%decomp_npools_leached(c,l) = 0._r8
         end do
         do j = 1, nlev
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               colnf%decomp_npools_leached(c,l) = &
                    colnf%decomp_npools_leached(c,l) + &
                    colnf%decomp_npools_transport_tendency(c,j,l) * dzsoi_decomp(j)

               colnf%bgc_npool_inputs(c,l) = colnf%bgc_npool_inputs(c,l) + &
                  (colnf%bgc_npool_ext_inputs_vr(c,j,l)-colnf%bgc_npool_ext_loss_vr(c,j,l))*dzsoi_decomp(j)
            end do
         end do
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            colnf%som_n_leached(c) = &
                 colnf%som_n_leached(c) + &
                 colnf%decomp_npools_leached(c,l)
         end do
      end do

      do fc = 1,num_soilc
         c = filter_soilc(fc)
         colnf%smin_no3_to_plant(c) = 0._r8
         colnf%smin_nh4_to_plant(c) = 0._r8
         colnf%plant_to_litter_nflux(c) = 0._r8
         colnf%plant_to_cwd_nflux(c) = 0._r8
         do j = 1, nlev
            colnf%plant_to_litter_nflux(c) = &
                 colnf%plant_to_litter_nflux(c)  + &
                 colnf%phenology_n_to_litr_met_n(c,j)* dzsoi_decomp(j) + &
                 colnf%phenology_n_to_litr_cel_n(c,j)* dzsoi_decomp(j) + &
                 colnf%phenology_n_to_litr_lig_n(c,j)* dzsoi_decomp(j) + &
                 colnf%gap_mortality_n_to_litr_met_n(c,j)* dzsoi_decomp(j) + &
                 colnf%gap_mortality_n_to_litr_cel_n(c,j)* dzsoi_decomp(j) + &
                 colnf%gap_mortality_n_to_litr_lig_n(c,j)* dzsoi_decomp(j) + &
                 colnf%m_n_to_litr_met_fire(c,j)* dzsoi_decomp(j) + &
                 colnf%m_n_to_litr_cel_fire(c,j)* dzsoi_decomp(j) + &
                 colnf%m_n_to_litr_lig_fire(c,j)* dzsoi_decomp(j)
            colnf%plant_to_cwd_nflux(c) = &
                 colnf%plant_to_cwd_nflux(c) + &
                 colnf%gap_mortality_n_to_cwdn(c,j)* dzsoi_decomp(j) + &
                 colnf%fire_mortality_n_to_cwdn(c,j)* dzsoi_decomp(j)
         end do
      end do

      if (use_nitrif_denitrif) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            do j = 1, nlev
               colnf%smin_no3_to_plant(c)= colnf%smin_no3_to_plant(c) + &
                    colnf%smin_no3_to_plant_vr(c,j) * dzsoi_decomp(j)
               colnf%smin_nh4_to_plant(c)= colnf%smin_nh4_to_plant(c) + &
                    colnf%smin_nh4_to_plant_vr(c,j) * dzsoi_decomp(j)
            enddo
         enddo
      endif

      ! bgc interface & pflotran
      if (use_clm_interface .and. (use_pflotran .and. pf_cmode)) then
          call colnf_SummaryInt_acc(colnf,bounds, num_soilc, filter_soilc, dtime)
      end if

    end subroutine colnf_summary_acc

    !------------------------------------------------------------------------
    subroutine colnf_summaryint_acc (colnf,bounds,num_soilc, filter_soilc, dtime)
      !
      ! !DESCRIPTION:
      ! Column-level nitrogen flux summary for PFLOTRAN interface
      !$acc routine seq
      ! !ARGUMENTS:
      implicit none
      type (column_nitrogen_flux)    :: colnf
      type(bounds_type) ,  intent(in) :: bounds
      integer,             intent(in) :: num_soilc       ! number of soil columns in filter
      integer,             intent(in) :: filter_soilc(:) ! filter for soil columns
      real(r8), intent(in    ):: dtime       ! radiation time step (seconds)

      !
      ! !LOCAL VARIABLES:
      integer :: c,j,l       ! indices
      integer :: fc          ! column filter indices
      !------------------------------------------------------------------------

      !#py dtime = real( get_step_size(), r8 )
      ! nitrification-denitrification rates (not yet passing out from PF, but will)
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         colnf%f_nit(c)   = 0._r8
         colnf%f_denit(c) = 0._r8
         do j = 1, nlevdecomp_full
            colnf%f_nit_vr(c,j) = 0._r8
            colnf%f_nit(c)  = colnf%f_nit(c) + &
                                 colnf%f_nit_vr(c,j)*dzsoi_decomp(j)

            colnf%f_denit_vr(c,j) = 0._r8
            colnf%f_denit(c) = colnf%f_denit(c) + &
                                 colnf%f_denit_vr(c,j)*dzsoi_decomp(j)

         end do
         colnf%denit(c)      = colnf%f_denit(c)
      end do

      ! the following are from pflotran bgc, and vertically down to 'nlevdecomp_full'
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         colnf%f_n2_soil(c)    = 0._r8
         colnf%f_n2o_soil(c)   = 0._r8
         colnf%f_ngas_decomp(c)= 0._r8
         colnf%f_ngas_nitri(c) = 0._r8
         colnf%f_ngas_denit(c) = 0._r8
         colnf%smin_no3_leached(c) = 0._r8
         colnf%smin_no3_runoff(c)  = 0._r8
         colnf%sminn_leached(c)    = 0._r8
         do j = 1, nlevdecomp_full
            ! all N2/N2O gas exchange between atm. and soil (i.e., dissolving - degassing)
            colnf%f_n2_soil(c)  = colnf%f_n2_soil(c) + &
                                  colnf%f_n2_soil_vr(c,j)*dzsoi_decomp(j)
            colnf%f_n2o_soil(c) = colnf%f_n2o_soil(c) + &
                                  colnf%f_n2o_soil_vr(c,j)*dzsoi_decomp(j)

            ! all N2/N2O production from soil bgc N processes (mineralization-nitrification-denitrification)
            ! note: those are directly dissolved into aq. gas species, which would be exchanging with atm.
            colnf%f_ngas_decomp(c) = colnf%f_ngas_decomp(c) + &
                                     colnf%f_ngas_decomp_vr(c,j)*dzsoi_decomp(j)
            colnf%f_ngas_nitri(c)  = colnf%f_ngas_nitri(c) + &
                                     colnf%f_ngas_nitri_vr(c,j)*dzsoi_decomp(j)
            colnf%f_ngas_denit(c)  = colnf%f_ngas_denit(c) + &
                                     colnf%f_ngas_denit_vr(c,j)*dzsoi_decomp(j)

            ! leaching/runoff fluxes summed vertically
            ! (1) if not hydroloy-coupled, advection from CLM-CN, plus diffusion from PF
            ! (2) if hydrology-coupled, all from PF (i.e. 'no3_net_transport_vr_col');
            colnf%smin_no3_leached(c) = colnf%smin_no3_leached(c) + &
                                        colnf%no3_net_transport_vr(c,j) * dzsoi_decomp(j)

            if(.not. pf_hmode) then ! this is from CLM-CN's leaching subroutine
                colnf%smin_no3_leached(c) = colnf%smin_no3_leached(c) + &
                                        colnf%smin_no3_leached_vr(c,j) * dzsoi_decomp(j)
                colnf%smin_no3_runoff(c)  = colnf%smin_no3_runoff(c) + &
                                        colnf%smin_no3_runoff_vr(c,j) * dzsoi_decomp(j)
            endif

            ! assign all no3-N leaching/runof,including diffusion from PF, to all mineral-N
            colnf%sminn_leached_vr(c,j) = colnf%smin_no3_leached_vr(c,j) + &
                                             colnf%smin_no3_runoff_vr(c,j) +  &
                                        colnf%nh4_net_transport_vr(c,j) * dzsoi_decomp(j)

            colnf%sminn_leached(c) = colnf%sminn_leached(c) + &
                                        colnf%sminn_leached_vr(c,j)*dzsoi_decomp(j)
         end do !j = 1, nlevdecomp_full
         ! for balance-checking
         colnf%denit(c)     = colnf%f_ngas_denit(c)
         colnf%f_n2o_nit(c) = colnf%f_ngas_decomp(c) + colnf%f_ngas_nitri(c)
      end do !fc = 1,num_soilc
      ! summarize at column-level vertically-resolved littering/removal for PFLOTRAN bgc input needs
      ! first it needs to save the total column-level N rate btw plant pool and decomposible pools at previous time step
      ! for adjusting difference when doing balance check

      do fc = 1,num_soilc
        c = filter_soilc(fc)
        colnf%externaln_to_decomp_delta(c) = 0._r8
        do j = 1, nlevdecomp_full
           do l = 1, ndecomp_pools
              colnf%externaln_to_decomp_delta(c) =    &
                 colnf%externaln_to_decomp_delta(c) + &
                   colnf%externaln_to_decomp_npools(c,j,l)*dzsoi_decomp(j)
           end do

        end do
      end do

      ! do the initialization for the following variable here.
      ! DON'T do so in the beginning of CLM-CN time-step (otherwise the above saved will not work)
      do fc = 1,num_soilc
           c = filter_soilc(fc)
           colnf%externaln_to_decomp_npools(c, 1:nlevdecomp_full, 1:ndecomp_pools) = 0._r8
      end do

      ! add up all vertically-resolved addition/removal rates (gC/m3/s) of decomp_pools
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         do j = 1, nlevdecomp_full
            do l = 1, ndecomp_pools
               ! for litter C pools
               if (l==i_met_lit) then
                  colnf%externaln_to_decomp_npools(c,j,l) =              &
                      colnf%externaln_to_decomp_npools(c,j,l)            &
                       + colnf%phenology_n_to_litr_met_n(c,j)            &
                       + colnf%dwt_frootn_to_litr_met_n(c,j)             &
                       + colnf%gap_mortality_n_to_litr_met_n(c,j)        &
                       + colnf%harvest_n_to_litr_met_n(c,j)              &
                       + colnf%m_n_to_litr_met_fire(c,j)

               elseif (l==i_cel_lit) then
                  colnf%externaln_to_decomp_npools(c,j,l) =              &
                      colnf%externaln_to_decomp_npools(c,j,l)            &
                       + colnf%phenology_n_to_litr_cel_n(c,j)            &
                       + colnf%dwt_frootn_to_litr_cel_n(c,j)             &
                       + colnf%gap_mortality_n_to_litr_cel_n(c,j)        &
                       + colnf%harvest_n_to_litr_cel_n(c,j)              &
                       + colnf%m_n_to_litr_cel_fire(c,j)

               elseif (l==i_lig_lit) then
                  colnf%externaln_to_decomp_npools(c,j,l) =              &
                      colnf%externaln_to_decomp_npools(c,j,l)            &
                       + colnf%phenology_n_to_litr_lig_n(c,j)            &
                       + colnf%dwt_frootn_to_litr_lig_n(c,j)             &
                       + colnf%gap_mortality_n_to_litr_lig_n(c,j)        &
                       + colnf%harvest_n_to_litr_lig_n(c,j)              &
                       + colnf%m_n_to_litr_lig_fire(c,j)

               ! for cwd
               elseif (l==i_cwd) then
                  colnf%externaln_to_decomp_npools(c,j,l) =              &
                      colnf%externaln_to_decomp_npools(c,j,l)            &
                       + colnf%dwt_livecrootn_to_cwdn(c,j)               &
                       + colnf%dwt_deadcrootn_to_cwdn(c,j)               &
                       + colnf%gap_mortality_n_to_cwdn(c,j)              &
                       + colnf%harvest_n_to_cwdn(c,j)                    &
                       + colnf%fire_mortality_n_to_cwdn(c,j)

               end if

               ! the following is the net changes of plant N to decompible N poools between time-step
               ! in pflotran, decomposible N pools increments ARE from previous time-step (saved above);
               ! while, in CLM-CN all plant N pools are updated with current N fluxes among plant and ground/soil.
               ! therefore, when do balance check it is needed to adjust the time-lag of changes.
               colnf%externaln_to_decomp_delta(c) =   &
                           colnf%externaln_to_decomp_delta(c) - &
                           colnf%externaln_to_decomp_npools(c,j,l)*dzsoi_decomp(j)

               if (abs(colnf%externaln_to_decomp_npools(c,j,l))<=1.e-21_r8) then
                   colnf%externaln_to_decomp_npools(c,j,l) = 0._r8
               end if
            end do !l = 1, ndecomp_pools
         end do !j = 1, nlevdecomp_full
      end do !fc = 1,num_soilc


      ! if pflotran hydrology NOT coupled, need to do:
      ! saving for (next time-step) possible including of RT mass-transfer in PFLOTRAN bgc coupling.
      ! (NOT USED anymore - 04/26/2017)
      if (.not. pf_hmode) then
         do j = 1, nlevdecomp_full
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               colnf%no3_net_transport_vr(c,j) = colnf%smin_no3_runoff_vr(c,j) + &
                                              colnf%smin_no3_leached_vr(c,j)
            end do
         end do
      else
         do j = 1, nlevdecomp_full
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               colnf%no3_net_transport_vr(c,j) = 0._r8
            end do
         end do
      end if

      ! change the sign so that it is the increments from the previous time-step (unit: g/m2/s)
      do fc = 1, num_soilc
         c = filter_soilc(fc)
         colnf%externaln_to_decomp_delta(c) = -colnf%externaln_to_decomp_delta(c)
      end do

    end subroutine colnf_summaryint_acc


    !-----------------------------------------------------------------------
    subroutine colpf_setvalues_acc( colpf, num_column, filter_column, value_column)
      !$acc routine seq
      ! !DESCRIPTION:
      ! Set phosphorus flux variables
      !
      use ColumnDataType, only : column_phosphorus_flux
      ! !ARGUMENTS:
      type(column_phosphorus_flux), intent(inout) :: colpf
      integer , intent(in) :: num_column
      integer , intent(in) :: filter_column(:)
      real(r8), intent(in) :: value_column
      !
      ! !LOCAL VARIABLES:
      integer :: fi,i,j,k,l     ! loop index
      !------------------------------------------------------------------------
      do j = 1, nlevdecomp_full
         do fi = 1,num_column
            i = filter_column(fi)

            ! phenology: litterfall and crop fluxes associated wit
            colpf%phenology_p_to_litr_met_p(i,j)        = value_column
            colpf%phenology_p_to_litr_cel_p(i,j)        = value_column
            colpf%phenology_p_to_litr_lig_p(i,j)        = value_column

            ! gap mortality
            colpf%gap_mortality_p_to_litr_met_p(i,j)    = value_column
            colpf%gap_mortality_p_to_litr_cel_p(i,j)    = value_column
            colpf%gap_mortality_p_to_litr_lig_p(i,j)    = value_column
            colpf%gap_mortality_p_to_cwdp(i,j)          = value_column

            ! fire
            colpf%fire_mortality_p_to_cwdp(i,j)         = value_column
            colpf%m_p_to_litr_met_fire(i,j)             = value_column
            colpf%m_p_to_litr_cel_fire(i,j)             = value_column
            colpf%m_p_to_litr_lig_fire(i,j)             = value_column

            ! harvest
            colpf%harvest_p_to_litr_met_p(i,j)          = value_column
            colpf%harvest_p_to_litr_cel_p(i,j)          = value_column
            colpf%harvest_p_to_litr_lig_p(i,j)          = value_column
            colpf%harvest_p_to_cwdp(i,j)                = value_column

            colpf%primp_to_labilep_vr(i,j)              = value_column
            colpf%labilep_to_secondp_vr(i,j)            = value_column
            colpf%secondp_to_labilep_vr(i,j)            = value_column
            colpf%secondp_to_occlp_vr(i,j)              = value_column

            colpf%sminp_leached_vr(i,j)                 = value_column

            colpf%potential_immob_p_vr(i,j)             = value_column
            colpf%actual_immob_p_vr(i,j)                = value_column
            colpf%sminp_to_plant_vr(i,j)                = value_column
            colpf%supplement_to_sminp_vr(i,j)           = value_column
            colpf%gross_pmin_vr(i,j)                    = value_column
            colpf%net_pmin_vr(i,j)                      = value_column
            colpf%biochem_pmin_vr(i,j)                  = value_column
            colpf%biochem_pmin_to_ecosysp_vr(i,j)       = value_column

            ! bgc interface & pflotran
            colpf%plant_pdemand_vr(i,j)                 = value_column
            colpf%adsorb_to_labilep_vr(i,j)             = value_column
            colpf%desorb_to_solutionp_vr(i,j)           = value_column
         end do
      end do

      do fi = 1,num_column
         i = filter_column(fi)

         colpf%pdep_to_sminp(i)             = value_column
         colpf%fert_p_to_sminp(i)           = value_column
         colpf%hrv_deadstemp_to_prod10p(i)  = value_column
         colpf%hrv_deadstemp_to_prod100p(i) = value_column
         colpf%hrv_cropp_to_prod1p(i)       = value_column
         colpf%prod10p_loss(i)              = value_column
         colpf%prod100p_loss(i)             = value_column
         colpf%product_ploss(i)             = value_column
         colpf%prod1p_loss(i)               = value_column
         colpf%potential_immob_p(i)         = value_column
         colpf%actual_immob_p(i)            = value_column
         colpf%sminp_to_plant(i)            = value_column
         colpf%supplement_to_sminp(i)       = value_column
         colpf%gross_pmin(i)                = value_column
         colpf%net_pmin(i)                  = value_column
         colpf%biochem_pmin(i)              = value_column
         colpf%primp_to_labilep(i)          = value_column
         colpf%labilep_to_secondp(i)        = value_column
         colpf%secondp_to_labilep(i)        = value_column
         colpf%secondp_to_occlp(i)          = value_column
         colpf%sminp_leached(i)             = value_column
         colpf%pinputs(i)                   = value_column
         colpf%poutputs(i)                  = value_column
         colpf%fire_ploss(i)                = value_column
         colpf%som_p_leached(i)             = value_column

         ! Zero p2c column fluxes
         colpf%fire_ploss(i)                = value_column
         colpf%wood_harvestp(i)             = value_column

         ! bgc-interface
         colpf%plant_pdemand(i)             = value_column

         colpf%fire_ploss(i)                = value_column
         colpf%wood_harvestp(i)             = value_column

         colpf%adsorb_to_labilep(i)         = value_column
         colpf%desorb_to_solutionp(i)       = value_column

      end do

      do k = 1, ndecomp_pools
         do fi = 1,num_column
            i = filter_column(fi)
            colpf%decomp_ppools_leached(i,k) = value_column
            colpf%m_decomp_ppools_to_fire(i,k) = value_column
         end do
      end do

      do k = 1, ndecomp_pools
         do j = 1, nlevdecomp_full
            do fi = 1,num_column
               i = filter_column(fi)
               colpf%m_decomp_ppools_to_fire_vr(i,j,k) = value_column
               colpf%decomp_ppools_transport_tendency(i,j,k) = value_column
            end do
         end do
      end do

      do l = 1, ndecomp_cascade_transitions
         do fi = 1,num_column
            i = filter_column(fi)
            colpf%decomp_cascade_ptransfer(i,l) = value_column
            colpf%decomp_cascade_sminp_flux(i,l) = value_column
         end do
      end do

      do l = 1, ndecomp_cascade_transitions
         do j = 1, nlevdecomp_full
            do fi = 1,num_column
               i = filter_column(fi)
               colpf%decomp_cascade_ptransfer_vr(i,j,l) = value_column
               colpf%decomp_cascade_sminp_flux_vr(i,j,l) = value_column
            end do
         end do
      end do

      do k = 1, ndecomp_pools
         do j = 1, nlevdecomp_full
            do fi = 1,num_column
               i = filter_column(fi)
               colpf%decomp_ppools_sourcesink(i,j,k) = value_column
               colpf%biochem_pmin_ppools_vr(i,j,k) = value_column
            end do
         end do
      end do

      ! pflotran
      do k = 1, ndecomp_pools
         do j = 1, nlevdecomp_full
            do fi = 1,num_column
               i = filter_column(fi)
               ! only initializing in the first time-step
               if ( colpf%externalp_to_decomp_ppools(i,j,k) == spval ) then
                  colpf%externalp_to_decomp_ppools(i,j,k) = value_column
               end if
            end do
         end do
      end do

      do j = 1, nlevdecomp_full
         do fi = 1,num_column
            i = filter_column(fi)
            ! only initializing in the first time-step
            if ( colpf%sminp_net_transport_vr(i,j) == spval ) then
               colpf%sminp_net_transport_vr(i,j) = value_column
            end if
         end do
      end do

      do fi = 1,num_column
         i = filter_column(fi)
         ! only initializing in the first time-step
         if ( colpf%externalp_to_decomp_delta(i) == spval ) then
            colpf%externalp_to_decomp_delta(i) = value_column
         end if
         if ( colpf%sminp_net_transport_delta(i) == spval ) then
            colpf%sminp_net_transport_delta(i)   = value_column
         end if
      end do

    end subroutine colpf_setvalues_acc

    !------------------------------------------------------------
    subroutine colcf_summary_for_ch4_acc( colcf, bounds, num_soilc, filter_soilc)
      !$acc routine seq
      ! !DESCRIPTION:
      ! summarize column-level fluxes for methane calculation
      !
      ! !USES:
      use ColumnDataType,   only : column_carbon_flux
      !
      ! !ARGUMENTS:
      type(column_carbon_flux), intent(inout)  :: colcf
      type(bounds_type), intent(in) :: bounds
      integer, intent(in) :: num_soilc
      integer, intent(in) :: filter_soilc(:)
      !
      ! !LOCAL VARIABLES
      integer :: fc, c
      integer :: j,k,l       ! indices
      !------------------------------------------------------------
      associate(&
           is_litter =>    decomp_cascade_con%is_litter , & ! Input:  [logical (:) ]  TRUE => pool is a litter pool
           is_soil   =>    decomp_cascade_con%is_soil   , & ! Input:  [logical (:) ]  TRUE => pool is a soil pool
           is_cwd    =>    decomp_cascade_con%is_cwd      &
      )
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         colcf%somhr(c)              = 0._r8
         colcf%lithr(c)              = 0._r8
         colcf%decomp_cascade_hr(c,1:ndecomp_cascade_transitions)= 0._r8
         if (.not. (use_pflotran .and. pf_cmode)) then
         ! pflotran has returned 'hr_vr(begc:endc,1:nlevdecomp)' to ALM before colcf subroutine is called in CNEcosystemDynNoLeaching2
         ! thus 'hr_vr_col' should NOT be set to 0
              colcf%hr_vr(c,1:nlevdecomp) = 0._r8
         end if
      enddo

      if ( (.not. is_active_betr_bgc           ) .and. &
           (.not. (use_pflotran .and. pf_cmode))) then
        ! vertically integrate HR and decomposition cascade fluxes
        do k = 1, ndecomp_cascade_transitions

         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               colcf%decomp_cascade_hr(c,k) = &
                  colcf%decomp_cascade_hr(c,k) + &
                  colcf%decomp_cascade_hr_vr(c,j,k) * dzsoi_decomp(j)

            end do
         end do
        end do

        ! litter heterotrophic respiration (LITHR)
        do k = 1, ndecomp_cascade_transitions
          if ( is_litter(decomp_cascade_con%cascade_donor_pool(k)) .or. is_cwd((decomp_cascade_con%cascade_donor_pool(k)))) then
            do fc = 1,num_soilc
              c = filter_soilc(fc)
              colcf%lithr(c) = &
                colcf%lithr(c) + &
                colcf%decomp_cascade_hr(c,k)
            end do
          end if
        end do

        ! soil organic matter heterotrophic respiration (SOMHR)
        do k = 1, ndecomp_cascade_transitions
          if ( is_soil(decomp_cascade_con%cascade_donor_pool(k)) ) then
            do fc = 1,num_soilc
              c = filter_soilc(fc)
              colcf%somhr(c) = &
                colcf%somhr(c) + &
                colcf%decomp_cascade_hr(c,k)
            end do
          end if
        end do

        ! total heterotrophic respiration, vertically resolved (HR)

        do k = 1, ndecomp_cascade_transitions
          do j = 1,nlevdecomp
            do fc = 1,num_soilc
              c = filter_soilc(fc)
              colcf%hr_vr(c,j) = &
                  colcf%hr_vr(c,j) + &
                  colcf%decomp_cascade_hr_vr(c,j,k)
            end do
          end do
        end do
      endif

      end associate

    end subroutine colcf_summary_for_ch4_acc

    !-----------------------------------------------------------------------
    subroutine colcf_summary_acc(colcf, bounds, num_soilc, filter_soilc, isotope, dtime)
      !
      ! !DESCRIPTION:
      ! column-level carbon flux summary calculations
      !
      !$acc routine seq
      ! !ARGUMENTS:
      type(column_carbon_flux)              :: colcf
      type(bounds_type)      , intent(in)    :: bounds
      integer                , intent(in)    :: num_soilc       ! number of soil columns in filter
      integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
      integer                , intent(in)    :: isotope  !c13 = 0, c14 = 1, bulk = 2
      real(r8)               , intent(in)    :: dtime
      !
      ! !LOCAL VARIABLES:
      real(r8) :: nfixlags ! temp variables for making lagged npp
      integer  :: c,p,j,k,l       ! indices
      integer  :: fc              ! lake filter indices
      real(r8) :: maxdepth        ! depth to integrate soil variables
      integer  :: nlev
      !-----------------------------------------------------------------------
      associate(&
           is_litter =>    decomp_cascade_con%is_litter , & ! Input:  [logical (:) ]  TRUE => pool is a litter pool
           is_soil   =>    decomp_cascade_con%is_soil   , & ! Input:  [logical (:) ]  TRUE => pool is a soil pool
           is_cwd    =>    decomp_cascade_con%is_cwd      & ! Input:  [logical (:) ]  TRUE => pool is a cwd pool
           )

      if (use_fates) return

      ! PET: retaining the following here during migration, but this is science code that should
      ! really be in the NDynamics module. Flag for relocation during ELM v2 code cleanup.
      if ( isotope == bulk) then
         if (nfix_timeconst > 0._r8 .and. nfix_timeconst < 500._r8 ) then

            ! this code is to calculate an exponentially-relaxed npp value for use in NDynamics code
            !#py dtime = get_step_size()
            nfixlags = nfix_timeconst * secspday

            do fc = 1,num_soilc
               c = filter_soilc(fc)
               if ( colcf%lag_npp(c) /= spval ) then
                  colcf%lag_npp(c) = &
                       colcf%lag_npp(c) * exp(-dtime/nfixlags) + &
                       colcf%npp(c) * (1._r8 - exp(-dtime/nfixlags))
               else
                  ! first timestep
                  colcf%lag_npp(c) = colcf%npp(c)
               endif
            end do
         endif
      endif
      nlev = nlevdecomp
      if (use_pflotran .and. pf_cmode) nlev = nlevdecomp_full

      ! some zeroing
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         colcf%cwdc_loss(c)          = 0._r8
         colcf%som_c_leached(c)      = 0._r8
      end do

      if ( (.not. is_active_betr_bgc           ) .and. &
           (.not. (use_pflotran .and. pf_cmode))) then

         ! vertically integrate HR and decomposition cascade fluxes
         do k = 1, ndecomp_cascade_transitions

         do j = 1,nlev
            do fc = 1,num_soilc
               c = filter_soilc(fc)

                  colcf%decomp_cascade_ctransfer(c,k) = &
                       colcf%decomp_cascade_ctransfer(c,k) + &
                       colcf%decomp_cascade_ctransfer_vr(c,j,k) * dzsoi_decomp(j)
               end do
            end do
         end do


         ! total heterotrophic respiration (HR)
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            colcf%hr(c) = &
                 colcf%lithr(c) + &
                 colcf%somhr(c)
         end do


      elseif (is_active_betr_bgc) then

         do fc = 1, num_soilc
            c = filter_soilc(fc)
            colcf%hr(c) = dot_product(colcf%hr_vr(c,1:nlevdecomp),dzsoi_decomp(1:nlevdecomp))
         enddo
      endif

      ! some zeroing
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         colcf%somhr(c)              = 0._r8
         colcf%lithr(c)              = 0._r8
         colcf%decomp_cascade_hr(c,1:ndecomp_cascade_transitions)= 0._r8
         if (.not. (use_pflotran .and. pf_cmode)) then
         ! pflotran has returned 'hr_vr(begc:endc,1:nlevdecomp)' to ALM before this subroutine is called in CNEcosystemDynNoLeaching2
         ! thus 'hr_vr_col' should NOT be set to 0
              colcf%hr_vr(c,1:nlevdecomp) = 0._r8
         end if
      enddo

      ! vertically integrate HR and decomposition cascade fluxes
      do k = 1, ndecomp_cascade_transitions

         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               colcf%decomp_cascade_hr(c,k) = &
                  colcf%decomp_cascade_hr(c,k) + &
                  colcf%decomp_cascade_hr_vr(c,j,k) * dzsoi_decomp(j)

            end do
         end do
      end do

      ! litter heterotrophic respiration (LITHR)
      do k = 1, ndecomp_cascade_transitions
         if ( is_litter(decomp_cascade_con%cascade_donor_pool(k)) .or. is_cwd((decomp_cascade_con%cascade_donor_pool(k)))) then
            do fc = 1,num_soilc
              c = filter_soilc(fc)
              colcf%lithr(c) = &
                colcf%lithr(c) + &
                colcf%decomp_cascade_hr(c,k)
            end do
         end if
      end do

      ! soil organic matter heterotrophic respiration (SOMHR)
      do k = 1, ndecomp_cascade_transitions
         if ( is_soil(decomp_cascade_con%cascade_donor_pool(k)) ) then
            do fc = 1,num_soilc
              c = filter_soilc(fc)
              colcf%somhr(c) = &
                colcf%somhr(c) + &
                colcf%decomp_cascade_hr(c,k)
            end do
         end if
      end do

      ! total heterotrophic respiration, vertically resolved (HR)

      do k = 1, ndecomp_cascade_transitions
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
              c = filter_soilc(fc)
              colcf%hr_vr(c,j) = &
                  colcf%hr_vr(c,j) + &
                  colcf%decomp_cascade_hr_vr(c,j,k)
            end do
         end do
      end do

      !----------------------------------------------------------------
      ! bgc interface & pflotran:
      !----------------------------------------------------------------
      if (use_clm_interface.and. (use_pflotran .and. pf_cmode)) then
          call colcf_summary_pf_acc(colcf, bounds, num_soilc, filter_soilc, dtime)
      end if
      !----------------------------------------------------------------

      do fc = 1,num_soilc
         c = filter_soilc(fc)
         ! total soil respiration, heterotrophic + root respiration (SR)
         colcf%sr(c) = &
              colcf%rr(c) + &
              colcf%hr(c)

         ! total ecosystem respiration, autotrophic + heterotrophic (ER)
         colcf%er(c) = &
              colcf%ar(c) + &
              colcf%hr(c)

         ! litter fire losses (LITFIRE)
         colcf%litfire(c) = 0._r8

         ! total product loss
         colcf%product_closs(c) = &
              colcf%prod10c_loss(c)  + &
              colcf%prod100c_loss(c) + &
              colcf%prod1c_loss(c)

         ! soil organic matter fire losses (SOMFIRE)
         colcf%somfire(c) = 0._r8

         ! total ecosystem fire losses (TOTFIRE)
         colcf%totfire(c) = &
              colcf%litfire(c) + &
              colcf%somfire(c) + &
              colcf%vegfire(c)
      end do

      ! vertically integrate column-level carbon fire losses
      do l = 1, ndecomp_pools
         do j = 1,nlev
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               colcf%m_decomp_cpools_to_fire(c,l) = &
                    colcf%m_decomp_cpools_to_fire(c,l) + &
                    colcf%m_decomp_cpools_to_fire_vr(c,j,l)*dzsoi_decomp(j)
            end do
         end do
      end do

      ! column-level carbon losses to fire, including pft losses
      do fc = 1,num_soilc
         c = filter_soilc(fc)

         colcf%fire_closs(c) = colcf%fire_closs_p2c(c)
         do l = 1, ndecomp_pools
            colcf%fire_closs(c) = &
                 colcf%fire_closs(c) + &
                 colcf%m_decomp_cpools_to_fire(c,l)
         end do

         ! column-level carbon losses due to landcover change
         colcf%dwt_closs(c) = &
              colcf%dwt_conv_cflux(c)

         ! net ecosystem production, excludes fire flux, landcover change, and loss from wood products, positive for sink (NEP)
         colcf%nep(c) = &
              colcf%gpp(c) - &
              colcf%er(c)

         ! net biome production of carbon, includes depletion from: fire flux, landcover change flux, and loss
         ! from wood products pools, positive for sink (NBP)
         colcf%nbp(c) =             &
              colcf%nep(c)        - &
              colcf%fire_closs(c) - &
              colcf%dwt_closs(c)  - &
              colcf%product_closs(c)

         ! net ecosystem exchange of carbon, includes fire flux, landcover change flux, loss
         ! from wood products pools, and hrv_xsmrpool flux, positive for source (NEE)
         colcf%nee(c) =                &
              -colcf%nep(c)           + &
              colcf%fire_closs(c)    + &
              colcf%dwt_closs(c)     + &
              colcf%product_closs(c) + &
              colcf%hrv_xsmrpool_to_atm(c)

         ! land use flux and land uptake
         colcf%landuseflux(c) = &
              colcf%dwt_closs(c) + &
              colcf%product_closs(c)

         colcf%landuptake(c) = &
              colcf%nee(c) - &
              colcf%landuseflux(c)
      end do

      if  (.not. is_active_betr_bgc) then

         ! (cWDC_HR) - coarse woody debris heterotrophic respiration
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            colcf%cwdc_hr(c) = 0._r8
         end do

         ! (cWDC_LOSS) - coarse woody debris C loss
         do l = 1, ndecomp_pools
            if ( is_cwd(l) ) then
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  colcf%cwdc_loss(c) = &
                       colcf%cwdc_loss(c) + &
                       colcf%m_decomp_cpools_to_fire(c,l)
               end do
            end if
         end do

         do k = 1, ndecomp_cascade_transitions
            if ( is_cwd(decomp_cascade_con%cascade_donor_pool(k)) ) then
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  colcf%cwdc_loss(c) = &
                       colcf%cwdc_loss(c) + &
                       colcf%decomp_cascade_ctransfer(c,k)
               end do
            end if
         end do

         if (.not.(use_pflotran .and. pf_cmode)) then
            ! (LITTERC_LOSS) - litter C loss
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               colcf%litterc_loss(c) = colcf%lithr(c)
            end do
         end if !(.not.(use_pflotran .and. pf_cmode))

         do l = 1, ndecomp_pools
            if ( is_litter(l) ) then
               do fc = 1,num_soilc
                   c = filter_soilc(fc)
                   colcf%litterc_loss(c) = &
                      colcf%litterc_loss(c) + &
                      colcf%m_decomp_cpools_to_fire(c,l)
               end do
            end if
         end do


         do k = 1, ndecomp_cascade_transitions
           if ( is_litter(decomp_cascade_con%cascade_donor_pool(k)) ) then
             do fc = 1,num_soilc
               c = filter_soilc(fc)
               colcf%litterc_loss(c) = &
                    colcf%litterc_loss(c) + &
                    colcf%decomp_cascade_ctransfer(c,k)
             end do
           end if
         end do

         if (use_pflotran .and. pf_cmode) then
            ! note: the follwoing should be useful to non-pflotran-coupled, but seems cause 1 BFB test unmatching.
            ! add up all vertical transport tendency terms and calculate total som leaching loss as the sum of these
            do l = 1, ndecomp_pools
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  colcf%decomp_cpools_leached(c,l) = 0._r8
               end do
               do j = 1, nlev
                  do fc = 1,num_soilc
                     c = filter_soilc(fc)
                     colcf%decomp_cpools_leached(c,l) = &
                       colcf%decomp_cpools_leached(c,l) + &
                       colcf%decomp_cpools_transport_tendency(c,j,l) * dzsoi_decomp(j)
                  end do
               end do
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  colcf%som_c_leached(c) = &
                     colcf%som_c_leached(c) + &
                     colcf%decomp_cpools_leached(c,l)
               end do
            end do
         end if

      end if ! .not. is_active_betr_bgc

      do fc = 1,num_soilc
          c = filter_soilc(fc)
          colcf%plant_to_litter_cflux(c) = 0._r8
          colcf%plant_to_cwd_cflux(c) = 0._r8
          do j = 1, nlev
              colcf%plant_to_litter_cflux(c) = &
                  colcf%plant_to_litter_cflux(c)  + &
                  colcf%phenology_c_to_litr_met_c(c,j)* dzsoi_decomp(j) + &
                  colcf%phenology_c_to_litr_cel_c(c,j)* dzsoi_decomp(j) + &
                  colcf%phenology_c_to_litr_lig_c(c,j)* dzsoi_decomp(j) + &
                  colcf%gap_mortality_c_to_litr_met_c(c,j)* dzsoi_decomp(j) + &
                  colcf%gap_mortality_c_to_litr_cel_c(c,j)* dzsoi_decomp(j) + &
                  colcf%gap_mortality_c_to_litr_lig_c(c,j)* dzsoi_decomp(j) + &
                  colcf%m_c_to_litr_met_fire(c,j)* dzsoi_decomp(j) + &
                  colcf%m_c_to_litr_cel_fire(c,j)* dzsoi_decomp(j) + &
                  colcf%m_c_to_litr_lig_fire(c,j)* dzsoi_decomp(j)
              colcf%plant_to_cwd_cflux(c) = &
                  colcf%plant_to_cwd_cflux(c) + &
                  colcf%gap_mortality_c_to_cwdc(c,j)* dzsoi_decomp(j) + &
                  colcf%fire_mortality_c_to_cwdc(c,j)* dzsoi_decomp(j)
          end do
      end do

      end associate

    end subroutine colcf_summary_acc

    !-------------------------------------------------------------------------------------------------
    subroutine colcf_summary_pf_acc(colcf, bounds, num_soilc, filter_soilc, dtime)
      !
      ! !DESCRIPTION:
      ! bgc interface & pflotran:
      ! On the radiation time step, perform column-level carbon
      ! summary calculations, which mainly from PFLOTRAN bgc
      !
      !$acc routine seq
      ! !ARGUMENTS:
      type(column_carbon_flux)       :: colcf
      type(bounds_type) ,  intent(in) :: bounds
      integer,             intent(in) :: num_soilc       ! number of soil columns in filter
      integer,             intent(in) :: filter_soilc(:) ! filter for soil columns
      real(r8), intent(in)            :: dtime                ! time-step (s)

      !
      ! !CALLED FROM:
      ! subroutine Summary (if plotran bgc coupled with CLM-CN
      !
      ! LOCAL VARIABLES:
      integer :: c,j,l                 ! indices
      integer :: fc                    ! column filter indices

      associate(&
          is_litter =>    decomp_cascade_con%is_litter , & ! Input:  [logical (:) ]  TRUE => pool is a litter pool
          is_soil   =>    decomp_cascade_con%is_soil   , & ! Input:  [logical (:) ]  TRUE => pool is a soil pool
          is_cwd    =>    decomp_cascade_con%is_cwd      & ! Input:  [logical (:) ]  TRUE => pool is a cwd pool
          )

      !#py dtime = get_step_size()
      ! total heterotrophic respiration (HR)
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         colcf%hr(c) = 0._r8
         do j = 1,nlevdecomp_full
            colcf%hr(c) = colcf%hr(c) + &
               colcf%hr_vr(c,j) * dzsoi_decomp(j)
         end do
      end do

      ! new variable to account for co2 exchange (not all HR goes to atm at current time-step)
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         colcf%f_co2_soil(c) = 0._r8
      end do
      do j = 1,nlevdecomp_full
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            colcf%f_co2_soil(c) = colcf%f_co2_soil(c) + &
               colcf%f_co2_soil_vr(c,j) * dzsoi_decomp(j)
         end do
      end do

      do fc = 1,num_soilc
         c = filter_soilc(fc)
         colcf%cwdc_hr(c)      = 0._r8
         colcf%cwdc_loss(c)    = 0._r8
         colcf%litterc_loss(c) = 0._r8
      end do

      do l = 1, ndecomp_pools
         if ( is_cwd(l) ) then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               do j = 1, nlevdecomp_full
                  colcf%cwdc_loss(c) = &
                     colcf%cwdc_loss(c) + &
                     colcf%decomp_cpools_sourcesink(c,j,l) / dtime
               end do
            end do
         end if

         if ( is_litter(l) ) then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               do j = 1, nlevdecomp_full
                  colcf%litterc_loss(c) = &
                     colcf%litterc_loss(c) + &
                     colcf%decomp_cpools_sourcesink(c,j,l) / dtime
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
          colcf%externalc_to_decomp_delta(c) = 0._r8
          do l = 1, ndecomp_pools
            do j = 1, nlevdecomp_full
              colcf%externalc_to_decomp_delta(c) = colcf%externalc_to_decomp_delta(c) + &
                                  colcf%externalc_to_decomp_cpools(c,j,l)*dzsoi_decomp(j)
            end do
          end do
      end do
      !
      ! do the initialization for the following variable here.
      ! DON'T do so in the beginning of CLM-CN time-step (otherwise the above saved will not work)

      do fc = 1,num_soilc
          c = filter_soilc(fc)
          colcf%externalc_to_decomp_cpools(c, 1:nlevdecomp_full, 1:ndecomp_pools) = 0._r8
      end do

      do fc = 1,num_soilc
         c = filter_soilc(fc)
         do l = 1, ndecomp_pools
            do j = 1, nlevdecomp_full
               ! for litter C pools
               if (l==i_met_lit) then
                  colcf%externalc_to_decomp_cpools(c,j,l) =                 &
                      colcf%externalc_to_decomp_cpools(c,j,l)               &
                          + colcf%phenology_c_to_litr_met_c(c,j)            &
                          + colcf%dwt_frootc_to_litr_met_c(c,j)             &
                          + colcf%gap_mortality_c_to_litr_met_c(c,j)        &
                          + colcf%harvest_c_to_litr_met_c(c,j)              &
                          + colcf%m_c_to_litr_met_fire(c,j)

               elseif (l==i_cel_lit) then
                  colcf%externalc_to_decomp_cpools(c,j,l) =                 &
                      colcf%externalc_to_decomp_cpools(c,j,l)               &
                          + colcf%phenology_c_to_litr_cel_c(c,j)            &
                          + colcf%dwt_frootc_to_litr_cel_c(c,j)             &
                          + colcf%gap_mortality_c_to_litr_cel_c(c,j)        &
                          + colcf%harvest_c_to_litr_cel_c(c,j)              &
                          + colcf%m_c_to_litr_cel_fire(c,j)

               elseif (l==i_lig_lit) then
                  colcf%externalc_to_decomp_cpools(c,j,l) =                 &
                      colcf%externalc_to_decomp_cpools(c,j,l)               &
                          + colcf%phenology_c_to_litr_lig_c(c,j)            &
                          + colcf%dwt_frootc_to_litr_lig_c(c,j)             &
                          + colcf%gap_mortality_c_to_litr_lig_c(c,j)        &
                          + colcf%harvest_c_to_litr_lig_c(c,j)              &
                          + colcf%m_c_to_litr_lig_fire(c,j)

               ! for cwd
               elseif (l==i_cwd) then
                  colcf%externalc_to_decomp_cpools(c,j,l) =                 &
                      colcf%externalc_to_decomp_cpools(c,j,l)               &
                          + colcf%dwt_livecrootc_to_cwdc(c,j)               &
                          + colcf%dwt_deadcrootc_to_cwdc(c,j)               &
                          + colcf%gap_mortality_c_to_cwdc(c,j)              &
                          + colcf%harvest_c_to_cwdc(c,j)                    &
                          + colcf%fire_mortality_c_to_cwdc(c,j)

               end if

               ! the following is the net changes of plant C to decompible C poools between time-step
               ! in pflotran, decomposible C pools increments ARE from previous time-step (saved above);
               ! while, in CLM-CN all plant C pools are updated with current C fluxes among plant and ground/soil.
               ! therefore, when do balance check it is needed to adjust the time-lag of changes.
               colcf%externalc_to_decomp_delta(c) = colcf%externalc_to_decomp_delta(c) - &
                                  colcf%externalc_to_decomp_cpools(c,j,l)*dzsoi_decomp(j)

               if (abs(colcf%externalc_to_decomp_cpools(c,j,l))<=1.e-20_r8) then
                   colcf%externalc_to_decomp_cpools(c,j,l) = 0._r8
               end if

            end do
         end do
      end do

      ! change the sign so that it is the increments from the previous time-step (unit: from g/m2/s)
      do fc = 1, num_soilc
         c = filter_soilc(fc)
         colcf%externalc_to_decomp_delta(c) = -colcf%externalc_to_decomp_delta(c)
      end do

      end associate

    end subroutine colcf_summary_pf_acc

    !-----------------------------------------------------------------------
    subroutine colps_setvalues_acc( colps , num_column, filter_column, value_column)
      !$acc routine seq
      ! !DESCRIPTION:
      ! Set phosphorus state variables, column-level
      !
      use ColumnDataType,  only : column_phosphorus_state
      ! !ARGUMENTS:
      type(column_phosphorus_state), intent(inout) :: colps
      integer , intent(in)            :: num_column
      integer , intent(in)            :: filter_column(:)
      real(r8), intent(in)            :: value_column
      !
      ! !LOCAL VARIABLES:
      integer :: fi,i     ! loop index
      integer :: j,k      ! indices
      !------------------------------------------------------------------------
      do fi = 1,num_column
         i = filter_column(fi)

         colps%sminp(i)       = value_column
         colps%solutionp(i)   = value_column
         colps%labilep(i)     = value_column
         colps%secondp(i)     = value_column
         colps%occlp(i)       = value_column
         colps%primp(i)       = value_column
         colps%ptrunc(i)      = value_column
         colps%cwdp(i)        = value_column
         colps%totlitp(i)     = value_column
         colps%totsomp(i)     = value_column
         colps%totecosysp(i)  = value_column
         colps%totcolp(i)     = value_column
         colps%totsomp_1m(i)  = value_column
         colps%totlitp_1m(i)  = value_column
      end do

      do j = 1,nlevdecomp_full
         do fi = 1,num_column
            i = filter_column(fi)
            colps%sminp_vr(i,j)       = value_column
            colps%solutionp_vr(i,j)   = value_column
            colps%labilep_vr(i,j)     = value_column
            colps%secondp_vr(i,j)     = value_column
            colps%occlp_vr(i,j)       = value_column
            colps%primp_vr(i,j)       = value_column
            colps%ptrunc_vr(i,j)      = value_column
         end do
      end do

      ! column and decomp_pools
      do k = 1, ndecomp_pools
         do fi = 1,num_column
            i = filter_column(fi)
            colps%decomp_ppools(i,k)    = value_column
            colps%decomp_ppools_1m(i,k) = value_column
         end do
      end do

      ! column levdecomp, and decomp_pools
      do j = 1,nlevdecomp_full
         do k = 1, ndecomp_pools
            do fi = 1,num_column
               i = filter_column(fi)
               colps%decomp_ppools_vr(i,j,k) = value_column
            end do
         end do
      end do

    end subroutine colps_setvalues_acc


    !-----------------------------------------------------------------------
    subroutine vegcf_setvalues_acc (vegcf, num_patch, filter_patch, value_patch)
      !$acc routine seq
      ! !DESCRIPTION:
      ! Set vegetation-level carbon fluxes
      !
      use VegetationDataType , only : vegetation_carbon_flux
      ! !ARGUMENTS:
      type(vegetation_carbon_flux), intent(inout) :: vegcf
      integer , intent(in) :: num_patch
      integer , intent(in) :: filter_patch(:)
      real(r8), intent(in) :: value_patch
      !
      ! !LOCAL VARIABLES:
      integer :: fi,i     ! loop index
      !------------------------------------------------------------------------

      if(.not.use_fates) then
         do fi = 1,num_patch
            i = filter_patch(fi)

            vegcf%m_leafc_to_litter(i)                   = value_patch
            vegcf%m_frootc_to_litter(i)                  = value_patch
            vegcf%m_leafc_storage_to_litter(i)           = value_patch
            vegcf%m_frootc_storage_to_litter(i)          = value_patch
            vegcf%m_livestemc_storage_to_litter(i)       = value_patch
            vegcf%m_deadstemc_storage_to_litter(i)       = value_patch
            vegcf%m_livecrootc_storage_to_litter(i)      = value_patch
            vegcf%m_deadcrootc_storage_to_litter(i)      = value_patch
            vegcf%m_leafc_xfer_to_litter(i)              = value_patch
            vegcf%m_frootc_xfer_to_litter(i)             = value_patch
            vegcf%m_livestemc_xfer_to_litter(i)          = value_patch
            vegcf%m_deadstemc_xfer_to_litter(i)          = value_patch
            vegcf%m_livecrootc_xfer_to_litter(i)         = value_patch
            vegcf%m_deadcrootc_xfer_to_litter(i)         = value_patch
            vegcf%m_livestemc_to_litter(i)               = value_patch
            vegcf%m_deadstemc_to_litter(i)               = value_patch
            vegcf%m_livecrootc_to_litter(i)              = value_patch
            vegcf%m_deadcrootc_to_litter(i)              = value_patch
            vegcf%m_gresp_storage_to_litter(i)           = value_patch
            vegcf%m_gresp_xfer_to_litter(i)              = value_patch
            vegcf%m_cpool_to_litter(i)                   = value_patch
            vegcf%hrv_leafc_to_litter(i)                 = value_patch
            vegcf%hrv_leafc_storage_to_litter(i)         = value_patch
            vegcf%hrv_leafc_xfer_to_litter(i)            = value_patch
            vegcf%hrv_frootc_to_litter(i)                = value_patch
            vegcf%hrv_frootc_storage_to_litter(i)        = value_patch
            vegcf%hrv_frootc_xfer_to_litter(i)           = value_patch
            vegcf%hrv_livestemc_to_litter(i)             = value_patch
            vegcf%hrv_livestemc_storage_to_litter(i)     = value_patch
            vegcf%hrv_livestemc_xfer_to_litter(i)        = value_patch
            vegcf%hrv_deadstemc_to_prod10c(i)            = value_patch
            vegcf%hrv_deadstemc_to_prod100c(i)           = value_patch
            vegcf%hrv_leafc_to_prod1c(i)                 = value_patch
            vegcf%hrv_livestemc_to_prod1c(i)             = value_patch
            vegcf%hrv_grainc_to_prod1c(i)                = value_patch
            vegcf%hrv_cropc_to_prod1c(i)                 = value_patch
            vegcf%hrv_deadstemc_storage_to_litter(i)     = value_patch
            vegcf%hrv_deadstemc_xfer_to_litter(i)        = value_patch
            vegcf%hrv_livecrootc_to_litter(i)            = value_patch
            vegcf%hrv_livecrootc_storage_to_litter(i)    = value_patch
            vegcf%hrv_livecrootc_xfer_to_litter(i)       = value_patch
            vegcf%hrv_deadcrootc_to_litter(i)            = value_patch
            vegcf%hrv_deadcrootc_storage_to_litter(i)    = value_patch
            vegcf%hrv_deadcrootc_xfer_to_litter(i)       = value_patch
            vegcf%hrv_gresp_storage_to_litter(i)         = value_patch
            vegcf%hrv_gresp_xfer_to_litter(i)            = value_patch
            vegcf%hrv_xsmrpool_to_atm(i)                 = value_patch
            vegcf%hrv_cpool_to_litter(i)                 = value_patch

            vegcf%m_leafc_to_fire(i)                     = value_patch
            vegcf%m_leafc_storage_to_fire(i)             = value_patch
            vegcf%m_leafc_xfer_to_fire(i)                = value_patch
            vegcf%m_livestemc_to_fire(i)                 = value_patch
            vegcf%m_livestemc_storage_to_fire(i)         = value_patch
            vegcf%m_livestemc_xfer_to_fire(i)            = value_patch
            vegcf%m_deadstemc_to_fire(i)                 = value_patch
            vegcf%m_deadstemc_storage_to_fire(i)         = value_patch
            vegcf%m_deadstemc_xfer_to_fire(i)            = value_patch
            vegcf%m_frootc_to_fire(i)                    = value_patch
            vegcf%m_frootc_storage_to_fire(i)            = value_patch
            vegcf%m_frootc_xfer_to_fire(i)               = value_patch
            vegcf%m_livecrootc_to_fire(i)                = value_patch
            vegcf%m_livecrootc_storage_to_fire(i)        = value_patch
            vegcf%m_livecrootc_xfer_to_fire(i)           = value_patch
            vegcf%m_deadcrootc_to_fire(i)                = value_patch
            vegcf%m_deadcrootc_storage_to_fire(i)        = value_patch
            vegcf%m_deadcrootc_xfer_to_fire(i)           = value_patch
            vegcf%m_gresp_storage_to_fire(i)             = value_patch
            vegcf%m_gresp_xfer_to_fire(i)                = value_patch
            vegcf%m_cpool_to_fire(i)                     = value_patch

            vegcf%m_leafc_to_litter_fire(i)              = value_patch
            vegcf%m_leafc_storage_to_litter_fire(i)      = value_patch
            vegcf%m_leafc_xfer_to_litter_fire(i)         = value_patch
            vegcf%m_livestemc_to_litter_fire(i)          = value_patch
            vegcf%m_livestemc_storage_to_litter_fire(i)  = value_patch
            vegcf%m_livestemc_xfer_to_litter_fire(i)     = value_patch
            vegcf%m_livestemc_to_deadstemc_fire(i)       = value_patch
            vegcf%m_deadstemc_to_litter_fire(i)          = value_patch
            vegcf%m_deadstemc_storage_to_litter_fire(i)  = value_patch
            vegcf%m_deadstemc_xfer_to_litter_fire(i)     = value_patch
            vegcf%m_frootc_to_litter_fire(i)             = value_patch
            vegcf%m_frootc_storage_to_litter_fire(i)     = value_patch
            vegcf%m_frootc_xfer_to_litter_fire(i)        = value_patch
            vegcf%m_livecrootc_to_litter_fire(i)         = value_patch
            vegcf%m_livecrootc_storage_to_litter_fire(i) = value_patch
            vegcf%m_livecrootc_xfer_to_litter_fire(i)    = value_patch
            vegcf%m_livecrootc_to_deadcrootc_fire(i)     = value_patch
            vegcf%m_deadcrootc_to_litter_fire(i)         = value_patch
            vegcf%m_deadcrootc_storage_to_litter_fire(i) = value_patch
            vegcf%m_deadcrootc_xfer_to_litter_fire(i)    = value_patch
            vegcf%m_gresp_storage_to_litter_fire(i)      = value_patch
            vegcf%m_gresp_xfer_to_litter_fire(i)         = value_patch
            vegcf%m_cpool_to_litter_fire(i)              = value_patch

            vegcf%leafc_xfer_to_leafc(i)                 = value_patch
            vegcf%frootc_xfer_to_frootc(i)               = value_patch
            vegcf%livestemc_xfer_to_livestemc(i)         = value_patch
            vegcf%deadstemc_xfer_to_deadstemc(i)         = value_patch
            vegcf%livecrootc_xfer_to_livecrootc(i)       = value_patch
            vegcf%deadcrootc_xfer_to_deadcrootc(i)       = value_patch
            vegcf%leafc_to_litter(i)                     = value_patch
            vegcf%frootc_to_litter(i)                    = value_patch
            vegcf%leaf_mr(i)                             = value_patch
            vegcf%froot_mr(i)                            = value_patch
            vegcf%livestem_mr(i)                         = value_patch
            vegcf%livecroot_mr(i)                        = value_patch
            vegcf%grain_mr(i)                            = value_patch
            vegcf%leaf_curmr(i)                          = value_patch
            vegcf%froot_curmr(i)                         = value_patch
            vegcf%livestem_curmr(i)                      = value_patch
            vegcf%livecroot_curmr(i)                     = value_patch
            vegcf%grain_curmr(i)                         = value_patch
            vegcf%leaf_xsmr(i)                           = value_patch
            vegcf%froot_xsmr(i)                          = value_patch
            vegcf%livestem_xsmr(i)                       = value_patch
            vegcf%livecroot_xsmr(i)                      = value_patch
            vegcf%grain_xsmr(i)                          = value_patch
            vegcf%xr(i)                                  = value_patch
            vegcf%psnsun_to_cpool(i)                     = value_patch
            vegcf%psnshade_to_cpool(i)                   = value_patch
            vegcf%cpool_to_xsmrpool(i)                   = value_patch
            vegcf%cpool_to_leafc(i)                      = value_patch
            vegcf%cpool_to_leafc_storage(i)              = value_patch
            vegcf%cpool_to_frootc(i)                     = value_patch
            vegcf%cpool_to_frootc_storage(i)             = value_patch
            vegcf%cpool_to_livestemc(i)                  = value_patch
            vegcf%cpool_to_livestemc_storage(i)          = value_patch
            vegcf%cpool_to_deadstemc(i)                  = value_patch
            vegcf%cpool_to_deadstemc_storage(i)          = value_patch
            vegcf%cpool_to_livecrootc(i)                 = value_patch
            vegcf%cpool_to_livecrootc_storage(i)         = value_patch
            vegcf%cpool_to_deadcrootc(i)                 = value_patch
            vegcf%cpool_to_deadcrootc_storage(i)         = value_patch
            vegcf%cpool_to_gresp_storage(i)              = value_patch
            vegcf%cpool_leaf_gr(i)                       = value_patch
            vegcf%cpool_leaf_storage_gr(i)               = value_patch
            vegcf%transfer_leaf_gr(i)                    = value_patch
            vegcf%cpool_froot_gr(i)                      = value_patch
            vegcf%cpool_froot_storage_gr(i)              = value_patch
            vegcf%transfer_froot_gr(i)                   = value_patch
            vegcf%cpool_livestem_gr(i)                   = value_patch
            vegcf%cpool_livestem_storage_gr(i)           = value_patch
            vegcf%transfer_livestem_gr(i)                = value_patch
            vegcf%cpool_deadstem_gr(i)                   = value_patch
            vegcf%cpool_deadstem_storage_gr(i)           = value_patch
            vegcf%transfer_deadstem_gr(i)                = value_patch
            vegcf%cpool_livecroot_gr(i)                  = value_patch
            vegcf%cpool_livecroot_storage_gr(i)          = value_patch
            vegcf%transfer_livecroot_gr(i)               = value_patch
            vegcf%cpool_deadcroot_gr(i)                  = value_patch
            vegcf%cpool_deadcroot_storage_gr(i)          = value_patch
            vegcf%transfer_deadcroot_gr(i)               = value_patch
            vegcf%leafc_storage_to_xfer(i)               = value_patch
            vegcf%frootc_storage_to_xfer(i)              = value_patch
            vegcf%livestemc_storage_to_xfer(i)           = value_patch
            vegcf%deadstemc_storage_to_xfer(i)           = value_patch
            vegcf%livecrootc_storage_to_xfer(i)          = value_patch
            vegcf%deadcrootc_storage_to_xfer(i)          = value_patch
            vegcf%gresp_storage_to_xfer(i)               = value_patch
            vegcf%livestemc_to_deadstemc(i)              = value_patch
            vegcf%livecrootc_to_deadcrootc(i)            = value_patch
            vegcf%gpp(i)                                 = value_patch
            vegcf%gpp_before_downreg(i)                  = value_patch
            vegcf%mr(i)                                  = value_patch
            vegcf%current_gr(i)                          = value_patch
            vegcf%transfer_gr(i)                         = value_patch
            vegcf%storage_gr(i)                          = value_patch
            vegcf%gr(i)                                  = value_patch
            vegcf%ar(i)                                  = value_patch
            vegcf%rr(i)                                  = value_patch
            vegcf%npp(i)                                 = value_patch
            vegcf%agnpp(i)                               = value_patch
            vegcf%bgnpp(i)                               = value_patch
            vegcf%agwdnpp(i)                               = value_patch
            vegcf%litfall(i)                             = value_patch
            vegcf%vegfire(i)                             = value_patch
            vegcf%wood_harvestc(i)                       = value_patch
            vegcf%cinputs(i)                             = value_patch
            vegcf%coutputs(i)                            = value_patch
            vegcf%fire_closs(i)                          = value_patch
            vegcf%frootc_alloc(i)                        = value_patch
            vegcf%frootc_loss(i)                         = value_patch
            vegcf%leafc_alloc(i)                         = value_patch
            vegcf%leafc_loss(i)                          = value_patch
            vegcf%woodc_alloc(i)                         = value_patch
            vegcf%woodc_loss(i)                          = value_patch
            vegcf%xsmrpool_turnover(i)                   = value_patch
         end do
      end if !(.not.use_fates)

      if ( crop_prog )then
         do fi = 1,num_patch
            i = filter_patch(fi)
            vegcf%xsmrpool_to_atm(i)         = value_patch
            vegcf%livestemc_to_litter(i)     = value_patch
            vegcf%grainc_to_food(i)          = value_patch
            vegcf%grainc_xfer_to_grainc(i)   = value_patch
            vegcf%cpool_to_grainc(i)         = value_patch
            vegcf%cpool_to_grainc_storage(i) = value_patch
            vegcf%cpool_grain_gr(i)          = value_patch
            vegcf%cpool_grain_storage_gr(i)  = value_patch
            vegcf%transfer_grain_gr(i)       = value_patch
            vegcf%grainc_storage_to_xfer(i)  = value_patch
            vegcf%crop_seedc_to_leaf(i)      = value_patch
         end do
      end if

    end subroutine vegcf_setvalues_acc

    !-----------------------------------------------------------------------
    subroutine vegnf_setvalues_acc ( vegnf, num_patch, filter_patch, value_patch)
      !$acc routine seq
      ! !DESCRIPTION:
      ! Set vegetation-level nitrogen fluxes
      !
      use VegetationDataType,  only : vegetation_nitrogen_flux
      ! !ARGUMENTS:
      type (vegetation_nitrogen_flux), intent(inout) :: vegnf
      integer , intent(in)             :: num_patch
      integer , intent(in)             :: filter_patch(:)
      real(r8), intent(in)             :: value_patch
      !
      ! !LOCAL VARIABLES:
      integer :: fi,i     ! loop index
      !------------------------------------------------------------------------

      do fi = 1,num_patch
         i=filter_patch(fi)

         vegnf%m_leafn_to_litter(i)                   = value_patch
         vegnf%m_frootn_to_litter(i)                  = value_patch
         vegnf%m_leafn_storage_to_litter(i)           = value_patch
         vegnf%m_frootn_storage_to_litter(i)          = value_patch
         vegnf%m_livestemn_storage_to_litter(i)       = value_patch
         vegnf%m_deadstemn_storage_to_litter(i)       = value_patch
         vegnf%m_livecrootn_storage_to_litter(i)      = value_patch
         vegnf%m_deadcrootn_storage_to_litter(i)      = value_patch
         vegnf%m_leafn_xfer_to_litter(i)              = value_patch
         vegnf%m_frootn_xfer_to_litter(i)             = value_patch
         vegnf%m_livestemn_xfer_to_litter(i)          = value_patch
         vegnf%m_deadstemn_xfer_to_litter(i)          = value_patch
         vegnf%m_livecrootn_xfer_to_litter(i)         = value_patch
         vegnf%m_deadcrootn_xfer_to_litter(i)         = value_patch
         vegnf%m_livestemn_to_litter(i)               = value_patch
         vegnf%m_deadstemn_to_litter(i)               = value_patch
         vegnf%m_livecrootn_to_litter(i)              = value_patch
         vegnf%m_deadcrootn_to_litter(i)              = value_patch
         vegnf%m_retransn_to_litter(i)                = value_patch
         vegnf%m_npool_to_litter(i)                   = value_patch
         vegnf%hrv_leafn_to_litter(i)                 = value_patch
         vegnf%hrv_frootn_to_litter(i)                = value_patch
         vegnf%hrv_leafn_storage_to_litter(i)         = value_patch
         vegnf%hrv_frootn_storage_to_litter(i)        = value_patch
         vegnf%hrv_livestemn_storage_to_litter(i)     = value_patch
         vegnf%hrv_deadstemn_storage_to_litter(i)     = value_patch
         vegnf%hrv_livecrootn_storage_to_litter(i)    = value_patch
         vegnf%hrv_deadcrootn_storage_to_litter(i)    = value_patch
         vegnf%hrv_leafn_xfer_to_litter(i)            = value_patch
         vegnf%hrv_frootn_xfer_to_litter(i)           = value_patch
         vegnf%hrv_livestemn_xfer_to_litter(i)        = value_patch
         vegnf%hrv_deadstemn_xfer_to_litter(i)        = value_patch
         vegnf%hrv_livecrootn_xfer_to_litter(i)       = value_patch
         vegnf%hrv_deadcrootn_xfer_to_litter(i)       = value_patch
         vegnf%hrv_livestemn_to_litter(i)             = value_patch
         vegnf%hrv_deadstemn_to_prod10n(i)            = value_patch
         vegnf%hrv_deadstemn_to_prod100n(i)           = value_patch
         vegnf%hrv_leafn_to_prod1n(i)                 = value_patch
         vegnf%hrv_livestemn_to_prod1n(i)             = value_patch
         vegnf%hrv_grainn_to_prod1n(i)                = value_patch
         vegnf%hrv_cropn_to_prod1n(i)                 = value_patch
         vegnf%hrv_livecrootn_to_litter(i)            = value_patch
         vegnf%hrv_deadcrootn_to_litter(i)            = value_patch
         vegnf%hrv_retransn_to_litter(i)              = value_patch
         vegnf%hrv_npool_to_litter(i)                 = value_patch

         vegnf%m_leafn_to_fire(i)                     = value_patch
         vegnf%m_leafn_storage_to_fire(i)             = value_patch
         vegnf%m_leafn_xfer_to_fire(i)                = value_patch
         vegnf%m_livestemn_to_fire(i)                 = value_patch
         vegnf%m_livestemn_storage_to_fire(i)         = value_patch
         vegnf%m_livestemn_xfer_to_fire(i)            = value_patch
         vegnf%m_deadstemn_to_fire(i)                 = value_patch
         vegnf%m_deadstemn_storage_to_fire(i)         = value_patch
         vegnf%m_deadstemn_xfer_to_fire(i)            = value_patch
         vegnf%m_frootn_to_fire(i)                    = value_patch
         vegnf%m_frootn_storage_to_fire(i)            = value_patch
         vegnf%m_frootn_xfer_to_fire(i)               = value_patch
         vegnf%m_livecrootn_to_fire(i)                = value_patch
         vegnf%m_livecrootn_storage_to_fire(i)        = value_patch
         vegnf%m_livecrootn_xfer_to_fire(i)           = value_patch
         vegnf%m_deadcrootn_to_fire(i)                = value_patch
         vegnf%m_deadcrootn_storage_to_fire(i)        = value_patch
         vegnf%m_deadcrootn_xfer_to_fire(i)           = value_patch
         vegnf%m_retransn_to_fire(i)                  = value_patch
         vegnf%m_npool_to_fire(i)                     = value_patch

         vegnf%m_leafn_to_litter_fire(i)              = value_patch
         vegnf%m_leafn_storage_to_litter_fire(i)      = value_patch
         vegnf%m_leafn_xfer_to_litter_fire(i)         = value_patch
         vegnf%m_livestemn_to_litter_fire(i)          = value_patch
         vegnf%m_livestemn_storage_to_litter_fire(i)  = value_patch
         vegnf%m_livestemn_xfer_to_litter_fire(i)     = value_patch
         vegnf%m_livestemn_to_deadstemn_fire(i)       = value_patch
         vegnf%m_deadstemn_to_litter_fire(i)          = value_patch
         vegnf%m_deadstemn_storage_to_litter_fire(i)  = value_patch
         vegnf%m_deadstemn_xfer_to_litter_fire(i)     = value_patch
         vegnf%m_frootn_to_litter_fire(i)             = value_patch
         vegnf%m_frootn_storage_to_litter_fire(i)     = value_patch
         vegnf%m_frootn_xfer_to_litter_fire(i)        = value_patch
         vegnf%m_livecrootn_to_litter_fire(i)         = value_patch
         vegnf%m_livecrootn_storage_to_litter_fire(i) = value_patch
         vegnf%m_livecrootn_xfer_to_litter_fire(i)    = value_patch
         vegnf%m_livecrootn_to_deadcrootn_fire(i)     = value_patch
         vegnf%m_deadcrootn_to_litter_fire(i)         = value_patch
         vegnf%m_deadcrootn_storage_to_litter_fire(i) = value_patch
         vegnf%m_deadcrootn_xfer_to_litter_fire(i)    = value_patch
         vegnf%m_retransn_to_litter_fire(i)           = value_patch
         vegnf%m_npool_to_litter_fire(i)              = value_patch

         vegnf%leafn_xfer_to_leafn(i)                 = value_patch
         vegnf%frootn_xfer_to_frootn(i)               = value_patch
         vegnf%livestemn_xfer_to_livestemn(i)         = value_patch
         vegnf%deadstemn_xfer_to_deadstemn(i)         = value_patch
         vegnf%livecrootn_xfer_to_livecrootn(i)       = value_patch
         vegnf%deadcrootn_xfer_to_deadcrootn(i)       = value_patch
         vegnf%leafn_to_litter(i)                     = value_patch
         vegnf%leafn_to_retransn(i)                   = value_patch
         vegnf%frootn_to_litter(i)                    = value_patch
         vegnf%retransn_to_npool(i)                   = value_patch
         vegnf%sminn_to_npool(i)                      = value_patch
         vegnf%npool_to_leafn(i)                      = value_patch
         vegnf%npool_to_leafn_storage(i)              = value_patch
         vegnf%npool_to_frootn(i)                     = value_patch
         vegnf%npool_to_frootn_storage(i)             = value_patch
         vegnf%npool_to_livestemn(i)                  = value_patch
         vegnf%npool_to_livestemn_storage(i)          = value_patch
         vegnf%npool_to_deadstemn(i)                  = value_patch
         vegnf%npool_to_deadstemn_storage(i)          = value_patch
         vegnf%npool_to_livecrootn(i)                 = value_patch
         vegnf%npool_to_livecrootn_storage(i)         = value_patch
         vegnf%npool_to_deadcrootn(i)                 = value_patch
         vegnf%npool_to_deadcrootn_storage(i)         = value_patch
         vegnf%leafn_storage_to_xfer(i)               = value_patch
         vegnf%frootn_storage_to_xfer(i)              = value_patch
         vegnf%livestemn_storage_to_xfer(i)           = value_patch
         vegnf%deadstemn_storage_to_xfer(i)           = value_patch
         vegnf%livecrootn_storage_to_xfer(i)          = value_patch
         vegnf%deadcrootn_storage_to_xfer(i)          = value_patch
         vegnf%livestemn_to_deadstemn(i)              = value_patch
         vegnf%livestemn_to_retransn(i)               = value_patch
         vegnf%livecrootn_to_deadcrootn(i)            = value_patch
         vegnf%livecrootn_to_retransn(i)              = value_patch
         vegnf%ndeploy(i)                             = value_patch
         vegnf%ninputs(i)                             = value_patch
         vegnf%noutputs(i)                            = value_patch
         vegnf%wood_harvestn(i)                       = value_patch
         vegnf%fire_nloss(i)                          = value_patch
         vegnf%nfix_to_plantn(i)                      = value_patch
         vegnf%gap_nloss_litter(i)                    = value_patch
         vegnf%fire_nloss_litter(i)                   = value_patch
         vegnf%hrv_nloss_litter(i)                    = value_patch
         vegnf%sen_nloss_litter(i)                    = value_patch
         vegnf%crop_seedn_to_leaf(i)                  = value_patch
         vegnf%livestemn_to_litter(i)                 = value_patch
      end do

      if ( crop_prog )then
         do fi = 1,num_patch
            i = filter_patch(fi)
            vegnf%grainn_to_food(i)                   = value_patch
            vegnf%grainn_xfer_to_grainn(i)            = value_patch
            vegnf%npool_to_grainn(i)                  = value_patch
            vegnf%npool_to_grainn_storage(i)          = value_patch
            vegnf%grainn_storage_to_xfer(i)           = value_patch
            vegnf%soyfixn(i)                          = value_patch
            vegnf%frootn_to_retransn(i)               = value_patch
         end do
      end if

    end subroutine vegnf_setvalues_acc
    !-----------------------------------------------------------------------
    subroutine vegpf_setvalues_acc( vegpf, num_patch, filter_patch, value_patch)
      !$acc routine seq
      ! !DESCRIPTION:
      ! Set phosphorus flux variables
      !
      use VegetationDataType, only : vegetation_phosphorus_flux
      ! !ARGUMENTS:
      type(vegetation_phosphorus_flux), intent(inout) :: vegpf
      integer , intent(in) :: num_patch
      integer , intent(in) :: filter_patch(:)
      real(r8), intent(in) :: value_patch
      !
      ! !LOCAL VARIABLES:
      integer :: fi,i     ! loop index
      !------------------------------------------------------------------------
      do fi = 1,num_patch
         i=filter_patch(fi)

         vegpf%m_leafp_to_litter(i)                   = value_patch
         vegpf%m_frootp_to_litter(i)                  = value_patch
         vegpf%m_leafp_storage_to_litter(i)           = value_patch
         vegpf%m_frootp_storage_to_litter(i)          = value_patch
         vegpf%m_livestemp_storage_to_litter(i)       = value_patch
         vegpf%m_deadstemp_storage_to_litter(i)       = value_patch
         vegpf%m_livecrootp_storage_to_litter(i)      = value_patch
         vegpf%m_deadcrootp_storage_to_litter(i)      = value_patch
         vegpf%m_leafp_xfer_to_litter(i)              = value_patch
         vegpf%m_frootp_xfer_to_litter(i)             = value_patch
         vegpf%m_livestemp_xfer_to_litter(i)          = value_patch
         vegpf%m_deadstemp_xfer_to_litter(i)          = value_patch
         vegpf%m_livecrootp_xfer_to_litter(i)         = value_patch
         vegpf%m_deadcrootp_xfer_to_litter(i)         = value_patch
         vegpf%m_livestemp_to_litter(i)               = value_patch
         vegpf%m_deadstemp_to_litter(i)               = value_patch
         vegpf%m_livecrootp_to_litter(i)              = value_patch
         vegpf%m_deadcrootp_to_litter(i)              = value_patch
         vegpf%m_retransp_to_litter(i)                = value_patch
         vegpf%m_ppool_to_litter(i)                   = value_patch
         vegpf%hrv_leafp_to_litter(i)                 = value_patch
         vegpf%hrv_frootp_to_litter(i)                = value_patch
         vegpf%hrv_leafp_storage_to_litter(i)         = value_patch
         vegpf%hrv_frootp_storage_to_litter(i)        = value_patch
         vegpf%hrv_livestemp_storage_to_litter(i)     = value_patch
         vegpf%hrv_deadstemp_storage_to_litter(i)     = value_patch
         vegpf%hrv_livecrootp_storage_to_litter(i)    = value_patch
         vegpf%hrv_deadcrootp_storage_to_litter(i)    = value_patch
         vegpf%hrv_leafp_xfer_to_litter(i)            = value_patch
         vegpf%hrv_frootp_xfer_to_litter(i)           = value_patch
         vegpf%hrv_livestemp_xfer_to_litter(i)        = value_patch
         vegpf%hrv_deadstemp_xfer_to_litter(i)        = value_patch
         vegpf%hrv_livecrootp_xfer_to_litter(i)       = value_patch
         vegpf%hrv_deadcrootp_xfer_to_litter(i)       = value_patch
         vegpf%hrv_livestemp_to_litter(i)             = value_patch
         vegpf%hrv_deadstemp_to_prod10p(i)            = value_patch
         vegpf%hrv_deadstemp_to_prod100p(i)           = value_patch
         vegpf%hrv_leafp_to_prod1p(i)                 = value_patch
         vegpf%hrv_livestemp_to_prod1p(i)             = value_patch
         vegpf%hrv_grainp_to_prod1p(i)                = value_patch
         vegpf%hrv_cropp_to_prod1p(i)                 = value_patch
         vegpf%hrv_livecrootp_to_litter(i)            = value_patch
         vegpf%hrv_deadcrootp_to_litter(i)            = value_patch
         vegpf%hrv_retransp_to_litter(i)              = value_patch
         vegpf%hrv_ppool_to_litter(i)                 = value_patch

         vegpf%m_leafp_to_fire(i)                     = value_patch
         vegpf%m_leafp_storage_to_fire(i)             = value_patch
         vegpf%m_leafp_xfer_to_fire(i)                = value_patch
         vegpf%m_livestemp_to_fire(i)                 = value_patch
         vegpf%m_livestemp_storage_to_fire(i)         = value_patch
         vegpf%m_livestemp_xfer_to_fire(i)            = value_patch
         vegpf%m_deadstemp_to_fire(i)                 = value_patch
         vegpf%m_deadstemp_storage_to_fire(i)         = value_patch
         vegpf%m_deadstemp_xfer_to_fire(i)            = value_patch
         vegpf%m_frootp_to_fire(i)                    = value_patch
         vegpf%m_frootp_storage_to_fire(i)            = value_patch
         vegpf%m_frootp_xfer_to_fire(i)               = value_patch
         vegpf%m_livecrootp_to_fire(i)                = value_patch
         vegpf%m_livecrootp_storage_to_fire(i)        = value_patch
         vegpf%m_livecrootp_xfer_to_fire(i)           = value_patch
         vegpf%m_deadcrootp_to_fire(i)                = value_patch
         vegpf%m_deadcrootp_storage_to_fire(i)        = value_patch
         vegpf%m_deadcrootp_xfer_to_fire(i)           = value_patch
         vegpf%m_retransp_to_fire(i)                  = value_patch
         vegpf%m_ppool_to_fire(i)                     = value_patch

         vegpf%m_leafp_to_litter_fire(i)              = value_patch
         vegpf%m_leafp_storage_to_litter_fire(i)      = value_patch
         vegpf%m_leafp_xfer_to_litter_fire(i)         = value_patch
         vegpf%m_livestemp_to_litter_fire(i)          = value_patch
         vegpf%m_livestemp_storage_to_litter_fire(i)  = value_patch
         vegpf%m_livestemp_xfer_to_litter_fire(i)     = value_patch
         vegpf%m_livestemp_to_deadstemp_fire(i)       = value_patch
         vegpf%m_deadstemp_to_litter_fire(i)          = value_patch
         vegpf%m_deadstemp_storage_to_litter_fire(i)  = value_patch
         vegpf%m_deadstemp_xfer_to_litter_fire(i)     = value_patch
         vegpf%m_frootp_to_litter_fire(i)             = value_patch
         vegpf%m_frootp_storage_to_litter_fire(i)     = value_patch
         vegpf%m_frootp_xfer_to_litter_fire(i)        = value_patch
         vegpf%m_livecrootp_to_litter_fire(i)         = value_patch
         vegpf%m_livecrootp_storage_to_litter_fire(i) = value_patch
         vegpf%m_livecrootp_xfer_to_litter_fire(i)    = value_patch
         vegpf%m_livecrootp_to_deadcrootp_fire(i)     = value_patch
         vegpf%m_deadcrootp_to_litter_fire(i)         = value_patch
         vegpf%m_deadcrootp_storage_to_litter_fire(i) = value_patch
         vegpf%m_deadcrootp_xfer_to_litter_fire(i)    = value_patch
         vegpf%m_retransp_to_litter_fire(i)           = value_patch
         vegpf%m_ppool_to_litter_fire(i)              = value_patch

         vegpf%leafp_xfer_to_leafp(i)                 = value_patch
         vegpf%frootp_xfer_to_frootp(i)               = value_patch
         vegpf%livestemp_xfer_to_livestemp(i)         = value_patch
         vegpf%deadstemp_xfer_to_deadstemp(i)         = value_patch
         vegpf%livecrootp_xfer_to_livecrootp(i)       = value_patch
         vegpf%deadcrootp_xfer_to_deadcrootp(i)       = value_patch
         vegpf%leafp_to_litter(i)                     = value_patch
         vegpf%leafp_to_retransp(i)                   = value_patch
         vegpf%frootp_to_litter(i)                    = value_patch
         vegpf%retransp_to_ppool(i)                   = value_patch
         vegpf%sminp_to_ppool(i)                      = value_patch
         vegpf%ppool_to_leafp(i)                      = value_patch
         vegpf%ppool_to_leafp_storage(i)              = value_patch
         vegpf%ppool_to_frootp(i)                     = value_patch
         vegpf%ppool_to_frootp_storage(i)             = value_patch
         vegpf%ppool_to_livestemp(i)                  = value_patch
         vegpf%ppool_to_livestemp_storage(i)          = value_patch
         vegpf%ppool_to_deadstemp(i)                  = value_patch
         vegpf%ppool_to_deadstemp_storage(i)          = value_patch
         vegpf%ppool_to_livecrootp(i)                 = value_patch
         vegpf%ppool_to_livecrootp_storage(i)         = value_patch
         vegpf%ppool_to_deadcrootp(i)                 = value_patch
         vegpf%ppool_to_deadcrootp_storage(i)         = value_patch
         vegpf%leafp_storage_to_xfer(i)               = value_patch
         vegpf%frootp_storage_to_xfer(i)              = value_patch
         vegpf%livestemp_storage_to_xfer(i)           = value_patch
         vegpf%deadstemp_storage_to_xfer(i)           = value_patch
         vegpf%livecrootp_storage_to_xfer(i)          = value_patch
         vegpf%deadcrootp_storage_to_xfer(i)          = value_patch
         vegpf%livestemp_to_deadstemp(i)              = value_patch
         vegpf%livestemp_to_retransp(i)               = value_patch
         vegpf%livecrootp_to_deadcrootp(i)            = value_patch
         vegpf%livecrootp_to_retransp(i)              = value_patch
         vegpf%pdeploy(i)                             = value_patch
         vegpf%pinputs(i)                             = value_patch
         vegpf%poutputs(i)                            = value_patch
         vegpf%wood_harvestp(i)                       = value_patch
         vegpf%fire_ploss(i)                          = value_patch
         vegpf%biochem_pmin_to_plant(i)               = value_patch
         vegpf%gap_ploss_litter(i)                    = value_patch
         vegpf%fire_ploss_litter(i)                   = value_patch
         vegpf%hrv_ploss_litter(i)                    = value_patch
         vegpf%sen_ploss_litter(i)                    = value_patch
         vegpf%livestemp_to_litter(i)                 = value_patch
      end do

      if ( crop_prog )then
         do fi = 1,num_patch
            i = filter_patch(fi)
            vegpf%grainp_to_food(i)                   = value_patch
            vegpf%grainp_xfer_to_grainp(i)            = value_patch
            vegpf%ppool_to_grainp(i)                  = value_patch
            vegpf%ppool_to_grainp_storage(i)          = value_patch
            vegpf%grainp_storage_to_xfer(i)           = value_patch
            vegpf%frootp_to_retransp(i)               = value_patch
            vegpf%crop_seedp_to_leaf(i)               = value_patch
         end do
      end if

    end subroutine vegpf_setvalues_acc

    !------------------------------------------------------------
    subroutine vegcf_summary_rr_acc(vegcf, bounds, num_soilp, filter_soilp, &
      num_soilc, filter_soilc, col_cf_input)
      !$acc routine seq
      ! !DESCRIPTION:
      ! summarize root respiration
      !
      ! !USES:
      !
      use VegetationDataType, only : vegetation_carbon_flux
      use ColumnDataType, only : column_carbon_flux
      ! !ARGUMENTS:
      type(vegetation_carbon_flux), intent(inout) :: vegcf
      type(bounds_type), intent(in) :: bounds
      integer, intent(in) :: num_soilp
      integer, intent(in) :: filter_soilp(:)
      integer, intent(in) :: num_soilc
      integer, intent(in) :: filter_soilc(:)
      type(column_carbon_flux), intent(inout) :: col_cf_input
      !
      ! !LOCAL VARIABLES
      integer :: fp, p
      !------------------------------------------------------------
      associate( &
        rr_patch => vegcf%rr, &
        rr_col   => col_cf_input%rr &
        )
      do fp = 1,num_soilp
        p = filter_soilp(fp)
        ! root respiration (RR)
        rr_patch(p) = &
        vegcf%froot_mr(p) + &
        vegcf%cpool_froot_gr(p) + &
        vegcf%cpool_livecroot_gr(p) + &
        vegcf%cpool_deadcroot_gr(p) + &
        vegcf%transfer_froot_gr(p) + &
        vegcf%transfer_livecroot_gr(p) + &
        vegcf%transfer_deadcroot_gr(p) + &
        vegcf%cpool_froot_storage_gr(p) + &
        vegcf%cpool_livecroot_storage_gr(p) + &
        vegcf%cpool_deadcroot_storage_gr(p)
      enddo
        call p2c_1d_filter(bounds, num_soilc, filter_soilc, &
                rr_patch(bounds%begp:bounds%endp), &
                rr_col(bounds%begc:bounds%endc))
      end associate
    end subroutine vegcf_summary_rr_acc

    !-----------------------------------------------------------------------
    subroutine vegcf_summary_acc(vegcf, bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, isotope, col_cf_input)
      !
      ! !DESCRIPTION:
      ! patch-level carbon flux summary calculations
      !
      ! !USES:
      !$acc routine seq
      ! !ARGUMENTS:
      type(vegetation_carbon_flux)                 :: vegcf
      type(bounds_type)      , intent(in)    :: bounds
      integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
      integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
      integer                , intent(in)    :: num_soilc       ! number of soil columns in filter
      integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
      integer       , intent(in)           :: isotope !c13 = 0, c14 = 1, bulk = 2
      type(column_carbon_flux), intent(inout):: col_cf_input    ! receives p2c output
      !
      ! !LOCAL VARIABLES:
      integer  :: p,j,k,l       ! indices
      integer  :: fp            ! lake filter indices
      !-----------------------------------------------------------------------
      associate( &
        gpp_patch => vegcf%gpp , &
        gpp_col   => col_cf_input%gpp , &
        ar_patch  => vegcf%ar , &
        ar_col    => col_cf_input%ar , &
        npp_patch => vegcf%npp , &
        npp_col   => col_cf_input%npp , &
        vegfire_patch => vegcf%vegfire , &
        vegfire_col   => col_cf_input%vegfire , &
        wood_harvestc_patch => vegcf%wood_harvestc , &
        wood_harvestc_col   => col_cf_input%wood_harvestc , &
        fire_closs_patch => vegcf%fire_closs , &
        fire_closs_col   => col_cf_input%fire_closs_p2c , &
        litfall_patch => vegcf%litfall , &
        litfall_col   => col_cf_input%litfall , &
        hrv_xsmrpool_to_atm_patch => vegcf%hrv_xsmrpool_to_atm , &
        hrv_xsmrpool_to_atm_col   => col_cf_input%hrv_xsmrpool_to_atm  &
        )
      if (use_fates) return

      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! calculate pft-level summary carbon fluxes and states

         ! gross primary production (GPP)
          gpp_patch(p) = &
              vegcf%psnsun_to_cpool(p) + &
              vegcf%psnshade_to_cpool(p)

         ! maintenance respiration (MR)
         if ( isotope == c13 .or. isotope == c14) then
            vegcf%leaf_mr(p)      = vegcf%leaf_curmr(p)      + vegcf%leaf_xsmr(p)
            vegcf%froot_mr(p)     = vegcf%froot_curmr(p)     + vegcf%froot_xsmr(p)
            vegcf%livestem_mr(p)  = vegcf%livestem_curmr(p)  + vegcf%livestem_xsmr(p)
            vegcf%livecroot_mr(p) = vegcf%livecroot_curmr(p) + vegcf%livecroot_xsmr(p)
         endif

         vegcf%mr(p)  = &
              vegcf%leaf_mr(p)     + &
              vegcf%froot_mr(p)    + &
              vegcf%livestem_mr(p) + &
              vegcf%livecroot_mr(p)

         ! growth respiration (GR)
         ! current GR is respired vegcf time step for new growth displayed in vegcf timestep
         vegcf%current_gr(p) = &
              vegcf%cpool_leaf_gr(p)      + &
              vegcf%cpool_froot_gr(p)     + &
              vegcf%cpool_livestem_gr(p)  + &
              vegcf%cpool_deadstem_gr(p)  + &
              vegcf%cpool_livecroot_gr(p) + &
              vegcf%cpool_deadcroot_gr(p)

         ! transfer GR is respired vegcf time step for transfer growth displayed in vegcf timestep
         vegcf%transfer_gr(p) = &
              vegcf%transfer_leaf_gr(p)      + &
              vegcf%transfer_froot_gr(p)     + &
              vegcf%transfer_livestem_gr(p)  + &
              vegcf%transfer_deadstem_gr(p)  + &
              vegcf%transfer_livecroot_gr(p) + &
              vegcf%transfer_deadcroot_gr(p)

         ! storage GR is respired vegcf time step for growth sent to storage for later display
         vegcf%storage_gr(p) = &
              vegcf%cpool_leaf_storage_gr(p)      + &
              vegcf%cpool_froot_storage_gr(p)     + &
              vegcf%cpool_livestem_storage_gr(p)  + &
              vegcf%cpool_deadstem_storage_gr(p)  + &
              vegcf%cpool_livecroot_storage_gr(p) + &
              vegcf%cpool_deadcroot_storage_gr(p)

         if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
            vegcf%mr(p) = &
                 vegcf%mr(p) + &
                 vegcf%grain_mr(p)

            vegcf%current_gr(p) = &
                 vegcf%current_gr(p) + &
                 vegcf%cpool_grain_gr(p)

            vegcf%transfer_gr(p) = &
                 vegcf%transfer_gr(p) + &
                 vegcf%transfer_grain_gr(p)

            vegcf%storage_gr(p) = &
                 vegcf%storage_gr(p) + &
                 vegcf%cpool_grain_storage_gr(p)
         end if

         ! GR is the sum of current + transfer + storage GR
         vegcf%gr(p) = &
              vegcf%current_gr(p)  + &
              vegcf%transfer_gr(p) + &
              vegcf%storage_gr(p)

         ! autotrophic respiration (AR)
         if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
            ar_patch(p) = &
                 vegcf%mr(p) + &
                 vegcf%gr(p) + &
                 vegcf%xr(p) + &
                 vegcf%xsmrpool_to_atm(p) ! xsmr... is -ve (slevis)
            if (nu_com .ne. 'RD' ) then
               ar_patch(p) = ar_patch(p) + &
                    vegcf%xsmrpool_turnover(p)
            end if
         else
            ar_patch(p) = &
                 vegcf%mr(p) + &
                 vegcf%gr(p) + &
                 vegcf%xr(p)
            if (nu_com .ne. 'RD' ) then
               ar_patch(p) = ar_patch(p) + &
                    vegcf%xsmrpool_turnover(p)
            end if
         end if

         ! net primary production (NPP)
         npp_patch(p) = &
              gpp_patch(p) - &
              ar_patch(p)

         ! update the annual NPP accumulator, for use in allocation code
         if (isotope == bulk) then
            vegcf%tempsum_npp(p) = &
                 vegcf%tempsum_npp(p) + &
                 npp_patch(p)
         end if

         ! litterfall (LITFALL)

         litfall_patch(p) = &
              vegcf%leafc_to_litter(p)                     + &
              vegcf%frootc_to_litter(p)                    + &
              vegcf%m_leafc_to_litter(p)                   + &
              vegcf%m_leafc_storage_to_litter(p)           + &
              vegcf%m_leafc_xfer_to_litter(p)              + &
              vegcf%m_frootc_to_litter(p)                  + &
              vegcf%m_frootc_storage_to_litter(p)          + &
              vegcf%m_frootc_xfer_to_litter(p)             + &
              vegcf%m_livestemc_to_litter(p)               + &
              vegcf%m_livestemc_storage_to_litter(p)       + &
              vegcf%m_livestemc_xfer_to_litter(p)          + &
              vegcf%m_deadstemc_to_litter(p)               + &
              vegcf%m_deadstemc_storage_to_litter(p)       + &
              vegcf%m_deadstemc_xfer_to_litter(p)          + &
              vegcf%m_livecrootc_to_litter(p)              + &
              vegcf%m_livecrootc_storage_to_litter(p)      + &
              vegcf%m_livecrootc_xfer_to_litter(p)         + &
              vegcf%m_deadcrootc_to_litter(p)              + &
              vegcf%m_deadcrootc_storage_to_litter(p)      + &
              vegcf%m_deadcrootc_xfer_to_litter(p)         + &
              vegcf%m_gresp_storage_to_litter(p)           + &
              vegcf%m_gresp_xfer_to_litter(p)              + &
              vegcf%m_leafc_to_litter_fire(p)              + &
              vegcf%m_leafc_storage_to_litter_fire(p)      + &
              vegcf%m_leafc_xfer_to_litter_fire(p)         + &
              vegcf%m_livestemc_to_litter_fire(p)          + &
              vegcf%m_livestemc_storage_to_litter_fire(p)  + &
              vegcf%m_livestemc_xfer_to_litter_fire(p)     + &
              vegcf%m_deadstemc_to_litter_fire(p)          + &
              vegcf%m_deadstemc_storage_to_litter_fire(p)  + &
              vegcf%m_deadstemc_xfer_to_litter_fire(p)     + &
              vegcf%m_frootc_to_litter_fire(p)             + &
              vegcf%m_frootc_storage_to_litter_fire(p)     + &
              vegcf%m_frootc_xfer_to_litter_fire(p)        + &
              vegcf%m_livecrootc_to_litter_fire(p)         + &
              vegcf%m_livecrootc_storage_to_litter_fire(p) + &
              vegcf%m_livecrootc_xfer_to_litter_fire(p)    + &
              vegcf%m_deadcrootc_to_litter_fire(p)         + &
              vegcf%m_deadcrootc_storage_to_litter_fire(p) + &
              vegcf%m_deadcrootc_xfer_to_litter_fire(p)    + &
              vegcf%m_gresp_storage_to_litter_fire(p)      + &
              vegcf%m_gresp_xfer_to_litter_fire(p)         + &
              vegcf%hrv_leafc_to_litter(p)                 + &
              vegcf%hrv_leafc_storage_to_litter(p)         + &
              vegcf%hrv_leafc_xfer_to_litter(p)            + &
              vegcf%hrv_frootc_to_litter(p)                + &
              vegcf%hrv_frootc_storage_to_litter(p)        + &
              vegcf%hrv_frootc_xfer_to_litter(p)           + &
              vegcf%hrv_livestemc_to_litter(p)             + &
              vegcf%hrv_livestemc_storage_to_litter(p)     + &
              vegcf%hrv_livestemc_xfer_to_litter(p)        + &
              vegcf%hrv_deadstemc_storage_to_litter(p)     + &
              vegcf%hrv_deadstemc_xfer_to_litter(p)        + &
              vegcf%hrv_livecrootc_to_litter(p)            + &
              vegcf%hrv_livecrootc_storage_to_litter(p)    + &
              vegcf%hrv_livecrootc_xfer_to_litter(p)       + &
              vegcf%hrv_deadcrootc_to_litter(p)            + &
              vegcf%hrv_deadcrootc_storage_to_litter(p)    + &
              vegcf%hrv_deadcrootc_xfer_to_litter(p)       + &
              vegcf%hrv_gresp_storage_to_litter(p)         + &
              vegcf%hrv_gresp_xfer_to_litter(p)            + &
              vegcf%hrv_cpool_to_litter(p)

         ! patch-level fire losses (VEGFIRE)
         vegfire_patch(p) = 0._r8

         ! patch-level wood harvest
         wood_harvestc_patch(p) = &
              vegcf%hrv_deadstemc_to_prod10c(p) + &
              vegcf%hrv_deadstemc_to_prod100c(p)
         if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
            wood_harvestc_patch(p) = &
                 wood_harvestc_patch(p) + &
                 vegcf%hrv_cropc_to_prod1c(p)
         end if

         ! patch-level carbon losses to fire changed by F. Li and S. Levis
         fire_closs_patch(p) = &
              vegcf%m_leafc_to_fire(p)                + &
              vegcf%m_leafc_storage_to_fire(p)        + &
              vegcf%m_leafc_xfer_to_fire(p)           + &
              vegcf%m_frootc_to_fire(p)               + &
              vegcf%m_frootc_storage_to_fire(p)       + &
              vegcf%m_frootc_xfer_to_fire(p)          + &
              vegcf%m_livestemc_to_fire(p)            + &
              vegcf%m_livestemc_storage_to_fire(p)    + &
              vegcf%m_livestemc_xfer_to_fire(p)       + &
              vegcf%m_deadstemc_to_fire(p)            + &
              vegcf%m_deadstemc_storage_to_fire(p)    + &
              vegcf%m_deadstemc_xfer_to_fire(p)       + &
              vegcf%m_livecrootc_to_fire(p)           + &
              vegcf%m_livecrootc_storage_to_fire(p)   + &
              vegcf%m_livecrootc_xfer_to_fire(p)      + &
              vegcf%m_deadcrootc_to_fire(p)           + &
              vegcf%m_deadcrootc_storage_to_fire(p)   + &
              vegcf%m_deadcrootc_xfer_to_fire(p)      + &
              vegcf%m_gresp_storage_to_fire(p)        + &
              vegcf%m_gresp_xfer_to_fire(p)           + &
              vegcf%m_cpool_to_fire(p)

         if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
            litfall_patch(p) =                  &
                 litfall_patch(p)             + &
                 vegcf%livestemc_to_litter(p) + &
                 vegcf%grainc_to_food(p)
         end if

         ! new summary variables for CLAMP

         ! (FROOTC_ALLOC) - fine root C allocation
         vegcf%frootc_alloc(p) = &
              vegcf%frootc_xfer_to_frootc(p)    + &
              vegcf%cpool_to_frootc(p)

         ! (FROOTC_LOSS) - fine root C loss changed by F. Li and S. Levis
         vegcf%frootc_loss(p) = &
              vegcf%m_frootc_to_litter(p)       + &
              vegcf%m_frootc_to_fire(p)         + &
              vegcf%m_frootc_to_litter_fire(p)  + &
              vegcf%hrv_frootc_to_litter(p)     + &
              vegcf%frootc_to_litter(p)

         ! (LEAFC_ALLOC) - leaf C allocation
         vegcf%leafc_alloc(p) = &
              vegcf%leafc_xfer_to_leafc(p)    + &
              vegcf%cpool_to_leafc(p)

         ! (LEAFC_LOSS) - leaf C loss changed by F. Li and S. Levis
         vegcf%leafc_loss(p) = &
              vegcf%m_leafc_to_litter(p)      + &
              vegcf%m_leafc_to_fire(p)        + &
              vegcf%m_leafc_to_litter_fire(p) + &
              vegcf%hrv_leafc_to_litter(p)    + &
              vegcf%leafc_to_litter(p)

         if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
            vegcf%leafc_loss(p) = &
                 vegcf%leafc_loss(p) + &
                 vegcf%hrv_leafc_to_prod1c(p)
         end if


         ! (WOODC_ALLOC) - wood C allocation
         vegcf%woodc_alloc(p) = &
              vegcf%livestemc_xfer_to_livestemc(p)   + &
              vegcf%deadstemc_xfer_to_deadstemc(p)   + &
              vegcf%livecrootc_xfer_to_livecrootc(p) + &
              vegcf%deadcrootc_xfer_to_deadcrootc(p) + &
              vegcf%cpool_to_livestemc(p)            + &
              vegcf%cpool_to_deadstemc(p)            + &
              vegcf%cpool_to_livecrootc(p)           + &
              vegcf%cpool_to_deadcrootc(p)

         ! (WOODC_LOSS) - wood C loss
         vegcf%woodc_loss(p) = &
              vegcf%m_livestemc_to_litter(p)            + &
              vegcf%m_deadstemc_to_litter(p)            + &
              vegcf%m_livecrootc_to_litter(p)           + &
              vegcf%m_deadcrootc_to_litter(p)           + &
              vegcf%m_livestemc_to_fire(p)              + &
              vegcf%m_deadstemc_to_fire(p)              + &
              vegcf%m_livecrootc_to_fire(p)             + &
              vegcf%m_deadcrootc_to_fire(p)             + &
              vegcf%hrv_livestemc_to_litter(p)          + &
              vegcf%hrv_livestemc_storage_to_litter(p)  + &
              vegcf%hrv_livestemc_xfer_to_litter(p)     + &
              vegcf%hrv_deadstemc_to_prod10c(p)         + &
              vegcf%hrv_deadstemc_to_prod100c(p)        + &
              vegcf%hrv_deadstemc_storage_to_litter(p)  + &
              vegcf%hrv_deadstemc_xfer_to_litter(p)     + &
              vegcf%hrv_livecrootc_to_litter(p)         + &
              vegcf%hrv_livecrootc_storage_to_litter(p) + &
              vegcf%hrv_livecrootc_xfer_to_litter(p)    + &
              vegcf%hrv_deadcrootc_to_litter(p)         + &
              vegcf%hrv_deadcrootc_storage_to_litter(p) + &
              vegcf%hrv_deadcrootc_xfer_to_litter(p)
         ! putting the harvested crop stem and grain in the wood loss bdrewniak
         if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
            vegcf%woodc_loss(p) = &
                 vegcf%woodc_loss(p) + &
                 vegcf%hrv_grainc_to_prod1c(p) + &
                 vegcf%hrv_livestemc_to_prod1c(p)
         end if
      end do  ! end of patches loop

      ! use p2c routine to get selected column-average patch-level fluxes and states
      call p2c(bounds, num_soilc, filter_soilc, &
              gpp_patch(bounds%begp:bounds%endp), &
              gpp_col(bounds%begc:bounds%endc))

      call p2c(bounds, num_soilc, filter_soilc, &
              ar_patch(bounds%begp:bounds%endp), &
              ar_col  (bounds%begc:bounds%endc))

      call p2c(bounds, num_soilc, filter_soilc, &
              npp_patch(bounds%begp:bounds%endp), &
              npp_col(bounds%begc:bounds%endc))

      call p2c(bounds, num_soilc, filter_soilc, &
              vegfire_patch(bounds%begp:bounds%endp), &
              vegfire_col(bounds%begc:bounds%endc))

      call p2c(bounds, num_soilc, filter_soilc, &
           wood_harvestc_patch(bounds%begp:bounds%endp), &
           wood_harvestc_col(bounds%begc:bounds%endc))

      call p2c(bounds, num_soilc, filter_soilc, &
           fire_closs_patch(bounds%begp:bounds%endp), &
           fire_closs_col(bounds%begc:bounds%endc))

      call p2c(bounds, num_soilc, filter_soilc, &
           litfall_patch(bounds%begp:bounds%endp), &
           litfall_col(bounds%begc:bounds%endc))

      call p2c(bounds, num_soilc, filter_soilc, &
           hrv_xsmrpool_to_atm_patch(bounds%begp:bounds%endp), &
           hrv_xsmrpool_to_atm_col(bounds%begc:bounds%endc))

      end associate

    end subroutine vegcf_summary_acc

    !------------------------------------------------------------
    subroutine vegcf_summary_for_ch4_acc( vegcf, bounds, num_soilp, filter_soilp)
      !$acc routine seq
      ! !DESCRIPTION:
      ! summarize vegetation-level fluxes for methane calculation
      !
      ! !USES:
      use VegetationDataType, only : vegetation_carbon_flux
      !
      ! !ARGUMENTS:
      type(vegetation_carbon_flux), intent(inout) :: vegcf
      type(bounds_type), intent(in) :: bounds
      integer, intent(in) :: num_soilp
      integer, intent(in) :: filter_soilp(:)
      !
      ! !LOCAL VARIABLES
      integer :: fp, p
      !------------------------------------------------------------

      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! aboveground NPP: leaf, live stem, dead stem (AGNPP)
         ! vegcf is supposed to correspond as closely as possible to
         ! field measurements of AGNPP, so it ignores the storage pools
         ! and only treats the fluxes into displayed pools.

         vegcf%agnpp(p) = &
              vegcf%cpool_to_leafc(p)                  + &
              vegcf%leafc_xfer_to_leafc(p)             + &
              vegcf%cpool_to_livestemc(p)              + &
              vegcf%livestemc_xfer_to_livestemc(p)     + &
              vegcf%cpool_to_deadstemc(p)              + &
              vegcf%deadstemc_xfer_to_deadstemc(p)

         if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
            vegcf%agnpp(p) =                    &
                 vegcf%agnpp(p)               + &
                 vegcf%cpool_to_grainc(p)     + &
                 vegcf%grainc_xfer_to_grainc(p)
         endif

         ! belowground NPP: fine root, live coarse root, dead coarse root (BGNPP)
         ! vegcf is supposed to correspond as closely as possible to
         ! field measurements of BGNPP, so it ignores the storage pools
         ! and only treats the fluxes into displayed pools.
         vegcf%bgnpp(p) = &
              vegcf%cpool_to_frootc(p)                   + &
              vegcf%frootc_xfer_to_frootc(p)             + &
              vegcf%cpool_to_livecrootc(p)               + &
              vegcf%livecrootc_xfer_to_livecrootc(p)     + &
              vegcf%cpool_to_deadcrootc(p)               + &
              vegcf%deadcrootc_xfer_to_deadcrootc(p)

         vegcf%agwdnpp(p) = &
              vegcf%cpool_to_livestemc(p)              + &
              vegcf%livestemc_xfer_to_livestemc(p)     + &
              vegcf%cpool_to_deadstemc(p)              + &
              vegcf%deadstemc_xfer_to_deadstemc(p)

      end do

    end subroutine vegcf_summary_for_ch4_acc

    !-----------------------------------------------------------------------
    subroutine vegps_setvalues_acc( vegps, num_patch, filter_patch, value_patch)
      !$acc routine seq
      ! !DESCRIPTION:
      ! Set phosphorus state variables, column-level
      !
      use VegetationDataType, only : vegetation_phosphorus_state
      ! !ARGUMENTS:
      type (vegetation_phosphorus_state), intent(inout) :: vegps
      integer , intent(in)                :: num_patch
      integer , intent(in)                :: filter_patch(:)
      real(r8), intent(in)                :: value_patch
      !
      ! !LOCAL VARIABLES:
      integer :: fi,i     ! loop index
      integer :: j,k      ! indices
      !------------------------------------------------------------------------
      do fi = 1,num_patch
         i = filter_patch(fi)

         vegps%leafp(i)              = value_patch
         vegps%leafp_storage(i)      = value_patch
         vegps%leafp_xfer(i)         = value_patch
         vegps%frootp(i)             = value_patch
         vegps%frootp_storage(i)     = value_patch
         vegps%frootp_xfer(i)        = value_patch
         vegps%livestemp(i)          = value_patch
         vegps%livestemp_storage(i)  = value_patch
         vegps%livestemp_xfer(i)     = value_patch
         vegps%deadstemp(i)          = value_patch
         vegps%deadstemp_storage(i)  = value_patch
         vegps%deadstemp_xfer(i)     = value_patch
         vegps%livecrootp(i)         = value_patch
         vegps%livecrootp_storage(i) = value_patch
         vegps%livecrootp_xfer(i)    = value_patch
         vegps%deadcrootp(i)         = value_patch
         vegps%deadcrootp_storage(i) = value_patch
         vegps%deadcrootp_xfer(i)    = value_patch
         vegps%retransp(i)           = value_patch
         vegps%ppool(i)              = value_patch
         vegps%ptrunc(i)             = value_patch
         vegps%dispvegp(i)           = value_patch
         vegps%storvegp(i)           = value_patch
         vegps%totvegp(i)            = value_patch
         vegps%totpftp(i)            = value_patch
      end do

      if ( crop_prog )then
         do fi = 1,num_patch
            i = filter_patch(fi)
            vegps%grainp(i)            = value_patch
            vegps%grainp_storage(i)    = value_patch
            vegps%grainp_xfer(i)       = value_patch
            vegps%cropseedp_deficit(i) = value_patch
         end do
      end if

    end subroutine vegps_setvalues_acc


    !-----------------------------------------------------------------------
    subroutine vegcs_summary_acc(vegcs, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_cs)
      !
      ! !DESCRIPTION:
      ! Vegetation-level carbon state summary calculations
      !$acc routine seq
      ! !ARGUMENTS:
      type(vegetation_carbon_state)             :: vegcs
      type(bounds_type)         , intent(in)    :: bounds
      integer                   , intent(in)    :: num_soilc       ! number of soil columns in filter
      integer                   , intent(in)    :: filter_soilc(:) ! filter for soil columns
      integer                   , intent(in)    :: num_soilp       ! number of soil patches in filter
      integer                   , intent(in)    :: filter_soilp(:) ! filter for soil patches
      type (column_carbon_state), intent(inout) :: col_cs          ! column-level state for p2c

      !
      ! !LOCAL VARIABLES:
      integer  :: c,p             ! indices
      integer  :: fp              ! filter indices
      real(r8) :: maxdepth        ! depth to integrate soil variables
      !-----------------------------------------------------------------------
      associate(&
        totpftc_patch  => vegcs%totpftc         , &
        totpftc_col    => col_cs%totpftc, &
        totvegc_patch  => vegcs%totvegc , &
        totvegc_col    => col_cs%totvegc, &
        totvegc_abg_patch  => vegcs%totvegc_abg , &
        totvegc_abg_col    => col_cs%totvegc_abg, &
        cropseedc_deficit_patch  => vegcs%cropseedc_deficit ,&
        cropseedc_deficit_col    => col_cs%cropseedc_deficit &
        )
      if (use_fates) return

      do fp = 1,num_soilp
         p = filter_soilp(fp)
         ! displayed vegetation carbon, excluding storage and cpool (DISPVEGC)
         vegcs%dispvegc(p) =        &
              vegcs%leafc(p)      + &
              vegcs%frootc(p)     + &
              vegcs%livestemc(p)  + &
              vegcs%deadstemc(p)  + &
              vegcs%livecrootc(p) + &
              vegcs%deadcrootc(p)

         ! stored vegetation carbon, excluding cpool (STORVEGC)
         vegcs%storvegc(p) =                &
              vegcs%cpool(p)              + &
              vegcs%leafc_storage(p)      + &
              vegcs%frootc_storage(p)     + &
              vegcs%livestemc_storage(p)  + &
              vegcs%deadstemc_storage(p)  + &
              vegcs%livecrootc_storage(p) + &
              vegcs%deadcrootc_storage(p) + &
              vegcs%leafc_xfer(p)         + &
              vegcs%frootc_xfer(p)        + &
              vegcs%livestemc_xfer(p)     + &
              vegcs%deadstemc_xfer(p)     + &
              vegcs%livecrootc_xfer(p)    + &
              vegcs%deadcrootc_xfer(p)    + &
              vegcs%gresp_storage(p)      + &
              vegcs%gresp_xfer(p)

         if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
            vegcs%storvegc(p) =            &
                 vegcs%storvegc(p)       + &
                 vegcs%grainc_storage(p) + &
                 vegcs%grainc_xfer(p)

            vegcs%dispvegc(p) =            &
                 vegcs%dispvegc(p)       + &
                 vegcs%grainc(p)
         end if

         ! total vegetation carbon, excluding cpool (TOTVEGC)
         totvegc_patch(p) = &
              vegcs%dispvegc(p) + &
              vegcs%storvegc(p)

         ! total pft-level carbon, including xsmrpool, ctrunc
         totpftc_patch(p) = &
              vegcs%totvegc(p) + &
              vegcs%xsmrpool(p) + &
              vegcs%ctrunc(p)

         ! (WOODC) - wood C
         vegcs%woodc(p) = &
              vegcs%deadstemc(p)    + &
              vegcs%livestemc(p)    + &
              vegcs%deadcrootc(p)   + &
              vegcs%livecrootc(p)

         totvegc_abg_patch(p) = &
                 vegcs%leafc(p)              + &
                 vegcs%leafc_storage(p)      + &
                 vegcs%leafc_xfer(p)         + &
                 vegcs%livestemc(p)          + &
                 vegcs%livestemc_storage(p)  + &
                 vegcs%livestemc_xfer(p)     + &
                 vegcs%deadstemc(p)          + &
                 vegcs%deadstemc_storage(p)  + &
                 vegcs%deadstemc_xfer(p)
      end do ! filtered veg list

      ! a few vegetation-to-column summaries
      call p2c(bounds, num_soilc, filter_soilc, &
           totpftc_patch(bounds%begp:bounds%endp) , &
           totpftc_col(bounds%begc:bounds%endc))

      call p2c(bounds, num_soilc, filter_soilc, &
           totvegc_patch(bounds%begp:bounds%endp) , &
           totvegc_col(bounds%begc:bounds%endc))

      call p2c(bounds, num_soilc, filter_soilc, &
           totvegc_abg_patch(bounds%begp:bounds%endp), &
           totvegc_abg_col(bounds%begc:bounds%endc))

      call p2c(bounds, num_soilc, filter_soilc, &
           cropseedc_deficit_patch(bounds%begp:bounds%endp), &
           cropseedc_deficit_col(bounds%begc:bounds%endc))
      end associate
    end subroutine vegcs_summary_acc

    !-----------------------------------------------------------------------
    subroutine colcs_summary_acc(colcs, bounds, num_soilc, filter_soilc)
      !
      ! !DESCRIPTION:
      ! Column-level carbon state summary calculations
      !$acc routine seq
      ! !ARGUMENTS:
      type(column_carbon_state) :: colcs
      type(bounds_type)      , intent(in)    :: bounds
      integer                , intent(in)    :: num_soilc       ! number of soil columns in filter
      integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
      !
      ! !LOCAL VARIABLES:
      real(r8) :: nfixlags, dtime ! temp variables for making lagged npp
      integer  :: c,p,j,k,l       ! indices
      integer  :: fp,fc           ! lake filter indices
      real(r8) :: maxdepth        ! depth to integrate soil variables
      integer  :: nlev
      !-----------------------------------------------------------------------

      if (use_fates) return

      nlev = nlevdecomp
      if (use_pflotran .and. pf_cmode) nlev = nlevdecomp_full

      ! vertically integrate each of the decomposing C pools
      do l = 1, ndecomp_pools
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            colcs%decomp_cpools(c,l) = 0._r8
         end do
      end do

      do l = 1, ndecomp_pools
         do j = 1, nlev
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               colcs%decomp_cpools(c,l) = &
                    colcs%decomp_cpools(c,l) + &
                    colcs%decomp_cpools_vr(c,j,l) * dzsoi_decomp(j)
            end do
         end do
      end do

      if ( nlevdecomp > 1) then
         ! vertically integrate each of the decomposing C pools to 1 meter
         maxdepth = 1._r8
         do l = 1, ndecomp_pools
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               colcs%decomp_cpools_1m(c,l) = 0._r8
            end do
         end do
         do l = 1, ndecomp_pools
            do j = 1, nlevdecomp
               if ( zisoi(j) <= maxdepth ) then
                  do fc = 1,num_soilc
                     c = filter_soilc(fc)
                     colcs%decomp_cpools_1m(c,l) = &
                          colcs%decomp_cpools_1m(c,l) + &
                          colcs%decomp_cpools_vr(c,j,l) * dzsoi_decomp(j)
                  end do
               elseif ( zisoi(j-1) < maxdepth ) then
                  do fc = 1,num_soilc
                     c = filter_soilc(fc)
                     colcs%decomp_cpools_1m(c,l) = &
                          colcs%decomp_cpools_1m(c,l) + &
                          colcs%decomp_cpools_vr(c,j,l) * (maxdepth - zisoi(j-1))
                  end do
               endif
            end do
         end do

         ! total litter carbon in the top meter (TOTLITC_1m)
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            colcs%totlitc_1m(c)         = 0._r8
         end do
         do l = 1, ndecomp_pools
            if ( decomp_cascade_con%is_litter(l) ) then
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  colcs%totlitc_1m(c) = &
                       colcs%totlitc_1m(c) + &
                       colcs%decomp_cpools_1m(c,l)
               end do
            endif
         end do

         ! total soil organic matter carbon in the top meter (TOTSOMC_1m)
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            colcs%totsomc_1m(c) = 0._r8
         end do
         do l = 1, ndecomp_pools
            if ( decomp_cascade_con%is_soil(l) ) then
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  colcs%totsomc_1m(c) = &
                       colcs%totsomc_1m(c) + &
                       colcs%decomp_cpools_1m(c,l)
               end do
            end if
         end do
      end if ! nlevdecomp>1

      ! total litter carbon (TOTLITC)
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         colcs%totlitc(c) = 0._r8
      end do
      do l = 1, ndecomp_pools
         if ( decomp_cascade_con%is_litter(l) ) then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               colcs%totlitc(c) = &
                    colcs%totlitc(c) + &
                    colcs%decomp_cpools(c,l)
            end do
         endif
      end do

      ! total soil organic matter carbon (TOTSOMC)
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         colcs%totsomc(c) = 0._r8
      end do
      do l = 1, ndecomp_pools
         if ( decomp_cascade_con%is_soil(l) ) then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               colcs%totsomc(c) = &
                    colcs%totsomc(c) + &
                    colcs%decomp_cpools(c,l)
            end do
         end if
      end do

      ! coarse woody debris carbon
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         colcs%cwdc(c) = 0._r8
      end do
      do l = 1, ndecomp_pools
         if ( decomp_cascade_con%is_cwd(l) ) then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               colcs%cwdc(c) = &
                    colcs%cwdc(c) + &
                    colcs%decomp_cpools(c,l)
            end do
         end if
      end do

      ! truncation carbon
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         colcs%ctrunc(c) = 0._r8
      end do
      do j = 1, nlev
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            colcs%ctrunc(c) = &
                 colcs%ctrunc(c) + &
                 colcs%ctrunc_vr(c,j) * dzsoi_decomp(j)
         end do
      end do

      do fc = 1,num_soilc
         c = filter_soilc(fc)

         ! total product carbon
         colcs%totprodc(c) = &
              colcs%prod10c(c)  + &
              colcs%prod100c(c) + &
              colcs%prod1c(c)

         ! total ecosystem carbon, including veg but excluding cpool (TOTECOSYSC)
         colcs%totecosysc(c) = &
              colcs%cwdc(c)     + &
              colcs%totlitc(c)  + &
              colcs%totsomc(c)  + &
              colcs%totprodc(c) + &
              colcs%totvegc(c)

         ! total column carbon, including veg and cpool (TOTCOLC)
         ! adding col_ctrunc, seedc
         colcs%totcolc(c) = &
              colcs%totpftc(c)  + &
              colcs%cwdc(c)     + &
              colcs%totlitc(c)  + &
              colcs%totsomc(c)  + &
              colcs%prod1c(c)   + &
              colcs%ctrunc(c)   + &
              colcs%cropseedc_deficit(c)
         colcs%totabgc(c) = &
              colcs%totpftc(c)  + &
              colcs%totprodc(c) + &
              colcs%seedc(c)    + &
              colcs%ctrunc(c)
      end do
    end subroutine colcs_summary_acc


    !-----------------------------------------------------------------------
    subroutine colns_summary_acc(colns, bounds, num_soilc, filter_soilc)
      !$acc routine seq
      ! !ARGUMENTS:
      type(column_nitrogen_state)  :: colns
      type(bounds_type) , intent(in) :: bounds
      integer           , intent(in) :: num_soilc       ! number of soil columns in filter
      integer           , intent(in) :: filter_soilc(:) ! filter for soil columns
      !
      ! !LOCAL VARIABLES:
      integer  :: c,p,j,k,l   ! indices
      integer  :: fp,fc       ! lake filter indices
      real(r8) :: maxdepth    ! depth to integrate soil variables
      integer  :: nlev
      !-----------------------------------------------------------------------

      ! vertically integrate NO3 NH4 N2O pools
      nlev = nlevdecomp
      if (use_pflotran .and. pf_cmode) nlev = nlevdecomp_full

      if (use_nitrif_denitrif .or. (use_pflotran .and. pf_cmode)) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            colns%smin_no3(c) = 0._r8
            colns%smin_nh4(c) = 0._r8
            if(use_pflotran .and. pf_cmode) then
               colns%smin_nh4sorb(c) = 0._r8
            end if
         end do
         do j = 1, nlev
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               colns%smin_no3(c) = &
                    colns%smin_no3(c) + &
                    colns%smin_no3_vr(c,j) * dzsoi_decomp(j)

               colns%smin_nh4(c) = &
                    colns%smin_nh4(c) + &
                    colns%smin_nh4_vr(c,j) * dzsoi_decomp(j)
               if(use_pflotran .and. pf_cmode) then
                  colns%smin_nh4sorb(c) = &
                    colns%smin_nh4sorb(c) + &
                    colns%smin_nh4sorb_vr(c,j) * dzsoi_decomp(j)
               end if
             end do
          end do

       end if

      ! vertically integrate each of the decomposing N pools
      do l = 1, ndecomp_pools
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            colns%decomp_npools(c,l) = 0._r8
         end do
         do j = 1, nlev
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               colns%decomp_npools(c,l) = &
                    colns%decomp_npools(c,l) + &
                    colns%decomp_npools_vr(c,j,l) * dzsoi_decomp(j)
            end do
         end do
      end do

      ! for vertically-resolved soil biogeochemistry, calculate some diagnostics of carbon pools to a given depth
      if ( nlevdecomp > 1) then

         do l = 1, ndecomp_pools
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               colns%decomp_npools_1m(c,l) = 0._r8
            end do
         end do

         ! vertically integrate each of the decomposing n pools to 1 meter
         maxdepth = 1._r8
         do l = 1, ndecomp_pools
            do j = 1, nlevdecomp
               if ( zisoi(j) <= maxdepth ) then
                  do fc = 1,num_soilc
                     c = filter_soilc(fc)
                     colns%decomp_npools_1m(c,l) = &
                          colns%decomp_npools_1m(c,l) + &
                          colns%decomp_npools_vr(c,j,l) * dzsoi_decomp(j)
                  end do
               elseif ( zisoi(j-1) < maxdepth ) then
                  do fc = 1,num_soilc
                     c = filter_soilc(fc)
                     colns%decomp_npools_1m(c,l) = &
                          colns%decomp_npools_1m(c,l) + &
                          colns%decomp_npools_vr(c,j,l) * (maxdepth - zisoi(j-1))
                  end do
               endif
            end do
         end do

         ! total litter nitrogen to 1 meter (TOTLITN_1m)
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            colns%totlitn_1m(c) = 0._r8
         end do
         do l = 1, ndecomp_pools
            if ( decomp_cascade_con%is_litter(l) ) then
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  colns%totlitn_1m(c) = &
                       colns%totlitn_1m(c) + &
                       colns%decomp_npools_1m(c,l)
               end do
            end if
         end do

         ! total soil organic matter nitrogen to 1 meter (TOTSOMN_1m)
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            colns%totsomn_1m(c) = 0._r8
         end do
         do l = 1, ndecomp_pools
            if ( decomp_cascade_con%is_soil(l) ) then
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  colns%totsomn_1m(c) = &
                       colns%totsomn_1m(c) + &
                       colns%decomp_npools_1m(c,l)
               end do
            end if
         end do

      endif

      ! total litter nitrogen (TOTLITN)
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         colns%totlitn(c)    = 0._r8
      end do
      do l = 1, ndecomp_pools
         if ( decomp_cascade_con%is_litter(l) ) then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               colns%totlitn(c) = &
                    colns%totlitn(c) + &
                    colns%decomp_npools(c,l)
            end do
         end if
      end do

      ! total soil organic matter nitrogen (TOTSOMN)
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         colns%totsomn(c)    = 0._r8
      end do
      do l = 1, ndecomp_pools
         if ( decomp_cascade_con%is_soil(l) ) then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               colns%totsomn(c) = &
                    colns%totsomn(c) + &
                    colns%decomp_npools(c,l)
            end do
         end if
      end do

      ! total cwdn
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         colns%cwdn(c) = 0._r8
      end do
      do l = 1, ndecomp_pools
         if ( decomp_cascade_con%is_cwd(l) ) then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               colns%cwdn(c) = &
                    colns%cwdn(c) + &
                    colns%decomp_npools(c,l)
            end do
         end if
      end do

      ! total sminn
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         colns%sminn(c)      = 0._r8
      end do
      do j = 1, nlev
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            colns%sminn(c) = &
                 colns%sminn(c) + &
                 colns%sminn_vr(c,j) * dzsoi_decomp(j)
         end do
      end do

      ! total col_ntrunc
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         colns%ntrunc(c) = 0._r8
      end do
      do j = 1, nlev
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            colns%ntrunc(c) = &
                 colns%ntrunc(c) + &
                 colns%ntrunc_vr(c,j) * dzsoi_decomp(j)
         end do
      end do

      do fc = 1,num_soilc
         c = filter_soilc(fc)

         ! total wood product nitrogen
         colns%totprodn(c) = &
              colns%prod1n(c) + &
              colns%prod10n(c) + &
              colns%prod100n(c)

         ! total ecosystem nitrogen, including veg (TOTECOSYSN)
         colns%totecosysn(c) = &
              colns%cwdn(c) + &
              colns%totlitn(c) + &
              colns%totsomn(c) + &
              colns%sminn(c) + &
              colns%totprodn(c) + &
              colns%totvegn(c)

         ! total column nitrogen, including pft (TOTCOLN)
         colns%totcoln(c) = &
              colns%totpftn(c) + &
              colns%cwdn(c) + &
              colns%totlitn(c) + &
              colns%totsomn(c) + &
              colns%sminn(c) + &
              colns%prod1n(c) + &
              colns%ntrunc(c)+ &
              colns%plant_n_buffer(c) + &
              colns%cropseedn_deficit(c)

         colns%totabgn (c) =  &
              colns%totpftn(c) + &
              colns%totprodn(c) + &
              colns%seedn(c) + &
              colns%ntrunc(c)+ &
              colns%plant_n_buffer(c)

         colns%totblgn(c) = &
              colns%cwdn(c) + &
              colns%totlitn(c) + &
              colns%totsomn(c) + &
              colns%sminn(c)
      end do

    end subroutine colns_summary_acc

    !-----------------------------------------------------------------------
    subroutine colpf_summary_acc(colpf, bounds, num_soilc, filter_soilc, dtime)
      !$acc routine seq
      ! !ARGUMENTS:
      type(column_phosphorus_flux) :: colpf
      type(bounds_type) , intent(in) :: bounds
      integer           , intent(in) :: num_soilc       ! number of soil columns in filter
      integer           , intent(in) :: filter_soilc(:) ! filter for soil columns
      real(r8)   , intent(in) :: dtime
      !
      ! !LOCAL VARIABLES:
      integer  :: c,j,k,l   ! indices
      integer  :: fc       ! lake filter indices
      real(r8) :: maxdepth    ! depth to integrate soil variables
      !-----------------------------------------------------------------------
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         colpf%supplement_to_sminp(c) = 0._r8
         colpf%som_p_leached(c)       = 0._r8
      end do

      ! pflotran
      if (.not.(use_pflotran .and. pf_cmode)) then
      ! vertically integrate decomposing P cascade fluxes and soil mineral P fluxes associated with decomposition cascade
      do k = 1, ndecomp_cascade_transitions
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               colpf%decomp_cascade_ptransfer(c,k) = &
                    colpf%decomp_cascade_ptransfer(c,k) + &
                    colpf%decomp_cascade_ptransfer_vr(c,j,k) * dzsoi_decomp(j)

               colpf%decomp_cascade_sminp_flux(c,k) = &
                    colpf%decomp_cascade_sminp_flux(c,k) + &
                    colpf%decomp_cascade_sminp_flux_vr(c,j,k) * dzsoi_decomp(j)
            end do
         end do
      end do
      end if !if (.not.(use_pflotran .and. pf_cmode))

      ! vertically integrate inorganic P flux
      do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            colpf%primp_to_labilep(c) = &
                 colpf%primp_to_labilep(c) + &
                 colpf%primp_to_labilep_vr(c,j) * dzsoi_decomp(j)
         end do
      end do

      do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            colpf%labilep_to_secondp(c) = &
                 colpf%labilep_to_secondp(c) + &
                 colpf%labilep_to_secondp_vr(c,j) * dzsoi_decomp(j)
         end do
      end do

      do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            colpf%secondp_to_labilep(c) = &
                 colpf%secondp_to_labilep(c) + &
                 colpf%secondp_to_labilep_vr(c,j) * dzsoi_decomp(j)
         end do
      end do

      do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            colpf%secondp_to_occlp(c) = &
                 colpf%secondp_to_occlp(c) + &
                 colpf%secondp_to_occlp_vr(c,j) * dzsoi_decomp(j)
         end do
      end do

      ! vertically integrate leaching flux
      do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            colpf%sminp_leached(c) = &
                 colpf%sminp_leached(c) + &
                 colpf%sminp_leached_vr(c,j) * dzsoi_decomp(j)
         end do
      end do

      ! vertically integrate column-level fire P losses
      do k = 1, ndecomp_pools
         do j = 1, nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               colpf%m_decomp_ppools_to_fire(c,k) = &
                    colpf%m_decomp_ppools_to_fire(c,k) + &
                    colpf%m_decomp_ppools_to_fire_vr(c,j,k) * dzsoi_decomp(j)
            end do
         end do
      end do

      ! total column-level fire P losses
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         colpf%fire_ploss(c) = colpf%fire_ploss_p2c(c)
      end do
      do k = 1, ndecomp_pools
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            colpf%fire_ploss(c) = &
                 colpf%fire_ploss(c) + &
                 colpf%m_decomp_ppools_to_fire(c,k)
         end do
      end do

      ! supplementary P supplement_to_sminp
      do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            colpf%supplement_to_sminp(c) = &
                 colpf%supplement_to_sminp(c) + &
                 colpf%supplement_to_sminp_vr(c,j) * dzsoi_decomp(j)
         end do
      end do

      do fc = 1,num_soilc
         c = filter_soilc(fc)

         ! column-level P losses due to landcover change
         colpf%dwt_ploss(c) = &
              colpf%dwt_conv_pflux(c)

         ! total wood product P loss
         colpf%product_ploss(c) = &
              colpf%prod10p_loss(c) + &
              colpf%prod100p_loss(c) + &
              colpf%prod1p_loss(c)
      end do

      ! add up all vertical transport tendency terms and calculate total som leaching loss as the sum of these
      do l = 1, ndecomp_pools
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            colpf%decomp_ppools_leached(c,l) = 0._r8
         end do

         do j = 1, nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               colpf%decomp_ppools_leached(c,l) = &
                    colpf%decomp_ppools_leached(c,l) + &
                    colpf%decomp_ppools_transport_tendency(c,j,l) * dzsoi_decomp(j)
            end do
         end do

         do fc = 1,num_soilc
            c = filter_soilc(fc)
            colpf%som_p_leached(c) = &
                 colpf%som_p_leached(c) + &
                 colpf%decomp_ppools_leached(c,l)
         end do
      end do

      do fc = 1,num_soilc
         c = filter_soilc(fc)
         colpf%biochem_pmin(c) = 0.0_r8
      end do

      do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            colpf%biochem_pmin(c) = colpf%biochem_pmin(c) + &
                 colpf%biochem_pmin_vr(c,j)* dzsoi_decomp(j)
         end do
      end do

      do fc = 1,num_soilc
         c = filter_soilc(fc)
         colpf%adsorb_to_labilep(c) = 0._r8
         colpf%desorb_to_solutionp(c) = 0._r8
      end do

      do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            colpf%adsorb_to_labilep(c) = colpf%adsorb_to_labilep(c) + &
                 colpf%adsorb_to_labilep_vr(c,j)* dzsoi_decomp(j)
            colpf%desorb_to_solutionp(c) = colpf%desorb_to_solutionp(c) + &
                 colpf%desorb_to_solutionp_vr(c,j)* dzsoi_decomp(j)
         end do
      end do

      do fc = 1,num_soilc
         c = filter_soilc(fc)
         colpf%actual_immob_p(c) = 0._r8
         colpf%smin_p_to_plant(c) = 0._r8
         colpf%plant_to_litter_pflux(c) = 0._r8
         colpf%plant_to_cwd_pflux(c) = 0._r8
         do j = 1, nlevdecomp
            colpf%actual_immob_p(c)= colpf%actual_immob_p(c) + &
                 colpf%actual_immob_p_vr(c,j) * dzsoi_decomp(j)
            colpf%smin_p_to_plant(c)= colpf%smin_p_to_plant(c) + &
                 colpf%sminp_to_plant_vr(c,j) * dzsoi_decomp(j)
            colpf%plant_to_litter_pflux(c) = &
                 colpf%plant_to_litter_pflux(c)  + &
                 colpf%phenology_p_to_litr_met_p(c,j)* dzsoi_decomp(j) + &
                 colpf%phenology_p_to_litr_cel_p(c,j)* dzsoi_decomp(j) + &
                 colpf%phenology_p_to_litr_lig_p(c,j)* dzsoi_decomp(j) + &
                 colpf%gap_mortality_p_to_litr_met_p(c,j)* dzsoi_decomp(j) + &
                 colpf%gap_mortality_p_to_litr_cel_p(c,j)* dzsoi_decomp(j) + &
                 colpf%gap_mortality_p_to_litr_lig_p(c,j)* dzsoi_decomp(j) + &
                 colpf%m_p_to_litr_met_fire(c,j)* dzsoi_decomp(j) + &
                 colpf%m_p_to_litr_cel_fire(c,j)* dzsoi_decomp(j) + &
                 colpf%m_p_to_litr_lig_fire(c,j)* dzsoi_decomp(j)
            colpf%plant_to_cwd_pflux(c) = &
                 colpf%plant_to_cwd_pflux(c) + &
                 colpf%gap_mortality_p_to_cwdp(c,j)* dzsoi_decomp(j) + &
                 colpf%fire_mortality_p_to_cwdp(c,j)* dzsoi_decomp(j)
         end do
      end do

      ! bgc interface & pflotran:
      if (use_clm_interface) then
          call colpf_SummaryInt_acc(colpf,bounds, num_soilc, filter_soilc, dtime)
      end if

    end subroutine colpf_summary_acc

    !-------------------------------------------------------------------------------------------------
    subroutine colpf_summaryint_acc(colpf,bounds,num_soilc, filter_soilc, dtime)
      !$acc routine seq
      ! !DESCRIPTION:
      ! bgc interface & pflotran:
      !
      ! !ARGUMENTS:
      type(column_phosphorus_flux)  :: colpf
      type(bounds_type) ,  intent(in) :: bounds
      integer,             intent(in) :: num_soilc       ! number of soil columns in filter
      integer,             intent(in) :: filter_soilc(:) ! filter for soil columns
      real(r8)          ,  intent(in) :: dtime       ! radiation time step (seconds)
      !
      ! !LOCAL VARIABLES:
      integer :: c,j,l       ! indices
      integer :: fc          ! column filter indices

      !-----------------------------------------------------------------------

      ! set time steps
      !#py dtime = real( get_step_size(), r8 )
      if (use_pflotran .and. pf_cmode) then
          ! TODO
      end if

      ! summarize at column-level vertically-resolved littering/removal for PFLOTRAN bgc input needs
      ! first it needs to save the total column-level N rate btw plant pool and decomposible pools at previous time step
      ! for adjusting difference when doing balance check

      do fc = 1,num_soilc
         c = filter_soilc(fc)
         colpf%externalp_to_decomp_delta(c) = 0._r8
         colpf%sminp_net_transport_delta(c)   = 0._r8
         do j = 1, nlevdecomp
            do l = 1, ndecomp_pools
               colpf%externalp_to_decomp_delta(c) =    &
                  colpf%externalp_to_decomp_delta(c) + &
                  colpf%externalp_to_decomp_ppools(c,j,l)*dzsoi_decomp(j)
            end do

            ! sminp leaching/runoff at previous time-step, which may be as source by PFLOTRAN
            colpf%sminp_net_transport_delta(c) = &
               colpf%sminp_net_transport_delta(c) + &
               colpf%sminp_net_transport_vr(c,j)*dzsoi_decomp(j)

         end do
      end do

      colpf%externalp_to_decomp_ppools(:,:,:) = 0._r8
      colpf%sminp_net_transport_vr(:,:) = 0._r8

      ! add up all vertically-resolved addition/removal rates (gC/m3/s) of decomp_pools
      do l = 1, ndecomp_pools
         do j = 1, nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               ! for litter C pools
               if (l==i_met_lit) then
                  colpf%externalp_to_decomp_ppools(c,j,l) =              &
                      colpf%externalp_to_decomp_ppools(c,j,l)            &
                       + colpf%phenology_p_to_litr_met_p(c,j)            &
                       + colpf%dwt_frootp_to_litr_met_p(c,j)             &
                       + colpf%gap_mortality_p_to_litr_met_p(c,j)        &
                       + colpf%harvest_p_to_litr_met_p(c,j)

               elseif (l==i_cel_lit) then
                  colpf%externalp_to_decomp_ppools(c,j,l) =              &
                      colpf%externalp_to_decomp_ppools(c,j,l)            &
                       + colpf%phenology_p_to_litr_cel_p(c,j)            &
                       + colpf%dwt_frootp_to_litr_cel_p(c,j)             &
                       + colpf%gap_mortality_p_to_litr_cel_p(c,j)        &
                       + colpf%harvest_p_to_litr_cel_p(c,j)

               elseif (l==i_lig_lit) then
                  colpf%externalp_to_decomp_ppools(c,j,l) =              &
                      colpf%externalp_to_decomp_ppools(c,j,l)            &
                       + colpf%phenology_p_to_litr_lig_p(c,j)            &
                       + colpf%dwt_frootp_to_litr_lig_p(c,j)             &
                       + colpf%gap_mortality_p_to_litr_lig_p(c,j)        &
                       + colpf%harvest_p_to_litr_lig_p(c,j)

               ! for cwd
               elseif (l==i_cwd) then
                  colpf%externalp_to_decomp_ppools(c,j,l) =              &
                      colpf%externalp_to_decomp_ppools(c,j,l)            &
                       + colpf%dwt_livecrootp_to_cwdp(c,j)               &
                       + colpf%dwt_deadcrootp_to_cwdp(c,j)               &
                       + colpf%gap_mortality_p_to_cwdp(c,j)              &
                       + colpf%harvest_p_to_cwdp(c,j)

               end if

               ! the following is the net changes of plant N to decompible N poools between time-step
               ! in pflotran, decomposible N pools increments ARE from previous time-step (saved above);
               ! while, in CLM-CN all plant N pools are updated with current N fluxes among plant and ground/soil.
               ! therefore, when do balance check it is needed to adjust the time-lag of changes.
               colpf%externalp_to_decomp_delta(c) =     &
                   colpf%externalp_to_decomp_delta(c) - &
                   colpf%externalp_to_decomp_ppools(c,j,l)*dzsoi_decomp(j)

               if (abs(colpf%externalp_to_decomp_ppools(c,j,l))<=1.e-21_r8) then
                   colpf%externalp_to_decomp_ppools(c,j,l) = 0._r8
               end if
            end do ! num_soilc
         end do ! nlevdecomp
      end do ! ndecomp_pools

      ! if pflotran hydrology NOT coupled, need to adjust for sminp leaching, for balance check
      if (.not. pf_hmode) then
         do j = 1, nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               colpf%sminp_net_transport_vr(c,j) = 0._r8

               colpf%sminp_net_transport_delta(c) = &
                   colpf%sminp_net_transport_delta(c) - &
                   colpf%sminp_net_transport_vr(c,j)*dzsoi_decomp(j)
            end do
         end do
      end if

      ! change the sign so that it is the increments from the previous time-step (unit: g/m2/s)
      do fc = 1, num_soilc
         c = filter_soilc(fc)
         colpf%externalp_to_decomp_delta(c)     = -colpf%externalp_to_decomp_delta(c)
         colpf%sminp_net_transport_delta(c)     = -colpf%sminp_net_transport_delta(c)
      end do

    end subroutine colpf_summaryint_acc

    !-----------------------------------------------------------------------
    subroutine colps_summary_acc(colps, bounds, num_soilc, filter_soilc)
      !$acc routine seq
      ! !ARGUMENTS:
      type(column_phosphorus_state) :: colps
      type(bounds_type) , intent(in)  :: bounds
      integer           , intent(in)  :: num_soilc       ! number of soil columns in filter
      integer           , intent(in)  :: filter_soilc(:) ! filter for soil columns
      !
      ! !LOCAL VARIABLES:
      integer  :: c,j,k,l   ! indices
      integer  :: fc       ! lake filter indices
      real(r8) :: maxdepth    ! depth to integrate soil variables
      !-----------------------------------------------------------------------
     do fc = 1,num_soilc
        c = filter_soilc(fc)
        colps%solutionp(c) = 0._r8
        colps%labilep(c)   = 0._r8
        colps%secondp(c)   = 0._r8
        colps%occlp(c)     = 0._r8
        colps%primp(c)     = 0._r8
     end do

     do j = 1, nlevdecomp
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           colps%solutionp(c) = &
                colps%solutionp(c) + &
                colps%solutionp_vr(c,j) * dzsoi_decomp(j)
           colps%labilep(c) = &
                colps%labilep(c) + &
                colps%labilep_vr(c,j) * dzsoi_decomp(j)
           colps%secondp(c) = &
                colps%secondp(c) + &
                colps%secondp_vr(c,j) * dzsoi_decomp(j)
           colps%occlp(c) = &
                colps%occlp(c) + &
                colps%occlp_vr(c,j) * dzsoi_decomp(j)
           colps%primp(c) = &
                colps%primp(c) + &
                colps%primp_vr(c,j) * dzsoi_decomp(j)
        end do
     end do

     ! vertically integrate each of the decomposing P pools
     do l = 1, ndecomp_pools
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           colps%decomp_ppools(c,l) = 0._r8
        end do
        do j = 1, nlevdecomp
           do fc = 1,num_soilc
              c = filter_soilc(fc)
              colps%decomp_ppools(c,l) = &
                   colps%decomp_ppools(c,l) + &
                   colps%decomp_ppools_vr(c,j,l) * dzsoi_decomp(j)
           end do
        end do
     end do

     ! for vertically-resolved soil biogeochemistry, calculate some diagnostics of carbon pools to a given depth
     if ( nlevdecomp > 1) then
        do l = 1, ndecomp_pools
           do fc = 1,num_soilc
              c = filter_soilc(fc)
              colps%decomp_ppools_1m(c,l) = 0._r8
           end do
        end do

        ! vertically integrate each of the decomposing n pools to 1 meter
        maxdepth = 1._r8
        do l = 1, ndecomp_pools
           do j = 1, nlevdecomp
              if ( zisoi(j) <= maxdepth ) then
                 do fc = 1,num_soilc
                    c = filter_soilc(fc)
                    colps%decomp_ppools_1m(c,l) = &
                         colps%decomp_ppools_1m(c,l) + &
                         colps%decomp_ppools_vr(c,j,l) * dzsoi_decomp(j)
                 end do
              elseif ( zisoi(j-1) < maxdepth ) then
                 do fc = 1,num_soilc
                    c = filter_soilc(fc)
                    colps%decomp_ppools_1m(c,l) = &
                         colps%decomp_ppools_1m(c,l) + &
                         colps%decomp_ppools_vr(c,j,l) * (maxdepth - zisoi(j-1))
                 end do
              endif
           end do
        end do

        ! total litter phosphorus to 1 meter (TOTLITN_1m)
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           colps%totlitp_1m(c) = 0._r8
        end do
        do l = 1, ndecomp_pools
           if ( decomp_cascade_con%is_litter(l) ) then
              do fc = 1,num_soilc
                 c = filter_soilc(fc)
                 colps%totlitp_1m(c) = &
                      colps%totlitp_1m(c) + &
                      colps%decomp_ppools_1m(c,l)
              end do
           end if
        end do

        ! total soil organic matter phosphorus to 1 meter (TOTSOMN_1m)
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           colps%totsomp_1m(c) = 0._r8
        end do
        do l = 1, ndecomp_pools
           if ( decomp_cascade_con%is_soil(l) ) then
              do fc = 1,num_soilc
                 c = filter_soilc(fc)
                 colps%totsomp_1m(c) = &
                      colps%totsomp_1m(c) + &
                      colps%decomp_ppools_1m(c,l)
              end do
           end if
        end do

     endif

     ! total litter phosphorus (TOTLITN)
     do fc = 1,num_soilc
        c = filter_soilc(fc)
        colps%totlitp(c)    = 0._r8
     end do
     do l = 1, ndecomp_pools
        if ( decomp_cascade_con%is_litter(l) ) then
           do fc = 1,num_soilc
              c = filter_soilc(fc)
              colps%totlitp(c) = &
                   colps%totlitp(c) + &
                   colps%decomp_ppools(c,l)
           end do
        end if
     end do

     ! total soil organic matter phosphorus (TOTSOMN)
     do fc = 1,num_soilc
        c = filter_soilc(fc)
        colps%totsomp(c)    = 0._r8
     end do
     do l = 1, ndecomp_pools
        if ( decomp_cascade_con%is_soil(l) ) then
           do fc = 1,num_soilc
              c = filter_soilc(fc)
              colps%totsomp(c) = &
                   colps%totsomp(c) + &
                   colps%decomp_ppools(c,l)
           end do
        end if
     end do

     ! total cwdn
     do fc = 1,num_soilc
        c = filter_soilc(fc)
        colps%cwdp(c) = 0._r8
     end do
     do l = 1, ndecomp_pools
        if ( decomp_cascade_con%is_cwd(l) ) then
           do fc = 1,num_soilc
              c = filter_soilc(fc)
              colps%cwdp(c) = &
                   colps%cwdp(c) + &
                   colps%decomp_ppools(c,l)
           end do
        end if
     end do

     ! total sminp
     do fc = 1,num_soilc
        c = filter_soilc(fc)
        colps%sminp(c)      = 0._r8
     end do
     do j = 1, nlevdecomp
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           colps%sminp_vr(c,j) = colps%solutionp_vr(c,j)+ &
                                    colps%labilep_vr(c,j)+ &
                                    colps%secondp_vr(c,j)
        end do
     end do
     do j = 1, nlevdecomp
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           colps%sminp(c) =  colps%sminp(c) + &
           colps%sminp_vr(c,j) * dzsoi_decomp(j)
        end do
     end do

     ! total col_ntrunc
     do fc = 1,num_soilc
        c = filter_soilc(fc)
        colps%ptrunc(c) = 0._r8
     end do
     do j = 1, nlevdecomp
        do fc = 1,num_soilc
           c = filter_soilc(fc)
           colps%ptrunc(c) = &
                colps%ptrunc(c) + &
                colps%ptrunc_vr(c,j) * dzsoi_decomp(j)
        end do
     end do

     do fc = 1,num_soilc
        c = filter_soilc(fc)

        ! total wood product phosphorus
        colps%totprodp(c) = &
             colps%prod1p(c) + &
             colps%prod10p(c) + &
             colps%prod100p(c)

        ! total ecosystem phosphorus, including veg (TOTECOSYSP)
        colps%totecosysp(c) = &
             colps%cwdp(c) + &
             colps%totlitp(c) + &
             colps%totsomp(c) + &
             colps%solutionp(c) + &
             colps%labilep(c) + &
             colps%secondp(c) + &
             colps%primp(c) + &
             colps%occlp(c) + &
             colps%totprodp(c) + &
             colps%totvegp(c)

        ! total column phosphorus, including pft (TOTCOLP)
        colps%totcolp(c) = &
             colps%totpftp(c) + &
             colps%cwdp(c) + &
             colps%totlitp(c) + &
             colps%totsomp(c) + &
             colps%prod1p(c) + &
             colps%solutionp(c) + &
             colps%labilep(c) + &
             colps%secondp(c) + &
             colps%ptrunc(c) + &
             colps%cropseedp_deficit(c)
     end do

   end subroutine colps_summary_acc



     !-----------------------------------------------------------------------
     subroutine vegnf_summary_acc(vegnf, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_nf)
       !
       ! !DESCRIPTION:
       ! Vegetation-level nitrogen flux summary calculations
       !$acc routine seq
       ! !ARGUMENTS:
       type(vegetation_nitrogen_flux)             :: vegnf
       type(bounds_type)           , intent(in)    :: bounds
       integer                     , intent(in)    :: num_soilc       ! number of soil columns in filter
       integer                     , intent(in)    :: filter_soilc(:) ! filter for soil columns
       integer                     , intent(in)    :: num_soilp       ! number of soil patches in filter
       integer                     , intent(in)    :: filter_soilp(:) ! filter for soil patches
       type (column_nitrogen_flux), intent(inout)  :: col_nf          ! column-level nitrogen state for p2c

       !
       ! !LOCAL VARIABLES:
       integer  :: c,p             ! indices
       integer  :: fp              ! filter indices
       !-----------------------------------------------------------------------
       associate(&
         fire_nloss_patch => vegnf%fire_nloss ,&
         fire_nloss_col   => col_nf%fire_nloss_p2c ,&
         wood_harvestn_patch =>  vegnf%wood_harvestn, &
         wood_harvestn_col   => col_nf%wood_harvestn &
         )
       do fp = 1,num_soilp
          p = filter_soilp(fp)

          ! total N deployment (from sminn and retranslocated N pool) (NDEPLOY)
          vegnf%ndeploy(p) = &
               vegnf%sminn_to_npool(p) + &
               vegnf%retransn_to_npool(p)

          ! pft-level wood harvest
          wood_harvestn_patch(p) = &
               vegnf%hrv_deadstemn_to_prod10n(p) + &
               vegnf%hrv_deadstemn_to_prod100n(p)
          if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
               vegnf%wood_harvestn(p) = &
               vegnf%wood_harvestn(p) + &
               vegnf%hrv_cropn_to_prod1n(p)
          end if

          ! total pft-level fire N losses
          fire_nloss_patch(p) = &
               vegnf%m_leafn_to_fire(p)               + &
               vegnf%m_leafn_storage_to_fire(p)       + &
               vegnf%m_leafn_xfer_to_fire(p)          + &
               vegnf%m_frootn_to_fire(p)              + &
               vegnf%m_frootn_storage_to_fire(p)      + &
               vegnf%m_frootn_xfer_to_fire(p)         + &
               vegnf%m_livestemn_to_fire(p)           + &
               vegnf%m_livestemn_storage_to_fire(p)   + &
               vegnf%m_livestemn_xfer_to_fire(p)      + &
               vegnf%m_deadstemn_to_fire(p)           + &
               vegnf%m_deadstemn_storage_to_fire(p)   + &
               vegnf%m_deadstemn_xfer_to_fire(p)      + &
               vegnf%m_livecrootn_to_fire(p)          + &
               vegnf%m_livecrootn_storage_to_fire(p)  + &
               vegnf%m_livecrootn_xfer_to_fire(p)     + &
               vegnf%m_deadcrootn_to_fire(p)          + &
               vegnf%m_deadcrootn_storage_to_fire(p)  + &
               vegnf%m_deadcrootn_xfer_to_fire(p)     + &
               vegnf%m_retransn_to_fire(p)            + &
               vegnf%m_npool_to_fire(p)

         vegnf%gap_nloss_litter(p) = &
              vegnf%m_leafn_to_litter(p)              + &
              vegnf%m_leafn_storage_to_litter(p)      + &
              vegnf%m_leafn_xfer_to_litter(p)         + &
              vegnf%m_frootn_to_litter(p)             + &
              vegnf%m_frootn_storage_to_litter(p)     + &
              vegnf%m_frootn_xfer_to_litter(p)        + &
              vegnf%m_livestemn_to_litter(p)          + &
              vegnf%m_livestemn_storage_to_litter(p)  + &
              vegnf%m_livestemn_xfer_to_litter(p)     + &
              vegnf%m_deadstemn_to_litter(p)          + &
              vegnf%m_deadstemn_storage_to_litter(p)  + &
              vegnf%m_deadstemn_xfer_to_litter(p)     + &
              vegnf%m_livecrootn_to_litter(p)         + &
              vegnf%m_livecrootn_storage_to_litter(p) + &
              vegnf%m_livecrootn_xfer_to_litter(p)    + &
              vegnf%m_deadcrootn_to_litter(p)         + &
              vegnf%m_deadcrootn_storage_to_litter(p) + &
              vegnf%m_deadcrootn_xfer_to_litter(p)    + &
              vegnf%m_retransn_to_litter(p)           + &
              vegnf%m_npool_to_litter(p)

         vegnf%fire_nloss_litter(p) = &
              vegnf%m_deadstemn_to_litter_fire(p)     + &
              vegnf%m_deadcrootn_to_litter_fire(p)    + &
              vegnf%m_retransn_to_litter_fire(p)      + &
              vegnf%m_npool_to_litter_fire(p)         + &
              vegnf%m_leafn_to_litter_fire(p)         + &
              vegnf%m_frootn_to_litter_fire(p)        + &
              vegnf%m_livestemn_to_litter_fire(p)     + &
              vegnf%m_livecrootn_to_litter_fire(p)    + &
              vegnf%m_leafn_storage_to_litter_fire(p) + &
              vegnf%m_frootn_storage_to_litter_fire(p)       + &
              vegnf%m_livestemn_storage_to_litter_fire(p)    + &
              vegnf%m_deadstemn_storage_to_litter_fire(p)    + &
              vegnf%m_livecrootn_storage_to_litter_fire(p)   + &
              vegnf%m_deadcrootn_storage_to_litter_fire(p)   + &
              vegnf%m_leafn_xfer_to_litter_fire(p)           + &
              vegnf%m_frootn_xfer_to_litter_fire(p)          + &
              vegnf%m_livestemn_xfer_to_litter_fire(p)       + &
              vegnf%m_deadstemn_xfer_to_litter_fire(p)       + &
              vegnf%m_livecrootn_xfer_to_litter_fire(p)      + &
              vegnf%m_deadcrootn_xfer_to_litter_fire(p)

         vegnf%hrv_nloss_litter(p) = &
              vegnf%hrv_retransn_to_litter(p)          + &
              vegnf%hrv_npool_to_litter(p)             + &
              vegnf%hrv_leafn_to_litter(p)             + &
              vegnf%hrv_leafn_storage_to_litter(p)     + &
              vegnf%hrv_leafn_xfer_to_litter(p)        + &
              vegnf%hrv_frootn_to_litter(p)            + &
              vegnf%hrv_frootn_storage_to_litter(p)    + &
              vegnf%hrv_frootn_xfer_to_litter(p)       + &
              vegnf%hrv_livestemn_to_litter(p)         + &
              vegnf%hrv_livestemn_storage_to_litter(p) + &
              vegnf%hrv_livestemn_xfer_to_litter(p)    + &
              vegnf%hrv_deadstemn_storage_to_litter(p) + &
              vegnf%hrv_deadstemn_xfer_to_litter(p)    + &
              vegnf%hrv_livecrootn_to_litter(p)        + &
              vegnf%hrv_livecrootn_storage_to_litter(p)+ &
              vegnf%hrv_livecrootn_xfer_to_litter(p)   + &
              vegnf%hrv_deadcrootn_to_litter(p)        + &
              vegnf%hrv_deadcrootn_storage_to_litter(p)+ &
              vegnf%hrv_deadcrootn_xfer_to_litter(p)

          vegnf%sen_nloss_litter(p) = &
              vegnf%livestemn_to_litter(p)            + &
              vegnf%leafn_to_litter(p)                + &
              vegnf%frootn_to_litter(p)
       end do

       call p2c(bounds, num_soilc, filter_soilc, &
            fire_nloss_patch(bounds%begp:bounds%endp)    , &
            fire_nloss_col(bounds%begc:bounds%endc))

       call p2c(bounds, num_soilc, filter_soilc, &
            wood_harvestn_patch(bounds%begp:bounds%endp) , &
            wood_harvestn_col(bounds%begc:bounds%endc))

      end associate
     end subroutine vegnf_summary_acc

     !-----------------------------------------------------------------------
     subroutine vegns_summary_acc(vegns, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_ns)
       !
       ! !DESCRIPTION:
       ! Vegetation-level nitrogen state summary calculations
       !$acc routine seq
       ! !ARGUMENTS:
       type(vegetation_nitrogen_state)            :: vegns
       type(bounds_type)           , intent(in)    :: bounds
       integer                     , intent(in)    :: num_soilc       ! number of soil columns in filter
       integer                     , intent(in)    :: filter_soilc(:) ! filter for soil columns
       integer                     , intent(in)    :: num_soilp       ! number of soil patches in filter
       integer                     , intent(in)    :: filter_soilp(:) ! filter for soil patches
       type (column_nitrogen_state), intent(inout) :: col_ns          ! column-level nitrogen state for p2c

       !
       ! !LOCAL VARIABLES:
       integer  :: c,p             ! indices
       integer  :: fp              ! filter indices
       !-----------------------------------------------------------------------
       associate( &
        plant_n_buffer_patch  => vegns%plant_n_buffer  , &
        plant_n_buffer_col    => col_ns%plant_n_buffer , &
        totvegn_patch  => vegns%totvegn   , &
        totvegn_col    => col_ns%totvegn, &
        totpftn_patch  => vegns%totpftn   , &
        totpftn_col    => col_ns%totpftn, &
        cropseedn_deficit_patch  => vegns%cropseedn_deficit , &
        cropseedn_deficit_col    => col_ns%cropseedn_deficit &
         )
       do fp = 1,num_soilp
          p = filter_soilp(fp)

          ! displayed vegetation nitrogen, excluding storage (DISPVEGN)
          vegns%dispvegn(p) = &
               vegns%leafn(p)      + &
               vegns%frootn(p)     + &
               vegns%livestemn(p)  + &
               vegns%deadstemn(p)  + &
               vegns%livecrootn(p) + &
               vegns%deadcrootn(p)

         ! stored vegetation nitrogen, including retranslocated N pool (STORVEGN)
         vegns%storvegn(p) = &
              vegns%leafn_storage(p)      + &
              vegns%frootn_storage(p)     + &
              vegns%livestemn_storage(p)  + &
              vegns%deadstemn_storage(p)  + &
              vegns%livecrootn_storage(p) + &
              vegns%deadcrootn_storage(p) + &
              vegns%leafn_xfer(p)         + &
              vegns%frootn_xfer(p)        + &
              vegns%livestemn_xfer(p)     + &
              vegns%deadstemn_xfer(p)     + &
              vegns%livecrootn_xfer(p)    + &
              vegns%deadcrootn_xfer(p)    + &
              vegns%npool(p)              + &
              vegns%retransn(p)

         if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
            vegns%dispvegn(p) = &
                 vegns%dispvegn(p) + &
                 vegns%grainn(p)

            vegns%storvegn(p) = &
                 vegns%storvegn(p) + &
                 vegns%grainn_storage(p)     + &
                 vegns%grainn_xfer(p)
         end if

         ! total vegetation nitrogen (TOTVEGN)
         totvegn_patch(p) = &
              vegns%dispvegn(p) + &
              vegns%storvegn(p)

         ! total pft-level carbon (add pft_ntrunc)
         totpftn_patch(p) = &
              totvegn_patch(p) + &
              vegns%ntrunc(p)
      end do ! filtered veg loop

      call p2c(bounds, num_soilc, filter_soilc, &
           plant_n_buffer_patch(bounds%begp:bounds%endp)  , &
           plant_n_buffer_col(bounds%begc:bounds%endc))

      call p2c(bounds, num_soilc, filter_soilc, &
           totvegn_patch(bounds%begp:bounds%endp) , &
           totvegn_col(bounds%begc:bounds%endc))

      call p2c(bounds, num_soilc, filter_soilc, &
           totpftn_patch(bounds%begp:bounds%endp) , &
           totpftn_col(bounds%begc:bounds%endc))

      call p2c(bounds, num_soilc, filter_soilc, &
           cropseedn_deficit_patch(bounds%begp:bounds%endp) , &
           cropseedn_deficit_col(bounds%begc:bounds%endc))

      end associate

     end subroutine vegns_summary_acc


     !-----------------------------------------------------------------------
     subroutine vegpf_summary_acc(vegpf, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_pf)
       !$acc routine seq
       ! !ARGUMENTS:
       type (vegetation_phosphorus_flux) :: vegpf
       type(bounds_type) , intent(in)     :: bounds
       integer           , intent(in)     :: num_soilc       ! number of soil columns in filter
       integer           , intent(in)     :: filter_soilc(:) ! filter for soil columns
       integer           , intent(in)     :: num_soilp       ! number of soil patches in filter
       integer           , intent(in)     :: filter_soilp(:) ! filter for soil patches
       type(column_phosphorus_flux), intent(inout) :: col_pf
       !
       ! !LOCAL VARIABLES:
       integer  :: p, fp   ! indices
       !-----------------------------------------------------------------------
       associate( &
         fire_ploss_patch => vegpf%fire_ploss      ,&
         fire_ploss_col   => col_pf%fire_ploss_p2c ,&
         wood_harvestp_patch  => vegpf%wood_harvestp ,&
         wood_harvestp_col => col_pf%wood_harvestp &
         )
       do fp = 1,num_soilp
          p = filter_soilp(fp)

          ! total P deployment (from sminn and retranslocated P pool) (PDEPLOY)
          vegpf%pdeploy(p) = &
               vegpf%sminp_to_ppool(p) + &
               vegpf%retransp_to_ppool(p)

          ! pft-level wood harvest
          wood_harvestp_patch(p) = &
               vegpf%hrv_deadstemp_to_prod10p(p) + &
               vegpf%hrv_deadstemp_to_prod100p(p)
          if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
               wood_harvestp_patch(p) = &
                wood_harvestp_patch(p) + &
                vegpf%hrv_cropp_to_prod1p(p)
          end if

          ! total pft-level fire P losses
          fire_ploss_patch(p) = &
               vegpf%m_leafp_to_fire(p)               + &
               vegpf%m_leafp_storage_to_fire(p)       + &
               vegpf%m_leafp_xfer_to_fire(p)          + &
               vegpf%m_frootp_to_fire(p)              + &
               vegpf%m_frootp_storage_to_fire(p)      + &
               vegpf%m_frootp_xfer_to_fire(p)         + &
               vegpf%m_livestemp_to_fire(p)           + &
               vegpf%m_livestemp_storage_to_fire(p)   + &
               vegpf%m_livestemp_xfer_to_fire(p)      + &
               vegpf%m_deadstemp_to_fire(p)           + &
               vegpf%m_deadstemp_storage_to_fire(p)   + &
               vegpf%m_deadstemp_xfer_to_fire(p)      + &
               vegpf%m_livecrootp_to_fire(p)          + &
               vegpf%m_livecrootp_storage_to_fire(p)  + &
               vegpf%m_livecrootp_xfer_to_fire(p)     + &
               vegpf%m_deadcrootp_to_fire(p)          + &
               vegpf%m_deadcrootp_storage_to_fire(p)  + &
               vegpf%m_deadcrootp_xfer_to_fire(p)     + &
               vegpf%m_retransp_to_fire(p)            + &
               vegpf%m_ppool_to_fire(p)

         vegpf%gap_ploss_litter(p) = &
              vegpf%m_leafp_to_litter(p)              + &
              vegpf%m_leafp_storage_to_litter(p)      + &
              vegpf%m_leafp_xfer_to_litter(p)         + &
              vegpf%m_frootp_to_litter(p)             + &
              vegpf%m_frootp_storage_to_litter(p)     + &
              vegpf%m_frootp_xfer_to_litter(p)        + &
              vegpf%m_livestemp_to_litter(p)          + &
              vegpf%m_livestemp_storage_to_litter(p)  + &
              vegpf%m_livestemp_xfer_to_litter(p)     + &
              vegpf%m_deadstemp_to_litter(p)          + &
              vegpf%m_deadstemp_storage_to_litter(p)  + &
              vegpf%m_deadstemp_xfer_to_litter(p)     + &
              vegpf%m_livecrootp_to_litter(p)         + &
              vegpf%m_livecrootp_storage_to_litter(p) + &
              vegpf%m_livecrootp_xfer_to_litter(p)    + &
              vegpf%m_deadcrootp_to_litter(p)         + &
              vegpf%m_deadcrootp_storage_to_litter(p) + &
              vegpf%m_deadcrootp_xfer_to_litter(p)    + &
              vegpf%m_retransp_to_litter(p)           + &
              vegpf%m_ppool_to_litter(p)

         vegpf%fire_ploss_litter(p) = &
              vegpf%m_deadstemp_to_litter_fire(p)     + &
              vegpf%m_deadcrootp_to_litter_fire(p)    + &
              vegpf%m_retransp_to_litter_fire(p)      + &
              vegpf%m_ppool_to_litter_fire(p)         + &
              vegpf%m_leafp_to_litter_fire(p)         + &
              vegpf%m_frootp_to_litter_fire(p)        + &
              vegpf%m_livestemp_to_litter_fire(p)     + &
              vegpf%m_livecrootp_to_litter_fire(p)    + &
              vegpf%m_leafp_storage_to_litter_fire(p) + &
              vegpf%m_frootp_storage_to_litter_fire(p)       + &
              vegpf%m_livestemp_storage_to_litter_fire(p)    + &
              vegpf%m_deadstemp_storage_to_litter_fire(p)    + &
              vegpf%m_livecrootp_storage_to_litter_fire(p)   + &
              vegpf%m_deadcrootp_storage_to_litter_fire(p)   + &
              vegpf%m_leafp_xfer_to_litter_fire(p)           + &
              vegpf%m_frootp_xfer_to_litter_fire(p)          + &
              vegpf%m_livestemp_xfer_to_litter_fire(p)       + &
              vegpf%m_deadstemp_xfer_to_litter_fire(p)       + &
              vegpf%m_livecrootp_xfer_to_litter_fire(p)      + &
              vegpf%m_deadcrootp_xfer_to_litter_fire(p)

         vegpf%hrv_ploss_litter(p) = &
              vegpf%hrv_retransp_to_litter(p)         + &
              vegpf%hrv_ppool_to_litter(p)            + &
              vegpf%hrv_leafp_to_litter(p)            + &
              vegpf%hrv_leafp_storage_to_litter(p)    + &
              vegpf%hrv_leafp_xfer_to_litter(p)       + &
              vegpf%hrv_frootp_to_litter(p)           + &
              vegpf%hrv_frootp_storage_to_litter(p)   + &
              vegpf%hrv_frootp_xfer_to_litter(p)      + &
              vegpf%hrv_livestemp_to_litter(p)        + &
              vegpf%hrv_livestemp_storage_to_litter(p)+ &
              vegpf%hrv_livestemp_xfer_to_litter(p)   + &
              vegpf%hrv_deadstemp_storage_to_litter(p)+ &
              vegpf%hrv_deadstemp_xfer_to_litter(p)   + &
              vegpf%hrv_livecrootp_to_litter(p)       + &
              vegpf%hrv_livecrootp_storage_to_litter(p)+ &
              vegpf%hrv_livecrootp_xfer_to_litter(p)  + &
              vegpf%hrv_deadcrootp_to_litter(p)       + &
              vegpf%hrv_deadcrootp_storage_to_litter(p)+ &
              vegpf%hrv_deadcrootp_xfer_to_litter(p)

          vegpf%sen_ploss_litter(p) = &
              vegpf%livestemp_to_litter(p)            + &
              vegpf%leafp_to_litter(p)                + &
              vegpf%frootp_to_litter(p)
       end do

       call p2c(bounds, num_soilc, filter_soilc, &
            fire_ploss_patch(bounds%begp:bounds%endp)     , &
            fire_ploss_col(bounds%begc:bounds%endc) )

       call p2c(bounds, num_soilc, filter_soilc, &
            wood_harvestp_patch(bounds%begp:bounds%endp) , &
            wood_harvestp_col(bounds%begc:bounds%endc))

       end associate

     end subroutine vegpf_summary_acc


     !-----------------------------------------------------------------------
     subroutine vegps_summary_acc (vegps, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_ps)
       !$acc routine seq
       ! !ARGUMENTS:
       type (vegetation_phosphorus_state) :: vegps
       type(bounds_type) , intent(in)      :: bounds
       integer           , intent(in)      :: num_soilc       ! number of soil columns in filter
       integer           , intent(in)      :: filter_soilc(:) ! filter for soil columns
       integer           , intent(in)      :: num_soilp       ! number of soil patches in filter
       integer           , intent(in)      :: filter_soilp(:) ! filter for soil patches
       type(column_phosphorus_state), intent(inout) :: col_ps
       !
       ! !LOCAL VARIABLES:
       integer  :: p        ! indices
       integer  :: fp       ! lake filter indices
       !-----------------------------------------------------------------------
       associate( &
        totvegp_patch  => vegps%totvegp   , &
        totvegp_col    => col_ps%totvegp, &
        totpftp_patch  => vegps%totpftp   , &
        totpftp_col    => col_ps%totpftp, &
        cropseedp_deficit_patch  => vegps%cropseedp_deficit , &
        cropseedp_deficit_col    => col_ps%cropseedp_deficit &
         )
       do fp = 1,num_soilp
          p = filter_soilp(fp)

          ! displayed vegetation phosphorus, excluding storage (DISPVEGN)
          vegps%dispvegp(p) = &
               vegps%leafp(p)      + &
               vegps%frootp(p)     + &
               vegps%livestemp(p)  + &
               vegps%deadstemp(p)  + &
               vegps%livecrootp(p) + &
               vegps%deadcrootp(p)

         ! stored vegetation phosphorus, including retranslocated N pool (STORVEGN)
         vegps%storvegp(p) = &
              vegps%leafp_storage(p)      + &
              vegps%frootp_storage(p)     + &
              vegps%livestemp_storage(p)  + &
              vegps%deadstemp_storage(p)  + &
              vegps%livecrootp_storage(p) + &
              vegps%deadcrootp_storage(p) + &
              vegps%leafp_xfer(p)         + &
              vegps%frootp_xfer(p)        + &
              vegps%livestemp_xfer(p)     + &
              vegps%deadstemp_xfer(p)     + &
              vegps%livecrootp_xfer(p)    + &
              vegps%deadcrootp_xfer(p)    + &
              vegps%ppool(p)              + &
              vegps%retransp(p)

         if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
            vegps%dispvegp(p) = &
                 vegps%dispvegp(p) + &
                 vegps%grainp(p)

            vegps%storvegp(p) = &
                 vegps%storvegp(p) + &
                 vegps%grainp_storage(p)     + &
                 vegps%grainp_xfer(p)
         end if

         ! total vegetation phosphorus (TOTVEGN)
         totvegp_patch(p) = &
              vegps%dispvegp(p) + &
              vegps%storvegp(p)

         ! total pft-level carbon (add pft_ntrunc)
         totpftp_patch(p) = &
              totvegp_patch(p) + &
              vegps%ptrunc(p)

      end do

      call p2c(bounds, num_soilc, filter_soilc, &
           totvegp_patch(bounds%begp:bounds%endp)  , &
           totvegp_col(bounds%begc:bounds%endc) )

      call p2c(bounds, num_soilc, filter_soilc, &
           totpftp_patch(bounds%begp:bounds%endp) , &
           totpftp_col(bounds%begc:bounds%endc) )

      call p2c(bounds, num_soilc, filter_soilc, &
           cropseedp_deficit_patch(bounds%begp:bounds%endp) , &
           cropseedp_deficit_col(bounds%begc:bounds%endc) )
      end associate

     end subroutine vegps_summary_acc

     subroutine elm_zero_fluxes(bounds)
       !$acc routine seq
       use shr_kind_mod , only : r8 => shr_kind_r8
       use clm_varctl   , only : use_cn, use_c13, use_c14
       use clm_varpar   , only : nlevdecomp_full
       use decompMod    , only : bounds_type
       use VegetationDataType, only : veg_cs, veg_ns, veg_ps
       use ColumnDataType    , only : col_cf, col_nf, col_pf
       use ColumnDataType    , only : c13_col_cf
       use ColumnDataType    , only : c14_col_cf
       use GridcellDataType  , only : grc_cf, grc_nf, grc_pf
       use GridcellDataType  , only : c13_grc_cf
       use GridcellDataType  , only : c14_grc_cf

       implicit none

       type(bounds_type), intent(in) :: bounds
       integer :: g, c,j,p


       if (use_cn) then
         do p = bounds%begp, bounds%endp
          veg_cs%dispvegc(p)  = 0._r8
          veg_cs%storvegc(p)  = 0._r8
          veg_cs%totpftc (p)  = 0._r8

          veg_ps%dispvegp(p)  = 0._r8
          veg_ps%storvegp(p)  = 0._r8
          veg_ps%totvegp (p)  = 0._r8
          veg_ps%totpftp (p)  = 0._r8

          veg_ns%dispvegn(p)  = 0._r8
          veg_ns%storvegn(p)  = 0._r8
          veg_ns%totvegn (p)  = 0._r8
          veg_ns%totpftn (p)  = 0._r8
        end do

          !gridcell
          do g = bounds%begg, bounds%endg
             grc_nf%dwt_seedn_to_leaf(g)     = 0._r8
             grc_nf%dwt_seedn_to_deadstem(g) = 0._r8
             grc_nf%dwt_conv_nflux(g)        = 0._r8
             grc_nf%dwt_seedn_to_npool(g)    = 0._r8
             grc_nf%dwt_prod10n_gain(g)      = 0._r8
             grc_nf%dwt_prod100n_gain(g)     = 0._r8

             grc_cf%dwt_seedc_to_leaf(g)         = 0._r8
             grc_cf%dwt_seedc_to_deadstem(g)     = 0._r8
             grc_cf%dwt_conv_cflux(g)            = 0._r8
             grc_cf%dwt_prod10c_gain(g)          = 0._r8
             grc_cf%dwt_prod100c_gain(g)         = 0._r8
             grc_cf%hrv_deadstemc_to_prod10c(g)  = 0._r8
             grc_cf%hrv_deadstemc_to_prod100c(g) = 0._r8

             grc_pf%dwt_seedp_to_leaf(g)     = 0._r8
             grc_pf%dwt_seedp_to_deadstem(g) = 0._r8
             grc_pf%dwt_conv_pflux(g)        = 0._r8
             grc_pf%dwt_seedp_to_ppool(g)    = 0._r8
             grc_pf%dwt_prod10p_gain(g)      = 0._r8
             grc_pf%dwt_prod100p_gain(g)     = 0._r8
          end do

          !ColumnDataTypes
          do c = bounds%begc,bounds%endc
             col_cf%dwt_conv_cflux(c)           = 0._r8
             col_cf%dwt_prod10c_gain(c)         = 0._r8
             col_cf%dwt_prod100c_gain(c)        = 0._r8
             col_cf%dwt_slash_cflux(c)          = 0._r8

             col_nf%dwt_conv_nflux(c)        = 0._r8
             col_nf%dwt_prod10n_gain(c)      = 0._r8
             col_nf%dwt_prod100n_gain(c)     = 0._r8
             col_nf%dwt_slash_nflux(c)       = 0._r8

             col_pf%dwt_conv_pflux(c)        = 0._r8
             col_pf%dwt_prod10p_gain(c)      = 0._r8
             col_pf%dwt_prod100p_gain(c)     = 0._r8
             col_pf%dwt_slash_pflux(c)       = 0._r8
          end do

          do j = 1, nlevdecomp_full
             do c = bounds%begc,bounds%endc
                col_cf%dwt_frootc_to_litr_met_c(c,j)    = 0._r8
                col_cf%dwt_frootc_to_litr_cel_c(c,j)    = 0._r8
                col_cf%dwt_frootc_to_litr_lig_c(c,j)    = 0._r8
                col_cf%dwt_livecrootc_to_cwdc(c,j)      = 0._r8
                col_cf%dwt_deadcrootc_to_cwdc(c,j)      = 0._r8

                col_nf%dwt_frootn_to_litr_met_n(c,j) = 0._r8
                col_nf%dwt_frootn_to_litr_cel_n(c,j) = 0._r8
                col_nf%dwt_frootn_to_litr_lig_n(c,j) = 0._r8
                col_nf%dwt_livecrootn_to_cwdn(c,j)   = 0._r8
                col_nf%dwt_deadcrootn_to_cwdn(c,j)   = 0._r8

                col_pf%dwt_frootp_to_litr_met_p(c,j) = 0._r8
                col_pf%dwt_frootp_to_litr_cel_p(c,j) = 0._r8
                col_pf%dwt_frootp_to_litr_lig_p(c,j) = 0._r8
                col_pf%dwt_livecrootp_to_cwdp(c,j)   = 0._r8
                col_pf%dwt_deadcrootp_to_cwdp(c,j)   = 0._r8
             end do
          end do

          if (use_c13) then
            do g = bounds%begg, bounds%endg
               c13_grc_cf%dwt_seedc_to_leaf(g)         = 0._r8
               c13_grc_cf%dwt_seedc_to_deadstem(g)     = 0._r8
               c13_grc_cf%dwt_conv_cflux(g)            = 0._r8
               c13_grc_cf%dwt_prod10c_gain(g)          = 0._r8
               c13_grc_cf%dwt_prod100c_gain(g)         = 0._r8
               c13_grc_cf%hrv_deadstemc_to_prod10c(g)  = 0._r8
               c13_grc_cf%hrv_deadstemc_to_prod100c(g) = 0._r8
            end do
            do c = bounds%begc,bounds%endc
               c13_col_cf%dwt_conv_cflux(c)           = 0._r8
               c13_col_cf%dwt_prod10c_gain(c)         = 0._r8
               c13_col_cf%dwt_prod100c_gain(c)        = 0._r8
               c13_col_cf%dwt_slash_cflux(c)          = 0._r8
            end do
            !!
            do j = 1, nlevdecomp_full
               do c = bounds%begc,bounds%endc
                  c13_col_cf%dwt_frootc_to_litr_met_c(c,j)    = 0._r8
                  c13_col_cf%dwt_frootc_to_litr_cel_c(c,j)    = 0._r8
                  c13_col_cf%dwt_frootc_to_litr_lig_c(c,j)    = 0._r8
                  c13_col_cf%dwt_livecrootc_to_cwdc(c,j)      = 0._r8
                  c13_col_cf%dwt_deadcrootc_to_cwdc(c,j)      = 0._r8
               end do
            end do
            !!
          end if

          if (use_c14) then
            do g = bounds%begg, bounds%endg
               c14_grc_cf%dwt_seedc_to_leaf(g)         = 0._r8
               c14_grc_cf%dwt_seedc_to_deadstem(g)     = 0._r8
               c14_grc_cf%dwt_conv_cflux(g)            = 0._r8
               c14_grc_cf%dwt_prod10c_gain(g)          = 0._r8
               c14_grc_cf%dwt_prod100c_gain(g)         = 0._r8
               c14_grc_cf%hrv_deadstemc_to_prod10c(g)  = 0._r8
               c14_grc_cf%hrv_deadstemc_to_prod100c(g) = 0._r8
            end do
            do c = bounds%begc,bounds%endc
               c14_col_cf%dwt_conv_cflux(c)           = 0._r8
               c14_col_cf%dwt_prod10c_gain(c)         = 0._r8
               c14_col_cf%dwt_prod100c_gain(c)        = 0._r8
               c14_col_cf%dwt_slash_cflux(c)          = 0._r8
            end do

            do j = 1, nlevdecomp_full
               do c = bounds%begc,bounds%endc
                  c14_col_cf%dwt_frootc_to_litr_met_c(c,j)    = 0._r8
                  c14_col_cf%dwt_frootc_to_litr_cel_c(c,j)    = 0._r8
                  c14_col_cf%dwt_frootc_to_litr_lig_c(c,j)    = 0._r8
                  c14_col_cf%dwt_livecrootc_to_cwdc(c,j)      = 0._r8
                  c14_col_cf%dwt_deadcrootc_to_cwdc(c,j)      = 0._r8
               end do
            end do

          end if

        end if

     end subroutine elm_zero_fluxes

     subroutine col_wf_Reset(col_wf,num_nolakec,filter_nolakec)
       !$acc routine seq

       implicit none
       type(column_water_flux) , intent(inout) :: col_wf
       integer, value, intent(in) :: num_nolakec
       integer, intent(in) :: filter_nolakec(:)

       integer :: fc, c

       do fc = 1, num_nolakec
         c = filter_nolakec(fc)
         col_wf%qflx_snow2topsoi     (c)   = 0._r8
         col_wf%qflx_h2osfc2topsoi   (c)   = 0._r8
       enddo

     end subroutine col_wf_Reset

     subroutine set_fracsno(bounds)
       !$acc routine seq
       use ColumnType  , only : col_pp
       use LandunitType , only : lun_pp
       use ColumnDataType , only : col_ws

       type(bounds_type), intent(in) :: bounds

       integer :: c, l

       ! ============================================================================
       ! ! Fraction of soil covered by snow (Z.-L. Yang U. Texas)
       ! ============================================================================

       do c = bounds%begc,bounds%endc
          l = col_pp%landunit(c)
          if (lun_pp%urbpoi(l)) then
             ! Urban landunit use Bonan 1996 (LSM Technical Note)
             col_ws%frac_sno(c) = min( col_ws%snow_depth(c)/0.05_r8, 1._r8)
          end if
       end do
     end subroutine set_fracsno

     subroutine crop_vars_CropIncrementYear(num_pcropp,filter_pcropp,crop_vars)

       !$acc routine seq
       use CropType,    only :  crop_type
       implicit none

       integer, value, intent(in) :: num_pcropp
       integer, intent(in)        :: filter_pcropp(:)
       type(crop_type), intent(inout)     :: crop_vars

       integer :: fp , p

       do fp = 1, num_pcropp
          p = filter_pcropp(fp)
          crop_vars%nyrs_crop_active_patch(p) = crop_vars%nyrs_crop_active_patch(p) + 1
       end do

     end subroutine crop_vars_CropIncrementYear

end module Method_procs_acc
