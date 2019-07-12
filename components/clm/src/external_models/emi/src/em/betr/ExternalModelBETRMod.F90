module ExternalModelBETRMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module provides wrapper for BeTR in ALM
  !
  use abortutils                   , only : endrun
  use shr_kind_mod                 , only : r8 => shr_kind_r8
  use shr_log_mod                  , only : errMsg => shr_log_errMsg
  use EMI_DataMod                  , only : emi_data_list, emi_data
  use ColumnType                   , only : col_pp
  use ColumnDataType               , only : col_cs, col_cf
  use VegetationType               , only : veg_pp
  use decompMod                    , only : bounds_type
  use clm_varctl                   , only : use_c13, use_c14
  use ExternalModelConstants
  use EMI_Atm2LndType_Constants
  use EMI_CanopyStateType_Constants
  use EMI_ColumnType_Constants
  use EMI_EnergyFluxType_Constants
  use EMI_Filter_Constants
  use EMI_Landunit_Constants
  use EMI_SoilHydrologyType_Constants
  use EMI_SoilStateType_Constants
  use EMI_TemperatureType_Constants
  use EMI_WaterFluxType_Constants
  use EMI_WaterStateType_Constants
  use EMI_KineticParType_Constants
  use ExternalModelBaseType        , only : em_base_type
  use BeTRSimulationALM            , only : betr_simulation_alm_type
  !
  implicit none
  private
  !
  !required for initialization
  type, public, extends(em_base_type) :: em_betr_type

    !biophysical state variables
    integer :: index_l2e_state_frac_h2osfc_1d
    integer :: index_l2e_state_finundated_1d
    integer :: index_l2e_state_h2osoi_liq_2d
    integer :: index_l2e_state_h2osoi_ice_2d
    integer :: index_l2e_state_h2osoi_liqvol_2d
    integer :: index_l2e_state_h2osoi_icevol_2d
    integer :: index_l2e_state_h2osoi_vol_2d
    integer :: index_l2e_state_air_vol_2d
    integer :: index_l2e_state_rho_vap_2d
    integer :: index_l2e_state_rhvap_soi_2d
    integer :: index_l2e_state_smp_l_2d
    integer :: index_l2e_state_tsoi_10cm_1d
    integer :: index_l2e_state_tsoisno_2d
    integer :: index_l2e_state_t_veg_1d
    integer :: index_l2e_state_fracice_2d
    integer :: index_l2e_state_forc_pbot_downscaled_1d
    integer :: index_l2e_state_t_downscaled_1d
    integer :: index_l2e_state_altmax_1d
    integer :: index_l2e_state_altmax_lastyear_1d
    integer :: index_l2e_state_lbl_rsc_h2o_veg_1d
    integer :: index_l2e_state_elai_veg_1d
    integer :: index_l2e_state_soil_pH_2d
    integer :: index_l2e_state_bsw_2d
    integer :: index_l2e_state_watsat_2d
    integer :: index_l2e_state_eff_porosity_2d
    integer :: index_l2e_state_cellorg_2d
    integer :: index_l2e_state_cellclay_2d
    integer :: index_l2e_state_cellsand_2d
    integer :: index_l2e_state_bd_2d
    integer :: index_l2e_state_watfc_2d
    integer :: index_l2e_state_suscsat_2d
    integer :: index_l2e_state_rootfr_veg_2d

    !bgc state
    integer :: index_l2e_bgc_coeff_t_scalar_2d
    integer :: index_l2e_bgc_coeff_w_scalar_2d
    integer :: index_l2e_bgc_coeff_decompk_3d

    !biophysical flux variables
    integer :: index_l2e_flux_annsum_npp_1d
    integer :: index_l2e_flux_agnpp_1d
    integer :: index_l2e_flux_bgnpp_1d
    integer :: index_l2e_flux_tempavg_agnpp_1d
    integer :: index_l2e_flux_tempavg_bgnpp_1d
    integer :: index_l2e_flux_annavg_agnpp_1d
    integer :: index_l2e_flux_annavg_bgnpp_1d
    integer :: index_l2e_flux_qflx_infl_1d
    integer :: index_l2e_flux_qflx_totdrain_1d
    integer :: index_l2e_flux_qflx_gross_evap_soil_1d
    integer :: index_l2e_flux_qflx_gross_infl_soil_1d
    integer :: index_l2e_flux_qflx_surf_1d
    integer :: index_l2e_flux_qflx_dew_ground_1d
    integer :: index_l2e_flux_qflx_dew_snow_1d
    integer :: index_l2e_flux_qflx_sub_snow_1d
    integer :: index_l2e_flux_qflx_sub_snow_vol_1d
    integer :: index_l2e_flux_qflx_h2osfc2topsoil_1d
    integer :: index_l2e_flux_qflx_snow2topsoil_1d
    integer :: index_l2e_flux_qflx_rootsoi_2d
    integer :: index_l2e_flux_qflx_runoff_1d
    integer :: index_l2e_flux_qflx_drain_2d
    integer :: index_l2e_flux_qflx_tran_veg_1d
    integer :: index_l2e_flux_qflx_rootsoi_vegfrac_2d
    integer :: index_l2e_flux_qflx_bot_1d
    integer :: index_e2l_flux_qflx_infl_1d
    integer :: index_e2l_flux_qflx_adv_2d
    integer :: index_e2l_flux_totdrain_1d
    integer :: index_e2l_flux_gross_evap_soil_1d
    integer :: index_e2l_flux_gross_infl_soil_1d
    integer :: index_e2l_flux_drain_2d

    !carbon states
    integer :: index_l2e_bgc_state_c12_decomp_pools_3d
    integer :: index_l2e_bgc_state_c13_decomp_pools_3d
    integer :: index_l2e_bgc_state_c14_decomp_pools_3d
    integer :: index_e2l_bgc_state_c12_decomp_pools_3d
    integer :: index_e2l_bgc_state_c13_decomp_pools_3d
    integer :: index_e2l_bgc_state_c14_decomp_pools_3d
    integer :: index_e2l_bgc_state_c12_cwd_1d
    integer :: index_e2l_bgc_state_c12_totlit_1d
    integer :: index_e2l_bgc_state_c12_totsom_1d
    integer :: index_e2l_bgc_state_c12_totlit_1m_1d
    integer :: index_e2l_bgc_state_c12_totsom_1m_1d
    integer :: index_e2l_bgc_state_c13_cwd_1d
    integer :: index_e2l_bgc_state_c13_totlit_1d
    integer :: index_e2l_bgc_state_c13_totsom_1d
    integer :: index_e2l_bgc_state_c13_totlit_1m_1d
    integer :: index_e2l_bgc_state_c13_totsom_1m_1d
    integer :: index_e2l_bgc_state_c14_cwd_1d
    integer :: index_e2l_bgc_state_c14_totlit_1d
    integer :: index_e2l_bgc_state_c14_totsom_1d
    integer :: index_e2l_bgc_state_c14_totlit_1m_1d
    integer :: index_e2l_bgc_state_c14_totsom_1m_1d

    !carbon fluxes
    integer :: index_l2e_bgc_flux_c12_litr_met_2d
    integer :: index_l2e_bgc_flux_c12_litr_cel_2d
    integer :: index_l2e_bgc_flux_c12_litr_lig_2d
    integer :: index_l2e_bgc_flux_c12_litr_cwd_2d
    integer :: index_l2e_bgc_flux_c13_litr_met_2d
    integer :: index_l2e_bgc_flux_c13_litr_cel_2d
    integer :: index_l2e_bgc_flux_c13_litr_lig_2d
    integer :: index_l2e_bgc_flux_c13_litr_cwd_2d
    integer :: index_l2e_bgc_flux_c14_litr_met_2d
    integer :: index_l2e_bgc_flux_c14_litr_cel_2d
    integer :: index_l2e_bgc_flux_c14_litr_lig_2d
    integer :: index_l2e_bgc_flux_c14_litr_cwd_2d
    integer :: index_e2l_bgc_flux_c12_hr_1d
    integer :: index_e2l_bgc_flux_c13_hr_1d
    integer :: index_e2l_bgc_flux_c14_hr_1d
    integer :: index_e2l_bgc_flux_c12_fire_decomp_loss_1d
    integer :: index_e2l_bgc_flux_c13_fire_decomp_loss_1d
    integer :: index_e2l_bgc_flux_c14_fire_decomp_loss_1d
    integer :: index_e2l_bgc_flux_c12_som_leached_1d
    integer :: index_e2l_bgc_flux_c12_som_runoff_1d
    integer :: index_e2l_bgc_flux_c12_dic_leached_1d
    integer :: index_e2l_bgc_flux_c12_dic_runoff_1d

    !nitrogen states
    integer :: index_l2e_bgc_state_n14_decomp_pools_3d
    integer :: index_e2l_bgc_state_n14_decomp_pools_3d
    integer :: index_l2e_bgc_state_n14_nh4_soil_2d
    integer :: index_l2e_bgc_state_n14_no3_soil_2d
    integer :: index_e2l_bgc_state_n14_nh4_soil_2d
    integer :: index_e2l_bgc_state_n14_no3_soil_2d
    integer :: index_e2l_bgc_state_n14_cwd_1d
    integer :: index_e2l_bgc_state_n14_totlit_1d
    integer :: index_e2l_bgc_state_n14_totsom_1d
    integer :: index_e2l_bgc_state_n14_totlit_1m_1d
    integer :: index_e2l_bgc_state_n14_totsom_1m_1d

    !nitrogen fluxes
    integer :: index_l2e_bgc_flux_n14_litr_met_2d
    integer :: index_l2e_bgc_flux_n14_litr_cel_2d
    integer :: index_l2e_bgc_flux_n14_litr_lig_2d
    integer :: index_l2e_bgc_flux_n14_litr_cwd_2d
    integer :: index_l2e_bgc_flux_n14_nh4_fix_soil_2d
    integer :: index_l2e_bgc_flux_n14_nh4_atmdep_soil_2d
    integer :: index_l2e_bgc_flux_n14_no3_atmdep_soil_2d
    integer :: index_l2e_bgc_flux_n14_nh4_fert_soil_2d
    integer :: index_l2e_bgc_flux_n14_no3_fert_soil_2d
    integer :: index_e2l_bgc_flux_n14_fire_decomp_loss_1d
    integer :: index_e2l_bgc_flux_n14_som_leached_1d
    integer :: index_e2l_bgc_flux_n14_som_runoff_1d
    integer :: index_e2l_bgc_flux_n14_din_leached_1d
    integer :: index_e2l_bgc_flux_n14_din_runoff_1d
    integer :: index_e2l_bgc_flux_n14_nh4_leached_1d
    integer :: index_e2l_bgc_flux_h14_nh4_runoff_1d
    integer :: index_e2l_bgc_flux_n14_no3_leached_1d
    integer :: index_e2l_bgc_flux_h14_no3_runoff_1d
    integer :: index_e2l_bgc_flux_n14_nh3_soil_emi_1d
    integer :: index_e2l_bgc_flux_n14_supplement_sminn_1d
    integer :: index_e2l_bgc_flux_n14_nitrif_1d
    integer :: index_e2l_bgc_flux_n14_denitrif_1d
    integer :: index_e2l_bgc_flux_n14_n2o_nitrif_1d
    integer :: index_e2l_bgc_flux_n14_n2o_denitrif_1d
    integer :: index_e2l_bgc_flux_n14_smin_nh4_to_plant_veg_1d
    integer :: index_e2l_bgc_flux_n14_smin_no3_to_plant_veg_1d
    integer :: index_e2l_bgc_flux_n14_sminn_to_plant_veg_1d

    !phosphorus states
    integer :: index_l2e_bgc_state_p31_decomp_pools_3d
    integer :: index_e2l_bgc_state_p31_decomp_pools_3d
    integer :: index_l2e_bgc_state_p31_solutionp_soil_2d
    integer :: index_e2l_bgc_state_p31_solutionp_soil_2d
    integer :: index_e2l_bgc_state_p31_cwd_1d
    integer :: index_e2l_bgc_state_p31_totlit_1d
    integer :: index_e2l_bgc_state_p31_totsom_1d
    integer :: index_e2l_bgc_state_p31_totlit_1m_1d
    integer :: index_e2l_bgc_state_p31_totsom_1m_1d
    integer :: index_e2l_bgc_flux_p31_sminp_to_plant_veg_1d
    integer :: index_e2l_bgc_flux_p31_col_plant_demand_2d
    integer :: index_e2l_bgc_flux_p31_adsorb_to_labile_2d

    !phosphorus fluxes
    integer :: index_l2e_bgc_flux_p31_litr_met_2d
    integer :: index_l2e_bgc_flux_p31_litr_cel_2d
    integer :: index_l2e_bgc_flux_p31_litr_lig_2d
    integer :: index_l2e_bgc_flux_p31_litr_cwd_2d
    integer :: index_l2e_bgc_flux_p31_sminp_fert_soil_2d
    integer :: index_l2e_bgc_flux_p31_sminp_atmdep_soil_2d
    integer :: index_e2l_bgc_flux_p31_fire_decomp_loss_1d
    integer :: index_e2l_bgc_flux_p31_som_leached_1d
    integer :: index_e2l_bgc_flux_p31_som_runoff_1d
    integer :: index_e2l_bgc_flux_p31_dip_leached_1d
    integer :: index_e2l_bgc_flux_p31_dip_runoff_1d
    integer :: index_e2l_bgc_flux_p31_supplement_sminp_1d

    !kinetic parameters for nutrient coupling
    integer :: index_l2e_bgc_par_plant_nh4_vmax_2d
    integer :: index_l2e_bgc_par_plant_no3_vmax_2d
    integer :: index_l2e_bgc_par_plant_p_vmax_2d
    integer :: index_l2e_bgc_par_plant_nh4_km_2d
    integer :: index_l2e_bgc_par_plant_no3_km_2d
    integer :: index_l2e_bgc_par_plant_p_km_2d
    integer :: index_l2e_bgc_par_plant_eff_ncompetb_2d
    integer :: index_l2e_bgc_par_plant_eff_pcompetb_2d
    integer :: index_l2e_bgc_par_plant_eff_frootc_2d
    integer :: index_l2e_bgc_par_minsurf_p_compet_2d
    integer :: index_l2e_bgc_par_minsurf_nh4_compet_2d
    integer :: index_l2e_bgc_par_km_minsurf_nh4_2d
    integer :: index_l2e_bgc_par_km_minsurf_p_2d

    integer :: index_l2e_bgc_par_decomp_eff_ncompetb_2d
    integer :: index_l2e_bgc_par_decomp_eff_pcompetb_2d
    integer :: index_l2e_bgc_par_nitrif_eff_ncompetb_2d
    integer :: index_l2e_bgc_par_denitrif_eff_ncompetb_2d
    integer :: index_l2e_bgc_par_km_nitrif_2d
    integer :: index_l2e_bgc_par_km_denitrif_2d
    integer :: index_l2e_bgc_par_dsolutionp_dt_2d
    integer :: index_l2e_bgc_par_vmax_minsurf_p_2d
    integer :: index_l2e_bgc_par_dlabp_dt_2d

    integer :: index_l2e_filter_nolakec
    integer :: index_l2e_filter_num_nolakec
    class(betr_simulation_alm_type), pointer :: ep_betr
  contains
     procedure, public :: Populate_L2E_Init_List  => EM_BETR_Populate_L2E_Init_List
     procedure, public :: Populate_L2E_List       => EM_BETR_Populate_L2E_List
     procedure, public :: Populate_E2L_List       => EM_BETR_Populate_E2L_List
     procedure, public :: Init                    => EM_BETR_Init
     procedure, public :: Solve                   => EM_BETR_Solve
  end type em_BETR_type

  ! EM_BeTR_CalcSmpL

contains

    subroutine EM_BETR_Init(this, l2e_init_list, e2l_init_list, iam, bounds_clump)

    implicit none
    !
    ! !ARGUMENTS:
    class(em_betr_type)                  :: this
    class(emi_data_list) , intent(in)    :: l2e_init_list
    class(emi_data_list) , intent(inout) :: e2l_init_list
    integer              , intent(in)    :: iam
    type(bounds_type)    , intent(in)    :: bounds_clump


    end subroutine EM_BETR_Init
    !------------------------------------------------------------------------
  subroutine EM_BETR_Solve(this, em_stage, dt, nstep, clump_rank, l2e_list, e2l_list, &
       bounds_clump)
    !
    ! !DESCRIPTION:
    !
    !
    ! !USES:
    use shr_kind_mod              , only : r8 => shr_kind_r8
    use abortutils                , only : endrun
    use shr_log_mod               , only : errMsg => shr_log_errMsg
    use clm_varctl                , only : iulog
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_betr_type)                  :: this
    integer              , intent(in)    :: em_stage
    real(r8)             , intent(in)    :: dt
    integer              , intent(in)    :: nstep
    integer              , intent(in)    :: clump_rank
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    type(bounds_type)    , intent (in)   :: bounds_clump

    select case(em_stage)

    case (EM_BETR_BEGIN_MASS_BALANCE_STAGE)
       call EM_BETR_BeginMassBalance_Solve(dt, nstep, bounds)

    case (EM_BETR_END_MASS_BALANCE_STAGE)
       call EM_BETR_EndMassBalance_Solve(bounds)

    case (EM_BETR_PRE_DIAG_WATER_FLUX_STAGE)
       call EM_BETR_PreDiagSoilColWaterFlux_Solve(bounds, l2e_list, e2l_list)

    case (EM_BETR_STEP_WITHOUT_DRAINAGE_STAGE)
       call EM_BETR_StepWithoutDraiange(bounds, l2e_list, e2l_list)

    case (EM_BETR_STEP_WITH_DRAINAGE_STAGE)
       call EM_BETR_StepWithDraiange(bounds, l2e_list, e2l_list)

    case (EM_BETR_OUTLOOP_SOILBGC_STAGE)
       call EM_BeTR_OutLoopSoilBGC()

    case (EM_BETR_CALC_DEWSUB_FLUX_STAGE)
       call EM_BeTR_CalcDewSubFLux()

    case default
       write(iulog,*)'EM_BETR_Solve: Unknown em_stage.'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine EM_BETR_Solve

  !------------------------------------------------------------------------
  subroutine EM_BETR_Populate_L2E_Init_List(this, l2e_init_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by BeTR from ALM
    !
    ! !USES:
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_betr_type)                 :: this
    class(emi_data_list), intent(inout) :: l2e_init_list
    !
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index


  end subroutine EM_BETR_Populate_L2E_Init_List
  !------------------------------------------------------------------------
  subroutine EM_BETR_Populate_L2E_List(this, l2e_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by BeTR from ALM
    !
    ! !USES:
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_betr_type)                 :: this
    class(emi_data_list), intent(inout) :: l2e_list
    !
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index

    number_em_stages = 5
    allocate(em_stages(number_em_stages))

    em_stages(1) = EM_BETR_PRE_DIAG_WATER_FLUX_STAGE

    em_stages(2) = EM_BETR_STEP_WITHOUT_DRAINAGE_STAGE

    em_stages(3) = EM_BETR_STEP_WITH_DRAINAGE_STAGE

    em_stages(4) = EM_BETR_OUTLOOP_SOILBGC_STAGE

    em_stages(5) = EM_BETR_CALC_DEWSUB_FLUX_STAGE

    id                            = L2E_STATE_FRAC_H2OSFC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_frac_h2osfc_1d   = index

    id                            = L2E_STATE_FRAC_INUNDATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_finundated_1d    = index

    id                            = L2E_STATE_H2OSOI_LIQ_NLEVSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_liq_2d    = index

    id                            = L2E_STATE_H2OSOI_ICE_NLEVSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_ice_2d    = index

    id                            = L2E_STATE_H2OSOI_LIQ_VOL_NLEVSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_liqvol_2d = index

    id                            = L2E_STATE_H2OSOI_ICE_VOL_NLEVSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_icevol_2d = index

    id                            = L2E_STATE_H2OSOI_VOL_NLEVSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_h2osoi_vol_2d    = index

    id                            = L2E_STATE_AIR_VOL_NLEVSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_air_vol_2d       = index

    id                            = L2E_STATE_RHO_VAP_NLEVSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_rho_vap_2d       = index

    id                            = L2E_STATE_RHVAP_SOI_NLEVSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_rhvap_soi_2d     = index

    id                            = L2E_STATE_SOIL_MATRIC_POTENTIAL_NLEVSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_smp_l_2d         = index

    id                            = L2E_FILTER_NOLAKEC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_filter_nolakec      = index

    id                            = L2E_FILTER_NUM_NOLAKEC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_filter_num_nolakec  = index

    id                            = L2E_STATE_TSOI10CM
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_tsoi_10cm_1d = index

    id                            = L2E_STATE_TSOIL_NLEVSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_tsoisno_2d = index

    id                            = L2E_STATE_TVEG
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_t_veg_1d = index

    id                            = L2E_STATE_FRACICE
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_fracice_2d = index

    id                            = L2E_STATE_FORC_PBOT_DOWNSCALED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_forc_pbot_downscaled_1d = index

    id                            = L2E_STATE_FORC_T_DOWNSCALED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_t_downscaled_1d = index

    id                            = L2E_STATE_ALTMAX
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_altmax_1d = index

    id                            = L2E_STATE_ALTMAX_LASTYEAR
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_altmax_lastyear_1d = index

    id                            = L2E_STATE_LBL_RSC_H2O
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_lbl_rsc_h2o_veg_1d = index

    id                            = L2E_STATE_ELAI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_elai_veg_1d = index

    id                            = L2E_STATE_SOIL_PH
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_soil_pH_2d = index

    id                            = L2E_PARAMETER_BSWC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_bsw_2d = index

    id                            = L2E_PARAMETER_WATSATC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_watsat_2d = index

    id                            = L2E_PARAMETER_EFFPOROSITYC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_eff_porosity_2d = index

    id                            = L2E_PARAMETER_CELLORG
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_cellorg_2d = index

    id                            = L2E_PARAMETER_CELLCLAY
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_cellclay_2d = index

    id                            = L2E_PARAMETER_CELLSAND
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_cellsand_2d = index

    id                            = L2E_PARAMETER_BD
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_bd_2d = index

    id                            = L2E_PARAMETER_WATFC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_watfc_2d = index

    id                            = L2E_PARAMETER_SUCSATC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_suscsat_2d = index

    id                            = L2E_PARAMETER_ROOTFR_PATCH
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_state_rootfr_veg_2d = index

    !biophysical flux variables
!    id                            = L2E_FLUX_ANNUAL_SUM_NPP_PATCH
!    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
!    this%index_l2e_flux_annsum_npp_1d = index

!    id                            =
!    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
!    this%index_l2e_flux_agnpp_1d = index

!    id                            =
!    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
!    this%index_l2e_flux_bgnpp_1d = index

!    id                            =
!    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
!    this%index_l2e_flux_tempavg_agnpp_1d = index

!    id                            =
!    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
!    this%index_l2e_flux_tempavg_bgnpp_1d = index

!    id                            =
!    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
!    this%index_l2e_flux_annavg_agnpp_1d = index

!    id                            =
!    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
!    this%index_l2e_flux_annavg_bgnpp_1d = index

    id                            = L2E_FLUX_INFL
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_qflx_infl_1d = index

    id                            = L2E_FLUX_TOTDRAIN
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_qflx_totdrain_1d = index

    id                            = L2E_FLUX_GROSS_EVAP_SOIL
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_qflx_gross_evap_soil_1d = index

    id                            = L2E_FLUX_GROSS_INFL_SOIL
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_qflx_gross_infl_soil_1d = index

    id                            = L2E_FLUX_SURF
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_qflx_surf_1d = index

    id                            = L2E_FLUX_DEW_GRND
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_qflx_dew_ground_1d = index

    id                            = L2E_FLUX_DEW_SNOW
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_qflx_dew_snow_1d = index

    id                            = L2E_FLUX_SUB_SNOW
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_qflx_sub_snow_1d = index

    id                            = L2E_FLUX_SUB_SNOW_VOL
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_qflx_sub_snow_vol_1d = index

    id                            = L2E_FLUX_H2OSFC2TOPSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_qflx_h2osfc2topsoil_1d = index

    id                            = L2E_FLUX_SNOW2TOPSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_qflx_snow2topsoil_1d = index

    id                            = L2E_FLUX_ROOTSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_qflx_rootsoi_2d = index

    id                            = L2E_FLUX_SURF
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_qflx_runoff_1d = index

    id                            = L2E_FLUX_DRAIN_VR
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_qflx_drain_2d = index

    id                            = L2E_FLUX_TRAN_VEG
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_qflx_tran_veg_1d = index

    id                            = L2E_FLUX_ROOTSOI_FRAC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_qflx_rootsoi_vegfrac_2d = index

    id                            = L2E_STATE_QCHARGE
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_flux_qflx_bot_1d = index

    id                            = L2E_FLUX_INFL
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_qflx_infl_1d = index

    id                            = L2E_FLUX_ADV
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_qflx_adv_2d = index

    id                            = L2E_FLUX_TOTDRAIN
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_totdrain_1d = index

    id                            = L2E_FLUX_GROSS_EVAP_SOIL
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_gross_evap_soil_1d = index

    id                            = L2E_FLUX_GROSS_INFL_SOIL
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_gross_infl_soil_1d = index

    id                            = L2E_FLUX_DRAIN_VR
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_flux_drain_2d = index

    !carbon states
    id                            = L2E_STATE_C12_CARBON_POOLS_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_state_c12_decomp_pools_3d = index

    if(use_c13)then
      id                            = L2E_STATE_C13_CARBON_POOLS_VERTICALLY_RESOLVED
      call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
      this%index_l2e_bgc_state_c13_decomp_pools_3d = index

    endif

    if(use_c14)then
      id                            = L2E_STATE_C14_CARBON_POOLS_VERTICALLY_RESOLVED
      call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
      this%index_l2e_bgc_state_c14_decomp_pools_3d = index
     endif

    id                            = E2L_STATE_C12_CARBON_POOLS_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_state_c12_decomp_pools_3d = index

    if(use_c13)then
      id                            = E2L_STATE_C13_CARBON_POOLS_VERTICALLY_RESOLVED
      call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
      this%index_e2l_bgc_state_c13_decomp_pools_3d = index
    endif

    if(use_c14)then
      id                            = E2L_STATE_C14_CARBON_POOLS_VERTICALLY_RESOLVED
      call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
      this%index_e2l_bgc_state_c14_decomp_pools_3d = index
    endif

    id                            = E2L_STATE_C12_CARBON_CWD_VERTICALLY_INTEGRATED_1M
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_state_c12_cwd_1d = index

    id                            = E2L_STATE_C12_CARBON_LITR_VERTICALLY_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_state_c12_totlit_1d = index

    id                            = E2L_STATE_C12_CARBON_SOM_VERTICALLY_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_state_c12_totsom_1d = index

    id                            = E2L_STATE_C12_CARBON_LITR_VERTICALLY_INTEGRATED_1M
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_state_c12_totlit_1m_1d = index

    id                            = E2L_STATE_C12_CARBON_SOM_VERTICALLY_INTEGRATED_1M
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_state_c12_totsom_1m_1d = index

    id                            = E2L_STATE_C13_CARBON_CWD_VERTICALLY_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_state_c13_cwd_1d = index

    id                            = E2L_STATE_C13_CARBON_LITR_VERTICALLY_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_state_c13_totlit_1d = index

    id                            = E2L_STATE_C13_CARBON_SOM_VERTICALLY_INTEGRATED_1M
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_state_c13_totsom_1d = index

    id                            = E2L_STATE_C13_CARBON_LITR_VERTICALLY_INTEGRATED_1M
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_state_c13_totlit_1m_1d = index

    id                            = E2L_STATE_C13_CARBON_SOM_VERTICALLY_INTEGRATED_1M
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_state_c13_totsom_1m_1d = index

    id                            = E2L_STATE_C14_CARBON_CWD_VERTICALLY_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_state_c14_cwd_1d = index

    id                            = E2L_STATE_C14_CARBON_LITR_VERTICALLY_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_state_c14_totlit_1d = index

    id                            = E2L_STATE_C14_CARBON_SOM_VERTICALLY_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_state_c14_totsom_1d = index

    id                            = E2L_STATE_C14_CARBON_LITR_VERTICALLY_INTEGRATED_1M
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_state_c14_totlit_1m_1d = index

    id                            = E2L_STATE_C14_CARBON_SOM_VERTICALLY_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_state_c14_totsom_1m_1d = index

    !carbon fluxes
    id                            = L2E_FLUX_CARBON_C12_LITR_MET_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_flux_c12_litr_met_2d = index

    id                            = L2E_FLUX_CARBON_C12_LITR_CEL_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_flux_c12_litr_cel_2d = index

    id                            = L2E_FLUX_CARBON_C12_LITR_LIG_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_flux_c12_litr_lig_2d = index

    id                            = L2E_FLUX_CARBON_C12_LITR_CWD_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_flux_c12_litr_cwd_2d = index

    id                            = L2E_FLUX_CARBON_C13_LITR_MET_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_flux_c13_litr_met_2d = index

    id                            = L2E_FLUX_CARBON_C13_LITR_CEL_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_flux_c13_litr_cel_2d = index

    id                            = L2E_FLUX_CARBON_C13_LITR_LIG_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_flux_c13_litr_lig_2d = index

    id                            = L2E_FLUX_CARBON_C13_LITR_CWD_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_flux_c13_litr_cwd_2d = index

    id                            = L2E_FLUX_CARBON_C14_LITR_MET_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_flux_c14_litr_met_2d = index

    id                            = L2E_FLUX_CARBON_C14_LITR_CEL_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_flux_c14_litr_cel_2d = index

    id                            = L2E_FLUX_CARBON_C14_LITR_LIG_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_flux_c14_litr_lig_2d = index

    id                            = L2E_FLUX_CARBON_C14_LITR_CWD_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_flux_c14_litr_cwd_2d = index

    id                            = E2L_FLUX_CARBON_C12_HR_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_c12_hr_1d = index

    id                            = E2L_FLUX_CARBON_C13_HR_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_c13_hr_1d = index

    id                            = E2L_FLUX_CARBON_C14_HR_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_c14_hr_1d = index

    id                            = E2L_FLUX_CARBON_C12_FIRE_DECOMP_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_c12_fire_decomp_loss_1d = index

    id                            = E2L_FLUX_CARBON_C13_FIRE_DECOMP_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_c13_fire_decomp_loss_1d = index

    id                            = E2L_FLUX_CARBON_C14_FIRE_DECOMP_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_c14_fire_decomp_loss_1d = index

    id                            = E2L_FLUX_CARBON_C12_SOM_LEACHED_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_c12_som_leached_1d = index

    id                            = E2L_FLUX_CARBON_C12_SOM_RUNOFF_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_c12_som_runoff_1d = index

    id                            = E2L_FLUX_CARBON_C12_DIC_LEACHED_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_c12_dic_leached_1d = index

    id                            = E2L_FLUX_CARBON_C12_DIC_RUNOFF_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_c12_dic_runoff_1d = index

    !nitrogen states
    id                            = L2E_STATE_NITROGEN_POOLS_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_state_n14_decomp_pools_3d = index

    id                            = E2L_STATE_NITROGEN_POOLS_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_state_n14_decomp_pools_3d = index

    id                            = L2E_STATE_SOIL_NH4_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_state_n14_nh4_soil_2d = index

    id                            = L2E_STATE_SOIL_NO3_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_state_n14_no3_soil_2d = index

    id                            = E2L_STATE_SOIL_NH4_VERTICALLY_RESOLVED
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_state_n14_nh4_soil_2d = index

    id                            = E2L_STATE_SOIL_NO3_VERTICALLY_RESOLVED
    call e2l_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_state_n14_no3_soil_2d = index

    id                            = E2L_STATE_NITROGEN_CWD_VERTICALLY_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_state_n14_cwd_1d = index

    id                            = E2L_STATE_NITROGEN_LITR_VERTICALLY_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_state_n14_totlit_1d = index

    id                            = E2L_STATE_NITROGEN_SOM_VERTICALLY_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_state_n14_totsom_1d = index

    id                            = E2L_STATE_NITROGEN_LITR_VERTICALLY_INTEGRATED_1M
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_state_n14_totlit_1m_1d = index

    id                            = E2L_STATE_NITROGEN_SOM_VERTICALLY_INTEGRATED_1M
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_state_n14_totsom_1m_1d = index

    !nitrogen fluxes
    id                            =  L2E_FLUX_NITROGEN_LITR_MET_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_flux_n14_litr_met_2d = index

    id                            = L2E_FLUX_NITROGEN_LITR_CEL_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_flux_n14_litr_cel_2d = index

    id                            = L2E_FLUX_NITROGEN_LITR_LIG_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_flux_n14_litr_lig_2d = index

    id                            = L2E_FLUX_NITROGEN_LITR_CWD_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_flux_n14_litr_cwd_2d = index

    id                            = L2E_FLUX_NITROGEN_NH4_FIX_SOIL_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_flux_n14_nh4_fix_soil_2d = index

    id                            = L2E_FLUX_NITROGEN_NH4_ATMDEP_SOIL_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_flux_n14_nh4_atmdep_soil_2d = index

    id                            = L2E_FLUX_NITROGEN_NO3_ATMDEP_SOIL_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_flux_n14_no3_atmdep_soil_2d = index

    id                            = L2E_FLUX_NITROGEN_NH4_FERT_SOIL_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_flux_n14_nh4_fert_soil_2d = index

    id                            = L2E_FLUX_NITROGEN_NO3_FERT_SOIL_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_flux_n14_no3_fert_soil_2d = index

    id                            = E2L_FLUX_NITROGEN_FIRE_DECOMP_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_n14_fire_decomp_loss_1d = index

    id                            = E2L_FLUX_NITROGEN_SOM_LEACHED_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_n14_som_leached_1d = index

    id                            = E2L_FLUX_NITROGEN_SOM_RUNOFF_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_n14_som_runoff_1d = index

    id                            = E2L_FLUX_NITROGEN_DIN_LEACHED_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_n14_din_leached_1d = index

    id                            = E2L_FLUX_NITROGEN_DIN_RUNOFF_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_n14_din_runoff_1d = index

    id                            = E2L_FLUX_NITROGEN_NH4_LEACHED_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_n14_nh4_leached_1d = index

    id                            = E2L_FLUX_NITROGEN_NH4_RUNOFF_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_h14_nh4_runoff_1d = index

    id                            = E2L_FLUX_NITROGEN_NO3_LEACHED_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_n14_no3_leached_1d = index

    id                            = E2L_FLUX_NITROGEN_NO3_RUNOFF_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_h14_no3_runoff_1d = index

    id                            = E2L_FLUX_NITROGEN_NH4_SOI_EMI_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_n14_nh3_soil_emi_1d = index

    id                            = E2L_FLUX_NITROGEN_SOI_SUPP_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_n14_supplement_sminn_1d = index

    id                            = E2L_FLUX_NITROGEN_NITRIF_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_n14_nitrif_1d = index

    id                            = E2L_FLUX_NITROGEN_DENITRIF_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_n14_denitrif_1d = index

    id                            = E2L_FLUX_NITROGEN_N2O_NITRIF_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_n14_n2o_nitrif_1d = index

    id                            = E2L_FLUX_NITROGEN_N2O_DENITRIF_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_n14_n2o_denitrif_1d = index

    id                            = E2L_FLUX_NITROGEN_NH4_TO_PLANT_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_n14_smin_nh4_to_plant_veg_1d = index

    id                            = E2L_FLUX_NITROGEN_NO3_TO_PLANT_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_n14_smin_no3_to_plant_veg_1d = index

    id                            = E2L_FLUX_NITROGEN_SMINN_TO_PLANT_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_n14_sminn_to_plant_veg_1d = index

    !phosphorus states
    id                            = L2E_STATE_PHOSPHORUS_POOLS_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_state_p31_decomp_pools_3d = index

    id                            = E2L_STATE_PHOSPHORUS_POOLS_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_state_p31_decomp_pools_3d = index

    id                            = L2E_STATE_PHOSPHORUS_SOLUTION_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_state_p31_solutionp_soil_2d = index

    id                            = E2L_STATE_PHOSPHORUS_SOLUTION_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_state_p31_solutionp_soil_2d = index

    id                            = E2L_STATE_PHOSPHORUS_CWD_VERTICALLY_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_state_p31_cwd_1d = index

    id                            = E2L_STATE_PHOSPHORUS_LITR_VERTICALLY_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_state_p31_totlit_1d = index

    id                            = E2L_STATE_PHOSPHORUS_SOM_VERTICALLY_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_state_p31_totsom_1d = index

    id                            = E2L_STATE_PHOSPHORUS_LITR_VERTICALLY_INTEGRATED_1M
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_state_p31_totlit_1m_1d = index

    id                            = E2L_STATE_PHOSPHORUS_SOM_VERTICALLY_INTEGRATED_1M
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_state_p31_totsom_1m_1d = index

    id                            = E2L_FLUX_Phosphorus_SMINP_TO_PLANT_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_p31_sminp_to_plant_veg_1d = index

    id                            = E2L_FLUX_Phosphorus_SMINP_PLANT_DEMAND_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_p31_col_plant_demand_2d = index

    id                            = E2L_FLUX_Phosphorus_SMINP_ADSORB_TO_LABILE_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_p31_adsorb_to_labile_2d = index

    !phosphorus fluxes
    id                            = L2E_FLUX_Phosphorus_LITR_MET_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_flux_p31_litr_met_2d = index

    id                            = L2E_FLUX_Phosphorus_LITR_CEL_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_flux_p31_litr_cel_2d = index

    id                            = L2E_FLUX_Phosphorus_LITR_LIG_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_flux_p31_litr_lig_2d = index

    id                            = L2E_FLUX_Phosphorus_LITR_CWD_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_flux_p31_litr_cwd_2d = index

    id                            = L2E_FLUX_Phosphorus_SMINP_FERT_SOIL_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_flux_p31_sminp_fert_soil_2d = index

    id                            = L2E_FLUX_Phosphorus_SMINP_ATMDEP_SOIL_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_flux_p31_sminp_atmdep_soil_2d = index

    id                            = E2L_FLUX_Phosphorus_FIRE_DECOMP_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_p31_fire_decomp_loss_1d = index

    id                            = E2L_FLUX_Phosphorus_SOM_LEACHED_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_p31_som_leached_1d = index

    id                            = E2L_FLUX_Phosphorus_SOM_RUNOFF_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_p31_som_runoff_1d = index

    id                            = E2L_FLUX_Phosphorus_DIP_LEACHED_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_p31_dip_leached_1d = index

    id                            = E2L_FLUX_Phosphorus_DIP_RUNOFF_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_p31_dip_runoff_1d = index

    id                            = E2L_FLUX_Phosphorus_SOI_SUPP_INTEGRATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_e2l_bgc_flux_p31_supplement_sminp_1d = index

    !kinetic parameters for nutrient coupling
    id                            = L2E_parameter_plant_nh4_vmax_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_par_plant_nh4_vmax_2d = index

    id                            = L2E_parameter_plant_no3_vmax_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_par_plant_no3_vmax_2d = index

    id                            = L2E_parameter_plant_p_vmax_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_par_plant_p_vmax_2d = index

    id                            = L2E_parameter_plant_nh4_km_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_par_plant_nh4_km_2d = index

    id                            = L2E_parameter_plant_no3_km_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_par_plant_no3_km_2d = index

    id                            = L2E_parameter_plant_km_p_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_par_plant_p_km_2d = index

    id                            = L2E_parameter_plant_eff_ncompet_b_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_par_plant_eff_ncompetb_2d = index

    id                            = L2E_parameter_plant_eff_pcompet_b_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_par_plant_eff_pcompetb_2d = index

    id                            = L2E_parameter_plant_eff_frootc_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_par_plant_eff_frootc_2d = index

    id                            = L2E_parameter_minsurf_pcompet_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_par_minsurf_p_compet_2d = index

    id                            = L2E_parameter_nminsuf_nh4compet_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_par_minsurf_nh4_compet_2d = index

    id                            = L2E_parameter_km_minsurf_nh4_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_par_km_minsurf_nh4_2d = index

    id                            = L2E_parameter_km_minsurf_p_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_par_km_minsurf_p_2d = index

    id                            = L2E_parameter_decomp_eff_ncompet_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_par_decomp_eff_ncompetb_2d = index

    id                            = L2E_parameter_decomp_eff_pcompet_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_par_decomp_eff_pcompetb_2d = index

    id                            = L2E_parameter_nit_eff_ncompet_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_par_nitrif_eff_ncompetb_2d = index

    id                            = L2E_parameter_den_eff_ncompet_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_par_denitrif_eff_ncompetb_2d = index

    id                            = L2E_parameter_km_nit_nh4_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_par_km_nitrif_2d = index

    id                            = L2E_parameter_km_den_no3_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_par_km_denitrif_2d = index

    id                            = L2E_parameter_dsolutionp_dt_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_par_dsolutionp_dt_2d = index

    id                            = L2E_parameter_vmax_minsurf_p_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_par_vmax_minsurf_p_2d = index

    id                            = L2E_parameter_dlabp_dt_VERTICALLY_RESOLVED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    this%index_l2e_bgc_par_dlabp_dt_2d = index
    deallocate(em_stages)

  end subroutine EM_BETR_Populate_L2E_List

  !------------------------------------------------------------------------
  subroutine EM_BETR_Populate_E2L_List(this, e2l_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables returned by BeTR to ALM
    !
    ! !USES:
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_betr_type)                 :: this
    class(emi_data_list), intent(inout) :: e2l_list
    !
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages


  end subroutine EM_BETR_Populate_E2L_List



    !------------------------------------------------------------------------
  subroutine EM_BETR_BeginMassBalance_Solve(this, dt, nstep, bounds_clump)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use shr_kind_mod              , only : r8 => shr_kind_r8
    use abortutils                , only : endrun
    use shr_log_mod               , only : errMsg => shr_log_errMsg
    use clm_varctl                , only : iulog

    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_betr_type)                  :: this
    real(r8)             , intent(in)    :: dt
    integer              , intent(in)    :: nstep
    type(bounds_type)    , intent(in)    :: bounds_clump

    call ep_betr%SetClock(dtime = dt, nelapstep = nstep)

    call ep_betr%BeginMassBalanceCheck(bounds_clump)


  end subroutine EM_BETR_BeginMassBalance_Solve
    !------------------------------------------------------------------------
  subroutine EM_BETR_EndMassBalance_Solve(this, bounds_clump)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use shr_kind_mod              , only : r8 => shr_kind_r8
    use abortutils                , only : endrun
    use shr_log_mod               , only : errMsg => shr_log_errMsg
    use clm_varctl                , only : iulog

    !
    implicit none
    !
    ! !ARGUMENTS:
    class(em_betr_type)                  :: this
    type(bounds_type)    , intent(in)    :: bounds_clump

    call ep_betr%MassBalanceCheck(bounds_clump)

  end subroutine EM_BETR_EndMassBalance_Solve
    !------------------------------------------------------------------------
  subroutine EM_BETR_PreDiagSoilColWaterFlux_Solve(this, bounds, l2e_list, e2l_list)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use shr_kind_mod              , only : r8 => shr_kind_r8
    use abortutils                , only : endrun
    use shr_log_mod               , only : errMsg => shr_log_errMsg
    use clm_varctl                , only : iulog
    use clm_varpar                , only : nlevsoi
    use BeTR_decompMod            , only : betr_bounds_type

    !
    implicit none
    !
    ! !ARGUMENTS:
    type(bounds_type)    , intent(in)    :: bounds
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    !
    ! !LOCAL VARIABLES:
    real(r8), pointer      :: l2e_frac_h2osfc(:)
    real(r8), pointer      :: l2e_finundated(:)
    real(r8), pointer      :: l2e_h2osoi_liq(:,:)
    real(r8), pointer      :: l2e_h2osoi_ice(:,:)
    real(r8), pointer      :: l2e_h2osoi_liqvol(:,:)
    real(r8), pointer      :: l2e_h2osoi_icevol(:,:)
    real(r8), pointer      :: l2e_h2osoi_vol(:,:)
    real(r8), pointer      :: l2e_air_vol(:,:)
    real(r8), pointer      :: l2e_rho_vap(:,:)
    real(r8), pointer      :: l2e_rhvap_soi(:,:)
    real(r8), pointer      :: l2e_smp_l(:,:)
    integer, pointer       :: l2e_filter_nolakec(:)
    integer                :: l2e_num_nolakec
    integer                :: cc, c, fc, lbj, ubj

    type(betr_bounds_type) :: betr_bounds

    call l2e_list%GetPointerToReal1D(this%index_l2e_state_frac_h2osfc_1d  , l2e_frac_h2osfc    )
    call l2e_list%GetPointerToReal1D(this%index_l2e_state_finundated_1d   , l2e_finundated     )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_liq_2d   , l2e_h2osoi_liq     )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_ice_2d   , l2e_h2osoi_ice     )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_liqvol_2d, l2e_h2osoi_liqvol  )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_icevol_2d, l2e_h2osoi_icevol  )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_h2osoi_vol_2d   , l2e_h2osoi_vol     )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_air_vol_2d      , l2e_air_vol        )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_rho_vap_2d      , l2e_rho_vap        )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_rhvap_soi_2d    , l2e_rhvap_soi      )
    call l2e_list%GetPointerToReal2D(this%index_l2e_state_smp_l_2d        , l2e_smp_l          )

    call l2e_list%GetPointerToInt1D(this%index_l2e_filter_nolakec      , l2e_filter_nolakec )
    call l2e_list%GetIntValue(this%index_l2e_filter_num_nolakec        , l2e_num_nolakec    )

    cc  = 1
    lbj = 1
    ubj = nlevsoi

    do c = bounds%begc, bounds%endc
      if (.not. ep_betr%active_col(c)) cycle

      !assign waterstate
      ep_betr%biophys_forc(c)%finundated_col    (cc)         = l2e_finundated   (c)
      ep_betr%biophys_forc(c)%frac_h2osfc_col   (cc)         = l2e_frac_h2osfc  (c)
      ep_betr%biophys_forc(c)%h2osoi_liq_col    (cc,lbj:ubj) = l2e_h2osoi_liq   (c,lbj:ubj)
      ep_betr%biophys_forc(c)%h2osoi_ice_col    (cc,lbj:ubj) = l2e_h2osoi_ice   (c,lbj:ubj)
      ep_betr%biophys_forc(c)%h2osoi_liqvol_col (cc,lbj:ubj) = l2e_h2osoi_liqvol(c,lbj:ubj)
      ep_betr%biophys_forc(c)%h2osoi_icevol_col (cc,lbj:ubj) = l2e_h2osoi_icevol(c,lbj:ubj)
      ep_betr%biophys_forc(c)%h2osoi_vol_col    (cc,lbj:ubj) = l2e_h2osoi_vol   (c,lbj:ubj)
      ep_betr%biophys_forc(c)%air_vol_col       (cc,lbj:ubj) = l2e_air_vol      (c,lbj:ubj)
      ep_betr%biophys_forc(c)%rho_vap           (cc,lbj:ubj) = l2e_rho_vap      (c,lbj:ubj)
      ep_betr%biophys_forc(c)%rhvap_soi         (cc,lbj:ubj) = l2e_rhvap_soi    (c,lbj:ubj)
      ep_betr%biophys_forc(c)%smp_l_col         (cc,lbj:ubj) = l2e_smp_l        (c,lbj:ubj)

    enddo

    call ep_betr%PreDiagSoilColWaterFlux(l2e_num_nolakec , l2e_filter_nolakec)

  end subroutine EM_BETR_PreDiagSoilColWaterFlux_Solve

    !------------------------------------------------------------------------
  subroutine EM_BETR_StepWithoutDraiange(bounds, l2e_list, e2l_list)

    !
    implicit none
    !
    ! !ARGUMENTS:
    type(bounds_type)    , intent(in)    :: bounds
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list


  call ep_betr%StepWithoutDrainage(bounds,  col_pp, veg_pp)

  end subroutine EM_BETR_StepWithoutDraiange

    !------------------------------------------------------------------------
  subroutine EM_BETR_StepWithDraiange(bounds, l2e_list, e2l_list)
  implicit none

    ! !ARGUMENTS:
    class(em_betr_type)                  :: this
    type(bounds_type)    , intent(in)    :: bounds
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list


    call ep_betr%StepWithDrainage(bounds,  col_pp)
  end subroutine EM_BETR_StepWithDraiange

    !------------------------------------------------------------------------
  subroutine EM_BeTR_OutLoopSoilBGC(bounds, l2e_list, e2l_list )
  implicit none

    ! !ARGUMENTS:
    type(bounds_type)    , intent(in)    :: bounds
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list

    real(r8), pointer :: t_scalar(:,:)
    real(r8), pointer :: w_scalar(:,:)
    real(r8), pointer :: smin_no3_vr(:,:)
    real(r8), pointer :: smin_nh4_vr(:,:)
    real(r8), pointer :: solutionp_vr(:,:)

    integer :: c_l, j, fc


    call l2e_list%GetPointerToInt1D(this%index_l2e_filter_nolakec      , l2e_filter_nolakec )
    call l2e_list%GetIntValue(this%index_l2e_filter_num_nolakec        , l2e_num_nolakec    )
    call l2e_list%GetPointerToReal2D(this%index_l2e_bgc_coeff_t_scalar_2d       , t_scalar )
    call l2e_list%GetPointerToReal2D(this%index_l2e_bgc_coeff_w_scalar_2d       , w_scalar )
    call l2e_list%GetPointerToReal3D(this%index_l2e_bgc_state_c12_decomp_pools_3d , c12_decomp_cpools_vr_col )
    call l2e_list%GetPointerToReal3D(this%index_l2e_bgc_state_n14_decomp_pools_3d , decomp_npools_vr_col )
    call l2e_list%GetPointerToReal3D(this%index_l2e_bgc_state_p31_decomp_pools_3d , decomp_ppools_vr_col )
    call l2e_list%GetPointerToReal3D(this%index_l2e_bgc_coeff_decompk_3d      , decomp_k )

    call l2e_list%GetPointerToReal2D(this%index_l2e_bgc_par_plant_nh4_vmax_2d, plant_nh4_vmax_vr_patch )
    call l2e_list%GetPointerToReal2D(this%index_l2e_bgc_par_plant_no3_vmax_2d, plant_no3_vmax_vr_patch )
    call l2e_list%GetPointerToReal2D(this%index_l2e_bgc_par_plant_p_vmax_2d, plant_p_vmax_vr_patch )
    call l2e_list%GetPointerToReal2D(this%index_l2e_bgc_par_plant_nh4_km_2d, plant_nh4_km_vr_patch )
    call l2e_list%GetPointerToReal2D(this%index_l2e_bgc_par_plant_no3_km_2d, plant_no3_km_vr_patch )
    call l2e_list%GetPointerToReal2D(this%index_l2e_bgc_par_plant_p_km_2d, plant_p_km_vr_patch )
    call l2e_list%GetPointerToReal2D(this%index_l2e_bgc_par_plant_eff_ncompetb_2d, plant_eff_ncompet_b_vr_patch )
    call l2e_list%GetPointerToReal2D(this%index_l2e_bgc_par_plant_eff_pcompetb_2d, plant_eff_pcompet_b_vr_patch )
    call l2e_list%GetPointerToReal2D(this%index_l2e_bgc_par_plant_eff_frootc_2d, plant_eff_frootc_vr_patch )

    c_l=1
    call this%BeTRSetBounds(betr_bounds)
    do j = 1,nlevtrc_soil
      do fc = 1, num_soilc
        c = filter_soilc(fc)
        ep_betr%biophys_forc(c)%c12flx%in_t_scalar(c_l,j) = t_scalar(c,j)
        ep_betr%biophys_forc(c)%c12flx%in_w_scalar(c_l,j) = w_scalar(c,j)
        ep_betr%biophys_forc(c)%n14flx%in_sminn_no3_vr_col(c_l,j) = smin_no3_vr(c,j)
        ep_betr%biophys_forc(c)%n14flx%in_sminn_nh4_vr_col(c_l,j) = smin_nh4_vr(c,j)
        ep_betr%biophys_forc(c)%p31flx%in_sminp_vr_col(c_l,j) = solutionp_vr(c,j)
      enddo
    enddo

    do kk = 1, 7
      do j = 1,nlevtrc_soil
        do fc = 1, num_soilc
          c = filter_soilc(fc)
          ep_betr%biophys_forc(c)%c12flx%in_decomp_cpools_vr_col(c_l,j,kk)=c12_decomp_cpools_vr_col(c,j,kk)
          ep_betr%biophys_forc(c)%n14flx%in_decomp_npools_vr_col(c_l,j,kk)=decomp_npools_vr_col(c,j,kk)
          ep_betr%biophys_forc(c)%p31flx%in_decomp_ppools_vr_col(c_l,j,kk)=decomp_ppools_vr_col(c,j,kk)
          ep_betr%biogeo_flux(c)%c12flux_vars%decomp_k(c_l,j,kk)=decomp_k(c,j,kk)
        enddo
      enddo
    enddo

    !set autotrophic respiration

    !set kinetic parameters
    do fc = 1, num_soilc
      c = filter_soilc(fc)
      pp = 0
      val=1._r8

      do pi = 1, betr_maxpatch_pft
        if (pi <= col%npfts(c)) then
          p = col%pfti(c) + pi - 1
          if (pft%active(p) .and. (pft%itype(p) .ne. noveg)) then
            pp = pp + 1
            do j =1, betr_bounds%ubj
              ep_betr%betr(c)%plantNutkinetics%plant_nh4_vmax_vr_patch(pp,j) = plant_nh4_vmax_vr_patch(p,j)
              ep_betr%betr(c)%plantNutkinetics%plant_no3_vmax_vr_patch(pp,j) = plant_no3_vmax_vr_patch(p,j)
              ep_betr%betr(c)%plantNutkinetics%plant_p_vmax_vr_patch(pp,j) = plant_p_vmax_vr_patch(p,j)
              ep_betr%betr(c)%plantNutkinetics%plant_nh4_km_vr_patch(pp,j) = plant_nh4_km_vr_patch(p,j)/natomw
              ep_betr%betr(c)%plantNutkinetics%plant_no3_km_vr_patch(pp,j) = plant_no3_km_vr_patch(p,j)/natomw
              ep_betr%betr(c)%plantNutkinetics%plant_p_km_vr_patch(pp,j) = plant_p_km_vr_patch(p,j)/patomw
              ep_betr%betr(c)%plantNutkinetics%plant_eff_ncompet_b_vr_patch(pp,j)=plant_eff_ncompet_b_vr_patch(p,j)/natomw
              ep_betr%betr(c)%plantNutkinetics%plant_eff_pcompet_b_vr_patch(pp,j)=plant_eff_pcompet_b_vr_patch(p,j)/patomw
              ep_betr%betr(c)%plantNutkinetics%plant_eff_frootc_vr_patch(pp,j) = plant_eff_frootc_vr_patch(p,j)
            enddo
          endif
        endif
      enddo
      ep_betr%betr(c)%nactpft = pp
      do j = 1, betr_bounds%ubj
        ep_betr%betr(c)%plantNutkinetics%minsurf_p_compet_vr_col(c_l,j) = minsurf_p_compet_vr_col(c,j)/patomw
        ep_betr%betr(c)%plantNutkinetics%minsurf_nh4_compet_vr_col(c_l,j) = minsurf_nh4_compet_vr_col(c,j)/natomw
      enddo
    enddo

    !the following parameters are specific to ECACNP, and I assume they are
    !grid specific as they currently used in alm-cnp.
    if(index(reaction_method,'ecacnp')/=0 .or. index(reaction_method, 'ch4soil')/=0 &
       .or. index(reaction_method, 'v1eca')/=0)then
      do j =1, betr_bounds%ubj
        do fc = 1, num_soilc
          c = filter_soilc(fc)
          ep_betr%betr(c)%plantNutkinetics%km_minsurf_p_vr_col(c_l,j)  = km_minsurf_p_vr_col(c,j)/patomw
          ep_betr%betr(c)%plantNutkinetics%km_minsurf_nh4_vr_col(c_l,j)= km_minsurf_nh4_vr_col(c,j)/natomw
        enddo
      enddo
      if(lbcalib)then
        do j =1, betr_bounds%ubj
          do fc = 1, num_soilc
            c = filter_soilc(fc)
            ep_betr%betr(c)%plantNutkinetics%km_decomp_p_vr_col(c_l,j) = PlantMicKinetics_vars%km_decomp_p_vr_col(c,j)/patomw
            ep_betr%betr(c)%plantNutkinetics%km_decomp_nh4_vr_col(c_l,j)=PlantMicKinetics_vars%km_decomp_nh4_vr_col(c,j)/natomw
            ep_betr%betr(c)%plantNutkinetics%km_decomp_no3_vr_col(c_l,j)=PlantMicKinetics_vars%km_decomp_no3_vr_col(c,j)/natomw
            ep_betr%betr(c)%plantNutkinetics%km_nit_nh4_vr_col(c_l,j) = PlantMicKinetics_vars%km_nit_nh4_vr_col(c,j)/natomw
            ep_betr%betr(c)%plantNutkinetics%km_den_no3_vr_col(c_l,j) = PlantMicKinetics_vars%km_den_no3_vr_col(c,j)/natomw
          enddo
        enddo
      endif
    endif

    if(index(reaction_method,'v1eca')/=0)then
      do j =1, betr_bounds%ubj
        do fc = 1, num_soilc
          c = filter_soilc(fc)
          ep_betr%betr(c)%plantNutkinetics%decomp_eff_ncompet_b_vr_col(c_l,j)= decomp_eff_ncompet_b_vr_col(c,j)/natomw
          ep_betr%betr(c)%plantNutkinetics%decomp_eff_pcompet_b_vr_col(c_l,j)= decomp_eff_pcompet_b_vr_col(c,j)/patomw
          ep_betr%betr(c)%plantNutkinetics%nit_eff_ncompet_b_vr_col(c_l,j)   = nit_eff_ncompet_b_vr_col(c,j)/natomw
          ep_betr%betr(c)%plantNutkinetics%den_eff_ncompet_b_vr_col(c_l,j)   = den_eff_ncompet_b_vr_col(c,j)/natomw
          ep_betr%betr(c)%plantNutkinetics%km_nit_nh4_vr_col(c_l,j) = PlantMicKinetics_vars%km_nit_nh4_vr_col(c,j)/natomw
          ep_betr%betr(c)%plantNutkinetics%km_den_no3_vr_col(c_l,j) = PlantMicKinetics_vars%km_den_no3_vr_col(c,j)/natomw
          ep_betr%betr(c)%plantNutkinetics%dsolutionp_dt_vr_col(c_l,j)       = dsolutionp_dt_vr_col(c,j)/patomw   ! g/m2/s
          ep_betr%betr(c)%plantNutkinetics%vmax_minsurf_p_vr_col(c_l,j)      = vmax_minsurf_p_vr_col(c,j)/patomw   ! g/m3
          ep_betr%betr(c)%plantNutkinetics%dlabp_dt_vr_col(c_l,j)            = dlabp_dt_vr_col(c,j)/patomw
        enddo
      enddo
    endif

    !execute betr calculations
    call ep_betr%OutLoopSoilBGC(bounds,  col_pp, veg_pp)

    !output variables after executation of betr

  c_l=1
  do kk = 1, 7
    do j = 1,nlevtrc_soil
      do c = bounds%begc, bounds%endc
        if(.not. this%active_col(c))cycle

        decomp_k(c,j,kk) = this%biogeo_flux(c)%c12flux_vars%decomp_k(c_l,j,kk)
        c12_decomp_cpools_vr_col(c,j,kk)=this%biogeo_state(c)%c12state_vars%decomp_cpools_vr(c_l,j,kk)
        decomp_npools_vr_col(c,j,kk) = this%biogeo_state(c)%n14state_vars%decomp_npools_vr(c_l,j,kk)
        decomp_ppools_vr_col(c,j,kk) = this%biogeo_state(c)%p31state_vars%decomp_ppools_vr(c_l,j,kk)
      enddo
    enddo
  enddo
  !extract plant nutrient uptake fluxes, soil respiration, denitrification, nitrification
  !
  do c = bounds%begc, bounds%endc
    if(.not. this%active_col(c))cycle
    pi = 0
    do p = col%pfti(c), col%pftf(c)
      if (pft%active(p) .and. (pft%itype(p) .ne. noveg)) then
        pi = pi + 1
        nitrogenflux_vars%smin_nh4_to_plant_patch(p) = this%biogeo_flux(c)%n14flux_vars%smin_nh4_to_plant_patch(pi)
        nitrogenflux_vars%smin_no3_to_plant_patch(p) = this%biogeo_flux(c)%n14flux_vars%smin_no3_to_plant_patch(pi)
        nitrogenflux_vars%sminn_to_plant_patch(p) = nitrogenflux_vars%smin_nh4_to_plant_patch(p) + nitrogenflux_vars%smin_no3_to_plant_patch(p)
        phosphorusflux_vars%sminp_to_plant_patch(p)  = this%biogeo_flux(c)%p31flux_vars%sminp_to_plant_patch(pi)
      else
        nitrogenflux_vars%smin_nh4_to_plant_patch(p) = 0._r8
        nitrogenflux_vars%smin_no3_to_plant_patch(p) = 0._r8
        phosphorusflux_vars%sminp_to_plant_patch(p) = 0._r8
        nitrogenflux_vars%sminn_to_plant_patch(p) = 0._r8
      endif
    enddo
    nitrogenflux_vars%actual_immob_col(c)=0._r8
!    print*,'immob',nitrogenflux_vars%actual_immob_col(c)
  enddo

  do j = 1,nlevtrc_soil
    do c = bounds%begc, bounds%endc
      if(.not. this%active_col(c))cycle
      nitrogenflux_vars%col_plant_pdemand_vr(c,j)  = this%biogeo_flux(c)%p31flux_vars%col_plant_pdemand_vr(c_l,j)
      nitrogenflux_vars%f_denit_vr_col(c,j)        = this%biogeo_flux(c)%n14flux_vars%f_denit_vr_col(c_l,j)
      nitrogenflux_vars%f_n2o_denit_vr_col(c,j)    = this%biogeo_flux(c)%n14flux_vars%f_n2o_denit_vr_col(c_l,j)
      nitrogenflux_vars%f_nit_vr_col(c,j)          = this%biogeo_flux(c)%n14flux_vars%f_nit_vr_col(c_l,j)
      nitrogenflux_vars%f_n2o_nit_vr_col(c,j)      = this%biogeo_flux(c)%n14flux_vars%f_n2o_nit_vr_col(c_l,j)
      phosphorusflux_vars%adsorb_to_labilep_vr(c,j)= this%biogeo_flux(c)%p31flux_vars%adsorb_to_labilep_vr_col(c_l,j)
      c12_cflx_vars%hr_vr_col(c,j)                 = this%biogeo_flux(c)%c12flux_vars%hr_vr_col(c_l,j)
      c12_cflx_vars%phr_vr_col(c,j)                = this%biogeo_flux(c)%c12flux_vars%phr_vr_col(c_l,j)
      c12_cflx_vars%o_scalar_col(c,j)              = this%biogeo_flux(c)%c12flux_vars%o_scalar_col(c_l,j)
      nitrogenstate_vars%smin_nh4_vr_col(c,j)      = this%biogeo_state(c)%n14state_vars%sminn_nh4_vr_col(c_l,j)
      nitrogenstate_vars%smin_no3_vr_col(c,j)      = this%biogeo_state(c)%n14state_vars%sminn_no3_vr_col(c_l,j)
      nitrogenflux_vars%supplement_to_sminn_vr_col(c,j) = this%biogeo_flux(c)%n14flux_vars%supplement_to_sminn_vr_col(c_l,j)
      nitrogenflux_vars%smin_nh4_to_plant_vr_col(c,j) = this%biogeo_flux(c)%n14flux_vars%smin_nh4_to_plant_vr_col(c_l,j)
      nitrogenflux_vars%smin_no3_to_plant_vr_col(c,j) = this%biogeo_flux(c)%n14flux_vars%smin_no3_to_plant_vr_col(c_l,j)
      phosphorusflux_vars%sminp_to_plant_vr_col(c,j) = this%biogeo_flux(c)%p31flux_vars%sminp_to_plant_vr_col(c_l,j)
      phosphorusflux_vars%supplement_to_sminp_vr_col(c,j) =this%biogeo_flux(c)%p31flux_vars%supplement_to_sminp_vr_col(c_l,j)
      phosphorusflux_vars%net_mineralization_p_vr_col(c,j) = this%biogeo_flux(c)%p31flux_vars%net_mineralization_p_vr_col(c_l,j)
      nitrogenflux_vars%actual_immob_col(c)= nitrogenflux_vars%actual_immob_col(c) + col%dz(c,j)*&
         (this%biogeo_flux(c)%n14flux_vars%smin_nh4_immob_vr_col(c_l,j) + this%biogeo_flux(c)%n14flux_vars%smin_no3_immob_vr_col(c_l,j))
      c12_cflx_vars%somhr_col(c) = c12_cflx_vars%somhr_col(c)  + col%dz(c,j)*&
         this%biogeo_flux(c)%c12flux_vars%somhr_vr_col(c_l,j)
      c12_cflx_vars%lithr_col(c) = c12_cflx_vars%lithr_col(c)  + col%dz(c,j)*&
         this%biogeo_flux(c)%c12flux_vars%lithr_vr_col(c_l,j)
    enddo
  enddo

  end subroutine EM_BeTR_OutLoopSoilBGC
end module ExternalModelBETRMod
