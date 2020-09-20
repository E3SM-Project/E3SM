module clm_interface_funcsMod
!=================================================================================================
! CLM Theraml-Hydrology (TH) & BioGeoChemistry (BGC) Interface: Modules
! created: 8/25/2015
! updated: June-2017
!=================================================================================================

#include "shr_assert.h"


  ! MODULE: clm_interface_funcsMod
  !--------------------------------------------------------------------------------------
  ! DESCRIPTION:
  ! Coupling of CLM with any specific Soil BGC module Consists of 3 STEPS:
  ! STEP-1:   clm vars             -> clm_interface_data (i.e. elm_interface_dataType)  ; pass clm vars to clm_interface_data
  ! STEP-2:   clm_interface_data   -> soil bgc module -> clm_interface_data
  !      2.1: clm_interface_data   -> soil bgc module
  !      2.2: run soil bgc module
  !      2.3: soil bgc module      -> clm_interface_data
  ! STEP-3:   clm_interface_data   -> clm vars
  !--------------------------------------------------------------------------------------


  !--------------------------------------------------------------------------------------
  ! !USES:
  ! clm g/l/c/p constants
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use GridcellType          , only : grc_pp
  use LandunitType          , only : lun_pp
  use ColumnType            , only : col_pp
  use ColumnDataType        , only : col_es, col_ef, col_ws, col_wf 
  use ColumnDataType        , only : col_cs, col_cf
  use ColumnDataType        , only : col_ns, col_nf
  use ColumnDataType        , only : col_ps, col_pf
  use VegetationType        , only : veg_pp

  use decompMod             , only : bounds_type

  ! (dummy) variable definitions
  use atm2lndType           , only : atm2lnd_type
  use SoilStateType         , only : soilstate_type
  use WaterStateType        , only : waterstate_type
  use WaterFluxType         , only : waterflux_type
  use SoilHydrologyType     , only : soilhydrology_type
  use TemperatureType       , only : temperature_type
  use EnergyFluxType        , only : energyflux_type

  use CNStateType           , only : cnstate_type
  use CNCarbonFluxType      , only : carbonflux_type
  use CNCarbonStateType     , only : carbonstate_type
  use CNNitrogenFluxType    , only : nitrogenflux_type
  use CNNitrogenStateType   , only : nitrogenstate_type

  use CH4Mod                , only : ch4_type

  use PhotosynthesisType    , only : photosyns_type
  use cropType              , only : crop_type
  use CanopyStateType       , only : canopystate_type
  use PhosphorusStateType   , only : phosphorusstate_type
  use PhosphorusFluxType    , only : phosphorusflux_type

  use SoilWaterRetentionCurveMod    , only : soil_water_retention_curve_type

  use elm_interface_dataType        , only : elm_interface_data_type
  use elm_interface_thType          , only : elm_interface_th_datatype
  use elm_interface_bgcType         , only : elm_interface_bgc_datatype

  ! most used constants in this module
  use clm_varpar            , only : nlevsoi, nlevsno, nlevgrnd, nlevdecomp_full
  use clm_varpar            , only : ndecomp_pools, ndecomp_cascade_transitions
  use clm_varpar            , only : max_patch_per_col
  use clm_varcon            , only : denh2o, denice, tfrz, dzsoi_decomp
  use landunit_varcon       , only : istsoil, istcrop

  ! misc.
  use abortutils            , only : endrun
  use clm_varctl            , only : nu_com
  !--------------------------------------------------------------------------------------

  implicit none

  save

  private    ! By default everything is private

  ! LOCAL VARIABLES:


  !--------------------------------------------------------------------------------------
  ! (1) GENERIC SUBROUTINES: used by any specific soil BGC module
  ! pass clm variables to clm_bgc_data
  public    :: get_clm_data                 ! STEP-1: clm vars -> clm_interface_data

  ! pass clm variables to clm_interface_data, called by get_clm_data
  private   :: get_clm_soil_property        ! STEP-1.1: soil properties
  private   :: get_clm_soil_th_state        ! STEP-1.2: thermohydrology (TH) state vars
  private   :: get_clm_bgc_state            ! STEP-1.3: state vars
  private   :: get_clm_bgc_flux             ! STEP-1.4: flux vars

  ! STEP-3.x: clm_interface_data -> clm vars
  ! update clm variables from clm_interface_data,
  ! e.g., called in 'update_bgc_data_clm2clm' and 'update_bgc_data_pf2clm'
  ! specific bgc-module (e.g., PFLOTRAN) requires certain combination of these subroutines
  private   :: update_bgc_state_decomp
  private   :: update_bgc_state_smin
  private   :: update_bgc_flux_decomp_sourcesink
  private   :: update_bgc_flux_decomp_cascade
  private   :: update_bgc_flux_smin
  private   :: update_bgc_flux_nitdenit
  private   :: update_bgc_flux_gas_pf
  private   :: update_soil_moisture
  private   :: update_soil_temperature

  !--------------------------------------------------------------------------------------
  ! (2) SPECIFIC SUBROUTINES: used by a specific soil BGC module
  ! (2.1) Specific Subroutines for running clm-bgc (CN or BGC) through interface
  ! if (use_clm_interface .and. use_clm_bgc)
  public    :: clm_bgc_run              ! STEP-2:   clm_interface_data  -> clm-bgc module -> clm_interface_data    ; called in clm_driver
  private   :: clm_bgc_get_data         ! STEP-2.1: clm_interface_data  -> clm-bgc module                          ; called in clm_bgc_run
                                        ! STEP-2.2: run clm-bgc module                                             ; see SoilLittDecompAlloc in SoilLittDecompMod
  private   :: clm_bgc_update_data      ! STEP-2.3: clm-bgc module-> clm_interface_data                            ; called in clm_bgc_run
  public    :: update_bgc_data_clm2clm  ! STEP-3:   clm_interface_data  -> clm vars                                ; called in clm_driver

  ! (2.2) Specific Subroutines for CLM-PFLOTRAN Coupling: update clm variables from pflotran
  ! if (use_clm_interface .and. use_pflotran)
  public    :: update_bgc_data_pf2clm   ! STEP-3:   clm_interface_data  -> clm vars                                ; called in clm_driver
                                        ! STEP-2:   see 'clm_pf_run' in clm_interface_pflotranMod

  public    :: update_th_data_pf2clm
  !--------------------------------------------------------------------------------------

contains


!--------------------------------------------------------------------------------------
  subroutine get_clm_data(clm_idata,                              &
           bounds, num_soilc, filter_soilc,                       &
           num_soilp, filter_soilp,                               &
           atm2lnd_vars, soilstate_vars,                          &
           waterstate_vars, waterflux_vars,                       &
           temperature_vars, energyflux_vars,                     &
           cnstate_vars, carbonflux_vars, carbonstate_vars,       &
           nitrogenflux_vars, nitrogenstate_vars,                 &
           phosphorusflux_vars, phosphorusstate_vars,             &
           ch4_vars                                               &
           )

    implicit none

    ! !ARGUMENTS:
    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                     , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                     , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                     , intent(in)    :: filter_soilp(:)   ! filter for soil patches
    type(atm2lnd_type)          , intent(in)    :: atm2lnd_vars
    type(soilstate_type)        , intent(in)    :: soilstate_vars

    type(waterstate_type)       , intent(in)    :: waterstate_vars
    type(waterflux_type)        , intent(in)    :: waterflux_vars
    type(temperature_type)      , intent(in)    :: temperature_vars
    type(energyflux_type)       , intent(in)    :: energyflux_vars

    type(cnstate_type)          , intent(in)    :: cnstate_vars
    type(carbonflux_type)       , intent(in)    :: carbonflux_vars
    type(carbonstate_type)      , intent(in)    :: carbonstate_vars
    type(nitrogenflux_type)     , intent(in)    :: nitrogenflux_vars
    type(nitrogenstate_type)    , intent(in)    :: nitrogenstate_vars
    type(phosphorusflux_type)   , intent(in)    :: phosphorusflux_vars
    type(phosphorusstate_type)  , intent(in)    :: phosphorusstate_vars
    type(ch4_type)              , intent(in)    :: ch4_vars

    type(elm_interface_data_type), intent(inout) :: clm_idata

    ! LOCAL
    !type(elm_interface_th_datatype) , pointer :: clm_idata_th
    !type(elm_interface_bgc_datatype), pointer :: clm_idata_bgc


    character(len=256) :: subname = "get_clm_data"
    !-----------------------------------------------------------------------

     associate ( &
      clm_idata_th  => clm_idata%th,  &
      clm_idata_bgc => clm_idata%bgc  &
     )

    call get_clm_soil_property(clm_idata,                   &
                    bounds, num_soilc, filter_soilc,        &
                    soilstate_vars, cnstate_vars)

    call get_clm_soil_th_state(clm_idata_th,                &
                   bounds, num_soilc, filter_soilc,         &
                   atm2lnd_vars, soilstate_vars,            &
                   waterstate_vars, temperature_vars)

    call get_clm_soil_th_flux(clm_idata_th,                 &
                       bounds, num_soilc, filter_soilc,     &
                       waterflux_vars, energyflux_vars)

    call get_clm_bgc_state(clm_idata_bgc,                   &
                    bounds, num_soilc, filter_soilc,        &
                    atm2lnd_vars, soilstate_vars,           &
                    carbonstate_vars, nitrogenstate_vars,   &
                    phosphorusstate_vars,                   &
                    ch4_vars)

    call get_clm_bgc_flux(clm_idata_bgc,                    &
                    bounds, num_soilc, filter_soilc,        &
                    cnstate_vars, carbonflux_vars,          &
                    nitrogenflux_vars, phosphorusflux_vars, &
                    ch4_vars)

    end associate
  end subroutine get_clm_data
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
  subroutine get_clm_soil_property(clm_idata,               &
                        bounds, num_soilc, filter_soilc,    &
                        soilstate_vars, cnstate_vars)

    !
    ! !DESCRIPTION:
    ! get soil column physical and biogeochemical properties
    !
    ! !USES:
    use CNDecompCascadeConType, only : decomp_cascade_con
    !
    ! !ARGUMENTS:

    implicit none

    ! LOCAL VARAIBLES:

    type(bounds_type)        , intent(in) :: bounds           ! bounds
    integer                  , intent(in) :: num_soilc        ! number of column soil points in column filter
    integer                  , intent(in) :: filter_soilc(:)  ! column filter for soil points
    type(soilstate_type)     , intent(in) :: soilstate_vars
    type(cnstate_type)       , intent(in) :: cnstate_vars

    type(elm_interface_data_type), intent(inout) :: clm_idata

    integer  :: fc, g, l, c, j, k      ! indices
    integer  :: gcount, cellcount

    character(len= 32) :: subname = 'get_clm_soil_property' ! subroutine name

    associate ( &
         ! Assign local pointer to derived subtypes components (column-level)
         z                  => col_pp%z                                                , & !  [real(r8) (:,:)]  layer depth (m)
         dz                 => col_pp%dz                                               , & !  [real(r8) (:,:)]  layer thickness depth (m)
         zi                 => col_pp%zi                                               , & !

         bd                 => soilstate_vars%bd_col                                , & !
         bsw                => soilstate_vars%bsw_col                               , & !  [real(r8) (:,:)]  Clapp and Hornberger "b" (nlevgrnd)
         hksat              => soilstate_vars%hksat_col                             , & !  [real(r8) (:,:)]  hydraulic conductivity at saturation (mm H2O /s) (nlevgrnd)
         sucsat             => soilstate_vars%sucsat_col                            , & !  [real(r8) (:,:)]  minimum soil suction (mm) (nlevgrnd)
         watsat             => soilstate_vars%watsat_col                            , & !  [real(r8) (:,:)]  volumetric soil water at saturation (porosity) (nlevgrnd)
         watfc              => soilstate_vars%watfc_col                             , & !  [real(r8) (:,:)]  volumetric soil water at saturation (porosity) (nlevgrnd)
         watmin             => soilstate_vars%watmin_col                            , & !   col minimum volumetric soil water (nlevsoi)
         sucmin             => soilstate_vars%sucmin_col                            , & !   col minimum allowable soil liquid suction pressure (mm) [Note: sucmin_col is a negative value, while sucsat_col is a positive quantity]
         !
         cellorg            => soilstate_vars%cellorg_col                           , & !  Input:  [real(r8) (:,:)  ]  column 3D org (kg/m3 organic matter) (nlevgrnd)
         !
!          porosity           => soilstate_vars%porosity_col                          , &
         eff_porosity       => soilstate_vars%eff_porosity_col                      , &
         !
         rootfr             => soilstate_vars%rootfr_col                            , & ! pft-level effective fraction of roots in each soil layer
         !
         initial_cn_ratio   => decomp_cascade_con%initial_cn_ratio                  , &
         initial_cp_ratio   => decomp_cascade_con%initial_cp_ratio                  , &
         !
         decomp_pool_name   => decomp_cascade_con%decomp_pool_name_history          , &
         floating_cn_ratio  => decomp_cascade_con%floating_cn_ratio_decomp_pools    , &
         floating_cp_ratio  => decomp_cascade_con%floating_cp_ratio_decomp_pools    , &
         decomp_k_pools     => decomp_cascade_con%decomp_k_pools                    , &
         adfactor_kd_pools  => decomp_cascade_con%spinup_factor                     , &
         rf_decomp_cascade       => cnstate_vars%rf_decomp_cascade_col              , &
         pathfrac_decomp_cascade => cnstate_vars%pathfrac_decomp_cascade_col          &

         )

!-------------------------------------------------------------------------------------
    ! constants:
    clm_idata%bgc%ndecomp_pools          = ndecomp_pools
    clm_idata%bgc%decomp_pool_name(:)    = decomp_pool_name(:)
    clm_idata%bgc%floating_cn_ratio(:)   = floating_cn_ratio(:)
    clm_idata%bgc%floating_cp_ratio(:)   = floating_cp_ratio(:)

    clm_idata%bgc%initial_cn_ratio(:)    = initial_cn_ratio(:)
    clm_idata%bgc%initial_cp_ratio(:)    = initial_cp_ratio(:)
    clm_idata%bgc%decomp_k_pools(:)      = decomp_k_pools(1:ndecomp_pools)
    clm_idata%bgc%adfactor_kd_pools(:)   = adfactor_kd_pools(1:ndecomp_pools)

    do fc = 1, num_soilc
        c = filter_soilc(fc)

        clm_idata%z(c,:)                 = z(c,:)
        clm_idata%zi(c,:)                = zi(c,:)
        clm_idata%dz(c,:)                = dz(c,:)
        clm_idata%bd_col(c,:)            = bd(c,:)
        clm_idata%bsw_col(c,:)           = bsw(c,:)
        clm_idata%hksat_col(c,:)         = hksat(c,:)
        clm_idata%sucsat_col(c,:)        = sucsat(c,:)
        clm_idata%watsat_col(c,:)        = watsat(c,:)
        clm_idata%watfc_col(c,:)         = watfc(c,:)
        clm_idata%watmin_col(c,:)        = watmin(c,:)
        clm_idata%sucmin_col(c,:)        = sucmin(c,:)

        clm_idata%porosity_col(c,:)      = watsat(c,:)
        clm_idata%eff_porosity_col(c,:)  = eff_porosity(c,:)

        clm_idata%cellorg_col(c,:)       = cellorg(c,:)

        clm_idata%rootfr_col(c,:)        = rootfr(c,:)

        !
        do k = 1, ndecomp_cascade_transitions
            clm_idata%bgc%rf_decomp_cascade_col(c,:,k)           = rf_decomp_cascade(c,:,k)
            clm_idata%bgc%pathfrac_decomp_cascade_col(c,:,k)     = pathfrac_decomp_cascade(c,:,k)
        end do

    end do

  end associate
  end subroutine get_clm_soil_property
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
  subroutine get_clm_soil_th_state(clm_idata_th,            &
                       bounds, num_soilc, filter_soilc,     &
                       atm2lnd_vars, soilstate_vars,        &
                       waterstate_vars, temperature_vars)
  !
  ! !DESCRIPTION:
  !  get soil temperature/saturation from CLM to soil BGC module
  !
  ! !USES:
    use clm_time_manager    , only : get_nstep
    use shr_const_mod       , only : SHR_CONST_G


  ! !ARGUMENTS:
    implicit none

    type(bounds_type)        , intent(in) :: bounds           ! bounds
    integer                  , intent(in) :: num_soilc        ! number of column soil points in column filter
    integer                  , intent(in) :: filter_soilc(:)  ! column filter for soil points
    type(atm2lnd_type)       , intent(in) :: atm2lnd_vars
    type(soilstate_type)     , intent(in) :: soilstate_vars
    type(waterstate_type)    , intent(in) :: waterstate_vars
    type(temperature_type)   , intent(in) :: temperature_vars

    type(elm_interface_th_datatype)       , intent(inout) :: clm_idata_th

  ! !LOCAL VARIABLES:
    integer  :: fc, c, j         ! indices

  !EOP
  !-----------------------------------------------------------------------
    associate ( &
      soilpsi               => soilstate_vars%soilpsi_col               , & !
      !
      frac_sno_eff          => col_ws%frac_sno_eff         , & !
      frac_h2osfc           => col_ws%frac_h2osfc          , & !
      h2osoi_vol            => col_ws%h2osoi_vol           , & ! [real(r8) (:,:)] volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevgrnd)
      h2osoi_liq            => col_ws%h2osoi_liq           , & ! [real(r8) (:,:)] liquid water (kg/m2) (-nlevsno+1:nlevgrnd)
      h2osoi_ice            => col_ws%h2osoi_ice           , & ! [real(r8) (:,:)] ice lens (kg/m2) (-nlevsno+1:nlevgrnd)
      !
      t_soisno              => col_es%t_soisno            , & ! [real(r8) (:,:)] snow-soil temperature (Kelvin) (-nlevsno+1:nlevgrnd)
      t_grnd                => col_es%t_grnd              , & ! [real(r8) (:)] ground (snow/soil1/surfwater-mixed) temperature (Kelvin)
      t_h2osfc              => col_es%t_h2osfc            , & ! [real(r8) (:)] surface water temperature (Kelvin)
      t_nearsurf            => col_es%t_nearsurf          , & ! [real(r8) (:)] near surface air temperature (Kelvin)
      !
      forc_pbot             => atm2lnd_vars%forc_pbot_not_downscaled_grc  & ! atmospheric pressure (Pa)
      )


    !--------------------------------------------------------------------------------------
    ! grid:
    clm_idata_th%forc_pbot_grc    = forc_pbot

    do fc = 1,num_soilc
        c = filter_soilc(fc)

        clm_idata_th%frac_sno_eff_col(c)         = frac_sno_eff(c)
        clm_idata_th%frac_h2osfc_col(c)          = frac_h2osfc(c)

        clm_idata_th%t_grnd_col(c)               = t_grnd(c)
        clm_idata_th%t_h2osfc_col(c)             = t_h2osfc(c)
        clm_idata_th%t_nearsurf_col(c)           = t_nearsurf(c)

        do j = -nlevsno+1,nlevgrnd
            if(j>=1) then
                clm_idata_th%soilpsi_col(c,j)        = soilpsi(c,j)
                clm_idata_th%h2osoi_vol_col(c,j)     = h2osoi_vol(c,j)
            endif

            clm_idata_th%h2osoi_liq_col(c,j)         = h2osoi_liq(c,j)
            clm_idata_th%h2osoi_ice_col(c,j)         = h2osoi_ice(c,j)
            clm_idata_th%t_soisno_col(c,j)           = t_soisno(c,j)
        end do

    end do

    end associate
  end subroutine get_clm_soil_th_state
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
  subroutine get_clm_soil_th_flux(clm_idata_th,             &
                       bounds, num_soilc, filter_soilc,     &
                       waterflux_vars, energyflux_vars)
  !
  ! !DESCRIPTION:
  !  get soil temperature/saturation from CLM to soil BGC module
  !
  ! !USES:
    use clm_time_manager    , only : get_nstep
    use shr_const_mod       , only : SHR_CONST_G


  ! !ARGUMENTS:
    implicit none

    type(bounds_type)        , intent(in) :: bounds           ! bounds
    integer                  , intent(in) :: num_soilc        ! number of column soil points in column filter
    integer                  , intent(in) :: filter_soilc(:)  ! column filter for soil points
    type(waterflux_type)     , intent(in) :: waterflux_vars
    type(energyflux_type)    , intent(in) :: energyflux_vars

    type(elm_interface_th_datatype)       , intent(inout) :: clm_idata_th

  ! !LOCAL VARIABLES:
    integer  :: fc, c, j         ! indices

  !EOP
  !-----------------------------------------------------------------------
    associate ( &
      qflx_top_soil     => col_wf%qflx_top_soil         , & ! [real(:,:)] net liq. water input into top of soil column (mmH2O/s)
      qflx_evap_soil    => col_wf%qflx_ev_soil          , & ! [real(:)] ! col soil surface evaporation (mm H2O/s) (+ = to atm)
      qflx_evap_h2osfc  => col_wf%qflx_ev_h2osfc        , & ! [real(:)] ! col water surface evaporation (mm H2O/s) (+ = to atm)
      qflx_evap_snow    => col_wf%qflx_ev_snow          , & ! [real(:)] ! col snow surface evaporation (mm H2O/s) (+ = to atm)
      qflx_subl_snow    => col_wf%qflx_sub_snow         , & ! [real(:)] ! col snow sublimation (mm H2O/s) (+ = to atm)
      qflx_tran_veg     => col_wf%qflx_tran_veg         , & ! [real(:)] ! col plant transpiration (mm H2O/s) (+ = to atm)
      qflx_rootsoil     => col_wf%qflx_rootsoi          , & ! [real(:,:)] ! col vertically-resolved root and soil water exchange [mm H2O/s] [+ into root]
      !
      htvp              => col_ef%htvp                 , & ! [real(:) ! latent heat of vapor of water (or sublimation) [j/kg]
      eflx_bot          => col_ef%eflx_bot             , & ! [real(:) ! col heat flux from beneath the soil or ice column (W/m**2)
      eflx_soil_grnd    => col_ef%eflx_soil_grnd       , & ! [real(:) ! col soil (ground) heat flux (W/m**2) [+ = into ground]
      eflx_fgr0_snow    => col_ef%eflx_fgr0_snow       , & ! [real(:) ! col ground heat flux from snow bottom to first soil layer (W/m**2) [+ = into soil]
      eflx_fgr0_h2osfc  => col_ef%eflx_fgr0_h2osfc     , & ! [real(:) ! col ground heat flux from surface water bottom to first soil layer (W/m**2) [+ = into soil]
      eflx_fgr0_soil    => col_ef%eflx_fgr0_soil       , & ! [real(:) ! col ground heat flux from near-surface air to first soil layer (W/m**2) [+ = into soil]
      eflx_rnet_soil    => col_ef%eflx_rnet_soil         & ! [real(:) ! net radiation flux between soil layer 1 and above-air, excluding SH and LE (i.e. radiation form only ) (W/m2) [+ = into soil]

    )

    ! a few notes:
    !   - 'qflx_evap_soil' appears for total soil surface, esp. bare soil; 'qflx_ev_soil/snow/h2osfc' are actually applied for in soil water modules
    !   - 'qflx_ev_snow' vs. 'qflx_sub_snow': the former is for total evap from both solid/liq., the latter is from solid snow pack (normally shall be same)
    !                        there is another variable 'qlfx_evap_grnd', which are those from liq. water when snow
    !--------------------------------------------------------------------------------------
!
    do fc = 1,num_soilc
        c = filter_soilc(fc)

        clm_idata_th%qflx_top_soil_col(c)        = qflx_top_soil(c)
        clm_idata_th%qflx_evap_soil_col(c)       = qflx_evap_soil(c)
        clm_idata_th%qflx_evap_h2osfc_col(c)     = qflx_evap_h2osfc(c)
        clm_idata_th%qflx_evap_snow_col(c)       = qflx_evap_snow(c)
        clm_idata_th%qflx_subl_snow_col(c)       = qflx_subl_snow(c)
        clm_idata_th%qflx_tran_veg_col(c)        = qflx_tran_veg(c)

        do j = 1,nlevgrnd
            clm_idata_th%qflx_rootsoil_col(c,j)  = qflx_rootsoil(c,j)
        end do

        clm_idata_th%htvp_col(c)                 = htvp(c)
        clm_idata_th%eflx_bot_col(c)             = eflx_bot(c)
        clm_idata_th%eflx_soil_grnd_col(c)       = eflx_soil_grnd(c)
        clm_idata_th%eflx_fgr0_snow_col(c)       = eflx_fgr0_snow(c)
        clm_idata_th%eflx_fgr0_h2osfc_col(c)     = eflx_fgr0_h2osfc(c)
        clm_idata_th%eflx_fgr0_soil_col(c)       = eflx_fgr0_soil(c)
        clm_idata_th%eflx_rnet_soil_col(c)       = eflx_rnet_soil(c)

    end do

    end associate
  end subroutine get_clm_soil_th_flux
!--------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------
  subroutine get_clm_bgc_state(clm_bgc_data,                    &
                        bounds, num_soilc, filter_soilc,        &
                        atm2lnd_vars, soilstate_vars,           &
                        carbonstate_vars, nitrogenstate_vars,   &
                        phosphorusstate_vars,                   &
                        ch4_vars)

    ! get clm bgc state variables
    implicit none

    type(bounds_type)           , intent(in) :: bounds
    integer                     , intent(in) :: num_soilc         ! number of soil columns in filter
    integer                     , intent(in) :: filter_soilc(:)   ! filter for soil columns

    type(atm2lnd_type)          , intent(in) :: atm2lnd_vars
    type(soilstate_type)        , intent(in) :: soilstate_vars
    type(carbonstate_type)      , intent(in) :: carbonstate_vars
    type(nitrogenstate_type)    , intent(in) :: nitrogenstate_vars
    type(phosphorusstate_type)  , intent(in) :: phosphorusstate_vars
    type(ch4_type)              , intent(in) :: ch4_vars          ! not yet used, but will be.

    type(elm_interface_bgc_datatype), intent(inout) :: clm_bgc_data

    character(len=256) :: subname = "get_clm_bgc_state"

    ! Local variables
    integer  :: fc, c, j, k

    !------------------------------------------------------------------------------------------
    !
    associate ( &
       decomp_cpools_vr=> col_cs%decomp_cpools_vr     , & ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
       decomp_npools_vr=> col_ns%decomp_npools_vr   , & ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
       decomp_ppools_vr=> col_ps%decomp_ppools_vr , & ! [real(r8) (:,:,:) ! col (gP/m3) vertically-resolved decomposing (litter, cwd, soil) P pools
       smin_no3_vr     => col_ns%smin_no3_vr        , & ! (gN/m3) vertically-resolved soil mineral NO3
       smin_nh4_vr     => col_ns%smin_nh4_vr        , & ! (gN/m3) vertically-resolved soil mineral NH4
       smin_nh4sorb_vr => col_ns%smin_nh4sorb_vr    , & ! (gN/m3) vertically-resolved soil mineral NH4 absorbed
       !
       solutionp_vr    => col_ps%solutionp_vr     , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil solution P
       labilep_vr      => col_ps%labilep_vr       , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil labile mineral P
       secondp_vr      => col_ps%secondp_vr       , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil secondary mineralP
       sminp_vr        => col_ps%sminp_vr         , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil mineral P = solutionp + labilep + secondp
       occlp_vr        => col_ps%occlp_vr         , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil occluded mineral P
       primp_vr        => col_ps%primp_vr         , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil primary mineral P
       !
       forc_pco2             => atm2lnd_vars%forc_pco2_grc               , & ! partial pressure co2 (Pa)
       forc_pch4             => atm2lnd_vars%forc_pch4_grc               , & ! partial pressure ch4 (Pa)
       !
       o2stress_sat          => ch4_vars%o2stress_sat_col                , & ! Input:  [real(r8) (:,:)   ]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
       o2stress_unsat        => ch4_vars%o2stress_unsat_col              , & ! Input:  [real(r8) (:,:)   ]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
       finundated            => ch4_vars%finundated_col                  , & ! Input:  [real(r8) (:)     ]  fractional inundated area (excluding dedicated wetland columns)
       o2_decomp_depth_unsat => ch4_vars%o2_decomp_depth_unsat_col       , & ! Input:  [real(r8) (:,:)  ]  O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
       conc_o2_unsat         => ch4_vars%conc_o2_unsat_col               , & ! Input:  [real(r8) (:,:)  ]  O2 conc in each soil layer (mol/m3) (nlevsoi)
       o2_decomp_depth_sat   => ch4_vars%o2_decomp_depth_sat_col         , & ! Input:  [real(r8) (:,:)  ]  O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
       conc_o2_sat           => ch4_vars%conc_o2_sat_col                  & ! Input:  [real(r8) (:,:)  ]  O2 conc in each soil layer (mol/m3) (nlevsoi)
    )
!

    clm_bgc_data%forc_pco2_grc(:)                   = forc_pco2(:)
    clm_bgc_data%forc_pch4_grc(:)                   = forc_pch4(:)

    do fc = 1, num_soilc
        c = filter_soilc(fc)
            do k = 1, ndecomp_pools
                clm_bgc_data%decomp_cpools_vr_col(c,:,k)    = decomp_cpools_vr(c,:,k)
                clm_bgc_data%decomp_npools_vr_col(c,:,k)    = decomp_npools_vr(c,:,k)
                clm_bgc_data%decomp_ppools_vr_col(c,:,k)    = decomp_ppools_vr(c,:,k)

            end do

            clm_bgc_data%smin_no3_vr_col(c,:)           = smin_no3_vr(c,:)
            clm_bgc_data%smin_nh4_vr_col(c,:)           = smin_nh4_vr(c,:)
            clm_bgc_data%smin_nh4sorb_vr_col(c,:)       = smin_nh4sorb_vr(c,:)

            clm_bgc_data%solutionp_vr_col(c,:)          = solutionp_vr(c,:)
            clm_bgc_data%labilep_vr_col(c,:)            = labilep_vr(c,:)
            clm_bgc_data%secondp_vr_col(c,:)            = secondp_vr(c,:)
            clm_bgc_data%sminp_vr_col(c,:)              = solutionp_vr(c,:) + labilep_vr(c,:) + secondp_vr(c,:)
            clm_bgc_data%occlp_vr_col(c,:)              = occlp_vr(c,:)
            clm_bgc_data%primp_vr_col(c,:)              = primp_vr(c,:)

            clm_bgc_data%finundated_col(c)              = finundated(c)
            clm_bgc_data%o2stress_unsat_col(c,:)        = o2stress_unsat(c,:)
            clm_bgc_data%o2stress_sat_col(c,:)          = o2stress_sat(c,:)
            clm_bgc_data%o2_decomp_depth_unsat_col(c,:) = o2_decomp_depth_unsat(c,:)
            clm_bgc_data%conc_o2_unsat_col(c,:)         = conc_o2_unsat(c,:)
            clm_bgc_data%o2_decomp_depth_sat_col(c,:)   = o2_decomp_depth_sat(c,:)
            clm_bgc_data%conc_o2_sat_col(c,:)           = conc_o2_sat(c,:)

    end do

!-----------------------------------------------------------------------------

  end associate
  end subroutine get_clm_bgc_state
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
  subroutine get_clm_bgc_flux(clm_bgc_data,                      &
                        bounds, num_soilc, filter_soilc,         &
                        cnstate_vars, carbonflux_vars,           &
                        nitrogenflux_vars, phosphorusflux_vars,  &
                        ch4_vars)

  !
  ! !DESCRIPTION:
  ! get clm bgc flux variables: external inputs to bgc state variables (pools)
  !
  ! !USES:
    use clm_time_manager      , only : get_curr_date
    use clm_varctl            , only : spinup_state
    use CNDecompCascadeConType, only : decomp_cascade_con

  ! !ARGUMENTS:
    implicit none

    type(bounds_type) , intent(in) :: bounds
    integer           , intent(in) :: num_soilc       ! number of soil columns in filter
    integer           , intent(in) :: filter_soilc(:) ! filter for soil columns

    type(cnstate_type)                  , intent(in)    :: cnstate_vars
    type(carbonflux_type)               , intent(in)    :: carbonflux_vars
    type(nitrogenflux_type)             , intent(in)    :: nitrogenflux_vars
    type(phosphorusflux_type)           , intent(in)    :: phosphorusflux_vars
    type(ch4_type)                      , intent(in) :: ch4_vars          ! not yet used, but will be.

    type(elm_interface_bgc_datatype)   , intent(inout) :: clm_bgc_data

    character(len=256) :: subname = "get_clm_bgc_flux"


 ! !LOCAL VARIABLES:
    integer  :: fc, c, g, j, k                         ! do loop indices

    real(r8) :: dtime                               ! land model time step (sec)
    integer  :: year, mon, day, sec

    ! ratios of NH4:NO3 in N deposition and fertilization (temporarily set here, will be as inputs)
    real(r8) :: r_nh4_no3_dep(bounds%begc:bounds%endc)
    real(r8) :: r_nh4_no3_fert(bounds%begc:bounds%endc)
    real(r8) :: fnh4_dep, fnh4_fert

    ! C/N source/sink rates as inputs for soil bgc

    !
    !---------------------------------------------------------------------------
    !
    associate ( &
      ! plant litering and removal + SOM/LIT vertical transport
      externalc_to_decomp_cpools_vr    => col_cf%externalc_to_decomp_cpools        , &
      externaln_to_decomp_npools_vr    => col_nf%externaln_to_decomp_npools      , &
      externalp_to_decomp_ppools_vr    => col_pf%externalp_to_decomp_ppools    , &
      t_scalar                         => col_cf%t_scalar                          , & ! Output: [real(r8) (:,:)   ]  soil temperature scalar for decomp
      w_scalar                         => col_cf%w_scalar                          , & ! Output: [real(r8) (:,:)   ]  soil water scalar for decomp
      o_scalar                         => col_cf%o_scalar                          , & ! Output: [real(r8) (:,:)   ]  fraction by which decomposition is limited by anoxia
      ! inorg. nitrogen source
      ndep_to_sminn                    => col_nf%ndep_to_sminn                   , &
      nfix_to_sminn                    => col_nf%nfix_to_sminn                   , &
      fert_to_sminn                    => col_nf%fert_to_sminn                   , &
      soyfixn_to_sminn                 => col_nf%soyfixn_to_sminn                , &
      supplement_to_sminn_vr           => col_nf%supplement_to_sminn_vr          , &
      !
      nfixation_prof                   => cnstate_vars%nfixation_prof_col                       , &
      ndep_prof                        => cnstate_vars%ndep_prof_col                            , &

      decomp_k_scalar                  => cnstate_vars%scalaravg_col                            , &

      no3_net_transport_vr             => col_nf%no3_net_transport_vr            , &
      nh4_net_transport_vr             => col_nf%nh4_net_transport_vr            , &
      col_plant_ndemand_vr             => col_nf%plant_ndemand_vr                , &

      plant_ndemand_col                => col_nf%plant_ndemand                   , &
      plant_pdemand_col                => col_pf%plant_pdemand                 , &

      col_plant_pdemand_vr             => col_pf%plant_pdemand_vr              , &
      pdep_to_sminp                    => col_pf%pdep_to_sminp                 , &
      ! assume pdep_prof = ndep_prof
      fert_p_to_sminp                  => col_pf%fert_p_to_sminp               , &
      supplement_to_sminp_vr           => col_pf%supplement_to_sminp_vr        , &
      sminp_net_transport_vr           => col_pf%sminp_net_transport_vr          &
    )

    !
    call get_curr_date(year, mon, day, sec)

!
    r_nh4_no3_dep(:)  = 1.0_r8      ! temporarily assuming half of N dep is in NH4 and another half in NO3
    r_nh4_no3_fert(:) = 1.0_r8      ! temporarily assiming half of N fertilization is in NH4 and another half in NO3

!
    do fc = 1,num_soilc
        c = filter_soilc(fc)

        clm_bgc_data%t_scalar_col(c,:) = t_scalar(c,:)
        clm_bgc_data%w_scalar_col(c,:) = w_scalar(c,:)
        clm_bgc_data%o_scalar_col(c,:) = o_scalar(c,:)


        clm_bgc_data%plant_ndemand_col(c)   = plant_ndemand_col(c)
        clm_bgc_data%plant_pdemand_col(c)   = plant_pdemand_col(c)

        fnh4_dep  = max(0._r8, min(1.0_r8, 1._r8/(r_nh4_no3_dep(c)+1._r8)))
        fnh4_fert = max(0._r8, min(1.0_r8, 1._r8/(r_nh4_no3_fert(c)+1._r8)))

        do k = 1, ndecomp_pools
            clm_bgc_data%externalc_to_decomp_cpools_col(c,:,k)  = externalc_to_decomp_cpools_vr(c,:,k)
            clm_bgc_data%externaln_to_decomp_npools_col(c,:,k)  = externaln_to_decomp_npools_vr(c,:,k)
            clm_bgc_data%externalp_to_decomp_ppools_col(c,:,k)  = externalp_to_decomp_ppools_vr(c,:,k)
       end do

       ! the following is for CTC ad-spinup.
       ! There is a 'time' control here, so MUST be called each time-step,
       ! and then better put the code here rather than in 'get_clm_bgc_state'
       if (spinup_state == 1 .and. year >= 40 .and. nu_com .eq. 'RD') then
            clm_bgc_data%sitefactor_kd_vr_col(c,:) = decomp_k_scalar(c,:)
       else
            clm_bgc_data%sitefactor_kd_vr_col(c,:) = 1.0_r8
       end if



       clm_bgc_data%externaln_to_nh4_col(c,:)          =   fnh4_dep*ndep_to_sminn(c) * ndep_prof(c,:) +  &
                                                           fnh4_fert*fert_to_sminn(c) * ndep_prof(c,:) + &
                                                           fnh4_fert*supplement_to_sminn_vr(c,:) +       &
                                                           nfix_to_sminn(c) * nfixation_prof(c,:) +      &
                                                           soyfixn_to_sminn(c) * nfixation_prof(c,:)

       clm_bgc_data%externaln_to_no3_col(c,:)          =   (1._r8-fnh4_dep)*ndep_to_sminn(c) * ndep_prof(c, :) +  &
                                                           (1._r8-fnh4_fert)*fert_to_sminn(c) * ndep_prof(c, :) + &
                                                           (1._r8-fnh4_fert)*supplement_to_sminn_vr(c,:)



       clm_bgc_data%externalp_to_primp_col(c,:)        =   pdep_to_sminp(c)*ndep_prof(c, :)
       clm_bgc_data%externalp_to_labilep_col(c,:)      =   fert_p_to_sminp(c)*ndep_prof(c, :)
       clm_bgc_data%externalp_to_solutionp_col(c,:)    =   supplement_to_sminp_vr(c,:)

       clm_bgc_data%no3_net_transport_vr_col(c,:)      = no3_net_transport_vr(c,:)
       clm_bgc_data%nh4_net_transport_vr_col(c,:)      = nh4_net_transport_vr(c,:)
       clm_bgc_data%sminp_net_transport_vr_col(c,:)    = sminp_net_transport_vr(c,:)  !from solutionp

       clm_bgc_data%plant_ndemand_vr_col(c,:)          = col_plant_ndemand_vr(c,:)
       clm_bgc_data%plant_pdemand_vr_col(c,:)          = col_plant_pdemand_vr(c,:)

    end do ! fc = 1,num_soilc

    end associate
  end subroutine get_clm_bgc_flux
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
  subroutine update_soil_moisture(clm_idata_th,     &
           bounds, num_soilc, filter_soilc,   &
           soilstate_vars, waterstate_vars)

  !
  ! !DESCRIPTION:
  !
  !
  ! !USES:

  ! !ARGUMENTS:
    implicit none

    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_soilc        ! number of column soil points in column filter
    integer, intent(in) :: filter_soilc(:)  ! column filter for soil points
    type(soilstate_type), intent(inout)  :: soilstate_vars
    type(waterstate_type), intent(inout) :: waterstate_vars

    type(elm_interface_th_datatype), intent(in) :: clm_idata_th

  ! !LOCAL VARIABLES:
    integer  :: fc, c, j, g, gcount      ! indices

  !EOP
  !-----------------------------------------------------------------------
    associate ( &
         soilpsi_col    =>  soilstate_vars%soilpsi_col          , &
         !
         h2osoi_liq_col =>  col_ws%h2osoi_liq      , &
         h2osoi_ice_col =>  col_ws%h2osoi_ice      , &
         h2osoi_vol_col =>  col_ws%h2osoi_vol        &
    )

    do fc = 1,num_soilc
        c = filter_soilc(fc)

        soilpsi_col(c,:)    =  clm_idata_th%soilpsi_col(c,:)

        h2osoi_liq_col(c,:) =  clm_idata_th%h2osoi_liq_col(c,:)
        h2osoi_ice_col(c,:) =  clm_idata_th%h2osoi_ice_col(c,:)
        h2osoi_vol_col(c,:) =  clm_idata_th%h2osoi_vol_col(c,:)
    end do

    end associate
  end subroutine update_soil_moisture
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
  subroutine update_soil_temperature(clm_idata_th,     &
           bounds, num_soilc, filter_soilc,            &
           temperature_vars)

  !
  ! !DESCRIPTION:
  !
  !
  ! !USES:

  ! !ARGUMENTS:
    implicit none

    type(bounds_type)                   , intent(in)    :: bounds
    integer                             , intent(in)    :: num_soilc        ! number of column soil points in column filter
    integer                             , intent(in)    :: filter_soilc(:)  ! column filter for soil points
    type(temperature_type)              , intent(inout) :: temperature_vars
    type(elm_interface_th_datatype)     , intent(in)    :: clm_idata_th

  ! !LOCAL VARIABLES:
    integer  :: fc, c, j, g, gcount      ! indices

  !EOP
  !-----------------------------------------------------------------------
    associate ( &
         t_soisno   => col_es%t_soisno   & ! snow-soil temperature (Kelvin)
    )

    do fc = 1,num_soilc
        c = filter_soilc(fc)
        t_soisno(c,:)   = clm_idata_th%t_soisno_col(c,:)
    end do

    end associate
  end subroutine update_soil_temperature
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
  subroutine update_th_data_pf2clm(clm_idata_th,           &
           bounds, num_soilc, filter_soilc,                &
           waterstate_vars, waterflux_vars,                &
           temperature_vars, energyflux_vars,              &
           soilstate_vars, soilhydrology_vars)

    ! USES
    use clm_varctl          , only : use_pflotran, pf_tmode, pf_hmode

    implicit none

    ! !ARGUMENTS:
    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                     , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    type(waterstate_type)       , intent(inout) :: waterstate_vars
    type(waterflux_type)        , intent(inout) :: waterflux_vars
    type(temperature_type)      , intent(inout) :: temperature_vars
    type(soilstate_type)        , intent(inout) :: soilstate_vars
    type(soilhydrology_type)    , intent(inout) :: soilhydrology_vars
    type(energyflux_type)       , intent(inout) :: energyflux_vars

    type(elm_interface_th_datatype), intent(in) :: clm_idata_th

    !-----------------------------------------------------------------------

    character(len=256) :: subname = "update_th_data_pf2clm"

    if (pf_tmode) then
        call update_soil_temperature(clm_idata_th,      &
                   bounds, num_soilc, filter_soilc,     &
                   temperature_vars)
    end if

    if (pf_hmode) then
        call update_soil_moisture(clm_idata_th,         &
                   bounds, num_soilc, filter_soilc,     &
                   soilstate_vars, waterstate_vars)
    end if

  end subroutine update_th_data_pf2clm
!--------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------
  subroutine update_bgc_state_decomp(clm_bgc_data,  &
           bounds, num_soilc, filter_soilc,         &
           carbonstate_vars, nitrogenstate_vars,    &
           phosphorusstate_vars                     &
           )

    implicit none

    type(bounds_type)                   , intent(in)    :: bounds
    integer                             , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                             , intent(in)    :: filter_soilc(:)   ! filter for soil columns

    type(carbonstate_type)              , intent(inout) :: carbonstate_vars
    type(nitrogenstate_type)            , intent(inout) :: nitrogenstate_vars
    type(phosphorusstate_type)          , intent(inout) :: phosphorusstate_vars

    type(elm_interface_bgc_datatype)    , intent(in)    :: clm_bgc_data

    character(len=256) :: subname = "update_soil_bgc_state"

    integer  :: fc,c,j,k

!------------------------------------------------------------------------------------
    !
    associate ( &
       decomp_cpools_vr             => col_cs%decomp_cpools_vr           , &
       decomp_npools_vr             => col_ns%decomp_npools_vr         , &
       decomp_ppools_vr             => col_ps%decomp_ppools_vr         &
    )
! ------------------------------------------------------------------------
!
    do fc = 1, num_soilc
        c = filter_soilc(fc)
            do k = 1, ndecomp_pools
                decomp_cpools_vr(c,:,k) = clm_bgc_data%decomp_cpools_vr_col(c,:,k)
                decomp_npools_vr(c,:,k) = clm_bgc_data%decomp_npools_vr_col(c,:,k)
                decomp_ppools_vr(c,:,k) = clm_bgc_data%decomp_ppools_vr_col(c,:,k)
            end do
    end do

    end associate
  end subroutine update_bgc_state_decomp
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
    subroutine update_bgc_state_smin(clm_bgc_data,      &
           bounds, num_soilc, filter_soilc,             &
           nitrogenstate_vars, phosphorusstate_vars)

    use CNDecompCascadeConType, only : decomp_cascade_con
    use clm_time_manager, only : get_step_size

    implicit none

    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                     , intent(in)    :: filter_soilc(:)   ! filter for soil columns

    type(nitrogenstate_type)    , intent(inout) :: nitrogenstate_vars
    type(phosphorusstate_type)  , intent(inout) :: phosphorusstate_vars

    type(elm_interface_bgc_datatype), intent(in) :: clm_bgc_data

    character(len=256) :: subname = "update_bgc_state_smin"

    integer  :: fc,c,j

!------------------------------------------------------------------------------------
     !
     associate ( &
     sminn_vr           => col_ns%sminn_vr           , &
     smin_no3_vr        => col_ns%smin_no3_vr        , &
     smin_nh4_vr        => col_ns%smin_nh4_vr        , &
     smin_nh4sorb_vr    => col_ns%smin_nh4sorb_vr    , &

     solutionp_vr       => col_ps%solutionp_vr     , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil solution P
     labilep_vr         => col_ps%labilep_vr       , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil labile mineral P
     secondp_vr         => col_ps%secondp_vr       , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil secondary mineralP
     sminp_vr           => col_ps%sminp_vr         , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil mineral P = solutionp + labilep + second
     occlp_vr           => col_ps%occlp_vr         , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil occluded mineral P
     primp_vr           => col_ps%primp_vr           & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil primary mineral P

     )
! ------------------------------------------------------------------------
!
    do fc = 1, num_soilc
        c = filter_soilc(fc)
            smin_no3_vr(c,:)        = clm_bgc_data%smin_no3_vr_col(c,:)
            smin_nh4_vr(c,:)        = clm_bgc_data%smin_nh4_vr_col(c,:)
            smin_nh4sorb_vr(c,:)    = clm_bgc_data%smin_nh4sorb_vr_col(c,:)
            sminn_vr(c,:)           = clm_bgc_data%sminn_vr_col(c,:)

            solutionp_vr(c,:)       = clm_bgc_data%solutionp_vr_col(c,:)
            labilep_vr(c,:)         = clm_bgc_data%labilep_vr_col(c,:)
            secondp_vr(c,:)         = clm_bgc_data%secondp_vr_col(c,:)
            sminp_vr(c,:)           = clm_bgc_data%sminp_vr_col(c,:)
            occlp_vr(c,:)           = clm_bgc_data%occlp_vr_col(c,:)
            primp_vr(c,:)           = clm_bgc_data%primp_vr_col(c,:)
    end do

    end associate
  end subroutine update_bgc_state_smin
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
    subroutine update_bgc_flux_decomp_sourcesink(clm_bgc_data,       &
           bounds, num_soilc, filter_soilc,                          &
           carbonflux_vars, nitrogenflux_vars,                       &
           phosphorusflux_vars)

    use CNDecompCascadeConType, only : decomp_cascade_con

    implicit none

    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                     , intent(in)    :: filter_soilc(:)   ! filter for soil columns

    type(carbonflux_type)       , intent(inout) :: carbonflux_vars
    type(nitrogenflux_type)     , intent(inout) :: nitrogenflux_vars
    type(phosphorusflux_type)   , intent(inout) :: phosphorusflux_vars

    type(elm_interface_bgc_datatype), intent(in):: clm_bgc_data

    integer :: fc, c, j, k
    character(len=256) :: subname = "update_soil_bgc_pf2clm"

    associate ( &
     decomp_cpools_sourcesink_vr  => col_cf%decomp_cpools_sourcesink    , &
     decomp_npools_sourcesink_vr  => col_nf%decomp_npools_sourcesink  , &
     decomp_ppools_sourcesink_vr  => col_pf%decomp_ppools_sourcesink  &
     )

    do fc = 1, num_soilc
        c = filter_soilc(fc)
            do k = 1, ndecomp_pools
                decomp_cpools_sourcesink_vr(c,:,k) = clm_bgc_data%decomp_cpools_sourcesink_col(c,:,k)
                decomp_npools_sourcesink_vr(c,:,k) = clm_bgc_data%decomp_npools_sourcesink_col(c,:,k)
                decomp_ppools_sourcesink_vr(c,:,k) = clm_bgc_data%decomp_ppools_sourcesink_col(c,:,k)
            end do
    end do
    end associate
    end subroutine update_bgc_flux_decomp_sourcesink
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
    subroutine update_bgc_flux_decomp_cascade(clm_bgc_data, &
           bounds, num_soilc, filter_soilc,                 &
           carbonflux_vars, nitrogenflux_vars,              &
           phosphorusflux_vars)

    use CNDecompCascadeConType, only : decomp_cascade_con

    implicit none

    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                     , intent(in)    :: filter_soilc(:)   ! filter for soil columns

    type(carbonflux_type)       , intent(inout) :: carbonflux_vars
    type(nitrogenflux_type)     , intent(inout) :: nitrogenflux_vars
    type(phosphorusflux_type)   , intent(inout) :: phosphorusflux_vars

    type(elm_interface_bgc_datatype), intent(in):: clm_bgc_data

    integer :: fc, c, j, k
    character(len=256) :: subname = "update_soil_bgc_pf2clm"

    associate ( &
     phr_vr                           => col_cf%phr_vr                                 , & ! Output: [real(r8) (:,:)   ]  potential HR (gC/m3/s)
     fphr                             => col_cf%fphr                                   , & ! Output: [real(r8) (:,:)   ]  fraction of potential SOM + LITTER heterotrophic

     decomp_cascade_hr_vr_col         => col_cf%decomp_cascade_hr_vr                   , &

     decomp_cascade_ctransfer_vr_col  => col_cf%decomp_cascade_ctransfer_vr            , &
     decomp_cascade_ntransfer_vr_col  => col_nf%decomp_cascade_ntransfer_vr          , &
     decomp_cascade_ptransfer_vr_col  => col_pf%decomp_cascade_ptransfer_vr        , &

     decomp_cascade_sminn_flux_vr_col => col_nf%decomp_cascade_sminn_flux_vr         , & ! Output: [real(r8) (:,:,:) ]  vert-res mineral N flux for transition along decomposition cascade (gN/m3/s)
     decomp_cascade_sminp_flux_vr_col => col_pf%decomp_cascade_sminp_flux_vr       , & ! Output: [real(r8) (:,:,:) ]  vert-res mineral P flux for transition along decomposition cascade (gP/m3/s)

     sminn_to_denit_decomp_cascade_vr_col => col_nf%sminn_to_denit_decomp_cascade_vr   & ! Output: [real(r8) (:,:,:) ]
     )

    do fc = 1, num_soilc
        c = filter_soilc(fc)
            phr_vr(c,:) = clm_bgc_data%phr_vr_col(c,:)
            fphr(c,:)   = clm_bgc_data%fphr_col(c,:)

            do k = 1, ndecomp_cascade_transitions
                decomp_cascade_hr_vr_col(c,:,k)         = clm_bgc_data%decomp_cascade_hr_vr_col(c,:,k)
                decomp_cascade_ctransfer_vr_col(c,:,k)  = clm_bgc_data%decomp_cascade_ctransfer_vr_col(c,:,k)
                decomp_cascade_ntransfer_vr_col(c,:,k)  = clm_bgc_data%decomp_cascade_ntransfer_vr_col(c,:,k)
                decomp_cascade_ptransfer_vr_col(c,:,k)  = clm_bgc_data%decomp_cascade_ptransfer_vr_col(c,:,k)

                decomp_cascade_sminn_flux_vr_col(c,:,k)     = clm_bgc_data%decomp_cascade_sminn_flux_vr_col(c,:,k)
                decomp_cascade_sminp_flux_vr_col(c,:,k)     = clm_bgc_data%decomp_cascade_sminp_flux_vr_col(c,:,k)
                sminn_to_denit_decomp_cascade_vr_col(c,:,k) = clm_bgc_data%sminn_to_denit_decomp_cascade_vr_col(c,:,k)
            end do
    end do
    end associate
    end subroutine update_bgc_flux_decomp_cascade
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
    subroutine update_bgc_flux_smin(clm_bgc_data,   &
           bounds, num_soilc, filter_soilc,         &
           cnstate_vars,                            &
           nitrogenflux_vars, phosphorusflux_vars)

    implicit none

    type(bounds_type)                   , intent(in)    :: bounds
    integer                             , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                             , intent(in)    :: filter_soilc(:)   ! filter for soil columns

    type(cnstate_type)                  , intent(inout) :: cnstate_vars

    type(nitrogenflux_type)             , intent(inout) :: nitrogenflux_vars
    type(phosphorusflux_type)           , intent(inout) :: phosphorusflux_vars
    type(elm_interface_bgc_datatype)    , intent(in)    :: clm_bgc_data

    integer :: fc, c, j
    character(len=256) :: subname = "update_bgc_flux_smin"

    associate ( &
     fpg                          => cnstate_vars%fpg_col                            , & ! Output: [real(r8) (:)   ]  fraction of potential gpp (no units)
     fpi                          => cnstate_vars%fpi_col                            , & ! Output: [real(r8) (:)   ]  fraction of potential immobilization (no units)
     fpi_vr                       => cnstate_vars%fpi_vr_col                         , & ! Output: [real(r8) (:,:) ]  fraction of potential immobilization (no units)
     fpg_p                        => cnstate_vars%fpg_p_col                          , & ! Output: [real(r8) (:)   ]  fraction of potential gpp (no units)
     fpi_p                        => cnstate_vars%fpi_p_col                          , & ! Output: [real(r8) (:)   ]  fraction of potential immobilization (no units)
     fpi_p_vr                     => cnstate_vars%fpi_p_vr_col                       , & ! Output: [real(r8) (:,:) ]  fraction of potential immobilization (no units)

     potential_immob              => col_nf%potential_immob           , & ! Output: [real(r8) (:)   ]
     actual_immob                 => col_nf%actual_immob              , & ! Output: [real(r8) (:)   ]
     sminn_to_plant               => col_nf%sminn_to_plant            , & ! Output: [real(r8) (:)   ]

     sminn_to_plant_vr            => col_nf%sminn_to_plant_vr         , &
     smin_no3_to_plant_vr         => col_nf%smin_no3_to_plant_vr      , &
     smin_nh4_to_plant_vr         => col_nf%smin_nh4_to_plant_vr      , &
     potential_immob_vr           => col_nf%potential_immob_vr        , &
     actual_immob_vr              => col_nf%actual_immob_vr           , &
     actual_immob_no3_vr          => col_nf%actual_immob_no3_vr       , & ! Output: [real(r8) (:,:) ]
     actual_immob_nh4_vr          => col_nf%actual_immob_nh4_vr       , & ! Output: [real(r8) (:,:) ]
     gross_nmin_vr                => col_nf%gross_nmin_vr             , &
     net_nmin_vr                  => col_nf%net_nmin_vr               , & ! Output: [real(r8) (:,:)   ]

     sminn_to_denit_excess_vr     => col_nf%sminn_to_denit_excess_vr  , & ! Output: [real(r8) (:,:) ]
     supplement_to_sminn_vr       => col_nf%supplement_to_sminn_vr    , & ! Output: [real(r8) (:,:) ]

     no3_net_transport_vr         => col_nf%no3_net_transport_vr      , & ! Output: updated from PF, if coupled
     nh4_net_transport_vr         => col_nf%nh4_net_transport_vr      , & ! Output: updated from PF, if coupled

     potential_immob_p            => col_pf%potential_immob_p       , & ! Output: [real(r8) (:)   ]
     actual_immob_p               => col_pf%actual_immob_p          , & ! Output: [real(r8) (:)   ]
     sminp_to_plant               => col_pf%sminp_to_plant          , & ! Output: [real(r8) (:)   ]

     supplement_to_sminp_vr       => col_pf%supplement_to_sminp_vr  , & ! Output: [real(r8) (:,:) ]

     sminp_to_plant_vr            => col_pf%sminp_to_plant_vr       , &
     potential_immob_p_vr         => col_pf%potential_immob_p_vr    , & ! Output: [real(r8) (:,:)   ]
     actual_immob_p_vr            => col_pf%actual_immob_p_vr       , &
     gross_pmin_vr                => col_pf%gross_pmin_vr           , & ! Output: [real(r8) (:,:)   ]
     net_pmin_vr                  => col_pf%net_pmin_vr               & ! Output: [real(r8) (:,:)   ]
     )

     do fc = 1, num_soilc
        c = filter_soilc(fc)

        fpg(c)              = clm_bgc_data%fpg_col(c)
        fpi(c)              = clm_bgc_data%fpi_col(c)
        fpg_p(c)            = clm_bgc_data%fpg_p_col(c)
        fpi_p(c)            = clm_bgc_data%fpi_p_col(c)

        potential_immob(c)  = clm_bgc_data%potential_immob_col(c)
        actual_immob(c)     = clm_bgc_data%actual_immob_col(c)
        sminn_to_plant(c)   = clm_bgc_data%sminn_to_plant_col(c)

        potential_immob_p(c)  = clm_bgc_data%potential_immob_p_col(c)
        actual_immob_p(c)     = clm_bgc_data%actual_immob_p_col(c)
        sminp_to_plant(c)     = clm_bgc_data%sminp_to_plant_col(c)

        fpi_vr(c,:)                 = clm_bgc_data%fpi_vr_col(c,:)
        fpi_p_vr(c,:)               = clm_bgc_data%fpi_p_vr_col(c,:)

        sminn_to_plant_vr(c,:)      = clm_bgc_data%sminn_to_plant_vr_col(c,:)
        smin_no3_to_plant_vr(c,:)   = clm_bgc_data%smin_no3_to_plant_vr_col(c,:)
        smin_nh4_to_plant_vr(c,:)   = clm_bgc_data%smin_nh4_to_plant_vr_col(c,:)

        potential_immob_vr(c,:)     = clm_bgc_data%potential_immob_vr_col(c,:)
        actual_immob_vr(c,:)        = clm_bgc_data%actual_immob_vr_col(c,:)
        actual_immob_no3_vr(c,:)    = clm_bgc_data%actual_immob_no3_vr_col(c,:)
        actual_immob_nh4_vr(c,:)    = clm_bgc_data%actual_immob_nh4_vr_col(c,:)
        gross_nmin_vr(c,:)          = clm_bgc_data%gross_nmin_vr_col(c,:)
        net_nmin_vr(c,:)            = clm_bgc_data%net_nmin_vr_col(c,:)

        sminn_to_denit_excess_vr(c,:)=clm_bgc_data%sminn_to_denit_excess_vr_col(c,:)
        supplement_to_sminn_vr(c,:) = clm_bgc_data%supplement_to_sminn_vr_col(c,:)

        no3_net_transport_vr(c,:)   = clm_bgc_data%no3_net_transport_vr_col(c,:)
        nh4_net_transport_vr(c,:)   = clm_bgc_data%nh4_net_transport_vr_col(c,:)

        supplement_to_sminp_vr(c,:) = clm_bgc_data%supplement_to_sminp_vr_col(c,:)

        sminp_to_plant_vr(c,:)      = clm_bgc_data%sminp_to_plant_vr_col(c,:)
        potential_immob_p_vr(c,:)   = clm_bgc_data%potential_immob_p_vr_col(c,:)
        actual_immob_p_vr(c,:)      = clm_bgc_data%actual_immob_p_vr_col(c,:)
        gross_pmin_vr(c,:)          = clm_bgc_data%gross_pmin_vr_col(c,:)
        net_pmin_vr(c,:)            = clm_bgc_data%net_pmin_vr_col(c,:)     !NOT available in PF

      end do
    end associate
    end subroutine update_bgc_flux_smin
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
    subroutine update_bgc_flux_nitdenit(clm_bgc_data,   &
           bounds, num_soilc, filter_soilc,             &
           nitrogenflux_vars, phosphorusflux_vars)

    implicit none

    type(bounds_type)                   , intent(in)    :: bounds
    integer                             , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                             , intent(in)    :: filter_soilc(:)   ! filter for soil columns

    type(nitrogenflux_type)             , intent(inout) :: nitrogenflux_vars
    type(phosphorusflux_type)           , intent(inout) :: phosphorusflux_vars
    type(elm_interface_bgc_datatype)    , intent(in)    :: clm_bgc_data

    integer :: fc, c, j
    character(len=256) :: subname = "update_bgc_flux_nitdenit"

    associate ( &
         pot_f_nit_vr                 => col_nf%pot_f_nit_vr                  , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) potential soil nitrification flux
         pot_f_denit_vr               => col_nf%pot_f_denit_vr                , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) potential soil denitrification flux
         f_nit_vr                     => col_nf%f_nit_vr                      , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) soil nitrification flux
         f_denit_vr                   => col_nf%f_denit_vr                    , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) soil denitrification flux
         n2_n2o_ratio_denit_vr        => col_nf%n2_n2o_ratio_denit_vr         , & ! Output: [real(r8) (:,:) ]  ratio of N2 to N2O production by denitrification [gN/gN]
         f_n2o_denit_vr               => col_nf%f_n2o_denit_vr                , & ! Output: [real(r8) (:,:) ]  flux of N2O from denitrification [gN/m3/s]
         f_n2o_nit_vr                 => col_nf%f_n2o_nit_vr                    & ! Output: [real(r8) (:,:) ]  flux of N2O from nitrification [gN/m3/s]
    )

    do fc = 1, num_soilc
        c = filter_soilc(fc)
        pot_f_nit_vr(c,:)           = clm_bgc_data%pot_f_nit_vr_col(c,:)
        pot_f_denit_vr(c,:)         = clm_bgc_data%pot_f_denit_vr_col(c,:)
        f_nit_vr(c,:)               = clm_bgc_data%f_nit_vr_col(c,:)
        f_denit_vr(c,:)             = clm_bgc_data%f_denit_vr_col(c,:)
        n2_n2o_ratio_denit_vr(c,:)  = clm_bgc_data%n2_n2o_ratio_denit_vr_col(c,:)
        f_n2o_denit_vr(c,:)         = clm_bgc_data%f_n2o_denit_vr_col(c,:)
        f_n2o_nit_vr(c,:)           = clm_bgc_data%f_n2o_nit_vr_col(c,:)
    end do
    end associate
    end subroutine update_bgc_flux_nitdenit
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
  subroutine update_bgc_flux_gas_pf(clm_bgc_data,  &
     bounds, num_soilc, filter_soilc,              &
     carbonflux_vars, nitrogenflux_vars)

     ! PFLOTRAN gas fluxes
     implicit none

     type(bounds_type)                  , intent(in)    :: bounds
     integer                            , intent(in)    :: num_soilc       ! number of soil columns in filter
     integer                            , intent(in)    :: filter_soilc(:) ! filter for soil columns

     type(carbonflux_type)              , intent(inout) :: carbonflux_vars
     type(nitrogenflux_type)            , intent(inout) :: nitrogenflux_vars
     type(elm_interface_bgc_datatype)   , intent(in)    :: clm_bgc_data

     !character(len=256) :: subname = "get_pf_bgc_gaslosses"

     integer  :: fc, c, g, j

!------------------------------------------------------------------------------------
    associate ( &
     hr_vr                        => col_cf%hr_vr              , &
     f_co2_soil_vr                => col_cf%f_co2_soil_vr      , &
     f_n2o_soil_vr                => col_nf%f_n2o_soil_vr    , &
     f_n2_soil_vr                 => col_nf%f_n2_soil_vr     , &
     f_ngas_decomp_vr             => col_nf%f_ngas_decomp_vr , &
     f_ngas_nitri_vr              => col_nf%f_ngas_nitri_vr  , &
     f_ngas_denit_vr              => col_nf%f_ngas_denit_vr    &
     )
! ------------------------------------------------------------------------
    do fc = 1,num_soilc
        c = filter_soilc(fc)
        f_co2_soil_vr(c,:)         = clm_bgc_data%f_co2_soil_vr_col(c,:)
        f_n2_soil_vr(c,:)          = clm_bgc_data%f_n2_soil_vr_col(c,:)
        f_n2o_soil_vr(c,:)         = clm_bgc_data%f_n2o_soil_vr_col(c,:)

        hr_vr(c,:)                 = clm_bgc_data%hr_vr_col(c,:)
        f_ngas_decomp_vr(c,:)      = clm_bgc_data%f_ngas_decomp_vr_col(c,:)
        f_ngas_nitri_vr(c,:)       = clm_bgc_data%f_ngas_nitri_vr_col(c,:)
        f_ngas_denit_vr(c,:)       = clm_bgc_data%f_ngas_denit_vr_col(c,:)

     enddo ! do c = begc, endc
!
    end associate
  end subroutine update_bgc_flux_gas_pf
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
  subroutine update_bgc_data_pf2clm(clm_bgc_data, bounds,         &
           num_soilc, filter_soilc,                               &
           num_soilp, filter_soilp,                               &
           cnstate_vars, carbonflux_vars, carbonstate_vars,       &
           nitrogenflux_vars, nitrogenstate_vars,                 &
           phosphorusflux_vars, phosphorusstate_vars,             &
           ch4_vars)
    ! USES
    use clm_varctl          , only : use_pflotran, pf_cmode

    implicit none

    ! !ARGUMENTS:
    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                     , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                     , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                     , intent(in)    :: filter_soilp(:)   ! filter for soil patches

    type(cnstate_type)          , intent(inout) :: cnstate_vars
    type(carbonflux_type)       , intent(inout) :: carbonflux_vars
    type(carbonstate_type)      , intent(inout) :: carbonstate_vars
    type(nitrogenflux_type)     , intent(inout) :: nitrogenflux_vars
    type(nitrogenstate_type)    , intent(inout) :: nitrogenstate_vars
    type(phosphorusflux_type)   , intent(inout) :: phosphorusflux_vars
    type(phosphorusstate_type)  , intent(inout) :: phosphorusstate_vars
    type(ch4_type)              , intent(inout) :: ch4_vars

    type(elm_interface_bgc_datatype), intent(in):: clm_bgc_data

    !-----------------------------------------------------------------------

    character(len=256) :: subname = "update_bgc_data_pf2clm"


    if (pf_cmode) then
        ! bgc_state_decomp is updated in CLM
        ! by passing bgc_flux_decomp_sourcesink into SoilLittVertTransp
        call update_bgc_flux_decomp_sourcesink(clm_bgc_data,    &
                    bounds, num_soilc, filter_soilc,            &
                    carbonflux_vars, nitrogenflux_vars,         &
                    phosphorusflux_vars)

        call update_bgc_state_smin(clm_bgc_data,                &
                    bounds, num_soilc, filter_soilc,            &
                    nitrogenstate_vars, phosphorusstate_vars)

        call update_bgc_flux_smin(clm_bgc_data,                 &
                    bounds, num_soilc, filter_soilc,            &
                    cnstate_vars,                               &
                    nitrogenflux_vars, phosphorusflux_vars)

        call update_bgc_flux_gas_pf(clm_bgc_data,               &
                    bounds, num_soilc, filter_soilc,            &
                    carbonflux_vars, nitrogenflux_vars)

    end if

  end subroutine update_bgc_data_pf2clm

!--------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------
! BEG of CLM-bgc through interface
!--------------------------------------------------------------------------------------
  ! !INTERFACE:
  subroutine clm_bgc_run(clm_interface_data, bounds,        &
                num_soilc, filter_soilc,                    &
                num_soilp, filter_soilp,                    &
                canopystate_vars, soilstate_vars,           &
                temperature_vars, waterstate_vars,          &
                cnstate_vars, ch4_vars,                     &
                carbonstate_vars, carbonflux_vars,          &
                nitrogenstate_vars, nitrogenflux_vars,      &
                phosphorusstate_vars,phosphorusflux_vars)

    ! USES:
    use SoilLittDecompMod          , only: SoilLittDecompAlloc

    ! ARGUMENTS:
    type(bounds_type)                   , intent(in)    :: bounds
    integer                             , intent(in)    :: num_soilc          ! number of soil columns in filter
    integer                             , intent(in)    :: filter_soilc(:)    ! filter for soil columns
    integer                             , intent(in)    :: num_soilp          ! number of soil patches in filter
    integer                             , intent(in)    :: filter_soilp(:)    ! filter for soil patches
    type(canopystate_type)              , intent(inout) :: canopystate_vars
    type(soilstate_type)                , intent(inout) :: soilstate_vars
    type(temperature_type)              , intent(inout) :: temperature_vars
    type(waterstate_type)               , intent(inout) :: waterstate_vars
    type(cnstate_type)                  , intent(inout) :: cnstate_vars
    type(ch4_type)                      , intent(inout) :: ch4_vars
    type(carbonstate_type)              , intent(inout) :: carbonstate_vars
    type(carbonflux_type)               , intent(inout) :: carbonflux_vars
    type(nitrogenstate_type)            , intent(inout) :: nitrogenstate_vars
    type(nitrogenflux_type)             , intent(inout) :: nitrogenflux_vars
    type(phosphorusstate_type)          , intent(inout) :: phosphorusstate_vars
    type(phosphorusflux_type)           , intent(inout) :: phosphorusflux_vars

    type(elm_interface_data_type)       , intent(inout) :: clm_interface_data

    !-------------------------------------------------------------
    ! STEP-2: (i) pass data from clm_bgc_data to SoilLittDecompAlloc
    call clm_bgc_get_data(clm_interface_data, bounds,       &
                num_soilc, filter_soilc,                    &
                canopystate_vars, soilstate_vars,           &
                temperature_vars, waterstate_vars,          &
                cnstate_vars, ch4_vars,                     &
                carbonstate_vars, carbonflux_vars,          &
                nitrogenstate_vars, nitrogenflux_vars,      &
                phosphorusstate_vars,phosphorusflux_vars)

    ! STEP-2: (ii) run SoilLittDecompAlloc
    call SoilLittDecompAlloc (bounds, num_soilc, filter_soilc,    &
               num_soilp, filter_soilp,                     &
               canopystate_vars, soilstate_vars,            &
               temperature_vars, waterstate_vars,           &
               cnstate_vars, ch4_vars,                      &
               carbonstate_vars, carbonflux_vars,           &
               nitrogenstate_vars, nitrogenflux_vars,       &
               phosphorusstate_vars,phosphorusflux_vars)

    ! STEP-2: (iii) update clm_bgc_data from SoilLittDecompAlloc
    call clm_bgc_update_data(clm_interface_data%bgc, bounds, &
                num_soilc, filter_soilc,                     &
                cnstate_vars, carbonflux_vars,               &
                nitrogenflux_vars, phosphorusflux_vars)

  end subroutine clm_bgc_run
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
  ! !INTERFACE:
  ! pass data from clm_bgc_data to clm original data-types that used by SoilLittDecompAlloc
  subroutine clm_bgc_get_data(clm_interface_data,       &
            bounds, num_soilc, filter_soilc,            &
            canopystate_vars, soilstate_vars,           &
            temperature_vars, waterstate_vars,          &
            cnstate_vars, ch4_vars,                     &
            carbonstate_vars, carbonflux_vars,          &
            nitrogenstate_vars, nitrogenflux_vars,      &
            phosphorusstate_vars,phosphorusflux_vars)

    ! USES:


    ! ARGUMENTS:
    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: num_soilc          ! number of soil columns in filter
    integer                     , intent(in)    :: filter_soilc(:)    ! filter for soil columns
    type(canopystate_type)      , intent(inout) :: canopystate_vars
    type(soilstate_type)        , intent(inout) :: soilstate_vars
    type(temperature_type)      , intent(inout) :: temperature_vars
    type(waterstate_type)       , intent(inout) :: waterstate_vars
    type(cnstate_type)          , intent(inout) :: cnstate_vars
    type(ch4_type)              , intent(inout) :: ch4_vars
    type(carbonstate_type)      , intent(inout) :: carbonstate_vars
    type(carbonflux_type)       , intent(inout) :: carbonflux_vars
    type(nitrogenstate_type)    , intent(inout) :: nitrogenstate_vars
    type(nitrogenflux_type)     , intent(inout) :: nitrogenflux_vars
    type(phosphorusstate_type)  , intent(inout) :: phosphorusstate_vars
    type(phosphorusflux_type)   , intent(inout) :: phosphorusflux_vars

    type(elm_interface_data_type), intent(in)   :: clm_interface_data

    ! LOCAL VARIABLES:
    integer :: fc, c, j, k
    !-----------------------------------------------------------------------

    associate(&
        decomp_cpools_vr        => col_cs%decomp_cpools_vr        , & ! Input:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
        decomp_npools_vr        => col_ns%decomp_npools_vr      , & ! Input:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
        decomp_ppools_vr        => col_ps%decomp_ppools_vr    , & ! Input:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) P pools

        smin_no3_vr             => col_ns%smin_no3_vr           , &      ! (gN/m3) vertically-resolved soil mineral NO3
        smin_nh4_vr             => col_ns%smin_nh4_vr           , &      ! (gN/m3) vertically-resolved soil mineral NH4
        smin_nh4sorb_vr         => col_ns%smin_nh4sorb_vr       , &      ! (gN/m3) vertically-resolved soil mineral NH4 absorbed

        solutionp_vr            => col_ps%solutionp_vr        , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil solution P
        labilep_vr              => col_ps%labilep_vr          , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil labile mineral P
        secondp_vr              => col_ps%secondp_vr          , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil secondary mineralP
        sminp_vr                => col_ps%sminp_vr            , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil mineral P = solutionp + labilep + secondp
        occlp_vr                => col_ps%occlp_vr            , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil occluded mineral P
        primp_vr                => col_ps%primp_vr            , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil primary mineral P
        !
        plant_ndemand_col       => col_nf%plant_ndemand          , &
        plant_pdemand_col       => col_pf%plant_pdemand        , &
        !
        !alt_indx                => canopystate_vars%alt_indx_col                , & ! Input:  [integer  (:)     ]  current depth of thaw
        !
        watsat                  => soilstate_vars%watsat_col                    , & ! Input:  [real(r8) (:,:)  ]  volumetric soil water at saturation (porosity) (nlevgrnd)
        bd                      => soilstate_vars%bd_col                        , & ! Input:  [real(r8) (:,:)  ]  bulk density of dry soil material [kg/m3]
        watfc                   => soilstate_vars%watfc_col                     , & ! Input:  [real(r8) (:,:)  ]  volumetric soil water at field capacity (nlevsoi)
        bsw                     => soilstate_vars%bsw_col                       , & ! Input:  [real(r8) (:,:)  ]  Clapp and Hornberger "b" (nlevgrnd)
        cellorg                 => soilstate_vars%cellorg_col                   , & ! Input:  [real(r8) (:,:)  ]  column 3D org (kg/m3 organic matter) (nlevgrnd)
        sucsat                  => soilstate_vars%sucsat_col                    , & ! Input:  [real(r8) (:,:)  ]  minimum soil suction (mm)
        soilpsi                 => soilstate_vars%soilpsi_col                   , & ! Input:  [real(r8) (:,:)  ]  soil water potential in each soil layer (MPa)

        h2osoi_vol              => col_ws%h2osoi_vol               , & ! Input:  [real(r8) (:,:)  ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevgrnd)
        h2osoi_liq              => col_ws%h2osoi_liq               , & ! Input:  [real(r8) (:,:)  ]  liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)

        t_soisno                => col_es%t_soisno                , & ! Input:  [real(r8) (:,:)  ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)

        o2_decomp_depth_unsat   => ch4_vars%o2_decomp_depth_unsat_col           , & ! Input:  [real(r8) (:,:)  ]  O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
        o2_decomp_depth_sat     => ch4_vars%o2_decomp_depth_sat_col             , & ! Input:  [real(r8) (:,:)  ]  O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
        conc_o2_unsat           => ch4_vars%conc_o2_unsat_col                   , & ! Input:  [real(r8) (:,:)  ]  O2 conc in each soil layer (mol/m3) (nlevsoi)
        conc_o2_sat             => ch4_vars%conc_o2_sat_col                     , & ! Input:  [real(r8) (:,:)  ]  O2 conc in each soil layer (mol/m3) (nlevsoi)
        o2stress_unsat          => ch4_vars%o2stress_unsat_col                  , & ! Input:  [real(r8) (:,:)   ]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
        o2stress_sat            => ch4_vars%o2stress_sat_col                    , & ! Input:  [real(r8) (:,:)   ]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)

        finundated              => ch4_vars%finundated_col                        & ! Input:  [real(r8) (:)     ]  fractional inundated area (excluding dedicated wetland columns)
    )

    ! soil properties & thermohydrology

    do fc = 1, num_soilc
        c = filter_soilc(fc)

        plant_ndemand_col(c)        = clm_interface_data%bgc%plant_ndemand_col(c)
        plant_pdemand_col(c)        = clm_interface_data%bgc%plant_pdemand_col(c)

        finundated(c)               = clm_interface_data%bgc%finundated_col(c)

        bd(c,:)                     = clm_interface_data%bd_col(c,:)
        watsat(c,:)                 = clm_interface_data%watsat_col(c,:)
        bsw(c,:)                    = clm_interface_data%bsw_col(c,:)
        sucsat(c,:)                 = clm_interface_data%sucsat_col(c,:)
        watfc(c,:)                  = clm_interface_data%watfc_col(c,:)
        cellorg(c,:)                = clm_interface_data%cellorg_col(c,:)

        soilpsi(c,:)                = clm_interface_data%th%soilpsi_col(c,:)
        h2osoi_vol(c,:)             = clm_interface_data%th%h2osoi_vol_col(c,:)
        h2osoi_liq(c,:)             = clm_interface_data%th%h2osoi_liq_col(c,:)

        t_soisno(c,:)               = clm_interface_data%th%t_soisno_col(c,:)

        o2stress_unsat(c,:)         = clm_interface_data%bgc%o2stress_unsat_col(c,:)
        o2stress_sat(c,:)           = clm_interface_data%bgc%o2stress_sat_col(c,:)
        o2_decomp_depth_unsat(c,:)  = clm_interface_data%bgc%o2_decomp_depth_unsat_col(c,:)
        conc_o2_unsat(c,:)          = clm_interface_data%bgc%conc_o2_unsat_col(c,:)
        o2_decomp_depth_sat(c,:)    = clm_interface_data%bgc%o2_decomp_depth_sat_col(c,:)
        conc_o2_sat(c,:)            = clm_interface_data%bgc%conc_o2_sat_col(c,:)

    end do

    !state variables
    do fc = 1, num_soilc
        c = filter_soilc(fc)
            do k = 1, ndecomp_pools
                decomp_cpools_vr(c,:,k) = clm_interface_data%bgc%decomp_cpools_vr_col(c,:,k)
                decomp_npools_vr(c,:,k) = clm_interface_data%bgc%decomp_npools_vr_col(c,:,k)
                decomp_ppools_vr(c,:,k) = clm_interface_data%bgc%decomp_ppools_vr_col(c,:,k)
            end do

            smin_no3_vr(c,:)        = clm_interface_data%bgc%smin_no3_vr_col(c,:)
            smin_nh4_vr(c,:)        = clm_interface_data%bgc%smin_nh4_vr_col(c,:)
            smin_nh4sorb_vr(c,:)    = clm_interface_data%bgc%smin_nh4sorb_vr_col(c,:)

            solutionp_vr(c,:)       = clm_interface_data%bgc%solutionp_vr_col(c,:)
            labilep_vr(c,:)         = clm_interface_data%bgc%labilep_vr_col(c,:)
            secondp_vr(c,:)         = clm_interface_data%bgc%secondp_vr_col(c,:)
            sminp_vr(c,:)           = clm_interface_data%bgc%sminp_vr_col(c,:)
            occlp_vr(c,:)           = clm_interface_data%bgc%occlp_vr_col(c,:)
            primp_vr(c,:)           = clm_interface_data%bgc%primp_vr_col(c,:)
    end do

    end associate
  end subroutine clm_bgc_get_data
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
  ! !INTERFACE:
  ! pass data from clm_bgc to clm_bgc_data
  subroutine clm_bgc_update_data(clm_bgc_data, bounds,  &
            num_soilc, filter_soilc,                    &
            cnstate_vars, carbonflux_vars,              &
            nitrogenflux_vars, phosphorusflux_vars)

    ! USES:

    ! ARGUMENTS:
    type(bounds_type)                   , intent(in)    :: bounds
    integer                             , intent(in)    :: num_soilc          ! number of soil columns in filter
    integer                             , intent(in)    :: filter_soilc(:)    ! filter for soil columns

    type(cnstate_type)                  , intent(in)    :: cnstate_vars
    type(carbonflux_type)               , intent(in)    :: carbonflux_vars
    type(nitrogenflux_type)             , intent(in)    :: nitrogenflux_vars
    type(phosphorusflux_type)           , intent(in)    :: phosphorusflux_vars

    type(elm_interface_bgc_datatype)    , intent(inout) :: clm_bgc_data

    ! LOCAL VARIABLES:
    integer :: fc, c, j, k

    !-----------------------------------------------------------------------

    associate(                                                                                            &
         fpg                          => cnstate_vars%fpg_col                                           , & ! Output: [real(r8) (:)   ]  fraction of potential gpp (no units)
         fpi                          => cnstate_vars%fpi_col                                           , & ! Output: [real(r8) (:)   ]  fraction of potential immobilization (no units)
         fpi_vr                       => cnstate_vars%fpi_vr_col                                        , & ! Output: [real(r8) (:,:) ]  fraction of potential immobilization (no units)
         fpg_p                        => cnstate_vars%fpg_p_col                                         , & ! Output: [real(r8) (:)   ]  fraction of potential gpp (no units)
         fpi_p                        => cnstate_vars%fpi_p_col                                         , & ! Output: [real(r8) (:)   ]  fraction of potential immobilization (no units)
         fpi_p_vr                     => cnstate_vars%fpi_p_vr_col                                      , & ! Output: [real(r8) (:,:) ]  fraction of potential immobilization (no units)

         potential_immob              => col_nf%potential_immob                          , & ! Output: [real(r8) (:)   ]
         actual_immob                 => col_nf%actual_immob                             , & ! Output: [real(r8) (:)   ]
         sminn_to_plant               => col_nf%sminn_to_plant                           , & ! Output: [real(r8) (:)   ]
         sminn_to_denit_excess_vr     => col_nf%sminn_to_denit_excess_vr                 , & ! Output: [real(r8) (:,:) ]
         pot_f_nit_vr                 => col_nf%pot_f_nit_vr                             , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) potential soil nitrification flux
         pot_f_denit_vr               => col_nf%pot_f_denit_vr                           , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) potential soil denitrification flux
         f_nit_vr                     => col_nf%f_nit_vr                                 , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) soil nitrification flux
         f_denit_vr                   => col_nf%f_denit_vr                               , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) soil denitrification flux
         actual_immob_no3_vr          => col_nf%actual_immob_no3_vr                      , & ! Output: [real(r8) (:,:) ]
         actual_immob_nh4_vr          => col_nf%actual_immob_nh4_vr                      , & ! Output: [real(r8) (:,:) ]
         smin_no3_to_plant_vr         => col_nf%smin_no3_to_plant_vr                     , & ! Output: [real(r8) (:,:) ]
         smin_nh4_to_plant_vr         => col_nf%smin_nh4_to_plant_vr                     , & ! Output: [real(r8) (:,:) ]
         n2_n2o_ratio_denit_vr        => col_nf%n2_n2o_ratio_denit_vr                    , & ! Output: [real(r8) (:,:) ]  ratio of N2 to N2O production by denitrification [gN/gN]
         f_n2o_denit_vr               => col_nf%f_n2o_denit_vr                           , & ! Output: [real(r8) (:,:) ]  flux of N2O from denitrification [gN/m3/s]
         f_n2o_nit_vr                 => col_nf%f_n2o_nit_vr                             , & ! Output: [real(r8) (:,:) ]  flux of N2O from nitrification [gN/m3/s]
         supplement_to_sminn_vr       => col_nf%supplement_to_sminn_vr                   , & ! Output: [real(r8) (:,:) ]
         sminn_to_plant_vr            => col_nf%sminn_to_plant_vr                        , & ! Output: [real(r8) (:,:) ]
         actual_immob_vr              => col_nf%actual_immob_vr                          , & ! Output: [real(r8) (:,:) ]

         potential_immob_p            => col_pf%potential_immob_p                      , & ! Output: [real(r8) (:)   ]
         actual_immob_p               => col_pf%actual_immob_p                         , & ! Output: [real(r8) (:)   ]
         sminp_to_plant               => col_pf%sminp_to_plant                         , & ! Output: [real(r8) (:)   ]
         supplement_to_sminp_vr       => col_pf%supplement_to_sminp_vr                 , & ! Output: [real(r8) (:,:) ]
         sminp_to_plant_vr            => col_pf%sminp_to_plant_vr                      , & ! Output: [real(r8) (:,:) ]
         actual_immob_p_vr            => col_pf%actual_immob_p_vr                      , & ! Output: [real(r8) (:,:) ]

         decomp_cascade_ntransfer_vr      =>    col_nf%decomp_cascade_ntransfer_vr       , & ! Output: [real(r8) (:,:,:) ]  vert-res transfer of N from donor to receiver pool along decomp. cascade (gN/m3/s)
         decomp_cascade_sminn_flux_vr     =>    col_nf%decomp_cascade_sminn_flux_vr      , & ! Output: [real(r8) (:,:,:) ]  vert-res mineral N flux for transition along decomposition cascade (gN/m3/s)
         potential_immob_vr               =>    col_nf%potential_immob_vr                , & ! Output: [real(r8) (:,:)   ]
         sminn_to_denit_decomp_cascade_vr =>    col_nf%sminn_to_denit_decomp_cascade_vr  , & ! Output: [real(r8) (:,:,:) ]
         gross_nmin_vr                    =>    col_nf%gross_nmin_vr                     , & ! Output: [real(r8) (:,:)   ]
         net_nmin_vr                      =>    col_nf%net_nmin_vr                       , & ! Output: [real(r8) (:,:)   ]
         gross_nmin                       =>    col_nf%gross_nmin                        , & ! Output: [real(r8) (:)     ]  gross rate of N mineralization (gN/m2/s)
         net_nmin                         =>    col_nf%net_nmin                          , & ! Output: [real(r8) (:)     ]  net rate of N mineralization (gN/m2/s)

         ! add phosphorus
         decomp_cascade_ptransfer_vr      =>    col_pf%decomp_cascade_ptransfer_vr     , & ! Output: [real(r8) (:,:,:) ]  vert-res transfer of P from donor to receiver pool along decomp. cascade (gP/m3/s)
         decomp_cascade_sminp_flux_vr     =>    col_pf%decomp_cascade_sminp_flux_vr    , & ! Output: [real(r8) (:,:,:) ]  vert-res mineral P flux for transition along decomposition cascade (gP/m3/s)
         potential_immob_p_vr             =>    col_pf%potential_immob_p_vr            , & ! Output: [real(r8) (:,:)   ]
         gross_pmin_vr                    =>    col_pf%gross_pmin_vr                   , & ! Output: [real(r8) (:,:)   ]
         net_pmin_vr                      =>    col_pf%net_pmin_vr                     , & ! Output: [real(r8) (:,:)   ]
         gross_pmin                       =>    col_pf%gross_pmin                      , & ! Output: [real(r8) (:)     ]  gross rate of P mineralization (gP/m2/s)
         net_pmin                         =>    col_pf%net_pmin                        , & ! Output: [real(r8) (:)     ]  net rate of P mineralization (gP/m2/s)

         decomp_cascade_hr_vr             =>    col_cf%decomp_cascade_hr_vr                , & ! Output: [real(r8) (:,:,:) ]  vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
         decomp_cascade_ctransfer_vr      =>    col_cf%decomp_cascade_ctransfer_vr         , & ! Output: [real(r8) (:,:,:) ]  vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
         phr_vr                           =>    col_cf%phr_vr                              , & ! Output: [real(r8) (:,:)   ]  potential HR (gC/m3/s)
         fphr                             =>    col_cf%fphr                                  & ! Output: [real(r8) (:,:)   ]  fraction of potential SOM + LITTER heterotrophic
         )

    !---------------------------------------------------------------------------
        do fc = 1,num_soilc
            c = filter_soilc(fc)
            clm_bgc_data%fpg_col(c)                                   = fpg(c)
            clm_bgc_data%fpi_col(c)                                   = fpi(c)
            clm_bgc_data%fpi_vr_col(c,:)                              = fpi_vr(c,:)
            clm_bgc_data%fpg_p_col(c)                                 = fpg_p(c)
            clm_bgc_data%fpi_p_col(c)                                 = fpi_p(c)
            clm_bgc_data%fpi_p_vr_col(c,:)                            = fpi_p_vr(c,:)
            clm_bgc_data%potential_immob_col(c)                       = potential_immob(c)
            clm_bgc_data%actual_immob_col(c)                          = actual_immob(c)
            clm_bgc_data%sminn_to_plant_col(c)                        = sminn_to_plant(c)
            clm_bgc_data%sminn_to_denit_excess_vr_col(c,:)            = sminn_to_denit_excess_vr(c,:)
            clm_bgc_data%pot_f_nit_vr_col(c,:)                        = pot_f_nit_vr(c,:)
            clm_bgc_data%pot_f_denit_vr_col(c,:)                      = pot_f_denit_vr(c,:)
            clm_bgc_data%f_nit_vr_col(c,:)                            = f_nit_vr(c,:)
            clm_bgc_data%f_denit_vr_col(c,:)                          = f_denit_vr(c,:)
            clm_bgc_data%actual_immob_no3_vr_col(c,:)                 = actual_immob_no3_vr(c,:)
            clm_bgc_data%actual_immob_nh4_vr_col(c,:)                 = actual_immob_nh4_vr(c,:)
            clm_bgc_data%smin_no3_to_plant_vr_col(c,:)                = smin_no3_to_plant_vr(c,:)
            clm_bgc_data%smin_nh4_to_plant_vr_col(c,:)                = smin_nh4_to_plant_vr(c,:)
            clm_bgc_data%n2_n2o_ratio_denit_vr_col(c,:)               = n2_n2o_ratio_denit_vr(c,:)
            clm_bgc_data%f_n2o_denit_vr_col(c,:)                      = f_n2o_denit_vr(c,:)
            clm_bgc_data%f_n2o_nit_vr_col(c,:)                        = f_n2o_nit_vr(c,:)
            clm_bgc_data%supplement_to_sminn_vr_col(c,:)              = supplement_to_sminn_vr(c,:)
            clm_bgc_data%sminn_to_plant_vr_col(c,:)                   = sminn_to_plant_vr(c,:)
            clm_bgc_data%potential_immob_vr_col(c,:)                  = potential_immob_vr(c,:)
            clm_bgc_data%actual_immob_vr_col(c,:)                     = actual_immob_vr(c,:)
            clm_bgc_data%potential_immob_p_col(c)                     = potential_immob_p(c)
            clm_bgc_data%actual_immob_p_col(c)                        = actual_immob_p(c)
            clm_bgc_data%sminp_to_plant_col(c)                        = sminp_to_plant(c)
            clm_bgc_data%supplement_to_sminp_vr_col(c,:)              = supplement_to_sminp_vr(c,:)
            clm_bgc_data%sminp_to_plant_vr_col(c,:)                   = sminp_to_plant_vr(c,:)
            clm_bgc_data%potential_immob_p_vr_col(c,:)                = potential_immob_p_vr(c,:)
            clm_bgc_data%actual_immob_p_vr_col(c,:)                   = actual_immob_p_vr(c,:)

            clm_bgc_data%decomp_cascade_ntransfer_vr_col(c,:,:)       = decomp_cascade_ntransfer_vr(c,:,:)
            clm_bgc_data%decomp_cascade_sminn_flux_vr_col(c,:,:)      = decomp_cascade_sminn_flux_vr(c,:,:)
            clm_bgc_data%sminn_to_denit_decomp_cascade_vr_col(c,:,:)  = sminn_to_denit_decomp_cascade_vr(c,:,:)
            clm_bgc_data%gross_nmin_vr_col(c,:)                       = gross_nmin_vr(c,:)
            clm_bgc_data%net_nmin_vr_col(c,:)                         = net_nmin_vr(c,:)

            ! phosphorus
            clm_bgc_data%decomp_cascade_ptransfer_vr_col(c,:,:)       = decomp_cascade_ptransfer_vr(c,:,:)
            clm_bgc_data%decomp_cascade_sminp_flux_vr_col(c,:,:)      = decomp_cascade_sminp_flux_vr(c,:,:)
            clm_bgc_data%potential_immob_p_vr_col(c,:)                = potential_immob_p_vr(c,:)
            clm_bgc_data%gross_pmin_vr_col(c,:)                       = gross_pmin_vr(c,:)
            clm_bgc_data%net_pmin_vr_col(c,:)                         = net_pmin_vr(c,:)

            clm_bgc_data%decomp_cascade_hr_vr_col(c,:,:)              = decomp_cascade_hr_vr(c,:,:)
            clm_bgc_data%decomp_cascade_ctransfer_vr_col(c,:,:)       = decomp_cascade_ctransfer_vr(c,:,:)

            clm_bgc_data%phr_vr_col(c,:)                              = phr_vr(c,:)
            clm_bgc_data%fphr_col(c,:)                                = fphr(c,:)
        end do


    end associate
  end subroutine clm_bgc_update_data
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
  subroutine update_bgc_data_clm2clm(clm_bgc_data,bounds,         &
           num_soilc, filter_soilc,                               &
           num_soilp, filter_soilp,                               &
           cnstate_vars, carbonflux_vars, carbonstate_vars,       &
           nitrogenflux_vars, nitrogenstate_vars,                 &
           phosphorusflux_vars, phosphorusstate_vars,             &
           ch4_vars)
    ! USES


    implicit none

    ! !ARGUMENTS:
    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                     , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                     , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                     , intent(in)    :: filter_soilp(:)   ! filter for soil patches

    type(cnstate_type)          , intent(inout) :: cnstate_vars
    type(carbonflux_type)       , intent(inout) :: carbonflux_vars
    type(carbonstate_type)      , intent(inout) :: carbonstate_vars
    type(nitrogenflux_type)     , intent(inout) :: nitrogenflux_vars
    type(nitrogenstate_type)    , intent(inout) :: nitrogenstate_vars
    type(phosphorusflux_type)   , intent(inout) :: phosphorusflux_vars
    type(phosphorusstate_type)  , intent(inout) :: phosphorusstate_vars
    type(ch4_type)              , intent(inout) :: ch4_vars

    type(elm_interface_bgc_datatype), intent(in):: clm_bgc_data

    !-----------------------------------------------------------------------

    character(len=256) :: subname = "update_bgc_data_clm2clm"

    ! bgc_state_decomp is updated in CLM
    ! by passing bgc_flux_decomp_sourcesink into SoilLittVertTransp
    call update_bgc_flux_decomp_cascade(clm_bgc_data,   &
                bounds, num_soilc, filter_soilc,        &
                carbonflux_vars, nitrogenflux_vars,     &
                phosphorusflux_vars)

    call update_bgc_flux_smin(clm_bgc_data,             &
                bounds, num_soilc, filter_soilc,        &
                cnstate_vars,                           &
                nitrogenflux_vars, phosphorusflux_vars)

    call update_bgc_flux_nitdenit(clm_bgc_data,         &
                bounds, num_soilc, filter_soilc,        &
                nitrogenflux_vars, phosphorusflux_vars)

  end subroutine update_bgc_data_clm2clm
!--------------------------------------------------------------------------------------
! END of CLM-bgc through interface
!--------------------------------------------------------------------------------------


end module clm_interface_funcsMod

