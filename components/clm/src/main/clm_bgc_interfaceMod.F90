module clm_bgc_interfaceMod
!!=================================================================================================
! CLM BioGeoChemistry (BGC) Interface
!
! Created by wgs @ ORNL
!
! date: 8/25/2015
!!=================================================================================================

#include "shr_assert.h"


  !! MODULE: clm_bgc_interfaceMod
  !!--------------------------------------------------------------------------------------
  !! DESCRIPTION:
  !! Coupling of CLM with any specific Soil BGC module Consists of 3 STEPS:
  !! STEP-1:   clm vars         -> clm_bgc_data (i.e. clm_bgc_interface_data_type)  ; pass clm vars to clm_bgc_data
  !! STEP-2:   clm_bgc_data     -> soil bgc module -> clm_bgc_data
  !!      2.1: clm_bgc_data     -> soil bgc module
  !!      2.2: run soil bgc module
  !!      2.3: soil bgc module  -> clm_bgc_data
  !! STEP-3:   clm_bgc_data     -> clm vars
  !!--------------------------------------------------------------------------------------


  !!--------------------------------------------------------------------------------------
  ! !USES:
  ! clm g/l/c/p constants
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use GridcellType          , only : grc_pp
  use LandunitType          , only : lun_pp
  use ColumnType            , only : col_pp 
  use PatchType             , only : pft_pp

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

  use ch4Mod                , only : ch4_type

  use PhotosynthesisType    , only : photosyns_type
  use cropType              , only : crop_type
  use CanopyStateType       , only : canopystate_type
  use PhosphorusStateType   , only : phosphorusstate_type
  use PhosphorusFluxType    , only : phosphorusflux_type

  use SoilWaterRetentionCurveMod    , only : soil_water_retention_curve_type
  use clm_bgc_interface_data        , only : clm_bgc_interface_data_type

  ! most used constants in this module
  use clm_varpar            , only : nlevsoi, nlevsno,nlevgrnd, nlevdecomp, nlevdecomp_full
  use clm_varpar            , only : ndecomp_pools
  use clm_varpar            , only : max_patch_per_col
  use clm_varcon            , only : denh2o, denice, tfrz, dzsoi_decomp
  use landunit_varcon       , only : istsoil, istcrop

  ! misc.
  use abortutils            , only : endrun
  !!--------------------------------------------------------------------------------------

  implicit none

  save

  private    ! By default everything is private

  !! LOCAL VARIABLES:


  !!--------------------------------------------------------------------------------------
  !! (1) GENERIC SUBROUTINES: used by any specific soil BGC module
  !! pass clm variables to clm_bgc_data
  public    :: get_clm_bgc_data             !! STEP-1: clm vars -> clm_bgc_data

  !! pass clm variables to clm_bgc_data, called by get_clm_bgc_data
  private   :: get_clm_soil_property        !! STEP-1.1: soil properties
  private   :: get_clm_soil_thermohydro     !! STEP-1.2: thermohydrology vars
  private   :: get_clm_bgc_state            !! STEP-1.3: state vars
  private   :: get_clm_bgc_flux             !! STEP-1.4: flux vars

  !! STEP-3.x: clm_bgc_data -> clm vars
  !! update clm variables from clm_bgc_data,
  !! e.g., called in 'update_bgc_data_clm2clm' and 'update_bgc_data_pf2clm'
  !! specific bgc-module (e.g., PFLOTRAN) requires certain combination of these subroutines
  private   :: update_bgc_state_decomp
  private   :: update_bgc_state_smin
  private   :: update_bgc_flux_decomp_sourcesink
  private   :: update_bgc_flux_decomp_cascade
  private   :: update_bgc_flux_smin
  private   :: update_bgc_flux_nitdenit
  private   :: update_bgc_flux_gas_pf
  private   :: update_soil_moisture
  private   :: update_soil_temperature

  !!--------------------------------------------------------------------------------------
  !! (2) SPECIFIC SUBROUTINES: used by a specific soil BGC module
  !! (2.1) Specific Subroutines for running clm-bgc (CN or BGC) through interface
  !! if (use_bgc_interface .and. use_clm_bgc)
  public    :: clm_bgc_run              !! STEP-2:   clm_bgc_data  -> clm-bgc module -> clm_bgc_data    ; called in clm_driver
  private   :: clm_bgc_get_data         !! STEP-2.1: clm_bgc_data  -> clm-bgc module                    ; called in clm_bgc_run
                                        !! STEP-2.2: run clm-bgc module                                 ; see CNDecompAlloc in CNDecompMod
  private   :: clm_bgc_update_data      !! STEP-2.3: clm-bgc module-> clm_bgc_data                      ; called in clm_bgc_run
  public    :: update_bgc_data_clm2clm  !! STEP-3:   clm_bgc_data  -> clm vars                          ; called in clm_driver

  !! (2.2) Specific Subroutines for CLM-PFLOTRAN Coupling: update clm variables from pflotran
  !! if (use_bgc_interface .and. use_pflotran)
  public    :: update_bgc_data_pf2clm   !! STEP-3:   clm_bgc_data  -> clm vars                          ; called in clm_driver
                                        !! STEP-2:   see 'clm_pf_run' in clm_pflotran_interfaceMod
  !!--------------------------------------------------------------------------------------

contains


!!--------------------------------------------------------------------------------------
  subroutine get_clm_bgc_data(clm_bgc_data,bounds,                &
           num_soilc, filter_soilc,                               &
           num_soilp, filter_soilp,                               &
           atm2lnd_vars,                                          &
           waterstate_vars, waterflux_vars,                       &
           soilstate_vars,  temperature_vars, energyflux_vars,    &
           soilhydrology_vars, soil_water_retention_curve,        &
           cnstate_vars, carbonflux_vars, carbonstate_vars,       &
           nitrogenflux_vars, nitrogenstate_vars,                 &
           phosphorusflux_vars, phosphorusstate_vars,             &
           canopystate_vars, ch4_vars                             &
           )

    implicit none

    ! !ARGUMENTS:
    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                     , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                     , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                     , intent(in)    :: filter_soilp(:)   ! filter for soil patches
    type(atm2lnd_type)          , intent(in)    :: atm2lnd_vars
    type(canopystate_type)      , intent(in)    :: canopystate_vars
    type(waterstate_type)       , intent(in)    :: waterstate_vars
    type(waterflux_type)        , intent(in)    :: waterflux_vars
    type(soilstate_type)        , intent(in)    :: soilstate_vars
    type(temperature_type)      , intent(in)    :: temperature_vars
    type(soilhydrology_type)    , intent(in)    :: soilhydrology_vars
    type(energyflux_type)       , intent(in)    :: energyflux_vars

    type(cnstate_type)          , intent(in)    :: cnstate_vars
    type(carbonflux_type)       , intent(in)    :: carbonflux_vars
    type(carbonstate_type)      , intent(in)    :: carbonstate_vars
    type(nitrogenflux_type)     , intent(in)    :: nitrogenflux_vars
    type(nitrogenstate_type)    , intent(in)    :: nitrogenstate_vars
    type(phosphorusflux_type)   , intent(in)    :: phosphorusflux_vars
    type(phosphorusstate_type)  , intent(in)    :: phosphorusstate_vars
    type(ch4_type)              , intent(in)    :: ch4_vars

    class(soil_water_retention_curve_type)  , intent(in)    :: soil_water_retention_curve
    type(clm_bgc_interface_data_type)       , intent(inout) :: clm_bgc_data

    !-----------------------------------------------------------------------

    character(len=256) :: subname = "get_clm_bgc_data"

    call get_clm_soil_property(clm_bgc_data,                &
                    bounds, num_soilc, filter_soilc,        &
                    soilstate_vars)

    call get_clm_soil_thermohydro(clm_bgc_data,             &
                   bounds, num_soilc, filter_soilc,         &
                   atm2lnd_vars, soilstate_vars,            &
                   waterstate_vars, waterflux_vars,         &
                   temperature_vars, energyflux_vars,       &
                   soil_water_retention_curve,              &
                   canopystate_vars, ch4_vars)

    call get_clm_bgc_state(clm_bgc_data,                    &
                    bounds, num_soilc, filter_soilc,        &
                    carbonstate_vars, nitrogenstate_vars,   &
                    phosphorusstate_vars)

    call get_clm_bgc_flux(clm_bgc_data,                     &
                    bounds, num_soilc, filter_soilc,        &
                    cnstate_vars, carbonflux_vars,          &
                    nitrogenflux_vars, phosphorusflux_vars)

  end subroutine get_clm_bgc_data
!!--------------------------------------------------------------------------------------

!!--------------------------------------------------------------------------------------
  subroutine get_clm_soil_property(clm_bgc_data,            &
                        bounds, num_soilc, filter_soilc,    &
                        soilstate_vars)
    !
    ! !DESCRIPTION:
    ! get soil column physical properties
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

    type(clm_bgc_interface_data_type), intent(inout) :: clm_bgc_data

    integer  :: fc, g, l, c, j      ! indices
    integer  :: gcount, cellcount

    character(len= 32) :: subname = 'get_clm_soil_property' ! subroutine name

    associate( &
         ! Assign local pointer to derived subtypes components (column-level)
         z                  => col_pp%z                                                , & !  [real(r8) (:,:)]  layer depth (m)
         dz                 => col_pp%dz                                               , & !  [real(r8) (:,:)]  layer thickness depth (m)

         bd                 => soilstate_vars%bd_col                                , & !
         bsw                => soilstate_vars%bsw_col                               , & !  [real(r8) (:,:)]  Clapp and Hornberger "b" (nlevgrnd)
         hksat              => soilstate_vars%hksat_col                             , & !  [real(r8) (:,:)]  hydraulic conductivity at saturation (mm H2O /s) (nlevgrnd)
         sucsat             => soilstate_vars%sucsat_col                            , & !  [real(r8) (:,:)]  minimum soil suction (mm) (nlevgrnd)
         watsat             => soilstate_vars%watsat_col                            , & !  [real(r8) (:,:)]  volumetric soil water at saturation (porosity) (nlevgrnd)
         watfc              => soilstate_vars%watfc_col                             , & !  [real(r8) (:,:)]  volumetric soil water at saturation (porosity) (nlevgrnd)

         cellorg            => soilstate_vars%cellorg_col                           , & ! Input:  [real(r8) (:,:)  ]  column 3D org (kg/m3 organic matter) (nlevgrnd)

         porosity           => soilstate_vars%porosity_col                          , &
         eff_porosity       => soilstate_vars%eff_porosity_col                      , &

         initial_cn_ratio   => decomp_cascade_con%initial_cn_ratio                  , &
         initial_cp_ratio   => decomp_cascade_con%initial_cp_ratio                  , &

         decomp_pool_name   => decomp_cascade_con%decomp_pool_name_history          , &
         floating_cn_ratio  => decomp_cascade_con%floating_cn_ratio_decomp_pools    , &
         floating_cp_ratio  => decomp_cascade_con%floating_cp_ratio_decomp_pools      &

         )

!-------------------------------------------------------------------------------------
    !! constants:
    clm_bgc_data%ndecomp_pools          = ndecomp_pools
    clm_bgc_data%decomp_pool_name(:)    = decomp_pool_name(:)
    clm_bgc_data%floating_cn_ratio(:)   = floating_cn_ratio(:)
    clm_bgc_data%floating_cp_ratio(:)   = floating_cp_ratio(:)

    clm_bgc_data%initial_cn_ratio(:)    = initial_cn_ratio(:)
    clm_bgc_data%initial_cp_ratio(:)    = initial_cp_ratio(:)



    do fc = 1, num_soilc
        c = filter_soilc(fc)
!        do j = 1,nlevsoi
            clm_bgc_data%z(c,:)                 = z(c,:)
            clm_bgc_data%dz(c,:)                = dz(c,:)
            clm_bgc_data%bd_col(c,:)            = bd(c,:)
            clm_bgc_data%bsw_col(c,:)           = bsw(c,:)
            clm_bgc_data%hksat_col(c,:)         = hksat(c,:)
            clm_bgc_data%sucsat_col(c,:)        = sucsat(c,:)
            clm_bgc_data%watsat_col(c,:)        = watsat(c,:)
            clm_bgc_data%watfc_col(c,:)         = watfc(c,:)

            clm_bgc_data%porosity_col(c,:)      = porosity(c,:)
            clm_bgc_data%eff_porosity_col(c,:)  = eff_porosity(c,:)

            clm_bgc_data%cellorg_col(c,:)       = cellorg(c,:)
!        end do
    end do

    end associate
  end subroutine get_clm_soil_property
!!--------------------------------------------------------------------------------------

!!--------------------------------------------------------------------------------------
  subroutine get_clm_soil_thermohydro(clm_bgc_data,         &
                       bounds, num_soilc, filter_soilc,     &
                       atm2lnd_vars, soilstate_vars,        &
                       waterstate_vars, waterflux_vars,     &
                       temperature_vars, energyflux_vars,   &
                       soil_water_retention_curve,          &
                       canopystate_vars, ch4_vars)
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
    type(waterflux_type)     , intent(in) :: waterflux_vars
    type(temperature_type)   , intent(in) :: temperature_vars
    type(energyflux_type)    , intent(in) :: energyflux_vars
    type(canopystate_type)   , intent(in) :: canopystate_vars
    type(ch4_type)           , intent(in) :: ch4_vars

    class(soil_water_retention_curve_type)  , intent(in)    :: soil_water_retention_curve
    type(clm_bgc_interface_data_type)       , intent(inout) :: clm_bgc_data

  ! !LOCAL VARIABLES:
    integer  :: fc, c, j      ! indices
    integer  :: pftindex, p
!    real(r8) :: sattmp, psitmp, itheta
!    real(r8) :: watmin(num_soilc, nlevsoi)
!    real(r8) :: sucmin(num_soilc, nlevsoi)

  !EOP
  !-----------------------------------------------------------------------
    associate ( &
      gridcell              => col_pp%gridcell                             , & ! column's gridcell
      wtgcell               => col_pp%wtgcell                              , & ! column's weight relative to gridcell
      cactive               => col_pp%active                               , & ! [logical (:)]  column active or not
      dz                    => col_pp%dz                                   , & ! layer thickness depth (m)
      zi                    => col_pp%zi                                   , & ! interface depth (m)

      soilpsi               => soilstate_vars%soilpsi_col               , & ! soil water matric potential in each soil layer (MPa)
      rootfr                => soilstate_vars%rootfr_col                , & ! pft-level effective fraction of roots in each soil layer

      h2osoi_vol            => waterstate_vars%h2osoi_vol_col           , & ! Input:  [real(r8) (:,:)  ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevgrnd)
      h2osoi_liq            => waterstate_vars%h2osoi_liq_col           , & ! liquid water (kg/m2)
      h2osoi_ice            => waterstate_vars%h2osoi_ice_col           , & ! ice lens (kg/m2)
      frac_sno              => waterstate_vars%frac_sno_eff_col         , & ! Input: fraction of ground covered by snow (0 to 1)
      frac_h2osfc           => waterstate_vars%frac_h2osfc_col          , & ! Input: fraction of ground covered by surface water (0 to 1)

      t_soisno              => temperature_vars%t_soisno_col            , & ! snow-soil temperature (Kelvin)
      t_grnd                => temperature_vars%t_grnd_col              , & ! Input:  [real(r8) (:)]  ground surface temperature [K]

      forc_pbot             => atm2lnd_vars%forc_pbot_not_downscaled_grc, &  ! atmospheric pressure (Pa)
      forc_pco2             => atm2lnd_vars%forc_pco2_grc               , & ! partial pressure co2 (Pa)
      forc_pch4             => atm2lnd_vars%forc_pch4_grc               , & ! partial pressure ch4 (Pa)

      alt_indx              => canopystate_vars%alt_indx_col            , & ! Input:  [integer  (:)     ]  current depth of thaw
      o2stress_sat          => ch4_vars%o2stress_sat_col                , & ! Input:  [real(r8) (:,:)   ]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
      o2stress_unsat        => ch4_vars%o2stress_unsat_col              , & ! Input:  [real(r8) (:,:)   ]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
      finundated            => ch4_vars%finundated_col                  , & ! Input:  [real(r8) (:)     ]  fractional inundated area (excluding dedicated wetland columns)
      o2_decomp_depth_unsat => ch4_vars%o2_decomp_depth_unsat_col       , & ! Input:  [real(r8) (:,:)  ]  O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
      conc_o2_unsat         => ch4_vars%conc_o2_unsat_col               , & ! Input:  [real(r8) (:,:)  ]  O2 conc in each soil layer (mol/m3) (nlevsoi)
      o2_decomp_depth_sat   => ch4_vars%o2_decomp_depth_sat_col         , & ! Input:  [real(r8) (:,:)  ]  O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
      conc_o2_sat           => ch4_vars%conc_o2_sat_col                 , & ! Input:  [real(r8) (:,:)  ]  O2 conc in each soil layer (mol/m3) (nlevsoi)

      htvp                  => energyflux_vars%htvp_col                 , & ! Input:  [real(r8) (:)]  latent heat of vapor of water (or sublimation) [j/kg]
      eflx_bot              => energyflux_vars%eflx_bot_col             , & ! heat flux from beneath column (W/m**2) [+ = upward]
      eflx_gnet_patch       => energyflux_vars%eflx_gnet_patch          , & ! net ground heat flux into the surface (W/m**2) per patch
      eflx_soil_grnd_patch  => energyflux_vars%eflx_soil_grnd_patch     , & ! soil heat flux (W/m**2) [+ = into soil]

      qflx_top_soil         => waterflux_vars%qflx_top_soil_col         , & ! Input: net water input into soil from top (mm/s)
      qflx_ev_h2osfc        => waterflux_vars%qflx_ev_h2osfc_col        , & ! Input: column-level evaporation flux from h2osfc (W/m2) [+ to atm] : checking unit
      qflx_evap_soi         => waterflux_vars%qflx_evap_soi_col         , & ! Input: column-level soil evaporation (mm H2O/s) (+ = to atm)
      qflx_sub_snow         => waterflux_vars%qflx_sub_snow_col         , & ! Input: column-level evaporation flux from snow (mm H2O/s) [+ to atm]
      qflx_tran_veg         => waterflux_vars%qflx_tran_veg_col           & ! Input: pft-level vegetation transpiration (mm H2O/s) (+ = to atm)

    )

    !--------------------------------------------------------------------------------------
!
!    watmin(:,:) = 0.01_r8
!    sucmin(:,:) = 1.e8_r8

    !! grid:
    clm_bgc_data%forc_pbot_not_downscaled_grc(:)    = forc_pbot(:)
    clm_bgc_data%forc_pco2_grc(:)                   = forc_pco2(:)
    clm_bgc_data%forc_pch4_grc(:)                   = forc_pch4(:)


    do fc = 1,num_soilc
        c = filter_soilc(fc)

        clm_bgc_data%frac_sno_eff_col(c)    = frac_sno(c)
        clm_bgc_data%frac_h2osfc_col(c)     = frac_h2osfc(c)

        clm_bgc_data%t_grnd_col(c)          = t_grnd(c)

        clm_bgc_data%alt_indx_col(c)        = alt_indx(c)
        clm_bgc_data%finundated_col(c)      = finundated(c)

        clm_bgc_data%qflx_top_soil_col(c)   = qflx_top_soil(c)
        clm_bgc_data%qflx_ev_h2osfc_col(c)  = qflx_ev_h2osfc(c)
        clm_bgc_data%qflx_evap_soi_col(c)   = qflx_evap_soi(c)
        clm_bgc_data%qflx_sub_snow_col(c)   = qflx_sub_snow(c)
        clm_bgc_data%qflx_tran_veg_col(c)   = qflx_tran_veg(c)

        clm_bgc_data%htvp_col(c)            = htvp(c)
        clm_bgc_data%eflx_bot_col(c)        = eflx_bot(c)

!        do j = 1, nlevsoi
            clm_bgc_data%soilpsi_col(c,:)           = soilpsi(c,:)
            clm_bgc_data%rootfr_col(c,:)            = rootfr(c,:)

            clm_bgc_data%h2osoi_vol_col(c,:)        = h2osoi_vol(c,:)
            clm_bgc_data%h2osoi_liq_col(c,:)        = h2osoi_liq(c,:)
            clm_bgc_data%h2osoi_ice_col(c,:)        = h2osoi_ice(c,:)

            clm_bgc_data%t_soisno_col(c,:)          = t_soisno(c,:)

            clm_bgc_data%o2stress_unsat_col(c,:)        = o2stress_unsat(c,:)
            clm_bgc_data%o2stress_sat_col(c,:)          = o2stress_sat(c,:)
            clm_bgc_data%o2_decomp_depth_unsat_col(c,:) = o2_decomp_depth_unsat(c,:)
            clm_bgc_data%conc_o2_unsat_col(c,:)         = conc_o2_unsat(c,:)
            clm_bgc_data%o2_decomp_depth_sat_col(c,:)   = o2_decomp_depth_sat(c,:)
            clm_bgc_data%conc_o2_sat_col(c,:)           = conc_o2_sat(c,:)

!        end do
    end do

    ! CLM appears NO column-level ground-heat-flux variable, instead by 'patch'
    do fc = 1, num_soilc
        c = filter_soilc(fc)
        clm_bgc_data%eflx_soil_grnd_col(c)  = 0._r8
        clm_bgc_data%eflx_gnet_col(c)       = 0._r8
        do pftindex = 1, max_patch_per_col
            if (pftindex <= col_pp%npfts(c)) then
                p = col_pp%pfti(c) + pftindex - 1
                clm_bgc_data%eflx_soil_grnd_col(c)  = clm_bgc_data%eflx_soil_grnd_col(c) &
                                                    + eflx_soil_grnd_patch(p) * pft_pp%wtcol(p)           ! W/m2
                clm_bgc_data%eflx_gnet_col(c)       = clm_bgc_data%eflx_gnet_col(c) &
                                                    + eflx_gnet_patch(p) * pft_pp%wtcol(p)
            end if
       end do
    end do
!
!write(*,'(A30,12E14.6)')">>>DEBUG | soillsat=", soillsat_clmp_loc(1:10)
!write(*,'(A30,12E14.6)')">>>DEBUG | gsoilpsi[Pa]=", soilpsi_clmp_loc(1:10)
!write(*,'(A30,12E14.6)')">>>DEBUG | soilt[oC]=", soilt_clmp_loc(1:10)


   end associate
  end subroutine get_clm_soil_thermohydro
!!--------------------------------------------------------------------------------------

!!--------------------------------------------------------------------------------------
  subroutine get_clm_bgc_state(clm_bgc_data,                    &
                        bounds, num_soilc, filter_soilc,        &
                        carbonstate_vars, nitrogenstate_vars,   &
                        phosphorusstate_vars)

    !! get clm bgc state variables
    implicit none

    type(bounds_type)           , intent(in) :: bounds
    integer                     , intent(in) :: num_soilc         ! number of soil columns in filter
    integer                     , intent(in) :: filter_soilc(:)   ! filter for soil columns

    type(carbonstate_type)      , intent(in) :: carbonstate_vars
    type(nitrogenstate_type)    , intent(in) :: nitrogenstate_vars
    type(phosphorusstate_type)  , intent(in) :: phosphorusstate_vars
!    type(ch4_type)           , intent(in) :: ch4_vars

    type(clm_bgc_interface_data_type), intent(inout) :: clm_bgc_data

    character(len=256) :: subname = "get_clm_bgc_state"

    ! Local variables
    integer  :: fc, c, j, k
!    integer  :: gcount, cellcount
!    real(r8) :: wtgcell, realc_gcell, realn_gcell


    !------------------------------------------------------------------------------------------
    !
    associate ( &
       decomp_cpools_vr=> carbonstate_vars%decomp_cpools_vr_col     , &      ! (gC/m3) vertically-resolved decomposing (litter, cwd, soil) c pools
       decomp_npools_vr=> nitrogenstate_vars%decomp_npools_vr_col   , &      ! (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
       decomp_ppools_vr=> phosphorusstate_vars%decomp_ppools_vr_col , & ! [real(r8) (:,:,:) ! col (gP/m3) vertically-resolved decomposing (litter, cwd, soil) P pools

       smin_no3_vr     => nitrogenstate_vars%smin_no3_vr_col        , &      ! (gN/m3) vertically-resolved soil mineral NO3
       smin_nh4_vr     => nitrogenstate_vars%smin_nh4_vr_col        , &      ! (gN/m3) vertically-resolved soil mineral NH4
       smin_nh4sorb_vr => nitrogenstate_vars%smin_nh4sorb_vr_col    , &      ! (gN/m3) vertically-resolved soil mineral NH4 absorbed

       solutionp_vr    => phosphorusstate_vars%solutionp_vr_col     , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil solution P
       labilep_vr      => phosphorusstate_vars%labilep_vr_col       , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil labile mineral P
       secondp_vr      => phosphorusstate_vars%secondp_vr_col       , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil secondary mineralP
       sminp_vr        => phosphorusstate_vars%sminp_vr_col         , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil mineral P = solutionp + labilep + secondp
       occlp_vr        => phosphorusstate_vars%occlp_vr_col         , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil occluded mineral P
       primp_vr        => phosphorusstate_vars%primp_vr_col           & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil primary mineral P
    )
!
    do fc = 1, num_soilc
        c = filter_soilc(fc)
!        do j = 1, nlevdecomp
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
!        end do
    end do


  end associate
  end subroutine get_clm_bgc_state
!!--------------------------------------------------------------------------------------

!!--------------------------------------------------------------------------------------
  subroutine get_clm_bgc_flux(clm_bgc_data,                                 &
                        bounds, num_soilc, filter_soilc,                    &
                        cnstate_vars, carbonflux_vars, nitrogenflux_vars,   &
                        phosphorusflux_vars)

  !
  ! !DESCRIPTION:
  !! get clm bgc flux variables: external inputs to bgc state variables (pools)
  !
  ! !USES:

!    use clm_time_manager, only : get_step_size, get_nstep

!    use clm_varpar,       only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd


  ! !ARGUMENTS:
    implicit none

    type(bounds_type) , intent(in) :: bounds
    integer           , intent(in) :: num_soilc       ! number of soil columns in filter
    integer           , intent(in) :: filter_soilc(:) ! filter for soil columns

    type(cnstate_type)                  , intent(in)    :: cnstate_vars
    type(carbonflux_type)               , intent(in)    :: carbonflux_vars
    type(nitrogenflux_type)             , intent(in)    :: nitrogenflux_vars
    type(phosphorusflux_type)           , intent(in)    :: phosphorusflux_vars
    type(clm_bgc_interface_data_type)   , intent(inout) :: clm_bgc_data

    character(len=256) :: subname = "get_clm_bgc_flux"


 ! !LOCAL VARIABLES:
    integer  :: fc, c, g, j, k                         ! do loop indices
!    integer  :: gcount, cellcount
!    real(r8) :: wtgcell, realc_gcell, realn_gcell

    real(r8) :: dtime                               ! land model time step (sec)

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
      externalc_to_decomp_cpools_vr    => carbonflux_vars%externalc_to_decomp_cpools_col        , &
      externaln_to_decomp_npools_vr    => nitrogenflux_vars%externaln_to_decomp_npools_col      , &
      externalp_to_decomp_ppools_vr    => phosphorusflux_vars%externalp_to_decomp_ppools_col    , &
      ! inorg. nitrogen source
      ndep_to_sminn                    => nitrogenflux_vars%ndep_to_sminn_col                   , &
      nfix_to_sminn                    => nitrogenflux_vars%nfix_to_sminn_col                   , &
      fert_to_sminn                    => nitrogenflux_vars%fert_to_sminn_col                   , &
      soyfixn_to_sminn                 => nitrogenflux_vars%soyfixn_to_sminn_col                , &
      supplement_to_sminn_vr           => nitrogenflux_vars%supplement_to_sminn_vr_col          , &
      !
      nfixation_prof                   => cnstate_vars%nfixation_prof_col                       , &
      ndep_prof                        => cnstate_vars%ndep_prof_col                            , &

      no3_net_transport_vr             => nitrogenflux_vars%no3_net_transport_vr_col            , &
      col_plant_ndemand_vr             => nitrogenflux_vars%plant_ndemand_vr_col                , &

      plant_ndemand_col                => nitrogenflux_vars%plant_ndemand_col                   , &
      plant_pdemand_col                => phosphorusflux_vars%plant_pdemand_col                 , &

      col_plant_pdemand_vr             => phosphorusflux_vars%plant_pdemand_vr_col              , &
      pdep_to_sminp                    => phosphorusflux_vars%pdep_to_sminp_col                 , &
      ! assume pdep_prof = ndep_prof
      fert_p_to_sminp                  => phosphorusflux_vars%fert_p_to_sminp_col               , &
      supplement_to_sminp_vr           => phosphorusflux_vars%supplement_to_sminp_vr_col        , &
      sminp_net_transport_vr           => phosphorusflux_vars%sminp_net_transport_vr_col          &
    )

!    dtime = get_step_size()


!
    r_nh4_no3_dep(:)  = 1.0_r8      ! temporarily assuming half of N dep is in NH4 and another half in NO3
    r_nh4_no3_fert(:) = 1.0_r8      ! temporarily assiming half of N fertilization is in NH4 and another half in NO3
!
    do fc = 1,num_soilc
        c = filter_soilc(fc)
        clm_bgc_data%plant_ndemand_col(c)   = plant_ndemand_col(c)
        clm_bgc_data%plant_pdemand_col(c)   = plant_pdemand_col(c)

        fnh4_dep  = max(0._r8, min(1.0_r8, 1._r8/(r_nh4_no3_dep(c)+1._r8)))
        fnh4_fert = max(0._r8, min(1.0_r8, 1._r8/(r_nh4_no3_fert(c)+1._r8)))

!        do j = 1, nlevdecomp
            do k = 1, ndecomp_pools
                clm_bgc_data%externalc_to_decomp_cpools_col(c,:,k)  = externalc_to_decomp_cpools_vr(c,:,k)
                clm_bgc_data%externaln_to_decomp_npools_col(c,:,k)  = externaln_to_decomp_npools_vr(c,:,k)
                clm_bgc_data%externalp_to_decomp_ppools_col(c,:,k)  = externalp_to_decomp_ppools_vr(c,:,k)
            end do

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

            !! net flux to no3 = externaln_to_no3_col(c,j) - no3_net_transport_vr_col(c,j)
            clm_bgc_data%no3_net_transport_vr_col(c,:)      = no3_net_transport_vr(c,:)
            clm_bgc_data%sminp_net_transport_vr_col(c,:)    = sminp_net_transport_vr(c,:)  !!from solutionp

            clm_bgc_data%plant_ndemand_vr_col(c,:)          = col_plant_ndemand_vr(c,:)
            clm_bgc_data%plant_pdemand_vr_col(c,:)          = col_plant_pdemand_vr(c,:)

!        end do
    end do

    end associate
  end subroutine get_clm_bgc_flux
!!--------------------------------------------------------------------------------------

!!--------------------------------------------------------------------------------------
  subroutine update_soil_moisture(clm_bgc_data,     &
           bounds, num_soilc, filter_soilc,   &
           waterstate_vars)

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
!    type(soilstate_type) , intent(in)    :: soilstate_vars
    type(waterstate_type), intent(inout) :: waterstate_vars
    type(clm_bgc_interface_data_type), intent(in) :: clm_bgc_data

  ! !LOCAL VARIABLES:
    integer  :: fc, c, j, g, gcount      ! indices

  !EOP
  !-----------------------------------------------------------------------
    associate ( &
         h2osoi_liq_col =>  waterstate_vars%h2osoi_liq_col      , &
         h2osoi_ice_col =>  waterstate_vars%h2osoi_ice_col      , &
         h2osoi_vol_col =>  waterstate_vars%h2osoi_vol_col        &
    )

    do fc = 1,num_soilc
        c = filter_soilc(fc)
!        do j = 1, nlevsoi
            h2osoi_liq_col(c,:) =  clm_bgc_data%h2osoi_liq_col(c,:)
            h2osoi_ice_col(c,:) =  clm_bgc_data%h2osoi_ice_col(c,:)
            h2osoi_vol_col(c,:) =  clm_bgc_data%h2osoi_vol_col(c,:)
!        end do
    end do

    end associate
  end subroutine update_soil_moisture
!!--------------------------------------------------------------------------------------

!!--------------------------------------------------------------------------------------
  subroutine update_soil_temperature(clm_bgc_data,     &
           bounds, num_soilc, filter_soilc,   &
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
    type(clm_bgc_interface_data_type)   , intent(in)    :: clm_bgc_data

  ! !LOCAL VARIABLES:
    integer  :: fc, c, j, g, gcount      ! indices

  !EOP
  !-----------------------------------------------------------------------
    associate ( &
         t_soisno   => temperature_vars%t_soisno_col   & ! snow-soil temperature (Kelvin)
    )
    do fc = 1,num_soilc
        c = filter_soilc(fc)
!        do j = 1, nlevsoi
            t_soisno(c,:)   = clm_bgc_data%t_soisno_col(c,:)
!        end do
    end do

    end associate
  end subroutine update_soil_temperature
!!--------------------------------------------------------------------------------------

!!--------------------------------------------------------------------------------------
  subroutine update_bgc_state_decomp(clm_bgc_data,    &
           bounds, num_soilc, filter_soilc,         &
           carbonstate_vars, nitrogenstate_vars,    &
           phosphorusstate_vars                     &
           )

!    use CNDecompCascadeConType, only : decomp_cascade_con
!    use clm_time_manager, only : get_step_size
!#ifndef FLEXIBLE_POOLS
!    use clm_varpar, only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
!#endif

    implicit none

    type(bounds_type)                   , intent(in)    :: bounds
    integer                             , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                             , intent(in)    :: filter_soilc(:)   ! filter for soil columns

    type(carbonstate_type)              , intent(inout) :: carbonstate_vars
    type(nitrogenstate_type)            , intent(inout) :: nitrogenstate_vars
    type(phosphorusstate_type)          , intent(inout) :: phosphorusstate_vars

    type(clm_bgc_interface_data_type)   , intent(in)    :: clm_bgc_data

    character(len=256) :: subname = "update_soil_bgc_state"

    integer  :: fc,c,j,k
!    integer  :: gcount, cellcount
!    real(r8) :: wtgcell

!    real(r8) :: dtime            ! land model time step (sec)

!------------------------------------------------------------------------------------
     !
     associate ( &
     decomp_cpools_vr             => carbonstate_vars%decomp_cpools_vr_col           , &
     decomp_npools_vr             => nitrogenstate_vars%decomp_npools_vr_col         , &
     decomp_ppools_vr             => phosphorusstate_vars%decomp_ppools_vr_col         &
     )
! ------------------------------------------------------------------------
!     dtime = get_step_size()
!
    do fc = 1, num_soilc
        c = filter_soilc(fc)
!        do j = 1, nlevdecomp
            do k = 1, ndecomp_pools
                decomp_cpools_vr(c,:,k) = clm_bgc_data%decomp_cpools_vr_col(c,:,k)
                decomp_npools_vr(c,:,k) = clm_bgc_data%decomp_npools_vr_col(c,:,k)
                decomp_ppools_vr(c,:,k) = clm_bgc_data%decomp_ppools_vr_col(c,:,k)
            end do
!        end do
    end do

    end associate
  end subroutine update_bgc_state_decomp
!!--------------------------------------------------------------------------------------

!!--------------------------------------------------------------------------------------
    subroutine update_bgc_state_smin(clm_bgc_data,      &
           bounds, num_soilc, filter_soilc,             &
           nitrogenstate_vars, phosphorusstate_vars)

    use CNDecompCascadeConType, only : decomp_cascade_con
    use clm_time_manager, only : get_step_size
!#ifndef FLEXIBLE_POOLS
!    use clm_varpar, only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
!#endif

    implicit none

    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                     , intent(in)    :: filter_soilc(:)   ! filter for soil columns

    type(nitrogenstate_type)    , intent(inout) :: nitrogenstate_vars
    type(phosphorusstate_type)  , intent(inout) :: phosphorusstate_vars

    type(clm_bgc_interface_data_type), intent(in) :: clm_bgc_data

    character(len=256) :: subname = "update_bgc_state_smin"

    integer  :: fc,c,j
!    integer  :: gcount, cellcount
!    real(r8) :: wtgcell
!
!    real(r8) :: dtime            ! land model time step (sec)

!------------------------------------------------------------------------------------
     !
     associate ( &
     sminn_vr           => nitrogenstate_vars%sminn_vr_col           , &
     smin_no3_vr        => nitrogenstate_vars%smin_no3_vr_col        , &
     smin_nh4_vr        => nitrogenstate_vars%smin_nh4_vr_col        , &
     smin_nh4sorb_vr    => nitrogenstate_vars%smin_nh4sorb_vr_col    , &

     solutionp_vr       => phosphorusstate_vars%solutionp_vr_col     , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil solution P
     labilep_vr         => phosphorusstate_vars%labilep_vr_col       , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil labile mineral P
     secondp_vr         => phosphorusstate_vars%secondp_vr_col       , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil secondary mineralP
     sminp_vr           => phosphorusstate_vars%sminp_vr_col         , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil mineral P = solutionp + labilep + second
     occlp_vr           => phosphorusstate_vars%occlp_vr_col         , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil occluded mineral P
     primp_vr           => phosphorusstate_vars%primp_vr_col           & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil primary mineral P

     )
! ------------------------------------------------------------------------
!     dtime = get_step_size()
!
    do fc = 1, num_soilc
        c = filter_soilc(fc)
!        do j = 1, nlevdecomp
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
!        end do
    end do
!write(*,'(A30,12E14.6)')"DEBUG | clm UPDATE no3=",smin_no3_vr(1,1:nlevdecomp)

    end associate
  end subroutine update_bgc_state_smin
!!--------------------------------------------------------------------------------------

!!--------------------------------------------------------------------------------------
    subroutine update_bgc_flux_decomp_sourcesink(clm_bgc_data,       &
           bounds, num_soilc, filter_soilc,              &
           carbonflux_vars, nitrogenflux_vars, &
           phosphorusflux_vars)

    use CNDecompCascadeConType, only : decomp_cascade_con
    use clm_time_manager, only : get_step_size
!#ifndef FLEXIBLE_POOLS
!    use clm_varpar, only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
!#endif

    implicit none

    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                     , intent(in)    :: filter_soilc(:)   ! filter for soil columns

    type(carbonflux_type)       , intent(inout) :: carbonflux_vars
    type(nitrogenflux_type)     , intent(inout) :: nitrogenflux_vars
    type(phosphorusflux_type)   , intent(inout) :: phosphorusflux_vars

    type(clm_bgc_interface_data_type), intent(in) :: clm_bgc_data

    integer :: fc, c, j, k
    character(len=256) :: subname = "update_soil_bgc_pf2clm"

    associate ( &
     decomp_cpools_sourcesink_vr  => carbonflux_vars%decomp_cpools_sourcesink_col    , &
     decomp_npools_sourcesink_vr  => nitrogenflux_vars%decomp_npools_sourcesink_col  , &
     decomp_ppools_sourcesink_vr  => phosphorusflux_vars%decomp_ppools_sourcesink_col  &
     )

    do fc = 1, num_soilc
        c = filter_soilc(fc)
!        do j = 1, nlevdecomp
            do k = 1, ndecomp_pools
                decomp_cpools_sourcesink_vr(c,:,k) = clm_bgc_data%decomp_cpools_sourcesink_col(c,:,k)
                decomp_npools_sourcesink_vr(c,:,k) = clm_bgc_data%decomp_npools_sourcesink_col(c,:,k)
                decomp_ppools_sourcesink_vr(c,:,k) = clm_bgc_data%decomp_ppools_sourcesink_col(c,:,k)
            end do
!        end do
    end do
    end associate
    end subroutine update_bgc_flux_decomp_sourcesink
!!--------------------------------------------------------------------------------------

!!--------------------------------------------------------------------------------------
    subroutine update_bgc_flux_decomp_cascade(clm_bgc_data, &
           bounds, num_soilc, filter_soilc,                 &
           carbonflux_vars, nitrogenflux_vars,              &
           phosphorusflux_vars)

    use CNDecompCascadeConType, only : decomp_cascade_con
    use clm_time_manager, only : get_step_size
!#ifndef FLEXIBLE_POOLS
!    use clm_varpar, only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
!#endif

    implicit none

    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                     , intent(in)    :: filter_soilc(:)   ! filter for soil columns

    type(carbonflux_type)       , intent(inout) :: carbonflux_vars
    type(nitrogenflux_type)     , intent(inout) :: nitrogenflux_vars
    type(phosphorusflux_type)   , intent(inout) :: phosphorusflux_vars

    type(clm_bgc_interface_data_type), intent(in) :: clm_bgc_data

    integer :: fc, c, j, k
    character(len=256) :: subname = "update_soil_bgc_pf2clm"

    associate ( &
     phr_vr                           => carbonflux_vars%phr_vr_col                                 , & ! Output: [real(r8) (:,:)   ]  potential HR (gC/m3/s)
     fphr                             => carbonflux_vars%fphr_col                                   , & ! Output: [real(r8) (:,:)   ]  fraction of potential SOM + LITTER heterotrophic

     decomp_cascade_hr_vr_col         => carbonflux_vars%decomp_cascade_hr_vr_col                   , &

     decomp_cascade_ctransfer_vr_col  => carbonflux_vars%decomp_cascade_ctransfer_vr_col            , &
     decomp_cascade_ntransfer_vr_col  => nitrogenflux_vars%decomp_cascade_ntransfer_vr_col          , &
     decomp_cascade_ptransfer_vr_col  => phosphorusflux_vars%decomp_cascade_ptransfer_vr_col        , &

     decomp_cascade_sminn_flux_vr_col => nitrogenflux_vars%decomp_cascade_sminn_flux_vr_col         , & ! Output: [real(r8) (:,:,:) ]  vert-res mineral N flux for transition along decomposition cascade (gN/m3/s)
     decomp_cascade_sminp_flux_vr_col => phosphorusflux_vars%decomp_cascade_sminp_flux_vr_col       , & ! Output: [real(r8) (:,:,:) ]  vert-res mineral P flux for transition along decomposition cascade (gP/m3/s)

     sminn_to_denit_decomp_cascade_vr_col => nitrogenflux_vars%sminn_to_denit_decomp_cascade_vr_col   & ! Output: [real(r8) (:,:,:) ]
     )

    do fc = 1, num_soilc
        c = filter_soilc(fc)
!        do j = 1, nlevdecomp

            phr_vr(c,:) = clm_bgc_data%phr_vr_col(c,:)
            fphr(c,:)   = clm_bgc_data%fphr_col(c,:)

            do k = 1, ndecomp_pools
                decomp_cascade_hr_vr_col(c,:,k)         = clm_bgc_data%decomp_cascade_hr_vr_col(c,:,k)
                decomp_cascade_ctransfer_vr_col(c,:,k)  = clm_bgc_data%decomp_cascade_ctransfer_vr_col(c,:,k)
                decomp_cascade_ntransfer_vr_col(c,:,k)  = clm_bgc_data%decomp_cascade_ntransfer_vr_col(c,:,k)
                decomp_cascade_ptransfer_vr_col(c,:,k)  = clm_bgc_data%decomp_cascade_ptransfer_vr_col(c,:,k)

                decomp_cascade_sminn_flux_vr_col(c,:,k)     = clm_bgc_data%decomp_cascade_sminn_flux_vr_col(c,:,k)
                decomp_cascade_sminp_flux_vr_col(c,:,k)     = clm_bgc_data%decomp_cascade_sminp_flux_vr_col(c,:,k)
                sminn_to_denit_decomp_cascade_vr_col(c,:,k) = clm_bgc_data%sminn_to_denit_decomp_cascade_vr_col(c,:,k)
            end do
!        end do
    end do
    end associate
    end subroutine update_bgc_flux_decomp_cascade
!!--------------------------------------------------------------------------------------

!!--------------------------------------------------------------------------------------
    subroutine update_bgc_flux_smin(clm_bgc_data,   &
           bounds, num_soilc, filter_soilc,         &
           cnstate_vars,                            &
           nitrogenflux_vars, phosphorusflux_vars)

!    use CNDecompCascadeConType, only : decomp_cascade_con
!    use clm_time_manager, only : get_step_size
!#ifndef FLEXIBLE_POOLS
!    use clm_varpar, only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
!#endif

    implicit none

    type(bounds_type)                   , intent(in)    :: bounds
    integer                             , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                             , intent(in)    :: filter_soilc(:)   ! filter for soil columns

    type(cnstate_type)                  , intent(inout) :: cnstate_vars

    type(nitrogenflux_type)             , intent(inout) :: nitrogenflux_vars
    type(phosphorusflux_type)           , intent(inout) :: phosphorusflux_vars
    type(clm_bgc_interface_data_type)   , intent(in)    :: clm_bgc_data

    integer :: fc, c, j
    character(len=256) :: subname = "update_bgc_flux_smin"

    associate ( &
     fpg                          => cnstate_vars%fpg_col                            , & ! Output: [real(r8) (:)   ]  fraction of potential gpp (no units)
     fpi                          => cnstate_vars%fpi_col                            , & ! Output: [real(r8) (:)   ]  fraction of potential immobilization (no units)
     fpi_vr                       => cnstate_vars%fpi_vr_col                         , & ! Output: [real(r8) (:,:) ]  fraction of potential immobilization (no units)
     fpg_p                        => cnstate_vars%fpg_p_col                          , & ! Output: [real(r8) (:)   ]  fraction of potential gpp (no units)
     fpi_p                        => cnstate_vars%fpi_p_col                          , & ! Output: [real(r8) (:)   ]  fraction of potential immobilization (no units)
     fpi_p_vr                     => cnstate_vars%fpi_p_vr_col                       , & ! Output: [real(r8) (:,:) ]  fraction of potential immobilization (no units)

     potential_immob              => nitrogenflux_vars%potential_immob_col           , & ! Output: [real(r8) (:)   ]
     actual_immob                 => nitrogenflux_vars%actual_immob_col              , & ! Output: [real(r8) (:)   ]
     sminn_to_plant               => nitrogenflux_vars%sminn_to_plant_col            , & ! Output: [real(r8) (:)   ]

     sminn_to_plant_vr            => nitrogenflux_vars%sminn_to_plant_vr_col         , &
     smin_no3_to_plant_vr         => nitrogenflux_vars%smin_no3_to_plant_vr_col      , &
     smin_nh4_to_plant_vr         => nitrogenflux_vars%smin_nh4_to_plant_vr_col      , &
     potential_immob_vr           => nitrogenflux_vars%potential_immob_vr_col        , &
     actual_immob_vr              => nitrogenflux_vars%actual_immob_vr_col           , &
     actual_immob_no3_vr          => nitrogenflux_vars%actual_immob_no3_vr_col       , & ! Output: [real(r8) (:,:) ]
     actual_immob_nh4_vr          => nitrogenflux_vars%actual_immob_nh4_vr_col       , & ! Output: [real(r8) (:,:) ]
     gross_nmin_vr                => nitrogenflux_vars%gross_nmin_vr_col             , &
     net_nmin_vr                  => nitrogenflux_vars%net_nmin_vr_col               , & ! Output: [real(r8) (:,:)   ]

     sminn_to_denit_excess_vr     => nitrogenflux_vars%sminn_to_denit_excess_vr_col  , & ! Output: [real(r8) (:,:) ]
     supplement_to_sminn_vr       => nitrogenflux_vars%supplement_to_sminn_vr_col    , & ! Output: [real(r8) (:,:) ]

     potential_immob_p            => phosphorusflux_vars%potential_immob_p_col       , & ! Output: [real(r8) (:)   ]
     actual_immob_p               => phosphorusflux_vars%actual_immob_p_col          , & ! Output: [real(r8) (:)   ]
     sminp_to_plant               => phosphorusflux_vars%sminp_to_plant_col          , & ! Output: [real(r8) (:)   ]

     supplement_to_sminp_vr       => phosphorusflux_vars%supplement_to_sminp_vr_col  , & ! Output: [real(r8) (:,:) ]

     sminp_to_plant_vr            => phosphorusflux_vars%sminp_to_plant_vr_col       , &
     potential_immob_p_vr         => phosphorusflux_vars%potential_immob_p_vr_col    , & ! Output: [real(r8) (:,:)   ]
     actual_immob_p_vr            => phosphorusflux_vars%actual_immob_p_vr_col       , &
     gross_pmin_vr                => phosphorusflux_vars%gross_pmin_vr_col           , & ! Output: [real(r8) (:,:)   ]
     net_pmin_vr                  => phosphorusflux_vars%net_pmin_vr_col               & ! Output: [real(r8) (:,:)   ]
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

!        do j = 1, nlevdecomp
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
            net_nmin_vr(c,:)            = clm_bgc_data%net_nmin_vr_col(c,:)     !!NOT available in PF

            sminn_to_denit_excess_vr(c,:)=clm_bgc_data%sminn_to_denit_excess_vr_col(c,:)
            supplement_to_sminn_vr(c,:) = clm_bgc_data%supplement_to_sminn_vr_col(c,:)

            supplement_to_sminp_vr(c,:) = clm_bgc_data%supplement_to_sminp_vr_col(c,:)

            sminp_to_plant_vr(c,:)      = clm_bgc_data%sminp_to_plant_vr_col(c,:)
            potential_immob_p_vr(c,:)   = clm_bgc_data%potential_immob_p_vr_col(c,:)
            actual_immob_p_vr(c,:)      = clm_bgc_data%actual_immob_p_vr_col(c,:)
            gross_pmin_vr(c,:)          = clm_bgc_data%gross_pmin_vr_col(c,:)
            net_pmin_vr(c,:)            = clm_bgc_data%net_pmin_vr_col(c,:)     !!NOT available in PF

!        end do
    end do
    end associate
    end subroutine update_bgc_flux_smin
!!--------------------------------------------------------------------------------------

!!--------------------------------------------------------------------------------------
    subroutine update_bgc_flux_nitdenit(clm_bgc_data,   &
           bounds, num_soilc, filter_soilc,         &
           nitrogenflux_vars, phosphorusflux_vars)

!    use CNDecompCascadeConType, only : decomp_cascade_con
!    use clm_time_manager, only : get_step_size
!#ifndef FLEXIBLE_POOLS
!    use clm_varpar, only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
!#endif

    implicit none

    type(bounds_type)                   , intent(in)    :: bounds
    integer                             , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                             , intent(in)    :: filter_soilc(:)   ! filter for soil columns

    type(nitrogenflux_type)             , intent(inout) :: nitrogenflux_vars
    type(phosphorusflux_type)           , intent(inout) :: phosphorusflux_vars
    type(clm_bgc_interface_data_type)   , intent(in)    :: clm_bgc_data

    integer :: fc, c, j
    character(len=256) :: subname = "update_bgc_flux_nitdenit"

    associate ( &
         pot_f_nit_vr                 => nitrogenflux_vars%pot_f_nit_vr_col                  , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) potential soil nitrification flux
         pot_f_denit_vr               => nitrogenflux_vars%pot_f_denit_vr_col                , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) potential soil denitrification flux
         f_nit_vr                     => nitrogenflux_vars%f_nit_vr_col                      , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) soil nitrification flux
         f_denit_vr                   => nitrogenflux_vars%f_denit_vr_col                    , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) soil denitrification flux
         n2_n2o_ratio_denit_vr        => nitrogenflux_vars%n2_n2o_ratio_denit_vr_col         , & ! Output: [real(r8) (:,:) ]  ratio of N2 to N2O production by denitrification [gN/gN]
         f_n2o_denit_vr               => nitrogenflux_vars%f_n2o_denit_vr_col                , & ! Output: [real(r8) (:,:) ]  flux of N2O from denitrification [gN/m3/s]
         f_n2o_nit_vr                 => nitrogenflux_vars%f_n2o_nit_vr_col                    & ! Output: [real(r8) (:,:) ]  flux of N2O from nitrification [gN/m3/s]
    )

    do fc = 1, num_soilc
        c = filter_soilc(fc)
!        do j = 1, nlevdecomp

            pot_f_nit_vr(c,:)           = clm_bgc_data%pot_f_nit_vr_col(c,:)
            pot_f_denit_vr(c,:)         = clm_bgc_data%pot_f_denit_vr_col(c,:)
            f_nit_vr(c,:)               = clm_bgc_data%f_nit_vr_col(c,:)
            f_denit_vr(c,:)             = clm_bgc_data%f_denit_vr_col(c,:)
            n2_n2o_ratio_denit_vr(c,:)  = clm_bgc_data%n2_n2o_ratio_denit_vr_col(c,:)
            f_n2o_denit_vr(c,:)         = clm_bgc_data%f_n2o_denit_vr_col(c,:)
            f_n2o_nit_vr(c,:)           = clm_bgc_data%f_n2o_nit_vr_col(c,:)
!        end do
    end do
    end associate
    end subroutine update_bgc_flux_nitdenit
!!--------------------------------------------------------------------------------------

!!--------------------------------------------------------------------------------------
  subroutine update_bgc_flux_gas_pf(clm_bgc_data,  &
     bounds, num_soilc, filter_soilc,           &
     carbonflux_vars, nitrogenflux_vars)

!     use clm_time_manager, only : get_step_size, get_nstep

     ! PFLOTRAN gas fluxes
     implicit none

     type(bounds_type)                  , intent(in)    :: bounds
     integer                            , intent(in)    :: num_soilc       ! number of soil columns in filter
     integer                            , intent(in)    :: filter_soilc(:) ! filter for soil columns

     type(carbonflux_type)              , intent(inout) :: carbonflux_vars
     type(nitrogenflux_type)            , intent(inout) :: nitrogenflux_vars
     type(clm_bgc_interface_data_type)  , intent(in)    :: clm_bgc_data

     !character(len=256) :: subname = "get_pf_bgc_gaslosses"

     integer  :: fc, c, g, j
!     integer  :: gcount, cellcount
!     real(r8) :: dtime            ! land model time step (sec)
!     integer  :: nstep


!------------------------------------------------------------------------------------
    associate ( &
     hr_vr                        => carbonflux_vars%hr_vr_col              , &
     f_co2_soil_vr                => carbonflux_vars%f_co2_soil_vr_col      , &
     f_n2o_soil_vr                => nitrogenflux_vars%f_n2o_soil_vr_col    , &
     f_n2_soil_vr                 => nitrogenflux_vars%f_n2_soil_vr_col     , &
     f_ngas_decomp_vr             => nitrogenflux_vars%f_ngas_decomp_vr_col , &
     f_ngas_nitri_vr              => nitrogenflux_vars%f_ngas_nitri_vr_col  , &
     f_ngas_denit_vr              => nitrogenflux_vars%f_ngas_denit_vr_col    &
     )
! ------------------------------------------------------------------------
    do fc = 1,num_soilc
        c = filter_soilc(fc)
!        do j = 1, nlevdecomp
!
              f_co2_soil_vr(c,:)         = clm_bgc_data%f_co2_soil_vr_col(c,:)
              f_n2_soil_vr(c,:)          = clm_bgc_data%f_n2_soil_vr_col(c,:)
              f_n2o_soil_vr(c,:)         = clm_bgc_data%f_n2o_soil_vr_col(c,:)

              hr_vr(c,:)                 = clm_bgc_data%hr_vr_col(c,:)
              f_ngas_decomp_vr(c,:)      = clm_bgc_data%f_ngas_decomp_vr_col(c,:)
              f_ngas_nitri_vr(c,:)       = clm_bgc_data%f_ngas_nitri_vr_col(c,:)
              f_ngas_denit_vr(c,:)       = clm_bgc_data%f_ngas_denit_vr_col(c,:)
!       enddo
     enddo ! do c = begc, endc
!
    end associate
  end subroutine update_bgc_flux_gas_pf
!!--------------------------------------------------------------------------------------

!!--------------------------------------------------------------------------------------
  subroutine update_bgc_data_pf2clm(clm_bgc_data, bounds,         &
           num_soilc, filter_soilc,                               &
           num_soilp, filter_soilp,                               &
           atm2lnd_vars,                                          &
           waterstate_vars, waterflux_vars,                       &
           soilstate_vars,  temperature_vars, energyflux_vars,    &
           soilhydrology_vars, soil_water_retention_curve,        &
           cnstate_vars, carbonflux_vars, carbonstate_vars,       &
           nitrogenflux_vars, nitrogenstate_vars,                 &
           phosphorusflux_vars, phosphorusstate_vars,             &
           ch4_vars)
    !! USES
    use clm_varctl          , only : use_pflotran, pf_tmode, pf_hmode, pf_cmode

    implicit none

    ! !ARGUMENTS:
    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                     , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                     , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                     , intent(in)    :: filter_soilp(:)   ! filter for soil patches
    type(atm2lnd_type)          , intent(in)    :: atm2lnd_vars
    type(waterstate_type)       , intent(inout) :: waterstate_vars
    type(waterflux_type)        , intent(inout) :: waterflux_vars
    type(soilstate_type)        , intent(inout) :: soilstate_vars
    type(temperature_type)      , intent(inout) :: temperature_vars
    type(soilhydrology_type)    , intent(inout) :: soilhydrology_vars
    type(energyflux_type)       , intent(inout) :: energyflux_vars

    type(cnstate_type)          , intent(inout) :: cnstate_vars
    type(carbonflux_type)       , intent(inout) :: carbonflux_vars
    type(carbonstate_type)      , intent(inout) :: carbonstate_vars
    type(nitrogenflux_type)     , intent(inout) :: nitrogenflux_vars
    type(nitrogenstate_type)    , intent(inout) :: nitrogenstate_vars
    type(phosphorusflux_type)   , intent(inout) :: phosphorusflux_vars
    type(phosphorusstate_type)  , intent(inout) :: phosphorusstate_vars
    type(ch4_type)              , intent(inout) :: ch4_vars

    class(soil_water_retention_curve_type)  , intent(in) :: soil_water_retention_curve
    type(clm_bgc_interface_data_type)       , intent(in) :: clm_bgc_data

    !-----------------------------------------------------------------------

    character(len=256) :: subname = "get_clm_bgc_data"


    if (pf_cmode) then
        !! bgc_state_decomp is updated in CLM
        !! by passing bgc_flux_decomp_sourcesink into CNSoilLittVertTransp
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

    if (pf_tmode) then
        call update_soil_temperature(clm_bgc_data,      &
                   bounds, num_soilc, filter_soilc,     &
                   temperature_vars)
    end if

    if (pf_hmode) then
        call update_soil_moisture(clm_bgc_data,         &
                   bounds, num_soilc, filter_soilc,     &
                   waterstate_vars)
    end if

  end subroutine update_bgc_data_pf2clm
!!--------------------------------------------------------------------------------------

!!--------------------------------------------------------------------------------------


!!--------------------------------------------------------------------------------------
! BEG of CLM-bgc through interface
!!--------------------------------------------------------------------------------------
  ! !INTERFACE:
  subroutine clm_bgc_run(clm_bgc_data, bounds,              &
                num_soilc, filter_soilc,                    &
                num_soilp, filter_soilp,                    &
                canopystate_vars, soilstate_vars,           &
                temperature_vars, waterstate_vars,          &
                cnstate_vars, ch4_vars,                     &
                carbonstate_vars, carbonflux_vars,          &
                nitrogenstate_vars, nitrogenflux_vars,      &
                phosphorusstate_vars,phosphorusflux_vars)

    !! USES:
    use CNDecompMod          , only: CNDecompAlloc

    !! ARGUMENTS:
    type(bounds_type)                   , intent(in)    :: bounds
    integer                             , intent(in)    :: num_soilc          ! number of soil columns in filter
    integer                             , intent(in)    :: filter_soilc(:)    ! filter for soil columns
    integer                             , intent(in)    :: num_soilp          ! number of soil patches in filter
    integer                             , intent(in)    :: filter_soilp(:)    ! filter for soil patches
!    type(photosyns_type)                , intent(in)    :: photosyns_vars
    type(canopystate_type)              , intent(inout) :: canopystate_vars
    type(soilstate_type)                , intent(inout) :: soilstate_vars
    type(temperature_type)              , intent(inout) :: temperature_vars
    type(waterstate_type)               , intent(inout) :: waterstate_vars
    type(cnstate_type)                  , intent(inout) :: cnstate_vars
    type(ch4_type)                      , intent(inout) :: ch4_vars
    type(carbonstate_type)              , intent(inout) :: carbonstate_vars
    type(carbonflux_type)               , intent(inout) :: carbonflux_vars
!    type(carbonflux_type)               , intent(inout) :: c13_carbonflux_vars
!    type(carbonflux_type)               , intent(inout) :: c14_carbonflux_vars
    type(nitrogenstate_type)            , intent(inout) :: nitrogenstate_vars
    type(nitrogenflux_type)             , intent(inout) :: nitrogenflux_vars
!    type(crop_type)                     , intent(inout) :: crop_vars
    type(phosphorusstate_type)          , intent(inout) :: phosphorusstate_vars
    type(phosphorusflux_type)           , intent(inout) :: phosphorusflux_vars

    type(clm_bgc_interface_data_type)   , intent(inout) :: clm_bgc_data

    !!-------------------------------------------------------------
    !! STEP-2: (i) pass data from clm_bgc_data to CNDecompAlloc
    call clm_bgc_get_data(clm_bgc_data, bounds,             &
                num_soilc, filter_soilc,                    &
                canopystate_vars, soilstate_vars,           &
                temperature_vars, waterstate_vars,          &
                cnstate_vars, ch4_vars,                     &
                carbonstate_vars, carbonflux_vars,          &
                nitrogenstate_vars, nitrogenflux_vars,      &
                phosphorusstate_vars,phosphorusflux_vars)

    !! STEP-2: (ii) run CNDecompAlloc
    call CNDecompAlloc (bounds, num_soilc, filter_soilc,    &
               num_soilp, filter_soilp,                     &
               canopystate_vars, soilstate_vars,            &
               temperature_vars, waterstate_vars,           &
               cnstate_vars, ch4_vars,                      &
               carbonstate_vars, carbonflux_vars,           &
               nitrogenstate_vars, nitrogenflux_vars,       &
               phosphorusstate_vars,phosphorusflux_vars)

    !! STEP-2: (iii) update clm_bgc_data from CNDecompAlloc
    call clm_bgc_update_data(clm_bgc_data, bounds,          &
                num_soilc, filter_soilc,                    &
                cnstate_vars, carbonflux_vars,              &
                nitrogenflux_vars, phosphorusflux_vars)

  end subroutine clm_bgc_run
!!--------------------------------------------------------------------------------------

!!--------------------------------------------------------------------------------------
  ! !INTERFACE:
  !! pass data from clm_bgc_data to clm original data-types that used by CNDecompAlloc
  subroutine clm_bgc_get_data(clm_bgc_data, bounds,     &
            num_soilc, filter_soilc,                    &
            canopystate_vars, soilstate_vars,           &
            temperature_vars, waterstate_vars,          &
            cnstate_vars, ch4_vars,                     &
            carbonstate_vars, carbonflux_vars,          &
            nitrogenstate_vars, nitrogenflux_vars,      &
            phosphorusstate_vars,phosphorusflux_vars)

    !! USES:


    !! ARGUMENTS:
    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: num_soilc          ! number of soil columns in filter
    integer                     , intent(in)    :: filter_soilc(:)    ! filter for soil columns
!    integer                     , intent(in)    :: num_soilp          ! number of soil patches in filter
!    integer                     , intent(in)    :: filter_soilp(:)    ! filter for soil patches
!    type(photosyns_type)        , intent(in)    :: photosyns_vars
    type(canopystate_type)      , intent(inout) :: canopystate_vars
    type(soilstate_type)        , intent(inout) :: soilstate_vars
    type(temperature_type)      , intent(inout) :: temperature_vars
    type(waterstate_type)       , intent(inout) :: waterstate_vars
    type(cnstate_type)          , intent(inout) :: cnstate_vars
    type(ch4_type)              , intent(inout) :: ch4_vars
    type(carbonstate_type)      , intent(inout) :: carbonstate_vars
    type(carbonflux_type)       , intent(inout) :: carbonflux_vars
!    type(carbonflux_type)       , intent(inout) :: c13_carbonflux_vars
!    type(carbonflux_type)       , intent(inout) :: c14_carbonflux_vars
    type(nitrogenstate_type)    , intent(inout) :: nitrogenstate_vars
    type(nitrogenflux_type)     , intent(inout) :: nitrogenflux_vars
!    type(crop_type)             , intent(in)    :: crop_vars
    type(phosphorusstate_type)  , intent(inout) :: phosphorusstate_vars
    type(phosphorusflux_type)   , intent(inout) :: phosphorusflux_vars

    type(clm_bgc_interface_data_type), intent(in) :: clm_bgc_data

    !! LOCAL VARIABLES:
    integer :: fc, c, j, k
    !-----------------------------------------------------------------------

    associate(&
!         initial_cn_ratio                 =>    decomp_cascade_con%initial_cn_ratio                    , & ! Input:  [real(r8) (:)     ]  c:n ratio for initialization of pools
!         initial_cp_ratio                 =>    decomp_cascade_con%initial_cp_ratio                    , & ! Input:  [real(r8) (:)     ]  c:p ratio for initialization of pools
        decomp_cpools_vr        => carbonstate_vars%decomp_cpools_vr_col        , & ! Input:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
        decomp_npools_vr        => nitrogenstate_vars%decomp_npools_vr_col      , & ! Input:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
        decomp_ppools_vr        => phosphorusstate_vars%decomp_ppools_vr_col    , & ! Input:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) P pools

        smin_no3_vr             => nitrogenstate_vars%smin_no3_vr_col           , &      ! (gN/m3) vertically-resolved soil mineral NO3
        smin_nh4_vr             => nitrogenstate_vars%smin_nh4_vr_col           , &      ! (gN/m3) vertically-resolved soil mineral NH4
        smin_nh4sorb_vr         => nitrogenstate_vars%smin_nh4sorb_vr_col       , &      ! (gN/m3) vertically-resolved soil mineral NH4 absorbed

        solutionp_vr            => phosphorusstate_vars%solutionp_vr_col        , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil solution P
        labilep_vr              => phosphorusstate_vars%labilep_vr_col          , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil labile mineral P
        secondp_vr              => phosphorusstate_vars%secondp_vr_col          , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil secondary mineralP
        sminp_vr                => phosphorusstate_vars%sminp_vr_col            , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil mineral P = solutionp + labilep + secondp
        occlp_vr                => phosphorusstate_vars%occlp_vr_col            , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil occluded mineral P
        primp_vr                => phosphorusstate_vars%primp_vr_col            , & ! [real(r8) (:,:)   ! col (gP/m3) vertically-resolved soil primary mineral P

        plant_ndemand_col       => nitrogenflux_vars%plant_ndemand_col          , &
        plant_pdemand_col       => phosphorusflux_vars%plant_pdemand_col        , &

        alt_indx                => canopystate_vars%alt_indx_col                , & ! Input:  [integer  (:)     ]  current depth of thaw

        watsat                  => soilstate_vars%watsat_col                    , & ! Input:  [real(r8) (:,:)  ]  volumetric soil water at saturation (porosity) (nlevgrnd)
        bd                      => soilstate_vars%bd_col                        , & ! Input:  [real(r8) (:,:)  ]  bulk density of dry soil material [kg/m3]
        watfc                   => soilstate_vars%watfc_col                     , & ! Input:  [real(r8) (:,:)  ]  volumetric soil water at field capacity (nlevsoi)
        bsw                     => soilstate_vars%bsw_col                       , & ! Input:  [real(r8) (:,:)  ]  Clapp and Hornberger "b" (nlevgrnd)
        cellorg                 => soilstate_vars%cellorg_col                   , & ! Input:  [real(r8) (:,:)  ]  column 3D org (kg/m3 organic matter) (nlevgrnd)
        sucsat                  => soilstate_vars%sucsat_col                    , & ! Input:  [real(r8) (:,:)  ]  minimum soil suction (mm)
        soilpsi                 => soilstate_vars%soilpsi_col                   , & ! Input:  [real(r8) (:,:)  ]  soil water potential in each soil layer (MPa)

        h2osoi_vol              => waterstate_vars%h2osoi_vol_col               , & ! Input:  [real(r8) (:,:)  ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevgrnd)
        h2osoi_liq              => waterstate_vars%h2osoi_liq_col               , & ! Input:  [real(r8) (:,:)  ]  liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)

        t_soisno                => temperature_vars%t_soisno_col                , & ! Input:  [real(r8) (:,:)  ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)

        o2_decomp_depth_unsat   => ch4_vars%o2_decomp_depth_unsat_col           , & ! Input:  [real(r8) (:,:)  ]  O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
        o2_decomp_depth_sat     => ch4_vars%o2_decomp_depth_sat_col             , & ! Input:  [real(r8) (:,:)  ]  O2 consumption during decomposition in each soil layer (nlevsoi) (mol/m3/s)
        conc_o2_unsat           => ch4_vars%conc_o2_unsat_col                   , & ! Input:  [real(r8) (:,:)  ]  O2 conc in each soil layer (mol/m3) (nlevsoi)
        conc_o2_sat             => ch4_vars%conc_o2_sat_col                     , & ! Input:  [real(r8) (:,:)  ]  O2 conc in each soil layer (mol/m3) (nlevsoi)
        o2stress_unsat          => ch4_vars%o2stress_unsat_col                  , & ! Input:  [real(r8) (:,:)   ]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
        o2stress_sat            => ch4_vars%o2stress_sat_col                    , & ! Input:  [real(r8) (:,:)   ]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)

        finundated              => ch4_vars%finundated_col                        & ! Input:  [real(r8) (:)     ]  fractional inundated area (excluding dedicated wetland columns)
         )

    !! soil properties & thermohydrology

    do fc = 1, num_soilc
        c = filter_soilc(fc)

        plant_ndemand_col(c) = clm_bgc_data%plant_ndemand_col(c)
        plant_pdemand_col(c) = clm_bgc_data%plant_pdemand_col(c)

        alt_indx(c)          = clm_bgc_data%alt_indx_col(c)
        finundated(c)        = clm_bgc_data%finundated_col(c)

!        do j = 1,nlevsoi
            bd(c,:)                     = clm_bgc_data%bd_col(c,:)
            watsat(c,:)                 = clm_bgc_data%watsat_col(c,:)
            bsw(c,:)                    = clm_bgc_data%bsw_col(c,:)
            sucsat(c,:)                 = clm_bgc_data%sucsat_col(c,:)
            watfc(c,:)                  = clm_bgc_data%watfc_col(c,:)
            cellorg(c,:)                = clm_bgc_data%cellorg_col(c,:)

            soilpsi(c,:)                = clm_bgc_data%soilpsi_col(c,:)
            h2osoi_vol(c,:)             = clm_bgc_data%h2osoi_vol_col(c,:)
            h2osoi_liq(c,:)             = clm_bgc_data%h2osoi_liq_col(c,:)

            t_soisno(c,:)               = clm_bgc_data%t_soisno_col(c,:)

            o2stress_unsat(c,:)         = clm_bgc_data%o2stress_unsat_col(c,:)
            o2stress_sat(c,:)           = clm_bgc_data%o2stress_sat_col(c,:)
            o2_decomp_depth_unsat(c,:)  = clm_bgc_data%o2_decomp_depth_unsat_col(c,:)
            conc_o2_unsat(c,:)          = clm_bgc_data%conc_o2_unsat_col(c,:)
            o2_decomp_depth_sat(c,:)    = clm_bgc_data%o2_decomp_depth_sat_col(c,:)
            conc_o2_sat(c,:)            = clm_bgc_data%conc_o2_sat_col(c,:)

!        end do
    end do

    !!state variables
    do fc = 1, num_soilc
        c = filter_soilc(fc)
!        do j = 1, nlevdecomp
            do k = 1, ndecomp_pools
                decomp_cpools_vr(c,:,k) = clm_bgc_data%decomp_cpools_vr_col(c,:,k)
                decomp_npools_vr(c,:,k) = clm_bgc_data%decomp_npools_vr_col(c,:,k)
                decomp_ppools_vr(c,:,k) = clm_bgc_data%decomp_ppools_vr_col(c,:,k)
            end do

            smin_no3_vr(c,:)        = clm_bgc_data%smin_no3_vr_col(c,:)
            smin_nh4_vr(c,:)        = clm_bgc_data%smin_nh4_vr_col(c,:)
            smin_nh4sorb_vr(c,:)    = clm_bgc_data%smin_nh4sorb_vr_col(c,:)

            solutionp_vr(c,:)       = clm_bgc_data%solutionp_vr_col(c,:)
            labilep_vr(c,:)         = clm_bgc_data%labilep_vr_col(c,:)
            secondp_vr(c,:)         = clm_bgc_data%secondp_vr_col(c,:)
            sminp_vr(c,:)           = clm_bgc_data%sminp_vr_col(c,:)
            occlp_vr(c,:)           = clm_bgc_data%occlp_vr_col(c,:)
            primp_vr(c,:)           = clm_bgc_data%primp_vr_col(c,:)
!        end do
    end do

    end associate
  end subroutine clm_bgc_get_data
!!--------------------------------------------------------------------------------------

!!--------------------------------------------------------------------------------------
  ! !INTERFACE:
  !! pass data from clm_bgc to clm_bgc_data
  subroutine clm_bgc_update_data(clm_bgc_data, bounds,  &
            num_soilc, filter_soilc,                    &
            cnstate_vars, carbonflux_vars,              &
            nitrogenflux_vars, phosphorusflux_vars)

    !! USES:

    !! ARGUMENTS:
    type(bounds_type)                   , intent(in)    :: bounds
    integer                             , intent(in)    :: num_soilc          ! number of soil columns in filter
    integer                             , intent(in)    :: filter_soilc(:)    ! filter for soil columns
!    integer                             , intent(in)    :: num_soilp          ! number of soil patches in filter
!    integer                             , intent(in)    :: filter_soilp(:)    ! filter for soil patches

    type(cnstate_type)                  , intent(in)    :: cnstate_vars
    type(carbonflux_type)               , intent(in)    :: carbonflux_vars
!    type(carbonflux_type)               , intent(in)    :: c13_carbonflux_vars
!    type(carbonflux_type)               , intent(in)    :: c14_carbonflux_vars
    type(nitrogenflux_type)             , intent(in)    :: nitrogenflux_vars
    type(phosphorusflux_type)           , intent(in)    :: phosphorusflux_vars

    type(clm_bgc_interface_data_type)   , intent(inout) :: clm_bgc_data

    !! LOCAL VARIABLES:
    integer :: fc, c, j, k

    !-----------------------------------------------------------------------

    associate(                                                                                            &
         fpg                          => cnstate_vars%fpg_col                                           , & ! Output: [real(r8) (:)   ]  fraction of potential gpp (no units)
         fpi                          => cnstate_vars%fpi_col                                           , & ! Output: [real(r8) (:)   ]  fraction of potential immobilization (no units)
         fpi_vr                       => cnstate_vars%fpi_vr_col                                        , & ! Output: [real(r8) (:,:) ]  fraction of potential immobilization (no units)
         fpg_p                        => cnstate_vars%fpg_p_col                                         , & ! Output: [real(r8) (:)   ]  fraction of potential gpp (no units)
         fpi_p                        => cnstate_vars%fpi_p_col                                         , & ! Output: [real(r8) (:)   ]  fraction of potential immobilization (no units)
         fpi_p_vr                     => cnstate_vars%fpi_p_vr_col                                      , & ! Output: [real(r8) (:,:) ]  fraction of potential immobilization (no units)

         potential_immob              => nitrogenflux_vars%potential_immob_col                          , & ! Output: [real(r8) (:)   ]
         actual_immob                 => nitrogenflux_vars%actual_immob_col                             , & ! Output: [real(r8) (:)   ]
         sminn_to_plant               => nitrogenflux_vars%sminn_to_plant_col                           , & ! Output: [real(r8) (:)   ]
         sminn_to_denit_excess_vr     => nitrogenflux_vars%sminn_to_denit_excess_vr_col                 , & ! Output: [real(r8) (:,:) ]
         pot_f_nit_vr                 => nitrogenflux_vars%pot_f_nit_vr_col                             , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) potential soil nitrification flux
         pot_f_denit_vr               => nitrogenflux_vars%pot_f_denit_vr_col                           , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) potential soil denitrification flux
         f_nit_vr                     => nitrogenflux_vars%f_nit_vr_col                                 , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) soil nitrification flux
         f_denit_vr                   => nitrogenflux_vars%f_denit_vr_col                               , & ! Output: [real(r8) (:,:) ]  (gN/m3/s) soil denitrification flux
         actual_immob_no3_vr          => nitrogenflux_vars%actual_immob_no3_vr_col                      , & ! Output: [real(r8) (:,:) ]
         actual_immob_nh4_vr          => nitrogenflux_vars%actual_immob_nh4_vr_col                      , & ! Output: [real(r8) (:,:) ]
         smin_no3_to_plant_vr         => nitrogenflux_vars%smin_no3_to_plant_vr_col                     , & ! Output: [real(r8) (:,:) ]
         smin_nh4_to_plant_vr         => nitrogenflux_vars%smin_nh4_to_plant_vr_col                     , & ! Output: [real(r8) (:,:) ]
         n2_n2o_ratio_denit_vr        => nitrogenflux_vars%n2_n2o_ratio_denit_vr_col                    , & ! Output: [real(r8) (:,:) ]  ratio of N2 to N2O production by denitrification [gN/gN]
         f_n2o_denit_vr               => nitrogenflux_vars%f_n2o_denit_vr_col                           , & ! Output: [real(r8) (:,:) ]  flux of N2O from denitrification [gN/m3/s]
         f_n2o_nit_vr                 => nitrogenflux_vars%f_n2o_nit_vr_col                             , & ! Output: [real(r8) (:,:) ]  flux of N2O from nitrification [gN/m3/s]
         supplement_to_sminn_vr       => nitrogenflux_vars%supplement_to_sminn_vr_col                   , & ! Output: [real(r8) (:,:) ]
         sminn_to_plant_vr            => nitrogenflux_vars%sminn_to_plant_vr_col                        , & ! Output: [real(r8) (:,:) ]
         actual_immob_vr              => nitrogenflux_vars%actual_immob_vr_col                          , & ! Output: [real(r8) (:,:) ]

         potential_immob_p            => phosphorusflux_vars%potential_immob_p_col                      , & ! Output: [real(r8) (:)   ]
         actual_immob_p               => phosphorusflux_vars%actual_immob_p_col                         , & ! Output: [real(r8) (:)   ]
         sminp_to_plant               => phosphorusflux_vars%sminp_to_plant_col                         , & ! Output: [real(r8) (:)   ]
         supplement_to_sminp_vr       => phosphorusflux_vars%supplement_to_sminp_vr_col                 , & ! Output: [real(r8) (:,:) ]
         sminp_to_plant_vr            => phosphorusflux_vars%sminp_to_plant_vr_col                      , & ! Output: [real(r8) (:,:) ]
         actual_immob_p_vr            => phosphorusflux_vars%actual_immob_p_vr_col                      , & ! Output: [real(r8) (:,:) ]

         decomp_cascade_ntransfer_vr      =>    nitrogenflux_vars%decomp_cascade_ntransfer_vr_col       , & ! Output: [real(r8) (:,:,:) ]  vert-res transfer of N from donor to receiver pool along decomp. cascade (gN/m3/s)
         decomp_cascade_sminn_flux_vr     =>    nitrogenflux_vars%decomp_cascade_sminn_flux_vr_col      , & ! Output: [real(r8) (:,:,:) ]  vert-res mineral N flux for transition along decomposition cascade (gN/m3/s)
         potential_immob_vr               =>    nitrogenflux_vars%potential_immob_vr_col                , & ! Output: [real(r8) (:,:)   ]
         sminn_to_denit_decomp_cascade_vr =>    nitrogenflux_vars%sminn_to_denit_decomp_cascade_vr_col  , & ! Output: [real(r8) (:,:,:) ]
         gross_nmin_vr                    =>    nitrogenflux_vars%gross_nmin_vr_col                     , & ! Output: [real(r8) (:,:)   ]
         net_nmin_vr                      =>    nitrogenflux_vars%net_nmin_vr_col                       , & ! Output: [real(r8) (:,:)   ]
         gross_nmin                       =>    nitrogenflux_vars%gross_nmin_col                        , & ! Output: [real(r8) (:)     ]  gross rate of N mineralization (gN/m2/s)
         net_nmin                         =>    nitrogenflux_vars%net_nmin_col                          , & ! Output: [real(r8) (:)     ]  net rate of N mineralization (gN/m2/s)
        !!! add phosphorus
         decomp_cascade_ptransfer_vr      =>    phosphorusflux_vars%decomp_cascade_ptransfer_vr_col     , & ! Output: [real(r8) (:,:,:) ]  vert-res transfer of P from donor to receiver pool along decomp. cascade (gP/m3/s)
         decomp_cascade_sminp_flux_vr     =>    phosphorusflux_vars%decomp_cascade_sminp_flux_vr_col    , & ! Output: [real(r8) (:,:,:) ]  vert-res mineral P flux for transition along decomposition cascade (gP/m3/s)
         potential_immob_p_vr             =>    phosphorusflux_vars%potential_immob_p_vr_col            , & ! Output: [real(r8) (:,:)   ]
         gross_pmin_vr                    =>    phosphorusflux_vars%gross_pmin_vr_col                   , & ! Output: [real(r8) (:,:)   ]
         net_pmin_vr                      =>    phosphorusflux_vars%net_pmin_vr_col                     , & ! Output: [real(r8) (:,:)   ]
         gross_pmin                       =>    phosphorusflux_vars%gross_pmin_col                      , & ! Output: [real(r8) (:)     ]  gross rate of P mineralization (gP/m2/s)
         net_pmin                         =>    phosphorusflux_vars%net_pmin_col                        , & ! Output: [real(r8) (:)     ]  net rate of P mineralization (gP/m2/s)

         decomp_cascade_hr_vr             =>    carbonflux_vars%decomp_cascade_hr_vr_col                , & ! Output: [real(r8) (:,:,:) ]  vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
         decomp_cascade_ctransfer_vr      =>    carbonflux_vars%decomp_cascade_ctransfer_vr_col         , & ! Output: [real(r8) (:,:,:) ]  vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
!         decomp_k                         =>    carbonflux_vars%decomp_k_col                           , & ! Output: [real(r8) (:,:,:) ]  rate constant for decomposition (1./sec)
         phr_vr                           =>    carbonflux_vars%phr_vr_col                              , & ! Output: [real(r8) (:,:)   ]  potential HR (gC/m3/s)
         fphr                             =>    carbonflux_vars%fphr_col                                  & ! Output: [real(r8) (:,:)   ]  fraction of potential SOM + LITTER heterotrophic
         )

    !---------------------------------------------------------------------------
        clm_bgc_data%fpg_col                          = fpg
        clm_bgc_data%fpi_col                          = fpi
        clm_bgc_data%fpi_vr_col                       = fpi_vr
        clm_bgc_data%fpg_p_col                        = fpg_p
        clm_bgc_data%fpi_p_col                        = fpi_p
        clm_bgc_data%fpi_p_vr_col                     = fpi_p_vr
        clm_bgc_data%potential_immob_col              = potential_immob
        clm_bgc_data%actual_immob_col                 = actual_immob
        clm_bgc_data%sminn_to_plant_col               = sminn_to_plant
        clm_bgc_data%sminn_to_denit_excess_vr_col     = sminn_to_denit_excess_vr
        clm_bgc_data%pot_f_nit_vr_col                 = pot_f_nit_vr
        clm_bgc_data%pot_f_denit_vr_col               = pot_f_denit_vr
        clm_bgc_data%f_nit_vr_col                     = f_nit_vr
        clm_bgc_data%f_denit_vr_col                   = f_denit_vr
        clm_bgc_data%actual_immob_no3_vr_col          = actual_immob_no3_vr
        clm_bgc_data%actual_immob_nh4_vr_col          = actual_immob_nh4_vr
        clm_bgc_data%smin_no3_to_plant_vr_col         = smin_no3_to_plant_vr
        clm_bgc_data%smin_nh4_to_plant_vr_col         = smin_nh4_to_plant_vr
        clm_bgc_data%n2_n2o_ratio_denit_vr_col        = n2_n2o_ratio_denit_vr
        clm_bgc_data%f_n2o_denit_vr_col               = f_n2o_denit_vr
        clm_bgc_data%f_n2o_nit_vr_col                 = f_n2o_nit_vr
        clm_bgc_data%supplement_to_sminn_vr_col       = supplement_to_sminn_vr
        clm_bgc_data%sminn_to_plant_vr_col            = sminn_to_plant_vr
        clm_bgc_data%potential_immob_vr_col           = potential_immob_vr
        clm_bgc_data%actual_immob_vr_col              = actual_immob_vr
        clm_bgc_data%potential_immob_p_col            = potential_immob_p
        clm_bgc_data%actual_immob_p_col               = actual_immob_p
        clm_bgc_data%sminp_to_plant_col               = sminp_to_plant
        clm_bgc_data%supplement_to_sminp_vr_col       = supplement_to_sminp_vr
        clm_bgc_data%sminp_to_plant_vr_col            = sminp_to_plant_vr
        clm_bgc_data%potential_immob_p_vr_col         = potential_immob_p_vr
        clm_bgc_data%actual_immob_p_vr_col            = actual_immob_p_vr

        clm_bgc_data%decomp_cascade_ntransfer_vr_col      = decomp_cascade_ntransfer_vr
        clm_bgc_data%decomp_cascade_sminn_flux_vr_col     = decomp_cascade_sminn_flux_vr
        clm_bgc_data%potential_immob_vr_col               = potential_immob_vr
        clm_bgc_data%sminn_to_denit_decomp_cascade_vr_col = sminn_to_denit_decomp_cascade_vr
        clm_bgc_data%gross_nmin_vr_col                    = gross_nmin_vr
        clm_bgc_data%net_nmin_vr_col                      = net_nmin_vr

        !! phosphorus
        clm_bgc_data%decomp_cascade_ptransfer_vr_col      = decomp_cascade_ptransfer_vr
        clm_bgc_data%decomp_cascade_sminp_flux_vr_col     = decomp_cascade_sminp_flux_vr
        clm_bgc_data%potential_immob_p_vr_col             = potential_immob_p_vr
        clm_bgc_data%gross_pmin_vr_col                    = gross_pmin_vr
        clm_bgc_data%net_pmin_vr_col                      = net_pmin_vr

        clm_bgc_data%decomp_cascade_hr_vr_col             = decomp_cascade_hr_vr
        clm_bgc_data%decomp_cascade_ctransfer_vr_col      = decomp_cascade_ctransfer_vr

        clm_bgc_data%phr_vr_col                           = phr_vr
        clm_bgc_data%fphr_col                             = fphr


    end associate
  end subroutine clm_bgc_update_data
!!--------------------------------------------------------------------------------------

!!--------------------------------------------------------------------------------------
  subroutine update_bgc_data_clm2clm(clm_bgc_data,bounds,         &
           num_soilc, filter_soilc,                               &
           num_soilp, filter_soilp,                               &
           atm2lnd_vars,                                          &
           waterstate_vars, waterflux_vars,                       &
           soilstate_vars,  temperature_vars, energyflux_vars,    &
           soilhydrology_vars, soil_water_retention_curve,        &
           cnstate_vars, carbonflux_vars, carbonstate_vars,       &
           nitrogenflux_vars, nitrogenstate_vars,                 &
           phosphorusflux_vars, phosphorusstate_vars,             &
           ch4_vars)
    !! USES


    implicit none

    ! !ARGUMENTS:
    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                     , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                     , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                     , intent(in)    :: filter_soilp(:)   ! filter for soil patches
    type(atm2lnd_type)          , intent(in)    :: atm2lnd_vars
    type(waterstate_type)       , intent(inout) :: waterstate_vars
    type(waterflux_type)        , intent(inout) :: waterflux_vars
    type(soilstate_type)        , intent(inout) :: soilstate_vars
    type(temperature_type)      , intent(inout) :: temperature_vars
    type(soilhydrology_type)    , intent(inout) :: soilhydrology_vars
    type(energyflux_type)       , intent(inout) :: energyflux_vars

    type(cnstate_type)          , intent(inout) :: cnstate_vars
    type(carbonflux_type)       , intent(inout) :: carbonflux_vars
    type(carbonstate_type)      , intent(inout) :: carbonstate_vars
    type(nitrogenflux_type)     , intent(inout) :: nitrogenflux_vars
    type(nitrogenstate_type)    , intent(inout) :: nitrogenstate_vars
    type(phosphorusflux_type)   , intent(inout) :: phosphorusflux_vars
    type(phosphorusstate_type)  , intent(inout) :: phosphorusstate_vars
    type(ch4_type)              , intent(inout) :: ch4_vars

    class(soil_water_retention_curve_type)  , intent(in) :: soil_water_retention_curve
    type(clm_bgc_interface_data_type)       , intent(in) :: clm_bgc_data

    !-----------------------------------------------------------------------

    character(len=256) :: subname = "get_clm_bgc_data"

    !! bgc_state_decomp is updated in CLM
    !! by passing bgc_flux_decomp_sourcesink into CNSoilLittVertTransp
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
!!--------------------------------------------------------------------------------------
! END of CLM-bgc through interface
!!--------------------------------------------------------------------------------------

end module clm_bgc_interfaceMod

