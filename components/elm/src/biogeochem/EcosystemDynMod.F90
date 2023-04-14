module EcosystemDynMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Ecosystem dynamics: phenology, vegetation
  !
  ! !USES:
  use dynSubgridControlMod, only : get_do_harvest
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_sys_mod         , only : shr_sys_flush
  use elm_varctl          , only : use_c13, use_c14, use_fates, use_dynroot
  use decompMod           , only : bounds_type
  use perf_mod            , only : t_startf, t_stopf
  use spmdMod             , only : masterproc
  use elm_varctl          , only : use_century_decomp
  use elm_varctl          , only : ero_ccycle 
  use CNStateType         , only : cnstate_type
  use CanopyStateType     , only : canopystate_type
  use SoilStateType       , only : soilstate_type
  use atm2lndType         , only : atm2lnd_type
  use SoilStateType       , only : soilstate_type
  use CanopyStateType     , only : canopystate_type
  use PhotosynthesisType  , only : photosyns_type
  use CH4Mod              , only : ch4_type
  use EnergyFluxType      , only : energyflux_type
  use SoilHydrologyType   , only : soilhydrology_type
  use FrictionVelocityType, only : frictionvel_type
  use SedFluxType         , only : sedflux_type
  use ColumnDataType      , only : col_cs, c13_col_cs, c14_col_cs
  use ColumnDataType      , only : col_cf, c13_col_cf, c14_col_cf
  use ColumnDataType      , only : col_ns, col_nf
  use ColumnDataType      , only : col_ps, col_pf
  use VegetationDataType  , only : veg_cs, c13_veg_cs, c14_veg_cs
  use VegetationDataType  , only : veg_cf, c13_veg_cf, c14_veg_cf
  use VegetationDataType  , only : veg_ns, veg_nf
  use VegetationDataType  , only : veg_ps, veg_pf
  use ELMFatesInterfaceMod  , only : hlm_fates_interface_type
  use elm_instMod  , only : alm_fates
  ! bgc interface & pflotran
  use elm_varctl         , only : use_elm_interface, use_elm_bgc, use_pflotran, pf_cmode, pf_hmode
  use VerticalProfileMod , only : decomp_vertprofiles
  use AllocationMod      , only : nu_com_nfix, nu_com_phosphatase
  use elm_varctl         , only : nu_com, use_pheno_flux_limiter
  use PhenologyFLuxLimitMod , only : phenology_flux_limiter, InitPhenoFluxLimiter

  use timeinfoMod
  use perfMod_GPU
  use VegetationDataType , only : veg_cf_summary, veg_cf_summary_for_ch4, veg_cf_summary_rr
  use VegetationDataType , only : veg_nf_summary, veg_ns_summary, veg_cs_Summary
  use VegetationDataType , only : veg_pf_summary, veg_ps_summary
  use VegetationDataType , only : veg_cf_setvalues, veg_nf_setvalues, veg_pf_setvalues

  use ColumnDataType , only : col_cf_summary, col_nf_summary, col_pf_Summary
  use ColumnDataType , only : col_cs_summary, col_ns_summary, col_ps_summary
  use ColumnDataType , only : col_cf_summary_for_ch4
  use ColumnDataType , only : col_cf_setvalues, col_nf_setvalues, col_pf_setvalues 


  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: EcosystemDynInit          ! Ecosystem dynamics initialization
  public :: EcosystemDynLeaching      ! Ecosystem dynamics: phenology, vegetation, doing N leaching
  !----------------------------------------------------------------------
  ! bgc&th interface & pflotran:
  ! EcosystemDynNoLeaching is divided into 2 subroutines:
  public :: EcosystemDynNoLeaching1   ! Ecosystem dynamics: phenology, vegetation, before doing soil_bgc
  public :: EcosystemDynNoLeaching2   ! Ecosystem dynamics: phenology, vegetation, after doing soil_bgc & before doing N leaching
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine EcosystemDynInit(bounds, elm_fates)
    !
    ! !DESCRIPTION:
    ! Initialzation of the CN Ecosystem dynamics.
    !
    ! !USES:
    use AllocationMod, only : AllocationInit
    use PhenologyMod , only : PhenologyInit
    use FireMod      , only : FireInit
    use C14DecayMod  , only : C14_init_BombSpike
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds
    type(hlm_fates_interface_type), intent(in) :: elm_fates
    !-----------------------------------------------------------------------

    call AllocationInit (bounds, elm_fates)

    if(use_fates) return

    call PhenologyInit  (bounds)
    call FireInit       (bounds)

    if ( use_c14 ) then
       call C14_init_BombSpike()
    end if

    if(use_pheno_flux_limiter)then
        call InitPhenoFluxLimiter()
    endif

  end subroutine EcosystemDynInit


  !-----------------------------------------------------------------------

  subroutine EcosystemDynLeaching(bounds, num_soilc, filter_soilc, &
       num_soilp, filter_soilp, num_pcropp, filter_pcropp, doalb, &
       cnstate_vars,  &
       frictionvel_vars, canopystate_vars )
    !
    ! !DESCRIPTION:
    ! The core CN code is executed here. Calculates fluxes for maintenance
    ! respiration, decomposition, allocation, phenology, and growth respiration.
    ! These routines happen on the radiation time step so that canopy structure
    ! stays synchronized with albedo calculations.
    !
    ! !USES:
    use spmdMod              , only: masterproc
    use PhosphorusDynamicsMod         , only: PhosphorusWeathering,PhosphorusAdsportion,PhosphorusDesoprtion,PhosphorusOcclusion
    use PhosphorusDynamicsMod         , only: PhosphorusBiochemMin,PhosphorusLeaching
    use NitrogenDynamicsMod       , only: NitrogenLeaching
    use NitrogenStateUpdate3Mod   , only: NitrogenStateUpdate3
    use PhosphorusStateUpdate3Mod     , only: PhosphorusStateUpdate3
    use PrecisionControlMod  , only: PrecisionControl
    use perf_mod             , only: t_startf, t_stopf
    use shr_sys_mod          , only: shr_sys_flush
    use PhosphorusDynamicsMod         , only: PhosphorusBiochemMin_balance

    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                  , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:)   ! filter for soil patches
    integer                  , intent(in)    :: num_pcropp        ! number of prog. crop patches in filter
    integer                  , intent(in)    :: filter_pcropp(:)  ! filter for prognostic crop patches
    logical                  , intent(in)    :: doalb             ! true = surface albedo calculation time step
    type(cnstate_type)       , intent(inout) :: cnstate_vars
    type(frictionvel_type)   , intent(in)    :: frictionvel_vars
    type(canopystate_type)   , intent(inout) :: canopystate_vars

    character(len=64) :: event
    real(r8) :: dt
    !-----------------------------------------------------------------------
    dt = dtime_mod;

    ! only do if ed is off
    event = 'PhosphorusWeathering'
    call t_start_lnd(event)
    call PhosphorusWeathering(num_soilc, filter_soilc, &
         cnstate_vars, dt)
    call t_stop_lnd(event)

    event = 'PhosphorusAdsportion'
    call t_start_lnd(event)
    call PhosphorusAdsportion(num_soilc, filter_soilc, &
         cnstate_vars, dt)
    call t_stop_lnd(event)

    event = 'PhosphorusDesoprtion'
    call t_start_lnd(event)
    call PhosphorusDesoprtion(num_soilc, filter_soilc, &
         cnstate_vars, dt)
    call t_stop_lnd(event)

    event = 'PhosphorusOcclusion'
    call t_start_lnd(event)
    call PhosphorusOcclusion(num_soilc, filter_soilc, &
         cnstate_vars, dt)
    call t_stop_lnd(event)

    if (.not. nu_com_phosphatase) then
      event = 'PhosphorusBiochemMin'
       call t_start_lnd(event)
       call PhosphorusBiochemMin(bounds,num_soilc, filter_soilc, &
            cnstate_vars, dt)
       call t_stop_lnd(event)
    else
      ! nu_com_phosphatase is true
      !call t_startf('PhosphorusBiochemMin')
      !call PhosphorusBiochemMin_balance(bounds,num_soilc, filter_soilc, &
      !     cnstate_vars,nitrogenstate_vars,phosphorusstate_vars,phosphorusflux_vars)
      !call t_stopf('PhosphorusBiochemMin')
    end if

    !-----------------------------------------------------------------------
    ! pflotran: when both 'pf-bgc' and 'pf-h' on, no need to call CLM-CN's N leaching module
    if (.not. (pf_cmode .and. pf_hmode)) then
     call NitrogenLeaching(bounds, num_soilc, filter_soilc, dt)

     call PhosphorusLeaching(bounds, num_soilc, filter_soilc, dt)
    end if !(.not. (pf_cmode .and. pf_hmode))
       !-----------------------------------------------------------------------

    event = 'CNUpdate3'
    call t_start_lnd(event)
    call NitrogenStateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp,dt)
    call t_stop_lnd(event)

    event = 'PUpdate3'
    call t_start_lnd(event)
    call PhosphorusStateUpdate3(bounds,num_soilc, filter_soilc, num_soilp, filter_soilp, &
        cnstate_vars, dt)
    call t_stop_lnd(event)

    event = 'CNPsum'
    call t_start_lnd(event)
    call PrecisionControl(num_soilc, filter_soilc, num_soilp, filter_soilp  )

    call col_cf_summary_for_ch4(col_cf,bounds, num_soilc, filter_soilc)

    ! Only update the veg_ data structures if we are using cn
    if(.not.use_fates) then
       call veg_cf_summary_for_ch4(veg_cf,bounds, num_soilp, filter_soilp)
       call veg_cf_Summary(veg_cf,bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, 'bulk', col_cf)
       if ( use_c13 ) then
          call veg_cf_Summary(veg_cf,bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, 'c13', c13_col_cf)
          call col_cf_Summary(col_cf,bounds, num_soilc, filter_soilc, 'c13')
       end if
       if ( use_c14 ) then
          call veg_cf_Summary(veg_cf,bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, 'c14', c14_col_cf)
          call col_cf_Summary(col_cf,bounds, num_soilc, filter_soilc, 'c14')
       end if
       call veg_cs_Summary(veg_cs,bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_cs)
       if ( use_c13 ) then
          call veg_cs_Summary(c13_veg_cs,bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, c13_col_cs)
          call col_cs_Summary(c13_col_cs,bounds, num_soilc, filter_soilc)
       end if
       if ( use_c14 ) then
          call veg_cs_Summary(c14_veg_cs,bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, c14_col_cs)
          call col_cs_Summary(c14_col_cs,bounds, num_soilc, filter_soilc)
       end if
       call veg_nf_Summary(veg_nf, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_nf)
       call veg_ns_Summary(veg_ns, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_ns)
       call veg_pf_Summary(veg_pf, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_pf)
       call veg_ps_Summary(veg_ps, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_ps)
    else
       ! In this scenario, we simply zero all of the
       ! column level variables that would had been upscaled
       ! in the veg summary with p2c
       call col_cs%ZeroForFates(bounds,num_soilc, filter_soilc)
       call col_ns%ZeroForFates(bounds,num_soilc, filter_soilc)
       call col_ps%ZeroForFates(bounds,num_soilc, filter_soilc)
       call col_cf%ZeroForFates(bounds,num_soilc, filter_soilc)
       call col_nf%ZeroForFates(bounds,num_soilc, filter_soilc)
       call col_pf%ZeroForFates(bounds,num_soilc, filter_soilc)
    end if

    call col_cf_Summary(col_cf,bounds, num_soilc, filter_soilc, 'bulk')
    call col_cs_Summary(col_cs,bounds, num_soilc, filter_soilc)

    call col_nf_Summary(col_nf,bounds, num_soilc, filter_soilc)
    call col_ns_Summary(col_ns,bounds, num_soilc, filter_soilc)

    call col_pf_Summary(col_pf,bounds, num_soilc, filter_soilc)
    call col_ps_Summary(col_ps,bounds, num_soilc, filter_soilc)

    call t_stop_lnd(event)

  end subroutine EcosystemDynLeaching


!-------------------------------------------------------------------------------------------------
  subroutine EcosystemDynNoLeaching1(bounds,                          &
       num_soilc, filter_soilc,                                         &
       num_soilp, filter_soilp,                                         &
       cnstate_vars,  atm2lnd_vars,           &
       canopystate_vars, soilstate_vars, crop_vars,   &
       ch4_vars, photosyns_vars     )
    !-------------------------------------------------------------------
    ! bgc interface
    ! Phase-1 of EcosystemDynNoLeaching
    ! call Allocation1_PlantNPDemand before soil_bgc
    !-------------------------------------------------------------------

    ! !DESCRIPTION:
    ! The core CN code is executed here. Calculates fluxes for maintenance
    ! respiration, decomposition, allocation, phenology, and growth respiration.
    ! These routines happen on the radiation time step so that canopy structure
    ! stays synchronized with albedo calculations.
    !
    ! !USES:
    use NitrogenDynamicsMod         , only: NitrogenDeposition,NitrogenFixation, NitrogenFert, CNSoyfix
    use PhosphorusDynamicsMod       , only: PhosphorusDeposition
    use MaintenanceRespMod          , only: MaintenanceResp
    use DecompCascadeBGCMod   , only: decomp_rate_constants_bgc
    use DecompCascadeCNMod    , only: decomp_rate_constants_cn
    use CropType              , only: crop_type
    use elm_varpar            , only: crop_prog
    use AllocationMod         , only: Allocation1_PlantNPDemand ! Phase-1 of CNAllocation
    use NitrogenDynamicsMod   , only: NitrogenLeaching
    use PhosphorusDynamicsMod       , only: PhosphorusLeaching
    use NitrogenDynamicsMod         , only: NitrogenFixation_balance
    use PhosphorusDynamicsMod       , only: PhosphorusWeathering,PhosphorusAdsportion,PhosphorusDesoprtion,PhosphorusOcclusion
    use PhosphorusDynamicsMod       , only: PhosphorusBiochemMin,PhosphorusBiochemMin_balance

    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                  , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:)   ! filter for soil patches
    type(cnstate_type)       , intent(inout) :: cnstate_vars
    type(atm2lnd_type)       , intent(in)    :: atm2lnd_vars
    type(canopystate_type)   , intent(in)    :: canopystate_vars
    type(soilstate_type)     , intent(inout) :: soilstate_vars
    type(crop_type)          , intent(in)    :: crop_vars
    type(ch4_type)           , intent(in)    :: ch4_vars
    type(photosyns_type)     , intent(in)    :: photosyns_vars

    character(len=64) :: event
    real(r8) :: dt, dayspyr
    integer  :: year, mon, day, sec
    !-----------------------------------------------------------------------
    dt = dtime_mod;
    year = year_curr; mon = mon_curr; day = day_curr; sec= secs_curr
    dayspyr = dayspyr_mod
    !-----------------------------------------------------------------------

    ! Call the main CN routines

    ! --------------------------------------------------
    ! zero the C and N fluxes
    ! --------------------------------------------------

    event = 'CNZero'
    call t_start_lnd(event)

    if(.not.use_fates) then
       call veg_cf_SetValues(veg_cf, num_soilp, filter_soilp, 0._r8)
       if ( use_c13 ) then
          call veg_cf_setvalues(c13_veg_cf,num_soilp, filter_soilp, 0._r8)
          call col_cf_SetValues(c13_col_cf,num_soilc, filter_soilc, 0._r8)
       end if
       if ( use_c14 ) then
          call veg_cf_setvalues(c14_veg_cf,num_soilp, filter_soilp, 0._r8)
          call col_cf_SetValues(c14_col_cf,num_soilc, filter_soilc, 0._r8)
       end if
       call veg_nf_setvalues (veg_nf,num_soilp, filter_soilp, 0._r8)
       call veg_pf_SetValues (veg_pf,num_soilp, filter_soilp, 0._r8)
    end if

    call col_cf_SetValues(col_cf,num_soilc, filter_soilc, 0._r8)
    call col_nf_SetValues(col_nf,num_soilc, filter_soilc, 0._r8)
    call col_pf_SetValues(col_pf,num_soilc, filter_soilc, 0._r8)

    call t_stop_lnd(event)

    ! --------------------------------------------------
    ! Nitrogen Deposition, Fixation and Respiration, phosphorus dynamics
    ! --------------------------------------------------

    event = 'CNDeposition'
    call t_start_lnd(event)
    call NitrogenDeposition(bounds, &
         atm2lnd_vars, dt )
    call t_stop_lnd(event)

    event = 'CNFixation'
    if ( (.not. nu_com_nfix) .or. use_fates) then
      call t_start_lnd(event)
       call NitrogenFixation( num_soilc, filter_soilc, dayspyr)
       call t_stop_lnd(event)
    else
       ! nu_com_nfix is true
       call t_start_lnd(event)
       call NitrogenFixation_balance( num_soilc, filter_soilc, cnstate_vars )
       call t_stop_lnd(event)
    end if

   if(.not.use_fates)then

       event = 'MaintenanceResp'
       call t_start_lnd(event)
       if (crop_prog) then
          call NitrogenFert(bounds, num_soilc,filter_soilc )

          call CNSoyfix(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
                         crop_vars, cnstate_vars )
       end if
       ! This is auto-trophic respiration, thus don't call this for FATES
       call MaintenanceResp(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            canopystate_vars, soilstate_vars,  photosyns_vars )
       call t_stop_lnd(event)

    end if


    if ( nu_com .ne. 'RD') then
       ! for P competition purpose, calculate P fluxes that will potentially increase solution P pool
       ! then competitors take up solution P
       event ='PhosphorusWeathering'
       call t_start_lnd(event)
       call PhosphorusWeathering(num_soilc, filter_soilc, cnstate_vars, dt)
       call t_stop_lnd(event)

       event = 'PhosphorusBiochemMin'
       if (.not. nu_com_phosphatase) then
           call t_start_lnd(event)
           call PhosphorusBiochemMin(bounds,num_soilc, filter_soilc, &
                cnstate_vars, dt)
           call t_stop_lnd(event)
       else
           ! nu_com_phosphatase is true
           call t_start_lnd(event)
           call PhosphorusBiochemMin_balance(bounds,num_soilc, filter_soilc, &
                cnstate_vars, dt)
           call t_stop_lnd(event)
       end if
    end if

    ! --------------------------------------------------
    ! Phosphorus Deposition ! X.SHI
    ! --------------------------------------------------

    event = 'PhosphorusDeposition'
    call t_start_lnd(event)
    call PhosphorusDeposition(bounds,  atm2lnd_vars )
    call t_stop_lnd(event)

    !-------------------------------------------------------------------------------------------------
    ! plfotran: 'decomp_rate_constants' must be calculated before entering "clm_interface"
    if (use_century_decomp) then
      call decomp_rate_constants_bgc(bounds, num_soilc, filter_soilc, &
           canopystate_vars, soilstate_vars, ch4_vars, cnstate_vars)
    else
      call decomp_rate_constants_cn(bounds, num_soilc, filter_soilc, &
           canopystate_vars, soilstate_vars,  ch4_vars, cnstate_vars)
    end if

    !-------------------------------------------------------------------------------------------------
    ! 'decomp_vertprofiles' (calc nfixation_prof) is moved from SoilLittDecompAlloc:
    ! ------------------------------------------------------------------------------------------------
    call decomp_vertprofiles(bounds,                        &
          num_soilc, filter_soilc, num_soilp, filter_soilp, &
          soilstate_vars, canopystate_vars, cnstate_vars)

    !-------------------------------------------------------------------------------------------------
    ! Allocation1 is always called (w/ or w/o use_elm_interface)
    ! pflotran: call 'Allocation1' to obtain potential N demand for support initial GPP
    if(.not.use_fates)then

       event = 'CNAllocation - phase-1'
       call t_start_lnd(event)
       call Allocation1_PlantNPDemand (bounds                             , &
                num_soilc, filter_soilc, num_soilp, filter_soilp            , &
                photosyns_vars, crop_vars, canopystate_vars, cnstate_vars   , &
               dt, year )
       call t_stop_lnd(event)
    end if


  end subroutine EcosystemDynNoLeaching1

!-------------------------------------------------------------------------------------------------
  subroutine EcosystemDynNoLeaching2(bounds,                                  &
       num_soilc, filter_soilc,                                                 &
       num_soilp, filter_soilp, num_pcropp, filter_pcropp, doalb,               &
       num_ppercropp, filter_ppercropp,                                         &
       cnstate_vars,                                              &
       atm2lnd_vars,               &
       canopystate_vars, soilstate_vars,  crop_vars, ch4_vars, &
       photosyns_vars, soilhydrology_vars, energyflux_vars,          &
       sedflux_vars)
    !-------------------------------------------------------------------
    ! bgc interface
    ! Phase-2 of EcosystemDynNoLeaching
    ! call SoilLittDecompAlloc (w/o bgc_interface) & SoilLittDecompAlloc2
    !-------------------------------------------------------------------

    ! !DESCRIPTION:
    ! The core CN code is executed here. Calculates fluxes for maintenance
    ! respiration, decomposition, allocation, phenology, and growth respiration.
    ! These routines happen on the radiation time step so that canopy structure
    ! stays synchronized with albedo calculations.
    !
    ! !USES:
    use PhenologyMod         , only: Phenology, CNLitterToColumn
    use GrowthRespMod             , only: GrowthResp
    use CarbonStateUpdate1Mod     , only: CarbonStateUpdate1,CarbonStateUpdate0
    use NitrogenStateUpdate1Mod     , only: NitrogenStateUpdate1
    use PhosphorusStateUpdate1Mod       , only: PhosphorusStateUpdate1
    use GapMortalityMod        , only: GapMortality
    use CarbonStateUpdate2Mod     , only: CarbonStateUpdate2, CarbonStateUpdate2h
    use NitrogenStateUpdate2Mod     , only: NitrogenStateUpdate2, NitrogenStateUpdate2h
    use PhosphorusStateUpdate2Mod       , only: PhosphorusStateUpdate2, PhosphorusStateUpdate2h
    use FireMod              , only: FireArea, FireFluxes
    use ErosionMod           , only: ErosionFluxes
    use CarbonStateUpdate3Mod     , only: CarbonStateUpdate3
    use CarbonIsoFluxMod          , only: CarbonIsoFlux1, CarbonIsoFlux2, CarbonIsoFlux2h, CarbonIsoFlux3
    use C14DecayMod          , only: C14Decay, C14BombSpike
    use WoodProductsMod      , only: WoodProducts
    use CropHarvestPoolsMod  , only: CropHarvestPools
    use SoilLittVertTranspMod, only: SoilLittVertTransp
    use CropType               , only: crop_type
    use dynHarvestMod          , only: CNHarvest
    use RootDynamicsMod        , only: RootDynamics
    use SoilLittDecompMod            , only: SoilLittDecompAlloc
    use SoilLittDecompMod            , only: SoilLittDecompAlloc2 !after SoilLittDecompAlloc

    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                  , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:)   ! filter for soil patches
    integer                  , intent(in)    :: num_pcropp        ! number of prog. crop patches in filter
    integer                  , intent(in)    :: filter_pcropp(:)  ! filter for prognostic crop patches
    integer                  , intent(in)    :: num_ppercropp     ! number of prog perennial crop patches in filter
    integer                  , intent(in)    :: filter_ppercropp(:) ! filter for prognostic perennial crop patches
    logical                  , intent(in)    :: doalb             ! true = surface albedo calculation time step
    type(cnstate_type)       , intent(inout) :: cnstate_vars
    type(atm2lnd_type)       , intent(in)    :: atm2lnd_vars
    type(canopystate_type)   , intent(in)    :: canopystate_vars
    type(soilstate_type)     , intent(inout) :: soilstate_vars
    type(crop_type)          , intent(inout) :: crop_vars
    type(ch4_type)           , intent(in)    :: ch4_vars
    type(photosyns_type)     , intent(in)    :: photosyns_vars
    type(soilhydrology_type) , intent(in)    :: soilhydrology_vars
    type(energyflux_type)    , intent(in)    :: energyflux_vars
!
    type(sedflux_type)       , intent(in)    :: sedflux_vars
    character(len=64) :: event
    real(r8) :: dt
    integer :: c13, c14
    c13 = 0
    c14 = 1
    !-----------------------------------------------------------------------
    dt = dtime_mod
    ! Call the main CN routines

    event = 'SoilLittDecompAlloc'
    call t_start_lnd(event)
    !----------------------------------------------------------------
    if(.not.use_elm_interface) then
       ! directly run elm-bgc
       ! if (use_elm_interface & use_elm_bgc), then CNDecomAlloc is called in elm_driver
       call SoilLittDecompAlloc (bounds, num_soilc, filter_soilc,    &
                  num_soilp, filter_soilp,                     &
                  canopystate_vars, soilstate_vars,            &
                  cnstate_vars, ch4_vars,                      &
                  dt)
    end if !if(.not.use_elm_interface)

    call t_stopf('SoilLittDecompAlloc')

    event = 'SoilLittDecompAlloc2'
    call t_start_lnd(event)
    !----------------------------------------------------------------
    ! SoilLittDecompAlloc2 is called by both elm-bgc & pflotran
    ! pflotran: call 'SoilLittDecompAlloc2' to calculate some diagnostic variables and 'fpg' for plant N uptake
    ! pflotran & elm-bgc : 'Allocation3_AG' and vertically integrate net and gross mineralization fluxes
    call SoilLittDecompAlloc2 (bounds, num_soilc, filter_soilc, num_soilp, filter_soilp,           &
             photosyns_vars, canopystate_vars, soilstate_vars,         &
             cnstate_vars, ch4_vars,                  &
             crop_vars, atm2lnd_vars,                 &
             dt )
    call t_stop_lnd(event)

    !----------------------------------------------------------------

    if(.not.use_fates)then

        !--------------------------------------------
        ! Phenology
        !--------------------------------------------

        ! Phenology needs to be called after SoilLittDecompAlloc, because it
        ! depends on current time-step fluxes to new growth on the last
        ! litterfall timestep in deciduous systems

        event = 'Phenology'
        call t_start_lnd(event)
        call Phenology(num_soilc, filter_soilc, num_soilp, filter_soilp, &
             num_pcropp, filter_pcropp, num_ppercropp, filter_ppercropp, doalb, atm2lnd_vars, &
             crop_vars, canopystate_vars, soilstate_vars, &
             cnstate_vars )
        call t_stop_lnd(event)

        !--------------------------------------------
        ! Growth respiration
        !--------------------------------------------

        event = 'GrowthResp'
        call t_start_lnd(event)
        call GrowthResp(num_soilp, filter_soilp )
        call t_stop_lnd(event)

        call veg_cf_summary_rr(veg_cf,bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, col_cf)
        if(use_c13) then
           call veg_cf_summary_rr(c13_veg_cf,bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, c13_col_cf)
        endif
        if(use_c14) then
           call veg_cf_summary_rr(c14_veg_cf,bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, c14_col_cf)
        endif


        !--------------------------------------------
        ! Dynamic Roots
        !--------------------------------------------

        if( use_dynroot ) then
            event = 'RootDynamics'
            call t_start_lnd(event)

            call RootDynamics(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
                  canopystate_vars,   &
                  cnstate_vars, crop_vars, energyflux_vars, soilstate_vars)
            call t_stop_lnd(event)
        end if

        !--------------------------------------------
        ! CNUpdate0
        !--------------------------------------------

       event = 'CNUpdate0'
       call t_start_lnd(event)
       call CarbonStateUpdate0(num_soilp, filter_soilp,veg_cs,veg_cf, dt)
       if ( use_c13 ) then
          call CarbonStateUpdate0(num_soilp, filter_soilp,c13_veg_cs,c13_veg_cf, dt)
       end if
       if ( use_c14 ) then
          call CarbonStateUpdate0(num_soilp, filter_soilp,c14_veg_cs,c14_veg_cf ,dt)
       end if
       call t_stop_lnd(event)

        !--------------------------------------------
        if(use_pheno_flux_limiter)then
          event = 'phenology_flux_limiter'
          call t_start_lnd(event)
          call phenology_flux_limiter(bounds, num_soilc, filter_soilc,&
            num_soilp, filter_soilp, crop_vars, cnstate_vars,  &
            veg_cf, veg_cs, &
            c13_veg_cf, c13_veg_cs, &
            c14_veg_cf, c14_veg_cs, &
            veg_nf, veg_ns, veg_pf, veg_ps)
          call t_stop_lnd(event)
        endif
        event = 'CNLitterToColumn'
        call t_start_lnd(event)
        call CNLitterToColumn(num_soilp, filter_soilp, cnstate_vars )
        call t_stop_lnd(event)

        !--------------------------------------------
        ! Update1
        !--------------------------------------------
        event = 'CNUpdate1'
        call t_start_lnd(event)

        if ( use_c13 ) then
          call CarbonIsoFlux1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, &
               isotope=c13, isocol_cs=c13_col_cs, isoveg_cs=c13_veg_cs, isocol_cf=c13_col_cf, isoveg_cf=c13_veg_cf)
        end if

        if ( use_c14 ) then
          call CarbonIsoFlux1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars , &
               isotope=c14, isocol_cs=c14_col_cs, isoveg_cs=c14_veg_cs, isocol_cf=c14_col_cf, isoveg_cf=c14_veg_cf)
         end if
         call t_stop_lnd(event) 

    end if  ! if(.not.use_fates)

    if(use_fates) then
       ! In this scenario, we simply zero all of the
       ! column level variables that would had been upscaled
       ! in the veg summary with p2c
       call col_cf%ZeroForFatesRR(bounds,num_soilc, filter_soilc)

       ! Transfer fates litter fluxes into ELM source arrays
       call alm_fates%UpdateLitterFluxes(bounds)
    end if
    
   event = 'CNUpdate1'
   call t_start_lnd(event)

   call CarbonStateUpdate1(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
         crop_vars, col_cs, veg_cs, col_cf, veg_cf, dt)

   if ( use_c13 ) then
      call CarbonStateUpdate1(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
           crop_vars, c13_col_cs, c13_veg_cs, c13_col_cf, c13_veg_cf,dt)
   end if
   if ( use_c14 ) then
      call CarbonStateUpdate1(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
           crop_vars, c14_col_cs, c14_veg_cs, c14_col_cf, c14_veg_cf,dt)
   end if

   call NitrogenStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
        cnstate_vars, dt)

   call PhosphorusStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
        cnstate_vars, dt)


   call t_stop_lnd(event)

   event = 'SoilLittVertTransp'
   call t_start_lnd(event)
   call SoilLittVertTransp(bounds, &
            num_soilc, filter_soilc, &
            canopystate_vars, cnstate_vars )
       call t_stop_lnd(event)
   if(.not.use_fates)then
       event = 'CNGapMortality'
       call t_start_lnd(event)
       call GapMortality( num_soilc, filter_soilc, num_soilp, filter_soilp,&
                              cnstate_vars, crop_vars)
       call t_stop_lnd(event)

       !--------------------------------------------
       ! Update2
       !--------------------------------------------
       event = 'CNUpdate2'
       call t_start_lnd(event)

       if ( use_c13 ) then
           call CarbonIsoFlux2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
                cnstate_vars,  &
                isotope=c13, isocol_cs=c13_col_cs, isoveg_cs=c13_veg_cs, isocol_cf=c13_col_cf, isoveg_cf=c13_veg_cf)
       end if

       if ( use_c14 ) then
           call CarbonIsoFlux2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
                cnstate_vars, &
                isotope=c14, isocol_cs=c14_col_cs, isoveg_cs=c14_veg_cs, isocol_cf=c14_col_cf, isoveg_cf=c14_veg_cf)
       end if



       call CarbonStateUpdate2( num_soilc, filter_soilc, num_soilp, filter_soilp, &
            col_cs, veg_cs, col_cf, veg_cf, dt)

       if ( use_c13 ) then
          call CarbonStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
                c13_col_cs, c13_veg_cs, c13_col_cf, c13_veg_cf,dt)
       end if
       if ( use_c14 ) then
          call CarbonStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c14_col_cs, c14_veg_cs, c14_col_cf, c14_veg_cf, dt)
       end if
       call NitrogenStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            dt )

       call PhosphorusStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            dt)

       if (get_do_harvest()) then
          call CNHarvest(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars )
       end if

       if ( use_c13 ) then
          call CarbonIsoFlux2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, &
               isotope=c13, isocol_cs=c13_col_cs, isoveg_cs=c13_veg_cs, isocol_cf=c13_col_cf, isoveg_cf=c13_veg_cf)
       end if
       if ( use_c14 ) then
          call CarbonIsoFlux2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars , &
               isotope=c14, isocol_cs=c14_col_cs, isoveg_cs=c14_veg_cs, isocol_cf=c14_col_cf, isoveg_cf=c14_veg_cf)
       end if



       call CarbonStateUpdate2h( num_soilc, filter_soilc,  num_soilp, filter_soilp, &
             col_cs, veg_cs, col_cf, veg_cf, dt)
       if ( use_c13 ) then
          call CarbonStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c13_col_cs, c13_veg_cs, c13_col_cf, c13_veg_cf, dt)
       end if
       if ( use_c14 ) then
          call CarbonStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c14_col_cs, c14_veg_cs, c14_col_cf, c14_veg_cf, dt)
       end if

       call NitrogenStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            dt)

       call PhosphorusStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            dt)

       call WoodProducts(num_soilc, filter_soilc )

       call CropHarvestPools(num_soilc, filter_soilc, dt)

       call FireArea(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            atm2lnd_vars, energyflux_vars, soilhydrology_vars, &
            cnstate_vars )

       call FireFluxes(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars)


       call t_stop_lnd(event)

   else
       call alm_fates%wrap_WoodProducts(bounds, num_soilc, filter_soilc)

       call WoodProducts(num_soilc, filter_soilc )

       call CropHarvestPools(num_soilc, filter_soilc, dt)
      
   end if

   if ( ero_ccycle ) then
       event = 'ErosionFluxes'
       call t_start_lnd(event)
       call ErosionFluxes(bounds, num_soilc, filter_soilc, soilstate_vars, sedflux_vars )
       call t_stop_lnd(event)
   end if

   !--------------------------------------------
   ! Update3
   !--------------------------------------------

   if(.not.use_fates)then
       if ( use_c13 ) then
       call CarbonIsoFlux3(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars, &
            isotope=c13, isocol_cs=c13_col_cs, isoveg_cs=c13_veg_cs, isocol_cf=c13_col_cf, isoveg_cf=c13_veg_cf)
       end if
       if ( use_c14 ) then
          call CarbonIsoFlux3(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars , &
               isotope=c14, isocol_cs=c14_col_cs, isoveg_cs=c14_veg_cs, isocol_cf=c14_col_cf, isoveg_cf=c14_veg_cf)
       end if

       call CarbonStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
            col_cs, veg_cs, col_cf, veg_cf, dt)

       if ( use_c13 ) then
          call CarbonStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c13_col_cs, c13_veg_cs, c13_col_cf, c13_veg_cf, dt)
       end if
       if ( use_c14 ) then
          call CarbonStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c14_col_cs, c14_veg_cs, c14_col_cf, c14_veg_cf, dt)
       end if


       if ( use_c14 ) then
          call C14Decay(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars )

          call C14BombSpike(num_soilp, filter_soilp, &
               cnstate_vars)
       end if

       call veg_cf_summary_for_ch4(veg_cf,bounds, num_soilp, filter_soilp)
       if( use_c13 ) then
          call col_cf_summary_for_ch4(c13_col_cf,bounds, num_soilc, filter_soilc)
          call veg_cf_summary_for_ch4(c13_veg_cf,bounds, num_soilp, filter_soilp)
       endif
       if( use_c14 ) then
          call col_cf_summary_for_ch4(c14_col_cf,bounds, num_soilc, filter_soilc)
          call veg_cf_summary_for_ch4(c14_veg_cf,bounds, num_soilp, filter_soilp)
       endif

  end if !end of if not use_fates block

  call col_cf_summary_for_ch4(col_cf,bounds, num_soilc, filter_soilc)



end subroutine EcosystemDynNoLeaching2


end  module EcosystemDynMod
