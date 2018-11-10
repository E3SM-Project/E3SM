module EcosystemDynMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Ecosystem dynamics: phenology, vegetation
  !
  ! !USES:
  use dynSubgridControlMod, only : get_do_harvest
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_sys_mod         , only : shr_sys_flush
  use clm_varctl          , only : use_c13, use_c14, use_fates, use_dynroot
  use decompMod           , only : bounds_type
  use perf_mod            , only : t_startf, t_stopf
  use spmdMod             , only : masterproc
  use clm_varctl          , only : use_century_decomp
  use CNStateType         , only : cnstate_type
  use CNCarbonFluxType    , only : carbonflux_type
  use CNCarbonStateType   , only : carbonstate_type
  use CNNitrogenFluxType  , only : nitrogenflux_type
  use CNNitrogenStateType , only : nitrogenstate_type
  use CNDVType            , only : dgvs_type
  use CanopyStateType     , only : canopystate_type
  use SoilStateType       , only : soilstate_type
  use TemperatureType     , only : temperature_type
  use WaterstateType      , only : waterstate_type
  use WaterfluxType       , only : waterflux_type
  use atm2lndType         , only : atm2lnd_type
  use SoilStateType       , only : soilstate_type
  use CanopyStateType     , only : canopystate_type
  use TemperatureType     , only : temperature_type 
  use PhotosynthesisType  , only : photosyns_type
  use CH4Mod              , only : ch4_type
  use EnergyFluxType      , only : energyflux_type
  use SoilHydrologyType   , only : soilhydrology_type
  use FrictionVelocityType, only : frictionvel_type
  use PhosphorusFluxType  , only : phosphorusflux_type
  use PhosphorusStateType , only : phosphorusstate_type
  ! bgc interface & pflotran
  use clm_varctl          , only : use_clm_interface, use_clm_bgc, use_pflotran, pf_cmode, pf_hmode
  use VerticalProfileMod   , only : decomp_vertprofiles
  use AllocationMod     , only : nu_com_nfix, nu_com_phosphatase
  use clm_varctl          , only : nu_com
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
  subroutine EcosystemDynInit(bounds)
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
    !-----------------------------------------------------------------------

    call AllocationInit (bounds)
    call PhenologyInit  (bounds)
    call FireInit       (bounds)
    
    if ( use_c14 ) then
       call C14_init_BombSpike()
    end if

  end subroutine EcosystemDynInit


  !-----------------------------------------------------------------------

  subroutine EcosystemDynLeaching(bounds, num_soilc, filter_soilc, &
       num_soilp, filter_soilp, num_pcropp, filter_pcropp, doalb, &
       cnstate_vars, carbonflux_vars, carbonstate_vars, &
       c13_carbonflux_vars, c13_carbonstate_vars, &
       c14_carbonflux_vars, c14_carbonstate_vars, dgvs_vars, &
       nitrogenflux_vars, nitrogenstate_vars, &
       waterstate_vars, waterflux_vars, frictionvel_vars, canopystate_vars,&
       phosphorusflux_vars,phosphorusstate_vars)
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
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
    type(carbonstate_type)   , intent(inout) :: carbonstate_vars
    type(carbonflux_type)    , intent(inout) :: c13_carbonflux_vars
    type(carbonstate_type)   , intent(inout) :: c13_carbonstate_vars
    type(carbonflux_type)    , intent(inout) :: c14_carbonflux_vars
    type(carbonstate_type)   , intent(inout) :: c14_carbonstate_vars
    type(dgvs_type)          , intent(in)    :: dgvs_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(waterstate_type)    , intent(in)    :: waterstate_vars
    type(waterflux_type)     , intent(in)    :: waterflux_vars
    type(frictionvel_type)   , intent(in)    :: frictionvel_vars
    type(canopystate_type)   , intent(inout) :: canopystate_vars

    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars

    !-----------------------------------------------------------------------
  
    ! only do if ed is off
    if( .not. use_fates) then
       !if(.not.(use_pflotran.and.pf_cmode)) then
             call t_startf('PhosphorusWeathering')
             call PhosphorusWeathering(num_soilc, filter_soilc, &
                  cnstate_vars,phosphorusstate_vars,phosphorusflux_vars)
             call t_stopf('PhosphorusWeathering')

             call t_startf('PhosphorusAdsportion')
             call PhosphorusAdsportion(num_soilc, filter_soilc, &
                  cnstate_vars,phosphorusstate_vars,phosphorusflux_vars)
             call t_stopf('PhosphorusAdsportion')

             call t_startf('PhosphorusDesoprtion')
             call PhosphorusDesoprtion(num_soilc, filter_soilc, &
                  cnstate_vars,phosphorusstate_vars,phosphorusflux_vars)
             call t_stopf('PhosphorusDesoprtion')

             call t_startf('PhosphorusOcclusion')
             call PhosphorusOcclusion(num_soilc, filter_soilc, &
                  cnstate_vars,phosphorusstate_vars,phosphorusflux_vars)
             call t_stopf('PhosphorusOcclusion')

             if (.not. nu_com_phosphatase) then
                call t_startf('PhosphorusBiochemMin')
                call PhosphorusBiochemMin(bounds,num_soilc, filter_soilc, &
                     cnstate_vars,phosphorusstate_vars,phosphorusflux_vars)
                call t_stopf('PhosphorusBiochemMin')
             else
                ! nu_com_phosphatase is true
                !call t_startf('PhosphorusBiochemMin')
                !call PhosphorusBiochemMin_balance(bounds,num_soilc, filter_soilc, &
                !     cnstate_vars,nitrogenstate_vars,phosphorusstate_vars,phosphorusflux_vars)
                !call t_stopf('PhosphorusBiochemMin')
             end if
       !end if
       
       !-----------------------------------------------------------------------
       ! pflotran: when both 'pf-bgc' and 'pf-h' on, no need to call CLM-CN's N leaching module
       if (.not. (pf_cmode .and. pf_hmode)) then
         call NitrogenLeaching(bounds, num_soilc, filter_soilc, &
            waterstate_vars, waterflux_vars, nitrogenstate_vars, nitrogenflux_vars)

         call PhosphorusLeaching(bounds, num_soilc, filter_soilc, &
            waterstate_vars, waterflux_vars, phosphorusstate_vars, phosphorusflux_vars)
       end if !(.not. (pf_cmode .and. pf_hmode))
       !-----------------------------------------------------------------------

       call t_startf('CNUpdate3')

       call NitrogenStateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            nitrogenflux_vars, nitrogenstate_vars)
       call t_stopf('CNUpdate3')


       call t_startf('PUpdate3')
       call PhosphorusStateUpdate3(bounds,num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars,phosphorusflux_vars, phosphorusstate_vars)
       call t_stopf('PUpdate3')

       call t_startf('CNPsum')
       call PrecisionControl(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, nitrogenstate_vars,phosphorusstate_vars)

       call carbonflux_vars%summary_cflux_for_ch4(bounds, num_soilp, filter_soilp, num_soilc, filter_soilc)
       call carbonflux_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, 'bulk')
       call carbonstate_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
       if ( use_c13 ) then
          call c13_carbonflux_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, 'c13')
          call c13_carbonstate_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
       end if
       if ( use_c14 ) then
          call c14_carbonflux_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, 'c14')
          call c14_carbonstate_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
       end if
       call nitrogenflux_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
       call nitrogenstate_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)

       call phosphorusflux_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
       call phosphorusstate_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)


       call t_stopf('CNPsum')

    end if !end of if not use_fates block

  end subroutine EcosystemDynLeaching


!-------------------------------------------------------------------------------------------------
  subroutine EcosystemDynNoLeaching1(bounds,                          &
       num_soilc, filter_soilc,                                         &
       num_soilp, filter_soilp,                                         &
       cnstate_vars, carbonflux_vars, carbonstate_vars,                 &
       c13_carbonflux_vars,                                             &
       c14_carbonflux_vars,                                             &
       nitrogenflux_vars, nitrogenstate_vars,                           &
       atm2lnd_vars, waterstate_vars, waterflux_vars,                   &
       canopystate_vars, soilstate_vars, temperature_vars, crop_vars,   &
       ch4_vars, photosyns_vars,                                        &
       phosphorusflux_vars,phosphorusstate_vars)
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
    use PhosphorusDynamicsMod           , only: PhosphorusDeposition   
    use MaintenanceRespMod             , only: MaintenanceResp
!    use SoilLittDecompMod            , only: SoilLittDecompAlloc
!    use PhenologyMod         , only: Phenology
!    use GrowthRespMod             , only: GrowthResp
!    use CarbonStateUpdate1Mod     , only: CarbonStateUpdate1,CarbonStateUpdate0
!    use NitrogenStateUpdate1Mod     , only: NitrogenStateUpdate1
!    use PhosphorusStateUpdate1Mod       , only: PhosphorusStateUpdate1
!    use GapMortalityMod      , only: GapMortality
!    use CarbonStateUpdate2Mod     , only: CarbonStateUpdate2, CarbonStateUpdate2h
!    use NitrogenStateUpdate2Mod     , only: NitrogenStateUpdate2, NitrogenStateUpdate2h
!    use PhosphorusStateUpdate2Mod       , only: PhosphorusStateUpdate2, PhosphorusStateUpdate2h
!    use FireMod              , only: FireArea, FireFluxes
!    use CarbonStateUpdate3Mod     , only: CarbonStateUpdate3
!    use CarbonIsoFluxMod          , only: CarbonIsoFlux1, CarbonIsoFlux2, CarbonIsoFlux2h, CarbonIsoFlux3
!    use C14DecayMod          , only: C14Decay, C14BombSpike
!    use WoodProductsMod      , only: WoodProducts
!    use SoilLittVertTranspMod, only: SoilLittVertTransp
    use DecompCascadeBGCMod  , only: decomp_rate_constants_bgc
    use DecompCascadeCNMod   , only: decomp_rate_constants_cn
    use CropType               , only: crop_type
!    use dynHarvestMod          , only: CNHarvest
    use clm_varpar             , only: crop_prog
    use AllocationMod        , only: Allocation1_PlantNPDemand ! Phase-1 of CNAllocation
!    use SoilLittDecompMod            , only: SoilLittDecompAlloc2
    use NitrogenDynamicsMod         , only: NitrogenLeaching
    use PhosphorusDynamicsMod           , only: PhosphorusLeaching
    use NitrogenDynamicsMod         , only: NitrogenFixation_balance
    use PhosphorusDynamicsMod           , only: PhosphorusWeathering,PhosphorusAdsportion,PhosphorusDesoprtion,PhosphorusOcclusion
    use PhosphorusDynamicsMod           , only: PhosphorusBiochemMin,PhosphorusBiochemMin_balance
  
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                  , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:)   ! filter for soil patches
!    integer                  , intent(in)    :: num_pcropp        ! number of prog. crop patches in filter
!    integer                  , intent(in)    :: filter_pcropp(:)  ! filter for prognostic crop patches
!    logical                  , intent(in)    :: doalb             ! true = surface albedo calculation time step
    type(cnstate_type)       , intent(inout) :: cnstate_vars
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
    type(carbonstate_type)   , intent(inout) :: carbonstate_vars
    type(carbonflux_type)    , intent(inout) :: c13_carbonflux_vars
!    type(carbonstate_type)   , intent(inout) :: c13_carbonstate_vars
    type(carbonflux_type)    , intent(inout) :: c14_carbonflux_vars
!    type(carbonstate_type)   , intent(inout) :: c14_carbonstate_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(atm2lnd_type)       , intent(in)    :: atm2lnd_vars
    type(waterstate_type)    , intent(in)    :: waterstate_vars
    type(waterflux_type)     , intent(in)    :: waterflux_vars
    type(canopystate_type)   , intent(in)    :: canopystate_vars
    type(soilstate_type)     , intent(inout) :: soilstate_vars
    type(temperature_type)   , intent(inout) :: temperature_vars
    type(crop_type)          , intent(in)    :: crop_vars
    type(ch4_type)           , intent(in)    :: ch4_vars
!    type(dgvs_type)          , intent(inout) :: dgvs_vars
    type(photosyns_type)     , intent(in)    :: photosyns_vars
!    type(soilhydrology_type) , intent(in)    :: soilhydrology_vars
!    type(energyflux_type)    , intent(in)    :: energyflux_vars
!
    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars

    !-----------------------------------------------------------------------

    ! Call the main CN routines

    ! only do if ed is off
    if( .not. use_fates ) then

       ! --------------------------------------------------
       ! zero the column-level C and N fluxes
       ! --------------------------------------------------

       call t_startf('CNZero')

       call carbonflux_vars%SetValues( &
            num_soilp, filter_soilp, 0._r8, &
            num_soilc, filter_soilc, 0._r8)

       if ( use_c13 ) then
          call c13_carbonflux_vars%SetValues( &
               num_soilp, filter_soilp, 0._r8, &
               num_soilc, filter_soilc, 0._r8)
       end if

       if ( use_c14 ) then
          call c14_carbonflux_vars%SetValues( &
               num_soilp, filter_soilp, 0._r8, &
               num_soilc, filter_soilc, 0._r8)
       end if
       call nitrogenflux_vars%SetValues( &
            num_soilp, filter_soilp, 0._r8, &
            num_soilc, filter_soilc, 0._r8)

       call phosphorusflux_vars%SetValues( &
            num_soilp, filter_soilp, 0._r8, &
            num_soilc, filter_soilc, 0._r8)

       call t_stopf('CNZero')

       ! --------------------------------------------------
       ! Nitrogen Deposition, Fixation and Respiration, phosphorus dynamics
       ! --------------------------------------------------

       call t_startf('CNDeposition')
       call NitrogenDeposition(bounds, &
            atm2lnd_vars, nitrogenflux_vars)
       call t_stopf('CNDeposition')

       if (.not. nu_com_nfix) then 
          call t_startf('CNFixation')
          call NitrogenFixation( num_soilc, filter_soilc, &
               waterflux_vars, carbonflux_vars, nitrogenflux_vars)
          call t_stopf('CNFixation')
       else
          ! nu_com_nfix is true
          call t_startf('CNFixation')
          call NitrogenFixation_balance( num_soilc, filter_soilc, &
               cnstate_vars, carbonflux_vars, nitrogenstate_vars, nitrogenflux_vars, &
               temperature_vars, waterstate_vars, carbonstate_vars, phosphorusstate_vars)
          call t_stopf('CNFixation')
       end if

       call t_startf('MaintenanceResp')
       if (crop_prog) then
          call NitrogenFert(bounds, num_soilc,filter_soilc, &
               nitrogenflux_vars)

          call CNSoyfix(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
               waterstate_vars, crop_vars, cnstate_vars, &
               nitrogenstate_vars, nitrogenflux_vars)
       end if
       call MaintenanceResp(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            canopystate_vars, soilstate_vars, temperature_vars, photosyns_vars, &
            carbonflux_vars, carbonstate_vars, nitrogenstate_vars)
       call t_stopf('MaintenanceResp')

       if ( nu_com .ne. 'RD') then
          ! for P competition purpose, calculate P fluxes that will potentially increase solution P pool
          ! then competitors take up solution P
          call t_startf('PhosphorusWeathering')
          call PhosphorusWeathering(num_soilc, filter_soilc, &
               cnstate_vars,phosphorusstate_vars,phosphorusflux_vars)
          call t_stopf('PhosphorusWeathering')

          if (.not. nu_com_phosphatase) then
             call t_startf('PhosphorusBiochemMin')
             call PhosphorusBiochemMin(bounds,num_soilc, filter_soilc, &
                  cnstate_vars,phosphorusstate_vars,phosphorusflux_vars)
             call t_stopf('PhosphorusBiochemMin')
          else
             ! nu_com_phosphatase is true
             call t_startf('PhosphorusBiochemMin')
             call PhosphorusBiochemMin_balance(bounds,num_soilc, filter_soilc, &
                  cnstate_vars,nitrogenstate_vars,phosphorusstate_vars,phosphorusflux_vars)
             call t_stopf('PhosphorusBiochemMin')
          end if
       end if

       ! --------------------------------------------------
       ! Phosphorus Deposition ! X.SHI
       ! --------------------------------------------------

       call t_startf('PhosphorusDeposition')
       call PhosphorusDeposition(bounds, &
            atm2lnd_vars, phosphorusflux_vars)
       call t_stopf('PhosphorusDeposition')

       !-------------------------------------------------------------------------------------------------
       ! plfotran: 'decomp_rate_constants' must be calculated before entering "clm_interface"
       if (use_century_decomp) then
          call decomp_rate_constants_bgc(bounds, num_soilc, filter_soilc, &
               canopystate_vars, soilstate_vars, temperature_vars, ch4_vars, carbonflux_vars, cnstate_vars)
       else
          call decomp_rate_constants_cn(bounds, num_soilc, filter_soilc, &
               canopystate_vars, soilstate_vars, temperature_vars, ch4_vars, carbonflux_vars, cnstate_vars)
       end if

       !-------------------------------------------------------------------------------------------------
       ! 'decomp_vertprofiles' (calc nfixation_prof) is moved from SoilLittDecompAlloc:
       ! 'nfixation_prof' is used to 'calc_nuptake_prof' & 'calc_puptake_prof', which are called in Allocation1,2,3
       call decomp_vertprofiles(bounds,                      &
           num_soilc, filter_soilc, num_soilp, filter_soilp, &
           soilstate_vars, canopystate_vars, cnstate_vars)
       !-------------------------------------------------------------------------------------------------
       ! Allocation1 is always called (w/ or w/o use_clm_interface)
       ! pflotran: call 'Allocation1' to obtain potential N demand for support initial GPP
       call t_startf('CNAllocation - phase-1')
       call Allocation1_PlantNPDemand (bounds                             , &
                num_soilc, filter_soilc, num_soilp, filter_soilp            , &
                photosyns_vars, crop_vars, canopystate_vars, cnstate_vars   , &
                carbonstate_vars, carbonflux_vars, c13_carbonflux_vars      , &
                c14_carbonflux_vars, nitrogenstate_vars, nitrogenflux_vars  , &
                phosphorusstate_vars, phosphorusflux_vars)

       call t_stopf('CNAllocation - phase-1')

    end if !end of if not use_fates block

  end subroutine EcosystemDynNoLeaching1

!-------------------------------------------------------------------------------------------------
  subroutine EcosystemDynNoLeaching2(bounds,                                  &
       num_soilc, filter_soilc,                                                 &
       num_soilp, filter_soilp, num_pcropp, filter_pcropp, doalb,               &
       cnstate_vars, carbonflux_vars, carbonstate_vars,                         &
       c13_carbonflux_vars, c13_carbonstate_vars,                               &
       c14_carbonflux_vars, c14_carbonstate_vars,                               &
       nitrogenflux_vars, nitrogenstate_vars,                                   &
       atm2lnd_vars, waterstate_vars, waterflux_vars,                           &
       canopystate_vars, soilstate_vars, temperature_vars, crop_vars, ch4_vars, &
       dgvs_vars, photosyns_vars, soilhydrology_vars, energyflux_vars,          &
       phosphorusflux_vars,phosphorusstate_vars)
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
!    use NitrogenDynamicsMod         , only: NitrogenDeposition,NitrogenFixation, NitrogenFert, CNSoyfix
!    use MaintenanceRespMod             , only: MaintenanceResp
!    use SoilLittDecompMod            , only: SoilLittDecompAlloc
    use PhenologyMod         , only: Phenology
    use GrowthRespMod             , only: GrowthResp
    use CarbonStateUpdate1Mod     , only: CarbonStateUpdate1,CarbonStateUpdate0
    use NitrogenStateUpdate1Mod     , only: NitrogenStateUpdate1
    use PhosphorusStateUpdate1Mod       , only: PhosphorusStateUpdate1
    use GapMortalityMod        , only: GapMortality
    use CarbonStateUpdate2Mod     , only: CarbonStateUpdate2, CarbonStateUpdate2h
    use NitrogenStateUpdate2Mod     , only: NitrogenStateUpdate2, NitrogenStateUpdate2h
    use PhosphorusStateUpdate2Mod       , only: PhosphorusStateUpdate2, PhosphorusStateUpdate2h
    use FireMod              , only: FireArea, FireFluxes
    use CarbonStateUpdate3Mod     , only: CarbonStateUpdate3
    use CarbonIsoFluxMod          , only: CarbonIsoFlux1, CarbonIsoFlux2, CarbonIsoFlux2h, CarbonIsoFlux3
    use C14DecayMod          , only: C14Decay, C14BombSpike
    use WoodProductsMod      , only: WoodProducts
    use CropHarvestPoolsMod    , only: CropHarvestPools
    use SoilLittVertTranspMod, only: SoilLittVertTransp
!    use DecompCascadeBGCMod  , only: decomp_rate_constants_bgc
!    use DecompCascadeCNMod   , only: decomp_rate_constants_cn
    use CropType               , only: crop_type
    use dynHarvestMod          , only: CNHarvest
    use RootDynamicsMod           , only: RootDynamics
!    use clm_varpar             , only: crop_prog

!    use AllocationMod        , only: cnallocation
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
    logical                  , intent(in)    :: doalb             ! true = surface albedo calculation time step
    type(cnstate_type)       , intent(inout) :: cnstate_vars
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
    type(carbonstate_type)   , intent(inout) :: carbonstate_vars
    type(carbonflux_type)    , intent(inout) :: c13_carbonflux_vars
    type(carbonstate_type)   , intent(inout) :: c13_carbonstate_vars
    type(carbonflux_type)    , intent(inout) :: c14_carbonflux_vars
    type(carbonstate_type)   , intent(inout) :: c14_carbonstate_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(atm2lnd_type)       , intent(in)    :: atm2lnd_vars
    type(waterstate_type)    , intent(in)    :: waterstate_vars
    type(waterflux_type)     , intent(in)    :: waterflux_vars
    type(canopystate_type)   , intent(in)    :: canopystate_vars
    type(soilstate_type)     , intent(inout) :: soilstate_vars
    type(temperature_type)   , intent(inout) :: temperature_vars
    type(crop_type)          , intent(inout) :: crop_vars
    type(ch4_type)           , intent(in)    :: ch4_vars
    type(dgvs_type)          , intent(inout) :: dgvs_vars
    type(photosyns_type)     , intent(in)    :: photosyns_vars
    type(soilhydrology_type) , intent(in)    :: soilhydrology_vars
    type(energyflux_type)    , intent(in)    :: energyflux_vars
!
    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars

    !-----------------------------------------------------------------------

    ! Call the main CN routines
    ! only do if ed is off
    if( .not. use_fates ) then

       call t_startf('SoilLittDecompAlloc')
       !----------------------------------------------------------------
       if(.not.use_clm_interface) then
            ! directly run clm-bgc
            ! if (use_clm_interface & use_clm_bgc), then CNDecomAlloc is called in clm_driver
            call SoilLittDecompAlloc (bounds, num_soilc, filter_soilc,    &
                       num_soilp, filter_soilp,                     &
                       canopystate_vars, soilstate_vars,            &
                       temperature_vars, waterstate_vars,           &
                       cnstate_vars, ch4_vars,                      &
                       carbonstate_vars, carbonflux_vars,           &
                       nitrogenstate_vars, nitrogenflux_vars,       &
                       phosphorusstate_vars,phosphorusflux_vars)
       end if !if(.not.use_clm_interface)
       !----------------------------------------------------------------
       ! SoilLittDecompAlloc2 is called by both clm-bgc & pflotran
       ! pflotran: call 'SoilLittDecompAlloc2' to calculate some diagnostic variables and 'fpg' for plant N uptake
       ! pflotran & clm-bgc : 'Allocation3_AG' and vertically integrate net and gross mineralization fluxes
       call SoilLittDecompAlloc2 (bounds, num_soilc, filter_soilc, num_soilp, filter_soilp,           &
                photosyns_vars, canopystate_vars, soilstate_vars, temperature_vars,             &
                waterstate_vars, cnstate_vars, ch4_vars,                                        &
                carbonstate_vars, carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars,    &
                nitrogenstate_vars, nitrogenflux_vars, crop_vars, atm2lnd_vars,                 &
                phosphorusstate_vars,phosphorusflux_vars)

       !----------------------------------------------------------------
       call t_stopf('SoilLittDecompAlloc')
       !----------------------------------------------------------------

       !--------------------------------------------
       ! Phenology
       !--------------------------------------------

       ! Phenology needs to be called after SoilLittDecompAlloc, because it
       ! depends on current time-step fluxes to new growth on the last
       ! litterfall timestep in deciduous systems

       call t_startf('Phenology')
       call Phenology(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            num_pcropp, filter_pcropp, doalb, atm2lnd_vars, &
            waterstate_vars, temperature_vars, crop_vars, canopystate_vars, soilstate_vars, &
            dgvs_vars, cnstate_vars, carbonstate_vars, carbonflux_vars, &
            nitrogenstate_vars, nitrogenflux_vars,&
            phosphorusstate_vars,phosphorusflux_vars)
       call t_stopf('Phenology')

       !--------------------------------------------
       ! Growth respiration
       !--------------------------------------------

       call t_startf('GrowthResp')
       call GrowthResp(num_soilp, filter_soilp, &
            carbonflux_vars)
       call t_stopf('GrowthResp')
       call carbonflux_vars%summary_rr(bounds, num_soilp, filter_soilp, num_soilc, filter_soilc)

       if(use_c13) then
         call c13_carbonflux_vars%summary_rr(bounds, num_soilp, filter_soilp, num_soilc, filter_soilc)
       endif

       if(use_c14) then
         call c14_carbonflux_vars%summary_rr(bounds, num_soilp, filter_soilp, num_soilc, filter_soilc)
       endif
       !--------------------------------------------
       ! Dynamic Roots
       !--------------------------------------------

       if( use_dynroot ) then
          call t_startf('RootDynamics')

          call RootDynamics(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
               canopystate_vars, carbonstate_vars, nitrogenstate_vars, carbonflux_vars,  &
               cnstate_vars, crop_vars, energyflux_vars, soilstate_vars)
          call t_stopf('RootDynamics')
       end if

       !--------------------------------------------
       ! CNUpdate0
       !--------------------------------------------

       call t_startf('CNUpdate0')
       call CarbonStateUpdate0(&
            num_soilp, filter_soilp, &
            carbonflux_vars, carbonstate_vars)

       if ( use_c13 ) then
          call CarbonStateUpdate0(&
               num_soilp, filter_soilp, &
               c13_carbonflux_vars, c13_carbonstate_vars)
       end if

       if ( use_c14 ) then
          call CarbonStateUpdate0(&
               num_soilp, filter_soilp, &
               c14_carbonflux_vars, c14_carbonstate_vars)
       end if
       call t_stopf('CNUpdate0')

       !--------------------------------------------
       ! Update1
       !--------------------------------------------

       call t_startf('CNUpdate1')

       if ( use_c13 ) then
          call CarbonIsoFlux1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c13_carbonflux_vars, isotopestate_vars=c13_carbonstate_vars, &
               isotope='c13')
       end if

       if ( use_c14 ) then
          call CarbonIsoFlux1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c14_carbonflux_vars, isotopestate_vars=c14_carbonstate_vars, &
               isotope='c14')
       end if

       call CarbonStateUpdate1(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            crop_vars, carbonflux_vars, carbonstate_vars)

       if ( use_c13 ) then
          call CarbonStateUpdate1(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
               crop_vars, c13_carbonflux_vars, c13_carbonstate_vars)
       end if
       if ( use_c14 ) then
          call CarbonStateUpdate1(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
               crop_vars, c14_carbonflux_vars, c14_carbonstate_vars)
       end if

       call NitrogenStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars, nitrogenflux_vars, nitrogenstate_vars)

       call PhosphorusStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars, phosphorusflux_vars, phosphorusstate_vars)

       call t_stopf('CNUpdate1')

       call t_startf('SoilLittVertTransp')
       call SoilLittVertTransp(bounds, &
            num_soilc, filter_soilc, &
            canopystate_vars, cnstate_vars,                               &
            carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, &
            carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars,    &
            nitrogenstate_vars, nitrogenflux_vars,&
            phosphorusstate_vars,phosphorusflux_vars)
       call t_stopf('SoilLittVertTransp')

       call t_startf('CNGapMortality')
       call GapMortality( num_soilc, filter_soilc, num_soilp, filter_soilp, &
            dgvs_vars, cnstate_vars, &
            carbonstate_vars, nitrogenstate_vars, carbonflux_vars, nitrogenflux_vars,&
            phosphorusstate_vars,phosphorusflux_vars )
       call t_stopf('CNGapMortality')

       !--------------------------------------------
       ! Update2
       !--------------------------------------------

       call t_startf('CNUpdate2')

       if ( use_c13 ) then
          call CarbonIsoFlux2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c13_carbonflux_vars, isotopestate_vars=c13_carbonstate_vars, &
               isotope='c13')
       end if

       if ( use_c14 ) then
          call CarbonIsoFlux2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c14_carbonflux_vars, isotopestate_vars=c14_carbonstate_vars, &
               isotope='c14')
       end if

       call CarbonStateUpdate2( num_soilc, filter_soilc, num_soilp, filter_soilp, &
            carbonflux_vars, carbonstate_vars)

       if ( use_c13 ) then
          call CarbonStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c13_carbonflux_vars, c13_carbonstate_vars)
       end if
       if ( use_c14 ) then
          call CarbonStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c14_carbonflux_vars, c14_carbonstate_vars)
       end if
       call NitrogenStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            nitrogenflux_vars, nitrogenstate_vars)

       call PhosphorusStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            phosphorusflux_vars, phosphorusstate_vars)

       if (get_do_harvest()) then
          call CNHarvest(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonstate_vars, nitrogenstate_vars, carbonflux_vars, nitrogenflux_vars,&
               phosphorusstate_vars, phosphorusflux_vars)
       end if

       if ( use_c13 ) then
          call CarbonIsoFlux2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c13_carbonflux_vars, isotopestate_vars=c13_carbonstate_vars, &
               isotope='c13')
       end if
       if ( use_c14 ) then
          call CarbonIsoFlux2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c14_carbonflux_vars, isotopestate_vars=c14_carbonstate_vars, &
               isotope='c14')
       end if

       call CarbonStateUpdate2h( num_soilc, filter_soilc,  num_soilp, filter_soilp, &
            carbonflux_vars, carbonstate_vars)
       if ( use_c13 ) then
          call CarbonStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c13_carbonflux_vars, c13_carbonstate_vars)
       end if
       if ( use_c14 ) then
          call CarbonStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c14_carbonflux_vars, c14_carbonstate_vars)
       end if

       call NitrogenStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            nitrogenflux_vars, nitrogenstate_vars)

       call PhosphorusStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            phosphorusflux_vars, phosphorusstate_vars)

       call WoodProducts(num_soilc, filter_soilc, &
            carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, nitrogenstate_vars, &
            carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars, nitrogenflux_vars,&
            phosphorusstate_vars,phosphorusflux_vars)

       call CropHarvestPools(num_soilc, filter_soilc, &
            carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, nitrogenstate_vars, &
            phosphorusstate_vars, carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars, &
            nitrogenflux_vars, phosphorusflux_vars)

       call FireArea(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            atm2lnd_vars, temperature_vars, energyflux_vars, soilhydrology_vars, waterstate_vars, &
            cnstate_vars, carbonstate_vars)

       call FireFluxes(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            dgvs_vars, cnstate_vars, carbonstate_vars, nitrogenstate_vars, &
            carbonflux_vars,nitrogenflux_vars,phosphorusstate_vars,phosphorusflux_vars)

       call t_stopf('CNUpdate2')

       !--------------------------------------------
       ! Update3
       !--------------------------------------------

       if ( use_c13 ) then
          call CarbonIsoFlux3(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c13_carbonflux_vars, isotopestate_vars=c13_carbonstate_vars, &
               isotope='c13')
       end if
       if ( use_c14 ) then
          call CarbonIsoFlux3(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c14_carbonflux_vars, isotopestate_vars=c14_carbonstate_vars, &
               isotope='c14')
       end if

       call CarbonStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
            carbonflux_vars, carbonstate_vars)

       if ( use_c13 ) then
          call CarbonStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c13_carbonflux_vars, c13_carbonstate_vars)
       end if
       if ( use_c14 ) then
          call CarbonStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c14_carbonflux_vars, c14_carbonstate_vars)
       end if


       if ( use_c14 ) then
          call C14Decay(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, c14_carbonstate_vars)

          call C14BombSpike(num_soilp, filter_soilp, &
               cnstate_vars)
       end if

       call carbonflux_vars%summary_cflux_for_ch4(bounds, num_soilp, filter_soilp, num_soilc, filter_soilc)
       if( use_c13 ) then
          call c13_carbonflux_vars%summary_cflux_for_ch4(bounds, num_soilp, filter_soilp, num_soilc, filter_soilc)
       endif
       if( use_c14 ) then
          call c14_carbonflux_vars%summary_cflux_for_ch4(bounds, num_soilp, filter_soilp, num_soilc, filter_soilc)
       endif    
    end if !end of if not use_fates block

  end subroutine EcosystemDynNoLeaching2

  
end  module EcosystemDynMod
