module EcosystemDynMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Ecosystem dynamics: phenology, vegetation
  !
  ! !USES:
  use dynSubgridControlMod, only : get_do_harvest
  use shr_kind_mod        , only : r8 => shr_kind_r8
   !#py use shr_sys_mod         , only : shr_sys_flush
  use clm_varctl          , only : use_c13, use_c14, use_fates, use_dynroot
  use decompMod           , only : bounds_type
   !#py use perf_mod            , only : t_startf, t_stopf
   !#py use spmdMod             , only : masterproc
  use clm_varctl          , only : use_century_decomp
  use clm_varctl          , only : use_erosion
  use CNStateType         , only : cnstate_type
  use CanopyStateType     , only : canopystate_type
  use SoilStateType       , only : soilstate_type

  use atm2lndType         , only : atm2lnd_type
  use PhotosynthesisType  , only : photosyns_type
  use CH4Mod              , only : ch4_type
  use EnergyFluxType      , only : energyflux_type
  use SoilHydrologyType   , only : soilhydrology_type
  use FrictionVelocityType, only : frictionvel_type
  !NEW FROM MASTER
  use SedFluxType         , only : sedflux_type
  use ColumnDataType      , only : col_cs, c13_col_cs, c14_col_cs
  use ColumnDataType      , only : col_cf, c13_col_cf, c14_col_cf
  use ColumnDataType      , only : col_ns, col_nf
  use ColumnDataType      , only : col_ps, col_pf
  use VegetationDataType  , only : veg_cs, c13_veg_cs, c14_veg_cs
  use VegetationDataType  , only : veg_cf, c13_veg_cf, c14_veg_cf
  use VegetationDataType  , only : veg_ns, veg_nf
  use VegetationDataType  , only : veg_ps, veg_pf

  ! bgc interface & pflotran
  use clm_varctl          , only : use_clm_interface, use_clm_bgc, use_pflotran, pf_cmode, pf_hmode
  use VerticalProfileMod   , only : decomp_vertprofiles
  use AllocationMod     , only : nu_com_nfix, nu_com_phosphatase
  use clm_varctl          , only : nu_com, use_pheno_flux_limiter
  use PhenologyFLuxLimitMod , only : phenology_flux_limiter, InitPhenoFluxLimiter
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

    if(use_pheno_flux_limiter)then
      call InitPhenoFluxLimiter()
    endif
  end subroutine EcosystemDynInit


  !-----------------------------------------------------------------------

   subroutine EcosystemDynLeaching(bounds, num_soilc, filter_soilc, &
        num_soilp, filter_soilp, num_pcropp, filter_pcropp, doalb, &
        cnstate_vars,frictionvel_vars, canopystate_vars, dt)
     !
     ! !DESCRIPTION:
     ! The core CN code is executed here. Calculates fluxes for maintenance
     ! respiration, decomposition, allocation, phenology, and growth respiration.
     ! These routines happen on the radiation time step so that canopy structure
     ! stays synchronized with albedo calculations.
     !
     ! !USES:
      !$acc routine seq
     use PhosphorusDynamicsMod   , only: PhosphorusWeathering,PhosphorusAdsportion,PhosphorusDesoprtion,PhosphorusOcclusion
     use PhosphorusDynamicsMod   , only: PhosphorusBiochemMin,PhosphorusLeaching
     use NitrogenDynamicsMod     , only: NitrogenLeaching
     use NitrogenStateUpdate3Mod , only: NitrogenStateUpdate3
     use PhosphorusStateUpdate3Mod     , only: PhosphorusStateUpdate3
     use PrecisionControlMod     , only: PrecisionControl
     use PhosphorusDynamicsMod   , only: PhosphorusBiochemMin_balance
     use Method_procs_acc        , only : vegnf_summary_acc,vegcf_summary_acc
     use Method_procs_acc , only : vegcf_summary_for_CH4_acc, colcf_summary_for_ch4_acc
     use Method_procs_acc , only : colcf_summary_acc, colcs_summary_acc, vegcs_summary_acc
     use Method_procs_acc , only : colnf_summary_acc, vegns_summary_acc,colns_summary_acc
     use Method_procs_acc , only : vegpf_summary_acc,colpf_summary_acc, vegps_summary_acc
     use Method_procs_acc , only : colps_summary_acc
     !use Method_procs_acc , only : bulk, c13, c14
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
     real(r8)                 , intent(in)    :: dt

     integer :: bulk, c13, c14
     c13 = 0; c14 = 1; bulk = 2;
     !-----------------------------------------------------------------------

     ! only do if ed is off
     if( .not. use_fates) then
        !if(.not.(use_pflotran.and.pf_cmode)) then
              !#py call t_startf('PhosphorusWeathering')
              call PhosphorusWeathering(num_soilc, filter_soilc, &
                   cnstate_vars, dt)
              !#py call t_stopf('PhosphorusWeathering')

              !#py call t_startf('PhosphorusAdsportion')
              call PhosphorusAdsportion(num_soilc, filter_soilc, &
                   cnstate_vars, dt)
              !#py call t_stopf('PhosphorusAdsportion')

              !#py call t_startf('PhosphorusDesoprtion')
              call PhosphorusDesoprtion(num_soilc, filter_soilc, &
                   cnstate_vars, dt)
              !#py call t_stopf('PhosphorusDesoprtion')

              !#py call t_startf('PhosphorusOcclusion')
              call PhosphorusOcclusion(num_soilc, filter_soilc, &
                   cnstate_vars, dt)
              !#py call t_stopf('PhosphorusOcclusion')

              if (.not. nu_com_phosphatase) then
                 !#py call t_startf('PhosphorusBiochemMin')
                 call PhosphorusBiochemMin(bounds,num_soilc, filter_soilc, &
                      cnstate_vars, dt)
                 !#py call t_stopf('PhosphorusBiochemMin')
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
          call NitrogenLeaching(bounds, num_soilc, filter_soilc, dt)

          call PhosphorusLeaching(bounds, num_soilc, filter_soilc, dt)
        end if !(.not. (pf_cmode .and. pf_hmode))
        !-----------------------------------------------------------------------

        !#py call t_startf('CNUpdate3')

        call NitrogenStateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp,dt)
        !#py call t_stopf('CNUpdate3')


        !#py call t_startf('PUpdate3')
        call PhosphorusStateUpdate3(bounds,num_soilc, filter_soilc, num_soilp, filter_soilp, &
             cnstate_vars, dt)
        !#py call t_stopf('PUpdate3')

        !#py call t_startf('CNPsum')
        call PrecisionControl(num_soilc, filter_soilc, num_soilp, filter_soilp )

        call colcf_summary_for_ch4_acc(col_cf, bounds, num_soilc, filter_soilc)
        call vegcf_Summary_for_CH4_acc(veg_cf,bounds, num_soilp, filter_soilp)

        call vegcf_summary_acc(veg_cf,bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, bulk, col_cf)
        call colcf_summary_acc(col_cf,bounds, num_soilc, filter_soilc, bulk, dt)
        if ( use_c13 ) then
           call vegcf_summary_acc(c13_veg_cf,bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, c13, c13_col_cf)
           call colcf_summary_acc(c13_col_cf,bounds, num_soilc, filter_soilc, c13, dt)
        end if
        if ( use_c14 ) then
           call vegcf_summary_acc(c14_veg_cf,bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, c14, c14_col_cf)
           call colcf_summary_acc(c14_col_cf,bounds, num_soilc, filter_soilc, c14, dt)
        end if

        call vegcs_summary_acc(veg_cs,bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_cs)
        call colcs_summary_acc(col_cs,bounds, num_soilc, filter_soilc)
        if ( use_c13 ) then
           call vegcs_summary_acc(c13_veg_cs,bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, c13_col_cs)
           call colcs_summary_acc(c13_col_cs,bounds, num_soilc, filter_soilc)
        end if
        if ( use_c14 ) then
           call vegcs_summary_acc(c14_veg_cs,bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, c14_col_cs)
           call colcs_summary_acc(c14_col_cs,bounds, num_soilc, filter_soilc)

        end if

        call vegnf_summary_acc(veg_nf,bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_nf)
        call colnf_summary_acc(col_nf,bounds, num_soilc, filter_soilc, dt)

        call vegns_summary_acc(veg_ns,bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_ns)
        call colns_summary_acc(col_ns,bounds, num_soilc, filter_soilc)

        call vegpf_summary_acc(veg_pf,bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_pf)
        call colpf_summary_acc(col_pf,bounds, num_soilc, filter_soilc, dt )

        call vegps_summary_acc(veg_ps,bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_ps)
        call colps_summary_acc(col_ps,bounds, num_soilc, filter_soilc)

        !#py call t_stopf('CNPsum')

     end if !end of if not use_fates block

   end subroutine EcosystemDynLeaching


!-------------------------------------------------------------------------------------------------
  subroutine EcosystemDynNoLeaching1(bounds,     &
       num_soilc, filter_soilc,                    &
       num_soilp, filter_soilp,                    &
       cnstate_vars,   &
       atm2lnd_vars, canopystate_vars, soilstate_vars, crop_vars,   &
       ch4_vars, photosyns_vars, dt, dayspyr,year, mon, day, sec)
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
      !$acc routine seq
    use NitrogenDynamicsMod         , only: NitrogenDeposition,NitrogenFixation, NitrogenFert, CNSoyfix
    use PhosphorusDynamicsMod           , only: PhosphorusDeposition
    use MaintenanceRespMod             , only: MaintenanceResp
    use DecompCascadeBGCMod  , only: decomp_rate_constants_bgc
    use DecompCascadeCNMod   , only: decomp_rate_constants_cn
    use CropType               , only: crop_type
    use clm_varpar             , only: crop_prog
    use AllocationMod        , only: Allocation1_PlantNPDemand ! Phase-1 of CNAllocation
    use NitrogenDynamicsMod         , only: NitrogenLeaching
    use PhosphorusDynamicsMod           , only: PhosphorusLeaching
    use NitrogenDynamicsMod         , only: NitrogenFixation_balance
    use PhosphorusDynamicsMod           , only: PhosphorusWeathering,PhosphorusAdsportion,PhosphorusDesoprtion,PhosphorusOcclusion
    use PhosphorusDynamicsMod           , only: PhosphorusBiochemMin,PhosphorusBiochemMin_balance
    use Method_procs_acc      , only : colcf_setvalues_acc, vegcf_setvalues_acc
    use Method_procs_acc      , only : colnf_setvalues_acc, vegnf_setvalues_acc
    use Method_procs_acc      , only : colpf_setvalues_acc, vegpf_setvalues_acc



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
    !type(temperature_type)   , intent(inout) :: temperature_vars
    type(crop_type)          , intent(in)    :: crop_vars
    type(ch4_type)           , intent(in)    :: ch4_vars
    type(photosyns_type)     , intent(in)    :: photosyns_vars
    real(r8), intent(in) :: dt, dayspyr
    integer, intent(in) :: year, mon, day, sec
    !-----------------------------------------------------------------------

    ! Call the main CN routines

    ! only do if ed is off
    if( .not. use_fates ) then

       ! --------------------------------------------------
       ! zero the C and N fluxes
       ! --------------------------------------------------

        !#py call t_startf('CNZero')

       call colcf_setvalues_acc(col_cf,num_soilc, filter_soilc, 0._r8)
       call vegcf_setvalues_acc(veg_cf,num_soilp, filter_soilp, 0._r8)
      if ( use_c13 ) then
         call colcf_setvalues_acc(c13_col_cf, num_soilc, filter_soilc, 0._r8)
         call vegcf_setvalues_acc(c13_veg_cf, num_soilp, filter_soilp, 0._r8)
      end if

      if ( use_c14 ) then
         call colcf_setvalues_acc(c14_col_cf,num_soilc, filter_soilc, 0._r8)
         call vegcf_setvalues_acc(c14_veg_cf,num_soilp, filter_soilp, 0._r8)
      end if

       call vegnf_SetValues_acc(veg_nf,num_soilp, filter_soilp, 0._r8)
       call colnf_SetValues_acc(col_nf,num_soilc, filter_soilc, 0._r8)

       call vegpf_SetValues_acc(veg_pf,num_soilp, filter_soilp, 0._r8)
       call colpf_SetValues_acc(col_pf,num_soilc, filter_soilc, 0._r8)

        !#py call t_stopf('CNZero')

       ! --------------------------------------------------
       ! Nitrogen Deposition, Fixation and Respiration, phosphorus dynamics
       ! --------------------------------------------------

        !#py call t_startf('CNDeposition')
       call NitrogenDeposition(bounds, &
            atm2lnd_vars, dt)
        !#py call t_stopf('CNDeposition')

       if (.not. nu_com_nfix) then
           !#py call t_startf('CNFixation')
          call NitrogenFixation( num_soilc, filter_soilc, dayspyr)
           !#py call t_stopf('CNFixation')
       else
          ! nu_com_nfix is true
           !#py call t_startf('CNFixation')
          call NitrogenFixation_balance( num_soilc, filter_soilc,cnstate_vars)
           !#py call t_stopf('CNFixation')
       end if

        !#py call t_startf('MaintenanceResp')
       if (crop_prog) then
          call NitrogenFert(bounds, num_soilc,filter_soilc)

          call CNSoyfix(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
               crop_vars, cnstate_vars)
       end if
       call MaintenanceResp(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            canopystate_vars, soilstate_vars, photosyns_vars)
        !#py call t_stopf('MaintenanceResp')

       if ( nu_com .ne. 'RD') then
          ! for P competition purpose, calculate P fluxes that will potentially increase solution P pool
          ! then competitors take up solution P
           !#py call t_startf('PhosphorusWeathering')
          call PhosphorusWeathering(num_soilc, filter_soilc, &
               cnstate_vars, dt)
           !#py call t_stopf('PhosphorusWeathering')

          if (.not. nu_com_phosphatase) then
              !#py call t_startf('PhosphorusBiochemMin')
             call PhosphorusBiochemMin(bounds,num_soilc, filter_soilc, &
                  cnstate_vars, dt)
              !#py call t_stopf('PhosphorusBiochemMin')
          else
             ! nu_com_phosphatase is true
              !#py call t_startf('PhosphorusBiochemMin')
             call PhosphorusBiochemMin_balance(bounds,num_soilc, filter_soilc, &
                  cnstate_vars, dt)
              !#py call t_stopf('PhosphorusBiochemMin')
          end if
       end if

       ! --------------------------------------------------
       ! Phosphorus Deposition ! X.SHI
       ! --------------------------------------------------

        !#py call t_startf('PhosphorusDeposition')
       call PhosphorusDeposition(bounds, &
            atm2lnd_vars)
        !#py call t_stopf('PhosphorusDeposition')

       !-------------------------------------------------------------------------------------------------
       ! plfotran: 'decomp_rate_constants' must be calculated before entering "clm_interface"
       if (use_century_decomp) then
          call decomp_rate_constants_bgc(bounds, num_soilc, filter_soilc, &
               canopystate_vars, soilstate_vars, ch4_vars, cnstate_vars,&
               dt, dayspyr,year, mon, day, sec)
       else
          call decomp_rate_constants_cn(bounds, num_soilc, filter_soilc, &
               canopystate_vars, soilstate_vars, ch4_vars, cnstate_vars, &
               dt ,year, mon, day ,sec)
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
        !#py call t_startf('CNAllocation - phase-1')
       call Allocation1_PlantNPDemand (bounds                             , &
                num_soilc, filter_soilc, num_soilp, filter_soilp            , &
                photosyns_vars, crop_vars, canopystate_vars, cnstate_vars   , &
                dt, year)

        !#py call t_stopf('CNAllocation - phase-1')

    end if !end of if not use_fates block

  end subroutine EcosystemDynNoLeaching1

!-------------------------------------------------------------------------------------------------
   subroutine EcosystemDynNoLeaching2(bounds,                         &
        num_soilc, filter_soilc,                                      &
        num_soilp, filter_soilp, num_pcropp, filter_pcropp, doalb,    &
        cnstate_vars, atm2lnd_vars,  &
        canopystate_vars, soilstate_vars, crop_vars, ch4_vars, &
        photosyns_vars, soilhydrology_vars, energyflux_vars,  sedflux_vars &
        year, mon, day, sec, tod, offset, dayspyr, dt, nstep )
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
     !#py use abortutils       , only : endrun
      !$acc routine seq
     use PhenologyMod     , only: Phenology, CNLitterToColumn
     use GrowthRespMod    , only: GrowthResp
     use CarbonStateUpdate1Mod     , only: CarbonStateUpdate1,CarbonStateUpdate0
     use NitrogenStateUpdate1Mod   , only: NitrogenStateUpdate1
     use PhosphorusStateUpdate1Mod , only: PhosphorusStateUpdate1
     use GapMortalityMod        , only: GapMortality
     use CarbonStateUpdate2Mod  , only: CarbonStateUpdate2, CarbonStateUpdate2h
     use NitrogenStateUpdate2Mod     , only: NitrogenStateUpdate2, NitrogenStateUpdate2h
     use PhosphorusStateUpdate2Mod   , only: PhosphorusStateUpdate2, PhosphorusStateUpdate2h
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
     use RootDynamicsMod           , only: RootDynamics
     use SoilLittDecompMod            , only: SoilLittDecompAlloc
     use SoilLittDecompMod            , only: SoilLittDecompAlloc2 !after SoilLittDecompAlloc
     use pftvarcon        , only : grperc, grpnow, npcropmin
     use VegetationPropertiesType   , only : veg_vp
     use CNCarbonFluxType , only : carbonflux_type
     use VegetationType        , only : veg_pp
     use VegetationDataType    , only : veg_cf
     use Method_procs_acc      , only : vegcf_summary_rr_acc
     use Method_procs_acc      , only : colcf_Summary_for_CH4_acc
     use Method_procs_acc      , only : vegcf_summary_for_CH4_acc

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
     type(atm2lnd_type)       , intent(in)    :: atm2lnd_vars
     type(canopystate_type)   , intent(in)    :: canopystate_vars
     type(soilstate_type)     , intent(inout) :: soilstate_vars
     type(crop_type)          , intent(inout) :: crop_vars
     type(ch4_type)           , intent(in)    :: ch4_vars
     type(photosyns_type)     , intent(in)    :: photosyns_vars
     type(soilhydrology_type) , intent(in)    :: soilhydrology_vars
     type(energyflux_type)    , intent(in)    :: energyflux_vars
     type(sedflux_type)       , intent(in)    :: sedflux_vars
     integer, intent(in)  :: year, mon, day,sec, tod, offset
     real(r8), intent(in) :: dayspyr ! days per year
     real(r8), intent(in) :: dt
     integer , intent(in) :: nstep

     integer :: c13, c14
      c13 = 0
      c14 = 1
     ! Call the main CN routines
     ! only do if ed is off
     if( .not. use_fates ) then

        !#py call t_startf('SoilLittDecompAlloc')
        !----------------------------------------------------------------
        if(.not.use_clm_interface) then
             ! directly run clm-bgc
             ! if (use_clm_interface & use_clm_bgc), then CNDecomAlloc is called in clm_driver
             call SoilLittDecompAlloc (bounds, num_soilc, filter_soilc,    &
                        num_soilp, filter_soilp,                     &
                        canopystate_vars, soilstate_vars,            &
                        cnstate_vars, ch4_vars,dt)
        end if !if(.not.use_clm_interface)
        !----------------------------------------------------------------
        ! SoilLittDecompAlloc2 is called by both clm-bgc & pflotran
        ! pflotran: call 'SoilLittDecompAlloc2' to calculate some diagnostic variables and 'fpg' for plant N uptake
        ! pflotran & clm-bgc : 'Allocation3_AG' and vertically integrate net and gross mineralization fluxes
        call SoilLittDecompAlloc2 (bounds, num_soilc, filter_soilc, num_soilp, filter_soilp,   &
                 photosyns_vars, canopystate_vars, soilstate_vars,     &
                 cnstate_vars, ch4_vars, crop_vars, atm2lnd_vars,dt)

         !----------------------------------------------------------------
        !#py call t_stopf('SoilLittDecompAlloc')
        !----------------------------------------------------------------

        !--------------------------------------------
        ! Phenology
        !--------------------------------------------

        ! Phenology needs to be called after SoilLittDecompAlloc, because it
        ! depends on current time-step fluxes to new growth on the last
        ! litterfall timestep in deciduous systems

        !#py call t_startf('Phenology')
        call Phenology(num_soilc, filter_soilc, num_soilp, filter_soilp, &
             num_pcropp, filter_pcropp, doalb, atm2lnd_vars, &
             crop_vars, canopystate_vars, soilstate_vars, &
             cnstate_vars)
        !#py call t_stopf('Phenology')

        !--------------------------------------------
        ! Growth respiration
        !--------------------------------------------

        !#py call t_startf('GrowthResp')

        call GrowthResp(num_soilp, filter_soilp)


        !#py call t_stopf('GrowthResp')

        call vegcf_summary_rr_acc(veg_cf,bounds, num_soilp, filter_soilp,&
          num_soilc, filter_soilc, col_cf)
        if(use_c13) then
          call vegcf_summary_rr_acc(c13_veg_cf,bounds, num_soilp, &
            filter_soilp, num_soilc, filter_soilc, c13_col_cf)
        endif
        if(use_c14) then
          call vegcf_summary_rr_acc(c14_veg_cf, bounds, num_soilp, filter_soilp , &
            num_soilc, filter_soilc, c14_col_cf)
        endif

        !--------------------------------------------
        ! Dynamic Roots
        !--------------------------------------------

        if( use_dynroot ) then
           !#py call t_startf('RootDynamics')

           call RootDynamics(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
                canopystate_vars,  &
                cnstate_vars, crop_vars, energyflux_vars, soilstate_vars)
           !#py call t_stopf('RootDynamics')
        end if

        !--------------------------------------------
        ! CNUpdate0
        !--------------------------------------------

        !#py call t_startf('CNUpdate0')
        call CarbonStateUpdate0(num_soilp, filter_soilp, veg_cs, veg_cf, dt)
        if ( use_c13 ) then
           call CarbonStateUpdate0(num_soilp, filter_soilp, c13_veg_cs, c13_veg_cf, dt)
        end if
        if ( use_c14 ) then
           call CarbonStateUpdate0(num_soilp, filter_soilp, c14_veg_cs, c14_veg_cf, dt)
        end if
        !#py call t_stopf('CNUpdate0')

        !--------------------------------------------
        if(use_pheno_flux_limiter)then
          !#py call t_startf('phenology_flux_limiter')
          call phenology_flux_limiter(bounds, num_soilc, filter_soilc,&
            num_soilp, filter_soilp, crop_vars, cnstate_vars,  &
            veg_cf, veg_cs, &
            c13_veg_cf, c13_veg_cs, &
            c14_veg_cf, c14_veg_cs, &
            veg_nf, veg_ns, veg_pf, veg_ps)
          !#py call t_stopf('phenology_flux_limiter')
        endif
        !#py call t_startf('CNLitterToColumn')
        call CNLitterToColumn(num_soilc, filter_soilc, &
                cnstate_vars)
        !#py call t_stopf('CNLitterToColumn')
        !--------------------------------------------
        ! Update1
        !--------------------------------------------

        !#py call t_startf('CNUpdate1')

        if ( use_c13 ) then
           call CarbonIsoFlux1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
                cnstate_vars, &
                isotope=c13, isocol_cs=c13_col_cs, isoveg_cs=c13_veg_cs, &
                isocol_cf=c13_col_cf, isoveg_cf=c13_veg_cf)
        end if

        if ( use_c14 ) then
           call CarbonIsoFlux1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
                cnstate_vars,  &
                isotope=c14, isocol_cs=c14_col_cs, isoveg_cs=c14_veg_cs,&
                 isocol_cf=c14_col_cf, isoveg_cf=c14_veg_cf)
        end if

        call CarbonStateUpdate1(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
             crop_vars, col_cs, veg_cs, col_cf, veg_cf, dt)

        if ( use_c13 ) then
           call CarbonStateUpdate1(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
                crop_vars, c13_col_cs, c13_veg_cs, c13_col_cf, c13_veg_cf, dt)
        end if
        if ( use_c14 ) then
           call CarbonStateUpdate1(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
                crop_vars, c14_col_cs, c14_veg_cs, c14_col_cf, c14_veg_cf, dt)
        end if

        call NitrogenStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
             cnstate_vars, dt)

        call PhosphorusStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
             cnstate_vars,dt)

        !#py call t_stopf('CNUpdate1')

        !#py call t_startf('SoilLittVertTransp')
        call SoilLittVertTransp(bounds, &
             num_soilc, filter_soilc, &
             canopystate_vars, cnstate_vars, &
             dt, year, mon, day, sec )
        !#py call t_stopf('SoilLittVertTransp')

        !#py call t_startf('CNGapMortality')
        call GapMortality( num_soilc, filter_soilc, num_soilp, filter_soilp, &
             cnstate_vars, dayspyr )
        !#py call t_stopf('CNGapMortality')

        !--------------------------------------------
        ! Update2
        !--------------------------------------------

        !#py call t_startf('CNUpdate2')

        if ( use_c13 ) then
           call CarbonIsoFlux2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
                cnstate_vars, &
                isotope=c13, isocol_cs=c13_col_cs, isoveg_cs=c13_veg_cs,&
                isocol_cf=c13_col_cf, isoveg_cf=c13_veg_cf)
        end if

        if ( use_c14 ) then
           call CarbonIsoFlux2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
                cnstate_vars, &
                isotope=c14, isocol_cs=c14_col_cs, isoveg_cs=c14_veg_cs, &
                isocol_cf=c14_col_cf, isoveg_cf=c14_veg_cf)
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
        call NitrogenStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, dt)

        call PhosphorusStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, dt)

        if (get_do_harvest()) then
           call CNHarvest(num_soilc, filter_soilc, num_soilp, filter_soilp, &
                cnstate_vars, dayspyr)
        end if

        if ( use_c13 ) then
           call CarbonIsoFlux2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
                cnstate_vars, &
                isotope=c13, isocol_cs=c13_col_cs, isoveg_cs=c13_veg_cs, &
                isocol_cf=c13_col_cf, isoveg_cf=c13_veg_cf)
        end if
        if ( use_c14 ) then
           call CarbonIsoFlux2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
                cnstate_vars,&
                isotope=c14, isocol_cs=c14_col_cs, &
                isoveg_cs=c14_veg_cs, isocol_cf=c14_col_cf, isoveg_cf=c14_veg_cf)
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

        call NitrogenStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, dt)

        call PhosphorusStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, dt)

        call WoodProducts(num_soilc, filter_soilc, dt)

        call CropHarvestPools(num_soilc, filter_soilc,dt)

        call FireArea(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
             atm2lnd_vars,  energyflux_vars, soilhydrology_vars, &
             cnstate_vars, dt, dayspyr, year, mon, day, sec, nstep)

        call FireFluxes(num_soilc, filter_soilc, num_soilp, filter_soilp, &
             cnstate_vars, dt, dayspyr, year, mon, day, sec)

        if ( use_erosion ) then
            call ErosionFluxes(bounds, num_soilc, filter_soilc, soilstate_vars, sedflux_vars, &
                 carbonstate_vars, nitrogenstate_vars, phosphorusstate_vars, carbonflux_vars, &
                 nitrogenflux_vars, phosphorusflux_vars)
        end if

        !#py call t_stopf('CNUpdate2')

        !--------------------------------------------
        ! Update3
        !--------------------------------------------

        if ( use_c13 ) then
           call CarbonIsoFlux3(num_soilc, filter_soilc, num_soilp, filter_soilp, &
                cnstate_vars, &
                isotope=c13, isocol_cs=c13_col_cs, isoveg_cs=c13_veg_cs,&
                isocol_cf=c13_col_cf, isoveg_cf=c13_veg_cf)
        end if
        if ( use_c14 ) then
           call CarbonIsoFlux3(num_soilc, filter_soilc, num_soilp, filter_soilp, &
                cnstate_vars, &
                isotope=c14, isocol_cs=c14_col_cs, isoveg_cs=c14_veg_cs,&
                isocol_cf=c14_col_cf, isoveg_cf=c14_veg_cf)
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
                cnstate_vars,dt, dayspyr, year, mon, day, tod, offset )

          call C14BombSpike(num_soilp, filter_soilp, &
               cnstate_vars,year, mon, day, tod, offset, dayspyr)
        end if

        call colcf_Summary_for_CH4_acc(col_cf,bounds, num_soilc, filter_soilc)
        call vegcf_summary_for_ch4_acc(veg_cf,bounds, num_soilp, filter_soilp)
        if( use_c13 ) then
           call colcf_Summary_for_CH4_acc(c13_col_cf,bounds, num_soilc, filter_soilc)
           call vegcf_summary_for_CH4_acc(c13_veg_cf,bounds, num_soilp, filter_soilp)
        endif
        if( use_c14 ) then
           call colcf_Summary_for_CH4_acc(c14_col_cf,bounds, num_soilc, filter_soilc)
           call vegcf_Summary_for_CH4_acc(c14_veg_cf,bounds, num_soilp, filter_soilp)
        endif

     end if !end of if not use_fates block


   end subroutine EcosystemDynNoLeaching2


end  module EcosystemDynMod
