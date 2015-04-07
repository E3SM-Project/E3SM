module CNDriverMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Ecosystem dynamics: phenology, vegetation
  !
  ! !USES:
  use shr_kind_mod                    , only : r8 => shr_kind_r8
  use clm_varctl                      , only : flanduse_timeseries, use_c13, use_c14, use_ed
  use decompMod                       , only : bounds_type
  use perf_mod                        , only : t_startf, t_stopf
  use clm_varctl                      , only : use_century_decomp, use_nitrif_denitrif
  use CNVegStateType                  , only : cnveg_state_type
  use CNVegCarbonStateType	      , only : cnveg_carbonstate_type
  use CNVegCarbonFluxType	      , only : cnveg_carbonflux_type
  use CNVegNitrogenStateType	      , only : cnveg_nitrogenstate_type
  use CNVegNitrogenFluxType	      , only : cnveg_nitrogenflux_type
  use SoilBiogeochemStateType         , only : soilbiogeochem_state_type
  use SoilBiogeochemCarbonStateType   , only : soilbiogeochem_carbonstate_type
  use SoilBiogeochemCarbonFluxType    , only : soilbiogeochem_carbonflux_type
  use SoilBiogeochemNitrogenStateType , only : soilbiogeochem_nitrogenstate_type
  use SoilBiogeochemNitrogenFluxType  , only : soilbiogeochem_nitrogenflux_type
  use CNDVType                        , only : dgvs_type
  use CanopyStateType                 , only : canopystate_type
  use SoilStateType                   , only : soilstate_type
  use TemperatureType                 , only : temperature_type
  use WaterstateType                  , only : waterstate_type
  use WaterfluxType                   , only : waterflux_type
  use atm2lndType                     , only : atm2lnd_type
  use SoilStateType                   , only : soilstate_type
  use CanopyStateType                 , only : canopystate_type
  use TemperatureType                 , only : temperature_type 
  use PhotosynthesisMod               , only : photosyns_type
  use ch4Mod                          , only : ch4_type
  use EnergyFluxType                  , only : energyflux_type
  use SoilHydrologyType               , only : soilhydrology_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CNDriverInit         ! Ecosystem dynamics: initialization
  public :: CNDriverNoLeaching   ! Ecosystem dynamics: phenology, vegetation, before doing N leaching
  public :: CNDriverLeaching     ! Ecosystem dynamics: phenology, vegetation, doing N leaching
  public :: CNDriverSummary      ! Ecosystem dynamics: summary
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CNDriverInit(bounds)
    !
    ! !DESCRIPTION:
    ! Initialzation of the CN Ecosystem dynamics.
    !
    ! !USES:
    use CNPhenologyMod              , only : CNPhenologyInit
    use CNFireMod                   , only : CNFireInit
    use SoilBiogeochemCompetitionMod, only : SoilBiogeochemCompetitionInit
    !
    ! !ARGUMENTS:
    type(bounds_type)                      , intent(in) :: bounds      
    !-----------------------------------------------------------------------

    call SoilBiogeochemCompetitionInit(bounds)
    call CNPhenologyInit(bounds)
    call CNFireInit(bounds)
    
  end subroutine CNDriverInit

  !-----------------------------------------------------------------------
  subroutine CNDriverNoLeaching(bounds,                                                    &
       num_soilc, filter_soilc, num_soilp, filter_soilp, num_pcropp, filter_pcropp, doalb, &
       cnveg_state_inst,                                                                   &
       cnveg_carbonflux_inst, cnveg_carbonstate_inst,                                      &
       c13_cnveg_carbonflux_inst, c13_cnveg_carbonstate_inst,                              &
       c14_cnveg_carbonflux_inst, c14_cnveg_carbonstate_inst,                              &
       cnveg_nitrogenflux_inst, cnveg_nitrogenstate_inst,                                  &
       soilbiogeochem_carbonflux_inst, soilbiogeochem_carbonstate_inst,                    &
       c13_soilbiogeochem_carbonflux_inst, c13_soilbiogeochem_carbonstate_inst,            &
       c14_soilbiogeochem_carbonflux_inst, c14_soilbiogeochem_carbonstate_inst,            &
       soilbiogeochem_state_inst,                                                          &
       soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst,                &
       atm2lnd_inst, waterstate_inst, waterflux_inst,                                      &
       canopystate_inst, soilstate_inst, temperature_inst, crop_inst, ch4_inst,            &
       dgvs_inst, photosyns_inst, soilhydrology_inst, energyflux_inst, nutrient_competition_method)
    !
    ! !DESCRIPTION:
    ! The core CN code is executed here. Calculates fluxes for maintenance
    ! respiration, decomposition, allocation, phenology, and growth respiration.
    ! These routines happen on the radiation time step so that canopy structure
    ! stays synchronized with albedo calculations.
    !
    ! !USES:
    use clm_varpar                        , only: crop_prog, nlevgrnd, nlevdecomp_full 
    use clm_varpar                        , only: nlevdecomp, ndecomp_cascade_transitions, ndecomp_pools
    use subgridAveMod                     , only: p2c
    use CropType                          , only: crop_type
    use CNNDynamicsMod                    , only: CNNDeposition,CNNFixation, CNNFert, CNSoyfix
    use CNMRespMod                        , only: CNMResp
    use CNPhenologyMod                    , only: CNPhenology
    use CNGRespMod                        , only: CNGResp
    use CNFireMod                         , only: CNFireArea, CNFireFluxes
    use CNCIsoFluxMod                     , only: CIsoFlux1, CIsoFlux2, CIsoFlux2h, CIsoFlux3
    use CNC14DecayMod                     , only: C14Decay
    use CNWoodProductsMod                 , only: CNWoodProducts
    use CNCStateUpdate1Mod                , only: CStateUpdate1,CStateUpdate0
    use CNCStateUpdate2Mod                , only: CStateUpdate2, CStateUpdate2h
    use CNCStateUpdate3Mod                , only: CStateUpdate3
    use CNNStateUpdate1Mod                , only: NStateUpdate1
    use CNNStateUpdate2Mod                , only: NStateUpdate2, NStateUpdate2h
    use CNGapMortalityMod                 , only: CNGapMortality
    use dynHarvestMod                     , only: CNHarvest
    use SoilBiogeochemDecompCascadeBGCMod , only: decomp_rate_constants_bgc
    use SoilBiogeochemDecompCascadeCNMod  , only: decomp_rate_constants_cn
    use SoilBiogeochemCompetitionMod      , only: SoilBiogeochemCompetition
    use SoilBiogeochemDecompMod           , only: SoilBiogeochemDecomp
    use SoilBiogeochemLittVertTranspMod   , only: SoilBiogeochemLittVertTransp
    use SoilBiogeochemPotentialMod        , only: SoilBiogeochemPotential 
    use SoilBiogeochemVerticalProfileMod  , only: SoilBiogeochemVerticalProfile
    use SoilBiogeochemNitrifDenitrifMod   , only: SoilBiogeochemNitrifDenitrif
    use SoilBiogeochemNStateUpdate1Mod    , only: SoilBiogeochemNStateUpdate1
    use NutrientCompetitionMethodMod      , only: nutrient_competition_method_type
    !
    ! !ARGUMENTS:
    type(bounds_type)                       , intent(in)    :: bounds  
    integer                                 , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                                 , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                                 , intent(in)    :: filter_soilp(:)   ! filter for soil patches
    integer                                 , intent(in)    :: num_pcropp        ! number of prog. crop patches in filter
    integer                                 , intent(in)    :: filter_pcropp(:)  ! filter for prognostic crop patches
    logical                                 , intent(in)    :: doalb             ! true = surface albedo calculation time step
    type(cnveg_state_type)                  , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonflux_type)             , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)            , intent(inout) :: cnveg_carbonstate_inst
    type(cnveg_carbonflux_type)             , intent(inout) :: c13_cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)            , intent(inout) :: c13_cnveg_carbonstate_inst
    type(cnveg_carbonflux_type)             , intent(inout) :: c14_cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)            , intent(inout) :: c14_cnveg_carbonstate_inst
    type(cnveg_nitrogenflux_type)           , intent(inout) :: cnveg_nitrogenflux_inst
    type(cnveg_nitrogenstate_type)          , intent(inout) :: cnveg_nitrogenstate_inst
    type(soilbiogeochem_state_type)         , intent(inout) :: soilbiogeochem_state_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: c13_soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: c13_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: c14_soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: c14_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst
    type(atm2lnd_type)                      , intent(in)    :: atm2lnd_inst 
    type(waterstate_type)                   , intent(in)    :: waterstate_inst
    type(waterflux_type)                    , intent(in)    :: waterflux_inst
    type(canopystate_type)                  , intent(in)    :: canopystate_inst
    type(soilstate_type)                    , intent(in)    :: soilstate_inst
    type(temperature_type)                  , intent(inout) :: temperature_inst
    type(crop_type)                         , intent(in)    :: crop_inst
    type(ch4_type)                          , intent(in)    :: ch4_inst
    type(dgvs_type)                         , intent(inout) :: dgvs_inst
    type(photosyns_type)                    , intent(in)    :: photosyns_inst
    type(soilhydrology_type)                , intent(in)    :: soilhydrology_inst
    type(energyflux_type)                   , intent(in)    :: energyflux_inst
    class(nutrient_competition_method_type) , intent(in)    :: nutrient_competition_method
    !
    ! !LOCAL VARIABLES:
    real(r8):: cn_decomp_pools(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_pools)
    real(r8):: p_decomp_cpool_loss(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_cascade_transitions) !potential C loss from one pool to another
    real(r8):: pmnf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_cascade_transitions) !potential mineral N flux, from one pool to another
    real(r8):: arepr(bounds%begp:bounds%endp) ! reproduction allocation coefficient (only used for crop_prog)
    real(r8):: aroot(bounds%begp:bounds%endp) ! root allocation coefficient (only used for crop_prog)
    integer :: begp,endp
    integer :: begc,endc
    !-----------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc

    !real(r8) , intent(in)    :: rootfr_patch(bounds%begp:, 1:)          
    !integer  , intent(in)    :: altmax_lastyear_indx_col(bounds%begc:)  ! frost table depth (m)

    associate(                                                                    &
         rootfr_patch              => soilstate_inst%rootfr_patch               , & ! fraction of roots in each soil layer  (nlevgrnd)
         altmax_lastyear_indx_col  => canopystate_inst%altmax_lastyear_indx_col , & ! frost table depth (m)
         laisun                    => canopystate_inst%laisun_patch             , & ! Input:  [real(r8) (:)   ]  sunlit projected leaf area index        
         laisha                    => canopystate_inst%laisha_patch             , & ! Input:  [real(r8) (:)   ]  shaded projected leaf area index        
         frac_veg_nosno            => canopystate_inst%frac_veg_nosno_patch     , & ! Input:  [integer  (:)   ]  fraction of vegetation not covered by snow (0 OR 1) [-]
         frac_veg_nosno_alb        => canopystate_inst%frac_veg_nosno_alb_patch , & ! Output: [integer  (:) ] frac of vegetation not covered by snow [-]         
         tlai                      => canopystate_inst%tlai_patch               , & ! Input:  [real(r8) (:) ]  one-sided leaf area index, no burying by snow     
         tsai                      => canopystate_inst%tsai_patch               , & ! Input:  [real(r8) (:)   ]  one-sided stem area index, no burying by snow     
         elai                      => canopystate_inst%elai_patch               , & ! Output: [real(r8) (:) ] one-sided leaf area index with burying by snow    
         esai                      => canopystate_inst%esai_patch               , & ! Output: [real(r8) (:) ] one-sided stem area index with burying by snow    
         htop                      => canopystate_inst%htop_patch               , & ! Output: [real(r8) (:) ] canopy top (m)                                     
         hbot                      => canopystate_inst%hbot_patch                 & ! Output: [real(r8) (:) ] canopy bottom (m)                                  
      )

    ! --------------------------------------------------
    ! zero the column-level C and N fluxes
    ! --------------------------------------------------
    
    call t_startf('CNZero')

    call cnveg_carbonflux_inst%SetValues( &
         num_soilp, filter_soilp, 0._r8, &
         num_soilc, filter_soilc, 0._r8)
    if ( use_c13 ) then
       call c13_cnveg_carbonflux_inst%SetValues( &
            num_soilp, filter_soilp, 0._r8, &
            num_soilc, filter_soilc, 0._r8)
    end if
    if ( use_c14 ) then
       call c14_cnveg_carbonflux_inst%SetValues( &
            num_soilp, filter_soilp, 0._r8, &
            num_soilc, filter_soilc, 0._r8)
    end if

    call soilbiogeochem_carbonflux_inst%SetValues( &
         num_soilc, filter_soilc, 0._r8)
    if ( use_c13 ) then
       call c13_soilbiogeochem_carbonflux_inst%SetValues( &
            num_soilc, filter_soilc, 0._r8)
    end if
    if ( use_c14 ) then
       call c14_soilbiogeochem_carbonflux_inst%SetValues( &
            num_soilc, filter_soilc, 0._r8)
    end if

    call cnveg_carbonflux_inst%SetValues( &
         num_soilp, filter_soilp, 0._r8, &
         num_soilc, filter_soilc, 0._r8)
    if ( use_c13 ) then
       call c13_cnveg_carbonflux_inst%SetValues( &
            num_soilp, filter_soilp, 0._r8, &
            num_soilc, filter_soilc, 0._r8)
    end if
    if ( use_c14 ) then
       call c14_cnveg_carbonflux_inst%SetValues( &
            num_soilp, filter_soilp, 0._r8, &
            num_soilc, filter_soilc, 0._r8)
    end if

    call cnveg_nitrogenflux_inst%SetValues( &
         num_soilp, filter_soilp, 0._r8, &
         num_soilc, filter_soilc, 0._r8)

    call soilbiogeochem_nitrogenflux_inst%SetValues( &
         num_soilc, filter_soilc, 0._r8)

    call t_stopf('CNZero')

    ! --------------------------------------------------
    ! Nitrogen Deposition, Fixation and Respiration
    ! --------------------------------------------------

    call t_startf('CNDeposition')
    call CNNDeposition(bounds, &
         atm2lnd_inst, soilbiogeochem_nitrogenflux_inst)
    call t_stopf('CNDeposition')

    call t_startf('CNFixation')
    call CNNFixation( num_soilc, filter_soilc, &
         cnveg_carbonflux_inst, soilbiogeochem_nitrogenflux_inst)
    call t_stopf('CNFixation')

    if (crop_prog) then
       call CNNFert(bounds, num_soilc,filter_soilc, &
            cnveg_nitrogenflux_inst, soilbiogeochem_nitrogenflux_inst)

       call  CNSoyfix (bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            waterstate_inst, crop_inst, cnveg_state_inst, cnveg_nitrogenflux_inst , &
            soilbiogeochem_state_inst, soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst)
    end if

    call t_startf('CNMResp')
    call CNMResp(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
         canopystate_inst, soilstate_inst, temperature_inst, photosyns_inst, &
         cnveg_carbonflux_inst, cnveg_nitrogenstate_inst)
    call t_stopf('CNMResp')

    !--------------------------------------------
    ! Soil Biogeochemistry
    !--------------------------------------------

    if (use_century_decomp) then
       call decomp_rate_constants_bgc(bounds, num_soilc, filter_soilc, &
            canopystate_inst, soilstate_inst, temperature_inst, ch4_inst, soilbiogeochem_carbonflux_inst)
    else
       call decomp_rate_constants_cn(bounds, num_soilc, filter_soilc, &
            canopystate_inst, soilstate_inst, temperature_inst, ch4_inst, soilbiogeochem_carbonflux_inst)
    end if

    ! calculate potential decomp rates and total immobilization demand (previously inlined in CNDecompAlloc)
    call SoilBiogeochemPotential (bounds, num_soilc, filter_soilc,                                                    &
         soilbiogeochem_state_inst, soilbiogeochem_carbonstate_inst, soilbiogeochem_carbonflux_inst,                  &
         soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst,                                         &
         cn_decomp_pools=cn_decomp_pools(begc:endc,1:nlevdecomp,1:ndecomp_pools), & 
         p_decomp_cpool_loss=p_decomp_cpool_loss(begc:endc,1:nlevdecomp,1:ndecomp_cascade_transitions), &
         pmnf_decomp_cascade=pmnf_decomp_cascade(begc:endc,1:nlevdecomp,1:ndecomp_cascade_transitions)) 

    ! calculate vertical profiles for distributing soil and litter C and N (previously subroutine decomp_vertprofiles called from CNDecompAlloc)
    call SoilBiogeochemVerticalProfile(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
         canopystate_inst, soilstate_inst,soilbiogeochem_state_inst)

    ! calculate nitrification and denitrification rates (previously subroutine nitrif_denitrif called from CNDecompAlloc)
    if (use_nitrif_denitrif) then 
       call SoilBiogeochemNitrifDenitrif(bounds, num_soilc, filter_soilc, &
            soilstate_inst, waterstate_inst, temperature_inst, ch4_inst, &
            soilbiogeochem_carbonflux_inst, soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst)
    end if

    !--------------------------------------------
    ! Resolve the competition between plants and soil heterotrophs 
    ! for available soil mineral N resource 
    !--------------------------------------------

    call t_startf('CNDecompAlloc')

    ! Jinyun Tang: at this stage, the plant_nutrient_demand only calculates the plant ntirgeon demand.
    ! Assume phosphorus dynamics will be included in the future. Also, I consider plant_nutrient_demand
    ! as a generic interface to call actual nutrient calculation from different aboveground plantbgc. 
    ! Right now it is assumed the plant nutrient demand is summarized into columnwise demand, and the 
    ! nutrient redistribution after uptake is done by the plant bgc accordingly. 
    ! When nutrient competition is required to be done at cohort level both plant_nutrient_demand and 
    ! do_nutrient_competition should be modified, but that modification should not significantly change 
    ! the current interface.

    call nutrient_competition_method%calc_plant_nutrient_demand (bounds,  &
         num_soilp, filter_soilp,                                         &
         photosyns_inst, crop_inst, canopystate_inst,                     &
         cnveg_state_inst, cnveg_carbonstate_inst, cnveg_carbonflux_inst, &
         c13_cnveg_carbonflux_inst, c14_cnveg_carbonflux_inst,            &
         cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst,               &
         aroot=aroot(begp:endp), arepr=arepr(begp:endp))

    ! get the column-averaged plant_ndemand (needed for following call to SoilBiogeochemCompetition)

    call p2c(bounds, num_soilc, filter_soilc,                    &
         cnveg_nitrogenflux_inst%plant_ndemand_patch(begp:endp), &
         soilbiogeochem_state_inst%plant_ndemand_col(begc:endc))

    ! resolve plant/heterotroph competition for mineral N 

    call SoilBiogeochemCompetition (bounds, num_soilc, filter_soilc, &
         soilbiogeochem_state_inst, soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst)

    ! distribute the available N between the competing patches  on the basis of 
    ! relative demand, and allocate C and N to new growth and storage

    call nutrient_competition_method%calc_plant_nutrient_competition (bounds,&
         num_soilp, filter_soilp,                                            &
         cnveg_state_inst, cnveg_carbonflux_inst, c13_cnveg_carbonflux_inst, &
         c14_cnveg_carbonflux_inst, cnveg_nitrogenflux_inst,                 &
         aroot=aroot(begp:endp),                                             &
         arepr=arepr(begp:endp),                                             &
         fpg_col=soilbiogeochem_state_inst%fpg_col(begc:endc))

    call t_stopf('CNDecompAlloc')

    !--------------------------------------------
    ! Calculate litter and soil decomposition rate
    !--------------------------------------------

    ! Calculation of actual immobilization and decomp rates, following
    ! resolution of plant/heterotroph  competition for mineral N (previously inlined in CNDecompAllocation in CNDecompMod)

    call t_startf('SoilBiogeochemDecomp')

    call SoilBiogeochemDecomp (bounds, num_soilc, filter_soilc,                                                       &
         soilbiogeochem_state_inst, soilbiogeochem_carbonstate_inst, soilbiogeochem_carbonflux_inst,                  &
         soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst,                                         &
         cn_decomp_pools=cn_decomp_pools(begc:endc,1:nlevdecomp,1:ndecomp_pools),                       & 
         p_decomp_cpool_loss=p_decomp_cpool_loss(begc:endc,1:nlevdecomp,1:ndecomp_cascade_transitions), &
         pmnf_decomp_cascade=pmnf_decomp_cascade(begc:endc,1:nlevdecomp,1:ndecomp_cascade_transitions)) 

    call t_stopf('SoilBiogeochemDecomp')

    !--------------------------------------------
    ! Phenology
    !--------------------------------------------

    ! CNphenology needs to be called after above calls, since it depends on current
    ! time-step fluxes to new growth on the lastlitterfall timestep in deciduous systems

    call t_startf('CNPhenology')

    call CNPhenology (bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, num_pcropp, filter_pcropp, &
         doalb, waterstate_inst, temperature_inst, crop_inst, canopystate_inst, soilstate_inst, dgvs_inst, &
         cnveg_state_inst, cnveg_carbonstate_inst, cnveg_carbonflux_inst,                                  &
         cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst,                                                &
         leaf_prof_patch=soilbiogeochem_state_inst%leaf_prof_patch(begp:endp,1:nlevdecomp_full), &
         froot_prof_patch=soilbiogeochem_state_inst%froot_prof_patch(begp:endp,1:nlevdecomp_full))

    call t_stopf('CNPhenology')

    !--------------------------------------------
    ! Growth respiration
    !--------------------------------------------

    call t_startf('CNGResp')

    call CNGResp(num_soilp, filter_soilp,&
         cnveg_carbonflux_inst)

    call t_stopf('CNGResp')

    !--------------------------------------------
    ! CNUpdate0
    !--------------------------------------------

    call t_startf('CNUpdate0')

    call CStateUpdate0(num_soilp, filter_soilp, &
         cnveg_carbonflux_inst, cnveg_carbonstate_inst)

    if ( use_c13 ) then
       call CStateUpdate0(num_soilp, filter_soilp, &
            c13_cnveg_carbonflux_inst, c13_cnveg_carbonstate_inst)
    end if

    if ( use_c14 ) then
       call CStateUpdate0(num_soilp, filter_soilp, &
            c14_cnveg_carbonflux_inst, c14_cnveg_carbonstate_inst)
    end if

    call t_stopf('CNUpdate0')

    !--------------------------------------------
    ! Update1
    !--------------------------------------------

    call t_startf('CNUpdate1')

    ! Set the carbon isotopic flux variables (except for gap-phase mortality and fire fluxes)
    if ( use_c13 ) then
       call CIsoFlux1(num_soilc, filter_soilc, num_soilp, filter_soilp,              &
            soilbiogeochem_state_inst,                                               &
            soilbiogeochem_carbonflux_inst,  soilbiogeochem_carbonstate_inst,        &
            cnveg_carbonflux_inst, cnveg_carbonstate_inst,                           &
            c13_soilbiogeochem_carbonflux_inst, c13_soilbiogeochem_carbonstate_inst, &
            c13_cnveg_carbonflux_inst, c13_cnveg_carbonstate_inst,                   &
            isotope='c13')
    end if
    if ( use_c14 ) then
       call CIsoFlux1(num_soilc, filter_soilc, num_soilp, filter_soilp,              &
            soilbiogeochem_state_inst,                                               &
            soilbiogeochem_carbonflux_inst,  soilbiogeochem_carbonstate_inst,        &
            cnveg_carbonflux_inst, cnveg_carbonstate_inst,                           &
            c14_soilbiogeochem_carbonflux_inst, c14_soilbiogeochem_carbonstate_inst, &
            c14_cnveg_carbonflux_inst, c14_cnveg_carbonstate_inst,                   &
            isotope='c14')
    end if

    ! Update all prognostic carbon state variables (except for gap-phase mortality and fire fluxes)
    call CStateUpdate1( num_soilc, filter_soilc, num_soilp, filter_soilp, &
         cnveg_state_inst, cnveg_carbonflux_inst, cnveg_carbonstate_inst, &
         soilbiogeochem_carbonflux_inst)
    if ( use_c13 ) then
       call CStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnveg_state_inst, c13_cnveg_carbonflux_inst, c13_cnveg_carbonstate_inst, &
            c13_soilbiogeochem_carbonflux_inst)
    end if
    if ( use_c14 ) then
       call CStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnveg_state_inst, c14_cnveg_carbonflux_inst, c14_cnveg_carbonstate_inst, &
            c14_soilbiogeochem_carbonflux_inst)
    end if

    ! Update all prognostic nitrogen state variables (except for gap-phase mortality and fire fluxes)
    call NStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
         cnveg_nitrogenflux_inst, cnveg_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst) 

    call SoilBiogeochemNStateUpdate1(num_soilc, filter_soilc,  &
         soilbiogeochem_state_inst, soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst)

    call t_stopf('CNUpdate1')

    !--------------------------------------------
    ! Calculate vertical mixing of soil and litter pools
    !--------------------------------------------

    call t_startf('SoilBiogeochemLittVertTransp')

    call SoilBiogeochemLittVertTransp(bounds, num_soilc, filter_soilc,            &
         canopystate_inst, soilbiogeochem_state_inst,                             &
         soilbiogeochem_carbonstate_inst, soilbiogeochem_carbonflux_inst,         &
         c13_soilbiogeochem_carbonstate_inst, c13_soilbiogeochem_carbonflux_inst, &
         c14_soilbiogeochem_carbonstate_inst, c14_soilbiogeochem_carbonflux_inst, &
         soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst)

    call t_stopf('SoilBiogeochemLittVertTransp')

    !--------------------------------------------
    ! Calculate the gap mortality carbon and nitrogen fluxes
    !--------------------------------------------

    call t_startf('CNGapMortality')

    call CNGapMortality (bounds, num_soilc, filter_soilc, num_soilp, filter_soilp,                                &
         dgvs_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst,                                             &
         cnveg_carbonflux_inst, cnveg_nitrogenflux_inst,                                                          &
         leaf_prof_patch=soilbiogeochem_state_inst%leaf_prof_patch(begp:endp, 1:nlevdecomp_full),   &
         froot_prof_patch=soilbiogeochem_state_inst%froot_prof_patch(begp:endp, 1:nlevdecomp_full), & 
         croot_prof_patch=soilbiogeochem_state_inst%croot_prof_patch(begp:endp, 1:nlevdecomp_full), &
         stem_prof_patch=soilbiogeochem_state_inst%stem_prof_patch(begp:endp, 1:nlevdecomp_full))

    call t_stopf('CNGapMortality')

    !--------------------------------------------
    ! Update2 (gap mortality)
    !--------------------------------------------

    call t_startf('CNUpdate2')

    ! Set the carbon isotopic fluxes for gap mortality
    if ( use_c13 ) then
       call CIsoFlux2(num_soilc, filter_soilc, num_soilp, filter_soilp,               &
            soilbiogeochem_state_inst, cnveg_carbonflux_inst, cnveg_carbonstate_inst, &
            iso_cnveg_carbonflux_inst=c13_cnveg_carbonflux_inst,                      &
            iso_cnveg_carbonstate_inst=c13_cnveg_carbonstate_inst,                    &
            isotope='c13')
    end if
    if ( use_c14 ) then
       call CIsoFlux2(num_soilc, filter_soilc, num_soilp, filter_soilp,               &
            soilbiogeochem_state_inst, cnveg_carbonflux_inst, cnveg_carbonstate_inst, &
            iso_cnveg_carbonflux_inst=c14_cnveg_carbonflux_inst,                      &
            iso_cnveg_carbonstate_inst=c14_cnveg_carbonstate_inst,                    &
            isotope='c14')
    end if

    ! Update all the prognostic carbon state variables affected by gap-phase mortality fluxes
    call CStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
         cnveg_carbonflux_inst, cnveg_carbonstate_inst, soilbiogeochem_carbonstate_inst)
    if ( use_c13 ) then
       call CStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            c13_cnveg_carbonflux_inst, c13_cnveg_carbonstate_inst, c13_soilbiogeochem_carbonstate_inst)
    end if
    if ( use_c14 ) then
       call CStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            c14_cnveg_carbonflux_inst, c14_cnveg_carbonstate_inst, c14_soilbiogeochem_carbonstate_inst)
    end if

    ! Update all the prognostic nitrogen state variables affected by gap-phase mortality fluxes
    call NStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
         cnveg_nitrogenflux_inst, cnveg_nitrogenstate_inst, soilbiogeochem_nitrogenstate_inst)

    !--------------------------------------------
    ! Update2h (harvest)
    !--------------------------------------------

    ! Set harvest mortality routine 
    if (flanduse_timeseries /= ' ') then
       call CNHarvest(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            soilbiogeochem_state_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, &
            cnveg_carbonflux_inst, cnveg_nitrogenflux_inst)
    end if

    if ( use_c13 ) then
       call CIsoFlux2h(num_soilc, filter_soilc, num_soilp, filter_soilp,   &
            soilbiogeochem_state_inst,                                     &
            cnveg_carbonflux_inst, cnveg_carbonstate_inst,                 &
            c13_cnveg_carbonflux_inst, c13_cnveg_carbonstate_inst,         &                         
            isotope='c13')
    end if
    if ( use_c14 ) then
       call CIsoFlux2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            soilbiogeochem_state_inst,                                     &
            cnveg_carbonflux_inst, cnveg_carbonstate_inst,                 &
            c14_cnveg_carbonflux_inst, c14_cnveg_carbonstate_inst,         &                         
            isotope='c14')
    end if

    call CStateUpdate2h( num_soilc, filter_soilc,  num_soilp, filter_soilp, &
         cnveg_carbonflux_inst, cnveg_carbonstate_inst, soilbiogeochem_carbonstate_inst)
    if ( use_c13 ) then
       call CStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            c13_cnveg_carbonflux_inst, c13_cnveg_carbonstate_inst, c13_soilbiogeochem_carbonstate_inst)
    end if
    if ( use_c14 ) then
       call CStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            c14_cnveg_carbonflux_inst, c14_cnveg_carbonstate_inst, c14_soilbiogeochem_carbonstate_inst)
    end if

    call NStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
         cnveg_nitrogenflux_inst, cnveg_nitrogenstate_inst, soilbiogeochem_nitrogenstate_inst)

    !--------------------------------------------
    ! Calculate loss fluxes from wood products pools
    ! and update product pool state variables
    !--------------------------------------------

    call CNWoodProducts(num_soilc, filter_soilc, &
         cnveg_carbonstate_inst, c13_cnveg_carbonstate_inst, c14_cnveg_carbonstate_inst, &
         cnveg_carbonflux_inst, c13_cnveg_carbonflux_inst, c14_cnveg_carbonflux_inst, &
         cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst)

    !--------------------------------------------
    ! Calculate fire area and fluxes
    !--------------------------------------------

    call CNFireArea(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
         atm2lnd_inst, energyflux_inst, soilhydrology_inst, waterstate_inst, &
         cnveg_state_inst, cnveg_carbonstate_inst, &
         totlitc_col=soilbiogeochem_carbonstate_inst%totlitc_col(begc:endc), &
         decomp_cpools_vr_col=soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools), &
         t_soi17cm_col=temperature_inst%t_soi17cm_col(begc:endc))

    call CNFireFluxes(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp,                                                    &
         dgvs_inst, cnveg_state_inst,                                                                                              &
         cnveg_carbonstate_inst, cnveg_carbonflux_inst, cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst,                         &
         leaf_prof_patch=soilbiogeochem_state_inst%leaf_prof_patch(begp:endp, 1:nlevdecomp_full),                                  &
         froot_prof_patch=soilbiogeochem_state_inst%froot_prof_patch(begp:endp, 1:nlevdecomp_full),                                & 
         croot_prof_patch=soilbiogeochem_state_inst%croot_prof_patch(begp:endp, 1:nlevdecomp_full),                                &
         stem_prof_patch=soilbiogeochem_state_inst%stem_prof_patch(begp:endp, 1:nlevdecomp_full),                                  &
         totsomc_col=soilbiogeochem_carbonstate_inst%totsomc_col(begc:endc),                                                       &
         decomp_cpools_vr_col=soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools),   &
         decomp_npools_vr_col=soilbiogeochem_nitrogenstate_inst%decomp_npools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools), &
         somc_fire_col=soilbiogeochem_carbonflux_inst%somc_fire_col(begc:endc))

    call t_stopf('CNUpdate2')

    !--------------------------------------------
    ! Update3
    !--------------------------------------------

    if ( use_c13 ) then
       call CIsoFlux3(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            soilbiogeochem_state_inst , soilbiogeochem_carbonstate_inst,         &
            cnveg_carbonflux_inst, cnveg_carbonstate_inst,                       &
            c13_cnveg_carbonflux_inst, c13_cnveg_carbonstate_inst,               &
            c13_soilbiogeochem_carbonstate_inst, &
            isotope='c13')
    end if
    if ( use_c14 ) then
       call CIsoFlux3(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            soilbiogeochem_state_inst , soilbiogeochem_carbonstate_inst,         &
            cnveg_carbonflux_inst, cnveg_carbonstate_inst,                       &
            c14_cnveg_carbonflux_inst, c14_cnveg_carbonstate_inst,               &
            c14_soilbiogeochem_carbonstate_inst, &
            isotope='c14')
    end if

    call CStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
         cnveg_carbonflux_inst, cnveg_carbonstate_inst, soilbiogeochem_carbonstate_inst)

    if ( use_c13 ) then
       call CStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
            c13_cnveg_carbonflux_inst, c13_cnveg_carbonstate_inst, c13_soilbiogeochem_carbonstate_inst)
    end if

    if ( use_c14 ) then
       call CStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
            c14_cnveg_carbonflux_inst, c14_cnveg_carbonstate_inst, c14_soilbiogeochem_carbonstate_inst)

       call C14Decay(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            c14_cnveg_carbonstate_inst, c14_soilbiogeochem_carbonstate_inst)
    end if

    end associate

  end subroutine CNDriverNoLeaching
  
  !-----------------------------------------------------------------------
  subroutine CNDriverLeaching(bounds, &
       num_soilc, filter_soilc, num_soilp, filter_soilp, &
       waterstate_inst, waterflux_inst, &
       cnveg_nitrogenflux_inst, cnveg_nitrogenstate_inst, &
       soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst)
    !
    ! !DESCRIPTION:
    ! Update the nitrogen leaching rate as a function of soluble mineral N and total soil water outflow.
    ! Also update nitrogen state variables         
    !
    ! !USES:
    use SoilBiogeochemNLeachingMod, only: SoilBiogeochemNLeaching
    use CNNStateUpdate3Mod   , only: NStateUpdate3
    !
    ! !ARGUMENTS:
    type(bounds_type)                       , intent(in)    :: bounds  
    integer                                 , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                                 , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                                 , intent(in)    :: filter_soilp(:)   ! filter for soil patches
    type(waterstate_type)                   , intent(in)    :: waterstate_inst
    type(waterflux_type)                    , intent(in)    :: waterflux_inst
    type(cnveg_nitrogenflux_type)           , intent(inout) :: cnveg_nitrogenflux_inst
    type(cnveg_nitrogenstate_type)          , intent(inout) :: cnveg_nitrogenstate_inst
    type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst
    !-----------------------------------------------------------------------
  
    ! Mineral nitrogen dynamics (deposition, fixation, leaching)
    
    call SoilBiogeochemNLeaching(bounds, num_soilc, filter_soilc, &
         waterstate_inst, waterflux_inst, soilbiogeochem_nitrogenstate_inst, &
         soilbiogeochem_nitrogenflux_inst)

    ! Nitrogen state variable update, mortality fluxes.

    call t_startf('CNUpdate3')

    call NstateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp, &
         cnveg_nitrogenflux_inst, cnveg_nitrogenstate_inst, &
         soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst)

    call t_stopf('CNUpdate3')

  end subroutine CNDriverLeaching

  !-----------------------------------------------------------------------
  subroutine CNDriverSummary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       cnveg_state_inst, cnveg_carbonflux_inst, cnveg_carbonstate_inst, &
       c13_cnveg_carbonflux_inst, c13_cnveg_carbonstate_inst, &
       c14_cnveg_carbonflux_inst, c14_cnveg_carbonstate_inst, &
       cnveg_nitrogenflux_inst, cnveg_nitrogenstate_inst, &
       soilbiogeochem_carbonflux_inst, soilbiogeochem_carbonstate_inst, &
       c13_soilbiogeochem_carbonflux_inst, c13_soilbiogeochem_carbonstate_inst, &
       c14_soilbiogeochem_carbonflux_inst, c14_soilbiogeochem_carbonstate_inst, &
       soilbiogeochem_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst)
    !
    ! !DESCRIPTION:
    ! Call to all CN and SoilBiogeochem summary routines
    !
    ! !USES:
    use clm_varpar                        , only: ndecomp_cascade_transitions
    use CNPrecisionControlMod             , only: CNPrecisionControl
    use SoilBiogeochemPrecisionControlMod , only: SoilBiogeochemPrecisionControl
    !
    ! !ARGUMENTS:
    type(bounds_type)                       , intent(in)    :: bounds  
    integer                                 , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                                 , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                                 , intent(in)    :: filter_soilp(:)   ! filter for soil patches
    type(cnveg_state_type)                  , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonflux_type)             , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)            , intent(inout) :: cnveg_carbonstate_inst
    type(cnveg_carbonflux_type)             , intent(inout) :: c13_cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)            , intent(inout) :: c13_cnveg_carbonstate_inst
    type(cnveg_carbonflux_type)             , intent(inout) :: c14_cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)            , intent(inout) :: c14_cnveg_carbonstate_inst
    type(cnveg_nitrogenflux_type)           , intent(inout) :: cnveg_nitrogenflux_inst
    type(cnveg_nitrogenstate_type)          , intent(inout) :: cnveg_nitrogenstate_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: c13_soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: c13_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: c14_soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_carbonstate_type)   , intent(inout) :: c14_soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst
    !
    ! !LOCAL VARIABLES:
    integer :: begc,endc
    !-----------------------------------------------------------------------
  
    begc = bounds%begc; endc= bounds%endc

    ! Call to all summary routines

    call t_startf('CNsum')

    ! Set controls on very low values in critical state variables 

    call CNPrecisionControl(num_soilp, filter_soilp, &
         cnveg_carbonstate_inst, c13_cnveg_carbonstate_inst, &
         c14_cnveg_carbonstate_inst, cnveg_nitrogenstate_inst)

    call SoilBiogeochemPrecisionControl(num_soilc, filter_soilc,  &
         soilbiogeochem_carbonstate_inst, c13_soilbiogeochem_carbonstate_inst, &
         c14_soilbiogeochem_carbonstate_inst,soilbiogeochem_nitrogenstate_inst)

    ! Note - all summary updates to cnveg_carbonstate_inst and cnveg_carbonflux_inst are done in 
    ! soilbiogeochem_carbonstate_inst%summary and CNVeg_carbonstate_inst%summary

    ! ----------------------------------------------
    ! soilbiogeochem carbon/nitrogen state summary
    ! ----------------------------------------------

    call soilbiogeochem_carbonstate_inst%summary(bounds, num_soilc, filter_soilc)
    if ( use_c13 ) then
       call c13_soilbiogeochem_carbonstate_inst%summary(bounds, num_soilc, filter_soilc)
    end if
    if ( use_c14 ) then
       call c14_soilbiogeochem_carbonstate_inst%summary(bounds, num_soilc, filter_soilc)
    end if
    call soilbiogeochem_nitrogenstate_inst%summary(bounds, num_soilc, filter_soilc)

    ! ----------------------------------------------
    ! soilbiogeochem carbon/nitrogen flux summary
    ! ----------------------------------------------

    call soilbiogeochem_carbonflux_inst%Summary(bounds, num_soilc, filter_soilc)
    if ( use_c13 ) then
       call c13_soilbiogeochem_carbonflux_inst%Summary(bounds, num_soilc, filter_soilc)
    end if
    if ( use_c14 ) then
       call c14_soilbiogeochem_carbonflux_inst%Summary(bounds, num_soilc, filter_soilc)
    end if
    call soilbiogeochem_nitrogenflux_inst%Summary(bounds, num_soilc, filter_soilc)

    ! ----------------------------------------------
    ! cnveg carbon/nitrogen state summary
    ! ----------------------------------------------

    call cnveg_carbonstate_inst%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
         soilbiogeochem_cwdc_col=soilbiogeochem_carbonstate_inst%cwdc_col(begc:endc), &
         soilbiogeochem_totlitc_col=soilbiogeochem_carbonstate_inst%totlitc_col(begc:endc), &
         soilbiogeochem_totsomc_col=soilbiogeochem_carbonstate_inst%totsomc_col(begc:endc), &
         soilbiogeochem_ctrunc_col=soilbiogeochem_carbonstate_inst%ctrunc_col(begc:endc))

    if ( use_c13 ) then
       call c13_cnveg_carbonstate_inst%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            soilbiogeochem_cwdc_col=c13_soilbiogeochem_carbonstate_inst%cwdc_col(begc:endc), &
            soilbiogeochem_totlitc_col=c13_soilbiogeochem_carbonstate_inst%totlitc_col(begc:endc), &
            soilbiogeochem_totsomc_col=c13_soilbiogeochem_carbonstate_inst%totsomc_col(begc:endc), &
            soilbiogeochem_ctrunc_col=c13_soilbiogeochem_carbonstate_inst%ctrunc_col(begc:endc))
    end if

    if ( use_c14 ) then
       call c14_cnveg_carbonstate_inst%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            soilbiogeochem_cwdc_col=c14_soilbiogeochem_carbonstate_inst%cwdc_col(begc:endc), &
            soilbiogeochem_totlitc_col=c14_soilbiogeochem_carbonstate_inst%totlitc_col(begc:endc), &
            soilbiogeochem_totsomc_col=c14_soilbiogeochem_carbonstate_inst%totsomc_col(begc:endc), &
            soilbiogeochem_ctrunc_col=c14_soilbiogeochem_carbonstate_inst%ctrunc_col(begc:endc))
    end if

    call cnveg_nitrogenstate_inst%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
         soilbiogeochem_nitrogenstate_inst)

    ! ----------------------------------------------
    ! cnveg carbon/nitrogen flux summary
    ! ----------------------------------------------

    call cnveg_carbonflux_inst%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
         isotope='bulk', &
         soilbiogeochem_hr_col=soilbiogeochem_carbonflux_inst%hr_col(begc:endc), &
         soilbiogeochem_lithr_col=soilbiogeochem_carbonflux_inst%lithr_col(begc:endc), &  
         soilbiogeochem_decomp_cascade_ctransfer_col=&
         soilbiogeochem_carbonflux_inst%decomp_cascade_ctransfer_col(begc:endc,1:ndecomp_cascade_transitions))

    if ( use_c13 ) then
       call c13_cnveg_carbonflux_inst%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            isotope='c13', &
            soilbiogeochem_hr_col=c13_soilbiogeochem_carbonflux_inst%hr_col(begc:endc), &
            soilbiogeochem_lithr_col=c13_soilbiogeochem_carbonflux_inst%lithr_col(begc:endc), &  
            soilbiogeochem_decomp_cascade_ctransfer_col=&
            c13_soilbiogeochem_carbonflux_inst%decomp_cascade_ctransfer_col(begc:endc,1:ndecomp_cascade_transitions))
    end if

    if ( use_c14 ) then
       call c14_cnveg_carbonflux_inst%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            isotope='c14', &
            soilbiogeochem_hr_col=c14_soilbiogeochem_carbonflux_inst%hr_col(begc:endc), &
            soilbiogeochem_lithr_col=c14_soilbiogeochem_carbonflux_inst%lithr_col(begc:endc), &  
            soilbiogeochem_decomp_cascade_ctransfer_col=&
            c14_soilbiogeochem_carbonflux_inst%decomp_cascade_ctransfer_col(begc:endc,1:ndecomp_cascade_transitions))
    end if

    call cnveg_nitrogenflux_inst%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)

    call t_stopf('CNsum')

  end subroutine CNDriverSummary

end  module CNDriverMod
