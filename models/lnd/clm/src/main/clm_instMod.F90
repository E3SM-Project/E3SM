module clm_instMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Instances and definitions of all data types
  !
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use decompMod       , only : bounds_type
  use clm_varpar      , only : crop_prog, ndecomp_pools, nlevdecomp_full
  use clm_varctl      , only : use_cn, use_c13, use_c14, use_lch4, use_cndv, use_ed, use_voc
  use clm_varctl      , only : use_century_decomp
  use clm_varcon      , only : h2osno_max, bdsno, c13ratio, c14ratio
  use landunit_varcon , only : istice, istice_mec, istsoil
  use perf_mod        , only : t_startf, t_stopf

  !-----------------------------------------
  ! Constants
  !-----------------------------------------

  use UrbanParamsType                    , only : urbanparams_type   ! Constants 
  use UrbanParamsType                    , only : IsSimpleBuildTemp, IsProgBuildTemp
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con
  use CNDVType                           , only : dgv_ecophyscon     ! Constants 
  use EDEcophysConType                   , only : EDecophyscon       ! ED Constants

  !-----------------------------------------
  ! Definition of component types 
  !-----------------------------------------

  use AerosolMod                      , only : aerosol_type
  use CanopyStateType                 , only : canopystate_type
  use ch4Mod                          , only : ch4_type
  use CNVegStateType                  , only : cnveg_state_type
  use CNVegCarbonFluxType             , only : cnveg_carbonflux_type
  use CNVegCarbonStateType            , only : cnveg_carbonstate_type
  use CNVegNitrogenFluxType           , only : cnveg_nitrogenflux_type
  use CNVegNitrogenStateType          , only : cnveg_nitrogenstate_type
  use SoilBiogeochemStateType         , only : soilbiogeochem_state_type
  use SoilBiogeochemCarbonFluxType    , only : soilbiogeochem_carbonflux_type
  use SoilBiogeochemCarbonStateType   , only : soilbiogeochem_carbonstate_type
  use SoilBiogeochemNitrogenFluxType  , only : soilbiogeochem_nitrogenflux_type
  use SoilBiogeochemNitrogenStateType , only : soilbiogeochem_nitrogenstate_type
  use CNDVType                        , only : dgvs_type
  use CropType                        , only : crop_type
  use DryDepVelocity                  , only : drydepvel_type
  use DUSTMod                         , only : dust_type
  use EnergyFluxType                  , only : energyflux_type
  use FrictionVelocityMod             , only : frictionvel_type
  use IrrigationMod                   , only : irrigation_type
  use LakeStateType                   , only : lakestate_type
  use PhotosynthesisMod               , only : photosyns_type
  use SoilHydrologyType               , only : soilhydrology_type  
  use SoilStateType                   , only : soilstate_type
  use SolarAbsorbedType               , only : solarabs_type
  use SurfaceRadiationMod             , only : surfrad_type
  use SurfaceAlbedoType               , only : surfalb_type
  use TemperatureType                 , only : temperature_type
  use WaterFluxType                   , only : waterflux_type
  use WaterStateType                  , only : waterstate_type
  use UrbanParamsType                 , only : urbanparams_type
  use HumanIndexMod                   , only : humanindex_type
  use VOCEmissionMod                  , only : vocemis_type
  use atm2lndType                     , only : atm2lnd_type
  use lnd2atmType                     , only : lnd2atm_type
  use lnd2glcMod                      , only : lnd2glc_type 
  use glc2lndMod                      , only : glc2lnd_type
  use glcDiagnosticsMod               , only : glc_diagnostics_type
  use GridcellType                    , only : grc                
  use LandunitType                    , only : lun                
  use ColumnType                      , only : col                
  use PatchType                       , only : patch                
  use EDTypesMod                      , only : ed_site_type
  use EDPhenologyType                 , only : ed_phenology_type
  use EDCLMLinkMod                    , only : ed_clm_type
  use SoilWaterRetentionCurveMod      , only : soil_water_retention_curve_type
  use NutrientCompetitionMethodMod    , only : nutrient_competition_method_type
  !
  use SoilStateInitTimeConstMod       , only : SoilStateInitTimeConst
  use SoilHydrologyInitTimeConstMod   , only : SoilHydrologyInitTimeConst
  use SurfaceAlbedoMod                , only : SurfaceAlbedoInitTimeConst 
  use LakeCon                         , only : LakeConInit 
  !
  implicit none
  public   ! By default everything is public 
  !
  !-----------------------------------------
  ! Instances of component types
  !-----------------------------------------

  ! Physics types 
  type(aerosol_type)                      :: aerosol_inst
  type(canopystate_type)                  :: canopystate_inst
  type(energyflux_type)                   :: energyflux_inst
  type(frictionvel_type)                  :: frictionvel_inst
  type(irrigation_type)                   :: irrigation_inst
  type(lakestate_type)                    :: lakestate_inst
  type(photosyns_type)                    :: photosyns_inst
  type(soilstate_type)                    :: soilstate_inst
  type(soilhydrology_type)                :: soilhydrology_inst
  type(solarabs_type)                     :: solarabs_inst
  type(surfalb_type)                      :: surfalb_inst
  type(surfrad_type)                      :: surfrad_inst
  type(temperature_type)                  :: temperature_inst
  type(urbanparams_type)                  :: urbanparams_inst
  type(humanindex_type)                   :: humanindex_inst
  type(waterflux_type)                    :: waterflux_inst
  type(waterstate_type)                   :: waterstate_inst
  type(atm2lnd_type)                      :: atm2lnd_inst
  type(glc2lnd_type)                      :: glc2lnd_inst
  type(lnd2atm_type)                      :: lnd2atm_inst
  type(lnd2glc_type)                      :: lnd2glc_inst
  type(glc_diagnostics_type)              :: glc_diagnostics_inst
  class(soil_water_retention_curve_type) , allocatable :: soil_water_retention_curve

  ! CN vegetation types  
  type(cnveg_state_type)                  :: cnveg_state_inst
  type(cnveg_carbonstate_type)            :: cnveg_carbonstate_inst
  type(cnveg_carbonstate_type)            :: c13_cnveg_carbonstate_inst
  type(cnveg_carbonstate_type)            :: c14_cnveg_carbonstate_inst
  type(cnveg_carbonflux_type)             :: cnveg_carbonflux_inst
  type(cnveg_carbonflux_type)             :: c13_cnveg_carbonflux_inst
  type(cnveg_carbonflux_type)             :: c14_cnveg_carbonflux_inst
  type(cnveg_nitrogenstate_type)          :: cnveg_nitrogenstate_inst
  type(cnveg_nitrogenflux_type)           :: cnveg_nitrogenflux_inst
  class(nutrient_competition_method_type), allocatable :: nutrient_competition_method

  ! Soil biogeochem types 
  type(soilbiogeochem_state_type)         :: soilbiogeochem_state_inst
  type(soilbiogeochem_carbonstate_type)   :: soilbiogeochem_carbonstate_inst
  type(soilbiogeochem_carbonstate_type)   :: c13_soilbiogeochem_carbonstate_inst
  type(soilbiogeochem_carbonstate_type)   :: c14_soilbiogeochem_carbonstate_inst
  type(soilbiogeochem_carbonflux_type)    :: soilbiogeochem_carbonflux_inst
  type(soilbiogeochem_carbonflux_type)    :: c13_soilbiogeochem_carbonflux_inst
  type(soilbiogeochem_carbonflux_type)    :: c14_soilbiogeochem_carbonflux_inst
  type(soilbiogeochem_nitrogenstate_type) :: soilbiogeochem_nitrogenstate_inst
  type(soilbiogeochem_nitrogenflux_type)  :: soilbiogeochem_nitrogenflux_inst

  ! General biogeochem types
  type(ch4_type)                          :: ch4_inst
  type(dgvs_type)                         :: dgvs_inst
  type(crop_type)                         :: crop_inst
  type(dust_type)                         :: dust_inst
  type(vocemis_type)                      :: vocemis_inst
  type(drydepvel_type)                    :: drydepvel_inst

  ! ED types passed in from top level
  type(ed_site_type), allocatable, target :: ed_allsites_inst(:)
  type(ed_phenology_type)                 :: ed_phenology_inst
  type(ed_clm_type)                       :: ed_clm_inst
  !
  public :: clm_instInit
  public :: clm_instRest
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine clm_instInit(bounds)
    !
    ! !USES: 
    use clm_varpar                         , only : nlevsno, numpft, crop_prog
    use controlMod                         , only : nlfilename, fsurdat
    use domainMod                          , only : ldomain
    use SoilBiogeochemDecompCascadeBGCMod  , only : init_decompcascade_bgc
    use SoilBiogeochemDecompCascadeCNMod   , only : init_decompcascade_cn
    use SoilBiogeochemDecompCascadeContype , only : init_decomp_cascade_constants
    use EDEcophysConType                   , only : EDecophysconInit 
    use EDPftVarcon                        , only : EDpftvarcon_inst
    use initVerticalMod                    , only : initVertical
    use accumulMod                         , only : print_accum_fields 
    use SoilWaterRetentionCurveFactoryMod  , only : create_soil_water_retention_curve
    !
    ! !ARGUMENTS    
    type(bounds_type), intent(in) :: bounds  ! processor bounds
    !
    ! !LOCAL VARIABLES:
    integer               :: c,l,g
    integer               :: begp, endp
    integer               :: begc, endc
    integer               :: begl, endl
    real(r8), allocatable :: h2osno_col(:)
    real(r8), allocatable :: snow_depth_col(:)
    !----------------------------------------------------------------------

    ! Note: h2osno_col and snow_depth_col are initialized as local variable 
    ! since they are needed to initialize vertical data structures  

    begp = bounds%begp; endp = bounds%endp 
    begc = bounds%begc; endc = bounds%endc 
    begl = bounds%begl; endl = bounds%endl 

    allocate (h2osno_col(begc:endc))
    allocate (snow_depth_col(begc:endc))

    ! snow water
    ! Note: Glacier_mec columns are initialized with half the maximum snow cover.
    ! This gives more realistic values of qflx_glcice sooner in the simulation
    ! for columns with net ablation, at the cost of delaying ice formation
    ! in columns with net accumulation.
    do c = begc,endc
       l = col%landunit(c)
       g = col%gridcell(c)

       if (lun%itype(l)==istice) then
          h2osno_col(c) = h2osno_max
       elseif (lun%itype(l)==istice_mec .or. &
              (lun%itype(l)==istsoil .and. abs(grc%latdeg(g)) >= 60._r8)) then 
          ! In order to speed equilibration of the snow pack, initialize non-zero snow
          ! thickness in some places. This is mainly of interest for glacier spinup.
          ! However, putting in an explicit dependence on glcmask is problematic, because
          ! that means that answers change simply due to changing glcmask (which may be
          ! done simply to have additional virtual columns for the sake of diagnostics).
          ! Thus, we apply this non-zero initialization at all high latitude soil points.
          h2osno_col(c) = 0.5_r8 * h2osno_max   ! 50 cm if h2osno_max = 1 m
       else
          h2osno_col(c) = 0._r8
       endif
       snow_depth_col(c)  = h2osno_col(c) / bdsno
    end do

    ! Initialize urban constants

    call urbanparams_inst%Init(bounds)
    call humanindex_inst%Init(bounds)

    ! Initialize vertical data components 

    call initVertical(bounds,               &
         snow_depth_col(begc:endc),              &
         urbanparams_inst%thick_wall(begl:endl), &
         urbanparams_inst%thick_roof(begl:endl))

    ! Initialize clm->drv and drv->clm data structures

    call atm2lnd_inst%Init( bounds )
    call lnd2atm_inst%Init( bounds )

    ! Initialize glc2lnd and lnd2glc even if running without create_glacier_mec_landunit,
    ! because at least some variables (such as the icemask) are referred to in code that
    ! is executed even when running without glc_mec.

    call glc2lnd_inst%Init( bounds )
    call lnd2glc_inst%Init( bounds )

    ! Initialization of public data types

    call temperature_inst%Init(bounds,           &
         urbanparams_inst%em_roof(begl:endl),    &
         urbanparams_inst%em_wall(begl:endl),    &
         urbanparams_inst%em_improad(begl:endl), &
         urbanparams_inst%em_perroad(begl:endl), &
         IsSimpleBuildTemp(), IsProgBuildTemp() )

    call canopystate_inst%Init(bounds)

    call soilstate_inst%Init(bounds)
    call SoilStateInitTimeConst(bounds, soilstate_inst, nlfilename) ! sets hydraulic and thermal soil properties

    call waterstate_inst%Init(bounds,         &
         h2osno_col(begc:endc),                    &
         snow_depth_col(begc:endc),                &
         soilstate_inst%watsat_col(begc:endc, 1:), &
         temperature_inst%t_soisno_col(begc:endc, -nlevsno+1:) )

    call waterflux_inst%Init(bounds)

    ! WJS (6-24-14): Without the following write statement, the assertion in
    ! energyflux_inst%Init fails with pgi 13.9 on yellowstone. So for now, I'm leaving
    ! this write statement in place as a workaround for this problem.
    call energyflux_inst%Init(bounds, temperature_inst%t_grnd_col(begc:endc), &
         IsSimpleBuildTemp(), IsProgBuildTemp() )

    call aerosol_inst%Init(bounds)

    call frictionvel_inst%Init(bounds)

    call lakestate_inst%Init(bounds)
    call LakeConInit()

    call photosyns_inst%Init(bounds)

    call soilhydrology_inst%Init(bounds, nlfilename)
    call SoilHydrologyInitTimeConst(bounds, soilhydrology_inst) ! sets time constant properties

    call solarabs_inst%Init(bounds)

    call surfalb_inst%Init(bounds)
    call SurfaceAlbedoInitTimeConst(bounds)

    call surfrad_inst%Init(bounds)

    call dust_inst%Init(bounds)

    call glc_diagnostics_inst%Init(bounds)

    ! Once namelist options are added to control the soil water retention curve method,
    ! we'll need to either pass the namelist file as an argument to this routine, or pass
    ! the namelist value itself (if the namelist is read elsewhere).

    allocate(soil_water_retention_curve, &
         source=create_soil_water_retention_curve())

    call irrigation_inst%init(bounds, soilstate_inst, soil_water_retention_curve)

    ! Note - always initialize the memory for ch4_inst
    call ch4_inst%Init(bounds, soilstate_inst%cellorg_col(begc:endc, 1:), fsurdat)

    ! Note - always initialize the memory for cnveg_state_inst (used in biogeophys/)
    call cnveg_state_inst%Init(bounds)

    if (use_voc ) then
       call vocemis_inst%Init(bounds)
    end if

    call drydepvel_inst%Init(bounds)

    if (use_cn) then

       ! Note - always initialize the memory for the c13_xxx_inst and
       ! c14_xxx_inst data structure so that they can be used in 
       ! associate statements (nag compiler complains otherwise)

       ! Note that SoillBiogeochem types must ALWAYS be allocated first - since CNVeg_xxxType
       ! can reference SoilBiogeochem types (for both carbon and nitrogen)

       ! Initialize soilbiogeochem_state_inst

       call soilbiogeochem_state_inst%Init(bounds)

       ! Initialize decompcascade constants
       ! Note that init_decompcascade_bgc and init_decompcascade_cn need 
       ! soilbiogeochem_state_inst to be initialized

       call init_decomp_cascade_constants()
       if (use_century_decomp) then
          call init_decompcascade_bgc(bounds, soilbiogeochem_state_inst, soilstate_inst)
       else 
          call init_decompcascade_cn(bounds, soilbiogeochem_state_inst)
       end if

       ! Initalize soilbiogeochem carbon and nitrogen types

       call soilbiogeochem_carbonstate_inst%Init(bounds, carbon_type='c12', ratio=1._r8)
       if (use_c13) then
          call c13_soilbiogeochem_carbonstate_inst%Init(bounds, carbon_type='c13', ratio=c13ratio, &
               c12_soilbiogeochem_carbonstate_inst=soilbiogeochem_carbonstate_inst)
       end if
       if (use_c14) then
          call c14_soilbiogeochem_carbonstate_inst%Init(bounds, carbon_type='c14', ratio=c14ratio, &
               c12_soilbiogeochem_carbonstate_inst=soilbiogeochem_carbonstate_inst)
       end if
       call soilbiogeochem_nitrogenstate_inst%Init(bounds, &
            soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools), &
            soilbiogeochem_carbonstate_inst%decomp_cpools_col(begc:endc,1:ndecomp_pools),  &
            soilbiogeochem_carbonstate_inst%decomp_cpools_1m_col(begc:endc, 1:ndecomp_pools))

       call soilbiogeochem_carbonflux_inst%Init(bounds, carbon_type='c12') 
       if (use_c13) then
          call c13_soilbiogeochem_carbonflux_inst%Init(bounds, carbon_type='c13')
       end if
       if (use_c14) then
          call c14_soilbiogeochem_carbonflux_inst%Init(bounds, carbon_type='c14')
       end if
       call soilbiogeochem_nitrogenflux_inst%Init(bounds) 

       ! Initalize cnveg carbon and nitrogen types

       call cnveg_carbonstate_inst%Init(bounds, carbon_type='c12', ratio=1._r8)
       if (use_c13) then
          call c13_cnveg_carbonstate_inst%Init(bounds, carbon_type='c13', ratio=c13ratio, &
               c12_cnveg_carbonstate_inst=cnveg_carbonstate_inst)
       end if
       if (use_c14) then
          call c14_cnveg_carbonstate_inst%Init(bounds, carbon_type='c14', ratio=c14ratio, &
               c12_cnveg_carbonstate_inst=cnveg_carbonstate_inst)
       end if
       call cnveg_carbonflux_inst%Init(bounds, carbon_type='c12')
       if (use_c13) then
          call c13_cnveg_carbonflux_inst%Init(bounds, carbon_type='c13')
       end if
       if (use_c14) then
          call c14_cnveg_carbonflux_inst%Init(bounds, carbon_type='c14')
       end if
       call cnveg_nitrogenstate_inst%Init(bounds,                  &
            cnveg_carbonstate_inst%leafc_patch(begp:endp),         &
            cnveg_carbonstate_inst%leafc_storage_patch(begp:endp), &
            cnveg_carbonstate_inst%deadstemc_patch(begp:endp))
       call cnveg_nitrogenflux_inst%Init(bounds) 

       ! Note - always initialize the memory for the dgvs_inst data structure so
       ! that it can be used in associate statements (nag compiler complains otherwise)

       call dgvs_inst%Init(bounds)

       call crop_inst%Init(bounds)
       
    end if ! end of if use_cn 

    ! NOTE (MV, 10-24-2014): because ed_allsites is currently passed as arguments to
    ! biogeophys routines in the present implementation - it needs to be allocated - 
    ! if use_ed is not set, then this will not contain any significant memory 
    ! if use_ed is true, then the actual memory for all of the ED data structures
    ! is allocated in the call to EDInitMod - called from clm_initialize

    allocate (ed_allsites_inst(bounds%begg:bounds%endg))
    if (use_ed) then
       call ed_clm_inst%Init(bounds)
       call ed_phenology_inst%Init(bounds)
       call EDecophysconInit( EDpftvarcon_inst, numpft)
    end if

    deallocate (h2osno_col)
    deallocate (snow_depth_col)

    ! ------------------------------------------------------------------------
    ! Initialize accumulated fields
    ! ------------------------------------------------------------------------

    ! The time manager needs to be initialized before this called is made, since
    ! the step size is needed. 

    call t_startf('init_accflds')

    call atm2lnd_inst%InitAccBuffer(bounds)

    call temperature_inst%InitAccBuffer(bounds)

    if (use_ed) then
       call ed_phenology_inst%initAccBuffer(bounds)
    endif

    call canopystate_inst%InitAccBuffer(bounds)

    if (use_cndv) then
       call dgvs_inst%InitAccBuffer(bounds)
    end if

    if (crop_prog) then
       call crop_inst%InitAccBuffer(bounds)
    end if

    call print_accum_fields()

    call t_stopf('init_accflds')

  end subroutine clm_instInit

  !-----------------------------------------------------------------------
  subroutine clm_instRest(bounds, ncid, flag)
    !
    ! !USES:
    use ncdio_pio       , only : file_desc_t
    use EDRestVectorMod , only : EDRest
    use UrbanParamsType , only : IsSimpleBuildTemp, IsProgBuildTemp
    !
    ! !DESCRIPTION:
    ! Define/write/read CLM restart file.
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)    :: bounds          
    type(file_desc_t) , intent(inout) :: ncid ! netcdf id
    character(len=*)  , intent(in)    :: flag ! 'define', 'write', 'read' 
    !-----------------------------------------------------------------------

    call atm2lnd_inst%restart (bounds, ncid, flag=flag)

    call canopystate_inst%restart (bounds, ncid, flag=flag)

    call energyflux_inst%restart (bounds, ncid, flag=flag, &
         is_simple_buildtemp=IsSimpleBuildTemp(), is_prog_buildtemp=IsProgBuildTemp())

    call frictionvel_inst% restart (bounds, ncid, flag=flag)

    call lakestate_inst%restart (bounds, ncid, flag=flag)

    call photosyns_inst%restart (bounds, ncid, flag=flag)

    call soilhydrology_inst%restart (bounds, ncid, flag=flag)

    call solarabs_inst%restart (bounds, ncid, flag=flag)

    call temperature_inst%restart (bounds, ncid, flag=flag, &
         is_simple_buildtemp=IsSimpleBuildTemp(), is_prog_buildtemp=IsProgBuildTemp())

    call waterflux_inst%restart (bounds, ncid, flag=flag)

    call waterstate_inst%restart (bounds, ncid, flag=flag, &
         watsat_col=soilstate_inst%watsat_col(bounds%begc:bounds%endc,:)) 

    call irrigation_inst%restart (bounds, ncid, flag=flag)

    call aerosol_inst%restart (bounds, ncid,  flag=flag, &
         h2osoi_ice_col=waterstate_inst%h2osoi_ice_col(bounds%begc:bounds%endc,:), &
         h2osoi_liq_col=waterstate_inst%h2osoi_liq_col(bounds%begc:bounds%endc,:))

    call surfalb_inst%restart (bounds, ncid, flag=flag, &
         tlai_patch=canopystate_inst%tlai_patch(bounds%begp:bounds%endp), &
         tsai_patch=canopystate_inst%tsai_patch(bounds%begp:bounds%endp))

    if (use_lch4) then
       call ch4_inst%restart(bounds, ncid, flag=flag)
    end if

    if (use_cn) then

       call soilbiogeochem_carbonstate_inst%restart(bounds, ncid, flag=flag, carbon_type='c12')
       if (use_c13) then
          call c13_soilbiogeochem_carbonstate_inst%restart(bounds, ncid, flag=flag, carbon_type='c13', &
               c12_soilbiogeochem_carbonstate_inst=soilbiogeochem_carbonstate_inst)
       end if
       if (use_c14) then
          call c14_soilbiogeochem_carbonstate_inst%restart(bounds, ncid, flag=flag, carbon_type='c14', &
               c12_soilbiogeochem_carbonstate_inst=soilbiogeochem_carbonstate_inst)
       end if
       call soilbiogeochem_carbonflux_inst%restart(bounds, ncid, flag=flag)
       call soilbiogeochem_nitrogenstate_inst%restart(bounds, ncid, flag=flag)
       call soilbiogeochem_nitrogenflux_inst%restart(bounds, ncid, flag=flag)

       call cnveg_state_inst%restart(bounds, ncid, flag=flag)
       call cnveg_carbonstate_inst%restart(bounds, ncid, flag=flag, carbon_type='c12')
       if (use_c13) then
          call c13_cnveg_carbonstate_inst%restart(bounds, ncid, flag=flag, carbon_type='c13', &
               c12_cnveg_carbonstate_inst=cnveg_carbonstate_inst)
       end if
       if (use_c14) then
          call c14_cnveg_carbonstate_inst%restart(bounds, ncid, flag=flag, carbon_type='c14', &
               c12_cnveg_carbonstate_inst=cnveg_carbonstate_inst)
       end if
       call cnveg_carbonflux_inst%restart(bounds, ncid, flag=flag)
       call cnveg_nitrogenstate_inst%restart(bounds, ncid, flag=flag)
       call cnveg_nitrogenflux_inst%restart(bounds, ncid, flag=flag)

    end if

    if (use_cndv) then
       call dgvs_inst%Restart(bounds, ncid, flag=flag)
    end if

    if (use_ed) then
       call EDRest ( bounds, ncid, flag, ed_allsites_inst(bounds%begg:bounds%endg), &
            ed_clm_inst, ed_phenology_inst, waterstate_inst, canopystate_inst )
    end if

  end subroutine clm_instRest

end module clm_instMod

