module clm_instMod
  !-----------------------------------------------------------------------
  ! initialize clm data types
  !
  use shr_kind_mod               , only : r8 => shr_kind_r8
  use shr_log_mod                , only : errMsg => shr_log_errMsg
  use decompMod                  , only : bounds_type, get_proc_bounds
  use clm_varctl                 , only : use_cn, use_voc, use_c13, use_c14, use_fates, use_betr
  !-----------------------------------------
  ! Definition of component types
  !-----------------------------------------
  use AerosolType                , only : aerosol_type
  use CanopyStateType            , only : canopystate_type
  use CH4Mod                     , only : ch4_type
  use CNCarbonFluxType           , only : carbonflux_type
  use CNCarbonStateType          , only : carbonstate_type
  use CNStateType                , only : cnstate_type
  use CNNitrogenFluxType         , only : nitrogenflux_type
  use CNNitrogenStateType        , only : nitrogenstate_type

  use PhosphorusFluxType         , only : phosphorusflux_type
  use PhosphorusStateType        , only : phosphorusstate_type

  use CropType                   , only : crop_type
  use DryDepVelocity             , only : drydepvel_type
  use DUSTMod                    , only : dust_type
  use EnergyFluxType             , only : energyflux_type
  use FrictionVelocityType       , only : frictionvel_type
  use LakeStateType              , only : lakestate_type
  use PhotosynthesisType         , only : photosyns_type
  use SedFluxType                , only : sedflux_type
  use SoilHydrologyType          , only : soilhydrology_type
  use SoilStateType              , only : soilstate_type
  use SolarAbsorbedType          , only : solarabs_type
  use SurfaceRadiationMod        , only : surfrad_type
  use SurfaceAlbedoMod           , only : SurfaceAlbedoInitTimeConst !TODO - can this be merged into the type?
  use SurfaceAlbedoType          , only : surfalb_type
  use TemperatureType            , only : temperature_type
  use WaterfluxType              , only : waterflux_type
  use WaterstateType             , only : waterstate_type
  use VOCEmissionMod             , only : vocemis_type
  use atm2lndType                , only : atm2lnd_type
  use lnd2atmType                , only : lnd2atm_type
  use lnd2glcMod                 , only : lnd2glc_type
  use glc2lndMod                 , only : glc2lnd_type
  use glcDiagnosticsMod          , only : glc_diagnostics_type
  use SoilWaterRetentionCurveMod , only : soil_water_retention_curve_type
  use VegetationPropertiesType   , only : veg_vp             ! Ecophysical Constants
  use SoilorderConType           , only : soilordercon         ! Constants

  use GridcellDataType           , only : grc_es, grc_ef, grc_ws, grc_wf
  use GridcellDataType           , only : grc_cs, c13_grc_cs, c14_grc_cs
  use GridcellDataType           , only : grc_cf, c13_grc_cf, c14_grc_cf
  use GridcellDataType           , only : grc_ns, grc_nf
  use GridcellDataType           , only : grc_ps, grc_pf
  use LandunitType               , only : lun_pp
  use LandunitDataType           , only : lun_es, lun_ef, lun_ws
  use ColumnType                 , only : col_pp
  use ColumnDataType             , only : col_es, col_ef, col_ws, col_wf
  use ColumnDataType             , only : col_cs, c13_col_cs, c14_col_cs
  use ColumnDataType             , only : col_cf, c13_col_cf, c14_col_cf
  use ColumnDataType             , only : col_ns, col_nf
  use ColumnDataType             , only : col_ps, col_pf
  use VegetationType             , only : veg_pp
  use VegetationDataType         , only : veg_es, veg_ef, veg_ws, veg_wf
  use VegetationDataType         , only : veg_cs, c13_veg_cs, c14_veg_cs
  use VegetationDataType         , only : veg_cf, c13_veg_cf, c14_veg_cf
  use VegetationDataType         , only : veg_ns, veg_nf
  use VegetationDataType         , only : veg_ps, veg_pf

  use clm_interface_dataType     , only : clm_interface_data_type
  use ChemStateType              , only : chemstate_type     ! structure for chemical indices of the soil, such as pH and Eh
  use BeTRSimulationALM          , only : betr_simulation_alm_type
  use PlantMicKineticsMod        , only : PlantMicKinetics_type
  use CLMFatesInterfaceMod       , only : hlm_fates_interface_type

  ! instances declared in their own modules
  use UrbanParamsType            , only : urbanparams_vars


  !
  implicit none
  save

 public   ! By default everything is public
  !
  !-----------------------------------------
  ! Instances of component types
  !-----------------------------------------
  !
  type(ch4_type)                                      :: ch4_vars
  type(carbonstate_type)                              :: carbonstate_vars
  type(carbonstate_type)                              :: c13_carbonstate_vars
  type(carbonstate_type)                              :: c14_carbonstate_vars
  type(carbonflux_type)                               :: carbonflux_vars
  type(carbonflux_type)                               :: c13_carbonflux_vars
  type(carbonflux_type)                               :: c14_carbonflux_vars
  type(nitrogenstate_type)                            :: nitrogenstate_vars
  type(nitrogenflux_type)                             :: nitrogenflux_vars
  type(crop_type)                                     :: crop_vars
  type(cnstate_type)                                  :: cnstate_vars
  type(dust_type)                                     :: dust_vars
  type(vocemis_type)                                  :: vocemis_vars
  type(drydepvel_type)                                :: drydepvel_vars
  type(aerosol_type)                                  :: aerosol_vars
  type(canopystate_type)                              :: canopystate_vars
  type(energyflux_type)                               :: energyflux_vars
  type(frictionvel_type)                              :: frictionvel_vars
  type(lakestate_type)                                :: lakestate_vars
  type(photosyns_type)                                :: photosyns_vars
  type(sedflux_type)                                  :: sedflux_vars
  type(soilstate_type)                                :: soilstate_vars
  type(soilhydrology_type)                            :: soilhydrology_vars
  type(solarabs_type)                                 :: solarabs_vars
  type(surfalb_type)                                  :: surfalb_vars
  type(surfrad_type)                                  :: surfrad_vars
  type(temperature_type)                              :: temperature_vars
  type(waterflux_type)                                :: waterflux_vars
  type(waterstate_type)                               :: waterstate_vars
  type(atm2lnd_type)                                  :: atm2lnd_vars
  type(glc2lnd_type)                                  :: glc2lnd_vars
  type(lnd2atm_type)                                  :: lnd2atm_vars
  type(lnd2glc_type)                                  :: lnd2glc_vars
  type(glc_diagnostics_type)                          :: glc_diagnostics_vars
  class(soil_water_retention_curve_type), allocatable :: soil_water_retention_curve
  type(phosphorusstate_type)                          :: phosphorusstate_vars
  type(phosphorusflux_type)                           :: phosphorusflux_vars
  type(clm_interface_data_type)                       :: clm_interface_data
  type(chemstate_type)                                :: chemstate_vars
  type(hlm_fates_interface_type)                      :: alm_fates
  class(betr_simulation_alm_type), pointer            :: ep_betr
  type(PlantMicKinetics_type)                         :: PlantMicKinetics_vars
  public :: clm_inst_biogeochem
  public :: clm_inst_biogeophys
  public :: alm_fates

contains


  !-----------------------------------------------------------------------
  subroutine clm_inst_biogeochem(bounds_proc)

    !
    ! DESCRIPTION
    ! initialize biogeochemical variables
    use clm_varcon            , only : c13ratio, c14ratio
    use histFileMod           , only : hist_printflds
    implicit none
    type(bounds_type), intent(in) :: bounds_proc

    integer               :: begp, endp
    integer               :: begc, endc
    integer               :: begl, endl
    integer               :: begg, endg

    begp = bounds_proc%begp; endp = bounds_proc%endp
    begc = bounds_proc%begc; endc = bounds_proc%endc
    begl = bounds_proc%begl; endl = bounds_proc%endl
    begg = bounds_proc%begg; endg = bounds_proc%endg


    if (use_voc ) then
       call vocemis_vars%Init(bounds_proc)
    end if
    if (use_cn .or. use_fates) then

       ! Note - always initialize the memory for the c13_cs and
       ! c14_cs data structure so that they can be used in
       ! associate statements (nag compiler complains otherwise)

       call grc_cs%Init(begg, endg, carbon_type='c12')
       call col_cs%Init(begc, endc, carbon_type='c12', ratio=1._r8)
       call veg_cs%Init(begp, endp, carbon_type='c12', ratio=1._r8)
       if (use_c13) then
          call c13_grc_cs%Init(begg, endg,carbon_type='c13')
          call c13_col_cs%Init(begc, endc, carbon_type='c13', ratio=c13ratio, &
               c12_carbonstate_vars=col_cs)
          call c13_veg_cs%Init(begp, endp, carbon_type='c13', ratio=c13ratio)
       end if
       if (use_c14) then
          call c14_grc_cs%Init(begg, endg,carbon_type='c14')
          call c14_col_cs%Init(begc, endc, carbon_type='c14', ratio=c14ratio, &
               c12_carbonstate_vars=col_cs)
          call c14_veg_cs%Init(begp, endp, carbon_type='c14', ratio=c14ratio)
       end if

       ! Note - always initialize the memory for the c13_carbonflux_vars and
       ! c14_carbonflux_vars data structure so that they can be used in
       ! associate statements (nag compiler complains otherwise)

       call grc_cf%Init(begg, endg, carbon_type='c12')
       call col_cf%Init(begc, endc, carbon_type='c12')
       call veg_cf%Init(begp, endp, carbon_type='c12')
       if (use_c13) then
          call c13_grc_cf%Init(begg, endg, carbon_type='c13')
          call c13_col_cf%Init(begc, endc, carbon_type='c13')
          call c13_veg_cf%Init(begp, endp, carbon_type='c13')
       end if
       if (use_c14) then
          call c14_grc_cf%Init(begg, endg, carbon_type='c14')
          call c14_col_cf%Init(begc, endc, carbon_type='c14')
          call c14_veg_cf%Init(begp, endp, carbon_type='c14')
       end if
    endif

    if (use_cn) then
       call grc_ns%Init(begg, endg)
       call col_ns%Init(begc, endc, col_cs)
       call veg_ns%Init(begp, endp, veg_cs)
       
       call grc_nf%Init(begg, endg)
       call col_nf%Init(begc, endc)
       call veg_nf%Init(begp, endp)

       call grc_ps%Init(begg, endg)
       call col_ps%Init(begc, endc, col_cs)
       call veg_ps%Init(begp, endp, veg_cs)

       call grc_pf%Init(begg, endg)
       call col_pf%Init(begc, endc)
       call veg_pf%Init(begp, endp)

       call crop_vars%Init(bounds_proc)

       if(use_betr)then
         call PlantMicKinetics_vars%Init(bounds_proc)
       endif
    end if
    
    ! Initialize the Functionaly Assembled Terrestrial Ecosystem Simulator (FATES)
    if (use_fates) then
       call alm_fates%Init(bounds_proc)
    end if
       
    call hist_printflds()

  end subroutine clm_inst_biogeochem


  !-----------------------------------------------------------------------

    subroutine clm_inst_biogeophys(bounds_proc)
    !
    ! DESCRIPTION
    ! initialize biogeophysical variables
    !
    use shr_scam_mod                      , only : shr_scam_getCloseLatLon
    use landunit_varcon                   , only : istice, istice_mec, istsoil
    use clm_varcon                        , only : h2osno_max, bdsno
    use domainMod                         , only : ldomain
    use clm_varpar                        , only : nlevsno, numpft
    use clm_varctl                        , only : single_column, fsurdat, scmlat, scmlon
    use controlMod                        , only : nlfilename
    use SoilWaterRetentionCurveFactoryMod , only : create_soil_water_retention_curve
    use fileutils                         , only : getfil
    use VegetationPropertiesType          , only : veg_vp
    use SoilorderConType                  , only : soilorderconInit
    use LakeCon                           , only : LakeConInit
    use initVerticalMod                   , only : initVertical
    ! !ARGUMENTS
    implicit none
    type(bounds_type), intent(in) :: bounds_proc
    ! LOCAL VARIABLES
    integer               :: c,i,g,j,k,l,p! indices
    integer               :: begp, endp
    integer               :: begc, endc
    integer               :: begl, endl
    integer               :: closelatidx,closelonidx
    real(r8)              :: closelat,closelon
    real(r8), allocatable :: h2osno_col(:)
    real(r8), allocatable :: snow_depth_col(:)
    character(len=256)    :: locfn        ! local file name


    ! Note: h2osno_col and snow_depth_col are initialized as local variable
    ! since they are needed to initialize vertical data structures
    begp = bounds_proc%begp; endp = bounds_proc%endp
    begc = bounds_proc%begc; endc = bounds_proc%endc
    begl = bounds_proc%begl; endl = bounds_proc%endl

    allocate (h2osno_col(begc:endc))
    allocate (snow_depth_col(begc:endc))

    ! snow water
    ! Note: Glacier_mec columns are initialized with half the maximum snow cover.
    ! This gives more realistic values of qflx_glcice sooner in the simulation
    ! for columns with net ablation, at the cost of delaying ice formation
    ! in columns with net accumulation.
    do c = begc,endc
       l = col_pp%landunit(c)
       g = col_pp%gridcell(c)

       if (lun_pp%itype(l)==istice) then
          h2osno_col(c) = h2osno_max
       elseif (lun_pp%itype(l)==istice_mec .or. &
              (lun_pp%itype(l)==istsoil .and. ldomain%glcmask(g) > 0._r8)) then
          ! Initialize a non-zero snow thickness where the ice sheet can/potentially operate.
          ! Using glcmask to capture all potential vegetated points around GrIS (ideally
          ! we would use icemask from CISM, but that isn't available until after initialization.)
          h2osno_col(c) = 1.0_r8 * h2osno_max   ! start with full snow column so +SMB can begin immediately
       else
          h2osno_col(c) = 0._r8
       endif
       snow_depth_col(c)  = h2osno_col(c) / bdsno
    end do

   ! Initialize urban constants

    call urbanparams_vars%Init(bounds_proc)

    ! Initialize ecophys constants

    call veg_vp%Init()

    ! Initialize soil order related constants

    call soilorderconInit()

    ! Initialize lake constants

    call LakeConInit()

    ! Initialize surface albedo constants

    call SurfaceAlbedoInitTimeConst(bounds_proc)

    ! Initialize vertical data components

    call initVertical(bounds_proc,               &
         snow_depth_col(begc:endc),              &
         urbanparams_vars%thick_wall(begl:endl), &
         urbanparams_vars%thick_roof(begl:endl))

    ! Initialize clm->drv and drv->clm data structures

    call atm2lnd_vars%Init( bounds_proc )
    call lnd2atm_vars%Init( bounds_proc )

    ! Initialize glc2lnd and lnd2glc even if running without create_glacier_mec_landunit,
    ! because at least some variables (such as the icemask) are referred to in code that
    ! is executed even when running without glc_mec.
    call glc2lnd_vars%Init( bounds_proc )
    call lnd2glc_vars%Init( bounds_proc )

    ! If single-column determine closest latitude and longitude

    if (single_column) then
       call getfil (fsurdat, locfn, 0)
       call shr_scam_getCloseLatLon(locfn, scmlat, scmlon, &
            closelat, closelon, closelatidx, closelonidx)
    end if

    ! Initialization of public data types

    call grc_es%Init(bounds_proc%begg_all, bounds_proc%endg_all)
    call lun_es%Init(bounds_proc%begl_all, bounds_proc%endl_all)
    call col_es%Init(bounds_proc%begc_all, bounds_proc%endc_all)
    call veg_es%Init(bounds_proc%begp_all, bounds_proc%endp_all)
    
    call canopystate_vars%init(bounds_proc)

    call soilstate_vars%init(bounds_proc)

    call waterstate_vars%init(bounds_proc,         &
         h2osno_col(begc:endc),                    &
         snow_depth_col(begc:endc),                &
         soilstate_vars%watsat_col(begc:endc, 1:), &
         col_es%t_soisno(begc:endc, -nlevsno+1:) )

    call grc_ws%Init(bounds_proc%begg_all, bounds_proc%endg_all)
    call lun_ws%Init(bounds_proc%begl_all, bounds_proc%endl_all)
    call col_ws%Init(bounds_proc%begc_all, bounds_proc%endc_all, &
         h2osno_col(begc:endc),                    &
         snow_depth_col(begc:endc),                &
         soilstate_vars%watsat_col(begc:endc, 1:))
    call veg_ws%Init(bounds_proc%begp_all, bounds_proc%endp_all)

    call waterflux_vars%init(bounds_proc)

    call grc_wf%Init(bounds_proc%begg_all, bounds_proc%endg_all, bounds_proc)
    call col_wf%Init(bounds_proc%begc_all, bounds_proc%endc_all)
    call veg_wf%Init(bounds_proc%begp_all, bounds_proc%endp_all)

    call chemstate_vars%Init(bounds_proc)
    ! WJS (6-24-14): Without the following write statement, the assertion in
    ! energyflux_vars%init fails with pgi 13.9 on yellowstone. So for now, I'm leaving
    ! this write statement in place as a workaround for this problem.
    call energyflux_vars%init(bounds_proc, col_es%t_grnd(begc:endc))

    call grc_ef%Init(bounds_proc%begg_all, bounds_proc%endg_all)
    call lun_ef%Init(bounds_proc%begl_all, bounds_proc%endl_all)
    call col_ef%Init(bounds_proc%begc_all, bounds_proc%endc_all)
    call veg_ef%Init(bounds_proc%begp_all, bounds_proc%endp_all)

    call drydepvel_vars%Init(bounds_proc)
    call aerosol_vars%Init(bounds_proc)

    call frictionvel_vars%Init(bounds_proc)

    call lakestate_vars%Init(bounds_proc)

    call photosyns_vars%Init(bounds_proc)

    call soilhydrology_vars%Init(bounds_proc, nlfilename)

    call solarabs_vars%Init(bounds_proc)

    call surfalb_vars%Init(bounds_proc)

    call surfrad_vars%Init(bounds_proc)

    call dust_vars%Init(bounds_proc)

    call glc_diagnostics_vars%Init(bounds_proc)

    ! Once namelist options are added to control the soil water retention curve method,
    ! we'll need to either pass the namelist file as an argument to this routine, or pass
    ! the namelist value itself (if the namelist is read elsewhere).
    allocate(soil_water_retention_curve, &
         source=create_soil_water_retention_curve())


    ! Note - always initialize the memory for ch4_vars
    call ch4_vars%Init(bounds_proc, soilstate_vars%cellorg_col(begc:endc, 1:))

    ! Note - always initialize the memory for cnstate_vars (used in biogeophys/)
    call cnstate_vars%Init(bounds_proc)

    call sedflux_vars%Init(bounds_proc)
    ! --------------------------------------------------------------
    ! Initialise the BeTR
    ! --------------------------------------------------------------

    deallocate (h2osno_col)
    deallocate (snow_depth_col)

    end subroutine clm_inst_biogeophys


end module clm_instMod

