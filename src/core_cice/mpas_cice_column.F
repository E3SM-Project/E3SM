#define COMMA ,
#define STREAM_DEBUG_WRITE(M) ! write(stderrUnit,*) M
#define STREAM_WARNING_WRITE(M) write(stderrUnit,*) 'WARNING: '//M
#define STREAM_ERROR_WRITE(M) write(stderrUnit,*) 'ERROR: '//M

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_column
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 12th January 2015
!> \details
!>  
!
!-----------------------------------------------------------------------

module cice_column

  use mpas_grid_types
  use mpas_timekeeping

  implicit none

  private
  save

  ! constants
  real(kind=RKIND), parameter :: puny = 1e-11_RKIND
       
  ! tracer object
  type, private :: ciceTracerObjectType

     ! length of tracer array
     integer :: nTracers

     ! number of base tracers
     integer, parameter :: nBaseTracers = 3

     ! maximum number of ancestor tracers
     integer, parameter :: nMaxAncestorTracers = 2

     ! tracer array
     real(kind=RKIND), dimension(:), allocatable :: tracerArray ! trcrn

     ! index of the parent tracer
     integer, dimension(:), allocatable :: parentIndex ! trcr_depend

     ! first ancestor type mask
     real(kind=RKIND), dimension(:,:), allocatable :: firstAncestorMask !trcr_base

     ! indices of ancestor tracers excluding base tracer
     integer, dimension(:,:), allocatable :: ancestorIndices ! nt_strata

     ! number of ancestor tracers excluding base tracer
     integer, dimension(:), allocatable :: ancestorNumber ! n_trcr_strata

  end type ciceTracerObjectType

  type(ciceTracerObjectType), private :: ciceTracerObject


contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  column_vertical_thermodynamics
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 20th January 2015
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine column_vertical_thermodynamics(domain, clock)

    use ice_colpkg, only: colpkg_step_therm1

    type(domain_type), intent(inout) :: domain

    type(MPAS_clock_type), intent(in) :: clock

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         icestate, &
         tracers, &
         velocity_solver, &
         atmos_coupling, &
         atmos_forcing, &
         ocean_coupling, &
         ocean_forcing, &
         drag, &
         melt_growth_rates, &
         atmos_fluxes, &
         ocean_fluxes, &
         shortwave, &
         ponds, &
         aerosol, &
         diagnostics

    ! configs
    real(kind=RKIND), pointer :: &
         config_dt, &
         config_ratio_pond_depth_to_fraction, &
         config_min_retained_meltwater_fraction, &
         config_max_retained_meltwater_fraction, &
         config_pond_flushing_timescale

    character(len=strKIND), pointer :: &
         config_thermodynamics_type, &
         config_atmosphere_boundary_method, &
         config_pond_refreezing_type

    logical, pointer :: &
         config_use_ice_age, &
         config_use_first_year_ice, &
         config_use_aerosols, &
         config_use_meltponds, &
         config_use_cesm_meltponds, &
         config_use_level_meltponds, &
         config_use_topo_meltponds, &
         config_calc_surface_temperature, &
         config_calc_surface_stresses, &
         config_use_form_drag

    ! dimensions
    integer, pointer :: &
         nCellsSolve, &
         nCategories, &
         nIceLayers, &
         nSnowLayers, &
         nAerosols

    ! variables
    real(kind=RKIND), dimension(:), pointer :: &
         latCell, &
         iceAreaCellInitial, &
         iceAreaCell, &
         iceVolumeCell, &
         snowVolumeCell, &
         uVelocity, &
         vvelocity, &
         uAirVelocity, &
         vAirVelocity, &
         windSpeed, &
         airLevelHeight, &
         airSpecificHumidity, &
         airDensity, &
         airTemperature, &
         airReferenceTemperature2m, &
         airReferenceHumidity2m, &
         airReferenceSpeed10m, &
         airOceanDragCoefficientRatio, &
         oceanDragCoefficient, &
         oceanDragCoefficientSkin, &
         oceanDragCoefficientFloe, &
         oceanDragCoefficientKeel, &
         airDragCoefficient, &
         airDragCoefficientSkin, &
         airDragCoefficientFloe, &
         airDragCoefficientPond, &
         airDragCoefficientRidge, &
         dragFreeboard, &
         dragIceSnowDraft, &
         dragRidgeHeight, &
         dragRidgeSeparation, &
         dragKeelDepth, &
         dragKeelSeparation, &
         dragFloeLength, &
         dragFloeSeparation, &
         airStressForcingU, &
         airStressForcingV, &
         airStressCellU, &
         airStressCellV, &
         airPotentialTemperature, &
         seaSurfaceTemperature, &
         seaSurfaceSalinity, &
         seaFreezingTemperature, &
         oceanStressCellU, &
         oceanStressCellV, &
         freezingMeltingPotential, &
         lateralIceMeltFraction, &
         snowfallRate, &
         rainfallRate, &
         pondFreshWaterFlux, &
         surfaceHeatFlux, &
         surfaceConductiveFlux, &
         absorbedShortwaveFlux, &
         longwaveUp, &
         longwaveDown, &
         solarZenithAngleCosine, & 
         sensibleHeatFlux, &
         latentHeatFlux, &
         evaporativeHeatFlux, &
         oceanFreshWaterFlux, &
         oceanSaltFlux, &
         oceanHeatFlux, &
         oceanShortwaveFlux, &
         surfaceIceMelt, &
         basalIceMelt, &
         lateralIceMelt, &
         snowMelt, &
         congelation, &
         snowiceFormation, &
         frazilFormation, &
         meltOnset, &
         freezeOnset

    real(kind=RKIND), dimension(:,:), pointer :: &
         iceAreaCategoryInitial, &
         iceVolumeCategoryInitial, &
         snowVolumeCategoryInitial, &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory, &
         surfaceTemperature, &
         levelIceArea, &
         levelIceVolume, &
         pondArea, &
         pondDepth, &
         pondLidThickness, &
         iceAge, &
         firstYearIceArea, &
         surfaceShortwaveFlux, &
         interiorShortwaveFlux, &
         penetratingShortwaveFlux, &
         latentHeatFluxCategory, &
         surfaceIceMeltCategory, &
         basalIceMeltCategory, &
         snowMeltCategory, &
         congelationCategory, &
         snowiceFormationCategory, &
         atmosAerosolFlux, &
         oceanAerosolFlux, &
         pondSnowDepthDifference, &
         pondLidMeltFluxFraction, &
         surfaceHeatFluxCategory, &
         surfaceConductiveFluxCategory, &
         latentHeatFluxCouple, &
         surfaceHeatFluxCouple, &
         surfaceConductiveFluxCouple, &
         snowThicknessChange

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         snowEnthalpy, &
         iceEnthalpy, &
         iceSalinity, &
         absorbedShortwaveIceLayer, &
         absorbedShortwaveSnowLayer, &
         snowScatteringAerosol, &
         snowBodyAerosol, &
         iceScatteringAerosol, &
         iceBodyAerosol

    ! local
    integer :: &
         iCell, &
         iCategory, &
         iAerosol

    real(kind=RKIND), dimension(:,:,:), allocatable :: &
         specificSnowAerosol, &
         specificIceAerosol

    logical :: &
         northernHemisphereMask, &
         abortFlag

    character(len=strKIND) :: &
         abortMessage

    real(kind=RKIND) :: &
         dayOfYear

    ! day of year
    call get_day_of_year(clock, dayOfYear)

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestate)
       call MPAS_pool_get_subpool(block % structs, "icestate", tracers)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocity_solver)
       call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmos_coupling)
       call MPAS_pool_get_subpool(block % structs, "atmos_forcing", atmos_forcing)
       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", ocean_coupling)
       call MPAS_pool_get_subpool(block % structs, "ocean_forcing", ocean_forcing)
       call MPAS_pool_get_subpool(block % structs, "drag", drag)
       call MPAS_pool_get_subpool(block % structs, "melt_growth_rates", melt_growth_rates)
       call MPAS_pool_get_subpool(block % structs, "atmos_fluxes", atmos_fluxes)
       call MPAS_pool_get_subpool(block % structs, "ocean_fluxes", ocean_fluxes)
       call MPAS_pool_get_subpool(block % structs, "shortwave", shortwave)
       call MPAS_pool_get_subpool(block % structs, "ponds", ponds)
       call MPAS_pool_get_subpool(block % structs, "aerosol", aerosol)
       call MPAS_pool_get_subpool(block % structs, "diagnostics", diagnostics)

       call MPAS_pool_get_config(block % configs, "config_dt", config_dt)
       call MPAS_pool_get_config(block % configs, "config_thermodynamics_type", config_thermodynamics_type)
       call MPAS_pool_get_config(block % configs, "config_use_ice_age", config_use_ice_age)
       call MPAS_pool_get_config(block % configs, "config_use_first_year_ice", config_use_first_year_ice)
       call MPAS_pool_get_config(block % configs, "config_use_aerosols", config_use_aerosols)
       call MPAS_pool_get_config(block % configs, "config_use_meltponds", config_use_meltponds)
       call MPAS_pool_get_config(block % configs, "config_use_cesm_meltponds", config_use_cesm_meltponds)
       call MPAS_pool_get_config(block % configs, "config_use_level_meltponds", config_use_level_meltponds)
       call MPAS_pool_get_config(block % configs, "config_use_topo_meltponds", config_use_topo_meltponds)
       call MPAS_pool_get_config(block % configs, "config_atmosphere_boundary_method", config_atmosphere_boundary_method)
       call MPAS_pool_get_config(block % configs, "config_calc_surface_temperature", config_calc_surface_temperature)
       call MPAS_pool_get_config(block % configs, "config_calc_surface_stresses", config_calc_surface_stresses)
       call MPAS_pool_get_config(block % configs, "config_use_form_drag", config_use_form_drag)
       call MPAS_pool_get_config(block % configs, "config_ratio_pond_depth_to_fraction", config_ratio_pond_depth_to_fraction)
       call MPAS_pool_get_config(block % configs, "config_min_retained_meltwater_fraction", config_min_retained_meltwater_fraction)
       call MPAS_pool_get_config(block % configs, "config_max_retained_meltwater_fraction", config_max_retained_meltwater_fraction)
       call MPAS_pool_get_config(block % configs, "config_pond_flushing_timescale", config_pond_flushing_timescale)
       call MPAS_pool_get_config(block % configs, "config_pond_refreezing_type", config_pond_refreezing_type)

       call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)
       call MPAS_pool_get_dimension(mesh, "nIceLayers", nIceLayers)
       call MPAS_pool_get_dimension(mesh, "nSnowLayers", nSnowLayers)
       call MPAS_pool_get_dimension(mesh, "nAerosols", nAerosols)

       call MPAS_pool_get_array(mesh, "latCell", latCell)

       call MPAS_pool_get_array(icestate, "iceAreaCellInitial", iceAreaCellInitial)
       call MPAS_pool_get_array(icestate, "iceAreaCategoryInitial", iceAreaCategoryInitial)
       call MPAS_pool_get_array(icestate, "iceVolumeCategoryInitial", iceVolumeCategoryInitial)
       call MPAS_pool_get_array(icestate, "snowVolumeCategoryInitial", snowVolumeCategoryInitial)
       call MPAS_pool_get_array(icestate, "iceAreaCell", iceAreaCell)
       call MPAS_pool_get_array(icestate, "iceVolumeCell", iceVolumeCell)
       call MPAS_pool_get_array(icestate, "snowVolumeCell", snowVolumeCell)

       call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory)
       call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory)
       call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory)
       call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature)
       call MPAS_pool_get_array(tracers, "snowEnthalpy", snowEnthalpy)
       call MPAS_pool_get_array(tracers, "iceEnthalpy", iceEnthalpy)
       call MPAS_pool_get_array(tracers, "iceSalinity", iceSalinity)
       call MPAS_pool_get_array(tracers, "levelIceArea", levelIceArea)
       call MPAS_pool_get_array(tracers, "levelIceVolume", levelIceVolume)
       call MPAS_pool_get_array(tracers, "pondArea", pondArea)
       call MPAS_pool_get_array(tracers, "pondDepth", pondDepth)
       call MPAS_pool_get_array(tracers, "pondLidThickness", pondLidThickness)
       call MPAS_pool_get_array(tracers, "iceAge", iceAge)
       call MPAS_pool_get_array(tracers, "firstYearIceArea", firstYearIceArea)
       call MPAS_pool_get_array(tracers, "snowScatteringAerosol", snowScatteringAerosol)
       call MPAS_pool_get_array(tracers, "snowBodyAerosol", snowBodyAerosol)
       call MPAS_pool_get_array(tracers, "iceScatteringAerosol", iceScatteringAerosol)
       call MPAS_pool_get_array(tracers, "iceBodyAerosol", iceBodyAerosol)

       call MPAS_pool_get_array(velocity_solver, "uVelocity", uVelocity)
       call MPAS_pool_get_array(velocity_solver, "vVelocity", vVelocity)
       call MPAS_pool_get_array(velocity_solver, "airStressCellU", airStressCellU)
       call MPAS_pool_get_array(velocity_solver, "airStressCellV", airStressCellV)
       call MPAS_pool_get_array(velocity_solver, "oceanStressCellU", oceanStressCellU)
       call MPAS_pool_get_array(velocity_solver, "oceanStressCellV", oceanStressCellV)

       call MPAS_pool_get_array(atmos_coupling, "uAirVelocity", uAirVelocity)
       call MPAS_pool_get_array(atmos_coupling, "vAirVelocity", vAirVelocity)
       call MPAS_pool_get_array(atmos_coupling, "airLevelHeight", airLevelHeight)
       call MPAS_pool_get_array(atmos_coupling, "airSpecificHumidity", airSpecificHumidity)
       call MPAS_pool_get_array(atmos_coupling, "airDensity", airDensity)
       call MPAS_pool_get_array(atmos_coupling, "airTemperature", airTemperature)
       call MPAS_pool_get_array(atmos_coupling, "airPotentialTemperature", airPotentialTemperature)
       call MPAS_pool_get_array(atmos_coupling, "snowfallRate", snowfallRate)
       call MPAS_pool_get_array(atmos_coupling, "rainfallRate", rainfallRate)
       call MPAS_pool_get_array(atmos_coupling, "longwaveDown", longwaveDown)
       call MPAS_pool_get_array(atmos_coupling, "latentHeatFluxCouple", latentHeatFluxCouple)
       call MPAS_pool_get_array(atmos_coupling, "surfaceHeatFluxCouple", surfaceHeatFluxCouple)
       call MPAS_pool_get_array(atmos_coupling, "surfaceConductiveFluxCouple", surfaceConductiveFluxCouple)

       call MPAS_pool_get_array(atmos_forcing, "windSpeed", windSpeed)
       call MPAS_pool_get_array(atmos_forcing, "airReferenceTemperature2m", airReferenceTemperature2m)
       call MPAS_pool_get_array(atmos_forcing, "airReferenceHumidity2m", airReferenceHumidity2m)
       call MPAS_pool_get_array(atmos_forcing, "airReferenceSpeed10m", airReferenceSpeed10m)
       call MPAS_pool_get_array(atmos_forcing, "airStressForcingU", airStressForcingU)
       call MPAS_pool_get_array(atmos_forcing, "airStressForcingV", airStressForcingV)

       call MPAS_pool_get_array(ocean_coupling, "seaSurfaceTemperature", seaSurfaceTemperature)
       call MPAS_pool_get_array(ocean_coupling, "seaSurfaceSalinity", seaSurfaceSalinity)
       call MPAS_pool_get_array(ocean_coupling, "freezingMeltingPotential", freezingMeltingPotential)

       call MPAS_pool_get_array(ocean_forcing, "seaFreezingTemperature", seaFreezingTemperature)

       call MPAS_pool_get_array(drag, "airOceanDragCoefficientRatio", airOceanDragCoefficientRatio)
       call MPAS_pool_get_array(drag, "oceanDragCoefficient", oceanDragCoefficient)
       call MPAS_pool_get_array(drag, "oceanDragCoefficientSkin", oceanDragCoefficientSkin)
       call MPAS_pool_get_array(drag, "oceanDragCoefficientFloe", oceanDragCoefficientFloe)
       call MPAS_pool_get_array(drag, "oceanDragCoefficientKeel", oceanDragCoefficientKeel)
       call MPAS_pool_get_array(drag, "airDragCoefficient", airDragCoefficient)
       call MPAS_pool_get_array(drag, "airDragCoefficientSkin", airDragCoefficientSkin)
       call MPAS_pool_get_array(drag, "airDragCoefficientFloe", airDragCoefficientFloe)
       call MPAS_pool_get_array(drag, "airDragCoefficientPond", airDragCoefficientPond)
       call MPAS_pool_get_array(drag, "airDragCoefficientRidge", airDragCoefficientRidge)
       call MPAS_pool_get_array(drag, "dragFreeboard", dragFreeboard)
       call MPAS_pool_get_array(drag, "dragIceSnowDraft", dragIceSnowDraft)
       call MPAS_pool_get_array(drag, "dragRidgeHeight", dragRidgeHeight)
       call MPAS_pool_get_array(drag, "dragRidgeSeparation", dragRidgeSeparation)
       call MPAS_pool_get_array(drag, "dragKeelDepth", dragKeelDepth)
       call MPAS_pool_get_array(drag, "dragKeelSeparation", dragKeelSeparation)
       call MPAS_pool_get_array(drag, "dragFloeLength", dragFloeLength)
       call MPAS_pool_get_array(drag, "dragFloeSeparation", dragFloeSeparation)

       call MPAS_pool_get_array(melt_growth_rates, "lateralIceMeltFraction", lateralIceMeltFraction)
       call MPAS_pool_get_array(melt_growth_rates, "surfaceIceMelt", surfaceIceMelt)
       call MPAS_pool_get_array(melt_growth_rates, "surfaceIceMeltCategory", surfaceIceMeltCategory)
       call MPAS_pool_get_array(melt_growth_rates, "basalIceMelt", basalIceMelt )
       call MPAS_pool_get_array(melt_growth_rates, "basalIceMeltCategory", basalIceMeltCategory)
       call MPAS_pool_get_array(melt_growth_rates, "lateralIceMelt", lateralIceMelt)
       call MPAS_pool_get_array(melt_growth_rates, "snowMelt", snowMelt)
       call MPAS_pool_get_array(melt_growth_rates, "snowMeltCategory", snowMeltCategory)
       call MPAS_pool_get_array(melt_growth_rates, "congelation", congelation)
       call MPAS_pool_get_array(melt_growth_rates, "congelationCategory", congelationCategory)
       call MPAS_pool_get_array(melt_growth_rates, "snowiceFormation", snowiceFormation)
       call MPAS_pool_get_array(melt_growth_rates, "snowiceFormationCategory", snowiceFormationCategory)
       call MPAS_pool_get_array(melt_growth_rates, "snowThicknessChange", snowThicknessChange)
       call MPAS_pool_get_array(melt_growth_rates, "frazilFormation", frazilFormation)

       call MPAS_pool_get_array(atmos_fluxes, "surfaceHeatFlux", surfaceHeatFlux)
       call MPAS_pool_get_array(atmos_fluxes, "surfaceHeatFluxCategory", surfaceHeatFluxCategory)
       call MPAS_pool_get_array(atmos_fluxes, "surfaceConductiveFlux", surfaceConductiveFlux)
       call MPAS_pool_get_array(atmos_fluxes, "surfaceConductiveFluxCategory", surfaceConductiveFluxCategory)
       call MPAS_pool_get_array(atmos_fluxes, "longwaveUp", longwaveUp)
       call MPAS_pool_get_array(atmos_fluxes, "sensibleHeatFlux", sensibleHeatFlux)
       call MPAS_pool_get_array(atmos_fluxes, "latentHeatFlux", latentHeatFlux)
       call MPAS_pool_get_array(atmos_fluxes, "latentHeatFluxCategory", latentHeatFluxCategory)
       call MPAS_pool_get_array(atmos_fluxes, "evaporativeHeatFlux", evaporativeHeatFlux)

       call MPAS_pool_get_array(ocean_fluxes, "oceanFreshWaterFlux", oceanFreshWaterFlux)
       call MPAS_pool_get_array(ocean_fluxes, "oceanSaltFlux", oceanSaltFlux)
       call MPAS_pool_get_array(ocean_fluxes, "oceanHeatFlux", oceanHeatFlux)
       call MPAS_pool_get_array(ocean_fluxes, "oceanShortwaveFlux", oceanShortwaveFlux)

       call MPAS_pool_get_array(shortwave, "surfaceShortwaveFlux", surfaceShortwaveFlux)
       call MPAS_pool_get_array(shortwave, "interiorShortwaveFlux", interiorShortwaveFlux)
       call MPAS_pool_get_array(shortwave, "penetratingShortwaveFlux", penetratingShortwaveFlux)
       call MPAS_pool_get_array(shortwave, "absorbedShortwaveFlux", absorbedShortwaveFlux)
       call MPAS_pool_get_array(shortwave, "absorbedShortwaveIceLayer", absorbedShortwaveIceLayer)
       call MPAS_pool_get_array(shortwave, "absorbedShortwaveSnowLayer", absorbedShortwaveSnowLayer)
       call MPAS_pool_get_array(shortwave, "solarZenithAngleCosine", solarZenithAngleCosine)

       call MPAS_pool_get_array(aerosol, "atmosAerosolFlux", atmosAerosolFlux)
       call MPAS_pool_get_array(aerosol, "oceanAerosolFlux", oceanAerosolFlux)

       call MPAS_pool_get_array(ponds, "pondFreshWaterFlux", pondFreshWaterFlux)
       call MPAS_pool_get_array(ponds, "pondSnowDepthDifference", pondSnowDepthDifference)
       call MPAS_pool_get_array(ponds, "pondLidMeltFluxFraction", pondLidMeltFluxFraction)

       call MPAS_pool_get_array(diagnostics, "meltOnset", meltOnset)
       call MPAS_pool_get_array(diagnostics, "freezeOnset", freezeOnset)

       do iCell = 1, nCellsSolve

          ! initial state values
          iceAreaCellInitial(iCell) = iceAreaCell(iCell)

          do iCategory = 1, nCategories

             iceAreaCategoryInitial(iCategory,iCell) = iceAreaCategory(iCategory,iCell)
             iceVolumeCategoryInitial(iCategory,iCell) = iceVolumeCategory(iCategory,iCell)
             snowVolumeCategoryInitial(iCategory,iCell) = snowVolumeCategory(iCategory,iCell)

          enddo ! iCategory             

          ! aerosol
          if (config_use_aerosols) then

             allocate(specificSnowAerosol(nAerosols, 2, nCategories))
             allocate(specificSnowAerosol(nAerosols, 2, nCategories))

             do iCategory = 1, nCategories
                do iAerosol = 1, nAerosols
                   
                   specificSnowAerosol(iAerosol, 1, iCategory) = snowScatteringAerosol(iAerosol,iCategory,iCell) * snowVolumeCategoryInitial(iCategory,iCell)
                   specificSnowAerosol(iAerosol, 2, iCategory) = snowBodyAerosol(iAerosol,iCategory,iCell)       * snowVolumeCategoryInitial(iCategory,iCell)
                   
                   specificIceAerosol(iAerosol, 1, iCategory) = iceScatteringAerosol(iAerosol,iCategory,iCell)   * iceVolumeCategoryInitial(iCategory,iCell)
                   specificIceAerosol(iAerosol, 2, iCategory) = iceBodyAerosol(iAerosol,iCategory,iCell)         * iceVolumeCategoryInitial(iCategory,iCell)
                   
                enddo ! iAerosol
             enddo ! iCategory

          end if

          ! hemisphere mask
          if (latCell(iCell) > 0.0_RKIND) then
             northernHemisphereMask = .true.
          else
             northernHemisphereMask = .false.
          endif
          
          call colpkg_step_therm1(&
               config_dt, &
               nCategories, &
               nIceLayers, &
               nSnowLayers, &
               config_cice_int("config_thermodynamics_type", config_thermodynamics_type), &
               config_use_ice_age, &
               config_use_first_year_ice, &
               config_use_aerosols, &
               config_use_meltponds, &
               config_use_cesm_meltponds, &
               config_use_level_meltponds, &
               config_use_topo_meltponds, &
               config_atmosphere_boundary_method, &
               config_calc_surface_temperature, &
               config_calc_surface_stresses, &
               config_use_form_drag, &
               iceAreaCategoryInitial(:,iCell), &
               iceVolumeCategoryInitial(:,iCell), &
               snowVolumeCategoryInitial(:,iCell), &
               iceAreaCell(iCell), &
               iceAreaCategory(:,iCell), &
               iceVolumeCell(iCell), &
               iceVolumeCategory(:,iCell), &
               snowVolumeCell(iCell), &
               snowVolumeCategory(:,iCell), &
               uVelocity(iCell), &
               vVelocity(iCell), &
               surfaceTemperature(:,iCell), &
               snowEnthalpy(:,:,iCell), &
               iceEnthalpy(:,:,iCell), &
               iceSalinity(:,:,iCell), &
               levelIceArea(:,iCell), &
               levelIceVolume(:,iCell), &
               pondArea(:,iCell), &
               pondDepth(:,iCell), &
               pondLidThickness(:,iCell), &
               iceAge(:,iCell), &
               firstYearIceArea(:,iCell), &
               specificSnowAerosol(:,:,:), &
               specificIceAerosol(:,:,:), &
               uAirVelocity(iCell), &
               vAirVelocity(iCell), &
               windSpeed(iCell), &
               airLevelHeight(iCell), &
               airSpecificHumidity(iCell), &
               airDensity(iCell), &
               airTemperature(iCell), &
               airReferenceTemperature2m(iCell), &
               airReferenceHumidity2m(iCell), &
               airReferenceSpeed10m(iCell), &
               airOceanDragCoefficientRatio(iCell), &
               oceanDragCoefficient(iCell), &
               oceanDragCoefficientSkin(iCell), &
               oceanDragCoefficientFloe(iCell), &
               oceanDragCoefficientKeel(iCell), &
               airDragCoefficient(iCell), &
               airDragCoefficientSkin(iCell), &
               airDragCoefficientFloe(iCell), &
               airDragCoefficientPond(iCell), &
               airDragCoefficientRidge(iCell), &
               dragFreeboard(iCell), &
               dragIceSnowDraft(iCell), &
               dragRidgeHeight(iCell), &
               dragRidgeSeparation(iCell), &
               dragKeelDepth(iCell), &
               dragKeelSeparation(iCell), &
               dragFloeLength(iCell), &
               dragFloeSeparation(iCell), &
               airStressForcingU(iCell), &
               airStressForcingV(iCell), &
               airStressCellU(iCell), &
               airStressCellV(iCell), &
               airPotentialTemperature(iCell), &
               seaSurfaceTemperature(iCell), &
               seaSurfaceSalinity(iCell), &
               seaFreezingTemperature(iCell), &
               oceanStressCellU(iCell), &
               oceanStressCellV(iCell), &
               freezingMeltingPotential(iCell), &
               lateralIceMeltFraction(iCell), &
               snowfallRate(iCell), &
               rainfallRate(iCell), &
               pondFreshWaterFlux(iCell), &
               surfaceHeatFlux(iCell), &
               surfaceHeatFluxCategory(:,iCell), &
               surfaceConductiveFlux(iCell), &
               surfaceConductiveFluxCategory(:,iCell), &
               surfaceShortwaveFlux(:,iCell), &
               interiorShortwaveFlux(:,iCell), &
               penetratingShortwaveFlux(:,iCell), &
               absorbedShortwaveFlux(iCell), &
               longwaveUp(iCell), &
               absorbedShortwaveIceLayer(:,:,iCell), &
               absorbedShortwaveSnowLayer(:,:,iCell), &
               longwaveDown(iCell), &
               solarZenithAngleCosine(iCell), & 
               sensibleHeatFlux(iCell), &
               latentHeatFlux(iCell), &
               latentHeatFluxCategory(:,iCell), &
               evaporativeHeatFlux(iCell), &
               oceanFreshWaterFlux(iCell), &
               oceanSaltFlux(iCell), &
               oceanHeatFlux(iCell), &
               oceanShortwaveFlux(iCell), &
               latentHeatFluxCouple(:,iCell), &
               surfaceHeatFluxCouple(:,iCell), &
               surfaceConductiveFluxCouple(:,iCell), &
               atmosAerosolFlux(:,iCell), &
               oceanAerosolFlux(:,iCell), &
               config_ratio_pond_depth_to_fraction, &
               config_min_retained_meltwater_fraction, &
               config_max_retained_meltwater_fraction, &
               config_pond_flushing_timescale, &
               config_pond_refreezing_type, &
               pondSnowDepthDifference(:,iCell), &
               pondLidMeltFluxFraction(:,iCell), &
               surfaceIceMelt(iCell), &
               surfaceIceMeltCategory(:,iCell), &
               basalIceMelt(iCell), &
               basalIceMeltCategory(:,iCell), &
               lateralIceMelt(iCell), &
               snowMelt(iCell), &
               snowMeltCategory(:,iCell), &
               congelation(iCell), &
               congelationCategory(:,iCell), &
               snowiceFormation(iCell), &
               snowiceFormationCategory(:,iCell), &
               snowThicknessChange(:,iCell), &
               frazilFormation(iCell), &
               northernHemisphereMask, &
               .not. northernHemisphereMask, &
               meltOnset(iCell), &
               freezeOnset(iCell), &
               dayOfYear, &
               abortFlag, &
               abortMessage, &
               stdoutUnit)

          ! code abort
          if (abortFlag) then
             STREAM_ERROR_WRITE("column_vertical_thermodynamics")
             STREAM_ERROR_WRITE(trim(abortMessage))
             stop
          endif

          ! aerosol
          if (config_use_aerosols) then

             do iCategory = 1, nCategories
                do iAerosol = 1, nAerosols

                   if (snowVolumeCategory(iCategory,iCell) > puny) &
                        specificSnowAerosol(iAerosol, :, iCategory) = specificSnowAerosol(iAerosol, :, iCategory) / snowVolumeCategory(iCategory,iCell)

                   if (iceVolumeCategory(iCategory,iCell) > puny) &
                        specificIceAerosol(iAerosol, :, iCategory)  = specificIceAerosol(iAerosol, :, iCategory)  / iceVolumeCategory(iCategory,iCell)

                   snowScatteringAerosol(iAerosol,iCategory,iCell) = specificSnowAerosol(iAerosol, 1, iCategory)
                   snowBodyAerosol(iAerosol,iCategory,iCell)       = specificSnowAerosol(iAerosol, 2, iCategory)
                   
                   iceScatteringAerosol(iAerosol,iCategory,iCell)  = specificIceAerosol(iAerosol, 1, iCategory)
                   iceBodyAerosol(iAerosol,iCategory,iCell)        = specificIceAerosol(iAerosol, 2, iCategory)

                enddo ! iAerosol
             enddo ! iCategory

             deallocate(specificSnowAerosol)
             deallocate(specificIceAerosol)

          endif

       enddo ! iCell

       block => block % next
    end do
    
  end subroutine column_vertical_thermodynamics

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  column_itd_thermodynamics
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 21th January 2015
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine column_itd_thermodynamics(domain)

    use ice_colpkg, only: colpkg_step_therm2

    type(domain_type), intent(inout) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         icestate, &
         tracers, &
         atmos_coupling, &
         ocean_coupling, &
         ocean_forcing, &
         melt_growth_rates, &
         ponds, &
         biogeochemistry, &
         diagnostics

    ! configs
    real(kind=RKIND), pointer :: &
         config_dt

    logical, pointer :: &
         config_use_aerosols, &
         config_use_topo_meltponds, &
         config_heat_capacity, &
         config_update_ocean_fluxes

    ! dimensions
    integer, pointer :: &
         nCellsSolve, &
         nCategories, &
         nIceLayers, &
         nSnowLayers, &
         nAerosols, &
         nZBGCTracers

    ! variables
    real(kind=RKIND), dimension(:), pointer :: &
         openWaterArea, &
         iceAreaCell, &
         seaFreezingTemperature, &
         seaSurfaceSalinity, &
         lateralIceMeltFraction, &
         lateralIceMelt, &
         freezingMeltingPotential, &
         frazilFormation, &
         rainfallRate, &
         pondFreshWaterFlux, &
         oceanFreshWaterFlux, &
         oceanSaltFlux, &
         oceanHeatFlux, &
         freezeOnset

    real(kind=RKIND), dimension(:,:), pointer :: &
         iceAreaCategoryInitial, &
         iceVolumeCategoryInitial, &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory, &
         oceanAerosolFlux, &
         newlyFormedIce, &
         oceanBioFluxes, &
         oceanBioConcentrations

    ! local
    integer :: &
         iCell

    logical :: &
         abortFlag

    character(len=strKIND) :: &
         abortMessage

    real(kind=RKIND) :: &
         dayOfYear

    ! day of year
    call get_day_of_year(clock, dayOfYear)

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestate)
       call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
       call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmos_coupling)
       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", ocean_coupling)
       call MPAS_pool_get_subpool(block % structs, "ocean_forcing", ocean_forcing)
       call MPAS_pool_get_subpool(block % structs, "ocean_fluxes", ocean_fluxes)
       call MPAS_pool_get_subpool(block % structs, "melt_growth_rates", melt_growth_rates)
       call MPAS_pool_get_subpool(block % structs, "ponds", ponds)
       call MPAS_pool_get_subpool(block % structs, "biogeochemistry", biogeochemistry)
       call MPAS_pool_get_subpool(block % structs, "diagnostics", diagnostics)

       call MPAS_pool_get_config(block % configs, "config_dt", config_dt)
       call MPAS_pool_get_config(block % configs, "config_use_aerosols", config_use_aerosols)
       call MPAS_pool_get_config(block % configs, "config_use_topo_meltponds", config_use_topo_meltponds)
       call MPAS_pool_get_config(block % configs, "config_heat_capacity", config_heat_capacity)
       call MPAS_pool_get_config(block % configs, "config_update_ocean_fluxes", config_update_ocean_fluxes)

       call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)
       call MPAS_pool_get_dimension(mesh, "nIceLayers", nIceLayers)
       call MPAS_pool_get_dimension(mesh, "nSnowLayers", nSnowLayers)
       call MPAS_pool_get_dimension(mesh, "nAerosols", nAerosols)
       call MPAS_pool_get_dimension(mesh, "nZBGCTracers", nZBGCTracers)

       call MPAS_pool_get_array(icestate, "iceAreaCategoryInitial", iceAreaCategoryInitial)
       call MPAS_pool_get_array(icestate, "iceVolumeCategoryInitial", iceVolumeCategoryInitial)
       call MPAS_pool_get_array(icestate, "openWaterArea", openWaterArea)
       call MPAS_pool_get_array(icestate, "iceAreaCell", iceAreaCell)

       call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory)
       call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory)
       call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory)

       call MPAS_pool_get_array(atmos_coupling, "rainfallRate", rainfallRate)

       call MPAS_pool_get_array(ocean_coupling, "freezingMeltingPotential", freezingMeltingPotential)

       call MPAS_pool_get_array(ocean_forcing, "seaFreezingTemperature", seaFreezingTemperature)
       call MPAS_pool_get_array(ocean_forcing, "seaSurfaceSalinity", seaSurfaceSalinity)

       call MPAS_pool_get_array(ocean_fluxes, "oceanFreshWaterFlux", oceanFreshWaterFlux)
       call MPAS_pool_get_array(ocean_fluxes, "oceanSaltFlux", oceanSaltFlux)
       call MPAS_pool_get_array(ocean_fluxes, "oceanHeatFlux", oceanHeatFlux)

       call MPAS_pool_get_array(melt_growth_rates, "lateralIceMeltFraction", lateralIceMeltFraction)
       call MPAS_pool_get_array(melt_growth_rates, "lateralIceMelt", lateralIceMelt)
       call MPAS_pool_get_array(melt_growth_rates, "frazilFormation", frazilFormation)

       call MPAS_pool_get_array(ponds, "pondFreshWaterFlux", pondFreshWaterFlux)

       call MPAS_pool_get_array(aerosol, "oceanAerosolFlux", oceanAerosolFlux)

       call MPAS_pool_get_array(biogeochemistry, "newlyFormedIce", newlyFormedIce)
       call MPAS_pool_get_array(biogeochemistry, "oceanBioFluxes", oceanBioFluxes
       call MPAS_pool_get_array(biogeochemistry, "oceanBioConcentrations", oceanBioConcentrations)

       call MPAS_pool_get_array(diagnostics, "freezeOnset", freezeOnset)

       do iCell = 1, nCellsSolve

          call colpkg_step_therm2(&
               config_dt, &
               nCategories, &
               nAerosols, &
               ciceTracerObject % nTracers, & !ntrcr
               nIcelayers, &
               nSnowLayers, &
               config_cice_int("config_itd_conversion_type", config_itd_conversion_type)
               hin_max, & !!!!
               .true., &
               iceAreaCategory(:,iCell), &
               iceVolumeCategory(:,iCell) &
               snowVolumeCategory(:,iCell), &
               iceAreaCategoryInitial(:,iCell), &
               iceVolumeCategoryInitial(:,iCell), &
               ciceTracerObject % tracerArray, & !trcrn
               openWaterArea(iCell), &
               iceAreaCell(iCell), &
               ciceTracerObject % parentIndex, & !trcr_depend
               ciceTracerObject % firstAncestor, & !trcr_base
               ciceTracerObject % ancestorIndices, & !n_trcr_strata
               ciceTracerObject % ancestorNumber, & !nt_strata
               seaFreezingTemperature(iCell), &
               seaSurfaceSalinity(iCell), &
               salinz, & !!!!
               phi_init, & !!!!
               dSin0_frazil, & !!!!
               lateralIceMeltFraction(iCell), &
               lateralIceMelt(iCell), &
               freezingMeltingPotential(iCell), &
               frazilFormation(iCell), &
               rainfallRate(iCell), &
               pondFreshWaterFlux(iCell), &
               oceanFreshWaterFlux(iCell), &
               oceanSaltFlux(iCell), &
               oceanHeatFlux(iCell), &
               config_update_ocean_fluxes, &
               oceanAerosolFlux(:,iCell), &
               config_use_aerosols, &
               config_use_topo_meltponds, &
               config_heat_capacity, &
               nZBGCTracers, &
               newlyFormedIce(:,iCell), &
               oceanBioFluxes(:,iCell), &
               oceanBioConcentrations(:,iCell), &
               abortFlag, &
               abortMessage, &
               stdoutUnit, &
               freezeOnset(iCell), &
               dayOfYear)

          ! code abort
          if (abortFlag) then
             STREAM_ERROR_WRITE("column_itd_thermodynamics")
             STREAM_ERROR_WRITE(trim(abortMessage))
             stop
          endif

       enddo ! iCell

       block => block % next
    end do

  end subroutine column_itd_thermodynamics

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  config_cice_int
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 20th January 2015
!> \details
!>  
!
!-----------------------------------------------------------------------

  function config_cice_int(configName, configValue) result(configValueCice)

    character(len=*), intent(in) :: &
         configName, &
         configValue

    integer :: configValueCice

    select case (trim(configName))

    ! ktherm
    case ("config_thermodynamics_type")
       
       select case (trim(configValue))
       case ("zero layer")
          configValueCice = 0
       case ("bl99")
          configValueCice = 1
       case ("mushy")
          configValueCice = 2
       end select

    ! kitd
    case ("config_itd_conversion_type")

       select case (trim(configValue))
       case ("delta function")
          configValueCice = 0
       case ("linear remap")
          configValueCice = 1
       end select

    end select

  end function config_cice_int

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  get_day_of_year
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 20th January 2015
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine get_day_of_year(clock, dayOfYear)

    type(MPAS_clock_type), intent(in) :: &
         clock

    real(kind=RKIND), intent(out) :: &
         dayOfYear

    dayOfYear = 0.0_RKIND

  end subroutine get_day_of_year

!-----------------------------------------------------------------------
! CICE tracer object
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_cice_tracer_object
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 22nd January 2015
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine init_cice_tracer_object(domain, tracerObject)

    type(domain_type), intent(in) :: &
         domain

    type(ciceTracerObjectType), intent(inout) :: &
         tracerObject

    ! get the number of CICE tracers in trcrn
    call init_cice_tracer_object_tracer_number(domain, tracerObject)

    ! allocate the tracer array
    allocate(tracerObject % tracerArray(tracerObject % nTracers))

    ! allocate other arrays
    allocate(tracerObject % parentIndex(tracerObject % nTracers))
    allocate(tracerObject % firstAncestorMask(tracerObject % nTracers, tracerObject % nBaseTracers))
    allocate(tracerObject % ancestorIndices(tracerObject % nTracers, tracerObject % nMaxAncestorTracers))
    allocate(tracerObject % ancestorNumber(tracerObject % nTracers))

    ! set the parent indices
    call init_cice_tracer_object_parent_indices()

  end subroutine init_cice_tracer_object

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_cice_tracer_object_tracer_number
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 22nd January 2015
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine init_cice_tracer_object_tracer_number(domain, tracerObject)

    type(domain_type), intent(in) :: &
         domain

    type(ciceTracerObjectType), intent(inout) :: &
         tracerObject

    logical, pointer :: &
         config_use_ice_age, &
         config_use_first_year_ice, &
         config_use_level_ice, &
         config_use_meltponds, &
         config_use_level_meltponds, &
         config_use_topo_meltponds, &
         config_use_aerosols

    integer, pointer :: &
         nIceLayers, &
         nSnowLayers, &
         nAerosols

    call MPAS_pool_get_config(domain % configs, "config_use_ice_age", config_use_ice_age)
    call MPAS_pool_get_config(domain % configs, "config_use_first_year_ice", config_use_first_year_ice)
    call MPAS_pool_get_config(domain % configs, "config_use_level_ice", config_use_level_ice)
    call MPAS_pool_get_config(domain % configs, "config_use_meltponds", config_use_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_level_meltponds", config_use_level_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_topo_meltponds", config_use_topo_meltponds)
    call MPAS_pool_get_config(domain % configs, "config_use_aerosols", config_use_aerosols)

    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nIceLayers", nIceLayers)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nSnowLayers", nSnowLayers)
    call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nAerosols", nAerosols)

    ! surfaceTemperature
    tracerObject % nTracers = 1

    ! iceEnthalpy
    tracerObject % nTracers = tracerObject % nTracers + nIceLayers

    ! snowEnthalpy
    tracerObject % nTracers = tracerObject % nTracers + nSnowLayers

    ! ice Salinity
    tracerObject % nTracers = tracerObject % nTracers + nIceLayers

    ! iceAge
    if (config_use_ice_age) &
         tracerObject % nTracers = tracerObject % nTracers + 1

    ! firstYearIceArea
    if (config_use_first_year_ice) &
         tracerObject % nTracers = tracerObject % nTracers + 1

    ! level ice tracers
    if (config_use_level_ice) &
         tracerObject % nTracers = tracerObject % nTracers + 2
       
    ! pond tracers
    if (config_use_meltponds) &
         tracerObject % nTracers = tracerObject % nTracers + 2

    ! level or topo ponds
    if (config_use_meltponds .and. &
       (config_use_level_meltponds .or. config_use_topo_meltponds)) &
         tracerObject % nTracers = tracerObject % nTracers + 1

    ! aerosols
    if (config_use_aerosols) &
         tracerObject % nTracers = tracerObject % nTracers + nAerosols*4

  end subroutine init_cice_tracer_object_tracer_number

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_cice_tracer_object_parent_indices
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 22nd January 2015
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine init_cice_tracer_object_parent_indices




      trcr_depend(nt_Tsfc) = 0 ! ice/snow surface temperature
      do k = 1, nilyr
         trcr_depend(nt_sice + k - 1) = 1 ! volume-weighted ice salinity
         trcr_depend(nt_qice + k - 1) = 1 ! volume-weighted ice enthalpy
      enddo
      do k = 1, nslyr
         trcr_depend(nt_qsno + k - 1) = 2 ! volume-weighted snow enthalpy
      enddo
      if (tr_iage) trcr_depend(nt_iage)  = 1   ! volume-weighted ice age
      if (tr_FY)   trcr_depend(nt_FY)    = 0   ! area-weighted first-year ice area
      if (tr_lvl)  trcr_depend(nt_alvl)  = 0   ! level ice area
      if (tr_lvl)  trcr_depend(nt_vlvl)  = 1   ! level ice volume
      if (tr_pond_cesm) then
                   trcr_depend(nt_apnd)  = 0           ! melt pond area
                   trcr_depend(nt_hpnd)  = 2+nt_apnd   ! melt pond depth
      endif
      if (tr_pond_lvl) then
                   trcr_depend(nt_apnd)  = 2+nt_alvl   ! melt pond area
                   trcr_depend(nt_hpnd)  = 2+nt_apnd   ! melt pond depth
                   trcr_depend(nt_ipnd)  = 2+nt_apnd   ! refrozen pond lid
      endif
      if (tr_pond_topo) then
                   trcr_depend(nt_apnd)  = 0           ! melt pond area
                   trcr_depend(nt_hpnd)  = 2+nt_apnd   ! melt pond depth
                   trcr_depend(nt_ipnd)  = 2+nt_apnd   ! refrozen pond lid
      endif
      if (tr_aero) then ! volume-weighted aerosols
         do it = 1, n_aero
            trcr_depend(nt_aero+(it-1)*4  ) = 2 ! snow
            trcr_depend(nt_aero+(it-1)*4+1) = 2 ! snow
            trcr_depend(nt_aero+(it-1)*4+2) = 1 ! ice
            trcr_depend(nt_aero+(it-1)*4+3) = 1 ! ice
         enddo
      endif


  end subroutine init_cice_tracer_object_parent_indices

!-----------------------------------------------------------------------

end module cice_column


